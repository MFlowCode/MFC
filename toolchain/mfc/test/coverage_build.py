"""
File-level gcov coverage-map builder for MFC (master-side only).

Build MFC once with gfortran --coverage, run all tests individually, record
which .fpp files each test executes, and write a coverage map keyed by
param_hash (stable across cosmetic cases.py edits).

Workflow:
    ./mfc.sh build --gcov -j 8                        # one-time: instrumented build
    ./mfc.sh test --build-coverage-map --gcov -j 8    # one-time: populate the map
"""

import io
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from .. import common
from ..build import POST_PROCESS, PRE_PROCESS, SIMULATION
from ..common import MFCException
from ..printer import cons
from .case import POST_PROCESS_3D_PARAMS, get_post_process_mods, input_bubbles_lagrange


def _get_gcov_version(gcov_binary: str) -> str:
    """Return the version string from gcov --version."""
    try:
        result = subprocess.run([gcov_binary, "--version"], capture_output=True, text=True, timeout=10, check=False)
        for line in result.stdout.splitlines():
            if line.strip():
                return line.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return "unknown"


def find_gcov_binary() -> str:
    """
    Find a GNU gcov binary compatible with the system gfortran.

    On macOS with Homebrew GCC, the binary is gcov-{major} (e.g. gcov-15).
    On Linux with system GCC, plain gcov is usually correct.
    Apple LLVM's /usr/bin/gcov is incompatible with gfortran .gcda files.
    """
    # Determine gfortran major version
    major = None
    try:
        result = subprocess.run(["gfortran", "--version"], capture_output=True, text=True, timeout=10, check=False)
        m = re.search(r"(\d+)\.\d+\.\d+", result.stdout)
        if m:
            major = m.group(1)
    except (OSError, subprocess.SubprocessError):
        pass

    # Try versioned binary first (Homebrew macOS), then plain gcov
    candidates = []
    if major:
        candidates.append(f"gcov-{major}")
    candidates.append("gcov")

    for candidate in candidates:
        path = shutil.which(candidate)
        if path is None:
            continue
        try:
            result = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=10, check=False)
            version_out = result.stdout
            if "Apple LLVM" in version_out or "Apple clang" in version_out:
                continue  # Apple's gcov cannot parse GCC-generated .gcda files
            if "GCC" in version_out or "GNU" in version_out:
                return path
        except (OSError, subprocess.SubprocessError):
            continue

    raise MFCException(
        "GNU gcov not found. gcov is required for the coverage map.\n"
        "  On macOS (Homebrew):  brew install gcc\n"
        "  On Linux (Debian/Ubuntu): apt install gcc\n"
        "  On Linux (RHEL/CentOS):  yum install gcc\n"
        "Apple's /usr/bin/gcov is incompatible with gfortran .gcda files."
    )


def find_gcno_files(root_dir: str) -> list:
    """
    Walk build/ and return all .gcno files (excluding venv paths).
    Raises if none found (indicates build was not done with --gcov).
    """
    build_dir = Path(root_dir) / "build"
    gcno_files = [p for p in build_dir.rglob("*.gcno") if "venv" not in p.parts]
    if not gcno_files:
        raise MFCException("No .gcno files found. Build with --gcov instrumentation first:\n  ./mfc.sh build --gcov -j 8")
    return gcno_files


def _parse_gcov_json_output(raw_bytes: bytes, root_dir: str) -> Optional[set]:
    """
    Parse gcov JSON output and return the set of .fpp file paths with coverage.
    Handles both gzip-compressed (gcov 13+) and raw JSON (gcov 12) formats.
    Handles concatenated JSON objects from batched gcov calls (multiple .gcno
    files passed to a single gcov invocation).
    Only .fpp files with at least one executed line are included.
    """
    import gzip
    import json

    try:
        text = gzip.decompress(raw_bytes).decode("utf-8", errors="replace")
    except (gzip.BadGzipFile, OSError):
        text = raw_bytes.decode("utf-8", errors="replace")

    result = set()
    real_root = os.path.realpath(root_dir)
    parsed_any = False

    # Parse potentially concatenated JSON objects (one per .gcno file).
    decoder = json.JSONDecoder()
    pos = 0
    while pos < len(text):
        while pos < len(text) and text[pos] in " \t\n\r":
            pos += 1
        if pos >= len(text):
            break
        try:
            data, end_pos = decoder.raw_decode(text, pos)
            pos = end_pos
            parsed_any = True
        except json.JSONDecodeError:
            remaining = len(text) - pos
            if remaining > 0:
                cons.print(f"[yellow]Warning: gcov JSON parse error at offset {pos} ({remaining} bytes remaining) — coverage for this test is untrustworthy, omitting from map.[/yellow]")
            # A mid-stream parse error means the JSON stream was truncated or
            # corrupted.  A partial coverage set is untrustworthy: a .fpp that
            # would have been recorded in the missing portion would be silently
            # skipped by select_tests on future runs.  Return None so the caller
            # omits this test from the map entirely (conservatively included).
            return None

        for file_entry in data.get("files", []):
            file_path = file_entry.get("file", "")
            if not file_path.endswith(".fpp"):
                continue
            if any(line.get("count", 0) > 0 for line in file_entry.get("lines", [])):
                try:
                    rel_path = os.path.relpath(os.path.realpath(file_path), real_root)
                except ValueError:
                    rel_path = file_path
                # Only keep src/ paths — build/staging/ artifacts from
                # case-optimized builds are auto-generated and never
                # appear in PR diffs.
                if rel_path.startswith("src/"):
                    result.add(rel_path)

    # If no JSON was parsed at all (complete garbage input), return None
    # so the caller omits this test from the map (conservatively included
    # on future runs).  An empty set after successful parsing means the test
    # genuinely covers no .fpp files.
    if not parsed_any:
        return None

    return result


def _compute_gcov_prefix_strip(root_dir: str) -> str:
    """
    Compute GCOV_PREFIX_STRIP so .gcda files preserve the build/ tree.

    GCOV_PREFIX_STRIP removes N leading path components from the compile-time
    absolute .gcda path.  We strip all components of the MFC root directory
    so the prefix tree starts with ``build/staging/...``.
    """
    real_root = os.path.realpath(root_dir)
    return str(len(Path(real_root).parts) - 1)  # -1 excludes root '/'


def _collect_single_test_coverage(
    uuid: str,
    test_gcda: str,
    root_dir: str,
    gcov_bin: str,
) -> tuple:
    """
    Collect file-level coverage for a single test, fully self-contained.

    Copies .gcno files from the real build tree into the test's isolated
    .gcda directory (alongside the .gcda files), runs a batched gcov call,
    then removes the .gcno copies.  Each test has its own directory, so
    this is safe to call concurrently without touching the shared build tree.
    """
    build_subdir = os.path.join(test_gcda, "build")
    if not os.path.isdir(build_subdir):
        # No .gcda files produced — test may not have run or GCOV_PREFIX
        # was misconfigured.  Return None so the test is omitted from the
        # map (conservatively included on future runs).  The sanity check
        # at the end of build_coverage_map will catch systemic failures.
        cons.print(f"[yellow]Warning: No .gcda directory for {uuid} — GCOV_PREFIX may be misconfigured.[/yellow]")
        return uuid, None

    gcno_copies = []

    for dirpath, _, filenames in os.walk(build_subdir):
        for fname in filenames:
            if not fname.endswith(".gcda"):
                continue
            # Derive matching .gcno path in the real build tree
            gcda_path = os.path.join(dirpath, fname)
            rel = os.path.relpath(gcda_path, test_gcda)
            gcno_rel = rel[:-5] + ".gcno"
            gcno_src = os.path.join(root_dir, gcno_rel)
            if os.path.isfile(gcno_src):
                # Copy .gcno alongside .gcda in the test's isolated dir.
                # Wrap in try/except for NFS TOCTOU races (file may vanish
                # between isfile() and copy on networked filesystems).
                gcno_dst = os.path.join(dirpath, fname[:-5] + ".gcno")
                try:
                    shutil.copy2(gcno_src, gcno_dst)
                except OSError:
                    continue
                gcno_copies.append(gcno_dst)

    if not gcno_copies:
        # Genuinely no matching .gcno files — return empty list (not None).
        # None means "collection failed, conservatively include"; empty list
        # means "test produced .gcda but no .gcno matched", which is a real
        # (if unusual) result that should be cached as-is.
        return uuid, []

    # Batch: single gcov call for all .gcno files in this test.
    # Run from root_dir so source path resolution works correctly.
    cmd = [gcov_bin, "--json-format", "--stdout"] + gcno_copies
    try:
        proc = subprocess.run(cmd, capture_output=True, cwd=root_dir, timeout=120, check=False)
    except (subprocess.TimeoutExpired, subprocess.SubprocessError, OSError) as exc:
        cons.print(f"[yellow]Warning: gcov failed for {uuid}: {exc}[/yellow]")
        return uuid, None
    finally:
        for g in gcno_copies:
            try:
                os.remove(g)
            except OSError:
                pass

    if proc.returncode != 0 or not proc.stdout:
        if proc.returncode != 0:
            cons.print(f"[yellow]Warning: gcov exited {proc.returncode} for {uuid}[/yellow]")
        return uuid, None

    coverage = _parse_gcov_json_output(proc.stdout, root_dir)
    if coverage is None:
        # Decode failure — return None so the caller omits this test from
        # the map (absent entries are conservatively included).
        return uuid, None
    return uuid, sorted(coverage)


def _run_single_test_direct(test_info: dict, gcda_dir: str, strip: str) -> tuple:
    """
    Run a single test by invoking Fortran executables directly.

    Bypasses ``./mfc.sh run`` entirely (no Python startup, no Mako template
    rendering, no shell script generation).  Input files and binary paths are
    pre-computed by the caller.

    Returns (uuid, test_gcda_path, failures).
    """
    uuid = test_info["uuid"]
    test_dir = test_info["dir"]
    binaries = test_info["binaries"]  # ordered list of (target_name, bin_path)
    ppn = test_info["ppn"]

    test_gcda = os.path.join(gcda_dir, uuid)
    os.makedirs(test_gcda, exist_ok=True)

    env = {**os.environ, "GCOV_PREFIX": test_gcda, "GCOV_PREFIX_STRIP": strip}

    # MPI-compiled binaries must be launched via an MPI launcher (even ppn=1).
    # Use --bind-to none to avoid binding issues with concurrent launches.
    if shutil.which("mpirun"):
        mpi_cmd = ["mpirun", "--bind-to", "none", "-np", str(ppn)]
    elif shutil.which("srun"):
        mpi_cmd = ["srun", "--ntasks", str(ppn)]
    else:
        raise MFCException("No MPI launcher found (mpirun or srun). MFC binaries require an MPI launcher.\n  On Ubuntu: sudo apt install openmpi-bin\n  On macOS:  brew install open-mpi")

    failures = []
    for target_name, bin_path in binaries:
        if not os.path.isfile(bin_path):
            # Record missing binary as a failure and stop: downstream targets
            # depend on outputs from earlier ones (e.g. simulation needs the
            # grid from pre_process), so running them without a predecessor
            # produces misleading init-only gcda files.
            failures.append((target_name, "missing-binary", f"binary not found: {bin_path}"))
            break

        # Verify .inp file exists before running (diagnostic for transient
        # filesystem issues where the file goes missing between phases).
        inp_file = os.path.join(test_dir, f"{target_name}.inp")
        if not os.path.isfile(inp_file):
            failures.append((target_name, "missing-inp", f"{inp_file} not found before launch"))
            break

        cmd = mpi_cmd + [bin_path]
        try:
            result = subprocess.run(cmd, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env, cwd=test_dir, timeout=600)
            if result.returncode != 0:
                # Save last lines of output for debugging.  Stop here: a
                # failed pre_process/simulation leaves no valid outputs for
                # the next target, and running it produces spurious coverage.
                tail = "\n".join(result.stdout.strip().splitlines()[-15:])
                failures.append((target_name, result.returncode, tail))
                break
        except subprocess.TimeoutExpired:
            failures.append((target_name, "timeout", ""))
            break
        except (subprocess.SubprocessError, OSError) as exc:
            failures.append((target_name, str(exc), ""))
            break

    return uuid, test_gcda, failures


def _prepare_test(case) -> dict:
    """
    Prepare a test for direct execution: create directory, generate .inp
    files, and resolve binary paths.  All Python/toolchain overhead happens
    here (single-threaded) so the parallel phase is pure subprocess calls.

    Temporarily sets modified params on the case object (needed by
    get_dirpath/to_input_file/get_inp), then restores the original
    params in a finally block so callers can safely reuse the case list.
    """
    try:
        case.delete_output()
        case.create_directory()
    except OSError as exc:
        cons.print(f"[yellow]Warning: Failed to prepare test directory for {case.get_uuid()}: {exc}[/yellow]")
        raise

    # Lagrange bubble tests need input files generated before running.
    if case.params.get("bubbles_lagrange", "F") == "T":
        try:
            input_bubbles_lagrange(case)
        except Exception as exc:
            cons.print(f"[yellow]Warning: Failed to generate Lagrange bubble input for {case.get_uuid()}: {exc}[/yellow]")
            raise

    # Work on a copy so we don't permanently mutate the case object.
    params = dict(case.params)

    # Apply post_process output params so simulation writes data files that
    # post_process reads.  Mirrors the generated case.py logic that normally
    # runs via ./mfc.sh run (see POST_PROCESS_OUTPUT_PARAMS in case.py).
    params.update(get_post_process_mods(params))

    # Run only one timestep: we only need to know which source files are
    # *touched*, not verify correctness.  A single step exercises the key
    # code paths across all three executables while preventing heavy 3D tests
    # from timing out under gcov instrumentation (~10x slowdown).
    params["t_step_stop"] = 1

    # Adaptive-dt tests: post_process computes n_save = int(t_stop/t_save)+1
    # and iterates over that many save indices.  But with small t_step_stop
    # the simulation produces far fewer saves.  Clamp t_stop so post_process
    # only reads saves that actually exist.
    if params.get("cfl_adap_dt", "F") == "T":
        t_save = float(params.get("t_save", 1.0))
        params["t_stop"] = t_save  # n_save = 2: indices 0 and 1

    # Heavy 3D tests: remove vorticity output (omega_wrt + fd_order) for
    # 3D QBMM tests.  The regular test suite runs post_process only under
    # --test-all (which CI sets); heavy 3D QBMM configs are known to crash
    # post_process (exit code 2) when vorticity FD is enabled on large grids
    # with many QBMM variables.  Strip those params here so the coverage
    # build does not fail on those tests.
    if int(params.get("p", 0)) > 0 and params.get("qbmm", "F") == "T":
        for key in POST_PROCESS_3D_PARAMS:
            params.pop(key, None)

    # Temporarily set mutated params on the case object for get_dirpath(),
    # to_input_file(), and get_inp().  Always restore the original params
    # so build_coverage_map callers can safely reuse the case list.
    orig_params = case.params
    case.params = params
    try:
        test_dir = case.get_dirpath()
        input_file = case.to_input_file()

        # Write .inp files directly (no subprocess, no Mako templates).
        # Suppress console output from get_inp() to avoid one message per (test, target) pair.
        # Run all three executables to capture coverage across the full pipeline
        # (pre_process: grid/IC generation; simulation: RHS/time-stepper; post_process: field I/O).
        targets = [PRE_PROCESS, SIMULATION, POST_PROCESS]
        binaries = []
        # NOTE: not thread-safe — Phase 1 must remain single-threaded.
        orig_file = cons.raw.file
        cons.raw.file = io.StringIO()
        try:
            for target in targets:
                inp_content = case.get_inp(target)
                common.file_write(os.path.join(test_dir, f"{target.name}.inp"), inp_content)
                bin_path = target.get_install_binpath(input_file)
                binaries.append((target.name, bin_path))
        finally:
            cons.raw.file = orig_file
    finally:
        case.params = orig_params

    return {
        "uuid": case.get_uuid(),
        "dir": test_dir,
        "binaries": binaries,
        "ppn": getattr(case, "ppn", 1),
    }


def build_coverage_map(
    root_dir: str,
    cases: list,
    n_jobs: Optional[int] = None,
) -> None:
    """
    Build the file-level coverage map by running tests in parallel.

    Phase 1 — Prepare all tests: generate .inp files and resolve binary paths.
    This happens single-threaded so the parallel phase has zero Python overhead.

    Phase 2 — Run all tests concurrently.  Each worker invokes Fortran binaries
    directly (no ``./mfc.sh run``, no shell scripts).  Each test's GCOV_PREFIX
    points to an isolated directory so .gcda files don't collide.

    Phase 3 — For each test, temporarily copy .gcno files from the real build tree
    into the test's isolated .gcda directory, run gcov to collect which .fpp files
    had coverage, then remove the .gcno copies.

    Requires a prior ``--gcov`` build: ``./mfc.sh build --gcov -j 8``
    """
    gcov_bin = find_gcov_binary()
    gcno_files = find_gcno_files(root_dir)
    strip = _compute_gcov_prefix_strip(root_dir)

    if n_jobs is None:
        # Caller should pass n_jobs explicitly on SLURM systems;
        # os.cpu_count() may exceed the SLURM allocation.
        n_jobs = max(os.cpu_count() or 1, 1)
    # Cap Phase 2 test parallelism: each test spawns gcov-instrumented MPI
    # processes (~2-5 GB each under gcov).  Too many concurrent tests cause OOM.
    # Phase 3 gcov workers run at full n_jobs (gcov is lightweight by comparison).
    phase2_jobs = min(n_jobs, 16)
    cons.print(f"[bold]Building coverage map for {len(cases)} tests ({phase2_jobs} test workers, {n_jobs} gcov workers)...[/bold]")
    cons.print(f"[dim]Using gcov binary: {gcov_bin}[/dim]")
    cons.print(f"[dim]Found {len(gcno_files)} .gcno files[/dim]")
    cons.print(f"[dim]GCOV_PREFIX_STRIP={strip}[/dim]")
    cons.print()

    # Phase 1: Prepare all tests (single-threaded; scales linearly with test count).
    cons.print("[bold]Phase 1/3: Preparing tests...[/bold]")
    test_infos = []
    for i, case in enumerate(cases):
        try:
            test_infos.append(_prepare_test(case))
        except Exception as exc:
            cons.print(f"  [yellow]Warning: skipping {case.get_uuid()} — prep failed: {exc}[/yellow]")
        if (i + 1) % 100 == 0 or (i + 1) == len(cases):
            cons.print(f"  [{i + 1:3d}/{len(cases):3d}] prepared")
    cons.print()

    gcda_dir = tempfile.mkdtemp(prefix="mfc_gcov_")
    try:
        # Phase 2: Run all tests in parallel via direct binary invocation.
        cons.print("[bold]Phase 2/3: Running tests...[/bold]")
        test_results: dict = {}
        all_failures: dict = {}
        with ThreadPoolExecutor(max_workers=phase2_jobs) as pool:
            futures = {pool.submit(_run_single_test_direct, info, gcda_dir, strip): info for info in test_infos}
            for i, future in enumerate(as_completed(futures)):
                try:
                    uuid, test_gcda, failures = future.result()
                except Exception as exc:
                    info = futures[future]
                    cons.print(f"  [yellow]Warning: {info['uuid']} failed to run: {exc}[/yellow]")
                    continue
                if failures:
                    # A test that crashed mid-pipeline produced only partial .gcda
                    # files (e.g. simulation failed after pre_process ran).  That
                    # truncated coverage is untrustworthy: a later .fpp change that
                    # ran only in the missing stage would be incorrectly skipped.
                    # Record the failure for the warning summary but do NOT add this
                    # test to test_results — absent entries are conservatively
                    # included by select_tests (rung 5), never skipped.
                    all_failures[uuid] = failures
                    continue
                test_results[uuid] = test_gcda
                if (i + 1) % 50 == 0 or (i + 1) == len(test_infos):
                    cons.print(f"  [{i + 1:3d}/{len(test_infos):3d}] tests completed")

        if all_failures:
            cons.print()
            cons.print(f"[bold yellow]Warning: {len(all_failures)} tests had target failures:[/bold yellow]")
            for uuid, fails in sorted(all_failures.items()):
                fail_str = ", ".join(f"{t}={rc}" for t, rc, _ in fails)
                cons.print(f"  [yellow]{uuid}[/yellow]: {fail_str}")
                for target_name, _rc, tail in fails:
                    if tail:
                        cons.print(f"    {target_name} output (last 15 lines):")
                        for line in tail.splitlines():
                            cons.print(f"      {line}")

        # Diagnostic: verify .gcda files exist for at least one test.
        sample_uuid = next(iter(test_results), None)
        if sample_uuid:
            sample_gcda = test_results[sample_uuid]
            sample_build = os.path.join(sample_gcda, "build")
            if os.path.isdir(sample_build):
                gcda_count = sum(1 for _, _, fns in os.walk(sample_build) for f in fns if f.endswith(".gcda"))
                cons.print(f"[dim]Sample test {sample_uuid}: {gcda_count} .gcda files in {sample_build}[/dim]")
            else:
                cons.print(f"[yellow]Sample test {sample_uuid}: no build/ dir in {sample_gcda}[/yellow]")

        # Phase 3: Collect gcov coverage from each test's isolated .gcda directory.
        # .gcno files are temporarily copied alongside .gcda files, then removed.
        cons.print()
        cons.print("[bold]Phase 3/3: Collecting coverage...[/bold]")
        # Internal collection keyed by uuid (gcov machinery references uuids).
        cache: dict = {}
        completed = 0
        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            futures = {
                pool.submit(
                    _collect_single_test_coverage,
                    uuid,
                    test_gcda,
                    root_dir,
                    gcov_bin,
                ): uuid
                for uuid, test_gcda in test_results.items()
            }
            for future in as_completed(futures):
                try:
                    uuid, coverage = future.result()
                except Exception as exc:
                    uuid = futures[future]
                    cons.print(f"  [yellow]Warning: {uuid} coverage failed: {exc}[/yellow]")
                    # Do NOT store entry — absent entries are conservatively
                    # included by select_tests, while [] means "covers no files"
                    # and would permanently skip the test.
                    continue
                if coverage is None:
                    # Decode or collection failure — omit from map so the
                    # test is conservatively included on future runs.
                    continue
                cache[uuid] = coverage
                completed += 1
                if completed % 50 == 0 or completed == len(test_results):
                    cons.print(f"  [{completed:3d}/{len(test_results):3d}] tests processed")
    finally:
        try:
            shutil.rmtree(gcda_dir)
        except OSError as exc:
            cons.print(f"[yellow]Warning: Failed to clean up temp directory {gcda_dir}: {exc}[/yellow]")

    # Translate internal uuid keys -> stable param_hash keys for the committed map.
    from .coverage import COVERAGE_MAP_PATH, save_map

    key_by_uuid = {c.get_uuid(): c.coverage_key() for c in cases}
    out = {}
    for uuid, cov in cache.items():
        if uuid == "_meta" or cov is None:
            continue
        out[key_by_uuid.get(uuid, uuid)] = cov
    n_tests = sum(1 for v in out.values() if v)
    if n_tests == 0:
        raise MFCException("Coverage build produced zero coverage. Check the --gcov build and gcov binary.")
    git_sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=root_dir, check=False).stdout.strip()
    save_map(COVERAGE_MAP_PATH, out, n_tests=n_tests, git_sha=git_sha, gfortran_version=_get_gcov_version(gcov_bin))
    cons.print(f"[bold green]Coverage map written to {COVERAGE_MAP_PATH}[/bold green] ({n_tests} tests)")

    # Clean up test output directories from Phase 1/2 (grid files, restart files,
    # silo output, etc.).  These live on NFS scratch and can total several GB for
    # the full test suite.  Leaving them behind creates I/O pressure for subsequent
    # test jobs that share the same scratch filesystem.
    cons.print("[dim]Cleaning up test output directories...[/dim]")
    for case in cases:
        try:
            case.delete_output()
        except OSError:
            pass  # Best-effort; NFS errors are non-fatal here
