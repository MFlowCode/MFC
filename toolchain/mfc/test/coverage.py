"""
File-level gcov coverage-based test pruning for MFC.

Build MFC once with gfortran --coverage, run all tests individually, record
which .fpp files each test executes, and cache that mapping.

When files change on a PR, intersect the changed .fpp files against each test's
covered file set. Only tests that touch at least one changed file run.

Workflow:
    ./mfc.sh test --build-coverage-cache --gcov -j 8  # one-time: populate the cache
    ./mfc.sh test --only-changes -j 8                 # fast: run only affected tests
"""

import datetime
import gzip
import hashlib
import io
import json
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

COVERAGE_CACHE_PATH = Path(common.MFC_ROOT_DIR) / "toolchain/mfc/test/test_coverage_cache.json.gz"

# Changes to these files trigger the full test suite.
# CPU coverage cannot tell us about GPU directive changes (macro files), and
# toolchain files define or change the set of tests themselves.
# TEMP: CMakeLists.txt is omitted to test the cache rebuild pipeline in CI.
# Add "CMakeLists.txt" to this set before merge.
ALWAYS_RUN_ALL = frozenset(
    [
        "src/common/include/parallel_macros.fpp",
        "src/common/include/acc_macros.fpp",
        "src/common/include/omp_macros.fpp",
        "src/common/include/shared_parallel_macros.fpp",
        "src/common/include/macros.fpp",
        "src/common/include/case.fpp",
        "toolchain/mfc/test/case.py",
        "toolchain/mfc/test/cases.py",
        "toolchain/mfc/test/coverage.py",
        "toolchain/mfc/params/definitions.py",
        "toolchain/mfc/run/input.py",
        "toolchain/mfc/case_validator.py",
    ]
)

# Directory prefixes: any changed file under these paths triggers full suite.
# Note: src/simulation/include/ (.fpp files like inline_riemann.fpp) is NOT
# listed here — Fypp line markers (--line-marker-format=gfortran5) correctly
# attribute included file paths, so gcov coverage tracks them accurately.
ALWAYS_RUN_ALL_PREFIXES = ("toolchain/cmake/",)


def _get_gcov_version(gcov_binary: str) -> str:
    """Return the version string from gcov --version."""
    try:
        result = subprocess.run([gcov_binary, "--version"], capture_output=True, text=True, timeout=10, check=False)
        for line in result.stdout.splitlines():
            if line.strip():
                return line.strip()
    except Exception:
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
        "GNU gcov not found. gcov is required for the coverage cache.\n"
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


def _parse_gcov_json_output(raw_bytes: bytes, root_dir: str) -> set:
    """
    Parse gcov JSON output and return the set of .fpp file paths with coverage.
    Handles both gzip-compressed (gcov 13+) and raw JSON (gcov 12) formats.
    Handles concatenated JSON objects from batched gcov calls (multiple .gcno
    files passed to a single gcov invocation).
    Only .fpp files with at least one executed line are included.
    """
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
                cons.print(f"[yellow]Warning: gcov JSON parse error at offset {pos} ({remaining} bytes remaining) — partial coverage recorded for this test.[/yellow]")
            break

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
    # so the caller omits this test from the cache (conservatively included
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


def _collect_single_test_coverage(  # pylint: disable=too-many-locals
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
        # cache (conservatively included on future runs).  The sanity check
        # at the end of build_coverage_cache will catch systemic failures.
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
        return uuid, []

    # Batch: single gcov call for all .gcno files in this test.
    # Run from root_dir so source path resolution works correctly.
    cmd = [gcov_bin, "--json-format", "--stdout"] + gcno_copies
    try:
        proc = subprocess.run(cmd, capture_output=True, cwd=root_dir, timeout=120, check=False)
    except (subprocess.TimeoutExpired, subprocess.SubprocessError, OSError) as exc:
        cons.print(f"[yellow]Warning: gcov failed for {uuid}: {exc}[/yellow]")
        return uuid, []
    finally:
        for g in gcno_copies:
            try:
                os.remove(g)
            except OSError:
                pass

    if proc.returncode != 0 or not proc.stdout:
        if proc.returncode != 0:
            cons.print(f"[yellow]Warning: gcov exited {proc.returncode} for {uuid}[/yellow]")
        return uuid, []

    coverage = _parse_gcov_json_output(proc.stdout, root_dir)
    if coverage is None:
        # Decode failure — return None so the caller omits this test from
        # the cache (absent entries are conservatively included).
        return uuid, None
    return uuid, sorted(coverage)


def _run_single_test_direct(test_info: dict, gcda_dir: str, strip: str) -> tuple:  # pylint: disable=too-many-locals
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


def _prepare_test(case, root_dir: str) -> dict:  # pylint: disable=unused-argument,too-many-locals
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
    # 3D QBMM tests.  Normal test execution never runs post_process (only
    # PRE_PROCESS + SIMULATION, never POST_PROCESS), so post_process on
    # heavy 3D configs is untested.  Vorticity FD computation on large grids
    # with many QBMM variables causes post_process to crash (exit code 2).
    if int(params.get("p", 0)) > 0 and params.get("qbmm", "F") == "T":
        for key in POST_PROCESS_3D_PARAMS:
            params.pop(key, None)

    # Temporarily set mutated params on the case object for get_dirpath(),
    # to_input_file(), and get_inp().  Always restore the original params
    # so build_coverage_cache callers can safely reuse the case list.
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


def build_coverage_cache(  # pylint: disable=too-many-locals,too-many-statements
    root_dir: str,
    cases: list,
    n_jobs: int = None,
) -> None:
    """
    Build the file-level coverage cache by running tests in parallel.

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
    cons.print(f"[bold]Building coverage cache for {len(cases)} tests ({phase2_jobs} test workers, {n_jobs} gcov workers)...[/bold]")
    cons.print(f"[dim]Using gcov binary: {gcov_bin}[/dim]")
    cons.print(f"[dim]Found {len(gcno_files)} .gcno files[/dim]")
    cons.print(f"[dim]GCOV_PREFIX_STRIP={strip}[/dim]")
    cons.print()

    # Phase 1: Prepare all tests (single-threaded; scales linearly with test count).
    cons.print("[bold]Phase 1/3: Preparing tests...[/bold]")
    test_infos = []
    for i, case in enumerate(cases):
        try:
            test_infos.append(_prepare_test(case, root_dir))
        except Exception as exc:  # pylint: disable=broad-except
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
                except Exception as exc:  # pylint: disable=broad-except
                    info = futures[future]
                    cons.print(f"  [yellow]Warning: {info['uuid']} failed to run: {exc}[/yellow]")
                    continue
                test_results[uuid] = test_gcda
                if failures:
                    all_failures[uuid] = failures
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
                except Exception as exc:  # pylint: disable=broad-except
                    uuid = futures[future]
                    cons.print(f"  [yellow]Warning: {uuid} coverage failed: {exc}[/yellow]")
                    # Do NOT store entry — absent entries are conservatively
                    # included by filter_tests_by_coverage, while [] means
                    # "covers no files" and would permanently skip the test.
                    continue
                if coverage is None:
                    # Decode or collection failure — omit from cache so the
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

    # Sanity check: at least some tests should have non-empty coverage.
    tests_with_coverage = sum(1 for v in cache.values() if v)
    if tests_with_coverage == 0:
        raise MFCException("Coverage cache build produced zero coverage for all tests. Check that the build was done with --gcov and gcov is working correctly.")
    if tests_with_coverage < len(cases) // 2:
        cons.print(f"[bold yellow]Warning: Only {tests_with_coverage}/{len(cases)} tests have coverage data. Cache may be incomplete.[/bold yellow]")

    cases_py_path = Path(root_dir) / "toolchain/mfc/test/cases.py"
    try:
        cases_hash = hashlib.sha256(cases_py_path.read_bytes()).hexdigest()
    except OSError as exc:
        raise MFCException(f"Failed to read {cases_py_path} for cache metadata: {exc}") from exc
    gcov_version = _get_gcov_version(gcov_bin)

    cache["_meta"] = {
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "cases_hash": cases_hash,
        "gcov_version": gcov_version,
    }

    try:
        with gzip.open(COVERAGE_CACHE_PATH, "wt", encoding="utf-8") as f:
            json.dump(cache, f, indent=2)
    except OSError as exc:
        raise MFCException(f"Failed to write coverage cache to {COVERAGE_CACHE_PATH}: {exc}\nCheck disk space and filesystem permissions.") from exc

    cons.print()
    cons.print(f"[bold green]Coverage cache written to {COVERAGE_CACHE_PATH}[/bold green]")
    cons.print(f"[dim]Cache has {len(cases)} test entries.[/dim]")

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


def _normalize_cache(cache: dict) -> dict:
    """Convert old line-level cache format to file-level if needed.

    Old format: {uuid: {file: [lines], ...}, ...}
    New format: {uuid: [file, ...], ...}
    """
    result = {}
    for k, v in cache.items():
        if k == "_meta":
            result[k] = v
        elif isinstance(v, dict):
            result[k] = sorted(v.keys())
        elif isinstance(v, list):
            result[k] = v
        else:
            cons.print(f"[yellow]Warning: unexpected cache value type for {k}: {type(v).__name__} — omitting (test will be conservatively included).[/yellow]")
    return result


def load_coverage_cache(root_dir: str) -> Optional[dict]:
    """
    Load the coverage cache, returning None if missing or stale.

    Staleness is detected by comparing the SHA256 of cases.py at cache-build time
    against the current cases.py. Auto-converts old line-level format if needed.
    """
    if not COVERAGE_CACHE_PATH.exists():
        return None

    try:
        with gzip.open(COVERAGE_CACHE_PATH, "rt", encoding="utf-8") as f:
            cache = json.load(f)
    except (OSError, gzip.BadGzipFile, json.JSONDecodeError, UnicodeDecodeError) as exc:
        cons.print(f"[yellow]Warning: Coverage cache is unreadable or corrupt: {exc}[/yellow]")
        return None

    if not isinstance(cache, dict):
        cons.print("[yellow]Warning: Coverage cache has unexpected format.[/yellow]")
        return None

    cases_py = Path(root_dir) / "toolchain/mfc/test/cases.py"
    try:
        current_hash = hashlib.sha256(cases_py.read_bytes()).hexdigest()
    except FileNotFoundError:
        cons.print("[yellow]Warning: cases.py not found; cannot verify cache staleness.[/yellow]")
        return None
    stored_hash = cache.get("_meta", {}).get("cases_hash", "")

    if current_hash != stored_hash:
        cons.print("[yellow]Warning: Coverage cache is stale (cases.py changed).[/yellow]")
        return None

    cache = _normalize_cache(cache)

    # A cache with only _meta and no test entries is effectively empty
    # (e.g., corrupt rebuild).  Treat as missing so tests fall back to
    # the full suite with a clear message instead of silently including
    # every test as "not in cache".
    test_entries = {k for k in cache if k != "_meta"}
    if not test_entries:
        cons.print("[yellow]Warning: Coverage cache has no test entries.[/yellow]")
        return None

    return cache


def _parse_diff_files(diff_text: str) -> set:
    """
    Parse ``git diff --name-only`` output and return the set of changed file paths.
    """
    return {f for f in diff_text.strip().splitlines() if f}


def get_changed_files(root_dir: str, compare_branch: str = "master") -> Optional[set]:
    """
    Return the set of files changed in this branch relative to the merge-base
    with compare_branch, or None on git failure.

    Uses merge-base (not master tip) so that unrelated master advances don't
    appear as "your changes."
    """
    try:
        # Try local branch first, then origin/ remote ref (CI shallow clones).
        for ref in [compare_branch, f"origin/{compare_branch}"]:
            merge_base_result = subprocess.run(["git", "merge-base", ref, "HEAD"], capture_output=True, text=True, cwd=root_dir, timeout=30, check=False)
            if merge_base_result.returncode == 0:
                break
        else:
            return None
        merge_base = merge_base_result.stdout.strip()
        if not merge_base:
            return None

        diff_result = subprocess.run(["git", "diff", merge_base, "HEAD", "--name-only", "--no-color"], capture_output=True, text=True, cwd=root_dir, timeout=30, check=False)
        if diff_result.returncode != 0:
            return None

        return _parse_diff_files(diff_result.stdout)
    except (subprocess.TimeoutExpired, OSError):
        return None


def should_run_all_tests(changed_files: set) -> bool:
    """
    Return True if any changed file is in ALWAYS_RUN_ALL or under
    ALWAYS_RUN_ALL_PREFIXES.

    GPU macro files, Fypp includes, and build system files cannot be
    correctly analyzed by CPU coverage — changes to them must always
    trigger the full test suite.
    """
    if changed_files & ALWAYS_RUN_ALL:
        return True
    return any(f.startswith(ALWAYS_RUN_ALL_PREFIXES) for f in changed_files)


def filter_tests_by_coverage(cases: list, coverage_cache: dict, changed_files: set) -> tuple:
    """
    Filter test cases to only those whose covered files overlap with changed files.

    Returns (cases_to_run, skipped_cases).

    Conservative behavior:
    - Test not in cache (newly added) -> include it
    - Changed .f90/.f files under src/ -> run all tests (not in coverage cache)
    - No changed .fpp files -> skip all tests (this branch is unreachable from
      test.py, which handles the no-changed-fpp case before calling this function;
      retained as a safe fallback for direct callers)
    - Test has incomplete coverage (no simulation files recorded but simulation
      files changed) -> include it (cache build likely failed for this test)
    """
    changed_fpp = {f for f in changed_files if f.endswith(".fpp")}
    # .f90/.f files under src/ are real source but not in the coverage cache
    # (cache only tracks .fpp).  If any changed, run all tests conservatively.
    changed_f90 = {f for f in changed_files if f.startswith("src/") and (f.endswith(".f90") or f.endswith(".f"))}
    if changed_f90:
        return list(cases), []
    if not changed_fpp:
        return [], list(cases)

    changed_sim = any(f.startswith("src/simulation/") for f in changed_fpp)

    to_run = []
    skipped = []
    n_not_in_cache = 0
    n_no_sim_coverage = 0

    for case in cases:
        uuid = case.get_uuid()
        test_files = coverage_cache.get(uuid)

        if test_files is None:
            # Test not in cache (e.g., newly added) -> conservative: include
            to_run.append(case)
            n_not_in_cache += 1
            continue

        test_file_set = set(test_files)

        # If simulation files changed but this test has no simulation coverage,
        # include it conservatively — the cache build likely failed for this test.
        if changed_sim and not any(f.startswith("src/simulation/") for f in test_file_set):
            to_run.append(case)
            n_no_sim_coverage += 1
            continue

        if test_file_set & changed_fpp:
            to_run.append(case)
        else:
            skipped.append(case)

    if n_not_in_cache:
        cons.print(f"[dim]  {n_not_in_cache} test(s) included conservatively (not in cache)[/dim]")
    if n_no_sim_coverage:
        cons.print(f"[dim]  {n_no_sim_coverage} test(s) included conservatively (missing sim coverage)[/dim]")

    return to_run, skipped
