"""
File-level gcov coverage-based test pruning for MFC.

Build MFC once with gfortran --coverage, run all tests individually, record
which .fpp files each test executes, and cache that mapping.

When files change on a PR, intersect the changed .fpp files against each test's
covered file set. Only tests that touch at least one changed file run.

Workflow:
    ./mfc.sh build --gcov -j 8              # one-time: build with coverage
    ./mfc.sh test --build-coverage-cache    # one-time: populate the cache
    ./mfc.sh test --only-changes -j 8       # fast: run only affected tests
"""

import io
import os
import re
import json
import gzip
import shutil
import hashlib
import tempfile
import subprocess
import datetime
from pathlib import Path
from typing import Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

from ..printer import cons
from .. import common
from ..common import MFCException
from ..build import PRE_PROCESS, SIMULATION, POST_PROCESS, SYSCHECK
from .case import input_bubbles_lagrange


COVERAGE_CACHE_PATH = Path(common.MFC_ROOT_DIR) / "toolchain/mfc/test/test_coverage_cache.json.gz"

# Changes to these files trigger the full test suite.
# CPU coverage cannot tell us about GPU directive changes (macro files), and
# toolchain files define or change the set of tests themselves.
ALWAYS_RUN_ALL = frozenset([
    "src/common/include/parallel_macros.fpp",
    "src/common/include/acc_macros.fpp",
    "src/common/include/omp_macros.fpp",
    "src/common/include/shared_parallel_macros.fpp",
    "src/common/include/macros.fpp",
    "toolchain/mfc/test/cases.py",
    "toolchain/mfc/test/case.py",
    "toolchain/mfc/params/definitions.py",
    "toolchain/mfc/run/input.py",
    "toolchain/mfc/case_validator.py",
    "CMakeLists.txt",
])


def _get_gcov_version(gcov_binary: str) -> str:
    """Return the version string from gcov --version."""
    try:
        result = subprocess.run(
            [gcov_binary, "--version"],
            capture_output=True, text=True, timeout=10, check=False
        )
        for line in result.stdout.splitlines():
            if line.strip():
                return line.strip()
    except Exception:
        pass
    return "unknown"


def find_gcov_binary(_root_dir: str = "") -> str:  # pylint: disable=unused-argument
    """
    Find a GNU gcov binary compatible with the system gfortran.

    On macOS with Homebrew GCC, the binary is gcov-{major} (e.g. gcov-15).
    On Linux with system GCC, plain gcov is usually correct.
    Apple LLVM's /usr/bin/gcov is incompatible with gfortran .gcda files.
    """
    # Determine gfortran major version
    major = None
    try:
        result = subprocess.run(
            ["gfortran", "--version"],
            capture_output=True, text=True, timeout=10, check=False
        )
        m = re.search(r'(\d+)\.\d+\.\d+', result.stdout)
        if m:
            major = m.group(1)
    except Exception:
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
            result = subprocess.run(
                [path, "--version"],
                capture_output=True, text=True, timeout=10, check=False
            )
            version_out = result.stdout
            if "Apple LLVM" in version_out or "Apple clang" in version_out:
                continue  # Apple's gcov cannot parse GCC-generated .gcda files
            if "GCC" in version_out or "GNU" in version_out:
                return path
        except Exception:
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
    gcno_files = [
        p for p in build_dir.rglob("*.gcno")
        if "venv" not in p.parts
    ]
    if not gcno_files:
        raise MFCException(
            "No .gcno files found. Build with --gcov instrumentation first:\n"
            "  ./mfc.sh build --gcov -j 8"
        )
    return gcno_files


def zero_gcda_files(root_dir: str) -> None:
    """
    Delete all .gcda files under build/ (excluding venv).
    Called before each test run during cache building to isolate per-test coverage.
    """
    build_dir = Path(root_dir) / "build"
    for gcda in build_dir.rglob("*.gcda"):
        if "venv" not in gcda.parts:
            try:
                gcda.unlink()
            except OSError:
                pass


def _parse_gcov_json_output(raw_bytes: bytes, root_dir: str) -> set:
    """
    Parse gcov JSON output and return the set of .fpp file paths with coverage.
    Handles both gzip-compressed (gcov 13+) and raw JSON (gcov 12) formats.
    Only .fpp files with at least one executed line are included.
    """
    try:
        data = json.loads(gzip.decompress(raw_bytes))
    except (gzip.BadGzipFile, OSError):
        try:
            data = json.loads(raw_bytes)
        except (json.JSONDecodeError, ValueError):
            return set()
    except Exception:
        return set()

    result = set()
    real_root = os.path.realpath(root_dir)
    for file_entry in data.get("files", []):
        file_path = file_entry.get("file", "")
        if not file_path.endswith(".fpp"):
            continue
        if any(line.get("count", 0) > 0 for line in file_entry.get("lines", [])):
            try:
                rel_path = os.path.relpath(os.path.realpath(file_path), real_root)
            except ValueError:
                rel_path = file_path
            result.add(rel_path)

    return result


def collect_coverage_for_test(gcno_files: list, root_dir: str, gcov_binary: str) -> set:
    """
    Run gcov on all .gcno files and return the set of .fpp files with coverage.

    Expects .gcda files to be in their normal locations next to the .gcno files.
    """
    merged = set()

    for gcno_file in gcno_files:
        try:
            cmd = [gcov_binary, "--json-format", "--stdout", str(gcno_file)]
            proc = subprocess.run(
                cmd, capture_output=True, cwd=root_dir, timeout=60,
                check=False
            )
        except subprocess.TimeoutExpired:
            continue
        except Exception:
            continue

        if proc.returncode != 0 or not proc.stdout:
            continue

        merged.update(_parse_gcov_json_output(proc.stdout, root_dir))

    return merged


def _find_matching_gcno(root_dir: str) -> list:
    """
    Find .gcno files that have a matching .gcda in the build tree.

    After installing a test's .gcda files, only .gcno files with a sibling
    .gcda need gcov processing.  This typically reduces 414 .gcno files
    to ~50, giving an ~8x speedup.
    """
    build_dir = Path(root_dir) / "build"
    matching = []
    for gcda in build_dir.rglob("*.gcda"):
        if "venv" in gcda.parts:
            continue
        gcno = gcda.with_suffix(".gcno")
        if gcno.exists():
            matching.append(gcno)
    return matching


def _gcda_path_to_fpp(gcda_rel: str) -> str:
    """
    Map a .gcda relative path to the corresponding .fpp source path.

    Build tree layout:
        CMakeFiles/<target>.dir/fypp/<target>/<name>.fpp.f90.gcda -> src/<target>/<name>.fpp

    Returns empty string for non-.fpp files (plain .f90, modules/, ltrans).
    """
    # Extract path after CMakeFiles/<target>.dir/
    m = re.match(r'.*?CMakeFiles/[^/]+\.dir/(.*)', gcda_rel)
    if not m:
        return ""
    inner = m.group(1)  # e.g. fypp/simulation/m_rhs.fpp.f90.gcda

    # Only .fpp files: inner must contain ".fpp.f90.gcda"
    if ".fpp.f90.gcda" not in inner:
        return ""

    # fypp/<target>/<name>.fpp.f90.gcda -> src/<target>/<name>.fpp
    path = inner.replace(".f90.gcda", "")  # fypp/<target>/<name>.fpp
    if path.startswith("fypp/"):
        path = "src/" + path[5:]  # src/<target>/<name>.fpp
    return path


def collect_coverage_from_gcda(prefix_dir: str) -> set:
    """
    Infer file-level coverage from .gcda file existence in a GCOV_PREFIX tree.

    This is much faster than running gcov (instant vs minutes) because we
    only need to list files and map paths.  A .gcda file is created by
    gfortran's runtime for each compilation unit that executed at least one
    function, so its existence implies the source file had coverage.
    """
    result = set()
    build_subdir = os.path.join(prefix_dir, "build")
    if not os.path.isdir(build_subdir):
        return result
    for dirpath, _dirnames, filenames in os.walk(build_subdir):
        for fname in filenames:
            if not fname.endswith(".gcda"):
                continue
            full = os.path.join(dirpath, fname)
            rel = os.path.relpath(full, prefix_dir)
            fpp = _gcda_path_to_fpp(rel)
            if fpp:
                result.add(fpp)
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


def _install_gcda_files(prefix_dir: str, root_dir: str) -> int:
    """
    Copy .gcda files from a GCOV_PREFIX tree into the build directory.

    The prefix tree mirrors the build layout (e.g. ``<prefix>/build/staging/…``).
    Returns the number of files copied.
    """
    build_subdir = os.path.join(prefix_dir, "build")
    if not os.path.isdir(build_subdir):
        return 0
    count = 0
    for dirpath, _dirnames, filenames in os.walk(build_subdir):
        for fname in filenames:
            if not fname.endswith(".gcda"):
                continue
            src = os.path.join(dirpath, fname)
            rel = os.path.relpath(src, prefix_dir)
            dst = os.path.join(root_dir, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)
            count += 1
    return count


def _run_single_test_direct(test_info: dict, gcda_dir: str, strip: str) -> tuple:
    """
    Run a single test by invoking Fortran executables directly.

    Bypasses ``./mfc.sh run`` entirely (no Python startup, no Mako template
    rendering, no shell script generation).  Input files and binary paths are
    pre-computed by the caller.

    Returns (uuid, test_gcda_path).
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
        mpi_cmd = []

    for _, bin_path in binaries:
        if not os.path.isfile(bin_path):
            continue
        cmd = mpi_cmd + [bin_path]
        try:
            subprocess.run(cmd, check=False, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           env=env, cwd=test_dir, timeout=300)
        except Exception:
            pass

    return uuid, test_gcda


def _prepare_test(case, root_dir: str) -> dict:  # pylint: disable=unused-argument
    """
    Prepare a test for direct execution: create directory, generate .inp
    files, and resolve binary paths.  All Python/toolchain overhead happens
    here (single-threaded) so the parallel phase is pure subprocess calls.
    """
    try:
        case.delete_output()
        case.create_directory()
    except Exception:
        pass

    # Lagrange bubble tests need input files generated before running.
    if case.params.get("bubbles_lagrange", 'F') == 'T':
        try:
            input_bubbles_lagrange(case)
        except Exception:
            pass

    test_dir = case.get_dirpath()
    input_file = case.to_input_file()

    # Write .inp files directly (no subprocess, no Mako templates).
    # Suppress console output from get_inp() to avoid 555×4 messages.
    targets = [SYSCHECK, PRE_PROCESS, SIMULATION, POST_PROCESS]
    binaries = []
    orig_file = cons.raw.file
    cons.raw.file = io.StringIO()
    try:
        for target in targets:
            inp_content = case.get_inp(target)
            common.file_write(os.path.join(test_dir, f"{target.name}.inp"),
                              inp_content)
            bin_path = target.get_install_binpath(input_file)
            binaries.append((target.name, bin_path))
    finally:
        cons.raw.file = orig_file

    return {
        "uuid":     case.get_uuid(),
        "dir":      test_dir,
        "binaries": binaries,
        "ppn":      getattr(case, 'ppn', 1),
    }


def build_coverage_cache(  # pylint: disable=unused-argument,too-many-locals
    root_dir: str, cases: list, extra_args: list = None, n_jobs: int = None,
) -> None:
    """
    Build the file-level coverage cache by running tests in parallel.

    Phase 0 — Prepare all tests: generate .inp files and resolve binary paths.
    This happens single-threaded so the parallel phase has zero Python overhead.

    Phase 1 — Run all tests concurrently.  Each worker invokes Fortran binaries
    directly (no ``./mfc.sh run``, no shell scripts).  Each test's GCOV_PREFIX
    points to an isolated directory so .gcda files don't collide.

    Phase 2 — For each test, copy its .gcda tree into the real build directory,
    run gcov to collect which .fpp files had coverage, then remove the .gcda files.

    Requires a prior ``--gcov`` build: ``./mfc.sh build --gcov -j 8``
    """
    gcov_bin = find_gcov_binary(root_dir)
    gcno_files = find_gcno_files(root_dir)
    strip = _compute_gcov_prefix_strip(root_dir)

    if n_jobs is None:
        n_jobs = max(os.cpu_count() or 1, 1)
    cons.print(f"[bold]Building coverage cache for {len(cases)} tests "
               f"({n_jobs} parallel)...[/bold]")
    cons.print(f"[dim]Using gcov binary: {gcov_bin}[/dim]")
    cons.print(f"[dim]Found {len(gcno_files)} .gcno files[/dim]")
    cons.print(f"[dim]GCOV_PREFIX_STRIP={strip}[/dim]")
    cons.print()

    # Phase 0: Prepare all tests (single-threaded, ~30s for 555 tests).
    cons.print("[bold]Phase 0/2: Preparing tests...[/bold]")
    test_infos = []
    for i, case in enumerate(cases):
        test_infos.append(_prepare_test(case, root_dir))
        if (i + 1) % 100 == 0 or (i + 1) == len(cases):
            cons.print(f"  [{i+1:3d}/{len(cases):3d}] prepared")
    cons.print()

    gcda_dir = tempfile.mkdtemp(prefix="mfc_gcov_")

    # Phase 1: Run all tests in parallel via direct binary invocation.
    cons.print("[bold]Phase 1/2: Running tests...[/bold]")
    test_results: dict = {}
    with ThreadPoolExecutor(max_workers=n_jobs) as pool:
        futures = {
            pool.submit(_run_single_test_direct, info, gcda_dir, strip): info
            for info in test_infos
        }
        for i, future in enumerate(as_completed(futures)):
            uuid, test_gcda = future.result()
            test_results[uuid] = test_gcda
            if (i + 1) % 50 == 0 or (i + 1) == len(cases):
                cons.print(f"  [{i+1:3d}/{len(cases):3d}] tests completed")

    # Phase 2: Collect gcov coverage from each test's isolated .gcda directory.
    # For each test, copy its .gcda files into the build tree, run gcov only
    # on matching .gcno files (not all 414), then clean up.  Targeting matching
    # .gcno files gives ~8x speedup over the full scan.
    cons.print()
    cons.print("[bold]Phase 2/2: Collecting coverage...[/bold]")
    cache: dict = {}
    for i, (uuid, test_gcda) in enumerate(sorted(test_results.items())):
        zero_gcda_files(root_dir)
        n_copied = _install_gcda_files(test_gcda, root_dir)

        if n_copied == 0:
            coverage = set()
        else:
            # Only run gcov on .gcno files that have a matching .gcda installed.
            matching = _find_matching_gcno(root_dir)
            coverage = collect_coverage_for_test(
                matching or gcno_files, root_dir, gcov_bin
            )

        cache[uuid] = sorted(coverage)
        if (i + 1) % 50 == 0 or (i + 1) == len(cases):
            cons.print(f"  [{i+1:3d}/{len(cases):3d}] tests processed")

    zero_gcda_files(root_dir)

    # Clean up temp directory.
    shutil.rmtree(gcda_dir, ignore_errors=True)

    cases_py_path = Path(root_dir) / "toolchain/mfc/test/cases.py"
    cases_hash = hashlib.sha256(cases_py_path.read_bytes()).hexdigest()
    gcov_version = _get_gcov_version(gcov_bin)

    cache["_meta"] = {
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "cases_hash": cases_hash,
        "gcov_version": gcov_version,
    }

    with gzip.open(COVERAGE_CACHE_PATH, "wt", encoding="utf-8") as f:
        json.dump(cache, f, indent=2)

    cons.print()
    cons.print(f"[bold green]Coverage cache written to {COVERAGE_CACHE_PATH}[/bold green]")
    cons.print(f"[dim]Cache has {len(cases)} test entries.[/dim]")


def _normalize_cache(cache: dict) -> dict:
    """Convert old line-level cache format to file-level if needed.

    Old format: {uuid: {file: [lines], ...}, ...}
    New format: {uuid: [file, ...], ...}
    """
    for key, value in cache.items():
        if key == "_meta":
            continue
        if isinstance(value, dict):
            cache[key] = sorted(value.keys())
    return cache


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
    except (OSError, gzip.BadGzipFile, json.JSONDecodeError, UnicodeDecodeError):
        cons.print("[yellow]Warning: Coverage cache is unreadable or corrupt.[/yellow]")
        return None

    cases_py = Path(root_dir) / "toolchain/mfc/test/cases.py"
    current_hash = hashlib.sha256(cases_py.read_bytes()).hexdigest()
    stored_hash = cache.get("_meta", {}).get("cases_hash", "")

    if current_hash != stored_hash:
        cons.print("[yellow]Warning: Coverage cache is stale (cases.py changed).[/yellow]")
        return None

    return _normalize_cache(cache)


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
    merge_base_result = subprocess.run(
        ["git", "merge-base", compare_branch, "HEAD"],
        capture_output=True, text=True, cwd=root_dir, timeout=30, check=False
    )
    if merge_base_result.returncode != 0:
        return None
    merge_base = merge_base_result.stdout.strip()
    if not merge_base:
        return None

    diff_result = subprocess.run(
        ["git", "diff", merge_base, "HEAD", "--name-only", "--no-color"],
        capture_output=True, text=True, cwd=root_dir, timeout=30, check=False
    )
    if diff_result.returncode != 0:
        return None

    return _parse_diff_files(diff_result.stdout)


def should_run_all_tests(changed_files: set) -> bool:
    """
    Return True if any changed file is in ALWAYS_RUN_ALL.

    GPU macro files and toolchain files cannot be correctly analyzed by CPU
    coverage — changes to them must always trigger the full test suite.
    """
    return bool(changed_files & ALWAYS_RUN_ALL)


def filter_tests_by_coverage(
    cases: list, coverage_cache: dict, changed_files: set
) -> tuple:
    """
    Filter test cases to only those whose covered files overlap with changed files.

    Returns (cases_to_run, skipped_cases).

    Conservative behavior:
    - Test not in cache (newly added) -> include it
    - No changed .fpp files -> skip all tests
    - Test has incomplete coverage (no simulation files recorded but simulation
      files changed) -> include it (cache build likely failed for this test)
    """
    changed_fpp = {f for f in changed_files if f.endswith(".fpp")}
    if not changed_fpp:
        return [], list(cases)

    changed_sim = any(f.startswith("src/simulation/") for f in changed_fpp)

    to_run = []
    skipped = []

    for case in cases:
        uuid = case.get_uuid()
        test_files = coverage_cache.get(uuid)

        if test_files is None:
            # Test not in cache (e.g., newly added) -> conservative: include
            to_run.append(case)
            continue

        test_file_set = set(test_files)

        # If simulation files changed but this test has no simulation coverage,
        # include it conservatively — the cache build likely failed for this test.
        if changed_sim and not any(f.startswith("src/simulation/") for f in test_file_set):
            to_run.append(case)
            continue

        if test_file_set & changed_fpp:
            to_run.append(case)
        else:
            skipped.append(case)

    return to_run, skipped
