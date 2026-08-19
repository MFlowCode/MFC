import itertools
import math
import os
import shutil
import struct
import sys
import threading
import time
import typing
from random import sample, seed

import rich
import rich.table
from rich.panel import Panel

from .. import common, sched
from ..build import HDF5, POST_PROCESS, PRE_PROCESS, SIMULATION, build
from ..common import MFCException, does_command_exist, format_list_to_string, get_program_output
from ..packer import packer
from ..packer import tol as packtol
from ..printer import cons
from ..state import ARG
from .case import TestCase
from .cases import list_cases

nFAIL = 0
nPASS = 0
nSKIP = 0
current_test_number = 0
total_test_count = 0
errors = []
failed_tests = []  # Track failed test details for summary
test_start_time = None  # Track overall test duration

# Early abort thresholds
MIN_CASES_BEFORE_ABORT = 20
FAILURE_RATE_THRESHOLD = 0.3

# Per-test timeout (1 hour)
TEST_TIMEOUT_SECONDS = 3600

# Global abort flag for thread-safe early termination
# This flag is set when the failure rate exceeds the threshold, signaling
# all worker threads to exit gracefully. This avoids raising exceptions
# from worker threads which could leave the scheduler in an undefined state.
abort_tests = threading.Event()


class TestTimeoutError(MFCException):
    pass


# ib_state_*.dat record layout, written by s_write_serial_ib_state / s_write_parallel_ib_state in
# src/simulation/m_data_output.fpp: NFIELDS_PER_IB reals per IB, with x/y/z_centroid at fields 17:19.
_NFIELDS_PER_IB = 20
_CENTROID_SLICE = slice(16, 19)


def _read_ib_state_records(filepath: str, single: bool):
    if not os.path.isfile(filepath):
        raise MFCException(f"Expected IB state file does not exist: {filepath}")

    # ib_buf is real(wp) written with mpi_p, so the on-disk field width follows the build's working
    # precision. wp is single only under --single (--mixed keeps wp double, narrowing only stp), so this
    # keys off the build flag alone, not the case's `precision` (which selects database output format).
    field_size = 4 if single else 8
    fmt = "<" + ("f" if single else "d") * _NFIELDS_PER_IB
    record_size = _NFIELDS_PER_IB * field_size

    with open(filepath, "rb") as state_file:
        data = state_file.read()

    if len(data) % record_size != 0:
        raise MFCException(f"IB state file size is not a multiple of one IB record: {filepath}")

    return [struct.unpack(fmt, data[offset : offset + record_size]) for offset in range(0, len(data), record_size)]


def _assert_particle_cloud_non_overlap(case: TestCase, cloud_idx: int, records, start: int, count: int, tol: float):
    radius = case.params[f"particle_cloud({cloud_idx})%radius"]
    min_spacing = case.params.get(f"particle_cloud({cloud_idx})%min_spacing", 0.0)
    min_dist = 2.0 * radius + min_spacing

    for i in range(start, start + count):
        xi, yi, zi = records[i][_CENTROID_SLICE]
        for j in range(i + 1, start + count):
            xj, yj, zj = records[j][_CENTROID_SLICE]
            dist = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            if dist < min_dist - tol:
                raise MFCException(f"particle_cloud({cloud_idx}) particles overlap in ib_state_0.dat")


def _assert_particle_cloud_ib_state(case: TestCase):
    num_particle_clouds = case.params.get("num_particle_clouds", 0) or 0
    if num_particle_clouds <= 0 or case.params.get("ib_state_wrt", "F") != "T":
        return
    # file_per_process writes a different path and record layout (per-rank files with leading
    # num_local_ibs / gbl_patch_id integers); this reader only handles the single-file layout.
    if case.params.get("file_per_process", "F") == "T":
        return

    single = ARG("single")
    # Single precision carries ~1e-7 absolute error at O(1) coordinates, four orders past a 1e-12 slack, so
    # scale the geometric tolerance with the working precision to avoid spurious single/CI failures.
    tol = 1.0e-6 if single else 1.0e-12
    records = _read_ib_state_records(os.path.join(case.get_dirpath(), "restart_data", "ib_state_0.dat"), single)
    start = case.params.get("num_ibs", 0) or 0
    num_dims = 3 if (case.params.get("p", 0) or 0) > 0 else 2 if (case.params.get("n", 0) or 0) > 0 else 1

    for cloud_idx in range(1, num_particle_clouds + 1):
        geometry = case.params.get(f"particle_cloud({cloud_idx})%cloud_geometry", 1)
        packing_method = case.params.get(f"particle_cloud({cloud_idx})%packing_method", 1)
        count = case.params.get(f"particle_cloud({cloud_idx})%num_particles", 0) or 0
        radius = case.params[f"particle_cloud({cloud_idx})%radius"]
        records_end = start + count
        if records_end > len(records):
            raise MFCException(f"particle_cloud({cloud_idx}) expected {count} IB state records, found {len(records) - start}")

        # Box containment only holds for rejection sampling; lattice packing can place sites past the
        # requested region (issue #1730), so the bounds assertion is gated on packing_method == 1.
        if geometry == 1 and packing_method == 1:
            xc = case.params[f"particle_cloud({cloud_idx})%x_centroid"]
            yc = case.params.get(f"particle_cloud({cloud_idx})%y_centroid", 0.0)
            zc = case.params.get(f"particle_cloud({cloud_idx})%z_centroid", 0.0)
            lx = case.params[f"particle_cloud({cloud_idx})%length_x"]
            ly = case.params.get(f"particle_cloud({cloud_idx})%length_y", 0.0)
            lz = case.params.get(f"particle_cloud({cloud_idx})%length_z", 0.0)
            bounds = [(xc - lx / 2.0, xc + lx / 2.0), (yc - ly / 2.0, yc + ly / 2.0), (zc - lz / 2.0, zc + lz / 2.0)]
            for record in records[start:records_end]:
                for axis, coord in enumerate(record[_CENTROID_SLICE]):
                    if axis >= num_dims:
                        continue
                    lo, hi = bounds[axis]
                    if coord < lo - tol or coord > hi + tol:
                        raise MFCException(f"particle_cloud({cloud_idx}) box particle lies outside its cloud bounds in ib_state_0.dat")
        elif geometry == 2:
            xc = case.params[f"particle_cloud({cloud_idx})%x_centroid"]
            yc = case.params.get(f"particle_cloud({cloud_idx})%y_centroid", 0.0)
            zc = case.params.get(f"particle_cloud({cloud_idx})%z_centroid", 0.0)
            r_inner = case.params[f"particle_cloud({cloud_idx})%shell_inner_radius"] + radius
            r_outer = case.params[f"particle_cloud({cloud_idx})%shell_outer_radius"] - radius
            for record in records[start:records_end]:
                x, y, z = record[_CENTROID_SLICE]
                if num_dims < 3:
                    radial_dist = math.sqrt((x - xc) ** 2 + (y - yc) ** 2)
                    plane_coord = y - yc
                else:
                    radial_dist = math.sqrt((x - xc) ** 2 + (y - yc) ** 2 + (z - zc) ** 2)
                    plane_coord = z - zc
                if radial_dist < r_inner - tol or radial_dist > r_outer + tol:
                    raise MFCException(f"particle_cloud({cloud_idx}) shell particle violates radial clearance in ib_state_0.dat")
                if plane_coord < radius - tol:
                    raise MFCException(f"particle_cloud({cloud_idx}) shell particle violates flat-plane clearance in ib_state_0.dat")

        _assert_particle_cloud_non_overlap(case, cloud_idx, records, start, count, tol)
        start = records_end


def _filter_only(cases, skipped_cases):
    """Filter cases by --only terms using AND for labels, OR for UUIDs.

    Labels (non-UUID terms): case must match ALL labels (AND logic).
    UUIDs (8-char hex terms): case must match ANY UUID (OR logic).
    Mixed: keep case if all labels match OR any UUID matches.
    """

    def is_uuid(term):
        return len(term) == 8 and all(c in "0123456789abcdefABCDEF" for c in term)

    uuids = [t for t in ARG("only") if is_uuid(t)]
    labels = [t for t in ARG("only") if not is_uuid(t)]

    for case in cases[:]:
        check = set(case.trace.split(" -> "))
        check.add(case.get_uuid())

        label_ok = all(label in check for label in labels) if labels else True
        uuid_ok = any(u in check for u in uuids) if uuids else True

        if labels and uuids:
            keep = label_ok or uuid_ok
        elif labels:
            keep = label_ok
        else:
            keep = uuid_ok

        if not keep:
            cases.remove(case)
            skipped_cases.append(case)

    return cases, skipped_cases


def __filter(cases_) -> typing.Tuple[typing.List[TestCase], typing.List[TestCase]]:
    cases = cases_[:]
    selected_cases = []
    skipped_cases = []

    # Check "--from" and "--to" exist and are in the right order
    bFoundFrom, bFoundTo = (False, False)
    from_i = -1
    for i, case in enumerate(cases):
        if case.get_uuid() == ARG("from"):
            from_i = i
            bFoundFrom = True
            # Do not "continue" because "--to" might be the same as "--from"
        if bFoundFrom and case.get_uuid() == ARG("to"):
            cases = cases[from_i : i + 1]
            skipped_cases = [case for case in cases_ if case not in cases]
            bFoundTo = True
            break

    if not bFoundTo:
        raise MFCException("Testing: Your specified range [--from,--to] is incorrect. Please ensure both IDs exist and are in the correct order.")

    if len(ARG("only")) > 0:
        cases, skipped_cases = _filter_only(cases, skipped_cases)

        if not cases:
            raise MFCException(f"--only filter matched zero test cases. Specified: {ARG('only')}. Check that UUIDs/names are valid.")

    # Convergence cases are slow (multiple resolutions × MPI ranks). Skip
    # unless the user explicitly opted in via --only "Convergence" or a
    # specific convergence UUID. A label like --only "2D" must not
    # accidentally pull in "Convergence -> 2D -> ..." cases.
    if not ARG("list"):

        def is_uuid(term):
            return len(term) == 8 and all(c in "0123456789abcdefABCDEF" for c in term)

        only_terms = ARG("only")
        only_labels = [t for t in only_terms if not is_uuid(t)]
        only_uuids = [t for t in only_terms if is_uuid(t)]

        convergence_uuids = {c.get_uuid() for c in cases if getattr(c, "kind", "golden") == "convergence"}
        user_wants_convergence = "Convergence" in only_labels or any(u in convergence_uuids for u in only_uuids)

        if not user_wants_convergence:
            convergence_cases = [c for c in cases if getattr(c, "kind", "golden") == "convergence"]
            for c in convergence_cases:
                cases.remove(c)
                skipped_cases.append(c)

    for case in cases[:]:
        if case.ppn > 1 and not ARG("mpi"):
            cases.remove(case)
            skipped_cases.append(case)

    for case in cases[:]:
        if ARG("single"):
            skip = ["low_Mach", "Hypoelasticity", "teno", "Chemistry", "Phase Change model 6", "Axisymmetric", "Transducer", "Transducer Array", "Cylindrical", "HLLD", "Example"]
            if any(label in case.trace for label in skip):
                cases.remove(case)
                skipped_cases.append(case)

    for case in cases[:]:
        if ARG("gpu"):
            skip = ["Gauss Seidel"]
            if any(label in case.trace for label in skip):
                cases.remove(case)

    # Skip tests that fail under nvfortran in Docker (pass natively/Apptainer):
    #  - 3D_rayleigh_taylor_muscl: segfaults with nvfortran+MPI (seccomp/mprotect)
    if os.environ.get("FC") == "nvfortran" and os.path.exists("/.dockerenv"):
        nvhpc_skip_uuids = {}
        nvhpc_skip_traces = {"rayleigh_taylor_muscl"}
        for case in cases[:]:
            if case.get_uuid() in nvhpc_skip_uuids or any(t in case.trace for t in nvhpc_skip_traces):
                cases.remove(case)
                skipped_cases.append(case)

    if ARG("no_examples"):
        example_cases = [case for case in cases if "Example" in case.trace]
        skipped_cases += example_cases
        cases = [case for case in cases if case not in example_cases]

    if ARG("shard") is not None:
        parts = ARG("shard").split("/")
        if len(parts) != 2 or not all(p.isdigit() for p in parts) or int(parts[1]) < 1 or not 1 <= int(parts[0]) <= int(parts[1]):
            raise MFCException(f"Invalid --shard '{ARG('shard')}': expected 'i/n' with 1 <= i <= n (e.g., '1/2').")
        shard_idx, shard_count = int(parts[0]), int(parts[1])
        skipped_cases += [c for i, c in enumerate(cases) if i % shard_count != shard_idx - 1]
        cases = [c for i, c in enumerate(cases) if i % shard_count == shard_idx - 1]

        if not cases:
            raise MFCException(f"--shard {ARG('shard')} matched zero test cases. Total cases before sharding may be less than shard count.")

    if ARG("only_changes") or ARG("select_enforce"):
        import datetime

        from .. import common
        from .coverage import COVERAGE_MAP_PATH, format_summary, get_changed_files, load_map, select_tests

        entries, meta = load_map(COVERAGE_MAP_PATH)
        if entries is None:
            cons.print("[yellow]Coverage selection: map missing/corrupt — running full suite.[/yellow]")
        else:
            changed = get_changed_files(common.MFC_ROOT_DIR, ARG("changes_branch"), explicit=ARG("changed_files"))
            to_run, to_skip, reason = select_tests(cases, entries, changed)
            cons.print(format_summary(ran=len(to_run), total=len(cases), reason=reason, meta=meta, now=datetime.datetime.now(datetime.timezone.utc).isoformat()))
            if ARG("select_enforce"):
                skipped_cases += to_skip
                cases = to_run
            else:
                cons.print("[dim](shadow mode: running full suite; pass --select-enforce to actually skip)[/dim]")

    if ARG("percent") == 100:
        return cases, skipped_cases

    seed(time.time())

    selected_cases = sample(cases, k=int(len(cases) * ARG("percent") / 100.0))
    skipped_cases += [item for item in cases if item not in selected_cases]

    return selected_cases, skipped_cases


def test():
    global nFAIL, nPASS, nSKIP, total_test_count  # noqa: PLW0603
    global errors, failed_tests, test_start_time  # noqa: PLW0603

    test_start_time = time.time()  # Start timing
    failed_uuids_path = os.path.join(common.MFC_TEST_DIR, "failed_uuids.txt")
    cases = list_cases()

    # Delete UUIDs that are not in the list of cases from tests/
    if ARG("remove_old_tests"):
        dir_uuids = set(os.listdir(common.MFC_TEST_DIR))
        new_uuids = {case.get_uuid() for case in cases}

        for old_uuid in dir_uuids - new_uuids:
            cons.print(f"[bold red]Deleting:[/bold red] {old_uuid}")
            common.delete_directory(f"{common.MFC_TEST_DIR}/{old_uuid}")

        return

    if ARG("build_coverage_map"):
        from .coverage_build import build_coverage_map

        # Convergence tests are order-of-accuracy checks driven by convergence.py,
        # which fills in grid resolution and patch geometry per refinement level at
        # runtime. Their base params are skeletons (m=n=p=0, no geometry) that cannot
        # run standalone, so the direct-invocation collector cannot map them. Exclude
        # them: absent from the map, select_tests conservatively always-runs them
        # (rung 5), which is the desired behavior for convergence checks anyway.
        convergence = [b for b in cases if getattr(b, "kind", "golden") == "convergence"]
        coverage_cases = [b for b in cases if getattr(b, "kind", "golden") != "convergence"]
        if convergence:
            cons.print(f"[yellow]Excluding {len(convergence)} convergence tests from the coverage map (always-run by design).[/yellow]")

        all_cases = [b.to_case() for b in coverage_cases]
        unique = set()
        for case, code in itertools.product(all_cases, [PRE_PROCESS, SIMULATION, POST_PROCESS]):
            slug = code.get_slug(case.to_input_file())
            if slug not in unique:
                build(code, case.to_input_file())
                unique.add(slug)
        build_coverage_map(common.MFC_ROOT_DIR, all_cases, n_jobs=int(ARG("jobs")))
        return

    cases, skipped_cases = __filter(cases)
    cases = [_.to_case() for _ in cases]
    total_test_count = len(cases)

    if ARG("list"):
        table = rich.table.Table(title="MFC Test Cases", box=rich.table.box.SIMPLE)

        table.add_column("UUID", style="bold magenta", justify="center")
        table.add_column("Trace")

        for case in cases:
            table.add_row(case.get_uuid(), case.trace)

        rich.print(table)

        return

    # Some cases require a specific build of MFC for features like Chemistry,
    # Analytically defined patches, and --case-optimization. Here, we build all
    # the unique versions of MFC we need to run cases.
    codes = [PRE_PROCESS, SIMULATION] + ([POST_PROCESS] if ARG("test_all") else [])
    unique_builds = set()
    for case, code in itertools.product(cases, codes):
        slug = code.get_slug(case.to_input_file())
        if slug not in unique_builds:
            build(code, case.to_input_file())
            unique_builds.add(slug)

    cons.print()

    range_str = f"from [bold magenta]{ARG('from')}[/bold magenta] to [bold magenta]{ARG('to')}[/bold magenta]"

    if len(ARG("only")) > 0:
        range_str = "Only " + format_list_to_string(ARG("only"), "bold magenta", "Nothing to run")

    cons.print(f"[bold]Test {format_list_to_string([x.name for x in codes], 'magenta')}[/bold] | {range_str} ({len(cases)} test{'s' if len(cases) != 1 else ''})")
    cons.indent()

    # Run cases with multiple threads (if available)
    cons.print()
    cons.print("  Progress      Test Name                                        Time(s)   UUID")
    cons.print()

    # Select the correct number of threads to use to launch test cases
    # We can't use ARG("jobs") when the --case-optimization option is set
    # because running a test case may cause it to rebuild, and thus
    # interfere with the other test cases. It is a niche feature so we won't
    # engineer around this issue (for now).
    sched.sched([sched.Task(ppn=case.ppn, func=handle_case, args=[case], load=case.get_cell_count()) for case in cases], ARG("jobs"), ARG("gpus"))

    # Check if we aborted due to high failure rate
    if abort_tests.is_set():
        # Clean up stale failed_uuids.txt so CI doesn't retry wrong tests
        try:
            if os.path.exists(failed_uuids_path):
                os.remove(failed_uuids_path)
        except OSError:
            pass

        total_completed = nFAIL + nPASS
        cons.print()
        cons.unindent()
        if total_completed > 0:
            raise MFCException(f"Excessive test failures: {nFAIL}/{total_completed} failed ({nFAIL / total_completed * 100:.1f}%)")
        raise MFCException(f"Excessive test failures: {nFAIL} failed, but no tests were completed.")

    nSKIP = len(skipped_cases)
    cons.print()
    cons.unindent()

    # Calculate total test duration
    total_duration = time.time() - test_start_time
    minutes = int(total_duration // 60)
    seconds = total_duration % 60

    # Build the summary report
    _print_test_summary(nPASS, nFAIL, nSKIP, minutes, seconds, failed_tests, skipped_cases)

    # Write failed UUIDs to file for CI retry logic
    if failed_tests:
        with open(failed_uuids_path, "w") as f:
            for test_info in failed_tests:
                f.write(test_info["uuid"] + "\n")
    elif os.path.exists(failed_uuids_path):
        os.remove(failed_uuids_path)

    sys.exit(nFAIL)


def _print_test_summary(passed: int, failed: int, skipped: int, minutes: int, seconds: float, failed_test_list: list, _skipped_cases: list):
    """Print a comprehensive test summary report."""
    total = passed + failed + skipped

    # Build summary header
    if failed == 0:
        status_icon = "[bold green]✓[/bold green]"
        status_text = "[bold green]ALL TESTS PASSED[/bold green]"
        border_style = "green"
    else:
        status_icon = "[bold red]✗[/bold red]"
        status_text = f"[bold red]{failed} TEST{'S' if failed != 1 else ''} FAILED[/bold red]"
        border_style = "red"

    # Format time string
    if minutes > 0:
        time_str = f"{minutes}m {seconds:.1f}s"
    else:
        time_str = f"{seconds:.1f}s"

    # Build summary content
    summary_lines = [
        f"{status_icon} {status_text}",
        "",
        f"  [bold green]{passed:4d}[/bold green] passed",
        f"  [bold red]{failed:4d}[/bold red] failed",
        f"  [bold yellow]{skipped:4d}[/bold yellow] skipped",
        f"  [dim]{'─' * 12}[/dim]",
        f"  [bold]{total:4d}[/bold] total",
        "",
        f"  [dim]Time: {time_str}[/dim]",
    ]

    # Add failed tests details if any
    if failed_test_list:
        summary_lines.append("")
        summary_lines.append("  [bold red]Failed Tests:[/bold red]")
        for test_info in failed_test_list[:10]:  # Limit to first 10
            trace = test_info.get("trace", "Unknown")
            uuid = test_info.get("uuid", "Unknown")
            error_type = test_info.get("error_type", "")
            if len(trace) > 40:
                trace = trace[:37] + "..."
            summary_lines.append(f"    [red]•[/red] {trace}")
            summary_lines.append(f"      [dim]UUID: {uuid}[/dim]")
            if error_type:
                summary_lines.append(f"      [dim]({error_type})[/dim]")
        if len(failed_test_list) > 10:
            summary_lines.append(f"    [dim]... and {len(failed_test_list) - 10} more[/dim]")

    # Add next steps for failures
    if failed > 0:
        summary_lines.append("")
        summary_lines.append("  [bold]Next Steps:[/bold]")
        summary_lines.append("    • Run with [cyan]--generate[/cyan] to update golden files (if changes are intentional)")
        summary_lines.append("    • Check individual test output in [cyan]tests/<UUID>/[/cyan]")
        summary_lines.append("    • Run specific test: [cyan]./mfc.sh test --only <UUID>[/cyan]")

    cons.print()
    cons.raw.print(Panel("\n".join(summary_lines), title="[bold]Test Summary[/bold]", border_style=border_style, padding=(1, 2)))
    cons.print()


def _process_silo_file(silo_filepath: str, case: TestCase, out_filepath: str):
    """Process a single silo file with h5dump and check for NaNs/Infinities."""
    h5dump = f"{HDF5.get_install_dirpath(case.to_input_file())}/bin/h5dump"

    if not os.path.exists(h5dump or ""):
        if not does_command_exist("h5dump"):
            raise MFCException("h5dump couldn't be found.")
        h5dump = shutil.which("h5dump")

    output, err = get_program_output([h5dump, silo_filepath])

    if err != 0:
        raise MFCException(f"Test {case}: Failed to run h5dump. You can find the run's output in {out_filepath}, and the case dictionary in {case.get_filepath()}.")

    if "nan," in output:
        raise MFCException(f"Test {case}: Post Process has detected a NaN. You can find the run's output in {out_filepath}, and the case dictionary in {case.get_filepath()}.")

    if "inf," in output:
        raise MFCException(f"Test {case}: Post Process has detected an Infinity. You can find the run's output in {out_filepath}, and the case dictionary in {case.get_filepath()}.")


def _handle_convergence_case(case: TestCase, start_time: float):
    """Dispatch convergence/order-of-accuracy cases through convergence.py."""
    from .convergence import run_case

    if ARG("dry_run"):
        trace_display = case.trace if len(case.trace) <= 50 else case.trace[:47] + "..."
        cons.print(f"  (dry-run)     {trace_display:50s}   SKIP    [magenta]{case.get_uuid()}[/magenta]")
        return

    passed, output = run_case(case.convergence_spec)

    log_dir = os.path.join(common.MFC_TEST_DIR, case.get_uuid())
    common.create_directory(log_dir)
    common.file_write(os.path.join(log_dir, "convergence.log"), output)

    duration = time.time() - start_time
    global current_test_number  # noqa: PLW0603
    current_test_number += 1
    progress_str = f"({current_test_number:3d}/{total_test_count:3d})"
    trace_display = case.trace if len(case.trace) <= 50 else case.trace[:47] + "..."
    cons.print(f"  {progress_str}    {trace_display:50s}  {duration:6.2f}    [magenta]{case.get_uuid()}[/magenta]")

    if not passed:
        raise MFCException(f"Test {case}: convergence rate check failed (see {log_dir}/convergence.log)")


def _handle_case(case: TestCase, devices: typing.Set[int]):
    global current_test_number  # noqa: PLW0603
    start_time = time.time()

    if getattr(case, "kind", "golden") == "convergence":
        _handle_convergence_case(case, start_time)
        return

    # Set timeout using threading.Timer (works in worker threads)
    # Note: we intentionally do not use signal.alarm() here because signals
    # only work in the main thread; sched.sched runs tests in worker threads.
    # threading.Timer works correctly in this threaded context.
    timeout_flag = threading.Event()
    timeout_timer = threading.Timer(TEST_TIMEOUT_SECONDS, timeout_flag.set)
    timeout_timer.start()

    tol = case.compute_tolerance()
    case.delete_output()
    case.create_directory()

    if ARG("dry_run"):
        # Truncate long traces for readability
        trace_display = case.trace if len(case.trace) <= 50 else case.trace[:47] + "..."
        cons.print(f"  (dry-run)     {trace_display:50s}   SKIP    [magenta]{case.get_uuid()}[/magenta]")
        timeout_timer.cancel()
        return

    try:
        # Check timeout before starting
        if timeout_flag.is_set():
            raise TestTimeoutError("Test case exceeded 1 hour timeout")
        cmd = case.run([PRE_PROCESS, SIMULATION], gpus=devices)

        # Check timeout after simulation
        if timeout_flag.is_set():
            raise TestTimeoutError("Test case exceeded 1 hour timeout")

        out_filepath = os.path.join(case.get_dirpath(), "out_pre_sim.txt")

        common.file_write(out_filepath, cmd.stdout)

        if cmd.returncode != 0:
            cons.print(cmd.stdout)
            raise MFCException(f"Test {case}: Failed to execute MFC.")

        _assert_particle_cloud_ib_state(case)

        pack, err = packer.pack(case.get_dirpath())
        if err is not None:
            raise MFCException(f"Test {case}: {err}")

        if pack.has_bad_values():
            raise MFCException(f"Test {case}: NaN or Inf detected in the case.")

        golden_filepath = os.path.join(case.get_dirpath(), "golden.txt")
        if ARG("generate"):
            common.delete_file(golden_filepath)
            pack.save(golden_filepath)
        else:
            if not os.path.isfile(golden_filepath):
                raise MFCException(f"Test {case}: The golden file does not exist! To generate golden files, use the '--generate' flag.")

            golden = packer.load(golden_filepath)

            if ARG("add_new_variables"):
                for pfilepath, pentry in list(pack.entries.items()):
                    if golden.find(pfilepath) is None:
                        golden.set(pentry)

                for gfilepath, gentry in list(golden.entries.items()):
                    if pack.find(gfilepath) is None:
                        golden.remove(gentry)

                golden.save(golden_filepath)
            else:
                err, msg = packtol.compare(pack, packer.load(golden_filepath), packtol.Tolerance(tol, tol))
                if msg is not None:
                    raise MFCException(f"Test {case}: {msg}")

        # Restart roundtrip verification: run to midpoint, restart,
        # and compare restarted output against the straight run.
        if case.restart_check and not ARG("add_new_variables") and not ARG("generate"):
            straight_pack = pack

            if timeout_flag.is_set():
                raise TestTimeoutError("Test case exceeded 1 hour timeout")

            restart_result = case.run_restart([PRE_PROCESS, SIMULATION], devices)

            if timeout_flag.is_set():
                raise TestTimeoutError("Test case exceeded 1 hour timeout")

            out_filepath_restart = os.path.join(case.get_dirpath(), "out_restart.txt")
            common.file_write(out_filepath_restart, restart_result.stdout)

            if restart_result.returncode != 0:
                cons.print(restart_result.stdout)
                raise MFCException(f"Test {case}: Restart roundtrip run failed.")

            restart_pack, restart_err = packer.pack(case.get_dirpath())
            if restart_err is not None:
                raise MFCException(f"Test {case}: Restart pack error: {restart_err}")

            if restart_pack.has_bad_values():
                raise MFCException(f"Test {case}: NaN or Inf detected in restarted output.")

            _, restart_msg = packtol.compare(restart_pack, straight_pack, packtol.Tolerance(tol, tol))
            if restart_msg is not None:
                raise MFCException(f"Test {case}: Restart roundtrip mismatch: {restart_msg}")

        if ARG("test_all"):
            case.delete_output()
            # Check timeout before launching the (potentially long) post-process run
            if timeout_flag.is_set():
                raise TestTimeoutError("Test case exceeded 1 hour timeout")
            cmd = case.run([PRE_PROCESS, SIMULATION, POST_PROCESS], gpus=devices)
            out_filepath = os.path.join(case.get_dirpath(), "out_post.txt")
            common.file_write(out_filepath, cmd.stdout)

            silo_dir = os.path.join(case.get_dirpath(), "silo_hdf5", "p0")
            if os.path.isdir(silo_dir):
                for silo_filename in os.listdir(silo_dir):
                    silo_filepath = os.path.join(silo_dir, silo_filename)
                    _process_silo_file(silo_filepath, case, out_filepath)

        case.delete_output()

        end_time = time.time()
        duration = end_time - start_time

        current_test_number += 1
        progress_str = f"({current_test_number:3d}/{total_test_count:3d})"
        # Truncate long traces for readability, showing test name prominently
        trace_display = case.trace if len(case.trace) <= 50 else case.trace[:47] + "..."
        cons.print(f"  {progress_str}    {trace_display:50s}  {duration:6.2f}    [magenta]{case.get_uuid()}[/magenta]")

    except TestTimeoutError as exc:
        log_path = os.path.join(case.get_dirpath(), "out_pre_sim.txt")
        if os.path.exists(log_path):
            log_msg = f"Check the log at: {log_path}"
        else:
            log_msg = f"Log file ({log_path}) may not exist if the timeout occurred early."
        raise MFCException(f"Test {case} exceeded 1 hour timeout. This may indicate a hung simulation or misconfigured case. {log_msg}") from exc
    finally:
        timeout_timer.cancel()  # Cancel timeout timer


def handle_case(case: TestCase, devices: typing.Set[int]):
    global nFAIL, nPASS, nSKIP  # noqa: PLW0603
    global errors, failed_tests  # noqa: PLW0603

    # Check if we should abort before processing this case
    if abort_tests.is_set():
        return  # Exit gracefully if abort was requested

    nAttempts = 0
    if ARG("single"):
        max_attempts = max(ARG("max_attempts"), 3)
    else:
        max_attempts = ARG("max_attempts")

    while True:
        nAttempts += 1

        try:
            _handle_case(case, devices)
            if ARG("dry_run"):
                nSKIP += 1
            else:
                nPASS += 1
        except Exception as exc:
            if nAttempts < max_attempts:
                continue
            nFAIL += 1

            # Enhanced real-time failure feedback
            trace_display = case.trace if len(case.trace) <= 50 else case.trace[:47] + "..."
            cons.print()
            cons.print(f"  [bold red]✗ FAILED:[/bold red] {trace_display}")
            cons.print(f"    UUID: [magenta]{case.get_uuid()}[/magenta]")
            cons.print(f"    Attempts: {nAttempts}")

            cons.print(f"    Error: {exc}")

            # Provide helpful hints based on error type
            exc_lower = str(exc).lower()
            if "tolerance" in exc_lower or "golden" in exc_lower or "mismatch" in exc_lower:
                cons.print("    [dim]Hint: Consider --generate to update golden files or check tolerances[/dim]")
            elif "timeout" in exc_lower:
                cons.print("    [dim]Hint: Test may be hanging - check case configuration[/dim]")
            elif "nan" in exc_lower:
                cons.print("    [dim]Hint: NaN detected - check numerical stability of the case[/dim]")
            elif "failed to execute" in exc_lower:
                cons.print("    [dim]Hint: Check build logs and case parameters[/dim]")
            cons.print()

            # Track failed test details for summary
            error_type = ""
            exc_lower = str(exc).lower()
            if "tolerance" in exc_lower or "golden" in exc_lower or "mismatch" in exc_lower:
                error_type = "tolerance mismatch"
            elif "timeout" in exc_lower:
                error_type = "timeout"
            elif "nan" in exc_lower:
                error_type = "NaN detected"
            elif "failed to execute" in exc_lower:
                error_type = "execution failed"

            failed_tests.append({"trace": case.trace, "uuid": case.get_uuid(), "error_type": error_type, "attempts": nAttempts})

            # Still collect for final summary
            errors.append(f"[bold red]Failed test {case} after {nAttempts} attempt(s).[/bold red]")
            errors.append(f"{exc}")

        # Check if we should abort early due to high failure rate
        # Skip this check during dry-run (only builds, doesn't run tests)
        if not ARG("dry_run"):
            total_completed = nFAIL + nPASS
            if total_completed >= MIN_CASES_BEFORE_ABORT:
                failure_rate = nFAIL / total_completed
                if failure_rate >= FAILURE_RATE_THRESHOLD:
                    cons.print(f"\n[bold red]CRITICAL: {failure_rate * 100:.1f}% failure rate detected after {total_completed} tests.[/bold red]")
                    cons.print("[bold red]This suggests a systemic issue (bad build, broken environment, etc.)[/bold red]")
                    cons.print("[bold red]Aborting remaining tests to fail fast.[/bold red]\n")
                    # Set abort flag instead of raising exception from worker thread
                    abort_tests.set()
                    return  # Exit gracefully

        return
