import os, typing, shutil, time, itertools, threading
from random import sample, seed

import rich, rich.table
from rich.panel import Panel

from ..printer import cons
from ..        import common
from ..state   import ARG
from .case     import TestCase
from .cases    import list_cases
from ..        import sched
from ..common  import MFCException, does_command_exist, format_list_to_string, get_program_output
from ..build   import build, HDF5, PRE_PROCESS, SIMULATION, POST_PROCESS

from ..packer import tol as packtol
from ..packer import packer


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

def _filter_only(cases, skipped_cases):
    """Filter cases by --only terms using AND for labels, OR for UUIDs.

    Labels (non-UUID terms): case must match ALL labels (AND logic).
    UUIDs (8-char hex terms): case must match ANY UUID (OR logic).
    Mixed: keep case if all labels match OR any UUID matches.
    """
    def is_uuid(term):
        return len(term) == 8 and all(c in '0123456789abcdefABCDEF' for c in term)

    uuids  = [t for t in ARG("only") if is_uuid(t)]
    labels = [t for t in ARG("only") if not is_uuid(t)]

    for case in cases[:]:
        check = set(case.trace.split(" -> "))
        check.add(case.get_uuid())

        label_ok = all(label in check for label in labels) if labels else True
        uuid_ok  = any(u in check for u in uuids)  if uuids  else True

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


# pylint: disable=too-many-branches,too-many-locals,too-many-statements,trailing-whitespace
def __filter(cases_) -> typing.List[TestCase]:
    cases = cases_[:]
    selected_cases = []
    skipped_cases  = []

    # Check "--from" and "--to" exist and are in the right order
    bFoundFrom, bFoundTo = (False, False)
    from_i = -1
    for i, case in enumerate(cases):
        if case.get_uuid() == ARG("from"):
            from_i     = i
            bFoundFrom = True
            # Do not "continue" because "--to" might be the same as "--from"
        if bFoundFrom and case.get_uuid() == ARG("to"):
            cases    = cases[from_i:i+1]
            skipped_cases = [case for case in cases_ if case not in cases]
            bFoundTo = True
            break

    if not bFoundTo:
        raise MFCException("Testing: Your specified range [--from,--to] is incorrect. Please ensure both IDs exist and are in the correct order.")

    if len(ARG("only")) > 0:
        cases, skipped_cases = _filter_only(cases, skipped_cases)

        if not cases:
            raise MFCException(
                f"--only filter matched zero test cases. "
                f"Specified: {ARG('only')}. Check that UUIDs/names are valid."
            )

    # --only-changes: filter based on file-level gcov coverage
    if ARG("only_changes"):
        from .coverage import (  # pylint: disable=import-outside-toplevel
            load_coverage_cache, get_changed_files,
            should_run_all_tests, filter_tests_by_coverage,
        )

        cache = load_coverage_cache(common.MFC_ROOT_DIR)
        if cache is None:
            cons.print("[yellow]Coverage cache missing or stale.[/yellow]")
            cons.print("[yellow]Run: ./mfc.sh build --gcov -j 8 && ./mfc.sh test --build-coverage-cache[/yellow]")
            cons.print("[yellow]Falling back to full test suite.[/yellow]")
        else:
            changed_files = get_changed_files(common.MFC_ROOT_DIR, ARG("changes_branch"))

            if changed_files is None:
                cons.print("[yellow]git diff failed — falling back to full test suite.[/yellow]")
            elif should_run_all_tests(changed_files):
                cons.print()
                cons.print("[bold cyan]Coverage Change Analysis[/bold cyan]")
                cons.print("-" * 50)
                cons.print("[yellow]Infrastructure or macro file changed — running full test suite.[/yellow]")
                cons.print("-" * 50)
            else:
                changed_fpp = {f for f in changed_files if f.endswith(".fpp")}
                if not changed_fpp:
                    cons.print()
                    cons.print("[bold cyan]Coverage Change Analysis[/bold cyan]")
                    cons.print("-" * 50)
                    cons.print("[green]No .fpp source changes detected — skipping all tests.[/green]")
                    cons.print("-" * 50)
                    cons.print()
                    skipped_cases += cases
                    cases = []
                else:
                    cons.print()
                    cons.print("[bold cyan]Coverage Change Analysis[/bold cyan]")
                    cons.print("-" * 50)
                    for fpp_file in sorted(changed_fpp):
                        cons.print(f"  [green]*[/green] {fpp_file}")

                    cases, new_skipped = filter_tests_by_coverage(cases, cache, changed_files)
                    skipped_cases += new_skipped
                    cons.print(f"\n[bold]Tests to run: {len(cases)} / {len(cases) + len(new_skipped)}[/bold]")
                    cons.print("-" * 50)
                    cons.print()

    for case in cases[:]:
        if case.ppn > 1 and not ARG("mpi"):
            cases.remove(case)
            skipped_cases.append(case)

    for case in cases[:]:
        if ARG("single"):
            skip = ['low_Mach', 'Hypoelasticity', 'teno', 'Chemistry', 'Phase Change model 6'
            ,'Axisymmetric', 'Transducer', 'Transducer Array', 'Cylindrical', 'HLLD', 'Example']
            if any(label in case.trace for label in skip):
                cases.remove(case)
                skipped_cases.append(case)

    for case in cases[:]:
        if ARG("gpu"):
            skip = ['Gauss Seidel']
            if any(label in case.trace for label in skip):
                cases.remove(case)

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
            raise MFCException(
                f"--shard {ARG('shard')} matched zero test cases. "
                f"Total cases before sharding may be less than shard count."
            )

    if ARG("percent") == 100:
        return cases, skipped_cases

    seed(time.time())

    selected_cases = sample(cases, k=int(len(cases)*ARG("percent")/100.0))
    skipped_cases += [item for item in cases if item not in selected_cases]

    return selected_cases, skipped_cases

def test():
    # pylint: disable=global-statement, global-variable-not-assigned, too-many-statements, too-many-locals
    global nFAIL, nPASS, nSKIP, total_test_count
    global errors, failed_tests, test_start_time

    test_start_time = time.time()  # Start timing
    failed_uuids_path = os.path.join(common.MFC_TEST_DIR, "failed_uuids.txt")
    cases = list_cases()

    # Delete UUIDs that are not in the list of cases from tests/
    if ARG("remove_old_tests"):
        dir_uuids = set(os.listdir(common.MFC_TEST_DIR))
        new_uuids = { case.get_uuid() for case in cases }

        for old_uuid in dir_uuids - new_uuids:
            cons.print(f"[bold red]Deleting:[/bold red] {old_uuid}")
            common.delete_directory(f"{common.MFC_TEST_DIR}/{old_uuid}")

        return

    if ARG("build_coverage_cache"):
        from .coverage import build_coverage_cache  # pylint: disable=import-outside-toplevel
        all_cases = [b.to_case() for b in cases]

        # Build all unique slugs (Chemistry, case-optimization, etc.) so every
        # test has a compatible binary when run with --no-build.
        codes = [PRE_PROCESS, SIMULATION, POST_PROCESS]
        unique_builds = set()
        for case, code in itertools.product(all_cases, codes):
            slug = code.get_slug(case.to_input_file())
            if slug not in unique_builds:
                build(code, case.to_input_file())
                unique_builds.add(slug)

        build_coverage_cache(common.MFC_ROOT_DIR, all_cases,
                             extra_args=ARG("--"), n_jobs=int(ARG("jobs")))
        return

    cases, skipped_cases = __filter(cases)
    cases = [ _.to_case() for _ in cases ]
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
    codes = [PRE_PROCESS, SIMULATION] + ([POST_PROCESS] if ARG('test_all') else [])
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

    cons.print(f"[bold]Test {format_list_to_string([ x.name for x in codes ], 'magenta')}[/bold] | {range_str} ({len(cases)} test{'s' if len(cases) != 1 else ''})")
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
    sched.sched(
        [ sched.Task(ppn=case.ppn, func=handle_case, args=[case], load=case.get_cell_count()) for case in cases ],
        ARG("jobs"), ARG("gpus"))

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
            raise MFCException(
                f"Excessive test failures: {nFAIL}/{total_completed} "
                f"failed ({nFAIL/total_completed*100:.1f}%)"
            )
        raise MFCException(
            f"Excessive test failures: {nFAIL} failed, but no tests were completed."
        )

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
                f.write(test_info['uuid'] + "\n")
    elif os.path.exists(failed_uuids_path):
        os.remove(failed_uuids_path)

    exit(nFAIL)


def _print_test_summary(passed: int, failed: int, skipped: int, minutes: int, seconds: float,
                        failed_test_list: list, _skipped_cases: list):
    # pylint: disable=too-many-arguments, too-many-positional-arguments, too-many-locals
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
            trace = test_info.get('trace', 'Unknown')
            uuid = test_info.get('uuid', 'Unknown')
            error_type = test_info.get('error_type', '')
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
    cons.raw.print(Panel(
        "\n".join(summary_lines),
        title="[bold]Test Summary[/bold]",
        border_style=border_style,
        padding=(1, 2)
    ))
    cons.print()


# pylint: disable=too-many-locals, too-many-branches, too-many-statements, trailing-whitespace
def _process_silo_file(silo_filepath: str, case: TestCase, out_filepath: str):
    """Process a single silo file with h5dump and check for NaNs/Infinities."""
    h5dump = f"{HDF5.get_install_dirpath(case.to_input_file())}/bin/h5dump"

    if not os.path.exists(h5dump or ""):
        if not does_command_exist("h5dump"):
            raise MFCException("h5dump couldn't be found.")
        h5dump = shutil.which("h5dump")

    output, err = get_program_output([h5dump, silo_filepath])

    if err != 0:
        raise MFCException(
            f"Test {case}: Failed to run h5dump. You can find the run's output in {out_filepath}, "
            f"and the case dictionary in {case.get_filepath()}."
        )

    if "nan," in output:
        raise MFCException(
            f"Test {case}: Post Process has detected a NaN. You can find the run's output in {out_filepath}, "
            f"and the case dictionary in {case.get_filepath()}."
        )

    if "inf," in output:
        raise MFCException(
            f"Test {case}: Post Process has detected an Infinity. You can find the run's output in {out_filepath}, "
            f"and the case dictionary in {case.get_filepath()}."
        )


def _handle_case(case: TestCase, devices: typing.Set[int]):
    # pylint: disable=global-statement, global-variable-not-assigned
    global current_test_number
    start_time = time.time()

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

        pack, err = packer.pack(case.get_dirpath())
        if err is not None:
            raise MFCException(f"Test {case}: {err}")

        if pack.has_NaNs():
            raise MFCException(f"Test {case}: NaNs detected in the case.")

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

        if ARG("test_all"):
            case.delete_output()
            # Check timeout before launching the (potentially long) post-process run
            if timeout_flag.is_set():
                raise TestTimeoutError("Test case exceeded 1 hour timeout")
            cmd = case.run([PRE_PROCESS, SIMULATION, POST_PROCESS], gpus=devices)
            out_filepath = os.path.join(case.get_dirpath(), "out_post.txt")
            common.file_write(out_filepath, cmd.stdout)

            silo_dir = os.path.join(case.get_dirpath(), 'silo_hdf5', 'p0')
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
        log_path = os.path.join(case.get_dirpath(), 'out_pre_sim.txt')
        if os.path.exists(log_path):
            log_msg = f"Check the log at: {log_path}"
        else:
            log_msg = (
                f"Log file ({log_path}) may not exist if the timeout occurred early."
            )
        raise MFCException(
            f"Test {case} exceeded 1 hour timeout. "
            f"This may indicate a hung simulation or misconfigured case. "
            f"{log_msg}"
        ) from exc
    finally:
        timeout_timer.cancel()  # Cancel timeout timer


def handle_case(case: TestCase, devices: typing.Set[int]):
    # pylint: disable=global-statement, global-variable-not-assigned
    global nFAIL, nPASS, nSKIP
    global errors, failed_tests

    # Check if we should abort before processing this case
    if abort_tests.is_set():
        return  # Exit gracefully if abort was requested

    nAttempts = 0
    if ARG('single'):
        max_attempts = max(ARG('max_attempts'), 3)
    else:
        max_attempts = ARG('max_attempts')

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

            # Show truncated error message
            exc_str = str(exc)
            if len(exc_str) > 300:
                exc_str = exc_str[:297] + "..."
            cons.print(f"    Error: {exc_str}")

            # Provide helpful hints based on error type
            exc_lower = str(exc).lower()
            if "tolerance" in exc_lower or "golden" in exc_lower or "mismatch" in exc_lower:
                cons.print(f"    [dim]Hint: Consider --generate to update golden files or check tolerances[/dim]")
            elif "timeout" in exc_lower:
                cons.print(f"    [dim]Hint: Test may be hanging - check case configuration[/dim]")
            elif "nan" in exc_lower:
                cons.print(f"    [dim]Hint: NaN detected - check numerical stability of the case[/dim]")
            elif "failed to execute" in exc_lower:
                cons.print(f"    [dim]Hint: Check build logs and case parameters[/dim]")
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

            failed_tests.append({
                'trace': case.trace,
                'uuid': case.get_uuid(),
                'error_type': error_type,
                'attempts': nAttempts
            })

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
                    cons.print(f"\n[bold red]CRITICAL: {failure_rate*100:.1f}% failure rate detected after {total_completed} tests.[/bold red]")
                    cons.print("[bold red]This suggests a systemic issue (bad build, broken environment, etc.)[/bold red]")
                    cons.print("[bold red]Aborting remaining tests to fail fast.[/bold red]\n")
                    # Set abort flag instead of raising exception from worker thread
                    abort_tests.set()
                    return  # Exit gracefully

        return
