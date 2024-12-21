import os, typing, shutil, time, itertools
from random import sample, seed

import rich, rich.table

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
errors = []

# pylint: disable=too-many-branches, trailing-whitespace
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
            bFoundTo = True
            break

    if not bFoundTo:
        raise MFCException("Testing: Your specified range [--from,--to] is incorrect. Please ensure both IDs exist and are in the correct order.")

    if len(ARG("only")) > 0:
        for case in cases[:]:
            case: TestCase

            checkCase = case.trace.split(" -> ")
            checkCase.append(case.get_uuid())
            if not set(ARG("only")).issubset(set(checkCase)):
                cases.remove(case)

    for case in cases[:]:
        if case.ppn > 1 and not ARG("mpi"):
            cases.remove(case)
            skipped_cases.append(case)
    
    for case in cases[:]:
        if ARG("single"):
            skip = ['low_Mach', 'Hypoelasticity', 'teno', 'Chemistry', 'Phase Change model 6'
            ,'Axisymmetric', 'Transducer', 'Transducer Array', 'Cylindrical', 'Example']
            if any(label in case.trace for label in skip):
                cases.remove(case)


    if ARG("no_examples"):
        cases = [case for case in cases if not "Example" in case.trace]

    if ARG("percent") == 100:
        return cases, skipped_cases

    seed(time.time())

    selected_cases = sample(cases, k=int(len(cases)*ARG("percent")/100.0))
    skipped_cases = [item for item in cases if item not in selected_cases]

    return selected_cases, skipped_cases

def test():
    # pylint: disable=global-statement, global-variable-not-assigned
    global nFAIL, nPASS, nSKIP
    global errors

    cases = list_cases()

    # Delete UUIDs that are not in the list of cases from tests/
    if ARG("remove_old_tests"):
        dir_uuids = set(os.listdir(common.MFC_TEST_DIR))
        new_uuids = { case.get_uuid() for case in cases }

        for old_uuid in dir_uuids - new_uuids:
            cons.print(f"[bold red]Deleting:[/bold red] {old_uuid}")
            common.delete_directory(f"{common.MFC_TEST_DIR}/{old_uuid}")

        return

    cases, skipped_cases = __filter(cases)
    cases = [ _.to_case() for _ in cases ]

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
    cons.print(" tests/[bold magenta]UUID[/bold magenta]     (s)      Summary")
    cons.print()

    # Select the correct number of threads to use to launch test cases
    # We can't use ARG("jobs") when the --case-optimization option is set
    # because running a test case may cause it to rebuild, and thus
    # interfere with the other test cases. It is a niche feature so we won't
    # engineer around this issue (for now).
    sched.sched(
        [ sched.Task(ppn=case.ppn, func=handle_case, args=[case], load=case.get_cell_count()) for case in cases ],
        ARG("jobs"), ARG("gpus"))

    nSKIP = len(skipped_cases)
    cons.print()
    cons.unindent()
    cons.print(f"\nTest Summary: [bold green]{nPASS}[/bold green] passed, [bold red]{nFAIL}[/bold red] failed, [bold yellow]{nSKIP}[/bold yellow] skipped.\n")

    # Print a summary of all errors at the end if errors exist
    if len(errors) != 0:
        cons.print(f"[bold red]Failed Cases[/bold red]\n")
        for e in errors:
            cons.print(e)

    # Print the list of skipped cases
    if len(skipped_cases) != 0:
        cons.print("[bold yellow]Skipped Cases[/bold yellow]\n")
        for c in skipped_cases:
            cons.print(f"[bold yellow]{c.trace}[/bold yellow]")

    exit(nFAIL)


# pylint: disable=too-many-locals, too-many-branches, too-many-statements, trailing-whitespace
def _handle_case(case: TestCase, devices: typing.Set[int]):
    # pylint: disable=global-statement, global-variable-not-assigned
    start_time = time.time()

    tol = case.compute_tolerance()
    case.delete_output()
    case.create_directory()

    if ARG("dry_run"):
        cons.print(f"  [bold magenta]{case.get_uuid()}[/bold magenta]     SKIP     {case.trace}")
        return

    cmd = case.run([PRE_PROCESS, SIMULATION], gpus=devices)
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
        cmd = case.run([PRE_PROCESS, SIMULATION, POST_PROCESS], gpus=devices)
        out_filepath = os.path.join(case.get_dirpath(), "out_post.txt")
        common.file_write(out_filepath, cmd.stdout)

        for silo_filepath in os.listdir(os.path.join(case.get_dirpath(), 'silo_hdf5', 'p0')):
            silo_filepath = os.path.join(case.get_dirpath(), 'silo_hdf5', 'p0', silo_filepath)
            h5dump        = f"{HDF5.get_install_dirpath(case.to_input_file())}/bin/h5dump"

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

    case.delete_output()

    end_time = time.time()
    duration = end_time - start_time

    cons.print(f"  [bold magenta]{case.get_uuid()}[/bold magenta]    {duration:6.2f}    {case.trace}")


def handle_case(case: TestCase, devices: typing.Set[int]):
    # pylint: disable=global-statement, global-variable-not-assigned
    global nFAIL, nPASS, nSKIP
    global errors

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
            cons.print(f"[bold red]Failed test {case} after {nAttempts} attempt(s).[/bold red]")
            errors.append(f"[bold red]Failed test {case} after {nAttempts} attempt(s).[/bold red]")
            errors.append(f"{exc}")

        return
