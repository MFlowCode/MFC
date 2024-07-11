import os, math, typing, shutil, time
from random import sample

import rich, rich.table

from ..printer import cons
from ..        import common
from ..state   import ARG
from .case     import TestCase
from .cases    import list_cases
from ..        import sched
from ..run.input import MFCInputFile
from ..common  import MFCException, does_command_exist, format_list_to_string, get_program_output
from ..build   import build, HDF5, PRE_PROCESS, SIMULATION, POST_PROCESS

from ..packer import tol as packtol
from ..packer import packer


nFAIL = 0

def __filter(cases_) -> typing.List[TestCase]:
    cases = cases_[:]

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

    if ARG("percent") == 100:
        return cases

    return sample(cases, k=int(len(cases)*ARG("percent")/100.0))


def test():
    # pylint: disable=global-statement, global-variable-not-assigned
    global nFAIL

    cases = [ _.to_case() for _ in list_cases() ]

    # Delete UUIDs that are not in the list of cases from tests/
    if ARG("remove_old_tests"):
        dir_uuids = set(os.listdir(common.MFC_TESTDIR))
        new_uuids = { case.get_uuid() for case in cases }

        for old_uuid in dir_uuids - new_uuids:
            cons.print(f"[bold red]Deleting:[/bold red] {old_uuid}")
            common.delete_directory(f"{common.MFC_TESTDIR}/{old_uuid}")

        return

    cases = __filter(cases)

    if ARG("list"):
        table = rich.table.Table(title="MFC Test Cases", box=rich.table.box.SIMPLE)

        table.add_column("UUID", style="bold magenta", justify="center")
        table.add_column("Trace")

        for case in cases:
            table.add_row(case.get_uuid(), case.trace)

        rich.print(table)

        return

    codes = [PRE_PROCESS, SIMULATION] + ([POST_PROCESS] if ARG('test_all') else [])
    if not ARG("case_optimization"):
        build(codes)

    for case in cases:
        if case.rebuild:
            build(codes, MFCInputFile(os.path.basename(case.get_dirpath()), case.get_dirpath(), case.params))

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

    cons.print()
    if nFAIL == 0:
        cons.print("Tested Simulation [bold green]✓[/bold green]")
    else:
        raise MFCException(f"Testing: Encountered [bold red]{nFAIL}[/bold red] failure(s).")

    cons.unindent()


# pylint: disable=too-many-locals, too-many-branches, too-many-statements
def _handle_case(case: TestCase, devices: typing.Set[int]):
    start_time = time.time()

    tol = case.compute_tolerance()

    case.delete_output()
    case.create_directory()

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
            for pfilepath, pentry in pack.entries.items():
                if golden.find(pfilepath) is None:
                    golden.set(pentry)

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

        for t_step in [ i*case["t_step_save"] for i in range(0, math.floor(case["t_step_stop"] / case["t_step_save"]) + 1) ]:
            silo_filepath = os.path.join(case.get_dirpath(), 'silo_hdf5', 'p0', f'{t_step}.silo')
            if not os.path.exists(silo_filepath):
                silo_filepath = os.path.join(case.get_dirpath(), 'silo_hdf5', 'p_all', 'p0', f'{t_step}.silo')

            h5dump = f"{HDF5.get_install_dirpath(MFCInputFile(os.path.basename(case.get_filepath()), case.get_dirpath(), case.get_parameters()))}/bin/h5dump"

            if ARG("sys_hdf5"):
                if not does_command_exist("h5dump"):
                    raise MFCException("--sys-hdf5 was specified and h5dump couldn't be found.")

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
    # pylint: disable=global-statement
    global nFAIL

    nAttempts = 0

    while True:
        nAttempts += 1

        try:
            _handle_case(case, devices)
        except Exception as exc:
            if nAttempts < ARG("max_attempts"):
                cons.print(f"[bold yellow] Attempt {nAttempts}: Failed test {case.get_uuid()}. Retrying...[/bold yellow]")
                continue

            nFAIL += 1

            cons.print(f"[bold red]Failed test {case} after {nAttempts} attempt(s).[/bold red]")

            if ARG("relentless"):
                cons.print(f"{exc}")
            else:
                raise exc

        return
