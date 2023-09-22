import os, math, typing, shutil, time

from random    import sample
from ..printer import cons
from ..        import common
from ..state   import ARG
from .case     import TestCase
from .cases    import generate_cases
from ..        import sched
from ..common  import MFCException, does_command_exist, format_list_to_string, get_program_output
from ..build   import build_targets, HDF5, PRE_PROCESS, SIMULATION, POST_PROCESS

from ..packer import tol as packtol
from ..packer import packer

import rich, rich.table


CASES = generate_cases()
nFAIL = 0

def __filter():
    global CASES
    
    # Check "--from" and "--to" exist and are in the right order
    bFoundFrom, bFoundTo = (False, False)
    from_i = -1
    for i, case in enumerate(CASES):
        if case.get_uuid() == ARG("from"):
            from_i     = i
            bFoundFrom = True
            # Do not "continue" because "--to" might be the same as "--from"
        if bFoundFrom and case.get_uuid() == ARG("to"):
            CASES    = CASES[from_i:i+1]
            bFoundTo = True
            break

    if not bFoundTo:
        raise MFCException("Testing: Your specified range [--from,--to] is incorrect. Please ensure both IDs exist and are in the correct order.")

    if len(ARG("only")) > 0:
        for i, case in enumerate(CASES[:]):
            case: TestCase

            doKeep = False
            for o in ARG("only"):
                if str(o) == case.get_uuid():
                    doKeep = True
                    break

            if not doKeep:
                CASES.remove(case)

    if not ARG("mpi"):
        for case in CASES[:]:
            if case.ppn > 1:
                CASES.remove(case)

    if ARG("percent") == 100:
        return

    CASES = sample(CASES, k=int(len(CASES)*ARG("percent")/100.0))


def test():
    global CASES, nFAIL
   
    # Delete UUIDs that are not in the list of CASES from tests/
    if ARG("remove_old_tests"):
        dir_uuids = set(os.listdir(common.MFC_TESTDIR))
        new_uuids = set([case.get_uuid() for case in CASES])

        for old_uuid in dir_uuids - new_uuids:
            cons.print(f"[bold red]Deleting:[/bold red] {old_uuid}")
            common.delete_directory(f"{common.MFC_TESTDIR}/{old_uuid}")
        
        return

    __filter()

    if ARG("list"):
        table = rich.table.Table(title="MFC Test Cases", box=rich.table.box.SIMPLE)

        table.add_column("UUID", style="bold magenta", justify="center")
        table.add_column("Trace")

        for case in CASES:
            table.add_row(case.get_uuid(), case.trace)

        rich.print(table)

        return

    codes = [PRE_PROCESS, SIMULATION] + ([POST_PROCESS] if ARG('test_all') else [])
    if not ARG("case_optimization"):
        build_targets(codes)

    cons.print()

    range_str = f"from [bold magenta]{ARG('from')}[/bold magenta] to [bold magenta]{ARG('to')}[/bold magenta]"

    if len(ARG("only")) > 0:
        range_str = "Only " + format_list_to_string(ARG("only"), "bold magenta", "Nothing to run")
    
    cons.print(f"[bold]Test {format_list_to_string([ x.name for x in codes ], 'magenta')}[/bold] | {range_str} ({len(CASES)} test{'s' if len(CASES) != 1 else ''})")
    cons.indent()

    # Run CASES with multiple threads (if available)
    cons.print()
    cons.print(f" tests/[bold magenta]UUID[/bold magenta]     (s)      Summary")
    cons.print()
    
    # Select the correct number of threads to use to launch test CASES
    # We can't use ARG("jobs") when the --case-optimization option is set
    # because running a test case may cause it to rebuild, and thus
    # interfere with the other test CASES. It is a niche feature so we won't
    # engineer around this issue (for now).
    sched.sched(
        [ sched.Task(ppn=case.ppn, func=handle_case, args=[case], load=case.get_cell_count()) for case in CASES ],
        ARG("jobs"), ARG("gpus"))

    cons.print()
    if nFAIL == 0:
        cons.print(f"Tested Simulation [bold green]✓[/bold green]")
    else:
        if nFAIL == 1:
            raise MFCException(f"Testing: There was [bold red]1[/bold red] failure.")
        else:
            raise MFCException(f"Testing: There were [bold red]{nFAIL}[/bold red] failures.")

    cons.unindent()


def _handle_case(test: TestCase, devices: typing.Set[int]):
    start_time = time.time()

    if test.params.get("qbmm", 'F') == 'T':
        tol = 1e-10
    elif test.params.get("bubbles", 'F') == 'T':
        tol = 1e-10
    elif test.params.get("hypoelasticity", 'F') == 'T':
        tol = 1e-7 
    elif test.params.get("relax", 'F') == 'T':
        tol = 5e-7 
    else:
        tol = 1e-12

    test.delete_output()
    test.create_directory()

    cmd = test.run([PRE_PROCESS, SIMULATION], gpus=devices)

    out_filepath = os.path.join(test.get_dirpath(), "out_pre_sim.txt")

    common.file_write(out_filepath, cmd.stdout)

    if cmd.returncode != 0:
        cons.print(cmd.stdout)
        raise MFCException(f"Test {test}: Failed to execute MFC.")

    pack, err = packer.pack(test.get_dirpath())
    if err is not None:
        raise MFCException(f"Test {test}: {err}")

    if pack.has_NaNs():
        raise MFCException(f"Test {test}: NaNs detected in the case.")

    golden_filepath = os.path.join(test.get_dirpath(), "golden.txt")
    golden_meta_filepath = os.path.join(test.get_dirpath(), "golden-metadata.txt")
    if ARG("generate"):
        common.delete_file(golden_filepath)
        pack.save(golden_filepath)
    else:
        if not os.path.isfile(golden_filepath):
            raise MFCException(f"Test {test}: The golden file does not exist! To generate golden files, use the '--generate' flag.")

        golden = packer.load(golden_filepath)

        if ARG("add_new_variables"):
            for pfilepath, pentry in pack.entries.items():
                if golden.find(pfilepath) is None:
                    golden.set(pentry)

            golden.save(golden_filepath)
        else:
            err, msg = packtol.compare(pack, packer.load(golden_filepath), packtol.Tolerance(tol, tol))
            if msg is not None:
                raise MFCException(f"Test {test}: {msg}")

    if ARG("test_all"):
        test.delete_output()
        cmd = test.run([PRE_PROCESS, SIMULATION, POST_PROCESS], gpus=devices)
        out_filepath = os.path.join(test.get_dirpath(), "out_post.txt")
        common.file_write(out_filepath, cmd.stdout)

        for t_step in [ i*test["t_step_save"] for i in range(0, math.floor(test["t_step_stop"] / test["t_step_save"]) + 1) ]:
            silo_filepath = os.path.join(test.get_dirpath(), 'silo_hdf5', 'p0', f'{t_step}.silo')
            if not os.path.exists(silo_filepath):
                silo_filepath = os.path.join(test.get_dirpath(), 'silo_hdf5', 'p_all', 'p0', f'{t_step}.silo')
    
            h5dump = f"{HDF5.get_install_dirpath()}/bin/h5dump"

            if ARG("no_hdf5"):
                if not does_command_exist("h5dump"):
                    raise MFCException("--no-hdf5 was specified and h5dump couldn't be found.")
                
                h5dump = shutil.which("h5dump")

            output, err = get_program_output([h5dump, silo_filepath])

            if err != 0:
                raise MFCException(f"""Test {test}: Failed to run h5dump. You can find the run's output in {out_filepath}, and the case dictionary in {os.path.join(test.get_dirpath(), "case.py")}.""")

            if "nan," in output:
                raise MFCException(f"""Test {test}: Post Process has detected a NaN. You can find the run's output in {out_filepath}, and the case dictionary in {os.path.join(test.get_dirpath(), "case.py")}.""")

            if "inf," in output:
                raise MFCException(f"""Test {test}: Post Process has detected an Infinity. You can find the run's output in {out_filepath}, and the case dictionary in {os.path.join(test.get_dirpath(), "case.py")}.""")

    test.delete_output()

    end_time = time.time()
    duration = end_time - start_time

    cons.print(f"  [bold magenta]{test.get_uuid()}[/bold magenta]    {duration:6.2f}    {test.trace}")


def handle_case(test: TestCase, devices: typing.Set[int]):
    global nFAIL
    
    nAttempts = 0

    while True:
        nAttempts += 1

        try:
            _handle_case(test, devices)
        except Exception as exc:
            if nAttempts < ARG("max_attempts"):
                cons.print(f"[bold yellow] Attempt {nAttempts}: Failed test {test.get_uuid()}. Retrying...[/bold yellow]")
                continue

            nFAIL += 1

            cons.print(f"[bold red]Failed test {test} after {nAttempts} attempt(s).[/bold red]")

            if ARG("relentless"):
                cons.print(f"{exc}")
            else:
                raise exc

        return
