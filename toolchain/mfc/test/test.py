import os

from ..printer import cons
from ..        import common
from ..state   import ARG
from .case     import TestCase
from .cases    import generate_cases
from ..        import sched
from ..common  import MFCException
from ..build   import build_targets
from .         import pack as packer

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


def test():
    global CASES, nFAIL
    
    # Delete UUIDs that are not in the list of CASES from tests/
    if ARG("generate"):
        dir_uuids = set([name for name in os.listdir(".") if os.path.isdir(name)])
        new_uuids = set([case.get_uuid() for case in CASES])

        for old_uuid in dir_uuids - new_uuids:
            common.delete_directory(f"{common.MFC_TESTDIR}/{old_uuid}")

    __filter()

    if ARG("list"):
        table = rich.table.Table(title="MFC Test Cases", box=rich.table.box.SIMPLE)

        table.add_column("UUID", style="bold magenta", justify="center")
        table.add_column("Trace")

        for case in CASES:
            table.add_row(case.get_uuid(), case.trace)

        rich.print(table)

        return

    build_targets(["pre_process", "simulation"])

    range_str = f"from [bold magenta]{ARG('from')}[/bold magenta] to [bold magenta]{ARG('to')}[/bold magenta]"

    if len(ARG("only")) > 0:
        range_str = "Only " + common.format_list_to_string([
            f"[bold magenta]{uuid}[/bold magenta]" for uuid in ARG("only")
        ], "Nothing to run")


    cons.print(f"[bold]Test[/bold] | {range_str} ({len(CASES)} test{'s' if len(CASES) != 1 else ''})")
    cons.indent()


    # Run CASES with multiple threads (if available)
    cons.print()
    cons.print(f" tests/[bold magenta]UUID[/bold magenta]    Summary")
    cons.print()
    
    # Select the correct number of threads to use to launch test CASES
    # We can't use ARG("jobs") when the --case-optimization option is set
    # because running a test case may cause it to rebuild, and thus
    # interfere with the other test CASES. It is a niche feature so we won't
    # engineer around this issue (for now).
    nThreads = ARG("jobs") if not ARG("case_optimization") else 1
    tasks    = [
        sched.Task(ppn=case.ppn, func=handle_case, args=[ case ]) for case in CASES
    ]
    sched.sched(tasks, nThreads)
    
    cons.print()
    if nFAIL == 0:
        cons.print(f"Tested [bold green]✓[/bold green]")
        cons.unindent()
    else:
        if nFAIL == 1:
            raise MFCException(f"Testing: There was [bold red]1[/bold red] failure.")
        else:
            raise MFCException(f"Testing: There were [bold red]{nFAIL}[/bold red] failures.")


def handle_case(test: TestCase):
    global nFAIL
    
    try:
        test.create_directory()

        if test.params.get("qbmm", 'F') == 'T':
            tol = 1e-10
        elif test.params.get("bubbles", 'F') == 'T':
            tol = 1e-10
        elif test.params.get("hypoelasticity", 'F') == 'T':
            tol = 1e-7
        else:
            tol = 1e-12

        cmd = test.run()

        out_filepath = os.path.join(test.get_dirpath(), "out.txt")

        common.file_write(out_filepath, cmd.stdout)

        if cmd.returncode != 0:
            cons.print(cmd.stdout)
            raise MFCException(f"""Test {test}: Failed to execute MFC. You can find the run's output in {out_filepath}, and the case dictionary in {os.path.join(test.get_dirpath(), "case.py")}.""")

        pack = packer.generate(test)
        pack.save(os.path.join(test.get_dirpath(), "pack.txt"))

        golden_filepath = os.path.join(test.get_dirpath(), "golden.txt")
        if ARG("generate"):
            common.delete_file(golden_filepath)
            pack.save(golden_filepath)
        else:
            if not os.path.isfile(golden_filepath):
                raise MFCException(f"Test {test}: Golden file doesn't exist! To generate golden files, use the '-g' flag.")

            packer.check_tolerance(test, pack, packer.load(golden_filepath), tol)

        cons.print(f"  [bold magenta]{test.get_uuid()}[/bold magenta]    {test.trace}")
    except Exception as exc:
        nFAIL = nFAIL + 1

        if not ARG("relentless"):
            raise exc

        cons.print(f"[bold red]Failed test {test}.[/bold red]")
        cons.print(f"{exc}")
