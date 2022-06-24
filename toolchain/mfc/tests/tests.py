#!/usr/bin/env python3

import os


import common

from tests.case    import Case
from tests.cases   import generate_cases
from tests.threads import MFCTestThreadManager

from common import MFCException

import tests.pack


import rich
import rich.table


class MFCTest:
    def __init__(self, mfc):
        self.mfc    = mfc
        self.sched  = MFCTestThreadManager(1)
        self.cases  = generate_cases()

    def __filter_tests(self):
        # Check "--from" and "--to" exist and are in the right order
        bFoundFrom, bFoundTo = (False, False)
        from_i = -1
        for i, case in enumerate(self.cases):
            if case.get_uuid() == self.mfc.args["from"]:
                from_i     = i
                bFoundFrom = True
                # Do not "continue" because "--to" might be the same as "--from"
            if bFoundFrom and case.get_uuid() == self.mfc.args["to"]:
                self.cases = self.cases[from_i:i+1]
                bFoundTo = True
                break
        
        if not bFoundTo:
            raise MFCException("Testing: Your specified range [--from,--to] is incorrect. Please ensure both IDs exist and are in the correct order.")

        if len(self.mfc.args["only"]) > 0:
            for i, case in enumerate(self.cases[:]):
                case: Case

                doKeep = False
                for o in self.mfc.args["only"]:
                    if str(o) == case.get_uuid():
                        doKeep = True
                        break

                if not doKeep:
                    self.cases.remove(case)


    def execute(self):
        self.__filter_tests()

        if self.mfc.args["list"]:
            table = rich.table.Table(title="MFC Test Cases", box=rich.table.box.SIMPLE)

            table.add_column("UUID", style="magenta", justify="center")
            table.add_column("Trace")

            for case in self.cases:
                table.add_row(case.get_uuid(), case.trace)
                
            rich.print(table)

            return

        range_str = f"from [bold magenta]{self.mfc.args['from']}[/bold magenta] to [bold magenta]{self.mfc.args['to']}[/bold magenta]"

        if len(self.mfc.args["only"]) > 0:
            range_str = common.format_list_to_string([
                f"[bold magenta]{uuid}[/bold magenta]" for uuid in self.mfc.args["only"]
            ], "Nothing to run")

        rich.print(f"[bold]Test[/bold] | {range_str} ({len(self.cases)} tests)")

        # Clear previous tests if we wish to (re)generate golden files
        if self.mfc.args["generate"]:
            common.delete_directory_recursive(common.MFC_TESTDIR)
            common.create_directory(common.MFC_TESTDIR)

        # Run cases with multiple threads (if available)
        rich.print(f"")
        rich.print(f"  tests/[bold magenta]UUID[/bold magenta]   Error RE     Tol.     Summary")
        rich.print(f"")
        self.sched.run(self.cases, self.handle_case)

        rich.print(f"> Tested [bold green]✓[/bold green]")

    def handle_case(self, test: Case):
        test.create_directory()

        if test.case.bubbles.qbmm:
            tol = 1e-7
        elif test.case.bubbles.bubbles:
            tol = 1e-10
        else:
            tol = 1e-12

        cmd = test.run(self.mfc.args)

        common.file_write(f"{test.get_dirpath()}/out.txt", cmd.stdout)

        if cmd.returncode != 0:
            rich.print(cmd.stdout)
            raise MFCException(f"""\
tests/{test.get_uuid()}: Failed to execute MFC [{test.trace}]. Above is the output of MFC.
You can find the output in {test.get_dirpath()}/out.txt, and teh case dictionary in {test.get_dirpath()}/case.py.""")

        pack = tests.pack.generate(test)
        pack.save(f"{test.get_dirpath()}/pack.txt")

        golden_filepath = f"{test.get_dirpath()}/golden.txt"

        if self.mfc.args["generate"]:
            common.delete_file(golden_filepath)
            pack.save(golden_filepath)

        if not os.path.isfile(golden_filepath):
            raise MFCException(f"tests/{test.get_uuid()}: Golden file doesn't exist! To generate golden files, use the '-g' flag. [{test.trace}]")

        error = tests.pack.check_tolerance(test, pack, tests.pack.load(golden_filepath), tol)
        rich.print(f"  [bold magenta]{test.get_uuid()}[/bold magenta]    {error.relative:+0.1E}   {tol:+0.1E}   {test.trace}")
