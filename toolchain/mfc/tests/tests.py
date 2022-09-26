#!/usr/bin/env python3

import os

import rich

from ..util.printer import cons

from ..util import common

from .case    import Case
from .cases   import generate_cases
from .threads import MFCTestThreadManager

from ..util.common import MFCException

from ..build import build_target

from . import pack as packer


import rich
import rich.table


class MFCTest:
    def __init__(self, mfc):
        self.mfc    = mfc
        self.sched  = MFCTestThreadManager(1)
        self.cases  = generate_cases()
        self.nFail  = 0

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
        # Select the correct number of threads to use to launch test cases
        # We can't use args["jobs"] when the --case-optimization option is set
        # because running a test case may cause it to rebuild, and thus
        # interfere with the other test cases. It is a niche feature so we won't
        # engineer around this issue (for now).
        self.sched.nAvailable = self.mfc.args["jobs"] if not self.mfc.args["case_optimization"] else 1

        # Delete UUIDs that are not in the list of cases from tests/
        if self.mfc.args["generate"]:
            dir_uuids = set([name for name in os.listdir(".") if os.path.isdir(name)])
            new_uuids = set([case.get_uuid() for case in self.cases])

            for old_uuid in dir_uuids - new_uuids:
                common.delete_directory(f"{common.MFC_TESTDIR}/{old_uuid}")

        self.__filter_tests()

        if self.mfc.args["list"]:
            table = rich.table.Table(title="MFC Test Cases", box=rich.table.box.SIMPLE)

            table.add_column("UUID", style="bold magenta", justify="center")
            table.add_column("Trace")

            for case in self.cases:
                table.add_row(case.get_uuid(), case.trace)

            rich.print(table)

            return

        build_target(self.mfc, "pre_process")
        build_target(self.mfc, "simulation")

        range_str = f"from [bold magenta]{self.mfc.args['from']}[/bold magenta] to [bold magenta]{self.mfc.args['to']}[/bold magenta]"

        if len(self.mfc.args["only"]) > 0:
            range_str = "Only " + common.format_list_to_string([
                f"[bold magenta]{uuid}[/bold magenta]" for uuid in self.mfc.args["only"]
            ], "Nothing to run")


        cons.print(f"[bold]Test[/bold] | {range_str} ({len(self.cases)} test{'s' if len(self.cases) != 1 else ''})")
        cons.indent()


        # Run cases with multiple threads (if available)
        cons.print()
        cons.print(f" tests/[bold magenta]UUID[/bold magenta]    Summary")
        cons.print()
        self.sched.run(self.cases, self.handle_case)

        cons.print()
        if self.nFail == 0:
            cons.print(f"Tested [bold green]✓[/bold green]")
            cons.unindent()
        else:
            if self.nFail == 1:
                raise MFCException(f"Testing: There was [bold red]1[/bold red] failure.")
            else:
                raise MFCException(f"Testing: There were [bold red]{self.nFail}[/bold red] failures.")


    def handle_case(self, test: Case):
        try:
            test.create_directory()

            if test.params.get("qbmm", 'F') == 'T':
                tol = 1e-10
            elif test.params.get("bubbles", 'F') == 'T':
                tol = 1e-10
            else:
                tol = 1e-12

            cmd = test.run(self.mfc.args)

            out_filepath = os.path.join(test.get_dirpath(), "out.txt")

            common.file_write(out_filepath, cmd.stdout)

            if cmd.returncode != 0:
                cons.print(cmd.stdout)
                raise MFCException(f"""Test {test}: Failed to execute MFC. You can find the run's output in {out_filepath}, and the case dictionary in {os.path.join(test.get_dirpath(), "case.py")}.""")

            pack = packer.generate(test)
            pack.save(os.path.join(test.get_dirpath(), "pack.txt"))

            golden_filepath = os.path.join(test.get_dirpath(), "golden.txt")

            if self.mfc.args["generate"]:
                common.delete_file(golden_filepath)
                pack.save(golden_filepath)
            else:
                if not os.path.isfile(golden_filepath):
                    raise MFCException(f"Test {test}: Golden file doesn't exist! To generate golden files, use the '-g' flag.")

                packer.check_tolerance(test, pack, packer.load(golden_filepath), tol)

            cons.print(f"  [bold magenta]{test.get_uuid()}[/bold magenta]    {test.trace}")
        except Exception as exc:
            self.nFail = self.nFail + 1

            if not self.mfc.args["relentless"]:
                raise exc

            cons.print(f"[bold red]Failed test {test}.[/bold red]")
            cons.print(f"{exc}")
