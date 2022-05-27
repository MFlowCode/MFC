#!/usr/bin/env python3

import os


import common

from tests.case    import Case
from tests.cases   import generate_filtered_cases
from tests.threads import MFCTestThreadManager

from common import MFCException

import tests.pack


import rich


class MFCTest:
    def __init__(self, mfc):
        self.mfc    = mfc
        self.sched  = MFCTestThreadManager(self.mfc.args["jobs"])

    def test(self):
        rich.print("[bold][u]Test:[/u][/bold] (in tests/)")

        # Clear previous tests if we wish to (re)generate golden files
        if self.mfc.args["generate"]:
            common.delete_directory_recursive(common.MFC_TESTDIR)
            common.create_directory(common.MFC_TESTDIR)

        # Build mfc if required
        if not self.mfc.build.is_built("mfc"):
            rich.print(f"> [bold cyan]mfc[/bold cyan] needs (re)building...")
            self.mfc.build.build_target(f"mfc", "> > ")

        # Run cases with multiple threads (if available)
        rich.print(f" |-+------------+----------+----------+---------+")
        rich.print(f" | | tests/[bold magenta]UUID[/bold magenta] | Error RE |   Tol.   | Summary |")
        rich.print(f" |-+------------+----------+----------+---------+")
        self.sched.run(generate_filtered_cases(self.mfc.args), self.handle_case)

        rich.print(f"> Tested [bold green]✓[/bold green]")

    def handle_case(self, test: Case):
        test.create_directory()

        if test.params.get('qbmm', 'F') == 'T':
            tol = 1e-7
        elif test.params.get('bubbles', 'F') == 'T':
            tol = 1e-10
        else:
            tol = 1e-12

        cmd = test.run(self.mfc.args)

        common.file_write(f"{test.get_dirpath()}/out.txt", cmd.stdout)

        if cmd.returncode != 0:
            raise MFCException(f"tests/{test.get_uuid()}: Failed to execute MFC [{test.trace}]")

        pack = tests.pack.generate(test)
        pack.save(f"{test.get_dirpath()}/pack.txt")

        golden_filepath = f"{test.get_dirpath()}/golden.txt"

        if self.mfc.args["generate"]:
            common.delete_file(golden_filepath)
            pack.save(golden_filepath)

        if not os.path.isfile(golden_filepath):
            raise MFCException(f"tests/{test.get_uuid()}: Golden file doesn't exist! To generate golden files, use the '-g' flag. [{test.trace}]")

        error = tests.pack.check_tolerance(test, pack, tests.pack.load(golden_filepath), tol)
        rich.print(f" |->  [bold magenta]{test.get_uuid()}[/bold magenta]  | {error.relative:+0.1E} | {tol:+0.1E} | {test.trace})")
