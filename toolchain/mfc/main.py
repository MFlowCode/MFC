#!/usr/bin/env python3

import os
import signal
import getpass
import platform
import itertools

from mfc.util.common  import MFC_LOGO, MFCException, quit, delete_directory, format_list_to_string
from mfc.util.printer import cons

import mfc.args
import mfc.build
import mfc.cfg.user
import mfc.cfg.lock

import mfc.run.run
import mfc.tests.tests


class MFCState:
    def __init__(self) -> None:
        self.user = mfc.cfg.user.MFCUser()
        self.lock = mfc.cfg.lock.MFCLock(self.user)
        self.test = mfc.tests.tests.MFCTest(self)
        self.args = mfc.args.parse(self)
        self.run  = mfc.run.run.MFCRun(self)

        self.__handle_mode()
        self.__print_greeting()
        self.__run()


    def __handle_mode(self):
        # Handle mode change
        if self.args["mode"] != self.lock.mode:
            cons.print(f"[bold yellow]Switching to [bold magenta]{self.args['mode']}[/bold magenta] mode from [bold magenta]{self.lock.mode}[/bold magenta] mode:[/bold yellow]")
            self.lock.mode = self.args["mode"]
            self.lock.write()

            for target_name in mfc.build.get_mfc_target_names():
                t = mfc.build.get_target(target_name)
                dirpath = mfc.build.get_build_dirpath(t)
                cons.print(f"[bold red] - Removing {os.path.relpath(dirpath)}[/bold red]")
                delete_directory(dirpath)


    def __print_greeting(self):
        MFC_LOGO_LINES       = MFC_LOGO.splitlines()
        max_logo_line_length = max([ len(line) for line in MFC_LOGO_LINES ])

        host_line = f"{getpass.getuser()}@{platform.node()} [{platform.system()}]"

        MFC_SIDEBAR_LINES = [
            "",
            f"[bold]{host_line}[/bold]",
            '-' * len(host_line),
            "",
            "",
            f"[bold]--jobs:    [magenta]{self.args['jobs']}[/magenta][/bold]",
            f"[bold]--mode:    [magenta]{self.lock.mode}[/magenta][/bold]",
            f"[bold]--targets: {format_list_to_string([ f'[magenta]{target}[/magenta]' for target in self.args['targets']], 'None')}[/bold]",
            "",
            "",
            "[yellow]$ ./mfc.sh \[run, test, clean] --help[/yellow]",
        ]


        for a, b in itertools.zip_longest(MFC_LOGO_LINES, MFC_SIDEBAR_LINES):
            lhs = a.ljust(max_logo_line_length)
            rhs = b if b is not None else ''
            cons.print(
                f"[bold blue] {lhs} [/bold blue]  {rhs}",
                highlight=False
            )

        cons.print()


    def __run(self):
        if self.args["command"] == "test":
            self.test.execute()
        elif self.args["command"] == "run":
            self.run.run()
        elif self.args["command"] == "build":
            for target in self.args["targets"]:
                mfc.build.build_target(self, target)
        elif self.args["command"] == "clean":
            for target in self.args["targets"]:
                mfc.build.clean_target(self, target)


FILE_ISSUE_MSG = f"""\
We apologize for the inconvenience. If you believe this is an issue with MFC, \
please visit https://github.com/MFlowCode/MFC-develop/issues to file an issue.\
"""

if __name__ == "__main__":
    try:
        MFCState()
    except MFCException as exc:
        cons.reset()
        cons.print(f"""
--- [bold red]FATAL MFC ERROR[/bold red] ---

{str(exc)}
{FILE_ISSUE_MSG}
""")
        quit(signal.SIGTERM)
    except KeyboardInterrupt as exc:
        quit(signal.SIGTERM)
    except Exception as exc:
        cons.reset()
        cons.print_exception()
        cons.print(f"""
--- [bold red]FATAL MFC ERROR[/bold red] ---

An unexpected exception occurred:
{FILE_ISSUE_MSG}
""")
        quit(signal.SIGTERM)
