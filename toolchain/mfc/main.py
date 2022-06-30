#!/usr/bin/env python3

from mfc.printer import cons

import args
import user
import lock
import build
import common

import os
import signal

import run.run
import tests.tests


class MFCState:
    def __init__(self) -> None:
        cons.print(common.MFC_HEADER)

        self.user = user.MFCUser()
        self.lock = lock.MFCLock(self.user)
        self.test = tests.tests.MFCTest(self)
        self.args = args.parse(self)
        self.run  = run.run.MFCRun(self)

        # Handle mode change
        if self.args["mode"] != self.lock.mode:
            cons.print(f"[bold yellow]Switching to [bold magenta]{self.args['mode']}[/bold magenta] mode from [bold magenta]{self.lock.mode}[/bold magenta] mode:[/bold yellow]")
            self.lock.mode = self.args["mode"]
            self.lock.write()

            for dep_name in build.get_target("mfc").requires:
                t = build.get_target(dep_name)
                dirpath = build.get_build_dirpath(t)
                cons.print(f"[bold red] - Removing {os.path.relpath(dirpath)}[/bold red]")
                common.delete_directory(dirpath)

        cons.print(f"[bold yellow]You are currently in [bold magenta]{self.lock.mode}[/bold magenta] mode.[/bold yellow]\n")

        if self.args["command"] == "test":
            self.test.execute()
        elif self.args["command"] == "run":
            self.run.run()
        elif self.args["command"] == "build":
            for target in self.args["targets"]:
                build.build_target(self, target)
        elif self.args["command"] == "clean":
            for target in self.args["targets"]:
                build.clean_target(self, target)


FILE_ISSUE_MSG = f"""\
We apologize for the inconvenience. If you believe this is an issue with MFC, \
please visit https://github.com/MFlowCode/MFC-develop/issues to file an issue.\
"""

if __name__ == "__main__":
    try:
        MFCState()
    except common.MFCException as exc:
        cons.reset()
        cons.print(f"""
--- [bold red]FATAL MFC ERROR[/bold red] ---

{str(exc)}
{FILE_ISSUE_MSG}
""")
        common.quit(signal.SIGTERM)
    except KeyboardInterrupt as exc:
        common.quit(signal.SIGTERM)
    except Exception as exc:
        cons.reset()
        cons.print_exception()
        cons.print(f"""
--- [bold red]FATAL MFC ERROR[/bold red] ---

An unexpected exception occurred:
{FILE_ISSUE_MSG}
""")
        common.quit(signal.SIGTERM)
