#!/usr/bin/env python3

import rich
import rich.console

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
        rich.print(common.MFC_HEADER)

        self.user = user.MFCUser()
        self.lock = lock.MFCLock(self.user)
        self.test = tests.tests.MFCTest(self)
        self.args = args.parse(self)
        self.run  = run.run.MFCRun(self)

        # Handle mode change
        if self.args["mode"] != self.lock.mode:
            rich.print(f"[bold yellow]Switching to [bold magenta]{self.args['mode']}[/bold magenta] mode from [bold magenta]{self.lock.mode}[/bold magenta] mode.[/bold yellow]")
            self.lock.mode = self.args["mode"]
            self.lock.write()

            rich.print(f"[bold red]Purging build files...[/bold red]")

            for dep_name in build.get_target("mfc").requires:
                t = build.get_target(dep_name)
                dirpath = build.get_build_dirpath(t)
                rich.print(f"[bold red] - {dep_name}: {os.path.relpath(dirpath)}[/bold red]")
                common.delete_directory_recursive(dirpath)
        else:
            rich.print(f"[bold yellow]You are currently in [bold magenta]{self.lock.mode}[/bold magenta] mode.[/bold yellow]")
        rich.print("")

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
        rich.print(f"""
 --- [bold red]FATAL MFC ERROR[/bold red] ---
{str(exc)}
{FILE_ISSUE_MSG}
""")
        common.quit(signal.SIGTERM)
    except KeyboardInterrupt as exc:
        common.quit(signal.SIGTERM)
    except Exception as exc:
        rich.console.Console().print_exception()
        rich.print(f"""
 --- [bold red]FATAL MFC ERROR[/bold red] ---
An unexpected exception occurred:
{FILE_ISSUE_MSG}
""")
        common.quit(signal.SIGTERM)
