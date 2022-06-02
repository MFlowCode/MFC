#!/usr/bin/env python3

import rich
import rich.console

import args
import user
import build
import common

import signal

import run.run
import tests.tests


class MFCState:
    def __init__(self) -> None:
        rich.print(common.MFC_HEADER)

        self.user = user.MFCUser()
        self.args = args.parse(self)
        self.test = tests.tests.MFCTest(self)
        self.run  = run.run.MFCRun(self)

        if self.args["command"] == "test":
            self.test.test()
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
