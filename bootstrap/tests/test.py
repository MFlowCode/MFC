#!/usr/bin/env python3

import common

import os
import re
import time
import signal
import threading
#import traceback
import subprocess
import dataclasses

from tests.case  import Case
from tests.cases import generate_filtered_cases

from common import MFCException

import tests.pack


import rich
import rich.progress


@dataclasses.dataclass
class TestThreadHolder:
    thread: threading.Thread
    ppn:    int


class MFCTest:
    def __init__(self, mfc):
        self.mfc = mfc

    def test(self):
        rich.print("[bold][u]Test:[/u][/bold]")

        # Generate (re)case-generation
        if self.mfc.args["generate"]:
            common.delete_directory_recursive(common.MFC_TESTDIR)
            common.create_directory(common.MFC_TESTDIR)

        # Build mfc if required
        if not self.mfc.build.is_built("mfc"):
            rich.print(f"> [bold cyan]mfc[/bold cyan] needs (re)building...")
            self.mfc.build.build_target(f"mfc", "> > ")

        cases: list = generate_filtered_cases(self.mfc.args)

        with rich.progress.Progress() as progress:
            queue_tracker    = progress.add_task("Queued   ", total=len(cases))
            complete_tracker = progress.add_task("Completed", total=len(cases))

            nAvailableThreads = self.mfc.args["jobs"]
            threads           = []

            # Queue Tests
            for i, test in enumerate(cases):
                test: Case

                ppn = test["ppn"]

                # Wait until there are threads available
                while nAvailableThreads < ppn:
                    # This is important if "-j 1" is used (the default) since there
                    # are test cases that require ppn=2
                    if ppn > self.mfc.args["jobs"] and nAvailableThreads > 0:
                        break

                    # Keep track of threads that are done
                    for threadID, threadHolder in enumerate(threads):
                        threadHolder: TestThreadHolder

                        if not threadHolder.thread.is_alive():
                            nAvailableThreads += threadHolder.ppn
                            progress.advance(complete_tracker)
                            del threads[threadID]
                            break
                    
                    # Do not overwhelm this core with this loop
                    time.sleep(0.2)

                rich.print(f" |-> {str(i+1).zfill(len(str(len(cases))))}/{len(cases)} ([bold magenta]{test.get_uuid()}[/bold magenta]): {test.trace}")
                progress.advance(queue_tracker)
                thread = threading.Thread(target=self.handle_case, args=(test,))
                thread.start()
                threads.append(TestThreadHolder(thread, ppn))
                nAvailableThreads -= ppn

            # Wait for the lasts tests to complete
            while len(threads) != 0:
                for threadID, threadHolder in enumerate(threads):
                    threadHolder: TestThreadHolder

                    if not threadHolder.thread.is_alive():
                        del threads[threadID]
                        progress.advance(complete_tracker)
                        break
            
            time.sleep(0.2)

        rich.print(f"> Tested [bold green]✓[/bold green]")

    def handle_case(self, test: Case):
        try:
            test.create_directory()

            if('qbmm' in test.params):
                tol = 1e-7
            elif('bubbles' in test.params):
                tol = 1e-10
            else:
                tol = 1e-12

            cmd = test.run(self.mfc.args)

            common.file_write(f"{test.get_dirpath()}/out.txt", cmd.stdout)

            if cmd.returncode != 0:
                raise MFCException(f"tests/{test.get_uuid()}: Failed to execute MFC.")

            pack = tests.pack.generate(test)
            pack.save(f"{test.get_dirpath()}/pack.txt")

            golden_filepath = f"{test.get_dirpath()}/golden.txt"

            if self.mfc.args["generate"]:
                common.delete_file(golden_filepath)
                pack.save(golden_filepath)

            if not os.path.isfile(golden_filepath):
                raise MFCException(f"tests/{test.get_uuid()}: Golden file doesn't exist! To generate golden files, use the '-g' flag.")

            tests.pack.check_tolerance(test.get_uuid(), pack, tests.pack.load(golden_filepath), tol)
            
        except BaseException as exc:
            print(exc)
            #traceback.print_exc()
            
            # Exit
            os.kill(os.getpid(), signal.SIGTERM)
