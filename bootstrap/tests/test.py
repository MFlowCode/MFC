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

from pathlib import Path

from tests.case import Case
from tests.cases import generate_filtered_cases

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

    def golden_file_compare_match(self, truth: str, candidate: str, tol):
        if truth.count('\n') != candidate.count('\n'):
            return (False, "Line count didn't match.")

        if "NaN" in truth:
            return (False, "NaN in golden file")

        if "NaN" in candidate:
            return (False, "NaN in packed file")

        for candidate_line in candidate.splitlines():
            if candidate_line == "":
                continue

            file_subpath: str = candidate_line.split(' ')[0]

            line_trusted: str = ""
            for l in truth.splitlines():
                if l.startswith(file_subpath):
                    line_trusted = l
                    break

            if len(line_trusted) == 0:
                continue

            numbers_cand  = [ float(x) for x in candidate_line.strip().split(' ')[1:] ]
            numbers_trust = [ float(x) for x in line_trusted.strip().split(' ')[1:]   ]

            # Different amount of spaces, means that there are more entires in one than in the other
            if len(numbers_cand) != len(numbers_trust):
                return (False, "Variable count didn't match.")

            # check values one by one
            for i in range(len(numbers_cand)):
                abs_delta = abs(numbers_cand[i]-numbers_trust[i])
                rel_diff  = abs(abs_delta/numbers_trust[i]) if numbers_trust[i] != 0 else 0
                if    (abs_delta > tol and rel_diff > tol):
                    percent_diff = rel_diff*100
                    return (False, f"Error margin is too high for the value #{i+1} in {file_subpath}: ~{round(percent_diff, 5)}% (~{round(abs_delta, 5)}).")

        # Both tests gave the same results within an acceptable tolerance
        return (True, "")

    def handle_case(self, test: Case):
        try:
            test.create_directory()

            if('qbmm' in test.params):
                tol = 1e-7
            elif('bubbles' in test.params):
                tol = 1e-10
            else:
                tol = 1e-12

            def on_test_errror(msg: str = "", term_out: str = ""):
                common.clear_line()

                uuid = test.get_uuid()
                rich.print(f"> Test #{uuid}: Failed!")
                if msg != "":
                    rich.print(f"> {msg}")

                filepath = f"{common.MFC_TESTDIR}/failed_{test.get_uuid()}.txt"
                common.file_write(filepath, f"""\
(1/3) Test #{uuid}:
- Test UUID: {uuid}
- Summary:   {test.trace}
- Location:  {test.get_dirpath()}
- Error:     {msg}
- When:      {common.get_datetime_str()}

(2/3) Test case:
{test.gen_json_dict_str()}

(3/3) Terminal output:
{term_out}
""")

                rich.print(f"> Please read {filepath} for more information.")
                raise common.MFCException("Testing failed (view above).")

            cmd = subprocess.run(
                f'./mfc.sh run "{test.get_dirpath()}/case.py" -m "{self.mfc.args["mode"]}" -c {test["ppn"]} -t pre_process simulation 2>&1',
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                universal_newlines=True, shell=True
            )
            
            common.file_write(f"{test.get_dirpath()}/out.txt", cmd.stdout)

            if cmd.returncode != 0:
                on_test_errror("MFC Execution Failed.", cmd.stdout)

            pack = self.pack_case_output(test)
            common.file_write(f"{test.get_dirpath()}/pack.txt", pack)

            golden_filepath = f"{test.get_dirpath()}/golden.txt"

            if self.mfc.args["generate"]:
                common.delete_file(golden_filepath)
                common.file_write(golden_filepath, pack)

            if not os.path.isfile(golden_filepath):
                common.clear_line()
                on_test_errror("Golden file doesn't exist! To generate golden files, use the '-g' flag.", cmd.stdout)

            golden_file_content = common.file_read(golden_filepath)
            bSuccess, errorMsg  = self.golden_file_compare_match(golden_file_content, pack, tol)
            
            if not bSuccess:
                on_test_errror(errorMsg, cmd.stdout)
            
        except BaseException as exc:
            print(exc)
            #traceback.print_exc()
            
            # Exit
            os.kill(os.getpid(), signal.SIGTERM)

    def pack_case_output(self, case: Case):
        result: str = ""

        case_dir = case.get_dirpath()
        D_dir    = f"{case_dir}/D/"

        for filepath in list(Path(D_dir).rglob("*.dat")):
            file_content   = common.file_read(filepath)
            short_filepath = str(filepath).replace(f'{case_dir}/', '')

            result += f"{short_filepath} " + re.sub(r' +', ' ', file_content.replace('\n', ' ')).strip() + '\n'

        return result
