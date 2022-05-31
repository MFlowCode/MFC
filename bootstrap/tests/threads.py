import time
import threading
import dataclasses


from tests.case import Case


import rich
import rich.progress


class TestThread(threading.Thread):
    def __init__(self, *args, **kwargs):
        self.exc = None

        threading.Thread.__init__(self, *args, **kwargs)

    def run(self):
        try:
            if self._target:
                self._target(*self._args, **self._kwargs)
        except Exception as exc:
            self.exc = exc


@dataclasses.dataclass
class TestThreadHolder:
    thread: threading.Thread
    ppn:    int


class MFCTestThreadManager:
    def __init__(self, nThreads: int) -> None:
        self.threads    = []
        self.nThreads   = nThreads
        self.nAvailable = self.nThreads

    def join_first_dead_thread(self, progress, complete_tracker) -> None:
        for threadID, threadHolder in enumerate(self.threads):
            threadHolder: TestThreadHolder

            if not threadHolder.thread.is_alive():
                if threadHolder.thread.exc != None:
                    raise threadHolder.thread.exc

                self.nAvailable += threadHolder.ppn
                progress.advance(complete_tracker)

                del self.threads[threadID]

                break

    def run(self, cases: list, handle_case) -> None:
        with rich.progress.Progress() as progress:
            queue_tracker    = progress.add_task("Queued   ", total=len(cases))
            complete_tracker = progress.add_task("Completed", total=len(cases))

            # Queue Tests
            for i, test in enumerate(cases):
                test: Case

                # Wait until there are threads available
                while self.nAvailable < test.ppn:
                    # This is important if "-j 1" is used (the default) since there
                    # are test cases that require test.ppn=2
                    if test.ppn > self.nThreads and self.nAvailable > 0:
                        break

                    # Keep track of threads that are done
                    self.join_first_dead_thread(progress, complete_tracker)

                    # Do not overwhelm this core with this loop
                    time.sleep(0.05)

                # Launch Thread
                progress.advance(queue_tracker)

                thread = TestThread(target=handle_case, args=(test,))
                thread.start()

                self.threads.append(TestThreadHolder(thread, test.ppn))
                self.nAvailable -= test.ppn

            # Wait for the lasts tests to complete
            while len(self.threads) != 0:
                # Keep track of threads that are done
                self.join_first_dead_thread(progress, complete_tracker)

                # Do not overwhelm this core with this loop
                time.sleep(0.05)
