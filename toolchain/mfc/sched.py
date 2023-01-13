import time, typing, threading, dataclasses

from .printer import cons

import rich, rich.progress


class WorkerThread(threading.Thread):
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
class WorkerThreadHolder:
    thread: threading.Thread
    ppn:    int


@dataclasses.dataclass
class Task:
    ppn:  int
    func: typing.Callable
    args: typing.List[typing.Any]


def sched(tasks: typing.List[Task], nThreads: int):
    nAvailable: int = nThreads
    threads:    typing.List[WorkerThreadHolder] = []


    def join_first_dead_thread(progress, complete_tracker) -> None:
        nonlocal threads, nAvailable
        
        for threadID, threadHolder in enumerate(threads):
            if not threadHolder.thread.is_alive():
                if threadHolder.thread.exc != None:
                    raise threadHolder.thread.exc

                nAvailable += threadHolder.ppn
                progress.advance(complete_tracker)

                del threads[threadID]

                break


    with rich.progress.Progress(console=cons.raw, transient=True) as progress:
        queue_tracker    = progress.add_task("Queued   ", total=len(tasks))
        complete_tracker = progress.add_task("Completed", total=len(tasks))

        # Queue Tests
        for task in tasks:
            # Wait until there are threads available
            while nAvailable < task.ppn:
                # This is important if "-j 1" is used (the default) since there
                # are test cases that require test.ppn=2
                if task.ppn > nThreads and nAvailable > 0:
                    break

                # Keep track of threads that are done
                join_first_dead_thread(progress, complete_tracker)

                # Do not overwhelm this core with this loop
                time.sleep(0.05)

            # Launch Thread
            progress.advance(queue_tracker)

            thread = WorkerThread(target=task.func, args=tuple(task.args))
            thread.start()

            threads.append(WorkerThreadHolder(thread, task.ppn))
            nAvailable -= task.ppn


        # Wait for the lasts tests to complete
        while len(threads) != 0:
            # Keep track of threads that are done
            join_first_dead_thread(progress, complete_tracker)

            # Do not overwhelm this core with this loop
            time.sleep(0.05)
