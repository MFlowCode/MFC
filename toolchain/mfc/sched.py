import time, typing, threading, dataclasses
import rich, rich.progress
import traceback

from .printer import cons

class WorkerThread(threading.Thread):
    def __init__(self, *args, **kwargs):
        self.exc = None
        self.exc_info = None  # Store full exception information for better debugging
        self.completed_successfully = False  # Track if the target function completed

        threading.Thread.__init__(self, *args, **kwargs)

    def run(self):
        try:
            if self._target:
                self._target(*self._args, **self._kwargs)
                self.completed_successfully = True  # Mark as completed successfully
        except Exception as exc:
            self.exc = exc
            # Store the full traceback for better error reporting
            self.exc_info = traceback.format_exc()


@dataclasses.dataclass
class WorkerThreadHolder:
    thread:  threading.Thread
    ppn:     int
    load:    float
    devices: typing.Set[int]


@dataclasses.dataclass
class Task:
    ppn:  int
    func: typing.Callable
    args: typing.List[typing.Any]
    load: float

def sched(tasks: typing.List[Task], nThreads: int, devices: typing.Set[int] = None) -> None:
    nAvailable: int = nThreads
    threads:    typing.List[WorkerThreadHolder] = []

    sched.LOAD = { id: 0.0 for id in devices or [] }

    def join_first_dead_thread(progress, complete_tracker) -> None:
        nonlocal threads, nAvailable

        for threadID, threadHolder in enumerate(threads):
            # Check if thread is not alive OR if it's been running for too long
            thread_not_alive = not threadHolder.thread.is_alive()

            if thread_not_alive:
                # Properly join the thread with timeout to prevent infinite hangs
                try:
                    threadHolder.thread.join(timeout=30.0)  # 30 second timeout

                    # Double-check that thread actually finished joining
                    if threadHolder.thread.is_alive():
                        # Thread didn't finish within timeout - this is a serious issue
                        raise RuntimeError(f"Thread {threadID} failed to join within 30 seconds timeout. "
                                         f"Thread may be hung or in an inconsistent state.")

                except Exception as join_exc:
                    # Handle join-specific exceptions with more context
                    raise RuntimeError(f"Failed to join thread {threadID}: {join_exc}. "
                                     f"This may indicate a system threading issue or hung test case.") from join_exc

                # Check for and propagate any exceptions that occurred in the worker thread
                # But only if the worker function didn't complete successfully
                # (This allows test failures to be handled gracefully by handle_case)
                if threadHolder.thread.exc is not None:
                    if threadHolder.thread.completed_successfully:
                        # Test framework handled the exception gracefully (e.g., test failure)
                        # Don't re-raise - this is expected behavior
                        pass
                    # Unhandled exception - this indicates a real problem
                    elif hasattr(threadHolder.thread, 'exc_info') and threadHolder.thread.exc_info:
                        error_msg = f"Worker thread {threadID} failed with unhandled exception:\n{threadHolder.thread.exc_info}"
                        raise RuntimeError(error_msg) from threadHolder.thread.exc
                    else:
                        raise threadHolder.thread.exc

                nAvailable += threadHolder.ppn
                for device in threadHolder.devices or set():
                    sched.LOAD[device] -= threadHolder.load / threadHolder.ppn

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

            use_devices = None
            # Use the least loaded devices
            if devices is not None:
                use_devices = set()
                for _ in range(task.ppn):
                    device = min(sched.LOAD.items(), key=lambda x: x[1])[0]
                    sched.LOAD[device] += task.load / task.ppn
                    use_devices.add(device)

            nAvailable -= task.ppn

            thread = WorkerThread(target=task.func, args=tuple(task.args) + (use_devices,))
            thread.start()

            threads.append(WorkerThreadHolder(thread, task.ppn, task.load, use_devices))

        # Wait for the last tests to complete (MOVED INSIDE CONTEXT)
        while len(threads) != 0:
            # Keep track of threads that are done
            join_first_dead_thread(progress, complete_tracker)

            # Do not overwhelm this core with this loop
            time.sleep(0.05)

sched.LOAD = {}
