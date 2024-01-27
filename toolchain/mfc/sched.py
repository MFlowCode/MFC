import time, typing, threading, dataclasses
import rich, rich.progress

from .printer import cons




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
            if not threadHolder.thread.is_alive():
                if threadHolder.thread.exc is not None:
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


        # Wait for the lasts tests to complete
        while len(threads) != 0:
            # Keep track of threads that are done
            join_first_dead_thread(progress, complete_tracker)

            # Do not overwhelm this core with this loop
            time.sleep(0.05)

sched.LOAD = {}
