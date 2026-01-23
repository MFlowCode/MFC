import time, typing, threading, dataclasses
import rich, rich.progress
import traceback

from .printer import cons

# Thresholds for long-running test notifications
# Interactive mode: dimension-aware thresholds
INTERACTIVE_THRESHOLDS = {
    1: 30.0,   # 1D: 30 seconds
    2: 60.0,   # 2D: 1 minute
    3: 120.0,  # 3D: 2 minutes
}

# Headless mode: fixed time-based thresholds (regardless of dimensionality)
HEADLESS_THRESHOLDS = (
  (2  * 60,  "[italic yellow]Still running[/italic yellow] (>2min)"),
  (10 * 60,  "[italic yellow]Still running[/italic yellow] (>10min)"),
  (30 * 60,  "[bold red]Still running[/bold red] (>30min, may be hanging)"),
)

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
class WorkerThreadHolder:  # pylint: disable=too-many-instance-attributes
    thread:  threading.Thread
    ppn:     int
    load:    float
    devices: typing.Optional[typing.Set[int]]
    task:    typing.Optional['Task'] = None
    start:   float = 0.0
    # Track which milestones we've already logged
    notified_interactive: bool = False  # First notification in interactive mode (time varies by dimensionality)
    notified_2m:  bool = False  # Headless mode: 2 minute milestone
    notified_10m: bool = False  # Headless mode: 10 minute milestone
    notified_30m: bool = False  # Headless mode: 30 minute milestone


@dataclasses.dataclass
class Task:
    ppn:  int
    func: typing.Callable
    args: typing.List[typing.Any]
    load: float

def sched(tasks: typing.List[Task], nThreads: int, devices: typing.Optional[typing.Set[int]] = None) -> None:  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    nAvailable: int = nThreads
    threads:    typing.List[WorkerThreadHolder] = []

    sched.LOAD = { id: 0.0 for id in devices or [] }

    def get_case_dimensionality(case: typing.Any) -> int:
        """
        Determine if a test case is 1D, 2D, or 3D based on grid parameters.
        
        Grid parameters (m, n, p) represent the number of cells in x, y, z directions.
        Returns 3 if p != 0, 2 if n != 0, otherwise 1. Defaults to 1D if params unavailable.
        """
        if not hasattr(case, 'params'):
            return 1  # Default to 1D if we can't determine

        params = case.params
        p = params.get('p', 0)
        n = params.get('n', 0)

        if p != 0:
            return 3  # 3D
        if n != 0:
            return 2  # 2D
        return 1  # 1D

    def get_threshold_for_case(case: typing.Any) -> float:
        """
        Get the dimension-aware time threshold (in seconds) for interactive mode notifications.
        
        Returns 30s for 1D, 60s for 2D, 120s for 3D tests.
        """
        dim = get_case_dimensionality(case)
        return INTERACTIVE_THRESHOLDS.get(dim, INTERACTIVE_THRESHOLDS[1])

    def notify_long_running_threads(  # pylint: disable=too-many-branches
        progress: rich.progress.Progress,
        running_tracker: typing.Optional[rich.progress.TaskID],
        interactive: bool
    ) -> None:
        """
        Monitor and notify about long-running tests.

        In interactive mode: prints once when a test crosses its dimension-aware threshold
        and updates the live progress bar. In headless mode: prints milestone notifications
        at 2, 10, and 30 minutes.
        """
        now = time.time()
        long_running_for_progress = []

        for holder in threads:
            if not holder.thread.is_alive():
                continue

            elapsed = now - holder.start
            case = holder.task.args[0] if holder.task and holder.task.args else None
            case_uuid  = case.get_uuid() if hasattr(case, "get_uuid") else "unknown"
            case_trace = getattr(case, "trace", "")

            # --- interactive: dimension-aware thresholds ---
            if interactive:
                threshold = get_threshold_for_case(case)

                if elapsed >= threshold:
                    long_running_for_progress.append((case_uuid, case_trace))

                    # Print explicit line once when crossing threshold
                    if not holder.notified_interactive:
                        dim = get_case_dimensionality(case)
                        dim_label = f"{dim}D"
                        time_label = f"{int(threshold)}s" if threshold < 60 else f"{threshold/60:.0f}min"
                        cons.print(
                            f"  [italic yellow]Still running[/italic yellow] ({dim_label}, >{time_label}) "
                            f"[bold magenta]{case_uuid}[/bold magenta]  {case_trace}"
                        )
                        holder.notified_interactive = True

            # --- headless: milestone notifications at 2, 10, 30 minutes ---
            else:
                # 2 minutes
                if (not holder.notified_2m) and elapsed >= 2 * 60:
                    cons.print(
                        f"  {HEADLESS_THRESHOLDS[0][1]} "
                        f"[bold magenta]{case_uuid}[/bold magenta]  {case_trace}"
                    )
                    holder.notified_2m = True

                # 10 minutes
                if (not holder.notified_10m) and elapsed >= 10 * 60:
                    cons.print(
                        f"  {HEADLESS_THRESHOLDS[1][1]} "
                        f"[bold magenta]{case_uuid}[/bold magenta]  {case_trace}"
                    )
                    holder.notified_10m = True

                # 30 minutes
                if (not holder.notified_30m) and elapsed >= 30 * 60:
                    cons.print(
                        f"  {HEADLESS_THRESHOLDS[2][1]} "
                        f"[bold magenta]{case_uuid}[/bold magenta]  {case_trace}"
                    )
                    holder.notified_30m = True

        # update the interactive "Running" row
        if interactive and running_tracker is not None:
            if long_running_for_progress:
                summary = ", ".join(uuid for uuid, _ in long_running_for_progress[:5])
                if len(long_running_for_progress) > 5:
                    summary += f", +{len(long_running_for_progress) - 5} more"
                progress.update(running_tracker, description=f"Running (long): {summary}")
            else:
                progress.update(running_tracker, description="Running (long): none")

    def join_first_dead_thread(progress, complete_tracker, interactive: bool) -> None:
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
                if threadHolder.thread.exc is not None:
                    # Unhandled exception - propagate with full traceback if available
                    if threadHolder.thread.exc_info:
                        error_msg = f"Worker thread {threadID} failed with unhandled exception:\n{threadHolder.thread.exc_info}"
                        raise RuntimeError(error_msg) from threadHolder.thread.exc
                    raise threadHolder.thread.exc

                # Print completion message for long-running tests in interactive mode
                if interactive and threadHolder.notified_interactive:
                    elapsed = time.time() - threadHolder.start
                    case = threadHolder.task.args[0] if threadHolder.task and threadHolder.task.args else None
                    case_uuid  = case.get_uuid() if hasattr(case, "get_uuid") else "unknown"
                    case_trace = getattr(case, "trace", "")
                    cons.print(
                        f"  [italic green]Completed[/italic green] (after {elapsed:.1f}s) "
                        f"[bold magenta]{case_uuid}[/bold magenta]  {case_trace}"
                    )

                nAvailable += threadHolder.ppn
                for device in threadHolder.devices or set():
                    sched.LOAD[device] -= threadHolder.load / threadHolder.ppn

                progress.advance(complete_tracker)

                del threads[threadID]

                break

    with rich.progress.Progress(console=cons.raw, transient=True) as progress:
        interactive      = cons.raw.is_terminal
        queue_tracker    = progress.add_task("Queued    ", total=len(tasks))
        complete_tracker = progress.add_task("Completed ", total=len(tasks))
        running_tracker  = progress.add_task("Running   ", total=None) if interactive else None

        # Queue Tests
        for task in tasks:
            # Wait until there are threads available
            while nAvailable < task.ppn:
                # This is important if "-j 1" is used (the default) since there
                # are test cases that require test.ppn=2
                if task.ppn > nThreads and nAvailable > 0:
                    break

                # Keep track of threads that are done
                join_first_dead_thread(progress, complete_tracker, interactive)

                # Notify about long-running threads
                notify_long_running_threads(progress, running_tracker, interactive)

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

            threads.append(
                WorkerThreadHolder(
                    thread=thread,
                    ppn=task.ppn,
                    load=task.load,
                    devices=use_devices,
                    task=task,
                    start=time.time(),
                )
            )

        # Wait for the last tests to complete (MOVED INSIDE CONTEXT)
        while len(threads) != 0:
            # Keep track of threads that are done
            join_first_dead_thread(progress, complete_tracker, interactive)

            # Notify about long-running threads
            notify_long_running_threads(progress, running_tracker, interactive)

            # Do not overwhelm this core with this loop
            time.sleep(0.05)

sched.LOAD = {}
