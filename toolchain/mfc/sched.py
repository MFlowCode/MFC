import time, typing, threading, dataclasses
import rich, rich.progress
import traceback

from .printer import cons
from .state import ARG

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

    # Debug logging setup
    gpu_mode = devices is not None and len(devices) > 0
    debug_enabled = ARG("sched_debug", False)  # Check for --sched-debug flag
    
    def debug_log(msg):
        if debug_enabled:
            cons.print(msg)
    
    if debug_enabled:
        debug_log(f"[SCHED DEBUG] Starting scheduler: {len(tasks)} tasks, {nThreads} threads, GPU mode: {gpu_mode}")
        if gpu_mode:
            debug_log(f"[SCHED DEBUG] GPU devices: {devices}")

    def join_first_dead_thread(progress, complete_tracker) -> None:
        nonlocal threads, nAvailable

        debug_log(f"[SCHED DEBUG] Checking {len(threads)} active threads for completion")
        
        for threadID, threadHolder in enumerate(threads):
            # Check if thread is not alive OR if it's been running for too long
            thread_not_alive = not threadHolder.thread.is_alive()
            
            debug_log(f"[SCHED DEBUG] Thread {threadID}: alive={threadHolder.thread.is_alive()}, devices={threadHolder.devices}")

            if thread_not_alive:
                debug_log(f"[SCHED DEBUG] Thread {threadID} detected as dead, attempting to join...")
                
                # Properly join the thread with timeout to prevent infinite hangs
                join_start_time = time.time()
                timeout_duration = 120.0 if gpu_mode else 30.0  # Longer timeout for GPU
                
                debug_log(f"[SCHED DEBUG] Joining thread {threadID} with {timeout_duration}s timeout...")
                
                try:
                    threadHolder.thread.join(timeout=timeout_duration)
                    join_end_time = time.time()
                    join_duration = join_end_time - join_start_time
                    
                    debug_log(f"[SCHED DEBUG] Thread {threadID} join completed in {join_duration:.2f}s")

                    # Double-check that thread actually finished joining
                    if threadHolder.thread.is_alive():
                        # Thread didn't finish within timeout - this is a serious issue
                        debug_log(f"[SCHED DEBUG] ERROR: Thread {threadID} still alive after {timeout_duration}s timeout!")
                        debug_log(f"[SCHED DEBUG] Thread {threadID} devices: {threadHolder.devices}")
                        debug_log(f"[SCHED DEBUG] Thread {threadID} exception: {threadHolder.thread.exc}")
                        raise RuntimeError(f"Thread {threadID} failed to join within {timeout_duration} seconds timeout. "
                                         f"Thread may be hung or in an inconsistent state. "
                                         f"GPU devices: {threadHolder.devices}")

                except Exception as join_exc:
                    # Handle join-specific exceptions with more context
                    debug_log(f"[SCHED DEBUG] Exception during thread {threadID} join: {join_exc}")
                    raise RuntimeError(f"Failed to join thread {threadID}: {join_exc}. "
                                     f"This may indicate a system threading issue or hung test case. "
                                     f"GPU devices: {threadHolder.devices}") from join_exc

                debug_log(f"[SCHED DEBUG] Thread {threadID} successfully joined")

                # Check for and propagate any exceptions that occurred in the worker thread
                # But only if the worker function didn't complete successfully
                # (This allows test failures to be handled gracefully by handle_case)
                if threadHolder.thread.exc is not None:
                    debug_log(f"[SCHED DEBUG] Thread {threadID} had exception: {threadHolder.thread.exc}")
                    debug_log(f"[SCHED DEBUG] Thread {threadID} completed successfully: {threadHolder.thread.completed_successfully}")
                    
                    if threadHolder.thread.completed_successfully:
                        # Test framework handled the exception gracefully (e.g., test failure)
                        # Don't re-raise - this is expected behavior
                        debug_log(f"[SCHED DEBUG] Thread {threadID} exception was handled gracefully by test framework")
                        pass
                    # Unhandled exception - this indicates a real problem
                    elif hasattr(threadHolder.thread, 'exc_info') and threadHolder.thread.exc_info:
                        error_msg = f"Worker thread {threadID} failed with unhandled exception:\n{threadHolder.thread.exc_info}"
                        debug_log(f"[SCHED DEBUG] Thread {threadID} had unhandled exception!")
                        raise RuntimeError(error_msg) from threadHolder.thread.exc
                    else:
                        debug_log(f"[SCHED DEBUG] Thread {threadID} had unhandled exception without details")
                        raise threadHolder.thread.exc

                # Update scheduler state
                nAvailable += threadHolder.ppn
                for device in threadHolder.devices or set():
                    old_load = sched.LOAD[device]
                    sched.LOAD[device] -= threadHolder.load / threadHolder.ppn
                    debug_log(f"[SCHED DEBUG] Device {device} load: {old_load:.3f} -> {sched.LOAD[device]:.3f}")

                progress.advance(complete_tracker)

                debug_log(f"[SCHED DEBUG] Thread {threadID} cleanup complete, removing from active threads")
                del threads[threadID]

                break
        
        debug_log(f"[SCHED DEBUG] join_first_dead_thread completed, {len(threads)} threads remaining")

    with rich.progress.Progress(console=cons.raw, transient=True) as progress:
        queue_tracker    = progress.add_task("Queued   ", total=len(tasks))
        complete_tracker = progress.add_task("Completed", total=len(tasks))

        debug_log(f"[SCHED DEBUG] Starting task queue processing...")

        # Queue Tests
        for task_idx, task in enumerate(tasks):
            debug_log(f"[SCHED DEBUG] Processing task {task_idx+1}/{len(tasks)}: ppn={task.ppn}, load={task.load}")
            
            # Wait until there are threads available
            while nAvailable < task.ppn:
                debug_log(f"[SCHED DEBUG] Waiting for resources: need {task.ppn}, have {nAvailable}")
                
                # This is important if "-j 1" is used (the default) since there
                # are test cases that require test.ppn=2
                if task.ppn > nThreads and nAvailable > 0:
                    debug_log(f"[SCHED DEBUG] Task requires more threads ({task.ppn}) than available ({nThreads}), but some are free ({nAvailable})")
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
                debug_log(f"[SCHED DEBUG] Assigning GPU devices for task {task_idx+1}")
                for device_idx in range(task.ppn):
                    device = min(sched.LOAD.items(), key=lambda x: x[1])[0]
                    sched.LOAD[device] += task.load / task.ppn
                    use_devices.add(device)
                    debug_log(f"[SCHED DEBUG] Assigned device {device} (load now: {sched.LOAD[device]:.3f})")

            nAvailable -= task.ppn

            debug_log(f"[SCHED DEBUG] Starting thread for task {task_idx+1}, devices: {use_devices}")
            thread = WorkerThread(target=task.func, args=tuple(task.args) + (use_devices,))
            thread.start()

            threads.append(WorkerThreadHolder(thread, task.ppn, task.load, use_devices))
            debug_log(f"[SCHED DEBUG] Thread started for task {task_idx+1}, {len(threads)} active threads")

        debug_log(f"[SCHED DEBUG] All tasks queued, waiting for completion...")

        # Wait for the last tests to complete (MOVED INSIDE CONTEXT)
        while len(threads) != 0:
            debug_log(f"[SCHED DEBUG] Waiting for {len(threads)} threads to complete...")
            
            # Keep track of threads that are done
            join_first_dead_thread(progress, complete_tracker)

            # Do not overwhelm this core with this loop
            time.sleep(0.05)

        debug_log(f"[SCHED DEBUG] Scheduler completed successfully!")

sched.LOAD = {}
