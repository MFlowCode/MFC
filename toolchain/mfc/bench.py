import os, json, time, typing, datetime, subprocess

import rich.table

from .printer import cons
from .state   import ARG
from .build   import PRE_PROCESS, SIMULATION, build_targets
from .common  import system, MFC_SUBDIR
from .        import sched

def bench():
    build_targets([PRE_PROCESS, SIMULATION])
    
    cons.print("[bold]Benchmarking [magenta]simulation[/magenta]:[/bold]")
    cons.indent()
    
    CASES   = ["1D_bubblescreen", "1D_exercise_WENO", "1D_kapilashocktube"]
    RESULTS = []
    
    table = rich.table.Table(show_lines=False, show_edge=False)
    table.add_column("Case")
    table.add_column("(Simulation) Runtime (s)")
    
    def __worker(case: str, devices: typing.Set[int]):
        nonlocal RESULTS
        
        system(["./mfc.sh", "run", f"examples/{case}/case.py", "--no-build", "-t", "pre_process"], stdout=subprocess.DEVNULL)
        start   = time.monotonic()
        system(["./mfc.sh", "run", f"examples/{case}/case.py", "--no-build", "-t", "simulation"], stdout=subprocess.DEVNULL)
        end     = time.monotonic()
        runtime = datetime.timedelta(seconds=end - start).total_seconds()

        RESULTS.append({
            "name":  f"Simulation: {case}",
            "unit":  "seconds",
            "value": runtime
        })
        
        table.add_row(case, str(runtime))
    
    tasks: typing.List[sched.Task] = [
        sched.Task(1, __worker, [ case ], 1) for case in CASES
    ]
    
    cons.print()
    nThreads = min(ARG('jobs'), len(ARG('gpus'))) if ARG("gpu") else ARG('jobs')
    if ARG('case_optimization'):
        nThreads = 1

    sched.sched(tasks, nThreads, ARG("gpus"))
    cons.print()
    cons.unindent()
    cons.print("[bold]Benchmark Results:[/bold]")
    cons.print()
    cons.raw.print(table)
    cons.print()
    
    filepath = os.path.join(MFC_SUBDIR, "bench.json")
    with open(filepath, "w") as f:
        json.dump(RESULTS, f)
    
    cons.print(f"[bold green]✓[/bold green] Saved results to [magenta]{filepath}[/magenta].")
