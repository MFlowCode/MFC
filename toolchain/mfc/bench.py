import os, sys, uuid, subprocess, dataclasses, typing, math, json

import rich.table

from .printer import cons
from .state   import ARG, CFG
from .build   import get_targets, DEFAULT_TARGETS, SIMULATION
from .common  import system, MFC_BENCH_FILEPATH, MFC_BUILD_DIR, format_list_to_string
from .common  import file_load_yaml, file_dump_yaml, create_directory
from .common  import MFCException


@dataclasses.dataclass
class BenchCase:
    slug: str
    path: str
    args: typing.List[str]

def bench(targets = None):
    if targets is None:
        targets = ARG("targets")

    targets = get_targets(targets)

    bench_dirpath = os.path.join(MFC_BUILD_DIR, "benchmarks", str(uuid.uuid4())[:4])
    create_directory(bench_dirpath)

    cons.print()
    cons.print(f"[bold]Benchmarking {format_list_to_string(ARG('targets'), 'magenta')} ([magenta]{os.path.relpath(bench_dirpath)}[/magenta]):[/bold]")
    cons.indent()
    cons.print()

    CASES = [ BenchCase(**case) for case in file_load_yaml(MFC_BENCH_FILEPATH) ]

    for case in CASES:
        case.args = case.args + ARG("--")
        case.path = os.path.abspath(case.path)

    results = {
        "metadata": {
            "invocation": sys.argv[1:],
            "lock":       dataclasses.asdict(CFG())
        },
        "cases": {},
    }

    for i, case in enumerate(CASES):
        summary_filepath = os.path.join(bench_dirpath, f"{case.slug}.yaml")
        log_filepath     = os.path.join(bench_dirpath, f"{case.slug}.out")

        cons.print(f"{str(i+1).zfill(len(CASES) // 10 + 1)}/{len(CASES)}: {case.slug} @ [bold]{os.path.relpath(case.path)}[/bold]")
        cons.indent()
        cons.print()
        cons.print(f"> Log:     [bold]{os.path.relpath(log_filepath)}[/bold]")
        cons.print(f"> Summary: [bold]{os.path.relpath(summary_filepath)}[/bold]")

        with open(log_filepath, "w") as log_file:
            system(
                ["./mfc.sh", "run", case.path, "--case-optimization"] +
                ["--targets"] + [t.name for t in targets] +
                ["--output-summary", summary_filepath] +
                case.args +
                ["--", "--gbpp", ARG('mem')],
                stdout=log_file,
                stderr=subprocess.STDOUT)

        results["cases"][case.slug] = {
            "description":    dataclasses.asdict(case),
            "output_summary": file_load_yaml(summary_filepath),
        }

    file_dump_yaml(ARG("output"), results)

    cons.print(f"Wrote results to [bold magenta]{os.path.relpath(ARG('output'))}[/bold magenta].")

    cons.unindent()


# TODO: This function is too long and not nicely written at all. Someone should
#       refactor it...
# pylint: disable=too-many-branches
def diff():
    lhs, rhs = file_load_yaml(ARG("lhs")), file_load_yaml(ARG("rhs"))
    cons.print(f"[bold]Comparing Benchmarks: Speedups from [magenta]{os.path.relpath(ARG('lhs'))}[/magenta] to [magenta]{os.path.relpath(ARG('rhs'))}[/magenta] are displayed below. Thus, numbers > 1 represent increases in performance.[/bold]")

    if lhs["metadata"] != rhs["metadata"]:
        _lock_to_str = lambda lock: ' '.join([f"{k}={v}" for k, v in lock.items()])

        cons.print(f"""\
[bold yellow]Warning[/bold yellow]: Metadata in lhs and rhs are not equal.
    This could mean that the benchmarks are not comparable (e.g. one was run on CPUs and the other on GPUs).
    lhs:
    * Invocation: [magenta]{' '.join(lhs['metadata']['invocation'])}[/magenta]
    * Modes:      {_lock_to_str(lhs['metadata']['lock'])}
    rhs:
    * Invocation: {' '.join(rhs['metadata']['invocation'])}
    * Modes:      [magenta]{_lock_to_str(rhs['metadata']['lock'])}[/magenta]
        """)

    slugs = set(lhs["cases"].keys()) & set(rhs["cases"].keys())
    if len(slugs) not in [len(lhs["cases"]), len(rhs["cases"])]:
        cons.print(f"""\
[bold yellow]Warning[/bold yellow]: Cases in lhs and rhs are not equal.
    * rhs cases: {', '.join(set(rhs['cases'].keys()) - slugs)}.
    * lhs cases: {', '.join(set(lhs['cases'].keys()) - slugs)}.
    Using intersection: {slugs} with {len(slugs)} elements.
        """)

    cb_stats = {}
    if ARG("file") is not None and ARG("name") is not None:
        try:
            with open(ARG("file"), 'r') as f:
                data_json = json.load(f)
            
            cb_test = ARG("name")
            if "entries" in data_json and cb_test in data_json["entries"]:
                benchmark_runs = data_json["entries"][cb_test]
                case_times = {}
                for run in benchmark_runs:
                    if "benches" not in run:
                        continue
                    for bench in run["benches"]:
                        case_name = bench.get("name")
                        grind_value = bench.get("value")
                        if case_name is None or grind_value is None:
                            continue
                        if case_name not in case_times:
                            case_times[case_name] = []
                        case_times[case_name].append(grind_value)
                for case_name, values in case_times.items():
                    if len(values) > 0:
                        avg = sum(values) / len(values)
                        cb_stats[case_name] = {"avg": avg, "count": len(values)}
                                
                cons.print(f"[bold]Loaded cb data for test: [bold]{cb_test}[/bold] ({len(cb_stats)} cases)[/bold]")
            else:
                cons.print(f"[bold yellow]Warning[/bold yellow]: Test '[bold]{cb_test}[/bold]' not found in data file.")
        except Exception as e:
            cons.print(f"[bold yellow]Warning[/bold yellow]: Could not load data file: {e}")

    table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
    table.add_column("[bold]Case[/bold]",    justify="left")
    table.add_column("[bold]Pre Process[/bold]", justify="right")
    table.add_column("[bold]Simulation[/bold]", justify="right")
    table.add_column("[bold]Post Process[/bold]", justify="right")
    if cb_stats:
        table.add_column("[bold] CB (Grind)[/bold]", justify="right")

    err = 0
    for slug in slugs:
        lhs_summary, rhs_summary = lhs["cases"][slug]["output_summary"], rhs["cases"][slug]["output_summary"]
        speedups = ['N/A', 'N/A', 'N/A']
        grind_comparison = 'N/A'

        for i, target in enumerate(sorted(DEFAULT_TARGETS, key=lambda t: t.runOrder)):
            if (target.name not in lhs_summary) or (target.name not in rhs_summary):
                cons.print(f"{target.name} not present in lhs_summary or rhs_summary - Case: {slug}")
                err = 1; continue

            if not math.isfinite(lhs_summary[target.name]["exec"]) or not math.isfinite(rhs_summary[target.name]["exec"]):
                err = 1
                cons.print(f"lhs_summary or rhs_summary reports non-real exec time for {target.name} - Case: {slug}")
            try:
                exec_time_value = lhs_summary[target.name]["exec"] / rhs_summary[target.name]["exec"]
                if exec_time_value < 0.9:
                    cons.print(f"[bold yellow]Warning[/bold yellow]: Exec time speedup for {target.name} is less than 0.9 - Case: {slug}")
                speedups[i] = f"Exec: {exec_time_value:.2f}"
                if target == SIMULATION:
                    if not math.isfinite(lhs_summary[target.name]["grind"]) or not math.isfinite(rhs_summary[target.name]["grind"]):
                        err = 1
                        cons.print(f"lhs_summary or rhs_summary reports non-real grind time for {target.name} - Case: {slug}")

                    grind_time_value = lhs_summary[target.name]["grind"] / rhs_summary[target.name]["grind"]
                    speedups[i] += f" & Grind: {grind_time_value:.2f}"
                    if grind_time_value <0.95:
                        cons.print(f"[bold red]Error[/bold red]: Benchmarking failed since grind time speedup for {target.name} below acceptable threshold (<0.95) - Case: {slug}")
                        err = 1
                    if slug in cb_stats:
                        rhs_grind = rhs_summary[target.name]["grind"]
                        stats = cb_stats[slug]
                        avg = stats["avg"]
                        offset_pct = ((rhs_grind - avg) / avg * 100)
                        color = "red" if (offset_pct) > 0 else "yellow" if abs(offset_pct) == 0 else "green"
                        grind_comparison = f"[{color}]{offset_pct:+.2f}%[/{color}] (avg: {avg:.2f})"

            except Exception as _:
                pass

        row = [f"[magenta]{slug}[/magenta]", *speedups]
        if cb_stats:
            row.append(grind_comparison)
        table.add_row(*row)

    cons.raw.print(table)
    if err:
        raise MFCException("Benchmarking failed")
        