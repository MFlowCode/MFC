import os, sys, uuid, subprocess, dataclasses, typing

import rich.table

from .printer import cons
from .state   import ARG, CFG
from .build   import get_targets, DEFAULT_TARGETS
from .common  import system, MFC_BENCH_FILEPATH, MFC_SUBDIR, format_list_to_string
from .common  import file_load_yaml, file_dump_yaml, create_directory


@dataclasses.dataclass
class BenchCase:
    slug: str
    path: str
    args: typing.List[str]


def bench(targets = None):
    if targets is None:
        targets = ARG("targets")

    targets = get_targets(targets)

    bench_dirpath = os.path.join(MFC_SUBDIR, "benchmarks", str(uuid.uuid4())[:4])
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
                ["./mfc.sh", "run", case.path, ARG('mem'), "--case-optimization"] +
                ["--targets"] + [t.name for t in targets] +
                ["--output-summary", summary_filepath] +
                case.args,
                stdout=log_file,
                stderr=subprocess.STDOUT)

        results["cases"][case.slug] = {
            "description":    dataclasses.asdict(case),
            "output_summary": file_load_yaml(summary_filepath),
        }

    file_dump_yaml(ARG("output"), results)

    cons.print(f"Wrote results to [bold magenta]{os.path.relpath(ARG('output'))}[/bold magenta].")

    cons.unindent()


def diff():
    lhs, rhs = file_load_yaml(ARG("lhs")), file_load_yaml(ARG("rhs"))

    cons.print(f"[bold]Comparing Bencharks: [magenta]{os.path.relpath(ARG('lhs'))}[/magenta] is x times slower than [magenta]{os.path.relpath(ARG('rhs'))}[/magenta].[/bold]")

    if lhs["metadata"] != rhs["metadata"]:
        def _lock_to_str(lock):
            return ' '.join([f"{k}={v}" for k, v in lock.items()])

        cons.print(f"[bold yellow]Warning[/bold yellow]: Metadata in lhs and rhs are not equal.")
        cons.print(f" This could mean that the benchmarks are not comparable (e.g. one was run on CPUs and the other on GPUs).")
        cons.print(f" lhs:")
        cons.print(f" * Invocation: [magenta]{' '.join(lhs['metadata']['invocation'])}[/magenta]")
        cons.print(f" * Modes:      {_lock_to_str(lhs['metadata']['lock'])}")
        cons.print(f" rhs:")
        cons.print(f" * Invocation: {' '.join(rhs['metadata']['invocation'])}")
        cons.print(f" * Modes:      [magenta]{_lock_to_str(rhs['metadata']['lock'])}[/magenta]")

    slugs = set(lhs["cases"].keys()) & set(rhs["cases"].keys())
    if len(slugs) not in [len(lhs["cases"]), len(rhs["cases"])]:
        cons.print(f"[bold yellow]Warning[/bold yellow]: Cases in lhs and rhs are not equal.")
        cons.print(f" * rhs cases: {', '.join(set(rhs['cases'].keys()) - slugs)}.")
        cons.print(f" * lhs cases: {', '.join(set(lhs['cases'].keys()) - slugs)}.")
        cons.print(f" Using intersection: {slugs} with {len(slugs)} elements.")

    table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
    table.add_column("[bold]Case[/bold]",    justify="left")
    table.add_column("[bold]Pre Process[/bold]", justify="right")
    table.add_column("[bold]Simulation[/bold]", justify="right")
    table.add_column("[bold]Post Process[/bold]", justify="right")

    for slug in slugs:
        lhs_summary = lhs["cases"][slug]["output_summary"]
        rhs_summary = rhs["cases"][slug]["output_summary"]

        speedups = ['N/A', 'N/A', 'N/A']

        for i, target in enumerate(sorted(DEFAULT_TARGETS, key=lambda t: t.runOrder)):
            if target.name not in lhs_summary or target.name not in rhs_summary:
                continue

            speedups[i] = f"{lhs_summary[target.name] / rhs_summary[target.name]:.2f}x"

        table.add_row(f"[magenta]{slug}[/magenta]", *speedups)

    cons.raw.print(table)
