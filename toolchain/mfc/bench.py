import os, sys, uuid, subprocess, dataclasses, typing, math
import rich.table
from .printer import cons
from .state import ARG, CFG
from .build import get_targets, DEFAULT_TARGETS, SIMULATION
from .common import system, MFC_BENCH_FILEPATH, MFC_BUILD_DIR, format_list_to_string
from .common import file_load_yaml, file_dump_yaml, create_directory
from .common import MFCException

SINGLE_PRECISION_SPEEDUP_THRESHOLD = 1.25  # Minimum speedup for PR double vs single precision

@dataclasses.dataclass
class BenchCase:
    slug: str
    path: str
    args: typing.List[str]

def bench(targets=None, precision="double"):
    """
    Benchmarks the provided targets in the specified precision mode (single or double).
    """
    if targets is None:
        targets = ARG("targets")

    if precision not in ["single", "double"]:
        raise ValueError("Precision must be 'single' or 'double'.")

    # Set precision_flag based on precision
    if precision == "single":
        precision_flag = ["--single"]
    else:
        precision_flag = []  # No flag needed for double precision

    targets = get_targets(targets)
    bench_dirpath = os.path.join(MFC_BUILD_DIR, "benchmarks", str(uuid.uuid4())[:4])
    create_directory(bench_dirpath)

    cons.print(f"[bold]Benchmarking {format_list_to_string(ARG('targets'), 'magenta')} in {precision} precision "
               f"([magenta]{os.path.relpath(bench_dirpath)}[/magenta]):[/bold]")
    cons.indent()
    cons.print()

    CASES = [BenchCase(**case) for case in file_load_yaml(MFC_BENCH_FILEPATH)]

    for case in CASES:
        case.args = case.args + ARG("--")
        case.path = os.path.abspath(case.path)

    results = {
        "metadata": {
            "invocation": sys.argv[1:],
            "lock": dataclasses.asdict(CFG()),
            "precision": precision,
        },
        "cases": {},
    }

    for i, case in enumerate(CASES):
        summary_filepath = os.path.join(bench_dirpath, f"{case.slug}-{precision}.yaml")
        log_filepath = os.path.join(bench_dirpath, f"{case.slug}-{precision}.out")

        cons.print(f"{str(i + 1).zfill(len(str(len(CASES))))}/{len(CASES)}: {case.slug} @ [bold]{os.path.relpath(case.path)}[/bold]")
        cons.indent()
        cons.print(f"> Log:     [bold]{os.path.relpath(log_filepath)}[/bold]")
        cons.print(f"> Summary: [bold]{os.path.relpath(summary_filepath)}[/bold]")

        with open(log_filepath, "w") as log_file:
            system(
                ["./mfc.sh", "run", case.path, "--case-optimization"] +
                ["--targets"] + [t.name for t in targets] +
                ["--output-summary", summary_filepath] +
                case.args +
                precision_flag,  # Use the precision_flag here
                stdout=log_file,
                stderr=subprocess.STDOUT
            )

        results["cases"][case.slug] = {
            "description": dataclasses.asdict(case),
            "output_summary": file_load_yaml(summary_filepath),
        }
        cons.unindent()

    file_dump_yaml(ARG("output"), results)
    cons.print(f"Wrote results to [bold magenta]{os.path.relpath(ARG('output'))}[/bold magenta].")
    cons.unindent()

def diff():
    """
    Compares the results between two benchmark YAML files (lhs vs rhs).
    Checks both PR vs master and PR single vs PR double precision.
    """
    lhs, rhs = file_load_yaml(ARG("lhs")), file_load_yaml(ARG("rhs"))

    lhs_precision = lhs["metadata"].get("precision", "double")
    rhs_precision = rhs["metadata"].get("precision", "double")

    is_pr_single_vs_double = lhs_precision == "double" and rhs_precision == "single"

    cons.print(f"[bold]Comparing Benchmarks: Speedups from [magenta]{os.path.relpath(ARG('lhs'))}[/magenta] to "
               f"[magenta]{os.path.relpath(ARG('rhs'))}[/magenta][/bold]")

    if lhs["metadata"] != rhs["metadata"] and not is_pr_single_vs_double:
        def _lock_to_str(lock):
            return ' '.join([f"{k}={v}" for k, v in lock.items()])

        cons.print(f"""\
[bold yellow]Warning[/bold yellow]: Metadata in lhs and rhs are not equal.
    lhs:
    * Invocation: [magenta]{' '.join(lhs['metadata']['invocation'])}[/magenta]
    * Modes:      {_lock_to_str(lhs['metadata']['lock'])}
    * Precision:  {lhs_precision}
    rhs:
    * Invocation: {' '.join(rhs['metadata']['invocation'])}
    * Modes:      [magenta]{_lock_to_str(rhs['metadata']['lock'])}[/magenta]
    * Precision:  {rhs_precision}
        """)

    slugs = set(lhs["cases"].keys()) & set(rhs["cases"].keys())
    if len(slugs) not in [len(lhs["cases"]), len(rhs["cases"])]:
        cons.print(f"""\
[bold yellow]Warning[/bold yellow]: Cases in lhs and rhs are not equal.
    Using intersection: {slugs} with {len(slugs)} elements.""")

    table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
    table.add_column("[bold]Case[/bold]", justify="left")
    table.add_column("[bold]Speedup (Exec)[/bold]", justify="right")
    table.add_column("[bold]Speedup (Grind)[/bold]", justify="right")

    err = 0

    for slug in slugs:
        lhs_summary = lhs["cases"][slug]["output_summary"]
        rhs_summary = rhs["cases"][slug]["output_summary"]

        try:
            exec_speedup = lhs_summary["exec"] / rhs_summary["exec"]
            grind_speedup = lhs_summary["grind"] / rhs_summary["grind"]

            if is_pr_single_vs_double and exec_speedup < SINGLE_PRECISION_SPEEDUP_THRESHOLD:
                cons.print(f"[bold red]Error[/bold red]: Case {slug} failed speedup requirement: "
                           f"Exec speedup {exec_speedup:.2f} < {SINGLE_PRECISION_SPEEDUP_THRESHOLD}.")
                err = 1

            table.add_row(slug, f"{exec_speedup:.2f}x", f"{grind_speedup:.2f}x")

        except KeyError as e:
            table.add_row(slug, "Error", "Error")
            cons.print(f"[bold yellow]Warning[/bold yellow]: Missing key {e} for case {slug}.")
        except ZeroDivisionError:
            table.add_row(slug, "Inf", "Inf")
            cons.print(f"[bold yellow]Warning[/bold yellow]: Zero execution time in case {slug}.")

    cons.raw.print(table)

    if err:
        raise MFCException("Benchmarking failed: Some cases did not meet the performance requirements.")

