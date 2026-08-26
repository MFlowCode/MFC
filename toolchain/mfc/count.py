import glob
import os
import typing

import rich.table

from .common import MFC_ROOT_DIR, file_write, format_list_to_string
from .printer import cons
from .state import ARG


def handle_dir(mfc_dir: str, srcdirname: str) -> typing.Tuple[typing.Dict[str, int], int]:
    files = {}
    total = 0

    # Untrusted tree in CI: a symlink must not reach off it.
    root = os.path.realpath(os.path.join(mfc_dir, "src")) + os.sep

    for filepath in glob.glob(os.path.join(mfc_dir, "src", srcdirname, "**", "*.*f*"), recursive=True):
        if not os.path.isfile(filepath) or not os.path.realpath(filepath).startswith(root):
            continue

        with open(filepath, errors="replace") as f:
            counter = sum(1 for line in f if line.strip() and not line.lstrip().startswith("!"))

        files[os.path.relpath(filepath, mfc_dir)] = counter
        total += counter

    return (files, total)


def count():
    target_str_list = format_list_to_string(ARG("targets"), "magenta")

    cons.print(f"[bold]Counting lines of code in {target_str_list}[/bold] (excluding whitespace and comment lines)")
    cons.indent()

    total = 0
    for codedir in ["common"] + ARG("targets"):
        dirfiles, dircount = handle_dir(MFC_ROOT_DIR, codedir)
        table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
        table.add_column(f"File (in [magenta]{codedir}[/magenta])", justify="left")
        table.add_column(f"Lines ([cyan]{dircount}[/cyan])", justify="right")

        for filepath, n in sorted(dirfiles.items()):
            table.add_row(filepath, f"[bold cyan]{n}[/bold cyan]")

        total += dircount

        cons.raw.print(table)

    cons.print(f"[bold]Total {target_str_list} lines: [bold cyan]{total}[/bold cyan].[/bold]")
    cons.print()
    cons.unindent()


def _write_markdown(filepath: str, files: list, dirs: list, total: int, total_diff: int):
    # Empty file: the CI commenter posts nothing.
    if not files:
        file_write(filepath, "")
        return

    lines = ["### Lines of Code", "", "| File | Lines | Diff |", "| :--- | ---: | ---: |"]
    lines += [f"| `{path}` | {n} | {diff:+d} |" for path, n, diff in files]
    lines += ["", "| Directory | Lines | Diff |", "| :--- | ---: | ---: |"]
    lines += [f"| {codedir} | {n} | {diff:+d} |" for codedir, n, diff in dirs if diff != 0]
    lines += [f"| **total** | **{total}** | **{total_diff:+d}** |", ""]

    file_write(filepath, "\n".join(lines))


def count_diff():
    base_dir, pr_dir = ARG("base"), ARG("pr")
    target_str_list = format_list_to_string(ARG("targets"), "magenta")

    cons.print(f"[bold]Counting lines of code in {target_str_list}[/bold] (excluding whitespace and comment lines)")
    cons.indent()

    files, dirs, total, total_diff = [], [], 0, 0
    for codedir in ["common"] + ARG("targets"):
        base_files, base_count = handle_dir(base_dir, codedir)
        pr_files, pr_count = handle_dir(pr_dir, codedir)

        for filepath in set(base_files) | set(pr_files):
            diff = pr_files.get(filepath, 0) - base_files.get(filepath, 0)
            if diff != 0:
                files.append((filepath, pr_files.get(filepath, 0), diff))

        dirs.append((codedir, pr_count, pr_count - base_count))
        total, total_diff = total + pr_count, total_diff + pr_count - base_count

    files.sort(key=lambda entry: (-abs(entry[2]), entry[0]))

    table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
    table.add_column("Changed file", justify="left")
    table.add_column("Lines", justify="right")
    table.add_column("Diff", justify="right")
    for path, n, diff in files:
        table.add_row(path, f"[bold cyan]{n}[/bold cyan]", f"[bold {'red' if diff > 0 else 'green'}]{diff:+d}[/bold {'red' if diff > 0 else 'green'}]")
    cons.raw.print(table)

    cons.print(f"[bold]Total {target_str_list} lines: [bold cyan]{total}[/bold cyan] ([bold cyan]{total_diff:+d}[/bold cyan]).[/bold]")
    cons.print()
    cons.unindent()

    if ARG("markdown") is not None:
        _write_markdown(ARG("markdown"), files, dirs, total, total_diff)
