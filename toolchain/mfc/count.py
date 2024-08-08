import os, glob, typing, typing
import rich.table

from .state   import ARG
from .common  import MFC_ROOT_DIR, format_list_to_string, MFCException
from .printer import cons

def handle_dir(mfc_dir: str, srcdirname: str) -> typing.Tuple[typing.Dict[str, int], int]:
    files = {}
    total = 0

    for filepath in glob.glob(os.path.join(mfc_dir, 'src', srcdirname, '*.*f*')):
        with open(filepath) as f:
            counter = 0
            for l in f.read().split('\n'):
                # Skip whitespace
                if l.isspace() or len(l) == 0:
                    continue
                # Skip comments but not !$acc ones!
                if l.lstrip().startswith("!") and not l.lstrip().startswith("!$acc"):
                    continue
                counter += 1

            files[os.path.relpath(filepath, mfc_dir)] = counter
            total += counter

    return (files, total)

def count():
    target_str_list = format_list_to_string(ARG('targets'), 'magenta')

    cons.print(f"[bold]Counting lines of code in {target_str_list}[/bold] (excluding whitespace lines)")
    cons.indent()

    total = 0
    for codedir in ['common'] + ARG("targets"):
        dirfiles, dircount = handle_dir(MFC_ROOT_DIR, codedir)
        table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
        table.add_column(f"File (in [magenta]{codedir}[/magenta])", justify="left")
        table.add_column(f"Lines ([cyan]{dircount}[/cyan])", justify="right")

        for filepath, n in dirfiles.items():
            table.add_row(os.path.basename(filepath), f"[bold cyan]{n}[/bold cyan]")

        total += dircount

        cons.raw.print(table)

    cons.print(f"[bold]Total {target_str_list} lines: [bold cyan]{total}[/bold cyan].[/bold]")
    cons.print()
    cons.unindent()

# pylint: disable=too-many-locals
def count_diff():
    target_str_list = format_list_to_string(ARG('targets'), 'magenta')
    cons.print(f"[bold]Counting lines of code in {target_str_list}[/bold] (excluding whitespace lines)")
    cons.indent()

    total = 0
    MFC_COMPARE_DIR=os.getenv('MFC_PR')
    if MFC_COMPARE_DIR is None:
        raise MFCException("MFC_PR is not in your environment.")

    print('compare dir', MFC_COMPARE_DIR)

    # MFC_COMPARE_DIR="/Users/spencer/Downloads/MFC-shbfork"
    for codedir in ['common'] + ARG("targets"):
        dirfiles_root, dircount_root = handle_dir(MFC_ROOT_DIR, codedir)
        dirfiles_pr, dircount_pr = handle_dir(MFC_COMPARE_DIR, codedir)
        table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
        table.add_column(f"File (in [magenta]{codedir}[/magenta])", justify="left")
        table.add_column(f"Lines [HEAD] ([cyan]{dircount_root}[/cyan])", justify="right")
        table.add_column(f"Lines [PR] ([cyan]{dircount_pr}[/cyan])", justify="right")
        table.add_column("", justify="right")
        table.add_column("Diff", justify="right")

        for filepath in set(dirfiles_root.keys()) | set(dirfiles_pr.keys()):
            dirfiles_root[filepath] = dirfiles_root.get(filepath, 0)
            dirfiles_pr[filepath] = dirfiles_pr.get(filepath, 0)

            PLUS  = "++ "
            MINUS = "-- "

            diff_count = dirfiles_pr[filepath] - dirfiles_root[filepath]
            mycolor = "red" if diff_count > 0 else "green"
            mysymbol = PLUS if diff_count > 0 else MINUS
            table.add_row(os.path.basename(filepath),
                          f"[bold cyan]{dirfiles_root[filepath]}[/bold cyan]",
                          f"[bold cyan]{dirfiles_pr[filepath]}[/bold cyan]",
                          mysymbol,
                          f"[bold {mycolor}]{diff_count}[/bold {mycolor}]")

        total += dircount_root

        cons.raw.print(table)

    cons.print(f"[bold]Total {target_str_list} lines: [bold cyan]{total}[/bold cyan].[/bold]")
    cons.print()
    cons.unindent()
