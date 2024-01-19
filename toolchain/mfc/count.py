import os, glob, typing, typing
import rich.table

from .state   import ARG
from .common  import MFC_ROOTDIR, format_list_to_string, MFCException
from .printer import cons

def handle_dir(dirpath: str) -> typing.Tuple[typing.List[typing.Tuple[str, int]], int]:
    files = []
    total = 0

    for filepath in glob.glob(os.path.join(dirpath, '*.*f*')):
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

            files.append((filepath, counter))
            total += counter

    files.sort(key=lambda x: x[1], reverse=True)

    return (files, total)

def count():
    target_str_list = format_list_to_string(ARG('targets'), 'magenta')

    cons.print(f"[bold]Counting lines of code in {target_str_list}[/bold] (excluding whitespace lines)")
    cons.indent()

    total = 0
    for codedir in ['common'] + ARG("targets"):
        dirfiles, dircount = handle_dir(os.path.join(MFC_ROOTDIR, 'src', codedir))
        table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
        table.add_column(f"File (in [magenta]{codedir}[/magenta])", justify="left")
        table.add_column(f"Lines ([cyan]{dircount}[/cyan])", justify="right")

        for filepath, n in dirfiles:
            table.add_row(f"{os.path.basename(filepath)}", f"[bold cyan]{n}[/bold cyan]")

        total += dircount

        cons.raw.print(table)

    cons.print(f"[bold]Total {target_str_list} lines: [bold cyan]{total}[/bold cyan].[/bold]")
    cons.print()
    cons.unindent()

def count_diff():
    target_str_list = format_list_to_string(ARG('targets'), 'magenta')
    cons.print(f"[bold]Counting lines of code in {target_str_list}[/bold] (excluding whitespace lines)")
    cons.indent()

    total = 0
    MFC_COMPAREDIR=os.getenv('MFC_PR')
    if MFC_COMPAREDIR is None: raise MFCException("MFC_PR is not in your enviroment.")

    # MFC_COMPAREDIR="/Users/spencer/Downloads/MFC-shbfork"
    for codedir in ['common'] + ARG("targets"):
        dirfiles_root, dircount = handle_dir(os.path.join(MFC_ROOTDIR, 'src', codedir))
        dirfiles_pr, dircount_pr = handle_dir(os.path.join(MFC_COMPAREDIR, 'src', codedir))
        table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
        table.add_column(f"File (in [magenta]{codedir}[/magenta])", justify="left")
        table.add_column(f"Lines [MASTER] ([cyan]{dircount}[/cyan])", justify="right")
        table.add_column(f"Lines [INCOMING] ([cyan]{dircount_pr}[/cyan])", justify="right")
        table.add_column("Diff.", justify="right")

        ii = 0
        files_pr = [os.path.basename(dirfiles_pr[i][0]) for i in range(len(dirfiles_pr))]
        files_root = [os.path.basename(dirfiles_root[i][0]) for i in range(len(dirfiles_root))]

        for filepath, n in dirfiles_pr:
            for filepath_root, q in dirfiles_root:
                if os.path.basename(dirfiles_pr[ii][0]) == os.path.basename(filepath_root):
                    diff_count = n - dirfiles_root[ii][1]
                    mycolor = "red" if diff_count > 0 else "green"
                    table.add_row(f"{os.path.basename(filepath)}", f"[bold cyan]{n}[/bold cyan]", f"[bold cyan]{dirfiles_root[ii][1]}[/bold cyan]", f"[bold {mycolor}]{diff_count}[/bold {mycolor}]")

            if files_pr[ii] not in files_root:
                table.add_row(f"{os.path.basename(files_pr[ii])}", "----", f"[bold green]{n}[/bold green]", f"[bold green]{n}[/bold green]")
            ii += 1

        total += dircount

        cons.raw.print(table)

    cons.print(f"[bold]Total {target_str_list} lines: [bold cyan]{total}[/bold cyan].[/bold]")
    cons.print()
    cons.unindent()
