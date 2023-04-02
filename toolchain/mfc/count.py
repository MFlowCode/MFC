import os, glob, typing, typing

from .common  import MFC_ROOTDIR
from .printer import cons

import rich.table

def handle_dir(dirpath: str) -> typing.Tuple[typing.List[typing.Tuple[str, int]], int]:
    files = []
    total = 0

    for filepath in glob.glob(os.path.join(dirpath, '*.*f*')):
        with open(filepath) as f:
            count = sum(1 if not l.isspace() else 0 for l in f.read().split('\n'))
            files.append((filepath, count))
            total += count

    files.sort(key=lambda x: x[1], reverse=True)

    return (files, total)

def count():
    cons.print("[bold]Counting lines of code in [magenta]MFC[/magenta][/bold] (excluding whitespace lines)")
    cons.print()
    cons.indent()
  
    total = 0
    for codedir in ['common', 'pre_process', 'simulation', 'post_process']:
        dirfiles, dircount = handle_dir(os.path.join(MFC_ROOTDIR, 'src', codedir))        
        table = rich.table.Table(show_header=True, box=rich.table.box.SIMPLE)
        table.add_column(f"File (in [magenta]{codedir}[/magenta])", justify="left")
        table.add_column(f"Lines ([cyan]{dircount}[/cyan])", justify="right")
        
        for filepath, count in dirfiles:
            table.add_row(f"{os.path.basename(filepath)}", f"[bold cyan]{count}[/bold cyan]")
        
        total += dircount

        cons.raw.print(table)

    cons.print(f"[bold]Total MFC lines: [bold cyan]{total}[/bold cyan].[/bold]")
    cons.print()
    cons.unindent()

