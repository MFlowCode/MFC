import os, json, time, typing, datetime, subprocess

import rich.table

from .printer import cons
from .state   import ARG
from .build   import build_targets
from .common  import system, MFC_SUBDIR
from .        import sched

import glob

def count():
    
    cons.print("[bold]Counting lines of code in [magenta]MFC[/magenta][/bold] (including blank ones)")
    cons.indent()
  
    dirs = ['common', 'pre_process', 'simulation', 'post_process']
    for mydir in dirs:
        cons.print(f"[bold]In [magenta]{mydir}[/magenta]:[/bold]")
        dir_path = f'./src/{mydir}/*.*f*'
        test_list = glob.glob(dir_path)
        result = str(' '.join(test_list))
        os.system("wc -l " + result + " | sort | sed \$d")
        os.system("wc -l " + result + " | sort | tail -n 1")
        cons.print()

    cons.print("[bold]Total MFC lines:[/bold]")
    dir_path = f'./src/*/*.*f*'
    test_list = glob.glob(dir_path)
    result = str(' '.join(test_list))
    os.system("wc -l " + result + " | sort | tail -n 1")




