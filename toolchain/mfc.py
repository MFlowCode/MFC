#!/usr/bin/env python3

import signal, getpass, platform, itertools, dataclasses

from mfc         import args, lock, build, bench, state
from mfc.state   import ARG
from mfc.run     import run
from mfc.test    import test
from mfc.common  import MFC_LOGO, MFCException, quit, format_list_to_string, does_command_exist
from mfc.printer import cons


def __print_greeting():
    MFC_LOGO_LINES       = MFC_LOGO.splitlines()
    max_logo_line_length = max([ len(line) for line in MFC_LOGO_LINES ])

    host_line    = f"{getpass.getuser()}@{platform.node()} [{platform.system()}]"
    targets_line = f"[bold]--targets {format_list_to_string([ f'[magenta]{target}[/magenta]' for target in ARG('targets')], 'None')}[/bold]"

    MFC_SIDEBAR_LINES = [
        "",
        f"[bold]{host_line}[/bold]",
        '-' * len(host_line),
        "",
        f"[bold]--jobs [magenta]{ARG('jobs')}[/magenta][/bold]"
    ] + [
        f"[bold]--{'' if getattr(state.gCFG, field.name) else 'no-'}{field.name}[/bold]" for field in dataclasses.fields(state.gCFG)
    ] + [
        targets_line if ARG("command") != "test" else "",
        "",
        "[yellow]$ ./mfc.sh \[build, run, test, clean] --help[/yellow]",
    ]

    for a, b in itertools.zip_longest(MFC_LOGO_LINES, MFC_SIDEBAR_LINES):
        lhs = a.ljust(max_logo_line_length)
        rhs = b if b is not None else ''
        cons.print(
            f"[bold blue] {lhs} [/bold blue]  {rhs}",
            highlight=False
        )

    cons.print()


def __checks():
    if not does_command_exist("cmake"):
        raise MFCException("CMake is required to build MFC but couldn't be located on your system. Please ensure it installed and discoverable (e.g in your system's $PATH).")


def __run():    
    {"test":  test.test,   "run":   run.run,    "build": build.build,
     "clean": build.clean, "bench": bench.bench
    }[ARG("command")]()


if __name__ == "__main__":
    try:
        lock.init()
        state.gARG = args.parse(state.gCFG)

        lock.switch(state.MFCConfig.from_dict(state.gARG))
        
        __print_greeting()
        __checks()
        __run()

    except MFCException as exc:
        cons.reset()
        cons.print(f"""\


[bold red]Error[/bold red]: {str(exc)}
""")
        quit(signal.SIGTERM)
    except KeyboardInterrupt as exc:
        quit(signal.SIGTERM)
    except Exception as exc:
        cons.reset()
        cons.print_exception()
        cons.print(f"""\


[bold red]ERROR[/bold red]: An unexpected exception occurred: {str(exc)}
""")

        quit(signal.SIGTERM)
