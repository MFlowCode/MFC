#!/usr/bin/env python3

import signal, getpass, platform, itertools

# Only import what's needed for startup - other modules are loaded lazily
from mfc         import args, lock, state
from mfc.state   import ARG
from mfc.common  import MFC_LOGO, MFCException, quit, format_list_to_string, does_command_exist, setup_debug_logging
from mfc.printer import cons

def __print_greeting():
    MFC_LOGO_LINES       = MFC_LOGO.splitlines()
    max_logo_line_length = max(len(line) for line in MFC_LOGO_LINES)

    host_line    = f"{getpass.getuser()}@{platform.node()} [{platform.system()}]"
    targets_line = f"[bold]--targets {format_list_to_string(ARG('targets'), 'magenta', 'None')}[/bold]"
    help_line    = "$ ./mfc.sh (build, run, test, clean, count, packer) --help"

    MFC_SIDEBAR_LINES = [
        f"[bold]{host_line}[/bold]",
        '-' * len(host_line),
        '',
        f"[bold]--jobs [magenta]{ARG('jobs')}[/magenta][/bold]",
        f"[bold]{' '.join(state.gCFG.make_options())}[/bold]",
        targets_line if ARG("command") != "test" else "",
        '',
        '-' * len(help_line),
        f"[yellow]{help_line}[/yellow]",
    ]

    for a, b in itertools.zip_longest(MFC_LOGO_LINES, MFC_SIDEBAR_LINES):
        a = a or ''
        lhs = a.ljust(max_logo_line_length)
        rhs = b or ''
        cons.print(
            f"[bold]{lhs}[/bold] | {rhs}",
            highlight=False
        )

    cons.print()


def __checks():
    if not does_command_exist("cmake"):
        raise MFCException("CMake is required to build MFC but couldn't be located on your system. Please ensure it installed and discoverable (e.g in your system's $PATH).")


def __run():
    # Lazy import modules only when needed for the specific command
    cmd = ARG("command")

    if cmd == "test":
        from mfc.test import test  # pylint: disable=import-outside-toplevel
        test.test()
    elif cmd == "run":
        from mfc.run import run  # pylint: disable=import-outside-toplevel
        run.run()
    elif cmd == "build":
        from mfc import build  # pylint: disable=import-outside-toplevel
        build.build()
    elif cmd == "bench":
        from mfc import bench  # pylint: disable=import-outside-toplevel
        bench.bench()
    elif cmd == "bench_diff":
        from mfc import bench  # pylint: disable=import-outside-toplevel
        bench.diff()
    elif cmd == "count":
        from mfc import count  # pylint: disable=import-outside-toplevel
        count.count()
    elif cmd == "count_diff":
        from mfc import count  # pylint: disable=import-outside-toplevel
        count.count_diff()
    elif cmd == "clean":
        from mfc import clean  # pylint: disable=import-outside-toplevel
        clean.clean()
    elif cmd == "packer":
        from mfc.packer import packer  # pylint: disable=import-outside-toplevel
        packer.packer()
    elif cmd == "validate":
        from mfc import validate  # pylint: disable=import-outside-toplevel
        validate.validate()
    elif cmd == "new":
        from mfc import init  # pylint: disable=import-outside-toplevel
        init.init()
    elif cmd == "interactive":
        from mfc.user_guide import interactive_mode  # pylint: disable=import-outside-toplevel
        interactive_mode()
    elif cmd == "completion":
        from mfc import completion  # pylint: disable=import-outside-toplevel
        completion.completion()


if __name__ == "__main__":
    try:
        lock.init()
        state.gARG = args.parse(state.gCFG)

        # Setup debug logging if requested
        setup_debug_logging(ARG("debug_log"))

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
