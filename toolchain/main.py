#!/usr/bin/env python3

import os, signal, getpass, platform, itertools

# Only import what's needed for startup - other modules are loaded lazily
from mfc         import args, lock, state
from mfc.state   import ARG
from mfc.common  import MFC_LOGO, MFC_ROOT_DIR, MFCException, quit, format_list_to_string, does_command_exist, setup_debug_logging
from mfc.printer import cons


def __do_regenerate(toolchain: str):
    """Perform the actual regeneration of completion scripts and schema."""
    import json  # pylint: disable=import-outside-toplevel
    from pathlib import Path  # pylint: disable=import-outside-toplevel
    from mfc.cli.commands import MFC_CLI_SCHEMA  # pylint: disable=import-outside-toplevel
    from mfc.cli.completion_gen import generate_bash_completion, generate_zsh_completion  # pylint: disable=import-outside-toplevel
    from mfc.params.generators.json_schema_gen import generate_json_schema  # pylint: disable=import-outside-toplevel

    cons.print("[dim]Auto-regenerating completion scripts...[/dim]")

    completions_dir = Path(toolchain) / "completions"
    completions_dir.mkdir(exist_ok=True)

    # Generate completion files
    (completions_dir / "mfc.bash").write_text(generate_bash_completion(MFC_CLI_SCHEMA))
    (completions_dir / "_mfc").write_text(generate_zsh_completion(MFC_CLI_SCHEMA))

    # Generate JSON schema
    schema = generate_json_schema(include_descriptions=True)
    with open(Path(toolchain) / "mfc-case-schema.json", 'w', encoding='utf-8') as f:
        json.dump(schema, f, indent=2)


def __ensure_generated_files():
    """Auto-regenerate completion scripts and docs if source files have changed."""
    toolchain = os.path.join(MFC_ROOT_DIR, "toolchain")

    # Source files that trigger regeneration
    sources = [
        os.path.join(toolchain, "mfc", "cli", "commands.py"),
        os.path.join(toolchain, "mfc", "cli", "completion_gen.py"),
        os.path.join(toolchain, "mfc", "params", "definitions.py"),
    ]

    # Generated files to check
    generated = [
        os.path.join(toolchain, "completions", "mfc.bash"),
        os.path.join(toolchain, "mfc-case-schema.json"),
    ]

    # Get max mtime of source files
    try:
        source_mtime = max(os.path.getmtime(s) for s in sources if os.path.exists(s))
    except ValueError:
        return  # No source files found, skip check

    # Check if any generated file is missing or older than sources
    needs_regen = any(
        not os.path.exists(g) or os.path.getmtime(g) < source_mtime
        for g in generated
    )

    if needs_regen:
        __do_regenerate(toolchain)

def __print_greeting():
    MFC_LOGO_LINES       = MFC_LOGO.splitlines()
    max_logo_line_length = max(len(line) for line in MFC_LOGO_LINES)

    host_line    = f"{getpass.getuser()}@{platform.node()} [{platform.system()}]"
    targets_line = f"[bold]--targets {format_list_to_string(ARG('targets'), 'magenta', 'None')}[/bold]"
    help_line    = "$ ./mfc.sh (build, run, test, clean, new, validate, params) --help"

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


def __run():  # pylint: disable=too-many-branches
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
    elif cmd == "generate":
        from mfc import generate  # pylint: disable=import-outside-toplevel
        generate.generate()
    elif cmd == "params":
        from mfc import params_cmd  # pylint: disable=import-outside-toplevel
        params_cmd.params()


if __name__ == "__main__":
    try:
        lock.init()
        state.gARG = args.parse(state.gCFG)

        # Setup debug logging if requested
        setup_debug_logging(ARG("debug_log", dflt=None))

        lock.switch(state.MFCConfig.from_dict(state.gARG))

        # Ensure IDE configuration is up to date (lightweight check)
        from mfc.ide import ensure_vscode_settings  # pylint: disable=import-outside-toplevel
        ensure_vscode_settings()

        # Auto-regenerate completion scripts if source files changed
        __ensure_generated_files()

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
