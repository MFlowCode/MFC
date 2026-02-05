"""
Generate CLI documentation from schema.

This module generates markdown documentation that is always
in sync with the CLI schema definitions.
"""

from typing import List
from .schema import CLISchema, Command, Argument


def _format_argument_row(arg: Argument) -> str:
    """Format a single argument as a markdown table row."""
    flags = []
    if arg.short:
        flags.append(f"`-{arg.short}`")
    flags.append(f"`--{arg.name}`")
    flag_str = ", ".join(flags)

    # Format default value
    default_str = _format_default(arg.default)

    # Escape pipes in help text
    help_text = arg.help.replace("|", "\\|")

    return f"| {flag_str} | {help_text} | {default_str} |"


def _format_default(default) -> str:
    """Format a default value for display."""
    if default is None:
        return "-"
    if isinstance(default, bool):
        return "`true`" if default else "`false`"
    if isinstance(default, list):
        if not default:
            return "`[]`"
        result = ", ".join(str(d) for d in default[:3])
        if len(default) > 3:
            result += "..."
        return f"`{result}`"
    # Handle empty strings
    if default == "":
        return "-"
    return f"`{default}`"


def _get_command_arguments(cmd: Command, schema: CLISchema) -> List[Argument]:
    """Get all arguments for a command including common sets."""
    all_args = []

    # Collect from common sets
    for set_name in cmd.include_common:
        common_set = schema.get_common_set(set_name)
        if common_set and not common_set.mfc_config_flags:
            all_args.extend(common_set.arguments)

    # Add command-specific arguments
    all_args.extend(cmd.arguments)

    # Add mutually exclusive arguments
    for meg in cmd.mutually_exclusive:
        all_args.extend(meg.arguments)

    return all_args


def _generate_options_table(cmd: Command, schema: CLISchema) -> List[str]:
    """Generate the options table for a command."""
    lines = []
    all_args = _get_command_arguments(cmd, schema)

    if all_args or "mfc_config" in cmd.include_common:
        lines.append("**Options:**")
        lines.append("")
        lines.append("| Option | Description | Default |")
        lines.append("|--------|-------------|---------|")

        for arg in all_args:
            lines.append(_format_argument_row(arg))

        # Add MFC config flags if included
        if "mfc_config" in cmd.include_common:
            lines.append("| `--mpi`, `--no-mpi` | Enable/disable MPI | `true` |")
            lines.append("| `--gpu [acc/mp]`, `--no-gpu` | Enable GPU (OpenACC/OpenMP) | `no` |")
            lines.append("| `--debug`, `--no-debug` | Build with debug compiler flags | `false` |")

        lines.append("")

    return lines


def _generate_subcommands_section(cmd: Command) -> List[str]:
    """Generate the subcommands section for a command."""
    lines = []
    if not cmd.subcommands:
        return lines

    lines.append("**Subcommands:**")
    lines.append("")

    for subcmd in cmd.subcommands:
        lines.append(f"#### {cmd.name} {subcmd.name}")
        lines.append("")
        lines.append(subcmd.help)
        lines.append("")

        if subcmd.positionals:
            lines.append("Arguments:")
            for pos in subcmd.positionals:
                lines.append(f"- `{pos.name.upper()}` - {pos.help}")
            lines.append("")

        if subcmd.arguments:
            lines.append("Options:")
            lines.append("")
            lines.append("| Option | Description | Default |")
            lines.append("|--------|-------------|---------|")
            for arg in subcmd.arguments:
                lines.append(_format_argument_row(arg))
            lines.append("")

    return lines


def _generate_command_section(cmd: Command, schema: CLISchema) -> List[str]:
    """Generate markdown section for a single command."""
    lines = []

    # Command header (no alias in heading to keep anchor simple)
    lines.append(f"### {cmd.name}")
    lines.append("")

    # Alias note if present
    if cmd.aliases:
        lines.append(f"**Alias:** `{cmd.aliases[0]}`")
        lines.append("")

    # Description
    description = cmd.description or cmd.help
    lines.append(description)
    lines.append("")

    # Usage
    usage_parts = ["./mfc.sh", cmd.name]
    for pos in cmd.positionals:
        if pos.nargs == "?":
            usage_parts.append(f"[{pos.name.upper()}]")
        else:
            usage_parts.append(pos.name.upper())
    usage_parts.append("[OPTIONS]")
    lines.append(f"**Usage:** `{' '.join(usage_parts)}`")
    lines.append("")

    # Positional arguments
    if cmd.positionals:
        lines.append("**Arguments:**")
        lines.append("")
        for pos in cmd.positionals:
            lines.append(f"- `{pos.name.upper()}` - {pos.help}")
        lines.append("")

    # Options table
    lines.extend(_generate_options_table(cmd, schema))

    # Examples
    if cmd.examples:
        lines.append("**Examples:**")
        lines.append("")
        lines.append("```bash")
        for example in cmd.examples:
            lines.append(f"# {example.description}")
            lines.append(example.command)
            lines.append("")
        lines.append("```")
        lines.append("")

    # Subcommands
    lines.extend(_generate_subcommands_section(cmd))

    lines.append("---")
    lines.append("")

    return lines


def _generate_commands_by_category(
    schema: CLISchema,
    category_commands: List[str],
    header: str
) -> List[str]:
    """Generate command sections for a category."""
    lines = []
    matching = [c for c in schema.commands if c.name in category_commands]
    if matching:
        lines.append(f"## {header}")
        lines.append("")
        for cmd in matching:
            lines.extend(_generate_command_section(cmd, schema))
    return lines


def generate_cli_reference(schema: CLISchema) -> str:
    """Generate complete CLI reference documentation in markdown."""
    lines = [
        "@page cli-reference Command Line Reference",
        "",
        "# Command Line Reference",
        "",
        "> **Auto-generated** from `toolchain/mfc/cli/commands.py`",
        "> ",
        "> Regenerate with: `./mfc.sh generate`",
        "",
        "## Overview",
        "",
        schema.description,
        "",
        "## Quick Reference",
        "",
        "| Command | Alias | Description |",
        "|---------|-------|-------------|",
    ]

    # Quick reference table
    for cmd in schema.commands:
        alias = f"`{cmd.aliases[0]}`" if cmd.aliases else "-"
        # Use HTML code tag inside link to avoid Doxygen markdown parsing issues
        lines.append(f"| [<code>{cmd.name}</code>](#{cmd.name}) | {alias} | {cmd.help} |")

    lines.append("")
    lines.append("## Commands")
    lines.append("")

    # Command categories
    core_commands = ["build", "run", "test", "clean", "validate"]
    utility_commands = ["new", "params", "packer", "completion", "generate", "help"]
    dev_commands = ["lint", "format", "spelling", "count", "count_diff"]
    ci_commands = ["bench", "bench_diff"]
    other_commands = ["load", "interactive"]

    # Core workflow commands first (no header, directly under Commands)
    for cmd in schema.commands:
        if cmd.name in core_commands:
            lines.extend(_generate_command_section(cmd, schema))

    # Category sections
    lines.extend(_generate_commands_by_category(schema, utility_commands, "Utility Commands"))
    lines.extend(_generate_commands_by_category(schema, dev_commands, "Development Commands"))
    lines.extend(_generate_commands_by_category(schema, ci_commands, "CI Commands"))
    lines.extend(_generate_commands_by_category(schema, other_commands, "Other Commands"))

    # Common options section
    lines.extend([
        "## Common Options",
        "",
        "Many commands share common option sets:",
        "",
        "### Target Selection (`-t, --targets`)",
        "",
        "Available targets:",
        "- `pre_process` - Pre-processor",
        "- `simulation` - Main simulation",
        "- `post_process` - Post-processor",
        "- `syscheck` - System check utility",
        "- `documentation` - Build documentation",
        "",
        "### Build Configuration Flags",
        "",
        "| Flag | Description |",
        "|------|-------------|",
        "| `--mpi` / `--no-mpi` | Enable/disable MPI support |",
        "| `--gpu [acc/mp]` / `--no-gpu` | Enable GPU with OpenACC or OpenMP |",
        "| `--debug` / `--no-debug` | Build with debug compiler flags |",
        "| `--gcov` / `--no-gcov` | Enable code coverage |",
        "| `--single` / `--no-single` | Single precision |",
        "| `--mixed` / `--no-mixed` | Mixed precision |",
        "",
        "### Verbosity (`-v, --verbose`)",
        "",
        "Controls output verbosity level:",
        "",
        "- `-v` - Basic verbose output",
        "- `-vv` - Show build commands",
        "- `-vvv` - Full verbose output including CMake details",
        "",
        "### Debug Logging (`-d, --debug-log`)",
        "",
        "Enables debug logging for the Python toolchain (mfc.sh internals).",
        "This is for troubleshooting the build system, not the MFC simulation code.",
        "",
    ])

    return "\n".join(lines)


def generate_command_summary(schema: CLISchema) -> str:
    """Generate a concise command summary for quick reference."""
    lines = [
        "# MFC Quick Reference",
        "",
        "## Commands",
        "",
    ]

    for cmd in schema.commands:
        alias_str = f" ({cmd.aliases[0]})" if cmd.aliases else ""
        lines.append(f"- **{cmd.name}**{alias_str}: {cmd.help}")

    lines.extend([
        "",
        "## Common Patterns",
        "",
        "```bash",
        "# Build MFC",
        "./mfc.sh build",
        "./mfc.sh build --gpu          # With GPU support",
        "./mfc.sh build -j 8           # Parallel build",
        "",
        "# Run a case",
        "./mfc.sh run case.py",
        "./mfc.sh run case.py -n 4     # 4 MPI ranks",
        "",
        "# Run tests",
        "./mfc.sh test",
        "./mfc.sh test -j 4            # Parallel tests",
        "",
        "# Validate a case",
        "./mfc.sh validate case.py",
        "```",
        "",
    ])

    return "\n".join(lines)
