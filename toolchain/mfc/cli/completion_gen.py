"""
Generate shell completion scripts from CLI schema.

This module generates bash and zsh completion scripts that are always
in sync with the CLI schema definitions.
"""

from typing import List, Set
from .schema import CLISchema, Command, CompletionType


# Mapping of completion types to bash completion expressions
_BASH_COMPLETION_MAP = {
    CompletionType.FILES_PY: 'COMPREPLY=( $(compgen -f -X "!*.py" -- "${cur}") $(compgen -d -- "${cur}") )',
    CompletionType.FILES_PACK: 'COMPREPLY=( $(compgen -f -X "!*.pack" -- "${cur}") $(compgen -d -- "${cur}") )',
    CompletionType.FILES_YAML: 'COMPREPLY=( $(compgen -f -X "!*.y*ml" -- "${cur}") $(compgen -d -- "${cur}") )',
    CompletionType.FILES: 'COMPREPLY=( $(compgen -f -- "${cur}") $(compgen -d -- "${cur}") )',
    CompletionType.DIRECTORIES: 'COMPREPLY=( $(compgen -d -- "${cur}") )',
}


def _collect_all_options(cmd: Command, schema: CLISchema) -> List[str]:
    """Collect all option flags for a command including common sets."""
    options: Set[str] = set()

    # Common arguments
    for set_name in cmd.include_common:
        common_set = schema.get_common_set(set_name)
        if common_set is None:
            continue

        # MFC config flags
        if common_set.mfc_config_flags:
            options.update([
                "--mpi", "--no-mpi",
                "--gpu", "--no-gpu",
                "--debug", "--no-debug",
                "--gcov", "--no-gcov",
                "--unified", "--no-unified",
                "--single", "--no-single",
                "--mixed", "--no-mixed",
                "--fastmath", "--no-fastmath",
            ])
        else:
            for arg in common_set.arguments:
                if arg.short:
                    options.add(f"-{arg.short}")
                options.add(f"--{arg.name}")

    # Command-specific arguments
    for arg in cmd.arguments:
        if arg.short:
            options.add(f"-{arg.short}")
        options.add(f"--{arg.name}")

    # Mutually exclusive arguments
    for meg in cmd.mutually_exclusive:
        for arg in meg.arguments:
            if arg.short:
                options.add(f"-{arg.short}")
            options.add(f"--{arg.name}")

    return sorted(options)


def _bash_completion_for_type(comp_type: CompletionType, choices: List[str] = None) -> str:
    """Generate bash completion expression for a completion type."""
    if comp_type == CompletionType.CHOICES and choices:
        return f'COMPREPLY=( $(compgen -W "{" ".join(choices)}" -- "${{cur}}") )'
    return _BASH_COMPLETION_MAP.get(comp_type, "")


def _collect_all_args(cmd: Command, schema: CLISchema) -> list:
    """Collect all arguments for a command including common sets."""
    all_args = list(cmd.arguments)
    for meg in cmd.mutually_exclusive:
        all_args.extend(meg.arguments)
    for set_name in cmd.include_common:
        common_set = schema.get_common_set(set_name)
        if common_set and not common_set.mfc_config_flags:
            all_args.extend(common_set.arguments)
    return all_args


def _generate_bash_prev_cases(cmd: Command, schema: CLISchema) -> List[str]:
    """Generate bash prev-based completion cases for a command."""
    lines = []
    has_prev_cases = False
    completable_types = (CompletionType.CHOICES, CompletionType.FILES_PY,
                         CompletionType.FILES_PACK, CompletionType.FILES,
                         CompletionType.DIRECTORIES, CompletionType.FILES_YAML)

    for arg in _collect_all_args(cmd, schema):
        if arg.completion.type not in completable_types:
            continue

        if not has_prev_cases:
            lines.append('            case "${prev}" in')
            has_prev_cases = True

        flags = [f'-{arg.short}'] if arg.short else []
        flags.append(f'--{arg.name}')

        lines.append(f'                {"|".join(flags)})')
        comp_choices = arg.completion.choices or arg.choices
        completion_code = _bash_completion_for_type(arg.completion.type, comp_choices)
        if completion_code:
            lines.append(f'                    {completion_code}')
        lines.append('                    return 0')
        lines.append('                    ;;')

    if has_prev_cases:
        lines.append('            esac')

    return lines


def _generate_bash_command_case(cmd: Command, schema: CLISchema) -> List[str]:
    """Generate bash completion case for a single command."""
    lines = []

    # Include aliases in case pattern
    patterns = [cmd.name] + cmd.aliases
    lines.append(f'        {"|".join(patterns)})')

    options = _collect_all_options(cmd, schema)

    # Handle subcommands (like packer pack, packer compare)
    if cmd.subcommands:
        lines.append('            if [[ ${COMP_CWORD} -eq 2 ]]; then')
        subcmd_names = [sc.name for sc in cmd.subcommands]
        lines.append(f'                COMPREPLY=( $(compgen -W "{" ".join(subcmd_names)}" -- "${{cur}}") )')
        lines.append('                return 0')
        lines.append('            fi')
        lines.append('            ;;')
        return lines

    # Generate prev-based completion
    lines.extend(_generate_bash_prev_cases(cmd, schema))

    # Default completion - options or positional files
    if options:
        lines.append(f'            local opts="{" ".join(options)}"')
        lines.append('            if [[ "${cur}" == -* ]]; then')
        lines.append('                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )')

        if cmd.positionals and cmd.positionals[0].completion.type != CompletionType.NONE:
            lines.append('            else')
            pos = cmd.positionals[0]
            comp_choices = pos.completion.choices or pos.choices
            completion_code = _bash_completion_for_type(pos.completion.type, comp_choices)
            if completion_code:
                lines.append(f'                {completion_code}')

        lines.append('            fi')
    elif cmd.positionals and cmd.positionals[0].completion.type != CompletionType.NONE:
        pos = cmd.positionals[0]
        comp_choices = pos.completion.choices or pos.choices
        completion_code = _bash_completion_for_type(pos.completion.type, comp_choices)
        if completion_code:
            lines.append(f'            {completion_code}')

    lines.append('            return 0')
    lines.append('            ;;')
    return lines


def generate_bash_completion(schema: CLISchema) -> str:
    """Generate bash completion script from schema."""
    commands = schema.get_all_command_names()

    lines = [
        '#!/usr/bin/env bash',
        '# AUTO-GENERATED from cli/commands.py - Do not edit manually',
        '# Regenerate with: ./mfc.sh generate',
        '',
        '_mfc_completions() {',
        '    local cur prev command',
        '    COMPREPLY=()',
        '    cur="${COMP_WORDS[COMP_CWORD]}"',
        '    prev="${COMP_WORDS[COMP_CWORD-1]}"',
        '',
        f'    local commands="{" ".join(sorted(commands))}"',
        '',
        '    # First argument - complete commands',
        '    if [[ ${COMP_CWORD} -eq 1 ]]; then',
        '        COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )',
        '        return 0',
        '    fi',
        '',
        '    local command="${COMP_WORDS[1]}"',
        '',
        '    case "${command}" in',
    ]

    for cmd in schema.commands:
        if not cmd.arguments and not cmd.positionals and not cmd.include_common and not cmd.subcommands:
            continue
        lines.extend(_generate_bash_command_case(cmd, schema))

    lines.extend([
        '    esac',
        '',
        '    return 0',
        '}',
        '',
        'complete -o filenames -o bashdefault -F _mfc_completions ./mfc.sh',
        'complete -o filenames -o bashdefault -F _mfc_completions mfc.sh',
        'complete -o filenames -o bashdefault -F _mfc_completions mfc',
    ])

    return '\n'.join(lines)


def _zsh_completion_for_positional(pos, index: int) -> str:
    """Generate zsh completion spec for a positional argument."""
    completion = ""
    if pos.completion.type == CompletionType.FILES_PY:
        completion = ':_files -g "*.py"'
    elif pos.completion.type == CompletionType.FILES_PACK:
        completion = ':_files -g "*.pack"'
    elif pos.completion.type == CompletionType.CHOICES:
        choices = pos.completion.choices or pos.choices or []
        completion = f':({" ".join(choices)})'
    elif pos.completion.type == CompletionType.DIRECTORIES:
        completion = ':_files -/'
    elif pos.completion.type == CompletionType.FILES:
        completion = ':_files'

    help_text = pos.help.replace("'", "").replace("[", "").replace("]", "")[:120]
    return f"'{index}:{help_text}{completion}'"


def _zsh_completion_for_arg(arg) -> str:
    """Generate zsh completion suffix for an argument."""
    if arg.completion.type == CompletionType.CHOICES:
        choices = arg.completion.choices or arg.choices or []
        return f':({" ".join(str(c) for c in choices)})'
    if arg.completion.type == CompletionType.FILES_PY:
        return ':_files -g "*.py"'
    if arg.completion.type == CompletionType.FILES_PACK:
        return ':_files -g "*.pack"'
    if arg.completion.type == CompletionType.FILES:
        return ':_files'
    if arg.completion.type == CompletionType.DIRECTORIES:
        return ':_files -/'
    return ""


def _generate_zsh_command_args(cmd: Command, schema: CLISchema) -> List[str]:
    """Generate zsh argument lines for a command."""
    arg_lines = []

    # Positionals
    for i, pos in enumerate(cmd.positionals):
        arg_lines.append(_zsh_completion_for_positional(pos, i + 1))

    # Options from common sets
    for set_name in cmd.include_common:
        common_set = schema.get_common_set(set_name)
        if common_set is None:
            continue

        if common_set.mfc_config_flags:
            arg_lines.extend([
                "'--mpi[Enable MPI]'",
                "'--no-mpi[Disable MPI]'",
                "'--gpu[Enable GPU]:mode:(acc mp)'",
                "'--no-gpu[Disable GPU]'",
                "'--debug[Enable debug mode]'",
                "'--no-debug[Disable debug mode]'",
            ])
        else:
            for arg in common_set.arguments:
                desc = arg.help.replace("'", "").replace("[", "").replace("]", "")[:120]
                completion = _zsh_completion_for_arg(arg)
                if arg.short:
                    arg_lines.append(f"'-{arg.short}[{desc}]{completion}'")
                arg_lines.append(f"'--{arg.name}[{desc}]{completion}'")

    # Command-specific arguments
    all_args = list(cmd.arguments)
    for meg in cmd.mutually_exclusive:
        all_args.extend(meg.arguments)

    for arg in all_args:
        desc = arg.help.replace("'", "").replace("[", "").replace("]", "")[:120]
        completion = _zsh_completion_for_arg(arg)
        if arg.short:
            arg_lines.append(f"'-{arg.short}[{desc}]{completion}'")
        arg_lines.append(f"'--{arg.name}[{desc}]{completion}'")

    return arg_lines


def generate_zsh_completion(schema: CLISchema) -> str:
    """Generate zsh completion script from schema."""
    lines = [
        '#compdef mfc.sh ./mfc.sh mfc',
        '# AUTO-GENERATED from cli/commands.py - Do not edit manually',
        '# Regenerate with: ./mfc.sh generate',
        '',
        '_mfc() {',
        '    local context state state_descr line',
        '    typeset -A opt_args',
        '',
        '    local -a commands',
        '    commands=(',
    ]

    # Commands with descriptions
    for cmd in schema.commands:
        desc = cmd.help.replace("'", "\\'").replace('"', '\\"')
        lines.append(f'        "{cmd.name}:{desc}"')
        for alias in cmd.aliases:
            lines.append(f'        "{alias}:Alias for {cmd.name}"')

    lines.extend([
        '    )',
        '',
        '    _arguments -C \\',
        "        '1: :->command' \\",
        "        '*:: :->args'",
        '',
        '    case $state in',
        '        command)',
        "            _describe -t commands 'mfc command' commands",
        '            ;;',
        '        args)',
        '            case $words[1] in',
    ])

    # Generate case for each command
    for cmd in schema.commands:
        all_names = [cmd.name] + cmd.aliases
        for name in all_names:
            lines.append(f'                {name})')
            arg_lines = _generate_zsh_command_args(cmd, schema)
            if arg_lines:
                lines.append('                    _arguments \\')
                lines.append('                        ' + ' \\\n                        '.join(arg_lines))
            lines.append('                    ;;')

    lines.extend([
        '            esac',
        '            ;;',
        '    esac',
        '}',
        '',
        '_mfc "$@"',
    ])

    return '\n'.join(lines)
