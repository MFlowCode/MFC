"""
Generate shell completion scripts from CLI schema.

This module generates bash and zsh completion scripts that are always
in sync with the CLI schema definitions.
"""

from typing import List, Set
from .schema import CLISchema, Command, CompletionType


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
    elif comp_type == CompletionType.FILES_PY:
        return 'COMPREPLY=( $(compgen -f -X "!*.py" -- "${cur}") $(compgen -d -- "${cur}") )'
    elif comp_type == CompletionType.FILES_PACK:
        return 'COMPREPLY=( $(compgen -f -X "!*.pack" -- "${cur}") $(compgen -d -- "${cur}") )'
    elif comp_type == CompletionType.FILES_YAML:
        return 'COMPREPLY=( $(compgen -f -X "!*.y*ml" -- "${cur}") $(compgen -d -- "${cur}") )'
    elif comp_type == CompletionType.FILES:
        return 'COMPREPLY=( $(compgen -f -- "${cur}") $(compgen -d -- "${cur}") )'
    elif comp_type == CompletionType.DIRECTORIES:
        return 'COMPREPLY=( $(compgen -d -- "${cur}") )'
    return ""


def generate_bash_completion(schema: CLISchema) -> str:
    """Generate bash completion script from schema."""
    # Collect all command names and aliases
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

    # Generate case for each command
    for cmd in schema.commands:
        if not cmd.arguments and not cmd.positionals and not cmd.include_common and not cmd.subcommands:
            continue  # Skip simple commands with no completable arguments

        # Include aliases in case pattern
        patterns = [cmd.name] + cmd.aliases
        lines.append(f'        {"|".join(patterns)})')

        # Collect options for this command
        options = _collect_all_options(cmd, schema)

        # Handle subcommands (like packer pack, packer compare)
        if cmd.subcommands:
            lines.append('            if [[ ${COMP_CWORD} -eq 2 ]]; then')
            subcmd_names = [sc.name for sc in cmd.subcommands]
            lines.append(f'                COMPREPLY=( $(compgen -W "{" ".join(subcmd_names)}" -- "${{cur}}") )')
            lines.append('                return 0')
            lines.append('            fi')
            lines.append('            ;;')
            continue

        # Generate prev-based completion for arguments with choices
        has_prev_cases = False
        all_args = list(cmd.arguments)
        for meg in cmd.mutually_exclusive:
            all_args.extend(meg.arguments)
        # Add common set arguments
        for set_name in cmd.include_common:
            common_set = schema.get_common_set(set_name)
            if common_set and not common_set.mfc_config_flags:
                all_args.extend(common_set.arguments)

        for arg in all_args:
            if arg.completion.type not in (CompletionType.CHOICES, CompletionType.FILES_PY,
                                           CompletionType.FILES_PACK, CompletionType.FILES,
                                           CompletionType.DIRECTORIES, CompletionType.FILES_YAML):
                continue

            if not has_prev_cases:
                lines.append('            case "${prev}" in')
                has_prev_cases = True

            flags = []
            if arg.short:
                flags.append(f'-{arg.short}')
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

        # Default completion - options or positional files
        if options:
            lines.append(f'            local opts="{" ".join(options)}"')
            lines.append('            if [[ "${cur}" == -* ]]; then')
            lines.append('                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )')

            # Check for file positionals
            if cmd.positionals:
                pos = cmd.positionals[0]
                if pos.completion.type != CompletionType.NONE:
                    lines.append('            else')
                    comp_choices = pos.completion.choices or pos.choices
                    completion_code = _bash_completion_for_type(pos.completion.type, comp_choices)
                    if completion_code:
                        lines.append(f'                {completion_code}')

            lines.append('            fi')
        elif cmd.positionals:
            # Only positionals, no options
            pos = cmd.positionals[0]
            if pos.completion.type != CompletionType.NONE:
                comp_choices = pos.completion.choices or pos.choices
                completion_code = _bash_completion_for_type(pos.completion.type, comp_choices)
                if completion_code:
                    lines.append(f'            {completion_code}')

        lines.append('            return 0')
        lines.append('            ;;')

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
        # Handle aliases by matching command name
        all_names = [cmd.name] + cmd.aliases
        for name in all_names:
            lines.append(f'                {name})')

            arg_lines = []

            # Positionals
            for i, pos in enumerate(cmd.positionals):
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

                help_text = pos.help.replace("'", "").replace("[", "").replace("]", "")[:40]
                arg_lines.append(f"'{i+1}:{help_text}{completion}'")

            # Options from common sets
            for set_name in cmd.include_common:
                common_set = schema.get_common_set(set_name)
                if common_set is None:
                    continue

                if common_set.mfc_config_flags:
                    # Add MFC config flags
                    arg_lines.append("'--mpi[Enable MPI]'")
                    arg_lines.append("'--no-mpi[Disable MPI]'")
                    arg_lines.append("'--gpu[Enable GPU]:mode:(acc mp)'")
                    arg_lines.append("'--no-gpu[Disable GPU]'")
                    arg_lines.append("'--debug[Enable debug mode]'")
                    arg_lines.append("'--no-debug[Disable debug mode]'")
                else:
                    for arg in common_set.arguments:
                        desc = arg.help.replace("'", "").replace("[", "").replace("]", "")[:40]
                        completion = ""
                        if arg.completion.type == CompletionType.CHOICES:
                            choices = arg.completion.choices or arg.choices or []
                            completion = f':({" ".join(str(c) for c in choices)})'

                        if arg.short:
                            arg_lines.append(f"'-{arg.short}[{desc}]{completion}'")
                        arg_lines.append(f"'--{arg.name}[{desc}]{completion}'")

            # Command-specific arguments
            all_args = list(cmd.arguments)
            for meg in cmd.mutually_exclusive:
                all_args.extend(meg.arguments)

            for arg in all_args:
                desc = arg.help.replace("'", "").replace("[", "").replace("]", "")[:40]
                completion = ""
                if arg.completion.type == CompletionType.CHOICES:
                    choices = arg.completion.choices or arg.choices or []
                    completion = f':({" ".join(str(c) for c in choices)})'
                elif arg.completion.type == CompletionType.FILES_PY:
                    completion = ':_files -g "*.py"'
                elif arg.completion.type == CompletionType.FILES_PACK:
                    completion = ':_files -g "*.pack"'
                elif arg.completion.type == CompletionType.FILES:
                    completion = ':_files'
                elif arg.completion.type == CompletionType.DIRECTORIES:
                    completion = ':_files -/'

                if arg.short:
                    arg_lines.append(f"'-{arg.short}[{desc}]{completion}'")
                arg_lines.append(f"'--{arg.name}[{desc}]{completion}'")

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
