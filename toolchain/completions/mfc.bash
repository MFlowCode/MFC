#!/usr/bin/env bash
# Bash completion for MFC (mfc.sh)
#
# To enable, add one of the following to your ~/.bashrc:
#   source /path/to/MFC/toolchain/completions/mfc.bash
# Or:
#   eval "$(./mfc.sh completion bash)"

_mfc_completions() {
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Main commands
    local commands="build run test clean validate init count packer load lint format spelling interactive completion"

    # Build targets
    local targets="pre_process simulation post_process"

    # Templates for init
    local templates="1D_minimal 2D_minimal 3D_minimal"

    # GPU options
    local gpu_options="acc mp"

    # Engines for run
    local engines="interactive batch"

    # First argument - complete commands
    if [[ ${COMP_CWORD} -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
        return 0
    fi

    # Get the command (first argument after mfc.sh)
    local command="${COMP_WORDS[1]}"

    case "${command}" in
        build)
            case "${prev}" in
                -t|--targets)
                    COMPREPLY=( $(compgen -W "${targets}" -- "${cur}") )
                    return 0
                    ;;
                -j|--jobs)
                    COMPREPLY=( $(compgen -W "1 2 4 8 16 32" -- "${cur}") )
                    return 0
                    ;;
                --gpu)
                    COMPREPLY=( $(compgen -W "${gpu_options}" -- "${cur}") )
                    return 0
                    ;;
                -i|--input)
                    COMPREPLY=( $(compgen -f -X '!*.py' -- "${cur}") )
                    return 0
                    ;;
                *)
                    local build_opts="-t --targets -j --jobs -v --verbose -d --debug-log --mpi --no-mpi --gpu --no-gpu --debug --no-debug -i --input --case-optimization"
                    COMPREPLY=( $(compgen -W "${build_opts}" -- "${cur}") )
                    return 0
                    ;;
            esac
            ;;
        run)
            case "${prev}" in
                -n|--nodes|-N|--tasks-per-node|-w|--walltime|-a|--account|-p|--partition|-e|--engine)
                    if [[ "${prev}" == "-e" || "${prev}" == "--engine" ]]; then
                        COMPREPLY=( $(compgen -W "${engines}" -- "${cur}") )
                    fi
                    return 0
                    ;;
                *)
                    # Complete case files
                    if [[ "${cur}" == -* ]]; then
                        local run_opts="-n --nodes -N --tasks-per-node -w --walltime -a --account -p --partition -e --engine -t --targets --mpi --no-mpi --gpu --no-gpu"
                        COMPREPLY=( $(compgen -W "${run_opts}" -- "${cur}") )
                    else
                        COMPREPLY=( $(compgen -f -X '!*.py' -- "${cur}") )
                    fi
                    return 0
                    ;;
            esac
            ;;
        test)
            case "${prev}" in
                -j|--jobs)
                    COMPREPLY=( $(compgen -W "1 2 4 8 16 32" -- "${cur}") )
                    return 0
                    ;;
                --only|--from|--to)
                    # Would need to dynamically get test UUIDs - leave empty for now
                    return 0
                    ;;
                *)
                    local test_opts="-j --jobs --only --from --to --generate --mpi --no-mpi --gpu --no-gpu"
                    COMPREPLY=( $(compgen -W "${test_opts}" -- "${cur}") )
                    return 0
                    ;;
            esac
            ;;
        validate)
            # Complete .py files
            if [[ "${cur}" == -* ]]; then
                local validate_opts="-d --debug-log"
                COMPREPLY=( $(compgen -W "${validate_opts}" -- "${cur}") )
            else
                COMPREPLY=( $(compgen -f -X '!*.py' -- "${cur}") )
            fi
            return 0
            ;;
        init)
            case "${prev}" in
                -t|--template)
                    COMPREPLY=( $(compgen -W "${templates}" -- "${cur}") )
                    return 0
                    ;;
                *)
                    if [[ "${cur}" == -* ]]; then
                        local init_opts="-t --template -l --list"
                        COMPREPLY=( $(compgen -W "${init_opts}" -- "${cur}") )
                    else
                        # Complete directories for case name
                        COMPREPLY=( $(compgen -d -- "${cur}") )
                    fi
                    return 0
                    ;;
            esac
            ;;
        packer)
            case "${prev}" in
                packer)
                    COMPREPLY=( $(compgen -W "pack compare" -- "${cur}") )
                    return 0
                    ;;
                pack)
                    COMPREPLY=( $(compgen -f -X '!*.py' -- "${cur}") )
                    return 0
                    ;;
                compare)
                    COMPREPLY=( $(compgen -f -X '!*.pack' -- "${cur}") )
                    return 0
                    ;;
            esac
            ;;
        load)
            case "${prev}" in
                -c)
                    # Cluster options
                    COMPREPLY=( $(compgen -W "p f s a" -- "${cur}") )
                    return 0
                    ;;
                -m)
                    # Mode options
                    COMPREPLY=( $(compgen -W "c g" -- "${cur}") )
                    return 0
                    ;;
                *)
                    COMPREPLY=( $(compgen -W "-c -m" -- "${cur}") )
                    return 0
                    ;;
            esac
            ;;
        clean|count|interactive)
            # No additional arguments
            return 0
            ;;
        completion)
            case "${prev}" in
                completion)
                    COMPREPLY=( $(compgen -W "install uninstall status" -- "${cur}") )
                    return 0
                    ;;
                install)
                    COMPREPLY=( $(compgen -W "bash zsh" -- "${cur}") )
                    return 0
                    ;;
                *)
                    return 0
                    ;;
            esac
            ;;
    esac

    return 0
}

# Register the completion function
complete -F _mfc_completions ./mfc.sh
complete -F _mfc_completions mfc.sh
complete -F _mfc_completions mfc
