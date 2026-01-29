#!/usr/bin/env bash
# AUTO-GENERATED from cli/commands.py - Do not edit manually
# Regenerate with: ./mfc.sh generate

_mfc_completions() {
    local cur prev command
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    local commands="b bench bench_diff build c clean completion count count_diff format generate help interactive lint load new packer r run spelling t test v validate"

    # First argument - complete commands
    if [[ ${COMP_CWORD} -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
        return 0
    fi

    local command="${COMP_WORDS[1]}"

    case "${command}" in
        build|b)
            case "${prev}" in
                -i|--input)
                    COMPREPLY=( $(compgen -f -X "!*.py" -- "${cur}") $(compgen -d -- "${cur}") )
                    return 0
                    ;;
                -t|--targets)
                    COMPREPLY=( $(compgen -W "fftw hdf5 silo lapack hipfort pre_process simulation post_process syscheck documentation" -- "${cur}") )
                    return 0
                    ;;
            esac
            local opts="--case-optimization --debug --debug-log --fastmath --gcov --gpu --input --jobs --mixed --mpi --no-debug --no-fastmath --no-gcov --no-gpu --no-mixed --no-mpi --no-single --no-unified --single --targets --unified --verbose -d -i -j -t -v"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            fi
            return 0
            ;;
        run|r)
            case "${prev}" in
                -e|--engine)
                    COMPREPLY=( $(compgen -W "interactive batch" -- "${cur}") )
                    return 0
                    ;;
                -b|--binary)
                    COMPREPLY=( $(compgen -W "mpirun jsrun srun mpiexec" -- "${cur}") )
                    return 0
                    ;;
                -c|--computer)
                    COMPREPLY=( $(compgen -W "bridges2 carpenter carpenter-cray default delta deltaai frontier hipergator nautilus oscar phoenix phoenix-bench santis tuo" -- "${cur}") )
                    return 0
                    ;;
                -t|--targets)
                    COMPREPLY=( $(compgen -W "fftw hdf5 silo lapack hipfort pre_process simulation post_process syscheck documentation" -- "${cur}") )
                    return 0
                    ;;
            esac
            local opts="-# --account --binary --case-optimization --clean --computer --debug --debug-log --dry-run --email --engine --fastmath --gcov --gpu --gpus --jobs --mixed --mpi --name --ncu --no-build --no-debug --no-fastmath --no-gcov --no-gpu --no-mixed --no-mpi --no-single --no-unified --nodes --nsys --output-summary --partition --quality_of_service --rcu --rsys --scratch --single --targets --tasks-per-node --unified --verbose --wait --walltime -@ -N -a -b -c -d -e -g -j -n -o -p -q -s -t -v -w"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            else
                COMPREPLY=( $(compgen -f -X "!*.py" -- "${cur}") $(compgen -d -- "${cur}") )
            fi
            return 0
            ;;
        test|t)
            local opts="-% --add-new-variables --case-optimization --debug --debug-log --dry-run --fastmath --from --gcov --generate --gpu --gpus --jobs --list --max-attempts --mixed --mpi --no-build --no-debug --no-examples --no-fastmath --no-gcov --no-gpu --no-mixed --no-mpi --no-single --no-unified --only --percent --rdma-mpi --remove-old-tests --single --test-all --to --unified --verbose -a -d -f -g -j -l -m -o -t -v"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            fi
            return 0
            ;;
        clean|c)
            case "${prev}" in
                -t|--targets)
                    COMPREPLY=( $(compgen -W "fftw hdf5 silo lapack hipfort pre_process simulation post_process syscheck documentation" -- "${cur}") )
                    return 0
                    ;;
            esac
            local opts="--debug --debug-log --fastmath --gcov --gpu --jobs --mixed --mpi --no-debug --no-fastmath --no-gcov --no-gpu --no-mixed --no-mpi --no-single --no-unified --single --targets --unified --verbose -d -j -t -v"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            fi
            return 0
            ;;
        validate|v)
            local opts="--debug-log -d"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            else
                COMPREPLY=( $(compgen -f -X "!*.py" -- "${cur}") $(compgen -d -- "${cur}") )
            fi
            return 0
            ;;
        new)
            case "${prev}" in
                -t|--template)
                    COMPREPLY=( $(compgen -W "1D_minimal 2D_minimal 3D_minimal" -- "${cur}") )
                    return 0
                    ;;
            esac
            local opts="--list --template -l -t"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            else
                COMPREPLY=( $(compgen -d -- "${cur}") )
            fi
            return 0
            ;;
        packer)
            if [[ ${COMP_CWORD} -eq 2 ]]; then
                COMPREPLY=( $(compgen -W "pack compare" -- "${cur}") )
                return 0
            fi
            ;;
        completion)
            COMPREPLY=( $(compgen -W "install uninstall status" -- "${cur}") )
            return 0
            ;;
        help)
            COMPREPLY=( $(compgen -W "gpu clusters batch debugging" -- "${cur}") )
            return 0
            ;;
        generate)
            local opts="--check"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            fi
            return 0
            ;;
        bench)
            case "${prev}" in
                -t|--targets)
                    COMPREPLY=( $(compgen -W "fftw hdf5 silo lapack hipfort pre_process simulation post_process syscheck documentation" -- "${cur}") )
                    return 0
                    ;;
            esac
            local opts="--debug --debug-log --fastmath --gcov --gpu --gpus --jobs --mem --mixed --mpi --no-debug --no-fastmath --no-gcov --no-gpu --no-mixed --no-mpi --no-single --no-unified --output --single --targets --unified --verbose -d -g -j -m -o -t -v"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            fi
            return 0
            ;;
        bench_diff)
            local opts="--debug --debug-log --fastmath --gcov --gpu --jobs --mixed --mpi --no-debug --no-fastmath --no-gcov --no-gpu --no-mixed --no-mpi --no-single --no-unified --single --unified --verbose -d -j -v"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            fi
            return 0
            ;;
        count)
            case "${prev}" in
                -t|--targets)
                    COMPREPLY=( $(compgen -W "fftw hdf5 silo lapack hipfort pre_process simulation post_process syscheck documentation" -- "${cur}") )
                    return 0
                    ;;
            esac
            local opts="--debug --debug-log --fastmath --gcov --gpu --jobs --mixed --mpi --no-debug --no-fastmath --no-gcov --no-gpu --no-mixed --no-mpi --no-single --no-unified --single --targets --unified --verbose -d -j -t -v"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            fi
            return 0
            ;;
        count_diff)
            case "${prev}" in
                -t|--targets)
                    COMPREPLY=( $(compgen -W "fftw hdf5 silo lapack hipfort pre_process simulation post_process syscheck documentation" -- "${cur}") )
                    return 0
                    ;;
            esac
            local opts="--debug --debug-log --fastmath --gcov --gpu --jobs --mixed --mpi --no-debug --no-fastmath --no-gcov --no-gpu --no-mixed --no-mpi --no-single --no-unified --single --targets --unified --verbose -d -j -t -v"
            if [[ "${cur}" == -* ]]; then
                COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
            fi
            return 0
            ;;
    esac

    return 0
}

complete -o filenames -o bashdefault -F _mfc_completions ./mfc.sh
complete -o filenames -o bashdefault -F _mfc_completions mfc.sh
complete -o filenames -o bashdefault -F _mfc_completions mfc