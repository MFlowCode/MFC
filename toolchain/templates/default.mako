#!/usr/bin/env bash

. "${rootdir}/toolchain/util.sh"

% if engine == 'batch':
    error "The$MAGENTA default$COLOR_RESET template does not support batch jobs. Please use a different template via the $MAGENTA--computer$COLOR_RESET option.\n"
    exit 1
% endif

<%include file="prologue.mako"/>

warn "This is the$MAGENTA default$COLOR_RESET template."
warn "It is not intended to support all systems and execution engines."
warn "Please use a different template via the $MAGENTA--computer$COLOR_RESET option."

% if mpi:
    # Find a suitable MPI launcher and store it in the variable "binary".
    for binary in ${binary or ''} jsrun srun mpirun mpiexec; do
        if command -v $binary > /dev/null; then
            break
        fi
    done

    if ! command -v $binary > /dev/null; then
        error ":( Could not find a suitable MPI launcher.\n"
        exit 1
    else
        ok ":) Selected MPI launcher $MAGENTA$binary$COLOR_RESET. Use$MAGENTA --binary$COLOR_RESET to override."
    fi
% endif

% for binpath in binpaths:
    ok ":) Running $MAGENTA${binpath}$COLOR_RESET:\n"

    % if not mpi:
        ${' '.join([f"'{x}'" for x in profiler ])} "${binpath}"
    % else:
        if [ "$binary" == "jsrun" ]; then
            ${' '.join([f"'{x}'" for x in profiler ])}            \
                jsrun --nrs          ${tasks_per_node*nodes}      \
                      --cpu_per_rs   1                            \
                      --gpu_per_rs   ${1 if gpu else 0}           \
                      --tasks_per_rs 1                            \
                      ${' '.join([f"'{x}'" for x in ARG('--') ])} \
                      "${binpath}"
        elif [ "$binary" == "srun" ]; then
            ${' '.join([f"'{x}'" for x in profiler ])}           \
                srun --ntasks-per-node ${tasks_per_node}         \
                     ${' '.join([f"'{x}'" for x in ARG('--') ])} \
                     "${binpath}"
        elif [ "$binary" == "mpirun" ] || [ "$binary" == "mpiexec" ]; then
            ${' '.join([f"'{x}'" for x in profiler ])}              \
                $binary -np ${nodes*tasks_per_node}                 \
                        ${' '.join([f"'{x}'" for x in ARG('--') ])} \
                        "${binpath}"
        fi
    % endif

    code=$?
    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA${binpath}$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi

    echo
% endfor

<%include file="epilogue.mako"/>
