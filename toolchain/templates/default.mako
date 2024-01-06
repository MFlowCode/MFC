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
echo

% for binpath in binpaths:
    echo -e ":) Running $MAGENTA${binpath}$COLOR_RESET:\n"

    % if not mpi:
        ${' '.join([f"'{x}'" for x in profiler ])} "${binpath}"
    % else:
        if command -v jsrun > /dev/null; then
            jsrun --nrs          ${tasks_per_node*nodes}     \
                  --cpu_per_rs   1                           \
                  --gpu_per_rs   ${1 if gpu else 0}          \
                  --tasks_per_rs 1                           \
                  ${' '.join([f"'{x}'" for x in profiler ])} \
                  "${binpath}"
        elif command -v srun > /dev/null; then
            srun --ntasks-per-node ${tasks_per_node}         \
                 ${' '.join([f"'{x}'" for x in profiler ])}  \
                 "${binpath}"
        elif command -v mpirun > /dev/null; then
            mpirun -np ${nodes*tasks_per_node}                \
                   ${' '.join([f"'{x}'" for x in profiler ])} \
                   "${binpath}"
        else
            echo -e "\n:( Could not find a suitable MPI launcher.\n"
            exit 1
        fi
    % endif

    code=$?
    if [ $code -ne 0 ]; then
        echo -e "\n:( $MAGENTA${binpath}$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET.\n"
        exit 1
    fi

    echo
% endfor

<%include file="epilogue.mako"/>
