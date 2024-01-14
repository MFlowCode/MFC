#!/usr/bin/env bash

% if engine == 'batch':
#BSUB -J {{{name}}}
#BSUB -nnodes {{{nodes}}}
#BSUB -W {{{walltime[:-3]}}}
#BSUB -N
% if account:
#BSUB -P {{{account}}}
% endif
% endif

<%include file="prologue.mako"/>

ok ":) Loading modules:\n"
cd "${rootdir}"
. ./mfc.sh load -c s -m ${'g' if gpu else 'c'}
cd - > /dev/null
echo

% for binpath in binpaths:
    ok ":) Running ${binpath}:\n"

    % if not mpi:
        ${' '.join([f"'{x}'" for x in profiler ])} "${binpath}"
    % else:
        ${' '.join([f"'{x}'" for x in profiler ])}          \
            jsrun                                           \
                ${'--smpiargs="-gpu"' if gpu else ''}       \
                --nrs          ${tasks_per_node*nodes}      \
                --cpu_per_rs   1                            \
                --gpu_per_rs   ${1 if gpu else 0}           \
                --tasks_per_rs 1                            \
                ${' '.join([f"'{x}'" for x in ARG('--') ])} \
                "${binpath}"
    % endif

    % if engine == 'interactive':
        code=$?
        if [ $code -ne 0 ]; then
            echo
            error ":( $MAGENTA${binpath}$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
            echo
            exit 1
        fi
    % endif

    echo
% endfor

<%include file="epilogue.mako"/>
