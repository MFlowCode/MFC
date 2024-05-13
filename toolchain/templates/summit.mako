#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#BSUB -J {{{name}}}
#BSUB -nnodes {{{nodes}}}
#BSUB -W {{{walltime[:-3]}}}
#BSUB -N
% if account:
#BSUB -P {{{account}}}
% endif
% endif

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOTDIR}"
. ./mfc.sh load -c s -m ${'g' if gpu else 'c'}
cd - > /dev/null
echo

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${' '.join([f"'{x}'" for x in profiler ])} "${target.get_install_binpath(case)}")
    % else:
        (set -x; ${' '.join([f"'{x}'" for x in profiler ])} \
            jsrun                                           \
                ${'--smpiargs="-gpu"' if gpu else ''}       \
                --nrs          ${tasks_per_node*nodes}      \
                --cpu_per_rs   1                            \
                --gpu_per_rs   ${1 if gpu else 0}           \
                --tasks_per_rs 1                            \
                ${' '.join([f"'{x}'" for x in ARG('--') ])} \
                "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
