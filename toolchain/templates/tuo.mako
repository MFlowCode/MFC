#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
# flux: -N ${nodes}
# flux: -n ${tasks_per_node*nodes}
# flux: --job-name="${name}"
# flux: --output="${name}.out"
# flux: --error="${name}.err"
# flux: --time=${walltime}
# flux: --exclusive
# flux: --setattr=thp=always
% if account:
# flux: --bank=${account}
% endif
% if partition:
# flux: --queue=${partition}
% endif
% endif

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOT_DIR}"
% if engine == 'batch':
. ./mfc.sh load -c tuo -m ${'g' if gpu else 'c'}
% endif
cd - > /dev/null
echo

% if gpu:
    export MPICH_GPU_SUPPORT_ENABLED=1
% else:
    export MPICH_GPU_SUPPORT_ENABLED=0
% endif

export HSA_XNACK=0

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        (set -x; flux run \
            --nodes=${nodes} --ntasks=${tasks_per_node * nodes} \
            -o spindle.level=off --exclusive \
            % if gpu:
                --gpus-per-task 1 \
            % endif
            ${profiler} "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
