#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#PBS -l select=${nodes}:ncpus=192:mpiprocs=${tasks_per_node}
#PBS -N "${name}"
#PBS -l walltime=${walltime}
% if partition:
#PBS -q ${partition}
% endif
% if account:
#PBS -A ${account}
% endif
% if email:
#PBS -M ${email}
#PBS -m abe
% endif
#PBS -o "${name}.out"
#PBS -e "${name}.err"
#PBS -V
% endif

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOT_DIR}"
. ./mfc.sh load -c c -m ${'g' if gpu_enabled else 'c'}
cd - > /dev/null
echo

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        (set -x; ${profiler}                              \
            mpirun -np ${nodes*tasks_per_node}            \
                   "${target.get_install_binpath(case)}")
    % endif
            
    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
