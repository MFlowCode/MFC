#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --job-name="${name}"
#SBATCH --output="${name}.out"
#SBATCH --time=${walltime}
% if account:
#SBATCH --account=${account}
% endif
% if partition:
#SBATCH --partition=${partition}
% endif
% if quality_of_service:
#SBATCH --qos=${quality_of_service}
% endif
% if email:
#SBATCH --mail-user=${email}
#SBATCH --mail-type="BEGIN, END, FAIL"
% endif
% endif

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOTDIR}"
% if engine == 'batch':
. ./mfc.sh load -c f -m ${'g' if gpu else 'c'}
% endif
cd - > /dev/null
echo

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${' '.join([f"'{x}'" for x in profiler ])} "${target.get_install_binpath()}")
    % else:
        (set -x; ${' '.join([f"'{x}'" for x in profiler ])}    \
            srun -N ${nodes}                 \
                   -n ${tasks_per_node}                              \
                   ${' '.join([f"'{x}'" for x in ARG('--') ])} \
                   "${target.get_install_binpath()}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
