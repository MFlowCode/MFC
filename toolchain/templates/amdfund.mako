#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --job-name="${name}"
#SBATCH --output="${name}.out"
#SBATCH --error="${name}.err"
#SBATCH --time=${walltime}
## hpcfund does not track GPUs through Slurm (GresTypes=null): the MI250X nodes
## hand out every GCD with the node, so there is no --gres / --gpu-bind to set.
#SBATCH --partition=${partition or 'mi2508x'}
% if account:
#SBATCH --account="${account}"
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
cd "${MFC_ROOT_DIR}"
. ./mfc.sh load -c amdfund -m ${'g' if gpu_enabled else 'c'}
cd - > /dev/null
echo

% if gpu_enabled:
# Abort loudly if OpenMP offload cannot reach a GPU (e.g. run on a login node).
export OMP_TARGET_OFFLOAD=MANDATORY
% endif

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
