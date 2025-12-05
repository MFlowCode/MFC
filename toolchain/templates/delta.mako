#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --cpus-per-task=1
#SBATCH --job-name="${name}"
#SBATCH --time=${walltime}
% if partition:
#SBATCH --partition=${partition}
% endif
% if account:
#SBATCH --account="${account}"
% endif
% if gpu_enabled:
#SBATCH --gpus-per-node=${tasks_per_node}
#SBATCH --mem=208G
#SBATCH --gpu-bind=closest
% endif
#SBATCH --output="${name}.out"
#SBATCH --error="${name}.err"
#SBATCH --export=ALL
% if email:
#SBATCH --mail-user=${email}
#SBATCH --mail-type="BEGIN, END, FAIL"
% endif
% endif

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOT_DIR}"
. ./mfc.sh load -c d -m ${'g' if gpu_enabled else 'c'}
cd - > /dev/null
echo

# Fixes Delta not being able to find core library file
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/spack/deltas11-2023-03/apps/linux-rhel8-zen3/nvhpc-24.1/openmpi-4.1.5-zkiklxi/lib/

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        (set -x; ${profiler}                                   \
            srun --ntasks ${nodes*tasks_per_node}                 \
                   "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
