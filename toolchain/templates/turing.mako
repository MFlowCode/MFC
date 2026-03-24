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
#SBATCH --gres=gpu:1
#SBATCH --mem=208G
#SBATCH -C "A30|A100"
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
. ./mfc.sh load -c t -m ${'g' if gpu_enabled else 'c'}
cd - > /dev/null
echo

% if gpu_enabled:
    export LD_LIBRARY_PATH=/cm/shared/spack/opt/spack/linux-ubuntu20.04-x86_64/gcc-13.2.0/cuda-12.3.0-vuydybqum6mloi2vvov7yn2juaurmtao/lib64:$LD_LIBRARY_PATH 
% endif

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        % if gpu_enabled:
            (set -x; ${profiler}                                    \
                srun --gres=gpu:1 -C "A30|A100"                     \
                $MPI_HOME/bin/mpirun --np ${nodes*tasks_per_node}   \
                "${target.get_install_binpath(case)}")
        % else:
            (set -x; ${profiler}                                    \
                srun  --mpi=pmi2   --ntasks=${nodes*tasks_per_node} \
                "${target.get_install_binpath(case)}")
        % endif
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
