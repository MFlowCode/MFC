#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --cpus-per-task=1
#SBATCH --job-name="${name}"
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --exclude=aswin-01
###SBATCH --nodelist=compute-3-03
% if account:
#SBATCH --account="${account}"
% endif
% if gpu:
#SBATCH --gpu-bind=verbose,closest
#SBATCH --gres=gpu:v100-16:${tasks_per_node}
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

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        (set -x; ${profiler}                              \
	    srun --mpi=pmi2     \
                     ${' '.join([f"'{x}'" for x in ARG('--') ])} \
                     "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
