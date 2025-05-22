#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --job-name="${name}"
#SBATCH --output="${name}.out"
#SBATCH --time=${walltime}
#SBATCH --cpus-per-task=7
% if gpu:
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
% endif
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
cd "${MFC_ROOT_DIR}"
% if engine == 'batch':
. ./mfc.sh load -c f -m ${'g' if gpu else 'c'}
% endif
cd - > /dev/null
echo

% if gpu:
    export MPICH_GPU_SUPPORT_ENABLED=1
% else:
    export MPICH_GPU_SUPPORT_ENABLED=0
% endif

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; \
            % if target.name == 'simulation':
            ${profiler} \
        % endif
            "${target.get_install_binpath(case)}")
    % else:
        (set -x; srun \
        % if engine == 'interactive':
                --nodes ${nodes} --ntasks-per-node ${tasks_per_node} \
                --cpus-per-task 7                                    \
            % if gpu:
                --gpus-per-task 1 --gpu-bind closest                 \
            % endif
        % endif
        % if target.name == 'simulation':
                ${profiler} \
        % endif
                "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
