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
% if gpu_enabled:
#SBATCH --gres=gpu:V100:${tasks_per_node}
#SBATCH --mem-per-gpu=16G\
% endif
% if email:
#SBATCH --mail-user=${email}
#SBATCH --mail-type="BEGIN, END, FAIL"
% endif
% endif

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOT_DIR}"
. ./mfc.sh load -c p -m ${'g' if gpu_enabled else 'c'}
cd - > /dev/null
echo

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        mkdir -p /storage/scratch1/6/sbryngelson3/mytmp
        chmod 777 /storage/scratch1/6/sbryngelson3/mytmp
        (set -x; ${profiler}    \
            mpirun --mca orte_tmpdir_base /storage/scratch1/6/sbryngelson3/mytmp \
                   --np ${nodes*tasks_per_node}           \
                   --bind-to none                         \
                   "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
