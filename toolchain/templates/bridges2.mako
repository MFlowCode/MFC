#!/usr/bin/env bash

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

<%include file="prologue.mako"/>

echo -e ":) Loading modules:\n"
cd "${rootdir}"
. ./mfc.sh load -c b -m ${'g' if gpu else 'c'}
cd - > /dev/null
echo

% for binpath in binpaths:
    echo -e ":) Running ${binpath.split('/')[-1]}:\n"

    % if not mpi:
        ${' '.join([f"'{x}'" for x in profiler ])} "${binpath}"
    % else:
        mpirun -np ${nodes*tasks_per_node}             \
            ${' '.join([f"'{x}'" for x in profiler ])} \
            "${binpath}"
    % endif

    % if engine == 'interactive':
        code=$?
        if [ $code -ne 0 ]; then
            echo -e "\n:( $MAGENTA${binpath}$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET.\n"
            exit 1
        fi
    % endif

    echo
% endfor

<%include file="epilogue.mako"/>
