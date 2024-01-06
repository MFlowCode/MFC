#!/usr/bin/env bash

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
% if gpu:
#SBATCH --gres=gpu:V100:${tasks_per_node}
#SBATCH --mem-per-gpu=16G\
% endif
% if email:
#SBATCH --mail-user=${email}
#SBATCH --mail-type="BEGIN, END, FAIL"
% endif
% endif

<%include file="prologue.mako"/>

echo -e ":) Loading modules:\n"
cd "${rootdir}" && . ./mfc.sh load -c p -m ${'g' if gpu else 'c'} && cd -
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
