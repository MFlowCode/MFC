#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#SBATCH --uenv=icon/25.2:v1
#SBATCH --nodes=${nodes}
#SBATCH --reservation=g183
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --job-name="${name}"
#SBATCH --output="${name}.out"
#SBATCH --error="${name}.err"
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

# NVHPC and CUDA env vars
export NV_ACC_USE_MALLOC=0                    # use cudaMallocManaged instead of malloc ( compiled using -gpu=mem:unified )
export NVCOMPILER_ACC_NO_MEMHINTS=1           # disable implicit compiler hints
#export CUDA_BUFFER_PAGE_IN_THRESHOLD_MS=0.001 # workaround for copying to/from unpopulated buffers on GH

# Cray MPICH
export MPICH_GPU_SUPPORT_ENABLED=1
export FI_CXI_RX_MATCH_MODE=software
export FI_MR_CACHE_MONITOR=disabled
export MPICH_NO_BUFFER_ALIAS_CHECK=1

# CUSTOM env vars to MFC
export NVIDIA_ALLOC_MODE=0                    # do nothing
export NVIDIA_MANUAL_GPU_HINTS=1              # prefloc GPU on some
export NVIDIA_IGR_TEMPS_ON_GPU=3              # jac, jac_rhs, and jac_old on GPU
export NVIDIA_VARS_ON_GPU=7                   # q_cons_ts(1)%vf%sf for j=1-7 on GPU

# NSYS
export NSYS=0                                 # enable nsys profiling
export NSYS_FILE=myreport.qdrep

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOT_DIR}"
% if engine == 'batch':
. ./mfc.sh load -c san -m ${'g' if gpu else 'c'}
% endif
cd - > /dev/null
echo

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        (set -x; srun --unbuffered \
                --ntasks=${nodes*tasks_per_node}                     \
                --cpus-per-task 1                                    \
                --cpu-bind=none                                      \
            % if gpu:
                --gpus-per-task 1                                    \
            % endif
                --wait 200 --bcast=/tmp/${target.name}               \
                "${target.get_home_dirpath(case)}/misc/nvidia_uvm/bind.sh" \
            #% if target.name == 'simulation':
                #"${target.get_home_dirpath(case)}/misc/nvidia_uvm/nsys.sh" \
            #% endif
                "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
