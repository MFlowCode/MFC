#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>

% if engine == 'batch':
#SBATCH --uenv=icon/25.2:v1@santis
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --cpus-per-task=72
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

# We compiled the code using -gpu=unified:managedalloc, hence we use cudaMallocManaged for the dynamic allocations.
# Using NV_ACC_USE_MALLOC we could change to malloc at runtime. We choose to not do that here and stick with cudaMallocManaged and 2MB page sizes.
# https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html#memory-model
# https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html#command-line-options-selecting-compiler-memory-modes
export NV_ACC_USE_MALLOC=0

# For NVIDIA CUDA devices, controls the use of automatic memory hints at data constructs in the managed and unified memory modes.
# Below is a breakdown of the permitted values (case insensitive):
# - DEFAULT: Use the default settings. On NVIDIA Grace Hopper systems, the default is currently ENABLE_ALL; on all other systems, the default is DISABLE.
# - DISABLE: Memory hints are disabled for all data constructs.
# - ENABLE_EXPLICIT: Memory hints are enabled for explicit data constructs only.
# - ENABLE_ALL: Memory hints are enabled for explicit and implicit data constructs.
# https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html#environment-variables-controlling-device-memory-management
# Here we disable the implicit compiler hints.
# Using NVCOMPILER_ACC_NO_MEMHINTS is the legacy way and is still supported, but users should prefer NVCOMPILER_ACC_MEMHINTS when using newer nvhpc compilers.
export NVCOMPILER_ACC_NO_MEMHINTS=1           # disable implicit compiler hints - legacy way
export NVCOMPILER_ACC_MEMHINTS=DISABLE        # disable implicit compiler hints - new way

# Cray MPICH
export MPICH_GPU_SUPPORT_ENABLED=1
export FI_CXI_RX_MATCH_MODE=software
export FI_MR_CACHE_MONITOR=disabled
export MPICH_NO_BUFFER_ALIAS_CHECK=1

# NSYS
export NSYS=0                                 # enable nsys profiling
export NSYS_FILE=myreport.qdrep

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOT_DIR}"
% if engine == 'batch':
. ./mfc.sh load -c san -m ${'g' if gpu_enabled else 'c'}
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
                --cpus-per-task 72                                   \
                --cpu-bind=none                                      \
            % if gpu_enabled:
                --gpus-per-task 1                                    \
            % endif
                --wait 200 --bcast=/tmp/${target.name}               \
                "${target.get_home_dirpath()}/misc/nvidia_uvm/bind.sh" \
            % if target.name == 'simulation':
                "${target.get_home_dirpath()}/misc/nvidia_uvm/nsys.sh" \
            % endif
                "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
