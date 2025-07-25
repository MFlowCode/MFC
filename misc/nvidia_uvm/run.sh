#!/usr/bin/env bash

# TODO: Modify accordingly
PATH_TO_BINARY=${SCRATCH}/projects/cfd/mfc/MFC-Wilfong/build/install/cdcd4e8762/bin/

# NVHPC and CUDA env vars
export NV_ACC_USE_MALLOC=1                    # use malloc instead of cudaMallocManaged ( compiled using -gpu=mem:unified )
export NVCOMPILER_ACC_NO_MEMHINTS=1           # disable implicit compiler hints
export CUDA_BUFFER_PAGE_IN_THRESHOLD_MS=0.001 # workaround for copying to/from unpopulated buffers on GH

# Cray MPICH
export MPICH_GPU_SUPPORT_ENABLED=1            # MPICH with GPU support
export FI_CXI_RX_MATCH_MODE=software
export FI_MR_CACHE_MONITOR=disabled

# CUSTOM env vars to MFC
export NVIDIA_ALLOC_MODE=2                    # default alloc to prefloc CPU
export NVIDIA_MANUAL_GPU_HINTS=1              # prefloc GPU on some
export NVIDIA_IGR_TEMPS_ON_GPU=1              # jac on GPU and jac_rhs on CPU       ( NOTE: good default, tune based on size )
export NVIDIA_VARS_ON_GPU=7                   # q_cons_ts(1)%vf%sf for j=1-7 on GPU ( NOTE: good default, tune based on size )

# NSYS
export NSYS=1                                 # enable nsys profiling
export NSYS_FILE=report_uvm_single_N-499_nGPUs-4_params-${NVIDIA_VARS_ON_GPU}-${NVIDIA_IGR_TEMPS_ON_GPU}.qdrep

# Run using --cpu-bind=none because we use our own binding script
srun --ntasks 4 --cpu-bind=none ./bind.sh ./nsys.sh ${PATH_TO_BINARY}/simulation
