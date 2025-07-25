#!/usr/bin/env bash

# -------------------------------- #
# Binding for a single Santis node #
# -------------------------------- #

# Local rank
export local_rank="${OMPI_COMM_WORLD_LOCAL_RANK:-$SLURM_LOCALID}"

# Bind to GPU
export CUDA_VISIBLE_DEVICES="$local_rank"

# Binding to NIC
export MPICH_OFI_NIC_POLICY=USER
export MPICH_OFI_NIC_MAPPING="0:0; 1:1; 2:2; 3:3"

# Bind to cores ( second core per socket )
physcores=(0 72 144 216)

#echo rank: $local_rank, cores: ${physcores[$local_rank]}, GPU: $CUDA_VISIBLE_DEVICES, NIC mapping: $MPICH_OFI_NIC_POLICY
numactl -l --all --physcpubind=${physcores[$local_rank]} "$@"
