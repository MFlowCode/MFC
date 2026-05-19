#!/bin/bash
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --gres=gpu:pvc:1
#SBATCH -t 1:30:00
#SBATCH --mem=32G
#SBATCH -o /scratch/user/u.sb27915/MFC-intel/build_intel_gpu.log
#SBATCH -e /scratch/user/u.sb27915/MFC-intel/build_intel_gpu.log
#SBATCH -J mfc-intel-gpu-build

source /etc/profile
module load iimpi/2025a imkl/2025.1.0 CMake/3.31.3 Python/3.13.1
export I_MPI_F90=ifx FC=mpif90
export UV_CACHE_DIR=/scratch/user/u.sb27915/.cache/uv
export RUSTUP_HOME=/scratch/user/u.sb27915/.rustup
export CARGO_HOME=/scratch/user/u.sb27915/.cargo
export PATH=$CARGO_HOME/bin:$PATH

cd /scratch/user/u.sb27915/MFC-intel
./mfc.sh build -t simulation --gpu mp --intel-aot -j 8
