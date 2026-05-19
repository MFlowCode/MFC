#!/bin/bash
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --gres=gpu:pvc:1
#SBATCH -t 0:20:00
#SBATCH -o /scratch/user/u.sb27915/MFC-intel/mini/cl4.log
#SBATCH -e /scratch/user/u.sb27915/MFC-intel/mini/cl4.log
#SBATCH -J cl4-ice-repro

source /etc/profile
module load iimpi/2025a imkl/2025.1.0
export I_MPI_F90=ifx

F90=/sw/eb/sw/impi/2021.15.0-intel-compilers-2025.1.1/mpi/2021.15/bin/mpif90
MINI=/scratch/user/u.sb27915/MFC-intel/mini
BASE="-free -fiopenmp -fopenmp-targets=spir64 -fp-model=precise -march=native -mno-avx512fp16"

echo "=== m_cl4_mini (single-module same-module declare-target + char struct) ==="
$F90 --version 2>&1 | head -1
echo ""

try() {
    local desc="$1"; shift
    local out
    rm -f /tmp/cl4_$$.o /tmp/cl4_$$.x m_cl4_mini.mod
    out=$($F90 $BASE "$@" "$MINI/m_cl4_mini.f90" -o /tmp/cl4_$$.x 2>&1)
    local rc=$?
    rm -f /tmp/cl4_$$.o /tmp/cl4_$$.x m_cl4_mini.mod
    if [ $rc -eq 0 ]; then
        echo "OK   : $desc"
    elif echo "$out" | grep -qiE "5623|5633|internal abort|segmentation"; then
        echo "ICE  : $desc"
        echo "     $(echo "$out" | grep -iE "5623|5633|internal abort|error #" | head -1)"
    else
        echo "ERR  : $desc"
        echo "     $(echo "$out" | grep -v "^$\|warning\|remark\|#5117\|preprocessor" | head -3)"
    fi
}

try "-O3 (baseline)"           -O3
try "-O2"                      -O2
try "-O1"                      -O1
try "-O3 -fno-inline"          -O3 -fno-inline
try "-O3 -fno-vectorize"       -O3 -fno-vectorize
try "-O3 -fno-unroll-loops"    -O3 -fno-unroll-loops

echo ""
echo "=== Done ==="
