#!/bin/bash
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --gres=gpu:pvc:1
#SBATCH -t 0:15:00
#SBATCH -o /scratch/user/u.sb27915/MFC-intel/mini/cl5.log
#SBATCH -e /scratch/user/u.sb27915/MFC-intel/mini/cl5.log
#SBATCH -J cl5-ice-reduce

source /etc/profile
module load iimpi/2025a imkl/2025.1.0
export I_MPI_F90=ifx

F90=/sw/eb/sw/impi/2021.15.0-intel-compilers-2025.1.1/mpi/2021.15/bin/mpif90
MINI=/scratch/user/u.sb27915/MFC-intel/mini
BASE="-free -fiopenmp -fopenmp-targets=spir64"

echo "=== m_cl5_mini (reduced: char+array struct + 10 target loops + 10 same-module decl-target subs) ==="
$F90 --version 2>&1 | head -1
echo ""

try() {
    local desc="$1"; shift
    local out
    rm -f /tmp/cl5_$$.x m_cl5.mod
    out=$($F90 $BASE "$@" "$MINI/m_cl5_mini.f90" -o /tmp/cl5_$$.x 2>&1)
    local rc=$?
    rm -f /tmp/cl5_$$.x m_cl5.mod
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
try "-O3 -fno-inline"          -O3 -fno-inline
try "-O2"                      -O2
try "-O1"                      -O1

echo ""
echo "=== Done ==="
