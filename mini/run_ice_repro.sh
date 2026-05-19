#!/bin/bash
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --gres=gpu:pvc:1
#SBATCH -t 0:15:00
#SBATCH -o /scratch/user/u.sb27915/MFC-intel/mini/ice_repro.log
#SBATCH -e /scratch/user/u.sb27915/MFC-intel/mini/ice_repro.log
#SBATCH -J ice-repro-final

source /etc/profile
module load iimpi/2025a imkl/2025.1.0
export I_MPI_F90=ifx

# Use ifx directly (not mpif90) since this is a single-file reproducer
IFX=/sw/eb/sw/intel-compilers/2025.1.1/compiler/2025.1/bin/ifx
MINI=/scratch/user/u.sb27915/MFC-intel/mini
BASE="-free -fiopenmp -fopenmp-targets=spir64"

echo "=== ifx #5633 minimal reproducer: matmul(derived_type%mat, vec) in declare-target sub ==="
$IFX --version 2>&1 | head -1
echo ""

try() {
    local desc="$1"; shift
    local out
    rm -f /tmp/ice_$$.x m_ice.mod
    out=$($IFX $BASE "$@" "$MINI/m_ice_repro.f90" -o /tmp/ice_$$.x 2>&1)
    local rc=$?
    rm -f /tmp/ice_$$.x m_ice.mod
    if [ $rc -eq 0 ]; then
        echo "OK   : $desc"
    elif echo "$out" | grep -qiE "5633|internal abort|segmentation"; then
        echo "ICE  : $desc"
        echo "     $(echo "$out" | grep -iE "5633|error #" | head -1)"
    else
        echo "ERR  : $desc"
        echo "     $(echo "$out" | grep -v "^$\|warning\|remark\|#5117" | head -2)"
    fi
}

try "-O3"                      -O3
try "-O2"                      -O2
try "-O1"                      -O1
try "-O3 -fno-inline (fix)"    -O3 -fno-inline
try "-O0"                      -O0

echo ""
echo "--- Variant: matmul on local array (no struct member) ---"
python3 - "$MINI/m_ice_repro.f90" > /tmp/local_mat_$$.f90 << 'PYEOF'
import sys
src = open(sys.argv[1]).read()
# Add local_mat to the declaration block (before first executable)
src = src.replace(
    '        real(wp), dimension(3) :: dv, ldv\n        real(wp) :: d',
    '        real(wp), dimension(3) :: dv, ldv\n        real(wp), dimension(3,3) :: local_mat\n        real(wp) :: d')
# Copy struct member to local before matmul
src = src.replace(
    '        ! matmul with a derived-type member (patches(id)%mat) is the key trigger\n        ldv = matmul(patches(id)%mat, dv)',
    '        local_mat = patches(id)%mat\n        ! matmul with LOCAL copy (not struct member directly)\n        ldv = matmul(local_mat, dv)')
open('/dev/stdout', 'w').write(src)
PYEOF
rm -f /tmp/ice_lm_$$.x m_ice.mod
out=$($IFX $BASE -O3 /tmp/local_mat_$$.f90 -o /tmp/ice_lm_$$.x 2>&1)
rc=$?
rm -f /tmp/ice_lm_$$.x m_ice.mod /tmp/local_mat_$$.f90
if [ $rc -eq 0 ]; then
    echo "OK   : -O3 with local_mat copy -- direct struct access IS the trigger"
elif echo "$out" | grep -qiE "5633|internal abort|segmentation"; then
    echo "ICE  : -O3 with local_mat copy -- copying to local doesn't help"
else
    echo "ERR  : $(echo "$out" | grep -v '^$' | head -1)"
fi

echo ""
echo "=== Done ==="
