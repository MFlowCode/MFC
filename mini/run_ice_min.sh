#!/bin/bash
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --gres=gpu:pvc:1
#SBATCH -t 0:15:00
#SBATCH -o /scratch/user/u.sb27915/MFC-intel/mini/ice_min.log
#SBATCH -e /scratch/user/u.sb27915/MFC-intel/mini/ice_min.log
#SBATCH -J ice-min

source /etc/profile
module load iimpi/2025a imkl/2025.1.0

IFX=/sw/eb/sw/intel-compilers/2025.1.1/compiler/2025.1/bin/ifx
MINI=/scratch/user/u.sb27915/MFC-intel/mini
BASE="-free -fiopenmp -fopenmp-targets=spir64"

echo "=== ifx #5633 absolute minimum: matmul in declare-target sub called from target teams loop ==="
$IFX --version 2>&1 | head -1
echo ""

try() {
    local desc="$1"; shift
    local out
    rm -f /tmp/icemin_$$.x m_min.mod
    out=$($IFX $BASE "$@" "$MINI/m_ice_min.f90" -o /tmp/icemin_$$.x 2>&1)
    local rc=$?
    rm -f /tmp/icemin_$$.x m_min.mod
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

try "-O3"                   -O3
try "-O2"                   -O2
try "-O1"                   -O1
try "-O0"                   -O0
try "-O3 -fno-inline"       -O3 -fno-inline

echo ""
echo "--- Variant: replace matmul with element-wise loop (no matmul intrinsic) ---"
python3 - "$MINI/m_ice_min.f90" > /tmp/no_matmul_$$.f90 << 'PYEOF'
import sys
src = open(sys.argv[1]).read()
src = src.replace(
    "        b = real(k, wp) * matmul(A, x)",
    "        integer :: ii, jj\n        do ii = 1, 3\n            b(ii) = 0._wp\n            do jj = 1, 3\n                b(ii) = b(ii) + A(ii,jj)*x(jj)\n            end do\n            b(ii) = b(ii) * real(k, wp)\n        end do")
open('/dev/stdout', 'w').write(src)
PYEOF
rm -f /tmp/icemin_nm_$$.x m_min.mod
out=$($IFX $BASE -O3 /tmp/no_matmul_$$.f90 -o /tmp/icemin_nm_$$.x 2>&1)
rc=$?
rm -f /tmp/icemin_nm_$$.x m_min.mod /tmp/no_matmul_$$.f90
if [ $rc -eq 0 ]; then
    echo "OK   : -O3 no matmul (manual loops) -- matmul intrinsic IS the trigger"
elif echo "$out" | grep -qiE "5633|segmentation"; then
    echo "ICE  : -O3 no matmul -- ICE even without matmul intrinsic"
else
    echo "ERR  : $(echo "$out" | grep -v '^$' | head -1)"
fi

echo ""
echo "=== Done ==="
