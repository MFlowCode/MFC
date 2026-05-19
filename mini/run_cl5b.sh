#!/bin/bash
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --gres=gpu:pvc:1
#SBATCH -t 0:20:00
#SBATCH -o /scratch/user/u.sb27915/MFC-intel/mini/cl5b.log
#SBATCH -e /scratch/user/u.sb27915/MFC-intel/mini/cl5b.log
#SBATCH -J cl5b-ice-bisect

source /etc/profile
module load iimpi/2025a imkl/2025.1.0
export I_MPI_F90=ifx

F90=/sw/eb/sw/impi/2021.15.0-intel-compilers-2025.1.1/mpi/2021.15/bin/mpif90
MINI=/scratch/user/u.sb27915/MFC-intel/mini
BASE="-free -fiopenmp -fopenmp-targets=spir64"
SRC="$MINI/m_cl5_mini.f90"

echo "=== m_cl5 bisection: what is the minimum ICE trigger? ==="
$F90 --version 2>&1 | head -1
echo ""

try_src() {
    local desc="$1"; local src="$2"; shift 2
    local out
    rm -f /tmp/cl5b_$$.x m_cl5.mod
    out=$($F90 $BASE "$@" "$src" -o /tmp/cl5b_$$.x 2>&1)
    local rc=$?
    rm -f /tmp/cl5b_$$.x m_cl5.mod
    if [ $rc -eq 0 ]; then
        echo "OK   : $desc"
    elif echo "$out" | grep -qiE "5623|5633|internal abort|segmentation"; then
        echo "ICE  : $desc"
    else
        echo "ERR  : $desc"
        echo "     $(echo "$out" | grep -v "^$\|warning\|remark\|#5117" | head -2)"
    fi
}

# Baseline confirms ICE
try_src "baseline (full cl5)" "$SRC" -O3

# --- Remove char field ---
sed '/character(LEN=200)/d; s/patches(i)%label.*=.*//g' "$SRC" > /tmp/nochar_$$.f90
try_src "no char(LEN=200) field" /tmp/nochar_$$.f90 -O3
rm -f /tmp/nochar_$$.f90

# --- Remove interp_coeffs(2,2,2) ---
sed '/interp_coeffs/d' "$SRC" > /tmp/no_ic_$$.f90
try_src "no interp_coeffs(2,2,2)" /tmp/no_ic_$$.f90 -O3
rm -f /tmp/no_ic_$$.f90

# --- Keep only 5 loops (geo1-5) ---
python3 - "$SRC" > /tmp/cl5_5loops_$$.f90 << 'PYEOF'
import sys, re
src = open(sys.argv[1]).read()
# Remove geo6-geo10 loops from s_dispatch
for g in range(6, 11):
    src = re.sub(
        r'\s*!\$omp target teams loop[^\n]*\n\s*do i = 1, n\n\s*if \(patches.*geometry == %d\) call s_geo%d.*\n\s*end do\n' % (g,g),
        '\n', src)
# Remove geo6-geo10 subroutines
for g in range(6, 11):
    src = re.sub(
        r'\n\s*subroutine s_geo%d\(gp\).*?end subroutine s_geo%d\n' % (g, g),
        '\n', src, flags=re.DOTALL)
open('/dev/stdout', 'w').write(src)
PYEOF
try_src "5 loops only (geo1-5)" /tmp/cl5_5loops_$$.f90 -O3
rm -f /tmp/cl5_5loops_$$.f90

# --- Remove rot/rot_inv (3x3 arrays) from patch_t ---
sed '/real(wp), dimension(1:3,1:3)/d; /patches(i)%rot/d; /patches(id)%rot/d; /matmul/d' "$SRC" > /tmp/norot_$$.f90
try_src "no rotation matrices" /tmp/norot_$$.f90 -O3
rm -f /tmp/norot_$$.f90

# --- Remove allocatable xc/yc/zc (use constants instead) ---
python3 - "$SRC" > /tmp/noxc_$$.f90 << 'PYEOF'
import sys, re
src = open(sys.argv[1]).read()
# Replace allocatable declarations with fixed arrays
src = src.replace('real(wp), allocatable :: xc(:), yc(:), zc(:)', 'real(wp) :: xc(64), yc(64), zc(64)')
# Remove allocate/deallocate
src = re.sub(r'\s*allocate\(xc.*?\)\n', '\n', src)
src = re.sub(r'\s*deallocate\(xc.*?\)\n', '\n', src)
src = src.replace('integer, intent(in) :: n\n        integer :: i\n        allocate', 'integer, intent(in) :: n\n        integer :: i\n        ! allocate')
open('/dev/stdout', 'w').write(src)
PYEOF
try_src "no allocatable xc/yc/zc" /tmp/noxc_$$.f90 -O3
rm -f /tmp/noxc_$$.f90

echo ""
echo "=== Done ==="
