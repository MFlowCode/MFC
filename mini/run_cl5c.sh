#!/bin/bash
#SBATCH -p pvc
#SBATCH -N 1
#SBATCH --gres=gpu:pvc:1
#SBATCH -t 0:20:00
#SBATCH -o /scratch/user/u.sb27915/MFC-intel/mini/cl5c.log
#SBATCH -e /scratch/user/u.sb27915/MFC-intel/mini/cl5c.log
#SBATCH -J cl5c-ice-matmul

source /etc/profile
module load iimpi/2025a imkl/2025.1.0
export I_MPI_F90=ifx

F90=/sw/eb/sw/impi/2021.15.0-intel-compilers-2025.1.1/mpi/2021.15/bin/mpif90
MINI=/scratch/user/u.sb27915/MFC-intel/mini
BASE="-free -fiopenmp -fopenmp-targets=spir64"
SRC="$MINI/m_cl5_mini.f90"

echo "=== matmul vs rot(3,3) field isolation test ==="
$F90 --version 2>&1 | head -1
echo ""

try_f90() {
    local desc="$1"; local src="$2"; shift 2
    local out
    rm -f /tmp/cl5c_$$.x m_cl5.mod
    out=$($F90 $BASE "$@" "$src" -o /tmp/cl5c_$$.x 2>&1)
    local rc=$?
    rm -f /tmp/cl5c_$$.x m_cl5.mod
    if [ $rc -eq 0 ]; then echo "OK   : $desc"
    elif echo "$out" | grep -qiE "5633|internal abort|segmentation"; then
        echo "ICE  : $desc"
    else
        echo "ERR  : $desc"
        echo "     $(echo "$out" | grep -v "^$\|warning\|remark\|#5117" | head -2)"
    fi
}

# A: Keep rot(3,3) field but replace matmul with element-wise multiply
python3 - "$SRC" > /tmp/no_matmul_$$.f90 << 'PYEOF'
import sys, re
src = open(sys.argv[1]).read()
# Replace matmul(X, Y) with sum of manually computed dot products (fake matmul)
# Just replace matmul calls with equivalent array expressions
src = re.sub(r'matmul\(patches\(id\)%rot_inv, dv\)', 'dv', src)
src = re.sub(r'matmul\(patches\(id\)%rot_inv, xyz_local\)', 'dv', src)
src = re.sub(r'matmul\(patches\(id\)%rot, dv\)', 'dv', src)
src = re.sub(r'matmul\(patches\(id\)%rot_inv, ldv\)', 'ldv', src)
open('/dev/stdout', 'w').write(src)
PYEOF
try_f90 "rot(3,3) field kept, matmul REMOVED (identity)" /tmp/no_matmul_$$.f90 -O3
rm -f /tmp/no_matmul_$$.f90

# B: Remove rot(3,3) fields but keep matmul via a local array
python3 - "$SRC" > /tmp/no_rotfield_$$.f90 << 'PYEOF'
import sys, re
src = open(sys.argv[1]).read()
# Remove rot/rot_inv struct fields
src = re.sub(r'\s*real\(wp\), dimension\(1:3,1:3\) :: rot, rot_inv\n', '\n', src)
# Replace struct accesses with local identity
src = re.sub(r'patches\(id\)%rot_inv', 'reshape([1._wp,0._wp,0._wp,0._wp,1._wp,0._wp,0._wp,0._wp,1._wp],[3,3])', src)
src = re.sub(r'patches\(id\)%rot\b', 'reshape([1._wp,0._wp,0._wp,0._wp,1._wp,0._wp,0._wp,0._wp,1._wp],[3,3])', src)
# Remove initialization of rot/rot_inv in s_init
src = re.sub(r'\s*patches\(i\)%rot.*\n', '\n', src)
open('/dev/stdout', 'w').write(src)
PYEOF
try_f90 "rot(3,3) field REMOVED, matmul kept (local identity)" /tmp/no_rotfield_$$.f90 -O3
rm -f /tmp/no_rotfield_$$.f90

# C: Both removed — equivalent to previous "no rotation" test
python3 - "$SRC" > /tmp/neither_$$.f90 << 'PYEOF'
import sys, re
src = open(sys.argv[1]).read()
src = re.sub(r'\s*real\(wp\), dimension\(1:3,1:3\) :: rot, rot_inv\n', '\n', src)
src = re.sub(r'matmul\(.*?\)', 'dv', src)
src = re.sub(r'\s*patches\(i\)%rot.*\n', '\n', src)
open('/dev/stdout', 'w').write(src)
PYEOF
try_f90 "both rot(3,3) AND matmul removed" /tmp/neither_$$.f90 -O3
rm -f /tmp/neither_$$.f90

# D: Original with just geo3 (which has matmul) — minimum loop count?
python3 - "$SRC" > /tmp/geo3only_$$.f90 << 'PYEOF'
import sys, re
src = open(sys.argv[1]).read()
# Keep only 1 loop with matmul (geo3)
for g in [1,2,4,5,6,7,8,9,10]:
    src = re.sub(
        r'\s*!\$omp target teams loop[^\n]*\n\s*do i = 1, n\n\s*if \(patches.*geometry == %d\) call s_geo%d.*\n\s*end do\n' % (g,g),
        '\n', src)
    src = re.sub(
        r'\n\s*subroutine s_geo%d\(gp\).*?end subroutine s_geo%d\n' % (g, g),
        '\n', src, flags=re.DOTALL)
open('/dev/stdout', 'w').write(src)
PYEOF
try_f90 "1 loop (geo3 only, has matmul)" /tmp/geo3only_$$.f90 -O3
rm -f /tmp/geo3only_$$.f90

# E: 2 loops with matmul (geo3 + geo4)
python3 - "$SRC" > /tmp/geo34_$$.f90 << 'PYEOF'
import sys, re
src = open(sys.argv[1]).read()
for g in [1,2,5,6,7,8,9,10]:
    src = re.sub(
        r'\s*!\$omp target teams loop[^\n]*\n\s*do i = 1, n\n\s*if \(patches.*geometry == %d\) call s_geo%d.*\n\s*end do\n' % (g,g),
        '\n', src)
    src = re.sub(
        r'\n\s*subroutine s_geo%d\(gp\).*?end subroutine s_geo%d\n' % (g, g),
        '\n', src, flags=re.DOTALL)
open('/dev/stdout', 'w').write(src)
PYEOF
try_f90 "2 loops with matmul (geo3 + geo4)" /tmp/geo34_$$.f90 -O3
rm -f /tmp/geo34_$$.f90

echo ""
echo "=== Done ==="
