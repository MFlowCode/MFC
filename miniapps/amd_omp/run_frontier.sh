#!/bin/bash
#SBATCH -A cfd154
#SBATCH -J amd_omp_miniapps
#SBATCH -o %x_%j.out
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -q hackathon

cd /lustre/orion/cfd154/scratch/sbryngelson/MFC-derivedtypes/miniapps/amd_omp

module load python cmake cpe/25.09 PrgEnv-amd amd/7.2.0 rocm/7.2.0

make clean
make -j6 2>&1 | grep -E "Error|error|warning|Warning|FAIL|PASS|^make" || true
echo ""
echo "=== Unfiltered build for tests 07-10 ==="
for t in test07_module_param_array test08_dt_private_array_elems test09_dt_arr7_arithmetic test10_hlld_flux; do
    echo "--- Building $t ---"
    ftn -fopenmp --offload-arch=gfx90a -fopenmp-target-fast \
        -fopenmp-assume-threads-oversubscription \
        -fopenmp-assume-teams-oversubscription \
        -o $t $t.f90 2>&1
    echo "exit: $?"
done

echo ""
echo "=== Running miniapps ==="
for t in test01_dt_elem_access test02_dt_array_constructor test03_dt_whole_array_ops \
          test04_dt_zero_init test05_device_routine_int_arg test06_dt_member_as_array_arg \
          test07_module_param_array \
          test08_dt_private_array_elems test09_dt_arr7_arithmetic test10_hlld_flux; do
    if [ -x "./$t" ]; then
        echo "--- $t ---"
        ./$t
    else
        echo "--- $t --- SKIPPED (did not build)"
    fi
done
echo "=== Done ==="
