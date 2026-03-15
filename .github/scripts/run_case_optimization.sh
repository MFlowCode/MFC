#!/bin/bash

# Case-optimization CI test script.
# Runs inside a SLURM job — expects $job_device and $job_interface from submit.sh.

set -e

source .github/scripts/detect-gpus.sh
source .github/scripts/gpu-opts.sh

# Default to 1 GPU if detection found none but we're in GPU mode
if [ "$job_device" = "gpu" ] && [ "$ngpus" -eq 0 ]; then
    ngpus=1
fi

# Verify the venv Python interpreter exists (created by ./mfc.sh build)
if [ ! -x build/venv/bin/python3 ]; then
    echo "ERROR: build/venv/bin/python3 not found."
    echo "The MFC build venv may not have been created. Was the pre-build step successful?"
    exit 1
fi

benchmarks=(
    benchmarks/5eq_rk3_weno3_hllc/case.py
    benchmarks/viscous_weno5_sgb_acoustic/case.py
    benchmarks/hypo_hll/case.py
    benchmarks/ibm/case.py
    benchmarks/igr/case.py
)

passed=0
failed=0
failed_cases=""

for case in "${benchmarks[@]}"; do
    case_dir="$(dirname "$case")"
    case_name="$(basename "$case_dir")"
    echo ""
    echo "========================================"
    echo "Case-optimization test: $case_name"
    echo "========================================"

    # Clean any previous output
    rm -rf "$case_dir/D" "$case_dir/p_all" "$case_dir/restart_data"

    # Build + run with --case-optimization, small grid, 10 timesteps
    if ./mfc.sh run "$case" --case-optimization $gpu_opts -n "$ngpus" -j 8 -- --gbpp 1 --steps 10; then
        # Validate output
        if build/venv/bin/python3 .github/scripts/check_case_optimization_output.py "$case_dir"; then
            echo "PASS: $case_name"
            passed=$((passed + 1))
        else
            echo "FAIL: $case_name (validation error)"
            failed=$((failed + 1))
            failed_cases="$failed_cases $case_name"
        fi
    else
        echo "FAIL: $case_name (build or run error)"
        failed=$((failed + 1))
        failed_cases="$failed_cases $case_name"
    fi

    # Clean up output between cases
    rm -rf "$case_dir/D" "$case_dir/p_all" "$case_dir/restart_data"
done

echo ""
echo "========================================"
echo "Case-optimization summary: $passed passed, $failed failed"
if [ $failed -gt 0 ]; then
    echo "Failed cases:$failed_cases"
fi
echo "========================================"

[ $failed -eq 0 ] && exit 0 || exit 1
