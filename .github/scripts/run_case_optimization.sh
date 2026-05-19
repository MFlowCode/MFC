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

benchmarks=(
    benchmarks/5eq_rk3_weno3_hllc/case.py
    benchmarks/viscous_weno5_sgb_acoustic/case.py
    benchmarks/hypo_hll/case.py
    benchmarks/ibm/case.py
    benchmarks/igr/case.py
)

# For Frontier/Frontier AMD: deps were fetched on the login node via --deps-only;
# build case-optimized binaries here on the compute node before running.
# For Phoenix and frontier_amd: prebuild-case-optimization.sh already built
# everything in a prior SLURM job (via --dry-run), so skip the build here.
#
# Clean stale MFC target staging before building. On self-hosted CI runners,
# corrupted intermediate files from a prior failed build (e.g. CCE optcg crash)
# can persist and poison subsequent builds. Each case-opt config gets its own
# hash-named staging dir, but install dirs and other artifacts may be stale.
if [ "$job_cluster" != "phoenix" ] && [ "$job_cluster" != "frontier_amd" ]; then
    # Clean stale MFC target dirs (hash-named) from prior builds, but
    # preserve dependency dirs (hipfort, fftw, etc.) since the compute
    # node has no internet to re-fetch them.
    echo "=== Cleaning stale MFC target staging/install ==="
    find build/staging -maxdepth 1 -regex '.*/[0-9a-f]+' -type d -exec rm -rf {} + 2>/dev/null || true
    find build/install -maxdepth 1 -regex '.*/[0-9a-f]+' -type d -exec rm -rf {} + 2>/dev/null || true

    echo "=== Building case-optimized binaries on compute node ==="
    for case in "${benchmarks[@]}"; do
        echo "--- Building: $case ---"
        ./mfc.sh build -i "$case" --case-optimization $gpu_opts -j 8
    done
    echo "=== All case-optimized binaries built ==="
fi

passed=0
failed=0
failed_cases=""

for case in "${benchmarks[@]}"; do
    case_dir="$(dirname "$case")"
    case_name="$(basename "$case_dir")"
    echo ""
    echo "===================="
    echo "Case-optimization test: $case_name"
    echo "===================="

    # Clean any previous output
    rm -rf "$case_dir/D" "$case_dir/p_all" "$case_dir/restart_data"

    # Build + run with --case-optimization, small grid, 10 timesteps
    if ./mfc.sh run "$case" --case-optimization $gpu_opts -n "$ngpus" -j 8 -c "$job_cluster" -- --gbpp 1 --steps 10; then
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
echo "===================="
echo "Case-optimization summary: $passed passed, $failed failed"
if [ $failed -gt 0 ]; then
    echo "Failed cases:$failed_cases"
fi
echo "===================="

[ $failed -eq 0 ] && exit 0 || exit 1
