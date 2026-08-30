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

# Benchmark list (single source of truth shared with the pre-build) + optional
# sharding: a sharded run executes only the subset it has pre-built binaries
# for. $job_shard mirrors the pre-build shard so the two stay aligned.
source .github/scripts/case-optimization-benchmarks.sh
caseopt_parse_shard

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
    find build/staging -maxdepth 1 -regex '.*/\(gpu-acc\|gpu-mp\|cpu\)-.*' -type d -exec rm -rf {} + 2>/dev/null || true
    find build/install -maxdepth 1 -regex '.*/\(gpu-acc\|gpu-mp\|cpu\)-.*' -type d -exec rm -rf {} + 2>/dev/null || true

    echo "=== Building case-optimized binaries on compute node ==="
    idx=0
    for case in "${benchmarks[@]}"; do
        idx=$((idx + 1))
        caseopt_case_in_shard "$idx" || continue
        echo "--- Building: $case ---"
        ./mfc.sh build -i "$case" --case-optimization $gpu_opts -j 8
    done
    echo "=== All case-optimized binaries built ==="
    build_opts=""
else
    # Nothing was built here, so `mfc.sh run` must not build either. Left to
    # itself it re-runs `cmake --build` and `cmake --install` for every target
    # on every case, and syscheck/pre_process/post_process hash to a single
    # slug shared by all cases, so every shard installs to the same
    # build/staging and build/install paths. The shards share this workspace
    # and leave the SLURM queue together, then march through an identical
    # startup, so those installs collide and the losing shard dies with
    # "file INSTALL cannot copy file ... No such file or directory".
    # The pre-build serializes its own shared-target build for this reason;
    # here there is nothing to serialize, because there is nothing to build.
    build_opts="--no-build"
fi

passed=0
failed=0
failed_cases=""

idx=0
for case in "${benchmarks[@]}"; do
    idx=$((idx + 1))
    caseopt_case_in_shard "$idx" || continue
    case_dir="$(dirname "$case")"
    case_name="$(basename "$case_dir")"
    echo ""
    echo "===================="
    echo "Case-optimization test: $case_name"
    echo "===================="

    # Clean any previous output
    rm -rf "$case_dir/D" "$case_dir/p_all" "$case_dir/restart_data"

    # Run with --case-optimization, small grid, 10 timesteps. On phoenix and
    # frontier_amd $build_opts is --no-build, so the run reuses the pre-build's
    # binaries. On phoenix (a single, unsharded job) a prebuilt binary can go
    # missing at run time (the pre-build and run SLURM jobs may be separated by
    # a long embers queue wait, or clean_build on a retry may have wiped
    # build/), which makes mpirun fail to launch it. Only that specific failure
    # triggers a rebuild-and-rerun, so a lost artifact degrades to a slower run
    # instead of a red CI; any other failure (a real crash, a NaN, an MPI fault)
    # is reported as-is and never masked by the retry. frontier_amd is excluded:
    # its run is sharded across concurrent jobs sharing one workspace, so a
    # fallback rebuild would race on the shared install paths (the collision the
    # --no-build guard above prevents).
    run_log="$(mktemp)"
    ./mfc.sh run "$case" --case-optimization $gpu_opts $build_opts -n "$ngpus" -j 8 -c "$job_cluster" -- --gbpp 1 --steps 10 2>&1 | tee "$run_log"
    run_rc=${PIPESTATUS[0]}
    if [ "$run_rc" -eq 0 ]; then
        run_ok=1
    elif [ "$job_cluster" = phoenix ] && grep -q "could not access or execute" "$run_log"; then
        echo "NOTE: $case_name rebuilding in-run (prebuilt binary was unavailable)"
        rm -rf "$case_dir/D" "$case_dir/p_all" "$case_dir/restart_data"
        if ./mfc.sh run "$case" --case-optimization $gpu_opts -n "$ngpus" -j 8 -c "$job_cluster" -- --gbpp 1 --steps 10; then
            run_ok=1
        else
            run_ok=0
        fi
    else
        run_ok=0
    fi
    rm -f "$run_log"

    if [ "$run_ok" = 1 ]; then
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
