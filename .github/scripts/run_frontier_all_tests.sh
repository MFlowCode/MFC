#!/bin/bash
# Orchestrate all Frontier test configs in one multi-node SLURM allocation.
# 1. Builds all configs on the login node (in parallel, different modules each)
# 2. Submits a single 5-node SLURM job running tests in parallel via ssh

set -euo pipefail

# Ignore SIGHUP to survive login node session drops
trap '' HUP

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Config table: cluster device interface
configs=(
    "frontier     gpu acc"
    "frontier     gpu omp"
    "frontier     cpu none"
    "frontier_amd gpu omp"
    "frontier_amd cpu none"
)
num_nodes=${#configs[@]}

echo "=========================================="
echo "Frontier consolidated tests: $num_nodes configs"
echo "=========================================="

# --- Phase 1: Create per-config source copies ---
# Build exclude list to prevent copying into self
excludes=""
for cfg in "${configs[@]}"; do
    read -r cluster device interface <<< "$cfg"
    excludes+=" --exclude=test-${cluster}-${device}-${interface}"
done

for cfg in "${configs[@]}"; do
    read -r cluster device interface <<< "$cfg"
    dir="test-${cluster}-${device}-${interface}"
    echo "Creating source copy: $dir"
    rsync -a --link-dest="$(pwd)" $excludes ./ "$dir/"
done

# --- Phase 2: Build all configs on login node in parallel ---
echo ""
echo "=========================================="
echo "Starting parallel builds (${num_nodes} configs)..."
echo "=========================================="

build_pids=()
for cfg in "${configs[@]}"; do
    read -r cluster device interface <<< "$cfg"
    dir="test-${cluster}-${device}-${interface}"
    log="build-${cluster}-${device}-${interface}.log"
    echo "  Starting: $cluster $device $interface"
    (
        cd "$dir"
        bash .github/workflows/${cluster}/build.sh "$device" "$interface"
    ) > "$log" 2>&1 &
    build_pids+=($!)
done

# Wait for all builds and report results
build_failed=0
build_exits=()
for i in "${!build_pids[@]}"; do
    read -r cluster device interface <<< "${configs[$i]}"
    if wait "${build_pids[$i]}"; then
        echo "  Build PASSED: $cluster $device $interface"
        build_exits+=(0)
    else
        code=$?
        echo "  Build FAILED: $cluster $device $interface (exit $code)"
        build_exits+=($code)
        build_failed=1
    fi
done

# On failure, print logs for failed builds and abort
if [ "$build_failed" -ne 0 ]; then
    echo ""
    echo "=========================================="
    echo "Build failures detected. Printing failed logs:"
    echo "=========================================="
    for i in "${!configs[@]}"; do
        if [ "${build_exits[$i]}" -ne 0 ]; then
            read -r cluster device interface <<< "${configs[$i]}"
            log="build-${cluster}-${device}-${interface}.log"
            echo "--- $log ---"
            cat "$log"
            echo "--- end $log ---"
        fi
    done
    exit 1
fi

echo ""
echo "=========================================="
echo "All builds complete. Submitting 5-node SLURM job..."
echo "=========================================="

# --- Phase 3: Submit one sbatch job with N nodes ---
output_file="test-frontier-all.out"

submit_output=$(sbatch <<'OUTER'
#!/bin/bash
#SBATCH -J MFC-frontier-all-tests
#SBATCH -A ENG160
#SBATCH -N 5
#SBATCH -t 05:59:00
#SBATCH -otest-frontier-all.out
#SBATCH -p extended

set -x

cd "$SLURM_SUBMIT_DIR"
echo "Running in $(pwd)"
echo "Allocated nodes: $SLURM_NODELIST"

# Get list of individual node hostnames
mapfile -t nodes < <(scontrol show hostnames "$SLURM_NODELIST")
echo "Nodes: ${nodes[*]}"

# Config table (must match the outer script)
configs=(
    "frontier     gpu acc"
    "frontier     gpu omp"
    "frontier     cpu none"
    "frontier_amd gpu omp"
    "frontier_amd cpu none"
)

pids=()

for i in "${!configs[@]}"; do
    read -r cluster device interface <<< "${configs[$i]}"
    node="${nodes[$i]}"
    dir="test-${cluster}-${device}-${interface}"
    outfile="test-${cluster}-${device}-${interface}.out"

    echo "[$node] Starting test: $cluster $device $interface in $dir"

    ssh -q -o StrictHostKeyChecking=no "$node" \
        "cd $SLURM_SUBMIT_DIR/$dir && bash .github/scripts/frontier_test_config.sh $cluster $device $interface" \
        > "$outfile" 2>&1 &
    pids+=($!)
done

echo "All test configs launched, waiting for completion..."

# Wait for all and collect exit codes
overall_exit=0
for i in "${!pids[@]}"; do
    read -r cluster device interface <<< "${configs[$i]}"
    pid=${pids[$i]}
    if wait "$pid"; then
        echo "PASSED: $cluster $device $interface (PID $pid)"
    else
        code=$?
        echo "FAILED: $cluster $device $interface (PID $pid, exit code $code)"
        overall_exit=1
    fi
done

# Print summary
echo ""
echo "=========================================="
echo "Test summary:"
for cfg in "${configs[@]}"; do
    read -r cluster device interface <<< "$cfg"
    outfile="test-${cluster}-${device}-${interface}.out"
    if [ -f "$outfile" ]; then
        echo "  $cluster $device $interface: $(tail -n 1 "$outfile")"
    else
        echo "  $cluster $device $interface: NO OUTPUT FILE"
    fi
done
echo "=========================================="

exit $overall_exit
OUTER
)

job_id=$(echo "$submit_output" | awk '/Submitted batch job/ {print $4}')
if [ -z "$job_id" ]; then
    echo "ERROR: Failed to submit job. sbatch output:"
    echo "$submit_output"
    exit 1
fi

echo "Submitted batch job $job_id ($num_nodes nodes)"

# Monitor the job
bash "$SCRIPT_DIR/monitor_slurm_job.sh" "$job_id" "$output_file"
