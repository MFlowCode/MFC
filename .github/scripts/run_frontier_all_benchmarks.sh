#!/bin/bash
# Orchestrate all Frontier benchmark configs in one multi-node SLURM allocation.
# 1. Builds all configs on the login node (PR and master, in parallel)
# 2. Submits a single SLURM job running benchmarks in parallel via ssh

set -euo pipefail

# Ignore SIGHUP to survive login node session drops
trap '' HUP

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# SLURM parameters
SLURM_ACCOUNT="ENG160"
SLURM_PARTITION="extended"
SLURM_WALLTIME="05:59:00"
CONFIG_TIMEOUT=7200  # 120 min per config

# Benchmark configs: version cluster device interface
# 6 total: 3 configs x 2 versions (PR + master)
configs=(
    "pr     frontier     gpu acc"
    "pr     frontier     gpu omp"
    "pr     frontier_amd gpu omp"
    "master frontier     gpu acc"
    "master frontier     gpu omp"
    "master frontier_amd gpu omp"
)
num_nodes=${#configs[@]}

echo "=========================================="
echo "Frontier consolidated benchmarks: $num_nodes configs on $num_nodes nodes"
echo "=========================================="

# Write config file for sbatch to read (single source of truth)
config_file="frontier-bench-configs.txt"
printf '%s\n' "${configs[@]}" > "$config_file"

# --- Phase 1: Create per-config source copies ---
for cfg in "${configs[@]}"; do
    read -r version cluster device interface <<< "$cfg"
    dir="${version}-${cluster}-${device}-${interface}"
    echo "Creating source copy: $dir from $version/"
    rm -rf "$dir"
    cp -al "$version" "$dir" 2>/dev/null || cp -r "$version" "$dir"
done

# --- Phase 2: Build all configs on login node in parallel ---
# Avoid setuptools_scm git conflicts in hardlink copies (shared .git/index)
export SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0

MAX_PARALLEL=2

echo ""
echo "=========================================="
echo "Starting parallel builds (${num_nodes} configs, max $MAX_PARALLEL concurrent)..."
echo "=========================================="
build_pids=()
running=()
for i in "${!configs[@]}"; do
    # Wait until a build slot is available
    throttle_timer=0
    while [ ${#running[@]} -ge $MAX_PARALLEL ]; do
        sleep 5
        throttle_timer=$((throttle_timer + 5))
        still_running=()
        for pid in "${running[@]}"; do
            kill -0 "$pid" 2>/dev/null && still_running+=("$pid")
        done
        running=("${still_running[@]}")
        if [ $throttle_timer -ge 120 ]; then
            throttle_timer=0
            echo "--- Build heartbeat ($(date +%H:%M:%S)) ---"
            for j in $(seq 0 $((i - 1))); do
                read -r v c d iface <<< "${configs[$j]}"
                if kill -0 "${build_pids[$j]}" 2>/dev/null; then
                    last=$(tail -n 1 "build-${v}-${c}-${d}-${iface}.log" 2>/dev/null | head -c 120 || echo "")
                    echo "  $v $c $d $iface: $last"
                fi
            done
        fi
    done

    read -r version cluster device interface <<< "${configs[$i]}"
    dir="${version}-${cluster}-${device}-${interface}"
    log="build-${version}-${cluster}-${device}-${interface}.log"
    echo "  Starting: $version $cluster $device $interface"
    (
        cd "$dir"
        bash .github/workflows/${cluster}/build.sh "$device" "$interface" bench
    ) > "$log" 2>&1 &
    build_pids[$i]=$!
    running+=($!)
done

# Periodic heartbeat while builds run
(
    while true; do
        sleep 120
        alive=0
        for pid in "${build_pids[@]}"; do
            kill -0 "$pid" 2>/dev/null && alive=$((alive + 1))
        done
        [ "$alive" -eq 0 ] && break
        echo "--- Build heartbeat ($(date +%H:%M:%S)): $alive/${#build_pids[@]} running ---"
        for i in "${!configs[@]}"; do
            read -r version cluster device interface <<< "${configs[$i]}"
            log="build-${version}-${cluster}-${device}-${interface}.log"
            if kill -0 "${build_pids[$i]}" 2>/dev/null; then
                size=$(stat -c%s "$log" 2>/dev/null || echo 0)
                last=$(tail -n 1 "$log" 2>/dev/null | head -c 120 || echo "")
                echo "  $version $cluster $device $interface: running (${size} bytes) $last"
            fi
        done
    done
) &
heartbeat_pid=$!

# Wait for all builds and report results
build_failed=0
build_exits=()
for i in "${!build_pids[@]}"; do
    read -r version cluster device interface <<< "${configs[$i]}"
    if wait "${build_pids[$i]}"; then
        build_exits+=(0)
    else
        code=$?
        build_exits+=($code)
        build_failed=1
    fi
done

# Stop heartbeat
kill "$heartbeat_pid" 2>/dev/null || true; wait "$heartbeat_pid" 2>/dev/null || true

# Print build logs: passed builds collapsed, failed builds in full
for i in "${!configs[@]}"; do
    read -r version cluster device interface <<< "${configs[$i]}"
    log="build-${version}-${cluster}-${device}-${interface}.log"
    if [ "${build_exits[$i]}" -eq 0 ]; then
        echo "::group::Build PASSED: $version $cluster $device $interface"
        cat "$log"
        echo "::endgroup::"
    else
        echo "=========================================="
        echo "Build FAILED: $version $cluster $device $interface (exit ${build_exits[$i]})"
        echo "=========================================="
        cat "$log"
    fi
done

# Abort on failure
if [ "$build_failed" -ne 0 ]; then
    echo ""
    echo "=========================================="
    echo "Build failures detected — see logs above."
    echo "=========================================="
    exit 1
fi

echo ""
echo "=========================================="
echo "All builds complete. Submitting ${num_nodes}-node SLURM job..."
echo "=========================================="

# --- Phase 3: Submit one sbatch job with N nodes ---
output_file="bench-frontier-all.out"

submit_output=$(sbatch <<OUTER
#!/bin/bash
#SBATCH -J MFC-frontier-all-bench
#SBATCH -A $SLURM_ACCOUNT
#SBATCH -N $num_nodes
#SBATCH -t $SLURM_WALLTIME
#SBATCH -o$output_file
#SBATCH -p $SLURM_PARTITION

set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in \$(pwd)"
echo "Allocated nodes: \$SLURM_NODELIST"

# Get list of individual node hostnames
mapfile -t nodes < <(scontrol show hostnames "\$SLURM_NODELIST")
echo "Nodes: \${nodes[*]}"

# Read config table from file (written by outer script, avoids duplication)
mapfile -t configs < "$config_file"

pids=()

cleanup() {
    echo "Cleaning up — killing all remote processes..."
    for pid in "\${pids[@]}"; do
        kill "\$pid" 2>/dev/null
    done
    wait
}
trap cleanup EXIT

for i in "\${!configs[@]}"; do
    read -r version cluster device interface <<< "\${configs[\$i]}"
    node="\${nodes[\$i]}"
    dir="\${version}-\${cluster}-\${device}-\${interface}"
    outfile="\${dir}/bench-\${device}-\${interface}.out"

    echo "[\$node] Starting bench: \$version \$cluster \$device \$interface in \$dir"

    timeout $CONFIG_TIMEOUT ssh -q -o StrictHostKeyChecking=no "\$node" \
        "cd \$SLURM_SUBMIT_DIR/\$dir && bash .github/scripts/frontier_bench_config.sh \$cluster \$device \$interface" \
        > "\$outfile" 2>&1 &
    pids+=(\$!)
done

echo "All bench configs launched, waiting for completion..."

# Wait for all and collect exit codes
overall_exit=0
for i in "\${!pids[@]}"; do
    read -r version cluster device interface <<< "\${configs[\$i]}"
    pid=\${pids[\$i]}
    if wait "\$pid"; then
        echo "PASSED: \$version \$cluster \$device \$interface (PID \$pid)"
    else
        code=\$?
        echo "FAILED: \$version \$cluster \$device \$interface (PID \$pid, exit code \$code)"
        overall_exit=1
    fi
done

# Print summary
echo ""
echo "=========================================="
echo "Benchmark summary:"
for cfg in "\${configs[@]}"; do
    read -r version cluster device interface <<< "\$cfg"
    dir="\${version}-\${cluster}-\${device}-\${interface}"
    yaml="\${dir}/bench-\${device}-\${interface}.yaml"
    if [ -f "\$yaml" ]; then
        echo "  \$version \$cluster \$device \$interface: OK (\$(stat -c%s "\$yaml" 2>/dev/null) bytes)"
    else
        echo "  \$version \$cluster \$device \$interface: MISSING YAML"
    fi
done
echo "=========================================="

exit \$overall_exit
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
