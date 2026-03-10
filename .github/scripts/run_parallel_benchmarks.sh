#!/bin/bash
# Run PR and master benchmarks in parallel and verify outputs
# Usage: run_parallel_benchmarks.sh <device> <interface> <cluster>

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <device> <interface> <cluster>"
    exit 1
fi

device="$1"
interface="$2"
cluster="$3"

# Get the directory where this script lives (pr/.github/scripts/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=========================================="
echo "Starting parallel benchmark jobs..."
echo "=========================================="

# For Phoenix GPU benchmarks, select a consistent GPU partition before launching
# both parallel jobs so PR and master always land on the same GPU type.
if [ "$device" = "gpu" ] && [ "$cluster" = "phoenix" ]; then
    echo "Selecting Phoenix GPU partition for benchmark consistency..."
    # Require 2 nodes so both PR and master jobs can run concurrently.
    GPU_PARTITION_MIN_NODES=2 source "${SCRIPT_DIR}/select-gpu-partition.sh"
    BENCH_GPU_PARTITION="$SELECTED_GPU_PARTITION"
    export BENCH_GPU_PARTITION
fi

# Run both jobs with monitoring using dedicated script from PR
# Use stdbuf for line-buffered output and prefix each line for clarity
(set -o pipefail; stdbuf -oL -eL bash "${SCRIPT_DIR}/submit_and_monitor_bench.sh" pr "$device" "$interface" "$cluster" 2>&1 | while IFS= read -r line; do echo "[PR] $line"; done) &
pr_pid=$!
echo "PR job started in background (PID: $pr_pid)"

(set -o pipefail; stdbuf -oL -eL bash "${SCRIPT_DIR}/submit_and_monitor_bench.sh" master "$device" "$interface" "$cluster" 2>&1 | while IFS= read -r line; do echo "[MASTER] $line"; done) &
master_pid=$!
echo "Master job started in background (PID: $master_pid)"

echo "Waiting for both jobs to complete..."

# Wait and capture exit codes reliably.
# Use `wait ... || exit=$?` to avoid set -e aborting on the first failure
# (which would orphan the second job).
pr_exit=0
master_exit=0

wait "$pr_pid" || pr_exit=$?
if [ "$pr_exit" -ne 0 ]; then
  echo "PR job exited with code: $pr_exit"
  echo "Last 50 lines of PR job log:"
  tail -n 50 "pr/bench-${device}-${interface}.out" 2>/dev/null || echo "  Could not read PR log"
else
  echo "PR job completed successfully"
fi

wait "$master_pid" || master_exit=$?
if [ "$master_exit" -ne 0 ]; then
  echo "Master job exited with code: $master_exit"
  echo "Last 50 lines of master job log:"
  tail -n 50 "master/bench-${device}-${interface}.out" 2>/dev/null || echo "  Could not read master log"
else
  echo "Master job completed successfully"
fi

# Warn if either job failed (partial results may still be usable)
if [ "${pr_exit}" -ne 0 ] || [ "${master_exit}" -ne 0 ]; then
  echo "WARNING: Benchmark jobs had failures: pr=${pr_exit}, master=${master_exit}"
  echo "Checking for partial results..."
else
  echo "=========================================="
  echo "Both benchmark jobs completed successfully!"
  echo "=========================================="
fi

# Final verification that output files exist before proceeding
pr_yaml="pr/bench-${device}-${interface}.yaml"
master_yaml="master/bench-${device}-${interface}.yaml"

if [ ! -f "$pr_yaml" ]; then
  echo "ERROR: PR benchmark output not found: $pr_yaml"
  ls -la pr/ || true
  echo ""
  echo "Last 100 lines of PR log:"
  tail -n 100 "pr/bench-${device}-${interface}.out" 2>/dev/null || echo "  Could not read PR log"
  exit 1
fi

if [ ! -f "$master_yaml" ]; then
  echo "ERROR: Master benchmark output not found: $master_yaml"
  ls -la master/ || true
  echo ""
  echo "Last 100 lines of master log:"
  tail -n 100 "master/bench-${device}-${interface}.out" 2>/dev/null || echo "  Could not read master log"
  exit 1
fi

echo "Verified both YAML files exist:"
echo "  - $pr_yaml"
echo "  - $master_yaml"

