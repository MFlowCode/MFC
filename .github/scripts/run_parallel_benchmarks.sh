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

# Run both jobs with monitoring using dedicated script from PR
(bash "${SCRIPT_DIR}/submit_and_monitor_bench.sh" pr "$device" "$interface" "$cluster") &
pr_pid=$!
echo "PR job started in background (PID: $pr_pid)"

(bash "${SCRIPT_DIR}/submit_and_monitor_bench.sh" master "$device" "$interface" "$cluster") &
master_pid=$!
echo "Master job started in background (PID: $master_pid)"

echo "Waiting for both jobs to complete..."

# Wait and capture exit codes reliably
pr_exit=0
master_exit=0

if ! wait "$pr_pid"; then
  pr_exit=$?
  echo "PR job exited with code: $pr_exit"
else
  echo "PR job completed successfully"
fi

if ! wait "$master_pid"; then
  master_exit=$?
  echo "Master job exited with code: $master_exit"
else
  echo "Master job completed successfully"
fi

# Check if either job failed
if [ "${pr_exit}" -ne 0 ] || [ "${master_exit}" -ne 0 ]; then
  echo "ERROR: One or both benchmark jobs failed: pr_exit=${pr_exit}, master_exit=${master_exit}"
  exit 1
fi

echo "=========================================="
echo "Both benchmark jobs completed successfully!"
echo "=========================================="

# Final verification that output files exist before proceeding
pr_yaml="pr/bench-${device}-${interface}.yaml"
master_yaml="master/bench-${device}-${interface}.yaml"

if [ ! -f "$pr_yaml" ]; then
  echo "ERROR: PR benchmark output not found: $pr_yaml"
  ls -la pr/ || true
  exit 1
fi

if [ ! -f "$master_yaml" ]; then
  echo "ERROR: Master benchmark output not found: $master_yaml"
  ls -la master/ || true
  exit 1
fi

echo "Verified both YAML files exist:"
echo "  - $pr_yaml"
echo "  - $master_yaml"

