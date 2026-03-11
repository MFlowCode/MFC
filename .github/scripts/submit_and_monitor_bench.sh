#!/bin/bash
# Submit and monitor a benchmark job on a SLURM cluster
# Usage: submit_and_monitor_bench.sh <dir> <device> <interface> <cluster>

set -euo pipefail

if [ $# -ne 4 ]; then
    echo "Usage: $0 <dir> <device> <interface> <cluster>"
    exit 1
fi

dir="$1"
device="$2"
interface="$3"
cluster="$4"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[$dir] Submitting benchmark for $device-$interface on $cluster..."
cd "$dir"

# Use the PR's submit-slurm-job.sh and bench script for both master and PR jobs.
# The bench script must come from the PR tree (master may not have common/bench.sh
# yet), and the script only orchestrates build+bench — the actual MFC code under
# test is the cwd's checkout (master/ or pr/).
PR_BENCH_SCRIPT="$(cd "${SCRIPT_DIR}/../workflows/common" && pwd)/bench.sh"
bash "${SCRIPT_DIR}/submit-slurm-job.sh" "$PR_BENCH_SCRIPT" "$device" "$interface" "$cluster"

# Verify the YAML output file was created
job_slug="bench-$device-$interface"
yaml_file="${job_slug}.yaml"
if [ ! -f "$yaml_file" ]; then
    echo "[$dir] ERROR: Expected output file not found: $yaml_file"
    echo "[$dir] Directory contents:"
    ls -la *.yaml 2>/dev/null || echo "  No YAML files found"
    echo ""
    output_file="${job_slug}.out"
    echo "[$dir] Last 100 lines of job output ($output_file):"
    echo "----------------------------------------"
    tail -n 100 "$output_file" 2>/dev/null || echo "  Could not read output file"
    echo "----------------------------------------"
    exit 1
fi

echo "[$dir] Verified output file exists: $yaml_file ($(stat -f%z "$yaml_file" 2>/dev/null || stat -c%s "$yaml_file" 2>/dev/null) bytes)"
