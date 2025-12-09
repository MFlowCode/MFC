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

echo "[$dir] Submitting benchmark for $device-$interface on $cluster..."
cd "$dir"

# Submit job
submit_output=$(bash .github/workflows/$cluster/submit-bench.sh \
  .github/workflows/$cluster/bench.sh "$device" "$interface" 2>&1)

job_id=$(echo "$submit_output" | sed -n 's/.*Submitted batch job \([0-9][0-9]*\).*/\1/p')
job_slug="bench-$device-$interface"
output_file="${job_slug}.out"

if [ -z "$job_id" ]; then
  echo "[$dir] ERROR: Failed to submit job"
  echo "$submit_output"
  exit 1
fi

echo "[$dir] Job ID: $job_id, monitoring output file: $output_file"

# Use the monitoring script
bash .github/scripts/monitor_slurm_job.sh "$job_id" "$output_file"

echo "[$dir] Monitoring complete for job $job_id"

# Verify the YAML output file was created
yaml_file="${job_slug}.yaml"
if [ ! -f "$yaml_file" ]; then
  echo "[$dir] ERROR: Expected output file not found: $yaml_file"
  echo "[$dir] Directory contents:"
  ls -la *.yaml 2>/dev/null || echo "  No YAML files found"
  exit 1
fi

echo "[$dir] Verified output file exists: $yaml_file ($(stat -f%z "$yaml_file" 2>/dev/null || stat -c%s "$yaml_file" 2>/dev/null) bytes)"

