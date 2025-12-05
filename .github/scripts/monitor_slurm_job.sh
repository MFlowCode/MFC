#!/bin/bash
# Monitor a SLURM job and stream its output in real-time
# Usage: monitor_slurm_job.sh <job_id> <output_file>

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <job_id> <output_file>"
    exit 1
fi

job_id="$1"
output_file="$2"

echo "Submitted batch job $job_id"
echo "Monitoring output file: $output_file"

# Wait for file to appear (check job status if it takes a while)
echo "Waiting for job to start..."
while [ ! -f "$output_file" ]; do
  # Check if job failed to start
  if ! squeue -j "$job_id" &>/dev/null && [ ! -f "$output_file" ]; then
    echo "ERROR: Job $job_id finished without creating output file"
    exit 1
  fi
  sleep 5
done

echo "=== Streaming output for job $job_id ==="
# Stream output while job runs
tail -f "$output_file" &
tail_pid=$!

# Wait for job to complete with retry logic for transient squeue failures
squeue_failures=0
while true; do
  if squeue -j "$job_id" &>/dev/null; then
    squeue_failures=0
  else
    squeue_failures=$((squeue_failures + 1))
    # Allow a few transient failures before concluding job is done
    if [ $squeue_failures -ge 3 ]; then
      break
    fi
  fi
  sleep 5
done

# Stop tailing
kill $tail_pid 2>/dev/null || true

echo ""
echo "=== Final output ==="
cat "$output_file"

# Check exit status
exit_code=$(scontrol show job "$job_id" 2>/dev/null | grep -oP 'ExitCode=\K[0-9]+:[0-9]+' || echo "0:0")
if [ "$exit_code" != "0:0" ]; then
  echo "ERROR: Job $job_id failed with exit code $exit_code"
  exit 1
fi

echo "Job $job_id completed successfully"
exit 0

