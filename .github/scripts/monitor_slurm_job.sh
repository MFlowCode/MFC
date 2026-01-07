#!/bin/bash
# Monitor a SLURM job and stream its output in real-time
# Usage: monitor_slurm_job.sh <job_id> <output_file>

set -euo pipefail

# Cleanup handler to prevent orphaned tail processes
cleanup() {
  if [ -n "${tail_pid:-}" ]; then
    kill "${tail_pid}" 2>/dev/null || true
  fi
}
trap cleanup EXIT

if [ $# -ne 2 ]; then
    echo "Usage: $0 <job_id> <output_file>"
    exit 1
fi

job_id="$1"
output_file="$2"

echo "Submitted batch job $job_id"
echo "Monitoring output file: $output_file"

# Wait for file to appear with retry logic for transient squeue failures
echo "Waiting for job to start..."
squeue_retries=0
max_squeue_retries=5
while [ ! -f "$output_file" ]; do
  # Check if job is still queued/running
  if squeue -j "$job_id" &>/dev/null; then
    squeue_retries=0  # Reset on success
    sleep 5
  else
    squeue_retries=$((squeue_retries + 1))
    if [ $squeue_retries -ge $max_squeue_retries ]; then
      # Job not in queue and output file doesn't exist
      if [ ! -f "$output_file" ]; then
        echo "ERROR: Job $job_id not in queue and output file not created"
        exit 1
      fi
      break
    fi
    # Exponential backoff
    sleep_time=$((2 ** squeue_retries))
    echo "Warning: squeue check failed, retrying in ${sleep_time}s..."
    sleep $sleep_time
  fi
done

echo "=== Streaming output for job $job_id ==="

# Start tail and redirect its output to file descriptor 3 for multiplexing
# This allows us to stream tail output while also printing heartbeat messages
exec 3< <(stdbuf -oL -eL tail -f "$output_file" 2>&1)
tail_pid=$!

# Monitor job status and stream output simultaneously
squeue_failures=0
last_heartbeat=$(date +%s)

while true; do
  # Try to read from tail output (non-blocking via timeout)
  # Read multiple lines if available to avoid falling behind
  lines_read=0
  while IFS= read -r -t 0.1 line <&3 2>/dev/null; do
    echo "$line"
    lines_read=$((lines_read + 1))
    last_heartbeat=$(date +%s)
    # Limit burst reads to avoid starving the status check
    if [ $lines_read -ge 100 ]; then
      break
    fi
  done
  
  # Check job status
  current_time=$(date +%s)
  if ! squeue -j "$job_id" &>/dev/null; then
    squeue_failures=$((squeue_failures + 1))
    # Check if job actually completed using sacct (if available)
    if [ $squeue_failures -ge 3 ]; then
      if command -v sacct >/dev/null 2>&1; then
        state=$(sacct -j "$job_id" --format=State --noheader 2>/dev/null | head -n1 | awk '{print $1}')
        # Consider job done only if it reached a terminal state
        case "$state" in
          COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY)
            echo "[$(date +%H:%M:%S)] Job $job_id reached terminal state: $state"
            break
            ;;
          *)
            # treat as transient failure, reset failures and continue polling
            squeue_failures=0
            ;;
        esac
      else
        # No sacct: assume job completed after 3 failures
        echo "[$(date +%H:%M:%S)] Job $job_id no longer in queue"
        break
      fi
    fi
  else
    squeue_failures=0
    # Print heartbeat if no output for 60 seconds
    if [ $((current_time - last_heartbeat)) -ge 60 ]; then
      echo "[$(date +%H:%M:%S)] Job $job_id still running (no new output for 60s)..."
      last_heartbeat=$current_time
    fi
  fi
  
  # Sleep briefly between status checks
  sleep 1
done

# Drain any remaining output from tail after job completes
echo "Draining remaining output..."
drain_count=0
while IFS= read -r -t 0.5 line <&3 2>/dev/null; do
  echo "$line"
  drain_count=$((drain_count + 1))
  # Safety limit to avoid infinite loop
  if [ $drain_count -ge 10000 ]; then
    echo "Warning: Truncating remaining output after 10000 lines"
    break
  fi
done

# Close the file descriptor and kill tail
exec 3<&-
kill "${tail_pid}" 2>/dev/null || true

# Wait for output file to finish growing (stabilize) before stopping tail
if [ -f "$output_file" ]; then
  last_size=-1
  same_count=0
  while true; do
    size=$(stat -c%s "$output_file" 2>/dev/null || echo -1)
    if [ "$size" -eq "$last_size" ] && [ "$size" -ge 0 ]; then
      same_count=$((same_count + 1))
    else
      same_count=0
      last_size=$size
    fi
    # two consecutive stable checks (~10s) implies file likely flushed
    if [ $same_count -ge 2 ]; then
      break
    fi
    sleep 5
  done
fi

# Stop tailing (trap will also handle this on exit)
kill "${tail_pid}" 2>/dev/null || true

echo ""
echo "=== Final output ==="
cat "$output_file"

# Check exit status with sacct fallback
exit_code=""

# Try scontrol first (works for recent jobs)
scontrol_output=$(scontrol show job "$job_id" 2>/dev/null || echo "")
if [ -n "$scontrol_output" ]; then
  exit_code=$(echo "$scontrol_output" | grep -oE 'ExitCode=[0-9]+:[0-9]+' | cut -d= -f2 || echo "")
fi

# If scontrol failed or returned invalid job, try sacct (for completed/aged-out jobs)
if [ -z "$exit_code" ]; then
  echo "Warning: scontrol failed to get exit code, trying sacct..."
  sacct_output=$(sacct -j "$job_id" --format=ExitCode --noheader --parsable2 2>/dev/null | head -n1 || echo "")
  if [ -n "$sacct_output" ]; then
    exit_code="$sacct_output"
  fi
fi

# If we still can't determine exit code, fail explicitly
if [ -z "$exit_code" ]; then
  echo "ERROR: Unable to determine exit status for job $job_id"
  echo "Both scontrol and sacct failed to return valid exit code"
  exit 1
fi

# Check if job succeeded
if [ "$exit_code" != "0:0" ]; then
  echo "ERROR: Job $job_id failed with exit code $exit_code"
  exit 1
fi

echo "Job $job_id completed successfully"
exit 0

