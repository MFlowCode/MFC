#!/bin/bash
# Monitor a SLURM job and stream its output in real-time
# Usage: monitor_slurm_job.sh <job_id> <output_file>

set -euo pipefail

# Cleanup handler to prevent orphaned tail processes and cancel orphaned jobs
cleanup() {
  if [ -n "${tail_pid:-}" ]; then
    kill "${tail_pid}" 2>/dev/null || true
  fi
  # Cancel the SLURM job if the monitor is exiting due to an error
  # (e.g., the CI runner is being killed). Don't cancel on success.
  if [ "${monitor_success:-0}" -ne 1 ] && [ -n "${job_id:-}" ]; then
    echo "Monitor exiting abnormally — cancelling SLURM job $job_id"
    scancel "$job_id" 2>/dev/null || true
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

# Robustly check SLURM job state using squeue with sacct fallback.
# Returns the state string (PENDING, RUNNING, COMPLETED, FAILED, etc.)
# or "UNKNOWN" if both commands fail.
get_job_state() {
  local jid="$1"
  local state

  # Try squeue first (fast, works for active jobs)
  state=$(squeue -j "$jid" -h -o '%T' 2>/dev/null | head -n1 | tr -d ' ' || true)
  if [ -n "$state" ]; then
    echo "$state"
    return
  fi

  # Fallback to sacct (works for completed/historical jobs)
  if command -v sacct >/dev/null 2>&1; then
    state=$(sacct -j "$jid" -n -X -P -o State 2>/dev/null | head -n1 | cut -d'|' -f1 || true)
    if [ -n "$state" ]; then
      echo "$state"
      return
    fi
  fi

  echo "UNKNOWN"
}

# Check if a state is terminal (job is done, for better or worse)
is_terminal_state() {
  case "$1" in
    COMPLETED|FAILED|CANCELLED|CANCELLED+|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|DEADLINE)
      return 0 ;;
    *)
      return 1 ;;
  esac
}

# Wait for file to appear, using robust state checking.
# Never give up due to transient squeue/sacct failures — the CI job timeout
# is the ultimate backstop.
echo "Waiting for job to start..."
unknown_count=0
while [ ! -f "$output_file" ]; do
  state=$(get_job_state "$job_id")

  case "$state" in
    PENDING|CONFIGURING)
      unknown_count=0
      sleep 5
      ;;
    RUNNING|COMPLETING)
      unknown_count=0
      # Job is running but output file not yet visible (NFS delay)
      sleep 2
      ;;
    UNKNOWN)
      unknown_count=$((unknown_count + 1))
      # Only print warning periodically to avoid log spam
      if [ $((unknown_count % 12)) -eq 1 ]; then
        echo "Warning: Could not query job $job_id state (SLURM may be temporarily unavailable)..."
      fi
      sleep 5
      ;;
    *)
      # Terminal state — job finished without creating output
      if is_terminal_state "$state"; then
        echo "ERROR: Job $job_id reached terminal state ($state) without creating output file"
        exit 1
      fi
      # Unrecognized state, keep waiting
      sleep 5
      ;;
  esac
done

echo "=== Streaming output for job $job_id ==="

# Start tail and redirect its output to file descriptor 3 for multiplexing
# This allows us to stream tail output while also printing heartbeat messages
exec 3< <(stdbuf -oL -eL tail -f "$output_file" 2>&1)
tail_pid=$!

# Monitor job status and stream output simultaneously
last_heartbeat=$(date +%s)

while true; do
  # Try to read from tail output (non-blocking via timeout)
  # Read multiple lines if available to avoid falling behind
  lines_read=0
  while IFS= read -r -t 1 line <&3 2>/dev/null; do
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
  state=$(get_job_state "$job_id")

  if is_terminal_state "$state"; then
    echo "[$(date +%H:%M:%S)] Job $job_id reached terminal state: $state"
    break
  else
    # Print heartbeat if no output for 60 seconds
    if [ $((current_time - last_heartbeat)) -ge 60 ]; then
      echo "[$(date +%H:%M:%S)] Job $job_id state=$state (no new output for 60s)..."
      last_heartbeat=$current_time
    fi
  fi

  # Sleep briefly between status checks
  sleep 1
done

# Drain any remaining output from tail after job completes
echo "Draining remaining output..."
drain_count=0
while IFS= read -r -t 1 line <&3 2>/dev/null; do
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
tail_pid=""

# Wait for output file to stabilize (NFS flush) before final read
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

monitor_success=1
echo "Job $job_id completed successfully"
exit 0
