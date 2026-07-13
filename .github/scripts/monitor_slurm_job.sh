#!/bin/bash
# Monitor a SLURM job and stream its output in real-time
# Usage: monitor_slurm_job.sh <job_id> <output_file>

set -euo pipefail

# Cleanup handler to prevent orphaned tail processes and cancel orphaned jobs
cleanup() {
  if [ -n "${tail_pid:-}" ]; then
    kill "${tail_pid}" 2>/dev/null || true
  fi
  # Cancel the SLURM job only if it is still active in the scheduler.
  # If the job already left the queue (squeue returns empty), it has finished
  # and run_monitored_slurm_job.sh will recover via sacct — don't cancel it.
  if [ "${monitor_success:-0}" -ne 1 ] && [ -n "${job_id:-}" ]; then
    active_state=$(squeue -j "$job_id" -h -o '%T' 2>/dev/null | head -n1 | tr -d ' ' || echo "")
    if [ -n "$active_state" ]; then
      echo "Monitor exiting abnormally — cancelling SLURM job $job_id (state: $active_state)"
      scancel "$job_id" 2>/dev/null || true
    else
      echo "Monitor exiting abnormally — SLURM job $job_id already left queue, not cancelling"
    fi
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
    # When a job is preempted+requeued, sacct -X reports PREEMPTED for the
    # original attempt while the requeued run may have completed.  Check all
    # records (without -X) for a terminal state that supersedes PREEMPTED.
    if [ "$state" = "PREEMPTED" ]; then
      requeue_state=$(sacct -j "$jid" -n -P -o State 2>/dev/null | grep -v PREEMPTED | head -n1 | cut -d'|' -f1 || true)
      if [ -n "$requeue_state" ]; then
        state="$requeue_state"
      fi
    fi
    if [ -n "$state" ]; then
      echo "$state"
      return
    fi
  fi

  echo "UNKNOWN"
}

# Check if a state is terminal (job is done, for better or worse)
# PREEMPTED is intentionally excluded: with --requeue the job restarts under
# the same job ID and we must keep monitoring rather than exiting early.
is_terminal_state() {
  case "$1" in
    COMPLETED|FAILED|CANCELLED|CANCELLED+|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|DEADLINE|REVOKED)
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
    PENDING|CONFIGURING|PREEMPTED)
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

# Stream the job's output to the step log with a plain backgrounded `tail`.
# This previously used a timed read (`read -t 1`) over a `tail -f` process
# substitution; when that pipe broke, bash could crash in the read-builtin's
# alarm/longjmp unwind path (SIGSEGV), and the EXIT trap would then cancel a
# still-healthy job. A backgrounded tail avoids that construct entirely, and
# the final `cat` below reprints the whole file so nothing is lost if tail is
# killed mid-flush.
stdbuf -oL -eL tail -f "$output_file" 2>&1 &
tail_pid=$!

# Poll job status until it reaches a terminal state; streaming happens
# independently in the background tail above.
last_heartbeat=$(date +%s)
while true; do
  state=$(get_job_state "$job_id")

  if is_terminal_state "$state"; then
    echo "[$(date +%H:%M:%S)] Job $job_id reached terminal state: $state"
    break
  fi

  # Periodic heartbeat so the CI log never looks stalled during quiet phases.
  current_time=$(date +%s)
  if [ $((current_time - last_heartbeat)) -ge 60 ]; then
    echo "[$(date +%H:%M:%S)] Job $job_id state=$state..."
    last_heartbeat=$current_time
  fi

  sleep 5
done

# Give tail a moment to flush the final lines, then stop streaming.
sleep 2
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
