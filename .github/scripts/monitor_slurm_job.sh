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

# Put the one thing a reader needs on the run's summary page. Without this,
# learning why a job failed means opening a log of tens of thousands of lines --
# and an infrastructure fault looks exactly like a test failure until you do.
# Silent when not running under Actions.
ci_summary() {
  [ -n "${GITHUB_STEP_SUMMARY:-}" ] || return 0
  printf '%b\n' "$1" >> "$GITHUB_STEP_SUMMARY"
}

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

# Check if a state is terminal (job is done, for better or worse).
# PREEMPTED is handled separately (below): Phoenix preempts 'embers' jobs with
# PreemptMode=CANCEL, not REQUEUE (verified via `scontrol show config`), so a
# preempted job is killed outright and never restarts under the same ID.
# --requeue is a no-op for it. It is surfaced via PREEMPT_EXIT so the submit
# wrapper can resubmit a fresh job instead of failing the CI step.
PREEMPT_EXIT=76
is_terminal_state() {
  case "$1" in
    COMPLETED|FAILED|CANCELLED|CANCELLED+|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|DEADLINE|REVOKED)
      return 0 ;;
    *)
      return 1 ;;
  esac
}

# Optionally bound how long a job may sit un-started in the queue. On the
# preemptible Phoenix 'embers' QOS a job routinely stays PENDING for hours and
# needs most of the job-level `timeout-minutes` (480m) window to backfill onto a
# free node; that job timeout is the real backstop. Default to 0 (wait
# indefinitely, up to the job timeout) so ordinary queue pressure does not turn
# otherwise-healthy jobs into red CI. Set SLURM_MAX_QUEUE_SECONDS>0 to opt into
# an earlier queue-starvation cutoff where the scheduler is not preemptible.
: "${SLURM_MAX_QUEUE_SECONDS:=0}"   # 0 = wait indefinitely (job timeout is the backstop)
# Reject a non-integer override rather than silently skipping the budget.
if ! [[ "$SLURM_MAX_QUEUE_SECONDS" =~ ^[0-9]+$ ]]; then
  echo "ERROR: SLURM_MAX_QUEUE_SECONDS must be a non-negative integer (seconds), got '$SLURM_MAX_QUEUE_SECONDS'" >&2
  exit 1
fi
# How long to wait between status polls and between output-stabilization
# checks. Overridable so tests can exercise this script without sleeping
# through it; CI leaves it at the default.
: "${MFC_MONITOR_POLL_SECONDS:=5}"

queue_start=$(date +%s)

abort_queue_starvation() {
  local waited="$1"
  echo "##[error]SLURM job $job_id did not start within ${waited}s (SLURM_MAX_QUEUE_SECONDS=$SLURM_MAX_QUEUE_SECONDS)."
  echo "QUEUE STARVATION: the cluster scheduler could not start this job in time."
  echo "This is an infrastructure / queue-availability problem, NOT a code or test failure."
  echo "Cancelling the queued job so it does not keep holding a CI runner slot."
  scancel "$job_id" 2>/dev/null || true
  exit 75   # EX_TEMPFAIL — distinguishes queue starvation from a real test failure
}

# Wait for file to appear, using robust state checking.
# Never give up due to transient squeue/sacct failures — the queue-wait budget
# above (or the CI job timeout) is the ultimate backstop.
echo "Waiting for job to start..."
unknown_count=0
while [ ! -f "$output_file" ]; do
  state=$(get_job_state "$job_id")

  # A started job (RUNNING/COMPLETING) whose output file is merely NFS-delayed
  # is exempt, so work in progress is never killed here.
  if [ "$SLURM_MAX_QUEUE_SECONDS" -gt 0 ]; then
    case "$state" in
      RUNNING|COMPLETING) ;;
      *)
        waited=$(( $(date +%s) - queue_start ))
        if [ "$waited" -ge "$SLURM_MAX_QUEUE_SECONDS" ]; then
          abort_queue_starvation "$waited"
        fi
        ;;
    esac
  fi

  case "$state" in
    PREEMPTED)
      # Preempted before producing output (embers, PreemptMode=CANCEL): the job
      # is dead and will not requeue. Signal the caller to resubmit a fresh job.
      echo "[$(date +%H:%M:%S)] Job $job_id PREEMPTED before start/output — signaling resubmit."
      exit "$PREEMPT_EXIT"
      ;;
    PENDING|CONFIGURING)
      unknown_count=0
      sleep "$MFC_MONITOR_POLL_SECONDS"
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
      sleep "$MFC_MONITOR_POLL_SECONDS"
      ;;
    *)
      # Terminal state — job finished without creating output
      if is_terminal_state "$state"; then
        echo "ERROR: Job $job_id reached terminal state ($state) without creating output file"
        exit 1
      fi
      # Unrecognized state, keep waiting
      sleep "$MFC_MONITOR_POLL_SECONDS"
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

  if [ "$state" = "PREEMPTED" ]; then
    # Preempted mid-run (embers, PreemptMode=CANCEL): dead, will not requeue.
    echo "[$(date +%H:%M:%S)] Job $job_id PREEMPTED mid-run — signaling resubmit."
    exit "$PREEMPT_EXIT"
  fi

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

  sleep "$MFC_MONITOR_POLL_SECONDS"
done

# Give tail a moment to flush the final lines, then stop streaming. Whether it
# was still alive decides how much needs reprinting below: if it streamed the
# whole job, printing the file again just doubles every log.
sleep 2
streamed_ok=0
if kill -0 "${tail_pid}" 2>/dev/null; then
  streamed_ok=1
fi
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
    sleep "$MFC_MONITOR_POLL_SECONDS"
  done
fi

# Reprint only what streaming may have missed. `tail -f` above already emitted
# the whole file as it was written, so cat'ing it again duplicated every job's
# output -- measured at 3 copies of each line on a GPU job, and 65,000 lines of
# offload diagnostics repeated for a single fault. The reprint exists solely as
# a safety net for a tail that died mid-job, so it is bounded when tail survived
# and complete only when it did not.
echo ""
if [ "${streamed_ok:-0}" -eq 1 ]; then
  echo "=== Final output (tail; the full log streamed above) ==="
  tail -n "${MFC_MONITOR_FINAL_LINES:-40}" "$output_file"
else
  echo "=== Final output (streaming stopped early; reprinting in full) ==="
  cat "$output_file"
fi

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

# The preflight's node-fault verdict comes back as the job's own exit code.
# Relay it verbatim: flattening it to 1 would leave the submit wrapper unable to
# tell "this node is unusable" (exclude it and try again) from "the tests
# failed" (report it).
faulted_node=$(grep -oE 'MFC_FAULT_NODE=[^ ]+' "$output_file" 2>/dev/null | tail -n1 | cut -d= -f2 || true)

case "$exit_code" in
  77:*)
    echo "Job $job_id failed preflight: the node is unusable — signaling caller to exclude it and resubmit."
    ci_summary "### :warning: Infrastructure fault — not a code or test failure\n\nNode \`${faulted_node:-unknown}\` could not run MFC (job \`$job_id\`). It is excluded and the job resubmitted elsewhere.\n"
    monitor_success=1
    exit 77
    ;;
esac

# Check if job succeeded
if [ "$exit_code" != "0:0" ]; then
  echo "ERROR: Job $job_id failed with exit code $exit_code"
  # A GPU memory fault explains itself in a block the test harness prints; lift
  # it onto the summary page so the faulting kernel and source line are visible
  # without opening the log at all.
  if grep -q 'GPU fault summary' "$output_file" 2>/dev/null; then
    ci_summary "### GPU memory fault\n\n\`\`\`\n$(grep -A6 'GPU fault summary' "$output_file" | head -8 | sed 's/`/'"'"'/g')\n\`\`\`\n"
  else
    ci_summary "### Job \`$job_id\` failed (exit $exit_code)\n\n\`\`\`\n$(tail -n 15 "$output_file" | sed 's/`/'"'"'/g')\n\`\`\`\n"
  fi
  exit 1
fi

monitor_success=1
echo "Job $job_id completed successfully"
exit 0
