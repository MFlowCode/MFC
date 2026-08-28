#!/bin/bash
# Run monitor_slurm_job.sh and recover if the monitor is killed (e.g. SIGKILL
# from the runner OS) before the SLURM job completes.  When the monitor exits
# non-zero, sacct is used to verify the job's actual final state; if the SLURM
# job succeeded we exit 0 so the CI step is not falsely marked as failed.
#
# Usage: run_monitored_slurm_job.sh <job_id> <output_file>

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <job_id> <output_file>"
    exit 1
fi

job_id="$1"
output_file="$2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

monitor_exit=0
bash "$SCRIPT_DIR/monitor_slurm_job.sh" "$job_id" "$output_file" || monitor_exit=$?

# Exit 76 = the monitor detected preemption (Phoenix 'embers', PreemptMode=CANCEL).
# Propagate it verbatim so the submit wrapper resubmits a fresh job rather than
# treating a killed-by-preemption job as a genuine test failure.
if [ "$monitor_exit" -eq 76 ]; then
    echo "Monitor reports SLURM job $job_id was PREEMPTED — signaling caller to resubmit."
    exit 76
fi

if [ "$monitor_exit" -ne 0 ]; then
    echo "Monitor exited with code $monitor_exit; re-checking SLURM job $job_id final state..."
    # Give the SLURM epilog time to finalize if the job just finished
    sleep 30
    final_state=$(sacct -j "$job_id" -n -X -P -o State 2>/dev/null | head -n1 | cut -d'|' -f1 | tr -d ' ' || true)
    final_state="${final_state:-UNKNOWN}"
    final_exit=$(sacct -j "$job_id" -X --format=ExitCode --noheader --parsable2 2>/dev/null | head -n1 | tr -d ' ' || true)
    final_exit="${final_exit:-}"
    echo "Final SLURM state=$final_state exit=$final_exit"
    if [ "$final_state" = "PREEMPTED" ]; then
        # Monitor was killed (e.g. SIGKILL) before it could classify this, but
        # the authoritative final state is PREEMPTED — signal a resubmit.
        echo "SLURM job $job_id final state PREEMPTED — signaling caller to resubmit."
        exit 76
    fi
    if [ "$final_state" = "COMPLETED" ] && [ "$final_exit" = "0:0" ]; then
        echo "SLURM job $job_id completed successfully despite monitor failure — continuing."
    else
        echo "ERROR: SLURM job $job_id did not complete successfully (state=$final_state exit=$final_exit)"
        exit 1
    fi
fi
