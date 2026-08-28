#!/bin/bash
# Run PR and master benchmarks and verify outputs.
# Both SLURM jobs are submitted up front so they run concurrently on
# compute nodes (fair comparison under the same cluster load), but
# monitoring happens sequentially to stay within the per-user cgroup
# memory limit on login nodes (4 GB on Phoenix shared by 7 runners).
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

echo "===================="
echo "Starting benchmark jobs..."
echo "===================="

# For Phoenix GPU benchmarks, select a consistent GPU partition so PR and
# master always land on the same GPU type.
if [ "$device" = "gpu" ] && [ "$cluster" = "phoenix" ]; then
    echo "Selecting Phoenix GPU partition for benchmark consistency..."
    # Require 2 nodes so both jobs can run concurrently on compute.
    GPU_PARTITION_MIN_NODES=2 source "${SCRIPT_DIR}/select-gpu-partition.sh"
    BENCH_GPU_PARTITION="$SELECTED_GPU_PARTITION"
    export BENCH_GPU_PARTITION
fi

# The bench script must come from the PR tree (master may not have it).
PR_BENCH_SCRIPT="$(cd "${SCRIPT_DIR}/../workflows/common" && pwd)/bench.sh"
# Must match the slug computed by submit-slurm-job.sh:
#   basename("bench.sh") → "bench" → "bench-${device}-${interface}"
job_slug="bench-${device}-${interface}"

# --- Phase 1: Submit both SLURM jobs (no monitoring yet) ---
echo "Submitting PR benchmark..."
(cd pr && SUBMIT_ONLY=1 bash "${SCRIPT_DIR}/submit-slurm-job.sh" "$PR_BENCH_SCRIPT" "$device" "$interface" "$cluster")
pr_job_id=$(cat "pr/${job_slug}.slurm_job_id")
echo "PR job submitted: $pr_job_id"

echo "Submitting master benchmark..."
(cd master && SUBMIT_ONLY=1 bash "${SCRIPT_DIR}/submit-slurm-job.sh" "$PR_BENCH_SCRIPT" "$device" "$interface" "$cluster")
master_job_id=$(cat "master/${job_slug}.slurm_job_id")
echo "Master job submitted: $master_job_id"

echo "Both SLURM jobs submitted — running concurrently on compute nodes."
echo "Monitoring sequentially to conserve login node memory."

# --- Phase 2: Monitor sequentially (one at a time on login node) ---
# On Phoenix 'embers' a long benchmark job can be preempted (PreemptMode=CANCEL,
# so it is killed rather than requeued). On preemption (run_monitored exit 76)
# resubmit a fresh job in the same tree and re-monitor, bounded by
# MAX_PREEMPT_RESUBMITS (the 480m job timeout is the real backstop). Note: a
# resubmitted job no longer overlaps its counterpart, slightly reducing
# same-load fairness -- still preferable to failing the run on an infra preempt.
: "${MAX_PREEMPT_RESUBMITS:=10}"
monitor_bench_with_resubmit() {  # arg: <dir> (pr|master); sets BENCH_MON_RC
    local dir="$1"
    local out="${dir}/${job_slug}.out"
    local jobid attempt=0 rc
    jobid=$(cat "${dir}/${job_slug}.slurm_job_id")
    while :; do
        rc=0
        bash "${SCRIPT_DIR}/run_monitored_slurm_job.sh" "$jobid" "$out" || rc=$?
        if [ "$rc" -ne 76 ]; then
            BENCH_MON_RC="$rc"
            return
        fi
        if [ "$attempt" -ge "$MAX_PREEMPT_RESUBMITS" ]; then
            echo "::error::${dir} benchmark preempted ${MAX_PREEMPT_RESUBMITS}x without completing; giving up."
            BENCH_MON_RC=1
            return
        fi
        attempt=$((attempt + 1))
        echo "::warning::${dir} benchmark job $jobid was preempted; resubmitting (attempt ${attempt}/${MAX_PREEMPT_RESUBMITS})."
        rm -f "$out"
        ( cd "$dir" && SUBMIT_ONLY=1 bash "${SCRIPT_DIR}/submit-slurm-job.sh" "$PR_BENCH_SCRIPT" "$device" "$interface" "$cluster" )
        jobid=$(cat "${dir}/${job_slug}.slurm_job_id")
        echo "${dir} benchmark resubmitted as job $jobid"
    done
}

echo ""
echo "=== Monitoring PR job $pr_job_id ==="
monitor_bench_with_resubmit pr
pr_exit=$BENCH_MON_RC
if [ "$pr_exit" -ne 0 ]; then
    echo "PR job exited with code: $pr_exit"
    tail -n 50 "pr/${job_slug}.out" 2>/dev/null || echo "  Could not read PR log"
    # The PR benchmark run genuinely failed (cases crashed/hung/SIGTERM'd, not a
    # monitor false-positive -- run_monitored_slurm_job.sh re-checks sacct). Fail
    # the job instead of falling through to the YAML-exists check, which would let
    # a broken PR pass green as long as a partial YAML was written. Scoped to PR
    # only: a master/baseline infra flake stays a warning and does not red-cross.
    exit 1
else
    echo "PR job completed successfully"
fi

echo ""
echo "=== Monitoring master job $master_job_id ==="
monitor_bench_with_resubmit master
master_exit=$BENCH_MON_RC
if [ "$master_exit" -ne 0 ]; then
    echo "Master job exited with code: $master_exit"
    tail -n 50 "master/${job_slug}.out" 2>/dev/null || echo "  Could not read master log"
else
    echo "Master job completed successfully"
fi

# --- Phase 3: Verify outputs ---
if [ "${pr_exit}" -ne 0 ] || [ "${master_exit}" -ne 0 ]; then
    echo "WARNING: Benchmark jobs had failures: pr=${pr_exit}, master=${master_exit}"
    echo "Checking for partial results..."
else
    echo "===================="
    echo "Both benchmark jobs completed successfully!"
    echo "===================="
fi

pr_yaml="pr/${job_slug}.yaml"
master_yaml="master/${job_slug}.yaml"

# Wait briefly for YAML files to appear on NFS.  When monitoring starts
# after a job has already completed (common for the second job), the
# recovery path in run_monitored_slurm_job.sh sleeps 30s, but NFS
# propagation can take longer under load.
for yaml in "$pr_yaml" "$master_yaml"; do
    attempts=0
    while [ ! -f "$yaml" ] && [ $attempts -lt 6 ]; do
        echo "Waiting for $yaml to appear (NFS propagation)..."
        sleep 5
        attempts=$((attempts + 1))
    done
done

if [ ! -f "$pr_yaml" ]; then
    echo "ERROR: PR benchmark output not found: $pr_yaml"
    ls -la pr/ || true
    echo ""
    tail -n 100 "pr/${job_slug}.out" 2>/dev/null || echo "  Could not read PR log"
    exit 1
fi

if [ ! -f "$master_yaml" ]; then
    echo "ERROR: Master benchmark output not found: $master_yaml"
    ls -la master/ || true
    echo ""
    tail -n 100 "master/${job_slug}.out" 2>/dev/null || echo "  Could not read master log"
    exit 1
fi

echo "Verified both YAML files exist:"
echo "  - $pr_yaml"
echo "  - $master_yaml"
