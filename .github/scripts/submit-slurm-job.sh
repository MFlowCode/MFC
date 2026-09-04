#!/bin/bash
# Unified SLURM job submission and monitoring for all clusters.
# Submits a script as a SLURM batch job, then monitors it until completion.
# Rerun-safe: cancels stale jobs from previous runs before resubmission.
#
# Usage: submit-slurm-job.sh <script.sh> <cpu|gpu> <none|acc|omp> <cluster> [shard] [variant]

set -euo pipefail

# Ignore SIGHUP to survive login node session drops
trap '' HUP

usage() {
    echo "Usage: $0 <script.sh> <cpu|gpu> <none|acc|omp> <cluster> [shard] [variant]"
}

script_path="${1:-}"
device="${2:-}"
interface="${3:-}"
cluster="${4:-}"
shard="${5:-}"
# Build variant ("base"/"chem"), letting one build job be split into concurrent
# per-variant jobs. Empty = build everything in this job (default).
variant="${6:-}"

if [ -z "$script_path" ] || [ -z "$device" ] || [ -z "$interface" ] || [ -z "$cluster" ]; then
    usage
    exit 1
fi

sbatch_script_contents=$(cat "$script_path")

# Nodes this job must not be scheduled onto. Seeded below per cluster with hosts
# already known to be bad, and grown at runtime when the in-allocation preflight
# reports a node fault (exit 77). Growing it automatically is the point: the two
# Phoenix nodes in the seed list were each found by hand, diagnosed, and
# committed, after one of them alone had eaten 25 jobs.
node_exclude=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Detect job type from submitted script basename
script_basename="$(basename "$script_path" .sh)"
case "$script_basename" in
    bench*)          job_type="bench" ;;
    build-and-test*) job_type="buildtest" ;;
    *)               job_type="test"  ;;
esac

# --- Cluster configuration ---
case "$cluster" in
    phoenix)
        compiler_flag="p"
        account="gts-sbryngelson3"
        job_prefix="shb"
        qos="embers"
        extra_sbatch="#SBATCH --requeue"
        test_time="03:00:00"
        # Combined build+test needs build headroom on top of the test budget;
        # kept modest to still backfill under 'embers'.
        buildtest_time="03:30:00"
        bench_time="04:00:00"
        gpu_partition_dynamic=true
        ;;
    frontier)
        compiler_flag="f"
        account="CFD154"
        job_prefix="MFC"
        # The hackathon QOS was a temporary grant and is no longer held by
        # CFD154; submitting under it now fails outright with "Invalid qos
        # specification". "normal" is the only QOS on this allocation without a
        # one-job-at-a-time cap, so it is the only one that can run the CI
        # matrix concurrently. Note that the g1 partition carries its own
        # partition QOS (also named "g1"), which slurmctld applies on its own
        # when a job lands there. Do not add --qos=g1: CFD154 has no
        # association with it and sbatch rejects the job outright.
        qos="normal"
        # Let each job's slurmstepd broker its own steps instead of routing
        # every srun through slurmctld. The in-job test suite launches ~1700+
        # srun steps per allocation, which congests the Frontier controller.
        extra_sbatch="#SBATCH --stepmgr"
        test_time="01:59:00"
        bench_time="01:59:00"
        gpu_partition_dynamic=false
        ;;
    frontier_amd)
        compiler_flag="famd"
        account="CFD154"
        job_prefix="MFC"
        qos="normal"
        extra_sbatch="#SBATCH --stepmgr"
        test_time="01:59:00"
        bench_time="01:59:00"
        gpu_partition_dynamic=false
        ;;
    *)
        echo "ERROR: Unknown cluster '$cluster'"
        exit 1
        ;;
esac

# --- Time limit ---
case "$job_type" in
    bench)     sbatch_time="#SBATCH -t $bench_time" ;;
    buildtest) sbatch_time="#SBATCH -t ${buildtest_time:-$test_time}" ;;
    *)         sbatch_time="#SBATCH -t $test_time" ;;
esac

# --- Device-specific SBATCH options ---
if [ "$device" = "cpu" ]; then
    case "$cluster" in
        phoenix)
            sbatch_device_opts="\
#SBATCH -p cpu-small,cpu-medium,cpu-large
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=8G"
            ;;
        frontier|frontier_amd)
            # g1 is a dedicated 64-node carve-out; its nodes are not in batch,
            # so CI starts promptly instead of queueing behind the machine.
            sbatch_device_opts="\
#SBATCH -n 32
#SBATCH -p g1"
            ;;
    esac
elif [ "$device" = "gpu" ]; then
    # Determine GPU partition
    gpu_partition="batch"
    if [ "$gpu_partition_dynamic" = "true" ]; then
        # Use pre-selected bench partition if available, otherwise query sinfo
        if [ -n "${BENCH_GPU_PARTITION:-}" ]; then
            gpu_partition="$BENCH_GPU_PARTITION"
            echo "Using pre-selected bench partition: $gpu_partition (PR/master consistency)"
        else
            source "${SCRIPT_DIR}/select-gpu-partition.sh"
            gpu_partition="$SELECTED_GPU_PARTITION"
        fi
    fi

    case "$cluster" in
        phoenix)
            # --exclude is rendered separately (see $node_exclude) so the
            # preflight can add a node to it and resubmit.
            sbatch_device_opts="\
#SBATCH -p $gpu_partition
#SBATCH --ntasks-per-node=4
#SBATCH -G2"
            node_exclude="atl1-1-03-007-29-0,atl1-1-03-007-31-0"
            ;;
        frontier|frontier_amd)
            sbatch_device_opts="\
#SBATCH -n 8
#SBATCH -p g1"
            # Seed, same as phoenix above: the preflight adds nodes to this at
            # run time. frontier10202 produced all 183 GPU memory-access faults
            # in run 33553417354 (43 distinct tests) while the same lanes passed
            # on eight other g1 nodes with none. Its faults are intermittent --
            # 379 of 382 tests still passed there -- so syscheck can clear it.
            node_exclude="frontier10202"
            ;;
    esac
else
    usage
    exit 1
fi

# --- Job slug ---
shard_suffix=""
if [ -n "$shard" ]; then
    shard_suffix="-$(echo "$shard" | sed 's|/|-of-|')"
fi
variant_suffix=""
if [ -n "$variant" ]; then
    # Sanitized like the rest of the slug: these become file names in the workspace.
    variant_suffix="-$(echo "$variant" | sed 's/[^a-zA-Z0-9]/-/g')"
fi
job_slug="$(basename "$script_path" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g')-${device}-${interface}${shard_suffix}${variant_suffix}"
output_file="$job_slug.out"
id_file="${job_slug}.slurm_job_id"

# --- Idempotency: cancel stale jobs from previous runs ---
if [ -f "$id_file" ]; then
    existing_id=$(cat "$id_file")
    state=$(sacct -j "$existing_id" -n -X -P -o State 2>/dev/null | head -n1 | cut -d'|' -f1 | tr -d ' ' || true)
    case "${state:-UNKNOWN}" in
        RUNNING|PENDING|REQUEUED|COMPLETING)
            echo "Cancelling stale SLURM job $existing_id (state=$state) before resubmission"
            scancel "$existing_id" 2>/dev/null || true
            ;;
        *)
            echo "Stale job $existing_id (state=${state:-UNKNOWN}) — submitting fresh"
            ;;
    esac
    rm -f "$id_file"
fi

# Remove stale output file so the monitor doesn't pick up old content
# (a previous SLURM job's epilog can write to the .out file after our
# stale-job check, polluting the new job's output stream).
rm -f "$output_file"

# --- Module load mode (short form) ---
module_mode=$([ "$device" = "gpu" ] && echo "g" || echo "c")

# --- Submit (with retries for transient SLURM errors) ---
source "${SCRIPT_DIR}/retry-sbatch.sh"
# Re-rendered before every submission so a node added to $node_exclude by a
# failed preflight actually takes effect on the retry.
render_sbatch_script() {
    local exclude_directive=""
    # An `if`, not `[ ... ] && ...`: under `set -e` the latter returns non-zero
    # whenever the list is empty and would abort the script.
    if [ -n "$node_exclude" ]; then
        exclude_directive="#SBATCH --exclude=${node_exclude}"
    fi
    cat <<EOT
#!/bin/bash
#SBATCH -J ${job_prefix}-${job_slug}
#SBATCH --account=${account}
#SBATCH -N 1
${sbatch_device_opts}
${exclude_directive}
${sbatch_time}
#SBATCH --qos=${qos}
${extra_sbatch}
#SBATCH -o ${output_file}

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in \$(pwd):"

# Exported so wrapper scripts (build-and-test.sh) run child scripts that inherit these.
export job_slug="$job_slug"
export job_device="$device"
export job_interface="$interface"
export job_shard="$shard"
export job_variant="$variant"
export job_cluster="$cluster"
export GITHUB_EVENT_NAME="$GITHUB_EVENT_NAME"

. ./mfc.sh load -c $compiler_flag -m $module_mode

$sbatch_script_contents

EOT
}

# --- Submit + monitor, resubmitting on preemption
# Phoenix preempts 'embers' jobs with PreemptMode=CANCEL (not REQUEUE), so a
# preempted job is killed outright and `--requeue` never restarts it. When the
# monitor reports preemption (exit 76), submit a fresh job and monitor again.
# Bounded by MAX_PREEMPT_RESUBMITS as a runaway guard; the job-level
# `timeout-minutes` (480m) remains the real backstop.
: "${MAX_PREEMPT_RESUBMITS:=10}"
# Node faults get a much tighter bound than preemption: preemption is routine on
# 'embers' and says nothing about the node, whereas hitting a second unusable
# node in a row means the problem is the cluster, not the draw.
#
# One, not two. Each attempt costs a node if the probe is wrong, and it has been
# wrong: a bounded loop faithfully condemned three healthy Phoenix nodes when the
# probe was given a binary that could not run there. Bad nodes are concentrated
# -- one accounted for 25 of 29 ECC failures -- so a single requeue captures
# nearly all of the benefit at half the blast radius.
: "${MFC_MAX_NODE_RESUBMITS:=1}"
preempt_attempt=0
node_attempt=0
while :; do
    _sbatch_script=$(render_sbatch_script)
    job_id=$(retry_sbatch "$_sbatch_script")
    echo "Submitted batch job $job_id"
    echo "$job_id" > "$id_file"
    echo "Job ID written to $id_file"

    # SUBMIT_ONLY=1 (parallel submission, e.g. benchmarks): the caller monitors
    # each job and handles its own preemption resubmits.
    if [ "${SUBMIT_ONLY:-0}" = "1" ]; then
        echo "SUBMIT_ONLY mode: skipping monitor (job_id=$job_id output=$output_file)"
        break
    fi

    monitor_rc=0
    bash "$SCRIPT_DIR/run_monitored_slurm_job.sh" "$job_id" "$output_file" || monitor_rc=$?
    if [ "$monitor_rc" -eq 0 ]; then
        break
    fi
    if [ "$monitor_rc" -eq 76 ]; then
        if [ "$preempt_attempt" -lt "$MAX_PREEMPT_RESUBMITS" ]; then
            preempt_attempt=$((preempt_attempt + 1))
            echo "::warning::SLURM job $job_id was preempted (Phoenix embers). Resubmitting a fresh job (attempt ${preempt_attempt}/${MAX_PREEMPT_RESUBMITS})."
            rm -f "$output_file"
            continue
        fi
        echo "::error::SLURM job preempted ${MAX_PREEMPT_RESUBMITS} times without completing; giving up."
        exit 1
    fi
    if [ "$monitor_rc" -eq 77 ]; then
        # The in-allocation preflight found this node unusable before any real
        # work started. Exclude it and draw another node.
        faulted_node=$(bash "$SCRIPT_DIR/node-exclude.sh" node-from "$output_file")
        if [ "$node_attempt" -lt "$MFC_MAX_NODE_RESUBMITS" ]; then
            node_attempt=$((node_attempt + 1))
            node_exclude=$(bash "$SCRIPT_DIR/node-exclude.sh" merge "$node_exclude" "$faulted_node")
            echo "::warning::SLURM job $job_id failed preflight on ${faulted_node:-an unidentified node}. Excluding it and resubmitting (attempt ${node_attempt}/${MFC_MAX_NODE_RESUBMITS}). Excluding: ${node_exclude}"
            rm -f "$output_file"
            continue
        fi
        echo "::error::Preflight failed on $((MFC_MAX_NODE_RESUBMITS + 1)) nodes in a row (excluded: ${node_exclude})."
        echo "That is a cluster-wide problem rather than a bad draw; not resubmitting."
        exit 1
    fi
    # Genuine failure (not preemption or infrastructure).
    exit "$monitor_rc"
done
unset _sbatch_script
