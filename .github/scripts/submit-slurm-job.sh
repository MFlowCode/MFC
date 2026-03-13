#!/bin/bash
# Unified SLURM job submission and monitoring for all clusters.
# Submits a script as a SLURM batch job, then monitors it until completion.
# Rerun-safe: cancels stale jobs from previous runs before resubmission.
#
# Usage: submit-slurm-job.sh <script.sh> <cpu|gpu> <none|acc|omp> <cluster> [shard]

set -euo pipefail

# Ignore SIGHUP to survive login node session drops
trap '' HUP

usage() {
    echo "Usage: $0 <script.sh> <cpu|gpu> <none|acc|omp> <cluster> [shard]"
}

script_path="${1:-}"
device="${2:-}"
interface="${3:-}"
cluster="${4:-}"
shard="${5:-}"

if [ -z "$script_path" ] || [ -z "$device" ] || [ -z "$interface" ] || [ -z "$cluster" ]; then
    usage
    exit 1
fi

sbatch_script_contents=$(cat "$script_path")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Detect job type from submitted script basename
script_basename="$(basename "$script_path" .sh)"
case "$script_basename" in
    bench*) job_type="bench" ;;
    *)      job_type="test"  ;;
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
        bench_time="04:00:00"
        gpu_partition_dynamic=true
        ;;
    frontier)
        compiler_flag="f"
        account="CFD154"
        job_prefix="MFC"
        qos="develop"
        extra_sbatch=""
        test_time="01:59:00"
        bench_time="01:59:00"
        gpu_partition_dynamic=false
        ;;
    frontier_amd)
        compiler_flag="famd"
        account="CFD154"
        job_prefix="MFC"
        qos="develop"
        extra_sbatch=""
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
if [ "$job_type" = "bench" ]; then
    sbatch_time="#SBATCH -t $bench_time"
else
    sbatch_time="#SBATCH -t $test_time"
fi

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
            sbatch_device_opts="\
#SBATCH -n 32
#SBATCH -p service"
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
            sbatch_device_opts="\
#SBATCH -p $gpu_partition
#SBATCH --ntasks-per-node=4
#SBATCH -G2
#SBATCH --exclude=atl1-1-03-002-29-0"
            ;;
        frontier|frontier_amd)
            sbatch_device_opts="\
#SBATCH -n 8
#SBATCH -p service"
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
job_slug="$(basename "$script_path" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g')-${device}-${interface}${shard_suffix}"
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
_sbatch_script=$(cat <<EOT
#!/bin/bash
#SBATCH -J ${job_prefix}-${job_slug}
#SBATCH --account=${account}
#SBATCH -N 1
${sbatch_device_opts}
${sbatch_time}
#SBATCH --qos=${qos}
${extra_sbatch}
#SBATCH -o ${output_file}

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in \$(pwd):"

job_slug="$job_slug"
job_device="$device"
job_interface="$interface"
job_shard="$shard"
job_cluster="$cluster"

. ./mfc.sh load -c $compiler_flag -m $module_mode

$sbatch_script_contents

EOT
)

job_id=$(retry_sbatch "$_sbatch_script")
unset _sbatch_script

echo "Submitted batch job $job_id"
echo "$job_id" > "$id_file"
echo "Job ID written to $id_file"

# --- Monitor ---
bash "$SCRIPT_DIR/run_monitored_slurm_job.sh" "$job_id" "$output_file"
