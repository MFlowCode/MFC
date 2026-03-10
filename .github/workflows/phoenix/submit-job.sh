#!/bin/bash
# Submit a SLURM job without waiting for it to complete.
# Writes the job ID to <job_slug>.slurm_job_id so a separate monitor step can wait.
# Idempotent: if a job for this slug is still RUNNING or PENDING, skip resubmission.
#
# Usage: submit-job.sh [script.sh] [cpu|gpu] [none|acc|omp]

set -euo pipefail

# Ignore SIGHUP to survive login node session drops
trap '' HUP

usage() {
    echo "Usage: $0 [script.sh] [cpu|gpu] [none|acc|omp]"
}

if [ -z "${1:-}" ]; then
    usage
    exit 1
fi

sbatch_script_contents=$(cat "$1")

# Detect job type from submitted script basename
script_basename="$(basename "$1" .sh)"
case "$script_basename" in
    bench*) job_type="bench" ;;
    *)      job_type="test"  ;;
esac

sbatch_cpu_opts="\
#SBATCH -p cpu-small               # partition
#SBATCH --ntasks-per-node=24       # Number of cores per node required
#SBATCH --mem-per-cpu=2G           # Memory per core\
"

if [ "$job_type" = "bench" ]; then
    # BENCH_GPU_PARTITION is pre-selected by run_parallel_benchmarks.sh so both
    # PR and master jobs land on the same GPU type for a fair comparison.
    gpu_partition="${BENCH_GPU_PARTITION:-gpu-rtx6000}"
    echo "Submitting bench GPU job to partition: $gpu_partition (BENCH_GPU_PARTITION=${BENCH_GPU_PARTITION:-<unset, using default>})"
    sbatch_time="#SBATCH -t 04:00:00"
else
    source "$(dirname "${BASH_SOURCE[0]}")/../../scripts/select-gpu-partition.sh"
    gpu_partition="$SELECTED_GPU_PARTITION"
    sbatch_time="#SBATCH -t 03:00:00"
fi

sbatch_gpu_opts="\
#SBATCH -p $gpu_partition
#SBATCH --ntasks-per-node=4       # Number of cores per node required
#SBATCH -G2\
"

if [ "$2" = "cpu" ]; then
    sbatch_device_opts="$sbatch_cpu_opts"
elif [ "$2" = "gpu" ]; then
    sbatch_device_opts="$sbatch_gpu_opts"
else
    usage
    exit 1
fi

job_slug="$(basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g')-$2-$3"
output_file="$job_slug.out"
id_file="${job_slug}.slurm_job_id"

# On rerun, cancel any existing job for this slug and submit a fresh one.
# If the job is still live (RUNNING/PENDING), scancel it first as a safety net
# in case the "Cancel SLURM Jobs" step did not fire (e.g. runner was SIGKILL'd).
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

submit_output=$(sbatch <<EOT
#!/bin/bash
#SBATCH -Jshb-$job_slug            # Job name
#SBATCH --account=gts-sbryngelson3 # charge account
#SBATCH -N1                        # Number of nodes required
$sbatch_device_opts
$sbatch_time
#SBATCH -q embers                  # QOS Name
#SBATCH --requeue                  # Auto-requeue on preemption
#SBATCH -o$output_file             # Combined output and error messages file

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in \$(pwd):"

job_slug="$job_slug"
job_device="$2"
job_interface="$3"

. ./mfc.sh load -c p -m $2

$sbatch_script_contents

EOT
)

job_id=$(echo "$submit_output" | grep -oE '[0-9]+')
if [ -z "$job_id" ]; then
    echo "ERROR: Failed to submit job. sbatch output:"
    echo "$submit_output"
    exit 1
fi

echo "Submitted batch job $job_id"
echo "$job_id" > "$id_file"
echo "Job ID written to $id_file"
