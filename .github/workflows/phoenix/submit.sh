#!/bin/bash

set -e

# Ignore SIGHUP to survive login node session drops
trap '' HUP

usage() {
    echo "Usage: $0 [script.sh] [cpu|gpu] [none|acc|omp]"
}

if [ ! -z "$1" ]; then
    sbatch_script_contents=`cat $1`
else
    usage
    exit 1
fi

# Detect job type from submitted script basename
script_basename="$(basename "$1" .sh)"
case "$script_basename" in
    bench*) job_type="bench" ;;
    *)      job_type="test"  ;;
esac

sbatch_cpu_opts="\
#SBATCH -p cpu-gnr                 # partition (full Granite Rapids node)
#SBATCH --exclusive                # exclusive access to all cores
#SBATCH -C graniterapids           # constrain to GNR architecture\
"

if [ "$job_type" = "bench" ]; then
    sbatch_gpu_opts="\
#SBATCH -CL40S
#SBATCH --ntasks-per-node=4       # Number of MPI tasks per node required
#SBATCH -G2\
"
    sbatch_time="#SBATCH -t 04:00:00"
else
    sbatch_gpu_opts="\
#SBATCH -p gpu-v100,gpu-a100,gpu-h100,gpu-l40s
#SBATCH --ntasks-per-node=4       # Number of MPI tasks per node required
#SBATCH -G2\
"
    sbatch_time="#SBATCH -t 03:00:00"
fi

if [ "$2" = "cpu" ]; then
    sbatch_device_opts="$sbatch_cpu_opts"
elif [ "$2" = "gpu" ]; then
    sbatch_device_opts="$sbatch_gpu_opts"
else
    usage
    exit 1
fi

job_slug="`basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g'`-$2-$3"
output_file="$job_slug.out"

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
echo "Running in $(pwd):"

job_slug="$job_slug"
job_device="$2"
job_interface="$3"
export GITHUB_EVENT_NAME="$GITHUB_EVENT_NAME"

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

# Use resilient monitoring instead of sbatch -W
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
monitor_exit=0
bash "$SCRIPT_DIR/../../scripts/monitor_slurm_job.sh" "$job_id" "$output_file" || monitor_exit=$?

if [ "$monitor_exit" -ne 0 ]; then
  echo "Monitor exited with code $monitor_exit; re-checking SLURM job $job_id final state..."
  # Give the SLURM epilog time to finalize if the job just finished
  sleep 30
  final_state=$(sacct -j "$job_id" -n -X -P -o State 2>/dev/null | head -n1 | cut -d'|' -f1 | tr -d ' ' || echo "UNKNOWN")
  final_exit=$(sacct -j "$job_id" --format=ExitCode --noheader --parsable2 2>/dev/null | head -n1 | tr -d ' ' || echo "")
  echo "Final SLURM state=$final_state exit=$final_exit"
  if [ "$final_state" = "COMPLETED" ] && [ "$final_exit" = "0:0" ]; then
    echo "SLURM job $job_id completed successfully despite monitor failure — continuing."
  else
    echo "ERROR: SLURM job $job_id did not complete successfully (state=$final_state exit=$final_exit)"
    exit 1
  fi
fi
