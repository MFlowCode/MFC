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

sbatch_cpu_opts="\
#SBATCH -p cpu-small               # partition
#SBATCH --ntasks-per-node=24       # Number of cores per node required
#SBATCH --mem-per-cpu=2G           # Memory per core\
"

sbatch_gpu_opts="\
#SBATCH --gres=gpu:H200:2
#SBATCH --ntasks-per-node=8       # Number of tasks (MPI ranks) per node\
"

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
#SBATCH -t 03:00:00                # Duration of the job (Ex: 15 mins)
#SBATCH -q embers                  # QOS Name
#SBATCH -o$output_file             # Combined output and error messages file

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in $(pwd):"

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

# Use resilient monitoring instead of sbatch -W
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "$SCRIPT_DIR/../../scripts/monitor_slurm_job.sh" "$job_id" "$output_file"
