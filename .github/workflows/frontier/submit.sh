#!/bin/bash

set -e

# Ignore SIGHUP to survive login node session drops
trap '' HUP

usage() {
    echo "Usage: $0 [script.sh] [cpu|gpu]"
}

if [ ! -z "$1" ]; then
    sbatch_script_contents=`cat $1`
else
    usage
    exit 1
fi

if [ "$2" = "cpu" ]; then
    sbatch_device_opts="\
#SBATCH -n 32                       # Number of cores required"
elif [ "$2" = "gpu" ]; then
    sbatch_device_opts="\
#SBATCH -n 8                       # Number of cores required"
else
    usage
    exit 1
fi


job_slug="`basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g'`-$2-$3"
output_file="$job_slug.out"

submit_output=$(sbatch <<EOT
#!/bin/bash
#SBATCH -J MFC-$job_slug            # Job name
#SBATCH -A ENG160                  # charge account
#SBATCH -N 1                       # Number of nodes required
$sbatch_device_opts
#SBATCH -t 05:59:00                # Duration of the job (Ex: 15 mins)
#SBATCH -o$output_file             # Combined output and error messages file
#SBATCH -p extended                # Extended partition for shorter queues

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in $(pwd):"

job_slug="$job_slug"
job_device="$2"
job_interface="$3"

. ./mfc.sh load -c f -m $([ "$2" = "gpu" ] && echo "g" || echo "c")

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
