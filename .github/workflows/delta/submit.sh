#!/bin/bash

set -e

usage() {
    echo "Usage: $0 [script.sh] [cpu|gpu]"
}

if [ ! -z "$1" ]; then
    sbatch_script_contents=`cat $1`
else
    usage
    exit 1
fi

sbatch_cpu_opts="\
#SBATCH -p cpu
#SBATCH --account=bdiy-delta-cpu
"

sbatch_gpu_opts="\
#SBATCH -p gpuA100x4,gpuA100x4-interactive
#SBATCH --account=bdiy-delta-gpu
#SBATCH --gpus-per-node=4
"

if [ "$2" == "cpu" ]; then
    sbatch_device_opts="$sbatch_cpu_opts"
elif [ "$2" == "gpu" ]; then
    sbatch_device_opts="$sbatch_gpu_opts"
else
    usage
    exit 1
fi

job_slug="`basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g'`-$2"

sbatch <<EOT
#!/bin/bash
#SBATCH -Jshb-$job_slug            # Job name
#SBATCH -N1                        # Number of nodes required
$sbatch_device_opts
#SBATCH -t 01:00:00                # Duration of the job (Ex: 15 mins)
#SBATCH -n 20
#SBATCH -o$job_slug.out            # Combined output and error messages file
#SBATCH --constraint="scratch"
#SBATCH -W                         # Do not exit until the submitted job terminates.

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in $(pwd):"

job_slug="$job_slug"
job_device="$2"

. ./mfc.sh load -c d -m $2

$sbatch_script_contents

EOT

