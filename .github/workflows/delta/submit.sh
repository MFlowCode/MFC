#!/bin/bash

set -e

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
#SBATCH -p cpu               # partition
#SBATCH --ntasks-per-node=24       # Number of cores per node required
#SBATCH --account=bdiy-delta-cpu # charge account
"

sbatch_gpu_opts="\
#SBATCH --account=bdiy-delta-gpu
#SBATCH -p gpu
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

job_slug="`basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g'`-$2-$3"

sbatch <<EOT
#!/bin/bash
#SBATCH -Jshb-$job_slug            # Job name
#SBATCH -N1                        # Number of nodes required
$sbatch_device_opts
#SBATCH -t 01:00:00                # Duration of the job (Ex: 15 mins)
#SBATCH -o$job_slug.out            # Combined output and error messages file
#SBATCH -W                         # Do not exit until the submitted job terminates.
#SBATCH --exclusive                         # Do not exit until the submitted job terminates.

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in $(pwd):"

echo "Job started on nodes: $SLURM_JOB_NODELIST"
echo "Number of tasks requested: $SLURM_NTASKS"
echo "CPUs per task requested: $SLURM_CPUS_PER_TASK"
echo "Total allocated CPU cores: $(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))"
echo "Cores available on this specific node: $SLURM_CPUS_ON_NODE"

job_slug="$job_slug"
job_device="$2"
job_interface="$3"

. ./mfc.sh load -c d -m $2

$sbatch_script_contents

EOT
