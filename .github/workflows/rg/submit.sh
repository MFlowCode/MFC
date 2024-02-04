#!/bin/bash

set -e

usage() {
    echo "Usage: $0 [script.sh] [gpu]"
}

if [ ! -z "$1" ]; then
    sbatch_script_contents=`cat $1`
else
    usage
    exit 1
fi

sbatch_gpu_opts="\
#SBATCH -G 1
#SBATCH --nodelist violet1 #node\
"

if [ "$2" == "gpu" ]; then
    sbatch_device_opts="$sbatch_gpu_opts"
else
    usage
    exit 1
fi

job_slug="`basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g'`-$2"

sbatch <<EOT
#!/bin/bash
#SBATCH -JMFC-rg-$job_slug              # Job name
#SBATCH -N1 --ntasks-per-node=4         # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=4G                # Memory per core
#SBATCH -t 02:00:00                     # Duration of the job (Ex: 30 mins)
#SBATCH -p rg-nextgen-hpc               # Partition Name
#SBATCH -o$jobslug.out                  # Combined output
$sbatch_device_opts
#SBATCH -W                  # Do not exit until the submitted job terminates.

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in $(pwd):"

source /tools/misc/.read_profile

job_slug="$job_slug"
job_device="$2"

module load nvhpc cmake

$sbatch_script_contents

EOT

