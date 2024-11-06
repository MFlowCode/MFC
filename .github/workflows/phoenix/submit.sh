#!/bin/bash

set -e

usage() {
    echo "Usage: $0 [script.sh] [cpu|gpu] [--single (optional)]"
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
#SBATCH -C V100-16GB
#SBATCH -G2\
"

if [ "$2" == "cpu" ]; then
    sbatch_device_opts="$sbatch_cpu_opts"
elif [ "$2" == "gpu" ]; then
    sbatch_device_opts="$sbatch_gpu_opts"
else
    usage
    exit 1
fi

# Check for the --single flag
single_flag=""
if [ "$3" == "--single" ]; then
    single_flag="--single"
fi

job_slug="`basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g'`-$2"

# **Add this block to adjust job_slug if --single is used**
if [ "$single_flag" == "--single" ]; then
    job_slug="${job_slug}-single"
fi

sbatch <<EOT
#!/bin/bash
#SBATCH -Jshb-$job_slug            # Job name
#SBATCH --account=gts-sbryngelson3 # Charge account
#SBATCH -N1                        # Number of nodes required
$sbatch_device_opts
#SBATCH -t 04:00:00                # Duration of the job
#SBATCH -q embers                  # QOS Name
#SBATCH -o $job_slug.out           # Output file
#SBATCH -W                         # Wait for the job to finish

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in \$(pwd):"

job_slug="$job_slug"
job_device="$2"
single_flag="$single_flag"

. ./mfc.sh load -c p -m $2

$sbatch_script_contents

EOT
