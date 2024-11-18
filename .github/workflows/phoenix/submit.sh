#!/bin/bash

set -e
usage() {
    echo "Usage: $0 [script.sh] [cpu|gpu] [single|double]"
}

if [ ! -z "$1" ]; then
    script_path="$1"
else
    usage
    exit 1
fi

if [ "$2" == "cpu" ]; then
    sbatch_device_opts="#SBATCH -p cpu-small
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=2G"
elif [ "$2" == "gpu" ]; then
    sbatch_device_opts="#SBATCH -C V100-16GB
#SBATCH -G2"
else
    usage
    exit 1
fi

if [ "$3" != "single" ] && [ "$3" != "double" ]; then
    usage
    exit 1
fi

job_device="$2"
job_precision="$3"
job_slug="$(basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g')-$job_device-$job_precision"

sbatch <<EOT
#!/bin/bash
#SBATCH -Jshb-$job_slug            # Job name
#SBATCH --account=gts-sbryngelson3 # Charge account
#SBATCH -N1                        # Number of nodes required
$sbatch_device_opts
#SBATCH -t 02:00:00                # Duration of the job
#SBATCH -q embers                  # QOS Name
#SBATCH -o $job_slug-$job_precision.out  # Include precision in the log file
#SBATCH -W                         # Do not exit until the submitted job terminates.

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in \$(pwd):"

# Load necessary modules
. ./mfc.sh load -c p -m \$job_device

# Execute the script with arguments
bash $script_path "\$job_device" "\$job_precision"
EOT

