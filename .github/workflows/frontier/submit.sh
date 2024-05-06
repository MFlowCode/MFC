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

job_slug="`basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g'`-$2"

sbatch <<EOT
#!/bin/bash
#SBATCH -JMFC-$job_slug            # Job name
#SBATCH -A CFD154                  # charge account
#SBATCH -N 1                       # Number of nodes required
#SBATCH -n 8                       # Number of cores required
#SBATCH -t 02:00:00                # Duration of the job (Ex: 15 mins)
#SBATCH -q debug                   # QOS Name
#SBATCH -o$job_slug.out            # Combined output and error messages file
#SBATCH -W                         # Do not exit until the submitted job terminates.

set -e
set -x

cd "\$SLURM_SUBMIT_DIR"
echo "Running in $(pwd):"

job_slug="$job_slug"
job_device="$2"

. ./mfc.sh load -c f -m g

$sbatch_script_contents

EOT

