#!/bin/bash
#SBATCH -Jshb-test-jobs            # Job name
#SBATCH --account=gts-sbryngelson3 # charge account
#SBATCH -N1                        # Number of nodes and cores per node required
#SBATCH -CV100-16GB
#SBATCH -G2
#SBATCH -t 02:00:00                # Duration of the job (Ex: 15 mins)
#SBATCH -q embers                  # QOS Name
#SBATCH -otest.out                 # Combined output and error messages file
#SBATCH -W                         # Do not exit until the submitted job terminates.

cd "$SLURM_SUBMIT_DIR"
echo "Running in $(pwd):"

set -x

. ./mfc.sh load -c p -m GPU

gpu_count=$(nvidia-smi -L | wc -l)        # number of GPUs on node
gpu_ids=$(seq -s ' ' 0 $(($gpu_count-1))) # 0,1,2,...,gpu_count-1

./mfc.sh test -a -j 2 --gpu -g $gpu_ids -- -c phoenix
