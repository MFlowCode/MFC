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

./mfc.sh test -a -b mpirun -j $(nproc) \
              --gpu -g $(seq -s ',' 0 $(($(nvidia-smi -L | wc -l)-1)))

