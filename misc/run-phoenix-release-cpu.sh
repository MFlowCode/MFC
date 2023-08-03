#!/bin/bash
#SBATCH -Jshb-test-jobs            # Job name
#SBATCH --account=gts-sbryngelson3 # charge account
#SBATCH -N1 --ntasks-per-node=12   # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=2G           # Memory per core
#SBATCH -t 04:00:00                # Duration of the job (Ex: 15 mins)
#SBATCH -q embers                  # QOS Name
#SBATCH -otest.out                 # Combined output and error messages file
#SBATCH -W                         # Do not exit until the submitted job terminates.

cd "$SLURM_SUBMIT_DIR"
echo "Running in $(pwd):"

. ./mfc.sh load -c p -m gpu
./mfc.sh test -j $(nproc) -b mpirun -a

