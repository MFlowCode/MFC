#!/bin/bash
#SBATCH -Jshb-test-jobs                         # Job name
#SBATCH --account=gts-sbryngelson3               # charge account
#SBATCH -N1                                     # Number of nodes and cores per node required
#SBATCH --gres=gpu:V100:2
#SBATCH -CV100-16GB
#SBATCH -G2
#SBATCH --tasks-per-node=2
#SBATCH -t 02:00:00                              # Duration of the job (Ex: 15 mins)
#SBATCH -q embers                               # QOS Name
#SBATCH -otest.out                               # Combined output and error messages file
#SBATCH -W                                      # Do not exit until the submitted job terminates.

cd $SLURM_SUBMIT_DIR                            # Change to working directory
echo $(pwd)
. ./mfc.sh load -c p -m gpu
./mfc.sh test -b mpirun -a --gpu
