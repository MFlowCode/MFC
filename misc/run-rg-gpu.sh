#!/bin/bash
#SBATCH -JMFC-ci                        # Job name
#SBATCH -N1 --ntasks-per-node=4         # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=4G                # Memory per core
#SBATCH -t 02:00:00                     # Duration of the job (Ex: 30 mins)
#SBATCH -p rg-nextgen-hpc               # Partition Name
#SBATCH -o /projects/ci-runners/MFC/mfc-%j.out   # Combined output
#SBATCH --nodelist violet1	# Specify a specific node
#SBATCH -G 1				# Request a GPU on that node
#SBATCH -W                  # Do not exit until the submitted job terminates.

# #SBATCH -o mfc-%j.out   # Combined output

##Add commands here to build and execute
cd $GITHUB_WORKSPACE
hostname

pwd

#This line allows the GH runner to use the module command on the targeted node
source /tools/misc/.read_profile

#Load NVHPC SDK, which includes the latest CUDA support
module load nvhpc

gpu_count=$(nvidia-smi -L | wc -l)        # number of GPUs on node
gpu_ids=$(seq -s ' ' 0 $(($gpu_count-1))) # 0,1,2,...,gpu_count-1

./mfc.sh test -j 1 --gpu -g $gpu_ids
