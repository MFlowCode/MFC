#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --account=bdwd-delta-cpu
#SBATCH --job-name="matlab"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=log/matlab.out%j
#SBATCH --error=log/matlab.err%j
#SBATCH --export=ALL

echo "#######################################"
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes"
echo "Running on $SLURM_NPROCS processors"
echo "#######################################"

/sw/external/matlab/2024a/bin/matlab -nodesktop -nodisplay -nosplash -r "run run_turbulence.m; run average_tke_over_self_similar.m; exit;"

echo "#######################################"
echo "Program finished with exit code $? at: `date`"
echo "#######################################"
