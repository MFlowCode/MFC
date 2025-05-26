#!/bin/sh
#PBS -A ONRDC53612738
#PBS -q standard
#PBS -l matlab=1
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=192:mpiprocs=2
#PBS -o log/matlab.out
#PBS -e log/matlab.err

echo "########################################"
echo "Starting at `date`"
echo "Running on hosts: $BC_NODE_TYPE"
echo "Running on $BC_NODE_ALLOC nodes."
echo "Running on $BC_MPI_TASKS_ALLOC processors."
echo "Current working directory is `pwd`"
echo "########################################"

/app/matlab/R2023b/bin/matlab -nodisplay -nosplash -nodesktop -r "run('/p/home/hyeoksu/MFC-post-matlab/turbulence/run_turbulence.m');exit;"

echo "########################################"
echo "Program finished with exit code $? at: `date`"
echo "########################################"
