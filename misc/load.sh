
echo -n "Would you like to load modules for Summit (s) or Ascent (a)? (s/a): "
read response

if [ "$response" == "s" ]; then # For Summit
    module load lsf-tools/2.0
    module load darshan-runtime/3.3.1-lite
    module load DefApps
    module load spectrum-mpi/10.4.0.3-20210112
    module load hsi/5.0.2.p5
    module load xalt/1.2.1
    module load nvhpc/21.9
    module load cuda/11.4.2
elif [ "$response" == "a" ]; then # For Ascent
    module purge
    module load nvhpc/21.9
    module load spectrum-mpi
    module load cuda/11.2.0
    module load nsight-compute
    module load nsight-systems
else
    echo "Error: Requested system is not supported"
    exit 1
fi
