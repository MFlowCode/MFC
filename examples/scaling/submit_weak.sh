#!/bin/bash

# Initialize default values for optional arguments
NODES=(16 128 1024 8192)
MEM=(64) # Approximate problem size per GCD in GB

# Mandatory argument
ACCOUNT=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --account)
            ACCOUNT="$2"
            shift
            shift
            ;;
        --nodes)
            # Accept a comma-separated list or space-separated
            IFS=',' read -r -a NODES <<< "$2"
            shift
            shift
            ;;
        --mem)
            IFS=',' read -r -a MEM <<< "$2"
            shift
            shift
            ;;
        *)
            echo "Unknown option $1"
            shift
            ;;
    esac
done

# Check mandatory argument
if [[ -z "$ACCOUNT" ]]; then
    echo "Error: --account is required"
    exit 1
fi

# Print results
echo "ACCOUNT = $ACCOUNT"
echo "NODES = ${NODES[@]}"
echo "MEM = ${MEM[@]}"

mkdir -p examples/scaling/logs

for N in "${NODES[@]}"; do
    for M in "${MEM[@]}"; do
        echo -n "N=$N: M=$M"
        sbatch <<EOT
#!/bin/bash
#SBATCH --job-name MFC-W-$N-$M
#SBATCH --account=$ACCOUNT
#SBATCH --nodes=$N
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:20:00
#SBATCH --output=MFC-W-$N-$M.out
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest

set -e
set -x

. ./mfc.sh load -c f -m g

cd "\$SLURM_SUBMIT_DIR"

echo "Hosts"

srun hostname

echo "Start @ $(date)"

slug=$N-$M
case_dir=examples/scaling/weak-\$slug
mkdir -p \$case_dir
cp examples/scaling/case.py \$case_dir/

if [ ! -d \$case_dir/restart_data ]; then
    echo "Running pre_process"

    # Note: `time` is not used for performance measurement, only for monitoring
    #       the job's progress.
    time ./mfc.sh run \$case_dir/case.py -c frontier -n 8 -N $N --clean \
            -t pre_process -# weak-\$slug-pre -- --scaling weak \
            --memory $M \
            > examples/scaling/logs/weak-\$slug-pre.out 2>&1
fi

for rdma_mpi in F T; do

    slug="$N-$M-\$rdma_mpi"
    echo "Running \$slug"

    # Note: `time` is not used for performance measurement, only for monitoring
    #       the job's progress.
    time ./mfc.sh run \$case_dir/case.py -c frontier -n 8 -N $N -t simulation \
        --case-optimization -# weak-\$slug-sim -- --scaling weak \
        --memory $M --rdma_mpi \$rdma_mpi --n-steps 20 \
        > examples/scaling/logs/weak-\$slug-sim.out 2>&1

done

rm -rf \$case_dir

done

echo "End @ $(date)"

EOT
done
done
