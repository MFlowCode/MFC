#!/bin/bash

# Initialize default values for optional arguments
MEM=(8 16 32 64) # Approximate memory per device in GB

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
echo "MEM = ${MEM[@]}"

mkdir -p examples/scaling/logs

for M in "${MEM[@]}"; do
    echo "M=$M"
        sbatch <<EOT
#!/bin/bash
#SBATCH --job-name MFC-G-$M
#SBATCH --account=$ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --output=MFC-G-$M.out
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

slug=1-$M
case_dir=examples/scaling/grind-\$slug
mkdir -p \$case_dir
cp examples/scaling/case.py \$case_dir/

if [ ! -d \$case_dir/restart_data ]; then
    echo "Running pre_process"

    # Note: `time` is not used for performance measurement, only for monitoring
    #       the job's progress.
    time ./mfc.sh run \$case_dir/case.py -c frontier -n 1 -N 1 --clean \
            -t pre_process -# grind-\$slug-pre -- --scaling weak \
            --memory $M \
            > examples/scaling/logs/grind-\$slug-pre.out 2>&1
fi

slug="1-$M-F"
echo "Running \$slug"

# Note: `time` is not used for performance measurement, only for monitoring
#       the job's progress.
time ./mfc.sh run \$case_dir/case.py -c frontier -n 1 -N 1 -t simulation \
    --case-optimization -# grind-\$slug-sim -- --scaling weak \
    --memory $M --n-steps 100 --n-save 20 \
    > examples/scaling/logs/grind-\$slug-sim.out 2>&1

rm -rf \$case_dir

done

echo "End @ $(date)"

EOT
done
