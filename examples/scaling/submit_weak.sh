#!/bin/bash

mkdir -p examples/scaling/logs

NODES=(16 48 128 384 1024 3072 8192)
GBS=(64)
ACCOUNT=cfd154

for N in $NODES; do
    echo -n "N=$N: "
    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name MFC-W-$N
#SBATCH --account=$ACCOUNT
#SBATCH --nodes=$N
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:20:00
#SBATCH --qos=debug
#SBATCH --output=MFC-W-$N.out
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

for M in $GBS; do

    slug=$N-\$M
    case_dir=examples/scaling/run-\$slug
    mkdir -p \$case_dir
    cp examples/scaling/case.py \$case_dir/

    if [ ! -d \$case_dir/restart_data ]; then
        echo "Running pre_process"

        # Note: `time` is not used for performance measurement, only for monitoring
        #       the job's progress.
        time ./mfc.sh run \$case_dir/case.py -c frontier -n 8 -N $N --clean   \
                -t pre_process -# run-\$slug-pre -- --scaling weak \
                --memory \$M \
                > examples/scaling/logs/run-\$slug-pre.out 2>&1
    fi

    for rdma_mpi in F T; do

        slug="$N-\$M-\$rdma_mpi"
        echo "Running \$slug"

        # Note: `time` is not used for performance measurement, only for monitoring
        #       the job's progress.
        time ./mfc.sh run \$case_dir/case.py -c frontier -n 8 -N $N -t simulation \
            --case-optimization -# run-\$slug-sim -- --scaling weak \
            --memory \$M --rdma_mpi \$rdma_mpi --n-steps 20 \
            > examples/scaling/logs/run-\$slug-sim.out 2>&1

    done

    rm -rf \$case_dir

done

echo "End @ $(date)"

EOT
done
