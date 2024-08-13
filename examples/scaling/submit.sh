#!/bin/bash

mkdir -p examples/scaling/logs

for N in 1 2 4 8 16 32; do
    echo -n "N=$N: "
    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name MFC-S-$N
#SBATCH --account=CFD154
#SBATCH --nodes=$N
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:20:00
#SBATCH --qos=debug
#SBATCH --output=MFC-S-$N.out
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
#SBATCH --mail-user="hberre3@gatech.edu"
#SBATCH --mail-type="BEGIN, END, FAIL"

set -e
set -x

. ./mfc.sh load -c f -m g

cd "\$SLURM_SUBMIT_DIR"

echo "Hosts"

srun hostname

echo "Start @ $(date)"

for M in 256 384 512; do

    slug=$N-\$M
    case_dir=examples/scaling/run-\$slug
    mkdir -p \$case_dir
    cp examples/scaling/case.py \$case_dir/

    if [ ! -d \$case_dir/restart_data ]; then
        echo "Running pre_process"

        # Note: `time` is not used for performance measurement, only for monitoring
        #       the job's progress.
        time ./mfc.sh run \$case_dir/case.py -c frontier -n 8 -N $N --clean     \
                -t pre_process --no-build -# run-\$slug-pre -- --scaling strong \
                --memory \$M > examples/scaling/logs/run-\$slug-pre.out 2>&1
    fi

    for cu_mpi in F T; do

        slug="$N-\$M-\$cu_mpi"
        echo "Running \$slug"

        # Note: `time` is not used for performance measurement, only for monitoring
        #       the job's progress.
        time ./mfc.sh run \$case_dir/case.py -c frontier -n 8 -N $N         \
            -t simulation --case-optimization --no-build -# run-\$slug-sim  \
            -- --scaling strong --memory \$M --cu_mpi \$cu_mpi --n-steps 20 \
            > examples/scaling/logs/run-\$slug-sim.out 2>&1

    done

    rm -rf \$case_dir

done

echo "End @ $(date)"

EOT
done