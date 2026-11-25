#!/bin/bash
# Weak and strong scaling of examples/scaling/benchmark.py on Frontier (AMD flang, OpenMP offload, GPU-aware MPI,
# case-optimized: that fixes the numerics at compile time, not the grid, so one build serves every size).
#
#   scaling.sh <log dir> [node counts...]        default: 1 8 64
#
# One SLURM job per node count, submitted together; each runs the weak case (SCALING_WEAK_GB per
# GPU) and the strong case (SCALING_STRONG_GB in all) for 20 steps and leaves
# <log dir>/<weak|strong>-<nodes>.out for scaling-record.py. A node count whose job fails is
# missing from the result; the script fails only when none finished.
set -o pipefail            # not -u: mfc.sh load reads variables it has not set

logs=${1:?log directory}; shift
nodes=("$@"); [ ${#nodes[@]} -gt 0 ] || nodes=(1 8 64)
weak_gb=${SCALING_WEAK_GB:-16}
strong_gb=${SCALING_STRONG_GB:-256}
account=${SCALING_ACCOUNT:-CFD154}

. ./mfc.sh load -c famd -m g
if [ -z "${SCALING_NO_BUILD:-}" ]; then
    ./mfc.sh build -j 8 --gpu mp --no-debug -t pre_process simulation --case-optimization \
        -i examples/scaling/benchmark.py -- --scaling weak --memory "$weak_gb" || exit 1
fi

rm -rf "$logs" scaling-run; mkdir -p "$logs" scaling-run
for N in "${nodes[@]}"; do
    sbatch --wait <<EOT > "$logs/submit-$N.out" 2>&1 &
#!/bin/bash
#SBATCH --job-name=MFC-scaling-$N
#SBATCH --account=$account
#SBATCH --nodes=$N
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
#SBATCH --time=00:40:00
#SBATCH --output=$logs/job-$N.out

cd "\$SLURM_SUBMIT_DIR"
. ./mfc.sh load -c famd -m g

for kind in weak strong; do
    mem=\$([ \$kind = weak ] && echo $weak_gb || echo $strong_gb)
    dir=scaling-run/\$kind-$N
    mkdir -p \$dir && cp examples/scaling/benchmark.py \$dir/case.py
    ./mfc.sh run \$dir/case.py -c frontier_amd -n 8 -N $N --no-build -t pre_process \\
        -- --scaling \$kind --memory \$mem > $logs/\$kind-$N-pre.out 2>&1 || continue
    ./mfc.sh run \$dir/case.py -c frontier_amd -n 8 -N $N --no-build -t simulation --case-optimization \\
        -- --scaling \$kind --memory \$mem --rdma_mpi T --n-steps 20 > $logs/\$kind-$N.out 2>&1
    rm -rf \$dir
done
EOT
done
wait

# What each SLURM job cost, in the form the CI history records.
for N in "${nodes[@]}"; do
    id=$(grep -o 'Submitted batch job [0-9]*' "$logs/submit-$N.out" | grep -o '[0-9]*$')
    [ -n "$id" ] && sacct -j "$id" -X -n -P -o JobID,JobName,NNodes,ElapsedRaw,State | sed 's/^/::notice title=scaling-usage::/'
done

grep -l "Performance:" "$logs"/weak-*.out "$logs"/strong-*.out > /dev/null 2>&1
