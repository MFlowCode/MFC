#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 [cpu|gpu]"
    exit 1
}

[[ $# -eq 1 ]] || usage

device="$1"
job_slug="test-$1"

# common SBATCH directives
sbatch_common_opts="\
#SBATCH -J MFC-test-$device    # job name
#SBATCH --account=gts-sbryngelson3              # account
#SBATCH -N1                                     # nodes
#SBATCH -t 03:00:00                             # walltime
#SBATCH -q embers                               # QOS
#SBATCH -o $job_slug.out                        # stdout+stderr
#SBATCH --mem-per-cpu=2G                        # default mem (overridden below)
"

# CPU vs GPU overrides
if [[ "$device" == "cpu" ]]; then
  sbatch_device_opts="\
#SBATCH -p cpu-small
#SBATCH --ntasks-per-node=24
"
elif [[ "$device" == "gpu" ]]; then
  sbatch_device_opts="\
#SBATCH -p gpu-v100,gpu-a100,gpu-h100,gpu-l40s
#SBATCH --ntasks-per-node=4
#SBATCH -G2
"
else
  usage
fi

# submit and capture the JobID
JOBID=$(sbatch <<-EOT | awk '{print $4}'
	#!/usr/bin/env bash
	${sbatch_common_opts}
	${sbatch_device_opts}

	export job_slug="${job_slug}"
	export device="${device}"

	echo "Job slug is:" $job_slug
 	echo "Device is:" $device

	set -e -x

	cd "\$SLURM_SUBMIT_DIR"
	echo "Running in \$(pwd):"

	# load your modules & env
	. ./mfc.sh load -c p -m $device

	# user script contents
    tmpbuild=/storage/scratch1/6/sbryngelson3/mytmp_build
    currentdir=$tmpbuild/run-$(( RANDOM % 900 ))
    mkdir -p $tmpbuild
    mkdir -p $currentdir
    export TMPDIR=$currentdir

    n_test_threads=8

    build_opts=""
    if [ "$device" = "gpu" ]; then
        build_opts="--gpu"
    fi
    echo "build_opts =" $build_opts

    if [[ "$device" == "cpu" ]]; then
        echo "CPU BUILD"
    elif [[ "$device" == "gpu" ]]; then
        echo "GPU BUILD"
    else
      exit 1
    fi

    exit 1

    ./mfc.sh test --dry-run -j $n_test_threads $build_opts

    if [ "$device" = "gpu" ]; then
        gpu_count=$(nvidia-smi -L | wc -l)        # number of GPUs on node
        gpu_ids=$(seq -s ' ' 0 $(($gpu_count-1))) # 0,1,2,...,gpu_count-1
        device_opts="-g $gpu_ids"
        n_test_threads=`expr $gpu_count \* 2`
    fi

    ./mfc.sh test --max-attempts 3 -a -j $n_test_threads $device_opts -- -c phoenix

    sleep 10
    rm -rf "$currentdir" || true

    unset TMPDIR
EOT
)

echo "Submitted: SLURM job $JOBID"

# if this wrapper is killed/canceled, make sure SLURM job is cleaned up
trap '[[ -n "${JOBID:-}" ]] && scancel "$JOBID" >/dev/null 2>&1 || :' EXIT

# ────────── Poll until SLURM job finishes ──────────
while :; do
  # Try sacct first
  STATE=$(sacct -j "$JOBID" --format=State --noheader --parsable2 | head -n1)

  # Fallback to squeue if sacct is empty
  if [[ -z "$STATE" ]]; then
    STATE=$(squeue -j "$JOBID" -h -o "%T" || echo "")
  fi

  # If it’s one of SLURM’s terminal states, break immediately
  case "$STATE" in
    COMPLETED|FAILED|CANCELLED|TIMEOUT|PREEMPTED)
      echo "Completed: SLURM job $JOBID reached terminal state: $STATE"
      break
      ;;
    "")
      echo "Completed: SLURM job $JOBID no longer in queue; assuming finished"
      break
      ;;
    *)
      echo "Waiting: SLURM job $JOBID state: $STATE"
      sleep 10
      ;;
  esac
done

# Now retrieve the exit code and exit with it
EXIT_CODE=$(sacct -j "$JOBID" --noheader --format=ExitCode | head -1 | cut -d: -f1)
echo "Completed: SLURM job $JOBID exit code: $EXIT_CODE"
exit "$EXIT_CODE"
