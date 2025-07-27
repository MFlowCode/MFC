#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 [script.sh] [cpu|gpu]"
    exit 1
}

[[ $# -eq 2 ]] || usage

sbatch_script="$1"

device="$2"
job_slug="`basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g'`-$2"

# read the body of the user script
sbatch_body=$(<"$sbatch_script")

# common SBATCH directives
sbatch_common_opts="\
#SBATCH -J shb-${sbatch_script%%.sh}-$device    # job name
#SBATCH --account=gts-sbryngelson3              # account
#SBATCH -N1                                     # nodes
#SBATCH -t 02:00:00                             # walltime
#SBATCH -q embers                               # QOS
#SBATCH -o $job_slug.out                        # stdout+stderr
#SBATCH --mem-per-cpu=2G                        # default mem (overridden below)
"

# CPU vs GPU overrides
if [[ "$device" == "cpu" ]]; then
  sbatch_device_opts="\
#SBATCH -p cpu-small               # partition
#SBATCH --ntasks-per-node=24       # Number of cores per node required
#SBATCH --mem-per-cpu=2G           # Memory per core\
"
elif [[ "$device" == "gpu" ]]; then
  sbatch_device_opts="\
#SBATCH -CL40S
#SBATCH --ntasks-per-node=4       # Number of cores per node required
#SBATCH -G2\
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
	${sbatch_body}
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
