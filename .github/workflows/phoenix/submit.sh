#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 [script.sh] [cpu|gpu]"
    exit 1
}

[[ $# -eq 2 ]] || usage

sbatch_script="$1"
device="$2"

# read the body of the user script
sbatch_body=$(<"$sbatch_script")

# common SBATCH directives
sbatch_common_opts="\
#SBATCH -J shb-${sbatch_script%%.sh}-$device   # job name
#SBATCH --account=gts-sbryngelson3             # account
#SBATCH -N1                                     # nodes
#SBATCH -t 03:00:00                             # walltime
#SBATCH -q embers                               # QOS
#SBATCH -o ${sbatch_script%%.sh}-$device.%j.out  # stdout+stderr
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

	set -e -x

	cd "\$SLURM_SUBMIT_DIR"
	echo "Running in \$(pwd):"

	# load your modules & env
	. ./mfc.sh load -c p -m $device

	# user script contents
	${sbatch_body}
EOT
)

echo "🚀 Submitted SLURM job $JOBID"

# if this wrapper is killed/canceled, make sure SLURM job is cleaned up
trap '[[ -n "${JOBID:-}" ]] && scancel "$JOBID" >/dev/null 2>&1 || :' EXIT

# ────────── Poll until SLURM job finishes ──────────
while :; do
  # First try sacct (may take a moment to appear)
  STATE=$(sacct -j "$JOBID" --format=State --noheader --parsable2 | head -n1)

  # If sacct didn’t return anything yet, check squeue
  if [[ -z "$STATE" ]]; then
    STATE=$(squeue -j "$JOBID" -h -o "%T" || echo "")
  fi

  # When both are empty, the job has left the queue
  if [[ -z "$STATE" ]]; then
    echo "✅ SLURM job $JOBID no longer in queue; assuming finished"
    break
  fi

  echo "⏳ SLURM job $JOBID state: $STATE"
  sleep 10
done

# show final report
sacct -j "$JOBID" --format=JobID,State,ExitCode,Elapsed,MaxRSS

# exit with the job's real code (left of the colon)
EXIT_CODE=$(sacct -j "$JOBID" --noheader --format=ExitCode | head -1 | cut -d: -f1)
exit "$EXIT_CODE"
