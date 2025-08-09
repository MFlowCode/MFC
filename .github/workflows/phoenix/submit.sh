#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 [cpu|gpu]"
  exit 1
}

[[ $# -eq 1 ]] || usage

device="$1"
job_slug="test-$device"

# Build sbatch arguments (use CLI args instead of #SBATCH lines)
sbatch_args=(
  -J "MFC-test-$device"
  --account=gts-sbryngelson3
  -N 1
  -t 03:00:00
  -q embers
  -o "${job_slug}.out"
  --mem-per-cpu=2G
  # Export variables for the job environment
  --export=ALL,job_slug="$job_slug",device="$device"
)

if [[ "$device" == "cpu" ]]; then
  sbatch_args+=(
    -p cpu-small
    --ntasks-per-node=24
  )
elif [[ "$device" == "gpu" ]]; then
  sbatch_args+=(
    -p gpu-v100,gpu-a100,gpu-h100,gpu-l40s
    --ntasks-per-node=4
    -G 2
  )
else
  usage
fi

# submit and capture the JobID
JOBID=$(
  sbatch "${sbatch_args[@]}" <<'EOT' | awk '{print $4}'
#!/usr/bin/env bash
set -euo pipefail
set -x

echo "Job slug is: $job_slug"
echo "Device is:   $device"

cd "$SLURM_SUBMIT_DIR"
echo "Running in $(pwd)"

# load your modules & env
. ./mfc.sh load -c p -m "$device"

# user script contents
tmpbuild=/storage/scratch1/6/sbryngelson3/mytmp_build
mkdir -p "$tmpbuild"
currentdir="$tmpbuild/run-$(( RANDOM % 900 ))"
mkdir -p "$currentdir"
export TMPDIR="$currentdir"

n_test_threads=8
build_opts=""
if [[ "$device" == "gpu" ]]; then
  build_opts="--gpu"
fi
echo "build_opts = $build_opts"

if [[ "$device" == "cpu" ]]; then
  echo "CPU BUILD"
elif [[ "$device" == "gpu" ]]; then
  echo "GPU BUILD"
else
  echo "Unknown device: $device" >&2
  exit 1
fi

# Dry run (kept from your original)
./mfc.sh test --dry-run -j "$n_test_threads" $build_opts

# GPU-specific runtime options
device_opts=""
if [[ "$device" == "gpu" ]]; then
  if command -v nvidia-smi >/dev/null 2>&1; then
    gpu_count=$(nvidia-smi -L | wc -l)
  else
    gpu_count=0
  fi

  if [[ "$gpu_count" -gt 0 ]]; then
    gpu_ids=$(seq -s ' ' 0 $(( gpu_count - 1 )))
    device_opts="-g $gpu_ids"
    n_test_threads=$(( gpu_count * 2 ))
  else
    echo "No GPUs detected; continuing without -g list"
    device_opts=""
  fi
fi

./mfc.sh test --max-attempts 3 -a -j "$n_test_threads" ${device_opts:-} -- -c phoenix

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
# (small grace period in case accounting lags)
sleep 2
EXIT_CODE=$(sacct -j "$JOBID" --noheader --format=ExitCode | head -1 | cut -d: -f1 || echo 1)
echo "Completed: SLURM job $JOBID exit code: $EXIT_CODE"
exit "$EXIT_CODE"
