#!/bin/bash
# Submit a SLURM job and wait for it to complete.
# Delegates submission (with idempotency) to submit-job.sh, then monitors.
#
# Usage: submit.sh [script.sh] [cpu|gpu] [none|acc|omp]

set -euo pipefail

# Ignore SIGHUP to survive login node session drops
trap '' HUP

usage() {
    echo "Usage: $0 [script.sh] [cpu|gpu] [none|acc|omp]"
}

if [ -z "${1:-}" ]; then
    usage
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Submit (idempotent — skips resubmission if a live job already exists)
bash "$SCRIPT_DIR/submit-job.sh" "$@"

# Derive the same job slug and file paths as submit-job.sh.
# NOTE: this sed pipeline must stay identical to the one in submit-job.sh —
# if they diverge the id-file will not be found and the monitor will fail.
job_slug="$(basename "$1" | sed 's/\.sh$//' | sed 's/[^a-zA-Z0-9]/-/g')-$2-$3"
output_file="$job_slug.out"
id_file="${job_slug}.slurm_job_id"

job_id=$(cat "$id_file")
bash "$SCRIPT_DIR/../../scripts/run_monitored_slurm_job.sh" "$job_id" "$output_file"
