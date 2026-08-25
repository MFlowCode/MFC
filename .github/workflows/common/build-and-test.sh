#!/bin/bash
# Combined build-then-test for a single SLURM allocation.
# Runs inside a SLURM job via submit-slurm-job.sh. Doing both phases in one
# allocation means the scheduler queue wait is paid once instead of twice —
# on a contended, preemptible QOS (Phoenix 'embers') the queue wait, not the
# work, is what makes CI slow and flaky.
#
# Relies on submit-slurm-job.sh exporting the job_* vars (job_device,
# job_interface, job_cluster, ...) so the child scripts inherit them. cwd is
# $SLURM_SUBMIT_DIR (the workspace root), so build/ persists from build.sh into
# test.sh's --no-build run.

set -euo pipefail

echo "=== [build-and-test] Build phase ==="
bash .github/workflows/common/build.sh

echo "=== [build-and-test] Test phase ==="
bash .github/workflows/common/test.sh
