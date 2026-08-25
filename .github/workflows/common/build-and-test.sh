#!/bin/bash
# Build then test in a single SLURM allocation so the scheduler queue wait is
# paid once instead of twice (Phoenix 'embers' is queue-bound, not work-bound).
# submit-slurm-job.sh exports the job_* vars so these child scripts inherit
# them; cwd is the workspace root, so build/ persists into test.sh's --no-build.

set -euo pipefail

echo "=== [build-and-test] Build phase ==="
bash .github/workflows/common/build.sh

echo "=== [build-and-test] Test phase ==="
bash .github/workflows/common/test.sh
