#!/bin/bash
# Restart a GitHub Actions runner on a specific Phoenix login node.
#
# Kills any existing instance of the runner, then starts a new one.
# Uses 'bash -l' so the runner inherits the full login PATH (which
# includes /opt/slurm/current/bin — required for sbatch/squeue/sacct).
# Uses 'setsid' + stdin close for full terminal detachment so the SSH
# session exits cleanly without waiting for the runner process.
#
# The runner binary lives on shared storage, so no files need to be
# copied — only the process needs to run on the target node.
#
# Usage: bash restart-runner.sh <node> <runner-dir>
# Example: bash restart-runner.sh login-phoenix-gnr-2 /storage/scratch1/6/sbryngelson3/mfc-runners/actions-runner-3

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

if [ $# -ne 2 ]; then
    echo "Usage: $0 <node> <runner-dir>"
    echo "Example: $0 login-phoenix-gnr-2 /storage/scratch1/6/sbryngelson3/mfc-runners/actions-runner-3"
    exit 1
fi

node="$1"
dir="$2"
name=$(basename "$dir")

echo "Restarting $name on $node..."

stop_runner "$node" "$dir"

if start_runner "$node" "$dir"; then
    pid=$(ssh -o ConnectTimeout=5 "$node" "pgrep -f 'Runner.Listener.*$dir' | head -1" 2>/dev/null || true)
    if [ -n "$pid" ] && check_slurm_path "$node" "$pid"; then
        echo "  OK: PID $pid, slurm in PATH"
    elif [ -n "$pid" ]; then
        echo "  WARNING: PID $pid started but slurm NOT in PATH!"
        echo "  The runner may not be able to submit SLURM jobs."
    fi
else
    echo "  ERROR: Runner failed to start on $node"
    echo "  Try manually: ssh $node 'cd $dir && setsid bash -lc \"nohup ./run.sh >> runner-nohup.log 2>&1 &\"'"
fi
