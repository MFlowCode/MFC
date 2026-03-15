#!/bin/bash
# Restart a GitHub Actions runner on a specific Phoenix login node.
#
# Kills any existing instance, then starts a new one with a login shell
# (for /opt/slurm PATH) and full terminal detachment.
#
# Usage: bash restart-runner.sh <node> <runner-dir>
# Example: bash restart-runner.sh login-phoenix-gnr-2 /storage/scratch1/6/sbryngelson3/mfc-runners/actions-runner-3

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

if [ $# -ne 2 ]; then
    echo "Usage: $0 <node> <runner-dir>"
    echo "Nodes: ${NODES[*]}"
    exit 1
fi

node="$1"
dir="$2"
name=$(get_runner_name "$dir" 2>/dev/null || basename "$dir")

echo "Restarting $name on $node..."
stop_runner "$node" "$dir"

if start_runner "$node" "$dir"; then
    pids=$(find_pids "$node" "$dir")
    pid=${pids%% *}
    if has_slurm "$node" "$pid"; then
        echo "  OK: PID $pid, slurm in PATH"
    else
        echo "  WARNING: PID $pid but slurm MISSING from PATH"
    fi
else
    echo "  ERROR: Failed to start on $node"
fi
