#!/usr/bin/env bash
# Restart a single GitHub Actions runner on a given node.
#
# Sourced by site wrappers (frontier/restart-runner.sh, phoenix/restart-runner.sh)
# after config.sh is loaded. Stops any existing process, starts fresh, and
# verifies slurm is in PATH.
#
# Usage: bash restart-runner.sh <node> <runner-dir>
set -euo pipefail

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
    echo "$node" > "$dir/runner.node"
    pids=$(find_pids "$node" "$dir")
    pid=${pids%% *}
    if has_slurm "$node" "$pid"; then
        echo "  OK: PID $pid, slurm in PATH"
    else
        echo "  WARNING: PID $pid but slurm MISSING from PATH"
    fi
else
    echo "  ERROR: Failed to start on $node"
    exit 1
fi
