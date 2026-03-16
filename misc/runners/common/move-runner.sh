#!/usr/bin/env bash
# Move a runner to a different login node.
#
# Sourced by site wrappers (frontier/move-runner.sh, phoenix/move-runner.sh)
# after config.sh is loaded. Finds the runner by name, stops it on its current
# node, and starts it on the target node. Retries start once after 5 seconds.
#
# Usage: bash move-runner.sh <runner-name> <target-node>
set -euo pipefail

RUNNER_NAME="${1:?Usage: $0 <runner-name> <target-node>}"
TARGET_NODE="${2:?Usage: $0 <runner-name> <target-node>}"

# Validate target node
valid=0
for node in "${NODES[@]}"; do
    [ "$node" = "$TARGET_NODE" ] && valid=1 && break
done
if [ "$valid" -eq 0 ]; then
    echo "ERROR: '$TARGET_NODE' is not a valid login node." >&2
    echo "       Valid nodes: ${NODES[*]}" >&2
    exit 1
fi

# Find runner directory by name
runner_dir=""
while IFS= read -r dir; do
    if [ "$(get_runner_name "$dir")" = "$RUNNER_NAME" ]; then
        runner_dir="$dir"
        break
    fi
done < <(find_runner_dirs)

if [ -z "$runner_dir" ]; then
    echo "ERROR: Runner '$RUNNER_NAME' not found in known runner directories." >&2
    exit 1
fi

declare -f sync_runner_nodes > /dev/null 2>&1 && {
    echo "==> Syncing runner node locations..."
    sync_runner_nodes
}

echo "==> Locating $RUNNER_NAME..."
current_node=$(find_node "$runner_dir")

if [ "$current_node" = "$TARGET_NODE" ]; then
    echo "==> $RUNNER_NAME is already running on $TARGET_NODE. Nothing to do."
    exit 0
fi

if [ "$current_node" != "offline" ]; then
    echo "==> Stopping $RUNNER_NAME on $current_node..."
    stop_runner "$current_node" "$runner_dir"
fi

echo "==> Starting $RUNNER_NAME on $TARGET_NODE..."
if start_runner "$TARGET_NODE" "$runner_dir"; then
    echo "$TARGET_NODE" > "$runner_dir/runner.node"
    echo "==> $RUNNER_NAME is now running on $TARGET_NODE."
else
    echo "    First start attempt failed. Retrying in 5 seconds..."
    sleep 5
    if start_runner "$TARGET_NODE" "$runner_dir"; then
        echo "$TARGET_NODE" > "$runner_dir/runner.node"
        echo "==> $RUNNER_NAME is now running on $TARGET_NODE."
    else
        echo "ERROR: $RUNNER_NAME failed to start on $TARGET_NODE after retry." >&2
        exit 1
    fi
fi
