#!/usr/bin/env bash
# Move a Frontier runner to a different login node.
#
# Stops the runner on its current node, updates runner.node, and starts it on
# the target node. Retries the start once after 5 seconds if the first attempt
# fails.
#
# Usage: move-runner.sh <runner-name> <target-node>
# Example: move-runner.sh frontier-1 login01
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

RUNNER_NAME="${1:?Usage: $0 <runner-name> <target-node>}"
TARGET_NODE="${2:?Usage: $0 <runner-name> <target-node>}"

RUNNER_DIR="${SHARED_DIR}/${RUNNER_NAME}"

# --- Validate runner directory ---
if [ ! -d "$RUNNER_DIR" ]; then
    echo "ERROR: Runner directory not found: ${RUNNER_DIR}" >&2
    exit 1
fi

# --- Validate target node is in the known node list ---
valid=0
for node in "${NODES[@]}"; do
    [ "$node" = "$TARGET_NODE" ] && valid=1 && break
done
if [ "$valid" -eq 0 ]; then
    echo "ERROR: '${TARGET_NODE}' is not a valid Frontier login node." >&2
    echo "       Valid nodes: ${NODES[*]}" >&2
    exit 1
fi

# --- Find current node ---
echo "==> Locating ${RUNNER_NAME}..."
current_node=$(find_node "$RUNNER_DIR")

if [ "$current_node" = "$TARGET_NODE" ]; then
    echo "==> ${RUNNER_NAME} is already running on ${TARGET_NODE}. Nothing to do."
    exit 0
fi

# --- Stop runner on current node (if running) ---
if [ "$current_node" != "offline" ]; then
    echo "==> Stopping ${RUNNER_NAME} on ${current_node}..."
    stop_runner "$current_node" "$RUNNER_DIR"
fi

# --- Update runner.node ---
echo "$TARGET_NODE" > "${RUNNER_DIR}/runner.node"

# --- Start runner on target node (with one retry) ---
echo "==> Starting ${RUNNER_NAME} on ${TARGET_NODE}..."
if start_runner "$TARGET_NODE" "$RUNNER_DIR"; then
    echo "==> ${RUNNER_NAME} is now running on ${TARGET_NODE}."
else
    echo "    First start attempt failed. Retrying in 5 seconds..."
    sleep 5
    if start_runner "$TARGET_NODE" "$RUNNER_DIR"; then
        echo "==> ${RUNNER_NAME} is now running on ${TARGET_NODE}."
    else
        echo "ERROR: ${RUNNER_NAME} failed to start on ${TARGET_NODE} after retry." >&2
        exit 1
    fi
fi
