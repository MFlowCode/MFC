#!/usr/bin/env bash
# Stop and deregister a GitHub Actions runner on Frontier.
# Usage: stop-runner.sh <runner-name>  (e.g. frontier-12)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

RUNNER_NAME="${1:?Usage: $0 <runner-name>}"
RUNNER_DIR="${SHARED_DIR}/${RUNNER_NAME}"

if [ ! -d "${RUNNER_DIR}" ]; then
    echo "Runner dir not found: ${RUNNER_DIR}" >&2
    exit 1
fi

# --- Locate and kill the runner process ---
echo "==> Locating ${RUNNER_NAME}..."
node=$(find_node "$RUNNER_DIR")

if [ "$node" != "offline" ]; then
    echo "==> Stopping ${RUNNER_NAME} on ${node}..."
    stop_runner "$node" "$RUNNER_DIR"
    echo "==> Process stopped."
else
    echo "==> ${RUNNER_NAME} is not running (already offline)."
fi

# --- Deregister from GitHub ---
echo "==> Fetching runner ID from GitHub..."
runner_id=$(gh_list_runners | while read -r id name status busy; do
    [ "$name" = "$RUNNER_NAME" ] && echo "$id" && break
done)

if [ -n "$runner_id" ]; then
    echo "==> Deregistering runner (ID ${runner_id})..."
    gh_remove_runner "$runner_id"
    echo "==> Done."
else
    echo "==> Runner not found in GitHub API (may already be deregistered)."
fi
