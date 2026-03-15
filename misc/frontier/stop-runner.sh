#!/usr/bin/env bash
# Stop and deregister a GitHub Actions runner.
# Usage: ./stop-runner.sh <runner-name>  (e.g. frontier-12)
set -euo pipefail

SSH_OPTS="-o StrictHostKeyChecking=no -o ConnectTimeout=30 -o ServerAliveInterval=10 -o ServerAliveCountMax=3"

RUNNER_NAME="${1:?Usage: $0 <runner-name>}"
SHARED_DIR="/lustre/orion/cfd154/proj-shared/runners"
RUNNER_DIR="${SHARED_DIR}/${RUNNER_NAME}"
ORG="MFlowCode"

if [ ! -d "${RUNNER_DIR}" ]; then
    echo "Runner dir not found: ${RUNNER_DIR}"
    exit 1
fi

# --- Kill the process (SSH to correct node if needed) ---
PID_FILE="${RUNNER_DIR}/runner.pid"
NODE_FILE="${RUNNER_DIR}/runner.node"
CURRENT_NODE=$(hostname -s)

if [ -f "${PID_FILE}" ]; then
    PID=$(cat "${PID_FILE}")
    TARGET_NODE="${CURRENT_NODE}"
    [ -f "${NODE_FILE}" ] && TARGET_NODE=$(cat "${NODE_FILE}")

    echo "==> Killing PID ${PID} on ${TARGET_NODE}..."
    if [ "${TARGET_NODE}" = "${CURRENT_NODE}" ]; then
        kill "${PID}" 2>/dev/null || true
    else
        ssh ${SSH_OPTS} "${TARGET_NODE}" "kill ${PID} 2>/dev/null || true"
    fi
    rm -f "${PID_FILE}"
fi

# --- Deregister from GitHub ---
echo "==> Fetching removal token..."
REMOVE_TOKEN=$(gh api \
    --method POST \
    -H "Accept: application/vnd.github+json" \
    "/orgs/${ORG}/actions/runners/remove-token" \
    --jq '.token')

echo "==> Deregistering runner..."
"${RUNNER_DIR}/config.sh" remove --token "${REMOVE_TOKEN}" || true
echo "==> Done."
