#!/usr/bin/env bash
# Create, configure, and start a single GitHub Actions runner.
# Usage: ./make-runner.sh <runner-number> [login-node]
#   runner-number  Sequential number for this runner (e.g. 12)
#   login-node     Node to run on (default: current host). Runner will be
#                  started there via SSH (or locally if it matches current host).
# Example: ./make-runner.sh 12 login03
set -euo pipefail

SSH_OPTS="-o StrictHostKeyChecking=no -o ConnectTimeout=30 -o ServerAliveInterval=10 -o ServerAliveCountMax=3"

RUNNER_NUM="${1:?Usage: $0 <runner-number> [login-node]}"
TARGET_NODE="${2:-$(hostname -s)}"
RUNNER_VERSION="2.332.0"
TARBALL="actions-runner-linux-x64-${RUNNER_VERSION}.tar.gz"
SHARED_DIR="/lustre/orion/cfd154/proj-shared/runners"
ORG="MFlowCode"
RUNNER_GROUP="phoenix"
EXTRA_LABELS="frontier"
RUNNER_NAME="frontier-${RUNNER_NUM}"
RUNNER_DIR="${SHARED_DIR}/${RUNNER_NAME}"

echo "==> Setting up runner: ${RUNNER_NAME} on ${TARGET_NODE}"

# --- Download tarball once to shared dir ---
if [ ! -f "${SHARED_DIR}/${TARBALL}" ]; then
    echo "==> Downloading runner v${RUNNER_VERSION}..."
    curl -fsSL \
        "https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/${TARBALL}" \
        -o "${SHARED_DIR}/${TARBALL}"
fi

# --- Extract (always run locally — filesystem is shared across all nodes) ---
mkdir -p "${RUNNER_DIR}"
echo "==> Extracting runner into ${RUNNER_DIR}..."
tar xzf "${SHARED_DIR}/${TARBALL}" -C "${RUNNER_DIR}"

# Verify extraction succeeded
if [ ! -f "${RUNNER_DIR}/run.sh" ]; then
    echo "ERROR: Extraction failed — run.sh not found in ${RUNNER_DIR}" >&2
    exit 1
fi

# --- Configure (always run locally — filesystem is shared across all nodes) ---
echo "==> Fetching registration token..."
REG_TOKEN=$(gh api \
    --method POST \
    -H "Accept: application/vnd.github+json" \
    "/orgs/${ORG}/actions/runners/registration-token" \
    --jq '.token')

echo "==> Configuring runner..."
"${RUNNER_DIR}/config.sh" \
    --url "https://github.com/${ORG}" \
    --token "${REG_TOKEN}" \
    --name "${RUNNER_NAME}" \
    --labels "${EXTRA_LABELS}" \
    --runnergroup "${RUNNER_GROUP}" \
    --work "_work" \
    --unattended \
    --replace

# --- Store which node this runner lives on ---
echo "${TARGET_NODE}" > "${RUNNER_DIR}/runner.node"

# --- Start run.sh on the target node ---
echo "==> Starting runner on ${TARGET_NODE}..."
CURRENT_NODE=$(hostname -s)
if [ "${TARGET_NODE}" = "${CURRENT_NODE}" ]; then
    nohup "${RUNNER_DIR}/run.sh" >> "${RUNNER_DIR}/runner.log" 2>&1 < /dev/null &
    echo $! > "${RUNNER_DIR}/runner.pid"
else
    PID=$(ssh ${SSH_OPTS} "${TARGET_NODE}" \
        "nohup ${RUNNER_DIR}/run.sh >> ${RUNNER_DIR}/runner.log 2>&1 < /dev/null & echo \$!")
    echo "${PID}" > "${RUNNER_DIR}/runner.pid"
fi

echo "==> Runner PID: $(cat ${RUNNER_DIR}/runner.pid)"
echo "==> Log: ${RUNNER_DIR}/runner.log"
echo "==> Done. Runner '${RUNNER_NAME}' is starting on ${TARGET_NODE}."
