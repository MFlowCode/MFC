#!/usr/bin/env bash
# Create, configure, and start a single GitHub Actions runner on Frontier.
# Usage: make-runner.sh <runner-number> [login-node]
#   runner-number  Sequential number for this runner (e.g. 12)
#   login-node     Node to run on (default: current host)
# Example: make-runner.sh 12 login03
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

RUNNER_VERSION="${RUNNER_VERSION:-$(gh_latest_runner_version 2>/dev/null || echo "2.332.0")}"

RUNNER_NUM="${1:?Usage: $0 <runner-number> [login-node]}"
TARGET_NODE="${2:-$(hostname -s)}"

TARBALL="actions-runner-linux-x64-${RUNNER_VERSION}.tar.gz"
RUNNER_NAME="frontier-${RUNNER_NUM}"
RUNNER_DIR="${SHARED_DIR}/${RUNNER_NAME}"

echo "==> Using runner version ${RUNNER_VERSION}"
echo "==> Setting up runner: ${RUNNER_NAME} on ${TARGET_NODE}"

# --- Download tarball once to shared dir ---
if [ ! -f "${SHARED_DIR}/${TARBALL}" ]; then
    echo "==> Downloading runner v${RUNNER_VERSION}..."
    tmp="${SHARED_DIR}/${TARBALL}.tmp.$$"
    curl -fsSL \
        "https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/${TARBALL}" \
        -o "$tmp"
    mv "$tmp" "${SHARED_DIR}/${TARBALL}"
fi

# --- Extract (filesystem is shared across all nodes) ---
mkdir -p "${RUNNER_DIR}"
echo "==> Extracting runner into ${RUNNER_DIR}..."
tar xzf "${SHARED_DIR}/${TARBALL}" -C "${RUNNER_DIR}"

if [ ! -f "${RUNNER_DIR}/run.sh" ]; then
    echo "ERROR: Extraction failed — run.sh not found in ${RUNNER_DIR}" >&2
    exit 1
fi

# --- Configure ---
echo "==> Fetching registration token..."
REG_TOKEN=$(gh_registration_token)

echo "==> Configuring runner..."
"${RUNNER_DIR}/config.sh" \
    --url "https://github.com/${ORG}" \
    --token "${REG_TOKEN}" \
    --name "${RUNNER_NAME}" \
    --labels "${RUNNER_LABEL}" \
    --runnergroup "${RUNNER_GROUP}" \
    --work "_work" \
    --unattended \
    --replace

# --- Store which node this runner lives on ---
echo "${TARGET_NODE}" > "${RUNNER_DIR}/runner.node"

# --- Start runner on target node ---
echo "==> Starting runner on ${TARGET_NODE}..."
if start_runner "${TARGET_NODE}" "${RUNNER_DIR}"; then
    echo "==> Runner '${RUNNER_NAME}' is running on ${TARGET_NODE}."
else
    echo "ERROR: Runner '${RUNNER_NAME}' did not start on ${TARGET_NODE}." >&2
    exit 1
fi

echo "==> Log: ${RUNNER_DIR}/runner.log"
