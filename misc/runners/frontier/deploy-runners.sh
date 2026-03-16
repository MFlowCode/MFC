#!/usr/bin/env bash
# Deploy one runner per login node in parallel.
# Usage: deploy-runners.sh <start-num> <node1> [node2 ...]
# Example: deploy-runners.sh 17 login08 login09 login10
#   Deploys frontier-17 on login08, frontier-18 on login09, frontier-19 on login10.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

START_NUM="${1:?Usage: $0 <start-num> <node1> [node2 ...]}"
shift
TARGET_NODES=("$@")

if [ ${#TARGET_NODES[@]} -eq 0 ]; then
    echo "Error: no login nodes specified." >&2
    echo "Usage: $0 <start-num> <node1> [node2 ...]" >&2
    exit 1
fi

# Pre-download the runner tarball once before spawning parallel make-runner.sh
# instances. Without this, all instances race to download the same file
# concurrently and corrupt it. The tmp+mv ensures an atomic final placement.
RUNNER_VERSION="${RUNNER_VERSION:-$(gh_latest_runner_version 2>/dev/null || echo "2.332.0")}"
TARBALL="actions-runner-linux-x64-${RUNNER_VERSION}.tar.gz"
if [ ! -f "${TARBALL_CACHE_DIR}/${TARBALL}" ]; then
    echo "==> Downloading runner v${RUNNER_VERSION}..."
    tmp="${TARBALL_CACHE_DIR}/${TARBALL}.tmp.$$"
    curl -fsSL \
        "https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/${TARBALL}" \
        -o "$tmp"
    mv "$tmp" "${TARBALL_CACHE_DIR}/${TARBALL}"
fi
export RUNNER_VERSION

declare -a pids=() nums=() deploy_nodes=()
for i in "${!TARGET_NODES[@]}"; do
    NODE="${TARGET_NODES[$i]}"
    NUM=$((START_NUM + i))
    echo "==> Deploying frontier-${NUM} on ${NODE}..."
    "$SCRIPT_DIR/make-runner.sh" "${NUM}" "${NODE}" &
    pids+=($!)
    nums+=("$NUM")
    deploy_nodes+=("$NODE")
done

failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        echo "ERROR: frontier-${nums[$i]} on ${deploy_nodes[$i]} failed." >&2
        failed=$((failed + 1))
    fi
done

if [ "$failed" -gt 0 ]; then
    echo "==> $failed runner(s) failed to deploy." >&2
    exit 1
fi
echo "==> All runners deployed."
