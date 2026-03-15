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

for i in "${!TARGET_NODES[@]}"; do
    NODE="${TARGET_NODES[$i]}"
    NUM=$((START_NUM + i))
    echo "==> Deploying frontier-${NUM} on ${NODE}..."
    "$SCRIPT_DIR/make-runner.sh" "${NUM}" "${NODE}" &
done

wait
echo "==> All runners deployed."
