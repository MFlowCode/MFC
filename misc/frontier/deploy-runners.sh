#!/usr/bin/env bash
# Deploy one runner per login node.
# Usage: ./deploy-runners.sh <start-num> <node1> [node2 ...]
# Example: ./deploy-runners.sh 17 login08 login09 login10
#   Deploys frontier-17 on login08, frontier-18 on login09, frontier-19 on login10.
set -euo pipefail

SHARED_DIR="/lustre/orion/cfd154/proj-shared/runners"

START_NUM="${1:?Usage: $0 <start-num> <node1> [node2 ...]}"
shift
NODES=("$@")

if [ ${#NODES[@]} -eq 0 ]; then
    echo "Error: no login nodes specified." >&2
    echo "Usage: $0 <start-num> <node1> [node2 ...]" >&2
    exit 1
fi

for i in "${!NODES[@]}"; do
    NODE="${NODES[$i]}"
    NUM=$((START_NUM + i))
    echo "==> Deploying frontier-${NUM} on ${NODE}..."
    bash "${SHARED_DIR}/make-runner.sh" "${NUM}" "${NODE}" &
done

wait
echo "==> All runners deployed."
