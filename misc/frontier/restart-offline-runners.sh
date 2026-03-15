#!/usr/bin/env bash
# Restart all offline frontier runners.
#
# Queries GitHub for offline frontier runners, locates each via CWD-based
# process discovery, stops any stale processes, then restarts in parallel.
# Falls back to runner.node for the target node if the runner is truly offline.
#
# Usage: bash restart-offline-runners.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

echo "==> Checking for offline frontier runners..."

# Collect offline runner names from GitHub API
mapfile -t OFFLINE_NAMES < <(
    gh_list_runners | while read -r id name status busy; do
        [ "$status" = "offline" ] && echo "$name"
    done
)

if [ ${#OFFLINE_NAMES[@]} -eq 0 ]; then
    echo "==> All frontier runners are online. Nothing to do."
    exit 0
fi

echo "==> Offline runners: ${OFFLINE_NAMES[*]}"

restart_one() {
    local runner_name="$1"
    local dir="${SHARED_DIR}/${runner_name}"

    if [ ! -d "$dir" ]; then
        echo "WARN: No directory for ${runner_name}, skipping."
        return
    fi

    # Check if it's actually already running somewhere (GitHub may lag)
    local actual_node
    actual_node=$(find_node "$dir")

    if [ "$actual_node" != "offline" ]; then
        echo "==> ${runner_name} appears running on ${actual_node} (GitHub may lag) — stopping first..."
        stop_runner "$actual_node" "$dir"
    fi

    # Determine target node from runner.node fallback
    local target_node
    if [ -f "${dir}/runner.node" ]; then
        target_node=$(cat "${dir}/runner.node")
    else
        echo "WARN: No runner.node for ${runner_name}, skipping."
        return
    fi

    echo "==> Starting ${runner_name} on ${target_node}..."
    if start_runner "$target_node" "$dir"; then
        echo "    ${runner_name}: started on ${target_node}."
    else
        echo "    ${runner_name}: ERROR — failed to start on ${target_node}." >&2
    fi
}

# Restart all offline runners in parallel
for name in "${OFFLINE_NAMES[@]}"; do
    restart_one "$name" &
done

wait

echo ""
echo "==> Final status:"
gh_list_runners | while read -r id name status busy; do
    printf "  %-30s %s\n" "$name" "$status"
done
