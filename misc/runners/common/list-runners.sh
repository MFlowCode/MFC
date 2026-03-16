#!/usr/bin/env bash
# List all runners combining GitHub API status with live node process info.
#
# Sourced by site wrappers (frontier/list-runners.sh, phoenix/list-runners.sh)
# after config.sh is loaded. Uses a parallel SSH sweep across all nodes
# simultaneously (one SSH per node regardless of runner count).
# Shows name, GitHub status, node, slurm availability, and RSS.
# If CGROUP_LIMIT > 0, also shows a per-node memory summary.
#
# Usage: bash list-runners.sh
set -euo pipefail

declare -f sync_runner_nodes > /dev/null 2>&1 && {
    echo "==> Syncing runner node locations..."
    sync_runner_nodes
}

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

sweep_all_nodes "$tmpdir"

# Parse sweep results into associative arrays
declare -A runner_node runner_rss runner_slurm
for node in "${NODES[@]}"; do
    while IFS= read -r line; do
        read -r _s sweep_node dir rss slurm_ok <<< "$line"
        runner_node["$dir"]="$sweep_node"
        runner_rss["$dir"]="$rss"
        runner_slurm["$dir"]="$slurm_ok"
    done < <(grep '^RUNNER ' "$tmpdir/$node.out" 2>/dev/null || true)
done

# Fetch GitHub API status
declare -A gh_status gh_busy
while read -r _id name status busy; do
    gh_status["$name"]="$status"
    gh_busy["$name"]="$busy"
done < <(gh_list_runners)

# Print table
printf "%-25s %-8s %-20s %-8s %s\n" "NAME" "GITHUB" "NODE" "SLURM" "RSS"
printf "%s\n" "$(printf '%.0s-' {1..70})"

while IFS= read -r dir; do
    name=$(get_runner_name "$dir")
    [ -z "$name" ] && continue

    [ "${gh_busy[$name]:-false}" = "true" ] && gh_col="BUSY" || gh_col="${gh_status[$name]:-unknown}"

    actual_node="${runner_node[$dir]:-}"
    rss="${runner_rss[$dir]:-—}"
    slurm="${runner_slurm[$dir]:-—}"

    if [ -z "$actual_node" ]; then
        printf "%-25s %-8s %-20s %-8s %s\n" "$name" "$gh_col" "offline" "—" "—"
        continue
    fi

    # Flag stale runner.node entries
    node_col="$actual_node"
    if [ -f "$dir/runner.node" ]; then
        recorded=$(cat "$dir/runner.node")
        [ "$actual_node" != "$recorded" ] && node_col="${actual_node} *(stale: ${recorded})"
    fi

    printf "%-25s %-8s %-20s %-8s %sMB\n" "$name" "$gh_col" "$node_col" "$slurm" "$rss"
done < <(find_runner_dirs)

# Per-node memory summary (only when site has a cgroup limit)
if [ "${CGROUP_LIMIT:-0}" -gt 0 ]; then
    echo ""
    echo "=== Per-node memory ==="
    for node in "${NODES[@]}"; do
        count=$(ssh $SSH_OPTS "$node" \
            "ps aux | grep Runner.Listener | grep -v grep | wc -l" 2>/dev/null || echo 0)
        rss=$(ssh $SSH_OPTS "$node" \
            "ps -u \$(whoami) -o rss= 2>/dev/null | awk '{sum+=\$1} END {printf \"%.0f\", sum/1024}'" \
            2>/dev/null || echo "?")
        [[ "$rss" =~ ^[0-9]+$ ]] || rss=0
        echo "  $node: $count runners, ${rss} MB / ${CGROUP_LIMIT} MB ($(( CGROUP_LIMIT - rss )) MB free)"
    done
fi
