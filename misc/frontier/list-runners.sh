#!/usr/bin/env bash
# List all Frontier runners, combining GitHub API status with login-node process info.
#
# Uses a parallel SSH sweep across all 11 login nodes simultaneously to avoid
# the overhead of serial per-runner node discovery. Each node is queried once;
# results are correlated with GitHub API status and the local runner directories.
#
# Usage: bash list-runners.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

echo "==> Syncing runner node locations..."
sync_runner_nodes

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# --- Parallel SSH sweep across all login nodes ---
# Each node prints lines in the format: RUNNER <node> <dir> <rss_mb>
# The RUNNER sentinel prefix allows stripping MOTD noise with grep.
for node in "${NODES[@]}"; do
    ssh $SSH_OPTS "$node" '
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || true)
            rss=$(ps -p $p -o rss= 2>/dev/null | awk "{printf \"%.0f\", \$1/1024}" || echo 0)
            [ -n "$cwd" ] && echo "RUNNER '"$node"' $cwd $rss"
        done
    ' 2>/dev/null > "$tmpdir/$node.out" &
done

wait

# --- Build associative arrays from sweep results ---
declare -A runner_node runner_rss
for node in "${NODES[@]}"; do
    while IFS= read -r line; do
        # Each line: RUNNER <node> <dir> <rss_mb>
        read -r _sentinel sweep_node dir rss <<< "$line"
        runner_node["$dir"]="$sweep_node"
        runner_rss["$dir"]="$rss"
    done < <(grep '^RUNNER ' "$tmpdir/$node.out" 2>/dev/null || true)
done

# --- Fetch GitHub API status ---
declare -A gh_status gh_busy
while read -r _id name status busy; do
    gh_status["$name"]="$status"
    gh_busy["$name"]="$busy"
done < <(gh_list_runners)

# --- Print table ---
printf "%-25s %-8s %-14s %s\n" "NAME" "GITHUB" "NODE" "RSS"
printf "%s\n" "$(printf '%.0s-' {1..60})"

for dir in $(find_runner_dirs); do
    name=$(get_runner_name "$dir")
    [ -z "$name" ] && continue

    # GitHub status column
    api_status="${gh_status[$name]:-unknown}"
    api_busy="${gh_busy[$name]:-false}"
    if [ "$api_busy" = "true" ]; then
        gh_col="BUSY"
    else
        gh_col="$api_status"
    fi

    # Node and RSS from parallel sweep
    actual_node="${runner_node[$dir]:-}"
    rss="${runner_rss[$dir]:-}"

    if [ -z "$actual_node" ]; then
        printf "%-25s %-8s %-14s %s\n" "$name" "$gh_col" "offline" "—"
        continue
    fi

    # Compare sweep result to recorded runner.node; flag stale entries
    node_col="$actual_node"
    if [ -f "$dir/runner.node" ]; then
        recorded=$(cat "$dir/runner.node")
        if [ "$actual_node" != "$recorded" ]; then
            node_col="${actual_node} *(stale: ${recorded})"
        fi
    fi

    printf "%-25s %-8s %-14s %sMB\n" "$name" "$gh_col" "$node_col" "$rss"
done
