#!/bin/bash
# List all Phoenix runners, combining GitHub API status with login-node process info.
#
# Shows both what GitHub thinks (online/offline/busy) and the actual process
# state on the login nodes (which node, slurm PATH, memory).
#
# Usage: bash list-runners.sh

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

printf "%-25s %-8s %-22s %-8s %6s  %s\n" \
    "NAME" "GITHUB" "NODE" "SLURM" "RSS" "DIRECTORY"
printf "%s\n" "$(printf '%.0s-' {1..100})"

# Get GitHub API status for all Phoenix runners
declare -A gh_status gh_busy
while read -r id name status busy; do
    gh_status[$name]="$status"
    gh_busy[$name]="$busy"
done <<< "$(gh_list_runners)"

# Walk local runner directories and cross-reference
for dir in $(find_runner_dirs); do
    name=$(get_runner_name "$dir")
    [ -z "$name" ] && continue

    # GitHub status
    api_status="${gh_status[$name]:-unknown}"
    api_busy="${gh_busy[$name]:-false}"
    if [ "$api_busy" = "true" ]; then
        gh_col="BUSY"
    else
        gh_col="$api_status"
    fi

    # Node status
    node=$(find_node "$dir")
    if [ "$node" = "offline" ]; then
        printf "%-25s %-8s %-22s %-8s %6s  %s\n" \
            "$name" "$gh_col" "—" "—" "—" "$dir"
        continue
    fi

    # Process details (one SSH call)
    info=$(ssh -o ConnectTimeout=5 "$node" '
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || true)
            if [ "$cwd" = "'"$dir"'" ]; then
                slurm=$(cat /proc/$p/environ 2>/dev/null | tr "\0" "\n" | grep -c /opt/slurm || echo 0)
                [ "$slurm" -gt 0 ] && s="ok" || s="MISSING"
                rss=$(ps -p $p -o rss= 2>/dev/null | awk "{printf \"%.0f\", \$1/1024}" || echo "?")
                echo "$s $rss"
                exit
            fi
        done
        echo "? ?"
    ' 2>/dev/null || echo "? ?")
    read -r slurm rss <<< "$info"

    printf "%-25s %-8s %-22s %-8s %5sMB  %s\n" \
        "$name" "$gh_col" "$node" "$slurm" "$rss" "$dir"
done

echo ""
echo "=== Per-node memory ==="
for node in "${NODES[@]}"; do
    count=$(ssh -o ConnectTimeout=5 "$node" "ps aux | grep Runner.Listener | grep -v grep | wc -l" 2>/dev/null || echo 0)
    rss=$(ssh -o ConnectTimeout=5 "$node" "ps -u \$(whoami) -o rss= 2>/dev/null | awk '{sum+=\$1} END {printf \"%.0f\", sum/1024}'" 2>/dev/null || echo "?")
    echo "  $node: $count runners, ${rss} MB / ${CGROUP_LIMIT} MB ($(( CGROUP_LIMIT - rss )) MB free)"
done
