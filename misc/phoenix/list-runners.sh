#!/bin/bash
# List all registered GitHub Actions runners, showing which node each is on,
# whether it's busy, its memory usage, and whether slurm is in PATH.
#
# Output columns:
#   Name       — GitHub runner name (from .runner config)
#   Node       — login node it's running on, or "offline"
#   Status     — idle, BUSY (has Worker process), or OFFLINE
#   Slurm      — ok or MISSING (can't submit SLURM jobs)
#   RSS        — memory usage in MB
#   Pool       — GitHub runner group/pool
#   Directory  — filesystem path
#
# Usage: bash list-runners.sh

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

discover_runners

# Header
printf "%-25s %-22s %-8s %-8s %5s  %-10s %s\n" \
    "NAME" "NODE" "STATUS" "SLURM" "RSS" "POOL" "DIRECTORY"
printf "%s\n" "$(printf '%.0s-' {1..120})"

for i in "${!RUNNER_DIRS[@]}"; do
    dir="${RUNNER_DIRS[$i]}"
    name="${RUNNER_NAMES[$i]}"
    pool="${RUNNER_POOLS[$i]}"

    node=$(find_runner_node "$dir")

    if [ "$node" = "offline" ]; then
        printf "%-25s %-22s %-8s %-8s %5s  %-10s %s\n" \
            "$name" "—" "OFFLINE" "—" "—" "$pool" "$dir"
        continue
    fi

    # Get all info in one SSH call to reduce latency
    info=$(ssh -o ConnectTimeout=5 "$node" '
        pid=""
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || true)
            [ "$cwd" = "'"$dir"'" ] && pid=$p && break
        done
        [ -z "$pid" ] && echo "? ? ? ?" && exit
        worker=$(ps aux | grep Runner.Worker | grep "'"$dir"'" | grep -v grep | head -1 || true)
        [ -n "$worker" ] && status="BUSY" || status="idle"
        has_slurm=$(cat /proc/$pid/environ 2>/dev/null | tr "\0" "\n" | grep -c /opt/slurm || echo 0)
        [ "$has_slurm" -gt 0 ] && slurm="ok" || slurm="MISSING"
        rss=$(ps -p $pid -o rss= 2>/dev/null | awk "{printf \"%.0f\", \$1/1024}" || echo "?")
        echo "$pid $status $slurm $rss"
    ' 2>/dev/null || echo "? ? ? ?")
    read -r pid status slurm rss <<< "$info"

    printf "%-25s %-22s %-8s %-8s %5s  %-10s %s\n" \
        "$name" "$node" "$status" "$slurm" "${rss}MB" "$pool" "$dir"
done

# Per-node summary
echo ""
echo "=== Per-node memory ==="
for node in "${NODES[@]}"; do
    count=$(ssh -o ConnectTimeout=5 "$node" "ps aux | grep Runner.Listener | grep -v grep | wc -l" 2>/dev/null || echo 0)
    rss=$(ssh -o ConnectTimeout=5 "$node" "ps -u \$(whoami) -o rss= 2>/dev/null | awk '{sum+=\$1} END {printf \"%.0f\", sum/1024}'" 2>/dev/null || echo "?")
    echo "  $node: $count runners, ${rss} MB / ${CGROUP_LIMIT} MB ($(( CGROUP_LIMIT - rss )) MB free)"
done
