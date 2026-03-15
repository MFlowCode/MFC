#!/bin/bash
# Quick health check for GitHub Actions runners across Phoenix login nodes.
#
# Lighter than list-runners.sh — doesn't query each runner individually,
# just shows per-node counts and memory.  Use list-runners.sh for details.
#
# Usage: bash check-runners.sh

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

for node in "${NODES[@]}"; do
    echo "=== $node ==="
    ssh -o ConnectTimeout=5 "$node" '
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || echo "???")
            has_slurm=$(cat /proc/$p/environ 2>/dev/null | tr "\0" "\n" | grep -c /opt/slurm || echo 0)
            worker=$(ps aux | grep "Runner.Worker" | grep "$cwd" | grep -v grep | awk "{print \$2}" | head -1)
            [ -n "$worker" ] && status="BUSY" || status="idle"
            rss=$(ps -p $p -o rss= 2>/dev/null | awk "{printf \"%.0f\", \$1/1024}" || echo "?")
            name=$(basename "$cwd")
            parent=$(basename $(dirname "$cwd"))
            slurm_ok="ok"
            [ "$has_slurm" -eq 0 ] && slurm_ok="MISSING"
            printf "  %-40s %5s  slurm=%-7s  %s MB\n" "$parent/$name" "$status" "$slurm_ok" "$rss"
        done
    ' 2>/dev/null || echo "  (unreachable)"

    rss=$(ssh -o ConnectTimeout=5 "$node" "ps -u \$(whoami) -o rss= 2>/dev/null | awk '{sum+=\$1} END {printf \"%.0f\", sum/1024}'" 2>/dev/null || echo "?")
    [[ "$rss" =~ ^[0-9]+$ ]] || rss=0
    echo "  --- Total: ${rss} MB / ${CGROUP_LIMIT} MB ($(( CGROUP_LIMIT - rss )) MB free) ---"
    echo ""
done
