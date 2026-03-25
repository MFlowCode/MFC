#!/usr/bin/env bash
# Check runner health across all login nodes.
#
# Sourced by site wrappers (frontier/check-runners.sh, phoenix/check-runners.sh)
# after config.sh is loaded. Shows Runner.Listener processes per node with
# name, busy/idle status, slurm availability, and RSS memory.
# If CGROUP_LIMIT > 0, also shows per-node total memory vs the cgroup limit.
#
# Usage: bash check-runners.sh
set -euo pipefail

declare -f sync_runner_nodes > /dev/null 2>&1 && {
    echo "==> Syncing runner node locations..."
    sync_runner_nodes
}

for node in "${NODES[@]}"; do
    echo "=== $node ==="
    ssh $SSH_OPTS "$node" '
        found=0
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            found=1
            exe=$(readlink -f /proc/$p/exe 2>/dev/null || echo "???")
            dir=$(dirname "$(dirname "$exe")" 2>/dev/null || echo "???")
            name=$(basename "$dir")
            worker=$(ps aux | grep "Runner.Worker" | grep "$dir" | grep -v grep | awk "{print \$2}" | head -1)
            [ -n "$worker" ] && status="BUSY" || status="idle"
            rss=$(ps -p $p -o rss= 2>/dev/null | awk "{printf \"%.0f\", \$1/1024}" || echo "?")
            slurm=$(tr "\0" "\n" < /proc/$p/environ 2>/dev/null | grep -c "^PATH=.*slurm" || echo 0)
            [ "$slurm" -gt 0 ] && slurm_ok="ok" || slurm_ok="MISSING"
            printf "  %-30s %5s  slurm=%-7s  %s MB\n" "$name" "$status" "$slurm_ok" "$rss"
        done
        [ "$found" -eq 0 ] && echo "  (no runners)"
    ' 2>/dev/null || echo "  (unreachable)"

    if [ "${CGROUP_LIMIT:-0}" -gt 0 ]; then
        rss=$(ssh $SSH_OPTS "$node" \
            "ps -u \$(whoami) -o rss= 2>/dev/null | awk '{sum+=\$1} END {printf \"%.0f\", sum/1024}'" \
            2>/dev/null || echo "?")
        [[ "$rss" =~ ^[0-9]+$ ]] || rss=0
        echo "  --- Total: ${rss} MB / ${CGROUP_LIMIT} MB ($(( CGROUP_LIMIT - rss )) MB free) ---"
    fi
    echo ""
done
