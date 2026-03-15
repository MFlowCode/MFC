#!/usr/bin/env bash
# Quick health check for GitHub Actions runners across Frontier login nodes.
#
# SSHes to each login node, finds Runner.Listener processes, and shows
# runner name, status (idle/BUSY), and memory usage.
#
# Usage: bash check-runners.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

for node in "${NODES[@]}"; do
    echo "=== $node ==="
    ssh $SSH_OPTS "$node" '
        found=0
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            found=1
            cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || echo "???")
            worker=$(ps aux | grep "Runner.Worker" | grep "$cwd" | grep -v grep | awk "{print \$2}" | head -1)
            [ -n "$worker" ] && status="BUSY" || status="idle"
            rss=$(ps -p $p -o rss= 2>/dev/null | awk "{printf \"%.0f\", \$1/1024}" || echo "?")
            name=$(basename "$cwd")
            printf "  %-30s %5s  %s MB\n" "$name" "$status" "$rss"
        done
        [ "$found" -eq 0 ] && echo "  (no runners)"
    ' 2>/dev/null || echo "  (unreachable)"
    echo ""
done
