#!/bin/bash
# Restart all GitHub Actions runners on their current nodes.
#
# Useful after a login node reboot or when runners need a fresh start
# (e.g. to pick up a new PATH or clear stale state).  Each runner is
# restarted in place — no rebalancing is done.
#
# Usage: bash restart-all.sh              # dry run (show what would restart)
#        APPLY=1 bash restart-all.sh      # restart all runners
#        APPLY=1 FORCE=1 bash restart-all.sh  # restart even busy runners

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

discover_runners

echo "=== Discovering runners ==="
declare -a restart_list=()  # "node dir name"

for i in "${!RUNNER_DIRS[@]}"; do
    dir="${RUNNER_DIRS[$i]}"
    name="${RUNNER_NAMES[$i]}"
    node=$(find_runner_node "$dir")

    if [ "$node" = "offline" ]; then
        echo "  $name: OFFLINE (skipping — use rebalance-runners.sh to place it)"
        continue
    fi

    worker=$(ssh -o ConnectTimeout=5 "$node" "ps aux | grep Runner.Worker | grep '$dir' | grep -v grep" 2>/dev/null || true)
    if [ -n "$worker" ]; then
        echo "  $name: BUSY on $node"
        if [ "${FORCE:-0}" != "1" ]; then
            echo "    Skipping busy runner. Set FORCE=1 to restart anyway."
            continue
        fi
    else
        echo "  $name: idle on $node"
    fi

    restart_list+=("$node $dir $name")
done

if [ ${#restart_list[@]} -eq 0 ]; then
    echo ""
    echo "Nothing to restart."
    exit 0
fi

echo ""
echo "${#restart_list[@]} runners will be restarted."

if [ "${APPLY:-0}" != "1" ]; then
    echo "Dry run — set APPLY=1 to execute."
    exit 0
fi

echo ""
echo "=== Restarting ==="
success=0
fail=0
for entry in "${restart_list[@]}"; do
    read -r node dir name <<< "$entry"
    echo "--- $name on $node ---"
    stop_runner "$node" "$dir"
    if start_runner "$node" "$dir"; then
        pid=$(ssh -o ConnectTimeout=5 "$node" "pgrep -f 'Runner.Listener.*$dir' | head -1" 2>/dev/null || true)
        if [ -n "$pid" ] && check_slurm_path "$node" "$pid"; then
            echo "  OK: PID $pid, slurm in PATH"
            success=$((success + 1))
        elif [ -n "$pid" ]; then
            echo "  WARNING: PID $pid but slurm MISSING from PATH"
            fail=$((fail + 1))
        fi
    else
        echo "  ERROR: Failed to start"
        fail=$((fail + 1))
    fi
done

echo ""
echo "=== Summary: $success succeeded, $fail failed ==="
echo ""
bash "$SCRIPT_DIR/list-runners.sh"
