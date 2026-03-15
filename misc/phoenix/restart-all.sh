#!/bin/bash
# Restart all Phoenix runners on their current nodes.
#
# Useful after a login node reboot or to pick up environment changes.
# Restarts in place — no rebalancing. Skips busy runners unless FORCE=1.
#
# Usage: bash restart-all.sh              # dry run
#        APPLY=1 bash restart-all.sh      # execute
#        APPLY=1 FORCE=1 bash restart-all.sh  # restart busy runners too

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

echo "=== Discovering runners ==="
declare -a restart_list=()

while IFS= read -r dir; do
    name=$(get_runner_name "$dir")
    [ -z "$name" ] && continue
    node=$(find_node "$dir")

    if [ "$node" = "offline" ]; then
        echo "  $name: OFFLINE (use rebalance-runners.sh to place)"
        continue
    fi

    worker=$(ssh $SSH_OPTS "$node" "ps aux | grep Runner.Worker | grep '$dir' | grep -v grep" 2>/dev/null || true)
    if [ -n "$worker" ]; then
        echo "  $name: BUSY on $node"
        if [ "${FORCE:-0}" != "1" ]; then
            echo "    Skipping. Set FORCE=1 to restart anyway."
            continue
        fi
    else
        echo "  $name: idle on $node"
    fi

    restart_list+=("$node $dir $name")
done < <(find_runner_dirs)

if [ ${#restart_list[@]} -eq 0 ]; then
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
success=0; fail=0
for entry in "${restart_list[@]}"; do
    read -r node dir name <<< "$entry"
    echo "--- $name on $node ---"
    stop_runner "$node" "$dir"
    if start_runner "$node" "$dir"; then
        pids=$(find_pids "$node" "$dir")
        pid=${pids%% *}
        if has_slurm "$node" "$pid"; then
            echo "  OK: PID $pid, slurm in PATH"
            success=$((success + 1))
        else
            echo "  WARNING: PID $pid but slurm MISSING"
            fail=$((fail + 1))
        fi
    else
        echo "  ERROR: Failed to start"
        fail=$((fail + 1))
    fi
done

echo ""
echo "=== Summary: $success succeeded, $fail failed ==="
