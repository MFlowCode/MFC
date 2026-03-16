#!/usr/bin/env bash
# Restart all runners in place on their current nodes.
#
# Sourced by site wrappers (frontier/restart-all.sh, phoenix/restart-all.sh)
# after config.sh is loaded. Useful after a login node reboot or to pick up
# environment changes. Skips busy runners unless FORCE=1. Dry run by default.
#
# Usage: bash restart-all.sh              # dry run
#        APPLY=1 bash restart-all.sh      # execute
#        APPLY=1 FORCE=1 bash restart-all.sh  # restart busy runners too
set -euo pipefail

declare -f sync_runner_nodes > /dev/null 2>&1 && {
    echo "==> Syncing runner node locations..."
    sync_runner_nodes
}

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
        echo "$node" > "$dir/runner.node"
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
