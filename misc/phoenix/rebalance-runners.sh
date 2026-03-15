#!/bin/bash
# Automatically rebalance GitHub Actions runners across Phoenix login nodes.
#
# Computes the optimal distribution (equal runners per node), determines
# which runners need to move, and executes the moves.  Prefers to move
# idle runners over busy ones.
#
# Each Phoenix login node has a 4 GB per-user cgroup memory limit.
# Target: ~3-4 runners per node to leave headroom for CI work.
#
# Usage: bash rebalance-runners.sh              # dry run
#        APPLY=1 bash rebalance-runners.sh      # execute moves
#        APPLY=1 FORCE=1 bash rebalance-runners.sh  # move busy runners too

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

discover_runners
num_nodes=${#NODES[@]}
num_runners=${#RUNNER_DIRS[@]}
target_per_node=$(( num_runners / num_nodes ))
remainder=$(( num_runners % num_nodes ))

echo "=== Current state ==="
echo "Runners: $num_runners across $num_nodes nodes"
echo "Target:  $target_per_node per node (+1 on first $remainder nodes)"
echo ""

# Build current assignment: node -> list of runner indices
declare -A node_runners  # node -> space-separated indices
declare -A runner_node   # index -> node
declare -A runner_busy   # index -> 1 if busy

for node in "${NODES[@]}"; do
    node_runners[$node]=""
done

for i in "${!RUNNER_DIRS[@]}"; do
    dir="${RUNNER_DIRS[$i]}"
    node=$(find_runner_node "$dir")
    runner_node[$i]="$node"

    if [ "$node" != "offline" ]; then
        node_runners[$node]="${node_runners[$node]:-} $i"
        worker=$(ssh -o ConnectTimeout=5 "$node" "ps aux | grep Runner.Worker | grep '$dir' | grep -v grep" 2>/dev/null || true)
        [ -n "$worker" ] && runner_busy[$i]=1 || runner_busy[$i]=0
    else
        runner_busy[$i]=0
    fi
done

# Show current distribution
for node in "${NODES[@]}"; do
    indices=(${node_runners[$node]:-})
    echo "$node: ${#indices[@]} runners"
    for i in "${indices[@]}"; do
        busy_marker=""
        [ "${runner_busy[$i]:-0}" = "1" ] && busy_marker=" (BUSY)"
        echo "  ${RUNNER_NAMES[$i]}$busy_marker"
    done
done

# Find offline runners
offline=()
for i in "${!RUNNER_DIRS[@]}"; do
    [ "${runner_node[$i]}" = "offline" ] && offline+=("$i")
done
if [ ${#offline[@]} -gt 0 ]; then
    echo ""
    echo "OFFLINE runners:"
    for i in "${offline[@]}"; do
        echo "  ${RUNNER_NAMES[$i]} (${RUNNER_DIRS[$i]})"
    done
fi

echo ""

# Compute target per node
declare -A node_target
n=0
for node in "${NODES[@]}"; do
    node_target[$node]=$target_per_node
    if [ $n -lt $remainder ]; then
        node_target[$node]=$(( target_per_node + 1 ))
    fi
    n=$((n + 1))
done

# Determine moves needed
# Phase 1: identify overloaded nodes and runners to move away
moves=()  # "source_node dest_node runner_index"

# Collect runners to move from overloaded nodes (prefer idle runners)
to_place=()  # indices of runners that need a new home
for node in "${NODES[@]}"; do
    indices=(${node_runners[$node]:-})
    excess=$(( ${#indices[@]} - ${node_target[$node]} ))
    if [ $excess -le 0 ]; then
        continue
    fi
    # Sort: move idle runners first, then busy
    idle_here=()
    busy_here=()
    for i in "${indices[@]}"; do
        if [ "${runner_busy[$i]:-0}" = "1" ]; then
            busy_here+=("$i")
        else
            idle_here+=("$i")
        fi
    done
    # Pick from idle first
    moved=0
    for i in "${idle_here[@]}"; do
        [ $moved -ge $excess ] && break
        to_place+=("$node $i")
        moved=$((moved + 1))
    done
    for i in "${busy_here[@]}"; do
        [ $moved -ge $excess ] && break
        to_place+=("$node $i")
        moved=$((moved + 1))
    done
done

# Phase 2: assign runners to underloaded nodes
# Also place offline runners
for i in "${offline[@]}"; do
    to_place+=("offline $i")
done

for entry in "${to_place[@]}"; do
    read -r src_node runner_idx <<< "$entry"
    # Find the most underloaded node
    best_node=""
    best_deficit=-999
    for node in "${NODES[@]}"; do
        current=(${node_runners[$node]:-})
        deficit=$(( ${node_target[$node]} - ${#current[@]} ))
        if [ $deficit -gt $best_deficit ]; then
            best_deficit=$deficit
            best_node=$node
        fi
    done

    if [ -z "$best_node" ] || [ "$best_deficit" -le 0 ]; then
        echo "WARNING: No underloaded node for ${RUNNER_NAMES[$runner_idx]}, skipping"
        continue
    fi

    moves+=("$src_node $best_node $runner_idx")
    # Update bookkeeping
    if [ "$src_node" != "offline" ]; then
        # Remove from source
        new_list=""
        for idx in ${node_runners[$src_node]}; do
            [ "$idx" != "$runner_idx" ] && new_list="$new_list $idx"
        done
        node_runners[$src_node]="$new_list"
    fi
    node_runners[$best_node]="${node_runners[$best_node]:-} $runner_idx"
done

if [ ${#moves[@]} -eq 0 ]; then
    echo "Already balanced — no moves needed."
    exit 0
fi

# Show plan
echo "=== Planned moves ==="
has_busy=false
for move in "${moves[@]}"; do
    read -r src dst idx <<< "$move"
    busy_marker=""
    if [ "${runner_busy[$idx]:-0}" = "1" ]; then
        busy_marker=" (BUSY!)"
        has_busy=true
    fi
    echo "  ${RUNNER_NAMES[$idx]}: $src -> $dst$busy_marker"
done

echo ""
echo "=== Target distribution ==="
for node in "${NODES[@]}"; do
    indices=(${node_runners[$node]:-})
    echo "  $node: ${#indices[@]} runners"
done

if [ "$has_busy" = true ] && [ "${FORCE:-0}" != "1" ]; then
    echo ""
    echo "Some runners to move have active jobs. Set FORCE=1 to move them."
    exit 1
fi

if [ "${APPLY:-0}" != "1" ]; then
    echo ""
    echo "Dry run — set APPLY=1 to execute. Add FORCE=1 to move busy runners."
    exit 0
fi

# Execute
echo ""
echo "=== Executing moves ==="
for move in "${moves[@]}"; do
    read -r src dst idx <<< "$move"
    dir="${RUNNER_DIRS[$idx]}"
    name="${RUNNER_NAMES[$idx]}"
    echo "Moving $name: $src -> $dst"

    if [ "$src" != "offline" ]; then
        stop_runner "$src" "$dir"
    fi

    if start_runner "$dst" "$dir"; then
        pid=$(ssh -o ConnectTimeout=5 "$dst" "pgrep -f 'Runner.Listener.*$dir' | head -1" 2>/dev/null || true)
        if [ -n "$pid" ] && check_slurm_path "$dst" "$pid"; then
            echo "  OK: PID $pid on $dst, slurm in PATH"
        elif [ -n "$pid" ]; then
            echo "  WARNING: PID $pid on $dst but slurm MISSING from PATH"
        else
            echo "  ERROR: Process not found after start"
        fi
    else
        echo "  ERROR: Failed to start on $dst"
    fi
done

echo ""
echo "=== Final state ==="
bash "$SCRIPT_DIR/list-runners.sh"
