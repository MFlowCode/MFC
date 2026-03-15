#!/usr/bin/env bash
# Automatically rebalance Frontier runners across login nodes.
#
# Discovers all runner directories, checks which node each is on,
# computes the optimal distribution, and moves runners to balance.
# Prefers moving idle runners over busy ones. Also places offline runners.
#
# Usage: bash rebalance-runners.sh              # dry run
#        APPLY=1 bash rebalance-runners.sh      # execute
#        APPLY=1 FORCE=1 bash rebalance-runners.sh  # move busy runners too
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

echo "==> Syncing runner node locations..."
sync_runner_nodes

# Discover runners
declare -a dirs=() names=()
while IFS= read -r dir; do
    name=$(get_runner_name "$dir")
    [ -z "$name" ] && continue
    dirs+=("$dir")
    names+=("$name")
done < <(find_runner_dirs)

num_nodes=${#NODES[@]}
num_runners=${#dirs[@]}
target=$(( num_runners / num_nodes ))
remainder=$(( num_runners % num_nodes ))

echo "=== Current state ==="
echo "Runners: $num_runners across $num_nodes nodes"
echo "Target:  $target per node (+1 on first $remainder nodes)"
echo ""

# Map runners to nodes and check busy status
declare -A node_runners
declare -A runner_node runner_busy

for node in "${NODES[@]}"; do node_runners[$node]=""; done

for i in "${!dirs[@]}"; do
    node=$(find_node "${dirs[$i]}")
    runner_node[$i]="$node"
    if [ "$node" != "offline" ]; then
        node_runners[$node]="${node_runners[$node]:-} $i"
        worker=$(ssh $SSH_OPTS "$node" "ps aux | grep Runner.Worker | grep '${dirs[$i]}' | grep -v grep" 2>/dev/null || true)
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
        busy=""
        [ "${runner_busy[$i]:-0}" = "1" ] && busy=" (BUSY)"
        echo "  ${names[$i]}$busy"
    done
done

offline=()
for i in "${!dirs[@]}"; do
    [ "${runner_node[$i]}" = "offline" ] && offline+=("$i")
done
if [ ${#offline[@]} -gt 0 ]; then
    echo ""
    echo "OFFLINE:"
    for i in "${offline[@]}"; do echo "  ${names[$i]}"; done
fi
echo ""

# Compute per-node targets
declare -A node_target
n=0
for node in "${NODES[@]}"; do
    node_target[$node]=$target
    [ $n -lt $remainder ] && node_target[$node]=$(( target + 1 ))
    n=$((n + 1))
done

# Plan moves: pull runners from overloaded nodes (idle first)
to_place=()
for node in "${NODES[@]}"; do
    indices=(${node_runners[$node]:-})
    excess=$(( ${#indices[@]} - ${node_target[$node]} ))
    [ $excess -le 0 ] && continue
    idle=() busy=()
    for i in "${indices[@]}"; do
        [ "${runner_busy[$i]:-0}" = "1" ] && busy+=("$i") || idle+=("$i")
    done
    moved=0
    for i in "${idle[@]}" "${busy[@]}"; do
        [ $moved -ge $excess ] && break
        to_place+=("$node $i")
        moved=$((moved + 1))
    done
done

# Add offline runners to be placed
for i in "${offline[@]}"; do to_place+=("offline $i"); done

# Assign to underloaded nodes
moves=()
for entry in "${to_place[@]}"; do
    read -r src idx <<< "$entry"
    best="" best_deficit=-999
    for node in "${NODES[@]}"; do
        cur=(${node_runners[$node]:-})
        deficit=$(( ${node_target[$node]} - ${#cur[@]} ))
        [ $deficit -gt $best_deficit ] && best_deficit=$deficit && best=$node
    done
    [ -z "$best" ] || [ "$best_deficit" -le 0 ] && continue
    moves+=("$src $best $idx")
    # Update bookkeeping so subsequent assignments reflect this move
    if [ "$src" != "offline" ]; then
        new=""
        for j in ${node_runners[$src]}; do [ "$j" != "$idx" ] && new="$new $j"; done
        node_runners[$src]="$new"
    fi
    node_runners[$best]="${node_runners[$best]:-} $idx"
done

if [ ${#moves[@]} -eq 0 ]; then
    echo "Already balanced."
    exit 0
fi

echo "=== Planned moves ==="
has_busy=false
for move in "${moves[@]}"; do
    read -r src dst idx <<< "$move"
    busy=""
    [ "${runner_busy[$idx]:-0}" = "1" ] && busy=" (BUSY!)" && has_busy=true
    echo "  ${names[$idx]}: $src -> $dst$busy"
done
echo ""
echo "=== Target ==="
for node in "${NODES[@]}"; do
    cur=(${node_runners[$node]:-})
    echo "  $node: ${#cur[@]} runners"
done

[ "$has_busy" = true ] && [ "${FORCE:-0}" != "1" ] && echo "" && echo "Set FORCE=1 to move busy runners." && exit 1
[ "${APPLY:-0}" != "1" ] && echo "" && echo "Dry run — set APPLY=1 to execute." && exit 0

echo ""
echo "=== Executing ==="
for move in "${moves[@]}"; do
    read -r src dst idx <<< "$move"
    echo "Moving ${names[$idx]}: $src -> $dst"
    [ "$src" != "offline" ] && stop_runner "$src" "${dirs[$idx]}"
    if start_runner "$dst" "${dirs[$idx]}"; then
        echo "$dst" > "${dirs[$idx]}/runner.node"
        echo "  OK: ${names[$idx]} started on $dst"
    else
        echo "  ERROR: Failed to start ${names[$idx]} on $dst"
    fi
done

echo ""
bash "$SCRIPT_DIR/check-runners.sh"
