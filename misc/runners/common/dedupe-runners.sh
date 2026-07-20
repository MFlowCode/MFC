#!/usr/bin/env bash
# Detect and remove duplicate GitHub Actions runner listeners.
#
# Sourced by site wrappers (frontier/dedupe-runners.sh via the dispatcher)
# after config.sh is loaded.
#
# Because runner binaries live on shared storage, a botched move/restart can
# leave the SAME runner listening on two login nodes at once. The two processes
# then fight over the single GitHub registration and share one _diag/ directory,
# which cancels in-progress jobs and produces "file already exists" errors in
# job setup. This sweeps every node, groups listeners by runner directory, and
# for any runner found on more than one node keeps a single listener (preferring
# the busy one, then the node recorded in runner.node) and stops the rest.
#
# Usage: bash dedupe-runners.sh              # dry run
#        APPLY=1 bash dedupe-runners.sh      # execute
set -euo pipefail

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

echo "=== Sweeping all nodes for runner listeners ==="
mapfile -t lines < <(sweep_all_listeners "$tmpdir")

# Group listeners by directory: dir -> "node:pid:busy node:pid:busy ..."
declare -A locs
for line in "${lines[@]}"; do
    read -r _tag node pid dir busy <<< "$line"
    locs["$dir"]="${locs["$dir"]:-} $node:$pid:$busy"
done

dupes=0
declare -a kill_nodes=() kill_pids=() kill_names=()
for dir in "${!locs[@]}"; do
    entries=(${locs["$dir"]})
    [ "${#entries[@]}" -le 1 ] && continue
    dupes=$((dupes + 1))
    name=$(get_runner_name "$dir" 2>/dev/null || basename "$dir")

    recorded=""
    [ -f "$dir/runner.node" ] && recorded=$(cat "$dir/runner.node" 2>/dev/null || echo "")

    echo ""
    echo "DUPLICATE: $name on ${#entries[@]} nodes"
    for e in "${entries[@]}"; do
        IFS=: read -r n p b <<< "$e"
        tag=""
        [ "$n" = "$recorded" ] && tag=" (runner.node)"
        echo "  $n pid=$p $b$tag"
    done

    # Choose the keeper: first busy listener (never kill a running job if we
    # can avoid it), else the one matching runner.node, else the first.
    keeper=""
    for e in "${entries[@]}"; do
        IFS=: read -r n p b <<< "$e"
        [ "$b" = "busy" ] && keeper="$e" && break
    done
    if [ -z "$keeper" ] && [ -n "$recorded" ]; then
        for e in "${entries[@]}"; do
            IFS=: read -r n p b <<< "$e"
            [ "$n" = "$recorded" ] && keeper="$e" && break
        done
    fi
    [ -z "$keeper" ] && keeper="${entries[0]}"

    IFS=: read -r kn kp kb <<< "$keeper"
    echo "  -> keep $kn pid=$kp ($kb); stop the rest"
    for e in "${entries[@]}"; do
        [ "$e" = "$keeper" ] && continue
        IFS=: read -r n p _b <<< "$e"
        kill_nodes+=("$n"); kill_pids+=("$p"); kill_names+=("$name")
    done
done

if [ "$dupes" -eq 0 ]; then
    echo ""
    echo "No duplicate runners found."
    exit 0
fi

echo ""
echo "${#kill_pids[@]} extra listener(s) to stop across $dupes runner(s)."
if [ "${APPLY:-0}" != "1" ]; then
    echo "Dry run — set APPLY=1 to execute."
    exit 0
fi

echo ""
echo "=== Stopping extra listeners ==="
for i in "${!kill_pids[@]}"; do
    n="${kill_nodes[$i]}"; p="${kill_pids[$i]}"; nm="${kill_names[$i]}"
    echo "Stopping duplicate $nm: $n pid=$p"
    ssh $SSH_OPTS "$n" "kill $p" 2>/dev/null || true
    sleep 2
    if ssh $SSH_OPTS "$n" "kill -0 $p" 2>/dev/null; then
        ssh $SSH_OPTS "$n" "kill -9 $p" 2>/dev/null || true
    fi
done

echo ""
echo "Done. Re-run without APPLY to confirm no duplicates remain."
