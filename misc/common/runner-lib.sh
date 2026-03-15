#!/usr/bin/env bash
# Shared GitHub Actions runner management library.
#
# Sourced by site-specific config.sh files (misc/frontier/config.sh,
# misc/phoenix/config.sh). Callers must define ORG, NODES, and SSH_OPTS
# before sourcing this file.

# --- GitHub API ---

# Get a registration token for new runners.
gh_registration_token() {
    gh api "orgs/$ORG/actions/runners/registration-token" --jq .token
}

# Get the latest runner binary version.
gh_latest_runner_version() {
    gh api repos/actions/runner/releases/latest --jq '.tag_name | ltrimstr("v")'
}

# Remove a runner registration from GitHub.
# Args: $1 = runner ID (numeric, from API)
gh_remove_runner() {
    gh api "orgs/$ORG/actions/runners/$1" -X DELETE
}

# --- Local filesystem ---

# Get the GitHub runner name from a .runner config file.
# Args: $1 = runner directory
get_runner_name() {
    python3 -c "
import json
d = json.loads(open('$1/.runner').read().lstrip('\ufeff'))
print(d.get('agentName', ''))
" 2>/dev/null
}

# --- Login-node process management ---

# Find PIDs of a runner on a node by matching its executable path.
# Matches /proc/$p/exe against $dir/bin/Runner.Listener — intrinsic to
# the binary, independent of CWD or how the process was launched.
# Output is filtered to numeric lines only to strip SSH MOTD noise.
# Args: $1 = node, $2 = runner directory
# Prints: space-separated PIDs, or empty.
find_pids() {
    ssh $SSH_OPTS "$1" '
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            exe=$(readlink -f /proc/$p/exe 2>/dev/null || true)
            [ "$exe" = "'"$2"'/bin/Runner.Listener" ] && echo "$p"
        done
    ' 2>/dev/null | grep -E '^[0-9]+$' | tr '\n' ' ' || true
}

# Find which login node a runner is on.
# Args: $1 = runner directory
# Prints: node hostname, or "offline".
find_node() {
    for node in "${NODES[@]}"; do
        [ -n "$(find_pids "$node" "$1")" ] && echo "$node" && return
    done
    echo "offline"
}

# Start a runner on a node.
# Uses a login shell (bash -lc) so site PATH (e.g. SLURM) is available.
# Args: $1 = node, $2 = runner directory
# Returns: 0 if running after start, 1 otherwise.
start_runner() {
    local node="$1" dir="$2"
    timeout 15 ssh $SSH_OPTS "$node" \
        "cd $dir && setsid bash -lc 'nohup ./run.sh >> runner.log 2>&1 < /dev/null &'" \
        </dev/null 2>/dev/null || true
    sleep 3
    [ -n "$(find_pids "$node" "$dir")" ]
}

# Stop a runner on a node (SIGTERM then SIGKILL).
# Args: $1 = node, $2 = runner directory
stop_runner() {
    local node="$1" dir="$2" pids
    pids=$(find_pids "$node" "$dir")
    [ -z "$pids" ] && return 0
    for pid in $pids; do
        ssh $SSH_OPTS "$node" "kill $pid" 2>/dev/null || true
    done
    sleep 3
    pids=$(find_pids "$node" "$dir")
    for pid in $pids; do
        ssh $SSH_OPTS "$node" "kill -9 $pid" 2>/dev/null || true
    done
    sleep 1
}
