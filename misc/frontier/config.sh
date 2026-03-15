#!/usr/bin/env bash
# Shared configuration for Frontier GitHub Actions runner management.
#
# Sourced by all other scripts. Provides constants, GitHub API helpers,
# and login-node process management functions.

# --- Frontier constants ---
ORG="MFlowCode"
RUNNER_GROUP="phoenix"
RUNNER_LABEL="frontier"
NODES=(login01 login02 login03 login04 login05 login06 login07 login08 login09 login10 login11)
SHARED_DIR="/lustre/orion/cfd154/proj-shared/runners"

SSH_OPTS="-o StrictHostKeyChecking=no -o ConnectTimeout=10 -o BatchMode=yes -o ServerAliveInterval=10 -o ServerAliveCountMax=3"

# --- GitHub API ---

# List Frontier runners from the GitHub API.
# Prints: id name status busy (one runner per line)
gh_list_runners() {
    gh api "orgs/$ORG/actions/runners" --paginate \
        --jq ".runners[]
              | select(.labels | map(.name) | index(\"$RUNNER_LABEL\"))
              | \"\(.id) \(.name) \(.status) \(.busy)\""
}

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

# Find all runner directories on shared storage.
# Prints: one directory path per line.
find_runner_dirs() {
    for conf in "$SHARED_DIR"/frontier-*/.runner; do
        [ -f "$conf" ] && dirname "$conf"
    done
}

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

# Find PIDs of a runner on a node by matching its CWD.
# (Runner.Listener's command line is just "Runner.Listener run" — no path.)
# Frontier SSH prints MOTD to stdout, so output is filtered to numeric lines only.
# Args: $1 = node, $2 = runner directory
# Prints: space-separated PIDs, or empty.
find_pids() {
    ssh $SSH_OPTS "$1" '
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || true)
            [ "$cwd" = "'"$2"'" ] && echo "$p"
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
# Args: $1 = node, $2 = runner directory
# Returns: 0 if running after start, 1 otherwise.
start_runner() {
    local node="$1" dir="$2"
    timeout 15 ssh $SSH_OPTS "$node" \
        "cd $dir && setsid nohup ./run.sh >> runner.log 2>&1 < /dev/null &" \
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

# Sweep all nodes in parallel and update runner.node for any runner
# found running on a different node than recorded. Called at the top of
# every primary script to ensure runner.node always reflects reality,
# even if a runner was manually restarted on a different node.
sync_runner_nodes() {
    local tmpdir
    tmpdir=$(mktemp -d)
    trap 'rm -rf "$tmpdir"' RETURN

    for node in "${NODES[@]}"; do
        (
            ssh $SSH_OPTS "$node" '
                for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
                    cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || true)
                    [ -n "$cwd" ] && echo "'"$node"' $cwd"
                done
            ' 2>/dev/null | grep -E '^[a-z0-9]+ /'
        ) > "$tmpdir/$node" &
    done
    wait

    while IFS=' ' read -r node dir; do
        [ -f "$dir/runner.node" ] || continue
        local recorded
        recorded=$(cat "$dir/runner.node" 2>/dev/null || echo "")
        if [ "$node" != "$recorded" ]; then
            echo "==> $(basename "$dir"): runner.node updated $recorded -> $node"
            echo "$node" > "$dir/runner.node"
        fi
    done < <(cat "$tmpdir"/*)
}
