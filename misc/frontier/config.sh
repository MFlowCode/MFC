#!/usr/bin/env bash
# Shared configuration for Frontier GitHub Actions runner management.
#
# Sourced by all other scripts. Provides Frontier constants, GitHub API
# helpers, and login-node process management functions.

# --- Frontier constants ---
ORG="MFlowCode"
RUNNER_GROUP="phoenix"
RUNNER_LABEL="frontier"
NODES=(login01 login02 login03 login04 login05 login06 login07 login08 login09 login10 login11)
SHARED_DIR="/lustre/orion/cfd154/proj-shared/runners"

SSH_OPTS="-o StrictHostKeyChecking=no -o ConnectTimeout=10 -o BatchMode=yes -o ServerAliveInterval=10 -o ServerAliveCountMax=3"

source "$(dirname "${BASH_SOURCE[0]}")/../common/runner-lib.sh"

# --- GitHub API ---

# List Frontier runners from the GitHub API.
# Prints: id name status busy (one runner per line)
gh_list_runners() {
    gh api "orgs/$ORG/actions/runners" --paginate \
        --jq ".runners[]
              | select(.labels | map(.name) | index(\"$RUNNER_LABEL\"))
              | \"\(.id) \(.name) \(.status) \(.busy)\""
}

# --- Local filesystem ---

# Find all runner directories on shared storage.
# Prints: one directory path per line.
find_runner_dirs() {
    for conf in "$SHARED_DIR"/frontier-*/.runner; do
        [ -f "$conf" ] && dirname "$conf"
    done
}

# --- Login-node process management ---

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
