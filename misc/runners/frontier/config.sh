#!/usr/bin/env bash
# Shared configuration for Frontier GitHub Actions runner management.
#
# Sourced by all other scripts. Provides Frontier constants and
# site-specific functions. Common functions live in ../common/runner-lib.sh.

# --- Frontier constants ---
ORG="MFlowCode"
RUNNER_GROUP="phoenix"  # Both sites share one GitHub runner group named "phoenix"
RUNNER_LABEL="frontier"
NODES=(login01 login02 login03 login04 login05 login06 login07 login08 login09 login10 login11)
SHARED_DIR="/lustre/orion/cfd154/proj-shared/runners"

SSH_OPTS="-o StrictHostKeyChecking=no -o ConnectTimeout=10 -o BatchMode=yes -o ServerAliveInterval=10 -o ServerAliveCountMax=3"

source "$(dirname "${BASH_SOURCE[0]}")/../common/runner-lib.sh"

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
    local tmpdir node
    tmpdir=$(mktemp -d)
    trap 'rm -rf "$tmpdir"' RETURN

    sweep_all_nodes "$tmpdir"

    for node in "${NODES[@]}"; do
        while IFS= read -r line; do
            local dir sweep_node
            read -r _s sweep_node dir _rss _slurm <<< "$line"
            [ -f "$dir/runner.node" ] || continue
            local recorded
            recorded=$(cat "$dir/runner.node" 2>/dev/null || echo "")
            if [ "$sweep_node" != "$recorded" ]; then
                echo "==> $(basename "$dir"): runner.node updated $recorded -> $sweep_node"
                echo "$sweep_node" > "$dir/runner.node"
            fi
        done < <(grep '^RUNNER ' "$tmpdir/$node.out" 2>/dev/null || true)
    done
}
