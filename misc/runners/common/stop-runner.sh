#!/usr/bin/env bash
# Stop and deregister a GitHub Actions runner.
#
# Sourced by site wrappers (frontier/stop-runner.sh, phoenix/stop-runner.sh)
# after config.sh is loaded. Finds the runner directory by name via
# find_runner_dirs(), stops the process, and removes the GitHub registration.
#
# Usage: bash stop-runner.sh <runner-name>
set -euo pipefail

RUNNER_NAME="${1:?Usage: $0 <runner-name>}"

# Find runner directory by name
runner_dir=""
while IFS= read -r dir; do
    if [ "$(get_runner_name "$dir")" = "$RUNNER_NAME" ]; then
        runner_dir="$dir"
        break
    fi
done < <(find_runner_dirs)

if [ -z "$runner_dir" ]; then
    echo "ERROR: Runner '$RUNNER_NAME' not found in known runner directories." >&2
    exit 1
fi

# Locate and stop the process
echo "==> Locating $RUNNER_NAME..."
node=$(find_node "$runner_dir")

if [ "$node" != "offline" ]; then
    echo "==> Stopping $RUNNER_NAME on $node..."
    stop_runner "$node" "$runner_dir"
    echo "==> Process stopped."
else
    echo "==> $RUNNER_NAME is not running (already offline)."
fi

# Deregister from GitHub
echo "==> Fetching runner ID from GitHub..."
runner_id=""
runner_list=$(gh_list_runners 2>/dev/null) || {
    echo "WARNING: GitHub API call failed; runner may still be registered on GitHub." >&2
    exit 0
}
while read -r id name _status _busy; do
    [ "$name" = "$RUNNER_NAME" ] && runner_id="$id" && break
done <<< "$runner_list"

if [ -n "$runner_id" ]; then
    echo "==> Deregistering runner (ID $runner_id)..."
    gh_remove_runner "$runner_id"
    echo "==> Done."
else
    echo "==> Runner not found in GitHub API (may already be deregistered)."
fi
