#!/usr/bin/env bash
# Dispatcher for GitHub Actions runner management scripts.
#
# Loads site-specific configuration then runs the requested command.
# Common commands live in misc/common/; site-specific commands live in
# misc/<site>/. All site-specific scripts source their own config, so
# the dispatcher only pre-loads config for common commands.
#
# Usage: bash misc/runners/runner.sh <site> <command> [args...]
#
# Sites:    frontier  phoenix
# Common:   check-runners  list-runners  move-runner  rebalance-runners
#           restart-all  restart-runner  stop-runner  rerun-failed
# Frontier: make-runner  deploy-runners
# Phoenix:  create-runner
#
# Examples:
#   bash misc/runners/runner.sh frontier check-runners
#   bash misc/runners/runner.sh phoenix list-runners
#   APPLY=1 bash misc/runners/runner.sh frontier rebalance-runners
#   bash misc/runners/runner.sh frontier restart-runner login01 /path/to/runner
#   bash misc/runners/runner.sh frontier make-runner 23 login01
#   bash misc/runners/runner.sh phoenix create-runner phoenix-11 login-phoenix-gnr-2

set -euo pipefail

MISC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SITE="${1:?Usage: $0 <site> <command> [args...]}"
CMD="${2:?Usage: $0 <site> <command> [args...]}"
shift 2

if [ ! -f "$MISC_DIR/$SITE/config.sh" ]; then
    echo "ERROR: Unknown site '$SITE'. Known sites: frontier, phoenix" >&2
    exit 1
fi

# Site-specific scripts are standalone — they source their own config.
if [ -f "$MISC_DIR/$SITE/$CMD.sh" ]; then
    exec bash "$MISC_DIR/$SITE/$CMD.sh" "$@"
fi

# Common scripts need the site config pre-loaded.
if [ -f "$MISC_DIR/common/$CMD.sh" ]; then
    SITE_SCRIPT_DIR="$MISC_DIR/$SITE"
    source "$MISC_DIR/$SITE/config.sh"
    source "$MISC_DIR/common/$CMD.sh" "$@"
    exit
fi

echo "ERROR: Unknown command '$CMD' for site '$SITE'." >&2
echo "Common:   check-runners list-runners move-runner rebalance-runners restart-all restart-runner stop-runner rerun-failed" >&2
echo "Frontier: make-runner deploy-runners" >&2
echo "Phoenix:  create-runner" >&2
exit 1
