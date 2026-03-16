#!/usr/bin/env bash
# Automatically rebalance Frontier runners across login nodes.
# Thin wrapper — see misc/common/rebalance-runners.sh for the implementation.
#
# Usage: bash rebalance-runners.sh              # dry run
#        APPLY=1 bash rebalance-runners.sh      # execute
#        APPLY=1 FORCE=1 bash rebalance-runners.sh  # move busy runners too
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE_SCRIPT_DIR="$SCRIPT_DIR"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/../common/rebalance-runners.sh"
