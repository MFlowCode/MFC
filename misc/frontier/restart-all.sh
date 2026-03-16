#!/usr/bin/env bash
# Restart all Frontier runners in place on their current nodes.
# Thin wrapper — see misc/common/restart-all.sh for the implementation.
#
# Usage: bash restart-all.sh              # dry run
#        APPLY=1 bash restart-all.sh      # execute
#        APPLY=1 FORCE=1 bash restart-all.sh  # restart busy runners too
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/../common/restart-all.sh"
