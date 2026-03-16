#!/usr/bin/env bash
# Check runner health across all Frontier login nodes.
# Thin wrapper — see misc/common/check-runners.sh for the implementation.
#
# Usage: bash check-runners.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/../common/check-runners.sh"
