#!/usr/bin/env bash
# Stop and deregister a Phoenix runner.
# Thin wrapper — see misc/common/stop-runner.sh for the implementation.
#
# Usage: bash stop-runner.sh <runner-name>
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/../common/stop-runner.sh"
