#!/bin/bash
# Decide whether a failed build was a cluster-wide dependency outage.
#
# MFC bootstraps its Python toolchain into build/venv on the first ./mfc.sh call
# of a job, pulling from pypi.org. On Phoenix clean_build has just moved build/
# aside, so that happens every time; on Frontier it happens in the login-node
# "Fetch Dependencies" step. When the index is unreachable the build fails for a
# reason no other node improves on, so it is worth recording once and skipping
# the rest of the matrix rather than having each job spend ~33 minutes
# rediscovering it (17 Frontier jobs did exactly that on 2026-08-28).
#
# Usage: classify-build-failure.sh <logfile> <cluster>
#
# Exit codes:
#   78  cluster-wide dependency outage; it has been recorded
#   0   ordinary build failure, caller should keep its own exit code

set -uo pipefail

log="${1:-}"
cluster="${2:-}"

if [ -z "$log" ] || [ -z "$cluster" ]; then
    echo "Usage: $0 <logfile> <cluster>" >&2
    exit 0
fi

# No log means nothing to classify. Never claim an outage on absent evidence.
[ -f "$log" ] || exit 0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The URL may or may not be wrapped (uv quotes it in backticks, plain pip does
# not), so do not require a character between the colon and the scheme.
if grep -qE "Failed to fetch:[^h]*https?://pypi|uv install failed|\(venv\) Installation failed" "$log"; then
    bash "$SCRIPT_DIR/ci-outage.sh" mark "$cluster" \
        "PyPI/uv dependency install failed during build"
    exit 78
fi

exit 0
