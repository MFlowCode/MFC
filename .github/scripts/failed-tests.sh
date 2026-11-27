#!/bin/bash
# The tests a `./mfc.sh test -v` log says failed, as a Markdown list.
# Usage: failed-tests.sh <log>
failed=$(grep -a -A3 'FAILED:' "$1" | grep -ao 'Error: Test tests/.*' | sed 's/^Error: Test /- /; s/: Failed to execute MFC\.$//' | sort -u | head -40)
echo "${failed:-The job failed before any test did.}"
