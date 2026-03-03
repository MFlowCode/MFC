#!/bin/bash
# Runs ./mfc.sh test with all provided arguments, then retries a small number
# of sporadic failures (up to 5). Exits non-zero on real failures.
# Usage: bash .github/scripts/run-tests-with-retry.sh [mfc test args...]

# Extract flags that should carry over to retries (retries build their own
# argument list with --only, so we capture passthrough flags here).
PASSTHROUGH=""
for arg in "$@"; do
    case "$arg" in
        --test-all) PASSTHROUGH="$PASSTHROUGH --test-all" ;;
    esac
done

rm -f tests/failed_uuids.txt
TEST_EXIT=0
/bin/bash mfc.sh test "$@" || TEST_EXIT=$?

# Retry only if a small number of tests failed (sporadic failures)
if [ -s tests/failed_uuids.txt ]; then
    NUM_FAILED=$(wc -l < tests/failed_uuids.txt)
    if [ "$NUM_FAILED" -le 5 ]; then
        FAILED=$(tr '\n' ' ' < tests/failed_uuids.txt)
        echo ""
        echo "=== Retrying $NUM_FAILED failed test(s): $FAILED ==="
        echo ""
        /bin/bash mfc.sh test -v --max-attempts 3 -j "$(nproc)" --only $FAILED $PASSTHROUGH || exit $?
    else
        echo "Too many failures ($NUM_FAILED) to retry — likely a real issue."
        exit 1
    fi
elif [ "$TEST_EXIT" -ne 0 ]; then
    exit $TEST_EXIT
fi
