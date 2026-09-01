#!/bin/bash
# Aggregate retry outcomes across recent CI runs (see issue #1798).
#
# The test suite emits one line per run:
#   MFC_RETRY_STATS rescued=N failed=N passed=N by_class=execution failed:2,...
#
# This collects those lines over the last N Test Suite runs and totals them, so
# "are retries worth their cost" is answered with numbers instead of argument.
# Reading a hundred logs by hand is what made the question unanswerable before.
#
# Usage: bash .github/scripts/harvest-retry-stats.sh [runs]   (default 40)
#
# Interpreting the result: retries cost roughly 3x the wall clock of the
# failures they act on. A class with a rescue count at or near zero across a few
# hundred retried cases is paying that cost for nothing and is a candidate for
# skipping the retry entirely.

set -uo pipefail

REPO="${MFC_REPO:-MFlowCode/MFC}"
WORKFLOW="${MFC_WORKFLOW:-test.yml}"
LIMIT="${1:-40}"

command -v gh >/dev/null 2>&1 || { echo "gh CLI is required" >&2; exit 1; }

echo "Harvesting MFC_RETRY_STATS from the last $LIMIT $WORKFLOW runs of $REPO..." >&2

runs=$(gh run list --repo "$REPO" --workflow "$WORKFLOW" --limit "$LIMIT" \
        --json databaseId --jq '.[].databaseId' 2>/dev/null)

if [ -z "$runs" ]; then
    echo "No runs found." >&2
    exit 1
fi

stats=$(
    for run in $runs; do
        for job in $(gh api "repos/$REPO/actions/runs/$run/jobs?per_page=100" \
                        --jq '.jobs[] | select(.conclusion != "skipped") | .id' 2>/dev/null); do
            gh api "repos/$REPO/actions/jobs/$job/logs" 2>/dev/null \
                | tr -d '\r' | grep -o 'MFC_RETRY_STATS .*' || true
        done
    done
)

if [ -z "$stats" ]; then
    echo "No MFC_RETRY_STATS lines found (runs may predate the instrumentation)." >&2
    exit 1
fi

printf '%s\n' "$stats" | awk '
{
    for (i = 1; i <= NF; i++) {
        if ($i ~ /^rescued=/)  { split($i, a, "="); rescued += a[2] }
        if ($i ~ /^failed=/)   { split($i, a, "="); failed  += a[2] }
        if ($i ~ /^passed=/)   { split($i, a, "="); passed  += a[2] }
    }
    if (match($0, /by_class=.*/)) {
        cls = substr($0, RSTART + 9)
        if (cls != "none") {
            n = split(cls, parts, ",")
            for (i = 1; i <= n; i++) {
                split(parts[i], kv, ":")
                by[kv[1]] += kv[2]
            }
        }
    }
    runs++
}
END {
    printf "\nruns sampled:        %d\n", runs
    printf "tests passed:        %d\n", passed
    printf "tests failed:        %d\n", failed
    printf "rescued by a retry:  %d\n", rescued
    if (rescued > 0) {
        print  "\nrescues by failure class:"
        for (k in by) printf "  %6d  %s\n", by[k], k
    } else {
        print  "\nNo retry rescued anything in this sample."
        print  "Every retried case failed all of its attempts, so the retries"
        print  "cost roughly 3x the wall clock of those failures for no benefit."
    }
}'
