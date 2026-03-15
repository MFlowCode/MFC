#!/usr/bin/env bash
# Rerun failed GitHub Actions workflows on open non-draft MFC PRs and master.
#
# Checks the 5 most recent workflow runs per branch.  Only the failed jobs
# within each run are rerun (via `gh run rerun --failed`), not the entire
# workflow.  Runs that are already in progress or queued are skipped by `gh`.
#
# Requires: gh CLI authenticated with access to MFlowCode/MFC
#
# Usage: bash rerun-failed.sh          # dry run (show what would be rerun)
#        APPLY=1 bash rerun-failed.sh  # actually rerun failed workflows

set -euo pipefail

REPO="MFlowCode/MFC"

echo "Checking open non-draft PRs on $REPO..."
prs=$(gh pr list --repo "$REPO" --state open --json number,title,isDraft --jq '.[] | select(.isDraft == false) | .number')

if [ -z "$prs" ]; then
    echo "No open non-draft PRs found."
    exit 0
fi

rerun_count=0
for pr in $prs; do
    title=$(gh pr view --repo "$REPO" "$pr" --json title --jq .title)
    branch=$(gh pr view --repo "$REPO" "$pr" --json headRefName --jq .headRefName)

    # Find failed workflow runs on this PR's branch
    failed_runs=$(gh run list --repo "$REPO" --branch "$branch" --limit 5 \
        --json databaseId,status,conclusion,name \
        --jq '.[] | select(.conclusion == "failure") | "\(.databaseId) \(.name)"')

    if [ -n "$failed_runs" ]; then
        echo ""
        echo "=== PR #$pr: $title ==="
        echo "$failed_runs" | while read -r run_id run_name; do
            # Check which jobs failed
            failed_jobs=$(gh run view --repo "$REPO" "$run_id" \
                --json jobs --jq '.jobs[] | select(.conclusion == "failure" or .conclusion == "cancelled") | .name')
            echo "  Run $run_id ($run_name):"
            echo "$failed_jobs" | while read -r job; do
                echo "    - $job"
            done

            if [ "${APPLY:-0}" = "1" ]; then
                echo "  Rerunning failed jobs..."
                gh run rerun --repo "$REPO" "$run_id" --failed || echo "  WARNING: rerun failed (may already be rerunning)"
                rerun_count=$((rerun_count + 1))
            fi
        done
    fi
done

# Also check master branch
echo ""
echo "=== master branch ==="
master_failed=$(gh run list --repo "$REPO" --branch master --limit 5 \
    --json databaseId,status,conclusion,name \
    --jq '.[] | select(.conclusion == "failure") | "\(.databaseId) \(.name)"')
if [ -n "$master_failed" ]; then
    echo "$master_failed" | while read -r run_id run_name; do
        echo "  Run $run_id ($run_name)"
        if [ "${APPLY:-0}" = "1" ]; then
            echo "  Rerunning failed jobs..."
            gh run rerun --repo "$REPO" "$run_id" --failed || echo "  WARNING: rerun failed"
        fi
    done
else
    echo "  All passing"
fi

if [ "${APPLY:-0}" != "1" ]; then
    echo ""
    echo "Dry run — set APPLY=1 to actually rerun failed workflows."
fi
