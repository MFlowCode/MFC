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
        while read -r run_id run_name; do
            # Check which jobs failed; skip if run has expired or been deleted
            failed_jobs=$(gh run view --repo "$REPO" "$run_id" \
                --json jobs --jq '.jobs[] | select(.conclusion == "failure" or .conclusion == "cancelled") | .name' \
                2>/dev/null) || { echo "  WARNING: could not fetch jobs for run $run_id, skipping"; continue; }
            echo "  Run $run_id ($run_name):"
            if [ -n "$failed_jobs" ]; then
                while read -r job; do
                    echo "    - $job"
                done <<< "$failed_jobs"
            fi

            if [ "${APPLY:-0}" = "1" ]; then
                echo "  Rerunning failed jobs..."
                gh run rerun --repo "$REPO" "$run_id" --failed || echo "  WARNING: rerun failed (may already be rerunning)"
            fi
        done < <(echo "$failed_runs")
    fi
done

# Also check master branch
echo ""
echo "=== master branch ==="
master_failed=$(gh run list --repo "$REPO" --branch master --limit 5 \
    --json databaseId,status,conclusion,name \
    --jq '.[] | select(.conclusion == "failure") | "\(.databaseId) \(.name)"')
if [ -n "$master_failed" ]; then
    while read -r run_id run_name; do
        failed_jobs=$(gh run view --repo "$REPO" "$run_id" \
            --json jobs --jq '.jobs[] | select(.conclusion == "failure" or .conclusion == "cancelled") | .name' \
            2>/dev/null) || { echo "  WARNING: could not fetch jobs for run $run_id, skipping"; continue; }
        echo "  Run $run_id ($run_name):"
        if [ -n "$failed_jobs" ]; then
            while read -r job; do
                echo "    - $job"
            done <<< "$failed_jobs"
        fi
        if [ "${APPLY:-0}" = "1" ]; then
            echo "  Rerunning failed jobs..."
            gh run rerun --repo "$REPO" "$run_id" --failed || echo "  WARNING: rerun failed"
        fi
    done < <(echo "$master_failed")
else
    echo "  All passing"
fi

if [ "${APPLY:-0}" != "1" ]; then
    echo ""
    echo "Dry run — set APPLY=1 to actually rerun failed workflows."
fi
