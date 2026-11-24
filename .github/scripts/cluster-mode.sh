#!/bin/bash
# Decide what the HPC clusters run for this workflow run. Outputs:
#   mode   none   nothing on the clusters
#          smoke  two lanes, Frontier CCE OpenACC and Phoenix NVHPC OpenMP: build, then up to 120
#                 of the tests the change selects on each
#          full   every cluster lane, the whole suite, case-optimization and benchmarks
#   base   scheduled runs: the commit the last successful scheduled run tested
#
# Env: EVENT, REPO, RUN_ID, HEAD_SHA, PR_NUMBER, PR_DRAFT, CHECKALL, GH_TOKEN;
#      WORKFLOW (default test.yml) is the workflow whose scheduled runs are compared.
set -euo pipefail

LABEL="ci:clusters"
WORKFLOW="${WORKFLOW:-test.yml}"

# Changes that only a real cluster run can vouch for: full run without waiting for the label.
RISKY='^(src/common/include/[^/]*macros\.fpp|src/[^/]+/m_mpi_(common|proxy)\.fpp|CMakeLists\.txt|cmake/|toolchain/cmake/|toolchain/modules$|toolchain/templates/|\.github/scripts/|\.github/workflows/(common|frontier|frontier_amd)/|\.github/workflows/(test|bench|bisect)\.yml$)'
# Changes no cluster job reads: a PR made only of these runs nothing there.
IRRELEVANT='(\.md$|^(docs|misc|\.claude)/|^toolchain/mfc/(lint_|test_)[^/]*\.py$|^\.github/workflows/(docs|formatting|spelling|lint-toolchain|line-count|cleanliness|claude|claude-code-review|homebrew|homebrew-release|deploy-tap|docker|coverage|coverage-health|coverage-refresh|phoenix-alert|cluster-label|ci-record|runner-watch|test-toolchain-compat)\.yml$)'

mode=none base=""
case "$EVENT" in
    workflow_dispatch) mode=full ;;
    schedule)
        # Full run unless this commit has already passed a scheduled run, or failed two: a run lost to
        # an outage or a cancellation is retried, a commit that is really broken is not run forever.
        runs=$(gh api "repos/$REPO/actions/workflows/$WORKFLOW/runs?event=schedule&branch=master&per_page=50" \
            --jq "[.workflow_runs[] | select(.id != $RUN_ID) | {sha: .head_sha, conclusion}]")
        passed=$(jq --arg sha "$HEAD_SHA" '[.[] | select(.sha == $sha and .conclusion == "success")] | length' <<<"$runs")
        failed=$(jq --arg sha "$HEAD_SHA" '[.[] | select(.sha == $sha and .conclusion == "failure")] | length' <<<"$runs")
        if [ "$passed" -eq 0 ] && [ "$failed" -lt 2 ]; then mode=full; fi
        # The commit the last successful scheduled run tested: what a benchmark compares against.
        base=$(jq -r --arg sha "$HEAD_SHA" '[.[] | select(.sha != $sha and .conclusion == "success")][0].sha // ""' <<<"$runs")
        ;;
    pull_request)
        if [ "$PR_DRAFT" != "true" ]; then
            # Both read live, not from the event payload: a re-run after labelling carries the old payload.
            files=$(gh api --paginate "repos/$REPO/pulls/$PR_NUMBER/files?per_page=100" --jq '.[].filename')
            labels=$(gh api "repos/$REPO/issues/$PR_NUMBER/labels" --jq '.[].name')
            if grep -qxF "$LABEL" <<<"$labels" || grep -qE "$RISKY" <<<"$files"; then
                mode=full
            elif [ -n "$files" ] && [ -z "$(grep -vE "$IRRELEVANT" <<<"$files" || true)" ]; then
                mode=none
            elif [ "$CHECKALL" = "true" ]; then
                mode=smoke
            fi
        fi
        ;;
esac

echo "Cluster mode: $mode (base=${base:-none})"
if [ -n "${GITHUB_OUTPUT:-}" ]; then
    { echo "mode=$mode"; echo "base=$base"; } >> "$GITHUB_OUTPUT"
fi
