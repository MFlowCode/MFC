#!/bin/bash
# Decide what the HPC clusters run for this workflow run. Outputs:
#   mode   none   nothing on the clusters
#          smoke  two lanes, Frontier CCE OpenACC and Phoenix NVHPC OpenMP: build, then up to 120
#                 of the tests the change selects on each
#          full   every cluster lane, the whole suite, case-optimization and benchmarks
#   base   scheduled runs: the commit of the last scheduled run that counts as passed
#
# Env: EVENT, REPO, RUN_ID, HEAD_SHA, PR_NUMBER, PR_DRAFT, CHECKALL, GH_TOKEN;
#      ACTION and PR_HEAD_SHA (bench.yml: a run started by a label does nothing when the pull
#      request's commit has a benchmark already, running or done);
#      WORKFLOW (default test.yml) is the workflow whose scheduled runs are compared;
#      PARTIAL=1 counts a run as passed once any cluster job passed in it (bench.yml: each lane
#      records its own result, and one broken lane must not hold the base back or repeat the rest).
set -euo pipefail

LABEL="ci:clusters"
WORKFLOW="${WORKFLOW:-test.yml}"

# Changes that only a real cluster run can vouch for: full run without waiting for the label.
# So is any source change that adds or removes a line with a GPU_* macro (see below).
RISKY='^(src/common/include/[^/]*macros\.fpp|src/[^/]+/m_mpi_(common|proxy)\.fpp|CMakeLists\.txt|cmake/|toolchain/cmake/|toolchain/modules$|toolchain/templates/|\.github/scripts/|\.github/workflows/(common|frontier|frontier_amd)/|\.github/workflows/(test|bench|bisect)\.yml$)'
# Changes no cluster job reads: a PR made only of these runs nothing there.
IRRELEVANT='(\.md$|^(docs|misc|\.claude)/|^toolchain/mfc/(lint_|test_)[^/]*\.py$|^\.github/workflows/(docs|formatting|spelling|lint-toolchain|line-count|cleanliness|claude|claude-code-review|homebrew|homebrew-release|deploy-tap|docker|coverage|coverage-health|coverage-refresh|phoenix-alert|cluster-label|ci-record|runner-watch|test-toolchain-compat)\.yml$)'

mode=none base=""
case "$EVENT" in
    workflow_dispatch) mode=full ;;
    schedule)
        # Full run unless this commit has already passed a scheduled run, or failed two: a run lost to
        # an outage or a cancellation is retried, a commit that is really broken is not run forever.
        # Judged by the cluster jobs, not the run: the Frontier lanes are continue-on-error, so a run
        # they failed in still concludes "success".
        runs=$(gh api "repos/$REPO/actions/workflows/$WORKFLOW/runs?event=schedule&branch=master&per_page=50" \
            --jq "[.workflow_runs[] | select(.id != $RUN_ID) | {id, sha: .head_sha}]")
        passed=0 failed=0
        for id in $(jq -r '.[].id' <<<"$runs"); do
            sha=$(jq -r --argjson id "$id" '.[] | select(.id == $id) | .sha' <<<"$runs")
            [ "$sha" = "$HEAD_SHA" ] || [ -z "$base" ] || continue
            verdict=$(gh api --paginate "repos/$REPO/actions/runs/$id/jobs?per_page=100" \
                --jq '[.jobs[] | select(.labels | any(. == "frontier" or . == "gt")) | .conclusion]' | jq -rs 'add
                    | if length == 0 or all(. == "skipped") then "none" elif all(. == "success" or . == "skipped") then "passed"
                      elif $partial == "1" and any(. == "success") then "passed"
                      elif any(. == "failure") then "failed" else "unfinished" end' --arg partial "${PARTIAL:-0}")
            if [ "$sha" = "$HEAD_SHA" ]; then
                case "$verdict" in passed) passed=$((passed + 1)) ;; failed) failed=$((failed + 1)) ;; esac
            elif [ "$verdict" = passed ]; then
                base=$sha       # the commit the last passing scheduled run tested: what a benchmark compares against
            fi
        done
        if [ "$passed" -eq 0 ] && [ "$failed" -lt 2 ]; then mode=full; fi
        ;;
    pull_request)
        if [ "$PR_DRAFT" != "true" ]; then
            # Both read live, not from the event payload: a re-run after labelling carries the old payload.
            changes=$(gh api --paginate "repos/$REPO/pulls/$PR_NUMBER/files?per_page=100" | jq -s 'add')
            files=$(jq -r '.[].filename' <<<"$changes")
            # Source files whose diff adds or removes a line with a GPU macro (GPU_PARALLEL_LOOP, GPU_ROUTINE,
            # GPU_UPDATE, GPU_DECLARE, ...): where the compilers differ most. A diff too large for the API to return counts, since it cannot be read.
            gpu=$(jq '[.[] | select(.filename | test("^src/.*\\.fpp$")) | select(.status != "removed")
                | select(.patch == null or (.patch | split("\n") | any(test("^[+-].*GPU_[A-Z_]+"))))] | length' <<<"$changes")
            labels=$(gh api "repos/$REPO/issues/$PR_NUMBER/labels" --jq '.[].name')
            # Tests whose own files changed all run, beyond the smoke cap: past this many, one unsharded job cannot hold them.
            goldens=$({ grep -oE '^tests/[0-9A-F]{8}/' <<<"$files" || true; } | sort -u | wc -l)
            if grep -qxF "$LABEL" <<<"$labels" || grep -qE "$RISKY" <<<"$files" || [ "$gpu" -gt 0 ] || [ "$goldens" -gt 100 ]; then
                mode=full
            elif [ -n "$files" ] && [ -z "$(grep -vE "$IRRELEVANT" <<<"$files" || true)" ]; then
                mode=none
            elif [ "$CHECKALL" = "true" ]; then
                mode=smoke
            fi
        fi
        ;;
esac

if [ "$mode" = full ] && [ "${ACTION:-}" = labeled ] && [ -n "${PR_HEAD_SHA:-}" ]; then
    for id in $(gh api "repos/$REPO/actions/workflows/$WORKFLOW/runs?event=pull_request&head_sha=$PR_HEAD_SHA&per_page=20" \
        --jq ".workflow_runs[] | select(.id != $RUN_ID) | .id"); do
        ran=$(gh api --paginate "repos/$REPO/actions/runs/$id/jobs?per_page=100" \
            --jq '[.jobs[] | select(.labels | any(. == "frontier" or . == "gt")) | select(.conclusion != "skipped" and .conclusion != "cancelled")] | length' | jq -s add)
        if [ "${ran:-0}" -gt 0 ]; then echo "Run $id already benchmarks $PR_HEAD_SHA."; mode=none; break; fi
    done
fi

echo "Cluster mode: $mode (base=${base:-none})"
if [ -n "${GITHUB_OUTPUT:-}" ]; then
    { echo "mode=$mode"; echo "base=$base"; } >> "$GITHUB_OUTPUT"
fi
