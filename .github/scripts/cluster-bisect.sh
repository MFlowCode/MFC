#!/bin/bash
# Find the commit that broke one cluster lane on master: git bisect between the last green
# scheduled run and the red one, building each candidate and running only the failed tests.
# Dispatched by phoenix-master-alert.sh through bisect.yml.
#
# The result is printed, and written to the step output "message" for bisect.yml to post.
#
# Env: CLUSTER (frontier|frontier_amd|phoenix), INTERFACE (acc|omp), GOOD, BAD, TESTS (UUIDs;
#      empty bisects the build).
#   cluster-bisect.sh        run the bisection
#   cluster-bisect.sh step   test the checked-out commit (called by git bisect run)
set -uo pipefail

if [ "${1:-}" = "step" ]; then
    lane=(gpu "$INTERFACE" "$CLUSTER")
    # 125 tells git bisect this commit cannot be judged.
    [ -f .github/scripts/cluster-lane.sh ] || exit 125
    rm -f ./*.out
    if [ "$CLUSTER" != "phoenix" ]; then
        bash ".github/workflows/$CLUSTER/build.sh" gpu "$INTERFACE" || exit 125
    fi
    rc=0
    bash .github/scripts/cluster-lane.sh build "${lane[@]}" || rc=$?
    # 78: the cluster is down, which says nothing about this commit.
    if [ "$rc" -eq 78 ]; then exit 125; fi
    if [ "$rc" -ne 0 ]; then
        if [ -n "$TESTS" ]; then exit 125; else exit 1; fi
    fi
    [ -n "$TESTS" ] || exit 0
    rc=0
    MFC_TEST_EXTRA_OPTS="--only $TESTS" bash .github/scripts/cluster-lane.sh test "${lane[@]}" || rc=$?
    [ "$rc" -ne 0 ] || exit 0
    # Bad only if a test failed. A job that failed without one (queue wait, time limit, a bad
    # node, the cluster down) says nothing about this commit.
    if [ -n "$(cat tests/failed_uuids*.txt 2>/dev/null)" ]; then exit 1; fi
    exit 125
fi

# The step script is run from a copy: bisect checks out other commits underneath it.
self=$(mktemp)
log=$(mktemp)
start=$(git rev-parse HEAD)
trap 'rm -f "$self" "$log"; git bisect reset >/dev/null 2>&1 || true; git checkout -q "$start" 2>/dev/null || true' EXIT
cp "$0" "$self"

report() {  # report <exit code> <text>: print it, hand it to the workflow, stop
    {
        printf 'Bisected %s gpu-%s between %s (green) and %s (red)' "$CLUSTER" "$INTERFACE" "${GOOD:0:9}" "${BAD:0:9}"
        if [ -n "$TESTS" ]; then printf ', on tests %s' "$TESTS"; else printf ', on the build'; fi
        printf '.\n\n%s\n' "$2"
    } | tee "$log.msg"
    if [ -n "${GITHUB_OUTPUT:-}" ]; then { echo 'message<<BISECT_EOF'; cat "$log.msg"; echo 'BISECT_EOF'; } >> "$GITHUB_OUTPUT"; fi
    rm -f "$log.msg"
    exit "$1"
}

# git bisect never tests the bad end itself, and names it when every candidate passes. Test it
# first: a flaky test or an infrastructure failure does not reproduce, and has no culprit.
git checkout -q "$BAD" || exit 1
rc=0
bash "$self" step || rc=$?
case "$rc" in
    0)   report 1 "The failure did not reproduce at ${BAD:0:9}: no commit to blame. A flaky test or an infrastructure failure looks like this." ;;
    125) report 1 "${BAD:0:9} could not be judged (the cluster was down, or the job failed without a failed test), so nothing was bisected." ;;
esac

git bisect start --first-parent "$BAD" "$GOOD" || exit 1
git bisect run bash "$self" step | tee "$log"
culprit=$(sed -n 's/^\([0-9a-f]\{40\}\) is the first bad commit$/\1/p' "$log" | head -1)
if [ -n "$culprit" ]; then
    report 0 "First bad commit: $culprit $(git log -1 --format=%s "$culprit")"
fi
report 1 "No single commit found: the candidates in between could not be judged."
