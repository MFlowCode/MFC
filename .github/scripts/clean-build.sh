#!/bin/bash
# Provides clean_build(): renames build/ aside and deletes it in the background.
# mv is a metadata-only operation that succeeds even with stale NFS file handles,
# unlike rm -rf which fails on ESTALE. The background delete is best-effort and
# scoped to this job's PID to avoid races with concurrent matrix jobs.
#
# Usage: source .github/scripts/clean-build.sh
#        clean_build

clean_build() {
    # Clean up leftover stale directories from previous runs before adding a new one.
    rm -rf build.stale.* 2>/dev/null || true
    mv build "build.stale.$$" 2>/dev/null || true
    rm -rf "build.stale.$$" 2>/dev/null & disown
}
