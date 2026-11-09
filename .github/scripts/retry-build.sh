#!/bin/bash
# Provides retry_build(): 2-attempt loop.
# On failure of attempt 1, removes the build/staging and build/install slugs that
# attempt touched (the toolchain lists them in $MFC_SLUG_LOG) before attempt 2.
# Nothing else goes: a compute node cannot reinstall build/venv or prebuilt
# dependencies (no route to PyPI or the sources, #1813), the venv's install marker
# build/pyproject.toml must survive or the next ./mfc.sh reinstalls, and Frontier's
# concurrent build variants share build/, so another job's slugs must survive too.
# If RETRY_VALIDATE_CMD is set, runs it after a successful build; a non-zero
# exit triggers the same clear-and-retry, catching e.g. SIGILL from binaries
# compiled on a different CPU architecture.
# Usage: source .github/scripts/retry-build.sh
#        retry_build ./mfc.sh build -j 8 --gpu acc
#        RETRY_VALIDATE_CMD='./syscheck' retry_build ./mfc.sh build -j 8

# Delay between build attempts. Overridable so tests can exercise the retry
# path without waiting on it; CI leaves it at the default.
: "${MFC_BUILD_RETRY_DELAY:=30}"

nuke_build() {
    local slug
    sort -u "$MFC_SLUG_LOG" | while read -r slug; do
        [ -n "$slug" ] && rm -rf -- "build/staging/$slug" "build/install/$slug"
    done
    : > "$MFC_SLUG_LOG"
}

retry_build() {
    local max_attempts=2
    local validate_cmd="${RETRY_VALIDATE_CMD:-}"
    local attempt=1
    export MFC_SLUG_LOG
    MFC_SLUG_LOG=$(mktemp)
    while [ $attempt -le $max_attempts ]; do
        echo "Build attempt $attempt of $max_attempts..."
        if "$@"; then
            if [ -n "$validate_cmd" ]; then
                if ! eval "$validate_cmd"; then
                    echo "Post-build validation failed on attempt $attempt."
                    if [ $attempt -lt $max_attempts ]; then
                        echo "  Clearing the slugs this attempt built before retry..."
                        nuke_build
                        sleep 5
                        attempt=$((attempt + 1))
                        continue
                    else
                        echo "Validation still failing after $max_attempts attempts."
                        return 1
                    fi
                fi
            fi
            echo "Build succeeded on attempt $attempt."
            return 0
        fi
        if [ $attempt -lt $max_attempts ]; then
            echo "  Build failed — clearing the slugs this attempt built before retry..."
            nuke_build
            sleep "$MFC_BUILD_RETRY_DELAY"
        else
            echo "Build failed after $max_attempts attempts."
            return 1
        fi
        attempt=$((attempt + 1))
    done
}
