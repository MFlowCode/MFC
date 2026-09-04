#!/bin/bash
# Provides retry_build(): 2-attempt loop.
# On failure of attempt 1, nukes the build directory before attempt 2, keeping
# build/venv: a compute node cannot reinstall it (no route to PyPI), so removing
# it made every retry fail on a dependency fetch that could not succeed (#1813).
# If RETRY_VALIDATE_CMD is set, runs it after a successful build; a non-zero
# exit triggers the same nuke-and-retry, catching e.g. SIGILL from binaries
# compiled on a different CPU architecture.
# Usage: source .github/scripts/retry-build.sh
#        retry_build ./mfc.sh build -j 8 --gpu acc
#        RETRY_VALIDATE_CMD='./syscheck' retry_build ./mfc.sh build -j 8

# Delay between build attempts. Overridable so tests can exercise the retry
# path without waiting on it; CI leaves it at the default.
: "${MFC_BUILD_RETRY_DELAY:=30}"

nuke_build() {
    find build -mindepth 1 -maxdepth 1 ! -name venv -exec rm -rf -- {} + 2>/dev/null || true
}

retry_build() {
    local max_attempts=2
    local validate_cmd="${RETRY_VALIDATE_CMD:-}"
    local attempt=1
    while [ $attempt -le $max_attempts ]; do
        echo "Build attempt $attempt of $max_attempts..."
        if "$@"; then
            if [ -n "$validate_cmd" ]; then
                if ! eval "$validate_cmd"; then
                    echo "Post-build validation failed on attempt $attempt."
                    if [ $attempt -lt $max_attempts ]; then
                        echo "  Clearing the build directory (keeping build/venv) before retry..."
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
            echo "  Build failed — clearing the build directory (keeping build/venv) before retry..."
            nuke_build
            sleep "$MFC_BUILD_RETRY_DELAY"
        else
            echo "Build failed after $max_attempts attempts."
            return 1
        fi
        attempt=$((attempt + 1))
    done
}
