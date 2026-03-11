#!/bin/bash
# Provides retry_build(): 2-attempt loop.
# On failure of attempt 1, nukes the entire build directory before attempt 2.
# If RETRY_VALIDATE_CMD is set, runs it after a successful build; a non-zero
# exit triggers the same nuke-and-retry, catching e.g. SIGILL from binaries
# compiled on a different CPU architecture.
# Usage: source .github/scripts/retry-build.sh
#        retry_build ./mfc.sh build -j 8 --gpu acc
#        RETRY_VALIDATE_CMD='./syscheck' retry_build ./mfc.sh build -j 8

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
                        echo "  Nuking build directory before retry..."
                        rm -rf build 2>/dev/null || true
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
            echo "  Build failed — nuking build directory before retry..."
            rm -rf build 2>/dev/null || true
            sleep 30
        else
            echo "Build failed after $max_attempts attempts."
            return 1
        fi
        attempt=$((attempt + 1))
    done
}
