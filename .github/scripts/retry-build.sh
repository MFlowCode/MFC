#!/bin/bash
# Provides retry_build(): 2-attempt loop.
# On failure of attempt 1, nukes the entire build directory before attempt 2.
# Usage: source .github/scripts/retry-build.sh
#        retry_build ./mfc.sh build -j 8 --gpu acc

retry_build() {
    local max_attempts=2
    local attempt=1
    while [ $attempt -le $max_attempts ]; do
        echo "Build attempt $attempt of $max_attempts..."
        if "$@"; then
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
