#!/bin/bash
# Provides retry_build(): 3-attempt loop with configurable cleanup.
# Set RETRY_CLEAN_CMD to override cleanup (default: rm -rf build/staging build/install build/lock.yaml).
# Usage: source .github/scripts/retry-build.sh
#        retry_build ./mfc.sh build -j 8 --gpu acc

retry_build() {
    local clean_cmd="${RETRY_CLEAN_CMD:-rm -rf build/staging build/install build/lock.yaml}"
    local max_attempts=3
    local attempt=1
    while [ $attempt -le $max_attempts ]; do
        echo "Build attempt $attempt of $max_attempts..."
        if "$@"; then
            echo "Build succeeded on attempt $attempt."
            return 0
        fi
        if [ $attempt -lt $max_attempts ]; then
            echo "Build failed on attempt $attempt. Retrying in 30s..."
            eval "$clean_cmd"
            sleep 30
        else
            echo "Build failed after $max_attempts attempts."
            return 1
        fi
        attempt=$((attempt + 1))
    done
}
