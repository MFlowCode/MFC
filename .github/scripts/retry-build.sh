#!/bin/bash
# Provides retry_build(): 3-attempt loop with configurable cleanup.
# Set RETRY_CLEAN_CMD to override cleanup (default: rm -rf build/staging build/install build/lock.yaml).
# Set RETRY_VALIDATE_CMD to run a post-build validation; failure triggers a retry.
# Usage: source .github/scripts/retry-build.sh
#        retry_build ./mfc.sh build -j 8 --gpu acc

_retry_clean() {
    local clean_cmd="$1"
    if eval "$clean_cmd" 2>/dev/null; then
        return 0
    fi
    echo "  Cleanup failed; falling back to best-effort rm."
    rm -rf build/staging build/install build/lock.yaml 2>/dev/null || true
}

retry_build() {
    local clean_cmd="${RETRY_CLEAN_CMD:-rm -rf build/staging build/install build/lock.yaml}"
    local validate_cmd="${RETRY_VALIDATE_CMD:-}"
    local max_attempts=3
    local attempt=1
    while [ $attempt -le $max_attempts ]; do
        echo "Build attempt $attempt of $max_attempts..."
        if "$@"; then
            if [ -n "$validate_cmd" ]; then
                if ! eval "$validate_cmd"; then
                    echo "Post-build validation failed on attempt $attempt."
                    if [ $attempt -lt $max_attempts ]; then
                        echo "Cleaning and retrying in 5s..."
                        _retry_clean "$clean_cmd"
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
            echo "Build failed on attempt $attempt. Retrying in 30s..."
            _retry_clean "$clean_cmd"
            sleep 30
        else
            echo "Build failed after $max_attempts attempts."
            return 1
        fi
        attempt=$((attempt + 1))
    done
}
