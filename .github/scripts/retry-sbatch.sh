#!/bin/bash
# Provides retry_sbatch(): submits a job script string via sbatch with retries.
# Only retries on known transient SLURM/infrastructure errors (socket timeouts,
# connection failures). Hard failures (bad account, invalid partition, QOS
# violations) are not retried.
#
# Usage: source .github/scripts/retry-sbatch.sh
#        job_id=$(retry_sbatch "$script_contents")

retry_sbatch() {
    local script_contents="$1"
    local max_attempts=3
    local attempt=1
    local submit_output job_id last_output=""

    while [ $attempt -le $max_attempts ]; do
        echo "sbatch attempt $attempt of $max_attempts..." >&2
        submit_output=$(printf '%s\n' "$script_contents" | sbatch 2>&1) || true
        job_id=$(echo "$submit_output" | grep -oE 'Submitted batch job ([0-9]+)' | grep -oE '[0-9]+$')
        if [ -n "$job_id" ]; then
            echo "$job_id"
            return 0
        fi
        last_output="$submit_output"
        echo "sbatch failed: $submit_output" >&2
        if ! echo "$submit_output" | grep -qiE "timed out|connection refused|connection reset|temporarily unavailable"; then
            echo "Non-transient sbatch failure — not retrying." >&2
            return 1
        fi
        if [ $attempt -lt $max_attempts ]; then
            echo "Transient error — retrying in 30s..." >&2
            sleep 30
        fi
        attempt=$((attempt + 1))
    done

    echo "sbatch failed after $max_attempts attempts. Last error: $last_output" >&2
    return 1
}

