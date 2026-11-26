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
    # A timed-out sbatch can still have queued the job; a second copy would then run in the same
    # directory and delete the first one's test inputs. Find it by name, directory and submit time.
    local name since
    name=$(printf '%s\n' "$script_contents" | sed -n 's/^#SBATCH -J //p' | head -n 1)
    since=$(date +%Y-%m-%dT%H:%M:%S)

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
            # standard: Frontier's SLURM_TIME_FORMAT=relative prints only the time for today's jobs
            job_id=$(SLURM_TIME_FORMAT=standard squeue -u "$USER" -h -n "$name" -t PD,R,CF,RQ -o '%i|%Z|%V' 2>/dev/null |
                awk -F'|' -v d="$PWD" -v t="$since" '$2 == d && $3 >= t { print $1; exit }')
            if [ -n "$name" ] && [ -n "$job_id" ]; then
                echo "sbatch timed out but queued job $job_id — using it." >&2
                echo "$job_id"
                return 0
            fi
        fi
        attempt=$((attempt + 1))
    done

    echo "sbatch failed after $max_attempts attempts. Last error: $last_output" >&2
    return 1
}

