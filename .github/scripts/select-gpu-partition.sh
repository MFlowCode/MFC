#!/bin/bash
# Select the best available Phoenix GPU partition using sinfo.
# Sources into caller: exports SELECTED_GPU_PARTITION.
#
# Priority order prefers partitions most likely to have availability.
# V100 is last due to slower performance near the test time limit.
# Falls back to gpu-l40s if no partition meets the idle node threshold.
# RTX 6000 nodes are excluded (too slow for the test suite time limit).
#
# Optional: set GPU_PARTITION_MIN_NODES before sourcing to require a minimum
# number of idle/mix nodes (e.g. GPU_PARTITION_MIN_NODES=2 for parallel bench jobs).
#
# Usage: source .github/scripts/select-gpu-partition.sh

_GPU_PARTITION_PRIORITY="gpu-l40s gpu-h200 gpu-h100 gpu-a100 gpu-v100"
_GPU_PARTITION_FALLBACK="gpu-l40s"
_GPU_PARTITION_MIN_NODES="${GPU_PARTITION_MIN_NODES:-1}"

SELECTED_GPU_PARTITION=""
for _part in $_GPU_PARTITION_PRIORITY; do
    _idle=$(sinfo -p "$_part" --noheader -o "%t" 2>/dev/null | grep -cE "^(idle|mix)" || true)
    if [ "${_idle:-0}" -ge "$_GPU_PARTITION_MIN_NODES" ]; then
        SELECTED_GPU_PARTITION="$_part"
        echo "Selected GPU partition: $SELECTED_GPU_PARTITION ($_idle idle/mix nodes)"
        break
    fi
done

if [ -z "$SELECTED_GPU_PARTITION" ]; then
    echo "WARNING: No idle GPU partition found; falling back to $_GPU_PARTITION_FALLBACK (may queue)"
    SELECTED_GPU_PARTITION="$_GPU_PARTITION_FALLBACK"
fi

export SELECTED_GPU_PARTITION
unset _GPU_PARTITION_PRIORITY _GPU_PARTITION_FALLBACK _GPU_PARTITION_MIN_NODES _part _idle
