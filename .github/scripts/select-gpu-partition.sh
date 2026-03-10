#!/bin/bash
# Select the best available Phoenix GPU partition using sinfo.
# Sources into caller: exports SELECTED_GPU_PARTITION.
#
# Priority order prefers smaller/older nodes to leave modern GPUs free
# for production workloads. Falls back to gpu-rtx6000 if nothing is idle.
#
# Usage: source .github/scripts/select-gpu-partition.sh

_GPU_PARTITION_PRIORITY="gpu-rtx6000 gpu-l40s gpu-v100 gpu-h200 gpu-h100 gpu-a100"
_GPU_PARTITION_FALLBACK="gpu-rtx6000"

SELECTED_GPU_PARTITION=""
for _part in $_GPU_PARTITION_PRIORITY; do
    _idle=$(sinfo -p "$_part" --noheader -o "%t" 2>/dev/null | grep -cE "^(idle|mix)" || true)
    if [ "${_idle:-0}" -gt 0 ]; then
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
unset _GPU_PARTITION_PRIORITY _GPU_PARTITION_FALLBACK _part _idle
