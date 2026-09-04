#!/bin/bash
# Select the best available Phoenix GPU partition using sinfo.
# Sources into caller: exports SELECTED_GPU_PARTITION.
#
# Priority order prefers partitions most likely to have availability.
# V100 is last due to slower performance near the test time limit.
# Falls back to gpu-a100 if no partition meets the idle node threshold.
#
# gpu-l40s is out of rotation: it has been failing jobs for weeks, and it was
# also the partition CI kept selecting and then starving on.
#
# Only fully idle nodes count. A "mix" node is partially allocated and may have
# no free GPU, so counting it overstates availability: run 33553417354 picked
# gpu-l40s on "1 idle/mix nodes", then sat in the queue until the 3-5.5h job
# timeout without ever starting.
#
# The match is anchored at both ends. sinfo's %t suffixes a state to flag it --
# "*" not responding, "$" reserved for maintenance, "~" powered down -- and this
# cluster does emit them (drain*, down*, alloc$, drain$ are all live right now).
# A bare "^idle" would count idle* and idle$ as available and starve the job on
# nodes that cannot take it.
# RTX 6000 nodes are excluded (too slow for the test suite time limit).
#
# Optional: set GPU_PARTITION_MIN_NODES before sourcing to require a minimum
# number of idle nodes (e.g. GPU_PARTITION_MIN_NODES=2 for parallel bench jobs).
#
# Usage: source .github/scripts/select-gpu-partition.sh

_GPU_PARTITION_PRIORITY="gpu-h200 gpu-h100 gpu-a100 gpu-v100"
_GPU_PARTITION_FALLBACK="gpu-a100"
_GPU_PARTITION_MIN_NODES="${GPU_PARTITION_MIN_NODES:-1}"

SELECTED_GPU_PARTITION=""
for _part in $_GPU_PARTITION_PRIORITY; do
    _idle=$(sinfo -p "$_part" --noheader -o "%t" 2>/dev/null | grep -cE "^idle$" || true)
    if [ "${_idle:-0}" -ge "$_GPU_PARTITION_MIN_NODES" ]; then
        SELECTED_GPU_PARTITION="$_part"
        echo "Selected GPU partition: $SELECTED_GPU_PARTITION ($_idle idle nodes)"
        break
    fi
done

if [ -z "$SELECTED_GPU_PARTITION" ]; then
    echo "WARNING: No idle GPU partition found; falling back to $_GPU_PARTITION_FALLBACK (may queue)"
    SELECTED_GPU_PARTITION="$_GPU_PARTITION_FALLBACK"
fi

export SELECTED_GPU_PARTITION
unset _GPU_PARTITION_PRIORITY _GPU_PARTITION_FALLBACK _GPU_PARTITION_MIN_NODES _part _idle
