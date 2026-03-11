#!/bin/bash
# Detects GPUs (NVIDIA or AMD), sets $ngpus and $gpu_ids.
# Usage: source .github/scripts/detect-gpus.sh

ngpus=0
gpu_ids=""
if command -v nvidia-smi &>/dev/null; then
    ngpus=$(nvidia-smi -L | wc -l)
    gpu_ids=$(seq -s ' ' 0 $((ngpus - 1)))
elif command -v rocm-smi &>/dev/null; then
    gpu_ids=$(rocm-smi --showid | awk '{print $1}' | grep -Eo '[0-9]+' | uniq | tr '\n' ' ')
    ngpus=$(echo "$gpu_ids" | wc -w)
fi
