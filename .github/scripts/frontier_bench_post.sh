#!/bin/bash
# Post-process all Frontier benchmark results after the SLURM job completes.
# Runs bench_diff for each config, comparing master vs PR YAML outputs.

set -euo pipefail

# Benchmark configs: cluster device interface flag
bench_configs=(
    "frontier:gpu:acc:f"
    "frontier:gpu:omp:f"
    "frontier_amd:gpu:omp:famd"
)

for cfg in "${bench_configs[@]}"; do
    IFS=':' read -r cluster device interface flag <<< "$cfg"
    pr_yaml="pr-${cluster}-${device}-${interface}/bench-${device}-${interface}.yaml"
    master_yaml="master-${cluster}-${device}-${interface}/bench-${device}-${interface}.yaml"

    echo "=========================================="
    echo "bench_diff: $cluster $device $interface"
    echo "  PR:     $pr_yaml"
    echo "  Master: $master_yaml"
    echo "=========================================="

    if [ ! -f "$pr_yaml" ]; then
        echo "ERROR: PR YAML not found: $pr_yaml"
        exit 1
    fi
    if [ ! -f "$master_yaml" ]; then
        echo "ERROR: Master YAML not found: $master_yaml"
        exit 1
    fi

    (cd pr && . ./mfc.sh load -c "$flag" -m g)
    (cd pr && ./mfc.sh bench_diff "../$master_yaml" "../$pr_yaml")
done
