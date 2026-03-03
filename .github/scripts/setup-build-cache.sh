#!/bin/bash
# Sets up a persistent build cache for self-hosted CI runners.
# Creates a symlink: ./build -> <cache_root>/<key>/build
#
# Each runner gets its own cache keyed by (cluster, device, interface, runner).
# This avoids cross-runner path issues entirely — CMake's absolute paths are
# always correct because the same runner always uses the same workspace path.
#
# Usage: source .github/scripts/setup-build-cache.sh <cluster> <device> <interface>

_cache_cluster="${1:?Usage: setup-build-cache.sh <cluster> <device> <interface>}"
_cache_device="${2:?}"
_cache_interface="${3:-none}"
_cache_runner="${RUNNER_NAME:?RUNNER_NAME not set}"

# Select cache root based on cluster (each HPC system has its own persistent storage).
case "$_cache_cluster" in
    phoenix)
        _cache_root="/storage/coda1/d-coc/0/sbryngelson3/.mfc-ci-cache" ;;
    frontier|frontier_amd)
        _cache_root="/lustre/orion/cfd154/scratch/sbryngelson/.mfc-ci-cache" ;;
    *)
        echo "=== Build Cache Setup ==="
        echo "  No cache root configured for cluster '$_cache_cluster' — skipping."
        echo "========================="
        return 0 2>/dev/null || exit 0 ;;
esac

_cache_key="${_cache_cluster}-${_cache_device}-${_cache_interface}-${_cache_runner}"
_cache_base="${_cache_root}/${_cache_key}/build"

# Check if the cache directory is healthy (readable, writable, no stale handles).
_cache_healthy() {
    local dir="$1"
    if ! ls "$dir" > /dev/null 2>&1; then
        echo "  Health check FAILED: cannot list $dir"
        return 1
    fi
    if [ -e "$dir/lock.yaml" ] && ! stat "$dir/lock.yaml" > /dev/null 2>&1; then
        echo "  Health check FAILED: cannot stat $dir/lock.yaml"
        return 1
    fi
    local probe="$dir/.nfs_probe.$$"
    if ! touch "$probe" 2>/dev/null || ! rm -f "$probe" 2>/dev/null; then
        echo "  Health check FAILED: cannot write/remove probe in $dir"
        rm -f "$probe" 2>/dev/null
        return 1
    fi
    return 0
}

# Nuclear recovery: rename stale cache out of the way and create a fresh one.
# Uses mv (operates on parent directory entry) which works even when children
# have stale file handles that prevent rm -rf from succeeding.
_cache_nuke() {
    local base="${1:-$_cache_base}"
    local stale_name="${base}.stale.$(date +%s)"
    echo "  NFS cache nuke: parking stale dir -> $stale_name"
    if mv "$base" "$stale_name" 2>/dev/null; then
        echo "  NFS cache nuke: renamed successfully"
    else
        echo "  NFS cache nuke: mv failed, trying rm -rf as fallback"
        rm -rf "$base" 2>/dev/null || true
    fi
    mkdir -p "$base"
    echo "  NFS cache nuke: fresh cache created at $base"
}

mkdir -p "$_cache_base"
_cache_dir="$(cd "$_cache_base" && pwd -P)"

echo "=== Build Cache Setup ==="
echo "  Cache key: $_cache_key"
echo "  Cache dir: $_cache_dir"

# Pre-flight: detect stale NFS handles before wasting a build attempt.
if ! _cache_healthy "$_cache_dir"; then
    echo "  Stale NFS cache detected — nuking and recreating."
    _cache_nuke "$_cache_base"
    _cache_dir="$(cd "$_cache_base" && pwd -P)"
fi

# Replace any existing build/ (real dir or stale symlink) with a symlink
# to our runner-specific cache directory.
# Use unlink for symlinks to avoid rm -rf following the link and deleting
# the shared cache contents (which another runner may be using).
if [ -L "build" ]; then
    unlink "build"
elif [ -e "build" ]; then
    rm -rf "build"
fi

ln -s "$_cache_dir" "build"

echo "  Symlink: build -> $_cache_dir"

# Garbage-collect stale cache dirs parked by _cache_nuke more than 7 days ago.
_cache_parent="$(dirname "$_cache_base")"
find "$_cache_parent" -maxdepth 1 -name "*.stale.*" -mtime +7 -exec rm -rf {} + 2>/dev/null || true

echo "========================="
