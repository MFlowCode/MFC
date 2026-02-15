#!/bin/bash
# Sets up a persistent build cache for self-hosted CI runners.
# Creates a symlink: ./build -> <resolved scratch path>/.mfc-ci-cache/<key>/build
#
# This ensures that every run of the same config (cluster/device/interface) finds
# cached build artifacts regardless of which runner instance picks up the job.
#
# Concurrent safety: uses flock to serialize access per cache directory. If
# multiple PRs trigger the same config simultaneously, the second job waits
# for the first to finish (up to 1 hour), then gets a warm cache. If the lock
# times out, falls back to a local build (same as no caching).
#
# Usage: source .github/scripts/setup-build-cache.sh <cluster> <device> <interface>

_cache_cluster="${1:?Usage: setup-build-cache.sh <cluster> <device> <interface>}"
_cache_device="${2:?}"
_cache_interface="${3:-none}"

_cache_key="${_cache_cluster}-${_cache_device}-${_cache_interface}"
_cache_base="$HOME/scratch/.mfc-ci-cache/${_cache_key}/build"

# Create the cache dir, then resolve to a physical path (no symlinks).
# $HOME/scratch is typically a symlink to a scratch filesystem — resolving
# it ensures the build symlink target remains valid even if intermediate
# symlinks change.
mkdir -p "$_cache_base"
_cache_dir="$(cd "$_cache_base" && pwd -P)"

echo "=== Build Cache Setup ==="
echo "  Cache key: $_cache_key"
echo "  Cache dir: $_cache_dir"

# Acquire an exclusive lock on the cache directory to prevent concurrent
# builds from corrupting it. The lock is fd-based (flock on fd 9), so it
# auto-releases when the calling process exits — no stale locks.
#
# Timeout: 1 hour. If another build holds the lock, we wait. This is fine
# because the waiting job will get a warm cache when it finally acquires.
# If the lock can't be acquired after 1 hour, something is wrong — fall
# back to a local build in the workspace.
_cache_locked=false
_lock_file="$_cache_dir/.cache.lock"
exec 9>"$_lock_file"
echo "  Acquiring cache lock..."
if flock --timeout 3600 9; then
    _cache_locked=true
    echo "  Cache lock acquired"
else
    echo "  WARNING: Cache lock timeout (1h), building locally without cache"
    exec 9>&-
    # Remove any existing symlink to the shared cache so we don't write
    # into it without the lock. Then create a real local directory.
    if [ -L "build" ]; then
        rm -f "build"
    fi
    mkdir -p "build"
    echo "========================="
    return 0 2>/dev/null || true
fi

# If build/ exists (real dir or stale symlink), remove it.
# rm -rf on a symlink removes the symlink, not the target — cache is safe.
if [ -e "build" ] || [ -L "build" ]; then
    rm -rf "build"
fi

ln -s "$_cache_dir" "build"

# Handle cross-runner workspace path changes.
# CMakeCache.txt stores absolute paths from whichever runner instance
# originally configured the build. If we're on a different runner, sed-replace
# the old workspace path with the current one so CMake can do incremental builds.
_workspace_marker="$_cache_dir/.workspace_path"
if [ -f "$_workspace_marker" ]; then
    _old_workspace=$(cat "$_workspace_marker")
    if [ "$_old_workspace" != "$(pwd)" ]; then
        echo "  Workspace path changed: $_old_workspace -> $(pwd)"
        echo "  Updating cached paths..."
        # Update CMake build files in staging/
        find "$_cache_dir/staging" -type f \
            \( -name "CMakeCache.txt" -o -name "*.cmake" \
               -o -name "*.make" -o -name "Makefile" \
               -o -name "build.ninja" \) \
            -exec sed -i "s|${_old_workspace}|$(pwd)|g" {} + 2>/dev/null || true
        # Compiled binaries have stale paths baked in — delete install/
        # so CMake rebuilds and re-installs them with correct paths.
        echo "  Clearing install/ to force rebuild of binaries..."
        rm -rf "$_cache_dir/install"
    fi
fi
echo "$(pwd)" > "$_workspace_marker"

echo "  Symlink: build -> $_cache_dir"
echo "========================="
