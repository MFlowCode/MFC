#!/bin/bash
# Sets up a persistent build cache for self-hosted CI runners.
# Creates a symlink: ./build -> <scratch>/.mfc-ci-cache/<key>/build
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

_cache_key="${_cache_cluster}-${_cache_device}-${_cache_interface}-${_cache_runner}"
_cache_base="$HOME/scratch/.mfc-ci-cache/${_cache_key}/build"

mkdir -p "$_cache_base"
_cache_dir="$(cd "$_cache_base" && pwd -P)"

echo "=== Build Cache Setup ==="
echo "  Cache key: $_cache_key"
echo "  Cache dir: $_cache_dir"

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
echo "========================="
