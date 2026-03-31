#!/bin/bash
# Simple build script that captures compiler warnings.
# Builds with --reldebug to enable warning flags, then dumps warnings.
# No retries, no validation, no tests — just build and report.

set -uo pipefail

source .github/scripts/gpu-opts.sh

# Phoenix needs TMPDIR setup
if [ "$job_cluster" = "phoenix" ]; then
    tmpbuild=/storage/project/r-sbryngelson3-0/sbryngelson3/mytmp_build
    currentdir=$tmpbuild/run-$(( RANDOM % 9000 ))
    mkdir -p $tmpbuild $currentdir
    export TMPDIR=$currentdir
    trap 'rm -rf "$currentdir" || true' EXIT
fi

# Clean stale builds on Phoenix
if [ "$job_cluster" = "phoenix" ]; then
    ./mfc.sh clean 2>/dev/null || true
fi

# Build with reldebug (enables warning flags)
./mfc.sh build -v --reldebug -j 8 $gpu_opts 2>&1 || true

echo ""
echo "=== COMPILER WARNINGS SUMMARY ==="
echo ""

# Rebuild each MFC target single-threaded to get clean warnings with file/line
for staging_dir in build/staging/*/; do
    [ -f "$staging_dir/CMakeCache.txt" ] || continue
    case "$staging_dir" in *fftw*|*hdf5*|*silo*|*lapack*) continue ;; esac
    target=$(find "$staging_dir/CMakeFiles/" -maxdepth 1 -name '*.dir' -type d 2>/dev/null | head -1 | xargs basename 2>/dev/null | sed 's/\.dir//')
    [ -z "$target" ] && continue

    echo "--- $target ---"
    cd "$staging_dir"
    find fypp/ -name '*.f90' -exec touch {} + 2>/dev/null
    make "$target" -j1 2>&1 | grep -E "Warning:|warning:" -B2 || echo "(no warnings)"
    cd - > /dev/null
    echo ""
done
