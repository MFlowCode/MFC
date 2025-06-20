#!/bin/bash

# Ultra-minimal script to run just one MFC test with maximum space isolation
# Designed to completely avoid "No space left on device" errors

set -e

echo "=== Single Test Runner (No Space Issues) ==="

# Container selection
if [ -f "mfc_cpu_optimized.sif" ]; then
    CONTAINER="mfc_cpu_optimized.sif"
else
    CONTAINER="mfc_cpu.sif"
fi

echo "Using container: $CONTAINER"

# Clean up any existing processes/files
echo "Initial cleanup..."
pkill -f apptainer 2>/dev/null || true
rm -rf /tmp/single-test-* 2>/dev/null || true
sync

# Create isolated test environment
TEST_ID="single-test-$$-$(date +%s)"
TEST_BASE="/tmp/$TEST_ID"
mkdir -p "$TEST_BASE"/{work,cache,output,tests}
chmod 755 "$TEST_BASE" "$TEST_BASE"/{work,cache,output,tests}

echo "Test environment: $TEST_BASE"

# Check space
SPACE_MB=$(df /tmp | awk 'NR==2 {print int($4/1024)}')
echo "Available space: ${SPACE_MB}MB"

if [ $SPACE_MB -lt 500 ]; then
    echo "ERROR: Need at least 500MB, have ${SPACE_MB}MB"
    rm -rf "$TEST_BASE"
    exit 1
fi

# Get first test
echo "Getting test list..."
FIRST_TEST=$(apptainer run --no-home --containall \
    --bind "$TEST_BASE/work:/tmp/work" \
    --env TMPDIR="/tmp/work" \
    "$CONTAINER" \
    test --list | grep -E '^  [A-F0-9]{8}' | awk '{print $1}' | head -1)

if [ -z "$FIRST_TEST" ]; then
    echo "ERROR: No tests found"
    rm -rf "$TEST_BASE"
    exit 1
fi

echo "Running test: $FIRST_TEST"

# Run single test with complete isolation (NO writable tmpfs)
echo "Executing test..."
if apptainer run \
    --no-home \
    --containall \
    --bind "$TEST_BASE/tests:/opt/MFC/tests" \
    --bind "$TEST_BASE/cache:/tmp/cache" \
    --bind "$TEST_BASE/output:/tmp/output" \
    --bind "$TEST_BASE/work:/tmp/work" \
    --env TMPDIR="/tmp/work" \
    --env TEMP="/tmp/work" \
    --env TMP="/tmp/work" \
    --env MFC_TESTDIR="/opt/MFC/tests" \
    --env APPTAINER_CACHEDIR="/tmp/cache" \
    --env SINGULARITY_CACHEDIR="/tmp/cache" \
    "$CONTAINER" \
    test --no-build -f $FIRST_TEST -t $FIRST_TEST; then
    
    echo "✓ SUCCESS: Test $FIRST_TEST completed without space errors!"
    
else
    echo "✗ FAILED: Test $FIRST_TEST failed"
    echo "Checking what happened..."
    ls -la "$TEST_BASE"/ 2>/dev/null || echo "Test directory gone"
fi

# Final space check
FINAL_SPACE_MB=$(df /tmp | awk 'NR==2 {print int($4/1024)}')
SPACE_USED=$((SPACE_MB - FINAL_SPACE_MB))
echo "Space used: ${SPACE_USED}MB"

# Cleanup
echo "Cleaning up..."
rm -rf "$TEST_BASE" 2>/dev/null || true
pkill -f apptainer 2>/dev/null || true

echo "Test completed!"