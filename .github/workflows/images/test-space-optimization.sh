#!/bin/bash

# Quick test script to validate space optimizations
set -e

echo "=== MFC Space Optimization Test ==="
echo "Testing space management with 2 tests only"
echo ""

# Use optimized container if available
if [ -f "mfc_cpu_optimized.sif" ]; then
    CONTAINER_PATH="mfc_cpu_optimized.sif"
    echo "Using optimized container: $CONTAINER_PATH"
else
    CONTAINER_PATH="mfc_cpu.sif"
    echo "Using standard container: $CONTAINER_PATH"
fi

# Check initial space
initial_space=$(df /tmp | awk 'NR==2 {print $4}')
initial_mb=$((initial_space / 1024))
echo "Initial space: ${initial_mb}MB"

# Enhanced cleanup function
cleanup_test() {
    echo "Performing cleanup..."
    pkill -f apptainer 2>/dev/null || true
    sleep 1
    rm -rf /tmp/mfc-test-* 2>/dev/null || true
    rm -rf /tmp/apptainer-* 2>/dev/null || true
    rm -rf /tmp/singularity-* 2>/dev/null || true
    sync
}

# Get 2 test UUIDs
echo "Getting test list..."
TEST_UUIDS=$(apptainer run --writable-tmpfs --memory 8G --tmpdir /tmp \
    "$CONTAINER_PATH" test --list | grep -E '^  [A-F0-9]{8}' | awk '{print $1}' | head -2)

TEST_ARRAY=($TEST_UUIDS)
echo "Found ${#TEST_ARRAY[@]} tests to run: ${TEST_ARRAY[0]} ${TEST_ARRAY[1]}"

# Run tests with space isolation
test_work_dir="/tmp/mfc-test-$$"
mkdir -p "$test_work_dir"/{cache,tests,output}
chmod 755 "$test_work_dir" "$test_work_dir"/{cache,tests,output}

echo ""
echo "Running tests with isolated working directory: $test_work_dir"

# Monitor space before
before_space=$(df /tmp | awk 'NR==2 {print $4}')
before_mb=$((before_space / 1024))
echo "Space before test: ${before_mb}MB"

# Run the tests
echo "Executing tests..."
apptainer run \
    --writable-tmpfs \
    --memory 8G \
    --tmpdir "$test_work_dir" \
    --bind "$test_work_dir:/tmp/mfc-isolated" \
    --env TMPDIR="/tmp/mfc-isolated" \
    --env TEMP="/tmp/mfc-isolated" \
    --env TMP="/tmp/mfc-isolated" \
    --env MFC_TESTDIR="/tmp/mfc-isolated/tests" \
    --env APPTAINER_CACHEDIR="/tmp/mfc-isolated/cache" \
    --env SINGULARITY_CACHEDIR="/tmp/mfc-isolated/cache" \
    --env MFC_NO_VERBOSE=1 \
    --env MFC_QUIET=1 \
    "$CONTAINER_PATH" \
    test --no-build -f ${TEST_ARRAY[0]} -t ${TEST_ARRAY[1]} || {
    echo "Tests failed"
    rm -rf "$test_work_dir" 2>/dev/null || true
    cleanup_test
    exit 1
}

echo "Tests completed successfully!"

# Monitor space after
after_space=$(df /tmp | awk 'NR==2 {print $4}')
after_mb=$((after_space / 1024))
echo "Space after test: ${after_mb}MB"

# Calculate space used
space_used=$((before_mb - after_mb))
echo "Space used during test: ${space_used}MB"

# Clean up test directory
echo "Cleaning up test directory..."
rm -rf "$test_work_dir" 2>/dev/null || true

# Final cleanup
cleanup_test

# Final space check
final_space=$(df /tmp | awk 'NR==2 {print $4}')
final_mb=$((final_space / 1024))
echo "Final space: ${final_mb}MB"

# Calculate cleanup effectiveness
recovered_space=$((final_mb - after_mb))
echo "Space recovered by cleanup: ${recovered_space}MB"

echo ""
echo "=== Test Summary ==="
echo "Initial space:  ${initial_mb}MB"
echo "Used by tests:  ${space_used}MB"
echo "Recovered:      ${recovered_space}MB"
echo "Final space:    ${final_mb}MB"

if [ $recovered_space -ge 0 ]; then
    echo "✓ Space optimization successful!"
    exit 0
else
    echo "⚠ Space recovery was incomplete"
    exit 1
fi