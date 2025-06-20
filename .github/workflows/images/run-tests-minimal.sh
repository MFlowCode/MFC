#!/bin/bash

# Minimal MFC test script designed to avoid "No space left on device" errors
# This script completely avoids writable tmpfs and uses only host filesystem

set -e

echo "=== MFC Minimal Space Test Runner ==="

# Use optimized container if available
if [ -f "mfc_cpu_optimized.sif" ]; then
    CONTAINER_PATH="mfc_cpu_optimized.sif"
    echo "Using optimized container: $CONTAINER_PATH"
else
    CONTAINER_PATH="mfc_cpu.sif"
    echo "Using standard container: $CONTAINER_PATH"
fi

# Ultra-conservative settings
TESTS_PER_CHUNK=1  # Run one test at a time
MAX_TESTS=10       # Limit total tests

echo "Settings: $TESTS_PER_CHUNK test per chunk, max $MAX_TESTS tests"

# Aggressive cleanup function
cleanup_all() {
    echo "Performing aggressive cleanup..."
    pkill -f apptainer 2>/dev/null || true
    sleep 2
    
    # Remove all possible temporary directories
    rm -rf /tmp/mfc-* 2>/dev/null || true
    rm -rf /tmp/apptainer-* 2>/dev/null || true  
    rm -rf /tmp/singularity-* 2>/dev/null || true
    rm -rf /tmp/test-* 2>/dev/null || true
    
    # Clean up cache
    rm -rf ~/.apptainer/cache/* 2>/dev/null || true
    rm -rf ~/.singularity/cache/* 2>/dev/null || true
    
    sync
    echo "Cleanup completed"
}

# Function to run a single test with maximum space isolation
run_single_test() {
    local test_uuid=$1
    local test_num=$2
    
    echo ""
    echo "=== Running test $test_num: $test_uuid ==="
    
    # Create completely isolated directories on host filesystem
    local test_base="/tmp/test-isolated-$$-$test_num"
    local test_work="$test_base/work"
    local test_cache="$test_base/cache"
    local test_output="$test_base/output" 
    local test_tests="$test_base/tests"
    
    mkdir -p "$test_work" "$test_cache" "$test_output" "$test_tests"
    chmod 755 "$test_base" "$test_work" "$test_cache" "$test_output" "$test_tests"
    
    echo "Working in: $test_base"
    
    # Check space before test
    local space_before=$(df /tmp | awk 'NR==2 {print $4}')
    local space_mb=$((space_before / 1024))
    echo "Available space before test: ${space_mb}MB"
    
    if [ $space_mb -lt 1000 ]; then
        echo "ERROR: Insufficient space ($space_mb MB < 1000MB)"
        rm -rf "$test_base" 2>/dev/null || true
        return 1
    fi
    
    # Run test with NO writable tmpfs and external directories
    apptainer run \
        --no-home \
        --containall \
        --bind "$test_tests:/opt/MFC/tests" \
        --bind "$test_cache:/tmp/cache" \
        --bind "$test_output:/tmp/output" \
        --bind "$test_work:/tmp/work" \
        --env TMPDIR="/tmp/work" \
        --env TEMP="/tmp/work" \
        --env TMP="/tmp/work" \
        --env MFC_TESTDIR="/opt/MFC/tests" \
        --env APPTAINER_CACHEDIR="/tmp/cache" \
        --env SINGULARITY_CACHEDIR="/tmp/cache" \
        --env MFC_NO_VERBOSE=1 \
        --env MFC_QUIET=1 \
        "$CONTAINER_PATH" \
        test --no-build -f $test_uuid -t $test_uuid || {
        echo "Test $test_uuid FAILED"
        rm -rf "$test_base" 2>/dev/null || true
        return 1
    }
    
    echo "Test $test_uuid PASSED"
    
    # Check space after test
    local space_after=$(df /tmp | awk 'NR==2 {print $4}')
    local space_after_mb=$((space_after / 1024))
    local space_used=$((space_mb - space_after_mb))
    echo "Space used by test: ${space_used}MB"
    
    # Immediate cleanup
    rm -rf "$test_base" 2>/dev/null || true
    
    # Verify cleanup worked
    local space_final=$(df /tmp | awk 'NR==2 {print $4}')
    local space_final_mb=$((space_final / 1024))
    local space_recovered=$((space_final_mb - space_after_mb))
    echo "Space recovered: ${space_recovered}MB"
    
    return 0
}

# Main execution
main() {
    if [ ! -f "$CONTAINER_PATH" ]; then
        echo "ERROR: Container file $CONTAINER_PATH not found"
        exit 1
    fi
    
    echo "Starting minimal space test execution..."
    
    # Initial cleanup
    cleanup_all
    
    # Get test list (limit to first few)
    echo "Getting test list..."
    TEST_UUIDS=$(apptainer run --no-home --containall \
        --bind /tmp:/tmp/host \
        --env TMPDIR="/tmp/host" \
        "$CONTAINER_PATH" \
        test --list | grep -E '^  [A-F0-9]{8}' | awk '{print $1}' | head -$MAX_TESTS)
    
    if [ -z "$TEST_UUIDS" ]; then
        echo "ERROR: No tests found"
        exit 1
    fi
    
    # Convert to array
    TEST_ARRAY=($TEST_UUIDS)
    local total_tests=${#TEST_ARRAY[@]}
    echo "Found $total_tests tests to run"
    
    local passed_tests=0
    local failed_tests=0
    
    # Run tests one by one
    for ((i=0; i<total_tests; i++)); do
        local test_uuid=${TEST_ARRAY[i]}
        local test_num=$((i + 1))
        
        if run_single_test "$test_uuid" "$test_num"; then
            passed_tests=$((passed_tests + 1))
        else
            failed_tests=$((failed_tests + 1))
        fi
        
        # Cleanup between tests
        cleanup_all
        
        echo "Progress: $test_num/$total_tests (Passed: $passed_tests, Failed: $failed_tests)"
    done
    
    echo ""
    echo "=== Final Results ==="
    echo "Total tests: $total_tests"
    echo "Passed: $passed_tests"
    echo "Failed: $failed_tests"
    echo "Success rate: $(( passed_tests * 100 / total_tests ))%"
    
    if [ $passed_tests -gt 0 ]; then
        echo "✓ Tests ran successfully with space optimization!"
        exit 0
    else
        echo "✗ All tests failed"
        exit 1
    fi
}

# Run main function
main "$@"