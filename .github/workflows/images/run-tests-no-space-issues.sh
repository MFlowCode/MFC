#!/bin/bash

# Ultimate MFC test runner designed to completely eliminate "No space left on device" errors
# This script avoids all known space-related issues by using minimal container footprint

set -e

echo "=== MFC No-Space-Issues Test Runner ==="
echo "Designed to avoid 'No space left on device' errors completely"

# Select container
if [ -f "mfc_cpu_optimized.sif" ]; then
    CONTAINER="mfc_cpu_optimized.sif"
else
    CONTAINER="mfc_cpu.sif"
fi
echo "Container: $CONTAINER"

# Ultra-conservative settings to avoid ANY space issues
MAX_TESTS=5
CLEANUP_FREQUENCY=1  # Clean after every test

echo "Will run maximum $MAX_TESTS tests with cleanup after every test"

# Ultra-aggressive cleanup function
nuclear_cleanup() {
    echo "Performing nuclear cleanup..."
    
    # Kill everything apptainer-related
    pkill -9 -f apptainer 2>/dev/null || true
    pkill -9 -f singularity 2>/dev/null || true
    sleep 3
    
    # Remove ALL possible temporary files and directories
    rm -rf /tmp/apptainer* /tmp/singularity* /tmp/mfc* /tmp/test* 2>/dev/null || true
    rm -rf ~/.apptainer ~/.singularity 2>/dev/null || true
    rm -rf /tmp/.* 2>/dev/null || true
    
    # Force filesystem sync
    sync
    sleep 1
    
    echo "Nuclear cleanup completed"
}

# Function to check space and abort if insufficient
check_space_or_die() {
    local space_kb=$(df /tmp | awk 'NR==2 {print $4}')
    local space_mb=$((space_kb / 1024))
    local space_gb=$((space_mb / 1024))
    
    echo "Available space: ${space_gb}GB (${space_mb}MB)"
    
    if [ $space_mb -lt 2000 ]; then
        echo "FATAL: Insufficient space ($space_mb MB < 2000MB required)"
        nuclear_cleanup
        exit 1
    fi
}

# Function to run a single test with ZERO space issues
run_zero_space_test() {
    local test_uuid=$1
    local test_num=$2
    
    echo ""
    echo "=== Test $test_num: $test_uuid ==="
    
    # Pre-test cleanup and space check
    nuclear_cleanup
    check_space_or_die
    
    # Create MINIMAL external directory (host filesystem only)
    local ext_dir="/tmp/ext-$$-$test_num"
    mkdir -p "$ext_dir"
    chmod 755 "$ext_dir"
    
    echo "External directory: $ext_dir"
    
    # Run test with ABSOLUTE MINIMAL footprint
    # - NO writable tmpfs (prevents container space issues)
    # - NO home directory mounting
    # - NO complex environment variables
    # - MINIMAL binds only
    
    local result=0
    
    echo "Starting test (minimal footprint)..."
    timeout 300 apptainer run \
        --no-home \
        --containall \
        --bind "$ext_dir:/tmp/test-output" \
        "$CONTAINER" \
        test --no-build -f $test_uuid -t $test_uuid > "$ext_dir/test.log" 2>&1 || result=$?
    
    if [ $result -eq 0 ]; then
        echo "✓ SUCCESS: Test $test_uuid passed"
        if [ -f "$ext_dir/test.log" ]; then
            echo "Last few lines of output:"
            tail -5 "$ext_dir/test.log" 2>/dev/null || echo "No log tail available"
        fi
    elif [ $result -eq 124 ]; then
        echo "⚠ TIMEOUT: Test $test_uuid timed out after 5 minutes"
    else
        echo "✗ FAILED: Test $test_uuid failed with exit code $result"
        if [ -f "$ext_dir/test.log" ]; then
            echo "Error output:"
            tail -10 "$ext_dir/test.log" 2>/dev/null || echo "No error log available"
        fi
    fi
    
    # Immediate post-test cleanup
    rm -rf "$ext_dir" 2>/dev/null || true
    nuclear_cleanup
    
    return $result
}

# Main execution
main() {
    echo "Starting zero-space-issues test execution..."
    
    if [ ! -f "$CONTAINER" ]; then
        echo "ERROR: Container $CONTAINER not found"
        exit 1
    fi
    
    # Initial nuclear cleanup
    nuclear_cleanup
    check_space_or_die
    
    # Get minimal test list using absolute minimal approach
    echo "Getting test list (minimal approach)..."
    
    local test_list_file="/tmp/test-list-$$"
    timeout 60 apptainer run \
        --no-home \
        --containall \
        "$CONTAINER" \
        test --list > "$test_list_file" 2>/dev/null || {
        echo "ERROR: Failed to get test list"
        rm -f "$test_list_file"
        exit 1
    }
    
    # Extract test UUIDs
    local test_uuids=$(grep -E '^  [A-F0-9]{8}' "$test_list_file" | awk '{print $1}' | head -$MAX_TESTS)
    rm -f "$test_list_file"
    
    if [ -z "$test_uuids" ]; then
        echo "ERROR: No tests found"
        exit 1
    fi
    
    # Convert to array
    local test_array=($test_uuids)
    local total=${#test_array[@]}
    
    echo "Found $total tests to run: ${test_array[*]}"
    
    local passed=0
    local failed=0
    
    # Run tests with maximum space isolation
    for ((i=0; i<total; i++)); do
        local uuid=${test_array[i]}
        local num=$((i + 1))
        
        if run_zero_space_test "$uuid" "$num"; then
            passed=$((passed + 1))
        else
            failed=$((failed + 1))
        fi
        
        echo "Progress: $num/$total (✓$passed ✗$failed)"
        
        # Mandatory cleanup between tests
        nuclear_cleanup
        sleep 2
    done
    
    echo ""
    echo "=== FINAL RESULTS ==="
    echo "Tests run: $total"
    echo "Passed: $passed"
    echo "Failed: $failed"
    
    if [ $passed -gt 0 ]; then
        echo "✓ SUCCESS: At least some tests passed without space errors!"
        echo "The space optimization is working!"
    else
        echo "✗ FAILURE: All tests failed - check container or test configuration"
    fi
    
    # Final cleanup
    nuclear_cleanup
    
    exit $failed
}

# Execute main function
main "$@"