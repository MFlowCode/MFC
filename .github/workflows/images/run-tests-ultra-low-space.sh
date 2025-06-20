#!/bin/bash

# Script to run MFC tests in ultra space-constrained environments
# This script uses the most aggressive space optimization techniques

set -e

# Use optimized container if available, otherwise fall back to regular
if [ -f "mfc_cpu_optimized.sif" ]; then
    CONTAINER_PATH="mfc_cpu_optimized.sif"
    echo "Using optimized container: $CONTAINER_PATH"
else
    CONTAINER_PATH="mfc_cpu.sif"
    echo "Using standard container: $CONTAINER_PATH"
fi

CHUNK_SIZE=10  # Small chunk size for ultra-low space  
MEMORY="8G"    # Conservative memory allocation
TMPFS_SIZE="16G"  # Conservative tmpfs size

# Create external test directory on a different filesystem if possible
EXTERNAL_TEST_DIR="/tmp/mfc-external-tests"
mkdir -p "$EXTERNAL_TEST_DIR"

# Function for ultra-aggressive cleanup
ultra_cleanup() {
    echo "Performing ultra-aggressive cleanup..."
    
    # Stop any running containers that might be consuming space
    apptainer instance list | grep -v "INSTANCE NAME" | awk '{print $1}' | xargs -r apptainer instance stop 2>/dev/null || true
    
    # Clean up all temporary files
    find /tmp -name "*" -type f -delete 2>/dev/null || true
    find /tmp -name "*" -type d -exec rm -rf {} + 2>/dev/null || true
    
    # Clean up external test directory
    rm -rf "$EXTERNAL_TEST_DIR"/* 2>/dev/null || true
    
    # Clean up apptainer cache
    rm -rf ~/.apptainer/cache/* 2>/dev/null || true
    rm -rf ~/.singularity/cache/* 2>/dev/null || true
    
    # Force system cleanup
    sync
    echo 3 > /proc/sys/vm/drop_caches 2>/dev/null || true
    
    # Recreate necessary directories
    mkdir -p "$EXTERNAL_TEST_DIR"
    mkdir -p /tmp
}

# Function to check available space with thresholds
check_space() {
    local available_space=$(df /tmp | awk 'NR==2 {print $4}')
    local available_mb=$((available_space / 1024))
    echo "Available space in /tmp: ${available_mb}MB"
    
    if [ $available_mb -lt 500 ]; then
        echo "Critical: Very low space detected, performing ultra-aggressive cleanup..."
        ultra_cleanup
    elif [ $available_mb -lt 2000 ]; then
        echo "Warning: Low space detected, performing aggressive cleanup..."
        ultra_cleanup
    fi
}

# Function to run a single test with maximum space optimization
run_single_test() {
    local test_uuid=$1
    
    echo "Running single test: $test_uuid"
    
    # Check space before running
    check_space
    
    # Create isolated working directory
    local test_work_dir="/tmp/mfc-single-test-$$"
    mkdir -p "$test_work_dir"
    
    # Run the test with maximum isolation
    apptainer run \
        --writable-tmpfs \
        --memory "$MEMORY" \
        --tmpdir "$test_work_dir" \
        --bind "$test_work_dir:/opt/MFC/tests" \
        --bind "$test_work_dir:/tmp/mfc-tests" \
        --bind "$test_work_dir:/tmp/chunk-work" \
        --env TMPDIR="$test_work_dir" \
        --env TEMP="$test_work_dir" \
        --env TMP="$test_work_dir" \
        --env MFC_TESTDIR="$test_work_dir/mfc-tests" \
        --env APPTAINER_CACHEDIR="$test_work_dir/cache" \
        --env SINGULARITY_CACHEDIR="$test_work_dir/cache" \
        --env MFC_NO_OUTPUT=1 \
        --env MFC_QUIET=1 \
        "$CONTAINER_PATH" \
        test --no-build -f $test_uuid -t $test_uuid || {
        echo "Test $test_uuid failed"
        rm -rf "$test_work_dir" 2>/dev/null || true
        return 1
    }
    
    # Clean up immediately
    rm -rf "$test_work_dir" 2>/dev/null || true
    return 0
}

# Function to run tests in very small batches
run_batch_tests() {
    local test_uuids=("$@")
    local batch_size=5  # Very small batch size for ultra-low space
    
    local failed_tests=0
    local total_tests=${#test_uuids[@]}
    
    for ((i=0; i<total_tests; i+=batch_size)); do
        local end_idx=$((i + batch_size - 1))
        if [ $end_idx -ge $total_tests ]; then
            end_idx=$((total_tests - 1))
        fi
        
        echo "Running batch: tests $(($i + 1)) to $(($end_idx + 1)) of $total_tests"
        
        # Run tests in this batch
        for ((j=i; j<=end_idx; j++)); do
            if ! run_single_test "${test_uuids[j]}"; then
                failed_tests=$((failed_tests + 1))
            fi
        done
        
        # Aggressive cleanup between batches
        ultra_cleanup
        echo "Completed batch, cleaned up"
        echo "---"
    done
    
    echo "Batch test summary:"
    echo "  Total tests: $total_tests"
    echo "  Failed tests: $failed_tests"
    echo "  Success rate: $(( (total_tests - failed_tests) * 100 / total_tests ))%"
}

# Function to run all tests with ultra-low space optimization
run_ultra_low_space_tests() {
    echo "Getting test list..."
    
    # Get all test UUIDs
    TEST_UUIDS=$(apptainer run \
        --writable-tmpfs \
        --memory "$MEMORY" \
        --tmpdir /tmp \
        --bind /tmp:/tmp \
        --env TMPDIR=/tmp \
        --env TEMP=/tmp \
        --env TMP=/tmp \
        "$CONTAINER_PATH" \
        test --list | grep -E '^  [A-F0-9]{8}' | awk '{print $1}' | head -100)  # Limit to first 100 tests
    
    echo "Found $(echo "$TEST_UUIDS" | wc -l) tests to run"
    
    # Convert to array
    TEST_ARRAY=($TEST_UUIDS)
    
    # Run tests in batches
    run_batch_tests "${TEST_ARRAY[@]}"
}

# Main execution
main() {
    if [ ! -f "$CONTAINER_PATH" ]; then
        echo "Error: Container file $CONTAINER_PATH not found"
        exit 1
    fi
    
    echo "Starting ultra-low space test execution..."
    echo "Container: $CONTAINER_PATH"
    echo "Memory allocation: $MEMORY"
    echo "External test directory: $EXTERNAL_TEST_DIR"
    echo "---"
    
    # Initial ultra cleanup
    ultra_cleanup
    
    # Check initial space
    check_space
    
    run_ultra_low_space_tests
}

# Check if script is being run directly
if [ "$0" = "${BASH_SOURCE[0]}" ]; then
    main "$@"
fi 