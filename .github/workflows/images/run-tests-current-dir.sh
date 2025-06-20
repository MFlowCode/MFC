#!/bin/bash

# Script to run MFC tests using current directory for test output
# This avoids external directory binding issues

set -e

# Use optimized container if available, otherwise fall back to regular
if [ -f "mfc_cpu_optimized.sif" ]; then
    CONTAINER_PATH="mfc_cpu_optimized.sif"
    echo "Using optimized container: $CONTAINER_PATH"
else
    CONTAINER_PATH="mfc_cpu.sif"
    echo "Using standard container: $CONTAINER_PATH"
fi

CHUNK_SIZE=18  # Reduced chunk size for better space management
MEMORY="16G"   # Reduced memory allocation to avoid space issues

# Get current directory
CURRENT_DIR=$(pwd)
TEST_OUTPUT_DIR="$CURRENT_DIR/mfc-test-output"
mkdir -p "$TEST_OUTPUT_DIR"

# Enhanced cleanup function for current directory approach
cleanup_between_chunks() {
    echo "Performing enhanced cleanup for current directory approach..."
    
    # Stop any lingering apptainer processes
    pkill -f apptainer 2>/dev/null || true
    sleep 1
    
    # Clean up any test output directories
    find /tmp -name "mfc-*" -type d -exec rm -rf {} + 2>/dev/null || true
    find /tmp -name "apptainer-*" -type d -exec rm -rf {} + 2>/dev/null || true
    find /tmp -name "singularity-*" -type d -exec rm -rf {} + 2>/dev/null || true
    
    # Clean up common test output file types
    find /tmp -name "*.dat" -type f -delete 2>/dev/null || true
    find /tmp -name "*.h5" -type f -delete 2>/dev/null || true
    find /tmp -name "*.hdf5" -type f -delete 2>/dev/null || true
    find /tmp -name "*.vtk" -type f -delete 2>/dev/null || true
    find /tmp -name "*.silo" -type f -delete 2>/dev/null || true
    find /tmp -name "*.log" -type f -delete 2>/dev/null || true
    find /tmp -name "*.out" -type f -delete 2>/dev/null || true
    find /tmp -name "*.err" -type f -delete 2>/dev/null || true
    find /tmp -name "*.tmp" -type f -delete 2>/dev/null || true
    find /tmp -name "core.*" -type f -delete 2>/dev/null || true
    
    # Clean up test output directory
    rm -rf "$TEST_OUTPUT_DIR"/* 2>/dev/null || true
    mkdir -p "$TEST_OUTPUT_DIR" 2>/dev/null || true
    
    # Clean up cache directories
    rm -rf ~/.apptainer/cache/* 2>/dev/null || true
    rm -rf ~/.singularity/cache/* 2>/dev/null || true
    rm -rf /tmp/.apptainer* 2>/dev/null || true
    rm -rf /tmp/.singularity* 2>/dev/null || true
    
    # Force sync
    sync
    echo "Enhanced cleanup completed"
}

# Function to check available space
check_space() {
    local available_space=$(df /tmp | awk 'NR==2 {print $4}')
    local available_mb=$((available_space / 1024))
    echo "Available space in /tmp: ${available_mb}MB"
    
    if [ $available_mb -lt 1000 ]; then
        echo "Warning: Low space detected, performing aggressive cleanup..."
        cleanup_between_chunks
    fi
}

# Function to run a chunk of tests using current directory
run_test_chunk() {
    local start_idx=$1
    local end_idx=$2
    
    echo "Running tests chunk: $start_idx to $end_idx"
    
    # Check space before running
    check_space
    
    # Create a temporary working directory for this chunk
    local chunk_work_dir="/tmp/mfc-chunk-$$"
    mkdir -p "$chunk_work_dir"
    chmod 755 "$chunk_work_dir"
    
    apptainer run \
        --writable-tmpfs \
        --memory "$MEMORY" \
        --tmpdir "$chunk_work_dir" \
        --bind "$CURRENT_DIR:/opt/MFC" \
        --bind /tmp:/tmp \
        --env TMPDIR="$chunk_work_dir" \
        --env TEMP="$chunk_work_dir" \
        --env TMP="$chunk_work_dir" \
        --env MFC_TESTDIR="$TEST_OUTPUT_DIR" \
        --env APPTAINER_CACHEDIR="$chunk_work_dir/cache" \
        --env SINGULARITY_CACHEDIR="$chunk_work_dir/cache" \
        "$CONTAINER_PATH" \
        test --no-build -f $start_idx -t $end_idx || {
        echo "Test chunk $start_idx-$end_idx failed, continuing..."
        rm -rf "$chunk_work_dir" 2>/dev/null || true
        return 1
    }
    
    # Clean up chunk working directory
    rm -rf "$chunk_work_dir" 2>/dev/null || true
    
    cleanup_between_chunks
}

# Function to run all tests in smaller chunks
run_chunked_tests() {
    echo "Getting test list..."
    
    # Get all test UUIDs
    TEST_UUIDS=$(apptainer run \
        --writable-tmpfs \
        --memory "$MEMORY" \
        --tmpdir /tmp \
        --bind "$CURRENT_DIR:/opt/MFC" \
        --bind /tmp:/tmp \
        --env TMPDIR=/tmp \
        --env TEMP=/tmp \
        --env TMP=/tmp \
        "$CONTAINER_PATH" \
        test --list | grep -E '^  [A-F0-9]{8}' | awk '{print $1}' | head -100)  # Limit to first 100 tests
    
    echo "Found $(echo "$TEST_UUIDS" | wc -l) tests to run"
    
    # Convert to array
    TEST_ARRAY=($TEST_UUIDS)
    
    local failed_chunks=0
    local total_chunks=0
    
    # Run tests in chunks
    for ((i=0; i<${#TEST_ARRAY[@]}; i+=CHUNK_SIZE)); do
        total_chunks=$((total_chunks + 1))
        local chunk_start=${TEST_ARRAY[i]}
        local chunk_end_idx=$((i + CHUNK_SIZE - 1))
        
        if [ $chunk_end_idx -ge ${#TEST_ARRAY[@]} ]; then
            chunk_end_idx=$((${#TEST_ARRAY[@]} - 1))
        fi
        
        local chunk_end=${TEST_ARRAY[chunk_end_idx]}
        
        echo "Running chunk $total_chunks: tests $(($i + 1)) to $(($chunk_end_idx + 1))"
        
        if ! run_test_chunk "$chunk_start" "$chunk_end"; then
            failed_chunks=$((failed_chunks + 1))
        fi
        
        echo "Completed chunk $total_chunks"
        echo "---"
    done
    
    echo "Test summary:"
    echo "  Total chunks: $total_chunks"
    echo "  Failed chunks: $failed_chunks"
    echo "  Success rate: $(( (total_chunks - failed_chunks) * 100 / total_chunks ))%"
}

# Main execution
main() {
    if [ ! -f "$CONTAINER_PATH" ]; then
        echo "Error: Container file $CONTAINER_PATH not found"
        exit 1
    fi
    
    echo "Starting chunked test execution with current directory mounting..."
    echo "Container: $CONTAINER_PATH"
    echo "Memory allocation: $MEMORY"
    echo "Chunk size: $CHUNK_SIZE tests"
    echo "Current directory: $CURRENT_DIR"
    echo "Test output directory: $TEST_OUTPUT_DIR"
    echo "---"
    
    # Initial cleanup
    cleanup_between_chunks
    
    # Check initial space
    check_space
    
    run_chunked_tests
}

# Check if script is being run directly
if [ "$0" = "${BASH_SOURCE[0]}" ]; then
    main "$@"
fi 