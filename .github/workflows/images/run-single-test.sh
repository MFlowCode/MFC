#!/bin/bash

# Script to run a single MFC test with optimized space usage
# Usage: ./run-single-test.sh [TEST_UUID]

set -e

CONTAINER_PATH="mfc_cpu.sif"
TEST_UUID=${1:-"D79C3E6F"}  # Default to first test
MEMORY="32G"

echo "Running single test: $TEST_UUID"
echo "Memory allocation: $MEMORY"

# Create a larger tmpfs and run the test
apptainer run \
    --writable-tmpfs \
    --memory "$MEMORY" \
    --tmpdir /tmp \
    --bind /tmp:/tmp \
    --env TMPDIR=/tmp \
    --env TMP=/tmp \
    --env TEMP=/tmp \
    "$CONTAINER_PATH" \
    test --no-build -f "$TEST_UUID" -t "$TEST_UUID"

echo "Test $TEST_UUID completed"