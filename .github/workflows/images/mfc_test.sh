#!/bin/bash

# Exit on any error
set -ex

# apptainer run --fakeroot MFC/mfc_cpu.sif test -a --no-build

CONTAINER_IMAGE="$1"
EXEC_COMMAND="$2"

if [[ ! -f "$CONTAINER_IMAGE" ]]; then
    echo "Error: Container image '$CONTAINER_IMAGE' not found!" >&2
    exit 1
fi

echo "Container: $CONTAINER_IMAGE"
echo "Command: $EXEC_COMMAND"

if apptainer exec --fakeroot --writable-tmpfs "$CONTAINER_IMAGE" /bin/bash -c "$EXEC_COMMAND"; then
    echo "Tests completed successfully!"
else
    echo "Tests failed with exit code $?" >&2
    exit 1
fi