#!/bin/bash
# Create a writable directory for test outputs
mkdir -p /tmp/mfc_tmp
cd /tmp/mfc_tmp

# Copy MFC to writable location
cp -r /opt/MFC/* .

# Run tests in writable location
./mfc.sh test -a --no-build