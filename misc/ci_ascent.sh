#!/usr/bin/env bash

# This script is executed by the CI Workflow for MFC
# running Ascent @ ORNL through GitLab. (see .gitlab-ci.yml)

# (#) Exit if a command fails & Show each command
set -e
set -o pipefail
set -o xtrace

# (#) Useful constants
# ->  Name of repository (MFC/MFC-develop/etc..)
CI_MFC_DIR_NAME="$1"
# ->  Location of the cloned MFC by GitLab
CI_MFC_INIT_PATH="$(pwd)/$CI_MFC_DIR_NAME"
# ->  New/Desired location of MFC in our project directory to use jsrun
CI_MFC_DEST_PATH="$MEMBERWORK/cfd154/.ci/$CI_MFC_DIR_NAME"

# (#) Make sure we are in a cloned MFC directory
if [[ ! "$CI_MFC_INIT_PATH" == *"MFC"* ]]; then
    exit 1
fi

# (#) Copy files from initial MFC location to new
# ->  Remove new install folder if it exists (from a previous run)
if [[ -d "$CI_MFC_DEST_PATH" ]]; then
    rm -rf "$CI_MFC_DEST_PATH"
fi

# ->  Create the new MFC location
mkdir -p "$CI_MFC_DEST_PATH"
# ->  Copy files over
rsync -a "$CI_MFC_INIT_PATH/" "$CI_MFC_DEST_PATH"

# (#) Create symbolic link from the original MFC clone
#     to the new location of MFC so CI can find artifacts
rm -rf "./$CI_MFC_DIR_NAME"
ln -s "$CI_MFC_DEST_PATH" .

# (#) Now run & test MFC (in the new MFC directory)
cd "$CI_MFC_DEST_PATH"
echo "[CI] Running in $(pwd):"
# ->  Load GPU modules for Ascent (used for CPU/GPU since NVTX-support is required)
# ->  Disable command tracing for sourced file
set +o xtrace
. ./misc/load.sh -c a -m g
set -o xtrace
# (##) CPU Mode
./mfc.sh build -c release-cpu -j 8
./mfc.sh test  -c release-cpu
# (##) GPU Mode
./mfc.sh build -c release-gpu -j 8
./mfc.sh test  -c release-gpu
