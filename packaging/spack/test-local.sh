#!/bin/bash
# Local test script for Spack package
# Mimics the CI workflow

set -e

echo "==> Installing Spack..."
if [ ! -d "spack-test" ]; then
    git clone --depth=1 https://github.com/spack/spack.git spack-test
fi

cd spack-test

echo "==> Setting up Spack environment..."
. share/spack/setup-env.sh

echo "==> Finding compilers..."
spack compiler find

echo "==> Installing MFC package into Spack..."
mkdir -p var/spack/repos/builtin/packages/mfc
cp ../package.py var/spack/repos/builtin/packages/mfc/

echo ""
echo "==> Running Spack style check..."
spack style var/spack/repos/builtin/packages/mfc/package.py

echo ""
echo "==> Running Spack audit..."
spack audit packages mfc

echo ""
echo "==> Showing package info..."
spack info mfc

echo ""
echo "==> Testing default spec..."
spack spec mfc

echo ""
echo "==> Testing minimal spec..."
spack spec mfc~mpi~post_process

echo ""
echo "✅ All tests passed!"

