#!/bin/bash
# Script to build and test Apptainer/Singularity images locally

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}Starting local Apptainer/Singularity image builds...${NC}"

# Initialize counter for successful builds
IMAGE_COUNT=0

# Check if apptainer is installed
if ! command -v apptainer &> /dev/null; then
    echo -e "${RED}Error: Apptainer is not installed.${NC}"
    echo "Please install Apptainer first: https://apptainer.org/docs/admin/main/installation.html"
    exit 1
fi

# Set up cache directories with space optimization
echo -e "${YELLOW}Setting up cache directories with space optimization...${NC}"
export APPTAINER_CACHEDIR="/tmp/apptainer-build-cache"
export SINGULARITY_CACHEDIR="/tmp/singularity-build-cache"
mkdir -p "$APPTAINER_CACHEDIR" "$SINGULARITY_CACHEDIR"

# Clean up any existing cache
rm -rf "$APPTAINER_CACHEDIR"/* 2>/dev/null || true
rm -rf "$SINGULARITY_CACHEDIR"/* 2>/dev/null || true

echo "Cache directories: $APPTAINER_CACHEDIR, $SINGULARITY_CACHEDIR"

# Enable fakeroot if not already enabled
echo -e "${YELLOW}Checking fakeroot configuration...${NC}"
if ! apptainer config fakeroot --list 2>/dev/null | grep -q $(whoami); then
    echo "Fakeroot not enabled for current user. Enabling..."
    if ! sudo apptainer config fakeroot --enable $(whoami); then
        echo -e "${RED}Warning: Could not enable fakeroot. Builds may fail.${NC}"
        echo "You may need to run manually:"
        echo "sudo apptainer config fakeroot --enable $(whoami)"
    else
        echo -e "${GREEN}✓ Fakeroot enabled for $(whoami)${NC}"
    fi
else
    echo -e "${GREEN}✓ Fakeroot already enabled for $(whoami)${NC}"
fi

# Build CPU image with space management
echo -e "${GREEN}Building CPU image...${NC}"
echo "This may take 10-20 minutes depending on your system..."

# Check available space before building
AVAILABLE_SPACE=$(df /tmp | awk 'NR==2 {print $4}')
AVAILABLE_GB=$((AVAILABLE_SPACE / 1024 / 1024))
echo "Available space: ${AVAILABLE_GB}GB"

if [ $AVAILABLE_GB -lt 5 ]; then
    echo -e "${RED}Error: Insufficient space for building. Need at least 5GB free.${NC}"
    exit 1
fi

# Build with automatic cleanup on failure
cleanup_build() {
    echo "Cleaning up build artifacts..."
    rm -rf "$APPTAINER_CACHEDIR"/* 2>/dev/null || true
    rm -rf "$SINGULARITY_CACHEDIR"/* 2>/dev/null || true
    rm -f mfc_cpu.sif 2>/dev/null || true
}

if apptainer build --tmpdir /tmp mfc_cpu.sif Singularity.cpu; then
    echo -e "${GREEN}✓ CPU image built successfully${NC}"
    IMAGE_COUNT=$((IMAGE_COUNT + 1))
else
    echo -e "${RED}✗ CPU image build failed${NC}"
    cleanup_build
    echo "Trying with --fakeroot flag..."
    if apptainer build --fakeroot --tmpdir /tmp mfc_cpu.sif Singularity.cpu; then
        echo -e "${GREEN}✓ CPU image built successfully with fakeroot${NC}"
        IMAGE_COUNT=$((IMAGE_COUNT + 1))
    else
        echo -e "${RED}✗ CPU image build failed even with fakeroot${NC}"
        cleanup_build
    fi
fi

# Build CPU benchmark image
echo -e "${GREEN}Building CPU benchmark image...${NC}"
if apptainer build mfc_cpu_bench.sif Singularity.cpu_bench; then
    echo -e "${GREEN}✓ CPU benchmark image built successfully${NC}"
    IMAGE_COUNT=$((IMAGE_COUNT + 1))
else
    echo -e "${RED}✗ CPU benchmark image build failed${NC}"
    echo "Trying with --fakeroot flag..."
    if apptainer build --fakeroot mfc_cpu_bench.sif Singularity.cpu_bench; then
        echo -e "${GREEN}✓ CPU benchmark image built successfully with fakeroot${NC}"
        IMAGE_COUNT=$((IMAGE_COUNT + 1))
    else
        echo -e "${RED}✗ CPU benchmark image build failed even with fakeroot${NC}"
    fi
fi

# Build GPU image
echo -e "${GREEN}Building GPU image...${NC}"
if apptainer build mfc_gpu.sif Singularity.gpu; then
    echo -e "${GREEN}✓ GPU image built successfully${NC}"
    IMAGE_COUNT=$((IMAGE_COUNT + 1))
else
    echo -e "${RED}✗ GPU image build failed${NC}"
    echo "Trying with --fakeroot flag..."
    if apptainer build --fakeroot mfc_gpu.sif Singularity.gpu; then
        echo -e "${GREEN}✓ GPU image built successfully with fakeroot${NC}"
        IMAGE_COUNT=$((IMAGE_COUNT + 1))
    else
        echo -e "${RED}✗ GPU image build failed even with fakeroot${NC}"
    fi
fi

# Build GPU benchmark image
echo -e "${GREEN}Building GPU benchmark image...${NC}"
if apptainer build mfc_gpu_bench.sif Singularity.gpu_bench; then
    echo -e "${GREEN}✓ GPU benchmark image built successfully${NC}"
    IMAGE_COUNT=$((IMAGE_COUNT + 1))
else
    echo -e "${RED}✗ GPU benchmark image build failed${NC}"
    echo "Trying with --fakeroot flag..."
    if apptainer build --fakeroot mfc_gpu_bench.sif Singularity.gpu_bench; then
        echo -e "${GREEN}✓ GPU benchmark image built successfully with fakeroot${NC}"
        IMAGE_COUNT=$((IMAGE_COUNT + 1))
    else
        echo -e "${RED}✗ GPU benchmark image build failed even with fakeroot${NC}"
    fi
fi

echo -e "${GREEN}Build process completed!${NC}"
echo ""
echo "Successfully built $IMAGE_COUNT out of 4 images."
echo ""
echo "Built images:"
ls -lh *.sif 2>/dev/null || echo "No .sif files found"

if [ $IMAGE_COUNT -gt 0 ]; then
    echo ""
    echo "To test an image, run:"
    echo "  apptainer run mfc_cpu.sif test -j 4"
    echo "  apptainer run --nv mfc_gpu.sif test -j 4 --gpu"
    echo ""
    echo "For help with any image:"
    echo "  apptainer run-help mfc_cpu.sif"
else
    echo -e "${RED}No images were built successfully. Check the errors above.${NC}"
    exit 1
fi