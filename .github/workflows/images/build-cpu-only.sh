#!/bin/bash
# Simple script to build just the CPU image for testing

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}Building MFC CPU container image...${NC}"

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

# Clean up any existing images
echo -e "${YELLOW}Cleaning up existing images...${NC}"
rm -f mfc_cpu.sif

# Build CPU image with enhanced space management
echo -e "${GREEN}Building CPU image (this may take 15-30 minutes)...${NC}"
echo "Progress will be shown below:"

# Check available space before building
AVAILABLE_SPACE=$(df /tmp | awk 'NR==2 {print $4}')
AVAILABLE_GB=$((AVAILABLE_SPACE / 1024 / 1024))
echo "Available space: ${AVAILABLE_GB}GB"

if [ $AVAILABLE_GB -lt 5 ]; then
    echo -e "${RED}Error: Insufficient space for building. Need at least 5GB free.${NC}"
    echo "Current available: ${AVAILABLE_GB}GB"
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
    echo -e "${GREEN}✓ CPU image built successfully!${NC}"
    
    echo ""
    echo "Image details:"
    ls -lh mfc_cpu.sif
    
    echo ""
    echo "Testing the built image..."
    if apptainer exec mfc_cpu.sif /opt/MFC/mfc.sh --help > /dev/null 2>&1; then
        echo -e "${GREEN}✓ Image test passed - MFC is accessible${NC}"
    else
        echo -e "${YELLOW}⚠ Warning: Image built but MFC may not be properly installed${NC}"
    fi
    
    # Clean up build caches after successful build
    rm -rf "$APPTAINER_CACHEDIR" "$SINGULARITY_CACHEDIR" 2>/dev/null || true
    
    echo ""
    echo "To run the image:"
    echo "  apptainer run mfc_cpu.sif run examples/2D_shockbubble/case.py -n 4"
    echo "  apptainer run mfc_cpu.sif test -j 4"
    echo ""
    echo "For help:"
    echo "  apptainer run-help mfc_cpu.sif"
    
else
    echo -e "${RED}✗ CPU image build failed${NC}"
    cleanup_build
    echo ""
    echo "Common solutions:"
    echo "1. Make sure you have enough disk space (need ~5GB free, have ${AVAILABLE_GB}GB)"
    echo "2. Check your internet connection"
    echo "3. Try running with sudo privileges if available"
    echo "4. Check if /tmp has sufficient space"
    exit 1
fi