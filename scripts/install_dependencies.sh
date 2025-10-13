#!/bin/bash
# A solution for consolidating dependency installation commands for different operating systems (for macOS, Debian/Ubuntu, RedHat/CentOS)
# Allows minimal repetition in command installations in multiple CI workflow files.
# CI logs on which commands are being run.

# Instructions:
# 1. Run this command: chmod +x ./scripts/install_dependencies.sh
# 2. Run the installer: ./scripts/install_dependencies.sh

set -euo pipefail  
set -x              

OS="$(uname -s)"
echo "Installing dependencies for $OS OS"

# macOS
if [[ "$OSTYPE" == "darwin"* ]]; then
    # Use Homebrew for macOS package management
    echo "Detected macOS."
    export PATH="/opt/homebrew/bin:/opt/homebrew/sbin:$PATH"
    brew update

    # Packages to install
    declare -A packages
    packages=(
        ["cmake"]="cmake"
        ["gcc"]="gcc"
        ["python3"]="python@3.10"
        ["boost"]="boost"
        ["protobuf"]="protobuf"
    )

    # Install missing packages
    for exe in "${!packages[@]}"; do
        formula="${packages[$exe]}"
        if ! command -v "$exe" >/dev/null 2>&1; then
            echo "Installing $formula..."
            brew install --verbose "$formula"
        else
            echo "$exe already installed"
        fi
    done

    # Verification
    echo "Verifying installations..."
    which cmake; cmake --version
    which gcc; gcc --version
    which python3; python3 --version
    brew list boost || echo "boost not found"
    brew list protobuf || echo "protobuf not found"

    echo "macOS dependencies installed successfully."


# Debian/Ubuntu
elif [[ -f /etc/debian_version ]]; then
    echo "Detected Debian/Ubuntu"
    sudo apt-get update -y

    # Install all required packages (this time with openmpi and fftw)
    pkgs="tar wget make cmake gcc g++ python3 python3-dev openmpi-bin libopenmpi-dev fftw3 libfftw3-dev protobuf-compiler libboost-all-dev"
    sudo apt-get install -y $pkgs

    # Verification
    echo "Verifying installations..."
    cmake --version || echo "cmake not found"
    gcc --version || echo "gcc not found"
    python3 --version || echo "python3 not found"
    mpirun --version || echo "mpirun not found"
    ldconfig -p | grep -q libfftw3 || echo "FFTW library not found"

    echo "Linux dependencies installed successfully."


# RedHat/CentOS
elif [[ -f /etc/redhat-release ]]; then
    echo "Detected RedHat/CentOS"
    sudo yum install -y epel-release

    # Install all required packages (this time with openmpi and fftw)
    pkgs="tar wget make cmake3 gcc gcc-c++ python3 openmpi fftw protobuf boost-devel"
    sudo yum install -y $pkgs

    # Verification
    echo "Verifying installations..."
    cmake3 --version || echo "cmake3 not found"
    gcc --version || echo "gcc not found"
    python3 --version || echo "python3 not found"
    mpirun --version || echo "mpirun not found"
    rpm -q fftw || echo "FFTW library not found"

    echo "RedHat/CentOS dependencies installed successfully."


# If Unsupported OS
else
    echo "Unsupported OS: $OSTYPE"
    exit 1
fi

echo "Dependencies installed successfully."
