#!/bin/bash
# A solution for consolidating dependency installation commands for different operating systems (for macOS, Debian/Ubuntu, RedHat/CentOS)
# Allows minimal repetition in command installations in multiple CI workflow files.
# CI logs on which commands are being run.

# Instructions:
# 1. Run this command: chmod +x ./scripts/install_dependencies.sh
# 2. Run the installer: ./scripts/install_dependencies.sh

set -euo pipefail  
set -x
export DEBIAN_FRONTEND=noninteractive

OS="$(uname -s)"
echo "Installing dependencies for $OS OS"

if command -v sudo >/dev/null 2>&1; then
    SUDO="sudo"
else
    SUDO=""
fi

# macOS
# if [[ "$OSTYPE" == "darwin"* ]]; then
#     # Use Homebrew for macOS package management
#     echo "Detected macOS."
#     export PATH="/opt/homebrew/bin:/opt/homebrew/sbin:$PATH"
#     brew update

#     # Packages to install
#     declare -A packages
#     packages=(
#         ["cmake"]="cmake"
#         ["gcc"]="gcc"
#         ["python3"]="python@3.10"
#         ["boost"]="boost"
#         ["protobuf"]="protobuf"
#     )

#     # Install missing packages
#     for exe in "${!packages[@]}"; do
#         formula="${packages[$exe]}"
#         if ! command -v "$exe" >/dev/null 2>&1; then
#             echo "Installing $formula..."
#             brew install --verbose "$formula"
#         else
#             echo "$exe already installed"
#         fi
#     done

#     # Verification
#     echo "Verifying installations..."
#     which cmake; cmake --version
#     which gcc; gcc --version
#     which python3; python3 --version
#     brew list boost || echo "boost not found"
#     brew list protobuf || echo "protobuf not found"

#     echo "macOS dependencies installed successfully."
# macOS
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Detected macOS."
    export PATH="/opt/homebrew/bin:/opt/homebrew/sbin:$PATH"
    brew update
    brew upgrade || true

    # Install required packages (mirrors workflow setup)
    pkgs=(coreutils python fftw hdf5 gcc@15 boost open-mpi lapack cmake protobuf)

    for pkg in "${pkgs[@]}"; do
        if ! brew list "$pkg" >/dev/null 2>&1; then
            echo "Installing $pkg..."
            brew install "$pkg" || brew reinstall "$pkg"
        else
            echo "$pkg already installed."
        fi
    done

    # Fix potential linking issues for Homebrew-installed libraries
    brew link --overwrite python || true
    brew link --overwrite boost || true
    brew link --overwrite hdf5 || true
    brew link --overwrite cmake || true
    brew link --overwrite protobuf || true

    # Environment setup
    echo "FC=gfortran-15" >> "$GITHUB_ENV"
    echo "BOOST_INCLUDE=/opt/homebrew/include/" >> "$GITHUB_ENV"

    echo "macOS dependencies installed successfully."

# Debian/Ubuntu
elif [[ -f /etc/debian_version ]]; then
    echo "Detected Debian/Ubuntu"

    # Install all required packages (this time with openmpi and fftw)
    if [[ "${CI_INTEL:-false}" == "true" ]]; then
        echo "Intel OneAPI environment detected. Installing Intel compilers and MPI..."
        wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
        $SUDO apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
        $SUDO add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
        $SUDO apt-get update
        $SUDO apt-get install -y intel-oneapi-compiler-fortran intel-oneapi-mpi intel-oneapi-mpi-devel
        source /opt/intel/oneapi/setvars.sh
        if [[ -n "${GITHUB_ENV:-}" ]]; then
            echo "source /opt/intel/oneapi/setvars.sh" >> "$GITHUB_ENV"
        fi
        printenv >> "$GITHUB_ENV"
        echo "Intel OneAPI installation completed."
    else
        echo "Standard GNU environment detected. Installing open-source toolchain..."
        $SUDO apt-get update -y
        pkgs="tar wget make cmake gcc g++ python3 python3-dev openmpi-bin libopenmpi-dev fftw3 libfftw3-dev protobuf-compiler libboost-all-dev"
        $SUDO apt-get install -y $pkgs

        # Verification
        echo "Verifying installations..."
        cmake --version || echo "cmake not found"
        gcc --version || echo "gcc not found"
        python3 --version || echo "python3 not found"
        mpirun --version || echo "mpirun not found"
        ldconfig -p | grep -q libfftw3 || echo "FFTW library not found"

        echo "Linux (GNU) dependencies installed successfully."
    fi

# RedHat/CentOS
elif [[ -f /etc/redhat-release ]]; then
    echo "Detected RedHat/CentOS"
    $SUDO yum install -y epel-release

    # Install all required packages (this time with openmpi and fftw)
    pkgs="tar wget make cmake3 gcc gcc-c++ python3 openmpi fftw protobuf boost-devel"
    $SUDO yum install -y $pkgs

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
