#!/bin/bash
# A solution for consolidating dependency installation commands for different operating systems (for macOS, Debian/Ubuntu, RedHat/CentOS)
# Allows minimal repetition in command installations in multiple CI workflow files.
# CI logs on which commands are being run.

# Instructions:
# 1. Run this command: chmod +x ./scripts/install_dependencies.sh
# 2. Run the installer: ./scripts/install_dependencies.sh

set -euo pipefail
set -x

SUDO=${SUDO:-sudo}
CI_INTEL=${CI_INTEL:-false}

log() { echo ">>> $*"; }

# Guarded append to GITHUB_ENV (safe when running locally)
append_to_github_env() {
    local line="$1"
    if [[ -n "${GITHUB_ENV:-}" ]]; then
        printf '%s\n' "$line" >> "$GITHUB_ENV"
    else
        printf '%s\n' "$line"
    fi
}

verify_tool() {
    local t="$1"
    if command -v "$t" >/dev/null 2>&1; then
        "$t" --version 2>&1 | head -n1 || true
    else
        echo "$t not found"
    fi
}

# macOS / Homebrew
install_macos() {
    log "Detected macOS"

    if ! command -v brew >/dev/null 2>&1; then
        log "Homebrew not found. Please install Homebrew first."
        return 1
    fi

    export PATH="/opt/homebrew/bin:/opt/homebrew/sbin:$PATH"
    brew update
    brew upgrade || true

    pkgs=(coreutils python fftw hdf5 gcc@15 boost open-mpi lapack cmake protobuf)

    for pkg in "${pkgs[@]}"; do
        # Use brew list to detect installed formulae (correct for lib-only formulae)
        if ! brew list --formula | grep -q "^${pkg%%@*}$"; then
            log "Installing ${pkg}..."
            brew install "${pkg}" || brew reinstall "${pkg}" || log "Failed to install ${pkg} (continuing)"
        else
            log "${pkg} already installed"
        fi
    done

    # Attempt to fix common Homebrew linking quirks
    brew link --overwrite python || true
    brew link --overwrite boost || true
    brew link --overwrite hdf5 || true
    brew link --overwrite cmake || true
    brew link --overwrite protobuf || true

    # Export useful environment vars in GH Actions only (guarded)
    append_to_github_env "FC=gfortran-15"
    append_to_github_env "BOOST_INCLUDE=/opt/homebrew/include/"

    log "Verifying installations..."
    # Use guarded verification so script doesn't exit on set -e for these checks
    if ! verify_tool cmake >/dev/null 2>&1; then echo "cmake not found"; fi
    if ! verify_tool gcc >/dev/null 2>&1; then echo "gcc not found"; fi
    if ! verify_tool gfortran >/dev/null 2>&1; then echo "gfortran not found"; fi
    if ! verify_tool python3 >/dev/null 2>&1; then echo "python3 not found"; fi
    if ! verify_tool mpirun >/dev/null 2>&1; then echo "mpirun not found"; fi

    if brew list --formula | grep -q "^fftw$"; then
        echo "FFTW library found"
    else
        echo "FFTW library not found"
    fi

    log "macOS dependencies installed successfully."
}

# Debian / Ubuntu
install_debian() {
    log "Detected Debian/Ubuntu"

    $SUDO apt-get update -y
    $SUDO apt-get install -y build-essential wget tar cmake python3 python3-dev pkg-config

    if [[ "${CI_INTEL}" == "true" ]]; then
        log "Intel OneAPI: adding repo with signed-by keyring (secure method)"
        # download key and install to /etc/apt/keyrings (modern secure approach)
        wget -q https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB -O /tmp/intel-pubkey.PUB
        $SUDO mkdir -p /etc/apt/keyrings
        $SUDO install -m 0644 /tmp/intel-pubkey.PUB /etc/apt/keyrings/intel-archive-keyring.gpg
        echo "deb [signed-by=/etc/apt/keyrings/intel-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | $SUDO tee /etc/apt/sources.list.d/oneapi.list > /dev/null
        $SUDO apt-get update -y
        $SUDO apt-get install -y intel-oneapi-compiler-fortran intel-oneapi-mpi intel-oneapi-mpi-devel || log "Intel packages failed; continuing"

        # capture environment changes made by Intel setvars.sh and append only deltas to GITHUB_ENV
        if [[ -f /opt/intel/oneapi/setvars.sh ]]; then
            env > /tmp/env_before
            source /opt/intel/oneapi/setvars.sh || log "sourcing Intel setvars.sh failed"
            if [[ -n "${GITHUB_ENV:-}" ]]; then
                # append only newly added / changed env vars
                diff --unchanged-line-format="" --old-line-format="" --new-line-format="%L" <(sort /tmp/env_before) <(env | sort) >> "$GITHUB_ENV" || true
            fi
            rm -f /tmp/env_before
        fi
    fi

    # Install common dev packages including fftw/openmpi/protobuf/boost
    $SUDO apt-get install -y libfftw3-dev libopenmpi-dev openmpi-bin libprotobuf-dev protobuf-compiler libboost-all-dev

    log "Verifying installations..."
    if ! verify_tool cmake >/dev/null 2>&1; then echo "cmake not found"; fi
    if ! verify_tool gfortran >/dev/null 2>&1; then echo "gfortran not found"; fi
    if ! verify_tool python3 >/dev/null 2>&1; then echo "python3 not found"; fi
    if ! verify_tool mpirun >/dev/null 2>&1; then echo "mpirun not found"; fi
    if ! ldconfig -p | grep -q libfftw3; then echo "FFTW library not found"; fi
}

# RHEL / CentOS
install_rhel() {
    log "Detected RHEL/CentOS"

    pkgs="tar wget make cmake3 gcc gcc-c++ python3-devel openmpi openmpi-devel fftw fftw-devel protobuf protobuf-devel boost-devel"
    $SUDO yum install -y $pkgs || $SUDO dnf install -y $pkgs || log "yum/dnf install failed"

    log "Verifying installations..."
    if command -v cmake3 >/dev/null 2>&1; then
        verify_tool cmake3
    elif command -v cmake >/dev/null 2>&1; then
        verify_tool cmake
    else
        echo "cmake not found"
    fi

    if ! verify_tool gfortran >/dev/null 2>&1; then echo "gfortran not found"; fi
    if ! verify_tool python3 >/dev/null 2>&1; then echo "python3 not found"; fi
    if ! verify_tool mpirun >/dev/null 2>&1; then echo "mpirun not found"; fi
    if ! ldconfig -p | grep -q libfftw3; then echo "FFTW library not found"; fi
}

# Boost compile check
boost_check() {
    log "Verifying Boost by compiling a tiny program..."
    cat <<'EOF' > /tmp/boost_check.cpp
#include <boost/version.hpp>
#include <iostream>
int main(){ std::cout << "Boost version: " << BOOST_LIB_VERSION << std::endl; return 0;}
EOF

    # Try pkg-config first (if available)
    if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists boost; then
        g++ /tmp/boost_check.cpp -o /tmp/boost_check $(pkg-config --cflags --libs boost) || true
    else
        # Try a plain compile (system include path), fallback to /opt/homebrew/include
        if g++ /tmp/boost_check.cpp -o /tmp/boost_check 2>/dev/null; then
            true
        else
            if [[ -d /opt/homebrew/include ]]; then
                g++ -I/opt/homebrew/include /tmp/boost_check.cpp -o /tmp/boost_check || true
            else
                log "Could not compile boost_check; boost include path unknown"
            fi
        fi
    fi

    if [[ -x /tmp/boost_check ]]; then
        /tmp/boost_check || true
    else
        log "boost_check binary not available"
    fi
}

# Main entry
main() {
    case "$(uname -s)" in
        Darwin)
            install_macos
            ;;
        *)
            if [[ -f /etc/debian_version ]]; then
                install_debian
            elif [[ -f /etc/redhat-release ]] || [[ -f /etc/centos-release ]]; then
                install_rhel
            else
                log "Unsupported OS: $OSTYPE"
                return 1
            fi
            ;;
    esac

    boost_check
    log "Done."
}

main "$@"