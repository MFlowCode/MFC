#!/bin/bash

# If we downloaded our own version of CMake, let us make it accessible to
# any command we run hereafter.

if [ -f "$(pwd)/build/cmake/bin/cmake" ]; then
    export PATH="$(pwd)/build/cmake/bin:$PATH"

    log "Found$MAGENTA CMake$COLOR_RESET: $MAGENTA$(pwd)/build/cmake/bin/cmake$COLOR_RESET."
fi

bShouldInstallCMake=false
if ! command -v cmake > /dev/null 2>&1; then
    # Not installed
    bShouldInstallCMake=true

    warn "$MAGENTA""CMake$COLOR_RESET is not installed."
else
    cmake_verstr=$(cmake --version | tr ' ' '\n' | sed -n 3p)
    cmake_major=$(echo $cmake_verstr | tr '.' '\n' | sed -n 1p)
    cmake_minor=$(echo $cmake_verstr | tr '.' '\n' | sed -n 2p)
    cmake_patch=$(echo $cmake_verstr | tr '.' '\n' | sed -n 3p)
    cmake_version="$(printf %05d%05d%05d $cmake_major $cmake_minor $cmake_patch)"

    MFC_CMAKE_MIN_VERSTR=$(cat CMakeLists.txt | grep cmake_minimum_required | head -n 1 | sed 's/[^0-9,.]*//g')
    MFC_CMAKE_MIN_MAJOR=$(echo $MFC_CMAKE_MIN_VERSTR | tr '.' '\n' | head -n 1)
    MFC_CMAKE_MIN_MINOR=$(echo $MFC_CMAKE_MIN_VERSTR | tr '.' '\n' | head -n 2 | tail -n 1)
    MFC_CMAKE_MIN_PATCH=0
    MFC_CMAKE_MIN_VERSION="$(printf %05d%05d%05d $MFC_CMAKE_MIN_MAJOR $MFC_CMAKE_MIN_MINOR $MFC_CMAKE_MIN_PATCH)"

    if [ "$cmake_version" -lt "$MFC_CMAKE_MIN_VERSION" ]; then
        # Out of date
        bShouldInstallCMake=true

        log "$MAGENTA CMake$YELLOW is out of date (current: $MAGENTA$cmake_verstr$YELLOW < minimum: $MAGENTA$MFC_CMAKE_MIN_VERSTR$YELLOW).$COLOR_RESET"
    fi
fi

if [ "$bShouldInstallCMake" = true ]; then
    version="3.25.2"
    arch="$(uname -m)"

    bErrorDoSelfDownload=false
    if [[ "$OSTYPE" == "darwin"* ]]; then
        error "Cannot fetch$MAGENTA CMake$COLOR_RESET for$MAGENTA macOS$COLOR_RESET."
        bErrorDoSelfDownload=true
    fi

    if ! [[ "$arch" == "x86_64" || "$arch" == "aarch64" ]]; then
        error "Cannot fetch$MAGENTA CMake$COLOR_RESET for architecture $MAGENTA$arch$COLOR_RESET."
        bErrorDoSelfDownload=true
    fi

    if [ "$bErrorDoSelfDownload" = true ]; then
        log "Please install$MAGENTA CMake$COLOR_RESET manually:"
        log " - via$MAGENTA Aptitude$COLOR_RESET:  $ sudo apt install cmake"
        log " - via$MAGENTA Homebrew$COLOR_RESET:  $ brew install cmake"
        log " - via$MAGENTA Spack$COLOR_RESET:     $ spack install cmake"
        log " - via$MAGENTA Pacman$COLOR_RESET:    $ sudo pacman -S cmake"
        log " - via$MAGENTA Lmod$COLOR_RESET:      $ module load cmake (if available)"
        log " - via$MAGENTA cmake.org$COLOR_RESET: https://cmake.org/download/"

        exit 1
    fi

    if [ -d "$(pwd)/build/cmake" ]; then
        rm -rf "$(pwd)/build/cmake"
    fi

    mkdir -p "$(pwd)/build/cmake"

    filename="cmake-$version-linux-$arch.sh"
    repository="https://github.com/Kitware/CMake"
    url="$repository/releases/download/v$version/$filename"

    log "Downloading$MAGENTA CMake v$version$COLOR_RESET for $MAGENTA$arch$COLOR_RESET from $CYAN$repository$COLOR_RESET."

    cmake_fatal_error() {
        log   "$MAGENTA""CMake$COLOR_RESET is not discoverable or is an older release, incompatible with$MAGENTA MFC$COLOR_RESET. Please download"
        log   "or install a recent version of$MAGENTA CMake$COLOR_RESET to get past this step. If you are currently on a"
        log   "managed system like a cluster, provided there is no suitable environment module, you can"
        log   "either build it from source, or get it via$MAGENTA Spack$COLOR_RESET."
        log   "- The minimum required version is currently$MAGENTA CMake$COLOR_RESET v$MFC_CMAKE_MIN_MAJOR.$MFC_CMAKE_MIN_MINOR.$MFC_CMAKE_MIN_PATCH."
        log   "- We attempted to download$MAGENTA CMake$COLOR_RESET v$version from $url."

        exit 1
    }

    if ! command -v wget > /dev/null 2>&1; then
        error "$MAGENTA""wget$COLOR_RESET is not installed but is necessary to download$MAGENTA CMake$COLOR_RESET."

        cmake_fatal_error
    fi

    if ! wget -P "$(pwd)/build/cmake" "$url"; then
        error "Failed to download a compatible version of$MAGENTA CMake$COLOR_RESET."

        cmake_fatal_error
    fi

    log "Installing$MAGENTA CMake$COLOR_RESET into $MAGENTA$(pwd)/build/cmake$COLOR_RESET."

    if ! $SHELL "$(pwd)/build/cmake/$filename" "--skip-license" "--prefix=$(pwd)/build/cmake"; then
        error "Failed to install a compatible version of CMake."

        cmake_fatal_error
    fi

    rm "$(pwd)/build/cmake/$filename"

    ok "Found$MAGENTA CMake$COLOR_RESET: $MAGENTA$(pwd)/build/cmake/bin/cmake$COLOR_RESET."
    export PATH="$(pwd)/build/cmake/bin:$PATH"
fi
