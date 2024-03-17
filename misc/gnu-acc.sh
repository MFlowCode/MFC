#!/bin/bash

# Author:      Henry Le Berre
# Description: Builds gcc/g++/gfortran with nvptx-none and amdgcn-amdhsa offloading
#              from the devel/omp/gcc-12 branch.

set -e

JOBS=1
INSTALL_PREFIX=""

while [[ $# -gt 0 ]]; do
    case $1 in
    --prefix)        INSTALL_PREFIX="$2"; shift; shift; ;;
    -a|--amd)        AMD=1;               shift;        ;;
    -n|--nvidia)     NVIDIA=1;            shift;        ;;
    -d|--temp-dir)   TEMP_DIR="$2";       shift; shift; ;;
    -j|--jobs)       JOBS="$2";           shift; shift; ;;
    -h|--help)       PRINT_HELP=1;        shift;        ;;
    -*|--*)
        echo "Unknown option $1"
        exit 1
        ;;
    esac
done

if [[ -n "$PRINT_HELP" ]]; then
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --prefix <path>  Installation prefix."
    echo "  -a, --amd        Offload for AMD GPUs."
    echo "  -n, --nvidia     Offload for NVIDIA GPUs."
    echo "  -d, --temp-dir   Temporary directory to use for building. [default: mktemp -d]"
    echo "  -j, --jobs <n>   Number of jobs to run in parallel [default: $JOBS]."
    echo "  -h, --help       Print this help message."
    exit 0
fi

if [[ -z "$INSTALL_PREFIX" ]]; then
    echo "Error: No installation prefix specified."
    exit 1
fi

if [[ -z "$TEMP_DIR" ]]; then
    TEMP_DIR="$(mktemp -d)"
fi

if [[ ! -z "$NVIDIA" ]]; then
    git clone --depth=1                                           \
              'https://github.com/MentorEmbedded/nvptx-tools.git' \
              "$TEMP_DIR/repositories/nvptx-tools"

    mkdir -p "$TEMP_DIR/staging/nvptx-tools"

    (cd "$TEMP_DIR/staging/nvptx-tools"; "$TEMP_DIR/repositories/nvptx-tools/configure" --prefix="$INSTALL_PREFIX")

    make install \
        -j $JOBS \
        -C "$TEMP_DIR/staging/nvptx-tools"

    rm -rf "$TEMP_DIR/staging/nvptx-tools"
fi

if [[ ! -z "$AMD" ]]; then
    mkdir -p "$TEMP_DIR/repositories/llvm-project"

    wget -O "$TEMP_DIR/repositories/llvmorg.tar.gz" \
            'https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-13.0.1.tar.gz'

    tar -C   "$TEMP_DIR/repositories/llvm-project"   \
        -vxf "$TEMP_DIR/repositories/llvmorg.tar.gz" \
        --strip-components 1

    rm "$TEMP_DIR/repositories/llvmorg.tar.gz"

    cmake                                     \
        -D LLVM_ENABLE_PROJECTS=lld           \
        -B "$TEMP_DIR/staging/llvm-project"   \
        -D LLVM_TARGETS_TO_BUILD="X86;AMDGPU" \
        -S "$TEMP_DIR/repositories/llvm-project/llvm"

    make -C "$TEMP_DIR/staging/llvm-project" \
         -j $JOBS

    mkdir -p "$INSTALL_PREFIX/amdgcn-amdhsa/bin"

    cp -a "$TEMP_DIR/staging/llvm-project/bin/llvm-ar" \
            "$INSTALL_PREFIX/amdgcn-amdhsa/bin/ar"
    cp -a "$TEMP_DIR/staging/llvm-project/bin/llvm-ar" \
            "$INSTALL_PREFIX/amdgcn-amdhsa/bin/ranlib"
    cp -a "$TEMP_DIR/staging/llvm-project/bin/llvm-mc" \
            "$INSTALL_PREFIX/amdgcn-amdhsa/bin/as"
    cp -a "$TEMP_DIR/staging/llvm-project/bin/llvm-nm" \
            "$INSTALL_PREFIX/amdgcn-amdhsa/bin/nm"
    cp -a "$TEMP_DIR/staging/llvm-project/bin/lld"     \
            "$INSTALL_PREFIX/amdgcn-amdhsa/bin/ld"

    rm -rf "$TEMP_DIR/repositories/llvm-project"
fi

git clone --depth=1                                    \
          --branch master                              \
          'git://sourceware.org/git/newlib-cygwin.git' \
          "$TEMP_DIR/repositories/newlib-cygwin"

git clone --depth=1                       \
          --branch devel/omp/gcc-12       \
          'git://gcc.gnu.org/git/gcc.git' \
          "$TEMP_DIR/repositories/gcc"

(cd "$TEMP_DIR/repositories/gcc"; "$TEMP_DIR/repositories/gcc/contrib/download_prerequisites")

ln -s "$TEMP_DIR/repositories/newlib-cygwin/newlib" \
      "$TEMP_DIR/repositories/gcc/newlib"

if [[ ! -z "$AMD" ]]; then
    mkdir -p "$TEMP_DIR/staging/amdgcn-amdhsa-gcc"

    (cd "$TEMP_DIR/staging/amdgcn-amdhsa-gcc"; "$TEMP_DIR/repositories/gcc/configure" \
        --prefix="$INSTALL_PREFIX"                                  \
        --target=amdgcn-amdhsa                                      \
        --with-newlib                                               \
        --enable-languages=c,c++,fortran,lto                        \
        --disable-bootstrap                                         \
        --disable-libquadmath                                       \
        --disable-sjlj-exceptions                                   \
        --with-build-time-tools="$INSTALL_PREFIX/amdgcn-amdhsa/bin" \
        --enable-as-accelerator-for=x86_64-pc-linux-gnu)

    make -j $JOBS \
         -C "$TEMP_DIR/staging/amdgcn-amdhsa-gcc"

    make install  \
         -j $JOBS \
         -C "$TEMP_DIR/staging/amdgcn-amdhsa-gcc"

    rm -rf "$TEMP_DIR/staging/amdgcn-amdhsa-gcc"
fi

if [[ ! -z "$NVIDIA" ]]; then
    mkdir -p "$TEMP_DIR/staging/nvptx-none-gcc"

    (cd "$TEMP_DIR/staging/nvptx-none-gcc"; "$TEMP_DIR/repositories/gcc/configure" \
        --prefix="$INSTALL_PREFIX"                               \
        --target=nvptx-none                                      \
        --with-newlib                                            \
        --enable-languages=c,c++,fortran,lto                     \
        --disable-bootstrap                                      \
        --disable-libquadmath                                    \
        --disable-sjlj-exceptions                                \
        --enable-newlib-io-long-long                             \
        --with-build-time-tools="$INSTALL_PREFIX/nvptx-none/bin" \
        --enable-as-accelerator-for=x86_64-pc-linux-gnu)

    make -j $JOBS \
         -C "$TEMP_DIR/staging/nvptx-none-gcc"

    make install  \
         -j $JOBS \
         -C "$TEMP_DIR/staging/nvptx-none-gcc"
fi

mkdir -p "$TEMP_DIR/staging/host-gcc"

if [[ ! -z "$AMD" ]]; then
    OFFLOAD_OPT="--enable-offload-targets=amdgcn-amdhsa=$INSTALL_PREFIX/amdgcn-amdhsa"
fi

if [[ ! -z "$NVIDIA" ]]; then
    OFFLOAD_OPT="--enable-offload-targets=nvptx-none=$INSTALL_PREFIX/nvptx-none"
fi

(cd "$TEMP_DIR/staging/host-gcc"; "$TEMP_DIR/repositories/gcc/configure" \
    --host=`uname -m`-pc-linux-gnu         \
    --build=`uname -m`-pc-linux-gnu        \
    --target=`uname -m`-pc-linux-gnu       \
    --prefix="$INSTALL_PREFIX"             \
    --disable-multilib                     \
    --enable-languages=c,c++,fortran,lto   \
    --disable-bootstrap                    \
    "$OFFLOAD_OPT")

make -j $JOBS \
     -C "$TEMP_DIR/staging/host-gcc"

make install  \
     -j $JOBS \
     -C "$TEMP_DIR/staging/host-gcc"

rm -rf "$TEMP_DIR"
