#!/bin/bash

# Authors:     Alex Woods (GaTech) and Henry Le Berre (GaTech).
# Description: Builds GNU 13 compilers with amdgcn-amdhsa offloading.

set -v

work_dir=$HOME/USERSCRATCH/build-gcc-amdgpu/
install_dir=$HOME/tools/gcc/13

mkdir -p $work_dir
cd $work_dir

mkdir builds
build_dir=$work_dir/builds

## Clone and build llvm 13
git clone -b release/13.x https://github.com/llvm/llvm-project.git
 
llvmsrcdir=$work_dir/llvm-project
cd $build_dir
mkdir llvm
cd llvm
cmake -D 'LLVM_TARGETS_TO_BUILD=X86;AMDGPU' -D LLVM_ENABLE_PROJECTS=lld $llvmsrcdir/llvm
make -j $(nproc)

mkdir -p $install_dir/amdgcn-amdhsa/bin
cp -a bin/llvm-ar $install_dir/amdgcn-amdhsa/bin/ar
cp -a bin/llvm-ar $install_dir/amdgcn-amdhsa/bin/ranlib
cp -a bin/llvm-mc $install_dir/amdgcn-amdhsa/bin/as
cp -a bin/llvm-nm $install_dir/amdgcn-amdhsa/bin/nm
cp -a bin/lld $install_dir/amdgcn-amdhsa/bin/ld

# Clone gcc, newlib and link newlib in
cd $work_dir
git clone git://sourceware.org/git/newlib-cygwin.git newlib
git clone git://gcc.gnu.org/git/gcc.git gcc
cd gcc
contrib/download_prerequisites
ln -s ../newlib/newlib newlib
cd ..
target=$(gcc/config.guess)

# Build offloading GCC
cd $build_dir
mkdir build-offload-gcc
cd build-offload-gcc
$work_dir/gcc/configure \
  --target=amdgcn-amdhsa --enable-languages=c,lto,fortran \
  --disable-sjlj-exceptions --with-newlib \
  --enable-as-accelerator-for=x86_64-pc-linux-gnu \
  --with-build-time-tools=$install_dir/amdgcn-amdhsa/bin \
  --disable-libquadmath --prefix=$install_dir 
make -j $(nproc) || exit 1
make install || exit 1
cd ..

rm $work_dir/gcc/newlib

# Build host GCC
mkdir build-host-gcc
cd  build-host-gcc
$work_dir/gcc/configure \
    --build=x86_64-pc-linux-gnu --host=x86_64-pc-linux-gnu \
    --target=x86_64-pc-linux-gnu \
    --enable-offload-targets=amdgcn-amdhsa=$install_dir/amdgcn-amdhsa \
    --enable-languages="c,c++,fortran,lto" \
    --disable-multilib --prefix=$install_dir

make -j $(nproc) || exit 1
make install || exit 1
cd ..
