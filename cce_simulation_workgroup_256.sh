#! /bin/bash

CCE_LLVM_PATH=${CRAY_CCE_CLANGSHARE}/../

WGSIZE=256

echo "Disassembling"
${CCE_LLVM_PATH}/bin/llvm-dis build/simulation/simulation-cce-openmp-pre-llc.bc

echo "Globally setting amdgpu-flat-work-group-size size to 1,$WGSIZE"
sed "s/\"amdgpu-flat-work-group-size\"\=\"1,1024\"/\"amdgpu-flat-work-group-size\"\=\"1,${WGSIZE}\"/g" build/simulation/simulation-cce-openmp-pre-llc.ll > build/simulation/simulation-cce-openmp-pre-llc-wg${WGSIZE}.ll

echo "Invoking LLC to compile"
${CCE_LLVM_PATH}/bin/llc -mtriple=amdgcn-amd-amdhsa -disable-promote-alloca-to-lds -mcpu=gfx90a -amdgpu-dump-hsa-metadata build/simulation/simulation-cce-openmp-pre-llc-wg${WGSIZE}.ll -filetype=obj -o build/simulation/simulation-cce-openmp__llc_wg${WGSIZE}.amdgpu

echo "Linking to a CCE Offload module"
${CCE_LLVM_PATH}/bin/lld  -flavor gnu --no-undefined -shared -o build/simulation/simulation-wg${WGSIZE}.lld.exe build/simulation/simulation-cce-openmp__llc_wg${WGSIZE}.amdgpu

echo "Now "
echo "export CRAY_ACC_MODULE=${PWD}/build/simulation/simulation-wg${WGSIZE}.lld.exe"
echo "to use the new GPU offload code."
echo "To use the original build"
echo "unset CRAY_ACC_MODULE"
