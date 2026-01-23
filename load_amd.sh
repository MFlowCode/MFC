module load python cmake
module load cpe/25.09
module load PrgEnv-amd

AFAR_UMS_BASEDIR="/sw/crusher/ums/compilers/afar"
AFAR_UMS_LATEST=$(ls -d --color=never ${AFAR_UMS_BASEDIR}/*/ | tail -n1)
export OLCF_AFAR_ROOT=${AFAR_UMS_BASEDIR}/"rocm-afar-8873-drop-22.2.0"

export PATH=${OLCF_AFAR_ROOT}/lib/llvm/bin:${PATH}
export LD_LIBRARY_PATH=${OLCF_AFAR_ROOT}/lib:${OLCF_AFAR_ROOT}/lib/llvm/lib:${LD_LIBRARY_PATH}

export CRAY_MPICH_INC="-I${OLCF_AFAR_ROOT}/include/mpich3.4a2"
export CRAY_HIPFORT_INC="-I${OLCF_AFAR_ROOT}/include/hipfort/amdgcn"
export CRAY_HIP_INC="-I${OLCF_AFAR_ROOT}/include/hip"
export CRAY_MPICH_LIB="-L${CRAY_MPICH_PREFIX}/lib \
                        ${CRAY_PMI_POST_LINK_OPTS} \
                        -lmpifort_amd -lmpi_amd -lmpi -lpmi -lpmi2"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${CRAY_LD_LIBRARY_PATH}"
export CMAKE_PREFIX_PATH="${OLCF_AFAR_ROOT}:${CMAKE_PREFIX_PATH}"
export FC="${OLCF_AFAR_ROOT}/bin/amdflang"

unset MPICH_GPU_SUPPORT_ENABLED

# module use /ccs/home/bcornille/afar-drops/modulefiles/Core/
# module load rocm-afar-drop mpich cray-python
# module load cmake

export OMPX_APU_MAPS=0
export HSA_XNACK=0
export LIBOMPTARGET_INFO=0
export AMD_LOG_LEVEL=0
export OFFLOAD_TRACK_ALLOCATION_TRACES=false
