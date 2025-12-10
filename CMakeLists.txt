# SPDX-License-Identifier: MIT

# Developer's note:
#  CMake has extensive documentation available online. Searching "cmake <variable>"
#  or "cmake <function_name>" will return a cmake.org page with more information
#  than you might expect to find or need. Note that any variable or symbol
#  prefixed by CMAKE_ is reserved by CMake and has special meaning.

cmake_minimum_required(VERSION 3.20)


# We include C as a language because - for some reason -
# FIND_LIBRARY_USE_LIB64_PATHS is otherwise ignored.

project(MFC LANGUAGES C CXX Fortran)

# Build options exposed to users and their default values.

option(MFC_MPI           "Build with MPI"                                     ON)
option(MFC_OpenACC       "Build with OpenACC"                                OFF)
option(MFC_OpenMP        "Build with OpenMP"                                 OFF)
option(MFC_GCov          "Build with GCov"                                   OFF)
option(MFC_Unified       "Build with unified CPU & GPU memory (GH-200 only)" OFF)
option(MFC_Fastmath      "Build with -gpu=fastmath on NV GPUs"               OFF)
option(MFC_PRE_PROCESS   "Build pre_process"                                 OFF)
option(MFC_SIMULATION    "Build simulation"                                  OFF)
option(MFC_POST_PROCESS  "Build post_process"                                OFF)
option(MFC_SYSCHECK      "Build syscheck"                                    OFF)
option(MFC_DOCUMENTATION "Build documentation"                               OFF)
option(MFC_ALL           "Build everything"                                  OFF)
option(MFC_SINGLE_PRECISION "Build single precision"                         OFF)
option(MFC_MIXED_PRECISION "Build mixed precision"                           OFF)

if (MFC_ALL)
    set(MFC_PRE_PROCESS   ON FORCE)
    set(MFC_SIMULATION    ON FORCE)
    set(MFC_POST_PROCESS  ON FORCE)
    set(MFC_DOCUMENTATION ON FORCE)
endif()

if (MFC_SINGLE_PRECISION)
    add_compile_definitions(MFC_SINGLE_PRECISION)
else()
    add_compile_definitions(MFC_DOUBLE_PRECISION)
endif()

if (MFC_MIXED_PRECISION)
    add_compile_definitions(MFC_MIXED_PRECISION)
endif()


# CMake Library Imports

include(GNUInstallDirs)
include(ExternalProject)
include(CheckIPOSupported)
include(CMakeParseArguments)
include(CheckFortranCompilerFlag)


# Check Compiler Support: Some compilers, especially old ones, are known not to
# be able to compile MFC. We throw our own error early - at configure time - in
# such cases.

set(__err_msg "\
CMake detected the ${CMAKE_Fortran_COMPILER_ID} Fortran compiler \
v${CMAKE_Fortran_COMPILER_VERSION}. If you intended to use a different \
compiler (or a different version thereof), please:\n\
 - Install the compiler or load its module. (e.g. module load gcc/10.1)\n\
 - Set/Export the C, CXX, and FC environment variables. (e.g. 'export CC=gcc', \
 'export CXX=g++', and 'export FC=gfortran'.\n\
 - If using mfc.sh, delete the build/<code name> directory and try again. (e.g. 'rm -rf build/pre_process')")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5)
        message(FATAL_ERROR "ERROR: GNU v5.0 or newer is required to build MFC.\n${__err_msg}")
    endif()
    if (MFC_OpenACC)
        message(FATAL_ERROR "ERROR: MFC with GPU processing is not currently compatible with GNU compilers. Please use NVIDIA or Cray compilers.\n${__err_msg}")
    endif()
elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 21.7)
        message(FATAL_ERROR "ERROR: When using NVHPC, v21.7 or newer is required to build MFC.\n${__err_msg}")
    endif()

    if ((CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 24.5) AND (CMAKE_BUILD_TYPE STREQUAL "Debug") AND MFC_OpenACC)
        message(FATAL_ERROR "ERROR: When using NVHPC, MFC with Debug and GPU options requires NVHPC v24.5 or newer.\n${__err_msg}")
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
    message(FATAL_ERROR "ERROR: MFC does not support the Apple Clang compilers. Please consult the documentation.\n${__err_msg}")
endif()

# Fypp is required to build MFC. We try to locate it. The path to the binary is
# returned in FYPP_EXE upon success. This is used later.

find_program(FYPP_EXE fypp REQUIRED)


# Miscellaneous Configuration:
# *  Explicitly link to -ldl (or system equivalent)
# *  Request that FIND_LIBRARY searches <prefix>/lib/ and <prefix>/lib64/
# *  Let FindXXX use custom scripts from toolchain/cmake/.

link_libraries("${CMAKE_DL_LIBS}")
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake/regular")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake/cce")
endif()


# Compiler Flags: Here, we specify our own compiler flags for both release and
# debug builds. These include optimization and debug flags, as well as some that
# are required for a successful build of MFC.

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
    )

    if (MFC_GCov)

        # Warning present in gcc versions >= 12 that is treated as an error
        # This flag doesn't exist in gcc versions < 12
        if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 12)
            add_compile_options(
                -Wno-error=coverage-invalid-line-number
            )
        endif()

        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-fprofile-arcs>
            $<$<COMPILE_LANGUAGE:Fortran>:-ftest-coverage>
            $<$<COMPILE_LANGUAGE:Fortran>:-O1>
	    )

        add_link_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-lgcov>
            $<$<COMPILE_LANGUAGE:Fortran>:--coverage>
        )
    endif()

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
            -Wall
            -Wextra
            -fcheck=all,no-array-temps
            -fbacktrace
            -fimplicit-none
            -fsignaling-nans
            -finit-real=snan
            -finit-integer=-99999999
            -Wintrinsic-shadow
            -Wunderflow
            -Wrealloc-lhs
            -Wsurprising
	    )
    endif()

    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-fallow-invalid-boz>
            $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>
            $<$<COMPILE_LANGUAGE:Fortran>:-fcheck=bounds>
        )
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    add_compile_options(
        "SHELL:-M 296,878,1391,1069,5025"
        "SHELL:-h static" "SHELL:-h keepfiles"
        "SHELL:-h acc_model=auto_async_none"
        "SHELL: -h acc_model=no_fast_addr"
        "SHELL: -h list=adm"
        "SHELL: -munsafe-fp-atomics" # Not unsafe for operations we do
    )

    add_link_options("SHELL:-hkeepfiles")

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
                "SHELL:-h acc_model=auto_async_none"
                "SHELL: -h acc_model=no_fast_addr"
                "SHELL: -K trap=fp" "SHELL: -G2"
        )
        add_link_options("SHELL: -K trap=fp" "SHELL: -G2")
    endif()

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>
        $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8>
    )
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-free>)

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -Og -traceback -debug)
    endif()
elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
        $<$<COMPILE_LANGUAGE:Fortran>:-Minfo=inline>
        $<$<COMPILE_LANGUAGE:Fortran>:-Minfo=accel>
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-O0>
            $<$<COMPILE_LANGUAGE:Fortran>:-C>
            $<$<COMPILE_LANGUAGE:Fortran>:-g>
            $<$<COMPILE_LANGUAGE:Fortran>:-traceback>
            $<$<COMPILE_LANGUAGE:Fortran>:-Minform=inform>
            $<$<COMPILE_LANGUAGE:Fortran>:-Mbounds>
        )
    endif()

    if (DEFINED ENV{MFC_CUDA_CC})
        string(REGEX MATCHALL "[0-9]+" MFC_CUDA_CC $ENV{MFC_CUDA_CC})
        message(STATUS "Found $MFC_CUDA_CC specified. GPU code will be generated for compute capability(ies) ${MFC_CUDA_CC}.")
    endif()
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    # Processor tuning: Check if we can target the host's native CPU's ISA.
    CHECK_FORTRAN_COMPILER_FLAG("-march=native" SUPPORTS_MARCH_NATIVE)
    if (SUPPORTS_MARCH_NATIVE)
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-march=native>)
    else()
    	CHECK_FORTRAN_COMPILER_FLAG("-mcpu=native" SUPPORTS_MCPU_NATIVE)
        if (SUPPORTS_MCPU_NATIVE)
            add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-mcpu=native>)
        endif()
    endif()

    # Enable LTO/IPO if supported
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
        if (MFC_Unified)
            message(STATUS "LTO/IPO is not available with NVHPC using Unified Memory")
        elseif(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "23.11")
            message(STATUS "LTO/IPO is not supported in NVHPC Version < 23.11. Use a newer version of NVHPC for best performance.")
        else()
            message(STATUS "Performing IPO using -Mextract followed by -Minline")
            set(NVHPC_USE_TWO_PASS_IPO TRUE)
        endif()
    else()
        CHECK_IPO_SUPPORTED(RESULT SUPPORTS_IPO OUTPUT IPO_ERROR)
        if (SUPPORTS_IPO)
            message(STATUS "Enabled IPO / LTO")
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
	    else()
            message(STATUS "IPO / LTO is NOT available")
        endif()
    endif()
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(MFC_DEBUG)
endif()



# HANDLE_SOURCES: Given a target (herein <target>):
#
# *  Locate all source files for <target> of the type
#
#               src/[<target>,common]/[.,include]/*.[f90,fpp].
#
# *  For each .fpp file found with filepath <dirpath>/<filename>.fpp, using a
#    custom command, instruct CMake how to generate a file with path
#
#                     src/<target>/fypp/<filename>.f90
#
#    by running Fypp on <dirpath>/<filename>.fpp. It is important to understand
#    that this does not actually run the pre-processor. Rather, it instructs
#    CMake what to do when it finds a src/<target>/fypp/<filename>.f90 path
#    in the source list for a target. Thus, an association is made from an .f90
#    file to its corresponding .fpp file (if applicable) even though the
#    generation is of the form .fpp -> .f90.
#
#    This design has one limitation: If an .fpp file depends on another, for
#    example if it '#:include's it and uses a macro defined in it, then the
#    dependency will not be tracked. A modification to the .fpp file it depends
#    on will not trigger a re-run of Fypp on the .fpp file that depends on it.
#    As a compromise, both in speed and complexity, all .f90 files generated
#    from .fpp files are re-generated not only when their corresponding .fpp
#    file is modified, but also when any file with filepath of the form
#
#                     src/[<target>,common]/include/*.fpp
#
#    is modified. This is a reasonable compromise as modifications to .fpp files
#    in the include directories will be rare - by design. Other approaches would
#    have required a more complex CMakeLists.txt file (perhaps parsing the .fpp
#    files to determine their dependencies) or the endurment of longer
#    compilation times (by way of re-running Fypp on all .fpp files every time
#    one of them is modified).
#
#    .fpp files in src/common are treated as if they were in src/<target> (not
#    pre-processed to src/common/fypp/) so as not to clash with other targets'
#    .fpp files (this has caused problems in the past).
#
# *  Export, in the variable <target>_SRCs, a list of all source files (.f90)
#    that would compile to produce <target>. If <target> includes .fpp files,
#    then the list will include the paths to the corresponding .f90 files that
#    Fypp would generate from the .fpp files.
#
# This design allows us to be flexible in our use of Fypp as we don't have to
# worry about running the pre-processor on .fpp files when we create executables
# and generate documentation. Instead, we can simply include the list of .f90
# files that will eventually be used to compile <target>.

macro(HANDLE_SOURCES target useCommon)
    set(${target}_DIR "${CMAKE_SOURCE_DIR}/src/${target}")
    set(common_DIR    "${CMAKE_SOURCE_DIR}/src/common")

    string(TOUPPER ${target} ${target}_UPPER)

    # Gather:
    # *          src/[<target>,(common)]/*.f90
    # * (if any) <build>/modules/<target>/*.f90
    file(GLOB ${target}_F90s CONFIGURE_DEPENDS "${${target}_DIR}/*.f90"
                                               "${CMAKE_BINARY_DIR}/modules/${target}/*.f90")
    set(${target}_SRCs ${${target}_F90s})
    if (${useCommon})
        file(GLOB common_F90s CONFIGURE_DEPENDS "${common_DIR}/*.f90")
        list(APPEND ${target}_SRCs ${common_F90s})
    endif()

    # Gather:
    # *          src/[<target>,(common)]/*.fpp]
    # * (if any) <build>/modules/<target>/*.fpp
    file(GLOB ${target}_FPPs CONFIGURE_DEPENDS "${${target}_DIR}/*.fpp"
                                               "${CMAKE_BINARY_DIR}/modules/${target}/*.fpp")
    if (${useCommon})
        file(GLOB common_FPPs CONFIGURE_DEPENDS "${common_DIR}/*.fpp")

        # If we're building post_process, exclude m_compute_levelset.fpp
        if("${target}" STREQUAL "post_process")
            list(FILTER common_FPPs EXCLUDE REGEX ".*/m_compute_levelset\.fpp$")
            list(FILTER common_FPPs EXCLUDE REGEX ".*/m_ib_patches\.fpp$")
        endif()

        list(APPEND ${target}_FPPs ${common_FPPs})
    endif()

    # Gather:
    # *          src/[<target>,common]/include/*.fpp
    # * (if any) <build>/include/<target>/*.fpp
    file(GLOB ${target}_incs CONFIGURE_DEPENDS "${${target}_DIR}/include/*.fpp"
                                               "${CMAKE_BINARY_DIR}/include/${target}/*.fpp")

    if (${useCommon})
        file(GLOB common_incs CONFIGURE_DEPENDS "${common_DIR}/include/*.fpp")
        list(APPEND ${target}_incs ${common_incs})
    endif()

    # /path/to/*.fpp (used by <target>) -> <build>/fypp/<target>/*.f90
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/fypp/${target}")
    foreach(fpp ${${target}_FPPs})
        cmake_path(GET fpp FILENAME fpp_filename)
        set(f90 "${CMAKE_BINARY_DIR}/fypp/${target}/${fpp_filename}.f90")

        add_custom_command(
            OUTPUT   ${f90}
            COMMAND  ${FYPP_EXE} -m re
                                 -I "${CMAKE_BINARY_DIR}/include/${target}"
                                 -I "${${target}_DIR}/include"
                                 -I "${common_DIR}/include"
                                 -I "${common_DIR}"
                                 -D MFC_${CMAKE_Fortran_COMPILER_ID}
                                 -D MFC_${${target}_UPPER}
                                 -D MFC_COMPILER="${CMAKE_Fortran_COMPILER_ID}"
                                 -D MFC_CASE_OPTIMIZATION=False
                                 -D chemistry=False
                                 --line-numbering
                                 --no-folding
								 --line-length=999
		 						 --line-numbering-mode=nocontlines
                                 "${fpp}" "${f90}"
            DEPENDS  "${fpp};${${target}_incs}"
            COMMENT  "Preprocessing (Fypp) ${fpp_filename}"
            VERBATIM
        )

        list(APPEND ${target}_SRCs ${f90})
    endforeach()
endmacro()


HANDLE_SOURCES(pre_process ON)
HANDLE_SOURCES(simulation ON)
HANDLE_SOURCES(post_process ON)
HANDLE_SOURCES(syscheck OFF)


# MFC_SETUP_TARGET: Given a target (herein <target>), this macro creates a new
# executable <target> with the appropriate sources, compiler definitions, and
# linked libraries (assuming HANDLE_SOURCES was called on <target>).
#
# Inputs:
#  * TARGET   Target name (<target>)
#  * SOURCES  List of source files (.f90)
#  * OpenACC  (optional) Can be compiled with OpenACC.
#  * MPI      (optional) Can be compiled with MPI.
#  * SILO     (optional) Should be linked with SILO.
#  * HDF5     (optional) Should be linked with HDF5.
#  * FFTW     (optional) Should be linked with an FFTW-like library (fftw/cufftw),
#             depending on whether OpenACC is enabled and which compiler is
#             being used.
#  * LAPACK   (optional) Should be linked with LAPACK

function(MFC_SETUP_TARGET)
    cmake_parse_arguments(ARGS "OpenACC;MPI;SILO;HDF5;FFTW;LAPACK;OpenMP" "TARGET" "SOURCES" ${ARGN})

    add_executable(${ARGS_TARGET} ${ARGS_SOURCES})
    set(IPO_TARGETS ${ARGS_TARGET})
    # Here we need to split into "library" and "executable" to perform IPO on the NVIDIA compiler.
    # A little hacky, but it *is* an edge-case for *one* compiler.
    if (NVHPC_USE_TWO_PASS_IPO AND NOT(MFC_OpenMP AND ARGS_OpenMP))
        add_library(${ARGS_TARGET}_lib OBJECT ${ARGS_SOURCES})
        target_compile_options(${ARGS_TARGET}_lib PRIVATE
		    $<$<COMPILE_LANGUAGE:Fortran>:-Mextract=lib:${ARGS_TARGET}_lib>
	        $<$<COMPILE_LANGUAGE:Fortran>:-Minline>
    	)
        add_dependencies(${ARGS_TARGET} ${ARGS_TARGET}_lib)
        target_compile_options(${ARGS_TARGET} PRIVATE -Minline=lib:${ARGS_TARGET}_lib)
        list(PREPEND IPO_TARGETS ${ARGS_TARGET}_lib)
    endif()

    foreach (a_target ${IPO_TARGETS})
        set_target_properties(${a_target} PROPERTIES Fortran_PREPROCESS ON)

        target_include_directories(${a_target} PRIVATE
            "${CMAKE_SOURCE_DIR}/src/common"
            "${CMAKE_SOURCE_DIR}/src/common/include"
            "${CMAKE_SOURCE_DIR}/src/${ARGS_TARGET}")

        if (EXISTS "${CMAKE_SOURCE_DIR}/src/${ARGS_TARGET}/include")
            target_include_directories(${a_target} PRIVATE
                "${CMAKE_SOURCE_DIR}/src/${ARGS_TARGET}/include")
        endif()

        string(TOUPPER "${ARGS_TARGET}" ${ARGS_TARGET}_UPPER)
        target_compile_definitions(
            ${a_target} PRIVATE MFC_${CMAKE_Fortran_COMPILER_ID}
                                MFC_${${ARGS_TARGET}_UPPER}
        )

        if (MFC_MPI AND ARGS_MPI)
            find_package(MPI COMPONENTS Fortran REQUIRED)

            target_compile_definitions(${a_target} PRIVATE MFC_MPI)
            target_link_libraries     (${a_target} PRIVATE MPI::MPI_Fortran)
        endif()

        if (ARGS_SILO)
            find_package(SILO REQUIRED)
            target_link_libraries(${a_target} PRIVATE SILO::SILO stdc++)
        endif()

        if (ARGS_HDF5)
            find_package(HDF5 REQUIRED)
            target_link_libraries(${a_target} PRIVATE HDF5::HDF5)
        endif()

        if (ARGS_FFTW)
            if ((MFC_OpenACC AND ARGS_OpenACC) OR (MFC_OpenMP AND ARGS_OpenMP))
                if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
                    find_package(CUDAToolkit REQUIRED)
                    target_link_libraries(${a_target} PRIVATE CUDA::cudart CUDA::cufft)
                else()
                    find_package(hipfort COMPONENTS hipfft CONFIG REQUIRED)
                    target_link_libraries(${a_target} PRIVATE hipfort::hipfft)
                endif()
            else()
                find_package(FFTW REQUIRED)
                target_link_libraries(${a_target} PRIVATE FFTW::FFTW)
            endif()
        endif()

        if (ARGS_LAPACK)
            find_package(LAPACK REQUIRED)
            target_link_libraries(${a_target} PRIVATE LAPACK::LAPACK)
        endif()

        if ((MFC_OpenACC AND ARGS_OpenACC) OR (MFC_OpenMP AND ARGS_OpenMP))
            if ((MFC_OpenACC AND ARGS_OpenACC))
                find_package(OpenACC)

                # This should be equivalent to if (NOT OpenACC_FC_FOUND)
                if (NOT TARGET OpenACC::OpenACC_Fortran)
                    message(FATAL_ERROR "OpenACC + Fortran is unsupported.")
                endif()

                target_link_libraries(${a_target} PRIVATE OpenACC::OpenACC_Fortran)
                target_compile_definitions(${a_target} PRIVATE MFC_OpenACC MFC_GPU)
            elseif((MFC_OpenMP AND ARGS_OpenMP))
                find_package(OpenMP)

                # This should be equivalent to if (NOT OpenACC_FC_FOUND)
                if (NOT TARGET OpenMP::OpenMP_Fortran)
                    message(FATAL_ERROR "OpenMP + Fortran is unsupported.")
                endif()
                set(ENV{OMP_TARGET_OFFLOAD} [MANDATORY])
                # target_link_libraries(${a_target} PRIVATE OpenMP::OpenMP_Fortran)
                target_compile_definitions(${a_target} PRIVATE MFC_OpenMP MFC_GPU)

                if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
                    target_compile_options(${a_target} PRIVATE "-mp=gpu" "-Minfo=mp")
                    target_link_options(${a_target} PRIVATE "-mp=gpu")
                    set_target_properties(${a_target} PROPERTIES Fortran_FLAGS "-mp=gpu -gpu=ccall")
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
                    target_compile_options(${a_target} PRIVATE -fopenmp -fopenmp-targets=spir64)
                    target_link_options(${a_target} PRIVATE -fopenmp -fopenmp-targets=spir64)
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
                    target_compile_options(${a_target} PRIVATE -fopenmp)
                    target_link_options(${a_target} PRIVATE -fopenmp)
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
                    target_compile_options(${a_target} PRIVATE -fopenmp --offload-arch=gfx90a)
                    target_link_options(${a_target} PRIVATE -fopenmp --offload-arch=gfx90a)
                endif()
            endif()

            if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
                # FIXME: This should work with other cards than gfx90a ones.
                target_compile_options(${a_target} PRIVATE
                    "-foffload=amdgcn-amdhsa='-march=gfx90a'"
                    "-foffload-options=-lgfortran\ -lm"
                    "-fno-exceptions")
                if (MFC_Fastmath)
                    message(WARNING "--fastmath has no effect with the GNU compiler")
                endif()
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
                foreach (cc ${MFC_CUDA_CC})
                    target_compile_options(${a_target}
                        PRIVATE -gpu=cc${cc}
                    )
                endforeach()

                target_compile_options(${a_target}
                    PRIVATE -gpu=keep,ptxinfo,lineinfo
                )

                if (MFC_Fastmath)
                    target_compile_options(${a_target}
                        PRIVATE -gpu=fastmath
                    )
                endif()

                # GH-200 Unified Memory Support
                if (MFC_Unified)
                    target_compile_options(${ARGS_TARGET}
                        PRIVATE -gpu=mem:unified:managedalloc -cuda
                    )
                    # "This option must appear in both the compile and link lines" -- NVHPC Docs
                    target_link_options(${ARGS_TARGET}
                        PRIVATE -gpu=mem:unified:managedalloc -cuda
                    )
                endif()

                if (CMAKE_BUILD_TYPE STREQUAL "Debug")
                    target_compile_options(${a_target}
                        PRIVATE -gpu=debug
                    )
                endif()
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
                # Frontier Unified Memory Support
                if (MFC_Unified)
                    target_compile_options(${ARGS_TARGET}
                        PRIVATE -DFRONTIER_UNIFIED)
                endif()

                if (MFC_Fastmath)
                    message(WARNING "--fastmath has no effect with the CCE")
                endif()

                find_package(hipfort COMPONENTS hip CONFIG REQUIRED)
                target_link_libraries(${a_target} PRIVATE hipfort::hip hipfort::hipfort-amdgcn)
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")

                if (MFC_Unified)
                    target_compile_options(${ARGS_TARGET}
                        PRIVATE -DFRONTIER_UNIFIED)
                endif()

    		find_package(hipfort COMPONENTS hip CONFIG REQUIRED)
                target_link_libraries(${a_target} PRIVATE hipfort::hip hipfort::hipfort-amdgcn flang_rt.hostdevice)
            endif()
        elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
            target_compile_options(${a_target} PRIVATE "SHELL:-h noacc" "SHELL:-x acc")
            if (MFC_Fastmath)
                message(WARNING "--fastmath has no effect with the CCE")
            endif()
        endif()

        if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
            find_package(CUDAToolkit REQUIRED)
            if (TARGET CUDA::nvToolsExt) # CUDA <= 12.8
                target_link_libraries(${a_target} PRIVATE CUDA::nvToolsExt)
            else() # CUDA >= 12.9
                target_link_libraries(${a_target} PRIVATE nvhpcwrapnvtx)
                target_link_options(${a_target} PRIVATE "-cudalib=nvtx3")
            endif()
        endif()
    endforeach()

    install(TARGETS ${ARGS_TARGET} RUNTIME DESTINATION bin)
endfunction()


if (MFC_PRE_PROCESS)
    MFC_SETUP_TARGET(TARGET  pre_process
                     SOURCES "${pre_process_SRCs}"
                     MPI)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        target_compile_options(pre_process PRIVATE -hfp0)
    endif()
endif()

if (MFC_SIMULATION)
    MFC_SETUP_TARGET(TARGET  simulation
                     SOURCES "${simulation_SRCs}"
                     MPI FFTW OpenACC OpenMP)
endif()

if (MFC_POST_PROCESS)
    MFC_SETUP_TARGET(TARGET  post_process
                     SOURCES "${post_process_SRCs}"
                     MPI SILO HDF5 FFTW LAPACK)

    # -O0 is in response to https://github.com/MFlowCode/MFC-develop/issues/95
    target_compile_options(post_process PRIVATE -O0)
endif()

if (MFC_SYSCHECK)
    MFC_SETUP_TARGET(TARGET  syscheck
                     SOURCES "${syscheck_SRCs}"
                     MPI OpenACC OpenMP)
endif()

if (MFC_DOCUMENTATION)
    # Files in examples/ are used to generate docs/documentation/examples.md
    file(GLOB_RECURSE examples_DOCs CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/examples/*")

    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/examples.md"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/examples.sh;${examples_DOCs}"
        COMMAND "bash" "${CMAKE_CURRENT_SOURCE_DIR}/docs/examples.sh"
                       "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Generating examples.md"
        VERBATIM
    )

    # Generate case_constraints.md from case_validator.py and examples/
    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/case_constraints.md"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/gen_case_constraints_docs.py"
                "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/case_validator.py"
                "${examples_DOCs}"
        COMMAND "bash" "${CMAKE_CURRENT_SOURCE_DIR}/docs/gen_constraints.sh"
                       "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Generating case_constraints.md"
        VERBATIM
    )

    file(GLOB common_DOCs CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/*")

    # GEN_DOCS: Given a target name (herein <target>), this macro sets up a
    # target, <target>_docs, that generates documentation for <target> using
    # Doxygen. It is then added as a dependency of the documentation target.
    # Doxygen outputs HTML to the binary (build) directory. It is installed into
    # its own directory: <prefix>/docs/mfc/<target>. We use configure_file to
    # generate the Doxyfile for each target using paths that only apply to that
    # target (or when relative paths are impractical), specified by CMake.

    macro(GEN_DOCS target name)
        set(DOXYGEN_PROJECT_NAME     "\"${name}\"")
        set(DOXYGEN_INPUT            "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/${target}\"")

        foreach (f90 ${${target}_SRCs})
            set(DOXYGEN_INPUT "${DOXYGEN_INPUT} \"${f90}\"")
        endforeach()

        set(DOXYGEN_HTML_HEADER      "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/header.html\"")
        set(DOXYGEN_HTML_FOOTER      "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/footer.html\"")
        set(DOXYGEN_HTML_OUTPUT      "\"${CMAKE_CURRENT_BINARY_DIR}/${target}\"")
        set(DOXYGEN_MATHJAX_CODEFILE "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/config.js\"")
        set(DOXYGEN_PROJECT_LOGO     "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/res/icon.ico\"")
        set(DOXYGEN_IMAGE_PATH       "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/res\"\
                                      \"${CMAKE_CURRENT_SOURCE_DIR}/docs/${target}\"")

        file(MAKE_DIRECTORY "${DOXYGEN_OUTPUT_DIRECTORY}")

        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in"
            "${CMAKE_CURRENT_BINARY_DIR}/${target}-Doxyfile" @ONLY)

        set(opt_example_dependency "")
        set(opt_constraints_dependency "")
        if (${target} STREQUAL documentation)
            set(opt_example_dependency "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/examples.md")
            set(opt_constraints_dependency "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/case_constraints.md")
        endif()

        file(GLOB_RECURSE ${target}_DOCs CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/${target}/*")
        list(APPEND ${target}_DOCs "${common_DOCs}")

        add_custom_command(
            OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/${target}/html/index.html"
            DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${target}-Doxyfile"
                    "${opt_example_dependency}"
                    "${opt_constraints_dependency}"
                    "${${target}_SRCs}" "${${target}_DOCs}"
            COMMAND "${DOXYGEN_EXECUTABLE}" "${target}-Doxyfile"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
            COMMENT "${target}: Generating documentation"
        )

        add_custom_target(
            "${target}_doxygen" ALL
            DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${target}/html/index.html"
        )

        install(DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/${target}"
                DESTINATION "docs/mfc")

        add_dependencies("${target}_doxygen" doxygen-awesome-css)
        add_dependencies(documentation "${target}_doxygen")
    endmacro()

    add_custom_target(documentation)

    find_package(Doxygen REQUIRED dot REQUIRED)

    # > Fetch CSS Theme
    ExternalProject_Add(doxygen-awesome-css
        PREFIX            doxygen-awesome-css
        GIT_REPOSITORY    "https://github.com/jothepro/doxygen-awesome-css"
        GIT_TAG           "df88fe4fdd97714fadfd3ef17de0b4401f804052"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     ""
        INSTALL_COMMAND   ""
    )

    set(theme_dirpath
        "${CMAKE_CURRENT_BINARY_DIR}/doxygen-awesome-css/src/doxygen-awesome-css")

    set(DOXYGEN_HTML_EXTRA_STYLESHEET
        "\"${theme_dirpath}/doxygen-awesome.css\"\
         \"${theme_dirpath}/doxygen-awesome-sidebar-only.css\"")

    # > Generate Documentation & Landing Page
    GEN_DOCS(pre_process   "MFC: Pre-Process")
    GEN_DOCS(simulation    "MFC: Simulation")
    GEN_DOCS(post_process  "MFC: Post-Process")
    GEN_DOCS(documentation "MFC")

    # > Copy Resources (main landing page & assets)
    install(DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/docs/res"
            DESTINATION "docs/mfc")

    install(FILES       "${CMAKE_CURRENT_SOURCE_DIR}/docs/robots.txt"
            DESTINATION "docs/mfc")

    install(FILES       "${CMAKE_CURRENT_SOURCE_DIR}/docs/index.html"
            DESTINATION "docs/mfc")
endif()

site_name(SITE_NAME)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake/configuration.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/configuration.txt")
