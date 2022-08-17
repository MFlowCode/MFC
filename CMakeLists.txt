CMAKE_MINIMUM_REQUIRED(VERSION 3.17)

PROJECT(MFC LANGUAGES C CXX Fortran)


# === Options
OPTION(MFC_WITH_MPI           "Build with MPI"                                   ON)
OPTION(MFC_WITH_OPEN_ACC      "Build with OpenACC"                              OFF)
OPTION(MFC_BUILD_PRE_PROCESS  "Build pre_process"                               OFF)
OPTION(MFC_BUILD_SIMULATION   "Build simulation"                                OFF)
OPTION(MFC_BUILD_POST_PROCESS "Build post_process"                              OFF)
OPTION(MFC_BUILD_ALL          "Build pre_process, simulation, and post_process" OFF)

IF (MFC_BUILD_ALL)
    SET(MFC_BUILD_PRE_PROCESS  ON FORCE)
    SET(MFC_BUILD_SIMULATION   ON FORCE)
    SET(MFC_BUILD_POST_PROCESS ON FORCE)
ENDIF()


# === Imports
INCLUDE(GNUInstallDirs)
INCLUDE(CheckIPOSupported)
INCLUDE(CheckFortranCompilerFlag)


# === Check Compiler Support & Tools
# === === Compiler Support
IF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    
    IF (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5)
        MESSAGE(FATAL_ERROR "GNU v5.0 or newer is required to build MFC. "
                            "[current: ${CMAKE_Fortran_COMPILER_VERSION}]")
    ENDIF()

ELSEIF ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR
        (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    
    IF (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 21.7)
        MESSAGE(FATAL_ERROR "NVHPC v21.7 or newer is required to build MFC. "
                            "[current: ${CMAKE_Fortran_COMPILER_VERSION}]")
    ENDIF()

ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "AppleClang")
    
    MESSAGE(FATAL_ERROR "MFC does not support the Apple Clang compilers. "
                        "Please consult the README for more details.")

ENDIF()

# === === Find Fypp
FIND_PROGRAM(FYPP_PATH fypp REQUIRED)


# === Miscellaneous Configuration
# === === Enable C-Preprocessor for Fortran files
SET(CMAKE_Fortran_PREPROCESS ON)   
# === === Explicitly link to -ldl (or system equivalent)
LINK_LIBRARIES("${CMAKE_DL_LIBS}")
# === === Request that FIND_LIBRARY searches lib/ and lib64/
SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)
# === === Let FindXXX use custom scripts from cmake/ 
LIST(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake")


# === Compiler Flags
IF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

    ADD_COMPILE_OPTIONS(
        -w                      -fdefault-real-8
        -ffree-line-length-none -fdefault-double-8
    )

    IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
        ADD_COMPILE_OPTIONS(
            -fcheck=all     -fbacktrace
            -fimplicit-none -ffpe-trap=invalid,zero,overflow,underflow
        )
    ENDIF()

    IF (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
        ADD_COMPILE_OPTIONS(-fallow-invalid-boz -fallow-argument-mismatch)
    ENDIF()

ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")

    ADD_COMPILE_OPTIONS(-Mfreeform -Mpreprocess -fdefault-real-8)

ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

    ADD_COMPILE_OPTIONS(-free)

ELSEIF ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR
        (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))

    ADD_COMPILE_OPTIONS(-r8 -cpp -Mfreeform -Minfo=accel -Mr8intrinsics)
    
    IF (MFC_WITH_OPEN_ACC)
        IF     (CMAKE_BUILD_TYPE STREQUAL "Release")
            ADD_COMPILE_OPTIONS(-gpu=keep,ptxinfo,lineinfo)
        ELSEIF (CMAKE_BUILD_TYPE STREQUAL "Debug")
            ADD_COMPILE_OPTIONS(-gpu=keep,ptxinfo,lineinfo,autocompare,debug)
        ENDIF()
    ENDIF()

ENDIF()


# === Release Optimizations
IF (CMAKE_BUILD_TYPE STREQUAL "Release")

    # === === Processor tuning
    CHECK_FORTRAN_COMPILER_FLAG("-march=native" SUPPORTS_MARCH_NATIVE)
    IF (SUPPORTS_MARCH_NATIVE)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
    ELSE()
    	CHECK_FORTRAN_COMPILER_FLAG("-mcpu=native" SUPPORTS_MCPU_NATIVE)
        IF (SUPPORTS_MCPU_NATIVE)
            SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcpu=native")
        ENDIF()
    ENDIF()


    # === === Enable LTO/IPO if supported
    CHECK_IPO_SUPPORTED(RESULT SUPPORTS_IPO OUTPUT IPO_ERROR)
    IF (SUPPORTS_IPO)
        MESSAGE(STATUS "IPO / LTO enabled")
        SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    ENDIF()

ENDIF()


# === Locate Libraries
# === === MPI
IF (MFC_WITH_MPI)
    FIND_PACKAGE(MPI COMPONENTS Fortran REQUIRED)
ENDIF()

# === === CUDAToolkit
IF (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    FIND_PACKAGE(CUDAToolkit REQUIRED)
ENDIF()


# === === Enable OpenACC
IF (MFC_WITH_OPEN_ACC)
    FIND_PACKAGE(OpenACC REQUIRED)
ENDIF()


# === === Dependencies (FFTW, HDF5, and SILO)
IF ((MFC_BUILD_SIMULATION AND (NOT MFC_WITH_OPEN_ACC)) OR MFC_BUILD_POST_PROCESS)
    FIND_PACKAGE(FFTW REQUIRED)
ENDIF()

IF (MFC_BUILD_POST_PROCESS)
    FIND_PACKAGE(HDF5 REQUIRED)
    FIND_PACKAGE(SILO REQUIRED)
ENDIF()


# This macro locates Fypp files in src/<target>_code/ and instructs CMake how to
# preprocess them to generate their .f90 counterparts, should they be requested
# by an executable or library. They are preprocessed again if they are modified
# or if src/common_code/case.fpp changes because some #:include it, and produce
# a different .f90 file depending on case.fpp's contents. This macro sets the
# <target>_srcs CMake cache variaible the full list of .f90 files, including the
# preprocessed Fypp files.
MACRO(HANDLE_FYPP target)
    FILE(GLOB ${target}_f90s "${CMAKE_CURRENT_SOURCE_DIR}/src/${target}_code/*.f90")
    FILE(GLOB ${target}_fpps "${CMAKE_CURRENT_SOURCE_DIR}/src/${target}_code/*.fpp")

    # Create src/<target>_code/autogen folder in which to generate Fypp'd .f90 files.
    FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/${target}_code/autogen")

    SET(${target}_srcs "${${target}_f90s}")

    FOREACH(fpp_filepath IN LISTS ${target}_fpps)
        STRING(REGEX REPLACE ".fpp\$" ".f90" f90_filepath "${fpp_filepath}")

        GET_FILENAME_COMPONENT(fpp_filename "${fpp_filepath}" NAME)
        GET_FILENAME_COMPONENT(f90_dirpath  "${f90_filepath}" DIRECTORY)
        GET_FILENAME_COMPONENT(f90_filename "${f90_filepath}" NAME)

        SET(f90_filepath "${f90_dirpath}/autogen/${f90_filename}")

        ADD_CUSTOM_COMMAND(
            OUTPUT   "${f90_filepath}"
            COMMAND  "fypp" "${fpp_filepath}" "${f90_filepath}"
            DEPENDS  "${fpp_filepath}" "${CMAKE_CURRENT_SOURCE_DIR}/src/common_code/case.fpp"
            COMMENT  "Preprocessing ${fpp_filename}"
            VERBATIM
        )

        LIST(APPEND ${target}_srcs "${f90_filepath}")
    ENDFOREACH(fpp_filepath)

ENDMACRO()


# === Handle src/ Folder
# === === src/common_code
HANDLE_FYPP(common)


# === === src/pre_process_code
IF (MFC_BUILD_PRE_PROCESS)
    HANDLE_FYPP(pre_process)

    ADD_EXECUTABLE(pre_process "${pre_process_srcs}" "${common_srcs}")


    IF (MFC_WITH_MPI)
        TARGET_COMPILE_DEFINITIONS(pre_process PRIVATE MFC_MPI)
        TARGET_LINK_LIBRARIES     (pre_process PRIVATE MPI::MPI_Fortran)
    ENDIF()


    TARGET_COMPILE_OPTIONS(pre_process PRIVATE -DMFC_PRE_PROCESS)


    INSTALL(TARGETS pre_process RUNTIME DESTINATION bin)
ENDIF()


# === === src/simulation_code
IF (MFC_BUILD_SIMULATION)
    HANDLE_FYPP(simulation)

    ADD_EXECUTABLE(simulation "${simulation_srcs}" "${common_srcs}")


    IF (NOT MFC_WITH_OPEN_ACC)
        TARGET_LINK_LIBRARIES(simulation PRIVATE FFTW::FFTW)
    ENDIF()


    IF (MFC_WITH_MPI)
        TARGET_COMPILE_DEFINITIONS(simulation PRIVATE MFC_MPI)
        TARGET_LINK_LIBRARIES     (simulation PRIVATE MPI::MPI_Fortran)
    ENDIF()


    IF (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
        TARGET_LINK_LIBRARIES(simulation PRIVATE CUDA::nvToolsExt)
    ENDIF()


    IF (MFC_WITH_OPEN_ACC)
        TARGET_LINK_LIBRARIES(simulation PRIVATE OpenACC::OpenACC_Fortran)

        IF (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
            TARGET_LINK_LIBRARIES(simulation PRIVATE CUDA::cudart CUDA::cufft cutensor)

            # In the case of the NVHPC SDK, cuTENSOR can be found at a path similar to:
            # /opt/nvidia/hpc_sdk/Linux_x86_64/xx.yy/math_libs/lib64/ on a regular installation.
            # CUDAToolkit_LIBRARY_ROOT would then point to /opt/nvidia/hpc_sdk/Linux_x86_64/xx.yy/cuda/zz.ww.
            # Since the linker doesn't usually consider math_libs/ (for it wasn't told to),
            # we attempt to add it as a link directory.

            SET(_warn_msg_begin "Can't automatically determine the location of cuTENSOR.")
            SET(_warn_msg_end   "If linking fails, please explicitly specify a link directory.")

            IF (EXISTS "${CUDAToolkit_LIBRARY_ROOT}/../../math_libs/lib")
                TARGET_LINK_DIRECTORIES(simulation PRIVATE "${CUDAToolkit_LIBRARY_ROOT}/../../math_libs/lib")

                MESSAGE(WARNING "${_warn_msg_begin} Gave linker ${CUDAToolkit_LIBRARY_ROOT}/../../math_libs/lib   as a potential search path for cuTENSOR (based on a regular NVHPC installation). ${_warn_msg_end} Fortunately, linking will likely succeed.")
            ELSEIF (EXISTS "${CUDAToolkit_LIBRARY_ROOT}/../../math_libs/lib64")
                TARGET_LINK_DIRECTORIES(simulation PRIVATE "${CUDAToolkit_LIBRARY_ROOT}/../../math_libs/lib64")

                MESSAGE(WARNING "${_warn_msg_begin} Gave linker ${CUDAToolkit_LIBRARY_ROOT}/../../math_libs/lib64 as a potential search path for cuTENSOR (based on a regular NVHPC installation). ${_warn_msg_end} Fortunately, linking will likely succeed.")
            ELSE()
                MESSAGE(WARNING "${_warn_msg_begin} No common potential search directories were found. ${_warn_msg_end}")
            ENDIF()
        ENDIF()
    ENDIF()


    TARGET_COMPILE_OPTIONS(simulation PRIVATE -DMFC_SIMULATION)


    INSTALL(TARGETS simulation RUNTIME DESTINATION bin)
ENDIF()


# === === src/post_process_code
IF (MFC_BUILD_POST_PROCESS)
    HANDLE_FYPP(post_process)

    ADD_EXECUTABLE(post_process "${post_process_srcs}" "${common_srcs}")

    TARGET_LINK_LIBRARIES(post_process PRIVATE SILO::SILO HDF5::HDF5 FFTW::FFTW)

    IF (MFC_WITH_MPI)
        TARGET_COMPILE_DEFINITIONS(post_process PRIVATE MFC_MPI)
        TARGET_LINK_LIBRARIES     (post_process PRIVATE MPI::MPI_Fortran)
    ENDIF()

    TARGET_COMPILE_OPTIONS(post_process PRIVATE -DMFC_POST_PROCESS)


    INSTALL(TARGETS post_process RUNTIME DESTINATION bin)
ENDIF()
