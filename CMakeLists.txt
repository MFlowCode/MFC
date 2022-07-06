CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

PROJECT(MFC LANGUAGES C CXX Fortran)

# Imports
INCLUDE(GNUInstallDirs)
INCLUDE(ExternalProject)
INCLUDE(CheckIPOSupported)
INCLUDE(CheckFortranCompilerFlag)

# Tell CMake to also check for lib64/
SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

# Add some (required) and cross-compiler compiler flags
SET(CMAKE_Fortran_PREPROCESS ON)

LINK_LIBRARIES(${CMAKE_DL_LIBS})

# Add compiler specific flags
IF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    ADD_COMPILE_OPTIONS(
        -w
        -fdefault-real-8
        -fdefault-double-8
        -ffree-line-length-none
    )

    IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
        ADD_COMPILE_OPTIONS(
            -g
            -fcheck=all
            -fbacktrace
            -fimplicit-none
            -ffpe-trap=invalid,zero,overflow,underflow
        )
    ENDIF()

    IF (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
        ADD_COMPILE_OPTIONS(
            -fallow-invalid-boz
            -fallow-argument-mismatch
        )
    ENDIF()
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    ADD_COMPILE_OPTIONS(
        -Mfreeform
        -Mpreprocess
        -fdefault-real-8
    )
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    ADD_COMPILE_OPTIONS(
        -free
    )
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    ADD_COMPILE_OPTIONS(
        -r8
        -cpp
        -Mfreeform
        -lcutensor
        -Minfo=accel
        -Mr8intrinsics
    )
ENDIF()


IF (CMAKE_BUILD_TYPE STREQUAL "Release")

    CHECK_FORTRAN_COMPILER_FLAG("-march=native" SUPPORTS_MARCH_NATIVE)
    IF (SUPPORTS_MARCH_NATIVE)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
    ELSE()
    	CHECK_FORTRAN_COMPILER_FLAG("-mcpu=native" SUPPORTS_MCPU_NATIVE)
        IF (SUPPORTS_MCPU_NATIVE)
            SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcpu=native")
        ENDIF()
    ENDIF()

    # Check IPO (LTO)
    CHECK_IPO_SUPPORTED(RESULT SUPPORTS_IPO OUTPUT IPO_ERROR)
    IF (SUPPORTS_IPO)
        MESSAGE(STATUS "IPO / LTO enabled")
        SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    ENDIF()

ENDIF()


# Options
OPTION(MFC_WITH_OPEN_ACC      "Build with OpenACC" OFF)
OPTION(MFC_BUILD_PRE_PROCESS  "Build pre_process"  OFF)
OPTION(MFC_BUILD_SIMULATION   "Build simulation"   OFF)
OPTION(MFC_BUILD_POST_PROCESS "Build post_process" OFF)
OPTION(MFC_BUILD_ALL          "Build pre_process, simulation, and post_process" OFF)

IF (MFC_BUILD_ALL)
    SET(MFC_BUILD_PRE_PROCESS  ON FORCE)
    SET(MFC_BUILD_SIMULATION   ON FORCE)
    SET(MFC_BUILD_POST_PROCESS ON FORCE)
ENDIF()

# Find and configure all executables for MPI
FIND_PACKAGE(MPI COMPONENTS Fortran REQUIRED)

IF (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    FIND_PACKAGE(CUDAToolkit REQUIRED)
ENDIF()


# Enable OpenACC
IF (MFC_WITH_OPEN_ACC)
    FIND_PACKAGE(OpenACC REQUIRED)
ENDIF()


# Check compiler vendor is not unsupported
IF (CMAKE_Fortran_COMPILER_ID STREQUAL "AppleClang")
    MESSAGE(FATAL_ERROR "MFC does not support the Apple Clang compilers. Please consult the README for more details.")
ENDIF()


# Check compiler version is not unsupported
IF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    IF (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5)
        MESSAGE(FATAL_ERROR "GNU v5.0 or newer is required to build MFC. [current: ${CMAKE_Fortran_COMPILER_VERSION}]")
    ENDIF()
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    IF (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 21.7)
        MESSAGE(FATAL_ERROR "NVHPC v21.7 or newer is required to build MFC. [current: ${CMAKE_Fortran_COMPILER_VERSION}]")
    ENDIF()
ENDIF()


# Find Dependencies
IF (MFC_BUILD_SIMULATION OR MFC_BUILD_POST_PROCESS)
    FIND_PACKAGE(FFTW3 REQUIRED)
ENDIF()

IF (MFC_BUILD_POST_PROCESS)
    LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/dependencies/")

    FIND_PACKAGE(HDF5 REQUIRED)
    FIND_PACKAGE(SILO REQUIRED)
ENDIF()


# Follow src/ Instructions
ADD_SUBDIRECTORY(src)

