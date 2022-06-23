# Initial configuration
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

PROJECT(MFC LANGUAGES C CXX Fortran)


# Imports
INCLUDE(GNUInstallDirs)
INCLUDE(ExternalProject)


# Tell CMake to also check for lib64/
SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

# Find and configure all executables for MPI
FIND_PACKAGE(MPI REQUIRED COMPONENTS Fortran)

ADD_DEFINITIONS(${MPI_Fortran_COMPILE_FLAGS})
INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_DIRS})
LINK_LIBRARIES(${MPI_Fortran_LIBRARIES})
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")


# Add some (required) and cross-compiler compiler flags
SET(CMAKE_Fortran_PREPROCESS ON)


# Add compiler specific flags
IF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    ADD_COMPILE_OPTIONS(
        -w
        -cpp
        -march=native
        -fdefault-real-8
        -fdefault-double-8
        -Wno-argument-mismatch
        -ffree-line-length-none
    )

    IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
        ADD_COMPILE_OPTIONS(
            -fcheck=all
            -fbacktrace
            -fimplicit-none
            -ffpe-trap=invalid,zero,overflow,underflow
        )
    ENDIF()
    
    IF (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
        ADD_COMPILE_OPTIONS(-fallow-invalid-boz)
    ENDIF()
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    ADD_COMPILE_OPTIONS(
        -Mfreeform 
        -Mpreprocess
        -fdefault-real-8
    )
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    ADD_COMPILE_OPTIONS(
        -r8
        -cpp
        -Mfreeform
        -Minfo=accel
        -lnvToolsExt
    )
ENDIF()


# Options
OPTION(MFC_WITH_OPEN_ACC      "Build with OpenACC" OFF)
OPTION(MFC_BUILD_PRE_PROCESS  "Build pre_process"  OFF)
OPTION(MFC_BUILD_SIMULATION   "Build simulation"   OFF)
OPTION(MFC_BUILD_POST_PROCESS "Build post_process" OFF)
OPTION(MFC_BUILD_ALL          "Build pre_process, simulation, and post_process" OFF)

IF(MFC_BUILD_ALL)
    SET(MFC_BUILD_PRE_PROCESS  ON FORCE)
    SET(MFC_BUILD_SIMULATION   ON FORCE)
    SET(MFC_BUILD_POST_PROCESS ON FORCE)
ENDIF()


# Enable OpenACC
IF(MFC_WITH_OPEN_ACC)
    FIND_PACKAGE(OpenACC REQUIRED)
    ADD_COMPILE_OPTIONS(OpenACC_Fortran_FLAGS)

    IF (NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
        MESSAGE(FATAL_ERROR "MFC only supports NVHPC as an OpenACC compiler.")
    ENDIF()
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

