# Initial configuration
CMAKE_MINIMUM_REQUIRED(VERSION 3.9)

PROJECT(MFC LANGUAGES Fortran)


# Imports
INCLUDE(GNUInstallDirs)
INCLUDE(ExternalProject)

SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install")

# MPI
FIND_PACKAGE(MPI REQUIRED COMPONENTS Fortran)

ADD_DEFINITIONS(${MPI_Fortran_COMPILE_FLAGS})
INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_DIRS})
LINK_DIRECTORIES(${MPI_Fortran_LIBRARIES})
LINK_LIBRARIES(${MPI_Fortran_LIBRARIES})
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")


# Add some (required) and cross-compiler compiler flags
SET(CMAKE_Fortran_PREPROCESS ON)
ADD_COMPILE_OPTIONS("-Wall")

# Add compiler specific flags
IF (CMAKE_FC_COMPILER_ID STREQUAL "NVHPC")
    ADD_COMPILE_OPTIONS("-Minfo=accel")
ENDIF()


# Options
OPTION(MFC_WITH_ACC   "Build with OpenACC" OFF)
OPTION(MFC_BUILD_PRE  "Build pre_process"  OFF)
OPTION(MFC_BUILD_SIM  "Build simulation"   OFF)
OPTION(MFC_BUILD_POST "Build post_process" OFF)
OPTION(MFC_BUILD_ALL  "Build pre_process, simulation, and post_process" OFF)


# Enable OpenACC
IF(MFC_WITH_ACC)
    FIND_PACKAGE(OpenACC REQUIRED)

    IF (CMAKE_FC_COMPILER_ID STREQUAL "NVHPC")
        ADD_COMPILE_OPTIONS(-acc)
    ELSE()
        MESSAGE(FATAL_ERROR "MFC only supports NVHPC as an OpenACC compiler.")
    ENDIF()
ENDIF()


# Check compiler vendor is supported
IF (CMAKE_FC_COMPILER_ID STREQUAL "AppleClang")
    MESSAGE(FATAL_ERROR "MFC does not support the Apple Clang compilers. Please consult the README for more details.")
ENDIF()


# Check compiler version is supported
IF (CMAKE_FC_COMPILER_ID STREQUAL "GNU")
    IF (CMAKE_FC_COMPILER_VERSION VERSION_LESS 5)
        MESSAGE(FATAL_ERROR "GNU v5.0 or newer is required to build MFC. [current=${CMAKE_FC_COMPILER_VERSION}]")
    ENDIF()
ELSEIF (CMAKE_FC_COMPILER_ID STREQUAL "NVHPC")
    IF (CMAKE_FC_COMPILER_VERSION VERSION_LESS 21.7)
        MESSAGE(FATAL_ERROR "NVHPC v21.7 or newer is required to build MFC. [current=${CMAKE_FC_COMPILER_VERSION}]")
    ENDIF()
ENDIF()


IF(MFC_BUILD_ALL)
    SET(MFC_BUILD_PRE  ON FORCE)
    SET(MFC_BUILD_SIM  ON FORCE)
    SET(MFC_BUILD_POST ON FORCE)
ENDIF()


ADD_SUBDIRECTORY(src)

