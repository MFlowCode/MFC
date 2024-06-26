# This is a superbuild CMake file
# See: https://cmake.org/cmake/help/latest/module/ExternalProject.html.

cmake_minimum_required(VERSION 3.12)

project(MFC-Dependencies LANGUAGES Fortran)

if (MFC_SILO OR NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    enable_language(C CXX)
endif()

set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

# Imports
include(GNUInstallDirs)
include(ExternalProject)


# Options
option(MFC_FFTW    "Build the FFTW3 dependency" OFF)
option(MFC_HDF5    "Build the HDF5  dependency" OFF)
option(MFC_SILO    "Build the SILO  dependency" OFF)
option(MFC_HIPFORT "Build the HIPFORT dependency" OFF)


# FFTW3
if (MFC_FFTW)
    if (NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        ExternalProject_Add(fftw
            URL        "http://www.fftw.org/fftw-3.3.10.tar.gz"
            CMAKE_ARGS -DBUILD_TESTS=OFF
                       -DBUILD_SHARED_LIBS=OFF
                       "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
        )
    else()
        message(WARNING "The Fortran compiler vendor is Cray so FFTW3 will not be built. We will use cray-fftw instead.")
        add_custom_target(fftw)
    endif()
endif()


# HDF5
if (MFC_HDF5)
    ExternalProject_Add(hdf5
        GIT_REPOSITORY "https://github.com/HDFGroup/hdf5"
        GIT_TAG        hdf5-1_12_2
        GIT_SHALLOW    ON
        GIT_PROGRESS   ON
        CMAKE_ARGS     "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                        -DBUILD_SHARED_LIBS=OFF
                        -DFORTRAN_LIBRARIES=ON
                        -DBUILD_TESTING=OFF
                        -DHDF5_BUILD_UTILS=OFF
                        -DHDF5_BUILD_TOOLS=ON
                        -DHDF5_BUILD_EXAMPLES=OFF
    )
endif()


# SILO
if (MFC_SILO)
    find_package(Git REQUIRED)

    # If we are using the CCE, HDF5 is not built, and we wish to find
    # the system's cray-hdf5. Otherwise, we point SILO to find HDF5 in
    # our common install directory using SILO_HDF5_DIR.
    ExternalProject_Add(silo
        GIT_REPOSITORY "https://github.com/LLNL/Silo"
        GIT_TAG        438477c80d32a3e1757d4584b993f382cace1535
        GIT_PROGRESS   ON
        PATCH_COMMAND  "${GIT_EXECUTABLE}" stash
                       && "${GIT_EXECUTABLE}" apply "${CMAKE_SOURCE_DIR}/Silo.patch"
                       && "${GIT_EXECUTABLE}" apply "${CMAKE_SOURCE_DIR}/Silo-GNU-13.patch"
        CMAKE_ARGS      -DSILO_ENABLE_SHARED=OFF
                        -DSILO_ENABLE_SILOCK=OFF
                        -DSILO_ENABLE_BROWSER=OFF
                        -DHDF5_USE_STATIC_LIBRARIES=ON
                       "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                       "-DSILO_HDF5_DIR=${CMAKE_INSTALL_PREFIX}"
                       "$<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},Cray>:-DCMAKE_MODULE_PATH=${CMAKE_SOURCE_DIR}/../cmake>"
    )

    if (MFC_HDF5)
        add_dependencies(silo hdf5)
    endif()
endif()

# HIPFORT
if (MFC_HIPFORT)
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        ExternalProject_Add(hipfort
            GIT_REPOSITORY "https://github.com/ROCmSoftwarePlatform/hipfort"
            GIT_TAG        rocm-6.0.2
            GIT_SHALLOW    ON
            GIT_PROGRESS   ON
            CMAKE_ARGS     "-DHIPFORT_COMPILER=${CMAKE_Fortran_COMPILER}"
                           "-DHIPFORT_AR=${CMAKE_AR}"
                           "-DHIPFORT_RANLIB=${CMAKE_RANLIB}"
                           "-DHIPFORT_COMPILER_FLAGS=-f free -e F -O0 -h ipa0"
                           "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
        )
    else()
        message(WARNING "The Fortran compiler vendor is not Cray so HIPFORT will not be built.")
        add_custom_target(hipfort)
    endif()
endif()
