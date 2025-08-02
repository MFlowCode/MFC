# This is a superbuild CMake file
# See: https://cmake.org/cmake/help/latest/module/ExternalProject.html.

cmake_minimum_required(VERSION 3.12)

project(MFC-Dependencies LANGUAGES C CXX Fortran)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/regular")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/cce")
endif()

set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

# Imports
include(GNUInstallDirs)
include(ExternalProject)


# Options
option(MFC_FFTW    "Build the FFTW3 dependency" OFF)
option(MFC_HDF5    "Build the HDF5  dependency" OFF)
option(MFC_SILO    "Build the SILO  dependency" OFF)
option(MFC_LAPACK  "Build the LAPACK  dependency" OFF)
option(MFC_HIPFORT "Build the HIPFORT dependency" OFF)


# FFTW3
if (MFC_FFTW)
    find_package(FFTW QUIET)
    if (FFTW_FOUND)
        message(STATUS "FFTW found.")
        add_custom_target(fftw)
    else()
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
endif()


# HDF5
if (MFC_HDF5)
    find_package(HDF5)
    if (HDF5_FOUND)
        message(STATUS "HDF5 found.")
        add_custom_target(hdf5)
    else()
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
endif()


# SILO
if (MFC_SILO)
    find_package(SILO)
    if (SILO_FOUND)
        message(STATUS "SILO found.")
        add_custom_target(silo)
    else()
        find_package(Git REQUIRED)

        ExternalProject_Add(silo
            GIT_REPOSITORY "https://github.com/LLNL/Silo"
            GIT_TAG        9af504ef4fb79153e1fbf3bdb75421b6b65f6dc4
            GIT_PROGRESS   ON
            PATCH_COMMAND  "${GIT_EXECUTABLE}" stash
                        && "${GIT_EXECUTABLE}" apply "${CMAKE_SOURCE_DIR}/Silo.patch"
            CMAKE_ARGS     -DSILO_ENABLE_SHARED=OFF
                           -DSILO_ENABLE_SILOCK=OFF
                           -DSILO_ENABLE_BROWSER=OFF
                           -DFIND_LIBRARY_USE_LIB64_PATHS=ON
                           "-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"
                           "-DCMAKE_FIND_ROOT_PATH=${CMAKE_FIND_ROOT_PATH}"
                           "-DCMAKE_FIND_PACKAGE_REDIRECTS_DIR=${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}"
                           "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                           "$<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},Cray>:-DCMAKE_MODULE_PATH=${CMAKE_SOURCE_DIR}/../cmake/cce>"
        )
    endif()
endif()

# LAPACK
 if (MFC_LAPACK)
     find_package(LAPACK QUIET)
     if (LAPACK_FOUND)
         message(STATUS "LAPACK found.")
         add_custom_target(lapack)
     else()
         # Use reference LAPACK from netlib
         find_package(Git REQUIRED)

         ExternalProject_Add(lapack
             GIT_REPOSITORY "https://github.com/Reference-LAPACK/lapack"
             GIT_TAG        v3.12.0
             GIT_SHALLOW    ON
             GIT_PROGRESS   ON
             CMAKE_ARGS     -DBUILD_SHARED_LIBS=OFF
                            -DBUILD_TESTING=OFF
                            -DCBLAS=OFF
                            -DLAPACKE=OFF
                            -DBUILD_DEPRECATED=OFF
                            "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                            "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
         )
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
