# This is a superbuild CMake file
# See: https://cmake.org/cmake/help/latest/module/ExternalProject.html.

cmake_minimum_required(VERSION 3.12)

project(MFC-Dependencies LANGUAGES Fortran)

if (MFC_SILO OR NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    enable_language(C CXX)
endif()

# Imports
include(GNUInstallDirs)
include(ExternalProject)


# Options
option(MFC_FFTW "Build the FFTW3 dependency" OFF)
option(MFC_HDF5 "Build the HDF5  dependency" OFF)
option(MFC_SILO "Build the SILO  dependency" OFF)


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
    if (NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
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
                           -DHDF5_BUILD_TOOLS=OFF
                           -DHDF5_BUILD_EXAMPLES=OFF
        )
    else()
        message(WARNING "The Fortran compiler vendor is Cray so HDF5 will not be built. We will use cray-hdf5 instead.")
        add_custom_target(hdf5)
    endif()
endif()


# SILO
if (MFC_SILO)
    ExternalProject_Add(silo
        GIT_REPOSITORY "https://github.com/henryleberre/Silo"
        GIT_TAG        964d993039b467e3d472d3a70afd03380fe2c320
        GIT_SHALLOW    ON
        GIT_PROGRESS   ON
        CMAKE_ARGS     -DSILO_ENABLE_SHARED=OFF
                       -DSILO_ENABLE_SILOCK=OFF
                       -DSILO_ENABLE_BROWSER=OFF
                       -DHDF5_USE_STATIC_LIBRARIES=ON
                       "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                       $<$<BOOL:${MFC_BUILD_HDF5}>:"-DSILO_HDF5_DIR=\"${CMAKE_INSTALL_PREFIX}\"">
    )

    if (MFC_HDF5 AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        add_dependencies(silo hdf5)
    endif()
endif()

