# This is a superbuild CMake file
# See: https://cmake.org/cmake/help/latest/module/ExternalProject.html.

CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

PROJECT(MFC-Dependencies)


# Imports
INCLUDE(GNUInstallDirs)
INCLUDE(ExternalProject)


# Options
OPTION(MFC_BUILD_FFTW3 "Build the FFTW3 dependency" OFF)
OPTION(MFC_BUILD_HDF5  "Build the HDF5  dependency" OFF)
OPTION(MFC_BUILD_SILO  "Build the SILO  dependency (implies MFC_BUILD_HDF5)" OFF)

IF (MFC_BUILD_SILO)
    SET(MFC_BUILD_HDF5 ON)
ENDIF()


# FFTW3
IF (MFC_BUILD_FFTW3)
    ExternalProject_Add(fftw3
        URL        "http://www.fftw.org/fftw-3.3.10.tar.gz"
        CMAKE_ARGS -DBUILD_TESTS=OFF
                   -DBUILD_SHARED_LIBS=OFF
                   "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
    )
ENDIF()


# HDF5
IF (MFC_BUILD_HDF5)
    ExternalProject_Add(hdf5
        GIT_REPOSITORY "https://github.com/HDFGroup/hdf5"
        GIT_TAG        hdf5-1_12_2
        CMAKE_ARGS     "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                       -DBUILD_SHARED_LIBS=OFF
                       -DFORTRAN_LIBRARIES=ON
                       -DBUILD_TESTING=OFF
                       -DHDF5_BUILD_UTILS=OFF
                       -DHDF5_BUILD_TOOLS=OFF
                       -DHDF5_BUILD_EXAMPLES=OFF
    )
ENDIF()


# SILO
IF (MFC_BUILD_SILO)
    ExternalProject_Add(silo
        GIT_REPOSITORY   "https://github.com/henryleberre/Silo"
        GIT_TAG          bd507bf1d739902c437c8b68807c423878cbace0
        CMAKE_ARGS       "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                         "-DSILO_HDF5_DIR=${CMAKE_INSTALL_PREFIX}"
                         -DSILO_ENABLE_SHARED=OFF
                         -DSILO_ENABLE_SILOCK=OFF
                         -DSILO_ENABLE_BROWSER=OFF                         
    )

    ADD_DEPENDENCIES(silo hdf5)
ENDIF()

