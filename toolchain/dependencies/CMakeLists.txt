# This is a superbuild CMake file
# See: https://cmake.org/cmake/help/latest/module/ExternalProject.html.

CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

PROJECT(MFC-Dependencies)


# Imports
INCLUDE(GNUInstallDirs)
INCLUDE(ExternalProject)


# Options
OPTION(MFC_BUILD_FFTW "Build the FFTW3 dependency" OFF)
OPTION(MFC_BUILD_HDF5 "Build the HDF5  dependency" OFF)
OPTION(MFC_BUILD_SILO "Build the SILO  dependency" OFF)


# FFTW3
IF (MFC_BUILD_FFTW)
    ExternalProject_Add(fftw
        URL        "http://www.fftw.org/fftw-3.3.10.tar.gz"
        CMAKE_ARGS -DBUILD_TESTS=OFF
                   -DBUILD_SHARED_LIBS=OFF
                   "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
    )
ENDIF()


# HDF5
FIND_PACKAGE(HDF5)
IF (MFC_BUILD_HDF5 OR (MFC_BUILD_SILO AND NOT HDF5_FOUND))
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
        GIT_TAG          964d993039b467e3d472d3a70afd03380fe2c320
        CMAKE_ARGS       -DSILO_ENABLE_SHARED=OFF
                         -DSILO_ENABLE_SILOCK=OFF
                         -DSILO_ENABLE_BROWSER=OFF
                         -DHDF5_USE_STATIC_LIBRARIES=ON     
                         "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}" 
                         $<$<BOOL:${MFC_BUILD_HDF5}>:"-DSILO_HDF5_DIR=\"${CMAKE_INSTALL_PREFIX}\"">
    )

    IF (MFC_BUILD_HDF5)
        ADD_DEPENDENCIES(silo hdf5)
    ENDIF()
ENDIF()
