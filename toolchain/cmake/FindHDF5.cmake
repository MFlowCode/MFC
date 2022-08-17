# Attempt to find HDF5
# URL:  https://www.hdfgroup.org/solutions/hdf5/
# DOCS: https://cmake.org/cmake/help/latest/command/find_library.html
#       https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html

INCLUDE(FindPackageHandleStandardArgs)


FIND_PATH(HDF5_INCLUDE_DIR
    NAMES         hdf5.h
    PATH_SUFFIXES hdf5
)

FIND_LIBRARY(HDF5_LIBRARY
    NAMES         libhdf5 hdf5
    PATH_SUFFIXES hdf5
    NAMES_PER_DIR
)


FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    HDF5
    REQUIRED_VARS
        HDF5_LIBRARY
        HDF5_INCLUDE_DIR
)


IF (HDF5_FOUND AND NOT TARGET HDF5::HDF5)
    SET(HDF5_LIBRARIES    "${HDF5_LIBRARY}")
    SET(HDF5_INCLUDE_DIRS "${HDF5_INCLUDE_DIR}")

    ADD_LIBRARY(HDF5::HDF5 UNKNOWN IMPORTED)

    SET_TARGET_PROPERTIES(HDF5::HDF5 PROPERTIES
        IMPORTED_LOCATION             "${HDF5_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
    )
ENDIF()
