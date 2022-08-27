# Attempt to find NVIDIA cuTENSOR
# URL:  https://developer.nvidia.com/cutensor
# DOCS: https://cmake.org/cmake/help/latest/command/find_library.html
#       https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html

cmake_minimum_required(VERSION 3.18)

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE(CUDAToolkit)

MESSAGE(STATUS "CUDAToolkit_LIBRARY_ROOT=${CUDAToolkit_LIBRARY_ROOT}")

FIND_LIBRARY(cuTENSOR_LIBRARY
    NAMES         cutensor
    HINTS         "${CUDAToolkit_LIBRARY_ROOT}/../../math_libs"
    NAMES_PER_DIR
)


FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    cuTENSOR
    REQUIRED_VARS
        cuTENSOR_LIBRARY
)

IF (cuTENSOR_FOUND AND NOT TARGET cuTENSOR::cuTENSOR)
    SET(cuTENSOR_LIBRARIES    "${cuTENSOR_LIBRARY}")
    SET(cuTENSOR_INCLUDE_DIRS "${cuTENSOR_INCLUDE_DIR}")

    ADD_LIBRARY(cuTENSOR::cuTENSOR UNKNOWN IMPORTED)

    SET_TARGET_PROPERTIES(cuTENSOR::cuTENSOR PROPERTIES
        IMPORTED_LOCATION             "${cuTENSOR_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${cuTENSOR_INCLUDE_DIRS}"
    )
ENDIF()
