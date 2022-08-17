# Attempt to find FFTW
# URL:  https://www.fftw.org/
# DOCS: https://cmake.org/cmake/help/latest/command/find_library.html
#       https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html

INCLUDE(FindPackageHandleStandardArgs)


FIND_PATH(FFTW_INCLUDE_DIR
    NAMES         fftw3.h
    PATH_SUFFIXES fftw fftw3
)

FIND_LIBRARY(FFTW_LIBRARY
    NAMES         fftw3
    PATH_SUFFIXES fftw fftw3
    NAMES_PER_DIR
)


FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    FFTW
    REQUIRED_VARS
        FFTW_LIBRARY
        FFTW_INCLUDE_DIR
)


IF (FFTW_FOUND AND NOT TARGET FFTW::FFTW)
    SET(FFTW_LIBRARIES    "${FFTW_LIBRARY}")
    SET(FFTW_INCLUDE_DIRS "${FFTW_INCLUDE_DIR}")

    ADD_LIBRARY(FFTW::FFTW UNKNOWN IMPORTED)

    SET_TARGET_PROPERTIES(FFTW::FFTW PROPERTIES
        IMPORTED_LOCATION             "${FFTW_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
    )
ENDIF()
