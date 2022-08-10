# Attempt to find FFTW
# URL:  https://www.fftw.org/
# DOCS: https://cmake.org/cmake/help/latest/command/find_library.html
#       https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html

INCLUDE(FindPackageHandleStandardArgs)


FIND_PATH(FFTW_INSTALL_PREFIX NAMES include/fftw3.h)

FIND_LIBRARY(
    FFTW_LIBRARIES
    NAMES "fftw3"
    HINTS "${FFTW_INSTALL_PREFIX}/lib"
          "${FFTW_INSTALL_PREFIX}/lib64"
)

SET(FFTW_INCLUDE_DIRS "${FFTW_INSTALL_PREFIX}/include")

FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    FFTW DEFAULT_MSG
    FFTW_LIBRARIES
    FFTW_INCLUDE_DIRS
)
