# Attempt to find LLNL's Silo
# URL:  https://github.com/LLNL/Silo 
# DOCS: https://cmake.org/cmake/help/latest/command/find_library.html
#       https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html

INCLUDE(FindPackageHandleStandardArgs)


FIND_PATH(SILO_INSTALL_PREFIX NAMES include/silo.h)

FIND_LIBRARY(
    SILO_LIBRARIES
    NAMES "siloh5"
          "silo"
    HINTS "${SILO_INSTALL_PREFIX}/lib"
          "${SILO_INSTALL_PREFIX}/lib64"
)

FIND_PATH(
    SILO_INCLUDE_DIRS
    NAMES silo.h
    HINTS ${SILO_INSTALL_PREFIX}/include
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    SILO DEFAULT_MSG
    SILO_LIBRARIES
    SILO_INCLUDE_DIRS
)

