# Attempt to find LLNL's Silo
# URL:  https://github.com/LLNL/Silo 
# DOCS: https://cmake.org/cmake/help/latest/command/find_library.html
#       https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html

INCLUDE(FindPackageHandleStandardArgs)


FIND_PATH(SILO_INCLUDE_DIR
    NAMES         silo.h
    PATH_SUFFIXES silo
)

FIND_LIBRARY(SILO_LIBRARY
    NAMES         siloh5 silo
    PATH_SUFFIXES silo
    NAMES_PER_DIR
)


FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    SILO
    REQUIRED_VARS
        SILO_LIBRARY
        SILO_INCLUDE_DIR
)


IF (SILO_FOUND AND NOT TARGET SILO::SILO)
    SET(SILO_LIBRARIES    "${SILO_LIBRARY}")
    SET(SILO_INCLUDE_DIRS "${SILO_INCLUDE_DIR}")

    ADD_LIBRARY(SILO::SILO UNKNOWN IMPORTED)

    SET_TARGET_PROPERTIES(SILO::SILO PROPERTIES
        IMPORTED_LOCATION             "${SILO_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${SILO_INCLUDE_DIRS}"
    )
ENDIF()
