# Attempt to find LLNL's Silo
# URL:  https://github.com/LLNL/Silo 
# DOCS: https://cmake.org/cmake/help/latest/command/find_library.html
#       https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html

include(FindPackageHandleStandardArgs)


find_path(SILO_INCLUDE_DIR
    NAMES         silo.h
    PATH_SUFFIXES silo
)

find_library(SILO_LIBRARY
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


if (SILO_FOUND AND NOT TARGET SILO::SILO)
    set(SILO_LIBRARIES    "${SILO_LIBRARY}")
    set(SILO_INCLUDE_DIRS "${SILO_INCLUDE_DIR}")

    add_library(SILO::SILO UNKNOWN IMPORTED)

    set_target_properties(SILO::SILO PROPERTIES
        IMPORTED_LOCATION             "${SILO_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${SILO_INCLUDE_DIRS}"
    )
endif()
