# Attempt to find LAPACK (Linear Algebra PACKage)
 # URL:  https://www.netlib.org/lapack/
 # DOCS: https://cmake.org/cmake/help/latest/command/find_library.html
 #       https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html

 include(FindPackageHandleStandardArgs)

 # Special handling for Cray systems which have optimized math libraries
 if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
     # On Cray systems, LAPACK is typically provided by the cray-libsci package
     find_library(LAPACK_LIBRARY
         NAMES sci_cray sci_gnu sci_intel sci_pgi sci
         NAMES_PER_DIR
     )
     set(BLAS_LIBRARY "")  # BLAS is included in the sci library
 else()
     # Find LAPACK library for other compilers
     find_library(LAPACK_LIBRARY
         NAMES         lapack
         PATH_SUFFIXES lapack
         NAMES_PER_DIR
     )

     # Find BLAS library (required by LAPACK)
     find_library(BLAS_LIBRARY
         NAMES         blas openblas
         PATH_SUFFIXES blas
         NAMES_PER_DIR
     )

     # Some LAPACK implementations include BLAS
     if (NOT BLAS_LIBRARY)
         set(BLAS_LIBRARY "")
     endif()
 endif()

 FIND_PACKAGE_HANDLE_STANDARD_ARGS(
     LAPACK
     REQUIRED_VARS
         LAPACK_LIBRARY
 )

 if (LAPACK_FOUND AND NOT TARGET LAPACK::LAPACK)
     set(LAPACK_LIBRARIES "${LAPACK_LIBRARY}")
     if (BLAS_LIBRARY)
         list(APPEND LAPACK_LIBRARIES "${BLAS_LIBRARY}")
     endif()

     add_library(LAPACK::LAPACK INTERFACE IMPORTED)

     set_target_properties(LAPACK::LAPACK PROPERTIES
         INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}"
     )

     # Add math library for linking (commonly needed)
     if (UNIX AND NOT APPLE)
         set_property(TARGET LAPACK::LAPACK APPEND PROPERTY
             INTERFACE_LINK_LIBRARIES "m")
     endif()
endif()