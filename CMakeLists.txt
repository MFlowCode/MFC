cmake_minimum_required(VERSION 3.16)
project(MFC LANGUAGES Fortran C CXX)

enable_language(Fortran)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Find MPI
find_package(MPI REQUIRED)

# Add subdirectories
add_subdirectory(src)

# Unit tests
enable_testing()
add_subdirectory(tests)

# Add examples
add_subdirectory(examples)

# Install configuration
install(DIRECTORY src DESTINATION .)

# --- Removed dependency_check executable ---
# Previously we had:
# find_package(Boost REQUIRED)
# add_executable(dependency_check examples/dependency_check/dependency_check.cpp)
# target_include_directories(dependency_check PRIVATE ${Boost_INCLUDE_DIRS})
# These were removed since dependency verification is now handled
# directly inside scripts/install_dependencies.sh.