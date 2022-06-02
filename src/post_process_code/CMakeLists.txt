
FILE(GLOB post_process_f90s "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
ADD_EXECUTABLE(post_process "${post_process_f90s}")

# Configure with dependencies
TARGET_INCLUDE_DIRECTORIES(
    post_process
    PRIVATE "${FFTW3_INCLUDE_DIRS}"
            "${HDF5_INCLUDE_DIRS}"
            "${SILO_INCLUDE_DIRS}"
)

TARGET_LINK_LIBRARIES(
    post_process
    PRIVATE "${FFTW3_LIBRARIES}"
            "${HDF5_LIBRARIES}"
            "${SILO_LIBRARIES}"
)

INSTALL(TARGETS post_process RUNTIME DESTINATION bin)
