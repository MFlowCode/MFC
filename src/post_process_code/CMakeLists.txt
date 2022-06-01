
FILE(GLOB post_process_f90s "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
ADD_EXECUTABLE(post_process "${post_process_f90s}")

TARGET_INCLUDE_DIRECTORIES(post_process PRIVATE "${CMAKE_INSTALL_PREFIX}/include")
TARGET_LINK_LIBRARIES     (post_process PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/libfftw3.a")
TARGET_LINK_LIBRARIES     (post_process PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/libsiloh5.so")

INSTALL(TARGETS post_process RUNTIME DESTINATION bin)
