
FILE(GLOB pre_process_f90s "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
ADD_EXECUTABLE(pre_process "${pre_process_f90s}")

INSTALL(TARGETS pre_process RUNTIME DESTINATION bin)
