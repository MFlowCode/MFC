include(../common/CMakeLists.txt)

FILE(GLOB post_process_f90s "${CMAKE_CURRENT_LIST_DIR}/*.f90")
FILE(GLOB post_process_fpps "${CMAKE_CURRENT_LIST_DIR}/*.fpp")

# Run FYPP
FOREACH(in_fpp_file IN LISTS post_process_fpps)
    STRING(REGEX REPLACE ".fpp\$" ".f90" out_f90_file "${in_fpp_file}")

    GET_FILENAME_COMPONENT(in_fpp_file_name     "${in_fpp_file}"  NAME)
    GET_FILENAME_COMPONENT(out_f90_file_dirpath "${out_f90_file}" DIRECTORY)
    GET_FILENAME_COMPONENT(mytest2_file_name    "${out_f90_file}" NAME)

    SET(out_f90_file "${out_f90_file_dirpath}/autogen/${mytest2_file_name}")

    ADD_CUSTOM_COMMAND(
        OUTPUT          "${out_f90_file}"
        COMMAND         "fypp" "${in_fpp_file}" "${out_f90_file}"
        DEPENDS         "${in_fpp_file}" "${CMAKE_CURRENT_LIST_DIR}/../common/case.fpp"
        COMMENT         "Preprocessing ${in_fpp_file_name}"
        VERBATIM
    )
    
    SET(post_process_f90s ${post_process_f90s} "${out_f90_file}")
ENDFOREACH(in_fpp_file)


ADD_EXECUTABLE(post_process "${common_f90s}" "${post_process_f90s}")


TARGET_INCLUDE_DIRECTORIES(
    post_process
    PRIVATE "${FFTW3_INCLUDE_DIRS}"
            "${SILO_INCLUDE_DIRS}"
            "${HDF5_INCLUDE_DIRS}"
)

TARGET_LINK_DIRECTORIES(
    post_process
    PRIVATE "${FFTW3_LIBRARY_DIRS}"
            "${SILO_LIBRARY_DIRS}"
            "${HDF5_LIBRARY_DIRS}"
)

TARGET_LINK_LIBRARIES(
    post_process
    PRIVATE "${FFTW3_LIBRARIES}"
            "${SILO_LIBRARIES}"
            "${HDF5_LIBRARIES}"
)


INSTALL(TARGETS post_process RUNTIME DESTINATION bin)
TARGET_LINK_LIBRARIES(post_process PRIVATE MPI::MPI_Fortran)
TARGET_COMPILE_OPTIONS(post_process PRIVATE -DMFC_POST_PROCESS)
