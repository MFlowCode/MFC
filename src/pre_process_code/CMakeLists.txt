include(../common/CMakeLists.txt)

FILE(GLOB pre_process_f90s "${CMAKE_CURRENT_LIST_DIR}/*.f90")
FILE(GLOB pre_process_fpps "${CMAKE_CURRENT_LIST_DIR}/*.fpp")

# Run FYPP
FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/autogen")

FOREACH(in_fpp_file IN LISTS pre_process_fpps)
    STRING(REGEX REPLACE ".fpp\$" ".f90" out_f90_file "${in_fpp_file}")

    GET_FILENAME_COMPONENT(out_f90_file_dirpath "${out_f90_file}" DIRECTORY)
    GET_FILENAME_COMPONENT(mytest2_file_name    "${out_f90_file}" NAME)

    SET(out_f90_file "${out_f90_file_dirpath}/autogen/${mytest2_file_name}")

    ADD_CUSTOM_COMMAND(
        OUTPUT          "${out_f90_file}"
        COMMAND         "fypp" "${in_fpp_file}" "${out_f90_file}"
        MAIN_DEPENDENCY "${in_fpp_file}"
        COMMENT         "Fypp ${in_fpp_file} -> ${out_f90_file}"
        VERBATIM
    )
    
    SET(pre_process_f90s ${pre_process_f90s} "${out_f90_file}")
ENDFOREACH(in_fpp_file)


ADD_EXECUTABLE(pre_process "${common_f90s}" "${pre_process_f90s}")

INSTALL(TARGETS pre_process RUNTIME DESTINATION bin)
TARGET_LINK_LIBRARIES(pre_process PRIVATE MPI::MPI_Fortran)
