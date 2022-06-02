# Run FYPP
FILE(GLOB common_f90s "${CMAKE_CURRENT_LIST_DIR}/*.f90")
FILE(GLOB common_fpps "${CMAKE_CURRENT_LIST_DIR}/*.fpp")

FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/autogen")

FOREACH(in_fpp_file IN LISTS common_fpps)
    STRING(REGEX REPLACE ".fpp\$" ".f90" out_fpp_file "${in_fpp_file}")

    GET_FILENAME_COMPONENT(in_fpp_file_name     "${in_fpp_file}"  NAME)

    GET_FILENAME_COMPONENT(out_fpp_file_dirpath "${out_fpp_file}" DIRECTORY)
    GET_FILENAME_COMPONENT(out_fpp_file_name    "${out_fpp_file}" NAME)

    SET(out_fpp_file "${out_fpp_file_dirpath}/autogen/${out_fpp_file_name}")

    ADD_CUSTOM_COMMAND(
        OUTPUT          "${out_fpp_file}"
        COMMAND         "fypp" "${in_fpp_file}" "${out_fpp_file}"
        DEPENDS         "${in_fpp_file}" "${CMAKE_CURRENT_LIST_DIR}/case.fpp"
        COMMENT         "Preprocessing ${in_fpp_file_name}"
        VERBATIM
    )
    
    SET(common_f90s ${common_f90s} "${out_fpp_file}")
ENDFOREACH(in_fpp_file)

