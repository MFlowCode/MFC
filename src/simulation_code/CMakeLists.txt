include(../common/CMakeLists.txt)

FILE(GLOB simulation_f90s "${CMAKE_CURRENT_LIST_DIR}/*.f90")
FILE(GLOB simulation_fpps "${CMAKE_CURRENT_LIST_DIR}/*.fpp")

# Run FYPP
FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/autogen")

FOREACH(in_fpp_file IN LISTS simulation_fpps)
    STRING(REGEX REPLACE ".fpp\$" ".f90" out_f90_file "${in_fpp_file}")

    GET_FILENAME_COMPONENT(in_fpp_file_name     "${in_fpp_file}"  NAME)
    GET_FILENAME_COMPONENT(out_f90_file_dirpath "${out_f90_file}" DIRECTORY)
    GET_FILENAME_COMPONENT(mytest2_file_name    "${out_f90_file}" NAME)

    SET(out_f90_file "${out_f90_file_dirpath}/autogen/${mytest2_file_name}")

    ADD_CUSTOM_COMMAND(
        OUTPUT   "${out_f90_file}"
        COMMAND  "fypp" "${in_fpp_file}" "${out_f90_file}"
        DEPENDS  "${in_fpp_file}" "${CMAKE_CURRENT_LIST_DIR}/../common/case.fpp"
        COMMENT  "Preprocessing ${in_fpp_file_name}"
        VERBATIM
    )

    SET(simulation_f90s ${simulation_f90s} "${out_f90_file}")
ENDFOREACH(in_fpp_file)


ADD_EXECUTABLE(simulation "${common_f90s}" "${simulation_f90s}")


TARGET_INCLUDE_DIRECTORIES(simulation PRIVATE "${FFTW3_INCLUDE_DIRS}")
TARGET_LINK_DIRECTORIES   (simulation PRIVATE "${FFTW3_LIBRARY_DIRS}")
TARGET_LINK_LIBRARIES     (simulation PRIVATE "${FFTW3_LIBRARIES}")


INSTALL(TARGETS simulation RUNTIME DESTINATION bin)
TARGET_LINK_LIBRARIES(simulation PRIVATE MPI::MPI_Fortran)


IF (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    TARGET_LINK_LIBRARIES(simulation PRIVATE CUDA::nvToolsExt)
ENDIF()


IF (MFC_WITH_OPEN_ACC)
    TARGET_LINK_LIBRARIES(simulation PRIVATE OpenACC::OpenACC_Fortran)

    IF (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
        TARGET_LINK_LIBRARIES(simulation PRIVATE CUDA::cufft)
    ENDIF()
ENDIF ()

TARGET_COMPILE_OPTIONS(simulation PRIVATE -DMFC_SIMULATION)
