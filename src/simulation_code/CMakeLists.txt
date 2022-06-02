
# Run FYPP
FILE(GLOB simulation_f90s "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
FILE(GLOB simulation_fpps "${CMAKE_CURRENT_SOURCE_DIR}/*.fpp")


FOREACH(in_fpp_file IN LISTS simulation_fpps)
    STRING(REGEX REPLACE ".fpp\$" ".f90" out_fpp_file "${in_fpp_file}")
    
    ADD_CUSTOM_COMMAND(
        OUTPUT          "${out_fpp_file}"
        COMMAND          fypp "${in_fpp_file}" "${out_fpp_file}"
        MAIN_DEPENDENCY "${in_fpp_file}"
        VERBATIM
    )
    
    SET(simulation_f90s ${simulation_f90s} "${out_fpp_file}")
ENDFOREACH(in_fpp_file)


ADD_EXECUTABLE(simulation "${simulation_f90s}")


TARGET_INCLUDE_DIRECTORIES(simulation PRIVATE "${FFTW3_INCLUDE_DIRS}")
TARGET_LINK_DIRECTORIES   (simulation PRIVATE "${FFTW3_LIBRARY_DIRS}")
TARGET_LINK_LIBRARIES     (simulation PRIVATE "${FFTW3_LIBRARIES}")


INSTALL(TARGETS simulation RUNTIME DESTINATION bin)
