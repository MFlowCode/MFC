# cmake/ParamsCodegen.cmake — params codegen stamp, gen-file lists, and build-time custom command

# Generate Fortran parameter namelist/decl includes into the per-target build
# include directories before HANDLE_SOURCES globs them for Fypp.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(_mfc_gen_stamp "${CMAKE_BINARY_DIR}/mfc_params_gen.stamp")
file(GLOB_RECURSE _mfc_gen_inputs
    "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/params/*.py"
)
set(_mfc_needs_regen FALSE)
if(NOT EXISTS "${_mfc_gen_stamp}")
    set(_mfc_needs_regen TRUE)
else()
    foreach(_input IN LISTS _mfc_gen_inputs)
        if("${_input}" IS_NEWER_THAN "${_mfc_gen_stamp}")
            set(_mfc_needs_regen TRUE)
            break()
        endif()
    endforeach()
endif()
if(_mfc_needs_regen)
    execute_process(
        COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/params/generators/cmake_gen.py"
                "${CMAKE_BINARY_DIR}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        RESULT_VARIABLE _mfc_gen_result
        ERROR_VARIABLE  _mfc_gen_error
    )
    if(NOT _mfc_gen_result EQUAL 0)
        message(FATAL_ERROR "Fortran param generation failed:\n${_mfc_gen_error}")
    endif()
    file(TOUCH "${_mfc_gen_stamp}")
endif()

# Enumerate the 10 generated .fpp files explicitly so ninja can track them as
# build-time outputs and so HANDLE_SOURCES does not need a configure-time GLOB
# of ${CMAKE_BINARY_DIR}/include/<target>/ (which fails when the dir is empty).
set(_mfc_gen_inc "${CMAKE_BINARY_DIR}/include")
set(_mfc_gen_files_pre_process
    "${_mfc_gen_inc}/pre_process/generated_namelist.fpp"
    "${_mfc_gen_inc}/pre_process/generated_decls.fpp"
    "${_mfc_gen_inc}/pre_process/generated_constants.fpp"
    "${_mfc_gen_inc}/pre_process/generated_bcast.fpp"
    "${_mfc_gen_inc}/pre_process/generated_case_opt_decls.fpp"
)
set(_mfc_gen_files_simulation
    "${_mfc_gen_inc}/simulation/generated_namelist.fpp"
    "${_mfc_gen_inc}/simulation/generated_decls.fpp"
    "${_mfc_gen_inc}/simulation/generated_constants.fpp"
    "${_mfc_gen_inc}/simulation/generated_bcast.fpp"
    "${_mfc_gen_inc}/simulation/generated_case_opt_decls.fpp"
)
set(_mfc_gen_files_post_process
    "${_mfc_gen_inc}/post_process/generated_namelist.fpp"
    "${_mfc_gen_inc}/post_process/generated_decls.fpp"
    "${_mfc_gen_inc}/post_process/generated_constants.fpp"
    "${_mfc_gen_inc}/post_process/generated_bcast.fpp"
    "${_mfc_gen_inc}/post_process/generated_case_opt_decls.fpp"
)
set(_mfc_gen_files_syscheck)
set(_mfc_all_gen_files
    ${_mfc_gen_files_pre_process}
    ${_mfc_gen_files_simulation}
    ${_mfc_gen_files_post_process}
)

# Build-time regeneration: ninja re-runs cmake_gen.py and re-preprocesses
# any .fpp that includes a generated file whenever a params .py changes.
add_custom_command(
    OUTPUT   ${_mfc_all_gen_files}
    COMMAND  "${Python3_EXECUTABLE}"
             "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/params/generators/cmake_gen.py"
             "${CMAKE_BINARY_DIR}"
    DEPENDS  ${_mfc_gen_inputs}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT  "Regenerating Fortran parameter includes"
    VERBATIM
)
add_custom_target(mfc_params_gen DEPENDS ${_mfc_all_gen_files})
