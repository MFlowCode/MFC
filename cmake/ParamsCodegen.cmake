# cmake/ParamsCodegen.cmake — params codegen gen-file lists and build-time custom command

# Generate Fortran parameter namelist/decl includes into the per-target build
# include directories before HANDLE_SOURCES globs them for Fypp.
# cmake_gen.py creates the include dirs itself (path.parent.mkdir), so no
# file(MAKE_DIRECTORY ...) is needed here.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
file(GLOB_RECURSE _mfc_gen_inputs
    "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/params/*.py"
)

# Enumerate the 15 generated .fpp files explicitly so ninja can track them as
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
