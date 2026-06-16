# cmake/Fypp.cmake — FYPP_EXE discovery and HANDLE_SOURCES macro

# Fypp is required to build MFC. We try to locate it. The path to the binary is
# returned in FYPP_EXE upon success. This is used later.

find_program(FYPP_EXE fypp REQUIRED)


# HANDLE_SOURCES: Given a target (herein <target>):
#
# *  Locate all source files for <target> of the type
#
#               src/[<target>,common]/[.,include]/*.[f90,fpp].
#
# *  For each .fpp file found with filepath <dirpath>/<filename>.fpp, using a
#    custom command, instruct CMake how to generate a file with path
#
#                     src/<target>/fypp/<filename>.f90
#
#    by running Fypp on <dirpath>/<filename>.fpp. It is important to understand
#    that this does not actually run the pre-processor. Rather, it instructs
#    CMake what to do when it finds a src/<target>/fypp/<filename>.f90 path
#    in the source list for a target. Thus, an association is made from an .f90
#    file to its corresponding .fpp file (if applicable) even though the
#    generation is of the form .fpp -> .f90.
#
#    This design has one limitation: If an .fpp file depends on another, for
#    example if it '#:include's it and uses a macro defined in it, then the
#    dependency will not be tracked. A modification to the .fpp file it depends
#    on will not trigger a re-run of Fypp on the .fpp file that depends on it.
#    As a compromise, both in speed and complexity, all .f90 files generated
#    from .fpp files are re-generated not only when their corresponding .fpp
#    file is modified, but also when any file with filepath of the form
#
#                     src/[<target>,common]/include/*.fpp
#
#    is modified. This is a reasonable compromise as modifications to .fpp files
#    in the include directories will be rare - by design. Other approaches would
#    have required a more complex CMakeLists.txt file (perhaps parsing the .fpp
#    files to determine their dependencies) or the endurment of longer
#    compilation times (by way of re-running Fypp on all .fpp files every time
#    one of them is modified).
#
#    .fpp files in src/common are treated as if they were in src/<target> (not
#    pre-processed to src/common/fypp/) so as not to clash with other targets'
#    .fpp files (this has caused problems in the past).
#
# *  Export, in the variable <target>_SRCs, a list of all source files (.f90)
#    that would compile to produce <target>. If <target> includes .fpp files,
#    then the list will include the paths to the corresponding .f90 files that
#    Fypp would generate from the .fpp files.
#
# This design allows us to be flexible in our use of Fypp as we don't have to
# worry about running the pre-processor on .fpp files when we create executables
# and generate documentation. Instead, we can simply include the list of .f90
# files that will eventually be used to compile <target>.

macro(HANDLE_SOURCES target useCommon)

    set(${target}_DIR "${CMAKE_SOURCE_DIR}/src/${target}")
    set(common_DIR    "${CMAKE_SOURCE_DIR}/src/common")

    string(TOUPPER ${target} ${target}_UPPER)

    # Gather:
    # *          src/[<target>,(common)]/*.f90
    # * (if any) <build>/modules/<target>/*.f90
    file(GLOB ${target}_F90s CONFIGURE_DEPENDS "${${target}_DIR}/*.f90"
                                               "${CMAKE_BINARY_DIR}/modules/${target}/*.f90")
    set(${target}_SRCs ${${target}_F90s})
    if (${useCommon})
        file(GLOB common_F90s CONFIGURE_DEPENDS "${common_DIR}/*.f90")
        list(APPEND ${target}_SRCs ${common_F90s})
    endif()

    # Gather:
    # *          src/[<target>,(common)]/*.fpp]
    # * (if any) <build>/modules/<target>/*.fpp
    file(GLOB ${target}_FPPs CONFIGURE_DEPENDS "${${target}_DIR}/*.fpp"
                                               "${CMAKE_BINARY_DIR}/modules/${target}/*.fpp")
    if (${useCommon})
        file(GLOB common_FPPs CONFIGURE_DEPENDS "${common_DIR}/*.fpp")

        # If we're building post_process, exclude m_compute_levelset.fpp
        if("${target}" STREQUAL "post_process")
            list(FILTER common_FPPs EXCLUDE REGEX ".*/m_compute_levelset\.fpp$")
            list(FILTER common_FPPs EXCLUDE REGEX ".*/m_ib_patches\.fpp$")
        endif()

        list(APPEND ${target}_FPPs ${common_FPPs})
    endif()

    # Gather:
    # *          src/[<target>,common]/include/*.fpp
    # * generated includes from build/include/<target>/ (explicit list, not a GLOB)
    file(GLOB ${target}_incs CONFIGURE_DEPENDS "${${target}_DIR}/include/*.fpp")
    list(APPEND ${target}_incs ${_mfc_gen_files_${target}})

    if (${useCommon})
        file(GLOB common_incs CONFIGURE_DEPENDS "${common_DIR}/include/*.fpp")
        list(APPEND ${target}_incs ${common_incs})
    endif()

    # /path/to/*.fpp (used by <target>) -> <build>/fypp/<target>/*.f90
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/fypp/${target}")
    foreach(fpp ${${target}_FPPs})
        cmake_path(GET fpp FILENAME fpp_filename)
        set(f90 "${CMAKE_BINARY_DIR}/fypp/${target}/${fpp_filename}.f90")

        add_custom_command(
            OUTPUT   ${f90}
            COMMAND  ${FYPP_EXE} -m re
                                 -I "${CMAKE_BINARY_DIR}/include/${target}"
                                 -I "${${target}_DIR}/include"
                                 -I "${common_DIR}/include"
                                 -I "${common_DIR}"
                                 -D MFC_${CMAKE_Fortran_COMPILER_ID}
                                 -D MFC_${${target}_UPPER}
                                 -D MFC_COMPILER="${CMAKE_Fortran_COMPILER_ID}"
				                 -D MFC_CASE_OPTIMIZATION=False
                                 -D chemistry=False
                                 --line-numbering
                                 --no-folding
								 --line-length=999
		 						 --line-numbering-mode=nocontlines
                                 ${FYPP_GCOV_OPTS}
                                 "${fpp}" "${f90}"
            DEPENDS  "${fpp};${${target}_incs}"
            COMMENT  "Preprocessing (Fypp) ${fpp_filename}"
            VERBATIM
        )

        list(APPEND ${target}_SRCs ${f90})
    endforeach()
endmacro()
