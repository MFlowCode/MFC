# SPDX-License-Identifier: MIT

# Intel ifx/icx CI support
# Developer's note:
#  CMake has extensive documentation available online. Searching "cmake <variable>"
#  or "cmake <function_name>" will return a cmake.org page with more information
#  than you might expect to find or need. Note that any variable or symbol
#  prefixed by CMAKE_ is reserved by CMake and has special meaning.

cmake_minimum_required(VERSION 3.20)


# We include C as a language because - for some reason
# FIND_LIBRARY_USE_LIB64_PATHS is otherwise ignored.

project(MFC LANGUAGES C CXX Fortran)

# Build options exposed to users and their default values.

option(MFC_MPI           "Build with MPI"                                     ON)
option(MFC_OpenACC       "Build with OpenACC"                                OFF)
option(MFC_OpenMP        "Build with OpenMP"                                 OFF)
option(MFC_GCov          "Build with GCov"                                   OFF)
option(MFC_Unified       "Build with unified CPU & GPU memory (GH-200 only)" OFF)
option(MFC_Fastmath      "Build with -gpu=fastmath on NV GPUs"               OFF)
option(MFC_PRE_PROCESS   "Build pre_process"                                 OFF)
option(MFC_SIMULATION    "Build simulation"                                  OFF)
option(MFC_POST_PROCESS  "Build post_process"                                OFF)
option(MFC_SYSCHECK      "Build syscheck"                                    OFF)
option(MFC_DOCUMENTATION "Build documentation"                               OFF)
option(MFC_ALL           "Build everything"                                  OFF)
option(MFC_SINGLE_PRECISION "Build single precision"                         OFF)
option(MFC_MIXED_PRECISION "Build mixed precision"                           OFF)

if (MFC_ALL)
    set(MFC_PRE_PROCESS   ON FORCE)
    set(MFC_SIMULATION    ON FORCE)
    set(MFC_POST_PROCESS  ON FORCE)
    set(MFC_DOCUMENTATION ON FORCE)
endif()

# Validate CMAKE_BUILD_TYPE to catch typos (CMake is case-sensitive).
set(_VALID_BUILD_TYPES "Debug" "Release" "RelDebug" "")
if (NOT CMAKE_BUILD_TYPE IN_LIST _VALID_BUILD_TYPES)
    message(FATAL_ERROR "Unknown CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}'. Valid: Debug, RelDebug, Release")
endif()

# RelDebug: a lighter debug mode for CI. Compiler-specific blocks below add the
# actual flags; these defaults just tell CMake it is a recognised build type.
set(CMAKE_C_FLAGS_RELDEBUG       "-g" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELDEBUG     "-g" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELDEBUG "-g" CACHE STRING "")

if (MFC_SINGLE_PRECISION)
    add_compile_definitions(MFC_SINGLE_PRECISION)
else()
    add_compile_definitions(MFC_DOUBLE_PRECISION)
endif()

if (MFC_MIXED_PRECISION)
    add_compile_definitions(MFC_MIXED_PRECISION)
endif()


# CMake Library Imports

include(GNUInstallDirs)
include(ExternalProject)
include(CheckIPOSupported)
include(CMakeParseArguments)
include(CheckFortranCompilerFlag)


# Check Compiler Support: Some compilers, especially old ones, are known not to
# be able to compile MFC. We throw our own error early - at configure time - in
# such cases.

set(__err_msg "\
CMake detected the ${CMAKE_Fortran_COMPILER_ID} Fortran compiler \
v${CMAKE_Fortran_COMPILER_VERSION}. If you intended to use a different \
compiler (or a different version thereof), please:\n\
 - Install the compiler or load its module. (e.g. module load gcc/10.1)\n\
 - Set/Export the C, CXX, and FC environment variables. (e.g. 'export CC=gcc', \
 'export CXX=g++', and 'export FC=gfortran'.\n\
 - If using mfc.sh, delete the build/<code name> directory and try again. (e.g. 'rm -rf build/pre_process')")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5)
        message(FATAL_ERROR "ERROR: GNU v5.0 or newer is required to build MFC.\n${__err_msg}")
    endif()
    if (MFC_OpenACC)
        message(FATAL_ERROR "ERROR: MFC with GPU processing is not currently compatible with GNU compilers. Please use NVIDIA or Cray compilers.\n${__err_msg}")
    endif()
elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 21.7)
        message(FATAL_ERROR "ERROR: When using NVHPC, v21.7 or newer is required to build MFC.\n${__err_msg}")
    endif()

    if ((CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 24.5) AND (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelDebug") AND MFC_OpenACC)
        message(FATAL_ERROR "ERROR: When using NVHPC, MFC with Debug and GPU options requires NVHPC v24.5 or newer.\n${__err_msg}")
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
    message(FATAL_ERROR "ERROR: MFC does not support the Apple Clang compilers. Please consult the documentation.\n${__err_msg}")
endif()


# Miscellaneous Configuration:
# *  Explicitly link to -ldl (or system equivalent)
# *  Request that FIND_LIBRARY searches <prefix>/lib/ and <prefix>/lib64/
# *  Let FindXXX use custom scripts from toolchain/cmake/.

link_libraries("${CMAKE_DL_LIBS}")
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake/regular")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake/cce")
endif()


# Compiler Flags: Here, we specify our own compiler flags for both release and
# debug builds. These include optimization and debug flags, as well as some that
# are required for a successful build of MFC.

include(cmake/GPU.cmake)

# Fypp discovery and HANDLE_SOURCES macro definition.
include(cmake/Fypp.cmake)

# Params codegen: explicit gen-file lists and the build-time
# add_custom_command/target.  The custom command is the SOLE generator;
# no configure-time execute_process bootstrap.
# NOTE: _mfc_gen_files_* lists must exist before HANDLE_SOURCES is called.
include(cmake/ParamsCodegen.cmake)

HANDLE_SOURCES(pre_process ON)
HANDLE_SOURCES(simulation ON)
HANDLE_SOURCES(post_process ON)
HANDLE_SOURCES(syscheck OFF)

# MFC_SETUP_TARGET function definition and per-target instantiation.
include(cmake/MFCTargets.cmake)

if (MFC_PRE_PROCESS)
    MFC_SETUP_TARGET(TARGET  pre_process
                     SOURCES "${pre_process_SRCs}"
                     MPI)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        target_compile_options(pre_process PRIVATE -hfp0)
        # CCE 19.0.0 IPA ICE: bring_routine_resident SIGSEGV in optcg during
        # ipa_processing. Nondeterministic — any file can become the crash site
        # via cross-file inlining. Safe to disable IPA for the whole target
        # (CPU-only, no GPU device-call requirements). See PR #1286.
        target_compile_options(pre_process PRIVATE -Oipa0)
    endif()
endif()

if (MFC_SIMULATION)
    MFC_SETUP_TARGET(TARGET  simulation
                     SOURCES "${simulation_SRCs}"
                     MPI FFTW OpenACC OpenMP)
    # CCE 19.0.0 optcg ICE: IPA inlining triggers two distinct crashes:
    #   m_bubbles_EL:  castIsValid assertion (InstCombine/foldIntegerTypedPHI)
    #   m_phase_change, m_mpi_common, m_start_up: bring_routine_resident SIGSEGV
    # GPU builds: apply -Oipa0 per-file to known crash sites. Whole-target
    # -Oipa0 is not safe for GPU builds — it breaks cross-file !$acc routine /
    # !DIR$ INLINEALWAYS device-call registration, causing runtime aborts.
    # CPU builds: disable IPA for the whole target (safe — no GPU requirements,
    # and the crash is nondeterministic so any file can become the crash site).
    # See PR #1286.
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        if (MFC_OpenACC OR MFC_OpenMP)
            set_source_files_properties(
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_bubbles_EL.fpp.f90"
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_phase_change.fpp.f90"
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_mpi_common.fpp.f90"
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_start_up.fpp.f90"
                PROPERTIES COMPILE_OPTIONS "-Oipa0"
            )
        else()
            target_compile_options(simulation PRIVATE -Oipa0)
        endif()
    endif()
endif()

if (MFC_POST_PROCESS)
    MFC_SETUP_TARGET(TARGET  post_process
                     SOURCES "${post_process_SRCs}"
                     MPI SILO HDF5 FFTW LAPACK)

    # -O0 is in response to https://github.com/MFlowCode/MFC-develop/issues/95
    target_compile_options(post_process PRIVATE -O0)

    # flang-23/LLD defaults to PIE; SILO and LAPACK static libs on Frontier are
    # non-PIC, producing R_X86_64_32 relocations that LLD rejects in PIE mode.
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
        target_link_options(post_process PRIVATE -no-pie)
    endif()
endif()

if (MFC_SYSCHECK)
    MFC_SETUP_TARGET(TARGET  syscheck
                     SOURCES "${syscheck_SRCs}"
                     MPI OpenACC OpenMP)
endif()

if (MFC_DOCUMENTATION)
    # Files in examples/ are used to generate docs/documentation/examples.md
    file(GLOB_RECURSE examples_DOCs CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/examples/*")

    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/examples.md"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/examples.sh;${examples_DOCs}"
        COMMAND "bash" "${CMAKE_CURRENT_SOURCE_DIR}/docs/examples.sh"
                       "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Generating examples.md"
        VERBATIM
    )

    # Generate case_constraints.md and physics_constraints.md together.
    # Both are produced by gen_constraints.sh, so a single command avoids races.
    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/case_constraints.md"
                "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/physics_constraints.md"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/gen_case_constraints_docs.py"
                "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/gen_physics_docs.py"
                "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/case_validator.py"
                "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/params/definitions.py"
                "${examples_DOCs}"
        COMMAND "bash" "${CMAKE_CURRENT_SOURCE_DIR}/docs/gen_constraints.sh"
                       "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Generating case_constraints.md and physics_constraints.md"
        VERBATIM
    )

    # Generate cli-reference.md from cli/commands.py
    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/cli-reference.md"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/cli/commands.py"
                "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/cli/docs_gen.py"
        COMMAND "bash" "${CMAKE_CURRENT_SOURCE_DIR}/docs/gen_cli_reference.sh"
                       "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Generating cli-reference.md"
        VERBATIM
    )

    # Generate parameters.md from parameter registry
    # docs_gen.py now AST-parses case_validator.py, so it must be a dependency
    file(GLOB_RECURSE params_SRCs CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/params/*.py")
    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/parameters.md"
        DEPENDS "${params_SRCs}"
                "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/mfc/case_validator.py"
        COMMAND "bash" "${CMAKE_CURRENT_SOURCE_DIR}/docs/gen_parameters.sh"
                       "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Generating parameters.md"
        VERBATIM
    )

    file(GLOB common_DOCs CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/*")

    # GEN_DOCS: Given a target name (herein <target>), this macro sets up a
    # target, <target>_docs, that generates documentation for <target> using
    # Doxygen. It is then added as a dependency of the documentation target.
    # Doxygen outputs HTML to the binary (build) directory. It is installed into
    # its own directory: <prefix>/docs/mfc/<target>. We use configure_file to
    # generate the Doxyfile for each target using paths that only apply to that
    # target (or when relative paths are impractical), specified by CMake.

    macro(GEN_DOCS target name)
        set(DOXYGEN_PROJECT_NAME     "\"${name}\"")
        set(DOXYGEN_INPUT            "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/${target}\"")

        foreach (f90 ${${target}_SRCs})
            set(DOXYGEN_INPUT "${DOXYGEN_INPUT} \"${f90}\"")
        endforeach()

        set(DOXYGEN_HTML_HEADER      "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/header.html\"")
        set(DOXYGEN_HTML_FOOTER      "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/footer.html\"")
        set(DOXYGEN_HTML_OUTPUT      "\"${CMAKE_CURRENT_BINARY_DIR}/${target}\"")
        set(DOXYGEN_MATHJAX_CODEFILE "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/config.js\"")
        set(DOXYGEN_PROJECT_LOGO     "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/res/icon.ico\"")
        set(DOXYGEN_CITE_BIB_FILES   "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/references.bib\"")
        set(DOXYGEN_IMAGE_PATH       "\"${CMAKE_CURRENT_SOURCE_DIR}/docs/res\"\
                                      \"${CMAKE_CURRENT_SOURCE_DIR}/docs/${target}\"")

        set(DOXYGEN_WARN_LOGFILE "\"${CMAKE_CURRENT_BINARY_DIR}/${target}-doxygen-warnings.log\"")

        file(MAKE_DIRECTORY "${DOXYGEN_OUTPUT_DIRECTORY}")

        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in"
            "${CMAKE_CURRENT_BINARY_DIR}/${target}-Doxyfile" @ONLY)

        set(opt_example_dependency "")
        set(opt_constraints_dependency "")
        set(opt_physics_dependency "")
        set(opt_cli_reference_dependency "")
        set(opt_parameters_dependency "")
        if (${target} STREQUAL documentation)
            set(opt_example_dependency "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/examples.md")
            set(opt_constraints_dependency "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/case_constraints.md")
            set(opt_physics_dependency "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/physics_constraints.md")
            set(opt_cli_reference_dependency "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/cli-reference.md")
            set(opt_parameters_dependency "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/parameters.md")
        endif()

        file(GLOB_RECURSE ${target}_DOCs CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/${target}/*")
        list(APPEND ${target}_DOCs "${common_DOCs}")

        add_custom_command(
            OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/${target}/html/index.html"
            DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${target}-Doxyfile"
                    "${opt_example_dependency}"
                    "${opt_constraints_dependency}"
                    "${opt_physics_dependency}"
                    "${opt_cli_reference_dependency}"
                    "${opt_parameters_dependency}"
                    "${${target}_SRCs}" "${${target}_DOCs}"
            COMMAND "${DOXYGEN_EXECUTABLE}" "${target}-Doxyfile"
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/docs/postprocess_citations.py"
                    "${CMAKE_CURRENT_BINARY_DIR}/${target}"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
            COMMENT "${target}: Generating documentation"
        )

        add_custom_target(
            "${target}_doxygen" ALL
            DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${target}/html/index.html"
        )

        install(DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/${target}"
                DESTINATION "docs/mfc")

        add_dependencies("${target}_doxygen" doxygen-awesome-css)
        add_dependencies(documentation "${target}_doxygen")
    endmacro()

    add_custom_target(documentation)

    find_package(Doxygen REQUIRED dot REQUIRED)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    # > Fetch CSS Theme
    ExternalProject_Add(doxygen-awesome-css
        PREFIX            doxygen-awesome-css
        GIT_REPOSITORY    "https://github.com/jothepro/doxygen-awesome-css"
        GIT_TAG           "v2.4.1"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     ""
        INSTALL_COMMAND   ""
    )

    set(theme_dirpath
        "${CMAKE_CURRENT_BINARY_DIR}/doxygen-awesome-css/src/doxygen-awesome-css")

    set(DOXYGEN_HTML_EXTRA_STYLESHEET
        "\"${theme_dirpath}/doxygen-awesome.css\"\
         \"${theme_dirpath}/doxygen-awesome-sidebar-only.css\"\
         \"${CMAKE_CURRENT_SOURCE_DIR}/docs/custom.css\"")

    # > Generate Documentation & Landing Page
    GEN_DOCS(pre_process   "MFC")
    GEN_DOCS(simulation    "MFC")
    GEN_DOCS(post_process  "MFC")
    GEN_DOCS(api           "MFC")
    GEN_DOCS(documentation "MFC")

    # Generate API landing pages for pre_process, simulation, post_process.
    # Scans src/{target}/*.fpp to produce module lists in docs/{target}/readme.md.
    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/gen-api-landing.stamp"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/gen_api_landing.py"
                ${pre_process_FPPs} ${pre_process_F90s}
                ${simulation_FPPs} ${simulation_F90s}
                ${post_process_FPPs} ${post_process_F90s}
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/docs/gen_api_landing.py"
                "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMAND "${CMAKE_COMMAND}" -E touch "${CMAKE_CURRENT_BINARY_DIR}/gen-api-landing.stamp"
        COMMENT "Generating API landing pages"
        VERBATIM
    )
    add_custom_target(gen_api_landing DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/gen-api-landing.stamp")
    add_dependencies(pre_process_doxygen gen_api_landing)
    add_dependencies(simulation_doxygen gen_api_landing)
    add_dependencies(post_process_doxygen gen_api_landing)
    add_dependencies(api_doxygen gen_api_landing)

    # Fix @file/@brief headers to match actual module/program declarations.
    # Handles mixed-case Fortran names and catches stale module renames.
    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/fix-file-briefs.stamp"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/fix_file_briefs.py"
                ${pre_process_FPPs} ${pre_process_F90s}
                ${simulation_FPPs} ${simulation_F90s}
                ${post_process_FPPs} ${post_process_F90s}
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/docs/fix_file_briefs.py"
                "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMAND "${CMAKE_COMMAND}" -E touch "${CMAKE_CURRENT_BINARY_DIR}/fix-file-briefs.stamp"
        COMMENT "Fixing @file brief headers"
        VERBATIM
    )
    add_custom_target(fix_file_briefs DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/fix-file-briefs.stamp")
    add_dependencies(pre_process_doxygen fix_file_briefs)
    add_dependencies(simulation_doxygen fix_file_briefs)
    add_dependencies(post_process_doxygen fix_file_briefs)

    # Generate architecture.md from template + module_categories.json + source briefs.
    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/gen-architecture.stamp"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/gen_architecture.py"
                "${CMAKE_CURRENT_SOURCE_DIR}/docs/module_categories.json"
                "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/architecture.md.in"
                ${pre_process_FPPs} ${pre_process_F90s}
                ${simulation_FPPs} ${simulation_F90s}
                ${post_process_FPPs} ${post_process_F90s}
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/docs/gen_architecture.py"
                "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMAND "${CMAKE_COMMAND}" -E touch "${CMAKE_CURRENT_BINARY_DIR}/gen-architecture.stamp"
        COMMENT "Generating architecture page"
        VERBATIM
    )
    add_custom_target(gen_architecture DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/gen-architecture.stamp")
    add_dependencies(documentation_doxygen gen_architecture)

    # Inject per-page last-updated dates into documentation markdown files.
    # Runs after auto-generated .md files exist, before Doxygen processes them.
    # Uses a stamp file so it only runs once per build.
    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/inject-dates.stamp"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/examples.md"
                "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/case_constraints.md"
                "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/physics_constraints.md"
                "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/cli-reference.md"
                "${CMAKE_CURRENT_SOURCE_DIR}/docs/documentation/parameters.md"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/docs/inject-dates.py"
                "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMAND "${CMAKE_COMMAND}" -E touch "${CMAKE_CURRENT_BINARY_DIR}/inject-dates.stamp"
        COMMENT "Injecting page dates into documentation"
        VERBATIM
    )
    add_custom_target(inject_page_dates DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/inject-dates.stamp")
    add_dependencies(documentation_doxygen inject_page_dates)

    # > Copy Resources (main landing page & assets)
    install(DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/docs/res"
            DESTINATION "docs/mfc")

    install(FILES       "${CMAKE_CURRENT_SOURCE_DIR}/docs/robots.txt"
            DESTINATION "docs/mfc")

    install(FILES       "${CMAKE_CURRENT_SOURCE_DIR}/docs/index.html"
                        "${CMAKE_CURRENT_SOURCE_DIR}/docs/404.html"
            DESTINATION "docs/mfc")
endif()

site_name(SITE_NAME)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake/configuration.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/configuration.txt")
