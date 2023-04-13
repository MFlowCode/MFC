cmake_minimum_required(VERSION 3.20)

# Note: If C is not included as a language,
# FIND_LIBRARY_USE_LIB64_PATHS is ignored.
project(MFC LANGUAGES C Fortran)


# === Options
option(MFC_MPI           "Build with MPI"       ON)
option(MFC_OpenACC       "Build with OpenACC"  OFF)
option(MFC_PRE_PROCESS   "Build pre_process"   OFF)
option(MFC_SIMULATION    "Build simulation"    OFF)
option(MFC_POST_PROCESS  "Build post_process"  OFF)
option(MFC_DOCUMENTATION "Build documentation" OFF)
option(MFC_ALL           "Build everything"    OFF)

if (MFC_BUILD_ALL)
    set(MFC_PRE_PROCESS   ON FORCE)
    set(MFC_SIMULATION    ON FORCE)
    set(MFC_POST_PROCESS  ON FORCE)
    set(MFC_DOCUMENTATION ON FORCE)
endif()


# === Imports
include(GNUInstallDirs)
include(ExternalProject)
include(CheckIPOSupported)
include(CMakeParseArguments)
include(CheckFortranCompilerFlag)


# === Check Compiler Support & Tools
# === === Compiler Support
set(__err_msg "\
CMake detected the ${CMAKE_Fortran_COMPILER_ID} Fortran compiler \
v${CMAKE_Fortran_COMPILER_VERSION}. If you intended to use a different \
compiler (or a different version thereof), please:\n\
 - Install the compiler or load its module. (e.g. module load gcc/10.1)\n\
 - Set/Export the C, CXX, and FC environment variables. (e.g. 'export CC=gcc', \
 'export CXX=g++', and 'export FC=gfortran'.\n\
 - If using mfc.sh, delete the build/<code name> directory and try again. (e.g. 'rm -rf build/pre_process')\n")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5)
        message(FATAL_ERROR "${__err_msg}ERROR: GNU v5.0 or newer is required to build MFC.")
    endif()
elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 21.7)
        message(FATAL_ERROR "${__err_msg}ERROR: NVHPC v21.7 or newer is required to build MFC.")
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
    message(FATAL_ERROR "${__err_msg}ERROR: MFC does not support the Apple Clang compilers. Please consult the documentation.")
endif()

# === === Find Fypp
find_program(FYPP_EXE fypp REQUIRED)


# === Miscellaneous Configuration
# === === Enable C-Preprocessor for Fortran files
set(CMAKE_Fortran_PREPROCESS ON)
# === === Explicitly link to -ldl (or system equivalent)
link_libraries("${CMAKE_DL_LIBS}")
# === === Request that FIND_LIBRARY searches lib/ and lib64/
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)
# === === Let FindXXX use custom scripts from cmake/
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake")


# === Compiler Flags
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
            -fcheck=all     -fbacktrace
            -fimplicit-none
        )
    endif()

    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-fallow-invalid-boz>
            $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>
        )
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    add_compile_options(
        -h nomessage=296:878:1391:1069
        -h msgs -h static -h keepfiles
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-e D -h func_trace)
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>
        $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8>
    )
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-free>)

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -O0 -traceback -debug)
    endif()
elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
        -Minfo=accel
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran:-minline>
        )
    endif()
    
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-C -g -O0 -traceback -Mchkptr -Mchkstk -Minform=inform -Mbounds)
    endif()
endif()


# === Release Optimizations
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    # === === Processor tuning
    CHECK_FORTRAN_COMPILER_FLAG("-march=native" SUPPORTS_MARCH_NATIVE)
    if (SUPPORTS_MARCH_NATIVE)
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-march=native>)
    else()
    	CHECK_FORTRAN_COMPILER_FLAG("-mcpu=native" SUPPORTS_MCPU_NATIVE)
        if (SUPPORTS_MCPU_NATIVE)
            add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-mcpu=native>)
        endif()
    endif()


    # === === Enable LTO/IPO if supported
    CHECK_IPO_SUPPORTED(RESULT SUPPORTS_IPO OUTPUT IPO_ERROR)
    if (SUPPORTS_IPO)
        message(STATUS "Enabled IPO / LTO")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(STATUS "IPO / LTO is NOT available")
    endif()
endif()

## === DEBUG

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(MFC_DEBUG)
endif()

## === HANDLE_SOURCES
# Gather F90 source files for a given target, including common/, and preprocessing .fpp -> .f90.
# Outputs:
#  - <target>_SRCs: List of F90 source filepaths for src/[<target>,common]/ (including .fpp -> .f90)
macro(HANDLE_SOURCES target)
    set(${target}_DIR "${CMAKE_SOURCE_DIR}/src/${target}")
    set(common_DIR    "${CMAKE_SOURCE_DIR}/src/common")

    # Gather src/[<target>,common]/*.f90
    file(GLOB ${target}_F90s CONFIGURE_DEPENDS "${${target}_DIR}/*.f90")
    file(GLOB    common_F90s CONFIGURE_DEPENDS    "${common_DIR}/*.f90")
    set(${target}_SRCs ${${target}_F90s} ${common_F90s})

    # src/[<DIR>,common]/*.fpp -> src/<DIR>/autogen/*.f90
    file(GLOB ${target}_FPPs CONFIGURE_DEPENDS "${${target}_DIR}/*.fpp")
    file(GLOB    common_FPPs CONFIGURE_DEPENDS    "${common_DIR}/*.fpp")

    # Locate src/[<DIR>,common]/include/*.fpp
    file(GLOB ${target}_incs CONFIGURE_DEPENDS "${${target}_DIR}/include/*.fpp")
    file(GLOB    common_incs CONFIGURE_DEPENDS    "${common_DIR}/include/*.fpp")

    file(MAKE_DIRECTORY "${${target}_DIR}/autogen")
    foreach(fpp ${${target}_FPPs} ${common_FPPs})
        cmake_path(GET fpp FILENAME fpp_filename)
        set(f90 "${${target}_DIR}/autogen/${fpp_filename}.f90")

        add_custom_command(
            OUTPUT   ${f90}
            COMMAND  ${FYPP_EXE} -I "${common_DIR}"
                                 -I "${common_DIR}/include"
                                 -I "${${target}_DIR}/include"
                                 --line-numbering
                                 "${fpp}" "${f90}"
            DEPENDS  "${fpp};${${target}_incs};${common_incs}"
            COMMENT  "Preprocessing (Fypp) ${fpp_filename}"
            VERBATIM
        )

        list(APPEND ${target}_SRCs ${f90})
    endforeach()
endmacro()

HANDLE_SOURCES(pre_process)
HANDLE_SOURCES(simulation)
HANDLE_SOURCES(post_process)


## === MFC_SETUP_TARGET
# Setup a target with the following options:
#  - TARGET:  Target name
#  - SOURCES: List of source files
#  - OpenACC: Can be compiled with OpenACC
#  - MPI:     Can be compiled with MPI
function(MFC_SETUP_TARGET)
    cmake_parse_arguments(ARGS "OpenACC;MPI" "TARGET" "SOURCES" ${ARGN})

    add_executable(${ARGS_TARGET} ${ARGS_SOURCES})

    if (MFC_MPI AND ARGS_MPI)
        find_package(MPI COMPONENTS Fortran REQUIRED)

        target_compile_definitions(${ARGS_TARGET} PRIVATE MFC_MPI)
        target_link_libraries     (${ARGS_TARGET} PRIVATE MPI::MPI_Fortran)
    endif()

    if (MFC_OpenACC AND ARGS_OpenACC)
        find_package(OpenACC)

        # This should be equivalent to if (NOT OpenACC_FC_FOUND)
        if (NOT TARGET OpenACC::OpenACC_Fortran)
            message(FATAL_ERROR "OpenACC + Fortran is unsupported.")
        endif()

        target_link_libraries(${ARGS_TARGET} PRIVATE OpenACC::OpenACC_Fortran)

        if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_compile_options(${ARGS_TARGET} PRIVATE
                "-foffload=amdgcn-amdhsa='-march=gfx90a'"
                "-foffload-options=-lgfortran\ -lm"
                "-fno-exceptions")
        elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
            find_package(cuTENSOR)
            if (NOT cuTENSOR_FOUND)
                message(WARNING
                    "Failed to locate the NVIDIA cuTENSOR library. MFC will be "
                    "built without support for it, disallowing the use of "
                    "cu_tensor=T. This can result in degraded performance.")
            else()
                target_link_libraries     (${ARGS_TARGET} PRIVATE cuTENSOR::cuTENSOR)
                target_compile_definitions(${ARGS_TARGET} PRIVATE MFC_cuTENSOR)
            endif()

            if     (CMAKE_BUILD_TYPE STREQUAL "Release")
                target_compile_options(${ARGS_TARGET}
                    PRIVATE -gpu=keep,ptxinfo,lineinfo
                )
            elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
                target_compile_options(${ARGS_TARGET}
                    PRIVATE -gpu=keep,ptxinfo,lineinfo,autocompare,debug -O0
                )
            endif()
        endif()
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        target_compile_options(${ARGS_TARGET} PRIVATE -h noacc -x acc)
    endif()

    string(TOUPPER "${ARGS_TARGET}" ARGS_TARGET_UPPER)

    target_compile_definitions(${ARGS_TARGET} PRIVATE "MFC_${ARGS_TARGET_UPPER}")

    install(TARGETS ${ARGS_TARGET} RUNTIME DESTINATION bin)
endfunction()


if (MFC_PRE_PROCESS)
    MFC_SETUP_TARGET(
        TARGET  pre_process
        SOURCES "${pre_process_SRCs}"
        MPI
    )
endif()


if (MFC_SIMULATION)
    MFC_SETUP_TARGET(
        TARGET simulation
        SOURCES "${simulation_SRCs}"
        MPI OpenACC
    )

    if (MFC_OpenACC AND (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
        find_package(CUDAToolkit REQUIRED)

        target_link_libraries(simulation PRIVATE CUDA::cudart CUDA::cufft)
    else()
        find_package(FFTW REQUIRED)

        target_link_libraries(simulation PRIVATE FFTW::FFTW)
    endif()

    if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
        find_package(CUDAToolkit REQUIRED)
        target_link_libraries(simulation PRIVATE CUDA::nvToolsExt)
    endif()
endif()


if (MFC_POST_PROCESS)
    MFC_SETUP_TARGET(
        TARGET  post_process
        SOURCES "${post_process_SRCs}"
        MPI
    )

    find_package(FFTW REQUIRED)
    find_package(SILO REQUIRED)
    find_package(HDF5 REQUIRED)

    target_link_libraries(post_process PRIVATE SILO::SILO HDF5::HDF5 FFTW::FFTW)

    # -O0 is in response to https://github.com/MFlowCode/MFC-develop/issues/95
    target_compile_options(post_process PRIVATE -O0)
endif()


add_custom_target(documentation)

macro(GEN_DOCS target name)
    set(DOXYGEN_PROJECT_NAME "${name}")
    set(DOXYGEN_INPUT        "${CMAKE_CURRENT_SOURCE_DIR}/src/${target}/ ${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/")

    file(COPY        "${CMAKE_CURRENT_SOURCE_DIR}/docs/${target}"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doxygen")

    file(COPY        "${CMAKE_CURRENT_SOURCE_DIR}/docs/doxygen"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}")

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in"
        "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/Doxyfile"
        @ONLY
    )

    file(COPY        "${CMAKE_CURRENT_SOURCE_DIR}/docs/res"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doxygen")

    file(COPY        "${CMAKE_CURRENT_SOURCE_DIR}/docs/config.js"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doxygen")

    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/html/index.html"
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/Doxyfile"
                "${${target}_SRCs}"
        COMMAND "${DOXYGEN_EXECUTABLE}"
                "Doxyfile"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}"
        COMMENT "${target}: Generating documentation"
    )

    add_custom_target(
        "${target}_doxygen" ALL
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/html/index.html"
    )

    add_dependencies("${target}_doxygen" doxygen-awesome-css)
    add_dependencies(documentation "${target}_doxygen")

    install(
        DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/html/"
        DESTINATION "docs/mfc/${target}"
    )

    install(
        DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/doxygen/res/"
        DESTINATION "docs/mfc/res"
    )

    install(
        FILES       "${CMAKE_CURRENT_SOURCE_DIR}/docs/index.html"
        DESTINATION "docs/mfc"
    )
endmacro()


# === Documentation
if (MFC_DOCUMENTATION)
    # === === Find Doxygen
    find_package(Doxygen REQUIRED dot REQUIRED)

    # === === Fetch CSS Theme
    ExternalProject_Add(doxygen-awesome-css
        PREFIX            doxygen-awesome-css
        GIT_REPOSITORY    "https://github.com/jothepro/doxygen-awesome-css"
        GIT_TAG           "8a34f701d3435309a21080f44fbd4a466979aa19"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     ""
        INSTALL_COMMAND   ""
    )

    set(theme_dirpath "${CMAKE_CURRENT_BINARY_DIR}/doxygen-awesome-css/src/doxygen-awesome-css/")

    set(DOXYGEN_HTML_EXTRA_STYLESHEET "\"${theme_dirpath}/doxygen-awesome.css\"")

    # === === Generate Documentation
    GEN_DOCS(pre_process  "MFC: Pre-Process")
    GEN_DOCS(simulation   "MFC: Simulation")
    GEN_DOCS(post_process "MFC: Post-Process")

    # === === Generate Landing Page
    GEN_DOCS(documentation "MFC")
endif()
