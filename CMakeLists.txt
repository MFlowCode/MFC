cmake_minimum_required(VERSION 3.18)

project(MFC LANGUAGES C CXX Fortran)


# === Options
option(MFC_WITH_MPI            "Build with MPI"       ON)
option(MFC_WITH_OPEN_ACC       "Build with OpenACC"  OFF)
option(MFC_BUILD_PRE_PROCESS   "Build pre_process"   OFF)
option(MFC_BUILD_SIMULATION    "Build simulation"    OFF)
option(MFC_BUILD_POST_PROCESS  "Build post_process"  OFF)
option(MFC_BUILD_DOC           "Build documentation" OFF)
option(MFC_BUILD_ALL           "Build everything"    OFF)

if (MFC_BUILD_ALL)
    set(MFC_BUILD_PRE_PROCESS  ON FORCE)
    set(MFC_BUILD_SIMULATION   ON FORCE)
    set(MFC_BUILD_POST_PROCESS ON FORCE)
    set(MFC_BUILD_DOC          ON FORCE)
endif()


# === Imports
include(GNUInstallDirs)
include(ExternalProject)
include(CheckIPOSupported)
include(CheckFortranCompilerFlag)


# === Check Compiler Support & Tools
# === === Compiler Support
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5)
        message(FATAL_ERROR "GNU v5.0 or newer is required to build MFC. "
                            "[current: ${CMAKE_Fortran_COMPILER_VERSION}]")
    endif()

elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR
        (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))

    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 21.7)
        message(FATAL_ERROR "NVHPC v21.7 or newer is required to build MFC. "
                            "[current: ${CMAKE_Fortran_COMPILER_VERSION}]")
    endif()

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "AppleClang")

    message(FATAL_ERROR "MFC does not support the Apple Clang compilers. "
                        "Please consult the README for more details.")

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
        -w
        $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
        $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-double-8>
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
            -fcheck=all     -fbacktrace
            -fimplicit-none -ffpe-trap=invalid,zero,overflow,underflow
        )
    endif()

    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-fallow-invalid-boz>
            -fallow-argument-mismatch
        )
    endif()

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")

    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>
        $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8>
    )

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-free>)

elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR
        (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))

    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-r8>
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-Mr8intrinsics>
        $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
        -Minfo=accel
    )

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
        message(STATUS "IPO / LTO enabled")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(STATUS "IPO / LTO not available")
    endif()

endif()


# This macro locates Fypp files in src/<target>/ and instructs CMake how to
# preprocess them to generate their .f90 counterparts, should they be requested
# by an executable or library. They are preprocessed again if they are modified
# or if src/common/case.fpp changes because some #:include it, and produce
# a different .f90 file depending on case.fpp's contents. This macro sets the
# <target>_srcs CMake cache variaible the full list of .f90 files, including the
# preprocessed Fypp files.
macro(HANDLE_FYPP target)
    file(GLOB ${target}_f90s CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/${target}/*.f90")
    file(GLOB ${target}_fpps CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/${target}/*.fpp")

    # Create src/<target>/autogen folder in which to generate Fypp'd .f90 files.
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/${target}/autogen")

    set(${target}_srcs "${${target}_f90s}")

    foreach(fpp_filepath IN LISTS ${target}_fpps)
        string(REGEX REPLACE ".fpp\$" ".f90" f90_filepath "${fpp_filepath}")

        get_filename_component(fpp_filename "${fpp_filepath}" NAME)
        get_filename_component(f90_dirpath  "${f90_filepath}" DIRECTORY)
        get_filename_component(f90_filename "${f90_filepath}" NAME)

        set(f90_filepath "${f90_dirpath}/autogen/${f90_filename}")

        add_custom_command(
            OUTPUT   "${f90_filepath}"
            COMMAND  "${FYPP_EXE}" "${fpp_filepath}" "${f90_filepath}"
            DEPENDS  "${fpp_filepath}" "${CMAKE_CURRENT_SOURCE_DIR}/src/${target}/case.fpp"
            COMMENT  "Preprocessing ${fpp_filename}"
            VERBATIM
        )

        list(APPEND ${target}_srcs "${f90_filepath}")
    endforeach(fpp_filepath)

endmacro()


# === Handle src/ Folder
# === === Fypp / Sources
HANDLE_FYPP(common)
HANDLE_FYPP(pre_process)
HANDLE_FYPP(simulation)
HANDLE_FYPP(post_process)


# === === src/pre_process
if (MFC_BUILD_PRE_PROCESS)
    add_executable(pre_process "${pre_process_srcs}" "${common_srcs}")


    if (MFC_WITH_MPI)
        find_package(MPI COMPONENTS Fortran REQUIRED)

        target_compile_definitions(pre_process PRIVATE MFC_MPI)
        target_link_libraries     (pre_process PRIVATE MPI::MPI_Fortran)
    endif()


    target_compile_options(pre_process PRIVATE -DMFC_PRE_PROCESS)


    install(TARGETS pre_process RUNTIME DESTINATION bin)
endif()


# === === src/simulation
if (MFC_BUILD_SIMULATION)
    add_executable(simulation "${simulation_srcs}" "${common_srcs}")


    if (MFC_WITH_OPEN_ACC)
        find_package(OpenACC REQUIRED)
        target_link_libraries(simulation PRIVATE OpenACC::OpenACC_Fortran)


        if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
            find_package(CUDAToolkit REQUIRED)

            target_link_libraries(simulation
                PRIVATE CUDA::cudart
                        CUDA::cufft
            )


            find_package(cuTENSOR)

            if (NOT cuTENSOR_FOUND)
                message(WARNING
                    "Failed to locate the NVIDIA cuTENSOR library. MFC will be "
                    "built without support for it, disallowing the use of "
                    "cu_tensor=T. This can result in degraded performance.")
            else()
                target_link_libraries(simulation PRIVATE cuTENSOR::cuTENSOR)

                target_compile_definitions(simulation PRIVATE MFC_cuTENSOR)
            endif()


            if     (CMAKE_BUILD_TYPE STREQUAL "Release")
                target_compile_options(simulation
                    PRIVATE -gpu=keep,ptxinfo,lineinfo
                )
            elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
                target_compile_options(simulation
                    PRIVATE -gpu=keep,ptxinfo,lineinfo,autocompare,debug
                )
            endif()

        endif()
    else()
        find_package(FFTW REQUIRED)

        target_link_libraries(simulation PRIVATE FFTW::FFTW)
    endif()


    if (MFC_WITH_MPI)
        find_package(MPI COMPONENTS Fortran REQUIRED)

        target_compile_definitions(simulation PRIVATE MFC_MPI)
        target_link_libraries     (simulation PRIVATE MPI::MPI_Fortran)
    endif()


    if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
        find_package(CUDAToolkit REQUIRED)
        
        target_link_libraries(simulation PRIVATE CUDA::nvToolsExt)
    endif()


    target_compile_options(simulation PRIVATE -DMFC_SIMULATION)


    install(TARGETS simulation RUNTIME DESTINATION bin)
endif()


# === === src/post_process
if (MFC_BUILD_POST_PROCESS)
    add_executable(post_process "${post_process_srcs}" "${common_srcs}")

    find_package(FFTW REQUIRED)
    find_package(SILO REQUIRED)
    find_package(HDF5 REQUIRED)

    target_link_libraries(post_process PRIVATE SILO::SILO HDF5::HDF5 FFTW::FFTW)

    if (MFC_WITH_MPI)
        find_package(MPI COMPONENTS Fortran REQUIRED)

        target_compile_definitions(post_process PRIVATE MFC_MPI)
        target_link_libraries     (post_process PRIVATE MPI::MPI_Fortran)
    endif()

    # -O0 is in response to https://github.com/MFlowCode/MFC-develop/issues/95
    target_compile_options(post_process PRIVATE -DMFC_POST_PROCESS -O0)


    install(TARGETS post_process RUNTIME DESTINATION bin)
endif()


add_custom_target("doc")

macro(GEN_DOCS target name)

    set(DOXYGEN_PROJECT_NAME "${name}")
    set(DOXYGEN_INPUT        "${CMAKE_CURRENT_SOURCE_DIR}/src/${target}/ ${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/")

    file(COPY        "${CMAKE_CURRENT_SOURCE_DIR}/doc/${target}"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doxygen")

    file(COPY        "${CMAKE_CURRENT_SOURCE_DIR}/doc/header.html"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}")

    file(COPY        "${CMAKE_CURRENT_SOURCE_DIR}/doc/footer.html"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}")

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in"
        "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/Doxyfile"
        @ONLY
    )

    file(COPY        "${CMAKE_CURRENT_SOURCE_DIR}/doc/res"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doxygen")

    file(COPY        "${CMAKE_CURRENT_SOURCE_DIR}/doc/config.js"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doxygen")

    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/html/index.html"
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/Doxyfile"
                "${${target}_srcs}"
        COMMAND "${DOXYGEN_EXECUTABLE}"
                "Doxyfile"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}"
        COMMENT "${target}: Generating docs"
    )

    add_custom_target(
        "${target}_doxygen" ALL
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/html/index.html"
    )

    add_dependencies("${target}_doxygen" doxygen-awesome-css)
    add_dependencies("doc" "${target}_doxygen")

    install(
        DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/doxygen/${target}/html/"
        DESTINATION "doc/mfc/${target}"
    )

    install(
        DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/doxygen/res/"
        DESTINATION "doc/mfc/res"
    )

    install(
        FILES       "${CMAKE_CURRENT_SOURCE_DIR}/doc/index.html"
        DESTINATION "doc/mfc"
    )

endmacro()


# === Documentation
if (MFC_BUILD_DOC)

    # === === Find Doxygen
    find_package(Doxygen REQUIRED dot REQUIRED)

    # === === Fetch CSS Theme
    ExternalProject_Add(doxygen-awesome-css
        PREFIX            doxygen-awesome-css
        GIT_REPOSITORY    "https://github.com/jothepro/doxygen-awesome-css"
        GIT_TAG           "a5efba07a3d4fd6317d95657b3095b97e134b791"
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
    GEN_DOCS(landing "MFC")

endif()
