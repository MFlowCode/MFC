# cmake/MFCTargets.cmake — MFC_SETUP_TARGET macro and related

# MFC_SETUP_TARGET: Given a target (herein <target>), this macro creates a new
# executable <target> with the appropriate sources, compiler definitions, and
# linked libraries (assuming HANDLE_SOURCES was called on <target>).
#
# Inputs:
#  * TARGET   Target name (<target>)
#  * SOURCES  List of source files (.f90)
#  * OpenACC  (optional) Can be compiled with OpenACC.
#  * MPI      (optional) Can be compiled with MPI.
#  * SILO     (optional) Should be linked with SILO.
#  * HDF5     (optional) Should be linked with HDF5.
#  * FFTW     (optional) Should be linked with an FFTW-like library (fftw/cufftw),
#             depending on whether OpenACC is enabled and which compiler is
#             being used.
#  * LAPACK   (optional) Should be linked with LAPACK

function(MFC_SETUP_TARGET)
    cmake_parse_arguments(ARGS "OpenACC;MPI;SILO;HDF5;FFTW;LAPACK;OpenMP" "TARGET" "SOURCES" ${ARGN})

    add_executable(${ARGS_TARGET} ${ARGS_SOURCES})
    set(IPO_TARGETS ${ARGS_TARGET})
    # Here we need to split into "library" and "executable" to perform IPO on the NVIDIA compiler.
    # A little hacky, but it *is* an edge-case for *one* compiler.
    if (NVHPC_USE_TWO_PASS_IPO AND NOT(MFC_OpenMP AND ARGS_OpenMP))
        # nvfortran -Mextract does not produce .o files, only inline library
        # data. An OBJECT library with -Mextract causes CMake to rebuild
        # everything on every build because the expected .o outputs never
        # exist. We use a wrapper script as RULE_LAUNCH_COMPILE that runs
        # the compiler and then touches the expected .o output file.
        set(_ipo_wrapper "${CMAKE_BINARY_DIR}/${ARGS_TARGET}_extract_wrapper.sh")
        file(WRITE "${_ipo_wrapper}" [=[#!/bin/sh
# Find the -o argument (the object file CMake expects)
out=
prev=
for arg do
  if [ "$prev" = "-o" ]; then out="$arg"; break; fi
  prev="$arg"
done
# Run the compiler; propagate its exit status on failure
"$@"
status=$?
[ "$status" -eq 0 ] || exit "$status"
# Touch the .o so CMake's dependency tracking sees it
[ -n "$out" ] && touch "$out"
exit 0
]=])
        file(CHMOD "${_ipo_wrapper}" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
        add_library(${ARGS_TARGET}_lib OBJECT ${ARGS_SOURCES})
        set_target_properties(${ARGS_TARGET}_lib PROPERTIES
            RULE_LAUNCH_COMPILE "${_ipo_wrapper}")
        target_compile_options(${ARGS_TARGET}_lib PRIVATE
		    $<$<COMPILE_LANGUAGE:Fortran>:-Mextract=lib:${ARGS_TARGET}_lib>
	        $<$<COMPILE_LANGUAGE:Fortran>:-Minline>
    	)
        add_dependencies(${ARGS_TARGET} ${ARGS_TARGET}_lib)
        target_compile_options(${ARGS_TARGET} PRIVATE -Minline=lib:${ARGS_TARGET}_lib,except:f_is_default,except:s_compute_dt,except:my_inquire,except:s_mpi_abort,except:s_mpi_barrier,except:s_prohibit_abort,except:s_int_to_str,except:s_associate_cbc_coefficients_pointers)

        # Exclude m_start_up and m_cbc from cross-file inlining: these are
        # initialization/boundary code that trigger NVHPC 25.x fort2 ICE when
        # too many functions are cross-inlined into them. GPU hot-path files
        # (m_rhs, m_riemann_solvers, m_viscous, m_weno, etc.) keep full IPO.
        foreach(_no_inline_file m_start_up m_cbc)
            set_source_files_properties(
                "${CMAKE_BINARY_DIR}/fypp/${ARGS_TARGET}/${_no_inline_file}.fpp.f90"
                TARGET_DIRECTORY ${ARGS_TARGET}
                PROPERTIES COMPILE_OPTIONS "-Mnoinline"
            )
        endforeach()

        list(PREPEND IPO_TARGETS ${ARGS_TARGET}_lib)
    endif()

    foreach (a_target ${IPO_TARGETS})
        set_target_properties(${a_target} PROPERTIES Fortran_PREPROCESS ON)

        target_include_directories(${a_target} PRIVATE
            "${CMAKE_SOURCE_DIR}/src/common"
            "${CMAKE_SOURCE_DIR}/src/common/include"
            "${CMAKE_SOURCE_DIR}/src/${ARGS_TARGET}")

        if (EXISTS "${CMAKE_SOURCE_DIR}/src/${ARGS_TARGET}/include")
            target_include_directories(${a_target} PRIVATE
                "${CMAKE_SOURCE_DIR}/src/${ARGS_TARGET}/include")
        endif()

        string(TOUPPER "${ARGS_TARGET}" ${ARGS_TARGET}_UPPER)
        target_compile_definitions(
            ${a_target} PRIVATE MFC_${CMAKE_Fortran_COMPILER_ID}
                                MFC_${${ARGS_TARGET}_UPPER}
        )

        if (MFC_MPI AND ARGS_MPI)
            find_package(MPI COMPONENTS Fortran REQUIRED)

            target_compile_definitions(${a_target} PRIVATE MFC_MPI)
            if(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang" AND
               DEFINED ENV{CRAY_MPICH_INC} AND NOT "$ENV{CRAY_MPICH_INC}" STREQUAL "")
                target_compile_options(${a_target} PRIVATE "$ENV{CRAY_MPICH_INC}")
                target_link_libraries(${a_target} PRIVATE $ENV{CRAY_MPICH_LIB})
            else()
                target_link_libraries(${a_target} PRIVATE MPI::MPI_Fortran)
            endif()
        endif()

        if (ARGS_SILO)
            find_package(SILO REQUIRED)
            target_link_libraries(${a_target} PRIVATE SILO::SILO stdc++)
        endif()

        if (ARGS_HDF5)
            find_package(HDF5 REQUIRED)
            target_link_libraries(${a_target} PRIVATE HDF5::HDF5)
        endif()

        if (ARGS_FFTW)
            if ((MFC_OpenACC AND ARGS_OpenACC) OR (MFC_OpenMP AND ARGS_OpenMP))
                if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
                    find_package(CUDAToolkit REQUIRED)
                    target_link_libraries(${a_target} PRIVATE CUDA::cudart CUDA::cufft)
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
                    if(DEFINED ENV{CRAY_HIPFORT_LIB} AND NOT "$ENV{CRAY_HIPFORT_LIB}" STREQUAL "")
                        target_link_libraries(${a_target} PRIVATE $ENV{CRAY_HIPFORT_LIB})
                    else()
                        find_library(HIPFFT_LIB hipfft
                            HINTS "$ENV{OLCF_AFAR_ROOT}/lib" REQUIRED)
                        target_link_libraries(${a_target} PRIVATE ${HIPFFT_LIB})
                    endif()
                else()
                    find_package(hipfort COMPONENTS hipfft CONFIG REQUIRED)
                    target_link_libraries(${a_target} PRIVATE hipfort::hipfft)
                endif()
            else()
                find_package(FFTW REQUIRED)
                target_link_libraries(${a_target} PRIVATE FFTW::FFTW)
            endif()
        endif()

        if (ARGS_LAPACK)
            find_package(LAPACK REQUIRED)
            target_link_libraries(${a_target} PRIVATE LAPACK::LAPACK)
        endif()

        if ((MFC_OpenACC AND ARGS_OpenACC) OR (MFC_OpenMP AND ARGS_OpenMP))
            if ((MFC_OpenACC AND ARGS_OpenACC))
                find_package(OpenACC)

                # This should be equivalent to if (NOT OpenACC_FC_FOUND)
                if (NOT TARGET OpenACC::OpenACC_Fortran)
                    message(FATAL_ERROR "OpenACC + Fortran is unsupported.")
                endif()

                target_link_libraries(${a_target} PRIVATE OpenACC::OpenACC_Fortran)
                target_compile_definitions(${a_target} PRIVATE MFC_OpenACC MFC_GPU)
            elseif((MFC_OpenMP AND ARGS_OpenMP))
                find_package(OpenMP)

                # This should be equivalent to if (NOT OpenACC_FC_FOUND)
                if (NOT TARGET OpenMP::OpenMP_Fortran)
                    message(FATAL_ERROR "OpenMP + Fortran is unsupported.")
                endif()
                # target_link_libraries(${a_target} PRIVATE OpenMP::OpenMP_Fortran)
                target_compile_definitions(${a_target} PRIVATE MFC_OpenMP MFC_GPU)

                if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
                    target_compile_options(${a_target} PRIVATE "-mp=gpu" "-Minfo=mp")
                    target_link_options(${a_target} PRIVATE "-mp=gpu")
                    set_target_properties(${a_target} PROPERTIES Fortran_FLAGS "-mp=gpu -gpu=ccall")
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
                    target_compile_options(${a_target} PRIVATE -fopenmp -fopenmp-targets=spir64)
                    target_link_options(${a_target} PRIVATE -fopenmp -fopenmp-targets=spir64)
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
                    target_compile_options(${a_target} PRIVATE -fopenmp)
                    target_link_options(${a_target} PRIVATE -fopenmp)
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
                    target_compile_options(${a_target} PRIVATE -fopenmp --offload-arch=gfx90a -O3 -fopenmp-assume-threads-oversubscription -fopenmp-assume-teams-oversubscription)
                    target_link_options(${a_target} PRIVATE -fopenmp --offload-arch=gfx90a -flto-partitions=${MFC_BUILD_JOBS})
                endif()
            endif()

            if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
                # FIXME: This should work with other cards than gfx90a ones.
                target_compile_options(${a_target} PRIVATE
                    "-foffload=amdgcn-amdhsa='-march=gfx90a'"
                    "-foffload-options=-lgfortran\ -lm"
                    "-fno-exceptions")
                if (MFC_Fastmath)
                    message(WARNING "--fastmath has no effect with the GNU compiler")
                endif()
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
                foreach (cc ${MFC_CUDA_CC})
                    target_compile_options(${a_target}
                        PRIVATE -gpu=cc${cc}
                    )
                endforeach()

                target_compile_options(${a_target}
                    PRIVATE -gpu=keep,ptxinfo,lineinfo
                )

                if (MFC_Fastmath)
                    target_compile_options(${a_target}
                        PRIVATE -gpu=fastmath
                    )
                endif()

                # GH-200 Unified Memory Support
                if (MFC_Unified)
                    target_compile_options(${a_target}
                        PRIVATE -gpu=mem:unified:managedalloc -cuda
                    )
                    # "This option must appear in both the compile and link lines" -- NVHPC Docs
                    target_link_options(${a_target}
                        PRIVATE -gpu=mem:unified:managedalloc -cuda
                    )
                endif()

                if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelDebug")
                    target_compile_options(${a_target}
                        PRIVATE -gpu=debug
                    )
                endif()
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
                # Frontier Unified Memory Support
                if (MFC_Unified)
                    target_compile_options(${a_target}
                        PRIVATE -DFRONTIER_UNIFIED)
                endif()

                if (MFC_Fastmath)
                    message(WARNING "--fastmath has no effect with the CCE")
                endif()

                find_package(hipfort COMPONENTS hip CONFIG REQUIRED)
                target_link_libraries(${a_target} PRIVATE hipfort::hip hipfort::hipfort-amdgcn)
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")

                if (MFC_Unified)
                    target_compile_options(${a_target}
                        PRIVATE -DFRONTIER_UNIFIED)
                endif()

		        find_library(HIP_LIB amdhip64
                    HINTS "$ENV{OLCF_AFAR_ROOT}/lib" REQUIRED)
                find_library(HIPFORT_AMDGCN_LIB hipfort-amdgcn
                    HINTS "$ENV{OLCF_AFAR_ROOT}/lib" REQUIRED)
                target_include_directories(${a_target} PRIVATE
                    "$ENV{OLCF_AFAR_ROOT}/include/hipfort/amdgcn")
                target_link_libraries(${a_target} PRIVATE
                    ${HIP_LIB} ${HIPFORT_AMDGCN_LIB})

            endif()
        elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
            target_compile_options(${a_target} PRIVATE "SHELL:-h noacc" "SHELL:-x acc")
            if (MFC_Fastmath)
                message(WARNING "--fastmath has no effect with the CCE")
            endif()
        endif()

        if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
            find_package(CUDAToolkit REQUIRED)
            if (TARGET CUDA::nvToolsExt) # CUDA <= 12.8
                target_link_libraries(${a_target} PRIVATE CUDA::nvToolsExt)
            else() # CUDA >= 12.9
                target_link_libraries(${a_target} PRIVATE nvhpcwrapnvtx)
                target_link_options(${a_target} PRIVATE "-cudalib=nvtx3")
            endif()
        endif()
    endforeach()

    install(TARGETS ${ARGS_TARGET} RUNTIME DESTINATION bin)
endfunction()
