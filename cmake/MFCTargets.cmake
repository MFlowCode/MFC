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

        # Exclude these files from cross-file inlining. GPU hot-path files
        # (m_rhs, m_riemann_solvers, m_viscous, m_weno, etc.) keep full IPO.
        #
        #   m_start_up, m_cbc      initialization/boundary code that trigger
        #                          NVHPC 25.x fort2 ICEs when too many functions
        #                          are cross-inlined into them.
        #   m_pressure_relaxation  worse than an ICE, because it compiles: inlining
        #                          the equation-of-state chain (s_phase_coefficients
        #                          -> s_eos_coefficients -> s_reference_curve,
        #                          s_phase_internal_energy, f_pressure) into the
        #                          six-equation relaxation kernel silently drops the
        #                          internal-energy update, leaving alpha and
        #                          alpha_rho correct and int_en zero.
        foreach(_no_inline_file m_start_up m_cbc m_pressure_relaxation)
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

        if (MFC_NVSHMEM AND ARGS_TARGET STREQUAL "simulation")
            target_compile_definitions(${a_target} PRIVATE MFC_NVSHMEM)
            target_include_directories(${a_target} PRIVATE "${MFC_NVSHMEM_LOWERINGS_DIR}")
        endif()

        if (MFC_LAGRANGIAN_DEVICE_CELL_LIST AND ARGS_TARGET STREQUAL "simulation")
            target_compile_definitions(${a_target} PRIVATE MFC_LAGRANGIAN_DEVICE_CELL_LIST)
        endif()

        if (MFC_MPI AND ARGS_MPI)
            find_package(MPI COMPONENTS Fortran REQUIRED)

            target_compile_definitions(${a_target} PRIVATE MFC_MPI)
            # amdflang (LLVMFlang) needs MPI include/link paths passed explicitly, as
            # its mpi.mod must match the compiler ABI. Prefer the compiler-neutral
            # MFC_FLANG_MPI_* vars; fall back to the legacy CRAY_MPICH_* names (Frontier).
            set(_flang_mpi_inc "$ENV{MFC_FLANG_MPI_INC}")
            set(_flang_mpi_lib "$ENV{MFC_FLANG_MPI_LIB}")
            if("${_flang_mpi_inc}" STREQUAL "")
                set(_flang_mpi_inc "$ENV{CRAY_MPICH_INC}")
                set(_flang_mpi_lib "$ENV{CRAY_MPICH_LIB}")
            endif()
            if(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang" AND NOT "${_flang_mpi_inc}" STREQUAL "")
                if("${_flang_mpi_lib}" STREQUAL "")
                    message(FATAL_ERROR "amdflang MPI include path is set but its link libraries are not. "
                        "Set both MFC_FLANG_MPI_INC and MFC_FLANG_MPI_LIB (or both CRAY_MPICH_INC and CRAY_MPICH_LIB).")
                endif()
                target_compile_options(${a_target} PRIVATE "${_flang_mpi_inc}")
                target_link_libraries(${a_target} PRIVATE ${_flang_mpi_lib})
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
                    # Prefer the neutral MFC_FLANG_HIPFORT_LIB; fall back to CRAY_HIPFORT_LIB (Frontier).
                    set(_flang_hipfort_lib "$ENV{MFC_FLANG_HIPFORT_LIB}")
                    if("${_flang_hipfort_lib}" STREQUAL "")
                        set(_flang_hipfort_lib "$ENV{CRAY_HIPFORT_LIB}")
                    endif()
                    if(NOT "${_flang_hipfort_lib}" STREQUAL "")
                        target_link_libraries(${a_target} PRIVATE ${_flang_hipfort_lib})
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
                    # GPU architecture: gfx90a (MI250X) unless MFC_AMD_GPU_ARCH overrides it.
                    set(_amd_gpu_arch "gfx90a")
                    if (DEFINED ENV{MFC_AMD_GPU_ARCH})
                        set(_amd_gpu_arch "$ENV{MFC_AMD_GPU_ARCH}")
                    endif()
                    # assume-no-nested-parallelism is for correctness, not speed: the two
                    # oversubscription flags are the trigger for llvm/llvm-project#198621
                    # (a suffix of loop iterations silently skipped), which -O3 happens to
                    # mask today. Codegen is unchanged on amdflang; it is a win on upstream
                    # flang.
                    # assume-no-thread-state: no kernel changes OpenMP ICVs on the device, so
                    # the per-kernel thread-state bookkeeping is dropped (~1.2x on AFAR 24.3).
                    target_compile_options(${a_target} PRIVATE
                        -fopenmp
                        --offload-arch=${_amd_gpu_arch}
                        -O3
                        -fopenmp-assume-threads-oversubscription
                        -fopenmp-assume-teams-oversubscription
                        -fopenmp-assume-no-nested-parallelism
                        -fopenmp-assume-no-thread-state)
                    # attributor-max-pi-accesses: amdflang generates device code for the WHOLE
                    # image at link time, and once the image carries enough target regions the
                    # device link's Attributor exceeds its AAPointerInfo access cap on a
                    # heavily-shared object. Pointer information then goes pessimistic and
                    # OpenMPOpt's __kmpc_parallel cleanup fails module-wide: UNTOUCHED kernels
                    # regenerate with 2.4-4.5x worse ISA (register spills, +512 B LDS in every
                    # kernel) whenever ANY kernel is added or removed anywhere in the code.
                    # Raising the cap restores full pointer precision for the whole image and
                    # makes kernel quality independent of unrelated edits, at the price of a
                    # longer device link. See docs/documentation/gpuParallelization.md
                    # ("AMD flang known issues") for the failure signature. Still needed on AFAR 24.3
                    # (cap unchanged; ROCm/llvm-project#4070).
                    target_link_options(${a_target} PRIVATE -fopenmp --offload-arch=${_amd_gpu_arch} -flto-partitions=${MFC_BUILD_JOBS}
                        "SHELL:-Xoffload-linker -mllvm -Xoffload-linker -attributor-max-pi-accesses=16384")

                    # Weak device stubs for symbols AFAR 24.3's device flang runtime
                    # leaves undefined. Built by the drop's clang so the bitcode
                    # matches the device LTO link.
                    get_filename_component(_flang_bin "${CMAKE_Fortran_COMPILER}" DIRECTORY)
                    find_program(MFC_AMDCLANG NAMES amdclang clang HINTS "${_flang_bin}" NO_DEFAULT_PATH REQUIRED)
                    set(_rt_stubs_src "${CMAKE_SOURCE_DIR}/cmake/amdflang_device_stubs.c")
                    set(_rt_stubs "${CMAKE_CURRENT_BINARY_DIR}/${a_target}_amdflang_device_stubs.o")
                    add_custom_command(OUTPUT "${_rt_stubs}"
                        COMMAND "${MFC_AMDCLANG}" -O2 -fopenmp --offload-arch=${_amd_gpu_arch} -c "${_rt_stubs_src}" -o "${_rt_stubs}"
                        DEPENDS "${_rt_stubs_src}")
                    target_sources(${a_target} PRIVATE "${_rt_stubs}")
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

                # Search the AFAR drop first: the ROCm module also ships
                # libhipfort-amdgcn.a, and HINTS lose to CMAKE_PREFIX_PATH, which
                # would pair the drop's .mod files with another flang's archive.
                find_library(HIP_LIB amdhip64
                    PATHS "$ENV{OLCF_AFAR_ROOT}/lib" "$ENV{OLCF_AFAR_ROOT}/lib/llvm/lib"
                    NO_DEFAULT_PATH)
                find_library(HIP_LIB amdhip64 REQUIRED)
                find_library(HIPFORT_AMDGCN_LIB hipfort-amdgcn
                    PATHS "$ENV{OLCF_AFAR_ROOT}/lib/llvm/lib/fortran/flang"
                          "$ENV{OLCF_AFAR_ROOT}/lib" "$ENV{OLCF_AFAR_ROOT}/lib/llvm/lib"
                    NO_DEFAULT_PATH)
                find_library(HIPFORT_AMDGCN_LIB hipfort-amdgcn REQUIRED)
                # The hipfort module dir moved to lib/llvm/include in AFAR 23.x and
                # to lib/llvm/include/fortran/flang in 24.3; keep all three layouts.
                target_include_directories(${a_target} PRIVATE
                    "$ENV{OLCF_AFAR_ROOT}/include/hipfort/amdgcn"
                    "$ENV{OLCF_AFAR_ROOT}/lib/llvm/include/hipfort/amdgcn"
                    "$ENV{OLCF_AFAR_ROOT}/lib/llvm/include/fortran/flang/hipfort/amdgcn")
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
                # Link the Fortran NVTX wrapper library under the name this NVHPC ships.
                get_filename_component(_nvhpc_lib_dir "${CMAKE_Fortran_COMPILER}" DIRECTORY)
                get_filename_component(_nvhpc_lib_dir "${_nvhpc_lib_dir}/../lib" ABSOLUTE)
                if (EXISTS "${_nvhpc_lib_dir}/libnvhpcnvtx.so" OR EXISTS "${_nvhpc_lib_dir}/libnvhpcnvtx.a")
                    set(_nvhpc_nvtx_lib nvhpcnvtx)
                else()
                    set(_nvhpc_nvtx_lib nvhpcwrapnvtx)
                endif()
                target_link_libraries(${a_target} PRIVATE "${_nvhpc_nvtx_lib}")
                target_link_options(${a_target} PRIVATE "-cudalib=nvtx3")
            endif()
        endif()
    endforeach()

    if (MFC_HALO_SINGLE_ROUND AND ARGS_TARGET STREQUAL "simulation")
        # m_mpi_common declares module state under this define, so the IPO extraction
        # library needs it too to keep the module layout identical.
        foreach(a_target ${IPO_TARGETS})
            target_compile_definitions(${a_target} PRIVATE MFC_HALO_SINGLE_ROUND)
        endforeach()
    endif()

    # Target-level feature defines reach both IPO passes: a define the -Mextract
    # library misses lets -Minline splice a differently-compiled body into the
    # executable. MFC_ASYNC_QUEUE stays ungated by language so the generated
    # C++ backend follows the Fortran queue's stream.
    foreach(_define MFC_CBC_NO_HOST_SYNC MFC_IGR_SHALLOW_HALO MFC_CBC_PARALLEL_SCATTER MFC_ASYNC_QUEUE)
        if (${_define} AND ARGS_TARGET STREQUAL "simulation")
            foreach(a_target ${IPO_TARGETS})
                target_compile_definitions(${a_target} PRIVATE ${_define})
            endforeach()
        endif()
    endforeach()

    if (MFC_DIRECT_RIEMANN_FLUX AND ARGS_TARGET STREQUAL "simulation")
        set_property(
            SOURCE
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_rhs.fpp.f90"
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_riemann_state.fpp.f90"
            TARGET_DIRECTORY ${IPO_TARGETS}
            APPEND PROPERTY COMPILE_DEFINITIONS "MFC_DIRECT_RIEMANN_FLUX"
        )
    endif()

    if (MFC_MUSCL_DIRECT_LOADS AND ARGS_TARGET STREQUAL "simulation")
        set_source_files_properties(
            "${CMAKE_BINARY_DIR}/fypp/simulation/m_muscl.fpp.f90"
            TARGET_DIRECTORY ${IPO_TARGETS}
            PROPERTIES COMPILE_DEFINITIONS "MFC_MUSCL_DIRECT_LOADS"
        )
    endif()

    if (MFC_HLLC_REGISTER_TUNING AND ARGS_TARGET STREQUAL "simulation")
        set_source_files_properties(
            "${CMAKE_BINARY_DIR}/fypp/simulation/m_riemann_solver_hllc.fpp.f90"
            TARGET_DIRECTORY ${IPO_TARGETS}
            PROPERTIES COMPILE_DEFINITIONS "MFC_HLLC_REGISTER_TUNING"
        )
    endif()

    # The Fortran preprocessor picks m_weno's portable variant unless a selector is set;
    # AMD compilers always take it (their OpenMP path keeps the reference module).
    if (ARGS_TARGET STREQUAL "simulation" AND NOT ((MFC_FUSED_RHS OR MFC_HLLC_REGISTER_TUNING) AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang"))
        set_source_files_properties(
            "${CMAKE_BINARY_DIR}/fypp/simulation/m_weno.fpp.f90"
            TARGET_DIRECTORY ${IPO_TARGETS}
            PROPERTIES COMPILE_DEFINITIONS "MFC_WENO_PORTABLE"
        )
    endif()

    if (MFC_REGISTER_CAPS AND ARGS_TARGET STREQUAL "simulation")
        # -gpu=maxregcount is a translation-unit option, so each entry names one solver module.
        if (NOT (MFC_OpenACC AND ARGS_OpenACC AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC"))
            message(FATAL_ERROR "MFC_REGISTER_CAPS is set but this build does not use NVHPC OpenACC.")
        endif()
        foreach(_register_cap ${MFC_REGISTER_CAPS})
            string(REPLACE "=" ";" _register_cap "${_register_cap}")
            list(GET _register_cap 0 _register_cap_module)
            list(GET _register_cap 1 _register_cap_value)
            set_source_files_properties(
                "${CMAKE_BINARY_DIR}/fypp/simulation/${_register_cap_module}.fpp.f90"
                TARGET_DIRECTORY ${ARGS_TARGET}
                PROPERTIES COMPILE_OPTIONS "-gpu=maxregcount:${_register_cap_value}"
            )
        endforeach()
    endif()

    if (MFC_WENO_DIRECT_LOADS AND ARGS_TARGET STREQUAL "simulation")
        # THINC reads the cell averages in place under this flag, so its module and both of its callers see the define.
        set_property(
            SOURCE
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_weno.fpp.f90"
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_thinc.fpp.f90"
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_muscl.fpp.f90"
            TARGET_DIRECTORY ${IPO_TARGETS}
            APPEND PROPERTY COMPILE_DEFINITIONS "MFC_WENO_DIRECT_LOADS"
        )
    endif()

    if (MFC_PERSISTENT_BC AND ARGS_TARGET STREQUAL "simulation")
        set_property(
            SOURCE "${CMAKE_BINARY_DIR}/fypp/simulation/m_boundary_common.fpp.f90"
            TARGET_DIRECTORY ${IPO_TARGETS}
            APPEND PROPERTY COMPILE_DEFINITIONS "MFC_PERSISTENT_BC"
        )
    endif()

    if (MFC_FUSED_BUBBLE_SOURCE AND ARGS_TARGET STREQUAL "simulation")
        foreach(a_target ${IPO_TARGETS})
            target_compile_definitions(${a_target} PRIVATE MFC_FUSED_BUBBLE_SOURCE)
        endforeach()
    endif()

    if (MFC_QBMM_DEVICE_INVERSION AND ARGS_TARGET STREQUAL "simulation")
        set_property(
            SOURCE
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_qbmm.fpp.f90"
                "${CMAKE_BINARY_DIR}/fypp/simulation/m_qbmm_openacc.fpp.f90"
            TARGET_DIRECTORY ${IPO_TARGETS}
            APPEND PROPERTY COMPILE_DEFINITIONS "MFC_QBMM_MOM_INV_DEVICE"
        )
    endif()

    if (MFC_ACCURATE_GRIND_TIME AND ARGS_TARGET STREQUAL "simulation")
        foreach (_grind_src m_time_steppers m_start_up)
            set_property(
                SOURCE "${CMAKE_BINARY_DIR}/fypp/simulation/${_grind_src}.fpp.f90"
                TARGET_DIRECTORY ${IPO_TARGETS}
                APPEND PROPERTY COMPILE_DEFINITIONS "MFC_ACCURATE_GRIND_TIME"
            )
        endforeach()
    endif()

    if (MFC_HALO_INDEX_CACHE AND ARGS_TARGET STREQUAL "simulation")
        # m_mpi_common declares module state under this define, as for MFC_HALO_SINGLE_ROUND.
        foreach(a_target ${IPO_TARGETS})
            target_compile_definitions(${a_target} PRIVATE MFC_HALO_INDEX_CACHE)
        endforeach()
    endif()

    if (MFC_ACTIVE_BOX AND ARGS_TARGET STREQUAL "simulation")
        foreach(a_target ${IPO_TARGETS})
            target_compile_definitions(${a_target} PRIVATE MFC_ACTIVE_BOX)
        endforeach()
    endif()

    if (MFC_FUSED_RHS_ENABLED AND ARGS_TARGET STREQUAL "simulation")
        foreach(a_target ${IPO_TARGETS})
            target_compile_definitions(${a_target} PRIVATE MFC_FUSED_RHS)
        endforeach()
        # NVHPC's extraction pass and final compilation see the same selector.
        if (MFC_HALO_SINGLE_ROUND AND MFC_ASYNC_QUEUE)
            set_property(
                SOURCE "${CMAKE_BINARY_DIR}/fypp/simulation/m_rhs.fpp.f90"
                TARGET_DIRECTORY ${IPO_TARGETS}
                APPEND PROPERTY COMPILE_DEFINITIONS "MFC_NATIVE_HALO_OVERLAP"
            )
        endif()
    endif()

    install(TARGETS ${ARGS_TARGET} RUNTIME DESTINATION bin)
endfunction()
