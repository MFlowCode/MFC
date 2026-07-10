# cmake/GPU.cmake — compiler flags and GPU configuration

set(FYPP_GCOV_OPTS "")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
    )

    if (MFC_GCov)

        # Warning present in gcc versions >= 12 that is treated as an error
        # This flag doesn't exist in gcc versions < 12
        if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 12)
            add_compile_options(
                -Wno-error=coverage-invalid-line-number
            )
        endif()

        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-fprofile-arcs>
            $<$<COMPILE_LANGUAGE:Fortran>:-ftest-coverage>
        )

        add_link_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-lgcov>
            $<$<COMPILE_LANGUAGE:Fortran>:--coverage>
        )

        # Override Release -O3 with -O1 for gcov: coverage instrumentation is
        # inaccurate at -O3, and aggressive codegen (e.g. AVX-512 FP16 on
        # Granite Rapids) can emit instructions that older assemblers reject.
        set(CMAKE_Fortran_FLAGS_RELEASE "-O1 -DNDEBUG" CACHE STRING "" FORCE)

        # Use gfortran5 line markers so gcov can map coverage to .fpp sources.
        set(FYPP_GCOV_OPTS "--line-marker-format=gfortran5")
    endif()

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
            -Wall
            -Wextra
            -fcheck=all,no-array-temps
            -fbacktrace
            -fimplicit-none
            -fsignaling-nans
            -finit-real=snan
            -finit-integer=-99999999
            -Wconversion
            -Wintrinsic-shadow
            -Wunderflow
            -Wrealloc-lhs
            -Wsurprising
	    )
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options(
            -Og
            -Wall
            -Wextra
            -fcheck=bounds,pointer
            -fbacktrace
            -fimplicit-none
            -fsignaling-nans
            -finit-real=snan
            -finit-integer=-99999999
            -Wconversion
            -Wintrinsic-shadow
            -Wunderflow
            -Wrealloc-lhs
            -Wsurprising
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
        "SHELL:-M 296,878,1391,1069,990,5025,7208,7212,7242"
        "SHELL:-h static" "SHELL:-h keepfiles"
        "SHELL:-h acc_model=auto_async_none"
        "SHELL: -h acc_model=no_fast_addr"
        "SHELL: -h list=adm"
        "SHELL: -munsafe-fp-atomics" # Not unsafe for operations we do
    )

    add_link_options("SHELL:-hkeepfiles")

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
                "SHELL:-h acc_model=auto_async_none"
                "SHELL: -h acc_model=no_fast_addr"
                "SHELL: -K trap=fp" "SHELL: -g" "SHELL: -O0"
        )
        add_link_options("SHELL: -K trap=fp" "SHELL: -g" "SHELL: -O0")
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options(
                "SHELL:-h acc_model=auto_async_none"
                "SHELL: -h acc_model=no_fast_addr"
                "SHELL: -K trap=fp" "SHELL: -g" "SHELL: -O1"
        )
        add_link_options("SHELL: -K trap=fp" "SHELL: -g" "SHELL: -O1")
    endif()

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>
        $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8>
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-O1> $<$<COMPILE_LANGUAGE:Fortran>:-g>)
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
    # AMD/LLVM flang (amdflang). These are CPU-build flags only; the OpenMP-offload GPU
    # build sets its own flags (MFCTargets.cmake) and links via the offload wrapper, so
    # it is left untouched here.
    if (NOT (MFC_OpenMP OR MFC_OpenACC))
        if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelDebug")
            # Ordinary (non-LTO) objects plus static deps built without -fPIC
            # (e.g. libfftw3) hit ld.lld's default-PIE "relocation R_X86_64_32
            # against local symbol". Link a non-PIE executable to accept them.
            # Release links fine via LTO's link-time PIC codegen, so it is left alone.
            add_link_options(-no-pie)
            add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-O1> $<$<COMPILE_LANGUAGE:Fortran>:-g>)
        endif()
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-free>)

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -Og -traceback -debug -check all)
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options(-g -Og -traceback -check bounds)
    endif()
elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
        $<$<COMPILE_LANGUAGE:Fortran>:-Minfo=inline>
        $<$<COMPILE_LANGUAGE:Fortran>:-Minfo=accel>
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-O0>
            $<$<COMPILE_LANGUAGE:Fortran>:-C>
            $<$<COMPILE_LANGUAGE:Fortran>:-g>
            $<$<COMPILE_LANGUAGE:Fortran>:-traceback>
            $<$<COMPILE_LANGUAGE:Fortran>:-Minform=inform>
            $<$<COMPILE_LANGUAGE:Fortran>:-Mbounds>
        )
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-O1>
            $<$<COMPILE_LANGUAGE:Fortran>:-g>
            $<$<COMPILE_LANGUAGE:Fortran>:-traceback>
            $<$<COMPILE_LANGUAGE:Fortran>:-Mbounds>
        )
    endif()

    if (DEFINED ENV{MFC_CUDA_CC})
        string(REGEX MATCHALL "[0-9]+" MFC_CUDA_CC $ENV{MFC_CUDA_CC})
        message(STATUS "Found $MFC_CUDA_CC specified. GPU code will be generated for compute capability(ies) ${MFC_CUDA_CC}.")
    endif()
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    # Processor tuning: Check if we can target the host's native CPU's ISA.
    # Skip for gcov builds — -march=native on newer CPUs (e.g. Granite Rapids)
    # can emit instructions the system assembler doesn't support.
    if (NOT MFC_GCov)
        CHECK_FORTRAN_COMPILER_FLAG("-march=native" SUPPORTS_MARCH_NATIVE)
        if (SUPPORTS_MARCH_NATIVE)
            add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-march=native>)
            # Disable AVX-512 FP16: gfortran >=12 emits vmovw instructions on
            # Granite Rapids CPUs, but binutils <2.38 cannot assemble them.
            # FP16 is unused in MFC's double-precision computations.
            CHECK_FORTRAN_COMPILER_FLAG("-mno-avx512fp16" SUPPORTS_MNO_AVX512FP16)
            if (SUPPORTS_MNO_AVX512FP16)
                add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-mno-avx512fp16>)
            endif()
        else()
            CHECK_FORTRAN_COMPILER_FLAG("-mcpu=native" SUPPORTS_MCPU_NATIVE)
            if (SUPPORTS_MCPU_NATIVE)
                add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-mcpu=native>)
            endif()
        endif()
    endif()

    # Enable LTO/IPO if supported (skip for gcov — LTO interferes with coverage
    # instrumentation and can trigger assembler errors on newer architectures).
    if (MFC_GCov)
        message(STATUS "LTO/IPO disabled for gcov build")
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
        if (MFC_Unified)
            message(STATUS "LTO/IPO is not available with NVHPC using Unified Memory")
        else()
            message(STATUS "Performing IPO using -Mextract followed by -Minline")
            set(NVHPC_USE_TWO_PASS_IPO TRUE)
        endif()
    else()
        CHECK_IPO_SUPPORTED(RESULT SUPPORTS_IPO OUTPUT IPO_ERROR)
        if (SUPPORTS_IPO)
            message(STATUS "Enabled IPO / LTO")
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
	    else()
            message(STATUS "IPO / LTO is NOT available")
        endif()
    endif()
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelDebug")
    add_compile_definitions(MFC_DEBUG)
endif()
