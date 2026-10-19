# Generate the NVSHMEM lowering from the developer-authored MPI/Fypp regions.
set(MFC_NVSHMEM_GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/nvshmem")
set(MFC_NVSHMEM_LOWERINGS_DIR "${CMAKE_CURRENT_BINARY_DIR}/include/nvshmem")
set(MFC_NVSHMEM_MANIFEST "${MFC_NVSHMEM_GENERATED_DIR}/transports.json")
set(MFC_NVSHMEM_HOST_SOURCE "${MFC_NVSHMEM_GENERATED_DIR}/m_nvshmem_backend.cpp")
set(MFC_NVSHMEM_DEVICE_SOURCE "${MFC_NVSHMEM_GENERATED_DIR}/m_nvshmem_device.cu")
file(GLOB MFC_NVSHMEM_MPI_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/common/*.fpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/simulation/*.fpp")
set(MFC_NVSHMEM_MPI_SOURCE_ARGS)
foreach(_mfc_mpi_source IN LISTS MFC_NVSHMEM_MPI_SOURCES)
    list(APPEND MFC_NVSHMEM_MPI_SOURCE_ARGS --mpi-source "${_mfc_mpi_source}")
endforeach()
set(MFC_NVSHMEM_LOWERINGS
    "${MFC_NVSHMEM_LOWERINGS_DIR}/nvshmem_halo_axis_exchange.fpp"
    "${MFC_NVSHMEM_LOWERINGS_DIR}/nvshmem_halo_axis_complete.fpp"
    "${MFC_NVSHMEM_LOWERINGS_DIR}/nvshmem_halo_26_exchange.fpp"
    "${MFC_NVSHMEM_LOWERINGS_DIR}/nvshmem_halo_26_complete.fpp"
    "${MFC_NVSHMEM_LOWERINGS_DIR}/nvshmem_halo_26_channel.fpp"
    "${MFC_NVSHMEM_LOWERINGS_DIR}/nvshmem_igr_exchange.fpp"
    "${MFC_NVSHMEM_LOWERINGS_DIR}/nvshmem_igr_complete.fpp")

add_custom_command(
    OUTPUT "${MFC_NVSHMEM_MANIFEST}" "${MFC_NVSHMEM_HOST_SOURCE}" "${MFC_NVSHMEM_DEVICE_SOURCE}" ${MFC_NVSHMEM_LOWERINGS}
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${MFC_NVSHMEM_GENERATED_DIR}" "${MFC_NVSHMEM_LOWERINGS_DIR}"
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/src/nvshmem/gen_nvshmem.py"
            ${MFC_NVSHMEM_MPI_SOURCE_ARGS}
            --manifest-output "${MFC_NVSHMEM_MANIFEST}"
            --lowerings-dir "${MFC_NVSHMEM_LOWERINGS_DIR}"
            --host-template "${CMAKE_CURRENT_SOURCE_DIR}/src/nvshmem/templates/nvshmem_backend.cpp.in"
            --device-template "${CMAKE_CURRENT_SOURCE_DIR}/src/nvshmem/templates/nvshmem_device.cu.in"
            --host-output "${MFC_NVSHMEM_HOST_SOURCE}"
            --device-output "${MFC_NVSHMEM_DEVICE_SOURCE}"
    DEPENDS
        ${MFC_NVSHMEM_MPI_SOURCES}
        "${CMAKE_CURRENT_SOURCE_DIR}/src/nvshmem/gen_nvshmem.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/nvshmem/nvshmem_ir.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/nvshmem/templates/nvshmem_backend.cpp.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/nvshmem/templates/nvshmem_device.cu.in"
    COMMENT "Lowering developer MPI/Fypp regions to NVSHMEM NBI backends"
    VERBATIM)

set(MFC_NVSHMEM_HEADERS)
if (MFC_NVSHMEM AND (MFC_SIMULATION OR MFC_ALL))
    foreach(_lowering IN LISTS MFC_NVSHMEM_LOWERINGS)
        string(REGEX REPLACE "\\.fpp$" ".inc" _header "${_lowering}")
        add_custom_command(OUTPUT "${_header}"
            COMMAND ${FYPP_EXE} -m re
                -I "${CMAKE_BINARY_DIR}/include/simulation"
                -I "${CMAKE_SOURCE_DIR}/src/common/include"
                -D MFC_SIMULATION -D MFC_COMPILER="${CMAKE_Fortran_COMPILER_ID}"
                -D MFC_CASE_OPTIMIZATION=False -D chemistry=False
                -D MFC_ASYNC_QUEUE=${FYPP_ASYNC_QUEUE}
                --no-folding --line-length=999 "${_lowering}" "${_header}"
            DEPENDS "${_lowering}" "${CMAKE_BINARY_DIR}/include/simulation/case.fpp"
                "${CMAKE_SOURCE_DIR}/src/common/include/macros.fpp"
                "${CMAKE_SOURCE_DIR}/src/common/include/shared_parallel_macros.fpp"
                "${CMAKE_SOURCE_DIR}/src/common/include/parallel_macros.fpp"
                "${CMAKE_SOURCE_DIR}/src/common/include/acc_macros.fpp"
                "${CMAKE_SOURCE_DIR}/src/common/include/omp_macros.fpp"
            VERBATIM)
        list(APPEND MFC_NVSHMEM_HEADERS "${_header}")
    endforeach()
endif()
add_custom_target(mfc_nvshmem_codegen
                  DEPENDS "${MFC_NVSHMEM_MANIFEST}" "${MFC_NVSHMEM_HOST_SOURCE}" "${MFC_NVSHMEM_DEVICE_SOURCE}" ${MFC_NVSHMEM_LOWERINGS} ${MFC_NVSHMEM_HEADERS})
