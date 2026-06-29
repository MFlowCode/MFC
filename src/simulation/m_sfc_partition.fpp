!>
!!@file
!!@brief Contains module m_sfc_partition

#:include 'macros.fpp'

!> @brief Analysis-only weighted space-filling-curve partitioner (Tier-2 sub-project B).
module m_sfc_partition

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_load_weight, only: load_weight, s_compute_load_weight

    implicit none

    private
    public :: s_initialize_sfc_partition_module, s_finalize_sfc_partition_module, s_compute_sfc_partition, s_report_sfc_partition

    integer               :: n_tiles_x, n_tiles_y, n_tiles_z, n_tiles
    real(wp), allocatable :: tile_weight(:)  !< global per-tile aggregated cost (linear index)
    integer, allocatable  :: tile_rank(:)    !< proposed owning rank per tile

contains

    impure subroutine s_initialize_sfc_partition_module

        if (.not. sfc_partition_wrt) return
        n_tiles_x = (m_glb + 1 + partition_tile_size - 1)/partition_tile_size
        n_tiles_y = (n_glb + 1 + partition_tile_size - 1)/partition_tile_size
        n_tiles_z = (p_glb + 1 + partition_tile_size - 1)/partition_tile_size
        n_tiles = n_tiles_x*n_tiles_y*n_tiles_z
        @:ALLOCATE(tile_weight(0:n_tiles - 1))
        @:ALLOCATE(tile_rank(0:n_tiles - 1))

    end subroutine s_initialize_sfc_partition_module

    impure subroutine s_finalize_sfc_partition_module

        if (.not. sfc_partition_wrt) return
        @:DEALLOCATE(tile_weight)
        @:DEALLOCATE(tile_rank)

    end subroutine s_finalize_sfc_partition_module

    !> Task 2-4 fill this. Stub: no-op.
    impure subroutine s_compute_sfc_partition(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf

    end subroutine s_compute_sfc_partition

    !> Task 5 fills this. Stub: no-op.
    impure subroutine s_report_sfc_partition

    end subroutine s_report_sfc_partition

end module m_sfc_partition
