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
    public :: s_initialize_sfc_partition_module, s_finalize_sfc_partition_module, s_compute_sfc_partition, &
        & s_report_sfc_partition, n_tiles_x, n_tiles_y, n_tiles_z, n_tiles, tile_weight, tile_rank

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

    !> Aggregate per-cell load weights into global per-tile weights via MPI_ALLREDUCE.
    impure subroutine s_compute_sfc_partition(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), allocatable                               :: tile_local(:)
        integer                                             :: j, k, l, gx, gy, gz, tx, ty, tz, t, ierr
        integer                                             :: sfc_start(3)  !< bounds-safe copy of start_idx (0 for inactive dims)

        if (.not. sfc_partition_wrt) return

        call s_compute_load_weight(q_cons_vf)
        $:GPU_UPDATE(host='[load_weight%sf]')

        ! Build a 3-element start offset safe to access regardless of num_dims.
        ! start_idx is allocated (1:num_dims) only; higher dims are always 0.
        sfc_start = 0
        sfc_start(1) = start_idx(1)
        if (num_dims >= 2) sfc_start(2) = start_idx(2)
        if (num_dims >= 3) sfc_start(3) = start_idx(3)

        allocate (tile_local(0:n_tiles - 1)); tile_local = 0._wp
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    gx = sfc_start(1) + j
                    gy = sfc_start(2) + k
                    gz = sfc_start(3) + l
                    tx = gx/partition_tile_size; ty = gy/partition_tile_size; tz = gz/partition_tile_size
                    t = (tz*n_tiles_y + ty)*n_tiles_x + tx
                    tile_local(t) = tile_local(t) + real(load_weight%sf(j, k, l), wp)
                end do
            end do
        end do
#ifdef MFC_MPI
        call MPI_ALLREDUCE(tile_local, tile_weight, n_tiles, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        tile_weight = tile_local
#endif
        deallocate (tile_local)

    end subroutine s_compute_sfc_partition

    !> Task 5 fills this. Stub: no-op.
    impure subroutine s_report_sfc_partition

    end subroutine s_report_sfc_partition

end module m_sfc_partition
