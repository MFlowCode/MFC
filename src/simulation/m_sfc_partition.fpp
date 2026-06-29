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
        & s_report_sfc_partition, n_tiles_x, n_tiles_y, n_tiles_z, n_tiles, tile_weight, tile_rank, sfc_order

    integer               :: n_tiles_x, n_tiles_y, n_tiles_z, n_tiles
    real(wp), allocatable :: tile_weight(:)  !< global per-tile aggregated cost (linear index)
    integer, allocatable  :: tile_rank(:)    !< proposed owning rank per tile
    integer, allocatable  :: sfc_order(:)    !< tile linear indices in Morton order

contains

    impure subroutine s_initialize_sfc_partition_module

        if (.not. sfc_partition_wrt) return
        n_tiles_x = (m_glb + 1 + partition_tile_size - 1)/partition_tile_size
        n_tiles_y = (n_glb + 1 + partition_tile_size - 1)/partition_tile_size
        n_tiles_z = (p_glb + 1 + partition_tile_size - 1)/partition_tile_size
        n_tiles = n_tiles_x*n_tiles_y*n_tiles_z
        @:ALLOCATE(tile_weight(0:n_tiles - 1))
        @:ALLOCATE(tile_rank(0:n_tiles - 1))
        @:ALLOCATE(sfc_order(1:n_tiles))

    end subroutine s_initialize_sfc_partition_module

    impure subroutine s_finalize_sfc_partition_module

        if (.not. sfc_partition_wrt) return
        @:DEALLOCATE(tile_weight)
        @:DEALLOCATE(tile_rank)
        @:DEALLOCATE(sfc_order)

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

        call s_build_sfc_order(sfc_order)

    end subroutine s_compute_sfc_partition

    !> Returns the 63-bit Morton code for tile coordinates (tx, ty, tz).
    pure function f_morton(tx, ty, tz) result(code)

        integer, intent(in) :: tx, ty, tz
        integer(kind=8)     :: code, x, y, z
        integer             :: b

        x = int(tx, 8); y = int(ty, 8); z = int(tz, 8); code = 0_8
        do b = 0, 20
            code = ior(code, ishft(iand(ishft(x, -b), 1_8), 3*b))
            code = ior(code, ishft(iand(ishft(y, -b), 1_8), 3*b + 1))
            code = ior(code, ishft(iand(ishft(z, -b), 1_8), 3*b + 2))
        end do

    end function f_morton

    !> Fills order(1:n_tiles) with tile linear indices sorted by Morton code.
    impure subroutine s_build_sfc_order(order)

        integer, intent(out)         :: order(:)
        integer(kind=8), allocatable :: code(:)
        integer                      :: tx, ty, tz, t, i, jmin
        integer(kind=8)              :: cmin
        logical, allocatable         :: used(:)

        allocate (code(0:n_tiles - 1), used(0:n_tiles - 1)); used = .false.
        do tz = 0, n_tiles_z - 1
            do ty = 0, n_tiles_y - 1
                do tx = 0, n_tiles_x - 1
                    t = (tz*n_tiles_y + ty)*n_tiles_x + tx
                    code(t) = f_morton(tx, ty, tz)
                end do
            end do
        end do
        ! selection by min code (n_tiles is modest; O(n_tiles^2) acceptable, or replace with a sort)
        do i = 1, n_tiles
            cmin = huge(0_8); jmin = -1
            do t = 0, n_tiles - 1
                if (.not. used(t) .and. code(t) < cmin) then; cmin = code(t); jmin = t; end if
            end do
            order(i) = jmin; used(jmin) = .true.
        end do
        deallocate (code, used)

    end subroutine s_build_sfc_order

    !> Task 5 fills this. Stub: no-op.
    impure subroutine s_report_sfc_partition

    end subroutine s_report_sfc_partition

end module m_sfc_partition
