!>
!! @file
!! @brief Contains module m_multigrid_levels

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Layout and storage of the multigrid hierarchy, shared by the cycle (m_multigrid), the ghost exchange
!! (m_multigrid_exchange) and the gathered bottom solve (m_multigrid_bottom)
module m_multigrid_levels

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy

    implicit none

    private :: f_seam_rank

    !> Multigrid hierarchy, flattened, with one ghost layer in each active direction; level lv starts after mg_lev_off(lv), x
    !! fastest. Every level coarsens by two (an odd size folds its last cell into the last coarse cell) down to the bottom level,
    !! the first with at most mg_bottom_max cells in all (or one cell per rank). PCG's search direction follows the levels in
    !! mg_corr, at mg_dir_off, laid out as level 1
    integer                             :: mg_num_lev, mg_ghost_x, mg_ghost_y, mg_ghost_z, mg_dir_off
    integer, dimension(mg_maxlev)       :: mg_nx, mg_ny, mg_nz, mg_lev_off
    real(wp), allocatable, dimension(:) :: mg_diag, mg_cond_x, mg_cond_y, mg_cond_z, mg_corr, mg_rhs, mg_res
    $:GPU_DECLARE(create='[mg_diag, mg_cond_x, mg_cond_y, mg_cond_z, mg_corr, mg_rhs, mg_res]')
    !> Flat index of cell (j, k, l) on level lv: mg_origin(lv) + l*mg_stride_z(lv) + k*mg_stride_y(lv) + j (f_mg_ix, m_multigrid,
    !! beside the kernels that call it, as a call across modules is not inlined)
    integer, dimension(mg_maxlev) :: mg_origin, mg_stride_y, mg_stride_z
    $:GPU_DECLARE(create='[mg_origin, mg_stride_y, mg_stride_z]')
    !> Neighbor across side q = 1..6, (x, y, z) x (low, high): a rank, this rank for a periodic seam it owns alone, or -1
    integer, dimension(6) :: mg_nbr_rank

contains

    !> Levels, their storage, and the neighbor across each side
    impure subroutine s_initialize_multigrid_levels()

        integer  :: lv, tot
        real(wp) :: rlev

        ! Halve every direction (rounding down) until each rank holds one cell; all ranks take the global level count, so a rank
        ! that reaches one cell early keeps it, and neighbors, which share their tangential sizes, stay aligned. The hierarchy then
        ! stops at the first level small enough for the dense bottom solve

        mg_ghost_x = 1; mg_ghost_y = merge(1, 0, n > 0); mg_ghost_z = merge(1, 0, p > 0)
        mg_nx(1) = m + 1; mg_ny(1) = n + 1; mg_nz(1) = p + 1
        lv = 1
        do while (max(mg_nx(lv), mg_ny(lv), mg_nz(lv)) > 1)
            mg_nx(lv + 1) = max(mg_nx(lv)/2, 1); mg_ny(lv + 1) = max(mg_ny(lv)/2, 1); mg_nz(lv + 1) = max(mg_nz(lv)/2, 1)
            lv = lv + 1
        end do
        call s_mpi_allreduce_max(real(lv, wp), rlev)
        mg_num_lev = nint(rlev)
        do lv = lv + 1, mg_num_lev
            mg_nx(lv) = 1; mg_ny(lv) = 1; mg_nz(lv) = 1
        end do
        do lv = 1, mg_num_lev
            call s_mpi_allreduce_max(real(mg_nx(lv)*mg_ny(lv)*mg_nz(lv), wp), rlev)
            if (rlev*num_procs <= max(mg_bottom_max, num_procs)) exit
        end do
        mg_num_lev = min(lv, mg_num_lev)
        tot = 0
        do lv = 1, mg_num_lev
            mg_lev_off(lv) = tot
            mg_stride_y(lv) = mg_nx(lv) + 2*mg_ghost_x; mg_stride_z(lv) = mg_stride_y(lv)*(mg_ny(lv) + 2*mg_ghost_y)
            mg_origin(lv) = tot + mg_ghost_z*mg_stride_z(lv) + mg_ghost_y*mg_stride_y(lv) + mg_ghost_x + 1
            tot = tot + (mg_nx(lv) + 2*mg_ghost_x)*(mg_ny(lv) + 2*mg_ghost_y)*(mg_nz(lv) + 2*mg_ghost_z)
        end do
        mg_dir_off = tot
        $:GPU_UPDATE(device='[mg_origin, mg_stride_y, mg_stride_z]')
        @:ALLOCATE(mg_diag(tot), mg_cond_x(tot), mg_cond_y(tot), mg_cond_z(tot), mg_rhs(tot), mg_res(tot))
        @:ALLOCATE(mg_corr(tot + (mg_nx(1) + 2*mg_ghost_x)*(mg_ny(1) + 2*mg_ghost_y)*(mg_nz(1) + 2*mg_ghost_z)))

        #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            mg_nbr_rank(2*${D}$ - 1) = f_seam_rank(bc_${XYZ}$%beg, num_dims >= ${D}$)
            mg_nbr_rank(2*${D}$) = f_seam_rank(bc_${XYZ}$%end, num_dims >= ${D}$)
        #:endfor

    end subroutine s_initialize_multigrid_levels

    !> Neighbor across a domain side from its bc value: a rank, this rank for a periodic seam it owns alone, or -1
    pure function f_seam_rank(bc, active) result(rank)

        integer, intent(in) :: bc
        logical, intent(in) :: active
        integer             :: rank

        rank = -1
        if (.not. active) return
        if (bc >= 0) rank = bc
        if (bc == BC_PERIODIC) rank = proc_rank

    end function f_seam_rank

    impure subroutine s_finalize_multigrid_levels()

        @:DEALLOCATE(mg_diag, mg_cond_x, mg_cond_y, mg_cond_z, mg_corr, mg_rhs, mg_res)

    end subroutine s_finalize_multigrid_levels

end module m_multigrid_levels
