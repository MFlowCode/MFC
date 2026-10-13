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

    !> Multigrid hierarchy, flattened, with one ghost layer in each active direction; level lv starts after mg_off(lv), x fastest.
    !! Every level coarsens by two (an odd size folds its last cell into the last coarse cell) down to the bottom level, the first
    !! with at most mg_bottom_max cells in all (or one cell per rank). PCG's search direction follows the levels in mg_e, at
    !! mg_poff, laid out as level 1
    integer                             :: mg_nlev, mg_gx, mg_gy, mg_gz, mg_poff
    integer, dimension(mg_maxlev)       :: mg_nx, mg_ny, mg_nz, mg_off
    real(wp), allocatable, dimension(:) :: mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r
    $:GPU_DECLARE(create='[mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r]')
    !> Flat index of cell (i, j, k) on level lv: mg_o(lv) + k*mg_sz(lv) + j*mg_sy(lv) + i (f_mg_ix, m_multigrid, beside the kernels
    !! that call it, as a call across modules is not inlined)
    integer, dimension(mg_maxlev) :: mg_o, mg_sy, mg_sz
    $:GPU_DECLARE(create='[mg_o, mg_sy, mg_sz]')
    !> Neighbor across side q = 1..6, (x, y, z) x (low, high): a rank, this rank for a periodic seam it owns alone, or -1
    integer, dimension(6) :: mg_nbr

contains

    !> Levels, their storage, and the neighbor across each side
    impure subroutine s_initialize_multigrid_levels()

        integer  :: lv, tot
        real(wp) :: rlev

        ! Halve every direction (rounding down) until each rank holds one cell; all ranks take the global level count, so a rank
        ! that reaches one cell early keeps it, and neighbors, which share their tangential sizes, stay aligned. The hierarchy then
        ! stops at the first level small enough for the dense bottom solve

        mg_gx = 1; mg_gy = merge(1, 0, n > 0); mg_gz = merge(1, 0, p > 0)
        mg_nx(1) = m + 1; mg_ny(1) = n + 1; mg_nz(1) = p + 1
        lv = 1
        do while (max(mg_nx(lv), mg_ny(lv), mg_nz(lv)) > 1)
            mg_nx(lv + 1) = max(mg_nx(lv)/2, 1); mg_ny(lv + 1) = max(mg_ny(lv)/2, 1); mg_nz(lv + 1) = max(mg_nz(lv)/2, 1)
            lv = lv + 1
        end do
        call s_mpi_allreduce_max(real(lv, wp), rlev)
        mg_nlev = nint(rlev)
        do lv = lv + 1, mg_nlev
            mg_nx(lv) = 1; mg_ny(lv) = 1; mg_nz(lv) = 1
        end do
        do lv = 1, mg_nlev
            call s_mpi_allreduce_max(real(mg_nx(lv)*mg_ny(lv)*mg_nz(lv), wp), rlev)
            if (rlev*num_procs <= max(mg_bottom_max, num_procs)) exit
        end do
        mg_nlev = min(lv, mg_nlev)
        tot = 0
        do lv = 1, mg_nlev
            mg_off(lv) = tot
            mg_sy(lv) = mg_nx(lv) + 2*mg_gx; mg_sz(lv) = mg_sy(lv)*(mg_ny(lv) + 2*mg_gy)
            mg_o(lv) = tot + mg_gz*mg_sz(lv) + mg_gy*mg_sy(lv) + mg_gx + 1
            tot = tot + (mg_nx(lv) + 2*mg_gx)*(mg_ny(lv) + 2*mg_gy)*(mg_nz(lv) + 2*mg_gz)
        end do
        mg_poff = tot
        $:GPU_UPDATE(device='[mg_o, mg_sy, mg_sz]')
        @:ALLOCATE(mg_d(tot), mg_kx(tot), mg_ky(tot), mg_kz(tot), mg_f(tot), mg_r(tot))
        @:ALLOCATE(mg_e(tot + (mg_nx(1) + 2*mg_gx)*(mg_ny(1) + 2*mg_gy)*(mg_nz(1) + 2*mg_gz)))

        #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            mg_nbr(2*${D}$ - 1) = f_seam_rank(bc_${XYZ}$%beg, num_dims >= ${D}$)
            mg_nbr(2*${D}$) = f_seam_rank(bc_${XYZ}$%end, num_dims >= ${D}$)
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

        @:DEALLOCATE(mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r)

    end subroutine s_finalize_multigrid_levels

end module m_multigrid_levels
