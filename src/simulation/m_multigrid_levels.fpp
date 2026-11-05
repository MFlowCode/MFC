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

    private; public :: s_initialize_multigrid_levels, s_finalize_multigrid_levels, mg_num_lev, mg_ghost_x, mg_ghost_y, &
        & mg_ghost_z, mg_nx, mg_ny, mg_nz, mg_lev, mg_nbr_rank

    !> Multigrid hierarchy: level lv holds mg_nx(lv) x mg_ny(lv) x mg_nz(lv) cells with one ghost layer in each active direction.
    !! Every level coarsens by two (an odd size folds its last cell into the last coarse cell) down to the bottom level, the first
    !! with at most mg_bottom_max cells in all (or one cell per rank)
    integer                              :: mg_num_lev, mg_ghost_x, mg_ghost_y, mg_ghost_z
    integer, dimension(mg_maxlev)        :: mg_nx, mg_ny, mg_nz
    type(mg_level), dimension(mg_maxlev) :: mg_lev
    $:GPU_DECLARE(create='[mg_lev]')
    !> Neighbor across side q = 1..6, (x, y, z) x (low, high): a rank, this rank for a periodic seam it owns alone, or -1
    integer, dimension(6) :: mg_nbr_rank

contains

    !> Levels, their storage (with the K-cycle vectors of the coarse levels if k_vectors), and the neighbor across each side
    impure subroutine s_initialize_multigrid_levels(k_vectors)

        logical, intent(in)   :: k_vectors
        integer               :: lv, cells
        integer, dimension(3) :: lo, hi

        ! Halve every direction (rounding down) until each rank holds one cell; all ranks take the global level count, so a rank
        ! that reaches one cell early keeps it, and neighbors, which share their tangential sizes, stay aligned. The hierarchy then
        ! stops at the first level small enough for the dense bottom solve

        mg_ghost_x = 1
        mg_ghost_y = merge(1, 0, n > 0)
        mg_ghost_z = merge(1, 0, p > 0)
        mg_nx(1) = m + 1
        mg_ny(1) = n + 1
        mg_nz(1) = p + 1
        lv = 1
        do while (max(mg_nx(lv), mg_ny(lv), mg_nz(lv)) > 1)
            mg_nx(lv + 1) = max(mg_nx(lv)/2, 1)
            mg_ny(lv + 1) = max(mg_ny(lv)/2, 1)
            mg_nz(lv + 1) = max(mg_nz(lv)/2, 1)
            lv = lv + 1
        end do
        call s_mpi_allreduce_integer_max(lv, mg_num_lev)
        do lv = lv + 1, mg_num_lev
            mg_nx(lv) = 1
            mg_ny(lv) = 1
            mg_nz(lv) = 1
        end do
        do lv = 1, mg_num_lev
            call s_mpi_allreduce_integer_max(mg_nx(lv)*mg_ny(lv)*mg_nz(lv), cells)
            if (cells <= max(mg_bottom_max, num_procs)/num_procs) exit
        end do
        mg_num_lev = min(lv, mg_num_lev)
        do lv = 1, mg_num_lev
            lo = -[mg_ghost_x, mg_ghost_y, mg_ghost_z]
            hi = [mg_nx(lv), mg_ny(lv), mg_nz(lv)] - 1 - lo
            #:for VAR in ['diag', 'cond_x', 'cond_y', 'cond_z', 'corr', 'rhs', 'res']
                @:ALLOCATE(mg_lev(lv)%${VAR}$(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
            #:endfor
            if (k_vectors .and. lv > 1) then
                #:for VAR in ['k_rhs', 'k_dir', 'k_Adir']
                    @:ALLOCATE(mg_lev(lv)%${VAR}$(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
                #:endfor
            end if
            @:ACC_SETUP_mg_levels(mg_lev(lv))
        end do

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

        integer :: lv

        do lv = 1, mg_num_lev
            #:for VAR in ['diag', 'cond_x', 'cond_y', 'cond_z', 'corr', 'rhs', 'res', 'k_rhs', 'k_dir', 'k_Adir']
                if (associated(mg_lev(lv)%${VAR}$)) then
                    @:DEALLOCATE(mg_lev(lv)%${VAR}$)
                end if
            #:endfor
        end do

    end subroutine s_finalize_multigrid_levels

end module m_multigrid_levels
