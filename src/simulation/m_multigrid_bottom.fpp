!>
!! @file
!! @brief Contains module m_multigrid_bottom

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Bottom solve of the multigrid hierarchy: the last device level, one cell per rank or a few cells in all, gathered to every
!! rank and solved there by a host hierarchy of its own
module m_multigrid_bottom

#ifdef MFC_MPI
    use mpi
#endif
    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_multigrid_levels

    implicit none

    private; public :: s_initialize_multigrid_bottom, s_mg_bottom_build, s_mg_bottom_solve, s_finalize_multigrid_bottom

    !> Bottom level, gathered to every rank: rank r's bot_rank_ncells(r + 1) cells (a bot_rank_dims(:, r + 1) block, x fastest)
    !! start at bot_rank_disp(r + 1) + 1. Past mg_bottom_max cells (one per rank) every rank coarsens it in bot_lev and runs one
    !! V-cycle, unless bot_exact; the last level is factored in bot_chol, or solved by CG (bot_cg) if too large
    integer                                     :: bot_nlev
    logical                                     :: bot_exact  !< bottom solved exactly (proj_mg_bottom, or by rank count)
    logical                                     :: bot_cg     !< exactly, by CG, as it is too large to factor densely
    type(mg_bottom_level), dimension(mg_maxlev) :: bot_lev
    integer, allocatable, dimension(:)          :: bot_rank_ncells, bot_rank_disp
    integer, allocatable, dimension(:,:)        :: bot_rank_dims
    real(wp), allocatable, dimension(:,:)       :: bot_chol

contains

    !> Every rank's bottom-level block size, for the global numbering of the bottom cells, and the bottom hierarchy's levels; exact
    !! keeps the bottom whole (bot_exact)
    impure subroutine s_initialize_multigrid_bottom(exact)

        logical, intent(in)              :: exact
        integer, dimension(3, num_procs) :: crd, crd_c  !< rank-grid coordinates of level lv's cells, and of level lv + 1's
        integer, dimension(3)            :: coords
        integer                          :: r, lv, i, n, ierr

        bot_exact = exact
        allocate (bot_rank_ncells(num_procs), bot_rank_disp(num_procs), bot_rank_dims(3, num_procs))
        bot_rank_dims(:,proc_rank + 1) = [mg_nx(mg_num_lev), mg_ny(mg_num_lev), mg_nz(mg_num_lev)]
#ifdef MFC_MPI
        call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, bot_rank_dims, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
        bot_lev(1)%ncells = 0
        do r = 1, num_procs
            bot_rank_ncells(r) = product(bot_rank_dims(:,r))
            bot_rank_disp(r) = bot_lev(1)%ncells
            bot_lev(1)%ncells = bot_lev(1)%ncells + bot_rank_ncells(r)
        end do

        bot_nlev = 1
        bot_lev(1)%dims = [num_procs_x, num_procs_y, num_procs_z]
        do while (bot_lev(bot_nlev)%ncells > mg_bottom_max .and. .not. bot_exact)
            bot_lev(bot_nlev + 1)%dims = (bot_lev(bot_nlev)%dims + 1)/2
            bot_lev(bot_nlev + 1)%ncells = product(bot_lev(bot_nlev + 1)%dims)
            bot_nlev = bot_nlev + 1
        end do
        do lv = 1, bot_nlev
            n = bot_lev(lv)%ncells
            allocate (bot_lev(lv)%diag(n), bot_lev(lv)%cond(6, n), bot_lev(lv)%nbr(6, n), bot_lev(lv)%parent(n), &
                      & bot_lev(lv)%x(n), bot_lev(lv)%rhs(n))
        end do
        bot_cg = bot_lev(bot_nlev)%ncells > mg_bottom_max
        if (.not. bot_cg) allocate (bot_chol(bot_lev(bot_nlev)%ncells, bot_lev(bot_nlev)%ncells))
        if (bot_nlev == 1) return

        if (any(bot_rank_ncells /= 1) .or. product(bot_lev(1)%dims) /= bot_lev(1)%ncells) then
            call s_mpi_abort('multigrid: a bottom level past mg_bottom_max cells must hold one cell per rank')
        end if
        coords = 0
        coords(1:num_dims) = proc_coords(1:num_dims)
        crd(:,proc_rank + 1) = coords
#ifdef MFC_MPI
        call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, crd, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
        do lv = 1, bot_nlev - 1
            do i = 1, bot_lev(lv)%ncells
                coords = crd(:,i)/2
                bot_lev(lv)%parent(i) = 1 + coords(1) + bot_lev(lv + 1)%dims(1)*(coords(2) + bot_lev(lv + 1)%dims(2)*coords(3))
                crd_c(:,bot_lev(lv)%parent(i)) = coords
            end do
            crd = crd_c
        end do

    end subroutine s_initialize_multigrid_bottom

    !> Global bottom index of the neighbor of this rank's bottom cell (i, j, k) across side q, or 0 where a boundary has no seam
    pure function f_bot_nbr(q, i, j, k) result(g)

        integer, intent(in)   :: q, i, j, k
        integer               :: g
        integer, dimension(3) :: cell
        integer               :: dir, r

        cell = [i, j, k]
        dir = (q + 1)/2
        r = proc_rank
        g = 0
        if (mod(q, 2) == 1) then
            cell(dir) = cell(dir) - 1
            if (cell(dir) < 0) then
                r = mg_nbr_rank(q)
                if (r < 0) return
                cell(dir) = bot_rank_dims(dir, r + 1) - 1
            end if
        else
            cell(dir) = cell(dir) + 1
            if (cell(dir) >= bot_rank_dims(dir, proc_rank + 1)) then
                r = mg_nbr_rank(q)
                if (r < 0) return
                cell(dir) = 0
            end if
        end if
        g = bot_rank_disp(r + 1) + (cell(3)*bot_rank_dims(2, r + 1) + cell(2))*bot_rank_dims(1, r + 1) + cell(1) + 1

    end function f_bot_nbr

    !> Gather every rank's bottom rows (cell_rows(:, c): diagonal sum, then the conductance across each side, x fastest) into the
    !! global stencil, coarsen it and factor the last level; a face with both sides in one cell couples it to itself and cancels
    impure subroutine s_mg_bottom_build(cell_rows)

        real(wp), dimension(7, bot_rank_ncells(proc_rank + 1)), intent(in) :: cell_rows
        !> Row 1 is the diagonal sum; row 2q the neighbor across side q, and 2q + 1 its conductance
        real(wp), dimension(13, bot_rank_ncells(proc_rank + 1)) :: rows_loc
        real(wp), dimension(13, bot_lev(1)%ncells) :: rows_glb
        integer :: i, j, k, c, q, g, n, col, ierr, lv, parent_i, parent_j

        c = 0
        do k = 0, bot_rank_dims(3, proc_rank + 1) - 1
            do j = 0, bot_rank_dims(2, proc_rank + 1) - 1
                do i = 0, bot_rank_dims(1, proc_rank + 1) - 1
                    c = c + 1
                    rows_loc(1, c) = cell_rows(1, c)
                    do q = 1, 6
                        rows_loc(2*q, c) = 0._wp
                        if (q <= 2*num_dims) rows_loc(2*q, c) = real(f_bot_nbr(q, i, j, k), wp)
                        rows_loc(2*q + 1, c) = cell_rows(q + 1, c)
                    end do
                end do
            end do
        end do

        rows_glb(:,1:c) = rows_loc
#ifdef MFC_MPI
        call MPI_ALLGATHERV(rows_loc, 13*c, mpi_p, rows_glb, 13*bot_rank_ncells, 13*bot_rank_disp, mpi_p, MPI_COMM_WORLD, ierr)
#endif

        do i = 1, bot_lev(1)%ncells
            bot_lev(1)%diag(i) = rows_glb(1, i)
            do q = 1, 6
                g = nint(rows_glb(2*q, i))
                ! No neighbour: a domain face, whose conductance (a pressure outlet's; zero on a wall) couples to a fixed ghost
                if (g == 0) bot_lev(1)%diag(i) = bot_lev(1)%diag(i) + rows_glb(2*q + 1, i)
                if (g == i) g = 0
                bot_lev(1)%nbr(q, i) = max(g, 0)
                bot_lev(1)%cond(q, i) = merge(rows_glb(2*q + 1, i), 0._wp, g > 0)
                bot_lev(1)%diag(i) = bot_lev(1)%diag(i) + bot_lev(1)%cond(q, i)
            end do
        end do

        do lv = 1, bot_nlev - 1
            bot_lev(lv + 1)%diag = 0._wp
            bot_lev(lv + 1)%cond = 0._wp
            bot_lev(lv + 1)%nbr = 0
            do i = 1, bot_lev(lv)%ncells
                parent_i = bot_lev(lv)%parent(i)
                bot_lev(lv + 1)%diag(parent_i) = bot_lev(lv + 1)%diag(parent_i) + bot_lev(lv)%diag(i)
                do q = 1, 6
                    if (bot_lev(lv)%nbr(q, i) == 0) cycle
                    parent_j = bot_lev(lv)%parent(bot_lev(lv)%nbr(q, i))
                    if (parent_j == parent_i) then
                        bot_lev(lv + 1)%diag(parent_i) = bot_lev(lv + 1)%diag(parent_i) - bot_lev(lv)%cond(q, i)
                    else
                        bot_lev(lv + 1)%nbr(q, parent_i) = parent_j
                        bot_lev(lv + 1)%cond(q, parent_i) = bot_lev(lv + 1)%cond(q, parent_i) + bot_lev(lv)%cond(q, i)
                    end if
                end do
            end do
        end do

        if (bot_cg) return
        n = bot_lev(bot_nlev)%ncells
        bot_chol = 0._wp
        do i = 1, n
            bot_chol(i, i) = bot_lev(bot_nlev)%diag(i)
            do q = 1, 6
                j = bot_lev(bot_nlev)%nbr(q, i)
                if (j > 0) bot_chol(i, j) = bot_chol(i, j) - bot_lev(bot_nlev)%cond(q, i)
            end do
        end do

        ! In-place Cholesky, lower triangle
        do col = 1, n
            bot_chol(col, col) = sqrt(bot_chol(col, col) - sum(bot_chol(col,1:col - 1)**2))
            do i = col + 1, n
                bot_chol(i, col) = (bot_chol(i, col) - sum(bot_chol(i,1:col - 1)*bot_chol(col,1:col - 1)))/bot_chol(col, col)
            end do
        end do

    end subroutine s_mg_bottom_build

    !> Bottom solve: gather the ranks' right-hand sides (v, this rank's cells, x fastest), run the bottom hierarchy's V-cycle on
    !! them, and return this rank's block of the solution in v
    impure subroutine s_mg_bottom_solve(v)

        real(wp), dimension(bot_rank_ncells(proc_rank + 1)), intent(inout) :: v
        integer                                                            :: c, ierr

        c = size(v)
        bot_lev(1)%rhs(1:c) = v
#ifdef MFC_MPI
        call MPI_ALLGATHERV(v, c, mpi_p, bot_lev(1)%rhs, bot_rank_ncells, bot_rank_disp, mpi_p, MPI_COMM_WORLD, ierr)
#endif
        if (bot_cg) then
            call s_bot_cg()
        else
            call s_bot_cycle(1)
        end if
        v = bot_lev(1)%x(bot_rank_disp(proc_rank + 1) + 1:bot_rank_disp(proc_rank + 1) + c)

    end subroutine s_mg_bottom_solve

    !> One V-cycle on bottom level lv, its rhs into its x: proj_mg_sweeps forward Gauss-Seidel sweeps, the coarser levels, the same
    !! sweeps backward; the last level is solved densely. It stays linear (no K-cycle), as a bottom solve that varies with its
    !! right-hand side slows the outer flexible CG markedly. Every rank runs it on the same data, so all agree
    recursive subroutine s_bot_cycle(lv)

        integer, intent(in) :: lv
        integer             :: i, s, n, parent_i

        n = bot_lev(lv)%ncells
        if (lv == bot_nlev) then
            do i = 1, n
                bot_lev(lv)%x(i) = (bot_lev(lv)%rhs(i) - sum(bot_chol(i,1:i - 1)*bot_lev(lv)%x(1:i - 1)))/bot_chol(i, i)
            end do
            do i = n, 1, -1
                bot_lev(lv)%x(i) = (bot_lev(lv)%x(i) - sum(bot_chol(i + 1:n,i)*bot_lev(lv)%x(i + 1:n)))/bot_chol(i, i)
            end do
            return
        end if
        bot_lev(lv)%x = 0._wp
        bot_lev(lv + 1)%rhs = 0._wp
        do s = 1, proj_mg_sweeps
            call s_bot_gs(lv, 1)
        end do
        do i = 1, n
            parent_i = bot_lev(lv)%parent(i)
            bot_lev(lv + 1)%rhs(parent_i) = bot_lev(lv + 1)%rhs(parent_i) + bot_lev(lv)%rhs(i) - f_bot_ax(lv, i)
        end do
        call s_bot_cycle(lv + 1)
        do i = 1, n
            bot_lev(lv)%x(i) = bot_lev(lv)%x(i) + bot_lev(lv + 1)%x(bot_lev(lv)%parent(i))
        end do
        do s = 1, proj_mg_sweeps
            call s_bot_gs(lv, -1)
        end do

    end subroutine s_bot_cycle

    !> Exact bottom solve past mg_bottom_max cells, too many to factor densely: Jacobi-preconditioned CG on the level-1 stencil,
    !! its rhs into its x, from zero to mg_crs_tol relative to the rhs. Every rank runs it on the same data, so all agree
    subroutine s_bot_cg()

        real(wp), dimension(bot_lev(1)%ncells) :: r, z, d, q
        real(wp)                               :: rz, rz_new, alpha, bnorm
        integer                                :: it, i, f

        bot_lev(1)%x = 0._wp
        r = bot_lev(1)%rhs
        bnorm = sqrt(sum(r*r))
        if (bnorm == 0._wp) return
        z = r/bot_lev(1)%diag
        d = z
        rz = sum(r*z)
        do it = 1, bot_lev(1)%ncells
            do i = 1, bot_lev(1)%ncells
                q(i) = bot_lev(1)%diag(i)*d(i)
                do f = 1, 6
                    if (bot_lev(1)%nbr(f, i) > 0) q(i) = q(i) - bot_lev(1)%cond(f, i)*d(bot_lev(1)%nbr(f, i))
                end do
            end do
            alpha = rz/sum(d*q)
            bot_lev(1)%x = bot_lev(1)%x + alpha*d
            r = r - alpha*q
            if (sqrt(sum(r*r)) <= mg_crs_tol*bnorm) exit
            z = r/bot_lev(1)%diag
            rz_new = sum(r*z)
            d = z + (rz_new/rz)*d
            rz = rz_new
        end do

    end subroutine s_bot_cg

    !> (A x)_i of cell i of bottom level lv
    pure function f_bot_ax(lv, i) result(ax)

        integer, intent(in) :: lv, i
        real(wp)            :: ax
        integer             :: q

        ax = bot_lev(lv)%diag(i)*bot_lev(lv)%x(i)
        do q = 1, 6
            if (bot_lev(lv)%nbr(q, i) > 0) ax = ax - bot_lev(lv)%cond(q, i)*bot_lev(lv)%x(bot_lev(lv)%nbr(q, i))
        end do

    end function f_bot_ax

    !> One Gauss-Seidel sweep of x on bottom level lv, in cell order (dir = 1) or reversed (dir = -1)
    subroutine s_bot_gs(lv, dir)

        integer, intent(in) :: lv, dir
        integer             :: i, q, i0, i1
        real(wp)            :: s

        i0 = 1
        i1 = bot_lev(lv)%ncells
        if (dir < 0) then
            i0 = i1
            i1 = 1
        end if
        do i = i0, i1, dir
            s = bot_lev(lv)%rhs(i)
            do q = 1, 6
                if (bot_lev(lv)%nbr(q, i) > 0) s = s + bot_lev(lv)%cond(q, i)*bot_lev(lv)%x(bot_lev(lv)%nbr(q, i))
            end do
            bot_lev(lv)%x(i) = s/bot_lev(lv)%diag(i)
        end do

    end subroutine s_bot_gs

    impure subroutine s_finalize_multigrid_bottom()

        integer :: lv

        if (allocated(bot_chol)) deallocate (bot_chol)
        do lv = 1, bot_nlev
            deallocate (bot_lev(lv)%diag, bot_lev(lv)%cond, bot_lev(lv)%nbr, bot_lev(lv)%parent, bot_lev(lv)%x, bot_lev(lv)%rhs)
        end do
        deallocate (bot_rank_ncells, bot_rank_disp, bot_rank_dims)

    end subroutine s_finalize_multigrid_bottom

end module m_multigrid_bottom
