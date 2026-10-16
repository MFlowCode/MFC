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

    !> Bottom level, gathered to every rank: bot_ncells cells in all, rank r's bot_rank_ncells(r + 1) of them (a bot_rank_dims(:, r
    !! + 1) block) numbered from bot_rank_disp(r + 1) + 1, x fastest. bot_ncells is at most mg_bottom_max, or num_procs past that,
    !! when it is the rank grid itself. Every rank continues the hierarchy on it: level lv holds bot_lev_ncells(lv) cells from
    !! bot_lev_off(lv) + 1, at rank-grid coordinates bot_coords, as a stencil (diagonal bot_diag, conductance bot_cond(q, i) to cell
    !! bot_nbr(q, i), 0 for none); a 2x2x2 box aggregates to cell bot_parent(i) of the next level, and the last level (at most
    !! mg_bottom_max cells) is factored densely in bot_chol. With one level the bottom solve is exact; past mg_bottom_max cells it
    !! is one V-cycle (an exact solve costs every rank more the more ranks there are), unless bot_exact keeps the bottom whole and
    !! CG solves it (bot_cg)
    integer                               :: bot_ncells, bot_nlev, bot_ncells_all
    logical                               :: bot_exact  !< bottom solved exactly (proj_mg_bottom, or by rank count)
    logical                               :: bot_cg     !< exactly, by CG, as it is too large to factor densely
    integer, dimension(mg_maxlev)         :: bot_lev_ncells, bot_lev_off
    integer, dimension(3, mg_maxlev)      :: bot_lev_dims
    integer, allocatable, dimension(:)    :: bot_rank_ncells, bot_rank_disp, bot_parent
    integer, allocatable, dimension(:,:)  :: bot_rank_dims, bot_nbr, bot_coords
    real(wp), allocatable, dimension(:,:) :: bot_chol, bot_cond
    real(wp), allocatable, dimension(:)   :: bot_diag, bot_x, bot_rhs

contains

    !> Every rank's bottom-level block size, for the global numbering of the bottom cells, and the bottom hierarchy's levels; exact
    !! keeps the bottom whole (bot_exact)
    impure subroutine s_initialize_multigrid_bottom(exact)

        logical, intent(in)   :: exact
        integer, dimension(3) :: coords
        integer               :: r, lv, i, ierr

        bot_exact = exact
        allocate (bot_rank_ncells(num_procs), bot_rank_disp(num_procs), bot_rank_dims(3, num_procs))
        bot_rank_dims(:,proc_rank + 1) = [mg_nx(mg_num_lev), mg_ny(mg_num_lev), mg_nz(mg_num_lev)]
#ifdef MFC_MPI
        call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, bot_rank_dims, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
        bot_ncells = 0
        do r = 1, num_procs
            bot_rank_ncells(r) = product(bot_rank_dims(:, &
                            & r)); bot_rank_disp(r) = bot_ncells; bot_ncells = bot_ncells + bot_rank_ncells(r)
        end do

        bot_nlev = 1; bot_lev_ncells(1) = bot_ncells; bot_lev_off(1) = 0; bot_ncells_all = bot_ncells
        bot_lev_dims(:,1) = [num_procs_x, num_procs_y, num_procs_z]
        do while (bot_lev_ncells(bot_nlev) > mg_bottom_max .and. .not. bot_exact)
            bot_lev_dims(:,bot_nlev + 1) = (bot_lev_dims(:,bot_nlev) + 1)/2
            bot_lev_off(bot_nlev + 1) = bot_ncells_all; bot_lev_ncells(bot_nlev + 1) = product(bot_lev_dims(:,bot_nlev + 1))
            bot_nlev = bot_nlev + 1; bot_ncells_all = bot_ncells_all + bot_lev_ncells(bot_nlev)
        end do
        allocate (bot_diag(bot_ncells_all), bot_cond(6, bot_ncells_all), bot_nbr(6, bot_ncells_all), bot_parent(bot_ncells_all), &
                  & bot_coords(3, bot_ncells_all))
        bot_cg = bot_lev_ncells(bot_nlev) > mg_bottom_max
        allocate (bot_x(bot_ncells_all), bot_rhs(bot_ncells_all))
        if (.not. bot_cg) allocate (bot_chol(bot_lev_ncells(bot_nlev), bot_lev_ncells(bot_nlev)))
        if (bot_nlev == 1) return

        if (any(bot_rank_ncells /= 1) .or. product(bot_lev_dims(:,1)) /= bot_ncells) then
            call s_mpi_abort('multigrid: a bottom level past mg_bottom_max cells must hold one cell per rank')
        end if
        coords = 0; coords(1:num_dims) = proc_coords(1:num_dims)
        bot_coords(:,proc_rank + 1) = coords
#ifdef MFC_MPI
        call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, bot_coords, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
        do lv = 1, bot_nlev - 1
            do i = bot_lev_off(lv) + 1, bot_lev_off(lv) + bot_lev_ncells(lv)
                coords = bot_coords(:,i)/2
                bot_parent(i) = bot_lev_off(lv + 1) + 1 + coords(1) + bot_lev_dims(1, lv + 1)*(coords(2) + bot_lev_dims(2, &
                           & lv + 1)*coords(3))
                bot_coords(:,bot_parent(i)) = coords
            end do
        end do

    end subroutine s_initialize_multigrid_bottom

    !> Global bottom index of the neighbor of this rank's bottom cell (i, j, k) across side q, or 0 where a boundary has no seam
    integer function f_bot_nbr(q, i, j, k) result(g)

        integer, intent(in)   :: q, i, j, k
        integer, dimension(3) :: cell
        integer               :: dir, r

        cell = [i, j, k]; dir = (q + 1)/2; r = proc_rank; g = 0
        if (mod(q, 2) == 1) then
            cell(dir) = cell(dir) - 1
            if (cell(dir) < 0) then
                r = mg_nbr_rank(q); if (r < 0) return
                cell(dir) = bot_rank_dims(dir, r + 1) - 1
            end if
        else
            cell(dir) = cell(dir) + 1
            if (cell(dir) >= bot_rank_dims(dir, proc_rank + 1)) then
                r = mg_nbr_rank(q); if (r < 0) return
                cell(dir) = 0
            end if
        end if
        g = bot_rank_disp(r + 1) + (cell(3)*bot_rank_dims(2, r + 1) + cell(2))*bot_rank_dims(1, r + 1) + cell(1) + 1

    end function f_bot_nbr

    !> Gather every bottom cell's row (cell_rows(:, c): diagonal sum and conductances of the x-low, x-high, y-low, y-high, z-low and
    !! z-high faces of this rank's cell c, x fastest; zero for an inactive direction) with the faces' global neighbors into the
    !! global bottom stencil, coarsen it as the device levels are, and factor the last level; a face whose two sides fall in one
    !! cell (a periodic seam across a one-cell width, or a face inside an aggregated box) couples the cell to itself and cancels
    impure subroutine s_mg_bottom_build(cell_rows)

        real(wp), dimension(7, bot_rank_ncells(proc_rank + 1)), intent(in) :: cell_rows
        real(wp), dimension(13, bot_rank_ncells(proc_rank + 1)) :: rows_loc
        real(wp), dimension(13, bot_ncells) :: rows_glb
        integer :: i, j, k, c, q, g, col, ierr, lv, lev_off, lev_n, parent_i, parent_j

        c = 0
        do k = 0, bot_rank_dims(3, proc_rank + 1) - 1
            do j = 0, bot_rank_dims(2, proc_rank + 1) - 1
                do i = 0, bot_rank_dims(1, proc_rank + 1) - 1
                    c = c + 1
                    rows_loc(1, c) = cell_rows(1, c); rows_loc(3:13:2,c) = cell_rows(2:7,c); rows_loc(2:12:2,c) = 0._wp
                    do q = 1, 2*num_dims
                        rows_loc(2*q, c) = real(f_bot_nbr(q, i, j, k), wp)
                    end do
                end do
            end do
        end do

        rows_glb(:,1:c) = rows_loc
#ifdef MFC_MPI
        call MPI_ALLGATHERV(rows_loc, 13*c, mpi_p, rows_glb, 13*bot_rank_ncells, 13*bot_rank_disp, mpi_p, MPI_COMM_WORLD, ierr)
#endif

        do i = 1, bot_ncells
            bot_diag(i) = rows_glb(1, i)
            do q = 1, 6
                g = nint(rows_glb(2*q, i))
                ! No neighbour: a domain face, whose conductance (a pressure outlet's; zero on a wall) couples to a fixed ghost
                if (g == 0) bot_diag(i) = bot_diag(i) + rows_glb(2*q + 1, i)
                if (g == i) g = 0
                bot_nbr(q, i) = max(g, 0)
                bot_cond(q, i) = merge(rows_glb(2*q + 1, i), 0._wp, g > 0)
                bot_diag(i) = bot_diag(i) + bot_cond(q, i)
            end do
        end do

        do lv = 1, bot_nlev - 1
            lev_off = bot_lev_off(lv + 1); lev_n = bot_lev_ncells(lv + 1)
            bot_diag(lev_off + 1:lev_off + lev_n) = 0._wp; bot_cond(:,lev_off + 1:lev_off + lev_n) = 0._wp; bot_nbr(:, &
                     & lev_off + 1:lev_off + lev_n) = 0
            do i = bot_lev_off(lv) + 1, bot_lev_off(lv) + bot_lev_ncells(lv)
                parent_i = bot_parent(i); bot_diag(parent_i) = bot_diag(parent_i) + bot_diag(i)
                do q = 1, 6
                    if (bot_nbr(q, i) == 0) cycle
                    parent_j = bot_parent(bot_nbr(q, i))
                    if (parent_j == parent_i) then
                        bot_diag(parent_i) = bot_diag(parent_i) - bot_cond(q, i)
                    else
                        bot_nbr(q, parent_i) = parent_j; bot_cond(q, parent_i) = bot_cond(q, parent_i) + bot_cond(q, i)
                    end if
                end do
            end do
        end do

        if (bot_cg) return
        lev_off = bot_lev_off(bot_nlev); lev_n = bot_lev_ncells(bot_nlev)
        bot_chol = 0._wp
        do i = 1, lev_n
            bot_chol(i, i) = bot_diag(lev_off + i)
            do q = 1, 6
                j = bot_nbr(q, lev_off + i)
                if (j > 0) bot_chol(i, j - lev_off) = bot_chol(i, j - lev_off) - bot_cond(q, lev_off + i)
            end do
        end do

        ! In-place Cholesky, lower triangle
        do col = 1, lev_n
            bot_chol(col, col) = sqrt(bot_chol(col, col) - sum(bot_chol(col,1:col - 1)**2))
            do i = col + 1, lev_n
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
        bot_rhs(1:c) = v
#ifdef MFC_MPI
        call MPI_ALLGATHERV(v, c, mpi_p, bot_rhs, bot_rank_ncells, bot_rank_disp, mpi_p, MPI_COMM_WORLD, ierr)
#endif
        if (bot_cg) then
            call s_bot_cg()
        else
            call s_bot_cycle(1)
        end if
        v = bot_x(bot_rank_disp(proc_rank + 1) + 1:bot_rank_disp(proc_rank + 1) + c)

    end subroutine s_mg_bottom_solve

    !> One V-cycle on bottom level lv, bot_rhs into bot_x: proj_mg_sweeps forward Gauss-Seidel sweeps, the coarser levels, the same
    !! sweeps backward; the last level is solved densely. It stays linear (no K-cycle), as a bottom solve that varies with its
    !! right-hand side slows the outer flexible CG markedly. Every rank runs it on the same data, so all agree
    recursive subroutine s_bot_cycle(lv)

        integer, intent(in) :: lv
        integer             :: i, s, lev_off, lev_n

        lev_off = bot_lev_off(lv); lev_n = bot_lev_ncells(lv)
        if (lv == bot_nlev) then
            do i = 1, lev_n
                bot_x(lev_off + i) = (bot_rhs(lev_off + i) - sum(bot_chol(i, &
                      & 1:i - 1)*bot_x(lev_off + 1:lev_off + i - 1)))/bot_chol(i, i)
            end do
            do i = lev_n, 1, -1
                bot_x(lev_off + i) = (bot_x(lev_off + i) - sum(bot_chol(i + 1:lev_n, &
                      & i)*bot_x(lev_off + i + 1:lev_off + lev_n)))/bot_chol(i, i)
            end do
            return
        end if
        bot_x(lev_off + 1:lev_off + lev_n) = 0._wp
        bot_rhs(bot_lev_off(lv + 1) + 1:bot_lev_off(lv + 1) + bot_lev_ncells(lv + 1)) = 0._wp
        do s = 1, proj_mg_sweeps
            call s_bot_gs(lv, 1)
        end do
        do i = lev_off + 1, lev_off + lev_n
            bot_rhs(bot_parent(i)) = bot_rhs(bot_parent(i)) + bot_rhs(i) - f_bot_ax(i)
        end do
        call s_bot_cycle(lv + 1)
        do i = lev_off + 1, lev_off + lev_n
            bot_x(i) = bot_x(i) + bot_x(bot_parent(i))
        end do
        do s = 1, proj_mg_sweeps
            call s_bot_gs(lv, -1)
        end do

    end subroutine s_bot_cycle

    !> Exact bottom solve past mg_bottom_max cells, too many to factor densely: Jacobi-preconditioned CG on the level-1 stencil,
    !! bot_rhs into bot_x, from zero to mg_crs_tol relative to bot_rhs. Every rank runs it on the same data, so all agree
    subroutine s_bot_cg()

        real(wp), dimension(bot_ncells) :: r, z, d, q
        real(wp)                        :: rz, rz_new, alpha, bnorm
        integer                         :: it, i, f

        bot_x(1:bot_ncells) = 0._wp; r = bot_rhs(1:bot_ncells)
        bnorm = sqrt(sum(r*r))
        if (bnorm == 0._wp) return
        z = r/bot_diag(1:bot_ncells); d = z; rz = sum(r*z)
        do it = 1, bot_ncells
            do i = 1, bot_ncells
                q(i) = bot_diag(i)*d(i)
                do f = 1, 6
                    if (bot_nbr(f, i) > 0) q(i) = q(i) - bot_cond(f, i)*d(bot_nbr(f, i))
                end do
            end do
            alpha = rz/sum(d*q)
            bot_x(1:bot_ncells) = bot_x(1:bot_ncells) + alpha*d; r = r - alpha*q
            if (sqrt(sum(r*r)) <= mg_crs_tol*bnorm) exit
            z = r/bot_diag(1:bot_ncells); rz_new = sum(r*z)
            d = z + (rz_new/rz)*d; rz = rz_new
        end do

    end subroutine s_bot_cg

    !> (A x)_i of bottom-hierarchy cell i
    pure real(wp) function f_bot_ax(i) result(ax)
        integer, intent(in) :: i
        integer             :: q

        ax = bot_diag(i)*bot_x(i)
        do q = 1, 6
            if (bot_nbr(q, i) > 0) ax = ax - bot_cond(q, i)*bot_x(bot_nbr(q, i))
        end do

    end function f_bot_ax

    !> One Gauss-Seidel sweep of bot_x on bottom level lv, in cell order (dir = 1) or reversed (dir = -1)
    subroutine s_bot_gs(lv, dir)

        integer, intent(in) :: lv, dir
        integer             :: i, q, i0, i1
        real(wp)            :: s

        i0 = bot_lev_off(lv) + 1; i1 = bot_lev_off(lv) + bot_lev_ncells(lv)
        if (dir < 0) then
            i0 = i1; i1 = bot_lev_off(lv) + 1
        end if
        do i = i0, i1, dir
            s = bot_rhs(i)
            do q = 1, 6
                if (bot_nbr(q, i) > 0) s = s + bot_cond(q, i)*bot_x(bot_nbr(q, i))
            end do
            bot_x(i) = s/bot_diag(i)
        end do

    end subroutine s_bot_gs

    impure subroutine s_finalize_multigrid_bottom()

        if (allocated(bot_chol)) deallocate (bot_chol)
        deallocate (bot_x, bot_rhs, bot_diag, bot_cond, bot_nbr, bot_parent, bot_coords, bot_rank_ncells, bot_rank_disp, &
                    & bot_rank_dims)

    end subroutine s_finalize_multigrid_bottom

end module m_multigrid_bottom
