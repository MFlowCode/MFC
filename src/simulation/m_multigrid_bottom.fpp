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

    !> Bottom level, gathered to every rank: crs_n cells in all, rank r's crs_cnt(r + 1) of them (a crs_sz(:, r + 1) block) numbered
    !! from crs_disp(r + 1) + 1, x fastest. crs_n is at most mg_bottom_max, or num_procs past that, when it is the rank grid itself.
    !! Every rank continues the hierarchy on it: level lv holds crs_ln(lv) cells from crs_loff(lv) + 1, at rank-grid coordinates
    !! crs_co, as a stencil (diagonal crs_ad, conductance crs_ak(q, i) to cell crs_aj(q, i), 0 for none); a 2x2x2 box aggregates to
    !! cell crs_agg(i) of the next level, and the last level (at most mg_bottom_max cells) is factored densely in crs_l. With one
    !! level the bottom solve is exact; past mg_bottom_max cells it is one V-cycle (an exact solve costs every rank more the more
    !! ranks there are), unless crs_exact keeps the bottom whole and CG solves it (crs_cg)
    integer                               :: crs_n, crs_nl, crs_tot
    logical                               :: crs_exact  !< bottom solved exactly (proj_mg_bottom, or by rank count)
    logical                               :: crs_cg     !< exactly, by CG, as it is too large to factor densely
    integer, dimension(mg_maxlev)         :: crs_ln, crs_loff
    integer, dimension(3, mg_maxlev)      :: crs_ld
    integer, allocatable, dimension(:)    :: crs_cnt, crs_disp, crs_agg
    integer, allocatable, dimension(:,:)  :: crs_sz, crs_aj, crs_co
    real(wp), allocatable, dimension(:,:) :: crs_l, crs_ak
    real(wp), allocatable, dimension(:)   :: crs_ad, crs_x, crs_b

contains

    !> Every rank's bottom-level block size, for the global numbering of the bottom cells, and the bottom hierarchy's levels; exact
    !! keeps the bottom whole (crs_exact)
    impure subroutine s_initialize_multigrid_bottom(exact)

        logical, intent(in)   :: exact
        integer, dimension(3) :: pc
        integer               :: r, lv, i, ierr

        crs_exact = exact
        allocate (crs_cnt(num_procs), crs_disp(num_procs), crs_sz(3, num_procs))
        crs_sz(:,proc_rank + 1) = [mg_nx(mg_nlev), mg_ny(mg_nlev), mg_nz(mg_nlev)]
#ifdef MFC_MPI
        call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, crs_sz, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
        crs_n = 0
        do r = 1, num_procs
            crs_cnt(r) = product(crs_sz(:,r)); crs_disp(r) = crs_n; crs_n = crs_n + crs_cnt(r)
        end do

        crs_nl = 1; crs_ln(1) = crs_n; crs_loff(1) = 0; crs_tot = crs_n
        crs_ld(:,1) = [num_procs_x, num_procs_y, num_procs_z]
        do while (crs_ln(crs_nl) > mg_bottom_max .and. .not. crs_exact)
            crs_ld(:,crs_nl + 1) = (crs_ld(:,crs_nl) + 1)/2
            crs_loff(crs_nl + 1) = crs_tot; crs_ln(crs_nl + 1) = product(crs_ld(:,crs_nl + 1))
            crs_nl = crs_nl + 1; crs_tot = crs_tot + crs_ln(crs_nl)
        end do
        allocate (crs_ad(crs_tot), crs_ak(6, crs_tot), crs_aj(6, crs_tot), crs_agg(crs_tot), crs_co(3, crs_tot))
        crs_cg = crs_ln(crs_nl) > mg_bottom_max
        allocate (crs_x(crs_tot), crs_b(crs_tot))
        if (.not. crs_cg) allocate (crs_l(crs_ln(crs_nl), crs_ln(crs_nl)))
        if (crs_nl == 1) return

        if (any(crs_cnt /= 1) .or. product(crs_ld(:,1)) /= crs_n) then
            call s_mpi_abort('multigrid: a bottom level past mg_bottom_max cells must hold one cell per rank')
        end if
        pc = 0; pc(1:num_dims) = proc_coords(1:num_dims)
        crs_co(:,proc_rank + 1) = pc
#ifdef MFC_MPI
        call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, crs_co, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
        do lv = 1, crs_nl - 1
            do i = crs_loff(lv) + 1, crs_loff(lv) + crs_ln(lv)
                pc = crs_co(:,i)/2
                crs_agg(i) = crs_loff(lv + 1) + 1 + pc(1) + crs_ld(1, lv + 1)*(pc(2) + crs_ld(2, lv + 1)*pc(3))
                crs_co(:,crs_agg(i)) = pc
            end do
        end do

    end subroutine s_initialize_multigrid_bottom

    !> Global bottom index of the neighbor of this rank's bottom cell (i, j, k) across side q, or 0 where a boundary has no seam
    integer function f_crs_nbr(q, i, j, k) result(g)

        integer, intent(in)   :: q, i, j, k
        integer, dimension(3) :: ic
        integer               :: nd, r

        ic = [i, j, k]; nd = (q + 1)/2; r = proc_rank; g = 0
        if (mod(q, 2) == 1) then
            ic(nd) = ic(nd) - 1
            if (ic(nd) < 0) then
                r = mg_nbr(q); if (r < 0) return
                ic(nd) = crs_sz(nd, r + 1) - 1
            end if
        else
            ic(nd) = ic(nd) + 1
            if (ic(nd) >= crs_sz(nd, proc_rank + 1)) then
                r = mg_nbr(q); if (r < 0) return
                ic(nd) = 0
            end if
        end if
        g = crs_disp(r + 1) + (ic(3)*crs_sz(2, r + 1) + ic(2))*crs_sz(1, r + 1) + ic(1) + 1

    end function f_crs_nbr

    !> Gather every bottom cell's row (dk(:, c): diagonal sum and conductances of the x-low, x-high, y-low, y-high, z-low and z-high
    !! faces of this rank's cell c, x fastest; zero for an inactive direction) with the faces' global neighbors into the global
    !! bottom stencil, coarsen it as the device levels are, and factor the last level; a face whose two sides fall in one cell (a
    !! periodic seam across a one-cell width, or a face inside an aggregated box) couples the cell to itself and cancels
    impure subroutine s_mg_bottom_build(dk)

        real(wp), dimension(7, crs_cnt(proc_rank + 1)), intent(in) :: dk
        real(wp), dimension(13, crs_cnt(proc_rank + 1))            :: row
        real(wp), dimension(13, crs_n)                             :: rows
        integer                                                    :: i, j, k, c, q, g, jr, ierr, lv, o, nn, ci, cj

        c = 0
        do k = 0, crs_sz(3, proc_rank + 1) - 1
            do j = 0, crs_sz(2, proc_rank + 1) - 1
                do i = 0, crs_sz(1, proc_rank + 1) - 1
                    c = c + 1
                    row(1, c) = dk(1, c); row(3:13:2,c) = dk(2:7,c); row(2:12:2,c) = 0._wp
                    do q = 1, 2*num_dims
                        row(2*q, c) = real(f_crs_nbr(q, i, j, k), wp)
                    end do
                end do
            end do
        end do

        rows(:,1:c) = row
#ifdef MFC_MPI
        call MPI_ALLGATHERV(row, 13*c, mpi_p, rows, 13*crs_cnt, 13*crs_disp, mpi_p, MPI_COMM_WORLD, ierr)
#endif

        do i = 1, crs_n
            crs_ad(i) = rows(1, i)
            do q = 1, 6
                g = nint(rows(2*q, i))
                ! No neighbour: a domain face, whose conductance (a pressure outlet's; zero on a wall) couples to a fixed ghost
                if (g == 0) crs_ad(i) = crs_ad(i) + rows(2*q + 1, i)
                if (g == i) g = 0
                crs_aj(q, i) = max(g, 0)
                crs_ak(q, i) = merge(rows(2*q + 1, i), 0._wp, g > 0)
                crs_ad(i) = crs_ad(i) + crs_ak(q, i)
            end do
        end do

        do lv = 1, crs_nl - 1
            o = crs_loff(lv + 1); nn = crs_ln(lv + 1)
            crs_ad(o + 1:o + nn) = 0._wp; crs_ak(:,o + 1:o + nn) = 0._wp; crs_aj(:,o + 1:o + nn) = 0
            do i = crs_loff(lv) + 1, crs_loff(lv) + crs_ln(lv)
                ci = crs_agg(i); crs_ad(ci) = crs_ad(ci) + crs_ad(i)
                do q = 1, 6
                    if (crs_aj(q, i) == 0) cycle
                    cj = crs_agg(crs_aj(q, i))
                    if (cj == ci) then
                        crs_ad(ci) = crs_ad(ci) - crs_ak(q, i)
                    else
                        crs_aj(q, ci) = cj; crs_ak(q, ci) = crs_ak(q, ci) + crs_ak(q, i)
                    end if
                end do
            end do
        end do

        if (crs_cg) return
        o = crs_loff(crs_nl); nn = crs_ln(crs_nl)
        crs_l = 0._wp
        do i = 1, nn
            crs_l(i, i) = crs_ad(o + i)
            do q = 1, 6
                j = crs_aj(q, o + i)
                if (j > 0) crs_l(i, j - o) = crs_l(i, j - o) - crs_ak(q, o + i)
            end do
        end do

        ! In-place Cholesky, lower triangle
        do jr = 1, nn
            crs_l(jr, jr) = sqrt(crs_l(jr, jr) - sum(crs_l(jr,1:jr - 1)**2))
            do i = jr + 1, nn
                crs_l(i, jr) = (crs_l(i, jr) - sum(crs_l(i,1:jr - 1)*crs_l(jr,1:jr - 1)))/crs_l(jr, jr)
            end do
        end do

    end subroutine s_mg_bottom_build

    !> Bottom solve: gather the ranks' right-hand sides (v, this rank's cells, x fastest), run the bottom hierarchy's V-cycle on
    !! them, and return this rank's block of the solution in v
    impure subroutine s_mg_bottom_solve(v)

        real(wp), dimension(crs_cnt(proc_rank + 1)), intent(inout) :: v
        integer                                                    :: c, ierr

        c = size(v)
        crs_b(1:c) = v
#ifdef MFC_MPI
        call MPI_ALLGATHERV(v, c, mpi_p, crs_b, crs_cnt, crs_disp, mpi_p, MPI_COMM_WORLD, ierr)
#endif
        if (crs_cg) then
            call s_crs_cg()
        else
            call s_crs_cycle(1)
        end if
        v = crs_x(crs_disp(proc_rank + 1) + 1:crs_disp(proc_rank + 1) + c)

    end subroutine s_mg_bottom_solve

    !> One V-cycle on bottom level lv, crs_b into crs_x: proj_mg_sweeps forward Gauss-Seidel sweeps, the coarser levels, the same
    !! sweeps backward; the last level is solved densely. It stays linear (no K-cycle), as a bottom solve that varies with its
    !! right-hand side slows the outer flexible CG markedly. Every rank runs it on the same data, so all agree
    recursive subroutine s_crs_cycle(lv)

        integer, intent(in) :: lv
        integer             :: i, s, o, nn

        o = crs_loff(lv); nn = crs_ln(lv)
        if (lv == crs_nl) then
            do i = 1, nn
                crs_x(o + i) = (crs_b(o + i) - sum(crs_l(i,1:i - 1)*crs_x(o + 1:o + i - 1)))/crs_l(i, i)
            end do
            do i = nn, 1, -1
                crs_x(o + i) = (crs_x(o + i) - sum(crs_l(i + 1:nn,i)*crs_x(o + i + 1:o + nn)))/crs_l(i, i)
            end do
            return
        end if
        crs_x(o + 1:o + nn) = 0._wp
        crs_b(crs_loff(lv + 1) + 1:crs_loff(lv + 1) + crs_ln(lv + 1)) = 0._wp
        do s = 1, proj_mg_sweeps
            call s_crs_gs(lv, 1)
        end do
        do i = o + 1, o + nn
            crs_b(crs_agg(i)) = crs_b(crs_agg(i)) + crs_b(i) - f_crs_ax(i)
        end do
        call s_crs_cycle(lv + 1)
        do i = o + 1, o + nn
            crs_x(i) = crs_x(i) + crs_x(crs_agg(i))
        end do
        do s = 1, proj_mg_sweeps
            call s_crs_gs(lv, -1)
        end do

    end subroutine s_crs_cycle

    !> Exact bottom solve past mg_bottom_max cells, too many to factor densely: Jacobi-preconditioned CG on the level-1 stencil,
    !! crs_b into crs_x, from zero to mg_crs_tol relative to crs_b. Every rank runs it on the same data, so all agree
    subroutine s_crs_cg()

        real(wp), dimension(crs_n) :: r, z, d, q
        real(wp)                   :: rz, rz_new, alpha, bnorm
        integer                    :: it, i, f

        crs_x(1:crs_n) = 0._wp; r = crs_b(1:crs_n)
        bnorm = sqrt(sum(r*r))
        if (bnorm == 0._wp) return
        z = r/crs_ad(1:crs_n); d = z; rz = sum(r*z)
        do it = 1, crs_n
            do i = 1, crs_n
                q(i) = crs_ad(i)*d(i)
                do f = 1, 6
                    if (crs_aj(f, i) > 0) q(i) = q(i) - crs_ak(f, i)*d(crs_aj(f, i))
                end do
            end do
            alpha = rz/sum(d*q)
            crs_x(1:crs_n) = crs_x(1:crs_n) + alpha*d; r = r - alpha*q
            if (sqrt(sum(r*r)) <= mg_crs_tol*bnorm) exit
            z = r/crs_ad(1:crs_n); rz_new = sum(r*z)
            d = z + (rz_new/rz)*d; rz = rz_new
        end do

    end subroutine s_crs_cg

    !> (A x)_i of bottom-hierarchy cell i
    pure real(wp) function f_crs_ax(i) result(ax)
        integer, intent(in) :: i
        integer             :: q

        ax = crs_ad(i)*crs_x(i)
        do q = 1, 6
            if (crs_aj(q, i) > 0) ax = ax - crs_ak(q, i)*crs_x(crs_aj(q, i))
        end do

    end function f_crs_ax

    !> One Gauss-Seidel sweep of crs_x on bottom level lv, in cell order (dir = 1) or reversed (dir = -1)
    subroutine s_crs_gs(lv, dir)

        integer, intent(in) :: lv, dir
        integer             :: i, q, i0, i1
        real(wp)            :: s

        i0 = crs_loff(lv) + 1; i1 = crs_loff(lv) + crs_ln(lv)
        if (dir < 0) then
            i0 = i1; i1 = crs_loff(lv) + 1
        end if
        do i = i0, i1, dir
            s = crs_b(i)
            do q = 1, 6
                if (crs_aj(q, i) > 0) s = s + crs_ak(q, i)*crs_x(crs_aj(q, i))
            end do
            crs_x(i) = s/crs_ad(i)
        end do

    end subroutine s_crs_gs

    impure subroutine s_finalize_multigrid_bottom()

        if (allocated(crs_l)) deallocate (crs_l)
        deallocate (crs_x, crs_b, crs_ad, crs_ak, crs_aj, crs_agg, crs_co, crs_cnt, crs_disp, crs_sz)

    end subroutine s_finalize_multigrid_bottom

end module m_multigrid_bottom
