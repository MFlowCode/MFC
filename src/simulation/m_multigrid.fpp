!>
!! @file
!! @brief Contains module m_multigrid

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Geometric multigrid for the projection's pressure system: a symmetric positive-definite seven-point (five, three) operator
!! A = D + sum over faces of a conductance times the difference across it, on this rank's Cartesian block. s_mg_set_operator takes D
!! and the face conductances and builds a Galerkin hierarchy by piecewise-constant aggregation; s_mg_solve runs flexible PCG,
!! preconditioned by one V- or K-cycle. Each level exchanges one ghost layer with its neighbors (one aggregated message per
!! neighbor) and coarsens to one cell per rank; every rank then continues the hierarchy on that gathered rank-level problem: solved
!! directly when small, else by one host V-cycle. The solve parameters are the projection's (proj_tol, proj_mg_*).
module m_multigrid

#ifdef MFC_MPI
    use mpi
#endif
    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_nvtx

    implicit none

    private; public :: s_initialize_multigrid_module, s_mg_set_operator, s_mg_solve, s_finalize_multigrid_module

    integer, parameter :: mg_maxlev = 24
    real(wp), parameter :: res_floor = 1.e2_wp*epsilon(1._wp)  !< residual round-off floor, relative to the right-hand side
    real(wp), allocatable, dimension(:,:,:) :: rs, zs, qs      !< PCG residual, preconditioned residual and operator image
    $:GPU_DECLARE(create='[rs, zs, qs]')

    !> Multigrid hierarchy, flattened, with one ghost layer in each active direction; level lv starts after mg_off(lv), x fastest.
    !! Every level coarsens by two (an odd size folds its last cell into the last coarse cell) down to the bottom level, the first
    !! with at most mg_bottom_max cells in all (or one cell per rank). PCG's search direction follows the levels in mg_e, at
    !! mg_poff, laid out as level 1
    integer, parameter                  :: mg_bottom_max = 128
    integer                             :: mg_nlev, mg_gx, mg_gy, mg_gz, mg_poff
    integer                             :: mg_blev             !< Level this solve stops at (s_mg_omega): mg_nlev, or the truncation
    real(wp), dimension(mg_maxlev)      :: mg_om               !< Coarse-correction scale into each level (s_mg_omega)
    real(wp), parameter                 :: mg_om_r0 = 0.02_wp  !< Helmholtz-to-Laplacian ratio that halves the scale's excess
    integer, dimension(mg_maxlev)       :: mg_nx, mg_ny, mg_nz, mg_off
    real(wp), allocatable, dimension(:) :: mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r
    $:GPU_DECLARE(create='[mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r]')
    !> Flat index of cell (i, j, k) on level lv: mg_o(lv) + k*mg_sz(lv) + j*mg_sy(lv) + i (f_mg_ix)
    integer, dimension(mg_maxlev) :: mg_o, mg_sy, mg_sz
    $:GPU_DECLARE(create='[mg_o, mg_sy, mg_sz]')
    !> K-cycle vectors of the coarse levels (indexed as mg_e; allocated with mg_kfull > 0): the right-hand side (kr), the first
    !! step's preconditioned direction (kc) and its image under the level operator (kv)
    real(wp), allocatable, dimension(:) :: mg_kr, mg_kc, mg_kv
    $:GPU_DECLARE(create='[mg_kr, mg_kc, mg_kv]')
    real(wp), parameter :: mg_kt = 0.25_wp  !< K-cycle skips its second step once the residual falls by this

    !> Ghost-layer exchange: side q = 1..6 is (x, y, z) x (low, high); mg_nbr(q) is the neighbor rank, this rank for a periodic
    !! seam, or -1. Each exchange sends one message per distinct neighbor, holding all its sides
    integer, dimension(6)               :: mg_nbr
    real(wp), allocatable, dimension(:) :: mg_sbuf, mg_rbuf
    $:GPU_DECLARE(create='[mg_sbuf, mg_rbuf]')
    !> Side q's action (0 none, 1 pack for another rank, 3 periodic self copy) and its offset in mg_sbuf (1) and mg_rbuf (2) by
    !! level; message i goes to rank mg_nlist(i), from mg_sbeg + 1 for mg_slen values, and arrives at mg_rbeg + 1 for mg_rlen
    integer                             :: mg_nmsg
    integer, dimension(6)               :: mg_smode, mg_nlist
    integer, dimension(6, mg_maxlev, 2) :: mg_boff
    integer, dimension(6, mg_maxlev)    :: mg_sbeg, mg_slen, mg_rbeg, mg_rlen
    $:GPU_DECLARE(create='[mg_smode, mg_boff]')

    !> Bottom level, gathered to every rank: crs_n cells in all, rank r's crs_cnt(r + 1) of them (a crs_sz(:, r + 1) block) numbered
    !! from crs_disp(r + 1) + 1, x fastest. crs_n is at most mg_bottom_max, or num_procs past that, when it is the rank grid itself.
    !! Every rank continues the hierarchy on it: level lv holds crs_ln(lv) cells from crs_loff(lv) + 1, at rank-grid coordinates
    !! crs_co, as a stencil (diagonal crs_ad, conductance crs_ak(q, i) to cell crs_aj(q, i), 0 for none); a 2x2x2 box aggregates to
    !! cell crs_agg(i) of the next level, and the last level (at most mg_bottom_max cells) is factored densely in crs_l. With one
    !! level the bottom solve is exact; past mg_bottom_max cells it is one V-cycle (an exact solve costs every rank more the more
    !! ranks there are), unless crs_exact keeps the bottom whole and CG solves it (crs_cg)
    integer :: crs_n, crs_nl, crs_tot
    logical :: crs_exact  !< bottom solved exactly (proj_mg_bottom, or by rank count)
    logical :: crs_cg  !< exactly, by CG, as it is too large to factor densely
    real(wp), parameter :: crs_tol = 1.e-12_wp  !< CG bottom tolerance, relative to its right-hand side
    integer :: mg_kfull  !< levels between K-cycle steps when a solve reaches the bottom
    integer :: mg_kspace  !< levels between K-cycle steps in this solve: mg_kfull, or V-cycles once truncated (s_mg_omega)
    integer, dimension(mg_maxlev) :: crs_ln, crs_loff
    integer, dimension(3, mg_maxlev) :: crs_ld
    integer, allocatable, dimension(:) :: crs_cnt, crs_disp, crs_agg
    integer, allocatable, dimension(:,:) :: crs_sz, crs_aj, crs_co
    real(wp), allocatable, dimension(:,:) :: crs_l, crs_ak
    real(wp), allocatable, dimension(:) :: crs_ad, crs_x, crs_b

contains

    impure subroutine s_initialize_multigrid_module()

        integer  :: lv, tot, kl, ku
        real(wp) :: rlev

        ! Automatic choices by rank count; the thresholds are machine-dependent (docs: "Projection method iterative solve tuning")

        mg_kfull = proj_mg_kcycle
        if (mg_kfull < 0) mg_kfull = merge(2, 0, num_procs >= proj_mg_k_ranks)
        mg_kspace = mg_kfull
        ! A K-cycle reaches the bottom several times per cycle, which multiplies the cost of an exact solve there
        crs_exact = proj_mg_bottom == 1 .or. (proj_mg_bottom == -1 .and. mg_kfull == 0 .and. num_procs <= proj_mg_cg_ranks)

        @:ALLOCATE(rs(0:m, 0:n, 0:p), zs(0:m, 0:n, 0:p), qs(0:m, 0:n, 0:p))

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
        kl = mg_off(min(2, mg_nlev)) + 1; ku = merge(mg_poff, kl - 1, mg_kfull > 0)
        @:ALLOCATE(mg_kr(kl:ku), mg_kc(kl:ku), mg_kv(kl:ku))
        tot = 2*((n + 1)*(p + 1) + (m + 1)*(p + 1) + (m + 1)*(n + 1))
        @:ALLOCATE(mg_sbuf(tot), mg_rbuf(tot))
        @:PIN_HOST(mg_sbuf, mg_rbuf)

        #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            mg_nbr(2*${D}$ - 1) = f_seam_rank(bc_${XYZ}$%beg, num_dims >= ${D}$)
            mg_nbr(2*${D}$) = f_seam_rank(bc_${XYZ}$%end, num_dims >= ${D}$)
        #:endfor
        call s_mg_exchange_layout()
        call s_mg_bottom_layout()

    end subroutine s_initialize_multigrid_module

    !> Solve A x = b, A the operator of s_mg_set_operator, by flexible PCG preconditioned by one multigrid cycle per iteration, from
    !! the initial guess in x to proj_tol relative to the initial residual (or a round-off floor relative to b), at most
    !! proj_max_iters iterations. Returns the iterations taken and the final residual relative to the initial one (0 for none)
    impure subroutine s_mg_solve(b, x, iters, res)

        real(wp), dimension(0:m,0:n,0:p), intent(in)    :: b
        real(wp), dimension(0:m,0:n,0:p), intent(inout) :: x
        integer, intent(out)                            :: iters
        real(wp), intent(out)                           :: res
        real(wp)                                        :: bnorm, rnorm, rtol, rz, rz_new, zq, dq, alpha, beta, r0
        integer                                         :: it, j, k, l, idx, gx, gy, gz, sh

        iters = 0; res = 0._wp
        gx = mg_gx; gy = mg_gy; gz = mg_gz; sh = mg_poff

        ! The search direction's wall ghosts stay zero: those faces carry no conductance, and the ghosts must not hold garbage
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -gz, p + gz
            do k = -gy, n + gy
                do j = -gx, m + gx
                    mg_e(f_mg_ix(1, j, k, l) + sh) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_set_dir(x)
        call s_apply_operator()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rs(j, k, l) = b(j, k, l) - qs(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Relative to the initial residual, i.e. to the pressure change being solved for; relative to b the tolerance would
        ! depend on the ambient pressure, which b carries in full
        bnorm = sqrt(f_dot(b, b))
        rnorm = sqrt(f_dot(rs, rs))
        rtol = max(proj_tol*rnorm, res_floor*bnorm)
        r0 = rnorm

        if (rnorm > rtol) then
            call nvtxStartRange("PROJ-VCYCLE")
            call s_mg_precond()
            call nvtxEndRange
            call s_set_dir(zs)
            rz = f_dot(rs, zs)

            do it = 1, proj_max_iters
                iters = iters + 1
                call s_apply_operator()
                dq = f_dot_dir(qs)
                alpha = rz/dq
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, idx]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            idx = f_mg_ix(1, j, k, l) + sh
                            x(j, k, l) = x(j, k, l) + alpha*mg_e(idx)
                            rs(j, k, l) = rs(j, k, l) - alpha*qs(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                rnorm = sqrt(f_dot(rs, rs))
                if (rnorm <= rtol) exit

                call nvtxStartRange("PROJ-VCYCLE")
                call s_mg_precond()
                call nvtxEndRange
                ! Flexible (Polak-Ribiere) beta, z.(r - r_old)/rz with r_old - r = alpha*q, as a K-cycle preconditioner varies
                call s_dot_rz_zq(rz_new, zq)
                beta = -alpha*zq/rz
                rz = rz_new
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, idx]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            idx = f_mg_ix(1, j, k, l) + sh
                            mg_e(idx) = zs(j, k, l) + beta*mg_e(idx)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            ! Inside the branch: a solve already below its floor takes no iterations and reports none
            res = rnorm/r0
        end if

    end subroutine s_mg_solve

    !> Search direction (interior) = v
    subroutine s_set_dir(v)

        real(wp), dimension(0:m,0:n,0:p), intent(in) :: v
        integer                                      :: j, k, l, sh

        sh = mg_poff
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mg_e(f_mg_ix(1, j, k, l) + sh) = v(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_set_dir

    !> qs = A d for the search direction d, held in level 1's layout: A is level 1's operator, whose boundary faces carry no
    !! conductance unless they are rank or periodic seams, so one exchanged ghost layer is all it needs
    impure subroutine s_apply_operator()

        integer  :: j, k, l, idx, sh
        real(wp) :: dg, nb

        call nvtxStartRange("PROJ-DIR-EXCHANGE")
        call s_mg_exchange(1, mg_poff)
        call nvtxEndRange
        sh = mg_poff
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, idx, dg, nb]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    idx = f_mg_ix(1, j, k, l)
                    call s_mg_row(1, idx, sh, dg, nb)
                    qs(j, k, l) = dg*mg_e(idx + sh) - nb
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_operator

    !> Global r.z and q.z over the interior, in one reduction
    impure subroutine s_dot_rz_zq(rz, zq)

        real(wp), intent(out)     :: rz, zq
        real(wp), dimension(2, 1) :: sl, sg
        real(wp)                  :: a, b
        integer                   :: j, k, l

        a = 0._wp; b = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[a, b]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    a = a + rs(j, k, l)*zs(j, k, l)
                    b = b + qs(j, k, l)*zs(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        sl(:,1) = [a, b]
        call nvtxStartRange("PROJ-ALLREDUCE")
        call s_mpi_allreduce_vectors_sum(sl, sg, 2, 1)
        call nvtxEndRange
        rz = sg(1, 1); zq = sg(2, 1)

    end subroutine s_dot_rz_zq

    !> Global inner product over the interior
    impure function f_dot(a, b) result(res)

        real(wp), dimension(0:m,0:n,0:p), intent(in) :: a, b
        real(wp)                                     :: res, loc
        integer                                      :: j, k, l

        loc = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[loc]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    loc = loc + a(j, k, l)*b(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call nvtxStartRange("PROJ-ALLREDUCE")
        call s_mpi_allreduce_sum(loc, res)
        call nvtxEndRange

    end function f_dot

    !> Global inner product of the search direction with a vector
    impure function f_dot_dir(a) result(res)

        real(wp), dimension(0:m,0:n,0:p), intent(in) :: a
        real(wp)                                     :: res, loc
        integer                                      :: j, k, l, sh

        sh = mg_poff
        loc = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[loc]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    loc = loc + mg_e(f_mg_ix(1, j, k, l) + sh)*a(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call nvtxStartRange("PROJ-ALLREDUCE")
        call s_mpi_allreduce_sum(loc, res)
        call nvtxEndRange

    end function f_dot_dir

    !> Flat index of cell (i, j, k) on level lv, ghosts at -1 and the level size
    pure function f_mg_ix(lv, xi, yi, zi) result(idx)

        $:GPU_ROUTINE(function_name='f_mg_ix', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: lv, xi, yi, zi
        integer             :: idx

        idx = mg_o(lv) + zi*mg_sz(lv) + yi*mg_sy(lv) + xi

    end function f_mg_ix

    !> Row of level lv's operator at flat index idx: diagonal dg and conductance-weighted neighbor sum nb, the neighbor values read
    !! sh past idx (a vector with the level's layout stored elsewhere in mg_e). Each cell stores the conductance of its low faces;
    !! boundary faces carry zero unless they are rank or periodic seams, whose neighbor values sit in the ghost layer
    pure subroutine s_mg_row(lv, idx, sh, dg, nb)

        $:GPU_ROUTINE(function_name='s_mg_row', parallelism='[seq]', cray_inline=True)

        integer, intent(in)   :: lv, idx, sh
        real(wp), intent(out) :: dg, nb
        integer               :: sy, sz

        sy = mg_sy(lv); sz = mg_sz(lv)
        dg = mg_d(idx) + mg_kx(idx) + mg_kx(idx + 1)
        nb = mg_kx(idx)*mg_e(idx + sh - 1) + mg_kx(idx + 1)*mg_e(idx + sh + 1)
        if (num_dims > 1) then
            dg = dg + mg_ky(idx) + mg_ky(idx + sy)
            nb = nb + mg_ky(idx)*mg_e(idx + sh - sy) + mg_ky(idx + sy)*mg_e(idx + sh + sy)
        end if
        if (num_dims > 2) then
            dg = dg + mg_kz(idx) + mg_kz(idx + sz)
            nb = nb + mg_kz(idx)*mg_e(idx + sh - sz) + mg_kz(idx + sz)*mg_e(idx + sh + sz)
        end if

    end subroutine s_mg_row

    !> First (lo) and last (hi) fine child of coarse index ic along a direction with nc coarse and nf fine cells: an odd last child
    !! folds into the last coarse cell, and ic = nc (the ghost past the end) maps to the fine high face
    pure subroutine s_mg_children(ic, nc, nf, lo, hi)

        $:GPU_ROUTINE(function_name='s_mg_children', parallelism='[seq]', cray_inline=True)

        integer, intent(in)  :: ic, nc, nf
        integer, intent(out) :: lo, hi

        lo = 2*ic
        hi = merge(nf - 1, 2*ic + 1, ic == nc - 1)
        if (ic >= nc) hi = 2*ic

    end subroutine s_mg_children

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

    !> Set the operator A of the solve and build its Galerkin hierarchy. Level 1 takes the diagonal d and each cell's low-face
    !! conductances (cx(j, k, l) couples cells j - 1 and j; the high boundary face sits at m + 1, and n + 1, p + 1); a boundary face
    !! with conductance couples to the ghost past it, which an exchange fills across rank and periodic seams and which is zero
    !! elsewhere. A coarse cell sums its children's diagonals, and a coarse face the fine faces on it: exact for any coefficient
    !! jump, which rediscretizing averaged coefficients would not be. The bottom level, one cell per rank, is gathered and factored.
    impure subroutine s_mg_set_operator(d, cx, cy, cz)

        real(wp), dimension(0:m,0:n,0:p), intent(in) :: d
        real(wp), dimension(0:m + 1,0:n,0:p), intent(in) :: cx
        real(wp), dimension(0:m,0:merge(n + 1, 0, n > 0),0:p), intent(in) :: cy
        real(wp), dimension(0:m,0:n,0:merge(p + 1, 0, p > 0)), intent(in) :: cz
        integer :: lv, nx, ny, nz, cnx, cny, cnz, ii, jj, kk, a, b, c, idx, cidx
        integer :: gx, gy, gz, fi, fj, fk, a0, a1, b0, b1, c0, c1
        real(wp) :: sd, skx, sky, skz

        gx = mg_gx; gy = mg_gy; gz = mg_gz

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx]')
        do kk = -gz, p + gz
            do jj = -gy, n + gy
                do ii = -gx, m + gx
                    idx = f_mg_ix(1, ii, jj, kk)
                    mg_d(idx) = 0._wp; mg_kx(idx) = 0._wp; mg_ky(idx) = 0._wp; mg_kz(idx) = 0._wp
                    mg_e(idx) = 0._wp; mg_f(idx) = 0._wp; mg_r(idx) = 0._wp
                    if (ii >= 0 .and. jj >= 0 .and. kk >= 0) then
                        if (ii <= m .and. jj <= n .and. kk <= p) mg_d(idx) = d(ii, jj, kk)
                        if (jj <= n .and. kk <= p) mg_kx(idx) = cx(ii, jj, kk)
                        if (num_dims > 1 .and. ii <= m .and. kk <= p) mg_ky(idx) = cy(ii, jj, kk)
                        if (num_dims > 2 .and. ii <= m .and. jj <= n) mg_kz(idx) = cz(ii, jj, kk)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        do lv = 1, mg_nlev - 1
            nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); cnx = mg_nx(lv + 1); cny = mg_ny(lv + 1); cnz = mg_nz(lv + 1)
            ! Coarse cells and faces, ghosts included (a ghost holds only its low face, the level's high boundary face)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, a, b, c, idx, cidx, sd, skx, sky, skz, fi, fj, fk, a0, a1, b0, &
                                & b1, c0, c1]')
            do kk = -gz, cnz - 1 + gz
                do jj = -gy, cny - 1 + gy
                    do ii = -gx, cnx - 1 + gx
                        sd = 0._wp; skx = 0._wp; sky = 0._wp; skz = 0._wp
                        if (ii >= 0 .and. jj >= 0 .and. kk >= 0) then
                            fi = merge(nx, 2*ii, ii == cnx); fj = merge(ny, 2*jj, jj == cny); fk = merge(nz, 2*kk, kk == cnz)
                            call s_mg_children(ii, cnx, nx, a0, a1)
                            call s_mg_children(jj, cny, ny, b0, b1)
                            call s_mg_children(kk, cnz, nz, c0, c1)
                            $:GPU_LOOP(parallelism='[seq]')
                            do c = c0, c1
                                $:GPU_LOOP(parallelism='[seq]')
                                do b = b0, b1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do a = a0, a1
                                        if (ii < cnx .and. jj < cny .and. kk < cnz) sd = sd + mg_d(f_mg_ix(lv, a, b, c))
                                        if (a == a0 .and. jj < cny .and. kk < cnz) skx = skx + mg_kx(f_mg_ix(lv, fi, b, c))
                                        if (b == b0 .and. ii < cnx .and. kk < cnz) sky = sky + mg_ky(f_mg_ix(lv, a, fj, c))
                                        if (c == c0 .and. ii < cnx .and. jj < cny) skz = skz + mg_kz(f_mg_ix(lv, a, b, fk))
                                    end do
                                end do
                            end do
                        end if
                        cidx = f_mg_ix(lv + 1, ii, jj, kk)
                        mg_d(cidx) = sd; mg_kx(cidx) = skx; mg_ky(cidx) = sky; mg_kz(cidx) = skz
                        mg_e(cidx) = 0._wp; mg_f(cidx) = 0._wp
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        call s_mg_omega()
        if (mg_blev == mg_nlev) call s_mg_bottom_build()

    end subroutine s_mg_set_operator

    !> Coarse-correction scale of each level. Piecewise-constant aggregation under-corrects the Laplacian but represents the
    !! Helmholtz diagonal exactly, so the correction from level lv + 1 is scaled by 1 + (proj_mg_omega - 1)/(1 + r/mg_om_r0), r
    !! being that level's Helmholtz-to-Laplacian diagonal ratio: the full scale for Poisson-like (low Mach) levels, none once
    !! compressibility dominates. The form and mg_om_r0 are fitted to TGV iteration counts from Mach 0.1 to 0.001
    impure subroutine s_mg_omega()

        real(wp), dimension(2, mg_maxlev) :: sums, sums_glb
        real(wp)                          :: sd, sk, dg, nb, r1
        integer                           :: lv, ii, jj, kk, idx, nx, ny, nz

        do lv = 1, mg_nlev
            nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv)
            sd = 0._wp; sk = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb]', reduction='[[sd, sk]]', reductionOp='[+]')
            do kk = 0, nz - 1
                do jj = 0, ny - 1
                    do ii = 0, nx - 1
                        idx = f_mg_ix(lv, ii, jj, kk)
                        call s_mg_row(lv, idx, 0, dg, nb)
                        sd = sd + mg_d(idx)
                        sk = sk + dg - mg_d(idx)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            sums(:,lv) = [sd, sk]
        end do
        call s_mpi_allreduce_vectors_sum(sums(:,1:mg_nlev), sums_glb(:,1:mg_nlev), 2, mg_nlev)
        do lv = 2, mg_nlev
            mg_om(lv - 1) = 1._wp + (proj_mg_omega - 1._wp)/(1._wp + sums_glb(1, lv)/(mg_om_r0*max(sums_glb(2, lv), tiny(1._wp))))
        end do

        ! Stop at the first level whose spacing spans proj_mg_trunc screening lengths. The pressure equation screens over
        ! l = h/sqrt(6 r1), r1 the fine level's Helmholtz-to-Laplacian diagonal ratio (l is about c*dt, the acoustic CFL in cells);
        ! error longer than l has eigenvalues near the Helmholtz diagonal, which CG resolves cheaply, so the coarser, rank-spanning
        ! levels and the gathered bottom add little. The coarse ratios cannot show this: summing aggregates doubles a level's
        ! ratio, where the spacing squared quadruples
        mg_blev = mg_nlev
        if (proj_mg_trunc > 0._wp) then
            r1 = sums_glb(1, 1)/max(sums_glb(2, 1), tiny(1._wp))
            do lv = 2, mg_nlev
                if (6._wp*r1*4._wp**(lv - 1) >= proj_mg_trunc**2) then
                    mg_blev = lv
                    exit
                end if
            end do
        end if

        ! K-cycles pay off by keeping iterations flat as the hierarchy deepens with the grid; a truncated hierarchy stops at a
        ! fixed depth set by the screening length, so with proj_mg_kcycle = -1 it keeps V-cycles, cheaper and fewer iterations
        mg_kspace = mg_kfull
        if (proj_mg_kcycle < 0 .and. mg_blev < mg_nlev) mg_kspace = 0

    end subroutine s_mg_omega

    !> Every rank's bottom-level block size, for the global numbering of the bottom cells, and the bottom hierarchy's levels
    impure subroutine s_mg_bottom_layout()

        integer, dimension(3) :: pc
        integer               :: r, lv, i, ierr

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

    end subroutine s_mg_bottom_layout

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

    !> Gather every bottom cell's row (its diagonal sum and its face conductances with their global neighbors) into the global
    !! bottom stencil, coarsen it as the device levels are, and factor the last level; a face whose two sides fall in one cell (a
    !! periodic seam across a one-cell width, or a face inside an aggregated box) couples the cell to itself and cancels
    impure subroutine s_mg_bottom_build()

        real(wp), dimension(13, crs_cnt(proc_rank + 1)) :: row
        real(wp), dimension(13, crs_n)                  :: rows
        integer                                         :: i, j, k, c, q, g, jr, idx, off, ex, ey, gx, gy, gz, nt, ierr
        integer                                         :: lv, o, nn, ci, cj

        off = mg_off(mg_nlev); gx = mg_gx; gy = mg_gy; gz = mg_gz
        ex = mg_nx(mg_nlev) + 2*gx; ey = mg_ny(mg_nlev) + 2*gy; nt = ex*ey*(mg_nz(mg_nlev) + 2*gz)
        $:GPU_UPDATE(host='[mg_d(off + 1:off + nt), mg_kx(off + 1:off + nt), mg_ky(off + 1:off + nt), mg_kz(off + 1:off + nt)]')
        c = 0
        do k = 0, mg_nz(mg_nlev) - 1
            do j = 0, mg_ny(mg_nlev) - 1
                do i = 0, mg_nx(mg_nlev) - 1
                    c = c + 1; idx = f_mg_ix(mg_nlev, i, j, k)
                    row(:,c) = 0._wp
                    row(1, c) = mg_d(idx)
                    row(3, c) = mg_kx(idx); row(5, c) = mg_kx(idx + 1)
                    if (num_dims > 1) then
                        row(7, c) = mg_ky(idx); row(9, c) = mg_ky(idx + ex)
                    end if
                    if (num_dims > 2) then
                        row(11, c) = mg_kz(idx); row(13, c) = mg_kz(idx + ex*ey)
                    end if
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

    !> Bottom solve: gather the ranks' right-hand sides, run the bottom hierarchy's V-cycle on them, keep this rank's block
    impure subroutine s_mg_bottom_solve()

        real(wp), dimension(crs_cnt(proc_rank + 1)) :: v
        integer                                     :: i, j, k, c, idx, off, ex, ey, gx, gy, gz, nt, ierr

        off = mg_off(mg_nlev); gx = mg_gx; gy = mg_gy; gz = mg_gz
        ex = mg_nx(mg_nlev) + 2*gx; ey = mg_ny(mg_nlev) + 2*gy; nt = ex*ey*(mg_nz(mg_nlev) + 2*gz)
        $:GPU_UPDATE(host='[mg_f(off + 1:off + nt)]')
        c = 0
        do k = 0, mg_nz(mg_nlev) - 1
            do j = 0, mg_ny(mg_nlev) - 1
                do i = 0, mg_nx(mg_nlev) - 1
                    c = c + 1; v(c) = mg_f(f_mg_ix(mg_nlev, i, j, k))
                end do
            end do
        end do
        crs_b(1:c) = v
#ifdef MFC_MPI
        call MPI_ALLGATHERV(v, c, mpi_p, crs_b, crs_cnt, crs_disp, mpi_p, MPI_COMM_WORLD, ierr)
#endif
        if (crs_cg) then
            call s_crs_cg()
        else
            call s_crs_cycle(1)
        end if

        ! The ghosts stay zero: nothing reads the bottom level's
        mg_e(off + 1:off + nt) = 0._wp
        c = 0
        do k = 0, mg_nz(mg_nlev) - 1
            do j = 0, mg_ny(mg_nlev) - 1
                do i = 0, mg_nx(mg_nlev) - 1
                    c = c + 1; mg_e(f_mg_ix(mg_nlev, i, j, k)) = crs_x(crs_disp(proc_rank + 1) + c)
                end do
            end do
        end do
        $:GPU_UPDATE(device='[mg_e(off + 1:off + nt)]')

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
    !! crs_b into crs_x, from zero to crs_tol relative to crs_b. Every rank runs it on the same data, so all agree
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
            if (sqrt(sum(r*r)) <= crs_tol*bnorm) exit
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

    !> Fill the ghost layer of the block of mg_e at offset off, laid out as level lv, from the neighbors: one message per distinct
    !! neighbor rank holding all the sides it shares with this one, and periodic seams a rank owns alone copied in place
    impure subroutine s_mg_exchange(lv, off)

        integer, intent(in)    :: lv, off
        integer, dimension(12) :: req
        integer                :: i, nt, ierr

        call s_mg_sides(lv, off, 1)
        if (mg_nmsg == 0) return
#ifdef MFC_MPI
        if (rdma_mpi) then
            $:GPU_WAIT()
            #:call GPU_HOST_DATA(use_device_addr='[mg_sbuf, mg_rbuf]')
                do i = 1, mg_nmsg
                    call MPI_IRECV(mg_rbuf(mg_rbeg(i, lv) + 1), mg_rlen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, req(i), &
                                   & ierr)
                    call MPI_ISEND(mg_sbuf(mg_sbeg(i, lv) + 1), mg_slen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, &
                                   & req(mg_nmsg + i), ierr)
                end do
                call MPI_WAITALL(2*mg_nmsg, req, MPI_STATUSES_IGNORE, ierr)
            #:endcall GPU_HOST_DATA
        else
            nt = mg_sbeg(mg_nmsg, lv) + mg_slen(mg_nmsg, lv)
            $:GPU_UPDATE(host='[mg_sbuf(1:nt)]')
            do i = 1, mg_nmsg
                call MPI_IRECV(mg_rbuf(mg_rbeg(i, lv) + 1), mg_rlen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, req(i), ierr)
                call MPI_ISEND(mg_sbuf(mg_sbeg(i, lv) + 1), mg_slen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, &
                               & req(mg_nmsg + i), ierr)
            end do
            call MPI_WAITALL(2*mg_nmsg, req, MPI_STATUSES_IGNORE, ierr)
            nt = mg_rbeg(mg_nmsg, lv) + mg_rlen(mg_nmsg, lv)
            $:GPU_UPDATE(device='[mg_rbuf(1:nt)]')
        end if
#endif
        call s_mg_sides(lv, off, 2)

    end subroutine s_mg_exchange

    !> Message layout of every level, fixed by the neighbor ranks: this rank sends its sides in its own side order, and receives
    !! each neighbor's sides in that neighbor's side order, so both ends agree without further communication
    subroutine s_mg_exchange_layout()

        integer :: lv, q, qq, i, tot

        mg_nmsg = 0
        do q = 1, 6
            mg_smode(q) = 0
            if (mg_nbr(q) < 0) cycle
            mg_smode(q) = merge(3, 1, mg_nbr(q) == proc_rank)
            if (mg_smode(q) == 3 .or. any(mg_nlist(1:mg_nmsg) == mg_nbr(q))) cycle
            mg_nmsg = mg_nmsg + 1; mg_nlist(mg_nmsg) = mg_nbr(q)
        end do
        mg_boff = 0
        do lv = 1, mg_nlev
            tot = 0
            do i = 1, mg_nmsg
                mg_sbeg(i, lv) = tot
                do q = 1, 6
                    if (mg_nbr(q) /= mg_nlist(i)) cycle
                    mg_boff(q, lv, 1) = tot; tot = tot + f_side_size(lv, q)
                end do
                mg_slen(i, lv) = tot - mg_sbeg(i, lv)
            end do
            tot = 0
            do i = 1, mg_nmsg
                mg_rbeg(i, lv) = tot
                do qq = 1, 6
                    q = f_opposite(qq)
                    if (mg_nbr(q) /= mg_nlist(i)) cycle
                    mg_boff(q, lv, 2) = tot; tot = tot + f_side_size(lv, q)
                end do
                mg_rlen(i, lv) = tot - mg_rbeg(i, lv)
            end do
        end do
        $:GPU_UPDATE(device='[mg_smode, mg_boff]')

    end subroutine s_mg_exchange_layout

    !> Side q of the opposite end of the same direction
    pure integer function f_opposite(q)

        integer, intent(in) :: q

        f_opposite = q + merge(1, -1, mod(q, 2) == 1)

    end function f_opposite

    !> Cells in the layer of side q on level lv
    pure integer function f_side_size(lv, q)

        integer, intent(in) :: lv, q

        select case ((q + 1)/2)
        case (1); f_side_size = mg_ny(lv)*mg_nz(lv)
        case (2); f_side_size = mg_nx(lv)*mg_nz(lv)
        case default; f_side_size = mg_nx(lv)*mg_ny(lv)
        end select

    end function f_side_size

    !> All six sides of the block at offset off laid out as level lv, in one kernel. Phase 1 copies each periodic seam a rank owns
    !! alone into its ghost layer and packs each side facing another rank into mg_sbuf; phase 2 unpacks mg_rbuf into those ghosts
    impure subroutine s_mg_sides(lv, off, phase)

        integer, intent(in) :: lv, off, phase
        integer             :: nx, ny, nz, ex, ey, base, nmax, q, a, b, md, n1, n2, s1, s2, st, nn, lay, gho, idx, bidx, lvl, ph

        ! Dummies may alias host variables, which a device kernel must not reference

        lvl = lv; ph = phase
        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); ex = nx + 2*mg_gx; ey = ny + 2*mg_gy
        base = off + (mg_gz*ey + mg_gy)*ex + mg_gx + 1
        nmax = max(nx, ny, nz)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[q, a, b, md, n1, n2, s1, s2, st, nn, lay, gho, idx, bidx]')
        do q = 1, 6
            do b = 0, nmax - 1
                do a = 0, nmax - 1
                    ! Phase 1 acts on self seams (3) and packs (1), phase 2 unpacks (becomes 2)
                    md = mg_smode(q)
                    if (ph == 2) md = merge(2, 0, md == 1)
                    if ((q + 1)/2 == 1) then
                        nn = nx; st = 1; n1 = ny; s1 = ex; n2 = nz; s2 = ex*ey
                    else if ((q + 1)/2 == 2) then
                        nn = ny; st = ex; n1 = nx; s1 = 1; n2 = nz; s2 = ex*ey
                    else
                        nn = nz; st = ex*ey; n1 = nx; s1 = 1; n2 = ny; s2 = ex
                    end if
                    if (md > 0 .and. a < n1 .and. b < n2) then
                        lay = merge(0, nn - 1, mod(q, 2) == 1)
                        gho = merge(-1, nn, mod(q, 2) == 1)
                        idx = base + a*s1 + b*s2
                        bidx = mg_boff(q, lvl, merge(1, 2, md == 1)) + b*n1 + a + 1
                        if (md == 1) then
                            mg_sbuf(bidx) = mg_e(idx + lay*st)
                        else if (md == 2) then
                            mg_e(idx + gho*st) = mg_rbuf(bidx)
                        else
                            mg_e(idx + gho*st) = mg_e(idx + (nn - 1 - lay)*st)
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_sides

    !> Multigrid preconditioner on rs, returned in zs, coupled across ranks
    impure subroutine s_mg_precond()

        integer :: j, k, l, gx, gy, gz

        gx = mg_gx; gy = mg_gy; gz = mg_gz
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -gz, p + gz
            do k = -gy, n + gy
                do j = -gx, m + gx
                    mg_e(f_mg_ix(1, j, k, l)) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mg_f(f_mg_ix(1, j, k, l)) = rs(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_mg_cycle(1)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    zs(j, k, l) = mg_e(f_mg_ix(1, j, k, l))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_precond

    !> One cycle on level lv, from mg_f(lv) into mg_e(lv), which must start at zero, ghosts included. Each smoothing phase freezes
    !! the ghost layer for all its sweeps; pre-sweeps run red then black and post-sweeps black then red, so a V-cycle is a symmetric
    !! operator. The coarse problem gets one cycle, or every mg_kspace levels two flexible-CG steps (s_mg_kstep)
    recursive impure subroutine s_mg_cycle(lv)

        integer, intent(in) :: lv
        integer             :: i

        if (lv == mg_nlev) then
            call nvtxStartRange("PROJ-MG-BOTTOM")
            call s_mg_bottom_solve()
            call nvtxEndRange
            return
        end if
        ! Two exchanges per level: pre-smoothing starts from zero, so its frozen ghosts are already current; the residual and the
        ! post-smoothing each take one exchange
        do i = 1, proj_mg_sweeps
            call s_mg_smooth(lv, 0); call s_mg_smooth(lv, 1)
        end do
        call nvtxStartRange("PROJ-MG-EXCHANGE")
        call s_mg_exchange(lv, mg_off(lv))
        call nvtxEndRange
        ! At a truncation level (mg_blev < mg_nlev) the cycle is its two smoothing phases alone, still symmetric; CG takes the
        ! screened error longer than the level's spacing
        if (lv < mg_blev) then
            call s_mg_restrict(lv)
            if (mg_kspace > 0 .and. lv + 1 < mg_blev .and. mod(lv, max(mg_kspace, 1)) == 0) then
                call s_mg_kstep(lv + 1)
            else
                call s_mg_cycle(lv + 1)
            end if
            call s_mg_prolong(lv)
            call nvtxStartRange("PROJ-MG-EXCHANGE")
            call s_mg_exchange(lv, mg_off(lv))
            call nvtxEndRange
        end if
        do i = 1, proj_mg_sweeps
            call s_mg_smooth(lv, 1); call s_mg_smooth(lv, 0)
        end do

    end subroutine s_mg_cycle

    !> K-cycle on coarse level c (Notay & Vassilevski, NLAA 15, 2008): two flexible-CG steps on its system, mg_f(c) into mg_e(c),
    !! each preconditioned by one cycle; the second is skipped once the first has cut the residual by mg_kt. The dot products of a
    !! step are fused into one reduction
    recursive impure subroutine s_mg_kstep(c)

        integer, intent(in)       :: c
        real(wp), dimension(5, 1) :: sl, sg
        real(wp)                  :: s1, s2, s3, s4, s5, a1, rho2, cc, cd, dg, nb
        integer                   :: ii, jj, kk, idx, nx, ny, nz, off, gx, gy, gz

        nx = mg_nx(c); ny = mg_ny(c); nz = mg_nz(c); off = mg_off(c); gx = mg_gx; gy = mg_gy; gz = mg_gz

        ! Step 1: c1 = cycle(r), v = A c1
        call s_mg_cycle(c)
        call s_mg_exchange(c, off)
        s1 = 0._wp; s2 = 0._wp; s3 = 0._wp; s4 = 0._wp; s5 = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb]', reduction='[[s1, s2, s3, s4, s5]]', reductionOp='[+]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    idx = f_mg_ix(c, ii, jj, kk)
                    call s_mg_row(c, idx, 0, dg, nb)
                    mg_kr(idx) = mg_f(idx); mg_kc(idx) = mg_e(idx); mg_kv(idx) = dg*mg_e(idx) - nb
                    s1 = s1 + mg_e(idx)*mg_kv(idx); s2 = s2 + mg_e(idx)*mg_f(idx); s3 = s3 + mg_f(idx)*mg_f(idx)
                    s4 = s4 + mg_kv(idx)*mg_f(idx); s5 = s5 + mg_kv(idx)*mg_kv(idx)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        sl(:,1) = [s1, s2, s3, s4, s5]
        call s_mpi_allreduce_vectors_sum(sl, sg, 5, 1)
        s1 = sg(1, 1); s2 = sg(2, 1); s3 = sg(3, 1); s4 = sg(4, 1); s5 = sg(5, 1)
        a1 = s2/s1
        cc = a1; cd = 0._wp

        ! Step 2, unless |r - a1 v| <= mg_kt |r| already: d = cycle(r - a1 v), from zero
        if (s3 - 2._wp*a1*s4 + a1*a1*s5 > mg_kt*mg_kt*s3) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx]')
            do kk = -gz, nz - 1 + gz
                do jj = -gy, ny - 1 + gy
                    do ii = -gx, nx - 1 + gx
                        idx = f_mg_ix(c, ii, jj, kk)
                        mg_e(idx) = 0._wp
                        if (ii >= 0 .and. ii < nx .and. jj >= 0 .and. jj < ny .and. kk >= 0 .and. kk < nz) then
                            mg_f(idx) = mg_kr(idx) - a1*mg_kv(idx)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_mg_cycle(c)
            call s_mg_exchange(c, off)
            s2 = 0._wp; s4 = 0._wp; s5 = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb]', reduction='[[s2, s4, s5]]', reductionOp='[+]')
            do kk = 0, nz - 1
                do jj = 0, ny - 1
                    do ii = 0, nx - 1
                        idx = f_mg_ix(c, ii, jj, kk)
                        call s_mg_row(c, idx, 0, dg, nb)
                        s2 = s2 + mg_e(idx)*mg_kv(idx); s4 = s4 + mg_e(idx)*(dg*mg_e(idx) - nb); s5 = s5 + mg_e(idx)*mg_f(idx)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            sl(1:3,1) = [s2, s4, s5]
            call s_mpi_allreduce_vectors_sum(sl(1:3,:), sg(1:3,:), 3, 1)
            ! gamma = d.v, beta = d.A d, alpha2 = d.(r - a1 v); rho2 = beta - gamma^2/rho1
            rho2 = sg(2, 1) - sg(1, 1)**2/s1
            if (rho2 > 0._wp) then
                cd = sg(3, 1)/rho2; cc = a1 - sg(1, 1)*cd/s1
            end if
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    idx = f_mg_ix(c, ii, jj, kk)
                    mg_e(idx) = cc*mg_kc(idx) + cd*mg_e(idx)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_kstep

    !> One red-black Gauss-Seidel half-sweep of the given color on level lv, one thread per cell of that color
    impure subroutine s_mg_smooth(lv, color)

        integer, intent(in) :: lv, color
        integer             :: nx, ny, nz, ih, ii, jj, kk, idx, cl
        real(wp)            :: dg, nb

        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); cl = color
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ih, ii, jj, kk, idx, dg, nb]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ih = 0, (nx - 1)/2
                    ii = 2*ih + mod(jj + kk + cl, 2)
                    if (ii < nx) then
                        idx = f_mg_ix(lv, ii, jj, kk)
                        call s_mg_row(lv, idx, 0, dg, nb)
                        mg_e(idx) = (mg_f(idx) + nb)/dg
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_smooth

    !> Residual of level lv (ghosts current) summed onto the coarse right-hand side; restriction is the transpose of prolongation
    impure subroutine s_mg_restrict(lv)

        integer, intent(in) :: lv
        integer             :: nx, ny, nz, gx, gy, gz, cnx, cny, cnz, ii, jj, kk, idx, a, b, c
        integer             :: a0, a1, b0, b1, c0, c1
        real(wp)            :: dg, nb

        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); cnx = mg_nx(lv + 1); cny = mg_ny(lv + 1); cnz = mg_nz(lv + 1)
        gx = mg_gx; gy = mg_gy; gz = mg_gz

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    idx = f_mg_ix(lv, ii, jj, kk)
                    call s_mg_row(lv, idx, 0, dg, nb)
                    mg_r(idx) = mg_f(idx) - (dg*mg_e(idx) - nb)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! The coarse correction starts from zero, ghosts included, so its first sweep needs no exchange
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk]')
        do kk = -gz, cnz - 1 + gz
            do jj = -gy, cny - 1 + gy
                do ii = -gx, cnx - 1 + gx
                    mg_e(f_mg_ix(lv + 1, ii, jj, kk)) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Gather per coarse cell rather than scatter with atomics, so the sum order is fixed and results are reproducible
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, a, b, c, dg, a0, a1, b0, b1, c0, c1]')
        do kk = 0, cnz - 1
            do jj = 0, cny - 1
                do ii = 0, cnx - 1
                    dg = 0._wp
                    call s_mg_children(ii, cnx, nx, a0, a1)
                    call s_mg_children(jj, cny, ny, b0, b1)
                    call s_mg_children(kk, cnz, nz, c0, c1)
                    $:GPU_LOOP(parallelism='[seq]')
                    do c = c0, c1
                        $:GPU_LOOP(parallelism='[seq]')
                        do b = b0, b1
                            $:GPU_LOOP(parallelism='[seq]')
                            do a = a0, a1
                                dg = dg + mg_r(f_mg_ix(lv, a, b, c))
                            end do
                        end do
                    end do
                    mg_f(f_mg_ix(lv + 1, ii, jj, kk)) = dg
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_restrict

    !> Add each coarse correction, scaled by mg_om(lv) (see s_mg_omega), to its children; a scale in (0, 2) keeps the V-cycle
    !! symmetric positive definite
    impure subroutine s_mg_prolong(lv)

        integer, intent(in) :: lv
        integer             :: nx, ny, nz, cnx, cny, cnz, ii, jj, kk, idx
        real(wp)            :: om

        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); cnx = mg_nx(lv + 1); cny = mg_ny(lv + 1); cnz = mg_nz(lv + 1) &
                   & ; om = mg_om(lv)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    idx = f_mg_ix(lv, ii, jj, kk)
                    mg_e(idx) = mg_e(idx) + om*mg_e(f_mg_ix(lv + 1, min(ii/2, cnx - 1), min(jj/2, cny - 1), min(kk/2, cnz - 1)))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_prolong

    impure subroutine s_finalize_multigrid_module()

        @:DEALLOCATE(rs, zs, qs)
        @:DEALLOCATE(mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r, mg_kr, mg_kc, mg_kv)
        if (allocated(crs_l)) deallocate (crs_l)
        deallocate (crs_x, crs_b, crs_ad, crs_ak, crs_aj, crs_agg, crs_co, crs_cnt, crs_disp, crs_sz)
        @:UNPIN_HOST(mg_sbuf, mg_rbuf)
        @:DEALLOCATE(mg_sbuf, mg_rbuf)

    end subroutine s_finalize_multigrid_module

end module m_multigrid
