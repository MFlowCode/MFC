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

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_nvtx
    use m_multigrid_levels
    use m_multigrid_exchange
    use m_multigrid_bottom

    implicit none

    private; public :: s_initialize_multigrid_module, s_mg_set_operator, s_mg_solve, s_finalize_multigrid_module

    real(wp), allocatable, dimension(:,:,:) :: rs, zs, qs  !< PCG residual, preconditioned residual and operator image
    $:GPU_DECLARE(create='[rs, zs, qs]')

    !> K-cycle vectors of the coarse levels (indexed as mg_e; allocated with mg_kfull > 0): the right-hand side (kr), the first
    !! step's preconditioned direction (kc) and its image under the level operator (kv)
    real(wp), allocatable, dimension(:) :: mg_kr, mg_kc, mg_kv
    $:GPU_DECLARE(create='[mg_kr, mg_kc, mg_kv]')
    integer :: mg_kfull  !< levels between K-cycle steps when a solve reaches the bottom
    integer :: mg_kspace  !< levels between K-cycle steps in this solve: mg_kfull, or V-cycles once truncated (s_mg_omega)
    integer :: mg_blev  !< Level this solve stops at (s_mg_omega): mg_nlev, or the truncation
    real(wp), dimension(mg_maxlev) :: mg_om  !< Coarse-correction scale into each level (s_mg_omega)

contains

    impure subroutine s_initialize_multigrid_module()

        integer :: kl, ku

        ! Automatic choices by rank count; the thresholds are machine-dependent (docs: "Projection method iterative solve tuning")

        mg_kfull = proj_mg_kcycle
        if (mg_kfull < 0) mg_kfull = merge(2, 0, num_procs >= proj_mg_k_ranks)
        mg_kspace = mg_kfull

        @:ALLOCATE(rs(0:m, 0:n, 0:p), zs(0:m, 0:n, 0:p), qs(0:m, 0:n, 0:p))
        call s_initialize_multigrid_levels()
        kl = mg_off(min(2, mg_nlev)) + 1; ku = merge(mg_poff, kl - 1, mg_kfull > 0)
        @:ALLOCATE(mg_kr(kl:ku), mg_kc(kl:ku), mg_kv(kl:ku))
        call s_initialize_multigrid_exchange()
        ! A K-cycle reaches the bottom several times per cycle, which multiplies the cost of an exact solve there
        call s_initialize_multigrid_bottom(proj_mg_bottom == 1 .or. (proj_mg_bottom == -1 .and. mg_kfull == 0 &
                                           & .and. num_procs <= proj_mg_cg_ranks))

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
        rtol = max(proj_tol*rnorm, mg_res_floor*bnorm)
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
        if (mg_blev == mg_nlev) call s_mg_bottom_set()

    end subroutine s_mg_set_operator

    !> Hand the bottom level's operator to m_multigrid_bottom: each cell's diagonal and face conductances, x fastest
    impure subroutine s_mg_bottom_set()

        real(wp), dimension(7, mg_nx(mg_nlev)*mg_ny(mg_nlev)*mg_nz(mg_nlev)) :: dk
        integer                                                              :: i, j, k, c, idx, lv, o, nt

        lv = mg_nlev; o = mg_off(lv); nt = mg_sz(lv)*(mg_nz(lv) + 2*mg_gz)
        $:GPU_UPDATE(host='[mg_d(o + 1:o + nt), mg_kx(o + 1:o + nt), mg_ky(o + 1:o + nt), mg_kz(o + 1:o + nt)]')
        dk = 0._wp; c = 0
        do k = 0, mg_nz(lv) - 1
            do j = 0, mg_ny(lv) - 1
                do i = 0, mg_nx(lv) - 1
                    c = c + 1; idx = f_mg_ix(lv, i, j, k)
                    dk(1:3,c) = [mg_d(idx), mg_kx(idx), mg_kx(idx + 1)]
                    if (num_dims > 1) dk(4:5,c) = [mg_ky(idx), mg_ky(idx + mg_sy(lv))]
                    if (num_dims > 2) dk(6:7,c) = [mg_kz(idx), mg_kz(idx + mg_sz(lv))]
                end do
            end do
        end do
        call s_mg_bottom_build(dk)

    end subroutine s_mg_bottom_set

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

    !> Bottom level, mg_f into mg_e, by the gathered solve of m_multigrid_bottom
    impure subroutine s_mg_bottom()

        real(wp), dimension(mg_nx(mg_nlev)*mg_ny(mg_nlev)*mg_nz(mg_nlev)) :: v
        integer                                                           :: i, j, k, c, lv, o, nt

        lv = mg_nlev; o = mg_off(lv); nt = mg_sz(lv)*(mg_nz(lv) + 2*mg_gz)
        $:GPU_UPDATE(host='[mg_f(o + 1:o + nt)]')
        c = 0
        do k = 0, mg_nz(lv) - 1
            do j = 0, mg_ny(lv) - 1
                do i = 0, mg_nx(lv) - 1
                    c = c + 1; v(c) = mg_f(f_mg_ix(lv, i, j, k))
                end do
            end do
        end do
        call s_mg_bottom_solve(v)
        ! The ghosts stay zero: nothing reads the bottom level's
        mg_e(o + 1:o + nt) = 0._wp
        c = 0
        do k = 0, mg_nz(lv) - 1
            do j = 0, mg_ny(lv) - 1
                do i = 0, mg_nx(lv) - 1
                    c = c + 1; mg_e(f_mg_ix(lv, i, j, k)) = v(c)
                end do
            end do
        end do
        $:GPU_UPDATE(device='[mg_e(o + 1:o + nt)]')

    end subroutine s_mg_bottom

    !> One cycle on level lv, from mg_f(lv) into mg_e(lv), which must start at zero, ghosts included. Each smoothing phase freezes
    !! the ghost layer for all its sweeps; pre-sweeps run red then black and post-sweeps black then red, so a V-cycle is a symmetric
    !! operator. The coarse problem gets one cycle, or every mg_kspace levels two flexible-CG steps (s_mg_kstep)
    recursive impure subroutine s_mg_cycle(lv)

        integer, intent(in) :: lv
        integer             :: i

        if (lv == mg_nlev) then
            call nvtxStartRange("PROJ-MG-BOTTOM")
            call s_mg_bottom()
            call nvtxEndRange
            return
        end if
        ! Two exchanges per level: pre-smoothing starts from zero, so its frozen ghosts are already current; the residual and the
        ! post-smoothing each take one exchange
        do i = 1, proj_mg_sweeps
            call s_mg_smooth(lv, 0)
            call s_mg_smooth(lv, 1)
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
            call s_mg_smooth(lv, 1)
            call s_mg_smooth(lv, 0)
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

        @:DEALLOCATE(rs, zs, qs, mg_kr, mg_kc, mg_kv)
        call s_finalize_multigrid_bottom()
        call s_finalize_multigrid_exchange()
        call s_finalize_multigrid_levels()

    end subroutine s_finalize_multigrid_module

end module m_multigrid
