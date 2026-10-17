!>
!! @file
!! @brief Contains module m_multigrid

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Flexible PCG for the projection's pressure system A = D + face conductances (seven-point, SPD), preconditioned by one V-
!! or K-cycle of a Galerkin multigrid hierarchy that coarsens to one cell per rank, then solves the gathered bottom
!! (m_multigrid_bottom)
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

    !> PCG residual, preconditioned residual and operator image
    real(wp), allocatable, dimension(:,:,:) :: pcg_res, pcg_prec, pcg_Ad
    $:GPU_DECLARE(create='[pcg_res, pcg_prec, pcg_Ad]')
    !> PCG search direction, with level 1's ghost layer
    real(wp), allocatable, dimension(:,:,:) :: mg_dir
    $:GPU_DECLARE(create='[mg_dir]')
    integer :: mg_k_spacing_dflt  !< levels between K-cycle steps when a solve reaches the bottom
    !> levels between K-cycle steps in this solve: mg_k_spacing_dflt, or V-cycles once truncated (s_mg_omega)
    integer                        :: mg_k_spacing
    integer                        :: mg_bottom_lev  !< Level this solve stops at (s_mg_omega): mg_num_lev, or the truncation
    real(wp), dimension(mg_maxlev) :: mg_cc_scale    !< Coarse-correction scale into each level (s_mg_omega)

contains

    impure subroutine s_initialize_multigrid_module()

        ! Automatic choices by rank count; the thresholds are machine-dependent (docs: "Projection method iterative solve tuning")

        mg_k_spacing_dflt = proj_mg_kcycle
        if (mg_k_spacing_dflt < 0) mg_k_spacing_dflt = merge(2, 0, num_procs >= proj_mg_k_ranks)
        mg_k_spacing = mg_k_spacing_dflt

        @:ALLOCATE(pcg_res(0:m, 0:n, 0:p), pcg_prec(0:m, 0:n, 0:p), pcg_Ad(0:m, 0:n, 0:p))
        call s_initialize_multigrid_levels(mg_k_spacing_dflt > 0)
        @:ALLOCATE(mg_dir(-mg_ghost_x:m + mg_ghost_x, -mg_ghost_y:n + mg_ghost_y, -mg_ghost_z:p + mg_ghost_z))
        call s_initialize_multigrid_exchange()
        ! A K-cycle reaches the bottom several times per cycle, which multiplies the cost of an exact solve there
        call s_initialize_multigrid_bottom(proj_mg_bottom == 1 .or. (proj_mg_bottom == -1 .and. mg_k_spacing_dflt == 0 &
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
        integer                                         :: it, j, k, l, gx, gy, gz

        iters = 0
        res = 0._wp
        gx = mg_ghost_x
        gy = mg_ghost_y
        gz = mg_ghost_z

        ! The search direction's wall ghosts stay zero: those faces carry no conductance, and the ghosts must not hold garbage
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -gz, p + gz
            do k = -gy, n + gy
                do j = -gx, m + gx
                    mg_dir(j, k, l) = 0._wp
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
                    pcg_res(j, k, l) = b(j, k, l) - pcg_Ad(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Relative to the initial residual, i.e. to the pressure change being solved for; relative to b the tolerance would
        ! depend on the ambient pressure, which b carries in full
        bnorm = sqrt(f_dot(b, b))
        rnorm = sqrt(f_dot(pcg_res, pcg_res))
        rtol = max(proj_tol*rnorm, mg_res_floor*bnorm)
        r0 = rnorm

        if (rnorm > rtol) then
            call nvtxStartRange("PROJ-VCYCLE")
            call s_mg_precond()
            call nvtxEndRange
            call s_set_dir(pcg_prec)
            rz = f_dot(pcg_res, pcg_prec)

            do it = 1, proj_max_iters
                iters = iters + 1
                call s_apply_operator()
                dq = f_dot_dir(pcg_Ad)
                alpha = rz/dq
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            x(j, k, l) = x(j, k, l) + alpha*mg_dir(j, k, l)
                            pcg_res(j, k, l) = pcg_res(j, k, l) - alpha*pcg_Ad(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                rnorm = sqrt(f_dot(pcg_res, pcg_res))
                if (rnorm <= rtol) exit

                call nvtxStartRange("PROJ-VCYCLE")
                call s_mg_precond()
                call nvtxEndRange
                ! Flexible (Polak-Ribiere) beta, z.(r - r_old)/rz with r_old - r = alpha*q, as a K-cycle preconditioner varies
                call s_dot_rz_zq(rz_new, zq)
                beta = -alpha*zq/rz
                rz = rz_new
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            mg_dir(j, k, l) = pcg_prec(j, k, l) + beta*mg_dir(j, k, l)
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
        integer                                      :: j, k, l

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mg_dir(j, k, l) = v(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_set_dir

    !> pcg_Ad = A d for the search direction d, held in level 1's layout: A is level 1's operator, whose boundary faces carry no
    !! conductance unless they are rank or periodic seams, so one exchanged ghost layer is all it needs
    impure subroutine s_apply_operator()

        integer  :: j, k, l
        real(wp) :: diag, nb_sum

        call nvtxStartRange("PROJ-DIR-EXCHANGE")
        call s_mg_exchange(1, mg_dir)
        call nvtxEndRange
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, diag, nb_sum]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    call s_mg_row_dir(1, j, k, l, diag, nb_sum)
                    pcg_Ad(j, k, l) = diag*mg_dir(j, k, l) - nb_sum
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_operator

    !> Global r.z and q.z over the interior, in one reduction
    impure subroutine s_dot_rz_zq(rz, zq)

        real(wp), intent(out)     :: rz, zq
        real(wp), dimension(2, 1) :: sums_loc, sums_glb
        real(wp)                  :: a, b
        integer                   :: j, k, l

        a = 0._wp
        b = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[a, b]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    a = a + pcg_res(j, k, l)*pcg_prec(j, k, l)
                    b = b + pcg_Ad(j, k, l)*pcg_prec(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        sums_loc(:,1) = [a, b]
        call nvtxStartRange("PROJ-ALLREDUCE")
        call s_mpi_allreduce_vectors_sum(sums_loc, sums_glb, 2, 1)
        call nvtxEndRange
        rz = sums_glb(1, 1)
        zq = sums_glb(2, 1)

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
        integer                                      :: j, k, l

        loc = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[loc]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    loc = loc + mg_dir(j, k, l)*a(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call nvtxStartRange("PROJ-ALLREDUCE")
        call s_mpi_allreduce_sum(loc, res)
        call nvtxEndRange

    end function f_dot_dir

    #:for NAME, VEC, WHAT in [('s_mg_row', 'mg_lev(lv)%corr', 'the correction'), ('s_mg_row_dir', 'mg_dir', 'the PCG direction')]
        !> Row of level lv's operator at cell (j, k, l): diagonal diag and conductance-weighted neighbor sum nb_sum of ${WHAT}$.
        !! Each cell stores the conductance of its low faces; boundary faces carry zero unless they are rank or periodic seams,
        !! whose neighbor values sit in the ghost layer
        pure subroutine ${NAME}$(lv, j, k, l, diag, nb_sum)

            $:GPU_ROUTINE(function_name=NAME, parallelism='[seq]', cray_inline=True)

            integer, intent(in)   :: lv, j, k, l
            real(wp), intent(out) :: diag, nb_sum

            diag = mg_lev(lv)%diag(j, k, l)
            nb_sum = 0._wp
            #:for DIM, XYZ, LO, HI in [(1, 'x', 'j - 1, k, l', 'j + 1, k, l'), (2, 'y', 'j, k - 1, l', 'j, k + 1, l'), (3, 'z', &
                                        & 'j, k, l - 1', 'j, k, l + 1')]
                if (num_dims >= ${DIM}$) then
                    diag = diag + mg_lev(lv)%cond_${XYZ}$(j, k, l) + mg_lev(lv)%cond_${XYZ}$(${HI}$)
                    nb_sum = nb_sum + mg_lev(lv)%cond_${XYZ}$(j, k, l)*${VEC}$(${LO}$)
                    nb_sum = nb_sum + mg_lev(lv)%cond_${XYZ}$(${HI}$)*${VEC}$(${HI}$)
                end if
            #:endfor

        end subroutine ${NAME}$
    #:endfor

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

    !> Set A from the diagonal d and face conductances (cx(j, k, l) couples cells j and j + 1, from j = -1; a boundary face couples
    !! to its ghost) and build the Galerkin hierarchy: a coarse cell sums its children's diagonals, a coarse face its fine faces
    impure subroutine s_mg_set_operator(d, cx, cy, cz)

        real(wp), dimension(0:m,0:n,0:p), intent(in)  :: d
        real(wp), dimension(-1:m,0:n,0:p), intent(in) :: cx
        real(wp), dimension(0:m,-1:n,0:p), intent(in) :: cy
        real(wp), dimension(0:m,0:n,-1:p), intent(in) :: cz
        integer                                       :: lv, nx, ny, nz, nx_c, ny_c, nz_c, j, k, l, a, b, c
        integer                                       :: gx, gy, gz, j_f, k_f, l_f, j_lo, j_hi, k_lo, k_hi, l_lo, l_hi
        real(wp)                                      :: sum_diag, sum_cond_x, sum_cond_y, sum_cond_z

        gx = mg_ghost_x
        gy = mg_ghost_y
        gz = mg_ghost_z

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -gz, p + gz
            do k = -gy, n + gy
                do j = -gx, m + gx
                    #:for VAR in ['diag', 'cond_x', 'cond_y', 'cond_z', 'corr', 'rhs', 'res']
                        mg_lev(1)%${VAR}$(j, k, l) = 0._wp
                    #:endfor
                    if (j >= 0 .and. k >= 0 .and. l >= 0) then
                        if (j <= m .and. k <= n .and. l <= p) mg_lev(1)%diag(j, k, l) = d(j, k, l)
                        ! Each cell stores the conductance of its low faces
                        if (k <= n .and. l <= p) mg_lev(1)%cond_x(j, k, l) = cx(j - 1, k, l)
                        if (num_dims > 1 .and. j <= m .and. l <= p) mg_lev(1)%cond_y(j, k, l) = cy(j, k - 1, l)
                        if (num_dims > 2 .and. j <= m .and. k <= n) mg_lev(1)%cond_z(j, k, l) = cz(j, k, l - 1)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        do lv = 1, mg_num_lev - 1
            nx = mg_nx(lv)
            ny = mg_ny(lv)
            nz = mg_nz(lv)
            nx_c = mg_nx(lv + 1)
            ny_c = mg_ny(lv + 1)
            nz_c = mg_nz(lv + 1)
            ! Coarse cells and faces, ghosts included (a ghost holds only its low face, the level's high boundary face)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, a, b, c, sum_diag, sum_cond_x, sum_cond_y, sum_cond_z, j_f, k_f, &
                                & l_f, j_lo, j_hi, k_lo, k_hi, l_lo, l_hi]')
            do l = -gz, nz_c - 1 + gz
                do k = -gy, ny_c - 1 + gy
                    do j = -gx, nx_c - 1 + gx
                        sum_diag = 0._wp
                        sum_cond_x = 0._wp
                        sum_cond_y = 0._wp
                        sum_cond_z = 0._wp
                        if (j >= 0 .and. k >= 0 .and. l >= 0) then
                            j_f = merge(nx, 2*j, j == nx_c)
                            k_f = merge(ny, 2*k, k == ny_c)
                            l_f = merge(nz, 2*l, l == nz_c)
                            call s_mg_children(j, nx_c, nx, j_lo, j_hi)
                            call s_mg_children(k, ny_c, ny, k_lo, k_hi)
                            call s_mg_children(l, nz_c, nz, l_lo, l_hi)
                            $:GPU_LOOP(parallelism='[seq]')
                            do c = l_lo, l_hi
                                $:GPU_LOOP(parallelism='[seq]')
                                do b = k_lo, k_hi
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do a = j_lo, j_hi
                                        if (j < nx_c .and. k < ny_c .and. l < nz_c) sum_diag = sum_diag + mg_lev(lv)%diag(a, b, c)
                                        if (a == j_lo .and. k < ny_c .and. l < nz_c) sum_cond_x = sum_cond_x &
                                            & + mg_lev(lv)%cond_x(j_f, b, c)
                                        if (b == k_lo .and. j < nx_c .and. l < nz_c) sum_cond_y = sum_cond_y &
                                            & + mg_lev(lv)%cond_y(a, k_f, c)
                                        if (c == l_lo .and. j < nx_c .and. k < ny_c) sum_cond_z = sum_cond_z &
                                            & + mg_lev(lv)%cond_z(a, b, l_f)
                                    end do
                                end do
                            end do
                        end if
                        mg_lev(lv + 1)%diag(j, k, l) = sum_diag
                        mg_lev(lv + 1)%cond_x(j, k, l) = sum_cond_x
                        mg_lev(lv + 1)%cond_y(j, k, l) = sum_cond_y
                        mg_lev(lv + 1)%cond_z(j, k, l) = sum_cond_z
                        mg_lev(lv + 1)%corr(j, k, l) = 0._wp
                        mg_lev(lv + 1)%rhs(j, k, l) = 0._wp
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        call s_mg_omega()
        if (mg_bottom_lev == mg_num_lev) call s_mg_bottom_set()

    end subroutine s_mg_set_operator

    !> Hand the bottom level's operator to m_multigrid_bottom: each cell's diagonal and face conductances, x fastest
    impure subroutine s_mg_bottom_set()

        real(wp), dimension(7, mg_nx(mg_num_lev)*mg_ny(mg_num_lev)*mg_nz(mg_num_lev)) :: cell_rows
        integer                                                                       :: i, j, k, c, lv

        lv = mg_num_lev
        $:GPU_UPDATE(host='[mg_lev(lv)%diag, mg_lev(lv)%cond_x, mg_lev(lv)%cond_y, mg_lev(lv)%cond_z]')
        cell_rows = 0._wp
        c = 0
        do k = 0, mg_nz(lv) - 1
            do j = 0, mg_ny(lv) - 1
                do i = 0, mg_nx(lv) - 1
                    c = c + 1
                    cell_rows(1:3,c) = [mg_lev(lv)%diag(i, j, k), mg_lev(lv)%cond_x(i, j, k), mg_lev(lv)%cond_x(i + 1, j, k)]
                    if (num_dims > 1) cell_rows(4:5,c) = [mg_lev(lv)%cond_y(i, j, k), mg_lev(lv)%cond_y(i, j + 1, k)]
                    if (num_dims > 2) cell_rows(6:7,c) = [mg_lev(lv)%cond_z(i, j, k), mg_lev(lv)%cond_z(i, j, k + 1)]
                end do
            end do
        end do
        call s_mg_bottom_build(cell_rows)

    end subroutine s_mg_bottom_set

    !> Coarse-correction scale 1 + (proj_mg_omega - 1)/(1 + r/mg_om_r0) of each level, r its Helmholtz-to-Laplacian diagonal ratio
    !! (aggregation under-corrects only the Laplacian), then the truncation level and the K-cycle spacing
    impure subroutine s_mg_omega()

        real(wp), dimension(2, mg_maxlev) :: sums, sums_glb
        real(wp)                          :: sum_diag, sum_cond, diag, nb_sum, r1
        integer                           :: lv, j, k, l, nx, ny, nz

        do lv = 1, mg_num_lev
            nx = mg_nx(lv)
            ny = mg_ny(lv)
            nz = mg_nz(lv)
            sum_diag = 0._wp
            sum_cond = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, diag, nb_sum]', reduction='[[sum_diag, sum_cond]]', &
                                & reductionOp='[+]')
            do l = 0, nz - 1
                do k = 0, ny - 1
                    do j = 0, nx - 1
                        call s_mg_row(lv, j, k, l, diag, nb_sum)
                        sum_diag = sum_diag + mg_lev(lv)%diag(j, k, l)
                        sum_cond = sum_cond + diag - mg_lev(lv)%diag(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            sums(:,lv) = [sum_diag, sum_cond]
        end do
        call s_mpi_allreduce_vectors_sum(sums(:,1:mg_num_lev), sums_glb(:,1:mg_num_lev), 2, mg_num_lev)
        do lv = 2, mg_num_lev
            mg_cc_scale(lv - 1) = 1._wp + (proj_mg_omega - 1._wp)/(1._wp + sums_glb(1, lv)/(mg_om_r0*max(sums_glb(2, lv), &
                        & tiny(1._wp))))
        end do

        ! Stop at the first level whose spacing spans proj_mg_trunc screening lengths h/sqrt(6 r1) (about c*dt); CG resolves
        ! longer error cheaply. r1 is the fine level's ratio, as aggregation doubles a level's ratio while h^2 quadruples
        mg_bottom_lev = mg_num_lev
        if (proj_mg_trunc > 0._wp) then
            r1 = sums_glb(1, 1)/max(sums_glb(2, 1), tiny(1._wp))
            do lv = 2, mg_num_lev
                if (6._wp*r1*4._wp**(lv - 1) >= proj_mg_trunc**2) then
                    mg_bottom_lev = lv
                    exit
                end if
            end do
        end if

        ! A truncated hierarchy has a fixed depth, so proj_mg_kcycle = -1 keeps V-cycles there
        mg_k_spacing = mg_k_spacing_dflt
        if (proj_mg_kcycle < 0 .and. mg_bottom_lev < mg_num_lev) mg_k_spacing = 0

    end subroutine s_mg_omega

    !> Multigrid preconditioner on pcg_res, returned in pcg_prec, coupled across ranks
    impure subroutine s_mg_precond()

        integer :: j, k, l, gx, gy, gz

        gx = mg_ghost_x
        gy = mg_ghost_y
        gz = mg_ghost_z
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -gz, p + gz
            do k = -gy, n + gy
                do j = -gx, m + gx
                    mg_lev(1)%corr(j, k, l) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mg_lev(1)%rhs(j, k, l) = pcg_res(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_mg_cycle(1)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    pcg_prec(j, k, l) = mg_lev(1)%corr(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_precond

    !> Bottom level, rhs into corr, by the gathered solve of m_multigrid_bottom
    impure subroutine s_mg_bottom()

        real(wp), dimension(mg_nx(mg_num_lev)*mg_ny(mg_num_lev)*mg_nz(mg_num_lev)) :: v
        integer                                                                    :: i, j, k, c, lv

        lv = mg_num_lev
        $:GPU_UPDATE(host='[mg_lev(lv)%rhs]')
        c = 0
        do k = 0, mg_nz(lv) - 1
            do j = 0, mg_ny(lv) - 1
                do i = 0, mg_nx(lv) - 1
                    c = c + 1
                    v(c) = mg_lev(lv)%rhs(i, j, k)
                end do
            end do
        end do
        call s_mg_bottom_solve(v)
        ! The ghosts stay zero: nothing reads the bottom level's
        mg_lev(lv)%corr = 0._wp
        c = 0
        do k = 0, mg_nz(lv) - 1
            do j = 0, mg_ny(lv) - 1
                do i = 0, mg_nx(lv) - 1
                    c = c + 1
                    mg_lev(lv)%corr(i, j, k) = v(c)
                end do
            end do
        end do
        $:GPU_UPDATE(device='[mg_lev(lv)%corr]')

    end subroutine s_mg_bottom

    !> One cycle on level lv, from mg_lev(lv)%rhs into mg_lev(lv)%corr, which must start at zero, ghosts included. Each smoothing
    !! phase freezes the ghost layer for all its sweeps; pre-sweeps run red then black and post-sweeps black then red, so a V-cycle
    !! is a symmetric operator. The coarse problem gets one cycle, or every mg_k_spacing levels two flexible-CG steps (s_mg_kstep)
    recursive impure subroutine s_mg_cycle(lv)

        integer, intent(in) :: lv
        integer             :: i

        if (lv == mg_num_lev) then
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
        call s_mg_exchange(lv, mg_lev(lv)%corr)
        call nvtxEndRange
        ! At a truncation level (mg_bottom_lev < mg_num_lev) the cycle is its two smoothing phases alone, still symmetric; CG takes
        ! the screened error longer than the level's spacing
        if (lv < mg_bottom_lev) then
            call s_mg_restrict(lv)
            if (mg_k_spacing > 0 .and. lv + 1 < mg_bottom_lev .and. mod(lv, max(mg_k_spacing, 1)) == 0) then
                call s_mg_kstep(lv + 1)
            else
                call s_mg_cycle(lv + 1)
            end if
            call s_mg_prolong(lv)
            call nvtxStartRange("PROJ-MG-EXCHANGE")
            call s_mg_exchange(lv, mg_lev(lv)%corr)
            call nvtxEndRange
        end if
        do i = 1, proj_mg_sweeps
            call s_mg_smooth(lv, 1)
            call s_mg_smooth(lv, 0)
        end do

    end subroutine s_mg_cycle

    !> K-cycle on coarse level c, Notay and Vassilevski NLAA (2008): two flexible-CG steps on its system, mg_lev(c)%rhs into
    !! mg_lev(c)%corr, each preconditioned by one cycle; the second is skipped once the first has cut the residual by mg_kt. The dot
    !! products of a step are fused into one reduction
    recursive impure subroutine s_mg_kstep(c)

        integer, intent(in) :: c
        !> Global dot products. Step 1: c1.v, c1.r, r.r, v.r, v.v; step 2 reuses 2:4 for d.v, d.A d, d.(r - alpha1 v)
        real(wp), dimension(5, 1) :: dots
        real(wp)                  :: s1, s2, s3, s4, s5  !< kernel reduction accumulators
        real(wp)                  :: alpha1, rho2, coef_c1, coef_d, diag, nb_sum
        integer                   :: j, k, l, nx, ny, nz, gx, gy, gz

        nx = mg_nx(c)
        ny = mg_ny(c)
        nz = mg_nz(c)
        gx = mg_ghost_x
        gy = mg_ghost_y
        gz = mg_ghost_z

        ! Step 1: c1 = cycle(r), v = A c1
        call s_mg_cycle(c)
        call s_mg_exchange(c, mg_lev(c)%corr)
        s1 = 0._wp
        s2 = 0._wp
        s3 = 0._wp
        s4 = 0._wp
        s5 = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, diag, nb_sum]', reduction='[[s1, s2, s3, s4, s5]]', reductionOp='[+]')
        do l = 0, nz - 1
            do k = 0, ny - 1
                do j = 0, nx - 1
                    call s_mg_row(c, j, k, l, diag, nb_sum)
                    mg_lev(c)%k_rhs(j, k, l) = mg_lev(c)%rhs(j, k, l)
                    mg_lev(c)%k_dir(j, k, l) = mg_lev(c)%corr(j, k, l)
                    mg_lev(c)%k_Adir(j, k, l) = diag*mg_lev(c)%corr(j, k, l) - nb_sum
                    s1 = s1 + mg_lev(c)%corr(j, k, l)*mg_lev(c)%k_Adir(j, k, l)
                    s2 = s2 + mg_lev(c)%corr(j, k, l)*mg_lev(c)%rhs(j, k, l)
                    s3 = s3 + mg_lev(c)%rhs(j, k, l)*mg_lev(c)%rhs(j, k, l)
                    s4 = s4 + mg_lev(c)%k_Adir(j, k, l)*mg_lev(c)%rhs(j, k, l)
                    s5 = s5 + mg_lev(c)%k_Adir(j, k, l)*mg_lev(c)%k_Adir(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        dots(:,1) = [s1, s2, s3, s4, s5]
        call s_mpi_allreduce_vectors_sum(dots, dots, 5, 1)
        alpha1 = dots(2, 1)/dots(1, 1)
        coef_c1 = alpha1
        coef_d = 0._wp

        ! Step 2, unless |r - alpha1 v| <= mg_kt |r| already: d = cycle(r - alpha1 v), from zero
        if (dots(3, 1) - 2._wp*alpha1*dots(4, 1) + alpha1*alpha1*dots(5, 1) > mg_kt*mg_kt*dots(3, 1)) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = -gz, nz - 1 + gz
                do k = -gy, ny - 1 + gy
                    do j = -gx, nx - 1 + gx
                        mg_lev(c)%corr(j, k, l) = 0._wp
                        if (j >= 0 .and. j < nx .and. k >= 0 .and. k < ny .and. l >= 0 .and. l < nz) then
                            mg_lev(c)%rhs(j, k, l) = mg_lev(c)%k_rhs(j, k, l) - alpha1*mg_lev(c)%k_Adir(j, k, l)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_mg_cycle(c)
            call s_mg_exchange(c, mg_lev(c)%corr)
            s1 = 0._wp
            s2 = 0._wp
            s3 = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, diag, nb_sum]', reduction='[[s1, s2, s3]]', reductionOp='[+]')
            do l = 0, nz - 1
                do k = 0, ny - 1
                    do j = 0, nx - 1
                        call s_mg_row(c, j, k, l, diag, nb_sum)
                        s1 = s1 + mg_lev(c)%corr(j, k, l)*mg_lev(c)%k_Adir(j, k, l)
                        s2 = s2 + mg_lev(c)%corr(j, k, l)*(diag*mg_lev(c)%corr(j, k, l) - nb_sum)
                        s3 = s3 + mg_lev(c)%corr(j, k, l)*mg_lev(c)%rhs(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            dots(2:4,1) = [s1, s2, s3]
            call s_mpi_allreduce_vectors_sum(dots(2:4,:), dots(2:4,:), 3, 1)
            ! gamma = d.v, beta = d.A d, alpha2 = d.(r - alpha1 v); rho2 = beta - gamma^2/rho1
            rho2 = dots(3, 1) - dots(2, 1)**2/dots(1, 1)
            if (rho2 > 0._wp) then
                coef_d = dots(4, 1)/rho2
                coef_c1 = alpha1 - dots(2, 1)*coef_d/dots(1, 1)
            end if
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, nz - 1
            do k = 0, ny - 1
                do j = 0, nx - 1
                    mg_lev(c)%corr(j, k, l) = coef_c1*mg_lev(c)%k_dir(j, k, l) + coef_d*mg_lev(c)%corr(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_kstep

    !> One red-black Gauss-Seidel half-sweep of the given color on level lv, one thread per cell of that color
    impure subroutine s_mg_smooth(lv, color)

        integer, intent(in) :: lv, color
        integer             :: nx, ny, nz, ih, j, k, l
        real(wp)            :: diag, nb_sum

        nx = mg_nx(lv)
        ny = mg_ny(lv)
        nz = mg_nz(lv)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ih, j, k, l, diag, nb_sum]')
        do l = 0, nz - 1
            do k = 0, ny - 1
                do ih = 0, (nx - 1)/2
                    j = 2*ih + mod(k + l + color, 2)
                    if (j < nx) then
                        call s_mg_row(lv, j, k, l, diag, nb_sum)
                        mg_lev(lv)%corr(j, k, l) = (mg_lev(lv)%rhs(j, k, l) + nb_sum)/diag
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_smooth

    !> Residual of level lv (ghosts current) summed onto the coarse right-hand side; restriction is the transpose of prolongation
    impure subroutine s_mg_restrict(lv)

        integer, intent(in) :: lv
        integer             :: nx, ny, nz, gx, gy, gz, nx_c, ny_c, nz_c, j, k, l, a, b, c
        integer             :: j_lo, j_hi, k_lo, k_hi, l_lo, l_hi
        real(wp)            :: diag, nb_sum

        nx = mg_nx(lv)
        ny = mg_ny(lv)
        nz = mg_nz(lv)
        nx_c = mg_nx(lv + 1)
        ny_c = mg_ny(lv + 1)
        nz_c = mg_nz(lv + 1)
        gx = mg_ghost_x
        gy = mg_ghost_y
        gz = mg_ghost_z

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, diag, nb_sum]')
        do l = 0, nz - 1
            do k = 0, ny - 1
                do j = 0, nx - 1
                    call s_mg_row(lv, j, k, l, diag, nb_sum)
                    mg_lev(lv)%res(j, k, l) = mg_lev(lv)%rhs(j, k, l) - (diag*mg_lev(lv)%corr(j, k, l) - nb_sum)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! The coarse correction starts from zero, ghosts included, so its first sweep needs no exchange
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -gz, nz_c - 1 + gz
            do k = -gy, ny_c - 1 + gy
                do j = -gx, nx_c - 1 + gx
                    mg_lev(lv + 1)%corr(j, k, l) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Gather per coarse cell rather than scatter with atomics, so the sum order is fixed and results are reproducible
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, a, b, c, diag, j_lo, j_hi, k_lo, k_hi, l_lo, l_hi]')
        do l = 0, nz_c - 1
            do k = 0, ny_c - 1
                do j = 0, nx_c - 1
                    diag = 0._wp
                    call s_mg_children(j, nx_c, nx, j_lo, j_hi)
                    call s_mg_children(k, ny_c, ny, k_lo, k_hi)
                    call s_mg_children(l, nz_c, nz, l_lo, l_hi)
                    $:GPU_LOOP(parallelism='[seq]')
                    do c = l_lo, l_hi
                        $:GPU_LOOP(parallelism='[seq]')
                        do b = k_lo, k_hi
                            $:GPU_LOOP(parallelism='[seq]')
                            do a = j_lo, j_hi
                                diag = diag + mg_lev(lv)%res(a, b, c)
                            end do
                        end do
                    end do
                    mg_lev(lv + 1)%rhs(j, k, l) = diag
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_restrict

    !> Add each coarse correction, scaled by mg_cc_scale(lv) (see s_mg_omega), to its children; a scale in (0, 2) keeps the V-cycle
    !! symmetric positive definite
    impure subroutine s_mg_prolong(lv)

        integer, intent(in) :: lv
        integer             :: nx, ny, nz, nx_c, ny_c, nz_c, j, k, l
        real(wp)            :: cc_scale

        nx = mg_nx(lv)
        ny = mg_ny(lv)
        nz = mg_nz(lv)
        nx_c = mg_nx(lv + 1)
        ny_c = mg_ny(lv + 1)
        nz_c = mg_nz(lv + 1)
        cc_scale = mg_cc_scale(lv)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, nz - 1
            do k = 0, ny - 1
                do j = 0, nx - 1
                    mg_lev(lv)%corr(j, k, l) = mg_lev(lv)%corr(j, k, l) + cc_scale*mg_lev(lv + 1)%corr(min(j/2, nx_c - 1), &
                           & min(k/2, ny_c - 1), min(l/2, nz_c - 1))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_prolong

    impure subroutine s_finalize_multigrid_module()

        @:DEALLOCATE(pcg_res, pcg_prec, pcg_Ad, mg_dir)
        call s_finalize_multigrid_bottom()
        call s_finalize_multigrid_exchange()
        call s_finalize_multigrid_levels()

    end subroutine s_finalize_multigrid_module

end module m_multigrid
