!>
!! @file
!! @brief Contains module m_projection

#:include 'case.fpp'
#:include 'macros.fpp'

#! Accumulate one face of a multigrid row: conductance kf (read from `karr` at `kidx`)
#! couples to the neighbor at `nidx`; skipped when `cond` fails (level boundary)
#:def MG_FACE(cond, karr, kidx, nidx)
    if (${cond}$) then
        kf = ${karr}$ (${kidx}$)
        dg = dg + kf
        nb = nb + kf*mg_e(${nidx}$)
    end if
#:enddef

#! Multigrid row at flat index idx on a level of size nx, ny, nz: dg = diagonal,
#! nb = sum of conductance-weighted neighbor values
#:def MG_ROW()
    dg = mg_d(idx)
    nb = 0._wp
    @:MG_FACE(ii > 0, mg_kx, idx, idx - 1)
    @:MG_FACE(ii < nx - 1, mg_kx, idx + 1, idx + 1)
    @:MG_FACE(jj > 0, mg_ky, idx, idx - nx)
    @:MG_FACE(jj < ny - 1, mg_ky, idx + nx, idx + nx)
    @:MG_FACE(kk > 0, mg_kz, idx, idx - nx*ny)
    @:MG_FACE(kk < nz - 1, mg_kz, idx + nx*ny, idx + nx*ny)
#:enddef

!> All-Mach pressure projection (Fuster & Popinet, JCP 374, 2018). Advection uses a persistent face velocity; the pressure then
!! solves p - rho*c^2*tau^2*div(rho_f^-1 grad p) = p_adv - rho*c^2*tau*div(u*_f) with div, grad and the Laplacian all taken on
!! faces, so they compose exactly and the projected face velocity satisfies the discrete pressure equation. The solve is PCG
!! preconditioned by a rank-local geometric multigrid V-cycle: the preconditioner may ignore rank seams because CG restores the
!! global coupling, which keeps the MPI cost to one fine-level halo exchange per iteration.
!> @brief All-Mach pressure projection
module m_projection

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_boundary_common
    use m_eos
    use m_body_forces, only: s_compute_acceleration
    use m_riemann_state, only: Re_avg_rsx_vf, vel_src_rsx_vf, s_compute_interface_reynolds
    use m_surface_tension, only: c_divs

    implicit none

    private; public :: s_initialize_projection_module, s_projection_rhs, s_projection_face_props, s_projection_apply, &
        & s_finalize_projection_module

    integer, parameter :: mg_maxlev = 24
    integer, parameter :: mg_nu = 2                            !< symmetric smoothing sweeps per level
    real(wp), parameter :: res_floor = 1.e2_wp*epsilon(1._wp)  !< residual round-off floor, relative to the right-hand side
    real(wp), allocatable, dimension(:,:,:,:) :: uf            !< face velocity; index j is the face between cells j and j+1
    real(wp), allocatable, dimension(:,:,:) :: divu, rhs_p     !< div of uf, and the pressure transport rate
    real(wp), allocatable, dimension(:,:,:) :: pflx            !< upwind pressure flux on the faces of one direction
    real(wp), allocatable, dimension(:,:,:) :: p_stage, p_step0
    real(wp), allocatable, dimension(:,:,:) :: rhoc            !< star density with one ghost layer
    real(wp), allocatable, dimension(:,:,:) :: dcoef, bvec     !< SPD system: D_c and right-hand side
    real(wp), allocatable, dimension(:,:,:) :: xs, rs, zs, qs  !< PCG vectors
    real(stp), allocatable, dimension(:,:,:), target :: pk     !< search direction, then solution, with ghosts for the halo
    type(scalar_field), dimension(1) :: pk_sf
    $:GPU_DECLARE(create='[pk_sf]')
    $:GPU_DECLARE(create='[uf, divu, rhs_p, pflx, p_stage, p_step0, rhoc, dcoef, bvec, xs, rs, zs, qs, pk]')
    !> Well-balanced surface tension, with one ghost layer: curvature (1) and |grad c| (2), both zero outside the interface band
    real(wp), allocatable, dimension(:,:,:,:) :: kap
    $:GPU_DECLARE(create='[kap]')

    !> Multigrid hierarchy, flattened: level lv occupies mg_off(lv)+1 .. mg_off(lv)+nx*ny*nz, x fastest
    integer                             :: mg_nlev
    integer, dimension(mg_maxlev)       :: mg_nx, mg_ny, mg_nz, mg_off, mg_sx, mg_sy, mg_sz
    real(wp), allocatable, dimension(:) :: mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r
    $:GPU_DECLARE(create='[mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r]')

    logical, dimension(3) :: wall_lo, wall_hi    !< this rank owns a solid wall face on that side
    logical               :: faces_ready         !< uf has been seeded from the cell velocities
    logical               :: wb_st               !< well-balanced surface tension
    integer               :: gk0, gk1, gl0, gl1  !< y and z extents including one ghost layer where those directions exist

contains

    impure subroutine s_initialize_projection_module()

        integer :: lv, tot

#ifdef MFC_MIXED_PRECISION
        call s_mpi_abort('proj_method needs stp = wp; mixed precision is not supported')
#endif

        @:ALLOCATE(uf(-1:m + 1, -1:n + 1, -1:p + 1, 1:num_dims))
        @:ALLOCATE(divu(0:m, 0:n, 0:p), rhs_p(0:m, 0:n, 0:p), p_stage(0:m, 0:n, 0:p), p_step0(0:m, 0:n, 0:p))
        @:ALLOCATE(pflx(-1:m + 1, -1:n + 1, -1:p + 1), rhoc(-1:m + 1, -1:n + 1, -1:p + 1))
        @:ALLOCATE(dcoef(0:m, 0:n, 0:p), bvec(0:m, 0:n, 0:p))
        @:ALLOCATE(xs(0:m, 0:n, 0:p), rs(0:m, 0:n, 0:p), zs(0:m, 0:n, 0:p), qs(0:m, 0:n, 0:p))
        @:ALLOCATE(pk(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))

        pk_sf(1)%sf => pk
        $:GPU_ENTER_DATA(copyin='[pk_sf(1)%sf]')
        $:GPU_ENTER_DATA(attach='[pk_sf(1)%sf]')

        ! Rank-local hierarchy: a direction halves while its local size is even
        mg_nx(1) = m + 1; mg_ny(1) = n + 1; mg_nz(1) = p + 1
        mg_nlev = 1
        do while (mg_nlev < mg_maxlev)
            lv = mg_nlev
            mg_sx(lv) = merge(2, 1, mod(mg_nx(lv), 2) == 0)
            mg_sy(lv) = merge(2, 1, mod(mg_ny(lv), 2) == 0)
            mg_sz(lv) = merge(2, 1, mod(mg_nz(lv), 2) == 0)
            if (mg_sx(lv)*mg_sy(lv)*mg_sz(lv) == 1) exit
            mg_nx(lv + 1) = mg_nx(lv)/mg_sx(lv); mg_ny(lv + 1) = mg_ny(lv)/mg_sy(lv); mg_nz(lv + 1) = mg_nz(lv)/mg_sz(lv)
            mg_nlev = lv + 1
        end do
        tot = 0
        do lv = 1, mg_nlev
            mg_off(lv) = tot
            tot = tot + mg_nx(lv)*mg_ny(lv)*mg_nz(lv)
        end do
        @:ALLOCATE(mg_d(tot), mg_kx(tot), mg_ky(tot), mg_kz(tot), mg_e(tot), mg_f(tot), mg_r(tot))

        faces_ready = .false.
        gk0 = merge(-1, 0, n > 0); gk1 = merge(n + 1, n, n > 0)
        gl0 = merge(-1, 0, p > 0); gl1 = merge(p + 1, p, p > 0)

        wb_st = surface_tension .and. surface_tension_model == surface_tension_model_well_balanced
        if (wb_st) then
            @:ALLOCATE(kap(-1:m + 1, gk0:gk1, gl0:gl1, 1:2))
        else
            @:ALLOCATE(kap(0:0, 0:0, 0:0, 1:2))
        end if

        #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            wall_lo(${D}$) = any(bc_${XYZ}$%beg == [BC_REFLECTIVE, BC_SLIP_WALL, BC_NO_SLIP_WALL])
            wall_hi(${D}$) = any(bc_${XYZ}$%end == [BC_REFLECTIVE, BC_SLIP_WALL, BC_NO_SLIP_WALL])
        #:endfor

    end subroutine s_initialize_projection_module

    !> Face velocity from the average of the cell velocities either side (primitive velocities, with ghosts)
    subroutine s_projection_init_faces(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: j, k, l

        #:for D, IP1, LB, KB, JB in [(1, 'j + 1, k, l', 0, 0, -1), (2, 'j, k + 1, l', 0, -1, 0), (3, 'j, k, l + 1', -1, 0, 0)]
            if (num_dims >= ${D}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            uf(j, k, l, ${D}$) = 0.5_wp*(real(q_prim_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(j, k, l), &
                               & wp) + real(q_prim_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(${IP1}$), wp))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor
        call s_zero_wall_faces()

    end subroutine s_projection_init_faces

    !> Normal velocity vanishes on solid walls
    subroutine s_zero_wall_faces()

        integer :: j, k, l

        #:set UB = {'j': 'm', 'k': 'n', 'l': 'p'}
        #:for D, NV, LO, HI in [(1, 'j', '-1, k, l', 'm, k, l'), (2, 'k', 'j, -1, l', 'j, n, l'), (3, 'l', 'j, k, -1', 'j, k, p')]
            #:set TV = [v for v in ['l', 'k', 'j'] if v != NV]
            #:for SIDE, FACE in [('lo', LO), ('hi', HI)]
                if (num_dims >= ${D}$) then
                    if (wall_${SIDE}$(${D}$)) then
                        $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                        do ${TV[0]}$ = 0, ${UB[TV[0]]}$
                            do ${TV[1]}$ = 0, ${UB[TV[1]]}$
                                uf(${FACE}$, ${D}$) = 0._wp
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                end if
            #:endfor
        #:endfor

    end subroutine s_zero_wall_faces

    !> Advective right-hand side of one direction sweep. Every quantity is carried by the projected face velocity and upwinded on
    !! its sign; the momentum flux is the summed partial-density flux times the upwind velocity, so mass and momentum move with one
    !! operator. Energy is left at zero here and rebuilt from the equation of state after the pressure solve.
    subroutine s_projection_rhs(id, qfl_rs, qfr_rs, q_prim_vf, flux_vf, rhs_vf)

        integer, intent(in)                                                                 :: id
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: qfl_rs, qfr_rs
        type(scalar_field), dimension(sys_size), intent(in)                                 :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout)                              :: flux_vf, rhs_vf
        real(wp)                                                                            :: vf, a_up, ar_up, fm
        logical                                                                             :: up_l
        integer                                                                             :: i, j, k, l

        if (id == 1) then
            if (.not. faces_ready) then
                call s_projection_init_faces(q_prim_vf)
                faces_ready = .true.
            end if
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu(j, k, l) = f_div_uf(j, k, l)
                        p_stage(j, k, l) = real(q_prim_vf(eqn_idx%E)%sf(j, k, l), wp)
                        ! Transport sources q*div(u), which make alpha and p advect rather than compress
                        rhs_p(j, k, l) = p_stage(j, k, l)*divu(j, k, l)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, sys_size
                            rhs_vf(i)%sf(j, k, l) = 0._stp
                        end do
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            rhs_vf(i)%sf(j, k, l) = real(real(q_prim_vf(i)%sf(j, k, l), wp)*divu(j, k, l), stp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        #:for D, SV, COORDS, JB, KB, LB, DXV in [(1, 'j', '{SI}, k, l', -1, 0, 0, 'dx'), &
            (2, 'k', 'j, {SI}, l', 0, -1, 0, 'dy'), (3, 'l', 'j, k, {SI}', 0, 0, -1, 'dz')]
            #:set SF = lambda offs: COORDS.format(SI=SV + offs)
            if (id == ${D}$) then
                ! Face fluxes. Left state of face j is the right edge of cell j, right state the left edge of cell j+1
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, vf, a_up, ar_up, fm, up_l]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            vf = uf(j, k, l, ${D}$)
                            up_l = vf >= 0._wp
                            fm = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                ! Reconstructed partial densities upwinded directly. Taking the phase density alpha_rho/alpha from
                                ! the upwind cell instead diverges where a phase is vanishing: both are round-off there, and their
                                ! ratio (seen at 6e10) times the face alpha flux injects mass
                                if (up_l) then
                                    a_up = qfl_rs(${SF('')}$, eqn_idx%adv%beg + i - 1)
                                    ar_up = qfl_rs(${SF('')}$, i)
                                else
                                    a_up = qfr_rs(${SF(' + 1')}$, eqn_idx%adv%beg + i - 1)
                                    ar_up = qfr_rs(${SF(' + 1')}$, i)
                                end if
                                flux_vf(eqn_idx%adv%beg + i - 1)%sf(${SF('')}$) = real(a_up*vf, stp)
                                flux_vf(i)%sf(${SF('')}$) = real(ar_up*vf, stp)
                                fm = fm + ar_up*vf
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                if (up_l) then
                                    flux_vf(eqn_idx%mom%beg + i - 1)%sf(${SF('')}$) = real(fm*qfl_rs(${SF('')}$, &
                                            & eqn_idx%mom%beg + i - 1), stp)
                                else
                                    flux_vf(eqn_idx%mom%beg + i - 1)%sf(${SF('')}$) = real(fm*qfr_rs(${SF(' + 1')}$, &
                                            & eqn_idx%mom%beg + i - 1), stp)
                                end if
                            end do
                            if (up_l) then
                                pflx(j, k, l) = vf*qfl_rs(${SF('')}$, eqn_idx%E)
                            else
                                pflx(j, k, l) = vf*qfr_rs(${SF(' + 1')}$, eqn_idx%E)
                            end if
                            ! The color function is not reconstructed; upwind its cell value
                            if (surface_tension) then
                                if (up_l) then
                                    flux_vf(eqn_idx%c)%sf(${SF('')}$) = real(vf*real(q_prim_vf(eqn_idx%c)%sf(${SF('')}$), wp), stp)
                                else
                                    flux_vf(eqn_idx%c)%sf(${SF('')}$) = real(vf*real(q_prim_vf(eqn_idx%c)%sf(${SF(' + 1')}$), &
                                            & wp), stp)
                                end if
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, sys_size
                                if (i /= eqn_idx%E) then
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + real((real(flux_vf(i)%sf(${SF(' - 1')}$), &
                                           & wp) - real(flux_vf(i)%sf(j, k, l), wp))/${DXV}$(${SV}$), stp)
                                end if
                            end do
                            rhs_p(j, k, l) = rhs_p(j, k, l) + (pflx(${SF(' - 1')}$) - pflx(j, k, l))/${DXV}$(${SV}$)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_rhs

    !> Face data that the viscous and capillary source fluxes otherwise take from a Riemann solve: interface Reynolds numbers, the
    !! mean face velocity (read only by their energy terms, which the projection rebuilds from the EOS), and for surface tension the
    !! face velocity whose divergence makes the color function advect
    subroutine s_projection_face_props(id, qfl_rs, qfr_rs, flux_src_vf)

        integer, intent(in)                                                                 :: id
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: qfl_rs, qfr_rs
        type(scalar_field), dimension(sys_size), intent(inout)                              :: flux_src_vf

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: al, ar
        #:else
            real(wp), dimension(num_fluids) :: al, ar
        #:endif
        real(wp), dimension(2) :: re_l, re_r
        integer                :: i, j, k, l, rs1, rs2

        rs1 = Re_size(1); rs2 = Re_size(2)

        #:for D, SV, COORDS, JB, KB, LB in [(1, 'j', '{SI}, k, l', -1, 0, 0), (2, 'k', 'j, {SI}, l', 0, -1, 0), &
            (3, 'l', 'j, k, {SI}', 0, 0, -1)]
            #:set SF = lambda offs: COORDS.format(SI=SV + offs)
            if (id == ${D}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, al, ar, re_l, re_r]', firstprivate='[rs1, rs2]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    al(i) = qfl_rs(${SF('')}$, eqn_idx%adv%beg + i - 1)
                                    ar(i) = qfr_rs(${SF(' + 1')}$, eqn_idx%adv%beg + i - 1)
                                end do
                                call s_compute_interface_reynolds(al, re_l, rs1, rs2)
                                call s_compute_interface_reynolds(ar, re_r, rs1, rs2)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_avg_rsx_vf(j, k, l, i) = 2._wp/(1._wp/re_l(i) + 1._wp/re_r(i))
                                end do
                            end if
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
                                vel_src_rsx_vf(j, k, l, i) = 0.5_wp*(qfl_rs(${SF('')}$, &
                                               & eqn_idx%mom%beg + i - 1) + qfr_rs(${SF(' + 1')}$, eqn_idx%mom%beg + i - 1))
                            end do
                            if (surface_tension) flux_src_vf(eqn_idx%adv%beg)%sf(j, k, l) = real(uf(j, k, l, ${D}$), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_face_props

    !> Divergence of the face velocity in cell (j, k, l): the one operator the transport sources and the pressure equation share
    function f_div_uf(j, k, l) result(dv)

        $:GPU_ROUTINE(function_name='f_div_uf', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: j, k, l
        real(wp)            :: dv

        dv = (uf(j, k, l, 1) - uf(j - 1, k, l, 1))/dx(j)
        if (num_dims > 1) dv = dv + (uf(j, k, l, 2) - uf(j, k - 1, l, 2))/dy(k)
        if (num_dims > 2) dv = dv + (uf(j, k, l, 3) - uf(j, k, l - 1, 3))/dz(l)

    end function f_div_uf

    !> Face conductance A_f/(rho_f*d_f) with the arithmetic face density: the inertia of a face volume straddling an interface,
    !! which a heavy phase resting on a light one needs to stay at rest
    pure function f_cond(ra, rb, area, dist) result(kf)

        $:GPU_ROUTINE(function_name='f_cond', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: ra, rb, area, dist
        real(wp)             :: kf

        kf = 2._wp*area/(max(ra + rb, sgm_eps)*dist)

    end function f_cond

    !> Well-balanced (Brackbill CSF) capillary acceleration of a face: sigma*kappa_f*(c_b - c_a)/(d_f*rho_f), with rho_f the same
    !! arithmetic face density as the pressure operator, so a constant curvature is balanced exactly by a pressure jump. kappa_f is
    !! the |grad c|-weighted mean of the adjacent cells' curvature (weights wa, wb)
    pure function f_capillary_accel(ka, kb, wa, wb, ca, cb, ra, rb, dist) result(acc)

        $:GPU_ROUTINE(function_name='f_capillary_accel', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: ka, kb, wa, wb, ca, cb, ra, rb, dist
        real(wp)             :: acc

        acc = 0._wp
        if (wa + wb > 0._wp) acc = 2._wp*sigma*(wa*ka + wb*kb)/(wa + wb)*(cb - ca)/(dist*max(ra + rb, sgm_eps))

    end function f_capillary_accel

    !> Interface normal component, grad_d(c)/|grad(c)|; zero outside the interface band
    pure function f_normal(gd, g) result(nd)

        $:GPU_ROUTINE(function_name='f_normal', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: gd, g
        real(wp)             :: nd

        nd = 0._wp
        if (g > capillary_cutoff) nd = gd/g

    end function f_normal

    !> Curvature kappa = -div(n) in the interface band, from the color-function gradient s_get_capillary left with filled ghosts
    subroutine s_compute_curvature()

        real(wp) :: kv
        integer  :: j, k, l, k0, k1, l0, l1

        k0 = gk0; k1 = gk1; l0 = gl0; l1 = gl1
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, kv]')
        do l = l0, l1
            do k = k0, k1
                do j = -1, m + 1
                    kap(j, k, l, 1) = 0._wp
                    kap(j, k, l, 2) = 0._wp
                    if (real(c_divs(num_dims + 1)%sf(j, k, l), wp) > capillary_cutoff) then
                        #:set NRM = lambda d, &
                            & idx: f"f_normal(real(c_divs({d})%sf({idx}), wp), real(c_divs(num_dims + 1)%sf({idx}), wp))"
                        kv = -(${NRM(1, 'j + 1, k, l')}$ - ${NRM(1, 'j - 1, k, l')}$)/(x_cc(j + 1) - x_cc(j - 1))
                        if (num_dims > 1) kv = kv - (${NRM(2, 'j, k + 1, l')}$ - ${NRM(2, 'j, k - 1, l')}$)/(y_cc(k + 1) - y_cc(k &
                            & - 1))
                        if (num_dims > 2) kv = kv - (${NRM(3, 'j, k, l + 1')}$ - ${NRM(3, 'j, k, l - 1')}$)/(z_cc(l + 1) - z_cc(l &
                            & - 1))
                        kap(j, k, l, 1) = kv
                        kap(j, k, l, 2) = real(c_divs(num_dims + 1)%sf(j, k, l), wp)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_curvature

    !> Pressure solve and correction on the blended (star) state of one RK stage
    impure subroutine s_projection_apply(q_cons_vf, bc_type, pb_in, mv_in, q_T_sf, rkc1, rkc2, rkc3, rkc4, stage)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(scalar_field), intent(inout) :: q_T_sf
        real(wp), intent(in) :: rkc1, rkc2, rkc3, rkc4
        integer, intent(in) :: stage
        real(wp) :: tau, rho, gam, pinf, qv, rc2, dv
        real(wp) :: vol, ke, ga, gf
        real(wp), dimension(3) :: acc
        logical :: wlo, whi, wbl

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: ar, al
        #:else
            real(wp), dimension(num_fluids) :: ar, al
        #:endif
        integer :: i, j, k, l, k0, k1, l0, l1

        tau = rkc3*dt/rkc4
        k0 = gk0; k1 = gk1; l0 = gl0; l1 = gl1

        call s_populate_variables_buffers(bc_type, q_cons_vf, pb_in, mv_in, q_T_sf)

        acc = 0._wp
        if (bodyForces) then
            call s_compute_acceleration(mytime)
            #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (bf_${XYZ}$) acc(${D}$) = accel_bf(${D}$)
            #:endfor
        end if

        ! Star density with its ghosts, and the face predictor from the star cell velocities
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
        do l = l0, l1
            do k = k0, k1
                do j = -1, m + 1
                    rhoc(j, k, l) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        rhoc(j, k, l) = rhoc(j, k, l) + real(q_cons_vf(i)%sf(j, k, l), wp)
                    end do
                    rhoc(j, k, l) = max(rhoc(j, k, l), sgm_eps)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (wb_st) call s_compute_curvature()
        wbl = wb_st

        ! Body forces and well-balanced surface tension enter on faces, where they meet the pressure gradient that balances them
        #:for D, DXV, SV, IP1, LB, KB, JB in [(1, 'dx', 'j', 'j + 1, k, l', 0, 0, -1), (2, 'dy', 'k', 'j, k + 1, l', 0, -1, 0), &
            (3, 'dz', 'l', 'j, k, l + 1', -1, 0, 0)]
            if (num_dims >= ${D}$) then
                ga = tau*acc(${D}$)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            uf(j, k, l, ${D}$) = 0.5_wp*(real(q_cons_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(j, k, l), wp)/rhoc(j, k, &
                               & l) + real(q_cons_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(${IP1}$), wp)/rhoc(${IP1}$)) + ga
                            if (wbl) uf(j, k, l, ${D}$) = uf(j, k, l, ${D}$) + tau*f_capillary_accel(kap(j, k, l, 1), &
                                & kap(${IP1}$, 1), kap(j, k, l, 2), kap(${IP1}$, 2), real(q_cons_vf(eqn_idx%c)%sf(j, k, l), wp), &
                                & real(q_cons_vf(eqn_idx%c)%sf(${IP1}$), wp), rhoc(j, k, l), rhoc(${IP1}$), &
                                & 0.5_wp*(${DXV}$(${SV}$) + ${DXV}$(${SV}$ + 1)))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor
        call s_zero_wall_faces()

        ! SPD system D_c p + sum_f K_f (p - p_nb) = b, the Helmholtz row scaled by V_c/(rho c^2 tau^2)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho, gam, pinf, qv, rc2, dv, vol, ar, al]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (stage == 1) p_step0(j, k, l) = p_stage(j, k, l)
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ar(i) = real(q_cons_vf(i)%sf(j, k, l), wp)
                        al(i) = real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                    end do
                    call s_compute_mixture_coefficients(ar, al, rho, gam, pinf, qv)
                    ! Allaire's model advects alpha, so gamma_mix and pi_inf_mix are advected and Dp/Dt = -K div(u) with K the
                    ! mixture bulk modulus (not Wood's, which belongs to the Kapila model)
                    rc2 = max(f_bulk_modulus(p_stage(j, k, l), gam, pinf), sgm_eps)
                    dv = f_div_uf(j, k, l)
                    vol = dx(j)
                    if (num_dims > 1) vol = vol*dy(k)
                    if (num_dims > 2) vol = vol*dz(l)
                    dcoef(j, k, l) = vol/(rc2*tau*tau)
                    bvec(j, k, l) = dcoef(j, k, l)*((rkc1*p_stage(j, k, l) + rkc2*p_step0(j, k, l) + rkc3*dt*rhs_p(j, k, &
                         & l))/rkc4 - rc2*tau*dv)
                    xs(j, k, l) = p_stage(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_pcg_solve(bc_type)

        ! Face correction with the operator's own conductance (area 1), so div(uf) matches the solved pressure exactly. Cells take
        ! the mean of their faces' net acceleration (body force less pressure gradient, zero on walls): a hydrostatic balance on
        ! the faces then leaves the cells at rest too, and for uniform density this is the centered pressure gradient
        #:for D, DXV, SV, UB, IP1, IM1, LB, KB, JB in [(1, 'dx', 'j', 'm', 'j + 1, k, l', 'j - 1, k, l', 0, 0, -1), &
            (2, 'dy', 'k', 'n', 'j, k + 1, l', 'j, k - 1, l', 0, -1, 0), (3, 'dz', 'l', 'p', 'j, k, l + 1', 'j, k, l - 1', -1, 0, &
             & 0)]
            if (num_dims >= ${D}$) then
                ga = tau*acc(${D}$)
                wlo = wall_lo(${D}$); whi = wall_hi(${D}$)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, gf]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            gf = tau*f_cond(rhoc(j, k, l), rhoc(${IP1}$), 1._wp, &
                                            & 0.5_wp*(${DXV}$(${SV}$) + ${DXV}$(${SV}$ + 1)))*(real(pk(${IP1}$), wp) - real(pk(j, &
                                            & k, l), wp))
                            uf(j, k, l, ${D}$) = uf(j, k, l, ${D}$) - gf
                            pflx(j, k, l) = ga - gf
                            if (wbl) pflx(j, k, l) = pflx(j, k, l) + tau*f_capillary_accel(kap(j, k, l, 1), kap(${IP1}$, 1), &
                                & kap(j, k, l, 2), kap(${IP1}$, 2), real(q_cons_vf(eqn_idx%c)%sf(j, k, l), wp), &
                                & real(q_cons_vf(eqn_idx%c)%sf(${IP1}$), wp), rhoc(j, k, l), rhoc(${IP1}$), &
                                & 0.5_wp*(${DXV}$(${SV}$) + ${DXV}$(${SV}$ + 1)))
                            if ((${SV}$ == -1 .and. wlo) .or. (${SV}$ == ${UB}$ .and. whi)) pflx(j, k, l) = 0._wp
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(j, k, &
                                      & l) = real(real(q_cons_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(j, k, l), wp) + rhoc(j, k, &
                                      & l)*0.5_wp*(pflx(${IM1}$) + pflx(j, k, l)), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor
        call s_zero_wall_faces()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho, gam, pinf, qv, ke, ar, al]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ar(i) = real(q_cons_vf(i)%sf(j, k, l), wp)
                        al(i) = real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                    end do
                    call s_compute_mixture_coefficients(ar, al, rho, gam, pinf, qv)
                    ke = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        ke = ke + 0.5_wp*real(q_cons_vf(i)%sf(j, k, l), wp)*(real(q_cons_vf(i)%sf(j, k, l), wp)/rho)
                    end do
                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = real(gam*real(pk(j, k, l), wp) + pinf + qv + ke, stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_projection_apply

    !> PCG on the SPD pressure system, preconditioned by one multigrid V-cycle. The solution is left in pk with filled ghosts.
    impure subroutine s_pcg_solve(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp)                                                   :: bnorm, rnorm, rtol, rz, rz_new, dq, alpha, beta
        integer                                                    :: it, j, k, l

        call s_mg_build()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    pk(j, k, l) = real(xs(j, k, l), stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_apply_operator(bc_type)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rs(j, k, l) = bvec(j, k, l) - qs(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Relative to the initial residual, i.e. to the pressure change being solved for; relative to b the tolerance would
        ! depend on the ambient pressure, which b carries in full
        bnorm = sqrt(f_dot(bvec, bvec))
        rnorm = sqrt(f_dot(rs, rs))
        rtol = max(proj_tol*rnorm, res_floor*bnorm)

        if (rnorm > rtol) then
            call s_mg_vcycle()
            call s_copy_to_pk(zs)
            rz = f_dot(rs, zs)

            do it = 1, proj_max_iters
                call s_apply_operator(bc_type)
                dq = f_dot_pk(qs)
                alpha = rz/dq
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            xs(j, k, l) = xs(j, k, l) + alpha*real(pk(j, k, l), wp)
                            rs(j, k, l) = rs(j, k, l) - alpha*qs(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                rnorm = sqrt(f_dot(rs, rs))
                if (rnorm <= rtol) exit

                call s_mg_vcycle()
                rz_new = f_dot(rs, zs)
                beta = rz_new/rz
                rz = rz_new
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            pk(j, k, l) = real(zs(j, k, l) + beta*real(pk(j, k, l), wp), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end if
        call s_copy_to_pk(xs)
        call s_populate_F_igr_buffers(bc_type, pk_sf)

    end subroutine s_pcg_solve

    subroutine s_copy_to_pk(v)

        real(wp), dimension(0:m,0:n,0:p), intent(in) :: v
        integer                                      :: j, k, l

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    pk(j, k, l) = real(v(j, k, l), stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_copy_to_pk

    !> qs = A pk, with the ghosts of pk filled first so periodic, MPI and wall neighbors all enter through the same stencil (a
    !! wall's even reflection makes its face term vanish, which is the Neumann condition)
    impure subroutine s_apply_operator(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp)                                                   :: s, pc, area
        integer                                                    :: j, k, l

        call s_populate_F_igr_buffers(bc_type, pk_sf)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, s, pc, area]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    pc = real(pk(j, k, l), wp)
                    s = dcoef(j, k, l)*pc
                    area = 1._wp
                    if (num_dims > 1) area = area*dy(k)
                    if (num_dims > 2) area = area*dz(l)
                    s = s + f_cond(rhoc(j, k, l), rhoc(j - 1, k, l), area, 0.5_wp*(dx(j - 1) + dx(j)))*(pc - real(pk(j - 1, k, &
                                   & l), wp)) + f_cond(rhoc(j, k, l), rhoc(j + 1, k, l), area, &
                                   & 0.5_wp*(dx(j) + dx(j + 1)))*(pc - real(pk(j + 1, k, l), wp))
                    if (num_dims > 1) then
                        area = dx(j)
                        if (num_dims > 2) area = area*dz(l)
                        s = s + f_cond(rhoc(j, k, l), rhoc(j, k - 1, l), area, 0.5_wp*(dy(k - 1) + dy(k)))*(pc - real(pk(j, &
                                       & k - 1, l), wp)) + f_cond(rhoc(j, k, l), rhoc(j, k + 1, l), area, &
                                       & 0.5_wp*(dy(k) + dy(k + 1)))*(pc - real(pk(j, k + 1, l), wp))
                    end if
                    if (num_dims > 2) then
                        area = dx(j)*dy(k)
                        s = s + f_cond(rhoc(j, k, l), rhoc(j, k, l - 1), area, 0.5_wp*(dz(l - 1) + dz(l)))*(pc - real(pk(j, k, &
                                       & l - 1), wp)) + f_cond(rhoc(j, k, l), rhoc(j, k, l + 1), area, &
                                       & 0.5_wp*(dz(l) + dz(l + 1)))*(pc - real(pk(j, k, l + 1), wp))
                    end if
                    qs(j, k, l) = s
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_operator

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
        call s_mpi_allreduce_sum(loc, res)

    end function f_dot

    !> Global inner product of the search direction with a vector
    impure function f_dot_pk(a) result(res)

        real(wp), dimension(0:m,0:n,0:p), intent(in) :: a
        real(wp)                                     :: res, loc
        integer                                      :: j, k, l

        loc = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[loc]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    loc = loc + real(pk(j, k, l), wp)*a(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_mpi_allreduce_sum(loc, res)

    end function f_dot_pk

    !> Level 1 from the fine system, with this rank's boundary faces dropped, then Galerkin coarsening for piecewise-constant
    !! aggregation: a coarse diagonal sums its children, a coarse face sums the fine faces lying on it. Exact for any coefficient
    !! jump, which rediscretizing an averaged density would not be.
    impure subroutine s_mg_build()

        integer  :: lv, nx, ny, nz, off, cnx, cny, cnz, coff, sx, sy, sz, ii, jj, kk, a, b, c, idx, cidx
        real(wp) :: area, sd, skx, sky, skz

        nx = mg_nx(1); ny = mg_ny(1)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, area]')
        do kk = 0, p
            do jj = 0, n
                do ii = 0, m
                    idx = (kk*ny + jj)*nx + ii + 1
                    mg_d(idx) = dcoef(ii, jj, kk)
                    area = 1._wp
                    if (num_dims > 1) area = dy(jj)
                    if (num_dims > 2) area = area*dz(kk)
                    mg_kx(idx) = merge(f_cond(rhoc(ii, jj, kk), rhoc(ii - 1, jj, kk), area, 0.5_wp*(dx(ii - 1) + dx(ii))), 0._wp, &
                          & ii > 0)
                    mg_ky(idx) = 0._wp
                    mg_kz(idx) = 0._wp
                    if (num_dims > 1 .and. jj > 0) then
                        area = dx(ii)
                        if (num_dims > 2) area = area*dz(kk)
                        mg_ky(idx) = f_cond(rhoc(ii, jj, kk), rhoc(ii, jj - 1, kk), area, 0.5_wp*(dy(jj - 1) + dy(jj)))
                    end if
                    if (num_dims > 2 .and. kk > 0) then
                        mg_kz(idx) = f_cond(rhoc(ii, jj, kk), rhoc(ii, jj, kk - 1), dx(ii)*dy(jj), 0.5_wp*(dz(kk - 1) + dz(kk)))
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        do lv = 1, mg_nlev - 1
            nx = mg_nx(lv); ny = mg_ny(lv); off = mg_off(lv)
            cnx = mg_nx(lv + 1); cny = mg_ny(lv + 1); cnz = mg_nz(lv + 1); coff = mg_off(lv + 1)
            sx = mg_sx(lv); sy = mg_sy(lv); sz = mg_sz(lv)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, a, b, c, idx, cidx, sd, skx, sky, skz]')
            do kk = 0, cnz - 1
                do jj = 0, cny - 1
                    do ii = 0, cnx - 1
                        sd = 0._wp; skx = 0._wp; sky = 0._wp; skz = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do c = 0, sz - 1
                            $:GPU_LOOP(parallelism='[seq]')
                            do b = 0, sy - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do a = 0, sx - 1
                                    idx = off + ((sz*kk + c)*ny + sy*jj + b)*nx + sx*ii + a + 1
                                    sd = sd + mg_d(idx)
                                    if (a == 0) skx = skx + mg_kx(idx)
                                    if (b == 0) sky = sky + mg_ky(idx)
                                    if (c == 0) skz = skz + mg_kz(idx)
                                end do
                            end do
                        end do
                        cidx = coff + (kk*cny + jj)*cnx + ii + 1
                        mg_d(cidx) = sd; mg_kx(cidx) = skx; mg_ky(cidx) = sky; mg_kz(cidx) = skz
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_mg_build

    !> One symmetric V-cycle on rs, returned in zs. Pre-smoothing runs red then black and post-smoothing the reverse, which is what
    !! makes the cycle a symmetric operator and so a valid CG preconditioner.
    impure subroutine s_mg_vcycle()

        integer :: lv, j, k, l, idx, nx, ny, i

        nx = mg_nx(1); ny = mg_ny(1)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, idx]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    idx = (l*ny + k)*nx + j + 1
                    mg_f(idx) = rs(j, k, l)
                    mg_e(idx) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        do lv = 1, mg_nlev - 1
            do i = 1, mg_nu
                call s_mg_smooth(lv, 0); call s_mg_smooth(lv, 1)
            end do
            call s_mg_restrict(lv)
        end do
        do i = 1, max(mg_nx(mg_nlev), mg_ny(mg_nlev), mg_nz(mg_nlev))
            call s_mg_smooth(mg_nlev, 0); call s_mg_smooth(mg_nlev, 1)
            call s_mg_smooth(mg_nlev, 1); call s_mg_smooth(mg_nlev, 0)
        end do
        do lv = mg_nlev - 1, 1, -1
            call s_mg_prolong(lv)
            do i = 1, mg_nu
                call s_mg_smooth(lv, 1); call s_mg_smooth(lv, 0)
            end do
        end do

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, idx]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    idx = (l*ny + k)*nx + j + 1
                    zs(j, k, l) = mg_e(idx)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_vcycle

    !> Gauss-Seidel update of one color on one level
    impure subroutine s_mg_smooth(lv, color)

        integer, intent(in) :: lv, color
        integer             :: nx, ny, nz, off, ii, jj, kk, idx
        real(wp)            :: dg, nb, kf

        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); off = mg_off(lv)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb, kf]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    if (mod(ii + jj + kk, 2) == color) then
                        idx = off + (kk*ny + jj)*nx + ii + 1
                        @:MG_ROW()
                        mg_e(idx) = (mg_f(idx) + nb)/dg
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_smooth

    !> Residual of level lv summed onto the coarse right-hand side (restriction is the transpose of the prolongation)
    impure subroutine s_mg_restrict(lv)

        integer, intent(in) :: lv
        integer             :: nx, ny, nz, off, cnx, cny, cnz, coff, sx, sy, sz, ii, jj, kk, idx, cidx, a, b, c
        real(wp)            :: dg, nb, kf

        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); off = mg_off(lv)
        cnx = mg_nx(lv + 1); cny = mg_ny(lv + 1); cnz = mg_nz(lv + 1); coff = mg_off(lv + 1)
        sx = mg_sx(lv); sy = mg_sy(lv); sz = mg_sz(lv)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb, kf]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    idx = off + (kk*ny + jj)*nx + ii + 1
                    @:MG_ROW()
                    mg_r(idx) = mg_f(idx) - (dg*mg_e(idx) - nb)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Gather per coarse cell rather than scatter with atomics, so the sum order is fixed and results are reproducible
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, a, b, c, cidx, dg]')
        do kk = 0, cnz - 1
            do jj = 0, cny - 1
                do ii = 0, cnx - 1
                    dg = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do c = 0, sz - 1
                        $:GPU_LOOP(parallelism='[seq]')
                        do b = 0, sy - 1
                            $:GPU_LOOP(parallelism='[seq]')
                            do a = 0, sx - 1
                                dg = dg + mg_r(off + ((sz*kk + c)*ny + sy*jj + b)*nx + sx*ii + a + 1)
                            end do
                        end do
                    end do
                    cidx = coff + (kk*cny + jj)*cnx + ii + 1
                    mg_f(cidx) = dg
                    mg_e(cidx) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_restrict

    !> Add each coarse correction to its children
    impure subroutine s_mg_prolong(lv)

        integer, intent(in) :: lv
        integer             :: nx, ny, nz, off, cnx, cny, coff, sx, sy, sz, ii, jj, kk, idx

        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); off = mg_off(lv)
        cnx = mg_nx(lv + 1); cny = mg_ny(lv + 1); coff = mg_off(lv + 1)
        sx = mg_sx(lv); sy = mg_sy(lv); sz = mg_sz(lv)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    idx = off + (kk*ny + jj)*nx + ii + 1
                    mg_e(idx) = mg_e(idx) + mg_e(coff + ((kk/sz)*cny + jj/sy)*cnx + ii/sx + 1)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_prolong

    impure subroutine s_finalize_projection_module()

        $:GPU_EXIT_DATA(detach='[pk_sf(1)%sf]')
        @:DEALLOCATE(uf, divu, rhs_p, p_stage, p_step0, pflx, rhoc, dcoef, bvec, xs, rs, zs, qs, pk, kap)
        @:DEALLOCATE(mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r)

    end subroutine s_finalize_projection_module

end module m_projection
