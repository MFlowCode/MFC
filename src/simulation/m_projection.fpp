!>
!! @file
!! @brief Contains module m_projection

#:include 'case.fpp'
#:include 'macros.fpp'

!> All-Mach pressure projection (Fuster & Popinet, JCP 374, 2018). Advection uses a persistent face velocity; the pressure then
!! solves p - rho*c^2*tau^2*div(rho_f^-1 grad p) = p_adv - rho*c^2*tau*div(u*_f) with div, grad and the Laplacian all taken on
!! faces, so they compose exactly and the projected face velocity satisfies the discrete pressure equation. The solve is PCG
!! preconditioned by geometric multigrid (m_multigrid).
!> @brief All-Mach pressure projection
module m_projection

#ifdef MFC_MPI
    use mpi
#endif
    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_boundary_common
    use m_eos
    use m_body_forces, only: s_compute_acceleration
    use m_riemann_state, only: Re_avg_rsx_vf, vel_src_rsx_vf, s_compute_interface_reynolds
    use m_projection_helpers
    use m_multigrid, only: s_initialize_multigrid_module, s_mg_set_operator, s_mg_solve, s_finalize_multigrid_module
    use m_nvtx

    implicit none

    private; public :: s_initialize_projection_module, s_projection_rhs, s_projection_face_props, s_projection_heat, &
        & s_projection_apply, s_finalize_projection_module

    real(wp), allocatable, dimension(:,:,:,:) :: uf           !< face velocity; index j is the face between cells j and j+1
    real(wp), allocatable, dimension(:,:,:,:) :: uf0          !< the face velocity that transported the current stage
    real(wp), allocatable, dimension(:,:,:)   :: divu, rhs_p  !< div of uf, and the pressure transport rate
    real(wp), allocatable, dimension(:,:,:)   :: pflx         !< upwind pressure flux on the faces of one direction
    real(wp), allocatable, dimension(:,:,:)   :: p_stage, p_step0
    real(wp), allocatable, dimension(:,:,:)   :: rhoc         !< star density with one ghost layer
    real(wp), allocatable, dimension(:,:,:)   :: dcoef, bvec  !< SPD system: D_c and right-hand side
    real(wp), allocatable, dimension(:,:,:)   :: xs           !< pressure-solve initial guess, then solution
    !> Face conductances A_f/(rho_f d_f) of the pressure operator, x, y and z faces; index j is the face between cells j - 1 and j
    real(wp), allocatable, dimension(:,:,:)          :: cfx, cfy, cfz
    real(stp), allocatable, dimension(:,:,:), target :: pk  !< search direction, then solution, with ghosts for the halo
    type(scalar_field), dimension(1)                 :: pk_sf
    $:GPU_DECLARE(create='[pk_sf]')
    $:GPU_DECLARE(create='[uf, uf0, divu, rhs_p, pflx, p_stage, p_step0, rhoc, dcoef, bvec, xs, cfx, cfy, cfz, pk]')
    !> Well-balanced surface tension, with one ghost layer: curvature (1) and |grad c| (2), both zero outside the interface band
    real(wp), allocatable, dimension(:,:,:,:) :: kap
    real(wp), allocatable, dimension(:,:,:,:) :: gnd  !< grad(alpha_1) (1:num_dims) and its magnitude (0), two ghost layers
    $:GPU_DECLARE(create='[kap, gnd]')

    !> Boundary code of each domain face, low (1) and high (2) side, copied from bc_type: walls and boundary patches per face
    integer, allocatable, dimension(:,:,:) :: fbx, fby, fbz
    $:GPU_DECLARE(create='[fbx, fby, fbz]')
    !> An extrapolation (-3) face is a pressure outlet held at bc_[x,y,z]%pres_out where that is set: its conductance couples the
    !! boundary cell to a ghost at pout, so flow leaves where the solve would otherwise see a wall
    logical, dimension(3)  :: has_pout
    real(wp), dimension(3) :: pout
    logical, dimension(3)  :: seam_lo, seam_hi    !< that side couples to another rank or, periodically, to this one
    logical                :: faces_ready         !< uf has been seeded from the cell velocities
    logical                :: wb_st               !< well-balanced surface tension
    integer                :: gk0, gk1, gl0, gl1  !< y and z extents including one ghost layer where those directions exist
    !> Loop bounds by direction d: face_bounds(:, d) spans the d faces (-1 to m, n, p along d); plane_bounds(:, d, ip) the plane
    !! normal to d at index -1 (ip = -1), 0 (ip = 0), the last cell (ip = 1) or the ghost past it (ip = 2)
    type(int_bounds_info), dimension(3, 3)      :: face_bounds
    type(int_bounds_info), dimension(3, 3,-1:2) :: plane_bounds
    $:GPU_DECLARE(create='[face_bounds, plane_bounds]')

contains

    impure subroutine s_initialize_projection_module()

        integer :: d

#ifdef MFC_MIXED_PRECISION
        call s_mpi_abort('proj_method needs stp = wp; mixed precision is not supported')
#endif

        @:ALLOCATE(uf(-1:m + 1, -1:n + 1, -1:p + 1, 1:num_dims), uf0(-1:m + 1, -1:n + 1, -1:p + 1, 1:num_dims))
        @:ALLOCATE(divu(0:m, 0:n, 0:p), rhs_p(0:m, 0:n, 0:p), p_stage(0:m, 0:n, 0:p), p_step0(0:m, 0:n, 0:p))
        @:ALLOCATE(pflx(-1:m + 1, -1:n + 1, -1:p + 1), rhoc(-1:m + 1, -1:n + 1, -1:p + 1))
        @:ALLOCATE(dcoef(0:m, 0:n, 0:p), bvec(0:m, 0:n, 0:p))
        @:ALLOCATE(xs(0:m, 0:n, 0:p), cfx(0:m + 1, 0:n, 0:p), cfy(0:m, 0:merge(n + 1, 0, n > 0), 0:p))
        @:ALLOCATE(cfz(0:m, 0:n, 0:merge(p + 1, 0, p > 0)))
        @:ALLOCATE(pk(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))

        pk_sf(1)%sf => pk
        $:GPU_ENTER_DATA(copyin='[pk_sf(1)%sf]')
        $:GPU_ENTER_DATA(attach='[pk_sf(1)%sf]')

        #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            seam_lo(${D}$) = bc_${XYZ}$%beg >= 0 .or. bc_${XYZ}$%beg == BC_PERIODIC
            seam_hi(${D}$) = bc_${XYZ}$%end >= 0 .or. bc_${XYZ}$%end == BC_PERIODIC
        #:endfor
        call s_initialize_multigrid_module()

        #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            has_pout(${D}$) = .not. f_is_default(bc_${XYZ}$%pres_out) .and. num_dims >= ${D}$
            pout(${D}$) = merge(bc_${XYZ}$%pres_out, 0._wp, has_pout(${D}$))
        #:endfor
        @:ALLOCATE(fbx(0:n, 0:p, 1:2), fby(0:merge(m, 0, n > 0), 0:p, 1:2), fbz(0:merge(m, 0, p > 0), 0:merge(n, 0, p > 0), 1:2))
        faces_ready = .false.
        do d = 1, num_dims
            call s_bounds(d, -1, f_last(d), face_bounds(:,d))
            call s_bounds(d, -1, -1, plane_bounds(:,d, -1))
            call s_bounds(d, 0, 0, plane_bounds(:,d, 0))
            call s_bounds(d, f_last(d), f_last(d), plane_bounds(:,d, 1))
            call s_bounds(d, f_last(d) + 1, f_last(d) + 1, plane_bounds(:,d, 2))
        end do
        $:GPU_UPDATE(device='[face_bounds, plane_bounds]')
        gk0 = merge(-1, 0, n > 0); gk1 = merge(n + 1, n, n > 0)
        gl0 = merge(-1, 0, p > 0); gl1 = merge(p + 1, p, p > 0)

        wb_st = surface_tension .and. surface_tension_model == surface_tension_model_well_balanced
        if (wb_st) then
            @:ALLOCATE(kap(-1:m + 1, gk0:gk1, gl0:gl1, 1:2))
            @:ALLOCATE(gnd(-2:m + 2, 2*gk0:gk1 - gk0, 2*gl0:gl1 - gl0, 0:3))
        else
            @:ALLOCATE(kap(0:0, 0:0, 0:0, 1:2), gnd(0:0, 0:0, 0:0, 0:3))
        end if

    end subroutine s_initialize_projection_module

    !> Boundary code of the side s (1 low, 2 high) domain face of direction d that (j, k, l) lies across from
    integer function f_bc_code(d, j, k, l, s)

        $:GPU_ROUTINE(function_name='f_bc_code', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: d, j, k, l, s

        if (d == 1) then
            f_bc_code = fbx(k, l, s)
        else if (d == 2) then
            f_bc_code = fby(j, l, s)
        else
            f_bc_code = fbz(j, k, s)
        end if

    end function f_bc_code

    !> Conductance of the direction-d face at index (j, k, l) (cfx, cfy, cfz)
    real(wp) function f_cf(d, j, k, l)
        $:GPU_ROUTINE(function_name='f_cf', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: d, j, k, l

        if (d == 1) then
            f_cf = cfx(j, k, l)
        else if (d == 2) then
            f_cf = cfy(j, k, l)
        else
            f_cf = cfz(j, k, l)
        end if

    end function f_cf

    !> Copy the boundary code of each domain face from bc_type, on the host where pre-process left it
    impure subroutine s_projection_face_codes(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer                                                    :: e

        do e = 1, 2
            fbx(:,:,e) = int(bc_type(1, e)%sf(0,:,:))
            if (num_dims > 1) fby(:,:,e) = int(bc_type(2, e)%sf(0:m,0,:))
            if (num_dims > 2) fbz(:,:,e) = int(bc_type(3, e)%sf(0:m,0:n,0))
        end do
        $:GPU_UPDATE(device='[fbx, fby, fbz]')

    end subroutine s_projection_face_codes

    !> Face velocity from the average of the cell velocities either side (primitive velocities, with ghosts)
    subroutine s_projection_init_faces(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: d, iv, j, k, l, oj, ok, ol

        do d = 1, num_dims
            call s_unit(d, oj, ok, ol)
            iv = eqn_idx%mom%beg + d - 1
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = face_bounds(3, d)%beg, face_bounds(3, d)%end
                do k = face_bounds(2, d)%beg, face_bounds(2, d)%end
                    do j = face_bounds(1, d)%beg, face_bounds(1, d)%end
                        uf(j, k, l, d) = 0.5_wp*(real(q_prim_vf(iv)%sf(j, k, l), wp) + real(q_prim_vf(iv)%sf(j + oj, k + ok, &
                           & l + ol), wp))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do
        call s_zero_wall_faces()

    end subroutine s_projection_init_faces

    !> A Dirichlet (-17) face, a whole boundary or a patch such as a nozzle set into a wall, carries the normal velocity of its
    !! ghost cell, which s_dirichlet sets from the boundary buffer and its inflow ramp: prescribed, it is neither predicted nor
    !! pressure-corrected (the pressure ghost is a zero-gradient copy) and the solve sees it as a source of volume
    subroutine s_dirichlet_faces(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: d, s, g, ip, iv, j, k, l, oj, ok, ol

        ! The low face and its ghost cell share index -1; the high face is m (n, p), its ghost (g = 1) one further

        do d = 1, num_dims
            call s_unit(d, oj, ok, ol)
            iv = eqn_idx%mom%beg + d - 1
            do s = 1, 2
                g = s - 1
                ip = merge(-1, 1, s == 1)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = plane_bounds(3, d, ip)%beg, plane_bounds(3, d, ip)%end
                    do k = plane_bounds(2, d, ip)%beg, plane_bounds(2, d, ip)%end
                        do j = plane_bounds(1, d, ip)%beg, plane_bounds(1, d, ip)%end
                            if (f_bc_code(d, j, k, l, s) == BC_DIRICHLET) uf(j, k, l, d) = real(q_prim_vf(iv)%sf(j + g*oj, &
                                & k + g*ok, l + g*ol), wp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

    end subroutine s_dirichlet_faces

    !> Normal velocity vanishes on solid walls
    subroutine s_zero_wall_faces()

        integer :: d, s, ip, j, k, l

        do d = 1, num_dims
            do s = 1, 2
                ip = merge(-1, 1, s == 1)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = plane_bounds(3, d, ip)%beg, plane_bounds(3, d, ip)%end
                    do k = plane_bounds(2, d, ip)%beg, plane_bounds(2, d, ip)%end
                        do j = plane_bounds(1, d, ip)%beg, plane_bounds(1, d, ip)%end
                            if (f_wall_face(f_bc_code(d, j, k, l, s))) uf(j, k, l, d) = 0._wp
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

    end subroutine s_zero_wall_faces

    !> Advective right-hand side of one direction sweep. Every quantity is carried by the projected face velocity and upwinded on
    !! its sign; the momentum flux is the summed partial-density flux times the upwind velocity, so mass and momentum move with one
    !! operator. Energy is left at zero here and rebuilt from the equation of state after the pressure solve.
    subroutine s_projection_rhs(id, qfl_rs, qfr_rs, q_prim_vf, flux_vf, rhs_vf, bc_type)

        integer, intent(in)                                                                 :: id
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: qfl_rs, qfr_rs
        type(scalar_field), dimension(sys_size), intent(in)                                 :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout)                              :: flux_vf, rhs_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                          :: bc_type
        real(wp)                                                                            :: vf, a_up, ar_up, fm
        logical                                                                             :: up_l
        integer                                                                             :: i, j, k, l, d, oj, ok, ol

        if (id == 1) then
            if (.not. faces_ready) then
                call s_projection_face_codes(bc_type)
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
            call s_dirichlet_faces(q_prim_vf)
        end if

        d = id
        call s_unit(d, oj, ok, ol)
        ! Face fluxes. Left state of face j is the right edge of cell j, right state the left edge of cell j+1
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, vf, a_up, ar_up, fm, up_l]')
        do l = face_bounds(3, d)%beg, face_bounds(3, d)%end
            do k = face_bounds(2, d)%beg, face_bounds(2, d)%end
                do j = face_bounds(1, d)%beg, face_bounds(1, d)%end
                    vf = uf(j, k, l, d)
                    up_l = vf >= 0._wp
                    fm = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ! Reconstructed partial densities upwinded directly. Taking the phase density alpha_rho/alpha from the
                        ! upwind cell instead diverges where a phase is vanishing: both are round-off there, and their ratio
                        ! (seen at 6e10) times the face alpha flux injects mass
                        if (up_l) then
                            a_up = qfl_rs(j, k, l, eqn_idx%adv%beg + i - 1)
                            ar_up = qfl_rs(j, k, l, i)
                        else
                            a_up = qfr_rs(j + oj, k + ok, l + ol, eqn_idx%adv%beg + i - 1)
                            ar_up = qfr_rs(j + oj, k + ok, l + ol, i)
                        end if
                        flux_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l) = real(a_up*vf, stp)
                        flux_vf(i)%sf(j, k, l) = real(ar_up*vf, stp)
                        fm = fm + ar_up*vf
                    end do
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_dims
                        if (up_l) then
                            flux_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = real(fm*qfl_rs(j, k, l, eqn_idx%mom%beg + i - 1), stp)
                        else
                            flux_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = real(fm*qfr_rs(j + oj, k + ok, l + ol, &
                                    & eqn_idx%mom%beg + i - 1), stp)
                        end if
                    end do
                    if (up_l) then
                        pflx(j, k, l) = vf*qfl_rs(j, k, l, eqn_idx%E)
                    else
                        pflx(j, k, l) = vf*qfr_rs(j + oj, k + ok, l + ol, eqn_idx%E)
                    end if
                    ! The color function is not reconstructed; upwind its cell value
                    if (surface_tension) then
                        if (up_l) then
                            flux_vf(eqn_idx%c)%sf(j, k, l) = real(vf*real(q_prim_vf(eqn_idx%c)%sf(j, k, l), wp), stp)
                        else
                            flux_vf(eqn_idx%c)%sf(j, k, l) = real(vf*real(q_prim_vf(eqn_idx%c)%sf(j + oj, k + ok, l + ol), wp), stp)
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
                            rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + real((real(flux_vf(i)%sf(j - oj, k - ok, l - ol), &
                                   & wp) - real(flux_vf(i)%sf(j, k, l), wp))/f_width(d, j, k, l), stp)
                        end if
                    end do
                    rhs_p(j, k, l) = rhs_p(j, k, l) + (pflx(j - oj, k - ok, l - ol) - pflx(j, k, l))/f_width(d, j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

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
        integer                :: i, j, k, l, rs1, rs2, d, oj, ok, ol

        rs1 = Re_size(1); rs2 = Re_size(2)

        d = id
        call s_unit(d, oj, ok, ol)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, al, ar, re_l, re_r]', firstprivate='[rs1, rs2]')
        do l = face_bounds(3, d)%beg, face_bounds(3, d)%end
            do k = face_bounds(2, d)%beg, face_bounds(2, d)%end
                do j = face_bounds(1, d)%beg, face_bounds(1, d)%end
                    if (viscous) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            al(i) = qfl_rs(j, k, l, eqn_idx%adv%beg + i - 1)
                            ar(i) = qfr_rs(j + oj, k + ok, l + ol, eqn_idx%adv%beg + i - 1)
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
                        vel_src_rsx_vf(j, k, l, i) = 0.5_wp*(qfl_rs(j, k, l, eqn_idx%mom%beg + i - 1) + qfr_rs(j + oj, k + ok, &
                                       & l + ol, eqn_idx%mom%beg + i - 1))
                    end do
                    if (surface_tension) flux_src_vf(eqn_idx%adv%beg)%sf(j, k, l) = real(uf(j, k, l, d), stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_projection_face_props

    !> Conduction heating of one direction sweep from the energy source flux fe (-k dT/dn on faces). At fixed alpha and phase
    !! densities d(rho e) = Gamma dp, so the pressure transport rate gains -div(q)/Gamma, and the energy rebuilt from p carries it
    subroutine s_projection_heat(id, q_prim_vf, fe)

        integer, intent(in)                                 :: id
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), intent(in)                      :: fe
        real(wp)                                            :: rho, gam, pinf, qv

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: ar, al
        #:else
            real(wp), dimension(num_fluids) :: ar, al
        #:endif
        integer :: i, j, k, l, d, oj, ok, ol

        d = id
        call s_unit(d, oj, ok, ol)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho, gam, pinf, qv, ar, al]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ar(i) = real(q_prim_vf(i)%sf(j, k, l), wp)
                        al(i) = real(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                    end do
                    call s_compute_mixture_coefficients(ar, al, rho, gam, pinf, qv)
                    rhs_p(j, k, l) = rhs_p(j, k, l) + (real(fe%sf(j - oj, k - ok, l - ol), wp) - real(fe%sf(j, k, l), &
                          & wp))/(f_width(d, j, k, l)*gam)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_projection_heat

    !> Divergence of the face velocity in cell (j, k, l): the one operator the transport sources and the pressure equation share
    function f_div_uf(j, k, l) result(dv)

        $:GPU_ROUTINE(function_name='f_div_uf', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: j, k, l
        real(wp)            :: dv

        dv = (uf(j, k, l, 1) - uf(j - 1, k, l, 1))/dx(j)
        if (num_dims > 1) dv = dv + (uf(j, k, l, 2) - uf(j, k - 1, l, 2))/dy(k)
        if (num_dims > 2) dv = dv + (uf(j, k, l, 3) - uf(j, k, l - 1, 3))/dz(l)

    end function f_div_uf

    !> Capillary acceleration of the direction-d face between (j, k, l) and the next cell, (oj, ok, ol) further on
    function f_face_capillary(d, j, k, l, oj, ok, ol, q_cons_vf) result(acc)

        $:GPU_ROUTINE(function_name='f_face_capillary', parallelism='[seq]', cray_inline=True)

        integer, intent(in)                                 :: d, j, k, l, oj, ok, ol
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp)                                            :: acc

        acc = f_capillary_accel(kap(j, k, l, 1), kap(j + oj, k + ok, l + ol, 1), kap(j, k, l, 2), kap(j + oj, k + ok, l + ol, 2), &
                                & real(q_cons_vf(eqn_idx%adv%beg)%sf(j, k, l), wp), real(q_cons_vf(eqn_idx%adv%beg)%sf(j + oj, &
                                & k + ok, l + ol), wp), rhoc(j, k, l), rhoc(j + oj, k + ok, l + ol), f_face_dist(d, j, k, l))

    end function f_face_capillary

    !> Interface normal component d at cell (j, k, l), from gnd
    function f_nrm(d, j, k, l) result(nd)

        $:GPU_ROUTINE(function_name='f_nrm', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: d, j, k, l
        real(wp)            :: nd

        nd = f_normal(gnd(j, k, l, d), gnd(j, k, l, 0))

    end function f_nrm

    !> Curvature kappa = -div(n) in the interface band, with n = grad(alpha_1)/|grad(alpha_1)|. The CSF force sigma*kappa*grad(c) is
    !! unchanged under c -> 1 - c, so the conservatively transported, compression-sharpened volume fraction serves as the indicator;
    !! the color function, upwinded for the stress-tensor model, smears and loses its interior value over time
    subroutine s_compute_curvature(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp)                                            :: kv
        integer                                             :: j, k, l, k0, k1, l0, l1, dk, dl, ia

        k0 = gk0; k1 = gk1; l0 = gl0; l1 = gl1
        dk = -k0; dl = -l0
        ia = eqn_idx%adv%beg
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = l0 - dl, l1 + dl
            do k = k0 - dk, k1 + dk
                do j = -2, m + 2
                    gnd(j, k, l, 1) = real(q_cons_vf(ia)%sf(j + 1, k, l) - q_cons_vf(ia)%sf(j - 1, k, l), &
                        & wp)/(x_cc(j + 1) - x_cc(j - 1))
                    gnd(j, k, l, 2) = 0._wp; gnd(j, k, l, 3) = 0._wp
                    if (num_dims > 1) gnd(j, k, l, 2) = real(q_cons_vf(ia)%sf(j, k + 1, l) - q_cons_vf(ia)%sf(j, k - 1, l), &
                        & wp)/(y_cc(k + 1) - y_cc(k - 1))
                    if (num_dims > 2) gnd(j, k, l, 3) = real(q_cons_vf(ia)%sf(j, k, l + 1) - q_cons_vf(ia)%sf(j, k, l - 1), &
                        & wp)/(z_cc(l + 1) - z_cc(l - 1))
                    gnd(j, k, l, 0) = sqrt(gnd(j, k, l, 1)**2 + gnd(j, k, l, 2)**2 + gnd(j, k, l, 3)**2)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, kv]')
        do l = l0, l1
            do k = k0, k1
                do j = -1, m + 1
                    kap(j, k, l, 1) = 0._wp
                    kap(j, k, l, 2) = 0._wp
                    if (gnd(j, k, l, 0) > capillary_cutoff) then
                        kv = -(f_nrm(1, j + 1, k, l) - f_nrm(1, j - 1, k, l))/(x_cc(j + 1) - x_cc(j - 1))
                        if (num_dims > 1) kv = kv - (f_nrm(2, j, k + 1, l) - f_nrm(2, j, k - 1, l))/(y_cc(k + 1) - y_cc(k - 1))
                        if (num_dims > 2) kv = kv - (f_nrm(3, j, k, l + 1) - f_nrm(3, j, k, l - 1))/(z_cc(l + 1) - z_cc(l - 1))
                        kap(j, k, l, 1) = kv
                        kap(j, k, l, 2) = gnd(j, k, l, 0)
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
        logical :: wbl, fix

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: ar, al
        #:else
            real(wp), dimension(num_fluids) :: ar, al
        #:endif
        integer :: i, j, k, l, k0, k1, l0, l1, d, s, ia, iv, oj, ok, ol

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

        ! Keep the face velocity this stage was transported with, for s_mass_lag below
        if (proj_lag_corr > 0) then
            $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
            do i = 1, num_dims
                do l = l0, l1
                    do k = k0, k1
                        do j = -1, m + 1
                            uf0(j, k, l, i) = uf(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
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

        if (wb_st) call s_compute_curvature(q_cons_vf)
        wbl = wb_st

        ! Body forces and well-balanced surface tension enter on faces, where they meet the pressure gradient that balances them.
        ! A Dirichlet face keeps its prescribed velocity
        do d = 1, num_dims
            call s_unit(d, oj, ok, ol)
            iv = eqn_idx%mom%beg + d - 1
            ga = tau*acc(d)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, fix, ia]')
            do l = face_bounds(3, d)%beg, face_bounds(3, d)%end
                do k = face_bounds(2, d)%beg, face_bounds(2, d)%end
                    do j = face_bounds(1, d)%beg, face_bounds(1, d)%end
                        ia = f_along(d, j, k, l)
                        fix = .false.
                        if (ia == -1) fix = f_bc_code(d, j, k, l, 1) == BC_DIRICHLET
                        if (ia == f_last(d)) fix = f_bc_code(d, j, k, l, 2) == BC_DIRICHLET
                        if (.not. fix) then
                            uf(j, k, l, d) = 0.5_wp*(real(q_cons_vf(iv)%sf(j, k, l), wp)/rhoc(j, k, &
                               & l) + real(q_cons_vf(iv)%sf(j + oj, k + ok, l + ol), wp)/rhoc(j + oj, k + ok, l + ol)) + ga
                            if (wbl) uf(j, k, l, d) = uf(j, k, l, d) + tau*f_face_capillary(d, j, k, l, oj, ok, ol, q_cons_vf)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do
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

        if (stage == 1) then
            proj_pcg_iters = 0
            proj_pcg_res = 0._wp
        end if
        call nvtxStartRange("PROJ-PCG")
        call s_pressure_solve(bc_type)
        call nvtxEndRange

        if (proj_lag_corr == 2) call s_mass_lag(q_cons_vf, tau, .true.)

        ! Face correction with the operator's own conductance (area 1), so div(uf) matches the solved pressure exactly. Cells take
        ! the mean of their faces' net acceleration (body force less pressure gradient, zero on walls): a hydrostatic balance on
        ! the faces then leaves the cells at rest too, and for uniform density this is the centered pressure gradient
        do d = 1, num_dims
            call s_unit(d, oj, ok, ol)
            iv = eqn_idx%mom%beg + d - 1
            ga = tau*acc(d)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, gf, ia, s]')
            do l = face_bounds(3, d)%beg, face_bounds(3, d)%end
                do k = face_bounds(2, d)%beg, face_bounds(2, d)%end
                    do j = face_bounds(1, d)%beg, face_bounds(1, d)%end
                        gf = f_face_gf(rhoc(j, k, l), rhoc(j + oj, k + ok, l + ol), f_face_dist(d, j, k, l), real(pk(j, k, l), &
                                       & wp), real(pk(j + oj, k + ok, l + ol), wp), tau)
                        uf(j, k, l, d) = uf(j, k, l, d) - gf
                        pflx(j, k, l) = ga - gf
                        if (wbl) pflx(j, k, l) = pflx(j, k, l) + tau*f_face_capillary(d, j, k, l, oj, ok, ol, q_cons_vf)
                        ia = f_along(d, j, k, l)
                        s = 0
                        if (ia == -1) s = 1
                        if (ia == f_last(d)) s = 2
                        if (s > 0) then
                            if (f_wall_face(f_bc_code(d, j, k, l, s)) .or. f_bc_code(d, j, k, l, s) == BC_DIRICHLET) then
                                pflx(j, k, l) = 0._wp
                            end if
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(iv)%sf(j, k, l) = real(real(q_cons_vf(iv)%sf(j, k, l), wp) + rhoc(j, k, l)*0.5_wp*(pflx(j - oj, &
                                  & k - ok, l - ol) + pflx(j, k, l)), stp)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do
        call s_zero_wall_faces()

        if (proj_lag_corr > 0) call s_mass_lag(q_cons_vf, tau, .false.)

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

    !> The stage transported the partial densities, volume fractions and momentum with uf0, while the pressure solve compressed p
    !! with the corrected face velocity: left alone, density and pressure disagree by tau*rho*div(u_f - uf0) every stage, an entropy
    !! error that accumulates wherever that divergence persists (small scales, walls, low acoustic CFL), and the moved mass would
    !! not carry its momentum. Move them by the difference too, conservatively (alpha in its advective form), upwinded on the
    !! corrected face velocity. Momentum (mom) is moved before the correction, which forms that velocity here as the correction
    !! will, so every upwind value, ghosts included, is the pre-correction state on both sides of a rank seam; the densities and
    !! alpha after it, where uf is final and the correction has read the transported alpha. Each component is summed into rhs_p,
    !! free once the pressure system is built, then applied.
    subroutine s_mass_lag(q_cons_vf, tau, mom)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(in)                                   :: tau
        logical, intent(in)                                    :: mom
        real(wp)                                               :: c, dd, du, vn, tl
        integer                                                :: iq, q, j, k, l, nalpha, i1, i2, d, oj, ok, ol
        logical                                                :: pre

        tl = tau; pre = mom  ! a dummy may alias a host variable, which a device kernel must not reference
        ! A lone fluid's alpha is 1, which the difference leaves unchanged
        nalpha = merge(0, num_fluids, num_fluids == 1)
        i1 = merge(num_fluids + nalpha + 1, 1, mom)
        i2 = merge(num_fluids + nalpha + num_dims, num_fluids + nalpha, mom)
        do iq = i1, i2
            if (iq <= num_fluids) then
                q = iq
            else if (iq <= num_fluids + nalpha) then
                q = eqn_idx%adv%beg + iq - num_fluids - 1
            else
                q = eqn_idx%mom%beg + iq - num_fluids - nalpha - 1
            end if
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, c, dd, du, vn, d, oj, ok, ol]', firstprivate='[q, tl, pre]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        c = 0._wp; dd = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do d = 1, num_dims
                            call s_unit(d, oj, ok, ol)
                            ! High face: between this cell and the next
                            vn = uf(j, k, l, d)
                            ! The correction
                            if (pre) vn = vn - f_face_gf(rhoc(j, k, l), rhoc(j + oj, k + ok, l + ol), f_face_dist(d, j, k, l), &
                                & real(pk(j, k, l), wp), real(pk(j + oj, k + ok, l + ol), wp), tl)
                            if (f_along(d, j, k, l) == f_last(d)) then
                                if (f_wall_face(f_bc_code(d, j, k, l, 2))) vn = 0._wp
                            end if
                            du = vn - uf0(j, k, l, d)
                            c = c + du*real(merge(q_cons_vf(q)%sf(j, k, l), q_cons_vf(q)%sf(j + oj, k + ok, l + ol), &
                                            & vn >= 0._wp), wp)/f_width(d, j, k, l)
                            dd = dd + du/f_width(d, j, k, l)
                            ! Low face: between the previous cell and this one
                            vn = uf(j - oj, k - ok, l - ol, d)
                            if (pre) vn = vn - f_face_gf(rhoc(j - oj, k - ok, l - ol), rhoc(j, k, l), f_face_dist(d, j - oj, &
                                & k - ok, l - ol), real(pk(j - oj, k - ok, l - ol), wp), real(pk(j, k, l), wp), tl)
                            if (f_along(d, j, k, l) == 0) then
                                if (f_wall_face(f_bc_code(d, j, k, l, 1))) vn = 0._wp
                            end if
                            du = vn - uf0(j - oj, k - ok, l - ol, d)
                            c = c - du*real(merge(q_cons_vf(q)%sf(j - oj, k - ok, l - ol), q_cons_vf(q)%sf(j, k, l), &
                                            & vn >= 0._wp), wp)/f_width(d, j, k, l)
                            dd = dd - du/f_width(d, j, k, l)
                        end do
                        if (iq > num_fluids .and. iq <= num_fluids + nalpha) c = c - real(q_cons_vf(q)%sf(j, k, l), wp)*dd
                        rhs_p(j, k, l) = c
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', firstprivate='[q, tl]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(q)%sf(j, k, l) = real(real(q_cons_vf(q)%sf(j, k, l), wp) - tl*rhs_p(j, k, l), stp)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_mass_lag

    !> Face conductances of the pressure operator, with the operator's own face density (f_cond). A boundary face carries one only
    !! where it couples to something: a rank or periodic seam, or a pressure outlet's ghost at pout; elsewhere it is a wall to the
    !! solve
    subroutine s_projection_conductances()

        integer  :: j, k, l
        real(wp) :: area
        logical  :: sl1, sh1, sl2, sh2, sl3, sh3, po1, po2, po3

        sl1 = seam_lo(1); sh1 = seam_hi(1); sl2 = seam_lo(2); sh2 = seam_hi(2); sl3 = seam_lo(3); sh3 = seam_hi(3)
        po1 = has_pout(1); po2 = has_pout(2); po3 = has_pout(3)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, area]', firstprivate='[sl1, sh1, sl2, sh2, sl3, sh3, po1, po2, po3]')
        do l = 0, merge(p + 1, p, p > 0)
            do k = 0, merge(n + 1, n, n > 0)
                do j = 0, m + 1
                    if (k <= n .and. l <= p) then
                        cfx(j, k, l) = 0._wp
                        if ((j > 0 .and. j <= m) .or. (j == 0 .and. (sl1 .or. (po1 .and. fbx(k, l, &
                            & 1) == BC_GHOST_EXTRAP))) .or. (j == m + 1 .and. (sh1 .or. (po1 .and. fbx(k, l, &
                            & 2) == BC_GHOST_EXTRAP)))) then
                            area = 1._wp
                            if (num_dims > 1) area = dy(k)
                            if (num_dims > 2) area = area*dz(l)
                            cfx(j, k, l) = f_cond(rhoc(j, k, l), rhoc(j - 1, k, l), area, 0.5_wp*(dx(j - 1) + dx(j)))
                        end if
                    end if
                    if (num_dims > 1 .and. j <= m .and. l <= p) then
                        cfy(j, k, l) = 0._wp
                        if ((k > 0 .and. k <= n) .or. (k == 0 .and. (sl2 .or. (po2 .and. fby(j, l, &
                            & 1) == BC_GHOST_EXTRAP))) .or. (k == n + 1 .and. (sh2 .or. (po2 .and. fby(j, l, &
                            & 2) == BC_GHOST_EXTRAP)))) then
                            area = dx(j)
                            if (num_dims > 2) area = area*dz(l)
                            cfy(j, k, l) = f_cond(rhoc(j, k, l), rhoc(j, k - 1, l), area, 0.5_wp*(dy(k - 1) + dy(k)))
                        end if
                    end if
                    if (num_dims > 2 .and. j <= m .and. k <= n) then
                        cfz(j, k, l) = 0._wp
                        if ((l > 0 .and. l <= p) .or. (l == 0 .and. (sl3 .or. (po3 .and. fbz(j, k, &
                            & 1) == BC_GHOST_EXTRAP))) .or. (l == p + 1 .and. (sh3 .or. (po3 .and. fbz(j, k, &
                            & 2) == BC_GHOST_EXTRAP)))) then
                            cfz(j, k, l) = f_cond(rhoc(j, k, l), rhoc(j, k, l - 1), dx(j)*dy(k), 0.5_wp*(dz(l - 1) + dz(l)))
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_projection_conductances

    !> Solve the pressure system for pk, with filled ghosts, from the initial guess in xs
    impure subroutine s_pressure_solve(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp)                                                   :: res
        integer                                                    :: its, j, k, l

        call nvtxStartRange("PROJ-MG-BUILD")
        call s_projection_conductances()
        call s_mg_set_operator(dcoef, cfx, cfy, cfz)
        call nvtxEndRange
        if (any(has_pout)) call s_pout_rhs()
        call s_mg_solve(bvec, xs, its, res)
        proj_pcg_iters = proj_pcg_iters + its
        proj_pcg_res = max(proj_pcg_res, res)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    pk(j, k, l) = real(xs(j, k, l), stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_populate_F_igr_buffers(bc_type, pk_sf)
        if (any(has_pout)) call s_pout_ghosts()

    end subroutine s_pressure_solve

    !> A pressure outlet's fixed ghost enters the right-hand side as its face's conductance times pout
    subroutine s_pout_rhs()

        integer  :: d, s, g, ip, j, k, l, oj, ok, ol
        real(wp) :: po

        ! The boundary cell is 0 or m (n, p); its outlet face the cell's own low face or (g = 1) the next one

        do d = 1, num_dims
            if (.not. has_pout(d)) cycle
            po = pout(d)
            call s_unit(d, oj, ok, ol)
            do s = 1, 2
                g = s - 1
                ip = g
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', firstprivate='[po]')
                do l = plane_bounds(3, d, ip)%beg, plane_bounds(3, d, ip)%end
                    do k = plane_bounds(2, d, ip)%beg, plane_bounds(2, d, ip)%end
                        do j = plane_bounds(1, d, ip)%beg, plane_bounds(1, d, ip)%end
                            if (f_bc_code(d, j, k, l, s) == BC_GHOST_EXTRAP) bvec(j, k, l) = bvec(j, k, l) + f_cf(d, j + g*oj, &
                                & k + g*ok, l + g*ol)*po
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

    end subroutine s_pout_rhs

    !> A pressure outlet's ghost holds pout, so the face correction pushes flow out against it
    subroutine s_pout_ghosts()

        integer  :: d, s, ip, j, k, l
        real(wp) :: po

        do d = 1, num_dims
            if (.not. has_pout(d)) cycle
            po = pout(d)
            do s = 1, 2
                ! The ghost past the boundary face: -1, or m + 1 (n + 1, p + 1)
                ip = merge(-1, 2, s == 1)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', firstprivate='[po]')
                do l = plane_bounds(3, d, ip)%beg, plane_bounds(3, d, ip)%end
                    do k = plane_bounds(2, d, ip)%beg, plane_bounds(2, d, ip)%end
                        do j = plane_bounds(1, d, ip)%beg, plane_bounds(1, d, ip)%end
                            if (f_bc_code(d, j, k, l, s) == BC_GHOST_EXTRAP) pk(j, k, l) = real(po, stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

    end subroutine s_pout_ghosts

    impure subroutine s_finalize_projection_module()

        $:GPU_EXIT_DATA(detach='[pk_sf(1)%sf]')
        @:DEALLOCATE(uf, uf0, divu, rhs_p, p_stage, p_step0, pflx, rhoc, dcoef, bvec, xs, cfx, cfy, cfz, pk, kap, gnd, fbx, fby, &
                     & fbz)
        call s_finalize_multigrid_module()

    end subroutine s_finalize_projection_module

end module m_projection
