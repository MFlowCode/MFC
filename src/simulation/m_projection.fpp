!>
!! @file
!! @brief Contains module m_projection

!> @brief All-Mach pressure projection (Fuster & Popinet, JCP 374, 2018). Advection uses a persistent face velocity; the pressure
!! then solves p - rho*c^2*tau^2*div(rho_f^-1 grad p) = p_adv - rho*c^2*tau*div(u*_f) with div, grad and the Laplacian all taken on
!! faces, so they compose exactly and the projected face velocity satisfies the discrete pressure equation. The solve is PCG
!! preconditioned by geometric multigrid (m_multigrid).

#:include 'case.fpp'
#:include 'macros.fpp'

module m_projection

#ifdef MFC_MPI
    use mpi
#endif
    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_boundary_common
    use m_eos
    use m_body_forces
    use m_riemann_state
    use m_projection_helpers
    use m_multigrid
    use m_nvtx

    implicit none

    private; public :: s_initialize_projection_module, s_projection_rhs, s_projection_face_props, s_projection_heat, &
        & s_projection_apply, s_finalize_projection_module

    real(wp), allocatable, dimension(:,:,:,:) :: u_face             !< face velocity; index j is the face between cells j and j+1
    real(wp), allocatable, dimension(:,:,:,:) :: u0_face            !< the face velocity that transported the current stage
    real(wp), allocatable, dimension(:,:,:)   :: divu_face, p_rate  !< div of u_face, and the pressure transport rate
    !> upwind pressure flux on the faces of one direction, then their acceleration in the correction
    real(wp), allocatable, dimension(:,:,:) :: face_flux
    real(wp), allocatable, dimension(:,:,:) :: p_stage, p_step0
    real(wp), allocatable, dimension(:,:,:) :: rhostar            !< star density with one ghost layer
    real(wp), allocatable, dimension(:,:,:) :: p_diag, p_sys_rhs  !< SPD system: D_c and right-hand side
    real(wp), allocatable, dimension(:,:,:) :: p_guess            !< pressure-solve initial guess, then solution
    !> Face conductances A_f/(rho_f d_f) of the pressure operator, x, y and z faces; indexed as u_face
    real(wp), allocatable, dimension(:,:,:)          :: cond_x, cond_y, cond_z
    real(stp), allocatable, dimension(:,:,:), target :: p_solved  !< solved pressure, with ghosts for the halo
    type(scalar_field), dimension(1)                 :: p_solved_sf
    $:GPU_DECLARE(create='[p_solved_sf]')
    $:GPU_DECLARE(create='[u_face, u0_face, divu_face, p_rate, face_flux, p_stage, p_step0, rhostar, p_diag, p_sys_rhs, p_guess]')
    $:GPU_DECLARE(create='[cond_x, cond_y, cond_z, p_solved]')
    !> Well-balanced surface tension, set over one ghost layer: curvature (1) and |grad c| (2), both zero outside the interface band
    real(wp), allocatable, dimension(:,:,:,:) :: curv
    !> grad(alpha_1) (1:num_dims) and its magnitude (0), set over two ghost layers
    real(wp), allocatable, dimension(:,:,:,:) :: grad_alpha
    $:GPU_DECLARE(create='[curv, grad_alpha]')

    !> Boundary code of each domain face, low (1) and high (2) side, copied from bc_type: walls and boundary patches per face
    integer, allocatable, dimension(:,:,:) :: bc_face_x, bc_face_y, bc_face_z
    $:GPU_DECLARE(create='[bc_face_x, bc_face_y, bc_face_z]')
    !> An extrapolation (-3) face is a pressure outlet held at bc_[x,y,z]%pres_out where that is set: its conductance couples the
    !! boundary cell to a ghost at pres_out, so flow leaves where the solve would otherwise see a wall
    logical, dimension(3)  :: has_pres_out
    real(wp), dimension(3) :: pres_out
    logical, dimension(3)  :: seam_lo, seam_hi  !< that side couples to another rank or, periodically, to this one
    $:GPU_DECLARE(create='[has_pres_out, pres_out, seam_lo, seam_hi]')
    logical :: faces_ready       !< u_face has been seeded from the cell velocities
    logical :: well_balanced_st  !< well-balanced surface tension
    $:GPU_DECLARE(create='[well_balanced_st]')
    !> Loop bounds by direction d: face_bounds(:, d) spans the d faces (-1 to m, n, p along d); plane_bounds(:, d, ip) the plane
    !! normal to d at index -1 (ip = -1), 0 (ip = 0), the last cell (ip = 1) or the ghost past it (ip = 2)
    type(int_bounds_info), dimension(3, 3)      :: face_bounds
    type(int_bounds_info), dimension(3, 3,-1:2) :: plane_bounds
    !> The cells and g ghost layers (g = 1, 2) in each active direction: ghost_bounds(:, g)
    type(int_bounds_info), dimension(3, 2) :: ghost_bounds
    $:GPU_DECLARE(create='[face_bounds, plane_bounds, ghost_bounds]')

contains

    impure subroutine s_initialize_projection_module()

        integer :: n_ghost

#ifdef MFC_MIXED_PRECISION
        call s_mpi_abort('proj_method needs stp = wp; mixed precision is not supported')
#endif

        @:ALLOCATE(u_face(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:num_dims))
        @:ALLOCATE(u0_face(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:num_dims))
        @:ALLOCATE(divu_face(0:m, 0:n, 0:p), p_rate(0:m, 0:n, 0:p), p_stage(0:m, 0:n, 0:p), p_step0(0:m, 0:n, 0:p))
        @:ALLOCATE(face_flux(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end), &
                   & rhostar(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ALLOCATE(p_diag(0:m, 0:n, 0:p), p_sys_rhs(0:m, 0:n, 0:p))
        @:ALLOCATE(p_guess(0:m, 0:n, 0:p), cond_x(-1:m, 0:n, 0:p), cond_y(0:m, -1:n, 0:p), cond_z(0:m, 0:n, -1:p))
        @:ALLOCATE(p_solved(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))

        ! Zeroed so that copying u_face over the whole buffer (s_projection_apply) never moves unset values
        u_face = 0._wp
        $:GPU_UPDATE(device='[u_face]')

        p_solved_sf(1)%sf => p_solved  ! Alias to reuse existing buffer population routine
        $:GPU_ENTER_DATA(copyin='[p_solved_sf(1)%sf]')
        $:GPU_ENTER_DATA(attach='[p_solved_sf(1)%sf]')

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            seam_lo(${NORM_DIR}$) = bc_${XYZ}$%beg >= 0 .or. bc_${XYZ}$%beg == BC_PERIODIC
            seam_hi(${NORM_DIR}$) = bc_${XYZ}$%end >= 0 .or. bc_${XYZ}$%end == BC_PERIODIC
            has_pres_out(${NORM_DIR}$) = .not. f_is_default(bc_${XYZ}$%pres_out) .and. num_dims >= ${NORM_DIR}$
            pres_out(${NORM_DIR}$) = merge(bc_${XYZ}$%pres_out, 0._wp, has_pres_out(${NORM_DIR}$))
        #:endfor
        $:GPU_UPDATE(device='[has_pres_out, pres_out, seam_lo, seam_hi]')
        call s_initialize_multigrid_module()

        @:ALLOCATE(bc_face_x(0:n, 0:p, 1:2), bc_face_y(0:merge(m, 0, n > 0), 0:p, 1:2), bc_face_z(0:merge(m, 0, p > 0), &
                   & 0:merge(n, 0, p > 0), 1:2))
        faces_ready = .false.
        do n_ghost = 1, 2
            ghost_bounds(:,n_ghost)%beg = 0
            ghost_bounds(1, n_ghost)%end = m
            ghost_bounds(2, n_ghost)%end = n
            ghost_bounds(3, n_ghost)%end = p
        end do
        #:for NORM_DIR, LAST in [(1, 'm'), (2, 'n'), (3, 'p')]
            if (num_dims >= ${NORM_DIR}$) then
                call s_set_bounds(${NORM_DIR}$, -1, ${LAST}$, face_bounds(:,${NORM_DIR}$))
                call s_set_bounds(${NORM_DIR}$, -1, -1, plane_bounds(:,${NORM_DIR}$, -1))
                call s_set_bounds(${NORM_DIR}$, 0, 0, plane_bounds(:,${NORM_DIR}$, 0))
                call s_set_bounds(${NORM_DIR}$, ${LAST}$, ${LAST}$, plane_bounds(:,${NORM_DIR}$, 1))
                call s_set_bounds(${NORM_DIR}$, ${LAST}$ + 1, ${LAST}$ + 1, plane_bounds(:,${NORM_DIR}$, 2))
                do n_ghost = 1, 2
                    ghost_bounds(${NORM_DIR}$, n_ghost)%beg = -n_ghost
                    ghost_bounds(${NORM_DIR}$, n_ghost)%end = ${LAST}$ + n_ghost
                end do
            end if
        #:endfor
        $:GPU_UPDATE(device='[face_bounds, plane_bounds, ghost_bounds]')

        well_balanced_st = surface_tension .and. surface_tension_model == surface_tension_model_well_balanced
        $:GPU_UPDATE(device='[well_balanced_st]')
        if (well_balanced_st) then
            @:ALLOCATE(curv(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:2))
            @:ALLOCATE(grad_alpha(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 0:3))
        else
            @:ALLOCATE(curv(0:0, 0:0, 0:0, 1:2), grad_alpha(0:0, 0:0, 0:0, 0:3))
        end if

    end subroutine s_initialize_projection_module

    !> Copy the boundary code of each domain face from bc_type, on the host where pre-process left it
    impure subroutine s_projection_face_codes(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer                                                    :: side

        do side = 1, 2
            bc_face_x(:,:,side) = int(bc_type(1, side)%sf(0,:,:))
            if (num_dims > 1) bc_face_y(:,:,side) = int(bc_type(2, side)%sf(0:m,0,:))
            if (num_dims > 2) bc_face_z(:,:,side) = int(bc_type(3, side)%sf(0:m,0:n,0))
        end do
        $:GPU_UPDATE(device='[bc_face_x, bc_face_y, bc_face_z]')

    end subroutine s_projection_face_codes

    !> Face velocity from the average of the cell velocities either side (primitive velocities, with ghosts)
    subroutine s_projection_init_faces(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: mom_idx, j, k, l

        #:for NORM_DIR, STENCIL_VAR, COORDS in &
                [(1, 'j', '{STENCIL_IDX}, k, l'), &
                 (2, 'k', 'j, {STENCIL_IDX}, l'), &
                 (3, 'l', 'j, k, {STENCIL_IDX}')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (num_dims >= ${NORM_DIR}$) then
                mom_idx = eqn_idx%mom%beg + ${NORM_DIR - 1}$
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = face_bounds(3, ${NORM_DIR}$)%beg, face_bounds(3, ${NORM_DIR}$)%end
                    do k = face_bounds(2, ${NORM_DIR}$)%beg, face_bounds(2, ${NORM_DIR}$)%end
                        do j = face_bounds(1, ${NORM_DIR}$)%beg, face_bounds(1, ${NORM_DIR}$)%end
                            u_face(j, k, l, ${NORM_DIR}$) = 0.5_wp*(real(q_prim_vf(mom_idx)%sf(j, k, l), &
                                   & wp) + real(q_prim_vf(mom_idx)%sf(${SF(' + 1')}$), wp))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor
        call s_zero_wall_faces()

    end subroutine s_projection_init_faces

    !> A Dirichlet (-17) face carries the normal velocity of its ghost cell, which s_dirichlet sets from the boundary buffer and its
    !! inflow ramp
    subroutine s_dirichlet_faces(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: mom_idx, j, k, l

        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS, BC_COORDS in &
                [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'k, l'), &
                 (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'j, l'), &
                 (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'j, k')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (num_dims >= ${NORM_DIR}$) then
                mom_idx = eqn_idx%mom%beg + ${NORM_DIR - 1}$
                #:for S, IP, GHOST in [(1, -1, ''), (2, 1, ' + 1')]
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = plane_bounds(3, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(3, ${NORM_DIR}$, ${IP}$)%end
                        do k = plane_bounds(2, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(2, ${NORM_DIR}$, ${IP}$)%end
                            do j = plane_bounds(1, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(1, ${NORM_DIR}$, ${IP}$)%end
                                if (bc_face_${XYZ}$(${BC_COORDS}$, ${S}$) == BC_DIRICHLET) then
                                    u_face(j, k, l, ${NORM_DIR}$) = real(q_prim_vf(mom_idx)%sf(${SF(GHOST)}$), wp)
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:endfor
            end if
        #:endfor

    end subroutine s_dirichlet_faces

    !> Normal velocity vanishes on solid walls
    subroutine s_zero_wall_faces()

        integer :: j, k, l

        #:for NORM_DIR, XYZ, BC_COORDS in [(1, 'x', 'k, l'), (2, 'y', 'j, l'), (3, 'z', 'j, k')]
            if (num_dims >= ${NORM_DIR}$) then
                #:for S, IP in [(1, -1), (2, 1)]
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = plane_bounds(3, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(3, ${NORM_DIR}$, ${IP}$)%end
                        do k = plane_bounds(2, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(2, ${NORM_DIR}$, ${IP}$)%end
                            do j = plane_bounds(1, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(1, ${NORM_DIR}$, ${IP}$)%end
                                if (f_wall_face(bc_face_${XYZ}$(${BC_COORDS}$, ${S}$))) then
                                    u_face(j, k, l, ${NORM_DIR}$) = 0._wp
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:endfor
            end if
        #:endfor

    end subroutine s_zero_wall_faces

    !> Advective right-hand side of one direction sweep. Every quantity is carried by the projected face velocity and upwinded on
    !! its sign; the momentum flux is the summed partial-density flux times the upwind velocity, so mass and momentum move with one
    !! operator. Energy is left at zero here and rebuilt from the equation of state after the pressure solve.
    subroutine s_projection_rhs(id, qfl_rs, qfr_rs, q_prim_vf, flux_vf, rhs_vf, bc_type)

        integer, intent(in) :: id
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: qfl_rs, qfr_rs
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, rhs_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp) :: u_f, alpha_up, alpha_rho_up, mass_flux
        logical :: upwind_l
        integer :: i, j, k, l

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
                        divu_face(j, k, l) = f_div_u_face(j, k, l)
                        p_stage(j, k, l) = real(q_prim_vf(eqn_idx%E)%sf(j, k, l), wp)
                        ! Transport sources q*div(u), which make alpha and p advect rather than compress
                        p_rate(j, k, l) = p_stage(j, k, l)*divu_face(j, k, l)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, sys_size
                            rhs_vf(i)%sf(j, k, l) = 0._stp
                        end do
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            rhs_vf(i)%sf(j, k, l) = real(real(q_prim_vf(i)%sf(j, k, l), wp)*divu_face(j, k, l), stp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_dirichlet_faces(q_prim_vf)
        end if

        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS in &
                [(1, 'x', 'j', '{STENCIL_IDX}, k, l'), &
                 (2, 'y', 'k', 'j, {STENCIL_IDX}, l'), &
                 (3, 'z', 'l', 'j, k, {STENCIL_IDX}')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (id == ${NORM_DIR}$) then
                ! Face fluxes. Left state of a face is the right edge of the cell before it, right state the left edge of the next
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, u_f, alpha_up, alpha_rho_up, mass_flux, upwind_l]')
                do l = face_bounds(3, ${NORM_DIR}$)%beg, face_bounds(3, ${NORM_DIR}$)%end
                    do k = face_bounds(2, ${NORM_DIR}$)%beg, face_bounds(2, ${NORM_DIR}$)%end
                        do j = face_bounds(1, ${NORM_DIR}$)%beg, face_bounds(1, ${NORM_DIR}$)%end
                            u_f = u_face(j, k, l, ${NORM_DIR}$)
                            upwind_l = u_f >= 0._wp
                            mass_flux = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                ! Reconstructed partial densities upwinded directly. Taking the phase density alpha_rho/alpha
                                ! from the upwind cell instead diverges where a phase is vanishing: both are round-off there,
                                ! and their ratio (seen at 6e10) times the face alpha flux injects mass
                                if (upwind_l) then
                                    alpha_up = qfl_rs(j, k, l, eqn_idx%adv%beg + i - 1)
                                    alpha_rho_up = qfl_rs(j, k, l, i)
                                else
                                    alpha_up = qfr_rs(${SF(' + 1')}$, eqn_idx%adv%beg + i - 1)
                                    alpha_rho_up = qfr_rs(${SF(' + 1')}$, i)
                                end if
                                flux_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l) = real(alpha_up*u_f, stp)
                                flux_vf(i)%sf(j, k, l) = real(alpha_rho_up*u_f, stp)
                                mass_flux = mass_flux + alpha_rho_up*u_f
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                if (upwind_l) then
                                    flux_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = real(mass_flux*qfl_rs(j, k, l, &
                                            & eqn_idx%mom%beg + i - 1), stp)
                                else
                                    flux_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = real(mass_flux*qfr_rs(${SF(' + 1')}$, &
                                            & eqn_idx%mom%beg + i - 1), stp)
                                end if
                            end do
                            if (upwind_l) then
                                face_flux(j, k, l) = u_f*qfl_rs(j, k, l, eqn_idx%E)
                            else
                                face_flux(j, k, l) = u_f*qfr_rs(${SF(' + 1')}$, eqn_idx%E)
                            end if
                            ! The color function is not reconstructed; upwind its cell value
                            if (surface_tension) then
                                if (upwind_l) then
                                    flux_vf(eqn_idx%c)%sf(j, k, l) = real(u_f*real(q_prim_vf(eqn_idx%c)%sf(j, k, l), wp), stp)
                                else
                                    flux_vf(eqn_idx%c)%sf(j, k, l) = real(u_f*real(q_prim_vf(eqn_idx%c)%sf(${SF(' + 1')}$), wp), &
                                            & stp)
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
                                           & wp) - real(flux_vf(i)%sf(j, k, l), wp))/d${XYZ}$(${SV}$), stp)
                                end if
                            end do
                            p_rate(j, k, l) = p_rate(j, k, l) + (face_flux(${SF(' - 1')}$) - face_flux(j, k, l))/d${XYZ}$(${SV}$)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_rhs

    !> Face data that the viscous and capillary source fluxes otherwise take from a Riemann solve: interface Reynolds numbers, the
    !! mean face velocity, and for surface tension the face velocity whose divergence makes the color function advect
    subroutine s_projection_face_props(id, qfl_rs, qfr_rs, flux_src_vf)

        integer, intent(in)                                                                 :: id
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: qfl_rs, qfr_rs
        type(scalar_field), dimension(sys_size), intent(inout)                              :: flux_src_vf

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_L, alpha_R
        #:else
            real(wp), dimension(num_fluids) :: alpha_L, alpha_R
        #:endif
        real(wp), dimension(2) :: Re_L, Re_R
        integer                :: i, j, k, l, Re_size_loc1, Re_size_loc2

        Re_size_loc1 = Re_size(1)
        Re_size_loc2 = Re_size(2)

        #:for NORM_DIR, STENCIL_VAR, COORDS in &
                [(1, 'j', '{STENCIL_IDX}, k, l'), &
                 (2, 'k', 'j, {STENCIL_IDX}, l'), &
                 (3, 'l', 'j, k, {STENCIL_IDX}')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (id == ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, alpha_L, alpha_R, Re_L, Re_R]', &
                                    & firstprivate='[Re_size_loc1, Re_size_loc2]')
                do l = face_bounds(3, ${NORM_DIR}$)%beg, face_bounds(3, ${NORM_DIR}$)%end
                    do k = face_bounds(2, ${NORM_DIR}$)%beg, face_bounds(2, ${NORM_DIR}$)%end
                        do j = face_bounds(1, ${NORM_DIR}$)%beg, face_bounds(1, ${NORM_DIR}$)%end
                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qfl_rs(j, k, l, eqn_idx%adv%beg + i - 1)
                                    alpha_R(i) = qfr_rs(${SF(' + 1')}$, eqn_idx%adv%beg + i - 1)
                                end do
                                call s_compute_interface_reynolds(alpha_L, Re_L, Re_size_loc1, Re_size_loc2)
                                call s_compute_interface_reynolds(alpha_R, Re_R, Re_size_loc1, Re_size_loc2)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_avg_rsx_vf(j, k, l, i) = 2._wp/(1._wp/Re_L(i) + 1._wp/Re_R(i))
                                end do
                            end if
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
                                vel_src_rsx_vf(j, k, l, i) = 0.5_wp*(qfl_rs(j, k, l, &
                                               & eqn_idx%mom%beg + i - 1) + qfr_rs(${SF(' + 1')}$, eqn_idx%mom%beg + i - 1))
                            end do
                            if (surface_tension) flux_src_vf(eqn_idx%adv%beg)%sf(j, k, l) = real(u_face(j, k, l, ${NORM_DIR}$), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_face_props

    !> Conduction heating of one direction sweep from the energy source flux fe (-k dT/dn on faces). At fixed alpha and phase
    !! densities d(rho e) = Gamma dp, so the pressure transport rate gains -div(q)/Gamma, and the energy rebuilt from p carries it
    subroutine s_projection_heat(id, q_prim_vf, fe)

        integer, intent(in)                                 :: id
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), intent(in)                      :: fe
        real(wp)                                            :: rho, gamma, pi_inf, qv

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho, alpha
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho, alpha
        #:endif
        integer :: i, j, k, l

        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS in &
                [(1, 'x', 'j', '{STENCIL_IDX}, k, l'), &
                 (2, 'y', 'k', 'j, {STENCIL_IDX}, l'), &
                 (3, 'z', 'l', 'j, k, {STENCIL_IDX}')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (id == ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho, gamma, pi_inf, qv, alpha_rho, alpha]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho(i) = real(q_prim_vf(i)%sf(j, k, l), wp)
                                alpha(i) = real(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                            end do
                            call s_compute_mixture_coefficients(alpha_rho, alpha, rho, gamma, pi_inf, qv)
                            p_rate(j, k, l) = p_rate(j, k, l) + (real(fe%sf(${SF(' - 1')}$), wp) - real(fe%sf(j, k, l), &
                                   & wp))/(d${XYZ}$(${SV}$)*gamma)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_heat

    !> Divergence of the face velocity in cell (j, k, l): the one operator the transport sources and the pressure equation share
    function f_div_u_face(j, k, l) result(div_u)

        $:GPU_ROUTINE(function_name='f_div_u_face', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: j, k, l
        real(wp)            :: div_u

        div_u = (u_face(j, k, l, 1) - u_face(j - 1, k, l, 1))/dx(j)
        if (num_dims > 1) div_u = div_u + (u_face(j, k, l, 2) - u_face(j, k - 1, l, 2))/dy(k)
        if (num_dims > 2) div_u = div_u + (u_face(j, k, l, 3) - u_face(j, k, l - 1, 3))/dz(l)

    end function f_div_u_face

    !> Capillary acceleration of the face between cells (j, k, l) and (jn, kn, ln), their centres dist apart
    function f_face_capillary(j, k, l, j2, k2, l2, dist, q_cons_vf) result(acc)

        $:GPU_ROUTINE(function_name='f_face_capillary', parallelism='[seq]', cray_inline=True)

        integer, intent(in)                                 :: j, k, l, j2, k2, l2
        real(wp), intent(in)                                :: dist
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp)                                            :: acc

        acc = f_capillary_accel(curv(j, k, l, 1), curv(j2, k2, l2, 1), curv(j, k, l, 2), curv(j2, k2, l2, 2), &
                                & real(q_cons_vf(eqn_idx%adv%beg)%sf(j, k, l), wp), real(q_cons_vf(eqn_idx%adv%beg)%sf(j2, k2, &
                                & l2), wp), rhostar(j, k, l), rhostar(j2, k2, l2), dist)

    end function f_face_capillary

    !> Interface normal component d at cell (j, k, l), from grad_alpha
    function f_interface_normal(d, j, k, l) result(n_d)

        $:GPU_ROUTINE(function_name='f_interface_normal', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: d, j, k, l
        real(wp)            :: n_d

        n_d = f_normal(grad_alpha(j, k, l, d), grad_alpha(j, k, l, 0))

    end function f_interface_normal

    !> Curvature kappa = -div(n) in the interface band, with n = grad(alpha_1)/|grad(alpha_1)|.
    subroutine s_compute_curvature(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp)                                            :: curv_c
        integer                                             :: j, k, l, adv_idx

        adv_idx = eqn_idx%adv%beg
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = ghost_bounds(3, 2)%beg, ghost_bounds(3, 2)%end
            do k = ghost_bounds(2, 2)%beg, ghost_bounds(2, 2)%end
                do j = ghost_bounds(1, 2)%beg, ghost_bounds(1, 2)%end
                    grad_alpha(j, k, l, 1) = real(q_cons_vf(adv_idx)%sf(j + 1, k, l) - q_cons_vf(adv_idx)%sf(j - 1, k, l), &
                               & wp)/(x_cc(j + 1) - x_cc(j - 1))
                    grad_alpha(j, k, l, 2) = 0._wp
                    grad_alpha(j, k, l, 3) = 0._wp
                    if (num_dims > 1) grad_alpha(j, k, l, 2) = real(q_cons_vf(adv_idx)%sf(j, k + 1, l) - q_cons_vf(adv_idx)%sf(j, &
                        & k - 1, l), wp)/(y_cc(k + 1) - y_cc(k - 1))
                    if (num_dims > 2) grad_alpha(j, k, l, 3) = real(q_cons_vf(adv_idx)%sf(j, k, l + 1) - q_cons_vf(adv_idx)%sf(j, &
                        & k, l - 1), wp)/(z_cc(l + 1) - z_cc(l - 1))
                    grad_alpha(j, k, l, 0) = sqrt(grad_alpha(j, k, l, 1)**2 + grad_alpha(j, k, l, 2)**2 + grad_alpha(j, k, l, 3)**2)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, curv_c]')
        do l = ghost_bounds(3, 1)%beg, ghost_bounds(3, 1)%end
            do k = ghost_bounds(2, 1)%beg, ghost_bounds(2, 1)%end
                do j = ghost_bounds(1, 1)%beg, ghost_bounds(1, 1)%end
                    curv(j, k, l, 1) = 0._wp
                    curv(j, k, l, 2) = 0._wp
                    if (grad_alpha(j, k, l, 0) > capillary_cutoff) then
                        curv_c = -(f_interface_normal(1, j + 1, k, l) - f_interface_normal(1, j - 1, k, &
                                   & l))/(x_cc(j + 1) - x_cc(j - 1))
                        if (num_dims > 1) curv_c = curv_c - (f_interface_normal(2, j, k + 1, l) - f_interface_normal(2, j, k - 1, &
                            & l))/(y_cc(k + 1) - y_cc(k - 1))
                        if (num_dims > 2) curv_c = curv_c - (f_interface_normal(3, j, k, l + 1) - f_interface_normal(3, j, k, &
                            & l - 1))/(z_cc(l + 1) - z_cc(l - 1))
                        curv(j, k, l, 1) = curv_c
                        curv(j, k, l, 2) = grad_alpha(j, k, l, 0)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_curvature

    !> Pressure solve and correction on the blended (star) state of one RK stage
    impure subroutine s_projection_apply(q_cons_vf, bc_type, pb_in, mv_in, q_T_sf, rkc, stage)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(scalar_field), intent(inout) :: q_T_sf
        real(wp), dimension(4), intent(in) :: rkc
        integer, intent(in) :: stage
        real(wp) :: tau, rho, gamma, pi_inf, qv, rho_c2, div_u
        real(wp) :: vol, dyn_pres, dvel_bf, dvel_p
        real(wp), dimension(3) :: accel
        logical :: dirichlet

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho, alpha
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho, alpha
        #:endif
        integer :: i, j, k, l, mom_idx, bc_code

        tau = rkc(3)*dt/rkc(4)

        call s_populate_variables_buffers(bc_type, q_cons_vf, pb_in, mv_in, q_T_sf)

        accel = 0._wp
        if (bodyForces) then
            call s_compute_acceleration(mytime)
            #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (bf_${XYZ}$) accel(${D}$) = accel_bf(${D}$)
            #:endfor
        end if

        ! Keep the face velocity this stage was transported with, for s_mass_lag below
        if (proj_lag_corr > 0) then
            $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
            do i = 1, num_dims
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            u0_face(j, k, l, i) = u_face(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! Star density with its ghosts, and the face predictor from the star cell velocities
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    rhostar(j, k, l) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        rhostar(j, k, l) = rhostar(j, k, l) + real(q_cons_vf(i)%sf(j, k, l), wp)
                    end do
                    rhostar(j, k, l) = max(rhostar(j, k, l), sgm_eps)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (well_balanced_st) call s_compute_curvature(q_cons_vf)

        ! Body forces and well-balanced surface tension enter on faces, where they meet the pressure gradient that balances them.
        ! A Dirichlet face keeps its prescribed velocity
        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS, BC_COORDS, LAST in &
                [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'k, l', 'm'), &
                 (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'j, l', 'n'), &
                 (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'j, k', 'p')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (num_dims >= ${NORM_DIR}$) then
                mom_idx = eqn_idx%mom%beg + ${NORM_DIR - 1}$
                dvel_bf = tau*accel(${NORM_DIR}$)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, dirichlet]')
                do l = face_bounds(3, ${NORM_DIR}$)%beg, face_bounds(3, ${NORM_DIR}$)%end
                    do k = face_bounds(2, ${NORM_DIR}$)%beg, face_bounds(2, ${NORM_DIR}$)%end
                        do j = face_bounds(1, ${NORM_DIR}$)%beg, face_bounds(1, ${NORM_DIR}$)%end
                            dirichlet = .false.
                            if (${SV}$ == -1) dirichlet = bc_face_${XYZ}$(${BC_COORDS}$, 1) == BC_DIRICHLET
                            if (${SV}$ == ${LAST}$) dirichlet = bc_face_${XYZ}$(${BC_COORDS}$, 2) == BC_DIRICHLET
                            if (.not. dirichlet) then
                                u_face(j, k, l, ${NORM_DIR}$) = 0.5_wp*(real(q_cons_vf(mom_idx)%sf(j, k, l), wp)/rhostar(j, k, &
                                       & l) + real(q_cons_vf(mom_idx)%sf(${SF(' + 1')}$), wp)/rhostar(${SF(' + 1')}$)) + dvel_bf
                                if (well_balanced_st) u_face(j, k, l, ${NORM_DIR}$) = u_face(j, k, l, &
                                    & ${NORM_DIR}$) + tau*f_face_capillary(j, k, l, ${SF(' + 1')}$, &
                                    & 0.5_wp*(d${XYZ}$(${SV}$) + d${XYZ}$(${SV}$ + 1)), q_cons_vf)
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor
        call s_zero_wall_faces()

        ! SPD system D_c p + sum_f K_f (p - p_nb) = b, the Helmholtz row scaled by V_c/(rho c^2 tau^2)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho, gamma, pi_inf, qv, rho_c2, div_u, vol, alpha_rho, alpha]', &
                            & copyin='[rkc]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (stage == 1) p_step0(j, k, l) = p_stage(j, k, l)
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_rho(i) = real(q_cons_vf(i)%sf(j, k, l), wp)
                        alpha(i) = real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                    end do
                    call s_compute_mixture_coefficients(alpha_rho, alpha, rho, gamma, pi_inf, qv)
                    rho_c2 = max(f_bulk_modulus(p_stage(j, k, l), gamma, pi_inf), sgm_eps)
                    div_u = f_div_u_face(j, k, l)
                    vol = dx(j)
                    if (num_dims > 1) vol = vol*dy(k)
                    if (num_dims > 2) vol = vol*dz(l)
                    p_diag(j, k, l) = vol/(rho_c2*tau*tau)
                    p_sys_rhs(j, k, l) = p_diag(j, k, l)*((rkc(1)*p_stage(j, k, l) + rkc(2)*p_step0(j, k, &
                              & l) + rkc(3)*dt*p_rate(j, k, l))/rkc(4) - rho_c2*tau*div_u)
                    p_guess(j, k, l) = p_stage(j, k, l)
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

        ! Face correction with the operator's own conductance (area 1), so div(u_face) matches the solved pressure exactly. Cells
        ! take the mean of their faces' net acceleration (body force less pressure gradient, zero on walls): a hydrostatic balance
        ! on the faces then leaves the cells at rest too, and for uniform density this is the centered pressure gradient
        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS, BC_COORDS, LAST in &
                [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'k, l', 'm'), &
                 (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'j, l', 'n'), &
                 (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'j, k', 'p')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (num_dims >= ${NORM_DIR}$) then
                mom_idx = eqn_idx%mom%beg + ${NORM_DIR - 1}$
                dvel_bf = tau*accel(${NORM_DIR}$)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, dvel_p, bc_code]')
                do l = face_bounds(3, ${NORM_DIR}$)%beg, face_bounds(3, ${NORM_DIR}$)%end
                    do k = face_bounds(2, ${NORM_DIR}$)%beg, face_bounds(2, ${NORM_DIR}$)%end
                        do j = face_bounds(1, ${NORM_DIR}$)%beg, face_bounds(1, ${NORM_DIR}$)%end
                            dvel_p = f_face_dvel_p(rhostar(j, k, l), rhostar(${SF(' + 1')}$), &
                                                   & 0.5_wp*(d${XYZ}$(${SV}$) + d${XYZ}$(${SV}$ + 1)), real(p_solved(j, k, l), &
                                                   & wp), real(p_solved(${SF(' + 1')}$), wp), tau)
                            u_face(j, k, l, ${NORM_DIR}$) = u_face(j, k, l, ${NORM_DIR}$) - dvel_p
                            face_flux(j, k, l) = dvel_bf - dvel_p
                            if (well_balanced_st) face_flux(j, k, l) = face_flux(j, k, l) + tau*f_face_capillary(j, k, l, &
                                & ${SF(' + 1')}$, 0.5_wp*(d${XYZ}$(${SV}$) + d${XYZ}$(${SV}$ + 1)), q_cons_vf)
                            #:for S, EDGE in [(1, '-1'), (2, LAST)]
                                if (${SV}$ == ${EDGE}$) then
                                    bc_code = bc_face_${XYZ}$(${BC_COORDS}$, ${S}$)
                                    if (f_wall_face(bc_code) .or. bc_code == BC_DIRICHLET) face_flux(j, k, l) = 0._wp
                                end if
                            #:endfor
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_vf(mom_idx)%sf(j, k, l) = real(real(q_cons_vf(mom_idx)%sf(j, k, l), wp) + rhostar(j, k, &
                                      & l)*0.5_wp*(face_flux(${SF(' - 1')}$) + face_flux(j, k, l)), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor
        call s_zero_wall_faces()

        if (proj_lag_corr > 0) call s_mass_lag(q_cons_vf, tau, .false.)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho, gamma, pi_inf, qv, dyn_pres, alpha_rho, alpha]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_rho(i) = real(q_cons_vf(i)%sf(j, k, l), wp)
                        alpha(i) = real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                    end do
                    call s_compute_mixture_coefficients(alpha_rho, alpha, rho, gamma, pi_inf, qv)
                    dyn_pres = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        dyn_pres = dyn_pres + 0.5_wp*real(q_cons_vf(i)%sf(j, k, l), wp)*(real(q_cons_vf(i)%sf(j, k, l), wp)/rho)
                    end do
                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = real(gamma*real(p_solved(j, k, l), wp) + pi_inf + qv + dyn_pres, stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_projection_apply

    !> The stage transported the partial densities, volume fractions and momentum with u0_face, while the pressure solve compressed
    !! p with the corrected face velocity: left alone, density and pressure disagree by tau*rho*div(u_f - u0_face) every stage, an
    !! entropy error that accumulates wherever that divergence persists (small scales, walls, low acoustic CFL), and the moved mass
    !! would not carry its momentum. Move them by the difference too, conservatively (alpha in its advective form), upwinded on the
    !! corrected face velocity. Momentum (mom) is moved before the correction, which forms that velocity here as the correction
    !! will, so every upwind value, ghosts included, is the pre-correction state on both sides of a rank seam; the densities and
    !! alpha after it, where u_face is final and the correction has read the transported alpha. Each component is summed into
    !! p_rate, free once the pressure system is built, then applied.
    subroutine s_mass_lag(q_cons_vf, tau, mom)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(in)                                   :: tau
        logical, intent(in)                                    :: mom
        real(wp)                                               :: lag_div_q, lag_div_u, du_face, u_f, tau_loc
        integer                                                :: lag_var, q_idx, j, k, l, n_alpha, var_beg, var_end
        logical                                                :: mom_loc

        ! A dummy may alias a host variable, which a device kernel must not reference
        tau_loc = tau
        mom_loc = mom
        ! A lone fluid's alpha is 1, which the difference leaves unchanged
        n_alpha = merge(0, num_fluids, num_fluids == 1)
        var_beg = merge(num_fluids + n_alpha + 1, 1, mom)
        var_end = merge(num_fluids + n_alpha + num_dims, num_fluids + n_alpha, mom)
        do lag_var = var_beg, var_end
            if (lag_var <= num_fluids) then
                q_idx = lag_var
            else if (lag_var <= num_fluids + n_alpha) then
                q_idx = eqn_idx%adv%beg + lag_var - num_fluids - 1
            else
                q_idx = eqn_idx%mom%beg + lag_var - num_fluids - n_alpha - 1
            end if
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, lag_div_q, lag_div_u, du_face, u_f]', firstprivate='[q_idx, &
                                & tau_loc, mom_loc]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        lag_div_q = 0._wp
                        lag_div_u = 0._wp
                        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS, BC_COORDS, LAST in &
                                [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'k, l', 'm'), &
                                 (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'j, l', 'n'), &
                                 (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'j, k', 'p')]
                            #:set SV = STENCIL_VAR
                            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
                            if (num_dims >= ${NORM_DIR}$) then
                                ! High face: between this cell and the next
                                u_f = u_face(j, k, l, ${NORM_DIR}$)
                                ! The correction
                                if (mom_loc) u_f = u_f - f_face_dvel_p(rhostar(j, k, l), rhostar(${SF(' + 1')}$), &
                                    & 0.5_wp*(d${XYZ}$(${SV}$) + d${XYZ}$(${SV}$ + 1)), real(p_solved(j, k, l), wp), &
                                    & real(p_solved(${SF(' + 1')}$), wp), tau_loc)
                                if (${SV}$ == ${LAST}$) then
                                    if (f_wall_face(bc_face_${XYZ}$(${BC_COORDS}$, 2))) u_f = 0._wp
                                end if
                                du_face = u_f - u0_face(j, k, l, ${NORM_DIR}$)
                                lag_div_q = lag_div_q + du_face*real(merge(q_cons_vf(q_idx)%sf(j, k, l), &
                                                                     & q_cons_vf(q_idx)%sf(${SF(' + 1')}$), u_f >= 0._wp), &
                                                                     & wp)/d${XYZ}$(${SV}$)
                                lag_div_u = lag_div_u + du_face/d${XYZ}$(${SV}$)
                                ! Low face: between the previous cell and this one
                                u_f = u_face(${SF(' - 1')}$, ${NORM_DIR}$)
                                if (mom_loc) u_f = u_f - f_face_dvel_p(rhostar(${SF(' - 1')}$), rhostar(j, k, l), &
                                    & 0.5_wp*(d${XYZ}$(${SV}$ - 1) + d${XYZ}$(${SV}$)), real(p_solved(${SF(' - 1')}$), wp), &
                                    & real(p_solved(j, k, l), wp), tau_loc)
                                if (${SV}$ == 0) then
                                    if (f_wall_face(bc_face_${XYZ}$(${BC_COORDS}$, 1))) u_f = 0._wp
                                end if
                                du_face = u_f - u0_face(${SF(' - 1')}$, ${NORM_DIR}$)
                                lag_div_q = lag_div_q - du_face*real(merge(q_cons_vf(q_idx)%sf(${SF(' - 1')}$), &
                                                                     & q_cons_vf(q_idx)%sf(j, k, l), u_f >= 0._wp), &
                                                                     & wp)/d${XYZ}$(${SV}$)
                                lag_div_u = lag_div_u - du_face/d${XYZ}$(${SV}$)
                            end if
                        #:endfor
                        if (lag_var > num_fluids .and. lag_var <= num_fluids + n_alpha) lag_div_q = lag_div_q &
                            & - real(q_cons_vf(q_idx)%sf(j, k, l), wp)*lag_div_u
                        p_rate(j, k, l) = lag_div_q
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', firstprivate='[q_idx, tau_loc]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(q_idx)%sf(j, k, l) = real(real(q_cons_vf(q_idx)%sf(j, k, l), wp) - tau_loc*p_rate(j, k, l), stp)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_mass_lag

    !> Face conductances of the pressure operator, with the operator's own face density (f_cond). A boundary face carries one only
    !! where it couples to something: a rank or periodic seam, or a pressure outlet's ghost at pres_out; elsewhere it is a wall
    subroutine s_projection_conductances()

        integer  :: j, k, l
        real(wp) :: area
        logical  :: coupled

        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS, BC_COORDS, LAST in &
                [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'k, l', 'm'), &
                 (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'j, l', 'n'), &
                 (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'j, k', 'p')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (num_dims >= ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, area, coupled]')
                do l = face_bounds(3, ${NORM_DIR}$)%beg, face_bounds(3, ${NORM_DIR}$)%end
                    do k = face_bounds(2, ${NORM_DIR}$)%beg, face_bounds(2, ${NORM_DIR}$)%end
                        do j = face_bounds(1, ${NORM_DIR}$)%beg, face_bounds(1, ${NORM_DIR}$)%end
                            coupled = .true.
                            if (${SV}$ == -1) coupled = seam_lo(${NORM_DIR}$) .or. (has_pres_out(${NORM_DIR}$) &
                                & .and. bc_face_${XYZ}$(${BC_COORDS}$, 1) == BC_GHOST_EXTRAP)
                            if (${SV}$ == ${LAST}$) coupled = seam_hi(${NORM_DIR}$) .or. (has_pres_out(${NORM_DIR}$) &
                                & .and. bc_face_${XYZ}$(${BC_COORDS}$, 2) == BC_GHOST_EXTRAP)
                            cond_${XYZ}$(j, k, l) = 0._wp
                            if (coupled) then
                                #:if NORM_DIR == 1
                                    area = 1._wp
                                    if (num_dims > 1) area = dy(k)
                                    if (num_dims > 2) area = area*dz(l)
                                #:elif NORM_DIR == 2
                                    area = dx(j)
                                    if (num_dims > 2) area = area*dz(l)
                                #:else
                                    area = dx(j)*dy(k)
                                #:endif
                                cond_${XYZ}$(j, k, l) = f_cond(rhostar(j, k, l), rhostar(${SF(' + 1')}$), area, &
                                             & 0.5_wp*(d${XYZ}$(${SV}$) + d${XYZ}$(${SV}$ + 1)))
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_conductances

    !> Solve the pressure system for p_solved, with filled ghosts, from the initial guess in p_guess
    impure subroutine s_pressure_solve(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp)                                                   :: res
        integer                                                    :: iters, j, k, l

        call nvtxStartRange("PROJ-MG-BUILD")
        call s_projection_conductances()
        call s_mg_set_operator(p_diag, cond_x, cond_y, cond_z)
        call nvtxEndRange
        if (any(has_pres_out)) call s_pout_rhs()
        call s_mg_solve(p_sys_rhs, p_guess, iters, res)
        proj_pcg_iters = proj_pcg_iters + iters
        proj_pcg_res = max(proj_pcg_res, res)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    p_solved(j, k, l) = real(p_guess(j, k, l), stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_populate_F_igr_buffers(bc_type, p_solved_sf)
        if (any(has_pres_out)) call s_pout_ghosts()

    end subroutine s_pressure_solve

    !> A pressure outlet's fixed ghost enters the right-hand side as its face's conductance times pres_out
    subroutine s_pout_rhs()

        integer :: j, k, l

        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS, BC_COORDS in &
                [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'k, l'), &
                 (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'j, l'), &
                 (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'j, k')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (has_pres_out(${NORM_DIR}$)) then
                #:for S, IP, FACE in [(1, 0, ' - 1'), (2, 1, '')]
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = plane_bounds(3, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(3, ${NORM_DIR}$, ${IP}$)%end
                        do k = plane_bounds(2, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(2, ${NORM_DIR}$, ${IP}$)%end
                            do j = plane_bounds(1, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(1, ${NORM_DIR}$, ${IP}$)%end
                                if (bc_face_${XYZ}$(${BC_COORDS}$, ${S}$) == BC_GHOST_EXTRAP) then
                                    p_sys_rhs(j, k, l) = p_sys_rhs(j, k, l) + cond_${XYZ}$(${SF(FACE)}$)*pres_out(${NORM_DIR}$)
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:endfor
            end if
        #:endfor

    end subroutine s_pout_rhs

    !> A pressure outlet's ghost holds pres_out, so the face correction pushes flow out against it
    subroutine s_pout_ghosts()

        integer :: j, k, l

        #:for NORM_DIR, XYZ, BC_COORDS in [(1, 'x', 'k, l'), (2, 'y', 'j, l'), (3, 'z', 'j, k')]
            if (has_pres_out(${NORM_DIR}$)) then
                #:for S, IP in [(1, -1), (2, 2)]
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = plane_bounds(3, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(3, ${NORM_DIR}$, ${IP}$)%end
                        do k = plane_bounds(2, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(2, ${NORM_DIR}$, ${IP}$)%end
                            do j = plane_bounds(1, ${NORM_DIR}$, ${IP}$)%beg, plane_bounds(1, ${NORM_DIR}$, ${IP}$)%end
                                if (bc_face_${XYZ}$(${BC_COORDS}$, ${S}$) == BC_GHOST_EXTRAP) p_solved(j, k, &
                                    & l) = real(pres_out(${NORM_DIR}$), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:endfor
            end if
        #:endfor

    end subroutine s_pout_ghosts

    impure subroutine s_finalize_projection_module()

        $:GPU_EXIT_DATA(detach='[p_solved_sf(1)%sf]')
        @:DEALLOCATE(u_face, u0_face, divu_face, p_rate, p_stage, p_step0, face_flux, rhostar, p_diag, p_sys_rhs, p_guess)
        @:DEALLOCATE(cond_x, cond_y, cond_z, p_solved, curv, grad_alpha, bc_face_x, bc_face_y, bc_face_z)
        call s_finalize_multigrid_module()

    end subroutine s_finalize_projection_module

end module m_projection
