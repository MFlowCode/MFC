!>
!! @file
!! @brief Contains module m_ibm

#:include 'macros.fpp'

!> @brief Ghost-node immersed boundary method: locates ghost/image points, computes interpolation coefficients, and corrects the
!! flow state
module m_ibm

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_eos
    use m_helper
    use m_helper_basic
    use m_constants
    use m_compute_levelset
    use m_ib_patches
    use m_model
    use m_patch_geometries
    use m_collisions
    use m_riemann_state, only: flux_rsx_vf, flux_src_rsx_vf
    use m_thermochem, only: num_species, gas_constant, molecular_weights, get_mixture_molecular_weight, get_mixture_energy_mass, &
        & get_mixture_thermal_conductivity_mixavg, get_species_mass_diffusivities_mixavg

    use m_surface_thermochem, only: get_surface_net_production_rates, get_surface_reaction_heat_flux
    use m_thermochem_state, only: get_species_enthalpies_mass

    implicit none

    private :: s_compute_image_points, s_compute_interpolation_coeffs, s_interpolate_image_point, s_find_ghost_points, &
        & s_find_num_ghost_points, s_compute_ghost_point_pressure, s_compute_ghost_point_velocity, s_interpolation_weights, &
        & s_wall_fill, f_tsc_apply, f_quad_fill
    ; public :: s_initialize_ibm_module, s_ibm_setup, s_ibm_correct_state, s_finalize_ibm_module, s_report_ibm_surface, &
        & s_write_ib_surface, s_ib_flux_wall_prepare, s_update_ib_temperatures

    !> Ghost points at which the reacting-surface Newton solve did not reach its tolerance, so the point fell back to a chemically
    !! inert wall (keeping a prescribed Twall). Counted because that fallback is otherwise indistinguishable from a surface
    !! mechanism that simply does nothing, and a run can look converged while the surface chemistry never engaged.
    integer :: n_surface_not_converged = 0

    !> Ghost points whose levelset distance was not positive, leaving the surface gradient (X_IP - X_s)/d undefined. Separate from
    !! the counter above because this is a grid/geometry degeneracy, not a kinetics failure.
    integer :: n_surface_ill_posed = 0

    !> Quadratic wall-normal fill of one ghost point (ib_wall_order = 2): TSC stencils of two probes along the surface normal and
    !! the Lagrange weights that carry the profile through them back to the ghost (see s_wall_fill).
    type :: wall_fill
        integer                   :: tier  !< 1 quadratic, 2 linear mirror (the outer probe lacks fluid), 3 the default treatment
        real(wp), dimension(3)    :: lw  !< Lagrange weights on (phi_w, phi_1, phi_2)
        real(wp), dimension(2)    :: llin  !< Linear mirror on (phi_w, phi_1), which bounds the quadratic
        integer, dimension(3)     :: ic1, ic2  !< TSC stencil centres of the probes
        real(wp), dimension(3, 3) :: wd1, wd2  !< TSC weights of the probes
    end type wall_fill

    !> Ghost-point updates per fill tier, summed over the run and reported at its end, because a tight packing can sit mostly in
    !! tier 2 and a second-order claim is then not true of that surface.
    integer, dimension(3)       :: n_fill_tier = 0
    type(integer_field), public :: ib_markers
    $:GPU_DECLARE(create='[ib_markers]')

    !> 1 at ghost points whose image point reaches the fluid, 0 elsewhere. Ghost points with a buried image point average over these
    !! neighbors, and never over each other.
    type(integer_field) :: corrected_gps
    $:GPU_DECLARE(create='[corrected_gps]')

    type(ghost_point), dimension(:), allocatable :: ghost_points

    !> New primitive state per ghost point, staged so it is published only after every image-point stencil has been read
    !! (MFlowCode/MFC#1901). Writing q_prim inside the same loop that reads it at neighboring image-point stencils makes a ghost's
    !! result depend on whether its neighbors were updated first. Inverse-distance weights mask solid cells, so the race is latent
    !! there; a stencil that reads ghost cells (ib_ip_interp = 2) needs it.
    real(wp), dimension(:,:), allocatable :: gp_prim_new
    $:GPU_DECLARE(create='[ghost_points]')

    !> Surface record per ghost point from the latest s_ibm_correct_state, kept when ib_surface_wrt: (1) area weight, (2) wall
    !! temperature, (3) gasified mass flux. Only ghost points within two cells of the surface carry a weight, so summing weight*flux
    !! over them integrates over the surface; see s_record_gp_surface.
    integer, parameter                    :: ib_surf_nvars = 3
    real(wp), allocatable, dimension(:,:) :: gp_surf
    $:GPU_DECLARE(create='[gp_surf]')

    !> ib_flux_wall: per ghost point with a fluid face neighbor of a reacting body, from its surface solve, the species mass fluxes
    !! into the gas (1:num_species), the energy flux (ibw_e), the blowing velocity (ibw_u:ibw_u + 2) and the wall area per face
    !! area, 1/sum|n_i| (ibw_a, 0 off a reacting wall); ib_wall_theta: each fluid cell's positivity factor on its wall faces' rates
    integer, parameter, public :: ibw_e = num_species + 1, ibw_u = num_species + 2
    integer, parameter, public :: ibw_a = num_species + 5
    ! Not allocatable: the AMD OpenMP kernels of m_riemann_solver_hllc assume every allocatable they name is present, and these
    ! are only filled under ib_flux_wall (their sf stay unassociated otherwise, as alf_sum's do in m_rhs)
    type(scalar_field), dimension(ibw_a), public :: ib_wall_src
    type(scalar_field), public                   :: ib_wall_theta
    $:GPU_DECLARE(create='[ib_wall_src, ib_wall_theta]')

    !> gp_surf is filled when it is written (ib_surface_wrt) or integrated (ib_exchange).
    logical         :: record_surface
    logical, public :: ib_exchange = .false.  !< ib_heat is accumulated: some IB has thermal_bc = 3, or ib_heat_wrt
    $:GPU_DECLARE(create='[record_surface]')

    !> Per global IB, for the thermal_bc = 3 update and ib_heat_wrt: (1) surface area and (4) gasification rate, from gp_surf; (2)
    !! energy into the solid and (3) mass out of it, both summed over the faces between the body's cells and fluid cells with the RK
    !! stage weights, so over a step they are exactly the energy and mass the fluid exchanged with the body.
    real(wp), allocatable, dimension(:,:) :: ib_heat
    $:GPU_DECLARE(create='[ib_heat]')

    integer :: num_gps  !< Number of ghost points
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[gp_layers, num_gps]')
#elif defined(MFC_OpenMP)
    $:GPU_DECLARE(create='[num_gps]')
#endif
    logical :: moving_immersed_boundary_flag
    !> Cells per direction in the image-point interpolation stencil: 2 (inverse distance, bilinear) or 3 (TSC)
    integer :: ib_ip_width
    logical :: ib_keep_ghost  !< The interpolation or the wall probes read solid cells, so ghost cells keep the last correction
    logical :: ib_wall_quad   !< ib_wall_order = 2: quadratic wall-normal ghost velocity at no-slip walls
    $:GPU_DECLARE(create='[ib_ip_width, ib_keep_ghost, ib_wall_quad]')
    logical :: centroid_offsets_active = .false.  !< some patch is a moving airfoil or STL, which carries a centroid_offset

    ! IB force reduction (s_communicate_ib_forces): distinct neighborhood ranks in ascending order, and buffers allocated once
    integer, parameter    :: ib_rec_len = 7  !< gbl_patch_id, force(1:3), torque(1:3)
    integer               :: n_ib_nbrs = 0
    integer, allocatable  :: ib_nbrs(:), ib_send_counts(:), ib_reqs(:), ib_stats(:,:)
    real(wp), allocatable :: ib_send_bufs(:,:), ib_recv_bufs(:,:)
    private :: ib_rec_len, n_ib_nbrs, ib_nbrs, ib_send_counts, ib_reqs, ib_stats, ib_send_bufs, ib_recv_bufs

    !> Force (1:3) and torque (4:6) on each local IB from the fluid-body face fluxes of the last RHS evaluation
    real(wp), allocatable, dimension(:,:) :: ib_face_ft
    $:GPU_DECLARE(create='[ib_face_ft]')

contains

    !> Allocates memory for the variables in the IBM module
    impure subroutine s_initialize_ibm_module()

        if (p > 0) then
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, -buff_size:p+buff_size))
            @:ALLOCATE(corrected_gps%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, -buff_size:p+buff_size))
        else
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, 0:0))
            @:ALLOCATE(corrected_gps%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, 0:0))
        end if

        @:ACC_SETUP_SFs(ib_markers)
        @:ACC_SETUP_SFs(corrected_gps)

        $:GPU_ENTER_DATA(copyin='[num_gps]')
        ib_ip_width = merge(3, 2, ib_ip_interp == 2)
        ib_wall_quad = ib_wall_order == 2
        ib_keep_ghost = ib_ip_interp == 2 .or. ib_wall_quad
        $:GPU_UPDATE(device='[ib_second_order_vel, ib_ip_min_dist, ib_ip_interp, ib_ip_width, ib_keep_ghost, ib_wall_quad]')

        if (collision_model > 0) call s_initialize_collisions_module()

    end subroutine s_initialize_ibm_module

    !> Initializes the values of various IBM variables, such as ghost points and image points.
    impure subroutine s_ibm_setup()

        integer         :: i, j, k, gid
        integer(kind=8) :: n_need
        real(wp)        :: t_init  !< initial time for prescribed kinematics
        real(wp)        :: max_num_gps_rank
        integer(kind=8) :: max_num_gps
        integer(kind=8) :: n_lumped

        call nvtxStartRange("SETUP-IBM-MODULE")

        ! GPU routines require updated cell centers
        $:GPU_UPDATE(device='[num_ibs, num_gbl_ibs, x_cc, y_cc, dx, dy, ib_bc_x%beg, ib_bc_x%end, ib_bc_y%beg, ib_bc_y%end]')
        if (p /= 0) then
            $:GPU_UPDATE(device='[z_cc, dz, ib_bc_z%beg, ib_bc_z%end]')
        end if
        $:GPU_UPDATE(device='[patch_ib(1:num_ibs), glb_bounds]')

        ! do all set up for moving immersed boundaries; prescribed kinematics are evaluated at the initial time so the
        ! first stage already sees the correct body state, velocity and angular velocity
        t_init = mytime
        $:GPU_PARALLEL_LOOP(private='[i]', copyin='[t_init]')
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) then
                if (patch_ib(i)%kin_model > 0) call s_prescribed_kinematics(i, t_init)
                call s_compute_moment_of_inertia(patch_ib(i), patch_ib(i)%angular_vel, patch_ib(i)%moment)
            end if
            call s_update_ib_rotation_matrix(i)
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_UPDATE(host='[patch_ib(1:num_ibs)]')

        ! allocate some arrays for MPI communication, if required by this simulation
#ifdef MFC_MPI
        if (num_procs > 1) call s_setup_ib_force_comm()
#endif

        @:ALLOCATE(ib_face_ft(6, num_local_ibs_max))
        ib_face_ft = 0._wp
        $:GPU_UPDATE(device='[ib_face_ft]')

        call s_update_ib_lookup()

        ! recompute the new ib_patch locations
        ib_markers%sf = 0._wp
        corrected_gps%sf = 0
        $:GPU_UPDATE(device='[ib_markers%sf, corrected_gps%sf]')
        call s_apply_ib_patches(ib_markers)
        $:GPU_UPDATE(host='[ib_markers%sf]')
        call s_check_every_patch_marked()
        ! One reduction decides whether any patch needs an offset, so cases without one (particle clouds have thousands of
        ! patches) skip the per-patch collectives. Then loop over global ids: ranks hold different patches.
        call s_mpi_allreduce_integer_sum(int(count([(f_needs_centroid_offset(patch_ib(i)), i=1, num_ibs)]), 8), n_need)
        centroid_offsets_active = n_need > 0_8
        if (centroid_offsets_active) then
            do gid = 1, num_gbl_ibs
                call s_get_neighborhood_idx(gid, i)
                call s_compute_centroid_offset(gid, i)
            end do
            call s_restore_centroid_offsets(t_init)
        end if
        do i = 1, num_ibs
            $:GPU_UPDATE(device='[patch_ib(i)]')
        end do

        ! find the number of ghost points and set them to be the maximum total across ranks
        call s_find_num_ghost_points(num_gps)
        if (moving_immersed_boundary_flag) then
            ! Every rank allocates the same size, which must hold any one rank's peak as the body moves: twice the global
            ! total, but at least eight times the largest per-rank count (a thin body crossing few of many ranks)
            call s_mpi_allreduce_integer_sum(int(num_gps, 8), max_num_gps)
            call s_mpi_allreduce_max(real(num_gps, wp), max_num_gps_rank)
            max_num_gps = max(max_num_gps*2_8, 8_8*int(max_num_gps_rank, 8))
            max_num_gps = min(max_num_gps, int(m + 1, 8)*int(n + 1, 8)*int(p + 1, 8))
        else
            max_num_gps = int(num_gps, 8)
        end if

        ! set the size of the ghost point arrays to be the amount of points total, plus a factor of 2 buffer
        $:GPU_UPDATE(device='[num_gps]')
        @:ALLOCATE(ghost_points(1:max_num_gps))
        @:ALLOCATE(gp_prim_new(1:sys_size, 1:max_num_gps))

        $:GPU_ENTER_DATA(copyin='[ghost_points, gp_prim_new]')
        ! Global, not per rank: s_update_ib_temperatures is collective, and a rank whose neighborhood holds no thermal_bc = 3
        ! body must still take part in its allreduce.
        call s_mpi_allreduce_integer_sum(int(count(patch_ib(1:num_ibs)%thermal_bc == 3), 8), n_lumped)
        ib_exchange = n_lumped > 0 .or. ib_heat_wrt
        record_surface = ib_surface_wrt .or. ib_exchange
        $:GPU_UPDATE(device='[record_surface]')
        if (record_surface) then
            @:ALLOCATE(gp_surf(ib_surf_nvars, 1:max_num_gps))
            gp_surf = 0._wp
            $:GPU_UPDATE(device='[gp_surf]')
        end if
        if (ib_exchange) then
            @:ALLOCATE(ib_heat(4, num_gbl_ibs))
            ib_heat = 0._wp
            $:GPU_UPDATE(device='[ib_heat]')
        end if
        ! Ghost-cell IBM, Tseng & Ferziger JCP (2003), Mittal & Iaccarino ARFM (2005)
        call s_find_ghost_points()
        call s_apply_levelset(ghost_points, num_gps)

        call s_compute_image_points()
        call s_compute_interpolation_coeffs()
        if (ib_flux_wall) call s_ib_flux_wall_setup()

        call nvtxEndRange

    end subroutine s_ibm_setup

    !> Ghost pressure from the wall-normal momentum balance, dp/dn = -rho (a_w . n) with n out of the body and a_w the acceleration
    !! of the wall point: the body's translation (force/mass for a two-way coupled body) and the centripetal acceleration omega x
    !! (omega x r) of its rotation about its centroid (a steady spin; angular acceleration is not included). A fluid held at a
    !! rotating wall's velocity needs this pressure gradient to stay on its circular path; mirroring the image-point pressure
    !! instead leaves the wall layer unbalanced. Clamped both ways: the linearization holds only while the correction is order one,
    !! and an unbounded one drives the ghost state to vacuum.
    subroutine s_compute_ghost_point_pressure(gp, gp_patch_id, radial_vector, rho, pres_IP, pres_GP)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(in)      :: gp
        integer, intent(in)                :: gp_patch_id
        real(wp), dimension(3), intent(in) :: radial_vector  !< Centroid to ghost point
        real(wp), intent(in)               :: rho, pres_IP
        real(wp), intent(out)              :: pres_GP
        real(wp), dimension(3)             :: accel, norm, r_w, w_x_r, w_x_w_x_r

        norm = gp%levelset_norm/max(sqrt(sum(gp%levelset_norm**2)), sgm_eps)
        accel = 0._wp
        if (patch_ib(gp_patch_id)%moving_ibm > 1) accel = patch_ib(gp_patch_id)%force/patch_ib(gp_patch_id)%mass
        r_w = radial_vector + abs(real(gp%levelset, kind=wp))*norm  ! the wall point, relative to the centroid
        call s_cross_product(patch_ib(gp_patch_id)%angular_vel, r_w, w_x_r)
        call s_cross_product(patch_ib(gp_patch_id)%angular_vel, w_x_r, w_x_w_x_r)

        pres_GP = pres_IP/min(max(1._wp - (abs(gp%levelset) + gp%ip_dist)*rho/pres_IP*dot_product(accel + w_x_w_x_r, norm), &
                              & 5.e-1_wp), 2._wp)

    end subroutine s_compute_ghost_point_pressure

    subroutine s_compute_ghost_point_velocity(gp, gp_patch_id, radial_vector, vel_IP, pres_IP, vel_GP, q_prim_vf, wf)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(in)                       :: gp
        integer, intent(in)                                 :: gp_patch_id
        real(wp), dimension(3), intent(in)                  :: radial_vector, vel_IP
        real(wp), intent(in)                                :: pres_IP
        real(wp), dimension(3), intent(out)                 :: vel_GP
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf  !< Probed by the quadratic fill
        type(wall_fill), intent(in)                         :: wf
        real(wp), dimension(3)                              :: norm, vel_norm_IP, rotation_velocity
        real(wp)                                            :: buf, v_blow_eff
        integer                                             :: q          !< Iterator variable

        ! Calculate velocity of ghost cell
        if (gp%slip) then
            norm(1:3) = gp%levelset_norm
            buf = sqrt(sum(norm**2))
            norm = norm/buf
            vel_norm_IP = sum(vel_IP*norm)*norm
            vel_GP = vel_IP - vel_norm_IP
            if (patch_ib(gp_patch_id)%moving_ibm /= 0) then
                ! compute the linear velocity of the ghost point due to rotation
                call s_cross_product(patch_ib(gp_patch_id)%angular_vel, radial_vector, rotation_velocity)

                ! add only the component of the IB's motion that is normal to the surface
                vel_GP = vel_GP + sum((patch_ib(gp_patch_id)%vel + rotation_velocity)*norm)*norm
            end if
        else
            if (patch_ib(gp_patch_id)%moving_ibm == 0) then
                ! we know the object is not moving if moving_ibm is 0 (false)
                vel_GP = 0._wp
            else
                ! convert the angular velocity from the inertial reference frame to the fluids frame, then convert to linear
                ! velocity
                call s_cross_product(patch_ib(gp_patch_id)%angular_vel, radial_vector, rotation_velocity)
                do q = 1, 3
                    ! if mibm is 1 or 2, then the boundary may be moving
                    vel_GP(q) = patch_ib(gp_patch_id)%vel(q)  ! add the linear velocity
                    vel_GP(q) = vel_GP(q) + rotation_velocity(q)  ! add the rotational velocity
                end do
            end if
        end if

        ! Burning/injecting surface: superimpose wall-normal (outward) blowing on the
        ! ghost velocity so the immersed surface transpires/injects gas into the flow.
        if (patch_ib(gp_patch_id)%v_blow > 0._wp) then
            v_blow_eff = patch_ib(gp_patch_id)%v_blow
            ! Pressure-coupled burn rate (Vieille's law r_dot ~ p^n)
            if (patch_ib(gp_patch_id)%burn_rate_pref > 0._wp) then
                ! max(pres_IP, 0) guards the fractional power against a transient negative
                v_blow_eff = v_blow_eff*(max(pres_IP, &
                                         & 0._wp)/patch_ib(gp_patch_id)%burn_rate_pref)**patch_ib(gp_patch_id)%burn_rate_exp
            end if
            norm(1:3) = gp%levelset_norm
            buf = sqrt(sum(norm**2))
            if (buf > 0._wp) vel_GP = vel_GP + v_blow_eff*norm/buf
        end if

        if (wf%tier < 3) then
            ! Quadratic profile through the wall value and the two probes. The ghost value is what enforces no-slip in the solve,
            ! so its order sets the order of the wall shear and the drag.
            do q = 1, num_dims
                vel_GP(q) = f_quad_fill(vel_GP(q), f_tsc_apply(q_prim_vf(eqn_idx%mom%beg + q - 1), wf%ic1, wf%wd1), &
                       & f_tsc_apply(q_prim_vf(eqn_idx%mom%beg + q - 1), wf%ic2, wf%wd2), wf%lw, wf%llin)
            end do
        else if (ib_second_order_vel) then
            ! extrapolate linearly through the boundary intercept, Mittal et al. (2008)
            vel_GP = vel_GP - abs(gp%levelset)/max(gp%ip_dist, sgm_eps)*(vel_IP - vel_GP)
        end if

    end subroutine s_compute_ghost_point_velocity

    !> Update the conservative variables at the ghost points
    subroutine s_ibm_correct_state(q_cons_vf, q_prim_vf, pb_in, mv_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf  !< Primitive Variables
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf  !< Primitive Variables
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        integer :: i, j, k, l, q, r  !< Iterator variables
        integer :: jj, kk, ll  !< Neighbor iterators
        integer :: rad, rad_z, num_nbrs  !< Neighbor stencil radius and population
        integer :: patch_id, patch_id_temp  !< Patch ID of ghost point
        real(wp) :: rho, gamma, pi_inf, dyn_pres  !< Mixture variables
        real(wp) :: vel_sum_g, E_ghost  !< Ghost-point velocity magnitude and energy
        real(wp), dimension(2) :: Re_K
        real(wp) :: G_K
        real(wp) :: qv_K
        real(wp) :: pres_IP, pres_GP
        real(wp), dimension(3) :: vel_IP
        real(wp) :: c_IP
        type(wall_fill) :: wf  !< Quadratic wall-normal fill of the current ghost point
        integer :: nf1, nf2, nf3  !< Per-call fill-tier tallies

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)  :: Gs
            real(wp), dimension(3)  :: alpha_rho_IP, alpha_IP, alpha_rho_GP
            real(wp), dimension(3)  :: r_IP, v_IP, pb_IP, mv_IP
            real(wp), dimension(18) :: nmom_IP
            real(wp), dimension(12) :: presb_IP, massv_IP
        #:else
            real(wp), dimension(num_fluids) :: Gs
            real(wp), dimension(num_fluids) :: alpha_rho_IP, alpha_IP, alpha_rho_GP
            real(wp), dimension(nb)         :: r_IP, v_IP, pb_IP, mv_IP
            real(wp), dimension(nb*nmom)    :: nmom_IP
            real(wp), dimension(nb*nnode)   :: presb_IP, massv_IP
        #:endif
        real(wp), dimension(${NUM_SPECIES}$) :: Ys_IP, Ys_g, Ys_s, W_species
        real(wp) :: alpha_q, alpha_rho_q, e_q
        real(wp) :: T_IP, mw_IP, e_IP  !< Image-point temperature, mixture MW, and mass-specific internal energy (chemistry)
        real(wp) :: T_s, T_g, mw_s, mw_g, rho_s, rho_g, mdot_s, v_stefan, d
        real(wp), dimension(${NUM_SPECIES}$) :: mdot_k, h_w  !< ib_flux_wall: species mass fluxes into the gas; species enthalpies
        real(wp) :: e_s, nrm_l  !< ib_flux_wall: energy flux into the gas; length of the ghost's normal
        logical :: surface_converged
        integer :: n_not_converged, n_ill_posed  !< Per-call reacting-surface failure tallies (see the module-level counters)
        !> Host copies; amdflang reads Re_size stale cross-TU
        integer                :: Re_size_loc1, Re_size_loc2
        real(wp), dimension(3) :: norm  !< Levelset normal at the ghost point, normalized below with buf
        ! Primitive variables at the image point associated with a ghost point, interpolated from surrounding fluid cells.

        real(wp), dimension(3) :: physical_loc   !< Physical loc of GP
        real(wp), dimension(3) :: vel_g          !< Velocity of GP
        real(wp), dimension(3) :: radial_vector  !< vector from centroid to ghost point
        real(wp)               :: nbub
        real(wp)               :: buf, buf_prim
        type(ghost_point)      :: gp

        ! set the Moving IBM interior conservative variables. A solid cell that is a ghost keeps the value of the last correction
        ! when the interpolation reads solid cells (ib_keep_ghost): resetting it to the body velocity, and without chemistry to a
        ! pressure of 1 Pa, would feed those values into the image-point state. Inverse-distance weights drop solid cells, so they
        ! never see it.
        $:GPU_PARALLEL_LOOP(private='[i, j, k, patch_id, rho, patch_id_temp]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    patch_id = ib_markers%sf(j, k, l)
                    if (patch_id /= 0) then
                        call s_decode_patch_periodicity(patch_id, patch_id_temp)
                        call s_get_neighborhood_idx(patch_id_temp, patch_id)
                        if (patch_id > 0 .and. .not. (ib_keep_ghost .and. corrected_gps%sf(j, k, l) /= 0)) then
                            ! skip pressure correction with chemistry to prevent unphysical pressure
                            if (.not. chemistry) q_prim_vf(eqn_idx%E)%sf(j, k, l) = 1._wp
                            rho = 0._wp
                            do i = 1, num_fluids
                                rho = rho + q_prim_vf(eqn_idx%cont%beg + i - 1)%sf(j, k, l)
                            end do

                            ! Sets the momentum
                            do i = 1, num_dims
                                q_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = patch_ib(patch_id)%vel(i)*rho
                                q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = patch_ib(patch_id)%vel(i)
                            end do
                        end if  ! patch_id > 0
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        n_not_converged = 0
        n_ill_posed = 0
        nf1 = 0; nf2 = 0; nf3 = 0
        Re_size_loc1 = Re_size(1); Re_size_loc2 = Re_size(2)

        if (num_gps > 0) then
            ! Seed the staging buffer, so a component no ghost point writes survives the flush unchanged.
            $:GPU_PARALLEL_LOOP(private='[i, q, j, k, l, gp]', present='[ghost_points, gp_prim_new]')
            do i = 1, num_gps
                gp = ghost_points(i)
                j = gp%loc(1); k = gp%loc(2); l = gp%loc(3)
                do q = 1, sys_size
                    gp_prim_new(q, i) = q_prim_vf(q)%sf(j, k, l)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            $:GPU_PARALLEL_LOOP(private='[i, physical_loc, dyn_pres, alpha_rho_IP, alpha_IP, alpha_rho_GP, pres_IP, pres_GP, &
                                & vel_IP, vel_g, r_IP, v_IP, pb_IP, mv_IP, nmom_IP, presb_IP, massv_IP, rho, gamma, pi_inf, Re_K, &
                                & G_K, Gs, gp, norm, buf, radial_vector, j, k, l, q, qv_K, c_IP, nbub, patch_id, Ys_IP, &
                                & W_species, T_IP, mw_IP, e_IP, Ys_g, Ys_s, T_s, T_g, mw_s, mw_g, rho_s, rho_g, mdot_s, v_stefan, &
                                & d, surface_converged, vel_sum_g, E_ghost, alpha_q, alpha_rho_q, e_q, mdot_k, h_w, e_s, nrm_l, &
                                & wf]', reduction='[[n_not_converged, n_ill_posed, nf1, nf2, nf3]]', reductionOp='[+]', &
                                & present='[ghost_points]', firstprivate='[Re_size_loc1, Re_size_loc2]')
            do i = 1, num_gps
                if (record_surface) then
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, ib_surf_nvars
                        gp_surf(q, i) = 0._wp
                    end do
                end if
                gp = ghost_points(i)
                if (.not. gp%interp_valid) cycle
                j = gp%loc(1)
                k = gp%loc(2)
                l = gp%loc(3)
                patch_id = ghost_points(i)%ib_patch_id

                ! Calculate physical location of GP
                physical_loc = [x_cc(j), y_cc(k), 0._wp]
                if (num_dims == 3) physical_loc(3) = z_cc(l)

                ! Interpolate primitive variables at image point associated w/ GP
                if (bubbles_euler .and. .not. qbmm) then
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, &
                                                   & pb_IP, mv_IP)
                else if (qbmm .and. polytropic) then
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, &
                                                   & pb_IP, mv_IP, nmom_IP)
                else if (qbmm .and. .not. polytropic) then
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, &
                                                   & pb_IP, mv_IP, nmom_IP, pb_in, mv_in, presb_IP, massv_IP)
                else if (chemistry) then
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, Ys_IP=Ys_IP)
                else
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP)
                end if

                ! Injecting (burning) surface: replace the mirrored ghost composition with pure
                ! injected fuel at the local pressure and the ambient (image-point) temperature.
                ! Setting a consistent injected density here (rather than reusing the heavy ambient
                ! rho) keeps the light fuel at a physical temperature and feeds the surface flame.
                if (chemistry .and. patch_ib(patch_id)%inj_species > 0) then
                    call get_mixture_molecular_weight(Ys_IP, mw_IP)
                    T_IP = pres_IP*mw_IP/(alpha_rho_IP(1)*gas_constant)
                    Ys_IP = 0._wp
                    Ys_IP(patch_ib(patch_id)%inj_species) = 1._wp
                    call get_mixture_molecular_weight(Ys_IP, mw_IP)
                    alpha_rho_IP(1) = pres_IP*mw_IP/(T_IP*gas_constant)
                end if

                ! Thermal and heterogeneous reacting-surface boundary conditions.
                v_stefan = 0._wp
                surface_converged = .false.

                if (chemistry .and. patch_ib(patch_id)%inj_species == 0) then
                    ! Intrinsic gas state at the image point:     rho_IP = (alpha*rho)_IP / alpha_IP
                    call get_mixture_molecular_weight(Ys_IP, mw_IP)
                    T_IP = pres_IP*mw_IP*alpha_IP(1)/(alpha_rho_IP(1)*gas_constant)

                    if (patch_ib(patch_id)%surface_reaction == 1) then
                        ! Heterogeneous reacting surface.
                        d = abs(real(gp%levelset, kind=wp))

                        W_species(1:num_species) = molecular_weights(:)

                        ! d is the gas-side gradient length for every surface flux, so a ghost point sitting exactly on the
                        ! immersed surface (grid-aligned rectangle patch, coincident STL facet) leaves the whole balance
                        ! undefined. Left alone the 1/d produces NaN residuals, and a NaN then defeats the pivot test in
                        ! s_solve_surface_linear_system -- NaN <= epsilon is false -- so the solve reports success and writes
                        ! NaN into the ghost state. Such a point is not solvable, so it is recorded and skipped rather than
                        ! having its distance rescaled to something it is not.
                        if (d > 0._wp) then
                            call s_solve_surface(pres_IP, T_IP, patch_ib(patch_id)%Twall, d, Ys_IP, W_species, merge(1, &
                                                 & patch_ib(patch_id)%thermal_bc, f_holds_twall(patch_ib(patch_id)%thermal_bc)), &
                                                 & Ys_s, T_s, mdot_s, surface_converged)
                            if (.not. surface_converged) n_not_converged = n_not_converged + 1
                        else
                            n_ill_posed = n_ill_posed + 1
                        end if

                        if (surface_converged) then
                            call get_mixture_molecular_weight(Ys_s, mw_s)

                            ! Intrinsic gas density at the reacting surface.
                            rho_s = pres_IP*mw_s/(gas_constant*T_s)
                            if (rho_s > 0._wp) v_stefan = mdot_s/rho_s
                        end if

                        ! ib_flux_wall: the rates this ghost's wall faces carry (zero where the solve failed: an inert, impermeable
                        ! face)
                        if (ib_flux_wall) then
                            if (ib_wall_src(ibw_a)%sf(j, k, l) > 0._wp) then
                                mdot_k = 0._wp
                                e_s = 0._wp
                                if (surface_converged .and. rho_s > 0._wp) then
                                    ! W_k omega_k and the energy they carry, sum_k mdot_k (h_k(T_s) + v_n^2/2), per unit wall area
                                    call get_surface_net_production_rates(rho_s, T_s, Ys_s, mdot_k)
                                    call get_species_enthalpies_mass(T_s, h_w)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, num_species
                                        mdot_k(q) = molecular_weights(q)*mdot_k(q)
                                        e_s = e_s + mdot_k(q)*(h_w(q) + 0.5_wp*v_stefan**2)
                                    end do
                                end if
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = 1, num_species
                                    ib_wall_src(q)%sf(j, k, l) = mdot_k(q)
                                end do
                                ib_wall_src(ibw_e)%sf(j, k, l) = e_s
                                nrm_l = max(sqrt(sum(gp%levelset_norm**2)), sgm_eps)
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = 1, 3
                                    ib_wall_src(ibw_u + q - 1)%sf(j, k, l) = v_stefan*gp%levelset_norm(q)/nrm_l
                                end do
                            end if
                        end if
                    end if

                    ! Inert surface, or a reacting surface whose solve failed: zero species flux, so the surface composition
                    ! is the image-point composition. Only the chemistry is dropped on failure; a prescribed Twall is kept,
                    ! and thermal_bc = 2, which has no wall temperature without the energy balance, takes a zero gradient.
                    if (.not. surface_converged) then
                        Ys_s(1:num_species) = Ys_IP(1:num_species)
                        T_s = T_IP
                        if (f_holds_twall(patch_ib(patch_id)%thermal_bc)) T_s = patch_ib(patch_id)%Twall
                    end if

                    if (record_surface) call s_record_gp_surface(i, gp, T_s, mdot_s, surface_converged)

                    call s_blend_ghost_temperature(T_IP, T_s, T_surface_min, T_surface_max, T_g)
                    call s_blend_ghost_species(Ys_IP, Ys_s, Ys_g)

                    call get_mixture_molecular_weight(Ys_g, mw_g)
                    alpha_rho_IP(1) = alpha_IP(1)*pres_IP*mw_g/(gas_constant*T_g)
                end if

                ! Isothermal wall without chemistry: keep p_IP, set the ghost density from the mirrored temperature.
                if (.not. chemistry .and. patch_ib(patch_id)%thermal_bc == 1) then
                    call s_phase_temperature(alpha_rho_IP(1)/alpha_IP(1), pres_IP, 1, T_IP)
                    call s_blend_ghost_temperature(T_IP, patch_ib(patch_id)%Twall, 0._wp, huge(1._wp), T_g)
                    call s_phase_density(pres_IP, T_g, 1, rho_g)
                    alpha_rho_IP(1) = alpha_IP(1)*rho_g
                end if

                ! If in simulation, use acc mixture subroutines
                if (hypoelasticity) then
                    call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv_K, alpha_IP, alpha_rho_IP, Re_K, &
                        & G_K, Gs, Re_size_loc1, Re_size_loc2)
                else
                    call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv_K, alpha_IP, alpha_rho_IP, Re_K, &
                        & Re_size_loc1=Re_size_loc1, Re_size_loc2=Re_size_loc2)
                end if

                if (surface_tension) gp_prim_new(eqn_idx%c, i) = c_IP

                ! get the vector that points from the centroid to the ghost
                radial_vector(1) = physical_loc(1) - (patch_ib(patch_id)%x_centroid + real(ghost_points(i)%x_periodicity, &
                              & wp)*(glb_bounds(1)%end - glb_bounds(1)%beg))
                radial_vector(2) = physical_loc(2) - (patch_ib(patch_id)%y_centroid + real(ghost_points(i)%y_periodicity, &
                              & wp)*(glb_bounds(2)%end - glb_bounds(2)%beg))
                radial_vector(3) = 0._wp
                if (num_dims == 3) radial_vector(3) = physical_loc(3) - (patch_ib(patch_id)%z_centroid &
                    & + real(ghost_points(i)%z_periodicity, wp)*(glb_bounds(3)%end - glb_bounds(3)%beg))

                ! set the pressure and density. A static body (its wall does not move whatever its angular_vel) and a prescribed
                ! body keep the image-point pressure unless the latter spins about its centroid; prescribed kinematics (kin_model >
                ! 0) rotate about a hinge, so the centroid's own acceleration would be needed
                if (patch_ib(patch_id)%moving_ibm == 0 .or. (patch_ib(patch_id)%moving_ibm == 1 &
                    & .and. (all(patch_ib(patch_id)%angular_vel == 0._wp) .or. patch_ib(patch_id)%kin_model > 0))) then
                    pres_GP = pres_IP
                    alpha_rho_GP = alpha_rho_IP
                else
                    call s_compute_ghost_point_pressure(gp, patch_id, radial_vector, rho, pres_IP, pres_GP)

                    ! The adiabatic wall condition T_GP = T_IP the correction is derived from also
                    ! fixes the ghost density: p + B = (n - 1)*cv*rho*T at both points under the one
                    ! temperature leaves each partial density carrying the pressure ratio. The volume
                    ! fractions are untouched, so only rho and qv move with it.
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, num_fluids
                        alpha_rho_GP(q) = alpha_rho_IP(q)*(pres_GP + isentrope_B(q))/(pres_IP + isentrope_B(q))
                    end do
                    call s_compute_mixture_coefficients(alpha_rho_GP, alpha_IP, rho, gamma, pi_inf, qv_K)
                end if
                gp_prim_new(eqn_idx%E, i) = pres_GP

                ! Set q_prim_vf params at GP
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    gp_prim_new(q, i) = alpha_rho_GP(q)
                    gp_prim_new(eqn_idx%adv%beg + q - 1, i) = alpha_IP(q)
                end do

                ! Calculate velocity of ghost cell
                wf%tier = 3
                if (ib_wall_quad .and. .not. gp%slip) then
                    call s_wall_fill(gp, wf)
                    if (wf%tier == 1) nf1 = nf1 + 1
                    if (wf%tier == 2) nf2 = nf2 + 1
                    if (wf%tier == 3) nf3 = nf3 + 1
                end if
                call s_compute_ghost_point_velocity(gp, patch_id, radial_vector, vel_IP, pres_IP, vel_g, q_prim_vf, wf)

                if (chemistry .and. patch_ib(patch_id)%inj_species == 0 .and. patch_ib(patch_id)%surface_reaction == 1 &
                    & .and. surface_converged) then
                    norm(1:3) = gp%levelset_norm
                    buf = sqrt(sum(norm**2))
                    if (buf > 0._wp) vel_g = vel_g + v_stefan*norm/buf
                end if

                ! Set momentum
                vel_sum_g = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = eqn_idx%mom%beg, eqn_idx%mom%end
                    q_cons_vf(q)%sf(j, k, l) = rho*vel_g(q - eqn_idx%mom%beg + 1)
                    gp_prim_new(q, i) = vel_g(q - eqn_idx%mom%beg + 1)
                    vel_sum_g = vel_sum_g + vel_g(q - eqn_idx%mom%beg + 1)**2._wp
                end do
                dyn_pres = 5.e-1_wp*rho*vel_sum_g

                ! Set continuity and adv vars
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    q_cons_vf(q)%sf(j, k, l) = alpha_rho_GP(q)
                    q_cons_vf(eqn_idx%adv%beg + q - 1)%sf(j, k, l) = alpha_IP(q)
                end do

                ! Set color function
                if (surface_tension) q_cons_vf(eqn_idx%c)%sf(j, k, l) = c_IP

                ! Set Energy
                if (chemistry) then
                    ! Use the reconstructed thermal/species ghost state for an inert
                    ! thermal surface or a converged heterogeneous reacting surface.
                    if (patch_ib(patch_id)%inj_species == 0 .and. (patch_ib(patch_id)%surface_reaction == 0 &
                        & .or. surface_converged)) then

                        call get_mixture_energy_mass(T_g, Ys_g, e_IP)
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_species
                            q_cons_vf(eqn_idx%species%beg + q - 1)%sf(j, k, l) = rho*Ys_g(q)
                        end do
                    else
                        ! Ordinary chemistry/injection, or fallback after a failed
                        ! heterogeneous surface solve: retain the image-point state.
                        mw_IP = 0._wp
                        call get_mixture_molecular_weight(Ys_IP, mw_IP)
                        T_IP = pres_IP*mw_IP/(rho*gas_constant)
                        call get_mixture_energy_mass(T_IP, Ys_IP, e_IP)
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_species
                            q_cons_vf(eqn_idx%species%beg + q - 1)%sf(j, k, l) = rho*Ys_IP(q)
                        end do
                    end if

                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = rho*e_IP + dyn_pres
                else
                    call s_compute_energy(pres_GP, alpha_rho_GP, alpha_IP, vel_sum_g, E_ghost)
                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = E_ghost
                end if
                ! Set bubble vars
                if (bubbles_euler .and. .not. qbmm) then
                    call s_comp_n_from_prim(alpha_IP(1), r_IP, nbub, weight)
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        q_cons_vf(eqn_idx%bub%beg + (q - 1)*2)%sf(j, k, l) = nbub*r_IP(q)
                        q_cons_vf(eqn_idx%bub%beg + (q - 1)*2 + 1)%sf(j, k, l) = nbub*v_IP(q)
                        if (.not. polytropic) then
                            q_cons_vf(eqn_idx%bub%beg + (q - 1)*4)%sf(j, k, l) = nbub*r_IP(q)
                            q_cons_vf(eqn_idx%bub%beg + (q - 1)*4 + 1)%sf(j, k, l) = nbub*v_IP(q)
                            q_cons_vf(eqn_idx%bub%beg + (q - 1)*4 + 2)%sf(j, k, l) = nbub*pb_IP(q)
                            q_cons_vf(eqn_idx%bub%beg + (q - 1)*4 + 3)%sf(j, k, l) = nbub*mv_IP(q)
                        end if
                    end do
                end if

                if (qbmm) then
                    nbub = nmom_IP(1)
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb*nmom
                        q_cons_vf(eqn_idx%bub%beg + q - 1)%sf(j, k, l) = nbub*nmom_IP(q)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        q_cons_vf(eqn_idx%bub%beg + (q - 1)*nmom)%sf(j, k, l) = nbub
                    end do

                    if (.not. polytropic) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = 1, nnode
                                pb_in(j, k, l, r, q) = presb_IP((q - 1)*nnode + r)
                                mv_in(j, k, l, r, q) = massv_IP((q - 1)*nnode + r)
                            end do
                        end do
                    end if
                end if

                if (model_eqns == model_eqns_6eq) then
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = eqn_idx%int_en%beg, eqn_idx%int_en%end
                        alpha_q = alpha_IP(q - eqn_idx%int_en%beg + 1)
                        alpha_rho_q = alpha_rho_GP(q - eqn_idx%int_en%beg + 1)
                        call s_phase_internal_energy(pres_GP, alpha_q, alpha_rho_q, q - eqn_idx%int_en%beg + 1, e_q)
                        q_cons_vf(q)%sf(j, k, l) = e_q
                    end do
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! Publish the staged state now that every stencil has been read, so the buried-point average below sees it.
            $:GPU_PARALLEL_LOOP(private='[i, q, j, k, l, gp]', present='[ghost_points, gp_prim_new]')
            do i = 1, num_gps
                gp = ghost_points(i)
                j = gp%loc(1); k = gp%loc(2); l = gp%loc(3)
                do q = 1, sys_size
                    q_prim_vf(q)%sf(j, k, l) = gp_prim_new(q, i)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! A ghost point whose image point is buried in a neighboring IB has no fluid to mirror, so
            ! it takes the average of the ghost points corrected above, growing the stencil until it
            ! reaches one. Those neighbors already carry the wall condition, so the average does too.
            $:GPU_PARALLEL_LOOP(private='[i, j, k, l, q, r, jj, kk, ll, gp, rad, rad_z, num_nbrs, buf, buf_prim]', &
                                & present='[ghost_points]')
            do i = 1, num_gps
                gp = ghost_points(i)
                if (gp%interp_valid) cycle
                j = gp%loc(1)
                k = gp%loc(2)
                l = gp%loc(3)

                num_nbrs = 0
                rad = 0
                rad_z = 0
                do while (num_nbrs == 0 .and. rad < gp_layers)
                    rad = rad + 1
                    if (p /= 0) rad_z = rad
                    do jj = j - rad, j + rad
                        do kk = k - rad, k + rad
                            do ll = l - rad_z, l + rad_z
                                num_nbrs = num_nbrs + corrected_gps%sf(jj, kk, ll)
                            end do
                        end do
                    end do
                end do

                ! Fully enclosed: nothing to average from, so leave the interior placeholder in place
                if (num_nbrs == 0) cycle

                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, sys_size
                    buf = 0._wp
                    buf_prim = 0._wp
                    do jj = j - rad, j + rad
                        do kk = k - rad, k + rad
                            do ll = l - rad_z, l + rad_z
                                if (corrected_gps%sf(jj, kk, ll) == 1) then
                                    buf = buf + q_cons_vf(q)%sf(jj, kk, ll)
                                    buf_prim = buf_prim + q_prim_vf(q)%sf(jj, kk, ll)
                                end if
                            end do
                        end do
                    end do
                    q_cons_vf(q)%sf(j, k, l) = buf/real(num_nbrs, wp)
                    q_prim_vf(q)%sf(j, k, l) = buf_prim/real(num_nbrs, wp)
                end do

                if (qbmm .and. .not. polytropic) then
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        $:GPU_LOOP(parallelism='[seq]')
                        do r = 1, nnode
                            buf = 0._wp
                            buf_prim = 0._wp
                            do jj = j - rad, j + rad
                                do kk = k - rad, k + rad
                                    do ll = l - rad_z, l + rad_z
                                        if (corrected_gps%sf(jj, kk, ll) == 1) then
                                            buf = buf + pb_in(jj, kk, ll, r, q)
                                            buf_prim = buf_prim + mv_in(jj, kk, ll, r, q)
                                        end if
                                    end do
                                end do
                            end do
                            pb_in(j, k, l, r, q) = buf/real(num_nbrs, wp)
                            mv_in(j, k, l, r, q) = buf_prim/real(num_nbrs, wp)
                        end do
                    end do
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        n_surface_not_converged = n_surface_not_converged + n_not_converged
        n_surface_ill_posed = n_surface_ill_posed + n_ill_posed
        n_fill_tier = n_fill_tier + [nf1, nf2, nf3]

    end subroutine s_ibm_correct_state

    !> One line at the end of a run if the reacting surface ever failed to solve. Silence means every ghost point on every step
    !! reached the surface-balance tolerance, i.e. the surface chemistry was actually applied everywhere it was asked for.
    impure subroutine s_report_ibm_surface

        integer               :: n_not_converged_glb, n_ill_posed_glb, q
        integer, dimension(3) :: ft  !< Global fill-tier tallies

        call s_mpi_reduce_int_sum(n_surface_not_converged, n_not_converged_glb)
        call s_mpi_reduce_int_sum(n_surface_ill_posed, n_ill_posed_glb)
        ! Every rank reduces before the root-only return below: a collective after it would deadlock rank 0
        do q = 1, 3
            call s_mpi_reduce_int_sum(n_fill_tier(q), ft(q))
        end do

        if (proc_rank /= 0) return

        if (n_not_converged_glb > 0) then
            print '(A,I0,A)', ' Immersed-boundary surface chemistry: the Newton solve did not converge at ', n_not_converged_glb, &
                & ' ghost-point updates, which fell back to a chemically inert wall (Twall kept for thermal_bc = 1 or 3, ' &
                & // 'zero temperature gradient for thermal_bc = 2).'
        end if

        if (n_ill_posed_glb > 0) then
            print '(A,I0,A)', ' Immersed-boundary surface chemistry: ', n_ill_posed_glb, &
                & ' ghost-point updates had a zero levelset distance, leaving the surface balance undefined; those points ' &
                & // 'fell back to a chemically inert wall (Twall kept for thermal_bc = 1 or 3, zero temperature gradient for ' &
                & // 'thermal_bc = 2). Check the immersed geometry against the grid.'
        end if

        if (ib_wall_quad .and. sum(ft) > 0) then
            print '(A,3(F6.1,A))', ' [WALLFILL] ghost fill tiers: quadratic ', 1.e2_wp*ft(1)/sum(ft), '%, linear ', &
                & 1.e2_wp*ft(2)/sum(ft), '%, default ', 1.e2_wp*ft(3)/sum(ft), '%'
        end if

    end subroutine s_report_ibm_surface

    !> Store ghost point i's surface record in gp_surf. The weight is the surface area the point stands for. Ghost points within two
    !! cell sizes h = dV^(1/d) of the surface fill a band of volume ~2hA, so dV/(2h) each sums to the area A (a length in 2D)
    !! whatever the surface's orientation to the grid. The band lies inside the solid, where a layer at depth s has area A (1 -
    !! s/R)^(d-1) on a curved surface; for circles, spheres and cylinders the weight is scaled back by (R/(R - s))^(d-1), which
    !! removes that bias (sampled on random lattice offsets: < 0.1% mean, 0.3% spread for a sphere at R = 12h). The factor is
    !! limited to R/(R - s) <= 2, so a body of only a few cells cannot give an infinite or negative weight.
    subroutine s_record_gp_surface(i, gp, T_s, mdot_s, reacting)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)           :: i
        type(ghost_point), intent(in) :: gp
        real(wp), intent(in)          :: T_s, mdot_s
        logical, intent(in)           :: reacting
        real(wp)                      :: d, dV, h, R, curv, axial
        real(wp), dimension(3)        :: axis

        d = abs(real(gp%levelset, kind=wp))
        dV = dx(gp%loc(1))*dy(gp%loc(2))
        if (num_dims == 3) dV = dV*dz(gp%loc(3))
        h = dV**(1._wp/real(num_dims, wp))
        if (.not. (d > 0._wp .and. d <= 2._wp*h)) return

        gp_surf(1, i) = dV/(2._wp*h)
        R = patch_ib(gp%ib_patch_id)%radius
        curv = 1._wp
        if (R > 0._wp) curv = R/max(R - d, 0.5_wp*R)
        select case (patch_ib(gp%ib_patch_id)%geometry)
        case (2, 8)
            gp_surf(1, i) = gp_surf(1, i)*curv**(num_dims - 1)
        case (10)
            ! A cylinder's axis is its body-frame x, which rotation_matrix maps to its first column in the lab frame. A normal
            ! along the axis is a flat cap, which needs no correction.
            axis(:) = patch_ib(gp%ib_patch_id)%rotation_matrix(:,1)
            axial = abs(dot_product(axis, gp%levelset_norm))
            if (axial < 0.5_wp) gp_surf(1, i) = gp_surf(1, i)*curv
        end select

        ! ib_flux_wall: a reacting body's records come only from the ghosts that set its face fluxes, with their faces' wall area
        if (ib_flux_wall) then
            if (patch_ib(gp%ib_patch_id)%surface_reaction == 1) gp_surf(1, i) = f_flux_wall_record_area(ib_wall_src, &
                & ib_wall_theta, ib_markers, gp%loc)
        end if

        gp_surf(2, i) = T_s
        gp_surf(3, i) = 0._wp
        if (reacting) gp_surf(3, i) = mdot_s

    end subroutine s_record_gp_surface

    !> Write this rank's surface records to D/ib_surface_<rank>_<save>.dat, one line per weighted ghost point.
    impure subroutine s_write_ib_surface(save_count)

        integer, intent(in)                  :: save_count
        character(LEN=path_len + 2*name_len) :: file_loc
        real(wp)                             :: x(3)
        integer                              :: i, g, unit
        integer, allocatable                 :: order(:)

        if (.not. allocated(gp_surf)) return

        if (num_gps > 0) then
            $:GPU_UPDATE(host='[gp_surf(:, 1:num_gps), ghost_points(1:num_gps)]')
        end if

        call s_sort_gps_by_location(order)

        write (file_loc, '(A,I0,A,I0,A)') trim(case_dir) // '/D/ib_surface_', proc_rank, '_', save_count, '.dat'
        open (newunit=unit, file=trim(file_loc), status='replace', action='write')
        do i = 1, num_gps
            g = order(i)
            if (.not. gp_surf(1, g) > 0._wp) cycle
            x = 0._wp
            x(1) = x_cc(ghost_points(g)%loc(1))
            x(2) = y_cc(ghost_points(g)%loc(2))
            if (num_dims == 3) x(3) = z_cc(ghost_points(g)%loc(3))
            write (unit, '(3ES16.8,1X,I0,3ES15.6,3ES16.8)') x, patch_ib(ghost_points(g)%ib_patch_id)%gbl_patch_id, &
                   & ghost_points(g)%levelset_norm, gp_surf(:,g)
        end do
        close (unit)

    end subroutine s_write_ib_surface

    !> The surface temperature is held at Twall: prescribed (thermal_bc = 1) or evolved by the heat balance (thermal_bc = 3).
    pure logical function f_holds_twall(thermal_bc)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: thermal_bc

        f_holds_twall = thermal_bc == 1 .or. thermal_bc == 3

    end function f_holds_twall

    !> Solid volume of a circle (per unit depth), sphere or cylinder IB.
    pure real(wp) function f_ib_volume(ib_patch)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(ib_patch_parameters), intent(in) :: ib_patch

        select case (ib_patch%geometry)
        case (2)
            f_ib_volume = pi*ib_patch%radius**2
        case (8)
            f_ib_volume = 4._wp/3._wp*pi*ib_patch%radius**3
        case default
            f_ib_volume = pi*ib_patch%radius**2*max(ib_patch%length_x, ib_patch%length_y, ib_patch%length_z)
        end select

    end function f_ib_volume

    !> Advance the temperature of each thermal_bc = 3 IB over one step, treating it as one lumped body (small Biot number). With
    !! E_solid = m c_s (T - T_ref), the step's face energy in from s_accumulate_ib_face_fluxes and the gasification rate from the
    !! surface solve, m c_s dT/dt = Q_in + c_s (T - T_ref) Mdot_out + heat_power - emissivity sigma A (T^4 - T_rad^4). Mdot_out is
    !! the solid's mass loss: the face mass sum also carries the ghost cells' non-conservation (a net inflow on a reacting
    !! cylinder). T_ref = 298.15 K is the reference of the formation enthalpies, so the gas energy the solid's lost carbon carries
    !! is consistent. The step is explicit: the body's thermal time m c_s/(h A) is many orders of magnitude longer than a flow step.
    !! For D/ib_heat.dat every body also keeps its gasification rate, face mass out and fluid heat in, Q_in (+ the c_s term for
    !! thermal_bc = 3), in patch_ib%mdot, %mdot_face and %qdot.
    impure subroutine s_update_ib_temperatures(dt_step)

        real(wp), intent(in) :: dt_step
        real(wp), parameter  :: sigma_sb = 5.670374419e-8_wp  !< Stefan-Boltzmann constant [W/m^2/K^4]
        real(wp), parameter  :: T_ref = 298.15_wp
        real(wp)             :: q_in
        integer              :: i, g

        ! Area and gasification rate from the surface records of the last stage

        $:GPU_PARALLEL_LOOP(private='[g]')
        do g = 1, num_gbl_ibs
            ib_heat(1, g) = 0._wp
            ib_heat(4, g) = 0._wp
        end do
        $:END_GPU_PARALLEL_LOOP()
        if (num_gps > 0) then
            $:GPU_PARALLEL_LOOP(private='[i, g]')
            do i = 1, num_gps
                if (gp_surf(1, i) > 0._wp) then
                    g = patch_ib(ghost_points(i)%ib_patch_id)%gbl_patch_id
                    $:GPU_ATOMIC(atomic='update')
                    ib_heat(1, g) = ib_heat(1, g) + gp_surf(1, i)
                    $:GPU_ATOMIC(atomic='update')
                    ib_heat(4, g) = ib_heat(4, g) + gp_surf(1, i)*gp_surf(3, i)
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        $:GPU_UPDATE(host='[ib_heat]')
        call s_mpi_allreduce_vectors_sum(ib_heat, ib_heat, 4, num_gbl_ibs)
        $:GPU_UPDATE(device='[ib_heat]')

        ! On the device, where patch_ib is current (moving IBs update it there). Twall advances only for the bodies this rank holds,
        ! which stays consistent because thermal_bc = 3 requires moving_ibm = 0 (case_validator)
        $:GPU_PARALLEL_LOOP(private='[i, g, q_in]', copyin='[dt_step]')
        do i = 1, num_ibs
            g = patch_ib(i)%gbl_patch_id
            q_in = ib_heat(2, g)
            if (patch_ib(i)%thermal_bc == 3) q_in = q_in + patch_ib(i)%cp_solid*(patch_ib(i)%Twall - T_ref)*ib_heat(4, g)
            patch_ib(i)%mdot = ib_heat(4, g)
            patch_ib(i)%mdot_face = ib_heat(3, g)
            patch_ib(i)%qdot = q_in
            if (patch_ib(i)%thermal_bc == 3) then
                q_in = q_in + patch_ib(i)%heat_power - patch_ib(i)%emissivity*sigma_sb*ib_heat(1, &
                                       & g)*(patch_ib(i)%Twall**4 - patch_ib(i)%T_rad**4)
                patch_ib(i)%Twall = min(max(patch_ib(i)%Twall &
                         & + dt_step*q_in/(patch_ib(i)%rho_solid*patch_ib(i)%cp_solid*f_ib_volume(patch_ib(i))), T_surface_min), &
                         & T_surface_max)
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! The next step accumulates afresh
        $:GPU_PARALLEL_LOOP(private='[g]')
        do g = 1, num_gbl_ibs
            ib_heat(2, g) = 0._wp
            ib_heat(3, g) = 0._wp
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_update_ib_temperatures

    !> Ghost-point indices ordered by grid location, so per-point output does not inherit the discovery order. s_find_ghost_points
    !! assigns slots from a GPU atomic, so that order is the thread schedule's and differs between CPU and GPU builds
    !! (MFlowCode/MFC#2039). Heapsort, since a GPU run's input is unordered and a particle array carries O(10^5) ghost points.
    pure subroutine s_sort_gps_by_location(order)

        integer, allocatable, intent(out) :: order(:)
        integer(8), allocatable           :: key(:)
        integer                           :: i, root, g
        integer(8)                        :: k

        allocate (order(max(num_gps, 1)), key(max(num_gps, 1)))
        do i = 1, num_gps
            order(i) = i
            ! Lexicographic (i, j, k) as one key; loc is a buffer-inclusive index, so offset it to stay non-negative.
            k = int(ghost_points(i)%loc(1) + buff_size, 8)
            k = k*int(n + 1 + 2*buff_size, 8) + int(ghost_points(i)%loc(2) + buff_size, 8)
            key(i) = k*int(p + 1 + 2*buff_size, 8) + int(ghost_points(i)%loc(3) + buff_size, 8)
        end do

        do root = num_gps/2, 1, -1
            call s_sift_gp(order, key, root, num_gps)
        end do
        do i = num_gps, 2, -1
            g = order(1); order(1) = order(i); order(i) = g
            call s_sift_gp(order, key, 1, i - 1)
        end do

    end subroutine s_sort_gps_by_location

    !> Restore the max-heap rooted at root over order(1:last), comparing by key.
    pure subroutine s_sift_gp(order, key, root, last)

        integer, intent(inout) :: order(:)
        integer(8), intent(in) :: key(:)
        integer, intent(in)    :: root, last
        integer                :: r, child, g

        r = root
        do
            child = 2*r
            if (child > last) exit
            if (child < last) then
                if (key(order(child + 1)) > key(order(child))) child = child + 1
            end if
            if (key(order(child)) <= key(order(r))) exit
            g = order(r); order(r) = order(child); order(child) = g
            r = child
        end do

    end subroutine s_sift_gp

    !> Compute the image points for each ghost point
    impure subroutine s_compute_image_points()

        real(wp)                        :: dist
        real(wp)                        :: min_cell_width
        real(wp)                        :: temp_loc
        real(wp), pointer, dimension(:) :: s_cc => null()
        integer                         :: bound
        type(ghost_point)               :: gp
        integer                         :: q, dim      !< Iterator variables
        integer                         :: i, j, k, l  !< Location indexes
        integer                         :: patch_id    !< IB Patch ID
        integer                         :: dir
        integer                         :: index
        logical                         :: bounds_error

        bounds_error = .false.

        $:GPU_PARALLEL_LOOP(private='[q, gp, i, j, k, patch_id, dist, min_cell_width, dim, bound, dir, index, temp_loc, s_cc]', &
                            & copy='[bounds_error]', present='[ghost_points]')
        do q = 1, num_gps
            gp = ghost_points(q)
            i = gp%loc(1)
            j = gp%loc(2)
            k = gp%loc(3)

            ! Calculate and store the precise location of the image point
            patch_id = gp%ib_patch_id
            dist = abs(real(gp%levelset, kind=wp))
            ghost_points(q)%ip_dist = dist
            ! keep the image point stencil off the wall
            if (ib_second_order_vel) then
                min_cell_width = min(dx(i), dy(j))
                if (p > 0) min_cell_width = min(min_cell_width, dz(k))
                ghost_points(q)%ip_dist = max(dist, ib_ip_min_dist*min_cell_width)
            end if
            ! Per-component scalar stores: the array-valued form of this assignment was intermittently overwritten by another
            ! thread on AMD GPUs (a few ip_loc(2) of ~4000 per call), which sent the image-point search past the halo
            ghost_points(q)%ip_loc(1) = x_cc(i) + (dist + ghost_points(q)%ip_dist)*gp%levelset_norm(1)
            ghost_points(q)%ip_loc(2) = y_cc(j) + (dist + ghost_points(q)%ip_dist)*gp%levelset_norm(2)
            ghost_points(q)%ip_loc(3) = 0._wp
            if (p > 0) ghost_points(q)%ip_loc(3) = z_cc(k) + (dist + ghost_points(q)%ip_dist)*gp%levelset_norm(3)

            ! Find the closest grid point to the image point
            do dim = 1, num_dims
                ! s_cc points to the dim array we need
                if (dim == 1) then
                    s_cc => x_cc
                    bound = m + buff_size - 1
                else if (dim == 2) then
                    s_cc => y_cc
                    bound = n + buff_size - 1
                else
                    s_cc => z_cc
                    bound = p + buff_size - 1
                end if

                if (f_approx_equal(gp%levelset_norm(dim), 0._wp)) then
                    ! if the ghost point is almost equal to a cell location, we set it equal and continue
                    ghost_points(q)%ip_grid(dim) = ghost_points(q)%loc(dim)
                else
                    if (gp%levelset_norm(dim) > 0) then
                        dir = 1
                    else
                        dir = -1
                    end if

                    index = ghost_points(q)%loc(dim)
                    temp_loc = ghost_points(q)%ip_loc(dim)
                    do while ((temp_loc < s_cc(index) .or. temp_loc > s_cc(index + 1)) .and. (.not. bounds_error))
                        index = index + dir
                        if (index < -buff_size .or. index > bound) then
#if !defined(MFC_OpenACC) && !defined(MFC_OpenMP)
                            print *, "A required image point is not located in this computational domain."
                            print *, "Ghost Point is located at :"
                            if (p == 0) then
                                print *, [x_cc(i), y_cc(j)]
                            else
                                print *, [x_cc(i), y_cc(j), z_cc(k)]
                            end if
                            print *, "We are searching in dimension ", dim, " for image point at ", ghost_points(q)%ip_loc(:)
                            print *, "Domain size: "
                            print *, "x: ", x_cc(-buff_size), " to: ", x_cc(m + buff_size - 1)
                            print *, "y: ", y_cc(-buff_size), " to: ", y_cc(n + buff_size - 1)
                            if (p /= 0) print *, "z: ", z_cc(-buff_size), " to: ", z_cc(p + buff_size - 1)
                            print *, "Image point is located approximately ", &
                                & (ghost_points(q)%loc(dim) - ghost_points(q) %ip_loc(dim))/(s_cc(1) - s_cc(0)), " grid cells away"
                            print *, "Levelset ", dist, " and Norm: ", gp%levelset_norm
                            print *, &
                                & "A short term fix may include increasing buff_size further in m_helper_basic (currently set to a minimum of 10)"
#endif
                            bounds_error = .true.
                        end if
                    end do

                    ghost_points(q)%ip_grid(dim) = index
                    if (ghost_points(q)%DB(dim) == -1) then
                        ghost_points(q)%ip_grid(dim) = ghost_points(q)%loc(dim) + 1
                    else if (ghost_points(q)%DB(dim) == 1) then
                        ghost_points(q)%ip_grid(dim) = ghost_points(q)%loc(dim) - 1
                    end if
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        @:PROHIBIT(bounds_error, "Ghost Point and Image Point on Different Processors. Exiting")

    end subroutine s_compute_image_points

    !> Count the number of ghost points for memory allocation
    subroutine s_find_num_ghost_points(num_gps_out)

        integer, intent(out) :: num_gps_out
        integer              :: i, j, k, ii, jj, kk, gp_layers_z  !< Iterator variables
        integer              :: num_gps_local                     !< local copies of the gp count to support GPU compute
        logical              :: is_gp

        num_gps_local = 0
        gp_layers_z = gp_layers
        if (p == 0) gp_layers_z = 0

        $:GPU_PARALLEL_LOOP(private='[i, j, k, ii, jj, kk, is_gp]', copy='[num_gps_local]', firstprivate='[gp_layers, &
                            & gp_layers_z]', collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (ib_markers%sf(i, j, k) /= 0) then
                        is_gp = .false.
                        marker_search: do ii = i - gp_layers, i + gp_layers
                            do jj = j - gp_layers, j + gp_layers
                                do kk = k - gp_layers_z, k + gp_layers_z
                                    if (ib_markers%sf(ii, jj, kk) == 0) then
                                        ! if any neighbors are not in the IB, it is a ghost point
                                        is_gp = .true.
                                        exit marker_search
                                    end if
                                end do
                            end do
                        end do marker_search

                        if (is_gp) then
                            $:GPU_ATOMIC(atomic='update')
                            num_gps_local = num_gps_local + 1
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        num_gps_out = num_gps_local

    end subroutine s_find_num_ghost_points

    !> Locate all ghost points in the domain
    subroutine s_find_ghost_points()

        integer :: i, j, k, ii, jj, kk, gp_layers_z  !< Iterator variables
        integer :: xp, yp, zp                        !< periodicities
        integer :: count, count_i, local_idx
        integer :: patch_id, encoded_patch_id, neighborhood_patch_id
        logical :: is_gp

        count = 0
        count_i = 0
        gp_layers_z = gp_layers
        if (p == 0) gp_layers_z = 0

        $:GPU_PARALLEL_LOOP(private='[i, j, k, ii, jj, kk, is_gp, local_idx, patch_id, encoded_patch_id, neighborhood_patch_id, &
                            & xp, yp, zp]', copyin='[count, count_i, glb_bounds]', firstprivate='[gp_layers, gp_layers_z]', &
                            & present='[ghost_points]', collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (ib_markers%sf(i, j, k) /= 0) then
                        is_gp = .false.
                        marker_search: do ii = i - gp_layers, i + gp_layers
                            do jj = j - gp_layers, j + gp_layers
                                do kk = k - gp_layers_z, k + gp_layers_z
                                    if (ib_markers%sf(ii, jj, kk) == 0) then
                                        ! if any neighbors are not in the IB, it is a ghost point
                                        is_gp = .true.
                                        exit marker_search
                                    end if
                                end do
                            end do
                        end do marker_search

                        if (is_gp) then
                            $:GPU_ATOMIC(atomic='capture')
                            count = count + 1
                            local_idx = count
                            $:END_GPU_ATOMIC_CAPTURE()
                            ghost_points(local_idx)%loc(1) = i
                            ghost_points(local_idx)%loc(2) = j
                            ghost_points(local_idx)%loc(3) = k
                            encoded_patch_id = ib_markers%sf(i, j, k)
                            call s_decode_patch_periodicity(encoded_patch_id, patch_id, xp, yp, zp)
                            call s_get_neighborhood_idx(patch_id, neighborhood_patch_id)
                            ghost_points(local_idx)%ib_patch_id = neighborhood_patch_id
                            ghost_points(local_idx)%x_periodicity = xp
                            ghost_points(local_idx)%y_periodicity = yp
                            ghost_points(local_idx)%z_periodicity = zp
                            ghost_points(local_idx)%slip = patch_ib(neighborhood_patch_id)%slip

                            #:for X, ID, IDX in [('x', 1, 'i'), ('y', 2, 'j'), ('z', 3, 'k')]
                                ghost_points(local_idx)%DB(${ID}$) = 0
                                if (${ID}$ <= num_dims) then  ! Onlyrun the numeher of dimensions present
                                    if (ib_bc_${X}$%beg /= BC_PERIODIC) then
                                        if ((${X}$_cc(${IDX}$) - d${X}$(${IDX}$)) < glb_bounds(${ID}$)%beg) then
                                            ! if the grid cell is in a wall on the "left"
                                            ghost_points(local_idx)%DB(${ID}$) = -1
                                        else if ((${X}$_cc(${IDX}$) + d${X}$(${IDX}$)) > glb_bounds(${ID}$)%end) then
                                            ! if the grid cell is in a wall on the "right"
                                            ghost_points(local_idx)%DB(${ID}$) = 1
                                        end if
                                    end if
                                end if
                            #:endfor
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_find_ghost_points

    !> Compute the interpolation coefficients for image points
    subroutine s_compute_interpolation_coeffs()

        real(wp), dimension(2, 2, 2) :: dist
        real(wp), dimension(2, 2, 2) :: fluid_mask
        real(wp), dimension(2, 2, 2) :: interp_coeffs
        real(wp), dimension(3, 3, 3) :: stencil_coeffs          !< Coefficients on the stencil that starts at ip_stencil
        real(wp), dimension(3, 3)    :: wd                      !< 1D TSC weights per direction
        real(wp)                     :: fluid_w                 !< Fraction of the TSC weight carried by fluid cells
        integer, dimension(3)        :: ic                      !< Cell centring the TSC stencil
        type(ghost_point)            :: gp
        integer                      :: q, i, j, k, ii, jj, kk  !< Grid indexes and iterators

        $:GPU_PARALLEL_LOOP(private='[q, i, j, k, ii, jj, kk, dist, gp, interp_coeffs, fluid_mask, stencil_coeffs, wd, fluid_w, &
                            & ic]', present='[ghost_points]')
        do q = 1, num_gps
            gp = ghost_points(q)
            stencil_coeffs = 0._wp
            ic = 0
            if (ib_ip_interp == 2) then
                call s_tsc_weights(ib_markers, gp%ip_loc, ic, wd, fluid_w)
                stencil_coeffs = f_tsc_stencil(wd)
                ghost_points(q)%interp_valid = fluid_w > 0._wp
            else
                ! Get the interpolation points
                i = gp%ip_grid(1)
                j = gp%ip_grid(2)
                if (p /= 0) then
                    k = gp%ip_grid(3)
                else
                    k = 0
                end if

                ! get the distance to a cell in each direction
                dist = 0._wp
                do ii = 0, 1
                    do jj = 0, 1
                        if (p == 0) then
                            dist(1 + ii, 1 + jj, 1) = sqrt((x_cc(i + ii) - gp%ip_loc(1))**2 + (y_cc(j + jj) - gp%ip_loc(2))**2)
                        else
                            do kk = 0, 1
                                dist(1 + ii, 1 + jj, &
                                     & 1 + kk) = sqrt((x_cc(i + ii) - gp%ip_loc(1))**2 + (y_cc(j + jj) - gp%ip_loc(2))**2 &
                                     & + (z_cc(k + kk) - gp%ip_loc(3))**2)
                            end do
                        end if
                    end do
                end do

                ! flag the fluid cells of the stencil
                fluid_mask = 0._wp
                do kk = 0, min(p, 1)
                    do jj = 0, 1
                        do ii = 0, 1
                            if (ib_markers%sf(i + ii, j + jj, k + kk) == 0) fluid_mask(ii + 1, jj + 1, kk + 1) = 1._wp
                        end do
                    end do
                end do

                call s_interpolation_weights(gp, i, j, k, dist, fluid_mask, interp_coeffs)

                stencil_coeffs(1:2,1:2,1:2) = interp_coeffs
                ghost_points(q)%interp_valid = any(interp_coeffs > 0._wp)
            end if

            ! An image point buried in a neighboring IB gets no weights at all: its ghost point is
            ! averaged from its own neighbors instead, in s_ibm_correct_state.
            do kk = 1, 3
                do jj = 1, 3
                    do ii = 1, 3
                        ghost_points(q)%interp_coeffs(ii, jj, kk) = stencil_coeffs(ii, jj, kk)
                    end do
                end do
                ! The stencil's lower corner: one cell below the TSC centre, or the image point's grid cell
                ghost_points(q)%ip_stencil(kk) = merge(ic(kk) - 1, gp%ip_grid(kk), ib_ip_interp == 2)
            end do
            if (p == 0) ghost_points(q)%ip_stencil(3) = 0
            if (ghost_points(q)%interp_valid) corrected_gps%sf(gp%loc(1), gp%loc(2), gp%loc(3)) = 1
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_interpolation_coeffs

    !> Triangular-shaped-cloud (quadratic B-spline) weights of a point on the 3x3(x3) stencil centred on the cell nearest to it
    !! (Ezra, Peles & Kozak 2025, Int. J. Heat Mass Transf. 252:127314), and the fraction of that weight carried by fluid cells. The
    !! kernel is a partition of unity and smooth across cell boundaries, so an interpolated wall quantity converges one order faster
    !! than from inverse-distance or bilinear weights, which are only continuous. The weights use the fractional position in the
    !! bracketing cell at integer offsets, which keeps the partition of unity exactly on a stretched grid. The stencil is kept
    !! inside the halo. Solid cells in it carry the ghost values of the last correction (see ib_keep_ghost); a point buried in a
    !! body has no fluid weight.
    subroutine s_tsc_weights(markers, loc, ic, wd, fluid_frac)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(integer_field), intent(in) :: markers    !< ib_markers, passed so the device routine reads the caller's present copy
        real(wp), dimension(3), intent(in) :: loc
        integer, dimension(3), intent(out) :: ic      !< Cell centring the stencil
        real(wp), dimension(3, 3), intent(out) :: wd  !< 1D weights on ic - 1, ic, ic + 1, per direction
        real(wp), intent(out) :: fluid_frac
        real(wp) :: ds, frac, wsum, wfluid, buf
        integer :: dim, ii, jj, kk

        ic = 0
        wd = 0._wp
        do dim = 1, merge(3, 2, p /= 0)
            if (dim == 1) then
                ic(1) = max(min(floor((loc(1) - x_cc(0))/(x_cc(1) - x_cc(0))), m + buff_size - 2), -buff_size + 2)
                ! the estimate assumes a uniform grid: one step brackets the point on a mildly stretched one
                if (loc(1) < x_cc(ic(1))) ic(1) = max(ic(1) - 1, -buff_size + 2)
                if (loc(1) >= x_cc(ic(1) + 1)) ic(1) = min(ic(1) + 1, m + buff_size - 2)
                ds = x_cc(ic(1) + 1) - x_cc(ic(1)); frac = (loc(1) - x_cc(ic(1)))/ds
            else if (dim == 2) then
                ic(2) = max(min(floor((loc(2) - y_cc(0))/(y_cc(1) - y_cc(0))), n + buff_size - 2), -buff_size + 2)
                if (loc(2) < y_cc(ic(2))) ic(2) = max(ic(2) - 1, -buff_size + 2)
                if (loc(2) >= y_cc(ic(2) + 1)) ic(2) = min(ic(2) + 1, n + buff_size - 2)
                ds = y_cc(ic(2) + 1) - y_cc(ic(2)); frac = (loc(2) - y_cc(ic(2)))/ds
            else
                ic(3) = max(min(floor((loc(3) - z_cc(0))/(z_cc(1) - z_cc(0))), p + buff_size - 2), -buff_size + 2)
                if (loc(3) < z_cc(ic(3))) ic(3) = max(ic(3) - 1, -buff_size + 2)
                if (loc(3) >= z_cc(ic(3) + 1)) ic(3) = min(ic(3) + 1, p + buff_size - 2)
                ds = z_cc(ic(3) + 1) - z_cc(ic(3)); frac = (loc(3) - z_cc(ic(3)))/ds
            end if
            if (frac > 0.5_wp) then
                ic(dim) = ic(dim) + 1
                frac = frac - 1._wp
            end if
            wd(1, dim) = 0.5_wp*(1.5_wp - abs(frac + 1._wp))**2
            wd(2, dim) = 0.75_wp - frac*frac
            wd(3, dim) = 0.5_wp*(1.5_wp - abs(frac - 1._wp))**2
        end do
        if (p == 0) wd(1, 3) = 1._wp  ! 2D: the consumers index this slot as k - k1 + 1 = 1

        wsum = 0._wp; wfluid = 0._wp
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, merge(3, 1, p /= 0)
                    buf = wd(ii, 1)*wd(jj, 2)*wd(kk, 3)
                    wsum = wsum + buf
                    if (markers%sf(ic(1) + ii - 2, ic(2) + jj - 2, ic(3) + (kk - 2)*min(p, 1)) == 0) wfluid = wfluid + buf
                end do
            end do
        end do
        fluid_frac = wfluid/wsum

    end subroutine s_tsc_weights

    !> Apply TSC weights from s_tsc_weights to one field.
    function f_tsc_apply(sf, ic, wd) result(val)

        $:GPU_ROUTINE(function_name='f_tsc_apply', parallelism='[seq]', cray_inline=True)
        type(scalar_field), intent(in)        :: sf
        integer, dimension(3), intent(in)     :: ic
        real(wp), dimension(3, 3), intent(in) :: wd
        real(wp)                              :: val
        integer                               :: ii, jj, kk

        val = 0._wp
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, merge(3, 1, p /= 0)
                    val = val + wd(ii, 1)*wd(jj, 2)*wd(kk, 3)*sf%sf(ic(1) + ii - 2, ic(2) + jj - 2, ic(3) + (kk - 2)*min(p, 1))
                end do
            end do
        end do

    end function f_tsc_apply

    !> Ghost value from the quadratic wall-normal profile, bounded by the linear mirror (minmod). The two extrapolations agree
    !! exactly for a linear profile, so this limits only the curvature correction, and to the size the linear mirror already
    !! produces. Unlimited, the quadratic amplifies that correction with ghost depth and can leave unphysical near-wall states.
    function f_quad_fill(phi_w, phi_1, phi_2, lw, llin) result(phi_g)

        $:GPU_ROUTINE(function_name='f_quad_fill', parallelism='[seq]', cray_inline=True)
        real(wp), intent(in)               :: phi_w, phi_1, phi_2
        real(wp), dimension(3), intent(in) :: lw
        real(wp), dimension(2), intent(in) :: llin
        real(wp)                           :: phi_g, d_quad, d_lin

        d_quad = lw(1)*phi_w + lw(2)*phi_1 + lw(3)*phi_2 - phi_w
        d_lin = llin(1)*phi_w + llin(2)*phi_1 - phi_w
        if (d_quad*d_lin <= 0._wp) then
            phi_g = phi_w + d_lin  ! opposite signs: the curvature term has reversed the extrapolation
        else if (abs(d_quad) <= abs(d_lin)) then
            phi_g = phi_w + d_quad
        else
            phi_g = phi_w + d_lin
        end if

    end function f_quad_fill

    !> Probes and Lagrange weights of the quadratic wall-normal fill of a ghost point. The wall is at s = 0, the probes at s1 and s2
    !! = 2 s1 along the surface normal from the wall foot, the ghost at s = -d (d = |levelset|):
    !!
    !!     phi_g = L0 phi_w + L1 phi_1 + L2 phi_2
    !!     L0 = (d + s1)(d + s2)/(s1 s2),  L1 = -d (d + s2)/(s1 (s2 - s1)),  L2 = d (d + s1)/(s2 (s2 - s1))
    !!
    !! L0 + L1 + L2 = 1 and at d = 0 they reduce to (1, 0, 0), so the wall value is exact on the wall. s1 = max(h, d) bounds the
    !! weights by (3, -3, 1) at any ghost depth. Tier 1 uses both probes; tier 2 is the linear mirror where the outer probe has too
    !! little fluid around it (a narrow gap, a rim, a neighbouring body); tier 3 leaves the default treatment.
    subroutine s_wall_fill(gp, wf)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(ghost_point), intent(in) :: gp
        type(wall_fill), intent(out)  :: wf
        real(wp), parameter           :: FLUID_FRAC_MIN = 0.5_wp  !< A probe needs this much fluid weight to be used
        real(wp), dimension(3)        :: norm, wall_loc
        real(wp)                      :: d, s1, s2, f1, f2, dV

        dV = dx(gp%loc(1))*dy(gp%loc(2))
        if (p > 0) dV = dV*dz(gp%loc(3))
        d = abs(real(gp%levelset, kind=wp))
        s1 = max(dV**(1._wp/real(num_dims, wp)), d)
        s2 = 2._wp*s1

        norm = gp%levelset_norm/max(sqrt(sum(gp%levelset_norm**2)), sgm_eps)
        wall_loc = 0._wp
        wall_loc(1) = x_cc(gp%loc(1)); wall_loc(2) = y_cc(gp%loc(2))
        if (p > 0) wall_loc(3) = z_cc(gp%loc(3))
        wall_loc = wall_loc + d*norm
        call s_tsc_weights(ib_markers, wall_loc + s1*norm, wf%ic1, wf%wd1, f1)
        call s_tsc_weights(ib_markers, wall_loc + s2*norm, wf%ic2, wf%wd2, f2)

        wf%llin = [1._wp + d/s1, -d/s1]
        wf%lw = [wf%llin(1), wf%llin(2), 0._wp]
        wf%tier = 3
        if (f1 < FLUID_FRAC_MIN) return
        wf%tier = 2
        if (f2 < FLUID_FRAC_MIN) return
        wf%lw = [(d + s1)*(d + s2)/(s1*s2), -d*(d + s2)/(s1*(s2 - s1)), d*(d + s1)/(s2*(s2 - s1))]
        wf%tier = 1

    end subroutine s_wall_fill

    !> Tensor product of the 1D TSC weights: the coefficient on each cell of the 3x3x3 stencil.
    function f_tsc_stencil(wd) result(c)

        $:GPU_ROUTINE(function_name='f_tsc_stencil', parallelism='[seq]', cray_inline=True)
        real(wp), dimension(3, 3), intent(in) :: wd
        real(wp), dimension(3, 3, 3)          :: c
        integer                               :: ii, jj, kk

        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, 3
                    c(ii, jj, kk) = wd(ii, 1)*wd(jj, 2)*wd(kk, 3)
                end do
            end do
        end do

    end function f_tsc_stencil

    !> Fluid-masked image-point weights: bilinear for second-order IBM, else inverse distance squared
    subroutine s_interpolation_weights(gp, i, j, k, dist, fluid_mask, interp_coeffs)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(in)             :: gp
        integer, intent(in)                       :: i, j, k      !< lower corner of the stencil
        real(wp), dimension(2, 2, 2), intent(in)  :: dist, fluid_mask
        real(wp), dimension(2, 2, 2), intent(out) :: interp_coeffs
        real(wp), dimension(2, 2, 2)              :: weights
        real(wp), dimension(2, 3)                 :: axis_weight  !< (stencil side, dimension) linear weights
        real(wp)                                  :: weight_sum
        integer                                   :: ii, jj, kk, num_z_layers

        num_z_layers = min(p, 1) + 1

        if (ib_second_order_vel) then
            ! bilinear weights, Mittal et al. (2008)
            axis_weight(2, 1) = (gp%ip_loc(1) - x_cc(i))/(x_cc(i + 1) - x_cc(i))
            axis_weight(2, 2) = (gp%ip_loc(2) - y_cc(j))/(y_cc(j + 1) - y_cc(j))
            axis_weight(2, 3) = 0._wp
            if (p /= 0) axis_weight(2, 3) = (gp%ip_loc(3) - z_cc(k))/(z_cc(k + 1) - z_cc(k))
            axis_weight(2,:) = min(max(axis_weight(2,:), 0._wp), 1._wp)  ! clamp to the cell edges
            axis_weight(1,:) = 1._wp - axis_weight(2,:)
            do kk = 1, num_z_layers
                do jj = 1, 2
                    do ii = 1, 2
                        weights(ii, jj, kk) = axis_weight(ii, 1)*axis_weight(jj, 2)*axis_weight(kk, 3)
                    end do
                end do
            end do
        else
            ! inverse distance squared, Chaudhuri et al. JCP (2011)
            if (minval(dist(:,:,1:num_z_layers)) <= 1.e-16_wp) then
                weights(:,:,1:num_z_layers) = merge(1._wp, 0._wp, dist(:,:,1:num_z_layers) <= 1.e-16_wp)
            else
                weights(:,:,1:num_z_layers) = 1._wp/dist(:,:,1:num_z_layers)**2
            end if
        end if

        interp_coeffs = 0._wp
        weights(:,:,1:num_z_layers) = fluid_mask(:,:,1:num_z_layers)*weights(:,:,1:num_z_layers)
        weight_sum = sum(weights(:,:,1:num_z_layers))
        if (weight_sum > 0._wp) interp_coeffs(:,:,1:num_z_layers) = weights(:,:,1:num_z_layers)/weight_sum

    end subroutine s_interpolation_weights

    !> Interpolate primitive variables to a ghost point's image point using bilinear or trilinear interpolation
    subroutine s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, pb_IP, mv_IP, &
                                         & nmom_IP, pb_in, mv_in, presb_IP, massv_IP, Ys_IP)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf  !< Primitive Variables
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_in, mv_in
        type(ghost_point), intent(in) :: gp
        real(wp), intent(inout) :: pres_IP
        real(wp), dimension(3), intent(inout) :: vel_IP
        real(wp), intent(inout) :: c_IP
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(inout) :: alpha_IP, alpha_rho_IP
        #:else
            real(wp), dimension(num_fluids), intent(inout) :: alpha_IP, alpha_rho_IP
        #:endif
        real(wp), optional, dimension(:), intent(inout) :: r_IP, v_IP, pb_IP, mv_IP
        real(wp), optional, dimension(:), intent(inout) :: nmom_IP
        real(wp), optional, dimension(:), intent(inout) :: presb_IP, massv_IP
        real(wp), optional, dimension(:), intent(inout) :: Ys_IP  !< Interpolated species mass fractions (chemistry)
        integer                                         :: i, j, k, l, q  !< Iterator variables
        integer                                         :: i1, i2, j1, j2, k1, k2  !< Iterator variables
        real(wp)                                        :: coeff

        i1 = gp%ip_stencil(1); i2 = i1 + ib_ip_width - 1
        j1 = gp%ip_stencil(2); j2 = j1 + ib_ip_width - 1
        k1 = gp%ip_stencil(3); k2 = k1 + ib_ip_width - 1

        if (p == 0) then
            k1 = 0
            k2 = 0
        end if

        alpha_rho_IP = 0._wp
        alpha_IP = 0._wp
        pres_IP = 0._wp
        vel_IP = 0._wp

        if (chemistry) Ys_IP = 0._wp

        if (surface_tension) c_IP = 0._wp

        if (bubbles_euler) then
            r_IP = 0._wp
            v_IP = 0._wp
            if (.not. polytropic) then
                mv_IP = 0._wp
                pb_IP = 0._wp
            end if
        end if

        if (qbmm) then
            nmom_IP = 0._wp
            if (.not. polytropic) then
                presb_IP = 0._wp
                massv_IP = 0._wp
            end if
        end if

        $:GPU_LOOP(parallelism='[seq]')
        do i = i1, i2
            $:GPU_LOOP(parallelism='[seq]')
            do j = j1, j2
                $:GPU_LOOP(parallelism='[seq]')
                do k = k1, k2
                    coeff = gp%interp_coeffs(i - i1 + 1, j - j1 + 1, k - k1 + 1)

                    pres_IP = pres_IP + coeff*q_prim_vf(eqn_idx%E)%sf(i, j, k)

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = eqn_idx%mom%beg, eqn_idx%mom%end
                        vel_IP(q + 1 - eqn_idx%mom%beg) = vel_IP(q + 1 - eqn_idx%mom%beg) + coeff*q_prim_vf(q)%sf(i, j, k)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do l = eqn_idx%cont%beg, eqn_idx%cont%end
                        alpha_rho_IP(l) = alpha_rho_IP(l) + coeff*q_prim_vf(l)%sf(i, j, k)
                        alpha_IP(l) = alpha_IP(l) + coeff*q_prim_vf(eqn_idx%adv%beg + l - 1)%sf(i, j, k)
                    end do

                    if (surface_tension) then
                        c_IP = c_IP + coeff*q_prim_vf(eqn_idx%c)%sf(i, j, k)
                    end if

                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_species
                            Ys_IP(q) = Ys_IP(q) + coeff*q_prim_vf(eqn_idx%species%beg + q - 1)%sf(i, j, k)
                        end do
                    end if

                    if (bubbles_euler .and. .not. qbmm) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do l = 1, nb
                            if (polytropic) then
                                r_IP(l) = r_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + (l - 1)*2)%sf(i, j, k)
                                v_IP(l) = v_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + 1 + (l - 1)*2)%sf(i, j, k)
                            else
                                r_IP(l) = r_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + (l - 1)*4)%sf(i, j, k)
                                v_IP(l) = v_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + 1 + (l - 1)*4)%sf(i, j, k)
                                pb_IP(l) = pb_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + 2 + (l - 1)*4)%sf(i, j, k)
                                mv_IP(l) = mv_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + 3 + (l - 1)*4)%sf(i, j, k)
                            end if
                        end do
                    end if

                    if (qbmm) then
                        do l = 1, nb*nmom
                            nmom_IP(l) = nmom_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg - 1 + l)%sf(i, j, k)
                        end do
                        if (.not. polytropic) then
                            do q = 1, nb
                                do l = 1, nnode
                                    presb_IP((q - 1)*nnode + l) = presb_IP((q - 1)*nnode + l) + coeff*real(pb_in(i, j, k, l, q), &
                                             & kind=wp)
                                    massv_IP((q - 1)*nnode + l) = massv_IP((q - 1)*nnode + l) + coeff*real(mv_in(i, j, k, l, q), &
                                             & kind=wp)
                                end do
                            end do
                        end if
                    end if
                end do
            end do
        end do

    end subroutine s_interpolate_image_point

    !> Resets the current indexes of immersed boundaries and replaces them after updating
    !> the position of each moving immersed boundary
    impure subroutine s_update_mib(num_ibs)

        integer, intent(in) :: num_ibs
        integer             :: i, j, k, z_buff_size

        call nvtxStartRange("UPDATE-MIBM")

        ! Clears the existing immersed boundary indices over the same halo extent s_get_bounding_indices draws them
        z_buff_size = 0; if (p /= 0) z_buff_size = buff_size
        $:GPU_PARALLEL_LOOP(private='[i, j, k]')
        do i = -buff_size, m + buff_size; do j = -buff_size, n + buff_size; do k = -z_buff_size, p + z_buff_size
            ib_markers%sf(i, j, k) = 0._wp
            corrected_gps%sf(i, j, k) = 0
        end do; end do; end do
        $:END_GPU_PARALLEL_LOOP()

        ! recalulcate the rotation matrix based upon the new angles
        $:GPU_PARALLEL_LOOP(private='[i]')
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) then
                call s_update_ib_rotation_matrix(i)
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! recompute the new ib_patch locations and broadcast them.
        call nvtxStartRange("APPLY-IB-PATCHES")
        call s_apply_ib_patches(ib_markers)
        call nvtxEndRange

        call nvtxStartRange("COMPUTE-GHOST-POINTS")
        ! recalculate the ghost point locations and coefficients
        call s_find_num_ghost_points(num_gps)
        ! overrunning ghost_points is a silent device fault; name the numbers instead
        if (int(num_gps, 8) > int(size(ghost_points), 8)) then
            print '(A,I0,A,I0,A,I0)', 'ghost points on rank ', proc_rank, ' grew to ', num_gps, ' but the array holds ', &
                & size(ghost_points)
            call s_mpi_abort('moving immersed boundary: ghost-point count exceeded its allocation')
        end if
        ! num_gps is a declare-target module variable and bounds every device loop over the ghost points; it is
        ! copied to the device once in s_ibm_setup, so without this refresh the kernels keep the setup-time count
        ! as a moving body's count changes: stale list entries beyond the current count are written as wall
        ! states into cells that are now fluid, new ghost cells are left uncorrected, and which entries those are
        ! depends on the atomic fill order -- nondeterministic results run to run (#1886).
        $:GPU_UPDATE(device='[num_gps]')
        call s_find_ghost_points()
        call nvtxEndRange

        call nvtxStartRange("COMPUTE-IMAGE-POINTS")
        call s_apply_levelset(ghost_points, num_gps)
        call s_compute_image_points()
        call s_compute_interpolation_coeffs()
        call nvtxEndRange

        call nvtxEndRange

    end subroutine s_update_mib

    !> log(cosh(x)) without overflow for large |x|
    pure function f_log_cosh(x) result(y)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: x
        real(wp)             :: y
        y = abs(x) + log(0.5_wp*(1._wp + exp(-2._wp*abs(x))))

    end function f_log_cosh

    !> Prescribed kinematics for patch i at time t: hinged flapping (kin_model = 1) or the Eldredge pitch ramp (kin_model = 2).
    !! Flapping rolls phi about the lab x axis through the hinge and pitches theta about the body spanwise (y) axis through it; the
    !! ramp holds phi at zero and drives theta alone. Both compose as R = Rx(phi) Ry(theta) (MFC's angle convention) and set the
    !! angles, centroid, velocity and lab-frame angular velocity; nothing is integrated, so restarts are exact.
    subroutine s_prescribed_kinematics(i, t)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)    :: i
        real(wp), intent(in)   :: t
        real(wp)               :: tau, amp, damp, omega, phi, theta, phid, thetad, arg_phi, arg_th, t_p, a_s
        real(wp), dimension(3) :: r, c, w, v

        tau = t - patch_ib(i)%kin_t0
        omega = 2._wp*pi*patch_ib(i)%kin_freq
        amp = 0._wp
        damp = 0._wp
        if (tau >= 0._wp) then
            amp = 1._wp
            if (patch_ib(i)%kin_ramp > 0._wp .and. tau < patch_ib(i)%kin_ramp) then
                amp = 0.5_wp*(1._wp - cos(pi*tau/patch_ib(i)%kin_ramp))
                damp = 0.5_wp*pi/patch_ib(i)%kin_ramp*sin(pi*tau/patch_ib(i)%kin_ramp)
            end if
        end if
        arg_phi = omega*tau
        arg_th = omega*tau + patch_ib(i)%kin_phase

        if (patch_ib(i)%kin_model == 2) then
            ! Eldredge smoothed linear pitch ramp and hold (AIAA canonical pitch-up): with a = kin_smooth, nominal rate
            ! Omega = kin_pitch_rate and pitch time t_p = theta0/Omega, G = log[cosh(a tau)/cosh(a (tau - t_p))] rises from
            ! -a t_p to +a t_p, so theta = mean + theta0/2 [1 + G/(a t_p)] and theta' = Omega/2 [tanh(a tau) - tanh(a (tau - t_p))]
            phi = 0._wp
            phid = 0._wp
            t_p = patch_ib(i)%kin_theta0/patch_ib(i)%kin_pitch_rate
            a_s = patch_ib(i)%kin_smooth
            theta = patch_ib(i)%kin_theta_mean + 0.5_wp*patch_ib(i)%kin_theta0*(1._wp + (f_log_cosh(a_s*tau) &
                             & - f_log_cosh(a_s*(tau - t_p)))/(a_s*t_p))
            thetad = 0.5_wp*patch_ib(i)%kin_pitch_rate*(tanh(a_s*tau) - tanh(a_s*(tau - t_p)))
        else
            phi = amp*patch_ib(i)%kin_phi0*sin(arg_phi)
            phid = damp*patch_ib(i)%kin_phi0*sin(arg_phi) + amp*patch_ib(i)%kin_phi0*omega*cos(arg_phi)
            theta = patch_ib(i)%kin_theta_mean + amp*patch_ib(i)%kin_theta0*sin(arg_th)
            thetad = damp*patch_ib(i)%kin_theta0*sin(arg_th) + amp*patch_ib(i)%kin_theta0*omega*cos(arg_th)
        end if

        ! centroid relative to the hinge: Rx(phi) Ry(theta) offset. kin_offset points to the case-file centroid; a
        ! centroid moved to the centre of mass (geometries 4, 5, 11, 12) sits centroid_offset short of it
        r = patch_ib(i)%kin_offset - patch_ib(i)%centroid_offset
        c(1) = cos(theta)*r(1) + sin(theta)*r(3)
        c(2) = r(2)
        c(3) = -sin(theta)*r(1) + cos(theta)*r(3)
        v(1) = c(1)
        v(2) = cos(phi)*c(2) - sin(phi)*c(3)
        v(3) = sin(phi)*c(2) + cos(phi)*c(3)
        c = v

        ! lab-frame angular velocity: phi' e_x + theta' Rx(phi) e_y
        w(1) = phid
        w(2) = thetad*cos(phi)
        w(3) = thetad*sin(phi)
        call s_cross_product(w, c, v)

        patch_ib(i)%angles(1) = phi
        patch_ib(i)%angles(2) = theta
        patch_ib(i)%angles(3) = 0._wp
        patch_ib(i)%x_centroid = patch_ib(i)%kin_hinge(1) + c(1)
        patch_ib(i)%y_centroid = patch_ib(i)%kin_hinge(2) + c(2)
        patch_ib(i)%z_centroid = patch_ib(i)%kin_hinge(3) + c(3)
        patch_ib(i)%vel = v
        patch_ib(i)%angular_vel = w

    end subroutine s_prescribed_kinematics

    !> Add direction idir's fluid-body face fluxes to each IB. These are the fluxes the RHS differences, Riemann (flux_rsx_vf) plus
    !! viscous, capillary, conductive and diffusive (flux_src_rsx_vf), so each body receives exactly what the fluid loses. A face
    !! counts on the rank owning its fluid cell. With IB forces: force (1:3) and torque about the centroid (4:6) of each local IB
    !! into ib_face_ft, from the last RHS evaluation (idir = 1 starts afresh). With ib_exchange: energy in and mass out of each
    !! lumped body (every body with ib_heat_wrt) into ib_heat(2:3), weighted so a step's stages sum to its exchange (SSP-RK: 1; 1/2,
    !! 1/2; 1/6, 1/6, 2/3).
    subroutine s_accumulate_ib_face_fluxes(idir, stage)

        integer, intent(in)               :: idir, stage
        real(wp), dimension(3), parameter :: w_rk3 = [1._wp/6._wp, 1._wp/6._wp, 2._wp/3._wp]
        real(wp), dimension(3)            :: f, r, tq
        real(wp)                          :: area, w, f_E, f_m
        integer                           :: i, j, k, l, jn, kn, ln, ci, cj, ck, side, ib_id, ib_tmp, xp, yp, zp, g
        integer                           :: j0, k0, l0
        logical                           :: has_src, has_src_E, do_ft, all_ibs

        if (num_ibs == 0) return
        do_ft = moving_immersed_boundary_flag .or. ib_state_wrt
        all_ibs = ib_heat_wrt

        if (do_ft .and. idir == 1) then
            $:GPU_PARALLEL_LOOP(private='[i, l]', collapse=2)
            do i = 1, num_ibs
                do l = 1, 6
                    ib_face_ft(l, i) = 0._wp
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! The conditions under which s_compute_additional_physics_rhs differences flux_src for the momentum and the energy
        has_src = viscous .or. surface_tension .or. heat_conduction
        has_src_E = has_src .or. chem_params%diffusion
        w = 1._wp
        if (time_stepper == time_stepper_rk2) w = 0.5_wp
        if (time_stepper == time_stepper_rk3) w = w_rk3(stage)
        j0 = 0; k0 = 0; l0 = 0
        if (idir == 1) j0 = -1
        if (idir == 2) k0 = -1
        if (idir == 3) l0 = -1

        ! Face (j, k, l) lies between that cell and the next one along idir
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, jn, kn, ln, ci, cj, ck, side, ib_id, ib_tmp, xp, yp, zp, g, area, &
                            & f, r, tq, f_E, f_m]')
        do l = l0, p
            do k = k0, n
                do j = j0, m
                    jn = j; kn = k; ln = l
                    if (idir == 1) jn = j + 1
                    if (idir == 2) kn = k + 1
                    if (idir == 3) ln = l + 1
                    ! side = +1: body on the far side, the face flux enters it; -1: body on the near side, the flux leaves it
                    side = 0
                    if (ib_markers%sf(j, k, l) == 0 .and. ib_markers%sf(jn, kn, ln) /= 0) then
                        if (j >= 0 .and. k >= 0 .and. l >= 0) side = 1
                    else if (ib_markers%sf(j, k, l) /= 0 .and. ib_markers%sf(jn, kn, ln) == 0) then
                        if (jn <= m .and. kn <= n .and. ln <= p) side = -1
                    end if
                    if (side /= 0) then
                        ci = jn; cj = kn; ck = ln
                        if (side == -1) then
                            ci = j; cj = k; ck = l
                        end if
                        call s_decode_patch_periodicity(ib_markers%sf(ci, cj, ck), ib_tmp, xp, yp, zp)
                        call s_get_neighborhood_idx(ib_tmp, ib_id)
                        if (ib_id > 0) then
                            if (idir == 1) area = dy(k)
                            if (idir == 2) area = dx(j)
                            if (idir == 3) area = dx(j)*dy(k)
                            if (num_dims == 3 .and. idir /= 3) area = area*dz(l)
                            if (do_ft) then
                                ! Lever arm: face centre minus the centroid of the body's periodic image
                                r = 0._wp
                                r(1) = merge(x_cb(j), x_cc(j), idir == 1) - (patch_ib(ib_id)%x_centroid + real(xp, &
                                  & wp)*(glb_bounds(1)%end - glb_bounds(1)%beg))
                                r(2) = merge(y_cb(k), y_cc(k), idir == 2) - (patch_ib(ib_id)%y_centroid + real(yp, &
                                  & wp)*(glb_bounds(2)%end - glb_bounds(2)%beg))
                                if (num_dims == 3) r(3) = merge(z_cb(l), z_cc(l), &
                                    & idir == 3) - (patch_ib(ib_id)%z_centroid + real(zp, &
                                    & wp)*(glb_bounds(3)%end - glb_bounds(3)%beg))
                                f = 0._wp
                                do i = 1, num_dims
                                    f(i) = flux_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                    if (has_src) f(i) = f(i) + flux_src_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                    f(i) = real(side, wp)*area*f(i)
                                end do
                                call s_cross_product(r, f, tq)
                                do i = 1, 3
                                    $:GPU_ATOMIC(atomic='update')
                                    ib_face_ft(i, ib_id) = ib_face_ft(i, ib_id) + f(i)
                                    $:GPU_ATOMIC(atomic='update')
                                    ib_face_ft(i + 3, ib_id) = ib_face_ft(i + 3, ib_id) + tq(i)
                                end do
                            end if
                            if (all_ibs .or. patch_ib(ib_id)%thermal_bc == 3) then
                                f_E = flux_rsx_vf(j, k, l, eqn_idx%E)
                                if (has_src_E) f_E = f_E + flux_src_rsx_vf(j, k, l, eqn_idx%E)
                                f_m = 0._wp
                                do i = eqn_idx%cont%beg, eqn_idx%cont%end
                                    f_m = f_m + flux_rsx_vf(j, k, l, i)
                                end do
                                g = patch_ib(ib_id)%gbl_patch_id
                                $:GPU_ATOMIC(atomic='update')
                                ib_heat(2, g) = ib_heat(2, g) + real(side, wp)*w*area*f_E
                                $:GPU_ATOMIC(atomic='update')
                                ib_heat(3, g) = ib_heat(3, g) - real(side, wp)*w*area*f_m
                            end if
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_accumulate_ib_face_fluxes

    !> Set each IB's force and torque: the fluid's share from the face fluxes of the last RHS evaluation
    !! (s_accumulate_ib_face_fluxes), plus collisions and body forces
    subroutine s_compute_ib_forces()

        integer                         :: i, l
        real(wp), dimension(num_ibs, 3) :: forces, torques

        call nvtxStartRange("COMPUTE-IB-FORCES")

        forces = 0._wp
        torques = 0._wp

        ! no kernel over zero-size mapped arrays: a rank can hold no patch after a handoff
        if (num_ibs > 0) then
            $:GPU_PARALLEL_LOOP(private='[i, l]', copy='[forces, torques]')
            do i = 1, num_ibs
                do l = 1, 3
                    forces(i, l) = ib_face_ft(l, i)
                    torques(i, l) = ib_face_ft(l + 3, i)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        call nvtxStartRange("APPLY-COLLISION-FORCES")
        call s_apply_collision_forces(ghost_points, num_gps, ib_markers, forces, torques)
        call nvtxEndRange

        ! reduce the forces across local neighborhood ranks
        call nvtxStartRange("COMMUNICATE-IB-FORCES")
        call s_communicate_ib_forces(forces, torques)
        call nvtxEndRange

        ! Body forces after reducing, to avoid double counting. Only a two-way body feels them; a fixed or prescribed body keeps
        ! the fluid and collision force, and its mass may be unset (-1e6)
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm >= 2) then
                if (bf_x) forces(i, 1) = forces(i, 1) + accel_bf(1)*patch_ib(i)%mass
                if (bf_y) forces(i, 2) = forces(i, 2) + accel_bf(2)*patch_ib(i)%mass
                if (bf_z) forces(i, 3) = forces(i, 3) + accel_bf(3)*patch_ib(i)%mass
            end if
        end do

        ! apply the summed forces
        if (num_ibs > 0) then
            $:GPU_PARALLEL_LOOP(private='[i, l]', copyin='[forces, torques]')
            do i = 1, num_ibs
                do l = 1, 3
                    patch_ib(i)%force(l) = forces(i, l)
                    patch_ib(i)%torque(l) = torques(i, l)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        call nvtxEndRange

    end subroutine s_compute_ib_forces

    !> Computes the center of mass for IB patch types where we are unable to determine their center of mass analytically.
    !> These patches include things like NACA airfoils and STL models
    !> Move a moving airfoil/STL patch's centroid to its centre of mass, keeping the difference in centroid_offset. Collective:
    !! every rank calls it for every global id, contributing zeros for patches it does not hold.
    subroutine s_compute_centroid_offset(gid, ib_marker)

        integer, intent(in)      :: gid        !< global patch id
        integer, intent(in)      :: ib_marker  !< local index on this rank; <= 0 if not held
        integer                  :: i, j, k, num_cells_local, decoded_gbl_id
        integer(kind=8)          :: num_cells, needs_loc, needs_glb
        real(wp), dimension(1:3) :: center_of_mass, center_of_mass_local

        needs_loc = 0_8
        if (ib_marker > 0) then
            if (f_needs_centroid_offset(patch_ib(ib_marker))) needs_loc = 1_8
        end if
        call s_mpi_allreduce_integer_sum(needs_loc, needs_glb)
        if (needs_glb == 0_8) then
            if (ib_marker > 0) patch_ib(ib_marker)%centroid_offset(:) = 0._wp
            return
        end if

        center_of_mass_local = [0._wp, 0._wp, 0._wp]
        num_cells_local = 0
        if (ib_marker > 0) then
            ! get the summed mass distribution and number of cells to divide by
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) /= 0) then
                            call s_decode_patch_periodicity(ib_markers%sf(i, j, k), decoded_gbl_id)
                            if (decoded_gbl_id == gid) then
                                num_cells_local = num_cells_local + 1
                                center_of_mass_local = center_of_mass_local + [x_cc(i), y_cc(j), 0._wp]
                                if (num_dims == 3) center_of_mass_local(3) = center_of_mass_local(3) + z_cc(k)
                            end if
                        end if
                    end do
                end do
            end do
        end if

        ! reduce the mass contribution over all MPI ranks and compute COM
        call s_mpi_allreduce_integer_sum(int(num_cells_local, 8), num_cells)
        call s_mpi_allreduce_sum(center_of_mass_local(1), center_of_mass(1))
        call s_mpi_allreduce_sum(center_of_mass_local(2), center_of_mass(2))
        call s_mpi_allreduce_sum(center_of_mass_local(3), center_of_mass(3))
        if (ib_marker <= 0) return
        if (num_cells == 0) then
            patch_ib(ib_marker)%centroid_offset = [0._wp, 0._wp, 0._wp]
            return
        end if
        center_of_mass = center_of_mass/real(num_cells, wp)

        ! assign the centroid offset as a vector pointing from the true COM to the "centroid" in the input file and replace the
        ! current centroid
        patch_ib(ib_marker)%centroid_offset = [patch_ib(ib_marker)%x_centroid, patch_ib(ib_marker)%y_centroid, &
                 & patch_ib(ib_marker)%z_centroid] - center_of_mass
        patch_ib(ib_marker)%x_centroid = center_of_mass(1)
        patch_ib(ib_marker)%y_centroid = center_of_mass(2)
        patch_ib(ib_marker)%z_centroid = center_of_mass(3)

        ! rotate the centroid offset back into the local coords of the IB
        patch_ib(ib_marker)%centroid_offset = matmul(patch_ib(ib_marker)%rotation_matrix_inverse, &
                 & patch_ib(ib_marker)%centroid_offset)

    end subroutine s_compute_centroid_offset

    !> A moving airfoil or STL (geometries 4, 5, 11, 12) moves its centroid to the centre of mass and keeps a centroid_offset
    pure logical function f_needs_centroid_offset(patch)

        type(ib_patch_parameters), intent(in) :: patch

        f_needs_centroid_offset = patch%moving_ibm /= 0 .and. any(patch%geometry == [4, 5, 11, 12])

    end function f_needs_centroid_offset

    !> On restart, replace the re-measured centroid offsets with the ones the run was using (restart_data/ib_offset_<step>.dat,
    !! written with each checkpoint), and re-place kinematics-driven bodies about them. No file: the measured offsets stand.
    impure subroutine s_restore_centroid_offsets(t_init)

        real(wp), intent(in)                 :: t_init
        character(len=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist
        integer                              :: gid, i, ios, file_unit, step
        real(wp), dimension(3)               :: off

        step = t_step_start
        if (cfl_dt) step = n_start
        if (step == 0) return
        write (file_loc, '(A,I0,A)') trim(case_dir) // '/restart_data/ib_offset_', step, '.dat'
        inquire (file=trim(file_loc), exist=file_exist)
        if (.not. file_exist) return
        open (newunit=file_unit, file=trim(file_loc), status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read (file_unit, *, iostat=ios) gid, off
            if (ios /= 0) exit
            call s_get_neighborhood_idx(gid, i)
            if (i > 0) then
                patch_ib(i)%centroid_offset = off
                if (patch_ib(i)%moving_ibm /= 0 .and. patch_ib(i)%kin_model > 0) call s_prescribed_kinematics(i, t_init)
            end if
        end do
        close (file_unit)

    end subroutine s_restore_centroid_offsets

    !> Computes the moment of inertia for an immersed boundary
    subroutine s_compute_moment_of_inertia(patch, axis, moment)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ib_patch_parameters), intent(in) :: patch
        real(wp), dimension(3), intent(in)    :: axis
        real(wp), intent(out)                 :: moment
        real(wp)                              :: distance_to_axis, cell_volume
        real(wp), dimension(3)                :: position, closest_point_along_axis, vector_to_axis, normal_axis
        integer                               :: i, j, k, count, ib_marker

        ! if the IB is in 2D or a 3D sphere, we can compute this exactly
        if (patch%geometry == 2) then  ! circle
            moment = 0.5_wp*patch%mass*(patch%radius)**2
        else if (patch%geometry == 3) then  ! rectangle
            moment = patch%mass*(patch%length_x**2 + patch %length_y**2)/6._wp
        else if (patch%geometry == 6) then  ! ellipse
            moment = 0.0625_wp*patch%mass*(patch%length_x**2 + patch %length_y**2)
        else if (patch%geometry == 8) then  ! sphere
            moment = 0.4*patch%mass*(patch%radius)**2
        else  ! we do not have an analytic moment of inertia calculation and need to approximate it directly via a sum
            count = 0
            cell_volume = (x_cc(1) - x_cc(0))*(y_cc(1) - y_cc(0))
            ! computed without grid stretching. Update in the loop to perform with stretching
            if (p /= 0) then
                cell_volume = cell_volume*(z_cc(1) - z_cc(0))
            end if

            ib_marker = patch%gbl_patch_id

            if (p == 0) then
                normal_axis = [0, 0, 1]
            else if (sqrt(sum(axis**2)) < sgm_eps) then
                ! if the object is not actually rotating at this time, return a dummy value and exit
                moment = 1._wp
                return
            else
                normal_axis = axis/sqrt(sum(axis**2))
            end if

            moment = 0._wp

            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) == ib_marker) then
                            count = count + 1  ! increment the count of total cells in the boundary

                            ! get the position in local coordinates so that the axis passes through 0, 0, 0
                            if (num_dims < 3) then
                                position = [x_cc(i), y_cc(j), 0._wp] - [patch%x_centroid, patch%y_centroid, 0._wp]
                            else
                                position = [x_cc(i), y_cc(j), z_cc(k)] - [patch%x_centroid, patch%y_centroid, patch%z_centroid]
                            end if

                            ! project the position along the axis to find the closest distance to the rotation axis
                            closest_point_along_axis = normal_axis*dot_product(normal_axis, position)
                            vector_to_axis = position - closest_point_along_axis
                            distance_to_axis = dot_product(vector_to_axis, vector_to_axis)  ! saves the distance to the axis squared

                            ! compute the position component of the moment
                            moment = moment + distance_to_axis
                        end if
                    end do
                end do
            end do

            ! write the final moment assuming the points are all uniform density
            moment = moment*patch%mass/(count*cell_volume)
        end if

    end subroutine s_compute_moment_of_inertia

    !> Wrap immersed boundary positions across periodic domain boundaries
    subroutine s_wrap_periodic_ibs()

        integer :: patch_id

        $:GPU_PARALLEL_LOOP(private='[patch_id]')
        do patch_id = 1, num_ibs
            ! check domain wraps in x, y,
            #:for X, DIR in [('x', 1), ('y', 2), ('z', 3)]
                if (num_dims >= ${DIR}$) then
                    ! check for periodicity
                    if (ib_bc_${X}$%beg == BC_PERIODIC) then
                        ! check if the boundary has left the domain, and then correct
                        if (patch_ib(patch_id)%${X}$_centroid < glb_bounds(${DIR}$)%beg) then
                            ! if the boundary exited "left", wrap it back around to the "right"
                            patch_ib(patch_id)%${X}$_centroid = patch_ib(patch_id)%${X}$_centroid + (glb_bounds(${DIR}$)%end &
                                     & - glb_bounds(${DIR}$)%beg)
                        else if (patch_ib(patch_id)%${X}$_centroid > glb_bounds(${DIR}$)%end) then
                            ! if the boundary exited "right", wrap it back around to the "left"
                            patch_ib(patch_id)%${X}$_centroid = patch_ib(patch_id)%${X}$_centroid - (glb_bounds(${DIR}$)%end &
                                     & - glb_bounds(${DIR}$)%beg)
                        end if
                    end if
                end if
            #:endfor
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_wrap_periodic_ibs

    !> Builds the ascending list of distinct neighborhood ranks used by s_communicate_ib_forces and allocates its buffers once. The
    !! ascending order is the fixed order in which an owner adds its neighbors' contributions.
    impure subroutine s_setup_ib_force_comm()

#ifdef MFC_MPI
        integer                            :: i, j, r, max_recs
        integer, allocatable, dimension(:) :: flat

        if (num_procs == 1 .or. allocated(ib_nbrs)) return
        @:PROHIBIT(num_gbl_ibs >= min(2_8**digits(0._wp), int(huge(0), 8)), &
                   & "Too many IBs to encode their ids exactly in real(wp) force messages")

        ! a rank can fill several table slots (periodicity, few ranks) or be its own neighbor; keep each distinct rank once
        flat = reshape(ib_neighbor_ranks, [size(ib_neighbor_ranks)])
        allocate (ib_nbrs(size(flat)))
        n_ib_nbrs = 0
        do i = 1, size(flat)
            r = flat(i)
            if (r < 0 .or. r == proc_rank) cycle
            if (any(ib_nbrs(1:n_ib_nbrs) == r)) cycle
            ! insert in ascending order
            j = n_ib_nbrs
            do while (j > 0)
                if (ib_nbrs(j) < r) exit
                ib_nbrs(j + 1) = ib_nbrs(j)
                j = j - 1
            end do
            ib_nbrs(j + 1) = r
            n_ib_nbrs = n_ib_nbrs + 1
        end do

        ! an owner never receives more records than it owns, and a broadcast never carries more than its sender owns
        max_recs = min(size(patch_ib), num_local_ibs_max)
        allocate (ib_send_counts(n_ib_nbrs), ib_reqs(2*n_ib_nbrs), ib_stats(MPI_STATUS_SIZE, 2*n_ib_nbrs))
        allocate (ib_send_bufs(ib_rec_len*max_recs, max(n_ib_nbrs, 1)), ib_recv_bufs(ib_rec_len*max_recs, max(n_ib_nbrs, 1)))
#endif

    end subroutine s_setup_ib_force_comm

    !> Position of rank in ib_nbrs
    function f_ib_nbr_slot(rank) result(slot)

        integer, intent(in) :: rank
        integer             :: slot

        do slot = 1, n_ib_nbrs
            if (ib_nbrs(slot) == rank) return
        end do
        call s_mpi_abort('IB owner is not a neighborhood rank')

    end function f_ib_nbr_slot

    !> Reduces forces and torques so every rank holding an IB ends with bit-identical values. Each IB has one owner, the rank whose
    !! subdomain holds its centroid (patch_ib%owner_rank, stamped at handoff). Phase 1: each rank sends its nonzero partials to the
    !! owner, which adds them to its own partial in ascending source-rank order once all have arrived. Phase 2: each owner sends its
    !! totals to every neighborhood rank, which overwrite their partials. One nonblocking message per neighbor per phase.
    subroutine s_communicate_ib_forces(forces, torques)

        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques

#ifdef MFC_MPI
        integer :: i, j, s, r, slot, pos, nvals, ierr

        if (num_procs == 1 .or. n_ib_nbrs == 0) return

        ! Phase 1: send each partial to the owner of its IB
        do s = 1, n_ib_nbrs
            call MPI_IRECV(ib_recv_bufs(:,s), size(ib_recv_bufs, 1), mpi_p, ib_nbrs(s), 610, MPI_COMM_WORLD, ib_reqs(s), ierr)
        end do

        ib_send_counts = 0
        do i = 1, num_ibs
            if (patch_ib(i)%owner_rank == proc_rank) cycle
            if (all(forces(i,:) == 0._wp) .and. all(torques(i,:) == 0._wp)) cycle
            slot = f_ib_nbr_slot(patch_ib(i)%owner_rank)
            pos = ib_rec_len*ib_send_counts(slot)
            ib_send_bufs(pos + 1, slot) = real(patch_ib(i)%gbl_patch_id, wp)
            ib_send_bufs(pos + 2:pos + 4,slot) = forces(i,1:3)
            ib_send_bufs(pos + 5:pos + 7,slot) = torques(i,1:3)
            ib_send_counts(slot) = ib_send_counts(slot) + 1
        end do

        ! every neighbor gets a message, possibly empty, so every posted receive completes
        do s = 1, n_ib_nbrs
            call MPI_ISEND(ib_send_bufs(:,s), ib_rec_len*ib_send_counts(s), mpi_p, ib_nbrs(s), 610, MPI_COMM_WORLD, &
                           & ib_reqs(n_ib_nbrs + s), ierr)
        end do
        call MPI_WAITALL(2*n_ib_nbrs, ib_reqs, ib_stats, ierr)

        ! owner sum: own partial first, then each neighbor in ascending rank order, whatever order the messages arrived in
        do s = 1, n_ib_nbrs
            call MPI_GET_COUNT(ib_stats(:,s), mpi_p, nvals, ierr)
            do r = 0, nvals/ib_rec_len - 1
                pos = ib_rec_len*r
                call s_get_neighborhood_idx(nint(ib_recv_bufs(pos + 1, s)), j)
                @:ASSERT(j > 0, 'IB force contribution for an IB this rank does not hold')
                @:ASSERT(patch_ib(j)%owner_rank == proc_rank, 'IB force contribution for an IB this rank does not own')
                forces(j,1:3) = forces(j,1:3) + ib_recv_bufs(pos + 2:pos + 4,s)
                torques(j,1:3) = torques(j,1:3) + ib_recv_bufs(pos + 5:pos + 7,s)
            end do
        end do

        ! Phase 2: each owner sends its totals to every neighbor
        do s = 1, n_ib_nbrs
            call MPI_IRECV(ib_recv_bufs(:,s), size(ib_recv_bufs, 1), mpi_p, ib_nbrs(s), 611, MPI_COMM_WORLD, ib_reqs(s), ierr)
        end do

        do i = 1, num_local_ibs
            j = local_ib_patch_ids(i)
            pos = ib_rec_len*(i - 1)
            ib_send_bufs(pos + 1, 1) = real(patch_ib(j)%gbl_patch_id, wp)
            ib_send_bufs(pos + 2:pos + 4,1) = forces(j,1:3)
            ib_send_bufs(pos + 5:pos + 7,1) = torques(j,1:3)
        end do

        do s = 1, n_ib_nbrs
            call MPI_ISEND(ib_send_bufs(:,1), ib_rec_len*num_local_ibs, mpi_p, ib_nbrs(s), 611, MPI_COMM_WORLD, &
                           & ib_reqs(n_ib_nbrs + s), ierr)
        end do
        call MPI_WAITALL(2*n_ib_nbrs, ib_reqs, ib_stats, ierr)

        ! holders take the owner's bits
        do s = 1, n_ib_nbrs
            call MPI_GET_COUNT(ib_stats(:,s), mpi_p, nvals, ierr)
            do r = 0, nvals/ib_rec_len - 1
                pos = ib_rec_len*r
                call s_get_neighborhood_idx(nint(ib_recv_bufs(pos + 1, s)), j)
                if (j <= 0) cycle
                @:ASSERT(patch_ib(j)%owner_rank == ib_nbrs(s), 'IB force total from a rank that does not own the IB')
                forces(j,1:3) = ib_recv_bufs(pos + 2:pos + 4,s)
                torques(j,1:3) = ib_recv_bufs(pos + 5:pos + 7,s)
            end do
        end do
#endif

    end subroutine s_communicate_ib_forces

    !> @brief Swaps ownership of IBs and passes ownership of IBs to neighbor processors
    subroutine s_handoff_ib_ownership()

        integer                                         :: i, j, k, output_idx, local_output_idx
        integer                                         :: old_num_local_ibs, num_ibs_old, num_ibs_pre
        integer                                         :: new_count, recv_count
        integer                                         :: pack_pos, unpack_pos, buf_size, patch_bytes
        integer                                         :: send_neighbor, recv_neighbor, ierr
        integer                                         :: dx, dy, dz, tag, nbr_idx, nreqs
        real(wp), dimension(3)                          :: centroid
        logical                                         :: is_new
        type(ib_patch_parameters)                       :: tmp_patch
        integer, dimension(num_local_ibs_max)           :: local_ib_idx_old
        integer, dimension(num_ib_patches_max_namelist) :: old_to_new  ! old patch_ib slot -> slot after compaction
        ! (2R+1)^num_dims - 1 neighbors in the radius-R neighborhood; each gets its own recv buffer
        integer                       :: max_nbrs
        character(len=1), allocatable :: send_buf(:), recv_bufs(:,:)
        integer, allocatable          :: requests(:), recv_neighbor_list(:)

#ifdef MFC_MPI
        call nvtxStartRange("HANDOFF-IB-OWNERSHIP")

        if (num_procs > 1) then
            ! save a copy of the local IB's global indices to cross-reference for later.
            local_ib_idx_old = 0
            old_num_local_ibs = num_local_ibs
            do i = 1, num_local_ibs
                local_ib_idx_old(i) = patch_ib(local_ib_patch_ids(i))%gbl_patch_id
            end do

            ! Sync GPU-updated fields (angles, angular_vel, centroids) to host before
            ! compaction and MPI packing, which read from host memory.
            $:GPU_UPDATE(host='[patch_ib]')

            ! delete any particles that no longer need to be tracked and coalesce the array
            num_ibs_old = num_ibs
            output_idx = 0
            local_output_idx = 0
            do i = 1, num_ibs
                old_to_new(i) = -1
                centroid = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, 0._wp]
                if (num_dims == 3) centroid(3) = patch_ib(i)%z_centroid

                ! delete if not in neighborhood
                if (f_neighborhood_ranks_own_location(centroid)) then
                    output_idx = output_idx + 1
                    old_to_new(i) = output_idx
                    if (i /= output_idx) then
                        patch_ib(output_idx) = patch_ib(i)
                    end if

                    ! check if in local domain
                    if (f_local_rank_owns_location(centroid, glb_bounds)) then
                        local_output_idx = local_output_idx + 1
                        @:PROHIBIT(local_output_idx > num_local_ibs_max, &
                                   & "Too many IBs on a single processor rank. Modify case file or increase limit of num_local_ibs_max to resolve.")
                        local_ib_patch_ids(local_output_idx) = output_idx
                        patch_ib(output_idx)%owner_rank = proc_rank
                    end if
                end if
            end do
            num_ibs = output_idx
            num_local_ibs = local_output_idx
            ! num_ibs shrinks here, so refresh it with patch_ib: s_get_neighborhood_idx bounds its search by it on the device
            $:GPU_UPDATE(device='[patch_ib, num_ibs]')
            call s_compact_ib_lookup(old_to_new, num_ibs_old)

            ! Broadcast newly-owned patches to all neighborhood neighbors
            patch_bytes = storage_size(tmp_patch)/8
            ! a rank can hand off at most every global patch, so size by num_gbl_ibs, not num_local_ibs_max
            buf_size = storage_size(0)/8 + patch_bytes*max(1, min(num_local_ibs_max, num_gbl_ibs))
            max_nbrs = (2*ib_neighborhood_radius + 1)**num_dims - 1
            allocate (send_buf(buf_size), recv_bufs(buf_size, max_nbrs), requests(2*max_nbrs), recv_neighbor_list(max_nbrs))

            ! Write placeholder count at position 0
            pack_pos = 0
            call MPI_PACK(0, 1, MPI_INTEGER, send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)

            ! pack new patches and count them
            new_count = 0
            do i = 1, num_local_ibs
                k = local_ib_patch_ids(i)
                is_new = .true.
                do j = 1, old_num_local_ibs
                    if (patch_ib(k)%gbl_patch_id == local_ib_idx_old(j)) then
                        is_new = .false.
                        exit
                    end if
                end do
                if (is_new) then
                    call MPI_PACK(patch_ib(k), patch_bytes, MPI_BYTE, send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)
                    new_count = new_count + 1
                end if
            end do

            ! Overwrite the placeholder with the real count
            pack_pos = 0
            call MPI_PACK(new_count, 1, MPI_INTEGER, send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)
            pack_pos = storage_size(0)/8 + new_count*patch_bytes

            ! Post all receives first, then sends
            nreqs = 0
            nbr_idx = 0
            do dz = merge(-ib_neighborhood_radius, 0, num_dims == 3), merge(ib_neighborhood_radius, 0, num_dims == 3)
                do dy = -ib_neighborhood_radius, ib_neighborhood_radius
                    do dx = -ib_neighborhood_radius, ib_neighborhood_radius
                        if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
                        nbr_idx = nbr_idx + 1
                        ! one tag per offset in the (2R+1)^3 neighbourhood: a radix-3 encoding collides once R > 1
                        tag = 200 + ((dx + ib_neighborhood_radius)*(2*ib_neighborhood_radius + 1) + (dy + ib_neighborhood_radius)) &
                                     & *(2*ib_neighborhood_radius + 1) + (dz + ib_neighborhood_radius)
                        recv_neighbor = ib_neighbor_ranks(-dx, -dy, -dz)
                        recv_neighbor_list(nbr_idx) = MPI_PROC_NULL
                        if (recv_neighbor < 0) cycle
                        recv_neighbor_list(nbr_idx) = recv_neighbor
                        nreqs = nreqs + 1
                        call MPI_IRECV(recv_bufs(:,nbr_idx), buf_size, MPI_PACKED, recv_neighbor, tag, MPI_COMM_WORLD, &
                                       & requests(nreqs), ierr)
                    end do
                end do
            end do

            do dz = merge(-ib_neighborhood_radius, 0, num_dims == 3), merge(ib_neighborhood_radius, 0, num_dims == 3)
                do dy = -ib_neighborhood_radius, ib_neighborhood_radius
                    do dx = -ib_neighborhood_radius, ib_neighborhood_radius
                        if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
                        ! one tag per offset in the (2R+1)^3 neighbourhood: a radix-3 encoding collides once R > 1
                        tag = 200 + ((dx + ib_neighborhood_radius)*(2*ib_neighborhood_radius + 1) + (dy + ib_neighborhood_radius)) &
                                     & *(2*ib_neighborhood_radius + 1) + (dz + ib_neighborhood_radius)
                        send_neighbor = ib_neighbor_ranks(dx, dy, dz)
                        if (send_neighbor < 0) cycle
                        nreqs = nreqs + 1
                        call MPI_ISEND(send_buf, pack_pos, MPI_PACKED, send_neighbor, tag, MPI_COMM_WORLD, requests(nreqs), ierr)
                    end do
                end do
            end do

            call MPI_WAITALL(nreqs, requests, MPI_STATUSES_IGNORE, ierr)

            ! Unpack all received buffers
            num_ibs_pre = num_ibs
            do nbr_idx = 1, max_nbrs
                if (recv_neighbor_list(nbr_idx) == MPI_PROC_NULL) cycle
                unpack_pos = 0
                call MPI_UNPACK(recv_bufs(:,nbr_idx), buf_size, unpack_pos, recv_count, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                do i = 1, recv_count
                    call MPI_UNPACK(recv_bufs(:,nbr_idx), buf_size, unpack_pos, tmp_patch, patch_bytes, MPI_BYTE, MPI_COMM_WORLD, &
                                    & ierr)
                    call s_get_neighborhood_idx(tmp_patch%gbl_patch_id, j, num_ibs_pre)
                    if (j < 0) then
                        num_ibs = num_ibs + 1
                        @:ASSERT(num_ibs <= size(patch_ib), 'patch_ib overflow in neighborhood handoff')
                        patch_ib(num_ibs) = tmp_patch
                    else
                        ! an IB we already hold changed owner: record the new owner, keep our own copy of its state
                        patch_ib(j)%owner_rank = tmp_patch%owner_rank
                    end if
                end do
            end do

            deallocate (send_buf, recv_bufs, requests, recv_neighbor_list)
            $:GPU_UPDATE(device='[patch_ib, num_ibs]')
            call s_merge_ib_lookup(num_ibs_pre)
        end if

        call nvtxEndRange()
#endif

    end subroutine s_handoff_ib_ownership

    !> Abort if any namelist immersed boundary marks no cell anywhere. The owning rank is chosen by the patch centroid, so an STL
    !! placed with model_translate while its centroid is left elsewhere is handed to ranks that do not hold it, and it vanishes
    !! silently. Cloud particles are skipped: a particle's centroid is its body, so it cannot vanish this way.
    impure subroutine s_check_every_patch_marked()

        integer(kind=8), allocatable :: cnt_loc(:), cnt_glb(:)
        integer                      :: gid, i, j, k, num_checked

        num_checked = num_gbl_ibs - sum(particle_cloud(1:num_particle_clouds)%num_particles)
        if (num_checked <= 0) return
        allocate (cnt_loc(num_checked), cnt_glb(num_checked))
        cnt_loc = 0_8
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    if (ib_markers%sf(i, j, k) /= 0) then
                        call s_decode_patch_periodicity(ib_markers%sf(i, j, k), gid)
                        if (gid <= num_checked) cnt_loc(gid) = cnt_loc(gid) + 1_8
                    end if
                end do
            end do
        end do
        call s_mpi_allreduce_integer_sum_vec(cnt_loc, cnt_glb)
        @:PROHIBIT(any(cnt_glb == 0_8), &
                   & "An immersed boundary marked no cell anywhere: its centroid decides which rank owns it, so a body placed " &
                   & // "elsewhere with model_translate is handed to a rank that does not hold it. Set patch_ib%x/y/z_centroid " &
                   & // "to where the body actually is.")
        deallocate (cnt_loc, cnt_glb)

    end subroutine s_check_every_patch_marked

    !> Finalize the IBM module
    impure subroutine s_finalize_ibm_module()

        integer :: i

        @:DEALLOCATE(ib_markers%sf)
        @:DEALLOCATE(corrected_gps%sf)
        if (associated(ib_wall_theta%sf)) then
            do i = 1, ibw_a
                @:DEALLOCATE(ib_wall_src(i)%sf)
            end do
            @:DEALLOCATE(ib_wall_theta%sf)
        end if
        @:DEALLOCATE(ib_face_ft)
        do i = 1, num_ib_airfoils_max
            if (allocated(ib_airfoil_grids(i)%upper)) then
                @:DEALLOCATE(ib_airfoil_grids(i)%upper)
                @:DEALLOCATE(ib_airfoil_grids(i)%lower)
            end if
        end do
        if (allocated(models)) then
            @:DEALLOCATE(models)
        end if
        if (allocated(ghost_points)) then
            @:DEALLOCATE(ghost_points)
            @:DEALLOCATE(gp_prim_new)
        end if
        if (allocated(gp_surf)) then
            @:DEALLOCATE(gp_surf)
        end if
        if (allocated(ib_heat)) then
            @:DEALLOCATE(ib_heat)
        end if
        if (collision_model > 0) call s_finalize_collisions_module()
#ifdef MFC_MPI
        if (allocated(ib_nbrs)) deallocate (ib_nbrs, ib_send_counts, ib_reqs, ib_stats, ib_send_bufs, ib_recv_bufs)
#endif

    end subroutine s_finalize_ibm_module

    !> ib_flux_wall: allocate the wall fields, give each reacting body's ghost points with a fluid face neighbor their wall area per
    !! face area, and report the wall area the faces carry. A planar wall of area A has sum|n_i| A of staircase faces, so each face
    !! carries 1/sum|n_i| of its own area.
    impure subroutine s_ib_flux_wall_setup()

        integer  :: i, j, k, l, q, d, sgn
        logical  :: wet
        real(wp) :: area, area_glb, vol

        do q = 1, ibw_a
            @:ALLOCATE(ib_wall_src(q)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                       & idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(ib_wall_src(q))
        end do
        @:ALLOCATE(ib_wall_theta%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(ib_wall_theta)

        $:GPU_PARALLEL_LOOP(private='[j, k, l, q]', collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    ib_wall_theta%sf(j, k, l) = 1._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, ibw_a
                        ib_wall_src(q)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, d, sgn, wet]', present='[ghost_points]')
        do i = 1, num_gps
            j = ghost_points(i)%loc(1); k = ghost_points(i)%loc(2); l = ghost_points(i)%loc(3)
            wet = .false.
            do d = 1, num_dims
                do sgn = -1, 1, 2
                    if (ib_markers%sf(j + sgn*merge(1, 0, d == 1), k + sgn*merge(1, 0, d == 2), l + sgn*merge(1, 0, &
                        & d == 3)) == 0) wet = .true.
                end do
            end do
            if (wet .and. patch_ib(ghost_points(i)%ib_patch_id)%surface_reaction == 1 &
                & .and. patch_ib(ghost_points(i)%ib_patch_id)%inj_species == 0 .and. sum(abs(ghost_points(i)%levelset_norm)) &
                & > 0._wp) then
                ib_wall_src(ibw_a)%sf(j, k, l) = sqrt(sum(ghost_points(i)%levelset_norm**2))/sum(abs(ghost_points(i)%levelset_norm))
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_ib_wall_src_buffers()

        area = 0._wp
        $:GPU_PARALLEL_LOOP(private='[j, k, l, vol]', collapse=3, reduction='[[area]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (ib_markers%sf(j, k, l) == 0) then
                        vol = dx(j)*dy(k)
                        if (num_dims == 3) vol = vol*dz(l)
                        area = area + vol*f_wall_face_sum(ib_wall_src, j, k, l, 0)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_mpi_allreduce_sum(area, area_glb)
        if (proc_rank == 0) print '(A, ES14.6)', ' ib_flux_wall: reacting wall area carried by the faces = ', area_glb

    end subroutine s_ib_flux_wall_setup

    !> ib_flux_wall, before each Riemann solve (primitive buffers filled): fill ib_wall_src's buffer cells and set the positivity
    !! factor of every fluid cell in the domain and its first buffer layer (a wall face on a rank seam is solved on both sides), the
    !! largest theta <= 1 for which its wall faces take at most ib_wall_safety of any species in one forward-Euler step of size dt
    !! (each SSP-RK stage is a convex combination of such steps). Conservative: the body books the scaled flux.
    impure subroutine s_ib_flux_wall_prepare(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), parameter                                 :: ib_wall_safety = 0.9_wp
        real(wp)                                            :: dt_loc, use, th, rhoY
        integer, dimension(3)                               :: lo, hi
        integer                                             :: j, k, l, q

        call s_ib_wall_src_buffers()
        dt_loc = dt
        lo = 0; hi = [m, n, p]
        do q = 1, num_dims
            lo(q) = -1; hi(q) = hi(q) + 1
        end do
        $:GPU_PARALLEL_LOOP(private='[j, k, l, q, use, th, rhoY]', collapse=3, copyin='[dt_loc, lo, hi]')
        do l = lo(3), hi(3)
            do k = lo(2), hi(2)
                do j = lo(1), hi(1)
                    th = 1._wp
                    if (ib_markers%sf(j, k, l) == 0) then
                        if (f_wall_face_sum(ib_wall_src, j, k, l, 0) > 0._wp) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, num_species
                                use = dt_loc*f_wall_face_sum(ib_wall_src, j, k, l, q)
                                rhoY = q_prim_vf(eqn_idx%cont%beg)%sf(j, k, l)*q_prim_vf(eqn_idx%species%beg + q - 1)%sf(j, k, l)
                                ! Divide only when the species is consumed (use > 0): a trace species with rhoY < 0 and no
                                ! consumption would give 0/0, and min/max with a NaN are processor dependent
                                if (use > 0._wp) then
                                    if (rhoY <= 0._wp) then
                                        th = 0._wp
                                    else if (use > ib_wall_safety*rhoY) then
                                        th = min(th, ib_wall_safety*rhoY/use)
                                    end if
                                end if
                            end do
                        end if
                    end if
                    if (.not. (th >= 0._wp)) th = 0._wp
                    ib_wall_theta%sf(j, k, l) = th
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_flux_wall_prepare

    !> ib_flux_wall: sum over the reacting-wall faces of fluid cell (j, k, l) of the wall area per cell volume, times species q's
    !! consumption rate max(-mdot_q, 0) (q > 0) or 1 (q = 0)
    function f_wall_face_sum(src, j, k, l, q) result(s)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(ibw_a), intent(in) :: src  !< ib_wall_src, passed so the device routine reads the present copy
        integer, intent(in)                              :: j, k, l, q
        real(wp)                                         :: s, w
        real(wp), dimension(3)                           :: h
        integer, dimension(3)                            :: c
        integer                                          :: d, sgn

        s = 0._wp
        h = [dx(j), dy(k), 1._wp]
        if (num_dims == 3) h(3) = dz(l)
        do d = 1, num_dims
            do sgn = -1, 1, 2
                c = [j, k, l]
                c(d) = c(d) + sgn
                w = 1._wp
                if (q > 0) w = max(-src(q)%sf(c(1), c(2), c(3)), 0._wp)
                s = s + src(ibw_a)%sf(c(1), c(2), c(3))*w/h(d)
            end do
        end do

    end function f_wall_face_sum

    !> ib_flux_wall: the surface area a reacting body's ghost point (j, k, l) stands for in its surface record, the wall area its
    !! fluid faces carry (a theta times their area, as the Riemann solver takes them); 0 off the faces
    function f_flux_wall_record_area(src, theta, marks, loc) result(w)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(ibw_a), intent(in) :: src
        type(scalar_field), intent(in)                   :: theta
        type(integer_field), intent(in)                  :: marks  !< ib_markers, passed like src
        integer, dimension(3), intent(in)                :: loc
        real(wp)                                         :: w
        real(wp), dimension(3)                           :: h
        integer, dimension(3)                            :: c
        integer                                          :: d, sgn

        w = 0._wp
        h = [dx(loc(1)), dy(loc(2)), 1._wp]
        if (num_dims == 3) h(3) = dz(loc(3))
        do d = 1, num_dims
            do sgn = -1, 1, 2
                c = loc
                c(d) = c(d) + sgn
                if (marks%sf(c(1), c(2), c(3)) == 0) w = w + src(ibw_a)%sf(loc(1), loc(2), loc(3))*theta%sf(c(1), c(2), &
                    & c(3))*h(1)*h(2)*h(3)/h(d)
            end do
        end do

    end function f_flux_wall_record_area

    !> ib_flux_wall: fill ib_wall_src's buffer cells from the neighbouring ranks or the periodic image; they stay zero (no reacting
    !! wall) at other domain boundaries
    impure subroutine s_ib_wall_src_buffers()

        integer, dimension(3) :: lo, hi, c, dst, src_c
        integer               :: d, loc, edge, nd, j, k, l, q

        do d = 1, num_dims
            do loc = -1, 1, 2
                select case (d)
                case (1)
                    edge = merge(bc_x%beg, bc_x%end, loc == -1); nd = m
                case (2)
                    edge = merge(bc_y%beg, bc_y%end, loc == -1); nd = n
                case default
                    edge = merge(bc_z%beg, bc_z%end, loc == -1); nd = p
                end select
                if (edge >= 0) then
                    call s_mpi_sendrecv_variables_buffers(ib_wall_src, d, loc, ibw_a)
                else if (edge == BC_PERIODIC) then
                    lo = [idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg]
                    hi = [idwbuff(1)%end, idwbuff(2)%end, idwbuff(3)%end]
                    lo(d) = 1; hi(d) = buff_size
                    $:GPU_PARALLEL_LOOP(private='[j, k, l, q, c, dst, src_c]', collapse=3, copyin='[lo, hi, d, loc, nd]')
                    do l = lo(3), hi(3)
                        do k = lo(2), hi(2)
                            do j = lo(1), hi(1)
                                c = [j, k, l]
                                dst = c; src_c = c
                                if (loc == -1) then
                                    dst(d) = -c(d); src_c(d) = nd + 1 - c(d)
                                else
                                    dst(d) = nd + c(d); src_c(d) = c(d) - 1
                                end if
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = 1, ibw_a
                                    ib_wall_src(q)%sf(dst(1), dst(2), dst(3)) = ib_wall_src(q)%sf(src_c(1), src_c(2), src_c(3))
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end do
        end do

    end subroutine s_ib_wall_src_buffers

    !> Species flux residual at a heterogeneous reacting surface.
    subroutine s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_s, d
        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species), W_species(num_species)
        real(wp), intent(out) :: R_species(num_species), omega_s(num_species), mdot_s
        real(wp)              :: mw_IP, mw_s, rho_s, sum_BG
        real(wp)              :: Xs_IP(num_species), Xs_s(num_species)
        real(wp)              :: D_s(num_species), B_s(num_species), G_s(num_species)
        integer               :: k

        call get_mixture_molecular_weight(Ys_IP, mw_IP)
        call get_mixture_molecular_weight(Ys_s, mw_s)
        rho_s = pres*mw_s/(gas_constant*T_s)

        do k = 1, num_species
            Xs_IP(k) = Ys_IP(k)*mw_IP/W_species(k)
            Xs_s(k) = Ys_s(k)*mw_s/W_species(k)
        end do

        call get_species_mass_diffusivities_mixavg(pres, T_s, Ys_s, D_s)
        call get_surface_net_production_rates(rho_s, T_s, Ys_s, omega_s)

        mdot_s = 0._wp
        do k = 1, num_species
            mdot_s = mdot_s + W_species(k)*omega_s(k)
        end do

        sum_BG = 0._wp
        do k = 1, num_species
            B_s(k) = rho_s*D_s(k)*W_species(k)/mw_s
            G_s(k) = (Xs_IP(k) - Xs_s(k))/d
            sum_BG = sum_BG + B_s(k)*G_s(k)
        end do

        do k = 1, num_species
            R_species(k) = -B_s(k)*G_s(k) + Ys_s(k)*(sum_BG + mdot_s) - W_species(k)*omega_s(k)
        end do

    end subroutine s_surface_species_residual

    !> Surface energy residual: gas-side conduction balances heterogeneous reaction heat. Radiation and solid-side conduction are
    !! omitted.
    subroutine s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_s, d
        real(wp), intent(in)  :: Ys_s(num_species)
        real(wp), intent(out) :: R_energy
        real(wp)              :: mw_s, rho_s, k_s, q_rxn

        call get_mixture_molecular_weight(Ys_s, mw_s)
        rho_s = pres*mw_s/(gas_constant*T_s)

        call get_mixture_thermal_conductivity_mixavg(T_s, Ys_s, k_s)
        call get_surface_reaction_heat_flux(rho_s, T_s, Ys_s, q_rxn)

        R_energy = k_s*(T_s - T_IP)/d - q_rxn

    end subroutine s_surface_energy_residual

    !> Assemble the Newton residual for Ns species, with temperature appended only when it is solved. k_bath selects which species
    !! balance the sum constraint displaces; see s_solve_surface for why it must be the most abundant one.
    subroutine s_surface_residual(pres, T_IP, T_s, d, Ys_IP, Ys_s, W_species, k_bath, solve_temperature, flux_scale, &
                                  & energy_scale, R, R_species, omega_s, mdot_s)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_s, d, flux_scale, energy_scale
        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species), W_species(num_species)
        integer, intent(in)   :: k_bath
        logical, intent(in)   :: solve_temperature
        real(wp), intent(out) :: R(num_species + 1), R_species(num_species), omega_s(num_species), mdot_s
        real(wp)              :: R_energy
        integer               :: k

        call s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)

        R = 0._wp
        do k = 1, num_species
            R(k) = R_species(k)/flux_scale
        end do
        R(k_bath) = sum(Ys_s) - 1._wp

        if (solve_temperature) then
            call s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy)
            R(num_species + 1) = R_energy/energy_scale
        end if

    end subroutine s_surface_residual

    !> Ghost state for an immersed surface, from the image-point state and the surface state.
    !!
    !! The natural closure is the linear mirror phi_g = 2*phi_s - phi_IP, which is what makes the midpoint of the ghost/image
    !! pair reproduce the surface value. Its -1 coefficient on phi_IP is unconditional, though, so it drives strictly positive
    !! quantities negative whenever the surface value sits far below the image-point value: a cold wall in hot gas gives a
    !! negative ghost temperature and hence a negative ghost density, and a species the surface consumes faster than half the
    !! free-stream value gives a negative mass fraction. Rather than clamp the result -- which would break sum(Y) = 1 and hide
    !! the excursion -- the ghost state is the largest convex blend of the mirror back toward the surface value that stays
    !! physical:
    !!
    !!     phi_g(theta) = phi_s + theta*(phi_s - phi_IP),   theta in [0, 1]
    !!
    !! theta = 1 is the full second-order mirror; theta = 0 is the first-order Dirichlet ghost phi_g = phi_s, which is
    !! Gibou et al. (JCP 176:205, 2002) Eq. 17 and is what they likewise fall back to where the linear form is ill-behaved.
    !! One theta is shared by every species, so since both endpoints satisfy sum(Y) = 1, so does every blend between them --
    !! exactly, with no renormalization. Temperature gets its own theta: nothing couples it to the composition, and sharing
    !! would let a trace radical the surface consumes to ~1e-9 (O and OH at a burning carbon wall) drag the thermal mirror
    !! down with it, turning a 2100 K ghost into 1205 K over a species whose own excursion was 1e-9. The cost of limiting
    !! is small: the mirror already delivers only first-order wall flux (Ezra et al., Int. J. Heat Mass Transfer, 2025),
    !! and a boundary closure one order below the interior scheme retains the interior convergence rate (Gustafsson,
    !! Math. Comp. 29:396, 1975).
    !!
    !! This bounds the ghost state; it does not make the boundary strictly conservative. Baskaya et al. (Computers & Fluids
    !! 270:106134, 2024) trace the same cold-wall, large-gradient regime in ablation to ghost-cell mass conservation error
    !! that surfaces as spurious blowing, and resolve it only by moving to a flux-based cut-cell boundary -- a different
    !! discretization from this one, not a tuning of it.
    subroutine s_blend_ghost_temperature(T_IP, T_s, T_lo, T_hi, T_g)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: T_IP, T_s, T_lo, T_hi
        real(wp), intent(out) :: T_g
        ! Stop short of the bound so roundoff cannot cross it.
        real(wp), parameter :: blend_safety = 0.9_wp
        real(wp)            :: theta_T

        ! Largest theta <= 1 keeping T_g in [T_lo, T_hi]: the NASA fit range with chemistry, positivity without.
        ! The bounds are compared before dividing, so T_hi = huge() cannot overflow.
        theta_T = 1._wp
        if (T_s < T_IP .and. blend_safety*(T_s - T_lo) < T_IP - T_s) then
            theta_T = blend_safety*(T_s - T_lo)/(T_IP - T_s)
        else if (T_s > T_IP .and. blend_safety*(T_hi - T_s) < T_s - T_IP) then
            theta_T = blend_safety*(T_hi - T_s)/(T_s - T_IP)
        end if
        theta_T = max(theta_T, 0._wp)

        T_g = T_s + theta_T*(T_s - T_IP)

    end subroutine s_blend_ghost_temperature

    !> Species half of the ghost-state blend: one theta shared by every species (see s_blend_ghost_temperature).
    subroutine s_blend_ghost_species(Ys_IP, Ys_s, Ys_g)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species)
        real(wp), intent(out) :: Ys_g(num_species)
        real(wp), parameter   :: blend_safety = 0.9_wp
        real(wp)              :: theta_Y
        integer               :: k

        theta_Y = 1._wp
        do k = 1, num_species
            if (Ys_s(k) < Ys_IP(k)) theta_Y = min(theta_Y, blend_safety*Ys_s(k)/(Ys_IP(k) - Ys_s(k)))
        end do

        do k = 1, num_species
            Ys_g(k) = Ys_s(k) + theta_Y*(Ys_s(k) - Ys_IP(k))
        end do

    end subroutine s_blend_ghost_species

    !> Index of the most abundant species, i.e. the balance the sum constraint displaces.
    !!
    !! Closing sum(Y) = 1 by dropping one species' flux balance is standard (Surface CHEMKIN), but the dropped species absorbs
    !! the roundoff of every other balance, so it must be the bath gas. Dropping a fixed index instead -- the last species in
    !! the mechanism -- lands that error on whatever the mechanism happens to list last, which for the reduced GRI mechanism
    !! shipped with the reacting-surface example is H2O2, a trace radical whose own surface balance then goes unenforced.
    !! Cantera's solveSP re-scans for the largest species each iteration (evalSurfLarge) for this reason.
    subroutine s_pick_bath_species(Ys_s, k_bath)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: Ys_s(num_species)
        integer, intent(out) :: k_bath
        integer              :: k

        k_bath = 1
        do k = 2, num_species
            if (Ys_s(k) > Ys_s(k_bath)) k_bath = k
        end do

    end subroutine s_pick_bath_species

    !> Newton solve for a reacting surface. thermal_bc=0: zero-normal-gradient T; 1: prescribed T; 2: energy balance.
    subroutine s_solve_surface(pres, T_IP, T_wall, d, Ys_IP, W_species, thermal_bc, Ys_s, T_s, mdot_s, converged)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_wall, d
        real(wp), intent(in)  :: Ys_IP(num_species), W_species(num_species)
        integer, intent(in)   :: thermal_bc
        real(wp), intent(out) :: Ys_s(num_species), T_s, mdot_s
        logical, intent(out)  :: converged
        integer, parameter    :: max_iter = 30, max_backtrack = 20
        real(wp), parameter   :: fd_eps_Y = 1.e-7_wp, fd_eps_T = 1.e-6_wp
        ! The Jacobian is a forward difference, so the Newton update carries O(sqrt(epsilon)) noise and the iteration turns
        ! linear near the root. The residual itself is evaluated exactly, so a tighter tolerance is not unreachable, but it is
        ! reached by stalling rather than converging -- which the caller then reports as a failed surface. Cantera's solveSP
        ! runs the same 1e-7 forward-difference step against a 1e-4 relative tolerance; 1e-6 here stays well inside that.
        real(wp), parameter :: tol = 1.e-6_wp
        real(wp), parameter :: Y_tol = 100._wp*epsilon(1._wp)
        real(wp)            :: A(num_species + 1, num_species + 1), rhs(num_species + 1), delta(num_species + 1)
        real(wp)            :: R(num_species + 1), R_pert(num_species + 1), R_trial(num_species + 1)
        real(wp)            :: R_species(num_species), R_species_pert(num_species), R_species_trial(num_species)
        real(wp)            :: omega_s(num_species), omega_pert(num_species), omega_trial(num_species)
        real(wp)            :: Ys_pert(num_species), Ys_trial(num_species)
        real(wp)            :: mdot_pert, mdot_trial, T_pert, T_trial
        real(wp)            :: flux_scale, energy_scale, R_energy, dx, lambda, norm_R, norm_trial
        logical             :: solve_temperature, linear_success, accepted
        integer             :: nsolve, iter, j, iback, k_bath, k_bath_prev

        converged = .false.
        Ys_s = Ys_IP

        select case (thermal_bc)
        case (0)
            T_s = T_IP
            solve_temperature = .false.
        case (1)
            T_s = T_wall
            solve_temperature = .false.
        case (2)
            T_s = min(max(T_IP, T_surface_min), T_surface_max)
            solve_temperature = .true.
        case default
            T_s = T_IP
            omega_s = 0._wp
            mdot_s = 0._wp
            return
        end select

        nsolve = num_species + merge(1, 0, solve_temperature)

        call s_pick_bath_species(Ys_s, k_bath)

        call s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)
        flux_scale = max(maxval(abs(R_species)), 1.e-12_wp)
        energy_scale = 1._wp
        if (solve_temperature) then
            call s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy)
            energy_scale = max(abs(R_energy), 1._wp)
        end if

        call s_surface_residual(pres, T_IP, T_s, d, Ys_IP, Ys_s, W_species, k_bath, solve_temperature, flux_scale, energy_scale, &
                                & R, R_species, omega_s, mdot_s)
        norm_R = maxval(abs(R(1:nsolve)))
        if (norm_R < tol) then
            converged = .true.
            return
        end if

        do iter = 1, max_iter
            do j = 1, num_species
                Ys_pert = Ys_s
                T_pert = T_s
                dx = fd_eps_Y*max(abs(Ys_s(j)), 1._wp)
                if (Ys_s(j) + dx > 1._wp) dx = -dx
                Ys_pert(j) = Ys_pert(j) + dx

                call s_surface_residual(pres, T_IP, T_pert, d, Ys_IP, Ys_pert, W_species, k_bath, solve_temperature, flux_scale, &
                                        & energy_scale, R_pert, R_species_pert, omega_pert, mdot_pert)
                A(1:nsolve,j) = (R_pert(1:nsolve) - R(1:nsolve))/dx
            end do

            if (solve_temperature) then
                Ys_pert = Ys_s
                dx = fd_eps_T*max(abs(T_s), 1._wp)
                T_pert = T_s + dx
                if (T_pert > T_surface_max) then
                    dx = -dx
                    T_pert = T_s + dx
                end if

                call s_surface_residual(pres, T_IP, T_pert, d, Ys_IP, Ys_pert, W_species, k_bath, solve_temperature, flux_scale, &
                                        & energy_scale, R_pert, R_species_pert, omega_pert, mdot_pert)
                A(1:nsolve,nsolve) = (R_pert(1:nsolve) - R(1:nsolve))/dx
            end if

            rhs(1:nsolve) = -R(1:nsolve)
            call s_solve_surface_linear_system(A, rhs, delta, nsolve, linear_success)
            if (.not. linear_success) return

            lambda = 1._wp
            accepted = .false.
            do iback = 1, max_backtrack
                Ys_trial = Ys_s + lambda*delta(1:num_species)
                T_trial = T_s
                if (solve_temperature) T_trial = T_s + lambda*delta(nsolve)

                if (minval(Ys_trial) < -Y_tol .or. maxval(Ys_trial) > 1._wp + Y_tol .or. T_trial < T_surface_min &
                    & .or. T_trial > T_surface_max) then
                    lambda = 0.5_wp*lambda
                    cycle
                end if

                where (Ys_trial < 0._wp) Ys_trial = 0._wp
                where (Ys_trial > 1._wp) Ys_trial = 1._wp

                call s_surface_residual(pres, T_IP, T_trial, d, Ys_IP, Ys_trial, W_species, k_bath, solve_temperature, &
                                        & flux_scale, energy_scale, R_trial, R_species_trial, omega_trial, mdot_trial)
                norm_trial = maxval(abs(R_trial(1:nsolve)))
                if (norm_trial < norm_R) then
                    accepted = .true.
                    exit
                end if
                lambda = 0.5_wp*lambda
            end do

            if (.not. accepted) return

            Ys_s = Ys_trial
            T_s = T_trial
            R = R_trial
            R_species = R_species_trial
            omega_s = omega_trial
            mdot_s = mdot_trial
            norm_R = norm_trial

            ! The composition has moved, so the species the sum constraint displaces may no longer be the most abundant one.
            ! Re-closing changes what R means, so the residual is rebuilt on the rare steps where the choice actually changes.
            k_bath_prev = k_bath
            call s_pick_bath_species(Ys_s, k_bath)
            if (k_bath /= k_bath_prev) then
                call s_surface_residual(pres, T_IP, T_s, d, Ys_IP, Ys_s, W_species, k_bath, solve_temperature, flux_scale, &
                                        & energy_scale, R, R_species, omega_s, mdot_s)
                norm_R = maxval(abs(R(1:nsolve)))
            end if

            if (norm_R < tol) then
                converged = .true.
                return
            end if
        end do

    end subroutine s_solve_surface

    !> Small dense linear solve with partial pivoting for the local surface Newton system.
    subroutine s_solve_surface_linear_system(A, b, x, nsolve, success)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)     :: nsolve
        real(wp), intent(inout) :: A(num_species + 1, num_species + 1)
        real(wp), intent(inout) :: b(num_species + 1)
        real(wp), intent(out)   :: x(num_species + 1)
        logical, intent(out)    :: success
        real(wp)                :: factor, pivot_value, tmp, row_tmp(num_species + 1)
        integer                 :: i, j, k, pivot

        success = .true.
        x = 0._wp

        do k = 1, nsolve - 1
            pivot = k
            pivot_value = abs(A(k, k))
            do i = k + 1, nsolve
                if (abs(A(i, k)) > pivot_value) then
                    pivot = i
                    pivot_value = abs(A(i, k))
                end if
            end do
            if (pivot_value <= epsilon(1._wp)) then
                success = .false.
                return
            end if

            if (pivot /= k) then
                row_tmp(1:nsolve) = A(k,1:nsolve)
                A(k,1:nsolve) = A(pivot,1:nsolve)
                A(pivot,1:nsolve) = row_tmp(1:nsolve)
                tmp = b(k)
                b(k) = b(pivot)
                b(pivot) = tmp
            end if

            do i = k + 1, nsolve
                factor = A(i, k)/A(k, k)
                A(i, k) = 0._wp
                do j = k + 1, nsolve
                    A(i, j) = A(i, j) - factor*A(k, j)
                end do
                b(i) = b(i) - factor*b(k)
            end do
        end do

        if (abs(A(nsolve, nsolve)) <= epsilon(1._wp)) then
            success = .false.
            return
        end if

        x(nsolve) = b(nsolve)/A(nsolve, nsolve)
        do i = nsolve - 1, 1, -1
            tmp = b(i)
            do j = i + 1, nsolve
                tmp = tmp - A(i, j)*x(j)
            end do
            if (abs(A(i, i)) <= epsilon(1._wp)) then
                success = .false.
                return
            end if
            x(i) = tmp/A(i, i)
        end do

    end subroutine s_solve_surface_linear_system

end module m_ibm
