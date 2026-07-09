!>
!! @file
!! @brief Contains module m_ibm

#:include 'macros.fpp'

!> @brief Per-block fine-grid immersed-boundary state for single-body AMR (SP20 static, SP21 prescribed-motion). This is a SEPARATE
!! module from m_ibm on purpose: declaring these non-declare-target derived-type allocatables inside m_ibm corrupts CCE's OpenMP
!! declare-target descriptor table for m_ibm's ghost_points and aborts its @:ALLOCATE with "lib-4425: Uninitialized descriptor for
!! ALLOCATE statement argument" (even for a plain non-AMR IBM run, which never touches this state). A distinct module keeps m_ibm's
!! compiled image identical to the pre-AMR-IB baseline. Do not fold these declarations back into m_ibm.
module m_ibm_fine

    use m_derived_types

    implicit none

    !> Per-block fine IB state: each AMR slot keeps a HOST-side copy of its fine-grid markers field and ghost-point list, computed
    !! from the geometry at fine resolution so the body is resolved on the fine block. Both are COPIED into m_ibm's declare-target
    !! module globals (ib_markers/ghost_points/num_gps) for the fine advance's setup, per-substep moving recompute, and
    !! correct-state, then copied back. markers%sf and gps are host-only parking storage (no device mapping): the declare-target
    !! ib_markers/ ghost_points are allocated once and NEVER reallocated/move_alloc'd/pointer-swapped, so the swap syncs data via
    !! GPU_UPDATE rather than churning the device present table (the detach/attach/move_alloc of a declared array corrupts it on
    !! Cray).
    type ib_fine_state
        type(integer_field)                          :: markers
        type(ghost_point), allocatable, dimension(:) :: gps
        integer                                      :: num_gps
    end type ib_fine_state
    type(ib_fine_state), allocatable :: ib_fine(:)
    integer                          :: num_gps_save

    !> The coarse ghost points AND coarse markers are parked (host copies) across a fine swap in an EXTRA ib_fine slot rather than
    !! dedicated module variables: adding ANY new module-level derived-type allocatable to this module compile-time corrupts a
    !! sibling allocatable's descriptor on CCE OpenMP-offload (the plain-IBM lib-4425 class), so reuse ib_fine's proven-safe
    !! storage. ib_coarse_slot indexes that extra slot (= nslots+1).
    integer :: ib_coarse_slot = 0

    !> Fine-block ghost-point capacity (= buffered fine-block cell count), set by s_ibm_alloc_fine before the coarse s_ibm_setup so
    !! the declare-target ghost_points is sized once to hold the larger of the coarse and fine lists. 0 when AMR-IB is inactive.
    integer(kind=8) :: fine_gps_cap = 0_8

    !> Bounds for the ib_markers marker field, sized once to enclose BOTH the coarse and fine blocks so the declare-target
    !! ib_markers holds either without a reallocation/pointer-swap. Set by s_ibm_alloc_fine, read by s_ibm_setup.
    integer :: mkr_lo(3) = 0, mkr_hi(3) = 0
end module m_ibm_fine

!> @brief Ghost-node immersed boundary method: locates ghost/image points, computes interpolation coefficients, and corrects the
!! flow state
module m_ibm

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_helper
    use m_helper_basic
    use m_constants
    use m_compute_levelset
    use m_ib_patches
    use m_viscous
    use m_model
    use m_patch_geometries
    use m_collisions

    ! Fine-IB AMR state (ib_fine/num_gps_save/mkr bounds) lives in a
    ! separate module so it is NOT part of m_ibm's compiled image: on CCE OpenMP-offload,
    ! declaring those derived-type allocatables here corrupts the declare-target descriptor
    ! for ghost_points and aborts its @:ALLOCATE with lib-4425. See m_ibm_fine.
    use m_ibm_fine

    implicit none

    private :: s_compute_image_points, s_compute_interpolation_coeffs, s_interpolate_image_point, s_find_ghost_points, &
        & s_find_num_ghost_points
    ; public :: ib_gbl_idx_lookup, s_initialize_ibm_module, s_ibm_setup, s_ibm_correct_state, s_finalize_ibm_module, &
        & s_ibm_alloc_fine, s_ibm_setup_fine, s_ibm_swap_to_fine, s_ibm_restore_from_fine, num_gps

    type(integer_field), public :: ib_markers
    $:GPU_DECLARE(create='[ib_markers]')

    type(ghost_point), dimension(:), allocatable :: ghost_points
    $:GPU_DECLARE(create='[ghost_points]')

    integer :: num_gps  !< Number of ghost points
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[gp_layers, num_gps]')
#elif defined(MFC_OpenMP)
    $:GPU_DECLARE(create='[num_gps]')
#endif
    logical :: moving_immersed_boundary_flag

    ! IB MPI buffers
    integer, allocatable  :: send_ids(:), recv_ids(:)
    real(wp), allocatable :: send_ft(:,:), recv_ft(:,:)
    real(wp), allocatable :: recv_forces_snap(:,:), recv_torques_snap(:,:)

contains

    !> Allocates memory for the variables in the IBM module
    impure subroutine s_initialize_ibm_module()

        if (p > 0) then
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, -buff_size:p+buff_size))
        else
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, 0:0))
        end if

        @:ACC_SETUP_SFs(ib_markers)

        $:GPU_ENTER_DATA(copyin='[num_gps]')

        if (collision_model > 0) call s_initialize_collisions_module()

    end subroutine s_initialize_ibm_module

    !> Initializes the values of various IBM variables, such as ghost points and image points.
    impure subroutine s_ibm_setup()

        integer         :: i, j, k
        integer(kind=8) :: max_num_gps

        call nvtxStartRange("SETUP-IBM-MODULE")

        ! GPU routines require updated cell centers
        $:GPU_UPDATE(device='[num_ibs, num_gbl_ibs, x_cc, y_cc, dx, dy, x_domain, y_domain, ib_bc_x%beg, ib_bc_y%beg]')
        if (p /= 0) then
            $:GPU_UPDATE(device='[z_cc, dz, z_domain, ib_bc_z%beg]')
        end if
        $:GPU_UPDATE(device='[patch_ib(1:num_ibs)]')

        ! do all set up for moving immersed boundaries
        $:GPU_PARALLEL_LOOP(private='[i]')
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) then
                call s_compute_moment_of_inertia(patch_ib(i), patch_ib(i)%angular_vel, patch_ib(i)%moment)
            end if
            call s_update_ib_rotation_matrix(i)
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_UPDATE(host='[patch_ib(1:num_ibs)]')

        ! allocate some arrays for MPI communication, if required by this simulation
#ifdef MFC_MPI
        if (num_procs > 1) then
            @:ALLOCATE(send_ids(size(patch_ib)), send_ft(6, size(patch_ib)))
            allocate (recv_forces_snap(size(patch_ib), 3), recv_torques_snap(size(patch_ib), 3), recv_ids(size(patch_ib)), &
                      & recv_ft(6, size(patch_ib)))
        end if
#endif

        call s_update_ib_lookup()

        ! recompute the new ib_patch locations
        ib_markers%sf = 0._wp
        $:GPU_UPDATE(device='[ib_markers%sf]')
        call s_apply_ib_patches(ib_markers)
        $:GPU_UPDATE(host='[ib_markers%sf]')
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) call s_compute_centroid_offset(i)  ! offsets are computed after IB markers are generated
            $:GPU_UPDATE(device='[patch_ib(i)]')
        end do

        ! find the number of ghost points and set them to be the maximum total across ranks
        call s_find_num_ghost_points(num_gps)
        if (moving_immersed_boundary_flag) then
            call s_mpi_allreduce_integer_sum(int(num_gps, 8), max_num_gps)
            max_num_gps = min(max_num_gps*2_8, int(m + 1, 8)*int(n + 1, 8)*int(p + 1, 8))
        else
            max_num_gps = int(num_gps, 8)
        end if
        ! AMR-IB: the fine blocks swap their (larger) ghost-point list into this same declare-target
        ! ghost_points, which is allocated ONCE here and never reallocated (a realloc/move_alloc of a
        ! declare-target allocatable corrupts the Cray present table). Size it to hold the fine list too.
        ! fine_gps_cap is set by s_ibm_alloc_fine, which runs before this routine.
        if (allocated(ib_fine)) max_num_gps = max(max_num_gps, fine_gps_cap)

        ! set the size of the ghost point arrays to be the amount of points total, plus a factor of 2 buffer
        $:GPU_UPDATE(device='[num_gps]')
        ! ghost_points is GPU_DECLARE'd and @:ALLOCATE establishes its device mapping; no
        ! explicit copyin (contents are written by the device pipeline below) - an extra
        ! dynamic map on top of the declared entry corrupts CCE-OMP's descriptor (lib-4425)
        @:ALLOCATE(ghost_points(1:max_num_gps))
        ! Ghost-cell IBM, Tseng & Ferziger JCP (2003), Mittal & Iaccarino ARFM (2005)
        call s_find_ghost_points()
        call s_apply_levelset(ghost_points, num_gps)

        call s_compute_image_points(ghost_points)
        call s_compute_interpolation_coeffs(ghost_points)

        call nvtxEndRange

    end subroutine s_ibm_setup

    !> Update the conservative variables at the ghost points
    subroutine s_ibm_correct_state(q_cons_vf, q_prim_vf, pb_in, mv_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf  !< Primitive Variables
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf  !< Primitive Variables
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), optional, intent(inout) :: pb_in, mv_in
        integer :: i, j, k, l, q, r                                          !< Iterator variables
        integer :: patch_id, patch_id_temp                                   !< Patch ID of ghost point
        real(wp) :: rho, gamma, pi_inf, dyn_pres                             !< Mixture variables
        real(wp), dimension(2) :: Re_K
        real(wp) :: G_K
        real(wp) :: qv_K
        real(wp) :: pres_IP
        real(wp), dimension(3) :: vel_IP, vel_norm_IP
        real(wp) :: c_IP

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)  :: Gs
            real(wp), dimension(3)  :: alpha_rho_IP, alpha_IP
            real(wp), dimension(3)  :: r_IP, v_IP, pb_IP, mv_IP
            real(wp), dimension(18) :: nmom_IP
            real(wp), dimension(12) :: presb_IP, massv_IP
        #:else
            real(wp), dimension(num_fluids) :: Gs
            real(wp), dimension(num_fluids) :: alpha_rho_IP, alpha_IP
            real(wp), dimension(nb)         :: r_IP, v_IP, pb_IP, mv_IP
            real(wp), dimension(nb*nmom)    :: nmom_IP
            real(wp), dimension(nb*nnode)   :: presb_IP, massv_IP
        #:endif
        ! Primitive variables at the image point associated with a ghost point, interpolated from surrounding fluid cells.

        real(wp), dimension(3) :: norm               !< Normal vector from GP to IP
        real(wp), dimension(3) :: physical_loc       !< Physical loc of GP
        real(wp), dimension(3) :: vel_g              !< Velocity of GP
        real(wp), dimension(3) :: radial_vector      !< vector from centroid to ghost point
        real(wp), dimension(3) :: rotation_velocity  !< speed of the ghost point due to rotation
        real(wp)               :: nbub
        real(wp)               :: buf
        type(ghost_point)      :: gp
        type(ghost_point)      :: innerp

        ! set the Moving IBM interior conservative variables
        $:GPU_PARALLEL_LOOP(private='[i, j, k, patch_id, rho]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    patch_id = ib_markers%sf(j, k, l)
                    if (patch_id /= 0) then
                        call s_decode_patch_periodicity(patch_id, patch_id_temp)
                        call s_get_neighborhood_idx(patch_id_temp, patch_id)
                        if (patch_id > 0) then
                            q_prim_vf(eqn_idx%E)%sf(j, k, l) = 1._wp
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

        if (num_gps > 0) then
            $:GPU_PARALLEL_LOOP(private='[i, physical_loc, dyn_pres, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, vel_g, vel_norm_IP, &
                                & r_IP, v_IP, pb_IP, mv_IP, nmom_IP, presb_IP, massv_IP, rho, gamma, pi_inf, Re_K, G_K, Gs, gp, &
                                & innerp, norm, buf, radial_vector, rotation_velocity, j, k, l, q, qv_K, c_IP, nbub, patch_id]')
            do i = 1, num_gps
                gp = ghost_points(i)
                j = gp%loc(1)
                k = gp%loc(2)
                l = gp%loc(3)
                patch_id = ghost_points(i)%ib_patch_id

                ! Calculate physical location of GP
                if (p > 0) then
                    physical_loc = [x_cc(j), y_cc(k), z_cc(l)]
                else
                    physical_loc = [x_cc(j), y_cc(k), 0._wp]
                end if

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
                else
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP)
                end if

                dyn_pres = 0._wp

                ! Set q_prim_vf params at GP so that mixture vars calculated properly
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    q_prim_vf(q)%sf(j, k, l) = alpha_rho_IP(q)
                    q_prim_vf(eqn_idx%adv%beg + q - 1)%sf(j, k, l) = alpha_IP(q)
                end do

                if (surface_tension) then
                    q_prim_vf(eqn_idx%c)%sf(j, k, l) = c_IP
                end if

                ! set the pressure
                if (patch_ib(patch_id)%moving_ibm <= 1) then
                    q_prim_vf(eqn_idx%E)%sf(j, k, l) = pres_IP
                else
                    q_prim_vf(eqn_idx%E)%sf(j, k, l) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, num_fluids
                        ! Pressure correction for moving IB: accounts for acceleration of IB surface
                        q_prim_vf(eqn_idx%E)%sf(j, k, l) = q_prim_vf(eqn_idx%E)%sf(j, k, &
                                  & l) + pres_IP/(1._wp - 2._wp*abs(gp%levelset*alpha_rho_IP(q)/pres_IP) &
                                  & *dot_product(patch_ib(patch_id)%force/patch_ib(patch_id)%mass, gp%levelset_norm))
                    end do
                end if

                if (model_eqns /= model_eqns_4eq) then
                    ! If in simulation, use acc mixture subroutines
                    if (elasticity) then
                        call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv_K, alpha_IP, alpha_rho_IP, Re_K, &
                            & G_K, Gs)
                    else
                        call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv_K, alpha_IP, alpha_rho_IP, Re_K)
                    end if
                end if

                if (patch_ib(patch_id)%moving_ibm /= 0) then
                    ! get the vector that points from the centroid to the ghost
                    radial_vector(1) = physical_loc(1) - (patch_ib(patch_id)%x_centroid + real(ghost_points(i)%x_periodicity, &
                                  & wp)*(x_domain%end - x_domain%beg))
                    radial_vector(2) = physical_loc(2) - (patch_ib(patch_id)%y_centroid + real(ghost_points(i)%y_periodicity, &
                                  & wp)*(y_domain%end - y_domain%beg))
                    radial_vector(3) = 0._wp
                    if (num_dims == 3) radial_vector(3) = physical_loc(3) - (patch_ib(patch_id)%z_centroid &
                        & + real(ghost_points(i)%z_periodicity, wp)*(z_domain%end - z_domain%beg))
                end if

                ! Calculate velocity of ghost cell
                if (gp%slip) then
                    norm(1:3) = gp%levelset_norm
                    buf = sqrt(sum(norm**2))
                    norm = norm/buf
                    vel_norm_IP = sum(vel_IP*norm)*norm
                    vel_g = vel_IP - vel_norm_IP
                    if (patch_ib(patch_id)%moving_ibm /= 0) then
                        ! compute the linear velocity of the ghost point due to rotation
                        call s_cross_product(patch_ib(patch_id)%angular_vel, radial_vector, rotation_velocity)

                        ! add only the component of the IB's motion that is normal to the surface
                        vel_g = vel_g + sum((patch_ib(patch_id)%vel + rotation_velocity)*norm)*norm
                    end if
                else
                    if (patch_ib(patch_id)%moving_ibm == 0) then
                        ! we know the object is not moving if moving_ibm is 0 (false)
                        vel_g = 0._wp
                    else
                        ! convert the angular velocity from the inertial reference frame to the fluids frame, then convert to linear
                        ! velocity
                        call s_cross_product(patch_ib(patch_id)%angular_vel, radial_vector, rotation_velocity)
                        do q = 1, 3
                            ! if mibm is 1 or 2, then the boundary may be moving
                            vel_g(q) = patch_ib(patch_id)%vel(q)  ! add the linear velocity
                            vel_g(q) = vel_g(q) + rotation_velocity(q)  ! add the rotational velocity
                        end do
                    end if
                end if

                ! Set momentum
                $:GPU_LOOP(parallelism='[seq]')
                do q = eqn_idx%mom%beg, eqn_idx%mom%end
                    q_cons_vf(q)%sf(j, k, l) = rho*vel_g(q - eqn_idx%mom%beg + 1)
                    dyn_pres = dyn_pres + q_cons_vf(q)%sf(j, k, l)*vel_g(q - eqn_idx%mom%beg + 1)/2._wp
                end do

                ! Set continuity and adv vars
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    q_cons_vf(q)%sf(j, k, l) = alpha_rho_IP(q)
                    q_cons_vf(eqn_idx%adv%beg + q - 1)%sf(j, k, l) = alpha_IP(q)
                end do

                ! Set color function
                if (surface_tension) then
                    q_cons_vf(eqn_idx%c)%sf(j, k, l) = c_IP
                end if

                ! Set Energy
                if (bubbles_euler) then
                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = (1 - alpha_IP(1))*(gamma*pres_IP + pi_inf + dyn_pres)
                else
                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = gamma*pres_IP + pi_inf + dyn_pres
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
                        q_cons_vf(q)%sf(j, k, &
                                  & l) = alpha_IP(q - eqn_idx%int_en%beg + 1)*(gammas(q - eqn_idx%int_en%beg + 1)*pres_IP &
                                  & + pi_infs(q - eqn_idx%int_en%beg + 1))
                    end do
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_ibm_correct_state

    !> Compute the image points for each ghost point
    impure subroutine s_compute_image_points(ghost_points_in)

        type(ghost_point), dimension(num_gps), intent(inout) :: ghost_points_in
        real(wp)                                             :: dist
        real(wp), dimension(3)                               :: norm
        real(wp), dimension(3)                               :: physical_loc
        real(wp)                                             :: temp_loc
        real(wp), pointer, dimension(:)                      :: s_cc => null()
        integer                                              :: bound
        type(ghost_point)                                    :: gp
        integer                                              :: q, dim      !< Iterator variables
        integer                                              :: i, j, k, l  !< Location indexes
        integer                                              :: patch_id    !< IB Patch ID
        integer                                              :: dir
        integer                                              :: index
        logical                                              :: bounds_error

        bounds_error = .false.

        $:GPU_PARALLEL_LOOP(private='[q, gp, i, j, k, physical_loc, patch_id, dist, norm, dim, bound, dir, index, temp_loc, &
                            & s_cc]', copy='[bounds_error]')
        do q = 1, num_gps
            gp = ghost_points_in(q)
            i = gp%loc(1)
            j = gp%loc(2)
            k = gp%loc(3)

            ! Calculate physical location of ghost point
            if (p > 0) then
                physical_loc = [x_cc(i), y_cc(j), z_cc(k)]
            else
                physical_loc = [x_cc(i), y_cc(j), 0._wp]
            end if

            ! Calculate and store the precise location of the image point
            patch_id = gp%ib_patch_id
            dist = abs(real(gp%levelset, kind=wp))
            norm(:) = gp%levelset_norm
            ghost_points_in(q)%ip_loc(:) = physical_loc(:) + 2*dist*norm(:)

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

                if (f_approx_equal(norm(dim), 0._wp)) then
                    ! if the ghost point is almost equal to a cell location, we set it equal and continue
                    ghost_points_in(q)%ip_grid(dim) = ghost_points_in(q)%loc(dim)
                else
                    if (norm(dim) > 0) then
                        dir = 1
                    else
                        dir = -1
                    end if

                    index = ghost_points_in(q)%loc(dim)
                    temp_loc = ghost_points_in(q)%ip_loc(dim)
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
                            print *, "We are searching in dimension ", dim, " for image point at ", ghost_points_in(q)%ip_loc(:)
                            print *, "Domain size: ", [x_cc(-buff_size), y_cc(-buff_size), z_cc(-buff_size)]
                            print *, "x: ", x_cc(-buff_size), " to: ", x_cc(m + buff_size - 1)
                            print *, "y: ", y_cc(-buff_size), " to: ", y_cc(n + buff_size - 1)
                            if (p /= 0) print *, "z: ", z_cc(-buff_size), " to: ", z_cc(p + buff_size - 1)
                            print *, "Image point is located approximately ", &
                                & (ghost_points_in(q)%loc(dim) - ghost_points_in(q) %ip_loc(dim))/(s_cc(1) - s_cc(0)), &
                                & " grid cells away"
                            print *, "Levelset ", dist, " and Norm: ", norm(:)
                            print *, &
                                & "A short term fix may include increasing buff_size further in m_helper_basic (currently set to a minimum of 10)"
#endif
                            bounds_error = .true.
                        end if
                    end do

                    ghost_points_in(q)%ip_grid(dim) = index
                    if (ghost_points_in(q)%DB(dim) == -1) then
                        ghost_points_in(q)%ip_grid(dim) = ghost_points_in(q)%loc(dim) + 1
                    else if (ghost_points_in(q)%DB(dim) == 1) then
                        ghost_points_in(q)%ip_grid(dim) = ghost_points_in(q)%loc(dim) - 1
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
                            & gp_layers_z]', copyin='[ib_markers%sf]', collapse=3)
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

        ! Operates on the declare-target module-global ghost_points directly (not via a dummy argument):
        ! the on-device parallel loop writes it and the on-device sort below reorders it, both under the
        ! SAME declare-target name and both on the device, so nothing ever crosses host<->device here.
        integer           :: i, j, k, ii, jj, kk, gp_layers_z  !< Iterator variables
        integer           :: xp, yp, zp                        !< periodicities
        integer           :: count, count_i, local_idx
        integer           :: patch_id, encoded_patch_id, neighborhood_patch_id
        logical           :: is_gp
        integer           :: a, b                              !< insertion-sort indices
        logical           :: less                              !< lexicographic comparison result
        type(ghost_point) :: tmp                               !< insertion-sort scratch element

        count = 0
        count_i = 0
        gp_layers_z = gp_layers
        if (p == 0) gp_layers_z = 0

        $:GPU_PARALLEL_LOOP(private='[i, j, k, ii, jj, kk, is_gp, local_idx, patch_id, encoded_patch_id, neighborhood_patch_id, &
                            & xp, yp, zp]', copyin='[count, count_i, x_domain, y_domain, z_domain, ib_markers%sf]', &
                            & firstprivate='[gp_layers, gp_layers_z]', collapse=3)
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

                            ghost_points(local_idx)%loc = [i, j, k]
                            encoded_patch_id = ib_markers%sf(i, j, k)
                            call s_decode_patch_periodicity(encoded_patch_id, patch_id, xp, yp, zp)
                            call s_get_neighborhood_idx(patch_id, neighborhood_patch_id)
                            ghost_points(local_idx)%ib_patch_id = neighborhood_patch_id
                            ghost_points(local_idx)%x_periodicity = xp
                            ghost_points(local_idx)%y_periodicity = yp
                            ghost_points(local_idx)%z_periodicity = zp
                            ghost_points(local_idx)%slip = patch_ib(neighborhood_patch_id)%slip

                            if ((x_cc(i) - dx(i)) < x_domain%beg) then
                                ghost_points(local_idx)%DB(1) = -1
                            else if ((x_cc(i) + dx(i)) > x_domain%end) then
                                ghost_points(local_idx)%DB(1) = 1
                            else
                                ghost_points(local_idx)%DB(1) = 0
                            end if

                            if ((y_cc(j) - dy(j)) < y_domain%beg) then
                                ghost_points(local_idx)%DB(2) = -1
                            else if ((y_cc(j) + dy(j)) > y_domain%end) then
                                ghost_points(local_idx)%DB(2) = 1
                            else
                                ghost_points(local_idx)%DB(2) = 0
                            end if

                            if (p /= 0) then
                                if ((z_cc(k) - dz(k)) < z_domain%beg) then
                                    ghost_points(local_idx)%DB(3) = -1
                                else if ((z_cc(k) + dz(k)) > z_domain%end) then
                                    ghost_points(local_idx)%DB(3) = 1
                                else
                                    ghost_points(local_idx)%DB(3) = 0
                                end if
                            end if
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! The atomic capture above assigns array slots in thread-completion order, so the ghost-point LIST
        ! order is nondeterministic on the GPU and differs from the CPU's serial (loop) order. Any
        ! order-sensitive consumer downstream (e.g. the surface-force reduction) then produces a
        ! backend-dependent result, which the discrete image-point stencil amplifies -> the moving AMR-IB
        ! golden diverges across backends. (Only moving AMR-IB rebuilds the list on-device per substep, so
        ! only it is affected.) Reorder into deterministic lexicographic (i,j,k) order with a single-thread
        ! ON-DEVICE insertion sort (num_gps ~ O(1e2); the single-trip outer loop pins it to one thread).
        ! Sorting on the device is deliberate: a host round-trip here needs a GPU_UPDATE of the declare-target
        ! ghost_points, which fails Cray OpenACC's present-table lookup and aborts CCE OpenMP-offload with
        ! lib-4425 in the AMR fine path (this routine runs mid-swap, see s_ibm_swap_to_fine).
        $:GPU_PARALLEL_LOOP(private='[a, b, tmp, less]')
        do local_idx = 1, 1
            do a = 2, num_gps
                tmp = ghost_points(a)
                b = a - 1
                do
                    if (b < 1) exit
                    less = tmp%loc(1) < ghost_points(b)%loc(1) .or. (tmp%loc(1) == ghost_points(b)%loc(1) .and. (tmp%loc(2) &
                                   & < ghost_points(b)%loc(2) .or. (tmp%loc(2) == ghost_points(b)%loc(2) .and. tmp%loc(3) &
                                   & < ghost_points(b)%loc(3))))
                    if (.not. less) exit
                    ghost_points(b + 1) = ghost_points(b)
                    b = b - 1
                end do
                ghost_points(b + 1) = tmp
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_find_ghost_points

    !> Compute the interpolation coefficients for image points
    subroutine s_compute_interpolation_coeffs(ghost_points_in)

        type(ghost_point), dimension(num_gps), intent(inout) :: ghost_points_in
        real(wp), dimension(2, 2, 2)                         :: dist
        real(wp), dimension(2, 2, 2)                         :: alpha
        real(wp), dimension(2, 2, 2)                         :: interp_coeffs
        real(wp)                                             :: buf
        real(wp), dimension(2, 2, 2)                         :: eta
        type(ghost_point)                                    :: gp
        integer                                              :: q, i, j, k, ii, jj, kk  !< Grid indexes and iterators
        integer                                              :: patch_id
        logical                                              :: is_cell_center

        $:GPU_PARALLEL_LOOP(private='[q, i, j, k, ii, jj, kk, dist, buf, gp, interp_coeffs, eta, alpha, patch_id, &
                            & is_cell_center]', copyin='[ib_markers%sf]')
        do q = 1, num_gps
            gp = ghost_points_in(q)
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
            buf = 1._wp
            do ii = 0, 1
                do jj = 0, 1
                    if (p == 0) then
                        dist(1 + ii, 1 + jj, 1) = sqrt((x_cc(i + ii) - gp%ip_loc(1))**2 + (y_cc(j + jj) - gp%ip_loc(2))**2)
                    else
                        do kk = 0, 1
                            dist(1 + ii, 1 + jj, &
                                 & 1 + kk) = sqrt((x_cc(i + ii) - gp%ip_loc(1))**2 + (y_cc(j + jj) - gp%ip_loc(2))**2 + (z_cc(k &
                                 & + kk) - gp%ip_loc(3))**2)
                        end do
                    end if
                end do
            end do

            ! check if we are arbitrarily close to a cell center
            interp_coeffs = 0._wp
            is_cell_center = .false.
            check_is_cell_center: do ii = 0, 1
                do jj = 0, 1
                    if (dist(ii + 1, jj + 1, 1) <= 1.e-16_wp) then
                        interp_coeffs(ii + 1, jj + 1, 1) = 1._wp
                        is_cell_center = .true.
                        exit check_is_cell_center
                    else
                        if (p /= 0) then
                            if (dist(ii + 1, jj + 1, 2) <= 1.e-16_wp) then
                                interp_coeffs(ii + 1, jj + 1, 2) = 1._wp
                                is_cell_center = .true.
                                exit check_is_cell_center
                            end if
                        end if
                    end if
                end do
            end do check_is_cell_center

            if (.not. is_cell_center) then
                ! if we are not arbitrarily close, interpolate
                alpha = 1._wp
                patch_id = gp%ib_patch_id
                if (ib_markers%sf(i, j, k) /= 0) alpha(1, 1, 1) = 0._wp
                if (ib_markers%sf(i + 1, j, k) /= 0) alpha(2, 1, 1) = 0._wp
                if (ib_markers%sf(i, j + 1, k) /= 0) alpha(1, 2, 1) = 0._wp
                if (ib_markers%sf(i + 1, j + 1, k) /= 0) alpha(2, 2, 1) = 0._wp

                if (p == 0) then
                    eta(:,:,1) = 1._wp/dist(:,:,1)**2
                    buf = sum(alpha(:,:,1)*eta(:,:,1))
                    if (buf > 0._wp) then
                        interp_coeffs(:,:,1) = alpha(:,:,1)*eta(:,:,1)/buf
                    else
                        buf = sum(eta(:,:,1))
                        interp_coeffs(:,:,1) = eta(:,:,1)/buf
                    end if
                else
                    if (ib_markers%sf(i, j, k + 1) /= 0) alpha(1, 1, 2) = 0._wp
                    if (ib_markers%sf(i + 1, j, k + 1) /= 0) alpha(2, 1, 2) = 0._wp
                    if (ib_markers%sf(i, j + 1, k + 1) /= 0) alpha(1, 2, 2) = 0._wp
                    if (ib_markers%sf(i + 1, j + 1, k + 1) /= 0) alpha(2, 2, 2) = 0._wp
                    eta = 1._wp/dist**2
                    buf = sum(alpha*eta)

                    if (buf > 0._wp) then
                        interp_coeffs = alpha*eta/buf
                    else
                        buf = sum(eta)
                        interp_coeffs = eta/buf
                    end if
                end if
            end if

            ghost_points_in(q)%interp_coeffs = interp_coeffs
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_interpolation_coeffs

    !> Interpolate primitive variables to a ghost point's image point using bilinear or trilinear interpolation
    subroutine s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, pb_IP, mv_IP, &
                                         & nmom_IP, pb_in, mv_in, presb_IP, massv_IP)

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
        integer                                         :: i, j, k, l, q           !< Iterator variables
        integer                                         :: i1, i2, j1, j2, k1, k2  !< Iterator variables
        real(wp)                                        :: coeff

        i1 = gp%ip_grid(1); i2 = i1 + 1
        j1 = gp%ip_grid(2); j2 = j1 + 1
        k1 = gp%ip_grid(3); k2 = k1 + 1

        if (p == 0) then
            k1 = 0
            k2 = 0
        end if

        alpha_rho_IP = 0._wp
        alpha_IP = 0._wp
        pres_IP = 0._wp
        vel_IP = 0._wp

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
    impure subroutine s_update_mib(num_ibs, th)

        integer, intent(in) :: num_ibs
        !> AMR subcycling: fine sub-time fraction in [0,1] of the coarse step. When present and >= 0 the moving body is snapshotted
        !! to the linear time interpolation between its coarse t^n position (step_*) and t^{n+1} position (current), matching the
        !! fluid-ghost lerp the subcycle applies, and restored afterwards. Absent/negative => current position.
        real(wp), intent(in), optional :: th
        integer :: i, j, k, z_gp_layers
        logical :: snap
        real(wp) :: sc(3, num_ibs), sa(3, num_ibs)  !< body centroids/angles saved across the sub-time snapshot

        call nvtxStartRange("UPDATE-MIBM")

        snap = .false.
        if (present(th)) then
            if (th >= 0._wp) snap = .true.
        end if
        if (snap) then
            ! The body position/angles were just updated on device by the RK body-motion loop
            ! (m_time_steppers); sync to host before the host-side sub-time interpolation reads them, else the
            ! fine block is built at the stale t^n position on GPU (host-current on CPU, so this only bites GPU).
            $:GPU_UPDATE(host='[patch_ib(1:num_ibs)]')
            do i = 1, num_ibs
                sc(1, i) = patch_ib(i)%x_centroid; sc(2, i) = patch_ib(i)%y_centroid; sc(3, i) = patch_ib(i)%z_centroid
                sa(:,i) = patch_ib(i)%angles
                if (patch_ib(i)%moving_ibm /= 0) then
                    patch_ib(i)%x_centroid = (1._wp - th)*patch_ib(i)%step_x_centroid + th*sc(1, i)
                    patch_ib(i)%y_centroid = (1._wp - th)*patch_ib(i)%step_y_centroid + th*sc(2, i)
                    patch_ib(i)%z_centroid = (1._wp - th)*patch_ib(i)%step_z_centroid + th*sc(3, i)
                    patch_ib(i)%angles = (1._wp - th)*patch_ib(i)%step_angles + th*sa(:,i)
                end if
            end do
            $:GPU_UPDATE(device='[patch_ib(1:num_ibs)]')
        end if

        ! Clears the existing immersed boundary indices
        z_gp_layers = 0; if (p /= 0) z_gp_layers = gp_layers + 1
        $:GPU_PARALLEL_LOOP(private='[i, j, k]')
        do i = -gp_layers - 1, m + gp_layers + 1; do j = -gp_layers - 1, n + gp_layers + 1; do k = -z_gp_layers, p + z_gp_layers
            ib_markers%sf(i, j, k) = 0._wp
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
        ! the ghost_points capacity (a setup-time heuristic) can be outgrown when the moving surface's
        ! discrete cell count increases (body entering the domain, bodies separating, rotating
        ! non-convex geometry); the fill below has no bound check, so overflow would be a silent
        ! device out-of-bounds write. size(ghost_points) is the ACTIVE array's capacity - the coarse
        ! list here, or the fine slot's own (larger) list when the AMR advance has swapped it in.
        @:PROHIBIT(num_gps > size(ghost_points), &
                   & "moving IB: the ghost-point count outgrew the ghost-point array capacity; the body's surface-cell count increased beyond the setup-time sizing")
        call s_find_ghost_points()
        call nvtxEndRange

        call nvtxStartRange("COMPUTE-IMAGE-POINTS")
        call s_apply_levelset(ghost_points, num_gps)
        call s_compute_image_points(ghost_points)
        call s_compute_interpolation_coeffs(ghost_points)
        call nvtxEndRange

        if (snap) then
            do i = 1, num_ibs
                patch_ib(i)%x_centroid = sc(1, i); patch_ib(i)%y_centroid = sc(2, i); patch_ib(i)%z_centroid = sc(3, i)
                patch_ib(i)%angles = sa(:,i)
                if (patch_ib(i)%moving_ibm /= 0) call s_update_ib_rotation_matrix(i)
            end do
            $:GPU_UPDATE(device='[patch_ib(1:num_ibs)]')
        end if

        call nvtxEndRange

    end subroutine s_update_mib

    !> Compute pressure and viscous forces and torques on immersed bodies via volume integration
    subroutine s_compute_ib_forces(q_prim_vf, fluid_pp)

        type(scalar_field), dimension(1:sys_size), intent(in) :: q_prim_vf
        type(physical_parameters), dimension(1:num_fluids), intent(in) :: fluid_pp
        integer :: i, j, k, l, encoded_ib_idx, xp, yp, zp, ib_idx, ib_idx_temp, fluid_idx
        real(wp), dimension(num_ibs, 3) :: forces, torques
        ! viscous stress tensor with temp vectors to hold divergence calculations
        real(wp), dimension(1:3,1:3) :: viscous_stress
        real(wp), dimension(1:3)     :: local_force_contribution, radial_vector, local_torque_contribution
        real(wp)                     :: cell_volume, dynamic_viscosity

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: dynamic_viscosities
        #:else
            real(wp), dimension(num_fluids) :: dynamic_viscosities
        #:endif

        call nvtxStartRange("COMPUTE-IB-FORCES")

        forces = 0._wp
        torques = 0._wp

        if (viscous) then
            do fluid_idx = 1, num_fluids
                if (fluid_pp(fluid_idx)%Re(1) > 0._wp) then
                    dynamic_viscosities(fluid_idx) = 1._wp/fluid_pp(fluid_idx)%Re(1)
                else
                    dynamic_viscosities(fluid_idx) = 0._wp
                end if
            end do
        end if

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, xp, yp, zp, ib_idx, ib_idx_temp, encoded_ib_idx, fluid_idx, radial_vector, &
                            & local_force_contribution, cell_volume, local_torque_contribution, dynamic_viscosity, &
                            & viscous_stress]', copy='[forces, torques]', copyin='[dynamic_viscosities]', collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    encoded_ib_idx = ib_markers%sf(i, j, k)
                    if (encoded_ib_idx /= 0) then
                        call s_decode_patch_periodicity(encoded_ib_idx, ib_idx_temp, xp, yp, zp)
                        call s_get_neighborhood_idx(ib_idx_temp, ib_idx)  ! global patch ID -> local index
                        if (ib_idx > 0) then
                            ! get the vector pointing to the grid cell from the IB centroid
                            radial_vector(1) = x_cc(i) - (patch_ib(ib_idx)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg))
                            radial_vector(2) = y_cc(j) - (patch_ib(ib_idx)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg))
                            radial_vector(3) = 0._wp
                            if (num_dims == 3) radial_vector(3) = z_cc(k) - (patch_ib(ib_idx)%z_centroid + real(zp, &
                                & wp)*(z_domain%end - z_domain%beg))

                            local_force_contribution(:) = 0._wp

                            ! compute the pressure force component, which is the negative pressure gradient
                            do l = -fd_number, fd_number
                                local_force_contribution(1) = local_force_contribution(1) - (fd_coeff_x(l, &
                                                         & i)*q_prim_vf(eqn_idx%E)%sf(i + l, j, k))
                                local_force_contribution(2) = local_force_contribution(2) - (fd_coeff_y(l, &
                                                         & j)*q_prim_vf(eqn_idx%E)%sf(i, j + l, k))
                                if (num_dims == 3) then
                                    local_force_contribution(3) = local_force_contribution(3) - (fd_coeff_z(l, &
                                                             & k)*q_prim_vf(eqn_idx%E)%sf(i, j, k + l))
                                end if
                            end do

                            ! get the viscous stress and add its contribution if that is considered
                            if (viscous) then
                                ! compute the volume-weighted local dynamic viscosity
                                dynamic_viscosity = 0._wp
                                do fluid_idx = 1, num_fluids
                                    ! local dynamic viscosity is the dynamic viscosity of the fluid times alpha of the fluid
                                    dynamic_viscosity = dynamic_viscosity + (q_prim_vf(fluid_idx + eqn_idx%adv%beg - 1)%sf(i, j, &
                                        & k)*dynamic_viscosities(fluid_idx))
                                end do

                                do l = -fd_number, fd_number
                                    call s_compute_viscous_stress_tensor(viscous_stress, q_prim_vf, dynamic_viscosity, i + l, j, k)
                                    local_force_contribution(1:3) = local_force_contribution(1:3) + fd_coeff_x(l, &
                                                             & i)*viscous_stress(1,1:3)

                                    call s_compute_viscous_stress_tensor(viscous_stress, q_prim_vf, dynamic_viscosity, i, j + l, k)
                                    local_force_contribution(1:3) = local_force_contribution(1:3) + fd_coeff_y(l, &
                                                             & j)*viscous_stress(2,1:3)

                                    if (num_dims == 3) then
                                        call s_compute_viscous_stress_tensor(viscous_stress, q_prim_vf, dynamic_viscosity, i, j, &
                                                                             & k + l)
                                        local_force_contribution(1:3) = local_force_contribution(1:3) + fd_coeff_z(l, &
                                                                 & k)*viscous_stress(3,1:3)
                                    end if
                                end do
                            end if

                            call s_cross_product(radial_vector, local_force_contribution, local_torque_contribution)

                            ! Update the force and torque values atomically to prevent race conditions
                            cell_volume = dx(i)*dy(j)
                            if (num_dims == 3) cell_volume = cell_volume*dz(k)
                            do l = 1, num_dims
                                $:GPU_ATOMIC(atomic='update')
                                forces(ib_idx, l) = forces(ib_idx, l) + (local_force_contribution(l)*cell_volume)
                                $:GPU_ATOMIC(atomic='update')
                                torques(ib_idx, l) = torques(ib_idx, l) + local_torque_contribution(l)*cell_volume
                            end do
                        end if  ! ib_idx > 0
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_apply_collision_forces(ghost_points, num_gps, ib_markers, forces, torques)

        ! reduce the forces across local neighborhood ranks
        call s_communicate_ib_forces(forces, torques)

        ! consider body forces after reducing to avoid double counting
        do i = 1, num_ibs
            if (bf_x) then
                forces(i, 1) = forces(i, 1) + accel_bf(1)*patch_ib(i)%mass
            end if
            if (bf_y) then
                forces(i, 2) = forces(i, 2) + accel_bf(2)*patch_ib(i)%mass
            end if
            if (bf_z) then
                forces(i, 3) = forces(i, 3) + accel_bf(3)*patch_ib(i)%mass
            end if
        end do

        ! apply the summed forces
        $:GPU_PARALLEL_LOOP(private='[i]', copyin='[forces, torques]')
        do i = 1, num_ibs
            patch_ib(i)%force(:) = forces(i,:)
            patch_ib(i)%torque(:) = torques(i,:)
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_compute_ib_forces

    !> Computes the center of mass for IB patch types where we are unable to determine their center of mass analytically.
    !> These patches include things like NACA airfoils and STL models
    subroutine s_compute_centroid_offset(ib_marker)

        integer, intent(in)      :: ib_marker
        integer                  :: i, j, k, num_cells_local, decoded_gbl_id
        integer(kind=8)          :: num_cells
        real(wp), dimension(1:3) :: center_of_mass, center_of_mass_local

        ! Offset only needs to be computes for specific geometries

        if (patch_ib(ib_marker)%geometry == 4 .or. patch_ib(ib_marker)%geometry == 5 .or. patch_ib(ib_marker)%geometry == 11 &
            & .or. patch_ib(ib_marker)%geometry == 12) then
            center_of_mass_local = [0._wp, 0._wp, 0._wp]
            num_cells_local = 0

            ! get the summed mass distribution and number of cells to divide by
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) /= 0) then
                            call s_decode_patch_periodicity(ib_markers%sf(i, j, k), decoded_gbl_id)
                            if (decoded_gbl_id == patch_ib(ib_marker)%gbl_patch_id) then
                                num_cells_local = num_cells_local + 1
                                center_of_mass_local = center_of_mass_local + [x_cc(i), y_cc(j), 0._wp]
                                if (num_dims == 3) center_of_mass_local(3) = center_of_mass_local(3) + z_cc(k)
                            end if
                        end if
                    end do
                end do
            end do

            ! reduce the mass contribution over all MPI ranks and compute COM
            call s_mpi_allreduce_integer_sum(int(num_cells_local, 8), num_cells)
            if (num_cells /= 0) then
                call s_mpi_allreduce_sum(center_of_mass_local(1), center_of_mass(1))
                call s_mpi_allreduce_sum(center_of_mass_local(2), center_of_mass(2))
                call s_mpi_allreduce_sum(center_of_mass_local(3), center_of_mass(3))
                center_of_mass = center_of_mass/real(num_cells, wp)
            else
                patch_ib(ib_marker)%centroid_offset = [0._wp, 0._wp, 0._wp]
                return
            end if

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
        else
            patch_ib(ib_marker)%centroid_offset(:) = [0._wp, 0._wp, 0._wp]
        end if

    end subroutine s_compute_centroid_offset

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
            #:for X, ID in [('x', 1), ('y', 2), ('z', 3)]
                if (num_dims >= ${ID}$) then
                    ! check for periodicity
                    if (ib_bc_${X}$%beg == BC_PERIODIC) then
                        ! check if the boundary has left the domain, and then correct
                        if (patch_ib(patch_id)%${X}$_centroid < ${X}$_domain%beg) then
                            ! if the boundary exited "left", wrap it back around to the "right"
                            patch_ib(patch_id)%${X}$_centroid = patch_ib(patch_id)%${X}$_centroid + (${X}$_domain%end &
                                     & - ${X}$_domain%beg)
                        else if (patch_ib(patch_id)%${X}$_centroid > ${X}$_domain%end) then
                            ! if the boundary exited "right", wrap it back around to the "left"
                            patch_ib(patch_id)%${X}$_centroid = patch_ib(patch_id)%${X}$_centroid - (${X}$_domain%end &
                                     & - ${X}$_domain%beg)
                        end if
                    end if
                end if
            #:endfor
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_wrap_periodic_ibs

    !> @brief Swaps ownership of IBs and passes ownership of IBs to neighbor processors
    !> Reduces forces and torques across the local neighborhood without a global allreduce. Accumulation phase: 2 passes per
    !! dimension receiving from the low-index (-X) neighbor. Pass 1: add received values; save what was received as recv_snap. Pass
    !! 2: send current (post-pass-1) values; add received; subtract recv_snap to remove double-counting of the direct contribution
    !! already added in pass 1. Back-propagation phase: 2 passes per dimension receiving from the high-index (+X) neighbor, each
    !! overwriting local forces with the neighbor's accumulated total.
    subroutine s_communicate_ib_forces(forces, torques)

        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques

#ifdef MFC_MPI
        integer                       :: i, j, k, pack_pos, unpack_pos, buf_size, ierr
        integer                       :: send_neighbor, recv_neighbor, recv_count, tag
        character(len=1), allocatable :: ib_force_send_buf(:), ib_force_recv_buf(:)

        if (num_procs == 1) return

        buf_size = storage_size(0)/8 + (storage_size(0)/8 + 6*storage_size(0._wp)/8)*size(patch_ib)
        allocate (ib_force_send_buf(buf_size), ib_force_recv_buf(buf_size))

        ! Accumulation phase: propagate contributions toward the high-index corner.
        #:for X, ID in [('x', 1), ('y', 2), ('z', 3)]
            if (num_dims >= ${ID}$) then
                send_neighbor = merge(bc_${X}$%end, MPI_PROC_NULL, bc_${X}$%end >= 0)
                recv_neighbor = merge(bc_${X}$%beg, MPI_PROC_NULL, bc_${X}$%beg >= 0)

                recv_forces_snap = 0._wp
                recv_torques_snap = 0._wp
                tag = 300

                do k = 1, min(2*ib_neighborhood_radius, num_procs_${X}$ - 1)
                    ! send forces to +${X}$ neighbor; receive from -${X}$ neighbor. Add received values then
                    pack_pos = 0
                    $:GPU_PARALLEL_LOOP(private='[i]', copyin='[forces, torques]')
                    do i = 1, num_ibs
                        send_ids(i) = patch_ib(i)%gbl_patch_id
                        send_ft(1:3,i) = forces(i,:)
                        send_ft(4:6,i) = torques(i,:)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_UPDATE(host='[send_ids, send_ft]')
                    call MPI_PACK(num_ibs, 1, MPI_INTEGER, ib_force_send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)
                    call MPI_PACK(send_ids, num_ibs, MPI_INTEGER, ib_force_send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)
                    call MPI_PACK(send_ft, 6*num_ibs, mpi_p, ib_force_send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)
                    call MPI_SENDRECV(ib_force_send_buf, pack_pos, MPI_PACKED, send_neighbor, tag, ib_force_recv_buf, buf_size, &
                                      & MPI_PACKED, recv_neighbor, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    if (recv_neighbor /= MPI_PROC_NULL) then
                        unpack_pos = 0
                        call MPI_UNPACK(ib_force_recv_buf, buf_size, unpack_pos, recv_count, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                        call MPI_UNPACK(ib_force_recv_buf, buf_size, unpack_pos, recv_ids, recv_count, MPI_INTEGER, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_UNPACK(ib_force_recv_buf, buf_size, unpack_pos, recv_ft, 6*recv_count, mpi_p, MPI_COMM_WORLD, ierr)
                        $:GPU_PARALLEL_LOOP(private='[i, j]', copyin='[recv_ft, recv_ids]', copy='[forces, torques, &
                                            & recv_forces_snap, recv_torques_snap]')
                        do i = 1, recv_count
                            call s_get_neighborhood_idx(recv_ids(i), j)
                            if (j > 0) then
                                ! add forces and subtract recv_snap prevent double-counting
                                forces(j,:) = forces(j,:) + recv_ft(1:3,i) - recv_forces_snap(j,:)
                                torques(j,:) = torques(j,:) + recv_ft(4:6,i) - recv_torques_snap(j,:)
                                recv_forces_snap(j,:) = recv_ft(1:3,i)
                                recv_torques_snap(j,:) = recv_ft(4:6,i)
                            end if
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                    tag = tag + 2
                end do
            end if
        #:endfor

        ! Send final sums back to neighbors in -X direction
        #:for X, ID in [('x', 1), ('y', 2), ('z', 3)]
            if (num_dims >= ${ID}$) then
                send_neighbor = merge(bc_${X}$%beg, MPI_PROC_NULL, bc_${X}$%beg >= 0)
                recv_neighbor = merge(bc_${X}$%end, MPI_PROC_NULL, bc_${X}$%end >= 0)

                do k = 1, min(2*ib_neighborhood_radius, num_procs_${X}$ - 1)
                    pack_pos = 0
                    $:GPU_PARALLEL_LOOP(private='[i]', copyin='[forces, torques]')
                    do i = 1, num_ibs
                        send_ids(i) = patch_ib(i)%gbl_patch_id
                        send_ft(1:3,i) = forces(i,:)
                        send_ft(4:6,i) = torques(i,:)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_UPDATE(host='[send_ids, send_ft]')
                    call MPI_PACK(num_ibs, 1, MPI_INTEGER, ib_force_send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)
                    call MPI_PACK(send_ids, num_ibs, MPI_INTEGER, ib_force_send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)
                    call MPI_PACK(send_ft, 6*num_ibs, mpi_p, ib_force_send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)
                    call MPI_SENDRECV(ib_force_send_buf, pack_pos, MPI_PACKED, send_neighbor, tag, ib_force_recv_buf, buf_size, &
                                      & MPI_PACKED, recv_neighbor, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                    if (recv_neighbor /= MPI_PROC_NULL) then
                        unpack_pos = 0
                        call MPI_UNPACK(ib_force_recv_buf, buf_size, unpack_pos, recv_count, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                        call MPI_UNPACK(ib_force_recv_buf, buf_size, unpack_pos, recv_ids, recv_count, MPI_INTEGER, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_UNPACK(ib_force_recv_buf, buf_size, unpack_pos, recv_ft, 6*recv_count, mpi_p, MPI_COMM_WORLD, ierr)
                        $:GPU_PARALLEL_LOOP(private='[i, j]', copyin='[recv_ft, recv_ids]', copy='[forces, torques]')
                        do i = 1, recv_count
                            call s_get_neighborhood_idx(recv_ids(i), j)
                            if (j > 0) then
                                forces(j,:) = recv_ft(1:3,i)
                                torques(j,:) = recv_ft(4:6,i)
                            end if
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                    tag = tag + 2
                end do
            end if
        #:endfor
#endif

    end subroutine s_communicate_ib_forces

    subroutine s_handoff_ib_ownership()

        integer                               :: i, j, k, output_idx, local_output_idx
        integer                               :: old_num_local_ibs
        integer                               :: new_count, recv_count
        integer                               :: pack_pos, unpack_pos, buf_size, patch_bytes
        integer                               :: send_neighbor, recv_neighbor, ierr
        integer                               :: dx, dy, dz, tag, nbr_idx, nreqs
        real(wp), dimension(3)                :: centroid
        logical                               :: is_new
        type(ib_patch_parameters)             :: tmp_patch
        integer, dimension(num_local_ibs_max) :: local_ib_idx_old
        ! 26 neighbors max in 3D (8 in 2D); each gets its own recv buffer
        integer, parameter             :: max_nbrs = 26
        character(len=1), allocatable  :: send_buf(:), recv_bufs(:,:)
        integer, dimension(2*max_nbrs) :: requests
        integer, dimension(max_nbrs)   :: recv_neighbor_list

#ifdef MFC_MPI
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
            output_idx = 0
            local_output_idx = 0
            do i = 1, num_ibs
                centroid = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, 0._wp]
                if (num_dims == 3) centroid(3) = patch_ib(i)%z_centroid

                ! delete if not in neighborhood
                if (f_neighborhood_ranks_own_location(centroid)) then
                    output_idx = output_idx + 1
                    if (i /= output_idx) then
                        patch_ib(output_idx) = patch_ib(i)
                    end if

                    ! check if in local domain
                    if (f_local_rank_owns_location(centroid)) then
                        local_output_idx = local_output_idx + 1
                        local_ib_patch_ids(local_output_idx) = output_idx
                    end if
                end if
            end do
            num_ibs = output_idx
            num_local_ibs = local_output_idx
            $:GPU_UPDATE(device='[patch_ib]')
            call s_update_ib_lookup()

            ! Broadcast newly-owned patches to all neighborhood neighbors
            patch_bytes = storage_size(tmp_patch)/8
            buf_size = storage_size(0)/8 + patch_bytes*num_local_ibs_max
            allocate (send_buf(buf_size), recv_bufs(buf_size, max_nbrs))

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
            do dz = merge(-1, 0, num_dims == 3), merge(1, 0, num_dims == 3)
                do dy = -1, 1
                    do dx = -1, 1
                        if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
                        nbr_idx = nbr_idx + 1
                        tag = 200 + (dx + 1)*9 + (dy + 1)*3 + (dz + 1)
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

            do dz = merge(-1, 0, num_dims == 3), merge(1, 0, num_dims == 3)
                do dy = -1, 1
                    do dx = -1, 1
                        if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
                        tag = 200 + (dx + 1)*9 + (dy + 1)*3 + (dz + 1)
                        send_neighbor = ib_neighbor_ranks(dx, dy, dz)
                        if (send_neighbor < 0) cycle
                        nreqs = nreqs + 1
                        call MPI_ISEND(send_buf, pack_pos, MPI_PACKED, send_neighbor, tag, MPI_COMM_WORLD, requests(nreqs), ierr)
                    end do
                end do
            end do

            call MPI_WAITALL(nreqs, requests, MPI_STATUSES_IGNORE, ierr)

            ! Unpack all received buffers
            do nbr_idx = 1, merge(26, 8, num_dims == 3)
                if (recv_neighbor_list(nbr_idx) == MPI_PROC_NULL) cycle
                unpack_pos = 0
                call MPI_UNPACK(recv_bufs(:,nbr_idx), buf_size, unpack_pos, recv_count, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                do i = 1, recv_count
                    call MPI_UNPACK(recv_bufs(:,nbr_idx), buf_size, unpack_pos, tmp_patch, patch_bytes, MPI_BYTE, MPI_COMM_WORLD, &
                                    & ierr)
                    call s_get_neighborhood_idx(tmp_patch%gbl_patch_id, j)
                    if (j < 0) then
                        num_ibs = num_ibs + 1
                        @:ASSERT(num_ibs <= size(patch_ib), 'patch_ib overflow in neighborhood handoff')
                        patch_ib(num_ibs) = tmp_patch
                    end if
                end do
            end do

            deallocate (send_buf, recv_bufs)
            $:GPU_UPDATE(device='[patch_ib]')
            call s_update_ib_lookup()
        end if
#endif

    end subroutine s_handoff_ib_ownership

    subroutine s_get_neighborhood_idx(gbl_idx, neighborhood_idx)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)  :: gbl_idx
        integer, intent(out) :: neighborhood_idx
        integer              :: i

        neighborhood_idx = ib_gbl_idx_lookup(gbl_idx)

    end subroutine s_get_neighborhood_idx

    subroutine s_update_ib_lookup()

        integer :: i

        ib_gbl_idx_lookup = -1
        $:GPU_UPDATE(device='[ib_gbl_idx_lookup]')

        $:GPU_PARALLEL_LOOP(private='[i]')
        do i = 1, num_ibs
            ib_gbl_idx_lookup(patch_ib(i)%gbl_patch_id) = i
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_UPDATE(host='[ib_gbl_idx_lookup]')

    end subroutine s_update_ib_lookup

    !> Allocate the per-slot fine-IB marker fields (static-body AMR). One integer field per AMR slot, sized to the max buffered fine
    !! extents (mirrors the coarse ib_markers bounds); ghost-point lists start empty and are filled by s_ibm_setup_fine. No-op
    !! unless amr .and. ib.
    impure subroutine s_ibm_alloc_fine(nslots, f1_lo, f1_hi, f2_lo, f2_hi, f3_lo, f3_hi)

        integer, intent(in) :: nslots, f1_lo, f1_hi, f2_lo, f2_hi, f3_lo, f3_hi
        integer             :: islot

        ! ghost-point capacity for any fine block = its buffered cell count (a block can hold no more
        ! ghost points than cells). s_ibm_setup reads this to size the shared declare-target ghost_points.

        fine_gps_cap = int(f1_hi - f1_lo + 1, 8)*int(f2_hi - f2_lo + 1, 8)*int(f3_hi - f3_lo + 1, 8)

        ! Marker-field bounds = the coarse ib_markers bounds (m,n,p are the coarse grid here - alloc_fine runs
        ! before any fine swap, and ib_markers was already allocated to these in s_initialize_ibm_module). The
        ! host park copies must match ib_markers for the whole-array swap copy. The fine block's local index
        ! range must fit inside these (it does for a body smaller than the coarse block); guarded below.
        mkr_lo(1) = -buff_size; mkr_hi(1) = m + buff_size
        mkr_lo(2) = -buff_size; mkr_hi(2) = n + buff_size
        if (p > 0) then
            mkr_lo(3) = -buff_size; mkr_hi(3) = p + buff_size
        else
            mkr_lo(3) = 0; mkr_hi(3) = 0
        end if
        @:PROHIBIT(f1_hi > mkr_hi(1) .or. f2_hi > mkr_hi(2) .or. (p > 0 .and. f3_hi > mkr_hi(3)), &
                   & "AMR fine IB: fine block extent exceeds the coarse ib_markers bounds; the copy-based fine-marker swap needs ib_markers sized to enclose the fine block")

        ! Extra slot (nslots+1) parks the coarse ghost points AND coarse markers during a fine swap - reusing
        ! an ib_fine slot avoids adding a new module-level derived-type allocatable (which corrupts descriptors
        ! on CCE OpenMP, lib-4425). markers%sf and gps are HOST-only park storage (no ACC_SETUP / device map);
        ! the declare-target ib_markers/ghost_points hold the active data and the swap copies to/from these.
        ib_coarse_slot = nslots + 1
        allocate (ib_fine(1:ib_coarse_slot))
        do islot = 1, ib_coarse_slot
            allocate (ib_fine(islot)%markers%sf(mkr_lo(1):mkr_hi(1),mkr_lo(2):mkr_hi(2),mkr_lo(3):mkr_hi(3)))
            ib_fine(islot)%markers%sf = 0
            ib_fine(islot)%num_gps = 0
            allocate (ib_fine(islot)%gps(1:fine_gps_cap))
        end do

    end subroutine s_ibm_alloc_fine

    !> Swap the module IB globals (ib_markers/ghost_points/num_gps) to fine slot islot's stored state; the coarse state parks in the
    !! save slot (host copies of markers and ghost points). MUST be paired with s_ibm_restore_from_fine. Grid globals must already
    !! be swapped to the fine block.
    impure subroutine s_ibm_swap_to_fine(islot, gps_on_device)

        integer, intent(in) :: islot
        logical, intent(in) :: gps_on_device  !< fine ghost points already present on device (per-stage correct path)
        integer             :: n_c

        ! ib_markers: declare-target field, allocated once and device-resident; NEVER pointer-swapped or
        ! detach-attach'd (that corrupts the Cray present table and leaves the device descriptor pointing at
        ! the stale coarse array). Park the coarse markers as a host copy, then on the correct/moving path
        ! copy this slot's fine markers in and push to device. On the setup path s_ibm_setup_fine rebuilds the
        ! fine markers directly in ib_markers.

        $:GPU_UPDATE(host='[ib_markers%sf]')
        ib_fine(ib_coarse_slot)%markers%sf = ib_markers%sf
        if (gps_on_device) then
            ib_markers%sf = ib_fine(islot)%markers%sf
            $:GPU_UPDATE(device='[ib_markers%sf]')
        end if

        ! ghost_points: the declare-target array is allocated once (s_ibm_setup) and stays device-resident;
        ! it is NEVER move_alloc'd/reallocated/detach-attach'd (that corrupts the Cray present table). Park
        ! the coarse list as a host copy, then on the correct/moving path copy this slot's fine list in and
        ! push it to device. On the setup path the fine list does not exist yet - s_ibm_setup_fine fills
        ! ghost_points in place next.
        ! Update only the active 1:num_gps slice (all the host code below touches): amdflang lowers a
        ! whole-array update of this array-of-derived-type to a per-element custom mapper that exhausts
        ! the ROCm HSA resource pool (spurious OUT_OF_RESOURCES abort); elements past num_gps are stale.
        $:GPU_UPDATE(host='[ghost_points(1:num_gps)]')
        n_c = num_gps
        @:PROHIBIT(int(n_c, 8) > size(ib_fine(ib_coarse_slot)%gps, kind=8), &
                   & "AMR fine IB: coarse ghost-point count exceeds the coarse park capacity")
        ib_fine(ib_coarse_slot)%gps(1:n_c) = ghost_points(1:n_c)
        num_gps_save = num_gps
        num_gps = ib_fine(islot)%num_gps
        if (gps_on_device) then
            ghost_points(1:num_gps) = ib_fine(islot)%gps(1:num_gps)
            $:GPU_UPDATE(device='[ghost_points(1:num_gps)]')
        end if
        $:GPU_UPDATE(device='[num_gps]')

    end subroutine s_ibm_swap_to_fine

    !> Restore the coarse IB globals saved by s_ibm_swap_to_fine, parking the (possibly updated) fine state back in slot islot.
    impure subroutine s_ibm_restore_from_fine(islot)

        integer, intent(in) :: islot

        ! Mirror s_ibm_swap_to_fine: save this slot's (freshly computed / motion-updated) fine markers and
        ! ghost points back to its host store, then copy the parked coarse markers and ghost points back into
        ! the device-resident ib_markers/ghost_points. No pointer-swap/detach/move_alloc of the declared arrays.

        $:GPU_UPDATE(host='[ib_markers%sf]')
        ib_fine(islot)%markers%sf = ib_markers%sf
        ib_markers%sf = ib_fine(ib_coarse_slot)%markers%sf
        $:GPU_UPDATE(device='[ib_markers%sf]')

        $:GPU_UPDATE(host='[ghost_points(1:num_gps)]')
        ib_fine(islot)%gps(1:num_gps) = ghost_points(1:num_gps)
        ib_fine(islot)%num_gps = num_gps
        num_gps = num_gps_save
        ghost_points(1:num_gps) = ib_fine(ib_coarse_slot)%gps(1:num_gps)
        $:GPU_UPDATE(device='[ghost_points(1:num_gps)]')
        $:GPU_UPDATE(device='[num_gps]')

    end subroutine s_ibm_restore_from_fine

    !> Compute the fine-grid IB state (markers, ghost points, levelset, image points, interpolation coeffs) for the current block.
    !! The grid globals must be swapped to the fine block AND the IB globals swapped to this slot (s_ibm_swap_to_fine) before the
    !! call: the pipeline writes into the module globals, which then hold this slot's fine state. Mirrors the static-body portion of
    !! s_ibm_setup at fine resolution.
    impure subroutine s_ibm_setup_fine()

        ib_markers%sf = 0
        $:GPU_UPDATE(device='[ib_markers%sf]')
        call s_apply_ib_patches(ib_markers)
        $:GPU_UPDATE(host='[ib_markers%sf]')

        call s_find_num_ghost_points(num_gps)
        $:GPU_UPDATE(device='[num_gps]')
        ! ghost_points is allocated once (s_ibm_setup, sized to the fine-block cell-count cap) and filled
        ! in place; it is never reallocated here (a realloc/move_alloc of the declare-target array corrupts
        ! the Cray present table). The cap bounds any block's ghost-point count, so this cannot overflow.
        @:PROHIBIT(int(num_gps, 8) > size(ghost_points, kind=8), &
                   & "AMR fine IB: ghost-point count exceeds the ghost_points capacity sized at s_ibm_setup")

        call s_find_ghost_points()
        call s_apply_levelset(ghost_points, num_gps)
        call s_compute_image_points(ghost_points)
        call s_compute_interpolation_coeffs(ghost_points)

        if (num_gps > 0) print '(A,I0,A,I0)', ' [amr] block ', amr_cur, ': fine IB ghost points = ', num_gps

    end subroutine s_ibm_setup_fine

    !> Finalize the IBM module
    impure subroutine s_finalize_ibm_module()

        integer :: i

        @:DEALLOCATE(ib_markers%sf)
        @:DEALLOCATE(ib_gbl_idx_lookup)
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
        end if
        if (allocated(ib_fine)) then
            do i = 1, size(ib_fine)
                ! markers%sf and gps are host-only parking storage (no device mapping) - plain deallocate
                if (associated(ib_fine(i)%markers%sf)) then
                    deallocate (ib_fine(i)%markers%sf)
                end if
                if (allocated(ib_fine(i)%gps)) then
                    deallocate (ib_fine(i)%gps)
                end if
            end do
            deallocate (ib_fine)
        end if
        if (collision_model > 0) call s_finalize_collisions_module()
#ifdef MFC_MPI
        if (num_procs > 1) then
            @:DEALLOCATE(send_ids, send_ft)
            deallocate (recv_forces_snap, recv_torques_snap, recv_ids, recv_ft)
        end if
#endif

    end subroutine s_finalize_ibm_module

end module m_ibm
