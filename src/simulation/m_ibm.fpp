!>
!! @file m_ibm.fpp
!! @brief Contains module m_ibm

#:include 'macros.fpp'

!> @brief This module is used to handle all operations related to immersed
!!              boundary methods (IBMs)
module m_ibm

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_constants

    use m_compute_levelset

    use m_ib_patches

    use m_viscous

    implicit none

    private :: s_compute_image_points, &
               s_compute_interpolation_coeffs, &
               s_interpolate_image_point, &
               s_find_ghost_points, &
               s_find_num_ghost_points
    ; public :: s_initialize_ibm_module, &
 s_ibm_setup, &
 s_ibm_correct_state, &
 s_finalize_ibm_module

    integer, allocatable, dimension(:, :, :) :: patch_id_fp
    type(integer_field), public :: ib_markers
    type(levelset_field), public :: levelset
    type(levelset_norm_field), public :: levelset_norm
    $:GPU_DECLARE(create='[ib_markers,levelset,levelset_norm]')

    type(ghost_point), dimension(:), allocatable :: ghost_points
    type(ghost_point), dimension(:), allocatable :: inner_points
    $:GPU_DECLARE(create='[ghost_points,inner_points]')

    integer :: num_gps !< Number of ghost points
    integer :: num_inner_gps !< Number of ghost points
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[gp_layers,num_gps,num_inner_gps]')
#elif defined(MFC_OpenMP)
    $:GPU_DECLARE(create='[num_gps,num_inner_gps]')
#endif
    logical :: moving_immersed_boundary_flag

contains

    !>  Allocates memory for the variables in the IBM module
    impure subroutine s_initialize_ibm_module()

        if (p > 0) then
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, &
                -buff_size:n+buff_size, -buff_size:p+buff_size))
            @:ALLOCATE(levelset%sf(-buff_size:m+buff_size, &
                -buff_size:n+buff_size, -buff_size:p+buff_size, 1:num_ibs))
            @:ALLOCATE(levelset_norm%sf(-buff_size:m+buff_size, &
                -buff_size:n+buff_size, -buff_size:p+buff_size, 1:num_ibs, 1:3))
        else
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, &
                -buff_size:n+buff_size, 0:0))
            @:ALLOCATE(levelset%sf(-buff_size:m+buff_size, &
                -buff_size:n+buff_size, 0:0, 1:num_ibs))
            @:ALLOCATE(levelset_norm%sf(-buff_size:m+buff_size, &
                -buff_size:n+buff_size, 0:0, 1:num_ibs, 1:3))
        end if

        @:ACC_SETUP_SFs(ib_markers)
        @:ACC_SETUP_SFs(levelset)
        @:ACC_SETUP_SFs(levelset_norm)

        $:GPU_ENTER_DATA(copyin='[num_gps,num_inner_gps]')

    end subroutine s_initialize_ibm_module

    !> Initializes the values of various IBM variables, such as ghost points and
    !! image points.
    impure subroutine s_ibm_setup()

        integer :: i, j, k
        integer :: max_num_gps, max_num_inner_gps

        ! do all set up for moving immersed boundaries
        moving_immersed_boundary_flag = .false.
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) then
                call s_compute_moment_of_inertia(i, patch_ib(i)%angular_vel)
                moving_immersed_boundary_flag = .true.
            end if
            call s_update_ib_rotation_matrix(i)
            call s_compute_centroid_offset(i)
        end do
        $:GPU_ENTER_DATA(copyin='[patch_ib]')

        ! Allocating the patch identities bookkeeping variable
        allocate (patch_id_fp(0:m, 0:n, 0:p))

        $:GPU_UPDATE(device='[ib_markers%sf]')
        $:GPU_UPDATE(device='[levelset%sf]')
        $:GPU_UPDATE(device='[levelset_norm%sf]')

        ! Get neighboring IB variables from other processors
        call s_populate_ib_buffers()

        $:GPU_UPDATE(host='[ib_markers%sf]')

        ! find the number of ghost points and set them to be the maximum total across ranks
        call s_find_num_ghost_points(num_gps, num_inner_gps)
        call s_mpi_allreduce_integer_sum(num_gps, max_num_gps)
        call s_mpi_allreduce_integer_sum(num_inner_gps, max_num_inner_gps)

        ! set the size of the ghost point arrays to be the amount of points total, plus a factor of 2 buffer
        $:GPU_UPDATE(device='[num_gps, num_inner_gps]')
        @:ALLOCATE(ghost_points(1:int((max_num_gps + max_num_inner_gps) * 2.0)))
        @:ALLOCATE(inner_points(1:int((max_num_gps + max_num_inner_gps) * 2.0)))

        $:GPU_ENTER_DATA(copyin='[ghost_points,inner_points]')

        call s_find_ghost_points(ghost_points, inner_points)
        $:GPU_UPDATE(device='[ghost_points, inner_points]')

        call s_compute_image_points(ghost_points, levelset, levelset_norm)
        $:GPU_UPDATE(device='[ghost_points]')

        call s_compute_interpolation_coeffs(ghost_points)
        $:GPU_UPDATE(device='[ghost_points]')

    end subroutine s_ibm_setup

    subroutine s_populate_ib_buffers()

        #:for DIRC, DIRI in [('x', 1), ('y', 2), ('z', 3)]
            #:for LOCC, LOCI in [('beg', -1), ('end', 1)]
                if (bc_${DIRC}$%${LOCC}$ >= 0) then
                    call s_mpi_sendrecv_ib_buffers(ib_markers, ${DIRI}$, ${LOCI}$)
                end if
            #:endfor
        #:endfor

    end subroutine s_populate_ib_buffers

    !>  Subroutine that updates the conservative variables at the ghost points
        !!  @param q_cons_vf Conservative Variables
        !!  @param q_prim_vf Primitive variables
        !!  @param pb Internal bubble pressure
        !!  @param mv Mass of vapor in bubble
    subroutine s_ibm_correct_state(q_cons_vf, q_prim_vf, pb_in, mv_in)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf !< Primitive Variables

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_prim_vf !< Primitive Variables

        real(stp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), optional, intent(INOUT) :: pb_in, mv_in

        integer :: i, j, k, l, q, r!< Iterator variables
        integer :: patch_id !< Patch ID of ghost point
        real(wp) :: rho, gamma, pi_inf, dyn_pres !< Mixture variables
        real(wp), dimension(2) :: Re_K
        real(wp) :: G_K
        real(wp) :: qv_K
        real(wp), dimension(num_fluids) :: Gs

        real(wp) :: pres_IP
        real(wp), dimension(3) :: vel_IP, vel_norm_IP
        real(wp) :: c_IP
        real(wp), dimension(num_fluids) :: alpha_rho_IP, alpha_IP
        real(wp), dimension(nb) :: r_IP, v_IP, pb_IP, mv_IP
        real(wp), dimension(nb*nmom) :: nmom_IP
        real(wp), dimension(nb*nnode) :: presb_IP, massv_IP
        !! Primitive variables at the image point associated with a ghost point,
        !! interpolated from surrounding fluid cells.

        real(wp), dimension(3) :: norm !< Normal vector from GP to IP
        real(wp), dimension(3) :: physical_loc !< Physical loc of GP
        real(wp), dimension(3) :: vel_g !< Velocity of GP
        real(wp), dimension(3) :: radial_vector !< vector from centroid to ghost point
        real(wp), dimension(3) :: rotation_velocity !< speed of the ghost point due to rotation

        real(wp) :: nbub
        real(wp) :: buf
        type(ghost_point) :: gp
        type(ghost_point) :: innerp

        ! set the Moving IBM interior Pressure Values
        $:GPU_PARALLEL_LOOP(private='[i,j,k,patch_id,rho]', copyin='[E_idx,momxb]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    patch_id = ib_markers%sf(j, k, l)
                    if (patch_id /= 0) then
                        q_prim_vf(E_idx)%sf(j, k, l) = 1._wp
                        if (patch_ib(patch_id)%moving_ibm > 0) then
                            rho = 0._wp
                            do i = 1, num_fluids
                                rho = rho + q_prim_vf(contxb + i - 1)%sf(j, k, l)
                            end do

                            ! Sets the momentum
                            do i = 1, num_dims
                                q_cons_vf(momxb + i - 1)%sf(j, k, l) = patch_ib(patch_id)%vel(i)*rho
                                q_prim_vf(momxb + i - 1)%sf(j, k, l) = patch_ib(patch_id)%vel(i)
                            end do
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (num_gps > 0) then
            $:GPU_PARALLEL_LOOP(private='[i,physical_loc,dyn_pres,alpha_rho_IP, alpha_IP,pres_IP,vel_IP,vel_g,vel_norm_IP,r_IP, v_IP,pb_IP,mv_IP,nmom_IP,presb_IP,massv_IP,rho, gamma,pi_inf,Re_K,G_K,Gs,gp,innerp,norm,buf, radial_vector, rotation_velocity, j,k,l,q,qv_K,c_IP,nbub,patch_id]')
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

                !Interpolate primitive variables at image point associated w/ GP
                if (bubbles_euler .and. .not. qbmm) then
                    call s_interpolate_image_point(q_prim_vf, gp, &
                                                   alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, &
                                                   r_IP, v_IP, pb_IP, mv_IP)
                else if (qbmm .and. polytropic) then
                    call s_interpolate_image_point(q_prim_vf, gp, &
                                                   alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, &
                                                   r_IP, v_IP, pb_IP, mv_IP, nmom_IP)
                else if (qbmm .and. .not. polytropic) then
                    call s_interpolate_image_point(q_prim_vf, gp, &
                                                   alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, &
                                                   r_IP, v_IP, pb_IP, mv_IP, nmom_IP, pb_in, mv_in, presb_IP, massv_IP)
                else
                    call s_interpolate_image_point(q_prim_vf, gp, &
                                                   alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP)
                end if

                dyn_pres = 0._wp

                ! Set q_prim_vf params at GP so that mixture vars calculated properly
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    q_prim_vf(q)%sf(j, k, l) = alpha_rho_IP(q)
                    q_prim_vf(advxb + q - 1)%sf(j, k, l) = alpha_IP(q)
                end do

                if (surface_tension) then
                    q_prim_vf(c_idx)%sf(j, k, l) = c_IP
                end if

                ! set the pressure
                if (patch_ib(patch_id)%moving_ibm <= 1) then
                    q_prim_vf(E_idx)%sf(j, k, l) = pres_IP
                else
                    q_prim_vf(E_idx)%sf(j, k, l) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, num_fluids
                        ! Se the pressure inside a moving immersed boundary based upon the pressure of the image point. acceleration, and normal vector direction
                        q_prim_vf(E_idx)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l) + pres_IP/(1._wp - 2._wp*abs(levelset%sf(j, k, l, patch_id)*alpha_rho_IP(q)/pres_IP)*dot_product(patch_ib(patch_id)%force/patch_ib(patch_id)%mass, levelset_norm%sf(j, k, l, patch_id, :)))
                    end do
                end if

                if (model_eqns /= 4) then
                    ! If in simulation, use acc mixture subroutines
                    if (elasticity) then
                        call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv_K, alpha_IP, &
                                                                        alpha_rho_IP, Re_K, G_K, Gs)
                    else
                        call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv_K, alpha_IP, &
                                                                        alpha_rho_IP, Re_K)
                    end if
                end if

                ! Calculate velocity of ghost cell
                if (gp%slip) then
                    norm(1:3) = levelset_norm%sf(gp%loc(1), gp%loc(2), gp%loc(3), gp%ib_patch_id, 1:3)
                    buf = sqrt(sum(norm**2))
                    norm = norm/buf
                    vel_norm_IP = sum(vel_IP*norm)*norm
                    vel_g = vel_IP - vel_norm_IP
                    if (patch_ib(patch_id)%moving_ibm /= 0) then
                        ! compute the linear velocity of the ghost point due to rotation
                        radial_vector = physical_loc - [patch_ib(patch_id)%x_centroid, &
                                                        patch_ib(patch_id)%y_centroid, patch_ib(patch_id)%z_centroid]
                        call s_cross_product(matmul(patch_ib(patch_id)%rotation_matrix, patch_ib(patch_id)%angular_vel), radial_vector, rotation_velocity)

                        ! add only the component of the IB's motion that is normal to the surface
                        vel_g = vel_g + sum((patch_ib(patch_id)%vel + rotation_velocity)*norm)*norm
                    end if
                else
                    if (patch_ib(patch_id)%moving_ibm == 0) then
                        ! we know the object is not moving if moving_ibm is 0 (false)
                        vel_g = 0._wp
                    else
                        ! get the vector that points from the centroid to the ghost
                        radial_vector = physical_loc - [patch_ib(patch_id)%x_centroid, &
                                                        patch_ib(patch_id)%y_centroid, patch_ib(patch_id)%z_centroid]
                        ! convert the angular velocity from the inertial reference frame to the fluids frame, then convert to linear velocity
                        call s_cross_product(matmul(patch_ib(patch_id)%rotation_matrix, patch_ib(patch_id)%angular_vel), radial_vector, rotation_velocity)
                        do q = 1, 3
                            ! if mibm is 1 or 2, then the boundary may be moving
                            vel_g(q) = patch_ib(patch_id)%vel(q) ! add the linear velocity
                            vel_g(q) = vel_g(q) + rotation_velocity(q) ! add the rotational velocity
                        end do
                    end if
                end if

                ! Set momentum
                $:GPU_LOOP(parallelism='[seq]')
                do q = momxb, momxe
                    q_cons_vf(q)%sf(j, k, l) = rho*vel_g(q - momxb + 1)
                    dyn_pres = dyn_pres + q_cons_vf(q)%sf(j, k, l)* &
                               vel_g(q - momxb + 1)/2._wp
                end do

                ! Set continuity and adv vars
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    q_cons_vf(q)%sf(j, k, l) = alpha_rho_IP(q)
                    q_cons_vf(advxb + q - 1)%sf(j, k, l) = alpha_IP(q)
                end do

                ! Set color function
                if (surface_tension) then
                    q_cons_vf(c_idx)%sf(j, k, l) = c_IP
                end if

                ! Set Energy
                if (bubbles_euler) then
                    q_cons_vf(E_idx)%sf(j, k, l) = (1 - alpha_IP(1))*(gamma*pres_IP + pi_inf + dyn_pres)
                else
                    q_cons_vf(E_idx)%sf(j, k, l) = gamma*pres_IP + pi_inf + dyn_pres
                end if
                ! Set bubble vars
                if (bubbles_euler .and. .not. qbmm) then
                    call s_comp_n_from_prim(alpha_IP(1), r_IP, nbub, weight)
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        q_cons_vf(bubxb + (q - 1)*2)%sf(j, k, l) = nbub*r_IP(q)
                        q_cons_vf(bubxb + (q - 1)*2 + 1)%sf(j, k, l) = nbub*v_IP(q)
                        if (.not. polytropic) then
                            q_cons_vf(bubxb + (q - 1)*4)%sf(j, k, l) = nbub*r_IP(q)
                            q_cons_vf(bubxb + (q - 1)*4 + 1)%sf(j, k, l) = nbub*v_IP(q)
                            q_cons_vf(bubxb + (q - 1)*4 + 2)%sf(j, k, l) = nbub*pb_IP(q)
                            q_cons_vf(bubxb + (q - 1)*4 + 3)%sf(j, k, l) = nbub*mv_IP(q)
                        end if
                    end do
                end if

                if (qbmm) then

                    nbub = nmom_IP(1)
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb*nmom
                        q_cons_vf(bubxb + q - 1)%sf(j, k, l) = nbub*nmom_IP(q)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        q_cons_vf(bubxb + (q - 1)*nmom)%sf(j, k, l) = nbub
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

                if (model_eqns == 3) then
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = intxb, intxe
                        q_cons_vf(q)%sf(j, k, l) = alpha_IP(q - intxb + 1)*(gammas(q - intxb + 1)*pres_IP &
                                                                            + pi_infs(q - intxb + 1))
                    end do
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        !Correct the state of the inner points in IBs
        if (num_inner_gps > 0) then
            $:GPU_PARALLEL_LOOP(private='[i,physical_loc,dyn_pres,alpha_rho_IP, alpha_IP,vel_g,rho,gamma,pi_inf,Re_K,innerp,j,k,l,q]')
            do i = 1, num_inner_gps

                innerp = inner_points(i)
                j = innerp%loc(1)
                k = innerp%loc(2)
                l = innerp%loc(3)

                $:GPU_LOOP(parallelism='[seq]')
                do q = momxb, momxe
                    q_cons_vf(q)%sf(j, k, l) = 0._wp
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_ibm_correct_state

    !>  Function that computes the image points for each ghost point
        !!  @param ghost_points Ghost Points
        !!  @param levelset Closest distance from each grid cell to IB
        !!  @param levelset_norm Vector pointing in the direction of the closest distance
    impure subroutine s_compute_image_points(ghost_points_in, levelset_in, levelset_norm_in)

        type(ghost_point), dimension(num_gps), intent(INOUT) :: ghost_points_in
        type(levelset_field), intent(IN) :: levelset_in
        type(levelset_norm_field), intent(IN) :: levelset_norm_in

        real(wp) :: dist
        real(wp), dimension(3) :: norm
        real(wp), dimension(3) :: physical_loc
        real(wp) :: temp_loc
        real(wp), pointer, dimension(:) :: s_cc => null()
        integer :: bound
        type(ghost_point) :: gp

        integer :: q, dim !< Iterator variables
        integer :: i, j, k, l !< Location indexes
        integer :: patch_id !< IB Patch ID
        integer :: dir
        integer :: index

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
            dist = abs(real(levelset_in%sf(i, j, k, patch_id), kind=wp))
            norm(:) = levelset_norm_in%sf(i, j, k, patch_id, :)
            ghost_points_in(q)%ip_loc(:) = physical_loc(:) + 2*dist*norm(:)

            ! Find the closest grid point to the image point
            do dim = 1, num_dims

                ! s_cc points to the dim array we need
                if (dim == 1) then
                    s_cc => x_cc
                    bound = m + buff_size - 1
                elseif (dim == 2) then
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
                    do while ((temp_loc < s_cc(index) &
                               .or. temp_loc > s_cc(index + 1)))
                        index = index + dir
                        if (index < -buff_size .or. index > bound) then
                            print *, "temp_loc=", temp_loc, " s_cc(index)=", s_cc(index), " s_cc(index+1)=", s_cc(index + 1)
                            print *, "Increase buff_size further in m_helper_basic (currently set to a minimum of 10)"
                            error stop "Increase buff_size"
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

    end subroutine s_compute_image_points

    !> Function that finds the number of ghost points, used for allocating
    !! memory.
    subroutine s_find_num_ghost_points(num_gps_out, num_inner_gps_out)

        integer, intent(out) :: num_gps_out
        integer, intent(out) :: num_inner_gps_out

        integer, dimension(2*gp_layers + 1, 2*gp_layers + 1) &
            :: subsection_2D
        integer, dimension(2*gp_layers + 1, 2*gp_layers + 1, 2*gp_layers + 1) &
            :: subsection_3D
        integer :: i, j, k!< Iterator variables

        num_gps_out = 0
        num_inner_gps_out = 0

        do i = 0, m
            do j = 0, n
                if (p == 0) then
                    if (ib_markers%sf(i, j, 0) /= 0) then
                        subsection_2D = ib_markers%sf( &
                                        i - gp_layers:i + gp_layers, &
                                        j - gp_layers:j + gp_layers, 0)
                        if (any(subsection_2D == 0)) then
                            num_gps_out = num_gps_out + 1
                        else
                            num_inner_gps_out = num_inner_gps_out + 1
                        end if
                    end if
                else
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) /= 0) then
                            subsection_3D = ib_markers%sf( &
                                            i - gp_layers:i + gp_layers, &
                                            j - gp_layers:j + gp_layers, &
                                            k - gp_layers:k + gp_layers)
                            if (any(subsection_3D == 0)) then
                                num_gps_out = num_gps_out + 1
                            else
                                num_inner_gps_out = num_inner_gps_out + 1
                            end if
                        end if
                    end do
                end if
            end do
        end do

    end subroutine s_find_num_ghost_points

    !> Function that finds the ghost points
    subroutine s_find_ghost_points(ghost_points_in, inner_points_in)

        type(ghost_point), dimension(num_gps), intent(INOUT) :: ghost_points_in
        type(ghost_point), dimension(num_inner_gps), intent(INOUT) :: inner_points_in
        integer, dimension(2*gp_layers + 1, 2*gp_layers + 1) &
            :: subsection_2D
        integer, dimension(2*gp_layers + 1, 2*gp_layers + 1, 2*gp_layers + 1) &
            :: subsection_3D
        integer :: i, j, k !< Iterator variables
        integer :: count, count_i
        integer :: patch_id

        count = 1
        count_i = 1

        do i = 0, m
            do j = 0, n
                if (p == 0) then
                    ! 2D
                    if (ib_markers%sf(i, j, 0) /= 0) then
                        subsection_2D = ib_markers%sf( &
                                        i - gp_layers:i + gp_layers, &
                                        j - gp_layers:j + gp_layers, 0)
                        if (any(subsection_2D == 0)) then
                            ghost_points_in(count)%loc = [i, j, 0]
                            patch_id = ib_markers%sf(i, j, 0)
                            ghost_points_in(count)%ib_patch_id = &
                                patch_id

                            ghost_points_in(count)%slip = patch_ib(patch_id)%slip
                            ! ghost_points(count)%rank = proc_rank

                            if ((x_cc(i) - dx(i)) < x_domain%beg) then
                                ghost_points_in(count)%DB(1) = -1
                            else if ((x_cc(i) + dx(i)) > x_domain%end) then
                                ghost_points_in(count)%DB(1) = 1
                            else
                                ghost_points_in(count)%DB(1) = 0
                            end if

                            if ((y_cc(j) - dy(j)) < y_domain%beg) then
                                ghost_points_in(count)%DB(2) = -1
                            else if ((y_cc(j) + dy(j)) > y_domain%end) then
                                ghost_points_in(count)%DB(2) = 1
                            else
                                ghost_points_in(count)%DB(2) = 0
                            end if

                            count = count + 1

                        else
                            inner_points_in(count_i)%loc = [i, j, 0]
                            patch_id = ib_markers%sf(i, j, 0)
                            inner_points_in(count_i)%ib_patch_id = &
                                patch_id
                            inner_points_in(count_i)%slip = patch_ib(patch_id)%slip
                            count_i = count_i + 1

                        end if
                    end if
                else
                    ! 3D
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) /= 0) then
                            subsection_3D = ib_markers%sf( &
                                            i - gp_layers:i + gp_layers, &
                                            j - gp_layers:j + gp_layers, &
                                            k - gp_layers:k + gp_layers)
                            if (any(subsection_3D == 0)) then
                                ghost_points_in(count)%loc = [i, j, k]
                                patch_id = ib_markers%sf(i, j, k)
                                ghost_points_in(count)%ib_patch_id = &
                                    ib_markers%sf(i, j, k)
                                ghost_points_in(count)%slip = patch_ib(patch_id)%slip

                                if ((x_cc(i) - dx(i)) < x_domain%beg) then
                                    ghost_points_in(count)%DB(1) = -1
                                else if ((x_cc(i) + dx(i)) > x_domain%end) then
                                    ghost_points_in(count)%DB(1) = 1
                                else
                                    ghost_points_in(count)%DB(1) = 0
                                end if

                                if ((y_cc(j) - dy(j)) < y_domain%beg) then
                                    ghost_points_in(count)%DB(2) = -1
                                else if ((y_cc(j) + dy(j)) > y_domain%end) then
                                    ghost_points_in(count)%DB(2) = 1
                                else
                                    ghost_points_in(count)%DB(2) = 0
                                end if

                                if ((z_cc(k) - dz(k)) < z_domain%beg) then
                                    ghost_points_in(count)%DB(3) = -1
                                else if ((z_cc(k) + dz(k)) > z_domain%end) then
                                    ghost_points_in(count)%DB(3) = 1
                                else
                                    ghost_points_in(count)%DB(3) = 0
                                end if

                                count = count + 1
                            else
                                inner_points_in(count_i)%loc = [i, j, k]
                                patch_id = ib_markers%sf(i, j, k)
                                inner_points_in(count_i)%ib_patch_id = &
                                    ib_markers%sf(i, j, k)
                                inner_points_in(count_i)%slip = patch_ib(patch_id)%slip

                                count_i = count_i + 1
                            end if
                        end if
                    end do
                end if
            end do
        end do

    end subroutine s_find_ghost_points

    !>  Function that computes the interpolation coefficients of image points
    subroutine s_compute_interpolation_coeffs(ghost_points_in)

        type(ghost_point), dimension(num_gps), intent(INOUT) :: ghost_points_in

        real(wp), dimension(2, 2, 2) :: dist
        real(wp), dimension(2, 2, 2) :: alpha
        real(wp), dimension(2, 2, 2) :: interp_coeffs
        real(wp) :: buf
        real(wp), dimension(2, 2, 2) :: eta
        type(ghost_point) :: gp
        integer :: i !< Iterator variables
        integer :: i1, i2, j1, j2, k1, k2 !< Grid indexes
        integer :: patch_id

        ! 2D
        if (p <= 0) then
            do i = 1, num_gps
                gp = ghost_points_in(i)
                ! Get the interpolation points
                i1 = gp%ip_grid(1); i2 = i1 + 1
                j1 = gp%ip_grid(2); j2 = j1 + 1

                dist = 0._wp
                buf = 1._wp
                dist(1, 1, 1) = sqrt( &
                                (x_cc(i1) - gp%ip_loc(1))**2 + &
                                (y_cc(j1) - gp%ip_loc(2))**2)
                dist(2, 1, 1) = sqrt( &
                                (x_cc(i2) - gp%ip_loc(1))**2 + &
                                (y_cc(j1) - gp%ip_loc(2))**2)
                dist(1, 2, 1) = sqrt( &
                                (x_cc(i1) - gp%ip_loc(1))**2 + &
                                (y_cc(j2) - gp%ip_loc(2))**2)
                dist(2, 2, 1) = sqrt( &
                                (x_cc(i2) - gp%ip_loc(1))**2 + &
                                (y_cc(j2) - gp%ip_loc(2))**2)

                interp_coeffs = 0._wp

                if (dist(1, 1, 1) <= 1.e-16_wp) then
                    interp_coeffs(1, 1, 1) = 1._wp
                else if (dist(2, 1, 1) <= 1.e-16_wp) then
                    interp_coeffs(2, 1, 1) = 1._wp
                else if (dist(1, 2, 1) <= 1.e-16_wp) then
                    interp_coeffs(1, 2, 1) = 1._wp
                else if (dist(2, 2, 1) <= 1.e-16_wp) then
                    interp_coeffs(2, 2, 1) = 1._wp
                else
                    eta(:, :, 1) = 1._wp/dist(:, :, 1)**2
                    alpha = 1._wp
                    patch_id = gp%ib_patch_id
                    if (ib_markers%sf(i1, j1, 0) /= 0) alpha(1, 1, 1) = 0._wp
                    if (ib_markers%sf(i2, j1, 0) /= 0) alpha(2, 1, 1) = 0._wp
                    if (ib_markers%sf(i1, j2, 0) /= 0) alpha(1, 2, 1) = 0._wp
                    if (ib_markers%sf(i2, j2, 0) /= 0) alpha(2, 2, 1) = 0._wp
                    buf = sum(alpha(:, :, 1)*eta(:, :, 1))
                    if (buf > 0._wp) then
                        interp_coeffs(:, :, 1) = alpha(:, :, 1)*eta(:, :, 1)/buf
                    else
                        buf = sum(eta(:, :, 1))
                        interp_coeffs(:, :, 1) = eta(:, :, 1)/buf
                    end if
                end if

                ghost_points_in(i)%interp_coeffs = interp_coeffs
            end do

        else
            do i = 1, num_gps
                gp = ghost_points_in(i)
                ! Get the interpolation points
                i1 = gp%ip_grid(1); i2 = i1 + 1
                j1 = gp%ip_grid(2); j2 = j1 + 1
                k1 = gp%ip_grid(3); k2 = k1 + 1

                ! Get interpolation weights (Chaudhuri et al. 2011, JCP)
                dist(1, 1, 1) = sqrt( &
                                (x_cc(i1) - gp%ip_loc(1))**2 + &
                                (y_cc(j1) - gp%ip_loc(2))**2 + &
                                (z_cc(k1) - gp%ip_loc(3))**2)
                dist(2, 1, 1) = sqrt( &
                                (x_cc(i2) - gp%ip_loc(1))**2 + &
                                (y_cc(j1) - gp%ip_loc(2))**2 + &
                                (z_cc(k1) - gp%ip_loc(3))**2)
                dist(1, 2, 1) = sqrt( &
                                (x_cc(i1) - gp%ip_loc(1))**2 + &
                                (y_cc(j2) - gp%ip_loc(2))**2 + &
                                (z_cc(k1) - gp%ip_loc(3))**2)
                dist(2, 2, 1) = sqrt( &
                                (x_cc(i2) - gp%ip_loc(1))**2 + &
                                (y_cc(j2) - gp%ip_loc(2))**2 + &
                                (z_cc(k1) - gp%ip_loc(3))**2)
                dist(1, 1, 2) = sqrt( &
                                (x_cc(i1) - gp%ip_loc(1))**2 + &
                                (y_cc(j1) - gp%ip_loc(2))**2 + &
                                (z_cc(k2) - gp%ip_loc(3))**2)
                dist(2, 1, 2) = sqrt( &
                                (x_cc(i2) - gp%ip_loc(1))**2 + &
                                (y_cc(j1) - gp%ip_loc(2))**2 + &
                                (z_cc(k2) - gp%ip_loc(3))**2)
                dist(1, 2, 2) = sqrt( &
                                (x_cc(i1) - gp%ip_loc(1))**2 + &
                                (y_cc(j2) - gp%ip_loc(2))**2 + &
                                (z_cc(k2) - gp%ip_loc(3))**2)
                dist(2, 2, 2) = sqrt( &
                                (x_cc(i2) - gp%ip_loc(1))**2 + &
                                (y_cc(j2) - gp%ip_loc(2))**2 + &
                                (z_cc(k2) - gp%ip_loc(3))**2)
                interp_coeffs = 0._wp
                buf = 1._wp
                if (dist(1, 1, 1) <= 1.e-16_wp) then
                    interp_coeffs(1, 1, 1) = 1._wp
                else if (dist(2, 1, 1) <= 1.e-16_wp) then
                    interp_coeffs(2, 1, 1) = 1._wp
                else if (dist(1, 2, 1) <= 1.e-16_wp) then
                    interp_coeffs(1, 2, 1) = 1._wp
                else if (dist(2, 2, 1) <= 1.e-16_wp) then
                    interp_coeffs(2, 2, 1) = 1._wp
                else if (dist(1, 1, 2) <= 1.e-16_wp) then
                    interp_coeffs(1, 1, 2) = 1._wp
                else if (dist(2, 1, 2) <= 1.e-16_wp) then
                    interp_coeffs(2, 1, 2) = 1._wp
                else if (dist(1, 2, 2) <= 1.e-16_wp) then
                    interp_coeffs(1, 2, 2) = 1._wp
                else if (dist(2, 2, 2) <= 1.e-16_wp) then
                    interp_coeffs(2, 2, 2) = 1._wp
                else
                    eta = 1._wp/dist**2
                    alpha = 1._wp
                    if (ib_markers%sf(i1, j1, k1) /= 0) alpha(1, 1, 1) = 0._wp
                    if (ib_markers%sf(i2, j1, k1) /= 0) alpha(2, 1, 1) = 0._wp
                    if (ib_markers%sf(i1, j2, k1) /= 0) alpha(1, 2, 1) = 0._wp
                    if (ib_markers%sf(i2, j2, k1) /= 0) alpha(2, 2, 1) = 0._wp
                    if (ib_markers%sf(i1, j1, k2) /= 0) alpha(1, 1, 2) = 0._wp
                    if (ib_markers%sf(i2, j1, k2) /= 0) alpha(2, 1, 2) = 0._wp
                    if (ib_markers%sf(i1, j2, k2) /= 0) alpha(1, 2, 2) = 0._wp
                    if (ib_markers%sf(i2, j2, k2) /= 0) alpha(2, 2, 2) = 0._wp
                    buf = sum(alpha*eta)
                    if (buf > 0._wp) then
                        interp_coeffs = alpha*eta/buf
                    else
                        buf = sum(eta)
                        interp_coeffs = eta/buf
                    end if
                end if

                ghost_points_in(i)%interp_coeffs = interp_coeffs
            end do
        end if

    end subroutine s_compute_interpolation_coeffs

    !> Function that uses the interpolation coefficients and the current state
    !! at the cell centers in order to estimate the state at the image point
    subroutine s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, &
                                         pres_IP, vel_IP, c_IP, r_IP, v_IP, pb_IP, &
                                         mv_IP, nmom_IP, pb_in, mv_in, presb_IP, massv_IP)
        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf !< Primitive Variables

        real(stp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(IN) :: pb_in, mv_in

        type(ghost_point), intent(IN) :: gp
        real(wp), intent(INOUT) :: pres_IP
        real(wp), dimension(3), intent(INOUT) :: vel_IP
        real(wp), intent(INOUT) :: c_IP
        real(wp), dimension(num_fluids), intent(INOUT) :: alpha_IP, alpha_rho_IP
        real(wp), optional, dimension(:), intent(INOUT) :: r_IP, v_IP, pb_IP, mv_IP
        real(wp), optional, dimension(:), intent(INOUT) :: nmom_IP
        real(wp), optional, dimension(:), intent(INOUT) :: presb_IP, massv_IP

        integer :: i, j, k, l, q !< Iterator variables
        integer :: i1, i2, j1, j2, k1, k2 !< Iterator variables
        real(wp) :: coeff

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

                    pres_IP = pres_IP + coeff* &
                              q_prim_vf(E_idx)%sf(i, j, k)

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = momxb, momxe
                        vel_IP(q + 1 - momxb) = vel_IP(q + 1 - momxb) + coeff* &
                                                q_prim_vf(q)%sf(i, j, k)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do l = contxb, contxe
                        alpha_rho_IP(l) = alpha_rho_IP(l) + coeff* &
                                          q_prim_vf(l)%sf(i, j, k)
                        alpha_IP(l) = alpha_IP(l) + coeff* &
                                      q_prim_vf(advxb + l - 1)%sf(i, j, k)
                    end do

                    if (surface_tension) then
                        c_IP = c_IP + coeff*q_prim_vf(c_idx)%sf(i, j, k)
                    end if

                    if (bubbles_euler .and. .not. qbmm) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do l = 1, nb
                            if (polytropic) then
                                r_IP(l) = r_IP(l) + coeff*q_prim_vf(bubxb + (l - 1)*2)%sf(i, j, k)
                                v_IP(l) = v_IP(l) + coeff*q_prim_vf(bubxb + 1 + (l - 1)*2)%sf(i, j, k)
                            else
                                r_IP(l) = r_IP(l) + coeff*q_prim_vf(bubxb + (l - 1)*4)%sf(i, j, k)
                                v_IP(l) = v_IP(l) + coeff*q_prim_vf(bubxb + 1 + (l - 1)*4)%sf(i, j, k)
                                pb_IP(l) = pb_IP(l) + coeff*q_prim_vf(bubxb + 2 + (l - 1)*4)%sf(i, j, k)
                                mv_IP(l) = mv_IP(l) + coeff*q_prim_vf(bubxb + 3 + (l - 1)*4)%sf(i, j, k)
                            end if
                        end do
                    end if

                    if (qbmm) then
                        do l = 1, nb*nmom
                            nmom_IP(l) = nmom_IP(l) + coeff*q_prim_vf(bubxb - 1 + l)%sf(i, j, k)
                        end do
                        if (.not. polytropic) then
                            do q = 1, nb
                                do l = 1, nnode
                                    presb_IP((q - 1)*nnode + l) = presb_IP((q - 1)*nnode + l) + &
                                                                  coeff*real(pb_in(i, j, k, l, q), kind=wp)
                                    massv_IP((q - 1)*nnode + l) = massv_IP((q - 1)*nnode + l) + &
                                                                  coeff*real(mv_in(i, j, k, l, q), kind=wp)
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
    impure subroutine s_update_mib(num_ibs, levelset, levelset_norm)

        integer, intent(in) :: num_ibs
        type(levelset_field), intent(inout) :: levelset
        type(levelset_norm_field), intent(inout) :: levelset_norm

        integer :: i, ierr

        ! Clears the existing immersed boundary indices
        ib_markers%sf = 0._wp
        levelset%sf = 0._wp
        levelset_norm%sf = 0._wp

        ! recalulcate the rotation matrix based upon the new angles
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) then
                call s_update_ib_rotation_matrix(i)
            end if
        end do

        $:GPU_UPDATE(device='[patch_ib]')

        ! recompute the new ib_patch locations and broadcast them.
        call s_apply_ib_patches(ib_markers%sf(0:m, 0:n, 0:p), levelset, levelset_norm)
        call s_populate_ib_buffers() ! transmits the new IB markers via MPI
        $:GPU_UPDATE(device='[ib_markers%sf]')
        $:GPU_UPDATE(host='[levelset%sf, levelset_norm%sf]')

        ! recalculate the ghost point locations and coefficients
        call s_find_num_ghost_points(num_gps, num_inner_gps)
        $:GPU_UPDATE(device='[num_gps, num_inner_gps]')

        call s_find_ghost_points(ghost_points, inner_points)
        $:GPU_UPDATE(device='[ghost_points, inner_points]')

        call s_compute_image_points(ghost_points, levelset, levelset_norm)
        $:GPU_UPDATE(device='[ghost_points]')

        call s_compute_interpolation_coeffs(ghost_points)
        $:GPU_UPDATE(device='[ghost_points]')

    end subroutine s_update_mib

    ! compute the surface integrals of the IB via a volume integraion method described in
    ! "A coupled IBM/Euler-Lagrange framework for simulating shock-induced particle size segregation"
    ! by Archana Sridhar and Jesse Capecelatro
    subroutine s_compute_ib_forces(q_prim_vf, fluid_pp)

        ! real(wp), dimension(idwbuff(1)%beg:idwbuff(1)%end, &
        !             idwbuff(2)%beg:idwbuff(2)%end, &
        !             idwbuff(3)%beg:idwbuff(3)%end), intent(in) :: pressure
        type(scalar_field), dimension(1:sys_size), intent(in) :: q_prim_vf
        type(physical_parameters), dimension(1:num_fluids), intent(in) :: fluid_pp

        integer :: gp_id, i, j, k, l, q, ib_idx, fluid_idx
        real(wp), dimension(num_ibs, 3) :: forces, torques
        real(wp), dimension(1:3, 1:3) :: viscous_stress_div, viscous_stress_div_1, viscous_stress_div_2, viscous_cross_1, viscous_cross_2 ! viscous stress tensor with temp vectors to hold divergence calculations
        real(wp), dimension(1:3) :: local_force_contribution, radial_vector, local_torque_contribution, vel
        real(wp) :: cell_volume, dx, dy, dz, dynamic_viscosity
        real(wp), dimension(1:num_fluids) :: dynamic_viscosities

        forces = 0._wp
        torques = 0._wp

        if (viscous) then
            do fluid_idx = 1, num_fluids
                if (fluid_pp(fluid_idx)%Re(1) /= 0._wp) then
                    dynamic_viscosities(fluid_idx) = 1._wp/fluid_pp(fluid_idx)%Re(1)
                else
                    dynamic_viscosities(fluid_idx) = 0._wp
                end if
            end do
        end if

        $:GPU_PARALLEL_LOOP(private='[ib_idx,fluid_idx, radial_vector,local_force_contribution,cell_volume,local_torque_contribution, dynamic_viscosity, viscous_stress_div, viscous_stress_div_1, viscous_stress_div_2, viscous_cross_1, viscous_cross_2, dx, dy, dz]', copy='[forces,torques]', copyin='[ib_markers,patch_ib,dynamic_viscosities]', collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    ib_idx = ib_markers%sf(i, j, k)
                    if (ib_idx /= 0) then
                        ! get the vector pointing to the grid cell from the IB centroid
                        if (num_dims == 3) then
                            radial_vector = [x_cc(i), y_cc(j), z_cc(k)] - [patch_ib(ib_idx)%x_centroid, patch_ib(ib_idx)%y_centroid, patch_ib(ib_idx)%z_centroid]
                        else
                            radial_vector = [x_cc(i), y_cc(j), 0._wp] - [patch_ib(ib_idx)%x_centroid, patch_ib(ib_idx)%y_centroid, 0._wp]
                        end if
                        dx = x_cc(i + 1) - x_cc(i)
                        dy = y_cc(j + 1) - y_cc(j)

                        local_force_contribution(:) = 0._wp
                        do fluid_idx = 0, num_fluids - 1
                            ! Get the pressure contribution to force via a finite difference to compute the 2D components of the gradient of the pressure and cell volume
                            local_force_contribution(1) = local_force_contribution(1) - (q_prim_vf(E_idx + fluid_idx)%sf(i + 1, j, k) - q_prim_vf(E_idx + fluid_idx)%sf(i - 1, j, k))/(2._wp*dx) ! force is the negative pressure gradient
                            local_force_contribution(2) = local_force_contribution(2) - (q_prim_vf(E_idx + fluid_idx)%sf(i, j + 1, k) - q_prim_vf(E_idx + fluid_idx)%sf(i, j - 1, k))/(2._wp*dy)
                            cell_volume = abs(dx*dy)
                            ! add the 3D component of the pressure gradient, if we are working in 3 dimensions
                            if (num_dims == 3) then
                                dz = z_cc(k + 1) - z_cc(k)
                                local_force_contribution(3) = local_force_contribution(3) - (q_prim_vf(E_idx + fluid_idx)%sf(i, j, k + 1) - q_prim_vf(E_idx + fluid_idx)%sf(i, j, k - 1))/(2._wp*dz)
                                cell_volume = abs(cell_volume*dz)
                            end if
                        end do

                        ! Update the force values atomically to prevent race conditions
                        call s_cross_product(radial_vector, local_force_contribution, local_torque_contribution)

                        ! get the viscous stress and add its contribution if that is considered
                        ! TODO :: This is really bad code
                        if (viscous) then
                            ! compute the volume-weighted local dynamic viscosity
                            dynamic_viscosity = 0._wp
                            do fluid_idx = 1, num_fluids
                                ! local dynamic viscosity is the dynamic viscosity of the fluid times alpha of the fluid
                                dynamic_viscosity = dynamic_viscosity + (q_prim_vf(fluid_idx + advxb - 1)%sf(i, j, k)*dynamic_viscosities(fluid_idx))
                            end do

                            ! get the linear force component first
                            call s_compute_viscous_stress_tensor(viscous_stress_div_1, q_prim_vf, dynamic_viscosity, i - 1, j, k)
                            call s_compute_viscous_stress_tensor(viscous_stress_div_2, q_prim_vf, dynamic_viscosity, i + 1, j, k)
                            viscous_stress_div = (viscous_stress_div_2 - viscous_stress_div_1)/(2._wp*dx) ! get the x derivative of the viscous stress tensor
                            local_force_contribution(1:3) = local_force_contribution(1:3) + viscous_stress_div(1, 1:3) ! add te x components of the derivative to the force
                            do l = 1, 3
                                ! take the cross products for the torque component
                                call s_cross_product(radial_vector, viscous_stress_div_1(l, 1:3), viscous_cross_1(l, 1:3))
                                call s_cross_product(radial_vector, viscous_stress_div_2(l, 1:3), viscous_cross_2(l, 1:3))
                            end do

                            viscous_stress_div = (viscous_cross_2 - viscous_cross_1)/(2._wp*dx) ! get the x derivative of the cross product
                            local_torque_contribution(1:3) = local_torque_contribution(1:3) + viscous_stress_div(1, 1:3) ! apply the cross product derivative to the torque

                            call s_compute_viscous_stress_tensor(viscous_stress_div_1, q_prim_vf, dynamic_viscosity, i, j - 1, k)
                            call s_compute_viscous_stress_tensor(viscous_stress_div_2, q_prim_vf, dynamic_viscosity, i, j + 1, k)
                            viscous_stress_div = (viscous_stress_div_2 - viscous_stress_div_1)/(2._wp*dy)
                            local_force_contribution(1:3) = local_force_contribution(1:3) + viscous_stress_div(2, 1:3)
                            do l = 1, 3
                                call s_cross_product(radial_vector, viscous_stress_div_1(l, 1:3), viscous_cross_1(l, 1:3))
                                call s_cross_product(radial_vector, viscous_stress_div_2(l, 1:3), viscous_cross_2(l, 1:3))
                            end do

                            viscous_stress_div = (viscous_cross_2 - viscous_cross_1)/(2._wp*dy)
                            local_torque_contribution(1:3) = local_torque_contribution(1:3) + viscous_stress_div(2, 1:3)

                            if (num_dims == 3) then
                                call s_compute_viscous_stress_tensor(viscous_stress_div_1, q_prim_vf, dynamic_viscosity, i, j, k - 1)
                                call s_compute_viscous_stress_tensor(viscous_stress_div_2, q_prim_vf, dynamic_viscosity, i, j, k + 1)
                                viscous_stress_div = (viscous_stress_div_2 - viscous_stress_div_1)/(2._wp*dz)
                                local_force_contribution(1:3) = local_force_contribution(1:3) + viscous_stress_div(3, 1:3)
                                do l = 1, 3
                                    call s_cross_product(radial_vector, viscous_stress_div_1(l, 1:3), viscous_cross_1(l, 1:3))
                                    call s_cross_product(radial_vector, viscous_stress_div_2(l, 1:3), viscous_cross_2(l, 1:3))
                                end do
                                viscous_stress_div = (viscous_cross_2 - viscous_cross_1)/(2._wp*dz)
                                local_torque_contribution(1:3) = local_torque_contribution(1:3) + viscous_stress_div(3, 1:3)
                            end if
                        end if

                        do l = 1, 3
                            $:GPU_ATOMIC(atomic='update')
                            forces(ib_idx, l) = forces(ib_idx, l) + (local_force_contribution(l)*cell_volume)
                            $:GPU_ATOMIC(atomic='update')
                            torques(ib_idx, l) = torques(ib_idx, l) + local_torque_contribution(l)*cell_volume
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! reduce the forces across all MPI ranks
        call s_mpi_allreduce_vectors_sum(forces, forces, num_ibs, 3)
        call s_mpi_allreduce_vectors_sum(torques, torques, num_ibs, 3)

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
        do i = 1, num_ibs
            patch_ib(i)%force(:) = forces(i, :)
            patch_ib(i)%torque(:) = matmul(patch_ib(i)%rotation_matrix_inverse, torques(i, :)) ! torques must be converted to the local coordinates of the IB
        end do

    end subroutine s_compute_ib_forces

    !> Subroutine to deallocate memory reserved for the IBM module
    impure subroutine s_finalize_ibm_module()

        @:DEALLOCATE(ib_markers%sf)
        @:DEALLOCATE(levelset%sf)
        @:DEALLOCATE(levelset_norm%sf)

    end subroutine s_finalize_ibm_module

    !> Computes the center of mass for IB patch types where we are unable to determine their center of mass analytically.
    !> These patches include things like NACA airfoils and STL models
    subroutine s_compute_centroid_offset(ib_marker)

        integer, intent(in) :: ib_marker

        integer :: i, j, k, num_cells, num_cells_local
        real(wp), dimension(1:3) :: center_of_mass, center_of_mass_local

        ! Offset only needs to be computes for specific geometries
        if (patch_ib(ib_marker)%geometry == 4 .or. &
            patch_ib(ib_marker)%geometry == 5 .or. &
            patch_ib(ib_marker)%geometry == 11 .or. &
            patch_ib(ib_marker)%geometry == 12) then

            center_of_mass_local = [0._wp, 0._wp, 0._wp]
            num_cells_local = 0

            ! get the summed mass distribution and number of cells to divide by
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) == ib_marker) then
                            num_cells_local = num_cells_local + 1
                            center_of_mass_local = center_of_mass_local + [x_cc(i), y_cc(j), 0._wp]
                            if (num_dims == 3) center_of_mass_local(3) = center_of_mass_local(3) + z_cc(k)
                        end if
                    end do
                end do
            end do

            ! reduce the mass contribution over all MPI ranks and compute COM
            call s_mpi_allreduce_integer_sum(num_cells_local, num_cells)
            if (num_cells /= 0) then
                call s_mpi_allreduce_sum(center_of_mass_local(1), center_of_mass(1))
                call s_mpi_allreduce_sum(center_of_mass_local(2), center_of_mass(2))
                call s_mpi_allreduce_sum(center_of_mass_local(3), center_of_mass(3))
                center_of_mass = center_of_mass/real(num_cells, wp)
            else
                patch_ib(ib_marker)%centroid_offset = [0._wp, 0._wp, 0._wp]
                return
            end if

            ! assign the centroid offset as a vector pointing from the true COM to the "centroid" in the input file and replace the current centroid
            patch_ib(ib_marker)%centroid_offset = [patch_ib(ib_marker)%x_centroid, patch_ib(ib_marker)%y_centroid, patch_ib(ib_marker)%z_centroid] &
                                                  - center_of_mass
            patch_ib(ib_marker)%x_centroid = center_of_mass(1)
            patch_ib(ib_marker)%y_centroid = center_of_mass(2)
            patch_ib(ib_marker)%z_centroid = center_of_mass(3)

            ! rotate the centroid offset back into the local coords of the IB
            patch_ib(ib_marker)%centroid_offset = matmul(patch_ib(ib_marker)%rotation_matrix_inverse, patch_ib(ib_marker)%centroid_offset)
        else
            patch_ib(ib_marker)%centroid_offset(:) = [0._wp, 0._wp, 0._wp]
        end if

    end subroutine s_compute_centroid_offset

    subroutine s_compute_moment_of_inertia(ib_marker, axis)

        real(wp), dimension(3), intent(in) :: axis !< the axis about which we compute the moment. Only required in 3D.
        integer, intent(in) :: ib_marker

        real(wp) :: moment, distance_to_axis, cell_volume
        real(wp), dimension(3) :: position, closest_point_along_axis, vector_to_axis, normal_axis
        integer :: i, j, k, count

        if (p == 0) then
            normal_axis = [0, 0, 1]
        else if (sqrt(sum(axis**2)) == 0) then
            ! if the object is not actually rotating at this time, return a dummy value and exit
            patch_ib(ib_marker)%moment = 1._wp
            return
        else
            normal_axis = axis/sqrt(sum(axis))
        end if

        ! if the IB is in 2D or a 3D sphere, we can compute this exactly
        if (patch_ib(ib_marker)%geometry == 2) then ! circle
            patch_ib(ib_marker)%moment = 0.5_wp*patch_ib(ib_marker)%mass*(patch_ib(ib_marker)%radius)**2
        elseif (patch_ib(ib_marker)%geometry == 3) then ! rectangle
            patch_ib(ib_marker)%moment = patch_ib(ib_marker)%mass*(patch_ib(ib_marker)%length_x**2 + patch_ib(ib_marker)%length_y**2)/6._wp
        elseif (patch_ib(ib_marker)%geometry == 6) then ! ellipse
            patch_ib(ib_marker)%moment = 0.0625_wp*patch_ib(ib_marker)%mass*(patch_ib(ib_marker)%length_x**2 + patch_ib(ib_marker)%length_y**2)
        elseif (patch_ib(ib_marker)%geometry == 8) then ! sphere
            patch_ib(ib_marker)%moment = 0.4*patch_ib(ib_marker)%mass*(patch_ib(ib_marker)%radius)**2

        else ! we do not have an analytic moment of inertia calculation and need to approximate it directly via a sum
            count = 0
            moment = 0._wp
            cell_volume = (x_cc(1) - x_cc(0))*(y_cc(1) - y_cc(0)) ! computed without grid stretching. Update in the loop to perform with stretching
            if (p /= 0) then
                cell_volume = cell_volume*(z_cc(1) - z_cc(0))
            end if

            $:GPU_PARALLEL_LOOP(private='[position,closest_point_along_axis,vector_to_axis,distance_to_axis]', copy='[moment,count]', copyin='[ib_marker,cell_volume,normal_axis]', collapse=3)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) == ib_marker) then
                            $:GPU_ATOMIC(atomic='update')
                            count = count + 1 ! increment the count of total cells in the boundary

                            ! get the position in local coordinates so that the axis passes through 0, 0, 0
                            if (p == 0) then
                                position = [x_cc(i), y_cc(j), 0._wp] - [patch_ib(ib_marker)%x_centroid, patch_ib(ib_marker)%y_centroid, 0._wp]
                            else
                                position = [x_cc(i), y_cc(j), z_cc(k)] - [patch_ib(ib_marker)%x_centroid, patch_ib(ib_marker)%y_centroid, patch_ib(ib_marker)%z_centroid]
                            end if

                            ! project the position along the axis to find the closest distance to the rotation axis
                            closest_point_along_axis = normal_axis*dot_product(normal_axis, position)
                            vector_to_axis = position - closest_point_along_axis
                            distance_to_axis = dot_product(vector_to_axis, vector_to_axis) ! saves the distance to the axis squared

                            ! compute the position component of the moment
                            $:GPU_ATOMIC(atomic='update')
                            moment = moment + distance_to_axis
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! write the final moment assuming the points are all uniform density
            patch_ib(ib_marker)%moment = moment*patch_ib(ib_marker)%mass/(count*cell_volume)
            $:GPU_UPDATE(device='[patch_ib(ib_marker)%moment]')
        end if

    end subroutine s_compute_moment_of_inertia

    subroutine s_cross_product(a, b, c)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: a(3), b(3)
        real(wp), intent(out) :: c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end subroutine s_cross_product

end module m_ibm
