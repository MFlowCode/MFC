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

    type(integer_field), public :: ib_markers
    type(levelset_field), public :: levelset
    type(levelset_norm_field), public :: levelset_norm
    $:GPU_DECLARE(create='[ib_markers,levelset,levelset_norm]')

    type(ghost_point), dimension(:), allocatable :: ghost_points
    type(ghost_point), dimension(:), allocatable :: inner_points
    $:GPU_DECLARE(create='[ghost_points,inner_points]')

    integer :: num_gps !< Number of ghost points
    integer :: num_inner_gps !< Number of ghost points
    $:GPU_DECLARE(create='[gp_layers,num_gps,num_inner_gps]')

contains

    !>  Allocates memory for the variables in the IBM module
    impure subroutine s_initialize_ibm_module()

        if (p > 0) then
            @:ALLOCATE(ib_markers%sf(-gp_layers:m+gp_layers, &
                -gp_layers:n+gp_layers, -gp_layers:p+gp_layers))
            @:ALLOCATE(levelset%sf(-gp_layers:m+gp_layers, &
                -gp_layers:n+gp_layers, -gp_layers:p+gp_layers, 1:num_ibs))
            @:ALLOCATE(levelset_norm%sf(-gp_layers:m+gp_layers, &
                -gp_layers:n+gp_layers, -gp_layers:p+gp_layers, 1:num_ibs, 1:3))
        else
            @:ALLOCATE(ib_markers%sf(-gp_layers:m+gp_layers, &
                -gp_layers:n+gp_layers, 0:0))
            @:ALLOCATE(levelset%sf(-gp_layers:m+gp_layers, &
                -gp_layers:n+gp_layers, 0:0, 1:num_ibs))
            @:ALLOCATE(levelset_norm%sf(-gp_layers:m+gp_layers, &
                -gp_layers:n+gp_layers, 0:0, 1:num_ibs, 1:3))
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

        $:GPU_UPDATE(device='[ib_markers%sf]')
        $:GPU_UPDATE(device='[levelset%sf]')
        $:GPU_UPDATE(device='[levelset_norm%sf]')

        ! Get neighboring IB variables from other processors
        call s_populate_ib_buffers()

        $:GPU_UPDATE(host='[ib_markers%sf]')

        call s_find_num_ghost_points(num_gps, num_inner_gps)

        $:GPU_UPDATE(device='[num_gps, num_inner_gps]')
        @:ALLOCATE(ghost_points(1:num_gps))
        @:ALLOCATE(inner_points(1:num_inner_gps))

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
    pure subroutine s_ibm_correct_state(q_cons_vf, q_prim_vf, pb_in, mv_in)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf !< Primitive Variables

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_prim_vf !< Primitive Variables

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), optional, intent(INOUT) :: pb_in, mv_in

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

        real(wp) :: nbub
        real(wp) :: buf
        type(ghost_point) :: gp
        type(ghost_point) :: innerp

        $:GPU_PARALLEL_LOOP(private='[physical_loc,dyn_pres,alpha_rho_IP, &
            & alpha_IP,pres_IP,vel_IP,vel_g,vel_norm_IP,r_IP, &
            & v_IP,pb_IP,mv_IP,nmom_IP,presb_IP,massv_IP,rho, &
            & gamma,pi_inf,Re_K,G_K,Gs,gp,innerp,norm,buf, &
            & j,k,l,q]')
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

            if (model_eqns /= 4) then
                ! If in simulation, use acc mixture subroutines
                if (elasticity) then
                    call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv_K, alpha_IP, &
                                                                    alpha_rho_IP, Re_K, G_K, Gs)
                else if (bubbles_euler) then
                    call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, qv_K, alpha_IP, &
                                                                            alpha_rho_IP, Re_K)
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
            else
                vel_g = 0._wp
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
                do q = 1, nb*nmom
                    q_cons_vf(bubxb + q - 1)%sf(j, k, l) = nbub*nmom_IP(q)
                end do
                do q = 1, nb
                    q_cons_vf(bubxb + (q - 1)*nmom)%sf(j, k, l) = nbub
                end do

                if (.not. polytropic) then
                    do q = 1, nb
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

        !Correct the state of the inner points in IBs
        $:GPU_PARALLEL_LOOP(private='[physical_loc,dyn_pres,alpha_rho_IP, &
            & alpha_IP,vel_g,rho,gamma,pi_inf,Re_K,innerp, &
            & j,k,l,q]')
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
        integer :: i, j, k !< Location indexes
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
            dist = abs(levelset_in%sf(i, j, k, patch_id))
            norm(:) = levelset_norm_in%sf(i, j, k, patch_id, :)
            ghost_points_in(q)%ip_loc(:) = physical_loc(:) + 2*dist*norm(:)

            ! Find the closest grid point to the image point
            do dim = 1, num_dims

                ! s_cc points to the dim array we need
                if (dim == 1) then
                    s_cc => x_cc
                    bound = m
                elseif (dim == 2) then
                    s_cc => y_cc
                    bound = n
                else
                    s_cc => z_cc
                    bound = p
                end if

                if (f_approx_equal(norm(dim), 0._wp)) then
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
                               .or. temp_loc > s_cc(index + 1)) &
                              .and. (index >= 0 .and. index <= bound))
                        index = index + dir
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
    pure subroutine s_find_num_ghost_points(num_gps_out, num_inner_gps_out)

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
    pure subroutine s_find_ghost_points(ghost_points_in, inner_points_in)

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
    pure subroutine s_compute_interpolation_coeffs(ghost_points_in)

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
    pure subroutine s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, pb_IP, mv_IP, nmom_IP, pb_in, mv_in, presb_IP, massv_IP)
        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf !< Primitive Variables

        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(INOUT) :: pb_in, mv_in

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
                                    presb_IP((q - 1)*nnode + l) = presb_IP((q - 1)*nnode + l) + coeff*pb_in(i, j, k, l, q)
                                    massv_IP((q - 1)*nnode + l) = massv_IP((q - 1)*nnode + l) + coeff*mv_in(i, j, k, l, q)
                                end do
                            end do
                        end if

                    end if

                end do
            end do
        end do

    end subroutine s_interpolate_image_point

    !> Subroutine to deallocate memory reserved for the IBM module
    impure subroutine s_finalize_ibm_module()

        @:DEALLOCATE(ib_markers%sf)
        @:DEALLOCATE(levelset%sf)
        @:DEALLOCATE(levelset_norm%sf)

    end subroutine s_finalize_ibm_module

end module m_ibm
