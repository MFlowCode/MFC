!>
!! @file m_ibm.fpp
!! @brief Contains module m_ibm

#:include 'macros.fpp'

!> @brief This module is used to handle all operations related to immersed
!!              boundary methods (IBMs)
module m_ibm
    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper

    use m_compute_levelset

    ! ==========================================================================

    implicit none

    private :: s_compute_image_points, &
               s_compute_interpolation_coeffs, &
               s_interpolate_image_point, &
               s_compute_levelset, &
               s_find_ghost_points, &
               s_find_num_ghost_points
    ; public :: s_initialize_ibm_module, &
 s_ibm_setup, &
 s_ibm_correct_state, &
 s_finalize_ibm_module

    type(integer_field), public :: ib_markers
!$acc declare create(ib_markers)

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), levelset)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :, :), levelset_norm)
    @:CRAY_DECLARE_GLOBAL(type(ghost_point), dimension(:), ghost_points)
    @:CRAY_DECLARE_GLOBAL(type(ghost_point), dimension(:), inner_points)

    !$acc declare link(levelset, levelset_norm, ghost_points, inner_points)
#else

    !! Marker for solid cells. 0 if liquid, the patch id of its IB if solid
    real(kind(0d0)), dimension(:, :, :, :), allocatable :: levelset
    !! Matrix of distance to IB
    real(kind(0d0)), dimension(:, :, :, :, :), allocatable :: levelset_norm
    !! Matrix of normal vector to IB
    type(ghost_point), dimension(:), allocatable :: ghost_points
    type(ghost_point), dimension(:), allocatable :: inner_points
    !! Matrix of normal vector to IB

    !$acc declare create(levelset, levelset_norm, ghost_points, inner_points)
#endif

    integer :: gp_layers !< Number of ghost point layers
    integer :: num_gps !< Number of ghost points
    integer :: num_inner_gps !< Number of ghost points
    !$acc declare create(gp_layers, num_gps, num_inner_gps)

contains

    !>  Initialize IBM module
    subroutine s_initialize_ibm_module

        gp_layers = 3

        if (p > 0) then
            @:ALLOCATE(ib_markers%sf(-gp_layers:m+gp_layers, &
                -gp_layers:n+gp_layers, -gp_layers:p+gp_layers))
        else
            @:ALLOCATE(ib_markers%sf(-gp_layers:m+gp_layers, &
                -gp_layers:n+gp_layers, 0:0))
        end if
        @:ACC_SETUP_SFs(ib_markers)

        ! @:ALLOCATE(ib_markers%sf(0:m, 0:n, 0:p))
        @:ALLOCATE_GLOBAL(levelset(0:m, 0:n, 0:p, num_ibs))
        @:ALLOCATE_GLOBAL(levelset_norm(0:m, 0:n, 0:p, num_ibs, 3))

        !$acc enter data copyin(gp_layers, num_gps, num_inner_gps)

    end subroutine s_initialize_ibm_module

    subroutine s_ibm_setup

        integer :: i, j, k

        !$acc update device(ib_markers%sf)

        ! Get neighboring IB variables from other processors
        call s_mpi_sendrecv_ib_buffers(ib_markers, gp_layers)

        call s_find_num_ghost_points()

        !$acc update device(num_gps, num_inner_gps)
        @:ALLOCATE_GLOBAL(ghost_points(num_gps))
        @:ALLOCATE_GLOBAL(inner_points(num_inner_gps))

        !$acc enter data copyin(ghost_points, inner_points)

        call s_find_ghost_points(ghost_points, inner_points)
        !$acc update device(ghost_points, inner_points)

        call s_compute_levelset(levelset, levelset_norm)
        !$acc update device(levelset, levelset_norm)

        call s_compute_image_points(ghost_points, levelset, levelset_norm)
        !$acc update device(ghost_points)

        call s_compute_interpolation_coeffs(ghost_points)
        !$acc update device(ghost_points)

    end subroutine s_ibm_setup

    !>  Subroutine that updates the conservative variables at the ghost points
        !!  @param q_cons_vf Conservative variables
        !!  @param q_prim_vf Primitive variables
    subroutine s_ibm_correct_state(q_cons_vf, q_prim_vf, pb, mv)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_cons_vf !< Conservative Variables

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_prim_vf !< Primitive Variables

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), optional, intent(inout) :: pb, mv

        integer :: i, j, k, l, q, r!< Iterator variables
        integer :: patch_id !< Patch ID of ghost point
        real(kind(0d0)) :: rho, gamma, pi_inf, dyn_pres !< Mixture variables
        real(kind(0d0)), dimension(2) :: Re_K
        real(kind(0d0)) :: G_K
        real(kind(0d0)) :: qv_K
        real(kind(0d0)), dimension(num_fluids) :: Gs

        real(kind(0d0)) :: pres_IP, coeff
        real(kind(0d0)), dimension(3) :: vel_IP, vel_norm_IP
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_IP, alpha_IP
        real(kind(0d0)), dimension(nb) :: r_IP, v_IP, pb_IP, mv_IP
        real(kind(0d0)), dimension(nb*nmom) :: nmom_IP
        real(kind(0d0)), dimension(nb*nnode) :: presb_IP, massv_IP
        !! Primitive variables at the image point associated with a ghost point,
        !! interpolated from surrounding fluid cells.

        real(kind(0d0)), dimension(3) :: norm !< Normal vector from GP to IP
        real(kind(0d0)), dimension(3) :: physical_loc !< Physical loc of GP
        real(kind(0d0)), dimension(3) :: vel_g !< Velocity of GP

        real(kind(0d0)) :: nbub
        real(kind(0d0)) :: buf
        type(ghost_point) :: gp
        type(ghost_point) :: innerp

        !$acc parallel loop gang vector private(physical_loc, dyn_pres, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, vel_g, vel_norm_IP, r_IP, v_IP, pb_IP, mv_IP, nmom_IP, presb_IP, massv_IP, rho, gamma, pi_inf, Re_K, G_K, Gs, gp, innerp, norm, buf, j, k, l, q, coeff)
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
                physical_loc = [x_cc(j), y_cc(k), 0d0]
            end if

            !Interpolate primitive variables at image point associated w/ GP
            if (bubbles .and. .not. qbmm) then
                call s_interpolate_image_point(q_prim_vf, gp, &
                                               alpha_rho_IP, alpha_IP, pres_IP, vel_IP, &
                                               r_IP, v_IP, pb_IP, mv_IP)
            else if (qbmm .and. polytropic) then
                call s_interpolate_image_point(q_prim_vf, gp, &
                                               alpha_rho_IP, alpha_IP, pres_IP, vel_IP, &
                                               r_IP, v_IP, pb_IP, mv_IP, nmom_IP)
            else if (qbmm .and. .not. polytropic) then
                call s_interpolate_image_point(q_prim_vf, gp, &
                                               alpha_rho_IP, alpha_IP, pres_IP, vel_IP, &
                                               r_IP, v_IP, pb_IP, mv_IP, nmom_IP, pb, mv, presb_IP, massv_IP)
            else
                call s_interpolate_image_point(q_prim_vf, gp, &
                                               alpha_rho_IP, alpha_IP, pres_IP, vel_IP)
            end if

            dyn_pres = 0d0

            ! Set q_prim_vf params at GP so that mixture vars calculated properly
            !$acc loop seq
            do q = 1, num_fluids
                q_prim_vf(q)%sf(j, k, l) = alpha_rho_IP(q)
                q_prim_vf(advxb + q - 1)%sf(j, k, l) = alpha_IP(q)
            end do

            if (model_eqns /= 4) then
                ! If in simulation, use acc mixture subroutines
                if (hypoelasticity) then
                    call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv_K, alpha_IP, &
                                                                    alpha_rho_IP, Re_K, j, k, l, G_K, Gs)
                else if (bubbles) then
                    call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, qv_K, alpha_IP, &
                                                                            alpha_rho_IP, Re_K, j, k, l)
                else
                    call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv_K, alpha_IP, &
                                                                    alpha_rho_IP, Re_K, j, k, l)
                end if
            end if

            ! Calculate velocity of ghost cell
            if (gp%slip) then
                norm = gp%ip_loc - physical_loc !
                buf = sqrt(sum(norm**2))
                norm = norm/buf
                vel_norm_IP = sum(vel_IP*norm)*norm
                vel_g = vel_IP - vel_norm_IP
            else
                vel_g = 0d0
            end if

            ! Set momentum
            !$acc loop seq
            do q = momxb, momxe
                q_cons_vf(q)%sf(j, k, l) = rho*vel_g(q - momxb + 1)
                dyn_pres = dyn_pres + q_cons_vf(q)%sf(j, k, l)* &
                           vel_g(q - momxb + 1)/2d0
            end do

            ! Set continuity and adv vars
            !$acc loop seq
            do q = 1, num_fluids
                q_cons_vf(q)%sf(j, k, l) = alpha_rho_IP(q)
                q_cons_vf(advxb + q - 1)%sf(j, k, l) = alpha_IP(q)
            end do

            ! Set Energy
            if (bubbles) then
                q_cons_vf(E_idx)%sf(j, k, l) = (1 - alpha_IP(1))*(gamma*pres_IP + pi_inf + dyn_pres)
            else
                q_cons_vf(E_idx)%sf(j, k, l) = gamma*pres_IP + pi_inf + dyn_pres
            end if

            ! Set bubble vars
            if (bubbles .and. .not. qbmm) then
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
                            pb(j, k, l, r, q) = presb_IP((q - 1)*nnode + r)
                            mv(j, k, l, r, q) = massv_IP((q - 1)*nnode + r)
                        end do
                    end do
                end if
            end if

            if (model_eqns == 3) then
                !$acc loop seq
                do q = intxb, intxe
                    q_cons_vf(q)%sf(j, k, l) = alpha_IP(q - intxb + 1)*(gammas(q - intxb + 1)*pres_IP &
                                                                        + pi_infs(q - intxb + 1))
                end do
            end if
        end do

        !Correct the state of the inner points in IBs
        !$acc parallel loop gang vector private(physical_loc, dyn_pres, alpha_rho_IP, alpha_IP, vel_g, rho, gamma, pi_inf, Re_K, innerp, j, k, l, q)
        do i = 1, num_inner_gps

            vel_g = 0d0
            innerp = inner_points(i)
            j = innerp%loc(1)
            k = innerp%loc(2)
            l = innerp%loc(3)
            patch_id = inner_points(i)%ib_patch_id

            ! Calculate physical location of GP
            if (p > 0) then
                physical_loc = [x_cc(j), y_cc(k), z_cc(l)]
            else
                physical_loc = [x_cc(j), y_cc(k), 0d0]
            end if

            !$acc loop seq
            do q = 1, num_fluids
                q_prim_vf(q)%sf(j, k, l) = alpha_rho_IP(q)
                q_prim_vf(advxb + q - 1)%sf(j, k, l) = alpha_IP(q)
            end do

            call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv_K, alpha_IP, &
                                                            alpha_rho_IP, Re_K, j, k, l)

            dyn_pres = 0d0

            !$acc loop seq
            do q = momxb, momxe
                q_cons_vf(q)%sf(j, k, l) = rho*vel_g(q - momxb + 1)
                dyn_pres = dyn_pres + q_cons_vf(q)%sf(j, k, l)* &
                           vel_g(q - momxb + 1)/2d0
            end do
        end do

    end subroutine s_ibm_correct_state

    !>  Subroutine that computes that bubble wall pressure for Gilmore bubbles
    subroutine s_compute_image_points(ghost_points, levelset, levelset_norm)

        type(ghost_point), dimension(num_gps), intent(inout) :: ghost_points
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(in) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(in) :: levelset_norm

        real(kind(0d0)) :: dist
        real(kind(0d0)), dimension(3) :: norm
        real(kind(0d0)), dimension(3) :: physical_loc
        real(kind(0d0)) :: temp_loc
        real(kind(0d0)), pointer, dimension(:) :: s_cc => null()
        integer :: bound
        type(ghost_point) :: gp

        integer :: q, dim !< Iterator variables
        integer :: i, j, k !< Location indexes
        integer :: patch_id !< IB Patch ID
        integer :: dir
        integer :: index

        do q = 1, num_gps
            gp = ghost_points(q)
            i = gp%loc(1)
            j = gp%loc(2)
            k = gp%loc(3)

            ! Calculate physical location of ghost point
            if (p > 0) then
                physical_loc = [x_cc(i), y_cc(j), z_cc(k)]
            else
                physical_loc = [x_cc(i), y_cc(j), 0d0]
            end if

            ! Calculate and store the precise location of the image point
            patch_id = gp%ib_patch_id
            dist = abs(levelset(i, j, k, patch_id))
            norm(:) = levelset_norm(i, j, k, patch_id, :)
            ghost_points(q)%ip_loc(:) = physical_loc(:) + 2*dist*norm(:)

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

                if (norm(dim) == 0) then
                    ghost_points(q)%ip_grid(dim) = ghost_points(q)%loc(dim)
                else
                    if (norm(dim) > 0) then
                        dir = 1
                    else
                        dir = -1
                    end if

                    index = ghost_points(q)%loc(dim)
                    temp_loc = ghost_points(q)%ip_loc(dim)
                    do while ((temp_loc < s_cc(index) &
                               .or. temp_loc > s_cc(index + 1)) &
                              .and. (index >= 0 .and. index <= bound))
                        index = index + dir
                    end do
                    ghost_points(q)%ip_grid(dim) = index
                    if (ghost_points(q)%DB(dim) == -1) then
                        ghost_points(q)%ip_grid(dim) = ghost_points(q)%loc(dim) + 1
                    else if (ghost_points(q)%DB(dim) == 1) then
                        ghost_points(q)%ip_grid(dim) = ghost_points(q)%loc(dim) - 1
                    end if
                end if
            end do

            ! print *, "GP Loc: ", ghost_points(q)%loc(:)
            ! print *, "Norm: ", norm(:)
            ! print *, "Dist: ", abs(dist)
            ! print *, "IP Loc: ", ghost_points(q)%ip_grid(:)
            ! print *, "------"
        end do

#if 0
        if (proc_rank == 0) then

            open (unit=10, file=trim(case_dir)//'/gp.txt', status='replace')
            do i = 1, num_gps
                write (10, '(3F36.12)') x_cc(ghost_points(i)%loc(1)), y_cc(ghost_points(i)%loc(2))
            end do
            close (10)

            open (unit=10, file=trim(case_dir)//'/ip.txt', status='replace')
            do i = 1, num_gps
                write (10, '(3F36.12)') ghost_points(i)%ip_loc(1), ghost_points(i)%ip_loc(2)
            end do
            close (10)

        end if

#endif

    end subroutine s_compute_image_points

    subroutine s_find_num_ghost_points
        integer, dimension(2*gp_layers + 1, 2*gp_layers + 1) &
            :: subsection_2D
        integer, dimension(2*gp_layers + 1, 2*gp_layers + 1, 2*gp_layers + 1) &
            :: subsection_3D
        integer :: i, j, k, l, q !< Iterator variables

        do i = 0, m
            do j = 0, n
                if (p == 0) then
                    if (ib_markers%sf(i, j, 0) /= 0) then
                        subsection_2D = ib_markers%sf( &
                                        i - gp_layers:i + gp_layers, &
                                        j - gp_layers:j + gp_layers, 0)
                        if (any(subsection_2D == 0)) then
                            num_gps = num_gps + 1
                        else
                            num_inner_gps = num_inner_gps + 1
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
                                num_gps = num_gps + 1
                            else
                                num_inner_gps = num_inner_gps + 1
                            end if
                        end if
                    end do
                end if
            end do
        end do

    end subroutine s_find_num_ghost_points

    subroutine s_find_ghost_points(ghost_points, inner_points)

        type(ghost_point), dimension(num_gps), intent(inout) :: ghost_points
        type(ghost_point), dimension(num_inner_gps), intent(inout) :: inner_points

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
                            ghost_points(count)%loc = [i, j, 0]
                            patch_id = ib_markers%sf(i, j, 0)
                            ghost_points(count)%ib_patch_id = &
                                patch_id
                            ghost_points(count)%slip = patch_ib(patch_id)%slip

                            if ((x_cc(i) - dx(i)) < x_domain%beg) then
                                ghost_points(count)%DB(1) = -1
                            else if ((x_cc(i) + dx(i)) > x_domain%end) then
                                ghost_points(count)%DB(1) = 1
                            else
                                ghost_points(count)%DB(1) = 0
                            end if

                            if ((y_cc(j) - dy(j)) < y_domain%beg) then
                                ghost_points(count)%DB(2) = -1
                            else if ((y_cc(j) + dy(j)) > y_domain%end) then
                                ghost_points(count)%DB(2) = 1
                            else
                                ghost_points(count)%DB(2) = 0
                            end if

                            count = count + 1

                        else
                            inner_points(count_i)%loc = [i, j, 0]
                            patch_id = ib_markers%sf(i, j, 0)
                            inner_points(count_i)%ib_patch_id = &
                                patch_id
                            inner_points(count_i)%slip = patch_ib(patch_id)%slip
                            if ((x_cc(i) - dx(i)) < x_domain%beg .or. &
                                (x_cc(i) + dx(i)) > x_domain%end) then
                                ghost_points(count)%DB(1) = 1
                            else
                                ghost_points(count)%DB(1) = 0
                            end if

                            if ((y_cc(j) - dy(j)) < y_domain%beg .or. &
                                (y_cc(j) + dy(j)) > y_domain%end) then
                                ghost_points(count)%DB(2) = 1
                            else
                                ghost_points(count)%DB(2) = 0
                            end if

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
                                ghost_points(count)%loc = [i, j, k]
                                patch_id = ib_markers%sf(i, j, k)
                                ghost_points(count)%ib_patch_id = &
                                    ib_markers%sf(i, j, k)
                                ghost_points(count)%slip = patch_ib(patch_id)%slip

                                if ((x_cc(i) - dx(i)) < x_domain%beg) then
                                    ghost_points(count)%DB(1) = -1
                                else if ((x_cc(i) + dx(i)) > x_domain%end) then
                                    ghost_points(count)%DB(1) = 1
                                else
                                    ghost_points(count)%DB(1) = 0
                                end if

                                if ((y_cc(j) - dy(j)) < y_domain%beg) then
                                    ghost_points(count)%DB(2) = -1
                                else if ((y_cc(j) + dy(j)) > y_domain%end) then
                                    ghost_points(count)%DB(2) = 1
                                else
                                    ghost_points(count)%DB(2) = 0
                                end if

                                if ((z_cc(k) - dz(k)) < z_domain%beg) then
                                    ghost_points(count)%DB(3) = -1
                                else if ((z_cc(k) + dz(k)) > z_domain%end) then
                                    ghost_points(count)%DB(3) = 1
                                else
                                    ghost_points(count)%DB(3) = 0
                                end if

                                count = count + 1
                            else
                                inner_points(count_i)%loc = [i, j, k]
                                patch_id = ib_markers%sf(i, j, k)
                                inner_points(count_i)%ib_patch_id = &
                                    ib_markers%sf(i, j, k)
                                inner_points(count_i)%slip = patch_ib(patch_id)%slip

                                if ((x_cc(i) - dx(i)) < x_domain%beg) then
                                    ghost_points(count)%DB(1) = -1
                                else if ((x_cc(i) + dx(i)) > x_domain%end) then
                                    ghost_points(count)%DB(1) = 1
                                else
                                    ghost_points(count)%DB(1) = 0
                                end if

                                if ((y_cc(j) - dy(j)) < y_domain%beg) then
                                    ghost_points(count)%DB(2) = -1
                                else if ((y_cc(j) + dy(j)) > y_domain%end) then
                                    ghost_points(count)%DB(2) = 1
                                else
                                    ghost_points(count)%DB(2) = 0
                                end if

                                if ((z_cc(k) - dz(k)) < z_domain%beg) then
                                    ghost_points(count)%DB(3) = -1
                                else if ((z_cc(k) + dz(k)) > z_domain%end) then
                                    ghost_points(count)%DB(3) = 1
                                else
                                    ghost_points(count)%DB(3) = 0
                                end if

                                count_i = count_i + 1
                            end if
                        end if
                    end do
                end if
            end do
        end do

    end subroutine s_find_ghost_points

    !>  Function that computes that bubble wall pressure for Gilmore bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    subroutine s_compute_interpolation_coeffs(ghost_points)

        type(ghost_point), dimension(num_gps), intent(inout) :: ghost_points

        real(kind(0d0)), dimension(2, 2, 2) :: dist
        real(kind(0d0)), dimension(2, 2, 2) :: alpha
        real(kind(0d0)), dimension(2, 2, 2) :: interp_coeffs
        real(kind(0d0)) :: buf
        real(kind(0d0)), dimension(2, 2, 2) :: eta
        type(ghost_point) :: gp
        integer :: i, j, k, l, q !< Iterator variables
        integer :: i1, i2, j1, j2, k1, k2 !< Grid indexes
        integer :: patch_id

        ! 2D
        if (p <= 0) then
            do i = 1, num_gps
                gp = ghost_points(i)
                ! Get the interpolation points
                i1 = gp%ip_grid(1); i2 = i1 + 1
                j1 = gp%ip_grid(2); j2 = j1 + 1

                dist = 0d0
                buf = 1d0
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

                interp_coeffs = 0d0

                if (dist(1, 1, 1) <= 1d-16) then
                    interp_coeffs(1, 1, 1) = 1d0
                else if (dist(2, 1, 1) <= 1d-16) then
                    interp_coeffs(2, 1, 1) = 1d0
                else if (dist(1, 2, 1) <= 1d-16) then
                    interp_coeffs(1, 2, 1) = 1d0
                else if (dist(2, 2, 1) <= 1d-16) then
                    interp_coeffs(2, 2, 1) = 1d0
                else
                    eta(:, :, 1) = 1d0/dist(:, :, 1)**2
                    alpha = 1d0
                    patch_id = gp%ib_patch_id
                    if (ib_markers%sf(i1, j1, 0) /= 0) alpha(1, 1, 1) = 0d0
                    if (ib_markers%sf(i2, j1, 0) /= 0) alpha(2, 1, 1) = 0d0
                    if (ib_markers%sf(i1, j2, 0) /= 0) alpha(1, 2, 1) = 0d0
                    if (ib_markers%sf(i2, j2, 0) /= 0) alpha(2, 2, 1) = 0d0
                    buf = sum(alpha(:, :, 1)*eta(:, :, 1))
                    if (buf > 0d0) then
                        interp_coeffs(:, :, 1) = alpha(:, :, 1)*eta(:, :, 1)/buf
                    else
                        buf = sum(eta(:, :, 1))
                        interp_coeffs(:, :, 1) = eta(:, :, 1)/buf
                    end if
                end if

                ghost_points(i)%interp_coeffs = interp_coeffs
            end do

        else
            do i = 1, num_gps
                gp = ghost_points(i)
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
                interp_coeffs = 0d0
                buf = 1d0
                if (dist(1, 1, 1) <= 1d-16) then
                    interp_coeffs(1, 1, 1) = 1d0
                else if (dist(2, 1, 1) <= 1d-16) then
                    interp_coeffs(2, 1, 1) = 1d0
                else if (dist(1, 2, 1) <= 1d-16) then
                    interp_coeffs(1, 2, 1) = 1d0
                else if (dist(2, 2, 1) <= 1d-16) then
                    interp_coeffs(2, 2, 1) = 1d0
                else if (dist(1, 1, 2) <= 1d-16) then
                    interp_coeffs(1, 1, 2) = 1d0
                else if (dist(2, 1, 2) <= 1d-16) then
                    interp_coeffs(2, 1, 2) = 1d0
                else if (dist(1, 2, 2) <= 1d-16) then
                    interp_coeffs(1, 2, 2) = 1d0
                else if (dist(2, 2, 2) <= 1d-16) then
                    interp_coeffs(2, 2, 2) = 1d0
                else
                    eta = 1d0/dist**2
                    alpha = 1d0
                    if (ib_markers%sf(i1, j1, k1) /= 0) alpha(1, 1, 1) = 0d0
                    if (ib_markers%sf(i2, j1, k1) /= 0) alpha(2, 1, 1) = 0d0
                    if (ib_markers%sf(i1, j2, k1) /= 0) alpha(1, 2, 1) = 0d0
                    if (ib_markers%sf(i2, j2, k1) /= 0) alpha(2, 2, 1) = 0d0
                    if (ib_markers%sf(i1, j1, k2) /= 0) alpha(1, 1, 2) = 0d0
                    if (ib_markers%sf(i2, j1, k2) /= 0) alpha(2, 1, 2) = 0d0
                    if (ib_markers%sf(i1, j2, k2) /= 0) alpha(1, 2, 2) = 0d0
                    if (ib_markers%sf(i2, j2, k2) /= 0) alpha(2, 2, 2) = 0d0
                    buf = sum(alpha*eta)
                    if (buf > 0d0) then
                        interp_coeffs = alpha*eta/buf
                    else
                        buf = sum(eta)
                        interp_coeffs = eta/buf
                    end if
                end if

                ghost_points(i)%interp_coeffs = interp_coeffs
            end do
        end if

    end subroutine s_compute_interpolation_coeffs

    subroutine s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, r_IP, v_IP, pb_IP, mv_IP, nmom_IP, pb, mv, presb_IP, massv_IP)
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf !< Primitive Variables
        type(ghost_point), intent(in) :: gp
        real(kind(0d0)), dimension(num_fluids), intent(inout) :: alpha_IP, alpha_rho_IP
        real(kind(0d0)), intent(inout) :: pres_IP
        real(kind(0d0)), dimension(3), intent(inout) :: vel_IP
        real(kind(0d0)), optional, dimension(:), intent(inout) :: r_IP, v_IP, pb_IP, mv_IP
        real(kind(0d0)), optional, dimension(:), intent(inout) :: nmom_IP
        real(kind(0d0)), optional, dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
        real(kind(0d0)), optional, dimension(:), intent(inout) :: presb_IP, massv_IP

        integer :: i, j, k, l, q !< Iterator variables
        integer :: i1, i2, j1, j2, k1, k2 !< Iterator variables
        real(kind(0d0)) :: coeff

        i1 = gp%ip_grid(1); i2 = i1 + 1
        j1 = gp%ip_grid(2); j2 = j1 + 1
        k1 = gp%ip_grid(3); k2 = k1 + 1

        if (p == 0) then
            k1 = 0
            k2 = 0
        end if

        alpha_rho_IP = 0d0
        alpha_IP = 0d0
        pres_IP = 0d0
        vel_IP = 0d0

        if (bubbles) then
            r_IP = 0d0
            v_IP = 0d0
            if (.not. polytropic) then
                mv_IP = 0d0
                pb_IP = 0d0
            end if
        end if

        if (qbmm) then
            nmom_IP = 0d0
            if (.not. polytropic) then
                presb_IP = 0d0
                massv_IP = 0d0
            end if
        end if

        !$acc loop seq
        do i = i1, i2
            !$acc loop seq
            do j = j1, j2
                !$acc loop seq
                do k = k1, k2

                    coeff = gp%interp_coeffs(i - i1 + 1, j - j1 + 1, k - k1 + 1)

                    pres_IP = pres_IP + coeff* &
                              q_prim_vf(E_idx)%sf(i, j, k)

                    !$acc loop seq
                    do q = momxb, momxe
                        vel_IP(q + 1 - momxb) = vel_IP(q + 1 - momxb) + coeff* &
                                                q_prim_vf(q)%sf(i, j, k)
                    end do

                    !$acc loop seq
                    do l = contxb, contxe
                        alpha_rho_IP(l) = alpha_rho_IP(l) + coeff* &
                                          q_prim_vf(l)%sf(i, j, k)
                        alpha_IP(l) = alpha_IP(l) + coeff* &
                                      q_prim_vf(advxb + l - 1)%sf(i, j, k)
                    end do

                    if (bubbles .and. .not. qbmm) then
                        !$acc loop seq
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
                                    presb_IP((q - 1)*nnode + l) = presb_IP((q - 1)*nnode + l) + coeff*pb(i, j, k, l, q)
                                    massv_IP((q - 1)*nnode + l) = massv_IP((q - 1)*nnode + l) + coeff*mv(i, j, k, l, q)
                                end do
                            end do
                        end if

                    end if

                end do
            end do
        end do

    end subroutine s_interpolate_image_point

    !>  Subroutine that computes that bubble wall pressure for Gilmore bubbles
    subroutine s_compute_levelset(levelset, levelset_norm)

        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(inout) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(inout) :: levelset_norm
        integer :: i !< Iterator variables
        integer :: geometry

        do i = 1, num_ibs
            geometry = patch_ib(i)%geometry
            if (geometry == 2) then
                call s_compute_circle_levelset(levelset, levelset_norm, i)
            else if (geometry == 3) then
                call s_compute_rectangle_levelset(levelset, levelset_norm, i)
            else if (geometry == 4) then
                call s_compute_airfoil_levelset(levelset, levelset_norm, i)
            else if (geometry == 8) then
                call s_compute_sphere_levelset(levelset, levelset_norm, i)
            else if (geometry == 10) then
                call s_compute_cylinder_levelset(levelset, levelset_norm, i)
            else if (geometry == 11) then
                call s_compute_3D_airfoil_levelset(levelset, levelset_norm, i)
            end if
        end do

    end subroutine s_compute_levelset

    !>  Subroutine that computes that bubble wall pressure for Gilmore bubbles
    subroutine s_finalize_ibm_module

        @:DEALLOCATE(ib_markers%sf)
        @:DEALLOCATE_GLOBAL(levelset)
        @:DEALLOCATE_GLOBAL(levelset_norm)
    end subroutine s_finalize_ibm_module

end module m_ibm
