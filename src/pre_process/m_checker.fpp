!>
!!@file m_checker.f90
!!@brief Contains module m_checker

!> @brief The purpose of the module is to check for compatible input files
module m_checker

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file.
        !! Used by the pre_process stage
    subroutine s_check_inputs

        call s_check_parallel_io
        call s_check_inputs_restart
        call s_check_inputs_grid_stretching
        call s_check_inputs_qbmm_and_polydisperse
        call s_check_inputs_perturb_density
        call s_check_inputs_misc

    end subroutine s_check_inputs

    !> Checks if mpi is enabled with parallel_io
    subroutine s_check_parallel_io
#ifndef MFC_MPI
        if (parallel_io) then
            print '(A)', 'MFC built with --no-mpi requires parallel_io=F. '// &
                'Exiting ...'
            call s_mpi_abort()
        end if
#endif
    end subroutine s_check_parallel_io

    !> Checks constraints on the restart parameters
        !! (old_grid, old_ic, etc.)
    subroutine s_check_inputs_restart
        if ((.not. old_grid) .and. old_ic) then
            call s_mpi_abort('old_ic cannot be enabled with old_grid disabled. '// &
                             'Exiting ...')

        elseif ((old_grid .or. old_ic) .and. t_step_old == dflt_int) then
            call s_mpi_abort('old_grid or old_ic enabled, but t_step_old not set. '// &
                             'Exiting ...')

            ! Constraints on domain boundaries locations in the x-direction
        elseif ((old_grid .and. x_domain%beg /= dflt_real) &
                .or. &
                ((.not. old_grid) .and. &
                 x_domain%beg == dflt_real)) then
            call s_mpi_abort('old_grid is not enabled but x_domain%beg is '// &
                             'not set. Exiting ...')
        elseif ((old_grid .and. x_domain%end /= dflt_real) &
                .or. &
                ((.not. old_grid) .and. &
                 x_domain%end == dflt_real)) then
            call s_mpi_abort('old_grid is not enabled but x_domain%end is '// &
                             'not set. Exiting ...')
        elseif ((.not. old_grid) &
                .and. &
                x_domain%beg >= x_domain%end) then
            call s_mpi_abort('x_domain%beg must be less than x_domain%end. '// &
                             'Exiting ...')
        end if

        if (cyl_coord) then ! Cartesian coordinates

            ! in case restart of a simulation
            if (old_grid .and. old_ic) then
                ! checking of there is any input to the domains
                if ((x_domain%beg /= dflt_real .or. x_domain%end /= dflt_real) &
                    .or. &
                    (y_domain%beg /= dflt_real .or. y_domain%end /= dflt_real) &
                    .or. &
                    (y_domain%beg /= dflt_real .or. y_domain%end /= dflt_real)) then
                    call s_mpi_abort('x_domain, y_domain, and/or z_domain '// &
                                     'are not supported in restart mode. '// &
                                     '(old_grid = T and old_ic = T). Exiting ...')
                elseif (m == dflt_int .or. n == dflt_int .or. p == dflt_int) then
                    call s_mpi_abort('m, n, and p must be set in restart mode. '// &
                                     '(old_grid = T and old_ic = T). Exiting ...')
                end if
                ! in case it is NOT restart
                ! Constraints on domain boundaries for cylindrical coordinates
            elseif (n == 0) then
                call s_mpi_abort('n must be positive for cylindrical coordinates. '// &
                                 'Exiting ...')
            elseif (y_domain%beg /= 0d0 &
                    .or. &
                    y_domain%end == dflt_real &
                    .or. &
                    y_domain%end < 0d0 &
                    .or. &
                    y_domain%beg >= y_domain%end) then
                call s_mpi_abort('y_domain%beg must be 0 and '// &
                                 'y_domain%end must be positive and '// &
                                 'greater than y_domain%beg for cylindrical '// &
                                 'coordinates. Exiting ...')
            elseif ((p == 0 .and. z_domain%beg /= dflt_real) &
                    .or. &
                    (p == 0 .and. z_domain%end /= dflt_real)) then
                call s_mpi_abort('z_domain%beg and z_domain%end '// &
                                 'are not supported for p = 0. Exiting ...')
            elseif (p > 0 .and. (z_domain%beg /= 0d0 &
                                 .or. &
                                 z_domain%end /= 2d0*pi)) then
                call s_mpi_abort('z_domain%beg must be 0 and '// &
                                 'z_domain%end must be 2*pi for 3D cylindrical '// &
                                 'coordinates. Exiting ...')
            end if

        else
            ! Constraints on domain boundaries locations in the y-direction
            if ((n == 0 .and. y_domain%beg /= dflt_real) .or. &
                (n > 0 .and. ((old_grid .and. y_domain%beg /= dflt_real) .or. &
                              (.not. old_grid .and. y_domain%beg == dflt_real)))) then
                call s_mpi_abort('y_domain%beg must not be set '// &
                                 'when n = 0 or when n > 0 and old_grid = F, and '// &
                                 'must be set otherwise. Exiting ...')
            elseif ((n == 0 .and. y_domain%end /= dflt_real) .or. &
                    (n > 0 .and. ((old_grid .and. y_domain%end /= dflt_real) .or. &
                                  (.not. old_grid .and. y_domain%end == dflt_real)))) then
                call s_mpi_abort('y_domain%end must not be set '// &
                                 'when n = 0 or when n > 0 and old_grid = F, and '// &
                                 'must be set otherwise. Exiting ...')
            elseif (n > 0 .and. .not. old_grid .and. y_domain%beg >= y_domain%end) then
                call s_mpi_abort('y_domain%beg must be less than y_domain%end '// &
                                 'when both are set. '// &
                                 'Exiting ...')

                ! Constraints on domain boundaries locations in the z-direction
            elseif ((p == 0 .and. z_domain%beg /= dflt_real) .or. &
                    (p > 0 .and. ((old_grid .and. z_domain%beg /= dflt_real) .or. &
                                  (.not. old_grid .and. z_domain%beg == dflt_real)))) then
                call s_mpi_abort('z_domain%beg must not be set '// &
                                 'when p = 0 or when p > 0 and old_grid = F, and '// &
                                 'must be set otherwise. Exiting ...')
            elseif ((p == 0 .and. z_domain%end /= dflt_real) .or. &
                    (p > 0 .and. ((old_grid .and. z_domain%end /= dflt_real) .or. &
                                  (.not. old_grid .and. z_domain%end == dflt_real)))) then
                call s_mpi_abort('z_domain%end must not be set '// &
                                 'when p = 0 or when p > 0 and old_grid = F, and '// &
                                 'must be set otherwise. Exiting ...')
            elseif (p > 0 .and. .not. old_grid .and. z_domain%beg >= z_domain%end) then
                call s_mpi_abort('z_domain%beg must be less than z_domain%end '// &
                                 'when both are set. '// &
                                 'Exiting ...')
            end if
        end if

        if (num_patches < 0 .or. &
            (num_patches == 0 .and. t_step_old == dflt_int)) then
            call s_mpi_abort('num_patches must be non-negative for the '// &
                             'non-restart case. Exiting ...')
        end if

    end subroutine s_check_inputs_restart

    !> Checks constraints on grid stretching parameters
        !! (loops_x[y,z], stretch_x[y,z], etc.)
    subroutine s_check_inputs_grid_stretching
        ! Constraints on loops for grid stretching
        if (loops_z < 1) then
            call s_mpi_abort('loops_z must be positive. Exiting ...')
        elseif (loops_y < 1) then
            call s_mpi_abort('loops_y must be positive. Exiting ...')
        end if

        ! Constraints on the grid stretching in the x-direction
        if (stretch_x) then
            if (old_grid) then
                call s_mpi_abort('old_grid and stretch_x are incompatible. '// &
                                 'Exiting ...')
            elseif (a_x == dflt_real) then
                call s_mpi_abort('a_x must be set if stretch_x = T. Exiting ...')
            elseif (x_a == dflt_real) then
                call s_mpi_abort('x_a must be set if stretch_x = T. Exiting ...')
            elseif (x_b == dflt_real) then
                call s_mpi_abort('x_b must be set if stretch_x = T. Exiting ...')
            elseif (x_a >= x_b) then
                call s_mpi_abort('x_a must be less than x_b if stretch_x = T. '// &
                                 'Exiting ...')
            elseif ((a_x + log(cosh(a_x*(x_domain%beg - x_a))) &
                     + log(cosh(a_x*(x_domain%beg - x_b))) &
                     - 2d0*log(cosh(0.5d0*a_x*(x_b - x_a))))/a_x <= 0d0) then
                call s_mpi_abort('x_domain%beg is too close to x_a and x_b '// &
                                 'for the given a_x. Exiting ...')
            elseif ((a_x + log(cosh(a_x*(x_domain%end - x_a))) &
                     + log(cosh(a_x*(x_domain%end - x_b))) &
                     - 2d0*log(cosh(0.5d0*a_x*(x_b - x_a))))/a_x <= 0d0) then
                call s_mpi_abort('x_domain%end is too close to x_a and x_b '// &
                                 'for the given a_x. Exiting ...')
            end if
        end if

        ! Constraints on the grid stretching in the y-direction
        if (stretch_y) then
            if (old_grid) then
                call s_mpi_abort('old_grid and stretch_y are incompatible. '// &
                                 'Exiting ...')
            elseif (n == 0) then
                call s_mpi_abort('n must be positive if stretch_y = T. Exiting ...')
            elseif (a_y == dflt_real) then
                call s_mpi_abort('a_y must be set if stretch_y = T. Exiting ...')
            elseif (y_a == dflt_real) then
                call s_mpi_abort('y_a must be set if stretch_y = T. Exiting ...')
            elseif (y_b == dflt_real) then
                call s_mpi_abort('y_b must be set if stretch_y = T. Exiting ...')
            elseif (y_a >= y_b) then
                call s_mpi_abort('y_a must be less than y_b if stretch_y = T. '// &
                                 'Exiting ...')
            elseif ((a_y + log(cosh(a_y*(y_domain%beg - y_a))) &
                     + log(cosh(a_y*(y_domain%beg - y_b))) &
                     - 2d0*log(cosh(0.5d0*a_y*(y_b - y_a))))/a_y <= 0d0) then
                call s_mpi_abort('y_domain%beg is too close to y_a and y_b '// &
                                 'for the given a_y. Exiting ...')
            elseif ((a_y + log(cosh(a_y*(y_domain%end - y_a))) &
                     + log(cosh(a_y*(y_domain%end - y_b))) &
                     - 2d0*log(cosh(0.5d0*a_y*(y_b - y_a))))/a_y <= 0d0) then
                call s_mpi_abort('y_domain%end is too close to y_a and y_b '// &
                                 'for the given a_y. Exiting ...')
            end if
        end if

        ! Constraints on the grid stretching in the z-direction
        if (stretch_z) then
            if (old_grid) then
                call s_mpi_abort('old_grid and stretch_z are incompatible. '// &
                                 'Exiting ...')
            elseif (cyl_coord) then
                call s_mpi_abort('stretch_z is not supported for cylindrical '// &
                                 'coordinates. Exiting ...')
            elseif (p == 0) then
                call s_mpi_abort('p must be positive if stretch_z = T. Exiting ...')
            elseif (a_z == dflt_real) then
                call s_mpi_abort('a_z must be set if stretch_z = T. Exiting ...')
            elseif (z_a == dflt_real) then
                call s_mpi_abort('z_a must be set if stretch_z = T. Exiting ...')
            elseif (z_b == dflt_real) then
                call s_mpi_abort('z_b must be set if stretch_z = T. Exiting ...')
            elseif (z_a >= z_b) then
                call s_mpi_abort('z_a must be less than z_b if stretch_z = T. '// &
                                 'Exiting ...')
            elseif ((a_z + log(cosh(a_z*(z_domain%beg - z_a))) &
                     + log(cosh(a_z*(z_domain%beg - z_b))) &
                     - 2d0*log(cosh(0.5d0*a_z*(z_b - z_a))))/a_z <= 0d0) then
                call s_mpi_abort('z_domain%beg is too close to z_a and z_b '// &
                                 'for the given a_z. Exiting ...')
            elseif ((a_z + log(cosh(a_z*(z_domain%end - z_a))) &
                     + log(cosh(a_z*(z_domain%end - z_b))) &
                     - 2d0*log(cosh(0.5d0*a_z*(z_b - z_a))))/a_z <= 0d0) then
                call s_mpi_abort('z_domain%end is too close to z_a and z_b '// &
                                 'for the given a_z. Exiting ...')
            end if
        end if
    end subroutine s_check_inputs_grid_stretching

    !> Checks constraints on the QBMM and polydisperse bubble parameters
        !! (qbmm, polydisperse, dist_type, rhoRV, and R0_type)
    subroutine s_check_inputs_qbmm_and_polydisperse
        if (qbmm .and. dist_type == dflt_int) then
            call s_mpi_abort('dist_type must be set if using QBMM. Exiting ...')
        else if (qbmm .and. (dist_type /= 1) .and. rhoRV > 0d0) then
            call s_mpi_abort('rhoRV cannot be used with dist_type != 1. Exiting ...')
        else if (polydisperse .and. R0_type == dflt_int) then
            call s_mpi_abort('R0 type must be set if using Polydisperse. Exiting ...')
        end if
    end subroutine s_check_inputs_qbmm_and_polydisperse

    !> Checks constraints on initial partial density perturbation
        !! (perturb_flow, perturb_flow_fluid, perturb_flow_mag, perturb_sph,
        !! perturb_sph_fluid, and fluid_rho)
    subroutine s_check_inputs_perturb_density
        character(len=5) :: iStr !< for int to string conversion
        integer :: i

        if (perturb_flow &
            .and. &
            (perturb_flow_fluid == dflt_int .or. perturb_flow_mag == dflt_real)) then
            call s_mpi_abort('perturb_flow_fluid and perturb_flow_mag '// &
                             'must be set with perturb_flow = T. Exiting ...')
        elseif ((.not. perturb_flow) &
                .and. &
                (perturb_flow_fluid /= dflt_int .or. perturb_flow_mag /= dflt_real)) then
            call s_mpi_abort('perturb_flow_fluid and perturb_flow_mag '// &
                             'must not be set with perturb_flow = F. Exiting ...')
        elseif ((perturb_flow_fluid > num_fluids) &
                .or. &
                (perturb_flow_fluid < 0 .and. perturb_flow_fluid /= dflt_int)) then
            call s_mpi_abort('perturb_flow_fluid must be between 0 and '// &
                             'num_fluids. Exiting ...')
        elseif (perturb_sph .and. perturb_sph_fluid == dflt_int) then
            call s_mpi_abort('perturb_sph_fluid must be set with perturb_sph = T. '// &
                             'Exiting ...')
        elseif (.not. perturb_sph .and. perturb_sph_fluid /= dflt_int) then
            call s_mpi_abort('perturb_sph_fluid must not be set with perturb_sph = F. '// &
                             'Exiting ...')
        elseif ((perturb_sph_fluid > num_fluids) &
                .or. &
                (perturb_sph_fluid < 0 .and. perturb_sph_fluid /= dflt_int)) then
            call s_mpi_abort('perturb_sph_fluid must be between 0 and '// &
                             'num_fluids. Exiting ...')
        elseif ((any(fluid_rho /= dflt_real)) .and. (.not. perturb_sph)) then
            call s_mpi_abort('fluid_rho must not be set with perturb_sph = F. '// &
                             'Exiting ...')
        elseif (perturb_sph) then
            do i = 1, num_fluids
                call s_int_to_str(i, iStr)
                if (fluid_rho(i) == dflt_real) then
                    call s_mpi_abort('fluid_rho('//trim(iStr)//') must be set '// &
                                     'if perturb_sph = T. Exiting ...')
                end if
            end do
        end if
    end subroutine s_check_inputs_perturb_density

    !> Checks miscellaneous constraints
        !! (vel_profile and instability_wave)
    subroutine s_check_inputs_misc
        ! Hypertangent velocity profile
        if (vel_profile .and. (n == 0)) then
            call s_mpi_abort('vel_profile requires n > 0. Exiting ...')
        end if
        ! Instability wave
        if (instability_wave .and. (n == 0)) then
            call s_mpi_abort('instability_wave requires n > 0. Exiting ...')
        end if
    end subroutine s_check_inputs_misc

end module m_checker
