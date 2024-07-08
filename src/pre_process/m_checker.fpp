!>
!!@file m_checker.f90
!!@brief Contains module m_checker

!> @brief The purpose of the module is to check for compatible input files
module m_checker

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

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
        logical :: skip_check !< Flag to skip the check when iterating over
        !! x, y, and z directions, for special treatment of cylindrical coordinates

        if ((.not. old_grid) .and. old_ic) then
            call s_mpi_abort('old_ic can only be enabled with old_grid enabled. '// &
                             'Exiting ...')
        end if

        if (old_grid) then
            if (t_step_old == dflt_int) then
                call s_mpi_abort('old_grid is enabled, but t_step_old not set. '// &
                                 'Exiting ...')
            elseif ((.not. f_is_default(x_domain%beg)) .or. (.not. f_is_default(x_domain%end)) &
                    .or. &
                    (.not. f_is_default(y_domain%beg)) .or. (.not. f_is_default(y_domain%end)) &
                    .or. &
                    (.not. f_is_default(z_domain%beg)) .or. (.not. f_is_default(z_domain%end))) then
                call s_mpi_abort('x_domain, y_domain, and/or z_domain '// &
                                 'are not supported with old_grid enabled. '// &
                                 'Exiting ...')
            end if
        end if

        #:for DIR, VAR in [('x', 'm'), ('y', 'n'), ('z', 'p')]
            ! For cylindrical coordinates, the y and z directions use a different check
            #:if (DIR == 'y') or (DIR == 'z')
                skip_check = cyl_coord
            #:else
                skip_check = .false.
            #:endif

            if (.not. skip_check) then
                #:for BOUND in ['beg', 'end']
                    if (${VAR}$ == 0) then
                        if (.not. f_is_default((${DIR}$_domain%${BOUND}$))) then
                            call s_mpi_abort('${DIR}$_domain%${BOUND}$ must not '// &
                                             'be set when ${VAR}$ = 0. Exiting ...')
                        end if
                    else ! ${VAR}$ > 0
                        if (old_grid .and. (.not. f_is_default(${DIR}$_domain%${BOUND}$))) then
                            call s_mpi_abort('${DIR}$_domain%${BOUND}$ must not '// &
                                             'be set when ${VAR}$ > 0 and '// &
                                             'old_grid = T. Exiting ...')
                        elseif (.not. old_grid .and. f_is_default(${DIR}$_domain%${BOUND}$)) then
                            call s_mpi_abort('${DIR}$_domain%${BOUND}$ must be '// &
                                             'set when ${VAR}$ > 0 and '// &
                                             'old_grid = F. Exiting ...')
                        elseif (${DIR}$_domain%beg >= ${DIR}$_domain%end) then
                            call s_mpi_abort('${DIR}$_domain%beg must be less '// &
                                             'than ${DIR}$_domain%end when '// &
                                             'both are set. Exiting ...')
                        end if
                    end if
                #:endfor
            end if

        #:endfor

        ! Check for y and z directions for cylindrical coordinates
        if (cyl_coord) then
            if (n == 0) then
                call s_mpi_abort('n must be positive for cylindrical '// &
                                 'coordinates. Exiting ...')
            elseif (f_is_default(y_domain%beg) .or. f_is_default(y_domain%end)) then
                call s_mpi_abort('y_domain%beg and y_domain%end '// &
                                 'must be set for n = 0 '// &
                                 '(2D cylindrical coordinates). Exiting ...')
            elseif (y_domain%beg /= 0d0 .or. y_domain%end <= 0d0) then
                call s_mpi_abort('y_domain%beg must be 0 and y_domain%end '// &
                                 'must be positive for cylindrical '// &
                                 'coordinates. Exiting ...')
            end if

            if (p == 0) then
                if ((.not. f_is_default(z_domain%beg)) &
                    .or. &
                    (.not. f_is_default(z_domain%end))) then
                    call s_mpi_abort('z_domain%beg and z_domain%end '// &
                                     'are not supported for p = 0 '// &
                                     '(2D cylindrical coordinates). Exiting ...')
                end if
            else if (p > 0) then
                if (z_domain%beg /= 0d0 .or. z_domain%end /= 2d0*pi) then
                    call s_mpi_abort('z_domain%beg must be 0 and z_domain%end '// &
                                     'must be 2*pi for 3D cylindrical '// &
                                     'coordinates. Exiting ...')
                end if
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

        ! Constraints specific to stretch_y
        if (stretch_y .and. n == 0) then
            call s_mpi_abort('n must be positive if stretch_y = T. Exiting ...')
        end if

        ! Constraints specific to stretch_z
        if (stretch_z) then
            if (cyl_coord) then
                call s_mpi_abort('stretch_z is not supported for '// &
                                 'cylindrical coordinates. Exiting ...')
            elseif (p == 0) then
                call s_mpi_abort('p must be positive if stretch_z = T. '// &
                                 'Exiting ...')
            end if
        end if

        ! Common checks for all directions (stretch_x, stretch_y, and stretch_z)
        #:for X in ['x', 'y', 'z']
            if (stretch_${X}$) then
                if (old_grid) then
                    call s_mpi_abort('old_grid and stretch_${X}$ are '// &
                                     'incompatible. Exiting ...')
                elseif (f_is_default(a_${X}$)) then
                    call s_mpi_abort('a_${X}$ must be set with stretch_${X}$ '// &
                                     'enabled. Exiting ...')
                elseif (f_is_default(${X}$_a)) then
                    call s_mpi_abort('${X}$_a must be set with stretch_${X}$ '// &
                                     'enabled. Exiting ...')
                elseif (f_is_default(${X}$_b)) then
                    call s_mpi_abort('${X}$_b must be set with stretch_${X}$ '// &
                                     'enabled. Exiting ...')
                elseif (${X}$_a >= ${X}$_b) then
                    call s_mpi_abort('${X}$_a must be less than ${X}$_b with '// &
                                     'stretch_${X}$ enabled. Exiting ...')
                end if
                #:for BOUND in ['beg', 'end']
                    ! Note: `!&` is used to prevent fprettify errors
                    if ((a_${X}$ + log(cosh(a_${X}$*(${X}$_domain%${BOUND}$ - ${X}$_a))) & !&
                                 + log(cosh(a_${X}$*(${X}$_domain%${BOUND}$ - ${X}$_b))) & !&
                                 - 2d0*log(cosh(0.5d0*a_${X}$*(${X}$_b - ${X}$_a)))) / a_${X}$ <= 0d0) then !&
                        call s_mpi_abort('${X}$_domain%${BOUND}$ is too close '// &
                                         'to ${X}$_a and ${X}$_b for the given '// &
                                         'a_${X}$. Exiting ...')
                    end if
                #:endfor
            end if
        #:endfor
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
            (perturb_flow_fluid == dflt_int .or. f_is_default(perturb_flow_mag))) then
            call s_mpi_abort('perturb_flow_fluid and perturb_flow_mag '// &
                             'must be set with perturb_flow = T. Exiting ...')
        elseif ((.not. perturb_flow) &
                .and. &
                (perturb_flow_fluid /= dflt_int .or. (.not. f_is_default(perturb_flow_mag)))) then
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
        elseif ((.not. perturb_sph) .and. (.not. f_all_default(fluid_rho))) then
            call s_mpi_abort('fluid_rho must not be set with perturb_sph = F. '// &
                             'Exiting ...')
        end if

        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (perturb_sph .and. f_is_default(fluid_rho(i))) then
                call s_mpi_abort('fluid_rho('//trim(iStr)//') must be set '// &
                                 'if perturb_sph = T. Exiting ...')
            end if
        end do
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
