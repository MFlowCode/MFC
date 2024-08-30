!>
!! @file m_checker.f90
!! @brief Contains module m_checker

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
        !! Used by the post_process stage
    subroutine s_check_inputs

        call s_check_inputs_output_format
        call s_check_inputs_partial_density
        call s_check_inputs_velocity
        call s_check_inputs_flux_limiter
        call s_check_inputs_volume_fraction
        call s_check_inputs_vorticity
        call s_check_inputs_schlieren
        call s_check_inputs_surface_tension
        call s_check_inputs_no_flow_variables

    end subroutine s_check_inputs

    !> Checks constraints on output format parameters
    subroutine s_check_inputs_output_format
        if (format /= 1 .and. format /= 2) then
            call s_mpi_abort('format must be 1 or 2. Exiting ...')
        elseif (precision /= 1 .and. precision /= 2) then
            call s_mpi_abort('precision must be 1 or 2. Exiting ...')
        end if
    end subroutine s_check_inputs_output_format

    !> Checks constraints on partial density parameters
    subroutine s_check_inputs_partial_density
        character(len=5) :: iStr
        integer :: i

        do i = 1, num_fluids
            if (alpha_rho_wrt(i)) then
                call s_int_to_str(i, iStr)
                if (model_eqns == 1) then
                    call s_mpi_abort('alpha_rho_wrt('//trim(iStr)//') is not '// &
                                     'supported for model_eqns = 1. Exiting ...')
                end if
                if (i > num_fluids) then
                    call s_mpi_abort('Index of alpha_rho_wrt('//trim(iStr)//') '// &
                                     'exceeds the total number of fluids. Exiting ...')
                end if
            end if
        end do
    end subroutine s_check_inputs_partial_density

    !> Checks constraints on momentum parameters
    subroutine s_check_inputs_momentum
        if (n == 0 .and. mom_wrt(2)) then
            call s_mpi_abort('mom_wrt(2) is not supported for n = 0. Exiting ...')
        elseif (p == 0 .and. mom_wrt(3)) then
            call s_mpi_abort('mom_wrt(3) is not supported for p = 0. Exiting ...')
        end if
    end subroutine s_check_inputs_momentum

    !> Checks constraints on velocity parameters
    subroutine s_check_inputs_velocity
        if (n == 0 .and. vel_wrt(2)) then
            call s_mpi_abort('vel_wrt(2) is not supported for n = 0. Exiting ...')
        elseif (p == 0 .and. vel_wrt(3)) then
            call s_mpi_abort('vel_wrt(3) is not supported for p = 0. Exiting ...')
        end if
    end subroutine s_check_inputs_velocity

    !> Checks constraints on flux limiter parameters
    subroutine s_check_inputs_flux_limiter
        if (n == 0 .and. flux_wrt(2)) then
            call s_mpi_abort('flux_wrt(2) is not supported for n = 0. Exiting ...')
        elseif (p == 0 .and. flux_wrt(3)) then
            call s_mpi_abort('flux_wrt(3) is not supported for p = 0. Exiting ...')
        elseif (all(flux_lim /= (/dflt_int, 1, 2, 3, 4, 5, 6, 7/))) then
            call s_mpi_abort('flux_lim must be between 1 and 7. Exiting ...')
        end if
    end subroutine s_check_inputs_flux_limiter

    !> Checks constraints on volume fraction parameters
    subroutine s_check_inputs_volume_fraction
        character(len=5) :: iStr
        integer :: i

        do i = 1, num_fluids
            if (alpha_wrt(i)) then
                call s_int_to_str(i, iStr)
                if (model_eqns == 1) then
                    call s_mpi_abort('alpha_wrt('//trim(iStr)//') is not '// &
                                     'supported for model_eqns = 1. Exiting ...')
                end if
                if (i > num_fluids) then
                    call s_mpi_abort('Index of alpha_wrt('//trim(iStr)//') '// &
                                     'exceeds the total number of fluids. Exiting ...')
                end if
            end if
        end do
    end subroutine s_check_inputs_volume_fraction

    !> Checks constraints on vorticity parameters
    subroutine s_check_inputs_vorticity
        if (n == 0 .and. any(omega_wrt)) then
            call s_mpi_abort('omega_wrt is not supported for n = 0. Exiting ...')
        elseif (p == 0 .and. (omega_wrt(1) .or. omega_wrt(2))) then
            call s_mpi_abort('omega_wrt(1) and omega_wrt(2) are not supported '// &
                             'for p = 0. Exiting ...')
        elseif (any(omega_wrt) .and. fd_order == dflt_int) then
            call s_mpi_abort('fd_order must be set for omega_wrt. Exiting ...')
        end if
    end subroutine s_check_inputs_vorticity

    !> Checks constraints on numerical Schlieren parameters
        !! (schlieren_wrt and schlieren_alpha)
    subroutine s_check_inputs_schlieren
        character(len=5) :: iStr
        integer :: i

        if (n == 0 .and. schlieren_wrt) then
            call s_mpi_abort('schlieren_wrt is not supported for n = 0. Exiting ...')
        elseif (schlieren_wrt .and. fd_order == dflt_int) then
            call s_mpi_abort('fd_order must be set for schlieren_wrt. Exiting ...')
        end if

        do i = 1, num_fluids
            if (.not. f_is_default(schlieren_alpha(i))) then
                call s_int_to_str(i, iStr)
                if (schlieren_alpha(i) <= 0d0) then
                    call s_mpi_abort('schlieren_alpha('//trim(iStr)//') must be '// &
                                     'greater than zero. Exiting ...')
                elseif (i > num_fluids) then
                    call s_mpi_abort('Index of schlieren_alpha('//trim(iStr)//') '// &
                                     'exceeds the total number of fluids. Exiting ...')
                elseif (.not. schlieren_wrt) then
                    call s_mpi_abort('schlieren_alpha('//trim(iStr)//') should '// &
                                     'be set only with schlieren_wrt enabled. Exiting ...')
                end if
            end if
        end do
    end subroutine s_check_inputs_schlieren

    !> Checks constraints on surface tension parameters (cf_wrt and sigma)
    subroutine s_check_inputs_surface_tension
        if (f_is_default(sigma) .and. cf_wrt) then
            call s_mpi_abort('cf_wrt can only be anabled if the surface'// &
                             'coefficient is set')
        end if
    end subroutine s_check_inputs_surface_tension

    !> Checks constraints on the absence of flow variables
    subroutine s_check_inputs_no_flow_variables
        if (.not. any([ &
                      (/rho_wrt, E_wrt, pres_wrt, &
                        gamma_wrt, heat_ratio_wrt, &
                        pi_inf_wrt, pres_inf_wrt, &
                        cons_vars_wrt, prim_vars_wrt, &
                        c_wrt, schlieren_wrt/), &
                      alpha_rho_wrt, mom_wrt, vel_wrt, flux_wrt, &
                      alpha_wrt, omega_wrt])) then
            call s_mpi_abort('None of the flow variables have been '// &
                             'selected for post-process. Exiting ...')
        end if
    end subroutine s_check_inputs_no_flow_variables

end module m_checker
