!>
!!@file m_checker.f90
!!@brief Contains module m_checker

#:include 'macros.fpp'

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
        @:PROHIBIT(format /= 1 .and. format /= 2)
        @:PROHIBIT(precision /= 1 .and. precision /= 2)
    end subroutine s_check_inputs_output_format

    !> Checks constraints on partial density parameters
    subroutine s_check_inputs_partial_density
        character(len=5) :: iStr
        integer :: i

        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            @:PROHIBIT(alpha_rho_wrt(i) .and. model_eqns == 1, "alpha_rho_wrt("//trim(iStr)//") is not supported for model_eqns = 1")
            @:PROHIBIT(alpha_rho_wrt(i) .and. i > num_fluids, "Index of alpha_rho_wrt("//trim(iStr)//") exceeds the total number of fluids")
        end do
    end subroutine s_check_inputs_partial_density

    !> Checks constraints on momentum parameters
    subroutine s_check_inputs_momentum
        @:PROHIBIT(n == 0 .and. mom_wrt(2))
        @:PROHIBIT(p == 0 .and. mom_wrt(3))
    end subroutine s_check_inputs_momentum

    !> Checks constraints on velocity parameters
    subroutine s_check_inputs_velocity
        @:PROHIBIT(n == 0 .and. vel_wrt(2))
        @:PROHIBIT(p == 0 .and. vel_wrt(3))
    end subroutine s_check_inputs_velocity

    !> Checks constraints on flux limiter parameters
    subroutine s_check_inputs_flux_limiter
        @:PROHIBIT(n == 0 .and. flux_wrt(2))
        @:PROHIBIT(p == 0 .and. flux_wrt(3))
        @:PROHIBIT(all(flux_lim /= (/dflt_int, 1, 2, 3, 4, 5, 6, 7/)), "flux_lim must be between 1 and 7")
    end subroutine s_check_inputs_flux_limiter

    !> Checks constraints on volume fraction parameters
    subroutine s_check_inputs_volume_fraction
        character(len=5) :: iStr
        integer :: i

        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            @:PROHIBIT(alpha_wrt(i) .and. model_eqns == 1, "alpha_wrt("//trim(iStr)//") is not supported for model_eqns = 1")
            @:PROHIBIT(alpha_wrt(i) .and. i > num_fluids, "Index of alpha_wrt("//trim(iStr)//") exceeds the total number of fluids")
        end do
    end subroutine s_check_inputs_volume_fraction

    !> Checks constraints on vorticity parameters
    subroutine s_check_inputs_vorticity
        @:PROHIBIT(n == 0 .and. any(omega_wrt))
        @:PROHIBIT(p == 0 .and. (omega_wrt(1) .or. omega_wrt(2)))
        @:PROHIBIT(any(omega_wrt) .and. fd_order == dflt_int, "fd_order must be set for omega_wrt")
    end subroutine s_check_inputs_vorticity

    !> Checks constraints on numerical Schlieren parameters
        !! (schlieren_wrt and schlieren_alpha)
    subroutine s_check_inputs_schlieren
        character(len=5) :: iStr
        integer :: i

        @:PROHIBIT(n == 0 .and. schlieren_wrt)
        @:PROHIBIT(schlieren_wrt .and. fd_order == dflt_int, "fd_order must be set for schlieren_wrt")

        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            @:PROHIBIT(.not. f_is_default(schlieren_alpha(i)) .and. schlieren_alpha(i) <= 0d0, &
                "schlieren_alpha("//trim(iStr)//") must be greater than zero")
            @:PROHIBIT(.not. f_is_default(schlieren_alpha(i)) .and. i > num_fluids, &
                "Index of schlieren_alpha("//trim(iStr)//") exceeds the total number of fluids")
            @:PROHIBIT(.not. f_is_default(schlieren_alpha(i)) .and. (.not. schlieren_wrt), &
                "schlieren_alpha("//trim(iStr)//") should be set only with schlieren_wrt enabled")
        end do
    end subroutine s_check_inputs_schlieren

    !> Checks constraints on surface tension parameters (cf_wrt and sigma)
    subroutine s_check_inputs_surface_tension
        @:PROHIBIT(cf_wrt .and. f_is_default(sigma), &
            "cf_wrt can only be enabled if the surface coefficient is set")
    end subroutine s_check_inputs_surface_tension

    !> Checks constraints on the absence of flow variables
    subroutine s_check_inputs_no_flow_variables
        @:PROHIBIT(.not. any([ &
            (/rho_wrt, E_wrt, pres_wrt, gamma_wrt, heat_ratio_wrt, pi_inf_wrt, &
              pres_inf_wrt, cons_vars_wrt, prim_vars_wrt, c_wrt, schlieren_wrt/), &
            alpha_rho_wrt, mom_wrt, vel_wrt, flux_wrt, alpha_wrt, omega_wrt]), &
            "None of the flow variables have been selected for post-process. Exiting ...")
    end subroutine s_check_inputs_no_flow_variables

end module m_checker
