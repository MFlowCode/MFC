!>
!! @file m_checker.f90
!! @brief Contains module m_checker

!> @brief The purpose of the module is to check for compatible input files
module m_checker

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper

    implicit none

    private; public :: s_check_inputs

contains

    subroutine s_check_inputs()

        character(len=5) :: iStr, numStr !< for int to string conversion
        integer :: i

        ! Constraints on the domain size
        if (nGlobal < 2**(min(1, m) + min(1, n) + min(1, p))*num_procs) then
            call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)
            call s_mpi_abort('Total number of cells must be at least '// &
                             '(2^[number of dimensions])*num_procs, which is currently '// &
                             trim(numStr)//'. Exiting ...')
        end if

        ! Constraints on the format of the formatted database file(s)
        if (format /= 1 .and. format /= 2) then
            call s_mpi_abort('format must be 1 or 2. Exiting ...')

            ! Constraints on the precision of the formatted database file(s)
        elseif (precision /= 1 .and. precision /= 2) then
            call s_mpi_abort('precision must be 1 or 2. Exiting ...')
        end if

        ! Constraints on the post-processing of the partial densities
        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (alpha_rho_wrt(i) .and. model_eqns == 1) then
                call s_mpi_abort('alpha_rho_wrt('//trim(iStr)//') is not '// &
                                 'supported for model_eqns = 1. Exiting ...')
            end if
            if (alpha_rho_wrt(i) .and. i > num_fluids) then
                call s_mpi_abort('Index of alpha_rho_wrt('//trim(iStr)//') '// &
                                 'exceeds the total number of fluids. Exiting ...')
            end if
        end do

        ! Constraints on the post-processing of the momentum
        if (n == 0 .and. mom_wrt(2)) then
            call s_mpi_abort('mom_wrt(2) is not supported for n = 0. Exiting ...')
        elseif (p == 0 .and. mom_wrt(3)) then
            call s_mpi_abort('mom_wrt(3) is not supported for p = 0. Exiting ...')

            ! Constraints on the post-processing of the velocity
        elseif (n == 0 .and. vel_wrt(2)) then
            call s_mpi_abort('vel_wrt(2) is not supported for n = 0. Exiting ...')
        elseif (p == 0 .and. vel_wrt(3)) then
            call s_mpi_abort('vel_wrt(3) is not supported for p = 0. Exiting ...')
        end if

        ! Constraints on the post-processing of the flux limiter function
        if (n == 0 .and. flux_wrt(2)) then
            call s_mpi_abort('flux_wrt(2) is not supported for n = 0. Exiting ...')
        elseif (p == 0 .and. flux_wrt(3)) then
            call s_mpi_abort('flux_wrt(3) is not supported for p = 0. Exiting ...')
        elseif (all(flux_lim /= (/dflt_int, 1, 2, 3, 4, 5, 6, 7/))) then
            call s_mpi_abort('flux_lim must be between 1 and 7. Exiting ...')
        end if

        ! Constraints on the post-processing of the volume fractions
        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (alpha_wrt(i) .and. model_eqns == 1) then
                call s_mpi_abort('alpha_wrt('//trim(iStr)//') is not '// &
                                 'supported for model_eqns = 1. Exiting ...')
            end if
            if (alpha_wrt(i) .and. i > num_fluids) then
                call s_mpi_abort('Index of alpha_wrt('//trim(iStr)//') '// &
                                 'exceeds the total number of fluids. Exiting ...')
            end if
        end do

        ! Constraints on the post-processing of the vorticity
        if (p == 0 .and. omega_wrt(1)) then
            call s_mpi_abort('omega_wrt(1) is not supported for p = 0. Exiting ...')
        elseif (p == 0 .and. omega_wrt(2)) then
            call s_mpi_abort('omega_wrt(2) is not supported for p = 0. Exiting ...')

            ! Constraints on post-processing of numerical Schlieren function
        elseif (n == 0 .and. schlieren_wrt) then
            call s_mpi_abort('schlieren_wrt is not supported for n = 0. Exiting ...')

            ! Constraints on post-processing combination of flow variables
        elseif ((any(alpha_rho_wrt) .neqv. .true.) &
                .and. &
                (any(mom_wrt) .neqv. .true.) &
                .and. &
                (any(vel_wrt) .neqv. .true.) &
                .and. &
                (any(flux_wrt) .neqv. .true.) &
                .and. &
                (any((/rho_wrt, E_wrt, pres_wrt, &
                       gamma_wrt, heat_ratio_wrt, &
                       pi_inf_wrt, pres_inf_wrt, &
                       cons_vars_wrt, &
                       prim_vars_wrt, &
                       c_wrt, schlieren_wrt/)) .neqv. .true.) &
                .and. &
                (any(alpha_wrt) .neqv. .true.) &
                .and. &
                (any(omega_wrt) .neqv. .true.)) then
            call s_mpi_abort('None of the flow variables have been '// &
                             'selected for post-process. Exiting ...')
        end if

        ! Constraints on the coefficients of numerical Schlieren function
        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (schlieren_alpha(i) /= dflt_real .and. schlieren_alpha(i) <= 0d0) then
                call s_mpi_abort('schlieren_alpha('//trim(iStr)//') must be '// &
                                 'greater than zero. Exiting ...')
            elseif (i > num_fluids .and. schlieren_alpha(i) /= dflt_real) then
                call s_mpi_abort('Index of schlieren_alpha('//trim(iStr)//') '// &
                                 'exceeds the total number of fluids. Exiting ...')
            elseif (.not. schlieren_wrt .and. schlieren_alpha(i) /= dflt_real) then
                call s_mpi_abort('schlieren_alpha('//trim(iStr)//') should '// &
                                 'be set only with schlieren_wrt enabled. Exiting ...')
            end if
        end do

        ! Constraints on the order of the finite difference scheme
        if ((any(omega_wrt) .or. schlieren_wrt) .and. fd_order == dflt_int) then
            call s_mpi_abort('fd_order must be set for omega_wrt or '// &
                             'schlieren_wrt. Exiting ...')
        end if

        ! Constraints on the surface tension model
        if (sigma == dflt_real .and. cf_wrt) then
            call s_mpi_abort('cf_wrt can only be anabled if the surface'// &
                             'coefficient is set')
        end if

    end subroutine s_check_inputs

end module m_checker
