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

        integer :: bub_fac
        integer :: i
        character(len=5) :: iStr

        bub_fac = 0; 
        if (bubbles .and. (num_fluids == 1)) bub_fac = 1

        ! Constraints on dimensionality and the number of cells for the grid
        if (m <= 0) then
            call s_mpi_abort('Unsupported choice for the value of m. '// &
                             'Exiting ...')
        elseif (n < 0) then
            call s_mpi_abort('Unsupported choice for the value of n. '// &
                             'Exiting ...')
        elseif (p < 0) then
            call s_mpi_abort('Unsupported choice for the value of p. '// &
                             'Exiting ...')
        elseif (cyl_coord .and. p > 0 .and. mod(p, 2) /= 1) then
            call s_mpi_abort('Unsupported choice for the value of p. '// &
                             'Exiting ...')
        elseif (n == 0 .and. p > 0) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and p. Exiting ...')
        elseif (nGlobal < 2**(min(1, m) + min(1, n) + min(1, p))*num_procs) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for num_procs, m, n and p. '// &
                             'Exiting ...')

            ! Constraints on the time-stepping parameters
        elseif (t_step_start < 0) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             't_step_start. Exiting ...')
        elseif (t_step_stop < t_step_start) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for t_step_start and t_step_stop. '// &
                             'Exiting ...')
        elseif (t_step_save > t_step_stop - t_step_start) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for t_step_start, t_step_stop and '// &
                             't_step_save. Exiting ...')

            ! Constraints on model equations and number of fluids in the flow
        elseif (all(model_eqns /= (/1, 2, 3, 4/))) then
            call s_mpi_abort('Unsupported value of model_eqns. Exiting ...')
        elseif (num_fluids /= dflt_int &
                .and. &
                (num_fluids < 1 .or. num_fluids > num_fluids)) then
            call s_mpi_abort('Unsupported value of num_fluids. Exiting ...')
        elseif ((model_eqns == 1 .and. num_fluids /= dflt_int) &
                .or. &
                (model_eqns == 2 .and. num_fluids == dflt_int) &
                .or. &
                (model_eqns == 3 .and. num_fluids == dflt_int)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'model_eqns and num_fluids. '// &
                             'Exiting ...')
        elseif (model_eqns == 1 .and. adv_alphan) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'model_eqns and adv_alphan. '// &
                             'Exiting ...')

            ! Constraints on the order of the WENO scheme
        elseif (weno_order /= 1 .and. weno_order /= 3 &
                .and. &
                weno_order /= 5) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'weno_order. Exiting ...')
        elseif (m + 1 < weno_order) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for m and weno_order. Exiting ...')
        elseif (n > 0 .and. n + 1 < weno_order) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and weno_order. Exiting ...')
        elseif (p > 0 .and. p + 1 < weno_order) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for p and weno_order. Exiting ...')
        elseif (nGlobal < weno_order**(min(1, m) + min(1, n) + min(1, p))*num_procs) &
            then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for num_procs, m, n, p and '// &
                             'weno_order. Exiting ...')

            ! Constraints on the boundary conditions in the x-direction
        elseif (bc_x%beg < -16 .or. bc_x%beg > -1 .or. bc_x%beg == -14) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'bc_x%beg. Exiting ...')
        elseif (bc_x%end < -16 .or. bc_x%end > -1 .or. bc_x%beg == -14) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'bc_x%end. Exiting ...')
        elseif ((bc_x%beg == -1 .and. bc_x%end /= -1) &
                .or. &
                (bc_x%end == -1 .and. bc_x%beg /= -1)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for bc_x%beg and bc_x%end. '// &
                             'Exiting ...')

            ! Constraints on the boundary conditions in the y-direction
        elseif (bc_y%beg /= dflt_int &
                .and. &
                ((((cyl_coord .neqv. .true.) &
                   .or. &
                   (cyl_coord .and. p == 0)) &
                  .and. &
                  (bc_y%beg < -16 .or. bc_y%beg > -1 .or. bc_y%beg == -14)) &
                 .or. &
                 (cyl_coord .and. p > 0 &
                  .and. &
                  (bc_y%beg < -16 .or. bc_y%beg > -1)))) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'bc_y%beg. Exiting ...')
        elseif (bc_y%end /= dflt_int &
                .and. &
                (bc_y%end < -16 .or. bc_y%end > -1 .or. bc_y%end == -14)) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'bc_y%end. Exiting ...')
        elseif ((n == 0 .and. bc_y%beg /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%beg == dflt_int)) then
            call s_mpi_abort('Unsupported choice for the value of n and '// &
                             'bc_y%beg. Exiting ...')
        elseif ((n == 0 .and. bc_y%end /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%end == dflt_int)) then
            call s_mpi_abort('Unsupported choice for the value of n and '// &
                             'bc_y%end. Exiting ...')
        elseif (n > 0 &
                .and. &
                ((bc_y%beg == -1 .and. bc_y%end /= -1) &
                 .or. &
                 (bc_y%end == -1 .and. bc_y%beg /= -1))) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n, bc_y%beg and bc_y%end. '// &
                             'Exiting ...')

            ! Constraints on the boundary conditions in the z-direction
        elseif (bc_z%beg /= dflt_int &
                .and. &
                (bc_z%beg < -16 .or. bc_z%beg > -1 .or. bc_z%beg == -14)) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'bc_z%beg. Exiting ...')
        elseif (bc_z%end /= dflt_int &
                .and. &
                (bc_z%end < -16 .or. bc_z%end > -1 .or. bc_z%end == -14)) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'bc_z%end. Exiting ...')
        elseif (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == -13)) then
            call s_mpi_abort('Unsupported choice of boundary condition -13')
        elseif ((p == 0 .and. bc_z%beg /= dflt_int) &
                .or. &
                (p > 0 .and. bc_z%beg == dflt_int)) then
            call s_mpi_abort('Unsupported choice for the value of p and '// &
                             'bc_z%beg. Exiting ...')
        elseif ((p == 0 .and. bc_z%end /= dflt_int) &
                .or. &
                (p > 0 .and. bc_z%end == dflt_int)) then
            call s_mpi_abort('Unsupported choice for the value of p and '// &
                             'bc_z%end. Exiting ...')
        elseif (p > 0 &
                .and. &
                ((bc_z%beg == -1 .and. bc_z%end /= -1) &
                 .or. &
                 (bc_z%end == -1 .and. bc_z%beg /= -1))) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for p, bc_z%beg and bc_z%end. '// &
                             'Exiting ...')
        end if

        ! Constraints on the stiffened equation of state fluids parameters
        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (fluid_pp(i)%gamma /= dflt_real &
                .and. &
                fluid_pp(i)%gamma <= 0d0) then
                call s_mpi_abort('Unsupported value of '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%gamma /= dflt_real) then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of model_eqns '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%gamma <= 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. fluid_pp(i)%gamma /= dflt_real)) &
                then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of num_fluids '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif (fluid_pp(i)%pi_inf /= dflt_real &
                    .and. &
                    fluid_pp(i)%pi_inf < 0d0) then
                call s_mpi_abort('Unsupported value of '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%pi_inf /= dflt_real) then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of model_eqns '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf. Exiting ...')
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%pi_inf < 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. fluid_pp(i)%pi_inf /= dflt_real)) &
                then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of num_fluids '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf. Exiting ...')
            end if

        end do

        ! Constraints on the format of the formatted database file(s)
        if (format /= 1 .and. format /= 2) then
            call s_mpi_abort('Unsupported choice for the value of format. '// &
                             'Exiting ...')

            ! Constraints on the precision of the formatted database file(s)
        elseif (precision /= 1 .and. precision /= 2) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'precision. Exiting ...')
        end if

        ! Constraints on the post-processing of the partial densities
        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (((i > num_fluids .or. model_eqns == 1) &
                 .and. &
                 alpha_rho_wrt(i)) &
                .or. &
                ((i <= num_fluids .and. model_eqns == 1) &
                 .and. &
                 alpha_rho_wrt(i))) then
                call s_mpi_abort('Unsupported choice of the '// &
                                 'combination of values for '// &
                                 'model_eqns, num_fluids and '// &
                                 'alpha_rho_wrt('//trim(iStr)//'). Exiting ...')
            end if
        end do

        ! Constraints on the post-processing of the momentum
        if (n == 0 .and. mom_wrt(2)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and mom_wrt(2). Exiting ...')
        elseif (n == 0 .and. mom_wrt(3)) then
            call s_mpi_abort('Unsupported cohice of the combination of '// &
                             'values for n and mom_wrt(3). Exiting ...')
        elseif (p == 0 .and. mom_wrt(3)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for p and mom_wrt(3). Exiting ...')

            ! Constraints on the post-processing of the velocity
        elseif (n == 0 .and. vel_wrt(2)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and vel_wrt(2). Exiting ...')
        elseif (n == 0 .and. vel_wrt(3)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and vel_wrt(3). Exiting ...')
        elseif (p == 0 .and. vel_wrt(3)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for p and vel_wrt(3). Exiting ...')
        end if

        ! Constraints on the post-processing of the flux limiter function
        if (n == 0 .and. flux_wrt(2)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and flux_wrt(2). Exiting ...')
        elseif (n == 0 .and. flux_wrt(3)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and flux_wrt(3). Exiting ...')
        elseif (p == 0 .and. flux_wrt(3)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for p and flux_wrt(3). Exiting ...')
        elseif (all(flux_lim /= (/dflt_int, 1, 2, 3, 4, 5, 6, 7/))) then
            call s_mpi_abort('Unsupported value of flux_lim. Exiting ...')
        end if

        ! Constraints on the post-processing of the volume fractions
        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (((i > num_fluids .or. model_eqns == 1) &
                 .and. &
                 alpha_wrt(i)) &
                .or. &
                ((i <= num_fluids .and. model_eqns == 1) &
                 .and. &
                 alpha_wrt(i))) then
                call s_mpi_abort('Unsupported choice of the '// &
                                 'combination of values for '// &
                                 'model_eqns, num_fluids and '// &
                                 'alpha_wrt('//trim(iStr)//'). Exiting ...')
            end if
        end do

        ! Constraints on the post-processing of the vorticity
        if (n == 0 .and. omega_wrt(1)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and omega_wrt(1). Exiting ...')
        elseif (n == 0 .and. omega_wrt(2)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and omega_wrt(2). Exiting ...')
        elseif (n == 0 .and. omega_wrt(3)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and omega_wrt(3). Exiting ...')
        elseif (p == 0 .and. omega_wrt(1)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for p and omega_wrt(1). Exiting ...')
        elseif (p == 0 .and. omega_wrt(2)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for p and omega_wrt(2). Exiting ...')

            ! Constraints on post-processing of numerical Schlieren function
        elseif (n == 0 .and. schlieren_wrt) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for n and schlieren_wrt. Exiting ...')

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
            if (schlieren_alpha(i) /= dflt_real &
                .and. &
                schlieren_alpha(i) <= 0d0) then
                call s_mpi_abort('Unsupported choice for the value of '// &
                                 'schlieren_alpha('//trim(iStr)//'). Exiting ...')
            elseif (((i > num_fluids .or. (schlieren_wrt .neqv. .true.)) &
                     .and. &
                     schlieren_alpha(i) /= dflt_real) &
                    .or. &
                    ((i <= num_fluids .and. schlieren_wrt) &
                     .and. &
                     schlieren_alpha(i) <= 0d0)) then
                call s_mpi_abort('Unsupported choice of the '// &
                                 'combination of values for '// &
                                 'num_fluids, schlieren_wrt and '// &
                                 'schlieren_alpha('//trim(iStr)//'). Exiting ...')
            end if
        end do

        ! Constraints on the order of the finite difference scheme
        if (fd_order /= dflt_int &
            .and. &
            fd_order /= 1 .and. fd_order /= 2 .and. fd_order /= 4) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'fd_order. Exiting ...')
        elseif ((any(omega_wrt) .or. schlieren_wrt) &
                .and. &
                fd_order == dflt_int) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for omega_wrt, schlieren_wrt and '// &
                             'fd_order. Exiting ...')
        end if

        ! Moving Boundaries Checks: x boundaries
        if (any((/bc_x%vb1, bc_x%vb2, bc_x%vb3/) /= 0d0)) then
            if (bc_x%beg == 15) then
                if (any((/bc_x%vb2, bc_x%vb3/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_x%beg and"// &
                                     "bc_x%vb2 or bc_x%vb3. Exiting ...")
                end if
            elseif (bc_x%beg /= -16) then
                call s_mpi_abort("Unsupported combination of bc_x%beg and"// &
                                 "bc_x%vb1, bc_x%vb2, or bc_x%vb3. Exiting...")
            end if
        end if

        if (any((/bc_x%ve1, bc_x%ve2, bc_x%ve3/) /= 0d0)) then
            if (bc_x%end == 15) then
                if (any((/bc_x%ve2, bc_x%ve3/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_x%end and"// &
                                     "bc_x%ve2 or bc_x%ve3. Exiting ...")
                end if
            elseif (bc_x%end /= -16) then
                call s_mpi_abort("Unsupported combination of bc_x%end and"// &
                                 "bc_x%ve1, bc_x%ve2, or bc_x%ve3. Exiting...")
            end if
        end if

        ! Moving Boundaries Checks: y boundaries
        if (any((/bc_y%vb1, bc_y%vb2, bc_y%vb3/) /= 0d0)) then
            if (bc_y%beg == 15) then
                if (any((/bc_y%vb1, bc_y%vb3/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_y%beg and"// &
                                     "bc_y%vb1 or bc_y%vb3. Exiting ...")
                end if
            elseif (bc_y%beg /= -16) then
                call s_mpi_abort("Unsupported combination of bc_y%beg and"// &
                                 "bc_y%vb1, bc_y%vb2, or bc_y%vb3. Exiting...")
            end if
        end if

        if (any((/bc_y%ve1, bc_y%ve2, bc_y%ve3/) /= 0d0)) then
            if (bc_y%end == 15) then
                if (any((/bc_y%ve1, bc_y%ve3/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_y%end and"// &
                                     "bc_y%ve1 or bc_y%ve3. Exiting ...")
                end if
            elseif (bc_y%end /= -16) then
                call s_mpi_abort("Unsupported combination of bc_y%end and"// &
                                 "bc_y%ve1, bc_y%ve2, or bc_y%ve3. Exiting...")
            end if
        end if

        ! Moving Boundaries Checks: z boundaries
        if (any((/bc_z%vb1, bc_z%vb2, bc_z%vb3/) /= 0d0)) then
            if (bc_z%beg == 15) then
                if (any((/bc_x%vb1, bc_x%vb2/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_z%beg and"// &
                                     "bc_x%vb1 or bc_x%vb1. Exiting ...")
                end if
            elseif (bc_z%beg /= -16) then
                call s_mpi_abort("Unsupported combination of bc_z%beg and"// &
                                 "bc_z%vb1, bc_z%vb2, or bc_z%vb3. Exiting...")
            end if
        end if

        if (any((/bc_z%ve1, bc_z%ve2, bc_z%ve3/) /= 0d0)) then
            if (bc_z%end == 15) then
                if (any((/bc_x%ve1, bc_x%ve2/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_z%end and"// &
                                     "bc_z%ve2 or bc_z%ve3. Exiting ...")
                end if
            elseif (bc_z%end /= -16) then
                call s_mpi_abort("Unsupported combination of bc_z%end and"// &
                                 "bc_z%ve1, bc_z%ve2, or bc_z%ve3. Exiting...")
            end if
        end if

        ! Constraints on the surface tension model
        if (sigma /= dflt_real .and. sigma < 0d0) then
            call s_mpi_abort('The surface tension coefficient must be'// &
                             'greater than or equal to zero. Exiting ...')
        elseif (sigma /= dflt_real .and. model_eqns /= 3) then
            call s_mpi_abort("The surface tension model requires"// &
                             'model_eqns=3. Exiting ...')
        elseif (sigma == dflt_real .and. cf_wrt) then
            call s_mpi_abort('cf_wrt can only be anabled if the surface'// &
                             'coefficient is set')
        end if

    end subroutine s_check_inputs

end module m_checker
