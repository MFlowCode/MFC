!>
!! @file m_start_up.f90
!! @brief  Contains module m_start_up
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module contains the subroutines that read in and check the
!!              consistency of the user provided inputs.
module m_start_up

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_compile_specific
    ! ==========================================================================

    implicit none

contains

    !>  Checking that the user inputs make sense, i.e. that the
        !!      individual choices are compatible with the code's options
        !!      and that the combination of these choices results into a
        !!      valid configuration for the post-process
    subroutine s_check_input_file() ! --------------------------------------

        character(LEN=len_trim(case_dir)) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: dir_check !<
            !! Logical variable used to test the existence of folders

        integer :: i  !< Generic loop iterator
        integer :: bub_fac

        bub_fac = 0; 
        if (bubbles .and. (num_fluids == 1)) bub_fac = 1

        ! Checking the existence of the case folder
        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir)//'/.'

        call my_inquire(file_loc, dir_check)

        ! Constraint on the location of the case directory
        if (dir_check .neqv. .true.) then
            print '(A)', 'Unsupported choice for the value of '// &
                'case_dir. Exiting ...'
            call s_mpi_abort()

            ! Constraints on dimensionality and the number of cells for the grid
        elseif (m <= 0) then
            print '(A)', 'Unsupported choice for the value of m. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (n < 0) then
            print '(A)', 'Unsupported choice for the value of n. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (p < 0) then
            print '(A)', 'Unsupported choice for the value of p. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (cyl_coord .and. p > 0 .and. mod(p, 2) /= 1) then
            print '(A)', 'Unsupported choice for the value of p. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. p > 0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and p. Exiting ...'
            call s_mpi_abort()
        elseif ((m + 1)*(n + 1)*(p + 1) &
                < &
                2**(min(1, m) + min(1, n) + min(1, p))*num_procs) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for num_procs, m, n and p. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on the time-stepping parameters
        elseif (t_step_start < 0) then
            print '(A)', 'Unsupported choice for the value of '// &
                't_step_start. Exiting ...'
            call s_mpi_abort()
        elseif (t_step_stop < t_step_start) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for t_step_start and t_step_stop. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (t_step_save > t_step_stop - t_step_start) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for t_step_start, t_step_stop and '// &
                't_step_save. Exiting ...'
            call s_mpi_abort()

            ! Constraints on model equations and number of fluids in the flow
        elseif (all(model_eqns /= (/1, 2, 3, 4/))) then
            print '(A)', 'Unsupported value of model_eqns. Exiting ...'
            call s_mpi_abort()
        elseif (num_fluids /= dflt_int &
                .and. &
                (num_fluids < 1 .or. num_fluids > num_fluids_alloc)) then
            print '(A)', 'Unsupported value of num_fluids. Exiting ...'
            call s_mpi_abort()
        elseif ((model_eqns == 1 .and. num_fluids /= dflt_int) &
                .or. &
                (model_eqns == 2 .and. num_fluids == dflt_int) &
                .or. &
                (model_eqns == 3 .and. num_fluids == dflt_int)) then
            print '(A)', 'Unsupported combination of values of '// &
                'model_eqns and num_fluids. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (model_eqns == 1 .and. adv_alphan) then
            print '(A)', 'Unsupported combination of values of '// &
                'model_eqns and adv_alphan. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on the order of the WENO scheme
        elseif (weno_order /= 1 .and. weno_order /= 3 &
                .and. &
                weno_order /= 5) then
            print '(A)', 'Unsupported choice for the value of '// &
                'weno_order. Exiting ...'
            call s_mpi_abort()
        elseif (m + 1 < weno_order) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for m and weno_order. Exiting ...'
            call s_mpi_abort()
        elseif (n > 0 .and. n + 1 < weno_order) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and weno_order. Exiting ...'
            call s_mpi_abort()
        elseif (p > 0 .and. p + 1 < weno_order) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for p and weno_order. Exiting ...'
            call s_mpi_abort()
        elseif ((m + 1)*(n + 1)*(p + 1) &
                < &
                weno_order**(min(1, m) + min(1, n) + min(1, p))*num_procs) &
            then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for num_procs, m, n, p and '// &
                'weno_order. Exiting ...'
            call s_mpi_abort()

            ! Constraints on the boundary conditions in the x-direction
        elseif (bc_x%beg < -12 .or. bc_x%beg > -1) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_x%beg. Exiting ...'
            call s_mpi_abort()
        elseif (bc_x%end < -12 .or. bc_x%end > -1) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_x%end. Exiting ...'
            call s_mpi_abort()
        elseif ((bc_x%beg == -1 .and. bc_x%end /= -1) &
                .or. &
                (bc_x%end == -1 .and. bc_x%beg /= -1)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for bc_x%beg and bc_x%end. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on the boundary conditions in the y-direction
        elseif (bc_y%beg /= dflt_int &
                .and. &
                ((((cyl_coord .neqv. .true.) &
                   .or. &
                   (cyl_coord .and. p == 0)) &
                  .and. &
                  (bc_y%beg < -12 .or. bc_y%beg > -1)) &
                 .or. &
                 (cyl_coord .and. p > 0 &
                  .and. &
                  (bc_y%beg < -13 .or. bc_y%beg > -1)))) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_y%beg. Exiting ...'
            call s_mpi_abort()
        elseif (bc_y%end /= dflt_int &
                .and. &
                (bc_y%end < -12 .or. bc_y%end > -1)) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_y%end. Exiting ...'
            call s_mpi_abort()
        elseif ((n == 0 .and. bc_y%beg /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%beg == dflt_int)) then
            print '(A)', 'Unsupported choice for the value of n and '// &
                'bc_y%beg. Exiting ...'
            call s_mpi_abort()
        elseif ((n == 0 .and. bc_y%end /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%end == dflt_int)) then
            print '(A)', 'Unsupported choice for the value of n and '// &
                'bc_y%end. Exiting ...'
            call s_mpi_abort()
        elseif (n > 0 &
                .and. &
                ((bc_y%beg == -1 .and. bc_y%end /= -1) &
                 .or. &
                 (bc_y%end == -1 .and. bc_y%beg /= -1))) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n, bc_y%beg and bc_y%end. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on the boundary conditions in the z-direction
        elseif (bc_z%beg /= dflt_int &
                .and. &
                (bc_z%beg < -12 .or. bc_z%beg > -1)) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_z%beg. Exiting ...'
            call s_mpi_abort()
        elseif (bc_z%end /= dflt_int &
                .and. &
                (bc_z%end < -12 .or. bc_z%end > -1)) then
            print '(A)', 'Unsupported choice for the value of '// &
                'bc_z%end. Exiting ...'
            call s_mpi_abort()
        elseif ((p == 0 .and. bc_z%beg /= dflt_int) &
                .or. &
                (p > 0 .and. bc_z%beg == dflt_int)) then
            print '(A)', 'Unsupported choice for the value of p and '// &
                'bc_z%beg. Exiting ...'
            call s_mpi_abort()
        elseif ((p == 0 .and. bc_z%end /= dflt_int) &
                .or. &
                (p > 0 .and. bc_z%end == dflt_int)) then
            print '(A)', 'Unsupported choice for the value of p and '// &
                'bc_z%end. Exiting ...'
            call s_mpi_abort()
        elseif (p > 0 &
                .and. &
                ((bc_z%beg == -1 .and. bc_z%end /= -1) &
                 .or. &
                 (bc_z%end == -1 .and. bc_z%beg /= -1))) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for p, bc_z%beg and bc_z%end. '// &
                'Exiting ...'
            call s_mpi_abort()
        end if

        ! Constraints on the stiffened equation of state fluids parameters
        do i = 1, num_fluids_alloc

            if (fluid_pp(i)%gamma /= dflt_real &
                .and. &
                fluid_pp(i)%gamma <= 0d0) then
                print '(A,I0,A)', 'Unsupported value of '// &
                    'fluid_pp(', i, ')%'// &
                    'gamma. Exiting ...'
                call s_mpi_abort()
            elseif (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%gamma /= dflt_real) then
                print '(A,I0,A)', 'Unsupported combination '// &
                    'of values of model_eqns '// &
                    'and fluid_pp(', i, ')%'// &
                    'gamma. Exiting ...'
                call s_mpi_abort()
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%gamma <= 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. fluid_pp(i)%gamma /= dflt_real)) &
                then
                print '(A,I0,A)', 'Unsupported combination '// &
                    'of values of num_fluids '// &
                    'and fluid_pp(', i, ')%'// &
                    'gamma. Exiting ...'
                call s_mpi_abort()
            elseif (fluid_pp(i)%pi_inf /= dflt_real &
                    .and. &
                    fluid_pp(i)%pi_inf < 0d0) then
                print '(A,I0,A)', 'Unsupported value of '// &
                    'fluid_pp(', i, ')%'// &
                    'pi_inf. Exiting ...'
                call s_mpi_abort()
            elseif (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%pi_inf /= dflt_real) then
                print '(A,I0,A)', 'Unsupported combination '// &
                    'of values of model_eqns '// &
                    'and fluid_pp(', i, ')%'// &
                    'pi_inf. Exiting ...'
                call s_mpi_abort()
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%pi_inf < 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. fluid_pp(i)%pi_inf /= dflt_real)) &
                then
                print '(A,I0,A)', 'Unsupported combination '// &
                    'of values of num_fluids '// &
                    'and fluid_pp(', i, ')%'// &
                    'pi_inf. Exiting ...'
                call s_mpi_abort()
            end if

        end do

        ! Constraints on the format of the formatted database file(s)
        if (format /= 1 .and. format /= 2) then
            print '(A)', 'Unsupported choice for the value of format. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif ((precision /= 2) .and. (parallel_io .neqv. .false.)) then
            print '(A)', 'Unsupported combination of precision and parallel IO. '// &
                'Please use precision == 2 when enabling parallel_io.  Exiting ...'
            call s_mpi_abort()

            ! Constraints on the precision of the formatted database file(s)
        elseif (precision /= 1 .and. precision /= 2) then
            print '(A)', 'Unsupported choice for the value of '// &
                'precision. Exiting ...'
            call s_mpi_abort()

            ! Constraints on the option to coarsen the formatted database files
        elseif (coarsen_silo .and. format /= 1) then
            print '(A)', 'Unsupported combination of values of format '// &
                'and coarsen_silo. Exiting ...'
            call s_mpi_abort()
        elseif (coarsen_silo .and. n == 0) then
            print '(A)', 'Unsupported combination of values of n '// &
                'and coarsen_silo. Exiting ...'
            call s_mpi_abort()
        end if

        ! Constraints on the post-processing of the partial densities
        do i = 1, num_fluids_alloc
            if (((i > num_fluids .or. model_eqns == 1) &
                 .and. &
                 alpha_rho_wrt(i)) &
                .or. &
                ((i <= num_fluids .and. model_eqns == 1) &
                 .and. &
                 alpha_rho_wrt(i))) then
                print '(A,I0,A)', 'Unsupported choice of the '// &
                    'combination of values for '// &
                    'model_eqns, num_fluids and '// &
                    'alpha_rho_wrt(', i, '). Exiting ...'
                call s_mpi_abort()
            end if
        end do

        ! Constraints on the post-processing of the momentum
        if (n == 0 .and. mom_wrt(2)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and mom_wrt(2). Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. mom_wrt(3)) then
            print '(A)', 'Unsupported cohice of the combination of '// &
                'values for n and mom_wrt(3). Exiting ...'
            call s_mpi_abort()
        elseif (p == 0 .and. mom_wrt(3)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for p and mom_wrt(3). Exiting ...'
            call s_mpi_abort()

            ! Constraints on the post-processing of the velocity
        elseif (n == 0 .and. vel_wrt(2)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and vel_wrt(2). Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. vel_wrt(3)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and vel_wrt(3). Exiting ...'
            call s_mpi_abort()
        elseif (p == 0 .and. vel_wrt(3)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for p and vel_wrt(3). Exiting ...'
            call s_mpi_abort()
        end if

        ! Constraints on the post-processing of the flux limiter function
        if (n == 0 .and. flux_wrt(2)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and flux_wrt(2). Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. flux_wrt(3)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and flux_wrt(3). Exiting ...'
            call s_mpi_abort()
        elseif (p == 0 .and. flux_wrt(3)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for p and flux_wrt(3). Exiting ...'
            call s_mpi_abort()
        elseif (all(flux_lim /= (/dflt_int, 1, 2, 3, 4, 5, 6, 7/))) then
            print '(A)', 'Unsupported value of flux_lim. Exiting ...'
            call s_mpi_abort()
        end if

        ! Constraints on the post-processing of the volume fractions
        do i = 1, num_fluids_alloc
            if (((i > num_fluids .or. model_eqns == 1) &
                 .and. &
                 alpha_wrt(i)) &
                .or. &
                ((i <= num_fluids .and. model_eqns == 1) &
                 .and. &
                 alpha_wrt(i))) then
                print '(A,I0,A)', 'Unsupported choice of the '// &
                    'combination of values for '// &
                    'model_eqns, num_fluids and '// &
                    'alpha_wrt(', i, '). Exiting ...'
                call s_mpi_abort()
            end if
        end do

        ! Constraints on the post-processing of the vorticity
        if (n == 0 .and. omega_wrt(1)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and omega_wrt(1). Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. omega_wrt(2)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and omega_wrt(2). Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. omega_wrt(3)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and omega_wrt(3). Exiting ...'
            call s_mpi_abort()
        elseif (p == 0 .and. omega_wrt(1)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for p and omega_wrt(1). Exiting ...'
            call s_mpi_abort()
        elseif (p == 0 .and. omega_wrt(2)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for p and omega_wrt(2). Exiting ...'
            call s_mpi_abort()

            ! Constraints on post-processing of numerical Schlieren function
        elseif (n == 0 .and. schlieren_wrt) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and schlieren_wrt. Exiting ...'
            call s_mpi_abort()

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
            print '(A)', 'None of the flow variables have been '// &
                'selected for post-process. Exiting ...'
            call s_mpi_abort()
        end if

        ! Constraints on the coefficients of numerical Schlieren function
        do i = 1, num_fluids_alloc
            if (schlieren_alpha(i) /= dflt_real &
                .and. &
                schlieren_alpha(i) <= 0d0) then
                print '(A,I0,A)', 'Unsupported choice for the value of '// &
                    'schlieren_alpha(', i, '). Exiting ...'
                call s_mpi_abort()
            elseif (((i > num_fluids .or. (schlieren_wrt .neqv. .true.)) &
                     .and. &
                     schlieren_alpha(i) /= dflt_real) &
                    .or. &
                    ((i <= num_fluids .and. schlieren_wrt) &
                     .and. &
                     schlieren_alpha(i) <= 0d0)) then
                print '(A,I0,A)', 'Unsupported choice of the '// &
                    'combination of values for '// &
                    'num_fluids, schlieren_wrt and '// &
                    'schlieren_alpha(', i, '). Exiting ...'
                call s_mpi_abort()
            end if
        end do

        ! Constraints on the order of the finite difference scheme
        if (fd_order /= dflt_int &
            .and. &
            fd_order /= 1 .and. fd_order /= 2 .and. fd_order /= 4) then
            print '(A)', 'Unsupported choice for the value of '// &
                'fd_order. Exiting ...'
            call s_mpi_abort()
            !          ELSEIF(               (omega_wrt(1) .NEQV. .TRUE.)            &
            !                                           .AND.                        &
            !                                (omega_wrt(2) .NEQV. .TRUE.)            &
            !                                           .AND.                        &
            !                                (omega_wrt(3) .NEQV. .TRUE.)            &
            !                                           .AND.                        &
            !                                !(schlieren_wrt .NEQV. .TRUE.)           &
            !                                !           .AND.                        &
            !                                   fd_order /= dflt_int              ) THEN
            !              PRINT '(A)', 'AA Unsupported choice of the combination of '    // &
            !                           'values for omega_wrt, schlieren_wrt and '     // &
            !                           'fd_order. Exiting ...'
            !              CALL s_mpi_abort()
        elseif ((any(omega_wrt) .or. schlieren_wrt) &
                .and. &
                fd_order == dflt_int) then
            print '(A)', 'BB Unsupported choice of the combination of '// &
                'values for omega_wrt, schlieren_wrt and '// &
                'fd_order. Exiting ...'
            call s_mpi_abort()
        end if

    end subroutine s_check_input_file ! ------------------------------------

end module m_start_up
