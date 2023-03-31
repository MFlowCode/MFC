!>
!! @file m_start_up.f90
!! @brief  Contains module m_start_up

!> @brief This module contains the subroutines that read in and check the
!!              consistency of the user provided inputs.
module m_start_up

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_compile_specific

    use m_checker

    use m_helper
    ! ==========================================================================

    implicit none

contains

    !>  Reads the configuration file post_process.inp, in order
        !!      to populate parameters in module m_global_parameters.f90
        !!      with the user provided inputs
    subroutine s_read_input_file() ! ---------------------------------------

        character(LEN=name_len) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location

        integer :: iostatus
            !! Integer to check iostat of file read

        ! Namelist for all of the parameters to be inputed by the user
        namelist /user_inputs/ case_dir, m, n, p, t_step_start, &
            t_step_stop, t_step_save, model_eqns, &
            num_fluids, mpp_lim, adv_alphan, &
            weno_order, bc_x, &
            bc_y, bc_z, fluid_pp, format, precision, &
            hypoelasticity, G, &
            alpha_rho_wrt, rho_wrt, mom_wrt, vel_wrt, &
            E_wrt, pres_wrt, alpha_wrt, gamma_wrt, &
            heat_ratio_wrt, pi_inf_wrt, pres_inf_wrt, &
            cons_vars_wrt, prim_vars_wrt, c_wrt, &
            omega_wrt, qm_wrt, schlieren_wrt, schlieren_alpha, &
            fd_order, mixture_err, alt_soundspeed, &
            flux_lim, flux_wrt, cyl_coord, &
            parallel_io, coarsen_silo, &
            rhoref, pref, bubbles, R0ref, nb, &
            polytropic, thermal, Ca, Web, Re_inv, &
            polydisperse, poly_sigma

        ! Inquiring the status of the post_process.inp file
        file_loc = 'post_process.inp'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! Checking whether the input file is there. If it is, the input file
        ! is read. If not, the program is terminated.
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='formatted', &
                  STATUS='old', ACTION='read')
            read (1, NML=user_inputs, iostat=iostatus)

            if (iostatus /= 0) then
                call s_mpi_abort('Invalid line in post_process.inp. It is '// &
                    'likely due to a datatype mismatch. Exiting ...')
            end if

            close (1)
            ! Store m,n,p into global m,n,p
            m_glb = m
            n_glb = n
            p_glb = p
        else
            call s_mpi_abort('File post_process.inp is missing. Exiting ...')
        end if

    end subroutine s_read_input_file ! -------------------------------------

    !>  Checking that the user inputs make sense, i.e. that the
        !!      individual choices are compatible with the code's options
        !!      and that the combination of these choices results into a
        !!      valid configuration for the post-process
    subroutine s_check_input_file() ! --------------------------------------

        character(LEN=len_trim(case_dir)) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: dir_check !<
            !! Logical variable used to test the existence of folders

        ! Checking the existence of the case folder
        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir)//'/.'

        call my_inquire(file_loc, dir_check)

        ! Constraint on the location of the case directory
        if (dir_check .neqv. .true.) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                'case_dir. Exiting ...')
        end if

        call s_check_inputs()

    end subroutine s_check_input_file ! ------------------------------------

end module m_start_up
