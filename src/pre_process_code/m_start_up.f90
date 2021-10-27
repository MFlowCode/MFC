!>
!! @file m_start_up.f90
!! @brief Contains module m_start_up
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module contains subroutines that read, and check consistency
!!              of, the user provided inputs and data.
module m_start_up

    ! Dependencies =============================================================
    use m_derived_types          !< Definitions of the derived types

    use m_global_parameters      !< Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_data_output            !< Procedures to write the grid data and the
                                 !! conservative variables to files

    use mpi                      !< Message passing interface (MPI) module

    use m_compile_specific
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_start_up_module, &
 s_read_input_file, &
 s_check_input_file, &
 s_read_grid_data_files, &
 s_read_ic_data_files, &
 s_read_serial_grid_data_files, &
 s_read_serial_ic_data_files, &
 s_read_parallel_grid_data_files, &
 s_read_parallel_ic_data_files, &
 s_check_grid_data_files, &
 s_finalize_start_up_module

    abstract interface ! ===================================================

        subroutine s_read_abstract_grid_data_files(dflt_int)! ----------

            integer, intent(IN) :: dflt_int

        end subroutine s_read_abstract_grid_data_files ! ---------------

        subroutine s_read_abstract_ic_data_files(q_cons_vf) ! -----------

            import :: scalar_field, sys_size

            ! Conservative variables
            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: q_cons_vf

        end subroutine s_read_abstract_ic_data_files ! -----------------

    end interface ! ========================================================

    character(LEN=path_len + name_len) :: proc_rank_dir !<
    !! Location of the folder associated with the rank of the local processor

    character(LEN=path_len + 2*name_len), private :: t_step_dir !<
    !! Possible location of time-step folder containing preexisting grid and/or
    !! conservative variables data to be used as starting point for pre-process

    procedure(s_read_abstract_grid_data_files), pointer :: s_read_grid_data_files => null()
    procedure(s_read_abstract_ic_data_files), pointer :: s_read_ic_data_files => null()

contains

    !>  Reads the configuration file pre_process.inp, in order to
        !!      populate the parameters in module m_global_parameters.f90
        !!      with the user provided inputs
    subroutine s_read_input_file() ! ---------------------------------------

        character(LEN=name_len) :: file_loc  !<
            !! Generic string used to store the address of a particular file

        logical :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location

        ! Namelist for all of the parameters to be inputed by the user
        namelist /user_inputs/ case_dir, old_grid, old_ic, t_step_old, m, &
            n, p, x_domain, y_domain, z_domain, &
            stretch_x, stretch_y, stretch_z, a_x, a_y, &
            a_z, x_a, y_a, z_a, x_b, y_b, z_b, &
            model_eqns, num_fluids, &
            adv_alphan, mpp_lim, &
            weno_order, bc_x, bc_y, bc_z, num_patches, &
            patch_icpp, fluid_pp, &
            precision, parallel_io, &
            perturb_flow, perturb_flow_fluid, &
            perturb_sph, perturb_sph_fluid, fluid_rho, &
            cyl_coord, loops_x, loops_y, loops_z, &
            rhoref, pref, bubbles, R0ref, nb, &
            polytropic, thermal, Ca, Web, Re_inv, &
            polydisperse, poly_sigma, qbmm, &
            nnode, sigR, sigV, dist_type, rhoRV, R0_type

        ! Inquiring the status of the pre_process.inp file
        file_loc = 'pre_process.inp'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! Checking whether the input file is there. If it is, the input file
        ! is read. If not, the program is terminated.
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='formatted', &
                  STATUS='old', ACTION='read')
            read (1, NML=user_inputs)
            close (1)
            ! Store m,n,p into global m,n,p
            m_glb = m
            n_glb = n
            p_glb = p
        else
            print '(A)', 'File pre_process.inp is missing. Exiting ...'
            call s_mpi_abort()
        end if

    end subroutine s_read_input_file ! -------------------------------------

    !>  Checking that the user inputs make sense, i.e. that the
        !!      individual choices are compatible with the code's options
        !!      and that the combination of these choices results into a
        !!      valid configuration for the pre-process
    subroutine s_check_input_file() ! --------------------------------------

        character(LEN=len_trim(case_dir)) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: dir_check !<
            !! Logical variable used to test the existence of folders

        integer :: i !<
            !! Generic loop iterator

        integer :: bub_fac !<
            !! For allowing an extra fluid_pp if there are subgrid bubbles

        bub_fac = 0
        if (bubbles .and. (num_fluids == 1)) bub_fac = 1

        ! Checking the existence of the case folder
        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir)//'/.'

        call my_inquire(file_loc, dir_check)

        ! Startup checks for bubbles and bubble variables
        if (bubbles .and. (model_eqns .ne. 4 .and. model_eqns .ne. 2)) then
            print '(A)', 'Unsupported combination of values of '// &
                'bubbles and model_eqns. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (bubbles .and. polydisperse .and. (nb == 1)) then
            print '(A)', 'Polydisperse bubble dynamics requires nb > 1 '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (bubbles .and. polydisperse .and. (mod(nb, 2) == 0)) then
            print '(A)', 'nb must be odd '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (model_eqns == 4 .and. (rhoref == dflt_real)) then
            print '(A)', 'Unsupported combination of values of '// &
                'bubbles and rhoref. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (model_eqns == 4 .and. (pref == dflt_real)) then
            print '(A)', 'Unsupported combination of values of '// &
                'bubbles and pref. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (model_eqns == 4 .and. (num_fluids > 1)) then
            print '(A)', 'Unsupported combination of values of '// &
                'model_eqns and num_fluids. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (bubbles .and. (R0ref == dflt_real)) then
            print '(A)', 'Unsupported combination of values of '// &
                'bubbles and R0ref. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (bubbles .and. (nb == dflt_int)) then
            print '(a)', 'unsupported combination of values of '// &
                'bubbles and nb. '// &
                'exiting ...'
            call s_mpi_abort()
        elseif (bubbles .and. (thermal > 3)) then
            print '(a)', 'unsupported combination of values of '// &
                'bubbles and thermal. '// &
                'exiting ...'
            call s_mpi_abort()
        end if

        ! Constraint on the location of the case directory
        if (dir_check .neqv. .true.) then
            print '(A)', 'Unsupported choice for the value of '// &
                'case_dir.'
            print '(A)', 'WARNING: Ensure that compiler flags/choices in Makefiles match your compiler! '
            print '(A)', 'WARNING: Ensure that preprocessor flags are enabled! '
            call s_mpi_abort()

            ! Constraints on the use of a preexisting grid and initial condition
        elseif ((old_grid .neqv. .true.) .and. old_ic) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and old_ic. Exiting ...'
            call s_mpi_abort()

        elseif ((old_grid .or. old_ic) .and. t_step_old == dflt_int) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and old_ic and t_step_old. Exiting ...'
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
                'Total number of cells in azimuthal direction '// &
                'must be an even number. Exiting ...'
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

            ! Constraints on domain boundaries locations in the x-direction
        elseif ((old_grid .and. x_domain%beg /= dflt_real) &
                .or. &
                ((old_grid .neqv. .true.) .and. &
                 x_domain%beg == dflt_real)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and x_domain%beg. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif ((old_grid .and. x_domain%end /= dflt_real) &
                .or. &
                ((old_grid .neqv. .true.) .and. &
                 x_domain%end == dflt_real)) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and x_domain%end. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif ((old_grid .neqv. .true.) &
                .and. &
                x_domain%beg >= x_domain%end) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid, x_domain%beg and '// &
                'x_domain%end. Exiting ...'
            call s_mpi_abort()
        else if (qbmm .and. dist_type == dflt_int) then
            print '(A)', 'Dist type must be set if using QBMM. Exiting ...'
            call s_mpi_abort()
        else if (qbmm .and. (dist_type .ne. 1) .and. rhoRV > 0d0) then
            print '(A)', 'rhoRV cannot be used with dist_type \ne 1. Exiting ...'
            call s_mpi_abort()
        else if (polydisperse .and. R0_type == dflt_int) then
            print '(A)', 'R0 type must be set if using Polydisperse. Exiting ...'
            call s_mpi_abort()
        end if

        if (cyl_coord .neqv. .true.) then ! Cartesian coordinates

            ! Constraints on domain boundaries locations in the y-direction
            if ((n == 0 .and. y_domain%beg /= dflt_real) &
                .or. &
                (n > 0 &
                 .and. &
                 ((old_grid .and. y_domain%beg /= dflt_real) &
                  .or. &
                  ((old_grid .neqv. .true.) .and. &
                   y_domain%beg == dflt_real)))) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'values for old_grid, n and y_domain%beg. '// &
                    'Exiting ...'
                call s_mpi_abort()
            elseif ((n == 0 .and. y_domain%end /= dflt_real) &
                    .or. &
                    (n > 0 &
                     .and. &
                     ((old_grid .and. y_domain%end /= dflt_real) &
                      .or. &
                      ((old_grid .neqv. .true.) .and. &
                       y_domain%end == dflt_real)))) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'values for old_grid, n and y_domain%end. '// &
                    'Exiting ...'
                call s_mpi_abort()
            elseif (n > 0 &
                    .and. &
                    (old_grid .neqv. .true.) &
                    .and. &
                    y_domain%beg >= y_domain%end) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'values for old_grid, n, y_domain%beg and '// &
                    'y_domain%end. Exiting ...'
                call s_mpi_abort()

                ! Constraints on domain boundaries locations in the z-direction
            elseif ((p == 0 .and. z_domain%beg /= dflt_real) &
                    .or. &
                    (p > 0 &
                     .and. &
                     ((old_grid .and. z_domain%beg /= dflt_real) &
                      .or. &
                      ((old_grid .neqv. .true.) .and. &
                       z_domain%beg == dflt_real)))) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'values for old_grid, p and z_domain%beg. '// &
                    'Exiting ...'
                call s_mpi_abort()
            elseif ((p == 0 .and. z_domain%end /= dflt_real) &
                    .or. &
                    (p > 0 &
                     .and. &
                     ((old_grid .and. z_domain%end /= dflt_real) &
                      .or. &
                      ((old_grid .neqv. .true.) .and. &
                       z_domain%end == dflt_real)))) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'values for old_grid, p and z_domain%end. '// &
                    'Exiting ...'
                call s_mpi_abort()
            elseif (p > 0 &
                    .and. &
                    (old_grid .neqv. .true.) &
                    .and. &
                    z_domain%beg >= z_domain%end) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'values for old_grid, p, z_domain%beg and '// &
                    'z_domain%end. Exiting ...'
                call s_mpi_abort()
            end if

        else ! Cylindrical coordinates

            ! Constraints on domain boundaries for cylindrical coordinates
            if (n == 0 &
                .or. &
                y_domain%beg /= 0d0 &
                .or. &
                y_domain%end == dflt_real &
                .or. &
                y_domain%end < 0d0 &
                .or. &
                y_domain%beg >= y_domain%end) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'cyl_coord and n, y_domain%beg, or         '// &
                    'y_domain%end. Exiting ...'
                call s_mpi_abort()
            elseif ((p == 0 .and. z_domain%beg /= dflt_real) &
                    .or. &
                    (p == 0 .and. z_domain%end /= dflt_real)) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'cyl_coord and p, z_domain%beg, or '// &
                    'z_domain%end. Exiting ...'
                call s_mpi_abort()
            elseif (p > 0 .and. (z_domain%beg /= 0d0 &
                                 .or. &
                                 z_domain%end /= 2d0*pi)) then
                print '(A)', 'Unsupported choice of the combination of '// &
                    'cyl_coord and p, z_domain%beg, or '// &
                    'z_domain%end. Exiting ...'
                call s_mpi_abort()
            end if

        end if

        ! Constraints on the grid stretching in the x-direction
        if (old_grid .and. stretch_x) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and stretch_x. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x .and. a_x == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_x and a_x. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x .and. x_a == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_x and x_a. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x .and. x_b == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_x and x_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x .and. x_a >= x_b) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_x, x_a and x_b. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x &
                .and. &
                (a_x + log(cosh(a_x*(x_domain%beg - x_a))) &
                 + log(cosh(a_x*(x_domain%beg - x_b))) &
                 - 2d0*log(cosh(0.5d0*a_x*(x_b - x_a))))/a_x <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for x_domain%beg, stretch_x, a_x, '// &
                'x_a, and x_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_x &
                .and. &
                (a_x + log(cosh(a_x*(x_domain%end - x_a))) &
                 + log(cosh(a_x*(x_domain%end - x_b))) &
                 - 2d0*log(cosh(0.5d0*a_x*(x_b - x_a))))/a_x <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for x_domain%end, stretch_x, a_x, '// &
                'x_a, and x_b. Exiting ...'
            call s_mpi_abort()
        elseif (loops_x < 1) then
            print '(A)', 'Unsupported choice for the value of loops_x. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on the grid stretching in the y-direction
        elseif (old_grid .and. stretch_y) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and stretch_y. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (n == 0 .and. stretch_y) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for n and stretch_y. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y .and. a_y == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_y and a_y. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y .and. y_a == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_y and y_a. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y .and. y_b == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_y and y_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y .and. y_a >= y_b) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_y, y_a and y_b. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y &
                .and. &
                (a_y + log(cosh(a_y*(y_domain%beg - y_a))) &
                 + log(cosh(a_y*(y_domain%beg - y_b))) &
                 - 2d0*log(cosh(0.5d0*a_y*(y_b - y_a))))/a_y <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for y_domain%beg, stretch_y, a_y, '// &
                'y_a, and y_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_y &
                .and. &
                (a_y + log(cosh(a_y*(y_domain%end - y_a))) &
                 + log(cosh(a_y*(y_domain%end - y_b))) &
                 - 2d0*log(cosh(0.5d0*a_y*(y_b - y_a))))/a_y <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for y_domain%end, stretch_y, a_y, '// &
                'y_a, and y_b. Exiting ...'
            call s_mpi_abort()
        elseif (loops_y < 1) then
            print '(A)', 'Unsupported choice for the value of loops_y. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on the grid stretching in the z-direction
        elseif (old_grid .and. stretch_z) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for old_grid and stretch_z. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (cyl_coord .and. stretch_z) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for cyl_coord and stretch_z. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (p == 0 .and. stretch_z) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for p and stretch_z. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_z .and. a_z == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_z and a_z. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_z .and. z_a == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_z and z_a. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_z .and. z_b == dflt_real) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_z and z_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_z .and. z_a >= z_b) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for stretch_z, z_a and z_b. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (stretch_z &
                .and. &
                (a_z + log(cosh(a_z*(z_domain%beg - z_a))) &
                 + log(cosh(a_z*(z_domain%beg - z_b))) &
                 - 2d0*log(cosh(0.5d0*a_z*(z_b - z_a))))/a_z <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for z_domain%beg, stretch_z, a_z, '// &
                'z_a, and z_b. Exiting ...'
            call s_mpi_abort()
        elseif (stretch_z &
                .and. &
                (a_z + log(cosh(a_z*(z_domain%end - z_a))) &
                 + log(cosh(a_z*(z_domain%end - z_b))) &
                 - 2d0*log(cosh(0.5d0*a_z*(z_b - z_a))))/a_z <= 0d0) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for z_domain%end, stretch_z, a_z, '// &
                'z_a, and z_b. Exiting ...'
            call s_mpi_abort()
        elseif (loops_z < 1) then
            print '(A)', 'Unsupported choice for the value of loops_z. '// &
                'Exiting ...'
            call s_mpi_abort()

            ! Constraints on model equations and number of fluids in the flow
        elseif (all(model_eqns /= (/1, 2, 3, 4/))) then
            print '(A)', 'Unsupported value of model_eqns. Exiting ...'
            call s_mpi_abort()
        elseif (num_fluids /= dflt_int &
                .and. &
                (num_fluids < 1 .or. num_fluids > num_fluids_max)) then
            print '(A)', 'Unsupported value of num_fluids. Exiting ...'
            call s_mpi_abort()
        elseif ((model_eqns == 1 .and. num_fluids /= dflt_int) &
                .or. &
                (model_eqns == 2 .and. num_fluids == dflt_int)) then
            print '(A)', 'Unsupported combination of values of '// &
                'model_eqns and num_fluids. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (model_eqns == 1 .and. adv_alphan) then
            print '(A)', 'Unsupported combination of values of '// &
                'model_eqns and adv_alphan. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif (num_fluids == 1 .and. (adv_alphan .neqv. .true.)) then
            print '(A)', 'Unsupported combination of values of '// &
                'num_fluids and adv_alphan. '// &
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
        end if

        if (cyl_coord .neqv. .true.) then ! Cartesian coordinates

            ! Constraints on the boundary conditions in the y-direction
            if (bc_y%beg /= dflt_real &
                .and. &
                (bc_y%beg < -12 .or. bc_y%beg > -1)) then
                print '(A)', 'Unsupported choice for the value of '// &
                    'bc_y%beg. Exiting ...'
                call s_mpi_abort()
            elseif (bc_y%end /= dflt_real &
                    .and. &
                    (bc_y%end < -12 .or. bc_y%end > -1)) then
                print '(A)', 'Unsupported choice for the value of '// &
                    'bc_y%end. Exiting ...'
                call s_mpi_abort()
            elseif ((n == 0 .and. bc_y%beg /= dflt_real) &
                    .or. &
                    (n > 0 .and. bc_y%beg == dflt_real)) then
                print '(A)', 'Unsupported choice for the value of n and '// &
                    'bc_y%beg. Exiting ...'
                call s_mpi_abort()
            elseif ((n == 0 .and. bc_y%end /= dflt_real) &
                    .or. &
                    (n > 0 .and. bc_y%end == dflt_real)) then
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
            elseif (bc_z%beg /= dflt_real &
                    .and. &
                    (bc_z%beg < -12 .or. bc_z%beg > -1)) then
                print '(A)', 'Unsupported choice for the value of '// &
                    'bc_z%beg. Exiting ...'
                call s_mpi_abort()
            elseif (bc_z%end /= dflt_real &
                    .and. &
                    (bc_z%end < -12 .or. bc_z%end > -1)) then
                print '(A)', 'Unsupported choice for the value of '// &
                    'bc_z%end. Exiting ...'
                call s_mpi_abort()
            elseif ((p == 0 .and. bc_z%beg /= dflt_real) &
                    .or. &
                    (p > 0 .and. bc_z%beg == dflt_real)) then
                print '(A)', 'Unsupported choice for the value of p and '// &
                    'bc_z%beg. Exiting ...'
                call s_mpi_abort()
            elseif ((p == 0 .and. bc_z%end /= dflt_real) &
                    .or. &
                    (p > 0 .and. bc_z%end == dflt_real)) then
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

        else ! Cylindrical coordinates

            ! Constraints on the boundary conditions in the r-direction
            if (bc_y%beg /= dflt_real &
                .and. &
                ((p > 0 .and. bc_y%beg /= -13) &
                 .or. &
                 (p == 0 .and. bc_y%beg /= -2))) then
                print '(A)', 'Unsupported choice for the value of '// &
                    'bc_y%beg. Exiting ...'
                call s_mpi_abort()
            elseif (bc_y%end /= dflt_real &
                    .and. &
                    (bc_y%end < -12 .or. bc_y%end > -1)) then
                print '(A)', 'Unsupported choice for the value of '// &
                    'bc_y%end. Exiting ...'
                call s_mpi_abort()
            elseif ((n > 0 .and. bc_y%beg == dflt_real)) then
                print '(A)', 'Unsupported choice for the value of n and '// &
                    'bc_y%beg. Exiting ...'
                call s_mpi_abort()
            elseif ((n > 0 .and. bc_y%end == dflt_real)) then
                print '(A)', 'Unsupported choice for the value of n and '// &
                    'bc_y%end. Exiting ...'
                call s_mpi_abort()

                ! Constraints on the boundary conditions in the theta-direction
            elseif (bc_z%beg /= dflt_real &
                    .and. &
                    (bc_z%beg /= -1 .and. bc_z%beg /= -2)) then
                print '(A)', 'Unsupported choice for the value of '// &
                    'bc_z%beg. Exiting ...'
                call s_mpi_abort()
            elseif (bc_z%end /= dflt_real &
                    .and. &
                    (bc_z%end /= -1 .and. bc_z%end /= -2)) then
                print '(A)', 'Unsupported choice for the value of '// &
                    'bc_z%end. Exiting ...'
                call s_mpi_abort()
            elseif ((p == 0 .and. bc_z%beg /= dflt_real) &
                    .or. &
                    (p > 0 .and. bc_z%beg == dflt_real)) then
                print '(A)', 'Unsupported choice for the value of p and '// &
                    'bc_z%beg. Exiting ...'
                call s_mpi_abort()
            elseif ((p == 0 .and. bc_z%end /= dflt_real) &
                    .or. &
                    (p > 0 .and. bc_z%end == dflt_real)) then
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

        end if

        ! Constraints on number of patches making up the initial condition
        if (num_patches < 0 .or. num_patches > num_patches_max .or. &
            (num_patches == 0 .and. t_step_old == dflt_int)) then
            print '(A)', 'Unsupported choice for the value of '// &
                'num_patches. Exiting ...'
            call s_mpi_abort()
            ! Constraints on perturbing the initial condition
        elseif ((perturb_flow .and. perturb_flow_fluid == dflt_int) &
                .or. &
                ((perturb_flow .neqv. .true.) .and. (perturb_flow_fluid /= dflt_int))) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for perturb_flow and perturb_flow_fluid. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif ((perturb_flow_fluid > num_fluids) &
                .or. &
                (perturb_flow_fluid < 0 .and. perturb_flow_fluid /= dflt_int)) then
            print '(A)', 'Unsupported choice for the value of '// &
                'perturb_flow_fluid. Exiting ...'
            call s_mpi_abort()
        elseif ((perturb_sph .and. perturb_sph_fluid == dflt_int) &
                .or. &
                ((perturb_sph .neqv. .true.) .and. (perturb_sph_fluid /= dflt_int))) then
            print '(A)', 'Unsupported choice of the combination of '// &
                'values for perturb_sph and perturb_sph_fluid. '// &
                'Exiting ...'
            call s_mpi_abort()
        elseif ((perturb_sph_fluid > num_fluids) &
                .or. &
                (perturb_sph_fluid < 0 .and. perturb_sph_fluid /= dflt_int)) then
            print '(A)', 'Unsupported choice for the value of '// &
                'perturb_sph_fluid. Exiting ...'
            call s_mpi_abort()
        elseif ((any(fluid_rho /= dflt_real)) .and. (perturb_sph .neqv. .true.)) then
            print '(A)', 'Unsupported choices for values of perturb_sph '// &
                'and fluid_rho. Exiting ...'
            call s_mpi_abort()
        end if

        if (perturb_sph) then
            do i = 1, num_fluids
                if (fluid_rho(i) == dflt_real) then
                    print '(A,I0,A)', 'Unsupported choice for value of fluid_rho(', &
                        i, '). Exiting ...'
                    call s_mpi_abort()
                end if
            end do
        end if

        ! Constraints on the geometric initial condition patch parameters
        do i = 1, num_patches_max
            if (i <= num_patches) then
                if (patch_icpp(i)%geometry == 1) then
                    call s_check_line_segment_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 2) then
                    call s_check_circle_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 3) then
                    call s_check_rectangle_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 4) then
                    call s_check_line_sweep_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 5) then
                    call s_check_ellipse_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 6) then
                    call s_check_isentropic_vortex_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 7) then
                    call s_check_2D_analytical_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 8) then
                    call s_check_sphere_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 9) then
                    call s_check_cuboid_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 10) then
                    call s_check_cylinder_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 11) then
                    call s_check_plane_sweep_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 12) then
                    call s_check_ellipsoid_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 13) then
                    call s_check_3D_analytical_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 14) then
                    call s_check_spherical_harmonic_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 15) then
                    call s_check_1d_analytical_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 16) then
                    print *, '1d pressure sinusoid'
                elseif (patch_icpp(i)%geometry == 17) then
                    print *, '2d spiral'
                elseif (patch_icpp(i)%geometry == 18) then
                    print *, '2d var circle'
                elseif (patch_icpp(i)%geometry == 19) then
                    print *, '3d var circle'
                else
                    print '(A,I0,A)', 'Unsupported choice of the '// &
                        'geometry of active patch ', i, &
                        ' detected. Exiting ...'
                    call s_mpi_abort()
                end if
            else
                if (patch_icpp(i)%geometry == dflt_int) then
                    call s_check_inactive_patch_geometry(i)
                else
                    print '(A,I0,A)', 'Unsupported choice of the '// &
                        'geometry of inactive patch ', i, &
                        ' detected. Exiting ...'
                    call s_mpi_abort()
                end if
            end if
        end do

        ! Constraints on overwrite rights initial condition patch parameters
        do i = 1, num_patches_max
            if (i <= num_patches) then
                call s_check_active_patch_alteration_rights(i)
            else
                call s_check_inactive_patch_alteration_rights(i)
            end if
        end do

        ! Constraints on smoothing initial condition patch parameters
        do i = 1, num_patches_max
            if (i > 1 .and. (patch_icpp(i)%geometry == 2 .or. &
                             patch_icpp(i)%geometry == 4 .or. &
                             patch_icpp(i)%geometry == 5 .or. &
                             patch_icpp(i)%geometry == 8 .or. &
                             patch_icpp(i)%geometry == 10 .or. &
                             patch_icpp(i)%geometry == 11 .or. &
                             patch_icpp(i)%geometry == 12)) then
                call s_check_supported_patch_smoothing(i)
            else
                call s_check_unsupported_patch_smoothing(i)
            end if
        end do

        ! Constraints on flow variables initial condition patch parameters
        do i = 1, num_patches_max
            if (i <= num_patches) then
                call s_check_active_patch_primitive_variables(i)
            else
                call s_check_inactive_patch_primitive_variables(i)
            end if
        end do

        ! Constraints on the stiffened equation of state fluids parameters
        do i = 1, num_fluids_max

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

    end subroutine s_check_input_file ! ------------------------------------

    !> This subroutine verifies that the geometric parameters of
        !!      the line segment patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_line_segment_patch_geometry(patch_id) ! -------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the line segment patch
        if (n > 0 .or. patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            cyl_coord) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of line segment '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_line_segment_patch_geometry ! -------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the circle patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_circle_patch_geometry(patch_id) ! -------------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the circle patch
        if (n == 0 .or. p > 0 .or. patch_icpp(patch_id)%radius <= 0d0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of circle '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_circle_patch_geometry ! -------------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the rectangle patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_rectangle_patch_geometry(patch_id) ! ----------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the rectangle patch
        if (n == 0 .or. p > 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_y <= 0d0) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of rectangle '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_rectangle_patch_geometry ! ----------------------

    !> This subroutine verifies that the geometric parameters of
        !!      the line sweep patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_line_sweep_patch_geometry(patch_id) ! ---------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the line sweep patch
        if (n == 0 .or. p > 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(1) == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(2) == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(3) /= dflt_real) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of line sweep '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_line_sweep_patch_geometry ! ---------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the ellipse patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_ellipse_patch_geometry(patch_id) ! ------------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the ellipse patch
        if (n == 0 .or. p > 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(1) == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(2) == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(3) /= dflt_real) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of ellipse '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_ellipse_patch_geometry ! ------------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the isentropic vortex patch have been entered by the user
        !!      consistently.
        !!  @param patch_id Patch identifier
    subroutine s_check_isentropic_vortex_patch_geometry(patch_id) ! --------

        integer, intent(IN) :: patch_id

        ! Constraints on the isentropic vortex patch geometric parameters
        if (n == 0 .or. p > 0 .or. model_eqns == 2 &
            .or. &
            patch_icpp(patch_id)%radius <= 0d0 &
            .or. &
            patch_icpp(patch_id)%epsilon <= 0d0 &
            .or. &
            patch_icpp(patch_id)%beta <= 0d0) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of isentropic '// &
                'vortex patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_isentropic_vortex_patch_geometry ! --------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the analytical patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_1D_analytical_patch_geometry(patch_id) ! ---------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the analytical patch
        if (n > 0 .or. p > 0 &
            .or. &
            (model_eqns .ne. 4 .and. model_eqns .ne. 2) &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of 1D analytical '// &
                'patch ', patch_id, '. Exiting...'
            call s_mpi_abort()
        end if
    end subroutine s_check_1D_analytical_patch_geometry ! ---------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the analytical patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_2D_analytical_patch_geometry(patch_id) ! ---------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the analytical patch
        if (n == 0 .or. p > 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_y <= 0d0) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of 2D analytical '// &
                'patch ', patch_id, '. Exiting...'
            call s_mpi_abort()
        end if
    end subroutine s_check_2D_analytical_patch_geometry ! ---------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the analytical patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_3D_analytical_patch_geometry(patch_id) ! ---------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the analytical patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_y <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_z <= 0d0) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of 3D analytical '// &
                'patch ', patch_id, '. Exiting...'
            call s_mpi_abort()
        end if
    end subroutine s_check_3D_analytical_patch_geometry ! ---------------------

    !> This subroutine verifies that the geometric parameters of
        !!      the sphere patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_sphere_patch_geometry(patch_id) ! -------------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the sphere patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%radius <= 0d0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of sphere '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_sphere_patch_geometry ! -------------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the spherical harmonic  patch have consistently been
        !!      inputted by the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_spherical_harmonic_patch_geometry(patch_id) ! -------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the spherical harmonic patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%radius <= 0d0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            all(patch_icpp(patch_id)%epsilon /= (/1d0, 2d0, 3d0, 4d0, 5d0/)) &
            .or. &
            patch_icpp(patch_id)%beta < 0d0 &
            .or. &
            patch_icpp(patch_id)%beta > patch_icpp(patch_id)%epsilon) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of spherical '// &
                'harmonic patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_spherical_harmonic_patch_geometry ! -------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the cuboid patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_cuboid_patch_geometry(patch_id) ! -------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the cuboid patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_y <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_z <= 0d0) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of cuboid '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_cuboid_patch_geometry ! -------------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the cylinder patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_cylinder_patch_geometry(patch_id) ! -----------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the cylinder patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            (patch_icpp(patch_id)%length_x <= 0d0 .and. &
             patch_icpp(patch_id)%length_y <= 0d0 .and. &
             patch_icpp(patch_id)%length_z <= 0d0) &
            .or. &
            (patch_icpp(patch_id)%length_x > 0d0 .and. &
             (patch_icpp(patch_id)%length_y /= dflt_real .or. &
              patch_icpp(patch_id)%length_z /= dflt_real)) &
            .or. &
            (patch_icpp(patch_id)%length_y > 0d0 .and. &
             (patch_icpp(patch_id)%length_x /= dflt_real .or. &
              patch_icpp(patch_id)%length_z /= dflt_real)) &
            .or. &
            (patch_icpp(patch_id)%length_z > 0d0 .and. &
             (patch_icpp(patch_id)%length_x /= dflt_real .or. &
              patch_icpp(patch_id)%length_y /= dflt_real)) &
            .or. &
            patch_icpp(patch_id)%radius <= 0d0) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of cylinder '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_cylinder_patch_geometry ! -----------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the plane sweep patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_plane_sweep_patch_geometry(patch_id) ! --------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the plane sweep patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(1) == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(2) == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(3) == dflt_real) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of plane sweep '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_plane_sweep_patch_geometry ! --------------------

    !> This subroutine verifies that the geometric parameters of
        !!      the ellipsoid patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_ellipsoid_patch_geometry(patch_id) ! ----------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the ellipsoid patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(1) == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(2) == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(3) == dflt_real) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of ellipsoid '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_ellipsoid_patch_geometry ! ----------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the inactive patch remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_inactive_patch_geometry(patch_id) ! -----------------

        integer, intent(IN) :: patch_id

        ! Constraints on the geometric parameters of the inactive patch
        if (patch_icpp(patch_id)%x_centroid /= dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid /= dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid /= dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x /= dflt_real &
            .or. &
            patch_icpp(patch_id)%length_y /= dflt_real &
            .or. &
            patch_icpp(patch_id)%length_z /= dflt_real &
            .or. &
            patch_icpp(patch_id)%radius /= dflt_real &
            .or. &
            patch_icpp(patch_id)%epsilon /= dflt_real &
            .or. &
            patch_icpp(patch_id)%beta /= dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(1) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(2) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(3) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(1) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(2) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(3) /= dflt_real) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'geometric parameters of inactive '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_inactive_patch_geometry ! -----------------------

    !>  This subroutine verifies that any rights granted to the
        !!      given active patch, to overwrite the preceding active
        !!      patches, were consistently inputted by the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_active_patch_alteration_rights(patch_id) ! ----------

        integer, intent(IN) :: patch_id

        ! Constraints on the alteration rights of an active patch
        if (patch_icpp(patch_id)%alter_patch(0) .eqv. .false. &
            .or. &
            any(patch_icpp(patch_id)%alter_patch(patch_id:))) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'alteration rights of active '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_active_patch_alteration_rights ! ----------------

    !>  This subroutine verifies that the rights of the given
        !!      inactive patch, to overwrite the preceding patches,
        !!      remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_inactive_patch_alteration_rights(patch_id) ! --------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Constraints on the alteration rights of an inactive patch
        if (patch_icpp(patch_id)%alter_patch(0) .eqv. .false. &
            .or. &
            any(patch_icpp(patch_id)%alter_patch(1:))) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'alteration rights of inactive '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_inactive_patch_alteration_rights ! --------------

    !> This subroutine verifies that the smoothing parameters of
        !!      the given patch, which supports the smoothing out of its
        !!      boundaries, have consistently been inputted by the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_supported_patch_smoothing(patch_id) ! ---------------

        integer, intent(IN) :: patch_id

        ! Constraints on the smoothing parameters of a supported patch
        if ((patch_icpp(patch_id)%smoothen &
             .and. &
             (patch_icpp(patch_id)%smooth_patch_id >= patch_id &
              .or. &
              patch_icpp(patch_id)%smooth_patch_id == 0 &
              .or. &
              patch_icpp(patch_id)%smooth_coeff <= 0d0)) &
            .or. &
            ((patch_icpp(patch_id)%smoothen .neqv. .true.) &
             .and. &
             (patch_icpp(patch_id)%smooth_patch_id /= patch_id &
              .or. &
              patch_icpp(patch_id)%smooth_coeff /= dflt_real))) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'smoothing parameters of supported '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_supported_patch_smoothing ! ---------------------

    !> This subroutine verifies that the smoothing parameters of
        !!      the given patch, which does not support the smoothing out
        !!          of its boundaries, remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_unsupported_patch_smoothing(patch_id) ! -------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Constraints on the smoothing parameters of an unsupported patch
        if (patch_icpp(patch_id)%smoothen &
            .or. &
            patch_icpp(patch_id)%smooth_patch_id /= patch_id &
            .or. &
            patch_icpp(patch_id)%smooth_coeff /= dflt_real) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'smoothing parameters of unsupported '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_unsupported_patch_smoothing ! -------------------

    !>  This subroutine verifies that the primitive variables
        !!      associated with the given active patch are physically
        !!      consistent.
        !!  @param patch_id Patch identifier
    subroutine s_check_active_patch_primitive_variables(patch_id) ! --------

        integer, intent(IN) :: patch_id

        ! Constraints on the primitive variables of an active patch
        if (patch_icpp(patch_id)%vel(1) == dflt_real &
            .or. &
            (n == 0 .and. patch_icpp(patch_id)%vel(2) /= dflt_real) &
            .or. &
            (n > 0 .and. patch_icpp(patch_id)%vel(2) == dflt_real) &
            .or. &
            (p == 0 .and. patch_icpp(patch_id)%vel(3) /= dflt_real) &
            .or. &
            (p > 0 .and. patch_icpp(patch_id)%vel(3) == dflt_real) &
            !                               .OR.                           &
            !                 patch_icpp(patch_id)%pres <= 0d0             &
            .or. &
            (model_eqns == 1 .and. &
             (patch_icpp(patch_id)%rho <= 0d0 .or. &
              patch_icpp(patch_id)%gamma <= 0d0 .or. &
              patch_icpp(patch_id)%pi_inf < 0d0)) &
            .or. &
            (patch_icpp(patch_id)%geometry == 5 &
             .and. &
             patch_icpp(patch_id)%pi_inf > 0) &
            .or. &
            (model_eqns == 2 &
             .and. &
             (any(patch_icpp(patch_id)%alpha_rho(1:num_fluids) < 0d0) &
              !                         .OR.                           &
              !ANY(patch_icpp(patch_id)%alpha(1: adv_idx%end           &
              !                                 -   E_idx    ) < 0d0)  &
              !                         .OR.                           &
              !SUM(patch_icpp(patch_id)%alpha(1: adv_idx%end           &
              !                                 -   E_idx    ))> 1d0)
              ))) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'primitive variables of active '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

        if (model_eqns == 2 .and. num_fluids < num_fluids_max) then

            if (any(patch_icpp(patch_id)%alpha_rho(num_fluids + 1:) &
                    /= dflt_real) &
                .or. &
                any(patch_icpp(patch_id)%alpha(num_fluids + 1:) &
                    /= dflt_real) &
                .or. &
                ((adv_alphan .neqv. .true.) &
                 .and. &
                 patch_icpp(patch_id)%alpha(num_fluids) /= dflt_real) &
                .or. &
                (adv_alphan &
                 .and. &
                 patch_icpp(patch_id)%alpha(num_fluids) == dflt_real)) then

                print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                    'primitive variables of active '// &
                    'patch ', patch_id, '. Exiting ...'
                !CALL s_mpi_abort()

            end if

        end if

    end subroutine s_check_active_patch_primitive_variables ! --------------

    !>  This subroutine verifies that the primitive variables
        !!      associated with the given inactive patch remain unaltered
        !!      by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_inactive_patch_primitive_variables(patch_id) ! ------

        integer, intent(IN) :: patch_id

        ! Constraints on the primitive variables of an inactive patch
        if (any(patch_icpp(patch_id)%alpha_rho /= dflt_real) &
            .or. &
            patch_icpp(patch_id)%rho /= dflt_real &
            .or. &
            any(patch_icpp(patch_id)%vel /= dflt_real) &
            .or. &
            patch_icpp(patch_id)%pres /= dflt_real &
            .or. &
            any(patch_icpp(patch_id)%alpha /= dflt_real) &
            .or. &
            patch_icpp(patch_id)%gamma /= dflt_real &
            .or. &
            patch_icpp(patch_id)%pi_inf /= dflt_real) then

            print '(A,I0,A)', 'Inconsistency(ies) detected in '// &
                'primitive variables of inactive '// &
                'patch ', patch_id, '. Exiting ...'
            call s_mpi_abort()

        end if

    end subroutine s_check_inactive_patch_primitive_variables ! ------------

    !> The goal of this subroutine is to read in any preexisting
        !!      grid data as well as based on the imported grid, complete
        !!      the necessary global computational domain parameters.
        !! @param dflt_int Default null integer
    subroutine s_read_serial_grid_data_files(dflt_int) ! ---

        integer, intent(IN) :: dflt_int

        ! Generic string used to store the address of a particular file
        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc

        ! Logical variable used to test the existence of folders
        logical :: dir_check

        ! Generic logical used for the purpose of asserting whether a file
        ! is or is not present in the designated location
        logical :: file_check

        ! Setting address of the local processor rank and time-step directory
        write (proc_rank_dir, '(A,I0)') '/p', proc_rank
        proc_rank_dir = trim(case_dir)//trim(proc_rank_dir)

        write (t_step_dir, '(A,I0)') '/', t_step_old
        t_step_dir = trim(proc_rank_dir)//trim(t_step_dir)

        ! Inquiring as to the existence of the time-step directory
        file_loc = trim(t_step_dir)//'/.'
        call my_inquire(file_loc, dir_check)

        ! If the time-step directory is missing, the pre-process exits
        if (dir_check .neqv. .true.) then
            print '(A)', 'Time-step folder '//trim(t_step_dir)// &
                ' is missing. Exiting ...'
            call s_mpi_abort()
        end if

        ! Reading the Grid Data File for the x-direction ===================

        ! Checking whether x_cb.dat exists
        file_loc = trim(t_step_dir)//'/x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! If it exists, x_cb.dat is read
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS='old', ACTION='read')
            read (1) x_cb
            close (1)
        else
            print '(A)', 'File x_cb.dat is missing in '// &
                trim(t_step_dir)//'. Exiting ...'
            call s_mpi_abort()
        end if

        ! Computing cell-center locations
        x_cc = (x_cb(0:m) + x_cb(-1:(m - 1)))/2d0

        ! Computing minimum cell-width
        dx = minval(x_cb(0:m) - x_cb(-1:m - 1))
        if (num_procs > 1) call s_mpi_reduce_min(dx)

        ! Setting locations of domain bounds
        x_domain%beg = x_cb(-1)
        x_domain%end = x_cb(m)

        ! ==================================================================

        ! Reading the Grid Data File for the y-direction ===================

        if (n > 0) then

            ! Checking whether y_cb.dat exists
            file_loc = trim(t_step_dir)//'/y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            ! If it exists, y_cb.dat is read
            if (file_check) then
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS='old', ACTION='read')
                read (1) y_cb
                close (1)
            else
                print '(A)', 'File y_cb.dat is missing in '// &
                    trim(t_step_dir)//'. Exiting ...'
                call s_mpi_abort()
            end if

            ! Computing cell-center locations
            y_cc = (y_cb(0:n) + y_cb(-1:(n - 1)))/2d0

            ! Computing minimum cell-width
            dy = minval(y_cb(0:n) - y_cb(-1:n - 1))
            if (num_procs > 1) call s_mpi_reduce_min(dy)

            ! Setting locations of domain bounds
            y_domain%beg = y_cb(-1)
            y_domain%end = y_cb(n)

            ! ==================================================================

            ! Reading the Grid Data File for the z-direction ===================

            if (p > 0) then

                ! Checking whether z_cb.dat exists
                file_loc = trim(t_step_dir)//'/z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_check)

                ! If it exists, z_cb.dat is read
                if (file_check) then
                    open (1, FILE=trim(file_loc), FORM='unformatted', &
                          STATUS='old', ACTION='read')
                    read (1) z_cb
                    close (1)
                else
                    print '(A)', 'File z_cb.dat is missing in '// &
                        trim(t_step_dir)//'. Exiting ...'
                    call s_mpi_abort()
                end if

                ! Computing cell-center locations
                z_cc = (z_cb(0:p) + z_cb(-1:(p - 1)))/2d0

                ! Computing minimum cell-width
                dz = minval(z_cb(0:p) - z_cb(-1:p - 1))
                if (num_procs > 1) call s_mpi_reduce_min(dz)

                ! Setting locations of domain bounds
                z_domain%beg = z_cb(-1)
                z_domain%end = z_cb(p)

            end if

        end if

        ! ==================================================================

        ! If only the preexisting grid data files are read in and there will
        ! not be any preexisting initial condition data files imported, then
        ! the directory associated with the rank of the local processor may
        ! be cleaned to make room for the new pre-process data. In addition,
        ! the time-step directory that will contain the new grid and initial
        ! condition data are also generated.
        if (old_ic .neqv. .true.) then
            call SYSTEM('rm -rf '//trim(proc_rank_dir)//'/*')
            call SYSTEM('mkdir -p '//trim(proc_rank_dir)//'/0')
        end if

    end subroutine s_read_serial_grid_data_files ! --------------------------------

    !> Cell-boundary data are checked for consistency by looking
        !!      at the (non-)uniform cell-width distributions for all the
        !!      active coordinate directions and making sure that all of
        !!      the cell-widths are positively valued
    subroutine s_check_grid_data_files() ! -----------------

        ! Cell-boundary Data Consistency Check in x-direction ==============

        if (any(x_cb(0:m) - x_cb(-1:m - 1) <= 0d0)) then
            print '(A)', 'x_cb.dat in '//trim(t_step_dir)// &
                ' contains non-positive cell-spacings. Exiting ...'
            call s_mpi_abort()
        end if

        ! ==================================================================

        ! Cell-boundary Data Consistency Check in y-direction ==============

        if (n > 0) then

            if (any(y_cb(0:n) - y_cb(-1:n - 1) <= 0d0)) then
                print '(A)', 'y_cb.dat in '//trim(t_step_dir)// &
                    ' contains non-positive cell-spacings. '// &
                    'Exiting ...'
                call s_mpi_abort()
            end if

            ! ==================================================================

            ! Cell-boundary Data Consistency Check in z-direction ==============

            if (p > 0) then

                if (any(z_cb(0:p) - z_cb(-1:p - 1) <= 0d0)) then
                    print '(A)', 'z_cb.dat in '//trim(t_step_dir)// &
                        ' contains non-positive cell-spacings'// &
                        ' .Exiting ...'
                    call s_mpi_abort()
                end if

            end if

        end if

        ! ==================================================================

    end subroutine s_check_grid_data_files ! -------------------------------

    !> The goal of this subroutine is to read in any preexisting
        !!      initial condition data files so that they may be used by
        !!      the pre-process as a starting point in the creation of an
        !!      all new initial condition.
        !! @param q_cons_vf Conservative variables
    subroutine s_read_serial_ic_data_files(q_cons_vf) ! ---------------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc !<
        ! Generic string used to store the address of a particular file

        character(LEN= &
                  int(floor(log10(real(sys_size, kind(0d0))))) + 1) :: file_num !<
            !! Used to store the variable position, in character form, of the
            !! currently manipulated conservative variable file

        logical :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location

        integer :: i !< Generic loop iterator

        ! Reading the Conservative Variables Data Files ====================
        do i = 1, sys_size

            ! Checking whether data file associated with variable position
            ! of the currently manipulated conservative variable exists
            write (file_num, '(I0)') i
            file_loc = trim(t_step_dir)//'/q_cons_vf'// &
                       trim(file_num)//'.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            ! If it exists, the data file is read
            if (file_check) then
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS='old', ACTION='read')
                read (1) q_cons_vf(i)%sf
                close (1)
            else
                print '(A)', 'File q_cons_vf'//trim(file_num)// &
                    '.dat is missing in '//trim(t_step_dir)// &
                    '. Exiting ...'
                call s_mpi_abort()
            end if

        end do

        ! ==================================================================

        ! Since the preexisting grid and initial condition data files have
        ! been read in, the directory associated with the rank of the local
        ! process may be cleaned out to make room for new pre-process data.
        ! In addition, the time-step folder that will contain the new grid
        ! and initial condition data are also generated.
        call SYSTEM('rm -rf '//trim(proc_rank_dir)//'/*')
        call SYSTEM('mkdir -p '//trim(proc_rank_dir)//'/0')

    end subroutine s_read_serial_ic_data_files ! ----------------------------------

    !> Cell-boundary data are checked for consistency by looking
        !!      at the (non-)uniform cell-width distributions for all the
        !!      active coordinate directions and making sure that all of
        !!      the cell-widths are positively valued
        !! @param dflt_int Default null integer
    subroutine s_read_parallel_grid_data_files(dflt_int)

        integer, intent(IN) :: dflt_int

        real(kind(0d0)), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        ! Read in cell boundary locations in x-direction
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ_ALL(ifile, x_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            print '(A)', 'File ', trim(file_loc), ' is missing. Exiting... '
            call s_mpi_abort()
        end if

        ! Assigning local cell boundary locations
        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        ! Computing cell center locations
        x_cc(0:m) = (x_cb(0:m) + x_cb(-1:(m - 1)))/2d0
        ! Computing minimum cell width
        dx = minval(x_cb(0:m) - x_cb(-1:(m - 1)))
        if (num_procs > 1) call s_mpi_reduce_min(dx)
        ! Setting locations of domain bounds
        x_domain%beg = x_cb(-1)
        x_domain%end = x_cb(m)

        if (n > 0) then
            ! Read in cell boundary locations in y-direction
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ_ALL(ifile, y_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                print '(A)', 'File ', trim(file_loc), ' is missing. Exiting... '
                call s_mpi_abort()
            end if

            ! Assigning local cell boundary locations
            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            ! Computing cell center locations
            y_cc = (y_cb(0:n) + y_cb(-1:(n - 1)))/2d0
            ! Computing minimum cell width
            dy = minval(y_cb(0:n) - y_cb(-1:(n - 1)))
            if (num_procs > 1) call s_mpi_reduce_min(dy)
            ! Setting locations of domain bounds
            y_domain%beg = y_cb(-1)
            y_domain%end = y_cb(n)

            if (p > 0) then
                ! Read in cell boundary locations in z-direction
                file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ_ALL(ifile, z_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    print '(A)', 'File ', trim(file_loc), ' is missing. Exiting... '
                    call s_mpi_abort()
                end if

                ! Assigning local cell boundary locations
                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                ! Computing cell center locations
                z_cc = (z_cb(0:p) + z_cb(-1:(p - 1)))/2d0
                ! Computing minimum cell width
                dz = minval(z_cb(0:p) - z_cb(-1:(p - 1)))
                if (num_procs > 1) call s_mpi_reduce_min(dz)
                ! Setting locations of domain bounds
                z_domain%beg = z_cb(-1)
                z_domain%end = z_cb(p)

            end if
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

    end subroutine s_read_parallel_grid_data_files ! -----------------------

    !> The goal of this subroutine is to read in any preexisting
        !!      initial condition data files so that they may be used by
        !!      the pre-process as a starting point in the creation of an
        !!      all new initial condition.
        !! @param q_cons_vf Conservative variables
    subroutine s_read_parallel_ic_data_files(q_cons_vf) ! ------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        integer :: i

        ! Open the file to read
        write (file_loc, '(I0,A)') t_step_old, '.dat'
        file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Initialize MPI data I/O
            call s_initialize_mpi_data(q_cons_vf)

            ! Size of local arrays
            data_size = (m + 1)*(n + 1)*(p + 1)

            ! Resize some integers so MPI can read even the biggest files
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8d0, MPI_OFFSET_KIND)
            MOK = int(1d0, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            ! Read the data for each variable
            do i = 1, adv_idx%end
                var_MOK = int(i, MPI_OFFSET_KIND)

                ! Initial displacement to skip at beginning of file
                disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                       'native', mpi_info_int, ierr)
                call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                   MPI_DOUBLE_PRECISION, status, ierr)
            end do

            call s_mpi_barrier()

            call MPI_FILE_CLOSE(ifile, ierr)

        else
            print '(A)', 'File ', trim(file_loc), ' is missing. Exiting... '
            call s_mpi_abort()
        end if
        call s_mpi_barrier()
        if (proc_rank == 0) call SYSTEM('rm -f '//trim(file_loc))

    end subroutine s_read_parallel_ic_data_files ! -------------------------

    subroutine s_initialize_start_up_module() !-----------------------------

        if (parallel_io .neqv. .true.) then
            s_read_grid_data_files => s_read_serial_grid_data_files
            s_read_ic_data_files => s_read_serial_ic_data_files
        else
            s_read_grid_data_files => s_read_parallel_grid_data_files
            s_read_ic_data_files => s_read_parallel_ic_data_files
        end if

    end subroutine s_initialize_start_up_module ! --------------------------

    subroutine s_finalize_start_up_module() ! ------------------------------

        s_read_grid_data_files => null()
        s_read_ic_data_files => null()

    end subroutine s_finalize_start_up_module ! ----------------------------

end module m_start_up
