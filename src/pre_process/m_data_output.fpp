!>
!! @file m_data_output.f90
!! @brief Contains module m_data_output

!> @brief This module takes care of writing the grid and initial condition
!!              data files into the "0" time-step directory located in the folder
!!              associated with the rank of the local processor, which is a sub-
!!              directory of the case folder specified by the user in the input
!!              file pre_process.inp.
module m_data_output

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_helper

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_compile_specific

    use m_variables_conversion

    use m_helper

    use m_delay_file_access

    use m_boundary_common

    use m_boundary_conditions

    use m_thermochem, only: species_names

    use m_helper

    implicit none

    private; 
    public :: s_write_serial_data_files, &
              s_write_parallel_data_files, &
              s_write_data_files, &
              s_initialize_data_output_module, &
              s_finalize_data_output_module

    type(scalar_field), allocatable, dimension(:) :: q_cons_temp

    abstract interface

        !>  Interface for the conservative data
        !! @param q_cons_vf Conservative variables
        !! @param ib_markers track if a cell is within the immersed boundary
        !! @param levelset closest distance from every cell to the IB
        !! @param levelset_norm normalized vector from every cell to the closest point to the IB
        impure subroutine s_write_abstract_data_files(q_cons_vf, q_prim_vf, bc_type, ib_markers, levelset, levelset_norm)

            import :: scalar_field, integer_field, sys_size, m, n, p, &
                pres_field, levelset_field, levelset_norm_field, num_dims

            ! Conservative variables
            type(scalar_field), &
                dimension(sys_size), &
                intent(inout) :: q_cons_vf, q_prim_vf

            type(integer_field), &
                dimension(1:num_dims, -1:1), &
                intent(in) :: bc_type

            ! IB markers
            type(integer_field), &
                intent(in), optional :: ib_markers

            ! Levelset
            type(levelset_field), &
                intent(IN), optional :: levelset

            ! Levelset Norm
            type(levelset_norm_field), &
                intent(IN), optional :: levelset_norm

        end subroutine s_write_abstract_data_files
    end interface

    character(LEN=path_len + 2*name_len), private :: t_step_dir !<
    !! Time-step folder into which grid and initial condition data will be placed

    character(LEN=path_len + 2*name_len), public :: restart_dir !<
    !! Restart data folder

    procedure(s_write_abstract_data_files), pointer :: s_write_data_files => null()

contains

    !>  Writes grid and initial condition data files to the "0"
        !!  time-step directory in the local processor rank folder
        !! @param q_cons_vf Conservative variables
        !! @param ib_markers track if a cell is within the immersed boundary
        !! @param levelset closest distance from every cell to the IB
        !! @param levelset_norm normalized vector from every cell to the closest point to the IB
    impure subroutine s_write_serial_data_files(q_cons_vf, q_prim_vf, bc_type, ib_markers, levelset, levelset_norm)
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_cons_vf, q_prim_vf

        ! BC types
        type(integer_field), &
            dimension(1:num_dims, -1:1), &
            intent(in) :: bc_type

        ! IB markers
        type(integer_field), &
            intent(in), optional :: ib_markers

        ! Levelset
        type(levelset_field), &
            intent(IN), optional :: levelset

        ! Levelset Norm
        type(levelset_norm_field), &
            intent(IN), optional :: levelset_norm

        logical :: file_exist !< checks if file exists

        character(LEN=15) :: FMT
        character(LEN=3) :: status

        character(LEN= &
                  int(floor(log10(real(sys_size, wp)))) + 1) :: file_num !< Used to store
            !! the number, in character form, of the currently
            !! manipulated conservative variable data file

        character(LEN=len_trim(t_step_dir) + name_len) :: file_loc !<
            !! Generic string used to store the address of a particular file

        integer :: i, j, k, l, r, c !< Generic loop iterator
        integer :: t_step

        real(wp), dimension(nb) :: nRtmp         !< Temporary bubble concentration
        real(wp) :: nbub                         !< Temporary bubble number density
        real(wp) :: gamma, lit_gamma, pi_inf, qv !< Temporary EOS params
        real(wp) :: rho                          !< Temporary density
        real(wp) :: pres, T                         !< Temporary pressure

        real(wp) :: rhoYks(1:num_species) !< Temporary species mass fractions

        real(wp) :: pres_mag

        pres_mag = 0._wp

        T = dflt_T_guess

        t_step = 0

        ! Outputting the Locations of the Cell-boundaries

        if (old_grid) then
            status = 'old'
        else
            status = 'new'
        end if

        if (bc_io) then
            if (igr) then
                call s_write_serial_boundary_condition_files(q_cons_vf, bc_type, t_step_dir, old_grid)
            else
                call s_write_serial_boundary_condition_files(q_prim_vf, bc_type, t_step_dir, old_grid)
            end if
        end if

        ! x-coordinate direction
        file_loc = trim(t_step_dir)//'/x_cb.dat'
        open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
        write (1) x_cb(-1:m)
        close (1)

        ! y- and z-coordinate directions
        if (n > 0) then
            ! y-coordinate direction
            file_loc = trim(t_step_dir)//'/y_cb.dat'
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS=status)
            write (1) y_cb(-1:n)
            close (1)

            ! z-coordinate direction
            if (p > 0) then
                file_loc = trim(t_step_dir)//'/z_cb.dat'
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS=status)
                write (1) z_cb(-1:p)
                close (1)
            end if
        end if

        if (ib) then

            ! Outputting IB Markers
            file_loc = trim(t_step_dir)//'/ib.dat'

            open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
            write (1) ib_markers%sf
            close (1)

            do i = 1, num_ibs
                if (patch_ib(i)%geometry == 4) then

                    file_loc = trim(t_step_dir)//'/airfoil_u.dat'

                    open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
                    write (1) airfoil_grid_u(1:Np)
                    close (1)

                    file_loc = trim(t_step_dir)//'/airfoil_l.dat'

                    open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
                    write (1) airfoil_grid_l(1:Np)
                    close (1)
                end if
            end do

            ! Outtputting Levelset Info
            file_loc = trim(t_step_dir)//'/levelset.dat'

            open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
            write (1) levelset%sf
            close (1)

            file_loc = trim(t_step_dir)//'/levelset_norm.dat'

            open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
            write (1) levelset_norm%sf
            close (1)
        end if

        ! Outputting Conservative Variables
        do i = 1, sys_size
            write (file_num, '(I0)') i
            file_loc = trim(t_step_dir)//'/q_cons_vf'//trim(file_num) &
                       //'.dat'
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS=status)
            write (1) q_cons_vf(i)%sf(0:m, 0:n, 0:p)
            close (1)
        end do

        !Outputting pb and mv for non-polytropic qbmm
        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do r = 1, nnode
                    write (file_num, '(I0)') r + (i - 1)*nnode + sys_size
                    file_loc = trim(t_step_dir)//'/pb'//trim(file_num) &
                               //'.dat'
                    open (1, FILE=trim(file_loc), FORM='unformatted', &
                          STATUS=status)
                    write (1) pb%sf(:, :, :, r, i)
                    close (1)
                end do
            end do

            do i = 1, nb
                do r = 1, nnode
                    write (file_num, '(I0)') r + (i - 1)*nnode + sys_size
                    file_loc = trim(t_step_dir)//'/mv'//trim(file_num) &
                               //'.dat'
                    open (1, FILE=trim(file_loc), FORM='unformatted', &
                          STATUS=status)
                    write (1) mv%sf(:, :, :, r, i)
                    close (1)
                end do
            end do
        end if

        gamma = fluid_pp(1)%gamma
        lit_gamma = 1._wp/fluid_pp(1)%gamma + 1._wp
        pi_inf = fluid_pp(1)%pi_inf
        qv = fluid_pp(1)%qv

        if (precision == 1) then
            FMT = "(2F30.3)"
        else
            FMT = "(2F40.14)"
        end if

        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir)//'/D'
        file_loc = trim(t_step_dir)//'/.'

        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (.not. file_exist) call s_create_directory(trim(t_step_dir))

        if (cfl_dt) t_step = n_start

        !1D
        if (n == 0 .and. p == 0) then
            if (model_eqns == 2) then
                do i = 1, sys_size
                    write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/prim.', i, '.', proc_rank, '.', t_step, '.dat'

                    open (2, FILE=trim(file_loc))
                    do j = 0, m

                        if (chemistry) then
                            do c = 1, num_species
                                rhoYks(c) = q_cons_vf(chemxb + c - 1)%sf(j, 0, 0)
                            end do
                        end if

                        call s_convert_to_mixture_variables(q_cons_vf, j, 0, 0, rho, gamma, pi_inf, qv)

                        lit_gamma = 1._wp/gamma + 1._wp

                        if ((i >= chemxb) .and. (i <= chemxe)) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)/rho
                        else if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                                 .or. &
                                 ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
                                 .or. &
                                 ((i >= chemxb) .and. (i <= chemxe)) &
                                 ) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        else if (i == mom_idx%beg) then !u
                            write (2, FMT) x_cb(j), q_cons_vf(mom_idx%beg)%sf(j, 0, 0)/rho
                        else if (i == stress_idx%beg) then !tau_e
                            write (2, FMT) x_cb(j), q_cons_vf(stress_idx%beg)%sf(j, 0, 0)/rho
                        else if (i == E_idx) then !p
                            if (mhd) then
                                pres_mag = 0.5_wp*(Bx0**2 + q_cons_vf(B_idx%beg)%sf(j, 0, 0)**2 + q_cons_vf(B_idx%beg + 1)%sf(j, 0, 0)**2)
                            end if

                            call s_compute_pressure( &
                                q_cons_vf(E_idx)%sf(j, 0, 0), &
                                q_cons_vf(alf_idx)%sf(j, 0, 0), &
                                0.5_wp*(q_cons_vf(mom_idx%beg)%sf(j, 0, 0)**2._wp)/rho, &
                                pi_inf, gamma, rho, qv, rhoYks, pres, T, pres_mag=pres_mag)
                            write (2, FMT) x_cb(j), pres
                        else if (mhd) then
                            if (i == mom_idx%beg + 1) then ! v
                                write (2, FMT) x_cb(j), q_cons_vf(mom_idx%beg + 1)%sf(j, 0, 0)/rho
                            else if (i == mom_idx%beg + 2) then ! w
                                write (2, FMT) x_cb(j), q_cons_vf(mom_idx%beg + 2)%sf(j, 0, 0)/rho
                            else if (i == B_idx%beg) then ! By
                                write (2, FMT) x_cb(j), q_cons_vf(B_idx%beg)%sf(j, 0, 0)/rho
                            else if (i == B_idx%beg + 1) then ! Bz
                                write (2, FMT) x_cb(j), q_cons_vf(B_idx%beg + 1)%sf(j, 0, 0)/rho
                            end if
                        else if ((i >= bub_idx%beg) .and. (i <= bub_idx%end) .and. bubbles_euler) then

                            if (qbmm) then
                                nbub = q_cons_vf(bubxb)%sf(j, 0, 0)
                            else
                                if (adv_n) then
                                    nbub = q_cons_vf(n_idx)%sf(j, 0, 0)
                                else
                                    do k = 1, nb
                                        nRtmp(k) = q_cons_vf(bub_idx%rs(k))%sf(j, 0, 0)
                                    end do

                                    call s_comp_n_from_cons(q_cons_vf(alf_idx)%sf(j, 0, 0), nRtmp, nbub, weight)
                                end if
                            end if
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)/nbub
                        else if (i == n_idx .and. adv_n .and. bubbles_euler) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        else if (i == damage_idx) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        end if
                    end do
                    close (2)
                end do
            end if

            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'

                open (2, FILE=trim(file_loc))
                do j = 0, m
                    write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                end do
                close (2)
            end do

            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/pres.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_loc))
                        do j = 0, m
                            write (2, FMT) x_cb(j), pb%sf(j, 0, 0, r, i)
                        end do
                        close (2)
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/mv.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_loc))
                        do j = 0, m
                            write (2, FMT) x_cb(j), mv%sf(j, 0, 0, r, i)
                        end do
                        close (2)
                    end do
                end do
            end if
        end if

        if (precision == 1) then
            FMT = "(3F30.7)"
        else
            FMT = "(3F40.14)"
        end if

        ! 2D
        if ((n > 0) .and. (p == 0)) then
            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_loc))
                do j = 0, m
                    do k = 0, n
                        write (2, FMT) x_cb(j), y_cb(k), q_cons_vf(i)%sf(j, k, 0)
                    end do
                    write (2, *)
                end do
                close (2)
            end do

            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/pres.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_loc))
                        do j = 0, m
                            do k = 0, n
                                write (2, FMT) x_cb(j), y_cb(k), pb%sf(j, k, 0, r, i)
                            end do
                        end do
                        close (2)
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/mv.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_loc))
                        do j = 0, m
                            do k = 0, n
                                write (2, FMT) x_cb(j), y_cb(k), mv%sf(j, k, 0, r, i)
                            end do
                        end do
                        close (2)
                    end do
                end do
            end if
        end if

        if (precision == 1) then
            FMT = "(4F30.7)"
        else
            FMT = "(4F40.14)"
        end if

        ! 3D
        if (p > 0) then
            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_loc))
                do j = 0, m
                    do k = 0, n
                        do l = 0, p
                            write (2, FMT) x_cb(j), y_cb(k), z_cb(l), q_cons_vf(i)%sf(j, k, l)
                        end do
                        write (2, *)
                    end do
                    write (2, *)
                end do
                close (2)
            end do

            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/pres.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_loc))
                        do j = 0, m
                            do k = 0, n
                                do l = 0, p
                                    write (2, FMT) x_cb(j), y_cb(k), z_cb(l), pb%sf(j, k, l, r, i)
                                end do
                            end do
                        end do
                        close (2)
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/mv.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_loc))
                        do j = 0, m
                            do k = 0, n
                                do l = 0, p
                                    write (2, FMT) x_cb(j), y_cb(k), z_cb(l), mv%sf(j, k, l, r, i)
                                end do
                            end do
                        end do
                        close (2)
                    end do
                end do
            end if
        end if

        if (ib) then

            ! Write IB Markers
            write (file_loc, '(A,I2.2,A)') trim(t_step_dir)//'/ib_markers.', proc_rank, '.dat'
            open (2, FILE=trim(file_loc))
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        if (p > 0) then
                            write (2, FMT) x_cc(j), y_cc(k), z_cc(l), real(ib_markers%sf(j, k, l))
                        else
                            write (2, FMT) x_cc(j), y_cc(k), real(ib_markers%sf(j, k, l))
                        end if
                    end do
                end do
            end do

            close (2)

            do i = 1, num_ibs
                if (patch_ib(i)%geometry == 4) then

                    write (file_loc, '(A,I2.2,A)') trim(t_step_dir)//'/airfoil_u.', proc_rank, '.dat'
                    open (2, FILE=trim(file_loc))
                    do j = 1, Np
                        write (2, FMT) airfoil_grid_u(j)%x, airfoil_grid_u(j)%y
                    end do
                    close (2)

                    write (file_loc, '(A,I2.2,A)') trim(t_step_dir)//'/airfoil_l.', proc_rank, '.dat'
                    open (2, FILE=trim(file_loc))
                    do j = 1, Np
                        write (2, FMT) airfoil_grid_l(j)%x, airfoil_grid_l(j)%y
                    end do
                    close (2)

                    print *, "Np", Np
                end if
            end do
        end if

    end subroutine s_write_serial_data_files

    !> Writes grid and initial condition data files in parallel to the "0"
        !!  time-step directory in the local processor rank folder
        !! @param q_cons_vf Conservative variables
        !! @param ib_markers track if a cell is within the immersed boundary
        !! @param levelset closest distance from every cell to the IB
        !! @param levelset_norm normalized vector from every cell to the closest point to the IB
    impure subroutine s_write_parallel_data_files(q_cons_vf, q_prim_vf, bc_type, ib_markers, levelset, levelset_norm)

        ! Conservative variables
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_cons_vf, q_prim_vf

        type(integer_field), &
            dimension(1:num_dims, -1:1), &
            intent(in) :: bc_type

        ! IB markers
        type(integer_field), &
            intent(in), optional :: ib_markers

        ! Levelset
        type(levelset_field), &
            intent(IN), optional :: levelset

        ! Levelset Norm
        type(levelset_norm_field), &
            intent(IN), optional :: levelset_norm

#ifdef MFC_MPI

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist, dir_check

        ! Generic loop iterators
        integer :: i, j, k, l
        real(wp) :: loc_violations, glb_violations

        ! Downsample variables
        integer :: m_ds, n_ds, p_ds
        integer :: m_glb_ds, n_glb_ds, p_glb_ds
        integer :: m_glb_save, n_glb_save, p_glb_save ! Size of array being saved

        if (down_sample) then
            if ((mod(m + 1, 3) > 0) .or. (mod(n + 1, 3) > 0) .or. (mod(p + 1, 3) > 0)) then
                loc_violations = 1._wp
            end if
            call s_mpi_allreduce_sum(loc_violations, glb_violations)
            if (proc_rank == 0 .and. nint(glb_violations) > 0) then
                print *, "WARNING: Attempting to downsample data but there are"// &
                    "processors with local problem sizes that are not divisible by 3."
            end if
            call s_populate_variables_buffers(bc_type, q_cons_vf)
            call s_downsample_data(q_cons_vf, q_cons_temp, &
                                   m_ds, n_ds, p_ds, m_glb_ds, n_glb_ds, p_glb_ds)
        end if

        if (file_per_process) then
            if (proc_rank == 0) then
                file_loc = trim(case_dir)//'/restart_data/lustre_0'
                call my_inquire(file_loc, dir_check)
                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(file_loc))
                end if
                call s_create_directory(trim(file_loc))
            end if
            call s_mpi_barrier()
            call DelayFileAccess(proc_rank)

            ! Initialize MPI data I/O
            if (down_sample) then
                call s_initialize_mpi_data_ds(q_cons_temp)
            else
                if (ib) then
                    call s_initialize_mpi_data(q_cons_vf, ib_markers, &
                                               levelset, levelset_norm)
                else
                    call s_initialize_mpi_data(q_cons_vf)
                end if
            end if

            ! Open the file to write all flow variables
            if (cfl_dt) then
                write (file_loc, '(I0,A,i7.7,A)') n_start, '_', proc_rank, '.dat'
            else
                write (file_loc, '(I0,A,i7.7,A)') t_step_start, '_', proc_rank, '.dat'
            end if
            file_loc = trim(restart_dir)//'/lustre_0'//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            if (file_exist) call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            if (down_sample) then
                ! Size of local arrays
                data_size = (m_ds + 3)*(n_ds + 3)*(p_ds + 3)
                m_glb_save = m_glb_ds + 3
                n_glb_save = n_glb_ds + 3
                p_glb_save = p_glb_ds + 3
            else
                ! Size of local arrays
                data_size = (m + 1)*(n + 1)*(p + 1)
                m_glb_save = m_glb + 1
                n_glb_save = n_glb + 1
                p_glb_save = p_glb + 1
            end if

            ! Resize some integers so MPI can write even the biggest files
            m_MOK = int(m_glb_save, MPI_OFFSET_KIND)
            n_MOK = int(n_glb_save, MPI_OFFSET_KIND)
            p_MOK = int(p_glb_save, MPI_OFFSET_KIND)
            WP_MOK = int(8._wp, MPI_OFFSET_KIND)
            MOK = int(1._wp, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            ! Write the data for each variable
            if (bubbles_euler) then
                do i = 1, sys_size! adv_idx%end
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            mpi_p, status, ierr)
                end do
                !Additional variables pb and mv for non-polytropic qbmm
                if (qbmm .and. .not. polytropic) then
                    do i = sys_size + 1, sys_size + 2*nb*nnode
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                                mpi_p, status, ierr)
                    end do
                end if
            else
                if (down_sample) then
                    do i = 1, sys_size !TODO: check if this is right
                        !            do i = 1, adv_idx%end
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, q_cons_temp(i)%sf, data_size, &
                                                mpi_p, status, ierr)
                    end do
                else
                    do i = 1, sys_size !TODO: check if this is right
                        !            do i = 1, adv_idx%end
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                                mpi_p, status, ierr)
                    end do
                end if
            end if

            call MPI_FILE_CLOSE(ifile, ierr)

        else
            ! Initialize MPI data I/O
            if (ib) then
                call s_initialize_mpi_data(q_cons_vf, ib_markers, &
                                           levelset, levelset_norm)
            else
                call s_initialize_mpi_data(q_cons_vf)
            end if

            ! Open the file to write all flow variables
            if (cfl_dt) then
                write (file_loc, '(I0,A)') n_start, '.dat'
            else
                write (file_loc, '(I0,A)') t_step_start, '.dat'
            end if
            file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Size of local arrays
            data_size = (m + 1)*(n + 1)*(p + 1)

            ! Resize some integers so MPI can write even the biggest files
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8._wp, MPI_OFFSET_KIND)
            MOK = int(1._wp, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            ! Write the data for each variable
            if (bubbles_euler) then
                do i = 1, sys_size! adv_idx%end
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            mpi_p, status, ierr)
                end do
                !Additional variables pb and mv for non-polytropic qbmm
                if (qbmm .and. .not. polytropic) then
                    do i = sys_size + 1, sys_size + 2*nb*nnode
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                                mpi_p, status, ierr)
                    end do
                end if
            else
                do i = 1, sys_size !TODO: check if this is right
                    !            do i = 1, adv_idx%end
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            mpi_p, status, ierr)
                end do

            end if

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        ! IB Markers
        if (ib) then

            write (file_loc, '(A)') 'ib.dat'
            file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Initial displacement to skip at beginning of file
            disp = 0

            call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_IO_IB_DATA%view, &
                                   'native', mpi_info_int, ierr)
            call MPI_FILE_WRITE_ALL(ifile, MPI_IO_IB_DATA%var%sf, data_size, &
                                    MPI_INTEGER, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)

            ! Levelset
            write (file_loc, '(A)') 'levelset.dat'
            file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Initial displacement to skip at beginning of file
            disp = 0

            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_levelset_DATA%view, &
                                   'native', mpi_info_int, ierr)
            call MPI_FILE_WRITE_ALL(ifile, MPI_IO_levelset_DATA%var%sf, data_size*num_ibs, &
                                    mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)

            ! Levelset Norm
            write (file_loc, '(A)') 'levelset_norm.dat'
            file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Initial displacement to skip at beginning of file
            disp = 0

            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_levelsetnorm_DATA%view, &
                                   'native', mpi_info_int, ierr)
            call MPI_FILE_WRITE_ALL(ifile, MPI_IO_levelsetnorm_DATA%var%sf, data_size*num_ibs*3, &
                                    mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        if (ib) then

            do i = 1, num_ibs

                if (patch_ib(i)%geometry == 4) then

                    write (file_loc, '(A)') 'airfoil_l.dat'
                    file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)
                    if (file_exist .and. proc_rank == 0) then
                        call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
                    end if
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                                       mpi_info_int, ifile, ierr)

                    ! Initial displacement to skip at beginning of file
                    disp = 0

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_airfoil_IB_DATA%view(1), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_airfoil_IB_DATA%var(1:Np), 3*Np, &
                                            mpi_p, status, ierr)

                    call MPI_FILE_CLOSE(ifile, ierr)

                    write (file_loc, '(A)') 'airfoil_u.dat'
                    file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)
                    if (file_exist .and. proc_rank == 0) then
                        call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
                    end if
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                                       mpi_info_int, ifile, ierr)

                    ! Initial displacement to skip at beginning of file
                    disp = 0

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_airfoil_IB_DATA%view(2), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_airfoil_IB_DATA%var(Np + 1:2*Np), 3*Np, &
                                            mpi_p, status, ierr)

                    call MPI_FILE_CLOSE(ifile, ierr)
                end if
            end do

        end if
#endif

        if (bc_io) then
            if (igr) then
                call s_write_parallel_boundary_condition_files(q_cons_vf, bc_type)
            else
                call s_write_parallel_boundary_condition_files(q_prim_vf, bc_type)
            end if
        end if

    end subroutine s_write_parallel_data_files

    !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
    impure subroutine s_initialize_data_output_module
        ! Generic string used to store the address of a particular file
        character(LEN=len_trim(case_dir) + 2*name_len) :: file_loc
        character(len=15) :: temp
        character(LEN=1), dimension(3), parameter :: coord = (/'x', 'y', 'z'/)

        ! Generic logical used to check the existence of directories
        logical :: dir_check
        integer :: i

        integer :: m_ds, n_ds, p_ds !< down sample dimensions

        if (parallel_io .neqv. .true.) then
            ! Setting the address of the time-step directory
            write (t_step_dir, '(A,I0,A)') '/p_all/p', proc_rank, '/0'
            t_step_dir = trim(case_dir)//trim(t_step_dir)

            ! Checking the existence of the time-step directory, removing it, if
            ! it exists, and creating a new copy. Note that if preexisting grid
            ! and/or initial condition data are to be read in from the very same
            ! location, then the above described steps are not executed here but
            ! rather in the module m_start_up.f90.
            if (old_grid .neqv. .true.) then

                file_loc = trim(t_step_dir)//'/'

                call my_inquire(file_loc, dir_check)

                if (dir_check) call s_delete_directory(trim(t_step_dir))

                call s_create_directory(trim(t_step_dir))

            end if

            s_write_data_files => s_write_serial_data_files
        else
            write (restart_dir, '(A)') '/restart_data'
            restart_dir = trim(case_dir)//trim(restart_dir)

            if ((old_grid .neqv. .true.) .and. (proc_rank == 0)) then

                file_loc = trim(restart_dir)//'/'
                call my_inquire(file_loc, dir_check)

                if (dir_check) call s_delete_directory(trim(restart_dir))
                call s_create_directory(trim(restart_dir))
            end if

            call s_mpi_barrier()

            s_write_data_files => s_write_parallel_data_files

        end if

        open (1, FILE='indices.dat', STATUS='unknown')

        write (1, '(A)') "Warning: The creation of file is currently experimental."
        write (1, '(A)') "This file may contain errors and not support all features."

        write (1, '(A3,A20,A20)') "#", "Conservative", "Primitive"
        write (1, '(A)') "    "
        do i = contxb, contxe
            write (temp, '(I0)') i - contxb + 1
            write (1, '(I3,A20,A20)') i, "\alpha_{"//trim(temp)//"} \rho_{"//trim(temp)//"}", "\alpha_{"//trim(temp)//"} \rho"
        end do
        do i = momxb, momxe
            write (1, '(I3,A20,A20)') i, "\rho u_"//coord(i - momxb + 1), "u_"//coord(i - momxb + 1)
        end do
        do i = E_idx, E_idx
            write (1, '(I3,A20,A20)') i, "\rho U", "p"
        end do
        do i = advxb, advxe
            write (temp, '(I0)') i - contxb + 1
            write (1, '(I3,A20,A20)') i, "\alpha_{"//trim(temp)//"}", "\alpha_{"//trim(temp)//"}"
        end do
        if (chemistry) then
            do i = 1, num_species
                write (1, '(I3,A20,A20)') chemxb + i - 1, "Y_{"//trim(species_names(i))//"} \rho", "Y_{"//trim(species_names(i))//"}"
            end do
        end if

        write (1, '(A)') ""
        if (momxb /= 0) write (1, '("[",I2,",",I2,"]",A)') momxb, momxe, " Momentum"
        if (E_idx /= 0) write (1, '("[",I2,",",I2,"]",A)') E_idx, E_idx, " Energy/Pressure"
        if (advxb /= 0) write (1, '("[",I2,",",I2,"]",A)') advxb, advxe, " Advection"
        if (contxb /= 0) write (1, '("[",I2,",",I2,"]",A)') contxb, contxe, " Continuity"
        if (bubxb /= 0) write (1, '("[",I2,",",I2,"]",A)') bubxb, bubxe, " Bubbles_euler"
        if (strxb /= 0) write (1, '("[",I2,",",I2,"]",A)') strxb, strxe, " Stress"
        if (intxb /= 0) write (1, '("[",I2,",",I2,"]",A)') intxb, intxe, " Internal Energies"
        if (chemxb /= 0) write (1, '("[",I2,",",I2,"]",A)') chemxb, chemxe, " Chemistry"

        close (1)

        if (down_sample) then
            m_ds = int((m + 1)/3) - 1
            n_ds = int((n + 1)/3) - 1
            p_ds = int((p + 1)/3) - 1

            allocate (q_cons_temp(1:sys_size))
            do i = 1, sys_size
                allocate (q_cons_temp(i)%sf(-1:m_ds + 1, -1:n_ds + 1, -1:p_ds + 1))
            end do
        end if

    end subroutine s_initialize_data_output_module

    !> Resets s_write_data_files pointer
    impure subroutine s_finalize_data_output_module

        integer :: i

        s_write_data_files => null()

        if (down_sample) then
            do i = 1, sys_size
                deallocate (q_cons_temp(i)%sf)
            end do
            deallocate (q_cons_temp)
        end if

    end subroutine s_finalize_data_output_module

end module m_data_output
