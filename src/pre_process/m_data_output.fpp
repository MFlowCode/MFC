!>
!! @file
!! @brief Contains module m_data_output

!> @brief Writes grid and initial condition data to serial or parallel output files
module m_data_output

    use m_derived_types
    use m_global_parameters
    use m_helper
    use m_mpi_proxy
#ifdef MFC_MPI
    use mpi
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

    private
    public :: s_write_serial_data_files, s_write_parallel_data_files, s_write_data_files, s_initialize_data_output_module, &
        & s_finalize_data_output_module

    type(scalar_field), allocatable, dimension(:) :: q_cons_temp

    abstract interface

        !> Interface for the conservative data
        impure subroutine s_write_abstract_data_files(q_cons_vf, q_prim_vf, bc_type, q_T_sf)

            import :: scalar_field, integer_field, sys_size, m, n, p, pres_field, num_dims

            type(scalar_field), dimension(sys_size), intent(inout)      :: q_cons_vf, q_prim_vf
            type(integer_field), dimension(1:num_dims,-1:1), intent(in) :: bc_type
            type(scalar_field), intent(inout), optional                 :: q_T_sf

        end subroutine s_write_abstract_data_files
    end interface

    !> Time-step folder into which grid and initial condition data will be placed
    character(LEN=path_len + 2*name_len), private :: t_step_dir
    character(LEN=path_len + 2*name_len), public  :: restart_dir  !< Restart data folder
    procedure(s_write_abstract_data_files), pointer :: s_write_data_files => null()

contains

    !> Writes grid and initial condition data files to the "0" time-step directory in the local processor rank folder
    impure subroutine s_write_serial_data_files(q_cons_vf, q_prim_vf, bc_type, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(inout)      :: q_cons_vf, q_prim_vf
        type(integer_field), dimension(1:num_dims,-1:1), intent(in) :: bc_type
        logical                                                     :: file_exist
        character(LEN=15)                                           :: FMT
        character(LEN=3)                                            :: status
        character(LEN=int(floor(log10(real(sys_size, wp)))) + 1)    :: file_num
        character(LEN=len_trim(t_step_dir) + name_len)              :: file_loc
        integer                                                     :: i, j, k, l, r, c
        integer                                                     :: t_step
        real(wp), dimension(nb)                                     :: nRtmp
        real(wp)                                                    :: nbub
        real(wp)                                                    :: gamma, lit_gamma, pi_inf, qv
        real(wp)                                                    :: rho
        real(wp)                                                    :: pres, T
        real(wp)                                                    :: rhoYks(1:num_species)
        real(wp)                                                    :: pres_mag
        type(scalar_field), intent(inout), optional                 :: q_T_sf

        pres_mag = 0._wp

        T = dflt_T_guess

        t_step = 0

        if (old_grid) then
            status = 'old'
        else
            status = 'new'
        end if

        if (bc_io) then
            if (igr) then
                call s_write_serial_boundary_condition_files(q_cons_vf, bc_type, t_step_dir, old_grid)
            else
                call s_write_serial_boundary_condition_files(q_prim_vf, bc_type, t_step_dir, old_grid, q_T_sf)
            end if
        end if

        file_loc = trim(t_step_dir) // '/x_cb.dat'
        open (1, FILE=trim(file_loc), form='unformatted', STATUS=status)
        write (1) x_cb(-1:m)
        close (1)

        if (n > 0) then
            file_loc = trim(t_step_dir) // '/y_cb.dat'
            open (1, FILE=trim(file_loc), form='unformatted', STATUS=status)
            write (1) y_cb(-1:n)
            close (1)

            if (p > 0) then
                file_loc = trim(t_step_dir) // '/z_cb.dat'
                open (1, FILE=trim(file_loc), form='unformatted', STATUS=status)
                write (1) z_cb(-1:p)
                close (1)
            end if
        end if

        do i = 1, sys_size
            write (file_num, '(I0)') i
            file_loc = trim(t_step_dir) // '/q_cons_vf' // trim(file_num) // '.dat'
            open (1, FILE=trim(file_loc), form='unformatted', STATUS=status)
            write (1) q_cons_vf(i)%sf(0:m,0:n,0:p)
            close (1)
        end do

        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do r = 1, nnode
                    write (file_num, '(I0)') r + (i - 1)*nnode + sys_size
                    file_loc = trim(t_step_dir) // '/pb' // trim(file_num) // '.dat'
                    open (1, FILE=trim(file_loc), form='unformatted', STATUS=status)
                    write (1) pb%sf(:,:,:,r, i)
                    close (1)
                end do
            end do

            do i = 1, nb
                do r = 1, nnode
                    write (file_num, '(I0)') r + (i - 1)*nnode + sys_size
                    file_loc = trim(t_step_dir) // '/mv' // trim(file_num) // '.dat'
                    open (1, FILE=trim(file_loc), form='unformatted', STATUS=status)
                    write (1) mv%sf(:,:,:,r, i)
                    close (1)
                end do
            end do
        end if

        gamma = gammas(1)
        lit_gamma = gs_min(1)
        pi_inf = pi_infs(1)
        qv = qvs(1)

        if (precision == 1) then
            FMT = "(2F30.3)"
        else
            FMT = "(2F40.14)"
        end if

        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir) // '/D'
        file_loc = trim(t_step_dir) // '/.'

        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (.not. file_exist) call s_create_directory(trim(t_step_dir))

        if (cfl_dt) t_step = n_start

        if (n == 0 .and. p == 0) then
            if (model_eqns == 2) then
                do i = 1, sys_size
                    write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/prim.', i, '.', proc_rank, '.', t_step, '.dat'

                    open (2, FILE=trim(file_loc))
                    do j = 0, m
                        if (chemistry) then
                            do c = 1, num_species
                                rhoYks(c) = q_cons_vf(eqn_idx%species%beg + c - 1)%sf(j, 0, 0)
                            end do
                        end if

                        call s_convert_to_mixture_variables(q_cons_vf, j, 0, 0, rho, gamma, pi_inf, qv)

                        lit_gamma = 1._wp/gamma + 1._wp

                        if ((i >= eqn_idx%species%beg) .and. (i <= eqn_idx%species%end)) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)/rho
                        else if (((i >= eqn_idx%cont%beg) .and. (i <= eqn_idx%cont%end)) .or. ((i >= eqn_idx%adv%beg) &
                                 & .and. (i <= eqn_idx%adv%end)) .or. ((i >= eqn_idx%species%beg) .and. (i <= eqn_idx%species%end) &
                                 & )) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        else if (i == eqn_idx%mom%beg) then  ! u
                            write (2, FMT) x_cb(j), q_cons_vf(eqn_idx%mom%beg)%sf(j, 0, 0)/rho
                        else if (i == eqn_idx%stress%beg) then  ! tau_e
                            write (2, FMT) x_cb(j), q_cons_vf(eqn_idx%stress%beg)%sf(j, 0, 0)/rho
                        else if (i == eqn_idx%E) then  ! p
                            if (mhd) then
                                pres_mag = 0.5_wp*(Bx0**2 + q_cons_vf(eqn_idx%B%beg)%sf(j, 0, &
                                                   & 0)**2 + q_cons_vf(eqn_idx%B%beg + 1)%sf(j, 0, 0)**2)
                            end if

                            call s_compute_pressure(q_cons_vf(eqn_idx%E)%sf(j, 0, 0), q_cons_vf(eqn_idx%alf)%sf(j, 0, 0), &
                                                    & 0.5_wp*(q_cons_vf(eqn_idx%mom%beg)%sf(j, 0, 0)**2._wp)/rho, pi_inf, gamma, &
                                                    & rho, qv, rhoYks, pres, T, pres_mag=pres_mag)
                            write (2, FMT) x_cb(j), pres
                        else if (mhd) then
                            if (i == eqn_idx%mom%beg + 1) then  ! v
                                write (2, FMT) x_cb(j), q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, 0, 0)/rho
                            else if (i == eqn_idx%mom%beg + 2) then  ! w
                                write (2, FMT) x_cb(j), q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, 0, 0)/rho
                            else if (i == eqn_idx%B%beg) then  ! By
                                write (2, FMT) x_cb(j), q_cons_vf(eqn_idx%B%beg)%sf(j, 0, 0)/rho
                            else if (i == eqn_idx%B%beg + 1) then  ! Bz
                                write (2, FMT) x_cb(j), q_cons_vf(eqn_idx%B%beg + 1)%sf(j, 0, 0)/rho
                            end if
                        else if ((i >= eqn_idx%bub%beg) .and. (i <= eqn_idx%bub%end) .and. bubbles_euler) then
                            if (qbmm) then
                                nbub = q_cons_vf(eqn_idx%bub%beg)%sf(j, 0, 0)
                            else
                                if (adv_n) then
                                    nbub = q_cons_vf(eqn_idx%n)%sf(j, 0, 0)
                                else
                                    do k = 1, nb
                                        nRtmp(k) = q_cons_vf(qbmm_idx%rs(k))%sf(j, 0, 0)
                                    end do

                                    call s_comp_n_from_cons(real(q_cons_vf(eqn_idx%alf)%sf(j, 0, 0), kind=wp), nRtmp, nbub, weight)
                                end if
                            end if
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)/nbub
                        else if (i == eqn_idx%n .and. adv_n .and. bubbles_euler) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        else if (i == eqn_idx%damage) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        end if
                    end do
                    close (2)
                end do
            end if

            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step, '.dat'

                open (2, FILE=trim(file_loc))
                do j = 0, m
                    write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                end do
                close (2)
            end do

            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/pres.', i, '.', r, '.', proc_rank, &
                               & '.', t_step, '.dat'

                        open (2, FILE=trim(file_loc))
                        do j = 0, m
                            write (2, FMT) x_cb(j), pb%sf(j, 0, 0, r, i)
                        end do
                        close (2)
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/mv.', i, '.', r, '.', proc_rank, &
                               & '.', t_step, '.dat'

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

        if ((n > 0) .and. (p == 0)) then
            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step, '.dat'
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
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/pres.', i, '.', r, '.', proc_rank, &
                               & '.', t_step, '.dat'

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
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/mv.', i, '.', r, '.', proc_rank, &
                               & '.', t_step, '.dat'

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

        if (p > 0) then
            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step, '.dat'
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
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/pres.', i, '.', r, '.', proc_rank, &
                               & '.', t_step, '.dat'

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
                        write (file_loc, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/mv.', i, '.', r, '.', proc_rank, &
                               & '.', t_step, '.dat'

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

    end subroutine s_write_serial_data_files

    !> Writes grid and initial condition data files in parallel to the "0" time-step directory in the local processor rank folder
    impure subroutine s_write_parallel_data_files(q_cons_vf, q_prim_vf, bc_type, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(inout)      :: q_cons_vf, q_prim_vf
        type(integer_field), dimension(1:num_dims,-1:1), intent(in) :: bc_type
        type(scalar_field), optional, intent(inout)                 :: q_T_sf

#ifdef MFC_MPI
        integer                              :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE)  :: status
        integer(KIND=MPI_OFFSET_KIND)        :: disp
        integer(KIND=MPI_OFFSET_KIND)        :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: MOK
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist, dir_check
        integer                              :: i, j, k, l
        real(wp)                             :: loc_violations, glb_violations
        integer                              :: m_ds, n_ds, p_ds
        integer                              :: m_glb_ds, n_glb_ds, p_glb_ds
        integer                              :: m_glb_save, n_glb_save, p_glb_save  !< Size of array being saved

        loc_violations = 0._wp

        if (down_sample) then
            if ((mod(m + 1, 3) > 0) .or. (mod(n + 1, 3) > 0) .or. (mod(p + 1, 3) > 0)) then
                loc_violations = 1._wp
            end if
            call s_mpi_allreduce_sum(loc_violations, glb_violations)
            if (proc_rank == 0 .and. nint(glb_violations) > 0) then
                print *, &
                    & "WARNING: Attempting to downsample data but there are" &
                    & // "processors with local problem sizes that are not divisible by 3."
            end if
            call s_populate_variables_buffers(bc_type, q_cons_vf)
            call s_downsample_data(q_cons_vf, q_cons_temp, m_ds, n_ds, p_ds, m_glb_ds, n_glb_ds, p_glb_ds)
        end if

        if (file_per_process) then
            if (proc_rank == 0) then
                file_loc = trim(case_dir) // '/restart_data/lustre_0'
                call my_inquire(file_loc, dir_check)
                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(file_loc))
                end if
                call s_create_directory(trim(file_loc))
            end if
            call s_mpi_barrier()
            call DelayFileAccess(proc_rank)

            if (down_sample) then
                call s_initialize_mpi_data_ds(q_cons_temp)
            else
                call s_initialize_mpi_data(q_cons_vf)
            end if

            if (cfl_dt) then
                write (file_loc, '(I0,A,i7.7,A)') n_start, '_', proc_rank, '.dat'
            else
                write (file_loc, '(I0,A,i7.7,A)') t_step_start, '_', proc_rank, '.dat'
            end if
            file_loc = trim(restart_dir) // '/lustre_0' // trim(mpiiofs) // trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            if (file_exist) call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)

            if (down_sample) then
                data_size = (m_ds + 3)*(n_ds + 3)*(p_ds + 3)
                m_glb_save = m_glb_ds + 3
                n_glb_save = n_glb_ds + 3
                p_glb_save = p_glb_ds + 3
            else
                data_size = (m + 1)*(n + 1)*(p + 1)
                m_glb_save = m_glb + 1
                n_glb_save = n_glb + 1
                p_glb_save = p_glb + 1
            end if

            ! Resize some integers so MPI can write even the biggest files
            m_MOK = int(m_glb_save, MPI_OFFSET_KIND)
            n_MOK = int(n_glb_save, MPI_OFFSET_KIND)
            p_MOK = int(p_glb_save, MPI_OFFSET_KIND)
            WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
            MOK = int(1._wp, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            if (bubbles_euler) then
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                end do
                if (qbmm .and. .not. polytropic) then
                    do i = sys_size + 1, sys_size + 2*nb*nnode
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                end if
            else
                if (down_sample) then
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, q_cons_temp(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                else
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                end if
            end if

            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_initialize_mpi_data(q_cons_vf)

            if (cfl_dt) then
                write (file_loc, '(I0,A)') n_start, '.dat'
            else
                write (file_loc, '(I0,A)') t_step_start, '.dat'
            end if
            file_loc = trim(restart_dir) // trim(mpiiofs) // trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)

            data_size = (m + 1)*(n + 1)*(p + 1)

            ! Resize some integers so MPI can write even the biggest files
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
            MOK = int(1._wp, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            if (bubbles_euler) then
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_io_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                end do
                if (qbmm .and. .not. polytropic) then
                    do i = sys_size + 1, sys_size + 2*nb*nnode
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_io_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                end if
            else
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_io_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                end do
            end if

            call MPI_FILE_CLOSE(ifile, ierr)
        end if
#endif

        if (bc_io) then
            if (igr) then
                call s_write_parallel_boundary_condition_files(q_cons_vf, bc_type)
            else
                call s_write_parallel_boundary_condition_files(q_prim_vf, bc_type, q_T_sf)
            end if
        end if

    end subroutine s_write_parallel_data_files

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_data_output_module

        character(LEN=len_trim(case_dir) + 2*name_len) :: file_loc
        character(len=15)                              :: temp
        character(LEN=1), dimension(3), parameter      :: coord = (/'x', 'y', 'z'/)
        logical                                        :: dir_check
        integer                                        :: i, iu
        integer                                        :: m_ds, n_ds, p_ds

        if (parallel_io .neqv. .true.) then
            write (t_step_dir, '(A,I0,A)') '/p_all/p', proc_rank, '/0'
            t_step_dir = trim(case_dir) // trim(t_step_dir)

            if (old_grid .neqv. .true.) then
                file_loc = trim(t_step_dir) // '/'

                call my_inquire(file_loc, dir_check)

                if (dir_check) call s_delete_directory(trim(t_step_dir))

                call s_create_directory(trim(t_step_dir))
            end if

            s_write_data_files => s_write_serial_data_files
        else
            write (restart_dir, '(A)') '/restart_data'
            restart_dir = trim(case_dir) // trim(restart_dir)

            if ((old_grid .neqv. .true.) .and. (proc_rank == 0)) then
                file_loc = trim(restart_dir) // '/'
                call my_inquire(file_loc, dir_check)

                if (dir_check) call s_delete_directory(trim(restart_dir))
                call s_create_directory(trim(restart_dir))
            end if

            call s_mpi_barrier()

            s_write_data_files => s_write_parallel_data_files
        end if

        open (newunit=iu, file='indices.dat', status='unknown')

        write (iu, '(A)') "Warning: The creation of file is currently experimental."
        write (iu, '(A)') "This file may contain errors and not support all features."

        write (iu, '(A3,A20,A20)') "#", "Conservative", "Primitive"
        write (iu, '(A)') "    "
        do i = eqn_idx%cont%beg, eqn_idx%cont%end
            write (temp, '(I0)') i - eqn_idx%cont%beg + 1
            write (iu, '(I3,A20,A20)') i, "\alpha_{" // trim(temp) // "} \rho_{" // trim(temp) // "}", &
                   & "\alpha_{" // trim(temp) // "} \rho"
        end do
        do i = eqn_idx%mom%beg, eqn_idx%mom%end
            write (iu, '(I3,A20,A20)') i, "\rho u_" // coord(i - eqn_idx%mom%beg + 1), "u_" // coord(i - eqn_idx%mom%beg + 1)
        end do
        if (eqn_idx%E /= 0) write (iu, '(I3,A20,A20)') eqn_idx%E, "\rho U", "p"
        do i = eqn_idx%adv%beg, eqn_idx%adv%end
            write (temp, '(I0)') i - eqn_idx%cont%beg + 1
            write (iu, '(I3,A20,A20)') i, "\alpha_{" // trim(temp) // "}", "\alpha_{" // trim(temp) // "}"
        end do
        if (chemistry) then
            do i = 1, num_species
                write (iu, '(I3,A20,A20)') eqn_idx%species%beg + i - 1, "Y_{" // trim(species_names(i)) // "} \rho", &
                       & "Y_{" // trim(species_names(i)) // "}"
            end do
        end if

        write (iu, '(A)') ""
        call write_range(eqn_idx%cont%beg, eqn_idx%cont%end, " Continuity")
        call write_range(eqn_idx%mom%beg, eqn_idx%mom%end, " Momentum")
        call write_range(eqn_idx%E, eqn_idx%E, " Energy/Pressure")
        call write_range(eqn_idx%adv%beg, eqn_idx%adv%end, " Advection")
        call write_range(eqn_idx%bub%beg, eqn_idx%bub%end, " Bubbles")
        call write_range(eqn_idx%stress%beg, eqn_idx%stress%end, " Stress")
        call write_range(eqn_idx%int_en%beg, eqn_idx%int_en%end, " Internal Energies")
        call write_range(eqn_idx%xi%beg, eqn_idx%xi%end, " Reference Map")
        call write_range(eqn_idx%B%beg, eqn_idx%B%end, " Magnetic Field")
        call write_range(eqn_idx%c, eqn_idx%c, " Color Function")
        call write_range(eqn_idx%species%beg, eqn_idx%species%end, " Chemistry")

        close (iu)

        if (down_sample) then
            m_ds = int((m + 1)/3) - 1
            n_ds = int((n + 1)/3) - 1
            p_ds = int((p + 1)/3) - 1

            allocate (q_cons_temp(1:sys_size))
            do i = 1, sys_size
                allocate (q_cons_temp(i)%sf(-1:m_ds + 1,-1:n_ds + 1,-1:p_ds + 1))
            end do
        end if

    contains

        subroutine write_range(beg, end, label)

            integer, intent(in)      :: beg, end
            character(*), intent(in) :: label

            if (beg /= 0) write (iu, '("[",I0,",",I0,"]",A)') beg, end, label

        end subroutine write_range

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
