!>
!! @file
!! @brief Contains module m_euler_lagrange

#:include 'macros.fpp'

!> @brief Utilities shared by the Euler-Lagrange bubble and solid-particle solvers: cell location, domain and cell-volume helpers,
!! the void-fraction history file, and the MPI-IO restart files of the Lagrangian state.
module m_euler_lagrange

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_compile_specific

    implicit none

    private
    public :: LAG_EVOL_ID, LAG_VOID_ID, s_locate_cell, particle_in_domain_physical, s_get_char_vol, s_open_void_evol, &
        & s_write_void_evol, s_close_void_evol, s_open_lag_evol, s_close_lag_evol, s_write_lag_restart, s_read_lag_restart, &
        & s_get_lag_restart_point, s_set_lag_comm_coords, s_create_D_dir, s_read_lag_input, s_count_lag_glb

    integer, parameter :: LAG_EVOL_ID = 11  !< File id for the per-rank Lagrangian evolution files
    integer, parameter :: LAG_VOID_ID = 13  !< File id for D/voidfraction.dat
    !> Lagrangian volume fraction above which a cell counts in the cloud statistics
    real(wp), parameter :: void_stat_min = 5.e-11_wp

contains

    !> Save index and time the run starts from: zero on a fresh start, the restart point otherwise.
    subroutine s_get_lag_restart_point(save_count, qtime)

        integer, intent(out)  :: save_count
        real(wp), intent(out) :: qtime

        if (cfl_dt) then
            save_count = n_start
            qtime = n_start*t_save
        else
            save_count = t_step_start
            qtime = t_step_start*dt
        end if

    end subroutine s_get_lag_restart_point

    !> Physical coordinates for triggering MPI communication of a Lagrangian entity.
    impure subroutine s_set_lag_comm_coords()

        pcomm_coords(1)%beg = x_cb(-1)
        pcomm_coords(1)%end = x_cb(m)
        $:GPU_UPDATE(device='[pcomm_coords(1)]')
        if (n > 0) then
            pcomm_coords(2)%beg = y_cb(-1)
            pcomm_coords(2)%end = y_cb(n)
            $:GPU_UPDATE(device='[pcomm_coords(2)]')
            if (p > 0) then
                pcomm_coords(3)%beg = z_cb(-1)
                pcomm_coords(3)%end = z_cb(p)
                $:GPU_UPDATE(device='[pcomm_coords(3)]')
            end if
        end if

    end subroutine s_set_lag_comm_coords

    !> Create case_dir/D on rank 0 if it does not exist.
    impure subroutine s_create_D_dir()

        character(LEN=path_len + 2*name_len) :: path_D_dir
        logical                              :: dir_exist

        if (proc_rank == 0) then
            path_D_dir = trim(case_dir) // '/D'
            call my_inquire(trim(path_D_dir), dir_exist)
            if (.not. dir_exist) call s_create_directory(trim(path_D_dir))
        end if
        call s_mpi_barrier()

    end subroutine s_create_D_dir

    !> Read the initial Lagrangian entities from path, one per line with ncols values, and return those in this rank's physical
    !! domain: rows(1:n_loc, :) and their global IDs ids(1:n_loc), which are the line numbers among the n_read lines read.
    impure subroutine s_read_lag_input(path, ncols, n_glb, rows, ids, n_loc, n_read)

        character(len=*), intent(in)                       :: path
        integer, intent(in)                                :: ncols, n_glb
        real(wp), allocatable, dimension(:,:), intent(out) :: rows
        integer, allocatable, dimension(:), intent(out)    :: ids
        integer, intent(out)                               :: n_loc, n_read
        real(wp), dimension(ncols)                         :: vals
        character(LEN=1024)                                :: line
        character(LEN=200)                                 :: abort_msg
        integer                                            :: ios, ios_line
        logical                                            :: file_exist

        call my_inquire(trim(path), file_exist)
        if (.not. file_exist) call s_mpi_abort('Lagrangian input file ' // trim(path) // ' does not exist')

        allocate (rows(max(1, n_glb), ncols), ids(max(1, n_glb)))
        n_loc = 0
        n_read = 0
        open (94, file=trim(path), form='formatted', iostat=ios)
        do while (ios == 0)
            ! Read line by line so a line with extra columns cannot shift the lines after it
            read (94, '(A)', iostat=ios) line
            if (ios /= 0 .or. len_trim(line) == 0) cycle
            read (line, *, iostat=ios_line) vals
            if (ios_line /= 0) then
                write (abort_msg, '(A,I0,A)') 'Each line of ' // trim(path) // ' needs ', ncols, ' columns'
                call s_mpi_abort(trim(abort_msg))
            end if
            n_read = n_read + 1
            if (n_read > n_glb) call s_mpi_abort('More entries in ' // trim(path) // ' than the global count (nBubs_glb or ' &
                & // 'nParticles_glb)')
            if (particle_in_domain_physical(vals(1:3))) then
                n_loc = n_loc + 1
                rows(n_loc,:) = vals
                ids(n_loc) = n_read
            end if
        end do
        close (94)

    end subroutine s_read_lag_input

    !> Sum the local count n_loc over the ranks into n_glb and abort if no entity lies in the domain.
    impure subroutine s_count_lag_glb(n_loc, n_glb, path)

        integer, intent(in)          :: n_loc
        integer, intent(out)         :: n_glb
        character(len=*), intent(in) :: path

        if (num_procs > 1) then
            call s_mpi_reduce_int_sum(n_loc, n_glb)
        else
            n_glb = n_loc
        end if
        if (proc_rank == 0 .and. n_glb == 0) call s_mpi_abort('No Lagrangian entities in the domain. Check ' // trim(path))

    end subroutine s_count_lag_glb

    !> Open this rank's evolution file D/lag_<lag_kind>_evol_<rank>.dat: on a fresh start replace it and write the column headers
    !! (first header width 14, the others width), on a restart append to it.
    impure subroutine s_open_lag_evol(lag_kind, width, headers)

        character(len=*), intent(in)               :: lag_kind
        integer, intent(in)                        :: width
        character(len=*), dimension(:), intent(in) :: headers
        character(LEN=path_len + 2*name_len)       :: file_loc
        character(LEN=32)                          :: FMT
        integer                                    :: save_count, i
        real(wp)                                   :: qtime

        write (file_loc, '(A,A,A,I0,A)') 'lag_', lag_kind, '_evol_', proc_rank, '.dat'
        file_loc = trim(case_dir) // '/D/' // trim(file_loc)

        call s_get_lag_restart_point(save_count, qtime)
        if (save_count > 0) then
            open (LAG_EVOL_ID, FILE=trim(file_loc), form='formatted', position='append')
        else
            open (LAG_EVOL_ID, FILE=trim(file_loc), form='formatted', position='rewind', status='replace')
            write (FMT, '(A,I0,A,I0,A,I0,A)') '(A', width, ',A14,', size(headers) - 2, 'A', width, ')'
            write (LAG_EVOL_ID, FMT) (trim(headers(i)), i=1, size(headers))
        end if

    end subroutine s_open_lag_evol

    !> Locate the cell holding pos, starting the search from cell, and its computational coordinates scoord (cell index plus the
    !! fraction of the cell width). The search stops at the edge of the buffer region.
    subroutine s_locate_cell(pos, cell, scoord)

        $:GPU_ROUTINE(function_name='s_locate_cell',parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in)   :: pos
        real(wp), dimension(3), intent(out)  :: scoord
        integer, dimension(3), intent(inout) :: cell

        do while (pos(1) < x_cb(cell(1) - 1) .and. cell(1) > -buff_size)
            cell(1) = cell(1) - 1
        end do
        do while (pos(1) >= x_cb(cell(1)) .and. cell(1) < m + buff_size)
            cell(1) = cell(1) + 1
        end do

        do while (pos(2) < y_cb(cell(2) - 1) .and. cell(2) > -buff_size)
            cell(2) = cell(2) - 1
        end do
        do while (pos(2) >= y_cb(cell(2)) .and. cell(2) < n + buff_size)
            cell(2) = cell(2) + 1
        end do

        if (p > 0) then
            do while (pos(3) < z_cb(cell(3) - 1) .and. cell(3) > -buff_size)
                cell(3) = cell(3) - 1
            end do
            do while (pos(3) >= z_cb(cell(3)) .and. cell(3) < p + buff_size)
                cell(3) = cell(3) + 1
            end do
        else
            cell(3) = 0
        end if

        ! Cell 0 starts at the domain boundary, so the center of cell i is at scoord = i + 1/2
        scoord(1) = cell(1) + (pos(1) - x_cb(cell(1) - 1))/dx(cell(1))
        scoord(2) = cell(2) + (pos(2) - y_cb(cell(2) - 1))/dy(cell(2))
        scoord(3) = 0._wp
        if (p > 0) scoord(3) = cell(3) + (pos(3) - z_cb(cell(3) - 1))/dz(cell(3))

    end subroutine s_locate_cell

    !> True if pos lies in this rank's physical domain (ghost cells excluded).
    function particle_in_domain_physical(pos_part)

        $:GPU_ROUTINE(parallelism='[seq]')

        logical                            :: particle_in_domain_physical
        real(wp), dimension(3), intent(in) :: pos_part

        particle_in_domain_physical = ((pos_part(1) < x_cb(m)) .and. (pos_part(1) >= x_cb(-1)) .and. (pos_part(2) < y_cb(n)) &
                                       & .and. (pos_part(2) >= y_cb(-1)))

        if (p > 0) then
            particle_in_domain_physical = (particle_in_domain_physical .and. (pos_part(3) < z_cb(p)) .and. (pos_part(3) &
                                           & >= z_cb(-1)))
        end if

    end function particle_in_domain_physical

    !> Volume of cell (cellx, celly, cellz); in planar 2D the cell has depth charwidth.
    subroutine s_get_char_vol(cellx, celly, cellz, charwidth, Charvol)

        $:GPU_ROUTINE(function_name='s_get_char_vol',parallelism='[seq]', cray_inline=True)

        integer, intent(in)   :: cellx, celly, cellz
        real(wp), intent(in)  :: charwidth
        real(wp), intent(out) :: Charvol

        if (p > 0) then
            Charvol = dx(cellx)*dy(celly)*dz(cellz)
        else
            if (cyl_coord) then
                Charvol = dx(cellx)*dy(celly)*y_cc(celly)*2._wp*pi
            else
                Charvol = dx(cellx)*dy(celly)*charwidth
            end if
        end if

    end subroutine s_get_char_vol

    !> Open D/voidfraction.dat on rank 0, appending if it exists.
    impure subroutine s_open_void_evol

        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist

        if (proc_rank == 0) then
            file_loc = trim(case_dir) // '/D/voidfraction.dat'
            call my_inquire(trim(file_loc), file_exist)
            if (.not. file_exist) then
                open (LAG_VOID_ID, FILE=trim(file_loc), form='formatted', position='rewind')
            else
                open (LAG_VOID_ID, FILE=trim(file_loc), form='formatted', position='append')
            end if
        end if

    end subroutine s_open_void_evol

    !> Append the cloud's average and maximum Lagrangian volume fraction 1 - alpha_f and its total volume at qtime to
    !! D/voidfraction.dat. The average is over the cells the cloud occupies only.
    impure subroutine s_write_void_evol(qtime, alpha_f, charwidth)

        real(wp), intent(in)                                                              :: qtime, charwidth
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:), intent(in) :: alpha_f
        real(wp)                                                                          :: volcell, voltot
        real(wp)                                                                          :: lag_void_max, lag_void_avg, lag_vol
        real(wp)                                                                          :: void_max_glb, void_avg_glb, vol_glb
        integer                                                                           :: i, j, k

        lag_void_max = 0._wp
        lag_void_avg = 0._wp
        lag_vol = 0._wp
        $:GPU_PARALLEL_LOOP(private='[volcell]', collapse=3, reduction='[[lag_vol, lag_void_avg], [lag_void_max]]', &
                            & reductionOp='[+, MAX]', copy='[lag_vol, lag_void_avg, lag_void_max]', copyin='[charwidth]')
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    lag_void_max = max(lag_void_max, 1._wp - alpha_f(i, j, k))
                    call s_get_char_vol(i, j, k, charwidth, volcell)
                    if ((1._wp - alpha_f(i, j, k)) > void_stat_min) then
                        lag_void_avg = lag_void_avg + (1._wp - alpha_f(i, j, k))*volcell
                        lag_vol = lag_vol + volcell
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

#ifdef MFC_MPI
        if (num_procs > 1) then
            call s_mpi_allreduce_max(lag_void_max, void_max_glb)
            lag_void_max = void_max_glb
            call s_mpi_allreduce_sum(lag_vol, vol_glb)
            lag_vol = vol_glb
            call s_mpi_allreduce_sum(lag_void_avg, void_avg_glb)
            lag_void_avg = void_avg_glb
        end if
#endif
        voltot = lag_void_avg
        if (lag_vol > 0._wp) lag_void_avg = lag_void_avg/lag_vol

        if (proc_rank == 0) then
            write (LAG_VOID_ID, '(6X,4e24.8)') qtime, lag_void_avg, lag_void_max, voltot
        end if

    end subroutine s_write_void_evol

    !> Close D/voidfraction.dat.
    impure subroutine s_close_void_evol

        if (proc_rank == 0) close (LAG_VOID_ID)

    end subroutine s_close_void_evol

    !> Close this rank's Lagrangian evolution file.
    impure subroutine s_close_lag_evol

        close (LAG_EVOL_ID)

    end subroutine s_close_lag_evol

    !> Write restart_data/lag_<lag_kind>_<t_step>.dat: a header (total count, time, dt, number of ranks, per-rank counts) followed
    !! by the lag_io_vars columns of every rank's entities, which the caller packs into io_data(1:n_loc, :).
    impure subroutine s_write_lag_restart(lag_kind, t_step, io_data, n_loc)

        character(len=*), intent(in)         :: lag_kind
        integer, intent(in)                  :: t_step, n_loc
        real(wp), dimension(:,:), intent(in) :: io_data

#ifdef MFC_MPI
        character(LEN=path_len + 2*name_len)   :: file_loc
        logical                                :: file_exist
        integer                                :: tot_part, ifile, ierr, view
        integer, dimension(MPI_STATUS_SIZE)    :: status
        integer(KIND=MPI_OFFSET_KIND)          :: disp
        integer, dimension(2)                  :: gsizes, lsizes, start_idx_part
        integer, dimension(num_procs)          :: proc_counts
        real(wp), dimension(1:1,1:lag_io_vars) :: dummy

        dummy = 0._wp
        lsizes(1) = n_loc
        lsizes(2) = lag_io_vars

        call MPI_ALLREDUCE(n_loc, tot_part, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLGATHER(n_loc, 1, MPI_INTEGER, proc_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

        ! Starting row of this rank's entities
        call MPI_EXSCAN(lsizes(1), start_idx_part(1), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (proc_rank == 0) start_idx_part(1) = 0
        start_idx_part(2) = 0

        gsizes(1) = tot_part
        gsizes(2) = lag_io_vars

        write (file_loc, '(A,A,A,I0,A)') 'lag_', lag_kind, '_', t_step, '.dat'
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)

        if (proc_rank == 0) then
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist) call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
            call s_check_mpi_file_open(ierr, file_loc)
            call MPI_FILE_WRITE(ifile, tot_part, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_WRITE(ifile, mytime, 1, mpi_p, status, ierr)
            call MPI_FILE_WRITE(ifile, dt, 1, mpi_p, status, ierr)
            call MPI_FILE_WRITE(ifile, num_procs, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_WRITE(ifile, proc_counts, num_procs, MPI_INTEGER, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if (n_loc > 0) then
            call MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, MPI_ORDER_FORTRAN, mpi_p, view, ierr)
        else
            call MPI_TYPE_CONTIGUOUS(0, mpi_p, view, ierr)
        end if
        call MPI_TYPE_COMMIT(view, ierr)

        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
        call s_check_mpi_file_open(ierr, file_loc)

        ! Skip the header written by rank 0
        disp = int(sizeof(tot_part) + 2*sizeof(mytime) + sizeof(num_procs) + num_procs*sizeof(proc_counts(1)), MPI_OFFSET_KIND)
        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

        if (n_loc > 0) then
            call MPI_FILE_WRITE_ALL(ifile, io_data(1:n_loc,:), lag_io_vars*n_loc, mpi_p, status, ierr)
        else
            call MPI_FILE_WRITE_ALL(ifile, dummy, 0, mpi_p, status, ierr)
        end if

        call MPI_FILE_CLOSE(ifile, ierr)
        call MPI_TYPE_FREE(view, ierr)
#endif

    end subroutine s_write_lag_restart

    !> Read restart_data/lag_<lag_kind>_<save_count>.dat written by s_write_lag_restart: set mytime and dt from the header and
    !! return this rank's n_loc rows in io_data(1:n_loc, 1:lag_io_vars). Each rank reads what the same rank wrote, so the run must
    !! use the same number of ranks as the one that wrote the file.
    impure subroutine s_read_lag_restart(lag_kind, save_count, io_data, n_loc)

        character(len=*), intent(in)                       :: lag_kind
        integer, intent(in)                                :: save_count
        real(wp), allocatable, dimension(:,:), intent(out) :: io_data
        integer, intent(out)                               :: n_loc
        character(LEN=path_len + 2*name_len)               :: file_loc
        character(len=200)                                 :: abort_msg
        logical                                            :: file_exist

#ifndef MFC_MPI
        n_loc = 0
        @:PROHIBIT(.true., "Lagrangian restart requires MPI (--mpi)")
#else
        real(wp)                               :: file_time, file_dt
        integer                                :: file_num_procs, file_tot_part, ifile, ierr, view, i
        integer, dimension(MPI_STATUS_SIZE)    :: status
        integer(kind=MPI_OFFSET_KIND)          :: disp
        integer, dimension(2)                  :: gsizes, lsizes, start_idx_part
        integer, dimension(:), allocatable     :: proc_counts
        real(wp), dimension(1:1,1:lag_io_vars) :: dummy

        dummy = 0._wp
        n_loc = 0

        write (file_loc, '(A,A,A,I0,A)') 'lag_', lag_kind, '_', save_count, '.dat'
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)

        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) call s_mpi_abort('Restart file ' // trim(file_loc) // ' does not exist!')

        if (.not. parallel_io) return

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call s_check_mpi_file_open(ierr, file_loc)
            call MPI_FILE_READ(ifile, file_tot_part, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_READ(ifile, file_time, 1, mpi_p, status, ierr)
            call MPI_FILE_READ(ifile, file_dt, 1, mpi_p, status, ierr)
            call MPI_FILE_READ(ifile, file_num_procs, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        call MPI_BCAST(file_tot_part, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(file_time, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(file_dt, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(file_num_procs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        if (file_num_procs /= num_procs) then
            write (abort_msg, &
                   & '(A,I0,A,I0)') 'Lagrangian restart needs the same number of ranks as the run that wrote the ' &
                   & // 'restart file. Ranks in the file: ', file_num_procs, ', ranks in this run: ', num_procs
            call s_mpi_abort(trim(abort_msg))
        end if

        allocate (proc_counts(file_num_procs))
        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call s_check_mpi_file_open(ierr, file_loc)
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs), MPI_OFFSET_KIND)
            call MPI_FILE_SEEK(ifile, disp, MPI_SEEK_SET, ierr)
            call MPI_FILE_READ(ifile, proc_counts, file_num_procs, MPI_INTEGER, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        end if
        call MPI_BCAST(proc_counts, file_num_procs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        mytime = file_time
        dt = file_dt

        n_loc = proc_counts(proc_rank + 1)
        start_idx_part(1) = 0
        do i = 1, proc_rank
            start_idx_part(1) = start_idx_part(1) + proc_counts(i)
        end do
        start_idx_part(2) = 0
        lsizes(1) = n_loc
        lsizes(2) = lag_io_vars
        gsizes(1) = file_tot_part
        gsizes(2) = lag_io_vars

        allocate (io_data(max(1, n_loc),1:lag_io_vars))

        if (n_loc > 0) then
            call MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, MPI_ORDER_FORTRAN, mpi_p, view, ierr)
        else
            call MPI_TYPE_CONTIGUOUS(0, mpi_p, view, ierr)
        end if
        call MPI_TYPE_COMMIT(view, ierr)

        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
        call s_check_mpi_file_open(ierr, file_loc)

        ! Skip the header
        disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs) + file_num_procs*sizeof(proc_counts(1)), &
                   & MPI_OFFSET_KIND)
        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

        if (n_loc > 0) then
            call MPI_FILE_READ_ALL(ifile, io_data, lag_io_vars*n_loc, mpi_p, status, ierr)
        else
            call MPI_FILE_READ_ALL(ifile, dummy, 0, mpi_p, status, ierr)
        end if

        call MPI_FILE_CLOSE(ifile, ierr)
        call MPI_TYPE_FREE(view, ierr)
        deallocate (proc_counts)

        if (proc_rank == 0) then
            write (*, '(A,I0,A,A,A,I0)') 'Read ', file_tot_part, ' Lagrangian ', lag_kind, ' from restart file at t_step = ', &
                   & save_count
            write (*, '(A,E15.7,A,E15.7)') 'Restart time = ', mytime, ', dt = ', dt
        end if
#endif

    end subroutine s_read_lag_restart

end module m_euler_lagrange
