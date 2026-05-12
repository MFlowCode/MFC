!>
!! @file
!! @brief Contains module m_data_output

!> @brief Writes post-processed grid and flow-variable data to Silo-HDF5 or binary database files
module m_data_output

    use m_derived_types
    use m_global_parameters
    use m_derived_variables
    use m_mpi_proxy
    use m_compile_specific
    use m_helper
    use m_variables_conversion

    implicit none

    private; public :: s_initialize_data_output_module, s_define_output_region, s_open_formatted_database_file, &
        & s_open_intf_data_file, s_open_energy_data_file, s_write_grid_to_formatted_database_file, &
        & s_write_variable_to_formatted_database_file, s_write_lag_bubbles_results_to_text, &
        & s_write_lag_bubbles_to_formatted_database_file, s_write_ib_state_files, s_write_intf_data_file, &
        & s_write_energy_data_file, s_write_ib_bodies_to_formatted_database_file, s_close_formatted_database_file, &
        & s_close_intf_data_file, s_close_energy_data_file, s_finalize_data_output_module

    ! Include Silo-HDF5 interface library
    include 'silo_f9x.inc'

    ! Flow variable storage; q_root_sf gathers to rank 0 in 1D parallel runs
    real(wp), allocatable, dimension(:,:,:), public :: q_sf
    real(wp), allocatable, dimension(:,:,:)         :: q_root_sf
    real(wp), allocatable, dimension(:,:,:)         :: cyl_q_sf

    ! Single precision storage for flow variables
    real(sp), allocatable, dimension(:,:,:), public :: q_sf_s
    real(sp), allocatable, dimension(:,:,:)         :: q_root_sf_s
    real(sp), allocatable, dimension(:,:,:)         :: cyl_q_sf_s

    ! Spatial and data extents for VisIt visualization
    real(wp), allocatable, dimension(:,:) :: spatial_extents
    real(wp), allocatable, dimension(:,:) :: data_extents

    ! Ghost zone layer sizes (lo/hi) for subdomain connectivity in VisIt
    integer, allocatable, dimension(:) :: lo_offset
    integer, allocatable, dimension(:) :: hi_offset

    ! Track cell-boundary count per active coordinate direction
    integer, allocatable, dimension(:) :: dims

    ! Locations of various folders in the case's directory tree, associated with the choice of the formatted database format. These
    ! include, in order, the location of the folder named after the selected formatted database format, and the locations of two
    ! sub-directories of the latter, the first of which is named after the local processor rank, while the second is named 'root'.
    ! The folder associated with the local processor rank contains only the data pertaining to the part of the domain taken care of
    ! by the local processor. The root directory, on the other hand, will contain either the information about the connectivity
    ! required to put the entire domain back together, or the actual data associated with the entire computational domain. This all
    ! depends on dimensionality and the choice of the formatted database format.
    character(LEN=path_len + name_len)   :: dbdir
    character(LEN=path_len + 2*name_len) :: proc_rank_dir
    character(LEN=path_len + 2*name_len) :: rootdir

    ! Handles of the formatted database master/root file, slave/local processor file and options list. The list of options is
    ! explicitly used in the Silo- HDF5 database format to provide additional details about the contents of a formatted database
    ! file, such as the previously described spatial and data extents.
    integer :: dbroot
    integer :: dbfile
    integer :: optlist

    ! The total number of flow variable(s) to be stored in a formatted database file. Note that this is only needed when using the
    ! Binary format.
    integer :: dbvars

    ! Generic error flags utilized in the handling, checking and the reporting of the input and output operations errors with a
    ! formatted database file
    integer, private :: err

contains

    !> Allocate storage arrays, configure output directories, and count flow variables for formatted database output.
    impure subroutine s_initialize_data_output_module()

        character(LEN=len_trim(case_dir) + 2*name_len) :: file_loc
        logical                                        :: dir_check
        integer                                        :: i

        allocate (q_sf(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end))
        if (grid_geometry == 3) then
            allocate (cyl_q_sf(-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end,-offset_x%beg:m + offset_x%end))
        end if

        if (precision == 1) then
            allocate (q_sf_s(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end))
            if (grid_geometry == 3) then
                allocate (cyl_q_sf_s(-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end,-offset_x%beg:m + offset_x%end))
            end if
        end if

        if (n == 0) then
            allocate (q_root_sf(0:m_root,0:0,0:0))
            if (precision == 1) then
                allocate (q_root_sf_s(0:m_root,0:0,0:0))
            end if
        end if

        ! Allocating the spatial and data extents and also the variables for the offsets and the one bookkeeping the number of
        ! cell-boundaries in each active coordinate direction. Note that all these variables are only needed by the Silo-HDF5 format
        ! for multidimensional data.
        if (format == 1) then
            allocate (data_extents(1:2,0:num_procs - 1))

            if (p > 0) then
                allocate (spatial_extents(1:6,0:num_procs - 1))
                allocate (lo_offset(1:3))
                allocate (hi_offset(1:3))
                allocate (dims(1:3))
            else if (n > 0) then
                allocate (spatial_extents(1:4,0:num_procs - 1))
                allocate (lo_offset(1:2))
                allocate (hi_offset(1:2))
                allocate (dims(1:2))
            else
                allocate (spatial_extents(1:2,0:num_procs - 1))
                allocate (lo_offset(1:1))
                allocate (hi_offset(1:1))
                allocate (dims(1:1))
            end if
        end if

        ! The size of the ghost zone layer in each of the active coordinate directions was set in the module m_mpi_proxy.f90. The
        ! results are now transferred to the local variables of this module when they are required by the Silo-HDF5 format, for
        ! multidimensional data sets. With the same, latter, requirements, the variables bookkeeping the number of cell-boundaries
        ! in each active coordinate direction are also set here.
        if (format == 1) then
            if (p > 0) then
                if (grid_geometry == 3) then
                    lo_offset(:) = (/offset_y%beg, offset_z%beg, offset_x%beg/)
                    hi_offset(:) = (/offset_y%end, offset_z%end, offset_x%end/)
                else
                    lo_offset(:) = (/offset_x%beg, offset_y%beg, offset_z%beg/)
                    hi_offset(:) = (/offset_x%end, offset_y%end, offset_z%end/)
                end if

                if (grid_geometry == 3) then
                    dims(:) = (/n + offset_y%beg + offset_y%end + 2, p + offset_z%beg + offset_z%end + 2, &
                         & m + offset_x%beg + offset_x%end + 2/)
                else
                    dims(:) = (/m + offset_x%beg + offset_x%end + 2, n + offset_y%beg + offset_y%end + 2, &
                         & p + offset_z%beg + offset_z%end + 2/)
                end if
            else if (n > 0) then
                lo_offset(:) = (/offset_x%beg, offset_y%beg/)
                hi_offset(:) = (/offset_x%end, offset_y%end/)

                dims(:) = (/m + offset_x%beg + offset_x%end + 2, n + offset_y%beg + offset_y%end + 2/)
            else
                lo_offset(:) = (/offset_x%beg/)
                hi_offset(:) = (/offset_x%end/)
                dims(:) = (/m + offset_x%beg + offset_x%end + 2/)
            end if
        end if

        if (format == 1) then
            dbdir = trim(case_dir) // '/silo_hdf5'

            write (proc_rank_dir, '(A,I0)') '/p', proc_rank

            proc_rank_dir = trim(dbdir) // trim(proc_rank_dir)

            file_loc = trim(proc_rank_dir) // '/.'

            call my_inquire(file_loc, dir_check)
            if (dir_check .neqv. .true.) then
                call s_create_directory(trim(proc_rank_dir))
            end if

            if (proc_rank == 0) then
                rootdir = trim(dbdir) // '/root'

                file_loc = trim(rootdir) // '/.'

                call my_inquire(file_loc, dir_check)
                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(rootdir))
                end if
            end if
        else
            dbdir = trim(case_dir) // '/binary'

            write (proc_rank_dir, '(A,I0)') '/p', proc_rank

            proc_rank_dir = trim(dbdir) // trim(proc_rank_dir)

            file_loc = trim(proc_rank_dir) // '/.'

            call my_inquire(file_loc, dir_check)

            if (dir_check .neqv. .true.) then
                call s_create_directory(trim(proc_rank_dir))
            end if

            if (n == 0 .and. proc_rank == 0) then
                rootdir = trim(dbdir) // '/root'

                file_loc = trim(rootdir) // '/.'

                call my_inquire(file_loc, dir_check)

                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(rootdir))
                end if
            end if
        end if

        if (bubbles_lagrange) then  ! Lagrangian solver
            if (lag_txt_wrt) then
                dbdir = trim(case_dir) // '/lag_bubbles_post_process'
                file_loc = trim(dbdir) // '/.'
                call my_inquire(file_loc, dir_check)

                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(dbdir))
                end if
            end if
        end if

        ! Contrary to the Silo-HDF5 database format, handles of the Binary database master/root and slave/local process files are
        ! perfectly static throughout post-process. Hence, they are set here so that they do not have to be repetitively computed in
        ! later procedures.
        if (format == 2) then
            if (n == 0 .and. proc_rank == 0) dbroot = 2
            dbfile = 1
        end if

        if (format == 2) then
            dbvars = 0

            if ((model_eqns == 2) .or. (model_eqns == 3)) then
                do i = 1, num_fluids
                    if (alpha_rho_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                        dbvars = dbvars + 1
                    end if
                end do
            end if

            if ((rho_wrt .or. (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) .and. (.not. relativity)) then
                dbvars = dbvars + 1
            end if

            if (relativity .and. (rho_wrt .or. prim_vars_wrt)) dbvars = dbvars + 1
            if (relativity .and. (rho_wrt .or. cons_vars_wrt)) dbvars = dbvars + 1

            do i = 1, eqn_idx%E - eqn_idx%mom%beg
                if (mom_wrt(i) .or. cons_vars_wrt) dbvars = dbvars + 1
            end do

            do i = 1, eqn_idx%E - eqn_idx%mom%beg
                if (vel_wrt(i) .or. prim_vars_wrt) dbvars = dbvars + 1
            end do

            do i = 1, eqn_idx%E - eqn_idx%mom%beg
                if (flux_wrt(i)) dbvars = dbvars + 1
            end do

            if (E_wrt .or. cons_vars_wrt) dbvars = dbvars + 1
            if (pres_wrt .or. prim_vars_wrt) dbvars = dbvars + 1
            if (hypoelasticity) dbvars = dbvars + (num_dims*(num_dims + 1))/2
            if (cont_damage) dbvars = dbvars + 1
            if (hyper_cleaning) dbvars = dbvars + 1

            if (mhd) then
                if (n == 0) then
                    dbvars = dbvars + 2
                else
                    dbvars = dbvars + 3
                end if
            end if

            if ((model_eqns == 2) .or. (model_eqns == 3)) then
                do i = 1, num_fluids - 1
                    if (alpha_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                        dbvars = dbvars + 1
                    end if
                end do

                if (alpha_wrt(num_fluids) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                    dbvars = dbvars + 1
                end if
            end if

            if (gamma_wrt .or. (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) then
                dbvars = dbvars + 1
            end if

            if (heat_ratio_wrt) dbvars = dbvars + 1

            if (pi_inf_wrt .or. (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) then
                dbvars = dbvars + 1
            end if

            if (pres_inf_wrt) dbvars = dbvars + 1
            if (c_wrt) dbvars = dbvars + 1

            if (p > 0) then
                do i = 1, num_vels
                    if (omega_wrt(i)) dbvars = dbvars + 1
                end do
            else if (n > 0) then
                do i = 1, num_vels
                    if (omega_wrt(i)) dbvars = dbvars + 1
                end do
            end if

            if (schlieren_wrt) dbvars = dbvars + 1
        end if

    end subroutine s_initialize_data_output_module

    !> Compute the cell-index bounds for the user-specified partial output domain in each coordinate direction.
    impure subroutine s_define_output_region

        integer :: i
        integer :: lower_bound, upper_bound

        #:for X, M in [('x', 'm'), ('y', 'n'), ('z', 'p')]
            if (${M}$ == 0) return  ! Early return for y or z if simulation is 1D or 2D

            lower_bound = -offset_${X}$%beg
            upper_bound = ${M}$ + offset_${X}$%end

            do i = lower_bound, upper_bound
                if (${X}$%cc(i) > ${X}$_output%beg) then
                    ${X}$_output_idx%beg = i + offset_${X}$%beg
                    exit
                end if
            end do

            do i = upper_bound, lower_bound, -1
                if (${X}$%cc(i) < ${X}$_output%end) then
                    ${X}$_output_idx%end = i + offset_${X}$%beg
                    exit
                end if
            end do

            ! If no grid points are within the output region
            if ((${X}$%cc(lower_bound) > ${X}$_output%end) .or. (${X}$%cc(upper_bound) < ${X}$_output%beg)) then
                ${X}$_output_idx%beg = 0
                ${X}$_output_idx%end = 0
            end if
        #:endfor

    end subroutine s_define_output_region

    !> Open (or create) the Silo-HDF5 or Binary formatted database slave and master files for a given time step.
    impure subroutine s_open_formatted_database_file(t_step)

        integer, intent(in)                            :: t_step
        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc
        integer                                        :: ierr

        if (format == 1) then
            write (file_loc, '(A,I0,A)') '/', t_step, '.silo'
            file_loc = trim(proc_rank_dir) // trim(file_loc)

            ierr = DBCREATE(trim(file_loc), len_trim(file_loc), DB_CLOBBER, DB_LOCAL, 'MFC v3.0', 8, DB_HDF5, dbfile)

            if (dbfile == -1) then
                call s_mpi_abort('Unable to create Silo-HDF5 database ' // 'slave file ' // trim(file_loc) // '. ' // 'Exiting.')
            end if

            if (proc_rank == 0) then
                write (file_loc, '(A,I0,A)') '/collection_', t_step, '.silo'
                file_loc = trim(rootdir) // trim(file_loc)

                ierr = DBCREATE(trim(file_loc), len_trim(file_loc), DB_CLOBBER, DB_LOCAL, 'MFC v3.0', 8, DB_HDF5, dbroot)

                if (dbroot == -1) then
                    call s_mpi_abort('Unable to create Silo-HDF5 database ' // 'master file ' // trim(file_loc) // '. ' &
                                     & // 'Exiting.')
                end if
            end if
        else
            write (file_loc, '(A,I0,A)') '/', t_step, '.dat'
            file_loc = trim(proc_rank_dir) // trim(file_loc)

            open (dbfile, IOSTAT=err, FILE=trim(file_loc), form='unformatted', STATUS='replace')

            if (err /= 0) then
                call s_mpi_abort('Unable to create Binary database slave ' // 'file ' // trim(file_loc) // '. Exiting.')
            end if

            if (output_partial_domain) then
                write (dbfile) x_output_idx%end - x_output_idx%beg, y_output_idx%end - y_output_idx%beg, &
                       & z_output_idx%end - z_output_idx%beg, dbvars
            else
                write (dbfile) m, n, p, dbvars
            end if

            if (n == 0 .and. proc_rank == 0) then
                write (file_loc, '(A,I0,A)') '/', t_step, '.dat'
                file_loc = trim(rootdir) // trim(file_loc)

                open (dbroot, IOSTAT=err, FILE=trim(file_loc), form='unformatted', STATUS='replace')

                if (err /= 0) then
                    call s_mpi_abort('Unable to create Binary database ' // 'master file ' // trim(file_loc) // '. Exiting.')
                end if

                if (output_partial_domain) then
                    write (dbroot) x_output_idx%end - x_output_idx%beg, 0, 0, dbvars
                else
                    write (dbroot) m_root, 0, 0, dbvars
                end if
            end if
        end if

    end subroutine s_open_formatted_database_file

    !> Open the interface data file for appending extracted interface coordinates.
    impure subroutine s_open_intf_data_file()

        character(LEN=path_len + 3*name_len) :: file_path

        write (file_path, '(A)') '/intf_data.dat'
        file_path = trim(case_dir) // trim(file_path)

        open (211, FILE=trim(file_path), form='formatted', POSITION='append', STATUS='unknown')

    end subroutine s_open_intf_data_file

    !> Open the energy data file for appending volume-integrated energy budget quantities.
    impure subroutine s_open_energy_data_file()

        character(LEN=path_len + 3*name_len) :: file_path

        write (file_path, '(A)') '/eng_data.dat'
        file_path = trim(case_dir) // trim(file_path)

        open (251, FILE=trim(file_path), form='formatted', POSITION='append', STATUS='unknown')

    end subroutine s_open_energy_data_file

    !> Write the computational grid (cell-boundary coordinates) to the formatted database slave and master files.
    impure subroutine s_write_grid_to_formatted_database_file(t_step)

        integer, intent(in) :: t_step

        ! NAG compiler requires these to be statically sized
        character(LEN=4*name_len), dimension(num_procs) :: meshnames
        integer, dimension(num_procs)                   :: meshtypes
        integer                                         :: i
        integer                                         :: ierr

        if (format == 1) then
            ! For multidimensional data sets, the spatial extents of all of the grid(s) handled by the local processor(s) are
            ! recorded so that they may be written, by root processor, to the formatted database master file.
            if (num_procs > 1) then
                call s_mpi_gather_spatial_extents(spatial_extents)
            else if (p > 0) then
                if (grid_geometry == 3) then
                    spatial_extents(:,0) = (/minval(y%cb), minval(z%cb), minval(x%cb), maxval(y%cb), maxval(z%cb), maxval(x%cb)/)
                else
                    spatial_extents(:,0) = (/minval(x%cb), minval(y%cb), minval(z%cb), maxval(x%cb), maxval(y%cb), maxval(z%cb)/)
                end if
            else if (n > 0) then
                spatial_extents(:,0) = (/minval(x%cb), minval(y%cb), maxval(x%cb), maxval(y%cb)/)
            else
                spatial_extents(:,0) = (/minval(x%cb), maxval(x%cb)/)
            end if

            ! Next, the root processor proceeds to record all of the spatial extents in the formatted database master file. In
            ! addition, it also records a sub-domain connectivity map so that the entire grid may be reassembled by looking at the
            ! master file.
            if (proc_rank == 0) then
                do i = 1, num_procs
                    write (meshnames(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:rectilinear_grid'
                end do

                meshtypes = DB_QUAD_RECT

                err = DBSET2DSTRLEN(len(meshnames(1)))
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIOPT(optlist, DBOPT_EXTENTS_SIZE, size(spatial_extents, 1))
                err = DBADDDOPT(optlist, DBOPT_EXTENTS, spatial_extents)
                err = DBPUTMMESH(dbroot, 'rectilinear_grid', 16, num_procs, meshnames, len_trim(meshnames), meshtypes, optlist, &
                                 & ierr)
                err = DBFREEOPTLIST(optlist)
            end if

            ! Finally, the local quadrilateral mesh, either 2D or 3D, along with its offsets that indicate the presence and size of
            ! ghost zone layer(s), are put in the formatted database slave file.

            if (p > 0) then
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIAOPT(optlist, DBOPT_LO_OFFSET, size(lo_offset), lo_offset)
                err = DBADDIAOPT(optlist, DBOPT_HI_OFFSET, size(hi_offset), hi_offset)
                if (grid_geometry == 3) then
                    err = DBPUTQM(dbfile, 'rectilinear_grid', 16, 'x', 1, 'y', 1, 'z', 1, y%cb, z%cb, x%cb, dims, 3, DB_DOUBLE, &
                                  & DB_COLLINEAR, optlist, ierr)
                else
                    err = DBPUTQM(dbfile, 'rectilinear_grid', 16, 'x', 1, 'y', 1, 'z', 1, x%cb, y%cb, z%cb, dims, 3, DB_DOUBLE, &
                                  & DB_COLLINEAR, optlist, ierr)
                end if
                err = DBFREEOPTLIST(optlist)
            else if (n > 0) then
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIAOPT(optlist, DBOPT_LO_OFFSET, size(lo_offset), lo_offset)
                err = DBADDIAOPT(optlist, DBOPT_HI_OFFSET, size(hi_offset), hi_offset)
                err = DBPUTQM(dbfile, 'rectilinear_grid', 16, 'x', 1, 'y', 1, 'z', 1, x%cb, y%cb, DB_F77NULL, dims, 2, DB_DOUBLE, &
                              & DB_COLLINEAR, optlist, ierr)
                err = DBFREEOPTLIST(optlist)
            else
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIAOPT(optlist, DBOPT_LO_OFFSET, size(lo_offset), lo_offset)
                err = DBADDIAOPT(optlist, DBOPT_HI_OFFSET, size(hi_offset), hi_offset)
                err = DBPUTQM(dbfile, 'rectilinear_grid', 16, 'x', 1, 'y', 1, 'z', 1, x%cb, DB_F77NULL, DB_F77NULL, dims, 1, &
                              & DB_DOUBLE, DB_COLLINEAR, optlist, ierr)
                err = DBFREEOPTLIST(optlist)
            end if
        else if (format == 2) then
            ! Multidimensional local grid data is written to the formatted database slave file. Recall that no master file to
            ! maintained in multidimensions.
            if (p > 0) then
                if (precision == 1) then
                    write (dbfile) real(x%cb, sp), real(y%cb, sp), real(z%cb, sp)
                else
                    if (output_partial_domain) then
                        write (dbfile) x%cb(x_output_idx%beg - 1:x_output_idx%end), y%cb(y_output_idx%beg - 1:y_output_idx%end), &
                               & z%cb(z_output_idx%beg - 1:z_output_idx%end)
                    else
                        write (dbfile) x%cb, y%cb, z%cb
                    end if
                end if
            else if (n > 0) then
                if (precision == 1) then
                    write (dbfile) real(x%cb, sp), real(y%cb, sp)
                else
                    if (output_partial_domain) then
                        write (dbfile) x%cb(x_output_idx%beg - 1:x_output_idx%end), y%cb(y_output_idx%beg - 1:y_output_idx%end)
                    else
                        write (dbfile) x%cb, y%cb
                    end if
                end if

                ! One-dimensional local grid data is written to the formatted database slave file. In addition, the local grid data
                ! is put together by the root process and written to the master file.
            else
                if (precision == 1) then
                    write (dbfile) real(x%cb, sp)
                else
                    if (output_partial_domain) then
                        write (dbfile) x%cb(x_output_idx%beg - 1:x_output_idx%end)
                    else
                        write (dbfile) x%cb
                    end if
                end if

                if (num_procs > 1) then
                    call s_mpi_defragment_1d_grid_variable()
                else
                    x_root_cb(:) = x%cb(:)
                end if

                if (proc_rank == 0) then
                    if (precision == 1) then
                        write (dbroot) real(x_root_cb, wp)
                    else
                        if (output_partial_domain) then
                            write (dbroot) x_root_cb(x_output_idx%beg - 1:x_output_idx%end)
                        else
                            write (dbroot) x_root_cb
                        end if
                    end if
                end if
            end if
        end if

    end subroutine s_write_grid_to_formatted_database_file

    !> Write a single flow variable field to the formatted database slave and master files for a given time step.
    impure subroutine s_write_variable_to_formatted_database_file(varname, t_step)

        character(LEN=*), intent(in) :: varname
        integer, intent(in)          :: t_step

        ! NAG compiler requires these to be statically sized
        character(LEN=4*name_len), dimension(num_procs) :: varnames
        integer, dimension(num_procs)                   :: vartypes
        integer                                         :: i, j, k
        integer                                         :: ierr

        if (format == 1) then
            ! Determining the extents of the flow variable on each local process and gathering all this information on root process
            if (num_procs > 1) then
                call s_mpi_gather_data_extents(q_sf, data_extents)
            else
                data_extents(:,0) = (/minval(q_sf), maxval(q_sf)/)
            end if

            if (proc_rank == 0) then
                do i = 1, num_procs
                    write (varnames(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:' // trim(varname)
                end do

                vartypes = DB_QUADVAR

                err = DBSET2DSTRLEN(len(varnames(1)))
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIOPT(optlist, DBOPT_EXTENTS_SIZE, 2)
                err = DBADDDOPT(optlist, DBOPT_EXTENTS, data_extents)
                err = DBPUTMVAR(dbroot, trim(varname), len_trim(varname), num_procs, varnames, len_trim(varnames), vartypes, &
                                & optlist, ierr)
                err = DBFREEOPTLIST(optlist)
            end if

            if (wp == dp) then
                if (precision == 1) then
                    do i = -offset_x%beg, m + offset_x%end
                        do j = -offset_y%beg, n + offset_y%end
                            do k = -offset_z%beg, p + offset_z%end
                                q_sf_s(i, j, k) = real(q_sf(i, j, k), sp)
                            end do
                        end do
                    end do
                    if (grid_geometry == 3) then
                        do i = -offset_x%beg, m + offset_x%end
                            do j = -offset_y%beg, n + offset_y%end
                                do k = -offset_z%beg, p + offset_z%end
                                    cyl_q_sf_s(j, k, i) = q_sf_s(i, j, k)
                                end do
                            end do
                        end do
                    end if
                else
                    if (grid_geometry == 3) then
                        do i = -offset_x%beg, m + offset_x%end
                            do j = -offset_y%beg, n + offset_y%end
                                do k = -offset_z%beg, p + offset_z%end
                                    cyl_q_sf(j, k, i) = q_sf(i, j, k)
                                end do
                            end do
                        end do
                    end if
                end if
            else if (wp == sp) then
                do i = -offset_x%beg, m + offset_x%end
                    do j = -offset_y%beg, n + offset_y%end
                        do k = -offset_z%beg, p + offset_z%end
                            q_sf_s(i, j, k) = q_sf(i, j, k)
                        end do
                    end do
                end do
                if (grid_geometry == 3) then
                    do i = -offset_x%beg, m + offset_x%end
                        do j = -offset_y%beg, n + offset_y%end
                            do k = -offset_z%beg, p + offset_z%end
                                cyl_q_sf_s(j, k, i) = q_sf_s(i, j, k)
                            end do
                        end do
                    end do
                end if
            end if

            #:for PRECISION, SFX, DBT in [(1,'_s','DB_FLOAT'),(2,'',"DB_DOUBLE")]
                if (precision == ${PRECISION}$) then
                    if (p > 0) then
                        if (grid_geometry == 3) then
                            err = DBPUTQV1(dbfile, trim(varname), len_trim(varname), 'rectilinear_grid', 16, cyl_q_sf${SFX}$, &
                                           & dims - 1, 3, DB_F77NULL, 0, ${DBT}$, DB_ZONECENT, DB_F77NULL, ierr)
                        else
                            err = DBPUTQV1(dbfile, trim(varname), len_trim(varname), 'rectilinear_grid', 16, q_sf${SFX}$, &
                                           & dims - 1, 3, DB_F77NULL, 0, ${DBT}$, DB_ZONECENT, DB_F77NULL, ierr)
                        end if
                    else if (n > 0) then
                        err = DBPUTQV1(dbfile, trim(varname), len_trim(varname), 'rectilinear_grid', 16, q_sf${SFX}$, dims - 1, &
                                       & 2, DB_F77NULL, 0, ${DBT}$, DB_ZONECENT, DB_F77NULL, ierr)
                    else
                        err = DBPUTQV1(dbfile, trim(varname), len_trim(varname), 'rectilinear_grid', 16, q_sf${SFX}$, dims - 1, &
                                       & 1, DB_F77NULL, 0, ${DBT}$, DB_ZONECENT, DB_F77NULL, ierr)
                    end if
                end if
            #:endfor
        else
            ! Writing the name of the flow variable and its data, associated with the local processor, to the formatted database
            ! slave file
            if (precision == 1) then
                write (dbfile) varname, real(q_sf, wp)
            else
                write (dbfile) varname, q_sf
            end if

            ! In 1D, the root process also takes care of gathering the flow variable data from all of the local processor(s) and
            ! writes it to the formatted database master file.
            if (n == 0) then
                if (num_procs > 1) then
                    call s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)
                else
                    q_root_sf(:,:,:) = q_sf(:,:,:)
                end if

                if (proc_rank == 0) then
                    if (precision == 1) then
                        write (dbroot) varname, real(q_root_sf, wp)
                    else
                        write (dbroot) varname, q_root_sf
                    end if
                end if
            end if
        end if

    end subroutine s_write_variable_to_formatted_database_file

    !> Write the post-processed results in the folder 'lag_bubbles_data'
    impure subroutine s_write_lag_bubbles_results_to_text(t_step)

        integer, intent(in)                            :: t_step
        character(len=len_trim(case_dir) + 3*name_len) :: file_loc
        integer                                        :: id

#ifdef MFC_MPI
        real(wp), dimension(20)                :: inputvals
        real(wp)                               :: time_real
        integer, dimension(MPI_STATUS_SIZE)    :: status
        integer(KIND=MPI_OFFSET_KIND)          :: disp
        integer                                :: view
        logical                                :: file_exist
        integer, dimension(2)                  :: gsizes, lsizes, start_idx_part
        integer                                :: ifile
        integer                                :: ierr
        real(wp)                               :: file_time, file_dt
        integer                                :: file_num_procs, file_tot_part
        integer                                :: i
        integer, dimension(:), allocatable     :: proc_bubble_counts
        real(wp), dimension(1:1,1:lag_io_vars) :: lag_io_null

        lag_io_null = 0._wp

        ! Construct file path
        write (file_loc, '(A,I0,A)') 'lag_bubbles_', t_step, '.dat'
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)

        ! Check if file exists
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort('Restart file ' // trim(file_loc) // ' does not exist!')
        end if

        if (.not. parallel_io) return

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

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
        time_real = file_time

        allocate (proc_bubble_counts(file_num_procs))

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Skip to processor counts position
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs), MPI_OFFSET_KIND)
            call MPI_FILE_SEEK(ifile, disp, MPI_SEEK_SET, ierr)
            call MPI_FILE_READ(ifile, proc_bubble_counts, file_num_procs, MPI_INTEGER, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        call MPI_BCAST(proc_bubble_counts, file_num_procs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        gsizes(1) = file_tot_part
        gsizes(2) = lag_io_vars
        lsizes(1) = file_tot_part
        lsizes(2) = lag_io_vars
        start_idx_part(1) = 0
        start_idx_part(2) = 0

        call MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, MPI_ORDER_FORTRAN, mpi_p, view, ierr)
        call MPI_TYPE_COMMIT(view, ierr)

        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

        disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs) &
                   & + file_num_procs*sizeof(proc_bubble_counts(1)), MPI_OFFSET_KIND)
        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_null, ierr)

        allocate (MPI_IO_DATA_lg_bubbles(file_tot_part,1:lag_io_vars))

        call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA_lg_bubbles, lag_io_vars*file_tot_part, mpi_p, status, ierr)

        write (file_loc, '(A,I0,A)') 'lag_bubbles_post_process_', t_step, '.dat'
        file_loc = trim(case_dir) // '/lag_bubbles_post_process/' // trim(file_loc)

        if (proc_rank == 0) then
            open (unit=29, file=file_loc, form='formatted', position='rewind')

            if (lag_header) then
                write (29, '(A)', advance='no')
                if (lag_id_wrt) write (29, '(A8)', advance='no') 'id, '
                if (lag_pos_wrt) write (29, '(3(A17))', advance='no') 'px, ', 'py, ', 'pz, '
                if (lag_pos_prev_wrt) write (29, '(3(A17))', advance='no') 'pvx, ', 'pvy, ', 'pvz, '
                if (lag_vel_wrt) write (29, '(3(A17))', advance='no') 'vx, ', 'vy, ', 'vz, '
                if (lag_rad_wrt) write (29, '(A17)', advance='no') 'radius, '
                if (lag_rvel_wrt) write (29, '(A17)', advance='no') 'rvel, '
                if (lag_r0_wrt) write (29, '(A17)', advance='no') 'r0, '
                if (lag_rmax_wrt) write (29, '(A17)', advance='no') 'rmax, '
                if (lag_rmin_wrt) write (29, '(A17)', advance='no') 'rmin, '
                if (lag_dphidt_wrt) write (29, '(A17)', advance='no') 'dphidt, '
                if (lag_pres_wrt) write (29, '(A17)', advance='no') 'pressure, '
                if (lag_mv_wrt) write (29, '(A17)', advance='no') 'mv, '
                if (lag_mg_wrt) write (29, '(A17)', advance='no') 'mg, '
                if (lag_betaT_wrt) write (29, '(A17)', advance='no') 'betaT, '
                if (lag_betaC_wrt) write (29, '(A17)', advance='no') 'betaC, '
                write (29, '(A15)') 'time'
            end if

            do i = 1, file_tot_part
                id = int(MPI_IO_DATA_lg_bubbles(i, 1))
                inputvals(1:20) = MPI_IO_DATA_lg_bubbles(i,2:21)
                if (id > 0) then
                    write (29, '(100(A))', advance='no') ''
                    if (lag_id_wrt) write (29, '(I6, A)', advance='no') id, ', '
                    if (lag_pos_wrt) write (29, '(3(E15.7, A))', advance='no') inputvals(1), ', ', inputvals(2), ', ', &
                        & inputvals(3), ', '
                    if (lag_pos_prev_wrt) write (29, '(3(E15.7, A))', advance='no') inputvals(4), ', ', inputvals(5), ', ', &
                        & inputvals(6), ', '
                    if (lag_vel_wrt) write (29, '(3(E15.7, A))', advance='no') inputvals(7), ', ', inputvals(8), ', ', &
                        & inputvals(9), ', '
                    if (lag_rad_wrt) write (29, '(E15.7, A)', advance='no') inputvals(10), ', '
                    if (lag_rvel_wrt) write (29, '(E15.7, A)', advance='no') inputvals(11), ', '
                    if (lag_r0_wrt) write (29, '(E15.7, A)', advance='no') inputvals(12), ', '
                    if (lag_rmax_wrt) write (29, '(E15.7, A)', advance='no') inputvals(13), ', '
                    if (lag_rmin_wrt) write (29, '(E15.7, A)', advance='no') inputvals(14), ', '
                    if (lag_dphidt_wrt) write (29, '(E15.7, A)', advance='no') inputvals(15), ', '
                    if (lag_pres_wrt) write (29, '(E15.7, A)', advance='no') inputvals(16), ', '
                    if (lag_mv_wrt) write (29, '(E15.7, A)', advance='no') inputvals(17), ', '
                    if (lag_mg_wrt) write (29, '(E15.7, A)', advance='no') inputvals(18), ', '
                    if (lag_betaT_wrt) write (29, '(E15.7, A)', advance='no') inputvals(19), ', '
                    if (lag_betaC_wrt) write (29, '(E15.7, A)', advance='no') inputvals(20), ', '
                    write (29, '(E15.7)') time_real
                end if
            end do
            close (29)
        end if

        deallocate (MPI_IO_DATA_lg_bubbles)

        call s_mpi_barrier()

        call MPI_FILE_CLOSE(ifile, ierr)
#endif

    end subroutine s_write_lag_bubbles_results_to_text

    !> Read Lagrangian bubble restart data and write bubble positions and scalar fields to the Silo database.
    impure subroutine s_write_lag_bubbles_to_formatted_database_file(t_step)

        integer, intent(in)                            :: t_step
        character(len=len_trim(case_dir) + 3*name_len) :: file_loc
        integer                                        :: id

#ifdef MFC_MPI
        real(wp)                                        :: time_real
        integer, dimension(MPI_STATUS_SIZE)             :: status
        integer(KIND=MPI_OFFSET_KIND)                   :: disp
        integer                                         :: view
        logical                                         :: file_exist
        integer, dimension(2)                           :: gsizes, lsizes, start_idx_part
        integer                                         :: ifile, ierr, nBub
        real(wp)                                        :: file_time, file_dt
        integer                                         :: file_num_procs, file_tot_part
        integer, dimension(:), allocatable              :: proc_bubble_counts
        real(wp), dimension(1:1,1:lag_io_vars)          :: dummy
        character(LEN=4*name_len), dimension(num_procs) :: meshnames
        integer, dimension(num_procs)                   :: meshtypes
        real(wp)                                        :: dummy_data
        integer                                         :: i
        real(wp), dimension(:), allocatable             :: bub_id
        real(wp), dimension(:), allocatable             :: px, py, pz, ppx, ppy, ppz, vx, vy, vz
        real(wp), dimension(:), allocatable             :: radius, rvel, rnot, rmax, rmin, dphidt
        real(wp), dimension(:), allocatable             :: pressure, mv, mg, betaT, betaC

        dummy = 0._wp
        dummy_data = 0._wp

        ! Construct file path
        write (file_loc, '(A,I0,A)') 'lag_bubbles_', t_step, '.dat'
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)

        ! Check if file exists
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort('Restart file ' // trim(file_loc) // ' does not exist!')
        end if

        if (.not. parallel_io) return

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

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
        time_real = file_time

        allocate (proc_bubble_counts(file_num_procs))

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Skip to processor counts position
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs), MPI_OFFSET_KIND)
            call MPI_FILE_SEEK(ifile, disp, MPI_SEEK_SET, ierr)
            call MPI_FILE_READ(ifile, proc_bubble_counts, file_num_procs, MPI_INTEGER, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        call MPI_BCAST(proc_bubble_counts, file_num_procs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! Set time variables from file

        nBub = proc_bubble_counts(proc_rank + 1)

        start_idx_part(1) = 0
        do i = 1, proc_rank
            start_idx_part(1) = start_idx_part(1) + proc_bubble_counts(i)
        end do

        start_idx_part(2) = 0
        lsizes(1) = nBub
        lsizes(2) = lag_io_vars

        gsizes(1) = file_tot_part
        gsizes(2) = lag_io_vars

        if (nBub > 0) then
            #:for VAR in ['bub_id', 'px', 'py', 'pz', 'ppx', 'ppy', 'ppz', 'vx', 'vy', 'vz', &
                'radius', 'rvel', 'rnot', 'rmax', 'rmin', 'dphidt', &
                'pressure', 'mv', 'mg', 'betaT', 'betaC']
                allocate (${VAR}$ (nBub))
            #:endfor
            allocate (MPI_IO_DATA_lg_bubbles(nBub,1:lag_io_vars))

            call MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, MPI_ORDER_FORTRAN, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Skip extended header
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs) &
                       & + file_num_procs*sizeof(proc_bubble_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA_lg_bubbles, lag_io_vars*nBub, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
            call MPI_TYPE_FREE(view, ierr)

            ! Extract data from MPI_IO_DATA_lg_bubbles array Adjust these indices based on your actual data layout
            #:for VAR, IDX in [('bub_id', 1), ('px', 2), ('py',3), ('pz',4), ('ppx',5), ('ppy',6), ('ppz',7), &
                ('vx',8), ('vy',9), ('vz',10), ('radius',11), ('rvel',12), &
                ('rnot',13), ('rmax',14), ('rmin',15), ('dphidt',16), &
                ('pressure',17), ('mv',18), ('mg',19), ('betaT',20), ('betaC',21)]
                ${VAR}$ (:) = MPI_IO_DATA_lg_bubbles(:,${IDX}$)
            #:endfor

            ! Next, the root processor proceeds to record all of the spatial extents in the formatted database master file. In
            ! addition, it also records a sub-domain connectivity map so that the entire grid may be reassembled by looking at the
            ! master file.
            if (proc_rank == 0) then
                do i = 1, num_procs
                    write (meshnames(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:lag_bubbles'
                    meshtypes(i) = DB_POINTMESH
                end do
                err = DBSET2DSTRLEN(len(meshnames(1)))
                err = DBPUTMMESH(dbroot, 'lag_bubbles', 16, num_procs, meshnames, len_trim(meshnames), meshtypes, DB_F77NULL, ierr)
            end if

            err = DBPUTPM(dbfile, 'lag_bubbles', 11, 3, px, py, pz, nBub, DB_DOUBLE, DB_F77NULL, ierr)

            if (lag_id_wrt) call s_write_lag_variable_to_formatted_database_file('part_id', t_step, bub_id, nBub)
            if (lag_vel_wrt) then
                call s_write_lag_variable_to_formatted_database_file('part_vel1', t_step, vx, nBub)
                call s_write_lag_variable_to_formatted_database_file('part_vel2', t_step, vy, nBub)
                if (p > 0) call s_write_lag_variable_to_formatted_database_file('part_vel3', t_step, vz, nBub)
            end if
            if (lag_rad_wrt) call s_write_lag_variable_to_formatted_database_file('part_radius', t_step, radius, nBub)
            if (lag_rvel_wrt) call s_write_lag_variable_to_formatted_database_file('part_rdot', t_step, rvel, nBub)
            if (lag_r0_wrt) call s_write_lag_variable_to_formatted_database_file('part_r0', t_step, rnot, nBub)
            if (lag_rmax_wrt) call s_write_lag_variable_to_formatted_database_file('part_rmax', t_step, rmax, nBub)
            if (lag_rmin_wrt) call s_write_lag_variable_to_formatted_database_file('part_rmin', t_step, rmin, nBub)
            if (lag_dphidt_wrt) call s_write_lag_variable_to_formatted_database_file('part_dphidt', t_step, dphidt, nBub)
            if (lag_pres_wrt) call s_write_lag_variable_to_formatted_database_file('part_pressure', t_step, pressure, nBub)
            if (lag_mv_wrt) call s_write_lag_variable_to_formatted_database_file('part_mv', t_step, mv, nBub)
            if (lag_mg_wrt) call s_write_lag_variable_to_formatted_database_file('part_mg', t_step, mg, nBub)
            if (lag_betaT_wrt) call s_write_lag_variable_to_formatted_database_file('part_betaT', t_step, betaT, nBub)
            if (lag_betaC_wrt) call s_write_lag_variable_to_formatted_database_file('part_betaC', t_step, betaC, nBub)

            deallocate (bub_id, px, py, pz, ppx, ppy, ppz, vx, vy, vz, radius, rvel, rnot, rmax, rmin, dphidt, pressure, mv, mg, &
                        & betaT, betaC)
            deallocate (MPI_IO_DATA_lg_bubbles)
        else
            call MPI_TYPE_CONTIGUOUS(0, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Skip extended header
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs) &
                       & + file_num_procs*sizeof(proc_bubble_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_READ_ALL(ifile, dummy, 0, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
            call MPI_TYPE_FREE(view, ierr)

            if (proc_rank == 0) then
                do i = 1, num_procs
                    write (meshnames(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:lag_bubbles'
                    meshtypes(i) = DB_POINTMESH
                end do
                err = DBSET2DSTRLEN(len(meshnames(1)))
                err = DBPUTMMESH(dbroot, 'lag_bubbles', 16, num_procs, meshnames, len_trim(meshnames), meshtypes, DB_F77NULL, ierr)
            end if

            err = DBSETEMPTYOK(1)
            err = DBPUTPM(dbfile, 'lag_bubbles', 11, 3, dummy_data, dummy_data, dummy_data, 0, DB_DOUBLE, DB_F77NULL, ierr)

            if (lag_id_wrt) call s_write_lag_variable_to_formatted_database_file('part_id', t_step)
            if (lag_vel_wrt) then
                call s_write_lag_variable_to_formatted_database_file('part_vel1', t_step)
                call s_write_lag_variable_to_formatted_database_file('part_vel2', t_step)
                if (p > 0) call s_write_lag_variable_to_formatted_database_file('part_vel3', t_step)
            end if
            if (lag_rad_wrt) call s_write_lag_variable_to_formatted_database_file('part_radius', t_step)
            if (lag_rvel_wrt) call s_write_lag_variable_to_formatted_database_file('part_rdot', t_step)
            if (lag_r0_wrt) call s_write_lag_variable_to_formatted_database_file('part_r0', t_step)
            if (lag_rmax_wrt) call s_write_lag_variable_to_formatted_database_file('part_rmax', t_step)
            if (lag_rmin_wrt) call s_write_lag_variable_to_formatted_database_file('part_rmin', t_step)
            if (lag_dphidt_wrt) call s_write_lag_variable_to_formatted_database_file('part_dphidt', t_step)
            if (lag_pres_wrt) call s_write_lag_variable_to_formatted_database_file('part_pressure', t_step)
            if (lag_mv_wrt) call s_write_lag_variable_to_formatted_database_file('part_mv', t_step)
            if (lag_mg_wrt) call s_write_lag_variable_to_formatted_database_file('part_mg', t_step)
            if (lag_betaT_wrt) call s_write_lag_variable_to_formatted_database_file('part_betaT', t_step)
            if (lag_betaC_wrt) call s_write_lag_variable_to_formatted_database_file('part_betaC', t_step)
        end if
#endif

    end subroutine s_write_lag_bubbles_to_formatted_database_file

    !> Write a single Lagrangian bubble point-variable to the Silo database slave and master files.
    subroutine s_write_lag_variable_to_formatted_database_file(varname, t_step, data, nBubs)

        character(len=*), intent(in)                  :: varname
        integer, intent(in)                           :: t_step
        real(wp), dimension(1:), intent(in), optional :: data
        integer, intent(in), optional                 :: nBubs
        character(len=64), dimension(num_procs)       :: var_names
        integer, dimension(num_procs)                 :: var_types
        real(wp)                                      :: dummy_data
        integer                                       :: ierr
        integer                                       :: i

        dummy_data = 0._wp

        if (present(nBubs) .and. present(data)) then
            if (proc_rank == 0) then
                do i = 1, num_procs
                    write (var_names(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:' // trim(varname)
                    var_types(i) = DB_POINTVAR
                end do
                err = DBSET2DSTRLEN(len(var_names(1)))
                err = DBPUTMVAR(dbroot, trim(varname), len_trim(varname), num_procs, var_names, len_trim(var_names), var_types, &
                                & DB_F77NULL, ierr)
            end if

            err = DBPUTPV1(dbfile, trim(varname), len_trim(varname), 'lag_bubbles', 11, data, nBubs, DB_DOUBLE, DB_F77NULL, ierr)
        else
            if (proc_rank == 0) then
                do i = 1, num_procs
                    write (var_names(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:' // trim(varname)
                    var_types(i) = DB_POINTVAR
                end do
                err = DBSET2DSTRLEN(len(var_names(1)))
                err = DBSETEMPTYOK(1)
                err = DBPUTMVAR(dbroot, trim(varname), len_trim(varname), num_procs, var_names, len_trim(var_names), var_types, &
                                & DB_F77NULL, ierr)
            end if

            err = DBSETEMPTYOK(1)
            err = DBPUTPV1(dbfile, trim(varname), len_trim(varname), 'lag_bubbles', 11, dummy_data, 0, DB_DOUBLE, DB_F77NULL, ierr)
        end if

    end subroutine s_write_lag_variable_to_formatted_database_file

    !> Convert the binary immersed-boundary state file to per-body formatted text files
    impure subroutine s_write_ib_state_files()

        character(len=len_trim(case_dir) + 4*name_len) :: in_file, out_file, file_loc
        integer                                        :: iu_in, ios, i, rec_id
        integer, allocatable, dimension(:)             :: iu_out
        real(wp)                                       :: rec_time
        real(wp), dimension(3)                         :: rec_force, rec_torque
        real(wp), dimension(3)                         :: rec_vel, rec_angular_vel
        real(wp), dimension(3)                         :: rec_angles, rec_centroid

        file_loc = trim(case_dir) // '/D'

        in_file = trim(file_loc) // '/ib_state.dat'
        open (newunit=iu_in, file=trim(in_file), form='unformatted', access='stream', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            call s_mpi_abort('Cannot open IB state input file: ' // trim(in_file))
        end if

        allocate (iu_out(num_ibs))
        do i = 1, num_ibs
            write (out_file, '(A,I0,A)') trim(file_loc) // '/ib_', i, '.txt'
            open (newunit=iu_out(i), file=trim(out_file), form='formatted', status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                call s_mpi_abort('Cannot open IB state output file: ' // trim(out_file))
            end if
            write (iu_out(i), &
                   & '(A)') 'mytime fx fy fz Tau_x Tau_y Tau_z vx vy vz omega_x omega_y omega_z angle_x angle_y angle_z x_c y_c z_c'
        end do

        do
            read (iu_in, iostat=ios) rec_time, rec_id, rec_force, rec_torque, rec_vel, rec_angular_vel, rec_angles, &
                  & rec_centroid(1), rec_centroid(2), rec_centroid(3)
            if (ios /= 0) exit

            if (rec_id >= 1 .and. rec_id <= num_ibs) then
                write (iu_out(rec_id), '(19(ES24.16E3,1X))') rec_time, rec_force(1), rec_force(2), rec_force(3), rec_torque(1), &
                       & rec_torque(2), rec_torque(3), rec_vel(1), rec_vel(2), rec_vel(3), rec_angular_vel(1), &
                       & rec_angular_vel(2), rec_angular_vel(3), rec_angles(1), rec_angles(2), rec_angles(3), rec_centroid(1), &
                       & rec_centroid(2), rec_centroid(3)
            end if
        end do

        close (iu_in)
        do i = 1, num_ibs
            close (iu_out(i))
        end do
        deallocate (iu_out)

    end subroutine s_write_ib_state_files

    !> Extract the volume-fraction interface contour from primitive fields and write the coordinates to the interface data file.
    impure subroutine s_write_intf_data_file(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer :: i, j, k, l, cent
        integer :: counter, root  !< number of data points extracted to fit shape to SH perturbations
        real(wp), allocatable :: x_td(:), y_td(:), x_d1(:), y_d1(:), y_d(:), x_d(:)
        real(wp) :: axp, axm, ayp, aym, tgp, euc_d, thres, maxalph_loc, maxalph_glb

        allocate (x_d1(m*n))
        allocate (y_d1(m*n))
        counter = 0
        maxalph_loc = 0._wp
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    if (q_prim_vf(eqn_idx%E + 2)%sf(i, j, k) > maxalph_loc) then
                        maxalph_loc = q_prim_vf(eqn_idx%E + 2)%sf(i, j, k)
                    end if
                end do
            end do
        end do

        call s_mpi_allreduce_max(maxalph_loc, maxalph_glb)
        if (p > 0) then
            do l = 0, p
                if (z%cc(l) < z%spacing(l) .and. z%cc(l) > 0) then
                    cent = l
                end if
            end do
        else
            cent = 0
        end if

        thres = 0.9_wp*maxalph_glb
        do k = 0, n
            do j = 0, m
                axp = q_prim_vf(eqn_idx%E + 2)%sf(j + 1, k, cent)
                axm = q_prim_vf(eqn_idx%E + 2)%sf(j, k, cent)
                ayp = q_prim_vf(eqn_idx%E + 2)%sf(j, k + 1, cent)
                aym = q_prim_vf(eqn_idx%E + 2)%sf(j, k, cent)
                if ((axp > thres .and. axm < thres) .or. (axp < thres .and. axm > thres) .or. (ayp > thres .and. aym < thres) &
                    & .or. (ayp < thres .and. aym > thres)) then
                    if (counter == 0) then
                        counter = counter + 1
                        x_d1(counter) = x%cc(j)
                        y_d1(counter) = y%cc(k)
                    else
                        tgp = sqrt(x%spacing(j)**2 + y%spacing(k)**2)
                        do i = 1, counter
                            euc_d = sqrt((x%cc(j) - x_d1(i))**2 + (y%cc(k) - y_d1(i))**2)
                            if (euc_d < tgp) then
                                exit
                            else if (i == counter) then
                                counter = counter + 1
                                x_d1(counter) = x%cc(j)
                                y_d1(counter) = y%cc(k)
                            end if
                        end do
                    end if
                end if
            end do
        end do

        allocate (x_d(counter), y_d(counter))

        do i = 1, counter
            y_d(i) = y_d1(i)
            x_d(i) = x_d1(i)
        end do
        root = 0

        call s_mpi_gather_data(x_d, counter, x_td, root)
        call s_mpi_gather_data(y_d, counter, y_td, root)
        if (proc_rank == 0) then
            do i = 1, size(x_td)
                if (i == size(x_td)) then
                    write (211, '(F12.9,1X,F12.9,1X,I4)') x_td(i), y_td(i), size(x_td)
                else
                    write (211, '(F12.9,1X,F12.9,1X,F3.1)') x_td(i), y_td(i), 0._wp
                end if
            end do
        end if

    end subroutine s_write_intf_data_file

    !> Compute volume-integrated kinetic, potential, and internal energies and write the energy budget to the energy data file.
    impure subroutine s_write_energy_data_file(q_prim_vf, q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf, q_cons_vf
        real(wp) :: Elk, Egk, Elp, Egint, Vb, Vl, pres_av, Et
        real(wp) :: rho, pres, dV, tmp, gamma, pi_inf, MaxMa, MaxMa_glb, maxvel, c, Ma, H, qv
        real(wp), dimension(num_vels) :: vel
        real(wp), dimension(num_fluids) :: adv
        integer :: i, j, k, l, s  !< looping indices

        Egk = 0._wp
        Elp = 0._wp
        Egint = 0._wp
        Vb = 0._wp
        maxvel = 0._wp
        MaxMa = 0._wp
        Vl = 0._wp
        Elk = 0._wp
        Et = 0._wp
        Vb = 0._wp
        dV = 0._wp
        pres_av = 0._wp
        pres = 0._wp
        c = 0._wp

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    pres = 0._wp
                    dV = x%spacing(i)*y%spacing(j)*z%spacing(k)
                    rho = 0._wp
                    gamma = 0._wp
                    pi_inf = 0._wp
                    qv = 0._wp
                    pres = q_prim_vf(eqn_idx%E)%sf(i, j, k)
                    Egint = Egint + q_prim_vf(eqn_idx%E + 2)%sf(i, j, k)*(gammas(2)*pres)*dV
                    do s = 1, num_vels
                        vel(s) = q_prim_vf(num_fluids + s)%sf(i, j, k)
                        Egk = Egk + 0.5_wp*q_prim_vf(eqn_idx%E + 2)%sf(i, j, k)*q_prim_vf(2)%sf(i, j, k)*vel(s)*vel(s)*dV
                        Elk = Elk + 0.5_wp*q_prim_vf(eqn_idx%E + 1)%sf(i, j, k)*q_prim_vf(1)%sf(i, j, k)*vel(s)*vel(s)*dV
                        if (abs(vel(s)) > maxvel) then
                            maxvel = abs(vel(s))
                        end if
                    end do
                    do l = 1, eqn_idx%adv%end - eqn_idx%E
                        adv(l) = q_prim_vf(eqn_idx%E + l)%sf(i, j, k)
                        gamma = gamma + adv(l)*gammas(l)
                        pi_inf = pi_inf + adv(l)*pi_infs(l)
                        rho = rho + adv(l)*q_prim_vf(l)%sf(i, j, k)
                        qv = qv + adv(l)*q_prim_vf(l)%sf(i, j, k)*qvs(l)
                    end do

                    H = ((gamma + 1._wp)*pres + pi_inf + qv)/rho

                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, 0._wp, 0._wp, c, qv)

                    Ma = maxvel/c
                    if (Ma > MaxMa .and. (adv(1) > (1.0_wp - 1.0e-10_wp))) then
                        MaxMa = Ma
                    end if
                    Vl = Vl + adv(1)*dV
                    Vb = Vb + adv(2)*dV
                    pres_av = pres_av + adv(1)*pres*dV
                    Et = Et + q_cons_vf(eqn_idx%E)%sf(i, j, k)*dV
                end do
            end do
        end do

        tmp = pres_av
        call s_mpi_allreduce_sum(tmp, pres_av)
        tmp = Vl
        call s_mpi_allreduce_sum(tmp, Vl)

        call s_mpi_allreduce_max(MaxMa, MaxMa_glb)
        tmp = Elk
        call s_mpi_allreduce_sum(tmp, Elk)
        tmp = Egint
        call s_mpi_allreduce_sum(tmp, Egint)
        tmp = Egk
        call s_mpi_allreduce_sum(tmp, Egk)
        tmp = Vb
        call s_mpi_allreduce_sum(tmp, Vb)
        tmp = Et
        call s_mpi_allreduce_sum(tmp, Et)

        Elp = pres_av/Vl*Vb
        if (proc_rank == 0) then
            write (251, '(10X, 8F24.8)') Elp, Egint, Elk, Egk, Et, Vb, Vl, MaxMa_glb
        end if

    end subroutine s_write_energy_data_file

    !> Read IB state and write a Silo point mesh with per-body scalar fields.
    impure subroutine s_write_ib_bodies_to_formatted_database_file(t_step)

        integer, intent(in)                            :: t_step
        character(len=len_trim(case_dir) + 3*name_len) :: file_loc

#ifdef MFC_MPI
        integer, parameter                              :: NFIELDS_PER_IB = 20
        real(wp)                                        :: ib_buf(NFIELDS_PER_IB)
        real(wp), dimension(:,:), allocatable           :: ib_data
        logical                                         :: file_exist
        character(LEN=4*name_len), dimension(num_procs) :: meshnames
        integer, dimension(num_procs)                   :: meshtypes
        integer                                         :: i, ios, file_unit
        integer                                         :: ierr, nBodies
        real(wp), dimension(:), allocatable             :: px, py, pz
        real(wp), dimension(:), allocatable             :: force_x, force_y, force_z
        real(wp), dimension(:), allocatable             :: torque_x, torque_y, torque_z
        real(wp), dimension(:), allocatable             :: vel_x, vel_y, vel_z
        real(wp), dimension(:), allocatable             :: omega_x, omega_y, omega_z
        real(wp), dimension(:), allocatable             :: angle_x, angle_y, angle_z
        real(wp), dimension(:), allocatable             :: ib_diameter

        ! Build path to per-timestep IB state file
        write (file_loc, '(A,I0,A)') '/restart_data/ib_state_', t_step, '.dat'
        file_loc = trim(case_dir) // trim(file_loc)

        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort('Restart file ' // trim(file_loc) // ' does not exist!')
        end if

        nBodies = num_ibs

        if (nBodies > 0) then
            allocate (ib_data(nBodies, NFIELDS_PER_IB))
            allocate (px(nBodies), py(nBodies), pz(nBodies))
            allocate (force_x(nBodies), force_y(nBodies), force_z(nBodies))
            allocate (torque_x(nBodies), torque_y(nBodies), torque_z(nBodies))
            allocate (vel_x(nBodies), vel_y(nBodies), vel_z(nBodies))
            allocate (omega_x(nBodies), omega_y(nBodies), omega_z(nBodies))
            allocate (angle_x(nBodies), angle_y(nBodies), angle_z(nBodies))
            allocate (ib_diameter(nBodies))

            if (proc_rank == 0) then
                open (newunit=file_unit, file=trim(file_loc), form='unformatted', access='stream', status='old', iostat=ios)
                if (ios /= 0) call s_mpi_abort('Cannot open IB state file: ' // trim(file_loc))

                do i = 1, nBodies
                    read (file_unit, iostat=ios) ib_buf
                    if (ios /= 0) call s_mpi_abort('Error reading IB state file')
                    ib_data(i,:) = ib_buf(:)
                end do

                close (file_unit)
            end if

            call MPI_BCAST(ib_data, nBodies*NFIELDS_PER_IB, mpi_p, 0, MPI_COMM_WORLD, ierr)

            do i = 1, nBodies
                force_x(i) = ib_data(i, 2); force_y(i) = ib_data(i, 3); force_z(i) = ib_data(i, 4)
                torque_x(i) = ib_data(i, 5); torque_y(i) = ib_data(i, 6); torque_z(i) = ib_data(i, 7)
                vel_x(i) = ib_data(i, 8); vel_y(i) = ib_data(i, 9); vel_z(i) = ib_data(i, 10)
                omega_x(i) = ib_data(i, 11); omega_y(i) = ib_data(i, 12); omega_z(i) = ib_data(i, 13)
                angle_x(i) = ib_data(i, 14); angle_y(i) = ib_data(i, 15); angle_z(i) = ib_data(i, 16)
                px(i) = ib_data(i, 17); py(i) = ib_data(i, 18); pz(i) = ib_data(i, 19)
                ib_diameter(i) = ib_data(i, 20)*2.0_wp
            end do

            if (proc_rank == 0) then
                do i = 1, num_procs
                    write (meshnames(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:ib_bodies'
                    meshtypes(i) = DB_POINTMESH
                end do
                err = DBSET2DSTRLEN(len(meshnames(1)))
                err = DBPUTMMESH(dbroot, 'ib_bodies', 16, num_procs, meshnames, len_trim(meshnames), meshtypes, DB_F77NULL, ierr)
            end if

            err = DBPUTPM(dbfile, 'ib_bodies', 9, 3, px, py, pz, nBodies, DB_DOUBLE, DB_F77NULL, ierr)

            call s_write_ib_variable('ib_force_x', t_step, force_x, nBodies)
            call s_write_ib_variable('ib_force_y', t_step, force_y, nBodies)
            call s_write_ib_variable('ib_force_z', t_step, force_z, nBodies)
            call s_write_ib_variable('ib_torque_x', t_step, torque_x, nBodies)
            call s_write_ib_variable('ib_torque_y', t_step, torque_y, nBodies)
            call s_write_ib_variable('ib_torque_z', t_step, torque_z, nBodies)
            call s_write_ib_variable('ib_vel_x', t_step, vel_x, nBodies)
            call s_write_ib_variable('ib_vel_y', t_step, vel_y, nBodies)
            call s_write_ib_variable('ib_vel_z', t_step, vel_z, nBodies)
            call s_write_ib_variable('ib_omega_x', t_step, omega_x, nBodies)
            call s_write_ib_variable('ib_omega_y', t_step, omega_y, nBodies)
            call s_write_ib_variable('ib_omega_z', t_step, omega_z, nBodies)
            call s_write_ib_variable('ib_angle_x', t_step, angle_x, nBodies)
            call s_write_ib_variable('ib_angle_y', t_step, angle_y, nBodies)
            call s_write_ib_variable('ib_angle_z', t_step, angle_z, nBodies)
            call s_write_ib_variable('ib_diameter', t_step, ib_diameter, nBodies)

            deallocate (ib_data, px, py, pz, force_x, force_y, force_z)
            deallocate (torque_x, torque_y, torque_z, vel_x, vel_y, vel_z)
            deallocate (omega_x, omega_y, omega_z, angle_x, angle_y, angle_z)
            deallocate (ib_diameter)
        end if
#endif

    end subroutine s_write_ib_bodies_to_formatted_database_file

    !> Write a single IB point-variable to the Silo database slave and master files.
    subroutine s_write_ib_variable(varname, t_step, data, nBodies)

        character(len=*), intent(in)                    :: varname
        integer, intent(in)                             :: t_step
        real(wp), dimension(:), intent(in)              :: data
        integer, intent(in)                             :: nBodies
        character(len=4*name_len), dimension(num_procs) :: var_names
        integer, dimension(num_procs)                   :: var_types
        integer                                         :: ierr, i

        if (proc_rank == 0) then
            do i = 1, num_procs
                write (var_names(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:' // trim(varname)
                var_types(i) = DB_POINTVAR
            end do
            err = DBSET2DSTRLEN(len(var_names(1)))
            err = DBPUTMVAR(dbroot, trim(varname), len_trim(varname), num_procs, var_names, len_trim(var_names), var_types, &
                            & DB_F77NULL, ierr)
        end if

        err = DBPUTPV1(dbfile, trim(varname), len_trim(varname), 'ib_bodies', 9, data, nBodies, DB_DOUBLE, DB_F77NULL, ierr)

    end subroutine s_write_ib_variable

    !> Close the formatted database slave file and, for the root process, the master file.
    impure subroutine s_close_formatted_database_file()

        integer :: ierr

        if (format == 1) then
            ierr = DBCLOSE(dbfile)
            if (proc_rank == 0) ierr = DBCLOSE(dbroot)
        else
            close (dbfile)
            if (n == 0 .and. proc_rank == 0) close (dbroot)
        end if

    end subroutine s_close_formatted_database_file

    !> Close the interface data file.
    impure subroutine s_close_intf_data_file()

        close (211)

    end subroutine s_close_intf_data_file

    !> Close the energy data file.
    impure subroutine s_close_energy_data_file()

        close (251)

    end subroutine s_close_energy_data_file

    !> Deallocate module arrays and release all data-output resources.
    impure subroutine s_finalize_data_output_module()

        deallocate (q_sf)
        if (n == 0) deallocate (q_root_sf)
        if (grid_geometry == 3) then
            deallocate (cyl_q_sf)
        end if

        ! Deallocating spatial and data extents and also the variables for the offsets and the one bookkeeping the number of
        ! cell-boundaries in each active coordinate direction. Note that all these variables were only needed by Silo-HDF5 format
        ! for multidimensional data.
        if (format == 1) then
            deallocate (spatial_extents)
            deallocate (data_extents)
            deallocate (lo_offset)
            deallocate (hi_offset)
            deallocate (dims)
        end if

    end subroutine s_finalize_data_output_module

end module m_data_output
