!>
!! @file m_data_output.f90
!! @brief Contains module m_data_output

!> @brief This module enables the restructuring of the raw simulation data
!!              file(s) into formatted database file(s). The formats that may be
!!              chosen from include Silo-HDF5 and Binary. Each of these database
!!              structures contains information about the grid as well as each of
!!              the flow variable(s) that were chosen by the user to be included.
module m_data_output

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_derived_variables     !< Procedures used to compute quantities derived

    use m_mpi_proxy             ! Message passing interface (MPI) module proxy

    use m_compile_specific

    use m_helper

    implicit none

    private; public :: s_initialize_data_output_module, &
 s_define_output_region, &
 s_open_formatted_database_file, &
 s_open_intf_data_file, &
 s_open_energy_data_file, &
 s_write_grid_to_formatted_database_file, &
 s_write_variable_to_formatted_database_file, &
 s_write_lag_bubbles_results, &
 s_write_intf_data_file, &
 s_write_energy_data_file, &
 s_close_formatted_database_file, &
 s_close_intf_data_file, &
 s_close_energy_data_file, &
 s_finalize_data_output_module

    ! Including the Silo Fortran interface library that features the subroutines
    ! and parameters that are required to write in the Silo-HDF5 database format
    ! INCLUDE 'silo.inc'
    include 'silo_f9x.inc'

    ! Generic storage for flow variable(s) that are to be written to formatted
    ! database file(s). Note that for 1D simulations, q_root_sf is employed to
    ! gather the flow variable(s) from all sub-domains on to the root process.
    ! If the run is not parallel, but serial, then q_root_sf is equal to q_sf.
    real(wp), allocatable, dimension(:, :, :), public :: q_sf
    real(wp), allocatable, dimension(:, :, :) :: q_root_sf
    real(wp), allocatable, dimension(:, :, :) :: cyl_q_sf
    ! Single precision storage for flow variables
    real(sp), allocatable, dimension(:, :, :), public :: q_sf_s
    real(sp), allocatable, dimension(:, :, :) :: q_root_sf_s
    real(sp), allocatable, dimension(:, :, :) :: cyl_q_sf_s

    ! The spatial and data extents array variables contain information about the
    ! minimum and maximum values of the grid and flow variable(s), respectively.
    ! The purpose of bookkeeping this information is to boost the visualization
    ! of the Silo-HDF5 database file(s) in VisIt.
    real(wp), allocatable, dimension(:, :) :: spatial_extents
    real(wp), allocatable, dimension(:, :) :: data_extents

    ! The size of the ghost zone layer at beginning of each coordinate direction
    ! (lo) and at end of each coordinate direction (hi). Adding this information
    ! to Silo-HDF5 database file(s) is recommended since it supplies VisIt with
    ! connectivity information between the sub-domains of a parallel data set.
    integer, allocatable, dimension(:) :: lo_offset
    integer, allocatable, dimension(:) :: hi_offset

    ! For Silo-HDF5 database format, this variable is used to keep track of the
    ! number of cell-boundaries, for the grid associated with the local process,
    ! in each of the active coordinate directions.
    integer, allocatable, dimension(:) :: dims

    ! Locations of various folders in the case's directory tree, associated with
    ! the choice of the formatted database format. These include, in order, the
    ! location of the folder named after the selected formatted database format,
    ! and the locations of two sub-directories of the latter, the first of which
    ! is named after the local processor rank, while the second is named 'root'.
    ! The folder associated with the local processor rank contains only the data
    ! pertaining to the part of the domain taken care of by the local processor.
    ! The root directory, on the other hand, will contain either the information
    ! about the connectivity required to put the entire domain back together, or
    ! the actual data associated with the entire computational domain. This all
    ! depends on dimensionality and the choice of the formatted database format.
    character(LEN=path_len + name_len) :: dbdir
    character(LEN=path_len + 2*name_len) :: proc_rank_dir
    character(LEN=path_len + 2*name_len) :: rootdir

    ! Handles of the formatted database master/root file, slave/local processor
    ! file and options list. The list of options is explicitly used in the Silo-
    ! HDF5 database format to provide additional details about the contents of a
    ! formatted database file, such as the previously described spatial and data
    ! extents.
    integer :: dbroot
    integer :: dbfile
    integer :: optlist

    ! The total number of flow variable(s) to be stored in a formatted database
    ! file. Note that this is only needed when using the Binary format.
    integer :: dbvars

    ! Generic error flags utilized in the handling, checking and the reporting
    ! of the input and output operations errors with a formatted database file
    integer, private :: err, ierr

contains

    subroutine s_initialize_data_output_module()
        ! Description: Computation of parameters, allocation procedures, and/or
        !              any other tasks needed to properly setup the module

        ! Generic string used to store the location of a particular file
        character(LEN=len_trim(case_dir) + 2*name_len) :: file_loc

        ! Generic logical used to test the existence of a particular folder
        logical :: dir_check

        integer :: i

        ! Allocating the generic storage for the flow variable(s) that are
        ! going to be written to the formatted database file(s). Note once
        ! more that the root variable is only required for 1D computations.
        allocate (q_sf(-offset_x%beg:m + offset_x%end, &
                       -offset_y%beg:n + offset_y%end, &
                       -offset_z%beg:p + offset_z%end))
        if (grid_geometry == 3) then
            allocate (cyl_q_sf(-offset_y%beg:n + offset_y%end, &
                               -offset_z%beg:p + offset_z%end, &
                               -offset_x%beg:m + offset_x%end))
        end if

        if (precision == 1) then
            allocate (q_sf_s(-offset_x%beg:m + offset_x%end, &
                             -offset_y%beg:n + offset_y%end, &
                             -offset_z%beg:p + offset_z%end))
            if (grid_geometry == 3) then
                allocate (cyl_q_sf_s(-offset_y%beg:n + offset_y%end, &
                                     -offset_z%beg:p + offset_z%end, &
                                     -offset_x%beg:m + offset_x%end))
            end if
        end if

        if (n == 0) then
            allocate (q_root_sf(0:m_root, 0:0, 0:0))
            if (precision == 1) then
                allocate (q_root_sf_s(0:m_root, 0:0, 0:0))
            end if
        end if

        ! Allocating the spatial and data extents and also the variables for
        ! the offsets and the one bookkeeping the number of cell-boundaries
        ! in each active coordinate direction. Note that all these variables
        ! are only needed by the Silo-HDF5 format for multidimensional data.
        if (format == 1 .and. n > 0) then

            allocate (data_extents(1:2, 0:num_procs - 1))

            if (p > 0) then
                allocate (spatial_extents(1:6, 0:num_procs - 1))
                allocate (lo_offset(1:3))
                allocate (hi_offset(1:3))
                allocate (dims(1:3))
            else
                allocate (spatial_extents(1:4, 0:num_procs - 1))
                allocate (lo_offset(1:2))
                allocate (hi_offset(1:2))
                allocate (dims(1:2))
            end if

        end if

        ! The size of the ghost zone layer in each of the active coordinate
        ! directions was set in the module m_mpi_proxy.f90. The results are
        ! now transferred to the local variables of this module when they are
        ! required by the Silo-HDF5 format, for multidimensional data sets.
        ! With the same, latter, requirements, the variables bookkeeping the
        ! number of cell-boundaries in each active coordinate direction are
        ! also set here.
        if (format == 1 .and. n > 0) then
            if (p > 0) then
                if (grid_geometry == 3) then
                    lo_offset = (/offset_y%beg, offset_z%beg, offset_x%beg/)
                    hi_offset = (/offset_y%end, offset_z%end, offset_x%end/)
                else
                    lo_offset = (/offset_x%beg, offset_y%beg, offset_z%beg/)
                    hi_offset = (/offset_x%end, offset_y%end, offset_z%end/)
                end if

                if (grid_geometry == 3) then
                    dims = (/n + offset_y%beg + offset_y%end + 2, &
                             p + offset_z%beg + offset_z%end + 2, &
                             m + offset_x%beg + offset_x%end + 2/)
                else
                    dims = (/m + offset_x%beg + offset_x%end + 2, &
                             n + offset_y%beg + offset_y%end + 2, &
                             p + offset_z%beg + offset_z%end + 2/)
                end if
            else
                lo_offset = (/offset_x%beg, offset_y%beg/)
                hi_offset = (/offset_x%end, offset_y%end/)

                dims = (/m + offset_x%beg + offset_x%end + 2, &
                         n + offset_y%beg + offset_y%end + 2/)
            end if
        end if

        ! Generating Silo-HDF5 Directory Tree
        if (format == 1) then

            ! Creating the directory associated with the local process
            dbdir = trim(case_dir)//'/silo_hdf5'

            write (proc_rank_dir, '(A,I0)') '/p', proc_rank

            proc_rank_dir = trim(dbdir)//trim(proc_rank_dir)

            file_loc = trim(proc_rank_dir)//'/.'

            !INQUIRE( DIRECTORY = TRIM(file_loc), & ! Intel compiler
            !EXIST     = dir_check       )
            ! INQUIRE( FILE      = TRIM(file_loc), & ! NAG/PGI/GCC compiler
            !           EXIST     = dir_check       )
            call my_inquire(file_loc, dir_check)
            if (dir_check .neqv. .true.) then
                call s_create_directory(trim(proc_rank_dir))
            end if

            ! Creating the directory associated with the root process
            if (proc_rank == 0) then

                rootdir = trim(dbdir)//'/root'

                file_loc = trim(rootdir)//'/.'

                !INQUIRE( DIRECTORY = TRIM(file_loc), & ! Intel compiler
                !        EXIST     = dir_check       )
                !  INQUIRE( FILE      = TRIM(file_loc), & ! NAG/PGI/GCC compiler
                !           EXIST     = dir_check       )
                call my_inquire(file_loc, dir_check)
                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(rootdir))
                end if

            end if

            ! Generating Binary Directory Tree

        else

            ! Creating the directory associated with the local process
            dbdir = trim(case_dir)//'/binary'

            write (proc_rank_dir, '(A,I0)') '/p', proc_rank

            proc_rank_dir = trim(dbdir)//trim(proc_rank_dir)

            file_loc = trim(proc_rank_dir)//'/.'

            !INQUIRE( DIRECTORY = TRIM(file_loc), & ! Intel compiler
            !       EXIST     = dir_check       )
            !  INQUIRE( FILE      = TRIM(file_loc), & ! NAG/PGI/GCC compiler
            !           EXIST     = dir_check       )
            call my_inquire(file_loc, dir_check)

            if (dir_check .neqv. .true.) then
                call s_create_directory(trim(proc_rank_dir))
            end if

            ! Creating the directory associated with the root process
            if (n == 0 .and. proc_rank == 0) then

                rootdir = trim(dbdir)//'/root'

                file_loc = trim(rootdir)//'/.'

                !INQUIRE( DIRECTORY = TRIM(file_loc), & ! Intel compiler
                !        EXIST     = dir_check       )
                !  INQUIRE( FILE      = TRIM(file_loc), & ! NAG/PGI/GCC compiler
                !        EXIST     = dir_check       )
                call my_inquire(file_loc, dir_check)

                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(rootdir))
                end if

            end if

        end if

        if (bubbles_lagrange) then !Lagrangian solver
            dbdir = trim(case_dir)//'/lag_bubbles_post_process'
            file_loc = trim(dbdir)//'/.'
            call my_inquire(file_loc, dir_check)

            if (dir_check .neqv. .true.) then
                call s_create_directory(trim(dbdir))
            end if
        end if

        ! Contrary to the Silo-HDF5 database format, handles of the Binary
        ! database master/root and slave/local process files are perfectly
        ! static throughout post-process. Hence, they are set here so that
        ! they do not have to be repetitively computed in later procedures.
        if (format == 2) then
            if (n == 0 .and. proc_rank == 0) dbroot = 2
            dbfile = 1
        end if

        ! Querying Number of Flow Variable(s) in Binary Output

        if (format == 2) then

            ! Initializing the counter of the number of flow variable(s) to
            ! be written to the formatted database file(s)
            dbvars = 0

            ! Partial densities
            if ((model_eqns == 2) .or. (model_eqns == 3)) then
                do i = 1, num_fluids
                    if (alpha_rho_wrt(i) &
                        .or. &
                        (cons_vars_wrt .or. prim_vars_wrt)) then
                        dbvars = dbvars + 1
                    end if
                end do
            end if

            ! Density
            if (rho_wrt &
                .or. &
                (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) &
                then
                dbvars = dbvars + 1
            end if

            ! Momentum
            do i = 1, E_idx - mom_idx%beg
                if (mom_wrt(i) .or. cons_vars_wrt) dbvars = dbvars + 1
            end do

            ! Velocity
            do i = 1, E_idx - mom_idx%beg
                if (vel_wrt(i) .or. prim_vars_wrt) dbvars = dbvars + 1
            end do

            ! Flux limiter function
            do i = 1, E_idx - mom_idx%beg
                if (flux_wrt(i)) dbvars = dbvars + 1
            end do

            ! Energy
            if (E_wrt .or. cons_vars_wrt) dbvars = dbvars + 1

            ! Pressure
            if (pres_wrt .or. prim_vars_wrt) dbvars = dbvars + 1

            ! Elastic stresses
            if (hypoelasticity) dbvars = dbvars + (num_dims*(num_dims + 1))/2

            ! Volume fraction(s)
            if ((model_eqns == 2) .or. (model_eqns == 3)) then

                do i = 1, num_fluids - 1
                    if (alpha_wrt(i) &
                        .or. &
                        (cons_vars_wrt .or. prim_vars_wrt)) then
                        dbvars = dbvars + 1
                    end if
                end do

                if (alpha_wrt(num_fluids) &
                    .or. &
                    (cons_vars_wrt .or. prim_vars_wrt)) &
                    then
                    dbvars = dbvars + 1
                end if

            end if

            ! Specific heat ratio function
            if (gamma_wrt &
                .or. &
                (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) &
                then
                dbvars = dbvars + 1
            end if

            ! Specific heat ratio
            if (heat_ratio_wrt) dbvars = dbvars + 1

            ! Liquid stiffness function
            if (pi_inf_wrt &
                .or. &
                (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) &
                then
                dbvars = dbvars + 1
            end if

            ! Liquid stiffness
            if (pres_inf_wrt) dbvars = dbvars + 1

            ! Speed of sound
            if (c_wrt) dbvars = dbvars + 1

            ! Vorticity
            if (p > 0) then
                do i = 1, E_idx - mom_idx%beg
                    if (omega_wrt(i)) dbvars = dbvars + 1
                end do
            elseif (n > 0) then
                do i = 1, E_idx - cont_idx%end
                    if (omega_wrt(i)) dbvars = dbvars + 1
                end do
            end if

            ! Numerical Schlieren function
            if (schlieren_wrt) dbvars = dbvars + 1

        end if

        ! END: Querying Number of Flow Variable(s) in Binary Output

    end subroutine s_initialize_data_output_module

    subroutine s_define_output_region

        integer :: i
        integer :: lower_bound, upper_bound

        #:for X, M in [('x', 'm'), ('y', 'n'), ('z', 'p')]

            if (${M}$ == 0) return ! Early return for y or z if simulation is 1D or 2D

            lower_bound = -offset_${X}$%beg
            upper_bound = ${M}$+offset_${X}$%end

            do i = lower_bound, upper_bound
                if (${X}$_cc(i) > ${X}$_output%beg) then
                    ${X}$_output_idx%beg = i + offset_${X}$%beg
                    exit
                end if
            end do

            do i = upper_bound, lower_bound, -1
                if (${X}$_cc(i) < ${X}$_output%end) then
                    ${X}$_output_idx%end = i + offset_${X}$%beg
                    exit
                end if
            end do

            ! If no grid points are within the output region
            if ((${X}$_cc(lower_bound) > ${X}$_output%end) .or. (${X}$_cc(upper_bound) < ${X}$_output%beg)) then
                ${X}$_output_idx%beg = 0
                ${X}$_output_idx%end = 0
            end if

        #:endfor

    end subroutine s_define_output_region

    subroutine s_open_formatted_database_file(t_step)
        ! Description: This subroutine opens a new formatted database file, or
        !              replaces an old one, and readies it for the data storage
        !              of the grid and the flow variable(s) associated with the
        !              current time-step, t_step. This is performed by all the
        !              local process(es). The root processor, in addition, must
        !              also generate a master formatted database file whose job
        !              will be to link, and thus combine, the data from all of
        !              the local process(es). Note that for the Binary format,
        !              this extra task that is assigned to the root process is
        !              not performed in multidimensions.

        ! Time-step that is currently being post-processed
        integer, intent(IN) :: t_step

        ! Generic string used to store the location of a particular file
        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc

        ! Silo-HDF5 Database Format

        if (format == 1) then

            ! Generating the relative path to the formatted database slave
            ! file, that is to be opened for the current time-step, t_step
            write (file_loc, '(A,I0,A)') '/', t_step, '.silo'
            file_loc = trim(proc_rank_dir)//trim(file_loc)

            ! Creating formatted database slave file at the above location
            ! and setting up the structure of the file and its header info
            ierr = DBCREATE(trim(file_loc), len_trim(file_loc), &
                            DB_CLOBBER, DB_LOCAL, 'MFC v3.0', 8, &
                            DB_HDF5, dbfile)

            ! Verifying that the creation and setup process of the formatted
            ! database slave file has been performed without errors. If this
            ! is not the case, the post-process exits.
            if (dbfile == -1) then
                call s_mpi_abort('Unable to create Silo-HDF5 database '// &
                                 'slave file '//trim(file_loc)//'. '// &
                                 'Exiting ...')
            end if

            ! Next, analogous steps to the ones above are carried out by the
            ! root process to create and setup the formatted database master
            ! file.
            if (proc_rank == 0) then

                write (file_loc, '(A,I0,A)') '/collection_', t_step, '.silo'
                file_loc = trim(rootdir)//trim(file_loc)

                ierr = DBCREATE(trim(file_loc), len_trim(file_loc), &
                                DB_CLOBBER, DB_LOCAL, 'MFC v3.0', 8, &
                                DB_HDF5, dbroot)

                if (dbroot == -1) then
                    call s_mpi_abort('Unable to create Silo-HDF5 database '// &
                                     'master file '//trim(file_loc)//'. '// &
                                     'Exiting ...')
                end if

            end if

            ! Binary Database Format

        else

            ! Generating the relative path to the formatted database slave
            ! file, that is to be opened for the current time-step, t_step
            write (file_loc, '(A,I0,A)') '/', t_step, '.dat'
            file_loc = trim(proc_rank_dir)//trim(file_loc)

            ! Creating the formatted database slave file, at the previously
            ! precised relative path location, and setting up its structure
            open (dbfile, IOSTAT=err, FILE=trim(file_loc), &
                  FORM='unformatted', STATUS='replace')

            ! Verifying that the creation and setup process of the formatted
            ! database slave file has been performed without errors. If this
            ! is not the case, the post-process exits.
            if (err /= 0) then
                call s_mpi_abort('Unable to create Binary database slave '// &
                                 'file '//trim(file_loc)//'. Exiting ...')
            end if

            ! Further defining the structure of the formatted database slave
            ! file by describing in it the dimensionality of post-processed
            ! data as well as the total number of flow variable(s) that will
            ! eventually be stored in it
            if (output_partial_domain) then
                write (dbfile) x_output_idx%end - x_output_idx%beg, &
                    y_output_idx%end - y_output_idx%beg, &
                    z_output_idx%end - z_output_idx%beg, &
                    dbvars
            else
                write (dbfile) m, n, p, dbvars
            end if

            ! Next, analogous steps to the ones above are carried out by the
            ! root process to create and setup the formatted database master
            ! file. Note that this is only done in multidimensional cases.
            if (n == 0 .and. proc_rank == 0) then

                write (file_loc, '(A,I0,A)') '/', t_step, '.dat'
                file_loc = trim(rootdir)//trim(file_loc)

                open (dbroot, IOSTAT=err, FILE=trim(file_loc), &
                      FORM='unformatted', STATUS='replace')

                if (err /= 0) then
                    call s_mpi_abort('Unable to create Binary database '// &
                                     'master file '//trim(file_loc)// &
                                     '. Exiting ...')
                end if

                if (output_partial_domain) then
                    write (dbroot) x_output_idx%end - x_output_idx%beg, 0, 0, dbvars
                else
                    write (dbroot) m_root, 0, 0, dbvars
                end if

            end if

        end if

    end subroutine s_open_formatted_database_file

    subroutine s_open_intf_data_file()

        character(LEN=path_len + 3*name_len) :: file_path !<
              !! Relative path to a file in the case directory

        write (file_path, '(A)') '/intf_data.dat'
        file_path = trim(case_dir)//trim(file_path)

        ! Opening the simulation data file
        open (211, FILE=trim(file_path), &
              FORM='formatted', &
              POSITION='append', &
              STATUS='unknown')

    end subroutine s_open_intf_data_file

    subroutine s_open_energy_data_file()

        character(LEN=path_len + 3*name_len) :: file_path !<
              !! Relative path to a file in the case directory

        write (file_path, '(A)') '/eng_data.dat'
        file_path = trim(case_dir)//trim(file_path)

        ! Opening the simulation data file
        open (251, FILE=trim(file_path), &
              FORM='formatted', &
              POSITION='append', &
              STATUS='unknown')

    end subroutine s_open_energy_data_file

    subroutine s_write_grid_to_formatted_database_file(t_step)
        ! Description: The general objective of this subroutine is to write the
        !              necessary grid data to the formatted database file, for
        !              the current time-step, t_step. The local processor will
        !              write the grid data of the domain segment that it is in
        !              charge of to the formatted database slave file. The root
        !              process will additionally take care of linking that grid
        !              data in the formatted database master file. In the Silo-
        !              HDF5 database format, the spatial extents of each local
        !              process grid are also written to the master file. In the
        !              Binary format, note that no master file is maintained in
        !              multidimensions. Finally, in 1D, no grid data is written
        !              within this subroutine for the Silo-HDF5 format because
        !              curve objects rather than quadrilateral meshes are used.
        !              For curve objects, in contrast to the quadrilateral mesh
        !              objects, the grid data is included side by side with the
        !              flow variable data. Then, in this case, we take care of
        !              writing both the grid and the flow variable data in the
        !              subroutine s_write_variable_to_formatted_database_file.

        ! Time-step that is currently being post-processed
        integer, intent(IN) :: t_step

        ! Bookkeeping variables storing the name and type of mesh that is
        ! handled by the local processor(s). Note that due to an internal
        ! NAG Fortran compiler problem, these two variables could not be
        ! allocated dynamically.
        character(LEN=4*name_len), dimension(num_procs) :: meshnames
        integer, dimension(num_procs) :: meshtypes

        ! Generic loop iterator
        integer :: i

        ! Silo-HDF5 Database Format

        if (format == 1 .and. n > 0) then

            ! For multidimensional data sets, the spatial extents of all of
            ! the grid(s) handled by the local processor(s) are recorded so
            ! that they may be written, by root processor, to the formatted
            ! database master file.
            if (num_procs > 1) then
                call s_mpi_gather_spatial_extents(spatial_extents)

            elseif (p > 0) then
                if (grid_geometry == 3) then
                    spatial_extents(:, 0) = (/minval(y_cb), minval(z_cb), &
                                              minval(x_cb), maxval(y_cb), &
                                              maxval(z_cb), maxval(x_cb)/)
                else
                    spatial_extents(:, 0) = (/minval(x_cb), minval(y_cb), &
                                              minval(z_cb), maxval(x_cb), &
                                              maxval(y_cb), maxval(z_cb)/)
                end if

            else
                spatial_extents(:, 0) = (/minval(x_cb), minval(y_cb), &
                                          maxval(x_cb), maxval(y_cb)/)

            end if

            ! Next, the root processor proceeds to record all of the spatial
            ! extents in the formatted database master file. In addition, it
            ! also records a sub-domain connectivity map so that the entire
            ! grid may be reassembled by looking at the master file.
            if (proc_rank == 0) then

                do i = 1, num_procs
                    write (meshnames(i), '(A,I0,A,I0,A)') '../p', i - 1, &
                        '/', t_step, '.silo:rectilinear_grid'
                end do

                meshtypes = DB_QUAD_RECT

                err = DBSET2DSTRLEN(len(meshnames(1)))
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIOPT(optlist, DBOPT_EXTENTS_SIZE, &
                                size(spatial_extents, 1))
                err = DBADDDOPT(optlist, DBOPT_EXTENTS, spatial_extents)
                err = DBPUTMMESH(dbroot, 'rectilinear_grid', 16, &
                                 num_procs, meshnames, &
                                 len_trim(meshnames), &
                                 meshtypes, optlist, ierr)
                err = DBFREEOPTLIST(optlist)

            end if

            ! Finally, the local quadrilateral mesh, either 2D or 3D, along
            ! with its offsets that indicate the presence and size of ghost
            ! zone layer(s), are put in the formatted database slave file.

            if (precision == 1) then
                if (p > 0) then
                    z_cb_s = real(z_cb, sp)
                end if
                x_cb_s = real(x_cb, sp)
                y_cb_s = real(y_cb, sp)
            end if

            #:for PRECISION, SFX, DBT in [(1,'_s','DB_FLOAT'),(2,'',"DB_DOUBLE")]
                if (precision == ${PRECISION}$) then
                    if (p > 0) then
                        err = DBMKOPTLIST(2, optlist)
                        err = DBADDIOPT(optlist, DBOPT_LO_OFFSET, lo_offset)
                        err = DBADDIOPT(optlist, DBOPT_HI_OFFSET, hi_offset)
                        if (grid_geometry == 3) then
                            err = DBPUTQM(dbfile, 'rectilinear_grid', 16, &
                                          'x', 1, 'y', 1, 'z', 1, &
                                          y_cb${SFX}$, z_cb${SFX}$, x_cb${SFX}$, dims, 3, &
                                          ${DBT}$, DB_COLLINEAR, &
                                          optlist, ierr)
                        else
                            err = DBPUTQM(dbfile, 'rectilinear_grid', 16, &
                                          'x', 1, 'y', 1, 'z', 1, &
                                          x_cb${SFX}$, y_cb${SFX}$, z_cb${SFX}$, dims, 3, &
                                          ${DBT}$, DB_COLLINEAR, &
                                          optlist, ierr)
                        end if
                        err = DBFREEOPTLIST(optlist)
                    else
                        err = DBMKOPTLIST(2, optlist)
                        err = DBADDIOPT(optlist, DBOPT_LO_OFFSET, lo_offset)
                        err = DBADDIOPT(optlist, DBOPT_HI_OFFSET, hi_offset)
                        err = DBPUTQM(dbfile, 'rectilinear_grid', 16, &
                                      'x', 1, 'y', 1, 'z', 1, &
                                      x_cb${SFX}$, y_cb${SFX}$, DB_F77NULL, dims, 2, &
                                      ${DBT}$, DB_COLLINEAR, &
                                      optlist, ierr)
                        err = DBFREEOPTLIST(optlist)
                    end if
                end if
            #:endfor

            ! END: Silo-HDF5 Database Format

            ! Binary Database Format

        elseif (format == 2) then

            ! Multidimensional local grid data is written to the formatted
            ! database slave file. Recall that no master file to maintained
            ! in multidimensions.
            if (p > 0) then
                if (precision == 1) then
                    write (dbfile) real(x_cb, sp), &
                        real(y_cb, sp), &
                        real(z_cb, sp)
                else
                    if (output_partial_domain) then
                        write (dbfile) x_cb(x_output_idx%beg - 1:x_output_idx%end), &
                            y_cb(y_output_idx%beg - 1:y_output_idx%end), &
                            z_cb(z_output_idx%beg - 1:z_output_idx%end)
                    else
                        write (dbfile) x_cb, y_cb, z_cb
                    end if
                end if

            elseif (n > 0) then
                if (precision == 1) then
                    write (dbfile) real(x_cb, sp), &
                        real(y_cb, sp)
                else
                    if (output_partial_domain) then
                        write (dbfile) x_cb(x_output_idx%beg - 1:x_output_idx%end), &
                            y_cb(y_output_idx%beg - 1:y_output_idx%end)
                    else
                        write (dbfile) x_cb, y_cb
                    end if
                end if

                ! One-dimensional local grid data is written to the formatted
                ! database slave file. In addition, the local grid data is put
                ! together by the root process and written to the master file.
            else

                if (precision == 1) then
                    write (dbfile) real(x_cb, sp)
                else
                    if (output_partial_domain) then
                        write (dbfile) x_cb(x_output_idx%beg - 1:x_output_idx%end)
                    else
                        write (dbfile) x_cb
                    end if
                end if

                if (num_procs > 1) then
                    call s_mpi_defragment_1d_grid_variable()
                else
                    x_root_cb = x_cb
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

    subroutine s_write_variable_to_formatted_database_file(varname, t_step)
        ! Description: The goal of this subroutine is to write to the formatted
        !              database file the flow variable at the current time-step,
        !              t_step. The local process(es) write the part of the flow
        !              variable that they handle to the formatted database slave
        !              file. The root process, on the other hand, will also take
        !              care of connecting all of the flow variable data in the
        !              formatted database master file. In the Silo-HDF5 database
        !              format, the extents of each local process flow variable
        !              are also written to the master file. Note that in Binary
        !              format, no master file is maintained in multidimensions.
        !              Finally note that in 1D, grid data is also written within
        !              this subroutine for Silo-HDF5 database format since curve
        !              and not the quadrilateral variable objects are used, see
        !              description of s_write_grid_to_formatted_database_file
        !              for more details on this topic.

        ! Name of the flow variable, which will be written to the formatted
        ! database file at the current time-step, t_step
        character(LEN=*), intent(IN) :: varname

        ! Time-step that is currently being post-processed
        integer, intent(IN) :: t_step

        ! Bookkeeping variables storing the name and type of flow variable
        ! that is about to be handled by the local processor(s). Note that
        ! due to an internal NAG Fortran compiler problem, these variables
        ! could not be allocated dynamically.
        character(LEN=4*name_len), dimension(num_procs) :: varnames
        integer, dimension(num_procs) :: vartypes

        ! Generic loop iterator
        integer :: i, j, k
        real(wp) :: start, finish

        ! Silo-HDF5 Database Format

        if (format == 1) then

            ! In 1D, a curve object, featuring the local processor grid and
            ! flow variable data, is written to the formatted database slave
            ! file. The root process, on the other hand, will also take care
            ! of gathering the entire grid and associated flow variable data
            ! and write it to the formatted database master file.
            if (n == 0) then

                if (precision == 1 .and. wp == dp) then
                    x_cc_s = real(x_cc, sp)
                    q_sf_s = real(q_sf, sp)
                elseif (precision == 1 .and. wp == sp) then
                    x_cc_s = x_cc
                    q_sf_s = q_sf
                end if

                ! Writing the curve object associated with the local process
                ! to the formatted database slave file
                #:for PRECISION, SFX, DBT in [(1,'_s','DB_FLOAT'),(2,'',"DB_DOUBLE")]
                    if (precision == ${PRECISION}$) then
                        err = DBPUTCURVE(dbfile, trim(varname), len_trim(varname), &
                                         x_cc${SFX}$ (0:m), q_sf${SFX}$, ${DBT}$, m + 1, &
                                         DB_F77NULL, ierr)
                    end if
                #:endfor

                ! Assembling the local grid and flow variable data for the
                ! entire computational domain on to the root process

                if (num_procs > 1) then
                    call s_mpi_defragment_1d_grid_variable()
                    call s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)

                    if (precision == 1) then
                        x_root_cc_s = real(x_root_cc, sp)
                        q_root_sf_s = real(q_root_sf, sp)
                    end if
                else
                    if (precision == 1) then
                        x_root_cc_s = real(x_cc, sp)
                        q_root_sf_s = real(q_sf, sp)
                    else
                        x_root_cc = x_cc
                        q_root_sf = q_sf
                    end if
                end if

                ! Writing the curve object associated with the root process
                ! to the formatted database master file
                if (proc_rank == 0) then
                    #:for PRECISION, SFX, DBT in [(1,'_s','DB_FLOAT'),(2,'',"DB_DOUBLE")]
                        if (precision == ${PRECISION}$) then
                            err = DBPUTCURVE(dbroot, trim(varname), &
                                             len_trim(varname), &
                                             x_root_cc${SFX}$, q_root_sf${SFX}$, &
                                             ${DBT}$, m_root + 1, &
                                             DB_F77NULL, ierr)
                        end if
                    #:endfor
                end if

                return

                ! In multidimensions, the local process(es) take care of writing
                ! the flow variable data they are in charge of to the formatted
                ! database slave file. The root processor, additionally, is also
                ! responsible in gathering the flow variable extents of each of
                ! the local processor(s) and writing them to formatted database
                ! master file.
            else

                ! Determining the extents of the flow variable on each local
                ! process and gathering all this information on root process
                if (num_procs > 1) then
                    call s_mpi_gather_data_extents(q_sf, data_extents)
                else
                    data_extents(:, 0) = (/minval(q_sf), maxval(q_sf)/)
                end if

                ! Next, the root process proceeds to write the gathered flow
                ! variable data extents to formatted database master file.
                if (proc_rank == 0) then

                    do i = 1, num_procs
                        write (varnames(i), '(A,I0,A,I0,A)') '../p', i - 1, &
                            '/', t_step, '.silo:'//trim(varname)
                    end do

                    vartypes = DB_QUADVAR

                    err = DBSET2DSTRLEN(len(varnames(1)))
                    err = DBMKOPTLIST(2, optlist)
                    err = DBADDIOPT(optlist, DBOPT_EXTENTS_SIZE, 2)
                    err = DBADDDOPT(optlist, DBOPT_EXTENTS, data_extents)
                    err = DBPUTMVAR(dbroot, trim(varname), &
                                    len_trim(varname), num_procs, &
                                    varnames, len_trim(varnames), &
                                    vartypes, optlist, ierr)
                    err = DBFREEOPTLIST(optlist)

                end if

                ! Finally, each of the local processor(s) proceeds to write
                ! the flow variable data that it is responsible for to the
                ! formatted database slave file.
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
                elseif (wp == dp) then
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
                                err = DBPUTQV1(dbfile, trim(varname), &
                                               len_trim(varname), &
                                               'rectilinear_grid', 16, &
                                               cyl_q_sf${SFX}$, dims - 1, 3, DB_F77NULL, &
                                               0, ${DBT}$, DB_ZONECENT, &
                                               DB_F77NULL, ierr)
                            else
                                err = DBPUTQV1(dbfile, trim(varname), &
                                               len_trim(varname), &
                                               'rectilinear_grid', 16, &
                                               q_sf${SFX}$, dims - 1, 3, DB_F77NULL, &
                                               0, ${DBT}$, DB_ZONECENT, &
                                               DB_F77NULL, ierr)
                            end if
                        else
                            err = DBPUTQV1(dbfile, trim(varname), &
                                           len_trim(varname), &
                                           'rectilinear_grid', 16, &
                                           q_sf${SFX}$, dims - 1, 2, DB_F77NULL, &
                                           0, ${DBT}$, DB_ZONECENT, &
                                           DB_F77NULL, ierr)
                        end if
                    end if
                #:endfor

            end if

            ! END: Silo-HDF5 Database Format

            ! Binary Database Format

        else

            ! Writing the name of the flow variable and its data, associated
            ! with the local processor, to the formatted database slave file
            if (precision == 1) then
                write (dbfile) varname, real(q_sf, wp)
            else
                write (dbfile) varname, q_sf
            end if

            ! In 1D, the root process also takes care of gathering the flow
            ! variable data from all of the local processor(s) and writes it
            ! to the formatted database master file.
            if (n == 0) then

                if (num_procs > 1) then
                    call s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)
                else
                    q_root_sf = q_sf
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

    !>  Subroutine that writes the post processed results in the folder 'lag_bubbles_data'
            !!  @param t_step Current time step
    subroutine s_write_lag_bubbles_results(t_step)

        integer, intent(in) :: t_step
        character(len=len_trim(case_dir) + 2*name_len) :: t_step_dir
        character(len=len_trim(case_dir) + 3*name_len) :: file_loc
        logical :: dir_check
        integer :: id, nlg_bubs

#ifdef MFC_MPI
        real(wp), dimension(20) :: inputvals
        real(wp) :: id_real, time_real
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer :: view

        integer, dimension(3) :: cell
        logical :: indomain, lg_bub_file, lg_bub_data, file_exist

        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer :: ifile, ireq, ierr, data_size, tot_data
        integer :: i

        write (file_loc, '(A,I0,A)') 'lag_bubbles_mpi_io_', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            if (proc_rank == 0) then
                open (9, FILE=trim(file_loc), FORM='unformatted', STATUS='unknown')
                read (9) tot_data, time_real
                close (9)
            end if
        else
            print '(A)', trim(file_loc)//' is missing. Exiting ...'
            call s_mpi_abort
        end if

        call MPI_BCAST(tot_data, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(time_real, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

        gsizes(1) = tot_data
        gsizes(2) = 21
        lsizes(1) = tot_data
        lsizes(2) = 21
        start_idx_part(1) = 0
        start_idx_part(2) = 0

        call MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                      MPI_ORDER_FORTRAN, mpi_p, view, ierr)
        call MPI_TYPE_COMMIT(view, ierr)

        write (file_loc, '(A,I0,A)') 'lag_bubbles_', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=lg_bub_file)

        if (lg_bub_file) then

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)

            disp = 0._wp
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, &
                                   'native', mpi_info_null, ierr)

            allocate (MPI_IO_DATA_lg_bubbles(tot_data, 1:21))

            call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA_lg_bubbles, 21*tot_data, &
                                   mpi_p, status, ierr)

            write (file_loc, '(A,I0,A)') 'lag_bubbles_post_process_', t_step, '.dat'
            file_loc = trim(case_dir)//'/lag_bubbles_post_process/'//trim(file_loc)

            if (proc_rank == 0) then
                open (unit=29, file=file_loc, form='formatted', position='rewind')
                !write(29,*) 'lg_bubID, x, y, z, xPrev, yPrev, zPrev, xVel, yVel, ',   &
                !            'zVel, radius, interfaceVelocity, equilibriumRadius',       &
                !            'Rmax, Rmin, dphidt, pressure, mv, mg, betaT, betaC, time'
                do i = 1, tot_data
                    id = int(MPI_IO_DATA_lg_bubbles(i, 1))
                    inputvals(1:20) = MPI_IO_DATA_lg_bubbles(i, 2:21)
                    if (id > 0) then
                        write (29, 6) int(id), inputvals(1), inputvals(2), &
                            inputvals(3), inputvals(4), inputvals(5), inputvals(6), inputvals(7), &
                            inputvals(8), inputvals(9), inputvals(10), inputvals(11), &
                            inputvals(12), inputvals(13), inputvals(14), inputvals(15), &
                            inputvals(16), inputvals(17), inputvals(18), inputvals(19), &
                            inputvals(20), time_real
6                       format(I6, 21(1x, E15.7))
                    end if
                end do
                close (29)
            end if

            deallocate (MPI_IO_DATA_lg_bubbles)

        end if

        call s_mpi_barrier()

        call MPI_FILE_CLOSE(ifile, ierr)

#endif

    end subroutine s_write_lag_bubbles_results
    subroutine s_write_intf_data_file(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        integer :: i, j, k, l, w, cent !< Generic loop iterators
        integer :: ierr, counter, root !< number of data points extracted to fit shape to SH perturbations
        real(wp), dimension(num_fluids) :: alpha, vol_fluid, xcom, ycom, zcom
        real(wp), parameter :: pi = 4._wp*tan(1._wp)
        real(wp), allocatable :: x_td(:), y_td(:), x_d1(:), y_d1(:), y_d(:), x_d(:)
        real(wp) :: axp, axm, ayp, aym, azm, azp, tgp, euc_d, thres, maxalph_loc, maxalph_glb

        allocate (x_d1(m*n))
        allocate (y_d1(m*n))
        counter = 0
        maxalph_loc = 0_wp
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    if (q_prim_vf(E_idx + 2)%sf(i, j, k) > maxalph_loc) then
                        maxalph_loc = q_prim_vf(E_idx + 2)%sf(i, j, k)
                    end if
                end do
            end do
        end do

        call s_mpi_allreduce_max(maxalph_loc, maxalph_glb)
        if (p > 0) then
            do l = 0, p
                if (z_cc(l) < dz(l) .and. z_cc(l) > 0) then
                    cent = l
                end if
            end do
        else
            cent = 0
        end if

        thres = 0.9_wp*maxalph_glb
        do k = 0, n
            do j = 0, m
                axp = q_prim_vf(E_idx + 2)%sf(j + 1, k, cent)
                axm = q_prim_vf(E_idx + 2)%sf(j, k, cent)
                ayp = q_prim_vf(E_idx + 2)%sf(j, k + 1, cent)
                aym = q_prim_vf(E_idx + 2)%sf(j, k, cent)
                if ((axp > thres .and. axm < thres) .or. (axp < thres .and. axm > thres) &
                    .or. (ayp > thres .and. aym < thres) .or. (ayp < thres .and. aym > thres)) then
                    if (counter == 0) then
                        counter = counter + 1
                        x_d1(counter) = x_cc(j)
                        y_d1(counter) = y_cc(k)
                        euc_d = sqrt((x_cc(j) - x_d1(i))**2 + (y_cc(k) - y_d1(i))**2)
                        tgp = sqrt(dx(j)**2 + dy(k)**2)
                    else
                        euc_d = sqrt((x_cc(j) - x_d1(i))**2 + (y_cc(k) - y_d1(i))**2)
                        tgp = sqrt(dx(j)**2 + dy(k)**2)
                        do i = 1, counter
                            if (euc_d < tgp) then
                                cycle
                            elseif (euc_d > tgp .and. i == counter) then
                                counter = counter + 1
                                x_d1(counter) = x_cc(j)
                                y_d1(counter) = y_cc(k)

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
                    write (211, '(F12.9,1X,F12.9,1X,I4)') &
                        x_td(i), y_td(i), size(x_td)
                else
                    write (211, '(F12.9,1X,F12.9,1X,F3.1)') &
                        x_td(i), y_td(i), 0_wp
                end if
            end do
        end if

    end subroutine s_write_intf_data_file

    subroutine s_write_energy_data_file(q_prim_vf, q_cons_vf)
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf, q_cons_vf
        real(wp) :: Elk, Egk, Elp, Egint, Vb, Vl, pres_av, Et
        real(wp) :: rho, pres, dV, tmp, gamma, pi_inf, MaxMa, MaxMa_glb, maxvel, c, Ma, H
        real(wp), dimension(num_dims) :: vel
        real(wp), dimension(num_fluids) :: gammas, pi_infs, adv
        integer :: i, j, k, l, s !looping indices
        integer :: ierr, counter, root !< number of data points extracted to fit shape to SH perturbations

        Egk = 0_wp
        Elp = 0_wp
        Egint = 0_wp
        Vb = 0_wp
        maxvel = 0_wp
        MaxMa = 0_wp
        Vl = 0_wp
        Elk = 0_wp
        Et = 0_wp
        Vb = 0_wp
        dV = 0_wp
        pres_av = 0_wp
        pres = 0_wp
        c = 0._wp

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    pres = 0_wp
                    dV = dx(i)*dy(j)*dz(k)
                    rho = 0_wp
                    gamma = 0_wp
                    pi_inf = 0_wp
                    pres = q_prim_vf(E_idx)%sf(i, j, k)
                    Egint = Egint + q_prim_vf(E_idx + 2)%sf(i, j, k)*(fluid_pp(2)%gamma*pres)*dV
                    do s = 1, num_dims
                        vel(s) = q_prim_vf(num_fluids + s)%sf(i, j, k)
                        Egk = Egk + 0.5_wp*q_prim_vf(E_idx + 2)%sf(i, j, k)*q_prim_vf(2)%sf(i, j, k)*vel(s)*vel(s)*dV
                        Elk = Elk + 0.5_wp*q_prim_vf(E_idx + 1)%sf(i, j, k)*q_prim_vf(1)%sf(i, j, k)*vel(s)*vel(s)*dV
                        if (abs(vel(s)) > maxvel) then
                            maxvel = abs(vel(s))
                        end if
                    end do
                    do l = 1, adv_idx%end - E_idx
                        adv(l) = q_prim_vf(E_idx + l)%sf(i, j, k)
                        gamma = gamma + adv(l)*fluid_pp(l)%gamma
                        pi_inf = pi_inf + adv(l)*fluid_pp(l)%pi_inf
                        rho = rho + adv(l)*q_prim_vf(l)%sf(i, j, k)
                    end do

                    H = ((gamma + 1_wp)*pres + pi_inf)/rho

                    call s_compute_speed_of_sound(pres, rho, &
                                                  gamma, pi_inf, &
                                                  H, adv, 0._wp, 0._wp, c)

                    Ma = maxvel/c
                    if (Ma > MaxMa .and. (adv(1) > (1.0_wp - 1.0e-10_wp))) then
                        MaxMa = Ma
                    end if
                    Vl = Vl + adv(1)*dV
                    Vb = Vb + adv(2)*dV
                    pres_av = pres_av + adv(1)*pres*dV
                    Et = Et + q_cons_vf(E_idx)%sf(i, j, k)*dV
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
            write (251, '(10X, 8F24.8)') &
                Elp, &
                Egint, &
                Elk, &
                Egk, &
                Et, &
                Vb, &
                Vl, &
                MaxMa_glb
        end if

    end subroutine s_write_energy_data_file

    subroutine s_close_formatted_database_file()
        ! Description: The purpose of this subroutine is to close any formatted
        !              database file(s) that may be opened at the time-step that
        !              is currently being post-processed. The root process must
        !              typically close two files, one associated with the local
        !              sub-domain and the other with the entire domain. The non-
        !              root process(es) must close one file, which is associated
        !              with the local sub-domain. Note that for the Binary data-
        !              base format and multidimensional data, the root process
        !              only has to close the file associated with the local sub-
        !              domain, because one associated with the entire domain is
        !              not generated.

        ! Silo-HDF5 database format
        if (format == 1) then
            ierr = DBCLOSE(dbfile)
            if (proc_rank == 0) ierr = DBCLOSE(dbroot)

            ! Binary database format
        else
            close (dbfile)
            if (n == 0 .and. proc_rank == 0) close (dbroot)

        end if

    end subroutine s_close_formatted_database_file

    subroutine s_close_intf_data_file()

        close (211)

    end subroutine s_close_intf_data_file

    subroutine s_close_energy_data_file()

        close (251)

    end subroutine s_close_energy_data_file

    subroutine s_finalize_data_output_module()
        ! Description: Deallocation procedures for the module

        ! Deallocating the generic storage employed for the flow variable(s)
        ! that were written to the formatted database file(s). Note that the
        ! root variable is only deallocated in the case of a 1D computation.
        deallocate (q_sf)
        if (n == 0) deallocate (q_root_sf)
        if (grid_geometry == 3) then
            deallocate (cyl_q_sf)
        end if

        ! Deallocating spatial and data extents and also the variables for
        ! the offsets and the one bookkeeping the number of cell-boundaries
        ! in each active coordinate direction. Note that all these variables
        ! were only needed by Silo-HDF5 format for multidimensional data.
        if (format == 1 .and. n > 0) then
            deallocate (spatial_extents)
            deallocate (data_extents)
            deallocate (lo_offset)
            deallocate (hi_offset)
            deallocate (dims)
        end if

    end subroutine s_finalize_data_output_module

end module m_data_output
