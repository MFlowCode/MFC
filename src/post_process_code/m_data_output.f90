!>
!! @file m_data_output.f90
!! @brief Contains module m_data_output
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module enables the restructuring of the raw simulation data
!!              file(s) into formatted database file(s). The formats that may be
!!              chosen from include Silo-HDF5 and Binary. Each of these database
!!              structures contains information about the grid as well as each of
!!              the flow variable(s) that were chosen by the user to be included.
MODULE m_data_output
    
    
    ! Dependencies =============================================================
    ! USE f90_unix_proc         ! NAG Compiler Library of UNIX system commands
    
    USE m_derived_types         ! Definitions of the derived types
    
    USE m_global_parameters     ! Global parameters for the code
    
    USE m_mpi_proxy             ! Message passing interface (MPI) module proxy

    USE m_compile_specific
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_data_output_module, &
                       s_open_formatted_database_file, &
                       s_write_grid_to_formatted_database_file, &
                       s_write_variable_to_formatted_database_file, &
                       s_close_formatted_database_file, &
                       s_finalize_data_output_module
    
    ! Including the Silo Fortran interface library that features the subroutines
    ! and parameters that are required to write in the Silo-HDF5 database format
    ! INCLUDE 'silo.inc'
    INCLUDE 'silo_f9x.inc'
    
    
    ! Generic storage for flow variable(s) that are to be written to formatted
    ! database file(s). Note that for 1D simulations, q_root_sf is employed to
    ! gather the flow variable(s) from all sub-domains on to the root process.
    ! If the run is not parallel, but serial, then q_root_sf is equal to q_sf.
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: q_sf
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: dft_q_sf
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: coarse_x_q_sf
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: coarse_xy_q_sf
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: coarse_xyz_q_sf
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: q_root_sf
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: cyl_q_sf
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: cyl_coarse_q_sf
    
    ! The spatial and data extents array variables contain information about the
    ! minimum and maximum values of the grid and flow variable(s), respectively.
    ! The purpose of bookkeeping this information is to boost the visualization
    ! of the Silo-HDF5 database file(s) in VisIt.
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: spatial_extents
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: data_extents
    
    ! The size of the ghost zone layer at beginning of each coordinate direction
    ! (lo) and at end of each coordinate direction (hi). Adding this information
    ! to Silo-HDF5 database file(s) is recommended since it supplies VisIt with
    ! connectivity information between the sub-domains of a parallel data set.
    INTEGER, ALLOCATABLE, DIMENSION(:) :: lo_offset
    INTEGER, ALLOCATABLE, DIMENSION(:) :: hi_offset
    
    ! For Silo-HDF5 database format, this variable is used to keep track of the
    ! number of cell-boundaries, for the grid associated with the local process,
    ! in each of the active coordinate directions.
    INTEGER, ALLOCATABLE, DIMENSION(:) :: dims
    
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
    CHARACTER(LEN = path_len +   name_len) :: dbdir
    CHARACTER(LEN = path_len + 2*name_len) :: proc_rank_dir
    CHARACTER(LEN = path_len + 2*name_len) :: rootdir
    
    ! Handles of the formatted database master/root file, slave/local processor
    ! file and options list. The list of options is explicitly used in the Silo-
    ! HDF5 database format to provide additional details about the contents of a
    ! formatted database file, such as the previously described spatial and data
    ! extents.
    INTEGER :: dbroot
    INTEGER :: dbfile
    INTEGER :: optlist
    
    ! The total number of flow variable(s) to be stored in a formatted database
    ! file. Note that this is only needed when using the Binary format.
    INTEGER :: dbvars
    
    ! Generic error flags utilized in the handling, checking and the reporting
    ! of the input and output operations errors with a formatted database file
    INTEGER, PRIVATE :: err, ierr
    
    
    CONTAINS
        
        
        
        
        
        SUBROUTINE s_initialize_data_output_module() ! ----------------------------
        ! Description: Computation of parameters, allocation procedures, and/or
        !              any other tasks needed to properly setup the module
            
            
            ! Generic string used to store the location of a particular file
            CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: file_loc
            
            ! Generic logical used to test the existence of a particular folder
            LOGICAL :: dir_check
            
            ! Generic loop iterator
            INTEGER :: i
            
            
            ! Allocating the generic storage for the flow variable(s) that are
            ! going to be written to the formatted database file(s). Note once
            ! more that the root variable is only required for 1D computations.
            ALLOCATE(q_sf( -offset_x%beg : m + offset_x%end, &
                           -offset_y%beg : n + offset_y%end, &
                           -offset_z%beg : p + offset_z%end ))
            IF (fourier_decomp) THEN
                ALLOCATE(dft_q_sf( -offset_x%beg : m + offset_x%end, &
                                   -offset_y%beg : n + offset_y%end, &
                                   -offset_z%beg : p + offset_z%end ))
            END IF
            IF (grid_geometry == 3) THEN
                IF (coarsen_silo) THEN
                    ALLOCATE(cyl_coarse_q_sf( -offset_y%beg : (n/2) + offset_y%end, &
                                              -offset_z%beg : (p/2) + offset_z%end, &
                                              -offset_x%beg : (m/2) + offset_x%end ))
                ELSE
                    ALLOCATE(cyl_q_sf( -offset_y%beg : n + offset_y%end, &
                                       -offset_z%beg : p + offset_z%end, &
                                       -offset_x%beg : m + offset_x%end ))
                END IF
            END IF
            IF (coarsen_silo) THEN
                ALLOCATE(coarse_x_q_sf(   -offset_x%beg : (m/2) + offset_x%end, &
                                          -offset_y%beg :   n   + offset_y%end, &
                                          -offset_z%beg :   p   + offset_z%end ))
                ALLOCATE(coarse_xy_q_sf(  -offset_x%beg : (m/2) + offset_x%end, &
                                          -offset_y%beg : (n/2) + offset_y%end, &
                                          -offset_z%beg :   p   + offset_z%end ))
                ALLOCATE(coarse_xyz_q_sf( -offset_x%beg : (m/2) + offset_x%end, &
                                          -offset_y%beg : (n/2) + offset_y%end, &
                                          -offset_z%beg : (p/2) + offset_z%end ))
            END IF
            
            IF(n == 0) ALLOCATE(q_root_sf(0:m_root,0:0,0:0))
            
            
            ! Allocating the spatial and data extents and also the variables for
            ! the offsets and the one bookkeeping the number of cell-boundaries
            ! in each active coordinate direction. Note that all these variables
            ! are only needed by the Silo-HDF5 format for multidimensional data.
            IF(format == 1 .AND. n > 0) THEN
                
                ALLOCATE(data_extents(1:2,0:num_procs-1))
                
                IF(p > 0) THEN
                    ALLOCATE(spatial_extents(1:6,0:num_procs-1))
                    ALLOCATE(lo_offset(1:3))
                    ALLOCATE(hi_offset(1:3))
                    ALLOCATE(dims(1:3))
                ELSE
                    ALLOCATE(spatial_extents(1:4,0:num_procs-1))
                    ALLOCATE(lo_offset(1:2))
                    ALLOCATE(hi_offset(1:2))
                    ALLOCATE(dims(1:2))
                END IF
                
            END IF
            
            
            ! The size of the ghost zone layer in each of the active coordinate
            ! directions was set in the module m_mpi_proxy.f90. The results are
            ! now transfered to the local variables of this module when they are
            ! required by the Silo-HDF5 format, for multidimensional data sets.
            ! With the same, latter, requirements, the variables bookkeeping the
            ! number of cell-boundaries in each active coordinate direction are
            ! also set here.
            IF(format == 1 .AND. n > 0) THEN
                IF(p > 0) THEN
                    IF (grid_geometry == 3) THEN
                        lo_offset = (/ offset_y%beg, offset_z%beg, offset_x%beg /)
                        hi_offset = (/ offset_y%end, offset_z%end, offset_x%end /)
                    ELSE
                        lo_offset = (/ offset_x%beg, offset_y%beg, offset_z%beg /)
                        hi_offset = (/ offset_x%end, offset_y%end, offset_z%end /)
                    END IF
                    IF (coarsen_silo) THEN
                        IF (grid_geometry == 3) THEN
                            dims      = (/ (n/2)+offset_y%beg + offset_y%end + 2, &
                                           (p/2)+offset_z%beg + offset_z%end + 2, &
                                           (m/2)+offset_x%beg + offset_x%end + 2 /)
                        ELSE
                            dims      = (/ (m/2)+offset_x%beg + offset_x%end + 2, &
                                           (n/2)+offset_y%beg + offset_y%end + 2, &
                                           (p/2)+offset_z%beg + offset_z%end + 2 /)
                        END IF
                    ELSE
                        IF (grid_geometry == 3) THEN
                            dims      = (/ n+offset_y%beg + offset_y%end + 2, &
                                           p+offset_z%beg + offset_z%end + 2, &
                                           m+offset_x%beg + offset_x%end + 2 /)
                        ELSE
                            dims      = (/ m+offset_x%beg + offset_x%end + 2, &
                                           n+offset_y%beg + offset_y%end + 2, &
                                           p+offset_z%beg + offset_z%end + 2 /)
                        END IF
                    END IF
                ELSE
                    lo_offset = (/ offset_x%beg, offset_y%beg /)
                    hi_offset = (/ offset_x%end, offset_y%end /)
                    IF (coarsen_silo) THEN
                        dims      = (/ (m/2)+offset_x%beg + offset_x%end + 2, &
                                       (n/2)+offset_y%beg + offset_y%end + 2 /)
                    ELSE
                        dims      = (/ m+offset_x%beg + offset_x%end + 2, &
                                       n+offset_y%beg + offset_y%end + 2 /)
                    END IF
                END IF
            END IF
            
            
            ! Generating Silo-HDF5 Directory Tree ==============================
            
            IF(format == 1) THEN
                
                ! Creating the directory associated with the local process
                IF (coarsen_silo) THEN
                    IF (fourier_decomp) THEN
                        dbdir = TRIM(case_dir) // '/coarse_fourier_silo_hdf5'
                    ELSE
                        dbdir = TRIM(case_dir) // '/coarse_silo_hdf5'
                    END IF
                ELSE
                    IF (fourier_decomp) THEN
                        dbdir = TRIM(case_dir) // '/fourier_silo_hdf5'
                    ELSE
                        dbdir = TRIM(case_dir) // '/silo_hdf5'
                    END IF
                END IF
                
                WRITE(proc_rank_dir, '(A,I0)') '/p', proc_rank
                
                proc_rank_dir = TRIM(dbdir) // TRIM(proc_rank_dir)
                
                file_loc = TRIM(proc_rank_dir) // '/.'
                
               !INQUIRE( DIRECTORY = TRIM(file_loc), & ! Intel compiler
               !EXIST     = dir_check       )
               ! INQUIRE( FILE      = TRIM(file_loc), & ! NAG/PGI/GCC compiler
               !           EXIST     = dir_check       )
                call my_inquire(file_loc,dir_check)
               IF(dir_check .NEQV. .TRUE.) THEN
                   CALL SYSTEM('mkdir -p ' // TRIM(proc_rank_dir))
               END IF
                
                ! Creating the directory associated with the root process
                IF(proc_rank == 0) THEN
                    
                    rootdir = TRIM(dbdir) // '/root'
                    
                    file_loc = TRIM(rootdir) // '/.'
                    
                   !INQUIRE( DIRECTORY = TRIM(file_loc), & ! Intel compiler
                   !        EXIST     = dir_check       )
                   !  INQUIRE( FILE      = TRIM(file_loc), & ! NAG/PGI/GCC compiler
                   !           EXIST     = dir_check       )
                    call my_inquire(file_loc,dir_check)
                   IF(dir_check .NEQV. .TRUE.) THEN
                       CALL SYSTEM('mkdir ' // TRIM(rootdir))
                   END IF
                    
                END IF
                
            ! ==================================================================
                
                
            ! Generating Binary Directory Tree =================================
                
            ELSE
                
                ! Creating the directory associated with the local process
                dbdir = TRIM(case_dir) // '/binary'
                
                WRITE(proc_rank_dir, '(A,I0)') '/p', proc_rank
                
                proc_rank_dir = TRIM(dbdir) // TRIM(proc_rank_dir)
                
                file_loc = TRIM(proc_rank_dir) // '/.'
                
                !INQUIRE( DIRECTORY = TRIM(file_loc), & ! Intel compiler
                !       EXIST     = dir_check       )
                !  INQUIRE( FILE      = TRIM(file_loc), & ! NAG/PGI/GCC compiler
                !           EXIST     = dir_check       )
                call my_inquire(file_loc,dir_check)

                IF(dir_check .NEQV. .TRUE.) THEN
                    CALL SYSTEM('mkdir -p ' // TRIM(proc_rank_dir))
                END IF
                
                ! Creating the directory associated with the root process
                IF(n == 0 .AND. proc_rank == 0) THEN
                    
                    rootdir = TRIM(dbdir) // '/root'
                    
                    file_loc = TRIM(rootdir) // '/.'
                    
                   !INQUIRE( DIRECTORY = TRIM(file_loc), & ! Intel compiler
                   !        EXIST     = dir_check       )
                   !  INQUIRE( FILE      = TRIM(file_loc), & ! NAG/PGI/GCC compiler
                   !        EXIST     = dir_check       )
                   call my_inquire(file_loc,dir_check)

                   IF(dir_check .NEQV. .TRUE.) THEN
                       CALL SYSTEM('mkdir ' // TRIM(rootdir))
                   END IF
                    
                END IF
                
            END IF
            
            ! ==================================================================
            
            
            ! Contrary to the Silo-HDF5 database format, handles of the Binary
            ! database master/root and slave/local process files are perfectly
            ! static throughout post-process. Hence, they are set here so that
            ! they do not have to be repetitively computed in later procedures.
            IF(format == 2) THEN
                IF(n == 0 .AND. proc_rank == 0) dbroot = 2
                dbfile = 1
            END IF
            
            
            ! Querying Number of Flow Variable(s) in Binary Output =============
            
            IF(format == 2) THEN
                
                ! Initializing the counter of the number of flow variable(s) to
                ! be written to the formatted database file(s)
                dbvars = 0
                
                ! Partial densities
                IF((model_eqns == 2) .OR. (model_eqns == 3)) THEN
                    DO i = 1, num_fluids
                        IF(                alpha_rho_wrt(i)               &
                                                .OR.                      &
                                 (cons_vars_wrt .OR. prim_vars_wrt)       ) THEN
                            dbvars = dbvars + 1
                        END IF
                    END DO
                END IF
                
                ! Density
                IF(                            rho_wrt                         &
                                                .OR.                           &
                   (model_eqns == 1 .AND. (cons_vars_wrt .OR. prim_vars_wrt))) &
                                                THEN
                    dbvars = dbvars + 1
                END IF
                
                ! Momentum
                DO i = 1, E_idx-mom_idx%beg
                    IF(mom_wrt(i) .OR. cons_vars_wrt) dbvars = dbvars + 1
                END DO
                
                ! Velocity
                DO i = 1, E_idx-mom_idx%beg
                    IF(vel_wrt(i) .OR. prim_vars_wrt) dbvars = dbvars + 1
                END DO
                
                ! Flux limiter function
                DO i = 1, E_idx-mom_idx%beg
                    IF (flux_wrt(i)) dbvars = dbvars + 1
                END DO

                ! Energy
                IF(E_wrt .OR. cons_vars_wrt) dbvars = dbvars + 1
                
                ! Pressure
                IF(pres_wrt .OR. prim_vars_wrt) dbvars = dbvars + 1
                
                ! Volume fraction(s)
                IF((model_eqns == 2) .OR. (model_eqns == 3)) THEN
                    
                    DO i = 1, num_fluids-1
                        IF(                  alpha_wrt(i)                 &
                                                .OR.                      &
                                 (cons_vars_wrt .OR. prim_vars_wrt)       ) THEN
                            dbvars = dbvars + 1
                        END IF
                    END DO
                    
                    IF(                 alpha_wrt(num_fluids)                  &
                                                .OR.                           &
                        (adv_alphan .AND. (cons_vars_wrt .OR. prim_vars_wrt))) &
                                                THEN
                        dbvars = dbvars + 1
                    END IF
                    
                END IF
                
                ! Specific heat ratio function
                IF(                           gamma_wrt                        &
                                                .OR.                           &
                   (model_eqns == 1 .AND. (cons_vars_wrt .OR. prim_vars_wrt))) &
                                                THEN
                    dbvars = dbvars + 1
                END IF
                
                ! Specific heat ratio
                IF(heat_ratio_wrt) dbvars = dbvars + 1
                
                ! Liquid stiffness function
                IF(                          pi_inf_wrt                        &
                                                .OR.                           &
                   (model_eqns == 1 .AND. (cons_vars_wrt .OR. prim_vars_wrt))) &
                                                THEN
                    dbvars = dbvars + 1
                END IF
                
                ! Liquid stiffness
                IF(pres_inf_wrt) dbvars = dbvars + 1
                
                ! Speed of sound
                IF(c_wrt) dbvars = dbvars + 1
                
                ! Vorticity
                IF (p > 0) THEN
                    DO i = 1, E_idx-mom_idx%beg
                        IF (omega_wrt(i)) dbvars = dbvars + 1
                    END DO
                ELSEIF (n > 0) THEN
                    DO i = 1, E_idx-cont_idx%end
                        IF (omega_wrt(i)) dbvars = dbvars + 1
                    END DO
                END IF
        
                ! Curvature
                DO i = 1, num_fluids
                    IF (kappa_wrt(i)) THEN
                        dbvars = dbvars + 1
                    END IF
                END DO
                
                ! Numerical Schlieren function
                IF(schlieren_wrt) dbvars = dbvars + 1
                
            END IF
            
            ! END: Querying Number of Flow Variable(s) in Binary Output ========
            
            
            ! Modifiying the value of the precision variable, which is used to
            ! indicate the floating point precision of the data that is stored
            ! in the formatted database file(s). This is performed so that this
            ! variable may be directly used as input to the functions that are
            ! in charge of writing the data. Only possible for Silo-HDF5 format.
             IF(format == 1) THEN
                IF(precision == 1) THEN          ! Single precision
                    precision = DB_FLOAT
                ELSE
                    precision = DB_DOUBLE        ! Double precision
                END IF
             END IF
             
             
        END SUBROUTINE s_initialize_data_output_module ! --------------------------
        
        
        
        
        
        SUBROUTINE s_open_formatted_database_file(t_step) ! --------------------
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
            INTEGER, INTENT(IN) :: t_step
            
            ! Generic string used to store the location of a particular file
            CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc
            
            
            ! Silo-HDF5 Database Format ========================================
            
            IF(format == 1) THEN
                
                ! Generating the relative path to the formatted database slave
                ! file, that is to be opened for the current time-step, t_step
                WRITE(file_loc, '(A,I0,A)') '/', t_step, '.silo'
                file_loc = TRIM(proc_rank_dir) // TRIM(file_loc)
                
                ! Creating formatted database slave file at the above location
                ! and setting up the structure of the file and its header info
                ierr = DBCREATE( TRIM(file_loc), LEN_TRIM(file_loc),  &
                                 DB_CLOBBER, DB_LOCAL, 'MFC v3.0', 8, &
                                 DB_HDF5, dbfile )
                
                ! Verifying that the creation and setup process of the formatted
                ! database slave file has been performed without errors. If this
                ! is not the case, the post-process exits.
                IF(dbfile == -1) THEN
                    PRINT '(A)', 'Unable to create Silo-HDF5 database ' // &
                                 'slave file '// TRIM(file_loc) // '. ' // &
                                 'Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
                ! Next, analogous steps to the ones above are carried out by the
                ! root process to create and setup the formatted database master
                ! file.
                IF(proc_rank == 0) THEN
                    
                    WRITE(file_loc, '(A,I0,A)') '/collection_', t_step, '.silo'
                    file_loc = TRIM(rootdir) // TRIM(file_loc)
                    
                    ierr = DBCREATE( TRIM(file_loc), LEN_TRIM(file_loc),  &
                                     DB_CLOBBER, DB_LOCAL, 'MFC v3.0', 8, &
                                     DB_HDF5, dbroot )
                    
                    IF(dbroot == -1) THEN
                        PRINT '(A)', 'Unable to create Silo-HDF5 database ' // &
                                     'master file '// TRIM(file_loc) //'. ' // &
                                     'Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                END IF
                
            ! ==================================================================
                
                
            ! Binary Database Format ===========================================
                
            ELSE
                
                ! Generating the relative path to the formatted database slave
                ! file, that is to be opened for the current time-step, t_step
                WRITE(file_loc, '(A,I0,A)') '/', t_step, '.dat'
                file_loc = TRIM(proc_rank_dir) // TRIM(file_loc)
                
                ! Creating the formatted database slave file, at the previously
                ! precised relative path location, and setting up its structure
                OPEN( dbfile, IOSTAT = err, FILE = TRIM(file_loc), &
                        FORM = 'unformatted', STATUS = 'replace'   )
                
                ! Verifying that the creation and setup process of the formatted
                ! database slave file has been performed without errors. If this
                ! is not the case, the post-process exits.
                IF(err /= 0) THEN
                    PRINT '(A)', 'Unable to create Binary database slave ' // &
                                 'file '// TRIM(file_loc)//'. Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
                ! Further defining the structure of the formatted database slave
                ! file by describing in it the dimensionality of post-processed
                ! data as well as the total number of flow variable(s) that will
                ! eventually be stored in it
                WRITE(dbfile) m, n, p, dbvars
                
                ! Next, analogous steps to the ones above are carried out by the
                ! root process to create and setup the formatted database master
                ! file. Note that this is only done in multidimensional cases.
                IF(n == 0 .AND. proc_rank == 0) THEN
                    
                    WRITE(file_loc, '(A,I0,A)') '/', t_step, '.dat'
                    file_loc = TRIM(rootdir) // TRIM(file_loc)
                    
                    OPEN( dbroot, IOSTAT = err, FILE = TRIM(file_loc), &
                            FORM = 'unformatted', STATUS = 'replace'   )
                     
                    IF(err /= 0) THEN
                        PRINT '(A)', 'Unable to create Binary database ' // &
                                     'master file ' // TRIM(file_loc)    // &
                                     '. Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    WRITE(dbroot) m_root, 0, 0, dbvars
                    
                END IF
                
            END IF
            
            ! END: Binary Database Format ======================================
            
            
        END SUBROUTINE s_open_formatted_database_file ! ------------------------
        
        
        
        
        
        SUBROUTINE s_write_grid_to_formatted_database_file(t_step) ! -----------
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
            INTEGER, INTENT(IN) :: t_step
            
            ! Bookkeeping variables storing the name and type of mesh that is
            ! handled by the local processor(s). Note that due to an internal
            ! NAG Fortran compiler problem, these two variables could not be
            ! allocated dynamically.
            CHARACTER(LEN = 4*name_len), DIMENSION(num_procs) :: meshnames
            INTEGER                    , DIMENSION(num_procs) :: meshtypes
            
            ! Generic loop iterator
            INTEGER :: i
            
            
            ! Silo-HDF5 Database Format ========================================
            
            IF(format == 1 .AND. n > 0) THEN
                
                ! For multidimensional data sets, the spatial extents of all of
                ! the grid(s) handled by the local processor(s) are recorded so
                ! that they may be written, by root processor, to the formatted
                ! database master file.
                IF(num_procs > 1) THEN
                    CALL s_mpi_gather_spatial_extents(spatial_extents)
                    
                ELSEIF(p > 0) THEN
                    IF (grid_geometry == 3) THEN
                        spatial_extents(:,0) = (/ MINVAL(y_cb), MINVAL(z_cb), &
                                                  MINVAL(x_cb), MAXVAL(y_cb), &
                                                  MAXVAL(z_cb), MAXVAL(x_cb) /)
                    ELSE
                        spatial_extents(:,0) = (/ MINVAL(x_cb), MINVAL(y_cb), &
                                                  MINVAL(z_cb), MAXVAL(x_cb), &
                                                  MAXVAL(y_cb), MAXVAL(z_cb) /)
                    END IF
                    
                ELSE
                    spatial_extents(:,0) = (/ MINVAL(x_cb), MINVAL(y_cb), &
                                              MAXVAL(x_cb), MAXVAL(y_cb) /)
                    
                END IF
                
                ! Next, the root processor proceeds to record all of the spatial
                ! extents in the formatted database master file. In addition, it
                ! also records a sub-domain connectivity map so that the entire
                ! grid may be reassembled by looking at the master file.
                IF(proc_rank == 0) THEN
                    
                    DO i = 1, num_procs
                        WRITE(meshnames(i), '(A,I0,A,I0,A)') '../p', i-1, &
                               '/', t_step, '.silo:rectilinear_grid'
                    END DO
                    
                    meshtypes = DB_QUAD_RECT
                    
                    err = DBSET2DSTRLEN(LEN(meshnames(1)))
                    err = DBMKOPTLIST(2, optlist)
                    err = DBADDIOPT( optlist, DBOPT_EXTENTS_SIZE,     &
                                              SIZE(spatial_extents,1) )
                    err = DBADDDOPT(optlist, DBOPT_EXTENTS, spatial_extents)
                    err = DBPUTMMESH( dbroot, 'rectilinear_grid', 16,   &
                                               num_procs, meshnames,    &
                                               LEN_TRIM(meshnames),     &
                                               meshtypes, optlist, ierr )
                    err = DBFREEOPTLIST(optlist)
                    
                END IF
                
                ! Finally, the local quadrilateral mesh, either 2D or 3D, along
                ! with its offsets that indicate the presence and size of ghost
                ! zone layer(s), are put in the formatted database slave file.

                IF (coarsen_silo) THEN
                    coarse_x_cb(-1-offset_x%beg : -1) = x_cb(-1-offset_x%beg : -1)
                    coarse_x_cb(m/2 : (m/2) + offset_x%end) = x_cb(m : m + offset_x%end)
                    DO i = 1, m, 2
                        coarse_x_cb((i-1)/2) = x_cb(i)
                    END DO

                    IF (n > 0) THEN
                        coarse_y_cb(-1-offset_y%beg : -1) = y_cb(-1-offset_y%beg : -1)
                        coarse_y_cb(n/2 : (n/2) + offset_y%end) = y_cb(n : n + offset_y%end)
                        DO i = 1, n, 2
                            coarse_y_cb((i-1)/2) = y_cb(i)
                        END DO

                        IF (p > 0) THEN
                            coarse_z_cb(-1-offset_z%beg : -1) = z_cb(-1-offset_z%beg : -1)
                            coarse_z_cb(p/2 : (p/2) + offset_z%end) = z_cb(p : p + offset_z%end)
                            DO i = 1, p, 2
                                coarse_z_cb((i-1)/2) = z_cb(i)
                            END DO
                        END IF
                    END IF
                END IF

                IF(p > 0) THEN
                    err = DBMKOPTLIST(2, optlist)
                    err = DBADDIOPT(optlist, DBOPT_LO_OFFSET, lo_offset)
                    err = DBADDIOPT(optlist, DBOPT_HI_OFFSET, hi_offset)
                    IF (coarsen_silo) THEN
                        IF (grid_geometry == 3) THEN
                            err = DBPUTQM( dbfile, 'rectilinear_grid', 16,       &
                                                      'x', 1, 'y', 1, 'z', 1,    &
                                                      coarse_y_cb, coarse_z_cb, coarse_x_cb, dims, 3, &
                                                      precision, DB_COLLINEAR,   &
                                                      optlist, ierr              )
                        ELSE
                            err = DBPUTQM( dbfile, 'rectilinear_grid', 16,       &
                                                      'x', 1, 'y', 1, 'z', 1,    &
                                                      coarse_x_cb, coarse_y_cb, coarse_z_cb, dims, 3, &
                                                      precision, DB_COLLINEAR,   &
                                                      optlist, ierr              )
                        END IF
                    ELSE
                        IF (grid_geometry == 3) THEN
                            err = DBPUTQM( dbfile, 'rectilinear_grid', 16,       &
                                                      'x', 1, 'y', 1, 'z', 1,    &
                                                      y_cb, z_cb, x_cb, dims, 3, &
                                                      precision, DB_COLLINEAR,   &
                                                      optlist, ierr              )
                        ELSE
                            err = DBPUTQM( dbfile, 'rectilinear_grid', 16,       &
                                                      'x', 1, 'y', 1, 'z', 1,    &
                                                      x_cb, y_cb, z_cb, dims, 3, &
                                                      precision, DB_COLLINEAR,   &
                                                      optlist, ierr              )
                        END IF
                    END IF
                    err = DBFREEOPTLIST(optlist)
                    
                ELSE
                    err = DBMKOPTLIST(2, optlist)
                    err = DBADDIOPT(optlist, DBOPT_LO_OFFSET, lo_offset)
                    err = DBADDIOPT(optlist, DBOPT_HI_OFFSET, hi_offset)
                    IF (coarsen_silo) THEN
                        err = DBPUTQM( dbfile, 'rectilinear_grid', 16,          &
                                               'x', 1, 'y', 1, 'z', 1,          &
                                               coarse_x_cb, coarse_y_cb, DB_F77NULL, dims, 2, &
                                               precision, DB_COLLINEAR,         &
                                               optlist, ierr                    )
                    ELSE
                        err = DBPUTQM( dbfile, 'rectilinear_grid', 16,          &
                                               'x', 1, 'y', 1, 'z', 1,          &
                                               x_cb, y_cb, DB_F77NULL, dims, 2, &
                                               precision, DB_COLLINEAR,         &
                                               optlist, ierr                    )
                    END IF
                    err = DBFREEOPTLIST(optlist)
                    
                END IF
                
            ! END: Silo-HDF5 Database Format ===================================
                
                
            ! Binary Database Format ===========================================
                
            ELSEIF(format == 2) THEN
                
                ! Multidimensional local grid data is written to the formatted
                ! database slave file. Recall that no master file to maintained
                ! in multidimensions.
                IF(p > 0) THEN
                    IF(precision == 1) THEN
                        WRITE(dbfile) REAL(x_cb, KIND(0.0)), &
                                      REAL(y_cb, KIND(0.0)), &
                                      REAL(z_cb, KIND(0.0))
                    ELSE
                        WRITE(dbfile) x_cb, y_cb, z_cb
                    END IF
                    
                ELSEIF(n > 0) THEN
                    IF(precision == 1) THEN
                        WRITE(dbfile) REAL(x_cb, KIND(0.0)), &
                                      REAL(y_cb, KIND(0.0))
                    ELSE
                        WRITE(dbfile) x_cb, y_cb
                    END IF
                    
                ! One-dimensional local grid data is written to the formatted
                ! database slave file. In addition, the local grid data is put
                ! together by the root process and written to the master file.
                ELSE
                    
                    IF(precision == 1) THEN
                        WRITE(dbfile) REAL(x_cb, KIND(0.0))
                    ELSE
                        WRITE(dbfile) x_cb
                    END IF
                    
                    IF(num_procs > 1) THEN
                        CALL s_mpi_defragment_1d_grid_variable()
                    ELSE
                        x_root_cb = x_cb
                    END IF
                    
                    IF(proc_rank == 0) THEN
                        IF(precision == 1) THEN
                            WRITE(dbroot) REAL(x_root_cb, KIND(0.0))
                        ELSE
                            WRITE(dbroot) x_root_cb
                        END IF
                    END IF
                    
                END IF
                
            END IF
            
            ! ==================================================================
            
            
        END SUBROUTINE s_write_grid_to_formatted_database_file ! ---------------
        
        
        
        
        
        SUBROUTINE s_write_variable_to_formatted_database_file(varname, t_step)
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
            CHARACTER(LEN = *), INTENT(IN) :: varname
            
            ! Time-step that is currently being post-processed
            INTEGER, INTENT(IN) :: t_step
            
            ! Bookkeeping variables storing the name and type of flow variable
            ! that is about to be handled by the local processor(s). Note that
            ! due to an internal NAG Fortran compiler problem, these variables
            ! could not be allocated dynamically.
            CHARACTER(LEN = 4*name_len), DIMENSION(num_procs) :: varnames
            INTEGER                    , DIMENSION(num_procs) :: vartypes
            
            ! Generic loop iterator
            INTEGER :: i,j,k
            
            ! Silo-HDF5 Database Format ========================================
            
            IF(format == 1) THEN
                
                ! In 1D, a curve object, featuring the local processor grid and
                ! flow variable data, is written to the formatted database slave
                ! file. The root process, on the other hand, will also take care
                ! of gathering the entire grid and associated flow variable data
                ! and write it to the formatted database master file.
                IF(n == 0) THEN
                    
                    ! Writing the curve object associated with the local process
                    ! to the formatted database slave file
                    err = DBPUTCURVE( dbfile, TRIM(varname), LEN_TRIM(varname),&
                                              x_cc(0:m), q_sf, precision, m+1, &
                                              DB_F77NULL, ierr                 )
                    
                    ! Assembling the local grid and flow variable data for the
                    ! entire computational domain on to the root process
                    IF(num_procs > 1) THEN
                       CALL s_mpi_defragment_1d_grid_variable()
                       CALL s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)
                    ELSE
                       x_root_cc = x_cc(0:m)
                       q_root_sf = q_sf
                    END IF
                    
                    ! Writing the curve object associated with the root process
                    ! to the formatted database master file
                    IF(proc_rank == 0) THEN
                        err = DBPUTCURVE( dbroot, TRIM(varname),        &
                                                  LEN_TRIM(varname),    &
                                                  x_root_cc, q_root_sf, &
                                                  precision, m_root+1,  &
                                                  DB_F77NULL, ierr      )
                    END IF
                    
                    RETURN
                    
                ! In multidimensions, the local process(es) take care of writing
                ! the flow variable data they are in charge of to the formatted
                ! database slave file. The root processor, additionally, is also
                ! responsible in gathering the flow variable extents of each of
                ! the local processor(s) and writing them to formatted database
                ! master file.
                ELSE
                    
                    ! Determining the extents of the flow variable on each local
                    ! process and gathering all this information on root process
                    IF(num_procs > 1) THEN
                        CALL s_mpi_gather_data_extents(q_sf, data_extents)
                    ELSE
                        data_extents(:,0) = (/ MINVAL(q_sf), MAXVAL(q_sf) /)
                    END IF
                    
                    ! Next, the root process proceeds to write the gathered flow
                    ! variable data extents to formatted database master file.
                    IF(proc_rank == 0) THEN
                        
                        DO i = 1, num_procs
                           WRITE(varnames(i), '(A,I0,A,I0,A)') '../p', i-1, &
                                  '/', t_step, '.silo:' // TRIM(varname)
                        END DO
                        
                        vartypes = DB_QUADVAR
                        
                        err = DBSET2DSTRLEN(LEN(varnames(1)))
                        err = DBMKOPTLIST(2, optlist)
                        err = DBADDIOPT(optlist, DBOPT_EXTENTS_SIZE, 2)
                        err = DBADDDOPT(optlist, DBOPT_EXTENTS, data_extents)
                        err = DBPUTMVAR( dbroot, TRIM(varname),                &
                                                 LEN_TRIM(varname), num_procs, &
                                                 varnames, LEN_TRIM(varnames), &
                                                 vartypes, optlist, ierr       )
                        err = DBFREEOPTLIST(optlist)
                        
                    END IF
                    
                    ! Finally, each of the local processor(s) proceeds to write
                    ! the flow variable data that it is responsible for to the
                    ! formatted database slave file.

                    IF (coarsen_silo) CALL s_coarsen_variable()

                    IF (grid_geometry == 3) THEN
                        IF (coarsen_silo) THEN
                            DO i = -offset_x%beg, (m/2) + offset_x%end
                                DO j = -offset_y%beg, (n/2) + offset_y%end
                                    DO k = -offset_z%beg, (p/2) + offset_z%end
                                        cyl_coarse_q_sf(j,k,i) = coarse_xyz_q_sf(i,j,k)
                                    END DO
                                END DO
                            END DO
                        ELSE
                            DO i = -offset_x%beg, m + offset_x%end
                                DO j = -offset_y%beg, n + offset_y%end
                                    DO k = -offset_z%beg, p + offset_z%end
                                        cyl_q_sf(j,k,i) = q_sf(i,j,k)
                                    END DO
                                END DO
                            END DO
                        END IF
                    END IF

                    IF(p > 0) THEN
                        IF (coarsen_silo) THEN
                            IF (grid_geometry == 3) THEN
                                err = DBPUTQV1( dbfile, TRIM(varname),               &
                                                        LEN_TRIM(varname),           &
                                                        'rectilinear_grid', 16,      &
                                                        cyl_coarse_q_sf, dims-1, 3, DB_F77NULL, &
                                                        0, precision, DB_ZONECENT,   &
                                                        DB_F77NULL, ierr             )
                            ELSE
                                err = DBPUTQV1( dbfile, TRIM(varname),               &
                                                        LEN_TRIM(varname),           &
                                                        'rectilinear_grid', 16,      &
                                                        coarse_xyz_q_sf, dims-1, 3, DB_F77NULL, &
                                                        0, precision, DB_ZONECENT,   &
                                                        DB_F77NULL, ierr             )
                            END IF
                        ELSE
                            IF (grid_geometry == 3) THEN
                                err = DBPUTQV1( dbfile, TRIM(varname),               &
                                                        LEN_TRIM(varname),           &
                                                        'rectilinear_grid', 16,      &
                                                        cyl_q_sf, dims-1, 3, DB_F77NULL, &
                                                        0, precision, DB_ZONECENT,   &
                                                        DB_F77NULL, ierr             )
                            ELSE
                                err = DBPUTQV1( dbfile, TRIM(varname),               &
                                                        LEN_TRIM(varname),           &
                                                        'rectilinear_grid', 16,      &
                                                        q_sf, dims-1, 3, DB_F77NULL, &
                                                        0, precision, DB_ZONECENT,   &
                                                        DB_F77NULL, ierr             )
                            END IF
                        END IF
                    ELSE
                        IF (coarsen_silo) THEN
                            err = DBPUTQV1( dbfile, TRIM(varname),               &
                                                    LEN_TRIM(varname),           &
                                                    'rectilinear_grid', 16,      &
                                                    coarse_xy_q_sf, dims-1, 2, DB_F77NULL, &
                                                    0, precision, DB_ZONECENT,   &
                                                    DB_F77NULL, ierr             )
                        ELSE
                            err = DBPUTQV1( dbfile, TRIM(varname),               &
                                                    LEN_TRIM(varname),           &
                                                    'rectilinear_grid', 16,      &
                                                    q_sf, dims-1, 2, DB_F77NULL, &
                                                    0, precision, DB_ZONECENT,   &
                                                    DB_F77NULL, ierr             )
                        END IF
                    END IF
                    
                END IF
                
            ! END: Silo-HDF5 Database Format ===================================
                
                
            ! Binary Database Format ===========================================
                
            ELSE
                
                ! Writing the name of the flow variable and its data, associated
                ! with the local processor, to the formatted database slave file
                IF(precision == 1) THEN
                    WRITE(dbfile) varname, REAL(q_sf, KIND(0.0))
                ELSE
                    WRITE(dbfile) varname, q_sf
                END IF
                
                ! In 1D, the root process also takes care of gathering the flow
                ! variable data from all of the local processor(s) and writes it
                ! to the formatted database master file.
                IF(n == 0) THEN
                    
                    IF(num_procs > 1) THEN
                       CALL s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)
                    ELSE
                       q_root_sf = q_sf
                    END IF
                    
                    IF(proc_rank == 0) THEN
                        IF(precision == 1) THEN
                            WRITE(dbroot) varname, REAL(q_root_sf, KIND(0.0))
                        ELSE
                            WRITE(dbroot) varname, q_root_sf
                        END IF
                    END IF
                    
                END IF
                
            END IF
            
            ! ==================================================================
            
            
        END SUBROUTINE s_write_variable_to_formatted_database_file ! -----------
        
        
        
        
        
        SUBROUTINE s_coarsen_variable() ! --------------------------------------
        ! Description: The purpose of this subroutine is to coarsen any variable
        ! that is to be written to the formatted database file by averaging every
        ! two cells together into a single value. This averaging is done separately
        ! in each dimension.

            ! Generic loop iterator
            INTEGER :: i,j,k

            ! Average q_sf onto coarser grid
            coarse_x_q_sf(-offset_x%beg : -1,:,:) = q_sf(-offset_x%beg : -1,:,:)
            coarse_x_q_sf((m/2)+1 : (m/2)+offset_x%end,:,:) = q_sf(m+1 : m+offset_x%end,:,:)
            DO i = 1, m, 2
                coarse_x_q_sf(   (i-1)/2,-offset_y%beg:n+offset_y%end,-offset_z%beg:p+offset_z%end)   &
                    = 5d-1*(q_sf(  i    ,-offset_y%beg:n+offset_y%end,-offset_z%beg:p+offset_z%end) + &
                            q_sf( i-1   ,-offset_y%beg:n+offset_y%end,-offset_z%beg:p+offset_z%end))
            END DO
            IF (MOD(m,2)==0) coarse_x_q_sf(m/2,:,:) = q_sf(m,:,:)

            IF (n > 0) THEN
                coarse_xy_q_sf(:,-offset_y%beg : -1,:) = coarse_x_q_sf(:,-offset_y%beg : -1,:)
                coarse_xy_q_sf(:,(n/2)+1 : (n/2)+offset_y%end,:) = coarse_x_q_sf(:,n+1 : n+offset_y%end,:)
                DO i = 1, n, 2
                    coarse_xy_q_sf(           -offset_x%beg:(m/2)+offset_x%end,(i-1)/2,-offset_z%beg:p+offset_z%end)   &
                        = 5d-1*(coarse_x_q_sf(-offset_x%beg:(m/2)+offset_x%end,  i    ,-offset_z%beg:p+offset_z%end) + &
                                coarse_x_q_sf(-offset_x%beg:(m/2)+offset_x%end, i-1   ,-offset_z%beg:p+offset_z%end))
                END DO
                IF (MOD(n,2)==0) coarse_xy_q_sf(:,n/2,:) = coarse_x_q_sf(:,n,:)

                IF (p > 0) THEN
                    coarse_xyz_q_sf(:,:,-offset_z%beg : -1) = coarse_xy_q_sf(:,:,-offset_z%beg : -1)
                    coarse_xyz_q_sf(:,:,(p/2)+1 : (p/2)+offset_z%end) = coarse_xy_q_sf(:,:,p+1 : p+offset_z%end)
                    DO i = 1, p, 2
                        coarse_xyz_q_sf(           -offset_x%beg:(m/2)+offset_x%end,-offset_y%beg:(n/2)+offset_y%end,(i-1)/2)   &
                            = 5d-1*(coarse_xy_q_sf(-offset_x%beg:(m/2)+offset_x%end,-offset_y%beg:(n/2)+offset_y%end,  i    ) + &
                                    coarse_xy_q_sf(-offset_x%beg:(m/2)+offset_x%end,-offset_y%beg:(n/2)+offset_y%end, i-1   ))
                    END DO
                    IF (MOD(p,2)==0) coarse_xyz_q_sf(:,:,p/2) = coarse_xy_q_sf(:,:,p)
                END IF
            END IF

        END SUBROUTINE s_coarsen_variable ! ------------------------------------





        SUBROUTINE s_close_formatted_database_file() ! -------------------------
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
            IF(format == 1) THEN
                ierr = DBCLOSE(dbfile)
                IF(proc_rank == 0) ierr = DBCLOSE(dbroot)
                
            ! Binary database format
            ELSE
                CLOSE(dbfile)
                IF(n == 0 .AND. proc_rank == 0) CLOSE(dbroot)
                
            END IF
            
            
        END SUBROUTINE s_close_formatted_database_file ! -----------------------
        
        
        
        
        
        SUBROUTINE s_finalize_data_output_module() ! -------------------------
        ! Description: Deallocation procedures for the module
            
            
            ! Deallocating the generic storage employed for the flow variable(s)
            ! that were written to the formatted database file(s). Note that the
            ! root variable is only deallocated in the case of a 1D computation.
            DEALLOCATE(q_sf)
            IF (fourier_decomp) DEALLOCATE(dft_q_sf)
            IF (coarsen_silo) DEALLOCATE(coarse_x_q_sf,coarse_xy_q_sf,coarse_xyz_q_sf)
            IF(n == 0) DEALLOCATE(q_root_sf)
            IF (grid_geometry == 3) THEN
                IF (coarsen_silo) THEN
                    DEALLOCATE(cyl_coarse_q_sf)
                ELSE
                    DEALLOCATE(cyl_q_sf)
                END IF
            END IF
            
            
            ! Deallocating spatial and data extents and also the variables for
            ! the offsets and the one bookkeeping the number of cell-boundaries
            ! in each active coordinate direction. Note that all these variables
            ! were only needed by Silo-HDF5 format for multidimensional data.
            IF(format == 1 .AND. n > 0) THEN
                DEALLOCATE(spatial_extents)
                DEALLOCATE(data_extents)
                DEALLOCATE(lo_offset)
                DEALLOCATE(hi_offset)
                DEALLOCATE(dims)
            END IF
            
            
        END SUBROUTINE s_finalize_data_output_module ! -----------------------
        
        
        
        
        
END MODULE m_data_output
