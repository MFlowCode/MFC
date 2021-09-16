!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!! Copyright 2021
!!
!! Permission is hereby granted, free of charge, to any person 
!! obtaining a copy of this software and associated documentation 
!! files (the "Software"), to deal in the Software without 
!! restriction, including without limitation the rights to use, 
!! copy, modify, merge, publish, distribute, sublicense, 
!! and/or sell copies of the Software, and to permit persons 
!! to whom the Software is furnished to do so, subject to the 
!! following conditions:
!!
!! The above copyright notice and this permission notice shall 
!! be included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
!! FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
!! OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
!! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
!! THE SOFTWARE.

!>
!! @file m_data_input.f90
!> @brief Contains module m_data_input
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module features procedures, which for a specific time-step,
!!             read in the raw simulation data for the grid and the conservative
!!             variables and fill out their buffer regions.
MODULE m_data_input
    
    
    ! Dependencies =============================================================
    USE mpi                     !< Message passing interface (MPI) module
    
    USE m_derived_types         !< Definitions of the derived types
    
    USE m_global_parameters     !< Global parameters for the code
    
    USE m_mpi_proxy             !< Message passing interface (MPI) module proxy

    USE m_compile_specific
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_data_input_module, &
                       s_read_data_files, &
                       s_read_serial_data_files, &
                       s_read_parallel_data_files, &
                       s_populate_grid_variables_buffer_regions, &
                       s_populate_conservative_variables_buffer_regions, &
                       s_finalize_data_input_module

    ABSTRACT INTERFACE ! ===================================================

        !> Subroutine for reading data files
        !!  @param t_step Current time-step to input
        SUBROUTINE s_read_abstract_data_files(t_step) ! ------------

            INTEGER, INTENT(IN) :: t_step

        END SUBROUTINE s_read_abstract_data_files ! ----------------

    END INTERFACE ! ========================================================


    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:), PUBLIC :: q_cons_vf !<
    !! Conservative variables
    
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:), PUBLIC :: q_prim_vf !<
    !! Primitive variables
    
    PROCEDURE(s_read_abstract_data_files), POINTER :: s_read_data_files => NULL()
    
    CONTAINS
        
        
        
        
        !>  This subroutine is called at each time-step that has to
        !!      be post-processed in order to read the raw data files
        !!      present in the corresponding time-step directory and to
        !!      populate the associated grid and conservative variables.        
        !!  @param t_step Current time-step
        SUBROUTINE s_read_serial_data_files(t_step) ! -----------------------------

            INTEGER, INTENT(IN) :: t_step
            
            CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir !<
            !! Location of the time-step directory associated with t_step
            
            CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc !<
            !! Generic string used to store the location of a particular file
            
            CHARACTER(LEN = &
            INT(FLOOR(LOG10(REAL(sys_size, KIND(0d0))))) + 1) :: file_num !<
            !! Used to store the variable position, in character form, of the
            !! currently manipulated conservative variable file
            
            LOGICAL :: dir_check !<
            !! Generic logical used to test the existence of a particular folder
            
            LOGICAL :: file_check  !<
            !! Generic logical used to test the existence of a particular file

            INTEGER :: i !< Generic loop iterator
            
            
            ! Setting location of time-step folder based on current time-step
            WRITE(t_step_dir, '(A,I0,A,I0)') '/p_all/p', proc_rank, '/', t_step
            t_step_dir = TRIM(case_dir) // TRIM(t_step_dir)
            
            
            ! Inquiring as to the existence of the time-step directory
            file_loc = TRIM(t_step_dir) // '/.'
            
            CALL my_inquire(file_loc,dir_check) 
            
            ! If the time-step directory is missing, the post-process exits.
            IF(dir_check .NEQV. .TRUE.) THEN
                PRINT '(A)', 'Time-step folder ' // TRIM(t_step_dir) // &
                             ' is missing. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
            ! Reading the Grid Data File for the x-direction ===================
            
            ! Checking whether x_cb.dat exists
            file_loc = TRIM(t_step_dir) // '/x_cb.dat'
            INQUIRE(FILE = TRIM(file_loc), EXIST = file_check)
            
            ! Reading x_cb.dat if it exists, exiting otherwise
            IF(file_check) THEN
                OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                        STATUS = 'old', ACTION = 'read')
                READ(1) x_cb(-1:m)
                CLOSE(1)
            ELSE
                PRINT '(A)', 'File x_cb.dat is missing in ' // &
                             TRIM(t_step_dir) // '. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            ! Computing the cell-width distribution
            dx(0:m) = x_cb(0:m) - x_cb(-1:m-1)
            
            ! Computing the cell-center locations
            x_cc(0:m) = x_cb(-1:m-1) + dx(0:m)/2d0
            
            ! ==================================================================
            
            
            ! Reading the Grid Data File for the y-direction ===================
            
            IF(n > 0) THEN
                
                ! Checking whether y_cb.dat exists
                file_loc = TRIM(t_step_dir) // '/y_cb.dat'
                INQUIRE(FILE = TRIM(file_loc), EXIST = file_check)
                
                ! Reading y_cb.dat if it exists, exiting otherwise
                IF(file_check) THEN
                    OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                            STATUS = 'old', ACTION = 'read')
                    READ(1) y_cb(-1:n)
                    CLOSE(1)
                ELSE
                    PRINT '(A)', 'File y_cb.dat is missing in ' // &
                                 TRIM(t_step_dir) // '. Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
                ! Computing the cell-width distribution
                dy(0:n) = y_cb(0:n) - y_cb(-1:n-1)
                
                ! Computing the cell-center locations
                y_cc(0:n) = y_cb(-1:n-1) + dy(0:n)/2d0
                
            ! ==================================================================
                
                
            ! Reading the Grid Data File for the z-direction ===================
                
                IF(p > 0) THEN
                    
                    ! Checking whether z_cb.dat exists
                    file_loc = TRIM(t_step_dir) // '/z_cb.dat'
                    INQUIRE(FILE = TRIM(file_loc), EXIST = file_check)
                    
                    ! Reading z_cb.dat if it exists, exiting otherwise
                    IF(file_check) THEN
                        OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                                STATUS = 'old', ACTION = 'read')
                        READ(1) z_cb(-1:p)
                        CLOSE(1)
                    ELSE
                        PRINT '(A)', 'File z_cb.dat is missing in ' // &
                                     TRIM(t_step_dir) // '. Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    ! Computing the cell-width distribution
                    dz(0:p) = z_cb(0:p) - z_cb(-1:p-1)
                    
                    ! Computing the cell-center locations
                    z_cc(0:p) = z_cb(-1:p-1) + dz(0:p)/2d0
                    
                END IF
                
            END IF
            
            ! ==================================================================
            
            
            ! Reading the Conservative Variables Data Files ====================
            DO i = 1, sys_size 
                
                ! Checking whether the data file associated with the variable
                ! position of currently manipulated conservative variable exists
                WRITE(file_num, '(I0)') i
                file_loc = TRIM(t_step_dir) // '/q_cons_vf' // &
                           TRIM(file_num) // '.dat'
                INQUIRE(FILE = TRIM(file_loc), EXIST = file_check)
                
                ! Reading the data file if it exists, exiting otherwise
                IF(file_check) THEN
                    OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                            STATUS = 'old', ACTION = 'read')
                    READ(1) q_cons_vf(i)%sf(0:m,0:n,0:p)
                    CLOSE(1)
                ELSE
                    PRINT '(A)', 'File q_cons_vf' // TRIM(file_num) // &
                                 '.dat is missing in ' // TRIM(t_step_dir) // &
                                 '. Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
            END DO
            
            ! ==================================================================
            
            
        END SUBROUTINE s_read_serial_data_files ! ---------------------------------
        
        
        
        !>  This subroutine is called at each time-step that has to
        !!      be post-processed in order to parallel-read the raw data files
        !!      present in the corresponding time-step directory and to
        !!      populate the associated grid and conservative variables.        
        !!  @param t_step Current time-step        
        SUBROUTINE s_read_parallel_data_files(t_step) ! ---------------------------

            INTEGER, INTENT(IN) :: t_step

            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: x_cb_glb, y_cb_glb, z_cb_glb

            INTEGER :: ifile, ierr, data_size
            INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
            INTEGER(KIND=MPI_OFFSET_KIND) :: disp
            INTEGER(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: MOK

            CHARACTER(LEN=path_len + 2*name_len) :: file_loc
            LOGICAL :: file_exist

            INTEGER :: i

            ALLOCATE(x_cb_glb(-1:m_glb))
            ALLOCATE(y_cb_glb(-1:n_glb))
            ALLOCATE(z_cb_glb(-1:p_glb))

            ! Read in cell boundary locations in x-direction
            file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'x_cb.dat'
            INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)

            IF (file_exist) THEN
                data_size = m_glb+2
                CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
                CALL MPI_FILE_READ(ifile,x_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                CALL MPI_FILE_CLOSE(ifile,ierr)
            ELSE
                PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting...'
                CALL s_mpi_abort()
            END IF

            ! Assigning local cell boundary locations
            x_cb(-1:m) = x_cb_glb((start_idx(1)-1):(start_idx(1)+m))
            ! Computing the cell width distribution
            dx(0:m) = x_cb(0:m) - x_cb(-1:m-1)
            ! Computing the cell center location
            x_cc(0:m) = x_cb(-1:m-1) + dx(0:m)/2d0

            IF (n > 0) THEN
                ! Read in cell boundary locations in y-direction
                file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'y_cb.dat'
                INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
        
                IF (file_exist) THEN
                    data_size = n_glb+2
                    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
                    CALL MPI_FILE_READ(ifile,y_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                    CALL MPI_FILE_CLOSE(ifile,ierr)
                ELSE
                    PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting...'
                    CALL s_mpi_abort()
                END IF
        
                ! Assigning local cell boundary locations
                y_cb(-1:n) = y_cb_glb((start_idx(2)-1):(start_idx(2)+n))
                ! Computing the cell width distribution
                dy(0:n) = y_cb(0:n) - y_cb(-1:n-1)
                ! Computing the cell center location
                y_cc(0:n) = y_cb(-1:n-1) + dy(0:n)/2d0

                IF (p > 0) THEN
                    ! Read in cell boundary locations in z-direction
                    file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'z_cb.dat'
                    INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
            
                    IF (file_exist) THEN
                        data_size = p_glb+2
                        CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
                        CALL MPI_FILE_READ(ifile,z_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                        CALL MPI_FILE_CLOSE(ifile,ierr)
                    ELSE
                        PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting...'
                        CALL s_mpi_abort()
                    END IF
            
                    ! Assigning local cell boundary locations
                    z_cb(-1:p) = z_cb_glb((start_idx(3)-1):(start_idx(3)+p))
                    ! Computing the cell width distribution
                    dz(0:p) = z_cb(0:p) - z_cb(-1:p-1)
                    ! Computing the cell center location
                    z_cc(0:p) = z_cb(-1:p-1) + dz(0:p)/2d0
                END IF
            END IF

            ! Open the file to read conservative variables
            WRITE(file_loc, '(I0,A)') t_step, '.dat'
            file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
            INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)

            IF (file_exist) THEN
                CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)

                ! Initialize MPI data I/O
                CALL s_initialize_mpi_data(q_cons_vf)

                ! Size of local arrays
                data_size = (m+1)*(n+1)*(p+1)

                ! Resize some integers so MPI can read even the biggest file
                m_MOK     = INT(m_glb+1,     MPI_OFFSET_KIND)
                n_MOK     = INT(n_glb+1,     MPI_OFFSET_KIND)
                p_MOK     = INT(p_glb+1,     MPI_OFFSET_KIND)
                WP_MOK    = INT(8d0,      MPI_OFFSET_KIND)
                MOK       = INT(1d0,      MPI_OFFSET_KIND)
                str_MOK   = INT(name_len, MPI_OFFSET_KIND)
                NVARS_MOK = INT(sys_size, MPI_OFFSET_KIND)

                ! Read the data for each variable
                IF (bubbles) THEN
                    DO i = 1,sys_size
                        var_MOK = INT(i, MPI_OFFSET_KIND)

                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*MAX(MOK,n_MOK)*MAX(MOK,p_MOK)*WP_MOK*(var_MOK-1)

                        CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i), &
                                    'native',mpi_info_int,ierr)
                        CALL MPI_FILE_READ(ifile,MPI_IO_DATA%var(i)%sf,data_size, &
                                    MPI_DOUBLE_PRECISION,status,ierr)
                    END DO
                ELSE
                    DO i = 1, adv_idx%end
                        var_MOK = INT(i, MPI_OFFSET_KIND)

                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*MAX(MOK,n_MOK)*MAX(MOK,p_MOK)*WP_MOK*(var_MOK-1)

                        CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i), &
                                    'native',mpi_info_int,ierr)
                        CALL MPI_FILE_READ(ifile,MPI_IO_DATA%var(i)%sf,data_size, &
                                    MPI_DOUBLE_PRECISION,status,ierr)
                    END DO
                END IF

                CALL s_mpi_barrier()

                CALL MPI_FILE_CLOSE(ifile,ierr)
            ELSE
                PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting...'
                CALL s_mpi_abort()
            END IF

            DEALLOCATE(x_cb_glb, y_cb_glb, z_cb_glb)
        END SUBROUTINE s_read_parallel_data_files ! -------------------------------




        !>  The following subroutine populates the buffer regions of
        !!      the cell-width spacings, the cell-boundary locations and
        !!      the cell-center locations. Note that the buffer regions
        !!      of the last two variables should be interpreted slightly
        !!      differently than usual. They are really ghost zones that
        !!      are used in aiding the multidimensional visualization of
        !!      Silo database files, in VisIt, when processor boundary
        !!      conditions are present.
        SUBROUTINE s_populate_grid_variables_buffer_regions() ! ----------------

            INTEGER :: i !< Generic loop iterator
            
            ! Populating Buffer Regions in the x-direction =====================
            
            ! Ghost-cell extrapolation BC at the beginning
            IF(bc_x%beg <= -3) THEN
                
                DO i = 1, buff_size
                    dx(-i) = dx(0)
                END DO
                
            ! Symmetry BC at the beginning
            ELSEIF(bc_x%beg == -2) THEN
                
                DO i = 1, buff_size
                    dx(-i) = dx(i-1)
                END DO
                
            ! Periodic BC at the beginning
            ELSEIF(bc_x%beg == -1) THEN
                
                DO i = 1, buff_size
                    dx(-i) = dx((m+1)-i)
                END DO
                
            ! Processor BC at the beginning
            ELSE
                
                CALL s_mpi_sendrecv_grid_vars_buffer_regions('beg','x')
                
                DO i = 1, offset_x%beg
                    x_cb(-1-i) = x_cb(-i) - dx(-i)
                END DO
                
                DO i = 1, buff_size
                    x_cc(-i) = x_cc(1-i) - (dx(1-i) + dx(-i))/2d0
                END DO
                
            END IF
            
            ! Ghost-cell extrapolation BC at the end
            IF(bc_x%end <= -3) THEN
                
                DO i = 1, buff_size
                    dx(m+i) = dx(m)
                END DO
                
            ! Symmetry BC at the end
            ELSEIF(bc_x%end == -2) THEN
                
                DO i = 1, buff_size
                    dx(m+i) = dx((m+1)-i)
                END DO
                
            ! Periodic BC at the end
            ELSEIF(bc_x%end == -1) THEN
                
                DO i = 1, buff_size
                    dx(m+i) = dx(i-1)
                END DO
                
            ! Processor BC at the end
            ELSE
                
                CALL s_mpi_sendrecv_grid_vars_buffer_regions('end','x')
                
                DO i = 1, offset_x%end
                    x_cb(m+i) = x_cb(m+(i-1)) + dx(m+i)
                END DO
                
                DO i = 1, buff_size
                    x_cc(m+i) = x_cc(m+(i-1)) + (dx(m+(i-1)) + dx(m+i))/2d0
                END DO
                
            END IF
            
            ! END: Populating Buffer Regions in the x-direction ================
            
            
            ! Populating Buffer Regions in the y-direction =====================
            
            IF(n > 0) THEN
                
                ! Ghost-cell extrapolation BC at the beginning
                IF(bc_y%beg <= -3 .AND. bc_y%beg /= -13) THEN
                    
                    DO i = 1, buff_size
                        dy(-i) = dy(0)
                    END DO
                    
                ! Symmetry BC at the beginning
                ELSEIF(bc_y%beg == -2 .OR. bc_y%beg == -13) THEN
                    
                    DO i = 1, buff_size
                        dy(-i) = dy(i-1)
                    END DO
                    
                ! Periodic BC at the beginning
                ELSEIF(bc_y%beg == -1) THEN
                    
                    DO i = 1, buff_size
                        dy(-i) = dy((n+1)-i)
                    END DO
                    
                ! Processor BC at the beginning
                ELSE
                    
                    CALL s_mpi_sendrecv_grid_vars_buffer_regions('beg','y')
                    
                    DO i = 1, offset_y%beg
                        y_cb(-1-i) = y_cb(-i) - dy(-i)
                    END DO
                    
                    DO i = 1, buff_size
                        y_cc(-i) = y_cc(1-i) - (dy(1-i) + dy(-i))/2d0
                    END DO
                    
                END IF
                
                ! Ghost-cell extrapolation BC at the end
                IF(bc_y%end <= -3) THEN
                    
                    DO i = 1, buff_size
                        dy(n+i) = dy(n)
                    END DO
                    
                ! Symmetry BC at the end
                ELSEIF(bc_y%end == -2) THEN
                    
                    DO i = 1, buff_size
                        dy(n+i) = dy((n+1)-i)
                    END DO
                    
                ! Periodic BC at the end
                ELSEIF(bc_y%end == -1) THEN
                    
                    DO i = 1, buff_size
                        dy(n+i) = dy(i-1)
                    END DO
                    
                ! Processor BC at the end
                ELSE
                    
                    CALL s_mpi_sendrecv_grid_vars_buffer_regions('end','y')
                    
                    DO i = 1, offset_y%end
                        y_cb(n+i) = y_cb(n+(i-1)) + dy(n+i)
                    END DO
                    
                    DO i = 1, buff_size
                        y_cc(n+i) = y_cc(n+(i-1)) + (dy(n+(i-1)) + dy(n+i))/2d0
                    END DO
                    
                END IF
                
            ! END: Populating Buffer Regions in the y-direction ================
                
                
            ! Populating Buffer Regions in the z-direction =====================
                
                IF(p > 0) THEN
                    
                    ! Ghost-cell extrapolation BC at the beginning
                    IF(bc_z%beg <= -3) THEN
                        
                        DO i = 1, buff_size
                            dz(-i) = dz(0)
                        END DO
                        
                    ! Symmetry BC at the beginning
                    ELSEIF(bc_z%beg == -2) THEN
                        
                        DO i = 1, buff_size
                            dz(-i) = dz(i-1)
                        END DO
                        
                    ! Periodic BC at the beginning
                    ELSEIF(bc_z%beg == -1) THEN
                        
                        DO i = 1, buff_size
                            dz(-i) = dz((p+1)-i)
                        END DO
                        
                    ! Processor BC at the beginning
                    ELSE
                        
                        CALL s_mpi_sendrecv_grid_vars_buffer_regions('beg','z')
                        
                        DO i = 1, offset_z%beg
                            z_cb(-1-i) = z_cb(-i) - dz(-i)
                        END DO
                        
                        DO i = 1, buff_size
                            z_cc(-i) = z_cc(1-i) - (dz(1-i) + dz(-i))/2d0
                        END DO
                        
                    END IF
                    
                    ! Ghost-cell extrapolation BC at the end
                    IF(bc_z%end <= -3) THEN
                        
                        DO i = 1, buff_size
                            dz(p+i) = dz(p)
                        END DO
                        
                    ! Symmetry BC at the end
                    ELSEIF(bc_z%end == -2) THEN
                        
                        DO i = 1, buff_size
                            dz(p+i) = dz((p+1)-i)
                        END DO
                        
                    ! Periodic BC at the end
                    ELSEIF(bc_z%end == -1) THEN
                        
                        DO i = 1, buff_size
                            dz(p+i) = dz(i-1)
                        END DO
                        
                    ! Processor BC at the end
                    ELSE
                        
                        CALL s_mpi_sendrecv_grid_vars_buffer_regions('end','z')
                        
                        DO i = 1, offset_z%end
                            z_cb(p+i) = z_cb(p+(i-1)) + dz(p+i)
                        END DO
                        
                        DO i = 1, buff_size
                            z_cc(p+i) = z_cc(p+(i-1)) + (dz(p+(i-1)) + dz(p+i))/2d0
                        END DO
                        
                    END IF
                    
                END IF
                
            END IF
            
            ! END: Populating Buffer Regions in the z-direction ================
            
            
        END SUBROUTINE s_populate_grid_variables_buffer_regions ! --------------
        
        
        
        
        !>  The purpose of this procedure is to populate the buffers
        !!      of the cell-average conservative variables, depending on
        !!      the boundary conditions.        
        SUBROUTINE s_populate_conservative_variables_buffer_regions() ! --------


            INTEGER :: i,j,k !< Generic loop iterators
            
            
            ! Populating Buffer Regions in the x-direction =====================
            
            ! Ghost-cell extrapolation BC at the beginning
            IF(bc_x%beg <= -3) THEN
                
                DO j = 1, buff_size
                    DO i = 1, sys_size
                        q_cons_vf(i)%sf(-j,0:n,0:p) = q_cons_vf(i)%sf(0,0:n,0:p)
                    END DO
                END DO
                
            ! Symmetry BC at the beginning
            ELSEIF(bc_x%beg == -2) THEN
                
                DO j = 1, buff_size
                    
                    ! Density or partial densities
                    DO i = 1, cont_idx%end
                        q_cons_vf(i)%sf(-j,0:n,0:p) = &
                                     q_cons_vf(i)%sf(j-1,0:n,0:p)
                    END DO
                    
                    ! x-component of momentum
                    q_cons_vf(mom_idx%beg)%sf(-j,0:n,0:p) = &
                                -q_cons_vf(mom_idx%beg)%sf(j-1,0:n,0:p)
                    
                    ! Remaining momentum component(s), if any, as well as the
                    ! energy and the variable(s) from advection equation(s)
                    DO i = mom_idx%beg+1, sys_size
                        q_cons_vf(i)%sf(-j,0:n,0:p) = &
                                     q_cons_vf(i)%sf(j-1,0:n,0:p)
                    END DO
                    
                END DO
                
            ! Periodic BC at the beginning
            ELSEIF(bc_x%beg == -1) THEN
                
                DO j = 1, buff_size
                    DO i = 1, sys_size
                        q_cons_vf(i)%sf(-j,0:n,0:p) = &
                                    q_cons_vf(i)%sf((m+1)-j,0:n,0:p)
                    END DO
                END DO
                
            ! Processor BC at the beginning
            ELSE
                
                CALL s_mpi_sendrecv_cons_vars_buffer_regions( q_cons_vf, &
                                                              'beg', 'x' )
                
            END IF
            
            ! Ghost-cell extrapolation BC at the end
            IF(bc_x%end <= -3) THEN
                
                DO j = 1, buff_size
                    DO i = 1, sys_size
                        q_cons_vf(i)%sf(m+j,0:n,0:p) = &
                                      q_cons_vf(i)%sf(m,0:n,0:p)
                    END DO
                END DO
                
            ! Symmetry BC at the end
            ELSEIF(bc_x%end == -2) THEN
                
                DO j = 1, buff_size
                    
                    ! Density or partial densities
                    DO i = 1, cont_idx%end
                        q_cons_vf(i)%sf(m+j,0:n,0:p) = &
                                    q_cons_vf(i)%sf((m+1)-j,0:n,0:p)
                    END DO
                    
                    ! x-component of momentum
                    q_cons_vf(mom_idx%beg)%sf(m+j,0:n,0:p) = &
                             -q_cons_vf(mom_idx%beg)%sf((m+1)-j,0:n,0:p)
                    
                    ! Remaining momentum component(s), if any, as well as the
                    ! energy and the variable(s) from advection equation(s)
                    DO i = mom_idx%beg+1, sys_size
                        q_cons_vf(i)%sf(m+j,0:n,0:p) = &
                                    q_cons_vf(i)%sf((m+1)-j,0:n,0:p)
                    END DO
                    
                END DO
                
            ! Perodic BC at the end
            ELSEIF(bc_x%end == -1) THEN
                
                DO j = 1, buff_size
                    DO i = 1, sys_size
                        q_cons_vf(i)%sf(m+j,0:n,0:p) = &
                                      q_cons_vf(i)%sf(j-1,0:n,0:p)
                    END DO
                END DO
                
            ! Processor BC at the end
            ELSE
                
                CALL s_mpi_sendrecv_cons_vars_buffer_regions( q_cons_vf, &
                                                              'end', 'x' )
                
            END IF
            
            ! END: Populating Buffer Regions in the x-direction ================
            
            
            ! Populating Buffer Regions in the y-direction =====================
            
            IF(n > 0) THEN
                
                ! Ghost-cell extrapolation BC at the beginning
                IF(bc_y%beg <= -3 .AND. bc_y%beg /= -13) THEN
                    
                    DO j = 1, buff_size
                        DO i = 1, sys_size
                            q_cons_vf(i)%sf(:,-j,0:p) = q_cons_vf(i)%sf(:,0,0:p)
                        END DO
                    END DO
                
                ! Axis BC at the beginning
                ELSEIF(bc_y%beg == -13) THEN

                    DO j = 1, buff_size
                        DO k = 0, p
                            IF (z_cc(k) < pi) THEN
                                DO i = 1, mom_idx%beg
                                    q_cons_vf(i)%sf(:,-j ,     k     ) = &
                                    q_cons_vf(i)%sf(:,j-1,k+((p+1)/2))
                                END DO
    
                                 q_cons_vf(mom_idx%beg+1)%sf(:,-j ,     k     ) = &
                                -q_cons_vf(mom_idx%beg+1)%sf(:,j-1,k+((p+1)/2))
                                
                                 q_cons_vf(mom_idx%end)%sf(:,-j ,     k     ) = &
                                -q_cons_vf(mom_idx%end)%sf(:,j-1,k+((p+1)/2))
                                
                                DO i = E_idx, sys_size
                                    q_cons_vf(i)%sf(:,-j ,     k     ) = &
                                    q_cons_vf(i)%sf(:,j-1,k+((p+1)/2))
                                END DO
                            ELSE
                                DO i = 1, mom_idx%beg
                                    q_cons_vf(i)%sf(:,-j ,     k     ) = &
                                    q_cons_vf(i)%sf(:,j-1,k-((p+1)/2))
                                END DO
    
                                 q_cons_vf(mom_idx%beg+1)%sf(:,-j ,     k     ) = &
                                -q_cons_vf(mom_idx%beg+1)%sf(:,j-1,k-((p+1)/2))
                                
                                 q_cons_vf(mom_idx%end)%sf(:,-j ,     k     ) = &
                                -q_cons_vf(mom_idx%end)%sf(:,j-1,k-((p+1)/2))
                                
                                DO i = E_idx, sys_size
                                    q_cons_vf(i)%sf(:,-j ,     k     ) = &
                                    q_cons_vf(i)%sf(:,j-1,k-((p+1)/2))
                                END DO
                            END IF
                        END DO
                    END DO

                ! Symmetry BC at the beginning
                ELSEIF(bc_y%beg == -2) THEN
                    
                    DO j = 1, buff_size
                        
                        ! Density or partial densities and x-momentum component
                        DO i = 1, mom_idx%beg
                            q_cons_vf(i)%sf(:,-j,0:p) = &
                                        q_cons_vf(i)%sf(:,j-1,0:p)
                        END DO
                        
                        ! y-component of momentum
                        q_cons_vf(mom_idx%beg+1)%sf(:,-j,0:p) = &
                                 -q_cons_vf(mom_idx%beg+1)%sf(:,j-1,0:p)
                        
                        ! Remaining z-momentum component, if any, as well as the
                        ! energy and variable(s) from advection equation(s)
                        DO i = mom_idx%beg+2, sys_size
                            q_cons_vf(i)%sf(:,-j,0:p) = &
                                        q_cons_vf(i)%sf(:,j-1,0:p)
                        END DO
                        
                    END DO
                    
                ! Periodic BC at the beginning
                ELSEIF(bc_y%beg == -1) THEN
                    
                    DO j = 1, buff_size
                        DO i = 1, sys_size
                            q_cons_vf(i)%sf(:,-j,0:p) = &
                                      q_cons_vf(i)%sf(:,(n+1)-j,0:p)
                        END DO
                    END DO
                    
                ! Processor BC at the beginning
                ELSE
                    
                    CALL s_mpi_sendrecv_cons_vars_buffer_regions( q_cons_vf, &
                                                                  'beg', 'y' )
                    
                END IF
                
                ! Ghost-cell extrapolation BC at the end
                IF(bc_y%end <= -3) THEN
                    
                    DO j = 1, buff_size
                        DO i = 1, sys_size
                            q_cons_vf(i)%sf(:,n+j,0:p) = &
                                          q_cons_vf(i)%sf(:,n,0:p)
                        END DO
                    END DO
                    
                ! Symmetry BC at the end
                ELSEIF(bc_y%end == -2) THEN
                    
                    DO j = 1, buff_size
                        
                        ! Density or partial densities and x-momentum component
                        DO i = 1, mom_idx%beg
                            q_cons_vf(i)%sf(:,n+j,0:p) = &
                                        q_cons_vf(i)%sf(:,(n+1)-j,0:p)
                        END DO
                        
                        ! y-component of momentum
                        q_cons_vf(mom_idx%beg+1)%sf(:,n+j,0:p) = &
                                 -q_cons_vf(mom_idx%beg+1)%sf(:,(n+1)-j,0:p)
                        
                        ! Remaining z-momentum component, if any, as well as the 
                        ! energy and variable(s) from advection equation(s)
                        DO i = mom_idx%beg+2, sys_size
                            q_cons_vf(i)%sf(:,n+j,0:p) = &
                                        q_cons_vf(i)%sf(:,(n+1)-j,0:p)
                        END DO
                        
                    END DO
                    
                ! Perodic BC at the end
                ELSEIF(bc_y%end == -1) THEN
                    
                    DO j = 1, buff_size
                        DO i = 1, sys_size
                            q_cons_vf(i)%sf(:,n+j,0:p) = &
                                          q_cons_vf(i)%sf(:,j-1,0:p)
                        END DO
                    END DO
                    
                ! Processor BC at the end
                ELSE
                    
                    CALL s_mpi_sendrecv_cons_vars_buffer_regions( q_cons_vf, &
                                                                  'end', 'y' )
                    
                END IF
                
            ! END: Populating Buffer Regions in the y-direction ================
                
                
            ! Populating Buffer Regions in the z-direction =====================
                
                IF(p > 0) THEN
                    
                    ! Ghost-cell extrapolation BC at the beginning
                    IF(bc_z%beg <= -3) THEN
                        
                        DO j = 1, buff_size
                            DO i = 1, sys_size
                                q_cons_vf(i)%sf(:,:,-j) = q_cons_vf(i)%sf(:,:,0)
                            END DO
                        END DO
                        
                    ! Symmetry BC at the beginning
                    ELSEIF(bc_z%beg == -2) THEN
                        
                        DO j = 1, buff_size
                            
                            ! Density or the partial densities and the momentum
                            ! components in x- and y-directions
                            DO i = 1, mom_idx%beg+1
                                q_cons_vf(i)%sf(:,:,-j) = &
                                            q_cons_vf(i)%sf(:,:,j-1)
                            END DO
                            
                            ! z-component of momentum
                            q_cons_vf(mom_idx%end)%sf(:,:,-j) = &
                                      -q_cons_vf(mom_idx%end)%sf(:,:,j-1)
                            
                            ! Energy and advection equation(s) variable(s)
                            DO i = E_idx, sys_size
                                q_cons_vf(i)%sf(:,:,-j) = &
                                            q_cons_vf(i)%sf(:,:,j-1)
                            END DO
                            
                        END DO
                        
                    ! Periodic BC at the beginning
                    ELSEIF(bc_z%beg == -1) THEN
                        
                        DO j = 1, buff_size
                            DO i = 1, sys_size
                                q_cons_vf(i)%sf(:,:,-j) = &
                                          q_cons_vf(i)%sf(:,:,(p+1)-j)
                            END DO
                        END DO
                        
                    ! Processor BC at the beginning
                    ELSE
                        
                        CALL s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf,&
                                                                     'beg', 'z')
                        
                    END IF
                    
                    ! Ghost-cell extrapolation BC at the end
                    IF(bc_z%end <= -3) THEN
                        
                        DO j = 1, buff_size
                            DO i = 1, sys_size
                                q_cons_vf(i)%sf(:,:,p+j) = &
                                             q_cons_vf(i)%sf(:,:,p)
                            END DO
                        END DO
                        
                    ! Symmetry BC at the end
                    ELSEIF(bc_z%end == -2) THEN
                        
                        DO j = 1, buff_size
                            
                            ! Density or the partial densities and the momentum
                            ! components in x- and y-directions
                            DO i = 1, mom_idx%beg+1
                                q_cons_vf(i)%sf(:,:,p+j) = &
                                          q_cons_vf(i)%sf(:,:,(p+1)-j)
                            END DO
                            
                            ! z-component of momentum
                            q_cons_vf(mom_idx%end)%sf(:,:,p+j) = &
                                     -q_cons_vf(mom_idx%end)%sf(:,:,(p+1)-j)
                            
                            ! Energy and advection equation(s) variable(s)
                            DO i = E_idx, sys_size
                                q_cons_vf(i)%sf(:,:,p+j) = &
                                          q_cons_vf(i)%sf(:,:,(p+1)-j)
                            END DO
                            
                        END DO
                        
                    ! Perodic BC at the end
                    ELSEIF(bc_z%end == -1) THEN
                        
                        DO j = 1, buff_size
                            DO i = 1, sys_size
                                q_cons_vf(i)%sf(:,:,p+j) = &
                                            q_cons_vf(i)%sf(:,:,j-1)
                            END DO
                        END DO
                        
                    ! Processor BC at the end
                    ELSE
                        
                        CALL s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf,&
                                                                     'end', 'z')
                        
                    END IF
                    
                END IF
                
            END IF
            
            ! END: Populating Buffer Regions in the z-direction ================
            
            
        END SUBROUTINE s_populate_conservative_variables_buffer_regions ! ------
        
        
        
        
        !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module       
        SUBROUTINE s_initialize_data_input_module() ! -----------------------------

            INTEGER :: i !< Generic loop iterator
            
            
            ! Allocating the parts of the conservative and primitive variables
            ! that do not require the direct knowledge of the dimensionality of
            ! the simulation
            ALLOCATE(q_cons_vf(1:sys_size))
            ALLOCATE(q_prim_vf(1:sys_size))
            
            
            ! Allocating the parts of the conservative and primitive variables
            ! that do require the direct knowledge of the dimensionality of the
            ! simulation
            
            ! Simulation is at least 2D
            IF(n > 0) THEN
                
                ! Simulation is 3D
                IF(p > 0) THEN
                    
                    DO i = 1, sys_size
                        ALLOCATE(q_cons_vf(i)%sf( -buff_size:m+buff_size, &
                                                  -buff_size:n+buff_size, &
                                                  -buff_size:p+buff_size ))
                        ALLOCATE(q_prim_vf(i)%sf( -buff_size:m+buff_size, &
                                                  -buff_size:n+buff_size, &
                                                  -buff_size:p+buff_size ))
                    END DO
                    
                ! Simulation is 2D
                ELSE
                    
                    DO i = 1, sys_size
                        ALLOCATE(q_cons_vf(i)%sf( -buff_size:m+buff_size, &
                                                  -buff_size:n+buff_size, &
                                                           0:0           ))
                        ALLOCATE(q_prim_vf(i)%sf( -buff_size:m+buff_size, &
                                                  -buff_size:n+buff_size, &
                                                           0:0           ))
                    END DO
                    
                END IF
                
            ! Simulation is 1D
            ELSE
                
                DO i = 1, sys_size
                    ALLOCATE(q_cons_vf(i)%sf( -buff_size:m+buff_size, &
                                                       0:0          , &
                                                       0:0           ))
                    ALLOCATE(q_prim_vf(i)%sf( -buff_size:m+buff_size, &
                                                       0:0          , &
                                                       0:0           ))
                END DO
                
            END IF
            
            IF (parallel_io .NEQV. .TRUE.) THEN
                s_read_data_files => s_read_serial_data_files
            ELSE
                s_read_data_files => s_read_parallel_data_files
            END IF
            
        END SUBROUTINE s_initialize_data_input_module ! ---------------------------
        
        
        
        
        
        !> Deallocation procedures for the module
        SUBROUTINE s_finalize_data_input_module() ! --------------------------


            INTEGER :: i !< Generic loop iterator
            
            
            ! Deallocating the conservative and primitive variables
            DO i = 1, sys_size
                DEALLOCATE(q_cons_vf(i)%sf)
                DEALLOCATE(q_prim_vf(i)%sf)
            END DO
            
            DEALLOCATE(q_cons_vf)
            DEALLOCATE(q_prim_vf)
            
            s_read_data_files => NULL()
            
        END SUBROUTINE s_finalize_data_input_module ! ------------------------
        
        
        
        
        
END MODULE m_data_input
