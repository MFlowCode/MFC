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
!! @file m_grid.f90
!! @brief Contains module m_grid
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief  This module takes care of creating the rectilinear grid on which
!!              the data for the initial condition will be laid out and on which
!!              the simulation will eventually be computed. The grid may either
!!              be uniform or non-uniform. Non-uniform grids are generated using
!!              the hyperbolic tangent function, see Johnsen (2007) for details.
!!              Alternatively to synthesizing a new grid, the user may select to
!!              read in a preexisting one. This is carried out through the module
!!              m_start_up.f90. In such a case, the responsibility of this module
!!              becomes only to allocate/deallocate the necessary grid variables
!!              for the cell-centers and cell-boundaries locations.
MODULE m_grid
    
    
    ! Dependencies =============================================================
    USE m_derived_types         ! Definitions of the derived types
    
    USE m_global_parameters     ! Global parameters for the code
    
    USE m_mpi_proxy             ! Message passing interface (MPI) module proxy
    
    USE mpi                     ! Message passing interface (MPI) module
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    
    PRIVATE; PUBLIC :: s_initialize_grid_module, &
                       s_generate_grid, &
                       s_generate_serial_grid, &
                       s_generate_parallel_grid, &
                       s_finalize_grid_module

    ABSTRACT INTERFACE ! ===================================================

        SUBROUTINE s_generate_abstract_grid(dflt_int) ! ------------------------

            INTEGER, INTENT(IN) :: dflt_int

        END SUBROUTINE s_generate_abstract_grid ! ----------------------

    END INTERFACE ! ========================================================

    PROCEDURE(s_generate_abstract_grid), POINTER :: s_generate_grid => NULL()
    
    CONTAINS
        
        
        !> The following subroutine generates either a uniform or
        !!              non-uniform rectilinear grid in serial, defined by the parameters
        !!              inputted by the user. The grid information is stored in
        !!              the grid variables containing coordinates of the cell-
        !!              centers and cell-boundaries.
        !! @param dflt_int is the default null integer used for checking initializationss  
        SUBROUTINE s_generate_serial_grid(dflt_int) ! -----------------------------------------
            
            INTEGER, INTENT(IN) :: dflt_int
            
            ! Generic loop iterator
            INTEGER :: i,j              !< generic loop operatorss
            REAL(KIND(0d0)) :: length   !< domain lengths
            
            ! Grid Generation in the x-direction ===============================
            dx = (x_domain%end - x_domain%beg) / REAL(m+1, KIND(0d0))
            
            DO i = 0, m
                x_cc( i ) = x_domain%beg + 5d-1*dx*REAL(2*i+1, KIND(0d0))
                x_cb(i-1) = x_domain%beg +      dx*REAL(  i  , KIND(0d0))
            END DO
            
            x_cb(m) = x_domain%end
            
            IF(stretch_x) THEN
                
                length = abs(x_cb(m)-x_cb(0))
                x_cb = x_cb / length
                x_a = x_a / length
                x_b = x_b / length

                DO j = 1, loops_x
                    DO i = -1, m
                        x_cb(i) = x_cb(i) / a_x * &
                                ( a_x + LOG(COSH(  a_x*(x_cb(i) - x_a)  )) &
                                      + LOG(COSH(  a_x*(x_cb(i) - x_b)  )) &
                                  - 2d0*LOG(COSH( 0.5d0*a_x*(x_b - x_a) )) )
                    END DO
                END DO
                x_cb = x_cb * length
                
                x_cc = (x_cb(0:m) + x_cb(-1:m-1))/2d0
                
                dx = MINVAL(x_cb(0:m) - x_cb(-1:m-1))
                print*, 'Stretched grid: min/max x grid: ', minval(x_cc(:)), maxval(x_cc(:))
                IF(num_procs > 1) CALL s_mpi_reduce_min(dx)
                
            END IF
            ! ==================================================================
            
            
            ! Grid Generation in the y-direction ===============================
            IF(n == 0) RETURN
            
            IF (grid_geometry == 2.AND.y_domain%beg.eq.0.0d0) THEN
            !IF (grid_geometry == 2) THEN

               dy = (y_domain%end - y_domain%beg) / REAL(2*n+1, KIND(0d0))
               
               y_cc(0)  = y_domain%beg + 5d-1*dy
               y_cb(-1) = y_domain%beg

               DO i = 1, n
                   y_cc( i ) = y_domain%beg + 2d0*dy*REAL(  i  , KIND(0d0))
                   y_cb(i-1) = y_domain%beg +     dy*REAL(2*i-1, KIND(0d0))
               END DO

            ELSE
             
               dy = (y_domain%end - y_domain%beg) / REAL(n+1, KIND(0d0))
            
               DO i = 0, n
                   y_cc( i ) = y_domain%beg + 5d-1*dy*REAL(2*i+1, KIND(0d0))
                   y_cb(i-1) = y_domain%beg +      dy*REAL(  i  , KIND(0d0))
               END DO

            END IF

            y_cb(n) = y_domain%end
            
            IF(stretch_y) THEN
                
                DO j = 1, loops_y
                    DO i = -1, n
                        y_cb(i) = y_cb(i) / a_y * &
                                ( a_y + LOG(COSH(  a_y*(y_cb(i) - y_a)  )) &
                                      + LOG(COSH(  a_y*(y_cb(i) - y_b)  )) &
                                  - 2d0*LOG(COSH( 0.5d0*a_y*(y_b - y_a) )) )
                    END DO
                END DO
                
                y_cc = (y_cb(0:n) + y_cb(-1:n-1))/2d0
                
                dy = MINVAL(y_cb(0:n) - y_cb(-1:n-1))
                
                IF(num_procs > 1) CALL s_mpi_reduce_min(dy)
                
            END IF
            ! ==================================================================
            
            
            ! Grid Generation in the z-direction ===============================
            IF(p == 0) RETURN
            
            dz = (z_domain%end - z_domain%beg) / REAL(p+1, KIND(0d0))
            
            DO i = 0, p
                z_cc( i ) = z_domain%beg + 5d-1*dz*REAL(2*i+1, KIND(0d0))
                z_cb(i-1) = z_domain%beg +      dz*REAL(  i  , KIND(0d0))
            END DO
            
            z_cb(p) = z_domain%end
            
            IF(stretch_z) THEN
                
                DO j = 1, loops_z
                    DO i = -1, p
                        z_cb(i) = z_cb(i) / a_z * &
                                ( a_z + LOG(COSH(  a_z*(z_cb(i) - z_a)  )) &
                                      + LOG(COSH(  a_z*(z_cb(i) - z_b)  )) &
                                  - 2d0*LOG(COSH( 0.5d0*a_z*(z_b - z_a) )) )
                    END DO
                END DO
                
                z_cc = (z_cb(0:p) + z_cb(-1:p-1))/2d0
                
                dz = MINVAL(z_cb(0:p) - z_cb(-1:p-1))
                
                IF(num_procs > 1) CALL s_mpi_reduce_min(dz)
                
            END IF
            ! ==================================================================
            
            
        END SUBROUTINE s_generate_serial_grid ! ---------------------------------------
        
        
        
        !> The following subroutine generates either a uniform or
        !!              non-uniform rectilinear grid in parallel, defined by the parameters
        !!              inputted by the user. The grid information is stored in
        !!              the grid variables containing coordinates of the cell-
        !!              centers and cell-boundaries.
        !! @param dflt_int is she default null integer used for checking initializationss         
        SUBROUTINE s_generate_parallel_grid(dflt_int) !-------------------------

            INTEGER, INTENT(IN) :: dflt_int

            ! Locations of cell boundaries
            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: x_cb_glb, y_cb_glb, z_cb_glb !<
            !! Locations of cell boundaries


            CHARACTER(LEN = path_len + name_len) :: file_loc !< 
            !! Generic string used to store the address of a file

            INTEGER :: ifile, ierr, data_size
            INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

            INTEGER :: i,j !< Generic loop integers

            ALLOCATE(x_cb_glb(-1:m_glb))
            ALLOCATE(y_cb_glb(-1:n_glb))
            ALLOCATE(z_cb_glb(-1:p_glb))

            ! Grid generation in the x-direction
            dx = (x_domain%end - x_domain%beg) / REAL(m_glb+1, KIND(0d0))
            DO i = 0, m_glb
                x_cb_glb(i-1) = x_domain%beg + dx*REAL(i,KIND(0d0))
            END DO
            x_cb_glb(m_glb) = x_domain%end
            IF (stretch_x) THEN
                DO j = 1, loops_x
                    DO i = -1, m_glb
                        x_cb_glb(i) = x_cb_glb(i) / a_x * &
                            (a_x + LOG(COSH( a_x*(x_cb_glb(i) - x_a))) & 
                                 + LOG(COSH( a_x*(x_cb_glb(i) - x_b))) & 
                             - 2d0*LOG(COSH(0.5d0*a_x*(x_b - x_a))) )
                    END DO
                END DO
            END IF

            ! Grid generation in the y-direction
            IF (n_glb > 0) THEN

                IF (grid_geometry == 2.AND.y_domain%beg.eq.0.0d0) THEN
                    dy = (y_domain%end - y_domain%beg) / REAL(2*n_glb+1, KIND(0d0))
                    y_cb_glb(-1) = y_domain%beg
                    DO i = 1, n_glb
                        y_cb_glb(i-1) = y_domain%beg +  dy*REAL(2*i-1, KIND(0d0))
                    END DO
                ELSE
                    dy = (y_domain%end - y_domain%beg) / REAL(n_glb+1, KIND(0d0))
                    DO i = 0, n_glb
                        y_cb_glb(i-1) = y_domain%beg + dy*REAL(i, KIND(0d0))
                    END DO
                END IF
                y_cb_glb(n_glb) = y_domain%end
                IF (stretch_y) THEN
                    DO j = 1, loops_y
                        DO i = -1, n_glb
                            y_cb_glb(i) = y_cb_glb(i) / a_y * &
                                (a_y + LOG(COSH( a_y*(y_cb_glb(i) - y_a))) & 
                                     + LOG(COSH( a_y*(y_cb_glb(i) - y_b))) & 
                                 - 2d0*LOG(COSH(0.5d0*a_y*(y_b - y_a))) )
                        END DO
                    END DO
                END IF

                ! Grid generation in the z-direction
                IF (p_glb > 0) THEN
                    dz = (z_domain%end - z_domain%beg) / REAL(p_glb+1, KIND(0d0))
                    DO i = 0, p_glb
                        z_cb_glb(i-1) = z_domain%beg + dz*REAL(i,KIND(0d0))
                    END DO
                    z_cb_glb(p_glb) = z_domain%end
                    IF (stretch_z) THEN
                        DO j = 1, loops_z
                            DO i = -1, p_glb
                                z_cb_glb(i) = z_cb_glb(i) / a_z * &
                                    (a_z + LOG(COSH( a_z*(z_cb_glb(i) - z_a))) & 
                                         + LOG(COSH( a_z*(z_cb_glb(i) - z_b))) & 
                                     - 2d0*LOG(COSH(0.5d0*a_z*(z_b - z_a))) )
                            END DO
                        END DO
                    END IF
                END IF
            END IF

            ! Write cell boundary locations to grid data files
            file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'x_cb.dat'
            data_size = m_glb+2
            CALL MPI_FILE_OPEN(MPI_COMM_SELF,file_loc,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                        mpi_info_int,ifile,ierr)
            CALL MPI_FILE_WRITE(ifile,x_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
            CALL MPI_FILE_CLOSE(ifile,ierr)

            IF (n > 0) THEN
                file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'y_cb.dat'
                data_size = n_glb+2
                CALL MPI_FILE_OPEN(MPI_COMM_SELF,file_loc,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                            mpi_info_int,ifile,ierr)
                CALL MPI_FILE_WRITE(ifile,y_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                CALL MPI_FILE_CLOSE(ifile,ierr)

                IF (p > 0) THEN
                    file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'z_cb.dat'
                    data_size = p_glb+2
                    CALL MPI_FILE_OPEN(MPI_COMM_SELF,file_loc,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                                mpi_info_int,ifile,ierr)
                    CALL MPI_FILE_WRITE(ifile,z_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                    CALL MPI_FILE_CLOSE(ifile,ierr)
                END IF
            END IF

            DEALLOCATE(x_cb_glb, y_cb_glb, z_cb_glb)

        END SUBROUTINE s_generate_parallel_grid ! ------------------------------

        !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
        SUBROUTINE s_initialize_grid_module() ! -----------------------------------
            
            IF (parallel_io .NEQV. .TRUE.) THEN
                s_generate_grid => s_generate_serial_grid
            ELSE
                s_generate_grid => s_generate_parallel_grid
            END IF
            
        END SUBROUTINE s_initialize_grid_module ! ---------------------------------
        
        !> Deallocation procedures for the module        
        SUBROUTINE s_finalize_grid_module() ! --------------------------------
            
            s_generate_grid => NULL()   
            
        END SUBROUTINE s_finalize_grid_module ! ------------------------------
        
END MODULE m_grid
