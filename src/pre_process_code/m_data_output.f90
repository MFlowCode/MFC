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
!! @file m_data_output.f90
!! @brief Contains module m_data_output
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module takes care of writing the grid and initial condition
!!              data files into the "0" time-step directory located in the folder
!!              associated with the rank of the local processor, which is a sub-
!!              directory of the case folder specified by the user in the input
!!              file pre_process.inp.
MODULE m_data_output
    
    ! Dependencies =============================================================
    USE m_derived_types         !< Definitions of the derived types
    
    USE m_global_parameters     !< Global parameters for the code
    
    USE m_mpi_proxy             !< Message passing interface (MPI) module proxy

    USE mpi                     !< Message passing interface (MPI) module
    
    USE m_compile_specific

    USE m_variables_conversion
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    
    PRIVATE; PUBLIC :: s_write_serial_data_files, &
                       s_write_parallel_data_files, &
                       s_write_data_files, &
                       s_initialize_data_output_module, &
                       s_finalize_data_output_module

    ABSTRACT INTERFACE ! ===================================================

        !>  Interface for the conservative data
        !! @param q_cons_vf The conservative variables
        SUBROUTINE s_write_abstract_data_files(q_cons_vf)

            IMPORT :: scalar_field, sys_size

            ! Conservative variables
            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf


        END SUBROUTINE s_write_abstract_data_files ! -------------------
    END INTERFACE ! ========================================================



    CHARACTER(LEN = path_len + 2*name_len), PRIVATE :: t_step_dir !<
    !! Time-step folder into which grid and initial condition data will be placed
    
    CHARACTER(LEN = path_len + 2*name_len), PUBLIC :: restart_dir !< 
    !! Restart data folder

    PROCEDURE(s_write_abstract_data_files), POINTER :: s_write_data_files => NULL()
    
    CONTAINS
        
        !>  Writes grid and initial condition data files to the "0"
        !!  time-step directory in the local processor rank folder
        !! @param q_cons_vf The conservative variables
        SUBROUTINE s_write_serial_data_files(q_cons_vf) ! -----------
            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf

            LOGICAL :: file_exist !< checks if file exists
            
            CHARACTER(LEN=15) :: FMT

            CHARACTER(LEN = &
            INT(FLOOR(LOG10(REAL(sys_size, KIND(0d0))))) + 1) :: file_num !< Used to store 
            !! the number, in character form, of the currently
            !! manipulated conservative variable data file

            CHARACTER(LEN = LEN_TRIM(t_step_dir) + name_len) :: file_loc !<
            !! Generic string used to store the address of a particular file
            
            INTEGER :: i,j,k,l !< Generic loop iterator
            INTEGER :: t_step
           
            REAL(KIND(0d0)), DIMENSION(nb) :: nRtmp         !< Temporary bubble concentration
            REAL(KIND(0d0)) :: nbub                         !< Temporary bubble number density
            REAL(KIND(0d0)) :: gamma, lit_gamma, pi_inf     !< Temporary EOS params
            REAL(KIND(0d0)) :: rho                          !< Temporary density

            t_step = 0

            ! Outputting the Locations of the Cell-boundaries ==================
            
            ! x-coordinate direction
            file_loc = TRIM(t_step_dir) // '/x_cb.dat'
            OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', STATUS = 'new')
            WRITE(1) x_cb
            CLOSE(1)
            
            ! y- and z-coordinate directions
            IF(n > 0) THEN
                ! y-coordinate direction
                file_loc = TRIM(t_step_dir) // '/y_cb.dat'
                OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                        STATUS = 'new')
                WRITE(1) y_cb
                CLOSE(1)
                
                ! z-coordinate direction
                IF(p > 0) THEN
                    file_loc = TRIM(t_step_dir) // '/z_cb.dat'
                    OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                            STATUS = 'new')
                    WRITE(1) z_cb
                    CLOSE(1)
                END IF
            END IF
            ! ==================================================================
            
            
            ! Outputting Conservative Variables ================================
            DO i = 1, sys_size
                WRITE(file_num, '(I0)') i
                file_loc = TRIM(t_step_dir) // '/q_cons_vf' // TRIM(file_num) &
                       // '.dat'
                OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                    STATUS = 'new')
                WRITE(1) q_cons_vf(i)%sf
                CLOSE(1)
            END DO
            ! ==================================================================
    
            gamma = fluid_pp(1)%gamma
            lit_gamma = 1d0/fluid_pp(1)%gamma + 1d0
            pi_inf = fluid_pp(1)%pi_inf

            IF (precision==1) THEN
                FMT="(2F30.7)"
            ELSE
                FMT="(2F40.14)"
            END IF
            
            WRITE(t_step_dir,'(A,I0,A,I0)') TRIM(case_dir) // '/D'
            file_loc = TRIM(t_step_dir) // '/.'
        
            INQUIRE( FILE = TRIM(file_loc), EXIST = file_exist )
        
            IF(.NOT.file_exist) CALL SYSTEM('mkdir -p ' // TRIM(t_step_dir))

            !1D
            IF (n ==0 .AND. p ==0) THEN
                IF (model_eqns == 2) THEN
                    DO i = 1, sys_size
                        WRITE(file_loc,'(A,I0,A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/prim.', i, '.', proc_rank, '.', t_step,'.dat'

                        OPEN(2,FILE= TRIM(file_loc) )
                            DO j=0,m
                                CALL s_convert_to_mixture_variables( q_cons_vf, j, 0, 0, rho, gamma, pi_inf)

                                lit_gamma = 1d0/gamma + 1d0
                                
                                IF ( ((i.ge.cont_idx%beg) .AND. (i.le.cont_idx%end))    &
                                                          .OR.                          &
                                     ((i.ge.adv_idx%beg)  .AND. (i.le.adv_idx%end))     &
                                                         ) THEN 
                                    WRITE(2,FMT) x_cb(j),q_cons_vf(i)%sf(j,0,0)
                                ELSE IF (i.eq.mom_idx%beg) THEN !u
                                    WRITE(2,FMT) x_cb(j),q_cons_vf(mom_idx%beg)%sf(j,0,0)/rho
                                ELSE IF (i.eq.stress_idx%beg) THEN !tau_e
                                    WRITE(2,FMT) x_cb(j),q_cons_vf(stress_idx%beg)%sf(j,0,0)/rho
                                ELSE IF (i.eq.E_idx) THEN !p
                                    IF (model_eqns == 4) THEN
                                        !Tait pressure from density
                                        WRITE(2,FMT) x_cb(j), &
                                            (pref + pi_inf) * (                     &
                                            ( q_cons_vf(1)%sf(j,0,0)/               &
                                            (rhoref*(1.d0-q_cons_vf(4)%sf(j,0,0)))  & 
                                            ) ** lit_gamma )                        &
                                            - pi_inf
                                    ELSE IF (model_eqns == 2 .AND. (bubbles .NEQV. .TRUE.)) THEN
                                        !Stiffened gas pressure from energy
                                        WRITE(2,FMT) x_cb(j), &
                                            (                                       & 
                                            q_cons_vf(E_idx)%sf(j,0,0)  -            &
                                            0.5d0*(q_cons_vf(mom_idx%beg)%sf(j,0,0)**2.d0)/rho - &
                                            pi_inf &
                                            ) / gamma
                                    ELSE
                                        !Stiffened gas pressure from energy with bubbles
                                        WRITE(2,FMT) x_cb(j), &
                                            (                                       & 
                                            (q_cons_vf(E_idx)%sf(j,0,0)  -            &
                                            0.5d0*(q_cons_vf(mom_idx%beg)%sf(j,0,0)**2.d0)/rho) / &
                                            (1.d0 - q_cons_vf(alf_idx)%sf(j,0,0)) - &
                                            pi_inf &
                                            ) / gamma
                                    END IF
                                ELSE IF ((i .GE. bub_idx%beg) .AND. (i .LE. bub_idx%end) .AND. bubbles) THEN
                                    DO k = 1,nb
                                        nRtmp(k) = q_cons_vf(bub_idx%rs(k))%sf(j,0,0)
                                    END DO
                                    CALL s_comp_n_from_cons( q_cons_vf(alf_idx)%sf(j,0,0), nRtmp, nbub)                                
                                    
                                    WRITE(2,FMT) x_cb(j),q_cons_vf(i)%sf(j,0,0)/nbub
                                END IF
                            END DO
                        CLOSE(2)
                    END DO
                END IF

                DO i = 1, sys_size    
                    WRITE(file_loc,'(A,I0,A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step,'.dat'

                    OPEN(2,FILE= TRIM(file_loc) )
                        DO j=0,m
                            WRITE(2,FMT) x_cb(j),q_cons_vf(i)%sf(j,0,0)
                        END DO
                    CLOSE(2)
                END DO
            END IF


            IF (precision==1) THEN
                FMT="(3F30.7)"
            ELSE
                FMT="(3F40.14)"
            END IF

            ! 2D
            IF ( (n>0) .AND. (p==0) ) THEN
                DO i = 1,sys_size
                    WRITE(file_loc,'(A,I0,A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step,'.dat'
                    OPEN(2,FILE= TRIM(file_loc) )
                        DO j=0,m
                        DO k=0,n
                            WRITE(2,FMT) x_cb(j),y_cb(k), q_cons_vf(i)%sf(j,k,0)
                        END DO
                        WRITE(2,*)
                        END DO
                    CLOSE(2)
                END DO
            END IF


            ! 3D
            IF ( p > 0) THEN
                DO i = 1,sys_size
                    WRITE(file_loc,'(A,I0,A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step,'.dat'
                    OPEN(2,FILE= TRIM(file_loc) )
                        DO j=0,m
                        DO k=0,n
                        DO l=0,p
                            WRITE(2,FMT) x_cb(j),y_cb(k),z_cb(l), q_cons_vf(i)%sf(j,k,l)
                        END DO
                        WRITE(2,*)
                        END DO
                        WRITE(2,*)
                        END DO
                    CLOSE(2)
                END DO
            END IF

        END SUBROUTINE s_write_serial_data_files ! ------------------------------------



        !> Writes grid and initial condition data files in parallel to the "0"
        !!  time-step directory in the local processor rank folder
        !! @param q_cons_vf The conservative variables
        SUBROUTINE s_write_parallel_data_files(q_cons_vf) ! --

            ! Conservative variables
            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf

            INTEGER :: ifile, ierr, data_size
            INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
            INTEGER(KIND=MPI_OFFSET_KIND) :: disp
            INTEGER(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: MOK

            CHARACTER(LEN=path_len + 2*name_len) :: file_loc 
            LOGICAL :: file_exist

            ! Generic loop iterator
            INTEGER :: i

            ! Initialize MPI data I/O
            CALL s_initialize_mpi_data(q_cons_vf)

            ! Open the file to write all flow variables
            WRITE(file_loc, '(A)') '0.dat'
            file_loc = TRIM(restart_dir) // TRIM(mpiiofs) // TRIM(file_loc)
            INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
            IF (file_exist .AND. proc_rank == 0) THEN
                CALL MPI_FILE_DELETE(file_loc,mpi_info_int,ierr)
            END IF
            CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                        mpi_info_int,ifile,ierr)

            ! Size of local arrays
            data_size = (m+1)*(n+1)*(p+1)

            ! Resize some integers so MPI can write even the biggest files
            m_MOK     = INT(m_glb+1,  MPI_OFFSET_KIND)
            n_MOK     = INT(n_glb+1,  MPI_OFFSET_KIND)
            p_MOK     = INT(p_glb+1,  MPI_OFFSET_KIND)
            WP_MOK    = INT(8d0,      MPI_OFFSET_KIND)
            MOK       = INT(1d0,      MPI_OFFSET_KIND)
            str_MOK   = INT(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = INT(sys_size, MPI_OFFSET_KIND)
            
            ! Write the data for each variable
            IF (bubbles) THEN
                DO i = 1, sys_size! adv_idx%end
                var_MOK = INT(i, MPI_OFFSET_KIND)

                ! Initial displacement to skip at beginning of file
                disp = m_MOK*MAX(MOK,n_MOK)*MAX(MOK,p_MOK)*WP_MOK*(var_MOK-1)

                CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i), &
                            'native',mpi_info_int,ierr)
                CALL MPI_FILE_WRITE_ALL(ifile,MPI_IO_DATA%var(i)%sf,data_size, &
                            MPI_DOUBLE_PRECISION,status,ierr)
                END DO
            ELSE
                DO i = 1, adv_idx%end
                var_MOK = INT(i, MPI_OFFSET_KIND)

                ! Initial displacement to skip at beginning of file
                disp = m_MOK*MAX(MOK,n_MOK)*MAX(MOK,p_MOK)*WP_MOK*(var_MOK-1)

                CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i), &
                            'native',mpi_info_int,ierr)
                CALL MPI_FILE_WRITE_ALL(ifile,MPI_IO_DATA%var(i)%sf,data_size, &
                            MPI_DOUBLE_PRECISION,status,ierr)
                END DO
            END IF

            CALL MPI_FILE_CLOSE(ifile,ierr)

        END SUBROUTINE s_write_parallel_data_files ! ---------------------------

        
        !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
        SUBROUTINE s_initialize_data_output_module() ! ----------------------------
            ! Generic string used to store the address of a particular file
            CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: file_loc
            
            ! Generic logical used to check the existence of directories
            LOGICAL :: dir_check
            
            
            IF (parallel_io .NEQV. .TRUE.) THEN
                ! Setting the address of the time-step directory
                WRITE(t_step_dir, '(A,I0,A)') '/p_all/p', proc_rank, '/0'
                t_step_dir = TRIM(case_dir) // TRIM(t_step_dir)
                
                
                ! Checking the existence of the time-step directory, removing it, if
                ! it exists, and creating a new copy. Note that if preexisting grid
                ! and/or initial condition data are to be read in from the very same
                ! location, then the above described steps are not executed here but
                ! rather in the module m_start_up.f90.
                IF(old_grid .NEQV. .TRUE.) THEN
                    
                    file_loc = TRIM(t_step_dir) // '/'
                    
                    CALL my_inquire(file_loc,dir_check)
        
                    IF(dir_check) CALL SYSTEM('rm -rf ' // TRIM(t_step_dir))
                    CALL SYSTEM('mkdir -p ' // TRIM(t_step_dir))
                    
                END IF
             
                s_write_data_files => s_write_serial_data_files
            ELSE
                WRITE(restart_dir, '(A)') '/restart_data'
                restart_dir = TRIM(case_dir) // TRIM(restart_dir)

                IF ((old_grid .NEQV. .TRUE.) .AND. (proc_rank == 0)) THEN

                    file_loc = TRIM(restart_dir) // '/'
                    CALL my_inquire(file_loc,dir_check)

                    IF (dir_check) CALL SYSTEM('rm -rf ' // TRIM(restart_dir))
                    CALL SYSTEM('mkdir -p ' // TRIM(restart_dir))
                END IF

                CALL s_mpi_barrier()

                s_write_data_files => s_write_parallel_data_files

            END IF
            
        END SUBROUTINE s_initialize_data_output_module ! --------------------------
        
        
        
        !> Resets s_write_data_files pointer
        SUBROUTINE s_finalize_data_output_module() ! ---------------------------

            s_write_data_files => NULL()

        END SUBROUTINE s_finalize_data_output_module ! -------------------------
        
        
        
        
        
END MODULE m_data_output
