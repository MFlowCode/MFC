! pMFC_v3.0 - Post-Process code: m_particles_output.f90
! Description: The module contains subroutines used to output particle variables
! Author: Kazuki Maeda
! Date: 01/01/17


MODULE m_particles_output

  ! Dependencies =============================================================
  USE m_global_parameters
  USE m_derived_types
  USE m_mpi_proxy
  USE m_particles_types
  ! ==========================================================================
 

  IMPLICIT NONE

         REAL(KIND(0.D0)), ALLOCATABLE, DIMENSION(:,:)     :: MPI_IO_DATA_particle
 
  CONTAINS


  SUBROUTINE write_particle_parallel (t_step,subflag)

    CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir
    CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc
    LOGICAL :: dir_check
    INTEGER :: id, nparticles

    INTEGER, INTENT(IN) :: t_step
    LOGICAL, OPTIONAL :: subflag
   
#ifdef MFC_MPI
    REAL(KIND(0.D0)), DIMENSION(20)   :: inputvals
    REAL(KIND(0.D0)) :: id_real, time_real
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER(KIND=MPI_OFFSET_KIND) :: disp
    INTEGER :: view
 
    TYPE(particledata)    , POINTER :: particleinfo
    TYPE(particleListinfo), POINTER :: particleListaux

    INTEGER, DIMENSION(3)   :: cell
    LOGICAL                 :: indomain, particle_file, particle_data, file_exist

    INTEGER, DIMENSION(2) :: gsizes, lsizes, start_idx_part
    INTEGER :: ifile, ireq, ierr, data_size, tot_data
    INTEGER :: i


    IF(num_procs >1 ) THEN
    PRINT '(A)', ' particle output does not support MPI. Exiting ...'
       CALL s_mpi_abort
    END IF
    
    IF(second_dir.AND. .NOT.PRESENT(subflag)) THEN
        WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io' , t_step, '.dat'
        file_loc = TRIM(case_dir) // '/restart_data2' // TRIM(mpiiofs) // TRIM(file_loc)
        INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
    ELSE
        IF(PRESENT(subflag)) THEN
           WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io_sub' , t_step, '.dat'
        ELSE
           WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io' , t_step, '.dat'
        END IF
        file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
        INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
    END IF

    IF (file_exist) THEN
      IF(proc_rank .EQ. 0) THEN
        OPEN(9, FILE = TRIM(file_loc), FORM = 'unformatted', STATUS = 'unknown')
          READ(9) tot_data, time_real
        CLOSE(9)
      END IF
    ELSE
      PRINT '(A)', TRIM(file_loc) // ' is missing. Exiting ...'
      CALL s_mpi_abort
    END IF
    
    CALL MPI_BCAST(tot_data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(time_real, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    gsizes(1)=tot_data
    gsizes(2)=21
    lsizes(1)=tot_data
    lsizes(2)=21
    start_idx_part(1)=0
    start_idx_part(2)=0

    CALL MPI_TYPE_CREATE_SUBARRAY(2,gsizes,lsizes,start_idx_part,&
           MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,view,ierr)
    CALL MPI_TYPE_COMMIT(view,ierr)

    ! Open the file to read particle variables
    IF(second_dir.AND. .NOT.PRESENT(subflag)) THEN
       WRITE(file_loc, '(A,I0,A)') 'particle' , t_step, '.dat'
       file_loc = TRIM(case_dir) // '/restart_data2' // TRIM(mpiiofs) // TRIM(file_loc)
       INQUIRE(FILE = TRIM(file_loc),EXIST = particle_file)
    ELSE
       IF(PRESENT(subflag)) THEN
          WRITE(file_loc, '(A,I0,A)') 'particle_sub' , t_step, '.dat'
       ELSE
          WRITE(file_loc, '(A,I0,A)') 'particle' , t_step, '.dat'
       END IF
       file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
       INQUIRE(FILE = TRIM(file_loc),EXIST = particle_file)
    END IF

    IF (particle_file) THEN

       CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY, &
                          mpi_info_int,ifile,ierr)
   
       disp = 0d0
       CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,view, &
                                          'native',mpi_info_null,ierr)
   
       ALLOCATE(MPI_IO_DATA_particle(tot_data,1:21))
   
       CALL MPI_FILE_READ_ALL(ifile,MPI_IO_DATA_particle,21*tot_data, &
                                    MPI_DOUBLE_PRECISION,status,ierr)

       IF(PRESENT(subflag)) THEN
          WRITE(file_loc, '(A,I0,A)') 'particles_data_sub' , t_step, '.dat'
       ELSE
          WRITE(file_loc, '(A,I0,A)') 'particles_data' , t_step, '.dat'
       END IF
       file_loc = TRIM(case_dir) // '/particles_data/' // TRIM(file_loc)

       IF(proc_rank == 0) OPEN(unit=29,file=file_loc,form='formatted')


       DO i = 1, tot_data
   
          id = INT(MPI_IO_DATA_particle(i,1))
          inputvals(1:20) = MPI_IO_DATA_particle(i,2:21)
      
        
           IF(id.GT.0) THEN
              WRITE(29, 6) INT(id), inputvals(1), inputvals(2), &
              inputvals(3), inputvals(4),inputvals(5), inputvals(6), inputvals(7), &
              inputvals(8),inputvals(9),inputvals(10), inputvals(11),&
              inputvals(12),inputvals(13),inputvals(14), inputvals(15),&
              inputvals(16),inputvals(17), inputvals(18), inputvals(19),&
              inputvals(20),time_real
              6 FORMAT(I6,21(1x,E15.7))
           END IF

       END DO

       CLOSE(29)

       DEALLOCATE(MPI_IO_DATA_particle)


    ENDIF
  
    CALL MPI_FILE_CLOSE(ifile,ierr)

#endif

  END SUBROUTINE write_particle_parallel


  SUBROUTINE get_tot_step (t_step,tot_step)

    CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir
    CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc
    LOGICAL :: dir_check
    INTEGER :: id, nparticles

    INTEGER, INTENT(IN) :: t_step
    INTEGER, INTENT(INOUT) :: tot_step

#ifdef MFC_MPI   
    REAL(KIND(0.D0)), DIMENSION(20)   :: inputvals
    REAL(KIND(0.D0)) :: id_real, time_real
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER(KIND=MPI_OFFSET_KIND) :: disp
    INTEGER :: view
 
    TYPE(particledata)    , POINTER :: particleinfo
    TYPE(particleListinfo), POINTER :: particleListaux

    INTEGER, DIMENSION(3)   :: cell
    LOGICAL                 :: indomain, particle_file, particle_data, file_exist

    INTEGER, DIMENSION(2) :: gsizes, lsizes, start_idx_part
    INTEGER :: ifile, ireq, ierr, data_size, tot_data
    INTEGER :: i

    INTEGER :: dmN1
    REAL(KIND(0.D0)) :: dmR1, dmR2

    IF(num_procs >1 ) THEN
    PRINT '(A)', ' particle output does not support MPI. Exiting ...'
       CALL s_mpi_abort
    END IF
    
    WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io' , t_step, '.dat'
    file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
    INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
 
    IF (file_exist) THEN
      IF(proc_rank .EQ. 0) THEN
        OPEN(9, FILE = TRIM(file_loc), FORM = 'unformatted', STATUS = 'unknown')
          READ(9) dmN1, dmR1, dmR2, tot_step
        CLOSE(9)
      END IF
    ELSE
      PRINT '(A)', TRIM(file_loc) // ' is missing. Exiting ...'
      CALL s_mpi_abort
    END IF

#endif

 END SUBROUTINE get_tot_step


END MODULE m_particles_output
