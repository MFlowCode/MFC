! pMFC_v3.0 - Simulation code: m_particles_output.f90
! Description: The module contains subroutines used to output particle variables
! Author: Kazuki Maeda
! Date: 01/01/17


MODULE m_particles_output

  USE m_global_parameters
  USE m_particles_types
  USE m_time_tmp
  USE m_mpi_particles

  IMPLICIT NONE

  !private; public :: write_samples

  REAL(KIND(0.D0)) :: Rmax=0.
  REAL(KIND(0.D0)) :: Rmin=1000000

  TYPE :: fileop
     CHARACTER(50) :: dir,prefix,ext
     INTEGER :: funit
  END TYPE fileop

  TYPE (fileop) :: particlefile     = fileop("./D/"  ,"particles",".bin",11), &
                   particlestatfile = fileop("./D/"  ,"particle_stats",".dat",12), &
                   samplesfile      = fileop("./D/"  ,"samples",".dat",13), &
                   voidfractionfile = fileop("./D/"  ,"voidfraction",".dat",14)
CONTAINS

  SUBROUTINE write_particles (qtime)

    TYPE(particlenode),POINTER     :: particle
    REAL(KIND(0.D0))               :: qtime, pinf
    INTEGER                        :: i
    INTEGER, DIMENSION(3)          :: cell
  
    IF (qtime.eq.0.d0) THEN
       CALL openfile(particlefile,append = .FALSE. )
    ELSE
       CALL openfile(particlefile,append = .TRUE. )
    END IF
  
    ! Cycle through list
      particle  => particlesubList%List%next
      DO WHILE(Associated(particle))
  
      WRITE(particlefile%funit) qtime, particle%data%id, &
      (particle%data%tmp%x(i), i=1,3), &
      particle%data%tmp%mv, particle%data%tmp%mv/(particle%data%tmp%mv + particle%data%mg), &
      particle%data%tmp%y(1), particle%data%tmp%y(2), particle%data%tmp%p
  
      particle => particle%next
      ENDDO
  
    CALL closefile(particlefile)

  END SUBROUTINE write_particles


  SUBROUTINE write_samples (qtime, q, q_prim)

    USE m_rhs
  
    REAL(KIND(0.D0))                   :: qtime
    INTEGER                            :: i
    TYPE(scalar_field), DIMENSION(sys_size)      :: q
    TYPE(scalar_field), DIMENSION(sys_size)      :: q_prim
   
   !type(vector_field) :: gm_alpha_qp
   !TYPE(bounds_info) :: ix,iy,iz
   !TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: q_particle
   

    IF (avgdensFlag) THEN
      !CALL s_convert_conservative_to_primitive_variables(q, q_prim, gm_alpha_qp(0,0,0)%vf, ix, iy, iz, q_particle(1))
      CALL s_convert_conservative_to_primitive_variables(q, q_prim, gm_alpha_qp%vf, ix, iy, iz, q_particle(1))
    ELSE
      !CALL s_convert_conservative_to_primitive_variables(q, q_prim,gm_alpha_qp(0,0,0)%vf, ix, iy, iz)
      CALL s_convert_conservative_to_primitive_variables(q, q_prim,gm_alpha_qp%vf, ix, iy, iz)
    ENDIF
  
    IF (qtime.eq.0.d0) THEN
       CALL openfile(samplesfile,append = .FALSE. )
    ELSE
       CALL openfile(samplesfile,append = .TRUE. )
    END IF
  
    DO i=1, nsamples,1
      WRITE(samplesfile%funit,*) prb(i)%id, qtime, q_prim(E_idx)%sf(prb(i)%ip(1),prb(i)%ip(2),prb(i)%ip(3))
    ENDDO 
  
    CALL closefile(samplesfile)
  
    CALL write_yz_averages(qtime, q, q_prim)

  END SUBROUTINE write_samples


  SUBROUTINE write_yz_averages(qtime, q, q_prim)

    USE m_rhs
    USE m_mpi_particles
  
    REAL(KIND(0.D0))                   :: qtime,intgrl,vol,totvol
    INTEGER                            :: i,j,k,l
    INTEGER, DIMENSION(3)              :: cell
    TYPE(scalar_field), DIMENSION(sys_size)        :: q
    TYPE(scalar_field), DIMENSION(sys_size)        :: q_prim
   
    IF (proc_rank.EQ.0) THEN
      IF (qtime.eq.0.d0) THEN
        OPEN(unit=2,file='integral_yz.dat', status = "unknown")
      ELSE
        OPEN(unit=2,file='integral_yz.dat', position = "append", status = "old")
      ENDIF
    ENDIF
 
    DO i=1,totsamples
      j=1
      intgrl = 0.0d0
      totvol = 0.0d0
      IF (nsamples.GT.0) THEN
      DO WHILE ((j.LE.nsamples).AND.(prb(j)%id.LT.i))
        j=j+1
      ENDDO
      IF (prb(j)%id.EQ.i) THEN
        vol    = 0.0d0
        cell(1) = prb(j)%ip(1)
        DO k=0,p
          DO l=0,n
            cell(2) = l
            cell(3) = k
            CALL get_char_vol (cell, vol)
            intgrl = intgrl + q_prim(E_idx)%sf(prb(j)%ip(1),l,k)*vol
            totvol = totvol + vol
          ENDDO
        ENDDO
      ENDIF
      ENDIF
  
      IF (num_procs > 1 ) THEN
        CALL get_sum (intgrl)
        CALL get_sum ( totvol )
      ENDIF
 
      IF (proc_rank.EQ.0) THEN
        !print*, qtime, proc_rank, intgrl/totvol, i, totsamples, nsamples
        IF (totvol.GT.0.0d0) WRITE(2,*) qtime, intgrl/totvol
      ENDIF
    ENDDO
  
    CLOSE(2)

  END SUBROUTINE write_yz_averages


  SUBROUTINE write_void_evol (qtime)

    USE m_rhs
    USE m_mpi_particles
  
    REAL(KIND(0.D0))                   :: qtime, voidmax,voidavg,vol,volcell,voltot
    INTEGER                            :: i,j,k
    INTEGER, DIMENSION(3)              :: cell
    LOGICAL                            :: prevfile
    !TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: q_particle

    IF (proc_rank.EQ.0) THEN
      IF (qtime.eq.0.d0) THEN
        INQUIRE (file='voidfraction_0.dat'   ,exist=prevfile )
        IF (prevfile) CALL SYSTEM('cp -r voidfraction_0.dat  voidfraction_0.old')
        CALL openfile(voidfractionfile,append = .FALSE. )
      ELSE
        CALL openfile(voidfractionfile,append = .TRUE. )
      END IF
    ENDIF
  
    voidmax = 0.0d0
    voidavg = 0.0d0
    vol     = 0.0d0
  
    DO  i=0,m
      DO j=0,n
        DO k=0,p
          cell(1)=i
          cell(2)=j
          cell(3)=k
          voidmax = max(voidmax, 1 - q_particle(1)%sf(i,j,k))
          CALL get_char_vol (cell, volcell)
          IF ((1. - q_particle(1)%sf(i,j,k)).GT.5e-11) THEN
            voidavg = voidavg + (1 - q_particle(1)%sf(i,j,k))*volcell
            vol     = vol + volcell
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  
    IF (num_procs > 1 ) THEN 
      CALL get_max ( voidmax )
      CALL get_sum ( vol )
      CALL get_sum ( voidavg )
    ENDIF
     
    voltot=voidavg 
    ! This voidavg value does not reflect the real void fraction in the cloud
    ! since the cell which does not have bubbles are not accounted
    IF (vol.gt.0.) voidavg = voidavg/vol
  
    IF (proc_rank.EQ.0) THEN
      WRITE(voidfractionfile%funit,*) qtime, voidavg, voidmax, voltot
      CALL closefile(voidfractionfile)
    ENDIF

  END SUBROUTINE write_void_evol


  SUBROUTINE write_restart_particles (t_step)

  ! Location of time-step folder corresponding to time-step, t_step
  CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir

  ! Generic string used to store the address of a particular file
  CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc

  ! Generic logical used for purpose of asserting whether a particular
  ! directory is or is not located in the designated location
  LOGICAL :: dir_check

  TYPE(particlenode),POINTER         :: particle
  INTEGER                            :: i, t_step

  ! Setting up the address of the time-directory
  IF(parallel_io .NEQV. .TRUE.) THEN
    WRITE(t_step_dir, '(A,I0,A,I0)') '/p', proc_rank, '/', t_step
  ELSE
    WRITE(t_step_dir, '(A,I0,A,I0)') 'restart_data'
  END IF

  t_step_dir = TRIM(case_dir) // TRIM(t_step_dir)

  file_loc = TRIM(t_step_dir) // '/.'

  IF(parallel_io .NEQV. .TRUE.) THEN
    WRITE(file_loc, '(I0)') proc_rank
    file_loc = TRIM(t_step_dir) // '/particles' // TRIM(file_loc) &
    // '.dat'
  ELSE
    WRITE(file_loc, '(I0)') proc_rank
    file_loc = TRIM(t_step_dir) // '/particles' // TRIM(t_step_dir) &
    // '.dat'
  END IF

  OPEN(2, FILE = TRIM(file_loc), FORM = 'unformatted', STATUS = 'unknown')

  ! Cycle through list
  particle  => particlesubList%List%next
  DO WHILE(Associated(particle))

    WRITE(2) particle%data%id, &
    (particle%data%x(i), i=1,3), &
    (particle%data%xprev(i), i=1,3), &
    (particle%data%u(i), i=1,3), &
    (particle%data%y(i), i=1,2), &
    particle%data%R0, particle%data%Rmax, particle%data%Rmin, &
    particle%data%dphidt, particle%data%p, particle%data%mv, &
    particle%data%mg, particle%data%betaT, particle%data%betaC
    particle => particle%next
  ENDDO
  
  CLOSE(2)

  END SUBROUTINE write_restart_particles


  SUBROUTINE write_restart_particles_parallel (t_step)

      ! Location of time-step folder corresponding to time-step, t_step
      CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir
    
      ! Generic string used to store the address of a particular file
      CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc
    
      ! Generic logical used for purpose of asserting whether a particular
      ! directory is or is not located in the designated location
      LOGICAL :: dir_check
      LOGICAL :: file_exist
    
      TYPE(particlenode), POINTER        :: particle
      INTEGER                            :: i, t_step
      INTEGER                            :: nparticles, tot_part, tot_part_wrtn, npart_wrtn
 
#ifdef MFC_MPI
      ! For Parallel I/O
      INTEGER :: ifile, ireq, ierr, data_size
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp
      INTEGER :: view
      INTEGER, DIMENSION(2) :: gsizes, lsizes, start_idx_part
      INTEGER, DIMENSION(num_procs) :: part_order, part_ord_mpi
    
      ! Number of particles in the node (fixme: better to broadcast just once)
    
      nparticles = 0d0
      IF(particlesubList%nb.NE.0) THEN
    
        particle  => particlesubList%List%next
        DO WHILE(Associated(particle))   
        
          IF(particle_in_domain_physical(particle%data%x(1:3))) THEN   
    
             nparticles = nparticles + 1
     
          END IF
     
          particle => particle%next 
     
        END DO
     
      END IF
       
      ! Total number of particles
      CALL MPI_ALLREDUCE(nparticles, tot_part, 1, MPI_INTEGER, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
    
      ! Total number of particles written so far
      CALL MPI_ALLREDUCE(npart_wrtn, tot_part_wrtn, 1, MPI_INTEGER, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
    
      ! lsizes(1) = MAX(1,nparticles)
      lsizes(1) = MAX(1,nparticles)
      lsizes(2) = 21
    
      ! If the partcle number is zero, put 1 since MPI cannot deal with writing
      ! zero particle
      part_order(:) = 1
      part_order(proc_rank+1) = MAX(1,nparticles)
    
      CALL MPI_ALLREDUCE(part_order, part_ord_mpi, num_procs, MPI_INTEGER, &
                            MPI_MAX, MPI_COMM_WORLD, ierr)
     
      gsizes(1) = SUM(part_ord_mpi(1:num_procs))
      gsizes(2) = 21
    
      start_idx_part(1) = SUM(part_ord_mpi(1:proc_rank+1))-part_ord_mpi(proc_rank+1)!-MAX(1,nparticles)
      start_idx_part(2) = 0
    
      WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io' , t_step, '.dat'
      file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
      INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
      IF (file_exist .AND. proc_rank == 0) THEN
          CALL MPI_FILE_DELETE(file_loc,mpi_info_int,ierr)
      END IF
    
      ! Writing down the total number of particles
      IF(proc_rank .EQ. 0) THEN
        OPEN(9, FILE = TRIM(file_loc), FORM = 'unformatted', STATUS = 'unknown')
          WRITE(9) gsizes(1), time_real, dt_next_inp, tot_step
        CLOSE(9)
      END IF
     
      CALL MPI_TYPE_CREATE_SUBARRAY(2,gsizes,lsizes,start_idx_part,&
               MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,view,ierr)
      CALL MPI_TYPE_COMMIT(view,ierr)
    
      ALLOCATE(MPI_IO_DATA_particle(1:MAX(1,nparticles),1:21))
       
      ! Open the file to write all flow variables
      WRITE(file_loc, '(A,I0,A)') 'particle' , t_step, '.dat'
      file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
      INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
      IF (file_exist .AND. proc_rank == 0) THEN
          CALL MPI_FILE_DELETE(file_loc,mpi_info_int,ierr)
      END IF
     
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE), &
                         mpi_info_int,ifile,ierr)
       
      disp = 0d0
    
      CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,view, &
                                         'native',mpi_info_null,ierr)
 
      ! Cycle through list
      i=1
  
      IF(nparticles.EQ.0) THEN
        MPI_IO_DATA_particle(1,1:21) = 0d0
      ELSE
  
        particle  => particlesubList%List%next
        DO WHILE(Associated(particle))
         
          IF(particle_in_domain_physical(particle%data%x(1:3))) THEN
  
             MPI_IO_DATA_particle(i,1) = REAL(particle%data%id)
             MPI_IO_DATA_particle(i,2:4) = particle%data%x(1:3)
             MPI_IO_DATA_particle(i,5:7) = particle%data%xprev(1:3)
             MPI_IO_DATA_particle(i,8:10) = particle%data%u(1:3)
             MPI_IO_DATA_particle(i,11:12) = particle%data%y(1:2)
             MPI_IO_DATA_particle(i,13) = particle%data%R0
             MPI_IO_DATA_particle(i,14) = particle%data%Rmax
             MPI_IO_DATA_particle(i,15) = particle%data%Rmin
             MPI_IO_DATA_particle(i,16) = particle%data%dphidt
             MPI_IO_DATA_particle(i,17) = particle%data%p
             MPI_IO_DATA_particle(i,18) = particle%data%mv
             MPI_IO_DATA_particle(i,19) = particle%data%mg
             MPI_IO_DATA_particle(i,20) = particle%data%betaT
             MPI_IO_DATA_particle(i,21) = particle%data%betaC
  
             i = i + 1
  
          END IF
  
          particle => particle%next 
  
       END DO
   
     END IF
 
     CALL MPI_FILE_WRITE_ALL(ifile,MPI_IO_DATA_particle,21*MAX(1,nparticles), &
                                   MPI_DOUBLE_PRECISION,status,ierr)
     
     CALL MPI_FILE_CLOSE(ifile,ierr)
   
     DEALLOCATE(MPI_IO_DATA_particle)

#endif

  END SUBROUTINE write_restart_particles_parallel


  SUBROUTINE write_restart_particles_parallel2 (t_step,subout)

      ! Location of time-step folder corresponding to time-step, t_step
      CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir
    
      ! Generic string used to store the address of a particular file
      CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc
    
      ! Generic logical used for purpose of asserting whether a particular
      ! directory is or is not located in the designated location
      LOGICAL :: dir_check
      LOGICAL :: file_exist
      
      ! flag for the sub-time step output of particle data
      LOGICAL, OPTIONAL :: subout
    
      TYPE(particlenode), POINTER        :: particle
      INTEGER                            :: i, t_step
      INTEGER                            :: nparticles, tot_part, tot_part_wrtn, npart_wrtn

#ifdef MFC_MPI

      ! For Parallel I/O
      INTEGER :: ifile, ireq, ierr, data_size
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp
      INTEGER :: view
      INTEGER, DIMENSION(2) :: gsizes, lsizes, start_idx_part
      INTEGER, DIMENSION(num_procs) :: part_order, part_ord_mpi
    
      ! Number of particles in the node (fixme: better to broadcast just once)
    
      nparticles = 0d0
      IF(particlesubList%nb.NE.0) THEN
    
        particle  => particlesubList%List%next
        DO WHILE(Associated(particle))   
        
          IF(particle_in_domain_physical(particle%data%x(1:3))) THEN   
    
             nparticles = nparticles + 1
     
          END IF
     
          particle => particle%next 
     
        END DO
     
      END IF
 
      ! Total number of particles
      CALL MPI_ALLREDUCE(nparticles, tot_part, 1, MPI_INTEGER, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
    
      ! Total number of particles written so far
      CALL MPI_ALLREDUCE(npart_wrtn, tot_part_wrtn, 1, MPI_INTEGER, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
    
      lsizes(1) = MAX(1,nparticles)
      lsizes(2) = 21
    
      ! If the partcle number is zero, put 1 since MPI cannot deal with writing
      ! zero particle
      part_order(:) = 1
      part_order(proc_rank+1) = MAX(1,nparticles)
    
      CALL MPI_ALLREDUCE(part_order, part_ord_mpi, num_procs, MPI_INTEGER, &
                            MPI_MAX, MPI_COMM_WORLD, ierr)
     
      gsizes(1) = SUM(part_ord_mpi(1:num_procs))
      gsizes(2) = 21
    
      start_idx_part(1) = SUM(part_ord_mpi(1:proc_rank+1))-part_ord_mpi(proc_rank+1)!-MAX(1,nparticles)
      start_idx_part(2) = 0
    
      ! Writing down the total number of particles
      IF(PRESENT(subout)) THEN
         WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io_sub' , tot_step, '.dat'
         file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
      ELSE
         WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io' , t_step, '.dat'
         file_loc = TRIM(case_dir) // '/restart_data2' // TRIM(mpiiofs) // TRIM(file_loc)
      END IF
      INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
      IF (file_exist .AND. proc_rank == 0) THEN
          CALL MPI_FILE_DELETE(file_loc,mpi_info_int,ierr)
      END IF
    
      IF(proc_rank .EQ. 0) THEN
         OPEN(9, FILE = TRIM(file_loc), FORM = 'unformatted', STATUS = 'unknown')
            WRITE(9) gsizes(1), time_real, dt_next_inp, tot_step
         CLOSE(9)
      END IF   
     
      CALL MPI_TYPE_CREATE_SUBARRAY(2,gsizes,lsizes,start_idx_part,&
               MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,view,ierr)
      CALL MPI_TYPE_COMMIT(view,ierr)
    
      ALLOCATE(MPI_IO_DATA_particle(1:MAX(1,nparticles),1:21))

      ! Open the file to write all flow variables
      IF(PRESENT(subout)) THEN
         WRITE(file_loc, '(A,I0,A)') 'particle_sub' , tot_step, '.dat'
         file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
      ELSE
         WRITE(file_loc, '(A,I0,A)') 'particle' , t_step, '.dat'
         file_loc = TRIM(case_dir) // '/restart_data2' // TRIM(mpiiofs) // TRIM(file_loc)
      END IF
      INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
      IF (file_exist .AND. proc_rank == 0) THEN
          CALL MPI_FILE_DELETE(file_loc,mpi_info_int,ierr)
      END IF
     
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE), &
                         mpi_info_int,ifile,ierr)
     
      disp = 0d0
    
      ! Define the view for each variable
      CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,view, &
                                         'native',mpi_info_null,ierr)

      ! Cycle through list
      i=1
  
      IF(nparticles.EQ.0) THEN
        MPI_IO_DATA_particle(1,1:21) = 0d0
      ELSE
  
        particle  => particlesubList%List%next
        DO WHILE(Associated(particle))
        
          IF(particle_in_domain_physical(particle%data%x(1:3))) THEN
  
             MPI_IO_DATA_particle(i,1) = REAL(particle%data%id)
             MPI_IO_DATA_particle(i,2:4) = particle%data%x(1:3)
             MPI_IO_DATA_particle(i,5:7) = particle%data%xprev(1:3)
             MPI_IO_DATA_particle(i,8:10) = particle%data%u(1:3)
             MPI_IO_DATA_particle(i,11:12) = particle%data%y(1:2)
             MPI_IO_DATA_particle(i,13) = particle%data%R0
             MPI_IO_DATA_particle(i,14) = particle%data%Rmax
             MPI_IO_DATA_particle(i,15) = particle%data%Rmin
             MPI_IO_DATA_particle(i,16) = particle%data%dphidt
             MPI_IO_DATA_particle(i,17) = particle%data%p
             MPI_IO_DATA_particle(i,18) = particle%data%mv
             MPI_IO_DATA_particle(i,19) = particle%data%mg
             MPI_IO_DATA_particle(i,20) = particle%data%betaT
             MPI_IO_DATA_particle(i,21) = particle%data%betaC
  
             i = i + 1
  
          END IF
  
          particle => particle%next 
  
       END DO
   
     END IF
 
     CALL MPI_FILE_WRITE_ALL(ifile,MPI_IO_DATA_particle,21*MAX(1,nparticles), &
                                   MPI_DOUBLE_PRECISION,status,ierr)
   
 
     CALL MPI_FILE_CLOSE(ifile,ierr)
  
     DEALLOCATE(MPI_IO_DATA_particle)

#endif

  END SUBROUTINE write_restart_particles_parallel2


  SUBROUTINE write_restart_particles_time(t_step)

  ! Location of time-step folder corresponding to time-step, t_step
  CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir

  ! Generic string used to store the address of a particular file
  CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc

  ! Generic logical used for purpose of asserting whether a particular
  ! directory is or is not located in the designated location
  LOGICAL :: dir_check
  LOGICAL :: file_exist

  INTEGER :: i, t_step
  INTEGER :: ierr

#ifdef MFC_MPI 
  WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io_time' , t_step, '.dat'
  file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
  INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
  IF (file_exist .AND. proc_rank == 0) THEN
      CALL MPI_FILE_DELETE(file_loc,mpi_info_int,ierr)
  END IF

  ! Writing down the total number of particles
  IF(proc_rank .EQ. 0) THEN
    OPEN(19, FILE = TRIM(file_loc), FORM = 'unformatted', STATUS = 'unknown')
      WRITE(19) time_real
    CLOSE(19)
  END IF
#endif
  END SUBROUTINE write_restart_particles_time


  SUBROUTINE read_restart_particles_time(t_step)

  ! Location of time-step folder corresponding to time-step, t_step
  CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir

  ! Generic string used to store the address of a particular file
  CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc

  ! Generic logical used for purpose of asserting whether a particular
  ! directory is or is not located in the designated location
  LOGICAL :: dir_check
  LOGICAL :: file_exist

  INTEGER :: i, t_step
  INTEGER :: ierr

#ifdef MFC_MPI 
  WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io' , t_step, '.dat'
  file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
  INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
  IF (file_exist .AND. proc_rank == 0) THEN
      CALL MPI_FILE_DELETE(file_loc,mpi_info_int,ierr)
  END IF

  ! Writing down the total number of particles
  IF(proc_rank .EQ. 0) THEN
    OPEN(19, FILE = TRIM(file_loc), FORM = 'unformatted', STATUS = 'unknown')
      WRITE(19) time_real
    CLOSE(19)
  END IF

#endif

  END SUBROUTINE read_restart_particles_time


  SUBROUTINE particle_stats ()

  TYPE(particlenode),POINTER                :: particle

  ! Cycle through list
    particle  => particlesubList%List%next
    DO WHILE(Associated(particle))
        Rmax=max(Rmax,particle%data%y(1)/particle%data%R0)
        Rmin=min(Rmin,particle%data%y(1)/particle%data%R0)
        particle%data%Rmax = max( particle%data%Rmax,particle%data%y(1)/particle%data%R0)
        particle%data%Rmin = min( particle%data%Rmin,particle%data%y(1)/particle%data%R0)
        particle => particle%next
    ENDDO

  END SUBROUTINE particle_stats


  SUBROUTINE write_particle_stats

  TYPE(particlenode),POINTER  :: particle
  INTEGER                     :: i

     CALL openfile(particlestatfile,append = .FALSE. )
     WRITE(particlestatfile%funit,*) '#',Rmax,Rmin
     WRITE(particlestatfile%funit,*) '#################'
     
         particle  => particlesubList%List%next
         DO WHILE(Associated(particle))
            WRITE(particlestatfile%funit,*) (particle%data%x(i), i=1,3),particle%data%Rmax,particle%data%Rmin
            particle => particle%next
         ENDDO
     
     CALL closefile(particlestatfile)

  END SUBROUTINE write_particle_stats
  

  SUBROUTINE openfile( thefile, version, append )

     INTEGER, OPTIONAL :: version
     LOGICAL, OPTIONAL :: append
     INTEGER           :: nochar
     TYPE (fileop)     :: thefile
     CHARACTER(6)      :: cversion, openversion
     CHARACTER(4)      :: crank
     CHARACTER(50)     :: theformat, rewindstatus, newness
     LOGICAL           :: do_warning

     CALL verno ( 4, proc_rank, crank, nochar )
     crank = crank(4-nochar+1:4)

    IF (thefile%ext(1:LEN_TRIM(thefile%ext)).eq.".bin" ) THEN
        theformat = "unformatted"
     ELSE
        theformat = "formatted"
     END IF

     IF (PRESENT(append)) THEN
       IF (append) THEN
        rewindstatus = "append"
       ELSE
        rewindstatus = "rewind"
       END IF
     ELSE
        rewindstatus = "rewind"
     END IF

     IF (PRESENT(version)) THEN
         newness ='unknown'
         CALL verno( 6, version, cversion, nochar )
         openversion = cversion(6-nochar+1:6)
     ELSE
         newness = 'unknown'
         openversion = ''
         INQUIRE( file=thefile%dir(1:LEN_TRIM(thefile%dir))       // &
                       thefile%prefix(1:LEN_TRIM(thefile%prefix)) // &
                       openversion(1:LEN_TRIM(openversion))       // &
                       '_'                                        // &
                       crank(1:LEN_TRIM(crank))                   // &
                       thefile%ext(1:LEN_TRIM(thefile%ext)),         &
                  exist=do_warning )
         IF ((do_warning).and.(rewindstatus.eq."rewind")) THEN
            if (proc_rank == 0) THEN  !print once
               WRITE(*,*) "Warning: the file you are opening"
               WRITE(*,*) thefile%dir(1:LEN_TRIM(thefile%dir))       // &
                          thefile%prefix(1:LEN_TRIM(thefile%prefix)) // &
                          openversion(1:LEN_TRIM(openversion))       // &
                          thefile%ext(1:LEN_TRIM(thefile%ext))
               WRITE(*,*) "with position=rewind already existed."
            endif
         END IF
     END IF

     OPEN(unit=thefile%funit,                                &
          file=thefile%dir(1:LEN_TRIM(thefile%dir))       // &
               thefile%prefix(1:LEN_TRIM(thefile%prefix)) // &
               openversion(1:LEN_TRIM(openversion))       // &
               '_'                                        // &
               crank(1:LEN_TRIM(crank))                   // &
               thefile%ext(1:LEN_TRIM(thefile%ext)),         &
          form = theformat(1:LEN_TRIM(theformat)),           &
          position = rewindstatus(1:LEN_TRIM(rewindstatus)), &
          status = newness(1:LEN_TRIM(newness)), RECL=5000000 )
 
  END SUBROUTINE openfile


  SUBROUTINE closefile( thefile )

    IMPLICIT NONE
    TYPE (fileop) :: thefile

    CLOSE( unit=thefile%funit )

  END SUBROUTINE closefile


  SUBROUTINE verno (n, index, str, len)

!
!  returns a character string for the integer "index" into "str"
!  it also returns len, the length of nonzero digits in str
!

    IMPLICIT NONE

    INTEGER      :: index, i, t, n, len
    CHARACTER(n) :: str

    DO i=1,n
        t = index/(10**(n-i))
        t = t - (t/10)*10
        str(i:i) = CHAR(48+t)
    END DO

    len = n

    lenloop : DO i=1,n
                  IF (str(i:i).eq.CHAR(48)) THEN
                     len = len - 1
                  ELSE
                     EXIT lenloop
                  END IF
              END DO lenloop

    IF (len.eq.0) len = 1

  END SUBROUTINE verno

END MODULE m_particles_output
