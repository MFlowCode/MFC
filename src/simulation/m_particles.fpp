! pMFC_v3.0 - Simulation code: m_mpi_particles.f90
! Description: The module contains subroutines for simulation of particle
! variables
! Author: Kazuki Maeda
! Date: 01/01/17


MODULE m_particles

    ! Dependencies =============================================================
    USE m_global_parameters
    USE m_derived_types
    USE m_particles_types
    USE m_rhs
    USE m_mpi_particles
    USE m_time_tmp
    use m_mpi_common
    use m_data_output
    use m_particles_output
    ! ==========================================================================


    IMPLICIT NONE

    CONTAINS

    SUBROUTINE s_initialize_lagrangian_solver(q_cons_vf, q_prim_vf, dtnext, dtdid, time_prev, tavg)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_prim_vf
        
        real(kind(0.d0)), intent(INOUT) :: dtnext, dtdid, time_prev, tavg

        real(kind(0.d0)) :: dtoutput, tend


        ! Time-step iterator to the first time-step
        dt0    = dt
        dtoutput = dt * t_step_save
        time_real = t_step_start * dt
        tend      = time_real + ( t_step_stop - t_step_start) * dt
        dt_next_inp = dt0

        ! Initializing particle module
        IF (particleflag) THEN
            CALL s_prepare_particle_module()
            CALL s_initialize_particles(q_cons_vf, q_prim_vf)
            IF(t_step_start.EQ.0) THEN
                IF (parallel_io .NEQV. .TRUE.) THEN
                    CALL s_write_initial_conc_serial(q_particle(1))
                ELSE
                    CALL s_write_data_files(q_cons_vf, q_prim_vf, t_step_start, q_particle(1))
                   !IF(second_dir) CALL s_write_parallel_data_files2(q_cons_ts(1)%vf, q_prim_vf, t_step_start, q_particle(1))
                END IF

                IF (parallel_io .NEQV. .TRUE.) THEN
                    CALL write_restart_particles (t_step_start)
                ELSE
                    CALL write_restart_particles_parallel (t_step_start)
                   !IF(second_dir) THEN
                       !CALL write_restart_particles_parallel2 (t_step_start)
                   !END IF
                   !CALL write_restart_particles_parallel2 (t_step_start,.TRUE.)
                END IF

                IF(avgdensFlag) CALL write_void_evol (time_real)
            END IF
            CALL s_populate_variables_buffers(q_cons_vf,q_particle)
        ELSE
            CALL s_populate_variables_buffers(q_cons_vf)
            !IF (We_size > 0 .AND. (We_riemann_flux .OR. We_rhs_flux)) THEN
                !CALL s_account_for_capillary_potential_energy(q_cons_ts(1)%vf)
            !END IF
        END IF

        IF(particleflag) THEN
            ! Recovering dtnext
            IF(t_step_start .NE. 0.0d0) THEN
                dtnext = dt_next_inp
            ELSE IF(dtmaxpart.GT.0.0d0) THEN
                dtnext = min(dt0, dtmaxpart)
            ELSE
                dtnext = dt0
            END IF
        ELSE 
        dtnext = dt0
        END IF

        ! Initialize samples particles
        INQUIRE (file='samples.inp', exist = do_sampling)
        IF(do_sampling .AND. particleflag) THEN
            CALL s_read_samples()
            IF(t_step_start.EQ.0) CALL write_samples (time_real, q_cons_vf, q_prim_vf)
        END IF

        IF(proc_rank == 0) THEN
            IF(particleflag.and.run_time_info) THEN
                OPEN(unit=2,file='times.dat')
                CLOSE(2)
            END IF

           !IF(second_dir) THEN
           !   OPEN(unit=4,file='second_dir_time.dat')
           !   CLOSE(4)
           !   IF(t_step_start.eq.0) THEN
           !     OPEN(unit=4, file='second_dir_time.dat', position = "append", status = "old")
           !       WRITE(4,*) time_real
           !     CLOSE(4)
           !   END IF
           !END IF

        ENDIF
        tavg = 0.d0

    END SUBROUTINE s_initialize_lagrangian_solver


    SUBROUTINE s_lagrangian_run_time_info(q_cons_vf, q_prim_vf, time_real, t_step, dtnext, tavg)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        real(kind(0d0)), intent(INOUT):: time_real, dtnext, tavg
        integer, intent(IN):: t_step

        IF(avgdensFlag) CALL write_void_evol(time_real)
        CALL system_clock(cpu_end)
        IF(proc_rank == 0) THEN
            if (t_step > 1) then
                tavg = tavg + REAL(cpu_end-cpu_start)/cpu_rate
                OPEN(unit=2, file='times.dat', position = "append", status = "old")
                WRITE(2,*) t_step, time_real,REAL(cpu_end-cpu_start)/cpu_rate,dtnext,tavg/real(t_step-1)
                CLOSE(2)
            end if
        END IF
        cpu_start = cpu_end
        IF(do_sampling) CALL write_samples(time_real, q_cons_vf, q_prim_vf)

    END SUBROUTINE s_lagrangian_run_time_info


    SUBROUTINE s_particles_nondimensionalize_inputs()
     
        !> Non-dimensionalize inputs
        if (particleflag) then
            epsilonb = 1
            Runiv = Runiv*Tini/(csonref**2)
            pvap = pvap/(rholiqref*csonref**2)
            cpgas = cpgas*Tini/(csonref**2)
            cpvapor = cpvapor*Tini/(csonref**2)
            kgas = kgas*Tini/(Lref*csonref**3*rholiqref)
            kvapor = kvapor*Tini/(Lref*csonref**3*rholiqref)
            diffcoefvap = diffcoefvap/(Lref*csonref)
            sigmabubble = sigmabubble/(rholiqref*csonref**2*Lref)
            viscref = viscref/(rholiqref*csonref*Lref)
            IF (.NOT.avgdensFlag) solverapproach = 0

        !> default values
        else
            avgdensFlag = .FALSE.
            particleoutFlag = .FALSE.
            particlestatFlag = .FALSE.
            RPflag = .FALSE.
            clusterflag = dflt_int
            stillparticlesflag= .FALSE.
            heatflag= dflt_int
            massflag= dflt_int
            csonref = dflt_real
            rholiqref = dflt_real
            Lref = dflt_real
            Tini = dflt_real
            Runiv = dflt_real
            gammagas = dflt_real
            gammavapor = dflt_real
            pvap = dflt_real
            cpgas = dflt_real
            cpvapor = dflt_real
            kgas = dflt_real
            kvapor = dflt_real
            MWgas = dflt_real
            MWvap = dflt_real
            diffcoefvap = dflt_real
            sigmabubble = dflt_real
            viscref = dflt_real
            RKeps = dflt_real
            ratiodt = dflt_int
            projectiontype = dflt_int
            smoothtype = dflt_int
            epsilonb= dflt_real
            coupledFlag = .FALSE.
            solverapproach = 2
            correctpresFlag = .FALSE.
            charwidth = dflt_real
            valmaxvoid = dflt_real
            dtmaxpart = dflt_real
            do_particles = .FALSE.
            bubblesources = .FALSE.
        end if

    END SUBROUTINE s_particles_nondimensionalize_inputs



  SUBROUTINE s_read_particle_input()

    !IMPLICIT NONE
    !OPEN(unit = 109, file = 'particles.inp', status='old')
    !READ(109,*) particleflag
    !READ(109,*) avgdensFlag
    !IF (.NOT.particleflag) avgdensFlag = .FALSE.
    !epsilonb = 1
    !READ(109,*) particleoutFlag
    !READ(109,*) particlestatFlag
    !READ(109,*) RPflag
    !READ(109,*) clusterflag
    !READ(109,*) stillparticlesflag
    !READ(109,*) heatflag
    !READ(109,*) massflag
    !READ(109,*) csonref
    !READ(109,*) rholiqref
    !READ(109,*) Lref
    !READ(109,*) Tini
    !READ(109,*) Runiv
    !Runiv = Runiv*Tini/(csonref**2)
    !READ(109,*) gammagas
    !READ(109,*) gammavapor
    !READ(109,*) pvap
    !pvap = pvap/(rholiqref*csonref**2)
    !READ(109,*) cpgas
    !cpgas = cpgas*Tini/(csonref**2)
    !READ(109,*) cpvapor
    !cpvapor = cpvapor*Tini/(csonref**2)
    !READ(109,*) kgas
    !kgas = kgas*Tini/(Lref*csonref**3*rholiqref)
    !READ(109,*) kvapor
    !kvapor = kvapor*Tini/(Lref*csonref**3*rholiqref) 
    !READ(109,*) MWgas
    !READ(109,*) MWvap
    !READ(109,*) diffcoefvap
    !diffcoefvap = diffcoefvap/(Lref*csonref)
    !READ(109,*) sigmabubble
    !sigmabubble = sigmabubble/(rholiqref*csonref**2*Lref)
    !READ(109,*) viscref
    !viscref = viscref/(rholiqref*csonref*Lref) 
    !READ(109,*) RKeps
    !READ(109,*) ratiodt
    !READ(109,*) projectiontype
    !READ(109,*) smoothtype
    !READ(109,*) epsilonb
    !READ(109,*) coupledFlag
    !READ(109,*) solverapproach
    !READ(109,*) correctpresFlag
    !READ(109,*) charwidth
    !READ(109, *, end = 110, err=110) valmaxvoid
    !READ(109, *, end = 110, err=110) dtmaxpart
!110  CLOSE(109)    

    !IF (.NOT.avgdensFlag) solverapproach = 0

  END SUBROUTINE s_read_particle_input



  SUBROUTINE s_read_samples()

    INTEGER :: i,j
    INTEGER, DIMENSION(3)            :: cell
    REAL(KIND(0.D0)), DIMENSION(3)   :: inputparticle
    LOGICAL  :: indomain

    nsamples=0
    cell = 0

    OPEN(unit=2,file='samples.inp',form='formatted')
112 READ(2,*, end=111) (inputparticle(i), i=1,3)
    indomain = particle_in_domain(inputparticle(1:3))
    IF (indomain) nsamples = nsamples+1
    GOTO 112
111 CLOSE(2)

    ALLOCATE(prb(nsamples))

    print *, 'samples in proc', proc_rank,':',nsamples
    IF (.NOT.particleflag) avgdensFlag = .FALSE.

    OPEN(unit=2,file='samples.inp',form='formatted')
    j=1
    totsamples=1
114  READ(2,*, end=113) (inputparticle(i), i=1,3)
    indomain = particle_in_domain(inputparticle(1:3))
    IF (indomain) THEN
      CALL locate_cell ( inputparticle, cell)
      prb(j)%id = totsamples
      prb(j)%ip = cell
      j=j+1
    ENDIF
    totsamples = totsamples+1
    GOTO 114

113 CONTINUE

  END SUBROUTINE s_read_samples                


  SUBROUTINE s_prepare_particle_module() ! ------------------------------------

    INTEGER :: i,imax

    IF((solverapproach.EQ.2).AND.avgdensFlag) THEN
      ! comp 1: (1 - beta)
      ! comp 2: dbetadt
      ! comp 3 - imax : auxiliary variables
      imax=4
      IF(clusterflag.GE.4) imax=10 !subgrid noise model
      bubblesources=.TRUE.
    ELSEIF ((solverapproach.EQ.0).AND.avgdensFlag) THEN
      imax = 3
      bubblesources=.FALSE.
    ELSE
      imax=2
      bubblesources=.FALSE.
    ENDIF

    ALLOCATE(q_particle(1:imax))
    DO i=1,imax
      IF(p > 0) THEN
        DIM=3
        ALLOCATE(q_particle(i)%sf(-buff_size : m+buff_size, &
        -buff_size : n+buff_size, &
        -buff_size : p+buff_size ))
      ELSE
        DIM=2
        ALLOCATE(q_particle(i)%sf(-buff_size : m+buff_size, &
        -buff_size : n+buff_size, &
        0 : 0             ))
      ENDIF
    ENDDO

    q_particle(1)%sf = 1.d0 !represents 1-beta

    DO i=2,imax
      q_particle(i)%sf = 0.d0
    ENDDO


  END SUBROUTINE s_prepare_particle_module



  SUBROUTINE s_deallocate_particles()

    INTEGER :: i,j,k,imax
    TYPE(particlenode),POINTER :: particle

    IF ((solverapproach.EQ.2).AND.avgdensFlag) THEN
      imax = 4
      IF(clusterflag.GE.4) imax=10 !subgrid noise model
    ELSEIF ((solverapproach.EQ.0).AND.avgdensFlag) THEN
      imax = 3
    ELSE
      imax = 2
    ENDIF

    DO i=1,imax
      DEALLOCATE(q_particle(i)%sf)
    ENDDO

    DEALLOCATE(q_particle)

    particle  => particlesubList%List%next
    DO WHILE(Associated(particle))
       CALL Remove_particle (particle,2)
    ENDDO

    DO k=-buff_size,p+buff_size
      DO j=-buff_size,n+buff_size
        DO i=-buff_size,m+buff_size
          DEALLOCATE ( qbl%fp(i,j,k)%List )
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE( qbl%fp )
    DEALLOCATE(cellwbList%List)
    DEALLOCATE(cellwbList)

  END SUBROUTINE s_deallocate_particles


  SUBROUTINE get_part_id
    id = id + 1 !do it in parallel?
  END SUBROUTINE get_part_id


  SUBROUTINE add_particle_rest (nparticles)

    CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir
    CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc
    LOGICAL :: dir_check
    INTEGER :: id, nparticles

    REAL(KIND(0.D0)), DIMENSION(20)   :: inputvals

    TYPE(particledata)    , POINTER :: particleinfo
    TYPE(particleListinfo), POINTER :: particleListaux

    INTEGER, DIMENSION(3)   :: cell
    LOGICAL                 :: indomain, particle_file

    ! Setting up the address of the time-directory
    WRITE(t_step_dir, '(A,I0,A,I0)') '/p', proc_rank, '/', t_step_start
    t_step_dir = TRIM(case_dir) // TRIM(t_step_dir)

    file_loc = TRIM(t_step_dir) // '/.'
    WRITE(file_loc, '(I0)') proc_rank
    file_loc = TRIM(t_step_dir) // '/particles' // TRIM(file_loc) // '.dat'

    INQUIRE (file= TRIM(file_loc) ,exist=particle_file   )
    IF (particle_file) THEN
    OPEN(2, FILE = TRIM(file_loc) ,FORM = "unformatted", position = "rewind" )

    921  READ(2,end=922) id, inputvals(1:20)

    indomain = particle_in_domain(inputvals(1:3))
    IF (indomain) THEN

      ALLOCATE(particleinfo)
      nparticles = nparticles + 1
      particleinfo%id        = id
      particleinfo%x(1:3)    = inputvals(1:3) 
      particleinfo%xprev(1:3)= inputvals(4:6) 
      particleinfo%u(1:3)    = inputvals(7:9) 
      particleinfo%y(1:2)    = inputvals(10:11) 
      particleinfo%R0        = inputvals(12) 
      particleinfo%Rmax      = inputvals(13) 
      particleinfo%Rmin      = inputvals(14) 
      particleinfo%dphidt    = inputvals(15) 
      particleinfo%p         = inputvals(16) 
      particleinfo%mv        = inputvals(17) 
      particleinfo%mg        = inputvals(18) 
      particleinfo%betaT     = inputvals(19) 
      particleinfo%betaC     = inputvals(20) 
      particleinfo%equilibrium = .FALSE.

      cell    = -buff_size

      CALL locate_cell ( particleinfo%x,  cell, particleinfo%tmp%s )

      particleListaux => qbl%fp(cell(1), cell(2), cell(3))
      IF (particleListaux%nb.EQ.0) CALL addtocell_list ( cell )
      CALL transfertotmp (particleinfo) 
      CALL addparticletolist (particleinfo,particleListaux)
    ENDIF
    GOTO 921
    922  close(2)

    ENDIF

  END SUBROUTINE add_particle_rest



  SUBROUTINE add_particle_rest_parallel (nparticles)

    CHARACTER(LEN = LEN_TRIM(case_dir) + 2*name_len) :: t_step_dir
    CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc
    LOGICAL :: dir_check
    INTEGER :: id, nparticles

#ifdef MFC_MPI

    REAL(KIND(0.D0)), DIMENSION(20)   :: inputvals
    REAL(KIND(0.D0)) :: id_real
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER(KIND=MPI_OFFSET_KIND) :: disp
    INTEGER :: view
 
    TYPE(particledata)    , POINTER :: particleinfo
    TYPE(particleListinfo), POINTER :: particleListaux

    INTEGER, DIMENSION(3)   :: cell
    LOGICAL                 :: indomain, particle_file, file_exist

    INTEGER, DIMENSION(2) :: gsizes, lsizes, start_idx_part
    INTEGER :: ifile, ireq, ierr, data_size, tot_data
    INTEGER :: i

    WRITE(file_loc, '(A,I0,A)') 'particle_mpi_io' , t_step_start, '.dat'
      file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
    INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)

    IF (file_exist) THEN
      IF(proc_rank .EQ. 0) THEN
        OPEN(9, FILE = TRIM(file_loc), FORM = 'unformatted', STATUS = 'unknown')
          READ(9) tot_data, time_real, dt_next_inp, tot_step
        CLOSE(9)
      END IF
    ELSE
      PRINT '(A)', TRIM(file_loc) // ' is missing. Exiting ...'
      CALL s_mpi_abort
    END IF
    
    CALL MPI_BCAST(tot_data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(time_real, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(dt_next_inp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    gsizes(1)=tot_data
    gsizes(2)=21
    lsizes(1)=tot_data
    lsizes(2)=21
    start_idx_part(1)=0
    start_idx_part(2)=0

    CALL MPI_TYPE_CREATE_SUBARRAY(2,gsizes,lsizes,start_idx_part,&
           MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,view,ierr)
    CALL MPI_TYPE_COMMIT(view,ierr)

    ! Open the file to write all flow variables
    WRITE(file_loc, '(A,I0,A)') 'particle' , t_step_start, '.dat'
    file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
    INQUIRE(FILE = TRIM(file_loc),EXIST = particle_file)
  
    IF (particle_file) THEN
  
        CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY, &
                           mpi_info_int,ifile,ierr)
      
        disp = 0d0
        CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,view, &
                                             'native',mpi_info_null,ierr)
    
        ALLOCATE(MPI_IO_DATA_particle(tot_data,1:21))
    
        CALL MPI_FILE_READ_ALL(ifile,MPI_IO_DATA_particle,21*tot_data, &
                                     MPI_DOUBLE_PRECISION,status,ierr)
    
        DO i = 1, tot_data
    
            id = INT(MPI_IO_DATA_particle(i,1))
            inputvals(1:20) = MPI_IO_DATA_particle(i,2:21)
        
            indomain = particle_in_domain(inputvals(1:3))
        
            IF (indomain.AND.(id.GT.0)) THEN
        
              ALLOCATE(particleinfo)
              nparticles = nparticles + 1
              particleinfo%id        = id
              particleinfo%x(1:3)    = inputvals(1:3) 
              particleinfo%xprev(1:3)= inputvals(4:6) 
              particleinfo%u(1:3)    = inputvals(7:9) 
              particleinfo%y(1:2)    = inputvals(10:11) 
              particleinfo%R0        = inputvals(12) 
              particleinfo%Rmax      = inputvals(13) 
              particleinfo%Rmin      = inputvals(14) 
              particleinfo%dphidt    = inputvals(15) 
              particleinfo%p         = inputvals(16) 
              particleinfo%mv        = inputvals(17) 
              particleinfo%mg        = inputvals(18) 
              particleinfo%betaT     = inputvals(19) 
              particleinfo%betaC     = inputvals(20) 
              particleinfo%equilibrium = .FALSE.
        
              cell    = -buff_size
        
              CALL locate_cell ( particleinfo%x,  cell, particleinfo%tmp%s )
        
              particleListaux => qbl%fp(cell(1), cell(2), cell(3))
              IF (particleListaux%nb.EQ.0) CALL addtocell_list ( cell )
              CALL transfertotmp (particleinfo) 
              CALL addparticletolist (particleinfo,particleListaux)
         
            ENDIF

        END DO
    
        DEALLOCATE(MPI_IO_DATA_particle)


    ENDIF
  
    CALL MPI_FILE_CLOSE(ifile,ierr)
#endif

  END SUBROUTINE add_particle_rest_parallel



  SUBROUTINE add_particle ( inputparticle, q_cons_vp, q_prim_vf )

    TYPE(particledata)    , POINTER :: particleinfo
    TYPE(particleListinfo), POINTER :: particleListaux
    TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_cons_vp
    TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
    REAL(KIND(0.D0)), DIMENSION(8),INTENT(IN)  :: inputparticle
    REAL(KIND(0.D0))                           :: pliq, volparticle, concvap, totalmass,&
                                                  kparticle,cpparticle,omegaN,PeG,PeT,rhol, cson, &
                                                  pcrit
    INTEGER, DIMENSION(3)                      :: cell
    REAL(KIND(0d0)), DIMENSION(2)              :: Re
    REAL(KIND(0d0)), DIMENSION( num_fluids, num_fluids ) :: We
 
    ALLOCATE(particleinfo)
    particleinfo%id        = id
    particleinfo%x(:)      = inputparticle(1:3)
    particleinfo%xprev(:)  = inputparticle(1:3)
    particleinfo%u(:)      = inputparticle(4:6)
    particleinfo%y(1)      = inputparticle(7)  ! particle radius
    particleinfo%R0        = inputparticle(7)
    particleinfo%y(2)      = inputparticle(8)  ! interface velocity
    particleinfo%Rmax      = 1.
    particleinfo%Rmin      = 1.
    particleinfo%dphidt    = 0.0d0

    IF(cyl_coord .AND. p.eq.0) THEN
        particleinfo%x(2) = DSQRT(particleinfo%x(2)**2d0+particleinfo%x(3)**2d0)
        !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
        particleinfo%x(3) = ATAN2(inputparticle(3),inputparticle(2))
        particleinfo%xprev = particleinfo%x
    ENDIF

    cell    = -buff_size
    CALL locate_cell ( particleinfo%x,  cell, particleinfo%tmp%s )

    pliq = Interpolate( particleinfo%tmp%s, q_prim_vf(E_idx))

    IF(pliq<0) PRINT *, "Negative pressure", proc_rank, q_cons_vp(E_idx)%sf(cell(1),cell(2),cell(3)),q_prim_vf(E_idx)%sf(cell(1),cell(2),cell(3)), cell

    CALL get_mixture_variables(q_prim_vf, q_cons_vp, pliq, cell(1), cell(2), cell(3), rhol, cson, Re, We)

    ! Intial particle pressure
    particleinfo%p =  pliq + 2.0*sigmabubble/particleinfo%R0
    IF(sigmabubble.NE.0.0d0) THEN
      pcrit = pvap - 4.0d0*sigmabubble/(3.d0*sqrt(3.0d0*particleinfo%p*particleinfo%R0**3/(2.0d0*sigmabubble)))
      ! fixme: For now I save only the last reference pressure (usually it is the same for
      ! all bubbles)
      pref  = particleinfo%p 
    ELSE
      pcrit = 0.0d0
    ENDIF

    particleinfo%equilibrium = .FALSE.

    ! Initial particle mass
    volparticle = 4.0d0/3.0d0*pi*particleinfo%R0**3 ! volume
    particleinfo%mv = pvap*volparticle*MWvap/Runiv*DBLE(massflag) ! vapermass
    particleinfo%mg = (particleinfo%p-pvap*DBLE(massflag))*volparticle*MWgas/Runiv ! gasmass
    IF (particleinfo%mg.LE.0.0d0) STOP 'the initial mass of gas inside the bubble is negative. Check your initial conditions'
    totalmass = particleinfo%mg + particleinfo%mv ! totalmass

    ! particle natural frequency
    concvap = particleinfo%mv/(particleinfo%mv+particleinfo%mg)
    omegaN = (3.0d0*(particleinfo%p-pvap*REAL(massflag))+4.0d0*sigmabubble/particleinfo%R0)/rhol
    IF (pvap*REAL(massflag).GT.particleinfo%p) THEN
      print *, 'Not allowed: bubble initially located in a region with pressure below the vapor pressure'
      print *, 'location:', particleinfo%x(1:3)
      stop
    ENDIF
    omegaN = DSQRT(omegaN/particleinfo%R0**2)

    cpparticle = concvap*cpvapor    + (1.0d0 - concvap)*cpgas
    kparticle  = concvap*kvapor     + (1.0d0 - concvap)*kgas
    ! Mass and heat transfer coefficients (based on Preston 2007)

    PeT   = totalmass/volparticle * cpparticle * particleinfo%R0**2*omegaN/kparticle
    particleinfo%betaT = transfercoeff(PeT,1.0d0)*REAL(heatflag)
    PeG   = particleinfo%R0**2*omegaN/diffcoefvap
    particleinfo%betaC = transfercoeff(PeG,1.0d0)*REAL(massflag)

    ! terms to work out directly the heat flux in getfluxes
    particleinfo%betaT = particleinfo%betaT*kparticle

    IF (particleinfo%mg.LE.0.0d0) STOP 'PROBLEM WITH THE MASS OF THE particle, CHECK PARTICLES ARE INSIDE THE DOMAIN'

    particleListaux => qbl%fp(cell(1),cell(2),cell(3))
    IF (particleListaux%nb.EQ.0) CALL addtocell_list ( cell )
    CALL transfertotmp (particleinfo) 
    CALL addparticletolist (particleinfo,particleListaux)

  END SUBROUTINE add_particle

  !------------------------------------------------------------------


  SUBROUTINE s_initialize_particles(q_cons_vp, q_prim_vf)

    USE m_particles_output

    ! Conservative variables
    TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vp
    TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_prim_vf
    TYPE(particlenode), POINTER           :: particle
    REAL(KIND(0.D0)), DIMENSION(8)   :: inputparticle
    REAL(KIND(0.D0))                 :: qtime
    INTEGER  :: i,j,k,l,nparticles
    LOGICAL  :: file_exist,indomain

    CALL s_populate_variables_buffers(q_cons_vp, q_particle)

    IF(model_eqns == 2 .AND. (adv_alphan .NEQV. .TRUE.)) THEN
             
        q_cons_vp(sys_size)%sf = 1d0
              
        DO i = adv_idx%beg, adv_idx%end
             q_cons_vp(sys_size)%sf = &
             q_cons_vp(sys_size)%sf - &
             q_cons_vp(i)%sf
        END DO
              
    END IF

    ix%beg = -buff_size; ix%end = m + buff_size
    iy%beg = -buff_size; iy%end = n + buff_size
    IF(p > 0) THEN
        iz%beg = -buff_size; iz%end = p + buff_size
    END IF

    ! To get the pressure to be used in add_particle()
    !CALL s_convert_conservative_to_primitive_variables(q_cons_vp, q_prim_vf, gm_alpha_qp(0,0,0)%vf, ix, iy, iz, q_particle(1))
    CALL s_convert_conservative_to_primitive_variables(q_cons_vp, q_prim_vf, gm_alpha_qp%vf, ix, iy, iz, q_particle(1))

    CALL compute_cell_centers ()

    nparticles = 0

    IF(num_procs > 1) THEN
      CALL initialize_particles_mpi()
    ENDIF

    CALL initialize_particlelists

    IF (t_step_start.eq.0) THEN
      INQUIRE (file='input/particles.dat'   ,exist=file_exist    )
      IF (file_exist) THEN
        OPEN(unit=85,file='input/particles.dat',form='formatted')
        101 READ(85 ,*,end=102)  (inputparticle(i), i=1,8)
        indomain = particle_in_domain(inputparticle(1:3))
        CALL get_part_id
        IF (indomain) THEN
          nparticles = nparticles + 1
          CALL add_particle(inputparticle, q_cons_vp, q_prim_vf)
        ENDIF
        GOTO 101
        102 CONTINUE
      ELSE
        STOP "If you include particles, you have to initializate them in input/particles.dat"
      END IF
    ELSE
      IF(parallel_io .NEQV. .TRUE.) THEN
        CALL add_particle_rest (nparticles)
      ELSE
        CALL add_particle_rest_parallel (nparticles)
      END IF
    ENDIF

    print *, " PARTICLES RUNNING, in proc", proc_rank, "number:", nparticles, "/",id
 
    CALL Create_sublists ()

    ! Apply density correction
    IF (avgdensFlag) THEN
      particle  => particlesubList%List%next
      DO WHILE(Associated(particle))
        CALL transfertotmp (particle%data)
        particle => particle%next
      ENDDO
      CALL voidfraction(q_cons_vp)
      IF ( solverapproach.EQ.1 ) THEN
       !Definition averaged quantities
        q_cons_vp(E_idx)%sf(0:m,0:n,0:p) = q_cons_vp(E_idx)%sf(0:m,0:n,0:p) * q_particle(1)%sf(0:m,0:n,0:p)
        DO i = 1, cont_idx%end 
          q_cons_vp(i)%sf(0:m,0:n,0:p)   = q_cons_vp(i)%sf(0:m,0:n,0:p)     * q_particle(1)%sf(0:m,0:n,0:p)
        END DO
        DO i = mom_idx%beg, mom_idx%end
          q_cons_vp(i)%sf(0:m,0:n,0:p)   = q_cons_vp(i)%sf(0:m,0:n,0:p)     * q_particle(1)%sf(0:m,0:n,0:p)
        ENDDO
      ENDIF
    ENDIF

    qtime = 0.0d0

    IF (particleoutFlag)  CALL write_particles (qtime)


    CALL s_populate_variables_buffers(q_cons_vp, q_particle)


    IF(model_eqns == 2 .AND. (adv_alphan .NEQV. .TRUE.)) THEN
             
        q_cons_vp(sys_size)%sf = 1d0
              
        DO i = adv_idx%beg, adv_idx%end
             q_cons_vp(sys_size)%sf = &
             q_cons_vp(sys_size)%sf - &
             q_cons_vp(i)%sf
        END DO
              
    END IF

    ix%beg = -buff_size; ix%end = m + buff_size
    iy%beg = -buff_size; iy%end = n + buff_size
    IF(p > 0) THEN
        iz%beg = -buff_size; iz%end = p + buff_size
    END IF

    !CALL s_convert_conservative_to_primitive_variables(q_cons_vp, q_prim_vf, gm_alpha_qp(0,0,0)%vf, &
    !                                                                      ix, iy, iz, q_particle(1))
    CALL s_convert_conservative_to_primitive_variables(q_cons_vp, q_prim_vf, gm_alpha_qp%vf, &
                                                                   ix, iy, iz, q_particle(1))
    
  END SUBROUTINE s_initialize_particles


  ! -------------------------------------------------------------------------------------------   
  FUNCTION transfercoeff (Pe,omegaN) 

    REAL(KIND(0.D0)) :: transfercoeff , Pe, omegaN
    COMPLEX          :: transferfunc, auxc

    transferfunc = CSQRT(CMPLX(0.0d0,Pe*omegaN))
    auxc = (cexp(-CMPLX(2.0d0,0.0d0)*transferfunc)+1.0d0)/(-cexp(-CMPLX(2.0d0,0.0d0)*transferfunc)+1.0d0)
    transferfunc = transferfunc*auxc-1.0d0
    transferfunc = 1.0d0/transferfunc - 3.0d0/CMPLX(0.0d0,Pe*omegaN)
    transfercoeff = REAL(1.0d0/transferfunc)

  END FUNCTION transfercoeff


  SUBROUTINE get_stddsv(cell, kernel, volpart, stddsv)

    REAL(KIND(0.D0))         :: chardist, volpart, rad, stddsv, charvol
    INTEGER, DIMENSION(3)    :: cell
    INTEGER                  :: kernel

    CALL get_char_dist(cell,chardist)
    CALL get_char_vol(cell,charvol)
    kernel = smoothtype
  
    IF (((volpart/charvol).GT.0.5d0*valmaxvoid).OR.(smoothtype.EQ.1)) THEN
      kernel = 1
      rad    = (3.0d0*volpart/(4.0d0*pi))**(1.0d0/3.0d0)
      stddsv = 1.0d0*epsilonb*MAX(chardist,rad)
    ELSE
      stddsv = 0.0d0
    ENDIF

  END SUBROUTINE get_stddsv


  SUBROUTINE voidfraction (q)

    USE m_kernel_functions

    TYPE(particlenode),POINTER       :: particle
    TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q
    REAL(KIND(0.D0))               :: volpart, totmass, stddsv, volpart2
    REAL(KIND(0.D0)), DIMENSION(3) :: nodecoord
    INTEGER, DIMENSION(3)          :: cell
    INTEGER                        :: i, j, k, l, kernel
    INTEGER, DIMENSION(3,2)        :: rangecells

    q_particle(1)%sf = 0.0d0
    q_particle(2)%sf = 0.0d0
    nodecoord(3)     = 0

    rangecells(1:3,2) = -buff_size
    rangecells(1,1) = m
    rangecells(2,1) = n
    rangecells(3,1) = p

    IF (projectiontype.EQ.0) THEN   

      particle  => particlesubList%List%next
      DO WHILE(Associated(particle))
        volpart = 4.0d0/3.0d0*pi*particle%data%tmp%y(1)**3
        cell = get_cell_from_s(particle%data%tmp%s)
        CALL update_rangecells (cell, rangecells)
        nodecoord(1) = particle%data%tmp%x(1)
        nodecoord(2) = particle%data%tmp%x(2)
        IF (p > 0) nodecoord(3) = particle%data%tmp%x(3)
        CALL get_stddsv(cell, kernel, volpart, stddsv)
        CALL smoothfunction ( q_particle(1), nodecoord, cell , volpart, kernel, stddsv)
        particle => particle%next
      ENDDO
	  
      subrange(1)%beg=max(rangecells(1,1), -buff_size)
      subrange(1)%end=min(rangecells(1,2),m+buff_size)
      subrange(2)%beg=max(rangecells(2,1), -buff_size)
      subrange(2)%end=min(rangecells(2,2),n+buff_size)
      subrange(3)%beg=max(rangecells(3,1), -buff_size)
      subrange(3)%end=min(rangecells(3,2),p+buff_size)

    ELSE
    
      particle  => particlesubList%List%next
      DO WHILE(Associated(particle))
        volpart = 4.0d0/3.0d0*pi*particle%data%tmp%y(1)**3
        CALL remeshdelta ( q_particle(2), particle%data, volpart, rangecells )     
        particle => particle%next
      ENDDO

      IF (bubblesources) q_particle(4)%sf = q_particle(2)%sf

      subrange(1)%beg=max(rangecells(1,1), -buff_size)
      subrange(1)%end=min(rangecells(1,2),m+buff_size)
      subrange(2)%beg=max(rangecells(2,1), -buff_size)
      subrange(2)%end=min(rangecells(2,2),n+buff_size)
      subrange(3)%beg=max(rangecells(3,1), -buff_size)
      subrange(3)%end=min(rangecells(3,2),p+buff_size)
  
      DO k=0,p
        DO j=0,n
          DO i=0,m
            cell(1) = i 
            cell(2) = j
            cell(3) = k
            nodecoord(1)=x_cc_lp(cell(1))
            nodecoord(2)=y_cc_lp(cell(2))
            IF (p.GT.0) nodecoord(3)=z_cc_lp(cell(3))

            IF (q_particle(2)%sf(i,j,k).NE.0.0d0) THEN
              CALL get_stddsv(cell, kernel, q_particle(2)%sf(i,j,k), stddsv)
              CALL smoothfunction ( q_particle(1), nodecoord, cell , q_particle(2)%sf(i,j,k), kernel, stddsv)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    !I store 1-beta (I should probably do it initializing
    !                q_particle(1)%fp=1 and using a negative volume)
    q_particle(1)%sf = 1. - q_particle(1)%sf

    subrange(1)%beg=max(subrange(1)%beg-CEILING(5*epsilonb), -buff_size*1.0d0)
    subrange(1)%end=min(subrange(1)%end+CEILING(5*epsilonb),(m+buff_size)*1.0d0)
    subrange(2)%beg=max(subrange(2)%beg-CEILING(5*epsilonb), -buff_size*1.0d0)
    subrange(2)%end=min(subrange(2)%end+CEILING(5*epsilonb),(n+buff_size)*1.0d0)

    IF (p > 3) THEN
      subrange(3)%beg=max(subrange(3)%beg-CEILING(5*epsilonb), -buff_size*1.0d0)
      subrange(3)%end=min(subrange(3)%end+CEILING(5*epsilonb),(p+buff_size)*1.0d0)
    ENDIF

    !=== dbetadt==========
    IF (bubblesources) THEN
      q_particle(2)%sf = 0.0d0
      q_particle(3)%sf = 0.0d0
      IF((clusterflag.GE.4)) THEN
         q_particle(5)%sf = 0.0d0
      END IF

      IF (projectiontype.EQ.0) THEN

        particle  => particlesubList%List%next
        DO WHILE(Associated(particle))
          volpart = 4.0d0/3.0d0*pi*particle%data%tmp%y(1)**3
          volpart2 = volpart
          cell = get_cell_from_s(particle%data%tmp%s)
          CALL get_stddsv(cell, kernel, volpart, stddsv)
          volpart = 4.0d0*pi*particle%data%tmp%y(1)**2*particle%data%tmp%y(2)
          nodecoord(1) = particle%data%tmp%x(1)
          nodecoord(2) = particle%data%tmp%x(2)
          IF (p > 0) nodecoord(3) = particle%data%tmp%x(3)
          IF(clusterflag.GE.4) THEN
            CALL smoothfunction ( q_particle(5), nodecoord, cell , volpart, kernel, stddsv, volpart2)
          END IF
          CALL smoothfunction ( q_particle(2), nodecoord, cell , volpart, kernel, stddsv)
          particle => particle%next
        ENDDO

      ELSE

        nodecoord(3)=0
        particle  => particlesubList%List%next
        DO WHILE(Associated(particle))
          volpart = 4.0d0*pi*particle%data%tmp%y(1)**2*particle%data%tmp%y(2)
          CALL remeshdelta ( q_particle(3), particle%data, volpart )     
          particle => particle%next
        ENDDO


        DO i=0,m
          DO j=0,n
            DO k=0,p
              cell(1) = i 
              cell(2) = j
              cell(3) = k
              nodecoord(1)=x_cc_lp(cell(1))
              nodecoord(2)=y_cc_lp(cell(2))
              IF (p > 0) nodecoord(3)=z_cc_lp(cell(3))

              IF (q_particle(3)%sf(i,j,k).NE.0.0d0) THEN
                CALL get_stddsv(cell, kernel, q_particle(4)%sf(i,j,k), stddsv)
                CALL smoothfunction( q_particle(2), nodecoord, cell , q_particle(3)%sf(i,j,k), kernel, stddsv)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF  
    ENDIF

    subrange(1)%beg=max(subrange(1)%beg,0*1.0d0)
    subrange(1)%end=min(subrange(1)%end,m*1.0d0)
    subrange(2)%beg=max(subrange(2)%beg,0*1.0d0)
    subrange(2)%end=min(subrange(2)%end,n*1.0d0)

    IF (p > 0) THEN
      subrange(3)%beg=max(subrange(3)%beg,0*1.0d0)
      subrange(3)%end=min(subrange(3)%end,p*1.0d0)
    ENDIF

    ! Limiting void fraction given max value
    DO k=0,p
      DO j=0,n
        DO i=0,m
          q_particle(1)%sf(i,j,k) = max(q_particle(1)%sf(i,j,k),1.d0-valmaxvoid)
        ENDDO
      ENDDO
    ENDDO
	

  END SUBROUTINE voidfraction


  SUBROUTINE add_sources(q,dq,q_prim)

    TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q
    TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: dq
    TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim
    INTEGER :: i,j,k,l

      DO k=0,p
        DO j=0,n
          DO i=0,m
            IF (q_particle(1)%sf(i,j,k).GT.(1.0d0-valmaxvoid)) THEN
              DO l=1,E_idx
                IF(clusterflag.GE.4) THEN
                   dq(l)%sf(i,j,k) = dq(l)%sf(i,j,k) + q(l)%sf(i,j,k)*(q_particle(2)%sf(i,j,k)+q_particle(5)%sf(i,j,k))
                ELSE
                   dq(l)%sf(i,j,k) = dq(l)%sf(i,j,k) + q(l)%sf(i,j,k)/q_particle(1)%sf(i,j,k)*q_particle(2)%sf(i,j,k)
                END IF
              ENDDO
            ENDIF
          ENDDO 
        ENDDO 
      ENDDO

      DO l=1,DIM
      CALL gradient_dir(q_prim(E_idx), q_particle(3),l)

      DO k=0,p
        DO j=0,n
          DO i=0,m
            IF (q_particle(1)%sf(i,j,k).GT.(1.0d0-valmaxvoid)) THEN
              dq(mom_idx%beg+l-1)%sf(i,j,k) = dq(mom_idx%beg+l-1)%sf(i,j,k) - (1.0d0-q_particle(1)%sf(i,j,k))/q_particle(1)%sf(i,j,k)*q_particle(3)%sf(i,j,k)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      !source in energy
      q_particle(3)%sf = q_prim(E_idx)%sf * q_prim(mom_idx%beg+l-1)%sf
      CALL gradient_dir(q_particle(3), q_particle(4),l)

      DO k=0,p
        DO j=0,n
          DO i=0,m
            IF (q_particle(1)%sf(i,j,k).GT.(1.0d0-valmaxvoid)) THEN
              dq(E_idx)%sf(i,j,k) = dq(E_idx)%sf(i,j,k) - q_particle(4)%sf(i,j,k)*(1.0d0-q_particle(1)%sf(i,j,k))/q_particle(1)%sf(i,j,k)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE add_sources


  !!----------------------------------------------------------------------------------------  
  !!----------------------------------------------------------------------------------------  
  !!                              BUBBLE DYNAMICS SUBROUTINES
  !!----------------------------------------------------------------------------------------  
  !!----------------------------------------------------------------------------------------  

  SUBROUTINE RKparticledyn (qtime,step,q,t_step,q_prim,dq,largestep)

  TYPE(particlenode),POINTER                  :: particle
  TYPE(scalar_field), DIMENSION(sys_size)   :: q
  TYPE(scalar_field), DIMENSION(sys_size)   :: q_prim
  TYPE(scalar_field), DIMENSION(sys_size), OPTIONAL :: dq
  REAL(KIND(0.D0)),DIMENSION(5)             :: intvalues
  REAL(KIND(0.D0)),DIMENSION(3)             :: totalforce, DupDt
  REAL(KIND(0.D0))                          :: gammaparticle,vaporflux,heatflux,qtime
  INTEGER,DIMENSION(3)                      :: cell
  INTEGER                                   :: i,step,j,k
  INTEGER,INTENT(IN)                        :: t_step
  LOGICAL,OPTIONAL                          :: largestep


  IF (avgdensFlag) THEN
    DO i=1, sys_size
      dq(i)%sf=0.0d0
    ENDDO
  ENDIF
  
  !update vbles
  IF (avgdensFlag) CALL voidfraction (q)
  DO i = 1, cont_idx%end
     q_prim(i)%sf => q(i)%sf
  END DO
  DO i = adv_idx%beg, sys_size
     q_prim(i)%sf => q(i)%sf
  END DO
  
  
  IF (coupledFlag) THEN
    !CALL s_compute_rhs(q, q_prim, dq, t_step=t_step, qtime=qtime, largestep=largestep)
    !if (proc_rank==0) print*,size(q), associate(q_prim), associate(dq), associate(t_step), associate(qtime)
    CALL s_compute_rhs(q, q_prim, dq, t_step=t_step, qtime=qtime)
    IF (num_procs > 1) THEN
         CALL bcst_largestep(largestep)
    ENDIF
    !IF (largestep) RETURN
  ENDIF
  

  IF (avgdensFlag) THEN
    
     CALL s_populate_variables_buffers(q,q_particle)
    
     IF(model_eqns == 2 .AND. (adv_alphan .NEQV. .TRUE.)) THEN
             
        q(sys_size)%sf = 1d0
              
        DO i = adv_idx%beg, adv_idx%end
             q(sys_size)%sf = &
             q(sys_size)%sf - &
             q(i)%sf
        END DO
              
     END IF
   
     ix%beg = -buff_size; ix%end = m + buff_size
     iy%beg = -buff_size; iy%end = n + buff_size
     IF(p > 0) THEN
       iz%beg = -buff_size; iz%end = p + buff_size
     END IF
     !CALL s_convert_conservative_to_primitive_variables(q, q_prim, gm_alpha_qp(0,0,0)%vf, ix, iy, iz, q_particle(1))
     CALL s_convert_conservative_to_primitive_variables(q, q_prim, gm_alpha_qp%vf, ix, iy, iz, q_particle(1))

     ! For conservative to primitive below
     DO i = 1, sys_size
       !q_cons_qp(0,0,0)%vf(i)%sf => q(i)%sf
       q_cons_qp%vf(i)%sf => q(i)%sf
     END DO

     IF ((clusterflag.GT.0).OR.(correctpresFlag)) CALL potentials (q_prim,q,qtime)

  ENDIF

  IF (bubblesources) CALL add_sources(q, dq, q_prim)

  particle  => particlesubList%List%next
  DO WHILE(Associated(particle))
    IF (.NOT.particle%data%equilibrium) THEN
      CALL getfluxes(particle%data, vaporflux, heatflux, gammaparticle)
      particle%data%dbdt(step)%dpbdt = deriv_gaspressure(particle%data, vaporflux, heatflux, gammaparticle )
      particle%data%dbdt(step)%dmvdt = 4.0d0*pi*particle%data%tmp%y(1)**2*vaporflux 
    ELSE
      particle%data%dbdt(step)%dpbdt = 0.0d0
      particle%data%dbdt(step)%dmvdt = 0.0d0
    ENDIF

    particle%data%dbdt(step)%dxdt(:) = 0.0d0
    particle%data%dbdt(step)%dudt(:) = 0.0d0

    particle => particle%next
  ENDDO

  !Radial motion
  IF (RPflag) THEN
    particle  => particlesubList%List%next
    DO WHILE(Associated(particle))
      IF (.NOT.particle%data%equilibrium) THEN
        CALL RPeq ( particle%data, step, q, q_prim, qtime)
      ELSE
        particle%data%dbdt(step)%dydt(2) = 0.
      ENDIF
      particle%data%dbdt(step)%dydt(1) = particle%data%tmp%y(2) !the derivative of the radius is the velocity
      particle => particle%next
    ENDDO
  ELSE
    particle  => particlesubList%List%next
    DO WHILE(Associated(particle))
      particle%data%dbdt(step)%dydt(2) = 0.0d0
      particle%data%dbdt(step)%dydt(1) = 0.0d0
      particle => particle%next
    ENDDO
  END IF
 
  DO i = 1, sys_size
    !NULLIFY(q_cons_qp(0,0,0)%vf(i)%sf)
    NULLIFY(q_cons_qp%vf(i)%sf)
  END DO

  END SUBROUTINE RKparticledyn


  FUNCTION pressureliq_int( pbubble, radius, bubblevel )

    REAL(KIND(0.D0)) :: pressureliq_int, radius,bubblevel, pbubble

    pressureliq_int=  pbubble -  2.0d0*sigmabubble/radius - 4.*viscref*bubblevel/radius

  END FUNCTION pressureliq_int


  SUBROUTINE RPeq (bubbletmp, step, q, q_prim, qtime)
  !subroutine containing the Rayleight--Plesset equation

    !TYPE( field_position ), DIMENSION(system_size) :: q
    TYPE(scalar_field), DIMENSION(sys_size)  :: q
    TYPE(scalar_field), DIMENSION(sys_size)  :: q_prim
    TYPE(particledata)                 :: bubbletmp
    TYPE(particlederivative)           :: dbdt
    REAL(KIND(0.D0))                   :: pliqint, pbubble, deltaP, pinf, termI, aux1, aux2, velint, &
                                          rhol, cson
    INTEGER, DIMENSION(3)              :: cell
    INTEGER                            :: step
    REAL(KIND(0d0)), DIMENSION(2)      :: Re
    REAL(KIND(0d0)), DIMENSION( num_fluids, num_fluids ) :: We
    REAL(KIND(0d0)), OPTIONAL          :: qtime

    REAL(KIND(0d0)) :: temp
         
    pbubble = min(bubbletmp%tmp%p, 1.d0) ! pres in the bubble
    pliqint = pressureliq_int(pbubble, bubbletmp%tmp%y(1), bubbletmp%tmp%y(2)) ! pres outside the bubble on the surface
    
    cell = get_cell_from_s(bubbletmp%tmp%s)

    ! getting p_inf
    pinf = get_pinf1(bubbletmp%tmp%s, q_prim(e_idx), 1, aux1, aux2)
     
    call get_mixture_variables(q_prim, q, pinf, cell(1), cell(2), cell(3), rhol, cson, re, we, q_particle(1))

    deltaP = pliqint - pinf
    termI = 0.0d0

    velint = bubbletmp%tmp%y(2) - bubbletmp%dbdt(step)%dmvdt/(4.0d0*pi*bubbletmp%tmp%y(1)**2*rhol)

    bubbletmp%dbdt(step)%dydt(2) =  ((1.0d0+velint/cson)*deltaP/rhol + termI &
                                    + bubbletmp%dbdt(step)%dpbdt*bubbletmp%tmp%y(1)/rhol/cson  &
                                    - velint**2*3.0d0/2.0d0*(1.0d0-velint/3.0d0/cson))         &
                                    / (bubbletmp%tmp%y(1)*(1.0d0-velint/cson))

  END SUBROUTINE RPeq


  FUNCTION get_pinf1(scoord, pres, ptype, preterm1, term2, Romega)

    TYPE(particlenode),POINTER         :: bubble
    TYPE(scalar_field)                 :: pres
    REAL(KIND(0.D0)), DIMENSION(3)     :: distance,center
    REAL(KIND(0.D0))                   :: get_pinf1, jac, dij, dpotjdt, dc, vol, aux,&
                                          volgas, term1, Rbeq, denom, Rmax, stddsv,&
                                          charvol, charpres, charvol2, charpres2
    REAL(KIND(0.D0)), DIMENSION(3)     :: scoord
    REAL(KIND(0.D0)), OPTIONAL         :: preterm1, term2, Romega
    INTEGER, DIMENSION(3)              :: cell, cellaux
    INTEGER, DIMENSION(3)              :: epsilonbaux
    INTEGER                            :: ptype !1=p at infinity, 2= averaged P at the bubble location
    INTEGER                            :: i, j, k, dir
    LOGICAL                            :: celloutside

    get_pinf1 = 0.0d0
    cell(:) = get_cell_from_s (scoord)
    
    IF((clusterflag.eq.0)) THEN
    !getting p_cell in terms of only the current cell by interpolation

       CALL get_char_vol (cell, vol)
   
       bubble  => qbl%fp(cell(1),cell(2),cell(3))%List%next
   
       Rmax=0.0d0
   
       DO WHILE(Associated(bubble))
         Rmax  = Rmax+bubble%data%tmp%y(1)**3
         bubble  => bubble%next
       ENDDO
   
       ! Surrogate bubble radius?
       Rmax = Rmax**(1.0d0/3.0d0)
 
       ! Getting the cell volulme as Omega
       CALL get_char_dist(cell,stddsv)
     
       !p_cell (interpolated)
       get_pinf1 = Interpolate( scoord, pres )
 
       !R_Omega
       dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)
   
    ELSE IF(clusterflag .eq. 1) THEN
    ! Just making the characteristic volume 3^3 times bigger than the cell volume
    ! DO NOT USE FOR PRODUCTION

       ! Getting Omega (not Omega_L)
       CALL get_char_vol (cell, vol)
       vol=vol*2.7d1 !Just multiplying by 3^3=27
   
       bubble  => qbl%fp(cell(1),cell(2),cell(3))%List%next
   
       Rmax=0.0d0
   
       DO WHILE(Associated(bubble))
         Rmax  = Rmax+bubble%data%tmp%y(1)**3
         bubble  => bubble%next
       ENDDO
   
       ! Surrogate bubble radius?
       Rmax = Rmax**(1.0d0/3.0d0)
   
       CALL get_char_dist(cell,stddsv)
     
       !p_cell to get interpolate
       get_pinf1 = Interpolate( scoord, pres )

       !R_Omega
       dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)
 
    ELSE IF(clusterflag .eq. 2) THEN
    !getting p_cell in terms of Omega, over which the kernel is applied

       PRINT '(A)', 'Bad cluster model. Exiting ...'
       CALL s_mpi_abort()

       ! Range of cells included in Omega
       IF(smoothtype.EQ.1) THEN
           epsilonbaux(:) = MIN(CEILING(epsilonb*3d0),buff_size)
       ELSE IF(projectiontype.EQ.1) THEN
           epsilonbaux(:) = 3
       END IF
       
       charvol = 0d0
       charpres = 0d0

 
       DO i=(cell(1)-epsilonbaux(1)),(cell(1)+epsilonbaux(1))
          DO j=(cell(2)-epsilonbaux(2)),(cell(2)+epsilonbaux(2))
             DO k=(cell(3)-epsilonbaux(3)),(cell(3)+epsilonbaux(3))
                !Excluding cellaux(:) outside the domain
                cellaux(1)=MIN(MAX(i,-buff_size),buff_size)
                cellaux(2)=MIN(MAX(j,-buff_size),buff_size)
                cellaux(3)=MIN(MAX(k,-buff_size),buff_size)
                CALL get_char_vol(cellaux,vol)
                charvol  = charvol + vol!*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
                charpres = charpres + pres%sf(cellaux(1),cellaux(2),cellaux(3)) &
                           * vol!*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
                charvol2  = charvol2 + vol*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
                charpres2 = charpres2 + pres%sf(cellaux(1),cellaux(2),cellaux(3)) &
                           * vol*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
             ENDDO
          ENDDO 
       ENDDO

       bubble  => qbl%fp(cell(1),cell(2),cell(3))%List%next
 
       CALL get_char_dist(cell,stddsv)
   
       get_pinf1 = charpres2/charvol2

       !R_Omega
       vol=charvol
       dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)
 
    ELSE IF(clusterflag .GE. 3 ) THEN

       ! Range of cells included in Omega
       IF(smoothtype.EQ.1) THEN
           epsilonbaux(:) = 3
       ELSE IF(projectiontype.EQ.1) THEN
           epsilonbaux(:) = 3
       END IF
 
       charvol = 0.d0
       charpres = 0.d0
       charvol2 = 0.d0
       charpres2 = 0.d0
       vol = 0.d0

       IF (DIM.EQ.3) THEN
          k = -epsilonbaux(3)
       ELSE
          k = 0
       ENDIF
      
       i=-epsilonbaux(1);j=-epsilonbaux(2)
  
3001   IF ((i.LE.epsilonbaux(1)).AND.(j.LE.epsilonbaux(2))) THEN
  
           celloutside = .FALSE.
  
           cellaux(1) = cell(1) + i
           cellaux(2) = cell(2) + j
           cellaux(3) = cell(3) + k
  
              IF (cellaux(1).LT.-buff_size) THEN
                  celloutside = .TRUE.
                  i = i+1
                  !i = i+abs(cellaux(1)-(-buff_size))
              ENDIF
  
              !check ghost part in y direction
              IF (cellaux(2).LT.-buff_size) THEN
                  celloutside = .TRUE.
                  j = j+1
                  !j = j+abs(cellaux(2)-(-buff_size))
              ENDIF
  
              IF(cyl_coord.AND.DIM.NE.3) THEN
                  IF (y_cc_lp(cellaux(2)).LT.0d0) THEN
                      celloutside = .TRUE.
                      j = j+1
                  ENDIF
              END IF
  
              ! Temp
              IF(DIM.eq.3) THEN
                  IF (cellaux(3).LT.-buff_size) THEN
                      celloutside = .TRUE.
                      k = k+1
                      !k = k+abs(cellaux(3)-(-buff_size))
                  ENDIF
              END IF

              if (cellaux(1).GT.m+buff_size) celloutside =.TRUE.
              if (cellaux(2).GT.n+buff_size) celloutside =.TRUE.
              if (cellaux(3).GT.p+buff_size) celloutside =.TRUE. 
              
              !print*, 'cellaux values x, y, z:', cellaux(1), cellaux(2), cellaux(3), celloutside, epsilonbaux(1), buff_size
              if (.not. celloutside) then
                 CALL get_char_vol(cellaux,vol)
                 charvol  = charvol + vol!*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
                 charpres = charpres + pres%sf(cellaux(1),cellaux(2),cellaux(3)) &
                         * vol!*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
                 charvol2  = charvol2 + vol*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
                 charpres2 = charpres2 + pres%sf(cellaux(1),cellaux(2),cellaux(3)) &
                         * vol*q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3))
              end if

           IF (j.LT.epsilonbaux(2)) THEN
               j = j+1
               GOTO 3001
           ENDIF
  
3002       j=-epsilonbaux(2)
           i = i+1
           GOTO 3001
  
       ENDIF
 
3003   IF ((DIM.EQ.3).AND.(k.LT.epsilonbaux(3))) THEN
           k = k+1
           i=-epsilonbaux(1);j=-epsilonbaux(2)
           GOTO 3001
       ENDIF

       get_pinf1 = charpres2/charvol2

       vol=charvol
       !dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)
       dc = (3.0d0*abs(vol)/(4.0d0*pi))**(1.0d0/3.0d0)  !positive volume
 
    ELSE

       PRINT '(A)', 'Check cluterflag. Exiting ...'
       CALL s_mpi_abort()

    END IF

    IF (correctpresFlag.AND.PRESENT(preterm1)) THEN

      dpotjdt = 0.0d0 !potential derivative contribution from other bubbles?
      volgas  = 0.0d0
      term1 = 0.0d0
      term2 = 0.0d0
      denom = 0.0d0
      bubble  => qbl%fp(cell(1),cell(2),cell(3))%List%next

      DO WHILE(Associated(bubble))
        volgas  = volgas + bubble%data%tmp%y(1)**3 !surrogate bubble volume
        denom   = denom + bubble%data%tmp%y(1)**2
        term1   = term1 + bubble%Data%dphidt*bubble%data%tmp%y(1)**2
        term2   = term2 + bubble%data%tmp%y(2)*bubble%data%tmp%y(1)**2
        bubble  => bubble%next
      ENDDO

      Rbeq = volgas**(1.0d0/3.0d0) !surrogate bubble radius
      aux = dc**3 - Rbeq**3
      term2 = term2/denom
      term2 = 3.0d0/2.0d0*term2**2*Rbeq**3*(1.0d0-Rbeq/dc)/aux
      preterm1  = 3.0d0/2.0d0*Rbeq*(dc**2 - Rbeq**2)/(aux*denom)

      !Control volume radius
      IF(PRESENT(Romega)) Romega = dc

      ! Getting p_inf
      IF (ptype.EQ.1) THEN
        get_pinf1 =  get_pinf1 + preterm1*term1 + term2
      ENDIF
 
    ENDIF

    ! Test if pinf is valid
    if ((get_pinf1>10) .OR. (get_pinf1<-10) .OR. (ieee_is_nan(get_pinf1))) THEN
        print*, 'Extreme pinf value of', get_pinf1,'charvol2', charvol2,&
                   'charpres2',charpres2,'dc',dc,'Rbeq',Rbeq,'vol', vol,&
                                 'charvol',charvol,'in rank', proc_rank,&
                                   'cell m n p', cell(1), cell(2), cell(3)
        if (charvol2==0) print*, 'charvol2 is equal to zero, q_particle%sf is ', &
            q_particle(1)%sf(cellaux(1),cellaux(2),cellaux(3)), cellaux(1), cellaux(2), cellaux(3)
        if (charpres2==0) print*, 'charpres2 is equal to zero, pres%sf is ', pres%sf(cellaux(1),cellaux(2),cellaux(3))
        call s_mpi_abort()
    endif


  END FUNCTION get_pinf1


  SUBROUTINE getfluxes(bubbletmp, vaporflux, heatflux, gammabubble)

      REAL(KIND(0.D0)) :: vaporflux, heatflux, concvapint,bubbleTemp,volbubble,kbubble,&
                          avgconc, Rmixt,gammabubble,rhogas
      TYPE(particledata) :: bubbletmp

      IF (massflag.EQ.0) THEN
        concvapint = 0.d0
      ELSE
        concvapint   = MWgas/MWvap*(bubbletmp%tmp%p/pvap-1.0d0)
        concvapint   = 1.0d0/(1.0d0+concvapint)
      ENDIF

      bubbleTemp   = (bubbletmp%mg/MWgas + bubbletmp%tmp%mv/MWvap)*Runiv
      volbubble    = 4.0d0/3.0d0*pi*bubbletmp%tmp%y(1)**3
      bubbleTemp   = bubbletmp%tmp%p*volbubble/bubbleTemp

      gammabubble  = concvapint*gammavapor + (1.0d0 - concvapint)*gammagas !needed later (deriv_gaspressure)
      heatflux     = -(gammabubble-1.0d0)/gammabubble*bubbletmp%betaT*(bubbleTemp-1.0d0)/bubbletmp%tmp%y(1)

      avgconc      = bubbletmp%tmp%mv/(bubbletmp%mg + bubbletmp%tmp%mv)
      Rmixt        = (concvapint/MWvap + (1.0d0 - concvapint)/MWgas)*Runiv

      concvapint   = min(concvapint, 0.99d0)
      vaporflux    = (1.0d0- concvapint)*bubbletmp%tmp%y(1)
      rhogas       = (bubbletmp%mg+bubbletmp%tmp%mv)/(4.0d0/3.0d0*pi*bubbletmp%tmp%y(1)**3)
      vaporflux    = -diffcoefvap*bubbletmp%betaC*(avgconc -concvapint)*rhogas/vaporflux

  END SUBROUTINE getfluxes



  FUNCTION deriv_gaspressure( bubbletmp, vaporflux, heatflux, gammabubble ) 

      TYPE(particledata) :: bubbletmp
      REAL(KIND(0.D0)) :: deriv_gaspressure, vaporflux,gammabubble, heatflux

      deriv_gaspressure = bubbletmp%tmp%p*bubbletmp%tmp%y(2) - heatflux - Runiv/MWvap*vaporflux
      deriv_gaspressure = -3.0d0*gammabubble/bubbletmp%tmp%y(1)*deriv_gaspressure 

  END FUNCTION deriv_gaspressure

  !--------------------------------------------------------------------------

  SUBROUTINE update(dt,RKstep,RKcoef, largestep,q,dq, q_prim )
      ! This subroutine update the particle variables

      TYPE(particlenode),POINTER                :: particle

      REAL(KIND(0.D0))                          :: dt
      REAL(KIND(0.D0)),DIMENSION(6),INTENT (IN) :: RKcoef
      INTEGER, INTENT (IN)                      :: RKstep
      INTEGER, DIMENSION(3)                     :: oldcell, newcell
      INTEGER                                   :: i,j
      LOGICAL                                   :: largestep,change,indomain
      TYPE(vector_field), DIMENSION(:), OPTIONAL :: q
      TYPE(vector_field), DIMENSION(:), OPTIONAL :: dq
      TYPE(scalar_field), DIMENSION(:), OPTIONAL :: q_prim
      TYPE (cellwb), POINTER        :: cellwbaux
      TYPE(particlenode),POINTER    :: nodeaux
      INTEGER,DIMENSION(3)          :: cell

      particle  => particlesubList%List%next

      DO WHILE(Associated(particle))

        oldcell =  get_cell_from_s(particle%data%tmp%s) 
        CALL transfertotmp (particle%data)

        DO i=1,RKstep
          particle%data%tmp%y(1:2) = particle%data%tmp%y(1:2) + dt*RKcoef(i)*particle%data%dbdt(i)%dydt(1:2)
          particle%data%tmp%x(1:3) = particle%data%tmp%x(1:3) + dt*RKcoef(i)*particle%data%dbdt(i)%dxdt(1:3)
          particle%data%tmp%u(1:3) = particle%data%tmp%u(1:3) + dt*RKcoef(i)*particle%data%dbdt(i)%dudt(1:3)
          particle%data%tmp%p      = particle%data%tmp%p      + dt*RKcoef(i)*particle%data%dbdt(i)%dpbdt
          particle%data%tmp%mv     = particle%data%tmp%mv     + dt*RKcoef(i)*particle%data%dbdt(i)%dmvdt
        ENDDO

        IF ((particle%data%tmp%y(1).LE.0.0d0).OR.(particle%data%tmp%x(1).NE.particle%data%tmp%x(1))) THEN
          !largestep=.TRUE.
          IF (dt.LT.2.d-15) THEN
            print *, 'warning large step',dt, particle%data%id
            CALL Remove_particle (particle,1)
            CALL Remove_particle (particle,2)
            GOTO 710
          ENDIF
          IF ((particle%data%tmp%x(1).NE.particle%data%tmp%x(1))) GOTO 711
        ENDIF

        indomain = particle_in_domain(particle%data%tmp%x)

        IF (.NOT.indomain) THEN 
          print *, 'not in domain', particle%data%id,particle%data%tmp%x(1),particle%data%xprev(1),x_cb(-buff_size-1),x_cb(n+buff_size)
          CALL Remove_particle (particle,2)
          GOTO 710
        ENDIF

        IF (particle%data%equilibrium) CALL equilibrium_state ( particle%data ,largestep, q_prim(E_idx) )

        particle => particle%next
        710 CONTINUE
      ENDDO

711   IF (num_procs > 1) THEN
         CALL bcst_largestep(largestep)
      ENDIF


      !IF (largestep) RETURN;
      !
      !   !update fluid variables
      IF (PRESENT(q)) THEN
        DO i = 1, sys_size
            q(2)%vf(i)%sf(0:m,0:n,0:p) = q(1)%vf(i)%sf(0:m,0:n,0:p) 
          DO j=1,RKstep
            q(2)%vf(i)%sf(0:m,0:n,0:p) = q(2)%vf(i)%sf(0:m,0:n,0:p) & 
                                                       + dt*RKcoef(j)*dq(j)%vf(i)%sf(0:m,0:n,0:p) 
          ENDDO
        ENDDO
      ENDIF

  END SUBROUTINE update


  SUBROUTINE equilibrium_state ( particle, largestep, pres )


      TYPE(particledata)     :: particle
      REAL(KIND(0.D0))       :: req(1)
      REAL(KIND(0.D0))       :: pinf, aux1, aux2
      INTEGER, DIMENSION(3)  :: cell
      LOGICAL                :: largestep,cond
      TYPE(scalar_field)     :: pres

      cell = get_cell_from_s(particle%tmp%s)
      pinf = Interpolate( particle%tmp%s, pres )
      
      req(1) = particle%tmp%y(1)
      cond = .FALSE.
      CALL get_equilibrium_radius (10000, req, 1.0d-10,1.0d-10, particle, pinf,cond)

      particle%tmp%y(1) = req(1)
      IF ((req(1).LE.0.0d0).OR.cond) THEN
        print *, 'released', particle%id
        particle%equilibrium = .FALSE.
        !largestep = .TRUE.
        RETURN
      ENDIF
      particle%tmp%p    = pinf - pvap + 2.0d0*sigmabubble/req(1)
      particle%tmp%mv   = pvap*4.0d0/3.0d0*pi*req(1)**3*MWvap/Runiv

  END SUBROUTINE equilibrium_state


      !-----------------------------------------------------------------------
  SUBROUTINE RKerror (timetmp, dt, RKcoef, errmax, largestep, t_step, q, q_prim, dq)
      ! This subroutine calculates the maximum error of the RK step

      TYPE(particlenode),POINTER                  :: particle
      REAL(KIND(0.D0)),INTENT (IN)              :: dt,timetmp
      REAL(KIND(0.D0)),DIMENSION(6),INTENT (IN) :: RKcoef
      REAL(KIND(0.D0))                          :: errmax,erraux,errb
      TYPE(vector_field), DIMENSION(:), OPTIONAL :: q
      TYPE(scalar_field), DIMENSION(:), OPTIONAL :: q_prim
      TYPE(vector_field), DIMENSION(:), OPTIONAL :: dq
      INTEGER :: i,j,k,l,l1,nb
      LOGICAL :: largestep
      INTEGER, INTENT(IN) :: t_step

      errmax = 0.0d0
      erraux = 0.0d0

      particle  => particlesubList%List%next
      DO WHILE(Associated(particle))
        errb   = 0.0d0
        IF (.NOT.particle%data%equilibrium) THEN
        !particle radius error
        DO i=1,6
          erraux = erraux + RKcoef(i)*particle%data%dbdt(i)%dydt(1)
        ENDDO
        errb=max(errb,abs(erraux)*dt/particle%data%R0)

        !interface velocity error
        erraux = 0.0d0
        DO i=1,6
          erraux = erraux + RKcoef(i)*particle%data%dbdt(i)%dydt(2)
        ENDDO
        errb=max(errb,abs(erraux)*dt)

        !particle velocity error
        DO j=1,3
          erraux = 0.0d0
          DO i=1,6
            erraux = erraux + RKcoef(i)*particle%data%dbdt(i)%dxdt(j)
          ENDDO
          errb=max(errb,abs(erraux)*dt/(abs(particle%data%tmp%u(j))+1.0d-4))
        ENDDO
        ENDIF
        errmax=max(errmax,errb)
        particle => particle%next
      ENDDO

      largestep = .FALSE.

      IF (PRESENT(q)) THEN 
        DO l1=1, cont_idx%end 
          DO k=0,p
            DO j=0,n
              DO i=0,m
                erraux = q(1)%vf(l1)%sf(i,j,k)
                DO l=1,6
                  erraux = erraux + dt*RKcoef(l)*dq(l)%vf(l1)%sf(i,j,k)
                ENDDO
                erraux = max(errmax,erraux)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO l1=mom_idx%beg, mom_idx%beg+DIM-1
          DO k=0,p
            DO j=0,n
              DO i=0,m
                erraux = q(1)%vf(l1)%sf(i,j,k)
                DO l=1,6
                  erraux = erraux + dt*RKcoef(l)*dq(l)%vf(l1)%sf(i,j,k)
                ENDDO
                erraux = max(errmax,erraux)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        !CALL s_compute_rhs(q(2)%vf,q_prim,dq(6)%vf,t_step=t_step,qtime=timetmp,largestep=largestep)
        CALL s_compute_rhs(q(2)%vf,q_prim,dq(6)%vf,t_step=t_step,qtime=timetmp)
        IF (num_procs > 1) CALL bcst_largestep(largestep)
      ENDIF

  END SUBROUTINE RKerror
      
      !-----------------------------------------------------------------------
 
  SUBROUTINE updateRK (q, update_fields, q_prim)

      TYPE(particlenode), POINTER          :: particle
      TYPE(vector_field), DIMENSION(:), OPTIONAL :: q
      TYPE(scalar_field), DIMENSION(:), OPTIONAL :: q_prim
      INTEGER                              :: i
      REAL(KIND(0.D0))                     :: pinf, pcrit
      LOGICAL                              :: update_fields, release_part

      particle  => particlesubList%List%next
      release_part = .FALSE.
      DO WHILE(Associated(particle)) 
          IF (.NOT.particle%data%equilibrium) THEN
            pinf = get_pinf1(particle%data%tmp%s, q_prim(E_idx), 1)
            pcrit = pvap - &
                    4.0d0*sigmabubble/(3.d0*sqrt(3.0d0*(pref+.0*sigmabubble/particle%data%R0)*particle%data%R0**3/(2.0d0*sigmabubble))) 
            pcrit = min(pcrit,-pref)
            IF (ABS((pcrit-pinf)/pcrit).LT.0.5d0) release_part = .TRUE.
          ENDIF
        CALL transfertodata (particle%data)
        particle => particle%next
      ENDDO

      ! releasing particles
      particle  => particlesubList%List%next
      IF (release_part) THEN
        DO WHILE(Associated(particle))
          particle%data%equilibrium = .FALSE.
          particle => particle%next
        ENDDO 
      ENDIF

      IF ((coupledflag.OR.bubblesources).AND.(update_fields)) THEN
        DO i=1,sys_size; q(1)%vf(i)%sf = q(2)%vf(i)%sf ;  ENDDO;
      ENDIF

      IF (avgdensFlag) THEN
        CALL voidfraction(q(1)%vf)
      ENDIF

  END SUBROUTINE updateRK

      !-----------------------------------------------------------------------


  SUBROUTINE potentials (q_prim,q,qtime)

    TYPE(scalar_field), DIMENSION(sys_size)   :: q_prim
    TYPE(scalar_field), DIMENSION(sys_size), OPTIONAL   :: q
    TYPE (cellwb), POINTER            :: cellwbaux
    INTEGER, DIMENSION(3)             :: cell 
    INTEGER                           :: info,i
    REAL(KIND(0.D0)),OPTIONAL         :: qtime

    cellwbaux => cellwbList%List%next
    DO WHILE(Associated(cellwbaux))
      cell = cellwbaux%data%coord
      IF(PRESENT(q)) THEN
         CALL solve_cell (cell,q_prim,q,qtime)
      ELSE
         CALL solve_cell (cell,q_prim)
      END IF
      cellwbaux => cellwbaux%next
    ENDDO

  END SUBROUTINE potentials


  SUBROUTINE solve_cell (cell,q_prim,q,qtime)

    TYPE(scalar_field), DIMENSION(sys_size)   :: q_prim
    TYPE(scalar_field), DIMENSION(sys_size), OPTIONAL   :: q
    INTEGER, DIMENSION(3)            :: cell
    REAL(KIND(0.D0))                 :: preterm1, term2, paux, pint, Romega, term1_fac, Rb
    REAL(KIND(0.D0)), DIMENSION(3)   :: scoord
    TYPE(particlenode), POINTER      :: particle
    REAL(KIND(0.D0)), OPTIONAL       :: qtime
    REAL(KIND(0.D0))                 :: rhol,cson
    REAL(KIND(0d0)), DIMENSION(2)              :: Re
    REAL(KIND(0d0)), DIMENSION( num_fluids, num_fluids ) :: We
 
    scoord(:) = cell(:) + 0.5d0
    paux = get_pinf1(scoord, q_prim(E_idx),2,preterm1,term2,Romega)
    particle => qbl%fp(cell(1),cell(2),cell(3))%List%next
  
    DO WHILE(Associated(particle))
      pint = pressureliq_int(particle%data%tmp%p,particle%data%tmp%y(1),particle%data%tmp%y(2)) + 0.5d0*particle%data%tmp%y(2)**2
      
      IF(clusterflag.EQ.3) THEN
 
         particle%Data%dphidt = (paux - pint) + term2
         
         ! Accouting for the potential induced by the bubble averaged over the control volume
         ! Note that this is based on the incompressible flow assumption near the bubble.
         Rb = particle%data%tmp%y(1)
         term1_fac=3.0d0/2.0d0*(Rb*(Romega**2d0-Rb**2d0))/(Romega**3d0-Rb**3d0)
         particle%Data%dphidt = particle%Data%dphidt/(1-term1_fac)

      END IF
  

      particle => particle%next

    ENDDO  

  END SUBROUTINE solve_cell

  SUBROUTINE s_populate_variables_buffers(v_vf, q_particle) ! ---------------
        ! Description: The purpose of this procedure is to populate the buffers
        !              of the conservative variables, depending on the selected
        !              boundary conditions. !Lagrangean particles

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: v_vf
            TYPE(scalar_field), DIMENSION(:), OPTIONAL :: q_particle

            ! Generic loop iterators
            INTEGER :: i,j

            ! Population of Buffers in x-direction =============================

            IF(bc_x%beg <= -3) THEN         ! Ghost-cell extrap. BC at beginning

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(-j,0:n,0:p) = &
                        v_vf(i)%sf( 0,0:n,0:p)
                    END DO
                END DO

            ELSEIF(bc_x%beg == -2) THEN     ! Symmetry BC at beginning

                DO j = 1, buff_size

                    DO i = 1, cont_idx%end
                        v_vf(i)%sf(-j ,0:n,0:p) = &
                        v_vf(i)%sf(j-1,0:n,0:p)
                    END DO

                    v_vf(mom_idx%beg)%sf(-j ,0:n,0:p) = &
                   -v_vf(mom_idx%beg)%sf(j-1,0:n,0:p)

                    DO i = mom_idx%beg+1, sys_size
                        v_vf(i)%sf(-j ,0:n,0:p) = &
                        v_vf(i)%sf(j-1,0:n,0:p)
                    END DO

                END DO

            ELSEIF(bc_x%beg == -1) THEN     ! Periodic BC at beginning

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(  -j   ,0:n,0:p) = &
                        v_vf(i)%sf(m-(j-1),0:n,0:p)
                    END DO
                END DO

            ELSE                            ! Processor BC at beginning

                IF(PRESENT(q_particle)) THEN
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=1, pbc_loc=-1, q_particle=q_particle )
                ELSE
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=1, pbc_loc=-1 )
                END IF

            END IF

            IF(bc_x%end <= -3) THEN         ! Ghost-cell extrap. BC at end

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(m+j,0:n,0:p) = &
                        v_vf(i)%sf( m ,0:n,0:p)
                    END DO
                END DO

            ELSEIF(bc_x%end == -2) THEN     ! Symmetry BC at end

                DO j = 1, buff_size

                    DO i = 1, cont_idx%end
                        v_vf(i)%sf(  m+j  ,0:n,0:p) = &
                        v_vf(i)%sf(m-(j-1),0:n,0:p)
                    END DO

                    v_vf(mom_idx%beg)%sf(  m+j  ,0:n,0:p) = &
                   -v_vf(mom_idx%beg)%sf(m-(j-1),0:n,0:p)

                    DO i = mom_idx%beg+1, sys_size
                        v_vf(i)%sf(  m+j  ,0:n,0:p) = &
                        v_vf(i)%sf(m-(j-1),0:n,0:p)
                    END DO

                END DO

            ELSEIF(bc_x%end == -1) THEN     ! Periodic BC at end

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(m+j,0:n,0:p) = &
                        v_vf(i)%sf(j-1,0:n,0:p)
                    END DO
                END DO

            ELSE                            ! Processor BC at end

                IF(PRESENT(q_particle)) THEN
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=1, pbc_loc=1, q_particle=q_particle )
                ELSE
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=1, pbc_loc=1 )
                END IF

            END IF

            ! END: Population of Buffers in x-direction ========================


            ! Population of Buffers in y-direction =============================

            IF(n == 0) THEN

                RETURN

            ELSEIF(bc_y%beg <= -3) THEN     ! Ghost-cell extrap. BC at beginning

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(:,-j,0:p) = &
                        v_vf(i)%sf(:, 0,0:p)
                    END DO
                END DO

            ELSEIF(bc_y%beg == -2) THEN     ! Symmetry BC at beginning

                DO j = 1, buff_size

                    DO i = 1, mom_idx%beg
                        v_vf(i)%sf(:,-j ,0:p) = &
                        v_vf(i)%sf(:,j-1,0:p)
                    END DO

                    v_vf(mom_idx%beg+1)%sf(:,-j ,0:p) = &
                   -v_vf(mom_idx%beg+1)%sf(:,j-1,0:p)

                    DO i = mom_idx%beg+2, sys_size
                        v_vf(i)%sf(:,-j ,0:p) = &
                        v_vf(i)%sf(:,j-1,0:p)
                    END DO

                END DO

            ELSEIF(bc_y%beg == -1) THEN     ! Periodic BC at beginning

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(:,  -j   ,0:p) = &
                        v_vf(i)%sf(:,n-(j-1),0:p)
                    END DO
                END DO

            ELSE                            ! Processor BC at beginning

                IF(PRESENT(q_particle)) THEN
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=2, pbc_loc=-1, q_particle=q_particle )
                ELSE
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=2, pbc_loc=-1 )
                END IF

            END IF

            IF(bc_y%end <= -3) THEN         ! Ghost-cell extrap. BC at end

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(:,n+j,0:p) = &
                        v_vf(i)%sf(:, n ,0:p)
                    END DO
                END DO

            ELSEIF(bc_y%end == -2) THEN     ! Symmetry BC at end

                DO j = 1, buff_size

                    DO i = 1, mom_idx%beg
                        v_vf(i)%sf(:,  n+j  ,0:p) = &
                        v_vf(i)%sf(:,n-(j-1),0:p)
                    END DO

                    v_vf(mom_idx%beg+1)%sf(:,  n+j  ,0:p) = &
                   -v_vf(mom_idx%beg+1)%sf(:,n-(j-1),0:p)

                    DO i = mom_idx%beg+2, sys_size
                        v_vf(i)%sf(:,  n+j  ,0:p) = &
                        v_vf(i)%sf(:,n-(j-1),0:p)
                    END DO

                END DO

            ELSEIF(bc_y%end == -1) THEN     ! Periodic BC at end

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(:,n+j,0:p) = &
                        v_vf(i)%sf(:,j-1,0:p)
                    END DO
                END DO

            ELSE                            ! Processor BC at end

                IF(PRESENT(q_particle)) THEN
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=2, pbc_loc=1, q_particle=q_particle )
                ELSE
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=2, pbc_loc=1 )
                END IF

            END IF

            ! END: Population of Buffers in y-direction ========================


            ! Population of Buffers in z-direction =============================

            IF(p == 0) THEN

                RETURN

            ELSEIF(bc_z%beg <= -3) THEN     ! Ghost-cell extrap. BC at beginning

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(:,:,-j) = &
                        v_vf(i)%sf(:,:, 0)
                    END DO
                END DO

            ELSEIF(bc_z%beg == -2) THEN     ! Symmetry BC at beginning

                DO j = 1, buff_size

                    DO i = 1, mom_idx%beg+1
                        v_vf(i)%sf(:,:,-j ) = &
                        v_vf(i)%sf(:,:,j-1)
                    END DO

                    v_vf(mom_idx%end)%sf(:,:,-j ) = &
                   -v_vf(mom_idx%end)%sf(:,:,j-1)

                    DO i = E_idx, sys_size
                        v_vf(i)%sf(:,:,-j ) = &
                        v_vf(i)%sf(:,:,j-1)
                    END DO

                END DO

            ELSEIF(bc_z%beg == -1) THEN     ! Periodic BC at beginning

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(:,:,  -j   ) = &
                        v_vf(i)%sf(:,:,p-(j-1))
                    END DO
                END DO

            ELSE                            ! Processor BC at beginning

                IF(PRESENT(q_particle)) THEN
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=3, pbc_loc=-1, q_particle=q_particle )
                ELSE
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=3, pbc_loc=-1 )
                END IF

            END IF

            IF(bc_z%end <= -3) THEN         ! Ghost-cell extrap. BC at end

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(:,:,p+j) = &
                        v_vf(i)%sf(:,:, p )
                    END DO
                END DO

            ELSEIF(bc_z%end == -2) THEN     ! Symmetry BC at end

                DO j = 1, buff_size

                    DO i = 1, mom_idx%beg+1
                        v_vf(i)%sf(:,:,  p+j  ) = &
                        v_vf(i)%sf(:,:,p-(j-1))
                    END DO

                    v_vf(mom_idx%end)%sf(:,:,  p+j  ) = &
                   -v_vf(mom_idx%end)%sf(:,:,p-(j-1))

                    DO i = E_idx, sys_size
                        v_vf(i)%sf(:,:,  p+j  ) = &
                        v_vf(i)%sf(:,:,p-(j-1))
                    END DO

                END DO

            ELSEIF(bc_z%end == -1) THEN     ! Periodic BC at end

                DO i = 1, sys_size
                    DO j = 1, buff_size
                        v_vf(i)%sf(:,:,p+j) = &
                        v_vf(i)%sf(:,:,j-1)
                    END DO
                END DO

            ELSE                            ! Processor BC at end

                IF(PRESENT(q_particle)) THEN
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=3, pbc_loc=1, q_particle=q_particle )
                ELSE
                    CALL s_mpi_sendrecv_conservative_variables_buffers( &
                                              v_vf, mpi_dir=3, pbc_loc=1 )
                END IF

            END IF

            ! END: Population of Buffers in z-direction ========================


        END SUBROUTINE s_populate_variables_buffers ! -------------


END MODULE m_particles
