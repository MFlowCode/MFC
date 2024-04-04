! pMFC_v3.0 - Simulation code: m_mpi_particles.f90
! Description: The module contains subroutines used for MPI for particles
! Author: Kazuki Maeda
! Date: 01/01/17

MODULE m_mpi_particles
     
     
        ! Dependencies =============================================================
#ifdef MFC_MPI
        USE mpi
#endif
        USE m_mpi_proxy                 ! Message passing interface (MPI) module
     
        USE m_global_parameters         ! Global parameters for the code

        USE m_particles_types           ! particle structures

        ! ==========================================================================

        IMPLICIT NONE 

        TYPE proc_Lists
            TYPE (particleListinfo),POINTER  :: List
            INTEGER                          :: target_proc
        END TYPE proc_Lists
        TYPE particletypetransf
            INTEGER                          :: yi
            REAL(KIND(0.D0)), DIMENSION(17)  :: yr
        END TYPE particletypetransf
        TYPE (proc_Lists),ALLOCATABLE,DIMENSION(:)   :: transf_part_Lists
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: id_List
        INTEGER :: n_neighbors
        INTEGER :: MPI_Req

        ! Variable to contain particle variables for Parallel I/O
        REAL(KIND(0.D0)), ALLOCATABLE, DIMENSION(:,:)     :: MPI_IO_DATA_particle
        ! Number of variables contained in MPI_IO_DATA_particle per particle
        INTEGER :: mpi_par_var

        integer, private :: ierr
        !INTEGER :: MPI_COMM_CART

CONTAINS

  SUBROUTINE reduce_dtnext( dtnext )
  
  REAL(KIND(0.D0)) :: dtnext, dtmin

#ifdef MFC_MPI
  !dtnext = min(dt0, dtnext) 
  dtnext = min(dt, dtnext)
  CALL MPI_ALLREDUCE (dtnext, dtmin, 1 , MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  dtnext = dtmin
#endif

  END SUBROUTINE reduce_dtnext


  SUBROUTINE s_mpi_bcast_user_particles()

#ifdef MFC_MPI
   CALL MPI_BCAST(do_particles, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(particleflag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(avgdensFlag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(particleoutFlag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(particlestatFlag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(RPflag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(clusterflag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(stillparticlesflag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
!   CALL MPI_BCAST(Ffluidflag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(heatflag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(massflag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(Runiv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(gammagas, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(gammavapor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(pvap, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(cpgas, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(cpvapor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(kgas, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(kvapor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(MWgas, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(MWvap, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(diffcoefvap, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(sigmabubble, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(viscref, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(RKeps, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(ratiodt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(projectiontype, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(smoothtype, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(epsilonb, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(coupledflag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(solverapproach, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(correctpresFlag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(charwidth, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(valmaxvoid, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(dtmaxpart, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(old_2D, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

  END SUBROUTINE s_mpi_bcast_user_particles


  SUBROUTINE initialize_particles_mpi

#ifdef MFC_MPI      
    INTEGER :: oldtypes(0:1), blockcounts(0:1), sizeint, sizedble
    INTEGER (kind=MPI_ADDRESS_KIND) :: offsets(0:1)
    INTEGER, ALLOCATABLE, DIMENSION(:) :: proc_coords
    INTEGER :: i,j,n_max

    !fixme: don't know if this definition of sizeint and sizedble are correct
    sizeint = 0
    sizedble = 0

    CALL MPI_TYPE_EXTENT(MPI_INTEGER, sizeint, ierr)

    offsets(0)     = 0
    oldtypes(0)    = MPI_INTEGER
    blockcounts(0) = 1

    CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, sizedble, ierr)
    offsets(1)     = sizeint
    oldtypes(1)    = MPI_DOUBLE_PRECISION
    blockcounts(1) = 17

    CALL MPI_TYPE_CREATE_STRUCT(2, blockcounts, offsets, oldtypes,particletype_mpi, ierr)
    CALL MPI_TYPE_COMMIT(particletype_mpi, ierr)

    ! getting the buffer size required by one particle
    CALL MPI_Pack_size (1, particletype_mpi, MPI_COMM_WORLD, unitbuffersize,ierr)

    !generating structures to transfer particles 
    IF (DIM.EQ.2) THEN

      ALLOCATE(proc_coords(1:2))
      proc_coords = 0
      n_max=8
      ! global position, position in the sorted list
      ALLOCATE(id_List(1:n_max,1:2))
      id_List(:,2) = -1
!     list order
!              4   3   2
!              5       1
!              6   7   8
!
      CALL MPI_CART_COORDS( MPI_COMM_CART, proc_rank, 3,proc_coords, ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(1,1), ierr)
      proc_coords(2) = proc_coords(2)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(2,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(3,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(4,1), ierr)
      proc_coords(2) = proc_coords(2)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(5,1), ierr)
      proc_coords(2) = proc_coords(2)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(6,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(7,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(8,1), ierr)
    ELSE  

      ALLOCATE(proc_coords(1:3))
      proc_coords = 0
      n_max=26
      ALLOCATE(id_List(1:n_max,1:2))
      id_List(:,2) = -1
      CALL MPI_CART_COORDS( MPI_COMM_CART, proc_rank, 3,proc_coords, ierr)
!     list order
!                9- 10 - 11      9  10 11          9 -10 - 11
!               /   /   /    <-  12 13 14         /   /   / |
!              4   3   2         15 16 17        4 - 3 - 2  |
!              5       1                        /   /   / | |
!              6   7   8                       18- 19- 20|/
!                                              21- 22- 23/
!                          <------------       24- 25- 26
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(1,1), ierr)
      proc_coords(2) = proc_coords(2)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(2,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(3,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(4,1), ierr)
      proc_coords(2) = proc_coords(2)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(5,1), ierr)
      proc_coords(2) = proc_coords(2)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(6,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(7,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(8,1), ierr)
      proc_coords(3) = proc_coords(3)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(26,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(25,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(24,1), ierr)
      proc_coords(2) = proc_coords(2)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(21,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(22,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(23,1), ierr)
      proc_coords(2) = proc_coords(2)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(20,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(19,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(18,1), ierr)
      proc_coords(3) = proc_coords(3)-2
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(9,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(10,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(11,1), ierr)
      proc_coords(2) = proc_coords(2)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(14,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(13,1), ierr)
      proc_coords(1) = proc_coords(1)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(12,1), ierr)
      proc_coords(2) = proc_coords(2)-1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(15,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(16,1), ierr)
      proc_coords(1) = proc_coords(1)+1
      CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, id_List(17,1), ierr)
    ENDIF


    DEALLOCATE(proc_coords)

    !count different neighbors
    n_neighbors = 0

    DO i=1,n_max
      j=1
      DO WHILE ((id_List(i,1).NE.id_List(j,1)))
        j = j + 1
      ENDDO
      IF (id_List(j,2).EQ.-1) THEN
        n_neighbors = n_neighbors + 1
        id_List(i,2)=n_neighbors
      ELSE
        id_List(i,2)=j
      ENDIF
    ENDDO

    
    ALLOCATE(transf_part_Lists(n_neighbors)) 
    DO i=1,n_neighbors
        j=1
        DO WHILE (i.NE.id_List(j,2))
          j = j +1
        ENDDO
        transf_part_Lists(i)%target_proc = id_List(j,1)
        ALLOCATE(transf_part_Lists(i)%List)
        CALL init_particle_list(transf_part_Lists(i)%List)
    ENDDO
#endif

  END SUBROUTINE initialize_particles_mpi


!  SUBROUTINE s_initialize_mpi_data_particle(particleinfo) ! --------------------------


            ! Particle variables
!            TYPE(particlenode), POINTER  :: particleinfo

            ! Number of particles on the MPI_process
!            INTEGER :: particle_num

            !INTEGER, DIMENSION(num_dims) :: gsizes, lsizes
            !INTEGER :: ierr, alt_sys

            ! Generic loop iterator
            !INTEGER :: i

!            MPI_IO_DATA_particle(1) = particleinfo%id
!            MPI_IO_DATA_particle(2:4) = particleinfo%x(1:3)
!            MPI_IO_DATA_particle(5:7) = particleinfo%xprev(1:3)
!            MPI_IO_DATA_particle(8:10) = particleinfo%u(1:3)
!            MPI_IO_DATA_particle(11:12) = particleinfo%y(1:2)
!            MPI_IO_DATA_particle(13) = particleinfo%R0
!            MPI_IO_DATA_particle(14) = particleinfo%Rmax
!            MPI_IO_DATA_particle(15) = particleinfo%Rmin
!            MPI_IO_DATA_particle(16) = particleinfo%dphidt
!            MPI_IO_DATA_particle(17) = particleinfo%p
!            MPI_IO_DATA_particle(18) = particleinfo%mv
!            MPI_IO_DATA_particle(19) = particleinfo%mg
!            MPI_IO_DATA_particle(20) = particleinfo%betaT
!            MPI_IO_DATA_particle(21) = particleinfo%betaC

           ! Define the view for each variable
           ! DO i = 1, alt_sys
           !     CALL MPI_TYPE_CREATE_SUBARRAY(num_dims,gsizes,lsizes,start_idx,&
           !         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i),ierr)
           !     CALL MPI_TYPE_COMMIT(MPI_IO_DATA%view(i),ierr)
           ! END DO


!  END SUBROUTINE s_initialize_mpi_data_particle ! ---------------------------------


  SUBROUTINE get_min( aux )

  REAL(KIND(0.D0)) :: aux,auxmin

#ifdef MFC_MPI
    CALL MPI_ALLREDUCE (aux, auxmin, 1 , MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    aux = auxmin
#endif

  END SUBROUTINE get_min
  
  SUBROUTINE get_max( aux )

  REAL(KIND(0.D0)) :: aux,auxmax

#ifdef MFC_MPI
    CALL MPI_ALLREDUCE (aux, auxmax, 1 , MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    aux = auxmax
#endif

  END SUBROUTINE get_max

  SUBROUTINE get_sum( aux )

  REAL(KIND(0.D0)) :: aux,auxmax

#ifdef MFC_MPI
    CALL MPI_ALLREDUCE (aux, auxmax, 1 , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    aux = auxmax
#endif

  END SUBROUTINE get_sum
  
  SUBROUTINE bcst_largestep(largestep)

  LOGICAL :: largestep,aux

#ifdef MFC_MPI
    CALL MPI_ALLREDUCE (largestep, aux, 1 , MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    largestep = aux
#endif

  END SUBROUTINE bcst_largestep
!--------------------------------------------------------------------------------------------

SUBROUTINE transfer_particles

#ifdef MFC_MPI
    IF (p > 0 ) THEN
      CALL add_to_transfer_list3D
    ELSE
      CALL add_to_transfer_list2D
    ENDIF
    CALL send_particles
    CALL rcv_particles ()
#endif

END SUBROUTINE transfer_particles

!--------------------------------------------------------------------------------------------

SUBROUTINE add_to_transfer_list2D

! list order
!              4   3   2
!              5       1
!              6   7   8
!
!
    TYPE (cellwb), POINTER          :: cellwbaux
    TYPE (particlenode),POINTER     :: node
    
    cellwbaux => cellwbList%List%next
    DO WHILE(Associated(cellwbaux))
       !print *, "mpi",proc_rank,cellwbList%nb
       node  => qbl%fp(cellwbaux%data%coord(1), cellwbaux%data%coord(2), cellwbaux%data%coord(3))%List%next
       DO WHILE(Associated(node))
          !x-right direction
          !print *, "mpi",proc_rank,node%data%xprev(1),node%data%x(1),x_cb(m-buff_size)
          IF ((node%data%xprev(1).LT.x_cb(m-buff_size)).AND.(node%data%x(1).GT.x_cb(m-buff_size))) THEN
!            print *, 'transfer x',id_List(1,2)
            CALL addparticletolist(node%data,transf_part_Lists(id_List(1,2))%List)

            IF ((node%data%x(2).GT.y_cb(n-buff_size)).AND.(node%data%xprev(2).LT.y_cb(n-buff_size))) THEN
!                print *, 'transfer y1', id_List(3,2)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(3,2))%List)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(2,2))%List)
            ELSEIF ((node%data%x(2).LT.y_cb(buff_size)).AND.(node%data%xprev(2).GT.y_cb(buff_size))) THEN
!                print *, 'transfer y2', id_List(7,2)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(7,2))%List)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(8,2))%List)
            ENDIF

            GOTO 525
          END IF

          !x-left direction
          IF ((node%data%xprev(1).GT.x_cb(buff_size)).AND.(node%data%x(1).LT.x_cb(buff_size))) THEN
!            print *, 'transfer x2', id_List(5,2)
            CALL addparticletolist(node%data,transf_part_Lists(id_List(5,2))%List)

            IF ((node%data%x(2).GT.y_cb(n-buff_size)).AND.(node%data%xprev(2).LT.y_cb(n-buff_size))) THEN
!                print *, 'transfer y3'
                CALL addparticletolist(node%data,transf_part_Lists(id_List(4,2))%List)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(3,2))%List)
            ELSEIF ((node%data%x(2).LT.y_cb(buff_size)).AND.(node%data%xprev(2).GT.y_cb(buff_size))) THEN
!                print *, 'transfer y4'
                CALL addparticletolist(node%data,transf_part_Lists(id_List(6,2))%List)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(7,2))%List)
            ENDIF

            GOTO 525
          END IF
          
          !y-up direction (sides already accounted for)
          IF ((node%data%xprev(2).LT.y_cb(n-buff_size)).AND.(node%data%x(2).GT.y_cb(n-buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(3,2))%List)
          ENDIF
          
          !y-down direction (sides already accounted for)
          IF ((node%data%xprev(2).GT.y_cb(buff_size)).AND.(node%data%x(2).LT.y_cb(buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(7,2))%List)
          ENDIF

525       CONTINUE
          node%data%xprev = node%data%x
          node => node%next
        ENDDO
        cellwbaux => cellwbaux%next
     ENDDO


END SUBROUTINE add_to_transfer_list2D

SUBROUTINE add_to_transfer_list3D
! list order

!                9- 10 - 11      9  10 11          9 -10 - 11
!               /   /   /    <-  12 13 14         /   /   / |
!              4   3   2         15 16 17        4 - 3 - 2  |
!              5       1                        /   /   / | |
!            / 6   7   8                       18- 19- 20|/
!           /         /                        21- 22- 23/
!              FRONT     <------------>        24- 25- 26

    TYPE (cellwb), POINTER          :: cellwbaux
    TYPE (particlenode),POINTER     :: node
    
    cellwbaux => cellwbList%List%next
    DO WHILE(Associated(cellwbaux))
       node  => qbl%fp(cellwbaux%data%coord(1), cellwbaux%data%coord(2), cellwbaux%data%coord(3))%List%next
       DO WHILE(Associated(node))
          !x-right direction
          IF ((node%data%xprev(1).LT.x_cb(m-buff_size)).AND.(node%data%x(1).GT.x_cb(m-buff_size))) THEN
            CALL addparticletolist(node%data,transf_part_Lists(id_List(1,2))%List)

            !z direction
            IF ((node%data%x(3).GT.z_cb(p-buff_size)).AND.(node%data%xprev(3).LT.z_cb(p-buff_size))) THEN
              CALL addparticletolist(node%data,transf_part_Lists(id_List(23,2))%List)
            ELSEIF ((node%data%x(3).LT.z_cb(buff_size)).AND.(node%data%xprev(3).GT.z_cb(buff_size))) THEN
              CALL addparticletolist(node%data,transf_part_Lists(id_List(14,2))%List)
            ENDIF

            IF ((node%data%x(2).GT.y_cb(n-buff_size)).AND.(node%data%xprev(2).LT.y_cb(n-buff_size))) THEN
              CALL addparticletolist(node%data,transf_part_Lists(id_List(2,2))%List)
              CALL addparticletolist(node%data,transf_part_Lists(id_List(3,2))%List)

              !z-direction
              IF ((node%data%x(3).GT.z_cb(p-buff_size)).AND.(node%data%xprev(3).LT.z_cb(p-buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(20,2))%List)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(19,2))%List)
              ELSEIF ((node%data%x(3).LT.z_cb(buff_size)).AND.(node%data%xprev(3).GT.z_cb(buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(11,2))%List)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(10,2))%List)
              ENDIF

            ENDIF

            IF ((node%data%x(2).LT.y_cb(buff_size)).AND.(node%data%xprev(2).GT.y_cb(buff_size))) THEN
              CALL addparticletolist(node%data,transf_part_Lists(id_List(7,2))%List)
              CALL addparticletolist(node%data,transf_part_Lists(id_List(8,2))%List)
                
              !z-direction
              IF ((node%data%x(3).GT.z_cb(p-buff_size)).AND.(node%data%xprev(3).LT.z_cb(p-buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(26,2))%List)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(25,2))%List)
              ELSEIF ((node%data%x(3).LT.z_cb(buff_size)).AND.(node%data%xprev(3).GT.z_cb(buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(16,2))%List)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(17,2))%List)
              ENDIF

            ENDIF
            GOTO 526
          END IF

          !x-left direction
          IF ((node%data%xprev(1).GT.x_cb(buff_size)).AND.(node%data%x(1).LT.x_cb(buff_size))) THEN
            CALL addparticletolist(node%data,transf_part_Lists(id_List(5,2))%List)

            !z-direction
            IF (node%data%x(3).GT.z_cb(p-buff_size).AND.(node%data%xprev(3).LT.z_cb(p-buff_size))) THEN
               CALL addparticletolist(node%data,transf_part_Lists(id_List(21,2))%List)
            ELSEIF ((node%data%x(3).LT.z_cb(buff_size)).AND.(node%data%xprev(3).GT.z_cb(buff_size))) THEN
               CALL addparticletolist(node%data,transf_part_Lists(id_List(12,2))%List)
            ENDIF

            IF ((node%data%x(2).GT.y_cb(n-buff_size)).AND.(node%data%xprev(2).LT.y_cb(n-buff_size))) THEN
               CALL addparticletolist(node%data,transf_part_Lists(id_List(3,2))%List)
               CALL addparticletolist(node%data,transf_part_Lists(id_List(4,2))%List)

                !z-direction
                IF (node%data%x(3).GT.z_cb(p-buff_size).AND.(node%data%xprev(3).LT.z_cb(p-buff_size))) THEN
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(18,2))%List)
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(19,2))%List)
                ELSEIF ((node%data%x(3).LT.z_cb(buff_size)).AND.(node%data%xprev(3).GT.z_cb(buff_size))) THEN
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(10,2))%List)
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(9,2))%List)
                ENDIF

             ELSEIF ((node%data%x(2).LT.y_cb(buff_size)).AND.(node%data%xprev(2).GT.y_cb(buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(7,2))%List)
                CALL addparticletolist(node%data,transf_part_Lists(id_List(6,2))%List)

                !z-direction
                IF ((node%data%x(3).GT.z_cb(p-buff_size)).AND.(node%data%xprev(3).LT.z_cb(p-buff_size))) THEN
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(24,2))%List)
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(25,2))%List)
                ELSEIF ((node%data%x(3).LT.z_cb(buff_size)).AND.(node%data%xprev(3).GT.z_cb(buff_size))) THEN
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(15,2))%List)
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(16,2))%List)
                ENDIF

            ENDIF
            GOTO 526
          END IF
          
          !y-up direction (sides already accounted for)
          IF ((node%data%xprev(2).LT.y_cb(n-buff_size)).AND.(node%data%x(2).GT.y_cb(n-buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(3,2))%List)

                !z-direction
                IF (node%data%x(3).GT.z_cb(p-buff_size).AND.(node%data%xprev(3).LT.z_cb(p-buff_size))) THEN
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(19,2))%List)
                ELSEIF ((node%data%x(3).LT.z_cb(buff_size)).AND.(node%data%xprev(3).GT.z_cb(buff_size))) THEN
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(10,2))%List)
                ENDIF

          ENDIF
          
          !y-down direction (sides already accounted for)
          IF ((node%data%xprev(2).GT.y_cb(buff_size)).AND.(node%data%x(2).LT.y_cb(buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(7,2))%List)

                !z-direction
                IF ((node%data%x(3).GT.z_cb(p-buff_size)).AND.(node%data%xprev(3).LT.z_cb(p-buff_size))) THEN
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(25,2))%List)
                ELSEIF ((node%data%x(3).LT.z_cb(buff_size)).AND.(node%data%xprev(3).GT.z_cb(buff_size))) THEN
                  CALL addparticletolist(node%data,transf_part_Lists(id_List(16,2))%List)
                ENDIF
          ENDIF
          
          !z-front
          IF ((node%data%xprev(3).LE.z_cb(p-buff_size)).AND.(node%data%x(3).GT.z_cb(p-buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(22,2))%List)
          !z-back
          ELSEIF ((node%data%xprev(3).GT.z_cb(buff_size)).AND.(node%data%x(3).LT.z_cb(buff_size))) THEN
                CALL addparticletolist(node%data,transf_part_Lists(id_List(13,2))%List)
          ENDIF

526       CONTINUE
          node%data%xprev = node%data%x
          node => node%next
        ENDDO
        cellwbaux => cellwbaux%next
     ENDDO
         

END SUBROUTINE add_to_transfer_list3D

SUBROUTINE send_particles

TYPE (particletypetransf), ALLOCATABLE    :: sendparticles(:)
TYPE(particletypetransf) :: senddata
TYPE(particlenode),POINTER :: particle
INTEGER :: tag, sendbytepos, buffersend, i

#ifdef MFC_MPI
tag = 0
!sender code
DO i=1,n_neighbors
  CALL MPI_ISEND ( transf_part_Lists(i)%List%nb, 1, MPI_INTEGER , transf_part_Lists(i)%target_proc, tag, MPI_COMM_WORLD, MPI_Req, ierr )
  !sending particles
  IF (transf_part_Lists(i)%List%nb.NE.0) THEN
    ALLOCATE(sendparticles(transf_part_Lists(i)%List%nb))
    buffersend=unitbuffersize*transf_part_Lists(i)%List%nb
    sendbytepos = 0
    particle  => transf_part_Lists(i)%List%List%next
    !packing
    DO WHILE(Associated(particle))
      !I do not understand why I have to transfer the info to senddata... but
      !it does not work if I try to transfer particleaux%P%Data directly (I think I have to create the proper type)
      CALL transfer2senddata(senddata,particle)

      CALL MPI_Pack (senddata, 1, particletype_mpi, sendparticles, buffersend, sendbytepos, MPI_COMM_WORLD, ierr)
      particle => particle%next
    ENDDO
    !sending
    CALL MPI_ISEND(sendparticles, sendbytepos, MPI_PACKED, transf_part_Lists(i)%target_proc, tag, MPI_COMM_WORLD, MPI_Req, ierr )
    DEALLOCATE (sendparticles)
    CALL Clear_particle_list(transf_part_Lists(i)%List)
  ENDIF
ENDDO
#endif

END SUBROUTINE send_particles

SUBROUTINE rcv_particles ( )

#ifdef MFC_MPI
TYPE (particletypetransf), ALLOCATABLE    :: recvparticles(:)
TYPE(particletypetransf) :: recvdata
INTEGER :: tag, nrcv_particles, recvbytepos, bufferrecv,bytesrecv,i
INTEGER :: stat(MPI_STATUS_SIZE)


tag = 0

!receiver code
DO i=1,n_neighbors
  !receiving particles
  CALL MPI_RECV  ( nrcv_particles, 1, MPI_INTEGER , transf_part_Lists(i)%target_proc, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
  IF (nrcv_particles.NE.0) THEN
   ALLOCATE(recvparticles(nrcv_particles))
  !receiver code
    bufferrecv=unitbuffersize*nrcv_particles
    CALL MPI_RECV(recvparticles, bufferrecv, MPI_PACKED, transf_part_Lists(i)%target_proc, tag,  MPI_COMM_WORLD, stat, ierr )
    CALL MPI_Get_Count (stat,MPI_PACKED,bytesrecv,ierr)
    recvbytepos=0
    !unpacking and adding to the list
    DO WHILE(recvbytepos.LT.bytesrecv)
          CALL Mpi_unpack ( recvparticles, bytesrecv, recvbytepos, recvdata, 1, particletype_mpi, MPI_COMM_WORLD,ierr )
          CALL  transfer2binfo(recvdata)
    ENDDO
    DEALLOCATE (recvparticles)
  ENDIF
ENDDO
#endif

END SUBROUTINE rcv_particles

!!---------------------------------------------------------
SUBROUTINE transfer2senddata(senddata,particle)

    TYPE(particletypetransf) :: senddata
    TYPE(particlenode),POINTER :: particle

         senddata%yi        = particle%Data%id
         senddata%yr(1:3)   = particle%Data%x(:)
         senddata%yr(4:6)   = particle%Data%u(:)
         senddata%yr(7)     = particle%Data%Rmax
         senddata%yr(8)     = particle%Data%Rmin
         senddata%yr(9)     = particle%Data%dphidt
         senddata%yr(10:11) = particle%Data%y(:)
         senddata%yr(12)    = particle%Data%R0
         senddata%yr(13)    = particle%Data%p
         senddata%yr(14)    = particle%Data%mg
         senddata%yr(15)    = particle%Data%mv
         senddata%yr(16)    = particle%Data%betaC
         senddata%yr(17)    = particle%Data%betaT

END SUBROUTINE transfer2senddata

SUBROUTINE transfer2binfo(recvdata)

    TYPE(particletypetransf) :: recvdata
    TYPE(particledata), POINTER :: particleinfo
    TYPE(particleListinfo), POINTER :: particleListaux
    INTEGER, DIMENSION(3) :: cell
    LOGICAL               :: indomain

         ALLOCATE(particleinfo)
         particleinfo%id      = recvdata%yi     
         particleinfo%x(:)    = recvdata%yr(1:3)
         particleinfo%xprev   = particleinfo%x
         particleinfo%u(:)    = recvdata%yr(4:6)   
         particleinfo%y(1:2)  = recvdata%yr(10:11)
         particleinfo%R0      = recvdata%yr(12)     
         particleinfo%p       = recvdata%yr(13)   
         particleinfo%mg      = recvdata%yr(14)    
         particleinfo%mv      = recvdata%yr(15)
         particleinfo%betaC   = recvdata%yr(16)
         particleinfo%betaT   = recvdata%yr(17)
         particleinfo%Rmax      = recvdata%yr(7)
         particleinfo%Rmin      = recvdata%yr(8)
         particleinfo%dphidt    = recvdata%yr(9)
         
         particleinfo%equilibrium = .FALSE.

         indomain = particle_in_domain(particleinfo%x)
  
         !relocating particle from periodic bc
         IF (.NOT.indomain) THEN
           IF (particleinfo%x(1).GT.x_cb(m+buff_size)) THEN
             particleinfo%x(1) = x_cb(-buff_size-1)
           ELSEIF (particleinfo%x(1).LT.x_cb(-buff_size-1)) THEN
            particleinfo%x(1) = x_cb(m+buff_size-1)
           ENDIF

           IF (particleinfo%x(2).GT.y_cb(n+buff_size)) THEN
             particleinfo%x(2) = y_cb(-buff_size-1)
           ELSEIF (particleinfo%x(2).LT.y_cb(-buff_size-1)) THEN
            particleinfo%x(2) = y_cb(n+buff_size-1)
          ENDIF

           IF (DIM.eq.3) THEN
             IF (particleinfo%x(3).GT.z_cb(p+buff_size)) THEN
               particleinfo%x(3) = z_cb(-buff_size-1)
             ELSEIF (particleinfo%x(3).LT.z_cb(-buff_size-1)) THEN 
              !particleinfo%x(3) = z_cb(p+buff_size)
              particleinfo%x(3) = z_cb(p+buff_size-1)
            ENDIF
           ENDIF
           particleinfo%xprev = particleinfo%x
         ENDIF

         cell(:) = -buff_size
         CALL locate_cell ( particleinfo%x,  cell, particleinfo%tmp%s )
         particleListaux => qbl%fp(cell(1),cell(2),cell(3))
         IF (particleListaux%nb.EQ.0) CALL addtocell_list ( cell )
         CALL transfertotmp (particleinfo) 
         CALL addparticletolist (particleinfo,particleListaux)
         CALL addparticletolist (particleinfo,particlesubList)

END SUBROUTINE transfer2binfo


END MODULE m_mpi_particles

