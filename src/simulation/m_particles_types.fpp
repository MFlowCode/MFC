! pMFC v3.0 - Simulation Code: m_particles_types.f90
! Description: This module features the definitions of all the custom-defined
!              derived types that are utilized throughout the simulation code
!              for particles.
! Author: Kazuki Maeda
! Date: 01/01/17

MODULE m_particles_types
  
  ! Dependencies =============================================================
  USE m_global_parameters
  USE m_derived_types
  ! ==========================================================================
 
  IMPLICIT NONE

  TYPE cellwbcoord
    INTEGER, DIMENSION(3) :: coord
  END TYPE cellwbcoord

  TYPE cellwb
      TYPE (cellwb), POINTER     ::  next, prev
      TYPE (cellwbcoord),POINTER ::  data
  END TYPE cellwb
    
  TYPE dirlist
      TYPE (dirlist), POINTER     ::  next
      INTEGER                     ::  dir
  END TYPE dirlist
  
  TYPE cellListinfo  
        TYPE(cellwb), POINTER :: List
        INTEGER               :: nb ! number of cells in the list. It is useful to construct the sublists 
  END TYPE cellListinfo

  TYPE particlederivative
       REAL(KIND(0.D0)), DIMENSION(3)   :: dxdt, dudt, dMdt
       REAL(KIND(0.D0)), DIMENSION(2)   :: dydt        
       REAL(KIND(0.D0))                 :: dpbdt, dmvdt, dphidt
  END TYPE particlederivative

  TYPE particletmp                        ! if this list is modified,transfertotmp has to be also modified 
     REAL(KIND(0.D0)), DIMENSION(3)   :: x, s, u ! x: real coord, s: comp coord, u: vel of the particle
     REAL(KIND(0.D0)), DIMENSION(2)   :: y ! y(1): radius, y(2): radial velocity
     REAL(KIND(0.D0))                 :: p, mv
  END TYPE particletmp

  TYPE particledata
       INTEGER                          :: id
       REAL(KIND(0.D0)), DIMENSION(3)   :: x, xprev      !physical and computational position
       REAL(KIND(0.D0)), DIMENSION(3)   :: u             !physical particle velocity
       REAL(KIND(0.D0)), DIMENSION(2)   :: y             !particle variables (rb,drbdt)
       REAL(KIND(0.D0))                 :: R0, p, mg, mv ! initial values (mg = mass of noncondensable gas, Cvap: mass vapor fraction)
       REAL(KIND(0.D0))                 :: betaC, betaT, dphidt
       REAL(KIND(0.D0))                 :: Rmax, Rmin !statistical data
       LOGICAL                          :: equilibrium
       TYPE (particletmp)               :: tmp        !temporal variable for intermediate steps
       TYPE (particlederivative)        :: dbdt(6)    !derivatives. It could be an allocable pointer
  END TYPE particledata
  
  TYPE particlenode   ! This structure is just required if we want to create different lists to the same elements
      TYPE (particlenode), POINTER     :: next, prev !pointer to the next element
      TYPE (particledata), POINTER     :: data
  ENDTYPE particlenode
  
  TYPE particleListinfo  
        TYPE(particlenode), POINTER     :: List
        INTEGER                         :: nb    !number of particles in the list. It is useful to construct the cell list
        TYPE(cellwb), POINTER           :: cellpointer !pointer to the list of the cells (to speed up the process of updating)
        TYPE(particleListinfo), POINTER :: next, prev
  END TYPE particleListinfo
   
  TYPE :: list3D
       TYPE(particleListinfo), ALLOCATABLE :: fp(:,:,:)
  END TYPE list3D

  TYPE(list3D),TARGET              :: qbl
  TYPE(cellListinfo), POINTER      :: cellwbList

  REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: x_cc_lp, y_cc_lp, z_cc_lp

  REAL(KIND(0.D0)) :: correctionvalue
  
  INTEGER               :: id, particletype_mpi, unitbuffersize
  INTEGER               :: ibmin,ibmax, jbmin, jbmax, kbmin, kbmax, DIM

  TYPE(bounds_info),DIMENSION(3) :: subrange

  TYPE (particleListinfo),   POINTER   :: particlesubList

  LOGICAL :: Ffluidflag
  LOGICAL :: old_2D
 

  CONTAINS


    SUBROUTINE compute_cell_centers ()
  
      INTEGER :: i
    
      ALLOCATE(x_cc_lp(-buff_size:m+buff_size))
      ALLOCATE(y_cc_lp(-buff_size:n+buff_size))
      IF (p.GT.0) ALLOCATE(z_cc_lp(-buff_size:p+buff_size))
    
      DO i=-buff_size,m+buff_size
        x_cc_lp(i) = (x_cb(i-1) + x_cb(i))*0.5d0
      ENDDO
      
      DO i=-buff_size,n+buff_size
        y_cc_lp(i) = (y_cb(i-1) + y_cb(i))*0.5d0
      ENDDO
      
      IF (p.GT.0) THEN
        DO i=-buff_size,p+buff_size
          z_cc_lp(i) = (z_cb(i-1) + z_cb(i))*0.5d0
        ENDDO
      ENDIF
  
    END SUBROUTINE compute_cell_centers
  
 
    FUNCTION get_cell_from_s (s)
  
      INTEGER, DIMENSION(3)          :: get_cell_from_s
      INTEGER                        :: i
      REAL(KIND(0.D0)), DIMENSION(3) :: s
    
      get_cell_from_s(:) = INT(s(:))
    
      DO i=1,DIM
        IF (s(i).LT.0.0d0) get_cell_from_s(i) = get_cell_from_s(i) - 1
      ENDDO
    
    END FUNCTION get_cell_from_s
 

    SUBROUTINE initialize_particlelists
  
      INTEGER :: i,j,k
    
      ALLOCATE( qbl%fp(-buff_size:m+buff_size,-buff_size:n+buff_size,-buff_size:p+buff_size)) 
      DO k=-buff_size,p+buff_size
        DO j=-buff_size,n+buff_size
          DO i=-buff_size,m+buff_size
            qbl%fp(i,j,k)%nb = 0
            ALLOCATE ( qbl%fp(i,j,k)%List )
            NULLIFY (qbl%fp(i,j,k)%List%next) !initialize lists
            NULLIFY (qbl%fp(i,j,k)%List%prev)
            NULLIFY (qbl%fp(i,j,k)%List%data)
          ENDDO
        ENDDO
      ENDDO
    
      CALL init_cell_list(cellwbList)
  
    END SUBROUTINE initialize_particlelists
  
   
    SUBROUTINE Remove_node_from_qbl (particle,nodeaux,cell)
    
      TYPE(particlenode),POINTER :: particle, nodeaux
      INTEGER, DIMENSION(3)      :: cell
      TYPE (cellwb), POINTER     :: cellwbaux
    
      801 IF (Associated(nodeaux%data,particle%data)) THEN
        CALL Remove_from_3Dlist(nodeaux,cell)  
        IF (qbl%fp(cell(1),cell(2),cell(3))%nb.EQ.0) THEN 
          cellwbaux => qbl%fp(cell(1),cell(2),cell(3))%cellpointer 
          CALL Remove_from_celllist(cellwbaux)
        ENDIF
      ELSE
        nodeaux => nodeaux%next
        GOTO 801
      ENDIF
    
    END SUBROUTINE Remove_node_from_qbl
    
   
    SUBROUTINE update_sublists(particle,oldcell)
    
      TYPE(particleListinfo), POINTER :: particleListaux 
      TYPE(particlenode),POINTER      :: particle, nodeaux
      TYPE (cellwb), POINTER        :: cellwbaux
      INTEGER,DIMENSION(3)          :: oldcell,cell
      INTEGER                       :: i
    
      cell = get_cell_from_s(particle%data%tmp%s)
      !Removing particle from the old cell
      nodeaux => qbl%fp(oldcell(1),oldcell(2),oldcell(3))%List%next
      CALL Remove_node_from_qbl (particle,nodeaux,cell)
    
      !Adding the particle to the new cell
      particleListaux => qbl%fp(cell(1),cell(2),cell(3))
      IF (particleListaux%nb.EQ.0) CALL addtocell_list ( cell )
      CALL addparticletolist (particle%data,particleListaux)
    
    END SUBROUTINE update_sublists
       
   
    SUBROUTINE addparticletolist(particle, particleList)
   
      !Adding particle to the head of the list
      TYPE(particledata),POINTER      :: particle
      TYPE(particlenode),POINTER      :: node
      TYPE(particleListinfo),POINTER  :: particleList
      
      ALLOCATE(node)
      particleList%nb   = particleList%nb + 1
      node%data => particle
      node%next => particleList%List%next
      nullify(node%prev)
      IF (ASSOCIATED(particleList%List%next)) particleList%List%next%prev => node
      particleList%List%next => node
    
    END SUBROUTINE addparticletolist
    
   
    SUBROUTINE Remove_particle (particle, originlist)
    
      TYPE(particlenode),POINTER :: sublistnode,particle, nodeaux
      TYPE (cellwb), POINTER     :: cellwbaux
      INTEGER, DIMENSION(3)      :: cell
      INTEGER                    :: originlist !1=from qbl, 2=SubList
    
      IF (originlist.EQ.1) THEN
        sublistnode => particlesubList%List%next
        DO WHILE (.NOT.Associated(sublistnode%data,particle%data))
          sublistnode => sublistnode%next
        ENDDO
      ELSE
        sublistnode => particle
      ENDIF
    
      !Removing element from the 3D list structure  
      cell = get_cell_from_s(particle%data%tmp%s)
      nodeaux => qbl%fp(cell(1),cell(2),cell(3))%List%next
      CALL Remove_node_from_qbl (particle,nodeaux,cell)
    
      !deallocating particle sublist node (if it is the first)
      !deallocating data
      IF (Associated(sublistnode%next)) THEN
          IF (Associated(sublistnode%prev)) THEN
            sublistnode%prev%next => sublistnode%next
            sublistnode%next%prev => sublistnode%prev
          ELSE
            nullify(sublistnode%next%prev)
          ENDIF
      ELSE
          IF (Associated(sublistnode%prev)) nullify(sublistnode%prev%next) !required?
      ENDIF
      
      IF (Associated(sublistnode%data,particlesubList%List%next%data)) THEN !first element 
        IF (Associated(sublistnode%next)) THEN
            particlesubList%List%next =>  sublistnode%next
        ELSE
            nullify(particlesubList%List%next)
        ENDIF  
      ENDIF
      particlesubList%nb = particlesubList%nb - 1
      
      !next particle
      IF (originlist.EQ.1) THEN
        particle => qbl%fp(cell(1),cell(2),cell(3))%List%next
      ELSE
        IF (Associated(particle%next)) THEN
          particle => particle%next
        ELSE
          nullify(particle)
        ENDIF
      ENDIF
    
      DEALLOCATE(sublistnode%data)  
      DEALLOCATE(sublistnode)  
    
    END SUBROUTINE Remove_particle
    
   
    SUBROUTINE Remove_from_3Dlist(nodeaux,oldcell)
    
      TYPE(particlenode),POINTER    :: oldnode,nodeaux
      INTEGER,DIMENSION(3)          :: oldcell
    
      oldnode   => nodeaux
      IF (Associated(nodeaux%next)) THEN
          IF (Associated(nodeaux%prev)) nodeaux%prev%next => nodeaux%next
          nodeaux%next%prev => nodeaux%prev
      ELSE
          IF (Associated(nodeaux%prev)) nullify(nodeaux%prev%next) !last element of the list
      ENDIF

      IF (Associated(nodeaux,qbl%fp(oldcell(1),oldcell(2),oldcell(3))%List%next)) THEN !first element
          IF (Associated(nodeaux%next)) THEN
              qbl%fp(oldcell(1),oldcell(2),oldcell(3))%List%next => nodeaux%next
          ELSE
              nullify(qbl%fp(oldcell(1),oldcell(2),oldcell(3))%List%next)
          ENDIF
      ENDIF
      DEALLOCATE(oldnode)
      qbl%fp(oldcell(1),oldcell(2),oldcell(3))%nb =   qbl%fp(oldcell(1),oldcell(2),oldcell(3))%nb - 1
 
    
    END SUBROUTINE Remove_from_3Dlist
    
   
    SUBROUTINE init_cell_list(cellList)
    
      TYPE(cellListinfo), POINTER :: cellList
  
      ALLOCATE(cellList)
      ALLOCATE(cellList%List)
      nullify (cellList%List%next)
      nullify (cellList%List%prev)
      nullify (cellList%List%data)
      cellList%nb   = 0
     
    END SUBROUTINE init_cell_list
    
   
    SUBROUTINE Remove_from_celllist(cellwbaux)
    
      TYPE (cellwb), POINTER        :: cellwbaux, oldcellnode
  
      oldcellnode  => cellwbaux
      IF (Associated(cellwbaux%next)) THEN
        IF (Associated(cellwbaux%prev)) cellwbaux%prev%next => cellwbaux%next
          cellwbaux%next%prev => cellwbaux%prev
        ELSE
          IF (Associated(cellwbaux%prev))  nullify(cellwbaux%prev%next) !last element of the list
      ENDIF
  
      IF (Associated(cellwbaux,cellwbList%List%next)) THEN !first element 
        IF (Associated(cellwbaux%next)) THEN
          cellwbList%List%next => cellwbaux%next
          ELSE
            nullify(cellwbList%List%next)
        ENDIF
      ENDIF
  
      cellwbList%nb = cellwbList%nb - 1
      DEALLOCATE(oldcellnode%data); DEALLOCATE(oldcellnode)
   
    END SUBROUTINE Remove_from_celllist
    
   
    SUBROUTINE init_particle_list(particleList)
    
      TYPE(particleListinfo), POINTER :: particleList
  
      ALLOCATE(particleList%List)
      nullify(particleList%next)
      nullify(particleList%prev)
      nullify (particleList%List%next)
      nullify (particleList%List%prev)
      nullify (particleList%List%data)
      particleList%nb = 0
  
    END SUBROUTINE init_particle_list
    

    SUBROUTINE addtocell_list( cell )

      !Adding bubble to the head of the list   
      TYPE(cellwb),POINTER          :: cellnode
      TYPE(cellwbcoord), POINTER    :: cellinfo
      INTEGER, DIMENSION(3)         :: cell     !cell 
  
      ALLOCATE(cellnode);  ALLOCATE(cellinfo)
      cellwbList%nb = cellwbList%nb + 1
      cellinfo%coord=cell
      cellnode%data => cellinfo
      cellnode%next => cellwbList%List%next
      nullify(cellnode%prev)
      qbl%fp(cell(1),cell(2),cell(3))%cellpointer => cellnode
      IF (ASSOCIATED(cellwbList%List%next)) cellwbList%List%next%prev => cellnode
      cellwbList%List%next => cellnode
    
    END SUBROUTINE addtocell_list
    
   
    SUBROUTINE Create_sublists ()
     
      ! It will be interesting to put this subroutine in structparticle.f90 and to
      ! pass the lists as arguments.
  
      TYPE (dirlist),POINTER    :: mdirList,dirnode
      TYPE (cellwb), POINTER    :: cellwbaux
      TYPE (particlenode),POINTER   :: node
      INTEGER,DIMENSION(3)      :: cell
      INTEGER                   :: i
  
      ALLOCATE(particlesubList)
      CALL init_particle_list (particlesubList)
      cellwbaux => cellwbList%List%next
   
      DO WHILE(Associated(cellwbaux))
        node  => qbl%fp(cellwbaux%data%coord(1), cellwbaux%data%coord(2), cellwbaux%data%coord(3))%List%next
        DO WHILE(Associated(node))
          CALL addparticletolist (node%data,particlesubList)
          node => node%next
        ENDDO
        cellwbaux => cellwbaux%next
      ENDDO
    
    END SUBROUTINE Create_sublists
    
   
    FUNCTION Interpolate (coord, q)
      
      !TYPE( field_position ),  INTENT(IN)         :: q
      TYPE(scalar_field), INTENT(IN)              :: q
      REAL(KIND(0.D0)), DIMENSION(3), INTENT(IN)  :: coord
      REAL(KIND(0.D0))                            :: Interpolate
  
      IF (p.eq.0) THEN
        Interpolate =  Interpolate2D (coord, q)
      ELSE
        Interpolate =  Interpolate3D (coord, q)
      ENDIF
    
    END FUNCTION Interpolate
 

    SUBROUTINE get_psi (coord,psi,cell)
      ! this subroutine returns the local psi coordinates 
      ! it also gives the cell with respect these coordinates are defined
    
      REAL(KIND(0.D0)), DIMENSION(3), INTENT(IN)  :: coord
      REAL(KIND(0.D0)), DIMENSION(3)              :: psi !local coordinates
      INTEGER, DIMENSION(3)                       :: cell
  
      cell = get_cell_from_s(coord)
  
      !obtain psi(1)
      psi(1) = (coord(1) - REAL(cell(1)))*dx(cell(1)) + x_cb(cell(1)-1)
      IF (cell(1).EQ.(m+buff_size)) THEN
        cell(1) = cell(1) - 1
        psi(1)    = 1.0d0
      ELSEIF (cell(1).EQ.(-buff_size)) THEN
        psi(1) = 0.0d0
      ELSE
        IF (psi(1).LT.x_cc_lp(cell(1))) cell(1) = cell(1) - 1
        psi(1) = ABS((psi(1) - x_cc_lp(cell(1)))/(x_cc_lp(cell(1)+1) - x_cc_lp(cell(1))))
      ENDIF
  
      !obtain psi(2)
      psi(2) = (coord(2) - REAL(cell(2)))*dy(cell(2)) + y_cb(cell(2)-1)
      IF (cell(2).EQ.(n+buff_size)) THEN
        cell(2) = cell(2) - 1
        psi(2)    = 1.0d0
      ELSEIF (cell(2).EQ.(-buff_size)) THEN
        psi(2) = 0.0d0
      ELSE
        IF (psi(2).LT.y_cc_lp(cell(2))) cell(2) = cell(2) - 1
        psi(2) = ABS((psi(2) - y_cc_lp(cell(2)))/(y_cc_lp(cell(2)+1) - y_cc_lp(cell(2))))
      ENDIF
  
      !obtain psi(3)
      IF (p.GT.0) THEN
        psi(3) = (coord(3) - REAL(cell(3)))*dz(cell(3)) + z_cb(cell(3)-1)
        IF (cell(3).EQ.(p+buff_size)) THEN
          cell(3) = cell(3) - 1
          psi(3)    = 1.0d0
          ELSEIF (cell(3).EQ.(-buff_size)) THEN
          psi(3) = 0.0d0
        ELSE
        IF (psi(3).LT.z_cc_lp(cell(3))) cell(3) = cell(3) - 1
          psi(3) = ABS((psi(3) - z_cc_lp(cell(3)))/(z_cc_lp(cell(3)+1) - z_cc_lp(cell(3))))
        ENDIF
      ELSE
        psi(3) = 0.0d0
      ENDIF
  
    END SUBROUTINE get_psi
 

    FUNCTION Interpolate2D (coord, q)

      !  !interpolation in the computational space
      !  fixme: This interpolation is not strictly correct for cell centered values
      !  Bilinear interpolation?
  
      !TYPE( field_position ),  INTENT(IN)         :: q
      TYPE(scalar_field), INTENT(IN)              :: q
      REAL(KIND(0.D0)), DIMENSION(3), INTENT(IN)  :: coord
      REAL(KIND(0.D0)),DIMENSION(3)               :: psi !local coordinates
      REAL(KIND(0.D0))                            :: Interpolate2D, tmp
      INTEGER, DIMENSION(3)                       :: cell
      INTEGER :: i
  
      CALL get_psi (coord,psi,cell)
  
      tmp =       q%sf(cell(1),cell(2),cell(3))      *(1.0d0-psi(1))*(1.0d0-psi(2))
      tmp = tmp + q%sf(cell(1)+1,cell(2),cell(3))    *psi(1)*(1.0d0-psi(2))
      tmp = tmp + q%sf(cell(1)+1,cell(2)+1,cell(3))  *psi(1)*psi(2)
      tmp = tmp + q%sf(cell(1),cell(2)+1,cell(3))    *(1.0d0-psi(1))*psi(2)
  
      Interpolate2D = tmp
  
    END FUNCTION Interpolate2D
   
   
    FUNCTION Interpolate3D (coord, q)
      ! Interpolation in the computational space
      !  fixme: This interpolation is not strictly correct for cell centered values
      ! Bilinear interpolation?
  
      !TYPE( field_position ),  INTENT(IN)         :: q
      TYPE(scalar_field), INTENT(IN)              :: q
      REAL(KIND(0.D0)), DIMENSION(3), INTENT(IN)  :: coord
      REAL(KIND(0.D0)), DIMENSION(3)              :: psi !local coordinates
      REAL(KIND(0.D0))                            :: Interpolate3D, tmp
      INTEGER, DIMENSION(3)                       :: cell
      INTEGER :: i
  
      CALL get_psi (coord, psi, cell)
  
      tmp =       q%sf(cell(1),cell(2),cell(3))      *(1.0d0-psi(1))*(1.0d0-psi(2))*(1.0d0-psi(3))
      tmp = tmp + q%sf(cell(1)+1,cell(2),cell(3))    *psi(1)*(1.0d0-psi(2))*(1.0d0-psi(3))
      tmp = tmp + q%sf(cell(1)+1,cell(2)+1,cell(3))  *psi(1)*psi(2)*(1.0d0-psi(3))
      tmp = tmp + q%sf(cell(1),cell(2)+1,cell(3))    *(1.0d0-psi(1))*psi(2)*(1.0d0-psi(3))
      tmp = tmp + q%sf(cell(1),cell(2),cell(3)+1)    *(1.0d0-psi(1))*(1.0d0-psi(2))*psi(3)
      tmp = tmp + q%sf(cell(1)+1,cell(2),cell(3)+1)  *psi(1)*(1.0d0-psi(2))*psi(3)
      tmp = tmp + q%sf(cell(1)+1,cell(2)+1,cell(3)+1)*psi(1)*psi(2)*psi(3)
      tmp = tmp + q%sf(cell(1),cell(2)+1,cell(3)+1)  *(1.0d0-psi(1))*psi(2)*psi(3)
  
      Interpolate3D = tmp
  
    END FUNCTION Interpolate3D
  
  
    SUBROUTINE locate_cell ( pos, cell, scoord )
      ! This subroutine gives back the computational coordinate of the pos into scoord.
      ! Computational coordinate has a unit of "cell".
      ! The numbering of the cell of which left boundary is the domain boundary is 0. 
      ! if comp.coord of the pos is s, the real coordinate of s is
      ! (the coordinate of the left boundary of the Floor(s)-th cell)
      ! + (s-(int(s))*(cell-width).
      ! In other words,  the coordinate of the center of the cell is x_cc_lp(cell).
  
      REAL(KIND(0.D0)), DIMENSION(3)   :: pos
      REAL(KIND(0.D0)), DIMENSION(3), OPTIONAL  :: scoord
      INTEGER, DIMENSION(3)            :: cell
      INTEGER :: i, j, k
    
      DO WHILE (pos(1).LT.x_cb(cell(1)-1))
        cell(1) = cell(1)-1
      ENDDO
  
      DO WHILE (pos(1).GT.x_cb(cell(1))) 
        cell(1) = cell(1)+1
      ENDDO
  
      DO WHILE (pos(2).LT.y_cb(cell(2)-1)) 
        cell(2) = cell(2)-1
      ENDDO
  
      DO WHILE (pos(2).GT.y_cb(cell(2))) 
        cell(2) = cell(2)+1
      ENDDO
   
      IF (p.gt.0) THEN
          DO WHILE (pos(3).LT.z_cb(cell(3)-1)) 
            cell(3) = cell(3)-1
          ENDDO
          DO WHILE (pos(3).GT.z_cb(cell(3))) 
            cell(3) = cell(3)+1
          ENDDO
      ENDIF
    
      !coordinates in computational space
      IF (PRESENT(scoord)) THEN
          scoord(1) = cell(1) + (pos(1) - x_cb(cell(1)-1))/dx(cell(1))
          scoord(2) = cell(2) + (pos(2) - y_cb(cell(2)-1))/dy(cell(2))
          scoord(3) = 0.0d0
          IF (p.GT.0) scoord(3) = cell(3) + (pos(3) - z_cb(cell(3)-1))/dz(cell(3))
          cell(:) = get_cell_from_s (scoord)
      ENDIF     
  
    END SUBROUTINE locate_cell
  
  
    SUBROUTINE get_epsilonbaux ( center, cell, stddsv, epsilonbaux )
  
      REAL(KIND(0.D0)), DIMENSION(3)  :: center,pos
      REAL(KIND(0.D0))                :: stddsv
      INTEGER                         :: i
      INTEGER, DIMENSION(3)           :: cell,newcell,epsilonbaux
      LOGICAL                         :: indomain
    
      DO i=1,DIM 
        pos(i)   = center(i) + stddsv
        indomain = particle_in_domain(pos) 
        IF (.NOT.indomain) THEN
          pos(i)   = center(i) - stddsv
          indomain = particle_in_domain(pos) 
        ENDIF
        IF (indomain) THEN
          newcell  = cell
          CALL locate_cell ( pos,  newcell)
          epsilonbaux(i) = 5*max(epsilonb,REAL(abs(newcell(i)-cell(i))))
        ELSE
          epsilonbaux(i) = 5*epsilonb
        ENDIF
      ENDDO
      
    END SUBROUTINE get_epsilonbaux
 

    SUBROUTINE transfertodata (particletmp)
    
      TYPE(particledata) :: particletmp
    
      particletmp%x = particletmp%tmp%x 
      particletmp%u = particletmp%tmp%u  
      particletmp%y = particletmp%tmp%y 
      particletmp%p = particletmp%tmp%p 
      particletmp%mv = particletmp%tmp%mv 
    
    END SUBROUTINE
    
 
    SUBROUTINE transfertotmp (particletmp)
  
      TYPE(particledata) :: particletmp
  
      particletmp%tmp%x  = particletmp%x
      particletmp%tmp%u  = particletmp%u
      particletmp%tmp%y  = particletmp%y
      particletmp%tmp%p  = particletmp%p
      particletmp%tmp%mv = particletmp%mv
  
    END SUBROUTINE
   
   
    SUBROUTINE clear_particle_list(particleList)
    
       TYPE(particleListinfo), POINTER :: particleList
       TYPE(particlenode),POINTER      :: particle, particleold
   
       particleList%nb = 0
   
       particle  => particleList%List%next
   
       DO WHILE(Associated(particle))
          particleold => particle
          particle    => particle%next
          DEALLOCATE(particleold)
          NULLIFY(particleold)
      ENDDO
   
      NULLIFY(particleList%List%next)

    END SUBROUTINE clear_particle_list
  
  
    SUBROUTINE get_char_dist ( cell, Chardist)
    
      REAL(KIND(0.D0))      :: Chardist
      INTEGER, DIMENSION(3) :: cell
      
      IF (p>0) THEN
         Chardist = (dx(cell(1))*dy(cell(2))*dz(cell(3)))**(1./3.)
      ELSE
         Chardist = sqrt(dx(cell(1))*dy(cell(2)))
      ENDIF
    
    END SUBROUTINE get_char_dist
    
    
    SUBROUTINE get_char_vol ( cell, Charvol)
    
      REAL(KIND(0.D0))      :: Charvol
      INTEGER, DIMENSION(3) :: cell
      
      IF (p>0) THEN
         Charvol = dx(cell(1))*dy(cell(2))*dz(cell(3))
      ELSE
        IF(cyl_coord) THEN
           Charvol = dx(cell(1))*dy(cell(2))*y_cc_lp(cell(2))*2d0*PI
        ELSE IF(old_2D .NEQV. .TRUE.) THEN
           Charvol = dx(cell(1))*dy(cell(2))*charwidth
        ELSE
           Charvol = dx(cell(1))*dy(cell(2))*charwidth
        END IF
      ENDIF
    
    END SUBROUTINE get_char_vol
    
    !-----------------------------------------------------------------------
    
    FUNCTION Modulus ( array )
    
      REAL(KIND(0.D0)),DIMENSION(3), INTENT(IN) :: array
      REAL(KIND(0.D0)) ::  modulus
    
      Modulus = sqrt(array(1)**2 + array(2)**2 + array(3)**2)
    
    END FUNCTION Modulus
    
    !-----------------------------------------------------------------------
    
    FUNCTION particle_in_domain (pos_part) 
    
      LOGICAL                        :: particle_in_domain
      REAL(KIND(0.D0)), DIMENSION(3) :: pos_part
    
      ! 2D
      IF(p.eq.0 .AND. cyl_coord.NEQV..TRUE.) THEN
        ! For the old 2D kernel (see Fuster and Colonius, JFM, 2011))
        IF(old_2D) THEN
           particle_in_domain = ((pos_part(1).LT.x_cb(m+buff_size)).AND.(pos_part(1).GE.x_cb(-buff_size-1)).AND. &
                                 (pos_part(2).LT.y_cb(n+buff_size)).AND.(pos_part(2).GE.y_cb(-buff_size-1)))
        ! For present 2D kernel
        ELSE
           ! Defining a virtual z-axis that has the same dimensions as y-axis
           ! defined in the input file
           particle_in_domain = ((pos_part(1).LT.x_cb(m+buff_size)).AND.(pos_part(1).GE.x_cb(-buff_size-1)).AND. &
                                 (pos_part(2).LT.y_cb(n+buff_size)).AND.(pos_part(2).GE.y_cb(-buff_size-1)).AND. &
                                 !(pos_part(3).LT.y_cb(n+buff_size)).AND.(pos_part(3).GE.y_cb(-buff_size-1)))
                                  (pos_part(3).LT.charwidth/2d0).AND.(pos_part(3).GE.-charwidth/2d0))
        END IF
      ELSE
        ! cyl_coord
        particle_in_domain = ((pos_part(1).LT.x_cb(m+buff_size)).AND.(pos_part(1).GE.x_cb(-buff_size-1)).AND. &
                              (ABS(pos_part(2)).LT.y_cb(n+buff_size)).AND.(ABS(pos_part(2)).GE.MAX(y_cb(-buff_size-1),0d0)))
      ENDIF   
   
      ! 3D
      IF (p > 0) THEN
        particle_in_domain = ((pos_part(1).LT.x_cb(m+buff_size)).AND.(pos_part(1).GE.x_cb(-buff_size-1)).AND. &
                             (pos_part(2).LT.y_cb(n+buff_size)).AND.(pos_part(2).GE.y_cb(-buff_size-1)).AND. &
                             (pos_part(3).LT.z_cb(p+buff_size)).AND.(pos_part(3).GE.z_cb(-buff_size-1)))
      ENDIF
    
      ! For symmetric boundary condition
      IF(bc_x%beg.eq.-2) THEN
        particle_in_domain = (particle_in_domain.AND.(pos_part(1).GE.x_cb(-1)))
      END IF
      IF(bc_x%end.eq.-2) THEN
        particle_in_domain = (particle_in_domain.AND.(pos_part(1).LT.x_cb(m)))
      END IF
      IF(bc_y%beg.eq.-2.AND.(.NOT.cyl_coord)) THEN
        particle_in_domain = (particle_in_domain.AND.(pos_part(2).GE.y_cb(-1)))
      END IF
      IF(bc_y%end.eq.-2.AND.(.NOT.cyl_coord)) THEN
        particle_in_domain = (particle_in_domain.AND.(pos_part(2).LT.y_cb(n)))
      END IF
    
      IF (p > 0) THEN
         IF(bc_z%beg.eq.-2) THEN
           particle_in_domain = (particle_in_domain.AND.(pos_part(3).GE.z_cb(-1)))
         END IF
         IF(bc_z%end.eq.-2) THEN
           particle_in_domain = (particle_in_domain.AND.(pos_part(3).LT.z_cb(p)))
         END IF
      END IF
    
    END FUNCTION particle_in_domain  
    
    
    FUNCTION particle_in_domain_physical (pos_part) 
      ! This subroutine is used for mpi_parallel_io
      LOGICAL                        :: particle_in_domain_physical
      REAL(KIND(0.D0)), DIMENSION(3) :: pos_part
    
      particle_in_domain_physical = ((pos_part(1).LT.x_cb(m)).AND.(pos_part(1).GE.x_cb(-1)).AND. &
                            (pos_part(2).LT.y_cb(n)).AND.(pos_part(2).GE.y_cb(-1)))
      !particle_in_domain = ((pos_part(1).LT.x_cb(m+buff_size)).AND.(pos_part(1).GE.x_cb(-buff_size-1)).AND. &
      !                      (pos_part(2).LT.y_cb(n+buff_size)).AND.(pos_part(2).GE.y_cb(-buff_size-1))).AND. &
      !                      (pos_part(3).LT.charwidth/2d0.AND.pos_part(3).GE.-charwidth/2d0)
    
    
      IF (p > 0) THEN
        particle_in_domain_physical = (particle_in_domain_physical.AND.(pos_part(3).LT.z_cb(p)).AND.(pos_part(3).GE.z_cb(-1)))
      ENDIF
    
    END FUNCTION particle_in_domain_physical  
    
    
    !-----------------------------------------------------------------------
    
    SUBROUTINE gradient_dir(q,dq,i)
        
      !TYPE( field_position )   :: q,dq
      TYPE(scalar_field)       :: q, dq
      INTEGER                  :: i
    
      IF (i.EQ.1) THEN
        CALL gradientx(q,dq)
      ELSE
        IF (i.EQ.2) THEN
          CALL gradienty(q,dq)
        ELSE
          CALL gradientz(q,dq)
        ENDIF
      ENDIF
    
    END SUBROUTINE gradient_dir
    
    SUBROUTINE gradientx(q,dq)
      
        !TYPE( field_position )   :: q,dq
        TYPE(scalar_field)       :: q, dq
        INTEGER                  :: i,j,k,l,lmax
        REAL(KIND(0.D0))         :: aux1,aux2 
       
        !DO k=subrange(3)%beg,subrange(3)%end
    !    DO j=subrange(2)%beg,subrange(2)%end
    !    DO i=max(subrange(1)%beg,-buffer_size+1),min(subrange(1)%end,m+buffer_size-1)
        DO k=0, p
        DO j=0, n
        DO i=0, m
            aux1= dx(i)+dx(i-1)
            aux2= dx(i)+dx(i+1)
            dq%sf(i,j,k) = q%sf(i,j,k)*(dx(i+1)-dx(i-1)) &
                         + q%sf(i+1,j,k)*aux1 &
                         - q%sf(i-1,j,k)*aux2
            dq%sf(i,j,k) = dq%sf(i,j,k)/(aux1*aux2)
        ENDDO
        ENDDO
        ENDDO
    
        !IF (subrange(1)%beg.EQ.-buffer_size) THEN
    !    DO k=-buffer_size, p+buffer_size
    !    DO j=-buffer_size, n+buffer_size
    !     !   DO j=subrange(2)%beg,subrange(2)%end
    !     !   DO k=subrange(3)%beg,subrange(3)%end
    !        dq%fp(-buffer_size,j,k) = (q%fp(-buffer_size+1,j,k) - q%fp(-buffer_size,j,k))/dx(-buffer_size)
    !        ENDDO
    !        ENDDO
    !    !ENDIF
    !    !IF (subrange(1)%end.EQ.m+buffer_size) THEN
    !    DO k=-buffer_size, p+buffer_size
    !    DO j=-buffer_size, n+buffer_size
    !     !   DO j=subrange(2)%beg,subrange(2)%end
    !     !   DO k=subrange(3)%beg,subrange(3)%end
    !        dq%fp(m+buffer_size,j,k) = (q%fp(m+buffer_size,j,k) - q%fp(m+buffer_size-1,j,k))/dx(m+buffer_size)
    !        ENDDO
    !        ENDDO
    !   ! ENDIF
    !        
    
    END SUBROUTINE gradientx
    
    !-----------------------------------------------------------------------
    
    SUBROUTINE gradienty(q,dq)
      
        !TYPE( field_position )   :: q,dq
        TYPE(scalar_field)       :: q, dq
        INTEGER                  :: i,j,k,l,lmax
        REAL(KIND(0.D0))         :: aux1,aux2 
        
        !DO k=subrange(3)%beg,subrange(3)%end
        !DO j=max(subrange(2)%beg,-buffer_size+1),min(subrange(2)%end,n+buffer_size-1)
        !DO i=subrange(1)%beg,subrange(1)%end
        DO k=0, p
        DO j=0, n
        DO i=0, m
            aux1= dy(j)+dy(j-1)
            aux2= dy(j)+dy(j+1)
            dq%sf(i,j,k) = q%sf(i,j,k)*(dy(j+1)-dy(j-1)) &
                         + q%sf(i,j+1,k)*aux1  &
                         - q%sf(i,j-1,k)*aux2  
            dq%sf(i,j,k) = dq%sf(i,j,k)/(aux1*aux2)
        ENDDO
        ENDDO
        ENDDO
    
    !    !IF (subrange(2)%beg.EQ.-buffer_size) THEN
    !        DO i=-buffer_size, m+buffer_size
    !        DO k=-buffer_size, p+buffer_size
    !        !DO i=subrange(1)%beg,subrange(1)%end
    !        !DO k=subrange(3)%beg,subrange(3)%end
    !          dq%fp(i,-buffer_size,k) = (q%fp(i,-buffer_size+1,k) - q%fp(i,-buffer_size,k))/dy(-buffer_size)
    !        ENDDO
    !        ENDDO
    !!    ENDIF
    !!    IF (subrange(2)%end.EQ.n+buffer_size) THEN
    !        DO i=-buffer_size, m+buffer_size
    !        DO k=-buffer_size, p+buffer_size
    !        !DO i=subrange(1)%beg,subrange(1)%end
    !        !DO k=subrange(3)%beg,subrange(3)%end
    !          dq%fp(i,n+buffer_size,k) = (q%fp(i,n+buffer_size,k) - q%fp(i,n+buffer_size-1,k))/dy(n+buffer_size)
    !        ENDDO
    !        ENDDO
    !!    ENDIF
    !        
    
    END SUBROUTINE gradienty
    
    !-----------------------------------------------------------------------
    
    SUBROUTINE gradientz(q,dq)
      
        !TYPE( field_position )   :: q,dq
        TYPE(scalar_field)       :: q, dq
        INTEGER                  :: i,j,k,l,lmax
        REAL(KIND(0.D0))         :: aux1,aux2 
        
        DO k=0, p
        DO j=0, n
        DO i=0, m
        !DO k=max(subrange(3)%beg,-buffer_size+1),min(subrange(3)%end,p+buffer_size-1)
        !DO j=subrange(2)%beg,subrange(2)%end
        !DO i=subrange(1)%beg,subrange(1)%end
            aux1= dz(k)+dz(k-1)
            aux2= dz(k)+dz(k+1)
            dq%sf(i,j,k) = q%sf(i,j,k)*(dz(k+1)-dz(k-1)) &
                         + q%sf(i,j,k+1)*aux1 &
                         - q%sf(i,j,k-1)*aux2 
            dq%sf(i,j,k) = dq%sf(i,j,k)/(aux1*aux2)
        ENDDO
        ENDDO
        ENDDO
    
    !    !IF (subrange(3)%beg.EQ.-buffer_size) THEN
    !    !    DO j=subrange(2)%beg,subrange(2)%end
    !    !    DO i=subrange(1)%beg,subrange(1)%end
    !    DO j=-buffer_size, n+buffer_size
    !    DO i=-buffer_size, m+buffer_size
    !          dq%fp(i,j,-buffer_size) = (q%fp(i,j,-buffer_size+1) - q%fp(i,j,-buffer_size))/dz(-buffer_size)
    !        ENDDO
    !        ENDDO
    !    !ENDIF
    !    !IF (subrange(3)%end.EQ.p+buffer_size) THEN
    !    DO j=-buffer_size, n+buffer_size
    !    DO i=-buffer_size, m+buffer_size
    !     !   DO j=subrange(2)%beg,subrange(2)%end
    !      !  DO i=subrange(1)%beg,subrange(1)%end
    !          dq%fp(i,j,p+buffer_size) = (q%fp(i,j,p+buffer_size) - q%fp(i,j,p+buffer_size-1))/dz(p+buffer_size)
    !        ENDDO
    !        ENDDO
    !!    ENDIF
    !        
    
    END SUBROUTINE gradientz
    
    !-----------------------------------------------------------------------
    
    SUBROUTINE get_equilibrium_radius (ntrial, x,  tolx, tolf, particle, pliq, largestep)
    
      ! uses lubksb, ludcmp usrfun
      !       Given an initial guess x for a root in n dimensions, take ntrial
      !       Newton-Rhapson steps to improve the root. Step if the root 
      !       converges in either summed absolute variable increments tolx or
      !       summed absolute function values tolf 
      
      INTEGER, PARAMETER  :: np=15,n=1
      INTEGER             :: ntrial,i,k,indx(np)
      TYPE(particledata)  :: particle
      REAL(KIND(0.D0))    :: tolf, tolx, x(n), pliq 
      REAL(KIND(0.D0))    :: d, errf,errx,fjac(np,np),fvec(np),p(np)
      LOGICAL             :: largestep
      
      DO k=1, ntrial
        ! function to supply the values of the function at x in fvec and the
        ! Jacobian Matrix at fjac  
        call eq_rad_jac(x,fvec,fjac,particle, pliq) 
        errf=0.
        DO i=1,n
          errf= errf+abs(fvec(i))
        ENDDO
        IF (errf.le.tolf) RETURN
        DO i=1,n
          p(i)=-fvec(i)
        ENDDO
        call ludcmp(fjac,indx,d,largestep)
        IF (largestep) RETURN
        call lubksb(fjac,indx,p)
        errx=0.
        DO i=1,n
          errx=errx+abs(p(i))
          x(i)=x(i)+p(i)
        ENDDO
        IF (errx.le.tolx) RETURN
      ENDDO

    END SUBROUTINE get_equilibrium_radius
    
    ! From Fuster's code
    ! functions to solve a nonlinear system of equations 
    ! In this case I only solve for the radius in the Laplace equation.
    SUBROUTINE ludcmp (a, indx, d,cond)
    
      INTEGER            :: indx(2)
      INTEGER, PARAMETER :: nmax=500,n=1, np=15
      REAL(KIND(0.D0)), PARAMETER :: tiny = 1.0d-20
      REAL(KIND(0.D0))   :: d,a(15,15)
      INTEGER            :: i,imax, j,k
      REAL(KIND(0.D0))   :: aamax , dum, sum, vv(nmax)
      LOGICAL            :: cond
      
      d=1.
      DO i=1,n
        aamax=0.
        DO j=1,n
          IF (abs(a(i,j)).GT.aamax) aamax=abs(a(i,j))
        ENDDO
        IF (aamax.eq.0.) THEN
          cond =.TRUE.
          RETURN
        ENDIF
        vv(i)=1./aamax
      ENDDO
      
      DO j=1,n
        DO i=1, j-1
          sum=a(i,j)
          DO k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          ENDDO
          a(i,j)=sum
        ENDDO
        aamax=0.
        DO i=j,n
          sum=a(i,j)
          DO k=1, j-1
            sum=sum-a(i,k)*a(k,j)
          ENDDO
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          IF (dum.ge.aamax) THEN
            imax=i
            aamax=dum
          ENDIF
        ENDDO
        IF (j.ne.imax) THEN
          DO k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          ENDDO
          d=-d
          vv(imax)=vv(j)
        ENDIF
        indx(j)=imax
        IF (a(j,j).EQ.0.) a(j,j)=tiny
        IF (j.ne.n) THEN
          dum = 1./a(j,j)
          DO i=j+1,n
            a(i,j)= a(i,j)*dum
          ENDDO
        ENDIF
      ENDDO
    
    END SUBROUTINE ludcmp


    SUBROUTINE lubksb(a,indx,b)
    
      INTEGER, PARAMETER :: np=15,n=1
      INTEGER            :: indx(n)
      REAL(KIND(0.D0))   :: a(np,np), b(n)
      INTEGER            :: i,ii,j,ll
      REAL(KIND(0.D0))   :: sum
      
      ii=0
      
      DO i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        IF (ii.ne.0) THEN
          DO j=ii,i-1
            sum=sum-a(i,j)*b(j)
          ENDDO
        ELSEIF (sum.ne.0.) THEN
          ii=i
        ENDIF
        b(i)=sum
      ENDDO
      
      DO i=n,1,-1
        sum=b(i)
        DO j=i+1,n
          sum=sum-a(i,j)*b(j)
        ENDDO
        b(i)=sum/a(i,i)
      ENDDO
    
    END SUBROUTINE lubksb


    SUBROUTINE eq_rad_jac (x,fvec,fjac,particle, pliq)
 
      ! Values are dimensionless such that T=1
      TYPE(particledata) :: particle
      INTEGER, PARAMETER :: np=15,n=1
      REAL(KIND(0.D0))   :: x(n), pliq
      REAL(KIND(0.D0))   :: fjac(np,np),fvec(np)
      
        fvec(1) = 4.0d0*pi*x(1)*(pliq-pvap) - 3.0d0*particle%mg*Runiv/MWgas/x(1)**2 + 8.0d0*pi*sigmabubble
        fjac(1,1) = 4.0d0*pi*(pliq-pvap) + 6.0d0*particle%mg*Runiv/MWgas/x(1)**3
      
    END SUBROUTINE eq_rad_jac


END MODULE m_particles_types
