! pMFC v3.0 - Post-Process Code: m_particles_types.f90
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
     REAL(KIND(0.D0)), DIMENSION(3)   :: x, s, u ! x: real eoord, s: comp coord, u: vel of the particle
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

  REAL(KIND(0.D0)) :: Tini, Runiv, gammagas, gammavapor, pvap, cpgas
  REAL(KIND(0.D0)) :: cpvapor, kgas, kvapor, MWgas, MWvap, diffcoefvap
  REAL(KIND(0.D0)) :: sigmabubble, viscref
  REAL(KIND(0.D0)) :: csonref, rholiqref, Lref
  REAL(KIND(0.D0)) :: correctionvalue, charwidth, valmaxvoid, dtmaxpart
  REAL(KIND(0.D0)) :: RKeps
  REAL(KIND(0.D0)), DIMENSION(3) :: gravity

  INTEGER               :: ratiodt
  INTEGER               :: id, particletype_mpi, unitbuffersize
  INTEGER               :: ibmin,ibmax, jbmin, jbmax, kbmin, kbmax, DIM
  INTEGER               :: clusterflag, massflag, heatflag
  INTEGER               :: smoothtype, projectiontype
  REAL(KIND(0.D0)) :: epsilonb

  TYPE(bounds_info),DIMENSION(3) :: subrange

  TYPE (particleListinfo),   POINTER   :: particlesubList
        
  LOGICAL :: particleoutFlag         ! Flag to activate bubble radius output
  LOGICAL :: RPflag, stillparticlesflag, Ffluidflag
  LOGICAL :: particlestatFlag,correctpresFlag
 
 !---------------------------------------------------------------------
 

END MODULE m_particles_types
