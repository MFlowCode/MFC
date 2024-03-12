! pMFC_v3.0 - Simulation code: m_kernel_functions.f90
! Description: The module contains kernel functions
! Author: Kazuki Maeda
! Date: 01/01/17

MODULE m_kernel_functions
  
  ! Dependencies =============================================================
  USE m_particles_types
  USE m_mpi_common
  ! ==========================================================================
  
  IMPLICIT NONE

  CONTAINS


  SUBROUTINE smoothfunction ( updatedvar, center, cell, strength, kernelused, stddsv, strength2 )
  
      !TYPE (field_position)           :: updatedvar
      TYPE(scalar_field)              :: updatedvar
      REAL(KIND(0.D0)), DIMENSION(3)  :: center
      REAL(KIND(0.D0))                :: strength, stddsv
      INTEGER, DIMENSION(3)           :: cell
      INTEGER                         :: kernelused
      REAL(KIND(0.D0)), OPTIONAL      :: strength2
  
      smoothfunc: SELECT CASE( kernelused )
      CASE (1)
          CALL gaussian ( updatedvar, center, cell, strength, stddsv, strength2 )
      CASE (2)
          CALL deltafunc ( updatedvar, cell, strength )
      END SELECT smoothfunc
  
  END SUBROUTINE smoothfunction
  


  SUBROUTINE smoothfunction_az ( beta,center,cell,volpart,kernelused,stddsv,dphidt,rdot,dpdx,dpdy,dudx,dudy)
  
      TYPE(scalar_field)              :: beta,dpdx,dpdy,dudx,dudy
      REAL(KIND(0.D0)), DIMENSION(3)  :: center
      REAL(KIND(0.D0))                :: stddsv,volpart,dphidt,rdot
      INTEGER, DIMENSION(3)           :: cell
      INTEGER                         :: kernelused
      INTEGER                         :: err_code, ierr
  
      smoothfunc: SELECT CASE( kernelused )
      CASE (1)
          CALL gaussian_az ( beta,center,cell,volpart,stddsv,dphidt,rdot,dpdx,dpdy,dudx,dudy )
      CASE (2)
          PRINT '(A)', "discrete kernel is not implemented for azimuthal subgrid modeling"
          !CALL MPI_ABORT(MPI_COMM_WORLD, err_code, ierr)
          call s_mpi_abort()
   
      END SELECT smoothfunc
  
  END SUBROUTINE smoothfunction_az


  
  SUBROUTINE gaussian ( updatedvar, center, cell, strength, stddsv, strength2)
  
      REAL(KIND(0.D0)), DIMENSION(3)  :: center
      REAL(KIND(0.D0))                :: strength, stddsv
      INTEGER                         :: i
      INTEGER, DIMENSION(3)           :: cell,epsilonbaux
      TYPE(scalar_field)              :: updatedvar
      REAL(KIND(0.D0)), OPTIONAL      :: strength2
     
      !fixme: epsilonbaux(:) must be strictly specified in terms of stddsv?
      !       however stddsv must be smaller than buffer?
      
      ! For now conservative for restart when epsilonbaux < buff_size
      epsilonbaux(:) = 3
    
      CALL applygaussian ( updatedvar, center, cell, strength, epsilonbaux, stddsv, strength2 )
      !For symmetric BC
      IF((bc_x%beg.eq.-2 .or. bc_x%end.eq.-2 .or. bc_y%beg.eq.-2 .or. bc_y%end.eq.-2 .or. bc_z%beg.eq.-2 .or. bc_z%end.eq.-2).OR.(bc_x%beg.eq.proc_rank .or. bc_x%end.eq.proc_rank .or. bc_y%beg.eq.proc_rank .or. bc_y%end.eq.proc_rank .or. bc_z%beg.eq.proc_rank .or. bc_z%end.eq.proc_rank)) THEN
          CALL gaussian_symmetric_bc (updatedvar, center, cell, strength, epsilonbaux, stddsv, strength2)
      END IF
  
  END SUBROUTINE gaussian
  


  SUBROUTINE gaussian_az (beta,center,cell,volpart,stddsv,dphidt,rdot,dpdx,dpdy,dudx,dudy)
  
      REAL(KIND(0.D0)), DIMENSION(3)  :: center
      REAL(KIND(0.D0))                :: strength,stddsv,volpart,dphidt,rdot
      INTEGER                         :: i
      INTEGER, DIMENSION(3)           :: cell,epsilonbaux
      TYPE(scalar_field)              :: beta, dpdx,dpdy,dudx,dudy
   
      ! For now conservative for restart when epsilonbaux < buff_size
      ! see SUBROUTINE gaussian
      epsilonbaux(:) = 3
     
      CALL applygaussian_az ( beta, center, cell, volpart, epsilonbaux, stddsv, dphidt,rdot,dpdx,dpdy,dudx,dudy )
  
      !For symmetric BC
      IF(bc_x%beg.eq.-2 .or. bc_x%end.eq.-2) THEN
          CALL gaussian_symmetric_bc_az (beta, center, cell, volpart, epsilonbaux, stddsv,dphidt,rdot,dpdx,dpdy,dudx,dudy)
      END IF
  
  END SUBROUTINE gaussian_az
  


  SUBROUTINE gaussian_symmetric_bc (updatedvar, center, cell, strength, epsilonbaux, stddsv, strength2)
  
      ! This subroutine smeares the strength of fictitious particles, across the
      ! boundary when symmetric boundary condition is used
      REAL(KIND(0.D0)), DIMENSION(3)  :: center,centertmp
      REAL(KIND(0.D0))                :: strength, stddsv
      INTEGER                         :: i
      INTEGER, DIMENSION(3)           :: cell,celltmp,epsilonbaux
      TYPE(scalar_field)              :: updatedvar
      REAL(KIND(0.D0)), OPTIONAL      :: strength2
   
      celltmp   = cell
      centertmp = center
  
      IF (((bc_x%beg.EQ.-2).OR.(bc_x%beg.EQ.proc_rank)) .AND. cell(1)+1.LE.epsilonbaux(1)) THEN
          
          CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,1)
          CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
          !Apply(x-d,y,z)
  
          IF (((bc_y%beg.EQ.-2).OR.(bc_y%beg.EQ.proc_rank)) .AND. cell(2)+1.LE.epsilonbaux(2) .AND.(.NOT.cyl_coord)) THEN
  
              CALL kernel_shift_bc(center,centertmp,cell,celltmp,2,1)
              CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
              !Apply(x-d,y-d,z)
  
              centertmp(1)= center(1)
              celltmp(1)= cell(1)
              CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
              !Apply(x,y-d,z)
  
              IF (p > 0) THEN
                  IF (((bc_z%beg.EQ.-2).OR.(bc_z%beg.EQ.proc_rank)) .AND. cell(3)+1.LE.epsilonbaux(3)) THEN
                  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y-d,z-d)
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y-d,z-d)
   
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y,z-d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z-d)
  
                  ELSEIF (((bc_z%end.EQ.-2).OR.(bc_z%end.EQ.proc_rank)) .AND. cell(3).GE.p+1-epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y-d,z+d)
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y-d,z+d)
   
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y,z+d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z+d)
  
                 ENDIF
  
              ENDIF
  
          ELSEIF (((bc_y%end.EQ.-2).OR.(bc_y%end.EQ.proc_rank)) .AND. cell(2).GE.n+1-epsilonbaux(2) .AND.(.NOT.cyl_coord)) THEN
  
              CALL kernel_shift_bc(center,centertmp,cell,celltmp,2,2)
              CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
              !Apply(x-d,y+d,z)
  
              centertmp(1)= center(1)
              celltmp(1)= cell(1)
              CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
              !Apply(x,y+d,z)
  
              IF (p > 0) THEN
                  IF (((bc_z%beg.EQ.-2).OR.(bc_z%beg.EQ.proc_rank)) .AND. cell(3)+1.LE.epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y+d,z-d)
                      
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y+d,z-d)
   
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y,z-d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z-d)
  
                    
                  ELSEIF (((bc_z%end.EQ.-2).OR.(bc_z%end.EQ.proc_rank)) .AND. cell(3).GE.p+1-epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y+d,z+d)
                      
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y+d,z+d)
   
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y,z+d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z+d)
  
  
                  ENDIF
              ENDIF
          ELSE
              IF (p > 0) THEN
                  IF (((bc_z%beg.EQ.-2).OR.(bc_z%beg.EQ.proc_rank)) .AND. cell(3)+1.LE.epsilonbaux(3)) THEN
   
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y,z-d)
   
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z-d)
   
                  ELSEIF (((bc_z%end.EQ.-2).OR.(bc_z%end.EQ.proc_rank)) .AND. cell(3)+1.GE.p+1-epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x-d,y,z+d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z+d)
   
                  ENDIF
              ENDIF
          ENDIF
      ELSEIF (((bc_x%end.EQ.-2).OR.(bc_x%end.EQ.proc_rank)) .AND. cell(1).GE.m+1-epsilonbaux(1)) THEN
  
          CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,2)
          CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
          !Apply(x+d,y,z)
  
          IF (((bc_y%beg.EQ.-2).OR.(bc_y%beg.EQ.proc_rank)) .AND. cell(2)+1.LE.epsilonbaux(2).AND.(.NOT.cyl_coord)) THEN
  
              CALL kernel_shift_bc(center,centertmp,cell,celltmp,2,1)
              CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
              !Apply(x+d,y-d,z)
  
              centertmp(1)= center(1)
              celltmp(1)= cell(1)
              CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
              !Apply(x,y-d,z)
  
              IF (p > 0) THEN
                  IF (((bc_z%beg.EQ.-2).OR.(bc_z%beg.EQ.proc_rank)) .AND. cell(3)+1.LE.epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y-d,z-d)
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y-d,z-d)
   
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y,z-d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z-d)
  
                  ELSEIF (((bc_z%end.EQ.-2).OR.(bc_z%end.eq.proc_rank)) .AND. cell(3).GE.p+1-epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y-d,z+d)
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y-d,z+d)
   
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y,z+d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z+d)
  
                  ENDIF
              ENDIF
          ELSEIF (((bc_y%end.EQ.-2).OR.(bc_y%end.EQ.proc_rank)) .AND. cell(2).GE.n+1-epsilonbaux(2).AND.(.NOT.cyl_coord)) THEN
  
                  CALL kernel_shift_bc(center,centertmp,cell,celltmp,2,2)
                  CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                  !Apply(x+d,y+d,z)
  
                  centertmp(1)= center(1)
                  celltmp(1)= cell(1)
                  CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                  !Apply(x,y+d,z)
  
              IF (p > 0) THEN
                  IF (((bc_z%beg.EQ.-2).OR.(bc_z%beg.EQ.proc_rank)) .AND. cell(3)+1.LE.epsilonbaux(3)) THEN
  
                      call kernel_shift_bc(center,centertmp,cell,celltmp,3,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y+d,z-d)
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y+d,z-d)
   
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y,z-d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z-d)
  
                  ELSEIF (((bc_z%end.EQ.-2).OR.(bc_z%end.EQ.proc_rank)) .AND. cell(3)+1.GE.p+1-epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y+d,z+d)
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,1,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y+d,z+d)
   
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y,z+d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z+d)
  
                  ENDIF
              ENDIF
          ELSE
              IF (p > 0) THEN
                  IF (((bc_z%beg.EQ.-2).OR.(bc_z%beg.EQ.proc_rank)) .AND. cell(3)+1.LE.epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y,z-d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z-d)
  
                  ELSEIF (((bc_z%end.EQ.-2).OR.(bc_z%end.EQ.proc_rank)) .AND. cell(3).GE.p+1-epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x+d,y,z+d)
  
                      centertmp(1)= center(1)
                      celltmp(1)= cell(1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z+d)
  
                  ENDIF
              ENDIF
          ENDIF
      ELSE
          IF (((bc_y%beg.EQ.-2).OR.(bc_y%beg.EQ.proc_rank)) .AND. cell(2)+1.LE.epsilonbaux(2).AND.(.NOT.cyl_coord)) THEN
  
              CALL kernel_shift_bc(center,centertmp,cell,celltmp,2,1)
              CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
              !Apply(x,y-d,z)
  
              IF(p > 0) THEN
                  IF (((bc_z%beg.EQ.-2).OR.(bc_z%beg.EQ.proc_rank)) .AND. cell(3)+1.LE.epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y-d,z-d)
  
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z-d)
  
                  ELSEIF (((bc_z%end.EQ.-2).OR.(bc_z%end.EQ.proc_rank)) .AND. cell(3).GE.p+1-epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y-d,z+d)
  
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z+d)
  
                  ENDIF
              ENDIF
          ELSEIF (((bc_y%end.eq.-2).OR.(bc_y%end.EQ.proc_rank)) .AND. cell(2).GE.n+1-epsilonbaux(2).AND.(.NOT.cyl_coord)) THEN
                  CALL kernel_shift_bc(center,centertmp,cell,celltmp,2,2)
                  CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                  !Apply(x,y+d,z)
  
              IF(p > 0) THEN
                  IF (((bc_z%beg.EQ.-2).OR.(bc_z%beg.EQ.proc_rank)) .AND. cell(3)+1.LE.epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y+d,z-d)
  
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z-d)
  
                  ELSEIF (((bc_z%end.EQ.-2).OR.(bc_z%end.EQ.proc_rank)) .AND. cell(3).GE.p+1-epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y+d,z+d)
  
                      centertmp(2)= center(2)
                      celltmp(2)= cell(2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z+d)
  
                  ENDIF
              ENDIF
          ELSE
              IF(p > 0) THEN
                  IF (((bc_z%beg.EQ.-2).OR.(bc_z%beg.EQ.proc_rank)) .AND. cell(3)+1.LE.epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,1)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z-d)
  
                  ELSEIF (((bc_z%end.EQ.-2).OR.(bc_z%end.EQ.proc_rank)) .AND. cell(3).GE.p+1-epsilonbaux(3)) THEN
  
                      CALL kernel_shift_bc(center,centertmp,cell,celltmp,3,2)
                      CALL applygaussian(updatedvar,centertmp,celltmp,strength,epsilonbaux,stddsv, strength2)
                      !Apply(x,y,z+d)
  
                  ENDIF
              ENDIF
          ENDIF
      ENDIF
  
  
  END SUBROUTINE gaussian_symmetric_bc


  
  SUBROUTINE gaussian_symmetric_bc_az (beta, center, cell, volpart, epsilonbaux, stddsv,dphidt,rdot,dpdx,dpdy,dudx,dudy)
  
      ! This subroutine smeares the strength of fictitious particles, across the
      ! boundary when symmetric boundary condition is used
      REAL(KIND(0.D0)), DIMENSION(3)  :: center,centertmp
      REAL(KIND(0.D0))                :: volpart, stddsv,dphidt,rdot
      INTEGER                         :: i
      INTEGER, DIMENSION(3)           :: cell,celltmp,epsilonbaux
      TYPE(scalar_field)              :: beta,dpdx,dpdy,dudx,dudy
  
      celltmp   = cell
      centertmp = centertmp
  
      ! For 2D axi-symmetric, we can consider only x boundaries.
      IF (bc_x%beg.eq.-2 .AND. cell(1)+1.LE.epsilonbaux(1)) THEN
          centertmp(1)= 2*x_cb(-1)-center(1)
          celltmp(1)= -1-cell(1)
          CALL applygaussian_az( beta, centertmp, celltmp, volpart, epsilonbaux, stddsv, dphidt,rdot,dpdx,dpdy,dudx,dudy )
      ELSEIF (bc_x%end.eq.-2 .AND. cell(1).GE.m+1-epsilonbaux(1)) THEN
          centertmp(1)= 2*x_cb(m)-center(1)
          celltmp(1)= 2*m+1-cell(1)
          CALL applygaussian_az( beta, centertmp, celltmp, volpart, epsilonbaux, stddsv, dphidt,rdot,dpdx,dpdy,dudx,dudy )
      ENDIF
  
  
  END SUBROUTINE gaussian_symmetric_bc_az


  
  SUBROUTINE applygaussian ( updatedvar, center, cell, strength, epsilonbaux, stddsv, strength2)
  
  
      REAL(KIND(0.D0)), DIMENSION(3)    :: nodecoord, center, auxvect
      REAL(KIND(0.D0))                  :: strength, func, distance, stddsv
      INTEGER, DIMENSION(3)             :: cell, cellaux, cellaux2
      INTEGER                           :: idir, idir2, i, j, k, ini,iend, jini,jend
      INTEGER, DIMENSION(3)             :: epsilonbaux
      TYPE(scalar_field)                :: updatedvar
      LOGICAL                           :: celloutside
      REAL(KIND(0.D0))                  :: theta, dtheta, L2, dz, Lz2
      INTEGER                           :: Nr, Nr_count
      REAL(KIND(0.D0)), OPTIONAL        :: strength2 ! for product of two smeared functions
  
      nodecoord = 0
  
      IF (DIM.EQ.3) THEN
          k = -epsilonbaux(3)
      ELSE
          k = 0
      ENDIF
      
       i=-epsilonbaux(1);j=-epsilonbaux(2)
  
  3001 IF ((i.LE.epsilonbaux(1)).AND.(j.LE.epsilonbaux(2))) THEN
  
           celloutside = .FALSE.
  
           cellaux(1) = cell(1) + i
           cellaux(2) = cell(2) + j
           cellaux(3) = cell(3) + k
  
              IF (cellaux(1).LT.-buff_size) THEN
                  celloutside = .TRUE.
                  i = i+1
              ENDIF
  
              !check ghost part in y direction
              IF (cellaux(2).LT.-buff_size) THEN
                  celloutside = .TRUE.
                  j = j+1
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
                  ENDIF
              END IF
  
              IF (celloutside) GOTO 3001
              ! Temp
              IF (cellaux(3).GT.p+buff_size) RETURN
              IF (cellaux(2).GT.n+buff_size) GOTO 3002
              IF (cellaux(1).GT.m+buff_size) GOTO 3003 !RETURN
  
              !z direction, periodic           
              nodecoord(1) = x_cc_lp(cellaux(1)) 
              nodecoord(2) = y_cc_lp(cellaux(2)) 
              IF (p > 0) nodecoord(3) = z_cc_lp(cellaux(3)) 
              auxvect(:) = center(:)-nodecoord(:)
              distance =  Modulus(auxvect)
              IF (DIM.EQ.3) THEN
                  func = exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3
              ELSE
                  
                  ! for 2D cylindrical coordinate we smear particles in the azimuthal
                  ! direction for given r
                  IF(cyl_coord) THEN
                      theta = 0d0
                      Nr = CEILING(2d0*PI*nodecoord(2)/(y_cb(cellaux(2))-y_cb(cellaux(2)-1)))
                      dtheta= 2d0*PI/Nr
                      L2 = center(2)**2d0+nodecoord(2)**2d0-2d0*center(2)*nodecoord(2)*COS(theta)
                      distance = DSQRT(auxvect(1)**2d0+L2)
                      ! Factor 2d0 is for symmetry (upper half of the
                      ! 2D field (+r) is considered)
                      func = dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3
                      Nr_count = 0
                      DO WHILE (Nr_count.LT.Nr-1)
                          Nr_count = Nr_count + 1
                          theta = Nr_count*dtheta
                          ! trigonometric relation
                          L2 = center(2)**2d0+nodecoord(2)**2d0-2d0*center(2)*nodecoord(2)*COS(theta)
                          distance = DSQRT(auxvect(1)**2d0+L2)
                          ! nodecoord(2)*dtheta is the azimuthal width of the cell
                          IF(PRESENT(strength2)) THEN
                             ! Product of two Gaussians
                             func = func + &
                             dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**6
                          ELSE
                             func = func + &
                             dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3
                          END IF
                      ENDDO
  
                  ! 2D with virtual depth
                  ELSE IF (old_2D.NEQV..TRUE.) THEN
                  
                      theta = 0d0
                      Nr = CEILING(charwidth/(y_cb(cellaux(2))-y_cb(cellaux(2)-1)))
                      Nr_count = 1-epsilonbaux(3)
                      dz = y_cb(cellaux(2)+1)-y_cb(cellaux(2))
                      Lz2 = (center(3)-(dz*(5d-1+Nr_count)-charwidth/2d0))**2d0
                      distance = DSQRT(auxvect(1)**2d0+auxvect(2)**2d0+Lz2)
                      func = dz/charwidth*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3
                      DO WHILE (Nr_count.LT.Nr-1+(epsilonbaux(3)-1))
                          Nr_count = Nr_count + 1
                          Lz2 = (center(3)-(dz*(5d-1+Nr_count)-charwidth/2d0))**2d0
                          distance = DSQRT(auxvect(1)**2d0+auxvect(2)**2d0+Lz2)
                          IF(PRESENT(strength2)) THEN
                             ! Product of two Gaussians
                             func = func + &
                             dz/charwidth*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**6
                          ELSE
                             func = func + &
                             dz/charwidth*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3
                          END IF
                      ENDDO
  
                  ELSE
  
                      func = exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3    
                      DO WHILE (distance.LT.5.0d0*stddsv)
                          auxvect(3) = auxvect(3) + charwidth
                          distance =  Modulus(auxvect)
                          func = func + 2.0d0*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3    
                      ENDDO
                  END IF
              END IF
  
              IF(PRESENT(strength2)) THEN
                  updatedvar%sf(cellaux(1),cellaux(2),cellaux(3)) = updatedvar%sf(cellaux(1),cellaux(2),cellaux(3)) + func*strength*strength2
              ELSE
                  updatedvar%sf(cellaux(1),cellaux(2),cellaux(3)) = updatedvar%sf(cellaux(1),cellaux(2),cellaux(3)) + func*strength
              END IF
  
           IF (j.LT.epsilonbaux(2)) THEN
               j = j+1
               GOTO 3001
           ENDIF
  
  3002     j=-epsilonbaux(2)
           i = i+1
           GOTO 3001
  
       ENDIF
  
  3003 IF ((DIM.EQ.3).AND.(k.LT.epsilonbaux(3))) THEN
           k = k+1
           i=-epsilonbaux(1);j=-epsilonbaux(2)
           GOTO 3001
       ENDIF
  
  END SUBROUTINE applygaussian


  
  SUBROUTINE applygaussian_az( beta, center, cell, volpart, epsilonbaux, stddsv, dphidt,rdot,dpdx,dpdy,dudx,dudy )
  
      REAL(KIND(0.D0)), DIMENSION(3)    :: nodecoord, center, auxvect
      REAL(KIND(0.D0)), DIMENSION(3)    :: node_rp,node_rm,node_zp,node_zm
      REAL(KIND(0.D0)), DIMENSION(3)    :: auxv_rp,auxv_rm,auxv_zp,auxv_zm
      REAL(KIND(0.D0))                  :: func, distance, stddsv,volpart,radpart
      REAL(KIND(0.D0))                  :: dist_rp,dist_rm,dist_zp,dist_zm
      REAL(KIND(0.D0))                  :: L2_rp,L2_rm,L2_zp,L2_zm
      REAL(KIND(0.D0))                  :: fac_d1_x,fac_d1_y,fac_d4_x,fac_d4_y
      REAL(KIND(0.D0))                  :: func_d1_x_beta,func_d1_y_beta,func_d4_x_beta,func_d4_y_beta
      INTEGER, DIMENSION(3)             :: cell, cellaux, cellaux2
      INTEGER                           :: idir, idir2, i, j, k, ini,iend, jini,jend
      INTEGER, DIMENSION(3)             :: epsilonbaux
      TYPE(scalar_field)                :: beta,dpdx,dpdy,dudx,dudy
      LOGICAL                           :: celloutside
      REAL(KIND(0.D0))                  :: theta, dtheta, L2
      INTEGER                           :: Nr, Nr_count
      REAL(KIND(0.D0))                  :: dphidt,rdot
  
      nodecoord = 0
      node_zp = 0
      node_zm = 0
      node_rp = 0
      node_rm = 0
  
      ! Getting the bubble radius from the volume
      radpart=(3.0d0/4.0d0/pi*volpart)**(1.0d0/3.0d0)
  
      k = 0
     
      i=-epsilonbaux(1);j=-epsilonbaux(2)
  
  3001 IF ((i.LE.epsilonbaux(1)).AND.(j.LE.epsilonbaux(2))) THEN
  
           celloutside = .FALSE.
  
           cellaux(1) = cell(1) + i
           cellaux(2) = cell(2) + j
  
           IF (cellaux(1).LT.-buff_size) THEN
               celloutside = .TRUE.
               i = i+1
           ENDIF
   
           !check ghost part in y direction
           IF (cellaux(2).LT.-buff_size) THEN
               celloutside = .TRUE.
               j = j+1
           ENDIF
   
           IF(cyl_coord) THEN
               IF (y_cc_lp(cellaux(2)).LT.0d0) THEN
                   celloutside = .TRUE.
                   j = j+1
               ENDIF
           END IF
   
           ! Temp
           IF (celloutside) GOTO 3001
           ! Temp
           IF (cellaux(2).GT.n+buff_size) GOTO 3002
           IF (cellaux(1).GT.m+buff_size) RETURN
   
           ! Coordinate of the node center
           nodecoord(1) = x_cc_lp(cellaux(1)) 
           nodecoord(2) = y_cc_lp(cellaux(2))
           auxvect(:) = center(:)-nodecoord(:)
           distance =  Modulus(auxvect)
   
           ! Coordinate of the node surface center (x+1/2,y)
           node_zp(1) = x_cb(cellaux(1))
           node_zp(2) = y_cc_lp(cellaux(2))
           auxv_zp(:) = center(:)-node_zp(:)
           dist_zp    = Modulus(auxv_zp)
   
           ! Coordinate of the node surface center (x-1/2,y)
           node_zm(1) = x_cb(cellaux(1)-1)
           node_zm(2) = y_cc_lp(cellaux(2))
           auxv_zm(:) = center(:)-node_zm(:)
           dist_zm    = Modulus(auxv_zm)
   
           ! Coordinate of the node surface center (x,y+1/2)
           node_rp(1) = x_cc_lp(cellaux(1))
           node_rp(2) = y_cb(cellaux(2))
           auxv_rp(:) = center(:)-node_rp(:)
           dist_rp    = Modulus(auxv_rp)
   
           ! Coordinate of the node surface center (x,y-1/2)
           node_rm(1) = x_cc_lp(cellaux(1))
           node_rm(2) = y_cb(cellaux(2)-1)
           auxv_rm(:) = center(:)-node_rm(:)
           dist_rm    = Modulus(auxv_rm)
              
           ! For 2D cylindrical coordinate we smear particles in the azimuthal
           ! direction for given r
           theta = 0d0
           Nr = CEILING(2d0*PI*nodecoord(2)/(y_cb(cellaux(2))-y_cb(cellaux(2)-1)))
           dtheta= 2d0*PI/Nr
  
           ! Node center for Volume
           L2 = center(2)**2d0+nodecoord(2)**2d0-2d0*center(2)*nodecoord(2)*COS(theta)
           distance = DSQRT(auxvect(1)**2d0+L2)
  
           ! Cell fases for dpdx and dudx
           L2_zp = center(2)**2d0+node_zp(2)**2d0-2d0*center(2)*node_zp(2)*COS(theta)
           dist_zp = DSQRT(auxv_zp(1)**2d0+L2_zp)
           L2_zm = center(2)**2d0+node_zm(2)**2d0-2d0*center(2)*node_zm(2)*COS(theta)
           dist_zm = DSQRT(auxv_zm(1)**2d0+L2_zm)
  
           ! Cell faces for dpdy and dudy
           L2_rp = center(2)**2d0+node_rp(2)**2d0-2d0*center(2)*node_rp(2)*COS(theta)
           dist_rp = DSQRT(auxv_rp(1)**2d0+L2_rp)
           L2_rm = center(2)**2d0+node_rm(2)**2d0-2d0*center(2)*node_rm(2)*COS(theta)
           dist_rm = DSQRT(auxv_rm(1)**2d0+L2_rm)
  
           ! Convolution kernels for nabla(1/dist) and nabla(1/dist^4) in x and y directions
           fac_d1_x  = (1.0d0/dist_zp-1.0d0/dist_zm)/(x_cb(cellaux(1))-x_cb(cellaux(1)-1))
           fac_d4_x  = ((1.0d0/dist_zp)**4.0d0-(1.0d0/dist_zm)**4.0d0)/(x_cb(cellaux(1))-x_cb(cellaux(1)-1))
           fac_d1_y  = (1.0d0/dist_rp-1.0d0/dist_rm)/(y_cb(cellaux(1))-y_cb(cellaux(1)-1))
           fac_d4_y  = ((1.0d0/dist_rp)**4.0d0-(1.0d0/dist_rm)**4.0d0)/(y_cb(cellaux(1))-y_cb(cellaux(1)-1))
  
           ! Getting the convolution at the cell
           func_d1_x_beta = dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3 &
                            * fac_d1_x
           func_d4_x_beta = dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3 &
                            * fac_d4_x
           func_d1_y_beta = dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3 &
                            * fac_d1_y
           func_d4_y_beta = dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3 &
                            * fac_d4_y
           
           ! Index for azimuthal cells
           Nr_count = 0
   
           DO WHILE (Nr_count.LT.Nr-1)
  
              Nr_count = Nr_count + 1
              theta = Nr_count*dtheta
  
              ! Void fraction
              ! trigonometric relation
              ! Cell faces for dpdx and dudx
              L2_zp = center(2)**2d0+node_zp(2)**2d0-2d0*center(2)*node_zp(2)*COS(theta)
              dist_zp = DSQRT(auxv_zp(1)**2d0+L2_zp)
              L2_zm = center(2)**2d0+node_zm(2)**2d0-2d0*center(2)*node_zm(2)*COS(theta)
              dist_zm = DSQRT(auxv_zm(1)**2d0+L2_zm)
     
              ! Cell faces for dpdy and dudy
              L2_rp = center(2)**2d0+node_rp(2)**2d0-2d0*center(2)*node_rp(2)*COS(theta)
              dist_rp = DSQRT(auxv_rp(1)**2d0+L2_rp)
              L2_rm = center(2)**2d0+node_rm(2)**2d0-2d0*center(2)*node_rm(2)*COS(theta)
              dist_rm = DSQRT(auxv_rm(1)**2d0+L2_rm)
     
              ! Convolution kernels for nabla(1/dist) and nabla(1/dist^4) in x and y directions
              fac_d1_x  = (1.0d0/dist_zp-1.0d0/dist_zm)/(x_cb(cellaux(1))-x_cb(cellaux(1)-1))
              fac_d4_x  = ((1.0d0/dist_zp)**4.0d0-(1.0d0/dist_zm)**4.0d0)/(x_cb(cellaux(1))-x_cb(cellaux(1)-1))
              fac_d1_y  = (1.0d0/dist_rp-1.0d0/dist_rm)/(y_cb(cellaux(1))-y_cb(cellaux(1)-1))
              fac_d4_y  = ((1.0d0/dist_rp)**4.0d0-(1.0d0/dist_rm)**4.0d0)/(y_cb(cellaux(1))-y_cb(cellaux(1)-1))
     
              ! Getting the convolution at the cell
              func_d1_x_beta = func_d1_x_beta &
                               + dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3 &
                               * fac_d1_x
              func_d4_x_beta = func_d4_x_beta &
                               + dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3 &
                               * fac_d4_x
              func_d1_y_beta = func_d1_y_beta &
                               + dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3 &
                               * fac_d1_y
              func_d4_y_beta = func_d4_y_beta &
                               + dtheta/2d0/PI*exp(-0.5d0*(distance/stddsv)**2)/(DSQRT(2.0d0*pi)*stddsv)**3 &
                               * fac_d4_y
  
           ENDDO
   
           dpdx%sf(cellaux(1),cellaux(2),cellaux(3)) = dpdx%sf(cellaux(1),cellaux(2),cellaux(3)) &
                                                       + volpart*(-func_d1_x_beta*radpart*dphidt - func_d4_x_beta*radpart**4.0d0*rdot**2.0d0/2.0d0)
           dpdy%sf(cellaux(1),cellaux(2),cellaux(3)) = dpdy%sf(cellaux(1),cellaux(2),cellaux(3)) &
                                                       + volpart*(-func_d1_y_beta*radpart*dphidt - func_d4_y_beta*radpart**4.0d0*rdot**2.0d0/2.0d0)
           dudx%sf(cellaux(1),cellaux(2),cellaux(3)) = dudx%sf(cellaux(1),cellaux(2),cellaux(3)) &
                                                       - volpart*(func_d1_x_beta*radpart**2.0d0*rdot)
           dudy%sf(cellaux(1),cellaux(2),cellaux(3)) = dudy%sf(cellaux(1),cellaux(2),cellaux(3)) &
                                                       - volpart*(func_d1_y_beta*radpart**2.0d0*rdot)
  
           IF (j.LT.epsilonbaux(2)) THEN
               j = j+1
               GOTO 3001
           ENDIF
  
  3002     j=-epsilonbaux(2)
           i = i+1
           GOTO 3001
  
       ENDIF
  
  END SUBROUTINE applygaussian_az
  
  
  
  SUBROUTINE deltafunc ( updatedvar, cell, strength)
  
      REAL(KIND(0.D0))                :: Dchar, strength,func,distance,Vol
      INTEGER, DIMENSION(3)           :: cell
      INTEGER                         :: idir,i,j
      !TYPE (field_position)           :: updatedvar
      TYPE(scalar_field)              :: updatedvar
   
      IF (DIM.EQ.2) THEN
          Vol    =  dx(cell(1))*dy(cell(2))*charwidth
          IF(cyl_coord) Vol = dx(cell(1))*dy(cell(2))*y_cc_lp(cell(2))*2d0*PI
      ELSE
          Vol    =  dx(cell(1))*dy(cell(2))*dz(cell(3))
      ENDIF
  
      updatedvar%sf(cell(1),cell(2),cell(3)) = updatedvar%sf(cell(1),cell(2),cell(3)) + strength/Vol
  
  END SUBROUTINE deltafunc
  
  !!---------------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE applydelta ( updatedvar, strength, cell, psi ) 
    !Discrete delta function?    
      
      IMPLICIT NONE
      
      !TYPE (field_position)           :: updatedvar
      TYPE(scalar_field)              :: updatedvar
      REAL(KIND(0.D0)),DIMENSION(2,2,2)   :: stencil !local coordinates
      REAL(KIND(0.D0)),DIMENSION(3)   :: psi !local coordinates
      REAL(KIND(0.D0))                :: strength
      INTEGER, DIMENSION(3)           :: cell, cellaux
      INTEGER                         :: idir,i,j,k
  
      stencil(1,1,1)  = (1.0d0-psi(1))*(1.0d0-psi(2))*(1.0d0-psi(3))
      stencil(2,1,1)  = psi(1)*(1.0d0-psi(2))*(1.0d0-psi(3))
      stencil(1,2,1)  = (1.0d0-psi(1))*psi(2)*(1.0d0-psi(3))
      stencil(2,2,1)  = psi(1)*psi(2)*(1.0d0-psi(3))
  
      IF ( p > 0 ) THEN
        stencil(1,1,2)  = (1.0d0-psi(1))*(1.0d0-psi(2))*psi(3)
        stencil(2,1,2)  = psi(1)*(1.0d0-psi(2))*psi(3)
        stencil(1,2,2)  = (1.0d0-psi(1))*psi(2)*psi(3)
        stencil(2,2,2)  = psi(1)*psi(2)*psi(3)
        DO i=0,1; DO j=0,1; DO k=0,1
          cellaux(1)   = min(max(cell(1)+i,-buff_size),m+buff_size)
          cellaux(2)   = min(max(cell(2)+j,-buff_size),n+buff_size)
          cellaux(3)   = min(max(cell(3)+k,-buff_size),p+buff_size)
          updatedvar%sf(cellaux(1),cellaux(2),cellaux(3)) = updatedvar%sf(cellaux(1),cellaux(2),cellaux(3)) + stencil(i+1,j+1,k+1)*strength
        ENDDO;ENDDO;ENDDO
      ELSE
        DO i=0,1; DO j=0,1
          cellaux(1)   = min(max(cell(1)+i,-buff_size),m+buff_size)
          cellaux(2)   = min(max(cell(2)+j,-buff_size),n+buff_size)
          updatedvar%sf(cellaux(1),cellaux(2),0) = updatedvar%sf(cellaux(1),cellaux(2),0) + stencil(i+1,j+1,1)*strength
        ENDDO;ENDDO
      ENDIF
  
  END SUBROUTINE applydelta
  !
  !!------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE remeshdelta ( updatedvar, node, strength, rangecells )
  !!Check in parallel!!!!!!!    
      
      REAL(KIND(0.D0))                :: strength,func,distance
      REAL(KIND(0.D0)),DIMENSION(3)   :: psi !local coordinates
      INTEGER, DIMENSION(3)           :: cell
      TYPE(particledata)              :: node
      !TYPE (field_position)           :: updatedvar
      TYPE(scalar_field)              :: updatedvar
      INTEGER                         :: i
      INTEGER, DIMENSION(3,2),OPTIONAL :: rangecells
  
      CALL get_psi (node%tmp%s,psi,cell)
      CALL applydelta ( updatedvar, strength, cell, psi ) 
      
      IF (PRESENT(rangecells)) CALL update_rangecells (cell, rangecells)
  
  END SUBROUTINE remeshdelta
  
  !!------------------------------------------------------------------------------------------------------
  
  SUBROUTINE update_rangecells (cell, rangecells)
  
    INTEGER, DIMENSION(3,2),OPTIONAL :: rangecells
    INTEGER, DIMENSION(3)            :: cell
    INTEGER                          :: i
  
    DO i=1,3
      rangecells(i,1) = max(min(rangecells(i,1),cell(i)-1),-buff_size)
      rangecells(i,2) = max(rangecells(i,2),cell(i)+1)
    ENDDO
    rangecells(1,2) =  min(rangecells(1,2),m+buff_size)
    rangecells(2,2) =  min(rangecells(2,2),n+buff_size)
    IF ( p.EQ.0 ) THEN
      rangecells(3,:) = 0
    ELSE 
      rangecells(3,2) = min(rangecells(3,2),p+buff_size)
    ENDIF
  
  END SUBROUTINE update_rangecells
  


  SUBROUTINE kernel_shift_bc(center,centertmp,cell,celltmp,i,j)
   
      REAL(KIND(0.D0)), DIMENSION(3)  :: center, centertmp
      REAL(KIND(0.D0))                :: strength, stddsv
      INTEGER, DIMENSION(3)           :: cell, celltmp
      ! Coordinate index(i=1,2,3 -> x,y,z) and boundary index (j=1,2 -> beg,end)
      INTEGER                         :: i, j
  
      IF(i == 1) THEN
  
          IF(j == 1) THEN
              IF(bc_x%beg.eq.proc_rank) THEN
                  centertmp(1) = center(1) + x_cb(m) - x_cb(-1)
                  celltmp(1)   = cell(1) + m+1
              ELSE
                  centertmp(1) = 2*x_cb(-1)-center(1)
                  celltmp(1)   = -1-cell(1)
              END IF
          ELSE
              IF(bc_x%end.eq.proc_rank) THEN
                  centertmp(1) = center(1) - x_cb(m) + x_cb(-1)
                  celltmp(1)   = cell(1) - m-1
              ELSE
                  centertmp(1) = 2*x_cb(m)-center(1)
                  celltmp(1)   = 2*m+1-cell(1)
              END IF
          END IF
      ELSE IF(i == 2) THEN
          IF(j == 1) THEN
              IF(bc_y%beg.eq.proc_rank) THEN
                  centertmp(2) = center(2) + y_cb(n) - y_cb(-1)
                  celltmp(2)   = cell(2) + n+1
              ELSE
                  centertmp(2) = 2*y_cb(-1)-center(2)
                  celltmp(2)   = -1-cell(2)
              END IF
          ELSE
              IF(bc_y%end.eq.proc_rank) THEN
                  centertmp(2) = center(2) - y_cb(n) + y_cb(-1)
                  celltmp(2)   = cell(2) - n-1
              ELSE
                  !print *, "yend"
                  centertmp(2)  = 2*y_cb(n)-center(2)
                  celltmp(2)   = 2*n+1-cell(2)
              END IF
          END IF
      ELSE IF(i == 3) THEN
          IF(j == 1) THEN
              IF(bc_z%beg.eq.proc_rank) THEN
                  centertmp(3) = center(3) + z_cb(p) - z_cb(-1)
                  celltmp(3)   = cell(3) + p+1
              ELSE
                  centertmp(3) = 2*z_cb(-1)-center(3)
                  celltmp(3)   = -1-cell(3)
              END IF
          ELSE
              IF(bc_z%end.eq.proc_rank) THEN
                  centertmp(3) = center(3) - z_cb(p) + z_cb(-1)
                  celltmp(3)   = cell(3) - p-1
              ELSE
                  centertmp(3) = 2*z_cb(p)-center(3)
                  celltmp(3)   = 2*p+1-cell(3)
              END IF
          END IF
      END IF
  
  END SUBROUTINE kernel_shift_bc


END MODULE m_kernel_functions
