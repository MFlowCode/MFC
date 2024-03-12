! pMFC_v3.0 - Simulation code: m_time_tmp.f90
! Description: The module contains a subroutine to store simulation time
! variables
! Author: Kazuki Maeda
! Date: 01/01/17

! Tentative treatment of simulation time
MODULE m_time_tmp

    IMPLICIT NONE

    INTEGER, PUBLIC :: time_tmp
    INTEGER, PUBLIC :: tot_step
    REAL(KIND(0d0)), PUBLIC :: time_real
    REAL(KIND(0d0)), PUBLIC :: dt_next_inp
    
END MODULE m_time_tmp

