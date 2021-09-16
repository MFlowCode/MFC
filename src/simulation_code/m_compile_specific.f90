!>
!! @file m_compile_specific.f90
!! @brief Contains module m_compile_specific
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module contains subroutines that are compiler specific
MODULE m_compile_specific

    IMPLICIT NONE

    CONTAINS

        !>  Inquires on the existence of a directory
        !!  @param fileloc File directory location
        !!  @param dircheck Switch that indicates if directory exists
        SUBROUTINE my_inquire(fileloc,dircheck)
            CHARACTER(LEN=*), INTENT(IN) :: fileloc
            LOGICAL, INTENT(INOUT) :: dircheck

#ifdef __INTEL_COMPILER
    INQUIRE(DIRECTORY=TRIM(fileloc),EXIST=dircheck)   !Intel
#else
    INQUIRE(FILE=TRIM(fileloc),EXIST=dircheck)        !GCC
#endif

        END SUBROUTINE my_inquire

END MODULE m_compile_specific
