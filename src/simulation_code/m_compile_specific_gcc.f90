!>
!! @file m_compile_specific_gcc.f90
!! @brief Contains module m_compile_specific
!! @author spencer
!! @version 1.1
!! @date 1/1/1

!> @brief This module contains subroutines that are compiler specific (GNU/INTEL)
MODULE m_compile_specific

    IMPLICIT NONE

    CONTAINS

        !>  Inquires on the existence of a directory
        !!  @param fileloc File directory location
        !!  @param dircheck Switch that indicates if directory exists
        SUBROUTINE my_inquire(fileloc,dircheck)
            CHARACTER(LEN=*), INTENT(IN) :: fileloc
            LOGICAL, INTENT(INOUT) :: dircheck

            INQUIRE(FILE=TRIM(fileloc),EXIST=dircheck)        !GCC
        END SUBROUTINE my_inquire

END MODULE m_compile_specific
