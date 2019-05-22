!>
!! @file m_compile_specific_intel.f90
!! @brief This module contains subroutines that are compiler specific (GNU/INTEL)
!! @author spencer
!! @version 1.1
!! @date 1/1/1
module m_compile_specific

    implicit none

    contains

        SUBROUTINE my_inquire(fileloc,dircheck)
            character(len=*), intent (in) :: fileloc
            logical, intent (inout) :: dircheck

            INQUIRE(DIRECTORY=TRIM(fileloc),EXIST=dircheck)   !Intel
        END SUBROUTINE my_inquire

end module m_compile_specific
