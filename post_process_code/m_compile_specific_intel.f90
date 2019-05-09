! Author: Spencer Bryngelson

module m_compile_specific

    implicit none

    contains

        SUBROUTINE my_inquire(fileloc,dircheck)
            character(len=*), intent (in) :: fileloc
            logical, intent (inout) :: dircheck

            INQUIRE(DIRECTORY=TRIM(fileloc),EXIST=dircheck)   !Intel
        END SUBROUTINE my_inquire

end module m_compile_specific
