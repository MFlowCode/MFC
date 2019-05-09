! Author: Spencer Bryngelson

module m_compile_specific

    implicit none

    contains

        SUBROUTINE my_inquire(fileloc,dircheck)
            character(len=*), intent (in) :: fileloc
            logical, intent (inout) :: dircheck

            INQUIRE(FILE=TRIM(fileloc),EXIST=dircheck)        !GCC
        END SUBROUTINE my_inquire

end module m_compile_specific
