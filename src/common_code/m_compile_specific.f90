!>
!! @file m_compile_specific.f90
!! @brief Contains module m_compile_specific
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module contains subroutines that are compiler specific
module m_compile_specific

    implicit none

contains

    !>  Creates a directory and all its parents if it does not exist
        !!  @param dir_name Directory path
    subroutine s_create_directory(dir_name)
        character(LEN=*), intent(IN) :: dir_name
        
#ifdef _WIN32
        call system('mkdir "'//dir_name//'" 2> NUL')
#else
        call system('mkdir -p "'//dir_name//'"')
#endif

    end subroutine s_create_directory

    subroutine s_delete_directory(dir_name)
        character(LEN=*), intent(IN) :: dir_name

#ifdef _WIN32
        call system('rmdir "'//dir_name//'" /s /q')
#else
        call system('rm -rf "'//dir_name//'"')
#endif

    end subroutine s_delete_directory

    !>  Inquires on the existence of a directory
        !!  @param fileloc File directory location
        !!  @param dircheck Switch that indicates if directory exists
    subroutine my_inquire(fileloc, dircheck)
        character(LEN=*), intent(IN) :: fileloc
        logical, intent(INOUT) :: dircheck

#ifdef __INTEL_COMPILER
        inquire (DIRECTORY=trim(fileloc), EXIST=dircheck)   !Intel
#else
        inquire (FILE=trim(fileloc), EXIST=dircheck)        !GCC
#endif

    end subroutine my_inquire

end module m_compile_specific
