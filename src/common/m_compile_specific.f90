!>
!! @file m_compile_specific.f90
!! @brief Contains module m_compile_specific

!> @brief This module contains subroutines that are compiler specific
module m_compile_specific

    ! Dependencies
    use m_mpi_proxy

    implicit none

contains

    !>  Creates a directory and all its parents if it does not exist
        !!  @param dir_name Directory path
    impure subroutine s_create_directory(dir_name)
        character(LEN=*), intent(in) :: dir_name

#ifdef _WIN32
        call system('mkdir "'//dir_name//'" 2> NUL')
#else
        call system('mkdir -p "'//dir_name//'"')
#endif

    end subroutine s_create_directory

    impure subroutine s_delete_file(filepath)
        character(LEN=*), intent(in) :: filepath

#ifdef _WIN32
        call system('del "'//filepath//'"')
#else
        call system('rm "'//filepath//'"')
#endif

    end subroutine s_delete_file

    impure subroutine s_delete_directory(dir_name)
        character(LEN=*), intent(in) :: dir_name

#ifdef _WIN32
        call system('rmdir "'//dir_name//'" /s /q')
#else
        call system('rm -r "'//dir_name//'"')
#endif

    end subroutine s_delete_directory

    !>  Inquires on the existence of a directory
        !!  @param fileloc File directory location
        !!  @param dircheck Switch that indicates if directory exists
    impure subroutine my_inquire(fileloc, dircheck)
        character(LEN=*), intent(in) :: fileloc
        logical, intent(inout) :: dircheck

#ifdef __INTEL_COMPILER
        inquire (DIRECTORY=trim(fileloc), EXIST=dircheck)   !Intel
#else
        inquire (FILE=trim(fileloc), EXIST=dircheck)        !GCC
#endif

    end subroutine my_inquire

    impure subroutine s_get_cwd(cwd)
        character(LEN=*), intent(out) :: cwd

        call GETCWD(cwd)
    end subroutine s_get_cwd

    impure subroutine s_get_basename(dirpath, basename)
        character(LEN=*), intent(in) :: dirpath
        character(LEN=*), intent(out) :: basename

        integer :: iUnit
        character(len=30) :: tmpfilepath

        write (tmpfilepath, '(A,I0)') 'basename_', proc_rank

#ifdef _WIN32
        call system('for /F %i in ("'//trim(dirpath)//'") do @echo %~ni > '//trim(tmpfilepath))
#else
        call system('basename "'//trim(dirpath)//'" > '//trim(tmpfilepath))
#endif

        open (newunit=iUnit, FILE=trim(tmpfilepath), FORM='formatted', STATUS='old')
        read (iUnit, '(A)') basename
        close (iUnit)

        call s_delete_file(trim(tmpfilepath))

    end subroutine s_get_basename

end module m_compile_specific
