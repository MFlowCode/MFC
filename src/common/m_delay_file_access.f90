!>
!! @file
!! @brief Contains module m_delay_file_access

!> @brief Rank-staggered file access delays to prevent I/O contention on parallel file systems
module m_delay_file_access

    use m_precision_select

    implicit none

    private

    public :: s_delay_file_access

    integer, private, parameter :: N_PROCESSES_FILE_ACCESS = 128, FILE_ACCESS_DELAY_UNIT = 10000

contains

    !> Introduce a rank-dependent busy-wait delay to stagger parallel file access and reduce I/O contention.
    impure subroutine s_delay_file_access(process_rank)

        integer, intent(in) :: process_rank
        integer             :: i, n_file_access_delay_iterations
        real(wp)            :: num, dummy

        n_file_access_delay_iterations = (process_rank/N_PROCESSES_FILE_ACCESS)*FILE_ACCESS_DELAY_UNIT

        do i = 1, n_file_access_delay_iterations
            call random_number(num)
            dummy = num*num
        end do

    end subroutine s_delay_file_access

end module m_delay_file_access
