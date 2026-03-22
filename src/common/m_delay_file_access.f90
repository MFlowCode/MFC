!>
!! @file
!! @brief Contains module m_delay_file_access

!> @brief Rank-staggered file access delays to prevent I/O contention on parallel file systems
module m_delay_file_access
    use m_precision_select
    implicit none
    private

    public :: DelayFileAccess

    integer, private, parameter :: &
        N_PROCESSES_FILE_ACCESS = 128, &
        FILE_ACCESS_DELAY_UNIT = 10000

contains

    !> @brief Introduces a rank-dependent busy-wait delay to stagger parallel file access and reduce I/O contention.
    impure subroutine DelayFileAccess(ProcessRank)
        integer, intent(in) :: ProcessRank

        integer :: iDelay, nFileAccessDelayIterations
        real(wp) :: Number, Dummy

        nFileAccessDelayIterations &
            = (ProcessRank/N_PROCESSES_FILE_ACCESS)*FILE_ACCESS_DELAY_UNIT

        do iDelay = 1, nFileAccessDelayIterations
            ! Wait my turn
            call random_number(Number)
            Dummy = Number*Number
        end do

    end subroutine DelayFileAccess

end module m_delay_file_access
