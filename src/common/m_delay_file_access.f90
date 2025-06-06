module m_delay_file_access
    use m_precision_select
    implicit none
    private

    public :: DelayFileAccess

    integer, private, parameter :: &
        N_PROCESSES_FILE_ACCESS = 128, &
        FILE_ACCESS_DELAY_UNIT = 10000

contains

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
