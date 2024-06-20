module m_delay_file_access
    implicit none
    private

    public :: DelayFileAccess

    integer, private, parameter :: &
        N_PROCESSES_FILE_ACCESS = 128, &
        FILE_ACCESS_DELAY_UNIT = 10000

contains

    subroutine DelayFileAccess(ProcessRank)
        integer, intent(in) :: ProcessRank

        integer :: iDelay, nFileAccessDelayIterations
        real(kind(0d0)) :: Number, Dummy

        nFileAccessDelayIterations &
            = (ProcessRank/N_PROCESSES_FILE_ACCESS)*FILE_ACCESS_DELAY_UNIT

        do iDelay = 1, nFileAccessDelayIterations
            !-- wait my turn
            call random_number(Number)
            Dummy = Number*Number
        end do

    end subroutine DelayFileAccess

end module m_delay_file_access
