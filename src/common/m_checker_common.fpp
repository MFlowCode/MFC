!>
!!@file
!!@brief Contains module m_checker_common

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Shared input validation checks for grid dimensions and AMD GPU compiler limits
module m_checker_common

    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic
    use m_helper

    implicit none

    private; public :: s_check_inputs_common

contains

    !> Checks compatibility of parameters in the input file. Used by all three stages
    impure subroutine s_check_inputs_common(check_total_cells, n_global)

        logical, intent(in)         :: check_total_cells
        integer(kind=8), intent(in) :: n_global

        if (check_total_cells) call s_check_total_cells(n_global)
        #:if MFC_FIXED_BOUNDS
            call s_check_fixed_bounds
        #:endif

    end subroutine s_check_inputs_common

    !> Verify that the total number of grid cells meets the minimum required by the number of dimensions and MPI ranks.
    impure subroutine s_check_total_cells(n_global)

        character(len=18)           :: numStr  !< for int to string conversion
        integer(kind=8)             :: min_cells
        integer(kind=8), intent(in) :: n_global

        min_cells = int(2, kind=8)**int(min(1, m) + min(1, n) + min(1, p), kind=8)*int(num_procs, kind=8)
        call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)

        @:PROHIBIT(n_global < min_cells, &
                   & "Total number of cells must be at least (2^[number of dimensions])*num_procs, " // "which is currently " &
                   & // trim(numStr))

    end subroutine s_check_total_cells

    !> Check that the case fits the fixed per-thread array bounds of a GPU simulation build.
    impure subroutine s_check_fixed_bounds

        #:if not MFC_CASE_OPTIMIZATION
            @:PROHIBIT(num_fluids > ${NUM_FLUIDS_MAX}$, &
                       & "num_fluids <= ${NUM_FLUIDS_MAX}$ in GPU builds; rebuild with --case-optimization")
            @:PROHIBIT((bubbles_euler .or. bubbles_lagrange) .and. nb > ${NB_MAX}$, &
                       & "nb <= ${NB_MAX}$ in GPU builds; rebuild with --case-optimization")
        #:endif

    end subroutine s_check_fixed_bounds

end module m_checker_common
