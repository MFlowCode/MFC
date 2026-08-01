!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'

!> @brief Checks pre-process input file parameters for compatibility and correctness
module m_checker

    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic
    use m_helper

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the pre_process stage
    impure subroutine s_check_inputs

        integer :: i

        do i = 1, num_patches
            @:PROHIBIT(patch_icpp(i)%rxn_val /= 0._wp .and. .not. jwl_reactive, "patch_icpp(i)%rxn_val requires jwl_reactive")
            @:PROHIBIT(patch_icpp(i)%rxn_val < 0._wp .or. patch_icpp(i)%rxn_val > 1._wp, "patch_icpp(i)%rxn_val must be in [0, 1]")
        end do

    end subroutine s_check_inputs

end module m_checker
