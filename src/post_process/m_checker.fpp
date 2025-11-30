!>
!!@file m_checker.f90
!!@brief Contains module m_checker

#:include 'macros.fpp'

!> @brief The purpose of the module is to check for compatible input files
module m_checker

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    implicit none

    private; public :: s_check_inputs, s_check_inputs_fft

contains

    !> Checks compatibility of parameters in the input file.
        !! Used by the post_process stage
    impure subroutine s_check_inputs

        call s_check_inputs_output_format

    end subroutine s_check_inputs

    !> Checks constraints on output format parameters
    impure subroutine s_check_inputs_output_format
        @:PROHIBIT(precision == 2 .and. wp == sp)
    end subroutine s_check_inputs_output_format

    !> Checks constraints on fft_wrt
    impure subroutine s_check_inputs_fft
        integer :: num_procs_y, num_procs_z

        @:PROHIBIT(fft_wrt .and. MOD(n_glb+1,n+1) /= 0, "FFT WRT requires n_glb to be divisible by num_procs_y")
        @:PROHIBIT(fft_wrt .and. MOD(p_glb+1,p+1) /= 0, "FFT WRT requires p_glb to be divisible by num_procs_z")
        num_procs_y = (n_glb + 1)/(n + 1)
        num_procs_z = (p_glb + 1)/(p + 1)
        @:PROHIBIT(fft_wrt .and. MOD(m_glb+1,num_procs_y) /= 0, "FFT WRT requires m_glb to be divisible by num_procs_y")
        @:PROHIBIT(fft_wrt .and. MOD(n_glb+1,num_procs_z) /= 0, "FFT WRT requires n_glb to be divisible by num_procs_z")
    end subroutine s_check_inputs_fft

end module m_checker
