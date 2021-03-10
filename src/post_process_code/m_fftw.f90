!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /
!!    / /  / / __/ / /___
!!   /_/  /_/_/    \____/
!!
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names
!!  are listed in the copyright file included with this source
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published
!!  by the Free Software Foundation, either version 3 of the license
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file m_fftw.f90
!! @brief Contains module m_fftw
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The module contains the subroutines for the FFT routines
module m_fftw

    ! Dependencies =============================================================
    use, intrinsic :: iso_c_binding

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_fftw_module, &
 s_apply_fourier_decomposition, &
 s_finalize_fftw_module

    include 'fftw3.f03'

    type(c_ptr) :: fwd_plan, bwd_plan
    type(c_ptr) :: fftw_real_data, fftw_cmplx_data, fftw_fltr_cmplx_data
    integer :: real_size, cmplx_size
    ! Real data
    real(c_double), pointer :: data_real(:)
    ! Complex data in Fourier space
    complex(c_double_complex), pointer :: data_cmplx(:)
    ! Filtered complex data in Fourier space
    complex(c_double_complex), pointer :: data_fltr_cmplx(:)

contains

    !>  The purpose of this subroutine is to create the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
    subroutine s_initialize_fftw_module() ! ----------------------------------

        ! Size of input array going into DFT
        real_size = p + 1
        ! Size of output array coming out of DFT
        cmplx_size = (p + 1)/2 + 1

        ! Allocate input and output DFT data sizes
        fftw_real_data = fftw_alloc_real(int(real_size, c_size_t))
        fftw_cmplx_data = fftw_alloc_complex(int(cmplx_size, c_size_t))
        fftw_fltr_cmplx_data = fftw_alloc_complex(int(cmplx_size, c_size_t))
        ! Associate input and output data pointers with allocated memory
        call c_f_pointer(fftw_real_data, data_real, [real_size])
        call c_f_pointer(fftw_cmplx_data, data_cmplx, [cmplx_size])
        call c_f_pointer(fftw_fltr_cmplx_data, data_fltr_cmplx, [cmplx_size])

        ! Generate plans for forward and backward DFTs
        fwd_plan = fftw_plan_dft_r2c_1d(real_size, data_real, data_cmplx, FFTW_ESTIMATE)
        bwd_plan = fftw_plan_dft_c2r_1d(real_size, data_fltr_cmplx, data_real, FFTW_ESTIMATE)

    end subroutine s_initialize_fftw_module ! ------------------------------

    !>  The purpose of this subroutine is to Fourier decompose
        !!      the flow field. Not done in most efficient manner since
        !!      subroutine is called for every mode, but can deal with
        !!      efficiency later.
        !!  @param q_sf Scalar field to transform
        !!  @param i Fourier component
    subroutine s_apply_fourier_decomposition(q_sf, i) ! -----------------------

        ! Variable to be Fourier decomposed
        real(kind(0d0)), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(INOUT) :: q_sf

        integer, intent(IN) :: i

        integer :: j, k

        do j = -offset_x%beg, m + offset_x%end
            do k = -offset_y%beg, n + offset_y%end
                data_fltr_cmplx(:) = (0d0, 0d0)
                data_real(1:p + 1) = q_sf(j, k, 0:p)
                call fftw_execute_dft_r2c(fwd_plan, data_real, data_cmplx)
                data_fltr_cmplx(i) = data_cmplx(i)
                call fftw_execute_dft_c2r(bwd_plan, data_fltr_cmplx, data_real)
                data_real(:) = data_real(:)/real(real_size, kind(0d0))
                q_sf(j, k, 0:p) = data_real(1:p + 1)
            end do
        end do

        ! Populate offset regions given that domain is azimuthally periodic
        do j = -offset_z%beg, -1
            q_sf(:, :, j) = q_sf(:, :, (p + 1) + j)
        end do
        do j = 1, offset_z%end
            q_sf(:, :, p + j) = q_sf(:, :, j - 1)
        end do

    end subroutine s_apply_fourier_decomposition ! -------------------------

    !>  The purpose of this subroutine is to destroy the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
    subroutine s_finalize_fftw_module() ! ------------------------------------

        call fftw_free(fftw_real_data)
        call fftw_free(fftw_cmplx_data)
        call fftw_free(fftw_fltr_cmplx_data)

        call fftw_destroy_plan(fwd_plan)
        call fftw_destroy_plan(bwd_plan)

    end subroutine s_finalize_fftw_module ! --------------------------------

end module
