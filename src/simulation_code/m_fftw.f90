!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!! Copyright 2021
!!
!! Permission is hereby granted, free of charge, to any person 
!! obtaining a copy of this software and associated documentation 
!! files (the "Software"), to deal in the Software without 
!! restriction, including without limitation the rights to use, 
!! copy, modify, merge, publish, distribute, sublicense, 
!! and/or sell copies of the Software, and to permit persons 
!! to whom the Software is furnished to do so, subject to the 
!! following conditions:
!!
!! The above copyright notice and this permission notice shall 
!! be included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
!! FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
!! OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
!! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
!! THE SOFTWARE.

!>
!! @file m_fftw.f90
!! @brief Contains module m_fftw
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The module contains the subroutines for the FFT routines
MODULE m_fftw

    ! Dependencies =============================================================
    USE, INTRINSIC :: ISO_C_BINDING

    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy
    ! ==========================================================================

    IMPLICIT NONE

    PRIVATE; PUBLIC :: s_initialize_fftw_module,   &
                       s_apply_fourier_filter,     &
                       s_finalize_fftw_module
    
    INCLUDE 'fftw3.f03'

    TYPE(C_PTR) :: fwd_plan, bwd_plan
    TYPE(C_PTR) :: fftw_real_data, fftw_cmplx_data, fftw_fltr_cmplx_data
    INTEGER :: real_size, cmplx_size

    REAL(C_DOUBLE), POINTER :: data_real(:) !< Real data

    COMPLEX(C_DOUBLE_COMPLEX), POINTER :: data_cmplx(:) !< 
    !! Complex data in Fourier space
    
    

    COMPLEX(C_DOUBLE_COMPLEX), POINTER :: data_fltr_cmplx(:) !< 
    !! Filtered complex data in Fourier space

    CONTAINS




        !>  The purpose of this subroutine is to create the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
        SUBROUTINE s_initialize_fftw_module() ! ----------------------------------



            ! Size of input array going into DFT
            real_size = p+1
            ! Size of output array coming out of DFT
            cmplx_size = (p+1)/2+1

            ! Allocate input and output DFT data sizes
            fftw_real_data       = fftw_alloc_real   (int( real_size, C_SIZE_T))
            fftw_cmplx_data      = fftw_alloc_complex(int(cmplx_size, C_SIZE_T))
            fftw_fltr_cmplx_data = fftw_alloc_complex(int(cmplx_size, C_SIZE_T))
            ! Associate input and output data pointers with allocated memory
            CALL c_f_pointer(fftw_real_data ,      data_real ,      [ real_size])
            CALL c_f_pointer(fftw_cmplx_data,      data_cmplx,      [cmplx_size])
            CALL c_f_pointer(fftw_fltr_cmplx_data, data_fltr_cmplx, [cmplx_size])

            ! Generate plans for forward and backward DFTs
            fwd_plan = fftw_plan_dft_r2c_1d(real_size, data_real      , data_cmplx, FFTW_ESTIMATE)
            bwd_plan = fftw_plan_dft_c2r_1d(real_size, data_fltr_cmplx, data_real , FFTW_ESTIMATE)

        END SUBROUTINE s_initialize_fftw_module ! ------------------------------




        !>  The purpose of this subroutine is to apply a Fourier low-
        !!      pass filter to the flow variables in the azimuthal direction
        !!      to remove the high-frequency content. This alleviates the 
        !!      restrictive CFL condition arising from cells near the axis.
        SUBROUTINE s_apply_fourier_filter(q_cons_vf) ! --------------------------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf

            INTEGER :: Nfq !< Number of kept modes

            INTEGER :: i,j,k !< Generic loop iterators

            ! Restrict filter to processors that have cells adjacent to axis
            IF (bc_y%beg >= 0) RETURN

            ! Keeping only the mean value for cells directly adjacent to axis
            Nfq = 3
            DO j = 0, m
                DO k = 1, sys_size
                    data_fltr_cmplx(:) = (0d0,0d0)
                    data_real(1:p+1) = q_cons_vf(k)%sf(j,0,0:p)
                    CALL fftw_execute_dft_r2c(fwd_plan, data_real, data_cmplx)
                    data_fltr_cmplx(1:Nfq) = data_cmplx(1:Nfq)
                    CALL fftw_execute_dft_c2r(bwd_plan, data_fltr_cmplx, data_real)
                    data_real(:) = data_real(:)/REAL(real_size,KIND(0d0))
                    q_cons_vf(k)%sf(j,0,0:p) = data_real(1:p+1)
                END DO
            END DO

            ! Apply Fourier filter to additional rings
            DO i = 1, fourier_rings
                Nfq = MIN(FLOOR(2d0*REAL(i,KIND(0d0))*pi),cmplx_size)
                DO j = 0, m
                    DO k = 1, sys_size
                        data_fltr_cmplx(:) = (0d0,0d0)
                        data_real(1:p+1) = q_cons_vf(k)%sf(j,i,0:p)
                        CALL fftw_execute_dft_r2c(fwd_plan, data_real, data_cmplx)
                        data_fltr_cmplx(1:Nfq) = data_cmplx(1:Nfq)
                        CALL fftw_execute_dft_c2r(bwd_plan, data_fltr_cmplx, data_real)
                        data_real(:) = data_real(:)/REAL(real_size,KIND(0d0))
                        q_cons_vf(k)%sf(j,i,0:p) = data_real(1:p+1)
                    END DO
                END DO
            END DO

        END SUBROUTINE s_apply_fourier_filter ! --------------------------------




        !>  The purpose of this subroutine is to destroy the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
        SUBROUTINE s_finalize_fftw_module() ! ------------------------------------
            CALL fftw_free(fftw_real_data)
            CALL fftw_free(fftw_cmplx_data)
            CALL fftw_free(fftw_fltr_cmplx_data)

            CALL fftw_destroy_plan(fwd_plan)
            CALL fftw_destroy_plan(bwd_plan)

        END SUBROUTINE s_finalize_fftw_module ! --------------------------------



END MODULE
