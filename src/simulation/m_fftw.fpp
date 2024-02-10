!>
!! @file m_fftw.f90
!! @brief Contains module m_fftw

#:include 'macros.fpp'

!> @brief The module contains the subroutines for the FFT routines
module m_fftw

    ! Dependencies =============================================================
    use, intrinsic :: iso_c_binding

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

#if defined(MFC_OpenACC) && defined(__PGI)
    use cufft
#endif
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_fftw_module, &
 s_apply_fourier_filter, &
 s_finalize_fftw_module

#if !(defined(MFC_OpenACC) && defined(__PGI))
    include 'fftw3.f03'
#endif

    type(c_ptr) :: fwd_plan, bwd_plan
    type(c_ptr) :: fftw_real_data, fftw_cmplx_data, fftw_fltr_cmplx_data
    integer :: real_size, cmplx_size, x_size, batch_size

    real(c_double), pointer :: data_real(:) !< Real data

    complex(c_double_complex), pointer :: data_cmplx(:) !<
    !! Complex data in Fourier space

    complex(c_double_complex), pointer :: data_fltr_cmplx(:) !<
    !! Filtered complex data in Fourier space

#if defined(MFC_OpenACC) && defined(__PGI)
    !$acc declare create(real_size, cmplx_size, x_size, batch_size)

    real(kind(0d0)), allocatable :: data_real_gpu(:)
    complex(kind(0d0)), allocatable :: data_cmplx_gpu(:)
    complex(kind(0d0)), allocatable :: data_fltr_cmplx_gpu(:)
    !$acc declare create(data_real_gpu, data_cmplx_gpu, data_fltr_cmplx_gpu)

    integer :: fwd_plan_gpu, bwd_plan_gpu, ierr

    integer, allocatable :: cufft_size(:), iembed(:), oembed(:)

    integer :: istride, ostride, idist, odist, rank
#endif

contains

    !>  The purpose of this subroutine is to create the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
    subroutine s_initialize_fftw_module() ! ----------------------------------

#if defined(MFC_OpenACC) && !defined(__PGI)

        print *, "The FFTW module is not supported when using OpenACC with a compiler other than NVHPC/PGI."
        stop 1

#endif

        ! Size of input array going into DFT
        real_size = p + 1
        ! Size of output array coming out of DFT
        cmplx_size = (p + 1)/2 + 1

        x_size = m + 1

        batch_size = x_size*sys_size

#if defined(MFC_OpenACC) && defined(__PGI)
        rank = 1; istride = 1; ostride = 1

        allocate (cufft_size(1:rank), iembed(1:rank), oembed(1:rank))

        cufft_size(1) = real_size; 
        iembed(1) = 0
        oembed(1) = 0

!$acc update device(real_size, cmplx_size, x_size, sys_size, batch_size)
#else
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
#endif

#if defined(MFC_OpenACC) && defined(__PGI)
        @:ALLOCATE(data_real_gpu(1:real_size*x_size*sys_size))
        @:ALLOCATE(data_cmplx_gpu(1:cmplx_size*x_size*sys_size))
        @:ALLOCATE(data_fltr_cmplx_gpu(1:cmplx_size*x_size*sys_size))

        ierr = cufftPlanMany(fwd_plan_gpu, rank, cufft_size, iembed, istride, real_size, oembed, ostride, cmplx_size, CUFFT_D2Z, batch_size)
        ierr = cufftPlanMany(bwd_plan_gpu, rank, cufft_size, iembed, istride, cmplx_size, oembed, ostride, real_size, CUFFT_Z2D, batch_size)
#endif

    end subroutine s_initialize_fftw_module ! ------------------------------

    !>  The purpose of this subroutine is to apply a Fourier low-
        !!      pass filter to the flow variables in the azimuthal direction
        !!      to remove the high-frequency content. This alleviates the
        !!      restrictive CFL condition arising from cells near the axis.
    subroutine s_apply_fourier_filter(q_cons_vf) ! --------------------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf

        integer :: Nfq !< Number of kept modes

        integer :: i, j, k, l !< Generic loop iterators

        ! Restrict filter to processors that have cells adjacent to axis
        if (bc_y%beg >= 0) return

#if defined(MFC_OpenACC) && defined(__PGI)

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 1, sys_size
            do j = 0, m
                do l = 1, cmplx_size
                    data_fltr_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size) = (0d0, 0d0)
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 1, sys_size
            do j = 0, m
                do l = 0, p
                    data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size) = q_cons_vf(k)%sf(j, 0, l)
                end do
            end do
        end do

        !$acc host_data use_device(data_real_gpu, data_cmplx_gpu)
        ierr = cufftExecD2Z(fwd_plan_gpu, data_real_gpu, data_cmplx_gpu)
        !$acc end host_data

        Nfq = 3

        !$acc parallel loop collapse(3) gang vector default(present) firstprivate(Nfq)
        do k = 1, sys_size
            do j = 0, m
                do l = 1, Nfq
                    data_fltr_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size) = data_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size)
                end do
            end do
        end do

        !$acc host_data use_device(data_real_gpu, data_fltr_cmplx_gpu)
        ierr = cufftExecZ2D(bwd_plan_gpu, data_fltr_cmplx_gpu, data_real_gpu)
        !$acc end host_data

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 1, sys_size
            do j = 0, m
                do l = 0, p
                    data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size) = data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size)/real(real_size, kind(0d0))
                    q_cons_vf(k)%sf(j, 0, l) = data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size)
                end do
            end do
        end do

        do i = 1, fourier_rings

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 1, sys_size
                do j = 0, m
                    do l = 1, cmplx_size
                        data_fltr_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size) = (0d0, 0d0)
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present) firstprivate(i)
            do k = 1, sys_size
                do j = 0, m
                    do l = 0, p
                        data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size) = q_cons_vf(k)%sf(j, i, l)
                    end do
                end do
            end do

            !$acc host_data use_device(data_real_gpu, data_cmplx_gpu)
            ierr = cufftExecD2Z(fwd_plan_gpu, data_real_gpu, data_cmplx_gpu)
            !$acc end host_data

            Nfq = min(floor(2d0*real(i, kind(0d0))*pi), cmplx_size)

            !$acc parallel loop collapse(3) gang vector default(present) firstprivate(Nfq)
            do k = 1, sys_size
                do j = 0, m
                    do l = 1, Nfq
                        data_fltr_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size) = data_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size)
                    end do
                end do
            end do

            !$acc host_data use_device(data_real_gpu, data_fltr_cmplx_gpu)
            ierr = cufftExecZ2D(bwd_plan_gpu, data_fltr_cmplx_gpu, data_real_gpu)
            !$acc end host_data

            !$acc parallel loop collapse(3) gang vector default(present) firstprivate(i)
            do k = 1, sys_size
                do j = 0, m
                    do l = 0, p
                        data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size) = data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size)/real(real_size, kind(0d0))
                        q_cons_vf(k)%sf(j, i, l) = data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size)
                    end do
                end do
            end do

        end do

#else
        Nfq = 3
        do j = 0, m
            do k = 1, sys_size
                data_fltr_cmplx(:) = (0d0, 0d0)
                data_real(1:p + 1) = q_cons_vf(k)%sf(j, 0, 0:p)
                call fftw_execute_dft_r2c(fwd_plan, data_real, data_cmplx)
                data_fltr_cmplx(1:Nfq) = data_cmplx(1:Nfq)
                call fftw_execute_dft_c2r(bwd_plan, data_fltr_cmplx, data_real)
                data_real(:) = data_real(:)/real(real_size, kind(0d0))
                q_cons_vf(k)%sf(j, 0, 0:p) = data_real(1:p + 1)
            end do
        end do

        ! Apply Fourier filter to additional rings
        do i = 1, fourier_rings
            Nfq = min(floor(2d0*real(i, kind(0d0))*pi), cmplx_size)
            do j = 0, m
                do k = 1, sys_size
                    data_fltr_cmplx(:) = (0d0, 0d0)
                    data_real(1:p + 1) = q_cons_vf(k)%sf(j, i, 0:p)
                    call fftw_execute_dft_r2c(fwd_plan, data_real, data_cmplx)
                    data_fltr_cmplx(1:Nfq) = data_cmplx(1:Nfq)
                    call fftw_execute_dft_c2r(bwd_plan, data_fltr_cmplx, data_real)
                    data_real(:) = data_real(:)/real(real_size, kind(0d0))
                    q_cons_vf(k)%sf(j, i, 0:p) = data_real(1:p + 1)
                end do
            end do
        end do
#endif

    end subroutine s_apply_fourier_filter ! --------------------------------

    !>  The purpose of this subroutine is to destroy the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
    subroutine s_finalize_fftw_module() ! ------------------------------------

#if defined(MFC_OpenACC) && defined(__PGI)
        @:DEALLOCATE(data_real_gpu, data_fltr_cmplx_gpu, data_cmplx_gpu)
        ierr = cufftDestroy(fwd_plan_gpu)
        ierr = cufftDestroy(bwd_plan_gpu)
#else
        call fftw_free(fftw_real_data)
        call fftw_free(fftw_cmplx_data)
        call fftw_free(fftw_fltr_cmplx_data)

        call fftw_destroy_plan(fwd_plan)
        call fftw_destroy_plan(bwd_plan)
#endif

    end subroutine s_finalize_fftw_module ! --------------------------------

end module
