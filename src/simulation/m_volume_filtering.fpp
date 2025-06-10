#:include 'macros.fpp'

module m_volume_filtering

    use, intrinsic :: iso_c_binding

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_ibm

    use m_boundary_common

#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

#if defined(MFC_OpenACC) && defined(__PGI)
    use cufft
#endif

    implicit none

    private; public :: s_initialize_fftw_explicit_filter_module, &
 s_initialize_filtering_kernel, s_initialize_fluid_indicator_function, s_initialize_filtered_fluid_indicator_function, & 
 s_finalize_fftw_explicit_filter_module, & 
 s_apply_fftw_filter_cons, s_apply_fftw_filter_tensor, s_apply_fftw_filter_scalarfield, &
 s_mpi_transpose_slabZ2Y, s_mpi_transpose_slabY2Z, s_mpi_FFT_fwd, s_mpi_FFT_bwd, &
 s_setup_terms_filtering, s_compute_pseudo_turbulent_reynolds_stress, s_compute_R_mu, s_compute_interphase_momentum_exchange_term

#if !defined(MFC_OpenACC)
    include 'fftw3.f03'
#endif

    integer :: ierr   

    ! fluid indicator function (1 = fluid, 0 = otherwise)
    type(scalar_field), public :: fluid_indicator_function_I

    !$acc declare create(fluid_indicator_function_I)

#if defined(MFC_OpenACC)
    ! GPU plans
    integer :: plan_x_fwd_gpu, plan_x_bwd_gpu, plan_y_gpu, plan_z_gpu
#else
    ! CPU plans
    type(c_ptr) :: plan_x_r2c_fwd, plan_x_c2r_bwd
    type(c_ptr) :: plan_y_c2c_fwd, plan_y_c2c_bwd 
    type(c_ptr) :: plan_z_c2c_fwd, plan_z_c2c_bwd
    type(c_ptr) :: plan_x_r2c_kernelG, plan_y_c2c_kernelG, plan_z_c2c_kernelG
#endif

    ! domain size information (global, complex, local)
    integer :: Nx, Ny, Nz, NxC, Nyloc, Nzloc

    ! 1D real and complex vectors for FFT routines
    real(c_double), allocatable :: data_real_in1d(:) 
    complex(c_double_complex), allocatable :: data_cmplx_out1d(:)
    complex(c_double_complex), allocatable :: data_cmplx_out1dy(:)

    ! 3D arrays for slab transposes
    complex(c_double_complex), allocatable :: data_cmplx_slabz(:, :, :), data_cmplx_slaby(:, :, :)

    ! input/output array for FFT routine
    real(c_double), allocatable :: data_real_3D_slabz(:, :, :)

    ! filtering kernel in physical space
    real(c_double), allocatable :: real_kernelG_in(:, :, :)

    ! FFT of filtering kernel
    complex(c_double_complex), allocatable :: cmplx_kernelG1d(:)

    !$acc declare create(Nx, Ny, Nz, NxC, Nyloc, Nzloc)
    !$acc declare create(data_real_in1d, data_cmplx_out1d, data_cmplx_out1dy, data_cmplx_slabz, data_cmplx_slaby, data_real_3D_slabz, real_kernelG_in, cmplx_kernelG1d)

contains

    !< create fft plans to be used for explicit filtering of data 
    subroutine s_initialize_fftw_explicit_filter_module
        integer :: size_n(1), inembed(1), onembed(1)

        !< global sizes 
        Nx = m_glb + 1
        Ny = n_glb + 1
        Nz = p_glb + 1

        !< complex size
        NxC = Nx/2 + 1

        !< local sizes on each processor
        Nyloc = Ny / num_procs
        Nzloc = p + 1

        !$acc update device(Nx, Ny, Nz, NxC, Nyloc, Nzloc)

        @:ALLOCATE(data_real_in1d(Nx*Ny*Nzloc))
        @:ALLOCATE(data_cmplx_out1d(NxC*Ny*Nz/num_procs))
        @:ALLOCATE(data_cmplx_out1dy(NxC*Ny*Nz/num_procs))
        @:ALLOCATE(cmplx_kernelG1d(NxC*Nyloc*Nz))
        @:ALLOCATE(real_kernelG_in(Nx, Ny, Nzloc))
        @:ALLOCATE(data_real_3D_slabz(Nx, Ny, Nzloc))
        @:ALLOCATE(data_cmplx_slabz(NxC, Ny, Nzloc))
        @:ALLOCATE(data_cmplx_slaby(NxC, Nyloc, Nz))

#if defined(MFC_OpenACC)
        !< GPU FFT plans
        !< X - plans
        size_n(1) = Nx
        inembed(1) = Nx
        onembed(1) = NxC
        ierr = cufftPlanMany(plan_x_fwd_gpu, 1, size_n, inembed, 1, Nx, onembed, 1, NxC, CUFFT_D2Z, Ny*Nzloc)
        size_n(1) = Nx
        inembed(1) = NxC
        onembed(1) = Nx  
        ierr = cufftPlanMany(plan_x_bwd_gpu, 1, size_n, inembed, 1, NxC, onembed, 1, Nx, CUFFT_Z2D, Ny*Nzloc)
        !< Y - plans
        size_n(1) = Ny
        inembed(1) = Ny
        onembed(1) = Ny
        ierr = cufftPlanMany(plan_y_gpu, 1, size_n, inembed, 1, Ny, onembed, 1, Ny, CUFFT_Z2Z, NxC*Nzloc)
        !< Z - plans
        size_n(1) = Nz 
        inembed(1) = Nz 
        onembed(1) = Nz 
        ierr = cufftPlanMany(plan_z_gpu, 1, size_n, inembed, 1, Nz, onembed, 1, Nz, CUFFT_Z2Z, NxC*Nyloc)
#else
        !< CPU FFT plans
        !< X - direction plans
        size_n(1) = Nx
        inembed(1) = Nx
        onembed(1) = NxC
        plan_x_r2c_fwd = fftw_plan_many_dft_r2c(1, size_n, Ny*Nzloc, &                  ! rank, n, howmany
                                                data_real_in1d, inembed, 1, Nx, &       ! in, inembed, istride, idist
                                                data_cmplx_out1d, onembed, 1, NxC, &    ! out, onembed, ostride, odist
                                                FFTW_MEASURE)                           ! sign, flags
        size_n(1) = Nx
        inembed(1) = NxC
        onembed(1) = Nx                                                         
        plan_x_c2r_bwd = fftw_plan_many_dft_c2r(1, size_n, Ny*Nzloc, & 
                                                data_cmplx_out1d, inembed, 1, NxC, & 
                                                data_real_in1d, onembed, 1, Nx, & 
                                                FFTW_MEASURE)
        !< Y - direction plans
        size_n(1) = Ny
        inembed(1) = Ny
        onembed(1) = Ny
        plan_y_c2c_fwd = fftw_plan_many_dft(1, size_n, NxC*Nzloc, & 
                                            data_cmplx_out1dy, inembed, 1, Ny, & 
                                            data_cmplx_out1dy, onembed, 1, Ny, & 
                                            FFTW_FORWARD, FFTW_MEASURE)
        plan_y_c2c_bwd = fftw_plan_many_dft(1, size_n, NxC*Nzloc, & 
                                            data_cmplx_out1dy, inembed, 1, Ny, & 
                                            data_cmplx_out1dy, onembed, 1, Ny, & 
                                            FFTW_BACKWARD, FFTW_MEASURE)
        !< Z - direction plans
        size_n(1) = Nz 
        inembed(1) = Nz 
        onembed(1) = Nz 
        plan_z_c2c_fwd = fftw_plan_many_dft(1, size_n, NxC*Nyloc, & 
                                            data_cmplx_out1d, inembed, 1, Nz, & 
                                            data_cmplx_out1d, onembed, 1, Nz, & 
                                            FFTW_FORWARD, FFTW_MEASURE)
        plan_z_c2c_bwd = fftw_plan_many_dft(1, size_n, NxC*Nyloc, & 
                                            data_cmplx_out1d, inembed, 1, Nz, &
                                            data_cmplx_out1d, onembed, 1, Nz, & 
                                            FFTW_BACKWARD, FFTW_MEASURE)
        ! forward plans for filtering kernel
        ! X kernel plan
        size_n(1) = Nx
        inembed(1) = Nx
        onembed(1) = NxC
        plan_x_r2c_kernelG = fftw_plan_many_dft_r2c(1, size_n, Ny*Nzloc, &                    
                                                    data_real_in1d, inembed, 1, Nx, &        
                                                    cmplx_kernelG1d, onembed, 1, NxC, &    
                                                    FFTW_MEASURE)          
        ! Y kernel plan                  
        size_n(1) = Ny
        inembed(1) = Ny
        onembed(1) = Ny
        plan_y_c2c_kernelG = fftw_plan_many_dft(1, size_n, NxC*Nzloc, & 
                                                data_cmplx_out1dy, inembed, 1, Ny, & 
                                                data_cmplx_out1dy, onembed, 1, Ny, & 
                                                FFTW_FORWARD, FFTW_MEASURE)
        ! Z kernel plan
        size_n(1) = Nz 
        inembed(1) = Nz 
        onembed(1) = Nz 
        plan_z_c2c_kernelG = fftw_plan_many_dft(1, size_n, NxC*Nyloc, & 
                                                cmplx_kernelG1d, inembed, 1, Nz, & 
                                                cmplx_kernelG1d, onembed, 1, Nz, & 
                                                FFTW_FORWARD, FFTW_MEASURE)
#endif
    end subroutine s_initialize_fftw_explicit_filter_module

    !< initialize the gaussian filtering kernel in real space and then compute its DFT
    subroutine s_initialize_filtering_kernel
        real(dp) :: sigma_stddev
        real(dp) :: Lx, Ly, Lz
        real(dp) :: x_r, y_r, z_r  
        real(dp) :: r2
        real(dp) :: G_norm_int, G_norm_int_glb
        integer :: i, j, k, idx

        ! gaussian filter
        sigma_stddev = 3.0_dp * 0.05_dp

        Lx = x_domain_end_glb - x_domain_beg_glb
        Ly = y_domain_end_glb - y_domain_beg_glb  
        Lz = z_domain_end_glb - z_domain_beg_glb    
        
        G_norm_int = 0.0_dp
   
        !$acc parallel loop collapse(3) gang vector default(present) reduction(+:G_norm_int) copyin(Lx, Ly, Lz, sigma_stddev) private(x_r, y_r, z_r, r2)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p 
                    x_r = min(abs(x_cc(i) - x_domain_beg_glb), Lx - abs(x_cc(i) - x_domain_beg_glb))
                    y_r = min(abs(y_cc(j) - y_domain_beg_glb), Ly - abs(y_cc(j) - y_domain_beg_glb))
                    z_r = min(abs(z_cc(k) - z_domain_beg_glb), Lz - abs(z_cc(k) - z_domain_beg_glb))

                    r2 = x_r**2 + y_r**2 + z_r**2

                    real_kernelG_in(i+1, j+1, k+1) = exp(-r2 / (2.0_dp*sigma_stddev**2))

                    G_norm_int = G_norm_int + real_kernelG_in(i+1, j+1, k+1)*dx(i)*dy(j)*dz(k)
                end do 
            end do
        end do

        call s_mpi_allreduce_sum(G_norm_int, G_norm_int_glb) 

        ! FFT of kernel
        ! normalize kernel
        !$acc parallel loop collapse(3) gang vector default(present) copyin(G_norm_int_glb)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_real_3D_slabz(i, j, k) = real_kernelG_in(i, j, k) / G_norm_int_glb
                end do 
            end do 
        end do 

        ! 3D z-slab -> 1D x, y, z
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_real_in1d(i + (j-1)*Nx + (k-1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                end do 
            end do 
        end do

        ! X FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, cmplx_kernelG1d)
#else
        call fftw_execute_dft_r2c(plan_x_r2c_kernelG, data_real_in1d, cmplx_kernelG1d)
#endif

        ! 1D x, y, z -> 1D y, x, z (CMPLX)
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC) = cmplx_kernelG1d(i + (j-1)*NxC + (k-1)*NxC*Ny)
                end do 
            end do 
        end do

        ! Y FFT 
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_y_c2c_kernelG, data_cmplx_out1dy, data_cmplx_out1dy)
#endif

        ! 1D y, x, z -> 3D z-slab
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_slabz(i, j, k) = data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC)
                end do 
            end do 
        end do 

        ! transpose z-slab to y-slab
        call s_mpi_transpose_slabZ2Y 

        ! 3D y-slab -> 1D z, x, y
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz
                    cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_slaby(i, j, k)
                end do 
            end do 
        end do

        ! Z FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_z_gpu, cmplx_kernelG1d, cmplx_kernelG1d, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_z_c2c_kernelG, cmplx_kernelG1d, cmplx_kernelG1d)
#endif

        ! normalize FFT 
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz
                    cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC) / (real(Nx*Ny*Nz, dp))
                end do 
            end do 
        end do

        ! return cmplx_kernelG1d: 1D z, x, y
    end subroutine s_initialize_filtering_kernel

    !< initialize fluid indicator function
    subroutine s_initialize_fluid_indicator_function 
        integer :: i, j, k 

        @:ALLOCATE(fluid_indicator_function_I%sf(0:m, 0:n, 0:p))
        @:ACC_SETUP_SFs(fluid_indicator_function_I)

        ! define fluid indicator function
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n 
                do k = 0, p
                    if (ib_markers%sf(i, j, k) == 0) then 
                        fluid_indicator_function_I%sf(i, j, k) = 1.0_dp
                    else 
                        fluid_indicator_function_I%sf(i, j, k) = 0.0_dp
                    end if
                end do
            end do
        end do

    end subroutine s_initialize_fluid_indicator_function

    !< compute the filtered fluid indicator function counterpart
    subroutine s_initialize_filtered_fluid_indicator_function(filtered_fluid_indicator_function)
        type(scalar_field) :: filtered_fluid_indicator_function

        integer :: i, j, k

        ! filter fluid indicator function -> stored in q_cons_vf(advxb)
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc 
                    data_real_3D_slabz(i, j, k) = fluid_indicator_function_I%sf(i-1, j-1, k-1)
                end do 
            end do 
        end do 

        call s_mpi_FFT_fwd 

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz
                    data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) * cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                end do
            end do 
        end do

        call s_mpi_FFT_bwd

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    filtered_fluid_indicator_function%sf(i-1, j-1, k-1) = data_real_3D_slabz(i, j, k) / (real(Nx*Ny*Nz, dp))
                end do 
            end do
        end do

    end subroutine s_initialize_filtered_fluid_indicator_function

    !< apply the gaussian filter to the conservative variables and compute their filtered components
    subroutine s_apply_fftw_filter_cons(q_cons_vf, q_cons_filtered)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_filtered

        integer :: l

        do l = 1, sys_size-1
            call s_apply_fftw_filter_scalarfield(q_cons_filtered(advxb), .true., q_cons_vf(l), q_cons_filtered(l))
        end do 

    end subroutine s_apply_fftw_filter_cons

    !< applies the gaussian filter to an arbitrary scalar field
    subroutine s_apply_fftw_filter_scalarfield(filtered_fluid_indicator_function, fluid_quantity, q_temp_in, q_temp_out)
        type(scalar_field), intent(in) :: filtered_fluid_indicator_function
        type(scalar_field), intent(inout) :: q_temp_in
        type(scalar_field), intent(inout), optional :: q_temp_out

        logical, intent(in) :: fluid_quantity !< whether or not convolution integral is over V_f or V_p^(i) - integral over fluid volume or particle volume

        integer :: i, j, k

        if (fluid_quantity) then 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p 
                        data_real_3D_slabz(i+1, j+1, k+1) = q_temp_in%sf(i, j, k) * fluid_indicator_function_I%sf(i, j, k)
                    end do 
                end do 
            end do
        else 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p 
                        data_real_3D_slabz(i+1, j+1, k+1) = q_temp_in%sf(i, j, k) * (1.0_dp - fluid_indicator_function_I%sf(i, j, k))
                    end do 
                end do 
            end do
        end if

        call s_mpi_FFT_fwd 

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz 
                    data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) * cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                end do 
            end do 
        end do

        call s_mpi_FFT_bwd

        if (present(q_temp_out)) then 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        q_temp_out%sf(i, j, k) = data_real_3D_slabz(i+1, j+1, k+1) / (real(Nx*Ny*Nz, dp) * filtered_fluid_indicator_function%sf(i, j, k))
                    end do 
                end do 
            end do
        else 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m
                do j = 0, n 
                    do k = 0, p 
                        q_temp_in%sf(i, j, k) = data_real_3D_slabz(i+1, j+1, k+1) / (real(Nx*Ny*Nz, dp) * filtered_fluid_indicator_function%sf(i, j, k))      
                    end do 
                end do 
            end do
        end if

    end subroutine s_apply_fftw_filter_scalarfield

    !< apply the gaussian filter to the requisite tensors to compute unclosed terms of interest
    subroutine s_apply_fftw_filter_tensor(pt_Re_stress, R_mu, q_cons_filtered, div_pres_visc_stress, pres_visc_stress_filtered)
        type(vector_field), dimension(1:num_dims), intent(inout) :: pt_Re_stress
        type(vector_field), dimension(1:num_dims), intent(inout) :: R_mu
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_filtered
        type(scalar_field), dimension(momxb:momxe), intent(inout) :: div_pres_visc_stress
        type(scalar_field), dimension(1:num_dims), intent(inout) :: pres_visc_stress_filtered

        integer :: i, j, k, l, q

        ! pseudo turbulent reynolds stress
        do l = 1, num_dims 
            do q = 1, num_dims
                call s_apply_fftw_filter_scalarfield(q_cons_filtered(advxb), .true., pt_Re_stress(l)%vf(q))
            end do
        end do 

        ! effective viscosity
        do l = 1, num_dims 
            do q = 1, num_dims
                call s_apply_fftw_filter_scalarfield(q_cons_filtered(advxb), .true., R_mu(l)%vf(q))
            end do
        end do 

        ! interphase momentum exchange
        do l = 1, num_dims
            call s_apply_fftw_filter_scalarfield(q_cons_filtered(advxb), .false., div_pres_visc_stress(momxb-1+l), pres_visc_stress_filtered(l))
        end do 

    end subroutine s_apply_fftw_filter_tensor

    !< transpose domain from z-slabs to y-slabs on each processor
    subroutine s_mpi_transpose_slabZ2Y
        complex(c_double_complex), allocatable :: sendbuf(:), recvbuf(:)
        integer :: dest_rank, src_rank
        integer :: i, j, k

        allocate(sendbuf(NxC*Nyloc*Nzloc*num_procs))
        allocate(recvbuf(NxC*Nyloc*Nzloc*num_procs))

        !$acc parallel loop collapse(4) gang vector default(present) copy(sendbuf)
        do dest_rank = 0, num_procs-1
            do k = 1, Nzloc 
                do j = 1, Nyloc
                    do i = 1, NxC
                        sendbuf(i + (j-1)*NxC + (k-1)*NxC*Nyloc + dest_rank*NxC*Nyloc*Nzloc) = data_cmplx_slabz(i, j+dest_rank*Nyloc, k)
                    end do 
                end do
            end do
        end do

        call MPI_Alltoall(sendbuf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, & 
                          recvbuf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop collapse(4) gang vector default(present) copy(recvbuf)
        do src_rank = 0, num_procs-1
            do k = 1, Nzloc 
                do j = 1, Nyloc
                    do i = 1, NxC
                        data_cmplx_slaby(i, j, k+src_rank*Nzloc) = recvbuf(i + (j-1)*NxC + (k-1)*NxC*Nyloc + src_rank*NxC*Nyloc*Nzloc)
                    end do 
                end do
            end do 
        end do

        deallocate(sendbuf, recvbuf)
    end subroutine s_mpi_transpose_slabZ2Y

    !< transpose domain from y-slabs to z-slabs on each processor
    subroutine s_mpi_transpose_slabY2Z 
        complex(c_double_complex), allocatable :: sendbuf(:), recvbuf(:)
        integer :: dest_rank, src_rank
        integer :: i, j, k

        allocate(sendbuf(NxC*Nyloc*Nzloc*num_procs))
        allocate(recvbuf(NxC*Nyloc*Nzloc*num_procs))

        !$acc parallel loop collapse(4) gang vector default(present) copy(sendbuf)
        do dest_rank = 0, num_procs-1
            do k = 1, Nzloc 
                do j = 1, Nyloc 
                    do i = 1, NxC 
                        sendbuf(i + (j-1)*NxC + (k-1)*NxC*Nyloc + dest_rank*NxC*Nyloc*Nzloc) = data_cmplx_slaby(i, j, k+dest_rank*Nzloc)
                    end do 
                end do 
            end do 
        end do

        call MPI_Alltoall(sendbuf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, & 
                          recvbuf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop collapse(4) gang vector default(present) copy(recvbuf) 
        do src_rank = 0, num_procs-1
            do k = 1, Nzloc
                do j = 1, Nyloc 
                    do i = 1, NxC 
                        data_cmplx_slabz(i, j+src_rank*Nyloc, k) = recvbuf(i + (j-1)*NxC + (k-1)*NxC*Nyloc + src_rank*NxC*Nyloc*Nzloc)
                    end do 
                end do
            end do 
        end do
        
        deallocate(sendbuf, recvbuf)
    end subroutine s_mpi_transpose_slabY2Z

    !< compute forward FFT, input: data_real_3D_slabz, output: data_cmplx_out1d
    subroutine s_mpi_FFT_fwd
        integer :: i, j, k

        ! 3D z-slab -> 1D x, y, z
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_real_in1d(i + (j-1)*Nx + (k-1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                end do 
            end do 
        end do

        ! X FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, data_cmplx_out1d)
#else
        call fftw_execute_dft_r2c(plan_x_r2c_fwd, data_real_in1d, data_cmplx_out1d)
#endif

        ! 1D x, y, z -> 1D y, x, z (CMPLX)
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC) = data_cmplx_out1d(i + (j-1)*NxC + (k-1)*NxC*Ny)
                end do 
            end do 
        end do

        ! Y FFT 
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_y_c2c_fwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif 

        ! 1D y, x, z -> 3D z-slab
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_slabz(i, j, k) = data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC)
                end do 
            end do 
        end do 

        ! transpose z-slab to y-slab
        call s_mpi_transpose_slabZ2Y 

        ! 3D y-slab -> 1D z, x, y
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz
                    data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_slaby(i, j, k)
                end do 
            end do 
        end do

        ! Z FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_z_c2c_fwd, data_cmplx_out1d, data_cmplx_out1d)
#endif

        ! return data_cmplx_out1d: 1D z, x, y
    end subroutine s_mpi_FFT_fwd

    !< compute inverse FFT, input: data_cmplx_out1d, output: data_real_3D_slabz
    subroutine s_mpi_FFT_bwd
        integer :: i, j, k

        ! Z inv FFT 
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_INVERSE)
#else
        call fftw_execute_dft(plan_z_c2c_bwd, data_cmplx_out1d, data_cmplx_out1d)
#endif

        ! 1D z, x, y -> 3D y-slab
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz 
                    data_cmplx_slaby(i, j, k) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                end do 
            end do 
        end do

        ! transpose y-slab to z-slab
        call s_mpi_transpose_slabY2Z

        ! 3D z-slab -> 1D y, x, z
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC) = data_cmplx_slabz(i, j, k)
                end do 
            end do 
        end do

        ! Y inv FFT 
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_INVERSE)
#else
        call fftw_execute_dft(plan_y_c2c_bwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif

        ! 1D y, x, z -> 1D x, y, z 
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_out1d(i + (j-1)*NxC + (k-1)*NxC*Ny) = data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC)
                end do 
            end do 
        end do

        ! X inv FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2D(plan_x_bwd_gpu, data_cmplx_out1d, data_real_in1d)
#else
        call fftw_execute_dft_c2r(plan_x_c2r_bwd, data_cmplx_out1d, data_real_in1d)
#endif

        ! 1D x, y, z -> 3D z-slab
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_real_3D_slabz(i, j, k) = data_real_in1d(i + (j-1)*Nx + (k-1)*Nx*Ny)
                end do 
            end do 
        end do

    end subroutine s_mpi_FFT_bwd

    !< setup for calculation of unclosed terms in volume filtered momentum eqn
    subroutine s_setup_terms_filtering(q_cons_vf, pt_Re_stress, R_mu)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(vector_field), dimension(1:num_dims), intent(inout) :: pt_Re_stress
        type(vector_field), dimension(1:num_dims), intent(inout) :: R_mu

        integer :: i, j, k, l, q

        ! pseudo turbulent reynolds stress setup
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    !$acc loop seq
                    do l = 1, num_dims
                        !$acc loop seq
                        do q = 1, num_dims
                            pt_Re_stress(l)%vf(q)%sf(i, j, k) = (q_cons_vf(momxb-1+l)%sf(i, j, k) * q_cons_vf(momxb-1+q)%sf(i, j, k)) / q_cons_vf(1)%sf(i, j, k) ! (rho*u x rho*u)/rho = rho*(u x u) 
                        end do
                    end do
                end do
            end do 
        end do

        ! set density and momentum buffers
#ifdef MFC_MPI
        do i = 1, momxe 
            call s_populate_scalarfield_buffers(q_cons_vf(i))
        end do
#else
        do i = 1, momxe
            q_cons_vf(i)%sf(-buff_size:-1, :, :) = q_cons_vf(i)%sf(m-buff_size+1:m, :, :)
            q_cons_vf(i)%sf(m+1:m+buff_size, :, :) = q_cons_vf(i)%sf(0:buff_size-1, :, :)

            q_cons_vf(i)%sf(:, -buff_size:-1, :) = q_cons_vf(i)%sf(:, n-buff_size+1:n, :)
            q_cons_vf(i)%sf(:, n+1:n+buff_size, :) = q_cons_vf(i)%sf(:, 0:buff_size-1, :)

            q_cons_vf(i)%sf(:, :, -buff_size:-1) = q_cons_vf(i)%sf(:, :, p-buff_size+1:p)
            q_cons_vf(i)%sf(:, :, p+1:p+buff_size) = q_cons_vf(i)%sf(:, :, 0:buff_size-1)
        end do
#endif
        
        ! R_mu setup
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    R_mu(1)%vf(1)%sf(i, j, k) = mu_visc * (2._wp*(q_cons_vf(momxb)%sf(i+1, j, k)/q_cons_vf(1)%sf(i+1, j, k) - q_cons_vf(momxb)%sf(i-1, j, k)/q_cons_vf(1)%sf(i-1, j, k))/(2._wp*dx(i)) & 
                                                - 2._wp/3._wp*((q_cons_vf(momxb)%sf(i+1, j, k)/q_cons_vf(1)%sf(i+1, j, k) - q_cons_vf(momxb)%sf(i-1, j, k)/q_cons_vf(1)%sf(i-1, j, k))/(2._wp*dx(i)) & 
                                                + (q_cons_vf(momxb+1)%sf(i, j+1, k)/q_cons_vf(1)%sf(i, j+1, k) - q_cons_vf(momxb+1)%sf(i, j-1, k)/q_cons_vf(1)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                                + (q_cons_vf(momxb+2)%sf(i, j, k+1)/q_cons_vf(1)%sf(i, j, k+1) - q_cons_vf(momxb+2)%sf(i, j, k-1)/q_cons_vf(1)%sf(i, j, k-1))/(2._wp*dz(k))))

                    R_mu(2)%vf(2)%sf(i, j, k) = mu_visc * (2._wp*(q_cons_vf(momxb+1)%sf(i, j+1, k)/q_cons_vf(1)%sf(i, j+1, k) - q_cons_vf(momxb+1)%sf(i, j-1, k)/q_cons_vf(1)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                                - 2._wp/3._wp*((q_cons_vf(momxb)%sf(i+1, j, k)/q_cons_vf(1)%sf(i+1, j, k) - q_cons_vf(momxb)%sf(i-1, j, k)/q_cons_vf(1)%sf(i-1, j, k))/(2._wp*dx(i)) & 
                                                + (q_cons_vf(momxb+1)%sf(i, j+1, k)/q_cons_vf(1)%sf(i, j+1, k) - q_cons_vf(momxb+1)%sf(i, j-1, k)/q_cons_vf(1)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                                + (q_cons_vf(momxb+2)%sf(i, j, k+1)/q_cons_vf(1)%sf(i, j, k+1) - q_cons_vf(momxb+2)%sf(i, j, k-1)/q_cons_vf(1)%sf(i, j, k-1))/(2._wp*dz(k))))

                    R_mu(3)%vf(3)%sf(i, j, k) = mu_visc * (2._wp*(q_cons_vf(momxb+2)%sf(i, j, k+1)/q_cons_vf(1)%sf(i, j, k+1) - q_cons_vf(momxb+2)%sf(i, j, k-1)/q_cons_vf(1)%sf(i, j, k-1))/(2._wp*dz(k)) & 
                                                - 2._wp/3._wp*((q_cons_vf(momxb)%sf(i+1, j, k)/q_cons_vf(1)%sf(i+1, j, k) - q_cons_vf(momxb)%sf(i-1, j, k)/q_cons_vf(1)%sf(i-1, j, k))/(2._wp*dx(i)) & 
                                                + (q_cons_vf(momxb+1)%sf(i, j+1, k)/q_cons_vf(1)%sf(i, j+1, k) - q_cons_vf(momxb+1)%sf(i, j-1, k)/q_cons_vf(1)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                                + (q_cons_vf(momxb+2)%sf(i, j, k+1)/q_cons_vf(1)%sf(i, j, k+1) - q_cons_vf(momxb+2)%sf(i, j, k-1)/q_cons_vf(1)%sf(i, j, k-1))/(2._wp*dz(k))))

                    R_mu(1)%vf(2)%sf(i, j, k) = mu_visc * ((q_cons_vf(momxb)%sf(i, j+1, k)/q_cons_vf(1)%sf(i, j+1, k) - q_cons_vf(momxb)%sf(i, j-1, k)/q_cons_vf(1)%sf(i, j-1, k))/(2._wp*dy(j))/q_cons_vf(1)%sf(i, j, k) & 
                                                + (q_cons_vf(momxb+1)%sf(i+1, j, k)/q_cons_vf(1)%sf(i+1, j, k) - q_cons_vf(momxb+1)%sf(i-1, j, k)/q_cons_vf(1)%sf(i-1, j, k))/(2._wp*dx(i))/q_cons_vf(1)%sf(i, j, k))
                                            
                    R_mu(2)%vf(1)%sf(i, j, k) = R_mu(1)%vf(2)%sf(i, j, k)

                    R_mu(1)%vf(3)%sf(i, j, k) = mu_visc * ((q_cons_vf(momxb)%sf(i, j, k+1)/q_cons_vf(1)%sf(i, j, k+1) - q_cons_vf(momxb)%sf(i, j, k-1)/q_cons_vf(1)%sf(i, j, k-1))/(2._wp*dz(k))/q_cons_vf(1)%sf(i, j, k) & 
                                                + (q_cons_vf(momxb+2)%sf(i+1, j, k)/q_cons_vf(1)%sf(i+1, j, k) - q_cons_vf(momxb+2)%sf(i-1, j, k)/q_cons_vf(1)%sf(i-1, j, k))/(2._wp*dx(i))/q_cons_vf(1)%sf(i, j, k))

                    R_mu(3)%vf(1)%sf(i, j, k) = R_mu(1)%vf(3)%sf(i, j, k)

                    R_mu(2)%vf(3)%sf(i, j, k) = mu_visc * ((q_cons_vf(momxb+1)%sf(i, j, k+1)/q_cons_vf(1)%sf(i, j, k+1) - q_cons_vf(momxb+1)%sf(i, j, k-1)/q_cons_vf(1)%sf(i, j, k-1))/(2._wp*dz(k))/q_cons_vf(1)%sf(i, j, k) & 
                                                + (q_cons_vf(momxb+2)%sf(i, j+1, k)/q_cons_vf(1)%sf(i, j+1, k) - q_cons_vf(momxb+2)%sf(i, j-1, k)/q_cons_vf(1)%sf(i, j-1, k))/(2._wp*dy(j))/q_cons_vf(1)%sf(i, j, k))

                    R_mu(3)%vf(2)%sf(i, j, k) = R_mu(2)%vf(3)%sf(i, j, k)
                end do
            end do
        end do

    end subroutine s_setup_terms_filtering

    subroutine s_compute_pseudo_turbulent_reynolds_stress(q_cons_filtered, pt_Re_stress, mag_div_Ru)
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_filtered
        type(vector_field), dimension(1:num_dims), intent(inout) :: pt_Re_stress
        type(scalar_field), intent(inout) :: mag_div_Ru
        real(wp), dimension(1:num_dims, 0:m, 0:n, 0:p) :: div_Ru
        integer :: i, j, k, l, q    

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p
                    !$acc loop seq
                    do l = 1, num_dims
                        !$acc loop seq
                        do q = 1, num_dims
                            pt_Re_stress(l)%vf(q)%sf(i, j, k) = pt_Re_stress(l)%vf(q)%sf(i, j, k) &
                                                              - (q_cons_filtered(momxb-1+l)%sf(i, j, k) * q_cons_filtered(momxb-1+q)%sf(i, j, k) / q_cons_filtered(1)%sf(i, j, k))
                        end do
                    end do
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p  
                    !$acc loop seq
                    do l = 1, num_dims
                        !$acc loop seq
                        do q = 1, num_dims
                            pt_Re_stress(l)%vf(q)%sf(i, j, k) = pt_Re_stress(l)%vf(q)%sf(i, j, k) * q_cons_filtered(advxb)%sf(i, j, k)
                        end do 
                    end do 
                end do
            end do 
        end do

        ! set boundary buffer zone values
#ifdef MFC_MPI
        do l = 1, num_dims 
            do q = 1, num_dims
                call s_populate_scalarfield_buffers(pt_Re_stress(l)%vf(q))
            end do 
        end do
#else
        do l = 1, num_dims
            do q = 1, num_dims
                pt_Re_stress(l)%vf(q)%sf(-buff_size:-1, :, :) = pt_Re_stress(l)%vf(q)%sf(m-buff_size+1:m, :, :)
                pt_Re_stress(l)%vf(q)%sf(m+1:m+buff_size, :, :) = pt_Re_stress(l)%vf(q)%sf(0:buff_size-1, :, :)

                pt_Re_stress(l)%vf(q)%sf(:, -buff_size:-1, :) = pt_Re_stress(l)%vf(q)%sf(:, n-buff_size+1:n, :)
                pt_Re_stress(l)%vf(q)%sf(:, n+1:n+buff_size, :) = pt_Re_stress(l)%vf(q)%sf(:, 0:buff_size-1, :)

                pt_Re_stress(l)%vf(q)%sf(:, :, -buff_size:-1) = pt_Re_stress(l)%vf(q)%sf(:, :, p-buff_size+1:p)
                pt_Re_stress(l)%vf(q)%sf(:, :, p+1:p+buff_size) = pt_Re_stress(l)%vf(q)%sf(:, :, 0:buff_size-1)
            end do
        end do
#endif

        ! div(Ru), using CD2 FD scheme 
        !$acc parallel loop collapse(3) gang vector default(present) copy(div_Ru)
        do i = 0, m
            do j = 0, n 
                do k = 0, p
                    !$acc loop seq
                    do l = 1, num_dims
                        div_Ru(l, i, j, k) = (pt_Re_stress(l)%vf(1)%sf(i+1, j, k) - pt_Re_stress(l)%vf(1)%sf(i-1, j, k))/(2._wp*dx(i)) &
                                           + (pt_Re_stress(l)%vf(2)%sf(i, j+1, k) - pt_Re_stress(l)%vf(2)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                           + (pt_Re_stress(l)%vf(3)%sf(i, j, k+1) - pt_Re_stress(l)%vf(3)%sf(i, j, k-1))/(2._wp*dz(k))
                    end do
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present) copyin(div_Ru)
        do i = 0, m
            do j = 0, n
                do k = 0, p 
                    mag_div_Ru%sf(i, j, k) = sqrt(div_Ru(1, i, j, k)**2 + div_Ru(2, i, j, k)**2 + div_Ru(3, i, j, k)**2)
                end do
            end do
        end do

    end subroutine s_compute_pseudo_turbulent_reynolds_stress

    subroutine s_compute_R_mu(q_cons_filtered, R_mu, mag_div_R_mu)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_filtered
        type(vector_field), dimension(1:num_dims), intent(inout) :: R_mu
        type(scalar_field), intent(inout) :: mag_div_R_mu
        real(wp), dimension(1:num_dims, 0:m, 0:n, 0:p) :: div_R_mu

        integer :: i, j, k, l, q

        ! set buffers for filtered momentum quantities and density
#ifdef MFC_MPI
        do i = 1, momxe 
            call s_populate_scalarfield_buffers(q_cons_filtered(i))
        end do
#else
        do i = 1, momxe
            q_cons_filtered(i)%sf(-buff_size:-1, :, :) = q_cons_filtered(i)%sf(m-buff_size+1:m, :, :)
            q_cons_filtered(i)%sf(m+1:m+buff_size, :, :) = q_cons_filtered(i)%sf(0:buff_size-1, :, :)

            q_cons_filtered(i)%sf(:, -buff_size:-1, :) = q_cons_filtered(i)%sf(:, n-buff_size+1:n, :)
            q_cons_filtered(i)%sf(:, n+1:n+buff_size, :) = q_cons_filtered(i)%sf(:, 0:buff_size-1, :)

            q_cons_filtered(i)%sf(:, :, -buff_size:-1) = q_cons_filtered(i)%sf(:, :, p-buff_size+1:p)
            q_cons_filtered(i)%sf(:, :, p+1:p+buff_size) = q_cons_filtered(i)%sf(:, :, 0:buff_size-1)
        end do
#endif

        ! calculate R_mu
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    R_mu(1)%vf(1)%sf(i, j, k) = R_mu(1)%vf(1)%sf(i, j, k) - mu_visc * (2._wp*(q_cons_filtered(momxb)%sf(i+1, j, k)/q_cons_filtered(1)%sf(i+1, j, k) - q_cons_filtered(momxb)%sf(i-1, j, k)/q_cons_filtered(1)%sf(i-1, j, k))/(2._wp*dx(i)) & 
                                            - 2._wp/3._wp*((q_cons_filtered(momxb)%sf(i+1, j, k)/q_cons_filtered(1)%sf(i+1, j, k) - q_cons_filtered(momxb)%sf(i-1, j, k)/q_cons_filtered(1)%sf(i-1, j, k))/(2._wp*dx(i)) & 
                                            + (q_cons_filtered(momxb+1)%sf(i, j+1, k)/q_cons_filtered(1)%sf(i, j+1, k) - q_cons_filtered(momxb+1)%sf(i, j-1, k)/q_cons_filtered(1)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                            + (q_cons_filtered(momxb+2)%sf(i, j, k+1)/q_cons_filtered(1)%sf(i, j, k+1) - q_cons_filtered(momxb+2)%sf(i, j, k-1)/q_cons_filtered(1)%sf(i, j, k-1))/(2._wp*dz(k))))

                    R_mu(2)%vf(2)%sf(i, j, k) = R_mu(2)%vf(2)%sf(i, j, k) - mu_visc * (2._wp*(q_cons_filtered(momxb+1)%sf(i, j+1, k)/q_cons_filtered(1)%sf(i, j+1, k) - q_cons_filtered(momxb+1)%sf(i, j-1, k)/q_cons_filtered(1)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                            - 2._wp/3._wp*((q_cons_filtered(momxb)%sf(i+1, j, k)/q_cons_filtered(1)%sf(i+1, j, k) - q_cons_filtered(momxb)%sf(i-1, j, k)/q_cons_filtered(1)%sf(i-1, j, k))/(2._wp*dx(i)) & 
                                            + (q_cons_filtered(momxb+1)%sf(i, j+1, k)/q_cons_filtered(1)%sf(i, j+1, k) - q_cons_filtered(momxb+1)%sf(i, j-1, k)/q_cons_filtered(1)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                            + (q_cons_filtered(momxb+2)%sf(i, j, k+1)/q_cons_filtered(1)%sf(i, j, k+1) - q_cons_filtered(momxb+2)%sf(i, j, k-1)/q_cons_filtered(1)%sf(i, j, k-1))/(2._wp*dz(k))))

                    R_mu(3)%vf(3)%sf(i, j, k) = R_mu(3)%vf(3)%sf(i, j, k) - mu_visc * (2._wp*(q_cons_filtered(momxb+2)%sf(i, j, k+1)/q_cons_filtered(1)%sf(i, j, k+1) - q_cons_filtered(momxb+2)%sf(i, j, k-1)/q_cons_filtered(1)%sf(i, j, k-1))/(2._wp*dz(k)) & 
                                            - 2._wp/3._wp*((q_cons_filtered(momxb)%sf(i+1, j, k)/q_cons_filtered(1)%sf(i+1, j, k) - q_cons_filtered(momxb)%sf(i-1, j, k)/q_cons_filtered(1)%sf(i-1, j, k))/(2._wp*dx(i)) & 
                                            + (q_cons_filtered(momxb+1)%sf(i, j+1, k)/q_cons_filtered(1)%sf(i, j+1, k) - q_cons_filtered(momxb+1)%sf(i, j-1, k)/q_cons_filtered(1)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                            + (q_cons_filtered(momxb+2)%sf(i, j, k+1)/q_cons_filtered(1)%sf(i, j, k+1) - q_cons_filtered(momxb+2)%sf(i, j, k-1)/q_cons_filtered(1)%sf(i, j, k-1))/(2._wp*dz(k))))

                    R_mu(1)%vf(2)%sf(i, j, k) = R_mu(1)%vf(2)%sf(i, j, k) - mu_visc * ((q_cons_filtered(momxb)%sf(i, j+1, k)/q_cons_filtered(1)%sf(i, j+1, k) - q_cons_filtered(momxb)%sf(i, j-1, k)/q_cons_filtered(1)%sf(i, j-1, k))/(2._wp*dy(j))/q_cons_filtered(1)%sf(i, j, k) & 
                                            + (q_cons_filtered(momxb+1)%sf(i+1, j, k)/q_cons_filtered(1)%sf(i+1, j, k) - q_cons_filtered(momxb+1)%sf(i-1, j, k)/q_cons_filtered(1)%sf(i-1, j, k))/(2._wp*dx(i))/q_cons_filtered(1)%sf(i, j, k))
                                        
                    R_mu(2)%vf(1)%sf(i, j, k) = R_mu(1)%vf(2)%sf(i, j, k)

                    R_mu(1)%vf(3)%sf(i, j, k) = R_mu(1)%vf(3)%sf(i, j, k) - mu_visc * ((q_cons_filtered(momxb)%sf(i, j, k+1)/q_cons_filtered(1)%sf(i, j, k+1) - q_cons_filtered(momxb)%sf(i, j, k-1)/q_cons_filtered(1)%sf(i, j, k-1))/(2._wp*dz(k))/q_cons_filtered(1)%sf(i, j, k) & 
                                            + (q_cons_filtered(momxb+2)%sf(i+1, j, k)/q_cons_filtered(1)%sf(i+1, j, k) - q_cons_filtered(momxb+2)%sf(i-1, j, k)/q_cons_filtered(1)%sf(i-1, j, k))/(2._wp*dx(i))/q_cons_filtered(1)%sf(i, j, k))

                    R_mu(3)%vf(1)%sf(i, j, k) = R_mu(1)%vf(3)%sf(i, j, k)

                    R_mu(2)%vf(3)%sf(i, j, k) = R_mu(2)%vf(3)%sf(i, j, k) - mu_visc * ((q_cons_filtered(momxb+1)%sf(i, j, k+1)/q_cons_filtered(1)%sf(i, j, k+1) - q_cons_filtered(momxb+1)%sf(i, j, k-1)/q_cons_filtered(1)%sf(i, j, k-1))/(2._wp*dz(k))/q_cons_filtered(1)%sf(i, j, k) & 
                                            + (q_cons_filtered(momxb+2)%sf(i, j+1, k)/q_cons_filtered(1)%sf(i, j+1, k) - q_cons_filtered(momxb+2)%sf(i, j-1, k)/q_cons_filtered(1)%sf(i, j-1, k))/(2._wp*dy(j))/q_cons_filtered(1)%sf(i, j, k))

                    R_mu(3)%vf(2)%sf(i, j, k) = R_mu(2)%vf(3)%sf(i, j, k)
                    
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p 
                    !$acc loop seq
                    do l = 1, num_dims
                        !$acc loop seq
                        do q = 1, num_dims
                            R_mu(l)%vf(q)%sf(i, j, k) = R_mu(l)%vf(q)%sf(i, j, k) * q_cons_filtered(advxb)%sf(i, j, k)
                        end do 
                    end do 
                end do
            end do 
        end do

        ! set boundary buffer zone values
#ifdef MFC_MPI
        do l = 1, num_dims
            do q = 1, num_dims
                call s_populate_scalarfield_buffers(R_mu(l)%vf(q))
            end do
        end do
#else
        do l = 1, num_dims
            do q = 1, num_dims
                R_mu(l)%vf(q)%sf(-buff_size:-1, :, :) = R_mu(l)%vf(q)%sf(m-buff_size+1:m, :, :)
                R_mu(l)%vf(q)%sf(m+1:m+buff_size, :, :) = R_mu(l)%vf(q)%sf(0:buff_size-1, :, :)

                R_mu(l)%vf(q)%sf(:, -buff_size:-1, :) = R_mu(l)%vf(q)%sf(:, n-buff_size+1:n, :)
                R_mu(l)%vf(q)%sf(:, n+1:n+buff_size, :) = R_mu(l)%vf(q)%sf(:, 0:buff_size-1, :)

                R_mu(l)%vf(q)%sf(:, :, -buff_size:-1) = R_mu(l)%vf(q)%sf(:, :, p-buff_size+1:p)
                R_mu(l)%vf(q)%sf(:, :, p+1:p+buff_size) = R_mu(l)%vf(q)%sf(:, :, 0:buff_size-1)
            end do
        end do
#endif

        ! div(R_mu), using CD2 FD scheme 
        !$acc parallel loop collapse(3) gang vector default(present) copy(div_R_mu)
        do i = 0, m
            do j = 0, n 
                do k = 0, p
                    !$acc loop seq
                    do l = 1, num_dims
                        div_R_mu(l, i, j, k) = (R_mu(l)%vf(1)%sf(i+1, j, k) - R_mu(l)%vf(1)%sf(i-1, j, k))/(2._wp*dx(i)) &
                                             + (R_mu(l)%vf(2)%sf(i, j+1, k) - R_mu(l)%vf(2)%sf(i, j-1, k))/(2._wp*dy(j)) & 
                                             + (R_mu(l)%vf(3)%sf(i, j, k+1) - R_mu(l)%vf(3)%sf(i, j, k-1))/(2._wp*dz(k))
                    end do
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present) copyin(div_R_mu)
        do i = 0, m
            do j = 0, n
                do k = 0, p 
                    mag_div_R_mu%sf(i, j, k) = sqrt(div_R_mu(1, i, j, k)**2 + div_R_mu(2, i, j, k)**2 + div_R_mu(3, i, j, k)**2)
                end do
            end do
        end do

    end subroutine s_compute_R_mu

    subroutine s_compute_interphase_momentum_exchange_term(pres_visc_stress_filtered, mag_F_IMET)
        type(scalar_field), dimension(1:num_dims), intent(in) :: pres_visc_stress_filtered
        type(scalar_field), intent(inout) :: mag_F_IMET

        integer :: i, j, k, l, q, ii

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p 
                    mag_F_IMET%sf(i, j, k) = sqrt(pres_visc_stress_filtered(1)%sf(i, j, k)**2 & 
                                                + pres_visc_stress_filtered(2)%sf(i, j, k)**2 & 
                                                + pres_visc_stress_filtered(3)%sf(i, j, k)**2)
                end do
            end do
        end do 

    end subroutine s_compute_interphase_momentum_exchange_term

    subroutine s_finalize_fftw_explicit_filter_module
        @:DEALLOCATE(fluid_indicator_function_I%sf)

        @:DEALLOCATE(data_real_in1d, data_cmplx_out1d, data_cmplx_out1dy)
        @:DEALLOCATE(cmplx_kernelG1d, real_kernelG_in)
        @:DEALLOCATE(data_real_3D_slabz, data_cmplx_slabz, data_cmplx_slaby)

#if defined(MFC_OpenACC)
        ierr = cufftDestroy(plan_x_fwd_gpu)
        ierr = cufftDestroy(plan_x_bwd_gpu) 
        ierr = cufftDestroy(plan_y_gpu)
        ierr = cufftDestroy(plan_z_gpu)
#else
        call fftw_destroy_plan(plan_x_r2c_fwd)
        call fftw_destroy_plan(plan_x_c2r_bwd)
        call fftw_destroy_plan(plan_y_c2c_fwd) 
        call fftw_destroy_plan(plan_y_c2c_bwd) 
        call fftw_destroy_plan(plan_z_c2c_fwd) 
        call fftw_destroy_plan(plan_z_c2c_bwd) 
        call fftw_destroy_plan(plan_x_r2c_kernelG)
        call fftw_destroy_plan(plan_y_c2c_kernelG)
        call fftw_destroy_plan(plan_z_c2c_kernelG)
#endif

    end subroutine s_finalize_fftw_explicit_filter_module

end module m_volume_filtering