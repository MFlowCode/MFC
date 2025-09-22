#:include 'macros.fpp'

module m_volume_filtering

    use, intrinsic :: iso_c_binding

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_ibm

    use m_boundary_common

    use m_nvtx

#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

#if defined(MFC_OpenACC) && defined(__PGI)
    use cufft
#endif

    implicit none

    private; public :: s_initialize_fftw_explicit_filter_module, &
 s_initialize_filtering_kernel, s_initialize_fluid_indicator_function, & 
 s_initialize_filtered_fluid_indicator_function, s_initialize_fluid_indicator_gradient, &
 s_finalize_fftw_explicit_filter_module, s_volume_filter_momentum_eqn, s_apply_fftw_filter_scalarfield, s_filter_tensor_field, &
 s_compute_viscous_stress_tensor, s_compute_stress_tensor, s_compute_divergence_stress_tensor, s_compute_particle_forces, &
 s_mpi_transpose_slabZ2Y, s_mpi_transpose_slabY2Z, s_mpi_transpose_slabZ2Y_tensor, s_mpi_transpose_slabY2Z_tensor, s_mpi_FFT_fwd, s_mpi_FFT_bwd, &
 s_setup_terms_filtering, s_compute_pseudo_turbulent_reynolds_stress, s_compute_effective_viscosity

#if !defined(MFC_OpenACC)
    include 'fftw3.f03'
#endif

    integer :: ierr   

    ! fluid indicator function (1 = fluid, 0 = otherwise)
    type(scalar_field), public :: fluid_indicator_function
    type(scalar_field), public :: filtered_fluid_indicator_function
    type(scalar_field), allocatable, dimension(:) :: grad_fluid_indicator

    ! volume filtered conservative variables
    type(scalar_field), allocatable, dimension(:), public :: q_cons_filtered
    type(scalar_field), public :: filtered_pressure

    ! viscous and pressure+viscous stress tensors
    type(vector_field), allocatable, dimension(:) :: visc_stress
    type(vector_field), allocatable, dimension(:) :: pres_visc_stress

    ! divergence of stress tensor
    type(scalar_field), allocatable, dimension(:) :: div_pres_visc_stress
    
    ! unclosed terms in volume filtered momentum equation
    type(vector_field), allocatable, dimension(:), public :: reynolds_stress
    type(vector_field), allocatable, dimension(:), public :: eff_visc
    type(scalar_field), allocatable, dimension(:), public :: int_mom_exch

    ! 1/mu
    real(wp), allocatable, dimension(:, :) :: Res

    ! x-,y-,z-direction forces on particles
    real(wp), allocatable, dimension(:, :) :: particle_forces

    !$acc declare create(fluid_indicator_function, filtered_fluid_indicator_function, grad_fluid_indicator)
    !$acc declare create(q_cons_filtered, filtered_pressure)
    !$acc declare create(visc_stress, pres_visc_stress, div_pres_visc_stress)
    !$acc declare create(reynolds_stress, eff_visc, int_mom_exch)
    !$acc declare create(Res, particle_forces)

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
    ! 3D arrays for slab transposes of tensor quantities
    complex(c_double_complex), allocatable :: data_cmplx_slabz_tensor(:, :, :, :), data_cmplx_slaby_tensor(:, :, :, :)

    ! input/output array for FFT routine
    real(c_double), allocatable :: data_real_3D_slabz(:, :, :)

    ! filtering kernel in physical space
    real(c_double), allocatable :: real_kernelG_in(:, :, :)

    ! FFT of filtering kernel
    complex(c_double_complex), allocatable :: cmplx_kernelG1d(:)

    !$acc declare create(Nx, Ny, Nz, NxC, Nyloc, Nzloc)
    !$acc declare create(data_real_in1d, data_cmplx_out1d, data_cmplx_out1dy)
    !$acc declare create(data_cmplx_slabz, data_cmplx_slaby, data_cmplx_slabz_tensor, data_cmplx_slaby_tensor, data_real_3D_slabz, real_kernelG_in, cmplx_kernelG1d)

    ! buffers for data transpose
    complex(c_double_complex), allocatable :: sendbuf_sf(:), recvbuf_sf(:)
    complex(c_double_complex), allocatable :: sendbuf_tensor(:), recvbuf_tensor(:)

contains

    !< create fft plans to be used for explicit filtering of data 
    subroutine s_initialize_fftw_explicit_filter_module
        integer :: i, j, k
        integer :: size_n(1), inembed(1), onembed(1)
        
        @:ALLOCATE(q_cons_filtered(1:sys_size-1))
        do i = 1, sys_size-1
            @:ALLOCATE(q_cons_filtered(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_cons_filtered(i))
        end do

        @:ALLOCATE(filtered_pressure%sf(idwbuff(1)%beg:idwbuff(1)%end, &
            idwbuff(2)%beg:idwbuff(2)%end, &
            idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(filtered_pressure)

        @:ALLOCATE(visc_stress(1:num_dims))
        do i = 1, num_dims
            @:ALLOCATE(visc_stress(i)%vf(1:num_dims))
        end do
        do i = 1, num_dims
            do j = 1, num_dims 
                @:ALLOCATE(visc_stress(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
            end do 
            @:ACC_SETUP_VFs(visc_stress(i))
        end do

        @:ALLOCATE(pres_visc_stress(1:num_dims))
        do i = 1, num_dims
            @:ALLOCATE(pres_visc_stress(i)%vf(1:num_dims))
        end do
        do i = 1, num_dims
            do j = 1, num_dims 
                @:ALLOCATE(pres_visc_stress(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
            end do 
            @:ACC_SETUP_VFs(pres_visc_stress(i))
        end do

        @:ALLOCATE(div_pres_visc_stress(1:num_dims))
        do i = 1, num_dims
            @:ALLOCATE(div_pres_visc_stress(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(div_pres_visc_stress(i))
        end do

        @:ALLOCATE(reynolds_stress(1:num_dims))
        do i = 1, num_dims
            @:ALLOCATE(reynolds_stress(i)%vf(1:num_dims))
        end do
        do i = 1, num_dims
            do j = 1, num_dims
                @:ALLOCATE(reynolds_stress(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
            end do
            @:ACC_SETUP_VFs(reynolds_stress(i))
        end do

        @:ALLOCATE(eff_visc(1:num_dims))
        do i = 1, num_dims
            @:ALLOCATE(eff_visc(i)%vf(1:num_dims))
        end do
        do i = 1, num_dims
            do j = 1, num_dims
                @:ALLOCATE(eff_visc(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
            end do
            @:ACC_SETUP_VFs(eff_visc(i))
        end do

        @:ALLOCATE(int_mom_exch(1:num_dims))
        do i = 1, num_dims
            @:ALLOCATE(int_mom_exch(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(int_mom_exch(i))
        end do

        if (viscous) then
            @:ALLOCATE(Res(1:2, 1:maxval(Re_size)))
        end if

        if (viscous) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            !$acc update device(Res)
        end if

        @:ALLOCATE(particle_forces(0:num_ibs, 3))

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
        @:ALLOCATE(data_cmplx_slabz_tensor(9, NxC, Ny, Nzloc))
        @:ALLOCATE(data_cmplx_slaby_tensor(9, NxC, Nyloc, Nz))

        allocate(sendbuf_sf(NxC*Nyloc*Nzloc*num_procs))
        allocate(recvbuf_sf(NxC*Nyloc*Nzloc*num_procs))
        allocate(sendbuf_tensor(9*NxC*Nyloc*Nzloc*num_procs))
        allocate(recvbuf_tensor(9*NxC*Nyloc*Nzloc*num_procs))

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

        ! file for particle forces
        if (proc_rank == 0) then
            open(unit=100, file='particle_force.bin', status='replace', form='unformatted', access='stream', action='write')
        end if

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
        sigma_stddev = filter_width

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

    !< initialize fluid indicator function and filtered fluid indicator function
    subroutine s_initialize_fluid_indicator_function 
        integer :: i, j, k 

        @:ALLOCATE(fluid_indicator_function%sf(-1:m+1, -1:n+1, -1:p+1))
        @:ACC_SETUP_SFs(fluid_indicator_function)

        ! define fluid indicator function
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = -1, m+1
            do j = -1, n+1 
                do k = -1, p+1
                    if (ib_markers%sf(i, j, k) == 0) then 
                        fluid_indicator_function%sf(i, j, k) = 1.0_dp
                    else 
                        fluid_indicator_function%sf(i, j, k) = 0.0_dp
                    end if
                end do
            end do
        end do
    
    end subroutine s_initialize_fluid_indicator_function

    subroutine s_initialize_filtered_fluid_indicator_function
        integer :: i, j, k
        
        @:ALLOCATE(filtered_fluid_indicator_function%sf(0:m, 0:n, 0:p))
        @:ACC_SETUP_SFs(filtered_fluid_indicator_function)

        ! filter fluid indicator function 
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc 
                    data_real_3D_slabz(i, j, k) = fluid_indicator_function%sf(i-1, j-1, k-1)
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


    subroutine s_initialize_fluid_indicator_gradient
        integer :: i, j, k

        @:ALLOCATE(grad_fluid_indicator(1:3))
        do i = 1, 3
            @:ALLOCATE(grad_fluid_indicator(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(grad_fluid_indicator(i))
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p 
                    grad_fluid_indicator(1)%sf(i, j, k) = (fluid_indicator_function%sf(i+1, j, k) - &
                                                           fluid_indicator_function%sf(i-1, j, k)) / & 
                                                           (x_cc(i+1) - x_cc(i-1))
                    grad_fluid_indicator(2)%sf(i, j, k) = (fluid_indicator_function%sf(i, j+1, k) - &
                                                           fluid_indicator_function%sf(i, j-1, k)) / & 
                                                           (y_cc(j+1) - y_cc(j-1))
                    grad_fluid_indicator(3)%sf(i, j, k) = (fluid_indicator_function%sf(i, j, k+1) - &
                                                           fluid_indicator_function%sf(i, j, k-1)) / & 
                                                           (z_cc(k+1) - z_cc(k-1))
                end do 
            end do 
        end do

    end subroutine s_initialize_fluid_indicator_gradient


    !< calculate the unclosed terms present in the volume filtered momentum equation
    subroutine s_volume_filter_momentum_eqn(q_cons_vf, q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k

        call nvtxStartRange("FILTER-CONS-VARS")
        do i = 1, sys_size-1
            call s_apply_fftw_filter_scalarfield(filtered_fluid_indicator_function, .true., q_cons_vf(i), q_cons_filtered(i))
        end do 
        call s_apply_fftw_filter_scalarfield(filtered_fluid_indicator_function, .true., q_prim_vf(E_idx), filtered_pressure)
        call nvtxEndRange

        call nvtxStartRange("COMPUTE-UNCLOSED-TERMS")
        call s_setup_terms_filtering(q_cons_vf, q_prim_vf, reynolds_stress, visc_stress, pres_visc_stress, div_pres_visc_stress)

        ! pseudo turbulent reynolds stress
        ! do i = 1, num_dims 
        !     do j = 1, num_dims
        !         call s_apply_fftw_filter_scalarfield(filtered_fluid_indicator_function, .true., reynolds_stress(i)%vf(j))
        !     end do
        ! end do 
        call s_filter_tensor_field(reynolds_stress)
        ! effective viscosity
        ! do i = 1, num_dims 
        !     do j = 1, num_dims
        !         call s_apply_fftw_filter_scalarfield(filtered_fluid_indicator_function, .true., visc_stress(i)%vf(j), eff_visc(i)%vf(j))
        !     end do
        ! end do 
        call s_filter_tensor_field(visc_stress, eff_visc)
        ! interphase momentum exchange
        call s_compute_interphase_momentum_exchange(filtered_fluid_indicator_function, grad_fluid_indicator, pres_visc_stress, int_mom_exch)

        call s_compute_pseudo_turbulent_reynolds_stress(q_cons_filtered, reynolds_stress)
        call s_compute_effective_viscosity(q_cons_filtered, eff_visc, visc_stress)
        call nvtxEndRange

    end subroutine s_volume_filter_momentum_eqn

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
                        data_real_3D_slabz(i+1, j+1, k+1) = q_temp_in%sf(i, j, k) * fluid_indicator_function%sf(i, j, k)
                    end do 
                end do 
            end do
        else 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p 
                        data_real_3D_slabz(i+1, j+1, k+1) = q_temp_in%sf(i, j, k) * (1.0_dp - fluid_indicator_function%sf(i, j, k))
                    end do 
                end do 
            end do
        end if

        call nvtxStartRange("FORWARD-3D-FFT")
        call s_mpi_FFT_fwd 
        call nvtxEndRange

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz 
                    data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) * cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                end do 
            end do 
        end do

        call nvtxStartRange("BACKWARD-3D-FFT")
        call s_mpi_FFT_bwd
        call nvtxEndRange

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

    ! compute viscous stress tensor
    subroutine s_compute_viscous_stress_tensor(visc_stress, q_prim_vf, q_cons_filtered)
        type(vector_field), dimension(num_dims), intent(inout) :: visc_stress 
        type(scalar_field), dimension(sys_size), intent(in), optional :: q_prim_vf
        type(scalar_field), dimension(sys_size-1), intent(in), optional :: q_cons_filtered
        real(wp) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz ! spatial velocity derivatives
        integer :: i, j, k 

        if (present(q_prim_vf)) then
            !$acc parallel loop collapse(3) gang vector default(present) private(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p
                        ! velocity gradients, local to each process
                        dudx = ( q_prim_vf(2)%sf(i+1, j, k) - q_prim_vf(2)%sf(i-1, j, k) ) / (dx(i-1) + dx(i+1))
                        dudy = ( q_prim_vf(2)%sf(i, j+1, k) - q_prim_vf(2)%sf(i, j-1, k) ) / (dy(j-1) + dy(j+1))
                        dudz = ( q_prim_vf(2)%sf(i, j, k+1) - q_prim_vf(2)%sf(i, j, k-1) ) / (dz(k-1) + dz(k+1))

                        dvdx = ( q_prim_vf(3)%sf(i+1, j, k) - q_prim_vf(3)%sf(i-1, j, k) ) / (dx(i-1) + dx(i+1))
                        dvdy = ( q_prim_vf(3)%sf(i, j+1, k) - q_prim_vf(3)%sf(i, j-1, k) ) / (dy(j-1) + dy(j+1))
                        dvdz = ( q_prim_vf(3)%sf(i, j, k+1) - q_prim_vf(3)%sf(i, j, k-1) ) / (dz(k-1) + dz(k+1))

                        dwdx = ( q_prim_vf(4)%sf(i+1, j, k) - q_prim_vf(4)%sf(i-1, j, k) ) / (dx(i-1) + dx(i+1))
                        dwdy = ( q_prim_vf(4)%sf(i, j+1, k) - q_prim_vf(4)%sf(i, j-1, k) ) / (dy(j-1) + dy(j+1))
                        dwdz = ( q_prim_vf(4)%sf(i, j, k+1) - q_prim_vf(4)%sf(i, j, k-1) ) / (dz(k-1) + dz(k+1))

                        ! viscous stress tensor, visc_stress(row, column)
                        visc_stress(1)%vf(1)%sf(i, j, k) = (4._wp/3._wp * dudx - 2._wp/3._wp * (dvdy + dwdz)) / Res(1, 1)
                        visc_stress(1)%vf(2)%sf(i, j, k) = (dudy + dvdx) / Res(1, 1)
                        visc_stress(1)%vf(3)%sf(i, j, k) = (dudz + dwdx) / Res(1, 1)
                        visc_stress(2)%vf(1)%sf(i, j, k) = (dvdx + dudy) / Res(1, 1)
                        visc_stress(2)%vf(2)%sf(i, j, k) = (4._wp/3._wp * dvdy - 2._wp/3._wp * (dudx + dwdz)) / Res(1, 1)
                        visc_stress(2)%vf(3)%sf(i, j, k) = (dvdz + dwdy) / Res(1, 1)
                        visc_stress(3)%vf(1)%sf(i, j, k) = (dwdx + dudz) / Res(1, 1)
                        visc_stress(3)%vf(2)%sf(i, j, k) = (dwdy + dvdz) / Res(1, 1)
                        visc_stress(3)%vf(3)%sf(i, j, k) = (4._wp/3._wp * dwdz - 2._wp/3._wp * (dudx + dvdy)) / Res(1, 1)
                    end do 
                end do 
            end do
        else if (present(q_cons_filtered)) then
            !$acc parallel loop collapse(3) gang vector default(present) private(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p
                        ! velocity gradients, local to each process
                        dudx = ( q_cons_filtered(2)%sf(i+1, j, k)/q_cons_filtered(1)%sf(i+1, j, k) - q_cons_filtered(2)%sf(i-1, j, k)/q_cons_filtered(1)%sf(i-1, j, k) ) / (dx(i-1) + dx(i+1))
                        dudy = ( q_cons_filtered(2)%sf(i, j+1, k)/q_cons_filtered(1)%sf(i, j+1, k) - q_cons_filtered(2)%sf(i, j-1, k)/q_cons_filtered(1)%sf(i, j-1, k) ) / (dy(j-1) + dy(j+1))
                        dudz = ( q_cons_filtered(2)%sf(i, j, k+1)/q_cons_filtered(1)%sf(i, j, k+1) - q_cons_filtered(2)%sf(i, j, k-1)/q_cons_filtered(1)%sf(i, j, k-1) ) / (dz(k-1) + dz(k+1))

                        dvdx = ( q_cons_filtered(3)%sf(i+1, j, k)/q_cons_filtered(1)%sf(i+1, j, k) - q_cons_filtered(3)%sf(i-1, j, k)/q_cons_filtered(1)%sf(i-1, j, k) ) / (dx(i-1) + dx(i+1))
                        dvdy = ( q_cons_filtered(3)%sf(i, j+1, k)/q_cons_filtered(1)%sf(i, j+1, k) - q_cons_filtered(3)%sf(i, j-1, k)/q_cons_filtered(1)%sf(i, j-1, k) ) / (dy(j-1) + dy(j+1))
                        dvdz = ( q_cons_filtered(3)%sf(i, j, k+1)/q_cons_filtered(1)%sf(i, j, k+1) - q_cons_filtered(3)%sf(i, j, k-1)/q_cons_filtered(1)%sf(i, j, k-1) ) / (dz(k-1) + dz(k+1))

                        dwdx = ( q_cons_filtered(4)%sf(i+1, j, k)/q_cons_filtered(1)%sf(i+1, j, k) - q_cons_filtered(4)%sf(i-1, j, k)/q_cons_filtered(1)%sf(i-1, j, k) ) / (dx(i-1) + dx(i+1))
                        dwdy = ( q_cons_filtered(4)%sf(i, j+1, k)/q_cons_filtered(1)%sf(i, j+1, k) - q_cons_filtered(4)%sf(i, j-1, k)/q_cons_filtered(1)%sf(i, j-1, k) ) / (dy(j-1) + dy(j+1))
                        dwdz = ( q_cons_filtered(4)%sf(i, j, k+1)/q_cons_filtered(1)%sf(i, j, k+1) - q_cons_filtered(4)%sf(i, j, k-1)/q_cons_filtered(1)%sf(i, j, k-1) ) / (dz(k-1) + dz(k+1))

                        ! viscous stress tensor, visc_stress(row, column)
                        visc_stress(1)%vf(1)%sf(i, j, k) = (4._wp/3._wp * dudx - 2._wp/3._wp * (dvdy + dwdz)) / Res(1, 1)
                        visc_stress(1)%vf(2)%sf(i, j, k) = (dudy + dvdx) / Res(1, 1)
                        visc_stress(1)%vf(3)%sf(i, j, k) = (dudz + dwdx) / Res(1, 1)
                        visc_stress(2)%vf(1)%sf(i, j, k) = (dvdx + dudy) / Res(1, 1)
                        visc_stress(2)%vf(2)%sf(i, j, k) = (4._wp/3._wp * dvdy - 2._wp/3._wp * (dudx + dwdz)) / Res(1, 1)
                        visc_stress(2)%vf(3)%sf(i, j, k) = (dvdz + dwdy) / Res(1, 1)
                        visc_stress(3)%vf(1)%sf(i, j, k) = (dwdx + dudz) / Res(1, 1)
                        visc_stress(3)%vf(2)%sf(i, j, k) = (dwdy + dvdz) / Res(1, 1)
                        visc_stress(3)%vf(3)%sf(i, j, k) = (4._wp/3._wp * dwdz - 2._wp/3._wp * (dudx + dvdy)) / Res(1, 1)
                    end do 
                end do 
            end do
        end if

    end subroutine s_compute_viscous_stress_tensor
    
    subroutine s_compute_stress_tensor(pres_visc_stress, visc_stress, q_cons_vf, q_prim_vf)
        type(vector_field), dimension(num_dims), intent(inout) :: pres_visc_stress
        type(vector_field), dimension(num_dims), intent(in) :: visc_stress
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp) :: pressure
        integer :: i, j, k

        !$acc parallel loop collapse(3) gang vector default(present) private(pressure)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p
                    pres_visc_stress(1)%vf(1)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k) - visc_stress(1)%vf(1)%sf(i, j, k)
                    pres_visc_stress(1)%vf(2)%sf(i, j, k) = - visc_stress(1)%vf(2)%sf(i, j, k) 
                    pres_visc_stress(1)%vf(3)%sf(i, j, k) = - visc_stress(1)%vf(3)%sf(i, j, k)
                    pres_visc_stress(2)%vf(1)%sf(i, j, k) = - visc_stress(2)%vf(1)%sf(i, j, k)
                    pres_visc_stress(2)%vf(2)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k) - visc_stress(2)%vf(2)%sf(i, j, k) 
                    pres_visc_stress(2)%vf(3)%sf(i, j, k) = - visc_stress(2)%vf(3)%sf(i, j, k)
                    pres_visc_stress(3)%vf(1)%sf(i, j, k) = - visc_stress(3)%vf(1)%sf(i, j, k)
                    pres_visc_stress(3)%vf(2)%sf(i, j, k) = - visc_stress(3)%vf(2)%sf(i, j, k)
                    pres_visc_stress(3)%vf(3)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k) - visc_stress(3)%vf(3)%sf(i, j, k)
                end do 
            end do 
        end do 

    end subroutine s_compute_stress_tensor

    !< compute the divergence of the pressure-viscous stress tensor
    subroutine s_compute_divergence_stress_tensor(div_stress_tensor, stress_tensor)
        type(scalar_field), dimension(num_dims), intent(inout) :: div_stress_tensor
        type(vector_field), dimension(num_dims), intent(in) :: stress_tensor
        integer :: i, j, k

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p
                    div_stress_tensor(1)%sf(i, j, k) = (stress_tensor(1)%vf(1)%sf(i+1, j, k) - stress_tensor(1)%vf(1)%sf(i-1, j, k)) / (dx(i-1) + dx(i+1)) &
                                                     + (stress_tensor(2)%vf(1)%sf(i, j+1, k) - stress_tensor(2)%vf(1)%sf(i, j-1, k)) / (dy(j-1) + dy(j+1)) &
                                                     + (stress_tensor(3)%vf(1)%sf(i, j, k+1) - stress_tensor(3)%vf(1)%sf(i, j, k-1)) / (dz(k-1) + dz(k+1))

                    div_stress_tensor(2)%sf(i, j, k) = (stress_tensor(1)%vf(2)%sf(i+1, j, k) - stress_tensor(1)%vf(2)%sf(i-1, j, k)) / (dx(i-1) + dx(i+1)) & 
                                                     + (stress_tensor(2)%vf(2)%sf(i, j+1, k) - stress_tensor(2)%vf(2)%sf(i, j-1, k)) / (dy(j-1) + dy(j+1)) & 
                                                     + (stress_tensor(3)%vf(2)%sf(i, j, k+1) - stress_tensor(3)%vf(2)%sf(i, j, k-1)) / (dz(k-1) + dz(k+1))

                    div_stress_tensor(3)%sf(i, j, k) = (stress_tensor(1)%vf(3)%sf(i+1, j, k) - stress_tensor(1)%vf(3)%sf(i-1, j, k)) / (dx(i-1) + dx(i+1)) & 
                                                     + (stress_tensor(2)%vf(3)%sf(i, j+1, k) - stress_tensor(2)%vf(3)%sf(i, j-1, k)) / (dy(j-1) + dy(j+1)) & 
                                                     + (stress_tensor(3)%vf(3)%sf(i, j, k+1) - stress_tensor(3)%vf(3)%sf(i, j ,k-1)) / (dz(k-1) + dz(k+1))
                end do 
            end do 
        end do

    end subroutine s_compute_divergence_stress_tensor

    !< setup for calculation of unclosed terms in volume filtered momentum eqn
    subroutine s_setup_terms_filtering(q_cons_vf, q_prim_vf, reynolds_stress, visc_stress, pres_visc_stress, div_pres_visc_stress)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(vector_field), dimension(1:num_dims), intent(inout) :: reynolds_stress
        type(vector_field), dimension(1:num_dims), intent(inout) :: visc_stress
        type(vector_field), dimension(1:num_dims), intent(inout) :: pres_visc_stress
        type(scalar_field), dimension(1:num_dims), intent(inout) :: div_pres_visc_stress

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
                            reynolds_stress(l)%vf(q)%sf(i, j, k) = q_cons_vf(1)%sf(i, j, k) * (q_prim_vf(momxb-1+l)%sf(i, j, k) * q_prim_vf(momxb-1+q)%sf(i, j, k)) ! rho*(u x u) 
                        end do
                    end do
                end do
            end do 
        end do

        ! set density and momentum buffers
#ifdef MFC_MPI
        do i = contxb, momxe 
            call s_populate_scalarfield_buffers(q_cons_vf(i))
        end do
#else
        do i = contxb, momxe
            q_cons_vf(i)%sf(-buff_size:-1, :, :) = q_cons_vf(i)%sf(m-buff_size+1:m, :, :)
            q_cons_vf(i)%sf(m+1:m+buff_size, :, :) = q_cons_vf(i)%sf(0:buff_size-1, :, :)

            q_cons_vf(i)%sf(:, -buff_size:-1, :) = q_cons_vf(i)%sf(:, n-buff_size+1:n, :)
            q_cons_vf(i)%sf(:, n+1:n+buff_size, :) = q_cons_vf(i)%sf(:, 0:buff_size-1, :)

            q_cons_vf(i)%sf(:, :, -buff_size:-1) = q_cons_vf(i)%sf(:, :, p-buff_size+1:p)
            q_cons_vf(i)%sf(:, :, p+1:p+buff_size) = q_cons_vf(i)%sf(:, :, 0:buff_size-1)
        end do
#endif
        
        ! effective viscosity setup, return viscous stress tensor
        call s_compute_viscous_stress_tensor(visc_stress, q_prim_vf=q_prim_vf)

        call s_compute_stress_tensor(pres_visc_stress, visc_stress, q_cons_vf, q_prim_vf)

        ! interphase momentum exchange term setup
        call s_compute_divergence_stress_tensor(div_pres_visc_stress, pres_visc_stress)

    end subroutine s_setup_terms_filtering

    subroutine s_compute_pseudo_turbulent_reynolds_stress(q_cons_filtered, reynolds_stress)
        type(scalar_field), dimension(sys_size-1), intent(in) :: q_cons_filtered
        type(vector_field), dimension(1:num_dims), intent(inout) :: reynolds_stress
        integer :: i, j, k, l, q    

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p
                    !$acc loop seq
                    do l = 1, num_dims
                        !$acc loop seq
                        do q = 1, num_dims
                            reynolds_stress(l)%vf(q)%sf(i, j, k) = reynolds_stress(l)%vf(q)%sf(i, j, k) &
                                - (q_cons_filtered(momxb-1+l)%sf(i, j, k) * q_cons_filtered(momxb-1+q)%sf(i, j, k) / q_cons_filtered(1)%sf(i, j, k))
                        end do
                    end do
                end do
            end do
        end do

    end subroutine s_compute_pseudo_turbulent_reynolds_stress

    subroutine s_compute_effective_viscosity(q_cons_filtered, eff_visc, visc_stress)
        type(scalar_field), dimension(1:sys_size-1), intent(inout) :: q_cons_filtered
        type(vector_field), dimension(1:num_dims), intent(inout) :: eff_visc
        type(vector_field), dimension(1:num_dims), intent(inout) :: visc_stress

        integer :: i, j, k, l, q

        ! set buffers for filtered momentum quantities and density
#ifdef MFC_MPI
        do i = contxb, momxe 
            call s_populate_scalarfield_buffers(q_cons_filtered(i))
        end do
#else
        do i = contxb, momxe
            q_cons_filtered(i)%sf(-buff_size:-1, :, :) = q_cons_filtered(i)%sf(m-buff_size+1:m, :, :)
            q_cons_filtered(i)%sf(m+1:m+buff_size, :, :) = q_cons_filtered(i)%sf(0:buff_size-1, :, :)

            q_cons_filtered(i)%sf(:, -buff_size:-1, :) = q_cons_filtered(i)%sf(:, n-buff_size+1:n, :)
            q_cons_filtered(i)%sf(:, n+1:n+buff_size, :) = q_cons_filtered(i)%sf(:, 0:buff_size-1, :)

            q_cons_filtered(i)%sf(:, :, -buff_size:-1) = q_cons_filtered(i)%sf(:, :, p-buff_size+1:p)
            q_cons_filtered(i)%sf(:, :, p+1:p+buff_size) = q_cons_filtered(i)%sf(:, :, 0:buff_size-1)
        end do
#endif

        ! calculate stress tensor with filtered quantities 
        call s_compute_viscous_stress_tensor(visc_stress, q_cons_filtered=q_cons_filtered)

        ! calculate eff_visc
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    !$acc loop seq
                    do l = 1, num_dims
                        !$acc loop seq
                        do q = 1, num_dims
                            eff_visc(l)%vf(q)%sf(i, j, k) = eff_visc(l)%vf(q)%sf(i, j, k) - visc_stress(l)%vf(q)%sf(i, j, k)
                        end do 
                    end do
                end do
            end do
        end do

    end subroutine s_compute_effective_viscosity

    subroutine s_compute_interphase_momentum_exchange(filtered_fluid_indicator_function, grad_fluid_indicator, pres_visc_stress, int_mom_exch)
        type(scalar_field), intent(in) :: filtered_fluid_indicator_function
        type(scalar_field), dimension(1:3), intent(in) :: grad_fluid_indicator
        type(vector_field), dimension(1:3), intent(in) :: pres_visc_stress
        type(scalar_field), dimension(1:3), intent(inout) :: int_mom_exch

        integer :: i, j, k, l

        ! x-, y-, z- component loop
        do l = 1, 3

            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p
                        data_real_3D_slabz(i+1, j+1, k+1) = pres_visc_stress(1)%vf(l)%sf(i, j, k) * grad_fluid_indicator(1)%sf(i, j, k) & 
                                                          + pres_visc_stress(2)%vf(l)%sf(i, j, k) * grad_fluid_indicator(2)%sf(i, j, k) & 
                                                          + pres_visc_stress(3)%vf(l)%sf(i, j, k) * grad_fluid_indicator(3)%sf(i, j, k)
                    end do 
                end do
            end do

            call nvtxStartRange("FORWARD-3D-FFT")
            call s_mpi_FFT_fwd 
            call nvtxEndRange

            ! convolution with filtering kernel
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, NxC 
                do j = 1, Nyloc 
                    do k = 1, Nz 
                        data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) * cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                    end do 
                end do 
            end do

            call nvtxStartRange("BACKWARD-3D-FFT")
            call s_mpi_FFT_bwd
            call nvtxEndRange

            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        int_mom_exch(l)%sf(i, j, k) = data_real_3D_slabz(i+1, j+1, k+1) / (real(Nx*Ny*Nz, dp))
                    end do 
                end do 
            end do
        end do ! end component loop

    end subroutine s_compute_interphase_momentum_exchange

    ! computes x-,y-,z-direction forces on particles
    subroutine s_compute_particle_forces
        real(wp), dimension(num_ibs, 3) :: force_glb
        real(wp) :: dvol
        integer :: i, j, k, l

        ! zero particle forces
        particle_forces = 0.0_wp
        !$acc update device(particle_forces)

        !$acc parallel loop collapse(3) gang vector default(present) private(dvol)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p
                    dvol = dx(i) * dy(j) * dz(k)
                    !$acc atomic
                    particle_forces(ib_markers%sf(i, j, k), 1) = particle_forces(ib_markers%sf(i, j, k), 1) - div_pres_visc_stress(1)%sf(i, j, k) * dvol
                    !$acc atomic
                    particle_forces(ib_markers%sf(i, j, k), 2) = particle_forces(ib_markers%sf(i, j, k), 2) - div_pres_visc_stress(2)%sf(i, j, k) * dvol
                    !$acc atomic
                    particle_forces(ib_markers%sf(i, j, k), 3) = particle_forces(ib_markers%sf(i, j, k), 3) - div_pres_visc_stress(3)%sf(i, j, k) * dvol
                end do 
            end do 
        end do

        !$acc update host(particle_forces)

        ! reduce particle forces across processors
        do i = 1, num_ibs
            call s_mpi_allreduce_sum(particle_forces(i, 1), force_glb(i, 1))
            call s_mpi_allreduce_sum(particle_forces(i, 2), force_glb(i, 2))
            call s_mpi_allreduce_sum(particle_forces(i, 3), force_glb(i, 3))
        end do

        ! if (proc_rank == 0) then
        !     print *, 'force', force_glb(1, 1)
        !     print *, 'C_D', 2._wp * force_glb(1, 1) / (rho_inf_ref * u_inf_ref**2 * pi * patch_ib(1)%radius**2)
        ! end if
        
        ! write particle forces to file
        if (proc_rank == 0) then
            write(100) force_glb
        end if
            
    end subroutine s_compute_particle_forces


    !< transpose domain from z-slabs to y-slabs on each processor
    subroutine s_mpi_transpose_slabZ2Y
        integer :: dest_rank, src_rank
        integer :: i, j, k

        !$acc parallel loop collapse(4) gang vector default(present) copy(sendbuf_sf)
        do dest_rank = 0, num_procs-1
            do k = 1, Nzloc 
                do j = 1, Nyloc
                    do i = 1, NxC
                        sendbuf_sf(i + (j-1)*NxC + (k-1)*NxC*Nyloc + dest_rank*NxC*Nyloc*Nzloc) = data_cmplx_slabz(i, j+dest_rank*Nyloc, k)
                    end do 
                end do
            end do
        end do

        call MPI_Alltoall(sendbuf_sf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, & 
                          recvbuf_sf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop collapse(4) gang vector default(present) copy(recvbuf_sf)
        do src_rank = 0, num_procs-1
            do k = 1, Nzloc 
                do j = 1, Nyloc
                    do i = 1, NxC
                        data_cmplx_slaby(i, j, k+src_rank*Nzloc) = recvbuf_sf(i + (j-1)*NxC + (k-1)*NxC*Nyloc + src_rank*NxC*Nyloc*Nzloc)
                    end do 
                end do
            end do 
        end do

    end subroutine s_mpi_transpose_slabZ2Y

    !< transpose domain from y-slabs to z-slabs on each processor
    subroutine s_mpi_transpose_slabY2Z 
        integer :: dest_rank, src_rank
        integer :: i, j, k

        !$acc parallel loop collapse(4) gang vector default(present) copy(sendbuf_sf)
        do dest_rank = 0, num_procs-1
            do k = 1, Nzloc 
                do j = 1, Nyloc 
                    do i = 1, NxC 
                        sendbuf_sf(i + (j-1)*NxC + (k-1)*NxC*Nyloc + dest_rank*NxC*Nyloc*Nzloc) = data_cmplx_slaby(i, j, k+dest_rank*Nzloc)
                    end do 
                end do 
            end do 
        end do

        call MPI_Alltoall(sendbuf_sf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, & 
                          recvbuf_sf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop collapse(4) gang vector default(present) copy(recvbuf_sf) 
        do src_rank = 0, num_procs-1
            do k = 1, Nzloc
                do j = 1, Nyloc 
                    do i = 1, NxC 
                        data_cmplx_slabz(i, j+src_rank*Nyloc, k) = recvbuf_sf(i + (j-1)*NxC + (k-1)*NxC*Nyloc + src_rank*NxC*Nyloc*Nzloc)
                    end do 
                end do
            end do 
        end do
        
    end subroutine s_mpi_transpose_slabY2Z

    !< transpose domain from z-slabs to y-slabs on each processor for batched 9 element tensors
    subroutine s_mpi_transpose_slabZ2Y_tensor
        integer :: dest_rank, src_rank
        integer :: i, j, k, l

        !$acc parallel loop collapse(5) gang vector default(present) copy(sendbuf_tensor)
        do dest_rank = 0, num_procs-1
            do k = 1, Nzloc 
                do j = 1, Nyloc
                    do i = 1, NxC
                        do l = 1, 9
                            sendbuf_tensor(l + (i-1)*9 + (j-1)*9*NxC + (k-1)*9*NxC*Nyloc + dest_rank*9*NxC*Nyloc*Nzloc) = data_cmplx_slabz_tensor(l, i, j+dest_rank*Nyloc, k)
                        end do 
                    end do
                end do
            end do
        end do 

        call MPI_Alltoall(sendbuf_tensor, 9*NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, & 
                          recvbuf_tensor, 9*NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop collapse(5) gang vector default(present) copy(recvbuf_tensor)
        do src_rank = 0, num_procs-1
            do k = 1, Nzloc 
                do j = 1, Nyloc
                    do i = 1, NxC
                        do l = 1, 9
                            data_cmplx_slaby_tensor(l, i, j, k+src_rank*Nzloc) = recvbuf_tensor(l + (i-1)*9 + (j-1)*9*NxC + (k-1)*9*NxC*Nyloc + src_rank*9*NxC*Nyloc*Nzloc)
                        end do 
                    end do
                end do 
            end do
        end do

    end subroutine s_mpi_transpose_slabZ2Y_tensor

    !< transpose domain from y-slabs to z-slabs on each processor for batched 9 element tensors
    subroutine s_mpi_transpose_slabY2Z_tensor
        integer :: dest_rank, src_rank
        integer :: i, j, k, l

        !$acc parallel loop collapse(5) gang vector default(present) copy(sendbuf_tensor)
        do dest_rank = 0, num_procs-1
            do k = 1, Nzloc 
                do j = 1, Nyloc 
                    do i = 1, NxC 
                        do l = 1, 9
                            sendbuf_tensor(l + (i-1)*9 + (j-1)*9*NxC + (k-1)*9*NxC*Nyloc + dest_rank*9*NxC*Nyloc*Nzloc) = data_cmplx_slaby_tensor(l, i, j, k+dest_rank*Nzloc)
                        end do 
                    end do 
                end do 
            end do
        end do

        call MPI_Alltoall(sendbuf_tensor, 9*NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, & 
                          recvbuf_tensor, 9*NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop collapse(5) gang vector default(present) copy(recvbuf_tensor) 
        do src_rank = 0, num_procs-1
            do k = 1, Nzloc
                do j = 1, Nyloc 
                    do i = 1, NxC 
                        do l = 1, 9
                            data_cmplx_slabz_tensor(l, i, j+src_rank*Nyloc, k) = recvbuf_tensor(l + (i-1)*9 + (j-1)*9*NxC + (k-1)*9*NxC*Nyloc + src_rank*9*NxC*Nyloc*Nzloc)
                        end do 
                    end do
                end do 
            end do
        end do
        
    end subroutine s_mpi_transpose_slabY2Z_tensor



    !< compute forward FFT, input: data_real_3D_slabz, output: data_cmplx_out1d
    subroutine s_filter_tensor_field(q_tensor_in, q_tensor_out)
        type(vector_field), dimension(3), intent(inout) :: q_tensor_in
        type(vector_field), dimension(3), intent(inout), optional :: q_tensor_out
        integer :: i, j, k, l, q

        ! ===== forward FFT =====
        ! outer tensor element loop
        do l = 1, 3
            do q = 1, 3

                !$acc parallel loop collapse(3)
                do i = 0, m 
                    do j = 0, n 
                        do k = 0, p 
                            data_real_3D_slabz(i+1, j+1, k+1) = q_tensor_in(l)%vf(q)%sf(i, j, k) * fluid_indicator_function%sf(i, j, k)
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
                            data_cmplx_slabz_tensor((l-1)*3 + q, i, j, k) = data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC)
                        end do 
                    end do 
                end do 
                ! pack data_cmplx_slabz_tensor for MPI tranpose
            end do
        end do 

        ! tensor MPI data transpose
        call s_mpi_transpose_slabZ2Y_tensor

        ! outer tensor element loop
        do l = 1, 3
            do q = 1, 3
                ! 3D y-slab -> 1D z, x, y
                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, NxC 
                    do j = 1, Nyloc 
                        do k = 1, Nz
                            data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_slaby_tensor((l-1)*3 + q, i, j, k)
                        end do 
                    end do 
                end do

                ! Z FFT
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_FORWARD)
#else
                call fftw_execute_dft(plan_z_c2c_fwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
                
                ! convolution with filtering kernel in Fourier space
                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, NxC 
                    do j = 1, Nyloc 
                        do k = 1, Nz 
                            data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) * cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                        end do 
                    end do 
                end do

                ! ===== begin backward FFT =====
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
                            data_cmplx_slaby_tensor((l-1)*3 + q, i, j, k) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                        end do 
                    end do 
                end do
                ! pack data_cmplx_slaby_tensor for MPI tranpose
            end do
        end do

        call s_mpi_transpose_slabY2Z_tensor

        ! outer tensor element loop
        do l = 1, 3
            do q = 1, 3
                
                ! 3D z-slab -> 1D y, x, z
                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, NxC 
                    do j = 1, Ny 
                        do k = 1, Nzloc
                            data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC) = data_cmplx_slabz_tensor((l-1)*3 + q, i, j, k)
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

                if (present(q_tensor_out)) then 
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = 0, m
                        do j = 0, n
                            do k = 0, p
                                q_tensor_out(l)%vf(q)%sf(i, j, k) = data_real_3D_slabz(i+1, j+1, k+1) / (real(Nx*Ny*Nz, dp) * filtered_fluid_indicator_function%sf(i, j, k))
                            end do 
                        end do 
                    end do
                else
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = 0, m
                        do j = 0, n
                            do k = 0, p
                                q_tensor_in(l)%vf(q)%sf(i, j, k) = data_real_3D_slabz(i+1, j+1, k+1) / (real(Nx*Ny*Nz, dp) * filtered_fluid_indicator_function%sf(i, j, k))
                            end do 
                        end do 
                    end do
                end if

            end do
        end do

    end subroutine s_filter_tensor_field



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
        call nvtxStartRange("SLAB-MPI-TRANSPOSE-Z2Y")
        call s_mpi_transpose_slabZ2Y 
        call nvtxEndRange

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
        call nvtxStartRange("SLAB-MPI-TRANSPOSE-Y2Z")
        call s_mpi_transpose_slabY2Z
        call nvtxEndRange

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

    subroutine s_finalize_fftw_explicit_filter_module
        integer :: i, j 

        @:DEALLOCATE(fluid_indicator_function%sf)
        @:DEALLOCATE(filtered_fluid_indicator_function%sf)
        do i = 1, 3 
            @:DEALLOCATE(grad_fluid_indicator(i)%sf)
        end do
        @:DEALLOCATE(grad_fluid_indicator)

        do i = 1, sys_size-1
            @:DEALLOCATE(q_cons_filtered(i)%sf)
        end do
        @:DEALLOCATE(q_cons_filtered)

        @:DEALLOCATE(filtered_pressure%sf)

        do i = 1, num_dims
            do j = 1, num_dims
                @:DEALLOCATE(visc_stress(i)%vf(j)%sf)
            end do 
            @:DEALLOCATE(visc_stress(i)%vf)
        end do
        @:DEALLOCATE(visc_stress)

        do i = 1, num_dims
            do j = 1, num_dims
                @:DEALLOCATE(pres_visc_stress(i)%vf(j)%sf)
            end do
            @:DEALLOCATE(pres_visc_stress(i)%vf)
        end do
        @:DEALLOCATE(pres_visc_stress)

        do i = 1, num_dims
            @:DEALLOCATE(div_pres_visc_stress(i)%sf)
        end do
        @:DEALLOCATE(div_pres_visc_stress)

        do i = 1, num_dims
            do j = 1, num_dims
                @:DEALLOCATE(reynolds_stress(i)%vf(j)%sf)
            end do
            @:DEALLOCATE(reynolds_stress(i)%vf)
        end do
        @:DEALLOCATE(reynolds_stress)

        do i = 1, num_dims
            do j = 1, num_dims
                @:DEALLOCATE(eff_visc(i)%vf(j)%sf)
            end do
            @:DEALLOCATE(eff_visc(i)%vf)
        end do
        @:DEALLOCATE(eff_visc)

        do i = 1, num_dims
            @:DEALLOCATE(int_mom_exch(i)%sf)
        end do
        @:DEALLOCATE(int_mom_exch)

        @:DEALLOCATE(Res)
        @:DEALLOCATE(particle_forces)

        @:DEALLOCATE(data_real_in1d, data_cmplx_out1d, data_cmplx_out1dy)
        @:DEALLOCATE(cmplx_kernelG1d, real_kernelG_in)
        @:DEALLOCATE(data_real_3D_slabz, data_cmplx_slabz, data_cmplx_slaby)
        @:DEALLOCATE(data_cmplx_slabz_tensor, data_cmplx_slaby_tensor)
        
        deallocate(sendbuf_sf, recvbuf_sf)
        deallocate(sendbuf_tensor, recvbuf_tensor)

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

        if (proc_rank == 0) then
            close(100)
        end if

    end subroutine s_finalize_fftw_explicit_filter_module

end module m_volume_filtering