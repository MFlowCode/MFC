#:include 'macros.fpp'

module m_volume_filtering

    use, intrinsic :: iso_c_binding

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_ibm

    use m_boundary_common

    use m_nvtx

    use m_variables_conversion

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
 s_finalize_fftw_explicit_filter_module, s_volume_filter_momentum_eqn, s_compute_particle_forces

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
    type(scalar_field), allocatable, dimension(:) :: q_prim_filtered
    type(scalar_field), public :: filtered_pressure

    ! viscous and pressure+viscous stress tensors
    type(scalar_field), allocatable, dimension(:, :) :: visc_stress
    type(scalar_field), allocatable, dimension(:, :) :: pres_visc_stress

    ! divergence of stress tensor
    type(scalar_field), allocatable, dimension(:) :: div_pres_visc_stress

    ! unclosed terms in volume filtered momentum equation
    type(scalar_field), allocatable, dimension(:, :), public :: reynolds_stress
    type(scalar_field), allocatable, dimension(:, :), public :: eff_visc
    type(scalar_field), allocatable, dimension(:), public :: int_mom_exch

    ! x-,y-,z-direction forces on particles
    real(wp), allocatable, dimension(:, :) :: particle_forces

    $:GPU_DECLARE(create='[fluid_indicator_function, filtered_fluid_indicator_function, grad_fluid_indicator]')
    $:GPU_DECLARE(create='[q_cons_filtered, filtered_pressure, visc_stress, pres_visc_stress, div_pres_visc_stress]')
    $:GPU_DECLARE(create='[reynolds_stress, eff_visc, int_mom_exch, particle_forces]')

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
    real(wp) :: fft_norm

    ! 1D real and complex vectors for FFT routines
    real(c_double), allocatable :: data_real_in1d(:)
    complex(c_double_complex), allocatable :: data_cmplx_out1d(:)
    complex(c_double_complex), allocatable :: data_cmplx_out1dy(:)

    ! 3D arrays for slab transposes
    complex(c_float_complex), allocatable :: data_cmplx_slabz(:, :, :), data_cmplx_slaby(:, :, :)
    ! 3D arrays for slab transposes of tensor quantities
    complex(c_float_complex), allocatable :: data_cmplx_slabz_batch(:, :, :, :), data_cmplx_slaby_batch(:, :, :, :)

    ! input/output array for FFT routine
    real(c_double), allocatable :: data_real_3D_slabz(:, :, :)

    ! filtering kernel in physical space
    real(c_double), allocatable :: real_kernelG_in(:, :, :)

    ! FFT of filtering kernel
    complex(c_double_complex), allocatable :: cmplx_kernelG1d(:)

    $:GPU_DECLARE(create='[Nx, Ny, Nz, NxC, Nyloc, Nzloc, fft_norm]')
    $:GPU_DECLARE(create='[data_real_in1d, data_cmplx_out1d, data_cmplx_out1dy, data_cmplx_slabz, data_cmplx_slaby, data_cmplx_slabz_batch, data_cmplx_slaby_batch]')
    $:GPU_DECLARE(create='[data_real_3D_slabz, real_kernelG_in, cmplx_kernelG1d]')

    ! buffers for data transpose
    complex(c_float_complex), allocatable :: sendbuf_sf(:), recvbuf_sf(:)
    complex(c_float_complex), allocatable :: sendbuf_batch(:), recvbuf_batch(:)

    $:GPU_DECLARE(create='[sendbuf_sf, recvbuf_sf, sendbuf_batch, recvbuf_batch]')

    ! total batch size for MPI_Alltoall used in 3D FFT
    integer :: fft_batch_size
    ! batch index locations
    integer :: reynolds_stress_idx, eff_visc_idx, int_mom_exch_idx

    $:GPU_DECLARE(create='[fft_batch_size, reynolds_stress_idx, eff_visc_idx, int_mom_exch_idx]')

contains

    !< create fft plans to be used for explicit filtering of data
    subroutine s_initialize_fftw_explicit_filter_module
        integer :: i, j, k
        integer :: size_n(1), inembed(1), onembed(1)

        @:ALLOCATE(q_cons_filtered(1:sys_size))
        do i = 1, sys_size
            @:ALLOCATE(q_cons_filtered(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_cons_filtered(i))
        end do

        @:ALLOCATE(q_prim_filtered(1:sys_size))
        do i = 1, sys_size
            @:ALLOCATE(q_prim_filtered(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_prim_filtered(i))
        end do

        @:ALLOCATE(filtered_pressure%sf(idwbuff(1)%beg:idwbuff(1)%end, &
            idwbuff(2)%beg:idwbuff(2)%end, &
            idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(filtered_pressure)

        @:ALLOCATE(visc_stress(1:num_dims, 1:num_dims))
        do i = 1, num_dims
            do j = 1, num_dims
                @:ALLOCATE(visc_stress(i, j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(visc_stress(i, j))
            end do
        end do

        @:ALLOCATE(pres_visc_stress(1:num_dims, 1:num_dims))
        do i = 1, num_dims
            do j = 1, num_dims
                @:ALLOCATE(pres_visc_stress(i, j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(pres_visc_stress(i, j))
            end do
        end do

        @:ALLOCATE(div_pres_visc_stress(1:num_dims))
        do i = 1, num_dims
            @:ALLOCATE(div_pres_visc_stress(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(div_pres_visc_stress(i))
        end do

        @:ALLOCATE(reynolds_stress(1:num_dims, 1:num_dims))
        do i = 1, num_dims
            do j = 1, num_dims
                @:ALLOCATE(reynolds_stress(i, j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(reynolds_stress(i, j))
            end do
        end do

        @:ALLOCATE(eff_visc(1:num_dims, 1:num_dims))
        do i = 1, num_dims
            do j = 1, num_dims
                @:ALLOCATE(eff_visc(i, j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(eff_visc(i, j))
            end do
        end do

        @:ALLOCATE(int_mom_exch(1:num_dims))
        do i = 1, num_dims
            @:ALLOCATE(int_mom_exch(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(int_mom_exch(i))
        end do

        @:ALLOCATE(particle_forces(0:num_ibs, 3))

        !< global sizes
        Nx = m_glb + 1
        Ny = n_glb + 1
        Nz = p_glb + 1
        fft_norm = 1._wp/real(Nx*Ny*Nz, wp)

        if (mod(Ny, num_procs) /= 0) then
            call s_mpi_abort('Volume filtering requires p to be divisible by the number of ranks')
        end if

        !< complex size
        NxC = Nx/2 + 1

        !< local sizes on each processor
        Nyloc = Ny/num_procs
        Nzloc = p + 1

        !< batch size used in MPI_Alltoall
        fft_batch_size = sys_size + 1 + 9 + 9 + 3 ! conservative vars, pressure, reynolds stress, viscous stress, interphase momentum exchange
        reynolds_stress_idx = sys_size + 1
        eff_visc_idx = reynolds_stress_idx + 9
        int_mom_exch_idx = eff_visc_idx + 9

        $:GPU_UPDATE(device='[Nx, Ny, Nz, fft_norm, NxC, Nyloc, Nzloc, fft_batch_size, reynolds_stress_idx, eff_visc_idx, int_mom_exch_idx]')

        @:ALLOCATE(data_real_in1d(Nx*Ny*Nzloc))
        @:ALLOCATE(data_cmplx_out1d(NxC*Ny*Nz/num_procs))
        @:ALLOCATE(data_cmplx_out1dy(NxC*Ny*Nz/num_procs))
        @:ALLOCATE(cmplx_kernelG1d(NxC*Nyloc*Nz))
        @:ALLOCATE(real_kernelG_in(Nx, Ny, Nzloc))
        @:ALLOCATE(data_real_3D_slabz(Nx, Ny, Nzloc))
        @:ALLOCATE(data_cmplx_slabz(NxC, Ny, Nzloc))
        @:ALLOCATE(data_cmplx_slaby(NxC, Nyloc, Nz))
        @:ALLOCATE(data_cmplx_slabz_batch(fft_batch_size, NxC, Ny, Nzloc))
        @:ALLOCATE(data_cmplx_slaby_batch(fft_batch_size, NxC, Nyloc, Nz))

        @:ALLOCATE(sendbuf_sf(NxC*Nyloc*Nzloc*num_procs))
        @:ALLOCATE(recvbuf_sf(NxC*Nyloc*Nzloc*num_procs))
        @:ALLOCATE(sendbuf_batch(fft_batch_size*NxC*Nyloc*Nzloc*num_procs))
        @:ALLOCATE(recvbuf_batch(fft_batch_size*NxC*Nyloc*Nzloc*num_procs))

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
        if (compute_particle_drag) then
            if (proc_rank == 0) then
                open (unit=100, file='particle_force.bin', status='replace', form='unformatted', access='stream', action='write')
            end if
        end if

    end subroutine s_initialize_fftw_explicit_filter_module

    !< initialize the gaussian filtering kernel in real space and then compute its DFT
    subroutine s_initialize_filtering_kernel
        real(wp) :: sigma_stddev
        real(wp) :: Lx, Ly, Lz
        real(wp) :: x_r, y_r, z_r
        real(wp) :: r2
        real(wp) :: G_norm_int, G_norm_int_glb
        integer :: i, j, k

        ! gaussian filter
        sigma_stddev = filter_width

        Lx = domain_glb(1, 2) - domain_glb(1, 1)
        Ly = domain_glb(2, 2) - domain_glb(2, 1)
        Lz = domain_glb(3, 2) - domain_glb(3, 1)

        G_norm_int = 0.0_wp

        $:GPU_PARALLEL_LOOP(collapse=3, reduction='[[G_norm_int]]', reductionOp='[+]', copyin='[Lx, Ly, Lz, sigma_stddev]', private='[x_r, y_r, z_r, r2]')
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    x_r = min(abs(x_cc(i) - domain_glb(1, 1)), Lx - abs(x_cc(i) - domain_glb(1, 1)))
                    y_r = min(abs(y_cc(j) - domain_glb(2, 1)), Ly - abs(y_cc(j) - domain_glb(2, 1)))
                    z_r = min(abs(z_cc(k) - domain_glb(3, 1)), Lz - abs(z_cc(k) - domain_glb(3, 1)))

                    r2 = x_r**2 + y_r**2 + z_r**2

                    real_kernelG_in(i + 1, j + 1, k + 1) = exp(-r2/(2.0_wp*sigma_stddev**2))

                    G_norm_int = G_norm_int + real_kernelG_in(i + 1, j + 1, k + 1)*dx(i)*dy(j)*dz(k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_mpi_allreduce_sum(G_norm_int, G_norm_int_glb)

        ! FFT of kernel
        ! normalize kernel
        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[G_norm_int_glb]')
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nzloc
                    data_real_3D_slabz(i, j, k) = real_kernelG_in(i, j, k)/G_norm_int_glb
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! 3D z-slab -> 1D x, y, z
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nzloc
                    data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! X FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, cmplx_kernelG1d)
#else
        call fftw_execute_dft_r2c(plan_x_r2c_kernelG, data_real_in1d, cmplx_kernelG1d)
#endif

        ! 1D x, y, z -> 1D y, x, z (CMPLX)
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = cmplx_kernelG1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Y FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_y_c2c_kernelG, data_cmplx_out1dy, data_cmplx_out1dy)
#endif

        ! 1D y, x, z -> 3D z-slab
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_slabz(i, j, k) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! transpose z-slab to y-slab
        call s_mpi_transpose_slabZ2Y

        ! 3D y-slab -> 1D z, x, y
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Nyloc
                do k = 1, Nz
                    cmplx_kernelG1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_slaby(i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Z FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_z_gpu, cmplx_kernelG1d, cmplx_kernelG1d, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_z_c2c_kernelG, cmplx_kernelG1d, cmplx_kernelG1d)
#endif

        ! normalize FFT
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Nyloc
                do k = 1, Nz
                    cmplx_kernelG1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = cmplx_kernelG1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)*fft_norm
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! return cmplx_kernelG1d: 1D z, x, y
    end subroutine s_initialize_filtering_kernel

    !< initialize fluid indicator function and filtered fluid indicator function
    subroutine s_initialize_fluid_indicator_function(bc_type)
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type
        integer :: i, j, k

        @:ALLOCATE(fluid_indicator_function%sf(idwbuff(1)%beg:idwbuff(1)%end, &
            idwbuff(2)%beg:idwbuff(2)%end, &
            idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(fluid_indicator_function)

        ! define fluid indicator function
        if (ib) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) == 0) then
                            fluid_indicator_function%sf(i, j, k) = 1.0_wp
                        else
                            fluid_indicator_function%sf(i, j, k) = 0.0_wp
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        fluid_indicator_function%sf(i, j, k) = 1.0_wp
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        call s_populate_scalarfield_buffers(bc_type, fluid_indicator_function)

    end subroutine s_initialize_fluid_indicator_function

    subroutine s_initialize_filtered_fluid_indicator_function
        integer :: i, j, k

        @:ALLOCATE(filtered_fluid_indicator_function%sf(0:m, 0:n, 0:p))
        @:ACC_SETUP_SFs(filtered_fluid_indicator_function)

        ! filter fluid indicator function
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nzloc
                    data_real_3D_slabz(i, j, k) = fluid_indicator_function%sf(i - 1, j - 1, k - 1)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_mpi_FFT_fwd

        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Nyloc
                do k = 1, Nz
                    data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)*cmplx_kernelG1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_mpi_FFT_bwd

        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nzloc
                    filtered_fluid_indicator_function%sf(i - 1, j - 1, k - 1) = data_real_3D_slabz(i, j, k)*fft_norm
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_initialize_filtered_fluid_indicator_function

    subroutine s_initialize_fluid_indicator_gradient
        real(wp), dimension(-4:4) :: fd_coeffs
        integer :: i, j, k, l

        fd_coeffs(-4) = 1._wp/280._wp
        fd_coeffs(-3) = -4._wp/105._wp
        fd_coeffs(-2) = 1._wp/5._wp
        fd_coeffs(-1) = -4._wp/5._wp
        fd_coeffs(0) = 0._wp
        fd_coeffs(1) = 4._wp/5._wp
        fd_coeffs(2) = -1._wp/5._wp
        fd_coeffs(3) = 4._wp/105._wp
        fd_coeffs(4) = -1._wp/280._wp

        @:ALLOCATE(grad_fluid_indicator(1:3))
        do i = 1, num_dims
            @:ALLOCATE(grad_fluid_indicator(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(grad_fluid_indicator(i))
            grad_fluid_indicator(i)%sf = 0._wp
            $:GPU_UPDATE(device='[grad_fluid_indicator(i)%sf]')
        end do

        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[fd_coeffs]')
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    $:GPU_LOOP(parallelism='[seq]')
                    do l = -4, 4
                        grad_fluid_indicator(1)%sf(i, j, k) = grad_fluid_indicator(1)%sf(i, j, k) + fd_coeffs(l)*fluid_indicator_function%sf(i + l, j, k)/dx(i)
                        grad_fluid_indicator(2)%sf(i, j, k) = grad_fluid_indicator(2)%sf(i, j, k) + fd_coeffs(l)*fluid_indicator_function%sf(i, j + l, k)/dy(j)
                        grad_fluid_indicator(3)%sf(i, j, k) = grad_fluid_indicator(3)%sf(i, j, k) + fd_coeffs(l)*fluid_indicator_function%sf(i, j, k + l)/dz(k)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_initialize_fluid_indicator_gradient

    !< calculate the unclosed terms present in the volume filtered momentum equation
    subroutine s_volume_filter_momentum_eqn(q_cons_vf, q_prim_vf, q_T_sf, dyn_visc, bc_type)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), intent(inout) :: q_T_sf
        real(wp), intent(in) :: dyn_visc
        type(integer_field), dimension(num_dims, -1:1), intent(in) :: bc_type
        integer :: i, j, k

        call s_setup_terms_filtering(q_cons_vf, q_prim_vf, dyn_visc, reynolds_stress, visc_stress, pres_visc_stress, div_pres_visc_stress, bc_type)

        call s_filter_batch(q_cons_vf, q_cons_filtered, q_prim_vf(E_idx), filtered_pressure, reynolds_stress, visc_stress, eff_visc, int_mom_exch)

        ! generate Favre filtered primitives
        call s_convert_conservative_to_primitive_variables(q_cons_filtered, q_T_sf, q_prim_filtered, idwint)

        call s_compute_pseudo_turbulent_reynolds_stress(q_cons_filtered, reynolds_stress)

        call s_compute_effective_viscosity(q_prim_filtered, eff_visc, visc_stress, dyn_visc, bc_type)

    end subroutine s_volume_filter_momentum_eqn

    ! compute viscous stress tensor
    subroutine s_generate_viscous_stress_tensor(visc_stress, q_prim_vf, dyn_visc)
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: visc_stress
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), intent(in) :: dyn_visc
        real(wp) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz ! spatial velocity derivatives
        integer :: i, j, k

        $:GPU_PARALLEL_LOOP(collapse=3, private='[dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz]')
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    ! velocity gradients, local to each process
                    dudx = (q_prim_vf(2)%sf(i + 1, j, k) - q_prim_vf(2)%sf(i - 1, j, k))/(2._wp*dx(i))
                    dudy = (q_prim_vf(2)%sf(i, j + 1, k) - q_prim_vf(2)%sf(i, j - 1, k))/(2._wp*dy(j))
                    dudz = (q_prim_vf(2)%sf(i, j, k + 1) - q_prim_vf(2)%sf(i, j, k - 1))/(2._wp*dz(k))

                    dvdx = (q_prim_vf(3)%sf(i + 1, j, k) - q_prim_vf(3)%sf(i - 1, j, k))/(2._wp*dx(i))
                    dvdy = (q_prim_vf(3)%sf(i, j + 1, k) - q_prim_vf(3)%sf(i, j - 1, k))/(2._wp*dy(j))
                    dvdz = (q_prim_vf(3)%sf(i, j, k + 1) - q_prim_vf(3)%sf(i, j, k - 1))/(2._wp*dz(k))

                    dwdx = (q_prim_vf(4)%sf(i + 1, j, k) - q_prim_vf(4)%sf(i - 1, j, k))/(2._wp*dx(i))
                    dwdy = (q_prim_vf(4)%sf(i, j + 1, k) - q_prim_vf(4)%sf(i, j - 1, k))/(2._wp*dy(j))
                    dwdz = (q_prim_vf(4)%sf(i, j, k + 1) - q_prim_vf(4)%sf(i, j, k - 1))/(2._wp*dz(k))

                    ! viscous stress tensor, visc_stress(row, column)
                    visc_stress(1, 1)%sf(i, j, k) = dyn_visc*(4._wp/3._wp*dudx - 2._wp/3._wp*(dvdy + dwdz))
                    visc_stress(1, 2)%sf(i, j, k) = dyn_visc*(dudy + dvdx)
                    visc_stress(1, 3)%sf(i, j, k) = dyn_visc*(dudz + dwdx)
                    visc_stress(2, 1)%sf(i, j, k) = dyn_visc*(dvdx + dudy)
                    visc_stress(2, 2)%sf(i, j, k) = dyn_visc*(4._wp/3._wp*dvdy - 2._wp/3._wp*(dudx + dwdz))
                    visc_stress(2, 3)%sf(i, j, k) = dyn_visc*(dvdz + dwdy)
                    visc_stress(3, 1)%sf(i, j, k) = dyn_visc*(dwdx + dudz)
                    visc_stress(3, 2)%sf(i, j, k) = dyn_visc*(dwdy + dvdz)
                    visc_stress(3, 3)%sf(i, j, k) = dyn_visc*(4._wp/3._wp*dwdz - 2._wp/3._wp*(dudx + dvdy))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_generate_viscous_stress_tensor

    subroutine s_compute_stress_tensor(pres_visc_stress, visc_stress, q_prim_vf)
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: pres_visc_stress
        type(scalar_field), dimension(num_dims, num_dims), intent(in) :: visc_stress
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer :: i, j, k

        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    pres_visc_stress(1, 1)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k) - visc_stress(1, 1)%sf(i, j, k)
                    pres_visc_stress(1, 2)%sf(i, j, k) = -visc_stress(1, 2)%sf(i, j, k)
                    pres_visc_stress(1, 3)%sf(i, j, k) = -visc_stress(1, 3)%sf(i, j, k)
                    pres_visc_stress(2, 1)%sf(i, j, k) = -visc_stress(2, 1)%sf(i, j, k)
                    pres_visc_stress(2, 2)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k) - visc_stress(2, 2)%sf(i, j, k)
                    pres_visc_stress(2, 3)%sf(i, j, k) = -visc_stress(2, 3)%sf(i, j, k)
                    pres_visc_stress(3, 1)%sf(i, j, k) = -visc_stress(3, 1)%sf(i, j, k)
                    pres_visc_stress(3, 2)%sf(i, j, k) = -visc_stress(3, 2)%sf(i, j, k)
                    pres_visc_stress(3, 3)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k) - visc_stress(3, 3)%sf(i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_stress_tensor

    !< compute the divergence of the pressure-viscous stress tensor
    subroutine s_compute_divergence_stress_tensor(div_stress_tensor, stress_tensor)
        type(scalar_field), dimension(num_dims), intent(inout) :: div_stress_tensor
        type(scalar_field), dimension(num_dims, num_dims), intent(in) :: stress_tensor
        integer :: i, j, k

        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    div_stress_tensor(1)%sf(i, j, k) = (stress_tensor(1, 1)%sf(i + 1, j, k) - stress_tensor(1, 1)%sf(i - 1, j, k))/(2._wp*dx(i)) &
                                                       + (stress_tensor(1, 2)%sf(i, j + 1, k) - stress_tensor(1, 2)%sf(i, j - 1, k))/(2._wp*dy(j)) &
                                                       + (stress_tensor(1, 3)%sf(i, j, k + 1) - stress_tensor(1, 3)%sf(i, j, k - 1))/(2._wp*dz(k))

                    div_stress_tensor(2)%sf(i, j, k) = (stress_tensor(2, 1)%sf(i + 1, j, k) - stress_tensor(2, 1)%sf(i - 1, j, k))/(2._wp*dx(i)) &
                                                       + (stress_tensor(2, 2)%sf(i, j + 1, k) - stress_tensor(2, 2)%sf(i, j - 1, k))/(2._wp*dy(j)) &
                                                       + (stress_tensor(2, 3)%sf(i, j, k + 1) - stress_tensor(2, 3)%sf(i, j, k - 1))/(2._wp*dz(k))

                    div_stress_tensor(3)%sf(i, j, k) = (stress_tensor(3, 1)%sf(i + 1, j, k) - stress_tensor(3, 1)%sf(i - 1, j, k))/(2._wp*dx(i)) &
                                                       + (stress_tensor(3, 2)%sf(i, j + 1, k) - stress_tensor(3, 2)%sf(i, j - 1, k))/(2._wp*dy(j)) &
                                                       + (stress_tensor(3, 3)%sf(i, j, k + 1) - stress_tensor(3, 3)%sf(i, j, k - 1))/(2._wp*dz(k))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_divergence_stress_tensor

    !< setup for calculation of unclosed terms in volume filtered momentum eqn
    subroutine s_setup_terms_filtering(q_cons_vf, q_prim_vf, dyn_visc, reynolds_stress, visc_stress, pres_visc_stress, div_pres_visc_stress, bc_type)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), intent(in) :: dyn_visc
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: reynolds_stress
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: visc_stress
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: pres_visc_stress
        type(scalar_field), dimension(num_dims), intent(inout) :: div_pres_visc_stress
        type(integer_field), dimension(num_dims, -1:1), intent(in) :: bc_type

        integer :: i, j, k, l, q

        ! pseudo turbulent reynolds stress setup
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    $:GPU_LOOP(parallelism='[seq]')
                    do l = 1, num_dims
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_dims
                            reynolds_stress(l, q)%sf(i, j, k) = q_cons_vf(1)%sf(i, j, k)*(q_prim_vf(momxb - 1 + l)%sf(i, j, k)*q_prim_vf(momxb - 1 + q)%sf(i, j, k)) ! rho*(u x u)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! set density and momentum buffers
        do i = contxb, momxe
            call s_populate_scalarfield_buffers(bc_type, q_prim_vf(i))
        end do

        ! effective viscosity setup, return viscous stress tensor
        call s_generate_viscous_stress_tensor(visc_stress, q_prim_vf, dyn_visc)

        call s_compute_stress_tensor(pres_visc_stress, visc_stress, q_prim_vf)

        ! set stress tensor buffers for taking divergence
        do i = 1, num_dims
            do j = 1, num_dims
                call s_populate_scalarfield_buffers(bc_type, pres_visc_stress(i, j))
            end do
        end do

        ! interphase momentum exchange term setup
        call s_compute_divergence_stress_tensor(div_pres_visc_stress, pres_visc_stress)

    end subroutine s_setup_terms_filtering

    subroutine s_compute_pseudo_turbulent_reynolds_stress(q_cons_filtered, reynolds_stress)
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_filtered
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: reynolds_stress
        integer :: i, j, k, l, q

        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    $:GPU_LOOP(parallelism='[seq]')
                    do l = 1, num_dims
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_dims
                            reynolds_stress(l, q)%sf(i, j, k) = reynolds_stress(l, q)%sf(i, j, k) &
                                                                - (q_cons_filtered(momxb - 1 + l)%sf(i, j, k)*q_cons_filtered(momxb - 1 + q)%sf(i, j, k)/q_cons_filtered(1)%sf(i, j, k))
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_pseudo_turbulent_reynolds_stress

    subroutine s_compute_effective_viscosity(q_prim_filtered, eff_visc, visc_stress, dyn_visc, bc_type)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_filtered
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: eff_visc
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: visc_stress
        real(wp), intent(in) :: dyn_visc
        type(integer_field), dimension(num_dims, -1:1), intent(in) :: bc_type

        integer :: i, j, k, l, q

        ! set buffers for filtered momentum quantities and density
        do i = momxb, momxe
            call s_populate_scalarfield_buffers(bc_type, q_prim_filtered(i))
        end do

        ! calculate stress tensor with filtered quantities
        call s_generate_viscous_stress_tensor(visc_stress, q_prim_filtered, dyn_visc)

        ! calculate eff_visc
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    $:GPU_LOOP(parallelism='[seq]')
                    do l = 1, num_dims
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_dims
                            eff_visc(l, q)%sf(i, j, k) = eff_visc(l, q)%sf(i, j, k) - visc_stress(l, q)%sf(i, j, k)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_effective_viscosity

    ! computes x-,y-,z-direction forces on particles
    subroutine s_compute_particle_forces
        real(wp), dimension(num_ibs, 3) :: force_glb
        real(wp) :: dvol
        integer :: i, j, k, l

        ! zero particle forces
        particle_forces = 0.0_wp
        $:GPU_UPDATE(device='[particle_forces]')

        $:GPU_PARALLEL_LOOP(collapse=3, private='[dvol]')
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    dvol = dx(i)*dy(j)*dz(k)
                    $:GPU_ATOMIC(atomic='update')
                    particle_forces(ib_markers%sf(i, j, k), 1) = particle_forces(ib_markers%sf(i, j, k), 1) - (div_pres_visc_stress(1)%sf(i, j, k)*dvol)
                    $:GPU_ATOMIC(atomic='update')
                    particle_forces(ib_markers%sf(i, j, k), 2) = particle_forces(ib_markers%sf(i, j, k), 2) - (div_pres_visc_stress(2)%sf(i, j, k)*dvol)
                    $:GPU_ATOMIC(atomic='update')
                    particle_forces(ib_markers%sf(i, j, k), 3) = particle_forces(ib_markers%sf(i, j, k), 3) - (div_pres_visc_stress(3)%sf(i, j, k)*dvol)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_UPDATE(host='[particle_forces]')

        ! reduce particle forces across processors
        do i = 1, num_ibs
            call s_mpi_allreduce_sum(particle_forces(i, 1), force_glb(i, 1))
            call s_mpi_allreduce_sum(particle_forces(i, 2), force_glb(i, 2))
            call s_mpi_allreduce_sum(particle_forces(i, 3), force_glb(i, 3))
        end do

        ! write particle forces to file
        if (proc_rank == 0) then
            write (100) force_glb
            flush (100)
            !print *, force_glb(1, 1) / (0.5_wp * rho_inf_ref * u_inf_ref**2 * pi * patch_ib(1)%radius**2)
        end if

    end subroutine s_compute_particle_forces

    !< transpose domain from z-slabs to y-slabs on each processor
    subroutine s_mpi_transpose_slabZ2Y
        integer :: dest_rank, src_rank
        integer :: i, j, k

        $:GPU_PARALLEL_LOOP(collapse=4)
        do dest_rank = 0, num_procs - 1
            do k = 1, Nzloc
                do j = 1, Nyloc
                    do i = 1, NxC
                        sendbuf_sf(i + (j - 1)*NxC + (k - 1)*NxC*Nyloc + dest_rank*NxC*Nyloc*Nzloc) = data_cmplx_slabz(i, j + dest_rank*Nyloc, k)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_UPDATE(host='[sendbuf_sf]')
#ifdef MFC_MPI
        call MPI_Alltoall(sendbuf_sf, NxC*Nyloc*Nzloc, MPI_COMPLEX, &
                          recvbuf_sf, NxC*Nyloc*Nzloc, MPI_COMPLEX, MPI_COMM_WORLD, ierr)
#endif
        $:GPU_UPDATE(device='[recvbuf_sf]')

        $:GPU_PARALLEL_LOOP(collapse=4)
        do src_rank = 0, num_procs - 1
            do k = 1, Nzloc
                do j = 1, Nyloc
                    do i = 1, NxC
                        data_cmplx_slaby(i, j, k + src_rank*Nzloc) = recvbuf_sf(i + (j - 1)*NxC + (k - 1)*NxC*Nyloc + src_rank*NxC*Nyloc*Nzloc)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mpi_transpose_slabZ2Y

    !< transpose domain from y-slabs to z-slabs on each processor
    subroutine s_mpi_transpose_slabY2Z
        integer :: dest_rank, src_rank
        integer :: i, j, k

        $:GPU_PARALLEL_LOOP(collapse=4)
        do dest_rank = 0, num_procs - 1
            do k = 1, Nzloc
                do j = 1, Nyloc
                    do i = 1, NxC
                        sendbuf_sf(i + (j - 1)*NxC + (k - 1)*NxC*Nyloc + dest_rank*NxC*Nyloc*Nzloc) = data_cmplx_slaby(i, j, k + dest_rank*Nzloc)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_UPDATE(host='[sendbuf_sf]')
#ifdef MFC_MPI
        call MPI_Alltoall(sendbuf_sf, NxC*Nyloc*Nzloc, MPI_COMPLEX, &
                          recvbuf_sf, NxC*Nyloc*Nzloc, MPI_COMPLEX, MPI_COMM_WORLD, ierr)
#endif
        $:GPU_UPDATE(device='[recvbuf_sf]')

        $:GPU_PARALLEL_LOOP(collapse=4)
        do src_rank = 0, num_procs - 1
            do k = 1, Nzloc
                do j = 1, Nyloc
                    do i = 1, NxC
                        data_cmplx_slabz(i, j + src_rank*Nyloc, k) = recvbuf_sf(i + (j - 1)*NxC + (k - 1)*NxC*Nyloc + src_rank*NxC*Nyloc*Nzloc)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mpi_transpose_slabY2Z

    !< transpose domain from z-slabs to y-slabs on each processor for batched fft_batch_size element tensors
    subroutine s_mpi_transpose_slabZ2Y_batch
        integer :: dest_rank, src_rank
        integer :: i, j, k, l

        $:GPU_PARALLEL_LOOP(collapse=5)
        do dest_rank = 0, num_procs - 1
            do k = 1, Nzloc
                do j = 1, Nyloc
                    do i = 1, NxC
                        do l = 1, fft_batch_size
                            sendbuf_batch(l + (i - 1)*fft_batch_size + (j - 1)*fft_batch_size*NxC + (k - 1)*fft_batch_size*NxC*Nyloc + dest_rank*fft_batch_size*NxC*Nyloc*Nzloc) = data_cmplx_slabz_batch(l, i, j + dest_rank*Nyloc, k)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_UPDATE(host='[sendbuf_batch]')
#ifdef MFC_MPI
        call MPI_Alltoall(sendbuf_batch, fft_batch_size*NxC*Nyloc*Nzloc, MPI_COMPLEX, &
                          recvbuf_batch, fft_batch_size*NxC*Nyloc*Nzloc, MPI_COMPLEX, MPI_COMM_WORLD, ierr)
#endif
        $:GPU_UPDATE(device='[recvbuf_batch]')

        $:GPU_PARALLEL_LOOP(collapse=5)
        do src_rank = 0, num_procs - 1
            do k = 1, Nzloc
                do j = 1, Nyloc
                    do i = 1, NxC
                        do l = 1, fft_batch_size
                            data_cmplx_slaby_batch(l, i, j, k + src_rank*Nzloc) = recvbuf_batch(l + (i - 1)*fft_batch_size + (j - 1)*fft_batch_size*NxC + (k - 1)*fft_batch_size*NxC*Nyloc + src_rank*fft_batch_size*NxC*Nyloc*Nzloc)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mpi_transpose_slabZ2Y_batch

    !< transpose domain from y-slabs to z-slabs on each processor for batched fft_batch_size element tensors
    subroutine s_mpi_transpose_slabY2Z_batch
        integer :: dest_rank, src_rank
        integer :: i, j, k, l

        $:GPU_PARALLEL_LOOP(collapse=5)
        do dest_rank = 0, num_procs - 1
            do k = 1, Nzloc
                do j = 1, Nyloc
                    do i = 1, NxC
                        do l = 1, fft_batch_size
                            sendbuf_batch(l + (i - 1)*fft_batch_size + (j - 1)*fft_batch_size*NxC + (k - 1)*fft_batch_size*NxC*Nyloc + dest_rank*fft_batch_size*NxC*Nyloc*Nzloc) = data_cmplx_slaby_batch(l, i, j, k + dest_rank*Nzloc)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_UPDATE(host='[sendbuf_batch]')
#ifdef MFC_MPI
        call MPI_Alltoall(sendbuf_batch, fft_batch_size*NxC*Nyloc*Nzloc, MPI_COMPLEX, &
                          recvbuf_batch, fft_batch_size*NxC*Nyloc*Nzloc, MPI_COMPLEX, MPI_COMM_WORLD, ierr)
#endif
        $:GPU_UPDATE(device='[recvbuf_batch]')

        $:GPU_PARALLEL_LOOP(collapse=5)
        do src_rank = 0, num_procs - 1
            do k = 1, Nzloc
                do j = 1, Nyloc
                    do i = 1, NxC
                        do l = 1, fft_batch_size
                            data_cmplx_slabz_batch(l, i, j + src_rank*Nyloc, k) = recvbuf_batch(l + (i - 1)*fft_batch_size + (j - 1)*fft_batch_size*NxC + (k - 1)*fft_batch_size*NxC*Nyloc + src_rank*fft_batch_size*NxC*Nyloc*Nzloc)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mpi_transpose_slabY2Z_batch

    !< compute forward FFT, input: data_real_3D_slabz, output: data_cmplx_out1d
    subroutine s_filter_batch(q_cons_vf, q_cons_filtered, pressure, filtered_pressure, reynolds_stress, visc_stress, eff_visc, int_mom_exch)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_filtered
        type(scalar_field), intent(inout) :: pressure
        type(scalar_field), intent(inout) :: filtered_pressure
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: reynolds_stress
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: visc_stress
        type(scalar_field), dimension(num_dims, num_dims), intent(inout) :: eff_visc
        type(scalar_field), dimension(num_dims), intent(inout) :: int_mom_exch
        integer :: i, j, k, l, q

        ! cons vars: X fwd FFT, Y fwd FFT
        do l = 1, sys_size
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        data_real_3D_slabz(i + 1, j + 1, k + 1) = q_cons_vf(l)%sf(i, j, k)*fluid_indicator_function%sf(i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, Nx
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, data_cmplx_out1d)
#else
            call fftw_execute_dft_r2c(plan_x_r2c_fwd, data_real_in1d, data_cmplx_out1d)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
            call fftw_execute_dft(plan_y_c2c_fwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_cmplx_slabz_batch(l, i, j, k) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        ! pressure: X fwd FFT, Y fwd FFT
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    data_real_3D_slabz(i + 1, j + 1, k + 1) = pressure%sf(i, j, k)*fluid_indicator_function%sf(i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nzloc
                    data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
        ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, data_cmplx_out1d)
#else
        call fftw_execute_dft_r2c(plan_x_r2c_fwd, data_real_in1d, data_cmplx_out1d)
#endif
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_y_c2c_fwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_slabz_batch(sys_size + 1, i, j, k) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! reynolds stress: X fwd FFT, Y fwd FFT
        do l = 1, num_dims
            do q = 1, num_dims
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 0, m
                    do j = 0, n
                        do k = 0, p
                            data_real_3D_slabz(i + 1, j + 1, k + 1) = reynolds_stress(l, q)%sf(i, j, k)*fluid_indicator_function%sf(i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, Nx
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, data_cmplx_out1d)
#else
                call fftw_execute_dft_r2c(plan_x_r2c_fwd, data_real_in1d, data_cmplx_out1d)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
                call fftw_execute_dft(plan_y_c2c_fwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_cmplx_slabz_batch(reynolds_stress_idx + 3*(l - 1) + q, i, j, k) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

        ! effective viscosity: X fwd FFT, Y fwd FFT
        do l = 1, num_dims
            do q = 1, num_dims
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 0, m
                    do j = 0, n
                        do k = 0, p
                            data_real_3D_slabz(i + 1, j + 1, k + 1) = visc_stress(l, q)%sf(i, j, k)*fluid_indicator_function%sf(i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, Nx
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, data_cmplx_out1d)
#else
                call fftw_execute_dft_r2c(plan_x_r2c_fwd, data_real_in1d, data_cmplx_out1d)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
                call fftw_execute_dft(plan_y_c2c_fwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_cmplx_slabz_batch(eff_visc_idx + 3*(l - 1) + q, i, j, k) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

        ! interphase momentum exchange: X fwd FFT, Y fwd FFT
        do l = 1, num_dims
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        data_real_3D_slabz(i + 1, j + 1, k + 1) = pres_visc_stress(l, 1)%sf(i, j, k)*grad_fluid_indicator(1)%sf(i, j, k) &
                                                                  + pres_visc_stress(l, 2)%sf(i, j, k)*grad_fluid_indicator(2)%sf(i, j, k) &
                                                                  + pres_visc_stress(l, 3)%sf(i, j, k)*grad_fluid_indicator(3)%sf(i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, Nx
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, data_cmplx_out1d)
#else
            call fftw_execute_dft_r2c(plan_x_r2c_fwd, data_real_in1d, data_cmplx_out1d)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
            call fftw_execute_dft(plan_y_c2c_fwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_cmplx_slabz_batch(int_mom_exch_idx + l, i, j, k) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        call s_mpi_transpose_slabZ2Y_batch

        ! cons vars: Z fwd FFT, convolution, Z bwd FFT
        do l = 1, sys_size
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Nyloc
                    do k = 1, Nz
                        data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_slaby_batch(l, i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_FORWARD)
#else
            call fftw_execute_dft(plan_z_c2c_fwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Nyloc
                    do k = 1, Nz
                        data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)*cmplx_kernelG1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_INVERSE)
#else
            call fftw_execute_dft(plan_z_c2c_bwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Nyloc
                    do k = 1, Nz
                        data_cmplx_slaby_batch(l, i, j, k) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        ! pressure: Z fwd FFT, convolution, Z bwd FFT
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Nyloc
                do k = 1, Nz
                    data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_slaby_batch(sys_size + 1, i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_z_c2c_fwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Nyloc
                do k = 1, Nz
                    data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)*cmplx_kernelG1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_INVERSE)
#else
        call fftw_execute_dft(plan_z_c2c_bwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Nyloc
                do k = 1, Nz
                    data_cmplx_slaby_batch(sys_size + 1, i, j, k) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! reynolds stress: Z fwd FFT, convolution, Z bwd FFT
        do l = 1, num_dims
            do q = 1, num_dims
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Nyloc
                        do k = 1, Nz
                            data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_slaby_batch(reynolds_stress_idx + 3*(l - 1) + q, i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_FORWARD)
#else
                call fftw_execute_dft(plan_z_c2c_fwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Nyloc
                        do k = 1, Nz
                            data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)*cmplx_kernelG1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_INVERSE)
#else
                call fftw_execute_dft(plan_z_c2c_bwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Nyloc
                        do k = 1, Nz
                            data_cmplx_slaby_batch(reynolds_stress_idx + 3*(l - 1) + q, i, j, k) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

        ! effective viscosity: Z fwd FFT, convolution, Z bwd FFT
        do l = 1, num_dims
            do q = 1, num_dims
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Nyloc
                        do k = 1, Nz
                            data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_slaby_batch(eff_visc_idx + 3*(l - 1) + q, i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_FORWARD)
#else
                call fftw_execute_dft(plan_z_c2c_fwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Nyloc
                        do k = 1, Nz
                            data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)*cmplx_kernelG1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_INVERSE)
#else
                call fftw_execute_dft(plan_z_c2c_bwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Nyloc
                        do k = 1, Nz
                            data_cmplx_slaby_batch(eff_visc_idx + 3*(l - 1) + q, i, j, k) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

        ! interphase momentum exchange: Z fwd FFT, convolution, Z bwd FFT
        do l = 1, num_dims
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Nyloc
                    do k = 1, Nz
                        data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_slaby_batch(int_mom_exch_idx + l, i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_FORWARD)
#else
            call fftw_execute_dft(plan_z_c2c_fwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Nyloc
                    do k = 1, Nz
                        data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)*cmplx_kernelG1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_INVERSE)
#else
            call fftw_execute_dft(plan_z_c2c_bwd, data_cmplx_out1d, data_cmplx_out1d)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Nyloc
                    do k = 1, Nz
                        data_cmplx_slaby_batch(int_mom_exch_idx + l, i, j, k) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        call s_mpi_transpose_slabY2Z_batch

        ! cons vars: Y bwd FFT, X bwd FFT
        do l = 1, sys_size
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_slabz_batch(l, i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_INVERSE)
#else
            call fftw_execute_dft(plan_y_c2c_bwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2D(plan_x_bwd_gpu, data_cmplx_out1d, data_real_in1d)
#else
            call fftw_execute_dft_c2r(plan_x_c2r_bwd, data_cmplx_out1d, data_real_in1d)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, Nx
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_real_3D_slabz(i, j, k) = data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        q_cons_filtered(l)%sf(i, j, k) = data_real_3D_slabz(i + 1, j + 1, k + 1)*fft_norm/filtered_fluid_indicator_function%sf(i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        ! pressure: Y bwd FFT, X bwd FFT
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_slabz_batch(sys_size + 1, i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_INVERSE)
#else
        call fftw_execute_dft(plan_y_c2c_bwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2D(plan_x_bwd_gpu, data_cmplx_out1d, data_real_in1d)
#else
        call fftw_execute_dft_c2r(plan_x_c2r_bwd, data_cmplx_out1d, data_real_in1d)
#endif
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nzloc
                    data_real_3D_slabz(i, j, k) = data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    filtered_pressure%sf(i, j, k) = data_real_3D_slabz(i + 1, j + 1, k + 1)*fft_norm/filtered_fluid_indicator_function%sf(i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! reynolds stress: Y bwd FFT, X bwd FFT
        do l = 1, num_dims
            do q = 1, num_dims
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_slabz_batch(reynolds_stress_idx + 3*(l - 1) + q, i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_INVERSE)
#else
                call fftw_execute_dft(plan_y_c2c_bwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2D(plan_x_bwd_gpu, data_cmplx_out1d, data_real_in1d)
#else
                call fftw_execute_dft_c2r(plan_x_c2r_bwd, data_cmplx_out1d, data_real_in1d)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, Nx
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_real_3D_slabz(i, j, k) = data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 0, m
                    do j = 0, n
                        do k = 0, p
                            reynolds_stress(l, q)%sf(i, j, k) = data_real_3D_slabz(i + 1, j + 1, k + 1)*fft_norm/filtered_fluid_indicator_function%sf(i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

        ! effective viscosity: Y bwd FFT, X bwd FFT
        do l = 1, num_dims
            do q = 1, num_dims
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_slabz_batch(eff_visc_idx + 3*(l - 1) + q, i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_INVERSE)
#else
                call fftw_execute_dft(plan_y_c2c_bwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, NxC
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
                ierr = cufftExecZ2D(plan_x_bwd_gpu, data_cmplx_out1d, data_real_in1d)
#else
                call fftw_execute_dft_c2r(plan_x_c2r_bwd, data_cmplx_out1d, data_real_in1d)
#endif
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, Nx
                    do j = 1, Ny
                        do k = 1, Nzloc
                            data_real_3D_slabz(i, j, k) = data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 0, m
                    do j = 0, n
                        do k = 0, p
                            eff_visc(l, q)%sf(i, j, k) = data_real_3D_slabz(i + 1, j + 1, k + 1)*fft_norm/filtered_fluid_indicator_function%sf(i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

        ! interphase momentum exchange: Y bwd FFT, X bwd FFT
        do l = 1, num_dims
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_slabz_batch(int_mom_exch_idx + l, i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_INVERSE)
#else
            call fftw_execute_dft(plan_y_c2c_bwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, NxC
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_OpenACC)
            ierr = cufftExecZ2D(plan_x_bwd_gpu, data_cmplx_out1d, data_real_in1d)
#else
            call fftw_execute_dft_c2r(plan_x_c2r_bwd, data_cmplx_out1d, data_real_in1d)
#endif
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 1, Nx
                do j = 1, Ny
                    do k = 1, Nzloc
                        data_real_3D_slabz(i, j, k) = data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            $:GPU_PARALLEL_LOOP(collapse=3)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        int_mom_exch(l)%sf(i, j, k) = data_real_3D_slabz(i + 1, j + 1, k + 1)*fft_norm
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_filter_batch

    !< compute forward FFT, input: data_real_3D_slabz, output: data_cmplx_out1d
    subroutine s_mpi_FFT_fwd
        integer :: i, j, k

        ! 3D z-slab -> 1D x, y, z
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nzloc
                    data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! X FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, data_cmplx_out1d)
#else
        call fftw_execute_dft_r2c(plan_x_r2c_fwd, data_real_in1d, data_cmplx_out1d)
#endif

        ! 1D x, y, z -> 1D y, x, z (CMPLX)
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Y FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_y_c2c_fwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif

        ! 1D y, x, z -> 3D z-slab
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_slabz(i, j, k) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! transpose z-slab to y-slab
        call nvtxStartRange("SLAB-MPI-TRANSPOSE-Z2Y")
        call s_mpi_transpose_slabZ2Y
        call nvtxEndRange

        ! 3D y-slab -> 1D z, x, y
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Nyloc
                do k = 1, Nz
                    data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC) = data_cmplx_slaby(i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

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
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Nyloc
                do k = 1, Nz
                    data_cmplx_slaby(i, j, k) = data_cmplx_out1d(k + (i - 1)*Nz + (j - 1)*Nz*NxC)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! transpose y-slab to z-slab
        call nvtxStartRange("SLAB-MPI-TRANSPOSE-Y2Z")
        call s_mpi_transpose_slabY2Z
        call nvtxEndRange

        ! 3D z-slab -> 1D y, x, z
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC) = data_cmplx_slabz(i, j, k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Y inv FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_INVERSE)
#else
        call fftw_execute_dft(plan_y_c2c_bwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif

        ! 1D y, x, z -> 1D x, y, z
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, NxC
            do j = 1, Ny
                do k = 1, Nzloc
                    data_cmplx_out1d(i + (j - 1)*NxC + (k - 1)*NxC*Ny) = data_cmplx_out1dy(j + (i - 1)*Ny + (k - 1)*Ny*NxC)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! X inv FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2D(plan_x_bwd_gpu, data_cmplx_out1d, data_real_in1d)
#else
        call fftw_execute_dft_c2r(plan_x_c2r_bwd, data_cmplx_out1d, data_real_in1d)
#endif

        ! 1D x, y, z -> 3D z-slab
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nzloc
                    data_real_3D_slabz(i, j, k) = data_real_in1d(i + (j - 1)*Nx + (k - 1)*Nx*Ny)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mpi_FFT_bwd

    subroutine s_finalize_fftw_explicit_filter_module
        integer :: i, j

        @:DEALLOCATE(fluid_indicator_function%sf)
        @:DEALLOCATE(filtered_fluid_indicator_function%sf)
        do i = 1, num_dims
            @:DEALLOCATE(grad_fluid_indicator(i)%sf)
        end do
        @:DEALLOCATE(grad_fluid_indicator)

        do i = 1, sys_size
            @:DEALLOCATE(q_cons_filtered(i)%sf)
        end do
        @:DEALLOCATE(q_cons_filtered)

        do i = 1, sys_size
            @:DEALLOCATE(q_prim_filtered(i)%sf)
        end do
        @:DEALLOCATE(q_prim_filtered)

        @:DEALLOCATE(filtered_pressure%sf)

        do i = 1, num_dims
            do j = 1, num_dims
                @:DEALLOCATE(visc_stress(i, j)%sf)
            end do
        end do
        @:DEALLOCATE(visc_stress)

        do i = 1, num_dims
            do j = 1, num_dims
                @:DEALLOCATE(pres_visc_stress(i, j)%sf)
            end do
        end do
        @:DEALLOCATE(pres_visc_stress)

        do i = 1, num_dims
            @:DEALLOCATE(div_pres_visc_stress(i)%sf)
        end do
        @:DEALLOCATE(div_pres_visc_stress)

        do i = 1, num_dims
            do j = 1, num_dims
                @:DEALLOCATE(reynolds_stress(i, j)%sf)
            end do
        end do
        @:DEALLOCATE(reynolds_stress)

        do i = 1, num_dims
            do j = 1, num_dims
                @:DEALLOCATE(eff_visc(i, j)%sf)
            end do
        end do
        @:DEALLOCATE(eff_visc)

        do i = 1, num_dims
            @:DEALLOCATE(int_mom_exch(i)%sf)
        end do
        @:DEALLOCATE(int_mom_exch)

        @:DEALLOCATE(particle_forces)

        @:DEALLOCATE(data_real_in1d, data_cmplx_out1d, data_cmplx_out1dy)
        @:DEALLOCATE(cmplx_kernelG1d, real_kernelG_in)
        @:DEALLOCATE(data_real_3D_slabz, data_cmplx_slabz, data_cmplx_slaby)
        @:DEALLOCATE(data_cmplx_slabz_batch, data_cmplx_slaby_batch)

        @:DEALLOCATE(sendbuf_sf, recvbuf_sf)
        @:DEALLOCATE(sendbuf_batch, recvbuf_batch)

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

        if (compute_particle_drag) then
            if (proc_rank == 0) then
                close (100)
            end if
        end if

    end subroutine s_finalize_fftw_explicit_filter_module

end module m_volume_filtering
