#:include 'case.fpp'
#:include 'macros.fpp'

module m_igr

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters

    use m_variables_conversion

    use m_mpi_proxy

    use m_helper

    use m_boundary_common

    implicit none

    private; public :: s_initialize_igr_module, &
 s_igr_iterative_solve, &
 s_igr_riemann_solver, &
 s_igr_sigma_x, &
 s_igr_flux_add, &
 s_finalize_igr_module

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
    integer, dimension(3) :: nv_uvm_temp_on_gpu
    real(wp), pointer, contiguous, dimension(:, :, :) :: jac, jac_rhs, jac_old
    real(wp), allocatable, dimension(:, :, :), pinned, target :: jac_host
    real(wp), allocatable, dimension(:, :, :), pinned, target :: jac_rhs_host
    real(wp), allocatable, dimension(:, :, :), pinned, target :: jac_old_host
#else
    real(wp), allocatable, target, dimension(:, :, :) :: jac
    real(wp), allocatable, dimension(:, :, :) :: jac_rhs, jac_old
    $:GPU_DECLARE(create='[jac, jac_rhs, jac_old]')
#endif
    type(scalar_field), dimension(1) :: jac_sf
    $:GPU_DECLARE(create='[jac_sf]')

    real(wp), allocatable, dimension(:, :) :: Res_igr
    $:GPU_DECLARE(create='[Res_igr]')

    real(wp) :: alf_igr
    $:GPU_DECLARE(create='[alf_igr]')

    #:if not MFC_CASE_OPTIMIZATION
        integer :: vidxb, vidxe
        $:GPU_DECLARE(create='[vidxb, vidxe]')

        integer :: offxL, offxR
        $:GPU_DECLARE(create='[offxL, offxR]')

        real(wp), allocatable, dimension(:) :: coeff_L, coeff_R
        $:GPU_DECLARE(create='[coeff_L, coeff_R]')
    #:else
        #:if igr_order == 5
            integer, parameter :: vidxb = -2
            integer, parameter :: vidxe = 3

            integer, parameter :: offxL = 2
            integer, parameter :: offxR = 3
#if defined(MFC_OpenMP)
            real(wp) :: coeff_L(1:5) = [ &
                        -3._wp/60._wp, &  ! Index -1
                        27._wp/60._wp, &  ! Index 0
                        47._wp/60._wp, &  ! Index 1
                        -13._wp/60._wp, &  ! Index 2
                        2._wp/60._wp &  ! Index 3
                        ]

            real(wp) :: coeff_R(1:5) = [ &
                        2._wp/60._wp, &  ! Index -2
                        -13._wp/60._wp, &  ! Index -1
                        47._wp/60._wp, &  ! Index 0
                        27._wp/60._wp, &  ! Index 1
                        -3._wp/60._wp &  ! Index 2
                        ]
#else
            real(wp), parameter :: coeff_L(1:5) = [ &
                                   -3._wp/60._wp, &  ! Index -1
                                   27._wp/60._wp, &  ! Index 0
                                   47._wp/60._wp, &  ! Index 1
                                   -13._wp/60._wp, &  ! Index 2
                                   2._wp/60._wp &  ! Index 3
                                   ]

            real(wp), parameter :: coeff_R(1:5) = [ &
                                   2._wp/60._wp, &  ! Index -2
                                   -13._wp/60._wp, &  ! Index -1
                                   47._wp/60._wp, &  ! Index 0
                                   27._wp/60._wp, &  ! Index 1
                                   -3._wp/60._wp &  ! Index 2
                                   ]
#endif
        #:elif igr_order == 3
            integer, parameter :: vidxb = -1
            integer, parameter :: vidxe = 2

            integer, parameter :: offxL = 1
            integer, parameter :: offxR = 2
#if defined(MFC_OpenMP)
            real(wp) :: coeff_L(1:3) = [ &
                        2._wp/6._wp, & ! Index 0
                        5._wp/6._wp, & ! Index 1
                        -1._wp/6._wp & ! Index 2
                        ]
            real(wp) :: coeff_R(1:3) = [ &
                        -1._wp/6._wp, & ! Index -1
                        5._wp/6._wp, & ! Index 0
                        2._wp/6._wp & ! Index 1
                        ]
#else
            real(wp), parameter :: coeff_L(1:3) = [ &
                                   2._wp/6._wp, & ! Index 0
                                   5._wp/6._wp, & ! Index 1
                                   -1._wp/6._wp & ! Index 2
                                   ]
            real(wp), parameter :: coeff_R(1:3) = [ &
                                   -1._wp/6._wp, & ! Index -1
                                   5._wp/6._wp, & ! Index 0
                                   2._wp/6._wp & ! Index 1
                                   ]
#endif
        #:endif

#if defined(MFC_OpenMP)
        $:GPU_DECLARE(create='[coeff_L, coeff_R]')
#endif
    #:endif

    integer(kind=8) :: i, j, k, l, q, r

contains

    subroutine s_initialize_igr_module()

        if (viscous) then
            @:ALLOCATE(Res_igr(1:2, 1:maxval(Re_size)))
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res_igr(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            $:GPU_UPDATE(device='[Res_igr, Re_idx, Re_size]')
            @:PREFER_GPU(Res_igr)
            @:PREFER_GPU(Re_idx)
        end if

#ifndef __NVCOMPILER_GPU_UNIFIED_MEM
        @:ALLOCATE(jac(idwbuff(1)%beg:idwbuff(1)%end, &
            idwbuff(2)%beg:idwbuff(2)%end, &
            idwbuff(3)%beg:idwbuff(3)%end))
        @:ALLOCATE(jac_rhs(-1:m,-1:n,-1:p))

        if (igr_iter_solver == 1) then ! Jacobi iteration
            @:ALLOCATE(jac_old(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
        end if
#else
        ! create map
        nv_uvm_temp_on_gpu(1:3) = 0
        nv_uvm_temp_on_gpu(1:nv_uvm_igr_temps_on_gpu) = 1

        if (nv_uvm_temp_on_gpu(1) == 1) then
            @:ALLOCATE(jac(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:PREFER_GPU(jac)
        else
            allocate (jac_host(idwbuff(1)%beg:idwbuff(1)%end, &
                               idwbuff(2)%beg:idwbuff(2)%end, &
                               idwbuff(3)%beg:idwbuff(3)%end))

            jac(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end) => jac_host(:, :, :)
        end if

        if (nv_uvm_temp_on_gpu(2) == 1) then
            @:ALLOCATE(jac_rhs(-1:m,-1:n,-1:p))
            @:PREFER_GPU(jac_rhs)
        else
            allocate (jac_rhs_host(-1:m, -1:n, -1:p))
            jac_rhs(-1:m, -1:n, -1:p) => jac_rhs_host(:, :, :)
        end if

        if (igr_iter_solver == 1) then ! Jacobi iteration
            if (nv_uvm_temp_on_gpu(3) == 1) then
                @:ALLOCATE(jac_old(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:PREFER_GPU(jac_old)
            else
                allocate (jac_old_host(idwbuff(1)%beg:idwbuff(1)%end, &
                                       idwbuff(2)%beg:idwbuff(2)%end, &
                                       idwbuff(3)%beg:idwbuff(3)%end))

                jac_old(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end) => jac_old_host(:, :, :)
            end if
        end if
#endif

        $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    jac(j, k, l) = 0._stp
                    if (igr_iter_solver == 1) jac_old(j, k, l) = 0._stp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (p == 0) then
            alf_igr = alf_factor*max(dx(1), dy(1))**2._wp
        else
            alf_igr = alf_factor*max(dx(1), dy(1), dz(1))**2._wp
        end if
        $:GPU_UPDATE(device='[alf_igr]')

        #:if not MFC_CASE_OPTIMIZATION
            if (igr_order == 3) then
                vidxb = -1; vidxe = 2; offxL = 1; offxR = 2
                $:GPU_UPDATE(device='[vidxb, vidxe, offxL, offxR]')

                @:ALLOCATE(coeff_L(1:3))
                coeff_L(1) = (2._wp/6._wp)
                coeff_L(2) = (5._wp/6._wp)
                coeff_L(3) = (-1._wp/6._wp)

                @:ALLOCATE(coeff_R(1:3))
                coeff_R(3) = (2._wp/6._wp)
                coeff_R(2) = (5._wp/6._wp)
                coeff_R(1) = (-1._wp/6._wp)

            elseif (igr_order == 5) then
                vidxb = -2; vidxe = 3; offxL = 2; offxR = 3
                $:GPU_UPDATE(device='[vidxb, vidxe, offxL, offxR]')

                @:ALLOCATE(coeff_L(1:5))
                coeff_L(1) = (-3._wp/60._wp)
                coeff_L(2) = (27._wp/60._wp)
                coeff_L(3) = (47._wp/60._wp)
                coeff_L(4) = (-13._wp/60._wp)
                coeff_L(5) = (2._wp/60._wp)

                @:ALLOCATE(coeff_R(1:5))
                coeff_R(5) = (-3._wp/60._wp)
                coeff_R(4) = (27._wp/60._wp)
                coeff_R(3) = (47._wp/60._wp)
                coeff_R(2) = (-13._wp/60._wp)
                coeff_R(1) = (2._wp/60._wp)
            end if

            $:GPU_UPDATE(device='[coeff_L]')
            $:GPU_UPDATE(device='[coeff_R]')
        #:else
#if defined(MFC_OpenMP)
            $:GPU_UPDATE(device='[coeff_L]')
            $:GPU_UPDATE(device='[coeff_R]')
#endif
        #:endif
        jac_sf(1)%sf => jac
        $:GPU_ENTER_DATA(copyin='[jac_sf(1)%sf]')
        $:GPU_ENTER_DATA(attach='[jac_sf(1)%sf]')

    end subroutine s_initialize_igr_module

    subroutine s_igr_iterative_solve(q_cons_vf, bc_type, t_step)
#ifdef _CRAYFTN
        !DIR$ OPTIMIZE (-haggress)
#endif
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type
        integer, intent(in) :: t_step

        real(wp) :: rho_rx, rho_ry, rho_rz, rho_lx, rho_ly, rho_lz
        real(wp) :: fd_coeff
        integer :: num_iters

        if (t_step == t_step_start) then
            num_iters = num_igr_warm_start_iters
        else
            num_iters = num_igr_iters
        end if

        do q = 1, num_iters
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j,k,l,rho_lx, rho_rx, rho_ly, rho_ry, rho_lz, rho_rz, fd_coeff]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rho_lx = 0._wp
                        rho_rx = 0._wp
                        rho_ly = 0._wp
                        rho_ry = 0._wp
                        rho_lz = 0._wp
                        rho_rz = 0._wp
                        fd_coeff = 0._wp

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            rho_lx = rho_lx + real(q_cons_vf(i)%sf(j, k, l) + q_cons_vf(i)%sf(j - 1, k, l), kind=wp)/2._wp
                            rho_rx = rho_rx + real(q_cons_vf(i)%sf(j, k, l) + q_cons_vf(i)%sf(j + 1, k, l), kind=wp)/2._wp
                            rho_ly = rho_ly + real(q_cons_vf(i)%sf(j, k, l) + q_cons_vf(i)%sf(j, k - 1, l), kind=wp)/2._wp
                            rho_ry = rho_ry + real(q_cons_vf(i)%sf(j, k, l) + q_cons_vf(i)%sf(j, k + 1, l), kind=wp)/2._wp
                            if (p > 0) then
                                rho_lz = rho_lz + real(q_cons_vf(i)%sf(j, k, l) + q_cons_vf(i)%sf(j, k, l - 1), kind=wp)/2._wp
                                rho_rz = rho_rz + real(q_cons_vf(i)%sf(j, k, l) + q_cons_vf(i)%sf(j, k, l + 1), kind=wp)/2._wp
                            end if
                            fd_coeff = fd_coeff + q_cons_vf(i)%sf(j, k, l)
                        end do

                        fd_coeff = 1._wp/fd_coeff + alf_igr* &
                                   ((1._wp/dx(j)**2._wp)*(1._wp/rho_lx + 1._wp/rho_rx) + &
                                    (1._wp/dy(k)**2._wp)*(1._wp/rho_ly + 1._wp/rho_ry))

                        if (num_dims == 3) then
                            fd_coeff = fd_coeff + alf_igr*(1._wp/dz(l)**2._wp)*(1._wp/rho_lz + 1._wp/rho_rz)
                        end if

                        if (igr_iter_solver == 1) then ! Jacobi iteration
                            if (num_dims == 3) then
                                jac(j, k, l) = real((alf_igr/fd_coeff)* &
                                                    ((1._wp/dx(j)**2._wp)*(jac_old(j - 1, k, l)/rho_lx + jac_old(j + 1, k, l)/rho_rx) + &
                                                     (1._wp/dy(k)**2._wp)*(jac_old(j, k - 1, l)/rho_ly + jac_old(j, k + 1, l)/rho_ry) + &
                                                     (1._wp/dz(l)**2._wp)*(jac_old(j, k, l - 1)/rho_lz + jac_old(j, k, l + 1)/rho_rz)) + &
                                                    real(jac_rhs(j, k, l), kind=wp)/fd_coeff, kind=stp)
                            else
                                jac(j, k, l) = real((alf_igr/fd_coeff)* &
                                                    ((1._wp/dx(j)**2._wp)*(real(jac_old(j - 1, k, l), kind=wp)/rho_lx + real(jac_old(j + 1, k, l), kind=wp)/rho_rx) + &
                                                     (1._wp/dy(k)**2._wp)*(real(jac_old(j, k - 1, l), kind=wp)/rho_ly + real(jac_old(j, k + 1, l), kind=wp)/rho_ry)) + &
                                                    real(jac_rhs(j, k, l), kind=wp)/fd_coeff, kind=stp)
                            end if
                        else ! Gauss Seidel iteration
                            if (num_dims == 3) then
                                jac(j, k, l) = real((alf_igr/fd_coeff)* &
                                                    ((1._wp/dx(j)**2._wp)*(jac(j - 1, k, l)/rho_lx + jac(j + 1, k, l)/rho_rx) + &
                                                     (1._wp/dy(k)**2._wp)*(jac(j, k - 1, l)/rho_ly + jac(j, k + 1, l)/rho_ry) + &
                                                     (1._wp/dz(l)**2._wp)*(jac(j, k, l - 1)/rho_lz + jac(j, k, l + 1)/rho_rz)) + &
                                                    real(jac_rhs(j, k, l), kind=wp)/fd_coeff, kind=stp)
                            else
                                jac(j, k, l) = real((alf_igr/fd_coeff)* &
                                                    ((1._wp/dx(j)**2._wp)*(jac(j - 1, k, l)/rho_lx + jac(j + 1, k, l)/rho_rx) + &
                                                     (1._wp/dy(k)**2._wp)*(jac(j, k - 1, l)/rho_ly + jac(j, k + 1, l)/rho_ry)) + &
                                                    real(jac_rhs(j, k, l), kind=wp)/fd_coeff, kind=stp)
                            end if
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            call s_populate_F_igr_buffers(bc_type, jac_sf)

            if (igr_iter_solver == 1) then ! Jacobi iteration
                $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            jac_old(j, k, l) = jac(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end do

    end subroutine s_igr_iterative_solve

    subroutine s_igr_sigma_x(q_cons_vf, rhs_vf)
#ifdef _CRAYFTN
        !DIR$ OPTIMIZE (-haggress)
#endif
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: rhs_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_cons_vf

        real(wp) :: F_L, vel_L, rho_L, F_R, vel_R, rho_R
        #:if not MFC_CASE_OPTIMIZATION
            real(wp), dimension(num_fluids_max) :: alpha_rho_L, alpha_rho_R
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho_L, alpha_rho_R
        #:endif

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j,k,l,F_L, vel_L, alpha_rho_L, F_R, vel_R, alpha_rho_R, rho_L, rho_R]')
        do l = 0, p
            do k = 0, n
                do j = -1, m

                    F_L = 0._wp; F_R = 0._wp
                    vel_L = 0._wp; vel_R = 0._wp
                    rho_L = 0._wp; rho_R = 0._wp

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_rho_L(i) = 0._wp
                        alpha_rho_R(i) = 0._wp
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = vidxb + 1, vidxe
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            alpha_rho_L(i) = alpha_rho_L(i) + coeff_L(q + offxL)*q_cons_vf(i)%sf(j + q, k, l)
                        end do

                        vel_L = vel_L + coeff_L(q + offxL)*q_cons_vf(momxb)%sf(j + q, k, l)
                        F_L = F_L + coeff_L(q + offxL)*jac(j + q, k, l)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = vidxb, vidxe - 1
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            alpha_rho_R(i) = alpha_rho_R(i) + coeff_R(q + offxR)*q_cons_vf(i)%sf(j + q, k, l)
                        end do

                        vel_R = vel_R + coeff_R(q + offxR)*q_cons_vf(momxb)%sf(j + q, k, l)
                        F_R = F_R + coeff_R(q + offxR)*jac(j + q, k, l)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        rho_L = rho_L + alpha_rho_L(i)
                        rho_R = rho_R + alpha_rho_R(i)
                    end do

                    vel_L = vel_L/rho_L
                    vel_R = vel_R/rho_R

                    #:for LR in ['L', 'R']
                        $:GPU_ATOMIC(atomic='update')
                        rhs_vf(momxb)%sf(j + 1, k, l) = rhs_vf(momxb)%sf(j + 1, k, l) + &
                                                        real(0.5_wp*dt*F_${LR}$*(1._wp/dx(j + 1)), kind=stp)
                        $:GPU_ATOMIC(atomic='update')
                        rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) + &
                                                        real(0.5_wp*dt*vel_${LR}$*F_${LR}$*(1._wp/dx(j + 1)), kind=stp)
                        $:GPU_ATOMIC(atomic='update')
                        rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                    real(0.5_wp*dt*F_${LR}$*(1._wp/dx(j)), kind=stp)
                        $:GPU_ATOMIC(atomic='update')
                        rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                    real(0.5_wp*dt*vel_${LR}$*F_${LR}$*(1._wp/dx(j)), kind=stp)
                    #:endfor
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_igr_sigma_x

    subroutine s_igr_riemann_solver(q_cons_vf, rhs_vf, idir)
#ifdef _CRAYFTN
        !DIR$ OPTIMIZE (-haggress)
#endif
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: rhs_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_cons_vf
        integer, intent(in) :: idir

        real(wp) :: cfl
        real(wp) :: rho_L, gamma_L, pi_inf_L, E_L, mu_L, F_L, pres_L
        real(wp) :: rho_R, gamma_R, pi_inf_R, E_R, mu_R, F_R, pres_R
        real(wp), dimension(3) :: vflux_L_arr, vflux_R_arr
        real(wp), dimension(-1:1) :: rho_sf_small
        #:if not MFC_CASE_OPTIMIZATION
            real(wp), dimension(num_fluids_max) :: alpha_rho_L, alpha_L, alpha_R, alpha_rho_R
            real(wp), dimension(3) :: vel_L, vel_R
            real(wp), dimension(3, 3) :: dvel
            real(wp), dimension(3) :: dvel_small
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho_L, alpha_L, alpha_R, alpha_rho_R
            real(wp), dimension(num_dims) :: vel_L, vel_R
            real(wp), dimension(num_dims, num_dims) :: dvel
            real(wp), dimension(num_dims) :: dvel_small
        #:endif

        if (idir == 1) then
            if (p == 0) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j,k,l,rho_L, rho_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, mu_L, mu_R, vel_L, vel_R, pres_L, pres_R, alpha_L, alpha_R, alpha_rho_L, alpha_rho_R, F_L, F_R, E_L, E_R, cfl, dvel, dvel_small, rho_sf_small, vflux_L_arr, vflux_R_arr]')
                do l = 0, p
                    do k = 0, n
                        do j = -1, m

                            vflux_L_arr = 0._wp
                            vflux_R_arr = 0._wp

                            #:if MFC_CASE_OPTIMIZATION
                                #:if igr_order == 5
                                    !DIR$ unroll 6
                                #:elif igr_order == 3
                                    !DIR$ unroll 4
                                #:endif
                            #:endif
                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb, vidxe
                                !x-direction contributions
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = -1, 1
                                    rho_L = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do r = 1, num_fluids
                                        rho_L = rho_L + q_cons_vf(r)%sf(j + i + q, k, l)
                                    end do
                                    rho_sf_small(i) = rho_L
                                end do

                                dvel_small(1) = (1/(2._wp*dx(j)))*( &
                                                1._wp*q_cons_vf(momxb)%sf(j + 1 + q, k, l)/rho_sf_small(1) - &
                                                1._wp*q_cons_vf(momxb)%sf(j - 1 + q, k, l)/rho_sf_small(-1))
                                dvel_small(2) = (1/(2._wp*dx(j)))*( &
                                                q_cons_vf(momxb + 1)%sf(j + 1 + q, k, l)/rho_sf_small(1) - &
                                                q_cons_vf(momxb + 1)%sf(j - 1 + q, k, l)/rho_sf_small(-1))

                                if (q == 0) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        dvel(i, 1) = dvel_small(i)
                                    end do
                                end if

                                if (q > vidxb) then
                                    vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(2))
                                    vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(4._wp*dvel_small(1))/3._wp
                                end if
                                if (q < vidxe) then
                                    vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(2))
                                    vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(4._wp*dvel_small(1))/3._wp
                                end if

                                !y-direction contributions
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = -1, 1
                                    rho_L = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do r = 1, num_fluids
                                        rho_L = rho_L + q_cons_vf(r)%sf(j + q, k + i, l)
                                    end do
                                    rho_sf_small(i) = rho_L
                                end do

                                dvel_small(1) = (1/(2._wp*dy(k)))*( &
                                                q_cons_vf(momxb)%sf(j + q, k + 1, l)/rho_sf_small(1) - &
                                                q_cons_vf(momxb)%sf(j + q, k - 1, l)/rho_sf_small(-1))
                                dvel_small(2) = (1/(2._wp*dy(k)))*( &
                                                q_cons_vf(momxb + 1)%sf(j + q, k + 1, l)/rho_sf_small(1) - &
                                                q_cons_vf(momxb + 1)%sf(j + q, k - 1, l)/rho_sf_small(-1))

                                if (q == 0) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        dvel(i, 2) = dvel_small(i)
                                    end do
                                end if

                                if (q > vidxb) then
                                    vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(1))
                                    vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(-2._wp*dvel_small(2))/3._wp
                                end if
                                if (q < vidxe) then
                                    vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(1))
                                    vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(-2._wp*dvel_small(2))/3._wp
                                end if

                                if (q == 0) then
                                    jac_rhs(j, k, l) = real(alf_igr*(2._wp*(dvel(1, 2)*dvel(2, 1)) &
                                                                     + dvel(1, 1)**2._wp + dvel(2, 2)**2._wp &
                                                                     + (dvel(1, 1) + dvel(2, 2))**2._wp), kind=stp)
                                end if
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_L(i) = 0._wp
                                alpha_rho_R(i) = 0._wp
                                alpha_L(i) = 0._wp
                                alpha_R(i) = 0._wp
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = 0._wp
                                vel_R(i) = 0._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb + 1, vidxe
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = alpha_rho_L(i) + coeff_L(q + offxL)*q_cons_vf(i)%sf(j + q, k, l)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        alpha_L(i) = alpha_L(i) + coeff_L(q + offxL)*q_cons_vf(E_idx + i)%sf(j + q, k, l)
                                    end do
                                else
                                    alpha_L(1) = 1._wp
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = vel_L(i) + coeff_L(q + offxL)*q_cons_vf(momxb + i - 1)%sf(j + q, k, l)
                                end do
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb, vidxe - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_R(i) = alpha_rho_R(i) + coeff_R(q + offxR)*q_cons_vf(i)%sf(j + q, k, l)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        alpha_R(i) = alpha_R(i) + coeff_R(q + offxR)*q_cons_vf(E_idx + i)%sf(j + q, k, l)
                                    end do
                                else
                                    alpha_R(1) = 1._wp
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_R(i) = vel_R(i) + coeff_R(q + offxR)*q_cons_vf(momxb + i - 1)%sf(j + q, k, l)
                                end do
                            end do

                            if (num_fluids > 1) then

                                alpha_L(num_fluids) = 1._wp
                                alpha_R(num_fluids) = 1._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    alpha_L(num_fluids) = alpha_L(num_fluids) - alpha_L(i)
                                    alpha_R(num_fluids) = alpha_R(num_fluids) - alpha_R(i)
                                end do
                            end if

                            rho_L = 0._wp; rho_R = 0._wp
                            gamma_L = 0._wp; gamma_R = 0._wp
                            pi_inf_L = 0._wp; pi_inf_R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho_L = rho_L + alpha_rho_L(i)
                                gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)

                                rho_R = rho_R + alpha_rho_R(i)
                                gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = vel_L(i)/rho_L
                                vel_R(i) = vel_R(i)/rho_R
                            end do

                            if (viscous) then
                                mu_L = 0._wp; mu_R = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    mu_L = alpha_L(i)/Res_igr(1, i) + mu_L
                                    mu_R = alpha_R(i)/Res_igr(1, i) + mu_R
                                end do

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j + 1, k, l) = rhs_vf(momxb + 1)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dx(j + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dx(j)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j + 1, k, l) = rhs_vf(momxb + 1)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dx(j + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dx(j)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j + 1, k, l) = rhs_vf(momxb)%sf(j + 1, k, l) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dx(j + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dx(j)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j + 1, k, l) = rhs_vf(momxb)%sf(j + 1, k, l) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dx(j + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dx(j)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j)), kind=stp)
                            end if

                            E_L = 0._wp; E_R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb + 1, vidxe
                                E_L = E_L + coeff_L(q + offxL)*q_cons_vf(E_idx)%sf(j + q, k, l)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb, vidxe - 1
                                E_R = E_R + coeff_R(q + offxR)*q_cons_vf(E_idx)%sf(j + q, k, l)
                            end do

                            call s_get_derived_states(E_L, gamma_L, pi_inf_L, rho_L, vel_L, &
                                                      E_R, gamma_R, pi_inf_R, rho_R, vel_R, &
                                                      pres_L, pres_R, cfl)

                            do i = 1, num_fluids
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j + 1, k, l) = rhs_vf(i)%sf(j + 1, k, l) + &
                                                            real((0.5_wp*dt*(alpha_rho_L(i)* &
                                                                             vel_L(1))*(1._wp/dx(j + 1)) - &
                                                                  0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                        real((0.5_wp*dt*(alpha_rho_L(i)* &
                                                                         vel_L(1))*(1._wp/dx(j)) - &
                                                              0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dx(j))), kind=stp)
                            end do

                            if (num_fluids > 1) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j + 1, k, l) = rhs_vf(advxb + i - 1)%sf(j + 1, k, l) + &
                                                                            real((0.5_wp*dt*(alpha_L(i)* &
                                                                                             vel_L(1))*(1._wp/dx(j + 1)) - &
                                                                                  0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dx(j + 1))), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j + 1, k, l) = rhs_vf(advxb + i - 1)%sf(j + 1, k, l) &
                                                                            - real((0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j + 1, k, l)*vel_L(1)*(1._wp/dx(j + 1))), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                        real((0.5_wp*dt*(alpha_L(i)* &
                                                                                         vel_L(1))*(1._wp/dx(j)) - &
                                                                              0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dx(j))), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                        + real((0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_L(1)*(1._wp/dx(j))), kind=stp)
                                end do
                            end if

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j + 1, k, l) = rhs_vf(momxb)%sf(j + 1, k, l) + &
                                                            real((0.5_wp*dt*(rho_L*(vel_L(1))**2.0 + &
                                                                             pres_L)*(1._wp/dx(j + 1)) - &
                                                                  0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dx(j + 1))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j + 1, k, l) = rhs_vf(momxb + 1)%sf(j + 1, k, l) + &
                                                                real((0.5_wp*dt*rho_L*vel_L(1)*vel_L(2)*(1._wp/dx(j + 1)) - &
                                                                      0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dx(j + 1))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) + &
                                                            real((0.5_wp*dt*(vel_L(1)*(E_L + &
                                                                                       pres_L))*(1._wp/dx(j + 1)) - &
                                                                  0.5_wp*dt*cfl*(E_L)*(1._wp/dx(j + 1))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                        real((0.5_wp*dt*(rho_L*(vel_L(1))**2.0 + &
                                                                         pres_L)*(1._wp/dx(j)) - &
                                                              0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dx(j))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                            real((0.5_wp*dt*rho_L*vel_L(1)*vel_L(2)*(1._wp/dx(j)) - &
                                                                  0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dx(j))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                        real((0.5_wp*dt*(vel_L(1)*(E_L + &
                                                                                   pres_L))*(1._wp/dx(j)) - &
                                                              0.5_wp*dt*cfl*(E_L)*(1._wp/dx(j))), kind=stp)

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j + 1, k, l) = rhs_vf(i)%sf(j + 1, k, l) + &
                                                            real((0.5_wp*dt*(alpha_rho_R(i)* &
                                                                             vel_R(1))*(1._wp/dx(j + 1)) + &
                                                                  0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                        real((0.5_wp*dt*(alpha_rho_R(i)* &
                                                                         vel_R(1))*(1._wp/dx(j)) + &
                                                              0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dx(j))), kind=stp)
                            end do

                            if (num_fluids > 1) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j + 1, k, l) = rhs_vf(advxb + i - 1)%sf(j + 1, k, l) + &
                                                                            real((0.5_wp*dt*(alpha_R(i)* &
                                                                                             vel_R(1))*(1._wp/dx(j + 1)) + &
                                                                                  0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dx(j + 1))), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j + 1, k, l) = rhs_vf(advxb + i - 1)%sf(j + 1, k, l) &
                                                                            - real((0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j + 1, k, l)*vel_R(1)*(1._wp/dx(j + 1))), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                        real((0.5_wp*dt*(alpha_R(i)* &
                                                                                         vel_R(1))*(1._wp/dx(j)) + &
                                                                              0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dx(j))), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                        + real((0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_R(1)*(1._wp/dx(j))), kind=stp)
                                end do
                            end if

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j + 1, k, l) = rhs_vf(momxb)%sf(j + 1, k, l) + &
                                                            real((0.5_wp*dt*(rho_R*(vel_R(1))**2.0 + &
                                                                             pres_R)*(1._wp/dx(j + 1)) + &
                                                                  0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dx(j + 1))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j + 1, k, l) = rhs_vf(momxb + 1)%sf(j + 1, k, l) + &
                                                                real((0.5_wp*dt*rho_R*vel_R(1)*vel_R(2)*(1._wp/dx(j + 1)) + &
                                                                      0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dx(j + 1))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) + &
                                                            real((0.5_wp*dt*(vel_R(1)*(E_R + &
                                                                                       pres_R))*(1._wp/dx(j + 1)) + &
                                                                  0.5_wp*dt*cfl*(E_R)*(1._wp/dx(j + 1))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                        real((0.5_wp*dt*(rho_R*(vel_R(1))**2.0 + &
                                                                         pres_R)*(1._wp/dx(j)) + &
                                                              0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dx(j))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                            real((0.5_wp*dt*rho_R*vel_R(1)*vel_R(2)*(1._wp/dx(j)) + &
                                                                  0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dx(j))), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                        real((0.5_wp*dt*(vel_R(1)*(E_R + &
                                                                                   pres_R))*(1._wp/dx(j)) + &
                                                              0.5_wp*dt*cfl*(E_R)*(1._wp/dx(j))), kind=stp)

                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else
                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j,k,l,rho_L, rho_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, mu_L, mu_R, vel_L, vel_R, pres_L, pres_R, alpha_L, alpha_R, alpha_rho_L, alpha_rho_R, F_L, F_R, E_L, E_R, cfl, dvel, dvel_small, rho_sf_small, vflux_L_arr, vflux_R_arr]')
                    do l = 0, p
                        do k = 0, n
                            do j = -1, m

                                vflux_L_arr = 0._wp
                                vflux_R_arr = 0._wp

                                #:if MFC_CASE_OPTIMIZATION
                                    #:if igr_order == 5
                                        !DIR$ unroll 6
                                    #:elif igr_order == 3
                                        !DIR$ unroll 4
                                    #:endif
                                #:endif
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb, vidxe

                                    !x-direction contributions
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_cons_vf(r)%sf(j + i + q, k, l)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dx(j)))*( &
                                                    q_cons_vf(momxb)%sf(j + 1 + q, k, l)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb)%sf(j - 1 + q, k, l)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dx(j)))*( &
                                                    q_cons_vf(momxb + 1)%sf(j + 1 + q, k, l)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 1)%sf(j - 1 + q, k, l)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dx(j)))*( &
                                                    q_cons_vf(momxb + 2)%sf(j + 1 + q, k, l)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 2)%sf(j - 1 + q, k, l)/rho_sf_small(-1))

                                    if (q == 0) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, num_dims
                                            dvel(i, 1) = dvel_small(i)
                                        end do
                                    end if

                                    if (q > vidxb) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(2))
                                        vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q + offxL)*(dvel_small(3))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(4._wp*dvel_small(1))/3._wp
                                    end if
                                    if (q < vidxe) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(2))
                                        vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q + offxR)*(dvel_small(3))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(4._wp*dvel_small(1))/3._wp
                                    end if

                                    !y-direction contributions
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_cons_vf(r)%sf(j + q, k + i, l)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dy(k)))*( &
                                                    q_cons_vf(momxb)%sf(j + q, k + 1, l)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb)%sf(j + q, k - 1, l)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dy(k)))*( &
                                                    q_cons_vf(momxb + 1)%sf(j + q, k + 1, l)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 1)%sf(j + q, k - 1, l)/rho_sf_small(-1))
                                    if (q == 0) dvel_small(3) = (1/(2._wp*dy(k)))*( &
                                                                q_cons_vf(momxb + 2)%sf(j + q, k + 1, l)/rho_sf_small(1) - &
                                                                q_cons_vf(momxb + 2)%sf(j + q, k - 1, l)/rho_sf_small(-1))
                                    if (q == 0) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, num_dims
                                            dvel(i, 2) = dvel_small(i)
                                        end do
                                    end if

                                    if (q > vidxb) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(1))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(-2._wp*dvel_small(2))/3._wp
                                    end if
                                    if (q < vidxe) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(1))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(-2._wp*dvel_small(2))/3._wp
                                    end if

                                    !z-direction contributions
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_cons_vf(r)%sf(j + q, k, l + i)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dz(l)))*( &
                                                    q_cons_vf(momxb)%sf(j + q, k, l + 1)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb)%sf(j + q, k, l - 1)/rho_sf_small(-1))
                                    if (q == 0) dvel_small(2) = (1/(2._wp*dz(l)))*( &
                                                                q_cons_vf(momxb + 1)%sf(j + q, k, l + 1)/rho_sf_small(1) - &
                                                                q_cons_vf(momxb + 1)%sf(j + q, k, l - 1)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dz(l)))*( &
                                                    q_cons_vf(momxb + 2)%sf(j + q, k, l + 1)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 2)%sf(j + q, k, l - 1)/rho_sf_small(-1))
                                    if (q == 0) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, num_dims
                                            dvel(i, 3) = dvel_small(i)
                                        end do
                                    end if

                                    if (q > vidxb) then
                                        vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q + offxL)*(dvel_small(1))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(-2._wp*dvel_small(3))/3._wp
                                    end if
                                    if (q < vidxe) then
                                        vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q + offxR)*(dvel_small(1))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(-2._wp*dvel_small(3))/3._wp
                                    end if

                                    if (q == 0) then
                                        jac_rhs(j, k, l) = real(alf_igr*(2._wp*(dvel(1, 2)*dvel(2, 1) &
                                                                                + dvel(1, 3)*dvel(3, 1) &
                                                                                + dvel(2, 3)*dvel(3, 2)) &
                                                                         + dvel(1, 1)**2._wp + dvel(2, 2)**2._wp &
                                                                         + dvel(3, 3)**2._wp &
                                                                         + (dvel(1, 1) + dvel(2, 2) + dvel(3, 3))**2._wp), kind=stp)
                                    end if
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = 0._wp
                                    alpha_rho_R(i) = 0._wp
                                    alpha_L(i) = 0._wp
                                    alpha_R(i) = 0._wp
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = 0._wp
                                    vel_R(i) = 0._wp
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb + 1, vidxe
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        alpha_rho_L(i) = alpha_rho_L(i) + coeff_L(q + offxL)*q_cons_vf(i)%sf(j + q, k, l)
                                    end do

                                    if (num_fluids > 1) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, num_fluids - 1
                                            alpha_L(i) = alpha_L(i) + coeff_L(q + offxL)*q_cons_vf(E_idx + i)%sf(j + q, k, l)
                                        end do
                                    else
                                        alpha_L(1) = 1._wp
                                    end if

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        vel_L(i) = vel_L(i) + coeff_L(q + offxL)*q_cons_vf(momxb + i - 1)%sf(j + q, k, l)
                                    end do
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb, vidxe - 1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        alpha_rho_R(i) = alpha_rho_R(i) + coeff_R(q + offxR)*q_cons_vf(i)%sf(j + q, k, l)
                                    end do

                                    if (num_fluids > 1) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, num_fluids - 1
                                            alpha_R(i) = alpha_R(i) + coeff_R(q + offxR)*q_cons_vf(E_idx + i)%sf(j + q, k, l)
                                        end do
                                    else
                                        alpha_R(1) = 1._wp
                                    end if

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        vel_R(i) = vel_R(i) + coeff_R(q + offxR)*q_cons_vf(momxb + i - 1)%sf(j + q, k, l)
                                    end do
                                end do

                                if (num_fluids > 1) then

                                    alpha_L(num_fluids) = 1._wp
                                    alpha_R(num_fluids) = 1._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        alpha_L(num_fluids) = alpha_L(num_fluids) - alpha_L(i)
                                        alpha_R(num_fluids) = alpha_R(num_fluids) - alpha_R(i)
                                    end do
                                end if

                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_L = rho_L + alpha_rho_L(i)
                                    gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                    pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)

                                    rho_R = rho_R + alpha_rho_R(i)
                                    gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                    pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = vel_L(i)/rho_L
                                    vel_R(i) = vel_R(i)/rho_R
                                end do

                                if (viscous) then
                                    mu_L = 0._wp
                                    mu_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        mu_L = alpha_L(i)/Res_igr(1, i) + mu_L
                                        mu_R = alpha_R(i)/Res_igr(1, i) + mu_R
                                    end do

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 1)%sf(j + 1, k, l) = rhs_vf(momxb + 1)%sf(j + 1, k, l) - &
                                                                        real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dx(j + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dx(j)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 1)%sf(j + 1, k, l) = rhs_vf(momxb + 1)%sf(j + 1, k, l) - &
                                                                        real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dx(j + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dx(j)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 2)%sf(j + 1, k, l) = rhs_vf(momxb + 2)%sf(j + 1, k, l) - &
                                                                        real(0.5_wp*dt*mu_L*vflux_L_arr(2)*(1._wp/dx(j + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dx(j + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) + &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(2)*(1._wp/dx(j)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dx(j)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 2)%sf(j + 1, k, l) = rhs_vf(momxb + 2)%sf(j + 1, k, l) - &
                                                                        real(0.5_wp*dt*mu_R*vflux_R_arr(2)*(1._wp/dx(j + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dx(j + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) + &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(2)*(1._wp/dx(j)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dx(j)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb)%sf(j + 1, k, l) = rhs_vf(momxb)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dx(j + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dx(j)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb)%sf(j + 1, k, l) = rhs_vf(momxb)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dx(j + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dx(j)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j)), kind=stp)
                                end if

                                E_L = 0._wp; E_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb + 1, vidxe
                                    E_L = E_L + coeff_L(q + offxL)*q_cons_vf(E_idx)%sf(j + q, k, l)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb, vidxe - 1
                                    E_R = E_R + coeff_R(q + offxR)*q_cons_vf(E_idx)%sf(j + q, k, l)
                                end do

                                call s_get_derived_states(E_L, gamma_L, pi_inf_L, rho_L, vel_L, &
                                                          E_R, gamma_R, pi_inf_R, rho_R, vel_R, &
                                                          pres_L, pres_R, cfl)

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(i)%sf(j + 1, k, l) = rhs_vf(i)%sf(j + 1, k, l) + &
                                                                real((0.5_wp*dt*(alpha_rho_L(i)* &
                                                                                 vel_L(1))*(1._wp/dx(j + 1)) - &
                                                                      0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dx(j + 1))), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                            real((0.5_wp*dt*(alpha_rho_L(i)* &
                                                                             vel_L(1))*(1._wp/dx(j)) - &
                                                                  0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dx(j))), kind=stp)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j + 1, k, l) = rhs_vf(advxb + i - 1)%sf(j + 1, k, l) + &
                                                                                real((0.5_wp*dt*(alpha_L(i)* &
                                                                                                 vel_L(1))*(1._wp/dx(j + 1)) - &
                                                                                      0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dx(j + 1))), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j + 1, k, l) = rhs_vf(advxb + i - 1)%sf(j + 1, k, l) &
                                                                                - real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j + 1, k, l)*vel_L(1)*(1._wp/dx(j + 1)), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                            real(0.5_wp*dt*(alpha_L(i)* &
                                                                                            vel_L(1))*(1._wp/dx(j)) - &
                                                                                 0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dx(j)), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                            + real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_L(1)*(1._wp/dx(j)), kind=stp)
                                    end do
                                end if

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j + 1, k, l) = rhs_vf(momxb)%sf(j + 1, k, l) + &
                                                                real((0.5_wp*dt*(rho_L*(vel_L(1))**2.0 + &
                                                                                 pres_L)*(1._wp/dx(j + 1)) - &
                                                                      0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j + 1, k, l) = rhs_vf(momxb + 1)%sf(j + 1, k, l) + &
                                                                    real((0.5_wp*dt*rho_L*vel_L(1)*vel_L(2)*(1._wp/dx(j + 1)) - &
                                                                          0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j + 1, k, l) = rhs_vf(momxb + 2)%sf(j + 1, k, l) + &
                                                                    real((0.5_wp*dt*rho_L*vel_L(1)*vel_L(3)*(1._wp/dx(j + 1)) - &
                                                                          0.5_wp*dt*cfl*(rho_L*vel_L(3))*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) + &
                                                                real((0.5_wp*dt*(vel_L(1)*(E_L + &
                                                                                           pres_L))*(1._wp/dx(j + 1)) - &
                                                                      0.5_wp*dt*cfl*(E_L)*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                            real((0.5_wp*dt*(rho_L*(vel_L(1))**2.0 + &
                                                                             pres_L)*(1._wp/dx(j)) - &
                                                                  0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dx(j))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                                real((0.5_wp*dt*rho_L*vel_L(1)*vel_L(2)*(1._wp/dx(j)) - &
                                                                      0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dx(j))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) - &
                                                                real((0.5_wp*dt*rho_L*vel_L(1)*vel_L(3)*(1._wp/dx(j)) - &
                                                                      0.5_wp*dt*cfl*(rho_L*vel_L(3))*(1._wp/dx(j))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                            real((0.5_wp*dt*(vel_L(1)*(E_L + &
                                                                                       pres_L))*(1._wp/dx(j)) - &
                                                                  0.5_wp*dt*cfl*(E_L)*(1._wp/dx(j))), kind=stp)

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(i)%sf(j + 1, k, l) = rhs_vf(i)%sf(j + 1, k, l) + &
                                                                real((0.5_wp*dt*(alpha_rho_R(i)* &
                                                                                 vel_R(1))*(1._wp/dx(j + 1)) + &
                                                                      0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dx(j + 1))), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                            real((0.5_wp*dt*(alpha_rho_R(i)* &
                                                                             vel_R(1))*(1._wp/dx(j)) + &
                                                                  0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dx(j))), kind=stp)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j + 1, k, l) = rhs_vf(advxb + i - 1)%sf(j + 1, k, l) + &
                                                                                real((0.5_wp*dt*(alpha_R(i)* &
                                                                                                 vel_R(1))*(1._wp/dx(j + 1)) + &
                                                                                      0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dx(j + 1))), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j + 1, k, l) = rhs_vf(advxb + i - 1)%sf(j + 1, k, l) &
                                                                                - real((0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j + 1, k, l)*vel_R(1)*(1._wp/dx(j + 1))), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                            real((0.5_wp*dt*(alpha_R(i)* &
                                                                                             vel_R(1))*(1._wp/dx(j)) + &
                                                                                  0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dx(j))), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                            + real((0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_R(1)*(1._wp/dx(j))), kind=stp)
                                    end do
                                end if

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j + 1, k, l) = rhs_vf(momxb)%sf(j + 1, k, l) + &
                                                                real((0.5_wp*dt*(rho_R*(vel_R(1))**2.0 + &
                                                                                 pres_R)*(1._wp/dx(j + 1)) + &
                                                                      0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j + 1, k, l) = rhs_vf(momxb + 1)%sf(j + 1, k, l) + &
                                                                    real((0.5_wp*dt*rho_R*vel_R(1)*vel_R(2)*(1._wp/dx(j + 1)) + &
                                                                          0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j + 1, k, l) = rhs_vf(momxb + 2)%sf(j + 1, k, l) + &
                                                                    real((0.5_wp*dt*rho_R*vel_R(1)*vel_R(3)*(1._wp/dx(j + 1)) + &
                                                                          0.5_wp*dt*cfl*(rho_R*vel_R(3))*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j + 1, k, l) = rhs_vf(E_idx)%sf(j + 1, k, l) + &
                                                                real((0.5_wp*dt*(vel_R(1)*(E_R + &
                                                                                           pres_R))*(1._wp/dx(j + 1)) + &
                                                                      0.5_wp*dt*cfl*(E_R)*(1._wp/dx(j + 1))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                            real((0.5_wp*dt*(rho_R*(vel_R(1))**2.0 + &
                                                                             pres_R)*(1._wp/dx(j)) + &
                                                                  0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dx(j))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                                real((0.5_wp*dt*rho_R*vel_R(1)*vel_R(2)*(1._wp/dx(j)) + &
                                                                      0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dx(j))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) - &
                                                                real((0.5_wp*dt*rho_R*vel_R(1)*vel_R(3)*(1._wp/dx(j)) + &
                                                                      0.5_wp*dt*cfl*(rho_R*vel_R(3))*(1._wp/dx(j))), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                            real((0.5_wp*dt*(vel_R(1)*(E_R + &
                                                                                       pres_R))*(1._wp/dx(j)) + &
                                                                  0.5_wp*dt*cfl*(E_R)*(1._wp/dx(j))), kind=stp)

                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:endif
            end if
        else if (idir == 2) then
            if (p == 0) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j,k,l,rho_L, rho_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, mu_L, mu_R, vel_L, vel_R, pres_L, pres_R, alpha_L, alpha_R, alpha_rho_L, alpha_rho_R, F_L, F_R, E_L, E_R, cfl, dvel_small, rho_sf_small, vflux_L_arr, vflux_R_arr]')
                do l = 0, p
                    do k = -1, n
                        do j = 0, m

                            if (viscous) then
                                vflux_L_arr = 0._wp
                                vflux_R_arr = 0._wp

                                #:if MFC_CASE_OPTIMIZATION
                                    #:if igr_order == 5
                                        !DIR$ unroll 6
                                    #:elif igr_order == 3
                                        !DIR$ unroll 4
                                    #:endif
                                #:endif
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb, vidxe

                                    !x-direction contributions
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_cons_vf(r)%sf(j + i, k + q, l)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dx(j)))*( &
                                                    q_cons_vf(momxb)%sf(j + 1, k + q, l)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb)%sf(j - 1, k + q, l)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dx(j)))*( &
                                                    q_cons_vf(momxb + 1)%sf(j + 1, k + q, l)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 1)%sf(j - 1, k + q, l)/rho_sf_small(-1))

                                    if (q > vidxb) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(2))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(-2._wp*dvel_small(1))/3._wp
                                    end if
                                    if (q < vidxe) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(2))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(-2._wp*dvel_small(1))/3._wp
                                    end if

                                    !y-direction contributions
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_cons_vf(r)%sf(j, k + i + q, l)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dy(k)))*( &
                                                    q_cons_vf(momxb)%sf(j, k + 1 + q, l)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb)%sf(j, k - 1 + q, l)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dy(k)))*( &
                                                    q_cons_vf(momxb + 1)%sf(j, k + 1 + q, l)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 1)%sf(j, k - 1 + q, l)/rho_sf_small(-1))

                                    if (q > vidxb) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(1))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(4._wp*dvel_small(2))/3._wp
                                    end if
                                    if (q < vidxe) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(1))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(4._wp*dvel_small(2))/3._wp
                                    end if
                                end do
                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_L(i) = 0._wp
                                alpha_rho_R(i) = 0._wp
                                alpha_L(i) = 0._wp
                                alpha_R(i) = 0._wp
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = 0._wp
                                vel_R(i) = 0._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb + 1, vidxe
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = alpha_rho_L(i) + coeff_L(q + offxL)*q_cons_vf(i)%sf(j, k + q, l)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        alpha_L(i) = alpha_L(i) + coeff_L(q + offxL)*q_cons_vf(E_idx + i)%sf(j, k + q, l)
                                    end do
                                else
                                    alpha_L(1) = 1._wp
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = vel_L(i) + coeff_L(q + offxL)*q_cons_vf(momxb + i - 1)%sf(j, k + q, l)
                                end do
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb, vidxe - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_R(i) = alpha_rho_R(i) + coeff_R(q + offxR)*q_cons_vf(i)%sf(j, k + q, l)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        alpha_R(i) = alpha_R(i) + coeff_R(q + offxR)*q_cons_vf(E_idx + i)%sf(j, k + q, l)
                                    end do
                                else
                                    alpha_R(1) = 1._wp
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_R(i) = vel_R(i) + coeff_R(q + offxR)*q_cons_vf(momxb + i - 1)%sf(j, k + q, l)
                                end do
                            end do

                            if (num_fluids > 1) then

                                alpha_L(num_fluids) = 1._wp
                                alpha_R(num_fluids) = 1._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    alpha_L(num_fluids) = alpha_L(num_fluids) - alpha_L(i)
                                    alpha_R(num_fluids) = alpha_R(num_fluids) - alpha_R(i)
                                end do
                            end if

                            rho_L = 0._wp; rho_R = 0._wp
                            gamma_L = 0._wp; gamma_R = 0._wp
                            pi_inf_L = 0._wp; pi_inf_R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho_L = rho_L + alpha_rho_L(i)
                                gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)

                                rho_R = rho_R + alpha_rho_R(i)
                                gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = vel_L(i)/rho_L
                                vel_R(i) = vel_R(i)/rho_R
                            end do

                            if (viscous) then
                                mu_L = 0._wp
                                mu_R = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    mu_L = alpha_L(i)/Res_igr(1, i) + mu_L
                                    mu_R = alpha_R(i)/Res_igr(1, i) + mu_R
                                end do

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k + 1, l) = rhs_vf(momxb)%sf(j, k + 1, l) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dy(k + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dy(k)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k + 1, l) = rhs_vf(momxb)%sf(j, k + 1, l) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dy(k + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dy(k)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k + 1, l) = rhs_vf(momxb + 1)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dy(k + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dy(k)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k + 1, l) = rhs_vf(momxb + 1)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dy(k + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dy(k)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k)), kind=stp)
                            end if

                            E_L = 0._wp; E_R = 0._wp
                            F_L = 0._wp; F_R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb + 1, vidxe
                                E_L = E_L + coeff_L(q + offxL)*q_cons_vf(E_idx)%sf(j, k + q, l)
                                F_L = F_L + coeff_L(q + offxL)*jac(j, k + q, l)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb, vidxe - 1
                                E_R = E_R + coeff_R(q + offxR)*q_cons_vf(E_idx)%sf(j, k + q, l)
                                F_R = F_R + coeff_R(q + offxR)*jac(j, k + q, l)
                            end do

                            call s_get_derived_states(E_L, gamma_L, pi_inf_L, rho_L, vel_L, &
                                                      E_R, gamma_R, pi_inf_R, rho_R, vel_R, &
                                                      pres_L, pres_R, cfl)

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k + 1, l) = rhs_vf(i)%sf(j, k + 1, l) + &
                                                            real(0.5_wp*dt*(alpha_rho_L(i)* &
                                                                            vel_L(2))*(1._wp/dy(k + 1)) - &
                                                                 0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*(alpha_rho_L(i)* &
                                                                        vel_L(2))*(1._wp/dy(k)) - &
                                                             0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dy(k)), kind=stp)
                            end do

                            if (num_fluids > 1) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k + 1, l) = rhs_vf(advxb + i - 1)%sf(j, k + 1, l) + &
                                                                            real(0.5_wp*dt*(alpha_L(i)* &
                                                                                            vel_L(2))*(1._wp/dy(k + 1)) - &
                                                                                 0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k + 1, l) = rhs_vf(advxb + i - 1)%sf(j, k + 1, l) &
                                                                            - real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k + 1, l)*vel_L(2)*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                        real(0.5_wp*dt*(alpha_L(i)* &
                                                                                        vel_L(2))*(1._wp/dy(k)) - &
                                                                             0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dy(k)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                        + real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_L(2)*(1._wp/dy(k)), kind=stp)
                                end do
                            end if

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k + 1, l) = rhs_vf(momxb + 1)%sf(j, k + 1, l) + &
                                                                real(0.5_wp*dt*(rho_L*(vel_L(2))**2.0 + &
                                                                                pres_L + F_L)*(1._wp/dy(k + 1)) - &
                                                                     0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dy(k + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k + 1, l) = rhs_vf(momxb)%sf(j, k + 1, l) + &
                                                            real(0.5_wp*dt*rho_L*vel_L(1)*vel_L(2)*(1._wp/dy(k + 1)) - &
                                                                 0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dy(k + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) + &
                                                            real(0.5_wp*dt*(vel_L(2)*(E_L + &
                                                                                      pres_L + F_L))*(1._wp/dy(k + 1)) - &
                                                                 0.5_wp*dt*cfl*(E_L)*(1._wp/dy(k + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*(rho_L*(vel_L(2))**2.0 + &
                                                                            pres_L + F_L)*(1._wp/dy(k)) - &
                                                                 0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dy(k)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*rho_L*vel_L(1)*vel_L(2)*(1._wp/dy(k)) - &
                                                             0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dy(k)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*(vel_L(2)*(E_L + &
                                                                                  pres_L + F_L))*(1._wp/dy(k)) - &
                                                             0.5_wp*dt*cfl*(E_L)*(1._wp/dy(k)), kind=stp)

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k + 1, l) = rhs_vf(i)%sf(j, k + 1, l) + &
                                                            real(0.5_wp*dt*(alpha_rho_R(i)* &
                                                                            vel_R(2))*(1._wp/dy(k + 1)) + &
                                                                 0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dy(k + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*(alpha_rho_R(i)* &
                                                                        vel_R(2))*(1._wp/dy(k)) + &
                                                             0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dy(k)), kind=stp)
                            end do

                            if (num_fluids > 1) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k + 1, l) = rhs_vf(advxb + i - 1)%sf(j, k + 1, l) + &
                                                                            real(0.5_wp*dt*(alpha_R(i)* &
                                                                                            vel_R(2))*(1._wp/dy(k + 1)) + &
                                                                                 0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k + 1, l) = rhs_vf(advxb + i - 1)%sf(j, k + 1, l) &
                                                                            - real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k + 1, l)*vel_R(2)*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                        real(0.5_wp*dt*(alpha_R(i)* &
                                                                                        vel_R(2))*(1._wp/dy(k)) + &
                                                                             0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dy(k)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                        + real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_R(2)*(1._wp/dy(k)), kind=stp)
                                end do
                            end if
                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k + 1, l) = rhs_vf(momxb + 1)%sf(j, k + 1, l) + &
                                                                real(0.5_wp*dt*(rho_R*(vel_R(2))**2.0 + &
                                                                                pres_R + F_R)*(1._wp/dy(k + 1)) + &
                                                                     0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dy(k + 1)), kind=stp)
                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k + 1, l) = rhs_vf(momxb)%sf(j, k + 1, l) + &
                                                            real(0.5_wp*dt*rho_R*vel_R(2)*vel_R(1)*(1._wp/dy(k + 1)) + &
                                                                 0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dy(k + 1)), kind=stp)
                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) + &
                                                            real(0.5_wp*dt*(vel_R(2)*(E_R + &
                                                                                      pres_R + F_R))*(1._wp/dy(k + 1)) + &
                                                                 0.5_wp*dt*cfl*(E_R)*(1._wp/dy(k + 1)), kind=stp)
                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*(rho_R*(vel_R(2))**2.0 + &
                                                                            pres_R + F_R)*(1._wp/dy(k)) + &
                                                                 0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dy(k)), kind=stp)
                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*rho_R*vel_R(2)*vel_R(1)*(1._wp/dy(k)) + &
                                                             0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dy(k)), kind=stp)
                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*(vel_R(2)*(E_R + &
                                                                                  pres_R + F_R))*(1._wp/dy(k)) + &
                                                             0.5_wp*dt*cfl*(E_R)*(1._wp/dy(k)), kind=stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else
                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j,k,l,rho_L, rho_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, mu_L, mu_R, vel_L, vel_R, pres_L, pres_R, alpha_L, alpha_R, alpha_rho_L, alpha_rho_R, F_L, F_R, E_L, E_R, cfl, dvel_small, rho_sf_small, vflux_L_arr, vflux_R_arr]')
                    do l = 0, p
                        do k = -1, n
                            do j = 0, m

                                if (viscous) then
                                    vflux_L_arr = 0._wp
                                    vflux_R_arr = 0._wp

                                    #:if MFC_CASE_OPTIMIZATION
                                        #:if igr_order == 5
                                            !DIR$ unroll 6
                                        #:elif igr_order == 3
                                            !DIR$ unroll 4
                                        #:endif
                                    #:endif
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = vidxb, vidxe

                                        !x-direction contributions
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = -1, 1
                                            rho_L = 0._wp
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do r = 1, num_fluids
                                                rho_L = rho_L + q_cons_vf(r)%sf(j + i, k + q, l)
                                            end do
                                            rho_sf_small(i) = rho_L
                                        end do

                                        dvel_small(1) = (1/(2._wp*dx(j)))*( &
                                                        q_cons_vf(momxb)%sf(j + 1, k + q, l)/rho_sf_small(1) - &
                                                        q_cons_vf(momxb)%sf(j - 1, k + q, l)/rho_sf_small(-1))
                                        dvel_small(2) = (1/(2._wp*dx(j)))*( &
                                                        q_cons_vf(momxb + 1)%sf(j + 1, k + q, l)/rho_sf_small(1) - &
                                                        q_cons_vf(momxb + 1)%sf(j - 1, k + q, l)/rho_sf_small(-1))

                                        if (q > vidxb) then
                                            vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(2))
                                            vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(-2._wp*dvel_small(1))/3._wp
                                        end if
                                        if (q < vidxe) then
                                            vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(2))
                                            vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(-2._wp*dvel_small(1))/3._wp
                                        end if

                                        !y-direction contributions
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = -1, 1
                                            rho_L = 0._wp
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do r = 1, num_fluids
                                                rho_L = rho_L + q_cons_vf(r)%sf(j, k + i + q, l)
                                            end do
                                            rho_sf_small(i) = rho_L
                                        end do

                                        dvel_small(1) = (1/(2._wp*dy(k)))*( &
                                                        q_cons_vf(momxb)%sf(j, k + 1 + q, l)/rho_sf_small(1) - &
                                                        q_cons_vf(momxb)%sf(j, k - 1 + q, l)/rho_sf_small(-1))
                                        dvel_small(2) = (1/(2._wp*dy(k)))*( &
                                                        q_cons_vf(momxb + 1)%sf(j, k + 1 + q, l)/rho_sf_small(1) - &
                                                        q_cons_vf(momxb + 1)%sf(j, k - 1 + q, l)/rho_sf_small(-1))
                                        dvel_small(3) = (1/(2._wp*dy(k)))*( &
                                                        q_cons_vf(momxb + 2)%sf(j, k + 1 + q, l)/rho_sf_small(1) - &
                                                        q_cons_vf(momxb + 2)%sf(j, k - 1 + q, l)/rho_sf_small(-1))

                                        if (q > vidxb) then
                                            vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(1))
                                            vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q + offxL)*(dvel_small(3))
                                            vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(4._wp*dvel_small(2))/3._wp
                                        end if
                                        if (q < vidxe) then
                                            vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(1))
                                            vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q + offxR)*(dvel_small(3))
                                            vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(4._wp*dvel_small(2))/3._wp
                                        end if

                                        !z-direction contributions
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = -1, 1
                                            rho_L = 0._wp
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do r = 1, num_fluids
                                                rho_L = rho_L + q_cons_vf(r)%sf(j, k + q, l + i)
                                            end do
                                            rho_sf_small(i) = rho_L
                                        end do

                                        dvel_small(2) = (1/(2._wp*dz(l)))*( &
                                                        q_cons_vf(momxb + 1)%sf(j, k + q, l + 1)/rho_sf_small(1) - &
                                                        q_cons_vf(momxb + 1)%sf(j, k + q, l - 1)/rho_sf_small(-1))
                                        dvel_small(3) = (1/(2._wp*dz(l)))*( &
                                                        q_cons_vf(momxb + 2)%sf(j, k + q, l + 1)/rho_sf_small(1) - &
                                                        q_cons_vf(momxb + 2)%sf(j, k + q, l - 1)/rho_sf_small(-1))
                                        if (q > vidxb) then
                                            vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q + offxL)*(dvel_small(2))
                                            vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(-2._wp*dvel_small(3))/3._wp
                                        end if
                                        if (q < vidxe) then
                                            vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q + offxR)*(dvel_small(2))
                                            vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(-2._wp*dvel_small(3))/3._wp
                                        end if
                                    end do
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = 0._wp
                                    alpha_rho_R(i) = 0._wp
                                    alpha_L(i) = 0._wp
                                    alpha_R(i) = 0._wp
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = 0._wp
                                    vel_R(i) = 0._wp
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb + 1, vidxe
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        alpha_rho_L(i) = alpha_rho_L(i) + coeff_L(q + offxL)*q_cons_vf(i)%sf(j, k + q, l)
                                    end do

                                    if (num_fluids > 1) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, num_fluids - 1
                                            alpha_L(i) = alpha_L(i) + coeff_L(q + offxL)*q_cons_vf(E_idx + i)%sf(j, k + q, l)
                                        end do
                                    else
                                        alpha_L(1) = 1._wp
                                    end if

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        vel_L(i) = vel_L(i) + coeff_L(q + offxL)*q_cons_vf(momxb + i - 1)%sf(j, k + q, l)
                                    end do
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb, vidxe - 1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        alpha_rho_R(i) = alpha_rho_R(i) + coeff_R(q + offxR)*q_cons_vf(i)%sf(j, k + q, l)
                                    end do

                                    if (num_fluids > 1) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, num_fluids - 1
                                            alpha_R(i) = alpha_R(i) + coeff_R(q + offxR)*q_cons_vf(E_idx + i)%sf(j, k + q, l)
                                        end do
                                    else
                                        alpha_R(1) = 1._wp
                                    end if

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        vel_R(i) = vel_R(i) + coeff_R(q + offxR)*q_cons_vf(momxb + i - 1)%sf(j, k + q, l)
                                    end do
                                end do

                                if (num_fluids > 1) then

                                    alpha_L(num_fluids) = 1._wp
                                    alpha_R(num_fluids) = 1._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        alpha_L(num_fluids) = alpha_L(num_fluids) - alpha_L(i)
                                        alpha_R(num_fluids) = alpha_R(num_fluids) - alpha_R(i)
                                    end do
                                end if

                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_L = rho_L + alpha_rho_L(i)
                                    gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                    pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)

                                    rho_R = rho_R + alpha_rho_R(i)
                                    gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                    pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = vel_L(i)/rho_L
                                    vel_R(i) = vel_R(i)/rho_R
                                end do

                                if (viscous) then
                                    mu_L = 0._wp
                                    mu_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        mu_L = alpha_L(i)/Res_igr(1, i) + mu_L
                                        mu_R = alpha_R(i)/Res_igr(1, i) + mu_R
                                    end do

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb)%sf(j, k + 1, l) = rhs_vf(momxb)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dy(k + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dy(k)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb)%sf(j, k + 1, l) = rhs_vf(momxb)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dy(k + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dy(k)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 2)%sf(j, k + 1, l) = rhs_vf(momxb + 2)%sf(j, k + 1, l) - &
                                                                        real(0.5_wp*dt*mu_L*vflux_L_arr(2)*(1._wp/dy(k + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) + &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(2)*(1._wp/dy(k)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dy(k)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 2)%sf(j, k + 1, l) = rhs_vf(momxb + 2)%sf(j, k + 1, l) - &
                                                                        real(0.5_wp*dt*mu_R*vflux_R_arr(2)*(1._wp/dy(k + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) + &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(2)*(1._wp/dy(k)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dy(k)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 1)%sf(j, k + 1, l) = rhs_vf(momxb + 1)%sf(j, k + 1, l) - &
                                                                        real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dy(k + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dy(k)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 1)%sf(j, k + 1, l) = rhs_vf(momxb + 1)%sf(j, k + 1, l) - &
                                                                        real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dy(k + 1)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dy(k)), kind=stp)
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k)), kind=stp)
                                end if

                                E_L = 0._wp; E_R = 0._wp
                                F_L = 0._wp; F_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb + 1, vidxe
                                    E_L = E_L + coeff_L(q + offxL)*q_cons_vf(E_idx)%sf(j, k + q, l)
                                    F_L = F_L + coeff_L(q + offxL)*jac(j, k + q, l)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb, vidxe - 1
                                    E_R = E_R + coeff_R(q + offxR)*q_cons_vf(E_idx)%sf(j, k + q, l)
                                    F_R = F_R + coeff_R(q + offxR)*jac(j, k + q, l)
                                end do

                                call s_get_derived_states(E_L, gamma_L, pi_inf_L, rho_L, vel_L, &
                                                          E_R, gamma_R, pi_inf_R, rho_R, vel_R, &
                                                          pres_L, pres_R, cfl)

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(i)%sf(j, k + 1, l) = rhs_vf(i)%sf(j, k + 1, l) + &
                                                                real(0.5_wp*dt*(alpha_rho_L(i)* &
                                                                                vel_L(2))*(1._wp/dy(k + 1)) - &
                                                                     0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*(alpha_rho_L(i)* &
                                                                            vel_L(2))*(1._wp/dy(k)) - &
                                                                 0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dy(k)), kind=stp)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k + 1, l) = rhs_vf(advxb + i - 1)%sf(j, k + 1, l) + &
                                                                                real(0.5_wp*dt*(alpha_L(i)* &
                                                                                                vel_L(2))*(1._wp/dy(k + 1)) - &
                                                                                     0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dy(k + 1)), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k + 1, l) = rhs_vf(advxb + i - 1)%sf(j, k + 1, l) &
                                                                                - real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k + 1, l)*vel_L(2)*(1._wp/dy(k + 1)), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                            real(0.5_wp*dt*(alpha_L(i)* &
                                                                                            vel_L(2))*(1._wp/dy(k)) - &
                                                                                 0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dy(k)), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                            + real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_L(2)*(1._wp/dy(k)), kind=stp)
                                    end do
                                end if

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k + 1, l) = rhs_vf(momxb + 1)%sf(j, k + 1, l) + &
                                                                    real(0.5_wp*dt*(rho_L*(vel_L(2))**2.0 + &
                                                                                    pres_L + F_L)*(1._wp/dy(k + 1)) - &
                                                                         0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k + 1, l) = rhs_vf(momxb)%sf(j, k + 1, l) + &
                                                                real(0.5_wp*dt*rho_L*vel_L(1)*vel_L(2)*(1._wp/dy(k + 1)) - &
                                                                     0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k + 1, l) = rhs_vf(momxb + 2)%sf(j, k + 1, l) + &
                                                                    real(0.5_wp*dt*rho_L*vel_L(3)*vel_L(2)*(1._wp/dy(k + 1)) - &
                                                                         0.5_wp*dt*cfl*(rho_L*vel_L(3))*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) + &
                                                                real(0.5_wp*dt*(vel_L(2)*(E_L + &
                                                                                          pres_L + F_L))*(1._wp/dy(k + 1)) - &
                                                                     0.5_wp*dt*cfl*(E_L)*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                                real(0.5_wp*dt*(rho_L*(vel_L(2))**2.0 + &
                                                                                pres_L + F_L)*(1._wp/dy(k)) - &
                                                                     0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dy(k)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*rho_L*vel_L(1)*vel_L(2)*(1._wp/dy(k)) - &
                                                                 0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dy(k)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) - &
                                                                real(0.5_wp*dt*rho_L*vel_L(3)*vel_L(2)*(1._wp/dy(k)) - &
                                                                     0.5_wp*dt*cfl*(rho_L*vel_L(3))*(1._wp/dy(k)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*(vel_L(2)*(E_L + &
                                                                                      pres_L + F_L))*(1._wp/dy(k)) - &
                                                                 0.5_wp*dt*cfl*(E_L)*(1._wp/dy(k)), kind=stp)

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(i)%sf(j, k + 1, l) = rhs_vf(i)%sf(j, k + 1, l) + &
                                                                real(0.5_wp*dt*(alpha_rho_R(i)* &
                                                                                vel_R(2))*(1._wp/dy(k + 1)) + &
                                                                     0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dy(k + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*(alpha_rho_R(i)* &
                                                                            vel_R(2))*(1._wp/dy(k)) + &
                                                                 0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dy(k)), kind=stp)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k + 1, l) = rhs_vf(advxb + i - 1)%sf(j, k + 1, l) + &
                                                                                real(0.5_wp*dt*(alpha_R(i)* &
                                                                                                vel_R(2))*(1._wp/dy(k + 1)) + &
                                                                                     0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dy(k + 1)), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k + 1, l) = rhs_vf(advxb + i - 1)%sf(j, k + 1, l) &
                                                                                - real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k + 1, l)*vel_R(2)*(1._wp/dy(k + 1)), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                            real(0.5_wp*dt*(alpha_R(i)* &
                                                                                            vel_R(2))*(1._wp/dy(k)) + &
                                                                                 0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dy(k)), kind=stp)

                                        $:GPU_ATOMIC(atomic='update')
                                        rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                            + real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_R(2)*(1._wp/dy(k)), kind=stp)
                                    end do
                                end if

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k + 1, l) = rhs_vf(momxb + 1)%sf(j, k + 1, l) + &
                                                                    real(0.5_wp*dt*(rho_R*(vel_R(2))**2.0 + &
                                                                                    pres_R + F_R)*(1._wp/dy(k + 1)) + &
                                                                         0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k + 1, l) = rhs_vf(momxb)%sf(j, k + 1, l) + &
                                                                real(0.5_wp*dt*rho_R*vel_R(2)*vel_R(1)*(1._wp/dy(k + 1)) + &
                                                                     0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k + 1, l) = rhs_vf(momxb + 2)%sf(j, k + 1, l) + &
                                                                    real(0.5_wp*dt*rho_R*vel_R(2)*vel_R(3)*(1._wp/dy(k + 1)) + &
                                                                         0.5_wp*dt*cfl*(rho_R*vel_R(3))*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k + 1, l) = rhs_vf(E_idx)%sf(j, k + 1, l) + &
                                                                real(0.5_wp*dt*(vel_R(2)*(E_R + &
                                                                                          pres_R + F_R))*(1._wp/dy(k + 1)) + &
                                                                     0.5_wp*dt*cfl*(E_R)*(1._wp/dy(k + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                                real(0.5_wp*dt*(rho_R*(vel_R(2))**2.0 + &
                                                                                pres_R + F_R)*(1._wp/dy(k)) + &
                                                                     0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dy(k)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*rho_R*vel_R(2)*vel_R(1)*(1._wp/dy(k)) + &
                                                                 0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dy(k)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) - &
                                                                real(0.5_wp*dt*rho_R*vel_R(2)*vel_R(3)*(1._wp/dy(k)) + &
                                                                     0.5_wp*dt*cfl*(rho_R*vel_R(3))*(1._wp/dy(k)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*(vel_R(2)*(E_R + &
                                                                                      pres_R + F_R))*(1._wp/dy(k)) + &
                                                                 0.5_wp*dt*cfl*(E_R)*(1._wp/dy(k)), kind=stp)

                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:endif
            end if
        elseif (idir == 3) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j,k,l,rho_L, rho_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, mu_L, mu_R, vel_L, vel_R, pres_L, pres_R, alpha_L, alpha_R, alpha_rho_L, alpha_rho_R, F_L, F_R, E_L, E_R, cfl, dvel_small, rho_sf_small, vflux_L_arr, vflux_R_arr]')
                do l = -1, p
                    do k = 0, n
                        do j = 0, m

                            if (viscous) then
                                vflux_L_arr = 0._wp
                                vflux_R_arr = 0._wp

                                #:if MFC_CASE_OPTIMIZATION
                                    #:if igr_order == 5
                                        !DIR$ unroll 6
                                    #:elif igr_order == 3
                                        !DIR$ unroll 4
                                    #:endif
                                #:endif
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = vidxb, vidxe

                                    !x-direction contributions
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_cons_vf(r)%sf(j + i, k, l + q)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dx(j)))*( &
                                                    q_cons_vf(momxb)%sf(j + 1, k, l + q)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb)%sf(j - 1, k, l + q)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dx(j)))*( &
                                                    q_cons_vf(momxb + 2)%sf(j + 1, k, l + q)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 2)%sf(j - 1, k, l + q)/rho_sf_small(-1))

                                    if (q > vidxb) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(3))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(-2._wp*dvel_small(1))/3._wp
                                    end if
                                    if (q < vidxe) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(3))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(-2._wp*dvel_small(1))/3._wp
                                    end if

                                    !y-direction contributions
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_cons_vf(r)%sf(j, k + i, l + q)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(2) = (1/(2._wp*dy(k)))*( &
                                                    q_cons_vf(momxb + 1)%sf(j, k + 1, l + q)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 1)%sf(j, k - 1, l + q)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dy(k)))*( &
                                                    q_cons_vf(momxb + 2)%sf(j, k + 1, l + q)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 2)%sf(j, k - 1, l + q)/rho_sf_small(-1))

                                    if (q > vidxb) then
                                        vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q + offxL)*(dvel_small(3))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(-2._wp*dvel_small(2))/3._wp
                                    end if
                                    if (q < vidxe) then
                                        vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q + offxR)*(dvel_small(3))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(-2._wp*dvel_small(2))/3._wp
                                    end if

                                    !z-direction contributions
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_cons_vf(r)%sf(j, k, l + i + q)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do
                                    dvel_small(1) = (1/(2._wp*dz(l)))*( &
                                                    q_cons_vf(momxb)%sf(j, k, l + 1 + q)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb)%sf(j, k, l - 1 + q)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dz(l)))*( &
                                                    q_cons_vf(momxb + 1)%sf(j, k, l + 1 + q)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 1)%sf(j, k, l - 1 + q)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dz(l)))*( &
                                                    q_cons_vf(momxb + 2)%sf(j, k, l + 1 + q)/rho_sf_small(1) - &
                                                    q_cons_vf(momxb + 2)%sf(j, k, l - 1 + q)/rho_sf_small(-1))
                                    if (q > vidxb) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q + offxL)*(dvel_small(1))
                                        vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q + offxL)*(dvel_small(2))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q + offxL)*(4._wp*dvel_small(3))/3._wp
                                    end if
                                    if (q < vidxe) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q + offxR)*(dvel_small(1))
                                        vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q + offxR)*(dvel_small(2))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q + offxR)*(4._wp*dvel_small(3))/3._wp
                                    end if
                                end do
                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_L(i) = 0._wp
                                alpha_rho_R(i) = 0._wp
                                alpha_L(i) = 0._wp
                                alpha_R(i) = 0._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = 0._wp
                                vel_R(i) = 0._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb + 1, vidxe
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = alpha_rho_L(i) + coeff_L(q + offxL)*q_cons_vf(i)%sf(j, k, l + q)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        alpha_L(i) = alpha_L(i) + coeff_L(q + offxL)*q_cons_vf(E_idx + i)%sf(j, k, l + q)
                                    end do
                                else
                                    alpha_L(1) = 1._wp
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = vel_L(i) + coeff_L(q + offxL)*q_cons_vf(momxb + i - 1)%sf(j, k, l + q)
                                end do
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb, vidxe - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_R(i) = alpha_rho_R(i) + coeff_R(q + offxR)*q_cons_vf(i)%sf(j, k, l + q)
                                end do

                                if (num_fluids > 1) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        alpha_R(i) = alpha_R(i) + coeff_R(q + offxR)*q_cons_vf(E_idx + i)%sf(j, k, l + q)
                                    end do
                                else
                                    alpha_R(1) = 1._wp
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_R(i) = vel_R(i) + coeff_R(q + offxR)*q_cons_vf(momxb + i - 1)%sf(j, k, l + q)
                                end do
                            end do

                            if (num_fluids > 1) then

                                alpha_L(num_fluids) = 1._wp
                                alpha_R(num_fluids) = 1._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    alpha_L(num_fluids) = alpha_L(num_fluids) - alpha_L(i)
                                    alpha_R(num_fluids) = alpha_R(num_fluids) - alpha_R(i)
                                end do
                            end if

                            rho_L = 0._wp; rho_R = 0._wp
                            gamma_L = 0._wp; gamma_R = 0._wp
                            pi_inf_L = 0._wp; pi_inf_R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho_L = rho_L + alpha_rho_L(i)
                                gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)

                                rho_R = rho_R + alpha_rho_R(i)
                                gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = vel_L(i)/rho_L
                                vel_R(i) = vel_R(i)/rho_R
                            end do

                            if (viscous) then
                                mu_L = 0._wp
                                mu_R = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    mu_L = alpha_L(i)/Res_igr(1, i) + mu_L
                                    mu_R = alpha_R(i)/Res_igr(1, i) + mu_R
                                end do

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l + 1) = rhs_vf(momxb)%sf(j, k, l + 1) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dz(l + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l + 1) = rhs_vf(E_idx)%sf(j, k, l + 1) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dz(l + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(1)*(1._wp/dz(l)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dz(l)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l + 1) = rhs_vf(momxb)%sf(j, k, l + 1) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dz(l + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l + 1) = rhs_vf(E_idx)%sf(j, k, l + 1) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dz(l + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(1)*(1._wp/dz(l)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dz(l)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l + 1) = rhs_vf(momxb + 1)%sf(j, k, l + 1) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(2)*(1._wp/dz(l + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l + 1) = rhs_vf(E_idx)%sf(j, k, l + 1) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(2)*vel_L(2)*(1._wp/dz(l + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(2)*(1._wp/dz(l)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(2)*vel_L(2)*(1._wp/dz(l)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l + 1) = rhs_vf(momxb + 1)%sf(j, k, l + 1) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(2)*(1._wp/dz(l + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l + 1) = rhs_vf(E_idx)%sf(j, k, l + 1) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(2)*vel_R(2)*(1._wp/dz(l + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(2)*(1._wp/dz(l)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(2)*vel_R(2)*(1._wp/dz(l)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k, l + 1) = rhs_vf(momxb + 2)%sf(j, k, l + 1) - &
                                                                    real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dz(l + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l + 1) = rhs_vf(E_idx)%sf(j, k, l + 1) - &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(3)*(1._wp/dz(l + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_L*vflux_L_arr(3)*(1._wp/dz(l)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_L*vflux_L_arr(3)*vel_L(3)*(1._wp/dz(l)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k, l + 1) = rhs_vf(momxb + 2)%sf(j, k, l + 1) - &
                                                                    real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dz(l + 1)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l + 1) = rhs_vf(E_idx)%sf(j, k, l + 1) - &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(3)*(1._wp/dz(l + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) + &
                                                                real(0.5_wp*dt*mu_R*vflux_R_arr(3)*(1._wp/dz(l)), kind=stp)
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                            real(0.5_wp*dt*mu_R*vflux_R_arr(3)*vel_R(3)*(1._wp/dz(l)), kind=stp)
                            end if

                            E_L = 0._wp; E_R = 0._wp
                            F_L = 0._wp; F_R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb + 1, vidxe
                                E_L = E_L + coeff_L(q + offxL)*q_cons_vf(E_idx)%sf(j, k, l + q)
                                F_L = F_L + coeff_L(q + offxL)*jac(j, k, l + q)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = vidxb, vidxe - 1
                                E_R = E_R + coeff_R(q + offxR)*q_cons_vf(E_idx)%sf(j, k, l + q)
                                F_R = F_R + coeff_R(q + offxR)*jac(j, k, l + q)
                            end do

                            call s_get_derived_states(E_L, gamma_L, pi_inf_L, rho_L, vel_L, &
                                                      E_R, gamma_R, pi_inf_R, rho_R, vel_R, &
                                                      pres_L, pres_R, cfl)

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k, l + 1) = rhs_vf(i)%sf(j, k, l + 1) + &
                                                            real(0.5_wp*dt*(alpha_rho_L(i)* &
                                                                            vel_L(3))*(1._wp/dz(l + 1)) - &
                                                                 0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dz(l + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*(alpha_rho_L(i)* &
                                                                        vel_L(3))*(1._wp/dz(l)) - &
                                                             0.5_wp*dt*cfl*(alpha_rho_L(i))*(1._wp/dz(l)), kind=stp)
                            end do

                            if (num_fluids > 1) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l + 1) = rhs_vf(advxb + i - 1)%sf(j, k, l + 1) + &
                                                                            real(0.5_wp*dt*(alpha_L(i)* &
                                                                                            vel_L(3))*(1._wp/dz(l + 1)) - &
                                                                                 0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dz(l + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l + 1) = rhs_vf(advxb + i - 1)%sf(j, k, l + 1) &
                                                                            - real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l + 1)*vel_L(3)*(1._wp/dz(l + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                        real(0.5_wp*dt*(alpha_L(i)* &
                                                                                        vel_L(3))*(1._wp/dz(l)) - &
                                                                             0.5_wp*dt*cfl*(alpha_L(i))*(1._wp/dz(l)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                        + real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_L(3)*(1._wp/dz(l)), kind=stp)
                                end do
                            end if

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 2)%sf(j, k, l + 1) = rhs_vf(momxb + 2)%sf(j, k, l + 1) + &
                                                                real(0.5_wp*dt*(rho_L*(vel_L(3))**2.0 + &
                                                                                pres_L + F_L)*(1._wp/dz(l + 1)) - &
                                                                     0.5_wp*dt*cfl*(rho_L*vel_L(3))*(1._wp/dz(l + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k, l + 1) = rhs_vf(momxb)%sf(j, k, l + 1) + &
                                                            real(0.5_wp*dt*rho_L*vel_L(1)*vel_L(3)*(1._wp/dz(l + 1)) - &
                                                                 0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dz(l + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k, l + 1) = rhs_vf(momxb + 1)%sf(j, k, l + 1) + &
                                                                real(0.5_wp*dt*rho_L*vel_L(2)*vel_L(3)*(1._wp/dz(l + 1)) - &
                                                                     0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dz(l + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k, l + 1) = rhs_vf(E_idx)%sf(j, k, l + 1) + &
                                                            real(0.5_wp*dt*(vel_L(3)*(E_L + &
                                                                                      pres_L + F_L))*(1._wp/dz(l + 1)) - &
                                                                 0.5_wp*dt*cfl*(E_L)*(1._wp/dz(l + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*(rho_L*(vel_L(3))**2.0 + &
                                                                            pres_L + F_L)*(1._wp/dz(l)) - &
                                                                 0.5_wp*dt*cfl*(rho_L*vel_L(3))*(1._wp/dz(l)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*rho_L*vel_L(1)*vel_L(3)*(1._wp/dz(l)) - &
                                                             0.5_wp*dt*cfl*(rho_L*vel_L(1))*(1._wp/dz(l)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*rho_L*vel_L(2)*vel_L(3)*(1._wp/dz(l)) - &
                                                                 0.5_wp*dt*cfl*(rho_L*vel_L(2))*(1._wp/dz(l)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*(vel_L(3)*(E_L + &
                                                                                  pres_L + F_L))*(1._wp/dz(l)) - &
                                                             0.5_wp*dt*cfl*(E_L)*(1._wp/dz(l)), kind=stp)

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k, l + 1) = rhs_vf(i)%sf(j, k, l + 1) + &
                                                            real(0.5_wp*dt*(alpha_rho_R(i)* &
                                                                            vel_R(3))*(1._wp/dz(l + 1)) + &
                                                                 0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dz(l + 1)), kind=stp)

                                $:GPU_ATOMIC(atomic='update')
                                rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*(alpha_rho_R(i)* &
                                                                        vel_R(3))*(1._wp/dz(l)) + &
                                                             0.5_wp*dt*cfl*(alpha_rho_R(i))*(1._wp/dz(l)), kind=stp)
                            end do

                            if (num_fluids > 1) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l + 1) = rhs_vf(advxb + i - 1)%sf(j, k, l + 1) + &
                                                                            real(0.5_wp*dt*(alpha_R(i)* &
                                                                                            vel_R(3))*(1._wp/dz(l + 1)) + &
                                                                                 0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dz(l + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l + 1) = rhs_vf(advxb + i - 1)%sf(j, k, l + 1) &
                                                                            - real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l + 1)*vel_R(3)*(1._wp/dz(l + 1)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) - &
                                                                        real(0.5_wp*dt*(alpha_R(i)* &
                                                                                        vel_R(3))*(1._wp/dz(l)) + &
                                                                             0.5_wp*dt*cfl*(alpha_R(i))*(1._wp/dz(l)), kind=stp)

                                    $:GPU_ATOMIC(atomic='update')
                                    rhs_vf(advxb + i - 1)%sf(j, k, l) = rhs_vf(advxb + i - 1)%sf(j, k, l) &
                                                                        + real(0.5_wp*dt*q_cons_vf(advxb + i - 1)%sf(j, k, l)*vel_R(3)*(1._wp/dz(l)), kind=stp)
                                end do
                            end if

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 2)%sf(j, k, l + 1) = rhs_vf(momxb + 2)%sf(j, k, l + 1) + &
                                                                real(0.5_wp*dt*(rho_R*(vel_R(3))**2.0 + &
                                                                                pres_R + F_R)*(1._wp/dz(l + 1)) + &
                                                                     0.5_wp*dt*cfl*(rho_R*vel_R(3))*(1._wp/dz(l + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k, l + 1) = rhs_vf(momxb)%sf(j, k, l + 1) + &
                                                            real(0.5_wp*dt*rho_R*vel_R(1)*vel_R(3)*(1._wp/dz(l + 1)) + &
                                                                 0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dz(l + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k, l + 1) = rhs_vf(momxb + 1)%sf(j, k, l + 1) + &
                                                                real(0.5_wp*dt*rho_R*vel_R(2)*vel_R(3)*(1._wp/dz(l + 1)) + &
                                                                     0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dz(l + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k, l + 1) = rhs_vf(E_idx)%sf(j, k, l + 1) + &
                                                            real(0.5_wp*dt*(vel_R(3)*(E_R + &
                                                                                      pres_R + F_R))*(1._wp/dz(l + 1)) + &
                                                                 0.5_wp*dt*cfl*(E_R)*(1._wp/dz(l + 1)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*(rho_R*(vel_R(3))**2.0 + &
                                                                            pres_R + F_R)*(1._wp/dz(l)) + &
                                                                 0.5_wp*dt*cfl*(rho_R*vel_R(3))*(1._wp/dz(l)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*rho_R*vel_R(1)*vel_R(3)*(1._wp/dz(l)) + &
                                                             0.5_wp*dt*cfl*(rho_R*vel_R(1))*(1._wp/dz(l)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) - &
                                                            real(0.5_wp*dt*rho_R*vel_R(2)*vel_R(3)*(1._wp/dz(l)) + &
                                                                 0.5_wp*dt*cfl*(rho_R*vel_R(2))*(1._wp/dz(l)), kind=stp)

                            $:GPU_ATOMIC(atomic='update')
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) - &
                                                        real(0.5_wp*dt*(vel_R(3)*(E_R + &
                                                                                  pres_R + F_R))*(1._wp/dz(l)) + &
                                                             0.5_wp*dt*cfl*(E_R)*(1._wp/dz(l)), kind=stp)

                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:endif
        end if

    end subroutine s_igr_riemann_solver

    subroutine s_get_derived_states(E_L, gamma_L, pi_inf_L, rho_L, vel_L, &
                                    E_R, gamma_R, pi_inf_R, rho_R, vel_R, &
                                    pres_L, pres_R, cfl)
        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: E_L, gamma_L, pi_inf_L, rho_L
        real(wp), intent(in) :: E_R, gamma_R, pi_inf_R, rho_R
        real(wp), dimension(num_dims), intent(in) :: vel_L, vel_R
        real(wp), intent(out) :: pres_L, pres_R, cfl

        real(wp) :: a_L, a_R

        if (num_dims == 2) then
            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1)**2._wp + vel_L(2)**2._wp))/gamma_L
            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1)**2._wp + vel_R(2)**2._wp))/gamma_R

            if (igr_pres_lim) then
                pres_L = max(pres_L, 0._wp)
                pres_R = max(pres_R, 0._wp)
            end if

            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L/gamma_L)/rho_L)
            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R/gamma_R)/rho_R)

            cfl = max(sqrt(vel_L(1)**2._wp + vel_L(2)**2._wp), &
                      sqrt(vel_R(1)**2._wp + vel_R(2)**2._wp)) + &
                  max(a_L, a_R)
        elseif (num_dims == 3) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1)**2._wp + vel_L(2)**2._wp + vel_L(3)**2._wp))/gamma_L
                pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1)**2._wp + vel_R(2)**2._wp + vel_R(3)**2._wp))/gamma_R

                if (igr_pres_lim) then
                    pres_L = max(pres_L, 0._wp)
                    pres_R = max(pres_R, 0._wp)
                end if

                a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L/gamma_L)/rho_L)
                a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R/gamma_R)/rho_R)

                cfl = max(sqrt(vel_L(1)**2._wp + vel_L(2)**2._wp + vel_L(3)**2._wp), &
                          sqrt(vel_R(1)**2._wp + vel_R(2)**2._wp + vel_R(3)**2._wp)) + &
                      max(a_L, a_R)
            #:endif
        end if

    end subroutine s_get_derived_states

    subroutine s_igr_flux_add(q_cons_vf, rhs_vf, flux_vf, idir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_cons_vf, flux_vf, rhs_vf

        integer, intent(in) :: idir

        if (idir == 1) then
            $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = 1._wp/dx(j)* &
                                                    (flux_vf(i)%sf(j - 1, k, l) &
                                                     - flux_vf(i)%sf(j, k, l))
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (idir == 2) then
            $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = &
                                rhs_vf(i)%sf(j, k, l) + 1._wp/dy(k)* &
                                (flux_vf(i)%sf(j, k - 1, l) &
                                 - flux_vf(i)%sf(j, k, l))
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (idir == 3) then
            $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = &
                                rhs_vf(i)%sf(j, k, l) + 1._wp/dz(l)* &
                                (flux_vf(i)%sf(j, k, l - 1) &
                                 - flux_vf(i)%sf(j, k, l))
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_igr_flux_add

    subroutine s_finalize_igr_module()

        if (viscous) then
            @:DEALLOCATE(Res_igr)
        end if

#ifndef __NVCOMPILER_GPU_UNIFIED_MEM
        @:DEALLOCATE(jac, jac_rhs)

        if (igr_iter_solver == 1) then ! Jacobi iteration
            @:DEALLOCATE(jac_old)
        end if
#else
        if (nv_uvm_temp_on_gpu(1) == 1) then
            @:DEALLOCATE(jac)
        else
            nullify (jac)
            deallocate (jac_host)
        end if

        if (nv_uvm_temp_on_gpu(2) == 1) then
            @:DEALLOCATE(jac_rhs)
        else
            nullify (jac_rhs)
            deallocate (jac_rhs_host)
        end if

        if (igr_iter_solver == 1) then ! Jacobi iteration
            if (nv_uvm_temp_on_gpu(3) == 1) then
                @:DEALLOCATE(jac_old)
            else
                nullify (jac_old)
                deallocate (jac_old_host)
            end if
        end if
#endif

        #:if not MFC_CASE_OPTIMIZATION
            @:DEALLOCATE(coeff_L, coeff_R)
        #:endif

    end subroutine s_finalize_igr_module

end module m_igr
