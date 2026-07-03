!>
!! @file
!! @brief Contains module m_riemann_state

!> @brief Shared Riemann-solver module state and the per-sweep setup, state-buffer population, viscous source flux, and finalization
!! helpers
#:include 'case.fpp'
#:include 'macros.fpp'

module m_riemann_state

    use m_derived_types
    use m_global_parameters
    use m_constants, only: riemann_solver_hll, riemann_solver_hlld
    use m_hb_function

    implicit none

    !> The cell-boundary values of the fluxes (src - source) that are computed through the chosen Riemann problem solver, and the
    !! direct evaluation of source terms, by using the left and right states given in qK_prim_rs_vf, dqK_prim_ds_vf where ds = dx,
    !! dy or dz.
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: flux_rsx_vf, flux_src_rsx_vf
    $:GPU_DECLARE(create='[flux_rsx_vf, flux_src_rsx_vf]')
    !> @}

    !> The cell-boundary values of the geometrical source flux that are computed through the chosen Riemann problem solver by using
    !! the left and right states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: flux_gsrc_rsx_vf
    $:GPU_DECLARE(create='[flux_gsrc_rsx_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: nc_iface_vel_rsx_vf
    $:GPU_DECLARE(create='[nc_iface_vel_rsx_vf]')

    !> Dual-pass HLLD second flux set: the hat_R-anchored fluxes (and, for axisymmetric runs, the hat_R interface velocities)
    !! written by the same fused solve that fills flux_rsx / nc_iface_vel_rsx with the hat_L-anchored values. Allocated only when
    !! hypo_nc_mode_dual_pass.
    real(wp), allocatable, dimension(:,:,:,:) :: flux_hatR_rsx_vf
    $:GPU_DECLARE(create='[flux_hatR_rsx_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: nc_iface_vel_hatR_rsx_vf
    $:GPU_DECLARE(create='[nc_iface_vel_hatR_rsx_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: flux_gsrc_hatR_rsx_vf
    $:GPU_DECLARE(create='[flux_gsrc_hatR_rsx_vf]')
    !> @}

    ! Cell-boundary velocity from Riemann solution; used for source flux

    real(wp), allocatable, dimension(:,:,:,:) :: vel_src_rsx_vf
    $:GPU_DECLARE(create='[vel_src_rsx_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: mom_sp_rsx_vf
    $:GPU_DECLARE(create='[mom_sp_rsx_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: Re_avg_rsx_vf
    $:GPU_DECLARE(create='[Re_avg_rsx_vf]')

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1, is2, is3
    type(int_bounds_info) :: isx, isy, isz
    !> @}

    $:GPU_DECLARE(create='[is1, is2, is3, isx, isy, isz]')

    real(wp), allocatable, dimension(:) :: Gs_rs
    $:GPU_DECLARE(create='[Gs_rs]')

    real(wp), allocatable, dimension(:,:) :: Res_gs
    $:GPU_DECLARE(create='[Res_gs]')

contains

    !> Dispatch to the subroutines that are utilized to compute the viscous source fluxes for either Cartesian or cylindrical
    !! geometries. For more information please refer to: 1) s_compute_cartesian_viscous_source_flux 2)
    !! s_compute_cylindrical_viscous_source_flux
    subroutine s_compute_viscous_source_flux(velL_vf, dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, velR_vf, dvelR_dx_vf, dvelR_dy_vf, &
        & dvelR_dz_vf, flux_src_vf, q_prim_vf, norm_dir, ix, iy, iz)

        type(scalar_field), dimension(num_vels), intent(in) :: velL_vf, velR_vf, dvelL_dx_vf, dvelR_dx_vf, dvelL_dy_vf, &
             & dvelR_dy_vf, dvelL_dz_vf, dvelR_dz_vf

        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        if (grid_geometry == 3) then
            call s_compute_cylindrical_viscous_source_flux(velL_vf, dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, velR_vf, dvelR_dx_vf, &
                & dvelR_dy_vf, dvelR_dz_vf, flux_src_vf, q_prim_vf, norm_dir, ix, iy, iz)
        else
            call s_compute_cartesian_viscous_source_flux(dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, dvelR_dx_vf, dvelR_dy_vf, &
                & dvelR_dz_vf, flux_src_vf, q_prim_vf, norm_dir)
        end if

    end subroutine s_compute_viscous_source_flux

    !> Populate the left and right Riemann state variable buffers based on boundary conditions
    subroutine s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, &
        & qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        integer, intent(in)               :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz
        integer                           :: i, j, k, l  !< Generic loop iterator

        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1._wp, 0._wp, 0._wp/)
        else if (norm_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
            dir_idx = (/2, 1, 3/); dir_flg = (/0._wp, 1._wp, 0._wp/)
        else
            is1 = iz; is2 = iy; is3 = ix
            dir_idx = (/3, 1, 2/); dir_flg = (/0._wp, 0._wp, 1._wp/)
        end if

        $:GPU_UPDATE(device='[is1, is2, is3]')

        if (elasticity) then
            ! dir_idx_tau(1:3) = (nn, nt, nt2): face-normal stress row for wave speeds and momentum flux. stress_perm(1:n_stress) =
            ! full tensor permutation mapping F_HLL local basis index -> physical storage index. Local order: (nn, nt, tt, nt2,
            ! t1t2, t2t2). In 2D only entries 1-3 are used.
            if (norm_dir == 1) then
                dir_idx_tau = (/1, 2, 4/)
                stress_perm = (/1, 2, 3, 4, 5, 6/)
            else if (norm_dir == 2) then
                dir_idx_tau = (/3, 2, 5/)
                stress_perm = (/3, 2, 1, 5, 4, 6/)
            else
                dir_idx_tau = (/6, 4, 5/)
                stress_perm = (/6, 4, 1, 5, 2, 3/)
            end if
        end if

        isx = ix; isy = iy; isz = iz
        ! for stuff in the same module
        $:GPU_UPDATE(device='[isx, isy, isz]')
        ! for stuff in different modules
        $:GPU_UPDATE(device='[dir_idx, dir_flg, dir_idx_tau, stress_perm]')

        ! Population of Buffers in x-direction
        if (norm_dir == 1) then
            if (bc_x%beg == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsx_vf(-1, k, l, i) = qR_prim_rsx_vf(0, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end
                                dqL_prim_dx_vf(i)%sf(-1, k, l) = dqR_prim_dx_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    if (n > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end
                                    dqL_prim_dy_vf(i)%sf(-1, k, l) = dqR_prim_dy_vf(i)%sf(0, k, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        if (p > 0) then
                            $:GPU_PARALLEL_LOOP(collapse=3)
                            do i = eqn_idx%mom%beg, eqn_idx%mom%end
                                do l = isz%beg, isz%end
                                    do k = isy%beg, isy%end
                                        dqL_prim_dz_vf(i)%sf(-1, k, l) = dqR_prim_dz_vf(i)%sf(0, k, l)
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    end if
                end if
            end if

            if (bc_x%end == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at end
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsx_vf(m + 1, k, l, i) = qL_prim_rsx_vf(m, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end
                                dqR_prim_dx_vf(i)%sf(m + 1, k, l) = dqL_prim_dx_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    if (n > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end
                                    dqR_prim_dy_vf(i)%sf(m + 1, k, l) = dqL_prim_dy_vf(i)%sf(m, k, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        if (p > 0) then
                            $:GPU_PARALLEL_LOOP(collapse=3)
                            do i = eqn_idx%mom%beg, eqn_idx%mom%end
                                do l = isz%beg, isz%end
                                    do k = isy%beg, isy%end
                                        dqR_prim_dz_vf(i)%sf(m + 1, k, l) = dqL_prim_dz_vf(i)%sf(m, k, l)
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    end if
                end if
            end if
            ! END: Population of Buffers in x-direction

            ! Population of Buffers in y-direction
        else if (norm_dir == 2) then
            if (bc_y%beg == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsx_vf(k, -1, l, i) = qR_prim_rsx_vf(k, 0, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, -1, l) = dqR_prim_dx_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, -1, l) = dqR_prim_dy_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    if (p > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            do l = isz%beg, isz%end
                                do j = isx%beg, isx%end
                                    dqL_prim_dz_vf(i)%sf(j, -1, l) = dqR_prim_dz_vf(i)%sf(j, 0, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                end if
            end if

            if (bc_y%end == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at end
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsx_vf(k, n + 1, l, i) = qL_prim_rsx_vf(k, n, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, n + 1, l) = dqL_prim_dx_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, n + 1, l) = dqL_prim_dy_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    if (p > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            do l = isz%beg, isz%end
                                do j = isx%beg, isx%end
                                    dqR_prim_dz_vf(i)%sf(j, n + 1, l) = dqL_prim_dz_vf(i)%sf(j, n, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                end if
            end if
            ! END: Population of Buffers in y-direction

            ! Population of Buffers in z-direction
        else
            if (bc_z%beg == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                            qL_prim_rsx_vf(l, k, -1, i) = qR_prim_rsx_vf(l, k, 0, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, k, -1) = dqR_prim_dx_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, k, -1) = dqR_prim_dy_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dz_vf(i)%sf(j, k, -1) = dqR_prim_dz_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if

            if (bc_z%end == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at end
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                            qR_prim_rsx_vf(l, k, p + 1, i) = qL_prim_rsx_vf(l, k, p, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, k, p + 1) = dqL_prim_dx_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, k, p + 1) = dqL_prim_dy_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dz_vf(i)%sf(j, k, p + 1) = dqL_prim_dz_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        end if
        ! END: Population of Buffers in z-direction

    end subroutine s_populate_riemann_states_variables_buffers

    !> Set up the chosen Riemann solver algorithm for the current direction
    subroutine s_initialize_riemann_solver(flux_src_vf, norm_dir)

        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        integer, intent(in)                                    :: norm_dir
        integer                                                :: i, j, k, l  !< Generic loop iterators

        ! Reshaping Inputted Data in x-direction

        if (norm_dir == 1) then
            if (viscous .or. (surface_tension)) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%mom%beg, eqn_idx%E
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                flux_src_vf(i)%sf(j, k, l) = 0._wp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%E, eqn_idx%species%end
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                if (i == eqn_idx%E .or. i >= eqn_idx%species%beg) then
                                    flux_src_vf(i)%sf(j, k, l) = 0._wp
                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (qbmm) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, 4
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end + 1
                                mom_sp_rsx_vf(j, k, l, i) = mom_sp(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            ! Reshaping Inputted Data in y-direction
        else if (norm_dir == 2) then
            if (viscous .or. (surface_tension)) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%mom%beg, eqn_idx%E
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_src_vf(i)%sf(k, j, l) = 0._wp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%E, eqn_idx%species%end
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                if (i == eqn_idx%E .or. i >= eqn_idx%species%beg) then
                                    flux_src_vf(i)%sf(k, j, l) = 0._wp
                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (qbmm) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, 4
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end + 1
                                mom_sp_rsx_vf(k, j, l, i) = mom_sp(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            ! Reshaping Inputted Data in z-direction
        else
            if (viscous .or. (surface_tension)) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%mom%beg, eqn_idx%E
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                flux_src_vf(i)%sf(l, k, j) = 0._wp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%E, eqn_idx%species%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                if (i == eqn_idx%E .or. i >= eqn_idx%species%beg) then
                                    flux_src_vf(i)%sf(l, k, j) = 0._wp
                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (qbmm) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, 4
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end + 1
                                mom_sp_rsx_vf(l, k, j, i) = mom_sp(i)%sf(l, k, j)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

    end subroutine s_initialize_riemann_solver

    !> Compute cylindrical viscous source flux contributions for momentum and energy
    subroutine s_compute_cylindrical_viscous_source_flux(velL_vf, dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, velR_vf, dvelR_dx_vf, &
        & dvelR_dy_vf, dvelR_dz_vf, flux_src_vf, q_prim_vf, norm_dir, ix, iy, iz)

        type(scalar_field), dimension(num_dims), intent(in)    :: velL_vf, velR_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dx_vf, dvelR_dx_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dy_vf, dvelR_dy_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dz_vf, dvelR_dz_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        ! Local variables

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: avg_v_int     !< Averaged interface velocity (\f$v_x, v_y, v_z\f$) (grid directions).
            real(wp), dimension(3) :: avg_dvdx_int  !< Averaged interface \f$\partial v_i/\partial x\f$ (grid dir 1).
            real(wp), dimension(3) :: avg_dvdy_int  !< Averaged interface \f$\partial v_i/\partial y\f$ (grid dir 2).
            real(wp), dimension(3) :: avg_dvdz_int  !< Averaged interface \f$\partial v_i/\partial z\f$ (grid dir 3).
            real(wp), dimension(3) :: vel_src_int   !< Interface velocity (\f$v_1,v_2,v_3\f$) (grid directions) for viscous work.

            !> Shear stress vector (\f$\sigma_{N1}, \sigma_{N2}, \sigma_{N3}\f$) on N-face (grid directions).
            real(wp), dimension(3) :: stress_vector_shear
        #:else
            real(wp), dimension(num_dims) :: avg_v_int     !< Averaged interface velocity (\f$v_x, v_y, v_z\f$) (grid directions).
            real(wp), dimension(num_dims) :: avg_dvdx_int  !< Averaged interface \f$\partial v_i/\partial x\f$ (grid dir 1).
            real(wp), dimension(num_dims) :: avg_dvdy_int  !< Averaged interface \f$\partial v_i/\partial y\f$ (grid dir 2).
            real(wp), dimension(num_dims) :: avg_dvdz_int  !< Averaged interface \f$\partial v_i/\partial z\f$ (grid dir 3).
            !> Interface velocity (\f$v_1,v_2,v_3\f$) (grid directions) for viscous work.
            real(wp), dimension(num_dims) :: vel_src_int
            !> Shear stress vector (\f$\sigma_{N1}, \sigma_{N2}, \sigma_{N3}\f$) on N-face (grid directions).
            real(wp), dimension(num_dims) :: stress_vector_shear
        #:endif
        real(wp) :: stress_normal_bulk  !< Normal bulk stress component \f$\sigma_{NN}\f$ on N-face.
        real(wp) :: Re_s, Re_b  !< Effective interface shear and bulk Reynolds numbers.
        real(wp) :: r_eff  !< Effective radius at interface for cylindrical terms.
        real(wp) :: div_v_term_const  !< Common term \f$-(2/3)(\nabla \cdot \mathbf{v}) / \text{Re}_s\f$ for shear stress diagonal.
        real(wp) :: divergence_cyl  !< Full divergence \f$\nabla \cdot \mathbf{v}\f$ in cylindrical coordinates.
        integer :: j, k, l  !< Loop iterators for \f$x, y, z\f$ grid directions.
        integer :: i_vel  !< Loop iterator for velocity components.
        integer :: idx_rp(3)  !< Indices \f$(j,k,l)\f$ of 'right' point for averaging.
        real(wp) :: gamma_dot, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz
        real(wp), dimension(2) :: Re_nn
        real(wp), dimension(num_fluids) :: alpha_avg
        integer :: fl

        $:GPU_PARALLEL_LOOP(collapse=3, private='[idx_rp, avg_v_int, avg_dvdx_int, avg_dvdy_int, avg_dvdz_int, Re_s, Re_b, &
                            & vel_src_int, r_eff, divergence_cyl, stress_vector_shear, stress_normal_bulk, div_v_term_const, &
                            & gamma_dot, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz, Re_nn, alpha_avg, fl]')
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end
                    ! Determine indices for the 'right' state for averaging across the interface
                    idx_rp = [j, k, l]
                    idx_rp(norm_dir) = idx_rp(norm_dir) + 1

                    ! Average velocities and their derivatives at the interface For cylindrical: x-dir ~ axial (z_cyl), y-dir ~
                    ! radial (r_cyl), z-dir ~ azimuthal (theta_cyl)
                    $:GPU_LOOP(parallelism='[seq]')
                    do i_vel = 1, num_dims
                        avg_v_int(i_vel) = 0.5_wp*(velL_vf(i_vel)%sf(j, k, l) + velR_vf(i_vel)%sf(idx_rp(1), idx_rp(2), idx_rp(3)))

                        avg_dvdx_int(i_vel) = 0.5_wp*(dvelL_dx_vf(i_vel)%sf(j, k, l) + dvelR_dx_vf(i_vel)%sf(idx_rp(1), &
                                     & idx_rp(2), idx_rp(3)))
                        if (num_dims > 1) then
                            avg_dvdy_int(i_vel) = 0.5_wp*(dvelL_dy_vf(i_vel)%sf(j, k, l) + dvelR_dy_vf(i_vel)%sf(idx_rp(1), &
                                         & idx_rp(2), idx_rp(3)))
                        else
                            avg_dvdy_int(i_vel) = 0.0_wp
                        end if
                        if (num_dims > 2) then
                            avg_dvdz_int(i_vel) = 0.5_wp*(dvelL_dz_vf(i_vel)%sf(j, k, l) + dvelR_dz_vf(i_vel)%sf(idx_rp(1), &
                                         & idx_rp(2), idx_rp(3)))
                        else
                            avg_dvdz_int(i_vel) = 0.0_wp
                        end if
                    end do

                    ! Non-Newtonian effective shear rate from grid-direction strain components.
                    ! NOTE: curvature corrections to gamma_dot (e.g. hoop strain) are not included
                    ! here - a documented first-version limitation. The Reynolds override and the
                    ! Newtonian path below are exact.
                    if (any_non_newtonian) then
                        D_xx = avg_dvdx_int(1); D_yy = 0._wp; D_zz = 0._wp
                        D_xy = 0._wp; D_xz = 0._wp; D_yz = 0._wp
                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                            if (num_dims > 1) then
                                D_yy = avg_dvdy_int(2)
                                D_xy = 0.5_wp*(avg_dvdy_int(1) + avg_dvdx_int(2))
                            end if
                        #:endif
                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                            if (num_dims > 2) then
                                D_zz = avg_dvdz_int(3)
                                D_xz = 0.5_wp*(avg_dvdz_int(1) + avg_dvdx_int(3))
                                D_yz = 0.5_wp*(avg_dvdz_int(2) + avg_dvdy_int(3))
                            end if
                        #:endif
                        gamma_dot = f_compute_shear_rate_from_components(D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)
                        do fl = 1, num_fluids
                            alpha_avg(fl) = 0.5_wp*(q_prim_vf(eqn_idx%adv%beg + fl - 1)%sf(j, k, &
                                      & l) + q_prim_vf(eqn_idx%adv%beg + fl - 1)%sf(idx_rp(1), idx_rp(2), idx_rp(3)))
                            ! Raw cell-centered alphas can under/overshoot near interfaces; clamp to [0,1]
                            alpha_avg(fl) = min(max(alpha_avg(fl), 0._wp), 1._wp)
                        end do
                        call s_compute_mixture_inv_re(alpha_avg, gamma_dot, Res_gs, Re_nn)
                    end if

                    ! Get Re numbers and interface velocity for viscous work
                    select case (norm_dir)
                    case (1)  ! x-face (axial face in z_cyl direction)
                        if (any_non_newtonian) then
                            Re_s = Re_nn(1)
                            Re_b = Re_nn(2)
                        else
                            Re_s = Re_avg_rsx_vf(j, k, l, 1)
                            Re_b = Re_avg_rsx_vf(j, k, l, 2)
                        end if
                        vel_src_int = vel_src_rsx_vf(j, k, l,1:num_dims)
                        r_eff = y_cc(k)
                    case (2)  ! y-face (radial face in r_cyl direction)
                        if (any_non_newtonian) then
                            Re_s = Re_nn(1)
                            Re_b = Re_nn(2)
                        else
                            Re_s = Re_avg_rsx_vf(j, k, l, 1)
                            Re_b = Re_avg_rsx_vf(j, k, l, 2)
                        end if
                        vel_src_int = vel_src_rsx_vf(j, k, l,1:num_dims)
                        r_eff = y_cb(k)
                    case (3)  ! z-face (azimuthal face in theta_cyl direction)
                        if (any_non_newtonian) then
                            Re_s = Re_nn(1)
                            Re_b = Re_nn(2)
                        else
                            Re_s = Re_avg_rsx_vf(j, k, l, 1)
                            Re_b = Re_avg_rsx_vf(j, k, l, 2)
                        end if
                        vel_src_int = vel_src_rsx_vf(j, k, l,1:num_dims)
                        r_eff = y_cc(k)
                    end select

                    ! Divergence in cylindrical coordinates (vx=vz_cyl, vy=vr_cyl, vz=vtheta_cyl)
                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                        divergence_cyl = avg_dvdx_int(1) + avg_dvdy_int(2) + avg_v_int(2)/r_eff
                        if (num_dims > 2) then
                            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                divergence_cyl = divergence_cyl + avg_dvdz_int(3)/r_eff
                            #:endif
                        end if
                    #:endif

                    stress_vector_shear = 0.0_wp
                    stress_normal_bulk = 0.0_wp

                    if (shear_stress) then
                        div_v_term_const = -(2.0_wp/3.0_wp)*divergence_cyl/Re_s

                        select case (norm_dir)
                        case (1)  ! X-face (axial normal, z_cyl)
                            stress_vector_shear(1) = (2.0_wp*avg_dvdx_int(1))/Re_s + div_v_term_const
                            if (num_dims > 1) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    stress_vector_shear(2) = (avg_dvdy_int(1) + avg_dvdx_int(2))/Re_s
                                #:endif
                            end if
                            if (num_dims > 2) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    stress_vector_shear(3) = (avg_dvdz_int(1)/r_eff + avg_dvdx_int(3))/Re_s
                                #:endif
                            end if
                        case (2)  ! Y-face (radial normal, r_cyl)
                            if (num_dims > 1) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    stress_vector_shear(1) = (avg_dvdy_int(1) + avg_dvdx_int(2))/Re_s
                                    stress_vector_shear(2) = (2.0_wp*avg_dvdy_int(2))/Re_s + div_v_term_const
                                    if (num_dims > 2) then
                                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                            stress_vector_shear(3) = (avg_dvdz_int(2)/r_eff - avg_v_int(3)/r_eff + avg_dvdy_int(3) &
                                                                & )/Re_s
                                        #:endif
                                    end if
                                #:endif
                            else
                                stress_vector_shear(1) = (2.0_wp*avg_dvdx_int(1))/Re_s + div_v_term_const
                            end if
                        case (3)  ! Z-face (azimuthal normal, theta_cyl)
                            if (num_dims > 2) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    stress_vector_shear(1) = (avg_dvdz_int(1)/r_eff + avg_dvdx_int(3))/Re_s
                                    stress_vector_shear(2) = (avg_dvdz_int(2)/r_eff - avg_v_int(3)/r_eff + avg_dvdy_int(3))/Re_s
                                    stress_vector_shear(3) = (2.0_wp*(avg_dvdz_int(3)/r_eff + avg_v_int(2)/r_eff))/Re_s &
                                                        & + div_v_term_const
                                #:endif
                            end if
                        end select

                        $:GPU_LOOP(parallelism='[seq]')
                        do i_vel = 1, num_dims
                            flux_src_vf(eqn_idx%mom%beg + i_vel - 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + i_vel - 1)%sf(j, &
                                        & k, l) - stress_vector_shear(i_vel)
                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                        & l) - vel_src_int(i_vel)*stress_vector_shear(i_vel)
                        end do
                    end if

                    if (bulk_stress) then
                        stress_normal_bulk = divergence_cyl/Re_b

                        flux_src_vf(eqn_idx%mom%beg + norm_dir - 1)%sf(j, k, &
                                    & l) = flux_src_vf(eqn_idx%mom%beg + norm_dir - 1)%sf(j, k, l) - stress_normal_bulk
                        flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                    & l) - vel_src_int(norm_dir)*stress_normal_bulk
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_cylindrical_viscous_source_flux

    !> Compute Cartesian viscous source flux contributions for momentum and energy
    subroutine s_compute_cartesian_viscous_source_flux(dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, dvelR_dx_vf, dvelR_dy_vf, &
        & dvelR_dz_vf, flux_src_vf, q_prim_vf, norm_dir)

        ! Arguments
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dx_vf, dvelR_dx_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dy_vf, dvelR_dy_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dz_vf, dvelR_dz_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        integer, intent(in)                                    :: norm_dir

        ! Local variables

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3, 3) :: vel_grad_avg          !< Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
            real(wp), dimension(3, 3) :: current_tau_shear     !< Current shear stress tensor.
            real(wp), dimension(3, 3) :: current_tau_bulk      !< Current bulk stress tensor.
            real(wp), dimension(3)    :: vel_src_at_interface  !< Interface velocities (u,v,w) for viscous work.
        #:else
            real(wp), dimension(num_dims, num_dims) :: vel_grad_avg  !< Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
            real(wp), dimension(num_dims, num_dims) :: current_tau_shear  !< Current shear stress tensor.
            real(wp), dimension(num_dims, num_dims) :: current_tau_bulk  !< Current bulk stress tensor.
            real(wp), dimension(num_dims)           :: vel_src_at_interface  !< Interface velocities (u,v,w) for viscous work.
        #:endif
        integer, dimension(3)           :: idx_right_phys  !< Physical (j,k,l) indices for right state.
        real(wp)                        :: Re_shear        !< Interface shear Reynolds number.
        real(wp)                        :: Re_bulk         !< Interface bulk Reynolds number.
        integer                         :: j_loop          !< Physical x-index loop iterator.
        integer                         :: k_loop          !< Physical y-index loop iterator.
        integer                         :: l_loop          !< Physical z-index loop iterator.
        integer                         :: i_dim           !< Generic dimension/component iterator.
        integer                         :: vel_comp_idx    !< Velocity component iterator (1=u, 2=v, 3=w).
        real(wp)                        :: divergence_v    !< Velocity divergence at interface.
        real(wp)                        :: gamma_dot, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz
        real(wp), dimension(2)          :: Re_nn
        real(wp), dimension(num_fluids) :: alpha_avg
        integer                         :: fl

        $:GPU_PARALLEL_LOOP(collapse=3, private='[idx_right_phys, vel_grad_avg, current_tau_shear, current_tau_bulk, &
                            & vel_src_at_interface, Re_shear, Re_bulk, divergence_v, i_dim, vel_comp_idx, gamma_dot, D_xx, D_yy, &
                            & D_zz, D_xy, D_xz, D_yz, Re_nn, alpha_avg, fl]')
        do l_loop = isz%beg, isz%end
            do k_loop = isy%beg, isy%end
                do j_loop = isx%beg, isx%end
                    idx_right_phys(1) = j_loop
                    idx_right_phys(2) = k_loop
                    idx_right_phys(3) = l_loop
                    idx_right_phys(norm_dir) = idx_right_phys(norm_dir) + 1

                    vel_grad_avg = 0.0_wp
                    do vel_comp_idx = 1, num_dims
                        vel_grad_avg(vel_comp_idx, 1) = 0.5_wp*(dvelL_dx_vf(vel_comp_idx)%sf(j_loop, k_loop, &
                                     & l_loop) + dvelR_dx_vf(vel_comp_idx)%sf(idx_right_phys(1), idx_right_phys(2), &
                                     & idx_right_phys(3)))
                        if (num_dims > 1) then
                            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                vel_grad_avg(vel_comp_idx, 2) = 0.5_wp*(dvelL_dy_vf(vel_comp_idx)%sf(j_loop, k_loop, &
                                             & l_loop) + dvelR_dy_vf(vel_comp_idx)%sf(idx_right_phys(1), idx_right_phys(2), &
                                             & idx_right_phys(3)))
                            #:endif
                        end if
                        if (num_dims > 2) then
                            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                vel_grad_avg(vel_comp_idx, 3) = 0.5_wp*(dvelL_dz_vf(vel_comp_idx)%sf(j_loop, k_loop, &
                                             & l_loop) + dvelR_dz_vf(vel_comp_idx)%sf(idx_right_phys(1), idx_right_phys(2), &
                                             & idx_right_phys(3)))
                            #:endif
                        end if
                    end do

                    if (any_non_newtonian) then
                        D_xx = vel_grad_avg(1, 1); D_yy = 0._wp; D_zz = 0._wp
                        D_xy = 0._wp; D_xz = 0._wp; D_yz = 0._wp
                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                            if (num_dims > 1) then
                                D_yy = vel_grad_avg(2, 2)
                                D_xy = 0.5_wp*(vel_grad_avg(1, 2) + vel_grad_avg(2, 1))
                            end if
                        #:endif
                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                            if (num_dims > 2) then
                                D_zz = vel_grad_avg(3, 3)
                                D_xz = 0.5_wp*(vel_grad_avg(1, 3) + vel_grad_avg(3, 1))
                                D_yz = 0.5_wp*(vel_grad_avg(2, 3) + vel_grad_avg(3, 2))
                            end if
                        #:endif
                        gamma_dot = f_compute_shear_rate_from_components(D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)
                        do fl = 1, num_fluids
                            alpha_avg(fl) = 0.5_wp*(q_prim_vf(eqn_idx%adv%beg + fl - 1)%sf(j_loop, k_loop, &
                                      & l_loop) + q_prim_vf(eqn_idx%adv%beg + fl - 1)%sf(idx_right_phys(1), idx_right_phys(2), &
                                      & idx_right_phys(3)))
                            ! Raw cell-centered alphas can under/overshoot near interfaces; clamp to [0,1]
                            alpha_avg(fl) = min(max(alpha_avg(fl), 0._wp), 1._wp)
                        end do
                        call s_compute_mixture_inv_re(alpha_avg, gamma_dot, Res_gs, Re_nn)
                    end if

                    divergence_v = 0.0_wp
                    do i_dim = 1, num_dims
                        divergence_v = divergence_v + vel_grad_avg(i_dim, i_dim)
                    end do

                    vel_src_at_interface = 0.0_wp
                    if (norm_dir == 1) then
                        if (any_non_newtonian) then
                            Re_shear = Re_nn(1)
                            Re_bulk = Re_nn(2)
                        else
                            Re_shear = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 1)
                            Re_bulk = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 2)
                        end if
                        do i_dim = 1, num_dims
                            vel_src_at_interface(i_dim) = vel_src_rsx_vf(j_loop, k_loop, l_loop, i_dim)
                        end do
                    else if (norm_dir == 2) then
                        if (any_non_newtonian) then
                            Re_shear = Re_nn(1)
                            Re_bulk = Re_nn(2)
                        else
                            Re_shear = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 1)
                            Re_bulk = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 2)
                        end if
                        do i_dim = 1, num_dims
                            vel_src_at_interface(i_dim) = vel_src_rsx_vf(j_loop, k_loop, l_loop, i_dim)
                        end do
                    else
                        if (any_non_newtonian) then
                            Re_shear = Re_nn(1)
                            Re_bulk = Re_nn(2)
                        else
                            Re_shear = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 1)
                            Re_bulk = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 2)
                        end if
                        do i_dim = 1, num_dims
                            vel_src_at_interface(i_dim) = vel_src_rsx_vf(j_loop, k_loop, l_loop, i_dim)
                        end do
                    end if

                    if (shear_stress) then
                        ! current_tau_shear = 0.0_wp
                        call s_calculate_shear_stress_tensor(vel_grad_avg, Re_shear, divergence_v, current_tau_shear)

                        do i_dim = 1, num_dims
                            flux_src_vf(eqn_idx%mom%beg + i_dim - 1)%sf(j_loop, k_loop, &
                                        & l_loop) = flux_src_vf(eqn_idx%mom%beg + i_dim - 1)%sf(j_loop, k_loop, &
                                        & l_loop) - current_tau_shear(norm_dir, i_dim)

                            flux_src_vf(eqn_idx%E)%sf(j_loop, k_loop, l_loop) = flux_src_vf(eqn_idx%E)%sf(j_loop, k_loop, &
                                        & l_loop) - vel_src_at_interface(i_dim)*current_tau_shear(norm_dir, i_dim)
                        end do
                    end if

                    if (bulk_stress) then
                        ! current_tau_bulk = 0.0_wp
                        call s_calculate_bulk_stress_tensor(Re_bulk, divergence_v, current_tau_bulk)

                        do i_dim = 1, num_dims
                            flux_src_vf(eqn_idx%mom%beg + i_dim - 1)%sf(j_loop, k_loop, &
                                        & l_loop) = flux_src_vf(eqn_idx%mom%beg + i_dim - 1)%sf(j_loop, k_loop, &
                                        & l_loop) - current_tau_bulk(norm_dir, i_dim)

                            flux_src_vf(eqn_idx%E)%sf(j_loop, k_loop, l_loop) = flux_src_vf(eqn_idx%E)%sf(j_loop, k_loop, &
                                        & l_loop) - vel_src_at_interface(i_dim)*current_tau_bulk(norm_dir, i_dim)
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_cartesian_viscous_source_flux

    !> Compute shear stress tensor components
    subroutine s_calculate_shear_stress_tensor(vel_grad_avg, Re_shear, divergence_v, tau_shear_out)

        $:GPU_ROUTINE(parallelism='[seq]')

        ! Arguments
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3, 3), intent(in)  :: vel_grad_avg
            real(wp), dimension(3, 3), intent(out) :: tau_shear_out
        #:else
            real(wp), dimension(num_dims, num_dims), intent(in)  :: vel_grad_avg
            real(wp), dimension(num_dims, num_dims), intent(out) :: tau_shear_out
        #:endif
        real(wp), intent(in) :: Re_shear
        real(wp), intent(in) :: divergence_v

        ! Local variables
        integer :: i_dim  !< Loop iterator for face normal.
        integer :: j_dim  !< Loop iterator for force component direction.
        tau_shear_out = 0.0_wp

        do i_dim = 1, num_dims
            do j_dim = 1, num_dims
                tau_shear_out(i_dim, j_dim) = (vel_grad_avg(j_dim, i_dim) + vel_grad_avg(i_dim, j_dim))/Re_shear
                if (i_dim == j_dim) then
                    tau_shear_out(i_dim, j_dim) = tau_shear_out(i_dim, j_dim) - (2.0_wp/3.0_wp)*divergence_v/Re_shear
                end if
            end do
        end do

    end subroutine s_calculate_shear_stress_tensor

    !> Compute bulk stress tensor components (diagonal only)
    subroutine s_calculate_bulk_stress_tensor(Re_bulk, divergence_v, tau_bulk_out)

        $:GPU_ROUTINE(parallelism='[seq]')

        ! Arguments
        real(wp), intent(in) :: Re_bulk
        real(wp), intent(in) :: divergence_v
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3, 3), intent(out) :: tau_bulk_out
        #:else
            real(wp), dimension(num_dims, num_dims), intent(out) :: tau_bulk_out
        #:endif

        ! Local variables
        integer :: i_dim  !< Loop iterator for diagonal components.
        tau_bulk_out = 0.0_wp

        do i_dim = 1, num_dims
            tau_bulk_out(i_dim, i_dim) = divergence_v/Re_bulk
        end do

    end subroutine s_calculate_bulk_stress_tensor

    !> Reshape and copy the Riemann-solver flux buffers back to the physical-space output arrays for the selected sweep direction,
    !! finalizing the Riemann solve. Two variants are emitted from one template so the shared unpermute logic cannot drift apart:
    !! the plain routine also copies the advection flux_src set and the grid_geometry==3 z-sweep geometric source flux, while the
    !! _hatR variant unpermutes the hat_R-anchored flux_hatR_rs* set of the fused dual-pass HLLD solve (called between the two RHS
    !! assemblies) and is a strict subset: flux_src is anchor-independent (already finalized with the hat_L set) and its geometric
    !! source flux only exists for the axisymmetric y-sweep.
    #:for SUFFIX in ['', '_hatR']
        #:set INFIX = 'hatR_' if SUFFIX else ''
        #:if SUFFIX == ''
            subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

                type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf

            #:else
                subroutine s_finalize_riemann_solver_hatR(flux_vf, flux_gsrc_vf, norm_dir)

                    type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_gsrc_vf

                #:endif
                integer, intent(in) :: norm_dir
                integer             :: i, j, k, l  !< Generic loop iterators
                ! Reshaping Outputted Data in y-direction

                if (norm_dir == 2) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = 1, sys_size
                        do l = is3%beg, is3%end
                            do j = is1%beg, is1%end
                                do k = is2%beg, is2%end
                                    flux_vf(i)%sf(k, j, l) = flux_${INFIX}$rsx_vf(k, j, l, i)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    if (cyl_coord) then
                        $:GPU_PARALLEL_LOOP(collapse=4)
                        do i = 1, sys_size
                            do l = is3%beg, is3%end
                                do j = is1%beg, is1%end
                                    do k = is2%beg, is2%end
                                        flux_gsrc_vf(i)%sf(k, j, l) = flux_gsrc_${INFIX}$rsx_vf(k, j, l, i)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if

                    #:if SUFFIX == ''
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do l = is3%beg, is3%end
                            do j = is1%beg, is1%end
                                do k = is2%beg, is2%end
                                    flux_src_vf(eqn_idx%adv%beg)%sf(k, j, l) = flux_src_rsx_vf(k, j, l, eqn_idx%adv%beg)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        ! Copy the per-fluid flux_src entries when they are structurally present. HLLD writes zeros here; these
                        ! entries
                        ! are kept for consistency.
                        if (adv_src_mode == adv_src_mode_alpha_iface .or. adv_src_mode == adv_src_mode_none) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                                do l = is3%beg, is3%end
                                    do j = is1%beg, is1%end
                                        do k = is2%beg, is2%end
                                            flux_src_vf(i)%sf(k, j, l) = flux_src_rsx_vf(k, j, l, i)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    #:endif
                    ! Reshaping Outputted Data in z-direction
                else if (norm_dir == 3) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = 1, sys_size
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                do l = is3%beg, is3%end
                                    flux_vf(i)%sf(l, k, j) = flux_${INFIX}$rsx_vf(l, k, j, i)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    #:if SUFFIX == ''
                        if (grid_geometry == 3) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = 1, sys_size
                                do j = is1%beg, is1%end
                                    do k = is2%beg, is2%end
                                        do l = is3%beg, is3%end
                                            flux_gsrc_vf(i)%sf(l, k, j) = flux_gsrc_rsx_vf(l, k, j, i)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if

                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                do l = is3%beg, is3%end
                                    flux_src_vf(eqn_idx%adv%beg)%sf(l, k, j) = flux_src_rsx_vf(l, k, j, eqn_idx%adv%beg)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        ! Copy the per-fluid flux_src entries when they are structurally present. HLLD writes zeros here; these
                        ! entries
                        ! are kept for consistency.
                        if (adv_src_mode == adv_src_mode_alpha_iface .or. adv_src_mode == adv_src_mode_none) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                                do j = is1%beg, is1%end
                                    do k = is2%beg, is2%end
                                        do l = is3%beg, is3%end
                                            flux_src_vf(i)%sf(l, k, j) = flux_src_rsx_vf(l, k, j, i)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    #:endif
                else if (norm_dir == 1) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = 1, sys_size
                        do l = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                do j = is1%beg, is1%end
                                    flux_vf(i)%sf(j, k, l) = flux_${INFIX}$rsx_vf(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    #:if SUFFIX == ''
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do l = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                do j = is1%beg, is1%end
                                    flux_src_vf(eqn_idx%adv%beg)%sf(j, k, l) = flux_src_rsx_vf(j, k, l, eqn_idx%adv%beg)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        ! Copy the per-fluid flux_src entries when they are structurally present. HLLD writes zeros here; these
                        ! entries
                        ! are kept for consistency.
                        if (adv_src_mode == adv_src_mode_alpha_iface .or. adv_src_mode == adv_src_mode_none) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                                do l = is3%beg, is3%end
                                    do k = is2%beg, is2%end
                                        do j = is1%beg, is1%end
                                            flux_src_vf(i)%sf(j, k, l) = flux_src_rsx_vf(j, k, l, i)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    #:endif
                end if

            end subroutine s_finalize_riemann_solver${SUFFIX}$
        #:endfor
    end module m_riemann_state
