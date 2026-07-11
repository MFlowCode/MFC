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
    use m_constants, only: riemann_solver_hll, riemann_solver_hlld, verysmall
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
            if (norm_dir == 1) then
                dir_idx_tau = (/1, 2, 4/)
            else if (norm_dir == 2) then
                dir_idx_tau = (/3, 2, 5/)
            else
                dir_idx_tau = (/6, 4, 5/)
            end if
        end if

        isx = ix; isy = iy; isz = iz
        ! for stuff in the same module
        $:GPU_UPDATE(device='[isx, isy, isz]')
        ! for stuff in different modules
        $:GPU_UPDATE(device='[dir_idx, dir_flg, dir_idx_tau]')

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

    !> Accumulate the mixture density, specific heat ratio function, liquid stiffness function, and internal energy reference of one
    !! Riemann state from its partial densities and volume fractions. The number of fluids is an explicit argument because the
    !! 5-equation bubble model accumulates over num_fluids - 1 fluids.
    subroutine s_accumulate_mixture_properties(nf, alpha_rho_K, alpha_K, rho_K, gamma_K, pi_inf_K, qv_K)

        $:GPU_ROUTINE(function_name='s_accumulate_mixture_properties', parallelism='[seq]', cray_inline=True)

        integer, intent(in)                 :: nf  !< Number of fluids to accumulate over
        real(wp), dimension(nf), intent(in) :: alpha_rho_K, alpha_K
        real(wp), intent(out)               :: rho_K, gamma_K, pi_inf_K, qv_K
        integer                             :: i   !< Loop iterator over fluids

        rho_K = 0._wp
        gamma_K = 0._wp
        pi_inf_K = 0._wp
        qv_K = 0._wp

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, nf
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
            qv_K = qv_K + alpha_rho_K(i)*qvs(i)
        end do

    end subroutine s_accumulate_mixture_properties

    !> Compute the shear and volume Reynolds numbers of one Riemann state by inverse-weighting the fluid Reynolds numbers with the
    !! volume fractions.
    subroutine s_compute_interface_reynolds(alpha_K, Re_K, Re_size_loc1, Re_size_loc2)

        $:GPU_ROUTINE(function_name='s_compute_interface_reynolds', parallelism='[seq]', cray_inline=True)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha_K
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha_K
        #:endif
        real(wp), dimension(2), intent(out) :: Re_K
        !> host copies of Re_size; amdflang reads the declare-target original stale cross-TU
        integer, intent(in) :: Re_size_loc1, Re_size_loc2
        integer             :: i, q  !< Loop iterators

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, 2
            Re_K(i) = dflt_real

            if (merge(Re_size_loc1, Re_size_loc2, i == 1) > 0) Re_K(i) = 0._wp

            $:GPU_LOOP(parallelism='[seq]')
            do q = 1, merge(Re_size_loc1, Re_size_loc2, i == 1)
                Re_K(i) = alpha_K(Re_idx(i, q))/Res_gs(i, q) + Re_K(i)
            end do

            Re_K(i) = 1._wp/max(Re_K(i), sgm_eps)
        end do

    end subroutine s_compute_interface_reynolds

    !> Accumulate the hypoelastic stress contribution to the energies of the left and right Riemann states: mix the shear modulus
    !! over the fluids, scale it by the continuum damage state when damage is modeled, and add the elastic energy of each stress
    !! component (doubled for the shear components) when both mixture moduli are non-negligible. The elastic shear stresses are
    !! loaded from the state buffers by the caller, which reuses them for the stress fluxes and elastic wave speeds. The G >
    !! verysmall gate is a deliberate maintainer ruling that replaces HLL's former hard-coded G > 1000 stability floor, retiring its
    !! "TODO take out if statement if stable without".
    subroutine s_compute_hypoelastic_interface_energy(nf, alpha_L, alpha_R, damage_L, damage_R, tau_e_L, tau_e_R, G_L, G_R, E_L, &
        & E_R)

        $:GPU_ROUTINE(function_name='s_compute_hypoelastic_interface_energy', parallelism='[seq]', cray_inline=True)

        integer, intent(in)                 :: nf                  !< Number of fluids to mix the shear modulus over
        real(wp), dimension(nf), intent(in) :: alpha_L, alpha_R    !< Left and right volume fractions
        real(wp), intent(in)                :: damage_L, damage_R  !< Continuum damage states (referenced only when cont_damage)
        real(wp), dimension(6), intent(in)  :: tau_e_L, tau_e_R    !< Left and right elastic shear stresses
        real(wp), intent(out)               :: G_L, G_R            !< Left and right mixture shear moduli
        real(wp), intent(inout)             :: E_L, E_R            !< Left and right state energies
        integer                             :: i                   !< Loop iterator
        real(wp)                            :: G_gate              !< Floor below which the elastic energy term is skipped

        G_L = 0._wp; G_R = 0._wp

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, nf
            G_L = G_L + alpha_L(i)*Gs_rs(i)
            G_R = G_R + alpha_R(i)*Gs_rs(i)
        end do

        if (cont_damage) then
            G_L = G_L*max((1._wp - damage_L), 0._wp)
            G_R = G_R*max((1._wp - damage_R), 0._wp)
        end if

        ! Under continuum damage a heavily-damaged interface can drive G -> 0 while the reconstructed stress does
        ! not relax with it, so tau^2/(4G) blows up. It stays finite (and negligible) on most backends but goes
        ! NaN under macOS gfortran's libm. Restore HLL's former stability floor for the damage case only; for
        ! undamaged states G >> this floor so master's verysmall gate is unchanged.
        G_gate = merge(1.e3_wp, verysmall, cont_damage)

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
            ! Elastic contribution to energy if G large enough
            if ((G_L > G_gate) .and. (G_R > G_gate)) then
                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                ! Double for shear stresses
                if (any(eqn_idx%stress%beg - 1 + i == shear_indices)) then
                    E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                    E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                end if
            end if
        end do

    end subroutine s_compute_hypoelastic_interface_energy

    !> Compute the advective part of the HLLC star-state momentum flux in the wave-normal direction (pressure excluded), used to
    !! assemble the geometrical source flux of the cylindrical and azimuthal sweeps.
    function f_compute_hllc_star_momentum_flux(rho_L, rho_R, vel_L_norm, vel_R_norm, s_M, s_P, s_S, xi_L, xi_R, xi_M, xi_P, &
        & dir_flg_norm) result(flux_mom)

        $:GPU_ROUTINE(function_name='f_compute_hllc_star_momentum_flux', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: rho_L, rho_R            !< Left and right densities
        real(wp), intent(in) :: vel_L_norm, vel_R_norm  !< Left and right wave-normal velocities
        real(wp), intent(in) :: s_M, s_P, s_S           !< Clamped left/right and contact wave speeds
        real(wp), intent(in) :: xi_L, xi_R, xi_M, xi_P  !< Star-state compression factors and upwind selectors
        real(wp), intent(in) :: dir_flg_norm            !< Direction flag of the wave-normal direction
        real(wp)             :: flux_mom

        flux_mom = xi_M*(rho_L*(vel_L_norm*vel_L_norm + s_M*(xi_L*(dir_flg_norm*s_S + (1._wp - dir_flg_norm)*vel_L_norm) &
                         & - vel_L_norm))) + xi_P*(rho_R*(vel_R_norm*vel_R_norm + s_P*(xi_R*(dir_flg_norm*s_S + (1._wp &
                         & - dir_flg_norm)*vel_R_norm) - vel_R_norm)))

    end function f_compute_hllc_star_momentum_flux

    !> Hybrid-Riemann smooth-face flux: overwrite the Riemann flux at a WENO-smooth face with either a central (hybrid_smooth_flux
    !! == 1) or a local Lax-Friedrichs / Rusanov (hybrid_smooth_flux == 2) flux, reusing the left/right states the calling solver
    !! already reconstructed. Shared by hll/hllc/lf so the smooth-flux path stays identical across solvers; the reshaped flux
    !! buffers are module-scope so only the (j,k,l) index and the small per-fluid/per-dim state arrays need to be passed.
    subroutine s_compute_hybrid_smooth_flux(j, k, l, alpha_rho_L, alpha_L, alpha_rho_R, alpha_R, vel_L, vel_R, c_L, c_R, rho_L, &
                                            & rho_R, pres_L, pres_R, E_L, E_R)

        $:GPU_ROUTINE(function_name='s_compute_hybrid_smooth_flux', parallelism='[seq]', cray_inline=True)

        integer, intent(in)                         :: j, k, l
        real(wp), dimension(num_fluids), intent(in) :: alpha_rho_L, alpha_L, alpha_rho_R, alpha_R
        real(wp), dimension(num_dims), intent(in)   :: vel_L, vel_R
        real(wp), intent(in)                        :: c_L, c_R, rho_L, rho_R, pres_L, pres_R, E_L, E_R
        real(wp)                                    :: lam, FL, FR, UL, UR
        integer                                     :: i

        ! Rusanov dissipation coefficient (dropped for the pure-central flux, hybrid_smooth_flux == 1)
        lam = max(abs(vel_L(dir_idx(1))) + c_L, abs(vel_R(dir_idx(1))) + c_R)

        ! Continuity
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, eqn_idx%cont%end
            FL = alpha_rho_L(i)*vel_L(dir_idx(1))
            FR = alpha_rho_R(i)*vel_R(dir_idx(1))
            flux_rsx_vf(j, k, l, i) = 0.5_wp*(FL + FR) - real(hybrid_smooth_flux - 1, &
                        & wp)*0.5_wp*lam*(alpha_rho_R(i) - alpha_rho_L(i))
        end do

        ! Momentum
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_dims
            FL = rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*pres_L
            FR = rho_R*vel_R(dir_idx(1))*vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*pres_R
            UL = rho_L*vel_L(dir_idx(i))
            UR = rho_R*vel_R(dir_idx(i))
            flux_rsx_vf(j, k, l, eqn_idx%cont%end + dir_idx(i)) = 0.5_wp*(FL + FR) - real(hybrid_smooth_flux - 1, &
                        & wp)*0.5_wp*lam*(UR - UL)
        end do

        ! Energy
        FL = vel_L(dir_idx(1))*(E_L + pres_L)
        FR = vel_R(dir_idx(1))*(E_R + pres_R)
        flux_rsx_vf(j, k, l, eqn_idx%E) = 0.5_wp*(FL + FR) - real(hybrid_smooth_flux - 1, wp)*0.5_wp*lam*(E_R - E_L)

        ! Volume fractions (advection)
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            FL = alpha_L(i)*vel_L(dir_idx(1))
            FR = alpha_R(i)*vel_R(dir_idx(1))
            flux_rsx_vf(j, k, l, eqn_idx%adv%beg + i - 1) = 0.5_wp*(FL + FR) - real(hybrid_smooth_flux - 1, &
                        & wp)*0.5_wp*lam*(alpha_R(i) - alpha_L(i))
        end do

        ! Internal energies (6-equation model only)
        if (model_eqns == model_eqns_6eq) then
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                UL = alpha_L(i)*(gammas(i)*pres_L + pi_infs(i)) + alpha_rho_L(i)*qvs(i)
                UR = alpha_R(i)*(gammas(i)*pres_R + pi_infs(i)) + alpha_rho_R(i)*qvs(i)
                FL = UL*vel_L(dir_idx(1))
                FR = UR*vel_R(dir_idx(1))
                flux_rsx_vf(j, k, l, eqn_idx%int_en%beg + i - 1) = 0.5_wp*(FL + FR) - real(hybrid_smooth_flux - 1, &
                            & wp)*0.5_wp*lam*(UR - UL)
            end do
        end if

        ! Non-conservative advection source: central interface velocity
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_dims
            vel_src_rsx_vf(j, k, l, dir_idx(i)) = 0.5_wp*(vel_L(dir_idx(i)) + vel_R(dir_idx(i)))
        end do
        flux_src_rsx_vf(j, k, l, eqn_idx%adv%beg) = vel_src_rsx_vf(j, k, l, dir_idx(1))

    end subroutine s_compute_hybrid_smooth_flux

    !> Deallocation and/or disassociation procedures that are needed to finalize the selected Riemann problem solver
    subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        integer, intent(in)                                    :: norm_dir
        integer                                                :: i, j, k, l  !< Generic loop iterators
        ! Reshaping Outputted Data in y-direction

        if (norm_dir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            flux_vf(i)%sf(k, j, l) = flux_rsx_vf(k, j, l, i)
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
                                flux_gsrc_vf(i)%sf(k, j, l) = flux_gsrc_rsx_vf(k, j, l, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3%beg, is3%end
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        flux_src_vf(eqn_idx%adv%beg)%sf(k, j, l) = flux_src_rsx_vf(k, j, l, eqn_idx%adv%beg)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) then
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
            ! Reshaping Outputted Data in z-direction
        else if (norm_dir == 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                            flux_vf(i)%sf(l, k, j) = flux_rsx_vf(l, k, j, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
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

            if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) then
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
        else if (norm_dir == 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            flux_vf(i)%sf(j, k, l) = flux_rsx_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = is1%beg, is1%end
                        flux_src_vf(eqn_idx%adv%beg)%sf(j, k, l) = flux_src_rsx_vf(j, k, l, eqn_idx%adv%beg)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) then
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
        end if

    end subroutine s_finalize_riemann_solver

end module m_riemann_state
