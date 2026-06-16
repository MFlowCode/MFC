!>
!! @file
!! @brief Contains module m_surface_tension

#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'inline_capillary.fpp'

!> @brief Computes capillary source fluxes and color-function gradients for the diffuse-interface surface tension model
module m_surface_tension

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_weno
    use m_muscl
    use m_helper
    use m_boundary_common

    implicit none

    private; public :: s_initialize_surface_tension_module, s_compute_capillary_source_flux, s_get_capillary, &
        & s_finalize_surface_tension_module

    !> @name color function gradient components and magnitude
    !> @{
    type(scalar_field), allocatable, dimension(:) :: c_divs
    !> @}
    $:GPU_DECLARE(create='[c_divs]')

    !> @name cell boundary reconstructed gradient components and magnitude
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: gL_x, gR_x
    !> @}
    $:GPU_DECLARE(create='[gL_x, gR_x]')

    type(int_bounds_info) :: is1, is2, is3, iv
    $:GPU_DECLARE(create='[is1, is2, is3, iv]')

contains

    !> Allocate and initialize surface tension module arrays
    impure subroutine s_initialize_surface_tension_module

        integer :: j

        @:ALLOCATE(c_divs(1:num_dims + 1))

        do j = 1, num_dims + 1
            @:ALLOCATE(c_divs(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(c_divs(j))
        end do

        @:ALLOCATE(gL_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))
        @:ALLOCATE(gR_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))

    end subroutine s_initialize_surface_tension_module

    !> Compute the capillary source flux from reconstructed color-gradient fields
    subroutine s_compute_capillary_source_flux(vSrc_rsx_vf, flux_src_vf, id, isx, isy, isz)

        real(wp), dimension(-1:,-1:,-1:,1:), intent(in)        :: vSrc_rsx_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        integer, intent(in)                                    :: id
        type(int_bounds_info), intent(in)                      :: isx, isy, isz

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3, 3) :: Omega
        #:else
            real(wp), dimension(num_dims, num_dims) :: Omega
        #:endif
        real(wp) :: w1L, w1R, w2L, w2R, w3L, w3R, w1, w2, w3
        real(wp) :: normWL, normWR, normW
        integer  :: j, k, l, i

        if (id == 1) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[Omega, w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW]')
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end
                        w1L = gL_x(j, k, l, 1)
                        w2L = gL_x(j, k, l, 2)
                        w3L = 0._wp
                        if (p > 0) w3L = gL_x(j, k, l, 3)

                        w1R = gR_x(j + 1, k, l, 1)
                        w2R = gR_x(j + 1, k, l, 2)
                        w3R = 0._wp
                        if (p > 0) w3R = gR_x(j + 1, k, l, 3)

                        normWL = gL_x(j, k, l, num_dims + 1)
                        normWR = gR_x(j + 1, k, l, num_dims + 1)

                        w1 = (w1L + w1R)/2._wp
                        w2 = (w2L + w2R)/2._wp
                        w3 = (w3L + w3R)/2._wp
                        normW = (normWL + normWR)/2._wp

                        if (normW > capillary_cutoff) then
                            @:compute_capillary_stress_tensor()

                            do i = 1, num_dims
                                flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, &
                                            & l) + Omega(1, i)

                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) + Omega(1, &
                                            & i)*vSrc_rsx_vf(j, k, l, i)
                            end do

                            ! Continuum surface force capillary stress, Schmidmayer et al. JCP (2017)
                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                        & l) + sigma*c_divs(num_dims + 1)%sf(j, k, l)*vSrc_rsx_vf(j, k, l, 1)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (id == 2) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                $:GPU_PARALLEL_LOOP(collapse=3, private='[Omega, w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW]')
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end
                            w1L = gL_x(j, k, l, 1)
                            w2L = gL_x(j, k, l, 2)
                            w3L = 0._wp
                            if (p > 0) w3L = gL_x(j, k, l, 3)

                            w1R = gR_x(j, k + 1, l, 1)
                            w2R = gR_x(j, k + 1, l, 2)
                            w3R = 0._wp
                            if (p > 0) w3R = gR_x(j, k + 1, l, 3)

                            normWL = gL_x(j, k, l, num_dims + 1)
                            normWR = gR_x(j, k + 1, l, num_dims + 1)

                            w1 = (w1L + w1R)/2._wp
                            w2 = (w2L + w2R)/2._wp
                            w3 = (w3L + w3R)/2._wp
                            normW = (normWL + normWR)/2._wp

                            if (normW > capillary_cutoff) then
                                @:compute_capillary_stress_tensor()

                                do i = 1, num_dims
                                    flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, &
                                                & k, l) + Omega(2, i)

                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) + Omega(2, &
                                                & i)*vSrc_rsx_vf(j, k, l, i)
                                end do

                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                            & l) + sigma*c_divs(num_dims + 1)%sf(j, k, l)*vSrc_rsx_vf(j, k, l, 2)
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:endif
        else if (id == 3) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                $:GPU_PARALLEL_LOOP(collapse=3, private='[Omega, w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW]')
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end
                            w1L = gL_x(j, k, l, 1)
                            w2L = gL_x(j, k, l, 2)
                            w3L = 0._wp
                            if (p > 0) w3L = gL_x(j, k, l, 3)

                            w1R = gR_x(j, k, l + 1, 1)
                            w2R = gR_x(j, k, l + 1, 2)
                            w3R = 0._wp
                            if (p > 0) w3R = gR_x(j, k, l + 1, 3)

                            normWL = gL_x(j, k, l, num_dims + 1)
                            normWR = gR_x(j, k, l + 1, num_dims + 1)

                            w1 = (w1L + w1R)/2._wp
                            w2 = (w2L + w2R)/2._wp
                            w3 = (w3L + w3R)/2._wp
                            normW = (normWL + normWR)/2._wp

                            if (normW > capillary_cutoff) then
                                @:compute_capillary_stress_tensor()

                                do i = 1, num_dims
                                    flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, &
                                                & k, l) + Omega(3, i)

                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) + Omega(3, &
                                                & i)*vSrc_rsx_vf(j, k, l, i)
                                end do

                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                            & l) + sigma*c_divs(num_dims + 1)%sf(j, k, l)*vSrc_rsx_vf(j, k, l, 3)
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:endif
        end if

    end subroutine s_compute_capillary_source_flux

    !> Compute color-function gradients and reconstruct them at cell boundaries
    impure subroutine s_get_capillary(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(in)        :: q_prim_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(int_bounds_info)                                      :: isx, isy, isz
        integer                                                    :: j, k, l, i

        isx%beg = -1; isy%beg = 0; isz%beg = 0

        if (m > 0) isy%beg = -1; if (p > 0) isz%beg = -1

        isx%end = m; isy%end = n; isz%end = p

        ! compute gradient components
        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs(1)%sf(j, k, l) = 1._wp/(x_cc(j + 1) - x_cc(j - 1))*(q_prim_vf(eqn_idx%c)%sf(j + 1, k, &
                           & l) - q_prim_vf(eqn_idx%c)%sf(j - 1, k, l))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs(2)%sf(j, k, l) = 1._wp/(y_cc(k + 1) - y_cc(k - 1))*(q_prim_vf(eqn_idx%c)%sf(j, k + 1, &
                           & l) - q_prim_vf(eqn_idx%c)%sf(j, k - 1, l))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (p > 0) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        c_divs(3)%sf(j, k, l) = 1._wp/(z_cc(l + 1) - z_cc(l - 1))*(q_prim_vf(eqn_idx%c)%sf(j, k, &
                               & l + 1) - q_prim_vf(eqn_idx%c)%sf(j, k, l - 1))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs(num_dims + 1)%sf(j, k, l) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_dims
                        c_divs(num_dims + 1)%sf(j, k, l) = c_divs(num_dims + 1)%sf(j, k, l) + c_divs(i)%sf(j, k, l)**2._wp
                    end do

                    c_divs(num_dims + 1)%sf(j, k, l) = sqrt(real(c_divs(num_dims + 1)%sf(j, k, l), kind=wp))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_populate_capillary_buffers(c_divs, bc_type, bc_xyz_info(bc_x, bc_y, bc_z))

        iv%beg = 1; iv%end = num_dims + 1

        ! reconstruct gradient components at cell boundaries
        call s_reconstruct_cell_boundary_values_capillary(c_divs, gL_x, gR_x, i)

    end subroutine s_get_capillary

    !> Reconstruct left and right cell-boundary values of capillary variables
    subroutine s_reconstruct_cell_boundary_values_capillary(v_vf, vL_x, vR_x, norm_dir)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in)                                  :: v_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,iv%beg:), intent(out) :: vL_x
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,iv%beg:), intent(out) :: vR_x
        integer, intent(in)                                                                       :: norm_dir
        integer                                                                                   :: i, j, k, l

        $:GPU_UPDATE(device='[iv]')

        $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
        do i = iv%beg, iv%end
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        vL_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                        vR_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_reconstruct_cell_boundary_values_capillary

    !> Finalize the surface tension module
    impure subroutine s_finalize_surface_tension_module

        integer :: j

        do j = 1, num_dims
            @:DEALLOCATE(c_divs(j)%sf)
        end do

        @:DEALLOCATE(c_divs)

        @:DEALLOCATE(gL_x, gR_x)

    end subroutine s_finalize_surface_tension_module

end module m_surface_tension
