#:include 'macros.fpp'
#:include 'inline_capillary.fpp'

!> @brief This module is used to compute source terms for surface tension model
module m_surface_tension

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion

    use m_weno

    use m_muscl                !< Monotonic Upstream-centered (MUSCL)
                               !! schemes for conservation laws

    use m_helper

    use m_boundary_common

    implicit none

    private; public :: s_initialize_surface_tension_module, &
 s_compute_capillary_source_flux, &
 s_get_capillary, &
 s_finalize_surface_tension_module

    !> @name color function gradient components and magnitude
    !> @{
    type(scalar_field), allocatable, dimension(:) :: c_divs
    !> @)
    $:GPU_DECLARE(create='[c_divs]')

    !> @name cell boundary reconstructed gradient components and magnitude
    !> @{
    real(wp), allocatable, dimension(:, :, :, :) :: gL_x, gR_x, gL_y, gR_y, gL_z, gR_z
    !> @}
    $:GPU_DECLARE(create='[gL_x,gR_x,gL_y,gR_y,gL_z,gR_z]')

    type(int_bounds_info) :: is1, is2, is3, iv
    $:GPU_DECLARE(create='[is1,is2,is3,iv]')

contains

    impure subroutine s_initialize_surface_tension_module

        integer :: j

        @:ALLOCATE(c_divs(1:num_dims + 1))

        do j = 1, num_dims + 1
            @:ALLOCATE(c_divs(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(c_divs(j))
        end do

        @:ALLOCATE(gL_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))
        @:ALLOCATE(gR_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))

        @:ALLOCATE(gL_y(idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))
        @:ALLOCATE(gR_y(idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))

        if (p > 0) then
            @:ALLOCATE(gL_z(idwbuff(3)%beg:idwbuff(3)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, num_dims + 1))
            @:ALLOCATE(gR_z(idwbuff(3)%beg:idwbuff(3)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, num_dims + 1))
        end if
    end subroutine s_initialize_surface_tension_module

    pure subroutine s_compute_capillary_source_flux( &
        vSrc_rsx_vf, vSrc_rsy_vf, vSrc_rsz_vf, &
        flux_src_vf, &
        id, isx, isy, isz)

        real(wp), dimension(-1:, 0:, 0:, 1:), intent(in) :: vSrc_rsx_vf
        real(wp), dimension(-1:, 0:, 0:, 1:), intent(in) :: vSrc_rsy_vf
        real(wp), dimension(-1:, 0:, 0:, 1:), intent(in) :: vSrc_rsz_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_src_vf
        integer, intent(in) :: id
        type(int_bounds_info), intent(in) :: isx, isy, isz

        real(wp), dimension(num_dims, num_dims) :: Omega
        real(wp) :: w1L, w1R, w2L, w2R, w3L, w3R, w1, w2, w3
        real(wp) :: normWL, normWR, normW
        integer :: j, k, l, i

        if (id == 1) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[Omega, w1L, w2L, w3L, &
                & w1R, w2R, w3R, w1, w2, w3, normWL, &
                & normWR, normW]')
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

                                flux_src_vf(momxb + i - 1)%sf(j, k, l) = &
                                    flux_src_vf(momxb + i - 1)%sf(j, k, l) + Omega(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                                 Omega(1, i)*vSrc_rsx_vf(j, k, l, i)

                            end do

                            flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                             sigma*c_divs(num_dims + 1)%sf(j, k, l)*vSrc_rsx_vf(j, k, l, 1)
                        end if
                    end do
                end do
            end do

        elseif (id == 2) then

            $:GPU_PARALLEL_LOOP(collapse=3, private='[Omega, w1L, w2L, w3L, &
                & w1R, w2R, w3R, w1, w2, w3, normWL, normWR, &
                & normW]')
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end

                        w1L = gL_y(k, j, l, 1)
                        w2L = gL_y(k, j, l, 2)
                        w3L = 0._wp
                        if (p > 0) w3L = gL_y(k, j, l, 3)

                        w1R = gR_y(k + 1, j, l, 1)
                        w2R = gR_y(k + 1, j, l, 2)
                        w3R = 0._wp
                        if (p > 0) w3R = gR_y(k + 1, j, l, 3)

                        normWL = gL_y(k, j, l, num_dims + 1)
                        normWR = gR_y(k + 1, j, l, num_dims + 1)

                        w1 = (w1L + w1R)/2._wp
                        w2 = (w2L + w2R)/2._wp
                        w3 = (w3L + w3R)/2._wp
                        normW = (normWL + normWR)/2._wp

                        if (normW > capillary_cutoff) then
                            @:compute_capillary_stress_tensor()

                            do i = 1, num_dims

                                flux_src_vf(momxb + i - 1)%sf(j, k, l) = &
                                    flux_src_vf(momxb + i - 1)%sf(j, k, l) + Omega(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                                 Omega(2, i)*vSrc_rsy_vf(k, j, l, i)

                            end do

                            flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                             sigma*c_divs(num_dims + 1)%sf(j, k, l)*vSrc_rsy_vf(k, j, l, 2)
                        end if
                    end do
                end do
            end do

        elseif (id == 3) then

            $:GPU_PARALLEL_LOOP(collapse=3, private='[Omega, w1L, w2L, w3L, &
                & w1R, w2R, w3R, w1, w2, w3, normWL, normWR, &
                & normW]')
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end

                        w1L = gL_z(l, k, j, 1)
                        w2L = gL_z(l, k, j, 2)
                        w3L = 0._wp
                        if (p > 0) w3L = gL_z(l, k, j, 3)

                        w1R = gR_z(l + 1, k, j, 1)
                        w2R = gR_z(l + 1, k, j, 2)
                        w3R = 0._wp
                        if (p > 0) w3R = gR_z(l + 1, k, j, 3)

                        normWL = gL_z(l, k, j, num_dims + 1)
                        normWR = gR_z(l + 1, k, j, num_dims + 1)

                        w1 = (w1L + w1R)/2._wp
                        w2 = (w2L + w2R)/2._wp
                        w3 = (w3L + w3R)/2._wp
                        normW = (normWL + normWR)/2._wp

                        if (normW > capillary_cutoff) then
                            @:compute_capillary_stress_tensor()

                            do i = 1, num_dims

                                flux_src_vf(momxb + i - 1)%sf(j, k, l) = &
                                    flux_src_vf(momxb + i - 1)%sf(j, k, l) + Omega(3, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                                 Omega(3, i)*vSrc_rsz_vf(l, k, j, i)

                            end do

                            flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                             sigma*c_divs(num_dims + 1)%sf(j, k, l)*vSrc_rsz_vf(l, k, j, 3)
                        end if
                    end do
                end do
            end do

        end if

    end subroutine s_compute_capillary_source_flux

    impure subroutine s_get_capillary(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        type(int_bounds_info) :: isx, isy, isz
        integer :: j, k, l, i

        isx%beg = -1; isy%beg = 0; isz%beg = 0

        if (m > 0) isy%beg = -1; if (p > 0) isz%beg = -1

        isx%end = m; isy%end = n; isz%end = p

        ! compute gradient components
        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs(1)%sf(j, k, l) = 1._wp/(x_cc(j + 1) - x_cc(j - 1))* &
                                            (q_prim_vf(c_idx)%sf(j + 1, k, l) - q_prim_vf(c_idx)%sf(j - 1, k, l))
                end do
            end do
        end do

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs(2)%sf(j, k, l) = 1._wp/(y_cc(k + 1) - y_cc(k - 1))* &
                                            (q_prim_vf(c_idx)%sf(j, k + 1, l) - q_prim_vf(c_idx)%sf(j, k - 1, l))
                end do
            end do
        end do

        if (p > 0) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        c_divs(3)%sf(j, k, l) = 1._wp/(z_cc(l + 1) - z_cc(l - 1))* &
                                                (q_prim_vf(c_idx)%sf(j, k, l + 1) - q_prim_vf(c_idx)%sf(j, k, l - 1))
                    end do
                end do
            end do
        end if

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs(num_dims + 1)%sf(j, k, l) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_dims
                        c_divs(num_dims + 1)%sf(j, k, l) = &
                            c_divs(num_dims + 1)%sf(j, k, l) + &
                            c_divs(i)%sf(j, k, l)**2._wp
                    end do
                    c_divs(num_dims + 1)%sf(j, k, l) = &
                        sqrt(c_divs(num_dims + 1)%sf(j, k, l))
                end do
            end do
        end do

        call s_populate_capillary_buffers(c_divs, bc_type)

        iv%beg = 1; iv%end = num_dims + 1

        ! reconstruct gradient components at cell boundaries
        do i = 1, num_dims
            call s_reconstruct_cell_boundary_values_capillary(c_divs, gL_x, gL_y, gL_z, gR_x, gR_y, gR_z, i)
        end do

    end subroutine s_get_capillary

    subroutine s_reconstruct_cell_boundary_values_capillary(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, &
                                                            norm_dir)
        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, iv%beg:), intent(out) :: vL_x, vL_y, vL_z
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, iv%beg:), intent(out) :: vR_x, vR_y, vR_z
        integer, intent(in) :: norm_dir

        integer :: recon_dir !< Coordinate direction of the reconstruction

        integer :: i, j, k, l

        #:for SCHEME, TYPE in [('weno', 'WENO_TYPE'),('muscl', 'MUSCL_TYPE')]
            if (recon_type == ${TYPE}$) then
                ! Reconstruction in s1-direction

                if (norm_dir == 1) then
                    is1 = idwbuff(1); is2 = idwbuff(2); is3 = idwbuff(3)
                    recon_dir = 1; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn

                elseif (norm_dir == 2) then
                    is1 = idwbuff(2); is2 = idwbuff(1); is3 = idwbuff(3)
                    recon_dir = 2; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn

                else
                    is1 = idwbuff(3); is2 = idwbuff(2); is3 = idwbuff(1)
                    recon_dir = 3; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn

                end if

                $:GPU_UPDATE(device='[is1,is2,is3,iv]')

                if (recon_dir == 1) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = iv%beg, iv%end
                        do l = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                do j = is1%beg, is1%end
                                    vL_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                                    vR_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                                end do
                            end do
                        end do
                    end do
                else if (recon_dir == 2) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = iv%beg, iv%end
                        do l = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                do j = is1%beg, is1%end
                                    vL_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                                    vR_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                                end do
                            end do
                        end do
                    end do
                else if (recon_dir == 3) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = iv%beg, iv%end
                        do l = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                do j = is1%beg, is1%end
                                    vL_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                                    vR_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                                end do
                            end do
                        end do
                    end do
                end if
            end if
        #:endfor

    end subroutine s_reconstruct_cell_boundary_values_capillary

    impure subroutine s_finalize_surface_tension_module
        integer :: j

        do j = 1, num_dims
            @:DEALLOCATE(c_divs(j)%sf)
        end do

        @:DEALLOCATE(c_divs)

        @:DEALLOCATE(gL_x, gR_x)

        @:DEALLOCATE(gL_y, gR_y)
        if (p > 0) then
            @:DEALLOCATE(gL_z, gR_z)
        end if

    end subroutine s_finalize_surface_tension_module

end module m_surface_tension
