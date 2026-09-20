!>
!! @file
!! @brief Contains module m_conduction

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Fourier heat conduction, div(k grad T), as a face-centered source flux on the energy equation. The term has no
!! cross-derivatives, so the direction-split face difference below is exact for it.
module m_conduction

    use m_global_parameters

    implicit none

    private; public :: s_compute_conduction_source_flux, s_compute_conduction_axis_source

    type(int_bounds_info) :: isc1, isc2, isc3
    $:GPU_DECLARE(create='[isc1, isc2, isc3]')
    integer, dimension(3) :: offsets_c
    $:GPU_DECLARE(create='[offsets_c]')

contains

    !> Accumulate -k*dT/dx_idir into the energy source flux at each idir-normal face.
    subroutine s_compute_conduction_source_flux(idir, q_prim_qp, q_T_sf, flux_src_vf, irx, iry, irz)

        integer, intent(in)                                    :: idir
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_qp
        type(scalar_field), intent(in)                         :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(int_bounds_info), intent(in)                      :: irx, iry, irz
        real(wp)                                               :: k_face, dT_dxi, grid_spacing, alpha_face
        integer                                                :: x, y, z, i

        isc1 = irx; isc2 = iry; isc3 = irz
        offsets_c = 0
        offsets_c(idir) = 1

        $:GPU_UPDATE(device='[isc1, isc2, isc3, offsets_c]')

        $:GPU_PARALLEL_LOOP(collapse=3, private='[k_face, dT_dxi, grid_spacing, alpha_face, i]')
        do z = isc3%beg, isc3%end
            do y = isc2%beg, isc2%end
                do x = isc1%beg, isc1%end
                    select case (idir)
                    case (1)
                        grid_spacing = x_cc(x + 1) - x_cc(x)
                    case (2)
                        grid_spacing = y_cc(y + 1) - y_cc(y)
                    case (3)
                        grid_spacing = z_cc(z + 1) - z_cc(z)
                        ! 3D cylindrical z is the azimuth: dz is in radians. The r**2 carries both metric
                        ! factors of (1/r**2) d2T/dtheta2, since m_rhs divides this flux by dz alone.
                        if (grid_geometry == 3) grid_spacing = y_cc(y)**2*grid_spacing
                    end select

                    ! Volume-fraction-weighted face conductivity. Raw cell-centered alphas over- and
                    ! undershoot near interfaces, so clamp the face average as the viscous path does.
                    k_face = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_face = 0.5_wp*(q_prim_qp(eqn_idx%adv%beg + i - 1)%sf(x, y, &
                                             & z) + q_prim_qp(eqn_idx%adv%beg + i - 1)%sf(x + offsets_c(1), y + offsets_c(2), &
                                             & z + offsets_c(3)))
                        alpha_face = min(max(alpha_face, 0._wp), 1._wp)
                        k_face = k_face + alpha_face*fluid_k_therm(i)
                    end do

                    dT_dxi = (q_T_sf%sf(x + offsets_c(1), y + offsets_c(2), z + offsets_c(3)) - q_T_sf%sf(x, y, z))/grid_spacing

                    flux_src_vf(eqn_idx%E)%sf(x, y, z) = flux_src_vf(eqn_idx%E)%sf(x, y, z) - k_face*dT_dxi
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_conduction_source_flux

    !> Cell-centered -k*dT/dr for the axis cell of a cylindrical grid. The generic geometric source in m_rhs uses face values of
    !! flux_src(E) in this role; the axis cell has no face pair, so it reads tau_Re_vf(E) instead. Accumulates: call after
    !! s_compute_viscous_stress_cylindrical_boundary, which zeroes tau_Re_vf(mom%beg:E); without viscosity this routine does that
    !! zeroing itself.
    subroutine s_compute_conduction_axis_source(q_prim_vf, q_T_sf, tau_Re_vf, ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(in)      :: q_prim_vf
        type(scalar_field), intent(in)                           :: q_T_sf
        type(scalar_field), dimension(1:sys_size), intent(inout) :: tau_Re_vf
        type(int_bounds_info), intent(in)                        :: ix, iy, iz
        real(wp)                                                 :: k_cell, dT_dr, dr_l, dr_r, alpha_cell
        integer                                                  :: j, k, l, i

        isc1 = ix; isc2 = iy; isc3 = iz
        $:GPU_UPDATE(device='[isc1, isc2, isc3]')

        if (.not. viscous) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i]')
            do l = isc3%beg, isc3%end
                do k = isc2%beg, isc2%end
                    do j = isc1%beg, isc1%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%mom%beg, eqn_idx%E
                            tau_Re_vf(i)%sf(j, k, l) = 0._wp
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[k_cell, dT_dr, dr_l, dr_r, alpha_cell, i]')
        do l = isc3%beg, isc3%end
            do k = isc2%beg + 1, isc2%end - 1
                do j = isc1%beg, isc1%end
                    k_cell = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_cell = min(max(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), 0._wp), 1._wp)
                        k_cell = k_cell + alpha_cell*fluid_k_therm(i)
                    end do

                    ! Distance-weighted central difference. The axisymmetric grid halves the cell at
                    ! r = 0, so the plain form would be only first order over exactly this stencil.
                    dr_l = y_cc(k) - y_cc(k - 1)
                    dr_r = y_cc(k + 1) - y_cc(k)
                    dT_dr = ((q_T_sf%sf(j, k, l) - q_T_sf%sf(j, k - 1, l))*dr_r/dr_l + (q_T_sf%sf(j, k + 1, l) - q_T_sf%sf(j, k, &
                             & l))*dr_l/dr_r)/(dr_l + dr_r)

                    tau_Re_vf(eqn_idx%E)%sf(j, k, l) = tau_Re_vf(eqn_idx%E)%sf(j, k, l) - k_cell*dT_dr
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_conduction_axis_source

end module m_conduction
