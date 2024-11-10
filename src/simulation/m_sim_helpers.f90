module m_sim_helpers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters

    use m_variables_conversion
    ! ==========================================================================

    implicit none

    private; public :: s_compute_enthalpy, &
 s_compute_stability_from_dt, &
 s_compute_dt_from_cfl

contains

    !> Computes enthalpy
        !! @param q_prim_vf cell centered primitive variables
        !! @param pres mixture pressure
        !! @param rho mixture density
        !! @param gamma mixture gamma
        !! @param pi_inf mixture pi_inf
        !! @param Re mixture reynolds number
        !! @param H mixture enthalpy
        !! @param alpha component alphas
        !! @param vel directional velocities
        !! @param vel_sum squard sum of velocity components
        !! @param j x index
        !! @param k y index
        !! @param l z index
    subroutine s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, j, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_compute_enthalpy
#else
        !$acc routine seq
#endif

        type(scalar_field), dimension(sys_size) :: q_prim_vf
        real(wp), dimension(num_fluids) :: alpha_rho
        real(wp), dimension(num_fluids) :: alpha
        real(wp), dimension(num_dims) :: vel
        real(wp) :: rho, gamma, pi_inf, qv, vel_sum, E, H, pres
        real(wp), dimension(2) :: Re
        integer :: i, j, k, l

        do i = 1, num_fluids
            alpha_rho(i) = q_prim_vf(i)%sf(j, k, l)
            alpha(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
        end do

        if (bubbles) then
            call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, j, k, l)
        else
            call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, j, k, l)
        end if

        do i = 1, num_dims
            vel(i) = q_prim_vf(contxe + i)%sf(j, k, l)
        end do

        vel_sum = 0._wp
        do i = 1, num_dims
            vel_sum = vel_sum + vel(i)**2._wp
        end do

        pres = q_prim_vf(E_idx)%sf(j, k, l)

        E = gamma*pres + pi_inf + 5e-1_wp*rho*vel_sum + qv

        H = (E + pres)/rho

    end subroutine s_compute_enthalpy

    !> Computes stability criterion for a specified dt
        !! @param vel directional velocities
        !! @param c mixture speed of sound
        !! @param Re_l mixture Reynolds number
        !! @param j x index
        !! @param k y index
        !! @param l z index
        !! @param icfl_sf cell centered inviscid cfl number
        !! @param vcfl_sf (optional) cell centered viscous cfl number
        !! @param Rc_sf (optional) cell centered Rc
    subroutine s_compute_stability_from_dt(vel, c, rho, Re_l, j, k, l, icfl_sf, vcfl_sf, Rc_sf)
        !$acc routine seq
        real(wp), dimension(num_dims) :: vel
        real(wp) :: c, rho
        real(wp), dimension(0:m, 0:n, 0:p) :: icfl_sf
        real(wp), dimension(0:m, 0:n, 0:p), optional :: vcfl_sf, Rc_sf
        real(wp) :: fltr_dtheta   !<
             !! Modified dtheta accounting for Fourier filtering in azimuthal direction.
        integer :: j, k, l
        integer :: Nfq
        real(wp), dimension(2) :: Re_l

        if (grid_geometry == 3) then
            if (k == 0) then
                fltr_dtheta = 2._wp*pi*y_cb(0)/3._wp
            elseif (k <= fourier_rings) then
                Nfq = min(floor(2._wp*real(k, wp)*pi), (p + 1)/2 + 1)
                fltr_dtheta = 2._wp*pi*y_cb(k - 1)/real(Nfq, wp)
            else
                fltr_dtheta = y_cb(k - 1)*dz(l)
            end if
        end if

        if (p > 0) then
            !3D
            if (grid_geometry == 3) then
                icfl_sf(j, k, l) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                          dy(k)/(abs(vel(2)) + c), &
                                          fltr_dtheta/(abs(vel(3)) + c))
            else
                icfl_sf(j, k, l) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                          dy(k)/(abs(vel(2)) + c), &
                                          dz(l)/(abs(vel(3)) + c))
            end if

            if (viscous) then

                if (grid_geometry == 3) then
                    vcfl_sf(j, k, l) = maxval(dt/Re_l/rho) &
                                       /min(dx(j), dy(k), fltr_dtheta)**2._wp

                    Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), &
                                         dy(k)*(abs(vel(2)) + c), &
                                         fltr_dtheta*(abs(vel(3)) + c)) &
                                     /maxval(1._wp/Re_l)
                else
                    vcfl_sf(j, k, l) = maxval(dt/Re_l/rho) &
                                       /min(dx(j), dy(k), dz(l))**2._wp

                    Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), &
                                         dy(k)*(abs(vel(2)) + c), &
                                         dz(l)*(abs(vel(3)) + c)) &
                                     /maxval(1._wp/Re_l)
                end if

            end if

        elseif (n > 0) then
            !2D
            icfl_sf(j, k, l) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                      dy(k)/(abs(vel(2)) + c))

            if (viscous) then

                vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/min(dx(j), dy(k))**2._wp

                Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), &
                                     dy(k)*(abs(vel(2)) + c)) &
                                 /maxval(1._wp/Re_l)

            end if

        else
            !1D
            icfl_sf(j, k, l) = (dt/dx(j))*(abs(vel(1)) + c)

            if (viscous) then

                vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/dx(j)**2._wp

                Rc_sf(j, k, l) = dx(j)*(abs(vel(1)) + c)/maxval(1._wp/Re_l)

            end if

        end if

    end subroutine s_compute_stability_from_dt

    !> Computes dt for a specified CFL number
        !! @param vel directional velocities
        !! @param max_dt cell centered maximum dt
        !! @param rho cell centered density
        !! @param Re_l cell centered Reynolds number
        !! @param j x coordinate
        !! @param k y coordinate
        !! @param l z coordinate
    subroutine s_compute_dt_from_cfl(vel, c, max_dt, rho, Re_l, j, k, l)
        !$acc routine seq
        real(wp), dimension(num_dims) :: vel
        real(wp) :: c, icfl_dt, vcfl_dt, rho
        real(wp), dimension(0:m, 0:n, 0:p) :: max_dt
        real(wp) :: fltr_dtheta   !<
             !! Modified dtheta accounting for Fourier filtering in azimuthal direction.
        integer :: j, k, l
        integer :: Nfq
        real(wp), dimension(2) :: Re_l

        if (grid_geometry == 3) then
            if (k == 0) then
                fltr_dtheta = 2._wp*pi*y_cb(0)/3._wp
            elseif (k <= fourier_rings) then
                Nfq = min(floor(2._wp*real(k, wp)*pi), (p + 1)/2 + 1)
                fltr_dtheta = 2._wp*pi*y_cb(k - 1)/real(Nfq, wp)
            else
                fltr_dtheta = y_cb(k - 1)*dz(l)
            end if
        end if

        if (p > 0) then
            !3D
            if (grid_geometry == 3) then
                icfl_dt = cfl_target*min(dx(j)/(abs(vel(1)) + c), &
                                         dy(k)/(abs(vel(2)) + c), &
                                         fltr_dtheta/(abs(vel(3)) + c))
            else
                icfl_dt = cfl_target*min(dx(j)/(abs(vel(1)) + c), &
                                         dy(k)/(abs(vel(2)) + c), &
                                         dz(l)/(abs(vel(3)) + c))
            end if

            if (viscous) then
                if (grid_geometry == 3) then
                    vcfl_dt = cfl_target*(min(dx(j), dy(k), fltr_dtheta)**2._wp) &
                              /minval(1/(rho*Re_l))
                else
                    vcfl_dt = cfl_target*(min(dx(j), dy(k), dz(l))**2._wp) &
                              /minval(1/(rho*Re_l))
                end if
            end if

        elseif (n > 0) then
            !2D
            icfl_dt = cfl_target*min(dx(j)/(abs(vel(1)) + c), &
                                     dy(k)/(abs(vel(2)) + c))

            if (viscous) then
                vcfl_dt = cfl_target*(min(dx(j), dy(k))**2._wp)/maxval((1/Re_l)/rho)
            end if

        else
            !1D
            icfl_dt = cfl_target*(dx(j)/(abs(vel(1)) + c))

            if (viscous) then
                vcfl_dt = cfl_target*(dx(j)**2._wp)/minval(1/(rho*Re_l))
            end if

        end if

        if (any(re_size > 0)) then
            max_dt(j, k, l) = min(icfl_dt, vcfl_dt)
        else
            max_dt(j, k, l) = icfl_dt
        end if

    end subroutine s_compute_dt_from_cfl

end module m_sim_helpers
