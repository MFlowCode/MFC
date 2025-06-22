module m_sim_helpers

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters

    use m_variables_conversion

    implicit none

    private; public :: s_compute_enthalpy, &
 s_compute_stability_from_dt, &
 s_compute_dt_from_cfl

contains

    !> Computes CFL-related terms (handles both Fourier filtering and inviscid CFL)
        !! This function consolidates the repeated CFL calculation patterns found in both
        !! s_compute_stability_from_dt and s_compute_dt_from_cfl subroutines.
        !! It handles Fourier filtering for cylindrical coordinates and computes inviscid
        !! CFL terms for 1D, 2D, and 3D cases.
        !! @param vel directional velocities
        !! @param c mixture speed of sound
        !! @param j x coordinate index
        !! @param k y coordinate index  
        !! @param l z coordinate index
        !! @param for_dt_calc logical flag: true for dt calculation, false for stability calculation
        !! @return cfl_terms computed CFL terms (with appropriate scaling applied)
    pure function f_compute_cfl_terms(vel, c, j, k, l, for_dt_calc) result(cfl_terms)
        !$acc routine seq
        real(wp), dimension(num_vels), intent(in) :: vel
        real(wp), intent(in) :: c
        integer, intent(in) :: j, k, l
        logical, intent(in) :: for_dt_calc
        real(wp) :: cfl_terms
        real(wp) :: fltr_dtheta
        integer :: Nfq

        ! Compute filtered dtheta for cylindrical coordinates
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

        ! Compute CFL terms based on dimensionality
        if (p > 0) then
            !3D
            if (grid_geometry == 3) then
                cfl_terms = min(dx(j)/(abs(vel(1)) + c), &
                               dy(k)/(abs(vel(2)) + c), &
                               fltr_dtheta/(abs(vel(3)) + c))
            else
                cfl_terms = min(dx(j)/(abs(vel(1)) + c), &
                               dy(k)/(abs(vel(2)) + c), &
                               dz(l)/(abs(vel(3)) + c))
            end if
        elseif (n > 0) then
            !2D
            cfl_terms = min(dx(j)/(abs(vel(1)) + c), &
                           dy(k)/(abs(vel(2)) + c))
        else
            !1D - special handling for different calculation types
            if (for_dt_calc) then
                cfl_terms = dx(j)/(abs(vel(1)) + c)
            else
                cfl_terms = (1._wp/dx(j))*(abs(vel(1)) + c)
            end if
        end if

        ! Apply appropriate factor for stability vs dt calculation
        if (.not. for_dt_calc .and. (p > 0 .or. n > 0)) then
            cfl_terms = 1._wp/cfl_terms
        end if
    end function f_compute_cfl_terms

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
    pure subroutine s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, j, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_compute_enthalpy
#else
        !$acc routine seq
#endif

        type(scalar_field), intent(in), dimension(sys_size) :: q_prim_vf
        real(wp), intent(inout), dimension(num_fluids) :: alpha
        real(wp), intent(inout), dimension(num_vels) :: vel
        real(wp), intent(inout) :: rho, gamma, pi_inf, vel_sum, H, pres
        integer, intent(in) :: j, k, l
        real(wp), dimension(2), intent(inout) :: Re

        real(wp), dimension(num_fluids) :: alpha_rho, Gs
        real(wp) :: qv, E, G

        integer :: i

        !$acc loop seq
        do i = 1, num_fluids
            alpha_rho(i) = q_prim_vf(i)%sf(j, k, l)
            alpha(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
        end do

        if (elasticity) then
            call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, alpha, &
                                                            alpha_rho, Re, G, Gs)
        elseif (bubbles_euler) then
            call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re)
        else
            call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re)
        end if

        !$acc loop seq
        do i = 1, num_vels
            vel(i) = q_prim_vf(contxe + i)%sf(j, k, l)
        end do

        vel_sum = 0._wp
        !$acc loop seq
        do i = 1, num_vels
            vel_sum = vel_sum + vel(i)**2._wp
        end do

        pres = q_prim_vf(E_idx)%sf(j, k, l)

        E = gamma*pres + pi_inf + 5.e-1_wp*rho*vel_sum + qv

        ! ENERGY ADJUSTMENTS FOR HYPERELASTIC ENERGY
        if (hyperelasticity) then
            E = E + G*q_prim_vf(xiend + 1)%sf(j, k, l)
        end if

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
    pure subroutine s_compute_stability_from_dt(vel, c, rho, Re_l, j, k, l, icfl_sf, vcfl_sf, Rc_sf)
        !$acc routine seq
        real(wp), intent(in), dimension(num_vels) :: vel
        real(wp), intent(in) :: c, rho
        real(wp), dimension(0:m, 0:n, 0:p), intent(inout) :: icfl_sf
        real(wp), dimension(0:m, 0:n, 0:p), intent(inout), optional :: vcfl_sf, Rc_sf
        real(wp), dimension(2), intent(in) :: Re_l
        integer, intent(in) :: j, k, l

        ! Inviscid CFL calculation
        icfl_sf(j, k, l) = dt*f_compute_cfl_terms(vel, c, j, k, l, .false.)

        ! Viscous calculations (simplified with common patterns)
        if (viscous) then
            if (p > 0) then
                !3D
                if (grid_geometry == 3) then
                    vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/min(dx(j), dy(k), &
                        2._wp*pi*y_cb(max(0,k-1))/max(3._wp, real(min(floor(2._wp*real(max(1,k), wp)*pi), (p + 1)/2 + 1), wp)))**2._wp
                    Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), dy(k)*(abs(vel(2)) + c), &
                        2._wp*pi*y_cb(max(0,k-1))/max(3._wp, real(min(floor(2._wp*real(max(1,k), wp)*pi), (p + 1)/2 + 1), wp))*(abs(vel(3)) + c))/maxval(1._wp/Re_l)
                else
                    vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/min(dx(j), dy(k), dz(l))**2._wp
                    Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), dy(k)*(abs(vel(2)) + c), dz(l)*(abs(vel(3)) + c))/maxval(1._wp/Re_l)
                end if
            elseif (n > 0) then
                !2D
                vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/min(dx(j), dy(k))**2._wp
                Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), dy(k)*(abs(vel(2)) + c))/maxval(1._wp/Re_l)
            else
                !1D
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
    pure subroutine s_compute_dt_from_cfl(vel, c, max_dt, rho, Re_l, j, k, l)
        !$acc routine seq
        real(wp), dimension(num_vels), intent(in) :: vel
        real(wp), intent(in) :: c, rho
        real(wp), dimension(0:m, 0:n, 0:p), intent(inout) :: max_dt
        real(wp), dimension(2), intent(in) :: Re_l
        integer, intent(in) :: j, k, l

        real(wp) :: icfl_dt, vcfl_dt

        ! Inviscid CFL calculation
        icfl_dt = cfl_target*f_compute_cfl_terms(vel, c, j, k, l, .true.)

        ! Viscous calculations (simplified with common patterns)
        if (viscous) then
            if (p > 0) then
                !3D
                if (grid_geometry == 3) then
                    vcfl_dt = cfl_target*min(dx(j), dy(k), &
                        2._wp*pi*y_cb(max(0,k-1))/max(3._wp, real(min(floor(2._wp*real(max(1,k), wp)*pi), (p + 1)/2 + 1), wp)))**2._wp/minval(1/(rho*Re_l))
                else
                    vcfl_dt = cfl_target*min(dx(j), dy(k), dz(l))**2._wp/minval(1/(rho*Re_l))
                end if
            elseif (n > 0) then
                !2D
                vcfl_dt = cfl_target*min(dx(j), dy(k))**2._wp/maxval((1/Re_l)/rho)
            else
                !1D
                vcfl_dt = cfl_target*dx(j)**2._wp/minval(1/(rho*Re_l))
            end if
        end if

        if (any(re_size > 0)) then
            max_dt(j, k, l) = min(icfl_dt, vcfl_dt)
        else
            max_dt(j, k, l) = icfl_dt
        end if

    end subroutine s_compute_dt_from_cfl

end module m_sim_helpers
