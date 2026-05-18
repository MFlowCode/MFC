!>
!! @file
!! @brief Contains module m_sim_helpers

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Simulation helper routines for enthalpy computation, CFL calculation, and stability checks
module m_sim_helpers

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion

    implicit none

    private; public :: s_compute_enthalpy, s_compute_stability_from_dt, s_compute_dt_from_cfl

contains

    !> Computes the modified dtheta for Fourier filtering in azimuthal direction
    function f_compute_filtered_dtheta(k, l) result(fltr_dtheta)

        $:GPU_ROUTINE(parallelism='[seq]')
        integer, intent(in) :: k, l
        real(wp)            :: fltr_dtheta
        integer             :: Nfq

        if (grid_geometry == 3) then
            if (k == 0) then
                fltr_dtheta = 2._wp*pi*y_cb(0)/3._wp
            else if (k <= fourier_rings) then
                Nfq = min(floor(2._wp*real(k, wp)*pi), (p + 1)/2 + 1)
                fltr_dtheta = 2._wp*pi*y_cb(k - 1)/real(Nfq, wp)
            else
                fltr_dtheta = y_cb(k - 1)*dz(l)
            end if
        else
            fltr_dtheta = 0._wp
        end if

    end function f_compute_filtered_dtheta

    !> Computes inviscid CFL terms for multi-dimensional cases (2D/3D only)
    function f_compute_multidim_cfl_terms(vel, c, j, k, l) result(cfl_terms)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(:), intent(in) :: vel
        real(wp), intent(in)               :: c
        integer, intent(in)                :: j, k, l
        real(wp)                           :: cfl_terms
        real(wp)                           :: fltr_dtheta

        fltr_dtheta = f_compute_filtered_dtheta(k, l)

        if (p > 0) then
            ! 3D
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                if (grid_geometry == 3) then
                    cfl_terms = min(dx(j)/(abs(vel(1)) + c), dy(k)/(abs(vel(2)) + c), fltr_dtheta/(abs(vel(3)) + c))
                else
                    cfl_terms = min(dx(j)/(abs(vel(1)) + c), dy(k)/(abs(vel(2)) + c), dz(l)/(abs(vel(3)) + c))
                end if
            #:endif
        else
            ! 2D
            cfl_terms = min(dx(j)/(abs(vel(1)) + c), dy(k)/(abs(vel(2)) + c))
        end if

    end function f_compute_multidim_cfl_terms

    !> Computes enthalpy
    subroutine s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, qv, j, k, l)

        $:GPU_ROUTINE(function_name='s_compute_enthalpy',parallelism='[seq]', cray_inline=True)

        type(scalar_field), intent(in), dimension(:) :: q_prim_vf
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), intent(inout), dimension(3) :: alpha
            real(wp), intent(inout), dimension(3) :: vel
        #:else
            real(wp), intent(inout), dimension(:) :: alpha
            real(wp), intent(inout), dimension(:) :: vel
        #:endif
        real(wp), intent(inout)               :: rho, gamma, pi_inf, vel_sum, H, pres
        real(wp), intent(out)                 :: qv
        integer, intent(in)                   :: j, k, l
        real(wp), dimension(2), intent(inout) :: Re
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho, Gs
        #:else
            real(wp), dimension(num_fluids_max) :: alpha_rho, Gs
        #:endif
        real(wp) :: E, G_local
        integer  :: i

        call s_compute_species_fraction(q_prim_vf, j, k, l, alpha_rho, alpha)

        if (elasticity) then
            call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, G_local, Gs)
        else
            call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re)
        end if

        if (igr) then
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_vels
                vel(i) = q_prim_vf(eqn_idx%cont%end + i)%sf(j, k, l)/rho
            end do
        else
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_vels
                vel(i) = q_prim_vf(eqn_idx%cont%end + i)%sf(j, k, l)
            end do
        end if

        vel_sum = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_vels
            vel_sum = vel_sum + vel(i)**2._wp
        end do

        if (igr) then
            E = q_prim_vf(eqn_idx%E)%sf(j, k, l)
            pres = (E - pi_inf - qv - 5.e-1_wp*rho*vel_sum)/gamma
        else
            pres = q_prim_vf(eqn_idx%E)%sf(j, k, l)
            E = gamma*pres + pi_inf + 5.e-1_wp*rho*vel_sum + qv
        end if

        ! Adjust energy for hyperelasticity
        if (hyperelasticity) then
            E = E + G_local*q_prim_vf(eqn_idx%xi%end + 1)%sf(j, k, l)
        end if

        H = (E + pres)/rho

    end subroutine s_compute_enthalpy

    !> Computes stability criterion for a specified dt
    subroutine s_compute_stability_from_dt(vel, c, rho, Re_l, j, k, l, icfl_sf, vcfl_sf, Rc_sf)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in), dimension(:)                        :: vel
        real(wp), intent(in)                                      :: c, rho
        real(wp), dimension(0:m,0:n,0:p), intent(inout)           :: icfl_sf
        real(wp), dimension(0:m,0:n,0:p), intent(inout), optional :: vcfl_sf, Rc_sf
        real(wp), dimension(2), intent(in)                        :: Re_l
        integer, intent(in)                                       :: j, k, l
        real(wp)                                                  :: fltr_dtheta

        ! Inviscid CFL calculation
        if (p > 0 .or. n > 0) then
            ! 2D/3D
            icfl_sf(j, k, l) = dt/f_compute_multidim_cfl_terms(vel, c, j, k, l)
        else
            ! 1D
            icfl_sf(j, k, l) = (dt/dx(j))*(abs(vel(1)) + c)
        end if

        ! Viscous calculations
        if (viscous) then
            if (p > 0) then
                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                    ! 3D
                    if (grid_geometry == 3) then
                        fltr_dtheta = f_compute_filtered_dtheta(k, l)
                        vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/min(dx(j), dy(k), fltr_dtheta)**2._wp
                        Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), dy(k)*(abs(vel(2)) + c), &
                              & fltr_dtheta*(abs(vel(3)) + c))/maxval(1._wp/Re_l)
                    else
                        vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/min(dx(j), dy(k), dz(l))**2._wp
                        Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), dy(k)*(abs(vel(2)) + c), &
                              & dz(l)*(abs(vel(3)) + c))/maxval(1._wp/Re_l)
                    end if
                #:endif
            else if (n > 0) then
                ! 2D
                vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/min(dx(j), dy(k))**2._wp
                Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), dy(k)*(abs(vel(2)) + c))/maxval(1._wp/Re_l)
            else
                ! 1D
                vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/dx(j)**2._wp
                Rc_sf(j, k, l) = dx(j)*(abs(vel(1)) + c)/maxval(1._wp/Re_l)
            end if
        end if

    end subroutine s_compute_stability_from_dt

    !> Computes dt for a specified CFL number
    subroutine s_compute_dt_from_cfl(vel, c, max_dt, rho, Re_l, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(:), intent(in)              :: vel
        real(wp), intent(in)                            :: c, rho
        real(wp), dimension(0:m,0:n,0:p), intent(inout) :: max_dt
        real(wp), dimension(2), intent(in)              :: Re_l
        integer, intent(in)                             :: j, k, l
        real(wp)                                        :: icfl_dt, vcfl_dt
        real(wp)                                        :: fltr_dtheta

        ! Inviscid CFL calculation
        if (p > 0 .or. n > 0) then
            ! 2D/3D cases
            icfl_dt = cfl_target*f_compute_multidim_cfl_terms(vel, c, j, k, l)
        else
            ! 1D case
            icfl_dt = cfl_target*(dx(j)/(abs(vel(1)) + c))
        end if

        ! Viscous calculations
        if (viscous) then
            if (p > 0) then
                ! 3D
                if (grid_geometry == 3) then
                    fltr_dtheta = f_compute_filtered_dtheta(k, l)
                    vcfl_dt = cfl_target*(min(dx(j), dy(k), fltr_dtheta)**2._wp)/maxval(1/(rho*Re_l))
                else
                    vcfl_dt = cfl_target*(min(dx(j), dy(k), dz(l))**2._wp)/maxval(1/(rho*Re_l))
                end if
            else if (n > 0) then
                ! 2D
                vcfl_dt = cfl_target*(min(dx(j), dy(k))**2._wp)/maxval((1/Re_l)/rho)
            else
                ! 1D
                vcfl_dt = cfl_target*(dx(j)**2._wp)/maxval(1/(rho*Re_l))
            end if
        end if

        if (any(Re_size > 0)) then
            max_dt(j, k, l) = min(icfl_dt, vcfl_dt)
        else
            max_dt(j, k, l) = icfl_dt
        end if

    end subroutine s_compute_dt_from_cfl

end module m_sim_helpers
