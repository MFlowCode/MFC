!>
!! @file m_hb_function.f90
!! @brief Herschel-Bulkley non-Newtonian viscosity: formula, shear rate, and mixture inverse-Re.

#:include 'macros.fpp'

!> @brief Herschel-Bulkley non-Newtonian viscosity: formula, shear rate, and mixture inverse-Re.
module m_hb_function

    use m_derived_types      !< Definitions of the derived types
    use m_global_parameters  !< Re_size, Re_idx, sgm_eps, dflt_real, any_non_newtonian, hb_* arrays
    use m_constants          !< verysmall

    implicit none

    private; public :: f_compute_hb_viscosity, f_compute_shear_rate_from_components, s_compute_mixture_inv_re

contains

    !> Papanastasiou-regularized Herschel-Bulkley viscosity.
    function f_compute_hb_viscosity(tau0, K_val, nn_val, mu_min_val, mu_max_val, shear_rate, hb_m_val) result(mu)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: tau0, K_val, nn_val
        real(wp), intent(in) :: mu_min_val, mu_max_val
        real(wp), intent(in) :: shear_rate, hb_m_val
        real(wp)             :: mu
        real(wp)             :: yield_term, power_law_term, g_eff

        g_eff = max(shear_rate, verysmall)
        if (shear_rate <= verysmall) then
            yield_term = tau0*hb_m_val
        else
            yield_term = tau0*(1._wp - exp(-hb_m_val*shear_rate))/shear_rate
        end if
        power_law_term = K_val*(g_eff**(nn_val - 1._wp))

        mu = yield_term + power_law_term
        if (mu_min_val > dflt_real) mu = max(mu, mu_min_val)
        if (mu_max_val > dflt_real) mu = min(mu, mu_max_val)

    end function f_compute_hb_viscosity

    !> Shear rate gamma_dot = sqrt(2 D_ij D_ij). Absent dims pass 0.
    function f_compute_shear_rate_from_components(D_xx, D_yy, D_zz, D_xy, D_xz, D_yz) result(shear_rate)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: D_xx, D_yy, D_zz, D_xy, D_xz, D_yz
        real(wp)             :: shear_rate

        shear_rate = sqrt(2._wp*(D_xx*D_xx + D_yy*D_yy + D_zz*D_zz + 2._wp*(D_xy*D_xy + D_xz*D_xz + D_yz*D_yz)))

    end function f_compute_shear_rate_from_components

    !> Mixture inverse Reynolds (= 1/mu_mix) per direction (1=shear, 2=bulk) at one state. Newtonian fluids accumulate the legacy
    !! per-fluid alpha/Res term. For non-Newtonian fluids only the shear direction (i==1) uses the Herschel-Bulkley viscosity at the
    !! given shear rate; the bulk direction (i==2) contributes zero, since non-Newtonian bulk viscosity is not supported.
    subroutine s_compute_mixture_inv_re(alpha, shear_rate, Res, Re_out)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(*), intent(in)   :: alpha
        real(wp), intent(in)                 :: shear_rate
        real(wp), dimension(2,*), intent(in) :: Res
        real(wp), dimension(2), intent(out)  :: Re_out
        integer                              :: i, q, fl
        real(wp)                             :: mu_q

        ! Plain serial loops: this is already a seq device routine, so no GPU_LOOP
        ! directives (they emit empty on Cray/AMD and are redundant here).
        do i = 1, 2
            Re_out(i) = dflt_real
            if (Re_size(i) > 0) Re_out(i) = 0._wp
            do q = 1, Re_size(i)
                fl = Re_idx(i, q)
                if (any_non_newtonian .and. is_non_newtonian(fl)) then
                    if (i == 1) then
                        mu_q = f_compute_hb_viscosity(hb_tau0(fl), hb_K(fl), hb_nn(fl), hb_mu_min(fl), hb_mu_max(fl), shear_rate, &
                                                      & hb_m_arr(fl))
                    else
                        mu_q = 0._wp
                        if (hb_mu_bulk(fl) > dflt_real) mu_q = hb_mu_bulk(fl)
                    end if
                    Re_out(i) = alpha(fl)*mu_q + Re_out(i)
                else
                    Re_out(i) = alpha(fl)/Res(i, q) + Re_out(i)
                end if
            end do
            Re_out(i) = 1._wp/max(Re_out(i), sgm_eps)
        end do

    end subroutine s_compute_mixture_inv_re

end module m_hb_function
