!>
!! @file m_arrhenius_viscosity.fpp
!! @brief Temperature-dependent (Arrhenius) shear viscosity mu = exp(C + D/T) and the mixture shear Reynolds number it gives.

#:include 'macros.fpp'

!> @brief Temperature-dependent (Arrhenius) shear viscosity mu = exp(C + D/T) and the mixture shear Reynolds number it gives.
module m_arrhenius_viscosity

    use m_derived_types      !< Definitions of the derived types
    use m_global_parameters  !< sgm_eps, visc_models, visc_cs, visc_ds, fluid_inv_re
    use m_constants, only: visc_exp_arg_cap
    use m_eos, only: f_mixture_temperature

    implicit none

    private; public :: f_compute_arrhenius_shear_re

contains

    !> Mixture shear Reynolds number 1/sum(alpha_i*mu_i) of one state at its thermal-equilibrium mixture temperature, the same
    !! closure heat conduction uses. Arrhenius fluids take mu_i = exp(C_i + D_i/T); the others keep their constant 1/Re_i(1).
    function f_compute_arrhenius_shear_re(alpha, alpha_rho, pres, gamma, pi_inf) result(Re_shear)

        $:GPU_ROUTINE(function_name='f_compute_arrhenius_shear_re', parallelism='[seq]', cray_inline=True)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha, alpha_rho
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha, alpha_rho
        #:endif
        real(wp), intent(in) :: pres, gamma, pi_inf
        real(wp)             :: Re_shear
        real(wp)             :: T, inv_re
        integer              :: fl

        T = max(f_mixture_temperature(alpha_rho, pres, gamma, pi_inf), sgm_eps)

        inv_re = 0._wp
        do fl = 1, num_fluids
            if (visc_models(fl) == visc_model_arrhenius) then
                inv_re = inv_re + alpha(fl)*exp(min(visc_cs(fl) + visc_ds(fl)/T, visc_exp_arg_cap))
            else
                inv_re = inv_re + alpha(fl)*fluid_inv_re(fl)
            end if
        end do

        Re_shear = 1._wp/max(inv_re, sgm_eps)

    end function f_compute_arrhenius_shear_re

end module m_arrhenius_viscosity
