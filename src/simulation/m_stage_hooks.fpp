!>
!!@file
!!@brief Contains module m_stage_hooks

#:include 'macros.fpp'

!> @brief The hooks applied to the conservative state after each RK stage update, as one ordered list shared by the coarse grid
!! (m_time_steppers) and the AMR fine blocks (m_amr_advance), so the two levels cannot run them in different orders.
module m_stage_hooks

    use m_derived_types
    use m_global_parameters
    use m_body_forces, only: s_apply_bodyforces, s_apply_synthetic_turbulence_force
    use m_fftw, only: s_apply_fourier_filter
    use m_pressure_relaxation, only: s_pressure_relaxation_procedure
    use m_bubbles_EE, only: s_comp_alpha_from_n
    use m_hypoelastic, only: s_enforce_cont_damage_bounds
    use m_species_bounds, only: s_clean_species
    use m_constants, only: model_eqns_6eq

    implicit none

    private
    public :: s_stage_hooks_pre_ib, s_stage_hooks_post_ib, f_stage_hooks_pre_ib_on_block, f_stage_hooks_post_ib

contains

    !> Sources and closures, before the IB correction. on_block: q_cons is an AMR fine block in its own frame, where the azimuthal
    !! filter does not apply (it is a whole-ring transform; case_validator rejects cylindrical grids with amr). q_prim is read by
    !! the synthetic turbulence forcing only.
    impure subroutine s_stage_hooks_pre_ib(q_cons, q_prim, rhs, ldt, on_block, q_T_sf)

        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_cons, rhs
        type(scalar_field), dimension(1:sys_size), intent(in)    :: q_prim
        real(wp), intent(in)                                     :: ldt  !< the stage's rk_coef(s, 3)*dt/rk_coef(s, 4)
        logical, intent(in)                                      :: on_block
        type(scalar_field), intent(inout), optional              :: q_T_sf

        if (bodyForces) call s_apply_bodyforces(q_cons, q_prim, rhs, ldt, q_T_sf)
        if (synthetic_turbulence) call s_apply_synthetic_turbulence_force(q_cons, q_prim, rhs, ldt)
        if (grid_geometry == 3 .and. .not. on_block) call s_apply_fourier_filter(q_cons)
        if (model_eqns == model_eqns_6eq .and. (.not. relax)) call s_pressure_relaxation_procedure(q_cons)
        if (adv_n) call s_comp_alpha_from_n(q_cons)

    end subroutine s_stage_hooks_pre_ib

    !> Clamps, after the IB correction. lvl: the AMR level of the frame (0 = the coarse grid), for the species cleanup's report.
    impure subroutine s_stage_hooks_post_ib(q_cons, lvl)

        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_cons
        integer, intent(in)                                      :: lvl

        if (cont_damage) call s_enforce_cont_damage_bounds(q_cons)
        if (chemistry) call s_clean_species(q_cons, lvl)

    end subroutine s_stage_hooks_post_ib

    !> Whether s_stage_hooks_pre_ib does anything on a fine block (it costs a frame swap there).
    pure logical function f_stage_hooks_pre_ib_on_block() result(on)

        on = bodyForces .or. synthetic_turbulence .or. (model_eqns == model_eqns_6eq .and. (.not. relax)) .or. adv_n

    end function f_stage_hooks_pre_ib_on_block

    !> Whether s_stage_hooks_post_ib does anything.
    pure logical function f_stage_hooks_post_ib() result(on)

        on = cont_damage .or. chemistry

    end function f_stage_hooks_post_ib

end module m_stage_hooks
