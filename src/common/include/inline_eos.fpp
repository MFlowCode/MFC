!>
!! @file inline_eos.fpp
!! @brief Helpers for modules that compile a private copy of the hot EOS routines (see inline_eos_procs.fpp).
!!
!! A kernel module does three things:
!!   1. $:EOS_INLINE_USE('m_eos') in place of its "use m_eos": the public routines are renamed to <name>_ext so they
!!      do not clash with the local copies.
!!   2. $:EOS_INLINE_PRIVATE() after "implicit none": the copies stay private, so they are never re-exported.
!!   3. $:EOS_INLINE_PROCS() right after "contains" (the macro comes from inline_eos_procs.fpp, included below).
!! Call sites stay as they are and resolve to the local copy.

#:include 'inline_eos_procs.fpp'

#:set EOS_INLINE_NAMES = ['f_isentrope_exponent', 'f_isentrope_pressure', 's_phase_coefficients', 'f_bulk_modulus', &
    & 'f_relativistic_enthalpy', 's_compute_mixture_coefficients', 's_compute_speed_of_sound', &
    & 's_compute_speed_of_sound_avg', 's_compute_energy', 'f_is_state_dependent', 'f_pressure', 's_phase_internal_energy', &
    & 's_phase_bulk_modulus', 's_phase_pressure_on_isentrope', 's_phase_temperature', 's_compute_mixture_coefficients_dt']

#:def EOS_INLINE_USE(module, only=None)
    #:if only is None
        use ${module}$, &
    #:else
        #! the copies' own dependencies, which a plain "use m_eos" would bring in
        use ${module}$, only: ${''.join(x + ', ' for x in only)}$s_eos_coefficients, f_has_isentropic_reference, s_reference_curve, &
        & s_rk4, &
    #:endif
    #:for n in EOS_INLINE_NAMES[:-1]
        & ${n}$_ext => ${n}$, &
    #:endfor
    & ${EOS_INLINE_NAMES[-1]}$_ext => ${EOS_INLINE_NAMES[-1]}$
    use m_constants, only: ode_isentrope, ode_reference_temperature
#:enddef

#:def EOS_INLINE_PRIVATE()
    #:for n in EOS_INLINE_NAMES
        private :: ${n}$
    #:endfor
#:enddef
