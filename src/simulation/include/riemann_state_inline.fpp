!>
!! @file riemann_state_inline.fpp
!! @brief Rename / private helpers for modules that compile a private copy of riemann_state_inline_procs.fpp.
!! Usage as in common/include/inline_eos.fpp.

#:include 'riemann_state_inline_procs.fpp'

#:set RS_INLINE_NAMES = ['f_low_Mach_zcoef', 'f_low_Mach_pcorr_hllc', 's_apply_low_Mach_velocity', 's_compute_average_state', &
    & 's_compute_hypoelastic_interface_energy', 'f_compute_hllc_star_momentum_flux']

#:def RS_INLINE_USE(module)
    use ${module}$, &
    #:for n in RS_INLINE_NAMES[:-1]
        & ${n}$_ext => ${n}$, &
    #:endfor
    & ${RS_INLINE_NAMES[-1]}$_ext => ${RS_INLINE_NAMES[-1]}$
#:enddef

#:def RS_INLINE_PRIVATE()
    #:for n in RS_INLINE_NAMES
        private :: ${n}$
    #:endfor
#:enddef
