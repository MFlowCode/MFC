!>
!! @file
!> @brief Fused OpenACC WENO/HLLC/RHS kernel for five- and six-equation mixtures

#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'hll_velocity.fpp'
#! MFC_FUSION_STAGING_DEFAULT comes from the toolchain's case.fpp; plain-CMake renders (docs) take the literal, mirroring tuning.py.
#:set MFC_FUSION_STAGING_BYTES = MFC_FUSION_STAGING_BYTES if defined('MFC_FUSION_STAGING_BYTES') &
    & else (MFC_FUSION_STAGING_DEFAULT if defined('MFC_FUSION_STAGING_DEFAULT') else 32*1024)
#:set MFC_FUSION_X_PADDED_CELLS = MFC_FUSION_X_PADDED_CELLS if defined('MFC_FUSION_X_PADDED_CELLS') else 16
#:set MFC_FUSION_TRANSVERSE_CELLS = MFC_FUSION_TRANSVERSE_CELLS if defined('MFC_FUSION_TRANSVERSE_CELLS') else 4
#:set MFC_FUSION_TRANSVERSE_WIDE_CELLS = MFC_FUSION_TRANSVERSE_WIDE_CELLS if defined('MFC_FUSION_TRANSVERSE_WIDE_CELLS') else 4
#:set MFC_FUSION_TRANSVERSE_LANES = MFC_FUSION_TRANSVERSE_LANES if defined('MFC_FUSION_TRANSVERSE_LANES') else 2
#:set MFC_FUSION_TRANSVERSE_WIDE_LANES = MFC_FUSION_TRANSVERSE_WIDE_LANES if defined('MFC_FUSION_TRANSVERSE_WIDE_LANES') else 2
#:set MFC_FUSION_COMPACT_LANES = MFC_FUSION_COMPACT_LANES if defined('MFC_FUSION_COMPACT_LANES') else 2
#:set MFC_FUSION_CHEM_LANES = MFC_FUSION_CHEM_LANES if defined('MFC_FUSION_CHEM_LANES') else (1, 2)
#:set MFC_FUSION_VECTOR_X = MFC_FUSION_VECTOR_X if defined('MFC_FUSION_VECTOR_X') else 128
#:set MFC_FUSION_VECTOR_TRANSVERSE = MFC_FUSION_VECTOR_TRANSVERSE if defined('MFC_FUSION_VECTOR_TRANSVERSE') else 64
#:set MFC_FUSION_VECTOR_COMPACT = MFC_FUSION_VECTOR_COMPACT if defined('MFC_FUSION_VECTOR_COMPACT') else 128
#:set MFC_FUSION_VECTOR_CHEM_PER_LANE = MFC_FUSION_VECTOR_CHEM_PER_LANE if defined('MFC_FUSION_VECTOR_CHEM_PER_LANE') else 32
#:include 'mhd_face.fpp'
#:include 'hllc_burn.fpp'
#:include 'hypoelastic_energy.fpp'
#:include 'chemistry_face.fpp'

#:def hllc_rk_update(FIELD)
    #:set DMG_ROW = NADV + STRESS_ROWS + 1 if cont_damage else 0
    #:if hllc_reactive_burn
        ! With rburn%substeps > 0, s_reactive_burn_substep integrates the burn after the flow update.
        if (rburn%substeps == 0) then
            $:hllc_burn_rate(f'real({FIELD("stage", ENERGY)}, wp)', f'real({FIELD("stage", ENERGY + 2)}, wp)', &
                             & f'real({FIELD("stage", 1)}, wp)', f'real({FIELD("stage", ENERGY + 1)}, wp)')
            burn_mass_rate = (real(${FIELD('stage', 1)}$, wp) + real(${FIELD('stage', 2)}$, wp))*burn_rate
            value1 = value1 - burn_mass_rate; value2 = value2 + burn_mass_rate
            value${ENERGY + 1}$ = value${ENERGY + 1}$ - burn_rate
            value${ENERGY + 2}$ = value${ENERGY + 2}$ + burn_rate
        end if
    #:endif
    if (save_rk_state) then
        #:for EQ in range(1, NE + 1)
            ${FIELD('stor', EQ)}$ = ${FIELD('cons', EQ)}$
        #:endfor
    end if
    #:for EQ in range(1, NE + 1)
        ${FIELD('output', EQ)}$ = real((rk_q_coef*real(${FIELD('cons', EQ)}$, wp) + rk_stor_coef*real(${FIELD('stor', EQ)}$, &
                & wp) + rk_rhs_coef*value${EQ}$)/rk_den_coef, stp)
    #:endfor
    #:for F in range(1, NF + 1)
        arho${F}$ = real(${FIELD('output', F)}$, wp)
        alpha${F}$ = real(${FIELD('output', ENERGY + F)}$, wp)
    #:endfor
    #:if hllc_mpp_lim
        #:for EPS in ['1.e-16_wp', 'sgm_eps']
            #:for F in range(1, NF + 1)
                arho${F}$ = max(0._wp, arho${F}$)
                alpha${F}$ = min(max(0._wp, alpha${F}$), 1._wp)
            #:endfor
            alpha_sum = ${' + '.join(f'alpha{f}' for f in range(1, NF + 1))}$
            #:for F in range(1, NF + 1)
                alpha${F}$ = alpha${F}$/max(alpha_sum, ${EPS}$)
            #:endfor
        #:endfor
    #:endif
    rho = max(${' + '.join(f'arho{f}' for f in range(1, NF + 1))}$, sgm_eps)
    #:if any_state_dependent_eos
        ! Each phase's coefficients at its own density, as in s_compute_mixture_coefficients.
        gamma_mix = 0._wp
        pi_mix = 0._wp
        #:for F in range(1, NF + 1)
            call s_phase_coefficients(arho${F}$, alpha${F}$, ${F}$, eos_rho, eos_gamma, eos_pi, eos_dpi, eos_dgamma)
            gamma_mix = gamma_mix + alpha${F}$*eos_gamma
            pi_mix = pi_mix + alpha${F}$*eos_pi
        #:endfor
    #:else
        gamma_mix = ${' + '.join(f'alpha{f}*gamma{f}' for f in range(1, NF + 1))}$
        pi_mix = ${' + '.join(f'alpha{f}*pi{f}' for f in range(1, NF + 1))}$
    #:endif
    #:if not qv_zero
        qv_mix = ${' + '.join(f'arho{f}*qv{f}' for f in range(1, NF + 1))}$
    #:endif
    inv_rho = 1._wp/rho
    #:for V in range(1, NV + 1)
        velocity${V}$ = real(${FIELD('output', NF + V)}$, wp)*inv_rho
    #:endfor
    dynamic_pressure = ${' + '.join(f'0.5_wp*real({FIELD("output", NF + v)}, wp)*velocity{v}' for v in range(1, NV + 1))}$
    #:if mhd
        magnetic_pressure = 0.5_wp*(${'magnetic_bx**2 + ' if ND == 1 else ''}$ ${' + '.join(f'real({FIELD("output", NADV + row)}, wp)**2' for row in range(1, (2 if ND == 1 else 3) + 1))}$)
    #:endif
    pressure = (real(${FIELD('output', ENERGY)}$, &
                & wp) - dynamic_pressure - pi_mix${' - qv_mix' if not qv_zero else ''}$${' - magnetic_pressure' if mhd else ''}$) &
                & /gamma_mix
    #:if mhd
        #:for ROW in range(1, MAG_ROWS + 1)
            ${FIELD('prim', NADV + ROW)}$ = ${FIELD('output', NADV + ROW)}$
        #:endfor
    #:endif
    #:if hypoelasticity
        ! The elastic energy uses the undamaged mixture modulus; damage scales the wave-speed moduli.
        shear_mix = max(${' + '.join(f'alpha{f}*Gs_rs({f})' for f in range(1, NF + 1))}$, 0._wp)
        #:for ROW in range(1, STRESS_ROWS + 1)
            tau = real(${FIELD('output', NADV + ROW)}$, wp)*inv_rho
            ${FIELD('prim', NADV + ROW)}$ = real(tau, stp)
            if (shear_mix > 1.e-12_wp) then
                elastic_pressure = tau*tau*(1._wp/(4._wp*shear_mix))*(1._wp/gamma_mix)
                pressure = pressure - elastic_pressure
                #:if ROW in (2, 4, 5) and ROW <= ND*(ND + 1)//2
                    pressure = pressure - elastic_pressure
                #:endif
            end if
        #:endfor
        #:if cont_damage
            ! D = U_D/m_s (damageable-solid partial mass), clamped to [0, 1].
            solid_partial_density = 0._wp
            #:for F in range(1, NF + 1)
                if (Gs_rs(${F}$) > verysmall) solid_partial_density = solid_partial_density + real(${FIELD('output', F)}$, wp)
            #:endfor
            ${FIELD('prim', DMG_ROW)}$ = real(min(max(real(${FIELD('output', DMG_ROW)}$, wp)/max(solid_partial_density, &
                    & verysmall), 0._wp), 1._wp), stp)
        #:endif
    #:endif
    #:for F in range(1, NF + 1)
        ${FIELD('prim', F)}$ = ${FIELD('output', F)}$
        ${FIELD('prim', ENERGY + F)}$ = ${FIELD('output', ENERGY + F)}$
    #:endfor
    #:for V in range(1, NV + 1)
        ${FIELD('prim', NF + V)}$ = real(velocity${V}$, stp)
    #:endfor
    ${FIELD('prim', ENERGY)}$ = real(pressure, stp)
    #:if surface_tension
        ${FIELD('prim', NE)}$ = ${FIELD('output', NE)}$
    #:endif
#:enddef

module m_fused_rhs_openacc

#ifdef MFC_FUSED_RHS
    use m_derived_types
    use m_global_parameters, only: m, n, p, idwbuff, dx, dy, dz, y_cc, y_cb, gammas, pi_infs, qvs, mixture_err, bc_x, bc_y, bc_z, &
        & buff_size
    use m_global_parameters, only: rburn, isentrope_B, isentrope_n, cvs, ic_eps, ic_beta
    #:if MFC_CASE_OPTIMIZATION and cont_damage
        use m_global_parameters, only: cont_damage
    #:endif
    #:if MFC_CASE_OPTIMIZATION and viscous
        use m_global_parameters, only: fluid_pp, shear_stress, bulk_stress
    #:endif
    #:if MFC_CASE_OPTIMIZATION and hllc_non_newtonian
        use m_global_parameters, only: hb_tau0, hb_K, hb_nn, hb_mu_min, hb_mu_max, hb_m_arr
        use m_hb_function, only: f_compute_hb_viscosity, f_compute_shear_rate_from_components
    #:endif
    #:if MFC_CASE_OPTIMIZATION and hllc_int_comp == 2 and num_dims > 1
        use m_thinc, only: mthinc_nhat, mthinc_d, gq3_pts, gq3_wts
    #:endif
    #:if MFC_CASE_OPTIMIZATION and bubbles_euler
        use m_global_parameters, only: R0, weight
        use m_bubbles, only: f_cpbw_KM
        #:if polytropic and not qbmm and not viscous
            use m_global_parameters, only: Ca, Web, Re_inv, gam, Eu
            use m_helper_basic, only: f_is_default
        #:endif
        use m_constants, only: pi, small_alf
        #:if qbmm
            use m_global_parameters, only: mom_sp
        #:endif
    #:endif
    #:if MFC_CASE_OPTIMIZATION and mhd
        use m_global_parameters, only: Bx0, hyper_cleaning_speed
    #:endif
    use m_constants, only: sgm_eps, verysmall, moncon_cutoff, BC_CHAR_SLIP_WALL, BC_CHAR_SUP_OUTFLOW, BC_RIEMANN_EXTRAP
    #:if MFC_CASE_OPTIMIZATION and hypoelasticity
        use m_riemann_state, only: Gs_rs
    #:endif
    #:if MFC_CASE_OPTIMIZATION and chemistry
        use m_thermochem, only: gas_constant
        #:if not hllc_scalar_thermo
            use m_thermochem, only: molecular_weights, get_mixture_molecular_weight, get_species_specific_heats_r, &
                & get_mixture_specific_heat_cp_mass, get_mixture_specific_heat_cv_mass, get_mixture_energy_mass, &
                & get_mixture_viscosity_mixavg
        #:endif
    #:endif
    use m_riemann_state, only: flux_rsx_vf, flux_src_rsx_vf
    #:if MFC_CASE_OPTIMIZATION and any_state_dependent_eos
        use m_eos, only: f_bulk_modulus, s_phase_coefficients, s_phase_internal_energy, s_phase_pressure_on_isentrope, &
            & s_phase_bulk_modulus, s_phase_temperature
    #:else
        use m_eos, only: f_bulk_modulus
    #:endif
    #:if MFC_CASE_OPTIMIZATION and ((viscous and not hllc_viscous_stress) or surface_tension)
        use m_riemann_state, only: vel_src_rsx_vf, Re_avg_rsx_vf
    #:endif
#ifdef MFC_ACTIVE_BOX
    use m_active_box, only: ab_x, ab_y, ab_z, ab_active, ab_nb, ab_box, ab_inner
#endif

    implicit none

    private

    #:if MFC_CASE_OPTIMIZATION and hllc_alt_soundspeed
        public :: f_hllc_kapila
    #:endif

    #:if MFC_CASE_OPTIMIZATION and 1 <= num_dims <= 3 and 1 <= num_fluids <= 4 and (model_eqns == 2 or (model_eqns == 3 and hllc_riemann == 2 and hllc_low_mach == 0 and not hllc_alt_soundspeed)) and (not bubbles_euler or hllc_bubbles) and (not chemistry or hllc_chemistry) and (not mhd or hllc_mhd) and (not relativity or hllc_mhd) and (not hypoelasticity or hllc_hypoelastic) and (not surface_tension or (hllc_riemann == 2 and not viscous))
        #:set NF = num_fluids
        #:set ND = num_dims
        #:set NV = 3 if mhd else ND
        #:set MAG_ROWS = (2 if ND == 1 else 3) + int(hyper_cleaning) if mhd else 0
        #:set ENERGY = NF + NV + 1
        #:set SIXEQ = model_eqns == 3
        #:set NADV = ENERGY + NF
        #:set BUB_STRIDE = 6 if qbmm else 2 if polytropic else 4
        #:set BUB_ROWS = nb*BUB_STRIDE if bubbles_euler else 0
        #:set BUB_SIZE = BUB_ROWS + (1 if bubbles_euler and adv_n else 0)
        #:set ELASTIC = hypoelasticity
        #:set STRESS_ROWS = ND*(ND + 1)//2 + int(cyl_coord) if ELASTIC else 0
        #! Continuum damage adds the advected row U_D = m_s*D behind the stresses.
        #:set DMG_ROWS = 1 if hypoelasticity and cont_damage else 0
        #:set STRESS_INDEX = lambda a, b: max(a, b)*(max(a, b) - 1)//2 + min(a, b)
        #:set CHEM_ROWS = hllc_num_species if chemistry else 0
        #:set NE = NADV + (NF if SIXEQ else 0) + (1 if surface_tension else 0) + BUB_SIZE + STRESS_ROWS &
                           & + DMG_ROWS + CHEM_ROWS + MAG_ROWS
        #! A state-dependent EOS reads per-phase coefficients, which the bubbly closures and the
        #! two-fluid Kapila term take from gammas and pi_infs.
        #:assert not (any_state_dependent_eos and bubbles_euler)
        #:assert not (any_state_dependent_eos and hllc_alt_soundspeed and NF != 2)
        #:set EARLY_CHEM = chemistry and not viscous and hllc_riemann in (1, 2, 5) and hllc_wave_speeds == 1 and hllc_low_mach == 0
        #:set EARLY_BUBBLE = bubbles_euler and polytropic and not qbmm and not viscous
        #:set EARLY_QBMM = bubbles_euler and qbmm and not viscous and hllc_riemann == 2 and hllc_wave_speeds == 1 and hllc_low_mach == 0 and not hllc_alt_soundspeed and not cyl_coord
        #:set NREC = NADV + (1 if surface_tension else 0) + BUB_SIZE + STRESS_ROWS &
                      & + DMG_ROWS + CHEM_ROWS &
                      & + MAG_ROWS
        #:set NRIEMANN = NADV + 1 if bubbles_euler and qbmm else NREC
        ! Colour follows the fractions in reconstruction, skipping six-equation internal energies.
        #:set REC_FIELD = lambda slot: NE if surface_tension and slot == NREC else slot
        #:set STATE = lambda side, &
            & eq: f'density_{side}{eq}' if eq <= NF else (f'velocity_{side}({eq - NF})' if eq < ENERGY else (f'pressure_{side}' if eq == ENERGY else f'alpha_{side}{eq - ENERGY}'))
        #:set HLL_FAMILY = hllc_riemann in (1, 5)
        #:set SOURCE_ROWS = range(1, NF + 1) if hllc_riemann == 1 else [1]
        #:set SOURCE = lambda f: f'source({f})' if hllc_riemann == 1 else 'source'
        public :: s_fused_rhs_openacc, s_weno_fused_rhs_openacc_kernel, s_fused_rhs_openacc_overlap_bounds
        #:if not SIXEQ
            public :: s_fused_rhs_openacc_boundary_rk, s_fused_rhs_openacc_update_rk
        #:endif
        #:if MFC_CASE_OPTIMIZATION and viscous and hllc_vel_gradients
            public :: ${', '.join(f's_fused_velocity_gradients_{a}' for a in 'xyz'[:num_dims])}$
        #:endif
    #:endif

contains

    #:if MFC_CASE_OPTIMIZATION and hllc_int_comp == 2 and num_dims > 1
        #:include 'thinc_integrals.fpp'
    #:endif

    #:if MFC_CASE_OPTIMIZATION and hllc_alt_soundspeed
        real(wp) function f_hllc_kapila(pressure, alpha_one, alpha_two, alpha_rho_one, alpha_rho_two, gamma_one, gamma_two, &
             & pi_one, pi_two) result(k_term)
            $:GPU_ROUTINE(function_name='f_hllc_kapila', parallelism='[seq]')

            real(wp), value, intent(in) :: pressure, alpha_one, alpha_two, alpha_rho_one, alpha_rho_two
            real(wp), value, intent(in) :: gamma_one, gamma_two, pi_one, pi_two
            real(wp)                    :: bulk_one, bulk_two

            #:if any_state_dependent_eos
                ! Each phase's bulk modulus at its own density.
                call s_phase_bulk_modulus(pressure, alpha_one, alpha_rho_one, 1, bulk_one)
                call s_phase_bulk_modulus(pressure, alpha_two, alpha_rho_two, 2, bulk_two)
            #:else
                bulk_one = ((gamma_one + 1._wp)*pressure + pi_one)/gamma_one
                bulk_two = ((gamma_two + 1._wp)*pressure + pi_two)/gamma_two
            #:endif
            k_term = alpha_one*alpha_two*(bulk_two - bulk_one)/(alpha_one*bulk_two + alpha_two*bulk_one)

        end function f_hllc_kapila
    #:endif

    #:if MFC_CASE_OPTIMIZATION and 1 <= num_dims <= 3 and 1 <= num_fluids <= 4 and (model_eqns == 2 or (model_eqns == 3 and hllc_riemann == 2 and hllc_low_mach == 0 and not hllc_alt_soundspeed)) and (not bubbles_euler or hllc_bubbles) and (not chemistry or hllc_chemistry) and (not mhd or hllc_mhd) and (not relativity or hllc_mhd) and (not hypoelasticity or hllc_hypoelastic) and (not surface_tension or (hllc_riemann == 2 and not viscous))
        #:set RK_FIELDS = ['cons', 'stor', 'output', 'prim']
        #:set RK_ARRAY_GROUPS = [', '.join(f'{field}{eq}' for eq in range(1, NE + 1)) for field in RK_FIELDS]
        #:set RK_SCALARS = 'save_rk_state, rk_q_coef, rk_stor_coef, rk_rhs_coef, rk_den_coef'
        #:set RK_GROUPS = RK_ARRAY_GROUPS + [RK_SCALARS]
        #:set RK_PRIVATE = [f'value{eq}' for eq in range(1, NE + 1)] + [f'{field}{f}' for field in ['arho', &
            & 'alpha'] for f in range(1, NF + 1)] + [f'velocity{v}' for v in range(1, NV + 1)] + ['alpha_sum', 'rho', 'inv_rho', &
            & 'gamma_mix', 'pi_mix', 'qv_mix', 'dynamic_pressure', 'pressure', 'burn_rate', 'burn_mass_rate', 'burn_drive', &
            & 'burn_temperature'] + (['eos_rho', 'eos_gamma', 'eos_pi', 'eos_dpi', 'eos_dgamma'] if any_state_dependent_eos else []) + (['shear_mix', 'tau', 'elastic_pressure'] + (['solid_partial_density'] if cont_damage else []) if hypoelasticity else []) + (['magnetic_pressure'] if mhd else [])
        #:include 'fused_weno.fpp'
        #:set X_COLUMN_CELLS = MFC_FUSION_X_PADDED_CELLS
        #:if chemistry
            #:set X_VISC_ROWS = ND + 1 if viscous and hllc_viscous_stress else 0
            #:set X_COLUMN_CELLS = min(X_COLUMN_CELLS, &
                                     & (MFC_FUSION_STAGING_BYTES//8 - (2*WENO_RADIUS + 4 if WENO_RADIUS > 0 else 2)*NREC &
                                     & - X_VISC_ROWS)//((3 if WENO_RADIUS > 0 else 2)*NREC + X_VISC_ROWS))
        #:endif

        subroutine s_fused_rhs_openacc_overlap_bounds(first, last)

            integer, intent(out) :: first, last

            first = ((buff_size + ${X_COLUMN_CELLS}$)/${X_COLUMN_CELLS}$)*${X_COLUMN_CELLS}$
            last = ((m - buff_size)/${X_COLUMN_CELLS}$)*${X_COLUMN_CELLS}$ - 1

        end subroutine s_fused_rhs_openacc_overlap_bounds

        #:set COMPRESS = NF == 2 and WENO_RADIUS > 0 and hllc_int_comp in (1, 2)
        #:set MTHINC = hllc_int_comp == 2 and ND > 1
        #:set PHYS_GROUPS = [', '.join(f'gamma{f}' for f in range(1, NF + 1)), ', '.join(f'pi{f}' for f in range(1, NF + 1)), &
                                       & ', '.join(f'qv{f}' for f in range(1, NF + 1)), 'mixture_error']
        ! Non-Newtonian Reynolds numbers come from the face strain rate, so the per-fluid
        ! Herschel-Bulkley parameters are launch arguments.
        #:set VISC_REAL_GROUPS = [', '.join(f're_{stress}{f}' for f in range(1, NF + 1)) for stress in ['shear', &
                                            & 'bulk']] if viscous else []
        #:set NN_GROUPS = [', '.join(f'hb_tau0_{f}, hb_K_{f}, hb_nn_{f}, hb_mu_min_{f}, hb_mu_max_{f}, ' &
                                     & + f'hb_m_{f}' for f in range(1, NF + 1))] if viscous and hllc_non_newtonian else []
        #:set MHD_GROUPS = ['magnetic_bx, cleaning_speed'] if mhd else []
        #:set PHYS_GROUPS = PHYS_GROUPS + MHD_GROUPS + VISC_REAL_GROUPS + (['shear_active, bulk_active'] &
            & if viscous else []) + NN_GROUPS
        #:set BUB_FLAGS = ', bubble_web_active, bubble_re_active' if EARLY_BUBBLE else ''
        #:set PHYS_GROUPS = PHYS_GROUPS + (['bubble_web_active, bubble_re_active'] if EARLY_BUBBLE else [])
        #:set BUB_HOST_FLAGS = ', .not. f_is_default(Web), .not. f_is_default(Re_inv)' if EARLY_BUBBLE else ''
        #:set QBMM_GROUPS = ['moment1, moment3, moment4', 'mi, mj'] if bubbles_euler and qbmm else []
        #:set MOMENT_ARGS = ', '.join('moment_' + side + str(i) for side in ['l', 'r'] for i in [1, 3, 4])
        #:set BUB_RESULTS = 'bubble_transport, ' if bubbles_euler and qbmm else ''
        #:set VISC_RESULTS = ('viscous_velocity, ' + ('re_interface, ' if viscous and not hllc_non_newtonian else '')) &
                              & if viscous or surface_tension else ''
        #:set VISC_HOST_ARGS = (', ' + ', '.join(f'fluid_pp({f})%Re({s})' for s in [1, 2] for f in range(1, &
                                & NF + 1)) + ', shear_stress, bulk_stress') if viscous else ''
        #:set NN_HOST_ARGS = (', ' + ', '.join(f'{arr}({f})' for f in range(1, NF + 1) for arr in ['hb_tau0', 'hb_K', &
                                & 'hb_nn', 'hb_mu_min', 'hb_mu_max', 'hb_m_arr'])) if viscous and hllc_non_newtonian else ''
        #:set IN_COLUMN_STRESS = viscous and hllc_viscous_stress
        #! 3D cylindrical stress follows s_compute_cylindrical_viscous_source_flux.
        #:set CYL_STRESS = cyl_coord and num_dims > 2
        #! weno_Re_flux reconstructs the face gradients in the column; other viscous cases read
        #! the face gradients the host passes.
        #:set GRAD_WENO_FLUX = hllc_grad_weno
        #:set GRAD_FIELDS = [f'grad_{a}_vf' for a in 'xyz'[:ND]] if IN_COLUMN_STRESS else []
        #:set GRAD_GROUPS = [', '.join(f'grad_{a}{v}' for a in 'xyz'[:ND] for v in range(1, ND + 1)), &
                                       & 'gi, gj'] if IN_COLUMN_STRESS else []
        #:if bubbles_euler
            subroutine s_hllc_bubble_state(q, rho, pressure, number_density, &
                                           & p_tilde${BUB_FLAGS}$${', moment_r3, moment_r3v2, moment_pbwr3' if qbmm else ''}$)
                $:GPU_ROUTINE(function_name='s_hllc_bubble_state', parallelism='[seq]')
                real(wp), intent(in)        :: q(${NRIEMANN}$)
                real(wp), value, intent(in) :: rho, pressure
                real(wp), intent(out)       :: number_density, p_tilde
                real(wp)                    :: r3_sum, r3v2_sum, pbwr3_sum
                #:if EARLY_BUBBLE
                    logical, value, intent(in) :: bubble_web_active, bubble_re_active
                    real(wp)                   :: radius_b, r0_b, velocity_b, ratio_b, ca_b, web_b, re_b
                #:endif
                #:if qbmm
                    real(wp), value, intent(in) :: moment_r3, moment_r3v2, moment_pbwr3
                    number_density = q(${NADV + 1}$)
                    r3_sum = moment_r3; r3v2_sum = moment_r3v2; pbwr3_sum = moment_pbwr3
                #:else
                    real(wp) :: radius, velocity, cubed, pbw
                    number_density = ${'q(' + str(NE) + ')' if adv_n else '0._wp'}$
                    r3_sum = 0._wp; r3v2_sum = 0._wp; pbwr3_sum = 0._wp
                    #:for BIN in range(nb)
                        radius = q(${NADV + 1 + BIN*BUB_STRIDE}$)
                        velocity = q(${NADV + 2 + BIN*BUB_STRIDE}$)
                        cubed = radius**3
                        #:if EARLY_BUBBLE
                            ! f_cpbw_KM with the invariant activation checks supplied by the host.
                            radius_b = real(radius, wp); r0_b = real(R0(${BIN + 1}$), wp)
                            velocity_b = real(velocity, wp)
                            ca_b = real(Ca, wp); web_b = real(Web, wp); re_b = real(Re_inv, wp)
                            ratio_b = (r0_b/radius_b)**(3._wp*real(gam, wp))
                            pbw = real(ca_b*ratio_b - ca_b + real(Eu, wp), wp)
                            if (bubble_web_active) then
                                pbw = real(real(pbw, wp) + (2._wp/(web_b*r0_b))*ratio_b, wp)
                                pbw = real(real(pbw, wp) - 2._wp/(radius_b*web_b), wp)
                            end if
                            if (bubble_re_active) pbw = real(real(pbw, wp) - 4._wp*re_b*velocity_b/radius_b, wp)
                        #:else
                            pbw = f_cpbw_KM(R0(${BIN + 1}$), radius, velocity, &
                                            & ${'0._wp' if polytropic else 'q(' + str(NADV + 3 + BIN*BUB_STRIDE) + ')'}$)
                        #:endif
                        #:if not adv_n
                            number_density = number_density + cubed*weight(${BIN + 1}$)
                        #:endif
                        r3_sum = r3_sum + cubed*weight(${BIN + 1}$)
                        r3v2_sum = r3v2_sum + cubed*velocity**2*weight(${BIN + 1}$)
                        pbwr3_sum = pbwr3_sum + pbw*cubed*weight(${BIN + 1}$)
                    #:endfor
                    #:if not adv_n
                        number_density = (3._wp/(4._wp*pi))*q(${NADV}$)/number_density
                    #:endif
                #:endif
                if (q(${NADV}$) < small_alf .or. r3_sum < small_alf) then
                    p_tilde = pressure - q(${NADV}$)*pressure
                else
                    p_tilde = pressure - q(${NADV}$)*(pressure - pbwr3_sum/r3_sum - rho*r3v2_sum/r3_sum)
                end if

            end subroutine s_hllc_bubble_state
        #:endif

        #:for DIR, AXIS in enumerate('xyz'[:ND], 1)
            #:if COMPRESS
                subroutine s_fused_thinc_${AXIS}$(alpha_center, &
                    #:if MTHINC
                        & normal_hat, normal_position, normal_index, &
                    #:else
                        & alpha_lower, alpha_upper, &
                    #:endif
                    & epsilon_ic, beta_ic, alpha_left, alpha_right, compressed)

                    $:GPU_ROUTINE(function_name='s_fused_thinc_${AXIS}$', parallelism='[seq]')

                    real(wp), value, intent(in) :: alpha_center, epsilon_ic, beta_ic
                    real(wp), intent(out)       :: alpha_left, alpha_right
                    logical, intent(out)        :: compressed
                    #:if MTHINC
                        real(wp), intent(in)       :: normal_hat(*), normal_position(*)
                        integer, value, intent(in) :: normal_index
                        real(wp)                   :: nh1, nh2, nh3, position
                    #:else
                        real(wp), value, intent(in) :: alpha_lower, alpha_upper
                        real(wp)                    :: alpha_min, alpha_span, profile_sign, profile, offset, tanh_beta
                    #:endif

                    compressed = .false.
                    if (.not. (alpha_center >= epsilon_ic .and. alpha_center <= 1._wp - epsilon_ic)) return
                    #:if MTHINC
                        nh1 = normal_hat(3*normal_index - 2)
                        nh2 = normal_hat(3*normal_index - 1)
                        nh3 = normal_hat(3*normal_index)
                        if (.not. (nh1*nh1 + nh2*nh2 + nh3*nh3 > 0.5_wp)) return
                        position = normal_position(normal_index)
                        alpha_left = f_mthinc_face_average(nh1, nh2, nh3, position, beta_ic, ${DIR}$, -0.5_wp, ${ND}$)
                        alpha_right = f_mthinc_face_average(nh1, nh2, nh3, position, beta_ic, ${DIR}$, 0.5_wp, ${ND}$)
                    #:else
                        if (.not. ((alpha_upper - alpha_center)*(alpha_center - alpha_lower) > moncon_cutoff)) return
                        profile_sign = merge(1._wp, -1._wp, alpha_upper - alpha_lower > 0._wp)
                        alpha_min = min(alpha_lower, alpha_upper)
                        alpha_span = max(alpha_lower, alpha_upper) - alpha_min
                        profile = (alpha_center - alpha_min + sgm_eps)/(alpha_span + sgm_eps)
                        profile = exp(profile_sign*beta_ic*(2._wp*profile - 1._wp))
                        tanh_beta = tanh(beta_ic)
                        offset = (profile/cosh(beta_ic) - 1._wp)/tanh_beta
                        alpha_left = alpha_min + 0.5_wp*alpha_span*(1._wp + profile_sign*offset)
                        alpha_right = alpha_min + 0.5_wp*alpha_span*(1._wp + profile_sign*(tanh_beta + offset)/(1._wp &
                            & + offset*tanh_beta))
                    #:endif
                    alpha_left = min(max(alpha_left, epsilon_ic), 1._wp - epsilon_ic)
                    alpha_right = min(max(alpha_right, epsilon_ic), 1._wp - epsilon_ic)
                    compressed = .true.

                end subroutine s_fused_thinc_${AXIS}$
            #:endif
            #:set GEO = cyl_coord and DIR > 1
            #:set FLUX_VALUE = lambda eq: f'flux({eq})'
            subroutine s_hllc_flux_${AXIS}$(q_left, q_right, flux, source, &
                                            & ${VISC_RESULTS}$${BUB_RESULTS}$ ${'geometric_flux, ' if GEO else ''}$${PHYS_GROUPS[0]}$, &
                                            & ${PHYS_GROUPS[1]}$, &
                                            & ${', '.join(PHYS_GROUPS[2:] + ([MOMENT_ARGS] if bubbles_euler and qbmm else []))}$)

                $:GPU_ROUTINE(function_name='s_hllc_flux_${AXIS}$', parallelism='[seq]')

                #:for GROUP in PHYS_GROUPS[:3]
                    real(wp), value, intent(in) :: ${GROUP}$
                #:endfor
                #:if mhd
                    real(wp), value, intent(in) :: magnetic_bx, cleaning_speed
                #:endif
                logical, value, intent(in) :: mixture_error
                #:if EARLY_BUBBLE
                    logical, value, intent(in) :: bubble_web_active, bubble_re_active
                #:endif
                #:if bubbles_euler
                    real(wp) :: nbub_l, nbub_r, p_tilde_l, p_tilde_r
                    #:if qbmm
                        real(wp), value, intent(in) :: ${MOMENT_ARGS}$
                        real(wp), intent(out)       :: bubble_transport(4)
                    #:endif
                #:endif
                real(wp), intent(in)  :: q_left(${NRIEMANN}$), q_right(${NRIEMANN}$)
                real(wp), intent(out) :: flux(${NE}$), source${f'({NF})' if hllc_riemann == 1 else ''}$
                #:if viscous
                    #:for GROUP in VISC_REAL_GROUPS
                        real(wp), value, intent(in) :: ${GROUP}$
                    #:endfor
                    logical, value, intent(in) :: shear_active, bulk_active
                    #:if hllc_non_newtonian
                        ! Carried with PHYS_GROUPS; the column computes the non-Newtonian Reynolds numbers.
                        #:for GROUP in NN_GROUPS
                            real(wp), value, intent(in) :: ${GROUP}$
                        #:endfor
                    #:endif
                    real(wp), intent(out) :: viscous_velocity(${ND}$)${', re_interface(2)' if not hllc_non_newtonian else ''}$
                    #:if not hllc_non_newtonian
                        real(wp) :: inverse_l, inverse_r, re_l, re_r
                    #:endif
                #:endif
                #:if surface_tension
                    real(wp), intent(out) :: viscous_velocity(${ND}$)
                #:endif
                real(wp) :: rho_l, rho_r, gamma_l, gamma_r, pi_l, pi_r, qv_l, qv_r, alpha_sum_l, alpha_sum_r
                real(wp) :: pressure_l, pressure_r, normal_l, normal_r, energy_l, energy_r, c_l, c_r
                real(wp) :: s_l, s_r, s_s, s_m, s_p, xi_l, xi_r, xi_l_m1, xi_r_m1, xi_m, xi_p
                real(wp) :: velocity_sq_l, velocity_sq_r, velocity_l(${NV}$), velocity_r(${NV}$)
                #:for SIDE in ['l', 'r']
                    real(wp) :: ${', '.join(f'density_{SIDE}{f}, alpha_{SIDE}{f}' for f in range(1, NF + 1))}$
                #:endfor
                #:if ELASTIC
                    real(wp) :: tau_l(6), tau_r(6)
                #:endif
                #:if hypoelasticity
                    real(wp) :: shear_l, shear_r, stress_power_l, stress_power_r
                #:endif
                #:if cont_damage
                    real(wp) :: damage_l, damage_r
                    real(wp) :: solid_partial_density_l, solid_partial_density_r
                #:endif
                #:if chemistry
                    #:set CHEM_MU = viscous and not hllc_non_newtonian
                    #:for SIDE in ['l', 'r']
                        #:if not hllc_scalar_thermo
                            real(wp) :: Ys_${SIDE}$(${CHEM_ROWS}$), Xs_${SIDE}$(${CHEM_ROWS}$), Cp_i${SIDE}$(${CHEM_ROWS}$), &
                                 & Gamma_i${SIDE}$(${CHEM_ROWS}$)
                        #:endif
                        real(wp) :: temperature_${SIDE}$, MW_${SIDE}$, R_gas_${SIDE}$, Cp_${SIDE}$, Cv_${SIDE}$, Gamm_${SIDE}$, &
                             & mu_${SIDE}$
                        #:if hllc_scalar_thermo and CHEM_MU
                            real(wp) :: log_T_${SIDE}$, wilke_phi_${SIDE}$, &
                                 & ${', '.join(f'mu_{SIDE}_{sp}' for sp in range(1, CHEM_ROWS + 1))}$
                        #:endif
                    #:endfor
                #:endif
                logical :: identical
                #:if GEO
                    real(wp), intent(out) :: geometric_flux
                #:endif

                #:if hllc_low_mach != 0
                    real(wp) :: zcoef, pressure_correction
                #:endif
                #:if SIXEQ
                    real(wp) :: xi_mp, xi_pp, energy_star, phasic_pressure, rho_star, normal_star
                    real(wp) :: interface_pressure, pi_over, exponent
                #:endif
                #:if SIXEQ or hllc_alt_soundspeed or any_state_dependent_eos
                    real(wp) :: sound_sum
                #:endif
                #:if SIXEQ and any_state_dependent_eos
                    real(wp) :: p_isen_l, p_isen_r, alpha_k_star, alpha_rho_k_star, e_k_star
                #:endif
                #:if any_state_dependent_eos
                    ! Per-phase EOS coefficients and frozen bulk modulus, per fluid and face side.
                    real(wp) &
                         & :: ${', '.join(f'{v}_{side}{f}' for side in ['l', 'r'] for f in range(1, NF + 1) for v in ['eos_gamma', 'eos_pi', 'eos_dpi', 'eos_dgamma', 'eos_rho', 'eos_blkmod'])}$
                #:endif

                #:if hllc_wave_speeds == 2 and hllc_riemann != 5
                    real(wp) :: rho_avg, c_avg, velocity_sq_avg, pressure_star, mach_l, mach_r
                #:endif

                #:if mhd
                    $:mhd_face_declarations()
                    $:mhd_face(ND, DIR, hllc_riemann, lambda side, row: f'q_{"left" if side == "l" else "right"}({row})', &
                               & FLUX_VALUE, 'gamma1', 'pi1', 'qv1', 'magnetic_bx', 'cleaning_speed', 'mixture_error')
                    source = 0._wp
                    #:if viscous
                        #:for V in range(1, ND + 1)
                            #:if hllc_riemann == 1
                                $:hll_interface_velocity('mf_sl', 'mf_sr', f'q_left({NF + V})', f'q_right({NF + V})', &
                                                         & f'viscous_velocity({V})')
                            #:else
                                viscous_velocity(${V}$) = 0.5_wp*(q_left(${NF + V}$) + q_right(${NF + V}$))
                            #:endif
                        #:endfor
                        #:if not hllc_non_newtonian
                            #:for S, STRESS in enumerate(['shear', 'bulk'], 1)
                                re_interface(${S}$) = 1._wp
                                if (${STRESS}$_active) then
                                    re_l = 1._wp/max(q_left(6)/re_${STRESS}$1, sgm_eps)
                                    re_r = 1._wp/max(q_right(6)/re_${STRESS}$1, sgm_eps)
                                    re_interface(${S}$) = 2._wp/(1._wp/re_l + 1._wp/re_r)
                                end if
                            #:endfor
                        #:endif
                    #:endif
                #:else
                    #:if HLL_FAMILY
                        real(wp) :: raw_alpha_l(${NF}$), raw_alpha_r(${NF}$)

                        #:for F in range(1, NF + 1)
                            raw_alpha_l(${F}$) = q_left(${ENERGY + F}$)
                            raw_alpha_r(${F}$) = q_right(${ENERGY + F}$)
                        #:endfor
                    #:endif

                    #:if EARLY_BUBBLE or EARLY_QBMM or EARLY_CHEM
                        identical = ${' .and. '.join(f'q_left({eq}) == q_right({eq})' for eq in range(1, NRIEMANN + 1))}$
                        #:if EARLY_QBMM
                            ! Equal QBMM pressure corrections also need equal moments.
                            #:for MOM in [1, 3, 4]
                                identical = identical .and. moment_l${MOM}$ == moment_r${MOM}$
                            #:endfor
                        #:endif
                    #:endif
                    #:for SIDE, Q in [('l', 'q_left'), ('r', 'q_right')]
                        #:for F in range(1, NF + 1)
                            density_${SIDE}$${F}$ = ${Q}$(${F}$)
                            alpha_${SIDE}$${F}$ = ${Q}$(${ENERGY + F}$)
                        #:endfor
                        velocity_${SIDE}$ = ${Q}$(${NF + 1}$:${NF + NV}$)
                        pressure_${SIDE}$ = ${Q}$(${ENERGY}$)
                        #:if hllc_mpp_lim
                            alpha_sum_${SIDE}$ = 0._wp
                            #:for F in range(1, NF + 1)
                                density_${SIDE}$${F}$ = max(0._wp, density_${SIDE}$${F}$)
                                alpha_${SIDE}$${F}$ = min(max(0._wp, alpha_${SIDE}$${F}$), 1._wp)
                                alpha_sum_${SIDE}$ = alpha_sum_${SIDE}$ + alpha_${SIDE}$${F}$
                            #:endfor
                            #:for F in range(1, NF + 1)
                                alpha_${SIDE}$${F}$ = alpha_${SIDE}$${F}$/max(alpha_sum_${SIDE}$, sgm_eps)
                            #:endfor
                        #:endif
                        rho_${SIDE}$ = 0._wp; gamma_${SIDE}$ = 0._wp; pi_${SIDE}$ = 0._wp; qv_${SIDE}$ = 0._wp
                        #:for F in range(1, NF + 1)
                            rho_${SIDE}$ = rho_${SIDE}$ + density_${SIDE}$${F}$
                            #:if any_state_dependent_eos
                                ! Mixture coefficients and frozen bulk modulus from each phase at its own density.
                                call s_phase_coefficients(density_${SIDE}$${F}$, alpha_${SIDE}$${F}$, ${F}$, &
                                                          & eos_rho_${SIDE}$${F}$, eos_gamma_${SIDE}$${F}$, eos_pi_${SIDE}$${F}$, &
                                                          & eos_dpi_${SIDE}$${F}$, eos_dgamma_${SIDE}$${F}$)
                                gamma_${SIDE}$ = gamma_${SIDE}$ + alpha_${SIDE}$${F}$*eos_gamma_${SIDE}$${F}$
                                pi_${SIDE}$ = pi_${SIDE}$ + alpha_${SIDE}$${F}$*eos_pi_${SIDE}$${F}$
                                eos_blkmod_${SIDE}$${F}$ = f_bulk_modulus(pressure_${SIDE}$, eos_gamma_${SIDE}$${F}$, &
                                    & eos_pi_${SIDE}$${F}$) - eos_rho_${SIDE}$${F}$*(eos_dpi_${SIDE}$${F}$ &
                                    & + pressure_${SIDE}$*eos_dgamma_${SIDE}$${F}$)/eos_gamma_${SIDE}$${F}$
                            #:else
                                gamma_${SIDE}$ = gamma_${SIDE}$ + alpha_${SIDE}$${F}$*gamma${F}$
                                pi_${SIDE}$ = pi_${SIDE}$ + alpha_${SIDE}$${F}$*pi${F}$
                            #:endif
                            #:if not qv_zero
                                qv_${SIDE}$ = qv_${SIDE}$ + density_${SIDE}$${F}$*qv${F}$
                            #:endif
                        #:endfor
                        #:if bubbles_euler
                            gamma_${SIDE}$ = gamma1
                            pi_${SIDE}$ = pi1
                            qv_${SIDE}$ = density_${SIDE}$1*qv1
                        #:endif
                        normal_${SIDE}$ = velocity_${SIDE}$(${DIR}$)
                        velocity_sq_${SIDE}$ = ${' + '.join(f'velocity_{SIDE}({v})**2' for v in range(1, NV + 1))}$
                        #:if bubbles_euler
                            ! Bubbly energy as in s_compute_energy, with the liquid internal energy diluted by (1 - alf).
                            energy_${SIDE}$ = gamma_${SIDE}$*pressure_${SIDE}$ + pi_${SIDE}$
                            energy_${SIDE}$ = energy_${SIDE}$*(1._wp - alpha_${SIDE}$${NF}$)
                            energy_${SIDE}$ = energy_${SIDE}$ + qv_${SIDE}$ + 5.e-1_wp*rho_${SIDE}$*velocity_sq_${SIDE}$
                        #:elif SIXEQ
                            energy_${SIDE}$ = gamma_${SIDE}$*pressure_${SIDE}$ + pi_${SIDE}$ &
                                & + 0.5_wp*rho_${SIDE}$*velocity_sq_${SIDE}$${f' + qv_{SIDE}' if not qv_zero else ''}$
                        #:else
                            energy_${SIDE}$ = gamma_${SIDE}$*pressure_${SIDE}$ + pi_${SIDE}$ + qv_${SIDE}$ &
                                & + 0.5_wp*rho_${SIDE}$*velocity_sq_${SIDE}$
                        #:endif
                        #:if chemistry
                            #:if not hllc_scalar_thermo
                                Ys_${SIDE}$ = ${Q}$(${NADV + 1}$:${NE}$)
                            #:endif
                            $:chemistry_face(SIDE, f'rho_{SIDE}', f'pressure_{SIDE}', f'velocity_sq_{SIDE}', f'gamma_{SIDE}', &
                                             & f'energy_{SIDE}', f'temperature_{SIDE}', f'MW_{SIDE}', &
                                             & GAMMA_GUARD='.not. identical' if EARLY_CHEM and HLL_FAMILY else None, &
                                             & FRACTION=(lambda sp: f'{Q}({NADV + sp})') if hllc_scalar_thermo else None, MU=f'mu_{SIDE}' if hllc_scalar_thermo and CHEM_MU else None)
                        #:endif
                        #:if viscous and not hllc_non_newtonian and SIDE == 'r'
                            #:for S, STRESS in enumerate(['shear', 'bulk'], 1)
                                re_interface(${S}$) = 1._wp
                                if (${STRESS}$_active) then
                                    inverse_l = 0._wp; inverse_r = 0._wp
                                    #:for F in range(1, NF + 1)
                                        if (re_${STRESS}$${F}$ > 0._wp) then
                                            inverse_l = ${'(1._wp - alpha_l' + str(F) + ')' if bubbles_euler else 'alpha_l' + str(F)}$/re_${STRESS}$${F}$ + inverse_l
                                            inverse_r = ${'(1._wp - alpha_r' + str(F) + ')' if bubbles_euler else 'alpha_r' + str(F)}$/re_${STRESS}$${F}$ + inverse_r
                                        end if
                                    #:endfor
                                    re_l = 1._wp/max(inverse_l, sgm_eps)
                                    re_r = 1._wp/max(inverse_r, sgm_eps)
                                    #:if chemistry and S == 1
                                        #:if not hllc_scalar_thermo
                                            call get_mixture_viscosity_mixavg(temperature_l, Ys_l, mu_l)
                                            call get_mixture_viscosity_mixavg(temperature_r, Ys_r, mu_r)
                                        #:endif
                                        ! The scalar closure supplies the Wilke mixture viscosity in mu_l/mu_r.
                                        re_l = 1._wp/mu_l; re_r = 1._wp/mu_r
                                    #:endif
                                    re_interface(${S}$) = 2._wp/(1._wp/re_l + 1._wp/re_r)
                                end if
                            #:endfor
                        #:endif
                        #:if bubbles_euler
                            #:if qbmm
                                ! Bubble pressure correction, inline as in s_hllc_riemann_solver.
                                nbub_${SIDE}$ = ${Q}$(${NADV + 1}$)
                                if (${Q}$(${NADV}$) < small_alf .or. moment_${SIDE}$1 < small_alf) then
                                    p_tilde_${SIDE}$ = pressure_${SIDE}$ - ${Q}$(${NADV}$)*pressure_${SIDE}$
                                else
                                    p_tilde_${SIDE}$ = pressure_${SIDE}$ - ${Q}$(${NADV}$)*(pressure_${SIDE}$ &
                                        & - moment_${SIDE}$4/moment_${SIDE}$1 - rho_${SIDE}$*moment_${SIDE}$3/moment_${SIDE}$1)
                                end if
                            #:else
                                call s_hllc_bubble_state(${Q}$, rho_${SIDE}$, pressure_${SIDE}$, nbub_${SIDE}$, &
                                                         & p_tilde_${SIDE}$${BUB_FLAGS}$)
                            #:endif
                        #:endif
                        #:if (not (bubbles_euler and qbmm) or EARLY_QBMM) and (not HLL_FAMILY or EARLY_CHEM) &
                              & and hllc_wave_speeds == 1 and hllc_low_mach == 0 and not hllc_alt_soundspeed and (not cyl_coord or SIXEQ) and SIDE == ('l' if EARLY_BUBBLE or EARLY_QBMM or EARLY_CHEM else 'r')
                            ! Identical faces take the physical flux of the shared state.
                            #:if not EARLY_BUBBLE and not EARLY_QBMM and not EARLY_CHEM
                                identical = ${' .and. '.join(f'{STATE("l", eq)} == {STATE("r", eq)}' for eq in range(1, NADV + 1))}$
                            #:endif
                            #:if bubbles_euler and not EARLY_BUBBLE and not EARLY_QBMM
                                identical = identical &
                                    & .and. ${' .and. '.join(f'q_left({eq}) == q_right({eq})' for eq in range(NADV + 1, NE + 1))}$
                            #:endif
                            #:if surface_tension
                                identical = identical .and. q_left(${NREC}$) == q_right(${NREC}$)
                            #:endif
                            #:if chemistry and not EARLY_CHEM
                                identical = identical &
                                    & .and. ${' .and. '.join(f'q_left({eq}) == q_right({eq})' for eq in range(NADV + 1, NE + 1))}$
                            #:endif
                            if (identical) then
                                #:if chemistry
                                    #:for EQ in range(NADV + 1, NE + 1)
                                        ${FLUX_VALUE(EQ)}$ = rho_l*q_left(${EQ}$)*normal_l
                                    #:endfor
                                #:endif
                                #:for F in range(1, NF + 1)
                                    ${FLUX_VALUE(F)}$ = ${STATE('l', F)}$*normal_l
                                    ${FLUX_VALUE(ENERGY + F)}$ = ${'0._wp' if hllc_riemann == 1 else STATE('l', ENERGY + F) + '*normal_l'}$
                                    #:if SIXEQ
                                        #:if any_state_dependent_eos
                                            ! Star ratios are one, so the phasic flux is e_phase*normal_l.
                                            call s_phase_internal_energy(pressure_l, alpha_l${F}$, density_l${F}$, ${F}$, e_k_star)
                                            ${FLUX_VALUE(NADV + F)}$ = e_k_star*normal_l
                                        #:else
                                            ${FLUX_VALUE(NADV + F)}$ = (alpha_l${F}$*(gamma${F}$*pressure_l + pi${F}$) &
                                                         & + density_l${F}$*qv${F}$)*normal_l
                                        #:endif
                                    #:endif
                                #:endfor
                                #:for V in range(1, NV + 1)
                                    ${FLUX_VALUE(NF + V)}$ = rho_l*normal_l*${STATE('l', NF + V)}$${' + p_tilde_l' if bubbles_euler and V == DIR else ' + pressure_l' if V == DIR else ''}$
                                #:endfor
                                ${FLUX_VALUE(ENERGY)}$ = normal_l*(energy_l + ${'p_tilde_l' if bubbles_euler else 'pressure_l'}$)
                                #:if hllc_riemann == 1
                                    #:for F in SOURCE_ROWS
                                        ${SOURCE(F)}$ = alpha_l${F}$
                                    #:endfor
                                #:else
                                    source = normal_l
                                #:endif
                                #:if EARLY_QBMM
                                    ${FLUX_VALUE(NADV + 1)}$ = nbub_l*normal_l
                                    bubble_transport(1) = (0.5_wp + sign(0.5_wp, normal_l))*nbub_l
                                    bubble_transport(2) = (0.5_wp - sign(0.5_wp, normal_l))*nbub_l
                                    bubble_transport(3) = normal_l
                                    bubble_transport(4) = normal_l
                                #:elif bubbles_euler
                                    #:for EQ in range(NADV + 1, NADV + BUB_ROWS + 1)
                                        ${FLUX_VALUE(EQ)}$ = nbub_l*q_left(${EQ}$)*normal_l
                                    #:endfor
                                    #:if adv_n
                                        ${FLUX_VALUE(NE)}$ = nbub_l*normal_l
                                    #:endif
                                #:endif
                                #:if surface_tension
                                    ${FLUX_VALUE(NE)}$ = q_left(${NREC}$)*normal_l
                                #:endif
                                #:if viscous or surface_tension
                                    viscous_velocity = velocity_l
                                #:endif
                                #:if GEO
                                    geometric_flux = ${FLUX_VALUE(NF + DIR)}$ - pressure_l
                                #:endif
                                return
                            end if
                        #:endif
                    #:endfor

                    #:if any_state_dependent_eos
                        ! Frozen-mixing sound speed of s_compute_speed_of_sound, Wood's law under alt_soundspeed.
                        #:for SIDE in ['l', 'r']
                            sound_sum = 0._wp
                            #:if hllc_alt_soundspeed
                                #:for F in range(1, NF + 1)
                                    sound_sum = sound_sum + alpha_${SIDE}$${F}$/eos_blkmod_${SIDE}$${F}$
                                #:endfor
                                c_${SIDE}$ = 1._wp/(rho_${SIDE}$*sound_sum)
                            #:else
                                #:for F in range(1, NF + 1)
                                    sound_sum = sound_sum + alpha_${SIDE}$${F}$*eos_blkmod_${SIDE}$${F}$
                                #:endfor
                                c_${SIDE}$ = sound_sum/rho_${SIDE}$
                            #:endif
                        #:endfor
                    #:elif chemistry
                        c_l = (1._wp + 1._wp/gamma_l)*pressure_l/rho_l
                        c_r = (1._wp + 1._wp/gamma_r)*pressure_r/rho_r
                    #:elif bubbles_euler
                        ! Operation order of s_compute_speed_of_sound.
                        c_l = f_bulk_modulus(pressure_l, gamma_l, pi_l)/rho_l
                        c_r = f_bulk_modulus(pressure_r, gamma_r, pi_r)/rho_r
                        #:if not (hllc_mpp_lim and num_fluids > 1)
                            c_l = c_l/(1._wp - q_left(${ENERGY + NF}$))
                            c_r = c_r/(1._wp - q_right(${ENERGY + NF}$))
                        #:endif
                    #:elif SIXEQ
                        #:for SIDE in ['l', 'r']
                            sound_sum = 0._wp
                            #:for F in range(1, NF + 1)
                                sound_sum = sound_sum + alpha_${SIDE}$${F}$*f_bulk_modulus(pressure_${SIDE}$, gamma${F}$, pi${F}$)
                            #:endfor
                            c_${SIDE}$ = sound_sum/rho_${SIDE}$
                        #:endfor
                    #:elif hllc_alt_soundspeed
                        sound_sum = ${STATE('l', ENERGY + 1)}$/f_bulk_modulus(pressure_l, gamma1, pi1)
                        sound_sum = sound_sum + ${STATE('l', ENERGY + 2)}$/f_bulk_modulus(pressure_l, gamma2, pi2)
                        c_l = 1._wp/(rho_l*sound_sum)
                        sound_sum = ${STATE('r', ENERGY + 1)}$/f_bulk_modulus(pressure_r, gamma1, pi1)
                        sound_sum = sound_sum + ${STATE('r', ENERGY + 2)}$/f_bulk_modulus(pressure_r, gamma2, pi2)
                        c_r = 1._wp/(rho_r*sound_sum)
                    #:else
                        c_l = f_bulk_modulus(pressure_l, gamma_l, pi_l)/rho_l
                        c_r = f_bulk_modulus(pressure_r, gamma_r, pi_r)/rho_r
                    #:endif
                    if (mixture_error .and. c_l < 0._wp) then
                        c_l = 100._wp*sgm_eps
                    else
                        c_l = sqrt(c_l)
                    end if
                    if (mixture_error .and. c_r < 0._wp) then
                        c_r = 100._wp*sgm_eps
                    else
                        c_r = sqrt(c_r)
                    end if
                    #:if hllc_wave_speeds == 2 and hllc_riemann != 5
                        rho_avg = 0.5_wp*(rho_l + rho_r)
                        #:if bubbles_euler
                            c_avg = (1._wp/(0.5_wp*(gamma_l + gamma_r)) + 1._wp)*(pressure_r + pi_r/(gamma_r + 1._wp)) &
                                     & /(rho_avg*(1._wp - alpha_r1))
                        #:elif SIXEQ
                            c_avg = sound_sum/rho_avg
                        #:elif hllc_alt_soundspeed
                            c_avg = 1._wp/(rho_avg*sound_sum)
                        #:else
                            velocity_sq_avg = 0._wp
                            #:for V in range(1, NV + 1)
                                velocity_sq_avg = velocity_sq_avg + (0.5_wp*(${STATE('l', NF + V)}$ + ${STATE('r', NF + V)}$))**2
                            #:endfor
                            c_avg = (0.5_wp*((energy_l + pressure_l)/rho_l + (energy_r + pressure_r)/rho_r) &
                                     & - 0.5_wp*velocity_sq_avg${' - 0.5_wp*(qv_l + qv_r)/rho_avg' if not qv_zero else ''}$) &
                                     & /(0.5_wp*(gamma_l + gamma_r))
                        #:endif
                        if (mixture_error .and. c_avg < 0._wp) then
                            c_avg = 100._wp*sgm_eps
                        else
                            c_avg = sqrt(c_avg)
                        end if
                    #:endif
                    #:if ELASTIC
                        tau_l = 0._wp; tau_r = 0._wp
                        #:for ROW in range(1, STRESS_ROWS + 1)
                            tau_l(${ROW}$) = q_left(${NADV + ROW}$)
                            tau_r(${ROW}$) = q_right(${NADV + ROW}$)
                        #:endfor
                        #:if hypoelasticity
                            #:if cont_damage
                                damage_l = q_left(${NADV + STRESS_ROWS + 1}$)
                                damage_r = q_right(${NADV + STRESS_ROWS + 1}$)
                                $:hypoelastic_interface_energy(lambda side, f: f'raw_alpha_{side.lower()}({f})', lambda side, &
                                                               & row: f'tau_{side.lower()}({row})', &
                                                               & lambda side: f'shear_{side.lower()}', &
                                                               & lambda side: f'energy_{side.lower()}', lambda side: f'damage_{side.lower()}')
                            #:else
                                $:hypoelastic_interface_energy(lambda side, f: f'raw_alpha_{side.lower()}({f})', lambda side, &
                                                               & row: f'tau_{side.lower()}({row})', &
                                                               & lambda side: f'shear_{side.lower()}', lambda side: f'energy_{side.lower()}')
                            #:endif
                        #:endif
                    #:endif
                    #:if hllc_low_mach != 0
                        zcoef = min(1._wp, max(sqrt(velocity_sq_l)/c_l, sqrt(velocity_sq_r)/c_r))
                    #:endif
                    #:if hllc_low_mach == 2
                        normal_l = 0.5_wp*(${STATE('l', NF + DIR)}$ + ${STATE('r', NF + DIR)}$ + zcoef*(${STATE('l', NF + DIR)}$ &
                                           & - ${STATE('r', NF + DIR)}$))
                        normal_r = 0.5_wp*(${STATE('l', NF + DIR)}$ + ${STATE('r', NF + DIR)}$ + zcoef*(${STATE('r', NF + DIR)}$ &
                                           & - ${STATE('l', NF + DIR)}$))
                        ${STATE('l', NF + DIR)}$ = normal_l; ${STATE('r', NF + DIR)}$ = normal_r
                    #:endif
                    #:if hllc_riemann == 5
                        s_r = max(sqrt(velocity_sq_l), sqrt(velocity_sq_r)) + max(c_l, c_r)
                        s_l = -s_r
                    #:elif hllc_wave_speeds == 2
                        pressure_star = 0.5_wp*(pressure_l + pressure_r + rho_avg*c_avg*(normal_l - normal_r))
                        mach_l = max(1._wp, &
                                     & sqrt(1._wp + ((0.5_wp + gamma_l)/(1._wp + gamma_l))*(pressure_star - pressure_l) &
                                     & /(pressure_l + pi_l/(1._wp + gamma_l))))
                        mach_r = max(1._wp, &
                                     & sqrt(1._wp + ((0.5_wp + gamma_r)/(1._wp + gamma_r))*(pressure_star - pressure_r) &
                                     & /(pressure_r + pi_r/(1._wp + gamma_r))))
                        s_l = normal_l - c_l*mach_l
                        s_r = normal_r + c_r*mach_r
                    #:elif hypoelasticity
                        ! Elastic signal speed of Rodriguez et al. JCP (2019), floored like
                        ! f_elastic_signal_speed in the classic solver.
                        c_l = sqrt(max(verysmall, c_l*c_l + (4._wp*shear_l/3._wp + tau_l(${STRESS_INDEX(DIR, DIR)}$))/rho_l))
                        c_r = sqrt(max(verysmall, c_r*c_r + (4._wp*shear_r/3._wp + tau_r(${STRESS_INDEX(DIR, DIR)}$))/rho_r))
                        s_l = min(normal_l - c_l, normal_r - c_r)
                        s_r = max(normal_r + c_r, normal_l + c_l)
                    #:else
                        s_l = min(normal_l - c_l, normal_r - c_r)
                        s_r = max(normal_r + c_r, normal_l + c_l)
                    #:endif
                    #:if HLL_FAMILY
                        s_m = min(0._wp, s_l); s_p = max(0._wp, s_r)
                        #:for F in range(1, NF + 1)
                            ${FLUX_VALUE(F)}$ = (s_m*${STATE('r', F)}$*normal_r - s_p*${STATE('l', F)}$*normal_l &
                                         & + s_m*s_p*(${STATE('l', F)}$ - ${STATE('r', F)}$))/(s_m - s_p)
                            #:if hllc_riemann == 5
                                ${FLUX_VALUE(ENERGY + F)}$ = (s_m*raw_alpha_r(${F}$)*normal_r - s_p*raw_alpha_l(${F}$)*normal_l &
                                             & + s_m*s_p*(raw_alpha_l(${F}$) - raw_alpha_r(${F}$)))/(s_m - s_p)
                            #:else
                                ${FLUX_VALUE(ENERGY + F)}$ = (raw_alpha_l(${F}$) - raw_alpha_r(${F}$))*s_m*s_p/(s_m - s_p)
                            #:endif
                        #:endfor
                        #:for V in range(1, NV + 1)
                            ${FLUX_VALUE(NF + V)}$ = (s_m*(rho_r*normal_r*${STATE('r', NF + V)}$${' + pressure_r' if V == DIR else ''}$${' - tau_r(' + str(STRESS_INDEX(DIR, V)) + ')' if hypoelasticity else ''}$) - s_p*(rho_l*normal_l*${STATE('l', NF + V)}$${' + pressure_l' if V == DIR else ''}$${' - tau_l(' + str(STRESS_INDEX(DIR, V)) + ')' if hypoelasticity else ''}$) + s_m*s_p*(rho_l*${STATE('l', NF + V)}$ - rho_r*${STATE('r', NF + V)}$))/(s_m - s_p)
                        #:endfor
                        #:if hypoelasticity
                            stress_power_l = ${' + '.join(f'tau_l({STRESS_INDEX(DIR, v)})*velocity_l({v})' for v in range(1, NV + 1))}$
                            stress_power_r = ${' + '.join(f'tau_r({STRESS_INDEX(DIR, v)})*velocity_r({v})' for v in range(1, NV + 1))}$
                            ${FLUX_VALUE(ENERGY)}$ = (s_m*(normal_r*(energy_r + pressure_r) - stress_power_r) &
                                         & - s_p*(normal_l*(energy_l + pressure_l) - stress_power_l) + s_m*s_p*(energy_l &
                                         & - energy_r))/(s_m - s_p)
                            #:for ROW in range(1, STRESS_ROWS + 1)
                                ${FLUX_VALUE(NADV + ROW)}$ = (s_m*(rho_r*normal_r*tau_r(${ROW}$)) &
                                             & - s_p*(rho_l*normal_l*tau_l(${ROW}$)) + s_m*s_p*(rho_l*tau_l(${ROW}$) &
                                             & - rho_r*tau_r(${ROW}$)))/(s_m - s_p)
                            #:endfor
                            #:if cont_damage
                                ! U_D = m_s*D advects with the HLL mass-flux pattern.
                                solid_partial_density_l = 0._wp
                                solid_partial_density_r = 0._wp
                                #:for F in range(1, NF + 1)
                                    if (Gs_rs(${F}$) > verysmall) then
                                        solid_partial_density_l = solid_partial_density_l + density_l${F}$
                                        solid_partial_density_r = solid_partial_density_r + density_r${F}$
                                    end if
                                #:endfor
                                ${FLUX_VALUE(NADV + STRESS_ROWS + 1)}$ = (s_m*(solid_partial_density_r*normal_r*damage_r) &
                                             & - s_p*(solid_partial_density_l*normal_l*damage_l) &
                                             & + s_m*s_p*(solid_partial_density_l*damage_l - solid_partial_density_r*damage_r)) &
                                             & /(s_m - s_p)
                            #:endif
                        #:else
                            ${FLUX_VALUE(ENERGY)}$ = (s_m*normal_r*(energy_r + pressure_r) - s_p*normal_l*(energy_l + pressure_l) &
                                         & + s_m*s_p*(energy_l - energy_r))/(s_m - s_p)
                        #:endif
                        #:if hllc_riemann == 5
                            source = (s_m*normal_r - s_p*normal_l)/(s_m - s_p)
                        #:else
                            #:for F in SOURCE_ROWS
                                ${SOURCE(F)}$ = (s_m*raw_alpha_r(${F}$) - s_p*raw_alpha_l(${F}$))/(s_m - s_p)
                            #:endfor
                        #:endif
                    #:else
                        #:if hllc_wave_speeds == 2
                            s_s = 0.5_wp*(normal_l + normal_r + (pressure_l - pressure_r)/(rho_avg*c_avg))
                        #:else
                            s_s = (pressure_r - pressure_l + rho_l*normal_l*(s_l - normal_l) - rho_r*normal_r*(s_r - normal_r)) &
                                   & /(rho_l*(s_l - normal_l) - rho_r*(s_r - normal_r))
                        #:endif
                        s_m = min(0._wp, s_l); s_p = max(0._wp, s_r)
                        xi_l = (s_l - normal_l)/min(s_l - s_s, -sgm_eps)
                        xi_r = (s_r - normal_r)/max(s_r - s_s, sgm_eps)
                        xi_l_m1 = (s_s - normal_l)/min(s_l - s_s, -sgm_eps)
                        xi_r_m1 = (s_s - normal_r)/max(s_r - s_s, sgm_eps)
                        xi_m = 0.5_wp + sign(0.5_wp, s_s)
                        xi_p = 0.5_wp - sign(0.5_wp, s_s)
                        #:if SIXEQ
                            xi_mp = -min(0._wp, sign(1._wp, s_l))
                            xi_pp = max(0._wp, sign(1._wp, s_r))
                            energy_star = xi_m*(energy_l + xi_mp*(xi_l*(energy_l + (s_s - normal_l)*(rho_l*s_s + pressure_l/(s_l &
                                                & - normal_l))) - energy_l)) + xi_p*(energy_r + xi_pp*(xi_r*(energy_r + (s_s &
                                                & - normal_r)*(rho_r*s_s + pressure_r/(s_r - normal_r))) - energy_r))
                            interface_pressure = xi_m*(pressure_l + xi_mp*(rho_l*(s_l - normal_l)*(s_s - normal_l))) &
                                                       & + xi_p*(pressure_r + xi_pp*(rho_r*(s_r - normal_r)*(s_s - normal_r)))
                            rho_star = xi_m*(rho_l*(xi_mp*xi_l + 1._wp - xi_mp)) + xi_p*(rho_r*(xi_pp*xi_r + 1._wp - xi_pp))
                            normal_star = normal_l*(1._wp - xi_mp) + xi_mp*normal_r + xi_mp*xi_pp*(s_s - normal_r)
                            #:for F in range(1, NF + 1)
                                ${FLUX_VALUE(F)}$ = xi_m*density_l${F}$*(normal_l + s_m*xi_l_m1) + xi_p*density_r${F}$*(normal_r &
                                             & + s_p*xi_r_m1)
                                ${FLUX_VALUE(ENERGY + F)}$ = xi_m*alpha_l${F}$*s_s + xi_p*alpha_r${F}$*s_s
                                #:if any_state_dependent_eos
                                    ! Six-equation internal-energy flux from the phasic isentrope
                                    ! pressure at the star density, as in s_hllc_riemann_solver.
                                    call s_phase_pressure_on_isentrope(pressure_l, density_l${F}$/max(alpha_l${F}$, sgm_eps), &
                                                                       & xi_l, ${F}$, p_isen_l)
                                    call s_phase_pressure_on_isentrope(pressure_r, density_r${F}$/max(alpha_r${F}$, sgm_eps), &
                                                                       & xi_r, ${F}$, p_isen_r)
                                    phasic_pressure = xi_m*(xi_mp*(p_isen_l - pressure_l) + pressure_l) + xi_p*(xi_pp*(p_isen_r &
                                                            & - pressure_r) + pressure_r)
                                    alpha_k_star = xi_m*alpha_l${F}$ + xi_p*alpha_r${F}$
                                    alpha_rho_k_star = (xi_m*density_l${F}$ + xi_p*density_r${F}$)*(1._wp + xi_m*xi_mp*(xi_l &
                                                        & - 1._wp) + xi_p*xi_pp*(xi_r - 1._wp))
                                    call s_phase_internal_energy(phasic_pressure, alpha_k_star, alpha_rho_k_star, ${F}$, e_k_star)
                                    ${FLUX_VALUE(NADV + F)}$ = e_k_star*normal_star
                                #:else
                                    pi_over = pi${F}$/(1._wp + gamma${F}$)
                                    exponent = 1._wp/gamma${F}$ + 1._wp
                                    phasic_pressure = xi_m*(xi_mp*((pressure_l + pi_over)*xi_l**exponent - pi_over - pressure_l) &
                                                            & + pressure_l) + xi_p*(xi_pp*((pressure_r + pi_over)*xi_r**exponent &
                                                            & - pi_over - pressure_r) + pressure_r)
                                    ${FLUX_VALUE(NADV + F)}$ = ((xi_m*alpha_l${F}$ + xi_p*alpha_r${F}$) &
                                                 & *(gamma${F}$*phasic_pressure + pi${F}$) + (xi_m*density_l${F}$ &
                                                 & + xi_p*density_r${F}$)*qv${F}$)*normal_star
                                #:endif
                            #:endfor
                            #:for V in range(1, NV + 1)
                                ${FLUX_VALUE(NF + V)}$ = rho_star*normal_star*${'normal_star + interface_pressure' if V == DIR else f'(xi_m*velocity_l({V}) + xi_p*velocity_r({V}))'}$
                            #:endfor
                            ${FLUX_VALUE(ENERGY)}$ = (energy_star + interface_pressure)*normal_star
                            source = xi_m*(normal_l + (s_s*(xi_mp*xi_l_m1 + 1._wp) - normal_l)) + xi_p*(normal_r &
                                           & + (s_s*(xi_pp*xi_r_m1 + 1._wp) - normal_r))
                        #:elif bubbles_euler
                            #:for EQ in list(range(1, NF + 1)) + list(range(ENERGY + 1, NADV + 1))
                                ${FLUX_VALUE(EQ)}$ = xi_m*${STATE('l', EQ)}$*(normal_l + s_m*xi_l_m1) &
                                             & + xi_p*${STATE('r', EQ)}$*(normal_r + s_p*xi_r_m1)
                            #:endfor
                            #:for V in range(1, NV + 1)
                                ${FLUX_VALUE(NF + V)}$ = xi_m*(rho_l*(normal_l*velocity_l(${V}$) &
                                             & + s_m*(xi_l*${'s_s' if V == DIR else 'velocity_l(' + str(V) + ')'}$ &
                                             & - velocity_l(${V}$)))${' + p_tilde_l' if V == DIR else ''}$) &
                                             & + xi_p*(rho_r*(normal_r*velocity_r(${V}$) &
                                             & + s_p*(xi_r*${'s_s' if V == DIR else 'velocity_r(' + str(V) + ')'}$ &
                                             & - velocity_r(${V}$)))${' + p_tilde_r' if V == DIR else ''}$)
                            #:endfor
                            ${FLUX_VALUE(ENERGY)}$ = xi_m*(normal_l*(energy_l + p_tilde_l) + s_m*(xi_l*(energy_l + (s_s &
                                         & - normal_l)*(rho_l*s_s + p_tilde_l/(s_l - normal_l))) - energy_l)) &
                                         & + xi_p*(normal_r*(energy_r + p_tilde_r) + s_p*(xi_r*(energy_r + (s_s - normal_r) &
                                         & *(rho_r*s_s + p_tilde_r/(s_r - normal_r))) - energy_r))
                            source = xi_m*(normal_l + s_m*xi_l_m1) + xi_p*(normal_r + s_p*xi_r_m1)
                            #:if qbmm
                                bubble_transport(1) = xi_m*nbub_l
                                bubble_transport(2) = xi_p*nbub_r
                                bubble_transport(3) = normal_l + s_m*xi_l_m1
                                bubble_transport(4) = normal_r + s_p*xi_r_m1
                            #:else
                                #:for EQ in range(NADV + 1, NADV + BUB_ROWS + 1)
                                    ${FLUX_VALUE(EQ)}$ = xi_m*nbub_l*q_left(${EQ}$)*(normal_l + s_m*xi_l_m1) &
                                                 & + xi_p*nbub_r*q_right(${EQ}$)*(normal_r + s_p*xi_r_m1)
                                #:endfor
                            #:endif
                            #:for EQ in ([NADV + 1] if qbmm else []) + ([NE] if adv_n else [])
                                ${FLUX_VALUE(EQ)}$ = xi_m*nbub_l*(normal_l + s_m*xi_l_m1) + xi_p*nbub_r*(normal_r + s_p*xi_r_m1)
                            #:endfor
                        #:else
                            #:for EQ in list(range(1, NF + 1)) + list(range(ENERGY + 1, NADV + 1))
                                ${FLUX_VALUE(EQ)}$ = xi_m*${STATE('l', EQ)}$*(normal_l + s_m*xi_l_m1) &
                                             & + xi_p*${STATE('r', EQ)}$*(normal_r + s_p*xi_r_m1)
                            #:endfor
                            #:for V in range(1, NV + 1)
                                #:if V == DIR
                                    ${FLUX_VALUE(NF + V)}$ = xi_m*(rho_l*(normal_l*${STATE('l', NF + V)}$ + s_m*s_l*xi_l_m1) &
                                                 & + pressure_l) + xi_p*(rho_r*(normal_r*${STATE('r', NF + V)}$ + s_p*s_r*xi_r_m1) &
                                                 & + pressure_r)
                                #:else
                                    ! Transverse momentum star flux, unfactored as in s_hllc_riemann_solver.
                                    ${FLUX_VALUE(NF + V)}$ = xi_m*(rho_l*(normal_l*${STATE('l', NF + V)}$ &
                                                 & + s_m*${STATE('l', NF + V)}$*xi_l_m1)) &
                                                 & + xi_p*(rho_r*(normal_r*${STATE('r', NF + V)}$ &
                                                 & + s_p*${STATE('r', NF + V)}$*xi_r_m1))
                                #:endif
                            #:endfor
                            ${FLUX_VALUE(ENERGY)}$ = xi_m*(normal_l*(energy_l + pressure_l) + s_m*(energy_l*xi_l_m1 + xi_l*(s_s &
                                         & - normal_l)*(rho_l*s_s + pressure_l/(s_l - normal_l)))) + xi_p*(normal_r*(energy_r &
                                         & + pressure_r) + s_p*(energy_r*xi_r_m1 + xi_r*(s_s - normal_r)*(rho_r*s_s &
                                         & + pressure_r/(s_r - normal_r))))
                            source = xi_m*(normal_l + s_m*xi_l_m1) + xi_p*(normal_r + s_p*xi_r_m1)
                        #:endif
                    #:endif

                    #:if chemistry
                        #:for EQ in range(NADV + 1, NE + 1)
                            #:if HLL_FAMILY
                                ${FLUX_VALUE(EQ)}$ = (s_m*q_right(${EQ}$)*rho_r*normal_r - s_p*q_left(${EQ}$)*rho_l*normal_l &
                                             & + s_m*s_p*(q_left(${EQ}$)*rho_l - q_right(${EQ}$)*rho_r))/(s_m - s_p)
                            #:else
                                ${FLUX_VALUE(EQ)}$ = xi_m*rho_l*q_left(${EQ}$)*(normal_l + s_m*xi_l_m1) &
                                             & + xi_p*rho_r*q_right(${EQ}$)*(normal_r + s_p*xi_r_m1)
                            #:endif
                        #:endfor
                    #:endif
                    #:if surface_tension
                        #:if SIXEQ
                            ${FLUX_VALUE(NE)}$ = (xi_m*q_left(${NREC}$) + xi_p*q_right(${NREC}$))*s_s
                        #:else
                            ${FLUX_VALUE(NE)}$ = xi_m*q_left(${NREC}$)*(normal_l + s_m*xi_l_m1) + xi_p*q_right(${NREC}$) &
                                         & *(normal_r + s_p*xi_r_m1)
                        #:endif
                    #:endif
                    #:if viscous or surface_tension
                        #:for V in range(1, NV + 1)
                            #:if hllc_riemann == 1
                                $:hll_interface_velocity('s_l', 's_r', f'velocity_l({V})', f'velocity_r({V})', &
                                                         & f'viscous_velocity({V})')
                            #:elif hllc_riemann == 5
                                viscous_velocity(${V}$) = 0.5_wp*(velocity_l(${V}$) + velocity_r(${V}$))
                            #:elif V != DIR
                                viscous_velocity(${V}$) = xi_m*velocity_l(${V}$) + xi_p*velocity_r(${V}$)
                            #:elif SIXEQ
                                viscous_velocity(${V}$) = xi_m*(velocity_l(${V}$) + (s_s*(xi_mp*xi_l_m1 + 1._wp) &
                                                 & - velocity_l(${V}$))) + xi_p*(velocity_r(${V}$) + (s_s*(xi_pp*xi_r_m1 + 1._wp) &
                                                 & - velocity_r(${V}$)))
                            #:else
                                viscous_velocity(${V}$) = xi_m*(velocity_l(${V}$) + s_m*xi_l_m1) + xi_p*(velocity_r(${V}$) &
                                                 & + s_p*xi_r_m1)
                            #:endif
                        #:endfor
                    #:endif

                    #:if hllc_low_mach == 1
                        #:if HLL_FAMILY
                            pressure_correction = -(s_p - s_m)*(rho_l + rho_r)/8._wp*(zcoef - 1._wp)
                        #:else
                            pressure_correction = rho_l*rho_r*(s_l - normal_l)*(s_r - normal_r)*(normal_r - normal_l)/(rho_r*(s_r &
                                                               & - normal_r) - rho_l*(s_l - normal_l))*(zcoef - 1._wp)
                        #:endif
                        pressure_correction = (s_m/s_l)*(s_p/s_r)*pressure_correction
                        #:if HLL_FAMILY
                            #:for V in range(1, NV + 1)
                                ${FLUX_VALUE(NF + V)}$ = ${FLUX_VALUE(NF + V)}$ + pressure_correction*(${STATE('r', NF + V)}$ &
                                             & - ${STATE('l', NF + V)}$)
                            #:endfor
                            ${FLUX_VALUE(ENERGY)}$ = ${FLUX_VALUE(ENERGY)}$ + pressure_correction*(velocity_sq_r - velocity_sq_l) &
                                         & /2._wp
                        #:else
                            ${FLUX_VALUE(NF + DIR)}$ = ${FLUX_VALUE(NF + DIR)}$ + pressure_correction
                            ${FLUX_VALUE(ENERGY)}$ = ${FLUX_VALUE(ENERGY)}$ + pressure_correction*s_s
                        #:endif
                    #:endif

                    #:if GEO
                        #:if SIXEQ
                            geometric_flux = ${FLUX_VALUE(NF + DIR)}$ - interface_pressure
                        #:elif HLL_FAMILY
                            geometric_flux = ${FLUX_VALUE(NF + DIR)}$ - (s_m*pressure_r - s_p*pressure_l)/(s_m - s_p)
                            #:if hypoelasticity and DIR == 2
                                geometric_flux = geometric_flux + (s_m*tau_r(${STRESS_ROWS}$) - s_p*tau_l(${STRESS_ROWS}$))/(s_m &
                                                                   & - s_p)
                            #:endif
                        #:else
                            geometric_flux = xi_m*(rho_l*(normal_l*normal_l + s_m*(xi_l*s_s - normal_l))) &
                                                   & + xi_p*(rho_r*(normal_r*normal_r + s_p*(xi_r*s_s - normal_r)))
                        #:endif
                    #:endif
                #:endif

            end subroutine s_hllc_flux_${AXIS}$
        #:endfor

        #:if not SIXEQ
            #:for BOUNDARY_ONLY in [True, False]
                #:set RK_NAME = 's_fused_rhs_openacc_boundary_rk' if BOUNDARY_ONLY else 's_fused_rhs_openacc_update_rk'
                subroutine ${RK_NAME}$(q_cons, q_stor, q_output, q_next_prim, q_stage_prim, rhs, ${RK_SCALARS}$)

                    type(scalar_field), intent(inout) :: q_cons(${NE}$), q_stor(${NE}$), q_output(${NE}$), q_next_prim(${NE}$)
                    type(scalar_field), intent(in)    :: q_stage_prim(${NE}$), rhs(${NE}$)
                    logical, intent(in)               :: save_rk_state
                    real(wp), intent(in)              :: rk_q_coef, rk_stor_coef, rk_rhs_coef, rk_den_coef
                    integer                           :: side, boundary, normal_beg, normal_end

                    #:set AXIS = 'xyz'[ND - 1]
                    #:set EXTENT = ['m', 'n', 'p'][ND - 1]
                    #:if BOUNDARY_ONLY
                        do side = 1, 2
                            if (side == 1) then
                                boundary = bc_${AXIS}$%beg
                                normal_beg = 0; normal_end = buff_size + 1
                            else
                                boundary = bc_${AXIS}$%end
                                normal_beg = ${EXTENT}$ - buff_size - 1; normal_end = ${EXTENT}$
                            end if
                            if (boundary > BC_CHAR_SLIP_WALL .or. boundary < BC_CHAR_SUP_OUTFLOW) cycle
                        #:else
                            normal_beg = 0; normal_end = ${EXTENT}$
                        #:endif
                        call s_fused_rhs_openacc_boundary_launch(&
                        #:for FIELD in ['q_stage_prim', 'rhs', 'q_cons', 'q_stor', 'q_output', 'q_next_prim']
                            #:for EQ in range(1, NE + 1)
                                & ${FIELD}$(${EQ}$)%sf${'' if FIELD == 'rhs' else '(idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg)'}$, &
                                            & &
                            #:endfor
                        #:endfor
                        & m, n, p, idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg, size(q_cons(1)%sf, 1), size(q_cons(1)%sf, 2), &
                                           & size(rhs(1)%sf, 1), size(rhs(1)%sf, 2), normal_beg, normal_end, &
                                           & ${', '.join(f'gammas({f})' for f in range(1, NF + 1))}$, &
                                           & ${', '.join(f'pi_infs({f})' for f in range(1, NF + 1))}$, &
                                           & ${', '.join(f'qvs({f})' for f in range(1, NF + 1))}$${', Bx0' if mhd else ''}$, &
                                           & ${RK_SCALARS}$)
                        #:if BOUNDARY_ONLY
                        end do
                    #:endif

                end subroutine ${RK_NAME}$
            #:endfor

            subroutine s_fused_rhs_openacc_boundary_launch(&
                #:for FIELD in ['q', 'rhs'] + RK_FIELDS
                    & ${', '.join(f'{FIELD}{eq}' for eq in range(1, NE + 1))}$, &
                #:endfor
                & m_in, n_in, p_in, i_lo, j_lo, k_lo, ni, nj, ri, rj, normal_beg, normal_end, &
                    & ${', '.join(PHYS_GROUPS[:3])}$${', magnetic_bx' if mhd else ''}$, ${RK_SCALARS}$)

                #:for FIELD in ['q', 'rhs'] + RK_FIELDS
                    real(stp), &
                         & intent(${'in' if FIELD in ['q', 'rhs', 'cons'] else 'inout'}$) &
                         & :: ${', '.join(f'{FIELD}{eq}(*)' for eq in range(1, NE + 1))}$
                #:endfor
                integer, value, intent(in) :: m_in, n_in, p_in, i_lo, j_lo, k_lo, ni, nj, ri, rj, normal_beg, normal_end
                #:for GROUP in PHYS_GROUPS[:3]
                    real(wp), value, intent(in) :: ${GROUP}$
                #:endfor
                #:if mhd
                    real(wp), value, intent(in) :: magnetic_bx
                #:endif
                logical, value, intent(in)  :: save_rk_state
                real(wp), value, intent(in) :: rk_q_coef, rk_stor_coef, rk_rhs_coef, rk_den_coef
                real(wp)                    :: ${', '.join(RK_PRIVATE)}$
                integer                     :: x, y, z, q_index, rhs_index
                integer                     :: xb, xe, yb, ye, zb, ze

                ! The boundary RK update stays inside the active-box window; the band bounds
                ! follow the swept axis, so they clip to that axis's window.
                #:set BAND_BOX = ['ab_x', 'ab_y', 'ab_z'][ND - 1]
                xb = 0; xe = m_in
                yb = 0; ye = n_in
                zb = normal_beg; ze = normal_end
#ifdef MFC_ACTIVE_BOX
                if (ab_active) then
                    xb = max(xb, ab_x%beg); xe = min(xe, ab_x%end)
                    yb = max(yb, ab_y%beg); ye = min(ye, ab_y%end)
                    zb = max(zb, ${BAND_BOX}$%beg); ze = min(ze, ${BAND_BOX}$%end)
                end if
#endif

                ! Empty windows skip the launch, since a zero-gang kernel deadlocks inside captured CUDA graphs.
                if (zb <= ze .and. yb <= ye .and. xb <= xe) then
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[x, y, z, q_index, rhs_index, ' + ', '.join(RK_PRIVATE) + ']')
                    do z = ${'zb, ze' if ND == 3 else '0, p_in'}$
                        do y = ${'zb, ze' if ND == 2 else 'yb, ye'}$
                            do x = ${'zb, ze' if ND == 1 else 'xb, xe'}$
                                q_index = 1 + x - i_lo + ni*((y - j_lo) + nj*(z - k_lo))
                                rhs_index = 1 + x + ri*(y + rj*z)
                                #:for EQ in range(1, NE + 1)
                                    value${EQ}$ = real(rhs${EQ}$(rhs_index), wp)
                                #:endfor
                                $:hllc_rk_update((lambda field, eq: f'{"q" if field == "stage" else field}{eq}(q_index)'))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if

            end subroutine s_fused_rhs_openacc_boundary_launch
        #:endif

        #:for WENO in [False, True]
            #:set NAME = 's_weno_fused_rhs_openacc_kernel' if WENO else 's_fused_rhs_openacc'
            #:set STENCIL = WENO and WENO_RADIUS > 0
            #:set COEFS = STENCIL and not MUSCL
            #:set COMP = STENCIL and COMPRESS
            #:set COMP_GROUPS = (['epsilon_ic, beta_ic'] + (['normal_hat, normal_position, normal_ni, normal_nj'] &
                                 & if MTHINC else [])) if COMP else []
            #:set NFACE = NREC if surface_tension or bubbles_euler or hypoelasticity or chemistry or mhd else NADV
            #:set EXTRA_GROUPS = ['epsilon_muscl'] if STENCIL and MUSCL else ((['d_l, d_r', 'coef_lo, epsilon_weno, teno_cutoff'] if TENO7 else ['poly_l, poly_r, d_l, d_r, beta_coef', &
                & 'coef_lo, coef_n, uniform, epsilon_weno' + (', teno_cutoff' if teno else '')]) if COEFS else ([] if WENO else ['qL, qR, qli, qlj, qln']))
            #:set BASE_ARG_GROUPS = [', '.join(f'q{eq}' for eq in range(1, NE + 1)), ', '.join(f'rhs{eq}' for eq in range(1, NE + 1)), &
                                          & 'm_in, n_in, p_in, buff_in, i_lo, j_lo, k_lo', &
                                          & 'ni, nj, ri, rj, ds, ds_lo, normal_beg, normal_end'] + (['launch_beg, launch_end'] if WENO else []) + (['radius_cc, radius_cb, radius_lo'] if cyl_coord else []) + PHYS_GROUPS + QBMM_GROUPS + GRAD_GROUPS + EXTRA_GROUPS + COMP_GROUPS
            subroutine ${NAME}$(q_prim, rhs, direction, &
                #:if IN_COLUMN_STRESS
                    & ${', '.join(GRAD_FIELDS)}$, &
                #:endif
                #:if WENO
                    #:if STENCIL and MUSCL
                        & epsilon_muscl, &
                    #:elif COEFS
                        #:if TENO7
                            & d_l, d_r, coef_lo, epsilon_weno, teno_cutoff, &
                        #:else
                            & poly_l, poly_r, d_l, d_r, beta_coef, coef_lo, coef_n, uniform, &
                                & epsilon_weno${', teno_cutoff' if teno else ''}$, &
                        #:endif
                    #:endif
                    & q_cons, q_stor, q_output, q_next_prim, fuse_final_rk, ${RK_SCALARS}$, overlap_region)
                #:else
                    & qL, qR)
                #:endif

                type(scalar_field), intent(in)    :: q_prim(${NE}$)
                type(scalar_field), intent(inout) :: rhs(${NE}$)
                integer, intent(in)               :: direction
                #:if IN_COLUMN_STRESS
                    type(scalar_field), intent(in) :: ${', '.join(f'{field}(:)' for field in GRAD_FIELDS)}$
                #:endif
                #:if WENO
                    type(scalar_field), intent(inout) :: q_cons(${NE}$), q_stor(${NE}$), q_output(${NE}$), q_next_prim(${NE}$)
                    integer, optional, intent(in)     :: overlap_region
                    logical, intent(in)               :: fuse_final_rk, save_rk_state
                    real(wp), intent(in)              :: rk_q_coef, rk_stor_coef, rk_rhs_coef, rk_den_coef
                    #:if STENCIL and MUSCL
                        real(wp), intent(in) :: epsilon_muscl
                    #:elif COEFS
                        #:if TENO7
                            integer, intent(in)          :: coef_lo
                            real(wp), target, intent(in) :: d_l(*), d_r(*)
                            real(wp), intent(in)         :: epsilon_weno, teno_cutoff
                        #:else
                            integer, intent(in)          :: coef_lo, coef_n
                            real(wp), target, intent(in) :: poly_l(*), poly_r(*), d_l(*), d_r(*), beta_coef(*)
                            logical, intent(in)          :: uniform
                            real(wp), intent(in)         :: epsilon_weno${', teno_cutoff' if teno else ''}$
                        #:endif
                    #:endif
                #:else
                    real(wp), intent(in) :: qL(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:)
                    real(wp), intent(in) :: qR(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:)
                #:endif
                integer :: normal_beg, normal_end
                #:if WENO
                    integer :: region, part, part_count, first, last, full_beg, full_end, launch_beg, launch_end
                    region = 0
                    if (present(overlap_region)) region = overlap_region
                    if (direction /= 1) region = 0
                    call s_fused_rhs_openacc_overlap_bounds(first, last)
                    part_count = merge(2, 1, region == 2)
                #:endif

#ifdef MFC_ACTIVE_BOX
                ! An empty active-box window freezes the whole rank, so no sweep launches.
                if (ab_active .and. (ab_x%beg > ab_x%end .or. ab_y%beg > ab_y%end .or. ab_z%beg > ab_z%end)) return
#endif

                select case (direction)
                    #:for DIR, AXIS in enumerate('xyz'[:ND], 1)
                        #:set EXTENT = ['m', 'n', 'p'][DIR - 1]
                    case (${DIR}$)
                        normal_beg = 0; normal_end = ${EXTENT}$
                        if (bc_${AXIS}$%beg <= BC_CHAR_SLIP_WALL .and. bc_${AXIS}$%beg >= BC_CHAR_SUP_OUTFLOW) &
                            & normal_beg = buff_size + 2
                        if (bc_${AXIS}$%end <= BC_CHAR_SLIP_WALL .and. bc_${AXIS}$%end >= BC_CHAR_SUP_OUTFLOW) &
                            & normal_end = ${EXTENT}$ - buff_size - 2
                        #:if hllc_cbc
                            ! The swept-axis divergence and RK update run inside the active-box window.
                            ! The hllc_cbc columns guard the divergence loop on the normal range.
#ifdef MFC_ACTIVE_BOX
                            if (ab_active) then
                                normal_beg = max(normal_beg, ab_${AXIS}$%beg)
                                normal_end = min(normal_end, ab_${AXIS}$%end)
                            end if
#endif
                        #:endif
                        #:if WENO
                            full_beg = normal_beg; full_end = normal_end
                            do part = 1, part_count
                                normal_beg = full_beg; normal_end = full_end
                                ! CBC consumes face fluxes from columns outside the native RHS range.
                                launch_beg = 0; launch_end = ${EXTENT}$
                                if (region == 1) then
                                    normal_beg = max(full_beg, first); normal_end = min(full_end, last)
                                    launch_beg = first; launch_end = last
                                else if (region == 2) then
                                    if (part == 1) then
                                        normal_end = min(full_end, first - 1)
                                        launch_end = first - 1
                                    else
                                        normal_beg = max(full_beg, last + 1)
                                        launch_beg = last + 1
                                    end if
                                end if
                                #:if hllc_cbc
                                    ! Clip the column launch range to the active-box window, keeping the
                                    ! CBC bands of the swept axis, whose columns publish face fluxes.
#ifdef MFC_ACTIVE_BOX
                                    if (ab_active) then
                                        block
                                            integer :: lo, hi, l0, l1
                                            l0 = launch_beg; l1 = launch_end
                                            lo = ${EXTENT}$ + 1; hi = -1
                                            if (ab_${AXIS}$%end >= l0 .and. ab_${AXIS}$%beg <= l1) then
                                                lo = min(lo, max(ab_${AXIS}$%beg, l0))
                                                hi = max(hi, min(ab_${AXIS}$%end, l1))
                                            end if
                                            if (l0 == 0 .and. bc_${AXIS}$%beg <= BC_CHAR_SLIP_WALL &
                                                & .and. bc_${AXIS}$%beg >= BC_CHAR_SUP_OUTFLOW) then
                                                lo = min(lo, 0); hi = max(hi, buff_size + 2)
                                            end if
                                            if (l1 == ${EXTENT}$ .and. bc_${AXIS}$%end <= BC_CHAR_SLIP_WALL &
                                                & .and. bc_${AXIS}$%end >= BC_CHAR_SUP_OUTFLOW) then
                                                lo = min(lo, ${EXTENT}$ - buff_size - 2); hi = max(hi, ${EXTENT}$)
                                            end if
                                            launch_beg = lo; launch_end = hi
                                        end block
                                    end if
#endif
                                #:endif
                                if (launch_beg > launch_end) cycle
                            #:endif
                            call ${NAME}$_launch_${AXIS}$(&
                            #:for EQ in range(1, NE + 1)
                                & q_prim(${EQ}$)%sf(idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg), &
                            #:endfor
                            #:for EQ in range(1, NE + 1)
                                & rhs(${EQ}$)%sf, &
                            #:endfor
                            & m, n, p, buff_size, idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg, &
                            & size(q_prim(1)%sf, 1), size(q_prim(1)%sf, 2), size(rhs(1)%sf, 1), size(rhs(1)%sf, 2), &
                            & d${AXIS}$(idwbuff(${DIR}$)%beg), idwbuff(${DIR}$)%beg, normal_beg, normal_end, &
                            #:if WENO
                                & launch_beg, launch_end, &
                            #:endif
                            #:if cyl_coord
                                & y_cc(idwbuff(2)%beg), y_cb(idwbuff(2)%beg), idwbuff(2)%beg, &
                            #:endif
                            & ${', '.join(f'gammas({f})' for f in range(1, NF + 1))}$, &
                            & ${', '.join(f'pi_infs({f})' for f in range(1, NF + 1))}$, &
                            & ${', '.join(f'qvs({f})' for f in range(1, NF + 1))}$, mixture_err${', Bx0, hyper_cleaning_speed' if mhd else ''}$${VISC_HOST_ARGS}$${NN_HOST_ARGS}$${BUB_HOST_FLAGS}$${',' if (bubbles_euler and qbmm) or IN_COLUMN_STRESS or not WENO or STENCIL else ''}$ &
                            #:if bubbles_euler and qbmm
                                #:for MOM in [1, 3, 4]
                                    & mom_sp(${MOM}$)%sf(idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg), &
                                #:endfor
                                & size(mom_sp(1)%sf, 1), size(mom_sp(1)%sf, &
                                       & 2)${',' if IN_COLUMN_STRESS or not WENO or STENCIL else ''}$ &
                            #:endif
                            #:if IN_COLUMN_STRESS
                                #:for FIELD in GRAD_FIELDS
                                    #:for V in range(1, ND + 1)
                                        & ${FIELD}$(${NF + V}$)%sf, &
                                    #:endfor
                                #:endfor
                                & size(grad_x_vf(${NF + 1}$)%sf, 1), size(grad_x_vf(${NF + 1}$)%sf, &
                                       & 2)${',' if not WENO or STENCIL else ''}$ &
                            #:endif
                            #:if WENO
                                #:if STENCIL and MUSCL
                                    & epsilon_muscl &
                                #:elif COEFS
                                    #:if TENO7
                                        & d_l, d_r, coef_lo, epsilon_weno, teno_cutoff &
                                    #:else
                                        & poly_l, poly_r, d_l, d_r, beta_coef, coef_lo, coef_n, uniform, &
                                            & epsilon_weno${', teno_cutoff' if teno else ''}$ &
                                    #:endif
                                #:endif
                                #:if COMP
                                    &, ic_eps, ic_beta &
                                    #:if MTHINC
                                        &, mthinc_nhat(1, idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg), &
                                    & mthinc_d(idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg), &
                                    & size(mthinc_d, 1), size(mthinc_d, 2) &
                                    #:endif
                                #:endif
                                #:if DIR == ND and not SIXEQ and (not viscous or IN_COLUMN_STRESS) &
                                                                  & and not surface_tension and not bubbles_euler and not hypoelasticity and not chemistry and not relativity and (not mhd or not hyper_cleaning) &
                                                                  & and not CYL_STRESS
                                    #:for FIELD in ['q_cons', 'q_stor', 'q_output', 'q_next_prim']
                                        #:for EQ in range(1, NE + 1)
                                            &, ${FIELD}$(${EQ}$)%sf(idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg) &
                                        #:endfor
                                    #:endfor
                                    &, ${RK_SCALARS}$, fuse_final_rk &
                                #:endif
                                &)
                            #:else
                                & qL, qR, size(qL, 1), size(qL, 2), size(qL, 1)*size(qL, 2)*size(qL, 3))
                            #:endif
                            #:if WENO
                            end do
                        #:endif
                    #:endfor
                end select

            end subroutine ${NAME}$

            #:for DIR, AXIS in enumerate('xyz'[:ND], 1)
                #:set FINAL = WENO and DIR == ND and not SIXEQ and (not viscous or IN_COLUMN_STRESS) &
                        & and not surface_tension and not bubbles_euler and not hypoelasticity and not chemistry and not relativity and (not mhd or not hyper_cleaning) &
                        & and not CYL_STRESS
                #! grid_geometry 3 runs no fused RK update (f_use_fused_rhs_openacc_final_rk),
                #! so CYL_STRESS columns omit the RK epilogue.
                #:set KEEP_VALUES = FINAL or IN_COLUMN_STRESS
                #:set UPDATE_PRIVATE = RK_PRIVATE if FINAL else ([f'value{eq}' for eq in range(1, &
                        & NE + 1)] if KEEP_VALUES else [])
                #:set GEO = cyl_coord and DIR > 1
                #:set CBC_GEO = GEO and DIR == 2 and hllc_cbc
                #:set ARG_GROUPS = BASE_ARG_GROUPS + (RK_GROUPS + ['fuse_final_rk'] if FINAL else [])
                #:set EXTENT = ['m_in', 'n_in', 'p_in'][DIR - 1]
                #:set WIDE_TRANSVERSE = SIXEQ and WENO and WENO_RADIUS == 1 and not MUSCL
                #:set CELLS = MFC_FUSION_X_PADDED_CELLS if DIR == 1 else (MFC_FUSION_TRANSVERSE_WIDE_CELLS if WIDE_TRANSVERSE else MFC_FUSION_TRANSVERSE_CELLS)
                #:if chemistry
                    ! Size the column so every species row fits the staging buffer.
                    #:set VISC_ROWS = ND + 1 if IN_COLUMN_STRESS else 0
                    #:set SHARED_SLOPE = (3 if STENCIL else 2)*NREC + VISC_ROWS
                    #:set SHARED_FIXED = (2*WENO_RADIUS + 4 if STENCIL else 2)*NREC + VISC_ROWS
                    #:set CELLS = min(CELLS, (MFC_FUSION_STAGING_BYTES//8 - SHARED_FIXED)//SHARED_SLOPE)
                    #:assert CELLS >= 1
                    #:set CHEM_LANE_BYTES = 8*(SHARED_SLOPE*CELLS + SHARED_FIXED)
                #:endif
                #:set BUBBLE_COMPACT = (bubbles_euler or hypoelasticity) &
                                            & and DIR > 1 and STENCIL and (NREC*(3*MFC_FUSION_TRANSVERSE_CELLS + 4 &
                                            & + 2*WENO_RADIUS) + ((ND + 1)*(MFC_FUSION_TRANSVERSE_CELLS + 1) &
                                            & if IN_COLUMN_STRESS else 0))*8*MFC_FUSION_TRANSVERSE_LANES > MFC_FUSION_STAGING_BYTES
                #! 3D weno_Re_flux sweeps use compact lanes to hold the gradient columns.
                #:set GRAD_COMPACT = GRAD_WENO_FLUX and STENCIL and not TENO7 and ND > 2
                #:set LANES = 1 if DIR == 1 else (MFC_FUSION_COMPACT_LANES if BUBBLE_COMPACT or GRAD_COMPACT else (MFC_FUSION_TRANSVERSE_WIDE_LANES if WIDE_TRANSVERSE else MFC_FUSION_TRANSVERSE_LANES))
                #:if chemistry and DIR > 1
                    #:set LANES = max(lanes for lanes in MFC_FUSION_CHEM_LANES if lanes*CHEM_LANE_BYTES <= MFC_FUSION_STAGING_BYTES)
                #:endif
                #:set LANE_DIM = f'{LANES}, '
                #:set LANE_INDEX = 'lane, '
                #:set TANGENT = 'p_in' if DIR < 3 else 'n_in'
                #:set CROSS = 'n_in' if DIR == 1 else f'm_in/{LANES}'
                #:set BLOCK_LOOPS = [('block_x', CROSS), ('block_n', f'{EXTENT}/{CELLS}')]
                #:set GANG_COUNT = f' num_gangs(({TANGENT} + 1)*({CROSS} + 1)*({EXTENT}/{CELLS} + 1))' if bubbles_euler and qbmm else ''
                #:set XYZ = [['normal', 'block_x', 'tangent'], ['x', 'normal', 'tangent'], ['x', 'tangent', &
                            & 'normal']][DIR - 1]
                #:set INDEX = ', '.join(XYZ)
                #:set QINDEX = f'1 + {XYZ[0]} - i_lo + ni*(({XYZ[1]} - j_lo) + nj*({XYZ[2]} - k_lo))'
                #:set MOM_INDEX = QINDEX.replace('ni', 'mi').replace('nj', 'mj')
                #:set MOM_VALUES = ', '.join('real(moment' + str(i) + '(' + MOM_INDEX + (' + ' + ['1', 'mi', &
                                                 & 'mi*mj'][DIR - 1] if side == 'r' else '') + '), wp)' for side in ['l', &
                                                 & 'r'] for i in [1, 3, 4])
                #:set QLINDEX = QINDEX.replace('ni', 'qli').replace('nj', 'qlj')
                #:set GRAD_INDEX = QINDEX.replace('ni', 'gi').replace('nj', 'gj')
                #:set NORMAL_INDEX = QINDEX.replace('ni', 'normal_ni').replace('nj', 'normal_nj')
                #:set RINDEX = f'1 + {XYZ[0]} + ri*({XYZ[1]} + rj*{XYZ[2]})'
                #:set FACE_FLUX = lambda eq: f'left{eq}'
                #:set FACE_SOURCE = lambda f: f'right{f}'
                #:set QLSTRIDE = ['1', 'qli', 'qli*qlj'][DIR - 1]
                #:set GRAD_WENO = GRAD_WENO_FLUX and STENCIL and not TENO7
                #! GREAD reads the staged face average under GRAD_WENO and the host face gradients otherwise.
                #:set SWEEP_STRIDE = ['1', 'ni', 'ni*nj'][DIR - 1]
                #:set GREAD = (lambda a, v: 'gleft_{a}{v}('.format(a=a, &
                                   & v=v) + LANE_INDEX + 'state)') if GRAD_WENO else (lambda a, &
                                   & v: 'real(grad_{a}{v}(grad_index), wp)'.format(a=a, v=v))
                #:for COLUMN in [False, True]
                    #:set SUFFIX = 'column' if COLUMN else 'launch'
                    #! Cell-centred alpha pair across the face, from the column cells in the stencil
                    #! variants and from the global state otherwise.
                    #:set NN_ALPHAS = [('column{e}('.format(e=ENERGY + f) + LANE_INDEX &
                                            & + 'state + {r}) + column{e}('.format(e=ENERGY + f, &
                                            & r=WENO_RADIUS) + LANE_INDEX + 'state + {r})'.format(r=WENO_RADIUS + 1)) &
                                            & if STENCIL and SUFFIX == 'column' &
                                            & else ('q{e}(q_index) + q{e}(q_index + {s})'.format(e=ENERGY + f, &
                                            & s=SWEEP_STRIDE) &
                                            & if WENO else 'qR(ql_index + qln*{m}) + qL(ql_index + {t} + qln*{m})'.format(e=ENERGY + f, &
                                            & m=ENERGY + f - 1, t=QLSTRIDE)) for f in range(1, NF + 1)]
                    #! Velocity pair for the cylindrical hoop terms, from the reconstructed faces
                    #! under weno_Re_flux and from the cells flanking the face otherwise.
                    #:set CYL_VEL = [('q_left({e}) + q_right({e})'.format(e=NF + v) if hllc_weno_re_flux &
                                          & else ('column{e}('.format(e=NF + v) + LANE_INDEX + 'state + {r}) &
                                          & + column{e}('.format(e=NF + v, r=WENO_RADIUS) + LANE_INDEX &
                                          & + 'state + {r})'.format(r=WENO_RADIUS + 1)) if STENCIL and SUFFIX == 'column' &
                                          & else ('q{e}(q_index) + q{e}(q_index + {s})'.format(e=NF + v, &
                                          & s=SWEEP_STRIDE) if WENO else 'qR(ql_index + qln*{m}) + qL(ql_index + {t} &
                                          & + qln*{m})'.format(m=NF + v - 1, t=QLSTRIDE))) for v in (2, 3)]
                    subroutine ${NAME}$_${SUFFIX}$_${AXIS}$(&
                        #:for GROUP in ARG_GROUPS[:-1]
                            & ${GROUP}$, &
                        #:endfor
                        & ${ARG_GROUPS[-1]}$${', block_n, block_x, tangent' if COLUMN else ''}$)

                        #:if COLUMN
                            $:GPU_ROUTINE(function_name='${NAME}$_column_${AXIS}$', parallelism='[vector]')
                            integer, value, intent(in) :: block_n, block_x, tangent
                        #:else
                            integer :: block_n, block_x, tangent
                            integer :: cb_lo, cb_hi, tg_lo, tg_hi
                            ! Active-box cover column spans, box index and swept-axis CBC band bounds.
                            integer :: ab_b, ab_nlo, ab_nhi, ab_clo, ab_chi, ab_tlo, ab_thi, ab_band1_hi, ab_band2_lo
                            logical :: ab_inb, ab_band1, ab_band2
                        #:endif
                        #:for EQ in range(1, NE + 1)
                            real(stp), intent(in)    :: q${EQ}$(*)
                            real(stp), intent(inout) :: rhs${EQ}$(*)
                        #:endfor
                        integer, value, intent(in) :: m_in, n_in, p_in, buff_in, i_lo, j_lo, k_lo, ni, nj, ri, rj, ds_lo
                        integer, value, intent(in) :: normal_beg, normal_end
                        #:if WENO
                            integer, value, intent(in) :: launch_beg, launch_end
                        #:endif
                        real(wp), intent(in) :: ds(*)
                        #:if cyl_coord
                            real(wp), intent(in)       :: radius_cc(*)
                            real(wp), intent(in)       :: radius_cb(*)
                            integer, value, intent(in) :: radius_lo
                        #:endif
                        #:for GROUP in PHYS_GROUPS[:3]
                            real(wp), value, intent(in) :: ${GROUP}$
                        #:endfor
                        #:if mhd
                            real(wp), value, intent(in) :: magnetic_bx, cleaning_speed
                        #:endif
                        logical, value, intent(in) :: mixture_error
                        #:if EARLY_BUBBLE
                            logical, value, intent(in) :: bubble_web_active, bubble_re_active
                        #:endif
                        #:if viscous
                            #:for GROUP in VISC_REAL_GROUPS
                                real(wp), value, intent(in) :: ${GROUP}$
                            #:endfor
                            logical, value, intent(in) :: shear_active, bulk_active
                            #:if hllc_non_newtonian
                                #:for GROUP in NN_GROUPS
                                    real(wp), value, intent(in) :: ${GROUP}$
                                #:endfor
                            #:endif
                        #:endif
                        #:if bubbles_euler and qbmm
                            real(stp), intent(in)      :: moment1(*), moment3(*), moment4(*)
                            integer, value, intent(in) :: mi, mj
                        #:endif
                        #:if IN_COLUMN_STRESS
                            real(stp), intent(in)      :: ${GRAD_GROUPS[0].replace(',', '(*),')}$ (*)
                            integer, value, intent(in) :: gi, gj
                        #:endif
                        #:if COMP
                            real(wp), value, intent(in) :: epsilon_ic, beta_ic
                            #:if MTHINC
                                real(wp), intent(in)       :: normal_hat(*), normal_position(*)
                                integer, value, intent(in) :: normal_ni, normal_nj
                            #:endif
                        #:endif
                        #:if FINAL
                            #:for GROUP in RK_ARRAY_GROUPS
                                real(stp), intent(inout) :: ${GROUP.replace(',', '(*),')}$ (*)
                            #:endfor
                            logical, value, intent(in)  :: fuse_final_rk, save_rk_state
                            real(wp), value, intent(in) :: rk_q_coef, rk_stor_coef, rk_rhs_coef, rk_den_coef
                        #:endif
                        #:if STENCIL and MUSCL
                            real(wp), value, intent(in) :: epsilon_muscl
                        #:elif COEFS
                            #:if TENO7
                                integer, value, intent(in)  :: coef_lo
                                real(wp), intent(in)        :: d_l(*), d_r(*)
                                real(wp), value, intent(in) :: epsilon_weno, teno_cutoff
                            #:else
                                integer, value, intent(in)  :: coef_lo, coef_n
                                real(wp), intent(in)        :: poly_l(*), poly_r(*), d_l(*), d_r(*), beta_coef(*)
                                logical, value, intent(in)  :: uniform
                                real(wp), value, intent(in) :: epsilon_weno${', teno_cutoff' if teno else ''}$
                            #:endif
                        #:elif not WENO
                            real(wp), intent(in)       :: qL(*), qR(*)
                            integer, value, intent(in) :: qli, qlj, qln
                        #:endif
                        #:if COLUMN
                            #:for EQ in range(1, NFACE + 1)
                                #:if STENCIL
                                    real(stp) :: column${EQ}$(${LANE_DIM}$${CELLS + 2*(WENO_RADIUS + 1)}$)
                                #:endif
                                real(wp) :: left${EQ}$(${LANE_DIM}$${CELLS + 1}$), right${EQ}$(${LANE_DIM}$${CELLS + 1}$)
                            #:endfor
                            #:if SIXEQ
                                #:for EQ in range(NFACE + 1 if surface_tension else NADV + 1, NE + 1)
                                    real(wp) :: left${EQ}$(${LANE_DIM}$${CELLS + 1}$)
                                #:endfor
                            #:endif
                            real(wp) :: q_left(${NRIEMANN}$), q_right(${NRIEMANN}$), flux(${NE}$), &
                                 & source${f'({NF})' if hllc_riemann == 1 else ''}$, inv_ds, value, state_left, state_right
                            #:if hllc_alt_soundspeed
                                real(wp) :: k_term
                            #:endif
                            #:if GEO
                                real(wp) :: geometric_flux, radial_coord, half_over_radius
                            #:endif
                            #:if bubbles_euler and qbmm
                                real(wp) :: bubble_transport(4), bubble_left, bubble_right
                            #:endif
                            integer :: lane, state, normal, x, q_index, rhs_index, ql_index
                            #:if IN_COLUMN_STRESS
                                #:for V in range(1, ND + 2)
                                    real(wp) :: viscous${V}$(${LANE_DIM}$${CELLS + 1}$)
                                #:endfor
                                real(wp) :: viscous_velocity(${ND}$), re_interface(2), viscous_flux(${ND + 1}$), divergence, tau
                                integer  :: grad_index
                                #:if hllc_non_newtonian
                                    real(wp) :: alpha_avg, gamma_dot, d_xx, d_yy, d_zz, d_xy, d_xz, d_yz
                                #:endif
                                #:if CYL_STRESS
                                    real(wp) :: r_eff, div_v_term_const, sigma1, sigma2, sigma3, avg_v2, avg_v3
                                #:endif
                                #:if GRAD_WENO
                                    #:for A in range(1, ND + 1)
                                        #:for V in range(1, NV + 1)
                                            real(stp) &
                                                 & :: gcolumn_${'xyz'[A - 1] + str(V)}$(${LANE_DIM}$${CELLS + 2*(WENO_RADIUS + 1)}$)
                                            real(wp) :: gleft_${'xyz'[A - 1] + str(V)}$(${LANE_DIM}$${CELLS + 1}$), &
                                                 & gright_${'xyz'[A - 1] + str(V)}$(${LANE_DIM}$${CELLS + 1}$)
                                        #:endfor
                                    #:endfor
                                #:endif
                            #:elif viscous or surface_tension
                                real(wp) :: viscous_velocity(${ND}$) &
                                     & ${', re_interface(2)' if viscous and not hllc_non_newtonian else ''}$
                            #:endif
                            #:if COMP
                                real(wp) :: compression_alpha, compression_rho1, compression_rho2
                                logical  :: compress_state
                            #:endif
                            #:if UPDATE_PRIVATE
                                real(wp) :: ${', '.join(UPDATE_PRIVATE)}$
                            #:endif
                        #:endif
                        #:if not COLUMN
                            ! Clip the transverse column ranges to the active-box window, so the
                            ! final-direction RK update stays inside it. Window bounds are lane-stride
                            ! aligned, and the per-cell m_in guard covers the partial edge block.
                            #:set CROSS_BOX = 'ab_y' if DIR == 1 else 'ab_x'
                            #:set TGT_BOX = 'ab_z' if DIR < 3 else 'ab_y'
                            cb_lo = 0; cb_hi = ${'n_in' if DIR == 1 else 'm_in'}$
                            tg_lo = 0; tg_hi = ${'p_in' if DIR < 3 else 'n_in'}$
#ifdef MFC_ACTIVE_BOX
                            if (ab_active) then
                                cb_lo = max(cb_lo, ${CROSS_BOX}$%beg); cb_hi = min(cb_hi, ${CROSS_BOX}$%end)
                                tg_lo = max(tg_lo, ${TGT_BOX}$%beg); tg_hi = min(tg_hi, ${TGT_BOX}$%end)
                            end if
                            #! The swept-axis CBC bands keep running under the active-box cover,
                            #! since their columns publish the face fluxes CBC consumes.
                            #:set AB_LO_BAND = ('bc_' + AXIS + '%beg <= BC_CHAR_SLIP_WALL .and. bc_' + AXIS &
                                                    & + '%beg >= BC_CHAR_SUP_OUTFLOW') if hllc_cbc else '.false.'
                            #:set AB_HI_BAND = ('bc_' + AXIS + '%end <= BC_CHAR_SLIP_WALL .and. bc_' + AXIS &
                                                    & + '%end >= BC_CHAR_SUP_OUTFLOW') if hllc_cbc else '.false.'
                            if (ab_active .and. (ab_nb > 1 .or. ab_inner > 0)) then
                                ab_band1 = ${AB_LO_BAND}$
                                ab_band2 = ${AB_HI_BAND}$
                                ab_band1_hi = buff_size + 2
                                ab_band2_lo = ${EXTENT}$ - buff_size - 2
                            end if
#endif
                            $:GPU_PARALLEL_LOOP(collapse=3, parallelism='[gang]', &
                                                & extraAccArgs=('vector_length(%d)' % (MFC_FUSION_VECTOR_X if DIR == 1 else (max(MFC_FUSION_VECTOR_X, &
                                                & MFC_FUSION_VECTOR_CHEM_PER_LANE*LANES) &
                                                & if chemistry else MFC_FUSION_VECTOR_COMPACT if BUBBLE_COMPACT else MFC_FUSION_VECTOR_TRANSVERSE))) + GANG_COUNT)
                            do tangent = tg_lo, tg_hi
                                #:for LOOP, BOUND in BLOCK_LOOPS
                                    do ${LOOP}$ = ${'launch_beg/' + str(CELLS) if WENO and LOOP == 'block_n' else ('cb_lo' + ('/' + str(LANES) if DIR > 1 else '') if LOOP == 'block_x' else '0')}$, &
                                                                        & ${'launch_end/' + str(CELLS) if WENO and LOOP == 'block_n' else ('cb_hi' + ('/' + str(LANES) if DIR > 1 else '') if LOOP == 'block_x' else BOUND)}$
                                    #:endfor
#ifdef MFC_ACTIVE_BOX
                                    ! Columns outside every active box and CBC band are skipped
                                    ! with their fused RK update; their cells stay frozen.
                                    if (ab_active .and. (ab_nb > 1 .or. ab_inner > 0)) then
                                        ab_nlo = block_n*${CELLS}$; ab_nhi = ab_nlo + ${CELLS - 1}$
                                        #:if DIR == 1
                                            ab_clo = block_x; ab_chi = block_x
                                        #:else
                                            ab_clo = block_x*${LANES}$; ab_chi = ab_clo + ${LANES - 1}$
                                        #:endif
                                        ab_tlo = tangent; ab_thi = tangent
                                        ab_inb = .false.
                                        do ab_b = 1, ab_nb
                                            if (ab_nlo <= ab_box(${2*DIR}$, ab_b) .and. ab_nhi >= ab_box(${2*DIR - 1}$, &
                                                & ab_b) .and. ab_clo <= ab_box(${4 if DIR == 1 else 2}$, &
                                                & ab_b) .and. ab_chi >= ab_box(${3 if DIR == 1 else 1}$, &
                                                & ab_b) .and. ab_tlo <= ab_box(${6 if DIR < 3 else 4}$, &
                                                & ab_b) .and. ab_thi >= ab_box(${5 if DIR < 3 else 3}$, ab_b)) then
                                                ab_inb = .true.
                                                exit
                                            end if
                                        end do
                                        if (.not. ab_inb .and. .not. (ab_band1 .and. ab_nlo <= ab_band1_hi) .and. .not. (ab_band2 &
                                            & .and. ab_nhi >= ab_band2_lo)) cycle
                                    end if
#endif
                                    call ${NAME}$_column_${AXIS}$(&
                                    #:for GROUP in ARG_GROUPS
                                        & ${GROUP}$, &
                                    #:endfor
                                    & block_n, block_x, tangent)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    #:endif
                    #:if COLUMN
                        #:if STENCIL
                            $:GPU_LOOP(collapse=2, parallelism='[vector]', private='[normal, x, q_index]')
                            do state = 1, ${CELLS + 2*(WENO_RADIUS + 1)}$
                                do lane = 1, ${LANES}$
                                    normal = block_n*${CELLS}$ - ${WENO_RADIUS + 2}$ + state
                                    x = block_x*${LANES}$ + lane - 1
                                    q_index = ${QINDEX}$
                                    if (normal <= ${EXTENT}$ + buff_in${' .and. x <= m_in' if DIR > 1 else ''}$) then
                                        #:for EQ in range(1, NFACE + 1)
                                            column${EQ}$(${LANE_INDEX}$state) = q${REC_FIELD(EQ)}$(q_index)
                                        #:endfor
                                        #:if GRAD_WENO
                                            #:for A in range(1, ND + 1)
                                                #:for V in range(1, NV + 1)
                                                    gcolumn_${'xyz'[A - 1] + str(V)}$(${LANE_INDEX}$state) &
                                                                                 & = grad_${'xyz'[A - 1] + str(V)}$(${GRAD_INDEX}$)
                                                #:endfor
                                            #:endfor
                                        #:endif
                                    end if
                                end do
                            end do
                            $:GPU_LOOP(collapse=2, parallelism='[vector]', &
                                       & private='[normal, x, q_index, state_left, state_right' &
                                       & + (', compression_alpha, compression_rho1, compression_rho2, compress_state' if COMP else '') + ']')
                            do state = 1, ${CELLS + 2}$
                                do lane = 1, ${LANES}$
                                    normal = block_n*${CELLS}$ - 2 + state
                                    x = block_x*${LANES}$ + lane - 1
                                    q_index = ${QINDEX}$
                                    if (normal <= ${EXTENT}$ + 1${' .and. x <= m_in' if DIR > 1 else ''}$) then
                                        #:for EQ in range(1, NFACE + 1)
                                            #:if surface_tension and EQ == ENERGY
                                                state_left = real(column${EQ}$(${LANE_INDEX}$state + ${WENO_RADIUS}$), wp)
                                                state_right = state_left
                                            #:elif MUSCL
                                                call s_fused_weno_state(column${EQ}$(1, 1), lane + ${LANES}$*state, ${LANES}$, &
                                                                        & epsilon_muscl, state_left, state_right)
                                            #:elif TENO7
                                                call s_fused_weno_state(column${EQ}$(1, 1), &
                                                                        & lane + ${LANES}$*(state + ${WENO_RADIUS - 1}$), normal, &
                                                                        & ${LANES}$, d_l, d_r, coef_lo, epsilon_weno, &
                                                                        & teno_cutoff, state_left, state_right)
                                            #:else
                                                call s_fused_weno_state(column${EQ}$(1, 1), &
                                                                        & lane + ${LANES}$*(state + ${WENO_RADIUS - 1}$), normal, &
                                                                        & ${LANES}$, poly_l, poly_r, d_l, d_r, beta_coef, &
                                                                        & coef_lo, coef_n, uniform, &
                                                                        & ${'.true.' if mapped_weno else '.false.'}$, &
                                                                        & epsilon_weno${', teno_cutoff' if teno else ''}$, &
                                                                        & state_left, state_right)
                                            #:endif
                                            if (state <= ${CELLS + 1}$) left${EQ}$(${LANE_INDEX}$state) = state_right
                                            if (state >= 2) right${EQ}$(${LANE_INDEX}$state - 1) = state_left
                                        #:endfor
                                        #:if GRAD_WENO
                                            ! Reconstruct the cell-centred gradients like the state;
                                            ! the pass below averages each face pair into gleft.
                                            #:for A in range(1, ND + 1)
                                                #:for V in range(1, NV + 1)
                                                    call s_fused_weno_state(gcolumn_${'xyz'[A - 1] + str(V)}$(1, 1), &
                                                                            & lane + ${LANES}$*(state + ${WENO_RADIUS - 1}$), &
                                                                            & normal, ${LANES}$, poly_l, poly_r, d_l, d_r, &
                                                                            & beta_coef, coef_lo, coef_n, uniform, &
                                                                            & ${'.true.' if mapped_weno else '.false.'}$, &
                                                                            & epsilon_weno${', teno_cutoff' if teno else ''}$, &
                                                                            & state_left, state_right)
                                                    if (state <= ${CELLS + 1}$) &
                                                        & gleft_${'xyz'[A - 1] + str(V)}$(${LANE_INDEX}$state) = state_right
                                                    if (state >= 2) gright_${'xyz'[A - 1] + str(V)}$(${LANE_INDEX}$state - 1) &
                                                        & = state_left
                                                #:endfor
                                            #:endfor
                                        #:endif
                                        #:if COMP
                                            compression_alpha = real(column${ENERGY + 1}$(${LANE_INDEX}$state + ${WENO_RADIUS}$), &
                                                                     & wp)
                                            call s_fused_thinc_${AXIS}$(compression_alpha, &
                                            #:if MTHINC
                                                & normal_hat, normal_position, ${NORMAL_INDEX}$, &
                                            #:else
                                                & real(column${ENERGY + 1}$(${LANE_INDEX}$state + ${WENO_RADIUS - 1}$), wp), &
                                                    & real(column${ENERGY + 1}$(${LANE_INDEX}$state + ${WENO_RADIUS + 1}$), wp), &
                                            #:endif
                                            & epsilon_ic, beta_ic, state_left, state_right, compress_state)
                                            if (compress_state) then
                                                compression_rho1 = real(column1(${LANE_INDEX}$state + ${WENO_RADIUS}$), &
                                                                        & wp)/compression_alpha
                                                compression_rho2 = real(column2(${LANE_INDEX}$state + ${WENO_RADIUS}$), &
                                                                        & wp)/(1._wp - compression_alpha)
                                                #:for SIDE, ALPHA, SLOT, GUARD in [('left', 'state_right', 'state', &
                                                                                        & f'state <= {CELLS + 1}'), ('right', &
                                                                                        & 'state_left', 'state - 1', 'state >= 2')]
                                                    if (${GUARD}$) then
                                                        ${SIDE}$1(${LANE_INDEX}$${SLOT}$) = compression_rho1*${ALPHA}$
                                                        ${SIDE}$2(${LANE_INDEX}$${SLOT}$) = compression_rho2*(1._wp - ${ALPHA}$)
                                                        ${SIDE}$${ENERGY + 1}$(${LANE_INDEX}$${SLOT}$) = ${ALPHA}$
                                                        ${SIDE}$${NADV}$(${LANE_INDEX}$${SLOT}$) = 1._wp - ${ALPHA}$
                                                    end if
                                                #:endfor
                                            end if
                                        #:endif
                                    end if
                                end do
                            end do
                        #:endif
                        #:if GRAD_WENO
                            ! Average the two one-sided gradients of each face, as in
                            ! s_average_fused_face_gradients. This pass reads the previous loop's stores.
                            $:GPU_LOOP(collapse=2, parallelism='[vector]')
                            do state = 2, ${CELLS + 2}$
                                do lane = 1, ${LANES}$
                                    #:for A in range(1, ND + 1)
                                        #:for V in range(1, NV + 1)
                                            gleft_${'xyz'[A - 1] + str(V)}$(${LANE_INDEX}$state - 1) &
                                                                       & = 5.e-1_wp*(gleft_${'xyz'[A - 1] + str(V)}$(${LANE_INDEX}$state - 1) + gright_${'xyz'[A - 1] + str(V)}$(${LANE_INDEX}$state - 1))
                                        #:endfor
                                    #:endfor
                                end do
                            end do
                        #:endif
                        $:GPU_LOOP(collapse=2, parallelism='[vector]', private='[normal, x, q_index, ql_index, q_left, q_right, &
                                   & flux, source' + (', bubble_transport, bubble_left, bubble_right' if bubbles_euler and qbmm else '') + (', geometric_flux' if GEO else '') + (', grad_index, viscous_velocity, re_interface, viscous_flux, divergence, tau' + (', alpha_avg, gamma_dot, d_xx, d_yy, d_zz, d_xy, d_xz, d_yz' if hllc_non_newtonian else '') + (', r_eff, div_v_term_const, sigma1, sigma2, sigma3, avg_v2, avg_v3' if CYL_STRESS else '') if IN_COLUMN_STRESS else ((', viscous_velocity' + (', re_interface' if viscous and not hllc_non_newtonian else '')) if viscous or surface_tension else '')) + ']')
                        do state = 1, ${CELLS + 1}$
                            do lane = 1, ${LANES}$
                                normal = block_n*${CELLS}$ - 2 + state
                                x = block_x*${LANES}$ + lane - 1
                                q_index = ${QINDEX}$
                                if (normal <= ${EXTENT}$${' .and. x <= m_in' if DIR > 1 else ''}$) then
                                    #:if STENCIL
                                        #:for EQ in range(1, min(NFACE, NRIEMANN) + 1)
                                            q_left(${EQ}$) = left${EQ}$(${LANE_INDEX}$state)
                                            q_right(${EQ}$) = right${EQ}$(${LANE_INDEX}$state)
                                        #:endfor
                                    #:elif WENO
                                        #:for EQ in range(1, NRIEMANN + 1)
                                            q_left(${EQ}$) = real(q${REC_FIELD(EQ)}$(q_index), wp)
                                            q_right(${EQ}$) = real(q${REC_FIELD(EQ)}$(q_index + ${['1', 'ni', 'ni*nj'][DIR - 1]}$) &
                                                    & , wp)
                                        #:endfor
                                    #:else
                                        ql_index = ${QLINDEX}$
                                        #:for EQ in range(1, NRIEMANN + 1)
                                            q_left(${EQ}$) = qR(ql_index + qln*(${REC_FIELD(EQ)}$ - 1))
                                            q_right(${EQ}$) = qL(ql_index + ${QLSTRIDE}$ + qln*(${REC_FIELD(EQ)}$ - 1))
                                        #:endfor
                                    #:endif
                                    #:if hllc_extrap
                                        if (normal == -1 .and. bc_${AXIS}$%beg == BC_RIEMANN_EXTRAP) q_left = q_right
                                        if (normal == ${EXTENT}$ .and. bc_${AXIS}$%end == BC_RIEMANN_EXTRAP) q_right = q_left
                                    #:endif
                                    call s_hllc_flux_${AXIS}$(q_left, q_right, flux, source, &
                                                              & ${VISC_RESULTS}$${BUB_RESULTS}$${'geometric_flux, ' if GEO else ''}$${PHYS_GROUPS[0]}$, &
                                                              & ${PHYS_GROUPS[1]}$, &
                                                              & ${', '.join(PHYS_GROUPS[2:])}$${', ' + MOM_VALUES if bubbles_euler and qbmm else ''}$)
                                    #:if bubbles_euler and qbmm
                                        ! Passive moments use the completed interface transport coefficients.
                                        #:for EQ in range(NADV + 2, NE + 1)
                                            #:if STENCIL
                                                bubble_left = left${EQ}$(${LANE_INDEX}$state)
                                                bubble_right = right${EQ}$(${LANE_INDEX}$state)
                                            #:elif WENO
                                                bubble_left = real(q${EQ}$(q_index), wp)
                                                bubble_right = real(q${EQ}$(q_index + ${['1', 'ni', 'ni*nj'][DIR - 1]}$), wp)
                                            #:else
                                                bubble_left = qR(ql_index + qln*(${EQ}$ - 1))
                                                bubble_right = qL(ql_index + ${QLSTRIDE}$ + qln*(${EQ}$ - 1))
                                            #:endif
                                            #:if hllc_extrap
                                                if (normal == -1 .and. bc_${AXIS}$%beg == BC_RIEMANN_EXTRAP) &
                                                    & bubble_left = bubble_right
                                                if (normal == ${EXTENT}$ .and. bc_${AXIS}$%end == BC_RIEMANN_EXTRAP) &
                                                    & bubble_right = bubble_left
                                            #:endif
                                            flux(${EQ}$) = bubble_transport(1)*bubble_left*bubble_transport(3) &
                                                 & + bubble_transport(2)*bubble_right*bubble_transport(4)
                                        #:endfor
                                    #:endif
                                    #:if IN_COLUMN_STRESS
                                        grad_index = ${GRAD_INDEX}$
                                        #:if hllc_non_newtonian
                                            ! Interface Reynolds numbers from the face strain rate and the
                                            ! averaged cell alphas, with the Herschel-Bulkley law per fluid.
                                            d_xx = ${GREAD('x', 1)}$
                                            d_yy = 0._wp; d_zz = 0._wp; d_xy = 0._wp; d_xz = 0._wp; d_yz = 0._wp
                                            #:if ND > 1
                                                d_yy = ${GREAD('y', 2)}$
                                                d_xy = 0.5_wp*(${GREAD('y', 1)}$ + ${GREAD('x', 2)}$)
                                            #:endif
                                            #:if ND > 2
                                                d_zz = ${GREAD('z', 3)}$
                                                d_xz = 0.5_wp*(${GREAD('z', 1)}$ + ${GREAD('x', 3)}$)
                                                d_yz = 0.5_wp*(${GREAD('z', 2)}$ + ${GREAD('y', 3)}$)
                                            #:endif
                                            gamma_dot = f_compute_shear_rate_from_components(d_xx, d_yy, d_zz, d_xy, d_xz, d_yz)
                                            re_interface = 0._wp
                                            #:for F in range(1, NF + 1)
                                                alpha_avg = 0.5_wp*(${NN_ALPHAS[F - 1]}$)
                                                alpha_avg = min(max(alpha_avg, 0._wp), 1._wp)
                                                #:if hllc_nn_fluids[F - 1]
                                                    re_interface(1) = re_interface(1) &
                                                                 & + alpha_avg*f_compute_hb_viscosity(hb_tau0_${F}$, hb_K_${F}$, &
                                                                 & hb_nn_${F}$, hb_mu_min_${F}$, hb_mu_max_${F}$, gamma_dot, &
                                                                 & hb_m_${F}$)
                                                #:else
                                                    if (re_shear${F}$ > 0._wp) re_interface(1) = re_interface(1) &
                                                        & + alpha_avg/re_shear${F}$
                                                    if (re_bulk${F}$ > 0._wp) re_interface(2) = re_interface(2) &
                                                        & + alpha_avg/re_bulk${F}$
                                                #:endif
                                            #:endfor
                                            re_interface(1) = 1._wp/max(re_interface(1), sgm_eps)
                                            re_interface(2) = 1._wp/max(re_interface(2), sgm_eps)
                                        #:endif
                                        #:if CYL_STRESS
                                            ! Cylindrical stress as in s_compute_cylindrical_viscous_source_flux,
                                            ! with r_eff the face radius on radial faces and the cell radius otherwise.
                                            #:if DIR == 2
                                                r_eff = radius_cb(normal - radius_lo + 1)
                                            #:else
                                                r_eff = radius_cc(${XYZ[1]}$ - radius_lo + 1)
                                            #:endif
                                            avg_v2 = 0.5_wp*(${CYL_VEL[0]}$)
                                            avg_v3 = 0.5_wp*(${CYL_VEL[1]}$)
                                            divergence = ${GREAD('x', 1)}$ + ${GREAD('y', 2)}$ + avg_v2/r_eff &
                                                                 & + ${GREAD('z', 3)}$/r_eff
                                            viscous_flux = 0._wp
                                            if (shear_active) then
                                                div_v_term_const = -(2._wp/3._wp)*divergence/re_interface(1)
                                                #:if DIR == 1
                                                    sigma1 = (2._wp*${GREAD('x', 1)}$)/re_interface(1) + div_v_term_const
                                                    sigma2 = (${GREAD('y', 1)}$ + ${GREAD('x', 2)}$)/re_interface(1)
                                                    sigma3 = (${GREAD('z', 1)}$/r_eff + ${GREAD('x', 3)}$)/re_interface(1)
                                                #:elif DIR == 2
                                                    sigma1 = (${GREAD('y', 1)}$ + ${GREAD('x', 2)}$)/re_interface(1)
                                                    sigma2 = (2._wp*${GREAD('y', 2)}$)/re_interface(1) + div_v_term_const
                                                    sigma3 = (${GREAD('z', 2)}$/r_eff - avg_v3/r_eff + ${GREAD('y', 3)}$) &
                                                              & /re_interface(1)
                                                #:else
                                                    sigma1 = (${GREAD('z', 1)}$/r_eff + ${GREAD('x', 3)}$)/re_interface(1)
                                                    sigma2 = (${GREAD('z', 2)}$/r_eff - avg_v3/r_eff + ${GREAD('y', 3)}$) &
                                                              & /re_interface(1)
                                                    sigma3 = (2._wp*(${GREAD('z', 3)}$/r_eff + avg_v2/r_eff))/re_interface(1) &
                                                              & + div_v_term_const
                                                #:endif
                                                viscous_flux(1) = viscous_flux(1) - sigma1
                                                viscous_flux(2) = viscous_flux(2) - sigma2
                                                viscous_flux(3) = viscous_flux(3) - sigma3
                                                viscous_flux(${ND + 1}$) = viscous_flux(${ND + 1}$) - viscous_velocity(1)*sigma1
                                                viscous_flux(${ND + 1}$) = viscous_flux(${ND + 1}$) - viscous_velocity(2)*sigma2
                                                viscous_flux(${ND + 1}$) = viscous_flux(${ND + 1}$) - viscous_velocity(3)*sigma3
                                            end if
                                            if (bulk_active) then
                                                viscous_flux(${DIR}$) = viscous_flux(${DIR}$) - divergence/re_interface(2)
                                                viscous_flux(${ND + 1}$) = viscous_flux(${ND + 1}$) - viscous_velocity(${DIR}$) &
                                                             & *divergence/re_interface(2)
                                            end if
                                        #:else
                                            divergence = ${' + '.join(GREAD(a, v) for v, a in enumerate('xyz'[:ND], 1))}$
                                            viscous_flux = 0._wp
                                            if (shear_active) then
                                                #:for V in range(1, ND + 1)
                                                    tau = (${GREAD(AXIS, V)}$ + ${GREAD('xyz'[V - 1], DIR)}$)/re_interface(1)
                                                    #:if V == DIR
                                                        tau = tau - (2._wp/3._wp)*divergence/re_interface(1)
                                                    #:endif
                                                    viscous_flux(${V}$) = viscous_flux(${V}$) - tau
                                                    viscous_flux(${ND + 1}$) = viscous_flux(${ND + 1}$) - viscous_velocity(${V}$) &
                                                                 & *tau
                                                #:endfor
                                            end if
                                            if (bulk_active) then
                                                #:for V in range(1, ND + 1)
                                                    tau = ${'divergence/re_interface(2)' if V == DIR else '0._wp'}$
                                                    viscous_flux(${V}$) = viscous_flux(${V}$) - tau
                                                    viscous_flux(${ND + 1}$) = viscous_flux(${ND + 1}$) - viscous_velocity(${V}$) &
                                                                 & *tau
                                                #:endfor
                                            end if
                                        #:endif
                                        #:for V in range(1, ND + 2)
                                            viscous${V}$(${LANE_INDEX}$state) = viscous_flux(${V}$)
                                        #:endfor
                                    #:endif
                                    #:if (viscous and not IN_COLUMN_STRESS) or surface_tension
                                        if (state > 1 .or. block_n == 0) then
                                            #:for V in range(1, ND + 1)
                                                vel_src_rsx_vf(${INDEX}$, ${V}$) = viscous_velocity(${V}$)
                                            #:endfor
                                            #:if viscous and not hllc_non_newtonian
                                                Re_avg_rsx_vf(${INDEX}$, 1) = re_interface(1)
                                                Re_avg_rsx_vf(${INDEX}$, 2) = re_interface(2)
                                            #:endif
                                            #:if surface_tension
                                                flux_src_rsx_vf(${INDEX}$, ${ENERGY + 1}$) = source
                                            #:endif
                                        end if
                                    #:endif
                                    #:if bubbles_euler and qbmm and not polytropic
                                        if (state > 1 .or. block_n == 0) then
                                            #:for BIN in range(nb)
                                                #:for EQ in [NADV + 1 + BIN*6, NADV + 2 + BIN*6, NADV + 4 + BIN*6]
                                                    flux_rsx_vf(${INDEX}$, ${EQ}$) = flux(${EQ}$)
                                                #:endfor
                                            #:endfor
                                        end if
                                    #:endif
                                    #:if hllc_cbc
                                        if ((state > 1 .or. block_n == 0) .and. ((normal_beg > 0 .and. normal < normal_beg) &
                                            & .or. (normal_end < ${EXTENT}$ .and. normal >= normal_end))) then
                                            #:for EQ in range(1, NE + 1)
                                                flux_rsx_vf(${INDEX}$, ${EQ}$) = flux(${EQ}$)
                                            #:endfor
                                            #:for F in SOURCE_ROWS
                                                flux_src_rsx_vf(${INDEX}$, ${ENERGY + F}$) = ${SOURCE(F)}$
                                            #:endfor
                                        end if
                                    #:endif
                                    #:for EQ in range(1, NE + 1)
                                        ${FACE_FLUX(EQ)}$(${LANE_INDEX}$state) = flux(${EQ}$)
                                    #:endfor
                                    #:for F in SOURCE_ROWS
                                        ${FACE_SOURCE(F)}$(${LANE_INDEX}$state) = ${SOURCE(F)}$
                                    #:endfor
                                    #:if GEO
                                        right${NF + 1}$(${LANE_INDEX}$state) = geometric_flux
                                    #:endif
                                end if
                            end do
                        end do
                        $:GPU_LOOP(collapse=2, parallelism='[vector]', &
                                   & private='[normal, x, q_index, rhs_index, inv_ds, value' &
                                   & + (', k_term' if hllc_alt_soundspeed else '') &
                                   & + (', radial_coord, half_over_radius' if GEO else '') + (', ' + ', '.join(UPDATE_PRIVATE) &
                                   & if UPDATE_PRIVATE else '') + ']')
                        do state = 1, ${CELLS}$
                            do lane = 1, ${LANES}$
                                normal = block_n*${CELLS}$ + state - 1
                                x = block_x*${LANES}$ + lane - 1
                                q_index = ${QINDEX}$
                                if (normal <= ${EXTENT}$${' .and. x <= m_in' if DIR > 1 else ''}$ ${'.and. normal >= normal_beg .and. normal <= normal_end' if hllc_cbc and not CBC_GEO else ''}$) then
                                    rhs_index = ${RINDEX}$
                                    inv_ds = 1._wp/ds(normal - ds_lo + 1)
                                    #:if GEO
                                        radial_coord = radius_cc(${XYZ[1]}$ - radius_lo + 1)
                                        half_over_radius = 0.5_wp/radial_coord
                                        #:if DIR == 3
                                            inv_ds = 1._wp/(ds(normal - ds_lo + 1)*radial_coord)
                                        #:endif
                                    #:endif
                                    #:if CBC_GEO
                                        if (normal >= normal_beg .and. normal <= normal_end) then
                                        #:endif
                                        #:if hllc_alt_soundspeed
                                            ! Kapila term of the first and last volume fractions, as in m_rhs.
                                            k_term = f_hllc_kapila(real(q${ENERGY}$(q_index), wp), real(q${ENERGY + 1}$(q_index), &
                                                                   & wp), real(q${NADV}$(q_index), wp), real(q${1}$(q_index), &
                                                                   & wp), real(q${NF}$(q_index), wp), gamma1, gamma2, pi1, pi2)
                                        #:endif
                                        #:for EQ in range(1, NE + 1)
                                            value = ${FACE_FLUX(EQ)}$(${LANE_INDEX}$state) &
                                                                & - ${FACE_FLUX(EQ)}$(${LANE_INDEX}$state + 1)
                                            #:if not SIXEQ and ENERGY < EQ <= NADV and (not hllc_alt_soundspeed or EQ in (ENERGY + 1, &
                                                    & NADV))
                                                #:if hllc_riemann == 1
                                                    value = value + real(q${NF + DIR}$(q_index), &
                                                                         & wp)*(${FACE_SOURCE(EQ - ENERGY)}$(${LANE_INDEX}$state) &
                                                                         & - ${FACE_SOURCE(EQ - ENERGY)}$(${LANE_INDEX}$state + 1))
                                                #:else
                                                    #:set ADVECTED_ALPHA = f'real(q{EQ}(q_index), wp)'
                                                    #:if hllc_alt_soundspeed
                                                        #:set ADVECTED_ALPHA = f'({ADVECTED_ALPHA} {"+" if EQ == ENERGY + 1 else "-"} k_term)'
                                                    #:endif
                                                    value = value + ${ADVECTED_ALPHA}$*(${FACE_SOURCE(1)}$(${LANE_INDEX}$ state &
                                                                                        & + 1) &
                                                                                        & - ${FACE_SOURCE(1)}$(${LANE_INDEX}$state))
                                                #:endif
                                            #:endif
                                            #:if SIXEQ
                                                value = real(${'real(rhs' + str(EQ) + '(rhs_index), wp) + ' if DIR > 1 else ''}$inv_ds*value, &
                                                             & stp)
                                                #:if ENERGY < EQ <= NADV
                                                    value = value + inv_ds*real(q${EQ}$(q_index), &
                                                                                & wp)*(${FACE_SOURCE(1)}$(${LANE_INDEX}$state + 1) &
                                                                                & - ${FACE_SOURCE(1)}$(${LANE_INDEX}$state))
                                                #:elif NADV < EQ <= NADV + NF
                                                    value = value - inv_ds*real(q${EQ - NF}$(q_index), &
                                                                                & wp)*real(q${ENERGY}$(q_index), &
                                                                                & wp)*(${FACE_SOURCE(1)}$(${LANE_INDEX}$state + 1) &
                                                                                & - ${FACE_SOURCE(1)}$(${LANE_INDEX}$state))
                                                    #:if GEO and DIR == 2
                                                        value = value - half_over_radius*real(q${EQ - NF}$(q_index), &
                                                                                              & wp)*real(q${ENERGY}$(q_index), &
                                                                                              & wp) &
                                                                                              & *(${FACE_SOURCE(1)}$(${LANE_INDEX}$state + 1) + ${FACE_SOURCE(1)}$(${LANE_INDEX}$state))
                                                    #:endif
                                                #:endif
                                            #:endif
                                            #:if GEO
                                                #:if not SIXEQ
                                                    value = real(real(rhs${EQ}$(rhs_index), wp) + inv_ds*value, stp)
                                                #:endif
                                                #:if DIR == 2 and (hllc_riemann == 1 or EQ <= ENERGY or (SIXEQ and NADV < EQ <= NADV + NF))
                                                    #:set GROW = f'right{NF + 1}' if EQ == NF + 2 else f'left{EQ}'
                                                    value = value - half_over_radius*(${GROW}$(${LANE_INDEX}$state) &
                                                                                      & + ${GROW}$(${LANE_INDEX}$state + 1))
                                                #:elif DIR == 3 and EQ in (NF + 2, NF + 3)
                                                    #:set GROW = f'right{NF + 1}' if EQ == NF + 2 else f'left{NF + 2}'
                                                    value = value ${'+' if EQ == NF + 2 else '-'}$ half_over_radius*(${GROW}$(${LANE_INDEX}$state) + ${GROW}$(${LANE_INDEX}$state + 1))
                                                #:endif
                                                #:if DIR == 2 and hllc_alt_soundspeed and EQ in (ENERGY + 1, NADV)
                                                    value = value ${'+' if EQ == ENERGY + 1 else '-'}$ (k_term/(2._wp*radial_coord))*(right1(${LANE_INDEX}$state + 1) + right1(${LANE_INDEX}$state))
                                                #:endif
                                                ${f'value{EQ}' if KEEP_VALUES else f'rhs{EQ}(rhs_index)'}$ = real(value, stp)
                                            #:elif SIXEQ
                                                ${f'value{EQ}' if KEEP_VALUES else f'rhs{EQ}(rhs_index)'}$ = real(value, stp)
                                            #:else
                                                ${f'value{EQ}' if KEEP_VALUES else f'rhs{EQ}(rhs_index)'}$ = real(${'real(rhs' + str(EQ) + '(rhs_index), wp) + ' if DIR > 1 else ''}$inv_ds*value, &
                                                    & stp)
                                            #:endif
                                        #:endfor
                                        #:if IN_COLUMN_STRESS
                                            #:if GEO and DIR == 2
                                                if (.not. (normal == 0 .and. (bc_y%beg == -2 .or. bc_y%beg == -14))) then
                                                #:endif
                                                #:for EQ in list(range(NF + 1, NF + ND + 1)) + [ENERGY]
                                                    value${EQ}$ = real(value${EQ}$ &
                                                                       & + ${'1._wp/ds(normal - ds_lo + 1)' if CYL_STRESS and DIR == 3 else 'inv_ds'}$*(viscous${ND + 1 if EQ == ENERGY else EQ - NF}$(${LANE_INDEX}$state) - viscous${ND + 1 if EQ == ENERGY else EQ - NF}$(${LANE_INDEX}$state + 1)), &
                                                                       & stp)
                                                    #:if GEO and DIR == 2
                                                        value${EQ}$ = real(value${EQ}$ &
                                                                           & - half_over_radius*(viscous${ND + 1 if EQ == ENERGY else EQ - NF}$(${LANE_INDEX}$state) + viscous${ND + 1 if EQ == ENERGY else EQ - NF}$(${LANE_INDEX}$state + 1)), &
                                                                           & stp)
                                                    #:endif
                                                #:endfor
                                                #:if CYL_STRESS and DIR == 3
                                                    ! Azimuthal viscous divergence with the plain 1/dz plus the
                                                    ! hoop-stress cross-momentum terms.
                                                    value${NF + 2}$ = real(value${NF + 2}$ &
                                                                           & + 5.e-1_wp*(viscous3(${LANE_INDEX}$state) &
                                                                           & + viscous3(${LANE_INDEX}$state + 1)), stp)
                                                    value${NF + 3}$ = real(value${NF + 3}$ &
                                                                           & - 5.e-1_wp*(viscous2(${LANE_INDEX}$state) &
                                                                           & + viscous2(${LANE_INDEX}$state + 1)), stp)
                                                #:endif
                                                #:if GEO and DIR == 2
                                                end if
                                            #:endif
                                        #:endif
                                        #:if FINAL
                                            if (fuse_final_rk) then
                                                $:hllc_rk_update((lambda field, &
                                                                 & eq: (f'column{eq}({LANE_INDEX}state + {WENO_RADIUS + 1})' if STENCIL and eq <= NFACE else f'q{eq}(q_index)') if field == 'stage' else f'{field}{eq}(q_index)'))
                                            else
                                            #:endif
                                            #:if KEEP_VALUES
                                                #:for EQ in range(1, NE + 1)
                                                    rhs${EQ}$(rhs_index) = real(value${EQ}$, stp)
                                                #:endfor
                                            #:endif
                                            #:if FINAL
                                            end if
                                        #:endif
                                        #:if CBC_GEO
                                        else
                                            ! CBC changes the face fluxes and keeps the original geometric source.
                                            #:for EQ in range(1, NE + 1)
                                                #:if hllc_riemann == 1 or EQ <= ENERGY or (SIXEQ and NADV < EQ <= NADV + NF)
                                                    #:set GROW = f'right{NF + 1}' if EQ == NF + 2 else f'left{EQ}'
                                                    rhs${EQ}$(rhs_index) = rhs${EQ}$(rhs_index) &
                                                              & - half_over_radius*(${GROW}$(${LANE_INDEX}$state) &
                                                              & + ${GROW}$(${LANE_INDEX}$state + 1))
                                                #:endif
                                            #:endfor
                                        end if
                                    #:endif
                                    #:if IN_COLUMN_STRESS and hllc_cbc and not CBC_GEO
                                        ! CBC band cells receive the viscous divergence here;
                                        ! s_update_fused_boundary_band adds the CBC-adjusted
                                        ! advection divergence on top.
                                    else if (normal <= ${EXTENT}$${' .and. x <= m_in' if DIR > 1 else ''}$) then
                                        rhs_index = ${RINDEX}$
                                        inv_ds = 1._wp/ds(normal - ds_lo + 1)
                                        #:for EQ in list(range(NF + 1, NF + ND + 1)) + [ENERGY]
                                            rhs${EQ}$(rhs_index) &
                                                      & = real(${'real(rhs' + str(EQ) + '(rhs_index), wp) + ' if DIR > 1 else ''}$inv_ds*(viscous${ND + 1 if EQ == ENERGY else EQ - NF}$(${LANE_INDEX}$state) - viscous${ND + 1 if EQ == ENERGY else EQ - NF}$(${LANE_INDEX}$state + 1)), &
                                                      & stp)
                                        #:endfor
                                    #:endif
                                end if
                            end do
                        end do
                    #:endif

                end subroutine ${NAME}$_${SUFFIX}$_${AXIS}$
            #:endfor
        #:endfor
    #:endfor

    #:if MFC_CASE_OPTIMIZATION and viscous and hllc_vel_gradients
        !> Cell-centred velocity gradients for the weno_Re_flux viscous stress, one kernel per direction in place of the
        !! s_get_viscous pass. It reconstructs with s_fused_weno_state and closes with the s_apply_scalar_divergence_theorem
        !! expression; faces outside the s_weno write range at the buffered ends of the swept direction contribute zero.
        #:for DIR, AXIS in enumerate('xyz'[:ND], 1)
            #:set EXTENT = ['m_in', 'n_in', 'p_in'][DIR - 1]
            #:set COORD = ['j', 'k', 'l'][DIR - 1]
            #:set STRIDE = ['1', 'ni', 'ni*nj'][DIR - 1]
            #:set QIDX = '1 + j - i_lo + ni*((k - j_lo) + nj*(l - k_lo))'
            #:set GIDX = '1 + j - i_lo + gi*((k - j_lo) + gj*(l - k_lo))'
            #:set MAPPED = '.true.' if mapped_weno else '.false.'
            subroutine s_fused_velocity_gradients_${AXIS}$(${', '.join(f'q_v{v}' for v in range(1, ND + 1))}$, &
                & ${', '.join(f'grad_v{v}' for v in range(1, ND + 1))}$, m_in, n_in, p_in, buff_in, i_lo, j_lo, k_lo, ni, nj, gi, &
                & gj, poly_l, poly_r, d_l, d_r, beta_coef, coef_lo, coef_n, uniform, &
                & epsilon_weno${', teno_cutoff' if teno else ''}$, wa, ds, ds_lo)

                real(stp), intent(in)       :: ${', '.join(f'q_v{v}(*)' for v in range(1, ND + 1))}$
                real(stp), intent(inout)    :: ${', '.join(f'grad_v{v}(*)' for v in range(1, ND + 1))}$
                integer, value, intent(in)  :: m_in, n_in, p_in, buff_in, i_lo, j_lo, k_lo, ni, nj, gi, gj, coef_lo, coef_n, ds_lo
                real(wp), intent(in)        :: poly_l(*), poly_r(*), d_l(*), d_r(*), beta_coef(*), ds(*)
                logical, value, intent(in)  :: uniform
                real(wp), value, intent(in) :: epsilon_weno${', teno_cutoff' if teno else ''}$, wa

                call s_fused_velocity_gradients_${AXIS}$_launch(${', '.join(f'q_v{v}' for v in range(1, ND + 1))}$, &
                    & ${', '.join(f'grad_v{v}' for v in range(1, ND + 1))}$, m_in, n_in, p_in, buff_in, i_lo, j_lo, k_lo, ni, nj, &
                    & gi, gj, poly_l, poly_r, d_l, d_r, beta_coef, coef_lo, coef_n, uniform, &
                    & epsilon_weno${', teno_cutoff' if teno else ''}$, wa, ds, ds_lo)

            end subroutine s_fused_velocity_gradients_${AXIS}$

            ! The launch receives the coefficient tables as assumed-size dummies, so the data region
            ! checks their presence by address.
            subroutine s_fused_velocity_gradients_${AXIS}$_launch(${', '.join(f'q_v{v}' for v in range(1, ND + 1))}$, &
                & ${', '.join(f'grad_v{v}' for v in range(1, ND + 1))}$, m_in, n_in, p_in, buff_in, i_lo, j_lo, k_lo, ni, nj, gi, &
                & gj, poly_l, poly_r, d_l, d_r, beta_coef, coef_lo, coef_n, uniform, &
                & epsilon_weno${', teno_cutoff' if teno else ''}$, wa, ds, ds_lo)

                real(stp), intent(in)       :: ${', '.join(f'q_v{v}(*)' for v in range(1, ND + 1))}$
                real(stp), intent(inout)    :: ${', '.join(f'grad_v{v}(*)' for v in range(1, ND + 1))}$
                integer, value, intent(in)  :: m_in, n_in, p_in, buff_in, i_lo, j_lo, k_lo, ni, nj, gi, gj, coef_lo, coef_n, ds_lo
                real(wp), intent(in)        :: poly_l(*), poly_r(*), d_l(*), d_r(*), beta_coef(*), ds(*)
                logical, value, intent(in)  :: uniform
                real(wp), value, intent(in) :: epsilon_weno${', teno_cutoff' if teno else ''}$, wa
                integer                     :: j, k, l
                integer                     :: jb, je, kb, ke, lb, le
                real(wp)                    :: vL_c, vR_c, vL_p, vR_m, q_left, q_right

                ! Loop bounds are plain locals; derived-type components in region bounds miscompile on nvfortran.

                jb = idwbuff(1)%beg; je = idwbuff(1)%end
                kb = idwbuff(2)%beg; ke = idwbuff(2)%end
                lb = idwbuff(3)%beg; le = idwbuff(3)%end

                $:GPU_PARALLEL_LOOP(collapse=3, private='[vL_c, vR_c, vL_p, vR_m, q_left, q_right]')
                do l = lb, le
                    do k = kb, ke
                        do j = jb, je
                            ! Write range of s_apply_scalar_divergence_theorem, one cell inside the buffer along the sweep.
                            if (${COORD}$ >= 1 - buff_in .and. ${COORD}$ <= ${EXTENT}$ + buff_in - 1) then
                                #:for V in range(1, ND + 1)
                                    ! Edge values of the cell inside the s_weno write range, q_left at the
                                    ! left boundary and q_right at the right.
                                    if (${COORD}$ >= ${weno_polyn}$ - buff_in .and. ${COORD}$ <= ${EXTENT}$ + buff_in &
                                        & - ${weno_polyn}$) then
                                        call s_fused_weno_state(q_v${V}$, ${QIDX}$, ${COORD}$, ${STRIDE}$, poly_l, poly_r, d_l, &
                                                                & d_r, beta_coef, coef_lo, coef_n, uniform, ${MAPPED}$, &
                                                                & epsilon_weno${', teno_cutoff' if teno else ''}$, q_left, q_right)
                                        vL_c = q_left
                                        vR_c = q_right
                                    else
                                        vL_c = 0._wp
                                        vR_c = 0._wp
                                    end if
                                    vL_p = 0._wp
                                    vR_m = 0._wp
                                    if (wa /= 0._wp) then
                                        if (${COORD}$ > ${weno_polyn}$ - buff_in) then
                                            call s_fused_weno_state(q_v${V}$, ${QIDX}$ - ${STRIDE}$, ${COORD}$ - 1, ${STRIDE}$, &
                                                                    & poly_l, poly_r, d_l, d_r, beta_coef, coef_lo, coef_n, &
                                                                    & uniform, ${MAPPED}$, &
                                                                    & epsilon_weno${', teno_cutoff' if teno else ''}$, q_left, &
                                                                    & q_right)
                                            vR_m = q_right
                                        end if
                                        if (${COORD}$ < ${EXTENT}$ + buff_in - ${weno_polyn}$) then
                                            call s_fused_weno_state(q_v${V}$, ${QIDX}$ + ${STRIDE}$, ${COORD}$ + 1, ${STRIDE}$, &
                                                                    & poly_l, poly_r, d_l, d_r, beta_coef, coef_lo, coef_n, &
                                                                    & uniform, ${MAPPED}$, &
                                                                    & epsilon_weno${', teno_cutoff' if teno else ''}$, q_left, &
                                                                    & q_right)
                                            vL_p = q_left
                                        end if
                                    end if
                                    grad_v${V}$(${GIDX}$) = 1._wp/((1._wp + wa)*ds(${COORD}$ - ds_lo + 1))*(wa*vL_p + vR_c - vL_c &
                                                & - wa*vR_m)
                                #:endfor
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

            end subroutine s_fused_velocity_gradients_${AXIS}$_launch
        #:endfor
    #:endif
#:endif
#endif
end module m_fused_rhs_openacc
