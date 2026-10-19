! Chemistry face closure shared by the fused column and the Riemann solvers.
#:if not defined('hllc_scalar_thermo')
    #:set hllc_scalar_thermo = False
#:endif
#:def chemistry_face(SIDE, RHO, PRESSURE, VELOCITY_SQ, GAMMA, ENERGY, TEMP, WEIGHT, ENTHALPY=None, GAMMA_GUARD=None, &
                     & FRACTION=None, MU=None)
    #:if MU is not None and not hllc_scalar_thermo
        #:assert False, 'chemistry_face: the inline mixture viscosity needs the scalar thermodynamic closure'
    #:endif
    #:if hllc_scalar_thermo
        #:set Y = lambda i: FRACTION(i) if FRACTION is not None else f'Ys_{SIDE}({i})'
        #:set LIT = lambda x: repr(float(x)) + '_wp'
        R_gas_${SIDE}$ = 0._wp
        #:for SP in range(hllc_num_species)
            R_gas_${SIDE}$ = R_gas_${SIDE}$ + ${Y(SP + 1)}$*${LIT(hllc_thermo_inv_w[SP])}$
        #:endfor
        ${WEIGHT}$ = 1._wp/R_gas_${SIDE}$
        R_gas_${SIDE}$ = gas_constant/${WEIGHT}$
        ${TEMP}$ = ${PRESSURE}$/${RHO}$/R_gas_${SIDE}$
        ${ENERGY}$ = 0._wp; Cv_${SIDE}$ = 0._wp
        #:for SP in range(hllc_num_species)
            if (${TEMP}$ > ${LIT(hllc_thermo_mid[SP])}$) then
                #:for RANGE in [1, 0]
                    #:if RANGE == 0
                    else
                    #:endif
                    #:set H = hllc_thermo_h[SP][RANGE]
                    #:set CP = hllc_thermo_cp[SP][RANGE]
                    Gamm_${SIDE}$ = ${TEMP}$**2*${LIT(H[2])}$
                    Gamm_${SIDE}$ = ${TEMP}$**3*${LIT(H[3])}$ + Gamm_${SIDE}$
                    Gamm_${SIDE}$ = ${TEMP}$**4*${LIT(H[4])}$ + Gamm_${SIDE}$
                    Gamm_${SIDE}$ = (1._wp/${TEMP}$)*${LIT(H[5])}$ + Gamm_${SIDE}$
                    Gamm_${SIDE}$ = (${TEMP}$*${LIT(H[1])}$ + ${LIT(H[0])}$) + Gamm_${SIDE}$
                    #:if GAMMA_GUARD is not None
                        if (${GAMMA_GUARD}$) then
                        #:endif
                        Cp_${SIDE}$ = ${TEMP}$*${LIT(CP[1])}$ + ${LIT(CP[0])}$
                        Cp_${SIDE}$ = ${TEMP}$**2*${LIT(CP[2])}$ + Cp_${SIDE}$
                        Cp_${SIDE}$ = ${TEMP}$**3*${LIT(CP[3])}$ + Cp_${SIDE}$
                        Cp_${SIDE}$ = ${TEMP}$**4*${LIT(CP[4])}$ + Cp_${SIDE}$
                        #:if GAMMA_GUARD is not None
                        end if
                    #:endif
                #:endfor
            end if
            ${ENERGY}$ = ${ENERGY}$ + (Gamm_${SIDE}$ - 1._wp)*(${Y(SP + 1)}$*${LIT(hllc_thermo_inv_w[SP])}$)
            #:if GAMMA_GUARD is not None
                if (${GAMMA_GUARD}$) then
                #:endif
                Cv_${SIDE}$ = Cv_${SIDE}$ + (Cp_${SIDE}$ - 1._wp)*(${Y(SP + 1)}$*${LIT(hllc_thermo_inv_w[SP])}$)
                #:if GAMMA_GUARD is not None
                end if
            #:endif
        #:endfor
        #:if GAMMA_GUARD is not None
            if (${GAMMA_GUARD}$) then
            #:endif
            ${GAMMA}$ = (Cv_${SIDE}$*gas_constant)/R_gas_${SIDE}$
            #:if GAMMA_GUARD is not None
            end if
        #:endif
        ${ENERGY}$ = ${TEMP}$*(${ENERGY}$*gas_constant)
        ${ENERGY}$ = ${RHO}$*${ENERGY}$ + 0.5_wp*${RHO}$*${VELOCITY_SQ}$
        #:if MU is not None
            ! Wilke mixture viscosity, in the evaluation order of get_mixture_viscosity_mixavg.
            log_T_${SIDE}$ = log(${TEMP}$)
            #:for SP in range(hllc_num_species)
                mu_${SIDE}$_${SP + 1}$ = sqrt(${TEMP}$)*(${LIT(hllc_thermo_mu[SP][0])}$ &
                                              & + ${LIT(hllc_thermo_mu[SP][1])}$*log_T_${SIDE}$ &
                                              & + ${LIT(hllc_thermo_mu[SP][2])}$*log_T_${SIDE}$**2 &
                                              & + ${LIT(hllc_thermo_mu[SP][3])}$*log_T_${SIDE}$**3 &
                                              & + ${LIT(hllc_thermo_mu[SP][4])}$*log_T_${SIDE}$**4)**2
            #:endfor
            ${MU}$ = 0._wp
            #:for I in range(hllc_num_species)
                wilke_phi_${SIDE}$ = 0._wp
                #:for J in range(hllc_num_species)
                    wilke_phi_${SIDE}$ = wilke_phi_${SIDE}$ + (${LIT(hllc_thermo_inv_w[J])}$*${Y(J + 1)}$*${WEIGHT}$)*(1._wp &
                        & + sqrt(mu_${SIDE}$_${I + 1}$/mu_${SIDE}$_${J + 1}$*${LIT(hllc_thermo_wilke_c[I][J])}$)) &
                        & **2/${LIT(hllc_thermo_wilke_d[I][J])}$
                #:endfor
                ${MU}$ = ${MU}$ + (${LIT(hllc_thermo_inv_w[I])}$*${Y(I + 1)}$*${WEIGHT}$)*mu_${SIDE}$_${I + 1}$/wilke_phi_${SIDE}$
            #:endfor
        #:endif
    #:else
        call get_mixture_molecular_weight(Ys_${SIDE}$, ${WEIGHT}$)
        R_gas_${SIDE}$ = gas_constant/${WEIGHT}$
        ${TEMP}$ = ${PRESSURE}$/${RHO}$/R_gas_${SIDE}$
        #:if GAMMA_GUARD is not None
            if (${GAMMA_GUARD}$) then
            #:endif
            #:if hllc_gamma_method == 1 and hllc_riemann in (1, 5)
                ! Sum X_i/(Gamma_i - 1) = Cv_mix/R_mix for ideal-gas species.
                call get_mixture_specific_heat_cv_mass(${TEMP}$, Ys_${SIDE}$, Cv_${SIDE}$)
                ${GAMMA}$ = Cv_${SIDE}$/R_gas_${SIDE}$
            #:elif hllc_gamma_method == 1
                Xs_${SIDE}$(:) = Ys_${SIDE}$(:)*${WEIGHT}$/molecular_weights(:)
                call get_species_specific_heats_r(${TEMP}$, Cp_i${SIDE}$)
                Gamma_i${SIDE}$ = Cp_i${SIDE}$/(Cp_i${SIDE}$ - 1._wp)
                ${GAMMA}$ = sum(Xs_${SIDE}$(:)/(Gamma_i${SIDE}$(:) - 1._wp))
            #:else
                call get_mixture_specific_heat_cp_mass(${TEMP}$, Ys_${SIDE}$, Cp_${SIDE}$)
                call get_mixture_specific_heat_cv_mass(${TEMP}$, Ys_${SIDE}$, Cv_${SIDE}$)
                Gamm_${SIDE}$ = Cp_${SIDE}$/Cv_${SIDE}$
                ${GAMMA}$ = 1._wp/(Gamm_${SIDE}$ - 1._wp)
            #:endif
            #:if GAMMA_GUARD is not None
            end if
        #:endif
        call get_mixture_energy_mass(${TEMP}$, Ys_${SIDE}$, ${ENERGY}$)
        ${ENERGY}$ = ${RHO}$*${ENERGY}$ + 0.5_wp*${RHO}$*${VELOCITY_SQ}$
    #:endif
    #:if ENTHALPY is not None
        ${ENTHALPY}$ = (${ENERGY}$ + ${PRESSURE}$)/${RHO}$
    #:endif
#:enddef
