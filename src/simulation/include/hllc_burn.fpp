!> Evaluate the programmed-burn rate.
#:def hllc_burn_rate(PRES, LAMBDA, ARHO, ALPHA)
    burn_rate = 0._wp
    burn_drive = (${PRES}$ - rburn%pign)/rburn%pref
    if (burn_drive > 0._wp .and. ${LAMBDA}$ < 1._wp) then
        burn_rate = rburn%k*(1._wp - ${LAMBDA}$)*burn_drive**rburn%n
        if (rburn%ta > 0._wp) then
            #:if any_state_dependent_eos
                ! Reactant phasic temperature from its own EOS, as in s_burn_rate.
                call s_phase_temperature(${ARHO}$/max(${ALPHA}$, sgm_eps), ${PRES}$, 1, burn_temperature)
            #:else
                burn_temperature = (${PRES}$ + isentrope_B(1))/((isentrope_n(1) - 1._wp)*cvs(1)*${ARHO}$/max(${ALPHA}$, sgm_eps))
            #:endif
            burn_rate = burn_rate*exp(-rburn%ta/burn_temperature)
        end if
    end if
#:enddef
