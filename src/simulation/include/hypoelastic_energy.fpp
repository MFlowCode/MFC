! Shared elastic interface energy for the case-specialized HLL paths.
! Each face retains its elastic energy when its own shear modulus is non-negligible.
#:def hypoelastic_interface_energy(ALPHA, TAU, MODULUS, ENERGY, DAMAGE=None)
    #:for SIDE in ['L', 'R']
        ${MODULUS(SIDE)}$ = 0._wp
        #:for F in range(1, num_fluids + 1)
            ${MODULUS(SIDE)}$ = ${MODULUS(SIDE)}$ + ${ALPHA(SIDE, F)}$*Gs_rs(${F}$)
        #:endfor
    #:endfor
    ! Elastic energy uses the undamaged modulus, so this loop precedes the damage scaling
    #:for SIDE in ['L', 'R']
        if (${MODULUS(SIDE)}$ > verysmall) then
            #:for ROW in range(1, num_dims*(num_dims + 1)//2 + int(cyl_coord) + 1)
                ${ENERGY(SIDE)}$ = ${ENERGY(SIDE)}$ + ${TAU(SIDE, ROW)}$**2/(4._wp*${MODULUS(SIDE)}$)
                #:if ROW in (2, 4, 5) and ROW <= num_dims*(num_dims + 1)//2
                    ${ENERGY(SIDE)}$ = ${ENERGY(SIDE)}$ + ${TAU(SIDE, ROW)}$**2/(4._wp*${MODULUS(SIDE)}$)
                #:endif
            #:endfor
        end if
    #:endfor
    ! Damage scales the returned moduli, which set the wave speeds
    #:if DAMAGE is not None
        if (cont_damage) then
            #:for SIDE in ['L', 'R']
                ${MODULUS(SIDE)}$ = ${MODULUS(SIDE)}$*max(1._wp - ${DAMAGE(SIDE)}$, 0._wp)
            #:endfor
        end if
    #:endif
#:enddef
