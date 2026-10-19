! HLL transport velocity shared by NC terms and viscous work.
#:def hll_interface_velocity(SL, SR, VL, VR, OUT)
    if (0._wp <= ${SL}$) then
        ${OUT}$ = ${VL}$
    else if (${SR}$ <= 0._wp) then
        ${OUT}$ = ${VR}$
    else
        ${OUT}$ = (${SR}$*${VL}$ - ${SL}$*${VR}$)/(${SR}$ - ${SL}$)
    end if
#:enddef
