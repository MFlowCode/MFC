!>  This procedure calculates the speed of sound
        !! @param pres Pressure
        !! @param rho Cell averaged density
        !! @param pi_inf Cell averaged liquid stiffness
        !! @param gamma Cell averaged specific heat ratio
        !! @param H Cell averaged enthalpy
        !! @param adv Advection Variables
        !! @param vel_sum Sum of all velocities
        !! @param q_prim_vf Primitive vars in 1 direction
        !! @param flg Helps determine which conditionals to be called.
            ! flg >= 2: Check all conditionals
            ! flg =  1: Check for alt_soundspeed, otherwise run the 3rd conditional block
            ! flg = 0: Check for alt_soundspeed, otherwise use enthalpy
        !! @param c Speed of sound
        
#:def compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, q_prim_vf, j, k, l, flg, c)

    if (alt_soundspeed .and. ${flg}$ >= 0) then 
        blkmod1 = ((gammas(1) + 1d0)*${pres}$ + & 
                    pi_infs(1))/gammas(1) 
        blkmod2 = ((gammas(2) + 1d0)*${pres}$ + & 
                    pi_infs(2))/gammas(2) 
        ${c}$ = (1d0/(${rho}$*(${adv}$(1)/blkmod1 + ${adv}$(2)/blkmod2))) 
    elseif (model_eqns == 3 .and. ${flg}$ >= 2) then 
        ${c}$ = 0d0 
!$acc loop seq 
        do q = 1, num_fluids 
            ${c}$ = ${c}$ + ${q_prim_vf}$(${j}$, ${k}$, ${l}$, q + advxb - 1)*(1d0/gammas(q) + 1d0)* & 
                (${q_prim_vf}$(${j}$, ${k}$, ${l}$, E_idx) + pi_infs(q)/(gammas(q) + 1d0)) 
        end do 
        ${c}$ = ${c}$/${rho}$

    elseif (((model_eqns == 4) .or. (model_eqns == 2 .and. bubbles) .or. ${flg}$ == 1) &
            .and. ${flg}$ >= 1) then
        ! Sound speed for bubble mmixture to order O(\alpha)

        if (mpp_lim .and. (num_fluids > 1)) then
            ${c}$ = (1d0/${gamma}$ + 1d0)* &
                  (${pres}$ + ${pi_inf}$)/${rho}$
        else
            ${c}$ = &
                (1d0/${gamma}$ + 1d0)* &    
                (${pres}$ + ${pi_inf}$)/ &
                (${rho}$*(1d0 - ${adv}$(num_fluids)))
        end if
    else 
        ${c}$ = ((${H}$ - 5d-1*${vel_sum}$)/${gamma}$) 
    end if 

    if (mixture_err .and. ${c}$ < 0d0) then
        ${c}$ = 100.d0*sgm_eps
    else
        ${c}$ = sqrt(${c}$)
    end if
        
#:enddef compute_speed_of_sound


    !> Computes the bubble number density n from the conservative variables
        !! @param vftmp is the void fraction
        !! @param nRtmp is the bubble number  density times the bubble radii
        !! @param ntmp is the output number bubble density
#:def s_comp_n_from_cons(vftmp, nRtmp, ntmp, weight)
    !$acc routine seq
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: nRtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: nR3
        real(kind(0.d0)), dimension(nb) :: weight

        call s_quad(nRtmp**3.d0, nR3, weight)  !returns itself if NR0 = 1
        ntmp = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)
        ! ntmp = 1d0

#:enddef s_comp_n_from_cons