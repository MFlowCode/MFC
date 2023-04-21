#:def s_compute_speed_of_sound()
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, c)

        real(wp), intent(IN) :: pres
        real(wp), intent(IN) :: rho, gamma, pi_inf
        real(wp), intent(IN) :: H
        real(wp), dimension(num_fluids), intent(IN) :: adv
        real(wp), intent(IN) :: vel_sum
        real(wp), intent(OUT) :: c

        real(wp) :: blkmod1, blkmod2

        integer :: q

        if (alt_soundspeed) then 
            blkmod1 = ((gammas(1) + 1d0)*pres + & 
                        pi_infs(1))/gammas(1)
            blkmod2 = ((gammas(2) + 1d0)*pres + & 
                        pi_infs(2))/gammas(2) 
            c = (1d0/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2))) 
        elseif (model_eqns == 3) then 
            c = 0d0 
!$acc loop seq 
            do q = 1, num_fluids 
                c = c + adv(q)*(1d0/gammas(q) + 1d0)* & 
                    (pres + pi_infs(q)/(gammas(q) + 1d0)) 
            end do 
            c = c/rho

        elseif (((model_eqns == 4) .or. (model_eqns == 2 .and. bubbles))) then
            ! Sound speed for bubble mmixture to order O(\alpha)

            if (mpp_lim .and. (num_fluids > 1)) then
                c = (1d0/gamma + 1d0)* &
                      (pres + pi_inf)/rho
            else
                c = &
                    (1d0/gamma + 1d0)* &    
                    (pres + pi_inf)/ &
                    (rho*(1d0 - adv(num_fluids)))
            end if
        else 
            c = ((H - 5d-1*vel_sum)/gamma) 
        end if 

        if (mixture_err .and. c < 0d0) then
            c = 100.d0*sgm_eps
        else
            c = sqrt(c)
        end if
    end subroutine s_compute_speed_of_sound
#:enddef

