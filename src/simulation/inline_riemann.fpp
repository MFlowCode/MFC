#:def arithmetic_avg()
    rho_avg = 5d-1*(rho_L + rho_R)
    vel_avg_rms = 0d0
    !$acc loop seq
    do i = 1, num_dims
        vel_avg_rms = vel_avg_rms + (5d-1*(vel_L(i) + vel_R(i)))**2d0
    end do

    H_avg = 5d-1*(H_L + H_R)
    gamma_avg = 5d-1*(gamma_L + gamma_R)



#:enddef arithmetic_avg


#:def roe_avg()
    rho_avg = sqrt(rho_L*rho_R)
    vel_avg_rms = 0d0
    !$acc loop seq
    do i = 1, num_dims
        vel_avg_rms = vel_avg_rms + (sqrt(rho_L)*vel_L(i) + sqrt(rho_R)*vel_R(i))**2d0/ &
                      (sqrt(rho_L) + sqrt(rho_R))**2d0
    end do

    H_avg = (sqrt(rho_L)*H_L + sqrt(rho_R)*H_R)/ &
            (sqrt(rho_L) + sqrt(rho_R))

    gamma_avg = (sqrt(rho_L)*gamma_L + sqrt(rho_R)*gamma_R)/ &
                (sqrt(rho_L) + sqrt(rho_R))

    rho_avg = sqrt(rho_L*rho_R)
    vel_avg_rms = (sqrt(rho_L)*vel_L(1) + sqrt(rho_R)*vel_R(1))**2d0/ &
                  (sqrt(rho_L) + sqrt(rho_R))**2d0

#:enddef roe_avg

#:def compute_average_state()

if (avg_state == 1) then
    @:roe_avg()
end if

if (avg_state == 2) then
    @:arithmetic_avg()
end if

#:enddef compute_average_state

! #:def compute_wave_speeds()

! if (wave_speeds == 1) then
!     s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
!     s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

!     s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
!            (s_L - vel_L(dir_idx(1))) - &
!            rho_R*vel_R(dir_idx(1))* &
!            (s_R - vel_R(dir_idx(1)))) &
!           /(rho_L*(s_L - vel_L(dir_idx(1))) - &
!             rho_R*(s_R - vel_R(dir_idx(1))))
! elseif (wave_speeds == 2) then
!     pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
!                     (vel_L(dir_idx(1)) - &
!                      vel_R(dir_idx(1))))

!     pres_SR = pres_SL

!     Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
!                          (pres_SL/pres_L - 1d0)*pres_L/ &
!                          ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
!     Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
!                          (pres_SR/pres_R - 1d0)*pres_R/ &
!                          ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

!     s_L = vel_L(dir_idx(1)) - c_L*Ms_L
!     s_R = vel_R(dir_idx(1)) + c_R*Ms_R

!     s_S = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
!                 (pres_L - pres_R)/ &
!                 (rho_avg*c_avg))
! end if


!         if (wave_speeds == 1) then
!             s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
!             s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

!             s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
!                    (s_L - vel_L(dir_idx(1))) - &
!                    rho_R*vel_R(dir_idx(1))* &
!                    (s_R - vel_R(dir_idx(1)))) &
!                   /(rho_L*(s_L - vel_L(dir_idx(1))) - &
!                     rho_R*(s_R - vel_R(dir_idx(1))))
!         elseif (wave_speeds == 2) then
!             pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
!                             (vel_L(dir_idx(1)) - &
!                              vel_R(dir_idx(1))))

!             pres_SR = pres_SL

!             Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
!                                  (pres_SL/pres_L - 1d0)*pres_L/ &
!                                  ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
!             Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
!                                  (pres_SR/pres_R - 1d0)*pres_R/ &
!                                  ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

!             s_L = vel_L(dir_idx(1)) - c_L*Ms_L
!             s_R = vel_R(dir_idx(1)) + c_R*Ms_R

!             s_S = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
!                         (pres_L - pres_R)/ &
!                         (rho_avg*c_avg))
!         end if

! #:enddef compute_wave_speeds


