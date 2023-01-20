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
