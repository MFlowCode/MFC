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

    vel_avg_rms = (sqrt(rho_L)*vel_L(1) + sqrt(rho_R)*vel_R(1))**2d0/ &
                  (sqrt(rho_L) + sqrt(rho_R))**2d0

    if (chemistry) then
        eps = 0.001d0
        call get_species_enthalpies_rt(T_L, h_iL)
        call get_species_enthalpies_rt(T_R, h_iR)

        h_iL = h_iL*gas_constant/mol_weights*T_L
        h_iR = h_iR*gas_constant/mol_weights*T_R
        call get_species_specific_heats_r(T_L, Cp_iL)
        call get_species_specific_heats_r(T_R, Cp_iR)

        h_avg_2 = (sqrt(rho_L)*h_iL + sqrt(rho_R)*h_iR)/(sqrt(rho_L) + sqrt(rho_R))
        Yi_avg = (sqrt(rho_L)*Ys_L + sqrt(rho_R)*Ys_R)/(sqrt(rho_L) + sqrt(rho_R))
        T_avg = (sqrt(rho_L)*T_L + sqrt(rho_R)*T_R)/(sqrt(rho_L) + sqrt(rho_R))

        if (abs(T_L - T_R) < eps) then
            ! Case when T_L and T_R are very close
            Cp_avg = sum(Yi_avg(:)*(0.5d0*Cp_iL(i) + 0.5d0*Cp_iR(:))*gas_constant/mol_weights(:))
            Cv_avg = sum(Yi_avg(:)*((0.5d0*Cp_iL(i) + 0.5d0*Cp_iR(:))*gas_constant/mol_weights(:) - gas_constant/mol_weights(:)))
        else
            ! Normal calculation when T_L and T_R are sufficiently different
            Cp_avg = sum(Yi_avg(:)*(h_iR(:) - h_iL(:))/(T_R - T_L))
            Cv_avg = sum(Yi_avg(:)*((h_iR(:) - h_iL(:))/(T_R - T_L) - gas_constant/mol_weights(:)))
        end if

        gamma_avg = Cp_avg/Cv_avg

        Phi_avg(:) = (gamma_avg - 1.d0)*(vel_avg_rms/2.0d0 - h_avg_2(:)) + gamma_avg*gas_constant/mol_weights(:)*T_avg
        c_sum_Yi_Phi = sum(Yi_avg(:)*Phi_avg(:))
    end if

#:enddef roe_avg

#:def compute_average_state()

    if (avg_state == 1) then
        @:roe_avg()
    end if

    if (avg_state == 2) then
        @:arithmetic_avg()
    end if

#:enddef compute_average_state

#:def compute_low_Mach_correction()

    zcoef = min(1d0, max(vel_L_rms**5d-1/c_L, vel_R_rms**5d-1/c_R))
    pcorr = 0d0

    if (low_Mach == 1) then
        pcorr = rho_L*rho_R* &
                (s_L - vel_L(dir_idx(1)))*(s_R - vel_R(dir_idx(1)))*(vel_R(dir_idx(1)) - vel_L(dir_idx(1)))/ &
                (rho_R*(s_R - vel_R(dir_idx(1))) - rho_L*(s_L - vel_L(dir_idx(1))))* &
                (zcoef - 1d0)
    else if (low_Mach == 2) then
        vel_L_tmp = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + zcoef*(vel_L(dir_idx(1)) - vel_R(dir_idx(1))))
        vel_R_tmp = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + zcoef*(vel_R(dir_idx(1)) - vel_L(dir_idx(1))))
        vel_L(dir_idx(1)) = vel_L_tmp
        vel_R(dir_idx(1)) = vel_R_tmp
    end if

#:enddef compute_low_Mach_correction
