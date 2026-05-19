#:def arithmetic_avg()
    rho_avg = 5.e-1_wp*(rho%L + rho%R)
    vel_avg_rms = 0._wp
    $:GPU_LOOP(parallelism='[seq]')
    do i = 1, num_vels
        vel_avg_rms = vel_avg_rms + (5.e-1_wp*(vel%L(i) + vel%R(i)))**2._wp
    end do

    H_avg = 5.e-1_wp*(H%L + H%R)
    gamma_avg = 5.e-1_wp*(gamma%L + gamma%R)
    qv_avg = 5.e-1_wp*(qv%L + qv%R)
#:enddef arithmetic_avg

#:def roe_avg()
    rho_avg = sqrt(rho%L*rho%R)

    vel_avg_rms = 0._wp

    $:GPU_LOOP(parallelism='[seq]')
    do i = 1, num_vels
        vel_avg_rms = vel_avg_rms + (sqrt(rho%L)*vel%L(i) + sqrt(rho%R)*vel%R(i))**2._wp/(sqrt(rho%L) + sqrt(rho%R))**2._wp
    end do

    H_avg = (sqrt(rho%L)*H%L + sqrt(rho%R)*H%R)/(sqrt(rho%L) + sqrt(rho%R))

    gamma_avg = (sqrt(rho%L)*gamma%L + sqrt(rho%R)*gamma%R)/(sqrt(rho%L) + sqrt(rho%R))

    vel_avg_rms = (sqrt(rho%L)*vel%L(1) + sqrt(rho%R)*vel%R(1))**2._wp/(sqrt(rho%L) + sqrt(rho%R))**2._wp

    qv_avg = (sqrt(rho%L)*qv%L + sqrt(rho%R)*qv%R)/(sqrt(rho%L) + sqrt(rho%R))

    if (chemistry) then
        eps = 0.001_wp
        call get_species_enthalpies_rt(T%L, h_iL)
        call get_species_enthalpies_rt(T%R, h_iR)
        h_iL = h_iL*gas_constant/molecular_weights*T%L
        h_iR = h_iR*gas_constant/molecular_weights*T%R
        call get_species_specific_heats_r(T%L, Cp_iL)
        call get_species_specific_heats_r(T%R, Cp_iR)

        h_avg_2 = (sqrt(rho%L)*h_iL + sqrt(rho%R)*h_iR)/(sqrt(rho%L) + sqrt(rho%R))
        Yi_avg = (sqrt(rho%L)*Ys_L + sqrt(rho%R)*Ys_R)/(sqrt(rho%L) + sqrt(rho%R))
        T_avg = (sqrt(rho%L)*T%L + sqrt(rho%R)*T%R)/(sqrt(rho%L) + sqrt(rho%R))
        if (abs(T%L - T%R) < eps) then
            ! Case when T%L and T%R are very close
            Cp_avg = sum(Yi_avg(:)*(0.5_wp*Cp_iL(:) + 0.5_wp*Cp_iR(:))*gas_constant/molecular_weights(:))
            Cv_avg = sum(Yi_avg(:)*((0.5_wp*Cp_iL(:) + 0.5_wp*Cp_iR(:))*gas_constant/molecular_weights(:) &
                         & - gas_constant/molecular_weights(:)))
        else
            ! Normal calculation when T%L and T%R are sufficiently different
            Cp_avg = sum(Yi_avg(:)*(h_iR(:) - h_iL(:))/(T%R - T%L))
            Cv_avg = sum(Yi_avg(:)*((h_iR(:) - h_iL(:))/(T%R - T%L) - gas_constant/molecular_weights(:)))
        end if
        gamma_avg = Cp_avg/Cv_avg

        Phi_avg(:) = (gamma_avg - 1._wp)*(vel_avg_rms/2.0_wp - h_avg_2(:)) + gamma_avg*gas_constant/molecular_weights(:)*T_avg
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
    if (riemann_solver == 1 .or. riemann_solver == 5) then
        zcoef = min(1._wp, max(vel_rms%L**5.e-1_wp/c%L, vel_rms%R**5.e-1_wp/c%R))
        pcorr = 0._wp

        if (low_Mach == 1) then
            pcorr = -(s_P - s_M)*(rho%L + rho%R)/8._wp*(zcoef - 1._wp)
        end if
    else if (riemann_solver == 2) then
        zcoef = min(1._wp, max(vel_rms%L**5.e-1_wp/c%L, vel_rms%R**5.e-1_wp/c%R))
        pcorr = 0._wp

        if (low_Mach == 1) then
            pcorr = rho%L*rho%R*(s%L - vel%L(dir_idx(1)))*(s%R - vel%R(dir_idx(1)))*(vel%R(dir_idx(1)) - vel%L(dir_idx(1))) &
                                 & /(rho%R*(s%R - vel%R(dir_idx(1))) - rho%L*(s%L - vel%L(dir_idx(1))))*(zcoef - 1._wp)
        else if (low_Mach == 2) then
            vel_tmp%L = 5.e-1_wp*((vel%L(dir_idx(1)) + vel%R(dir_idx(1))) + zcoef*(vel%L(dir_idx(1)) - vel%R(dir_idx(1))))
            vel_tmp%R = 5.e-1_wp*((vel%L(dir_idx(1)) + vel%R(dir_idx(1))) + zcoef*(vel%R(dir_idx(1)) - vel%L(dir_idx(1))))
            vel%L(dir_idx(1)) = vel_tmp%L
            vel%R(dir_idx(1)) = vel_tmp%R
        end if
    end if
#:enddef compute_low_Mach_correction
