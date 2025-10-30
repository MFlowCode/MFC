#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph, Mach

    real(wp) :: eps

    ! Case 302 - IGR Jet
    real(wp) :: r, ux_th, ux_am, p_th, p_am, rho_th, rho_am, y_th, z_th, r_th, eps_smooth

    eps = 1e-9_wp

#:enddef

#:def Hardcoded3D()

    select case (patch_icpp(patch_id)%hcid)
    case (300) ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.025_wp/wl

        intH = amp*(sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + sin(2._wp*pi*z_cc(k)/lam - pi/2._wp)) + h

        alph = 5.e-1_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, k) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, k) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if

    case (301) ! (3D lung geometry in X direction, |sin(*)+sin(*)|)
        h = 0.0_wp
        lam = 1.0_wp
        amp = patch_icpp(patch_id)%a(2)
        intH = amp*abs((sin(2*pi*y_cc(j)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h)
        if (x_cc(i) > intH) then
            q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, k) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, k) = patch_icpp(1)%alpha(2)
        end if

    case (302) ! 3D Jet with IGR

        ux_th = 10*sqrt(1.4*0.4)
        ux_am = 0.0*sqrt(1.4)
        p_th = 2.0_wp
        p_am = 1.0_wp
        rho_th = 1._wp
        rho_am = 1._wp
        y_th = 0.0_wp
        z_th = 0.0_wp
        r_th = 1._wp
        eps_smooth = 1._wp
        eps = 1e-6

        r = sqrt((y_cc(j) - y_th)**2._wp + (z_cc(k) - z_th)**2._wp)

        q_prim_vf(momxb)%sf(i, j, k) = ux_th*f_cut_on(r - r_th, eps_smooth)*f_cut_on(x_cc(i), eps_smooth) + ux_am
        q_prim_vf(momxb + 1)%sf(i, j, k) = 0._wp
        q_prim_vf(momxe)%sf(i, j, k) = 0._wp

        q_prim_vf(advxb)%sf(i, j, k) = (1._wp - 2._wp*eps)*f_cut_on(r - r_th, eps_smooth)*f_cut_on(x_cc(i), eps_smooth) + eps

        q_prim_vf(contxb)%sf(i, j, k) = q_prim_vf(advxb)%sf(i, j, k)*rho_am
        q_prim_vf(contxe)%sf(i, j, k) = (1._wp - q_prim_vf(advxb)%sf(i, j, k))*rho_th

        q_prim_vf(E_idx)%sf(i, j, k) = p_th*f_cut_on(r - r_th, eps_smooth)*f_cut_on(x_cc(i), eps_smooth) + p_am

    case (370)
        ! This hardcoded case extrudes a 2D profile to initialize a 3D simulation domain
        @: HardcodedReadValues()

    case (380)
        ! This is patch is hard-coded for test suite optimization used in the
        ! 3D_TaylorGreenVortex case:
        ! This analytic patch used geometry 9
        Mach = 0.1
        if (patch_id == 1) then
            q_prim_vf(E_idx)%sf(i, j, k) = 101325 + (Mach**2*376.636429464809**2/16)*(cos(2*x_cc(i)/1) + cos(2*y_cc(j)/1))*(cos(2*z_cc(k)/1) + 2)
            q_prim_vf(momxb + 0)%sf(i, j, k) = Mach*376.636429464809*sin(x_cc(i)/1)*cos(y_cc(j)/1)*sin(z_cc(k)/1)
            q_prim_vf(momxb + 1)%sf(i, j, k) = -Mach*376.636429464809*cos(x_cc(i)/1)*sin(y_cc(j)/1)*sin(z_cc(k)/1)
        end if

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
