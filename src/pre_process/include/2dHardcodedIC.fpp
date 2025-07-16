#:def Hardcoded2DVariables()
    ! Place any declaration of intermediate variables here
    real(wp) :: eps
    real(wp) :: r, rmax, gam, umax, p0
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, intL, alph
    real(wp) :: factor

    eps = 1.e-9_wp
#:enddef

#:def Hardcoded2D()

    select case (patch_icpp(patch_id)%hcid) ! 2D_hardcoded_ic example case

    case (200)
        if (y_cc(j) <= (-x_cc(i)**3 + 1)**(1._wp/3._wp)) then
            ! Volume Fractions
            q_prim_vf(advxb)%sf(i, j, 0) = eps
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - eps
            ! Denssities
            q_prim_vf(contxb)%sf(i, j, 0) = eps*1000._wp
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - eps)*1._wp
            ! Pressure
            q_prim_vf(E_idx)%sf(i, j, 0) = 1000._wp
        end if
    case (202) ! Gresho vortex (Gouasmi et al 2022 JCP)
        r = ((x_cc(i) - 0.5_wp)**2 + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2._wp/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4*(1 - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*(-2 + 4*log(2._wp))
        end if
    case (203) ! Gresho vortex (Gouasmi et al 2022 JCP) with density correction
        r = ((x_cc(i) - 0.5_wp)**2._wp + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2._wp/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4._wp*(1._wp - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2._wp*(-2._wp + 4*log(2._wp))
        end if

        q_prim_vf(contxb)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)**(1._wp/gam)
    case (204) ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.05_wp/wl

        intH = amp*sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + h

        alph = 0.5_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, 0) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, 0) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if

    case (205) ! 2D lung wave interaction problem
        h = 0.0_wp           !non dim origin y
        lam = 1.0_wp         !non dim lambda
        amp = patch_icpp(patch_id)%a(2)         !to be changed later!       !non dim amplitude

        intH = amp*sin(2*pi*x_cc(i)/lam - pi/2) + h

        if (y_cc(j) > intH) then
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if

    case (206) ! 2D lung wave interaction problem - horizontal domain
        h = 0.0_wp           !non dim origin y
        lam = 1.0_wp         !non dim lambda
        amp = patch_icpp(patch_id)%a(2)

        intL = amp*sin(2*pi*y_cc(j)/lam - pi/2) + h

        if (x_cc(i) > intL) then        !this is the liquid
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if

    case (250) ! MHD Orszag-Tang vortex
        ! gamma = 5/3
        !   rho = 25/(36*pi)
        !     p = 5/(12*pi)
        !     v = (-sin(2*pi*y), sin(2*pi*x), 0)
        !     B = (-sin(2*pi*y)/sqrt(4*pi), sin(4*pi*x)/sqrt(4*pi), 0)

        q_prim_vf(momxb)%sf(i, j, 0) = -sin(2._wp*pi*y_cc(j))
        q_prim_vf(momxb + 1)%sf(i, j, 0) = sin(2._wp*pi*x_cc(i))

        q_prim_vf(B_idx%beg)%sf(i, j, 0) = -sin(2._wp*pi*y_cc(j))/sqrt(4._wp*pi)
        q_prim_vf(B_idx%beg + 1)%sf(i, j, 0) = sin(4._wp*pi*x_cc(i))/sqrt(4._wp*pi)

    case (251) ! RMHD Cylindrical Blast Wave [Mignone, 2006: Section 4.3.1]
        if (x_cc(i)**2 + y_cc(j)**2 < 0.08_wp**2) then
            q_prim_vf(contxb)%sf(i, j, 0) = 0.01
            q_prim_vf(E_idx)%sf(i, j, 0) = 1.0
        elseif (x_cc(i)**2 + y_cc(j)**2 <= 1._wp**2) then
            ! Linear interpolation between r=0.08 and r=1.0
            factor = (1.0_wp - sqrt(x_cc(i)**2 + y_cc(j)**2))/(1.0_wp - 0.08_wp)
            q_prim_vf(contxb)%sf(i, j, 0) = 0.01_wp*factor + 1.e-4_wp*(1.0_wp - factor)
            q_prim_vf(E_idx)%sf(i, j, 0) = 1.0_wp*factor + 3.e-5_wp*(1.0_wp - factor)
        else
            q_prim_vf(contxb)%sf(i, j, 0) = 1.e-4_wp
            q_prim_vf(E_idx)%sf(i, j, 0) = 3.e-5_wp
        end if

    case (270)
        ! This hardcoded case extrudes a 1D profile to initialize a 2D simulation domain
        @: HardcodedReadValues()

    case (280)
        ! This is patch is hard-coded for test suite optimization used in the
        ! 2D_isentropicvortex case:
        ! This analytic patch uses geometry 2
        if (patch_id == 1) then
            q_prim_vf(E_idx)%sf(i, j, 0) = 1.0*(1.0 - (1.0/1.0)*(5.0/(2.0*pi))*(5.0/(8.0*1.0*(1.4 + 1.0)*pi))*exp(2.0*1.0*(1.0 - (x_cc(i) - patch_icpp(1)%x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0)))**(1.4 + 1.0)
            q_prim_vf(contxb + 0)%sf(i, j, 0) = 1.0*(1.0 - (1.0/1.0)*(5.0/(2.0*pi))*(5.0/(8.0*1.0*(1.4 + 1.0)*pi))*exp(2.0*1.0*(1.0 - (x_cc(i) - patch_icpp(1)%x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0)))**1.4
            q_prim_vf(momxb + 0)%sf(i, j, 0) = 0.0 + (y_cc(j) - patch_icpp(1)%y_centroid)*(5.0/(2.0*pi))*exp(1.0*(1.0 - (x_cc(i) - patch_icpp(1)%x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0))
            q_prim_vf(momxb + 1)%sf(i, j, 0) = 0.0 - (x_cc(i) - patch_icpp(1)%x_centroid)*(5.0/(2.0*pi))*exp(1.0*(1.0 - (x_cc(i) - patch_icpp(1)%x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0))
        end if

    case (281)
        ! This is patch is hard-coded for test suite optimization used in the
        ! 2D_acoustic_pulse case:
        ! This analytic patch uses geometry 2
        if (patch_id == 2) then
            q_prim_vf(E_idx)%sf(i, j, 0) = 101325*(1 - 0.5*(1.4 - 1)*(0.4)**2*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2))))**(1.4/(1.4 - 1))
            q_prim_vf(contxb + 0)%sf(i, j, 0) = 1*(1 - 0.5*(1.4 - 1)*(0.4)**2*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2))))**(1/(1.4 - 1))
        end if

    case (282)
        ! This is patch is hard-coded for test suite optimization used in the
        ! 2D_zero_circ_vortex case:
        ! This analytic patch uses geometry 2
        if (patch_id == 2) then
            q_prim_vf(E_idx)%sf(i, j, 0) = 101325*(1 - 0.5*(1.4 - 1)*(0.1/0.3)**2*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2))))**(1.4/(1.4 - 1))
            q_prim_vf(contxb + 0)%sf(i, j, 0) = 1*(1 - 0.5*(1.4 - 1)*(0.1/0.3)**2*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2))))**(1/(1.4 - 1))
            q_prim_vf(momxb + 0)%sf(i, j, 0) = 112.99092883944267*(1 - (0.1/0.3))*y_cc(j)*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2)))
            q_prim_vf(momxb + 1)%sf(i, j, 0) = 112.99092883944267*((0.1/0.3))*x_cc(i)*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2)))
        end if

    case default
        if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
        end if

    end select

#:enddef
