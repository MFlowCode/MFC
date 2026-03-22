#:def Hardcoded2DVariables()
    ! Place any declaration of intermediate variables here
    real(wp) :: eps, eps_mhd, C_mhd
    real(wp) :: r, rmax, gam, umax, p0
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, intL, alph
    real(wp) :: factor
    real(wp) :: r0, alpha, r2
    real(wp) :: sinA, cosA

    real(wp) :: r_sq

    ! # 207
    real(wp) :: sigma, gauss1, gauss2
    ! # 208
    real(wp) :: ei, d, fsm, alpha_air, alpha_sf6

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

    case (207) ! Kelvin Helmholtz Instability
        sigma = 0.05_wp/sqrt(2.0_wp)
        gauss1 = exp(-(y_cc(j) - 0.75_wp)**2/(2.0_wp*sigma**2))
        gauss2 = exp(-(y_cc(j) - 0.25_wp)**2/(2.0_wp*sigma**2))
        q_prim_vf(momxb + 1)%sf(i, j, 0) = &
            0.1_wp*sin(4.0_wp*pi*x_cc(i))*(gauss1 + gauss2)

    case (208) ! Richtmeyer Meshkov Instability
        lam = 1.0_wp
        eps = 1.0e-6_wp
        ei = 5.0_wp
        ! Smoothening function to smooth out sharp discontinuity in the interface
        if (x_cc(i) <= 0.7_wp*lam) then
            d = x_cc(i) - lam*(0.4_wp - 0.1_wp*sin(2.0_wp*pi*(y_cc(j)/lam + 0.25_wp)))
            fsm = 0.5_wp*(1.0_wp + erf(d/(ei*sqrt(dx*dy))))
            alpha_air = eps + (1.0_wp - 2.0_wp*eps)*fsm
            alpha_sf6 = 1.0_wp - alpha_air
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_sf6*5.04_wp
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_air*1.0_wp
            q_prim_vf(advxb)%sf(i, j, 0) = alpha_sf6
            q_prim_vf(advxe)%sf(i, j, 0) = alpha_air
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

        ! case 252 is for the 2D MHD Rotor problem
    case (252) ! 2D MHD Rotor Problem
        ! Ambient conditions are set in the JSON file.
        ! This case imposes the dense, rotating cylinder.
        !
        ! gamma = 1.4
        ! Ambient medium (r > 0.1):
        !   rho = 1, p = 1, v = 0, B = (1,0,0)
        ! Rotor (r <= 0.1):
        !   rho = 10, p = 1
        !   v has angular velocity w=20, giving v_tan=2 at r=0.1

        ! Calculate distance squared from the center
        r_sq = (x_cc(i) - 0.5_wp)**2 + (y_cc(j) - 0.5_wp)**2

        ! inner radius of 0.1
        if (r_sq <= 0.1**2) then
            ! -- Inside the rotor --
            ! Set density uniformly to 10
            q_prim_vf(contxb)%sf(i, j, 0) = 10._wp

            ! Set vup constant rotation of rate v=2
            ! v_x = -omega * (y - y_c)
            ! v_y =  omega * (x - x_c)
            q_prim_vf(momxb)%sf(i, j, 0) = -20._wp*(y_cc(j) - 0.5_wp)
            q_prim_vf(momxb + 1)%sf(i, j, 0) = 20._wp*(x_cc(i) - 0.5_wp)

            ! taper width of 0.015
        else if (r_sq <= 0.115**2) then
            ! linearly smooth the function between r = 0.1 and 0.115
            q_prim_vf(contxb)%sf(i, j, 0) = 1._wp + 9._wp*(0.115_wp - sqrt(r_sq))/(0.015_wp)

            q_prim_vf(momxb)%sf(i, j, 0) = -(2._wp/sqrt(r_sq))*(y_cc(j) - 0.5_wp)*(0.115_wp - sqrt(r_sq))/(0.015_wp)
            q_prim_vf(momxb + 1)%sf(i, j, 0) = (2._wp/sqrt(r_sq))*(x_cc(i) - 0.5_wp)*(0.115_wp - sqrt(r_sq))/(0.015_wp)
        end if

    case (253) ! MHD Smooth Magnetic Vortex
        ! Section 5.2 of
        ! Implicit hybridized discontinuous Galerkin methods for compressible magnetohydrodynamics
        ! C. Ciuca, P. Fernandez, A. Christophe, N.C. Nguyen, J. Peraire

        ! velocity
        q_prim_vf(momxb)%sf(i, j, 0) = 1._wp - (y_cc(j)*exp(1 - (x_cc(i)**2 + y_cc(j)**2))/(2.*pi))
        q_prim_vf(momxb + 1)%sf(i, j, 0) = 1._wp + (x_cc(i)*exp(1 - (x_cc(i)**2 + y_cc(j)**2))/(2.*pi))

        ! magnetic field
        q_prim_vf(B_idx%beg)%sf(i, j, 0) = -y_cc(j)*exp(1 - (x_cc(i)**2 + y_cc(j)**2))/(2.*pi)
        q_prim_vf(B_idx%beg + 1)%sf(i, j, 0) = x_cc(i)*exp(1 - (x_cc(i)**2 + y_cc(j)**2))/(2.*pi)

        ! pressure
        q_prim_vf(E_idx)%sf(i, j, 0) = 1._wp + (1 - 2._wp*(x_cc(i)**2 + y_cc(j)**2))*exp(1 - (x_cc(i)**2 + y_cc(j)**2))/((2._wp*pi)**3)

    case (260)  ! Gaussian Divergence Pulse
        !  Bx(x) = 1 + C * erf((x-0.5)/σ)
        !  ⇒  ∂Bx/∂x = C * (2/√π) * exp[-((x-0.5)/σ)**2] * (1/σ)
        !  Choose C = ε * σ * √π / 2  ⇒ ∂Bx/∂x = ε * exp[-((x-0.5)/σ)**2]
        !  ψ is initialized to zero everywhere.

        eps_mhd = patch_icpp(patch_id)%a(2)
        sigma = patch_icpp(patch_id)%a(3)
        C_mhd = eps_mhd*sigma*sqrt(pi)*0.5_wp

        ! B-field
        q_prim_vf(B_idx%beg)%sf(i, j, 0) = 1._wp + C_mhd*erf((x_cc(i) - 0.5_wp)/sigma)

    case (261)  ! Blob
        r0 = 1._wp/sqrt(8._wp)
        r2 = x_cc(i)**2 + y_cc(j)**2
        r = sqrt(r2)
        alpha = r/r0
        if (alpha < 1) then
            q_prim_vf(B_idx%beg)%sf(i, j, 0) = 1._wp/sqrt(4._wp*pi)*(alpha**8 - 2._wp*alpha**4 + 1._wp)
            ! q_prim_vf(B_idx%beg)%sf(i,j,0) = 1._wp/sqrt(4000._wp*pi) * (4096._wp*r2**4 - 128._wp*r2**2 + 1._wp)
            ! q_prim_vf(B_idx%beg)%sf(i,j,0) = 1._wp/(4._wp*pi) * (alpha**8 - 2._wp*alpha**4 + 1._wp)
            ! q_prim_vf(E_idx)%sf(i,j,0) = 6._wp - q_prim_vf(B_idx%beg)%sf(i,j,0)**2/2._wp
        end if

    case (262)  ! Tilted 2D MHD shock‐tube at α = arctan2 (≈63.4°)
        ! rotate by α = atan(2)
        alpha = atan(2._wp)
        cosA = cos(alpha)
        sinA = sin(alpha)
        ! projection along shock normal
        r = x_cc(i)*cosA + y_cc(j)*sinA

        if (r <= 0.5_wp) then
            ! LEFT state: ρ=1, v∥=+10, v⊥=0, p=20, B∥=B⊥=5/√(4π)
            q_prim_vf(contxb)%sf(i, j, 0) = 1._wp
            q_prim_vf(momxb)%sf(i, j, 0) = 10._wp*cosA
            q_prim_vf(momxb + 1)%sf(i, j, 0) = 10._wp*sinA
            q_prim_vf(E_idx)%sf(i, j, 0) = 20._wp
            q_prim_vf(B_idx%beg)%sf(i, j, 0) = (5._wp/sqrt(4._wp*pi))*cosA &
                                               - (5._wp/sqrt(4._wp*pi))*sinA
            q_prim_vf(B_idx%beg + 1)%sf(i, j, 0) = (5._wp/sqrt(4._wp*pi))*sinA &
                                                   + (5._wp/sqrt(4._wp*pi))*cosA
        else
            ! RIGHT state: ρ=1, v∥=−10, v⊥=0, p=1, B∥=B⊥=5/√(4π)
            q_prim_vf(contxb)%sf(i, j, 0) = 1._wp
            q_prim_vf(momxb)%sf(i, j, 0) = -10._wp*cosA
            q_prim_vf(momxb + 1)%sf(i, j, 0) = -10._wp*sinA
            q_prim_vf(E_idx)%sf(i, j, 0) = 1._wp
            q_prim_vf(B_idx%beg)%sf(i, j, 0) = (5._wp/sqrt(4._wp*pi))*cosA &
                                               - (5._wp/sqrt(4._wp*pi))*sinA
            q_prim_vf(B_idx%beg + 1)%sf(i, j, 0) = (5._wp/sqrt(4._wp*pi))*sinA &
                                                   + (5._wp/sqrt(4._wp*pi))*cosA
        end if
        ! v^z and B^z remain zero by default

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
