#:def Hardcoded2DVariables()

    real(kind(0d0)) :: eps
    real(kind(0d0)) :: r, rmax, gam, umax, p0
    real(kind(0d0)) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, intL, alph

    eps = 1e-9

#:enddef

#:def Hardcoded2D()

    select case (patch_icpp(patch_id)%hcid) ! 2D_hardcoded_ic example case
    
    case (200)
        if (y_cc(j) <= (-x_cc(i)**3 + 1)**(1d0/3d0)) then
            ! Volume Fractions
            q_prim_vf(advxb)%sf(i, j, 0) = eps
            q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - eps
            ! Denssities
            q_prim_vf(contxb)%sf(i, j, 0) = eps*1000d0
            q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - eps)*1d0
            ! Pressure
            q_prim_vf(E_idx)%sf(i, j, 0) = 1000d0
        end if
    case (202) ! Gresho vortex (Gouasmi et al 2022 JCP)
        r = ((x_cc(i) - 0.5d0)**2 + (y_cc(j) - 0.5d0)**2)**0.5d0
        rmax = 0.2

        gam = 1d0 + 1d0/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1d0/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5d0)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5d0)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5d0)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2d0)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5d0)/r)*umax*(2d0 - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5d0)/r)*umax*(2d0 - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2d0 + 4*(1 - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0d0
            q_prim_vf(momxe)%sf(i, j, 0) = 0d0
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*(-2 + 4*log(2.))
        end if
    case (203) ! Gresho vortex (Gouasmi et al 2022 JCP) with density correction
        r = ((x_cc(i) - 0.5d0)**2 + (y_cc(j) - 0.5d0)**2)**0.5d0
        rmax = 0.2

        gam = 1d0 + 1d0/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1d0/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5d0)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5d0)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5d0)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2d0)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5d0)/r)*umax*(2d0 - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5d0)/r)*umax*(2d0 - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2d0 + 4*(1 - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0d0
            q_prim_vf(momxe)%sf(i, j, 0) = 0d0
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*(-2 + 4*log(2.))
        end if

        q_prim_vf(contxb)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)**(1d0/gam)

 case (204) ! Rayleigh-taylor problem
        rhoH = 3
        rhoL = 1
        pRef = 1e5
        pInt = pRef
        h = 0.7
        lam = 0.2
        wl = 2*pi/lam
        amp = 0.05/wl

        intH = amp*sin(2*pi*x_cc(i)/lam - pi/2) + h

        alph = 5d-1*(1 + tanh((y_cc(j) - intH)/2.5e-3))

        if (alph < eps) alph = eps
        if (alph > 1 - eps) alph = 1 - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1 - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1 - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, 0) = pref + rhoH*9.81*(1.2 - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1 - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1 - alph)*rhoL
            pInt = pref + rhoH*9.81*(1.2 - intH)
            q_prim_vf(E_idx)%sf(i, j, 0) = pInt + rhoL*9.81*(intH - y_cc(j))
        end if

    case (205) ! 2D lung wave interaction problem
        h = 0.0           !non dim origin y
        lam = 1.0         !non dim lambda
        amp =  patch_icpp(patch_id)%a2         !to be changed later!       !non dim amplitude       

        intH = amp*sin(2*pi*x_cc(i)/lam - pi/2)+h

       if (y_cc(j) > intH) then       
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
       end if
       
     case (206) ! 2D lung wave interaction problem - horizontal domain
        h = 0.0           !non dim origin y
        lam = 1.0         !non dim lambda
        amp =  patch_icpp(patch_id)%a2        
        
        intL = amp*sin(2*pi*y_cc(j)/lam - pi/2)+h

       if (x_cc(i) > intL) then        !this is the liquid
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
       end if

     case (207) ! Bumps for the patch geometry of the lung
        h = 0.0
        lam = 1.0
        amp =  patch_icpp(patch_id)%a2 
        
        
       
    case default
       if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
        end if
        
    end select

#:enddef
