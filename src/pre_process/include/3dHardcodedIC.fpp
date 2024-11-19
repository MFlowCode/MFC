#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here

    real(kind(0d0)) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

    real(kind(0d0)) :: eps
   
    real(kind(0d0)) :: rcoord, theta, phi, xi_sph, x_bcen, y_bcen, z_bcen, Rinit
    real(kind(0d0)) :: x_ccs, y_ccs, z_ccs
    real(kind(0d0)), dimension(num_dims) :: xi_cart
    integer :: l

    eps = 1e-9
#:enddef

#:def Hardcoded3D()

    select case (patch_icpp(patch_id)%hcid)
    case (300) ! Rayleigh-Taylor instability
        rhoH = 3
        rhoL = 1
        pRef = 1e5
        pInt = pRef
        h = 0.7
        lam = 0.2
        wl = 2*pi/lam
        amp = 0.025/wl

        intH = amp*(sin(2*pi*x_cc(i)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h

        alph = 5d-1*(1 + tanh((y_cc(j) - intH)/2.5e-3))

        if (alph < eps) alph = eps
        if (alph > 1 - eps) alph = 1 - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1 - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1 - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, k) = pref + rhoH*9.81*(1.2 - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1 - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1 - alph)*rhoL
            pInt = pref + rhoH*9.81*(1.2 - intH)
            q_prim_vf(E_idx)%sf(i, j, k) = pInt + rhoL*9.81*(intH - y_cc(j))
        end if

    case (301) ! (3D lung geometry in X direction, |sin(*)+sin(*)|)
        h = 0.0
        lam = 1.0
        amp = patch_icpp(patch_id)%a2
        intH = amp*abs((sin(2*pi*y_cc(j)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h)
        if (x_cc(i) > intH) then
            q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, k) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, k) = patch_icpp(1)%alpha(2)
        end if
    
    case (302) ! (3D lung geometry in X direction - axisym, with smoothing)
        lam = 1.0*200.E-06
        amp = patch_icpp(patch_id)%a2
        h = 0.125*amp
        
        intH = amp/2*(sin(2*pi*y_cc(j)/lam + pi/2) + sin(2*pi*z_cc(k)/lam + pi/2))
        
        alph = patch_icpp(2)%alpha(1) + (patch_icpp(1)%alpha(1)-patch_icpp(2)%alpha(1))/(h)*(x_cc(i)-(intH-h/2))
        
        if (x_cc(i) > intH + h/2) then
        
            q_prim_vf(advxb)%sf(i, j, k) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, k) = patch_icpp(1)%alpha(2)
            q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(1)%pres
        
        else if ((x_cc(i) .le. intH + h/2) .and. (x_cc(i) .ge. intH - h/2)) then
        
            q_prim_vf(advxb)%sf(i, j, k) = alph !0.5
            q_prim_vf(advxe)%sf(i, j, k) = 1- alph !0.5
            q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(1)%alpha_rho(1)/patch_icpp(1)%alpha(1)*alph!0.5
            q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(2)%alpha_rho(2)/patch_icpp(2)%alpha(2)*(1-alph)!0.5
            q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(1)%pres
        
        end if
    
    case (303) ! pre_stress for hyperelasticity, bubble in material

        R0ref = 30E-6    ! equilibrium radius
        Rinit = patch_icpp(3)%radius ! initial radius
        x_bcen = patch_icpp(3)%x_centroid
        y_bcen = patch_icpp(3)%y_centroid
        z_bcen = patch_icpp(3)%z_centroid
        x_ccs = x_cc(i) - x_bcen
        y_ccs = y_cc(j) - y_bcen
        z_ccs = z_cc(k) - z_bcen
        rcoord = sqrt(x_ccs**2 + y_ccs**2 + z_ccs**2)
        phi = atan2(y_ccs, x_ccs)
        theta = atan2(sqrt(x_ccs**2 + y_ccs**2), z_ccs)
        !spherical coord, assuming Rmax=1
        xi_sph = (rcoord**3 - R0ref**3 + Rinit**3)**(1d0/3d0) 
        xi_cart(1) = xi_sph*sin(theta)*cos(phi)
        xi_cart(2) = xi_sph*sin(theta)*sin(phi)
        xi_cart(3) = xi_sph*cos(theta)
        ! shift back
        xi_cart(1) = xi_cart(1) + x_bcen
        xi_cart(2) = xi_cart(2) + y_bcen
        xi_cart(3) = xi_cart(3) + z_bcen
!        print *, 'xi_cart(1) ::', xi_cart(1), 'xi_cart(2) ::', xi_cart(2), 'xi_cart(3) ::', xi_cart(3)
        ! assigning the reference map to the q_prim vector field
        do l = 1, 3
            q_prim_vf(l + xibeg - 1)%sf(i, j, k) = xi_cart(l)
        end do

        ! Put your variable assignments here
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
