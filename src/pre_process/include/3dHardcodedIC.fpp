#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here

    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

    real(wp) :: eps

    real(wp) :: pres
    real(wp), dimension(0:m_glb, 0:p_glb) :: ih
    integer :: i, j, pos, start, end
    character(len=10000) :: line
    character(len=25) :: value

    if (interface_file /= '.') then
        open(unit=10, file=trim(interface_file), status="old", action="read")
        do i = 0, m_glb
            read(10, '(A)') line  ! Read a full line as a string
            start = 1

            do j = 0, p_glb
                end = index(line(start:), ',')  ! Find the next comma
                if (end == 0) then
                    value = trim(adjustl(line(start:)))  ! Last value in the line
                else
                    value = trim(adjustl(line(start:start+end-2)))  ! Extract substring
                    start = start + end  ! Move to next value
                end if
                read(value, *) ih(i, j)  ! Convert string to numeric value
                if (.not. f_is_default(normMag)) ih(i,j )= ih(i,j) * normMag
                if (.not. f_is_default(normFac)) ih(i,j) = ih(i,j) / normFac
            end do
        end do
        close(10)

        print*, "Interface file "//trim(interface_file)//" read"
    end if

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

        alph = 5e-1_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

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

    case (302) ! (3D Perlin Noise Interface)

        alph = 0.5_wp * (1 + (1._wp - 2._wp * eps) * &
                    tanh((ih(start_idx(1) + i,start_idx(3) + k)  - y_cc(j))*100._wp))

        q_prim_vf(advxb)%sf(i,j,k) = alph
        q_prim_vf(advxe)%sf(i,j,k) = 1._wp - alph

        q_prim_vf(contxb)%sf(i,j,k) = q_prim_vf(advxb)%sf(i,j,k) * 1._wp
        q_prim_vf(contxe)%sf(i,j,k) = q_prim_vf(advxe)%sf(i,j,k) * (1._wp / 950._wp)

        q_prim_vf(E_idx)%sf(i,j,k) = p0 + &
            (q_prim_vf(contxb)%sf(i,j,k) + q_prim_vf(contxe)%sf(i,j,k)) * g0 * &
            (ih(start_idx(1) + i, start_idx(3) + k) - y_cc(j))

        if (surface_tension) q_prim_vf(c_idx)%sf(i,j,k) = alph

        ! Put your variable assignments here
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
