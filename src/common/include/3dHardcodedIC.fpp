#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph, Mach
    real(wp) :: eps

    ! IGR Jets
    ! Arrays to stor position and radii of jets from input file
    real(wp), dimension(:), allocatable :: y_th_arr, z_th_arr, r_th_arr
    ! Variables to describe initial condition of jet
    real(wp) :: r, ux_th, ux_am, p_th, p_am, rho_th, rho_am, y_th, z_th, r_th, eps_smooth
    real(wp) :: rcut, xcut ! Intermediate variables for creating smooth initial condition

    real(wp), dimension(0:n, 0:p) :: rcut_arr
    integer :: l, q, s ! Iterators for reading input files
    integer :: start, end ! Ints to keep track of position in file
    character(len=1000) :: line ! String to store line in ile
    character(len=25) :: value ! String to store value in line
    integer :: NJet ! Number of jets

    eps = 1e-9_wp

    if (patch_icpp(patch_id)%hcid == 303) then
        eps_smooth = 3._wp
        open (unit=10, file="njet.txt", status="old", action="read")
        read (10, *) NJet
        close (10)

        allocate (y_th_arr(0:NJet - 1))
        allocate (z_th_arr(0:NJet - 1))
        allocate (r_th_arr(0:NJet - 1))

        open (unit=10, file="jets.csv", status="old", action="read")
        do q = 0, NJet - 1
            read (10, '(A)') line  ! Read a full line as a string
            start = 1

            do l = 0, 2
                end = index(line(start:), ',')  ! Find the next comma
                if (end == 0) then
                    value = trim(adjustl(line(start:)))  ! Last value in the line
                else
                    value = trim(adjustl(line(start:start + end - 2)))  ! Extract substring
                    start = start + end  ! Move to next value
                end if
                if (l == 0) then
                    read (value, *) y_th_arr(q)  ! Convert string to numeric value
                elseif (l == 1) then
                    read (value, *) z_th_arr(q)
                else
                    read (value, *) r_th_arr(q)
                end if
            end do
        end do
        close (10)

        do q = 0, p
            do l = 0, n
                rcut = 0._wp
                do s = 0, NJet - 1
                    r = sqrt((y_cc(l) - y_th_arr(s))**2._wp + (z_cc(q) - z_th_arr(s))**2._wp)
                    rcut = rcut + f_cut_on(r - r_th_arr(s), eps_smooth)
                end do
                rcut_arr(l, q) = rcut
            end do
        end do
    end if

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
        rcut = f_cut_on(r - r_th, eps_smooth)
        xcut = f_cut_on(x_cc(i), eps_smooth)

        q_prim_vf(momxb)%sf(i, j, k) = ux_th*rcut*xcut + ux_am
        q_prim_vf(momxb + 1)%sf(i, j, k) = 0._wp
        q_prim_vf(momxe)%sf(i, j, k) = 0._wp

        if (num_fluids == 1) then
            q_prim_vf(contxb)%sf(i, j, k) = (rho_th - rho_am)*rcut*xcut + rho_am
        else
            q_prim_vf(advxb)%sf(i, j, k) = (1._wp - 2._wp*eps)*rcut*xcut + eps
            q_prim_vf(contxb)%sf(i, j, k) = rho_th*q_prim_vf(advxb)%sf(i, j, k)
            q_prim_vf(contxe)%sf(i, j, k) = rho_am*(1._wp - q_prim_vf(advxb)%sf(i, j, k))
        end if

        q_prim_vf(E_idx)%sf(i, j, k) = p_th*rcut*xcut + p_am

    case (303) ! 3D Multijet

        eps_smooth = 3.0_wp
        ux_th = 10*sqrt(1.4*0.4)
        ux_am = 2.5*sqrt(1.4*0.4)
        p_th = 0.8_wp
        p_am = 0.4_wp
        rho_th = 1._wp
        rho_am = 1._wp
        eps = 1e-6

        rcut = rcut_arr(j, k)
        xcut = f_cut_on(x_cc(i), eps_smooth)

        q_prim_vf(momxb)%sf(i, j, k) = ux_th*rcut*xcut + ux_am
        q_prim_vf(momxb + 1)%sf(i, j, k) = 0._wp
        q_prim_vf(momxe)%sf(i, j, k) = 0._wp

        if (num_fluids == 1) then
            q_prim_vf(contxb)%sf(i, j, k) = (rho_th - rho_am)*rcut*xcut + rho_am
        else
            q_prim_vf(advxb)%sf(i, j, k) = (1._wp - 2._wp*eps)*rcut*xcut + eps
            q_prim_vf(contxb)%sf(i, j, k) = rho_th*q_prim_vf(advxb)%sf(i, j, k)
            q_prim_vf(contxe)%sf(i, j, k) = rho_am*(1._wp - q_prim_vf(advxb)%sf(i, j, k))
        end if

        q_prim_vf(E_idx)%sf(i, j, k) = p_th*rcut*xcut + p_am

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
