#:def Hardcoded1DVariables()
    ! Place any declaration of intermediate variables here
    real(wp) :: x_mid_diffu, width_sq, profile_shape, temp, molar_mass_inv, y1, y2, y3, y4
#:enddef

#:def Hardcoded1D()
    select case (patch_icpp(patch_id)%hcid)
    case (150)  ! 1D Smooth Alfven Case for MHD
        ! velocity
        q_prim_vf(momxb + 1)%sf(i, 0, 0) = 0.1_wp*sin(2._wp*pi*x_cc(i))
        q_prim_vf(momxb + 2)%sf(i, 0, 0) = 0.1_wp*cos(2._wp*pi*x_cc(i))

        ! magnetic field
        q_prim_vf(B_idx%end - 1)%sf(i, 0, 0) = 0.1_wp*sin(2._wp*pi*x_cc(i))
        q_prim_vf(B_idx%end)%sf(i, 0, 0) = 0.1_wp*cos(2._wp*pi*x_cc(i))
    case (170)  ! 1D profile from external data (e.g. Cantera, SDtoolbox)
        ! This hardcoded case can be used to start a simulation with initial conditions given from a known 1D profile (e.g. Cantera,
        ! SDtoolbox)
        @: HardcodedReadValues()
    case (180)  ! Shu-Osher problem
        ! This is patch is hard-coded for test suite optimization used in the 1D_shuoser cases: "patch_icpp(2)%alpha_rho(1)": "1 +
        ! 0.2*sin(5*x)"
        if (patch_id == 2) then
            q_prim_vf(contxb + 0)%sf(i, 0, 0) = 1 + 0.2*sin(5*x_cc(i))
        end if
    case (181)  ! Titarev-Torro problem
        ! This is patch is hard-coded for test suite optimization used in the 1D_titarevtorro cases: "patch_icpp(2)%alpha_rho(1)":
        ! "1 + 0.1*sin(20*x*pi)"
        q_prim_vf(contxb + 0)%sf(i, 0, 0) = 1 + 0.1*sin(20*x_cc(i)*pi)
    case (182)  ! Multi-component diffusion
        ! This patch is a hard-coded for test suite optimization (multiple component diffusion)
        x_mid_diffu = 0.05_wp/2.0_wp
        width_sq = (2.5_wp*10.0_wp**(-3.0_wp))**2
        profile_shape = 1.0_wp - 0.5_wp*exp(-(x_cc(i) - x_mid_diffu)**2/width_sq)
        q_prim_vf(momxb)%sf(i, 0, 0) = 0.0_wp
        q_prim_vf(E_idx)%sf(i, 0, 0) = 1.01325_wp*(10.0_wp)**5
        q_prim_vf(advxb)%sf(i, 0, 0) = 1.0_wp

        y1 = (0.195_wp - 0.142_wp)*profile_shape + 0.142_wp
        y2 = (0.0_wp - 0.1_wp)*profile_shape + 0.1_wp
        y3 = (0.214_wp - 0.0_wp)*profile_shape + 0.0_wp
        y4 = (0.591_wp - 0.758_wp)*profile_shape + 0.758_wp

        q_prim_vf(chemxb)%sf(i, 0, 0) = y1
        q_prim_vf(chemxb + 1)%sf(i, 0, 0) = y2
        q_prim_vf(chemxb + 2)%sf(i, 0, 0) = y3
        q_prim_vf(chemxb + 3)%sf(i, 0, 0) = y4

        temp = (320.0_wp - 1350.0_wp)*profile_shape + 1350.0_wp

        molar_mass_inv = y1/31.998_wp + y2/18.01508_wp + y3/16.04256_wp + y4/28.0134_wp

        q_prim_vf(contxb)%sf(i, 0, 0) = 1.01325_wp*(10.0_wp)**5/(temp*8.3144626_wp*1000.0_wp*molar_mass_inv)

    case(191)
    ! 1. Constant Pressure (1 atm = 101325 Pa)
        q_prim_vf(E_idx)%sf(i, 0, 0) = 101325.0_wp

        ! 2. Zero Velocity everywhere (1D means only momxb is active)
        q_prim_vf(momxb)%sf(i, 0, 0) = 0.0_wp

        ! 3. Pure H2 Mass Fraction (H2 is the first species in h2o2.yaml)
        q_prim_vf(chemxb)%sf(i, 0, 0) = 1.0_wp

        ! 4. Piecewise "Hat" Temperature profile
        ! Domain is 0.0 to 0.05. Midpoint is 0.025.
        if (x_cc(i) <= 0.025_wp) then
            ! Left half: 600 K at x=0.0 up to 1000 K at x=0.025
            temp = 700.0_wp + ((1000.0_wp - 700.0_wp) / 0.025_wp) * x_cc(i)
        else
            ! Right half: 1000 K at x=0.025 down to 900 K at x=0.05
            temp = 1200.0_wp + ((900.0_wp - 1000.0_wp) / 0.025_wp) * (x_cc(i) - 0.025_wp)
        end if

        ! 5. Pure H2 Molar Mass Inverse (~2.016 g/mol)
        molar_mass_inv = 1.0_wp / 2.01588_wp

        ! 6. Calculate Ideal Gas Density: rho = P / (R_spec * T)
        ! R_univ = 8.3144626 J/(mol K). Factor of 1000 converts g/mol to kg/mol.
        q_prim_vf(contxb)%sf(i, 0, 0) = 101325.0_wp / (temp * 8.3144626_wp * 1000.0_wp * molar_mass_inv)

   case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
    end select
#:enddef
