#:def Hardcoded1DVariables()
    ! Place any declaration of intermediate variables here
    real(wp) :: x_mid_diffu, width_sq, profile_shape, temp, molar_mass_inv, y1, y2, y3, y4
#:enddef

#:def Hardcoded1D()
    select case (patch_icpp(patch_id)%hcid)
    case (170)
        ! This hardcoded case can be used to start a simulation with initial conditions given from a known 1D profile (e.g. Cantera, SDtoolbox)
        @: HardcodedReadValues()

    case (180)
        ! This is patch is hard-coded for test suite optimization used in the
        ! 1D_shuoser cases: "patch_icpp(2)%alpha_rho(1)": "1 + 0.2*sin(5*x)"
        if (patch_id == 2) then
            q_prim_vf(contxb + 0)%sf(i, 0, 0) = 1 + 0.2*sin(5*x_cc(i))
        end if

    case (181)
        ! This is patch is hard-coded for test suite optimization used in the
        ! 1D_titarevtorro cases: "patch_icpp(2)%alpha_rho(1)": "1 + 0.1*sin(20*x*pi)"
        q_prim_vf(contxb + 0)%sf(i, 0, 0) = 1 + 0.1*sin(20*x_cc(i)*pi)

    case (182)
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

        molar_mass_inv = y1/31.998_wp + &
                         y2/18.01508_wp + &
                         y3/16.04256_wp + &
                         y4/28.0134_wp

        q_prim_vf(contxb)%sf(i, 0, 0) = 1.01325_wp*(10.0_wp)**5/(temp*8.3144626_wp*1000.0_wp*molar_mass_inv)

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
