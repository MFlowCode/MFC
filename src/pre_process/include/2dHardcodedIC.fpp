#:def Hardcoded2DVariables()

    real(kind(0d0)) :: eps
    real(kind(0d0)) :: r, rmax, gam, umax, p0

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

    case default
        if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
        end if
    end select

#:enddef
