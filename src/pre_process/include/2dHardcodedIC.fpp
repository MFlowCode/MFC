#:def Hardcoded2DVariables()

    real(kind(0d0)) :: eps

    eps = 1e-9

#:enddef

#:def Hardcoded2D()
    
    select case(patch_icpp(patch_id)%hcid) ! 2D_hardcoded_ic example case
        case(200)
            if (y_cc(j) <= (-x_cc(i)**3 + 1)**(1d0/3d0)) then
                ! Volume Fractions
                q_prim_vf(advxb)%sf(i, j, 0) = eps
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0-eps
                ! Denssities
                q_prim_vf(contxb)%sf(i, j, 0) = eps*1000d0
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0-eps)*1d0
                ! Pressure
                q_prim_vf(E_idx)%sf(i, j, 0) = 1000d0
            end if
        case default
            if (proc_rank == 0) then
                call s_int_to_str(patch_id, iStr)
                call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
            end if
    end select

#:enddef
