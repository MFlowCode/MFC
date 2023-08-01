#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here
    

#:enddef

#:def Hardcoded3D()

    select case(patch_icpp(patch_id)%hcid)
        case(300)
            ! Put your variable assignments here
        case default
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
    end select

#:enddef
