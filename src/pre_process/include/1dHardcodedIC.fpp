#:def Hardcoded1DVariables()
    ! Place any declaration of intermediate variables here
#:enddef

#:def Hardcoded1D()

    select case (patch_icpp(patch_id)%hcid)
    case (170)
        ! This hardcoded case can be used to start a simulation with initial conditions given from a known 1D profile (e.g. Cantera, SDtoolbox)
        @: HardcodedReadValues()

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
