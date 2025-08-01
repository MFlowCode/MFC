#:include 'parallel_macros.fpp'

#:def LOG(expr)
#ifdef MFC_DEBUG
    block
        use iso_fortran_env, only: output_unit

        print *, '${_FILE_.split('/')[-1]}$:${_LINE_}$: ', ${expr}$

        call flush (output_unit)
    end block
#endif
#:enddef

#:def PREFER_GPU(*args)
#ifdef MFC_SIMULATION
#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
    block
    use cudafor, gpu_sum => sum, gpu_maxval => maxval, gpu_minval => minval
    integer :: istat
    integer :: prefer_gpu_mode
    character(len=10) :: prefer_gpu_mode_str

    ! environment variable
    call get_environment_variable("NVIDIA_MANUAL_GPU_HINTS", prefer_gpu_mode_str)
    if (trim(prefer_gpu_mode_str) == "0") then ! OFF
        prefer_gpu_mode = 0
    elseif (trim(prefer_gpu_mode_str) == "1") then ! ON
        prefer_gpu_mode = 1
    else ! default
        prefer_gpu_mode = 0
    endif

    if (prefer_gpu_mode .eq. 1) then
    #:for arg in args
        !print*, "Moving ${arg}$ to GPU => ", SHAPE(${arg}$)
        ! unset
        istat = cudaMemAdvise( c_devloc(${arg}$), SIZEOF(${arg}$), cudaMemAdviseUnSetPreferredLocation, cudaCpuDeviceId )
        if (istat /= cudaSuccess) then
            write(*,"('Error code: ',I0, ': ')") istat
            write(*,*) cudaGetErrorString(istat)
        endif
        ! set
        istat = cudaMemAdvise( c_devloc(${arg}$), SIZEOF(${arg}$), cudaMemAdviseSetPreferredLocation, 0 )
        if (istat /= cudaSuccess) then
            write(*,"('Error code: ',I0, ': ')") istat
            write(*,*) cudaGetErrorString(istat)
        endif
    #:endfor
    end if
    end block
#endif
#endif
#:enddef

#:def ALLOCATE(*args)
    @:LOG({'@:ALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    #:set allocated_variables = ', '.join(args)
    allocate (${allocated_variables}$)
    $:GPU_ENTER_DATA(create=('[' + allocated_variables + ']'))
#:enddef ALLOCATE

#:def DEALLOCATE(*args)
    @:LOG({'@:DEALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    #:set allocated_variables = ', '.join(args)
    $:GPU_EXIT_DATA(delete=('[' + allocated_variables + ']'))
    deallocate (${allocated_variables}$)
#:enddef DEALLOCATE

#:def ACC_SETUP_VFs(*args)
#ifdef _CRAYFTN
    block
        integer :: macros_setup_vfs_i

        @:LOG({'@:ACC_SETUP_VFs(${', '.join(args)}$)'})

        #:for arg in args
            $:GPU_ENTER_DATA(copyin=('[' + arg + ']'))
            $:GPU_ENTER_DATA(copyin=('[' + arg + '%vf]'))
            if (allocated(${arg}$%vf)) then
                do macros_setup_vfs_i = lbound(${arg}$%vf, 1), ubound(${arg}$%vf, 1)
                    if (associated(${arg}$%vf(macros_setup_vfs_i)%sf)) then
                        $:GPU_ENTER_DATA(copyin=('[' + arg + '%vf(macros_setup_vfs_i)]'))
                        $:GPU_ENTER_DATA(copyin=('[' + arg + '%vf(macros_setup_vfs_i)%sf]'))
                    end if
                end do
            end if
        #:endfor
    end block
#endif
#:enddef

#:def ACC_SETUP_SFs(*args)
#ifdef _CRAYFTN
    block

        @:LOG({'@:ACC_SETUP_SFs(${', '.join(args)}$)'})

        #:for arg in args
            $:GPU_ENTER_DATA(copyin=('[' + arg + ']'))
            if (associated(${arg}$%sf)) then
                $:GPU_ENTER_DATA(copyin=('[' + arg + '%sf]'))
            end if
        #:endfor
    end block
#endif
#:enddef

#:def ACC_SETUP_source_spatials(*args)
#ifdef _CRAYFTN
    block

        @:LOG({'@:ACC_SETUP_source_spatials(${', '.join(args)}$)'})

        #:for arg in args
            $:GPU_ENTER_DATA(copyin=('[' + arg + ']'))
            if (allocated(${arg}$%coord)) then
                $:GPU_ENTER_DATA(copyin=('[' + arg + '%coord]'))
            end if
            if (allocated(${arg}$%val)) then
                $:GPU_ENTER_DATA(copyin=('[' + arg + '%val]'))
            end if
            if (allocated(${arg}$%angle)) then
                $:GPU_ENTER_DATA(copyin=('[' + arg + '%angle]'))
            end if
            if (allocated(${arg}$%xyz_to_r_ratios)) then
                $:GPU_ENTER_DATA(copyin=('[' + arg + '%xyz_to_r_ratios]'))
            end if
        #:endfor
    end block
#endif
#:enddef

#:def PROHIBIT(condition, message = None)
    if (${condition}$) then
        call s_prohibit_abort("${condition}$", ${message or '""'}$)
    end if
#:enddef

#:def ASSERT(predicate, message = None)
    if (.not. (${predicate}$)) then
        call s_mpi_abort("${_FILE_.split('/')[-1]}$:${_LINE_}$: "// &
                         "Assertion failed: ${predicate}$. " &
                         //${message or '"No error description."'}$)
    end if
#:enddef
! New line at end of file is required for FYPP
