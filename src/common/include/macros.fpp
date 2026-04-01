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

! Caution: This macro requires the use of a binding script to set CUDA_VISIBLE_DEVICES, such that we have one GPU device per MPI
! rank. That's because for both cudaMemAdvise (preferred location) and cudaMemPrefetchAsync we use location = device_id = 0. For an
! example see misc/nvidia_uvm/bind.sh. NVIDIA unified memory page placement hint
#:def PREFER_GPU(*args)
#ifdef MFC_SIMULATION
#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
    block
        ! NVIDIA CUDA Fortran 25.3+: uses submodules (cuda_runtime_api, gpu_reductions, sort) See
        ! https://docs.nvidia.com/hpc-sdk/compilers/cuda-fortran-prog-guide/index.html#fortran-host-modules
#if __NVCOMPILER_MAJOR__ < 25 || (__NVCOMPILER_MAJOR__ == 25 && __NVCOMPILER_MINOR__ < 3)
        use cudafor, gpu_sum => sum, gpu_maxval => maxval, gpu_minval => minval
#else
        use cuda_runtime_api
#endif
        integer :: istat

        if (nv_uvm_pref_gpu) then
            #:for arg in args
                ! print*, "Moving ${arg}$ to GPU => ", SHAPE(${arg}$) set preferred location GPU
                istat = cudaMemAdvise(c_devloc(${arg}$), SIZEOF(${arg}$), cudaMemAdviseSetPreferredLocation, 0)
                if (istat /= cudaSuccess) then
                    write (*, "('Error code: ',I0, ': ')") istat
                    ! write(*,*) cudaGetErrorString(istat)
                end if
                ! set accessed by CPU
                istat = cudaMemAdvise(c_devloc(${arg}$), SIZEOF(${arg}$), cudaMemAdviseSetAccessedBy, cudaCpuDeviceId)
                if (istat /= cudaSuccess) then
                    write (*, "('Error code: ',I0, ': ')") istat
                    ! write(*,*) cudaGetErrorString(istat)
                end if
                ! prefetch to GPU - physically populate memory pages
                istat = cudaMemPrefetchAsync(c_devloc(${arg}$), SIZEOF(${arg}$), 0, 0)
                if (istat /= cudaSuccess) then
                    write (*, "('Error code: ',I0, ': ')") istat
                    ! write(*,*) cudaGetErrorString(istat)
                end if
            #:endfor
        end if
    end block
#endif
#endif
#:enddef

! Allocate and create GPU device memory
#:def ALLOCATE(*args)
    @:LOG({'@:ALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    #:set allocated_variables = ', '.join(args)
    allocate (${allocated_variables}$)
    #:set cleaned = []
    #:for a in args
        #:set s = a.rstrip()
        #:if s.endswith(')')
            #:set rev = s[::-1]
            #:set pos = next(i for i, ch, d in ( (j, c, sum(1 if t==')' else -1 if t=='(' else 0 for t in rev[:j+1])) for j, &
                             & c in enumerate(rev) ) if ch == '(' and d == 0 )
            #:set s = s[:len(s)-1-pos]
        #:endif
        $:cleaned.append(s)
    #:endfor
    #:set joined = ', '.join(cleaned)
    $:GPU_ENTER_DATA(create='[' + joined + ']')
#:enddef ALLOCATE

! Free GPU device memory and deallocate
#:def DEALLOCATE(*args)
    @:LOG({'@:DEALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    #:set allocated_variables = ', '.join(args)
    $:GPU_EXIT_DATA(delete=('[' + allocated_variables + ']'))
    deallocate (${allocated_variables}$)
#:enddef DEALLOCATE

! Cray-specific GPU pointer setup for vector fields
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

! Cray-specific GPU pointer setup for scalar fields
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

! Cray-specific GPU pointer setup for acoustic source spatials
#:def ACC_SETUP_source_spatials(*args)
#ifdef _CRAYFTN
    block
        @:LOG({'@:ACC_SETUP_source_spatials(${', '.join(args)}$)'})

        #:for arg in args
            $:GPU_ENTER_DATA(copyin=('[' + arg + ']'))
            if (associated(${arg}$%coord)) then
                $:GPU_ENTER_DATA(copyin=('[' + arg + '%coord]'))
            end if
            if (associated(${arg}$%val)) then
                $:GPU_ENTER_DATA(copyin=('[' + arg + '%val]'))
            end if
            if (associated(${arg}$%angle)) then
                $:GPU_ENTER_DATA(copyin=('[' + arg + '%angle]'))
            end if
            if (associated(${arg}$%xyz_to_r_ratios)) then
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
        call s_mpi_abort("${_FILE_.split('/')[-1]}$:${_LINE_}$: " // "Assertion failed: ${predicate}$. " &
                         & // ${message or '"No error description."'}$)
    end if
#:enddef
! New line at end of file is required for FYPP
