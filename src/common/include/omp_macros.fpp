#:include 'shared_parallel_macros.fpp'

#:set NVIDIA_COMPILER_ID="NVHPC"
#:set PGI_COMPILER_ID="PGI"
#:set INTEL_COMPILER_ID="Intel"
#:set CCE_COMPILER_ID="Cray"
#:set AMD_COMPILER_ID="LLVMFlang"

#:def OMP_MAP_STR(map_type, var_list)
    #:assert map_type is not None
    #:assert isinstance(map_type, str)
    #:if var_list is not None
        #:set map_clause = 'map(' + map_type + ':'
        #:set map_val = GEN_CLAUSE(map_clause, var_list)
    #:else
        #:set map_val = ''
    #:endif
    $:map_val
#:enddef

#:def OMP_DEFAULT_STR(default)
    #:if default is not None
        #:assert isinstance(default, str)
        #:assert (default == 'present' or default == 'none')
        #:if default == 'present'
            #:if MFC_COMPILER == NVIDIA_COMPILER_ID or MFC_COMPILER == PGI_COMPILER_ID
                #:set default_val = 'defaultmap(tofrom:aggregate) defaultmap(tofrom:allocatable) defaultmap(tofrom:pointer) '
            #:elif MFC_COMPILER == CCE_COMPILER_ID
                #:set default_val = 'defaultmap(tofrom:aggregate) defaultmap(present:allocatable) defaultmap(present:pointer) '
            #:elif MFC_COMPILER == AMD_COMPILER_ID
                #:set default_val = ''
            #:else
                #:set default_val = 'defaultmap(tofrom:aggregate) defaultmap(tofrom:allocatable) defaultmap(tofrom:pointer) '
            #:endif
        #:elif default == 'none'
            #:stop 'Not Supported Yet'
        #:endif
    #:else
        #:set default_val = ''
    #:endif
    $:default_val
#:enddef

#:def OMP_COPY_STR(copy)
    #:set copy_val = OMP_MAP_STR('tofrom', copy)
    $:copy_val
#:enddef

#:def OMP_COPYIN_STR(copyin)
    #:set copyin_val = OMP_MAP_STR('to', copyin)
    $:copyin_val
#:enddef

#:def OMP_COPYOUT_STR(copyout)
    #:set copyout_val = OMP_MAP_STR('from', copyout)
    $:copyout_val
#:enddef

#:def OMP_CREATE_STR(create)
    #:set create_val = OMP_MAP_STR('always,alloc', create)
    $:create_val
#:enddef

#:def OMP_DELETE_STR(delete)
    #:set create_val = OMP_MAP_STR('release', delete)
    $:create_val
#:enddef

#:def OMP_NOCREATE_STR(no_create)
    #! OpenMP has no no_create equivalent, so this is a documented no-op: the variable is
    #! left to whatever mapping the enclosing region already applies to it. Note that is
    #! NOT equivalent to no_create on most targets -- OMP_DEFAULT_STR emits nothing unless
    #! the caller passes default='present', and even then only CCE maps allocatables and
    #! pointers as present; NVHPC/PGI and the fallback map them tofrom, which copies rather
    #! than reuses, and LLVMFlang emits nothing at all.
    #! Do NOT #:stop here: GPU_DATA expands both backends before #if selects one, so
    #! aborting would break OpenACC builds, where no_create is supported natively.
    #:set no_create_val = ''
    $:no_create_val
#:enddef

#:def OMP_PRESENT_STR(present)
    #:set present_val = OMP_MAP_STR('present,alloc', present)
    $:present_val
#:enddef

#:def OMP_DEVICEPTR_STR(deviceptr)
    #:set deviceptr_val = GEN_PARENTHESES_CLAUSE('use_device_ptr', deviceptr)
    $:deviceptr_val
#:enddef

#! #:def OMP_ATTACH_STR(attach)
#! #:set attach_val = OMP_MAP_STR('always,to', attach)
#! $:attach_val
#! #:enddef

#! #:def OMP_DETACH_STR(detach)
#! #:set detach_val = OMP_MAP_STR('always,from', detach)
#! $:detach_val
#! #:enddef

#:def OMP_TO_STR(to)
    #:set to_val = GEN_PARENTHESES_CLAUSE('to', to)
    $:to_val
#:enddef

#:def OMP_FROM_STR(to)
    #:set from_val = GEN_PARENTHESES_CLAUSE('from', to)
    $:from_val
#:enddef

#:def OMP_PARALLELISM_STR(parallelism)
    #:set temp = ''
    $:temp
#:enddef

#:def OMP_USE_DEVICE_ADDR_STR(use_device_addr)
    #:set use_device_addr_val = GEN_PARENTHESES_CLAUSE('use_device_addr', use_device_addr)
    $:use_device_addr_val
#:enddef

#:def OMP_USE_DEVICE_PTR_STR(use_device_ptr)
    #:set use_device_ptr_val = GEN_PARENTHESES_CLAUSE('use_device_ptr', use_device_ptr)
    $:use_device_ptr_val
#:enddef

#:def OMP_PARALLEL(code, private=None, default='present', firstprivate=None, reduction=None, reductionOp=None, &
    & copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, &
    & no_create=None, present=None, deviceptr=None, attach=None, extraOmpArgs=None)
    #:set default_val = OMP_DEFAULT_STR(default)
    #:set private_val = GEN_PRIVATE_STR(private, False).strip('\n') + GEN_PRIVATE_STR(firstprivate, True).strip('\n')
    #:set reduction_val = GEN_REDUCTION_STR(reduction, reductionOp)
    #:set copy_val = OMP_COPY_STR(copy)
    #:set copyin_val = OMP_COPYIN_STR(copyin).strip('\n') + OMP_COPYIN_STR(copyinReadOnly).strip('\n')
    #:set copyout_val = OMP_COPYOUT_STR(copyout)
    #:set create_val = OMP_CREATE_STR(create)
    #:set no_create_val = OMP_NOCREATE_STR(no_create)
    #:set present_val = OMP_PRESENT_STR(present)
    #:set deviceptr_val = OMP_DEVICEPTR_STR(deviceptr)
    #:set attach_val = OMP_MAP_STR('always,tofrom', attach)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set omp_clause_val = default_val.strip('\n') + private_val.strip('\n') + reduction_val.strip('\n') + &
        & copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + &
        & deviceptr_val.strip('\n') + attach_val.strip('\n')

    #:set omp_clause_val = omp_clause_val.strip('\n')
    #:set omp_directive = FOLD_DIRECTIVE('!$omp target teams ' + omp_clause_val + extraOmpArgs_val.strip('\n'), '!$omp').strip('\n')

    #:set omp_end_directive = '!$omp end target teams'
    $:omp_directive
    $:code
    $:omp_end_directive
#:enddef

#:def OMP_PARALLEL_LOOP(collapse=None, private=None, parallelism='[gang, vector]', &
    & default='present', firstprivate=None, reduction=None, reductionOp=None, &
    & copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, &
    & no_create=None, present=None, deviceptr=None, attach=None, extraOmpArgs=None)

    #:set collapse_val = GEN_COLLAPSE_STR(collapse)
    #:set parallelism_val = OMP_PARALLELISM_STR(parallelism)
    #:set default_val = OMP_DEFAULT_STR(default)
    #:set private_val = GEN_PRIVATE_STR(private, False).strip('\n') + GEN_PRIVATE_STR(firstprivate, True).strip('\n')
    #:set reduction_val = GEN_REDUCTION_STR(reduction, reductionOp)
    #:set copy_val = OMP_COPY_STR(copy)
    #:set copyin_val = OMP_COPYIN_STR(copyin).strip('\n') + OMP_COPYIN_STR(copyinReadOnly).strip('\n')
    #:set copyout_val = OMP_COPYOUT_STR(copyout)
    #:set create_val = OMP_CREATE_STR(create)
    #:set no_create_val = OMP_NOCREATE_STR(no_create)
    #:set present_val = OMP_PRESENT_STR(present)
    #:set deviceptr_val = OMP_DEVICEPTR_STR(deviceptr)
    #:set attach_val = OMP_MAP_STR('always,tofrom', attach)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = collapse_val.strip('\n') + parallelism_val.strip('\n') + &
        & default_val.strip('\n') + private_val.strip('\n') + reduction_val.strip('\n') + &
        & copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + &
        & deviceptr_val.strip('\n') + attach_val.strip('\n')
    #! Hardcoding the parallelism for now

    #:if MFC_COMPILER == NVIDIA_COMPILER_ID or MFC_COMPILER == PGI_COMPILER_ID
        #:set omp_start_directive = '!$omp target teams loop defaultmap(firstprivate:scalar) bind(teams,parallel) '
    #:elif MFC_COMPILER == CCE_COMPILER_ID
        #:set omp_start_directive = '!$omp target teams distribute parallel do defaultmap(firstprivate:scalar) '
    #:elif MFC_COMPILER == AMD_COMPILER_ID
        #:set omp_start_directive = '!$omp target teams distribute parallel do '
    #:else
        #:set omp_start_directive = '!$omp target teams loop defaultmap(firstprivate:scalar) bind(teams,parallel) '
    #:endif

    #:set omp_directive = FOLD_DIRECTIVE(omp_start_directive + clause_val + extraOmpArgs_val.strip('\n'), '!$omp').strip('\n')
    $:omp_directive
#:enddef

#:def END_OMP_PARALLEL_LOOP()
    #:if MFC_COMPILER == NVIDIA_COMPILER_ID or MFC_COMPILER == PGI_COMPILER_ID
        #:set omp_end_directive = '!$omp end target teams loop'
    #:elif MFC_COMPILER == CCE_COMPILER_ID
        #:set omp_end_directive = '!$omp end target teams distribute parallel do'
    #:elif MFC_COMPILER == AMD_COMPILER_ID
        #:set omp_end_directive = '!$omp end target teams distribute parallel do'
    #:else
        #:set omp_end_directive = '!$omp end target teams loop'
    #:endif

    $:omp_end_directive
#:enddef

#:def OMP_ROUTINE(function_name, nohost, extraOmpArgs)
    #:assert isinstance(nohost, bool)
    #:if nohost == True
        #:set nohost_val = 'device_type(nohost) '
    #:else
        #:set nohost_val = 'device_type(any) '
    #:endif
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:if function_name is not None
        #:set function_name_val = '(' + function_name + ') '
    #:else
        #:set function_name_val = ''
    #:endif

    #:if MFC_COMPILER == AMD_COMPILER_ID
        #:set clause_val = ''
    #:else
        #:set clause_val = nohost_val.strip('\n')
    #:endif

    #:set omp_directive = '!$omp declare target ' + &
        & clause_val + extraOmpArgs_val.strip('\n')
    $:omp_directive
#:enddef

#:def OMP_DECLARE(copyin=None, copyinReadOnly=None, create=None, link=None, extraOmpArgs=None)
    #:set copyin_val = OMP_TO_STR(copyin).strip('\n') + OMP_TO_STR(copyinReadOnly).strip('\n')
    #:set create_val = GEN_CLAUSE('(', create)
    #:set link_val = GEN_LINK_STR(link)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = copyin_val.strip('\n') + &
        & create_val.strip('\n') + link_val.strip('\n')
    #:set omp_directive = '!$omp declare target ' + clause_val + extraOmpArgs_val.strip('\n')
    $:omp_directive
#:enddef

#! Not fully implemented yet (ignores most args right now)
#:def OMP_LOOP(collapse=None, parallelism=None, data_dependency=None, reduction=None, reductionOp=None, private=None, &
               & extraOmpArgs=None)
    #:if MFC_COMPILER == NVIDIA_COMPILER_ID or MFC_COMPILER == PGI_COMPILER_ID
        #:set omp_directive = '!$omp loop bind(thread)'
    #:elif MFC_COMPILER == CCE_COMPILER_ID or MFC_COMPILER == AMD_COMPILER_ID
        #:set omp_directive = ''
    #:else
        #:set omp_directive = ''
    #:endif
    $:omp_directive
#:enddef

#:def OMP_DATA(code, copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, no_create=None, present=None, &
               & deviceptr=None, attach=None, default=None, extraOmpArgs=None)
    #:assert code is not None
    #:assert isinstance(code, str)
    #:if code == '' or code.isspace()
        #:stop 'GPU_DATA macro has no effect on the code as it is not surrounding any code'
    #:endif
    #:set copy_val = OMP_COPY_STR(copy)
    #:set copyin_val = OMP_COPYIN_STR(copyin).strip('\n') + OMP_COPYIN_STR(copyinReadOnly).strip('\n')
    #:set copyout_val = OMP_COPYOUT_STR(copyout)
    #:set create_val = OMP_CREATE_STR(create)
    #:set no_create_val = OMP_NOCREATE_STR(no_create)
    #:set present_val = OMP_PRESENT_STR(present)
    #:set deviceptr_val = OMP_DEVICEPTR_STR(deviceptr)
    #:set attach_val = OMP_MAP_STR('always,tofrom', attach)
    #:set default_val = OMP_DEFAULT_STR(default)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + &
        & deviceptr_val.strip('\n') + attach_val.strip('\n') + &
        & default_val.strip('\n')
    #! An OpenMP `target data` region must carry at least one map, use_device_ptr or
    #! use_device_addr clause. no_create contributes none (it is a no-op here), and
    #! `default` is None by default, so GPU_DATA(no_create=x) with no other clause left
    #! clause_val empty and emitted a bare `!$omp target data` that no compiler accepts.
    #! Emit the body alone in that case -- a data region with nothing to map has nothing
    #! to do, which is what GPU_DATA's own #else branch already does when neither backend
    #! is enabled. A #:stop is not an option here, for the reason in OMP_NOCREATE_STR.
    #! Only no_create earns that silence. If the caller asked for some other clause and
    #! the backend produced nothing for it -- default='present' on LLVMFlang, where
    #! OMP_DEFAULT_STR returns an empty string -- keep emitting the bare directive so the
    #! build fails loudly rather than silently discarding a region that was asked for.
    #:set other_clause_requested = copy is not None or copyin is not None or &
        & copyinReadOnly is not None or copyout is not None or create is not None or &
        & present is not None or deviceptr is not None or attach is not None or &
        & default is not None
    #:set emitted_clause = clause_val.strip() != '' or extraOmpArgs_val.strip() != ''
    #:set has_clauses = emitted_clause or other_clause_requested
    #:set omp_directive = '!$omp target data ' + clause_val + extraOmpArgs_val.strip('\n')
    #:set end_omp_directive = '!$omp end target data'
    #:if has_clauses
        $:omp_directive
    #:endif
    $:code
    #:if has_clauses
        $:end_omp_directive
    #:endif
#:enddef

#:def OMP_ENTER_DATA(copyin=None, copyinReadOnly=None, create=None, attach=None, extraOmpArgs=None)
    #:set copyin_val = OMP_COPYIN_STR(copyin).strip('\n') + OMP_COPYIN_STR(copyinReadOnly).strip('\n')
    #:set create_val = OMP_CREATE_STR(create)
    #:set attach_val = OMP_MAP_STR('always,to', attach)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set omp_clause_val = copyin_val.strip('\n') + create_val.strip('\n') + attach_val.strip('\n')
    #:set omp_directive = '!$omp target enter data ' + omp_clause_val + extraOmpArgs_val.strip('\n')
    $:omp_directive
#:enddef

#:def OMP_EXIT_DATA(copyout=None, delete=None, detach=None, extraOmpArgs=None)
    #:set copyout_val = OMP_COPYOUT_STR(copyout)
    #:set delete_val = OMP_DELETE_STR(delete)
    #:set detach_val = OMP_MAP_STR('always,from', detach)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = copyout_val.strip('\n') + delete_val.strip('\n') + detach_val.strip('\n')
    #:set omp_directive = '!$omp target exit data ' + clause_val + extraOmpArgs_val.strip('\n')
    $:omp_directive
#:enddef

#:def OMP_UPDATE(host=None, device=None, extraOmpArgs=None)
    #:set host_val = OMP_FROM_STR(host)
    #:set device_val = OMP_TO_STR(device)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = host_val.strip('\n') + device_val.strip('\n')
    #:set omp_directive = '!$omp target update ' + clause_val + extraOmpArgs_val.strip('\n')
    $:omp_directive
#:enddef

#:def OMP_HOST_DATA(code, use_device_addr, use_device_ptr, extraOmpArgs)
    #:assert code is not None
    #:assert isinstance(code, str)
    #:if code == '' or code.isspace()
        #:stop 'GPU_HOST_DATA macro has no effect on the code as it is not surrounding any code'
    #:endif
    #:set use_device_addr_val = OMP_USE_DEVICE_ADDR_STR(use_device_addr)
    #:set use_device_ptr_val = OMP_USE_DEVICE_PTR_STR(use_device_ptr)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = use_device_addr_val.strip('\n') + use_device_ptr_val.strip('\n')
    #:set omp_directive = '!$omp target data ' + clause_val + extraOmpArgs_val.strip('\n')
    #:set omp_end_directive = '!$omp end target data'
    $:omp_directive
    $:code
    $:omp_end_directive
#:enddef

#:def OMP_ATOMIC(atomic, extraOmpArgs=None)
    #:assert isinstance(atomic, str)
    #:assert (atomic == 'read' or atomic == 'write' or atomic == 'update' or atomic == 'capture')
    #:set atomic_val = atomic
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = atomic_val.strip('\n')
    #:set omp_directive = '!$omp atomic ' + clause_val + extraOmpArgs_val.strip('\n')
    $:omp_directive
#:enddef

#:def OMP_WAIT(extraOmpArgs=None)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = ''
    #:set omp_directive = '!$omp barrier ' + clause_val + extraOmpArgs_val.strip('\n')
    $:omp_directive
#:enddef

#:def UNDEF_AMD(code)
    #:if MFC_COMPILER != AMD_COMPILER_ID
        $:code
    #:endif
#:enddef

#:def DEF_AMD(code)
    #:if MFC_COMPILER == AMD_COMPILER_ID
        $:code
    #:endif
#:enddef

#:def UNDEF_CCE(code)
    #:if MFC_COMPILER != CCE_COMPILER_ID
        $:code
    #:endif
#:enddef

#:def DEF_CCE(code)
    #:if MFC_COMPILER == CCE_COMPILER_ID
        $:code
    #:endif
#:enddef

#:def UNDEF_NVIDIA(code)
    #:if MFC_COMPILER != NVIDIA_COMPILER_ID and MFC_COMPILER != PGI_COMPILER_ID
        $:code
    #:endif
#:enddef
! New line at end of file is required for FYPP
