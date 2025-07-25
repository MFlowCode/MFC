#:include 'shared_parallel_macros.fpp'

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
            #! #:set default_val = 'defaultmap(present:aggregate) defaultmap(present:allocatable) defaultmap(present:pointer) '
            #:set default_val = 'defaultmap(tofrom:aggregate) defaultmap(tofrom:allocatable) defaultmap(tofrom:pointer) '
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
    #:if no_create is not None
        #:stop 'no_create is not supported yet'
    #:endif
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

#:def OMP_ATTACH_STR(attach)
    #! #:if attach is not None
        #! #:stop 'attach is not supported yet'
    #! #:endif
    #:set attach_val = ''
    $:attach_val
#:enddef

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
    #:set attach_val = OMP_MAP_STR('tofrom', attach)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set omp_clause_val = default_val.strip('\n') + private_val.strip('\n') + reduction_val.strip('\n') + &
        & copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + &
        & deviceptr_val.strip('\n') + attach_val.strip('\n')
    
    #:set omp_clause_val = omp_clause_val.strip('\n')
    #:set omp_directive = '!$omp target teams ' + omp_clause_val + extraOmpArgs_val.strip('\n')

    #:set omp_end_directive = '!$omp end target teams'
    $:omp_directive
    $:code
    $:omp_end_directive
#:enddef

#:def OMP_PARALLEL_LOOP(code, collapse=None, private=None, parallelism='[gang, vector]', &
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
    #:set attach_val = OMP_MAP_STR('tofrom', attach)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = collapse_val.strip('\n') + parallelism_val.strip('\n') + &
        & default_val.strip('\n') + private_val.strip('\n') + reduction_val.strip('\n') + &
        & copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + &
        & deviceptr_val.strip('\n') + attach_val.strip('\n')
    #! Hardcoding the parallelism for now
    #:set omp_directive = '!$omp target teams loop ' + &
        & clause_val + extraOmpArgs_val.strip('\n')
    #:set omp_end_directive = '!$omp end target teams loop'
    $:omp_directive
    $:code
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
    #:set clause_val = nohost_val.strip('\n')
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

#! Not implemented yet
#:def OMP_LOOP(collapse=None, parallelism=None, data_dependency=None, reduction=None, reductionOp=None, private=None, extraOmpArgs=None)
    #! loop is going to be ignored since all loops right now are seq
    #:set temp = '!$omp loop bind(thread)'
    $:temp
#:enddef

#:def OMP_DATA(code, copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, no_create=None, present=None, deviceptr=None, attach=None, default=None, extraOmpArgs=None)
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
    #:set attach_val = OMP_MAP_STR('tofrom', attach)
    #:set default_val = OMP_DEFAULT_STR(default)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set clause_val = copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + & 
        & deviceptr_val.strip('\n') + attach_val.strip('\n') + &
        & default_val.strip('\n')
    #:set omp_directive = '!$omp target data ' + clause_val + extraOmpArgs_val.strip('\n')
    #:set end_omp_directive = '!$omp end target data'
    $:omp_directive
    $:code
    $:end_omp_directive
#:enddef

#:def OMP_ENTER_DATA(copyin=None, copyinReadOnly=None, create=None, attach=None, extraOmpArgs=None)
    #:set copyin_val = OMP_COPYIN_STR(copyin).strip('\n') + OMP_COPYIN_STR(copyinReadOnly).strip('\n')
    #:set create_val = OMP_CREATE_STR(create)
    #:set attach_val = OMP_MAP_STR('to', attach)
    #:set extraOmpArgs_val = GEN_EXTRA_ARGS_STR(extraOmpArgs)
    #:set omp_clause_val = copyin_val.strip('\n') + create_val.strip('\n') + attach_val.strip('\n')
    #:set omp_directive = '!$omp target enter data ' + omp_clause_val + extraOmpArgs_val.strip('\n')
    $:omp_directive
#:enddef

#:def OMP_EXIT_DATA(copyout=None, delete=None, detach=None, extraOmpArgs=None)
    #:set copyout_val = OMP_COPYOUT_STR(copyout)
    #:set delete_val = OMP_DELETE_STR(delete)
    #:set detach_val = OMP_MAP_STR('from', detach)
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
! New line at end of file is required for FYPP