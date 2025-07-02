#:mute

#:def ASSERT_LIST(data, datatype)
    #:assert data is not None
    #:assert isinstance(data, list)
    #:assert len(data) != 0
    #:assert all(isinstance(element, datatype) for element in data)
#:enddef

#:def GEN_PARENTHESES_CLAUSE(clause_name, clause_str)
    #:set clause_regex = re.compile(',(?![^(]*\\))')
    #:assert isinstance(clause_name, str)
    #:if clause_str is not None
        #:set count = 0
        #:assert isinstance(clause_str, str)
        #:assert clause_str[0] == '[' and clause_str[-1] == ']'
        #:for c in clause_str
            #:if c == '('
                #:set count = count + 1
            #:elif c == ')'
                #:set count = count - 1
            #:endif
            #:if c == ',' and count > 1
                #:stop 'Nested parentheses with comma inside is not supported. Incorrect clause: {}'.format(clause_str)
            #:elif count < 0
                #:stop 'Missing parentheses. Incorrect clause: {}'.format(clause_str)
            #:endif
        #:endfor
        #:set clause_str = re.sub(clause_regex, ';', clause_str)
        #:set clause_list = [x.strip() for x in clause_str.strip('[]').split(';')]
        $:ASSERT_LIST(clause_list, str)
        #:set clause_str = clause_name + '(' + ', '.join(clause_list) + ') '
    #:else
        #:set clause_str = ''
    #:endif
    $:clause_str
#:enddef

#:def GEN_PRIVATE_STR(private, initialized_values)
    #:assert isinstance(initialized_values, bool)
    #:if initialized_values == True
        #:set private_val = GEN_PARENTHESES_CLAUSE('firstprivate', private)
    #:else
        #:set private_val = GEN_PARENTHESES_CLAUSE('private', private)
    #:endif
    $:private_val
#:enddef

#:def GEN_COPY_STR(copy)
    #:set copy_val = GEN_PARENTHESES_CLAUSE('copy', copy)
    $:copy_val
#:enddef

#:def GEN_COPYIN_STR(copyin, readonly)
    #:assert isinstance(readonly, bool)
    #:set copyin_val = GEN_PARENTHESES_CLAUSE('copyin', copyin)
    #:if copyin is not None and readonly == True
        #:set index = copyin_val.find('copyin(') + len('copyin(')
        #:set copyin_val = copyin_val[:index] + 'readonly:' + copyin_val[index:]
    #:endif
    $:copyin_val
#:enddef

#:def GEN_COPYOUT_STR(copyout)
    #:set copyout_val = GEN_PARENTHESES_CLAUSE('copyout', copyout)
    $:copyout_val
#:enddef

#:def GEN_CREATE_STR(create)
    #:set create_val = GEN_PARENTHESES_CLAUSE('create', create)
    $:create_val
#:enddef

#:def GEN_NOCREATE_STR(no_create)
    #:set nocreate_val = GEN_PARENTHESES_CLAUSE('no_create', no_create)
    $:nocreate_val
#:enddef

#:def GEN_DELETE_STR(delete)
    #:set delete_val = GEN_PARENTHESES_CLAUSE('delete', delete)
    $:delete_val
#:enddef

#:def GEN_PRESENT_STR(present)
    #:set present_val = GEN_PARENTHESES_CLAUSE('present', present)
    $:present_val
#:enddef

#:def GEN_DEVICEPTR_STR(deviceptr)
    #:set deviceptr_val = GEN_PARENTHESES_CLAUSE('deviceptr', deviceptr)
    $:deviceptr_val
#:enddef

#:def GEN_ATTACH_STR(attach)
    #:set attach_val = GEN_PARENTHESES_CLAUSE('attach', attach)
    $:attach_val
#:enddef

#:def GEN_DETACH_STR(detach)
    #:set detach_val = GEN_PARENTHESES_CLAUSE('detach', detach)
    $:detach_val
#:enddef

#:def GEN_LINK_STR(link)
    #:set link_val = GEN_PARENTHESES_CLAUSE('link', link)
    $:link_val
#:enddef

#:def GEN_EXTRA_ARGS_STR(extraArgs)
    #:if extraArgs is not None
        #:assert isinstance(extraArgs, str)
        #:set extraArgs_val = extraArgs
    #:else
        #:set extraArgs_val = ''
    #:endif
    $:extraArgs_val
#:enddef

#:def GEN_PARALLELISM_STR(parallelism)
    #:if parallelism is not None
        #:assert isinstance(parallelism, str)
        #:assert parallelism[0] == '[' and parallelism[-1] == ']'
        #:set parallelism_list = [x.strip() for x in parallelism.strip('[]').split(',')]
        $:ASSERT_LIST(parallelism_list, str)
        #:assert all((element == 'gang' or element == 'worker' or &
            & element == 'vector' or element == 'seq') for element in parallelism_list)
        #:set parallelism_val = ' '.join(parallelism_list) + ' '
    #:else
        #:set parallelism_val = ''
    #:endif
    $:parallelism_val
#:enddef

#:def GEN_COLLAPSE_STR(collapse)
    #:if collapse is not None
        #:set collapse = int(collapse)
        #:assert isinstance(collapse, int)
        #:assert collapse > 1
        #:set collapse_val = 'collapse(' + str(collapse) + ') '
    #:else
        #:set collapse_val = ''
    #:endif
    $:collapse_val
#:enddef

#:def GEN_DEFAULT_STR(default)
    #:if default is not None
        #:assert isinstance(default, str)
        #:assert (default == 'present' or default == 'none')
        #:set default_val = 'default(' + default + ') '
    #:else
        #:set default_val = ''
    #:endif
    $:default_val
#:enddef

#:def GEN_REDUCTION_STR(reduction, reductionOp)
    #:if reduction is not None and reductionOp is not None
        #:assert isinstance(reduction, str)
        #:assert isinstance(reductionOp, str)
        #:assert reduction[0] == '[' and reduction[-1] == ']'
        #:assert reductionOp[0] == '[' and reductionOp[-1] == ']'
        #:set reduction = reduction.replace(' ', '')
        #:set reduction = reduction[1:-1]
        #:set reduction_list = reduction.split('],')
        #:set reduction_list = [str + ']' for str in reduction_list]
        #:assert all(str[0] == '[' and str[-1] == ']' for str in reduction_list)

        #:set reductionOp_list = [x.strip() for x in reductionOp.strip('[]').split(',')]
        $:ASSERT_LIST(reduction_list, str)
        $:ASSERT_LIST(reductionOp_list, str)
        #:assert len(reduction_list) == len(reductionOp_list)
        #:set reduction_val = ''
        #:for i in range(len(reduction_list))
            #:set temp_clause = GEN_PARENTHESES_CLAUSE('reduction', reduction_list[i]).strip('\n')
            #:set ind = temp_clause.find('reduction(') + len('reduction(')
            #:set reduction_val = reduction_val.strip('\n') + temp_clause[:ind] + reductionOp_list[i] + ':' + temp_clause[ind:]
        #:endfor
    #:elif reduction is not None or reductionOp is not None
        #:stop 'Cannot set the reduction list or reduction operation without setting the other'
    #:else
        #:set reduction_val = ''
    #:endif
    $:reduction_val
#:enddef

#:def GEN_HOST_STR(host)
    #:set host_val = GEN_PARENTHESES_CLAUSE('host', host)
    $:host_val
#:enddef

#:def GEN_DEVICE_STR(device)
    #:set device_val = GEN_PARENTHESES_CLAUSE('device', device)
    $:device_val
#:enddef

#:def GEN_USE_DEVICE_STR(use_device)
    #:set use_device_val = GEN_PARENTHESES_CLAUSE('use_device', use_device)
    $:use_device_val
#:enddef

#:def GPU_PARALLEL(code, private=None, default='present', firstprivate=None, reduction=None, reductionOp=None, &
    & copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, &
    & no_create=None, present=None, deviceptr=None, attach=None, extraAccArgs=None)
    #:set default_val = GEN_DEFAULT_STR(default)
    #:set private_val = GEN_PRIVATE_STR(private, False).strip('\n') + GEN_PRIVATE_STR(firstprivate, True).strip('\n')
    #:set reduction_val = GEN_REDUCTION_STR(reduction, reductionOp)
    #:set copy_val = GEN_COPY_STR(copy)
    #:set copyin_val = GEN_COPYIN_STR(copyin, False).strip('\n') + GEN_COPYIN_STR(copyinReadOnly, True).strip('\n')
    #:set copyout_val = GEN_COPYOUT_STR(copyout)
    #:set create_val = GEN_CREATE_STR(create)
    #:set no_create_val = GEN_NOCREATE_STR(no_create)
    #:set present_val = GEN_PRESENT_STR(present)
    #:set deviceptr_val = GEN_DEVICEPTR_STR(deviceptr)
    #:set attach_val = GEN_ATTACH_STR(attach)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = default_val.strip('\n') + private_val.strip('\n') + reduction_val.strip('\n') + &
        & copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + &
        & deviceptr_val.strip('\n') + attach_val.strip('\n')
    #:set acc_directive = '!$acc parallel ' + &
        & clause_val + extraAccArgs_val.strip('\n')
    #:set end_acc_directive = '!$acc end parallel'
    $:acc_directive
    $:code
    $:end_acc_directive
#:enddef


#:def GPU_PARALLEL_LOOP(collapse=None, private=None, parallelism='[gang, vector]', &
    & default='present', firstprivate=None, reduction=None, reductionOp=None, &
    & copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, &
    & no_create=None, present=None, deviceptr=None, attach=None, extraAccArgs=None)
    #:set collapse_val = GEN_COLLAPSE_STR(collapse)
    #:set parallelism_val = GEN_PARALLELISM_STR(parallelism)
    #:set default_val = GEN_DEFAULT_STR(default)
    #:set private_val = GEN_PRIVATE_STR(private, False).strip('\n') + GEN_PRIVATE_STR(firstprivate, True).strip('\n')
    #:set reduction_val = GEN_REDUCTION_STR(reduction, reductionOp)
    #:set copy_val = GEN_COPY_STR(copy)
    #:set copyin_val = GEN_COPYIN_STR(copyin, False).strip('\n') + GEN_COPYIN_STR(copyinReadOnly, True).strip('\n')
    #:set copyout_val = GEN_COPYOUT_STR(copyout)
    #:set create_val = GEN_CREATE_STR(create)
    #:set no_create_val = GEN_NOCREATE_STR(no_create)
    #:set present_val = GEN_PRESENT_STR(present)
    #:set deviceptr_val = GEN_DEVICEPTR_STR(deviceptr)
    #:set attach_val = GEN_ATTACH_STR(attach)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = collapse_val.strip('\n') + parallelism_val.strip('\n') + &
        & default_val.strip('\n') + private_val.strip('\n') + reduction_val.strip('\n') + &
        & copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + &
        & deviceptr_val.strip('\n') + attach_val.strip('\n')
    #:set acc_directive = '!$acc parallel loop ' + &
        & clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def GPU_ROUTINE(function_name=None, parallelism=None, nohost=False, cray_inline=False, extraAccArgs=None)
    #:assert isinstance(cray_inline, bool)
    #:set parallelism_val = GEN_PARALLELISM_STR(parallelism)
    #:assert isinstance(nohost, bool)
    #:if nohost == True
        #:set nohost_val = 'nohost'
    #:else
        #:set nohost_val = ''
    #:endif
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = parallelism_val.strip('\n') + nohost_val.strip('\n')
    #:set acc_directive = '!$acc routine ' + &
        & clause_val + extraAccArgs_val.strip('\n')
    #:if cray_inline == True
        #:if not isinstance(function_name, str)
            #:stop "When inlining for Cray Compiler, function name must be given and given as a string"
        #:endif
        #:set cray_directive = ('!DIR$ INLINEALWAYS ' + function_name).strip('\n')
#ifdef _CRAYFTN
        $:cray_directive
#else
        $:acc_directive
#endif
    #:else
        $:acc_directive
    #:endif
#:enddef

#:def GPU_DECLARE(copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, present=None, deviceptr=None, link=None, extraAccArgs=None)
    #:set copy_val = GEN_COPY_STR(copy)
    #:set copyin_val = GEN_COPYIN_STR(copyin, False).strip('\n') + GEN_COPYIN_STR(copyinReadOnly, True).strip('\n')
    #:set copyout_val = GEN_COPYOUT_STR(copyout)
    #:set create_val = GEN_CREATE_STR(create)
    #:set present_val = GEN_PRESENT_STR(present)
    #:set deviceptr_val = GEN_DEVICEPTR_STR(deviceptr)
    #:set link_val = GEN_LINK_STR(link)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & present_val.strip('\n') + deviceptr_val.strip('\n') + &
        & link_val.strip('\n')
    #:set acc_directive = '!$acc declare ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def GPU_LOOP(collapse=None, parallelism=None, data_dependency=None, reduction=None, reductionOp=None, private=None, extraAccArgs=None)
    #:set collapse_val = GEN_COLLAPSE_STR(collapse)
    #:set parallelism_val = GEN_PARALLELISM_STR(parallelism)
    #:if data_dependency is not None
        #:assert isinstance(data_dependency, str)
        #:assert (data_dependency == 'auto' or data_dependency == 'independent')
        #:set data_dependency_val = data_dependency 
    #:else
        #:set data_dependency_val = ''
    #:endif
    #:set private_val = GEN_PRIVATE_STR(private, False)
    #:set reduction_val = GEN_REDUCTION_STR(reduction, reductionOp)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = collapse_val.strip('\n') + parallelism_val.strip('\n') + &
        & data_dependency_val.strip('\n') + private_val.strip('\n') + &
        & reduction_val.strip('\n')
    #:set acc_directive = '!$acc loop ' + &
        & clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def GPU_DATA(code, copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, no_create=None, present=None, deviceptr=None, attach=None, default=None, extraAccArgs=None)
    #:assert code is not None
    #:assert isinstance(code, str)
    #:if code == '' or code.isspace()
        #:stop 'GPU_DATA macro has no effect on the code as it is not surrounding any code'
    #:endif
    #:set copy_val = GEN_COPY_STR(copy)
    #:set copyin_val = GEN_COPYIN_STR(copyin, False).strip('\n') + GEN_COPYIN_STR(copyinReadOnly, True).strip('\n')
    #:set copyout_val = GEN_COPYOUT_STR(copyout)
    #:set create_val = GEN_CREATE_STR(create)
    #:set no_create_val = GEN_NOCREATE_STR(no_create)
    #:set present_val = GEN_PRESENT_STR(present)
    #:set deviceptr_val = GEN_DEVICEPTR_STR(deviceptr)
    #:set attach_val = GEN_ATTACH_STR(attach)
    #:set default_val = GEN_DEFAULT_STR(default)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + & 
        & deviceptr_val.strip('\n') + attach_val.strip('\n') + &
        & default_val.strip('\n')
    #:set acc_directive = '!$acc data ' + clause_val + extraAccArgs_val.strip('\n')
    #:set end_acc_directive = '!$acc end data'
    $:acc_directive
    $:code
    $:end_acc_directive
#:enddef

#:def GPU_HOST_DATA(code, use_device=None, extraAccArgs=None)
    #:assert code is not None
    #:assert isinstance(code, str)
    #:if code == '' or code.isspace()
        #:stop 'GPU_HOST_DATA macro has no effect on the code as it is not surrounding any code'
    #:endif
    #:set use_device_val = GEN_USE_DEVICE_STR(use_device)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = use_device_val.strip('\n')
    #:set acc_directive = '!$acc host_data ' + clause_val + extraAccArgs_val.strip('\n')
    #:set end_acc_directive = '!$acc end host_data'
    $:acc_directive
    $:code
    $:end_acc_directive
#:enddef

#:def GPU_ENTER_DATA(copyin=None, copyinReadOnly=None, create=None, attach=None, extraAccArgs=None)
    #:set copyin_val = GEN_COPYIN_STR(copyin, False).strip('\n') + GEN_COPYIN_STR(copyinReadOnly, True).strip('\n')
    #:set create_val = GEN_CREATE_STR(create)
    #:set attach_val = GEN_ATTACH_STR(attach)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = copyin_val.strip('\n') + create_val.strip('\n') + attach_val.strip('\n')
    #:set acc_directive = '!$acc enter data ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def GPU_EXIT_DATA(copyout=None, delete=None, detach=None, extraAccArgs=None)
    #:set copyout_val = GEN_COPYOUT_STR(copyout)
    #:set delete_val = GEN_DELETE_STR(delete)
    #:set detach_val = GEN_DETACH_STR(detach)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = copyout_val.strip('\n') + delete_val.strip('\n') + detach_val.strip('\n')
    #:set acc_directive = '!$acc exit data ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def GPU_CACHE(cache, extraAccArgs=None)
    #:set cache_val = GEN_PARENTHESES_CLAUSE('cache', cache)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = cache_val.strip('\n')
    #:set acc_directive = '!$acc ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def GPU_ATOMIC(atomic, extraAccArgs=None)
    #:assert isinstance(atomic, str)
    #:assert (atomic == 'read' or atomic == 'write' or atomic == 'update' or atomic == 'capture')
    #:set atomic_val = atomic
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = atomic_val.strip('\n')
    #:set acc_directive = '!$acc atomic ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def GPU_UPDATE(host=None, device=None, extraAccArgs=None)
    #:set host_val = GEN_HOST_STR(host)
    #:set device_val = GEN_DEVICE_STR(device)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = host_val.strip('\n') + device_val.strip('\n')
    #:set acc_directive = '!$acc update ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def GPU_WAIT(extraAccArgs=None)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = ''
    #:set acc_directive = '!$acc wait ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:endmute
! New line at end of file is required for FYPP
