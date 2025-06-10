#:def LOG(expr)
#ifdef MFC_DEBUG
    block
        use iso_fortran_env, only: output_unit

        print *, '${_FILE_.split('/')[-1]}$:${_LINE_}$: ', ${expr}$

        call flush (output_unit)
    end block
#endif
#:enddef

#:def ALLOCATE(*args)
    @:LOG({'@:ALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    allocate (${', '.join(args)}$)
    !$acc enter data create(${', '.join(args)}$)
#:enddef ALLOCATE

#:def DEALLOCATE(*args)
    @:LOG({'@:DEALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    deallocate (${', '.join(args)}$)
    !$acc exit data delete(${', '.join(args)}$)
#:enddef DEALLOCATE

#:def ACC_SETUP_VFs(*args)
#ifdef _CRAYFTN
    block
        integer :: macros_setup_vfs_i

        @:LOG({'@:ACC_SETUP_VFs(${', '.join(args)}$)'})

        #:for arg in args
            !$acc enter data copyin(${arg}$)
            !$acc enter data copyin(${arg}$%vf)
            if (allocated(${arg}$%vf)) then
                do macros_setup_vfs_i = lbound(${arg}$%vf, 1), ubound(${arg}$%vf, 1)
                    if (associated(${arg}$%vf(macros_setup_vfs_i)%sf)) then
                        !$acc enter data copyin(${arg}$%vf(macros_setup_vfs_i))
                        !$acc enter data create(${arg}$%vf(macros_setup_vfs_i)%sf)
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
            !$acc enter data copyin(${arg}$)
            if (associated(${arg}$%sf)) then
                !$acc enter data create(${arg}$%sf)
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
            !$acc enter data copyin(${arg}$)
            if (allocated(${arg}$%coord)) then
                !$acc enter data create(${arg}$%coord)
            end if
            if (allocated(${arg}$%val)) then
                !$acc enter data create(${arg}$%val)
            end if
            if (allocated(${arg}$%angle)) then
                !$acc enter data create(${arg}$%angle)
            end if
            if (allocated(${arg}$%xyz_to_r_ratios)) then
                !$acc enter data create(${arg}$%xyz_to_r_ratios)
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

#define t_vec3   real(wp), dimension(1:3)
#define t_mat4x4 real(wp), dimension(1:4,1:4)

#:def ASSERT(predicate, message = None)
    if (.not. (${predicate}$)) then
        call s_mpi_abort("${_FILE_.split('/')[-1]}$:${_LINE_}$: "// &
                         "Assertion failed: ${predicate}$. " &
                         //${message or '"No error description."'}$)
    end if
#:enddef

#:def ASSERT_LIST(data, datatype)
    #:assert data is not None
    #:assert isinstance(data, list)
    #:assert len(data) != 0
    #:assert all(isinstance(element, datatype) for element in data)
#:enddef

#:def GEN_PARENTHESES_CLAUSE(clause_name, clause_list)
    #:assert isinstance(clause_name, str)
    #:if clause_list is not None
        $:ASSERT_LIST(clause_list, str)
        #:set clause_str = clause_name + '(' + ', '.join(clause_list) + ') '
    #:else
        #:set clause_str = ''
    #:endif
    $:clause_str
#:enddef

#:def GEN_PRIVATE_STR(private, initalized_values)
    #:assert isinstance(initalized_values, bool)
    #:if initalized_values == True
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
    #:if copyin is not None
        $:ASSERT_LIST(copyin, str)
        #:if readonly == True
            #:set copyin_val = 'copyin(readonly:' + ', '.join(copyin) + ') '
        #:else
            #:set copyin_val = 'copyin(' + ', '.join(copyin) + ') '
        #:endif
    #:else
        #:set copyin_val = ''
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
        $:ASSERT_LIST(parallelism, str)
        #:assert all((element == 'gang' or element == 'worker' or &
            & element == 'vector' or element == 'seq') for element in parallelism)
        #:set parallelism_val = ' '.join(parallelism) + ' '
    #:else
        #:set parallelism_val = ''
    #:endif
    $:parallelism_val
#:enddef

#:def GEN_COLLAPSE_STR(collapse)
    #:if collapse is not None
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
        #:if isinstance(reduction, list) and isinstance(reductionOp, list)
            $:ASSERT_LIST(reduction, list)
            $:ASSERT_LIST(reductionOp, str)
            #:assert all(len(element) != 0 for element in reduction)
            #:assert all(isinstance(element, str) for sublist in reduction for element in sublist)
            #:assert len(reduction) == len(reductionOp)
            #:set reduction_list = ['reduction(' + op + ':' + ', '.join(red) + ') ' for (red, op) in zip(reduction, reductionOp)]
            #:set reduction_val = ' '. join(reduction_list) + ' '
        #:elif isinstance(reduction, list) and isinstance(reductionOp, str)
            $:ASSERT_LIST(reduction, str)
            #:assert isinstance(reductionOp, str)
            #:set reduction_val = 'reduction(' + reductionOp + ':' + ', '.join(reduction) + ') '
        #:else
            #:stop 'Invalid datatypes for reduction or reductionOp. Must be list of lists and lists or list and str respectively'
        #:endif
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


#:def PARALLEL_LOOP(collapse=None, private=None, parallelism=['gang', 'vector'], &
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

#:def ROUTINE(parallelism=['seq'], nohost=False, extraAccArgs=None)
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
    $:acc_directive
#:enddef

#:def DECLARE(copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, present=None, deviceptr=None, link=None, extraAccArgs=None)
    
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

#:def LOOP(collapse=None, parallelism=["seq"], data_dependency=None, reduction=None, reductionOp=None, private=None, extraAccArgs=None)

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

#:def DATA(copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, no_create=None, present=None, deviceptr=None, attach=None, default=None, extraAccArgs=None)
    
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
    $:acc_directive
#:enddef

#:def ENTER_DATA(copyin=None, copyinReadOnly=None, create=None, attach=None, extraAccArgs=None)
    #:set copyin_val = GEN_COPYIN_STR(copyin, False).strip('\n') + GEN_COPYIN_STR(copyinReadOnly, True).strip('\n')

    #:set create_val = GEN_CREATE_STR(create)

    #:set attach_val = GEN_ATTACH_STR(attach)

    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)

    #:set clause_val = copyin_val.strip('\n') + create_val.strip('\n') + attach_val.strip('\n')
    #:set acc_directive = '!$acc enter data ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def EXIT_DATA(copyout=None, delete=None, detach=None, extraAccArgs=None)
    #:set copyout_val = GEN_COPYOUT_STR(copyout)

    #:set delete_val = GEN_DELETE_STR(delete)

    #:set detach_val = GEN_DETACH_STR(detach)

    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)

    #:set clause_val = copyout_val.strip('\n') + delete_val.strip('\n') + detach_val.strip('\n')
    #:set acc_directive = '!$acc exit data ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def CACHE(cache, extraAccArgs=None)
    #:set cache_val = GEN_PARENTHESES_CLAUSE('cache', cache)

    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)

    #:set clause_val = cache_val.strip('\n')
    #:set acc_directive = '!$acc ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def ATOMIC(atomic='update', extraAccArgs=None)
    #:assert isinstance(atomic, str)
    #:assert (atomic == 'read' or atomic == 'write' or atomic == 'update' or atomic == 'capture')

    #:set atomic_val = atomic

    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)

    #:set clause_val = atomic_val.strip('\n')
    #:set acc_directive = '!$acc atomic ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def UPDATE(host=None, device=None, extraAccArgs=None)
    #:set host_val = GEN_HOST_STR(host)

    #:set device_val = GEN_DEVICE_STR(device)

    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)

    #:set clause_val = host_val.strip('\n') + device_val.strip('\n')
    #:set acc_directive = '!$acc update ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def WAIT(host=None, device=None, extraAccArgs=None)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)

    #:set clause_val = host_val.strip('\n') + device_val.strip('\n')
    #:set acc_directive = '!$acc update ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef
