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

#:def GEN_COPY_STR(copy)
    #:set copy_val = GEN_PARENTHESES_CLAUSE('copy', copy)
    $:copy_val
#:enddef

#:def GEN_COPYIN_STR(copyin, copyinReadOnly)
    #:if copyin is not None
        $:ASSERT_LIST(copyin, str)
        #:assert isinstance(copyinReadOnly, bool)
        #:if copyinReadOnly == True
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
    #:set parallelism_val = parallelism_val.strip('\n')
    $:parallelism_val
#:enddef

#:def PARALLEL_LOOP(collapse=None, private=None, parallelism=['gang', 'vector'], &
    & default='present', firstprivate=None, reduction=None, reductionOp=None, &
    & copy=None, copyin=None, copyinReadOnly=False, copyout=None, create=None, &
    & extraAccArgs=None)
    #:if collapse is not None
        #:assert isinstance(collapse, int)
        #:assert collapse > 1
        #:set collapse_val = 'collapse(' + str(collapse) + ') '
    #:else
        #:set collapse_val = ''
    #:endif

    #:set private_val = GEN_PARENTHESES_CLAUSE('private', private)

    #:if default is not None
        #:assert isinstance(default, str)
        #:assert (default == 'present' or default == 'none')
        #:set default_val = 'default(' + default + ') '
    #:else
        #:set default_val = ''
    #:endif

    #:set parallelism_val = GEN_PARALLELISM_STR(parallelism)

    #:set firstprivate_val = GEN_PARENTHESES_CLAUSE('firstprivate', firstprivate)

    #:if reduction is not None and reductionOp is not None
        #:if isinstance(reduction, list) and isinstance(reductionOp, list)
            $:ASSERT_LIST(reduction, list)
            $:ASSERT_LIST(reductionOp, str)
            #:assert all(len(element) != 0 for element in reduction)
            #:assert all(type(element) == str for sublist in reduction for element in sublist)
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

    #:set copy_val = GEN_COPY_STR(copy)

    #:set copyin_val = GEN_COPYIN_STR(copyin, copyinReadOnly)

    #:set copyout_val = GEN_COPYOUT_STR(copyout)

    #:set create_val = GEN_CREATE_STR(create)

    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)

    #:set clause_val = collapse_val.strip('\n') + parallelism_val.strip('\n') + &
        & default_val.strip('\n') + private_val.strip('\n') + &
        & firstprivate_val.strip('\n') + reduction_val.strip('\n') + &
        & copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n')
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

#:def DECLARE(copy=None, copyin=None, copyinReadOnly=False, copyout=None, create=None, present=None, deviceptr=None, link=None, extraAccArgs=None)
    #:set copy_val = GEN_COPY_STR(copy)

    #:set copyin_val = GEN_COPYIN_STR(copyin, copyinReadOnly)

    #:set copyout_val = GEN_COPYOUT_STR(copyout)

    #:set create_val = GEN_CREATE_STR(create)

    #:set present_val = GEN_PARENTHESES_CLAUSE('present', present)

    #:set deviceptr_val = GEN_PARENTHESES_CLAUSE('deviceptr', deviceptr)

    #:set link_val = GEN_PARENTHESES_CLAUSE('link', link)

    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)

    #:set clause_val = copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & present_val.strip('\n') + deviceptr_val.strip('\n') + &
        & link_val.strip('\n')
    #:set acc_directive = '!$acc declare ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def LOOP(collapse=None, parallelism=["seq"], data_dependency=None, reduction_list=None, reductionOp=None, private=None, extraAccArgs=None)
    #:if collapse is not None
        #:assert isinstance(collapse, int)
        #:assert collapse > 1
        #:set collapse_val = 'collapse(' + str(collapse) + ') '
    #:else
        #:set collapse_val = ''
    #:endif

    #:set parallelism_val = GEN_PARALLELISM_STR(parallelism)


    #:if data_dependency is not None
        #:assert isinstance(data_dependency, str)
        #:assert (data_dependency == 'auto' or data_dependency == 'independent')
        #:set data_dependency_val = data_dependency 
    #:else
        #:set data_dependency_val = ''
    #:endif
    
    #:set private_val = GEN_PARENTHESES_CLAUSE('private', private)

    #:if reduction is not None and reductionOp is not None
        #:if isinstance(reduction, list) and isinstance(reductionOp, list)
            $:ASSERT_LIST(reduction, list)
            $:ASSERT_LIST(reductionOp, str)
            #:assert all(len(element) != 0 for element in reduction)
            #:assert all(type(element) == str for sublist in reduction for element in sublist)
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

    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)

    #:set clause_val = collapse_val.strip('\n') + parallelism_val.strip('\n') + &
        & data_dependency_val.strip('\n') + private_val.strip('\n') + &
        & reduction_val.strip('\n')
    #:set acc_directive = '!$acc serial ' + &
        & clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef