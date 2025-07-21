#:include 'shared_parallel_macros.fpp'
#:include 'omp_macros.fpp'
#:include 'acc_macros.fpp'

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
    & no_create=None, present=None, deviceptr=None, attach=None, extraAccArgs=None, extraOmpArgs=None)

    #:set acc_code = ACC_PARALLEL(code, private, default, firstprivate, reduction, reductionOp, copy, copyin, copyinReadOnly, copyout, create, no_create, present, deviceptr, attach, extraAccArgs)
    #:set omp_code = OMP_PARALLEL(code, private, default, firstprivate, reduction, reductionOp, copy, copyin, copyinReadOnly, copyout, create, no_create, present, deviceptr, attach, extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#else
    $:code
#endif

#:enddef

#:def GPU_PARALLEL_LOOP(code, collapse=None, private=None, parallelism='[gang, vector]', &
    & default='present', firstprivate=None, reduction=None, reductionOp=None, &
    & copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, &
    & no_create=None, present=None, deviceptr=None, attach=None, extraAccArgs=None, extraOmpArgs=None)

    #:set acc_code = ACC_PARALLEL_LOOP(code, collapse, private, parallelism, default, firstprivate, reduction, reductionOp, copy, copyin, copyinReadOnly, copyout, create, no_create, present, deviceptr, attach, extraAccArgs)
    #:set omp_code = OMP_PARALLEL_LOOP(code, collapse, private, parallelism, default, firstprivate, reduction, reductionOp, copy, copyin, copyinReadOnly, copyout, create, no_create, present, deviceptr, attach, extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#endif
#:enddef

#:def GPU_ROUTINE(function_name=None, parallelism=None, nohost=False, cray_inline=False, extraAccArgs=None, extraOmpArgs=None)
    #:assert isinstance(cray_inline, bool)
    #:set acc_directive = ACC_ROUTINE(function_name=function_name, parallelism=parallelism, nohost=nohost, extraAccArgs=extraAccArgs)
    #:set omp_directive = OMP_ROUTINE(function_name=function_name, nohost=nohost, extraOmpArgs=extraOmpArgs)

    #:if cray_inline == True
        #:if not isinstance(function_name, str)
            #:stop "When inlining for Cray Compiler, function name must be given and given as a string"
        #:endif
        #:set cray_directive = ('!DIR$ INLINEALWAYS ' + function_name).strip('\n')
#ifdef _CRAYFTN
        $:cray_directive
#elif MFC_OpenACC
        $:acc_directive
#elif MFC_OpenMP
        $:omp_directive
#endif
    #:else
#if MFC_OpenACC
        $:acc_directive
#elif MFC_OpenMP
        $:omp_directive
#endif
    #:endif
#:enddef

#:def GPU_DECLARE(copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, present=None, deviceptr=None, link=None, extraAccArgs=None, extraOmpArgs=None)
    #:set acc_code = ACC_DECLARE(copy=copy, copyin=copyin, copyinReadOnly=copyinReadOnly, copyout=copyout, create=create, present=present, deviceptr=deviceptr, link=link, extraAccArgs=None)
    #:assert copyout is None
    #:assert present is None
    #:assert deviceptr is None
    #:assert copy is None
    #:set omp_code = OMP_DECLARE(copyin=copyin, copyinReadOnly=copyinReadOnly, create=create, link=link, extraOmpArgs=extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#endif
#:enddef

#:def GPU_LOOP(collapse=None, parallelism=None, data_dependency=None, reduction=None, reductionOp=None, private=None, extraAccArgs=None, extraOmpArgs=None)
    #:set acc_code = ACC_LOOP(collapse=collapse, parallelism=parallelism, data_dependency=data_dependency, reduction=reduction, reductionOp=reductionOp, private=private, extraAccArgs=extraAccArgs)
    #:set omp_code = OMP_LOOP(collapse=collapse, parallelism=parallelism, data_dependency=data_dependency, reduction=reduction, reductionOp=reductionOp, private=private, extraOmpArgs=extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#endif
#:enddef

#:def GPU_DATA(code, copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, no_create=None, present=None, deviceptr=None, attach=None, default=None, extraAccArgs=None, extraOmpArgs=None)
    #:set acc_code = ACC_DATA(code=code, copy=copy, copyin=copyin, copyinReadOnly=copyinReadOnly, copyout=copyout, create=create, no_create=no_create, present=present, deviceptr=deviceptr, attach=attach, default=default, extraAccArgs=extraAccArgs)
    #:set omp_code = OMP_DATA(code=code, copy=copy, copyin=copyin, copyinReadOnly=copyinReadOnly, copyout=copyout, create=create, no_create=no_create, present=present, deviceptr=deviceptr, attach=attach, default=default, extraOmpArgs=extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#else
    $:code
#endif
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

#:def GPU_ENTER_DATA(copyin=None, copyinReadOnly=None, create=None, attach=None, extraAccArgs=None, extraOmpArgs=None)
    #:set acc_code = ACC_ENTER_DATA(copyin=copyin, copyinReadOnly=copyinReadOnly, create=create, attach=attach, extraAccArgs=extraAccArgs)
    #:set omp_code = OMP_ENTER_DATA(copyin=copyin, copyinReadOnly=copyinReadOnly, create=create, attach=attach, extraOmpArgs=extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#endif
#:enddef

#:def GPU_EXIT_DATA(copyout=None, delete=None, detach=None, extraAccArgs=None, extraOmpArgs=None)
    #:set acc_code = ACC_EXIT_DATA(copyout=copyout, delete=delete, detach=detach, extraAccArgs=extraAccArgs)
    #:set omp_code = OMP_EXIT_DATA(copyout=copyout, delete=delete, detach=detach, extraOmpArgs=extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#endif
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

#:def GPU_UPDATE(host=None, device=None, extraAccArgs=None, extraOmpArgs=None)
    #:set acc_code = ACC_UPDATE(host=host, device=device, extraAccArgs=extraAccArgs)
    #:set omp_code = OMP_UPDATE(host=host, device=device, extraOmpArgs=extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#endif
#:enddef

#:def GPU_WAIT(extraAccArgs=None)
    #:set extraAccArgs_val = GEN_EXTRA_ARGS_STR(extraAccArgs)
    #:set clause_val = ''
    #:set acc_directive = '!$acc wait ' + clause_val + extraAccArgs_val.strip('\n')
    $:acc_directive
#:enddef

#:def USE_GPU_MODULE()

#if defined(MFC_OpenACC)
    use openacc
#elif defined(MFC_OpenMP)
    use omp_lib
#endif

#:enddef
! New line at end of file is required for FYPP
