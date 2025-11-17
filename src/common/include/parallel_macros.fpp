#:include 'shared_parallel_macros.fpp'
#:include 'omp_macros.fpp'
#:include 'acc_macros.fpp'

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

#:def GPU_PARALLEL_LOOP_OLD(code, collapse=None, private=None, parallelism='[gang, vector]', &
    & default='present', firstprivate=None, reduction=None, reductionOp=None, &
    & copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, &
    & no_create=None, present=None, deviceptr=None, attach=None, extraAccArgs=None, extraOmpArgs=None)

    #:set acc_code = ACC_PARALLEL_LOOP_OLD(code, collapse, private, parallelism, default, firstprivate, reduction, reductionOp, copy, copyin, copyinReadOnly, copyout, create, no_create, present, deviceptr, attach, extraAccArgs)
    #:set omp_code = OMP_PARALLEL_LOOP_OLD(code, collapse, private, parallelism, default, firstprivate, reduction, reductionOp, copy, copyin, copyinReadOnly, copyout, create, no_create, present, deviceptr, attach, extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#else
    $:code
#endif
#:enddef

#:def GPU_PARALLEL_LOOP(collapse=None, private=None, parallelism='[gang, vector]', &
    & default='present', firstprivate=None, reduction=None, reductionOp=None, &
    & copy=None, copyin=None, copyinReadOnly=None, copyout=None, create=None, &
    & no_create=None, present=None, deviceptr=None, attach=None, extraAccArgs=None, extraOmpArgs=None)

    #:set acc_directive = ACC_PARALLEL_LOOP(collapse, private, parallelism, default, firstprivate, reduction, reductionOp, copy, copyin, copyinReadOnly, copyout, create, no_create, present, deviceptr, attach, extraAccArgs)
    #:set omp_directive = OMP_PARALLEL_LOOP(collapse, private, parallelism, default, firstprivate, reduction, reductionOp, copy, copyin, copyinReadOnly, copyout, create, no_create, present, deviceptr, attach, extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_directive
#elif defined(MFC_OpenMP)
    $:omp_directive
#endif

#:enddef

#:def END_GPU_PARALLEL_LOOP()

    #:set acc_end_directive = '!$acc end parallel loop'
    #:set omp_end_directive = END_OMP_PARALLEL_LOOP()

#if defined(MFC_OpenACC)
    $:acc_end_directive
#elif defined(MFC_OpenMP)
    $:omp_end_directive
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

#:def GPU_HOST_DATA(code, use_device_addr=None, use_device_ptr=None, extraAccArgs=None, extraOmpArgs=None)

    #:if use_device_addr is not None and use_device_ptr is not None
        #:set use_device_addr_end_index = len(use_device_addr) - 1
        #:set use_device = use_device_addr + use_device_ptr
        $:use_device[use_device_addr_end_index] = ','
        $:use_device[use_device_addr_end_index + 1] = ' '
    #:elif use_device_addr is not None or use_device_ptr is not None
        #:if use_device_addr is not None
            #:set use_device = use_device_addr
        #:else
            #:set use_device = use_device_ptr
        #:endif
    #:else
        #:set use_device = None
    #:endif
    #:set acc_code = ACC_HOST_DATA(code=code, use_device=use_device, extraAccArgs=extraAccArgs)
    #:set omp_code = OMP_HOST_DATA(code=code, use_device_addr=use_device_addr, use_device_ptr=use_device_ptr, extraOmpArgs=extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#else
    $:code
#endif
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

#:def GPU_ATOMIC(atomic, extraAccArgs=None, extraOmpArgs=None)
    #:set acc_code = ACC_ATOMIC(atomic=atomic, extraAccArgs=extraAccArgs)
    #:set omp_code = OMP_ATOMIC(atomic=atomic, extraOmpArgs=extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#endif
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

#:def GPU_WAIT(extraAccArgs=None, extraOmpArgs=None)
    #:set acc_code = ACC_WAIT(extraAccArgs=extraAccArgs)
    #:set omp_code = OMP_WAIT(extraOmpArgs=extraOmpArgs)

#if defined(MFC_OpenACC)
    $:acc_code
#elif defined(MFC_OpenMP)
    $:omp_code
#endif
#:enddef

#:def USE_GPU_MODULE()

#if defined(MFC_OpenACC)
    use openacc
#elif defined(MFC_OpenMP)
    use omp_lib
#endif

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

#:set USING_NVHPC = (MFC_COMPILER == NVIDIA_COMPILER_ID or MFC_COMPILER == PGI_COMPILER_ID)
#:set USING_CCE = (MFC_COMPILER == CCE_COMPILER_ID)
! New line at end of file is required for FYPP
