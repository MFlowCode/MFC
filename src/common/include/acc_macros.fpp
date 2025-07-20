#:include 'shared_parallel_macros.fpp'

#:def ACC_PARALLEL(code, private=None, default='present', firstprivate=None, reduction=None, reductionOp=None, &
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
    #:set acc_clause_val = default_val.strip('\n') + private_val.strip('\n') + reduction_val.strip('\n') + &
        & copy_val.strip('\n') + copyin_val.strip('\n') + &
        & copyout_val.strip('\n') + create_val.strip('\n') + &
        & no_create_val.strip('\n') + present_val.strip('\n') + &
        & deviceptr_val.strip('\n') + attach_val.strip('\n')
    #:set acc_directive = '!$acc parallel ' + &
        & acc_clause_val + extraAccArgs_val.strip('\n')
    #:set end_acc_directive = '!$acc end parallel'
    $:acc_directive
    $:code
    $:end_acc_directive
#:enddef

#:def ACC_PARALLEL_LOOP(collapse=None, private=None, parallelism='[gang, vector]', &
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
! New line at end of file is required for FYPP