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

#:def parallel_loop(collapse=None, private=None, parallelism=["gang", "vector"], default="present", firstprivate=None, reduction=None, reductionOp=None, copy=None, copyin=None, copyinReadOnly=False, copyout=None, create=None)
    #:if collapse is not None
        #:assert isinstance(collapse, int)
        #:assert collapse > 0
        #:set collapse_val = 'collapse(' + str(collapse) + ') '
    #:else
        #:set collapse_val = ""
    #:endif

    #:if private is not None
        #:assert isinstance(private, list)
        #:assert len(private) != 0
        #:assert all(type(element) == str for element in private)
        #:set private_val = 'private(' + ', '.join(private) + ') '
    #:else
        #:set private_val = ""
    #:endif

    #:if default is not None
        #:assert isinstance(default, str)
        #:assert (default == "present" or default == "none")
        #:set default_val = 'default(' + default + ') '
    #:else
        #:set default_val = ""
    #:endif

    #:if parallelism is not None
        #:assert isinstance(parallelism, list)
        #:assert len(parallelism) != 0
        #:assert all(type(element) == str for element in parallelism)
        #:assert all((element == "gang" or element == "worker" or element == "vector") for element in parallelism)
        #:set parallelism_val = " ".join(parallelism) + " "
    #:else
        #:set parallelism_val = ""
    #:endif

    #:if firstprivate is not None
        #:assert isinstance(firstprivate, list)
        #:assert len(firstprivate) != 0
        #:assert all(type(element) == str for element in firstprivate)
        #:set firstprivate_val = 'firstprivate(' + ', '.join(firstprivate) + ') '
    #:else
        #:set firstprivate_val = ""
    #:endif

    #:if reduction is not None and reductionOp is not None
        #:assert isinstance(reduction, list)
        #:assert len(reduction) != 0
        #:assert all(type(element) == str for element in reduction)
        #:assert isintance(reductionOp, str)
        #:set reduction_val = 'reduction(' + reductionOp + ':' + ', '.join(reduction) + ') '
    #:elif reduction is not None or reductionOp is not None
        #:stop "Cannot set the reduction list or reduction operation without setting the other"
    #:else
        #:set reduction_val = ""
    #:endif

    #:if copy is not None
        #:assert isinstance(copy, list)
        #:assert len(copy) != 0
        #:assert all(type(element) == str for element in copy)
        #:set copy_val = 'copy(' + ', '.join(copy) + ') '
    #:else
        #:set copy_val = ""
    #:endif
    
    #:if copyin is not None
        #:assert isinstance(copyin, list)
        #:assert len(copyin) != 0
        #:assert all(type(element) == str for element in copyin)
        #:assert isinstance(copyinReadOnly, bool)
        #:if copyinReadOnly == True
            #:set copyin_val = 'copyin(readonly:' + ', '.join(copyin) + ') '
        #:else
            #:set copyin_val = 'copyin(' + ', '.join(copyin) + ') '
        #:endif
    #:else
        #:set copyin_val = ""
    #:endif

    #:if copyout is not None
        #:assert isinstance(copyout, list)
        #:assert len(copyout) != 0
        #:assert all(type(element) == str for element in copyout)
        #:set copyout_val = 'copyout(' + ', '.join(copyout) + ') '
    #:else
        #:set copyout_val = ""
    #:endif

    #:if create is not None
        #:assert isinstance(create, list)
        #:assert len(create) != 0
        #:assert all(type(element) == str for element in create)
        #:set create_val = 'create(' + ', '.join(create) + ') '
    #:else
        #:set create_val = ""
    #:endif

    #:set clause_val = collapse_val + parallelism_val + default_val + private_val + firstprivate_val + reduction_val + copy_val + copyin_val + copyout_val + create_val
    !$acc parallel loop ${clause_val}$

#:enddef
