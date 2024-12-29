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
