#:def LOG(expr)
#ifdef MFC_DEBUG
    block
        use iso_fortran_env, only: output_unit

        print *, '${_FILE_.split('/')[-1]}$:${_LINE_}$: ', ${expr}$
        call flush(output_unit)
    end block
#endif
#:enddef

#:def ALLOCATE(*args)
    @:LOG({'@:ALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    allocate(${', '.join(args)}$)
    !$acc enter data create(${', '.join(args)}$)
#:enddef ALLOCATE

#:def DEALLOCATE(*args)
    @:LOG({'@:DEALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    deallocate(${', '.join(args)}$)
    !$acc exit data delete(${', '.join(args)}$)
#:enddef DEALLOCATE

#:def ALLOCATE_GLOBAL(*args)
    @:LOG({'@:ALLOCATE_GLOBAL(${re.sub(' +', ' ', ', '.join(args))}$)'})
#ifdef _CRAYFTN
    allocate(${', '.join(('p_' + arg.strip() for arg in args))}$)
    #:for arg in args
        ${re.sub('\(.*\)','',arg)}$ => ${ 'p_' + re.sub('\(.*\)','',arg.strip()) }$
    #:endfor
    !$acc enter data create(${', '.join(('p_' + re.sub('\(.*\)','',arg.strip()) for arg in args))}$) &
    !$acc& attach(${', '.join(map(lambda x: re.sub('\(.*\)','',x), args))}$)
#else
    allocate(${', '.join(args)}$)
    !$acc enter data create(${', '.join(args)}$)
#endif

#:enddef ALLOCATE_GLOBAL

#:def DEALLOCATE_GLOBAL(*args)
    @:LOG({'@:DEALLOCATE_GLOBAL(${re.sub(' +', ' ', ', '.join(args))}$)'})
#ifdef _CRAYFTN
    !$acc exit data delete(${', '.join(('p_' + arg.strip() for arg in args))}$) &
    !$acc& detach(${', '.join(args)}$)
    #:for arg in args
        nullify( ${arg}$ )
    #:endfor
    deallocate(${', '.join(('p_' + arg.strip() for arg in args))}$)
#else
    deallocate(${', '.join(args)}$)
    !$acc exit data delete(${', '.join(args)}$)
#endif

#:enddef DEALLOCATE_GLOBAL

#:def CRAY_DECLARE_GLOBAL(intype, dim, *args)
    ${intype}$, ${dim}$, allocatable, target :: ${', '.join(('p_' + arg.strip() for arg in args))}$
    ${intype}$, ${dim}$, pointer :: ${', '.join(args)}$
#:enddef CRAY_DECLARE_GLOBAL

#:def ACC_SETUP_VFs(*args)
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
#:enddef
