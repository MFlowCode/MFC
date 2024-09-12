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
#ifndef CRAY_ACC_WAR
!$acc enter data create(${', '.join(args)}$)
#endif
#:enddef ALLOCATE

#:def DEALLOCATE(*args)
    @:LOG({'@:DEALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    deallocate (${', '.join(args)}$)
#ifndef CRAY_ACC_WAR
!$acc exit data delete(${', '.join(args)}$)
#endif
#:enddef DEALLOCATE

#:def ALLOCATE_GLOBAL(*args)
    @:LOG({'@:ALLOCATE_GLOBAL(${re.sub(' +', ' ', ', '.join(args))}$)'})
#ifdef CRAY_ACC_WAR
    allocate (${', '.join(('p_' + arg.strip() for arg in args))}$)
    #:for arg in args
        ${re.sub('\\(.*\\)','',arg)}$ => ${ 'p_' + re.sub('\\(.*\\)','',arg.strip()) }$
    #:endfor
    !$acc enter data create(${', '.join(('p_' + re.sub('\\(.*\\)','',arg.strip()) for arg in args))}$) &
    !$acc& attach(${', '.join(map(lambda x: re.sub('\\(.*\\)','',x), args))}$)
#else
    allocate (${', '.join(args)}$)
    !$acc enter data create(${', '.join(args)}$)
#endif

#:enddef ALLOCATE_GLOBAL

#:def DEALLOCATE_GLOBAL(*args)
    @:LOG({'@:DEALLOCATE_GLOBAL(${re.sub(' +', ' ', ', '.join(args))}$)'})
#ifdef CRAY_ACC_WAR
    !$acc exit data delete(${', '.join(('p_' + arg.strip() for arg in args))}$) &
    !$acc& detach(${', '.join(args)}$)
    #:for arg in args
        nullify (${arg}$)
    #:endfor
    deallocate (${', '.join(('p_' + arg.strip() for arg in args))}$)
#else
    deallocate (${', '.join(args)}$)
    !$acc exit data delete(${', '.join(args)}$)
#endif

#:enddef DEALLOCATE_GLOBAL

#:def CRAY_DECLARE_GLOBAL(intype, dim, *args)
#ifdef CRAY_ACC_WAR
    ${intype}$, ${dim}$, allocatable, target :: ${', '.join(('p_' + arg.strip() for arg in args))}$
    ${intype}$, ${dim}$, pointer :: ${', '.join(args)}$
#else
    ${intype}$, ${dim}$, allocatable :: ${', '.join(args)}$
#endif
#:enddef CRAY_DECLARE_GLOBAL

#:def CRAY_DECLARE_GLOBAL_SCALAR(intype, *args)
#ifdef CRAY_ACC_WAR
    ${intype}$, target :: ${', '.join(('p_' + arg.strip() for arg in args))}$
    ${intype}$, pointer :: ${', '.join(args)}$
#else
    ${intype}$::${', '.join(args)}$
#endif
#:enddef CRAY_DECLARE_GLOBAL_SCALAR

#:def ACC_SETUP_VFs(*args)
#ifdef CRAY_ACC_WAR
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
#ifdef CRAY_ACC_WAR
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
#ifdef CRAY_ACC_WAR
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

#:def PROHIBIT(*args)
    #:set condition = args[0]
    #:if len(args) == 1
        #:set message = '""'
    #:else
        #:set message = args[1]
    #:endif
    if (${condition}$) then
        call s_prohibit_abort("${condition}$", ${message}$)
    end if
#:enddef

#define t_vec3   real(kind(0d0)), dimension(1:3)
#define t_mat4x4 real(kind(0d0)), dimension(1:4,1:4)
