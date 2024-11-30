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

#define t_vec3   real(kind(0d0)), dimension(1:3)
#define t_mat4x4 real(kind(0d0)), dimension(1:4,1:4)

#:def ASSERT(predicate, message = None)
    if (.not. (${predicate}$)) then
        call s_mpi_abort("${_FILE_.split('/')[-1]}$:${_LINE_}$: "// &
                         "Assertion failed: ${predicate}$. " &
                         //${message or '"No error description."'}$)
    end if
#:enddef

#:def MAKE_LOOPS(impl, loops, gpu = None)
    #:if gpu != False
        !$acc parallel loop collapse(${len(loops)}$) gang vector default(present)
    #:endif
    #:for index, lbound, hbound in (loops or [])
        do ${index}$ = ${lbound}$, ${hbound}$
    #:endfor

    $:impl

    #:for i in range(len(loops or []))
        end do
        #:endfor
    #:enddef

    #:def ITERATE_OVER_BUFFER_REGION_SIDED(impl, dir = None, loc = None, thickness = None, gpu = None, outer_loops = None, inner_loops = None, pack_v_size = None)
        #:set thickness   = thickness or 'buff_size'
        #:set outer_loops = outer_loops or []
        #:set inner_loops = inner_loops or []
        #:set pack_v_size = pack_v_size or 'v_size'

        if (${dir}$ <= num_dims) then

            if ((${dir}$) == +1 .and. (${loc}$) == -1) then

                #:block MAKE_LOOPS(loops = outer_loops + [("l", 0, "p"), ("k", 0, "n"), ("j", 1, thickness)] + inner_loops, gpu=gpu)
                    x = -j; y = k; z = l ! Regular
                    sx = j - 1; sy = k; sz = l ! Symmetry
                    px = m - (j - 1); py = k; pz = l ! Periodic
                    ex = 0; ey = k; ez = l ! Extrapolation
                    locx = -1; locy = 0; locz = 0

                    exlhs = ex; eylhs = ey; ezlhs = ez
                    exlhs = f_clamp(exlhs, 0, m); eylhs = f_clamp(eylhs, 0, n); ezlhs = f_clamp(ezlhs, 0, p)

#ifndef MFC_PRE_PROCESS
                    pack_idr = (i - 1) + ${pack_v_size}$*((j - 1) + ${thickness}$*(k + (n + 1)*l))
#endif

                    $:impl
                #:endblock

            else if ((${dir}$) == +1 .and. (${loc}$) == +1) then

                #:block MAKE_LOOPS(loops = outer_loops + [("l", 0, "p"), ("k", 0, "n"), ("j", 1, thickness)] + inner_loops, gpu=gpu)
                    x = m + j; y = k; z = l ! Regular
                    sx = m - (j - 1); sy = k; sz = l ! Symmetry
                    px = j - 1; py = k; pz = l ! Periodic
                    ex = m; ey = k; ez = l ! Extrapolation
                    locx = +1; locy = 0; locz = 0

                    exlhs = 0; eylhs = ey; ezlhs = ez
                    exlhs = f_clamp(exlhs, 0, m); eylhs = f_clamp(eylhs, 0, n); ezlhs = f_clamp(ezlhs, 0, p)

#ifndef MFC_PRE_PROCESS
                    pack_idr = (i - 1) + ${pack_v_size}$*((j - 1) + ${thickness}$*(k + (n + 1)*l))
#endif

                    $:impl
                #:endblock

            else if ((${dir}$) == +2 .and. (${loc}$) == -1) then

                #:block MAKE_LOOPS(loops = outer_loops + [("k", 0, "p"), ("j", 1, thickness), ("l", f"-({thickness})", f"m + ({thickness})")] + inner_loops, gpu=gpu)
                    x = l; y = -j; z = k ! Regular
                    sx = l; sy = j - 1; sz = k ! Symmetry
                    px = l; py = n - (j - 1); pz = k ! Periodic
                    ex = l; ey = 0; ez = k ! Extrapolation
                    locx = 0; locy = -1; locz = 0

                    exlhs = ex; eylhs = ey; ezlhs = ez
                    exlhs = f_clamp(exlhs, 0, m); eylhs = f_clamp(eylhs, 0, n); ezlhs = f_clamp(ezlhs, 0, p)

#ifndef MFC_PRE_PROCESS
                    pack_idr = (i - 1) + ${pack_v_size}$*((l + ${thickness}$) + (m + 2*${thickness}$+1)*(j - 1 + ${thickness}$*k))
#endif

                    $:impl
                #:endblock

            else if ((${dir}$) == +2 .and. (${loc}$) == +1) then

                #:block MAKE_LOOPS(loops = outer_loops + [("k", 0, "p"), ("j", 1, thickness), ("l", f"-({thickness})", f"m + ({thickness})")] + inner_loops, gpu=gpu)
                    x = l; y = n + j; z = k ! Regular
                    sx = l; sy = n - (j - 1); sz = k ! Symmetry
                    px = l; py = j - 1; pz = k ! Periodic
                    ex = l; ey = n; ez = k ! Extrapolation
                    locx = 0; locy = +1; locz = 0

                    exlhs = ex; eylhs = 0; ezlhs = ez
                    exlhs = f_clamp(exlhs, 0, m); eylhs = f_clamp(eylhs, 0, n); ezlhs = f_clamp(ezlhs, 0, p)

#ifndef MFC_PRE_PROCESS
                    pack_idr = (i - 1) + ${pack_v_size}$*((l + ${thickness}$) + (m + 2*${thickness}$+1)*(j - 1 + ${thickness}$*k))
#endif

                    $:impl
                #:endblock

            else if ((${dir}$) == +3 .and. (${loc}$) == -1) then

                #:block MAKE_LOOPS(loops = outer_loops + [("j", 1, thickness), ("l", f"-({thickness})", f"n + ({thickness})"), ("k", f"-({thickness})", f"m + ({thickness})")] + inner_loops, gpu=gpu)
                    x = k; y = l; z = -j          ! Regular
                    sx = k; sy = l; sz = j - 1       ! Symmetry
                    px = k; py = l; pz = p - (j - 1) ! Periodic
                    ex = k; ey = l; ez = 0           ! Extrapolation
                    locx = 0; locy = 0; locz = -1

                    exlhs = ex; eylhs = ey; ezlhs = ez
                    exlhs = f_clamp(exlhs, 0, m); eylhs = f_clamp(eylhs, 0, n); ezlhs = f_clamp(ezlhs, 0, p)

#ifndef MFC_PRE_PROCESS
                    pack_idr = (i - 1) + ${pack_v_size}$*((k + ${thickness}$) + (m + 2*${thickness}$+1)*((l + ${thickness}$) + (n + 2*${thickness}$+1)*(j - 1)))
#endif

                    $:impl
                #:endblock

            else if ((${dir}$) == +3 .and. (${loc}$) == +1) then

                #:block MAKE_LOOPS(loops = outer_loops + [("j", 1, thickness), ("l", f"-({thickness})", f"n + ({thickness})"), ("k", f"-({thickness})", f"m + ({thickness})")] + inner_loops, gpu=gpu)
                    x = k; y = l; z = p + j       ! Regular
                    sx = k; sy = l; sz = p - (j - 1) ! Symmetry
                    px = k; py = l; pz = j - 1       ! Periodic
                    ex = k; ey = l; ez = p           ! Extrapolation
                    locx = 0; locy = 0; locz = +1

                    exlhs = ex; eylhs = ey; ezlhs = 0
                    exlhs = f_clamp(exlhs, 0, m); eylhs = f_clamp(eylhs, 0, n); ezlhs = f_clamp(ezlhs, 0, p)

#ifndef MFC_PRE_PROCESS
                    pack_idr = (i - 1) + ${pack_v_size}$*((k + ${thickness}$) + (m + 2*${thickness}$+1)*((l + ${thickness}$) + (n + 2*${thickness}$+1)*(j - 1)))
#endif

                    $:impl
                #:endblock

            else

                stop "Invalid boundary condition direction or location."

            end if

        end if
    #:enddef ITERATE_OVER_BUFFER_REGION_SIDED

    #:def BOUNDARY_CONDITION_INTEGER_DECLARATIONS()
        integer :: i, j, k, l, q

        integer :: x, y, z
        integer :: sx, sy, sz
        integer :: px, py, pz
        integer :: ex, ey, ez
        integer :: locx, locy, locz
        integer :: exlhs, eylhs, ezlhs
        integer :: pack_idr
    #:enddef
