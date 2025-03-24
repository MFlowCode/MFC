#:def PRIM_GHOST_CELL_EXTRAPOLATION_BC(DEST, SRC)
    do i = 1, sys_size
        do j = 1, buff_size
            q_prim_vf(i)%sf(${DEST}$) = q_prim_vf(i)%sf(${SRC}$)
        end do
    end do
#:enddef

#:def PRIM_SYMMETRY_BC(DIR,DEST,SRC)
    do j = 1, buff_size
        do i = 1, momxb+${DIR}$-2
            q_prim_vf(i)%sf(${DEST}$) = q_prim_vf(i)%sf(${SRC}$)
        end do

        q_prim_vf(i)%sf(${DEST}$) = -q_prim_vf(i)%sf(${SRC}$)

        do i = momxb+${DIR}$, sys_size
            q_prim_vf(i)%sf(${DEST}$) = q_prim_vf(i)%sf(${SRC}$)
        end do

        if (hyperelasticity) then
           q_prim_vf(xibeg + ${DIR}$ - 1)%sf(${DEST}$) = &
                -q_prim_vf(xibeg + ${DIR}$ - 1)%sf(l, j - 1, k)
        end if
    end do
#:enddef

#:def PRIM_PERIODIC_BC(DEST, SRC)
    do i = 1, sys_size
        do j = 1, buff_size
            q_prim_vf(i)%sf(${DEST}$) = q_prim_vf(i)%sf(${SRC}$)
        end do
    end do
#:enddef

#:def PRIM_SLIP_WALL_BC(DIR,LOC)
    #:if DIR == "x"
        #:if LOC == "L"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb) then
                        !q_prim_vf(i)%sf(-j, k, l) = &
                            !-q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb1
                        q_prim_vf(i)%sf(-j, k, l) = 0
                    else
                        q_prim_vf(i)%sf(-j, k, l) = &
                            q_prim_vf(i)%sf(0, k, l)
                    end if
                end do
            end do
        #:elif LOC == "R"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb) then
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve1
                    else
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            q_prim_vf(i)%sf(m, k, l)
                    end if
                end do
            end do
        #:endif
    #:elif DIR == "y"
        #:if LOC == "L"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb + 1) then
                        q_prim_vf(i)%sf(k, -j, l) = &
                            -q_prim_vf(i)%sf(k, j - 1, l) + 2._wp*bc_y%vb2
                    else
                        q_prim_vf(i)%sf(k, -j, l) = &
                            q_prim_vf(i)%sf(k, 0, l)
                    end if
                end do
            end do
        #:elif LOC == "R"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb + 1) then
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            -q_prim_vf(i)%sf(k, n - (j - 1), l) + 2._wp*bc_y%ve2
                    else
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            q_prim_vf(i)%sf(k, n, l)
                    end if
                end do
            end do
        #:endif
    #:elif DIR == "z"
        #:if LOC == "L"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxe) then
                        q_prim_vf(i)%sf(k, l, -j) = &
                            -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb3
                    else
                        q_prim_vf(i)%sf(k, l, -j) = &
                            q_prim_vf(i)%sf(k, l, 0)
                    end if
                end do
            end do
        #:elif LOC == "R"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxe) then
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve3
                    else
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            q_prim_vf(i)%sf(k, l, p)
                    end if
                end do
            end do
        #:endif
    #:endif
#:enddef

#:def PRIM_NO_SLIP_WALL_BC(DIR,LOC)
    #:if DIR == "x"
        #:if LOC == "L"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb) then
                        !q_prim_vf(i)%sf(-j, k, l) = &
                            !-q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb1
                        q_prim_vf(i)%sf(-j, k, l) = 0._wp
                    elseif (i == momxb + 1 .and. num_dims > 1) then
                        !q_prim_vf(i)%sf(-j, k, l) = &
                            !-q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb2
                        q_prim_vf(i)%sf(-j, k, l) = 0._wp
                    elseif (i == momxb + 2 .and. num_dims > 2) then
                        !q_prim_vf(i)%sf(-j, k, l) = &
                            !-q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb3
                        q_prim_vf(i)%sf(-j, k, l) = 0._wp
                    else
                        q_prim_vf(i)%sf(-j, k, l) = &
                            q_prim_vf(i)%sf(0, k, l)
                    end if
                end do
            end do
        #:elif LOC == "R"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb) then
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve1
                    elseif (i == momxb + 1 .and. num_dims > 1) then
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve2
                    elseif (i == momxb + 2 .and. num_dims > 2) then
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve3
                    else
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            q_prim_vf(i)%sf(m, k, l)
                    end if
                end do
            end do
        #:endif
    #:elif DIR == "y"
        #:if LOC == "L"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb) then
                        q_prim_vf(i)%sf(k, -j, l) = &
                            -q_prim_vf(i)%sf(k, j - 1, l) + 2._wp*bc_y%vb1
                    elseif (i == momxb + 1 .and. num_dims > 1) then
                        q_prim_vf(i)%sf(k, -j, l) = &
                            -q_prim_vf(i)%sf(k, j - 1, l) + 2._wp*bc_y%vb2
                    elseif (i == momxb + 2 .and. num_dims > 2) then
                        q_prim_vf(i)%sf(k, -j, l) = &
                            -q_prim_vf(i)%sf(k, j - 1, l) + 2._wp*bc_y%vb3
                    else
                        q_prim_vf(i)%sf(k, -j, l) = &
                            q_prim_vf(i)%sf(k, 0, l)
                    end if
                end do
            end do
        #:elif LOC == "R"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb) then
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            -q_prim_vf(i)%sf(k, n - (j - 1), l) + 2._wp*bc_y%ve1
                    elseif (i == momxb + 1 .and. num_dims > 1) then
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            -q_prim_vf(i)%sf(k, n - (j - 1), l) + 2._wp*bc_y%ve2
                    elseif (i == momxb + 2 .and. num_dims > 2) then
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            -q_prim_vf(i)%sf(k, n - (j - 1), l) + 2._wp*bc_y%ve3
                    else
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            q_prim_vf(i)%sf(k, n, l)
                    end if
                end do
            end do
        #:endif
    #:elif DIR == "z"
        #:if LOC == "L"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb) then
                        q_prim_vf(i)%sf(k, l, -j) = &
                            -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb1
                    elseif (i == momxb + 1 .and. num_dims > 1) then
                        q_prim_vf(i)%sf(k, l, -j) = &
                            -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb2
                    elseif (i == momxb + 2 .and. num_dims > 2) then
                        q_prim_vf(i)%sf(k, l, -j) = &
                            -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb3
                    else
                        q_prim_vf(i)%sf(k, l, -j) = &
                            q_prim_vf(i)%sf(k, l, 0)
                    end if
                end do
            end do
        #:elif LOC == "R"
            do i = 1, sys_size
                do j = 1, buff_size
                    if (i == momxb) then
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve1
                    elseif (i == momxb + 1 .and. num_dims > 1) then
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve2
                    elseif (i == momxb + 2 .and. num_dims > 2) then
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve3
                    else
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            q_prim_vf(i)%sf(k, l, p)
                    end if
                end do
            end do
        #:endif
    #:endif
#:enddef

#:def PRIM_DIRICHLET_BC(DIR, LOC, DEST, SRC)
    do i = 1, sys_size
        do j = 1, buff_size
            q_prim_vf(i)%sf(${DEST}$) = bc_buffers(${DIR}$,${LOC}$)%sf(${SRC}$)
        end do
    end do
#:enddef

#:def QBMM_BC(DEST, SRC)
    do i = 1, nb
        do q = 1, nnode
            do j = 1, buff_size
                pb(${DEST}$) = pb(${SRC}$)
                mv(${DEST}$) = mv(${SRC}$)
            end do
        end do
    end do
#:enddef

