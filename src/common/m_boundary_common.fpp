!>
!! @file m_boundary_conditions_common.fpp
!! @brief Contains module m_boundary_conditions_common

!> @brief The purpose of the module is to apply noncharacteristic and processor
!! boundary condiitons

#:include 'macros.fpp'
#:include 'inline_boundary_conditions.fpp'

module m_boundary_common

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy

    use m_constants

    implicit none

    type(scalar_field), dimension(:, :), allocatable :: bc_buffers
    !$acc declare create(bc_buffers)

    real(wp) :: bcxb, bcxe, bcyb, bcye, bczb, bcze

#ifdef MFC_MPI
    integer, dimension(1:3, -1:1) :: MPI_BC_TYPE_TYPE, MPI_BC_BUFFER_TYPE
#endif

#ifdef MFC_SIMULATION
    private; public :: s_initialize_boundary_common_module, &
 s_populate_variables_buffers, &
 s_create_mpi_types, &
 s_populate_capillary_buffers, &
 s_finalize_boundary_common_module
#else
    private; public :: s_initialize_boundary_common_module, &
 s_populate_variables_buffers, &
 s_create_mpi_types, &
 s_finalize_boundary_common_module
#endif

    public :: bc_buffers, bcxb, bcxe, bcyb, bcye, bczb, bcze

#ifdef MFC_MPI
    public :: MPI_BC_TYPE_TYPE, MPI_BC_BUFFER_TYPE
#endif

contains

    subroutine s_initialize_boundary_common_module()

        bcxb = bc_x%beg; bcxe = bc_x%end; bcyb = bc_y%beg; bcye = bc_y%end; bczb = bc_z%beg; bcze = bc_z%end

        @:ALLOCATE(bc_buffers(1:num_dims, -1:1))

#ifndef MFC_POST_PROCESS

#ifdef MFC_PRE_PROCESS
        if (save_bc) then
#elif MFC_SIMULATION
        if (read_bc) then
#endif
            @:ALLOCATE(bc_buffers(1, -1)%sf(1:sys_size, 0:n, 0:p))
            @:ALLOCATE(bc_buffers(1, 1)%sf(1:sys_size, 0:n, 0:p))
            @:ACC_SETUP_SFs(bc_buffers(1,-1), bc_buffers(1,1))
            if (n > 0) then
                @:ALLOCATE(bc_buffers(2,-1)%sf(-buff_size:m+buff_size,1:sys_size,0:p))
                @:ALLOCATE(bc_buffers(2,1)%sf(-buff_size:m+buff_size,1:sys_size,0:p))
                @:ACC_SETUP_SFs(bc_buffers(2,-1), bc_buffers(2,1))
                if (p > 0) then
                    @:ALLOCATE(bc_buffers(3,-1)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,1:sys_size))
                    @:ALLOCATE(bc_buffers(3,1)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,1:sys_size))
                    @:ACC_SETUP_SFs(bc_buffers(3,-1), bc_buffers(3,1))
                end if
            end if
        end if
#endif

    end subroutine s_initialize_boundary_common_module

    subroutine s_populate_variables_buffers(q_prim_vf, pb, mv, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv

        integer :: i, j, k, l, q

        !< x-direction
        if (bcxb >= 0) then
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 1, -1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    if (bc_type(1, -1)%sf(0, k, l) >= -13 .and. bc_type(1, -1)%sf(0, k, l) <= -3) then
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("-j,k,l","0,k,l")}$
                    elseif (bc_type(1, -1)%sf(0, k, l) == -2) then
                        ${PRIM_SYMMETRY_BC(1,"-j,k,l","j-1,k,l")}$
                    elseif (bc_type(1, -1)%sf(0, k, l) == -1) then
                        ${PRIM_PERIODIC_BC("-j,k,l","m-(j-1),k,l")}$
                    elseif (bc_type(1, -1)%sf(0, k, l) == -15) then
                        ${PRIM_SLIP_WALL_BC("x","L")}$
                    elseif (bc_type(1, -1)%sf(0, k, l) == -16) then
                        ${PRIM_NO_SLIP_WALL_BC("x","L")}$
                    elseif (bc_type(1, -1)%sf(0, k, l) == -17) then
#ifdef MFC_PRE_PROCESS
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("-j,k,l","0,k,l")}$
#else
                        ${PRIM_DIRICHLET_BC(1,-1,"-j,k,l","i,k,l")}$
#endif
                    end if
                end do
            end do
        end if

        if (bcxe >= 0) then
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 1, 1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    if (bc_type(1, 1)%sf(0, k, l) >= -13 .and. bc_type(1, 1)%sf(0, k, l) <= -3) then
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("m+j,k,l","m,k,l")}$
                    elseif (bc_type(1, 1)%sf(0, k, l) == -2) then
                        ${PRIM_SYMMETRY_BC(1,"m+j,k,l","m - (j-1),k,l")}$
                    elseif (bc_type(1, 1)%sf(0, k, l) == -1) then
                        ${PRIM_PERIODIC_BC("m+j,k,l","j-1,k,l")}$
                    elseif (bc_type(1, 1)%sf(0, k, l) == -15) then
                        ${PRIM_SLIP_WALL_BC("x","R")}$
                    elseif (bc_type(1, 1)%sf(0, k, l) == -16) then
                        ${PRIM_NO_SLIP_WALL_BC("x","R")}$
                    elseif (bc_type(1, 1)%sf(0, k, l) == -17) then
#ifdef MFC_PRE_PROCESS
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("m+j,k,l","m,k,l")}$
#else
                        ${PRIM_DIRICHLET_BC(1,1,"m+j,k,l","i,k,l")}$
#endif
                    end if
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            if (bcxb < 0) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = 0, p
                    do k = 0, n
                        if (bc_type(1, -1)%sf(0, k, l) >= -13 .and. bc_type(1, -1)%sf(0, k, l) <= -3) then
                            ${QBMM_BC("-j,k,l,q,i","0,k,l,q,i")}$
                        elseif (bc_type(1, -1)%sf(0, k, l) == -2) then
                            ${QBMM_BC("-j,k,l,q,i","j-1,k,l,q,i")}$
                        elseif (bc_type(1, -1)%sf(0, k, l) == -1) then
                            ${QBMM_BC("-j,k,l,q,i","m - (j-1),k,l,q,i")}$
                        end if
                    end do
                end do
            end if

            if (bcxe < 0) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = 0, p
                    do k = 0, n
                        if (bc_type(1, 1)%sf(0, k, l) >= -13 .and. bc_type(1, 1)%sf(0, k, l) <= -3) then
                            ${QBMM_BC("m+j,k,l,q,i","m,k,l,q,i")}$
                        elseif (bc_type(1, 1)%sf(0, k, l) == -2) then
                            ${QBMM_BC("m+j,k,l,q,i","m - (j-1),k,l,q,i")}$
                        elseif (bc_type(1, 1)%sf(0, k, l) == -1) then
                            ${QBMM_BC("m+j,k,l,q,i","j-1,k,l,q,i")}$
                        end if
                    end do
                end do
            end if
        end if

        if (n == 0) return

        !< y-direction
        if (bcyb >= 0) then
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 2, -1)
        elseif (bcyb == -14) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
                        if (z_cc(k) < pi) then
                            !$acc loop seq
                            do i = 1, momxb
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k + ((p + 1)/2))
                            end do

                            q_prim_vf(momxb + 1)%sf(l, -j, k) = &
                                -q_prim_vf(momxb + 1)%sf(l, j - 1, k + ((p + 1)/2))

                            q_prim_vf(momxe)%sf(l, -j, k) = &
                                -q_prim_vf(momxe)%sf(l, j - 1, k + ((p + 1)/2))

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k + ((p + 1)/2))
                            end do
                        else
                            !$acc loop seq
                            do i = 1, momxb
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k - ((p + 1)/2))
                            end do

                            q_prim_vf(momxb + 1)%sf(l, -j, k) = &
                                -q_prim_vf(momxb + 1)%sf(l, j - 1, k - ((p + 1)/2))

                            q_prim_vf(momxe)%sf(l, -j, k) = &
                                -q_prim_vf(momxe)%sf(l, j - 1, k - ((p + 1)/2))

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k - ((p + 1)/2))
                            end do
                        end if
                    end do
                end do
            end do
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = -buff_size, m + buff_size
                    if (bc_type(2, -1)%sf(k, 0, l) >= -13 .and. bc_type(2, -1)%sf(k, 0, l) <= -3) then
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("k,-j,l","k,0,l")}$
                    elseif (bc_type(2, -1)%sf(k, 0, l) == -2) then
                        ${PRIM_SYMMETRY_BC(2,"k,-j,l","k,j-1,l")}$
                    elseif (bc_type(2, -1)%sf(k, 0, l) == -1) then
                        ${PRIM_PERIODIC_BC("k,-j,l","k,n-(j-1),l")}$
                    elseif (bc_type(2, -1)%sf(k, 0, l) == -15) then
                        ${PRIM_SLIP_WALL_BC("y","L")}$
                    elseif (bc_type(2, -1)%sf(k, 0, l) == -16) then
                        ${PRIM_NO_SLIP_WALL_BC("y","L")}$
                    elseif (bc_type(2, -1)%sf(k, 0, l) == -17) then
#ifdef MFC_PRE_PROCESS
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("k,-j,l","k,0,l")}$
#else
                        ${PRIM_DIRICHLET_BC(2,-1,"k,-j,l","k,i,l")}$
#endif
                    end if
                end do
            end do
        end if

        if (bcye >= 0) then
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 2, 1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = -buff_size, m + buff_size
                    if (bc_type(2, 1)%sf(k, 0, l) >= -13 .and. bc_type(2, 1)%sf(k, 0, l) <= -3) then
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("k,n+j,l","k,n,l")}$
                    elseif (bc_type(2, 1)%sf(k, 0, l) == -2) then
                        ${PRIM_SYMMETRY_BC(2,"k,n+j,l","k,n - (j-1),l")}$
                    elseif (bc_type(2, 1)%sf(k, 0, l) == -1) then
                        ${PRIM_PERIODIC_BC("k,n+j,l","k,j-1,l")}$
                    elseif (bc_type(2, 1)%sf(k, 0, l) == -15) then
                        ${PRIM_SLIP_WALL_BC("y","R")}$
                    elseif (bc_type(2, 1)%sf(k, 0, l) == -16) then
                        ${PRIM_NO_SLIP_WALL_BC("y","R")}$
                    elseif (bc_type(2, 1)%sf(k, 0, l) == -17) then
#ifdef FMC_PRE_PROCESS
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("k,n+j,l","k,n,l")}$
#else
                        ${PRIM_DIRICHLET_BC(2,1,"k,n+j,l","k,i,l")}$
#endif
                    end if
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            if (bcyb < 0) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = 0, p
                    do k = -buff_size, m + buff_size
                        if (bc_type(2, -1)%sf(k, 0, l) >= -13 .and. bc_type(2, -1)%sf(k, 0, l) <= -3) then
                            ${QBMM_BC("k,-j,l,q,i","k,0,l,q,i")}$
                        elseif (bc_type(2, -1)%sf(k, 0, l) == -2) then
                            ${QBMM_BC("k,-j,l,q,i","k,j-1,l,q,i")}$
                        elseif (bc_type(2, -1)%sf(k, 0, l) == -1) then
                            ${QBMM_BC("k,-j,l,q,i","k,n - (j-1),l,q,i")}$
                        end if
                    end do
                end do
            end if

            if (bcye < 0) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = 0, p
                    do k = -buff_size, m + buff_size
                        if (bc_type(2, 1)%sf(k, 0, l) >= -13 .and. bc_type(2, 1)%sf(k, 0, l) <= -3) then
                            ${QBMM_BC("k,n+j,l,q,i","k,n,l,q,i")}$
                        elseif (bc_type(2, 1)%sf(k, 0, l) == -2) then
                            ${QBMM_BC("k,n+j,l,q,i","k,n - (j-1),l,q,i")}$
                        elseif (bc_type(2, 1)%sf(k, 0, l) == -1) then
                            ${QBMM_BC("k,n+j,l,q,i","k,j-1,k,q,i")}$
                        end if
                    end do
                end do
            end if
        end if

        if (p == 0) return

        !< z-direction
        if (bczb >= 0) then
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 3, -1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = -buff_size, n + buff_size
                do k = -buff_size, m + buff_size
                    if (bc_type(3, -1)%sf(k, l, 0) >= -13 .and. bc_type(3, -1)%sf(k, l, 0) <= -3) then
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("k,l,-j","k,l,0")}$
                    elseif (bc_type(3, -1)%sf(k, l, 0) == -2) then
                        ${PRIM_SYMMETRY_BC(3,"k,l,-j","k,l,j-1")}$
                    elseif (bc_type(3, -1)%sf(k, l, 0) == -1) then
                        ${PRIM_PERIODIC_BC("k,l,-j","k,l,p-(j-1)")}$
                    elseif (bc_type(3, -1)%sf(k, l, 0) == -15) then
                        ${PRIM_SLIP_WALL_BC("z","L")}$
                    elseif (bc_type(3, -1)%sf(k, l, 0) == -16) then
                        ${PRIM_NO_SLIP_WALL_BC("z","L")}$
                    elseif (bc_type(3, -1)%sf(k, l, 0) == -17) then
#ifdef MFC_PRE_PROCESS
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("k,l,-j","k,l,0")}$
#else
                        ${PRIM_DIRICHLET_BC(3,-1,"k,l,-j","k,l,i")}$
#endif
                    end if
                end do
            end do
        end if

        if (bcze >= 0) then
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 3, 1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = -buff_size, n + buff_size
                do k = -buff_size, m + buff_size
                    if (bc_type(3, 1)%sf(k, l, 0) >= -13 .and. bc_type(3, 1)%sf(k, l, 0) <= -3) then
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("k,l,p+j","k,l,p")}$
                    elseif (bc_type(3, 1)%sf(k, l, 0) == -2) then
                        ${PRIM_SYMMETRY_BC(3,"k,l,p+j","k,l,p - (j-1)")}$
                    elseif (bc_type(3, 1)%sf(k, l, 0) == -1) then
                        ${PRIM_PERIODIC_BC("k,l,p+j","k,l,j-1")}$
                    elseif (bc_type(3, 1)%sf(k, l, 0) == -15) then
                        ${PRIM_SLIP_WALL_BC("z","R")}$
                    elseif (bc_type(3, 1)%sf(k, l, 0) == -16) then
                        ${PRIM_NO_SLIP_WALL_BC("z","R")}$
                    elseif (bc_type(3, 1)%sf(k, l, 0) == -17) then
#ifdef MFC_PRE_PROCESS
                        ${PRIM_GHOST_CELL_EXTRAPOLATION_BC("k,l,p+j","k,l,p")}$
#else
                        ${PRIM_DIRICHLET_BC(3,1,"k,l,p+j","k,l,i")}$
#endif
                    end if
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            if (bczb < 0) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = -buff_size, n + buff_size
                    do k = -buff_size, m + buff_size
                        if (bc_type(3, -1)%sf(k, l, 0) >= -13 .and. bc_type(3, -1)%sf(k, l, 0) <= -3) then
                            ${QBMM_BC("k,l,-j,q,i","k,l,0,q,i")}$
                        elseif (bc_type(3, -1)%sf(k, l, 0) == -2) then
                            ${QBMM_BC("k,l,-j,q,i","k,l,j-1,q,i")}$
                        elseif (bc_type(3, -1)%sf(k, l, 0) == -1) then
                            ${QBMM_BC("k,l,-j,q,i","k,l,p - (j-1),q,i")}$
                        end if
                    end do
                end do
            end if

            if (bcze < 0) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = -buff_size, n + buff_size
                    do k = -buff_size, m + buff_size
                        if (bc_type(3, 1)%sf(k, l, 0) >= -13 .and. bc_type(3, 1)%sf(k, l, 0) <= -3) then
                            ${QBMM_BC("k,l,p+j,q,i","k,l,p,q,i")}$
                        elseif (bc_type(3, 1)%sf(k, l, 0) == -2) then
                            ${QBMM_BC("k,l,p+j,q,i","k,l,p - (j-1),q,i")}$
                        elseif (bc_type(3, 1)%sf(k, l, 0) == -1) then
                            ${QBMM_BC("k,l,p+j,q,i","k,l,j-1,q,i")}$
                        end if
                    end do
                end do
            end if
        end if

    end subroutine s_populate_variables_buffers

#ifdef MFC_SIMULATION
    subroutine s_populate_capillary_buffers(c_divs, bc_type)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        integer :: i, j, k, l

        !< x-direction
        if (bcxb >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 1, -1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    if (bc_type(1, -1)%sf(0, k, l) == -1) then
                        ${COLOR_FUNC_EXTRAPOLATION("-j,k,l","m - (j-1),k,l")}$
                    elseif (bc_type(1, -1)%sf(0, k, l) == -2) then
                        ${COLOR_FUNC_SLIP_WALL_BC(1,"-j,k,l","j-1,k,l")}$
                    else
                        ${COLOR_FUNC_EXTRAPOLATION("-j,k,l","0,k,l")}$
                    end if
                end do
            end do
        end if

        if (bcxe >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 1, 1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    if (bc_type(1, 1)%sf(0, k, l) == -1) then
                        ${COLOR_FUNC_EXTRAPOLATION("m+j,k,l","j-1,k,l")}$
                    elseif (bc_type(1, 1)%sf(0, k, l) == -2) then
                        ${COLOR_FUNC_SLIP_WALL_BC(1,"m+j,k,l","m - (j-1),k,l")}$
                    else
                        ${COLOR_FUNC_EXTRAPOLATION("m+j,k,l","m,k,l")}$
                    end if
                end do
            end do
        end if

        if (n == 0) return

        !< y-direction
        if (bcyb >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 2, -1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = -buff_size, m + buff_size
                    if (bc_type(2, -1)%sf(k, 0, l) == -1) then
                        ${COLOR_FUNC_EXTRAPOLATION("k,-j,l","k,n - (j-1),l")}$
                    elseif (bc_type(2, -1)%sf(k, 0, l) == -2) then
                        ${COLOR_FUNC_SLIP_WALL_BC(2,"k,-j,l","k,j-1,l")}$
                    else
                        ${COLOR_FUNC_EXTRAPOLATION("k,-j,l","k,0,l")}$
                    end if
                end do
            end do
        end if

        if (bcye >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 2, 1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = -buff_size, m + buff_size
                    if (bc_type(2, 1)%sf(k, 0, l) == -1) then
                        ${COLOR_FUNC_EXTRAPOLATION("k,n+j,l","k,j-1,l")}$
                    elseif (bc_type(2, 1)%sf(k, 0, l) == -2) then
                        ${COLOR_FUNC_SLIP_WALL_BC(2,"k,n+j,l","k,n - (j-1),l")}$
                    else
                        ${COLOR_FUNC_EXTRAPOLATION("k,n+j,l","k,n,l")}$
                    end if
                end do
            end do
        end if

        if (p == 0) return

        !< z-direction
        if (bczb >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 3, -1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = -buff_size, n + buff_size
                do k = -buff_size, m + buff_size
                    if (bc_type(3, -1)%sf(k, l, 0) == -1) then
                        ${COLOR_FUNC_EXTRAPOLATION("k,l,-j","k,l,p - (j-1)")}$
                    elseif (bc_type(3, -1)%sf(k, l, 0) == -2) then
                        ${COLOR_FUNC_SLIP_WALL_BC(3,"k,l,-j","k,l,j-1")}$
                    else
                        ${COLOR_FUNC_EXTRAPOLATION("k,l,-j","k,l,0")}$
                    end if
                end do
            end do
        end if

        if (bcze >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 3, 1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = -buff_size, n + buff_size
                do k = -buff_size, m + buff_size
                    if (bc_type(3, 1)%sf(k, l, 0) == -1) then
                        ${COLOR_FUNC_EXTRAPOLATION("k,l,p+j","k,l,j-1")}$
                    elseif (bc_type(3, 1)%sf(k, l, 0) == -2) then
                        ${COLOR_FUNC_SLIP_WALL_BC(3,"k,l,p+j","k,l,p - (j-1)")}$
                    else
                        ${COLOR_FUNC_EXTRAPOLATION("k,l,p+j","k,l,p")}$
                    end if
                end do
            end do
        end if

    end subroutine s_populate_capillary_buffers
#endif

    subroutine s_create_mpi_types(bc_type)

        type(integer_field), dimension(1:num_dims, -1:1) :: bc_type

#ifdef MFC_MPI
        integer :: dir, loc
        integer, dimension(3) :: sf_start_idx, sf_extents_loc
        integer :: ifile, ierr, data_size

        do dir = 1, num_dims
            do loc = -1, 1, 2
                sf_start_idx = (/0, 0, 0/)
                sf_extents_loc = shape(bc_type(dir, loc)%sf)

                call MPI_TYPE_CREATE_SUBARRAY(num_dims, sf_extents_loc, sf_extents_loc, sf_start_idx, &
                                              MPI_ORDER_FORTRAN, MPI_INTEGER, MPI_BC_TYPE_TYPE(dir, loc), ierr)
                call MPI_TYPE_COMMIT(MPI_BC_TYPE_TYPE(dir, loc), ierr)
            end do
        end do

        do dir = 1, num_dims
            do loc = -1, 1, 2
                sf_start_idx = (/0, 0, 0/)
                sf_extents_loc = shape(bc_buffers(dir, loc)%sf)

                call MPI_TYPE_CREATE_SUBARRAY(num_dims, sf_extents_loc, sf_extents_loc, sf_start_idx, &
                                              MPI_ORDER_FORTRAN, mpi_p, MPI_BC_BUFFER_TYPE(dir, loc), ierr)
                call MPI_TYPE_COMMIT(MPI_BC_BUFFER_TYPE(dir, loc), ierr)
            end do
        end do
#endif
    end subroutine s_create_mpi_types

    subroutine s_finalize_boundary_common_module()

#ifndef MFC_POST_PROCESS

#ifdef MFC_PRE_PROCESS
        if (save_bc) then
#elif MFC_SIMULATION
        if (read_bc) then
#endif
            deallocate (bc_buffers(1, -1)%sf)
            deallocate (bc_buffers(1, 1)%sf)
            if (n > 0) then
                deallocate (bc_buffers(2, -1)%sf)
                deallocate (bc_buffers(2, 1)%sf)
                if (p > 0) then
                    deallocate (bc_buffers(3, -1)%sf)
                    deallocate (bc_buffers(3, 1)%sf)
                end if
            end if
        end if

        deallocate (bc_buffers)
#endif

    end subroutine s_finalize_boundary_common_module

end module m_boundary_common
