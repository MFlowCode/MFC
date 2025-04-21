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

            !>  The purpose of this procedure is to populate the buffers
    !!      of the primitive variables, depending on the selected
    !!      boundary conditions.
            subroutine s_populate_variables_buffers(q_prim_vf, pb, mv, bc_type)

                type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
                type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type
                real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv

                integer :: bc_loc, bc_dir
                integer :: l, k

                ! Population of Buffers in x-direction
                !$acc parallel loop collapse(2) gang vector default(present)
                if (bcxb >= 0) then
                    call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 1, -1)
                else
                    do l = 0, p
                        do k = 0, n
                            select case (int(bc_type(1, -1)%sf(0, k, l)))
                            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, -1, l, k)
                            case (-2)     ! Symmetry BC at beginning
                                call s_symmetry(q_prim_vf, pb, mv, 1, -1, l, k)
                            case (-1)     ! Periodic BC at beginning
                                call s_periodic(q_prim_vf, pb, mv, 1, -1, l, k)
                            case (-15)    ! Slip wall BC at beginning
                                call s_slip_wall(q_prim_vf, pb, mv, 1, -1, l, k)
                            case (-16)    ! No-slip wall BC at beginning
                                call s_no_slip_wall(q_prim_vf, pb, mv, 1, -1, l, k)
                            case (-17)
#ifdef MFC_PRE_PROCESS
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, -1, l, k)
#else
                                call s_dirichlet(q_prim_vf, pb, mv, 1, -1, l, k)
#endif
                            end select
                        end do
                    end do
                end if

                if (bcxe >= 0) then
                    call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 1, 1)
                else
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            select case (int(bc_type(1, 1)%sf(0, k, l)))
                            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, 1, l, k)
                            case (-2)     ! Symmetry BC at beginning
                                call s_symmetry(q_prim_vf, pb, mv, 1, 1, l, k)
                            case (-1)     ! Periodic BC at beginning
                                call s_periodic(q_prim_vf, pb, mv, 1, 1, l, k)
                            case (-15)    ! Slip wall BC at beginning
                                call s_slip_wall(q_prim_vf, pb, mv, 1, 1, l, k)
                            case (-16)    ! No-slip wall BC at beginning
                                call s_no_slip_wall(q_prim_vf, pb, mv, 1, 1, l, k)
                            case (-17)
                                print *, "here1"
#ifdef MFC_PRE_PROCESS
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, 1, l, k)
#else
                                call s_dirichlet(q_prim_vf, pb, mv, 1, 1, l, k)
#endif
                                print *, "here"
                            end select
                        end do
                    end do
                end if

                ! Population of Buffers in y-direction

                if (n == 0) return

                if (bcyb >= 0) then
                    call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 2, -1)
                else
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do l = 0, p
                        do k = -buff_size, m + buff_size
                            select case (int(bc_type(2, -1)%sf(k, 0, l)))
                            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, -1, l, k)
                            case (-2)     ! Symmetry BC at beginning
                                call s_axis(q_prim_vf, pb, mv, 2, -1, l, k)
                            case (-1)     ! Periodic BC at beginning
                                call s_periodic(q_prim_vf, pb, mv, 2, -1, l, k)
                            case (-15)    ! Slip wall BC at beginning
                                call s_slip_wall(q_prim_vf, pb, mv, 2, -1, l, k)
                            case (-16)    ! No-slip wall BC at beginning
                                call s_no_slip_wall(q_prim_vf, pb, mv, 2, -1, l, k)
                            case (-17)
#ifdef MFC_PRE_PROCESS
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, -1, l, k)
#else
                                call s_dirichlet(q_prim_vf, pb, mv, 2, -1, l, k)
#endif
                            case default ! Processor BC at beginning
                                call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 2, -1)
                            end select
                        end do
                    end do
                end if

                if (bcye >= 0) then
                    call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 2, 1)
                else
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do l = 0, p
                        do k = -buff_size, m + buff_size
                            select case (int(bc_type(2, 1)%sf(k, 0, l)))
                            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, 1, l, k)
                            case (-2)     ! Symmetry BC at beginning
                                call s_symmetry(q_prim_vf, pb, mv, 2, 1, l, k)
                            case (-1)     ! Periodic BC at beginning
                                call s_periodic(q_prim_vf, pb, mv, 2, 1, l, k)
                            case (-15)    ! Slip wall BC at beginning
                                call s_slip_wall(q_prim_vf, pb, mv, 2, 1, l, k)
                            case (-16)    ! No-slip wall BC at beginning
                                call s_no_slip_wall(q_prim_vf, pb, mv, 2, 1, l, k)
                            case (-17)
#ifdef MFC_PRE_PROCESS
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, 1, l, k)
#else
                                call s_dirichlet(q_prim_vf, pb, mv, 2, 1, l, k)
#endif
                            end select
                        end do
                    end do
                end if

                ! Population of Buffers in z-direction

                if (p == 0) return
                if (bczb >= 0) then
                    call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 3, -1)
                else
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            select case (int(bc_type(3, -1)%sf(k, l, 0)))
                            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, -1, l, k)
                            case (-2)     ! Symmetry BC at beginning
                                call s_symmetry(q_prim_vf, pb, mv, 3, -1, l, k)
                            case (-1)     ! Periodic BC at beginning
                                call s_periodic(q_prim_vf, pb, mv, 3, -1, l, k)
                            case (-15)    ! Slip wall BC at beginning
                                call s_slip_wall(q_prim_vf, pb, mv, 3, -1, l, k)
                            case (-16)    ! No-slip wall BC at beginning
                                call s_no_slip_wall(q_prim_vf, pb, mv, 3, -1, l, k)
                            case (-17)
#ifdef MFC_PRE_PROCESS
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, -1, l, k)
#else
                                call s_dirichlet(q_prim_vf, pb, mv, 3, -1, l, k)
#endif
                            end select
                        end do
                    end do
                end if

                if (bcze >= 0) then
                    call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 3, 1)
                else
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            select case (int(bc_type(3, 1)%sf(k, l, 0)))
                            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, 1, l, k)
                            case (-2)     ! Symmetry BC at beginning
                                call s_symmetry(q_prim_vf, pb, mv, 3, 1, l, k)
                            case (-1)     ! Periodic BC at beginning
                                call s_periodic(q_prim_vf, pb, mv, 3, 1, l, k)
                            case (-15)    ! Slip wall BC at beginning
                                call s_slip_wall(q_prim_vf, pb, mv, 3, 1, l, k)
                            case (-16)    ! No-slip wall BC at beginning
                                call s_no_slip_wall(q_prim_vf, pb, mv, 3, 1, l, k)
                            case (-17)
#ifdef MFC_PRE_PROCESS
                                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, 1, l, k)
#else
                                call s_dirichlet(q_prim_vf, pb, mv, 3, 1, l, k)
#endif
                            end select
                        end do
                    end do
                end if
                ! END: Population of Buffers in z-direction

            end subroutine s_populate_variables_buffers

            subroutine s_ghost_cell_extrapolation(q_prim_vf, pb, mv, bc_dir, bc_loc, l, k)

                type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
                real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
                integer, intent(in) :: bc_dir, bc_loc
                integer, intent(in) :: l, k

                integer :: j, q, i

                !< x-direction
                if (bc_dir == 1) then !< x-direction

                    if (bc_loc == -1) then !bc_x%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(0, k, l)
                            end do
                        end do

                    else !< bc_x%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(m, k, l)
                            end do
                        end do

                    end if

                    !< y-direction
                elseif (bc_dir == 2) then !< y-direction

                    if (bc_loc == -1) then !< bc_y%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, -j, l) = &
                                    q_prim_vf(i)%sf(k, 0, l)
                            end do
                        end do

                    else !< bc_y%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, n + j, l) = &
                                    q_prim_vf(i)%sf(k, n, l)
                            end do
                        end do

                    end if

                    !< z-direction
                elseif (bc_dir == 3) then !< z-direction

                    if (bc_loc == -1) then !< bc_z%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, 0)
                            end do
                        end do

                    else !< bc_z%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, p)
                            end do
                        end do

                    end if

                end if

            end subroutine s_ghost_cell_extrapolation

            subroutine s_symmetry(q_prim_vf, pb, mv, bc_dir, bc_loc, l, k)

                type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
                real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
                integer, intent(in) :: bc_dir, bc_loc
                integer, intent(in) :: l, k

                integer :: j, q, i

                !< x-direction
                if (bc_dir == 1) then

                    if (bc_loc == -1) then !< bc_x%beg

                        !$acc parallel loop gang vector default(present)
                        do j = 1, buff_size
                            !$acc loop seq
                            do i = 1, contxe
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(j - 1, k, l)
                            end do

                            q_prim_vf(momxb)%sf(-j, k, l) = &
                                -q_prim_vf(momxb)%sf(j - 1, k, l)

                            !$acc loop seq
                            do i = momxb + 1, sys_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(j - 1, k, l)
                            end do

                            if (elasticity) then
                                do i = 1, shear_BC_flip_num
                                    q_prim_vf(shear_BC_flip_indices(1, i))%sf(-j, k, l) = &
                                        -q_prim_vf(shear_BC_flip_indices(1, i))%sf(j - 1, k, l)
                                end do
                            end if

                            if (hyperelasticity) then
                                q_prim_vf(xibeg)%sf(-j, k, l) = &
                                    -q_prim_vf(xibeg)%sf(j - 1, k, l)
                            end if

                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(-j, k, l, q, i) = &
                                            pb(j - 1, k, l, q, i)
                                        mv(-j, k, l, q, i) = &
                                            mv(j - 1, k, l, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif

                    else !< bc_x%end

                        !$acc parallel loop default(present)
                        do j = 1, buff_size
                            !$acc loop seq
                            do i = 1, contxe
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(m - (j - 1), k, l)
                            end do

                            q_prim_vf(momxb)%sf(m + j, k, l) = &
                                -q_prim_vf(momxb)%sf(m - (j - 1), k, l)

                            !$acc loop seq
                            do i = momxb + 1, sys_size
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(m - (j - 1), k, l)
                            end do

                            if (elasticity) then
                                do i = 1, shear_BC_flip_num
                                    q_prim_vf(shear_BC_flip_indices(1, i))%sf(m + j, k, l) = &
                                        -q_prim_vf(shear_BC_flip_indices(1, i))%sf(m - (j - 1), k, l)
                                end do
                            end if

                            if (hyperelasticity) then
                                q_prim_vf(xibeg)%sf(m + j, k, l) = &
                                    -q_prim_vf(xibeg)%sf(m - (j - 1), k, l)
                            end if

                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(m + j, k, l, q, i) = &
                                            pb(m - (j - 1), k, l, q, i)
                                        mv(m + j, k, l, q, i) = &
                                            mv(m - (j - 1), k, l, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    end if

                    !< y-direction
                elseif (bc_dir == 2) then

                    if (bc_loc == -1) then !< bc_y%beg

                        !$acc parallel loop gang vector default(present)
                        do j = 1, buff_size
                            !$acc loop seq
                            do i = 1, momxb
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k)
                            end do

                            q_prim_vf(momxb + 1)%sf(l, -j, k) = &
                                -q_prim_vf(momxb + 1)%sf(l, j - 1, k)

                            !$acc loop seq
                            do i = momxb + 2, sys_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k)
                            end do

                            if (elasticity) then
                                do i = 1, shear_BC_flip_num
                                    q_prim_vf(shear_BC_flip_indices(2, i))%sf(l, -j, k) = &
                                        -q_prim_vf(shear_BC_flip_indices(2, i))%sf(l, j - 1, k)
                                end do
                            end if

                            if (hyperelasticity) then
                                q_prim_vf(xibeg + 1)%sf(l, -j, k) = &
                                    -q_prim_vf(xibeg + 1)%sf(l, j - 1, k)
                            end if
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(l, -j, k, q, i) = &
                                            pb(l, j - 1, k, q, i)
                                        mv(l, -j, k, q, i) = &
                                            mv(l, j - 1, k, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    else !< bc_y%end

                        !$acc parallel loop gang vector default(present)
                        do j = 1, buff_size
                            !$acc loop seq
                            do i = 1, momxb
                                q_prim_vf(i)%sf(l, n + j, k) = &
                                    q_prim_vf(i)%sf(l, n - (j - 1), k)
                            end do

                            q_prim_vf(momxb + 1)%sf(l, n + j, k) = &
                                -q_prim_vf(momxb + 1)%sf(l, n - (j - 1), k)

                            !$acc loop seq
                            do i = momxb + 2, sys_size
                                q_prim_vf(i)%sf(l, n + j, k) = &
                                    q_prim_vf(i)%sf(l, n - (j - 1), k)
                            end do

                            if (elasticity) then
                                do i = 1, shear_BC_flip_num
                                    q_prim_vf(shear_BC_flip_indices(2, i))%sf(l, n + j, k) = &
                                        -q_prim_vf(shear_BC_flip_indices(2, i))%sf(l, n - (j - 1), k)
                                end do
                            end if

                            if (hyperelasticity) then
                                q_prim_vf(xibeg + 1)%sf(l, n + j, k) = &
                                    -q_prim_vf(xibeg + 1)%sf(l, n - (j - 1), k)
                            end if
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(l, n + j, k, q, i) = &
                                            pb(l, n - (j - 1), k, q, i)
                                        mv(l, n + j, k, q, i) = &
                                            mv(l, n - (j - 1), k, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    end if

                    !< z-direction
                elseif (bc_dir == 3) then

                    if (bc_loc == -1) then !< bc_z%beg

                        !$acc parallel loop gang vector default(present)
                        do j = 1, buff_size
                            !$acc loop seq
                            do i = 1, momxb + 1
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, j - 1)
                            end do

                            q_prim_vf(momxe)%sf(k, l, -j) = &
                                -q_prim_vf(momxe)%sf(k, l, j - 1)

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, j - 1)
                            end do

                            if (elasticity) then
                                do i = 1, shear_BC_flip_num
                                    q_prim_vf(shear_BC_flip_indices(3, i))%sf(k, l, -j) = &
                                        -q_prim_vf(shear_BC_flip_indices(3, i))%sf(k, l, j - 1)
                                end do
                            end if

                            if (hyperelasticity) then
                                q_prim_vf(xiend)%sf(k, l, -j) = &
                                    -q_prim_vf(xiend)%sf(k, l, j - 1)
                            end if
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(k, l, -j, q, i) = &
                                            pb(k, l, j - 1, q, i)
                                        mv(k, l, -j, q, i) = &
                                            mv(k, l, j - 1, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    else !< bc_z%end

                        !$acc parallel loop gang vector default(present)
                        do j = 1, buff_size
                            !$acc loop seq
                            do i = 1, momxb + 1
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, p - (j - 1))
                            end do

                            q_prim_vf(momxe)%sf(k, l, p + j) = &
                                -q_prim_vf(momxe)%sf(k, l, p - (j - 1))

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, p - (j - 1))
                            end do

                            if (elasticity) then
                                do i = 1, shear_BC_flip_num
                                    q_prim_vf(shear_BC_flip_indices(3, i))%sf(k, l, p + j) = &
                                        -q_prim_vf(shear_BC_flip_indices(3, i))%sf(k, l, p - (j - 1))
                                end do
                            end if

                            if (hyperelasticity) then
                                q_prim_vf(xiend)%sf(k, l, p + j) = &
                                    -q_prim_vf(xiend)%sf(k, l, p - (j - 1))
                            end if
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(k, l, p + j, q, i) = &
                                            pb(k, l, p - (j - 1), q, i)
                                        mv(k, l, p + j, q, i) = &
                                            mv(k, l, p - (j - 1), q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    end if

                end if

            end subroutine s_symmetry

            subroutine s_periodic(q_prim_vf, pb, mv, bc_dir, bc_loc, l, k)

                type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
                real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
                integer, intent(in) :: bc_dir, bc_loc
                integer, intent(in) :: l, k

                integer :: j, q, i

                !< x-direction
                if (bc_dir == 1) then

                    if (bc_loc == -1) then !< bc_x%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(m - (j - 1), k, l)
                            end do
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(-j, k, l, q, i) = &
                                            pb(m - (j - 1), k, l, q, i)
                                        mv(-j, k, l, q, i) = &
                                            mv(m - (j - 1), k, l, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    else !< bc_x%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(j - 1, k, l)
                            end do
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(m + j, k, l, q, i) = &
                                            pb(j - 1, k, l, q, i)
                                        mv(m + j, k, l, q, i) = &
                                            mv(j - 1, k, l, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    end if

                    !< y-direction
                elseif (bc_dir == 2) then

                    if (bc_loc == -1) then !< bc_y%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, n - (j - 1), k)
                            end do
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(l, -j, k, q, i) = &
                                            pb(l, n - (j - 1), k, q, i)
                                        mv(l, -j, k, q, i) = &
                                            mv(l, n - (j - 1), k, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    else !< bc_y%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(l, n + j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k)
                            end do
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(l, n + j, k, q, i) = &
                                            pb(l, (j - 1), k, q, i)
                                        mv(l, n + j, k, q, i) = &
                                            mv(l, (j - 1), k, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    end if

                    !< z-direction
                elseif (bc_dir == 3) then

                    if (bc_loc == -1) then !< bc_z%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, p - (j - 1))
                            end do
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(k, l, -j, q, i) = &
                                            pb(k, l, p - (j - 1), q, i)
                                        mv(k, l, -j, q, i) = &
                                            mv(k, l, p - (j - 1), q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    else !< bc_z%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, j - 1)
                            end do
                        end do
#ifdef MFC_SIMULATION
                        if (qbmm .and. .not. polytropic) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = 1, nb
                                do q = 1, nnode
                                    do j = 1, buff_size
                                        pb(k, l, p + j, q, i) = &
                                            pb(k, l, j - 1, q, i)
                                        mv(k, l, p + j, q, i) = &
                                            mv(k, l, j - 1, q, i)
                                    end do
                                end do
                            end do
                        end if
#endif
                    end if

                end if

            end subroutine s_periodic

            subroutine s_axis(q_prim_vf, pb, mv, bc_dir, bc_loc, l, k)

                type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
                real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
                integer, intent(in) :: bc_dir, bc_loc
                integer, intent(in) :: l, k

                integer :: j, q, i

                !$acc parallel loop gang vector default(present)
                do j = 1, buff_size
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
#ifdef MFC_SIMULATION
                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb(l, -j, k, q, i) = &
                                    pb(l, j - 1, k - ((p + 1)/2), q, i)
                                mv(l, -j, k, q, i) = &
                                    mv(l, j - 1, k - ((p + 1)/2), q, i)
                            end do
                        end do
                    end do
                end if
#endif
            end subroutine s_axis

            subroutine s_slip_wall(q_prim_vf, pb, mv, bc_dir, bc_loc, l, k)

                type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
                real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
                integer, intent(in) :: bc_dir, bc_loc
                integer, intent(in) :: l, k

                integer :: j, q, i

                !< x-direction
                if (bc_dir == 1) then

                    if (bc_loc == -1) then !< bc_x%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb1
                                else
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        q_prim_vf(i)%sf(0, k, l)
                                end if
                            end do
                        end do

                    else !< bc_x%end

                        !$acc parallel loop collapse(2) gang vector default(present)
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

                    end if

                    !< y-direction
                elseif (bc_dir == 2) then

                    if (bc_loc == -1) then !< bc_y%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                if (i == momxb + 1) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2._wp*bc_y%vb2
                                else
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        q_prim_vf(i)%sf(l, 0, k)
                                end if
                            end do
                        end do

                    else !< bc_y%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                if (i == momxb + 1) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2._wp*bc_y%ve2
                                else
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        q_prim_vf(i)%sf(l, n, k)
                                end if
                            end do
                        end do

                    end if

                    !< z-direction
                elseif (bc_dir == 3) then

                    if (bc_loc == -1) then !< bc_z%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
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

                    else !< bc_z%end

                        !$acc parallel loop collapse(4) gang vector default(present)
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

                    end if

                end if

            end subroutine s_slip_wall

            subroutine s_no_slip_wall(q_prim_vf, pb, mv, bc_dir, bc_loc, l, k)

                type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
                real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
                integer, intent(in) :: bc_dir, bc_loc
                integer, intent(in) :: l, k

                integer :: j, q, i

                !< x-direction
                if (bc_dir == 1) then

                    if (bc_loc == -1) then !< bc_x%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb3
                                else
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        q_prim_vf(i)%sf(0, k, l)
                                end if
                            end do
                        end do

                    else !< bc_x%end

                        !$acc parallel loop collapse(2) gang vector default(present)
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

                    end if

                    !< y-direction
                elseif (bc_dir == 2) then

                    if (bc_loc == -1) then !< bc_y%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2._wp*bc_y%vb1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2._wp*bc_y%vb2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2._wp*bc_y%vb3
                                else
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        q_prim_vf(i)%sf(l, 0, k)
                                end if
                            end do
                        end do

                    else !< bc_y%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2._wp*bc_y%ve1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2._wp*bc_y%ve2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2._wp*bc_y%ve3
                                else
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        q_prim_vf(i)%sf(l, n, k)
                                end if
                            end do
                        end do

                    end if

                    !< z-direction
                elseif (bc_dir == 3) then

                    if (bc_loc == -1) then !< bc_z%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
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

                    else !< bc_z%end

                        !$acc parallel loop collapse(2) gang vector default(present)
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

                    end if

                end if

            end subroutine s_no_slip_wall

            subroutine s_qbmm_extrapolation(pb, mv, bc_dir, bc_loc, l, k)

                real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
                integer, intent(in) :: bc_dir, bc_loc
                integer, intent(in) :: l, k

                integer :: j, q, i

                !< x-direction
                if (bc_dir == 1) then

                    if (bc_loc == -1) then !< bc_x%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, nb
                            do q = 1, nnode
                                do j = 1, buff_size
                                    pb(-j, k, l, q, i) = &
                                        pb(0, k, l, q, i)
                                    mv(-j, k, l, q, i) = &
                                        mv(0, k, l, q, i)
                                end do
                            end do
                        end do

                    else !< bc_x%end

                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = 1, nb
                            do q = 1, nnode
                                do j = 1, buff_size
                                    pb(m + j, k, l, q, i) = &
                                        pb(m, k, l, q, i)
                                    mv(m + j, k, l, q, i) = &
                                        mv(m, k, l, q, i)
                                end do
                            end do
                        end do

                    end if

                    !< y-direction
                elseif (bc_dir == 2) then

                    if (bc_loc == -1) then !< bc_y%beg

                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = 1, nb
                            do q = 1, nnode
                                do j = 1, buff_size
                                    pb(l, -j, k, q, i) = &
                                        pb(l, 0, k, q, i)
                                    mv(l, -j, k, q, i) = &
                                        mv(l, 0, k, q, i)
                                end do
                            end do
                        end do

                    else !< bc_y%end

                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = 1, nb
                            do q = 1, nnode
                                do j = 1, buff_size
                                    pb(l, n + j, k, q, i) = &
                                        pb(l, n, k, q, i)
                                    mv(l, n + j, k, q, i) = &
                                        mv(l, n, k, q, i)
                                end do
                            end do
                        end do

                    end if

                    !< z-direction
                elseif (bc_dir == 3) then

                    if (bc_loc == -1) then !< bc_z%beg

                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = 1, nb
                            do q = 1, nnode
                                do j = 1, buff_size
                                    pb(k, l, -j, q, i) = &
                                        pb(k, l, 0, q, i)
                                    mv(k, l, -j, q, i) = &
                                        mv(k, l, 0, q, i)
                                end do
                            end do
                        end do

                    else !< bc_z%end

                        !$acc parallel loop collapse(5) gang vector default(present)
                        do i = 1, nb
                            do q = 1, nnode
                                do j = 1, buff_size
                                    pb(k, l, p + j, q, i) = &
                                        pb(k, l, p, q, i)
                                    mv(k, l, p + j, q, i) = &
                                        mv(k, l, p, q, i)
                                end do
                            end do
                        end do

                    end if

                end if

            end subroutine s_qbmm_extrapolation

            subroutine s_dirichlet(q_prim_vf, pb, mv, bc_dir, bc_loc, l, k)

                type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
                real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
                integer, intent(in) :: bc_dir, bc_loc
                integer, intent(in) :: l, k

                integer :: j, i, q

                !< x-direction
                if (bc_dir == 1) then !< x-direction

                    if (bc_loc == -1) then !bc_x%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    bc_buffers(1, -1)%sf(i, k, l)
                            end do
                        end do

                    else !< bc_x%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    bc_buffers(1, 1)%sf(i, k, l)
                            end do
                        end do

                    end if

                    !< y-direction
                elseif (bc_dir == 2) then !< y-direction

                    if (bc_loc == -1) then !< bc_y%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, -j, l) = &
                                    bc_buffers(2, -1)%sf(k, i, l)
                            end do
                        end do

                    else !< bc_y%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, n + j, l) = &
                                    bc_buffers(2, 1)%sf(k, i, l)
                            end do
                        end do

                    end if

                    !< z-direction
                elseif (bc_dir == 3) then !< z-direction

                    if (bc_loc == -1) then !< bc_z%beg

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    bc_buffers(3, -1)%sf(k, l, i)
                            end do
                        end do

                    else !< bc_z%end

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do i = 1, sys_size
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    bc_buffers(3, 1)%sf(k, l, i)
                            end do
                        end do

                    end if

                end if
            end subroutine s_dirichlet

#ifdef MFC_SIMULATION
            subroutine s_populate_capillary_buffers(c_divs, bc_type)

                type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
                type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type
                integer :: i, j, k, l

                ! x - direction
                if (bc_x%beg <= -3) then !< ghost cell extrapolation
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    c_divs(i)%sf(-j, k, l) = &
                                        c_divs(i)%sf(0, k, l)
                                end do
                            end do
                        end do
                    end do
                elseif (bc_x%beg == -2) then !< slip wall or reflective
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    if (i == 1) then
                                        c_divs(i)%sf(-j, k, l) = &
                                            -c_divs(i)%sf(j - 1, k, l)
                                    else
                                        c_divs(i)%sf(-j, k, l) = &
                                            c_divs(i)%sf(j - 1, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do
                elseif (bc_x%beg == -1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    c_divs(i)%sf(-j, k, l) = &
                                        c_divs(i)%sf(m - (j - 1), k, l)
                                end do
                            end do
                        end do
                    end do
                else
                    call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 1, -1)
                end if

                if (bc_x%end <= -3) then !< ghost-cell extrapolation
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    c_divs(i)%sf(m + j, k, l) = &
                                        c_divs(i)%sf(m, k, l)
                                end do
                            end do
                        end do
                    end do
                elseif (bc_x%end == -2) then
                    !$acc parallel loop collapse(4) default(present)
                    do i = 1, num_dims + 1
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    if (i == 1) then
                                        c_divs(i)%sf(m + j, k, l) = &
                                            -c_divs(i)%sf(m - (j - 1), k, l)
                                    else
                                        c_divs(i)%sf(m + j, k, l) = &
                                            c_divs(i)%sf(m - (j - 1), k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do
                else if (bc_x%end == -1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    c_divs(i)%sf(m + j, k, l) = &
                                        c_divs(i)%sf(j - 1, k, l)
                                end do
                            end do
                        end do
                    end do
                else
                    call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 1, 1)
                end if

                if (n == 0) then
                    return
                elseif (bc_y%beg <= -3) then !< ghost-cell extrapolation
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    c_divs(i)%sf(l, -j, k) = &
                                        c_divs(i)%sf(l, 0, k)
                                end do
                            end do
                        end do
                    end do
                elseif (bc_y%beg == -2) then !< slip wall or reflective
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    if (i == 2) then
                                        c_divs(i)%sf(l, -j, k) = &
                                            -c_divs(i)%sf(l, j - 1, k)
                                    else
                                        c_divs(i)%sf(l, -j, k) = &
                                            c_divs(i)%sf(l, j - 1, k)
                                    end if
                                end do
                            end do
                        end do
                    end do
                elseif (bc_y%beg == -1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    c_divs(i)%sf(l, -j, k) = &
                                        c_divs(i)%sf(l, n - (j - 1), k)
                                end do
                            end do
                        end do
                    end do
                else
                    call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 2, -1)
                end if

                if (bc_y%end <= -3) then !< ghost-cell extrapolation
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    c_divs(i)%sf(l, n + j, k) = &
                                        c_divs(i)%sf(l, n, k)
                                end do
                            end do
                        end do
                    end do
                elseif (bc_y%end == -2) then !< slip wall or reflective
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    if (i == 2) then
                                        c_divs(i)%sf(l, n + j, k) = &
                                            -c_divs(i)%sf(l, n - (j - 1), k)
                                    else
                                        c_divs(i)%sf(l, n + j, k) = &
                                            c_divs(i)%sf(l, n - (j - 1), k)
                                    end if
                                end do
                            end do
                        end do
                    end do
                elseif (bc_y%end == -1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    c_divs(i)%sf(l, n + j, k) = &
                                        c_divs(i)%sf(l, j - 1, k)
                                end do
                            end do
                        end do
                    end do
                else
                    call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 2, 1)
                end if

                if (p == 0) then
                    return
                elseif (bc_z%beg <= -3) then !< ghost-cell extrapolation
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    c_divs(i)%sf(k, l, -j) = &
                                        c_divs(i)%sf(k, l, 0)
                                end do
                            end do
                        end do
                    end do
                elseif (bc_z%beg == -2) then !< symmetry
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    if (i == 3) then
                                        c_divs(i)%sf(k, l, -j) = &
                                            -c_divs(i)%sf(k, l, j - 1)
                                    else
                                        c_divs(i)%sf(k, l, -j) = &
                                            c_divs(i)%sf(k, l, j - 1)
                                    end if
                                end do
                            end do
                        end do
                    end do
                elseif (bc_z%beg == -1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    c_divs(i)%sf(k, l, -j) = &
                                        c_divs(i)%sf(k, l, p - (j - 1))
                                end do
                            end do
                        end do
                    end do
                else
                    call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 3, -1)
                end if

                if (bc_z%end <= -3) then !< ghost-cell extrapolation
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    c_divs(i)%sf(k, l, p + j) = &
                                        c_divs(i)%sf(k, l, p)
                                end do
                            end do
                        end do
                    end do
                elseif (bc_z%end == -2) then !< symmetry
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    if (i == 3) then
                                        c_divs(i)%sf(k, l, p + j) = &
                                            -c_divs(i)%sf(k, l, p - (j - 1))
                                    else
                                        c_divs(i)%sf(k, l, p + j) = &
                                            c_divs(i)%sf(k, l, p - (j - 1))
                                    end if
                                end do
                            end do
                        end do
                    end do
                elseif (bc_z%end == -1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, num_dims + 1
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    c_divs(i)%sf(k, l, p + j) = &
                                        c_divs(i)%sf(k, l, j - 1)
                                end do
                            end do
                        end do
                    end do
                else
                    call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 3, 1)
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
