!>
!! @file m_boundary_conditions_common.fpp
!! @brief Contains module m_boundary_conditions_common

!> @brief The purpose of the module is to apply noncharacteristic and processor
!! boundary condiitons

#:include 'macros.fpp'

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

    private; public :: s_initialize_boundary_common_module, &
 s_populate_variables_buffers, &
 s_create_mpi_types, &
 s_populate_capillary_buffers, &
 s_finalize_boundary_common_module

    public :: bc_buffers, bcxb, bcxe, bcyb, bcye, bczb, bcze

#ifdef MFC_MPI
    public :: MPI_BC_TYPE_TYPE, MPI_BC_BUFFER_TYPE
#endif

contains

    impure subroutine s_initialize_boundary_common_module()

        bcxb = bc_x%beg; bcxe = bc_x%end; bcyb = bc_y%beg; bcye = bc_y%end; bczb = bc_z%beg; bcze = bc_z%end

        @:ALLOCATE(bc_buffers(1:num_dims, -1:1))

#ifndef MFC_POST_PROCESS
        if (bc_io) then
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
    impure subroutine s_populate_variables_buffers(q_prim_vf, pb, mv, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        integer :: k, l

        ! Population of Buffers in x-direction
        if (bcxb >= 0) then
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, pb, mv, 1, -1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    select case (int(bc_type(1, -1)%sf(0, k, l)))
                    case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP)
                        call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, -1, k, l)
                    case (BC_REFLECTIVE)
                        call s_symmetry(q_prim_vf, pb, mv, 1, -1, k, l)
                    case (BC_PERIODIC)
                        call s_periodic(q_prim_vf, pb, mv, 1, -1, k, l)
                    case (BC_SLIP_WALL)
                        call s_slip_wall(q_prim_vf, pb, mv, 1, -1, k, l)
                    case (BC_NO_SLIP_WALL)
                        call s_no_slip_wall(q_prim_vf, pb, mv, 1, -1, k, l)
                    case (BC_DIRICHLET)
                        call s_dirichlet(q_prim_vf, pb, mv, 1, -1, k, l)
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
                    case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP) ! Ghost-cell extrap. BC at end
                        call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, 1, k, l)
                    case (BC_REFLECTIVE)
                        call s_symmetry(q_prim_vf, pb, mv, 1, 1, k, l)
                    case (BC_PERIODIC)
                        call s_periodic(q_prim_vf, pb, mv, 1, 1, k, l)
                    case (BC_SLIP_WALL)
                        call s_slip_wall(q_prim_vf, pb, mv, 1, 1, k, l)
                    case (BC_NO_SLIP_WALL)
                        call s_no_slip_wall(q_prim_vf, pb, mv, 1, 1, k, l)
                    case (BC_DIRICHLET)
                        call s_dirichlet(q_prim_vf, pb, mv, 1, 1, k, l)
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
                    case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP)
                        call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, -1, k, l)
                    case (BC_AXIS)
                        call s_axis(q_prim_vf, pb, mv, 2, -1, k, l)
                    case (BC_REFLECTIVE)
                        call s_symmetry(q_prim_vf, pb, mv, 2, -1, k, l)
                    case (BC_PERIODIC)
                        call s_periodic(q_prim_vf, pb, mv, 2, -1, k, l)
                    case (BC_SLIP_WALL)
                        call s_slip_wall(q_prim_vf, pb, mv, 2, -1, k, l)
                    case (BC_NO_SLIP_WALL)
                        call s_no_slip_wall(q_prim_vf, pb, mv, 2, -1, k, l)
                    case (BC_DIRICHLET)
                        call s_dirichlet(q_prim_vf, pb, mv, 2, -1, k, l)
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
                    case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP)
                        call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, 1, k, l)
                    case (BC_REFLECTIVE)
                        call s_symmetry(q_prim_vf, pb, mv, 2, 1, k, l)
                    case (BC_PERIODIC)
                        call s_periodic(q_prim_vf, pb, mv, 2, 1, k, l)
                    case (BC_SLIP_WALL)
                        call s_slip_wall(q_prim_vf, pb, mv, 2, 1, k, l)
                    case (BC_NO_SLIP_WALL)
                        call s_no_slip_wall(q_prim_vf, pb, mv, 2, 1, k, l)
                    case (BC_DIRICHLET)
                        call s_dirichlet(q_prim_vf, pb, mv, 2, 1, k, l)
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
                    case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP)
                        call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, -1, k, l)
                    case (BC_REFLECTIVE)
                        call s_symmetry(q_prim_vf, pb, mv, 3, -1, k, l)
                    case (BC_PERIODIC)
                        call s_periodic(q_prim_vf, pb, mv, 3, -1, k, l)
                    case (BC_SLIP_WALL)
                        call s_slip_wall(q_prim_vf, pb, mv, 3, -1, k, l)
                    case (BC_NO_SLIP_WALL)
                        call s_no_slip_wall(q_prim_vf, pb, mv, 3, -1, k, l)
                    case (BC_DIRICHLET)
                        call s_dirichlet(q_prim_vf, pb, mv, 3, -1, k, l)
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
                    case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP)
                        call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, 1, k, l)
                    case (BC_REFLECTIVE)
                        call s_symmetry(q_prim_vf, pb, mv, 3, 1, k, l)
                    case (BC_PERIODIC)
                        call s_periodic(q_prim_vf, pb, mv, 3, 1, k, l)
                    case (BC_SlIP_WALL)
                        call s_slip_wall(q_prim_vf, pb, mv, 3, 1, k, l)
                    case (BC_NO_SLIP_WALL)
                        call s_no_slip_wall(q_prim_vf, pb, mv, 3, 1, k, l)
                    case (BC_DIRICHLET)
                        call s_dirichlet(q_prim_vf, pb, mv, 3, 1, k, l)
                    end select
                end do
            end do
        end if
        ! END: Population of Buffers in z-direction

    end subroutine s_populate_variables_buffers

    pure subroutine s_ghost_cell_extrapolation(q_prim_vf, pb, mv, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_ghost_cell_extrapolation
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, i

        if (qbmm .and. .not. polytropic) then
            call s_qbmm_extrapolation(pb, mv, bc_dir, bc_loc, k, l)
        end if

        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !bc_x%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(-j, k, l) = &
                            q_prim_vf(i)%sf(0, k, l)
                    end do
                end do
            else !< bc_x%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            q_prim_vf(i)%sf(m, k, l)
                    end do
                end do
            end if
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, -j, l) = &
                            q_prim_vf(i)%sf(k, 0, l)
                    end do
                end do
            else !< bc_y%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            q_prim_vf(i)%sf(k, n, l)
                    end do
                end do
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, -j) = &
                            q_prim_vf(i)%sf(k, l, 0)
                    end do
                end do
            else !< bc_z%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            q_prim_vf(i)%sf(k, l, p)
                    end do
                end do
            end if
        end if

    end subroutine s_ghost_cell_extrapolation

    pure subroutine s_symmetry(q_prim_vf, pb, mv, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_symmetry
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, q, i

        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !< bc_x%beg
                do j = 1, buff_size
                    do i = 1, contxe
                        q_prim_vf(i)%sf(-j, k, l) = &
                            q_prim_vf(i)%sf(j - 1, k, l)
                    end do

                    q_prim_vf(momxb)%sf(-j, k, l) = &
                        -q_prim_vf(momxb)%sf(j - 1, k, l)

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

                if (qbmm .and. .not. polytropic) then
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
            else !< bc_x%end
                do j = 1, buff_size
                    do i = 1, contxe
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            q_prim_vf(i)%sf(m - (j - 1), k, l)
                    end do

                    q_prim_vf(momxb)%sf(m + j, k, l) = &
                        -q_prim_vf(momxb)%sf(m - (j - 1), k, l)

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
                if (qbmm .and. .not. polytropic) then

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
            end if
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
                do j = 1, buff_size
                    do i = 1, momxb
                        q_prim_vf(i)%sf(k, -j, l) = &
                            q_prim_vf(i)%sf(k, j - 1, l)
                    end do

                    q_prim_vf(momxb + 1)%sf(k, -j, l) = &
                        -q_prim_vf(momxb + 1)%sf(k, j - 1, l)

                    do i = momxb + 2, sys_size
                        q_prim_vf(i)%sf(k, -j, l) = &
                            q_prim_vf(i)%sf(k, j - 1, l)
                    end do

                    if (elasticity) then
                        do i = 1, shear_BC_flip_num
                            q_prim_vf(shear_BC_flip_indices(2, i))%sf(k, -j, l) = &
                                -q_prim_vf(shear_BC_flip_indices(2, i))%sf(k, j - 1, l)
                        end do
                    end if

                    if (hyperelasticity) then
                        q_prim_vf(xibeg + 1)%sf(k, -j, l) = &
                            -q_prim_vf(xibeg + 1)%sf(k, j - 1, l)
                    end if
                end do

                if (qbmm .and. .not. polytropic) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb(k, -j, l, q, i) = &
                                    pb(k, j - 1, l, q, i)
                                mv(k, -j, l, q, i) = &
                                    mv(k, j - 1, l, q, i)
                            end do
                        end do
                    end do
                end if
            else !< bc_y%end
                do j = 1, buff_size
                    do i = 1, momxb
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            q_prim_vf(i)%sf(k, n - (j - 1), l)
                    end do

                    q_prim_vf(momxb + 1)%sf(k, n + j, l) = &
                        -q_prim_vf(momxb + 1)%sf(k, n - (j - 1), l)

                    do i = momxb + 2, sys_size
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            q_prim_vf(i)%sf(k, n - (j - 1), l)
                    end do

                    if (elasticity) then
                        do i = 1, shear_BC_flip_num
                            q_prim_vf(shear_BC_flip_indices(2, i))%sf(k, n + j, l) = &
                                -q_prim_vf(shear_BC_flip_indices(2, i))%sf(k, n - (j - 1), l)
                        end do
                    end if

                    if (hyperelasticity) then
                        q_prim_vf(xibeg + 1)%sf(k, n + j, l) = &
                            -q_prim_vf(xibeg + 1)%sf(k, n - (j - 1), l)
                    end if
                end do

                if (qbmm .and. .not. polytropic) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb(k, n + j, l, q, i) = &
                                    pb(k, n - (j - 1), l, q, i)
                                mv(k, n + j, l, q, i) = &
                                    mv(k, n - (j - 1), l, q, i)
                            end do
                        end do
                    end do
                end if
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
                do j = 1, buff_size
                    do i = 1, momxb + 1
                        q_prim_vf(i)%sf(k, l, -j) = &
                            q_prim_vf(i)%sf(k, l, j - 1)
                    end do

                    q_prim_vf(momxe)%sf(k, l, -j) = &
                        -q_prim_vf(momxe)%sf(k, l, j - 1)

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

                if (qbmm .and. .not. polytropic) then
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
            else !< bc_z%end
                do j = 1, buff_size
                    do i = 1, momxb + 1
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            q_prim_vf(i)%sf(k, l, p - (j - 1))
                    end do

                    q_prim_vf(momxe)%sf(k, l, p + j) = &
                        -q_prim_vf(momxe)%sf(k, l, p - (j - 1))

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

                if (qbmm .and. .not. polytropic) then
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
            end if
        end if

    end subroutine s_symmetry

    pure subroutine s_periodic(q_prim_vf, pb, mv, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_periodic
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, q, i

        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !< bc_x%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(-j, k, l) = &
                            q_prim_vf(i)%sf(m - (j - 1), k, l)
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
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
            else !< bc_x%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            q_prim_vf(i)%sf(j - 1, k, l)
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
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
            end if
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, -j, l) = &
                            q_prim_vf(i)%sf(k, n - (j - 1), l)
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb(k, -j, l, q, i) = &
                                    pb(k, n - (j - 1), l, q, i)
                                mv(k, -j, l, q, i) = &
                                    mv(k, n - (j - 1), l, q, i)
                            end do
                        end do
                    end do
                end if
            else !< bc_y%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            q_prim_vf(i)%sf(k, j - 1, l)
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb(k, n + j, l, q, i) = &
                                    pb(k, (j - 1), l, q, i)
                                mv(k, n + j, l, q, i) = &
                                    mv(k, (j - 1), l, q, i)
                            end do
                        end do
                    end do
                end if
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, -j) = &
                            q_prim_vf(i)%sf(k, l, p - (j - 1))
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
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
            else !< bc_z%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            q_prim_vf(i)%sf(k, l, j - 1)
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
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
            end if
        end if

    end subroutine s_periodic

    pure subroutine s_axis(q_prim_vf, pb, mv, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_axis
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, q, i

        do j = 1, buff_size
            if (z_cc(l) < pi) then
                do i = 1, momxb
                    q_prim_vf(i)%sf(k, -j, l) = &
                        q_prim_vf(i)%sf(k, j - 1, l + ((p + 1)/2))
                end do

                q_prim_vf(momxb + 1)%sf(k, -j, l) = &
                    -q_prim_vf(momxb + 1)%sf(k, j - 1, l + ((p + 1)/2))

                q_prim_vf(momxe)%sf(k, -j, l) = &
                    -q_prim_vf(momxe)%sf(k, j - 1, l + ((p + 1)/2))

                do i = E_idx, sys_size
                    q_prim_vf(i)%sf(k, -j, l) = &
                        q_prim_vf(i)%sf(k, j - 1, l + ((p + 1)/2))
                end do
            else
                do i = 1, momxb
                    q_prim_vf(i)%sf(k, -j, l) = &
                        q_prim_vf(i)%sf(k, j - 1, l - ((p + 1)/2))
                end do

                q_prim_vf(momxb + 1)%sf(k, -j, l) = &
                    -q_prim_vf(momxb + 1)%sf(k, j - 1, l - ((p + 1)/2))

                q_prim_vf(momxe)%sf(k, -j, l) = &
                    -q_prim_vf(momxe)%sf(k, j - 1, l - ((p + 1)/2))

                do i = E_idx, sys_size
                    q_prim_vf(i)%sf(k, -j, l) = &
                        q_prim_vf(i)%sf(k, j - 1, l - ((p + 1)/2))
                end do
            end if
        end do

        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do q = 1, nnode
                    do j = 1, buff_size
                        pb(k, -j, l, q, i) = &
                            pb(k, j - 1, l - ((p + 1)/2), q, i)
                        mv(k, -j, l, q, i) = &
                            mv(k, j - 1, l - ((p + 1)/2), q, i)
                    end do
                end do
            end do
        end if

    end subroutine s_axis

    pure subroutine s_slip_wall(q_prim_vf, pb, mv, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_slip_wall
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, i

        if (qbmm .and. .not. polytropic) then
            call s_qbmm_extrapolation(pb, mv, bc_dir, bc_loc, k, l)
        end if

        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !< bc_x%beg
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
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
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
            else !< bc_y%end
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
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
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

    pure subroutine s_no_slip_wall(q_prim_vf, pb, mv, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_no_slip_wall
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, i

        if (qbmm .and. .not. polytropic) then
            call s_qbmm_extrapolation(pb, mv, bc_dir, bc_loc, k, l)
        end if

        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !< bc_x%beg
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
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
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
            else !< bc_y%end
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
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
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

    pure subroutine s_dirichlet(q_prim_vf, pb, mv, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_dirichlet
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, i

#ifdef MFC_PRE_PROCESS
        call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, -1, k, l)
#else
        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !bc_x%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(-j, k, l) = &
                            bc_buffers(1, -1)%sf(i, k, l)
                    end do
                end do
            else !< bc_x%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(m + j, k, l) = &
                            bc_buffers(1, 1)%sf(i, k, l)
                    end do
                end do
            end if
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, -j, l) = &
                            bc_buffers(2, -1)%sf(k, i, l)
                    end do
                end do
            else !< bc_y%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, n + j, l) = &
                            bc_buffers(2, 1)%sf(k, i, l)
                    end do
                end do
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, -j) = &
                            bc_buffers(3, -1)%sf(k, l, i)
                    end do
                end do
            else !< bc_z%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, p + j) = &
                            bc_buffers(3, 1)%sf(k, l, i)
                    end do
                end do
            end if
        end if
#endif

    end subroutine s_dirichlet

    pure subroutine s_qbmm_extrapolation(pb, mv, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_qbmm_extrapolation
#else
        !$acc routine seq
#endif
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, q, i

        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !bc_x%beg
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb(-j, k, l, q, i) = pb(0, k, l, q, i)
                            mv(-j, k, l, q, i) = mv(0, k, l, q, i)
                        end do
                    end do
                end do
            else !< bc_x%end
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb(m + j, k, l, q, i) = pb(m, k, l, q, i)
                            mv(m + j, k, l, q, i) = mv(m, k, l, q, i)
                        end do
                    end do
                end do
            end if
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb(k, -j, l, q, i) = pb(k, 0, l, q, i)
                            mv(k, -j, l, q, i) = mv(k, 0, l, q, i)
                        end do
                    end do
                end do
            else !< bc_y%end
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb(k, n + j, l, q, i) = pb(k, n, l, q, i)
                            mv(k, n + j, l, q, i) = mv(k, n, l, q, i)
                        end do
                    end do
                end do
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb(k, l, -j, q, i) = pb(k, l, 0, q, i)
                            mv(k, l, -j, q, i) = mv(k, l, 0, q, i)
                        end do
                    end do
                end do
            else !< bc_z%end
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb(k, l, p + j, q, i) = pb(k, l, p, q, i)
                            mv(k, l, p + j, q, i) = mv(k, l, p, q, i)
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_qbmm_extrapolation

    impure subroutine s_populate_capillary_buffers(c_divs, bc_type)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        integer :: k, l

        !< x-direction
        if (bcxb >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 1, -1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    select case (bc_type(1, -1)%sf(0, k, l))
                    case (BC_PERIODIC)
                        call s_color_function_periodic(c_divs, 1, -1, k, l)
                    case (BC_REFLECTIVE)
                        call s_color_function_reflective(c_divs, 1, -1, k, l)
                    case default
                        call s_color_function_ghost_cell_extrapolation(c_divs, 1, -1, k, l)
                    end select
                end do
            end do
        end if

        if (bcxe >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 1, 1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    select case (bc_type(1, 1)%sf(0, k, l))
                    case (BC_PERIODIC)
                        call s_color_function_periodic(c_divs, 1, 1, k, l)
                    case (BC_REFLECTIVE)
                        call s_color_function_reflective(c_divs, 1, 1, k, l)
                    case default
                        call s_color_function_ghost_cell_extrapolation(c_divs, 1, 1, k, l)
                    end select
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
                    select case (bc_type(2, -1)%sf(k, 0, l))
                    case (BC_PERIODIC)
                        call s_color_function_periodic(c_divs, 2, -1, k, l)
                    case (BC_REFLECTIVE)
                        call s_color_function_reflective(c_divs, 2, -1, k, l)
                    case default
                        call s_color_function_ghost_cell_extrapolation(c_divs, 2, -1, k, l)
                    end select
                end do
            end do
        end if

        if (bcye >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 2, 1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, p
                do k = -buff_size, m + buff_size
                    select case (bc_type(2, 1)%sf(k, 0, l))
                    case (BC_PERIODIC)
                        call s_color_function_periodic(c_divs, 2, 1, k, l)
                    case (BC_REFLECTIVE)
                        call s_color_function_reflective(c_divs, 2, 1, k, l)
                    case default
                        call s_color_function_ghost_cell_extrapolation(c_divs, 2, 1, k, l)
                    end select
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
                    select case (bc_type(3, -1)%sf(k, l, 0))
                    case (BC_PERIODIC)
                        call s_color_function_periodic(c_divs, 3, -1, k, l)
                    case (BC_REFLECTIVE)
                        call s_color_function_reflective(c_divs, 3, -1, k, l)
                    case default
                        call s_color_function_ghost_cell_extrapolation(c_divs, 3, -1, k, l)
                    end select
                end do
            end do
        end if

        if (bcze >= 0) then
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs, 3, 1)
        else
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = -buff_size, n + buff_size
                do k = -buff_size, m + buff_size
                    select case (bc_type(3, 1)%sf(k, l, 0))
                    case (BC_PERIODIC)
                        call s_color_function_periodic(c_divs, 3, 1, k, l)
                    case (BC_REFLECTIVE)
                        call s_color_function_reflective(c_divs, 3, 1, k, l)
                    case default
                        call s_color_function_ghost_cell_extrapolation(c_divs, 3, 1, k, l)
                    end select
                end do
            end do
        end if
    end subroutine s_populate_capillary_buffers

    pure subroutine s_color_function_periodic(c_divs, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_color_function_periodic
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, i

        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !bc_x%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(-j, k, l) = c_divs(i)%sf(m - (j - 1), k, l)
                    end do
                end do
            else !< bc_x%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(m + j, k, l) = c_divs(i)%sf(j - 1, k, l)
                    end do
                end do
            end if
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, -j, l) = c_divs(i)%sf(k, n - (j - 1), l)
                    end do
                end do
            else !< bc_y%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, n + j, l) = c_divs(i)%sf(k, j - 1, l)
                    end do
                end do
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, -j) = c_divs(i)%sf(k, l, p - (j - 1))
                    end do
                end do
            else !< bc_z%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, p + j) = c_divs(i)%sf(k, l, j - 1)
                    end do
                end do
            end if
        end if

    end subroutine s_color_function_periodic

    pure subroutine s_color_function_reflective(c_divs, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_color_function_reflective
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, i

        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !bc_x%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(-j, k, l) = -c_divs(i)%sf(j - 1, k, l)
                        else
                            c_divs(i)%sf(-j, k, l) = c_divs(i)%sf(j - 1, k, l)
                        end if
                    end do
                end do
            else !< bc_x%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(m + j, k, l) = -c_divs(i)%sf(m - (j - 1), k, l)
                        else
                            c_divs(i)%sf(m + j, k, l) = c_divs(i)%sf(m - (j - 1), k, l)
                        end if
                    end do
                end do
            end if
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, -j, l) = -c_divs(i)%sf(k, j - 1, l)
                        else
                            c_divs(i)%sf(k, -j, l) = c_divs(i)%sf(k, j - 1, l)
                        end if
                    end do
                end do
            else !< bc_y%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, n + j, l) = -c_divs(i)%sf(k, n - (j - 1), l)
                        else
                            c_divs(i)%sf(k, n + j, l) = c_divs(i)%sf(k, n - (j - 1), l)
                        end if
                    end do
                end do
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, l, -j) = -c_divs(i)%sf(k, l, j - 1)
                        else
                            c_divs(i)%sf(k, l, -j) = c_divs(i)%sf(k, l, j - 1)
                        end if
                    end do
                end do
            else !< bc_z%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, l, p + j) = -c_divs(i)%sf(k, l, p - (j - 1))
                        else
                            c_divs(i)%sf(k, l, p + j) = c_divs(i)%sf(k, l, p - (j - 1))
                        end if
                    end do
                end do
            end if
        end if

    end subroutine s_color_function_reflective

    pure subroutine s_color_function_ghost_cell_extrapolation(c_divs, bc_dir, bc_loc, k, l)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_color_function_ghost_cell_extrapolation
#else
        !$acc routine seq
#endif
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer, intent(in) :: bc_dir, bc_loc
        integer, intent(in) :: k, l

        integer :: j, i

        if (bc_dir == 1) then !< x-direction
            if (bc_loc == -1) then !bc_x%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(-j, k, l) = c_divs(i)%sf(0, k, l)
                    end do
                end do
            else !< bc_x%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(m + j, k, l) = c_divs(i)%sf(m, k, l)
                    end do
                end do
            end if
        elseif (bc_dir == 2) then !< y-direction
            if (bc_loc == -1) then !< bc_y%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, -j, l) = c_divs(i)%sf(k, 0, l)
                    end do
                end do
            else !< bc_y%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, n + j, l) = c_divs(i)%sf(k, n, l)
                    end do
                end do
            end if
        elseif (bc_dir == 3) then !< z-direction
            if (bc_loc == -1) then !< bc_z%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, -j) = c_divs(i)%sf(k, l, 0)
                    end do
                end do
            else !< bc_z%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, p + j) = c_divs(i)%sf(k, l, p)
                    end do
                end do
            end if
        end if

    end subroutine s_color_function_ghost_cell_extrapolation

    impure subroutine s_create_mpi_types(bc_type)

        type(integer_field), dimension(1:num_dims, -1:1) :: bc_type

#ifdef MFC_MPI
        integer :: dir, loc
        integer, dimension(3) :: sf_start_idx, sf_extents_loc
        integer :: ierr

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

    impure subroutine s_finalize_boundary_common_module()

#ifndef MFC_POST_PROCESS
        if (bc_io) then
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
#endif
        deallocate (bc_buffers)

    end subroutine s_finalize_boundary_common_module

end module m_boundary_common
