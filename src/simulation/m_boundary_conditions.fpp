!>
!! @file m_boundary_conditions.fpp
!! @brief Contains module m_boundary_conditions

!> @brief The purpose of the module is to apply noncharacteristic and processor
!! boundary condiitons
module m_boundary_conditions

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy

    use m_constants
    ! ==========================================================================

    implicit none

    private; 
    public :: s_populate_primitive_variables_buffers, &
              s_populate_capillary_buffers

contains

    !>  The purpose of this procedure is to populate the buffers
    !!      of the primitive variables, depending on the selected
    !!      boundary conditions.
    !! @param q_prim_vf Primitive variable
    subroutine s_populate_primitive_variables_buffers(q_prim_vf, pb, mv)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv

        integer :: bc_loc, bc_dir

        ! Population of Buffers in x-direction =============================

        select case (bc_x%beg)
        case (-13:-3) ! Ghost-cell extrap. BC at beginning
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, -1)
        case (-2)     ! Symmetry BC at beginning
            call s_symmetry(q_prim_vf, pb, mv, 1, -1)
        case (-1)     ! Periodic BC at beginning
            call s_periodic(q_prim_vf, pb, mv, 1, -1)
        case (-15)    ! Slip wall BC at beginning
            call s_slip_wall(q_prim_vf, pb, mv, 1, -1)
        case (-16)    ! No-slip wall BC at beginning
            call s_no_slip_wall(q_prim_vf, pb, mv, 1, -1)
        case default ! Processor BC at beginning
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 1, -1)
        end select

        select case (bc_x%end)
        case (-13:-3) ! Ghost-cell extrap. BC at end
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 1, 1)
        case (-2)     ! Symmetry BC at end
            call s_symmetry(q_prim_vf, pb, mv, 1, 1)
        case (-1)     ! Periodic BC at end
            call s_periodic(q_prim_vf, pb, mv, 1, 1)
        case (-15)    ! Slip wall BC at end
            call s_slip_wall(q_prim_vf, pb, mv, 1, 1)
        case (-16)    ! No-slip wall bc at end
            call s_no_slip_wall(q_prim_vf, pb, mv, 1, 1)
        case default ! Processor BC at end
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 1, 1)
        end select

        if (qbmm .and. .not. polytropic) then
            select case (bc_x%beg)
            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                call s_qbmm_extrapolation(pb, mv, 1, -1)
            case (-15)    ! Slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 1, -1)
            case (-16)    ! No-slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 1, -1)
            end select

            select case (bc_x%end)
            case (-13:-3) ! Ghost-cell extrap. BC at end
                call s_qbmm_extrapolation(pb, mv, 1, 1)
            case (-15)    ! Slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 1, 1)
            case (-16)    ! No-slip wall bc at end
                call s_qbmm_extrapolation(pb, mv, 1, 1)
            end select
        end if

        ! END: Population of Buffers in x-direction ========================

        ! Population of Buffers in y-direction =============================

        if (n == 0) return

        select case (bc_y%beg)
        case (-13:-3) ! Ghost-cell extrap. BC at beginning
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, -1)
        case (-14)    ! Axis BC at beginning
            call s_axis(q_prim_vf, pb, mv, 2, -1)
        case (-2)     ! Symmetry BC at beginning
            call s_symmetry(q_prim_vf, pb, mv, 2, -1)
        case (-1)     ! Periodic BC at beginning
            call s_periodic(q_prim_vf, pb, mv, 2, -1)
        case (-15)    ! Slip wall BC at beginning
            call s_slip_wall(q_prim_vf, pb, mv, 2, -1)
        case (-16)    ! No-slip wall BC at beginning
            call s_no_slip_wall(q_prim_vf, pb, mv, 2, -1)
        case default ! Processor BC at beginning
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 2, -1)
        end select

        select case (bc_y%end)
        case (-13:-3) ! Ghost-cell extrap. BC at end
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 2, 1)
        case (-2)     ! Symmetry BC at end
            call s_symmetry(q_prim_vf, pb, mv, 2, 1)
        case (-1)     ! Periodic BC at end
            call s_periodic(q_prim_vf, pb, mv, 2, 1)
        case (-15)    ! Slip wall BC at end
            call s_slip_wall(q_prim_vf, pb, mv, 2, 1)
        case (-16)    ! No-slip wall BC at end
            call s_no_slip_wall(q_prim_vf, pb, mv, 2, 1)
        case default ! Processor BC at end
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 2, 1)
        end select

        if (qbmm .and. .not. polytropic) then

            select case (bc_y%beg)
            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                call s_qbmm_extrapolation(pb, mv, 2, -1)
            case (-15)    ! Slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 2, -1)
            case (-16)    ! No-slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 2, -1)
            end select

            select case (bc_y%end)
            case (-13:-3) ! Ghost-cell extrap. BC at end
                call s_qbmm_extrapolation(pb, mv, 2, 1)
            case (-15)    ! Slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 2, 1)
            case (-16)    ! No-slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 2, 1)
            end select

        end if

        ! END: Population of Buffers in y-direction ========================

        ! Population of Buffers in z-direction =============================

        if (p == 0) return

        select case (bc_z%beg)
        case (-13:-3) ! Ghost-cell extrap. BC at beginning
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, -1)
        case (-2)     ! Symmetry BC at beginning
            call s_symmetry(q_prim_vf, pb, mv, 3, -1)
        case (-1)     ! Periodic BC at beginning
            call s_periodic(q_prim_vf, pb, mv, 3, -1)
        case (-15)    ! Slip wall BC at beginning
            call s_slip_wall(q_prim_vf, pb, mv, 3, -1)
        case (-16)    ! No-slip wall BC at beginning
            call s_no_slip_wall(q_prim_vf, pb, mv, 3, -1)
        case default ! Processor BC at beginning
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 3, -1)
        end select

        select case (bc_z%end)
        case (-13:-3) ! Ghost-cell extrap. BC at end
            call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, 3, 1)
        case (-2)     ! Symmetry BC at end
            call s_symmetry(q_prim_vf, pb, mv, 3, 1)
        case (-1)     ! Periodic BC at end
            call s_periodic(q_prim_vf, pb, mv, 3, 1)
        case (-15)    ! Slip wall BC at end
            call s_slip_wall(q_prim_vf, pb, mv, 3, 1)
        case (-16)    ! No-slip wall BC at end
            call s_no_slip_wall(q_prim_vf, pb, mv, 3, 1)
        case default ! Processor BC at end
            call s_mpi_sendrecv_variables_buffers( &
                q_prim_vf, pb, mv, 3, 1)
        end select

        if (qbmm .and. .not. polytropic) then

            select case (bc_z%beg)
            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                call s_qbmm_extrapolation(pb, mv, 3, -1)
            case (-15)    ! Slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 3, -1)
            case (-16)    ! No-slip wall BC at beginning
                call s_qbmm_extrapolation(pb, mv, 3, -1)
            end select

            select case (bc_z%end)
            case (-13:-3) ! Ghost-cell extrap. BC at end
                call s_qbmm_extrapolation(pb, mv, 3, 1)
            case (-15)    ! Slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 3, 1)
            case (-16)    ! No-slip wall BC at end
                call s_qbmm_extrapolation(pb, mv, 3, 1)
            end select

        end if

        ! END: Population of Buffers in z-direction ========================

    end subroutine s_populate_primitive_variables_buffers

    subroutine s_ghost_cell_extrapolation(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc
        integer :: j, k, l, q, i

        !< x-direction =========================================================
        if (bc_dir == 1) then !< x-direction

            if (bc_loc == -1) then !bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do
                end do

            else !< bc_x%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do
                end do

            end if

            !< y-direction =========================================================
        elseif (bc_dir == 2) then !< y-direction

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, 0, k)
                            end do
                        end do
                    end do
                end do

            else !< bc_y%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(l, n + j, k) = &
                                    q_prim_vf(i)%sf(l, n, k)
                            end do
                        end do
                    end do
                end do

            end if

            !< z-direction =========================================================
        elseif (bc_dir == 3) then !< z-direction

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, 0)
                            end do
                        end do
                    end do
                end do

            else !< bc_z%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, p)
                            end do
                        end do
                    end do
                end do

            end if

        end if
        !< =====================================================================

    end subroutine s_ghost_cell_extrapolation

    subroutine s_symmetry(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction =========================================================
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, p
                    do k = 0, n
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
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 1, buff_size
                                        pb(-j, k, l, q, i) = &
                                            pb(j - 1, k, l, q, i)
                                        mv(-j, k, l, q, i) = &
                                            mv(j - 1, k, l, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            else !< bc_x%end

                !$acc parallel loop collapse(3) default(present)
                do l = 0, p
                    do k = 0, n
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

                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 1, buff_size
                                        pb(m + j, k, l, q, i) = &
                                            pb(m - (j - 1), k, l, q, i)
                                        mv(m + j, k, l, q, i) = &
                                            mv(m - (j - 1), k, l, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            end if

            !< y-direction =========================================================
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(3) gang vector default(present)
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
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
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do k = 0, p
                                do j = 1, buff_size
                                    do l = -buff_size, m + buff_size
                                        pb(l, -j, k, q, i) = &
                                            pb(l, j - 1, k, q, i)
                                        mv(l, -j, k, q, i) = &
                                            mv(l, j - 1, k, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            else !< bc_y%end

                !$acc parallel loop collapse(3) gang vector default(present)
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
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
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do k = 0, p
                                do j = 1, buff_size
                                    do l = -buff_size, m + buff_size
                                        pb(l, n + j, k, q, i) = &
                                            pb(l, n - (j - 1), k, q, i)
                                        mv(l, n + j, k, q, i) = &
                                            mv(l, n - (j - 1), k, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            end if

            !< z-direction =========================================================
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(3) gang vector default(present)
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
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
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                do l = -buff_size, n + buff_size
                                    do k = -buff_size, m + buff_size
                                        pb(k, l, -j, q, i) = &
                                            pb(k, l, j - 1, q, i)
                                        mv(k, l, -j, q, i) = &
                                            mv(k, l, j - 1, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            else !< bc_z%end

                !$acc parallel loop collapse(3) gang vector default(present)
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
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
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                do l = -buff_size, n + buff_size
                                    do k = -buff_size, m + buff_size
                                        pb(k, l, p + j, q, i) = &
                                            pb(k, l, p - (j - 1), q, i)
                                        mv(k, l, p + j, q, i) = &
                                            mv(k, l, p - (j - 1), q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            end if

        end if
        !< =====================================================================

    end subroutine s_symmetry

    subroutine s_periodic(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction =========================================================
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(-j, k, l) = &
                                    q_prim_vf(i)%sf(m - (j - 1), k, l)
                            end do
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 1, buff_size
                                        pb(-j, k, l, q, i) = &
                                            pb(m - (j - 1), k, l, q, i)
                                        mv(-j, k, l, q, i) = &
                                            mv(m - (j - 1), k, l, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            else !< bc_x%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_prim_vf(i)%sf(m + j, k, l) = &
                                    q_prim_vf(i)%sf(j - 1, k, l)
                            end do
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 1, buff_size
                                        pb(m + j, k, l, q, i) = &
                                            pb(j - 1, k, l, q, i)
                                        mv(m + j, k, l, q, i) = &
                                            mv(j - 1, k, l, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            end if

            !< y-direction =========================================================
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(l, -j, k) = &
                                    q_prim_vf(i)%sf(l, n - (j - 1), k)
                            end do
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do k = 0, p
                                do j = 1, buff_size
                                    do l = -buff_size, m + buff_size
                                        pb(l, -j, k, q, i) = &
                                            pb(l, n - (j - 1), k, q, i)
                                        mv(l, -j, k, q, i) = &
                                            mv(l, n - (j - 1), k, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            else !< bc_y%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(l, n + j, k) = &
                                    q_prim_vf(i)%sf(l, j - 1, k)
                            end do
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do k = 0, p
                                do j = 1, buff_size
                                    do l = -buff_size, m + buff_size
                                        pb(l, n + j, k, q, i) = &
                                            pb(l, (j - 1), k, q, i)
                                        mv(l, n + j, k, q, i) = &
                                            mv(l, (j - 1), k, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            end if

            !< z-direction =========================================================
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(k, l, -j) = &
                                    q_prim_vf(i)%sf(k, l, p - (j - 1))
                            end do
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                do l = -buff_size, n + buff_size
                                    do k = -buff_size, m + buff_size
                                        pb(k, l, -j, q, i) = &
                                            pb(k, l, p - (j - 1), q, i)
                                        mv(k, l, -j, q, i) = &
                                            mv(k, l, p - (j - 1), q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            else !< bc_z%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                q_prim_vf(i)%sf(k, l, p + j) = &
                                    q_prim_vf(i)%sf(k, l, j - 1)
                            end do
                        end do
                    end do
                end do

                if (qbmm .and. .not. polytropic) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                do l = -buff_size, n + buff_size
                                    do k = -buff_size, m + buff_size
                                        pb(k, l, p + j, q, i) = &
                                            pb(k, l, j - 1, q, i)
                                        mv(k, l, p + j, q, i) = &
                                            mv(k, l, j - 1, q, i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if

            end if

        end if
        !< =====================================================================

    end subroutine s_periodic

    subroutine s_axis(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

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

        if (qbmm .and. .not. polytropic) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do q = 1, nnode
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                pb(l, -j, k, q, i) = &
                                    pb(l, j - 1, k - ((p + 1)/2), q, i)
                                mv(l, -j, k, q, i) = &
                                    mv(l, j - 1, k - ((p + 1)/2), q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_axis

    subroutine s_slip_wall(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction =========================================================
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2d0*bc_x%vb1
                                else
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        q_prim_vf(i)%sf(0, k, l)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_x%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2d0*bc_x%ve1
                                else
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        q_prim_vf(i)%sf(m, k, l)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

            !< y-direction =========================================================
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                if (i == momxb + 1) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2d0*bc_y%vb2
                                else
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        q_prim_vf(i)%sf(l, 0, k)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_y%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                if (i == momxb + 1) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2d0*bc_y%ve2
                                else
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        q_prim_vf(i)%sf(l, n, k)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

            !< z-direction =========================================================
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                if (i == momxe) then
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        -q_prim_vf(i)%sf(k, l, j - 1) + 2d0*bc_z%vb3
                                else
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        q_prim_vf(i)%sf(k, l, 0)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_z%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                if (i == momxe) then
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2d0*bc_z%ve3
                                else
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        q_prim_vf(i)%sf(k, l, p)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

        end if
        !< =====================================================================

    end subroutine s_slip_wall

    subroutine s_no_slip_wall(q_prim_vf, pb, mv, bc_dir, bc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction =========================================================
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2d0*bc_x%vb1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2d0*bc_x%vb2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        -q_prim_vf(i)%sf(j - 1, k, l) + 2d0*bc_x%vb3
                                else
                                    q_prim_vf(i)%sf(-j, k, l) = &
                                        q_prim_vf(i)%sf(0, k, l)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_x%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2d0*bc_x%ve1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2d0*bc_x%ve2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2d0*bc_x%ve3
                                else
                                    q_prim_vf(i)%sf(m + j, k, l) = &
                                        q_prim_vf(i)%sf(m, k, l)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

            !< y-direction =========================================================
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2d0*bc_y%vb1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2d0*bc_y%vb2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        -q_prim_vf(i)%sf(l, j - 1, k) + 2d0*bc_y%vb3
                                else
                                    q_prim_vf(i)%sf(l, -j, k) = &
                                        q_prim_vf(i)%sf(l, 0, k)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_y%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do k = 0, p
                        do j = 1, buff_size
                            do l = -buff_size, m + buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2d0*bc_y%ve1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2d0*bc_y%ve2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        -q_prim_vf(i)%sf(l, n - (j - 1), k) + 2d0*bc_y%ve3
                                else
                                    q_prim_vf(i)%sf(l, n + j, k) = &
                                        q_prim_vf(i)%sf(l, n, k)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

            !< z-direction =========================================================
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        -q_prim_vf(i)%sf(k, l, j - 1) + 2d0*bc_z%vb1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        -q_prim_vf(i)%sf(k, l, j - 1) + 2d0*bc_z%vb2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        -q_prim_vf(i)%sf(k, l, j - 1) + 2d0*bc_z%vb3
                                else
                                    q_prim_vf(i)%sf(k, l, -j) = &
                                        q_prim_vf(i)%sf(k, l, 0)
                                end if
                            end do
                        end do
                    end do
                end do

            else !< bc_z%end

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = 1, buff_size
                        do l = -buff_size, n + buff_size
                            do k = -buff_size, m + buff_size
                                if (i == momxb) then
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2d0*bc_z%ve1
                                elseif (i == momxb + 1 .and. num_dims > 1) then
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2d0*bc_z%ve2
                                elseif (i == momxb + 2 .and. num_dims > 2) then
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2d0*bc_z%ve3
                                else
                                    q_prim_vf(i)%sf(k, l, p + j) = &
                                        q_prim_vf(i)%sf(k, l, p)
                                end if
                            end do
                        end do
                    end do
                end do

            end if

        end if
        !< =====================================================================

    end subroutine s_no_slip_wall

    subroutine s_qbmm_extrapolation(pb, mv, bc_dir, bc_loc)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: bc_dir, bc_loc

        integer :: j, k, l, q, i

        !< x-direction =========================================================
        if (bc_dir == 1) then

            if (bc_loc == -1) then !< bc_x%beg

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    pb(-j, k, l, q, i) = &
                                        pb(0, k, l, q, i)
                                    mv(-j, k, l, q, i) = &
                                        mv(0, k, l, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            else !< bc_x%end

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do l = 0, p
                            do k = 0, n
                                do j = 1, buff_size
                                    pb(m + j, k, l, q, i) = &
                                        pb(m, k, l, q, i)
                                    mv(m + j, k, l, q, i) = &
                                        mv(m, k, l, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            end if

            !< y-direction =========================================================
        elseif (bc_dir == 2) then

            if (bc_loc == -1) then !< bc_y%beg

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    pb(l, -j, k, q, i) = &
                                        pb(l, 0, k, q, i)
                                    mv(l, -j, k, q, i) = &
                                        mv(l, 0, k, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            else !< bc_y%end

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do k = 0, p
                            do j = 1, buff_size
                                do l = -buff_size, m + buff_size
                                    pb(l, n + j, k, q, i) = &
                                        pb(l, n, k, q, i)
                                    mv(l, n + j, k, q, i) = &
                                        mv(l, n, k, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            end if

            !< z-direction =========================================================
        elseif (bc_dir == 3) then

            if (bc_loc == -1) then !< bc_z%beg

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    pb(k, l, -j, q, i) = &
                                        pb(k, l, 0, q, i)
                                    mv(k, l, -j, q, i) = &
                                        mv(k, l, 0, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            else !< bc_z%end

                !$acc parallel loop collapse(5) gang vector default(present)
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            do l = -buff_size, n + buff_size
                                do k = -buff_size, m + buff_size
                                    pb(k, l, p + j, q, i) = &
                                        pb(k, l, p, q, i)
                                    mv(k, l, p + j, q, i) = &
                                        mv(k, l, p, q, i)
                                end do
                            end do
                        end do
                    end do
                end do

            end if

        end if

    end subroutine

    subroutine s_populate_capillary_buffers(c_divs)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
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

end module m_boundary_conditions
