!>
!! @file
!! @brief Contains module m_boundary_primitives

!> @brief Per-cell noncharacteristic boundary condition primitives applied in the ghost cells
#:include 'case.fpp'
#:include 'macros.fpp'

module m_boundary_primitives

    use m_derived_types
    use m_global_parameters
    use m_constants

    implicit none

    type(scalar_field), dimension(:,:), allocatable :: bc_buffers
    $:GPU_DECLARE(create='[bc_buffers]')

contains

    !> Fill ghost cells by copying the nearest boundary cell value along the specified direction.
    subroutine s_ghost_cell_extrapolation(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)

        $:GPU_ROUTINE(function_name='s_ghost_cell_extrapolation', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in)                                    :: bc_dir, bc_loc
        integer, intent(in)                                    :: k, l
        integer                                                :: j, i
        type(scalar_field), optional, intent(inout)            :: q_T_sf

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(-j, k, l) = q_prim_vf(i)%sf(0, k, l)
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(-j, k, l) = q_T_sf%sf(0, k, l)
                    end do
                end if
            else  !< bc_x%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(m + j, k, l) = q_prim_vf(i)%sf(m, k, l)
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(m + j, k, l) = q_T_sf%sf(m, k, l)
                    end do
                end if
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, 0, l)
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(k, -j, l) = q_T_sf%sf(k, 0, l)
                    end do
                end if
            else  !< bc_y%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, n + j, l) = q_prim_vf(i)%sf(k, n, l)
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(k, n + j, l) = q_T_sf%sf(k, n, l)
                    end do
                end if
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, -j) = q_prim_vf(i)%sf(k, l, 0)
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(k, l, -j) = q_T_sf%sf(k, l, 0)
                    end do
                end if
            else  !< bc_z%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, p + j) = q_prim_vf(i)%sf(k, l, p)
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(k, l, p + j) = q_T_sf%sf(k, l, p)
                    end do
                end if
            end if
        end if

    end subroutine s_ghost_cell_extrapolation

    !> Apply reflective (symmetry) boundary conditions by mirroring primitive variables and flipping the normal velocity component.
    subroutine s_symmetry(q_prim_vf, bc_dir, bc_loc, k, l, pb_in, mv_in, q_T_sf)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout)                                               :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        integer, intent(in)                                                                                  :: bc_dir, bc_loc
        integer, intent(in)                                                                                  :: k, l
        integer                                                                                              :: j, q, i
        type(scalar_field), optional, intent(inout)                                                          :: q_T_sf

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  !< bc_x%beg
                do j = 1, buff_size
                    do i = 1, eqn_idx%cont%end
                        q_prim_vf(i)%sf(-j, k, l) = q_prim_vf(i)%sf(j - 1, k, l)
                    end do

                    q_prim_vf(eqn_idx%mom%beg)%sf(-j, k, l) = -q_prim_vf(eqn_idx%mom%beg)%sf(j - 1, k, l)

                    do i = eqn_idx%mom%beg + 1, sys_size
                        q_prim_vf(i)%sf(-j, k, l) = q_prim_vf(i)%sf(j - 1, k, l)
                    end do

                    if (chemistry .and. present(q_T_sf)) then
                        q_T_sf%sf(-j, k, l) = q_T_sf%sf(j - 1, k, l)
                    end if

                    if (elasticity) then
                        do i = 1, shear_BC_flip_num
                            q_prim_vf(shear_BC_flip_indices(1, i))%sf(-j, k, l) = -q_prim_vf(shear_BC_flip_indices(1, &
                                      & i))%sf(j - 1, k, l)
                        end do
                    end if

                    if (hyperelasticity) then
                        q_prim_vf(eqn_idx%xi%beg)%sf(-j, k, l) = -q_prim_vf(eqn_idx%xi%beg)%sf(j - 1, k, l)
                    end if
                end do

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(-j, k, l, q, i) = pb_in(j - 1, k, l, q, i)
                                mv_in(-j, k, l, q, i) = mv_in(j - 1, k, l, q, i)
                            end do
                        end do
                    end do
                end if
            else  !< bc_x%end
                do j = 1, buff_size
                    do i = 1, eqn_idx%cont%end
                        q_prim_vf(i)%sf(m + j, k, l) = q_prim_vf(i)%sf(m - (j - 1), k, l)
                    end do

                    q_prim_vf(eqn_idx%mom%beg)%sf(m + j, k, l) = -q_prim_vf(eqn_idx%mom%beg)%sf(m - (j - 1), k, l)

                    do i = eqn_idx%mom%beg + 1, sys_size
                        q_prim_vf(i)%sf(m + j, k, l) = q_prim_vf(i)%sf(m - (j - 1), k, l)
                    end do

                    if (chemistry .and. present(q_T_sf)) then
                        q_T_sf%sf(m + j, k, l) = q_T_sf%sf(m - (j - 1), k, l)
                    end if

                    if (elasticity) then
                        do i = 1, shear_BC_flip_num
                            q_prim_vf(shear_BC_flip_indices(1, i))%sf(m + j, k, l) = -q_prim_vf(shear_BC_flip_indices(1, &
                                      & i))%sf(m - (j - 1), k, l)
                        end do
                    end if

                    if (hyperelasticity) then
                        q_prim_vf(eqn_idx%xi%beg)%sf(m + j, k, l) = -q_prim_vf(eqn_idx%xi%beg)%sf(m - (j - 1), k, l)
                    end if
                end do
                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(m + j, k, l, q, i) = pb_in(m - (j - 1), k, l, q, i)
                                mv_in(m + j, k, l, q, i) = mv_in(m - (j - 1), k, l, q, i)
                            end do
                        end do
                    end do
                end if
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do j = 1, buff_size
                    do i = 1, eqn_idx%mom%beg
                        q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, j - 1, l)
                    end do

                    q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, -j, l) = -q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, j - 1, l)

                    do i = eqn_idx%mom%beg + 2, sys_size
                        q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, j - 1, l)
                    end do

                    if (chemistry .and. present(q_T_sf)) then
                        q_T_sf%sf(k, -j, l) = q_T_sf%sf(k, j - 1, l)
                    end if

                    if (elasticity) then
                        do i = 1, shear_BC_flip_num
                            q_prim_vf(shear_BC_flip_indices(2, i))%sf(k, -j, l) = -q_prim_vf(shear_BC_flip_indices(2, i))%sf(k, &
                                      & j - 1, l)
                        end do
                    end if

                    if (hyperelasticity) then
                        q_prim_vf(eqn_idx%xi%beg + 1)%sf(k, -j, l) = -q_prim_vf(eqn_idx%xi%beg + 1)%sf(k, j - 1, l)
                    end if
                end do

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(k, -j, l, q, i) = pb_in(k, j - 1, l, q, i)
                                mv_in(k, -j, l, q, i) = mv_in(k, j - 1, l, q, i)
                            end do
                        end do
                    end do
                end if
            else  !< bc_y%end
                do j = 1, buff_size
                    do i = 1, eqn_idx%mom%beg
                        q_prim_vf(i)%sf(k, n + j, l) = q_prim_vf(i)%sf(k, n - (j - 1), l)
                    end do

                    q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, n + j, l) = -q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, n - (j - 1), l)

                    do i = eqn_idx%mom%beg + 2, sys_size
                        q_prim_vf(i)%sf(k, n + j, l) = q_prim_vf(i)%sf(k, n - (j - 1), l)
                    end do

                    if (chemistry .and. present(q_T_sf)) then
                        q_T_sf%sf(k, n + j, l) = q_T_sf%sf(k, n - (j - 1), l)
                    end if

                    if (elasticity) then
                        do i = 1, shear_BC_flip_num
                            q_prim_vf(shear_BC_flip_indices(2, i))%sf(k, n + j, l) = -q_prim_vf(shear_BC_flip_indices(2, &
                                      & i))%sf(k, n - (j - 1), l)
                        end do
                    end if

                    if (hyperelasticity) then
                        q_prim_vf(eqn_idx%xi%beg + 1)%sf(k, n + j, l) = -q_prim_vf(eqn_idx%xi%beg + 1)%sf(k, n - (j - 1), l)
                    end if
                end do

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(k, n + j, l, q, i) = pb_in(k, n - (j - 1), l, q, i)
                                mv_in(k, n + j, l, q, i) = mv_in(k, n - (j - 1), l, q, i)
                            end do
                        end do
                    end do
                end if
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do j = 1, buff_size
                    do i = 1, eqn_idx%mom%beg + 1
                        q_prim_vf(i)%sf(k, l, -j) = q_prim_vf(i)%sf(k, l, j - 1)
                    end do

                    q_prim_vf(eqn_idx%mom%end)%sf(k, l, -j) = -q_prim_vf(eqn_idx%mom%end)%sf(k, l, j - 1)

                    do i = eqn_idx%E, sys_size
                        q_prim_vf(i)%sf(k, l, -j) = q_prim_vf(i)%sf(k, l, j - 1)
                    end do

                    if (chemistry .and. present(q_T_sf)) then
                        q_T_sf%sf(k, l, -j) = q_T_sf%sf(k, l, j - 1)
                    end if

                    if (elasticity) then
                        do i = 1, shear_BC_flip_num
                            q_prim_vf(shear_BC_flip_indices(3, i))%sf(k, l, -j) = -q_prim_vf(shear_BC_flip_indices(3, i))%sf(k, &
                                      & l, j - 1)
                        end do
                    end if

                    if (hyperelasticity) then
                        q_prim_vf(eqn_idx%xi%end)%sf(k, l, -j) = -q_prim_vf(eqn_idx%xi%end)%sf(k, l, j - 1)
                    end if
                end do

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(k, l, -j, q, i) = pb_in(k, l, j - 1, q, i)
                                mv_in(k, l, -j, q, i) = mv_in(k, l, j - 1, q, i)
                            end do
                        end do
                    end do
                end if
            else  !< bc_z%end
                do j = 1, buff_size
                    do i = 1, eqn_idx%mom%beg + 1
                        q_prim_vf(i)%sf(k, l, p + j) = q_prim_vf(i)%sf(k, l, p - (j - 1))
                    end do

                    q_prim_vf(eqn_idx%mom%end)%sf(k, l, p + j) = -q_prim_vf(eqn_idx%mom%end)%sf(k, l, p - (j - 1))

                    do i = eqn_idx%E, sys_size
                        q_prim_vf(i)%sf(k, l, p + j) = q_prim_vf(i)%sf(k, l, p - (j - 1))
                    end do

                    if (chemistry .and. present(q_T_sf)) then
                        q_T_sf%sf(k, l, p + j) = q_T_sf%sf(k, l, p - (j - 1))
                    end if

                    if (elasticity) then
                        do i = 1, shear_BC_flip_num
                            q_prim_vf(shear_BC_flip_indices(3, i))%sf(k, l, p + j) = -q_prim_vf(shear_BC_flip_indices(3, &
                                      & i))%sf(k, l, p - (j - 1))
                        end do
                    end if

                    if (hyperelasticity) then
                        q_prim_vf(eqn_idx%xi%end)%sf(k, l, p + j) = -q_prim_vf(eqn_idx%xi%end)%sf(k, l, p - (j - 1))
                    end if
                end do

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(k, l, p + j, q, i) = pb_in(k, l, p - (j - 1), q, i)
                                mv_in(k, l, p + j, q, i) = mv_in(k, l, p - (j - 1), q, i)
                            end do
                        end do
                    end do
                end if
            end if
        end if

    end subroutine s_symmetry

    !> Apply periodic boundary conditions by copying values from the opposite domain boundary.
    subroutine s_periodic(q_prim_vf, bc_dir, bc_loc, k, l, pb_in, mv_in, q_T_sf)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout)                                               :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        integer, intent(in)                                                                                  :: bc_dir, bc_loc
        integer, intent(in)                                                                                  :: k, l
        integer                                                                                              :: j, q, i
        type(scalar_field), optional, intent(inout)                                                          :: q_T_sf

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  !< bc_x%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(-j, k, l) = q_prim_vf(i)%sf(m - (j - 1), k, l)
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(-j, k, l) = q_T_sf%sf(m - (j - 1), k, l)
                    end do
                end if

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(-j, k, l, q, i) = pb_in(m - (j - 1), k, l, q, i)
                                mv_in(-j, k, l, q, i) = mv_in(m - (j - 1), k, l, q, i)
                            end do
                        end do
                    end do
                end if
            else  !< bc_x%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(m + j, k, l) = q_prim_vf(i)%sf(j - 1, k, l)
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(m + j, k, l) = q_T_sf%sf(j - 1, k, l)
                    end do
                end if

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(m + j, k, l, q, i) = pb_in(j - 1, k, l, q, i)
                                mv_in(m + j, k, l, q, i) = mv_in(j - 1, k, l, q, i)
                            end do
                        end do
                    end do
                end if
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, n - (j - 1), l)
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(k, -j, l) = q_T_sf%sf(k, n - (j - 1), l)
                    end do
                end if

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(k, -j, l, q, i) = pb_in(k, n - (j - 1), l, q, i)
                                mv_in(k, -j, l, q, i) = mv_in(k, n - (j - 1), l, q, i)
                            end do
                        end do
                    end do
                end if
            else  !< bc_y%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, n + j, l) = q_prim_vf(i)%sf(k, j - 1, l)
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(k, n + j, l) = q_T_sf%sf(k, j - 1, l)
                    end do
                end if

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(k, n + j, l, q, i) = pb_in(k, (j - 1), l, q, i)
                                mv_in(k, n + j, l, q, i) = mv_in(k, (j - 1), l, q, i)
                            end do
                        end do
                    end do
                end if
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, -j) = q_prim_vf(i)%sf(k, l, p - (j - 1))
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(k, l, -j) = q_T_sf%sf(k, l, p - (j - 1))
                    end do
                end if

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(k, l, -j, q, i) = pb_in(k, l, p - (j - 1), q, i)
                                mv_in(k, l, -j, q, i) = mv_in(k, l, p - (j - 1), q, i)
                            end do
                        end do
                    end do
                end if
            else  !< bc_z%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(k, l, p + j) = q_prim_vf(i)%sf(k, l, j - 1)
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(k, l, p + j) = q_T_sf%sf(k, l, j - 1)
                    end do
                end if

                if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
                    do i = 1, nb
                        do q = 1, nnode
                            do j = 1, buff_size
                                pb_in(k, l, p + j, q, i) = pb_in(k, l, j - 1, q, i)
                                mv_in(k, l, p + j, q, i) = mv_in(k, l, j - 1, q, i)
                            end do
                        end do
                    end do
                end if
            end if
        end if

    end subroutine s_periodic

    !> Apply axis boundary conditions for cylindrical coordinates by reflecting values across the axis with azimuthal phase shift.
    subroutine s_axis(q_prim_vf, pb_in, mv_in, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout)                                               :: q_prim_vf
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), optional, intent(inout) :: pb_in, mv_in
        integer, intent(in)                                                                                  :: k, l
        integer                                                                                              :: j, q, i

        do j = 1, buff_size
            if (z_cc(l) < pi) then
                do i = 1, eqn_idx%mom%beg
                    q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, j - 1, l + ((p + 1)/2))
                end do

                q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, -j, l) = -q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, j - 1, l + ((p + 1)/2))

                q_prim_vf(eqn_idx%mom%end)%sf(k, -j, l) = -q_prim_vf(eqn_idx%mom%end)%sf(k, j - 1, l + ((p + 1)/2))

                do i = eqn_idx%E, sys_size
                    q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, j - 1, l + ((p + 1)/2))
                end do
            else
                do i = 1, eqn_idx%mom%beg
                    q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, j - 1, l - ((p + 1)/2))
                end do

                q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, -j, l) = -q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, j - 1, l - ((p + 1)/2))

                q_prim_vf(eqn_idx%mom%end)%sf(k, -j, l) = -q_prim_vf(eqn_idx%mom%end)%sf(k, j - 1, l - ((p + 1)/2))

                do i = eqn_idx%E, sys_size
                    q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, j - 1, l - ((p + 1)/2))
                end do
            end if
        end do

        if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
            do i = 1, nb
                do q = 1, nnode
                    do j = 1, buff_size
                        if (z_cc(l) < pi) then
                            pb_in(k, -j, l, q, i) = pb_in(k, j - 1, l + ((p + 1)/2), q, i)
                            mv_in(k, -j, l, q, i) = mv_in(k, j - 1, l + ((p + 1)/2), q, i)
                        else
                            pb_in(k, -j, l, q, i) = pb_in(k, j - 1, l - ((p + 1)/2), q, i)
                            mv_in(k, -j, l, q, i) = mv_in(k, j - 1, l - ((p + 1)/2), q, i)
                        end if
                    end do
                end do
            end do
        end if

    end subroutine s_axis

    !> Apply slip wall boundary conditions by extrapolating scalars and reflecting the wall-normal velocity component.
    subroutine s_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)

        $:GPU_ROUTINE(function_name='s_slip_wall',parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in)                                    :: bc_dir, bc_loc
        integer, intent(in)                                    :: k, l
        integer                                                :: j, i
        type(scalar_field), optional, intent(inout)            :: q_T_sf

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  !< bc_x%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg) then
                            q_prim_vf(i)%sf(-j, k, l) = -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb1
                        else
                            q_prim_vf(i)%sf(-j, k, l) = q_prim_vf(i)%sf(0, k, l)
                        end if
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    if (bc_x%isothermal_in) then
                        do j = 1, buff_size
                            q_T_sf%sf(-j, k, l) = 2._wp*bc_x%Twall_in - q_T_sf%sf(j - 1, k, l)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(-j, k, l) = q_T_sf%sf(0, k, l)
                        end do
                    end if
                end if
            else  !< bc_x%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg) then
                            q_prim_vf(i)%sf(m + j, k, l) = -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve1
                        else
                            q_prim_vf(i)%sf(m + j, k, l) = q_prim_vf(i)%sf(m, k, l)
                        end if
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    if (bc_x%isothermal_out) then
                        do j = 1, buff_size
                            q_T_sf%sf(m + j, k, l) = 2._wp*bc_x%Twall_out - q_T_sf%sf(m - (j - 1), k, l)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(m + j, k, l) = q_T_sf%sf(m, k, l)
                        end do
                    end if
                end if
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg + 1) then
                            q_prim_vf(i)%sf(k, -j, l) = -q_prim_vf(i)%sf(k, j - 1, l) + 2._wp*bc_y%vb2
                        else
                            q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, 0, l)
                        end if
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    if (bc_y%isothermal_in) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, -j, l) = 2._wp*bc_y%Twall_in - q_T_sf%sf(k, j - 1, l)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(k, -j, l) = q_T_sf%sf(k, 0, l)
                        end do
                    end if
                end if
            else  !< bc_y%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg + 1) then
                            q_prim_vf(i)%sf(k, n + j, l) = -q_prim_vf(i)%sf(k, n - (j - 1), l) + 2._wp*bc_y%ve2
                        else
                            q_prim_vf(i)%sf(k, n + j, l) = q_prim_vf(i)%sf(k, n, l)
                        end if
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    if (bc_y%isothermal_out) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, n + j, l) = 2._wp*bc_y%Twall_out - q_T_sf%sf(k, n - (j - 1), l)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(k, n + j, l) = q_T_sf%sf(k, n, l)
                        end do
                    end if
                end if
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%end) then
                            q_prim_vf(i)%sf(k, l, -j) = -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb3
                        else
                            q_prim_vf(i)%sf(k, l, -j) = q_prim_vf(i)%sf(k, l, 0)
                        end if
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    if (bc_z%isothermal_in) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, -j) = 2._wp*bc_z%Twall_in - q_T_sf%sf(k, l, j - 1)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, -j) = q_T_sf%sf(k, l, 0)
                        end do
                    end if
                end if
            else  !< bc_z%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%end) then
                            q_prim_vf(i)%sf(k, l, p + j) = -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve3
                        else
                            q_prim_vf(i)%sf(k, l, p + j) = q_prim_vf(i)%sf(k, l, p)
                        end if
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    if (bc_z%isothermal_out) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, p + j) = 2._wp*bc_z%Twall_out - q_T_sf%sf(k, l, p - (j - 1))
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, p + j) = q_T_sf%sf(k, l, p)
                        end do
                    end if
                end if
            end if
        end if

    end subroutine s_slip_wall

    !> Apply no-slip wall boundary conditions by reflecting and negating all velocity components at the wall.
    subroutine s_no_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)

        $:GPU_ROUTINE(function_name='s_no_slip_wall',parallelism='[seq]', cray_inline=True)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in)                                    :: bc_dir, bc_loc
        integer, intent(in)                                    :: k, l
        integer                                                :: j, i
        type(scalar_field), optional, intent(inout)            :: q_T_sf

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  !< bc_x%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg) then
                            q_prim_vf(i)%sf(-j, k, l) = -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb1
                        else if (i == eqn_idx%mom%beg + 1 .and. num_dims > 1) then
                            q_prim_vf(i)%sf(-j, k, l) = -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb2
                        else if (i == eqn_idx%mom%beg + 2 .and. num_dims > 2) then
                            q_prim_vf(i)%sf(-j, k, l) = -q_prim_vf(i)%sf(j - 1, k, l) + 2._wp*bc_x%vb3
                        else
                            q_prim_vf(i)%sf(-j, k, l) = q_prim_vf(i)%sf(0, k, l)
                        end if
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    if (bc_x%isothermal_in) then
                        do j = 1, buff_size
                            q_T_sf%sf(-j, k, l) = 2._wp*bc_x%Twall_in - q_T_sf%sf(j - 1, k, l)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(-j, k, l) = q_T_sf%sf(0, k, l)
                        end do
                    end if
                end if
            else  !< bc_x%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg) then
                            q_prim_vf(i)%sf(m + j, k, l) = -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve1
                        else if (i == eqn_idx%mom%beg + 1 .and. num_dims > 1) then
                            q_prim_vf(i)%sf(m + j, k, l) = -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve2
                        else if (i == eqn_idx%mom%beg + 2 .and. num_dims > 2) then
                            q_prim_vf(i)%sf(m + j, k, l) = -q_prim_vf(i)%sf(m - (j - 1), k, l) + 2._wp*bc_x%ve3
                        else
                            q_prim_vf(i)%sf(m + j, k, l) = q_prim_vf(i)%sf(m, k, l)
                        end if
                    end do
                end do

                if (chemistry .and. present(q_T_sf)) then
                    if (bc_x%isothermal_out) then
                        do j = 1, buff_size
                            q_T_sf%sf(m + j, k, l) = 2._wp*bc_x%Twall_out - q_T_sf%sf(m - (j - 1), k, l)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(m + j, k, l) = q_T_sf%sf(m, k, l)
                        end do
                    end if
                end if
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg) then
                            q_prim_vf(i)%sf(k, -j, l) = -q_prim_vf(i)%sf(k, j - 1, l) + 2._wp*bc_y%vb1
                        else if (i == eqn_idx%mom%beg + 1 .and. num_dims > 1) then
                            q_prim_vf(i)%sf(k, -j, l) = -q_prim_vf(i)%sf(k, j - 1, l) + 2._wp*bc_y%vb2
                        else if (i == eqn_idx%mom%beg + 2 .and. num_dims > 2) then
                            q_prim_vf(i)%sf(k, -j, l) = -q_prim_vf(i)%sf(k, j - 1, l) + 2._wp*bc_y%vb3
                        else
                            q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, 0, l)
                        end if
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    if (bc_y%isothermal_in) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, -j, l) = 2._wp*bc_y%Twall_in - q_T_sf%sf(k, j - 1, l)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(k, -j, l) = q_T_sf%sf(k, 0, l)
                        end do
                    end if
                end if
            else  !< bc_y%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg) then
                            q_prim_vf(i)%sf(k, n + j, l) = -q_prim_vf(i)%sf(k, n - (j - 1), l) + 2._wp*bc_y%ve1
                        else if (i == eqn_idx%mom%beg + 1 .and. num_dims > 1) then
                            q_prim_vf(i)%sf(k, n + j, l) = -q_prim_vf(i)%sf(k, n - (j - 1), l) + 2._wp*bc_y%ve2
                        else if (i == eqn_idx%mom%beg + 2 .and. num_dims > 2) then
                            q_prim_vf(i)%sf(k, n + j, l) = -q_prim_vf(i)%sf(k, n - (j - 1), l) + 2._wp*bc_y%ve3
                        else
                            q_prim_vf(i)%sf(k, n + j, l) = q_prim_vf(i)%sf(k, n, l)
                        end if
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    if (bc_y%isothermal_out) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, n + j, l) = 2._wp*bc_y%Twall_out - q_T_sf%sf(k, n - (j - 1), l)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(k, n + j, l) = q_T_sf%sf(k, n, l)
                        end do
                    end if
                end if
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg) then
                            q_prim_vf(i)%sf(k, l, -j) = -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb1
                        else if (i == eqn_idx%mom%beg + 1 .and. num_dims > 1) then
                            q_prim_vf(i)%sf(k, l, -j) = -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb2
                        else if (i == eqn_idx%mom%beg + 2 .and. num_dims > 2) then
                            q_prim_vf(i)%sf(k, l, -j) = -q_prim_vf(i)%sf(k, l, j - 1) + 2._wp*bc_z%vb3
                        else
                            q_prim_vf(i)%sf(k, l, -j) = q_prim_vf(i)%sf(k, l, 0)
                        end if
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    if (bc_z%isothermal_in) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, -j) = 2._wp*bc_z%Twall_in - q_T_sf%sf(k, l, j - 1)
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, -j) = q_T_sf%sf(k, l, 0)
                        end do
                    end if
                end if
            else  !< bc_z%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        if (i == eqn_idx%mom%beg) then
                            q_prim_vf(i)%sf(k, l, p + j) = -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve1
                        else if (i == eqn_idx%mom%beg + 1 .and. num_dims > 1) then
                            q_prim_vf(i)%sf(k, l, p + j) = -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve2
                        else if (i == eqn_idx%mom%beg + 2 .and. num_dims > 2) then
                            q_prim_vf(i)%sf(k, l, p + j) = -q_prim_vf(i)%sf(k, l, p - (j - 1)) + 2._wp*bc_z%ve3
                        else
                            q_prim_vf(i)%sf(k, l, p + j) = q_prim_vf(i)%sf(k, l, p)
                        end if
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    if (bc_z%isothermal_out) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, p + j) = 2._wp*bc_z%Twall_out - q_T_sf%sf(k, l, p - (j - 1))
                        end do
                    else
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, p + j) = q_T_sf%sf(k, l, p)
                        end do
                    end if
                end if
            end if
        end if

    end subroutine s_no_slip_wall

    !> Apply Dirichlet boundary conditions by prescribing ghost cell values from stored boundary buffers.
    subroutine s_dirichlet(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)

        $:GPU_ROUTINE(function_name='s_dirichlet',parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in)                                    :: bc_dir, bc_loc
        integer, intent(in)                                    :: k, l
        integer                                                :: j, i
        type(scalar_field), optional, intent(inout)            :: q_T_sf

#ifdef MFC_SIMULATION
        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(-j, k, l) = bc_buffers(1, 1)%sf(i, k, l)
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(-j, k, l) = bc_buffers(1, 1)%sf(sys_size + 1, k, l)
                    end do
                end if
            else  !< bc_x%end
                do i = 1, sys_size
                    do j = 1, buff_size
                        q_prim_vf(i)%sf(m + j, k, l) = bc_buffers(1, 2)%sf(i, k, l)
                    end do
                end do
                if (chemistry .and. present(q_T_sf)) then
                    do j = 1, buff_size
                        q_T_sf%sf(m + j, k, l) = bc_buffers(1, 2)%sf(sys_size + 1, k, l)
                    end do
                end if
            end if
        else if (bc_dir == 2) then  !< y-direction
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                if (bc_loc == -1) then  !< bc_y%beg
                    do i = 1, sys_size
                        do j = 1, buff_size
                            q_prim_vf(i)%sf(k, -j, l) = bc_buffers(2, 1)%sf(k, i, l)
                        end do
                    end do
                    if (chemistry .and. present(q_T_sf)) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, -j, l) = bc_buffers(2, 1)%sf(k, sys_size + 1, l)
                        end do
                    end if
                else  !< bc_y%end
                    do i = 1, sys_size
                        do j = 1, buff_size
                            q_prim_vf(i)%sf(k, n + j, l) = bc_buffers(2, 2)%sf(k, i, l)
                        end do
                    end do
                    if (chemistry .and. present(q_T_sf)) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, n + j, l) = bc_buffers(2, 2)%sf(k, sys_size + 1, l)
                        end do
                    end if
                end if
            #:endif
        else if (bc_dir == 3) then  !< z-direction
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                if (bc_loc == -1) then  !< bc_z%beg
                    do i = 1, sys_size
                        do j = 1, buff_size
                            q_prim_vf(i)%sf(k, l, -j) = bc_buffers(3, 1)%sf(k, l, i)
                        end do
                    end do
                    if (chemistry .and. present(q_T_sf)) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, -j) = bc_buffers(3, 1)%sf(k, l, sys_size + 1)
                        end do
                    end if
                else  !< bc_z%end
                    do i = 1, sys_size
                        do j = 1, buff_size
                            q_prim_vf(i)%sf(k, l, p + j) = bc_buffers(3, 2)%sf(k, l, i)
                        end do
                    end do
                    if (chemistry .and. present(q_T_sf)) then
                        do j = 1, buff_size
                            q_T_sf%sf(k, l, p + j) = bc_buffers(3, 2)%sf(k, l, sys_size + 1)
                        end do
                    end if
                end if
            #:endif
        end if
#else
        call s_ghost_cell_extrapolation(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)
#endif

    end subroutine s_dirichlet

    !> Extrapolate QBMM bubble pressure and mass-vapor variables into ghost cells by copying boundary values.
    subroutine s_qbmm_extrapolation(bc_dir, bc_loc, k, l, pb_in, mv_in)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        integer, intent(in)                                                                                  :: bc_dir, bc_loc
        integer, intent(in)                                                                                  :: k, l
        integer                                                                                              :: j, q, i

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb_in(-j, k, l, q, i) = pb_in(0, k, l, q, i)
                            mv_in(-j, k, l, q, i) = mv_in(0, k, l, q, i)
                        end do
                    end do
                end do
            else  !< bc_x%end
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb_in(m + j, k, l, q, i) = pb_in(m, k, l, q, i)
                            mv_in(m + j, k, l, q, i) = mv_in(m, k, l, q, i)
                        end do
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb_in(k, -j, l, q, i) = pb_in(k, 0, l, q, i)
                            mv_in(k, -j, l, q, i) = mv_in(k, 0, l, q, i)
                        end do
                    end do
                end do
            else  !< bc_y%end
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb_in(k, n + j, l, q, i) = pb_in(k, n, l, q, i)
                            mv_in(k, n + j, l, q, i) = mv_in(k, n, l, q, i)
                        end do
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb_in(k, l, -j, q, i) = pb_in(k, l, 0, q, i)
                            mv_in(k, l, -j, q, i) = mv_in(k, l, 0, q, i)
                        end do
                    end do
                end do
            else  !< bc_z%end
                do i = 1, nb
                    do q = 1, nnode
                        do j = 1, buff_size
                            pb_in(k, l, p + j, q, i) = pb_in(k, l, p, q, i)
                            mv_in(k, l, p + j, q, i) = mv_in(k, l, p, q, i)
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_qbmm_extrapolation

    !> Apply periodic boundary conditions to the color function and its divergence fields.
    subroutine s_color_function_periodic(c_divs, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_color_function_periodic', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer, intent(in)                                        :: bc_dir, bc_loc
        integer, intent(in)                                        :: k, l
        integer                                                    :: j, i

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(-j, k, l) = c_divs(i)%sf(m - (j - 1), k, l)
                    end do
                end do
            else  !< bc_x%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(m + j, k, l) = c_divs(i)%sf(j - 1, k, l)
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, -j, l) = c_divs(i)%sf(k, n - (j - 1), l)
                    end do
                end do
            else  !< bc_y%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, n + j, l) = c_divs(i)%sf(k, j - 1, l)
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, -j) = c_divs(i)%sf(k, l, p - (j - 1))
                    end do
                end do
            else  !< bc_z%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, p + j) = c_divs(i)%sf(k, l, j - 1)
                    end do
                end do
            end if
        end if

    end subroutine s_color_function_periodic

    !> Apply reflective boundary conditions to the color function and its divergence fields.
    subroutine s_color_function_reflective(c_divs, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_color_function_reflective', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer, intent(in)                                        :: bc_dir, bc_loc
        integer, intent(in)                                        :: k, l
        integer                                                    :: j, i

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(-j, k, l) = -c_divs(i)%sf(j - 1, k, l)
                        else
                            c_divs(i)%sf(-j, k, l) = c_divs(i)%sf(j - 1, k, l)
                        end if
                    end do
                end do
            else  !< bc_x%end
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
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, -j, l) = -c_divs(i)%sf(k, j - 1, l)
                        else
                            c_divs(i)%sf(k, -j, l) = c_divs(i)%sf(k, j - 1, l)
                        end if
                    end do
                end do
            else  !< bc_y%end
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
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, l, -j) = -c_divs(i)%sf(k, l, j - 1)
                        else
                            c_divs(i)%sf(k, l, -j) = c_divs(i)%sf(k, l, j - 1)
                        end if
                    end do
                end do
            else  !< bc_z%end
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

    !> Extrapolate the color function and its divergence into ghost cells by copying boundary values.
    subroutine s_color_function_ghost_cell_extrapolation(c_divs, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_color_function_ghost_cell_extrapolation', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer, intent(in)                                        :: bc_dir, bc_loc
        integer, intent(in)                                        :: k, l
        integer                                                    :: j, i

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(-j, k, l) = c_divs(i)%sf(0, k, l)
                    end do
                end do
            else  !< bc_x%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(m + j, k, l) = c_divs(i)%sf(m, k, l)
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, -j, l) = c_divs(i)%sf(k, 0, l)
                    end do
                end do
            else  !< bc_y%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, n + j, l) = c_divs(i)%sf(k, n, l)
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, -j) = c_divs(i)%sf(k, l, 0)
                    end do
                end do
            else  !< bc_z%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, p + j) = c_divs(i)%sf(k, l, p)
                    end do
                end do
            end if
        end if

    end subroutine s_color_function_ghost_cell_extrapolation

    !> Apply periodic boundary conditions to the IGR Jacobian field by copying values from the opposite domain boundary.
    subroutine s_F_igr_periodic(jac_sf, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_F_igr_periodic', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(1:), intent(inout) :: jac_sf
        integer, intent(in)                              :: bc_dir, bc_loc
        integer, intent(in)                              :: k, l
        integer                                          :: j

        if (bc_dir == 1) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(-j, k, l) = jac_sf(1)%sf(m - j + 1, k, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(m + j, k, l) = jac_sf(1)%sf(j - 1, k, l)
                end do
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, -j, l) = jac_sf(1)%sf(k, n - j + 1, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, n + j, l) = jac_sf(1)%sf(k, j - 1, l)
                end do
            end if
        else
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, -j) = jac_sf(1)%sf(k, l, p - j + 1)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, p + j) = jac_sf(1)%sf(k, l, j - 1)
                end do
            end if
        end if

    end subroutine s_F_igr_periodic

    !> Apply reflective boundary conditions to the IGR Jacobian field by mirroring values across the boundary.
    subroutine s_F_igr_reflective(jac_sf, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_F_igr_reflective', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(1:), intent(inout) :: jac_sf
        integer, intent(in)                              :: bc_dir, bc_loc
        integer, intent(in)                              :: k, l
        integer                                          :: j

        if (bc_dir == 1) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(-j, k, l) = jac_sf(1)%sf(j - 1, k, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(m + j, k, l) = jac_sf(1)%sf(m - (j - 1), k, l)
                end do
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, -j, l) = jac_sf(1)%sf(k, j - 1, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, n + j, l) = jac_sf(1)%sf(k, n - (j - 1), l)
                end do
            end if
        else
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, -j) = jac_sf(1)%sf(k, l, j - 1)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, p + j) = jac_sf(1)%sf(k, l, p - (j - 1))
                end do
            end if
        end if

    end subroutine s_F_igr_reflective

    !> Extrapolate the IGR Jacobian field into ghost cells by copying boundary values.
    subroutine s_F_igr_ghost_cell_extrapolation(jac_sf, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_F_igr_ghost_cell_extrapolation', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(1:), intent(inout) :: jac_sf
        integer, intent(in)                              :: bc_dir, bc_loc
        integer, intent(in)                              :: k, l
        integer                                          :: j

        if (bc_dir == 1) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(-j, k, l) = jac_sf(1)%sf(0, k, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(m + j, k, l) = jac_sf(1)%sf(m, k, l)
                end do
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, -j, l) = jac_sf(1)%sf(k, 0, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, n + j, l) = jac_sf(1)%sf(k, n, l)
                end do
            end if
        else
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, -j) = jac_sf(1)%sf(k, l, 0)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, p + j) = jac_sf(1)%sf(k, l, p)
                end do
            end if
        end if

    end subroutine s_F_igr_ghost_cell_extrapolation

#ifndef MFC_PRE_PROCESS
    !> Apply periodic boundary conditions to grid variables by copying cell widths from opposite domain boundary.
    subroutine s_grid_periodic_bc(bc_dir, bc_loc, offset_dir)

        integer, intent(in)               :: bc_dir, bc_loc
        type(int_bounds_info), intent(in) :: offset_dir
        integer                           :: i

        if (bc_dir == 1) then
            if (bc_loc == -1) then
                do i = 1, buff_size
                    dx(-i) = dx(m - (i - 1))
                end do
                do i = 1, offset_dir%beg
                    x_cb(-1 - i) = x_cb(-i) - dx(-i)
                end do
                do i = 1, buff_size
                    x_cc(-i) = x_cc(1 - i) - (dx(1 - i) + dx(-i))/2._wp
                end do
            else
                do i = 1, buff_size
                    dx(m + i) = dx(i - 1)
                end do
                do i = 1, offset_dir%end
                    x_cb(m + i) = x_cb(m + (i - 1)) + dx(m + i)
                end do
                do i = 1, buff_size
                    x_cc(m + i) = x_cc(m + (i - 1)) + (dx(m + (i - 1)) + dx(m + i))/2._wp
                end do
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                do i = 1, buff_size
                    dy(-i) = dy(n - (i - 1))
                end do
                do i = 1, offset_dir%beg
                    y_cb(-1 - i) = y_cb(-i) - dy(-i)
                end do
                do i = 1, buff_size
                    y_cc(-i) = y_cc(1 - i) - (dy(1 - i) + dy(-i))/2._wp
                end do
            else
                do i = 1, buff_size
                    dy(n + i) = dy(i - 1)
                end do
                do i = 1, offset_dir%end
                    y_cb(n + i) = y_cb(n + (i - 1)) + dy(n + i)
                end do
                do i = 1, buff_size
                    y_cc(n + i) = y_cc(n + (i - 1)) + (dy(n + (i - 1)) + dy(n + i))/2._wp
                end do
            end if
        else
            if (bc_loc == -1) then
                do i = 1, buff_size
                    dz(-i) = dz(p - (i - 1))
                end do
                do i = 1, offset_dir%beg
                    z_cb(-1 - i) = z_cb(-i) - dz(-i)
                end do
                do i = 1, buff_size
                    z_cc(-i) = z_cc(1 - i) - (dz(1 - i) + dz(-i))/2._wp
                end do
            else
                do i = 1, buff_size
                    dz(p + i) = dz(i - 1)
                end do
                do i = 1, offset_dir%end
                    z_cb(p + i) = z_cb(p + (i - 1)) + dz(p + i)
                end do
                do i = 1, buff_size
                    z_cc(p + i) = z_cc(p + (i - 1)) + (dz(p + (i - 1)) + dz(p + i))/2._wp
                end do
            end if
        end if

    end subroutine s_grid_periodic_bc

    !> Apply reflective boundary conditions to grid variables by mirroring cell widths across the boundary.
    subroutine s_grid_reflective_bc(bc_dir, bc_loc, offset_dir)

        integer, intent(in)               :: bc_dir, bc_loc
        type(int_bounds_info), intent(in) :: offset_dir
        integer                           :: i

        if (bc_dir == 1) then
            if (bc_loc == -1) then
                do i = 1, buff_size
                    dx(-i) = dx(i - 1)
                end do
                do i = 1, offset_dir%beg
                    x_cb(-1 - i) = x_cb(-i) - dx(-i)
                end do
                do i = 1, buff_size
                    x_cc(-i) = x_cc(1 - i) - (dx(1 - i) + dx(-i))/2._wp
                end do
            else
                do i = 1, buff_size
                    dx(m + i) = dx(m - (i - 1))
                end do
                do i = 1, offset_dir%end
                    x_cb(m + i) = x_cb(m + (i - 1)) + dx(m + i)
                end do
                do i = 1, buff_size
                    x_cc(m + i) = x_cc(m + (i - 1)) + (dx(m + (i - 1)) + dx(m + i))/2._wp
                end do
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                do i = 1, buff_size
                    dy(-i) = dy(i - 1)
                end do
                do i = 1, offset_dir%beg
                    y_cb(-1 - i) = y_cb(-i) - dy(-i)
                end do
                do i = 1, buff_size
                    y_cc(-i) = y_cc(1 - i) - (dy(1 - i) + dy(-i))/2._wp
                end do
            else
                do i = 1, buff_size
                    dy(n + i) = dy(n - (i - 1))
                end do
                do i = 1, offset_dir%end
                    y_cb(n + i) = y_cb(n + (i - 1)) + dy(n + i)
                end do
                do i = 1, buff_size
                    y_cc(n + i) = y_cc(n + (i - 1)) + (dy(n + (i - 1)) + dy(n + i))/2._wp
                end do
            end if
        else
            if (bc_loc == -1) then
                do i = 1, buff_size
                    dz(-i) = dz(i - 1)
                end do
                do i = 1, offset_dir%beg
                    z_cb(-1 - i) = z_cb(-i) - dz(-i)
                end do
                do i = 1, buff_size
                    z_cc(-i) = z_cc(1 - i) - (dz(1 - i) + dz(-i))/2._wp
                end do
            else
                do i = 1, buff_size
                    dz(p + i) = dz(p - (i - 1))
                end do
                do i = 1, offset_dir%end
                    z_cb(p + i) = z_cb(p + (i - 1)) + dz(p + i)
                end do
                do i = 1, buff_size
                    z_cc(p + i) = z_cc(p + (i - 1)) + (dz(p + (i - 1)) + dz(p + i))/2._wp
                end do
            end if
        end if

    end subroutine s_grid_reflective_bc

    !> Extrapolate grid variables by copying boundary cell width into ghost cells.
    subroutine s_grid_ghost_cell_extrapolation_bc(bc_dir, bc_loc, offset_dir)

        integer, intent(in)               :: bc_dir, bc_loc
        type(int_bounds_info), intent(in) :: offset_dir
        integer                           :: i

        if (bc_dir == 1) then
            if (bc_loc == -1) then
                do i = 1, buff_size
                    dx(-i) = dx(0)
                end do
                do i = 1, offset_dir%beg
                    x_cb(-1 - i) = x_cb(-i) - dx(-i)
                end do
                do i = 1, buff_size
                    x_cc(-i) = x_cc(1 - i) - (dx(1 - i) + dx(-i))/2._wp
                end do
            else
                do i = 1, buff_size
                    dx(m + i) = dx(m)
                end do
                do i = 1, offset_dir%end
                    x_cb(m + i) = x_cb(m + (i - 1)) + dx(m + i)
                end do
                do i = 1, buff_size
                    x_cc(m + i) = x_cc(m + (i - 1)) + (dx(m + (i - 1)) + dx(m + i))/2._wp
                end do
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                do i = 1, buff_size
                    dy(-i) = dy(0)
                end do
                do i = 1, offset_dir%beg
                    y_cb(-1 - i) = y_cb(-i) - dy(-i)
                end do
                do i = 1, buff_size
                    y_cc(-i) = y_cc(1 - i) - (dy(1 - i) + dy(-i))/2._wp
                end do
            else
                do i = 1, buff_size
                    dy(n + i) = dy(n)
                end do
                do i = 1, offset_dir%end
                    y_cb(n + i) = y_cb(n + (i - 1)) + dy(n + i)
                end do
                do i = 1, buff_size
                    y_cc(n + i) = y_cc(n + (i - 1)) + (dy(n + (i - 1)) + dy(n + i))/2._wp
                end do
            end if
        else
            if (bc_loc == -1) then
                do i = 1, buff_size
                    dz(-i) = dz(0)
                end do
                do i = 1, offset_dir%beg
                    z_cb(-1 - i) = z_cb(-i) - dz(-i)
                end do
                do i = 1, buff_size
                    z_cc(-i) = z_cc(1 - i) - (dz(1 - i) + dz(-i))/2._wp
                end do
            else
                do i = 1, buff_size
                    dz(p + i) = dz(p)
                end do
                do i = 1, offset_dir%end
                    z_cb(p + i) = z_cb(p + (i - 1)) + dz(p + i)
                end do
                do i = 1, buff_size
                    z_cc(p + i) = z_cc(p + (i - 1)) + (dz(p + (i - 1)) + dz(p + i))/2._wp
                end do
            end if
        end if

    end subroutine s_grid_ghost_cell_extrapolation_bc

    !> Apply axis boundary conditions to grid variables for cylindrical coordinates.
    subroutine s_grid_axis_bc(bc_loc, offset_dir)

        integer, intent(in)               :: bc_loc
        type(int_bounds_info), intent(in) :: offset_dir
        integer                           :: i

        if (bc_loc == -1) then
            do i = 1, buff_size
                dy(-i) = dy(i - 1)
            end do
            do i = 1, offset_dir%beg
                y_cb(-1 - i) = y_cb(-i) - dy(-i)
            end do
            do i = 1, buff_size
                y_cc(-i) = y_cc(1 - i) - (dy(1 - i) + dy(-i))/2._wp
            end do
        end if

    end subroutine s_grid_axis_bc
#endif
end module m_boundary_primitives
