!>
!! @file m_boundary_conditions.fpp
!! @brief Contains module m_boundary_conditions

!> @brief This module contains
module m_boundary_conditions

    use m_derived_types

    use m_global_parameters
#ifdef MFC_MPI
    use mpi
#endif
    use m_delay_file_access

    use m_compile_specific

    use m_boundary_common

    implicit none

    real(wp) :: x_centroid, y_centroid, z_centroid
    real(wp) :: length_x, length_y, length_z
    real(wp) :: radius
    type(bounds_info) :: x_boundary, y_boundary, z_boundary  !<

    private; public :: s_apply_boundary_patches

contains
    impure subroutine s_line_segment_bc(patch_id, bc_type)

        type(integer_field), dimension(1:num_dims, 1:2), intent(inout) :: bc_type
        integer, intent(in) :: patch_id

        integer :: j

        ! Patch is a vertical line at x_beg or x_end
        if (patch_bc(patch_id)%dir == 1) then
            y_centroid = patch_bc(patch_id)%centroid(2)
            length_y = patch_bc(patch_id)%length(2)

            y_boundary%beg = y_centroid - 0.5_wp*length_y
            y_boundary%end = y_centroid + 0.5_wp*length_y

            ! Patch is a vertical line at x_beg and x_beg is a domain boundary
            #:for BOUND, X, LOC, IDX in [('beg', '-i', -1, 1), ('end', 'm+i', 1, 2)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_x%${BOUND}$ < 0) then
                    do j = 0, n
                        if (y_cc(j) > y_boundary%beg .and. y_cc(j) < y_boundary%end) then
                            bc_type(1, ${IDX}$)%sf(0, j, 0) = patch_bc(patch_id)%type
                        end if
                    end do
                end if
            #:endfor
        end if

        ! Patch is a vertical line at y_beg or y_end
        if (patch_bc(patch_id)%dir == 2) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            length_x = patch_bc(patch_id)%length(1)

            x_boundary%beg = x_centroid - 0.5_wp*length_x
            x_boundary%end = x_centroid + 0.5_wp*length_x

            ! Patch is a vertical line at x_beg and x_beg is a domain boundary
            #:for BOUND, Y, LOC, IDX in [('beg', '-i', -1, 1), ('end', 'n+i', 1, 2)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_y%${BOUND}$ < 0) then
                    do j = 0, m
                        if (x_cc(j) > x_boundary%beg .and. x_cc(j) < x_boundary%end) then
                            bc_type(2, ${IDX}$)%sf(j, 0, 0) = patch_bc(patch_id)%type
                        end if
                    end do
                end if
            #:endfor
        end if

    end subroutine s_line_segment_bc

    impure subroutine s_circle_bc(patch_id, bc_type)

        type(integer_field), dimension(1:num_dims, 1:2), intent(inout) :: bc_type

        integer, intent(in) :: patch_id

        integer :: j, k
        if (patch_bc(patch_id)%dir == 1) then
            y_centroid = patch_bc(patch_id)%centroid(2)
            z_centroid = patch_bc(patch_id)%centroid(3)
            radius = patch_bc(patch_id)%radius
            ! Patch is a circle at x_beg and x_beg is a domain boundary
            #:for BOUND, X, LOC, IDX in [('beg', '-i', -1, 1), ('end', 'm+i', 1, 2)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_x%${BOUND}$ < 0) then
                    do k = 0, p
                        do j = 0, n
                            if ((z_cc(k) - z_centroid)**2._wp + &
                                (y_cc(j) - y_centroid)**2._wp <= radius**2._wp) then
                                bc_type(1, ${IDX}$)%sf(0, j, k) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if
        if (patch_bc(patch_id)%dir == 2) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            z_centroid = patch_bc(patch_id)%centroid(3)
            radius = patch_bc(patch_id)%radius
            ! Patch is a circle at y_beg and y_beg is a domain boundary
            #:for BOUND, Y, LOC, IDX in [('beg', '-i', -1, 1), ('end', 'n+i', 1, 2)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_y%${BOUND}$ < 0) then
                    do k = 0, p
                        do j = 0, m
                            if ((z_cc(k) - z_centroid)**2._wp + &
                                (x_cc(j) - x_centroid)**2._wp <= radius**2._wp) then
                                bc_type(2, ${IDX}$)%sf(j, 0, k) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if
        if (patch_bc(patch_id)%dir == 3) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            y_centroid = patch_bc(patch_id)%centroid(2)
            radius = patch_bc(patch_id)%radius
            #:for BOUND, Z, LOC, IDX in [('beg', '-i', -1, 1), ('end', 'p+i', 1, 2)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_z%${BOUND}$ < 0) then
                    do k = 0, n
                        do j = 0, m
                            if ((y_cc(k) - y_centroid)**2._wp + &
                                (x_cc(j) - x_centroid)**2._wp <= radius**2._wp) then
                                bc_type(3, ${IDX}$)%sf(j, k, 0) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if

    end subroutine s_circle_bc

    impure subroutine s_rectangle_bc(patch_id, bc_type)

        type(integer_field), dimension(1:num_dims, 1:2), intent(inout) :: bc_type

        integer, intent(in) :: patch_id
        integer :: j, k
        if (patch_bc(patch_id)%dir == 1) then
            y_centroid = patch_bc(patch_id)%centroid(2)
            z_centroid = patch_bc(patch_id)%centroid(3)
            length_y = patch_bc(patch_id)%length(2)
            length_z = patch_bc(patch_id)%length(3)

            y_boundary%beg = y_centroid - 0.5_wp*length_y
            y_boundary%end = y_centroid + 0.5_wp*length_y

            z_boundary%beg = z_centroid - 0.5_wp*length_z
            z_boundary%end = z_centroid + 0.5_wp*length_z
            ! Patch is a circle at x_beg and x_beg is a domain boundary
            #:for BOUND, X, LOC, IDX in [('beg', '-i', -1, 1), ('end', 'm+i', 1, 2)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_x%${BOUND}$ < 0) then
                    do k = 0, p
                        do j = 0, n
                            if (y_boundary%beg <= y_cc(j) .and. &
                                y_boundary%end >= y_cc(j) .and. &
                                z_boundary%beg <= z_cc(k) .and. &
                                z_boundary%end >= z_cc(k)) then
                                bc_type(1, ${IDX}$)%sf(0, j, k) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if
        if (patch_bc(patch_id)%dir == 2) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            z_centroid = patch_bc(patch_id)%centroid(3)
            length_x = patch_bc(patch_id)%length(1)
            length_z = patch_bc(patch_id)%length(3)

            x_boundary%beg = x_centroid - 0.5_wp*length_x
            x_boundary%end = x_centroid + 0.5_wp*length_x

            z_boundary%beg = z_centroid - 0.5_wp*length_z
            z_boundary%end = z_centroid + 0.5_wp*length_z
            ! Patch is a circle at y_beg and y_beg is a domain boundary
            #:for BOUND, Y, LOC, IDX in [('beg', '-i', -1, 1), ('end', 'n+i', 1, 2)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_y%${BOUND}$ < 0) then
                    do k = 0, p
                        do j = 0, m
                            if (x_boundary%beg <= x_cc(j) .and. &
                                x_boundary%end >= x_cc(j) .and. &
                                z_boundary%beg <= z_cc(k) .and. &
                                z_boundary%end >= z_cc(k)) then
                                bc_type(2, ${IDX}$)%sf(j, 0, k) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if
        if (patch_bc(patch_id)%dir == 3) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            y_centroid = patch_bc(patch_id)%centroid(2)
            length_x = patch_bc(patch_id)%length(1)
            length_y = patch_bc(patch_id)%length(2)

            x_boundary%beg = x_centroid - 0.5_wp*length_x
            x_boundary%end = x_centroid + 0.5_wp*length_x

            y_boundary%beg = y_centroid - 0.5_wp*length_y
            y_boundary%end = y_centroid + 0.5_wp*length_y
            #:for BOUND, Z, LOC, IDX in [('beg', '-i', -1, 1), ('end', 'p+i', 1, 2)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_z%${BOUND}$ < 0) then
                    do k = 0, n
                        do j = 0, m
                            if (x_boundary%beg <= x_cc(j) .and. &
                                x_boundary%end >= x_cc(j) .and. &
                                y_boundary%beg <= y_cc(k) .and. &
                                y_boundary%end >= y_cc(k)) then
                                bc_type(3, ${IDX}$)%sf(j, k, 0) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if

    end subroutine s_rectangle_bc

    impure subroutine s_apply_boundary_patches(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, 1:2) :: bc_type
        integer :: i

        !< Apply 2D patches to 3D domain
        if (p > 0) then
            do i = 1, num_bc_patches
                if (proc_rank == 0) then
                    print *, 'Processing boundary condition patch', i
                end if

                if (patch_bc(i)%geometry == 2) then
                    call s_circle_bc(i, bc_type)
                elseif (patch_bc(i)%geometry == 3) then
                    call s_rectangle_bc(i, bc_type)
                end if
            end do
            !< Apply 1D patches to 2D domain
        elseif (n > 0) then
            do i = 1, num_bc_patches
                if (proc_rank == 0) then
                    print *, 'Processing boundary condition patch', i
                end if

                if (patch_bc(i)%geometry == 1) then
                    call s_line_segment_bc(i, bc_type)
                end if
            end do
        end if

    end subroutine s_apply_boundary_patches

end module m_boundary_conditions
