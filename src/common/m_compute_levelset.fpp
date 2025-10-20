!>
!! @file m_compute_levelset.fpp
!! @brief Contains module m_compute_levelset

#:include 'macros.fpp'

!> @brief This module is used to handle all operations related to immersed
!!              boundary methods (IBMs)
module m_compute_levelset

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    implicit none

    private; public :: s_cylinder_levelset, s_circle_levelset, &
 s_airfoil_levelset, &
 s_3D_airfoil_levelset, &
 s_rectangle_levelset, &
 s_cuboid_levelset, &
 s_sphere_levelset

contains

    pure subroutine s_circle_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(wp) :: radius, dist
        real(wp) :: x_centroid, y_centroid
        real(wp), dimension(3) :: dist_vec

        integer :: i, j !< Loop index variables

        radius = patch_ib(ib_patch_id)%radius
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid

        do i = 0, m
            do j = 0, n

                dist_vec(1) = x_cc(i) - x_centroid
                dist_vec(2) = y_cc(j) - y_centroid
                dist_vec(3) = 0
                dist = sqrt(sum(dist_vec**2))
                levelset%sf(i, j, 0, ib_patch_id) = dist - radius
                if (f_approx_equal(dist, 0._wp)) then
                    levelset_norm%sf(i, j, 0, ib_patch_id, :) = 0
                else
                    levelset_norm%sf(i, j, 0, ib_patch_id, :) = &
                        dist_vec(:)/dist
                end if

            end do
        end do

    end subroutine s_circle_levelset

    pure subroutine s_airfoil_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(inout), optional :: levelset
        type(levelset_norm_field), intent(inout), optional :: levelset_norm
        integer, intent(in) :: ib_patch_id

        real(wp) :: dist, global_dist
        integer :: global_id
        real(wp) :: x_centroid, y_centroid
        real(wp), dimension(3) :: dist_vec

        real(wp), dimension(1:3) :: xy_local !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        integer :: i, j, k !< Loop index variables

        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)

        do i = 0, m
            do j = 0, n
                xy_local = [x_cc(i) - x_centroid, y_cc(j) - y_centroid, 0._wp] ! get coordinate frame centered on IB
                xy_local = matmul(inverse_rotation, xy_local) ! rotate the frame into the IB's coordinate

                if (xy_local(2) >= 0._wp) then
                    ! finds the location on the airfoil grid with the minimum distance (closest)
                    do k = 1, Np
                        dist_vec(1) = xy_local(1) - airfoil_grid_u(k)%x
                        dist_vec(2) = xy_local(2) - airfoil_grid_u(k)%y
                        dist_vec(3) = 0._wp
                        dist = sqrt(sum(dist_vec**2))
                        if (k == 1) then
                            global_dist = dist
                            global_id = k
                        else
                            if (dist < global_dist) then
                                global_dist = dist
                                global_id = k
                            end if
                        end if
                    end do
                    dist_vec(1) = xy_local(1) - airfoil_grid_u(global_id)%x
                    dist_vec(2) = xy_local(2) - airfoil_grid_u(global_id)%y
                    dist_vec(3) = 0
                    dist = global_dist
                else
                    ! TODO :: This looks identical to the above code but using the lower array. Refactor this.
                    do k = 1, Np
                        dist_vec(1) = xy_local(1) - airfoil_grid_l(k)%x
                        dist_vec(2) = xy_local(2) - airfoil_grid_l(k)%y
                        dist_vec(3) = 0
                        dist = sqrt(sum(dist_vec**2))
                        if (k == 1) then
                            global_dist = dist
                            global_id = k
                        else
                            if (dist < global_dist) then
                                global_dist = dist
                                global_id = k
                            end if
                        end if
                    end do
                    dist_vec(1) = xy_local(1) - airfoil_grid_l(global_id)%x
                    dist_vec(2) = xy_local(2) - airfoil_grid_l(global_id)%y
                    dist_vec(3) = 0
                    dist = global_dist
                end if

                levelset%sf(i, j, 0, ib_patch_id) = dist
                if (f_approx_equal(dist, 0._wp)) then
                    levelset_norm%sf(i, j, 0, ib_patch_id, :) = 0._wp
                else
                    levelset_norm%sf(i, j, 0, ib_patch_id, :) = &
                        matmul(rotation, dist_vec(:))/dist ! convert the normal vector back to global grid coordinates
                end if

            end do
        end do

    end subroutine s_airfoil_levelset

    pure subroutine s_3D_airfoil_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(wp) :: dist, dist_surf, dist_side, global_dist
        integer :: global_id
        real(wp) :: x_centroid, y_centroid, z_centroid, lz, z_max, z_min
        real(wp), dimension(3) :: dist_vec

        real(wp), dimension(1:3) :: xyz_local !< x, y, z coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        real(wp) :: length_z

        integer :: i, j, k, l !< Loop index variables

        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid
        lz = patch_ib(ib_patch_id)%length_z
        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)

        z_max = z_centroid + lz/2
        z_min = z_centroid - lz/2

        do l = 0, p
            do j = 0, n
                do i = 0, m

                    xyz_local = [x_cc(i) - x_centroid, y_cc(j) - y_centroid, z_cc(l) - z_centroid] ! get coordinate frame centered on IB
                    xyz_local = matmul(inverse_rotation, xyz_local) ! rotate the frame into the IB's coordinates

                    if (xyz_local(2) >= y_centroid) then
                        do k = 1, Np
                            dist_vec(1) = xyz_local(1) - airfoil_grid_u(k)%x
                            dist_vec(2) = xyz_local(2) - airfoil_grid_u(k)%y
                            dist_vec(3) = 0
                            dist_surf = sqrt(sum(dist_vec**2))
                            if (k == 1) then
                                global_dist = dist_surf
                                global_id = k
                            else
                                if (dist_surf < global_dist) then
                                    global_dist = dist_surf
                                    global_id = k
                                end if
                            end if
                        end do
                        dist_vec(1) = xyz_local(1) - airfoil_grid_u(global_id)%x
                        dist_vec(2) = xyz_local(2) - airfoil_grid_u(global_id)%y
                        dist_vec(3) = 0
                        dist_surf = global_dist
                    else
                        do k = 1, Np
                            dist_vec(1) = xyz_local(1) - airfoil_grid_l(k)%x
                            dist_vec(2) = xyz_local(2) - airfoil_grid_l(k)%y
                            dist_vec(3) = 0
                            dist_surf = sqrt(sum(dist_vec**2))
                            if (k == 1) then
                                global_dist = dist_surf
                                global_id = k
                            else
                                if (dist_surf < global_dist) then
                                    global_dist = dist_surf
                                    global_id = k
                                end if
                            end if
                        end do
                        dist_vec(1) = xyz_local(1) - airfoil_grid_l(global_id)%x
                        dist_vec(2) = xyz_local(2) - airfoil_grid_l(global_id)%y
                        dist_vec(3) = 0
                        dist_surf = global_dist
                    end if

                    dist_side = min(abs(z_cc(l) - z_min), abs(z_max - z_cc(l)))

                    if (dist_side < dist_surf) then
                        levelset%sf(i, j, l, ib_patch_id) = dist_side
                        if (f_approx_equal(dist_side, abs(z_cc(l) - z_min))) then
                            levelset_norm%sf(i, j, l, ib_patch_id, :) = (/0, 0, -1/)
                        else
                            levelset_norm%sf(i, j, l, ib_patch_id, :) = (/0, 0, 1/)
                        end if
                        levelset_norm%sf(i, j, l, ib_patch_id, :) = &
                            matmul(rotation, levelset_norm%sf(i, j, l, ib_patch_id, :)/dist_surf)
                    else
                        levelset%sf(i, j, l, ib_patch_id) = dist_surf
                        if (f_approx_equal(dist_surf, 0._wp)) then
                            levelset_norm%sf(i, j, l, ib_patch_id, :) = 0
                        else
                            levelset_norm%sf(i, j, l, ib_patch_id, :) = &
                                matmul(rotation, dist_vec(:)/dist_surf)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine s_3D_airfoil_levelset

    !>  Initialize IBM module
    pure subroutine s_rectangle_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm

        integer, intent(in) :: ib_patch_id
        real(wp) :: top_right(2), bottom_left(2)
        real(wp) :: min_dist
        real(wp) :: side_dists(4)

        real(wp) :: x_centroid, y_centroid
        real(wp) :: length_x, length_y
        real(wp), dimension(1:3) :: xy_local !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        integer :: i, j, k !< Loop index variables
        integer :: idx !< Shortest path direction indicator

        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)

        top_right(1) = length_x/2
        top_right(2) = length_y/2
        bottom_left(1) = -length_x/2
        bottom_left(2) = -length_y/2

        do i = 0, m
            do j = 0, n
                xy_local = [x_cc(i) - x_centroid, y_cc(j) - y_centroid, 0._wp]
                xy_local = matmul(inverse_rotation, xy_local)

                if ((xy_local(1) > bottom_left(1) .and. xy_local(1) < top_right(1)) .or. &
                    (xy_local(2) > bottom_left(2) .and. xy_local(2) < top_right(2))) then

                    side_dists(1) = bottom_left(1) - xy_local(1)
                    side_dists(2) = top_right(1) - xy_local(1)
                    side_dists(3) = bottom_left(2) - xy_local(2)
                    side_dists(4) = top_right(2) - xy_local(2)
                    min_dist = initial_distance_buffer
                    idx = 1

                    do k = 1, 4
                        if (abs(side_dists(k)) < abs(min_dist)) then
                            idx = k
                            min_dist = side_dists(idx)
                        end if
                    end do

                    levelset%sf(i, j, 0, ib_patch_id) = side_dists(idx)
                    levelset_norm%sf(i, j, 0, ib_patch_id, :) = 0._wp
                    if (.not. f_approx_equal(side_dists(idx), 0._wp)) then
                        if (idx == 1 .or. idx == 2) then
                            ! vector points along the x axis
                            levelset_norm%sf(i, j, 0, ib_patch_id, 1) = side_dists(idx)/ &
                                                                        abs(side_dists(idx))
                        else
                            ! vector points along the y axis
                            levelset_norm%sf(i, j, 0, ib_patch_id, 2) = side_dists(idx)/ &
                                                                        abs(side_dists(idx))
                        end if
                        ! convert the normal vector back into the global coordinate system
                        levelset_norm%sf(i, j, 0, ib_patch_id, :) = &
                            matmul(rotation, levelset_norm%sf(i, j, 0, ib_patch_id, :))
                    else
                        levelset_norm%sf(i, j, 0, ib_patch_id, :) = 0._wp
                    end if
                end if
            end do
        end do

    end subroutine s_rectangle_levelset

    pure subroutine s_cuboid_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm

        integer, intent(IN) :: ib_patch_id
        real(wp) :: Right, Left, Bottom, Top, Front, Back
        real(wp) :: x, y, z, min_dist
        real(wp) :: side_dists(6)

        real(wp) :: x_centroid, y_centroid, z_centroid
        real(wp) :: length_x, length_y, length_z
        real(wp), dimension(1:3) :: xyz_local !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        integer :: i, j, k, l, idx !< Loop index variables

        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        length_z = patch_ib(ib_patch_id)%length_z

        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid

        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)

        Right = length_x/2
        Left = -length_x/2

        Top = length_y/2
        Bottom = -length_y/2

        Front = length_z/2
        Back = -length_z/2

        do i = 0, m
            do j = 0, n
                do k = 0, p

                    xyz_local = [x_cc(i) - x_centroid, y_cc(j) - y_centroid, z_cc(k) - z_centroid] ! get coordinate frame centered on IB
                    xyz_local = matmul(inverse_rotation, xyz_local) ! rotate the frame into the IB's coordinate

                    if ((xyz_local(1) > Left .and. xyz_local(1) < Right) .or. &
                        (xyz_local(2) > Bottom .and. xyz_local(2) < Top) .or. &
                        (xyz_local(3) > Back .and. xyz_local(3) < Front)) then

                        side_dists(1) = Left - xyz_local(1)
                        side_dists(2) = xyz_local(1) - Right
                        side_dists(3) = Bottom - xyz_local(2)
                        side_dists(4) = xyz_local(2) - Top
                        side_dists(5) = Back - xyz_local(3)
                        side_dists(6) = xyz_local(3) - Front
                        min_dist = minval(abs(side_dists))

                        ! TODO :: The way that this is written, it looks like we will
                        ! trigger at the first size that is close to the minimum distance,
                        ! meaning corners where side_dists are the same will
                        ! trigger on what may not actually be the minimum,
                        ! leading to undesired behavior. This should be resolved
                        ! and this code should be cleaned up. I verified this behavior
                        ! with tests.
                        levelset_norm%sf(i, j, k, ib_patch_id, :) = 0._wp
                        if (f_approx_equal(min_dist, abs(side_dists(1)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(1)
                            if (.not. f_approx_equal(side_dists(1), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 1) = side_dists(1)/ &
                                                                            abs(side_dists(1))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(2)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(2)
                            if (.not. f_approx_equal(side_dists(2), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 1) = -side_dists(2)/ &
                                                                            abs(side_dists(2))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(3)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(3)
                            if (.not. f_approx_equal(side_dists(3), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 2) = side_dists(3)/ &
                                                                            abs(side_dists(3))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(4)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(4)
                            if (.not. f_approx_equal(side_dists(4), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 2) = -side_dists(4)/ &
                                                                            abs(side_dists(4))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(5)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(5)
                            if (.not. f_approx_equal(side_dists(5), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 3) = side_dists(5)/ &
                                                                            abs(side_dists(5))
                            end if

                        else if (f_approx_equal(min_dist, abs(side_dists(6)))) then
                            levelset%sf(i, j, k, ib_patch_id) = side_dists(6)
                            if (.not. f_approx_equal(side_dists(6), 0._wp)) then
                                levelset_norm%sf(i, j, k, ib_patch_id, 3) = -side_dists(6)/ &
                                                                            abs(side_dists(6))
                            end if
                        end if
                        levelset_norm%sf(i, j, 0, ib_patch_id, :) = &
                            matmul(rotation, levelset_norm%sf(i, j, 0, ib_patch_id, :))
                    end if
                end do
            end do
        end do

    end subroutine s_cuboid_levelset

    pure subroutine s_sphere_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(wp) :: radius, dist
        real(wp) :: x_centroid, y_centroid, z_centroid
        real(wp), dimension(3) :: dist_vec

        integer :: i, j, k !< Loop index variables

        radius = patch_ib(ib_patch_id)%radius
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid

        do i = 0, m
            do j = 0, n
                do k = 0, p
                    dist_vec(1) = x_cc(i) - x_centroid
                    dist_vec(2) = y_cc(j) - y_centroid
                    dist_vec(3) = z_cc(k) - z_centroid
                    dist = sqrt(sum(dist_vec**2))
                    levelset%sf(i, j, k, ib_patch_id) = dist - radius
                    if (f_approx_equal(dist, 0._wp)) then
                        levelset_norm%sf(i, j, k, ib_patch_id, :) = (/1, 0, 0/)
                    else
                        levelset_norm%sf(i, j, k, ib_patch_id, :) = &
                            dist_vec(:)/dist
                    end if
                end do
            end do
        end do

    end subroutine s_sphere_levelset

    pure subroutine s_cylinder_levelset(ib_patch_id, levelset, levelset_norm)

        type(levelset_field), intent(INOUT), optional :: levelset
        type(levelset_norm_field), intent(INOUT), optional :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(wp) :: radius
        real(wp) :: x_centroid, y_centroid, z_centroid
        real(wp) :: length_x, length_y, length_z
        real(wp), dimension(3) :: centroid_vec, dist_sides_vec, dist_surface_vec
        real(wp) :: dist_side, dist_surface, side_pos
        type(bounds_info) :: boundary
        integer :: i, j, k !< Loop index variables

        real(wp), dimension(1:3) :: xyz_local !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        radius = patch_ib(ib_patch_id)%radius
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid
        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        length_z = patch_ib(ib_patch_id)%length_z

        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)

        if (.not. f_approx_equal(length_x, 0._wp)) then
            boundary%beg = -0.5_wp*length_x
            boundary%end = 0.5_wp*length_x
            dist_sides_vec = (/1, 0, 0/)
            dist_surface_vec = (/0, 1, 1/)
        else if (.not. f_approx_equal(length_y, 0._wp)) then
            boundary%beg = -0.5_wp*length_y
            boundary%end = 0.5_wp*length_y
            dist_sides_vec = (/0, 1, 0/)
            dist_surface_vec = (/1, 0, 1/)
        else if (.not. f_approx_equal(length_z, 0._wp)) then
            boundary%beg = -0.5_wp*length_z
            boundary%end = 0.5_wp*length_z
            dist_sides_vec = (/0, 0, 1/)
            dist_surface_vec = (/1, 1, 0/)
        end if

        do i = 0, m
            do j = 0, n
                do k = 0, p
                    xyz_local = [x_cc(i) - x_centroid, y_cc(j) - y_centroid, z_cc(k) - z_centroid] ! get coordinate frame centered on IB
                    xyz_local = matmul(inverse_rotation, xyz_local) ! rotate the frame into the IB's coordinates

                    ! get distance to flat edge of cylinder
                    side_pos = dot_product(xyz_local, dist_sides_vec)
                    dist_side = min(abs(side_pos - boundary%beg), &
                                    abs(boundary%end - side_pos))
                    ! get distance to curved side of cylinder
                    dist_surface = norm2(xyz_local*dist_surface_vec) &
                                   - radius

                    if (dist_side < abs(dist_surface)) then
                        ! if the closest edge is flat
                        levelset%sf(i, j, k, ib_patch_id) = -dist_side
                        if (f_approx_equal(dist_side, abs(side_pos - boundary%beg))) then
                            levelset_norm%sf(i, j, k, ib_patch_id, :) = -dist_sides_vec
                        else
                            levelset_norm%sf(i, j, k, ib_patch_id, :) = dist_sides_vec
                        end if
                    else
                        levelset%sf(i, j, k, ib_patch_id) = dist_surface

                        levelset_norm%sf(i, j, k, ib_patch_id, :) = &
                            xyz_local*dist_surface_vec
                        levelset_norm%sf(i, j, k, ib_patch_id, :) = &
                            levelset_norm%sf(i, j, k, ib_patch_id, :)/ &
                            norm2(levelset_norm%sf(i, j, k, ib_patch_id, :))
                    end if
                    levelset_norm%sf(i, j, k, ib_patch_id, :) = &
                        matmul(rotation, levelset_norm%sf(i, j, k, ib_patch_id, :))
                end do
            end do
        end do

    end subroutine s_cylinder_levelset

end module m_compute_levelset
