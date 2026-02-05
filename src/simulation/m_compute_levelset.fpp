!>
!! @file m_compute_levelset.fpp
!! @brief Contains module m_compute_levelset

#:include 'macros.fpp'

!> @brief This module is used to handle all operations related to immersed
!!              boundary methods (IBMs)
module m_compute_levelset

    use m_ib_patches           !< The IB patch parameters

    use m_model                !< Subroutine(s) related to STL files

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    implicit none

    private; public :: s_apply_levelset

contains

    impure subroutine s_apply_levelset(gps, num_gps)

        type(ghost_point), dimension(:), intent(inout) :: gps
        integer, intent(in) :: num_gps

        integer :: i, patch_id, patch_geometry

        !  3D Patch Geometries
        if (p > 0) then

            $:GPU_PARALLEL_LOOP(private='[i]', copy='[gps]', copyin='[patch_ib,airfoil_grid_u,airfoil_grid_l]')
            do i = 1, num_gps

                patch_id = gps(i)%ib_patch_id
                patch_geometry = patch_ib(patch_id)%geometry

                if (patch_geometry == 8) then
                    call s_sphere_levelset(gps(i))
                elseif (patch_geometry == 9) then
                    call s_cuboid_levelset(gps(i))
                elseif (patch_geometry == 10) then
                    call s_cylinder_levelset(gps(i))
                elseif (patch_geometry == 11) then
                    call s_3d_airfoil_levelset(gps(i))
                    ! STL+IBM patch
                elseif (patch_geometry == 12) then
                    call s_model_levelset(gps(i))
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
            !> @}

            ! 2D Patch Geometries
        elseif (n > 0) then

            $:GPU_PARALLEL_LOOP(private='[i]', copy='[gps]', copyin='[patch_ib]')
            do i = 1, num_gps

                patch_id = gps(i)%ib_patch_id
                patch_geometry = patch_ib(patch_id)%geometry

                if (patch_geometry == 2) then
                    call s_circle_levelset(gps(i))
                elseif (patch_geometry == 3) then
                    call s_rectangle_levelset(gps(i))
                elseif (patch_geometry == 4) then
                    call s_airfoil_levelset(gps(i))
                    ! STL+IBM patch
                elseif (patch_geometry == 5) then
                    call s_model_levelset(gps(i))
                elseif (patch_geometry == 6) then
                    call s_ellipse_levelset(gps(i))
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
            !> @}

        end if

    end subroutine s_apply_levelset

    subroutine s_circle_levelset(gp)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(inout) :: gp

        real(wp) :: radius, dist
        real(wp), dimension(2) :: center
        real(wp), dimension(3) :: dist_vec

        integer :: i, j, ib_patch_id !< Loop index variables

        ib_patch_id = gp%ib_patch_id
        i = gp%loc(1)
        j = gp%loc(2)

        radius = patch_ib(ib_patch_id)%radius

        dist_vec(1) = x_cc(i) - patch_ib(ib_patch_id)%x_centroid
        dist_vec(2) = y_cc(j) - patch_ib(ib_patch_id)%y_centroid
        dist_vec(3) = 0._wp
        dist = sqrt(sum(dist_vec**2))

        gp%levelset = dist - radius
        if (f_approx_equal(dist, 0._wp)) then
            gp%levelset_norm = 0
        else
            gp%levelset_norm = dist_vec(:)/dist
        end if

    end subroutine s_circle_levelset

    subroutine s_airfoil_levelset(gp)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(inout) :: gp

        real(wp) :: dist, global_dist
        integer :: global_id
        real(wp), dimension(3) :: dist_vec

        real(wp), dimension(1:3) :: xy_local, offset !< x and y coordinates in local IB frame
        real(wp), dimension(1:2) :: center
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        integer :: i, j, k, ib_patch_id !< Loop index variables

        ib_patch_id = gp%ib_patch_id
        i = gp%loc(1)
        j = gp%loc(2)

        center(1) = patch_ib(ib_patch_id)%x_centroid
        center(2) = patch_ib(ib_patch_id)%y_centroid
        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)
        offset(:) = patch_ib(ib_patch_id)%centroid_offset(:)

        xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp] ! get coordinate frame centered on IB
        xy_local = matmul(inverse_rotation, xy_local) ! rotate the frame into the IB's coordinate
        xy_local = xy_local - offset ! airfoils are a patch that require a centroid offset

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

        gp%levelset = dist
        if (f_approx_equal(dist, 0._wp)) then
            gp%levelset_norm = 0._wp
        else
            gp%levelset_norm = matmul(rotation, dist_vec(:))/dist ! convert the normal vector back to global grid coordinates
        end if

    end subroutine s_airfoil_levelset

    subroutine s_3d_airfoil_levelset(gp)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(inout) :: gp

        real(wp) :: dist, dist_surf, dist_side, global_dist
        integer :: global_id
        real(wp) :: lz, z_max, z_min
        real(wp), dimension(3) :: dist_vec

        real(wp), dimension(1:3) :: xyz_local, center, offset !< x, y, z coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        real(wp) :: length_z

        integer :: i, j, k, l, ib_patch_id !< Loop index variables

        ib_patch_id = gp%ib_patch_id
        i = gp%loc(1)
        j = gp%loc(2)
        l = gp%loc(3)

        center(1) = patch_ib(ib_patch_id)%x_centroid
        center(2) = patch_ib(ib_patch_id)%y_centroid
        center(3) = patch_ib(ib_patch_id)%z_centroid
        lz = patch_ib(ib_patch_id)%length_z
        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)
        offset(:) = patch_ib(ib_patch_id)%centroid_offset(:)

        z_max = center(3) + lz/2
        z_min = center(3) - lz/2

        xyz_local = [x_cc(i) - center(1), y_cc(j) - center(2), z_cc(l) - center(3)] ! get coordinate frame centered on IB
        xyz_local = matmul(inverse_rotation, xyz_local) ! rotate the frame into the IB's coordinates
        xyz_local = xyz_local - offset ! airfoils are a patch that require a centroid offset

        if (xyz_local(2) >= center(2)) then
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
            gp%levelset = dist_side
            if (f_approx_equal(dist_side, abs(z_cc(l) - z_min))) then
                gp%levelset_norm = (/0._wp, 0._wp, -1._wp/)
            else
                gp%levelset_norm = (/0._wp, 0._wp, 1._wp/)
            end if
            gp%levelset_norm = matmul(rotation, gp%levelset_norm)
        else
            gp%levelset = dist_surf
            if (f_approx_equal(dist_surf, 0._wp)) then
                gp%levelset_norm = 0._wp
            else
                gp%levelset_norm = matmul(rotation, dist_vec(:)/dist_surf)
            end if
        end if

    end subroutine s_3d_airfoil_levelset

    !>  Initialize IBM module
    subroutine s_rectangle_levelset(gp)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(inout) :: gp

        real(wp) :: top_right(2), bottom_left(2)
        real(wp) :: min_dist
        real(wp) :: side_dists(4)

        real(wp) :: length_x, length_y
        real(wp), dimension(1:3) :: xy_local, dist_vec !< x and y coordinates in local IB frame
        real(wp), dimension(2) :: center !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        integer :: i, j, k !< Loop index variables
        integer :: idx !< Shortest path direction indicator
        integer :: ib_patch_id !< patch ID

        ib_patch_id = gp%ib_patch_id
        i = gp%loc(1)
        j = gp%loc(2)

        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        center(1) = patch_ib(ib_patch_id)%x_centroid
        center(2) = patch_ib(ib_patch_id)%y_centroid
        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)

        top_right(1) = length_x/2
        top_right(2) = length_y/2
        bottom_left(1) = -length_x/2
        bottom_left(2) = -length_y/2

        ! convert grid to local coordinates
        xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
        xy_local = matmul(inverse_rotation, xy_local)

        side_dists(1) = bottom_left(1) - xy_local(1)
        side_dists(2) = top_right(1) - xy_local(1)
        side_dists(3) = bottom_left(2) - xy_local(2)
        side_dists(4) = top_right(2) - xy_local(2)
        min_dist = side_dists(1)
        idx = 1

        do k = 2, 4
            if (abs(side_dists(k)) < abs(min_dist)) then
                idx = k
                min_dist = side_dists(idx)
            end if
        end do

        gp%levelset = side_dists(idx)
        dist_vec = 0._wp
        if (.not. f_approx_equal(side_dists(idx), 0._wp)) then
            if (idx == 1 .or. idx == 2) then
                ! vector points along the x axis
                dist_vec(1) = side_dists(idx)/abs(side_dists(idx))
            else
                ! vector points along the y axis
                dist_vec(2) = side_dists(idx)/abs(side_dists(idx))
            end if
            ! convert the normal vector back into the global coordinate system
            gp%levelset_norm = matmul(rotation, dist_vec)
        else
            gp%levelset_norm = 0._wp
        end if

    end subroutine s_rectangle_levelset

    subroutine s_ellipse_levelset(gp)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(inout) :: gp

        real(wp) :: ellipse_coeffs(2) ! a and b in the ellipse equation
        real(wp) :: quadratic_coeffs(3) ! A, B, C in the quadratic equation to compute levelset

        real(wp) :: length_x, length_y
        real(wp), dimension(1:3) :: xy_local, normal_vector !< x and y coordinates in local IB frame
        real(wp), dimension(2) :: center !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        integer :: i, j, k !< Loop index variables
        integer :: idx !< Shortest path direction indicator
        integer :: ib_patch_id !< patch ID

        ib_patch_id = gp%ib_patch_id
        i = gp%loc(1)
        j = gp%loc(2)

        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        center(1) = patch_ib(ib_patch_id)%x_centroid
        center(2) = patch_ib(ib_patch_id)%y_centroid
        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)

        ellipse_coeffs(1) = 0.5_wp*length_x
        ellipse_coeffs(2) = 0.5_wp*length_y

        xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
        xy_local = matmul(inverse_rotation, xy_local)

        ! we will get NaNs in the levelset if we compute this outside the ellipse
        if ((xy_local(1)/ellipse_coeffs(1))**2 + (xy_local(2)/ellipse_coeffs(2))**2 <= 1._wp) then

            normal_vector = xy_local
            normal_vector(2) = normal_vector(2)*(ellipse_coeffs(1)/ellipse_coeffs(2))**2._wp ! get the normal direction via the coordinate transformation method
            normal_vector = normal_vector/sqrt(dot_product(normal_vector, normal_vector)) ! normalize the vector
            gp%levelset_norm = matmul(rotation, normal_vector) ! save after rotating the vector to the global frame

            ! use the normal vector to set up the quadratic equation for the levelset, using A, B, and C in indices 1, 2, and 3
            quadratic_coeffs(1) = (normal_vector(1)/ellipse_coeffs(1))**2 + (normal_vector(2)/ellipse_coeffs(2))**2
            quadratic_coeffs(2) = 2._wp*((xy_local(1)*normal_vector(1)/(ellipse_coeffs(1)**2)) + (xy_local(2)*normal_vector(2)/(ellipse_coeffs(2)**2)))
            quadratic_coeffs(3) = (xy_local(1)/ellipse_coeffs(1))**2._wp + (xy_local(2)/ellipse_coeffs(2))**2._wp - 1._wp

            ! compute the levelset with the quadratic equation [ -B + sqrt(B^2 - 4AC) ] / 2A
            gp%levelset = -0.5_wp*(-quadratic_coeffs(2) + sqrt(quadratic_coeffs(2)**2._wp - 4._wp*quadratic_coeffs(1)*quadratic_coeffs(3)))/quadratic_coeffs(1)
        end if

    end subroutine s_ellipse_levelset

    subroutine s_cuboid_levelset(gp)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(inout) :: gp

        real(wp) :: Right, Left, Bottom, Top, Front, Back
        real(wp) :: min_dist
        real(wp) :: side_dists(6)

        real(wp), dimension(3) :: center
        real(wp) :: length_x, length_y, length_z
        real(wp), dimension(1:3) :: xyz_local, dist_vec !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        integer :: i, j, k !< Loop index variables
        integer :: ib_patch_id !< patch ID

        ib_patch_id = gp%ib_patch_id
        i = gp%loc(1)
        j = gp%loc(2)
        k = gp%loc(3)

        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        length_z = patch_ib(ib_patch_id)%length_z

        center(1) = patch_ib(ib_patch_id)%x_centroid
        center(2) = patch_ib(ib_patch_id)%y_centroid
        center(3) = patch_ib(ib_patch_id)%z_centroid

        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)

        Right = length_x/2
        Left = -length_x/2
        Top = length_y/2
        Bottom = -length_y/2
        Front = length_z/2
        Back = -length_z/2

        xyz_local = [x_cc(i), y_cc(j), z_cc(k)] - center ! get coordinate frame centered on IB
        xyz_local = matmul(inverse_rotation, xyz_local) ! rotate the frame into the IB's coordinate

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
        ! and this code should be cleaned up. It also means that
        ! rotating the box 90 degrees will cause tests to fail.
        dist_vec = 0._wp
        if (f_approx_equal(min_dist, abs(side_dists(1)))) then
            gp%levelset = side_dists(1)
            if (.not. f_approx_equal(side_dists(1), 0._wp)) then
                dist_vec(1) = side_dists(1)/abs(side_dists(1))
            end if

        else if (f_approx_equal(min_dist, abs(side_dists(2)))) then
            gp%levelset = side_dists(2)
            if (.not. f_approx_equal(side_dists(2), 0._wp)) then
                dist_vec(1) = -side_dists(2)/abs(side_dists(2))
            end if

        else if (f_approx_equal(min_dist, abs(side_dists(3)))) then
            gp%levelset = side_dists(3)
            if (.not. f_approx_equal(side_dists(3), 0._wp)) then
                dist_vec(2) = side_dists(3)/abs(side_dists(3))
            end if

        else if (f_approx_equal(min_dist, abs(side_dists(4)))) then
            gp%levelset = side_dists(4)
            if (.not. f_approx_equal(side_dists(4), 0._wp)) then
                dist_vec(2) = -side_dists(4)/abs(side_dists(4))
            end if

        else if (f_approx_equal(min_dist, abs(side_dists(5)))) then
            gp%levelset = side_dists(5)
            if (.not. f_approx_equal(side_dists(5), 0._wp)) then
                dist_vec(3) = side_dists(5)/abs(side_dists(5))
            end if

        else if (f_approx_equal(min_dist, abs(side_dists(6)))) then
            gp%levelset = side_dists(6)
            if (.not. f_approx_equal(side_dists(6), 0._wp)) then
                dist_vec(3) = -side_dists(6)/abs(side_dists(6))
            end if
        end if
        gp%levelset_norm = matmul(rotation, dist_vec)

    end subroutine s_cuboid_levelset

    subroutine s_sphere_levelset(gp)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(inout) :: gp

        real(wp) :: radius, dist
        real(wp), dimension(3) :: dist_vec, center

        integer :: i, j, k, ib_patch_id !< Loop index variables

        ib_patch_id = gp%ib_patch_id
        i = gp%loc(1)
        j = gp%loc(2)
        k = gp%loc(3)

        radius = patch_ib(ib_patch_id)%radius
        center(1) = patch_ib(ib_patch_id)%x_centroid
        center(2) = patch_ib(ib_patch_id)%y_centroid
        center(3) = patch_ib(ib_patch_id)%z_centroid

        dist_vec(1) = x_cc(i) - center(1)
        dist_vec(2) = y_cc(j) - center(2)
        dist_vec(3) = z_cc(k) - center(3)
        dist = sqrt(sum(dist_vec**2))
        gp%levelset = dist - radius
        if (f_approx_equal(dist, 0._wp)) then
            gp%levelset_norm = (/1, 0, 0/)
        else
            gp%levelset_norm = dist_vec(:)/dist
        end if

    end subroutine s_sphere_levelset

    subroutine s_cylinder_levelset(gp)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(inout) :: gp

        real(wp) :: radius
        real(wp), dimension(3) :: dist_sides_vec, dist_surface_vec, length
        real(wp), dimension(2) :: boundary
        real(wp) :: dist_side, dist_surface, side_pos
        integer :: i, j, k !< Loop index variables
        integer :: ib_patch_id !< patch ID

        real(wp), dimension(1:3) :: xyz_local, center !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: rotation, inverse_rotation

        ib_patch_id = gp%ib_patch_id
        i = gp%loc(1)
        j = gp%loc(2)
        k = gp%loc(3)

        radius = patch_ib(ib_patch_id)%radius
        center(1) = patch_ib(ib_patch_id)%x_centroid
        center(2) = patch_ib(ib_patch_id)%y_centroid
        center(3) = patch_ib(ib_patch_id)%z_centroid
        length(1) = patch_ib(ib_patch_id)%length_x
        length(2) = patch_ib(ib_patch_id)%length_y
        length(3) = patch_ib(ib_patch_id)%length_z

        inverse_rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix_inverse(:, :)
        rotation(:, :) = patch_ib(ib_patch_id)%rotation_matrix(:, :)

        if (.not. f_approx_equal(length(1), 0._wp)) then
            boundary(1) = -0.5_wp*length(1)
            boundary(2) = 0.5_wp*length(1)
            dist_sides_vec = (/1, 0, 0/)
            dist_surface_vec = (/0, 1, 1/)
        else if (.not. f_approx_equal(length(2), 0._wp)) then
            boundary(1) = -0.5_wp*length(2)
            boundary(2) = 0.5_wp*length(2)
            dist_sides_vec = (/0, 1, 0/)
            dist_surface_vec = (/1, 0, 1/)
        else if (.not. f_approx_equal(length(3), 0._wp)) then
            boundary(1) = -0.5_wp*length(3)
            boundary(2) = 0.5_wp*length(3)
            dist_sides_vec = (/0, 0, 1/)
            dist_surface_vec = (/1, 1, 0/)
        end if

        xyz_local = [x_cc(i), y_cc(j), z_cc(k)] - center ! get coordinate frame centered on IB
        xyz_local = matmul(inverse_rotation, xyz_local) ! rotate the frame into the IB's coordinates

        ! get distance to flat edge of cylinder
        side_pos = dot_product(xyz_local, dist_sides_vec)
        dist_side = min(abs(side_pos - boundary(1)), &
                        abs(boundary(2) - side_pos))
        ! get distance to curved side of cylinder
        dist_surface = norm2(xyz_local*dist_surface_vec) &
                       - radius

        if (dist_side < abs(dist_surface)) then
            ! if the closest edge is flat
            gp%levelset = -dist_side
            if (f_approx_equal(dist_side, abs(side_pos - boundary(1)))) then
                gp%levelset_norm = matmul(rotation, -dist_sides_vec)
            else
                gp%levelset_norm = matmul(rotation, dist_sides_vec)
            end if
        else
            gp%levelset = dist_surface
            xyz_local = xyz_local*dist_surface_vec
            xyz_local = xyz_local/max(norm2(xyz_local), sgm_eps)
            gp%levelset_norm = matmul(rotation, xyz_local)
        end if

    end subroutine s_cylinder_levelset

    !> The STL patch is a 2/3D geometry that is imported from an STL file.
    !! @param patch_id is the patch identifier
    !! @param STL_levelset STL levelset
    !! @param STL_levelset_norm STL levelset normals
    subroutine s_model_levelset(gp)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(inout) :: gp

        integer :: i, j, k, patch_id, boundary_edge_count, total_vertices
        type(t_model), pointer :: model
        real(wp), pointer, dimension(:, :, :) :: boundary_v
        real(wp), pointer, dimension(:, :) :: interpolated_boundary_v
        logical :: interpolate
        real(wp), dimension(1:3) :: point
        real(wp) :: normals(1:3) !< Boundary normal buffer
        real(wp) :: distance

        patch_id = gp%ib_patch_id
        i = gp%loc(1)
        j = gp%loc(2)
        k = gp%loc(3)

        ! load in model values
        model => models(patch_id)%model
        interpolate = models(patch_id)%interpolate
        boundary_edge_count = models(patch_id)%boundary_edge_count
        total_vertices = models(patch_id)%total_vertices
        boundary_v => models(patch_id)%boundary_v
        interpolated_boundary_v => models(patch_id)%interpolated_boundary_v

        ! determine where we are located in space
        point = (/x_cc(i), y_cc(j), 0._wp/)
        if (p > 0) then
            point(3) = z_cc(k)
        end if

        if (grid_geometry == 3) then
            point = f_convert_cyl_to_cart(point)
        end if

        ! 3D models
        if (p > 0) then

            ! Get the boundary normals and shortest distance between the cell center and the model boundary
            call f_distance_normals_3D(model, point, normals, distance)

            ! Get the shortest distance between the cell center and the interpolated model boundary
            if (interpolate) then
                gp%levelset = f_interpolated_distance(interpolated_boundary_v, total_vertices, point)
            else
                gp%levelset = distance
            end if

            ! Correct the sign of the levelset
            gp%levelset = -abs(gp%levelset)

            ! Assign the levelset_norm
            gp%levelset_norm = normals(1:3)
        else
            ! 2D models
            if (interpolate) then
                ! Get the shortest distance between the cell center and the model boundary
                gp%levelset = f_interpolated_distance(interpolated_boundary_v, total_vertices, point)
            else
                ! Get the shortest distance between the cell center and the interpolated model boundary
                gp%levelset = f_distance(boundary_v, boundary_edge_count, point)
            end if

            ! Correct the sign of the levelset
            gp%levelset = -abs(gp%levelset)

            ! Get the boundary normals
            call f_normals(boundary_v, &
                           boundary_edge_count, &
                           point, &
                           normals)

            ! Assign the levelset_norm
            gp%levelset_norm = normals(1:3)

        end if

    end subroutine s_model_levelset

end module m_compute_levelset
