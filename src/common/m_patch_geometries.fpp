!>
!! @file
!! @brief Contains module m_ibm

#:include 'macros.fpp'

!> @brief Contains helper functions specific to various patch gemoetries for determining if a grid cell lies inside of or outside of
!! a patch geometry
module m_patch_geometries

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_helper
    use m_helper_basic
    use m_constants
    use m_model

    implicit none

    public :: f_is_inside_sphere, f_is_inside_cylinder, f_is_inside_cuboid, f_is_inside_airfoil, f_is_inside_ellipse, &
        & s_compute_rotation_matrix

contains

    !> Check if the x, y, and z coordinates would be located inside a sphere with the patch_id's radius
    function f_is_inside_sphere(x, y, z, radius) result(is_inside)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: radius, x, y, z
        logical              :: is_inside

        is_inside = x**2 + y**2 + z**2 <= radius**2

    end function f_is_inside_sphere

    !> Check which length of the cylinder is not default. Use that direction as the height and the other two coordinate
    ! values as the radius check
    function f_is_inside_cylinder(polar_x, polar_y, height, radius, length) result(is_inside)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: polar_x, polar_y, height, radius, length
        logical              :: is_inside

        ! check if the circular component of the cylinder is correct
        is_inside = polar_x**2 + polar_y**2 <= radius**2

        ! in 3D, also check the length of the cylinder
        if (num_dims == 3) is_inside = is_inside .and. -0.5_wp*length <= height .and. 0.5_wp*length >= height

    end function f_is_inside_cylinder

    !> Check if the x, y, and possibly z coordinates would be located inside a cuboid with the patch_id's lengths
    function f_is_inside_cuboid(x, y, z, length) result(is_inside)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)               :: x, y, z
        real(wp), dimension(3), intent(in) :: length
        logical                            :: is_inside

        ! check if x and y are inside the rectangle plane at z=0
        is_inside = -0.5_wp*length(1) <= x .and. 0.5_wp*length(1) >= x .and. -0.5_wp*length(2) <= y .and. 0.5_wp*length(2) >= y

        ! if we are in 3D, this is a cuboid and so we must also check the z axis
        if (num_dims == 3) is_inside = is_inside .and. -0.5_wp*length(3) <= z .and. 0.5_wp*length(3) >= z

    end function f_is_inside_cuboid

    !> Check if the x, y, are bounded by a NACA airfoil. Check if the z coordinate is inside the left and right edges of the
    !! airfoil, if set.
    function f_is_inside_airfoil(x, y, z, length, airfoil_id) result(is_inside)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: x, y, z, length
        integer, intent(in)  :: airfoil_id
        logical              :: is_inside
        integer              :: k
        real(wp)             :: f

        is_inside = .false.

        ! check the initial x bounds of the grid cell
        if (.not. (x >= 0._wp .and. x <= ib_airfoil(airfoil_id)%c)) return

        ! if we are in 3D, we must also check the z axis
        if (num_dims == 3 .and. (.not. (-0.5_wp*length <= z .and. 0.5_wp*length >= z))) return

        ! our check branches for the upper and lower half of the airfoil
        if (y >= 0._wp) then
            ! increment the iterator so we know where in the airfoil arrays to look
            k = 1
            do while (ib_airfoil_grids(airfoil_id)%upper(k)%x < x)
                k = k + 1
            end do

            ! If the values are approximately equivalent, skip the next check
            if (f_approx_equal(ib_airfoil_grids(airfoil_id)%upper(k)%x, x)) then
                if (y <= ib_airfoil_grids(airfoil_id)%upper(k)%y) is_inside = .true.
            else
                ! check if the y value is below the upper edge of the airfoil
                f = (ib_airfoil_grids(airfoil_id)%upper(k)%x - x)/(ib_airfoil_grids(airfoil_id)%upper(k)%x &
                     & - ib_airfoil_grids(airfoil_id)%upper(k - 1)%x)
                if (y <= ((1._wp - f)*ib_airfoil_grids(airfoil_id)%upper(k)%y + f*ib_airfoil_grids(airfoil_id)%upper(k - 1)%y)) &
                    & is_inside = .true.
            end if
        else
            ! increment the iterator so we know where in the airfoil arrays to look
            k = 1
            do while (ib_airfoil_grids(airfoil_id)%lower(k)%x < x)
                k = k + 1
            end do

            ! If the values are approximately equivalent, skip the next check
            if (f_approx_equal(ib_airfoil_grids(airfoil_id)%lower(k)%x, x)) then
                if (y >= ib_airfoil_grids(airfoil_id)%lower(k)%y) is_inside = .true.
            else
                ! check if the y value is above the lower edge of the airfoil
                f = (ib_airfoil_grids(airfoil_id)%lower(k)%x - x)/(ib_airfoil_grids(airfoil_id)%lower(k)%x &
                     & - ib_airfoil_grids(airfoil_id)%lower(k - 1)%x)
                if (y >= ((1._wp - f)*ib_airfoil_grids(airfoil_id)%lower(k)%y + f*ib_airfoil_grids(airfoil_id)%lower(k - 1)%y)) &
                    & is_inside = .true.
            end if
        end if

    end function f_is_inside_airfoil

    function f_is_inside_ellipse(x, y, length) result(is_inside)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)               :: x, y
        real(wp), dimension(3), intent(in) :: length
        logical                            :: is_inside

        ! Ellipse condition (x/a)^2 + (y/b)^2 <= 1
        is_inside = (x/(0.5_wp*length(1)))**2 + (y/(0.5_wp*length(2)))**2 <= 1._wp

    end function f_is_inside_ellipse

    !> Compute a rotation matrix for converting to the rotating frame of the boundary
    subroutine s_compute_rotation_matrix(angles, rotation_matrix, rotation_matrix_inverse)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(1:3), intent(in)      :: angles
        real(wp), dimension(1:3,1:3), intent(out) :: rotation_matrix
        real(wp), dimension(1:3,1:3), intent(out) :: rotation_matrix_inverse
        real(wp), dimension(3, 3, 3)              :: rotation
        real(wp)                                  :: angle

        ! construct the x, y, and z rotation matrices

        if (num_dims == 3) then
            ! also compute the x and y axes in 3D
            angle = angles(1)
            rotation(1, 1,:) = [1._wp, 0._wp, 0._wp]
            rotation(1, 2,:) = [0._wp, cos(angle), -sin(angle)]
            rotation(1, 3,:) = [0._wp, sin(angle), cos(angle)]

            angle = angles(2)
            rotation(2, 1,:) = [cos(angle), 0._wp, sin(angle)]
            rotation(2, 2,:) = [0._wp, 1._wp, 0._wp]
            rotation(2, 3,:) = [-sin(angle), 0._wp, cos(angle)]

            ! apply the y rotation to the x rotation
            rotation_matrix(:,:) = matmul(rotation(1,:,:), rotation(2,:,:))
            rotation_matrix_inverse(:,:) = matmul(transpose(rotation(2,:,:)), transpose(rotation(1,:,:)))
        end if

        ! z component first, since it applies in 2D and 3D
        angle = angles(3)
        rotation(3, 1,:) = [cos(angle), -sin(angle), 0._wp]
        rotation(3, 2,:) = [sin(angle), cos(angle), 0._wp]
        rotation(3, 3,:) = [0._wp, 0._wp, 1._wp]

        if (num_dims == 3) then
            ! apply the z rotation to the xy rotation in 3D
            rotation_matrix(:,:) = matmul(rotation_matrix(:,:), rotation(3,:,:))
            rotation_matrix_inverse(:,:) = matmul(transpose(rotation(3,:,:)), rotation_matrix_inverse(:,:))
        else
            ! write out only the z rotation in 2D
            rotation_matrix(:,:) = rotation(3,:,:)
            rotation_matrix_inverse(:,:) = transpose(rotation(3,:,:))
        end if

    end subroutine s_compute_rotation_matrix

end module m_patch_geometries
