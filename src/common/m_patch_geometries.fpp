!>
!! @file
!! @brief Contains module m_ibm

#:include 'macros.fpp'

!> @brief Contains helper functions specific to various patch gemoetries for determining if a grid cell lies inside of or outside of a patch geometry
module m_patch_geometries

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_helper
    use m_helper_basic
    use m_constants
    use m_model

    implicit none

    public :: f_is_inside_circle, f_is_inside_sphere, f_is_inside_cylinder, f_is_inside_rectangle, f_is_inside_cuboid, &
        & f_is_inside_airfoil, f_is_inside_airfoil_3D, f_is_inside_ellipse

contains

    !> Check if the x, and y coordinates would be located inside a circle with the patch_id's radius
    function f_is_inside_circle(patch_id, x, y) result(is_inside)

      $:GPU_ROUTINE(parallelism='[seq]')

      integer, intent(in) :: patch_id
      real(wp), intent(in) :: x, y
      logical :: is_inside

      is_inside = x**2 + y**2 <= patch_ib(patch_id)%radius**2

    end function f_is_inside_circle

    !> Check if the x, y, and z coordinates would be located inside a sphere with the patch_id's radius
    function f_is_inside_sphere(patch_id, x, y, z) result(is_inside)

      $:GPU_ROUTINE(parallelism='[seq]')

      integer, intent(in) :: patch_id
      real(wp), intent(in) :: x, y, z
      logical :: is_inside

      is_inside = x**2 + y**2 + z**2 <= patch_ib(patch_id)%radius**2

    end function f_is_inside_sphere

    !> Check which length of the cylinder is not default. Use that direction as the height and the other two coordinate
    ! values as the radius check
    function f_is_inside_cylinder(patch_id, x, y, z) result(is_inside)

      $:GPU_ROUTINE(parallelism='[seq]')

      integer, intent(in) :: patch_id
      real(wp), intent(in) :: x, y, z
      logical :: is_inside

      ! check if the cylinder is extended in the x direction
      is_inside = (.not. f_is_default(patch_ib(patch_id)%length_x)) .and. y**2 + z &
      & **2 <= patch_ib(patch_id)%radius**2 .and. -0.5_wp*patch_ib(patch_id)%length_x <= x &
      & .and. 0.5_wp*patch_ib(patch_id)%length_x >= x

      ! check if the cylinder is extended in the z direction
      is_inside = is_inside .or. ((.not. f_is_default(patch_ib(patch_id)%length_y)) .and. x**2 + z &
      & **2 <= patch_ib(patch_id)%radius**2 .and. -0.5_wp*patch_ib(patch_id)%length_y <= y &
      & .and. 0.5_wp*patch_ib(patch_id)%length_y >= y)

      ! check if the cylinder is extended in the z direction
      is_inside = is_inside .or. ((.not. f_is_default(patch_ib(patch_id)%length_z)) .and. x**2 + y &
      & **2 <= patch_ib(patch_id)%radius**2 .and. -0.5_wp*patch_ib(patch_id)%length_z <= z &
      & .and. 0.5_wp*patch_ib(patch_id)%length_z >= z)

      ! TODO :: This can easily reuse the f_is_inside_circle function with an additional length requirement

    end function f_is_inside_cylinder

    !> Check if the x, y, and possibly z coordinates would be located inside a cuboid with the patch_id's lengths
    function f_is_inside_cuboid(patch_id, x, y, z) result(is_inside)

      $:GPU_ROUTINE(parallelism='[seq]')

      integer, intent(in) :: patch_id
      real(wp), intent(in) :: x, y, z
      logical :: is_inside

      ! check if x and y are inside the rectangle plane at z=0
      is_inside = -0.5_wp*patch_ib(patch_id)%length_x <= x .and. 0.5_wp*patch_ib(patch_id)%length_x >= x &
      is_inside = is_inside .and. -0.5_wp*patch_ib(patch_id)%length_y <= y.and. 0.5_wp*patch_ib(patch_id)%length_y >= y

      ! if we are in 3D, this is a cuboid and so we must also check the z axis
      if (num_dims == 3) is_inside = is_inside .and. -0.5_wp*patch_ib(patch_id)%length_z <= z .and. 0.5_wp*patch_ib(patch_id)%length_z >= z

    end function f_is_inside_cuboid

    !> Check if the x, y, are bounded by a NACA airfoil. Check if the z coordinate is inside the left and right edges of the airfoil, if set.
    function f_is_inside_airfoil(patch_id, x, y, z) result(is_inside)

      $:GPU_ROUTINE(parallelism='[seq]')

      integer, intent(in) :: patch_id
      real(wp), intent(in) :: x, y, z
      logical :: is_inside

      integer :: k, airfoil_id
      real(wp) :: f

      airfoil_id = patch_ib(patch_id)%airfoil_id

      is_inside = .false.

      ! check the initial x bounds of the grid cell
      if (.not. (x >= 0._wp .and. x <= ib_airfoil(airfoil_id)%c)) return

      ! if we are in 3D, we must also check the z axis
      if (num_dims == 3 .and. (.not. (-0.5_wp*patch_ib(patch_id)%length_z <= z .and. 0.5_wp*patch_ib(patch_id)%length_z >= z))) return

      ! our check branches for the upper and lower half of the airfoil
      if (y >= 0._wp) then

          ! incriment the iterator so we know where in the airfoil arrays to look
          k = 1
          do while (ib_airfoil_grids(airfoil_id)%upper(k)%x < x)
              k = k + 1
          end do

          ! If the values are approximately equivilent, skip the next check
          if (f_approx_equal(ib_airfoil_grids(airfoil_id)%upper(k)%x, x)) then
              if (y <= ib_airfoil_grids(airfoil_id)%upper(k)%y) is_inside = .true.
          else

              ! check if the y value is below the upper edge of the airfoil
              f = (ib_airfoil_grids(airfoil_id)%upper(k)%x - x) &
                    & /(ib_airfoil_grids(airfoil_id)%upper(k)%x - ib_airfoil_grids(airfoil_id)%upper(k - 1)%x)
              if (y <= ((1._wp - f)*ib_airfoil_grids(airfoil_id)%upper(k)%y &
                  & + f*ib_airfoil_grids(airfoil_id)%upper(k - 1)%y)) is_inside = .true.

          end if
      else

          ! incriment the iterator so we know where in the airfoil arrays to look
          k = 1
          do while (ib_airfoil_grids(airfoil_id)%lower(k)%x < x)
              k = k + 1
          end do

          ! If the values are approximately equivilent, skip the next check
          if (f_approx_equal(ib_airfoil_grids(airfoil_id)%lower(k)%x, x)) then
              if (y >= ib_airfoil_grids(airfoil_id)%lower(k)%y) is_inside = .true.
          else

              ! check if the y value is above the lower edge of the airfoil
              f = (ib_airfoil_grids(airfoil_id)%lower(k)%x - x) &
                    & /(ib_airfoil_grids(airfoil_id)%lower(k)%x - ib_airfoil_grids(airfoil_id)%lower(k - 1)%x)
              if (y >= ((1._wp - f)*ib_airfoil_grids(airfoil_id)%lower(k)%y &
                  & + f*ib_airfoil_grids(airfoil_id)%lower(k - 1)%y)) is_inside = .true.

          end if
      end if

    end function f_is_inside_airfoil

    function f_is_inside_ellipse(patch_id, x, y) result(is_inside)

      $:GPU_ROUTINE(parallelism='[seq]')

      integer, intent(in) :: patch_id
      real(wp), intent(in) :: x, y
      logical :: is_inside

    end function f_is_inside_ellipse

  end module m_patch_geometries
