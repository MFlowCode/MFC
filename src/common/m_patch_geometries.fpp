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

    public :: f_is_inside_circle, f_is_inside_sphere, f_is_inside_cylinder, f_is_inside_cuboid, &
        & f_is_inside_airfoil, f_is_inside_ellipse

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

    subroutine s_get_bounding_box_corner_distance(patch_id, bounding_distance)

      $:GPU_ROUTINE(parallelism='[seq]')

      integer, intent(in) :: patch_id
      real(wp), intent(out) :: bounding_distance

      ! TODO :: FILL ME IN

    end subroutine s_get_bounding_box_corner_distance

        !> Initialize the NACA surface grids for all airfoil IB patches. Must be called after the grid is established (so dx is valid)
    !! and before s_apply_ib_patches or s_apply_levelset.
    subroutine s_initialize_ib_airfoils()

      integer  :: i, j, airfoil_id
      integer  :: Np, Np1, Np2
      real(wp) :: ca_in, pa, ma, ta
      real(wp) :: xc, xa, yc, dycdxc, yt, xu, yu, xl, yl, sin_c, cos_c

      do i = 1, num_ibs
          if (patch_ib(i)%geometry /= 4 .and. patch_ib(i)%geometry /= 11) cycle

          airfoil_id = patch_ib(i)%airfoil_id
          ca_in = ib_airfoil(airfoil_id)%c
          pa = ib_airfoil(airfoil_id)%p
          ma = ib_airfoil(airfoil_id)%m
          ta = ib_airfoil(airfoil_id)%t

          Np1 = int((pa*ca_in/dx(0))*20)
          Np2 = int(((ca_in - pa*ca_in)/dx(0))*20)
          Np = Np1 + Np2 + 1
          ib_airfoil_grids(airfoil_id)%Np = Np
          $:GPU_UPDATE(device='[ib_airfoil_grids(airfoil_id)%Np]')

          if (.not. allocated(ib_airfoil_grids(airfoil_id)%upper)) then
              @:ALLOCATE(ib_airfoil_grids(airfoil_id)%upper(1:Np))
              @:ALLOCATE(ib_airfoil_grids(airfoil_id)%lower(1:Np))

              ib_airfoil_grids(airfoil_id)%upper(1)%x = 0._wp
              ib_airfoil_grids(airfoil_id)%upper(1)%y = 0._wp
              ib_airfoil_grids(airfoil_id)%lower(1)%x = 0._wp
              ib_airfoil_grids(airfoil_id)%lower(1)%y = 0._wp

              do j = 1, Np1 + Np2 - 1
                  if (j <= Np1) then
                      xc = j*(pa*ca_in/Np1)
                      xa = xc/ca_in
                      yc = (ma/pa**2)*(2*pa*xa - xa**2)
                      dycdxc = (2*ma/pa**2)*(pa - xa)
                  else
                      xc = pa*ca_in + (j - Np1)*((ca_in - pa*ca_in)/Np2)
                      xa = xc/ca_in
                      yc = (ma/(1 - pa)**2)*(1 - 2*pa + 2*pa*xa - xa**2)
                      dycdxc = (2*ma/(1 - pa)**2)*(pa - xa)
                  end if

                  yt = (5._wp*ta)*(0.2969_wp*xa**0.5_wp - 0.126_wp*xa - 0.3516_wp*xa**2._wp + 0.2843_wp*xa**3 - 0.1015_wp*xa**4)
                  sin_c = dycdxc/(1 + dycdxc**2)**0.5_wp
                  cos_c = 1/(1 + dycdxc**2)**0.5_wp

                  xu = (xa - yt*sin_c)*ca_in
                  yu = (yc + yt*cos_c)*ca_in
                  xl = (xa + yt*sin_c)*ca_in
                  yl = (yc - yt*cos_c)*ca_in

                  ib_airfoil_grids(airfoil_id)%upper(j + 1)%x = xu
                  ib_airfoil_grids(airfoil_id)%upper(j + 1)%y = yu
                  ib_airfoil_grids(airfoil_id)%lower(j + 1)%x = xl
                  ib_airfoil_grids(airfoil_id)%lower(j + 1)%y = yl
              end do

              ib_airfoil_grids(airfoil_id)%upper(Np)%x = ca_in
              ib_airfoil_grids(airfoil_id)%upper(Np)%y = 0._wp
              ib_airfoil_grids(airfoil_id)%lower(Np)%x = ca_in
              ib_airfoil_grids(airfoil_id)%lower(Np)%y = 0._wp

              $:GPU_UPDATE(device='[ib_airfoil_grids(airfoil_id)%upper, ib_airfoil_grids(airfoil_id)%lower]')
          end if
      end do

  end subroutine s_initialize_ib_airfoils

  end module m_patch_geometries
