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
        & f_is_inside_ellipsoid, s_build_airfoil_grid

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
        real(wp)             :: segment_fraction, surface_height
        real(wp)             :: vertex_1(1:2), vertex_2(1:2)

        is_inside = .false.

        ! check the initial x bounds of the grid cell
        if (.not. (x >= 0._wp .and. x <= ib_airfoil(airfoil_id)%c)) return

        ! if we are in 3D, we must also check the z axis
        if (num_dims == 3 .and. (.not. (-0.5_wp*length <= z .and. 0.5_wp*length >= z))) return

        ! find the segment spanning x on the surface bounding this half
        k = 1
        if (y >= 0._wp) then
            do while (ib_airfoil_grids(airfoil_id)%upper(k, 2, 1) < x)
                k = k + 1
            end do
            vertex_1 = ib_airfoil_grids(airfoil_id)%upper(k, 1,:)
            vertex_2 = ib_airfoil_grids(airfoil_id)%upper(k, 2,:)
        else
            do while (ib_airfoil_grids(airfoil_id)%lower(k, 2, 1) < x)
                k = k + 1
            end do
            vertex_1 = ib_airfoil_grids(airfoil_id)%lower(k, 1,:)
            vertex_2 = ib_airfoil_grids(airfoil_id)%lower(k, 2,:)
        end if

        ! inside if y is between the chord side and the surface height interpolated along the segment
        segment_fraction = (x - vertex_1(1))/(vertex_2(1) - vertex_1(1))
        surface_height = vertex_1(2) + segment_fraction*(vertex_2(2) - vertex_1(2))
        is_inside = merge(y <= surface_height, y >= surface_height, y >= 0._wp)

    end function f_is_inside_airfoil

    function f_is_inside_ellipse(x, y, length) result(is_inside)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)               :: x, y
        real(wp), dimension(3), intent(in) :: length
        logical                            :: is_inside

        ! Ellipse condition (x/a)^2 + (y/b)^2 <= 1
        is_inside = (x/(0.5_wp*length(1)))**2 + (y/(0.5_wp*length(2)))**2 <= 1._wp

    end function f_is_inside_ellipse

    !> Check if x, y, z lie inside the ellipsoid with full axis lengths length(1:3)
    function f_is_inside_ellipsoid(x, y, z, length) result(is_inside)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)               :: x, y, z
        real(wp), dimension(3), intent(in) :: length
        logical                            :: is_inside

        is_inside = (x/(0.5_wp*length(1)))**2 + (y/(0.5_wp*length(2)))**2 + (z/(0.5_wp*length(3)))**2 <= 1._wp

    end function f_is_inside_ellipsoid

    !> Build the NACA surface segments of ib_airfoil(airfoil_id) once, with 20 surface points per min_grid_spacing (the same on
    !! every rank). Shared by airfoil IB patches and airfoil initial-condition patches.
    subroutine s_build_airfoil_grid(airfoil_id, min_grid_spacing)

        integer, intent(in)  :: airfoil_id
        real(wp), intent(in) :: min_grid_spacing
        integer              :: j, Np
        real(wp)             :: ca_in, pa, ma, ta
        real(wp)             :: xa, yc, dycdxc, yt, xu, yu, xl, yl, sin_c, cos_c
        real(wp)             :: edge(1:2)  !< segment vector from vertex 1 to vertex 2

        ca_in = ib_airfoil(airfoil_id)%c
        pa = ib_airfoil(airfoil_id)%p
        ma = ib_airfoil(airfoil_id)%m
        ta = ib_airfoil(airfoil_id)%t

        Np = int(10._wp*ca_in/min_grid_spacing) + 1  ! 20 surface points per smallest cell
        ib_airfoil_grids(airfoil_id)%Np = Np
        $:GPU_UPDATE(device='[ib_airfoil_grids(airfoil_id)%Np]')

        if (.not. allocated(ib_airfoil_grids(airfoil_id)%upper)) then
            @:ALLOCATE(ib_airfoil_grids(airfoil_id)%upper(1:Np - 1, 1:3, 1:2))
            @:ALLOCATE(ib_airfoil_grids(airfoil_id)%lower(1:Np - 1, 1:3, 1:2))

            ! segments run from the leading edge to the trailing edge
            ib_airfoil_grids(airfoil_id)%upper(1, 1,:) = 0._wp
            ib_airfoil_grids(airfoil_id)%lower(1, 1,:) = 0._wp

            do j = 1, Np - 2
                ! cosine spacing clusters points at the leading and trailing edges, Katz & Plotkin (2001)
                ! "Low-Speed Aerodynamics", 2nd ed., Cambridge University Press
                xa = 0.5_wp*(1._wp - cos(pi*real(j, wp)/real(Np - 1, wp)))
                if (xa <= pa) then
                    yc = (ma/pa**2)*(2*pa*xa - xa**2)
                    dycdxc = (2*ma/pa**2)*(pa - xa)
                else
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

                ! each surface point ends segment j and starts segment j + 1
                ib_airfoil_grids(airfoil_id)%upper(j, 2,:) = [xu, yu]
                ib_airfoil_grids(airfoil_id)%upper(j + 1, 1,:) = [xu, yu]
                ib_airfoil_grids(airfoil_id)%lower(j, 2,:) = [xl, yl]
                ib_airfoil_grids(airfoil_id)%lower(j + 1, 1,:) = [xl, yl]
            end do

            ib_airfoil_grids(airfoil_id)%upper(Np - 1, 2,:) = [ca_in, 0._wp]
            ib_airfoil_grids(airfoil_id)%lower(Np - 1, 2,:) = [ca_in, 0._wp]

            ! outward unit normals: left of the edge on the upper surface, right of it on the lower
            do j = 1, Np - 1
                edge = ib_airfoil_grids(airfoil_id)%upper(j, 2,:) - ib_airfoil_grids(airfoil_id)%upper(j, 1,:)
                ib_airfoil_grids(airfoil_id)%upper(j, 3,:) = [-edge(2), edge(1)]/norm2(edge)
                edge = ib_airfoil_grids(airfoil_id)%lower(j, 2,:) - ib_airfoil_grids(airfoil_id)%lower(j, 1,:)
                ib_airfoil_grids(airfoil_id)%lower(j, 3,:) = [edge(2), -edge(1)]/norm2(edge)
            end do

            $:GPU_UPDATE(device='[ib_airfoil_grids(airfoil_id)%upper, ib_airfoil_grids(airfoil_id)%lower]')
        end if

    end subroutine s_build_airfoil_grid

end module m_patch_geometries
