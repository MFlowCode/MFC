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
        & f_is_inside_patch, s_update_patch_rotation_matrix, s_get_ib_bound, s_get_patch_bounding_box, s_initialize_airfoil_grid

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

    !> Is xyz_local, a point already shifted to the patch centroid and rotated into the patch frame, inside patch? The single
    !! geometry dispatch over the IB patch numbering (2D: 2 circle, 3 rectangle, 4 airfoil, 5 STL, 6 ellipse; 3D: 8 sphere, 9
    !! cuboid, 10 cylinder, 11 airfoil, 12 STL) for callers that place things inside a patch shape.
    function f_is_inside_patch(patch, xyz_local) result(is_inside)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ib_patch_parameters), intent(in) :: patch
        real(wp), dimension(3), intent(in)    :: xyz_local
        logical                               :: is_inside
        real(wp), dimension(3)                :: length, offset_local

        length = [patch%length_x, patch%length_y, patch%length_z]
        offset_local = xyz_local - patch%centroid_offset

        select case (patch%geometry)
        case (2)  ! circle
            is_inside = f_is_inside_cylinder(xyz_local(1), xyz_local(2), 0._wp, patch%radius, 0._wp)
        case (3, 9)  ! rectangle, cuboid
            is_inside = f_is_inside_cuboid(xyz_local(1), xyz_local(2), xyz_local(3), length)
        case (4)  ! 2D airfoil
            is_inside = f_is_inside_airfoil(offset_local(1), offset_local(2), 0._wp, 0._wp, patch%airfoil_id)
        case (11)  ! 3D airfoil
            is_inside = f_is_inside_airfoil(offset_local(1), offset_local(2), offset_local(3), patch%length_z, patch%airfoil_id)
        case (5, 12)  ! STL model
            is_inside = f_model_is_inside(gpu_ntrs(patch%model_id), patch%model_id, &
                                          & offset_local) > stl_models(patch%model_id)%model_threshold
        case (6)  ! ellipse
            is_inside = f_is_inside_ellipse(xyz_local(1), xyz_local(2), length)
        case (8)  ! sphere
            is_inside = f_is_inside_sphere(xyz_local(1), xyz_local(2), xyz_local(3), patch%radius)
        case (10)  ! cylinder along x
            is_inside = f_is_inside_cylinder(xyz_local(2), xyz_local(3), xyz_local(1), patch%radius, patch%length_x)
        case default
            is_inside = .false.
        end select

    end function f_is_inside_patch

    !> Compute a patch's rotation matrix (patch frame to fluid frame) and its inverse from its angles, z then y then x in 3D
    subroutine s_update_patch_rotation_matrix(patch)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ib_patch_parameters), intent(inout) :: patch
        real(wp), dimension(3, 3, 3)             :: rotation
        real(wp)                                 :: angle

        ! construct the x, y, and z rotation matrices

        if (num_dims == 3) then
            ! also compute the x and y axes in 3D
            angle = patch%angles(1)
            rotation(1, 1,:) = [1._wp, 0._wp, 0._wp]
            rotation(1, 2,:) = [0._wp, cos(angle), -sin(angle)]
            rotation(1, 3,:) = [0._wp, sin(angle), cos(angle)]

            angle = patch%angles(2)
            rotation(2, 1,:) = [cos(angle), 0._wp, sin(angle)]
            rotation(2, 2,:) = [0._wp, 1._wp, 0._wp]
            rotation(2, 3,:) = [-sin(angle), 0._wp, cos(angle)]

            ! apply the y rotation to the x rotation
            patch%rotation_matrix(:,:) = matmul(rotation(1,:,:), rotation(2,:,:))
            patch%rotation_matrix_inverse(:,:) = matmul(transpose(rotation(2,:,:)), transpose(rotation(1,:,:)))
        end if

        ! z component first, since it applies in 2D and 3D
        angle = patch%angles(3)
        rotation(3, 1,:) = [cos(angle), -sin(angle), 0._wp]
        rotation(3, 2,:) = [sin(angle), cos(angle), 0._wp]
        rotation(3, 3,:) = [0._wp, 0._wp, 1._wp]

        if (num_dims == 3) then
            ! apply the z rotation to the xy rotation in 3D
            patch%rotation_matrix(:,:) = matmul(patch%rotation_matrix(:,:), rotation(3,:,:))
            patch%rotation_matrix_inverse(:,:) = matmul(transpose(rotation(3,:,:)), patch%rotation_matrix_inverse(:,:))
        else
            ! write out only the z rotation in 2D
            patch%rotation_matrix(:,:) = rotation(3,:,:)
            patch%rotation_matrix_inverse(:,:) = transpose(rotation(3,:,:))
        end if

    end subroutine s_update_patch_rotation_matrix

    subroutine s_get_ib_bound(patch, bound)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ib_patch_parameters), intent(in) :: patch
        real(wp), intent(out)                 :: bound
        real(wp), dimension(2)                :: lx, ly, lz

        if (patch%geometry == 2 .or. patch%geometry == 8) then
            ! circle and sphere geometries
            bound = patch%radius
        else if (patch%geometry == 3) then
            bound = 0.5_wp*sqrt(patch%length_x**2 + patch%length_y**2)
        else if (patch%geometry == 4 .or. patch%geometry == 11) then
            ! rectangular geometries
            bound = ib_airfoil(patch%airfoil_id)%c
        else if (patch%geometry == 5) then
            ! STL model geometry
            lx(1) = stl_bounding_boxes(patch%model_id, 1, 1)
            lx(2) = stl_bounding_boxes(patch%model_id, 1, 3)
            ly(1) = stl_bounding_boxes(patch%model_id, 2, 1)
            ly(2) = stl_bounding_boxes(patch%model_id, 2, 3)

            bound = 0.5_wp*sqrt((lx(2) - lx(1))**2 + (ly(2) - ly(1))**2)
        else if (patch%geometry == 6) then
            ! ellipse geometry
            bound = 0.5_wp*max(patch%length_x, patch%length_y)
        else if (patch%geometry == 9) then
            ! cuboid geometries
            bound = 0.5_wp*sqrt(patch%length_x**2 + patch%length_y**2 + patch%length_z**2)
        else if (patch%geometry == 10) then
            ! cylinder geometry
            bound = sqrt(patch%radius**2 + patch%length_x**2)
        else if (patch%geometry == 12) then
            ! Local-space bounding box extents (min=1, max=2 in the third index)
            lx(1) = stl_bounding_boxes(patch%model_id, 1, 1) + patch%centroid_offset(1)
            lx(2) = stl_bounding_boxes(patch%model_id, 1, 3) + patch%centroid_offset(1)
            ly(1) = stl_bounding_boxes(patch%model_id, 2, 1) + patch%centroid_offset(2)
            ly(2) = stl_bounding_boxes(patch%model_id, 2, 3) + patch%centroid_offset(2)
            lz(1) = stl_bounding_boxes(patch%model_id, 3, 1) + patch%centroid_offset(3)
            lz(2) = stl_bounding_boxes(patch%model_id, 3, 3) + patch%centroid_offset(3)

            bound = 0.5_wp*sqrt((lx(2) - lx(1))**2 + (ly(2) - ly(1))**2 + (lz(2) - lz(1))**2)
        end if

    end subroutine s_get_ib_bound

    !> World-space axis-aligned bounding box of a patch centred at center: the rotated corners of an STL model's own box, otherwise
    !! a cube of half-width s_get_ib_bound about the centre.
    subroutine s_get_patch_bounding_box(patch, center, bbox_min, bbox_max)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ib_patch_parameters), intent(in) :: patch
        real(wp), dimension(3), intent(in)    :: center
        real(wp), dimension(3), intent(out)   :: bbox_min, bbox_max
        real(wp), dimension(3)                :: local_corner, world_corner
        real(wp), dimension(2)                :: lx, ly, lz
        real(wp)                              :: bound
        integer                               :: cx, cy, cz

        if (patch%geometry == 5) then
            ! STL model geometry
            lx(1) = stl_bounding_boxes(patch%model_id, 1, 1) + patch%centroid_offset(1)
            lx(2) = stl_bounding_boxes(patch%model_id, 1, 3) + patch%centroid_offset(1)
            ly(1) = stl_bounding_boxes(patch%model_id, 2, 1) + patch%centroid_offset(2)
            ly(2) = stl_bounding_boxes(patch%model_id, 2, 3) + patch%centroid_offset(2)

            bbox_min = 1e12
            bbox_max = -1e12
            ! Enumerate all 4 corners of the local bounding box, rotate to world space, track world-space AABB
            do cx = 1, 2
                do cy = 1, 2
                    local_corner = [lx(cx), ly(cy), 0._wp]
                    world_corner = matmul(patch%rotation_matrix, local_corner) + center
                    bbox_min(1) = min(bbox_min(1), world_corner(1))
                    bbox_min(2) = min(bbox_min(2), world_corner(2))
                    bbox_max(1) = max(bbox_max(1), world_corner(1))
                    bbox_max(2) = max(bbox_max(2), world_corner(2))
                end do
            end do
        else if (patch%geometry == 12) then
            ! Local-space bounding box extents (min=1, max=2 in the third index)
            lx(1) = stl_bounding_boxes(patch%model_id, 1, 1) + patch%centroid_offset(1)
            lx(2) = stl_bounding_boxes(patch%model_id, 1, 3) + patch%centroid_offset(1)
            ly(1) = stl_bounding_boxes(patch%model_id, 2, 1) + patch%centroid_offset(2)
            ly(2) = stl_bounding_boxes(patch%model_id, 2, 3) + patch%centroid_offset(2)
            lz(1) = stl_bounding_boxes(patch%model_id, 3, 1) + patch%centroid_offset(3)
            lz(2) = stl_bounding_boxes(patch%model_id, 3, 3) + patch%centroid_offset(3)

            bbox_min = 1e12
            bbox_max = -1e12
            ! Enumerate all 8 corners of the local bounding box, rotate to world space, track world-space AABB
            do cx = 1, 2
                do cy = 1, 2
                    do cz = 1, 2
                        local_corner = [lx(cx), ly(cy), lz(cz)]
                        world_corner = matmul(patch%rotation_matrix, local_corner) + center
                        bbox_min(1) = min(bbox_min(1), world_corner(1))
                        bbox_min(2) = min(bbox_min(2), world_corner(2))
                        bbox_min(3) = min(bbox_min(3), world_corner(3))
                        bbox_max(1) = max(bbox_max(1), world_corner(1))
                        bbox_max(2) = max(bbox_max(2), world_corner(2))
                        bbox_max(3) = max(bbox_max(3), world_corner(3))
                    end do
                end do
            end do
        else
            ! All other IBs
            call s_get_ib_bound(patch, bound)
            bbox_min = center - bound
            bbox_max = center + bound
        end if

    end subroutine s_get_patch_bounding_box

    !> Builds the NACA surface segments of airfoil airfoil_id on the host, once, with about 10 surface points per surface_spacing
    !! along the chord (IB patches pass the smallest grid cell). f_is_inside_airfoil reads them; a caller that tests on the device
    !! copies them there itself, since pre_process keeps airfoils on the host.
    subroutine s_initialize_airfoil_grid(airfoil_id, surface_spacing)

        integer, intent(in)  :: airfoil_id
        real(wp), intent(in) :: surface_spacing
        integer              :: j, Np
        real(wp)             :: ca_in, pa, ma, ta
        real(wp)             :: xa, yc, dycdxc, yt, xu, yu, xl, yl, sin_c, cos_c
        real(wp)             :: edge(1:2)  !< segment vector from vertex 1 to vertex 2

        ca_in = ib_airfoil(airfoil_id)%c
        pa = ib_airfoil(airfoil_id)%p
        ma = ib_airfoil(airfoil_id)%m
        ta = ib_airfoil(airfoil_id)%t

        Np = int(10._wp*ca_in/surface_spacing) + 1
        ib_airfoil_grids(airfoil_id)%Np = Np

        if (.not. allocated(ib_airfoil_grids(airfoil_id)%upper)) then
            allocate (ib_airfoil_grids(airfoil_id)%upper(1:Np - 1,1:3,1:2))
            allocate (ib_airfoil_grids(airfoil_id)%lower(1:Np - 1,1:3,1:2))

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
        end if

    end subroutine s_initialize_airfoil_grid

end module m_patch_geometries
