!>
!! @file
!! @brief Contains module m_patches

#:include 'case.fpp'
#:include 'ExtrusionHardcodedIC.fpp'
#:include '1dHardcodedIC.fpp'
#:include '2dHardcodedIC.fpp'
#:include '3dHardcodedIC.fpp'
#:include 'macros.fpp'

module m_ib_patches

    use m_model                 ! Subroutine(s) related to STL files

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     !< Definitions of the global parameters

    use m_helper_basic          !< Functions to compare floating point numbers

    use m_helper

    use m_mpi_common

    implicit none

    private; public :: s_apply_ib_patches, s_update_ib_rotation_matrix, s_instantiate_STL_models, models, f_convert_cyl_to_cart

    real(wp) :: x_centroid, y_centroid, z_centroid
    real(wp) :: length_x, length_y, length_z
    $:GPU_DECLARE(create='[x_centroid, y_centroid, z_centroid]')
    $:GPU_DECLARE(create='[length_x, length_y, length_z]')

    integer :: smooth_patch_id
    real(wp) :: smooth_coeff
    $:GPU_DECLARE(create='[smooth_patch_id, smooth_coeff]')
    !! These variables are analogous in both meaning and use to the similarly
    !! named components in the ic_patch_parameters type (see m_derived_types.f90
    !! for additional details). They are employed as a means to more concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    real(wp) :: cart_x, cart_y, cart_z
    real(wp) :: sph_phi !<
    $:GPU_DECLARE(create='[cart_x, cart_y, cart_z, sph_phi]')
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates

    type(bounds_info) :: x_boundary, y_boundary, z_boundary
    $:GPU_DECLARE(create='[x_boundary, y_boundary, z_boundary]')
    !! These variables combine the centroid and length parameters associated with
    !! a particular patch to yield the locations of the patch boundaries in the
    !! x-, y- and z-coordinate directions. They are used as a means to concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    character(len=5) :: istr ! string to store int to string result for error checking

    type(t_model_array), allocatable, target :: models(:)
    $:GPU_DECLARE(create='[models]')
    !! array of STL models that can be allocated and then used in IB marker and levelset compute

contains

    impure subroutine s_apply_ib_patches(ib_markers_sf)

        integer, dimension(:, :, :), intent(inout), optional :: ib_markers_sf

        integer :: i

        !  3D Patch Geometries
        if (p > 0) then

            !> IB Patches
            !> @{
            ! Spherical patch
            do i = 1, num_ibs

                if (patch_ib(i)%geometry == 8) then
                    call s_ib_sphere(i, ib_markers_sf)
                elseif (patch_ib(i)%geometry == 9) then
                    call s_ib_cuboid(i, ib_markers_sf)
                elseif (patch_ib(i)%geometry == 10) then
                    call s_ib_cylinder(i, ib_markers_sf)
                elseif (patch_ib(i)%geometry == 11) then
                    call s_ib_3D_airfoil(i, ib_markers_sf)
                    ! STL+IBM patch
                elseif (patch_ib(i)%geometry == 12) then
                    call s_ib_model(i, ib_markers_sf)
                end if
            end do
            !> @}

            ! 2D Patch Geometries
        elseif (n > 0) then

            !> IB Patches
            !> @{
            do i = 1, num_ibs
                if (patch_ib(i)%geometry == 2) then
                    call s_ib_circle(i, ib_markers_sf)
                elseif (patch_ib(i)%geometry == 3) then
                    call s_ib_rectangle(i, ib_markers_sf)
                elseif (patch_ib(i)%geometry == 4) then
                    call s_ib_airfoil(i, ib_markers_sf)
                    ! STL+IBM patch
                elseif (patch_ib(i)%geometry == 5) then
                    call s_ib_model(i, ib_markers_sf)
                elseif (patch_ib(i)%geometry == 6) then
                    call s_ib_ellipse(i, ib_markers_sf)
                end if
            end do
            !> @}

        end if

    end subroutine s_apply_ib_patches

    subroutine s_instantiate_STL_models()

        ! Variables for IBM+STL
        real(wp) :: normals(1:3) !< Boundary normal buffer
        integer :: boundary_vertex_count, boundary_edge_count, total_vertices !< Boundary vertex
        real(wp), allocatable, dimension(:, :, :) :: boundary_v !< Boundary vertex buffer
        real(wp), allocatable, dimension(:, :) :: interpolated_boundary_v !< Interpolated vertex buffer
        real(wp) :: distance !< Levelset distance buffer
        logical :: interpolate !< Logical variable to determine whether or not the model should be interpolated

        integer :: i, j, k !< Generic loop iterators
        integer :: patch_id

        type(t_bbox) :: bbox, bbox_old
        type(t_model) :: model
        type(ic_model_parameters) :: params

        real(wp) :: eta
        real(wp), dimension(1:3) :: point, model_center
        real(wp) :: grid_mm(1:3, 1:2)

        real(wp), dimension(1:4, 1:4) :: transform, transform_n

        call random_seed()

        do patch_id = 1, num_ibs
            if (patch_ib(patch_id)%geometry == 5 .or. patch_ib(patch_id)%geometry == 12) then
                allocate (models(patch_id)%model)
                print *, " * Reading model: "//trim(patch_ib(patch_id)%model_filepath)

                model = f_model_read(patch_ib(patch_id)%model_filepath)
                params%scale(:) = patch_ib(patch_id)%model_scale(:)
                params%translate(:) = patch_ib(patch_id)%model_translate(:)
                params%rotate(:) = patch_ib(patch_id)%model_rotate(:)
                params%spc = patch_ib(patch_id)%model_spc
                params%threshold = patch_ib(patch_id)%model_threshold

                if (f_approx_equal(dot_product(params%scale, params%scale), 0._wp)) then
                    params%scale(:) = 1._wp
                end if

                if (proc_rank == 0) then
                    print *, " * Transforming model."
                end if

                ! Get the model center before transforming the model
                bbox_old = f_create_bbox(model)
                model_center(1:3) = (bbox_old%min(1:3) + bbox_old%max(1:3))/2._wp

                ! Compute the transform matrices for vertices and normals
                transform = f_create_transform_matrix(params, model_center)
                transform_n = f_create_transform_matrix(params)

                call s_transform_model(model, transform, transform_n)

                ! Recreate the bounding box after transformation
                bbox = f_create_bbox(model)

                ! Show the number of vertices in the original STL model
                if (proc_rank == 0) then
                    print *, ' * Number of input model vertices:', 3*model%ntrs
                end if

                call f_check_boundary(model, boundary_v, boundary_vertex_count, boundary_edge_count)

                ! Check if the model needs interpolation
                if (p > 0) then
                    call f_check_interpolation_3D(model, (/dx, dy, dz/), interpolate)
                else
                    call f_check_interpolation_2D(boundary_v, boundary_edge_count, (/minval(dx), minval(dy), 0._wp/), interpolate)
                end if

                ! Show the number of edges and boundary edges in 2D STL models
                if (proc_rank == 0 .and. p == 0) then
                    print *, ' * Number of 2D model boundary edges:', boundary_edge_count
                end if

                ! Interpolate the STL model along the edges (2D) and on triangle facets (3D)
                if (interpolate) then
                    if (proc_rank == 0) then
                        print *, ' * Interpolating STL vertices.'
                    end if

                    if (p > 0) then
                        call f_interpolate_3D(model, (/dx, dy, dz/), interpolated_boundary_v, total_vertices)
                    else
                        call f_interpolate_2D(boundary_v, boundary_edge_count, (/dx, dy, dz/), interpolated_boundary_v, total_vertices)
                    end if

                    if (proc_rank == 0) then
                        print *, ' * Total number of interpolated boundary vertices:', total_vertices
                    end if
                end if

                if (proc_rank == 0) then
                    write (*, "(A, 3(2X, F20.10))") "    > Model:  Min:", bbox%min(1:3)
                    write (*, "(A, 3(2X, F20.10))") "    >         Cen:", (bbox%min(1:3) + bbox%max(1:3))/2._wp
                    write (*, "(A, 3(2X, F20.10))") "    >         Max:", bbox%max(1:3)

                    grid_mm(1, :) = (/minval(x_cc(0:m)) - 0.e5_wp*dx(0), maxval(x_cc(0:m)) + 0.e5_wp*dx(m)/)
                    grid_mm(2, :) = (/minval(y_cc(0:n)) - 0.e5_wp*dy(0), maxval(y_cc(0:n)) + 0.e5_wp*dy(n)/)

                    if (p > 0) then
                        grid_mm(3, :) = (/minval(z_cc(0:p)) - 0.e5_wp*dz(0), maxval(z_cc(0:p)) + 0.e5_wp*dz(p)/)
                    else
                        grid_mm(3, :) = (/0._wp, 0._wp/)
                    end if

                    write (*, "(A, 3(2X, F20.10))") "    > Domain: Min:", grid_mm(:, 1)
                    write (*, "(A, 3(2X, F20.10))") "    >         Cen:", (grid_mm(:, 1) + grid_mm(:, 2))/2._wp
                    write (*, "(A, 3(2X, F20.10))") "    >         Max:", grid_mm(:, 2)
                end if

                models(patch_id)%model = model
                models(patch_id)%boundary_v = boundary_v
                models(patch_id)%boundary_edge_count = boundary_edge_count
                models(patch_id)%interpolate = interpolate
                if (interpolate) then
                    models(patch_id)%interpolated_boundary_v = interpolated_boundary_v
                    models(patch_id)%total_vertices = total_vertices
                end if

            end if
        end do

    end subroutine s_instantiate_STL_models

    !> The circular patch is a 2D geometry that may be used, for
        !!              example, in creating a bubble or a droplet. The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !! @param patch_id is the patch identifier
        !! @param ib_markers_sf Array to track patch ids
        !! @param ib True if this patch is an immersed boundary
    subroutine s_ib_circle(patch_id, ib_markers_sf)

        integer, intent(in) :: patch_id
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: ib_markers_sf

        real(wp), dimension(1:2) :: center
        real(wp) :: radius

        integer :: i, j, k !< Generic loop iterators

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information

        center(1) = patch_ib(patch_id)%x_centroid
        center(2) = patch_ib(patch_id)%y_centroid
        radius = patch_ib(patch_id)%radius

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.

        $:GPU_PARALLEL_LOOP(private='[i,j]', copy='[ib_markers_sf]',&
                  & copyin='[patch_id,center,radius]', collapse=2)
        do j = 0, n
            do i = 0, m
                if ((x_cc(i) - center(1))**2 &
                    + (y_cc(j) - center(2))**2 <= radius**2) &
                    then
                    ib_markers_sf(i, j, 0) = patch_id
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_circle

    !! @param patch_id is the patch identifier
    !! @param ib_markers_sf Array to track patch ids
    subroutine s_ib_airfoil(patch_id, ib_markers_sf)

        integer, intent(in) :: patch_id
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: ib_markers_sf

        real(wp) :: f, ca_in, pa, ma, ta
        real(wp) :: xa, yt, xu, yu, xl, yl, xc, yc, dycdxc, sin_c, cos_c
        integer :: i, j, k
        integer :: Np1, Np2

        real(wp), dimension(1:3) :: xy_local, offset !< x and y coordinates in local IB frame
        real(wp), dimension(1:2) :: center !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: inverse_rotation

        center(1) = patch_ib(patch_id)%x_centroid
        center(2) = patch_ib(patch_id)%y_centroid
        ca_in = patch_ib(patch_id)%c
        pa = patch_ib(patch_id)%p
        ma = patch_ib(patch_id)%m
        ta = patch_ib(patch_id)%t
        inverse_rotation(:, :) = patch_ib(patch_id)%rotation_matrix_inverse(:, :)
        offset(:) = patch_ib(patch_id)%centroid_offset(:)

        Np1 = int((pa*ca_in/dx(0))*20)
        Np2 = int(((ca_in - pa*ca_in)/dx(0))*20)
        Np = Np1 + Np2 + 1
        $:GPU_UPDATE(device='[Np]')

        if (.not. allocated(airfoil_grid_u)) then
            @:ALLOCATE(airfoil_grid_u(1:Np))
            @:ALLOCATE(airfoil_grid_l(1:Np))

            ! TODO :: The below instantiations are already handles by the loop below
            airfoil_grid_u(1)%x = 0._wp
            airfoil_grid_u(1)%y = 0._wp

            airfoil_grid_l(1)%x = 0._wp
            airfoil_grid_l(1)%y = 0._wp

            do i = 1, Np1 + Np2 - 1
                ! TODO :: This allocated the upper and lower airfoil arrays, and does not need to be performed each time the IB markers are updated. Place this as a separate subroutine.
                if (i <= Np1) then
                    xc = i*(pa*ca_in/Np1)
                    xa = xc/ca_in
                    yc = (ma/pa**2)*(2*pa*xa - xa**2)
                    dycdxc = (2*ma/pa**2)*(pa - xa)
                else
                    xc = pa*ca_in + (i - Np1)*((ca_in - pa*ca_in)/Np2)
                    xa = xc/ca_in
                    yc = (ma/(1 - pa)**2)*(1 - 2*pa + 2*pa*xa - xa**2)
                    dycdxc = (2*ma/(1 - pa)**2)*(pa - xa)
                end if

                yt = (5._wp*ta)*(0.2969_wp*xa**0.5_wp - 0.126_wp*xa - 0.3516_wp*xa**2._wp + 0.2843_wp*xa**3 - 0.1015_wp*xa**4)
                sin_c = dycdxc/(1 + dycdxc**2)**0.5_wp
                cos_c = 1/(1 + dycdxc**2)**0.5_wp

                xu = xa - yt*sin_c
                yu = yc + yt*cos_c

                xl = xa + yt*sin_c
                yl = yc - yt*cos_c

                xu = xu*ca_in
                yu = yu*ca_in

                xl = xl*ca_in
                yl = yl*ca_in

                airfoil_grid_u(i + 1)%x = xu
                airfoil_grid_u(i + 1)%y = yu

                airfoil_grid_l(i + 1)%x = xl
                airfoil_grid_l(i + 1)%y = yl

            end do

            airfoil_grid_u(Np)%x = ca_in
            airfoil_grid_u(Np)%y = 0._wp

            airfoil_grid_l(Np)%x = ca_in
            airfoil_grid_l(Np)%y = 0._wp

            $:GPU_UPDATE(device='[airfoil_grid_l,airfoil_grid_u]')

        end if

        $:GPU_PARALLEL_LOOP(private='[i,j,xy_local,k,f]', copy='[ib_markers_sf]',&
                  & copyin='[patch_id,center,inverse_rotation,offset,ma,ca_in,airfoil_grid_u,airfoil_grid_l]', collapse=2)
        do j = 0, n
            do i = 0, m
                xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp] ! get coordinate frame centered on IB
                xy_local = matmul(inverse_rotation, xy_local) ! rotate the frame into the IB's coordinates
                xy_local = xy_local - offset ! airfoils are a patch that require a centroid offset

                if (xy_local(1) >= 0._wp .and. xy_local(1) <= ca_in) then
                    xa = xy_local(1)/ca_in
                    if (xa <= pa) then
                        yc = (ma/pa**2)*(2*pa*xa - xa**2)
                        dycdxc = (2*ma/pa**2)*(pa - xa)
                    else
                        yc = (ma/(1 - pa)**2)*(1 - 2*pa + 2*pa*xa - xa**2)
                        dycdxc = (2*ma/(1 - pa)**2)*(pa - xa)
                    end if
                    if (xy_local(2) >= 0._wp) then
                        k = 1
                        do while (airfoil_grid_u(k)%x < xy_local(1) .and. k <= Np)
                            k = k + 1
                        end do
                        if (f_approx_equal(airfoil_grid_u(k)%x, xy_local(1))) then
                            if (xy_local(2) <= airfoil_grid_u(k)%y) then
                                !!IB
                                ib_markers_sf(i, j, 0) = patch_id
                            end if
                        else
                            f = (airfoil_grid_u(k)%x - xy_local(1))/(airfoil_grid_u(k)%x - airfoil_grid_u(k - 1)%x)
                            if (xy_local(2) <= ((1._wp - f)*airfoil_grid_u(k)%y + f*airfoil_grid_u(k - 1)%y)) then
                                !!IB
                                ib_markers_sf(i, j, 0) = patch_id
                            end if
                        end if
                    else
                        k = 1
                        do while (airfoil_grid_l(k)%x < xy_local(1))
                            k = k + 1
                        end do
                        if (f_approx_equal(airfoil_grid_l(k)%x, xy_local(1))) then
                            if (xy_local(2) >= airfoil_grid_l(k)%y) then
                                !!IB
                                ib_markers_sf(i, j, 0) = patch_id
                            end if
                        else
                            f = (airfoil_grid_l(k)%x - xy_local(1))/(airfoil_grid_l(k)%x - airfoil_grid_l(k - 1)%x)

                            if (xy_local(2) >= ((1._wp - f)*airfoil_grid_l(k)%y + f*airfoil_grid_l(k - 1)%y)) then
                                !!IB
                                ib_markers_sf(i, j, 0) = patch_id
                            end if
                        end if
                    end if
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_airfoil

    !! @param patch_id is the patch identifier
    !! @param ib_markers_sf Array to track patch ids
    !! @param ib True if this patch is an immersed boundary
    subroutine s_ib_3D_airfoil(patch_id, ib_markers_sf)

        integer, intent(in) :: patch_id
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: ib_markers_sf

        real(wp) :: lz, z_max, z_min, f, ca_in, pa, ma, ta, xa, yt, xu, yu, xl, yl, xc, yc, dycdxc, sin_c, cos_c
        integer :: i, j, k, l
        integer :: Np1, Np2

        real(wp), dimension(1:3) :: xyz_local, center, offset !< x, y, z coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: inverse_rotation

        center(1) = patch_ib(patch_id)%x_centroid
        center(2) = patch_ib(patch_id)%y_centroid
        center(3) = patch_ib(patch_id)%z_centroid
        lz = patch_ib(patch_id)%length_z
        ca_in = patch_ib(patch_id)%c
        pa = patch_ib(patch_id)%p
        ma = patch_ib(patch_id)%m
        ta = patch_ib(patch_id)%t
        inverse_rotation(:, :) = patch_ib(patch_id)%rotation_matrix_inverse(:, :)
        offset(:) = patch_ib(patch_id)%centroid_offset(:)

        z_max = lz/2
        z_min = -lz/2

        Np1 = int((pa*ca_in/dx(0))*20)
        Np2 = int(((ca_in - pa*ca_in)/dx(0))*20)
        Np = Np1 + Np2 + 1
        $:GPU_UPDATE(device='[Np]')

        if (.not. allocated(airfoil_grid_u)) then

            @:ALLOCATE(airfoil_grid_u(1:Np))
            @:ALLOCATE(airfoil_grid_l(1:Np))

            airfoil_grid_u(1)%x = 0._wp
            airfoil_grid_u(1)%y = 0._wp

            airfoil_grid_l(1)%x = 0._wp
            airfoil_grid_l(1)%y = 0._wp

            do i = 1, Np1 + Np2 - 1
                if (i <= Np1) then
                    xc = i*(pa*ca_in/Np1)
                    xa = xc/ca_in
                    yc = (ma/pa**2)*(2*pa*xa - xa**2)
                    dycdxc = (2*ma/pa**2)*(pa - xa)
                else
                    xc = pa*ca_in + (i - Np1)*((ca_in - pa*ca_in)/Np2)
                    xa = xc/ca_in
                    yc = (ma/(1 - pa)**2)*(1 - 2*pa + 2*pa*xa - xa**2)
                    dycdxc = (2*ma/(1 - pa)**2)*(pa - xa)
                end if

                yt = (5._wp*ta)*(0.2969_wp*xa**0.5_wp - 0.126_wp*xa - 0.3516_wp*xa**2._wp + 0.2843_wp*xa**3 - 0.1015_wp*xa**4)
                sin_c = dycdxc/(1 + dycdxc**2)**0.5_wp
                cos_c = 1/(1 + dycdxc**2)**0.5_wp

                xu = xa - yt*sin_c
                yu = yc + yt*cos_c

                xl = xa + yt*sin_c
                yl = yc - yt*cos_c

                xu = xu*ca_in
                yu = yu*ca_in

                xl = xl*ca_in
                yl = yl*ca_in

                airfoil_grid_u(i + 1)%x = xu
                airfoil_grid_u(i + 1)%y = yu

                airfoil_grid_l(i + 1)%x = xl
                airfoil_grid_l(i + 1)%y = yl

            end do

            airfoil_grid_u(Np)%x = ca_in
            airfoil_grid_u(Np)%y = 0._wp

            airfoil_grid_l(Np)%x = ca_in
            airfoil_grid_l(Np)%y = 0._wp

            $:GPU_UPDATE(device='[airfoil_grid_l,airfoil_grid_u]')
        end if

        $:GPU_PARALLEL_LOOP(private='[i,j,l,xyz_local,k,f]', copy='[ib_markers_sf]',&
                  & copyin='[patch_id,center,inverse_rotation,offset,ma,ca_in,airfoil_grid_u,airfoil_grid_l,z_min,z_max]', collapse=3)
        do l = 0, p
            do j = 0, n
                do i = 0, m
                    xyz_local = [x_cc(i) - center(1), y_cc(j) - center(2), z_cc(l) - center(3)] ! get coordinate frame centered on IB
                    xyz_local = matmul(inverse_rotation, xyz_local) ! rotate the frame into the IB's coordinates
                    xyz_local = xyz_local - offset ! airfoils are a patch that require a centroid offset

                    if (xyz_local(3) >= z_min .and. xyz_local(3) <= z_max) then

                        if (xyz_local(1) >= 0._wp .and. xyz_local(1) <= ca_in) then
                            if (xyz_local(2) >= 0._wp) then
                                k = 1
                                do while (airfoil_grid_u(k)%x < xyz_local(1))
                                    k = k + 1
                                end do
                                if (f_approx_equal(airfoil_grid_u(k)%x, xyz_local(1))) then
                                    if (xyz_local(2) <= airfoil_grid_u(k)%y) then
                                        !IB
                                        ib_markers_sf(i, j, l) = patch_id
                                    end if
                                else
                                    f = (airfoil_grid_u(k)%x - xyz_local(1))/(airfoil_grid_u(k)%x - airfoil_grid_u(k - 1)%x)
                                    if (xyz_local(2) <= ((1._wp - f)*airfoil_grid_u(k)%y + f*airfoil_grid_u(k - 1)%y)) then
                                        !!IB
                                        ib_markers_sf(i, j, l) = patch_id
                                    end if
                                end if
                            else
                                k = 1
                                do while (airfoil_grid_l(k)%x < xyz_local(1))
                                    k = k + 1
                                end do
                                if (f_approx_equal(airfoil_grid_l(k)%x, xyz_local(1))) then
                                    if (xyz_local(2) >= airfoil_grid_l(k)%y) then
                                        !!IB
                                        ib_markers_sf(i, j, l) = patch_id
                                    end if
                                else
                                    f = (airfoil_grid_l(k)%x - xyz_local(1))/(airfoil_grid_l(k)%x - airfoil_grid_l(k - 1)%x)

                                    if (xyz_local(2) >= ((1._wp - f)*airfoil_grid_l(k)%y + f*airfoil_grid_l(k - 1)%y)) then
                                        !!IB
                                        ib_markers_sf(i, j, l) = patch_id
                                    end if
                                end if
                            end if
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_3D_airfoil

    !> The rectangular patch is a 2D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, in alignment with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x- and y-
        !!              coordinate directions are provided. Please note that the
        !!              rectangular patch DOES NOT allow for the smoothing of its
        !!              boundaries.
        !! @param patch_id is the patch identifier
        !! @param ib_markers_sf Array to track patch ids
        !! @param ib True if this patch is an immersed boundary
    subroutine s_ib_rectangle(patch_id, ib_markers_sf)

        integer, intent(in) :: patch_id
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: ib_markers_sf

        integer :: i, j, k !< generic loop iterators
        real(wp) :: pi_inf, gamma, lit_gamma !< Equation of state parameters
        real(wp), dimension(1:3) :: xy_local !< x and y coordinates in local IB frame
        real(wp), dimension(1:2) :: length, center !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: inverse_rotation

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1._wp + gamma)/gamma

        ! Transferring the rectangle's centroid and length information
        center(1) = patch_ib(patch_id)%x_centroid
        center(2) = patch_ib(patch_id)%y_centroid
        length(1) = patch_ib(patch_id)%length_x
        length(2) = patch_ib(patch_id)%length_y
        inverse_rotation(:, :) = patch_ib(patch_id)%rotation_matrix_inverse(:, :)

        ! Checking whether the rectangle covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        $:GPU_PARALLEL_LOOP(private='[i,j, xy_local]', copy='[ib_markers_sf]',&
                  & copyin='[patch_id,center,length,inverse_rotation,x_cc,y_cc]', collapse=2)
        do j = 0, n
            do i = 0, m
                ! get the x and y coordinates in the local IB frame
                xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
                xy_local = matmul(inverse_rotation, xy_local)

                if (-0.5_wp*length(1) <= xy_local(1) .and. &
                    0.5_wp*length(1) >= xy_local(1) .and. &
                    -0.5_wp*length(2) <= xy_local(2) .and. &
                    0.5_wp*length(2) >= xy_local(2)) then

                    ! Updating the patch identities bookkeeping variable
                    ib_markers_sf(i, j, 0) = patch_id

                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_rectangle

    !>          The spherical patch is a 3D geometry that may be used,
        !!              for example, in creating a bubble or a droplet. The patch
        !!              geometry is well-defined when its centroid and radius are
        !!              provided. Please note that the spherical patch DOES allow
        !!              for the smoothing of its boundary.
        !! @param patch_id is the patch identifier
        !! @param ib_markers_sf Array to track patch ids
        !! @param ib True if this patch is an immersed boundary
    subroutine s_ib_sphere(patch_id, ib_markers_sf)

        integer, intent(in) :: patch_id
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: ib_markers_sf

        ! Generic loop iterators
        integer :: i, j, k
        real(wp) :: radius
        real(wp), dimension(1:3) :: center

        !! Variables to initialize the pressure field that corresponds to the
            !! bubble-collapse test case found in Tiwari et al. (2013)

        ! Transferring spherical patch's radius, centroid, smoothing patch
        ! identity and smoothing coefficient information
        center(1) = patch_ib(patch_id)%x_centroid
        center(2) = patch_ib(patch_id)%y_centroid
        center(3) = patch_ib(patch_id)%z_centroid
        radius = patch_ib(patch_id)%radius

        ! Checking whether the sphere covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        $:GPU_PARALLEL_LOOP(private='[i,j,k,cart_y,cart_z]', copy='[ib_markers_sf]',&
                  & copyin='[patch_id,center,radius]', collapse=3)
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if
                    ! Updating the patch identities bookkeeping variable
                    if (((x_cc(i) - center(1))**2 &
                         + (cart_y - center(2))**2 &
                         + (cart_z - center(3))**2 <= radius**2)) then
                        ib_markers_sf(i, j, k) = patch_id
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_sphere

    !> The cuboidal patch is a 3D geometry that may be used, for
        !!              example, in creating a solid boundary, or pre-/post-shock
        !!              region, which is aligned with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x-, y- and
        !!              z-coordinate directions are provided. Please notice that
        !!              the cuboidal patch DOES NOT allow for the smearing of its
        !!              boundaries.
        !! @param patch_id is the patch identifier
        !! @param ib_markers_sf Array to track patch ids
    subroutine s_ib_cuboid(patch_id, ib_markers_sf)

        integer, intent(in) :: patch_id
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: ib_markers_sf

        integer :: i, j, k !< Generic loop iterators
        real(wp), dimension(1:3) :: xyz_local, center, length !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: inverse_rotation

        ! Transferring the cuboid's centroid and length information
        center(1) = patch_ib(patch_id)%x_centroid
        center(2) = patch_ib(patch_id)%y_centroid
        center(3) = patch_ib(patch_id)%z_centroid
        length(1) = patch_ib(patch_id)%length_x
        length(2) = patch_ib(patch_id)%length_y
        length(3) = patch_ib(patch_id)%length_z
        inverse_rotation(:, :) = patch_ib(patch_id)%rotation_matrix_inverse(:, :)

        ! Checking whether the cuboid covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! to that cell. If both queries check out, the primitive variables
        ! of the current patch are assigned to this cell.
        $:GPU_PARALLEL_LOOP(private='[i,j,k,xyz_local,cart_y,cart_z]', copy='[ib_markers_sf]',&
                  & copyin='[patch_id,center,length,inverse_rotation]', collapse=3)
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        ! TODO :: This does not work and is not covered by any tests. This should be fixed
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if
                    xyz_local = [x_cc(i), cart_y, cart_z] - center ! get coordinate frame centered on IB
                    xyz_local = matmul(inverse_rotation, xyz_local) ! rotate the frame into the IB's coordinates

                    if (-0.5*length(1) <= xyz_local(1) .and. &
                        0.5*length(1) >= xyz_local(1) .and. &
                        -0.5*length(2) <= xyz_local(2) .and. &
                        0.5*length(2) >= xyz_local(2) .and. &
                        -0.5*length(3) <= xyz_local(3) .and. &
                        0.5*length(3) >= xyz_local(3)) then

                        ! Updating the patch identities bookkeeping variable
                        ib_markers_sf(i, j, k) = patch_id
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_cuboid

    !> The cylindrical patch is a 3D geometry that may be used,
        !!              for example, in setting up a cylindrical solid boundary
        !!              confinement, like a blood vessel. The geometry of this
        !!              patch is well-defined when the centroid, the radius and
        !!              the length along the cylinder's axis, parallel to the x-,
        !!              y- or z-coordinate direction, are provided. Please note
        !!              that the cylindrical patch DOES allow for the smoothing
        !!              of its lateral boundary.
        !! @param patch_id is the patch identifier
        !! @param ib_markers_sf Array to track patch ids
        !! @param ib True if this patch is an immersed boundary
    subroutine s_ib_cylinder(patch_id, ib_markers_sf)

        integer, intent(in) :: patch_id
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: ib_markers_sf

        integer :: i, j, k !< Generic loop iterators
        real(wp) :: radius
        real(wp), dimension(1:3) :: xyz_local, center, length !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: inverse_rotation

        ! Transferring the cylindrical patch's centroid, length, radius,
        ! smoothing patch identity and smoothing coefficient information

        center(1) = patch_ib(patch_id)%x_centroid
        center(2) = patch_ib(patch_id)%y_centroid
        center(3) = patch_ib(patch_id)%z_centroid
        length(1) = patch_ib(patch_id)%length_x
        length(2) = patch_ib(patch_id)%length_y
        length(3) = patch_ib(patch_id)%length_z
        radius = patch_ib(patch_id)%radius
        inverse_rotation(:, :) = patch_ib(patch_id)%rotation_matrix_inverse(:, :)

        ! Checking whether the cylinder covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        $:GPU_PARALLEL_LOOP(private='[i,j,k,xyz_local,cart_y,cart_z]', copy='[ib_markers_sf]',&
                  & copyin='[patch_id,center,length,radius,inverse_rotation]', collapse=3)
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if
                    xyz_local = [x_cc(i), cart_y, cart_z] - center ! get coordinate frame centered on IB
                    xyz_local = matmul(inverse_rotation, xyz_local) ! rotate the frame into the IB's coordinates

                    if (((.not. f_is_default(length(1)) .and. &
                          xyz_local(2)**2 &
                          + xyz_local(3)**2 <= radius**2 .and. &
                          -0.5_wp*length(1) <= xyz_local(1) .and. &
                          0.5_wp*length(1) >= xyz_local(1)) &
                         .or. &
                         (.not. f_is_default(length(2)) .and. &
                          xyz_local(1)**2 &
                          + xyz_local(3)**2 <= radius**2 .and. &
                          -0.5_wp*length(2) <= xyz_local(2) .and. &
                          0.5_wp*length(2) >= xyz_local(2)) &
                         .or. &
                         (.not. f_is_default(length(3)) .and. &
                          xyz_local(1)**2 &
                          + xyz_local(2)**2 <= radius**2 .and. &
                          -0.5_wp*length(3) <= xyz_local(3) .and. &
                          0.5_wp*length(3) >= xyz_local(3)))) then

                        ! Updating the patch identities bookkeeping variable
                        ib_markers_sf(i, j, k) = patch_id
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_cylinder

    subroutine s_ib_ellipse(patch_id, ib_markers_sf)

        integer, intent(in) :: patch_id
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: ib_markers_sf

        integer :: i, j, k !< generic loop iterators
        real(wp), dimension(1:3) :: xy_local !< x and y coordinates in local IB frame
        real(wp), dimension(1:2) :: ellipse_coeffs !< a and b in the ellipse coefficients
        real(wp), dimension(1:2) :: center !< x and y coordinates in local IB frame
        real(wp), dimension(1:3, 1:3) :: inverse_rotation

        ! Transferring the ellipse's centroid and length information
        center(1) = patch_ib(patch_id)%x_centroid
        center(2) = patch_ib(patch_id)%y_centroid
        ellipse_coeffs(1) = 0.5_wp*patch_ib(patch_id)%length_x
        ellipse_coeffs(2) = 0.5_wp*patch_ib(patch_id)%length_y
        inverse_rotation(:, :) = patch_ib(patch_id)%rotation_matrix_inverse(:, :)

        ! Checking whether the ellipse covers a particular cell in the
        ! domain
        $:GPU_PARALLEL_LOOP(private='[i,j, xy_local]', copy='[ib_markers_sf]',&
                  & copyin='[patch_id,center,ellipse_coeffs,inverse_rotation,x_cc,y_cc]', collapse=2)
        do j = 0, n
            do i = 0, m
                ! get the x and y coordinates in the local IB frame
                xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
                xy_local = matmul(inverse_rotation, xy_local)

                ! Ellipse condition (x/a)^2 + (y/b)^2 <= 1
                if ((xy_local(1)/ellipse_coeffs(1))**2 + (xy_local(2)/ellipse_coeffs(2))**2 <= 1._wp) then
                    ! Updating the patch identities bookkeeping variable
                    ib_markers_sf(i, j, 0) = patch_id
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_ellipse

    !> The STL patch is a 2/3D geometry that is imported from an STL file.
    !! @param patch_id is the patch identifier
    !! @param ib_markers_sf Array to track patch ids
    !! @param STL_levelset STL levelset
    !! @param STL_levelset_norm STL levelset normals
    subroutine s_ib_model(patch_id, ib_markers_sf)

        integer, intent(in) :: patch_id
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: ib_markers_sf

        integer :: i, j, k !< Generic loop iterators

        type(t_model), pointer :: model

        real(wp) :: eta
        real(wp), dimension(1:3) :: point, local_point, offset
        real(wp), dimension(1:3) :: center, xyz_local
        real(wp), dimension(1:3, 1:3) :: inverse_rotation

        model => models(patch_id)%model
        center = 0._wp
        center(1) = patch_ib(patch_id)%x_centroid
        center(2) = patch_ib(patch_id)%y_centroid
        if (p > 0) center(3) = patch_ib(patch_id)%z_centroid
        inverse_rotation(:, :) = patch_ib(patch_id)%rotation_matrix_inverse(:, :)
        offset(:) = patch_ib(patch_id)%centroid_offset(:)

        do i = 0, m
            do j = 0, n
                do k = 0, p

                    xyz_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
                    if (p > 0) then
                        xyz_local(3) = z_cc(k) - center(3)
                    end if
                    xyz_local = matmul(inverse_rotation, xyz_local)
                    xyz_local = xyz_local - offset

                    if (grid_geometry == 3) then
                        xyz_local = f_convert_cyl_to_cart(xyz_local)
                    end if

                    if (p == 0) then
                        eta = f_model_is_inside(model, xyz_local, (/dx(i), dy(j), 0._wp/), patch_ib(patch_id)%model_spc)
                    else
                        eta = f_model_is_inside(model, xyz_local, (/dx(i), dy(j), dz(k)/), patch_ib(patch_id)%model_spc)
                    end if

                    ! Reading STL boundary vertices and compute the levelset and levelset_norm
                    if (eta > patch_ib(patch_id)%model_threshold) then
                        ib_markers_sf(i, j, k) = patch_id
                    end if

                end do
            end do
        end do

    end subroutine s_ib_model

    !> Subroutine that computes a rotation matrix for converting to the rotating frame of the boundary
    subroutine s_update_ib_rotation_matrix(patch_id)

        integer, intent(in) :: patch_id
        integer :: i

        real(wp), dimension(3, 3, 3) :: rotation
        real(wp) :: angle

        ! construct the x, y, and z rotation matrices
        if (num_dims == 3) then
            ! also compute the x and y axes in 3D
            angle = patch_ib(patch_id)%angles(1)
            rotation(1, 1, :) = [1._wp, 0._wp, 0._wp]
            rotation(1, 2, :) = [0._wp, cos(angle), -sin(angle)]
            rotation(1, 3, :) = [0._wp, sin(angle), cos(angle)]

            angle = patch_ib(patch_id)%angles(2)
            rotation(2, 1, :) = [cos(angle), 0._wp, sin(angle)]
            rotation(2, 2, :) = [0._wp, 1._wp, 0._wp]
            rotation(2, 3, :) = [-sin(angle), 0._wp, cos(angle)]

            ! apply the y rotation to the x rotation
            patch_ib(patch_id)%rotation_matrix(:, :) = matmul(rotation(1, :, :), rotation(2, :, :))
            patch_ib(patch_id)%rotation_matrix_inverse(:, :) = matmul(transpose(rotation(2, :, :)), transpose(rotation(1, :, :)))
        end if

        ! z component first, since it applies in 2D and 3D
        angle = patch_ib(patch_id)%angles(3)
        rotation(3, 1, :) = [cos(angle), -sin(angle), 0._wp]
        rotation(3, 2, :) = [sin(angle), cos(angle), 0._wp]
        rotation(3, 3, :) = [0._wp, 0._wp, 1._wp]

        if (num_dims == 3) then
            ! apply the z rotation to the xy rotation in 3D
            patch_ib(patch_id)%rotation_matrix(:, :) = matmul(patch_ib(patch_id)%rotation_matrix(:, :), rotation(3, :, :))
            patch_ib(patch_id)%rotation_matrix_inverse(:, :) = matmul(transpose(rotation(3, :, :)), patch_ib(patch_id)%rotation_matrix_inverse(:, :))
        else
            ! write out only the z rotation in 2D
            patch_ib(patch_id)%rotation_matrix(:, :) = rotation(3, :, :)
            patch_ib(patch_id)%rotation_matrix_inverse(:, :) = transpose(rotation(3, :, :))
        end if

    end subroutine s_update_ib_rotation_matrix

    subroutine s_convert_cylindrical_to_cartesian_coord(cyl_y, cyl_z)
        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: cyl_y, cyl_z

        cart_y = cyl_y*sin(cyl_z)
        cart_z = cyl_y*cos(cyl_z)

    end subroutine s_convert_cylindrical_to_cartesian_coord

    pure function f_convert_cyl_to_cart(cyl) result(cart)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(1:3), intent(in) :: cyl
        real(wp), dimension(1:3) :: cart

        cart = (/cyl(1), &
                 cyl(2)*sin(cyl(3)), &
                 cyl(2)*cos(cyl(3))/)

    end function f_convert_cyl_to_cart

    subroutine s_convert_cylindrical_to_spherical_coord(cyl_x, cyl_y)
        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(IN) :: cyl_x, cyl_y

        sph_phi = atan(cyl_y/cyl_x)

    end subroutine s_convert_cylindrical_to_spherical_coord

    !> Archimedes spiral function
    !! @param myth Angle
    !! @param offset Thickness
    !! @param a Starting position
    pure elemental function f_r(myth, offset, a)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: myth, offset, a
        real(wp) :: b
        real(wp) :: f_r

        !r(th) = a + b*th

        b = 2._wp*a/(2._wp*pi)
        f_r = a + b*myth + offset
    end function f_r

end module m_ib_patches
