!>
!! @file
!! @brief Contains module m_ib_patches

#:include 'case.fpp'
#:include 'ExtrusionHardcodedIC.fpp'
#:include '1dHardcodedIC.fpp'
#:include '2dHardcodedIC.fpp'
#:include '3dHardcodedIC.fpp'
#:include 'macros.fpp'

!> @brief Immersed boundary patch geometry constructors for 2D and 3D shapes
module m_ib_patches

    use m_patch_geometries
    use m_model  ! Subroutine(s) related to STL files
    use m_derived_types  ! Definitions of the derived types
    use m_global_parameters
    use m_helper_basic
    use m_helper
    use m_mpi_common

    implicit none

    private; public :: s_apply_ib_patches, s_update_ib_rotation_matrix, s_instantiate_STL_models, s_decode_patch_periodicity, &
        & s_initialize_ib_airfoils

contains

    !> Apply all immersed boundary patch geometries to mark interior cells in the IB marker array
    impure subroutine s_apply_ib_patches(ib_markers)

        type(integer_field), intent(inout) :: ib_markers

        if (many_ib_patch_parallelism) then
            call s_apply_ib_patches_ib_parallelism(ib_markers)
        else
            call s_apply_ib_patches_grid_cell_parallelism(ib_markers)
        end if

    end subroutine s_apply_ib_patches

    subroutine s_apply_ib_patches_grid_cell_parallelism(ib_markers)

        type(integer_field), intent(inout) :: ib_markers
        integer :: patch_id, airfoil_id, model_id, encoded_patch_id, i, j, k, il, ir, jl, jr, kl, kr, xp, yp, zp  !< iterators
        integer :: xp_lower, xp_upper, yp_lower, yp_upper, zp_lower, zp_upper                                     !< periodic bounds
        real(wp), dimension(3) :: center, xyz_local, length
        real(wp) :: radius

        !  3D Patch Geometries

        if (num_dims == 3) then
            call s_get_periodicities(xp_lower, xp_upper, yp_lower, yp_upper, zp_lower, zp_upper)
            do xp = xp_lower, xp_upper
                do yp = yp_lower, yp_upper
                    do zp = zp_lower, zp_upper
                        do patch_id = 1, num_ibs
                            center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
                            center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
                            center(3) = patch_ib(patch_id)%z_centroid + real(zp, wp)*(z_domain%end - z_domain%beg)

                            ! encode the periodicity information into the patch_id
                            call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, zp, encoded_patch_id)

                            ! find the indices to the left and right of the IB in i, j, k
                            call get_bounding_indices(center, patch_id, il, ir, jl, jr, kl, kr)

                            $:GPU_PARALLEL_LOOP(private='[i, j, k, xyz_local, length, radius, airfoil_id]', copyin='[patch_id, &
                                                & encoded_patch_id, center]', collapse=3)
                            do k = kl, kr
                                do j = jl, jr
                                    do i = il, ir
                                        ! get coordinate frame centered on IB
                                        xyz_local = [x_cc(i) - center(1), y_cc(j) - center(2), z_cc(k) - center(3)]
                                        ! rotate the frame into the IB's coordinates
                                        xyz_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse, xyz_local)

                                        ! perform the interior check for the patch geometry of this IB
                                        if (patch_ib(patch_id)%geometry == 8) then
                                            ! sphere geometry
                                            radius = patch_ib(patch_id)%radius
                                            if (f_is_inside_sphere(xyz_local(1), xyz_local(2), xyz_local(3), &
                                                & radius)) ib_markers%sf(i, j, k) = encoded_patch_id
                                        else if (patch_ib(patch_id)%geometry == 9) then
                                            ! cuboid geometry
                                            length = [patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y, &
                                                               & patch_ib(patch_id)%length_z]
                                            if (f_is_inside_cuboid(xyz_local(1), xyz_local(2), xyz_local(3), &
                                                & length)) ib_markers%sf(i, j, k) = encoded_patch_id
                                        else if (patch_ib(patch_id)%geometry == 10) then
                                            ! cylinder geometry
                                            radius = patch_ib(patch_id)%radius
                                            if (f_is_inside_cylinder(xyz_local(2), xyz_local(3), xyz_local(1), radius, &
                                                & patch_ib(patch_id)%length_x)) ib_markers%sf(i, j, k) = encoded_patch_id
                                        else if (patch_ib(patch_id)%geometry == 11) then
                                            ! 3D airfoil geometry
                                            airfoil_id = patch_ib(patch_id)%airfoil_id
                                            xyz_local = xyz_local - patch_ib(patch_id)%centroid_offset
                                            if (f_is_inside_airfoil(xyz_local(1), xyz_local(2), xyz_local(3), &
                                                & patch_ib(patch_id)%length_z, airfoil_id)) ib_markers%sf(i, j, &
                                                & k) = encoded_patch_id
                                        else if (patch_ib(patch_id)%geometry == 12) then
                                            ! STL model geometry
                                            xyz_local = xyz_local - patch_ib(patch_id)%centroid_offset
                                            ! if (f_is_inside_model(patch_id, x_cc(i), y_cc(j), z_cc(k))) ib_markers%sf(i, j, &
                                            !     & k) = encoded_patch_id
                                        end if
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end do
                    end do
                end do
            end do

            ! 2D Patch Geometries
        else if (num_dims == 2) then
            call s_get_periodicities(xp_lower, xp_upper, yp_lower, yp_upper)
            do xp = xp_lower, xp_upper
                do yp = yp_lower, yp_upper
                    do patch_id = 1, num_ibs
                        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
                        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
                        center(3) = 0._wp

                        ! encode the periodicity information into the patch_id
                        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, 0, encoded_patch_id)

                        ! find the indices to the left and right of the IB in i, j, k
                        call get_bounding_indices(center, patch_id, il, ir, jl, jr, kl, kr)

                        $:GPU_PARALLEL_LOOP(private='[i, j, xyz_local, airfoil_id]', copyin='[patch_id, encoded_patch_id, &
                                            & center]', collapse=2)
                        do j = jl, jr
                            do i = il, ir
                                ! get coordinate frame centered on IB
                                xyz_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
                                ! rotate the frame into the IB's coordinates
                                xyz_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse, xyz_local)

                                ! perform the interior check for the patch geometry of this IB
                                if (patch_ib(patch_id)%geometry == 2) then
                                    ! circular geometries
                                    radius = patch_ib(patch_id)%radius
                                    if (f_is_inside_cylinder(xyz_local(1), xyz_local(2), 0._wp, radius, 0._wp)) ib_markers%sf(i, &
                                        & j, k) = encoded_patch_id
                                else if (patch_ib(patch_id)%geometry == 3) then
                                    ! rectangular geometries
                                    length = [patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y, 0._wp]
                                    if (f_is_inside_cuboid(xyz_local(1), xyz_local(2), xyz_local(3), length)) ib_markers%sf(i, j, &
                                        & 0) = encoded_patch_id
                                else if (patch_ib(patch_id)%geometry == 4) then
                                    ! 2D airfoil geometry
                                    airfoil_id = patch_ib(patch_id)%airfoil_id
                                    xyz_local = xyz_local - patch_ib(patch_id)%centroid_offset
                                    if (f_is_inside_airfoil(xyz_local(1), xyz_local(2), 0._wp, 0._wp, &
                                        & airfoil_id)) ib_markers%sf(i, j, 0) = encoded_patch_id
                                else if (patch_ib(patch_id)%geometry == 5) then
                                    ! STL model geometry
                                    xyz_local = xyz_local - patch_ib(patch_id)%centroid_offset
                                    ! if (f_is_inside_model(patch_id, x_cc(i), y_cc(j), 0._wp)) ib_markers%sf(i, j, &
                                    !     & 0) = encoded_patch_id
                                else if (patch_ib(patch_id)%geometry == 6) then
                                    ! ellipse geometry
                                    length = [patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y, 0._wp]
                                    if (f_is_inside_ellipse(xyz_local(1), xyz_local(2), length)) ib_markers%sf(i, j, &
                                        & 0) = encoded_patch_id
                                end if
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end do
                end do
            end do
        end if

    end subroutine s_apply_ib_patches_grid_cell_parallelism

    subroutine s_apply_ib_patches_ib_parallelism(ib_markers)

        type(integer_field), intent(inout) :: ib_markers
        integer :: patch_id, airfoil_id, model_id, encoded_patch_id, i, j, k, il, ir, jl, jr, kl, kr, xp, yp, zp  !< iterators
        integer :: xp_lower, xp_upper, yp_lower, yp_upper, zp_lower, zp_upper                                     !< periodic bounds
        real(wp), dimension(3) :: center, xyz_local, length
        real(wp) :: radius

        if (num_dims == 3) then
            ! get the periodicities
            call s_get_periodicities(xp_lower, xp_upper, yp_lower, yp_upper, zp_lower, zp_upper)

            do xp = xp_lower, xp_upper
                do yp = yp_lower, yp_upper
                    do zp = zp_lower, zp_upper
                        $:GPU_PARALLEL_LOOP(private='[xp, yp, zp, i, il, ir, j, jl, jr, k, kl, kr, xyz_local, length, radius, &
                                            & patch_id, airfoil_id, model_id, encoded_patch_id, center]')
                        do patch_id = 1, num_ibs
                            center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
                            center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
                            center(3) = patch_ib(patch_id)%z_centroid + real(zp, wp)*(z_domain%end - z_domain%beg)

                            ! encode the periodicity information into the patch_id
                            call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, zp, encoded_patch_id)

                            ! find the indices to the left and right of the IB in i, j, k
                            call get_bounding_indices(center, patch_id, il, ir, jl, jr, kl, kr)

                            do k = kl, kr
                                do j = jl, jr
                                    do i = il, ir
                                        ! get coordinate frame centered on IB
                                        xyz_local = [x_cc(i) - center(1), y_cc(j) - center(2), z_cc(k) - center(3)]
                                        ! rotate the frame into the IB's coordinates
                                        xyz_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse, xyz_local)

                                        ! perform the interior check for the patch geometry of this IB
                                        if (patch_ib(patch_id)%geometry == 8) then
                                            ! sphere geometry
                                            radius = patch_ib(patch_id)%radius
                                            if (f_is_inside_sphere(xyz_local(1), xyz_local(2), xyz_local(3), &
                                                & radius)) ib_markers%sf(i, j, k) = encoded_patch_id
                                        else if (patch_ib(patch_id)%geometry == 9) then
                                            ! cuboid geometry
                                            length = [patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y, &
                                                               & patch_ib(patch_id)%length_z]
                                            if (f_is_inside_cuboid(xyz_local(1), xyz_local(2), xyz_local(3), &
                                                & length)) ib_markers%sf(i, j, k) = encoded_patch_id
                                        else if (patch_ib(patch_id)%geometry == 10) then
                                            ! cylinder geometry
                                            radius = patch_ib(patch_id)%radius
                                            if (f_is_inside_cylinder(xyz_local(2), xyz_local(3), xyz_local(1), radius, &
                                                & patch_ib(patch_id)%length_x)) ib_markers%sf(i, j, k) = encoded_patch_id
                                        else if (patch_ib(patch_id)%geometry == 11) then
                                            ! 3D airfoil geometry
                                            airfoil_id = patch_ib(patch_id)%airfoil_id
                                            xyz_local = xyz_local - patch_ib(patch_id)%centroid_offset
                                            if (f_is_inside_airfoil(xyz_local(1), xyz_local(2), xyz_local(3), &
                                                & patch_ib(patch_id)%length_z, airfoil_id)) ib_markers%sf(i, j, &
                                                & k) = encoded_patch_id
                                        else if (patch_ib(patch_id)%geometry == 12) then
                                            ! STL model geometry
                                            xyz_local = xyz_local - patch_ib(patch_id)%centroid_offset
                                            ! if (f_is_inside_model(patch_id, x_cc(i), y_cc(j), z_cc(k))) ib_markers%sf(i, j, &
                                            !     & k) = encoded_patch_id
                                        end if
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end do
                end do
            end do
        else if (num_dims == 2) then
            ! get the periodicities
            call s_get_periodicities(xp_lower, xp_upper, yp_lower, yp_upper)

            $:GPU_PARALLEL_LOOP(private='[xp, yp, patch_id, center]', collapse=3)
            do xp = xp_lower, xp_upper
                do yp = yp_lower, yp_upper
                    $:GPU_PARALLEL_LOOP(private='[xp, yp, i, il, ir, j, jl, jr, xyz_local, length, radius, &
                                        & bounding_box_corner_distance, patch_id, airfoil_id, model_id, encoded_patch_id, center]')
                    do patch_id = 1, num_ibs
                        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
                        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
                        center(3) = 0._wp

                        ! encode the periodicity information into the patch_id
                        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, 0, encoded_patch_id)

                        ! find the indices to the left and right of the IB in i, j, k
                        ! find the indices to the left and right of the IB in i, j, k
                        call get_bounding_indices(center, patch_id, il, ir, jl, jr, kl, kr)

                        do j = jl, jr
                            do i = il, ir
                                ! get coordinate frame centered on IB
                                xyz_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
                                ! rotate the frame into the IB's coordinates
                                xyz_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse, xyz_local)

                                ! perform the interior check for the patch geometry of this IB
                                if (patch_ib(patch_id)%geometry == 2) then
                                    ! circular geometries
                                    radius = patch_ib(patch_id)%radius
                                    if (f_is_inside_cylinder(xyz_local(1), xyz_local(2), 0._wp, radius, 0._wp)) ib_markers%sf(i, &
                                        & j, k) = encoded_patch_id
                                else if (patch_ib(patch_id)%geometry == 3) then
                                    ! rectangular geometries
                                    length = [patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y, 0._wp]
                                    if (f_is_inside_cuboid(xyz_local(1), xyz_local(2), xyz_local(3), length)) ib_markers%sf(i, j, &
                                        & k) = encoded_patch_id
                                else if (patch_ib(patch_id)%geometry == 4) then
                                    ! 2D airfoil geometry
                                    airfoil_id = patch_ib(patch_id)%airfoil_id
                                    xyz_local = xyz_local - patch_ib(patch_id)%centroid_offset
                                    if (f_is_inside_airfoil(xyz_local(1), xyz_local(2), 0._wp, 0._wp, &
                                        & airfoil_id)) ib_markers%sf(i, j, k) = encoded_patch_id
                                else if (patch_ib(patch_id)%geometry == 5) then
                                    ! STL model geometry
                                    model_id = patch_ib(patch_id)%model_id
                                    xyz_local = xyz_local - patch_ib(patch_id)%centroid_offset
                                    ! if (f_is_inside_model(patch_id, x_cc(i), y_cc(j), 0._wp)) ib_markers%sf(i, j, &
                                    !     & k) = encoded_patch_id
                                else if (patch_ib(patch_id)%geometry == 6) then
                                    ! ellipse geometry
                                    length = [patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y, 0._wp]
                                    if (f_is_inside_ellipse(xyz_local(1), xyz_local(2), length)) ib_markers%sf(i, j, &
                                        & k) = encoded_patch_id
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end do
        end if

    end subroutine s_apply_ib_patches_ib_parallelism

    !> The STL patch is a 2D geometry that is imported from an STL file.
    subroutine s_ib_model(patch_id, ib_markers, xp, yp)

        integer, intent(in)                :: patch_id
        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in)                :: xp, yp                !< integers containing the periodicity projection information
        integer                            :: i, j, il, ir, jl, jr  !< Generic loop iterators
        integer                            :: model_id, encoded_patch_id
        integer                            :: cx, cy
        real(wp)                           :: lx(2), ly(2)
        real(wp), dimension(1:2)           :: bbox_min, bbox_max
        real(wp), dimension(1:3)           :: local_corner, world_corner
        real(wp)                           :: eta, threshold
        real(wp), dimension(1:3)           :: offset
        real(wp), dimension(1:3)           :: center, xy_local
        real(wp), dimension(1:3,1:3)       :: inverse_rotation, rotation

        center = 0._wp
        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
        inverse_rotation(:,:) = patch_ib(patch_id)%rotation_matrix_inverse(:,:)
        rotation(:,:) = patch_ib(patch_id)%rotation_matrix(:,:)
        offset(:) = patch_ib(patch_id)%centroid_offset(:)
        model_id = patch_ib(patch_id)%model_id
        threshold = stl_models(model_id)%model_threshold

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, 0, encoded_patch_id)

        il = -gp_layers - 1
        jl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        $:GPU_PARALLEL_LOOP(private='[i, j, xy_local, eta]', copyin='[patch_id, model_id, encoded_patch_id, center, &
                            & inverse_rotation, offset, threshold]', collapse=2)
        do i = il, ir
            do j = jl, jr
                xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
                xy_local = matmul(inverse_rotation, xy_local)
                xy_local = xy_local - offset

                eta = f_model_is_inside_flat(gpu_ntrs(model_id), model_id, xy_local)

                ! Reading STL boundary vertices and compute the levelset and levelset_norm
                if (eta > threshold) then
                    ib_markers%sf(i, j, 0) = encoded_patch_id
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_model

    !> The STL patch is a 3D geometry that is imported from an STL file.
    subroutine s_ib_3d_model(patch_id, ib_markers, xp, yp, zp)

        integer, intent(in)                :: patch_id
        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in)                :: xp, yp, zp  !< integers containing the periodicity projection information
        integer                            :: i, j, k, il, ir, jl, jr, kl, kr  !< Generic loop iterators
        integer                            :: model_id, encoded_patch_id
        real(wp)                           :: eta, threshold
        real(wp), dimension(1:3)           :: offset
        real(wp), dimension(1:3)           :: center, xyz_local
        real(wp), dimension(1:3,1:3)       :: inverse_rotation, rotation
        integer                            :: cx, cy, cz
        real(wp)                           :: lx(2), ly(2), lz(2)
        real(wp), dimension(1:3)           :: bbox_min, bbox_max, local_corner, world_corner

        center = 0._wp
        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
        center(3) = patch_ib(patch_id)%z_centroid + real(zp, wp)*(z_domain%end - z_domain%beg)
        inverse_rotation(:,:) = patch_ib(patch_id)%rotation_matrix_inverse(:,:)
        offset(:) = patch_ib(patch_id)%centroid_offset(:)
        model_id = patch_ib(patch_id)%model_id
        threshold = stl_models(model_id)%model_threshold
        rotation(:,:) = patch_ib(patch_id)%rotation_matrix(:,:)

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, zp, encoded_patch_id)

        il = -gp_layers - 1
        jl = -gp_layers - 1
        kl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        kr = p + gp_layers + 1

        $:GPU_PARALLEL_LOOP(private='[i, j, k, xyz_local, eta]', copyin='[patch_id, model_id, encoded_patch_id, center, &
                            & inverse_rotation, offset, threshold]', collapse=3)
        do i = il, ir
            do j = jl, jr
                do k = kl, kr
                    xyz_local = [x_cc(i) - center(1), y_cc(j) - center(2), z_cc(k) - center(3)]
                    xyz_local = matmul(inverse_rotation, xyz_local)
                    xyz_local = xyz_local - offset

                    eta = f_model_is_inside_flat(gpu_ntrs(model_id), model_id, xyz_local)

                    if (eta > threshold) then
                        ib_markers%sf(i, j, k) = encoded_patch_id
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_3d_model

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

    !> Compute a rotation matrix for converting to the rotating frame of the boundary
    subroutine s_update_ib_rotation_matrix(patch_id)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)          :: patch_id
        real(wp), dimension(3, 3, 3) :: rotation
        real(wp)                     :: angle

        ! construct the x, y, and z rotation matrices

        if (num_dims == 3) then
            ! also compute the x and y axes in 3D
            angle = patch_ib(patch_id)%angles(1)
            rotation(1, 1,:) = [1._wp, 0._wp, 0._wp]
            rotation(1, 2,:) = [0._wp, cos(angle), -sin(angle)]
            rotation(1, 3,:) = [0._wp, sin(angle), cos(angle)]

            angle = patch_ib(patch_id)%angles(2)
            rotation(2, 1,:) = [cos(angle), 0._wp, sin(angle)]
            rotation(2, 2,:) = [0._wp, 1._wp, 0._wp]
            rotation(2, 3,:) = [-sin(angle), 0._wp, cos(angle)]

            ! apply the y rotation to the x rotation
            patch_ib(patch_id)%rotation_matrix(:,:) = matmul(rotation(1,:,:), rotation(2,:,:))
            patch_ib(patch_id)%rotation_matrix_inverse(:,:) = matmul(transpose(rotation(2,:,:)), transpose(rotation(1,:,:)))
        end if

        ! z component first, since it applies in 2D and 3D
        angle = patch_ib(patch_id)%angles(3)
        rotation(3, 1,:) = [cos(angle), -sin(angle), 0._wp]
        rotation(3, 2,:) = [sin(angle), cos(angle), 0._wp]
        rotation(3, 3,:) = [0._wp, 0._wp, 1._wp]

        if (num_dims == 3) then
            ! apply the z rotation to the xy rotation in 3D
            patch_ib(patch_id)%rotation_matrix(:,:) = matmul(patch_ib(patch_id)%rotation_matrix(:,:), rotation(3,:,:))
            patch_ib(patch_id)%rotation_matrix_inverse(:,:) = matmul(transpose(rotation(3,:,:)), &
                     & patch_ib(patch_id)%rotation_matrix_inverse(:,:))
        else
            ! write out only the z rotation in 2D
            patch_ib(patch_id)%rotation_matrix(:,:) = rotation(3,:,:)
            patch_ib(patch_id)%rotation_matrix_inverse(:,:) = transpose(rotation(3,:,:))
        end if

    end subroutine s_update_ib_rotation_matrix

    subroutine get_bounding_indices(center, patch_id, il, ir, jl, jr, kl, kr)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: center
        integer, intent(in)                :: patch_id
        integer, intent(out)               :: il, ir, jl, jr, kl, kr
        real(wp), dimension(3)             :: bbox_min, bbox_max, local_corner, world_corner
        real(wp), dimension(2)             :: lx, ly, lz
        integer                            :: cx, cy, cz

        if (patch_ib(patch_id)%geometry == 2 .or. patch_ib(patch_id)%geometry == 8) then
            ! circle and sphere geometries
            bbox_min = center - patch_ib(patch_id)%radius
            bbox_max = center + patch_ib(patch_id)%radius
        else if (patch_ib(patch_id)%geometry == 3) then
            ! rectangular geometries
            bbox_min = center - 0.5_wp*sqrt(patch_ib(patch_id)%length_x**2 + patch_ib(patch_id)%length_y**2)
            bbox_max = center + 0.5_wp*sqrt(patch_ib(patch_id)%length_x**2 + patch_ib(patch_id)%length_y**2)
        else if (patch_ib(patch_id)%geometry == 4 .or. patch_ib(patch_id)%geometry == 11) then
            ! airfoil geometries TODO :: This can be better optimized, since airfoils are typically very long in one dimension
            bbox_min = center - ib_airfoil(patch_ib(patch_id)%airfoil_id)%c
            bbox_max = center + ib_airfoil(patch_ib(patch_id)%airfoil_id)%c
        else if (patch_ib(patch_id)%geometry == 5) then
            ! STL model geometry
            lx(1) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 1, 1) + patch_ib(patch_id)%centroid_offset(1)
            lx(2) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 1, 3) + patch_ib(patch_id)%centroid_offset(1)
            ly(1) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 2, 1) + patch_ib(patch_id)%centroid_offset(2)
            ly(2) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 2, 3) + patch_ib(patch_id)%centroid_offset(2)

            bbox_min = 1e12
            bbox_max = -1e12
            ! Enumerate all 4 corners of the local bounding box, rotate to world space, track world-space AABB
            do cx = 1, 2
                do cy = 1, 2
                    local_corner = [lx(cx), ly(cy), 0._wp]
                    world_corner = matmul(patch_ib(patch_id)%rotation_matrix, local_corner) + center
                    bbox_min(1) = min(bbox_min(1), world_corner(1))
                    bbox_min(2) = min(bbox_min(2), world_corner(2))
                    bbox_max(1) = max(bbox_max(1), world_corner(1))
                    bbox_max(2) = max(bbox_max(2), world_corner(2))
                end do
            end do
        else if (patch_ib(patch_id)%geometry == 6) then
            ! ellipse geometry
            bbox_min = center - 0.5_wp*max(patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y)
            bbox_max = center + 0.5_wp*max(patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y)
        else if (patch_ib(patch_id)%geometry == 9) then
            ! cuboid geometries
            bbox_min = center - 0.5_wp*sqrt(patch_ib(patch_id)%length_x**2 + patch_ib(patch_id)%length_y**2 &
                                            & + patch_ib(patch_id)%length_z**2)
            bbox_max = center + 0.5_wp*sqrt(patch_ib(patch_id)%length_x**2 + patch_ib(patch_id)%length_y**2 &
                                            & + patch_ib(patch_id)%length_z**2)
        else if (patch_ib(patch_id)%geometry == 10) then
            ! cylinder geometry
            bbox_min = center - sqrt(patch_ib(patch_id)%radius**2 + patch_ib(patch_id)%length_x**2)
            bbox_max = center + sqrt(patch_ib(patch_id)%radius**2 + patch_ib(patch_id)%length_x**2)
        else if (patch_ib(patch_id)%geometry == 12) then
            ! Local-space bounding box extents (min=1, max=2 in the third index)
            lx(1) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 1, 1) + patch_ib(patch_id)%centroid_offset(1)
            lx(2) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 1, 3) + patch_ib(patch_id)%centroid_offset(1)
            ly(1) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 2, 1) + patch_ib(patch_id)%centroid_offset(2)
            ly(2) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 2, 3) + patch_ib(patch_id)%centroid_offset(2)
            lz(1) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 3, 1) + patch_ib(patch_id)%centroid_offset(3)
            lz(2) = stl_bounding_boxes(patch_ib(patch_id)%model_id, 3, 3) + patch_ib(patch_id)%centroid_offset(3)

            bbox_min = 1e12
            bbox_max = -1e12
            ! Enumerate all 8 corners of the local bounding box, rotate to world space, track world-space AABB
            do cx = 1, 2
                do cy = 1, 2
                    do cz = 1, 2
                        local_corner = [lx(cx), ly(cy), lz(cz)]
                        world_corner = matmul(patch_ib(patch_id)%rotation_matrix, local_corner) + center
                        bbox_min(1) = min(bbox_min(1), world_corner(1))
                        bbox_min(2) = min(bbox_min(2), world_corner(2))
                        bbox_min(3) = min(bbox_min(3), world_corner(3))
                        bbox_max(1) = max(bbox_max(1), world_corner(1))
                        bbox_max(2) = max(bbox_max(2), world_corner(2))
                        bbox_max(3) = max(bbox_max(3), world_corner(3))
                    end do
                end do
            end do
        end if

        il = -gp_layers - 1
        jl = -gp_layers - 1
        kl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        kr = p + gp_layers + 1
        call get_indices_from_bounds(bbox_min(1), bbox_max(1), x_cc, il, ir)
        call get_indices_from_bounds(bbox_min(2), bbox_max(2), y_cc, jl, jr)
        if (num_dims == 3) call get_indices_from_bounds(bbox_min(3), bbox_max(3), z_cc, kl, kr)

    end subroutine get_bounding_indices

    subroutine get_indices_from_bounds(left_bound, right_bound, cell_centers, left_index, right_index)

        real(wp), intent(in)                         :: left_bound, right_bound
        integer, intent(inout)                       :: left_index, right_index
        real(wp), dimension(-buff_size:), intent(in) :: cell_centers
        integer                                      :: itr_left, itr_middle, itr_right

        itr_left = left_index
        itr_right = right_index

        do while (itr_left + 1 < itr_right)
            itr_middle = (itr_left + itr_right)/2
            if (cell_centers(itr_middle) < left_bound) then
                itr_left = itr_middle
            else if (cell_centers(itr_middle) > left_bound) then
                itr_right = itr_middle
            else
                itr_left = itr_middle
                exit
            end if
        end do
        left_index = itr_left

        itr_right = right_index
        do while (itr_left + 1 < itr_right)
            itr_middle = (itr_left + itr_right)/2
            if (cell_centers(itr_middle) < right_bound) then
                itr_left = itr_middle
            else if (cell_centers(itr_middle) > right_bound) then
                itr_right = itr_middle
            else
                itr_right = itr_middle
                exit
            end if
        end do
        right_index = itr_right

    end subroutine get_indices_from_bounds

    !> Encode the patch ID with a unique offset containing periodicity information
    subroutine s_encode_patch_periodicity(patch_id, x_periodicity, y_periodicity, z_periodicity, encoded_patch_id)

        integer, intent(in)  :: patch_id, x_periodicity, y_periodicity, z_periodicity
        integer, intent(out) :: encoded_patch_id
        integer              :: temp_x_per, temp_y_per, temp_z_per, offset

        encoded_patch_id = patch_id

        temp_x_per = x_periodicity; if (x_periodicity == -1) temp_x_per = 2
        temp_y_per = y_periodicity; if (y_periodicity == -1) temp_y_per = 2
        temp_z_per = z_periodicity; if (z_periodicity == -1) temp_z_per = 2

        offset = (num_gbl_ibs + 1)*temp_x_per + 3*(num_gbl_ibs + 1)*temp_y_per + 9*(num_gbl_ibs + 1)*temp_z_per
        encoded_patch_id = patch_id + offset

    end subroutine s_encode_patch_periodicity

    !> Decode the encoded ID to recover the original patch ID and periodicity
    subroutine s_decode_patch_periodicity(encoded_patch_id, patch_id, x_periodicity, y_periodicity, z_periodicity)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)            :: encoded_patch_id
        integer, intent(out)           :: patch_id
        integer, intent(out), optional :: x_periodicity, y_periodicity, z_periodicity
        integer                        :: offset, remainder, xp, yp, zp, base

        base = num_gbl_ibs + 1

        patch_id = mod(encoded_patch_id - 1, base) + 1
        offset = (encoded_patch_id - patch_id)/base

        xp = mod(offset, 3)
        remainder = offset/3
        yp = mod(remainder, 3)
        zp = remainder/3

        ! Reverse map: 2 -> -1, 0 -> 0, 1 -> 1
        if (present(x_periodicity) .and. present(y_periodicity) .and. present(z_periodicity)) then
            x_periodicity = xp; if (xp == 2) x_periodicity = -1
            y_periodicity = yp; if (yp == 2) y_periodicity = -1
            z_periodicity = zp; if (zp == 2) z_periodicity = -1
        end if

    end subroutine s_decode_patch_periodicity

    !> Determine the periodic wrapping bounds in each direction
    subroutine s_get_periodicities(xp_lower, xp_upper, yp_lower, yp_upper, zp_lower, zp_upper)

        integer, intent(out)           :: xp_lower, xp_upper, yp_lower, yp_upper
        integer, intent(out), optional :: zp_lower, zp_upper

        ! check domain wraps in x, y

        #:for X in [('x'), ('y')]
            ! check for periodicity
            if (ib_bc_${X}$%beg == BC_PERIODIC) then
                ${X}$p_lower = -1
                ${X}$p_upper = 1
            else
                ! if it is not periodic, then both elements are 0
                ${X}$p_lower = 0
                ${X}$p_upper = 0
            end if
        #:endfor

        ! z only if 3D
        if (present(zp_lower) .and. p /= 0) then
            if (ib_bc_z%beg == BC_PERIODIC) then
                zp_lower = -1
                zp_upper = 1
            else
                zp_lower = 0
                zp_upper = 0
            end if
        end if

    end subroutine s_get_periodicities

    !> Archimedes spiral function
    pure elemental function f_r(myth, offset, a)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: myth, offset, a
        real(wp)             :: b
        real(wp)             :: f_r

        ! r(th) = a + b*th

        b = 2._wp*a/(2._wp*pi)
        f_r = a + b*myth + offset

    end function f_r

end module m_ib_patches
