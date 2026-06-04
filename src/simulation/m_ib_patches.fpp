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

    !> Apply all immersed boundary patch geometries to mark interior cells in the IB marker array
    impure subroutine s_apply_ib_patches(ib_markers)

        type(integer_field), intent(inout) :: ib_markers
        integer                            :: i, xp, yp, zp                                               !< iterators
        integer                            :: xp_lower, xp_upper, yp_lower, yp_upper, zp_lower, zp_upper  !< periodic bounds

        !  3D Patch Geometries

        if (p > 0) then
            !> IB Patches
            !> @{
            call s_get_periodicities(xp_lower, xp_upper, yp_lower, yp_upper, zp_lower, zp_upper)
            do xp = xp_lower, xp_upper
                do yp = yp_lower, yp_upper
                    do zp = zp_lower, zp_upper
                        do i = 1, num_ibs
                            if (patch_ib(i)%geometry == 8) then
                                call s_ib_sphere(i, ib_markers, xp, yp, zp)
                            else if (patch_ib(i)%geometry == 9) then
                                call s_ib_cuboid(i, ib_markers, xp, yp, zp)
                            else if (patch_ib(i)%geometry == 10) then
                                call s_ib_cylinder(i, ib_markers, xp, yp, zp)
                            else if (patch_ib(i)%geometry == 11) then
                                call s_ib_3D_airfoil(i, ib_markers, xp, yp, zp)
                            else if (patch_ib(i)%geometry == 12) then
                                call s_ib_3d_model(i, ib_markers, xp, yp, zp)
                            end if
                        end do
                    end do
                end do
            end do
            !> @}

            ! 2D Patch Geometries
        else if (n > 0) then
            !> IB Patches
            !> @{
            call s_get_periodicities(xp_lower, xp_upper, yp_lower, yp_upper)
            do xp = xp_lower, xp_upper
                do yp = yp_lower, yp_upper
                    do i = 1, num_ibs
                        if (patch_ib(i)%geometry == 2) then
                            call s_ib_circle(i, ib_markers, xp, yp)
                        else if (patch_ib(i)%geometry == 3) then
                            call s_ib_rectangle(i, ib_markers, xp, yp)
                        else if (patch_ib(i)%geometry == 4) then
                            call s_ib_airfoil(i, ib_markers, xp, yp)
                        else if (patch_ib(i)%geometry == 5) then
                            call s_ib_model(i, ib_markers, xp, yp)
                        else if (patch_ib(i)%geometry == 6) then
                            call s_ib_ellipse(i, ib_markers, xp, yp)
                        end if
                    end do
                end do
            end do
            !> @}
        end if

    end subroutine s_apply_ib_patches

    !> Mark cells inside a circular immersed boundary
    subroutine s_ib_circle(patch_id, ib_markers, xp, yp)

        integer, intent(in)                :: patch_id
        integer, intent(in)                :: xp, yp                !< integers containing the periodicity projection information
        type(integer_field), intent(inout) :: ib_markers
        real(wp), dimension(1:2)           :: center
        real(wp)                           :: radius
        integer                            :: i, j, il, ir, jl, jr  !< Generic loop iterators
        integer                            :: encoded_patch_id

        ! Transferring the circular patch's radius, centroid, smearing patch identity and smearing coefficient information

        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
        radius = patch_ib(patch_id)%radius

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, 0, encoded_patch_id)

        ! find the indices to the left and right of the IB in i, j, k
        il = -gp_layers - 1
        jl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        call get_bounding_indices(center(1) - radius, center(1) + radius, x_cc, il, ir)
        call get_bounding_indices(center(2) - radius, center(2) + radius, y_cc, jl, jr)

        ! Assign primitive variables if circle covers cell and patch has write permission

        $:GPU_PARALLEL_LOOP(private='[i, j]', copyin='[encoded_patch_id, center, radius]', collapse=2)
        do j = jl, jr
            do i = il, ir
                if ((x_cc(i) - center(1))**2 + (y_cc(j) - center(2))**2 <= radius**2) then
                    ib_markers%sf(i, j, 0) = encoded_patch_id
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_circle

    !> Mark cells inside a 2D NACA 4-digit airfoil immersed boundary
    subroutine s_ib_airfoil(patch_id, ib_markers, xp, yp)

        integer, intent(in)                :: patch_id
        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in)                :: xp, yp            !< integers containing the periodicity projection information
        real(wp)                           :: f, ca_in
        integer                            :: i, j, k, il, ir, jl, jr
        integer                            :: Np_local, airfoil_id
        integer                            :: encoded_patch_id
        real(wp), dimension(1:3)           :: xy_local, offset  !< x and y coordinates in local IB frame
        real(wp), dimension(1:2)           :: center            !< x and y coordinates in local IB frame

        airfoil_id = patch_ib(patch_id)%airfoil_id
        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
        ca_in = ib_airfoil(airfoil_id)%c
        Np_local = ib_airfoil_grids(airfoil_id)%Np
        offset(:) = patch_ib(patch_id)%centroid_offset(:)

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, 0, encoded_patch_id)

        ! find the indices to the left and right of the IB in i, j, k
        il = -gp_layers - 1
        jl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        ! maximum distance any marker can be from the center is the length of the airfoil
        call get_bounding_indices(center(1) - ca_in, center(1) + ca_in, x_cc, il, ir)
        call get_bounding_indices(center(2) - ca_in, center(2) + ca_in, y_cc, jl, jr)

        $:GPU_PARALLEL_LOOP(private='[i, j, xy_local, k, f]', copyin='[encoded_patch_id, center, offset, ca_in, airfoil_id, &
                            & Np_local, ib_airfoil_grids(airfoil_id)%upper, ib_airfoil_grids(airfoil_id)%lower]', collapse=2)
        do j = jl, jr
            do i = il, ir
                xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]  ! get coordinate frame centered on IB
                xy_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse, xy_local)  ! rotate the frame into the IB's coordinates
                xy_local = xy_local - offset  ! airfoils are a patch that require a centroid offset

                if (xy_local(1) >= 0._wp .and. xy_local(1) <= ca_in) then
                    if (xy_local(2) >= 0._wp) then
                        k = 1
                        do while (ib_airfoil_grids(airfoil_id)%upper(k)%x < xy_local(1) .and. k <= Np_local)
                            k = k + 1
                        end do
                        if (f_approx_equal(ib_airfoil_grids(airfoil_id)%upper(k)%x, xy_local(1))) then
                            if (xy_local(2) <= ib_airfoil_grids(airfoil_id)%upper(k)%y) then
                                ib_markers%sf(i, j, 0) = encoded_patch_id
                            end if
                        else
                            f = (ib_airfoil_grids(airfoil_id)%upper(k)%x - xy_local(1))/(ib_airfoil_grids(airfoil_id)%upper(k)%x &
                                 & - ib_airfoil_grids(airfoil_id)%upper(k - 1)%x)
                            if (xy_local(2) <= ((1._wp - f)*ib_airfoil_grids(airfoil_id)%upper(k)%y &
                                & + f*ib_airfoil_grids(airfoil_id)%upper(k - 1)%y)) then
                                ib_markers%sf(i, j, 0) = encoded_patch_id
                            end if
                        end if
                    else
                        k = 1
                        do while (ib_airfoil_grids(airfoil_id)%lower(k)%x < xy_local(1))
                            k = k + 1
                        end do
                        if (f_approx_equal(ib_airfoil_grids(airfoil_id)%lower(k)%x, xy_local(1))) then
                            if (xy_local(2) >= ib_airfoil_grids(airfoil_id)%lower(k)%y) then
                                ib_markers%sf(i, j, 0) = encoded_patch_id
                            end if
                        else
                            f = (ib_airfoil_grids(airfoil_id)%lower(k)%x - xy_local(1))/(ib_airfoil_grids(airfoil_id)%lower(k)%x &
                                 & - ib_airfoil_grids(airfoil_id)%lower(k - 1)%x)
                            if (xy_local(2) >= ((1._wp - f)*ib_airfoil_grids(airfoil_id)%lower(k)%y &
                                & + f*ib_airfoil_grids(airfoil_id)%lower(k - 1)%y)) then
                                ib_markers%sf(i, j, 0) = encoded_patch_id
                            end if
                        end if
                    end if
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_airfoil

    !> Mark cells inside a 3D extruded NACA 4-digit airfoil immersed boundary with finite span
    subroutine s_ib_3D_airfoil(patch_id, ib_markers, xp, yp, zp)

        integer, intent(in)                :: patch_id
        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in)                :: xp, yp, zp  !< integers containing the periodicity projection information
        real(wp)                           :: lz, z_max, z_min, f, ca_in
        integer                            :: i, j, k, l, il, ir, jl, jr, ll, lr
        integer                            :: airfoil_id
        integer                            :: encoded_patch_id
        real(wp), dimension(1:3)           :: xyz_local, center, offset  !< x, y, z coordinates in local IB frame

        airfoil_id = patch_ib(patch_id)%airfoil_id
        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
        center(3) = patch_ib(patch_id)%z_centroid + real(zp, wp)*(z_domain%end - z_domain%beg)
        lz = patch_ib(patch_id)%length_z
        ca_in = ib_airfoil(airfoil_id)%c
        offset(:) = patch_ib(patch_id)%centroid_offset(:)

        z_max = lz/2
        z_min = -lz/2

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, zp, encoded_patch_id)

        ! find the indices to the left and right of the IB in i, j, k
        il = -gp_layers - 1
        jl = -gp_layers - 1
        ll = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        lr = p + gp_layers + 1
        ! maximum distance any marker can be from the center is the length of the airfoil
        call get_bounding_indices(center(1) - ca_in, center(1) + ca_in, x_cc, il, ir)
        call get_bounding_indices(center(2) - ca_in, center(2) + ca_in, y_cc, jl, jr)
        call get_bounding_indices(center(3) - ca_in, center(3) + ca_in, z_cc, ll, lr)

        $:GPU_PARALLEL_LOOP(private='[i, j, l, xyz_local, k, f]', copyin='[encoded_patch_id, center, offset, ca_in, airfoil_id, &
                            & ib_airfoil_grids(airfoil_id)%upper, ib_airfoil_grids(airfoil_id)%lower, z_min, z_max]', collapse=3)
        do l = ll, lr
            do j = jl, jr
                do i = il, ir
                    ! get coordinate frame centered on IB
                    xyz_local = [x_cc(i) - center(1), y_cc(j) - center(2), z_cc(l) - center(3)]
                    ! rotate the frame into the IB's coordinates
                    xyz_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse, xyz_local)
                    xyz_local = xyz_local - offset  ! airfoils are a patch that require a centroid offset

                    if (xyz_local(3) >= z_min .and. xyz_local(3) <= z_max) then
                        if (xyz_local(1) >= 0._wp .and. xyz_local(1) <= ca_in) then
                            if (xyz_local(2) >= 0._wp) then
                                k = 1
                                do while (ib_airfoil_grids(airfoil_id)%upper(k)%x < xyz_local(1))
                                    k = k + 1
                                end do
                                if (f_approx_equal(ib_airfoil_grids(airfoil_id)%upper(k)%x, xyz_local(1))) then
                                    if (xyz_local(2) <= ib_airfoil_grids(airfoil_id)%upper(k)%y) then
                                        ! IB
                                        ib_markers%sf(i, j, l) = encoded_patch_id
                                    end if
                                else
                                    f = (ib_airfoil_grids(airfoil_id)%upper(k)%x - xyz_local(1)) &
                                         & /(ib_airfoil_grids(airfoil_id)%upper(k)%x - ib_airfoil_grids(airfoil_id)%upper(k - 1)%x)
                                    if (xyz_local(2) <= ((1._wp - f)*ib_airfoil_grids(airfoil_id)%upper(k)%y &
                                        & + f*ib_airfoil_grids(airfoil_id)%upper(k - 1)%y)) then
                                        ib_markers%sf(i, j, l) = encoded_patch_id
                                    end if
                                end if
                            else
                                k = 1
                                do while (ib_airfoil_grids(airfoil_id)%lower(k)%x < xyz_local(1))
                                    k = k + 1
                                end do
                                if (f_approx_equal(ib_airfoil_grids(airfoil_id)%lower(k)%x, xyz_local(1))) then
                                    if (xyz_local(2) >= ib_airfoil_grids(airfoil_id)%lower(k)%y) then
                                        ib_markers%sf(i, j, l) = encoded_patch_id
                                    end if
                                else
                                    f = (ib_airfoil_grids(airfoil_id)%lower(k)%x - xyz_local(1)) &
                                         & /(ib_airfoil_grids(airfoil_id)%lower(k)%x - ib_airfoil_grids(airfoil_id)%lower(k - 1)%x)
                                    if (xyz_local(2) >= ((1._wp - f)*ib_airfoil_grids(airfoil_id)%lower(k)%y &
                                        & + f*ib_airfoil_grids(airfoil_id)%lower(k - 1)%y)) then
                                        ib_markers%sf(i, j, l) = encoded_patch_id
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

    !> Mark cells inside a rectangular immersed boundary
    subroutine s_ib_rectangle(patch_id, ib_markers, xp, yp)

        integer, intent(in)                :: patch_id
        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in)                :: xp, yp                !< integers containing the periodicity projection information
        integer                            :: i, j, il, ir, jl, jr  !< generic loop iterators
        integer                            :: encoded_patch_id
        real(wp)                           :: corner_distance       !< Equation of state parameters
        real(wp), dimension(1:3)           :: xy_local              !< x and y coordinates in local IB frame
        real(wp), dimension(1:2)           :: length, center        !< x and y coordinates in local IB frame
        real(wp), dimension(1:3,1:3)       :: inverse_rotation

        ! Transferring the rectangle's centroid and length information

        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, 0, encoded_patch_id)

        ! find the indices to the left and right of the IB in i, j, k
        il = -gp_layers - 1
        jl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        ! maximum distance any marker can be from the center
        corner_distance = 0.5_wp*sqrt(patch_ib(patch_id)%length_x**2 + patch_ib(patch_id)%length_y**2)
        call get_bounding_indices(center(1) - corner_distance, center(1) + corner_distance, x_cc, il, ir)
        call get_bounding_indices(center(2) - corner_distance, center(2) + corner_distance, y_cc, jl, jr)

        ! Assign primitive variables if rectangle covers cell and patch has write permission
        $:GPU_PARALLEL_LOOP(private='[i, j, xy_local]', copyin='[encoded_patch_id, center]', collapse=2)
        do j = jl, jr
            do i = il, ir
                ! get the x and y coordinates in the local IB frame
                xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
                xy_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse, xy_local)

                if (-0.5_wp*patch_ib(patch_id)%length_x <= xy_local(1) .and. 0.5_wp*patch_ib(patch_id)%length_x >= xy_local(1) &
                    & .and. -0.5_wp*patch_ib(patch_id)%length_y <= xy_local(2) &
                    & .and. 0.5_wp*patch_ib(patch_id)%length_y >= xy_local(2)) then
                    ! Updating the patch identities bookkeeping variable
                    ib_markers%sf(i, j, 0) = encoded_patch_id
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_rectangle

    !> Mark cells inside a spherical immersed boundary
    subroutine s_ib_sphere(patch_id, ib_markers, xp, yp, zp)

        integer, intent(in)                :: patch_id
        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in)                :: xp, yp, zp  !< integers containing the periodicity projection information
        ! Generic loop iterators
        integer                  :: i, j, k
        integer                  :: il, ir, jl, jr, kl, kr
        integer                  :: encoded_patch_id
        real(wp)                 :: radius
        real(wp), dimension(1:3) :: center

        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
        center(3) = patch_ib(patch_id)%z_centroid + real(zp, wp)*(z_domain%end - z_domain%beg)
        radius = patch_ib(patch_id)%radius

        ! completely skip particles not in the domain
        if (center(1) - radius > x_cc(m + gp_layers + 1) .or. center(1) + radius < x_cc(-gp_layers - 1) .or. center(2) &
            & - radius > y_cc(n + gp_layers + 1) .or. center(2) + radius < y_cc(-gp_layers - 1) .or. center(3) - radius > z_cc(p &
            & + gp_layers + 1) .or. center(3) + radius < z_cc(-gp_layers - 1)) then
            return
        end if

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, zp, encoded_patch_id)

        ! find the indices to the left and right of the IB in i, j, k
        il = -gp_layers - 1
        jl = -gp_layers - 1
        kl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        kr = p + gp_layers + 1
        call get_bounding_indices(center(1) - radius, center(1) + radius, x_cc, il, ir)
        call get_bounding_indices(center(2) - radius, center(2) + radius, y_cc, jl, jr)
        call get_bounding_indices(center(3) - radius, center(3) + radius, z_cc, kl, kr)

        ! Checking whether the sphere covers a particular cell in the domain and verifying whether the current patch has permission
        ! to write to that cell. If both queries check out, the primitive variables of the current patch are assigned to this cell.
        $:GPU_PARALLEL_LOOP(private='[i, j, k]', copyin='[encoded_patch_id, center, radius]', collapse=3)
        do k = kl, kr
            do j = jl, jr
                do i = il, ir
                    ! Updating the patch identities bookkeeping variable
                    if (((x_cc(i) - center(1))**2 + (y_cc(j) - center(2))**2 + (z_cc(k) - center(3))**2 <= radius**2)) then
                        ib_markers%sf(i, j, k) = encoded_patch_id
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_sphere

    !> Mark cells inside a cuboidal immersed boundary
    subroutine s_ib_cuboid(patch_id, ib_markers, xp, yp, zp)

        integer, intent(in)                :: patch_id
        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in)                :: xp, yp, zp  !< integers containing the periodicity projection information
        integer                            :: i, j, k, ir, il, jr, jl, kr, kl  !< Generic loop iterators
        integer                            :: encoded_patch_id
        real(wp), dimension(1:3)           :: xyz_local, center  !< x and y coordinates in local IB frame
        real(wp)                           :: corner_distance

        ! Transferring the cuboid's centroid

        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
        center(3) = patch_ib(patch_id)%z_centroid + real(zp, wp)*(z_domain%end - z_domain%beg)

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, zp, encoded_patch_id)

        ! find the indices to the left and right of the IB in i, j, k
        il = -gp_layers - 1
        jl = -gp_layers - 1
        kl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        kr = p + gp_layers + 1
        corner_distance = 0.5_wp*sqrt(patch_ib(patch_id)%length_x**2 + patch_ib(patch_id)%length_y**2 &
                                      & + patch_ib(patch_id)%length_z**2)  ! maximum distance any marker can be from the center
        call get_bounding_indices(center(1) - corner_distance, center(1) + corner_distance, x_cc, il, ir)
        call get_bounding_indices(center(2) - corner_distance, center(2) + corner_distance, y_cc, jl, jr)
        call get_bounding_indices(center(3) - corner_distance, center(3) + corner_distance, z_cc, kl, kr)

        ! Checking whether the cuboid covers a particular cell in the domain and verifying whether the current patch has permission
        ! to write to to that cell. If both queries check out, the primitive variables of the current patch are assigned to this
        ! cell.
        $:GPU_PARALLEL_LOOP(private='[i, j, k, xyz_local]', copyin='[encoded_patch_id, center]', collapse=3)
        do k = kl, kr
            do j = jl, jr
                do i = il, ir
                    xyz_local = [x_cc(i), y_cc(j), z_cc(k)] - center  ! get coordinate frame centered on IB
                    ! rotate the frame into the IB's coordinates
                    xyz_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse, xyz_local)

                    if (-0.5_wp*patch_ib(patch_id)%length_x <= xyz_local(1) &
                        & .and. 0.5_wp*patch_ib(patch_id)%length_x >= xyz_local(1) .and. &
                        & -0.5_wp*patch_ib(patch_id)%length_y <= xyz_local(2) &
                        & .and. 0.5_wp*patch_ib(patch_id)%length_y >= xyz_local(2) .and. &
                        & -0.5_wp*patch_ib(patch_id)%length_z <= xyz_local(3) &
                        & .and. 0.5_wp*patch_ib(patch_id)%length_z >= xyz_local(3)) then
                        ! Updating the patch identities bookkeeping variable
                        ib_markers%sf(i, j, k) = encoded_patch_id
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_cuboid

    !> Mark cells inside a cylindrical immersed boundary
    subroutine s_ib_cylinder(patch_id, ib_markers, xp, yp, zp)

        integer, intent(in)                :: patch_id
        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in)                :: xp, yp, zp  !< integers containing the periodicity projection information
        integer                            :: i, j, k, il, ir, jl, jr, kl, kr  !< Generic loop iterators
        integer                            :: encoded_patch_id
        real(wp), dimension(1:3)           :: xyz_local, center, length  !< x and y coordinates in local IB frame
        real(wp), dimension(1:3,1:3)       :: inverse_rotation
        real(wp)                           :: corner_distance

        ! Transferring the cylindrical patch's centroid

        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)
        center(3) = patch_ib(patch_id)%z_centroid + real(zp, wp)*(z_domain%end - z_domain%beg)
        length = [patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y, patch_ib(patch_id)%length_z]

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, zp, encoded_patch_id)

        il = -gp_layers - 1
        jl = -gp_layers - 1
        kl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        kr = p + gp_layers + 1
        corner_distance = sqrt(patch_ib(patch_id)%radius**2 + maxval(length)**2)  ! distance to rim of cylinder
        call get_bounding_indices(center(1) - corner_distance, center(1) + corner_distance, x_cc, il, ir)
        call get_bounding_indices(center(2) - corner_distance, center(2) + corner_distance, y_cc, jl, jr)
        call get_bounding_indices(center(3) - corner_distance, center(3) + corner_distance, z_cc, kl, kr)

        ! Checking whether the cylinder covers a particular cell in the domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out, the primitive variables of the current patch are assigned to
        ! this cell.
        $:GPU_PARALLEL_LOOP(private='[i, j, k, xyz_local]', copyin='[encoded_patch_id, center]', collapse=3)
        do k = kl, kr
            do j = jl, jr
                do i = il, ir
                    xyz_local = [x_cc(i), y_cc(j), z_cc(k)] - center  ! get coordinate frame centered on IB
                    ! rotate the frame into the IB's coordinates
                    xyz_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse, xyz_local)

                    if (((.not. f_is_default(patch_ib(patch_id)%length_x) .and. xyz_local(2)**2 + xyz_local(3) &
                        & **2 <= patch_ib(patch_id)%radius**2 .and. -0.5_wp*patch_ib(patch_id)%length_x <= xyz_local(1) &
                        & .and. 0.5_wp*patch_ib(patch_id)%length_x >= xyz_local(1)) &
                        & .or. (.not. f_is_default(patch_ib(patch_id)%length_y) .and. xyz_local(1)**2 + xyz_local(3) &
                        & **2 <= patch_ib(patch_id)%radius**2 .and. -0.5_wp*patch_ib(patch_id)%length_y <= xyz_local(2) &
                        & .and. 0.5_wp*patch_ib(patch_id)%length_y >= xyz_local(2)) &
                        & .or. (.not. f_is_default(patch_ib(patch_id)%length_z) .and. xyz_local(1)**2 + xyz_local(2) &
                        & **2 <= patch_ib(patch_id)%radius**2 .and. -0.5_wp*patch_ib(patch_id)%length_z <= xyz_local(3) &
                        & .and. 0.5_wp*patch_ib(patch_id)%length_z >= xyz_local(3)))) then
                        ! Updating the patch identities bookkeeping variable
                        ib_markers%sf(i, j, k) = encoded_patch_id
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_cylinder

    !> Mark cells inside a 2D elliptical immersed boundary
    subroutine s_ib_ellipse(patch_id, ib_markers, xp, yp)

        integer, intent(in)                :: patch_id
        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in)                :: xp, yp                !< integers containing the periodicity projection information
        integer                            :: i, j, il, ir, jl, jr  !< Generic loop iterators
        integer                            :: encoded_patch_id
        real(wp), dimension(1:3)           :: xy_local              !< x and y coordinates in local IB frame
        real(wp), dimension(1:2)           :: center                !< x and y coordinates in local IB frame
        real(wp)                           :: bounding_radius

        ! Transferring the ellipse's centroid

        center(1) = patch_ib(patch_id)%x_centroid + real(xp, wp)*(x_domain%end - x_domain%beg)
        center(2) = patch_ib(patch_id)%y_centroid + real(yp, wp)*(y_domain%end - y_domain%beg)

        ! encode the periodicity information into the patch_id
        call s_encode_patch_periodicity(patch_ib(patch_id)%gbl_patch_id, xp, yp, 0, encoded_patch_id)

        ! find the indices to the left and right of the IB in i, j, k
        bounding_radius = 0.5_wp*max(patch_ib(patch_id)%length_x, patch_ib(patch_id)%length_y)
        il = -gp_layers - 1
        jl = -gp_layers - 1
        ir = m + gp_layers + 1
        jr = n + gp_layers + 1
        call get_bounding_indices(center(1) - bounding_radius*2._wp, center(1) + bounding_radius*2._wp, x_cc, il, ir)
        call get_bounding_indices(center(2) - bounding_radius*2._wp, center(2) + bounding_radius*2._wp, y_cc, jl, jr)

        ! Checking whether the ellipse covers a particular cell in the domain
        $:GPU_PARALLEL_LOOP(private='[i, j, xy_local]', copyin='[encoded_patch_id, center]', collapse=2)
        do j = jl, jr
            do i = il, ir
                ! get the x and y coordinates in the local IB frame
                xy_local = [x_cc(i) - center(1), y_cc(j) - center(2), 0._wp]
                xy_local = matmul(patch_ib(patch_id)%rotation_matrix_inverse(:,:), xy_local)

                ! Ellipse condition (x/a)^2 + (y/b)^2 <= 1
                if ((xy_local(1)/(0.5_wp*patch_ib(patch_id)%length_x))**2 + (xy_local(2)/(0.5_wp*patch_ib(patch_id)%length_y)) &
                    & **2 <= 1._wp) then
                    ! Updating the patch identities bookkeeping variable
                    ib_markers%sf(i, j, 0) = encoded_patch_id
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_ellipse

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

        ! Local-space bounding box extents (min=1, max=2 in the third index)
        lx(1) = stl_bounding_boxes(model_id, 1, 1) + offset(1)
        lx(2) = stl_bounding_boxes(model_id, 1, 3) + offset(1)
        ly(1) = stl_bounding_boxes(model_id, 2, 1) + offset(2)
        ly(2) = stl_bounding_boxes(model_id, 2, 3) + offset(2)

        bbox_min = 1e12
        bbox_max = -1e12
        ! Enumerate all 8 corners of the local bounding box, rotate to world space, track world-space AABB
        do cx = 1, 2
            do cy = 1, 2
                local_corner = [lx(cx), ly(cy), 0._wp]
                world_corner = matmul(rotation, local_corner) + center
                bbox_min(1) = min(bbox_min(1), world_corner(1))
                bbox_min(2) = min(bbox_min(2), world_corner(2))
                bbox_max(1) = max(bbox_max(1), world_corner(1))
                bbox_max(2) = max(bbox_max(2), world_corner(2))
            end do
        end do

        call get_bounding_indices(bbox_min(1), bbox_max(1), x_cc, il, ir)
        call get_bounding_indices(bbox_min(2), bbox_max(2), y_cc, jl, jr)

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

        ! Local-space bounding box extents (min=1, max=2 in the third index)
        lx(1) = stl_bounding_boxes(model_id, 1, 1) + offset(1)
        lx(2) = stl_bounding_boxes(model_id, 1, 3) + offset(1)
        ly(1) = stl_bounding_boxes(model_id, 2, 1) + offset(2)
        ly(2) = stl_bounding_boxes(model_id, 2, 3) + offset(2)
        lz(1) = stl_bounding_boxes(model_id, 3, 1) + offset(3)
        lz(2) = stl_bounding_boxes(model_id, 3, 3) + offset(3)

        bbox_min = 1e12
        bbox_max = -1e12
        ! Enumerate all 8 corners of the local bounding box, rotate to world space, track world-space AABB
        do cx = 1, 2
            do cy = 1, 2
                do cz = 1, 2
                    local_corner = [lx(cx), ly(cy), lz(cz)]
                    world_corner = matmul(rotation, local_corner) + center
                    bbox_min(1) = min(bbox_min(1), world_corner(1))
                    bbox_min(2) = min(bbox_min(2), world_corner(2))
                    bbox_min(3) = min(bbox_min(3), world_corner(3))
                    bbox_max(1) = max(bbox_max(1), world_corner(1))
                    bbox_max(2) = max(bbox_max(2), world_corner(2))
                    bbox_max(3) = max(bbox_max(3), world_corner(3))
                end do
            end do
        end do

        call get_bounding_indices(bbox_min(1), bbox_max(1), x_cc, il, ir)
        call get_bounding_indices(bbox_min(2), bbox_max(2), y_cc, jl, jr)
        call get_bounding_indices(bbox_min(3), bbox_max(3), z_cc, kl, kr)

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

    subroutine get_bounding_indices(left_bound, right_bound, cell_centers, left_index, right_index)

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

    end subroutine get_bounding_indices

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
