!>
!! @file m_patches.fpp
!! @brief Contains module m_patches

#:include 'case.fpp'
#:include 'ExtrusionHardcodedIC.fpp'
#:include '1dHardcodedIC.fpp'
#:include '2dHardcodedIC.fpp'
#:include '3dHardcodedIC.fpp'
#:include 'macros.fpp'

module m_icpp_patches

    use m_model                 ! Subroutine(s) related to STL files

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     !< Definitions of the global parameters

    use m_helper_basic          !< Functions to compare floating point numbers

    use m_helper

    use m_compute_levelset      ! Subroutines to calculate levelsets for IBs

    use m_mpi_common

    use m_assign_variables

    use m_mpi_common

    use m_ib_patches

    use m_variables_conversion

    implicit none

    private; public :: s_apply_icpp_patches

    real(wp) :: x_centroid, y_centroid, z_centroid
    real(wp) :: length_x, length_y, length_z

    integer :: smooth_patch_id
    real(wp) :: smooth_coeff !<
    !! These variables are analogous in both meaning and use to the similarly
    !! named components in the ic_patch_parameters type (see m_derived_types.f90
    !! for additional details). They are employed as a means to more concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    real(wp) :: eta !<
    !! In the case that smoothing of patch boundaries is enabled and the boundary
    !! between two adjacent patches is to be smeared out, this variable's purpose
    !! is to act as a pseudo volume fraction to indicate the contribution of each
    !! patch toward the composition of a cell's fluid state.

    real(wp) :: cart_x, cart_y, cart_z
    real(wp) :: sph_phi !<
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates

    type(bounds_info) :: x_boundary, y_boundary, z_boundary  !<
    !! These variables combine the centroid and length parameters associated with
    !! a particular patch to yield the locations of the patch boundaries in the
    !! x-, y- and z-coordinate directions. They are used as a means to concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    character(len=5) :: istr ! string to store int to string result for error checking

contains

    impure subroutine s_apply_icpp_patches(patch_id_fp, q_prim_vf)

        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        integer :: i

        !  3D Patch Geometries
        if (p > 0) then

            do i = 1, num_patches

                if (proc_rank == 0) then
                    print *, 'Processing patch', i
                end if

                !> ICPP Patches
                !> @{
                ! Spherical patch
                if (patch_icpp(i)%geometry == 8) then
                    call s_icpp_sphere(i, patch_id_fp, q_prim_vf)
                    ! Cuboidal patch
                elseif (patch_icpp(i)%geometry == 9) then
                    call s_icpp_cuboid(i, patch_id_fp, q_prim_vf)
                    ! Cylindrical patch
                elseif (patch_icpp(i)%geometry == 10) then
                    call s_icpp_cylinder(i, patch_id_fp, q_prim_vf)
                    ! Swept plane patch
                elseif (patch_icpp(i)%geometry == 11) then
                    call s_icpp_sweep_plane(i, patch_id_fp, q_prim_vf)
                    ! Ellipsoidal patch
                elseif (patch_icpp(i)%geometry == 12) then
                    call s_icpp_ellipsoid(i, patch_id_fp, q_prim_vf)
                    ! Spherical harmonic patch
                elseif (patch_icpp(i)%geometry == 14) then
                    call s_icpp_spherical_harmonic(i, patch_id_fp, q_prim_vf)
                    ! 3D Modified circular patch
                elseif (patch_icpp(i)%geometry == 19) then
                    call s_icpp_3dvarcircle(i, patch_id_fp, q_prim_vf)
                    ! 3D STL patch
                elseif (patch_icpp(i)%geometry == 21) then
                    call s_icpp_model(i, patch_id_fp, q_prim_vf)
                end if
            end do
            !> @}

            ! 2D Patch Geometries
        elseif (n > 0) then

            do i = 1, num_patches

                if (proc_rank == 0) then
                    print *, 'Processing patch', i
                end if

                !> ICPP Patches
                !> @{
                ! Circular patch
                if (patch_icpp(i)%geometry == 2) then
                    call s_icpp_circle(i, patch_id_fp, q_prim_vf)
                    ! Rectangular patch
                elseif (patch_icpp(i)%geometry == 3) then
                    call s_icpp_rectangle(i, patch_id_fp, q_prim_vf)
                    ! Swept line patch
                elseif (patch_icpp(i)%geometry == 4) then
                    call s_icpp_sweep_line(i, patch_id_fp, q_prim_vf)
                    ! Elliptical patch
                elseif (patch_icpp(i)%geometry == 5) then
                    call s_icpp_ellipse(i, patch_id_fp, q_prim_vf)
                    ! Unimplemented patch (formerly isentropic vortex)
                elseif (patch_icpp(i)%geometry == 6) then
                    call s_mpi_abort('This used to be the isentropic vortex patch, '// &
                                     'which no longer exists. See Examples. Exiting.')
                    ! Spherical Harmonic Patch
                elseif (patch_icpp(i)%geometry == 14) then
                    call s_icpp_spherical_harmonic(i, patch_id_fp, q_prim_vf)
                    ! Spiral patch
                elseif (patch_icpp(i)%geometry == 17) then
                    call s_icpp_spiral(i, patch_id_fp, q_prim_vf)
                    ! Modified circular patch
                elseif (patch_icpp(i)%geometry == 18) then
                    call s_icpp_varcircle(i, patch_id_fp, q_prim_vf)
                    ! TaylorGreen vortex patch
                elseif (patch_icpp(i)%geometry == 20) then
                    call s_icpp_2D_TaylorGreen_vortex(i, patch_id_fp, q_prim_vf)
                    ! STL patch
                elseif (patch_icpp(i)%geometry == 21) then
                    call s_icpp_model(i, patch_id_fp, q_prim_vf)
                end if
                !> @}
            end do

            ! 1D Patch Geometries
        else

            do i = 1, num_patches

                if (proc_rank == 0) then
                    print *, 'Processing patch', i
                end if

                ! Line segment patch
                if (patch_icpp(i)%geometry == 1) then
                    call s_icpp_line_segment(i, patch_id_fp, q_prim_vf)
                    ! 1d analytical
                elseif (patch_icpp(i)%geometry == 16) then
                    call s_icpp_1d_bubble_pulse(i, patch_id_fp, q_prim_vf)
                end if
            end do

        end if

    end subroutine s_apply_icpp_patches

    !>          The line segment patch is a 1D geometry that may be used,
    !!              for example, in creating a Riemann problem. The geometry
    !!              of the patch is well-defined when its centroid and length
    !!              in the x-coordinate direction are provided. Note that the
    !!              line segment patch DOES NOT allow for the smearing of its
    !!              boundaries.
    !! @param patch_id patch identifier
    !! @param patch_id_fp Array to track patch ids
    !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_line_segment(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        ! Generic loop iterators
        integer :: i, j, k

        ! Placeholders for the cell boundary values
        real(wp) :: pi_inf, gamma, lit_gamma
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded1DVariables()

        pi_inf = pi_infs(1)
        gamma = gammas(1)
        lit_gamma = gs_min(1)
        j = 0
        k = 0

        ! Transferring the line segment's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and end x-coordinates of the line segment
        ! based on its centroid and length
        x_boundary%beg = x_centroid - 0.5_wp*length_x
        x_boundary%end = x_centroid + 0.5_wp*length_x

        ! Since the line segment patch does not allow for its boundaries to
        ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1._wp

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0, 0))) then

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0, &
                                                        eta, q_prim_vf, patch_id_fp)

                @:analytical()

                ! check if this should load a hardcoded patch
                if (patch_icpp(patch_id)%hcid /= dflt_int) then
                    @:Hardcoded1D()
                end if

                ! Updating the patch identities bookkeeping variable
                if (1._wp - eta < sgm_eps) patch_id_fp(i, 0, 0) = patch_id

            end if
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_line_segment

    !>  The spiral patch is a 2D geometry that may be used, The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !! @param patch_id patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    impure subroutine s_icpp_spiral(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k !< Generic loop iterators
        real(wp) :: th, thickness, nturns, mya
        real(wp) :: spiral_x_min, spiral_x_max, spiral_y_min, spiral_y_max
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded2DVariables()

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        mya = patch_icpp(patch_id)%radius
        thickness = patch_icpp(patch_id)%length_x
        nturns = patch_icpp(patch_id)%length_y

        !
        logic_grid = 0
        do k = 0, int(m*91*nturns)
            th = k/real(int(m*91._wp*nturns))*nturns*2._wp*pi

            spiral_x_min = minval((/f_r(th, 0.0_wp, mya)*cos(th), &
                                    f_r(th, thickness, mya)*cos(th)/))
            spiral_y_min = minval((/f_r(th, 0.0_wp, mya)*sin(th), &
                                    f_r(th, thickness, mya)*sin(th)/))

            spiral_x_max = maxval((/f_r(th, 0.0_wp, mya)*cos(th), &
                                    f_r(th, thickness, mya)*cos(th)/))
            spiral_y_max = maxval((/f_r(th, 0.0_wp, mya)*sin(th), &
                                    f_r(th, thickness, mya)*sin(th)/))

            do j = 0, n; do i = 0, m; 
                    if ((x_cc(i) > spiral_x_min) .and. (x_cc(i) < spiral_x_max) .and. &
                        (y_cc(j) > spiral_y_min) .and. (y_cc(j) < spiral_y_max)) then
                        logic_grid(i, j, 0) = 1
                    end if
                end do; end do
        end do

        do j = 0, n
            do i = 0, m
                if ((logic_grid(i, j, 0) == 1)) then
                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()
                    if (patch_icpp(patch_id)%hcid /= dflt_int) then
                        @:Hardcoded2D()
                    end if

                    ! Updating the patch identities bookkeeping variable
                    if (1._wp - eta < sgm_eps) patch_id_fp(i, j, 0) = patch_id
                end if
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_spiral

    !> The circular patch is a 2D geometry that may be used, for
        !!              example, in creating a bubble or a droplet. The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_circle(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        real(wp) :: radius

        integer :: i, j, k !< Generic loop iterators
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded2DVariables()

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information

        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the circular patch's boundary is enabled.
        eta = 1._wp

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.

        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then

                    eta = tanh(smooth_coeff/min(dx, dy)* &
                               (sqrt((x_cc(i) - x_centroid)**2 &
                                     + (y_cc(j) - y_centroid)**2) &
                                - radius))*(-0.5_wp) + 0.5_wp

                end if

                if (((x_cc(i) - x_centroid)**2 &
                     + (y_cc(j) - y_centroid)**2 <= radius**2 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    .or. &
                    patch_id_fp(i, j, 0) == smooth_patch_id) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()
                    if (patch_icpp(patch_id)%hcid /= dflt_int) then
                        @:Hardcoded2D()
                    end if

                end if
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_circle

    !>  The varcircle patch is a 2D geometry that may be used
        !!             . It  generatres an annulus
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_varcircle(patch_id, patch_id_fp, q_prim_vf)

        ! Patch identifier
        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        ! Generic loop iterators
        integer :: i, j, k
        real(wp) :: radius, myr, thickness
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded2DVariables()

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff
        thickness = patch_icpp(patch_id)%epsilon

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the circular patch's boundary is enabled.
        eta = 1._wp

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m
                myr = sqrt((x_cc(i) - x_centroid)**2 &
                           + (y_cc(j) - y_centroid)**2)

                if (myr <= radius + thickness/2._wp .and. &
                    myr >= radius - thickness/2._wp .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()
                    if (patch_icpp(patch_id)%hcid /= dflt_int) then
                        @:Hardcoded2D()
                    end if

                    ! Updating the patch identities bookkeeping variable
                    if (1._wp - eta < sgm_eps) patch_id_fp(i, j, 0) = patch_id

                    q_prim_vf(alf_idx)%sf(i, j, 0) = patch_icpp(patch_id)%alpha(1)* &
                                                     exp(-0.5_wp*((myr - radius)**2._wp)/(thickness/3._wp)**2._wp)
                end if

            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_varcircle

    !! @param patch_id is the patch identifier
    !! @param patch_id_fp Array to track patch ids
    !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_3dvarcircle(patch_id, patch_id_fp, q_prim_vf)

        ! Patch identifier
        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        ! Generic loop iterators
        integer :: i, j, k
        real(wp) :: radius, myr, thickness
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded3DVariables()

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        length_z = patch_icpp(patch_id)%length_z
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff
        thickness = patch_icpp(patch_id)%epsilon

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the circular patch's boundary is enabled.
        eta = 1._wp

        ! write for all z

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    myr = sqrt((x_cc(i) - x_centroid)**2 &
                               + (y_cc(j) - y_centroid)**2)

                    if (myr <= radius + thickness/2._wp .and. &
                        myr >= radius - thickness/2._wp .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()
                        if (patch_icpp(patch_id)%hcid /= dflt_int) then
                            @:Hardcoded3D()
                        end if

                        ! Updating the patch identities bookkeeping variable
                        if (1._wp - eta < sgm_eps) patch_id_fp(i, j, k) = patch_id

                        q_prim_vf(alf_idx)%sf(i, j, k) = patch_icpp(patch_id)%alpha(1)* &
                                                         exp(-0.5_wp*((myr - radius)**2._wp)/(thickness/3._wp)**2._wp)
                    end if

                end do
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_3dvarcircle

    !> The elliptical patch is a 2D geometry. The geometry of
        !!      the patch is well-defined when its centroid and radii
        !!      are provided. Note that the elliptical patch DOES allow
        !!      for the smoothing of its boundary
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_ellipse(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k !< Generic loop operators
        real(wp) :: a, b
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded2DVariables()

        ! Transferring the elliptical patch's radii, centroid, smearing
        ! patch identity, and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        a = patch_icpp(patch_id)%radii(1)
        b = patch_icpp(patch_id)%radii(2)
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value
        ! be modified as the patch is laid out on the grid, but only in
        ! the case that smoothing of the elliptical patch's boundary is
        ! enabled.
        eta = 1._wp

        ! Checking whether the ellipse covers a particular cell in the
        ! domain and verifying whether the current patch has permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then
                    eta = tanh(smooth_coeff/min(dx, dy)* &
                               (sqrt(((x_cc(i) - x_centroid)/a)**2 + &
                                     ((y_cc(j) - y_centroid)/b)**2) &
                                - 1._wp))*(-0.5_wp) + 0.5_wp
                end if

                if ((((x_cc(i) - x_centroid)/a)**2 + &
                     ((y_cc(j) - y_centroid)/b)**2 <= 1._wp &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    .or. &
                    patch_id_fp(i, j, 0) == smooth_patch_id) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()
                    if (patch_icpp(patch_id)%hcid /= dflt_int) then
                        @:Hardcoded2D()
                    end if

                    ! Updating the patch identities bookkeeping variable
                    if (1._wp - eta < sgm_eps) patch_id_fp(i, j, 0) = patch_id
                end if
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_ellipse

    !> The ellipsoidal patch is a 3D geometry. The geometry of
        !!       the patch is well-defined when its centroid and radii
        !!       are provided. Note that the ellipsoidal patch DOES allow
        !!       for the smoothing of its boundary
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_ellipsoid(patch_id, patch_id_fp, q_prim_vf)

        ! Patch identifier
        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        ! Generic loop iterators
        integer :: i, j, k
        real(wp) :: a, b, c
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded3DVariables()

        ! Transferring the ellipsoidal patch's radii, centroid, smearing
        ! patch identity, and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        a = patch_icpp(patch_id)%radii(1)
        b = patch_icpp(patch_id)%radii(2)
        c = patch_icpp(patch_id)%radii(3)
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value
        ! be modified as the patch is laid out on the grid, but only in
        ! the case that smoothing of the ellipsoidal patch's boundary is
        ! enabled.
        eta = 1._wp

        ! Checking whether the ellipsoid covers a particular cell in the
        ! domain and verifying whether the current patch has permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (patch_icpp(patch_id)%smoothen) then
                        eta = tanh(smooth_coeff/min(dx, dy, dz)* &
                                   (sqrt(((x_cc(i) - x_centroid)/a)**2 + &
                                         ((cart_y - y_centroid)/b)**2 + &
                                         ((cart_z - z_centroid)/c)**2) &
                                    - 1._wp))*(-0.5_wp) + 0.5_wp
                    end if

                    if ((((x_cc(i) - x_centroid)/a)**2 + &
                         ((cart_y - y_centroid)/b)**2 + &
                         ((cart_z - z_centroid)/c)**2 <= 1._wp &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()
                        if (patch_icpp(patch_id)%hcid /= dflt_int) then
                            @:Hardcoded3D()
                        end if

                        ! Updating the patch identities bookkeeping variable
                        if (1._wp - eta < sgm_eps) patch_id_fp(i, j, k) = patch_id
                    end if
                end do
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_ellipsoid

    !> The rectangular patch is a 2D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, in alignment with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x- and y-
        !!              coordinate directions are provided. Please note that the
        !!              rectangular patch DOES NOT allow for the smoothing of its
        !!              boundaries.
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_rectangle(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k !< generic loop iterators
        real(wp) :: pi_inf, gamma, lit_gamma !< Equation of state parameters
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded2DVariables()

        pi_inf = pi_infs(1)
        gamma = gammas(1)
        lit_gamma = gs_min(1)

        ! Transferring the rectangle's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y

        ! Computing the beginning and the end x- and y-coordinates of the
        ! rectangle based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5_wp*length_x
        x_boundary%end = x_centroid + 0.5_wp*length_x
        y_boundary%beg = y_centroid - 0.5_wp*length_y
        y_boundary%end = y_centroid + 0.5_wp*length_y

        ! Since the rectangular patch does not allow for its boundaries to
        ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1._wp

        ! Checking whether the rectangle covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m
                if (x_boundary%beg <= x_cc(i) .and. &
                    x_boundary%end >= x_cc(i) .and. &
                    y_boundary%beg <= y_cc(j) .and. &
                    y_boundary%end >= y_cc(j)) then
                    if (patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()

                        if (patch_icpp(patch_id)%hcid /= dflt_int) then
                            @:Hardcoded2D()
                        end if

                        if ((q_prim_vf(1)%sf(i, j, 0) < 1.e-10) .and. (model_eqns == 4)) then
                            !zero density, reassign according to Tait EOS
                            q_prim_vf(1)%sf(i, j, 0) = &
                                (((q_prim_vf(E_idx)%sf(i, j, 0) + pi_inf)/(pref + pi_inf))**(1._wp/lit_gamma))* &
                                rhoref*(1._wp - q_prim_vf(alf_idx)%sf(i, j, 0))
                        end if

                        ! Updating the patch identities bookkeeping variable
                        if (1._wp - eta < sgm_eps) patch_id_fp(i, j, 0) = patch_id
                    end if
                end if
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_rectangle

    !> The swept line patch is a 2D geometry that may be used,
        !!      for example, in creating a solid boundary, or pre-/post-
        !!      shock region, at an angle with respect to the axes of the
        !!      Cartesian coordinate system. The geometry of the patch is
        !!      well-defined when its centroid and normal vector, aimed
        !!      in the sweep direction, are provided. Note that the sweep
        !!      line patch DOES allow the smoothing of its boundary.
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_sweep_line(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k !< Generic loop operators
        real(wp) :: a, b, c
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded3DVariables()

        ! Transferring the centroid information of the line to be swept
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Obtaining coefficients of the equation describing the sweep line
        a = patch_icpp(patch_id)%normal(1)
        b = patch_icpp(patch_id)%normal(2)
        c = -a*x_centroid - b*y_centroid

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the sweep line patch's boundary is enabled.
        eta = 1._wp

        ! Checking whether the region swept by the line covers a particular
        ! cell in the domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out, the
        ! primitive variables of the current patch are written to this cell.
        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then
                    eta = 5.e-1_wp + 5.e-1_wp*tanh(smooth_coeff/min(dx, dy) &
                                                   *(a*x_cc(i) + b*y_cc(j) + c) &
                                                   /sqrt(a**2 + b**2))
                end if

                if ((a*x_cc(i) + b*y_cc(j) + c >= 0._wp &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    .or. &
                    patch_id_fp(i, j, 0) == smooth_patch_id) &
                    then
                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()
                    if (patch_icpp(patch_id)%hcid /= dflt_int) then
                        @:Hardcoded3D()
                    end if

                    ! Updating the patch identities bookkeeping variable
                    if (1._wp - eta < sgm_eps) patch_id_fp(i, j, 0) = patch_id
                end if

            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_sweep_line

    !> The Taylor Green vortex is 2D decaying vortex that may be used,
        !!              for example, to verify the effects of viscous attenuation.
        !!              Geometry of the patch is well-defined when its centroid
        !!              are provided.
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_2D_TaylorGreen_Vortex(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k !< generic loop iterators
        real(wp) :: pi_inf, gamma, lit_gamma !< equation of state parameters
        real(wp) :: L0, U0 !< Taylor Green Vortex parameters
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded2DVariables()

        pi_inf = pi_infs(1)
        gamma = gammas(1)
        lit_gamma = gs_min(1)

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5_wp*length_x
        x_boundary%end = x_centroid + 0.5_wp*length_x
        y_boundary%beg = y_centroid - 0.5_wp*length_y
        y_boundary%end = y_centroid + 0.5_wp*length_y

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1._wp
        ! U0 is the characteristic velocity of the vortex
        U0 = patch_icpp(patch_id)%vel(1)
        ! L0 is the characteristic length of the vortex
        L0 = patch_icpp(patch_id)%vel(2)
        ! Checking whether the patch covers a particular cell in the
        ! domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out,
        ! the primitive variables of the current patch are assigned
        ! to this cell.
        do j = 0, n
            do i = 0, m
                if (x_boundary%beg <= x_cc(i) .and. &
                    x_boundary%end >= x_cc(i) .and. &
                    y_boundary%beg <= y_cc(j) .and. &
                    y_boundary%end >= y_cc(j) .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()
                    if (patch_icpp(patch_id)%hcid /= dflt_int) then
                        @:Hardcoded2D()
                    end if

                    ! Updating the patch identities bookkeeping variable
                    if (1._wp - eta < sgm_eps) patch_id_fp(i, j, 0) = patch_id

                    ! Assign Parameters
                    q_prim_vf(mom_idx%beg)%sf(i, j, 0) = U0*sin(x_cc(i)/L0)*cos(y_cc(j)/L0)
                    q_prim_vf(mom_idx%end)%sf(i, j, 0) = -U0*cos(x_cc(i)/L0)*sin(y_cc(j)/L0)
                    q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(patch_id)%pres + (cos(2*x_cc(i))/L0 + &
                                                                                cos(2*y_cc(j))/L0)* &
                                                   (q_prim_vf(1)%sf(i, j, 0)*U0*U0)/16
                end if
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_2D_TaylorGreen_Vortex

        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_1d_bubble_pulse(patch_id, patch_id_fp, q_prim_vf)
        ! Description: This patch assigns the primitive variables as analytical
        !       functions such that the code can be verified.

        ! Patch identifier
        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        ! Generic loop iterators
        integer :: i, j, k
        ! Placeholders for the cell boundary values
        real(wp) :: pi_inf, gamma, lit_gamma
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded1DVariables()

        pi_inf = pi_infs(1)
        gamma = gammas(1)
        lit_gamma = gs_min(1)

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5_wp*length_x
        x_boundary%end = x_centroid + 0.5_wp*length_x

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1._wp

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0, 0))) then

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0, &
                                                        eta, q_prim_vf, patch_id_fp)

                @:analytical()
                if (patch_icpp(patch_id)%hcid /= dflt_int) then
                    @:Hardcoded1D()
                end if

            end if
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_1D_bubble_pulse

    !> This patch generates the shape of the spherical harmonics
        !!      as a perturbation to a perfect sphere
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_spherical_harmonic(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(IN) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        real(wp) :: r, x_p, eps, phi
        real(wp), dimension(2:9) :: as, Ps
        real(wp) :: radius, x_centroid_local, y_centroid_local, z_centroid_local, eta_local, smooth_coeff_local
        logical :: non_axis_sym_in

        integer :: i, j, k !< generic loop iterators

        ! Transferring the patch's centroid and radius information
        x_centroid_local = patch_icpp(patch_id)%x_centroid
        y_centroid_local = patch_icpp(patch_id)%y_centroid
        z_centroid_local = patch_icpp(patch_id)%z_centroid
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff_local = patch_icpp(patch_id)%smooth_coeff
        radius = patch_icpp(patch_id)%radius
        as(2) = patch_icpp(patch_id)%a(2)
        as(3) = patch_icpp(patch_id)%a(3)
        as(4) = patch_icpp(patch_id)%a(4)
        as(5) = patch_icpp(patch_id)%a(5)
        as(6) = patch_icpp(patch_id)%a(6)
        as(7) = patch_icpp(patch_id)%a(7)
        as(8) = patch_icpp(patch_id)%a(8)
        as(9) = patch_icpp(patch_id)%a(9)
        non_axis_sym_in = patch_icpp(patch_id)%non_axis_sym

        ! Since the analytical patch does not allow for its boundaries to get
        ! smoothed out, the pseudo volume fraction is set to 1 to make sure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta_local = 1._wp
        eps = 1.e-32_wp

        ! Checking whether the patch covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! to that cell. If both queries check out, the primitive variables
        ! of the current patch are assigned to this cell.
        if (p > 0 .and. .not. non_axis_sym_in) then
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        if (grid_geometry == 3) then
                            call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                        else
                            cart_y = y_cc(j)
                            cart_z = z_cc(k)
                        end if

                        r = sqrt((x_cc(i) - x_centroid_local)**2 + (cart_y - y_centroid_local)**2 + (cart_z - z_centroid_local)**2) + eps
                        if (x_cc(i) - x_centroid_local <= 0) then
                            x_p = -1._wp*abs(x_cc(i) - x_centroid_local + eps)/r
                        else
                            x_p = abs(x_cc(i) - x_centroid_local + eps)/r
                        end if

                        Ps(2) = unassociated_legendre(x_p, 2)
                        Ps(3) = unassociated_legendre(x_p, 3)
                        Ps(4) = unassociated_legendre(x_p, 4)
                        Ps(5) = unassociated_legendre(x_p, 5)
                        Ps(6) = unassociated_legendre(x_p, 6)
                        Ps(7) = unassociated_legendre(x_p, 7)
                        if ((x_cc(i) - x_centroid_local >= 0 &
                             .and. &
                             r - as(2)*Ps(2) - as(3)*Ps(3) - as(4)*Ps(4) - as(5)*Ps(5) - as(6)*Ps(6) - as(7)*Ps(7) <= radius &
                             .and. &
                             patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) .or. &
                            (patch_id_fp(i, j, k) == smooth_patch_id)) &
                            then
                            if (patch_icpp(patch_id)%smoothen) then
                                eta_local = tanh(smooth_coeff_local/min(dx, dy, dz)* &
                                                 ((r - as(2)*Ps(2) - as(3)*Ps(3) - as(4)*Ps(4) - as(5)*Ps(5) - as(6)*Ps(6) - as(7)*Ps(7)) &
                                                  - radius))*(-0.5_wp) + 0.5_wp
                            end if

                            call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                    eta_local, q_prim_vf, patch_id_fp)
                        end if

                    end do
                end do
            end do

        else if (p == 0) then
            do j = 0, n
                do i = 0, m

                    if (non_axis_sym_in) then
                        phi = atan(((y_cc(j) - y_centroid_local) + eps)/((x_cc(i) - x_centroid_local) + eps))
                        r = sqrt((x_cc(i) - x_centroid_local)**2._wp + (y_cc(j) - y_centroid_local)**2._wp) + eps
                        x_p = (eps)/r
                        Ps(2) = spherical_harmonic_func(x_p, phi, 2, 2)
                        Ps(3) = spherical_harmonic_func(x_p, phi, 3, 3)
                        Ps(4) = spherical_harmonic_func(x_p, phi, 4, 4)
                        Ps(5) = spherical_harmonic_func(x_p, phi, 5, 5)
                        Ps(6) = spherical_harmonic_func(x_p, phi, 6, 6)
                        Ps(7) = spherical_harmonic_func(x_p, phi, 7, 7)
                        Ps(8) = spherical_harmonic_func(x_p, phi, 8, 8)
                        Ps(9) = spherical_harmonic_func(x_p, phi, 9, 9)
                    else
                        r = sqrt((x_cc(i) - x_centroid_local)**2._wp + (y_cc(j) - y_centroid_local)**2._wp) + eps
                        x_p = abs(x_cc(i) - x_centroid_local + eps)/r
                        Ps(2) = unassociated_legendre(x_p, 2)
                        Ps(3) = unassociated_legendre(x_p, 3)
                        Ps(4) = unassociated_legendre(x_p, 4)
                        Ps(5) = unassociated_legendre(x_p, 5)
                        Ps(6) = unassociated_legendre(x_p, 6)
                        Ps(7) = unassociated_legendre(x_p, 7)
                        Ps(8) = unassociated_legendre(x_p, 8)
                        Ps(9) = unassociated_legendre(x_p, 9)
                    end if

                    if (x_cc(i) - x_centroid_local >= 0 &
                        .and. &
                        r - as(2)*Ps(2) - as(3)*Ps(3) - as(4)*Ps(4) - as(5)*Ps(5) - as(6)*Ps(6) - as(7)*Ps(7) - as(8)*Ps(8) - as(9)*Ps(9) <= radius .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                        then
                        call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                                eta_local, q_prim_vf, patch_id_fp)

                    elseif (x_cc(i) - x_centroid_local < 0 &
                            .and. &
                            r - as(2)*Ps(2) + as(3)*Ps(3) - as(4)*Ps(4) + as(5)*Ps(5) - as(6)*Ps(6) + as(7)*Ps(7) - as(8)*Ps(8) + as(9)*Ps(9) <= radius &
                            .and. &
                            patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                        then
                        call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                                eta_local, q_prim_vf, patch_id_fp)

                    end if
                end do
            end do
        end if

    end subroutine s_icpp_spherical_harmonic

    !>          The spherical patch is a 3D geometry that may be used,
        !!              for example, in creating a bubble or a droplet. The patch
        !!              geometry is well-defined when its centroid and radius are
        !!              provided. Please note that the spherical patch DOES allow
        !!              for the smoothing of its boundary.
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_sphere(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        ! Generic loop iterators
        integer :: i, j, k
        real(wp) :: radius
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded3DVariables()

        !! Variables to initialize the pressure field that corresponds to the
            !! bubble-collapse test case found in Tiwari et al. (2013)

        ! Transferring spherical patch's radius, centroid, smoothing patch
        ! identity and smoothing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the spherical patch's boundary is enabled.
        eta = 1._wp

        ! Checking whether the sphere covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (patch_icpp(patch_id)%smoothen) then
                        eta = tanh(smooth_coeff/min(dx, dy, dz)* &
                                   (sqrt((x_cc(i) - x_centroid)**2 &
                                         + (cart_y - y_centroid)**2 &
                                         + (cart_z - z_centroid)**2) &
                                    - radius))*(-0.5_wp) + 0.5_wp
                    end if

                    if ((((x_cc(i) - x_centroid)**2 &
                          + (cart_y - y_centroid)**2 &
                          + (cart_z - z_centroid)**2 <= radius**2) .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()
                        if (patch_icpp(patch_id)%hcid /= dflt_int) then
                            @:Hardcoded3D()
                        end if

                    end if
                end do
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_sphere

    !> The cuboidal patch is a 3D geometry that may be used, for
        !!              example, in creating a solid boundary, or pre-/post-shock
        !!              region, which is aligned with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x-, y- and
        !!              z-coordinate directions are provided. Please notice that
        !!              the cuboidal patch DOES NOT allow for the smearing of its
        !!              boundaries.
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_cuboid(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k !< Generic loop iterators
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded3DVariables()

        ! Transferring the cuboid's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y
        length_z = patch_icpp(patch_id)%length_z

        ! Computing the beginning and the end x-, y- and z-coordinates of
        ! the cuboid based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5_wp*length_x
        x_boundary%end = x_centroid + 0.5_wp*length_x
        y_boundary%beg = y_centroid - 0.5_wp*length_y
        y_boundary%end = y_centroid + 0.5_wp*length_y
        z_boundary%beg = z_centroid - 0.5_wp*length_z
        z_boundary%end = z_centroid + 0.5_wp*length_z

        ! Since the cuboidal patch does not allow for its boundaries to get
        ! smoothed out, the pseudo volume fraction is set to 1 to make sure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1._wp

        ! Checking whether the cuboid covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! to that cell. If both queries check out, the primitive variables
        ! of the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (x_boundary%beg <= x_cc(i) .and. &
                        x_boundary%end >= x_cc(i) .and. &
                        y_boundary%beg <= cart_y .and. &
                        y_boundary%end >= cart_y .and. &
                        z_boundary%beg <= cart_z .and. &
                        z_boundary%end >= cart_z) then

                        if (patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) then

                            call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                    eta, q_prim_vf, patch_id_fp)

                            @:analytical()
                            if (patch_icpp(patch_id)%hcid /= dflt_int) then
                                @:Hardcoded3D()
                            end if

                            ! Updating the patch identities bookkeeping variable
                            if (1._wp - eta < sgm_eps) patch_id_fp(i, j, k) = patch_id

                        end if
                    end if
                end do
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_cuboid

    !> The cylindrical patch is a 3D geometry that may be used,
        !!              for example, in setting up a cylindrical solid boundary
        !!              confinement, like a blood vessel. The geometry of this
        !!              patch is well-defined when the centroid, the radius and
        !!              the length along the cylinder's axis, parallel to the x-,
        !!              y- or z-coordinate direction, are provided. Please note
        !!              that the cylindrical patch DOES allow for the smoothing
        !!              of its lateral boundary.
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Array of primitive variables
    subroutine s_icpp_cylinder(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k !< Generic loop iterators
        real(wp) :: radius
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded3DVariables()

        ! Transferring the cylindrical patch's centroid, length, radius,
        ! smoothing patch identity and smoothing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y
        length_z = patch_icpp(patch_id)%length_z
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Computing the beginning and the end x-, y- and z-coordinates of
        ! the cylinder based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5_wp*length_x
        x_boundary%end = x_centroid + 0.5_wp*length_x
        y_boundary%beg = y_centroid - 0.5_wp*length_y
        y_boundary%end = y_centroid + 0.5_wp*length_y
        z_boundary%beg = z_centroid - 0.5_wp*length_z
        z_boundary%end = z_centroid + 0.5_wp*length_z

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smearing of the cylindrical patch's boundary is enabled.
        eta = 1._wp

        ! Checking whether the cylinder covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (patch_icpp(patch_id)%smoothen) then
                        if (.not. f_is_default(length_x)) then
                            eta = tanh(smooth_coeff/min(dy, dz)* &
                                       (sqrt((cart_y - y_centroid)**2 &
                                             + (cart_z - z_centroid)**2) &
                                        - radius))*(-0.5_wp) + 0.5_wp
                        elseif (.not. f_is_default(length_y)) then
                            eta = tanh(smooth_coeff/min(dx, dz)* &
                                       (sqrt((x_cc(i) - x_centroid)**2 &
                                             + (cart_z - z_centroid)**2) &
                                        - radius))*(-0.5_wp) + 0.5_wp
                        else
                            eta = tanh(smooth_coeff/min(dx, dy)* &
                                       (sqrt((x_cc(i) - x_centroid)**2 &
                                             + (cart_y - y_centroid)**2) &
                                        - radius))*(-0.5_wp) + 0.5_wp
                        end if
                    end if

                    if (((.not. f_is_default(length_x) .and. &
                          (cart_y - y_centroid)**2 &
                          + (cart_z - z_centroid)**2 <= radius**2 .and. &
                          x_boundary%beg <= x_cc(i) .and. &
                          x_boundary%end >= x_cc(i)) &
                         .or. &
                         (.not. f_is_default(length_y) .and. &
                          (x_cc(i) - x_centroid)**2 &
                          + (cart_z - z_centroid)**2 <= radius**2 .and. &
                          y_boundary%beg <= cart_y .and. &
                          y_boundary%end >= cart_y) &
                         .or. &
                         (.not. f_is_default(length_z) .and. &
                          (x_cc(i) - x_centroid)**2 &
                          + (cart_y - y_centroid)**2 <= radius**2 .and. &
                          z_boundary%beg <= cart_z .and. &
                          z_boundary%end >= cart_z) .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()
                        if (patch_icpp(patch_id)%hcid /= dflt_int) then
                            @:Hardcoded3D()
                        end if

                        ! Updating the patch identities bookkeeping variable
                        if (1._wp - eta < sgm_eps) patch_id_fp(i, j, k) = patch_id
                    end if
                end do
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_cylinder

    !>      The swept plane patch is a 3D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, at an angle with respect to the axes of the
        !!              Cartesian coordinate system. The geometry of the patch is
        !!              well-defined when its centroid and normal vector, aimed
        !!              in the sweep direction, are provided. Note that the sweep
        !!              plane patch DOES allow the smoothing of its boundary.
        !! @param patch_id is the patch identifier
        !! @param patch_id_fp Array to track patch ids
        !! @param q_prim_vf Primitive variables
    subroutine s_icpp_sweep_plane(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k !< Generic loop iterators
        real(wp) :: a, b, c, d
        @:HardcodedDimensionsExtrusion()
        @:Hardcoded3DVariables()

        ! Transferring the centroid information of the plane to be swept
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Obtaining coefficients of the equation describing the sweep plane
        a = patch_icpp(patch_id)%normal(1)
        b = patch_icpp(patch_id)%normal(2)
        c = patch_icpp(patch_id)%normal(3)
        d = -a*x_centroid - b*y_centroid - c*z_centroid

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smearing of the sweep plane patch's boundary is enabled.
        eta = 1._wp

        ! Checking whether the region swept by the plane covers a particular
        ! cell in the domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out, the
        ! primitive variables of the current patch are written to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (patch_icpp(patch_id)%smoothen) then
                        eta = 5.e-1_wp + 5.e-1_wp*tanh(smooth_coeff/min(dx, dy, dz) &
                                                       *(a*x_cc(i) + &
                                                         b*cart_y + &
                                                         c*cart_z + d) &
                                                       /sqrt(a**2 + b**2 + c**2))
                    end if

                    if ((a*x_cc(i) + b*cart_y + c*cart_z + d >= 0._wp &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()
                        if (patch_icpp(patch_id)%hcid /= dflt_int) then
                            @:Hardcoded3D()
                        end if

                        ! Updating the patch identities bookkeeping variable
                        if (1._wp - eta < sgm_eps) patch_id_fp(i, j, k) = patch_id
                    end if

                end do
            end do
        end do
        @:HardcodedDellacation()

    end subroutine s_icpp_sweep_plane

    !> The STL patch is a 2/3D geometry that is imported from an STL file.
    !! @param patch_id is the patch identifier
    !! @param patch_id_fp Array to track patch ids
    !! @param q_prim_vf Primitive variables
    !! @param STL_levelset STL levelset
    !! @param STL_levelset_norm STL levelset normals
    subroutine s_icpp_model(patch_id, patch_id_fp, q_prim_vf)

        integer, intent(in) :: patch_id
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        ! Variables for IBM+STL
        real(wp) :: normals(1:3) !< Boundary normal buffer
        integer :: boundary_vertex_count, boundary_edge_count, total_vertices !< Boundary vertex
        real(wp), allocatable, dimension(:, :, :) :: boundary_v !< Boundary vertex buffer
        real(wp), allocatable, dimension(:, :) :: interpolated_boundary_v !< Interpolated vertex buffer
        real(wp) :: distance !< Levelset distance buffer
        logical :: interpolate !< Logical variable to determine whether or not the model should be interpolated

        integer :: i, j, k !< Generic loop iterators

        type(t_bbox) :: bbox, bbox_old
        type(t_model) :: model
        type(ic_model_parameters) :: params

        real(wp), dimension(1:3) :: point, model_center

        real(wp) :: grid_mm(1:3, 1:2)

        integer :: cell_num
        integer :: ncells

        real(wp), dimension(1:4, 1:4) :: transform, transform_n

        if (proc_rank == 0) then
            print *, " * Reading model: "//trim(patch_icpp(patch_id)%model_filepath)
        end if

        model = f_model_read(patch_icpp(patch_id)%model_filepath)
        params%scale(:) = patch_icpp(patch_id)%model_scale(:)
        params%translate(:) = patch_icpp(patch_id)%model_translate(:)
        params%rotate(:) = patch_icpp(patch_id)%model_rotate(:)
        params%spc = patch_icpp(patch_id)%model_spc
        params%threshold = patch_icpp(patch_id)%model_threshold

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
            call f_check_interpolation_2D(boundary_v, boundary_edge_count, (/dx, dy, dz/), interpolate)
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

            !call s_model_write("__out__.stl", model)
            !call s_model_write("__out__.obj", model)

            grid_mm(1, :) = (/minval(x_cc) - 0.e5_wp*dx, maxval(x_cc) + 0.e5_wp*dx/)
            grid_mm(2, :) = (/minval(y_cc) - 0.e5_wp*dy, maxval(y_cc) + 0.e5_wp*dy/)

            if (p > 0) then
                grid_mm(3, :) = (/minval(z_cc) - 0.e5_wp*dz, maxval(z_cc) + 0.e5_wp*dz/)
            else
                grid_mm(3, :) = (/0._wp, 0._wp/)
            end if

            write (*, "(A, 3(2X, F20.10))") "    > Domain: Min:", grid_mm(:, 1)
            write (*, "(A, 3(2X, F20.10))") "    >         Cen:", (grid_mm(:, 1) + grid_mm(:, 2))/2._wp
            write (*, "(A, 3(2X, F20.10))") "    >         Max:", grid_mm(:, 2)
        end if

        ncells = (m + 1)*(n + 1)*(p + 1)
        do i = 0, m; do j = 0, n; do k = 0, p

                    cell_num = i*(n + 1)*(p + 1) + j*(p + 1) + (k + 1)
                    if (proc_rank == 0 .and. mod(cell_num, ncells/100) == 0) then
                        write (*, "(A, I3, A)", advance="no") &
                            char(13)//"  * Generating grid: ", &
                            nint(100*real(cell_num)/ncells), "%"
                    end if

                    point = (/x_cc(i), y_cc(j), 0._wp/)
                    if (p > 0) then
                        point(3) = z_cc(k)
                    end if

                    if (grid_geometry == 3) then
                        point = f_convert_cyl_to_cart(point)
                    end if

                    eta = f_model_is_inside(model, point, (/dx, dy, dz/), patch_icpp(patch_id)%model_spc)

                    if (eta > patch_icpp(patch_id)%model_threshold) then
                        eta = 1._wp
                    else if (.not. patch_icpp(patch_id)%smoothen) then
                        eta = 0._wp
                    end if

                    call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                            eta, q_prim_vf, patch_id_fp)

                    ! Note: Should probably use *eta* to compute primitive variables
                    ! if defining them analytically.
                    @:analytical()
                end do; end do; end do

        if (proc_rank == 0) then
            print *, ""
            print *, " * Cleaning up."
        end if

        call s_model_free(model)

    end subroutine s_icpp_model

    subroutine s_convert_cylindrical_to_cartesian_coord(cyl_y, cyl_z)
        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: cyl_y, cyl_z

        cart_y = cyl_y*sin(cyl_z)
        cart_z = cyl_y*cos(cyl_z)

    end subroutine s_convert_cylindrical_to_cartesian_coord

    function f_convert_cyl_to_cart(cyl) result(cart)

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
    elemental function f_r(myth, offset, a)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: myth, offset, a
        real(wp) :: b
        real(wp) :: f_r

        !r(th) = a + b*th

        b = 2._wp*a/(2._wp*pi)
        f_r = a + b*myth + offset
    end function f_r

end module m_icpp_patches
