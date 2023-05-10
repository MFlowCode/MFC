!>
!! @file m_patches.f90
!! @brief Contains module m_patches
module m_patches

    ! Dependencies =============================================================
    
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_helper

    use m_assign_variables
    ! ==========================================================================

    implicit none

    private; public :: s_line_segment, &
        s_spiral, &
        s_circle, &
        s_varcircle, &
        s_3dvarcircle, &
        s_ellipse, &
        s_ellipsoid, &
        s_rectangle, &
        s_sweep_line, &
        s_isentropic_vortex, &
        s_2D_TaylorGreen_vortex, &
        s_1D_analytical, &
        s_1d_bubble_pulse, &
        s_2D_analytical, &
        s_3D_analytical, &
        s_spherical_harmonic, &
        s_sphere, &
        s_cuboid, &
        s_cylinder, &
        s_sweep_plane


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

    real(wp) :: cart_y, cart_z
    real(wp) :: sph_phi !<
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates

    type(bounds_info) :: x_boundary, y_boundary, z_boundary  !<
    !! These variables combine the centroid and length parameters associated with
    !! a particular patch to yield the locations of the patch boundaries in the
    !! x-, y- and z-coordinate directions. They are used as a means to concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

contains
    !>          The line segment patch is a 1D geometry that may be used,
    !!              for example, in creating a Riemann problem. The geometry
    !!              of the patch is well-defined when its centroid and length
    !!              in the x-coordinate direction are provided. Note that the
    !!              line segment patch DOES NOT allow for the smearing of its
    !!              boundaries.
    !! @param patch_id patch identifier
    subroutine s_line_segment(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        real(wp) :: pi_inf, gamma, lit_gamma

        integer :: i, j  !< Generic loop operators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1._wp + gamma)/gamma

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

                !IF ( (q_prim_vf(1)%sf(i,0,0) < 1.e-12) .AND. (model_eqns .NE. 4)) THEN
                !    !zero density, reassign according to Tait EOS
                !    q_prim_vf(1)%sf(i,0,0) = &
                !        (((q_prim_vf(E_idx)%sf(i,0,0) + pi_inf)/(pref + pi_inf))**(1._wp/lit_gamma)) * &
                !        rhoref*(1._wp-q_prim_vf(alf_idx)%sf(i,0,0))
                !END IF
            end if
        end do

    end subroutine s_line_segment ! ----------------------------------------

    !>  The spiral patch is a 2D geometry that may be used, The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !!  @param patch_id patch identifier
    subroutine s_spiral(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        integer :: i, j, k !< Generic loop iterators
        real(wp) :: th, thickness, nturns, mya
        real(wp) :: spiral_x_min, spiral_x_max, spiral_y_min, spiral_y_max

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
                end if
            end do
        end do

    end subroutine s_spiral ! ----------------------------------------------

    !> The circular patch is a 2D geometry that may be used, for
        !!              example, in creating a bubble or a droplet. The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_circle(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: radius

        integer :: i, j !< Generic loop iterators

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
                end if

            end do
        end do

    end subroutine s_circle ! ----------------------------------------------

    !>             The varcircle patch is a 2D geometry that may be used
        !!             . It  generatres an annulus
        !! @param patch_id is the patch identifier
    subroutine s_varcircle(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: radius

        ! Generic loop iterators
        integer :: i, j

        real(wp) :: myr, thickness

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

                    q_prim_vf(alf_idx)%sf(i, j, 0) = patch_icpp(patch_id)%alpha(1)* &
                                                    exp(-0.5_wp*((myr - radius)**2._wp)/(thickness/3._wp)**2._wp)
                end if

            end do
        end do

    end subroutine s_varcircle ! ----------------------------------------------

    subroutine s_3dvarcircle(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: radius

        ! Generic loop iterators
        integer :: i, j, k

        real(wp) :: myr, thickness

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

                        q_prim_vf(alf_idx)%sf(i, j, k) = patch_icpp(patch_id)%alpha(1)* &
                                                        exp(-0.5_wp*((myr - radius)**2._wp)/(thickness/3._wp)**2._wp)
                    end if

                end do
            end do
        end do

    end subroutine s_3dvarcircle ! ----------------------------------------------

    !>      The elliptical patch is a 2D geometry. The geometry of
        !!      the patch is well-defined when its centroid and radii
        !!      are provided. Note that the elliptical patch DOES allow
        !!      for the smoothing of its boundary
        !! @param patch_id is the patch identifier
    subroutine s_ellipse(patch_id, patch_id_fp, q_prim_vf) ! ---------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: a, b

        integer :: i, j !< Generic loop operators

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
                end if
            end do
        end do

    end subroutine s_ellipse ! ---------------------------------------------

    !>      The ellipsoidal patch is a 3D geometry. The geometry of
        !!       the patch is well-defined when its centroid and radii
        !!       are provided. Note that the ellipsoidal patch DOES allow
        !!       for the smoothing of its boundary
        !! @param patch_id is the patch identifier
    subroutine s_ellipsoid(patch_id, patch_id_fp, q_prim_vf) ! -------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: a, b, c

        ! Generic loop iterators
        integer :: i, j, k

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
                    end if
                end do
            end do
        end do

    end subroutine s_ellipsoid ! -------------------------------------------

    !>      The rectangular patch is a 2D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, in alignment with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x- and y-
        !!              coordinate directions are provided. Please note that the
        !!              rectangular patch DOES NOT allow for the smoothing of its
        !!              boundaries.
        !! @param patch_id is the patch identifier
    subroutine s_rectangle(patch_id, patch_id_fp, q_prim_vf) ! -------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        real(wp) :: pi_inf, gamma, lit_gamma !< Equation of state parameters

        integer :: i, j !< generic loop iterators
        
        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1._wp + gamma)/gamma

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
                    y_boundary%end >= y_cc(j) &
                    .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                    eta, q_prim_vf, patch_id_fp)

                    if ((q_prim_vf(1)%sf(i, j, 0) < 1.e-10) .and. (model_eqns == 4)) then
                        !zero density, reassign according to Tait EOS
                        q_prim_vf(1)%sf(i, j, 0) = &
                            (((q_prim_vf(E_idx)%sf(i, j, 0) + pi_inf)/(pref + pi_inf))**(1._wp/lit_gamma))* &
                            rhoref*(1._wp - q_prim_vf(alf_idx)%sf(i, j, 0))
                    end if
                end if
            end do
        end do

    end subroutine s_rectangle ! -------------------------------------------

    !>  The swept line patch is a 2D geometry that may be used,
        !!      for example, in creating a solid boundary, or pre-/post-
        !!      shock region, at an angle with respect to the axes of the
        !!      Cartesian coordinate system. The geometry of the patch is
        !!      well-defined when its centroid and normal vector, aimed
        !!      in the sweep direction, are provided. Note that the sweep
        !!      line patch DOES allow the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_sweep_line(patch_id, patch_id_fp, q_prim_vf) ! ------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: a, b, c

        integer :: i, j !< Generic loop operators

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
                    eta = (5._wp * (10._wp ** -(1))) + (5._wp * (10._wp ** -(1)))*tanh(smooth_coeff/min(dx, dy) &
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

                end if

            end do
        end do

    end subroutine s_sweep_line ! ------------------------------------------

    !> The isentropic vortex is a 2D geometry that may be used,
        !!              for example, to generate an isentropic flow disturbance.
        !!              Geometry of the patch is well-defined when its centroid
        !!              and radius are provided. Notice that the patch DOES NOT
        !!              allow for the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_isentropic_vortex(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        ! Generic loop iterators
        integer :: i, j
        real(wp) :: radius

        ! Transferring isentropic vortex patch's centroid and radius info
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius

        ! Since the isentropic vortex patch does not allow for its boundary
        ! to get smoothed, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1._wp

        ! Verifying whether the isentropic vortex includes a particular cell
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries work out the primitive variables of the
        ! the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m

                if ((x_cc(i) - x_centroid)**2 &
                    + (y_cc(j) - y_centroid)**2 <= radius**2 &
                    .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, &
                                                            i, j, 0, &
                                            eta, q_prim_vf, patch_id_fp)

                end if

            end do
        end do

    end subroutine s_isentropic_vortex ! -----------------------------------
    
    !> The Taylor Green vortex is 2D decaying vortex that may be used,
        !!              for example, to verify the effects of viscous attenuation.
        !!              Geometry of the patch is well-defined when its centroid
        !!              are provided.
        !! @param patch_id is the patch identifier
    subroutine s_2D_TaylorGreen_Vortex(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< equation of state parameters
        real(kind(0d0)) :: L0, U0 !< Taylor Green Vortex parameters

        integer :: i, j !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1d0
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

                    ! Assign Parameters =========================================================
                    q_prim_vf(mom_idx%beg  )%sf(i,j,0) = U0*SIN(x_cc(i)/L0)*COS(y_cc(j)/L0)
                    q_prim_vf(mom_idx%end  )%sf(i,j,0) = -U0*COS(x_cc(i)/L0)*SIN(y_cc(j)/L0)
                    q_prim_vf(E_idx        )%sf(i,j,0) = patch_icpp(patch_id)%pres + (COS(2*x_cc(i))/L0 + &
                                                            COS(2*y_cc(j))/L0)* &
                                                            (q_prim_vf(1)%sf(i,j,0)*U0*U0)/16
                    ! ================================================================================

                end if
            end do
        end do

    end subroutine s_2D_TaylorGreen_Vortex ! -----------------------------------

    !>  This patch assigns the primitive variables as analytical
        !!  functions such that the code can be verified.
        !!  @param patch_id is the patch identifier
    subroutine s_1D_analytical(patch_id, patch_id_fp, q_prim_vf) ! ---------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        ! Placeholders for the cell boundary values
        real(wp) :: a, b, c, d, pi_inf, gamma, lit_gamma

        ! Generic loop iterators
        integer :: i, j

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1._wp + gamma)/gamma

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

                !what variables to alter
                !bump in pressure
                q_prim_vf(E_idx)%sf(i, 0, 0) = q_prim_vf(E_idx)%sf(i, 0, 0)* &
                                            (1._wp + 0.2_wp*exp(-1._wp*((x_cb(i) - x_centroid)**2._wp)/(2._wp*0.005_wp)))

                !bump in void fraction
                !q_prim_vf(adv_idx%beg)%sf(i,0,0) = q_prim_vf(adv_idx%beg)%sf(i,0,0) * &
                !    ( 1._wp + 0.2_wp*exp(-1._wp*((x_cb(i)-x_centroid)**2._wp)/(2._wp*0.005_wp)) )

                !bump in R(x)
                !q_prim_vf(adv_idx%end+1)%sf(i,0,0) = q_prim_vf(adv_idx%end+1)%sf(i,0,0) * &
                !    ( 1._wp + 0.2_wp*exp(-1._wp*((x_cb(i)-x_centroid)**2._wp)/(2._wp*0.005_wp)) )

                !IF (model_eqns == 4) THEN
                !reassign density
                !IF (num_fluids == 1) THEN
                ! q_prim_vf(1)%sf(i, 0, 0) = &
                !     (((q_prim_vf(E_idx)%sf(i, 0, 0) + pi_inf)/(pref + pi_inf))**(1._wp/lit_gamma))* &
                !     rhoref*(1._wp - q_prim_vf(alf_idx)%sf(i, 0, 0))
                !END IF
                !ELSE IF (model_eqns == 2) THEN
                !can manually adjust density here
                !q_prim_vf(1)%sf(i,0,0) = q_prim_vf(1)%sf(i,0,0) * &
                !    ( 1._wp + 0.2_wp*exp(-1._wp*((x_cb(i)-x_centroid)**2._wp)/(2._wp*0.005_wp)) )
                !END IF
            end if
        end do

    end subroutine s_1D_analytical ! ---------------------------------------

    subroutine s_1d_bubble_pulse(patch_id, patch_id_fp, q_prim_vf) ! ---------------------------------
        ! Description: This patch assigns the primitive variables as analytical
        !       functions such that the code can be verified.

        ! Patch identifier
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        ! Placeholders for the cell boundary values
        real(wp) :: fac, a, b, c, d, pi_inf, gamma, lit_gamma

        ! Generic loop iterators
        integer :: i, j

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1._wp + gamma)/gamma

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

                !what variables to alter
                !sinusoid in pressure
                q_prim_vf(E_idx)%sf(i, 0, 0) = q_prim_vf(E_idx)%sf(i, 0, 0)* &
                                            (1._wp + 0.1_wp*sin(-1._wp*(x_cb(i) - x_centroid)*2._wp*pi/length_x))

                !bump in void fraction
                !q_prim_vf(adv_idx%beg)%sf(i,0,0) = q_prim_vf(adv_idx%beg)%sf(i,0,0) * &
                !    ( 1._wp + 0.2_wp*exp(-1._wp*((x_cb(i)-x_centroid)**2._wp)/(2._wp*0.005_wp)) )

                !bump in R(x)
                !q_prim_vf(adv_idx%end+1)%sf(i,0,0) = q_prim_vf(adv_idx%end+1)%sf(i,0,0) * &
                !    ( 1._wp + 0.2_wp*exp(-1._wp*((x_cb(i)-x_centroid)**2._wp)/(2._wp*0.005_wp)) )

                !IF (model_eqns == 4) THEN
                !reassign density
                !IF (num_fluids == 1) THEN
                q_prim_vf(1)%sf(i, 0, 0) = &
                    (((q_prim_vf(E_idx)%sf(i, 0, 0) + pi_inf)/(pref + pi_inf))**(1._wp/lit_gamma))* &
                    rhoref*(1._wp - q_prim_vf(alf_idx)%sf(i, 0, 0))
                !END IF
                !ELSE IF (model_eqns == 2) THEN
                !can manually adjust density here
                !q_prim_vf(1)%sf(i,0,0) = q_prim_vf(1)%sf(i,0,0) * &
                !    ( 1._wp + 0.2_wp*exp(-1._wp*((x_cb(i)-x_centroid)**2._wp)/(2._wp*0.005_wp)) )
                !END IF
            end if
        end do

    end subroutine s_1D_bubble_pulse ! ---------------------------------------

    !>  This patch assigns the primitive variables as analytical
        !!  functions such that the code can be verified.
        !!  @param patch_id is the patch identifier
    subroutine s_2D_analytical(patch_id, patch_id_fp, q_prim_vf) ! ---------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        real(wp) :: a, b, c, d !< placeholderrs for the cell boundary values
        real(wp) :: pi_inf, gamma, lit_gamma !< equation of state parameters
        real(wp) :: a, b, c, d !< placeholderrs for the cell boundary values
        real(wp) :: pi_inf, gamma, lit_gamma !< equation of state parameters
        real(wp) :: l, U0 !< Taylor Green Vortex parameters

        integer :: i, j !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1._wp + gamma)/gamma

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
        l = 1._wp
        U0 = 10._wp ** -1
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

                    !what variables to alter
                    !x-y bump in pressure
                    q_prim_vf(E_idx)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)* &
                            (1._wp + 0.2_wp*exp(-1._wp*((x_cb(i) - x_centroid)**2._wp + (y_cb(j) - y_centroid)**2._wp)/(2._wp*0.005_wp)))

                    !x-bump
                    !q_prim_vf(E_idx)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)* &
                    !(1._wp + 0.2_wp*exp(-1._wp*((x_cb(i) - x_centroid)**2._wp)/(2._wp*0.005_wp)))

                    !bump in void fraction
                    !q_prim_vf(adv_idx%beg)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * &
                    !    ( 1._wp + 0.2_wp*exp(-1._wp*((x_cb(i)-x_centroid)**2._wp + (y_cb(j)-y_centroid)**2._wp)/(2._wp*0.005_wp)) )

                    !bump in R(x)
                    !q_prim_vf(adv_idx%end+1)%sf(i,j,0) = q_prim_vf(adv_idx%end+1)%sf(i,j,0) * &
                    !    ( 1._wp + 0.2_wp*exp(-1._wp*((x_cb(i)-x_centroid)**2._wp + (y_cb(j)-y_centroid)**2._wp)/(2._wp*0.005_wp)) )

                    !reassign density
                    !q_prim_vf(1)%sf(i, j, 0) = &
                    !(((q_prim_vf(E_idx)%sf(i, j, 0) + pi_inf)/(pref + pi_inf))**(1._wp/lit_gamma))* &
                    !rhoref*(1._wp - q_prim_vf(alf_idx)%sf(i, j, 0))

                    ! ================================================================================

                    ! Sinusoidal initial condition for all flow variables =============================

                    ! Cell-center values
    !                        a = 0._wp
    !                        b = 0._wp
    !                        c = 0._wp
    !                        d = 0._wp
    !                        q_prim_vf(adv_idx%beg)%sf(i,j,0) = SIN(x_cc(i)) * SIN(y_cc(j))
    !                        q_prim_vf(1)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * 1._wp
    !                        q_prim_vf(cont_idx%end)%sf(i,j,0) = (1._wp - q_prim_vf(adv_idx%beg)%sf(i,j,0)) * 1._wp
    !                        q_prim_vf(mom_idx%beg)%sf(i,j,0) = SIN(x_cc(i))
    !                        q_prim_vf(mom_idx%end)%sf(i,j,0) = SIN(y_cc(j))
    !                        q_prim_vf(E_idx)%sf(i,j,0) = 1._wp

                    ! Cell-average values
    !                       a = x_cc(i) - (5._wp * (10._wp ** -(1)))*dx ! x-beg
    !                       b = x_cc(i) + (5._wp * (10._wp ** -(1)))*dx ! x-end
    !                       c = y_cc(j) - (5._wp * (10._wp ** -(1)))*dy ! y-beg
    !                       d = y_cc(j) + (5._wp * (10._wp ** -(1)))*dy ! y-end
    !                       q_prim_vf(adv_idx%beg)%sf(i,j,0) = 1._wp/((b-a)*(d-c)) * &
    !                               (COS(a)*COS(c) - COS(a)*COS(d) - COS(b)*COS(c) + COS(b)*COS(d))
    !                       q_prim_vf(1)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * 1._wp
    !                       q_prim_vf(cont_idx%end)%sf(i,j,0) = (1._wp - q_prim_vf(adv_idx%beg)%sf(i,j,0)) * 1._wp
    !                       q_prim_vf(mom_idx%beg)%sf(i,j,0) = (COS(a) - COS(b))/(b-a)
    !                       q_prim_vf(mom_idx%end)%sf(i,j,0) = (COS(c) - COS(d))/(d-c)
    !                       q_prim_vf(E_idx)%sf(i,j,0) = 1._wp
                    ! ================================================================================
                   
                    ! Initial pressure profile smearing for bubble collapse case of Tiwari (2013) ====
                    !IF((       (x_cc(i))**2                     &
                    !         + (y_cc(j))**2 <= 1._wp**2)) THEN
                    !         q_prim_vf(E_idx)%sf(i,j,0) = (1._wp * (10._wp ** 5)) / 25._wp
                    !ELSE
                    !    q_prim_vf(E_idx)%sf(i,j,0) = (1._wp * (10._wp ** 5)) + 1._wp/SQRT(x_cc(i)**2+y_cc(j)**2) &
                    !                                    * (((1._wp * (10._wp ** 5))/25._wp) - (1._wp * (10._wp ** 5)))
                    !END IF
                    ! ================================================================================

                    ! Taylor Green Vortex===============================================
                    ! q_prim_vf(mom_idx%beg  )%sf(i,j,0) = U0*SIN(x_cc(i)/l)*COS(y_cc(j)/l)
                    ! q_prim_vf(mom_idx%end  )%sf(i,j,0) = -U0*COS(x_cc(i)/l)*SIN(y_cc(j)/l)
                    ! q_prim_vf(E_idx        )%sf(i,j,0) = 100000d0 + (COS(2*x_cc(i))/l + COS(2*y_cc(j))/l)* &
                    !                                         (q_prim_vf(1)%sf(i,j,0)*U0*U0)/16
                    ! ================================================================================

                end if
            end do
        end do

    end subroutine s_2D_analytical ! ---------------------------------------

    !>      This patch assigns the primitive variables as analytical
        !!      functions such that the code can be verified.
        !!      @param patch_id is the patch identifier
    subroutine s_3D_analytical(patch_id, patch_id_fp, q_prim_vf) ! ---------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: pi_inf, gamma, lit_gamma !< equation of state parameters

        integer :: i, j, k !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1._wp + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y
        length_z = patch_icpp(patch_id)%length_z

        ! Computing the beginning and the end x-, y- and z-coordinates of
        ! the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5_wp*length_x
        x_boundary%end = x_centroid + 0.5_wp*length_x
        y_boundary%beg = y_centroid - 0.5_wp*length_y
        y_boundary%end = y_centroid + 0.5_wp*length_y
        z_boundary%beg = z_centroid - 0.5_wp*length_z
        z_boundary%end = z_centroid + 0.5_wp*length_z

        ! Since the analytical patch does not allow for its boundaries to get
        ! smoothed out, the pseudo volume fraction is set to 1 to make sure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1._wp

        ! Checking whether the patch covers a particular cell in the domain
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
                        z_boundary%end >= cart_z &
                        .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                    eta, q_prim_vf, patch_id_fp)

                        !gaussian ball
                        !what variables to alter
                        !bump in pressure
                        q_prim_vf(E_idx)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k)* &
                                                    (1._wp + 0.2_wp*exp(-1._wp* &
                                                                        ((x_cb(i) - x_centroid)**2._wp + &
                                                                        (y_cb(j) - y_centroid)**2._wp + &
                                                                        (z_cb(k) - z_centroid)**2._wp) &
                                                                        /(2._wp*0.5_wp)))

                        !bump in void fraction
                        !                       q_prim_vf(adv_idx%beg)%sf(i, j, k) = q_prim_vf(adv_idx%beg)%sf(i, j, k)* &
                        !                                                           (1._wp + 0.2_wp*exp(-1._wp* &
                        !                                                                           ((x_cb(i) - x_centroid)**2._wp + (y_cb(j) - y_centroid)**2._wp + (z_cb(k) - z_centroid)**2._wp) &
                        !                                                                          /(2._wp*0.005_wp)))

                        !bump in R(x)
                        !       q_prim_vf(adv_idx%end + 1)%sf(i, j, k) = q_prim_vf(adv_idx%end + 1)%sf(i, j, k)* &
                        !                                               (1._wp + 0.2_wp*exp(-1._wp* &
                        !                                                                               ((x_cb(i) - x_centroid)**2._wp + (y_cb(j) - y_centroid)**2._wp + (z_cb(k) - z_centroid)**2._wp) &
    !                                                                                  /(2._wp*0.005_wp)))

                        !reassign density
                        !         q_prim_vf(1)%sf(i, j, k) = &
                        !             (((q_prim_vf(E_idx)%sf(i, j, k) + pi_inf)/(pref + pi_inf))**(1._wp/lit_gamma))* &
                        !            rhoref*(1._wp - q_prim_vf(E_idx + 1)%sf(i, j, k))

                        ! ================================================================================

                        ! Constant x-velocity in cylindrical grid ========================================
    !                        q_prim_vf(cont_idx%beg )%sf(i,j,k) = 1._wp
    !                        q_prim_vf(cont_idx%end )%sf(i,j,k) = 0._wp
    !                        q_prim_vf(mom_idx%beg  )%sf(i,j,k) = 0._wp
    !                        q_prim_vf(mom_idx%beg+1)%sf(i,j,k) = COS(z_cc(k))
    !                        q_prim_vf(mom_idx%end  )%sf(i,j,k) = -SIN(z_cc(k))
    !                        q_prim_vf(E_idx        )%sf(i,j,k) = 1._wp
    !                        q_prim_vf(adv_idx%beg  )%sf(i,j,k) = 1._wp
                        ! ================================================================================

                        ! Couette flow in cylindrical grid ===============================================
                        !q_prim_vf(cont_idx%beg )%sf(i,j,k) = 1._wp
                        !q_prim_vf(cont_idx%end )%sf(i,j,k) = 0._wp
                        !q_prim_vf(mom_idx%beg  )%sf(i,j,k) = 0._wp
                        !q_prim_vf(mom_idx%beg+1)%sf(i,j,k) = y_cc(j)*COS(z_cc(k))*SIN(z_cc(k))
                        !q_prim_vf(mom_idx%end  )%sf(i,j,k) = -y_cc(j)*SIN(z_cc(k))**2
                        !q_prim_vf(E_idx        )%sf(i,j,k) = 1._wp
                        !q_prim_vf(adv_idx%beg  )%sf(i,j,k) = 1._wp
                        ! ================================================================================

                    end if

                end do
            end do
        end do

    end subroutine s_3D_analytical ! ---------------------------------------

    !>      This patch generates the shape of the spherical harmonics
        !!      as a perturbation to a perfect sphere
        !!      @param patch_id is the patch identifier
    subroutine s_spherical_harmonic(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        real(wp) :: epsilon, beta
        real(wp) :: radius

        integer :: i, j, k !< generic loop iterators

        complex(wp) :: cmplx_i = (0._wp, 1._wp)
        complex(wp) :: H

        ! Transferring the patch's centroid and radius information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        radius = patch_icpp(patch_id)%radius
        epsilon = patch_icpp(patch_id)%epsilon
        beta = patch_icpp(patch_id)%beta

        ! Since the analytical patch does not allow for its boundaries to get
        ! smoothed out, the pseudo volume fraction is set to 1 to make sure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1._wp

        ! Checking whether the patch covers a particular cell in the domain
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

                    if (((x_cc(i) - x_centroid)**2 &
                        + (cart_y - y_centroid)**2 &
                        + (cart_z - z_centroid)**2 <= radius**2 &
                        .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k)))) &
                        then

                        call s_convert_cylindrical_to_spherical_coord(x_cc(i), y_cc(j))

                        if (epsilon == 1._wp) then
                            if (beta == 0._wp) then
                                H = (5._wp * (10._wp ** -(1)))*sqrt(3._wp/pi)*cos(sph_phi)
                            elseif (beta == 1._wp) then
                                H = (-5._wp * (10._wp ** -(1)))*sqrt(3._wp/(2._wp*pi))*exp(cmplx_i*z_cc(k))*sin(sph_phi)
                            end if
                        elseif (epsilon == 2._wp) then
                            if (beta == 0._wp) then
                                H = (25._wp * (10._wp ** -(2)))*sqrt(5._wp/pi)*(3._wp*cos(sph_phi)**2 - 1._wp)
                            elseif (beta == 1._wp) then
                                H = (-5._wp * (10._wp ** -(1)))*sqrt(15._wp/(2._wp*pi))*exp(cmplx_i*z_cc(k))*sin(sph_phi)*cos(sph_phi)
                            elseif (beta == 2._wp) then
                                H = (25._wp * (10._wp ** -(2)))*sqrt(15._wp/(2._wp*pi))*exp(2._wp*cmplx_i*z_cc(k))*sin(sph_phi)**2
                            end if
                        elseif (epsilon == 3._wp) then
                            if (beta == 0._wp) then
                                H = (25._wp * (10._wp ** -(2)))*sqrt(7._wp/pi)*(5._wp*cos(sph_phi)**3._wp - 3._wp*cos(sph_phi))
                            elseif (beta == 1._wp) then
                                H = (-125._wp * (10._wp ** -(3)))*sqrt(21._wp/pi)*exp(cmplx_i*z_cc(k))*sin(sph_phi)* &
                                    (5._wp*cos(sph_phi)**2 - 1._wp)
                            elseif (beta == 2._wp) then
                                H = (25._wp * (10._wp ** -(2)))*sqrt(105._wp/(2._wp*pi))*exp(2._wp*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**2*cos(sph_phi)
                            elseif (beta == 3._wp) then
                                H = (-125._wp * (10._wp ** -(3)))*sqrt(35._wp/pi)*exp(3._wp*cmplx_i*z_cc(k))*sin(sph_phi)**3._wp
                            end if
                        elseif (epsilon == 4._wp) then
                            if (beta == 0._wp) then
                                H = 3._wp/16._wp*sqrt(1._wp/pi)*(35._wp*cos(sph_phi)**4._wp - &
                                                        (3._wp * (10._wp ** 1))*cos(sph_phi)**2 + 3._wp)
                            elseif (beta == 1._wp) then
                                H = -3._wp/8._wp*sqrt(5._wp/pi)*exp(cmplx_i*z_cc(k))* &
                                    sin(sph_phi)*(7._wp*cos(sph_phi)**3._wp - 3._wp*cos(sph_phi))
                            elseif (beta == 2._wp) then
                                H = 3._wp/8._wp*sqrt(5._wp/(2._wp*pi))*exp(2._wp*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**2*(7._wp*cos(sph_phi)**2 - 1._wp)
                            elseif (beta == 3._wp) then
                                H = -3._wp/8._wp*sqrt(35._wp/pi)*exp(3._wp*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**3._wp*cos(sph_phi)
                            elseif (beta == 4._wp) then
                                H = 3._wp/16._wp*sqrt(35._wp/(2._wp*pi))*exp(4._wp*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**4._wp
                            end if
                        elseif (epsilon == 5._wp) then
                            if (beta == 0._wp) then
                                H = 1._wp/16._wp*sqrt(11._wp/pi)*(63._wp*cos(sph_phi)**5._wp - &
                                                            (7._wp * (10._wp ** 1))*cos(sph_phi)**3._wp + 15._wp*cos(sph_phi))
                            elseif (beta == 1._wp) then
                                H = -1._wp/16._wp*sqrt(165._wp/(2._wp*pi))*exp(cmplx_i*z_cc(k))* &
                                    sin(sph_phi)*(21._wp*cos(sph_phi)**4._wp - 14._wp*cos(sph_phi)**2 + 1._wp)
                            elseif (beta == 2._wp) then
                                H = (125._wp * (10._wp ** -(3)))*sqrt(1155._wp/(2._wp*pi))*exp(2._wp*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**2*(3._wp*cos(sph_phi)**3._wp - cos(sph_phi))
                            elseif (beta == 3._wp) then
                                H = -1._wp/32._wp*sqrt(385._wp/pi)*exp(3._wp*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**3._wp*(9._wp*cos(sph_phi)**2 - 1._wp)
                            elseif (beta == 4._wp) then
                                H = 3._wp/16._wp*sqrt(385._wp/(2._wp*pi))*exp(4._wp*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**4._wp*cos(sph_phi)
                            elseif (beta == 5._wp) then
                                H = -3._wp/32._wp*sqrt(77._wp/pi)*exp(5._wp*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**5._wp
                            end if
                        end if

                        q_prim_vf(adv_idx%beg)%sf(i, j, k) = 1._wp - abs(real(H, wp))

                    end if

                end do
            end do
        end do

    end subroutine s_spherical_harmonic ! ----------------------------------

    !>          The spherical patch is a 3D geometry that may be used,
        !!              for example, in creating a bubble or a droplet. The patch
        !!              geometry is well-defined when its centroid and radius are
        !!              provided. Please note that the spherical patch DOES allow
        !!              for the smoothing of its boundary.
        !!      @param patch_id is the patch identifier
    subroutine s_sphere(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: radius

        ! Generic loop iterators
        integer :: i, j, k !< generic loop iterators

        real(wp) :: radius_pressure, pressure_bubble, pressure_inf !<
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

                    if (((x_cc(i) - x_centroid)**2 &
                        + (cart_y - y_centroid)**2 &
                        + (cart_z - z_centroid)**2 <= radius**2 &
                        .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                    eta, q_prim_vf, patch_id_fp)

                    end if

                    ! Initialization of the pressure field that corresponds to the bubble-collapse
                    !! test case found in Tiwari et al. (2013)
                    ! radius_pressure = SQRT(x_cc(i)**2) ! 1D
                    ! radius_pressure = SQRT(x_cc(i)**2 + cart_y**2) ! 2D
                    ! radius_pressure = SQRT(x_cc(i)**2 + cart_y**2 + cart_z**2) ! 3D
                    ! pressure_bubble = 1.E+04
                    ! pressure_inf    = 1.E+05
                    ! q_prim_vf(E_idx)%sf(i,j,k) = pressure_inf + radius / radius_pressure * (pressure_bubble - pressure_inf)
                    !
                    ! IF(       ((  x_cc(i) - x_centroid)**2                    &
                    !          + (   cart_y - y_centroid)**2                    &
                    !          + (   cart_z - z_centroid)**2) <= radius**2)   &
                    !                               THEN
                    !
                    !    q_prim_vf(E_idx)%sf(i,j,k) = pressure_bubble
                    !
                    ! END IF

                end do
            end do
        end do

    end subroutine s_sphere ! ----------------------------------------------

    !>      The cuboidal patch is a 3D geometry that may be used, for
        !!              example, in creating a solid boundary, or pre-/post-shock
        !!              region, which is aligned with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x-, y- and
        !!              z-coordinate directions are provided. Please notice that
        !!              the cuboidal patch DOES NOT allow for the smearing of its
        !!              boundaries.
        !!      @param patch_id is the patch identifier
    subroutine s_cuboid(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        integer :: i, j, k !< Generic loop iterators

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
                        z_boundary%end >= cart_z &
                        .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                    eta, q_prim_vf, patch_id_fp)

                    end if
                end do
            end do
        end do

    end subroutine s_cuboid ! ----------------------------------------------

    !>              The cylindrical patch is a 3D geometry that may be used,
        !!              for example, in setting up a cylindrical solid boundary
        !!              confinement, like a blood vessel. The geometry of this
        !!              patch is well-defined when the centroid, the radius and
        !!              the length along the cylinder's axis, parallel to the x-,
        !!              y- or z-coordinate direction, are provided. Please note
        !!              that the cylindrical patch DOES allow for the smoothing
        !!              of its lateral boundary.
        !!      @param patch_id is the patch identifier
    subroutine s_cylinder(patch_id, patch_id_fp, q_prim_vf) ! --------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: radius

        integer :: i, j, k !< Generic loop iterators

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

                        if (length_x /= dflt_real) then
                            eta = tanh(smooth_coeff/min(dy, dz)* &
                                    (sqrt((cart_y - y_centroid)**2 &
                                            + (cart_z - z_centroid)**2) &
                                        - radius))*(-0.5_wp) + 0.5_wp
                        elseif (length_y /= dflt_real) then
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

                    if ((((length_x /= dflt_real .and. &
                        (cart_y - y_centroid)**2 &
                        + (cart_z - z_centroid)**2 <= radius**2 .and. &
                        x_boundary%beg <= x_cc(i) .and. &
                        x_boundary%end >= x_cc(i)) &
                        .or. &
                        (length_y /= dflt_real .and. &
                        (x_cc(i) - x_centroid)**2 &
                        + (cart_z - z_centroid)**2 <= radius**2 .and. &
                        y_boundary%beg <= cart_y .and. &
                        y_boundary%end >= cart_y) &
                        .or. &
                        (length_z /= dflt_real .and. &
                        (x_cc(i) - x_centroid)**2 &
                        + (cart_y - y_centroid)**2 <= radius**2 .and. &
                        z_boundary%beg <= cart_z .and. &
                        z_boundary%end >= cart_z)) &
                        .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                    eta, q_prim_vf, patch_id_fp)

                    end if

                end do
            end do
        end do

    end subroutine s_cylinder ! --------------------------------------------

    !>      The swept plane patch is a 3D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, at an angle with respect to the axes of the
        !!              Cartesian coordinate system. The geometry of the patch is
        !!              well-defined when its centroid and normal vector, aimed
        !!              in the sweep direction, are provided. Note that the sweep
        !!              plane patch DOES allow the smoothing of its boundary.
        !!      @param patch_id is the patch identifier
    subroutine s_sweep_plane(patch_id, patch_id_fp, q_prim_vf) ! -----------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: a, b, c, d

        integer :: i, j, k !< Generic loop iterators

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
                        eta = (5._wp * (10._wp ** -(1))) + (5._wp * (10._wp ** -(1)))*tanh(smooth_coeff/min(dx, dy, dz) &
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

                    end if

                end do
            end do
        end do

    end subroutine s_sweep_plane ! -----------------------------------------

    subroutine s_convert_cylindrical_to_cartesian_coord(cyl_y, cyl_z)
        !$acc routine seq
        real(wp), intent(IN) :: cyl_y, cyl_z

        cart_y = cyl_y*sin(cyl_z)
        cart_z = cyl_y*cos(cyl_z)

    end subroutine s_convert_cylindrical_to_cartesian_coord ! --------------

    subroutine s_convert_cylindrical_to_spherical_coord(cyl_x, cyl_y)
        !$acc routine seq
        real(wp), intent(IN) :: cyl_x, cyl_y

        sph_phi = atan(cyl_y/cyl_x)

    end subroutine s_convert_cylindrical_to_spherical_coord ! --------------

    !> Archimedes spiral function
    !! @param myth Angle
    !! @param offset Thickness
    !! @param a Starting position
    function f_r(myth, offset, a)
        !$acc routine seq
        real(wp), intent(IN) :: myth, offset, a
        real(wp) :: b
        real(wp) :: f_r

        !r(th) = a + b*th

        b = 2._wp*a/(2._wp*pi)
        f_r = a + b*myth + offset
    end function f_r

end module m_patches
