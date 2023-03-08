module m_create_patches

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    use m_assign_patches
    ! ==========================================================================

    implicit none


    real(kind(0d0)) :: radius
    real(kind(0d0)) :: smooth_coeff !<
    !! These variables are analogous in both meaning and use to the similarly
    !! named components in the ic_patch_parameters type (see m_derived_types.f90
    !! for additional details). They are employed as a means to more concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    type(bounds_info) :: x_boundary, y_boundary, z_boundary  !<
    !! These variables combine the centroid and length parameters associated with
    !! a particular patch to yield the locations of the patch boundaries in the
    !! x-, y- and z-coordinate directions. They are used as a means to concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    real(kind(0d0)) :: length_x, length_y, length_z

    real(kind(0d0)) :: a, b, c, d !<
    !! When a line or a plane sweep patch geometry is employed, these variables
    !! represent the coefficients associated with the equation describing the
    !! said line or plane.

    real(kind(0d0)) :: cart_y, cart_z
    real(kind(0d0)) :: sph_phi !<
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates
    
contains

    subroutine s_convert_cylindrical_to_cartesian_coord(cyl_y, cyl_z)

        real(kind(0d0)), intent(IN) :: cyl_y, cyl_z

        cart_y = cyl_y*sin(cyl_z)
        cart_z = cyl_y*cos(cyl_z)

    end subroutine s_convert_cylindrical_to_cartesian_coord ! --------------

    subroutine s_convert_cylindrical_to_spherical_coord(cyl_x, cyl_y)

        real(kind(0d0)), intent(IN) :: cyl_x, cyl_y

        sph_phi = atan(cyl_y/cyl_x)

    end subroutine s_convert_cylindrical_to_spherical_coord ! --------------

    subroutine s_perturb_sphere() ! ----------------------------------------

        integer :: i, j, k, l !< generic loop operators

        real(kind(0d0)) :: perturb_alpha
        real(kind(0d0)) :: alpha_unadv
        real(kind(0d0)) :: rand_real
        call random_seed()

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    call random_number(rand_real)

                    perturb_alpha = q_prim_vf(E_idx + perturb_sph_fluid)%sf(i, j, k)

                    ! Perturb partial density fields to match perturbed volume fraction fields
!                        IF ((perturb_alpha >= 25d-2) .AND. (perturb_alpha <= 75d-2)) THEN
                    if ((perturb_alpha /= 0d0) .and. (perturb_alpha /= 1d0)) then

                        ! Derive new partial densities
                        do l = 1, num_fluids
                            q_prim_vf(l)%sf(i, j, k) = q_prim_vf(E_idx + l)%sf(i, j, k)*fluid_rho(l)
                        end do

                    end if
                end do
            end do
        end do

    end subroutine s_perturb_sphere ! --------------------------------------

    subroutine s_perturb_surrounding_flow() ! ------------------------------

        integer :: i, j, k, l !<  generic loop iterators

        real(kind(0d0)) :: perturb_alpha
        real(kind(0d0)) :: rand_real
        call random_seed()

        ! Perturb partial density or velocity of surrounding flow by some random small amount of noise
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    perturb_alpha = q_prim_vf(E_idx + perturb_flow_fluid)%sf(i, j, k)
                    ! IF (perturb_alpha == 1d0) THEN
                    ! Perturb partial density
!                            CALL RANDOM_NUMBER(rand_real)
!                            rand_real = rand_real / 1d2 / 1d3
!                            q_prim_vf(perturb_flow_fluid)%sf(i,j,k) = q_prim_vf(perturb_flow_fluid)%sf(i,j,k) + rand_real
                    ! Perturb velocity
                    call random_number(rand_real)
                    rand_real = rand_real*1.d-2
                    q_prim_vf(mom_idx%beg)%sf(i, j, k) = (1.d0 + rand_real)*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    q_prim_vf(mom_idx%end)%sf(i, j, k) = rand_real*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    if (bubbles) then
                        q_prim_vf(alf_idx)%sf(i, j, k) = (1.d0 + rand_real)*q_prim_vf(alf_idx)%sf(i, j, k)
                    end if
                    ! END IF
                end do
            end do
        end do

    end subroutine s_perturb_surrounding_flow ! ----------------------------

    !>          The line segment patch is a 1D geometry that may be used,
        !!              for example, in creating a Riemann problem. The geometry
        !!              of the patch is well-defined when its centroid and length
        !!              in the x-coordinate direction are provided. Note that the
        !!              line segment patch DOES NOT allow for the smearing of its
        !!              boundaries.
        !! @param patch_id patch identifier
    subroutine s_line_segment(patch_id) ! ----------------------------------

        integer, intent(IN) :: patch_id

        real(kind(0d0)) :: pi_inf, gamma, lit_gamma

        integer :: i, j  !< Generic loop operators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the line segment's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and end x-coordinates of the line segment
        ! based on its centroid and length
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x

        ! Since the line segment patch does not allow for its boundaries to
        ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0, 0))) then

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0)

                !IF ( (q_prim_vf(1)%sf(i,0,0) < 1.e-12) .AND. (model_eqns .NE. 4)) THEN
                !    !zero density, reassign according to Tait EOS
                !    q_prim_vf(1)%sf(i,0,0) = &
                !        (((q_prim_vf(E_idx)%sf(i,0,0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma)) * &
                !        rhoref*(1d0-q_prim_vf(alf_idx)%sf(i,0,0))
                !END IF
            end if
        end do

    end subroutine s_line_segment ! ----------------------------------------

    !>  The spiral patch is a 2D geometry that may be used, The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !!  @param patch_id patch identifier
    subroutine s_spiral(patch_id) ! ----------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j, k !< Generic loop iterators
        real(kind(0d0)) :: th, thickness, nturns, mya
        real(kind(0d0)) :: spiral_x_min, spiral_x_max, spiral_y_min, spiral_y_max

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
            th = k/real(int(m*91d0*nturns))*nturns*2.d0*pi

            spiral_x_min = minval((/f_r(th, 0.0d0, mya)*cos(th), &
                                    f_r(th, thickness, mya)*cos(th)/))
            spiral_y_min = minval((/f_r(th, 0.0d0, mya)*sin(th), &
                                    f_r(th, thickness, mya)*sin(th)/))

            spiral_x_max = maxval((/f_r(th, 0.0d0, mya)*cos(th), &
                                    f_r(th, thickness, mya)*cos(th)/))
            spiral_y_max = maxval((/f_r(th, 0.0d0, mya)*sin(th), &
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
                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)
                end if
            end do
        end do

    end subroutine s_spiral ! ----------------------------------------------

    !> Archimedes spiral function
        !! @param myth Angle
        !! @param offset Thickness
        !! @param a Starting position
    function f_r(myth, offset, a)
        real(kind(0d0)), intent(IN) :: myth, offset, a
        real(kind(0d0)) :: b
        real(kind(0d0)) :: f_r

        !r(th) = a + b*th

        b = 2.d0*a/(2.d0*pi)
        f_r = a + b*myth + offset
    end function f_r

    !> The circular patch is a 2D geometry that may be used, for
        !!              example, in creating a bubble or a droplet. The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_circle(patch_id) ! ----------------------------------------

        integer, intent(IN) :: patch_id

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
        eta = 1d0

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
                                - radius))*(-0.5d0) + 0.5d0

                end if

                if (((x_cc(i) - x_centroid)**2 &
                     + (y_cc(j) - y_centroid)**2 <= radius**2 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    .or. &
                    patch_id_fp(i, j, 0) == smooth_patch_id) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

                end if

            end do
        end do

    end subroutine s_circle ! ----------------------------------------------

    !>             The varcircle patch is a 2D geometry that may be used
        !!             . It  generatres an annulus
        !! @param patch_id is the patch identifier
    subroutine s_varcircle(patch_id) ! ----------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j

        real(kind(0d0)) :: myr, thickness

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
        eta = 1d0

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m
                myr = dsqrt((x_cc(i) - x_centroid)**2 &
                            + (y_cc(j) - y_centroid)**2)

                if (myr <= radius + thickness/2.d0 .and. &
                    myr >= radius - thickness/2.d0 .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

                    q_prim_vf(alf_idx)%sf(i, j, 0) = patch_icpp(patch_id)%alpha(1)* &
                                                     dexp(-0.5d0*((myr - radius)**2.d0)/(thickness/3.d0)**2.d0)
                end if

            end do
        end do

    end subroutine s_varcircle ! ----------------------------------------------

    subroutine s_3dvarcircle(patch_id) ! ----------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j, k

        real(kind(0d0)) :: myr, thickness

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
        eta = 1d0

        ! write for all z

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    myr = dsqrt((x_cc(i) - x_centroid)**2 &
                                + (y_cc(j) - y_centroid)**2)

                    if (myr <= radius + thickness/2.d0 .and. &
                        myr >= radius - thickness/2.d0 .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

                        q_prim_vf(alf_idx)%sf(i, j, k) = patch_icpp(patch_id)%alpha(1)* &
                                                         dexp(-0.5d0*((myr - radius)**2.d0)/(thickness/3.d0)**2.d0)
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
    subroutine s_ellipse(patch_id) ! ---------------------------------------

        integer, intent(IN) :: patch_id

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
        eta = 1d0

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
                                - 1d0))*(-0.5d0) + 0.5d0
                end if

                if ((((x_cc(i) - x_centroid)/a)**2 + &
                     ((y_cc(j) - y_centroid)/b)**2 <= 1d0 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    .or. &
                    patch_id_fp(i, j, 0) == smooth_patch_id) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)
                end if
            end do
        end do

    end subroutine s_ellipse ! ---------------------------------------------

    !>      The ellipsoidal patch is a 3D geometry. The geometry of
        !!       the patch is well-defined when its centroid and radii
        !!       are provided. Note that the ellipsoidal patch DOES allow
        !!       for the smoothing of its boundary
        !! @param patch_id is the patch identifier
    subroutine s_ellipsoid(patch_id) ! -------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

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
        eta = 1d0

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
                                    - 1d0))*(-0.5d0) + 0.5d0
                    end if

                    if ((((x_cc(i) - x_centroid)/a)**2 + &
                         ((cart_y - y_centroid)/b)**2 + &
                         ((cart_z - z_centroid)/c)**2 <= 1d0 &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)
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
    subroutine s_rectangle(patch_id) ! -------------------------------------

        integer, intent(IN) :: patch_id

        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< Equation of state parameters

        integer :: i, j !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the rectangle's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y

        ! Computing the beginning and the end x- and y-coordinates of the
        ! rectangle based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y

        ! Since the rectangular patch does not allow for its boundaries to
        ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

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

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

                    if ((q_prim_vf(1)%sf(i, j, 0) < 1.e-10) .and. (model_eqns == 4)) then
                        !zero density, reassign according to Tait EOS
                        q_prim_vf(1)%sf(i, j, 0) = &
                            (((q_prim_vf(E_idx)%sf(i, j, 0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                            rhoref*(1d0 - q_prim_vf(alf_idx)%sf(i, j, 0))
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
    subroutine s_sweep_line(patch_id) ! ------------------------------------

        integer, intent(IN) :: patch_id

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
        eta = 1d0

        ! Checking whether the region swept by the line covers a particular
        ! cell in the domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out, the
        ! primitive variables of the current patch are written to this cell.
        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then
                    eta = 5d-1 + 5d-1*tanh(smooth_coeff/min(dx, dy) &
                                           *(a*x_cc(i) + b*y_cc(j) + c) &
                                           /sqrt(a**2 + b**2))
                end if

                if ((a*x_cc(i) + b*y_cc(j) + c >= 0d0 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    .or. &
                    patch_id_fp(i, j, 0) == smooth_patch_id) &
                    then
                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

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
    subroutine s_isentropic_vortex(patch_id) ! ----------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j

        ! Transferring isentropic vortex patch's centroid and radius info
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius

        ! Since the isentropic vortex patch does not allow for its boundary
        ! to get smoothed, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

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
                                                            i, j, 0)

                end if

            end do
        end do

    end subroutine s_isentropic_vortex ! -----------------------------------

    !>  This patch assigns the primitive variables as analytical
        !!  functions such that the code can be verified.
        !!  @param patch_id is the patch identifier
    subroutine s_1D_analytical(patch_id) ! ---------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Placeholders for the cell boundary values
        real(kind(0d0)) :: a, b, c, d, pi_inf, gamma, lit_gamma

        ! Generic loop iterators
        integer :: i, j

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1d0

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0, 0))) then

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0)

                !what variables to alter
                !bump in pressure
                q_prim_vf(E_idx)%sf(i, 0, 0) = q_prim_vf(E_idx)%sf(i, 0, 0)* &
                                               (1d0 + 0.2d0*dexp(-1d0*((x_cb(i) - x_centroid)**2.d0)/(2.d0*0.005d0)))

                !bump in void fraction
                !q_prim_vf(adv_idx%beg)%sf(i,0,0) = q_prim_vf(adv_idx%beg)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                !bump in R(x)
                !q_prim_vf(adv_idx%end+1)%sf(i,0,0) = q_prim_vf(adv_idx%end+1)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                !IF (model_eqns == 4) THEN
                !reassign density
                !IF (num_fluids == 1) THEN
                ! q_prim_vf(1)%sf(i, 0, 0) = &
                !     (((q_prim_vf(E_idx)%sf(i, 0, 0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                !     rhoref*(1d0 - q_prim_vf(alf_idx)%sf(i, 0, 0))
                !END IF
                !ELSE IF (model_eqns == 2) THEN
                !can manually adjust density here
                !q_prim_vf(1)%sf(i,0,0) = q_prim_vf(1)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )
                !END IF
            end if
        end do

    end subroutine s_1D_analytical ! ---------------------------------------

    subroutine s_1d_bubble_pulse(patch_id) ! ---------------------------------
        ! Description: This patch assigns the primitive variables as analytical
        !       functions such that the code can be verified.

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Placeholders for the cell boundary values
        real(kind(0d0)) :: fac, a, b, c, d, pi_inf, gamma, lit_gamma

        ! Generic loop iterators
        integer :: i, j

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1d0

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0, 0))) then

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0)

                !what variables to alter
                !sinusoid in pressure
                q_prim_vf(E_idx)%sf(i, 0, 0) = q_prim_vf(E_idx)%sf(i, 0, 0)* &
                                               (1d0 + 0.1d0*sin(-1d0*(x_cb(i) - x_centroid)*2d0*pi/length_x))

                !bump in void fraction
                !q_prim_vf(adv_idx%beg)%sf(i,0,0) = q_prim_vf(adv_idx%beg)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                !bump in R(x)
                !q_prim_vf(adv_idx%end+1)%sf(i,0,0) = q_prim_vf(adv_idx%end+1)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                !IF (model_eqns == 4) THEN
                !reassign density
                !IF (num_fluids == 1) THEN
                q_prim_vf(1)%sf(i, 0, 0) = &
                    (((q_prim_vf(E_idx)%sf(i, 0, 0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                    rhoref*(1d0 - q_prim_vf(alf_idx)%sf(i, 0, 0))
                !END IF
                !ELSE IF (model_eqns == 2) THEN
                !can manually adjust density here
                !q_prim_vf(1)%sf(i,0,0) = q_prim_vf(1)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )
                !END IF
            end if
        end do

    end subroutine s_1D_bubble_pulse ! ---------------------------------------

    !>  This patch assigns the primitive variables as analytical
        !!  functions such that the code can be verified.
        !!  @param patch_id is the patch identifier
    subroutine s_2D_analytical(patch_id) ! ---------------------------------

        integer, intent(IN) :: patch_id

        real(kind(0d0)) :: a, b, c, d !< placeholderrs for the cell boundary values
        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< equation of state parameters

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

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

                    !what variables to alter
                    !x-y bump in pressure
                    q_prim_vf(E_idx)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)* &
                               (1d0 + 0.2d0*dexp(-1d0*((x_cb(i) - x_centroid)**2.d0 + (y_cb(j) - y_centroid)**2.d0)/(2.d0*0.005d0)))

                    !x-bump
                    !q_prim_vf(E_idx)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)* &
                    !(1d0 + 0.2d0*dexp(-1d0*((x_cb(i) - x_centroid)**2.d0)/(2.d0*0.005d0)))

                    !bump in void fraction
                    !q_prim_vf(adv_idx%beg)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * &
                    !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0)/(2.d0*0.005d0)) )

                    !bump in R(x)
                    !q_prim_vf(adv_idx%end+1)%sf(i,j,0) = q_prim_vf(adv_idx%end+1)%sf(i,j,0) * &
                    !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0)/(2.d0*0.005d0)) )

                    !reassign density
                    !q_prim_vf(1)%sf(i, j, 0) = &
                    !(((q_prim_vf(E_idx)%sf(i, j, 0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                    !rhoref*(1d0 - q_prim_vf(alf_idx)%sf(i, j, 0))

                    ! ================================================================================

                    ! Sinusoidal initial condition for all flow variables =============================

                    ! Cell-center values
!                        a = 0d0
!                        b = 0d0
!                        c = 0d0
!                        d = 0d0
!                        q_prim_vf(adv_idx%beg)%sf(i,j,0) = SIN(x_cc(i)) * SIN(y_cc(j))
!                        q_prim_vf(1)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * 1d0
!                        q_prim_vf(cont_idx%end)%sf(i,j,0) = (1d0 - q_prim_vf(adv_idx%beg)%sf(i,j,0)) * 1d0
!                        q_prim_vf(mom_idx%beg)%sf(i,j,0) = SIN(x_cc(i))
!                        q_prim_vf(mom_idx%end)%sf(i,j,0) = SIN(y_cc(j))
!                        q_prim_vf(E_idx)%sf(i,j,0) = 1d0

                    ! Cell-average values
!                       a = x_cc(i) - 5d-1*dx ! x-beg
!                       b = x_cc(i) + 5d-1*dx ! x-end
!                       c = y_cc(j) - 5d-1*dy ! y-beg
!                       d = y_cc(j) + 5d-1*dy ! y-end
!                       q_prim_vf(adv_idx%beg)%sf(i,j,0) = 1d0/((b-a)*(d-c)) * &
!                               (COS(a)*COS(c) - COS(a)*COS(d) - COS(b)*COS(c) + COS(b)*COS(d))
!                       q_prim_vf(1)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * 1d0
!                       q_prim_vf(cont_idx%end)%sf(i,j,0) = (1d0 - q_prim_vf(adv_idx%beg)%sf(i,j,0)) * 1d0
!                       q_prim_vf(mom_idx%beg)%sf(i,j,0) = (COS(a) - COS(b))/(b-a)
!                       q_prim_vf(mom_idx%end)%sf(i,j,0) = (COS(c) - COS(d))/(d-c)
!                       q_prim_vf(E_idx)%sf(i,j,0) = 1d0
                    ! ================================================================================

                    ! Initial pressure profile smearing for bubble collapse case of Tiwari (2013) ====
                    !IF((       (x_cc(i))**2                     &
                    !         + (y_cc(j))**2 <= 1d0**2)) THEN
                    !         q_prim_vf(E_idx)%sf(i,j,0) = 1d5 / 25d0
                    !ELSE
                    !    q_prim_vf(E_idx)%sf(i,j,0) = 1d5 + 1d0/SQRT(x_cc(i)**2+y_cc(j)**2) &
                    !                                    * ((1d5/25d0) - 1d5)
                    !END IF
                    ! ================================================================================

                end if
            end do
        end do

    end subroutine s_2D_analytical ! ---------------------------------------

    !>      This patch assigns the primitive variables as analytical
        !!      functions such that the code can be verified.
        !!      @param patch_id is the patch identifier
    subroutine s_3D_analytical(patch_id) ! ---------------------------------

        integer, intent(IN) :: patch_id
        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< equation of state parameters

        integer :: i, j, k !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y
        length_z = patch_icpp(patch_id)%length_z

        ! Computing the beginning and the end x-, y- and z-coordinates of
        ! the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y
        z_boundary%beg = z_centroid - 0.5d0*length_z
        z_boundary%end = z_centroid + 0.5d0*length_z

        ! Since the analytical patch does not allow for its boundaries to get
        ! smoothed out, the pseudo volume fraction is set to 1 to make sure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

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

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

                        !gaussian ball
                        !what variables to alter
                        !bump in pressure
                        q_prim_vf(E_idx)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k)* &
                                                       (1d0 + 0.2d0*exp(-1d0* &
                                                                        ((x_cb(i) - x_centroid)**2.d0 + &
                                                                         (y_cb(j) - y_centroid)**2.d0 + &
                                                                         (z_cb(k) - z_centroid)**2.d0) &
                                                                        /(2.d0*0.5d0)))

                        !bump in void fraction
                        !                       q_prim_vf(adv_idx%beg)%sf(i, j, k) = q_prim_vf(adv_idx%beg)%sf(i, j, k)* &
                        !                                                           (1d0 + 0.2d0*exp(-1d0* &
                        !                                                                           ((x_cb(i) - x_centroid)**2.d0 + (y_cb(j) - y_centroid)**2.d0 + (z_cb(k) - z_centroid)**2.d0) &
                        !                                                                          /(2.d0*0.005d0)))

                        !bump in R(x)
                        !       q_prim_vf(adv_idx%end + 1)%sf(i, j, k) = q_prim_vf(adv_idx%end + 1)%sf(i, j, k)* &
                        !                                               (1d0 + 0.2d0*exp(-1d0* &
                        !                                                                               ((x_cb(i) - x_centroid)**2.d0 + (y_cb(j) - y_centroid)**2.d0 + (z_cb(k) - z_centroid)**2.d0) &
!                                                                                  /(2.d0*0.005d0)))

                        !reassign density
                        !         q_prim_vf(1)%sf(i, j, k) = &
                        !             (((q_prim_vf(E_idx)%sf(i, j, k) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                        !            rhoref*(1d0 - q_prim_vf(E_idx + 1)%sf(i, j, k))

                        ! ================================================================================

                        ! Constant x-velocity in cylindrical grid ========================================
!                        q_prim_vf(cont_idx%beg )%sf(i,j,k) = 1d0
!                        q_prim_vf(cont_idx%end )%sf(i,j,k) = 0d0
!                        q_prim_vf(mom_idx%beg  )%sf(i,j,k) = 0d0
!                        q_prim_vf(mom_idx%beg+1)%sf(i,j,k) = COS(z_cc(k))
!                        q_prim_vf(mom_idx%end  )%sf(i,j,k) = -SIN(z_cc(k))
!                        q_prim_vf(E_idx        )%sf(i,j,k) = 1d0
!                        q_prim_vf(adv_idx%beg  )%sf(i,j,k) = 1d0
                        ! ================================================================================

                        ! Couette flow in cylindrical grid ===============================================
                        !q_prim_vf(cont_idx%beg )%sf(i,j,k) = 1d0
                        !q_prim_vf(cont_idx%end )%sf(i,j,k) = 0d0
                        !q_prim_vf(mom_idx%beg  )%sf(i,j,k) = 0d0
                        !q_prim_vf(mom_idx%beg+1)%sf(i,j,k) = y_cc(j)*COS(z_cc(k))*SIN(z_cc(k))
                        !q_prim_vf(mom_idx%end  )%sf(i,j,k) = -y_cc(j)*SIN(z_cc(k))**2
                        !q_prim_vf(E_idx        )%sf(i,j,k) = 1d0
                        !q_prim_vf(adv_idx%beg  )%sf(i,j,k) = 1d0
                        ! ================================================================================

                    end if

                end do
            end do
        end do

    end subroutine s_3D_analytical ! ---------------------------------------

    !>      This patch generates the shape of the spherical harmonics
        !!      as a perturbation to a perfect sphere
        !!      @param patch_id is the patch identifier
    subroutine s_spherical_harmonic(patch_id) ! ----------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j, k !< generic loop iterators

        complex(kind(0d0)) :: cmplx_i = (0d0, 1d0)
        complex(kind(0d0)) :: H

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
        eta = 1d0

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

                        if (epsilon == 1d0) then
                            if (beta == 0d0) then
                                H = 5d-1*sqrt(3d0/pi)*cos(sph_phi)
                            elseif (beta == 1d0) then
                                H = -5d-1*sqrt(3d0/(2d0*pi))*exp(cmplx_i*z_cc(k))*sin(sph_phi)
                            end if
                        elseif (epsilon == 2d0) then
                            if (beta == 0d0) then
                                H = 25d-2*sqrt(5d0/pi)*(3d0*cos(sph_phi)**2 - 1d0)
                            elseif (beta == 1d0) then
                                H = -5d-1*sqrt(15d0/(2d0*pi))*exp(cmplx_i*z_cc(k))*sin(sph_phi)*cos(sph_phi)
                            elseif (beta == 2d0) then
                                H = 25d-2*sqrt(15d0/(2d0*pi))*exp(2d0*cmplx_i*z_cc(k))*sin(sph_phi)**2
                            end if
                        elseif (epsilon == 3d0) then
                            if (beta == 0d0) then
                                H = 25d-2*sqrt(7d0/pi)*(5d0*cos(sph_phi)**3d0 - 3d0*cos(sph_phi))
                            elseif (beta == 1d0) then
                                H = -125d-3*sqrt(21d0/pi)*exp(cmplx_i*z_cc(k))*sin(sph_phi)* &
                                    (5d0*cos(sph_phi)**2 - 1d0)
                            elseif (beta == 2d0) then
                                H = 25d-2*sqrt(105d0/(2d0*pi))*exp(2d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**2*cos(sph_phi)
                            elseif (beta == 3d0) then
                                H = -125d-3*sqrt(35d0/pi)*exp(3d0*cmplx_i*z_cc(k))*sin(sph_phi)**3d0
                            end if
                        elseif (epsilon == 4d0) then
                            if (beta == 0d0) then
                                H = 3d0/16d0*sqrt(1d0/pi)*(35d0*cos(sph_phi)**4d0 - &
                                                           3d1*cos(sph_phi)**2 + 3d0)
                            elseif (beta == 1d0) then
                                H = -3d0/8d0*sqrt(5d0/pi)*exp(cmplx_i*z_cc(k))* &
                                    sin(sph_phi)*(7d0*cos(sph_phi)**3d0 - 3d0*cos(sph_phi))
                            elseif (beta == 2d0) then
                                H = 3d0/8d0*sqrt(5d0/(2d0*pi))*exp(2d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**2*(7d0*cos(sph_phi)**2 - 1d0)
                            elseif (beta == 3d0) then
                                H = -3d0/8d0*sqrt(35d0/pi)*exp(3d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**3d0*cos(sph_phi)
                            elseif (beta == 4d0) then
                                H = 3d0/16d0*sqrt(35d0/(2d0*pi))*exp(4d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**4d0
                            end if
                        elseif (epsilon == 5d0) then
                            if (beta == 0d0) then
                                H = 1d0/16d0*sqrt(11d0/pi)*(63d0*cos(sph_phi)**5d0 - &
                                                            7d1*cos(sph_phi)**3d0 + 15d0*cos(sph_phi))
                            elseif (beta == 1d0) then
                                H = -1d0/16d0*sqrt(165d0/(2d0*pi))*exp(cmplx_i*z_cc(k))* &
                                    sin(sph_phi)*(21d0*cos(sph_phi)**4d0 - 14d0*cos(sph_phi)**2 + 1d0)
                            elseif (beta == 2d0) then
                                H = 125d-3*sqrt(1155d0/(2d0*pi))*exp(2d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**2*(3d0*cos(sph_phi)**3d0 - cos(sph_phi))
                            elseif (beta == 3d0) then
                                H = -1d0/32d0*sqrt(385d0/pi)*exp(3d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**3d0*(9d0*cos(sph_phi)**2 - 1d0)
                            elseif (beta == 4d0) then
                                H = 3d0/16d0*sqrt(385d0/(2d0*pi))*exp(4d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**4d0*cos(sph_phi)
                            elseif (beta == 5d0) then
                                H = -3d0/32d0*sqrt(77d0/pi)*exp(5d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**5d0
                            end if
                        end if

                        q_prim_vf(adv_idx%beg)%sf(i, j, k) = 1d0 - abs(real(H, kind(0d0)))

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
    subroutine s_sphere(patch_id) ! ----------------------------------------

        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j, k !< generic loop iterators

        real(kind(0d0)) :: radius_pressure, pressure_bubble, pressure_inf !<
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
        eta = 1d0

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
                                    - radius))*(-0.5d0) + 0.5d0

                    end if

                    if (((x_cc(i) - x_centroid)**2 &
                         + (cart_y - y_centroid)**2 &
                         + (cart_z - z_centroid)**2 <= radius**2 &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

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
    subroutine s_cuboid(patch_id) ! ----------------------------------------

        integer, intent(IN) :: patch_id

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
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y
        z_boundary%beg = z_centroid - 0.5d0*length_z
        z_boundary%end = z_centroid + 0.5d0*length_z

        ! Since the cuboidal patch does not allow for its boundaries to get
        ! smoothed out, the pseudo volume fraction is set to 1 to make sure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

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

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

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
    subroutine s_cylinder(patch_id) ! --------------------------------------

        integer, intent(IN) :: patch_id

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
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y
        z_boundary%beg = z_centroid - 0.5d0*length_z
        z_boundary%end = z_centroid + 0.5d0*length_z

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smearing of the cylindrical patch's boundary is enabled.
        eta = 1d0

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
                                        - radius))*(-0.5d0) + 0.5d0
                        elseif (length_y /= dflt_real) then
                            eta = tanh(smooth_coeff/min(dx, dz)* &
                                       (sqrt((x_cc(i) - x_centroid)**2 &
                                             + (cart_z - z_centroid)**2) &
                                        - radius))*(-0.5d0) + 0.5d0
                        else
                            eta = tanh(smooth_coeff/min(dx, dy)* &
                                       (sqrt((x_cc(i) - x_centroid)**2 &
                                             + (cart_y - y_centroid)**2) &
                                        - radius))*(-0.5d0) + 0.5d0
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

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

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
    subroutine s_sweep_plane(patch_id) ! -----------------------------------

        integer, intent(IN) :: patch_id

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
        eta = 1d0

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
                        eta = 5d-1 + 5d-1*tanh(smooth_coeff/min(dx, dy, dz) &
                                               *(a*x_cc(i) + &
                                                 b*cart_y + &
                                                 c*cart_z + d) &
                                               /sqrt(a**2 + b**2 + c**2))
                    end if

                    if ((a*x_cc(i) + b*cart_y + c*cart_z + d >= 0d0 &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

                    end if

                end do
            end do
        end do

    end subroutine s_sweep_plane ! -----------------------------------------

end module m_create_patches
