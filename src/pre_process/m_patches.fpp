!>
!! @file m_patches.fpp
!! @brief Contains module m_patches

#:include 'case.fpp'
#:include '1dHardcodedIC.fpp'
#:include '2dHardcodedIC.fpp'
#:include '3dHardcodedIC.fpp'
#:include 'macros.fpp'

module m_patches

    ! Dependencies =============================================================
    use m_model                 ! Subroutine(s) related to STL files

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_helper

    use m_mpi_common

    use m_assign_variables

    use m_mpi_common
    ! ==========================================================================

    implicit none

    private; public :: s_line_segment, &
 s_spiral, &
 s_circle, &
 s_airfoil, &
 s_3D_airfoil, &
 s_varcircle, &
 s_3dvarcircle, &
 s_ellipse, &
 s_ellipsoid, &
 s_rectangle, &
 s_sweep_line, &
 s_2D_TaylorGreen_vortex, &
 s_1D_analytical, &
 s_1d_bubble_pulse, &
 s_2D_analytical, &
 s_3D_analytical, &
 s_spherical_harmonic, &
 s_sphere, &
 s_cuboid, &
 s_cylinder, &
 s_sweep_plane, &
 s_model

    real(kind(0d0)) :: x_centroid, y_centroid, z_centroid
    real(kind(0d0)) :: length_x, length_y, length_z

    integer :: smooth_patch_id
    real(kind(0d0)) :: smooth_coeff !<
    !! These variables are analogous in both meaning and use to the similarly
    !! named components in the ic_patch_parameters type (see m_derived_types.f90
    !! for additional details). They are employed as a means to more concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    real(kind(0d0)) :: eta !<
    !! In the case that smoothing of patch boundaries is enabled and the boundary
    !! between two adjacent patches is to be smeared out, this variable's purpose
    !! is to act as a pseudo volume fraction to indicate the contribution of each
    !! patch toward the composition of a cell's fluid state.

    real(kind(0d0)) :: cart_y, cart_z
    real(kind(0d0)) :: sph_phi !<
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates

    type(bounds_info) :: x_boundary, y_boundary, z_boundary  !<
    !! These variables combine the centroid and length parameters associated with
    !! a particular patch to yield the locations of the patch boundaries in the
    !! x-, y- and z-coordinate directions. They are used as a means to concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    character(len=5) :: istr ! string to store int to string result for error checking

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

        real(kind(0d0)) :: pi_inf, gamma, lit_gamma

        integer :: i, j, k !< Generic loop operators

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

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0, &
                                                        eta, q_prim_vf, patch_id_fp)

                @:analytical()

                ! Updating the patch identities bookkeeping variable
                if (1d0 - eta < 1d-16) patch_id_fp(i, 0, 0) = patch_id

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
                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()

                    ! Updating the patch identities bookkeeping variable
                    if (1d0 - eta < 1d-16) patch_id_fp(i, j, 0) = patch_id
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
    subroutine s_circle(patch_id, patch_id_fp, q_prim_vf, ib) ! ----------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        logical, intent(IN) :: ib   !< True if this patch is an immersed boundary
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: radius

        integer :: i, j, k !< Generic loop iterators

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information

        if (.not. ib) then
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            radius = patch_icpp(patch_id)%radius
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff = patch_icpp(patch_id)%smooth_coeff
        else
            x_centroid = patch_ib(patch_id)%x_centroid
            y_centroid = patch_ib(patch_id)%y_centroid
            radius = patch_ib(patch_id)%radius
        end if

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

                if (.not. ib .and. patch_icpp(patch_id)%smoothen) then

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
                    (.not. ib .and. patch_id_fp(i, j, 0) == smooth_patch_id)) &
                    then
                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()
                end if

            end do
        end do

    end subroutine s_circle ! ----------------------------------------------

    subroutine s_airfoil(patch_id, patch_id_fp, q_prim_vf, ib)

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        logical, intent(IN) :: ib   !< True if this patch is an immersed boundary
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: x0, y0, f, x_act, y_act, ca, pa, ma, ta, theta, xa, ya, yt, xu, yu, xl, yl, xc, yc, dycdxc, sin_c, cos_c
        integer :: i, j, k, l
        integer :: Np1, Np2

        if (.not. ib) return
        x0 = patch_ib(patch_id)%x_centroid
        y0 = patch_ib(patch_id)%y_centroid
        ca = patch_ib(patch_id)%c
        pa = patch_ib(patch_id)%p
        ma = patch_ib(patch_id)%m
        ta = patch_ib(patch_id)%t
        theta = pi*patch_ib(patch_id)%theta/180d0

        Np1 = int((pa*ca/dx)*20)
        Np2 = int(((ca - pa*ca)/dx)*20)
        Np = Np1 + Np2 + 1

        allocate (airfoil_grid_u(1:Np))
        allocate (airfoil_grid_l(1:Np))

        airfoil_grid_u(1)%x = x0
        airfoil_grid_u(1)%y = y0

        airfoil_grid_l(1)%x = x0
        airfoil_grid_l(1)%y = y0

        eta = 1d0

        do i = 1, Np1 + Np2 - 1
            if (i <= Np1) then
                xc = x0 + i*(pa*ca/Np1)
                xa = (xc - x0)/ca
                yc = (ma/pa**2)*(2*pa*xa - xa**2)
                dycdxc = (2*ma/pa**2)*(pa - xa)
            else
                xc = x0 + pa*ca + (i - Np1)*((ca - pa*ca)/Np2)
                xa = (xc - x0)/ca
                yc = (ma/(1 - pa)**2)*(1 - 2*pa + 2*pa*xa - xa**2)
                dycdxc = (2*ma/(1 - pa)**2)*(pa - xa)
            end if

            yt = (5d0*ta)*(0.2969*xa**0.5d0 - 0.126*xa - 0.3516*xa**2d0 + 0.2843*xa**3 - 0.1015*xa**4)
            sin_c = dycdxc/(1 + dycdxc**2)**0.5
            cos_c = 1/(1 + dycdxc**2)**0.5

            xu = xa - yt*sin_c
            yu = yc + yt*cos_c

            xl = xa + yt*sin_c
            yl = yc - yt*cos_c

            xu = xu*ca + x0
            yu = yu*ca + y0

            xl = xl*ca + x0
            yl = yl*ca + y0

            airfoil_grid_u(i + 1)%x = xu
            airfoil_grid_u(i + 1)%y = yu

            airfoil_grid_l(i + 1)%x = xl
            airfoil_grid_l(i + 1)%y = yl

        end do

        airfoil_grid_u(Np)%x = x0 + ca
        airfoil_grid_u(Np)%y = y0

        airfoil_grid_l(Np)%x = x0 + ca
        airfoil_grid_l(Np)%y = y0

        do j = 0, n
            do i = 0, m

                if (patch_ib(patch_id)%theta /= dflt_real) then
                    x_act = (x_cc(i) - x0)*cos(theta) - (y_cc(j) - y0)*sin(theta) + x0
                    y_act = (x_cc(i) - x0)*sin(theta) + (y_cc(j) - y0)*cos(theta) + y0
                else
                    x_act = x_cc(i)
                    y_act = y_cc(j)
                end if

                if (x_act >= x0 .and. x_act <= x0 + ca) then
                    xa = (x_act - x0)/ca
                    if (xa <= pa) then
                        yc = (ma/pa**2)*(2*pa*xa - xa**2)
                        dycdxc = (2*ma/pa**2)*(pa - xa)
                    else
                        yc = (ma/(1 - pa)**2)*(1 - 2*pa + 2*pa*xa - xa**2)
                        dycdxc = (2*ma/(1 - pa)**2)*(pa - xa)
                    end if
                    if (y_act >= y0) then
                        k = 1
                        do while (airfoil_grid_u(k)%x < x_act)
                            k = k + 1
                        end do
                        if (airfoil_grid_u(k)%x == x_act) then
                            if (y_act <= airfoil_grid_u(k)%y) then
                                !!IB
                                !call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                !eta, q_prim_vf, patch_id_fp)
                                patch_id_fp(i, j, 0) = patch_id
                            end if
                        else
                            f = (airfoil_grid_u(k)%x - x_act)/(airfoil_grid_u(k)%x - airfoil_grid_u(k - 1)%x)
                            if (y_act <= ((1d0 - f)*airfoil_grid_u(k)%y + f*airfoil_grid_u(k - 1)%y)) then
                                !!IB
                                !call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                !eta, q_prim_vf, patch_id_fp)
                                patch_id_fp(i, j, 0) = patch_id
                            end if
                        end if
                    else
                        k = 1
                        do while (airfoil_grid_l(k)%x < x_act)
                            k = k + 1
                        end do
                        if (airfoil_grid_l(k)%x == x_act) then
                            if (y_act >= airfoil_grid_l(k)%y) then
                                !!IB
                                !call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                !eta, q_prim_vf, patch_id_fp)
                                patch_id_fp(i, j, 0) = patch_id
                            end if
                        else
                            f = (airfoil_grid_l(k)%x - x_act)/(airfoil_grid_l(k)%x - airfoil_grid_l(k - 1)%x)

                            if (y_act >= ((1d0 - f)*airfoil_grid_l(k)%y + f*airfoil_grid_l(k - 1)%y)) then
                                   !!IB
                                !call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                !eta, q_prim_vf, patch_id_fp)
                                patch_id_fp(i, j, 0) = patch_id
                            end if
                        end if
                    end if
                end if
            end do
        end do

        if (patch_ib(patch_id)%theta /= dflt_real) then
            do i = 1, Np
                airfoil_grid_l(i)%x = (airfoil_grid_l(i)%x - x0)*cos(theta) + (airfoil_grid_l(i)%y - y0)*sin(theta) + x0
                airfoil_grid_l(i)%y = -1d0*(airfoil_grid_l(i)%x - x0)*sin(theta) + (airfoil_grid_l(i)%y - y0)*cos(theta) + y0

                airfoil_grid_u(i)%x = (airfoil_grid_u(i)%x - x0)*cos(theta) + (airfoil_grid_u(i)%y - y0)*sin(theta) + x0
                airfoil_grid_u(i)%y = -1d0*(airfoil_grid_u(i)%x - x0)*sin(theta) + (airfoil_grid_u(i)%y - y0)*cos(theta) + y0
            end do
        end if

    end subroutine s_airfoil

    subroutine s_3D_airfoil(patch_id, patch_id_fp, q_prim_vf, ib)

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        logical, intent(IN) :: ib   !< True if this patch is an immersed boundary
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: x0, y0, z0, lz, z_max, z_min, f, x_act, y_act, ca, pa, ma, ta, theta, xa, ya, yt, xu, yu, xl, yl, xc, yc, dycdxc, sin_c, cos_c
        integer :: i, j, k, l
        integer :: Np1, Np2

        if (.not. ib) return
        x0 = patch_ib(patch_id)%x_centroid
        y0 = patch_ib(patch_id)%y_centroid
        z0 = patch_ib(patch_id)%z_centroid
        lz = patch_ib(patch_id)%length_z
        ca = patch_ib(patch_id)%c
        pa = patch_ib(patch_id)%p
        ma = patch_ib(patch_id)%m
        ta = patch_ib(patch_id)%t
        theta = pi*patch_ib(patch_id)%theta/180d0

        Np1 = int((pa*ca/dx)*20)
        Np2 = int(((ca - pa*ca)/dx)*20)
        Np = Np1 + Np2 + 1

        allocate (airfoil_grid_u(1:Np))
        allocate (airfoil_grid_l(1:Np))

        airfoil_grid_u(1)%x = x0
        airfoil_grid_u(1)%y = y0

        airfoil_grid_l(1)%x = x0
        airfoil_grid_l(1)%y = y0

        z_max = z0 + lz/2
        z_min = z0 - lz/2

        eta = 1d0

        do i = 1, Np1 + Np2 - 1
            if (i <= Np1) then
                xc = x0 + i*(pa*ca/Np1)
                xa = (xc - x0)/ca
                yc = (ma/pa**2)*(2*pa*xa - xa**2)
                dycdxc = (2*ma/pa**2)*(pa - xa)
            else
                xc = x0 + pa*ca + (i - Np1)*((ca - pa*ca)/Np2)
                xa = (xc - x0)/ca
                yc = (ma/(1 - pa)**2)*(1 - 2*pa + 2*pa*xa - xa**2)
                dycdxc = (2*ma/(1 - pa)**2)*(pa - xa)
            end if

            yt = (5d0*ta)*(0.2969*xa**0.5d0 - 0.126*xa - 0.3516*xa**2d0 + 0.2843*xa**3 - 0.1015*xa**4)
            sin_c = dycdxc/(1 + dycdxc**2)**0.5
            cos_c = 1/(1 + dycdxc**2)**0.5

            xu = xa - yt*sin_c
            yu = yc + yt*cos_c

            xl = xa + yt*sin_c
            yl = yc - yt*cos_c

            xu = xu*ca + x0
            yu = yu*ca + y0

            xl = xl*ca + x0
            yl = yl*ca + y0

            airfoil_grid_u(i + 1)%x = xu
            airfoil_grid_u(i + 1)%y = yu

            airfoil_grid_l(i + 1)%x = xl
            airfoil_grid_l(i + 1)%y = yl

        end do

        airfoil_grid_u(Np)%x = x0 + ca
        airfoil_grid_u(Np)%y = y0

        airfoil_grid_l(Np)%x = x0 + ca
        airfoil_grid_l(Np)%y = y0

        do l = 0, p
            if (z_cc(l) >= z_min .and. z_cc(l) <= z_max) then
                do j = 0, n
                    do i = 0, m

                        if (patch_ib(patch_id)%theta /= dflt_real) then
                            x_act = (x_cc(i) - x0)*cos(theta) - (y_cc(j) - y0)*sin(theta) + x0
                            y_act = (x_cc(i) - x0)*sin(theta) + (y_cc(j) - y0)*cos(theta) + y0
                        else
                            x_act = x_cc(i)
                            y_act = y_cc(j)
                        end if

                        if (x_act >= x0 .and. x_act <= x0 + ca) then
                            xa = (x_act - x0)/ca
                            if (xa <= pa) then
                                yc = (ma/pa**2)*(2*pa*xa - xa**2)
                                dycdxc = (2*ma/pa**2)*(pa - xa)
                            else
                                yc = (ma/(1 - pa)**2)*(1 - 2*pa + 2*pa*xa - xa**2)
                                dycdxc = (2*ma/(1 - pa)**2)*(pa - xa)
                            end if
                            if (y_act >= y0) then
                                k = 1
                                do while (airfoil_grid_u(k)%x < x_act)
                                    k = k + 1
                                end do
                                if (airfoil_grid_u(k)%x == x_act) then
                                    if (y_act <= airfoil_grid_u(k)%y) then
                                        !!IB
                                        !call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                        !eta, q_prim_vf, patch_id_fp)
                                        patch_id_fp(i, j, l) = patch_id
                                    end if
                                else
                                    f = (airfoil_grid_u(k)%x - x_act)/(airfoil_grid_u(k)%x - airfoil_grid_u(k - 1)%x)
                                    if (y_act <= ((1d0 - f)*airfoil_grid_u(k)%y + f*airfoil_grid_u(k - 1)%y)) then
                                        !!IB
                                        !call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                        !eta, q_prim_vf, patch_id_fp)
                                        patch_id_fp(i, j, l) = patch_id
                                    end if
                                end if
                            else
                                k = 1
                                do while (airfoil_grid_l(k)%x < x_act)
                                    k = k + 1
                                end do
                                if (airfoil_grid_l(k)%x == x_act) then
                                    if (y_act >= airfoil_grid_l(k)%y) then
                                        !!IB
                                        !call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                        !eta, q_prim_vf, patch_id_fp)
                                        patch_id_fp(i, j, l) = patch_id
                                    end if
                                else
                                    f = (airfoil_grid_l(k)%x - x_act)/(airfoil_grid_l(k)%x - airfoil_grid_l(k - 1)%x)

                                    if (y_act >= ((1d0 - f)*airfoil_grid_l(k)%y + f*airfoil_grid_l(k - 1)%y)) then
                                           !!IB
                                        !call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                        !eta, q_prim_vf, patch_id_fp)
                                        patch_id_fp(i, j, l) = patch_id
                                    end if
                                end if
                            end if
                        end if
                    end do
                end do
            end if
        end do

        if (patch_ib(patch_id)%theta /= dflt_real) then
            do i = 1, Np
                airfoil_grid_l(i)%x = (airfoil_grid_l(i)%x - x0)*cos(theta) + (airfoil_grid_l(i)%y - y0)*sin(theta) + x0
                airfoil_grid_l(i)%y = -1d0*(airfoil_grid_l(i)%x - x0)*sin(theta) + (airfoil_grid_l(i)%y - y0)*cos(theta) + y0

                airfoil_grid_u(i)%x = (airfoil_grid_u(i)%x - x0)*cos(theta) + (airfoil_grid_u(i)%y - y0)*sin(theta) + x0
                airfoil_grid_u(i)%y = -1d0*(airfoil_grid_u(i)%x - x0)*sin(theta) + (airfoil_grid_u(i)%y - y0)*cos(theta) + y0
            end do
        end if

    end subroutine s_3D_airfoil

    !>             The varcircle patch is a 2D geometry that may be used
        !!             . It  generatres an annulus
        !! @param patch_id is the patch identifier
    subroutine s_varcircle(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: radius

        ! Generic loop iterators
        integer :: i, j, k

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

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()

                    ! Updating the patch identities bookkeeping variable
                    if (1d0 - eta < 1d-16) patch_id_fp(i, j, 0) = patch_id

                    q_prim_vf(alf_idx)%sf(i, j, 0) = patch_icpp(patch_id)%alpha(1)* &
                                                     dexp(-0.5d0*((myr - radius)**2.d0)/(thickness/3.d0)**2.d0)
                end if

            end do
        end do

    end subroutine s_varcircle ! ----------------------------------------------

    subroutine s_3dvarcircle(patch_id, patch_id_fp, q_prim_vf) ! ----------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: radius

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

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()

                        ! Updating the patch identities bookkeeping variable
                        if (1d0 - eta < 1d-16) patch_id_fp(i, j, k) = patch_id

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
    subroutine s_ellipse(patch_id, patch_id_fp, q_prim_vf) ! ---------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: a, b

        integer :: i, j, k !< Generic loop operators

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

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()

                    ! Updating the patch identities bookkeeping variable
                    if (1d0 - eta < 1d-16) patch_id_fp(i, j, 0) = patch_id
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
        real(kind(0d0)) :: a, b, c

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

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()

                        ! Updating the patch identities bookkeeping variable
                        if (1d0 - eta < 1d-16) patch_id_fp(i, j, k) = patch_id
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
    subroutine s_rectangle(patch_id, patch_id_fp, q_prim_vf, ib) ! -------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< Equation of state parameters
        logical :: ib !< True if this patch is an immersed boundary

        integer :: i, j, k !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the rectangle's centroid and length information
        if (.not. ib) then
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            length_x = patch_icpp(patch_id)%length_x
            length_y = patch_icpp(patch_id)%length_y
        else
            x_centroid = patch_ib(patch_id)%x_centroid
            y_centroid = patch_ib(patch_id)%y_centroid
            length_x = patch_ib(patch_id)%length_x
            length_y = patch_ib(patch_id)%length_y
        end if

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

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    if ((q_prim_vf(1)%sf(i, j, 0) < 1.e-10) .and. (model_eqns == 4)) then
                        !zero density, reassign according to Tait EOS
                        q_prim_vf(1)%sf(i, j, 0) = &
                            (((q_prim_vf(E_idx)%sf(i, j, 0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                            rhoref*(1d0 - q_prim_vf(alf_idx)%sf(i, j, 0))
                    end if

                    ! Updating the patch identities bookkeeping variable
                    if (1d0 - eta < 1d-16) patch_id_fp(i, j, 0) = patch_id

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
        real(kind(0d0)) :: a, b, c

        integer :: i, j, k !< Generic loop operators

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
                    call s_assign_patch_primitive_variables(patch_id, i, j, 0, &
                                                            eta, q_prim_vf, patch_id_fp)

                    @:analytical()

                    ! Updating the patch identities bookkeeping variable
                    if (1d0 - eta < 1d-16) patch_id_fp(i, j, 0) = patch_id
                end if

            end do
        end do

    end subroutine s_sweep_line ! ------------------------------------------

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

        integer :: i, j, k !< generic loop iterators

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

                    @:analytical()

                    ! Updating the patch identities bookkeeping variable
                    if (1d0 - eta < 1d-16) patch_id_fp(i, j, 0) = patch_id

                    ! Assign Parameters =========================================================
                    q_prim_vf(mom_idx%beg)%sf(i, j, 0) = U0*sin(x_cc(i)/L0)*cos(y_cc(j)/L0)
                    q_prim_vf(mom_idx%end)%sf(i, j, 0) = -U0*cos(x_cc(i)/L0)*sin(y_cc(j)/L0)
                    q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(patch_id)%pres + (cos(2*x_cc(i))/L0 + &
                                                                                cos(2*y_cc(j))/L0)* &
                                                   (q_prim_vf(1)%sf(i, j, 0)*U0*U0)/16
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
        real(kind(0d0)) :: a, b, c, d, pi_inf, gamma, lit_gamma

        ! Generic loop iterators
        integer :: i, j, k

        @:Hardcoded1DVariables()

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

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0, &
                                                        eta, q_prim_vf, patch_id_fp)

                @:Hardcoded1D()

                ! Updating the patch identities bookkeeping variable
                if (1d0 - eta < 1d-16) patch_id_fp(i, 0, 0) = patch_id

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
        real(kind(0d0)) :: fac, a, b, c, d, pi_inf, gamma, lit_gamma

        ! Generic loop iterators
        integer :: i, j, k

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

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0, &
                                                        eta, q_prim_vf, patch_id_fp)

                @:analytical()

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

        real(kind(0d0)) :: a, b, c, d !< placeholderrs for the cell boundary values
        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< equation of state parameters
        real(kind(0d0)) :: l, U0 !< Taylor Green Vortex parameters

        integer :: i, j, k !< generic loop iterators

        @:Hardcoded2DVariables()

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
        l = 1d0
        U0 = 0.1
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

                    @:Hardcoded2D()
                    ! Updating the patch identities bookkeeping variable
                    if (1d0 - eta < 1d-16) patch_id_fp(i, j, 0) = patch_id

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
        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< equation of state parameters

        integer :: i, j, k !< generic loop iterators

        @:Hardcoded3DVariables()

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

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:Hardcoded3D()

                        ! Updating the patch identities bookkeeping variable
                        if (1d0 - eta < 1d-16) patch_id_fp(i, j, k) = patch_id

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

        real(kind(0d0)) :: epsilon, beta
        real(kind(0d0)) :: radius

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
    subroutine s_sphere(patch_id, patch_id_fp, q_prim_vf, ib) ! ----------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        logical, intent(IN) :: ib   !< True if this patch is an immersed boundary
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: radius

        ! Generic loop iterators
        integer :: i, j, k !< generic loop iterators

        real(kind(0d0)) :: radius_pressure, pressure_bubble, pressure_inf !<
            !! Variables to initialize the pressure field that corresponds to the
            !! bubble-collapse test case found in Tiwari et al. (2013)

        ! Transferring spherical patch's radius, centroid, smoothing patch
        ! identity and smoothing coefficient information
        if (.not. ib) then
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            z_centroid = patch_icpp(patch_id)%z_centroid
            radius = patch_icpp(patch_id)%radius
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff = patch_icpp(patch_id)%smooth_coeff
        else
            x_centroid = patch_ib(patch_id)%x_centroid
            y_centroid = patch_ib(patch_id)%y_centroid
            z_centroid = patch_ib(patch_id)%z_centroid
            radius = patch_ib(patch_id)%radius
        end if

        print *, x_centroid, y_centroid, z_centroid, radius

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

                    if (.not. ib .and. patch_icpp(patch_id)%smoothen) then

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
                        (.not. ib .and. patch_id_fp(i, j, k) == smooth_patch_id)) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()
                    end if

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

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()

                        ! Updating the patch identities bookkeeping variable
                        if (1d0 - eta < 1d-16) patch_id_fp(i, j, k) = patch_id

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
    subroutine s_cylinder(patch_id, patch_id_fp, q_prim_vf, ib) ! --------------------------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        logical, intent(IN) :: ib   !< True if this patch is an immersed boundary
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: radius

        integer :: i, j, k !< Generic loop iterators

        ! Transferring the cylindrical patch's centroid, length, radius,
        ! smoothing patch identity and smoothing coefficient information

        if (.not. ib) then
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            z_centroid = patch_icpp(patch_id)%z_centroid
            length_x = patch_icpp(patch_id)%length_x
            length_y = patch_icpp(patch_id)%length_y
            length_z = patch_icpp(patch_id)%length_z
            radius = patch_icpp(patch_id)%radius
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff = patch_icpp(patch_id)%smooth_coeff
        else
            x_centroid = patch_ib(patch_id)%x_centroid
            y_centroid = patch_ib(patch_id)%y_centroid
            z_centroid = patch_ib(patch_id)%z_centroid
            length_x = patch_ib(patch_id)%length_x
            length_y = patch_ib(patch_id)%length_y
            length_z = patch_ib(patch_id)%length_z
            radius = patch_ib(patch_id)%radius
        end if

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

                    if (.not. ib .and. patch_icpp(patch_id)%smoothen) then

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
                        (.not. ib .and. patch_id_fp(i, j, k) == smooth_patch_id)) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()

                        ! Updating the patch identities bookkeeping variable
                        if (1d0 - eta < 1d-16) patch_id_fp(i, j, k) = patch_id
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
        real(kind(0d0)) :: a, b, c, d

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

                        call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                                eta, q_prim_vf, patch_id_fp)

                        @:analytical()

                        ! Updating the patch identities bookkeeping variable
                        if (1d0 - eta < 1d-16) patch_id_fp(i, j, k) = patch_id
                    end if

                end do
            end do
        end do

    end subroutine s_sweep_plane ! -----------------------------------------

    !> The STL patch is a 2/3D geometry that is imported from an STL file.
    !! @param patch_id is the patch identifier
    subroutine s_model(patch_id, patch_id_fp, q_prim_vf) ! ---------------------

        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf

        integer :: i, j, k !< Generic loop iterators

        type(t_bbox) :: bbox
        type(t_model) :: model
        type(ic_model_parameters) :: params

        t_vec3 :: point

        real(kind(0d0)) :: grid_mm(1:3, 1:2)

        integer :: cell_num
        integer :: ncells

        t_mat4x4 :: transform

        if (proc_rank == 0) then
            print *, " * Reading model: "//trim(patch_icpp(patch_id)%model%filepath)
        end if
        model = f_model_read(patch_icpp(patch_id)%model%filepath)

        if (proc_rank == 0) then
            print *, " * Transforming model..."
        end if
        transform = f_create_transform_matrix(patch_icpp(patch_id)%model)
        call s_transform_model(model, transform)

        bbox = f_create_bbox(model)

        if (proc_rank == 0) then
            write (*, "(A, 3(2X, F20.10))") "    > Model:  Min:", bbox%min(1:3)
            write (*, "(A, 3(2X, F20.10))") "    >         Cen:", (bbox%min(1:3) + bbox%max(1:3))/2d0
            write (*, "(A, 3(2X, F20.10))") "    >         Max:", bbox%max(1:3)

            !call s_model_write("__out__.stl", model)
            !call s_model_write("__out__.obj", model)

            grid_mm(1, :) = (/minval(x_cc) - 0d5*dx, maxval(x_cc) + 0d5*dx/)
            grid_mm(2, :) = (/minval(y_cc) - 0d5*dy, maxval(y_cc) + 0d5*dy/)

            if (p > 0) then
                grid_mm(3, :) = (/minval(z_cc) - 0d5*dz, maxval(z_cc) + 0d5*dz/)
            else
                grid_mm(3, :) = (/0d0, 0d0/)
            end if

            write (*, "(A, 3(2X, F20.10))") "    > Domain: Min:", grid_mm(:, 1)
            write (*, "(A, 3(2X, F20.10))") "    >         Cen:", (grid_mm(:, 1) + grid_mm(:, 2))/2d0
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

                    point = (/x_cc(i), y_cc(j), 0d0/)
                    if (p > 0) then
                        point(3) = z_cc(k)
                    end if

                    if (grid_geometry == 3) then
                        point = f_convert_cyl_to_cart(point)
                    end if

                    eta = f_model_is_inside(model, point, (/dx, dy, dz/), patch_icpp(patch_id)%model%spc)

                    call s_assign_patch_primitive_variables(patch_id, i, j, k, &
                                                            eta, q_prim_vf, patch_id_fp)

                    ! Note: Should probably use *eta* to compute primitive variables
                    ! if defining them analytically.
                    @:analytical()

                end do; end do; end do

        if (proc_rank == 0) then
            print *, ""
            print *, " * Cleaning up..."
        end if

        call s_model_free(model)

    end subroutine s_model ! ---------------------------------------------------

    subroutine s_convert_cylindrical_to_cartesian_coord(cyl_y, cyl_z)
        !$acc routine seq

        real(kind(0d0)), intent(IN) :: cyl_y, cyl_z

        cart_y = cyl_y*sin(cyl_z)
        cart_z = cyl_y*cos(cyl_z)

    end subroutine s_convert_cylindrical_to_cartesian_coord ! --------------

    function f_convert_cyl_to_cart(cyl) result(cart)

        !$acc routine seq

        t_vec3, intent(in) :: cyl
        t_vec3 :: cart

        cart = (/cyl(1), &
                 cyl(2)*sin(cyl(3)), &
                 cyl(2)*cos(cyl(3))/)

    end function f_convert_cyl_to_cart

    subroutine s_convert_cylindrical_to_spherical_coord(cyl_x, cyl_y)
        !$acc routine seq

        real(kind(0d0)), intent(IN) :: cyl_x, cyl_y

        sph_phi = atan(cyl_y/cyl_x)

    end subroutine s_convert_cylindrical_to_spherical_coord ! --------------

    !> Archimedes spiral function
    !! @param myth Angle
    !! @param offset Thickness
    !! @param a Starting position
    function f_r(myth, offset, a)
        !$acc routine seq
        real(kind(0d0)), intent(IN) :: myth, offset, a
        real(kind(0d0)) :: b
        real(kind(0d0)) :: f_r

        !r(th) = a + b*th

        b = 2.d0*a/(2.d0*pi)
        f_r = a + b*myth + offset
    end function f_r

end module m_patches
