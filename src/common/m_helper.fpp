#:include 'macros.fpp'
!>
!! @file m_helper.f90
!! @brief Contains module m_helper

module m_helper

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_common           !< MPI modules

    use ieee_arithmetic        !< For checking NaN

    ! ==========================================================================

    implicit none

    private; public :: s_compute_finite_difference_coefficients, &
 s_comp_n_from_prim, &
 s_comp_n_from_cons, &
 s_initialize_nonpoly, &
 s_simpson, &
 s_transcoeff, &
 s_int_to_str, &
 s_transform_vec, &
 s_transform_triangle, &
 s_transform_model, &
 s_swap, &
 f_cross, &
 f_create_transform_matrix, &
 f_create_bbox, &
 s_print_2D_array

contains

    !>  The purpose of this subroutine is to compute the finite-
        !!      difference coefficients for the centered schemes utilized
        !!      in computations of first order spatial derivatives in the
        !!      s-coordinate direction. The s-coordinate direction refers
        !!      to the x-, y- or z-coordinate direction, depending on the
        !!      subroutine's inputs. Note that coefficients of up to 4th
        !!      order accuracy are available.
        !!  @param q Number of cells in the s-coordinate direction
        !!  @param s_cc Locations of the cell-centers in the s-coordinate direction
        !!  @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction
    subroutine s_compute_finite_difference_coefficients(q, s_cc, fd_coeff_s, buff_size, &
                                                        fd_number_in, fd_order_in, offset_s)

        integer :: lB, lE !< loop bounds
        integer, intent(IN) :: q
        integer, intent(IN) :: buff_size, fd_number_in, fd_order_in
        type(int_bounds_info), optional, intent(IN) :: offset_s
        real(kind(0d0)), allocatable, dimension(:, :), intent(INOUT) :: fd_coeff_s

        real(kind(0d0)), &
            dimension(-buff_size:q + buff_size), &
            intent(IN) :: s_cc

        integer :: i !< Generic loop iterator

        if (present(offset_s)) then
            lB = -offset_s%beg
            lE = q + offset_s%end
        else
            lB = 0
            lE = q
        end if

        if (allocated(fd_coeff_s)) deallocate (fd_coeff_s)
        allocate (fd_coeff_s(-fd_number_in:fd_number_in, lb:lE))

        ! Computing the 1st order finite-difference coefficients
        if (fd_order_in == 1) then
            do i = lB, lE
                fd_coeff_s(-1, i) = 0d0
                fd_coeff_s(0, i) = -1d0/(s_cc(i + 1) - s_cc(i))
                fd_coeff_s(1, i) = -fd_coeff_s(0, i)
            end do

            ! Computing the 2nd order finite-difference coefficients
        elseif (fd_order_in == 2) then
            do i = lB, lE
                fd_coeff_s(-1, i) = -1d0/(s_cc(i + 1) - s_cc(i - 1))
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
            end do

            ! Computing the 4th order finite-difference coefficients
        else
            do i = lB, lE
                fd_coeff_s(-2, i) = 1d0/(s_cc(i - 2) - 8d0*s_cc(i - 1) - s_cc(i + 2) + 8d0*s_cc(i + 1))
                fd_coeff_s(-1, i) = -8d0*fd_coeff_s(-2, i)
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
                fd_coeff_s(2, i) = -fd_coeff_s(-2, i)
            end do

        end if

    end subroutine s_compute_finite_difference_coefficients ! --------------

    !> Computes the bubble number density n from the primitive variables
        !! @param vftmp is the void fraction
        !! @param Rtmp is the  bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_prim(vftmp, Rtmp, ntmp, weights)
        !$acc routine seq
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: Rtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: R3
        real(kind(0.d0)), dimension(nb) :: weights

        R3 = dot_product(weights, Rtmp**3.d0)
        ntmp = (3.d0/(4.d0*pi))*vftmp/R3

    end subroutine s_comp_n_from_prim

    subroutine s_comp_n_from_cons(vftmp, nRtmp, ntmp, weights)
        !$acc routine seq
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: nRtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: nR3
        real(kind(0.d0)), dimension(nb) :: weights

        nR3 = dot_product(weights, nRtmp**3.d0)
        ntmp = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)
        !ntmp = (3.d0/(4.d0*pi))*0.00001

        !print *, "nbub", ntmp

    end subroutine s_comp_n_from_cons

    subroutine s_print_2D_array(A, div)

        real(kind(0d0)), dimension(:, :) :: A
        integer :: i, j
        integer :: m, n
        real :: c
        real, optional :: div

        m = size(A, 1)
        n = size(A, 2)

        if (present(div)) then
            c = div
        else
            c = 1
        end if

        print *, m, n

        do i = 1, m
            do j = 1, n
                write (*, fmt="(F12.4)", advance="no") A(i, j)/c
            end do
            write (*, fmt="(A1)") " "
        end do
        write (*, fmt="(A1)") " "

    end subroutine

    !> Initializes non-polydisperse bubble modeling
    subroutine s_initialize_nonpoly
        integer :: ir
        real(kind(0.d0)) :: rhol0
        real(kind(0.d0)) :: pl0
        real(kind(0.d0)) :: uu
        real(kind(0.d0)) :: D_m
        real(kind(0.d0)) :: temp
        real(kind(0.d0)) :: omega_ref
        real(kind(0.d0)), dimension(Nb) :: chi_vw0
        real(kind(0.d0)), dimension(Nb) :: cp_m0
        real(kind(0.d0)), dimension(Nb) :: k_m0
        real(kind(0.d0)), dimension(Nb) :: rho_m0
        real(kind(0.d0)), dimension(Nb) :: x_vw

        real(kind(0.d0)), parameter :: k_poly = 1.d0 !<
            !! polytropic index used to compute isothermal natural frequency

        real(kind(0.d0)), parameter :: Ru = 8314.d0 !<
            !! universal gas constant

        rhol0 = rhoref
        pl0 = pref
#ifdef MFC_SIMULATION
        @:ALLOCATE_GLOBAL(pb0(nb), mass_n0(nb), mass_v0(nb), Pe_T(nb))
        @:ALLOCATE_GLOBAL(k_n(nb), k_v(nb), omegaN(nb))
        @:ALLOCATE_GLOBAL(Re_trans_T(nb), Re_trans_c(nb), Im_trans_T(nb), Im_trans_c(nb))
#else
        @:ALLOCATE(pb0(nb), mass_n0(nb), mass_v0(nb), Pe_T(nb))
        @:ALLOCATE(k_n(nb), k_v(nb), omegaN(nb))
        @:ALLOCATE(Re_trans_T(nb), Re_trans_c(nb), Im_trans_T(nb), Im_trans_c(nb))
#endif

        pb0(:) = dflt_real
        mass_n0(:) = dflt_real
        mass_v0(:) = dflt_real
        Pe_T(:) = dflt_real
        omegaN(:) = dflt_real

        mul0 = fluid_pp(1)%mul0
        ss = fluid_pp(1)%ss
        pv = fluid_pp(1)%pv
        gamma_v = fluid_pp(1)%gamma_v
        M_v = fluid_pp(1)%M_v
        mu_v = fluid_pp(1)%mu_v
        k_v(:) = fluid_pp(1)%k_v

        gamma_n = fluid_pp(2)%gamma_v
        M_n = fluid_pp(2)%M_v
        mu_n = fluid_pp(2)%mu_v
        k_n(:) = fluid_pp(2)%k_v

        gamma_m = gamma_n
        if (thermal == 2) gamma_m = 1.d0

        temp = 293.15d0
        D_m = 0.242d-4
        uu = DSQRT(pl0/rhol0)

        omega_ref = 3.d0*k_poly*Ca + 2.d0*(3.d0*k_poly - 1.d0)/Web

            !!! thermal properties !!!
        ! gas constants
        R_n = Ru/M_n
        R_v = Ru/M_v
        ! phi_vn & phi_nv (phi_nn = phi_vv = 1)
        phi_vn = (1.d0 + DSQRT(mu_v/mu_n)*(M_n/M_v)**(0.25d0))**2 &
                 /(DSQRT(8.d0)*DSQRT(1.d0 + M_v/M_n))
        phi_nv = (1.d0 + DSQRT(mu_n/mu_v)*(M_v/M_n)**(0.25d0))**2 &
                 /(DSQRT(8.d0)*DSQRT(1.d0 + M_n/M_v))
        ! internal bubble pressure
        pb0 = pl0 + 2.d0*ss/(R0ref*R0)

        ! mass fraction of vapor
        chi_vw0 = 1.d0/(1.d0 + R_v/R_n*(pb0/pv - 1.d0))
        ! specific heat for gas/vapor mixture
        cp_m0 = chi_vw0*R_v*gamma_v/(gamma_v - 1.d0) &
                + (1.d0 - chi_vw0)*R_n*gamma_n/(gamma_n - 1.d0)
        ! mole fraction of vapor
        x_vw = M_n*chi_vw0/(M_v + (M_n - M_v)*chi_vw0)
        ! thermal conductivity for gas/vapor mixture
        k_m0 = x_vw*k_v/(x_vw + (1.d0 - x_vw)*phi_vn) &
               + (1.d0 - x_vw)*k_n/(x_vw*phi_nv + 1.d0 - x_vw)
        ! mixture density
        rho_m0 = pv/(chi_vw0*R_v*temp)

        ! mass of gas/vapor computed using dimensional quantities
        mass_n0 = 4.d0*(pb0 - pv)*pi/(3.d0*R_n*temp*rhol0)*R0**3
        mass_v0 = 4.d0*pv*pi/(3.d0*R_v*temp*rhol0)*R0**3
        ! Peclet numbers
        Pe_T = rho_m0*cp_m0*uu*R0ref/k_m0
        Pe_c = uu*R0ref/D_m

        Tw = temp

        ! nondimensional properties
        !if(.not. qbmm) then
        R_n = rhol0*R_n*temp/pl0
        R_v = rhol0*R_v*temp/pl0
        k_n = k_n/k_m0
        k_v = k_v/k_m0
        pb0 = pb0/pl0
        pv = pv/pl0
        Tw = 1.d0
        pl0 = 1.d0

        rhoref = 1.d0
        pref = 1.d0
        !end if

        ! natural frequencies
        omegaN = DSQRT(3.d0*k_poly*Ca + 2.d0*(3.d0*k_poly - 1.d0)/(Web*R0))/R0
        do ir = 1, Nb
            call s_transcoeff(omegaN(ir)*R0(ir), Pe_T(ir)*R0(ir), &
                              Re_trans_T(ir), Im_trans_T(ir))
            call s_transcoeff(omegaN(ir)*R0(ir), Pe_c*R0(ir), &
                              Re_trans_c(ir), Im_trans_c(ir))
        end do
        Im_trans_T = 0d0

    end subroutine s_initialize_nonpoly

    !> Computes the transfer coefficient for the non-polytropic bubble compression process
        !! @param omega natural frqeuencies
        !! @param peclet Peclet number
        !! @param Re_trans Real part of the transport coefficients
        !! @param Im_trans Imaginary part of the transport coefficients
    subroutine s_transcoeff(omega, peclet, Re_trans, Im_trans)

        real(kind(0.d0)), intent(IN) :: omega
        real(kind(0.d0)), intent(IN) :: peclet
        real(kind(0.d0)), intent(OUT) :: Re_trans
        real(kind(0.d0)), intent(OUT) :: Im_trans
        complex :: trans, c1, c2, c3
        complex :: imag = (0., 1.)
        real(kind(0.d0)) :: f_transcoeff

        c1 = imag*omega*peclet
        c2 = CSQRT(c1)
        c3 = (CEXP(c2) - CEXP(-c2))/(CEXP(c2) + CEXP(-c2)) ! TANH(c2)
        trans = ((c2/c3 - 1.d0)**(-1) - 3.d0/c1)**(-1) ! transfer function

        Re_trans = dble(trans)
        Im_trans = aimag(trans)

    end subroutine s_transcoeff

    subroutine s_int_to_str(i, res)
        character(len=*) :: res
        integer, intent(in) :: i
        write (res, '(I0)') i
        res = trim(res)
    end subroutine

    !> Computes the Simpson weights for quadrature
    subroutine s_simpson

        integer :: ir
        real(kind(0.d0)) :: R0mn
        real(kind(0.d0)) :: R0mx
        real(kind(0.d0)) :: dphi
        real(kind(0.d0)) :: tmp
        real(kind(0.d0)) :: sd
        real(kind(0.d0)), dimension(nb) :: phi

        ! nondiml. min. & max. initial radii for numerical quadrature
        !sd   = 0.05D0
        !R0mn = 0.75D0
        !R0mx = 1.3D0

        !sd   = 0.3D0
        !R0mn = 0.3D0
        !R0mx = 6.D0

        !sd   = 0.7D0
        !R0mn = 0.12D0
        !R0mx = 150.D0

        sd = poly_sigma
        R0mn = 0.8d0*DEXP(-2.8d0*sd)
        R0mx = 0.2d0*DEXP(9.5d0*sd) + 1.d0

        ! phi = ln( R0 ) & return R0
        do ir = 1, nb
            phi(ir) = DLOG(R0mn) &
                      + dble(ir - 1)*DLOG(R0mx/R0mn)/dble(nb - 1)
            R0(ir) = DEXP(phi(ir))
        end do
        dphi = phi(2) - phi(1)

        ! weights for quadrature using Simpson's rule
        do ir = 2, nb - 1
            ! Gaussian
            tmp = DEXP(-0.5d0*(phi(ir)/sd)**2)/DSQRT(2.d0*pi)/sd
            if (mod(ir, 2) == 0) then
                weight(ir) = tmp*4.d0*dphi/3.d0
            else
                weight(ir) = tmp*2.d0*dphi/3.d0
            end if
        end do
        tmp = DEXP(-0.5d0*(phi(1)/sd)**2)/DSQRT(2.d0*pi)/sd
        weight(1) = tmp*dphi/3.d0
        tmp = DEXP(-0.5d0*(phi(nb)/sd)**2)/DSQRT(2.d0*pi)/sd
        weight(nb) = tmp*dphi/3.d0
    end subroutine s_simpson

    !> This procedure computes the cross product of two vectors.
    !! @param a First vector.
    !! @param b Second vector.
    !! @return The cross product of the two vectors.
    function f_cross(a, b) result(c)
        real(kind(0d0)), dimension(3), intent(in) :: a, b
        real(kind(0d0)), dimension(3) :: c

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function f_cross

    !> This procedure swaps two real numbers.
    !! @param lhs Left-hand side.
    !! @param rhs Right-hand side.
    subroutine s_swap(lhs, rhs)
        real(kind(0d0)), intent(inout) :: lhs, rhs
        real(kind(0d0)) :: ltemp

        ltemp = lhs
        lhs = rhs
        rhs = ltemp
    end subroutine s_swap

    !> This procedure creates a transformation matrix.
    !! @param  p Parameters for the transformation.
    !! @return Transformation matrix.
    function f_create_transform_matrix(p) result(out_matrix)

        type(ic_model_parameters) :: p

        t_mat4x4 :: sc, rz, rx, ry, tr, out_matrix

        sc = transpose(reshape([ &
                               p%scale(1), 0d0, 0d0, 0d0, &
                               0d0, p%scale(2), 0d0, 0d0, &
                               0d0, 0d0, p%scale(3), 0d0, &
                               0d0, 0d0, 0d0, 1d0], shape(sc)))

        rz = transpose(reshape([ &
                               cos(p%rotate(3)), -sin(p%rotate(3)), 0d0, 0d0, &
                               sin(p%rotate(3)), cos(p%rotate(3)), 0d0, 0d0, &
                               0d0, 0d0, 1d0, 0d0, &
                               0d0, 0d0, 0d0, 1d0], shape(rz)))

        rx = transpose(reshape([ &
                               1d0, 0d0, 0d0, 0d0, &
                               0d0, cos(p%rotate(1)), -sin(p%rotate(1)), 0d0, &
                               0d0, sin(p%rotate(1)), cos(p%rotate(1)), 0d0, &
                               0d0, 0d0, 0d0, 1d0], shape(rx)))

        ry = transpose(reshape([ &
                               cos(p%rotate(2)), 0d0, sin(p%rotate(2)), 0d0, &
                               0d0, 1d0, 0d0, 0d0, &
                               -sin(p%rotate(2)), 0d0, cos(p%rotate(2)), 0d0, &
                               0d0, 0d0, 0d0, 1d0], shape(ry)))

        tr = transpose(reshape([ &
                               1d0, 0d0, 0d0, p%translate(1), &
                               0d0, 1d0, 0d0, p%translate(2), &
                               0d0, 0d0, 1d0, p%translate(3), &
                               0d0, 0d0, 0d0, 1d0], shape(tr)))

        out_matrix = matmul(tr, matmul(ry, matmul(rx, matmul(rz, sc))))

    end function f_create_transform_matrix

    !> This procedure transforms a vector by a matrix.
    !! @param vec Vector to transform.
    !! @param matrix Transformation matrix.
    subroutine s_transform_vec(vec, matrix)

        t_vec3, intent(inout) :: vec
        t_mat4x4, intent(in) :: matrix

        real(kind(0d0)), dimension(1:4) :: tmp

        tmp = matmul(matrix, [vec(1), vec(2), vec(3), 1d0])
        vec = tmp(1:3)

    end subroutine s_transform_vec

    !> This procedure transforms a triangle by a matrix, one vertex at a time.
    !! @param triangle Triangle to transform.
    !! @param matrix   Transformation matrix.
    subroutine s_transform_triangle(triangle, matrix)

        type(t_triangle), intent(inout) :: triangle
        t_mat4x4, intent(in) :: matrix

        integer :: i

        real(kind(0d0)), dimension(1:4) :: tmp

        do i = 1, 3
            call s_transform_vec(triangle%v(i, :), matrix)
        end do

    end subroutine s_transform_triangle

    !> This procedure transforms a model by a matrix, one triangle at a time.
    !! @param model  Model to transform.
    !! @param matrix Transformation matrix.
    subroutine s_transform_model(model, matrix)

        type(t_model), intent(inout) :: model
        t_mat4x4, intent(in) :: matrix

        integer :: i

        do i = 1, size(model%trs)
            call s_transform_triangle(model%trs(i), matrix)
        end do

    end subroutine s_transform_model

    !> This procedure creates a bounding box for a model.
    !! @param model Model to create bounding box for.
    !! @return Bounding box.
    function f_create_bbox(model) result(bbox)

        type(t_model), intent(in) :: model
        type(t_bbox) :: bbox

        integer :: i, j

        if (size(model%trs) == 0) then
            bbox%min = 0d0
            bbox%max = 0d0
            return
        end if

        bbox%min = model%trs(1)%v(1, :)
        bbox%max = model%trs(1)%v(1, :)

        do i = 1, size(model%trs)
            do j = 1, 3
                bbox%min = min(bbox%min, model%trs(i)%v(j, :))
                bbox%max = max(bbox%max, model%trs(i)%v(j, :))
            end do
        end do

    end function f_create_bbox

end module m_helper
