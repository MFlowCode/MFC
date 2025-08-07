#:include 'macros.fpp'

!>
!! @file m_helper.f90
!! @brief Contains module m_helper

module m_helper

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use ieee_arithmetic        !< For checking NaN

    implicit none

    private; 
    public :: s_comp_n_from_prim, &
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
              s_print_2D_array, &
              f_xor, &
              f_logical_to_int, &
              unassociated_legendre, &
              associated_legendre, &
              spherical_harmonic_func, &
              double_factorial, &
              factorial, &
              f_cut_on, &
              f_cut_off, &
              s_downsample_data, &
              s_upsample_data

contains

    !> Computes the bubble number density n from the primitive variables
        !! @param vftmp is the void fraction
        !! @param Rtmp is the  bubble radii
        !! @param ntmp is the output number bubble density
    pure subroutine s_comp_n_from_prim(vftmp, Rtmp, ntmp, weights)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: vftmp
        real(wp), dimension(nb), intent(in) :: Rtmp
        real(wp), intent(out) :: ntmp
        real(wp), dimension(nb), intent(in) :: weights

        real(wp) :: R3

        R3 = dot_product(weights, Rtmp**3._wp)
        ntmp = (3._wp/(4._wp*pi))*vftmp/R3

    end subroutine s_comp_n_from_prim

    pure subroutine s_comp_n_from_cons(vftmp, nRtmp, ntmp, weights)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: vftmp
        real(wp), dimension(nb), intent(in) :: nRtmp
        real(wp), intent(out) :: ntmp
        real(wp), dimension(nb), intent(in) :: weights

        real(wp) :: nR3

        nR3 = dot_product(weights, nRtmp**3._wp)
        ntmp = sqrt((4._wp*pi/3._wp)*nR3/vftmp)

    end subroutine s_comp_n_from_cons

    impure subroutine s_print_2D_array(A, div)

        real(wp), dimension(:, :), intent(in) :: A
        real(wp), optional, intent(in) :: div

        integer :: i, j
        integer :: local_m, local_n
        real(wp) :: c

        local_m = size(A, 1)
        local_n = size(A, 2)

        if (present(div)) then
            c = div
        else
            c = 1._wp
        end if

        print *, local_m, local_n

        do i = 1, local_m
            do j = 1, local_n
                write (*, fmt="(F12.4)", advance="no") A(i, j)/c
            end do
            write (*, fmt="(A1)") " "
        end do
        write (*, fmt="(A1)") " "

    end subroutine s_print_2D_array

    !> Initializes non-polydisperse bubble modeling
    impure subroutine s_initialize_nonpoly

        integer :: ir
        real(wp) :: rhol0, pl0, uu, D_m, temp, omega_ref
        real(wp), dimension(Nb) :: chi_vw0, cp_m0, k_m0, rho_m0, x_vw

        real(wp), parameter :: k_poly = 1._wp !<
            !! polytropic index used to compute isothermal natural frequency

        real(wp), parameter :: Ru = 8314._wp !<
            !! universal gas constant

        rhol0 = rhoref
        pl0 = pref
#ifdef MFC_SIMULATION
        @:ALLOCATE(pb0(nb), mass_n0(nb), mass_v0(nb), Pe_T(nb))
        @:ALLOCATE(k_n(nb), k_v(nb), omegaN(nb))
        @:ALLOCATE(Re_trans_T(nb), Re_trans_c(nb), Im_trans_T(nb), Im_trans_c(nb))
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
        if (thermal == 2) gamma_m = 1._wp

        temp = 293.15_wp
        D_m = 0.242e-4_wp
        uu = sqrt(pl0/rhol0)

        omega_ref = 3._wp*k_poly*Ca + 2._wp*(3._wp*k_poly - 1._wp)/Web

            !!! thermal properties !!!
        ! gas constants
        R_n = Ru/M_n
        R_v = Ru/M_v
        ! phi_vn & phi_nv (phi_nn = phi_vv = 1)
        phi_vn = (1._wp + sqrt(mu_v/mu_n)*(M_n/M_v)**(0.25_wp))**2 &
                 /(sqrt(8._wp)*sqrt(1._wp + M_v/M_n))
        phi_nv = (1._wp + sqrt(mu_n/mu_v)*(M_v/M_n)**(0.25_wp))**2 &
                 /(sqrt(8._wp)*sqrt(1._wp + M_n/M_v))
        ! internal bubble pressure
        pb0(:) = pl0 + 2._wp*ss/(R0ref*R0(:))

        ! mass fraction of vapor
        chi_vw0 = 1._wp/(1._wp + R_v/R_n*(pb0/pv - 1._wp))
        ! specific heat for gas/vapor mixture
        cp_m0 = chi_vw0*R_v*gamma_v/(gamma_v - 1._wp) &
                + (1._wp - chi_vw0)*R_n*gamma_n/(gamma_n - 1._wp)
        ! mole fraction of vapor
        x_vw = M_n*chi_vw0/(M_v + (M_n - M_v)*chi_vw0)
        ! thermal conductivity for gas/vapor mixture
        k_m0 = x_vw*k_v/(x_vw + (1._wp - x_vw)*phi_vn) &
               + (1._wp - x_vw)*k_n/(x_vw*phi_nv + 1._wp - x_vw)
        ! mixture density
        rho_m0 = pv/(chi_vw0*R_v*temp)

        ! mass of gas/vapor computed using dimensional quantities
        mass_n0(:) = 4._wp*(pb0(:) - pv)*pi/(3._wp*R_n*temp*rhol0)*R0(:)**3
        mass_v0(:) = 4._wp*pv*pi/(3._wp*R_v*temp*rhol0)*R0(:)**3
        ! Peclet numbers
        Pe_T(:) = rho_m0*cp_m0(:)*uu*R0ref/k_m0(:)
        Pe_c = uu*R0ref/D_m

        Tw = temp

        ! nondimensional properties
        !if(.not. qbmm) then
        R_n = rhol0*R_n*temp/pl0
        R_v = rhol0*R_v*temp/pl0
        k_n(:) = k_n(:)/k_m0(:)
        k_v(:) = k_v(:)/k_m0(:)
        pb0 = pb0/pl0
        pv = pv/pl0
        Tw = 1._wp
        pl0 = 1._wp

        rhoref = 1._wp
        pref = 1._wp
        !end if

        ! natural frequencies
        omegaN(:) = sqrt(3._wp*k_poly*Ca + 2._wp*(3._wp*k_poly - 1._wp)/(Web*R0))/R0
        do ir = 1, Nb
            call s_transcoeff(omegaN(ir)*R0(ir), Pe_T(ir)*R0(ir), &
                              Re_trans_T(ir), Im_trans_T(ir))
            call s_transcoeff(omegaN(ir)*R0(ir), Pe_c*R0(ir), &
                              Re_trans_c(ir), Im_trans_c(ir))
        end do
        Im_trans_T = 0._wp

    end subroutine s_initialize_nonpoly

    !> Computes the transfer coefficient for the non-polytropic bubble compression process
        !! @param omega natural frequencies
        !! @param peclet Peclet number
        !! @param Re_trans Real part of the transport coefficients
        !! @param Im_trans Imaginary part of the transport coefficients
    pure elemental subroutine s_transcoeff(omega, peclet, Re_trans, Im_trans)

        real(wp), intent(in) :: omega, peclet
        real(wp), intent(out) :: Re_trans, Im_trans

        complex(wp) :: imag, trans, c1, c2, c3

        imag = (0._wp, 1._wp)

        c1 = imag*omega*peclet
        c2 = sqrt(c1)
        c3 = (exp(c2) - exp(-c2))/(exp(c2) + exp(-c2)) ! TANH(c2)
        trans = ((c2/c3 - 1._wp)**(-1) - 3._wp/c1)**(-1) ! transfer function

        Re_trans = trans
        Im_trans = aimag(trans)

    end subroutine s_transcoeff

    pure elemental subroutine s_int_to_str(i, res)

        integer, intent(in) :: i
        character(len=*), intent(inout) :: res

        write (res, '(I0)') i
        res = trim(res)
    end subroutine s_int_to_str

    !> Computes the Simpson weights for quadrature
    subroutine s_simpson(local_weight, local_R0)

        real(wp), dimension(:), intent(inout) :: local_weight
        real(wp), dimension(:), intent(inout) :: local_R0

        integer :: ir
        real(wp) :: R0mn, R0mx, dphi, tmp, sd
        real(wp), dimension(nb) :: phi

        sd = poly_sigma
        R0mn = 0.8_wp*exp(-2.8_wp*sd)
        R0mx = 0.2_wp*exp(9.5_wp*sd) + 1._wp

        ! phi = ln( R0 ) & return R0
        do ir = 1, nb
            phi(ir) = log(R0mn) &
                      + (ir - 1._wp)*log(R0mx/R0mn)/(nb - 1._wp)
            local_R0(ir) = exp(phi(ir))
        end do
        dphi = phi(2) - phi(1)

        ! weights for quadrature using Simpson's rule
        do ir = 2, nb - 1
            ! Gaussian
            tmp = exp(-0.5_wp*(phi(ir)/sd)**2)/sqrt(2._wp*pi)/sd
            if (mod(ir, 2) == 0) then
                local_weight(ir) = tmp*4._wp*dphi/3._wp
            else
                local_weight(ir) = tmp*2._wp*dphi/3._wp
            end if
        end do
        tmp = exp(-0.5_wp*(phi(1)/sd)**2)/sqrt(2._wp*pi)/sd
        local_weight(1) = tmp*dphi/3._wp
        tmp = exp(-0.5_wp*(phi(nb)/sd)**2)/sqrt(2._wp*pi)/sd
        local_weight(nb) = tmp*dphi/3._wp
    end subroutine s_simpson

    !> This procedure computes the cross product of two vectors.
    !! @param a First vector.
    !! @param b Second vector.
    !! @return The cross product of the two vectors.
    pure function f_cross(a, b) result(c)

        real(wp), dimension(3), intent(in) :: a, b
        real(wp), dimension(3) :: c

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function f_cross

    !> This procedure swaps two real numbers.
    !! @param lhs Left-hand side.
    !! @param rhs Right-hand side.
    pure elemental subroutine s_swap(lhs, rhs)

        real(wp), intent(inout) :: lhs, rhs
        real(wp) :: ltemp

        ltemp = lhs
        lhs = rhs
        rhs = ltemp
    end subroutine s_swap

    !> This procedure creates a transformation matrix.
    !! @param  p Parameters for the transformation.
    !! @return Transformation matrix.
    pure function f_create_transform_matrix(param, center) result(out_matrix)

        type(ic_model_parameters), intent(in) :: param
        real(wp), dimension(1:3), optional, intent(in) :: center
        real(wp), dimension(1:4, 1:4) :: sc, rz, rx, ry, tr, t_back, t_to_origin, out_matrix

        sc = transpose(reshape([ &
                               param%scale(1), 0._wp, 0._wp, 0._wp, &
                               0._wp, param%scale(2), 0._wp, 0._wp, &
                               0._wp, 0._wp, param%scale(3), 0._wp, &
                               0._wp, 0._wp, 0._wp, 1._wp], shape(sc)))

        rz = transpose(reshape([ &
                               cos(param%rotate(3)), -sin(param%rotate(3)), 0._wp, 0._wp, &
                               sin(param%rotate(3)), cos(param%rotate(3)), 0._wp, 0._wp, &
                               0._wp, 0._wp, 1._wp, 0._wp, &
                               0._wp, 0._wp, 0._wp, 1._wp], shape(rz)))

        rx = transpose(reshape([ &
                               1._wp, 0._wp, 0._wp, 0._wp, &
                               0._wp, cos(param%rotate(1)), -sin(param%rotate(1)), 0._wp, &
                               0._wp, sin(param%rotate(1)), cos(param%rotate(1)), 0._wp, &
                               0._wp, 0._wp, 0._wp, 1._wp], shape(rx)))

        ry = transpose(reshape([ &
                               cos(param%rotate(2)), 0._wp, sin(param%rotate(2)), 0._wp, &
                               0._wp, 1._wp, 0._wp, 0._wp, &
                               -sin(param%rotate(2)), 0._wp, cos(param%rotate(2)), 0._wp, &
                               0._wp, 0._wp, 0._wp, 1._wp], shape(ry)))

        tr = transpose(reshape([ &
                               1._wp, 0._wp, 0._wp, param%translate(1), &
                               0._wp, 1._wp, 0._wp, param%translate(2), &
                               0._wp, 0._wp, 1._wp, param%translate(3), &
                               0._wp, 0._wp, 0._wp, 1._wp], shape(tr)))

        if (present(center)) then
            ! Translation matrix to move center to the origin
            t_to_origin = transpose(reshape([ &
                                            1._wp, 0._wp, 0._wp, -center(1), &
                                            0._wp, 1._wp, 0._wp, -center(2), &
                                            0._wp, 0._wp, 1._wp, -center(3), &
                                            0._wp, 0._wp, 0._wp, 1._wp], shape(tr)))

            ! Translation matrix to move center back to original position
            t_back = transpose(reshape([ &
                                       1._wp, 0._wp, 0._wp, center(1), &
                                       0._wp, 1._wp, 0._wp, center(2), &
                                       0._wp, 0._wp, 1._wp, center(3), &
                                       0._wp, 0._wp, 0._wp, 1._wp], shape(tr)))

            out_matrix = matmul(tr, matmul(t_back, matmul(ry, matmul(rx, matmul(rz, matmul(sc, t_to_origin))))))
        else
            out_matrix = matmul(ry, matmul(rx, rz))
        end if

    end function f_create_transform_matrix

    !> This procedure transforms a vector by a matrix.
    !! @param vec Vector to transform.
    !! @param matrix Transformation matrix.
    pure subroutine s_transform_vec(vec, matrix)

        real(wp), dimension(1:3), intent(inout) :: vec
        real(wp), dimension(1:4, 1:4), intent(in) :: matrix

        real(wp), dimension(1:4) :: tmp

        tmp = matmul(matrix, [vec(1), vec(2), vec(3), 1._wp])
        vec = tmp(1:3)

    end subroutine s_transform_vec

    !> This procedure transforms a triangle by a matrix, one vertex at a time.
    !! @param triangle Triangle to transform.
    !! @param matrix   Transformation matrix.
    pure subroutine s_transform_triangle(triangle, matrix, matrix_n)

        type(t_triangle), intent(inout) :: triangle
        real(wp), dimension(1:4, 1:4), intent(in) :: matrix, matrix_n

        integer :: i

        do i = 1, 3
            call s_transform_vec(triangle%v(i, :), matrix)
        end do

        call s_transform_vec(triangle%n(1:3), matrix_n)

    end subroutine s_transform_triangle

    !> This procedure transforms a model by a matrix, one triangle at a time.
    !! @param model  Model to transform.
    !! @param matrix Transformation matrix.
    pure subroutine s_transform_model(model, matrix, matrix_n)

        type(t_model), intent(inout) :: model
        real(wp), dimension(1:4, 1:4), intent(in) :: matrix, matrix_n

        integer :: i

        do i = 1, size(model%trs)
            call s_transform_triangle(model%trs(i), matrix, matrix_n)
        end do

    end subroutine s_transform_model

    !> This procedure creates a bounding box for a model.
    !! @param model Model to create bounding box for.
    !! @return Bounding box.
    pure function f_create_bbox(model) result(bbox)

        type(t_model), intent(in) :: model
        type(t_bbox) :: bbox

        integer :: i, j

        if (size(model%trs) == 0) then
            bbox%min = 0._wp
            bbox%max = 0._wp
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

    !> This procedure performs xor on lhs and rhs.
    !! @param lhs logical input.
    !! @param rhs other logical input.
    !! @return xored result.
    pure elemental function f_xor(lhs, rhs) result(res)

        logical, intent(in) :: lhs, rhs
        logical :: res

        res = (lhs .and. .not. rhs) .or. (.not. lhs .and. rhs)
    end function f_xor

    !> This procedure converts logical to 1 or 0.
    !! @param perdicate A Logical argument.
    !! @return 1 if .true., 0 if .false..
    pure elemental function f_logical_to_int(predicate) result(int)

        logical, intent(in) :: predicate
        integer :: int

        if (predicate) then
            int = 1
        else
            int = 0
        end if
    end function f_logical_to_int

    !> This function generates the unassociated legendre poynomials
    !! @param x is the input value
    !! @param l is the degree
    !! @return P is the unassociated legendre polynomial evaluated at x
    pure recursive function unassociated_legendre(x, l) result(result_P)

        integer, intent(in) :: l
        real(wp), intent(in) :: x
        real(wp) :: result_P

        if (l == 0) then
            result_P = 1._wp
        else if (l == 1) then
            result_P = x
        else
            result_P = ((2*l - 1)*x*unassociated_legendre(x, l - 1) - (l - 1)*unassociated_legendre(x, l - 2))/l
        end if

    end function unassociated_legendre

    !> This function calculates the spherical harmonic function evaluated at x and phi
    !! @param x is the x coordinate
    !! @param phi is the phi coordinate
    !! @param l is the degree
    !! @param m_order is the order
    !! @return Y is the spherical harmonic function evaluated at x and phi
    pure recursive function spherical_harmonic_func(x, phi, l, m_order) result(Y)

        integer, intent(in) :: l, m_order
        real(wp), intent(in) :: x, phi
        real(wp) :: Y, prefactor, local_pi

        local_pi = acos(-1._wp)
        prefactor = sqrt((2*l + 1)/(4*local_pi)*factorial(l - m_order)/factorial(l + m_order)); 
        if (m_order == 0) then
            Y = prefactor*associated_legendre(x, l, m_order); 
        elseif (m_order > 0) then
            Y = (-1._wp)**m_order*sqrt(2._wp)*prefactor*associated_legendre(x, l, m_order)*cos(m_order*phi); 
        end if

    end function spherical_harmonic_func

    !> This function generates the associated legendre polynomials evaluated
    !! at x with inputs l and m
    !! @param x is the input value
    !! @param l is the degree
    !! @param m_order is the order
    !! @return P is the associated legendre polynomial evaluated at x
    pure recursive function associated_legendre(x, l, m_order) result(result_P)

        integer, intent(in) :: l, m_order
        real(wp), intent(in) :: x
        real(wp) :: result_P

        if (m_order <= 0 .and. l <= 0) then
            result_P = 1; 
        elseif (l == 1 .and. m_order <= 0) then
            result_P = x; 
        elseif (l == 1 .and. m_order == 1) then
            result_P = -(1 - x**2)**(1._wp/2._wp); 
        elseif (m_order == l) then
            result_P = (-1)**l*double_factorial(2*l - 1)*(1 - x**2)**(l/2); 
        elseif (m_order == l - 1) then
            result_P = x*(2*l - 1)*associated_legendre(x, l - 1, l - 1); 
        else
            result_P = ((2*l - 1)*x*associated_legendre(x, l - 1, m_order) - (l + m_order - 1)*associated_legendre(x, l - 2, m_order))/(l - m_order); 
        end if

    end function associated_legendre

    !> This function calculates the double factorial value of an integer
    !! @param n_in is the input integer
    !! @return R is the double factorial value of n
    pure elemental function double_factorial(n_in) result(R_result)

        integer, intent(in) :: n_in
        integer, parameter :: int64_kind = selected_int_kind(18) ! 18 bytes for 64-bit integer
        integer(kind=int64_kind) :: R_result
        integer :: i

        R_result = product((/(i, i=n_in, 1, -2)/))

    end function double_factorial

    !> The following function calculates the factorial value of an integer
    !! @param n_in is the input integer
    !! @return R is the factorial value of n
    pure elemental function factorial(n_in) result(R_result)

        integer, intent(in) :: n_in
        integer, parameter :: int64_kind = selected_int_kind(18) ! 18 bytes for 64-bit integer
        integer(kind=int64_kind) :: R_result

        integer :: i

        R_result = product((/(i, i=n_in, 1, -1)/))

    end function factorial

    !> This function calculates a smooth cut-on function that is zero for x values
    !! smaller than zero and goes to one. It can be used for generating smooth
    !! initial conditions
    !! @param x is the input value
    !! @param eps is the smoothing parameter
    !! @return fx is the cut-on function evaluated at x
    function f_cut_on(x, eps) result(fx)

        real(wp), intent(in) :: x, eps
        real(wp) :: fx

        fx = 1 - f_gx(x/eps)/(f_gx(x/eps) + f_gx(1 - x/eps))

    end function f_cut_on

    !> This function calculates a smooth cut-off function that is one for x values
    !! smaller than zero and goes to zero. It can be used for generating smooth
    !! initial conditions
    !! @param x is the input value
    !! @param eps is the smoothing parameter
    !! @return fx is the cut-ff function evaluated at x
    function f_cut_off(x, eps) result(fx)

        real(wp), intent(in) :: x, eps
        real(wp) :: fx

        fx = f_gx(x/eps)/(f_gx(x/eps) + f_gx(1 - x/eps))

    end function f_cut_off

    !> This function is a helper function for the functions f_cut_on and f_cut_off
    !! @param x is the input value
    !! @return gx is the result
    function f_gx(x) result(gx)

        real(wp), intent(in) :: x
        real(wp) :: gx

        if (x > 0) then
            gx = exp(-1._wp/x)
        else
            gx = 0._wp
        end if

    end function f_gx

    subroutine s_downsample_data(q_cons_vf, q_cons_temp, m_ds, n_ds, p_ds, m_glb_ds, n_glb_ds, p_glb_ds)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf, q_cons_temp

        ! Down sampling variables
        integer :: i, j, k, l
        integer :: ix, iy, iz, x_id, y_id, z_id
        integer, intent(inout) :: m_ds, n_ds, p_ds, m_glb_ds, n_glb_ds, p_glb_ds

        m_ds = int((m + 1)/3) - 1
        n_ds = int((n + 1)/3) - 1
        p_ds = int((p + 1)/3) - 1

        m_glb_ds = int((m_glb + 1)/3) - 1
        n_glb_ds = int((n_glb + 1)/3) - 1
        p_glb_ds = int((p_glb + 1)/3) - 1

        do i = 1, sys_size
            $:GPU_UPDATE(host='[q_cons_vf(i)%sf]')
        end do

        do l = -1, p_ds + 1
            do k = -1, n_ds + 1
                do j = -1, m_ds + 1
                    x_id = 3*j + 1
                    y_id = 3*k + 1
                    z_id = 3*l + 1
                    do i = 1, sys_size
                        q_cons_temp(i)%sf(j, k, l) = 0

                        do iz = -1, 1
                            do iy = -1, 1
                                do ix = -1, 1
                                    q_cons_temp(i)%sf(j, k, l) = q_cons_temp(i)%sf(j, k, l) &
                                                                 + (1._wp/27._wp)*q_cons_vf(i)%sf(x_id + ix, y_id + iy, z_id + iz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine s_downsample_data

    subroutine s_upsample_data(q_cons_vf, q_cons_temp)

        type(scalar_field), intent(inout), dimension(sys_size) :: q_cons_vf, q_cons_temp
        integer :: i, j, k, l
        integer :: ix, iy, iz
        integer :: x_id, y_id, z_id
        real(wp), dimension(4) :: temp

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    do i = 1, sys_size

                        ix = int(j/3._wp)
                        iy = int(k/3._wp)
                        iz = int(l/3._wp)

                        x_id = j - int(3*ix) - 1
                        y_id = k - int(3*iy) - 1
                        z_id = l - int(3*iz) - 1

                        temp(1) = (2._wp/3._wp)*q_cons_temp(i)%sf(ix, iy, iz) + (1._wp/3._wp)*q_cons_temp(i)%sf(ix + x_id, iy, iz)
                        temp(2) = (2._wp/3._wp)*q_cons_temp(i)%sf(ix, iy + y_id, iz) + (1._wp/3._wp)*q_cons_temp(i)%sf(ix + x_id, iy + y_id, iz)
                        temp(3) = (2._wp/3._wp)*temp(1) + (1._wp/3._wp)*temp(2)

                        temp(1) = (2._wp/3._wp)*q_cons_temp(i)%sf(ix, iy, iz + z_id) + (1._wp/3._wp)*q_cons_temp(i)%sf(ix + x_id, iy, iz + z_id)
                        temp(2) = (2._wp/3._wp)*q_cons_temp(i)%sf(ix, iy + y_id, iz + z_id) + (1._wp/3._wp)*q_cons_temp(i)%sf(ix + x_id, iy + y_id, iz + z_id)
                        temp(4) = (2._wp/3._wp)*temp(1) + (1._wp/3._wp)*temp(2)

                        q_cons_vf(i)%sf(j, k, l) = (2._wp/3._wp)*temp(3) + (1._wp/3._wp)*temp(4)

                    end do
                end do
            end do
        end do

    end subroutine s_upsample_data

end module m_helper
