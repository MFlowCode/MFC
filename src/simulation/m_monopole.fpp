!>
!! @file m_monopole.f90
!! @brief Contains module m_monopole

#:include 'macros.fpp'

!> @brief The module contains the subroutines used to create a monopole pressure source term
module m_monopole

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_bubbles              !< Bubble dynamic routines

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper_basic         !< Functions to compare floating point numbers
    ! ==========================================================================
    implicit none
    private; public :: s_initialize_monopole_module, s_monopole_calculations

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(integer, dimension(:), pulse, support)
    !$acc declare link(pulse, support)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :), loc_mono)
    !$acc declare link(loc_mono)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), mag, length, wavelength, frequency, gauss_sigma_dist, gauss_sigma_time, npulse, dir, delay)
    !$acc declare link(mag, length, wavelength, frequency, gauss_sigma_dist, gauss_sigma_time, npulse, dir, delay)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), foc_length, aperture)
    !$acc declare link(foc_length, aperture)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), element_spacing_angle, element_polygon_ratio, rotate_angle)
    !$acc declare link(element_spacing_angle, element_polygon_ratio, rotate_angle)

    @:CRAY_DECLARE_GLOBAL(integer, dimension(:), num_elements, element_on)
    !$acc declare link(num_elements, element_on)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), mono_mass_src, mono_e_src)
    !$acc declare link(mono_mass_src, mono_e_src)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), mono_mom_src)
    !$acc declare link(mono_mom_src)

#else
    integer, allocatable, dimension(:) :: pulse, support
    !$acc declare create(pulse, support)

    real(kind(0d0)), allocatable, target, dimension(:, :) :: loc_mono
    !$acc declare create(loc_mono)

    real(kind(0d0)), allocatable, dimension(:) :: mag, length, wavelength, frequency, gauss_sigma_dist, gauss_sigma_time, npulse, dir, delay
    !$acc declare create(mag, length, wavelength, frequency, gauss_sigma_dist, gauss_sigma_time, npulse, dir, delay)

    real(kind(0d0)), allocatable, dimension(:) :: foc_length, aperture
    !$acc declare create(foc_length, aperture)

    real(kind(0d0)), allocatable, dimension(:) :: element_spacing_angle, element_polygon_ratio, rotate_angle
    !$acc declare create(element_spacing_angle, element_polygon_ratio, rotate_angle)

    integer, allocatable, dimension(:) :: num_elements, element_on
    !$acc declare create(num_elements, element_on)

    !> @name Monopole source terms
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mono_mass_src, mono_e_src
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: mono_mom_src
    !> @}
    !$acc declare create(mono_mass_src, mono_e_src, mono_mom_src)

#endif

contains

    subroutine s_initialize_monopole_module
        integer :: i, j !< generic loop variables

        @:ALLOCATE_GLOBAL(loc_mono(1:3, 1:num_mono), mag(1:num_mono), support(1:num_mono), length(1:num_mono), wavelength(1:num_mono), frequency(1:num_mono), gauss_sigma_dist(1:num_mono), gauss_sigma_time(1:num_mono), foc_length(1:num_mono), aperture(1:num_mono), npulse(1:num_mono), pulse(1:num_mono), dir(1:num_mono), delay(1:num_mono), element_polygon_ratio(1:num_mono), rotate_angle(1:num_mono), element_spacing_angle(1:num_mono), num_elements(1:num_mono), element_on(1:num_mono))
        do i = 1, num_mono
            do j = 1, 3
                loc_mono(j, i) = mono(i)%loc(j)
            end do
            mag(i) = mono(i)%mag
            support(i) = mono(i)%support
            length(i) = mono(i)%length
            wavelength(i) = mono(i)%wavelength ! TODO EITHER wavelength OR frequency
            frequency(i) = mono(i)%frequency
            gauss_sigma_dist(i) = mono(i)%gauss_sigma_dist ! TODO EITHER gauss_sigma_dist OR gauss_sigma_time
            gauss_sigma_time(i) = mono(i)%gauss_sigma_time
            foc_length(i) = mono(i)%foc_length
            aperture(i) = mono(i)%aperture
            if (mono(i)%npulse == dflt_int) then
                npulse(i) = 1
            else
                npulse(i) = mono(i)%npulse
            end if
            if (mono(i)%pulse == dflt_int) then
                pulse(i) = 1
            else
                pulse(i) = mono(i)%pulse
            end if
            if (f_is_default(mono(i)%dir)) then
                dir(i) = 1d0
            else
                dir(i) = mono(i)%dir
            end if
            if (f_is_default(mono(i)%delay)) then
                delay(i) = 0d0
            else
                delay(i) = mono(i)%delay
            end if
            element_spacing_angle(i) = mono(i)%element_spacing_angle
            element_polygon_ratio(i) = mono(i)%element_polygon_ratio
            rotate_angle(i) = mono(i)%rotate_angle
            num_elements(i) = mono(i)%num_elements
            element_on(i) = mono(i)%element_on
        end do
        !$acc update device(loc_mono, mag, support, length, wavelength, frequency, gauss_sigma_dist, gauss_sigma_time, foc_length, aperture, npulse, pulse, dir, delay, element_polygon_ratio, rotate_angle, element_spacing_angle, num_elements, element_on)

        @:ALLOCATE_GLOBAL(mono_mass_src(0:m, 0:n, 0:p))
        @:ALLOCATE_GLOBAL(mono_mom_src(1:num_dims, 0:m, 0:n, 0:p))
        @:ALLOCATE_GLOBAL(mono_E_src(0:m, 0:n, 0:p))

    end subroutine

    subroutine s_monopole_calculations(q_cons_vf, &
                                       q_prim_vf, t_step, id, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf !<
        !! This variable contains the WENO-reconstructed values of the cell-average
        !! conservative variables, which are located in q_cons_vf, at cell-interior
        !! Gaussian quadrature points (QP).

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf !<
        !! The primitive variables at cell-interior Gaussian quadrature points. These
        !! are calculated from the conservative variables and gradient magnitude (GM)
        !! of the volume fractions, q_cons_qp and gm_alpha_qp, respectively.

        integer, intent(in) :: t_step, id

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(kind(0d0)) :: myR, myV, alf, myP, myRho, R2Vav

        integer :: i, j, k, l, q, ii !< generic loop variables
        integer :: term_index

        real(kind(0d0)), dimension(num_fluids) :: myalpha_rho, myalpha

        real(kind(0d0)) :: n_tait, B_tait, angle, ratio_x_r, ratio_y_r, ratio_z_r

        integer :: ndirs

        real(kind(0d0)) :: the_time, sound
        real(kind(0d0)) :: s1, s2

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mono_mass_src(j, k, l) = 0d0; mono_mom_src(1, j, k, l) = 0d0; mono_e_src(j, k, l) = 0d0; 
                    if (n > 0) then
                        mono_mom_src(2, j, k, l) = 0d0
                    end if
                    if (p > 0) then
                        mono_mom_src(3, j, k, l) = 0d0
                    end if
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present) private(myalpha_rho, myalpha)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    !$acc loop seq
                    do q = 1, num_mono

                        the_time = t_step*dt ! TODO move outside of the loop?
                        if (the_time < delay(q) .and. pulse(q) == 1) cycle ! need to explicitly set delay to 0 if default; now it hinges on the fact default is a negative value

                        !$acc loop seq
                        do ii = 1, num_fluids
                            myalpha_rho(ii) = q_cons_vf(ii)%sf(j, k, l)
                            myalpha(ii) = q_cons_vf(advxb + ii - 1)%sf(j, k, l)
                        end do

                        myRho = 0d0
                        n_tait = 0d0
                        B_tait = 0d0

                        if (bubbles) then
                            if (mpp_lim .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do ii = 1, num_fluids
                                    myRho = myRho + myalpha_rho(ii)
                                    n_tait = n_tait + myalpha(ii)*gammas(ii)
                                    B_tait = B_tait + myalpha(ii)*pi_infs(ii)
                                end do
                            else if (num_fluids > 2) then
                                !$acc loop seq
                                do ii = 1, num_fluids - 1
                                    myRho = myRho + myalpha_rho(ii)
                                    n_tait = n_tait + myalpha(ii)*gammas(ii)
                                    B_tait = B_tait + myalpha(ii)*pi_infs(ii)
                                end do
                            else
                                myRho = myalpha_rho(1)
                                n_tait = gammas(1)
                                B_tait = pi_infs(1)
                            end if
                        else
                            !$acc loop seq
                            do ii = 1, num_fluids
                                myRho = myRho + myalpha_rho(ii)
                                n_tait = n_tait + myalpha(ii)*gammas(ii)
                                B_tait = B_tait + myalpha(ii)*pi_infs(ii)
                            end do
                        end if
                        n_tait = 1d0/n_tait + 1d0 ! The usual little 'gamma'

                        sound = dsqrt(n_tait*(q_prim_vf(E_idx)%sf(j, k, l) + ((n_tait - 1d0)/n_tait)*B_tait)/myRho)

                        term_index = 2

                        angle = 0d0
                        ratio_x_r = 0d0
                        ratio_y_r = 0d0
                        ratio_z_r = 0d0

                        s2 = f_g(the_time, sound, q, term_index)* &
                             f_delta(j, k, l, loc_mono(:, q), length(q), q, angle, ratio_x_r, ratio_y_r, ratio_z_r)

                        if (support(q) == 5 .or. support(q) == 7) then
                            term_index = 1
                            s1 = f_g(the_time, sound, q, term_index)* &
                                 f_delta(j, k, l, loc_mono(:, q), length(q), q, angle, ratio_x_r, ratio_y_r, ratio_z_r)
                            mono_mass_src(j, k, l) = mono_mass_src(j, k, l) + s1
                        else
                            mono_mass_src(j, k, l) = mono_mass_src(j, k, l) + s2/sound
                        end if

                        if (n == 0) then ! 1D
                            if (dir(q) < 0d0) then ! Left-going wave
                                mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) - s2
                            else ! Right-going wave
                                mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2
                            end if
                        else if (p == 0) then ! 2D
                            if (.not. f_is_default(dir(q))) then
                                if (support(q) == 5 .or. support(q) == 7) then
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(angle)
                                    mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*sin(angle)
                                else
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(dir(q))
                                    mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*sin(dir(q))
                                end if
                            end if
                        else ! 3D
                            if (.not. f_is_default(dir(q))) then
                                if (support(q) == 5 .or. support(q) == 7) then
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*ratio_x_r
                                    mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*ratio_y_r
                                    mono_mom_src(3, j, k, l) = mono_mom_src(3, j, k, l) + s2*ratio_z_r
                                else
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(dir(q))
                                    mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*sin(dir(q))
                                end if
                            end if
                        end if

                        if (model_eqns /= 4) then
                            if (any(support(q) == (/5, 7/))) then
                                mono_E_src(j, k, l) = mono_E_src(j, k, l) + s1*sound**2d0/(n_tait - 1d0)
                            else
                                mono_E_src(j, k, l) = mono_E_src(j, k, l) + s2*sound/(n_tait - 1d0)
                            end if
                        end if

                    end do
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    !$acc loop seq
                    do q = contxb, contxe
                        rhs_vf(q)%sf(j, k, l) = rhs_vf(q)%sf(j, k, l) + mono_mass_src(j, k, l)
                    end do
                    !$acc loop seq
                    do q = momxb, momxe
                        rhs_vf(q)%sf(j, k, l) = rhs_vf(q)%sf(j, k, l) + mono_mom_src(q - contxe, j, k, l)
                    end do
                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + mono_e_src(j, k, l)
                end do
            end do
        end do

    end subroutine

    !> This function gives the temporally varying amplitude of the pulse
        !! @param the_time Simulation time
        !! @param sos Sound speed
    function f_g(the_time, sos, nm, term_index)
        !$acc routine seq
        real(kind(0d0)), intent(in) :: the_time, sos
        integer, intent(in) :: nm, term_index
        real(kind(0d0)) :: omega ! angular frequency
        real(kind(0d0)) :: f_g

        real(kind(0d0)) :: foc_length_factor ! Scale amplitude with radius for spherical support
        ! i.e. Spherical support -> 1/r scaling; Cylindrical support -> 1/sqrt(r) [^-0.5 -> ^-0.85] (empirical correction)

        if (n == 0) then
            foc_length_factor = 1d0
        else if (p == 0 .and. (.not. cyl_coord)) then ! 2D axisymmetric case is physically 3D
            foc_length_factor = foc_length(nm)**(-0.85d0); 
        else
            foc_length_factor = 1/foc_length(nm); 
        end if

        f_g = 0d0

        if (pulse(nm) == 1) then
            ! Sine wave
            if (f_is_default(frequency(nm)) .and. f_is_default(wavelength(nm))) then
                wavelength(nm) = length(nm) ! For CI test - TODO remove, add frequency to test case files, and regenerate tests
            end if
            if (f_is_default(frequency(nm))) then
                frequency(nm) = sos/wavelength(nm) ! TODO CHANGE
            end if
            if ((the_time - delay(nm))*frequency(nm) > npulse(nm)) return
            omega = 2d0*pi*frequency(nm)
            if (term_index == 1) then
                f_g = mag(nm)*sin((the_time - delay(nm))*omega)/sos &
                      + foc_length_factor*mag(nm)*(cos((the_time - delay(nm))*omega) - 1d0)/omega
            else
                f_g = mag(nm)*sin((the_time - delay(nm))*omega)
            end if
        else if (pulse(nm) == 2) then
            ! Gaussian pulse
            if (f_is_default(gauss_sigma_time(nm))) then
                gauss_sigma_time(nm) = sos/gauss_sigma_dist(nm) ! TODO CHANGE
            end if
            if (term_index == 1) then
                f_g = mag(nm)*dexp(-0.5d0*((the_time - delay(nm))**2d0)/(gauss_sigma_time(nm)**2d0))/sos - &
                      foc_length_factor*mag(nm)*dsqrt(pi/2)*gauss_sigma_time(nm)* &
                      (erf((the_time - delay(nm))/(dsqrt(2d0)*gauss_sigma_time(nm))) + 1)
            else
                f_g = mag(nm)*dexp(-0.5d0*((the_time - delay(nm))**2d0)/(gauss_sigma_time(nm)**2d0))
            end if
        else if (pulse(nm) == 3) then
            ! Square wave
            if (f_is_default(frequency(nm))) then
                frequency(nm) = sos/wavelength(nm) ! TODO CHANGE
            end if
            if ((the_time - delay(nm))*frequency(nm) > npulse(nm)) return
            omega = 2d0*pi*frequency(nm)
            f_g = mag(nm)*sign(1d0, sin((the_time - delay(nm))*omega))
        end if

    end function f_g

    !> This function give the spatial support of the acoustic source
        !! @param j First coordinate-direction location index
        !! @param k Second coordinate-direction location index
        !! @param l Third coordinate-direction location index
        !! @param mono_loc Nominal source term location
        !! @param mono_leng Length of source term in space
    function f_delta(j, k, l, mono_loc, mono_leng, nm, angle, ratio_x_r, ratio_y_r, ratio_z_r)

        !$acc routine seq
        real(kind(0d0)), dimension(3), intent(in) :: mono_loc
        integer, intent(in) :: nm
        real(kind(0d0)), intent(in) :: mono_leng
        integer, intent(in) :: j, k, l

        integer :: q
        real(kind(0d0)) :: h, hx, hy, hz
        real(kind(0d0)) :: hx_cyl, hy_cyl, hz_cyl
        real(kind(0d0)) :: hxnew, hynew
        real(kind(0d0)) :: hxnew_cyl, hynew_cyl
        real(kind(0d0)) :: sig
        real(kind(0d0)) :: f_delta
        real(kind(0d0)) :: angle, ratio_x_r, ratio_y_r, ratio_z_r

        integer :: elem, elem_min, elem_max
        real(kind(0d0)) :: current_angle, angle_per_elem, angle_half_aperture
        real(kind(0d0)) :: angle_min, angle_max
        real(kind(0d0)) :: aperture_element_3D, poly_side_length, angle_elem, alpha
        real(kind(0d0)) :: dist_interp_to_elem_center
        real(kind(0d0)) :: hy_elem_center, hz_elem_center
        real(kind(0d0)) :: x2, y2, z2, x3, y3, z3, C, f, R, theta
        real(kind(0d0)) :: norm

        if (n == 0) then
            sig = dx(j)
        else if (p == 0) then
            sig = maxval((/dx(j), dy(k)/))
        else
            sig = maxval((/dx(j), dy(k), dz(l)/))
        end if
        sig = sig*acoustic_spatial_support_width

        if (n == 0) then ! 1D
            if (support(nm) == 1) then
                ! 1D delta function
                hx = abs(mono_loc(1) - x_cc(j))

                f_delta = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                          dexp(-0.5d0*(hx/(sig/2d0))**2d0)
            else if (support(nm) == 0) then
                ! Support for all x
                f_delta = 1d0
            end if
        else if (p == 0) then ! 2D
            hx = mono_loc(1) - x_cc(j)
            hy = mono_loc(2) - y_cc(k)
            if (support(nm) == 1) then
                ! 2D delta function
                sig = mono_leng/20.d0 ! For CI test - TODO remove & regenerate tests
                h = dsqrt(hx**2d0 + hy**2d0)

                f_delta = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                          dexp(-0.5d0*((h/(sig/2d0))**2d0))
            else if (support(nm) == 2) then
                !only support for y \pm some value
                if (abs(hy) < length(nm)) then
                    f_delta = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                              dexp(-0.5d0*(hx/(sig/2d0))**2d0)
                else
                    f_delta = 0d0
                end if
            else if (support(nm) == 3) then
                ! Only support along some line
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)

                ! Rotate actual point by -theta
                hxnew = cos(dir(nm))*hx + sin(dir(nm))*hy
                hynew = -1d0*sin(dir(nm))*hx + cos(dir(nm))*hy
                if (abs(hynew) < mono_loc(3)/2d0) then
                    f_delta = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)
                else
                    f_delta = 0d0
                end if
            else if (support(nm) == 4) then
                ! Support for all y
                f_delta = 1d0/(dsqrt(2d0*pi)*sig)* &
                          dexp(-0.5d0*(hx/sig)**2d0)
            else if (support(nm) == 5) then
                ! Support along transducer in 2D or 2D axisymmetric coordinate
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)

                current_angle = -atan(hy/(foc_length(nm) - hx))
                angle_half_aperture = asin((aperture(nm)/2d0)/(foc_length(nm)))

                f_delta = 0d0
                if (abs(current_angle) < angle_half_aperture .and. hx < foc_length(nm)) then
                    hxnew = foc_length(nm) - dsqrt(hy**2d0 + (foc_length(nm) - hx)**2d0)
                    f_delta = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)
                    angle = -atan(hy/(foc_length(nm) - hx))
                end if
            else if (support(nm) == 7) then
                ! Support along transducer array in 2D
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)

                current_angle = -atan(hy/(foc_length(nm) - hx))
                angle_half_aperture = asin((aperture(nm)/2d0)/(foc_length(nm)))
                angle_per_elem = (2d0*angle_half_aperture - (num_elements(nm) - 1d0)*element_spacing_angle(nm))/num_elements(nm)
                hxnew = foc_length(nm) - dsqrt(hy**2d0 + (foc_length(nm) - hx)**2d0)

                if (element_on(nm) == 0) then ! Full transducer
                    elem_min = 1
                    elem_max = num_elements(nm)
                else ! Transducer element specified
                    elem_min = element_on(nm)
                    elem_max = element_on(nm)
                end if

                f_delta = 0d0 ! If not affected by any element
                do elem = elem_min, elem_max
                    angle_max = angle_half_aperture - (element_spacing_angle(nm) + angle_per_elem)*(elem - 1d0)
                    angle_min = angle_max - angle_per_elem

                    if (current_angle > angle_min .and. current_angle < angle_max .and. hx < foc_length(nm)) then
                        f_delta = dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)/(dsqrt(2d0*pi)*sig/2d0)
                        angle = current_angle
                        exit ! Assume elements don't overlap
                    end if
                end do
            end if

        else ! 3D
            hx = x_cc(j) - mono_loc(1)
            hy = y_cc(k) - mono_loc(2)
            hz = z_cc(l) - mono_loc(3)
            if (support(nm) == 3) then

                ! Rotate actual point by -theta
                hxnew = cos(dir(nm))*hx + sin(dir(nm))*hy
                hynew = -1d0*sin(dir(nm))*hx + cos(dir(nm))*hy

                if (abs(hynew) < length(nm)/2. .and. &
                    abs(hz) < length(nm)/2.) then
                    f_delta = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)
                else
                    f_delta = 0d0
                end if
            else if (support(nm) == 4) then
                ! Support for all x,y
                f_delta = 1d0/(dsqrt(2d0*pi)*sig)* &
                          dexp(-0.5d0*(hz/sig)**2d0)
            else if (support(nm) == 5) then
                ! Support along transducer in 3D
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)
                hz = z_cc(l) - mono_loc(3)

                current_angle = -atan(dsqrt(hy**2 + hz**2)/(foc_length(nm) - hx))
                angle_half_aperture = asin((aperture(nm)/2d0)/(foc_length(nm)))

                f_delta = 0d0
                if (abs(current_angle) < angle_half_aperture .and. hx < foc_length(nm)) then

                    hxnew = foc_length(nm) - dsqrt(hy**2d0 + hz**2d0 + (foc_length(nm) - hx)**2d0)
                    f_delta = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)

                    norm = dsqrt(hy**2d0 + hz**2d0 + (foc_length(nm) - hx)**2d0)
                    ratio_x_r = -(hx - foc_length(nm))/norm
                    ratio_y_r = -hy/norm
                    ratio_z_r = -hz/norm
                end if
            else if (support(nm) == -1) then ! Not working - disabled for now
                ! Support for transducer in 3D cylindrical coordinate
                sig = maxval((/dx(j), dy(k)*sin(dz(l)), dz(l)*cos(dz(l))/))*acoustic_spatial_support_width
                hx_cyl = x_cc(j) - mono_loc(1)
                hy_cyl = y_cc(k)*sin(z_cc(l)) - mono_loc(2)
                hz_cyl = y_cc(k)*cos(z_cc(l)) - mono_loc(3)

                ! Rotate actual point by -theta
                hxnew_cyl = cos(dir(nm))*hx_cyl + sin(dir(nm))*hy_cyl
                hynew_cyl = -1d0*sin(dir(nm))*hx_cyl + cos(dir(nm))*hy_cyl

                f_delta = 0d0
                if (abs(hynew_cyl) < length(nm)/2. .and. &
                    abs(hz_cyl) < length(nm)/2.) then
                    f_delta = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                              dexp(-0.5d0*(hxnew_cyl/(sig/2d0))**2d0)
                end if
            else if (support(nm) == 7) then
                ! Support along transducer ring in 3D

                poly_side_length = aperture(nm)*sin(pi/num_elements(nm))
                aperture_element_3D = poly_side_length*element_polygon_ratio(nm)

                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)
                hz = z_cc(l) - mono_loc(3)

                f = foc_length(nm)
                R = aperture(nm)/2d0

                if (element_on(nm) == 0) then ! Full transducer
                    elem_min = 1
                    elem_max = num_elements(nm)
                else ! Transducer element specified
                    elem_min = element_on(nm)
                    elem_max = element_on(nm)
                end if

                f_delta = 0d0 ! If not affected by any element
                do elem = elem_min, elem_max
                    angle_elem = 2d0*pi*real(elem, kind(0d0))/real(num_elements(nm), kind(0d0)) + rotate_angle(nm)

                    ! Point 2 is the elem center
                    x2 = f - dsqrt(f**2 - R**2); 
                    y2 = R*cos(angle_elem)
                    z2 = R*sin(angle_elem)

                    ! Construct a plane normal to the line from the focal point to the elem center,
                    ! Point 3 is the intercept of the plane and the line from the focal point to the current location
                    C = f**2d0/((hx - f)*(x2 - f) + hy*y2 + hz*z2) ! Intermediate step constant without much physical meaning
                    x3 = C*(hx - f) + f
                    y3 = C*hy
                    z3 = C*hz

                    dist_interp_to_elem_center = dsqrt((x2 - x3)**2d0 + (y2 - y3)**2d0 + (z2 - z3)**2d0)
                    if ((dist_interp_to_elem_center < aperture_element_3D/2d0) .and. (hx < f)) then
                        hxnew = dsqrt((x3 - hx)**2d0 + (y3 - hy)**2d0 + (z3 - hz)**2d0)
                        ! Note: hxnew here is not related to hx
                        !       It's the distance of the point of interest to the element plane
                        f_delta = dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)/(dsqrt(2d0*pi)*sig/2d0)
                        norm = dsqrt(hy**2d0 + hz**2d0 + (foc_length(nm) - hx)**2d0)
                        ratio_x_r = -(hx - foc_length(nm))/norm
                        ratio_y_r = -hy/norm
                        ratio_z_r = -hz/norm
                    end if
                end do

            end if
        end if

    end function f_delta

end module m_monopole
