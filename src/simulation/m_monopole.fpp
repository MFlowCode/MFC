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
            npulse(i) = mono(i)%npulse
            pulse(i) = mono(i)%pulse
            dir(i) = mono(i)%dir
            element_spacing_angle(i) = mono(i)%element_spacing_angle
            element_polygon_ratio(i) = mono(i)%element_polygon_ratio
            rotate_angle(i) = mono(i)%rotate_angle
            num_elements(i) = mono(i)%num_elements
            element_on(i) = mono(i)%element_on
            if (f_is_default(mono(i)%delay)) then ! m_checker guarantees mono(i)%delay is set for pulse = 2 (Gaussian)
                delay(i) = 0d0 ! Defaults to zero for sine and square waves
            else
                delay(i) = mono(i)%delay
            end if
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

        real(kind(0d0)), dimension(sys_size) :: q_cons_elements
        real(kind(0d0)) :: q_prim_element

        integer, intent(in) :: t_step, id

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        integer :: i, j, k, l, q !< generic loop variables
        integer :: mi !< monopole index
        integer :: term_index

        real(kind(0d0)) :: sim_time, c, small_gamma
        real(kind(0d0)) :: mass_src_diff, mom_src_diff
        real(kind(0d0)) :: angle, ratio_x_r, ratio_y_r, ratio_z_r

        sim_time = t_step*dt

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

        ! Monopoles are looped through sequentially because they can have very different computational costs
        !$acc parallel loop collapse(3) gang vector default(present) private(q_cons_elements, q_prim_element, c, small_gamma, angle, ratio_x_r, ratio_y_r, ratio_z_r)
        !$acc loop seq
        do mi = 1, num_mono
            if (sim_time < delay(mi) .and. (pulse(mi) == 1 .or. pulse(mi) == 3)) cycle
            do l = 0, p
                do k = 0, n
                    do j = 0, m

                        !$acc loop seq
                        do q = 1, sys_size
                            q_cons_elements(q) = q_cons_vf(q)%sf(j, k, l)
                        end do
                        q_prim_element = q_prim_vf(E_idx)%sf(j, k, l)

                        call s_compute_speed_of_sound_monopole(q_cons_elements, q_prim_element, c, small_gamma)

                        term_index = 2

                        mom_src_diff = f_source_temporal(sim_time, c, mi, term_index)* &
                                       f_source_spatial(j, k, l, loc_mono(:, mi), mi, angle, ratio_x_r, ratio_y_r, ratio_z_r)

                        if (support(mi) == 5 .or. support(mi) == 7) then
                            term_index = 1
                            mass_src_diff = f_source_temporal(sim_time, c, mi, term_index)* &
                                            f_source_spatial(j, k, l, loc_mono(:, mi), mi, angle, ratio_x_r, ratio_y_r, ratio_z_r)
                            mono_mass_src(j, k, l) = mono_mass_src(j, k, l) + mass_src_diff
                        else
                            mono_mass_src(j, k, l) = mono_mass_src(j, k, l) + mom_src_diff/c
                        end if

                        if (n == 0) then ! 1D
                            if (dir(mi) < 0d0) then ! Left-going wave
                                mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) - mom_src_diff
                            else ! Right-going wave
                                mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + mom_src_diff
                            end if
                        elseif (p == 0) then ! 2D
                            if (.not. f_is_default(dir(mi))) then
                                if (support(mi) == 5 .or. support(mi) == 7) then
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + mom_src_diff*cos(angle)
                                    mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + mom_src_diff*sin(angle)
                                else
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + mom_src_diff*cos(dir(mi))
                                    mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + mom_src_diff*sin(dir(mi))
                                end if
                            end if
                        else ! 3D
                            if (.not. f_is_default(dir(mi))) then
                                if (support(mi) == 5 .or. support(mi) == 7) then
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + mom_src_diff*ratio_x_r
                                    mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + mom_src_diff*ratio_y_r
                                    mono_mom_src(3, j, k, l) = mono_mom_src(3, j, k, l) + mom_src_diff*ratio_z_r
                                else
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + mom_src_diff*cos(dir(mi))
                                    mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + mom_src_diff*sin(dir(mi))
                                end if
                            end if
                        end if

                        if (model_eqns /= 4) then
                            if (any(support(mi) == (/5, 7/))) then
                                mono_E_src(j, k, l) = mono_E_src(j, k, l) + mass_src_diff*c**2d0/(small_gamma - 1d0)
                            else
                                mono_E_src(j, k, l) = mono_E_src(j, k, l) + mom_src_diff*c/(small_gamma - 1d0)
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
        !! @param sim_time Simulation time
        !! @param c Sound speed
        !! @param mi Monopole index
    function f_source_temporal(sim_time, c, mi, term_index) result(g)
        !$acc routine seq
        real(kind(0d0)), intent(in) :: sim_time, c
        integer, intent(in) :: mi, term_index
        real(kind(0d0)) :: omega ! angular frequency
        real(kind(0d0)) :: g

        real(kind(0d0)) :: foc_length_factor ! Scale amplitude with radius for spherical support
        ! i.e. Spherical support -> 1/r scaling; Cylindrical support -> 1/sqrt(r) [^-0.5 -> ^-0.85] (empirical correction)

        if (n == 0) then
            foc_length_factor = 1d0
        elseif (p == 0 .and. (.not. cyl_coord)) then ! 2D axisymmetric case is physically 3D
            foc_length_factor = foc_length(mi)**(-0.85d0); 
        else
            foc_length_factor = 1/foc_length(mi); 
        end if

        g = 0d0

        if (pulse(mi) == 1) then
            ! Sine wave
            if (f_is_default(frequency(mi)) .and. f_is_default(wavelength(mi))) then
                wavelength(mi) = length(mi) ! For CI test - TODO remove, add frequency to test case files, and regenerate tests
            end if
            if (f_is_default(frequency(mi))) then
                frequency(mi) = c/wavelength(mi) ! TODO CHANGE
            end if
            if ((sim_time - delay(mi))*frequency(mi) > npulse(mi)) return
            omega = 2d0*pi*frequency(mi)
            if (term_index == 1) then
                g = mag(mi)*sin((sim_time - delay(mi))*omega)/c &
                    + foc_length_factor*mag(mi)*(cos((sim_time - delay(mi))*omega) - 1d0)/omega
            else
                g = mag(mi)*sin((sim_time - delay(mi))*omega)
            end if
        elseif (pulse(mi) == 2) then
            ! Gaussian pulse
            if (f_is_default(gauss_sigma_time(mi))) then
                gauss_sigma_time(mi) = c/gauss_sigma_dist(mi) ! TODO CHANGE
            end if
            if (term_index == 1) then
                g = mag(mi)*dexp(-0.5d0*((sim_time - delay(mi))**2d0)/(gauss_sigma_time(mi)**2d0))/c - &
                    foc_length_factor*mag(mi)*dsqrt(pi/2)*gauss_sigma_time(mi)* &
                    (erf((sim_time - delay(mi))/(dsqrt(2d0)*gauss_sigma_time(mi))) + 1)
            else
                g = mag(mi)*dexp(-0.5d0*((sim_time - delay(mi))**2d0)/(gauss_sigma_time(mi)**2d0))
            end if
        elseif (pulse(mi) == 3) then
            ! Square wave
            if (f_is_default(frequency(mi))) then
                frequency(mi) = c/wavelength(mi) ! TODO CHANGE
            end if
            if ((sim_time - delay(mi))*frequency(mi) > npulse(mi)) return
            omega = 2d0*pi*frequency(mi)
            g = mag(mi)*sign(1d0, sin((sim_time - delay(mi))*omega))
        end if
    end function f_source_temporal

    !> This function give the spatial support of the acoustic source
        !! @param j First coordinate-direction location index
        !! @param k Second coordinate-direction location index
        !! @param l Third coordinate-direction location index
        !! @param mono_loc Nominal source term location
        !! @param mi Monopole index
        !! @param angle Angle of the source term with respect to the x-axis (2D or 2D axisymmetric)
        !! @param ratio_x_r Ratio of the x-component of the source term to the magnitude (3D)
        !! @param ratio_y_r Ratio of the y-component of the source term to the magnitude (3D)
        !! @param ratio_z_r Ratio of the z-component of the source term to the magnitude (3D)
    function f_source_spatial(j, k, l, mono_loc, mi, angle, ratio_x_r, ratio_y_r, ratio_z_r)

        !$acc routine seq
        real(kind(0d0)), dimension(3), intent(in) :: mono_loc
        integer, intent(in) :: mi
        integer, intent(in) :: j, k, l
        real(kind(0d0)), intent(out) :: angle, ratio_x_r, ratio_y_r, ratio_z_r

        integer :: q
        real(kind(0d0)) :: h, hx, hy, hz
        real(kind(0d0)) :: hx_cyl, hy_cyl, hz_cyl
        real(kind(0d0)) :: hxnew, hynew
        real(kind(0d0)) :: hxnew_cyl, hynew_cyl
        real(kind(0d0)) :: sig
        real(kind(0d0)) :: f_source_spatial

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
        elseif (p == 0) then
            sig = maxval((/dx(j), dy(k)/))
        else
            sig = maxval((/dx(j), dy(k), dz(l)/))
        end if
        sig = sig*acoustic_spatial_support_width

        if (n == 0) then ! 1D
            if (support(mi) == 1) then
                ! 1D delta function
                hx = abs(mono_loc(1) - x_cc(j))

                f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                                   dexp(-0.5d0*(hx/(sig/2d0))**2d0)
            elseif (support(mi) == 0) then
                ! Support for all x
                f_source_spatial = 1d0
            end if
        elseif (p == 0) then ! 2D
            hx = mono_loc(1) - x_cc(j)
            hy = mono_loc(2) - y_cc(k)
            if (support(mi) == 1) then
                ! 2D delta function
                h = dsqrt(hx**2d0 + hy**2d0)

                f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                                   dexp(-0.5d0*((h/(sig/2d0))**2d0))
            elseif (support(mi) == 2) then
                !only support for y \pm some value
                if (abs(hy) < length(mi)) then
                    f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                                       dexp(-0.5d0*(hx/(sig/2d0))**2d0)
                else
                    f_source_spatial = 0d0
                end if
            elseif (support(mi) == 3) then
                ! Only support along some line
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)

                ! Rotate actual point by -theta
                hxnew = cos(dir(mi))*hx + sin(dir(mi))*hy
                hynew = -1d0*sin(dir(mi))*hx + cos(dir(mi))*hy
                if (abs(hynew) < mono_loc(3)/2d0) then
                    f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                                       dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)
                else
                    f_source_spatial = 0d0
                end if
            elseif (support(mi) == 4) then
                ! Support for all y
                f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig)* &
                                   dexp(-0.5d0*(hx/sig)**2d0)
            elseif (support(mi) == 5) then
                ! Support along transducer in 2D or 2D axisymmetric coordinate
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)

                current_angle = -atan(hy/(foc_length(mi) - hx))
                angle_half_aperture = asin((aperture(mi)/2d0)/(foc_length(mi)))

                f_source_spatial = 0d0
                if (abs(current_angle) < angle_half_aperture .and. hx < foc_length(mi)) then
                    hxnew = foc_length(mi) - dsqrt(hy**2d0 + (foc_length(mi) - hx)**2d0)
                    f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                                       dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)
                    angle = -atan(hy/(foc_length(mi) - hx))
                end if
            elseif (support(mi) == 7) then
                ! Support along transducer array in 2D
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)

                current_angle = -atan(hy/(foc_length(mi) - hx))
                angle_half_aperture = asin((aperture(mi)/2d0)/(foc_length(mi)))
                angle_per_elem = (2d0*angle_half_aperture - (num_elements(mi) - 1d0)*element_spacing_angle(mi))/num_elements(mi)
                hxnew = foc_length(mi) - dsqrt(hy**2d0 + (foc_length(mi) - hx)**2d0)

                if (element_on(mi) == 0) then ! Full transducer
                    elem_min = 1
                    elem_max = num_elements(mi)
                else ! Transducer element specified
                    elem_min = element_on(mi)
                    elem_max = element_on(mi)
                end if

                f_source_spatial = 0d0 ! If not affected by any element
                do elem = elem_min, elem_max
                    angle_max = angle_half_aperture - (element_spacing_angle(mi) + angle_per_elem)*(elem - 1d0)
                    angle_min = angle_max - angle_per_elem

                    if (current_angle > angle_min .and. current_angle < angle_max .and. hx < foc_length(mi)) then
                        f_source_spatial = dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)/(dsqrt(2d0*pi)*sig/2d0)
                        angle = current_angle
                        exit ! Assume elements don't overlap
                    end if
                end do
            end if

        else ! 3D
            hx = x_cc(j) - mono_loc(1)
            hy = y_cc(k) - mono_loc(2)
            hz = z_cc(l) - mono_loc(3)
            if (support(mi) == 3) then

                ! Rotate actual point by -theta
                hxnew = cos(dir(mi))*hx + sin(dir(mi))*hy
                hynew = -1d0*sin(dir(mi))*hx + cos(dir(mi))*hy

                if (abs(hynew) < length(mi)/2. .and. &
                    abs(hz) < length(mi)/2.) then
                    f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                                       dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)
                else
                    f_source_spatial = 0d0
                end if
            elseif (support(mi) == 4) then
                ! Support for all x,y
                f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig)* &
                                   dexp(-0.5d0*(hz/sig)**2d0)
            elseif (support(mi) == 5) then
                ! Support along transducer in 3D
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)
                hz = z_cc(l) - mono_loc(3)

                current_angle = -atan(dsqrt(hy**2 + hz**2)/(foc_length(mi) - hx))
                angle_half_aperture = asin((aperture(mi)/2d0)/(foc_length(mi)))

                f_source_spatial = 0d0
                if (abs(current_angle) < angle_half_aperture .and. hx < foc_length(mi)) then

                    hxnew = foc_length(mi) - dsqrt(hy**2d0 + hz**2d0 + (foc_length(mi) - hx)**2d0)
                    f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                                       dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)

                    norm = dsqrt(hy**2d0 + hz**2d0 + (foc_length(mi) - hx)**2d0)
                    ratio_x_r = -(hx - foc_length(mi))/norm
                    ratio_y_r = -hy/norm
                    ratio_z_r = -hz/norm
                end if
            elseif (support(mi) == -1) then ! Not working - disabled for now
                ! Support for transducer in 3D cylindrical coordinate
                sig = maxval((/dx(j), dy(k)*sin(dz(l)), dz(l)*cos(dz(l))/))*acoustic_spatial_support_width
                hx_cyl = x_cc(j) - mono_loc(1)
                hy_cyl = y_cc(k)*sin(z_cc(l)) - mono_loc(2)
                hz_cyl = y_cc(k)*cos(z_cc(l)) - mono_loc(3)

                ! Rotate actual point by -theta
                hxnew_cyl = cos(dir(mi))*hx_cyl + sin(dir(mi))*hy_cyl
                hynew_cyl = -1d0*sin(dir(mi))*hx_cyl + cos(dir(mi))*hy_cyl

                f_source_spatial = 0d0
                if (abs(hynew_cyl) < length(mi)/2. .and. &
                    abs(hz_cyl) < length(mi)/2.) then
                    f_source_spatial = 1d0/(dsqrt(2d0*pi)*sig/2d0)* &
                                       dexp(-0.5d0*(hxnew_cyl/(sig/2d0))**2d0)
                end if
            elseif (support(mi) == 7) then
                ! Support along transducer ring in 3D

                poly_side_length = aperture(mi)*sin(pi/num_elements(mi))
                aperture_element_3D = poly_side_length*element_polygon_ratio(mi)

                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)
                hz = z_cc(l) - mono_loc(3)

                f = foc_length(mi)
                R = aperture(mi)/2d0

                if (element_on(mi) == 0) then ! Full transducer
                    elem_min = 1
                    elem_max = num_elements(mi)
                else ! Transducer element specified
                    elem_min = element_on(mi)
                    elem_max = element_on(mi)
                end if

                f_source_spatial = 0d0 ! If not affected by any element
                do elem = elem_min, elem_max
                    angle_elem = 2d0*pi*real(elem, kind(0d0))/real(num_elements(mi), kind(0d0)) + rotate_angle(mi)

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
                        f_source_spatial = dexp(-0.5d0*(hxnew/(sig/2d0))**2d0)/(dsqrt(2d0*pi)*sig/2d0)
                        norm = dsqrt(hy**2d0 + hz**2d0 + (foc_length(mi) - hx)**2d0)
                        ratio_x_r = -(hx - foc_length(mi))/norm
                        ratio_y_r = -hy/norm
                        ratio_z_r = -hz/norm
                    end if
                end do

            end if
        end if

    end function f_source_spatial

    subroutine s_compute_speed_of_sound_monopole(q_cons_elements, q_prim_element, sos, n_tait)
        real(kind(0d0)), dimension(sys_size), intent(in) :: q_cons_elements
        real(kind(0d0)), intent(in) :: q_prim_element
        real(kind(0d0)), intent(out) :: sos, n_tait

        real(kind(0d0)), dimension(num_fluids) :: myalpha_rho, myalpha
        real(kind(0d0)) :: myRho, B_tait
        integer :: q

        myRho = 0d0
        n_tait = 0d0
        B_tait = 0d0

        !$acc loop seq
        do q = 1, num_fluids
            myalpha_rho(q) = q_cons_elements(q)
            myalpha(q) = q_cons_elements(advxb + q - 1)
        end do

        if (bubbles) then
            if (mpp_lim .and. (num_fluids > 2)) then
                !$acc loop seq
                do q = 1, num_fluids
                    myRho = myRho + myalpha_rho(q)
                    n_tait = n_tait + myalpha(q)*gammas(q)
                    B_tait = B_tait + myalpha(q)*pi_infs(q)
                end do
            elseif (num_fluids > 2) then
                !$acc loop seq
                do q = 1, num_fluids - 1
                    myRho = myRho + myalpha_rho(q)
                    n_tait = n_tait + myalpha(q)*gammas(q)
                    B_tait = B_tait + myalpha(q)*pi_infs(q)
                end do
            else
                myRho = myalpha_rho(1)
                n_tait = gammas(1)
                B_tait = pi_infs(1)
            end if
        else
            !$acc loop seq
            do q = 1, num_fluids
                myRho = myRho + myalpha_rho(q)
                n_tait = n_tait + myalpha(q)*gammas(q)
                B_tait = B_tait + myalpha(q)*pi_infs(q)
            end do
        end if

        n_tait = 1d0/n_tait + 1d0

        sos = dsqrt(n_tait*(q_prim_element + ((n_tait - 1d0)/n_tait)*B_tait)/myRho)

    end subroutine s_compute_speed_of_sound_monopole

end module m_monopole
