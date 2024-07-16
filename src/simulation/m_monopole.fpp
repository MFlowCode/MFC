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

    !> This subroutine initializes the monopole module
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
            wavelength(i) = mono(i)%wavelength
            frequency(i) = mono(i)%frequency
            gauss_sigma_dist(i) = mono(i)%gauss_sigma_dist
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

    !> This subroutine updates the rhs by computing the mass, mom, energy sources
    !! @param q_cons_vf Conservative variables
    !! @param q_prim_vf Primitive variables
    !! @param t_step Current time step
    !! @param rhs_vf rhs variables
    subroutine s_monopole_calculations(q_cons_vf, q_prim_vf, t_step, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf !<
        !! This variable contains the WENO-reconstructed values of the cell-average
        !! conservative variables, which are located in q_cons_vf, at cell-interior
        !! Gaussian quadrature points (QP).

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf !<
        !! The primitive variables at cell-interior Gaussian quadrature points. These
        !! are calculated from the conservative variables and gradient magnitude (GM)
        !! of the volume fractions, q_cons_qp and gm_alpha_qp, respectively.

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        integer, intent(in) :: t_step

        real(kind(0d0)) :: q_prim_local, q_cons_local(sys_size)
        real(kind(0d0)) :: sim_time, c, small_gamma
        real(kind(0d0)) :: frequency_local, gauss_sigma_time_local
        real(kind(0d0)) :: mass_src_diff, mom_src_diff
        real(kind(0d0)) :: angle ! Angle with x-axis for mom source term vector
        real(kind(0d0)) :: xyz_to_r_ratios(3) ! [xyz]/r for mom source term vector
        real(kind(0d0)) :: source_temporal, source_spatial

        integer :: i, j, k, l, q !< generic loop variables
        integer :: mi !< monopole index

        logical :: freq_conv_flag, gauss_conv_flag

        integer, parameter :: mass_label = 1, mom_label = 2

        sim_time = t_step*dt

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mono_mass_src(j, k, l) = 0d0
                    mono_mom_src(1, j, k, l) = 0d0
                    mono_e_src(j, k, l) = 0d0
                    if (n > 0) mono_mom_src(2, j, k, l) = 0d0
                    if (p > 0) mono_mom_src(3, j, k, l) = 0d0
                end do
            end do
        end do

        ! Monopoles are looped through sequentially because they can have very different computational costs
        !$acc parallel loop collapse(3) gang vector default(present) private(q_cons_local, q_prim_local, freq_conv_flag, gauss_conv_flag, c, small_gamma, angle, xyz_to_r_ratios)
        !$acc loop seq
        do mi = 1, num_mono
            ! Skip if the pulse has not started yet for sine and square waves
            if (sim_time < delay(mi) .and. (pulse(mi) == 1 .or. pulse(mi) == 3)) cycle

            ! Decide if frequency need to be converted from wavelength
            frequency_local = frequency(mi)
            gauss_sigma_time_local = gauss_sigma_time(mi)
            freq_conv_flag = f_is_default(frequency(mi))
            gauss_conv_flag = f_is_default(gauss_sigma_time(mi))

            do l = 0, p
                do k = 0, n
                    do j = 0, m

                        ! Compute speed of sound
                        !$acc loop seq
                        do q = 1, sys_size
                            q_cons_local(q) = q_cons_vf(q)%sf(j, k, l)
                        end do
                        q_prim_local = q_prim_vf(E_idx)%sf(j, k, l)
                        call s_compute_speed_of_sound_monopole(q_cons_local, q_prim_local, c, small_gamma)

                        ! Wavelength to frequency conversion
                        if (freq_conv_flag .and. (pulse(mi) == 1 .or. pulse(mi) == 3)) then
                            frequency_local = c/wavelength(mi)
                        elseif (gauss_conv_flag .and. pulse(mi) == 2) then
                            gauss_sigma_time_local = c/gauss_sigma_dist(mi)
                        end if

                        ! Update momentum source term
                        call s_source_temporal(sim_time, c, mi, mom_label, frequency_local, gauss_sigma_time_local, source_temporal)
                        call s_source_spatial(j, k, l, loc_mono(:, mi), mi, source_spatial, angle, xyz_to_r_ratios)
                        mom_src_diff = source_temporal*source_spatial

                        if (n == 0) then ! 1D
                            mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + mom_src_diff*sign(1d0, dir(mi)) ! Left or right-going wave

                        elseif (p == 0) then ! 2D
                            if (support(mi) < 5) then ! Planar
                                angle = dir(mi)
                            end if
                            mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + mom_src_diff*cos(angle)
                            mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + mom_src_diff*sin(angle)

                        else ! 3D
                            if (support(mi) < 5) then ! Planar
                                mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + mom_src_diff*cos(dir(mi))
                                mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + mom_src_diff*sin(dir(mi))
                            else
                                mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + mom_src_diff*xyz_to_r_ratios(1)
                                mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + mom_src_diff*xyz_to_r_ratios(2)
                                mono_mom_src(3, j, k, l) = mono_mom_src(3, j, k, l) + mom_src_diff*xyz_to_r_ratios(3)
                            end if
                        end if

                        ! Update mass source term
                        if (support(mi) < 5) then ! Planar
                            mass_src_diff = mom_src_diff/c
                        else
                            call s_source_temporal(sim_time, c, mi, mass_label, frequency_local, gauss_sigma_time_local, source_temporal)
                            call s_source_spatial(j, k, l, loc_mono(:, mi), mi, source_spatial, angle, xyz_to_r_ratios)
                            mass_src_diff = source_temporal*source_spatial
                        end if
                        mono_mass_src(j, k, l) = mono_mass_src(j, k, l) + mass_src_diff

                        ! Update energy source term
                        if (model_eqns /= 4) then
                            mono_E_src(j, k, l) = mono_E_src(j, k, l) + mass_src_diff*c**2d0/(small_gamma - 1d0)
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

    !> This subroutine gives the temporally varying amplitude of the pulse
    !! @param sim_time Simulation time
    !! @param c Sound speed
    !! @param mi Monopole index
    !! @param term_index Index of the term to be calculated (1: mass source, 2: momentum source)
    !! @param frequency_local Frequency at the spatial location for sine and square waves
    !! @param gauss_sigma_time_local sigma in time for Gaussian pulse
    !! @param source Source term amplitude
    subroutine s_source_temporal(sim_time, c, mi, term_index, frequency_local, gauss_sigma_time_local, source)
        !$acc routine seq
        integer, intent(in) :: mi, term_index
        real(kind(0d0)), intent(in) :: sim_time, c
        real(kind(0d0)), intent(in) :: frequency_local, gauss_sigma_time_local
        real(kind(0d0)), intent(out) :: source

        real(kind(0d0)) :: omega ! angular frequency
        real(kind(0d0)) :: foc_length_factor ! Scale amplitude with radius for spherical support
        ! i.e. Spherical support -> 1/r scaling; Cylindrical support -> 1/sqrt(r) [empirical correction: ^-0.5 -> ^-0.85]
        integer, parameter :: mass_label = 1

        if (n == 0) then
            foc_length_factor = 1d0
        elseif (p == 0 .and. (.not. cyl_coord)) then ! 2D axisymmetric case is physically 3D
            foc_length_factor = foc_length(mi)**(-0.85d0); ! Empirical correction
        else
            foc_length_factor = 1/foc_length(mi); 
        end if

        source = 0d0

        if (pulse(mi) == 1) then ! Sine wave
            if ((sim_time - delay(mi))*frequency_local > npulse(mi)) return

            omega = 2d0*pi*frequency_local
            source = mag(mi)*sin((sim_time - delay(mi))*omega)

            if (term_index == mass_label) then
                source = source/c + foc_length_factor*mag(mi)*(cos((sim_time - delay(mi))*omega) - 1d0)/omega
            end if

        elseif (pulse(mi) == 2) then ! Gaussian pulse
            source = mag(mi)*dexp(-0.5d0*((sim_time - delay(mi))**2d0)/(gauss_sigma_time_local**2d0))

            if (term_index == mass_label) then
                source = source/c - &
                         foc_length_factor*mag(mi)*dsqrt(pi/2)*gauss_sigma_time_local* &
                         (erf((sim_time - delay(mi))/(dsqrt(2d0)*gauss_sigma_time_local)) + 1)
            end if

        elseif (pulse(mi) == 3) then ! Square wave
            if ((sim_time - delay(mi))*frequency_local > npulse(mi)) return

            omega = 2d0*pi*frequency_local
            source = mag(mi)*sign(1d0, sin((sim_time - delay(mi))*omega))

        end if
    end subroutine s_source_temporal

    !> This subroutine gives the spatial support of the acoustic source
    !! @param j x-index
    !! @param k y-index
    !! @param l z-index
    !! @param mono_loc Nominal source term location
    !! @param mi Monopole index
    !! @param source Source term amplitude
    !! @param angle Angle of the source term with respect to the x-axis (for 2D or 2D axisymmetric)
    !! @param xyz_to_r_ratios Ratio of the [xyz]-component of the source term to the magnitude (for 3D)
    subroutine s_source_spatial(j, k, l, mono_loc, mi, source, angle, xyz_to_r_ratios)
        !$acc routine seq
        integer, intent(in) :: j, k, l, mi
        real(kind(0d0)), dimension(3), intent(in) :: mono_loc
        real(kind(0d0)), intent(out) :: source, angle, xyz_to_r_ratios(3)

        real(kind(0d0)) :: sig, r(3)

        ! Calculate sig spatial support width
        if (n == 0) then
            sig = dx(j)
        elseif (p == 0) then
            sig = maxval((/dx(j), dy(k)/))
        else
            sig = maxval((/dx(j), dy(k), dz(l)/))
        end if
        sig = sig*acoustic_spatial_support_width

        ! Calculate displacement from monopole location
        r(1) = x_cc(j) - mono_loc(1)
        if (n /= 0) r(2) = y_cc(k) - mono_loc(2)
        if (p /= 0) r(3) = z_cc(l) - mono_loc(3)

        if (any(support(mi) == (/1, 2, 3, 4/))) then
            call s_source_spatial_planar(mi, sig, r, source)
        elseif (any(support(mi) == (/5, 6, 7/))) then
            call s_source_spatial_transducer(mi, sig, r, source, angle, xyz_to_r_ratios)
        elseif (any(support(mi) == (/9, 10, 11/))) then
            call s_source_spatial_transducer_array(mi, sig, r, source, angle, xyz_to_r_ratios)
        end if
    end subroutine s_source_spatial

    !> This subroutine calculates the spatial support for planar acoustic sources in 1D, 2D, and 3D
    !! @param mi Monopole index
    !! @param sig Sigma value for the Gaussian distribution
    !! @param r Displacement from source to current point
    !! @param source Source term amplitude
    subroutine s_source_spatial_planar(mi, sig, r, source)
        !$acc routine seq
        integer, intent(in) :: mi
        real(kind(0d0)), intent(in) :: sig, r(3)
        real(kind(0d0)), intent(out) :: source

        real(kind(0d0)) :: dist, rxnew, rynew

        source = 0d0

        if (n == 0) then ! 1D - Only support 1 is allowed
            source = 1d0/(dsqrt(2d0*pi)*sig/2d0)*dexp(-0.5d0*(r(1)/(sig/2d0))**2d0)

        elseif (p == 0) then ! 2D
            if (support(mi) == 1) then
                dist = dsqrt(r(1)**2d0 + r(2)**2d0)
                source = 1d0/(dsqrt(2d0*pi)*sig/2d0)*dexp(-0.5d0*((dist/(sig/2d0))**2d0))

            elseif (support(mi) == 2) then ! Only support for y \pm some value
                if (abs(r(2)) < length(mi)) then
                    source = 1d0/(dsqrt(2d0*pi)*sig/2d0)*dexp(-0.5d0*(r(1)/(sig/2d0))**2d0)
                end if

                ! elseif (support(mi) == 3) then
                !     rxnew = cos(dir(mi))*r(1) + sin(dir(mi))*r(2)
                !     rynew = -1d0*sin(dir(mi))*r(1) + cos(dir(mi))*r(2)
                !     if (abs(rynew) < mono_loc(3)/2d0) then
                !         source = 1d0/(dsqrt(2d0*pi)*sig/2d0)*dexp(-0.5d0*(rxnew/(sig/2d0))**2d0)
                !     end if

            elseif (support(mi) == 4) then ! Support for all y
                source = 1d0/(dsqrt(2d0*pi)*sig)*dexp(-0.5d0*(r(1)/sig)**2d0)
            end if

        else ! 3D
            if (support(mi) == 4) then ! Support for all x and y
                source = 1d0/(dsqrt(2d0*pi)*sig)*dexp(-0.5d0*(r(3)/sig)**2d0)
            end if

        end if
    end subroutine s_source_spatial_planar

    !> This subroutine calculates the spatial support for a single transducer in 2D, 2D axisymmetric, and 3D
    !! @param mi Monopole index
    !! @param sig Sigma value for the Gaussian distribution
    !! @param r Displacement from source to current point
    !! @param source Source term amplitude
    !! @param angle Angle of the source term with respect to the x-axis (for 2D or 2D axisymmetric)
    !! @param xyz_to_r_ratios Ratio of the [xyz]-component of the source term to the magnitude (for 3D)
    subroutine s_source_spatial_transducer(mi, sig, r, source, angle, xyz_to_r_ratios)
        !$acc routine seq
        integer, intent(in) :: mi
        real(kind(0d0)), intent(in) :: sig, r(3)
        real(kind(0d0)), intent(out) :: source, angle, xyz_to_r_ratios(3)

        real(kind(0d0)) :: current_angle, angle_half_aperture, dist, norm

        source = 0d0

        if (support(mi) == 5 .or. support(mi) == 6) then ! 2D or 2D axisymmetric
            current_angle = -atan(r(2)/(foc_length(mi) - r(1)))
            angle_half_aperture = asin((aperture(mi)/2d0)/(foc_length(mi)))

            if (abs(current_angle) < angle_half_aperture .and. r(1) < foc_length(mi)) then
                dist = foc_length(mi) - dsqrt(r(2)**2d0 + (foc_length(mi) - r(1))**2d0)
                source = 1d0/(dsqrt(2d0*pi)*sig/2d0)*dexp(-0.5d0*(dist/(sig/2d0))**2d0)
                angle = -atan(r(2)/(foc_length(mi) - r(1)))
            end if

        elseif (support(mi) == 7) then ! 3D
            current_angle = -atan(dsqrt(r(2)**2 + r(3)**2)/(foc_length(mi) - r(1)))
            angle_half_aperture = asin((aperture(mi)/2d0)/(foc_length(mi)))

            if (abs(current_angle) < angle_half_aperture .and. r(1) < foc_length(mi)) then
                dist = foc_length(mi) - dsqrt(r(2)**2d0 + r(3)**2d0 + (foc_length(mi) - r(1))**2d0)
                source = 1d0/(dsqrt(2d0*pi)*sig/2d0)*dexp(-0.5d0*(dist/(sig/2d0))**2d0)

                norm = dsqrt(r(2)**2d0 + r(3)**2d0 + (foc_length(mi) - r(1))**2d0)
                xyz_to_r_ratios(1) = -(r(1) - foc_length(mi))/norm
                xyz_to_r_ratios(2) = -r(2)/norm
                xyz_to_r_ratios(3) = -r(3)/norm
            end if

        end if
    end subroutine s_source_spatial_transducer

    !> This subroutine calculates the spatial support for multiple transducers in 2D, 2D axisymmetric, and 3D
    !! @param mi Monopole index
    !! @param sig Sigma value for the Gaussian distribution
    !! @param r Displacement from source to current point
    !! @param source Source term amplitude
    !! @param angle Angle of the source term with respect to the x-axis (for 2D or 2D axisymmetric)
    !! @param xyz_to_r_ratios Ratio of the [xyz]-component of the source term to the magnitude (for 3D)
    subroutine s_source_spatial_transducer_array(mi, sig, r, source, angle, xyz_to_r_ratios)
        !$acc routine seq
        integer, intent(in) :: mi
        real(kind(0d0)), intent(in) :: sig, r(3)
        real(kind(0d0)), intent(out) :: source, angle, xyz_to_r_ratios(3)

        integer :: elem, elem_min, elem_max
        real(kind(0d0)) :: current_angle, angle_half_aperture, angle_per_elem, dist
        real(kind(0d0)) :: angle_min, angle_max, norm
        real(kind(0d0)) :: poly_side_length, aperture_element_3D, angle_elem
        real(kind(0d0)) :: x2, y2, z2, x3, y3, z3, C, f, half_apert, dist_interp_to_elem_center

        if (element_on(mi) == 0) then ! Full transducer
            elem_min = 1
            elem_max = num_elements(mi)
        else ! Transducer element specified
            elem_min = element_on(mi)
            elem_max = element_on(mi)
        end if

        source = 0d0 ! If not affected by any element

        if (support(mi) == 9 .or. support(mi) == 10) then ! 2D or 2D axisymmetric
            current_angle = -atan(r(2)/(foc_length(mi) - r(1)))
            angle_half_aperture = asin((aperture(mi)/2d0)/(foc_length(mi)))
            angle_per_elem = (2d0*angle_half_aperture - (num_elements(mi) - 1d0)*element_spacing_angle(mi))/num_elements(mi)
            dist = foc_length(mi) - dsqrt(r(2)**2d0 + (foc_length(mi) - r(1))**2d0)

            do elem = elem_min, elem_max
                angle_max = angle_half_aperture - (element_spacing_angle(mi) + angle_per_elem)*(elem - 1d0)
                angle_min = angle_max - angle_per_elem

                if (current_angle > angle_min .and. current_angle < angle_max .and. r(1) < foc_length(mi)) then
                    source = dexp(-0.5d0*(dist/(sig/2d0))**2d0)/(dsqrt(2d0*pi)*sig/2d0)
                    angle = current_angle
                    exit ! Assume elements don't overlap
                end if
            end do

        elseif (support(mi) == 11) then ! 3D
            poly_side_length = aperture(mi)*sin(pi/num_elements(mi))
            aperture_element_3D = poly_side_length*element_polygon_ratio(mi)
            f = foc_length(mi)
            half_apert = aperture(mi)/2d0

            do elem = elem_min, elem_max
                angle_elem = 2d0*pi*real(elem, kind(0d0))/real(num_elements(mi), kind(0d0)) + rotate_angle(mi)

                ! Point 2 is the elem center
                x2 = f - dsqrt(f**2 - half_apert**2)
                y2 = half_apert*cos(angle_elem)
                z2 = half_apert*sin(angle_elem)

                ! Construct a plane normal to the line from the focal point to the elem center,
                ! Point 3 is the intercept of the plane and the line from the focal point to the current location
                C = f**2d0/((r(1) - f)*(x2 - f) + r(2)*y2 + r(3)*z2) ! Constant for intermediate step
                x3 = C*(r(1) - f) + f
                y3 = C*r(2)
                z3 = C*r(3)

                dist_interp_to_elem_center = dsqrt((x2 - x3)**2d0 + (y2 - y3)**2d0 + (z2 - z3)**2d0)
                if ((dist_interp_to_elem_center < aperture_element_3D/2d0) .and. (r(1) < f)) then
                    dist = dsqrt((x3 - r(1))**2d0 + (y3 - r(2))**2d0 + (z3 - r(3))**2d0)
                    source = dexp(-0.5d0*(dist/(sig/2d0))**2d0)/(dsqrt(2d0*pi)*sig/2d0)

                    norm = dsqrt(r(2)**2d0 + r(3)**2d0 + (f - r(1))**2d0)
                    xyz_to_r_ratios(1) = -(r(1) - f)/norm
                    xyz_to_r_ratios(2) = -r(2)/norm
                    xyz_to_r_ratios(3) = -r(3)/norm
                end if

            end do

        end if
    end subroutine s_source_spatial_transducer_array

    !> This subroutine calculates the speed of sound within the monopole module
    !! @param q_cons_local Conserved quantities of the elements
    !! @param q_prim_local Primitive quantity of the element
    !! @param c Speed of sound
    !! @param n_tait Gas constant (small gamma)
    subroutine s_compute_speed_of_sound_monopole(q_cons_local, q_prim_local, c, n_tait)
        real(kind(0d0)), dimension(sys_size), intent(in) :: q_cons_local
        real(kind(0d0)), intent(in) :: q_prim_local
        real(kind(0d0)), intent(out) :: c, n_tait

        real(kind(0d0)), dimension(num_fluids) :: myalpha_rho, myalpha
        real(kind(0d0)) :: myRho, B_tait
        integer :: q

        myRho = 0d0
        n_tait = 0d0
        B_tait = 0d0

        !$acc loop seq
        do q = 1, num_fluids
            myalpha_rho(q) = q_cons_local(q)
            myalpha(q) = q_cons_local(advxb + q - 1)
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

        c = dsqrt(n_tait*(q_prim_local + ((n_tait - 1d0)/n_tait)*B_tait)/myRho)

    end subroutine s_compute_speed_of_sound_monopole

end module m_monopole
