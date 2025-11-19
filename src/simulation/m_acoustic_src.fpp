!>
!! @file m_acoustic_src.fpp
!! @brief Contains module m_acoustic_src

#:include 'macros.fpp'

!> @brief The module contains the subroutines used to create a acoustic source pressure source term
module m_acoustic_src

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_bubbles              !< Bubble dynamic routines

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_constants            !< Definitions of the constants

    implicit none
    private; public :: s_initialize_acoustic_src, s_precalculate_acoustic_spatial_sources, s_acoustic_src_calculations

    integer, allocatable, dimension(:) :: pulse, support
    $:GPU_DECLARE(create='[pulse,support]')

    logical, allocatable, dimension(:) :: dipole
    $:GPU_DECLARE(create='[dipole]')

    real(wp), allocatable, target, dimension(:, :) :: loc_acoustic
    $:GPU_DECLARE(create='[loc_acoustic]')

    real(wp), allocatable, dimension(:) :: mag, length, height, wavelength, frequency
    real(wp), allocatable, dimension(:) :: gauss_sigma_dist, gauss_sigma_time, npulse, dir, delay
    $:GPU_DECLARE(create='[mag,length,height,wavelength,frequency]')
    $:GPU_DECLARE(create='[gauss_sigma_dist,gauss_sigma_time,npulse,dir,delay]')

    real(wp), allocatable, dimension(:) :: foc_length, aperture
    $:GPU_DECLARE(create='[foc_length,aperture]')

    real(wp), allocatable, dimension(:) :: element_spacing_angle, element_polygon_ratio, rotate_angle
    $:GPU_DECLARE(create='[element_spacing_angle,element_polygon_ratio,rotate_angle]')

    real(wp), allocatable, dimension(:) :: bb_bandwidth, bb_lowest_freq
    $:GPU_DECLARE(create='[bb_bandwidth,bb_lowest_freq]')

    integer, allocatable, dimension(:) :: num_elements, element_on, bb_num_freq
    $:GPU_DECLARE(create='[num_elements,element_on,bb_num_freq]')

    !> @name Acoustic source terms
    !> @{
    real(wp), allocatable, dimension(:, :, :) :: mass_src, e_src
    real(wp), allocatable, dimension(:, :, :, :) :: mom_src
    !> @}
    $:GPU_DECLARE(create='[mass_src,e_src,mom_src]')

    integer, dimension(:), allocatable :: source_spatials_num_points !< Number of non-zero source grid points for each source
    $:GPU_DECLARE(create='[source_spatials_num_points]')

    type(source_spatial_type), dimension(:), allocatable :: source_spatials !< Data of non-zero source grid points for each source
    $:GPU_DECLARE(create='[source_spatials]')

contains

    !> This subroutine initializes the acoustic source module
    impure subroutine s_initialize_acoustic_src
        integer :: i, j !< generic loop variables

        @:ALLOCATE(loc_acoustic(1:3, 1:num_source), mag(1:num_source), dipole(1:num_source), support(1:num_source), length(1:num_source), height(1:num_source), wavelength(1:num_source), frequency(1:num_source), gauss_sigma_dist(1:num_source), gauss_sigma_time(1:num_source), foc_length(1:num_source), aperture(1:num_source), npulse(1:num_source), pulse(1:num_source), dir(1:num_source), delay(1:num_source), element_polygon_ratio(1:num_source), rotate_angle(1:num_source), element_spacing_angle(1:num_source), num_elements(1:num_source), element_on(1:num_source), bb_num_freq(1:num_source), bb_bandwidth(1:num_source), bb_lowest_freq(1:num_source))

        do i = 1, num_source
            do j = 1, 3
                loc_acoustic(j, i) = acoustic(i)%loc(j)
            end do
            mag(i) = acoustic(i)%mag
            dipole(i) = acoustic(i)%dipole
            support(i) = acoustic(i)%support
            length(i) = acoustic(i)%length
            height(i) = acoustic(i)%height
            wavelength(i) = acoustic(i)%wavelength
            frequency(i) = acoustic(i)%frequency
            gauss_sigma_dist(i) = acoustic(i)%gauss_sigma_dist
            gauss_sigma_time(i) = acoustic(i)%gauss_sigma_time
            foc_length(i) = acoustic(i)%foc_length
            aperture(i) = acoustic(i)%aperture
            npulse(i) = acoustic(i)%npulse
            pulse(i) = acoustic(i)%pulse
            dir(i) = acoustic(i)%dir
            element_spacing_angle(i) = acoustic(i)%element_spacing_angle
            element_polygon_ratio(i) = acoustic(i)%element_polygon_ratio
            num_elements(i) = acoustic(i)%num_elements
            bb_num_freq(i) = acoustic(i)%bb_num_freq
            bb_bandwidth(i) = acoustic(i)%bb_bandwidth
            bb_lowest_freq(i) = acoustic(i)%bb_lowest_freq

            if (acoustic(i)%element_on == dflt_int) then
                element_on(i) = 0
            else
                element_on(i) = acoustic(i)%element_on
            end if
            if (f_is_default(acoustic(i)%rotate_angle)) then
                rotate_angle(i) = 0._wp
            else
                rotate_angle(i) = acoustic(i)%rotate_angle
            end if
            if (f_is_default(acoustic(i)%delay)) then ! m_checker guarantees acoustic(i)%delay is set for pulse = 2 (Gaussian)
                delay(i) = 0._wp ! Defaults to zero for sine and square waves
            else
                delay(i) = acoustic(i)%delay
            end if
        end do
        $:GPU_UPDATE(device='[loc_acoustic,mag,dipole,support,length, &
            & height,wavelength,frequency,gauss_sigma_dist, &
            & gauss_sigma_time,foc_length,aperture,npulse,pulse, &
            & dir,delay,element_polygon_ratio,rotate_angle, &
            & element_spacing_angle,num_elements,element_on, &
            & bb_num_freq,bb_bandwidth,bb_lowest_freq]')

        @:ALLOCATE(mass_src(0:m, 0:n, 0:p))
        @:ALLOCATE(mom_src(1:num_vels, 0:m, 0:n, 0:p))
        @:ALLOCATE(E_src(0:m, 0:n, 0:p))

    end subroutine s_initialize_acoustic_src

    !> This subroutine updates the rhs by computing the mass, mom, energy sources
    !! @param q_cons_vf Conservative variables
    !! @param q_prim_vf Primitive variables
    !! @param t_step Current time step
    !! @param rhs_vf rhs variables
    impure subroutine s_acoustic_src_calculations(q_cons_vf, q_prim_vf, t_step, rhs_vf)

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

        real(wp), dimension(num_fluids) :: myalpha, myalpha_rho
        real(wp) :: myRho, B_tait
        real(wp) :: sim_time, c, small_gamma
        real(wp) :: frequency_local, gauss_sigma_time_local
        real(wp) :: mass_src_diff, mom_src_diff
        real(wp) :: source_temporal
        real(wp) :: period_BB !< period of each sine wave in broadband source
        real(wp) :: sl_BB !< spectral level at each frequency
        real(wp) :: ffre_BB !< source term corresponding to each frequency
        real(wp) :: sum_BB !< total source term for the broadband wave
        real(wp), allocatable, dimension(:) :: phi_rn !< random phase shift for each frequency

        integer :: i, j, k, l, q !< generic loop variables
        integer :: ai !< acoustic source index
        integer :: num_points

        logical :: freq_conv_flag, gauss_conv_flag

        integer, parameter :: mass_label = 1, mom_label = 2

        sim_time = t_step*dt

        $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mass_src(j, k, l) = 0._wp
                    mom_src(1, j, k, l) = 0._wp
                    e_src(j, k, l) = 0._wp
                    if (n > 0) mom_src(2, j, k, l) = 0._wp
                    if (p > 0) mom_src(3, j, k, l) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Keep outer loop sequel because different sources can have very different number of points
        do ai = 1, num_source
            ! Skip if the pulse has not started yet for sine and square waves
            if (.not. (sim_time < delay(ai) .and. (pulse(ai) == 1 .or. pulse(ai) == 3))) then

                ! Decide if frequency need to be converted from wavelength
                freq_conv_flag = f_is_default(frequency(ai))
                gauss_conv_flag = f_is_default(gauss_sigma_time(ai))

                num_points = source_spatials_num_points(ai) ! Use scalar to force firstprivate to prevent GPU bug

                ! Calculate the broadband source
                period_BB = 0._wp
                sl_BB = 0._wp
                ffre_BB = 0._wp
                sum_BB = 0._wp

                ! Allocate buffers for random phase shift
                allocate (phi_rn(1:bb_num_freq(ai)))
                phi_rn(1:bb_num_freq(ai)) = 0._wp

                if (pulse(ai) == 4) then
                    call random_number(phi_rn(1:bb_num_freq(ai)))
                    ! Ensure all the ranks have the same random phase shift
                    call s_mpi_send_random_number(phi_rn, bb_num_freq(ai))
                end if

                do k = 1, bb_num_freq(ai)
                    ! Acoustic period of the wave at each discrete frequency
                    period_BB = 1._wp/(bb_lowest_freq(ai) + k*bb_bandwidth(ai))
                    ! Spectral level at each frequency
                    sl_BB = broadband_spectral_level_constant*mag(ai) + k*mag(ai)/broadband_spectral_level_growth_rate
                    ! Source term corresponding to each frequencies
                    ffre_BB = sqrt((2._wp*sl_BB*bb_bandwidth(ai)))*cos((sim_time)*2._wp*pi/period_BB + 2._wp*pi*phi_rn(k))
                    ! Sum up the source term of each frequency to obtain the total source term for broadband wave
                    sum_BB = sum_BB + ffre_BB
                end do

                deallocate (phi_rn)

                $:GPU_PARALLEL_LOOP(private='[myalpha,myalpha_rho, myRho, B_tait,c,  small_gamma, frequency_local, gauss_sigma_time_local, mass_src_diff, mom_src_diff, source_temporal, j, k, l, q ]', copyin = '[sum_BB, freq_conv_flag, gauss_conv_flag, sim_time]')
                do i = 1, num_points
                    j = source_spatials(ai)%coord(1, i)
                    k = source_spatials(ai)%coord(2, i)
                    l = source_spatials(ai)%coord(3, i)

                    ! Compute speed of sound
                    myRho = 0._wp
                    B_tait = 0._wp
                    small_gamma = 0._wp

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, num_fluids
                        myalpha_rho(q) = q_cons_vf(q)%sf(j, k, l)
                        myalpha(q) = q_cons_vf(advxb + q - 1)%sf(j, k, l)
                    end do

                    if (bubbles_euler) then
                        if (num_fluids > 2) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, num_fluids - 1
                                myRho = myRho + myalpha_rho(q)
                                B_tait = B_tait + myalpha(q)*pi_infs(q)
                                small_gamma = small_gamma + myalpha(q)*gammas(q)
                            end do
                        else
                            myRho = myalpha_rho(1)
                            B_tait = pi_infs(1)
                            small_gamma = gammas(1)
                        end if
                    end if

                    if ((.not. bubbles_euler) .or. (mpp_lim .and. (num_fluids > 2))) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_fluids
                            myRho = myRho + myalpha_rho(q)
                            B_tait = B_tait + myalpha(q)*pi_infs(q)
                            small_gamma = small_gamma + myalpha(q)*gammas(q)
                        end do
                    end if

                    small_gamma = 1._wp/small_gamma + 1._wp
                    c = sqrt(small_gamma*(q_prim_vf(E_idx)%sf(j, k, l) + ((small_gamma - 1._wp)/small_gamma)*B_tait)/myRho)

                    ! Wavelength to frequency conversion
                    if (pulse(ai) == 1 .or. pulse(ai) == 3) frequency_local = f_frequency_local(freq_conv_flag, ai, c)
                    if (pulse(ai) == 2) gauss_sigma_time_local = f_gauss_sigma_time_local(gauss_conv_flag, ai, c)

                    ! Update momentum source term
                    call s_source_temporal(sim_time, c, ai, mom_label, frequency_local, gauss_sigma_time_local, source_temporal, sum_BB)
                    mom_src_diff = source_temporal*source_spatials(ai)%val(i)

                    if (dipole(ai)) then ! Double amplitude & No momentum source term (only works for Planar)
                        mass_src(j, k, l) = mass_src(j, k, l) + 2._wp*mom_src_diff/c
                        if (model_eqns /= 4) E_src(j, k, l) = E_src(j, k, l) + 2._wp*mom_src_diff*c/(small_gamma - 1._wp)
                        cycle
                    end if

                    if (n == 0) then ! 1D
                        mom_src(1, j, k, l) = mom_src(1, j, k, l) + mom_src_diff*sign(1._wp, dir(ai)) ! Left or right-going wave

                    elseif (p == 0) then ! 2D
                        if (support(ai) < 5) then ! Planar
                            mom_src(1, j, k, l) = mom_src(1, j, k, l) + mom_src_diff*cos(dir(ai))
                            mom_src(2, j, k, l) = mom_src(2, j, k, l) + mom_src_diff*sin(dir(ai))
                        else
                            mom_src(1, j, k, l) = mom_src(1, j, k, l) + mom_src_diff*cos(source_spatials(ai)%angle(i))
                            mom_src(2, j, k, l) = mom_src(2, j, k, l) + mom_src_diff*sin(source_spatials(ai)%angle(i))
                        end if

                    else ! 3D
                        if (support(ai) < 5) then ! Planar
                            mom_src(1, j, k, l) = mom_src(1, j, k, l) + mom_src_diff*cos(dir(ai))
                            mom_src(2, j, k, l) = mom_src(2, j, k, l) + mom_src_diff*sin(dir(ai))
                        else
                            mom_src(1, j, k, l) = mom_src(1, j, k, l) + mom_src_diff*source_spatials(ai)%xyz_to_r_ratios(1, i)
                            mom_src(2, j, k, l) = mom_src(2, j, k, l) + mom_src_diff*source_spatials(ai)%xyz_to_r_ratios(2, i)
                            mom_src(3, j, k, l) = mom_src(3, j, k, l) + mom_src_diff*source_spatials(ai)%xyz_to_r_ratios(3, i)
                        end if
                    end if

                    ! Update mass source term
                    if (support(ai) < 5) then ! Planar
                        mass_src_diff = mom_src_diff/c
                    else ! Spherical or cylindrical support
                        ! Mass source term must be calculated differently using a correction term for spherical and cylindrical support
                        call s_source_temporal(sim_time, c, ai, mass_label, frequency_local, gauss_sigma_time_local, source_temporal, sum_BB)
                        mass_src_diff = source_temporal*source_spatials(ai)%val(i)
                    end if
                    mass_src(j, k, l) = mass_src(j, k, l) + mass_src_diff

                    ! Update energy source term
                    if (model_eqns /= 4) then
                        E_src(j, k, l) = E_src(j, k, l) + mass_src_diff*c**2._wp/(small_gamma - 1._wp)
                    end if

                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end do

        ! Update the rhs variables
        $:GPU_PARALLEL_LOOP(private='[j,k,l]',collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = contxb, contxe
                        rhs_vf(q)%sf(j, k, l) = rhs_vf(q)%sf(j, k, l) + mass_src(j, k, l)
                    end do
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = momxb, momxe
                        rhs_vf(q)%sf(j, k, l) = rhs_vf(q)%sf(j, k, l) + mom_src(q - contxe, j, k, l)
                    end do
                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + e_src(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
    end subroutine s_acoustic_src_calculations

    !> This subroutine gives the temporally varying amplitude of the pulse
    !! @param sim_time Simulation time
    !! @param c Sound speed
    !! @param ai Acoustic source index
    !! @param term_index Index of the term to be calculated (1: mass source, 2: momentum source)
    !! @param frequency_local Frequency at the spatial location for sine and square waves
    !! @param gauss_sigma_time_local sigma in time for Gaussian pulse
    !! @param source Source term amplitude
    elemental subroutine s_source_temporal(sim_time, c, ai, term_index, frequency_local, gauss_sigma_time_local, source, sum_BB)
        $:GPU_ROUTINE(parallelism='[seq]')
        integer, intent(in) :: ai, term_index
        real(wp), intent(in) :: sim_time, c, sum_BB
        real(wp), intent(in) :: frequency_local, gauss_sigma_time_local
        real(wp), intent(out) :: source

        real(wp) :: omega ! angular frequency
        real(wp) :: sine_wave ! sine function for square wave
        real(wp) :: foc_length_factor ! Scale amplitude with radius for spherical support
        ! i.e. Spherical support -> 1/r scaling; Cylindrical support -> 1/sqrt(r) [empirical correction: ^-0.5 -> ^-0.85]
        integer, parameter :: mass_label = 1

        if (n == 0) then
            foc_length_factor = 1._wp
        elseif (p == 0 .and. (.not. cyl_coord)) then ! 2D axisymmetric case is physically 3D
            foc_length_factor = foc_length(ai)**(-0.85_wp); ! Empirical correction
        else
            foc_length_factor = 1/foc_length(ai); 
        end if

        source = 0._wp

        if (pulse(ai) == 1) then ! Sine wave
            if ((sim_time - delay(ai))*frequency_local > npulse(ai)) return

            omega = 2._wp*pi*frequency_local
            source = mag(ai)*sin((sim_time - delay(ai))*omega)

            if (term_index == mass_label) then
                source = source/c + foc_length_factor*mag(ai)*(cos((sim_time - delay(ai))*omega) - 1._wp)/omega
            end if

        elseif (pulse(ai) == 2) then ! Gaussian pulse
            source = mag(ai)*exp(-0.5_wp*((sim_time - delay(ai))**2._wp)/(gauss_sigma_time_local**2._wp))

            if (term_index == mass_label) then
                source = source/c - &
                         foc_length_factor*mag(ai)*sqrt(pi/2)*gauss_sigma_time_local* &
                         (erf((sim_time - delay(ai))/(sqrt(2._wp)*gauss_sigma_time_local)) + 1)
            end if

        elseif (pulse(ai) == 3) then ! Square wave
            if ((sim_time - delay(ai))*frequency_local > npulse(ai)) return

            omega = 2._wp*pi*frequency_local
            sine_wave = sin((sim_time - delay(ai))*omega)
            source = mag(ai)*sign(1._wp, sine_wave)

            ! Prevent max-norm differences due to compilers to pass CI
            if (abs(sine_wave) < 1.e-2_wp) then
                source = mag(ai)*sine_wave*1.e2_wp
            end if

        elseif (pulse(ai) == 4) then ! Broadband wave
            source = sum_BB
        end if
    end subroutine s_source_temporal

    !> This subroutine identifies and precalculates the non-zero acoustic spatial sources before time-stepping
    impure subroutine s_precalculate_acoustic_spatial_sources
        integer :: j, k, l, ai
        integer :: count
        integer :: dim
        real(wp) :: source_spatial, angle, xyz_to_r_ratios(3)
        real(wp), parameter :: threshold = 1.e-10_wp

        if (n == 0) then
            dim = 1
        elseif (p == 0) then
            dim = 2
        else
            dim = 3
        end if

        @:ALLOCATE(source_spatials_num_points(1:num_source))
        @:ALLOCATE(source_spatials(1:num_source))

        do ai = 1, num_source
            ! First pass: Count the number of points for each source
            count = 0
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        call s_source_spatial(j, k, l, loc_acoustic(:, ai), ai, source_spatial, angle, xyz_to_r_ratios)
                        if (abs(source_spatial) < threshold) cycle
                        count = count + 1
                    end do
                end do
            end do
            source_spatials_num_points(ai) = count

            ! Allocate arrays with the correct size

            @:ALLOCATE(source_spatials(ai)%coord(1:3, 1:count))
            @:ALLOCATE(source_spatials(ai)%val(1:count))
            @:ALLOCATE(source_spatials(ai)%angle(1:count))
            @:ALLOCATE(source_spatials(ai)%xyz_to_r_ratios(1:3, 1:count))

            @:ACC_SETUP_source_spatials(source_spatials(ai))

            ! Second pass: Store the values
            count = 0 ! Reset counter
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        call s_source_spatial(j, k, l, loc_acoustic(:, ai), ai, source_spatial, angle, xyz_to_r_ratios)
                        if (abs(source_spatial) < threshold) cycle
                        count = count + 1
                        source_spatials(ai)%coord(1, count) = j
                        source_spatials(ai)%coord(2, count) = k
                        source_spatials(ai)%coord(3, count) = l
                        source_spatials(ai)%val(count) = source_spatial
                        if (support(ai) >= 5) then
                            if (dim == 2) source_spatials(ai)%angle(count) = angle
                            if (dim == 3) source_spatials(ai)%xyz_to_r_ratios(1:3, count) = xyz_to_r_ratios
                        end if
                    end do
                end do
            end do

            if (source_spatials_num_points(ai) /= count) then
                call s_mpi_abort('Fatal Error: Inconsistent allocation of source_spatials')
            end if

            $:GPU_UPDATE(device='[source_spatials(ai)%coord]')
            $:GPU_UPDATE(device='[source_spatials(ai)%val]')
            if (support(ai) >= 5) then
                if (dim == 2) then
                    $:GPU_UPDATE(device='[source_spatials(ai)%angle]')
                end if
                if (dim == 3) then
                    $:GPU_UPDATE(device='[source_spatials(ai)%xyz_to_r_ratios]')
                end if
            end if

        end do

#ifdef MFC_DEBUG
        do ai = 1, num_source
            write (*, '(A,I2,A,I8,A)') 'Acoustic source ', ai, ' has ', source_spatials_num_points(ai), &
                ' grid points with non-zero source term'
        end do
#endif

    end subroutine s_precalculate_acoustic_spatial_sources

    !> This subroutine gives the spatial support of the acoustic source
    !! @param j x-index
    !! @param k y-index
    !! @param l z-index
    !! @param loc Nominal source term location
    !! @param ai Acoustic source index
    !! @param source Source term amplitude
    !! @param angle Angle of the source term with respect to the x-axis (for 2D or 2D axisymmetric)
    !! @param xyz_to_r_ratios Ratios of the [xyz]-component of the source term to the magnitude (for 3D)
    subroutine s_source_spatial(j, k, l, loc, ai, source, angle, xyz_to_r_ratios)
        integer, intent(in) :: j, k, l, ai
        real(wp), dimension(3), intent(in) :: loc
        real(wp), intent(out) :: source, angle, xyz_to_r_ratios(3)

        real(wp) :: sig, r(3)

        ! Calculate sig spatial support width
        if (n == 0) then
            sig = dx(j)
        elseif (p == 0) then
            sig = maxval((/dx(j), dy(k)/))
        else
            sig = maxval((/dx(j), dy(k), dz(l)/))
        end if
        sig = sig*acoustic_spatial_support_width

        ! Calculate displacement from acoustic source location
        r(1) = x_cc(j) - loc(1)
        if (n /= 0) r(2) = y_cc(k) - loc(2)
        if (p /= 0) r(3) = z_cc(l) - loc(3)

        if (any(support(ai) == (/1, 2, 3, 4/))) then
            call s_source_spatial_planar(ai, sig, r, source)
        elseif (any(support(ai) == (/5, 6, 7/))) then
            call s_source_spatial_transducer(ai, sig, r, source, angle, xyz_to_r_ratios)
        elseif (any(support(ai) == (/9, 10, 11/))) then
            call s_source_spatial_transducer_array(ai, sig, r, source, angle, xyz_to_r_ratios)
        end if
    end subroutine s_source_spatial

    !> This subroutine calculates the spatial support for planar acoustic sources in 1D, 2D, and 3D
    !! @param ai Acoustic source index
    !! @param sig Sigma value for the Gaussian distribution
    !! @param r Displacement from source to current point
    !! @param source Source term amplitude
    subroutine s_source_spatial_planar(ai, sig, r, source)
        integer, intent(in) :: ai
        real(wp), intent(in) :: sig, r(3)
        real(wp), intent(out) :: source

        real(wp) :: dist

        source = 0._wp

        if (support(ai) == 1) then ! 1D
            source = 1._wp/(sqrt(2._wp*pi)*sig/2._wp)*exp(-0.5_wp*(r(1)/(sig/2._wp))**2._wp)

        elseif (support(ai) == 2 .or. support(ai) == 3) then ! 2D or 3D
            ! If we let unit vector e = (cos(dir), sin(dir)),
            dist = r(1)*cos(dir(ai)) + r(2)*sin(dir(ai)) ! dot(r,e)
            if ((r(1) - dist*cos(dir(ai)))**2._wp + (r(2) - dist*sin(dir(ai)))**2._wp < 0.25_wp*length(ai)**2._wp) then ! |r - dist*e| < length/2
                if (support(ai) /= 3 .or. abs(r(3)) < 0.25_wp*height(ai)) then ! additional height constraint for 3D
                    source = 1._wp/(sqrt(2._wp*pi)*sig/2._wp)*exp(-0.5_wp*(dist/(sig/2._wp))**2._wp)
                end if
            end if
        end if
    end subroutine s_source_spatial_planar

    !> This subroutine calculates the spatial support for a single transducer in 2D, 2D axisymmetric, and 3D
    !! @param ai Acoustic source index
    !! @param sig Sigma value for the Gaussian distribution
    !! @param r Displacement from source to current point
    !! @param source Source term amplitude
    !! @param angle Angle of the source term with respect to the x-axis (for 2D or 2D axisymmetric)
    !! @param xyz_to_r_ratios Ratios of the [xyz]-component of the source term to the magnitude (for 3D)
    subroutine s_source_spatial_transducer(ai, sig, r, source, angle, xyz_to_r_ratios)
        integer, intent(in) :: ai
        real(wp), intent(in) :: sig, r(3)
        real(wp), intent(out) :: source, angle, xyz_to_r_ratios(3)

        real(wp) :: current_angle, angle_half_aperture, dist, norm

        source = 0._wp ! If not affected by transducer
        angle = 0._wp
        xyz_to_r_ratios = 0._wp

        if (support(ai) == 5 .or. support(ai) == 6) then ! 2D or 2D axisymmetric
            current_angle = -atan(r(2)/(foc_length(ai) - r(1)))
            angle_half_aperture = asin((aperture(ai)/2._wp)/(foc_length(ai)))

            if (abs(current_angle) < angle_half_aperture .and. r(1) < foc_length(ai)) then
                dist = foc_length(ai) - sqrt(r(2)**2._wp + (foc_length(ai) - r(1))**2._wp)
                source = 1._wp/(sqrt(2._wp*pi)*sig/2._wp)*exp(-0.5_wp*(dist/(sig/2._wp))**2._wp)
                angle = -atan(r(2)/(foc_length(ai) - r(1)))
            end if

        elseif (support(ai) == 7) then ! 3D
            current_angle = -atan(sqrt(r(2)**2 + r(3)**2)/(foc_length(ai) - r(1)))
            angle_half_aperture = asin((aperture(ai)/2._wp)/(foc_length(ai)))

            if (abs(current_angle) < angle_half_aperture .and. r(1) < foc_length(ai)) then
                dist = foc_length(ai) - sqrt(r(2)**2._wp + r(3)**2._wp + (foc_length(ai) - r(1))**2._wp)
                source = 1._wp/(sqrt(2._wp*pi)*sig/2._wp)*exp(-0.5_wp*(dist/(sig/2._wp))**2._wp)

                norm = sqrt(r(2)**2._wp + r(3)**2._wp + (foc_length(ai) - r(1))**2._wp)
                xyz_to_r_ratios(1) = -(r(1) - foc_length(ai))/norm
                xyz_to_r_ratios(2) = -r(2)/norm
                xyz_to_r_ratios(3) = -r(3)/norm
            end if

        end if
    end subroutine s_source_spatial_transducer

    !> This subroutine calculates the spatial support for multiple transducers in 2D, 2D axisymmetric, and 3D
    !! @param ai Acoustic source index
    !! @param sig Sigma value for the Gaussian distribution
    !! @param r Displacement from source to current point
    !! @param source Source term amplitude
    !! @param angle Angle of the source term with respect to the x-axis (for 2D or 2D axisymmetric)
    !! @param xyz_to_r_ratios Ratios of the [xyz]-component of the source term to the magnitude (for 3D)
    subroutine s_source_spatial_transducer_array(ai, sig, r, source, angle, xyz_to_r_ratios)
        integer, intent(in) :: ai
        real(wp), intent(in) :: sig, r(3)
        real(wp), intent(out) :: source, angle, xyz_to_r_ratios(3)

        integer :: elem, elem_min, elem_max
        real(wp) :: current_angle, angle_half_aperture, angle_per_elem, dist
        real(wp) :: angle_min, angle_max, norm
        real(wp) :: poly_side_length, aperture_element_3D, angle_elem
        real(wp) :: x2, y2, z2, x3, y3, z3, C, f, half_apert, dist_interp_to_elem_center

        if (element_on(ai) == 0) then ! Full transducer
            elem_min = 1
            elem_max = num_elements(ai)
        else ! Transducer element specified
            elem_min = element_on(ai)
            elem_max = element_on(ai)
        end if

        source = 0._wp ! If not affected by any transducer element
        angle = 0._wp
        xyz_to_r_ratios = 0._wp

        if (support(ai) == 9 .or. support(ai) == 10) then ! 2D or 2D axisymmetric
            current_angle = -atan(r(2)/(foc_length(ai) - r(1)))
            angle_half_aperture = asin((aperture(ai)/2._wp)/(foc_length(ai)))
            angle_per_elem = (2._wp*angle_half_aperture - (num_elements(ai) - 1._wp)*element_spacing_angle(ai))/num_elements(ai)
            dist = foc_length(ai) - sqrt(r(2)**2._wp + (foc_length(ai) - r(1))**2._wp)

            do elem = elem_min, elem_max
                angle_max = angle_half_aperture - (element_spacing_angle(ai) + angle_per_elem)*(elem - 1._wp)
                angle_min = angle_max - angle_per_elem

                if (current_angle > angle_min .and. current_angle < angle_max .and. r(1) < foc_length(ai)) then
                    source = exp(-0.5_wp*(dist/(sig/2._wp))**2._wp)/(sqrt(2._wp*pi)*sig/2._wp)
                    angle = current_angle
                    exit ! Assume elements don't overlap
                end if
            end do

        elseif (support(ai) == 11) then ! 3D
            poly_side_length = aperture(ai)*sin(pi/num_elements(ai))
            aperture_element_3D = poly_side_length*element_polygon_ratio(ai)
            f = foc_length(ai)
            half_apert = aperture(ai)/2._wp

            do elem = elem_min, elem_max
                angle_elem = 2._wp*pi*real(elem, wp)/real(num_elements(ai), wp) + rotate_angle(ai)

                ! Point 2 is the elem center
                x2 = f - sqrt(f**2 - half_apert**2)
                y2 = half_apert*cos(angle_elem)
                z2 = half_apert*sin(angle_elem)

                ! Construct a plane normal to the line from the focal point to the elem center,
                ! Point 3 is the intercept of the plane and the line from the focal point to the current location
                C = f**2._wp/((r(1) - f)*(x2 - f) + r(2)*y2 + r(3)*z2) ! Constant for intermediate step
                x3 = C*(r(1) - f) + f
                y3 = C*r(2)
                z3 = C*r(3)

                dist_interp_to_elem_center = sqrt((x2 - x3)**2._wp + (y2 - y3)**2._wp + (z2 - z3)**2._wp)
                if ((dist_interp_to_elem_center < aperture_element_3D/2._wp) .and. (r(1) < f)) then
                    dist = sqrt((x3 - r(1))**2._wp + (y3 - r(2))**2._wp + (z3 - r(3))**2._wp)
                    source = exp(-0.5_wp*(dist/(sig/2._wp))**2._wp)/(sqrt(2._wp*pi)*sig/2._wp)

                    norm = sqrt(r(2)**2._wp + r(3)**2._wp + (f - r(1))**2._wp)
                    xyz_to_r_ratios(1) = -(r(1) - f)/norm
                    xyz_to_r_ratios(2) = -r(2)/norm
                    xyz_to_r_ratios(3) = -r(3)/norm
                end if

            end do

        end if
    end subroutine s_source_spatial_transducer_array

    !> This function performs wavelength to frequency conversion
    !! @param freq_conv_flag Determines if frequency is given or wavelength
    !! @param ai Acoustic source index
    !! @param c Speed of sound
    !! @return frequency_local Converted frequency
    elemental function f_frequency_local(freq_conv_flag, ai, c)
        $:GPU_ROUTINE(parallelism='[seq]')
        logical, intent(in) :: freq_conv_flag
        integer, intent(in) :: ai
        real(wp), intent(in) :: c
        real(wp) :: f_frequency_local

        if (freq_conv_flag) then
            f_frequency_local = c/wavelength(ai)
        else
            f_frequency_local = frequency(ai)
        end if
    end function f_frequency_local

    !> This function performs Gaussian sigma dist to time conversion
    !! @param gauss_conv_flag Determines if sigma_dist is given or sigma_time
    !! @param c Speed of sound
    !! @param ai Acoustic source index
    !! @return gauss_sigma_time_local Converted Gaussian sigma time
    function f_gauss_sigma_time_local(gauss_conv_flag, ai, c)
        $:GPU_ROUTINE(parallelism='[seq]')
        logical, intent(in) :: gauss_conv_flag
        integer, intent(in) :: ai
        real(wp), intent(in) :: c
        real(wp) :: f_gauss_sigma_time_local

        if (gauss_conv_flag) then
            f_gauss_sigma_time_local = gauss_sigma_dist(ai)/c
        else
            f_gauss_sigma_time_local = gauss_sigma_time(ai)
        end if
    end function f_gauss_sigma_time_local

end module m_acoustic_src
