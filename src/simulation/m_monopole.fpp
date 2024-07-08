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

    use m_helper_basic           !< Functions to compare floating point numbers
    ! ==========================================================================
    implicit none
    private; public :: s_initialize_monopole_module, s_monopole_calculations

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(integer, dimension(:), pulse, support)
    !$acc declare link(pulse, support)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :), loc_mono)
    !$acc declare link(loc_mono)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), foc_length, aperture, support_width)
    !$acc declare link(foc_length, aperture, support_width)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), mag, length, npulse, dir, delay)
    !$acc declare link(mag, length, npulse, dir, delay)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), mono_mass_src, mono_e_src)
    !$acc declare link(mono_mass_src, mono_e_src)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), mono_mom_src)
    !$acc declare link(mono_mom_src)

#else
    integer, allocatable, dimension(:) :: pulse, support
    !$acc declare create(pulse, support)

    real(kind(0d0)), allocatable, target, dimension(:, :) :: loc_mono
    !$acc declare create(loc_mono)

    real(kind(0d0)), allocatable, dimension(:) :: foc_length, aperture, support_width
    !$acc declare create(foc_length, aperture, support_width)

    real(kind(0d0)), allocatable, dimension(:) :: mag, length, npulse, dir, delay
    !$acc declare create(mag, length, npulse, dir, delay)

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

        @:ALLOCATE_GLOBAL(mag(1:num_mono), support(1:num_mono), length(1:num_mono), npulse(1:num_mono), pulse(1:num_mono), dir(1:num_mono), delay(1:num_mono), loc_mono(1:3, 1:num_mono), foc_length(1:num_mono), aperture(1:num_mono), support_width(1:num_mono))

        do i = 1, num_mono
            do j = 1, 3
                loc_mono(j, i) = mono(i)%loc(j)
            end do
            mag(i) = mono(i)%mag
            support(i) = mono(i)%support
            length(i) = mono(i)%length
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
            if (f_is_default(mono(i)%support_width)) then
                support_width(i) = 2.5d0
            else
                support_width(i) = mono(i)%support_width
            end if
        end do
        !$acc update device(mag, support, length, npulse, pulse, dir, delay, foc_length, aperture, loc_mono, support_width)

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

        real(kind(0d0)) :: n_tait, B_tait, angle, angle_z

        integer :: ndirs

        real(kind(0d0)) :: the_time, sound
        real(kind(0d0)) :: s2, const_sos, s1

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

                        the_time = t_step*dt
                        if ((the_time >= delay(q)) .or. f_is_default(delay(q))) then
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
                            n_tait = 1.d0/n_tait + 1.d0 !make this the usual little 'gamma'

                            sound = n_tait*(q_prim_vf(E_idx)%sf(j, k, l) + ((n_tait - 1d0)/n_tait)*B_tait)/myRho
                            sound = dsqrt(sound)

                            const_sos = n_tait*(1.01d5 + ((n_tait - 1d0)/n_tait)*B_tait)/myRho
                            const_sos = dsqrt(const_sos)

                            term_index = 2

                            angle = 0.d0
                            angle_z = 0.d0

                            s2 = f_g(the_time, sound, const_sos, q, term_index)* &
                                 f_delta(j, k, l, loc_mono(:, q), length(q), q, angle, angle_z)
                            !s2 = 1d0

                            if (support(q) == 5) then
                                term_index = 1
                                s1 = f_g(the_time, sound, const_sos, q, term_index)* &
                                     f_delta(j, k, l, loc_mono(:, q), length(q), q, angle, angle_z)
                            end if

                            mono_mass_src(j, k, l) = mono_mass_src(j, k, l) + s2/sound

                            if (n == 0) then
                                ! 1D
                                if (dir(q) < -0.1d0) then
                                    !left-going wave
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) - s2
                                else
                                    !right-going wave
                                    mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2
                                end if
                            else if (p == 0) then
                                if (.not. f_is_default(dir(q))) then
                                    ! 2d
                                    !mono_mom_src(1,j,k,l) = s2
                                    !mono_mom_src(2,j,k,l) = s2
                                    if (support(q) == 5) then
                                        mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(angle)
                                        mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*sin(angle)
                                    else
                                        mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(dir(q))
                                        mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*sin(dir(q))
                                    end if
                                end if
                            else
                                ! 3D
                                if (.not. f_is_default(dir(q))) then
                                    if (support(q) == 5) then
                                        mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(angle)
                                        mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*sin(angle)
                                    else if (support(q) == 6) then
                                        ! Cylindrical Coordinate
                                        mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(dir(q))
                                        mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l)
                                    else
                                        mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(dir(q))
                                        mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*sin(dir(q))
                                    end if
                                end if
                            end if

                            if (model_eqns /= 4) then
                                if (support(q) == 5) then
                                    mono_E_src(j, k, l) = mono_E_src(j, k, l) + s1*const_sos**2.d0/(n_tait - 1.d0)
                                else
                                    mono_E_src(j, k, l) = mono_E_src(j, k, l) + s2*sound/(n_tait - 1.d0)
                                end if
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
        !! @param mysos Alternative speed of sound for testing
    function f_g(the_time, sos, mysos, nm, term_index)
        !$acc routine seq
        real(kind(0d0)), intent(in) :: the_time, sos, mysos
        integer, intent(in) :: nm
        integer, intent(in) :: term_index

        real(kind(0d0)) :: period, t0, sigt, pa
        real(kind(0d0)) :: offset
        real(kind(0d0)) :: f_g

        offset = 0d0
        if (.not. f_is_default(delay(nm))) offset = delay(nm)

        if (pulse(nm) == 1) then
            ! Sine wave
            period = length(nm)/sos
            f_g = 0d0
            if (term_index == 1) then
                f_g = mag(nm)*sin((the_time)*2.d0*pi/period)/mysos &
                      + mag(nm)/foc_length(nm)*(1.d0/(2.d0*pi/period)*cos((the_time)*2.d0*pi/period) &
                                                - 1.d0/(2.d0*pi/period))
            elseif (the_time <= (npulse(nm)*period + offset)) then
                f_g = mag(nm)*sin((the_time + offset)*2.d0*pi/period)
            end if
        else if (pulse(nm) == 2) then
            ! Gaussian pulse
            sigt = length(nm)/sos/7.d0
            t0 = 3.5d0*sigt
            f_g = mag(nm)/(dsqrt(2.d0*pi)*sigt)* &
                  dexp(-0.5d0*((the_time - t0)**2.d0)/(sigt**2.d0))
        else if (pulse(nm) == 3) then
            ! Square wave
            sigt = length(nm)/sos
            t0 = 0d0; f_g = 0d0
            if (the_time > t0 .and. the_time < sigt) then
                f_g = mag(nm)
            end if
        end if

    end function f_g

    !> This function give the spatial support of the acoustic source
        !! @param j First coordinate-direction location index
        !! @param k Second coordinate-direction location index
        !! @param l Third coordinate-direction location index
        !! @param mono_loc Nominal source term location
        !! @param mono_leng Length of source term in space
    function f_delta(j, k, l, mono_loc, mono_leng, nm, angle, angle_z)

        !$acc routine seq
        integer, intent(in) :: j, k, l
        real(kind(0d0)), dimension(3), intent(in) :: mono_loc
        real(kind(0d0)), intent(in) :: mono_leng
        integer, intent(in) :: nm
        real(kind(0d0)), intent(out) :: angle
        real(kind(0d0)), intent(out) :: angle_z

        integer :: q
        real(kind(0d0)) :: h, hx, hy, hz
        real(kind(0d0)) :: hx_cyl, hy_cyl, hz_cyl
        real(kind(0d0)) :: hxnew, hynew
        real(kind(0d0)) :: hxnew_cyl, hynew_cyl
        real(kind(0d0)) :: sig
        real(kind(0d0)) :: f_delta

        if (n == 0) then
            sig = dx(j)
            sig = sig*support_width(nm)
        else if (p == 0) then
            sig = maxval((/dx(j), dy(k)/))
            sig = sig*support_width(nm)
        else
            sig = maxval((/dx(j), dy(k), dz(l)/))
            sig = sig*support_width(nm)
        end if

        if (n == 0) then      !1D
            if (support(nm) == 1) then
                ! 1D delta function
                hx = abs(mono_loc(1) - x_cc(j))

                f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                          dexp(-0.5d0*(hx/(sig/2.d0))**2.d0)
            else if (support(nm) == 0) then
                ! Support for all x
                f_delta = 1.d0
            end if
        else if (p == 0) then !2D
            hx = mono_loc(1) - x_cc(j)
            hy = mono_loc(2) - y_cc(k)
            if (support(nm) == 1) then
                ! 2D delta function
                sig = mono_leng/20.d0
                h = dsqrt(hx**2.d0 + hy**2.d0)

                f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                          dexp(-0.5d0*((h/(sig/2.d0))**2.d0))
            else if (support(nm) == 2) then
                !only support for y \pm some value
                if (abs(hy) < length(nm)) then
                    f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                              dexp(-0.5d0*(hx/(sig/2.d0))**2.d0)
                else
                    f_delta = 0d0
                end if
            else if (support(nm) == 3) then
                ! Only support along some line
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)

                ! Rotate actual point by -theta
                hxnew = cos(dir(nm))*hx + sin(dir(nm))*hy
                hynew = -1.d0*sin(dir(nm))*hx + cos(dir(nm))*hy
                if (abs(hynew) < mono_loc(3)/2.d0) then
                    f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2.d0))**2.d0)
                else
                    f_delta = 0d0
                end if
            else if (support(nm) == 4) then
                ! Support for all y
                f_delta = 1.d0/(dsqrt(2.d0*pi)*sig)* &
                          dexp(-0.5d0*(hx/sig)**2.d0)
            else if (support(nm) == 5) then
                ! Support along 'transducer'
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)

                hxnew = foc_length(nm) - dsqrt(hy**2.d0 + (foc_length(nm) - hx)**2.d0)
                if ((abs(hy) < aperture(nm)/2.d0) .and. (hx < foc_length(nm))) then
                    f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2.d0))**2.d0)
                    angle = -atan(hy/(foc_length(nm) - hx))
                else
                    f_delta = 0d0
                end if
            end if
        else !3D

            hx = x_cc(j) - mono_loc(1)
            hy = y_cc(k) - mono_loc(2)
            hz = z_cc(l) - mono_loc(3)
            if (support(nm) == 3) then

                ! Rotate actual point by -theta
                hxnew = cos(dir(nm))*hx + sin(dir(nm))*hy
                hynew = -1.d0*sin(dir(nm))*hx + cos(dir(nm))*hy

                if (abs(hynew) < length(nm)/2. .and. &
                    abs(hz) < length(nm)/2.) then
                    f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2.d0))**2.d0)
                else
                    f_delta = 0d0
                end if
            else if (support(nm) == 4) then
                ! Support for all x,y
                f_delta = 1.d0/(dsqrt(2.d0*pi)*sig)* &
                          dexp(-0.5d0*(hz/sig)**2.d0)
            else if (support(nm) == 5) then
                ! Support along 'transducer'
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)
                hz = z_cc(l) - mono_loc(3)

                hxnew = foc_length(nm) - dsqrt(hy**2.d0 + hz**2.d0 + (foc_length(nm) - hx)**2.d0)
                if ((dsqrt(hy**2.d0 + hz**2.d0) < aperture(nm)/2.d0) .and. &
                    (hx < foc_length(nm))) then

                    f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2.d0))**2.d0)

                    angle = -atan(hy/(foc_length(nm) - hx))
                    angle_z = -atan(hz/(foc_length(nm) - hx))
                else
                    f_delta = 0d0
                end if
            else if (support(nm) == 6) then
                ! Support for cylindrical coordinate system
                !sig = maxval((/dx(j), dy(k)*sin(dz(l)), dz(l)*cos(dz(l))/))
                sig = dx(j)
                sig = sig*support_width(nm)
                hx_cyl = x_cc(j) - mono_loc(1)
                hy_cyl = y_cc(k)*sin(z_cc(l)) - mono_loc(2)
                hz_cyl = y_cc(k)*cos(z_cc(l)) - mono_loc(3)

                ! Rotate actual point by -theta
                hxnew_cyl = cos(dir(nm))*hx_cyl + sin(dir(nm))*hy_cyl
                hynew_cyl = -1.d0*sin(dir(nm))*hx_cyl + cos(dir(nm))*hy_cyl

                if (abs(hynew_cyl) < length(nm)/2. .and. &
                    abs(hz_cyl) < length(nm)/2.) then
                    f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                              dexp(-0.5d0*(hxnew_cyl/(sig/2.d0))**2.d0)
                else
                    f_delta = 0d0
                end if

            end if
        end if

    end function f_delta

end module m_monopole
