!>
!! @file m_derived_variables.f90
!! @brief Contains module m_derived_variables

!> @brief This module features subroutines that allow for the derivation of
!!              numerous flow variables from the conservative and primitive ones.
!!              Currently, the available derived variables include the unadvected
!!              volume fraction, specific heat ratio, liquid stiffness, speed of
!!              sound, vorticity and the numerical Schlieren function.
module m_derived_variables

    ! Dependencies =============================================================
    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_helper

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_data_output           !< Data output module

    use m_time_steppers         !< Time-stepping algorithms
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_derived_variables_module, &
 s_initialize_derived_variables, &
 s_compute_derived_variables, &
 s_finalize_derived_variables_module

    !> @name Finite-difference coefficients
    !! Finite-difference (fd) coefficients in x-, y- and z-coordinate directions.
    !! Note that because sufficient boundary information is available for all the
    !! active coordinate directions, the centered family of the finite-difference
    !! schemes is used.
    !> @{
    real(kind(0d0)), public, allocatable, dimension(:, :) :: fd_coeff_x
    real(kind(0d0)), public, allocatable, dimension(:, :) :: fd_coeff_y
    real(kind(0d0)), public, allocatable, dimension(:, :) :: fd_coeff_z
    !> @}

contains

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    subroutine s_initialize_derived_variables_module() ! ----------------------

        ! Allocating the variables which will store the coefficients of the
        ! centered family of finite-difference schemes. Note that sufficient
        ! space is allocated so that the coefficients up to any chosen order
        ! of accuracy may be bookkept. However, if higher than fourth-order
        ! accuracy coefficients are wanted, the formulae required to compute
        ! these coefficients will have to be implemented in the subroutine
        ! s_compute_finite_difference_coefficients.

        ! Allocating centered finite-difference coefficients
        if (probe_wrt) then
            allocate (fd_coeff_x(-fd_number:fd_number, 0:m))
            if (n > 0) then
                allocate (fd_coeff_y(-fd_number:fd_number, 0:n))
                if (p > 0) then
                    allocate (fd_coeff_z(-fd_number:fd_number, 0:p))
                end if
            end if
        end if

    end subroutine s_initialize_derived_variables_module ! --------------------

    !> Allocate and open derived variables. Computing FD coefficients.
    subroutine s_initialize_derived_variables() ! -----------------------------

        ! Opening and writing header of CoM and flow probe files
        if (proc_rank == 0) then
            if (any(com_wrt)) then
                call s_open_com_files()
            end if
            if (any(cb_wrt)) then
                call s_open_cb_files()
            end if
            if (probe_wrt) then
                call s_open_probe_files()
            end if
        end if

        ! Computing centered finite difference coefficients
        if (probe_wrt) then
            call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, &
                                                             fd_number, fd_order)
            if (n > 0) then
                call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, &
                                                                 fd_number, fd_order)
                if (p > 0) then
                    call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z, buff_size, &
                                                                     fd_number, fd_order)
                end if
            end if
        end if

    end subroutine s_initialize_derived_variables ! -----------------------------

    !> Writes coherent body information, communication files, and probes.
        !!  @param t_step Current time-step
    subroutine s_compute_derived_variables(t_step) ! -----------------------

        integer, intent(IN) :: t_step

        integer :: i, j, k !< Generic loop iterators

        ! IF ((ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) .AND. (t_step > t_step_start + 2)) THEN
        if ((any(com_wrt) .or. any(cb_wrt) .or. probe_wrt)) then
            if (any(com_wrt)) then
                call s_derive_center_of_mass(q_prim_ts(0)%vf, &
                                             q_prim_ts(1)%vf, &
                                             q_prim_ts(2)%vf, &
                                             q_prim_ts(3)%vf, &
                                             q_com)
                call s_derive_higher_moments(q_prim_ts(0)%vf, moments)
                call s_write_com_files(t_step, q_com, moments)
            end if

            if (any(cb_wrt)) then
                call s_derive_fluid_bounds(q_prim_ts(0)%vf, bounds)
                call s_derive_coherent_body(q_prim_ts(0)%vf, cb_mass)
                call s_derive_centerline(q_prim_ts(0)%vf, cntrline)
                call s_write_cb_files(t_step, cb_mass, bounds, cntrline)
            end if

            if (probe_wrt) then
                call s_derive_acceleration_component(1, q_prim_ts(0)%vf, &
                                                     q_prim_ts(1)%vf, &
                                                     q_prim_ts(2)%vf, &
                                                     q_prim_ts(3)%vf, &
                                                     x_accel)
                if (n > 0) then
                    call s_derive_acceleration_component(2, q_prim_ts(0)%vf, &
                                                         q_prim_ts(1)%vf, &
                                                         q_prim_ts(2)%vf, &
                                                         q_prim_ts(3)%vf, &
                                                         y_accel)
                    if (p > 0) then
                        call s_derive_acceleration_component(3, q_prim_ts(0)%vf, &
                                                             q_prim_ts(1)%vf, &
                                                             q_prim_ts(2)%vf, &
                                                             q_prim_ts(3)%vf, &
                                                             z_accel)
                    end if
                end if

                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            if (p > 0) then
                                accel_mag(i, j, k) = sqrt(x_accel(i, j, k)**2d0 + &
                                                          y_accel(i, j, k)**2d0 + &
                                                          z_accel(i, j, k)**2d0)
                            elseif (n > 0) then
                                accel_mag(i, j, k) = sqrt(x_accel(i, j, k)**2d0 + &
                                                          y_accel(i, j, k)**2d0)
                            else
                                accel_mag(i, j, k) = x_accel(i, j, k)
                            end if
                        end do
                    end do
                end do

                call s_write_probe_files(t_step, q_cons_ts(1)%vf, accel_mag)
            end if
        end if

    end subroutine s_compute_derived_variables ! ---------------------------

    !> This subroutine receives as inputs the indicator of the
        !!      component of the acceleration that should be outputted and
        !!      the primitive variables. From those inputs, it proceeds
        !!      to calculate values of the desired acceleration component,
        !!      which are subsequently stored in derived flow quantity
        !!      storage variable, q_sf.
        !!  @param i Acceleration component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_prim_vf1 Primitive variables
        !!  @param q_prim_vf2 Primitive variables
        !!  @param q_prim_vf3 Primitive variables
        !!  @param q_sf Acceleration component
    subroutine s_derive_acceleration_component(i, q_prim_vf, q_prim_vf1, &
                                               q_prim_vf2, q_prim_vf3, q_sf) ! ----------

        integer, intent(IN) :: i

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf1
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf2
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf3

        real(kind(0d0)), dimension(0:m, 0:n, 0:p), intent(OUT) :: q_sf

        integer :: j, k, l, r !< Generic loop iterators

        ! Computing the acceleration component in the x-coordinate direction
        if (i == 1) then
            do l = 0, p
                do k = 0, n
                    do j = 0, m

                        q_sf(j, k, l) = (11d0*q_prim_vf(mom_idx%beg)%sf(j, k, l) &
                                         - 18d0*q_prim_vf1(mom_idx%beg)%sf(j, k, l) &
                                         + 9d0*q_prim_vf2(mom_idx%beg)%sf(j, k, l) &
                                         - 2d0*q_prim_vf3(mom_idx%beg)%sf(j, k, l))/(6d0*dt)

                        do r = -fd_number, fd_number
                            if (n == 0) then ! 1D simulation
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf(mom_idx%beg)%sf(r + j, k, l)
                            elseif (p == 0) then ! 2D simulation
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf(mom_idx%beg)%sf(r + j, k, l) &
                                                + q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                q_prim_vf(mom_idx%beg)%sf(j, r + k, l)
                            else ! 3D simulation
                                if (grid_geometry == 3) then
                                    q_sf(j, k, l) = q_sf(j, k, l) &
                                                    + q_prim_vf(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                    q_prim_vf(mom_idx%beg)%sf(r + j, k, l) &
                                                    + q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                    q_prim_vf(mom_idx%beg)%sf(j, r + k, l) &
                                                    + q_prim_vf(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                    q_prim_vf(mom_idx%beg)%sf(j, k, r + l)/y_cc(k)
                                else
                                    q_sf(j, k, l) = q_sf(j, k, l) &
                                                    + q_prim_vf(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                    q_prim_vf(mom_idx%beg)%sf(r + j, k, l) &
                                                    + q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                    q_prim_vf(mom_idx%beg)%sf(j, r + k, l) &
                                                    + q_prim_vf(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                    q_prim_vf(mom_idx%beg)%sf(j, k, r + l)
                                end if
                            end if
                        end do
                    end do
                end do
            end do

            ! Computing the acceleration component in the y-coordinate direction
        elseif (i == 2) then
            do l = 0, p
                do k = 0, n
                    do j = 0, m

                        q_sf(j, k, l) = (11d0*q_prim_vf(mom_idx%beg + 1)%sf(j, k, l) &
                                         - 18d0*q_prim_vf1(mom_idx%beg + 1)%sf(j, k, l) &
                                         + 9d0*q_prim_vf2(mom_idx%beg + 1)%sf(j, k, l) &
                                         - 2d0*q_prim_vf3(mom_idx%beg + 1)%sf(j, k, l))/(6d0*dt)

                        do r = -fd_number, fd_number
                            if (p == 0) then ! 2D simulation
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf(mom_idx%beg + 1)%sf(r + j, k, l) &
                                                + q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                q_prim_vf(mom_idx%beg + 1)%sf(j, r + k, l)
                            else ! 3D simulation
                                if (grid_geometry == 3) then
                                    q_sf(j, k, l) = q_sf(j, k, l) &
                                                    + q_prim_vf(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                    q_prim_vf(mom_idx%beg + 1)%sf(r + j, k, l) &
                                                    + q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                    q_prim_vf(mom_idx%beg + 1)%sf(j, r + k, l) &
                                                    + q_prim_vf(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                    q_prim_vf(mom_idx%beg + 1)%sf(j, k, r + l)/y_cc(k) &
                                                    - (q_prim_vf(mom_idx%end)%sf(j, k, l)**2d0)/y_cc(k)
                                else
                                    q_sf(j, k, l) = q_sf(j, k, l) &
                                                    + q_prim_vf(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                    q_prim_vf(mom_idx%beg + 1)%sf(r + j, k, l) &
                                                    + q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                    q_prim_vf(mom_idx%beg + 1)%sf(j, r + k, l) &
                                                    + q_prim_vf(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                    q_prim_vf(mom_idx%beg + 1)%sf(j, k, r + l)
                                end if
                            end if
                        end do
                    end do
                end do
            end do

            ! Computing the acceleration component in the z-coordinate direction
        else
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_sf(j, k, l) = (11d0*q_prim_vf(mom_idx%end)%sf(j, k, l) &
                                         - 18d0*q_prim_vf1(mom_idx%end)%sf(j, k, l) &
                                         + 9d0*q_prim_vf2(mom_idx%end)%sf(j, k, l) &
                                         - 2d0*q_prim_vf3(mom_idx%end)%sf(j, k, l))/(6d0*dt)

                        do r = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf(mom_idx%end)%sf(r + j, k, l) &
                                                + q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                q_prim_vf(mom_idx%end)%sf(j, r + k, l) &
                                                + q_prim_vf(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                q_prim_vf(mom_idx%end)%sf(j, k, r + l)/y_cc(k) &
                                                + (q_prim_vf(mom_idx%end)%sf(j, k, l)* &
                                                   q_prim_vf(mom_idx%beg + 1)%sf(j, k, l))/y_cc(k)
                            else
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf(mom_idx%end)%sf(r + j, k, l) &
                                                + q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                q_prim_vf(mom_idx%end)%sf(j, r + k, l) &
                                                + q_prim_vf(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                q_prim_vf(mom_idx%end)%sf(j, k, r + l)
                            end if
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_derive_acceleration_component ! --------------------------

    !> This subroutine is used together with the volume fraction
        !!      model and when called upon, it computes the location of
        !!      of the center of mass for each fluid from the inputted
        !!      primitive variables, q_prim_vf. The computed location
        !!      is then written to a formatted data file by the root
        !!      process.
        !!  @param q_prim_vf Primitive variables
        !!  @param q_prim_vf1 Primitive variables
        !!  @param q_prim_vf2 Primitive variables
        !!  @param q_prim_vf3 Primitive variables
        !!  @param q_com Mass,x-location,y-location,z-location,x-velocity,y-velocity,z-velocity,
        !!  x-acceleration, y-acceleration, z-acceleration, weighted
    subroutine s_derive_center_of_mass(q_prim_vf, q_prim_vf1, q_prim_vf2, q_prim_vf3, q_com)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf1
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf2
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf3
        real(kind(0d0)), dimension(num_fluids, 10), intent(INOUT) :: q_com

        real(kind(0d0)) :: xbeg, xend, ybeg, yend, zbeg, zend !<
            !! Maximum and minimum values of cell boundaries in each direction used in check for
            !! reflective BC in computation of center of mass

        integer :: i, j, k, l !< Generic loop iterators

        real(kind(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce

        real(kind(0d0)) :: dV !< Discrete cell volume

        real(kind(0d0)) :: cart_u_x, cart_u_y, &
                           cart_u_x1, cart_u_y1, &
                           cart_u_x2, cart_u_y2, &
                           cart_u_x3, cart_u_y3 !<
            !! Cartesian velocities

        if (n == 0) then !1D simulation

            do i = 1, num_fluids !Loop over individual fluids
                if (com_wrt(i)) then
                    q_com(i, :) = 0d0
                    do l = 0, p !Loop over grid
                        do k = 0, n
                            do j = 0, m

                                dV = dx(j)

                                ! Mass
                                q_com(i, 1) = q_com(i, 1) + q_prim_vf(i)%sf(j, k, l)*dV
                                ! x-location weighted
                                q_com(i, 2) = q_com(i, 2) + q_prim_vf(i)%sf(j, k, l)*dV*x_cc(j)
                                ! x-velocity weighted
                                q_com(i, 5) = q_com(i, 5) + q_prim_vf(i)%sf(j, k, l)*dV*q_prim_vf(mom_idx%beg)%sf(j, k, l)
                                ! x-acceleration weighted
                                q_com(i, 8) = q_com(i, 8) + dV*(11d0*(q_prim_vf(i)%sf(j, k, l) &
                                                                      *q_prim_vf(mom_idx%beg)%sf(j, k, l)) &
                                                            - 18d0*(q_prim_vf1(i)%sf(j, k, l)*q_prim_vf1(mom_idx%beg)%sf(j, k, l)) &
                                                             + 9d0*(q_prim_vf2(i)%sf(j, k, l)*q_prim_vf2(mom_idx%beg)%sf(j, k, l)) &
                                                     - 2d0*(q_prim_vf3(i)%sf(j, k, l)*q_prim_vf3(mom_idx%beg)%sf(j, k, l)))/(6d0*dt)
                            end do
                        end do
                    end do
                    ! Sum all components across all processors using MPI_ALLREDUCE
                    if (num_procs > 1) then
                        tmp = q_com(i, 1)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 1))
                        tmp = q_com(i, 2)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 2))
                        tmp = q_com(i, 5)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 5))
                        tmp = q_com(i, 8)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 8))
                    end if

                    ! Compute quotients
                    q_com(i, 2) = q_com(i, 2)/q_com(i, 1)
                    q_com(i, 5) = q_com(i, 5)/q_com(i, 1)
                    q_com(i, 8) = q_com(i, 8)/q_com(i, 1)
                end if
            end do

        elseif (p == 0) then !2D simulation

            do i = 1, num_fluids !Loop over individual fluids
                if (com_wrt(i)) then
                    q_com(i, :) = 0d0
                    do l = 0, p !Loop over grid
                        do k = 0, n
                            do j = 0, m

                                dV = dx(j)*dy(k)

                                ! Mass
                                q_com(i, 1) = q_com(i, 1) + q_prim_vf(i)%sf(j, k, l)*dV
                                ! x-location weighted
                                q_com(i, 2) = q_com(i, 2) + q_prim_vf(i)%sf(j, k, l)*dV*x_cc(j)
                                ! y-location weighted
                                q_com(i, 3) = q_com(i, 3) + q_prim_vf(i)%sf(j, k, l)*dV*y_cc(k)
                                ! x-velocity weighted
                                q_com(i, 5) = q_com(i, 5) + q_prim_vf(i)%sf(j, k, l)*dV*q_prim_vf(mom_idx%beg)%sf(j, k, l)
                                ! y-velocity weighted
                                q_com(i, 6) = q_com(i, 6) + q_prim_vf(i)%sf(j, k, l)*dV*q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)
                                ! x-acceleration weighted
                                q_com(i, 8) = q_com(i, 8) + dV* &
                                              (11d0*(q_prim_vf(i)%sf(j, k, l)*q_prim_vf(mom_idx%beg)%sf(j, k, l)) &
                                               - 18d0*(q_prim_vf1(i)%sf(j, k, l)*q_prim_vf1(mom_idx%beg)%sf(j, k, l)) &
                                               + 9d0*(q_prim_vf2(i)%sf(j, k, l)*q_prim_vf2(mom_idx%beg)%sf(j, k, l)) &
                                               - 2d0*(q_prim_vf3(i)%sf(j, k, l)*q_prim_vf3(mom_idx%beg)%sf(j, k, l)))/(6d0*dt)
                                ! y-acceleration weighted
                                q_com(i, 9) = q_com(i, 9) + dV* &
                                              (11d0*(q_prim_vf(i)%sf(j, k, l)*q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)) &
                                               - 18d0*(q_prim_vf1(i)%sf(j, k, l)*q_prim_vf1(mom_idx%beg + 1)%sf(j, k, l)) &
                                               + 9d0*(q_prim_vf2(i)%sf(j, k, l)*q_prim_vf2(mom_idx%beg + 1)%sf(j, k, l)) &
                                               - 2d0*(q_prim_vf3(i)%sf(j, k, l)*q_prim_vf3(mom_idx%beg + 1)%sf(j, k, l)))/(6d0*dt)
                            end do
                        end do
                    end do
                    ! Sum all components across all processors using MPI_ALLREDUCE
                    if (num_procs > 1) then
                        tmp = q_com(i, 1)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 1))
                        tmp = q_com(i, 2)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 2))
                        tmp = q_com(i, 3)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 3))
                        tmp = q_com(i, 5)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 5))
                        tmp = q_com(i, 6)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 6))
                        tmp = q_com(i, 8)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 8))
                        tmp = q_com(i, 9)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 9))
                    end if

                    ! Compute quotients
                    q_com(i, 2) = q_com(i, 2)/q_com(i, 1)
                    q_com(i, 3) = q_com(i, 3)/q_com(i, 1)
                    q_com(i, 5) = q_com(i, 5)/q_com(i, 1)
                    q_com(i, 6) = q_com(i, 6)/q_com(i, 1)
                    q_com(i, 8) = q_com(i, 8)/q_com(i, 1)
                    q_com(i, 9) = q_com(i, 9)/q_com(i, 1)
                end if
            end do

        else !3D simulation

            do i = 1, num_fluids !Loop over individual fluids
                if (com_wrt(i)) then
                    q_com(i, :) = 0d0
                    do l = 0, p !Loop over grid
                        do k = 0, n
                            do j = 0, m
                                if (grid_geometry == 3) then

                                    dV = (2d0*y_cb(k - 1)*dy(k) + dy(k)**2d0)/2d0*dx(j)*dz(l)
                                    cart_u_x = q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*cos(z_cc(l)) - &
                                               q_prim_vf(mom_idx%end)%sf(j, k, l)*sin(z_cc(l))
                                    cart_u_y = q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*sin(z_cc(l)) + &
                                               q_prim_vf(mom_idx%end)%sf(j, k, l)*cos(z_cc(l))
                                    cart_u_x1 = q_prim_vf1(mom_idx%beg + 1)%sf(j, k, l)*cos(z_cc(l)) - &
                                                q_prim_vf1(mom_idx%end)%sf(j, k, l)*sin(z_cc(l))
                                    cart_u_y1 = q_prim_vf1(mom_idx%beg + 1)%sf(j, k, l)*sin(z_cc(l)) + &
                                                q_prim_vf1(mom_idx%end)%sf(j, k, l)*cos(z_cc(l))
                                    cart_u_x2 = q_prim_vf2(mom_idx%beg + 1)%sf(j, k, l)*cos(z_cc(l)) - &
                                                q_prim_vf2(mom_idx%end)%sf(j, k, l)*sin(z_cc(l))
                                    cart_u_y2 = q_prim_vf2(mom_idx%beg + 1)%sf(j, k, l)*sin(z_cc(l)) + &
                                                q_prim_vf2(mom_idx%end)%sf(j, k, l)*cos(z_cc(l))
                                    cart_u_x3 = q_prim_vf3(mom_idx%beg + 1)%sf(j, k, l)*cos(z_cc(l)) - &
                                                q_prim_vf3(mom_idx%end)%sf(j, k, l)*sin(z_cc(l))
                                    cart_u_y3 = q_prim_vf3(mom_idx%beg + 1)%sf(j, k, l)*sin(z_cc(l)) + &
                                                q_prim_vf3(mom_idx%end)%sf(j, k, l)*cos(z_cc(l))

                                    ! Mass
                                    q_com(i, 1) = q_com(i, 1) + q_prim_vf(i)%sf(j, k, l)*dV
                                    ! x-location weighted
                                    q_com(i, 2) = q_com(i, 2) + q_prim_vf(i)%sf(j, k, l)*dV*y_cc(k)*cos(z_cc(l))
                                    ! y-location weighted
                                    q_com(i, 3) = q_com(i, 3) + q_prim_vf(i)%sf(j, k, l)*dV*y_cc(k)*sin(z_cc(l))
                                    ! z-location weighted
                                    q_com(i, 4) = q_com(i, 4) + q_prim_vf(i)%sf(j, k, l)*dV*x_cc(j)
                                    ! x-velocity weighted
                                    q_com(i, 5) = q_com(i, 5) + q_prim_vf(i)%sf(j, k, l)*dV*cart_u_x
                                    ! y-velocity weighted
                                    q_com(i, 6) = q_com(i, 6) + q_prim_vf(i)%sf(j, k, l)*dV*cart_u_y
                                    ! z-velocity weighted
                                    q_com(i, 7) = q_com(i, 7) + q_prim_vf(i)%sf(j, k, l)*dV*q_prim_vf(mom_idx%beg)%sf(j, k, l)
                                    ! x-acceleration weighted
                                    q_com(i, 8) = q_com(i, 8) + dV* &
                                                  (11d0*(q_prim_vf(i)%sf(j, k, l)*cart_u_x) &
                                                   - 18d0*(q_prim_vf1(i)%sf(j, k, l)*cart_u_x1) &
                                                   + 9d0*(q_prim_vf2(i)%sf(j, k, l)*cart_u_x2) &
                                                   - 2d0*(q_prim_vf3(i)%sf(j, k, l)*cart_u_x3))/(6d0*dt)
                                    ! y-acceleration weighted
                                    q_com(i, 9) = q_com(i, 9) + dV* &
                                                  (11d0*(q_prim_vf(i)%sf(j, k, l)*cart_u_y) &
                                                   - 18d0*(q_prim_vf1(i)%sf(j, k, l)*cart_u_y1) &
                                                   + 9d0*(q_prim_vf2(i)%sf(j, k, l)*cart_u_y2) &
                                                   - 2d0*(q_prim_vf3(i)%sf(j, k, l)*cart_u_y3))/(6d0*dt)
                                    ! z-acceleration weighted
                                    q_com(i, 10) = q_com(i, 10) + dV* &
                                                   (11d0*(q_prim_vf(i)%sf(j, k, l)*q_prim_vf(mom_idx%beg)%sf(j, k, l)) &
                                                    - 18d0*(q_prim_vf1(i)%sf(j, k, l)*q_prim_vf1(mom_idx%beg)%sf(j, k, l)) &
                                                    + 9d0*(q_prim_vf2(i)%sf(j, k, l)*q_prim_vf2(mom_idx%beg)%sf(j, k, l)) &
                                                    - 2d0*(q_prim_vf3(i)%sf(j, k, l)*q_prim_vf3(mom_idx%beg)%sf(j, k, l)))/(6d0*dt)
                                else

                                    dV = dx(j)*dy(k)*dz(l)

                                    ! Mass
                                    q_com(i, 1) = q_com(i, 1) + q_prim_vf(i)%sf(j, k, l)*dV
                                    ! x-location weighted
                                    q_com(i, 2) = q_com(i, 2) + q_prim_vf(i)%sf(j, k, l)*dV*x_cc(j)
                                    ! y-location weighted
                                    q_com(i, 3) = q_com(i, 3) + q_prim_vf(i)%sf(j, k, l)*dV*y_cc(k)
                                    ! z-location weighted
                                    q_com(i, 4) = q_com(i, 4) + q_prim_vf(i)%sf(j, k, l)*dV*z_cc(l)
                                    ! x-velocity weighted
                                    q_com(i, 5) = q_com(i, 5) + q_prim_vf(i)%sf(j, k, l)*dV*q_prim_vf(mom_idx%beg)%sf(j, k, l)
                                    ! y-velocity weighted
                                    q_com(i, 6) = q_com(i, 6) + q_prim_vf(i)%sf(j, k, l)*dV*q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)
                                    ! z-velocity weighted
                                    q_com(i, 7) = q_com(i, 7) + q_prim_vf(i)%sf(j, k, l)*dV*q_prim_vf(mom_idx%end)%sf(j, k, l)
                                    ! x-acceleration weighted
                                    q_com(i, 8) = q_com(i, 8) + dV* &
                                                  (11d0*(q_prim_vf(i)%sf(j, k, l)*q_prim_vf(mom_idx%beg)%sf(j, k, l)) &
                                                   - 18d0*(q_prim_vf1(i)%sf(j, k, l)*q_prim_vf1(mom_idx%beg)%sf(j, k, l)) &
                                                   + 9d0*(q_prim_vf2(i)%sf(j, k, l)*q_prim_vf2(mom_idx%beg)%sf(j, k, l)) &
                                                   - 2d0*(q_prim_vf3(i)%sf(j, k, l)*q_prim_vf3(mom_idx%beg)%sf(j, k, l)))/(6d0*dt)
                                    ! y-acceleration weighted
                                    q_com(i, 9) = q_com(i, 9) + dV* &
                                                  (11d0*(q_prim_vf(i)%sf(j, k, l)*q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)) &
                                                   - 18d0*(q_prim_vf1(i)%sf(j, k, l)*q_prim_vf1(mom_idx%beg + 1)%sf(j, k, l)) &
                                                   + 9d0*(q_prim_vf2(i)%sf(j, k, l)*q_prim_vf2(mom_idx%beg + 1)%sf(j, k, l)) &
                                                 - 2d0*(q_prim_vf3(i)%sf(j, k, l)*q_prim_vf3(mom_idx%beg + 1)%sf(j, k, l)))/(6d0*dt)
                                    ! z-acceleration weighted
                                    q_com(i, 10) = q_com(i, 10) + dV* &
                                                   (11d0*(q_prim_vf(i)%sf(j, k, l)*q_prim_vf(mom_idx%end)%sf(j, k, l)) &
                                                    - 18d0*(q_prim_vf1(i)%sf(j, k, l)*q_prim_vf1(mom_idx%end)%sf(j, k, l)) &
                                                    + 9d0*(q_prim_vf2(i)%sf(j, k, l)*q_prim_vf2(mom_idx%end)%sf(j, k, l)) &
                                                    - 2d0*(q_prim_vf3(i)%sf(j, k, l)*q_prim_vf3(mom_idx%end)%sf(j, k, l)))/(6d0*dt)
                                end if
                            end do
                        end do
                    end do
                    ! Sum all components across all processors using MPI_ALLREDUCE
                    if (num_procs > 1) then
                        tmp = q_com(i, 1)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 1))
                        tmp = q_com(i, 2)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 2))
                        tmp = q_com(i, 3)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 3))
                        tmp = q_com(i, 4)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 4))
                        tmp = q_com(i, 5)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 5))
                        tmp = q_com(i, 6)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 6))
                        tmp = q_com(i, 7)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 7))
                        tmp = q_com(i, 8)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 8))
                        tmp = q_com(i, 9)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 9))
                        tmp = q_com(i, 10)
                        call s_mpi_allreduce_sum(tmp, q_com(i, 10))
                    end if

                    ! Compute quotients
                    q_com(i, 2) = q_com(i, 2)/q_com(i, 1)
                    q_com(i, 3) = q_com(i, 3)/q_com(i, 1)
                    q_com(i, 4) = q_com(i, 4)/q_com(i, 1)
                    q_com(i, 5) = q_com(i, 5)/q_com(i, 1)
                    q_com(i, 6) = q_com(i, 6)/q_com(i, 1)
                    q_com(i, 7) = q_com(i, 7)/q_com(i, 1)
                    q_com(i, 8) = q_com(i, 8)/q_com(i, 1)
                    q_com(i, 9) = q_com(i, 9)/q_com(i, 1)
                    q_com(i, 10) = q_com(i, 10)/q_com(i, 1)
                end if
            end do

        end if

        ! Find computational domain boundaries
        if (num_procs > 1) then
            call s_mpi_allreduce_min(minval(x_cb(-1:m)), xbeg)
            call s_mpi_allreduce_max(maxval(x_cb(-1:m)), xend)
            if (n > 0) then
                call s_mpi_allreduce_min(minval(y_cb(-1:n)), ybeg)
                call s_mpi_allreduce_max(maxval(y_cb(-1:n)), yend)
                if (p > 0) then
                    call s_mpi_allreduce_min(minval(z_cb(-1:p)), zbeg)
                    call s_mpi_allreduce_max(maxval(z_cb(-1:p)), zend)
                end if
            end if
        else
            xbeg = minval(x_cb(-1:m))
            xend = maxval(x_cb(-1:m))
            if (n > 0) then
                ybeg = minval(y_cb(-1:n))
                yend = maxval(y_cb(-1:n))
                if (p > 0) then
                    zbeg = minval(z_cb(-1:p))
                    zend = maxval(z_cb(-1:p))
                end if
            end if
        end if

        do i = 1, num_fluids
            if (com_wrt(i)) then
                ! Check for reflective BC in x-direction
                if (bc_x_glb%beg == -2) then
                    q_com(i, 1) = q_com(i, 1)*2d0
                    q_com(i, 2) = xbeg
                    q_com(i, 5) = 0d0
                    q_com(i, 8) = 0d0
                elseif (bc_x_glb%end == -2) then
                    q_com(i, 1) = q_com(i, 1)*2d0
                    q_com(i, 2) = xend
                    q_com(i, 5) = 0d0
                    q_com(i, 8) = 0d0
                end if
                if (n > 0) then
                    ! Check for reflective BC in y-direction
                    if (bc_y_glb%beg == -2) then
                        q_com(i, 1) = q_com(i, 1)*2d0
                        q_com(i, 3) = ybeg
                        q_com(i, 6) = 0d0
                        q_com(i, 9) = 0d0
                    elseif (bc_y_glb%end == -2) then
                        q_com(i, 1) = q_com(i, 1)*2d0
                        q_com(i, 3) = yend
                        q_com(i, 6) = 0d0
                        q_com(i, 9) = 0d0
                    end if
                    if (p > 0) then
                        ! Check for reflective BC in z-direction
                        if (bc_z_glb%beg == -2) then
                            q_com(i, 1) = q_com(i, 1)*2d0
                            q_com(i, 4) = zbeg
                            q_com(i, 7) = 0d0
                            q_com(i, 10) = 0d0
                        elseif (bc_z_glb%end == -2) then
                            q_com(i, 1) = q_com(i, 1)*2d0
                            q_com(i, 4) = zend
                            q_com(i, 7) = 0d0
                            q_com(i, 10) = 0d0
                        end if

                    end if
                end if
            end if
        end do

    end subroutine s_derive_center_of_mass ! ----------------------------------

    !>  Subroutine to compute the higher moments in an attempt to find
        !!      the maximal size of the droplet
        !!  @param q_prim_vf Primitive variables
        !!  @param moments Higher moments (2 lateral directions, 5 moment orders)
    subroutine s_derive_higher_moments(q_prim_vf, moments)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_fluids, 2, 5), intent(INOUT) :: moments

        integer :: i, r !< Generic loop iterators

        ! Using the global boundary conditions, determine method of computing
        ! higher moments for y-direction
        if (n > 0) then
            if ((bc_y_glb%beg /= -2) .and. (bc_y_glb%end /= -2)) then
                ! Non-symmetric moments
                call s_non_symmetric_moments(q_prim_vf, moments, 1)
            elseif (((bc_y_glb%beg == -2) .and. (bc_y_glb%end == -2)) &
                    .or. &
                    ((bc_y_glb%beg == -1) .and. (bc_y_glb%end == -1))) then
                print '(A)', 'Periodic boundary conditions in y-direction. '// &
                    'Cannot compute higher moments. Exiting...'
                call s_mpi_abort()
            else
                call s_symmetric_moments(q_prim_vf, moments, 1)
            end if

            if (p > 0) then
                if ((bc_z_glb%beg /= -2) .and. (bc_z_glb%end /= -2)) then
                    ! Non-symmetric moments
                    call s_non_symmetric_moments(q_prim_vf, moments, 2)
                elseif (((bc_z_glb%beg == -2) .and. (bc_z_glb%end == -2)) &
                        .or. &
                        ((bc_z_glb%beg == -1) .and. (bc_z_glb%end == -1))) then
                    print '(A)', 'Periodic boundary conditions in z-direction. '// &
                        'Cannot compute higher moments. Exiting...'
                    call s_mpi_abort()
                else
                    call s_symmetric_moments(q_prim_vf, moments, 2)
                end if
            end if
        else !1D simulation
            do i = 1, num_fluids !Loop over individual fluids
                if (com_wrt(i)) then
                    do r = 1, 5
                        if (moment_order(r) /= dflt_int) then
                            moments(i, :, r) = 0d0
                        else
                            moments(i, :, r) = dflt_real
                        end if
                    end do
                end if
            end do
        end if

    end subroutine s_derive_higher_moments ! -----------------------------------------

    !> Compute non-symmetric moments
        !! @param q_prim_vf Primitive variables
        !! @param moments Higher moments(2 lateral directions, 5 moment orders)
        !! @param dir Current lateral direction
    subroutine s_non_symmetric_moments(q_prim_vf, moments, dir) ! ---------------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_fluids, 2, 5), intent(INOUT) :: moments
        integer, intent(IN) :: dir

        real(kind(0d0)), dimension(num_fluids, 5) :: pos_numer, neg_numer, pos_denom, neg_denom !<
            !! Numerator and denominator place holders for computation

        real(kind(0d0)) :: numer_weight     !< Numerator weight
        real(kind(0d0)) :: main_term        !< Constant term in both numerator and denominator
        real(kind(0d0)) :: dV               !< Discrete cell volume
        real(kind(0d0)) :: cart_x, cart_y   !< Cartesian x- and y-locations
        integer :: i, j, k, l, r    !< Generic loop iterators
        real(kind(0d0)) :: tmp  !< Temporary variable to store quantity for mpi_allreduce

        do i = 1, num_fluids
            if (com_wrt(i)) then
                pos_numer(i, :) = 0d0
                neg_numer(i, :) = 0d0
                pos_denom(i, :) = 0d0
                neg_denom(i, :) = 0d0
                do r = 1, 5
                    if (moment_order(r) /= dflt_int) then
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    if (q_prim_vf(i)%sf(j, k, l) > 5d-1) then
                                        if (grid_geometry == 3) then
                                            dV = (2d0*y_cb(k - 1)*dy(k) + dy(k)**2d0)/2d0*dx(j)*dz(l)
                                            cart_x = y_cc(k)*cos(z_cc(l))
                                            cart_y = y_cc(k)*sin(z_cc(l))
                                            main_term = q_prim_vf(i + E_idx)%sf(j, k, l)* &
                                                        (1d0 - q_prim_vf(i + E_idx)%sf(j, k, l))*dV
                                            if ((dir == 1) .and. (cart_x >= 0d0)) then
                                                numer_weight = cart_x**moment_order(r)

                                                pos_numer(i, r) = pos_numer(i, r) + numer_weight*main_term
                                                pos_denom(i, r) = pos_denom(i, r) + main_term
                                            elseif ((dir == 1) .and. (cart_x < 0d0)) then
                                                numer_weight = cart_x**moment_order(r)

                                                neg_numer(i, r) = neg_numer(i, r) + numer_weight*main_term
                                                neg_denom(i, r) = neg_denom(i, r) + main_term
                                            elseif ((dir == 2) .and. (cart_y >= 0d0)) then
                                                numer_weight = cart_y**moment_order(r)

                                                pos_numer(i, r) = pos_numer(i, r) + numer_weight*main_term
                                                pos_denom(i, r) = pos_denom(i, r) + main_term
                                            elseif ((dir == 2) .and. (cart_y < 0d0)) then
                                                numer_weight = cart_y**moment_order(r)

                                                neg_numer(i, r) = neg_numer(i, r) + numer_weight*main_term
                                                neg_denom(i, r) = neg_denom(i, r) + main_term
                                            end if
                                        else
                                            if (n > 0) then
                                                main_term = q_prim_vf(i + E_idx)%sf(j, k, l)* &
                                                            (1d0 - q_prim_vf(i + E_idx)%sf(j, k, l))* &
                                                            dx(j)*dy(k)
                                                if (p > 0) then
                                                    main_term = main_term*dz(l)
                                                end if
                                            end if
                                            if ((dir == 1) .and. (y_cc(k) >= 0d0)) then
                                                numer_weight = y_cc(k)**moment_order(r)

                                                pos_numer(i, r) = pos_numer(i, r) + numer_weight*main_term
                                                pos_denom(i, r) = pos_denom(i, r) + main_term
                                            elseif ((dir == 1) .and. (y_cc(k) < 0d0)) then
                                                numer_weight = y_cc(k)**moment_order(r)

                                                neg_numer(i, r) = neg_numer(i, r) + numer_weight*main_term
                                                neg_denom(i, r) = neg_denom(i, r) + main_term
                                            elseif ((dir == 2) .and. (z_cc(l) >= 0d0)) then
                                                numer_weight = z_cc(l)**moment_order(r)

                                                pos_numer(i, r) = pos_numer(i, r) + numer_weight*main_term
                                                pos_denom(i, r) = pos_denom(i, r) + main_term
                                            elseif ((dir == 2) .and. (z_cc(l) < 0d0)) then
                                                numer_weight = z_cc(l)**moment_order(r)

                                                neg_numer(i, r) = neg_numer(i, r) + numer_weight*main_term
                                                neg_denom(i, r) = neg_denom(i, r) + main_term
                                            end if
                                        end if
                                    end if
                                end do
                            end do
                        end do
                        ! Sum all components across all procs using MPI_ALLREDUCE
                        if (num_procs > 1) then
                            tmp = pos_numer(i, r)
                            call s_mpi_allreduce_sum(tmp, pos_numer(i, r))
                            tmp = neg_numer(i, r)
                            call s_mpi_allreduce_sum(tmp, neg_numer(i, r))
                            tmp = pos_denom(i, r)
                            call s_mpi_allreduce_sum(tmp, pos_denom(i, r))
                            tmp = neg_denom(i, r)
                            call s_mpi_allreduce_sum(tmp, neg_denom(i, r))
                        end if
                        ! Compute quotients and sum to get total moment
                        moments(i, dir, r) = (pos_numer(i, r)/pos_denom(i, r))**(1d0/moment_order(r)) + &
                                             (neg_numer(i, r)/neg_denom(i, r))**(1d0/moment_order(r))
                    else
                        moments(i, dir, r) = dflt_real
                    end if
                end do
            end if
        end do

    end subroutine s_non_symmetric_moments ! ------------------------------------------

    !> Compute symmetric moments
        !! @param q_prim_vf Primitive variables
        !! @param moments Higher moments(2 lateral directions, 5 moment orders)
        !! @param dir Current lateral direction
    subroutine s_symmetric_moments(q_prim_vf, moments, dir) ! ---------------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_fluids, 2, 5), intent(INOUT) :: moments
        integer, intent(IN) :: dir

        real(kind(0d0)), dimension(num_fluids, 5) :: numer, denom !<
            !! Numerator and denominator place holders for computation

        real(kind(0d0)) :: numer_weight     !< Numerator weight
        real(kind(0d0)) :: main_term        !< Constant term in both numerator and denominator
        real(kind(0d0)) :: dV               !< Discrete cell volume
        real(kind(0d0)) :: cart_x, cart_y   !< Cartesian x- and y-locations
        real(kind(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce

        integer :: i, j, k, l, r !< Generic loop iterators

        do i = 1, num_fluids
            if (com_wrt(i)) then
                numer(i, :) = 0d0
                denom(i, :) = 0d0
                do r = 1, 5
                    if (moment_order(r) /= dflt_int) then
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    if (q_prim_vf(i)%sf(j, k, l) > 5d-1) then
                                        if (grid_geometry == 3) then
                                            dV = (2d0*y_cb(k - 1)*dy(k) + dy(k)**2d0)/2d0*dx(j)*dz(l)
                                            cart_x = y_cc(k)*cos(z_cc(l))
                                            cart_y = y_cc(k)*sin(z_cc(l))
                                            main_term = q_prim_vf(i + E_idx)%sf(j, k, l)* &
                                                        (1d0 - q_prim_vf(i + E_idx)%sf(j, k, l))*dV
                                            if (dir == 1) then
                                                numer_weight = cart_x**moment_order(r)

                                                numer(i, r) = numer(i, r) + numer_weight*main_term
                                                denom(i, r) = denom(i, r) + main_term
                                            elseif (dir == 2) then
                                                numer_weight = cart_y**moment_order(r)

                                                numer(i, r) = numer(i, r) + numer_weight*main_term
                                                denom(i, r) = denom(i, r) + main_term
                                            end if
                                        else
                                            if (n > 0) then
                                                main_term = q_prim_vf(i + E_idx)%sf(j, k, l)* &
                                                            (1d0 - q_prim_vf(i + E_idx)%sf(j, k, l))* &
                                                            dx(j)*dy(k)
                                                if (p > 0) then
                                                    main_term = main_term*dz(l)
                                                end if
                                            end if
                                            if (dir == 1) then
                                                numer_weight = y_cc(k)**moment_order(r)

                                                numer(i, r) = numer(i, r) + numer_weight*main_term
                                                denom(i, r) = denom(i, r) + main_term
                                            elseif (dir == 2) then
                                                numer_weight = z_cc(l)**moment_order(r)

                                                numer(i, r) = numer(i, r) + numer_weight*main_term
                                                denom(i, r) = denom(i, r) + main_term
                                            end if
                                        end if
                                    end if
                                end do
                            end do
                        end do
                        ! Sum all components across all procs using MPI_ALLREDUCE
                        if (num_procs > 1) then
                            tmp = numer(i, r)
                            call s_mpi_allreduce_sum(tmp, numer(i, r))
                            tmp = denom(i, r)
                            call s_mpi_allreduce_sum(tmp, denom(i, r))
                        end if
                        ! Compute quotients and sum to get total moment
                        moments(i, dir, r) = (numer(i, r)/denom(i, r))**(1d0/moment_order(r))
                    else
                        moments(i, dir, r) = dflt_real
                    end if
                end do
            end if
        end do

    end subroutine s_symmetric_moments ! ------------------------------------------

    !>  This subroutine is used together with the volume fraction model
        !!      and when called upon, it computes the min and max bounds  of the
        !!      fluid in each direction in the domain.
        !!  @param q_prim_vf Primitive variables
        !!  @param bounds Variables storing the min and max bounds  of the fluids
    subroutine s_derive_fluid_bounds(q_prim_vf, bounds) ! -----------------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_fluids, 5, 6), intent(INOUT) :: bounds

        real(kind(0d0)) :: cart_x, cart_y, cart_z !< Cartesian x,y,z-locations
        real(kind(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce

        integer :: i, j, k, l, r !< Generic loop iterators

        if (n == 0) then ! 1D simulation
            do i = 1, num_fluids !Loop over individual fluids
                if (cb_wrt(i)) then
                    bounds(i, :, 1) = -1d0*dflt_real ! 1d6
                    bounds(i, :, 2) = dflt_real ! -1d6
                    do r = 1, 5
                        if (threshold_mf(r) /= dflt_real) then
                            do l = 0, p !Loop over grid
                                do k = 0, n
                                    do j = 0, m
                                        if ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                            .and. (x_cb(j - 1) <= bounds(i, r, 1))) then
                                            bounds(i, r, 1) = x_cb(j - 1)
                                        elseif ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                .and. (x_cb(j) >= bounds(i, r, 2))) then
                                            bounds(i, r, 2) = x_cb(j)
                                        end if
                                    end do
                                end do
                            end do

                            if (num_procs > 1) then
                                tmp = bounds(i, r, 1)
                                call s_mpi_allreduce_min(tmp, bounds(i, r, 1))
                                tmp = bounds(i, r, 2)
                                call s_mpi_allreduce_max(tmp, bounds(i, r, 2))
                            end if
                        else
                            bounds(i, r, 1) = dflt_real
                            bounds(i, r, 2) = dflt_real
                        end if
                    end do
                end if
            end do
        elseif (p == 0) then ! 2D simulation
            do i = 1, num_fluids !Loop over individual fluids
                if (cb_wrt(i)) then
                    bounds(i, :, 1) = -1d0*dflt_real ! 1d6
                    bounds(i, :, 2) = dflt_real ! -1d6
                    bounds(i, :, 3) = -1d0*dflt_real ! 1d6
                    bounds(i, :, 4) = dflt_real ! -1d6
                    do r = 1, 5
                        if (threshold_mf(r) /= dflt_real) then
                            do l = 0, p ! Loop over grid
                                do k = 0, n
                                    do j = 0, m
                                        if ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                            .and. (x_cb(j - 1) <= bounds(i, r, 1))) then
                                            bounds(i, r, 1) = x_cb(j - 1)
                                        elseif ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                .and. (x_cb(j) >= bounds(i, r, 2))) then
                                            bounds(i, r, 2) = x_cb(j)
                                        end if
                                        if ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                            .and. (y_cb(k - 1) <= bounds(i, r, 3))) then
                                            bounds(i, r, 3) = y_cb(k - 1)
                                        elseif ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                .and. (y_cb(k) >= bounds(i, r, 4))) then
                                            bounds(i, r, 4) = y_cb(k)
                                        end if
                                    end do
                                end do
                            end do

                            if (num_procs > 1) then
                                tmp = bounds(i, r, 1)
                                call s_mpi_allreduce_min(tmp, bounds(i, r, 1))
                                tmp = bounds(i, r, 2)
                                call s_mpi_allreduce_max(tmp, bounds(i, r, 2))
                                tmp = bounds(i, r, 3)
                                call s_mpi_allreduce_min(tmp, bounds(i, r, 3))
                                tmp = bounds(i, r, 4)
                                call s_mpi_allreduce_max(tmp, bounds(i, r, 4))
                            end if
                        else
                            bounds(i, r, 1) = dflt_real
                            bounds(i, r, 2) = dflt_real
                            bounds(i, r, 3) = dflt_real
                            bounds(i, r, 4) = dflt_real
                        end if
                    end do
                end if
            end do
        else ! 3D simulation
            do i = 1, num_fluids !Loop over individual fluids
                if (cb_wrt(i)) then
                    bounds(i, :, 1) = -1d0*dflt_real ! 1d6
                    bounds(i, :, 2) = dflt_real ! -1d6
                    bounds(i, :, 3) = -1d0*dflt_real ! 1d6
                    bounds(i, :, 4) = dflt_real ! -1d6
                    bounds(i, :, 5) = -1d0*dflt_real ! 1d6
                    bounds(i, :, 6) = dflt_real ! -1d6
                    do r = 1, 5
                        if (threshold_mf(r) /= dflt_real) then
                            do l = 0, p ! Loop over grid
                                do k = 0, n
                                    do j = 0, m
                                        if (grid_geometry == 3) then
                                            cart_x = y_cc(k)*cos(z_cc(l))
                                            cart_y = y_cc(k)*sin(z_cc(l))
                                            cart_z = x_cc(j)
                                            if ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                .and. (cart_x <= bounds(i, r, 1))) then
                                                bounds(i, r, 1) = cart_x
                                            elseif ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                    .and. (cart_x >= bounds(i, r, 2))) then
                                                bounds(i, r, 2) = cart_x
                                            end if
                                            if ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                .and. (cart_y <= bounds(i, r, 3))) then
                                                bounds(i, r, 3) = cart_y
                                            elseif ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                    .and. (cart_y >= bounds(i, r, 4))) then
                                                bounds(i, r, 4) = cart_y
                                            end if
                                            if ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                .and. (cart_z <= bounds(i, r, 5))) then
                                                bounds(i, r, 5) = cart_z
                                            elseif ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                    .and. (cart_z >= bounds(i, r, 6))) then
                                                bounds(i, r, 6) = cart_z
                                            end if
                                        else
                                            if ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                .and. (x_cb(j - 1) <= bounds(i, r, 1))) then
                                                bounds(i, r, 1) = x_cb(j - 1)
                                            elseif ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                    .and. (x_cb(j) >= bounds(i, r, 2))) then
                                                bounds(i, r, 2) = x_cb(j)
                                            end if
                                            if ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                .and. (y_cb(k - 1) <= bounds(i, r, 3))) then
                                                bounds(i, r, 3) = y_cb(k - 1)
                                            elseif ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                    .and. (y_cb(k) >= bounds(i, r, 4))) then
                                                bounds(i, r, 4) = y_cb(k)
                                            end if
                                            if ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                .and. (z_cb(l - 1) <= bounds(i, r, 5))) then
                                                bounds(i, r, 5) = z_cb(l - 1)
                                            elseif ((q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) &
                                                    .and. (z_cb(l) >= bounds(i, r, 6))) then
                                                bounds(i, r, 6) = z_cb(l)
                                            end if
                                        end if
                                    end do
                                end do
                            end do

                            if (num_procs > 1) then
                                tmp = bounds(i, r, 1)
                                call s_mpi_allreduce_min(tmp, bounds(i, r, 1))
                                tmp = bounds(i, r, 2)
                                call s_mpi_allreduce_max(tmp, bounds(i, r, 2))
                                tmp = bounds(i, r, 3)
                                call s_mpi_allreduce_min(tmp, bounds(i, r, 3))
                                tmp = bounds(i, r, 4)
                                call s_mpi_allreduce_max(tmp, bounds(i, r, 4))
                                tmp = bounds(i, r, 5)
                                call s_mpi_allreduce_min(tmp, bounds(i, r, 5))
                                tmp = bounds(i, r, 6)
                                call s_mpi_allreduce_max(tmp, bounds(i, r, 6))
                            end if
                        else
                            bounds(i, r, 1) = dflt_real
                            bounds(i, r, 2) = dflt_real
                            bounds(i, r, 3) = dflt_real
                            bounds(i, r, 4) = dflt_real
                            bounds(i, r, 5) = dflt_real
                            bounds(i, r, 6) = dflt_real
                        end if
                    end do
                end if
            end do
        end if

    end subroutine s_derive_fluid_bounds ! -------------------------------------------

    !>  This subroutine is used together with the volume fraction model
        !!      and when called upon, it computes the total mass of a fluid in the
        !!      entire domain for which the volume fraction is greater than a
        !!      threshold value. This gives the mass of the coherent body.
        !!  @param q_prim_vf Primitive variables
        !!  @param cb_mass Coherent body mass
    subroutine s_derive_coherent_body(q_prim_vf, cb_mass) ! --------------------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_fluids, 10), intent(INOUT) :: cb_mass

        real(kind(0d0)) :: dV  !< Discrete cell volume
        integer :: i, j, k, l, r   !< Generic loop iterators
        real(kind(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce

        do i = 1, num_fluids !Loop over individual fluids
            if (cb_wrt(i)) then
                cb_mass(i, :) = 0d0
                do r = 1, 5 ! Volume fraction threshold values
                    if (threshold_mf(r) /= dflt_real) then
                        do l = 0, p !Loop over grid
                            do k = 0, n
                                do j = 0, m
                                    if (q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) then
                                        if (n == 0) then
                                            dV = dx(j)
                                        elseif (p == 0) then
                                            dV = dx(j)*dy(k)
                                        else
                                            if (grid_geometry == 3) then
                                                dV = (2d0*y_cb(k - 1)*dy(k) + dy(k)**2d0)/2d0*dx(j)*dz(l)
                                            else
                                                dV = dx(j)*dy(k)*dz(l)
                                            end if
                                        end if
                                        cb_mass(i, r) = cb_mass(i, r) + q_prim_vf(i)%sf(j, k, l)*dV ! Mass
                                        cb_mass(i, r + 5) = cb_mass(i, r + 5) + dV ! Volume
                                    end if
                                end do
                            end do
                        end do

                        if (num_procs > 1) then
                            tmp = cb_mass(i, r)
                            call s_mpi_allreduce_sum(tmp, cb_mass(i, r))
                            tmp = cb_mass(i, r + 5)
                            call s_mpi_allreduce_sum(tmp, cb_mass(i, r + 5))
                        end if
                    else
                        cb_mass(i, r) = dflt_real
                        cb_mass(i, r + 5) = dflt_real
                    end if
                end do
            end if
        end do

        do i = 1, num_fluids
            if (cb_wrt(i)) then
                do r = 1, 5
                    if (threshold_mf(r) /= dflt_real) then
                        ! Check for reflective BC in x-direction
                        if (bc_x_glb%beg == -2) then
                            cb_mass(i, r) = cb_mass(i, r)*2d0
                            cb_mass(i, r + 5) = cb_mass(i, r + 5)*2d0
                        elseif (bc_x_glb%end == -2) then
                            cb_mass(i, r) = cb_mass(i, r)*2d0
                            cb_mass(i, r + 5) = cb_mass(i, r + 5)*2d0
                        end if
                        if (n > 0) then
                            ! Check for reflective BC in y-direction
                            if (bc_y_glb%beg == -2) then
                                cb_mass(i, r) = cb_mass(i, r)*2d0
                                cb_mass(i, r + 5) = cb_mass(i, r + 5)*2d0
                            elseif (bc_y_glb%end == -2) then
                                cb_mass(i, r) = cb_mass(i, r)*2d0
                                cb_mass(i, r + 5) = cb_mass(i, r + 5)*2d0
                            end if
                            if (p > 0) then
                                ! Check for reflective BC in z-direction
                                if (bc_z_glb%beg == -2) then
                                    cb_mass(i, r) = cb_mass(i, r)*2d0
                                    cb_mass(i, r + 5) = cb_mass(i, r + 5)*2d0
                                elseif (bc_z_glb%end == -2) then
                                    cb_mass(i, r) = cb_mass(i, r)*2d0
                                    cb_mass(i, r + 5) = cb_mass(i, r + 5)*2d0
                                end if

                            end if
                        end if
                    end if
                end do
            end if
        end do

    end subroutine s_derive_coherent_body ! ------------------------------

    !>  This subroutine is used together with the volume fraction model
        !!      and when called upon, it computes the centerline length of the
        !!      fluid.
        !!  @param q_prim_vf Primitive variables
        !!  @param cntrline Variables storing the centerline length of the fluids
    subroutine s_derive_centerline(q_prim_vf, cntrline) ! --------------------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_fluids, 5), intent(INOUT) :: cntrline

        real(kind(0d0)), dimension(5) :: cntrmin, cntrmax !< Placeholders
        real(kind(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce

        integer :: i, j, k, l, r !< Generic loop iterators

        if (n == 0) then ! 1D simulation
            do i = 1, num_fluids
                if (cb_wrt(i)) then
                    do r = 1, 5
                        if (threshold_mf(r) /= dflt_real) then
                            cntrline(i, r) = 0d0
                        else
                            cntrline(i, r) = dflt_real
                        end if
                    end do
                end if
            end do
        elseif ((p == 0) .or. (grid_geometry == 3)) then ! 2D simulation
            do i = 1, num_fluids
                if (cb_wrt(i)) then
                    cntrmin(:) = -1d0*dflt_real
                    cntrmax(:) = dflt_real
                    do r = 1, 5
                        if (threshold_mf(r) /= dflt_real) then
                            do l = 0, p
                                do k = 0, n
                                    do j = 0, m
                                        if ((y_cb(k - 1) == 0d0) .and. &
                                            (q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) .and. &
                                            (x_cb(j - 1) <= cntrmin(r))) then
                                            cntrmin(r) = x_cb(j - 1)
                                        elseif ((y_cb(k - 1) == 0d0) .and. &
                                                (q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) .and. &
                                                (x_cb(j) >= cntrmax(r))) then
                                            cntrmax(r) = x_cb(j)
                                        end if
                                    end do
                                end do
                            end do

                            if (num_procs > 1) then
                                tmp = cntrmin(r)
                                call s_mpi_allreduce_min(tmp, cntrmin(r))
                                tmp = cntrmax(r)
                                call s_mpi_allreduce_max(tmp, cntrmax(r))
                            end if

                            cntrline(i, r) = cntrmax(r) - cntrmin(r)
                        else
                            cntrline(i, r) = dflt_real
                        end if
                    end do
                end if
            end do
        else ! 3D simulation
            do i = 1, num_fluids
                if (cb_wrt(i)) then
                    cntrmin(:) = -1d0*dflt_real
                    cntrmax(:) = dflt_real
                    do r = 1, 5
                        if (threshold_mf(r) /= dflt_real) then
                            do l = 0, p
                                do k = 0, n
                                    do j = 0, m
                                        if ((y_cb(k - 1) == 0d0) .and. &
                                            (z_cb(l - 1) == 0d0) .and. &
                                            (q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) .and. &
                                            (x_cb(j - 1) <= cntrmin(r))) then
                                            cntrmin(r) = x_cb(j - 1)
                                        elseif ((y_cb(k - 1) == 0d0) .and. &
                                                (z_cb(l - 1) == 0d0) .and. &
                                                (q_prim_vf(i + E_idx)%sf(j, k, l) >= threshold_mf(r)) .and. &
                                                (x_cb(j) >= cntrmax(r))) then
                                            cntrmax(r) = x_cb(j)
                                        end if
                                    end do
                                end do
                            end do

                            if (num_procs > 1) then
                                tmp = cntrmin(r)
                                call s_mpi_allreduce_min(tmp, cntrmin(r))
                                tmp = cntrmax(r)
                                call s_mpi_allreduce_max(tmp, cntrmax(r))
                            end if

                            cntrline(i, r) = cntrmax(r) - cntrmin(r)
                        else
                            cntrline(i, r) = dflt_real
                        end if
                    end do
                end if
            end do
        end if

    end subroutine s_derive_centerline ! ----------------------------------

    !> Deallocation procedures for the module
    subroutine s_finalize_derived_variables_module() ! -------------------

        ! Closing CoM and flow probe files
        if (proc_rank == 0) then
            if (any(com_wrt)) then
                call s_close_com_files()
            end if
            if (any(cb_wrt)) then
                call s_close_cb_files()
            end if
            if (probe_wrt) then
                call s_close_probe_files()
            end if
        end if

        ! Deallocating the variables that might have been used to bookkeep
        ! the finite-difference coefficients in the x-, y- and z-directions
        if (allocated(fd_coeff_x)) deallocate (fd_coeff_x)
        if (allocated(fd_coeff_y)) deallocate (fd_coeff_y)
        if (allocated(fd_coeff_z)) deallocate (fd_coeff_z)

    end subroutine s_finalize_derived_variables_module ! -----------------

end module m_derived_variables
