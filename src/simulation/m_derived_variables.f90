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

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_data_output           !< Data output module

    use m_time_steppers         !< Time-stepping algorithms

    use m_helper
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

    ! @name Variables for computing acceleration
    !> @{
    real(kind(0d0)), public, allocatable, dimension(:, :, :) :: accel_mag
    real(kind(0d0)), public, allocatable, dimension(:, :, :) :: x_accel, y_accel, z_accel
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
            end if
            if (p > 0) then
                allocate (fd_coeff_z(-fd_number:fd_number, 0:p))
            end if

            allocate (accel_mag(0:m, 0:n, 0:p))
            allocate (x_accel(0:m, 0:n, 0:p))
            if (n > 0) then
                allocate (y_accel(0:m, 0:n, 0:p))
                if (p > 0) then
                    allocate (z_accel(0:m, 0:n, 0:p))
                end if
            end if
        end if

    end subroutine s_initialize_derived_variables_module ! --------------------

    !> Allocate and open derived variables. Computing FD coefficients.
    subroutine s_initialize_derived_variables() ! -----------------------------

        if (probe_wrt) then
            ! Opening and writing header of flow probe files
            if (proc_rank == 0) then
                call s_open_probe_files()
            end if

            ! Computing centered finite difference coefficients
            call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, &
                                                          fd_number, fd_order)

            if (n > 0) then
                call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, &
                                                              fd_number, fd_order)
            end if

            if (p > 0) then
                call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z, buff_size, &
                                                              fd_number, fd_order)
            end if

        end if

    end subroutine s_initialize_derived_variables ! -----------------------------

    !> Writes coherent body information, communication files, and probes.
        !!  @param t_step Current time-step
    subroutine s_compute_derived_variables(t_step) ! -----------------------

        integer, intent(IN) :: t_step

        integer :: i, j, k !< Generic loop iterators

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
            end if

            if (p > 0) then
                call s_derive_acceleration_component(3, q_prim_ts(0)%vf, &
                                                     q_prim_ts(1)%vf, &
                                                     q_prim_ts(2)%vf, &
                                                     q_prim_ts(3)%vf, &
                                                     z_accel)
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

    end subroutine s_compute_derived_variables ! ---------------------------

    !> This subroutine receives as inputs the indicator of the
        !!      component of the acceleration that should be outputted and
        !!      the primitive variables. From those inputs, it proceeds
        !!      to calculate values of the desired acceleration component,
        !!      which are subsequently stored in derived flow quantity
        !!      storage variable, q_sf.
        !!  @param i Acceleration component indicator
        !!  @param q_prim_vf0 Primitive variables
        !!  @param q_prim_vf1 Primitive variables
        !!  @param q_prim_vf2 Primitive variables
        !!  @param q_prim_vf3 Primitive variables
        !!  @param q_sf Acceleration component
    subroutine s_derive_acceleration_component(i, q_prim_vf0, q_prim_vf1, &
                                               q_prim_vf2, q_prim_vf3, q_sf) ! ----------
        !DIR$ INLINEALWAYS s_derive_acceleration_component
        integer, intent(IN) :: i

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf0
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

                        q_sf(j, k, l) = (11d0*q_prim_vf0(mom_idx%beg)%sf(j, k, l) &
                                         - 18d0*q_prim_vf1(mom_idx%beg)%sf(j, k, l) &
                                         + 9d0*q_prim_vf2(mom_idx%beg)%sf(j, k, l) &
                                         - 2d0*q_prim_vf3(mom_idx%beg)%sf(j, k, l))/(6d0*dt)

                        do r = -fd_number, fd_number
                            if (n == 0) then ! 1D simulation
                                print *, q_sf(j, k, l), q_prim_vf0(mom_idx%beg)%sf(j, k, l), fd_coeff_x(r, j), q_prim_vf0(mom_idx%beg)%sf(r + j, k, l)
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf0(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf0(mom_idx%beg)%sf(r + j, k, l)
                            elseif (p == 0) then ! 2D simulation
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf0(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf0(mom_idx%beg)%sf(r + j, k, l) &
                                                + q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                q_prim_vf0(mom_idx%beg)%sf(j, r + k, l)
                            else ! 3D simulation
                                if (grid_geometry == 3) then
                                    q_sf(j, k, l) = q_sf(j, k, l) &
                                                    + q_prim_vf0(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                    q_prim_vf0(mom_idx%beg)%sf(r + j, k, l) &
                                                    + q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                    q_prim_vf0(mom_idx%beg)%sf(j, r + k, l) &
                                                    + q_prim_vf0(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                    q_prim_vf0(mom_idx%beg)%sf(j, k, r + l)/y_cc(k)
                                else
                                    q_sf(j, k, l) = q_sf(j, k, l) &
                                                    + q_prim_vf0(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                    q_prim_vf0(mom_idx%beg)%sf(r + j, k, l) &
                                                    + q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                    q_prim_vf0(mom_idx%beg)%sf(j, r + k, l) &
                                                    + q_prim_vf0(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                    q_prim_vf0(mom_idx%beg)%sf(j, k, r + l)
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

                        q_sf(j, k, l) = (11d0*q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l) &
                                         - 18d0*q_prim_vf1(mom_idx%beg + 1)%sf(j, k, l) &
                                         + 9d0*q_prim_vf2(mom_idx%beg + 1)%sf(j, k, l) &
                                         - 2d0*q_prim_vf3(mom_idx%beg + 1)%sf(j, k, l))/(6d0*dt)

                        do r = -fd_number, fd_number
                            if (p == 0) then ! 2D simulation
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf0(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf0(mom_idx%beg + 1)%sf(r + j, k, l) &
                                                + q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                q_prim_vf0(mom_idx%beg + 1)%sf(j, r + k, l)
                            else ! 3D simulation
                                if (grid_geometry == 3) then
                                    q_sf(j, k, l) = q_sf(j, k, l) &
                                                    + q_prim_vf0(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                    q_prim_vf0(mom_idx%beg + 1)%sf(r + j, k, l) &
                                                    + q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                    q_prim_vf0(mom_idx%beg + 1)%sf(j, r + k, l) &
                                                    + q_prim_vf0(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                    q_prim_vf0(mom_idx%beg + 1)%sf(j, k, r + l)/y_cc(k) &
                                                    - (q_prim_vf0(mom_idx%end)%sf(j, k, l)**2d0)/y_cc(k)
                                else
                                    q_sf(j, k, l) = q_sf(j, k, l) &
                                                    + q_prim_vf0(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                    q_prim_vf0(mom_idx%beg + 1)%sf(r + j, k, l) &
                                                    + q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                    q_prim_vf0(mom_idx%beg + 1)%sf(j, r + k, l) &
                                                    + q_prim_vf0(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                    q_prim_vf0(mom_idx%beg + 1)%sf(j, k, r + l)
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
                        q_sf(j, k, l) = (11d0*q_prim_vf0(mom_idx%end)%sf(j, k, l) &
                                         - 18d0*q_prim_vf1(mom_idx%end)%sf(j, k, l) &
                                         + 9d0*q_prim_vf2(mom_idx%end)%sf(j, k, l) &
                                         - 2d0*q_prim_vf3(mom_idx%end)%sf(j, k, l))/(6d0*dt)

                        do r = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf0(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf0(mom_idx%end)%sf(r + j, k, l) &
                                                + q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                q_prim_vf0(mom_idx%end)%sf(j, r + k, l) &
                                                + q_prim_vf0(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                q_prim_vf0(mom_idx%end)%sf(j, k, r + l)/y_cc(k) &
                                                + (q_prim_vf0(mom_idx%end)%sf(j, k, l)* &
                                                   q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l))/y_cc(k)
                            else
                                q_sf(j, k, l) = q_sf(j, k, l) &
                                                + q_prim_vf0(mom_idx%beg)%sf(j, k, l)*fd_coeff_x(r, j)* &
                                                q_prim_vf0(mom_idx%end)%sf(r + j, k, l) &
                                                + q_prim_vf0(mom_idx%beg + 1)%sf(j, k, l)*fd_coeff_y(r, k)* &
                                                q_prim_vf0(mom_idx%end)%sf(j, r + k, l) &
                                                + q_prim_vf0(mom_idx%end)%sf(j, k, l)*fd_coeff_z(r, l)* &
                                                q_prim_vf0(mom_idx%end)%sf(j, k, r + l)
                            end if
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_derive_acceleration_component ! --------------------------

    !> Deallocation procedures for the module
    subroutine s_finalize_derived_variables_module() ! -------------------

        ! Closing CoM and flow probe files
        if (proc_rank == 0) then
            if (probe_wrt) then
                call s_close_probe_files()
            end if
        end if

        if (probe_wrt) then
            deallocate (accel_mag, x_accel)
            if (n > 0) then
                deallocate (y_accel)
                if (p > 0) then
                    deallocate (z_accel)
                end if
            end if
        end if

        ! Deallocating the variables that might have been used to bookkeep
        ! the finite-difference coefficients in the x-, y- and z-directions
        if (allocated(fd_coeff_x)) deallocate (fd_coeff_x)
        if (allocated(fd_coeff_y)) deallocate (fd_coeff_y)
        if (allocated(fd_coeff_z)) deallocate (fd_coeff_z)

    end subroutine s_finalize_derived_variables_module ! -----------------

end module m_derived_variables
