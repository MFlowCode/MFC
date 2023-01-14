!>
!! @file m_derived_variables.f90
!! @brief Contains module m_derived_variables

!> @brief This module features subroutines that allow for the derivation of
!!      numerous flow variables from the conservative and primitive ones.
!!      Currently, the available derived variables include the unadvected
!!      volume fraction, specific heat ratio, liquid stiffness, speed of
!!      sound, vorticity and the numerical Schlieren function.
module m_derived_variables

    ! Dependencies =============================================================
    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_derived_variables_module, &
 s_compute_finite_difference_coefficients, &
 s_derive_specific_heat_ratio, &
 s_derive_liquid_stiffness, &
 s_derive_sound_speed, &
 s_derive_flux_limiter, &
 s_derive_vorticity_component, &
 s_derive_numerical_schlieren_function, &
 s_finalize_derived_variables_module

    real(kind(0d0)), allocatable, dimension(:, :, :) :: gm_rho_sf !<
    !! Gradient magnitude (gm) of the density for each cell of the computational
    !! sub-domain. This variable is employed in the calculation of the numerical
    !! Schlieren function.

    !> @name Finite-difference (fd) coefficients in x-, y- and z-coordinate directions.
    !! Note that because sufficient boundary information is available for all the
    !! active coordinate directions, the centered family of the finite-difference
    !! schemes is used.
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :), public :: fd_coeff_x
    real(kind(0d0)), allocatable, dimension(:, :), public :: fd_coeff_y
    real(kind(0d0)), allocatable, dimension(:, :), public :: fd_coeff_z
    !> @}

    integer, private :: flg  !<
    !! Flagging (flg) variable used to annotate the dimensionality of the dataset
    !! that is undergoing the post-process. A flag value of 1 indicates that the
    !! dataset is 3D, while a flag value of 0 indicates that it is not. This flg
    !! variable is necessary to avoid cycling through the third dimension of the
    !! flow variable(s) when the simulation is not 3D and the size of the buffer
    !! is non-zero. Note that a similar procedure does not have to be applied to
    !! the second dimension since in 1D, the buffer size is always zero.

contains

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    subroutine s_initialize_derived_variables_module() ! ----------------------

        ! Allocating the gradient magnitude of the density variable provided
        ! that numerical Schlieren function is outputted during post-process
        if (schlieren_wrt) then
            allocate (gm_rho_sf(-offset_x%beg:m + offset_x%end, &
                                -offset_y%beg:n + offset_y%end, &
                                -offset_z%beg:p + offset_z%end))
        end if

        ! Allocating the variables which will store the coefficients of the
        ! centered family of finite-difference schemes. Note that sufficient
        ! space is allocated so that the coefficients up to any chosen order
        ! of accuracy may be bookkept. However, if higher than fourth-order
        ! accuracy coefficients are wanted, the formulae required to compute
        ! these coefficients will have to be implemented in the subroutine
        ! s_compute_finite_difference_coefficients.

        ! Allocating centered finite-difference coefficients in x-direction
        if (omega_wrt(2) .or. omega_wrt(3) .or. schlieren_wrt) then
            allocate (fd_coeff_x(-fd_number:fd_number, &
                                 -offset_x%beg:m + offset_x%end))
        end if

        ! Allocating centered finite-difference coefficients in y-direction
        if (omega_wrt(1) .or. omega_wrt(3) &
            .or. &
            (n > 0 .and. schlieren_wrt)) then
            allocate (fd_coeff_y(-fd_number:fd_number, &
                                 -offset_y%beg:n + offset_y%end))
        end if

        ! Allocating centered finite-difference coefficients in z-direction
        if (omega_wrt(1) .or. omega_wrt(2) &
            .or. &
            (p > 0 .and. schlieren_wrt)) then
            allocate (fd_coeff_z(-fd_number:fd_number, &
                                 -offset_z%beg:p + offset_z%end))
        end if

        ! Annotating the dimensionality of the dataset undergoing the post-
        ! process. A flag value of 1 indicates that the dataset is 3D, while
        ! a flag value of 0 indicates that it is not.
        if (p > 0) then
            flg = 1
        else
            flg = 0
        end if

    end subroutine s_initialize_derived_variables_module ! --------------------

    !> @name The purpose of this subroutine is to compute the finite-
        !!      difference coefficients for the centered schemes utilized
        !!      in computations of first order spatial derivatives in the
        !!      s-coordinate direction. The s-coordinate direction refers
        !!      to the x-, y- or z-coordinate direction, depending on the
        !!      subroutine's inputs. Note that coefficients of up to 4th
        !!      order accuracy are available.
        !!  @param q Number of cells in the s-coordinate direction
        !!  @param offset_s  Size of the ghost zone layer in the s-coordinate direction
        !!  @param s_cc Locations of the cell-centers in the s-coordinate direction
        !!  @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction
    subroutine s_compute_finite_difference_coefficients(q, offset_s, &
                                                        s_cc, fd_coeff_s)

        integer, intent(IN) :: q
        type(int_bounds_info), intent(IN) :: offset_s

        real(kind(0d0)), &
            dimension(-buff_size:q + buff_size), &
            intent(IN) :: s_cc

        real(kind(0d0)), &
            dimension(-fd_number:fd_number, -offset_s%beg:q + offset_s%end), &
            intent(INOUT) :: fd_coeff_s

        integer :: i !< Generic loop iterator

        ! Computing the 1st order finite-difference coefficients
        if (fd_order == 1) then
            do i = -offset_s%beg, q + offset_s%end
                fd_coeff_s(-1, i) = 0d0
                fd_coeff_s(0, i) = -1d0/(s_cc(i + 1) - s_cc(i))
                fd_coeff_s(1, i) = -fd_coeff_s(0, i)
            end do

            ! Computing the 2nd order finite-difference coefficients
        elseif (fd_order == 2) then
            do i = -offset_s%beg, q + offset_s%end
                fd_coeff_s(-1, i) = -1d0/(s_cc(i + 1) - s_cc(i - 1))
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
            end do

            ! Computing the 4th order finite-difference coefficients
        else
            do i = -offset_s%beg, q + offset_s%end
                fd_coeff_s(-2, i) = 1d0/(s_cc(i - 2) - 8d0*s_cc(i - 1) &
                                         - s_cc(i + 2) + 8d0*s_cc(i + 1))
                fd_coeff_s(-1, i) = -8d0*fd_coeff_s(-2, i)
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
                fd_coeff_s(2, i) = -fd_coeff_s(-2, i)
            end do

        end if

    end subroutine s_compute_finite_difference_coefficients ! --------------


    !>  This subroutine receives as input the specific heat ratio
        !!      function, gamma_sf, and derives from it the specific heat
        !!      ratio. The latter is stored in the derived flow quantity
        !!      storage variable, q_sf.
        !!  @param gamma_sf Specific heat ratio function
        !!  @param q_sf Specific heat ratio
    subroutine s_derive_specific_heat_ratio(gamma_sf, q_sf) ! --------------

        real(kind(0d0)), &
            dimension(-buff_size:m + buff_size, &
                      -buff_size:n + buff_size, &
                      -buff_size*flg:(p + buff_size)*flg), &
            intent(IN) :: gamma_sf

        real(kind(0d0)), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(INOUT) :: q_sf

        integer :: i, j, k !< Generic loop iterators

        ! Computing specific heat ratio from specific heat ratio function
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end
                    q_sf(i, j, k) = 1d0 + 1d0/gamma_sf(i, j, k)
                end do
            end do
        end do

    end subroutine s_derive_specific_heat_ratio ! --------------------------

    !>  This subroutine admits as inputs the specific heat ratio
        !!      function and the liquid stiffness function, gamma_sf and
        !!      pi_inf_sf, respectively. These are used to calculate the
        !!      values of the liquid stiffness, which are stored in the
        !!      derived flow quantity storage variable, q_sf.
        !!  @param gamma_sf Specific heat ratio
        !!  @param pi_inf_sf Liquid stiffness function
        !!  @param q_sf Liquid stiffness
    subroutine s_derive_liquid_stiffness(gamma_sf, pi_inf_sf, q_sf) ! ------

        real(kind(0d0)), &
            dimension(-buff_size:m + buff_size, &
                      -buff_size:n + buff_size, &
                      -buff_size*flg:(p + buff_size)*flg), &
            intent(IN) :: gamma_sf, pi_inf_sf

        real(kind(0d0)), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(INOUT) :: q_sf

        integer :: i, j, k !< Generic loop iterators

        ! Calculating the values of the liquid stiffness from those of the
        ! specific heat ratio function and the liquid stiffness function
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end
                    q_sf(i, j, k) = pi_inf_sf(i, j, k)/(gamma_sf(i, j, k) + 1d0)
                end do
            end do
        end do

    end subroutine s_derive_liquid_stiffness ! -----------------------------

    !> This subroutine admits as inputs the primitive variables,
        !!      the density, the specific heat ratio function and liquid
        !!      stiffness function. It then computes from those variables
        !!      the values of the speed of sound, which are stored in the
        !!      derived flow quantity storage variable, q_sf.
        !! @param q_prim_vf Primitive variables
        !! @param rho_sf Density
        !! @param gamma_sf Specific heat ratio function
        !! @param pi_inf_sf Liquid stiffness function
        !! @param q_sf Speed of sound
    subroutine s_derive_sound_speed(q_prim_vf, rho_sf, gamma_sf, & ! ------
                                    pi_inf_sf, q_sf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        real(kind(0d0)), &
            dimension(-buff_size:m + buff_size, &
                      -buff_size:n + buff_size, &
                      -buff_size*flg:(p + buff_size)*flg), &
            intent(IN) :: rho_sf, gamma_sf, pi_inf_sf

        real(kind(0d0)), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(INOUT) :: q_sf

        integer :: i, j, k !< Generic loop iterators

        ! Fluid bulk modulus for alternate sound speed
        real(kind(0d0)) :: blkmod1, blkmod2

        ! Computing speed of sound values from those of pressure, density,
        ! specific heat ratio function and the liquid stiffness function
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end

                    ! Compute mixture sound speed
                    if (alt_soundspeed .neqv. .true.) then
                        q_sf(i, j, k) = (((gamma_sf(i, j, k) + 1d0)* &
                                          q_prim_vf(E_idx)%sf(i, j, k) + &
                                          pi_inf_sf(i, j, k))/(gamma_sf(i, j, k)* &
                                                               rho_sf(i, j, k)))
                    else
                        blkmod1 = ((fluid_pp(1)%gamma + 1d0)*q_prim_vf(E_idx)%sf(i, j, k) + &
                                   fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                        blkmod2 = ((fluid_pp(2)%gamma + 1d0)*q_prim_vf(E_idx)%sf(i, j, k) + &
                                   fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                        q_sf(i, j, k) = (1d0/(rho_sf(i, j, k)*(q_prim_vf(adv_idx%beg)%sf(i, j, k)/blkmod1 + &
                                                               (1d0 - q_prim_vf(adv_idx%beg)%sf(i, j, k))/blkmod2)))
                    end if

                    if (mixture_err .and. q_sf(i, j, k) < 0d0) then
                        q_sf(i, j, k) = 1d-16
                    else
                        q_sf(i, j, k) = sqrt(q_sf(i, j, k))
                    end if
                end do
            end do
        end do

    end subroutine s_derive_sound_speed ! ----------------------------------

    !>  This subroutine derives the flux_limiter at cell boundary
        !!      i+1/2. This is an approximation because the velocity used
        !!      to determine the upwind direction is the velocity at the
        !!      cell center i instead of the contact velocity at the cell
        !!      boundary from the Riemann solver.
        !!  @param i Component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Flux limiter
    subroutine s_derive_flux_limiter(i, q_prim_vf, q_sf) ! -----------------

        integer, intent(IN) :: i

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        real(kind(0d0)), dimension(-offset_x%beg:m + offset_x%end, &
                                   -offset_y%beg:n + offset_y%end, &
                                   -offset_z%beg:p + offset_z%end), &
            intent(INOUT) :: q_sf

        real(kind(0d0)) :: top, bottom, slope !< Flux limiter calcs
        integer :: j, k, l !< Generic loop iterators

        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end
                    if (i == 1) then
                        if (q_prim_vf(cont_idx%end + i)%sf(j, k, l) >= 0d0) then
                            top = q_prim_vf(adv_idx%beg)%sf(j, k, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j - 1, k, l)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j + 1, k, l) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        else
                            top = q_prim_vf(adv_idx%beg)%sf(j + 2, k, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j + 1, k, l)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j + 1, k, l) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        end if
                    elseif (i == 2) then
                        if (q_prim_vf(cont_idx%end + i)%sf(j, k, l) >= 0d0) then
                            top = q_prim_vf(adv_idx%beg)%sf(j, k, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j, k - 1, l)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j, k + 1, l) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        else
                            top = q_prim_vf(adv_idx%beg)%sf(j, k + 2, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j, k + 1, l)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j, k + 1, l) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        end if
                    else
                        if (q_prim_vf(cont_idx%end + i)%sf(j, k, l) >= 0d0) then
                            top = q_prim_vf(adv_idx%beg)%sf(j, k, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j, k, l - 1)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j, k, l + 1) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        else
                            top = q_prim_vf(adv_idx%beg)%sf(j, k, l + 2) - &
                                  q_prim_vf(adv_idx%beg)%sf(j, k, l + 1)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j, k, l + 1) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        end if
                    end if

                    if (abs(top) < 1d-8) top = 0d0
                    if (abs(bottom) < 1d-8) bottom = 0d0

                    if (top == bottom) then
                        slope = 1d0
                        !       ELSEIF((top == 0d0 .AND. bottom /= 0d0) &
                        !               .OR.            &
                        !           (bottom == 0d0 .AND. top /= 0d0)) THEN
                        !           slope = 0d0
                    else
                        slope = (top*bottom)/(bottom**2d0 + 1d-16)
                    end if

                    ! Flux limiter function
                    if (flux_lim == 1) then ! MINMOD (MM)
                        q_sf(j, k, l) = max(0d0, min(1d0, slope))
                    elseif (flux_lim == 2) then ! MUSCL (MC)
                        q_sf(j, k, l) = max(0d0, min(2d0*slope, 5d-1*(1d0 + slope), 2d0))
                    elseif (flux_lim == 3) then ! OSPRE (OP)
                        q_sf(j, k, l) = (15d-1*(slope**2d0 + slope))/(slope**2d0 + slope + 1d0)
                    elseif (flux_lim == 4) then ! SUPERBEE (SB)
                        q_sf(j, k, l) = max(0d0, min(1d0, 2d0*slope), min(slope, 2d0))
                    elseif (flux_lim == 5) then ! SWEBY (SW) (beta = 1.5)
                        q_sf(j, k, l) = max(0d0, min(15d-1*slope, 1d0), min(slope, 15d-1))
                    elseif (flux_lim == 6) then ! VAN ALBADA (VA)
                        q_sf(j, k, l) = (slope**2d0 + slope)/(slope**2d0 + 1d0)
                    elseif (flux_lim == 7) then ! VAN LEER (VL)
                        q_sf(j, k, l) = (abs(slope) + slope)/(1d0 + abs(slope))
                    end if
                end do
            end do
        end do
    end subroutine s_derive_flux_limiter ! ---------------------------------


    !>  Computes the solution to the linear system Ax=b w/ sol = x
        !!  @param A Input matrix
        !!  @param b right-hand-side
        !!  @param sol Solution
        !!  @param ndim Problem size
    subroutine s_solve_linear_system(A, b, sol, ndim)

        integer, intent(IN) :: ndim
        real(kind(0d0)), dimension(ndim, ndim), intent(INOUT) :: A
        real(kind(0d0)), dimension(ndim), intent(INOUT) :: b
        real(kind(0d0)), dimension(ndim), intent(OUT) :: sol
        integer, dimension(ndim) :: ipiv

        integer :: nrhs, lda, ldb, info
        !EXTERNAL DGESV

        integer :: i, j, k

        ! Solve linear system using Intel MKL (Hooke)
!               nrhs = 1
!               lda = ndim
!               ldb = ndim
!
!               CALL DGESV(ndim, nrhs, A, lda, ipiv, b, ldb, info)
!
!               DO i = 1, ndim
!                   sol(i) = b(i)
!               END DO
!
!               IF (info /= 0) THEN
!                   PRINT '(A)', 'Trouble solving linear system'
!                   CALL s_mpi_abort()
!               END IF

        ! Solve linear system using own linear solver (Thomson/Darter/Comet/Stampede)
        ! Forward elimination
        do i = 1, ndim
            ! Pivoting
            j = i - 1 + maxloc(abs(A(i:ndim, i)), 1)
            sol = A(i, :)
            A(i, :) = A(j, :)
            A(j, :) = sol
            sol(1) = b(i)
            b(i) = b(j)
            b(j) = sol(1)
            ! Elimination
            b(i) = b(i)/A(i, i)
            A(i, :) = A(i, :)/A(i, i)
            do k = i + 1, ndim
                b(k) = b(k) - A(k, i)*b(i)
                A(k, :) = A(k, :) - A(k, i)*A(i, :)
            end do
        end do

        ! Backward substitution
        do i = ndim, 1, -1
            sol(i) = b(i)
            do k = i + 1, ndim
                sol(i) = sol(i) - A(i, k)*sol(k)
            end do
        end do

    end subroutine s_solve_linear_system ! -------------------------------------

    !>  This subroutine receives as inputs the indicator of the
        !!      component of the vorticity that should be outputted and
        !!      the primitive variables. From those inputs, it proceeds
        !!      to calculate values of the desired vorticity component,
        !!      which are subsequently stored in derived flow quantity
        !!      storage variable, q_sf.
        !!  @param i Vorticity component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Vorticity component
    subroutine s_derive_vorticity_component(i, q_prim_vf, q_sf) ! ----------

        integer, intent(IN) :: i

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        real(kind(0d0)), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(INOUT) :: q_sf

        integer :: j, k, l, r !< Generic loop iterators

        ! Computing the vorticity component in the x-coordinate direction
        if (i == 1) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k, l) = 0d0

                        do r = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                q_sf(j, k, l) = &
                                    q_sf(j, k, l) + 1d0/y_cc(k)* &
                                    (fd_coeff_y(r, k)*y_cc(r + k)* &
                                     q_prim_vf(mom_idx%end)%sf(j, r + k, l) &
                                     - fd_coeff_z(r, l)* &
                                     q_prim_vf(mom_idx%beg + 1)%sf(j, k, r + l))
                            else
                                q_sf(j, k, l) = &
                                    q_sf(j, k, l) + fd_coeff_y(r, k)* &
                                    q_prim_vf(mom_idx%end)%sf(j, r + k, l) &
                                    - fd_coeff_z(r, l)* &
                                    q_prim_vf(mom_idx%beg + 1)%sf(j, k, r + l)
                            end if
                        end do

                    end do
                end do
            end do

            ! Computing the vorticity component in the y-coordinate direction
        elseif (i == 2) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k, l) = 0d0

                        do r = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                q_sf(j, k, l) = &
                                    q_sf(j, k, l) + fd_coeff_z(r, l)/y_cc(k)* &
                                    q_prim_vf(mom_idx%beg)%sf(j, k, r + l) &
                                    - fd_coeff_x(r, j)* &
                                    q_prim_vf(mom_idx%end)%sf(r + j, k, l)
                            else
                                q_sf(j, k, l) = &
                                    q_sf(j, k, l) + fd_coeff_z(r, l)* &
                                    q_prim_vf(mom_idx%beg)%sf(j, k, r + l) &
                                    - fd_coeff_x(r, j)* &
                                    q_prim_vf(mom_idx%end)%sf(r + j, k, l)
                            end if
                        end do

                    end do
                end do
            end do

            ! Computing the vorticity component in the z-coordinate direction
        else
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k, l) = 0d0

                        do r = -fd_number, fd_number
                            q_sf(j, k, l) = &
                                q_sf(j, k, l) + fd_coeff_x(r, j)* &
                                q_prim_vf(mom_idx%beg + 1)%sf(r + j, k, l) &
                                - fd_coeff_y(r, k)* &
                                q_prim_vf(mom_idx%beg)%sf(j, r + k, l)
                        end do

                    end do
                end do
            end do
        end if

    end subroutine s_derive_vorticity_component ! --------------------------

    !>  This subroutine gets as inputs the conservative variables
        !!      and density. From those inputs, it proceeds to calculate
        !!      the values of the numerical Schlieren function, which are
        !!      subsequently stored in the derived flow quantity storage
        !!      variable, q_sf.
        !!  @param q_cons_vf Conservative variables
        !!  @param rho_sf Density
        !!  @param q_sf Numerical Schlieren function
    subroutine s_derive_numerical_schlieren_function(q_cons_vf, rho_sf, q_sf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        real(kind(0d0)), &
            dimension(-buff_size:m + buff_size, &
                      -buff_size:n + buff_size, &
                      -buff_size*flg:(p + buff_size)*flg), &
            intent(IN) :: rho_sf

        real(kind(0d0)), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(INOUT) :: q_sf

        real(kind(0d0)) :: drho_dx, drho_dy, drho_dz !<
            !! Spatial derivatives of the density in the x-, y- and z-directions

        real(kind(0d0)), dimension(2) :: gm_rho_max !<
            !! Maximum value of the gradient magnitude (gm) of the density field
            !! in entire computational domain and not just the local sub-domain.
            !! The first position in the variable contains the maximum value and
            !! the second contains the rank of the processor on which it occured.

        real(kind(0d0)) :: alpha_unadv !< Unadvected volume fraction

        integer :: i, j, k, l !< Generic loop iterators

        ! Computing Gradient Magnitude of Density ==========================

        ! Contributions from the x- and y-coordinate directions
        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end

                    drho_dx = 0d0
                    drho_dy = 0d0

                    do i = -fd_number, fd_number
                        drho_dx = drho_dx + fd_coeff_x(i, j)*rho_sf(i + j, k, l)
                        drho_dy = drho_dy + fd_coeff_y(i, k)*rho_sf(j, i + k, l)
                    end do

                    gm_rho_sf(j, k, l) = drho_dx*drho_dx + drho_dy*drho_dy

                end do
            end do
        end do

        ! Contribution from the z-coordinate direction
        if (p > 0) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        drho_dz = 0d0

                        do i = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                drho_dz = drho_dz + fd_coeff_z(i, l)/y_cc(k)* &
                                          rho_sf(j, k, i + l)
                            else
                                drho_dz = drho_dz + fd_coeff_z(i, l)* &
                                          rho_sf(j, k, i + l)
                            end if
                        end do

                        gm_rho_sf(j, k, l) = gm_rho_sf(j, k, l) &
                                             + drho_dz*drho_dz

                    end do
                end do
            end do
        end if

        ! Up until now, only the dot product of the gradient of the density
        ! field has been calculated and stored in the gradient magnitude of
        ! density variable. So now we proceed to take the square-root as to
        ! complete the desired calculation.
        gm_rho_sf = sqrt(gm_rho_sf)

        ! ==================================================================

        ! Determining the local maximum of the gradient magnitude of density
        ! and bookkeeping the result, along with rank of the local processor
        gm_rho_max = (/maxval(gm_rho_sf), real(proc_rank, kind(0d0))/)

        ! Comparing the local maximum gradient magnitude of the density on
        ! this processor to the those computed on the remaining processors.
        ! This allows for the global maximum to be computed and the rank of
        ! the processor on which it has occured to be recorded.
        if (num_procs > 1) call s_mpi_reduce_maxloc(gm_rho_max)

        ! Computing Numerical Schlieren Function ===========================

        ! The form of the numerical Schlieren function depends on the choice
        ! of the multicomponent flow model. For the gamma/pi_inf model, the
        ! exponential of the negative, normalized, gradient magnitude of the
        ! density is computed. For the volume fraction model, the amplitude
        ! of the exponential's inside is also modulated with respect to the
        ! identity of the fluid in which the function is evaluated. For more
        ! information, refer to Marquina and Mulet (2003).

        if (model_eqns == 1) then                    ! Gamma/pi_inf model
            q_sf = -gm_rho_sf/gm_rho_max(1)

        else                                        ! Volume fraction model
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k, l) = 0d0

                        do i = 1, adv_idx%end - E_idx
                            q_sf(j, k, l) = &
                                q_sf(j, k, l) - schlieren_alpha(i)* &
                                q_cons_vf(i + E_idx)%sf(j, k, l)* &
                                gm_rho_sf(j, k, l)/gm_rho_max(1)
                        end do

                        if (adv_alphan .neqv. .true.) then

                            alpha_unadv = 1d0

                            do i = 1, num_fluids - 1
                                alpha_unadv = alpha_unadv &
                                              - q_cons_vf(i + E_idx)%sf(j, k, l)
                            end do

                            q_sf(j, k, l) = q_sf(j, k, l) &
                                            - schlieren_alpha(num_fluids)* &
                                            alpha_unadv*gm_rho_sf(j, k, l)/ &
                                            gm_rho_max(1)

                        end if

                    end do
                end do
            end do
        end if

        ! Up until now, only the inside of the exponential of the numerical
        ! Schlieren function has been evaluated and stored. Then, to finish
        ! the computation, the exponential of the inside quantity is taken.
        q_sf = exp(q_sf)

        ! ==================================================================

    end subroutine s_derive_numerical_schlieren_function ! -----------------

    !>  Deallocation procedures for the module
    subroutine s_finalize_derived_variables_module() ! -------------------

        ! Deallocating the variable containing the gradient magnitude of the
        ! density field provided that the numerical Schlieren function was
        ! was outputted during the post-process
        if (schlieren_wrt) deallocate (gm_rho_sf)

        ! Deallocating the variables that might have been used to bookkeep
        ! the finite-difference coefficients in the x-, y- and z-directions
        if (allocated(fd_coeff_x)) deallocate (fd_coeff_x)
        if (allocated(fd_coeff_y)) deallocate (fd_coeff_y)
        if (allocated(fd_coeff_z)) deallocate (fd_coeff_z)

    end subroutine s_finalize_derived_variables_module ! -----------------

end module m_derived_variables
