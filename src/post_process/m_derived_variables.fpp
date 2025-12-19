!>
!! @file m_derived_variables.f90
!! @brief Contains module m_derived_variables

!> @brief This module features subroutines that allow for the derivation of
!!      numerous flow variables from the conservative and primitive ones.
!!      Currently, the available derived variables include the unadvected
!!      volume fraction, specific heat ratio, liquid stiffness, speed of
!!      sound, vorticity and the numerical Schlieren function.

module m_derived_variables

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_variables_conversion

    implicit none

    private; public :: s_initialize_derived_variables_module, &
 s_derive_specific_heat_ratio, &
 s_derive_liquid_stiffness, &
 s_derive_sound_speed, &
 s_derive_flux_limiter, &
 s_derive_vorticity_component, &
 s_derive_qm, &
 s_derive_liutex, &
 s_derive_numerical_schlieren_function, &
 s_compute_speed_of_sound, &
 s_finalize_derived_variables_module

    real(wp), allocatable, dimension(:, :, :) :: gm_rho_sf !<
    !! Gradient magnitude (gm) of the density for each cell of the computational
    !! sub-domain. This variable is employed in the calculation of the numerical
    !! Schlieren function.

    !> @name Finite-difference (fd) coefficients in x-, y- and z-coordinate directions.
    !! Note that because sufficient boundary information is available for all the
    !! active coordinate directions, the centered family of the finite-difference
    !! schemes is used.
    !> @{
    real(wp), allocatable, dimension(:, :), public :: fd_coeff_x
    real(wp), allocatable, dimension(:, :), public :: fd_coeff_y
    real(wp), allocatable, dimension(:, :), public :: fd_coeff_z
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
    impure subroutine s_initialize_derived_variables_module

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
        if (omega_wrt(2) .or. omega_wrt(3) .or. schlieren_wrt .or. liutex_wrt) then
            allocate (fd_coeff_x(-fd_number:fd_number, &
                                 -offset_x%beg:m + offset_x%end))
        end if

        ! Allocating centered finite-difference coefficients in y-direction
        if (omega_wrt(1) .or. omega_wrt(3) .or. liutex_wrt &
            .or. &
            (n > 0 .and. schlieren_wrt)) then
            allocate (fd_coeff_y(-fd_number:fd_number, &
                                 -offset_y%beg:n + offset_y%end))
        end if

        ! Allocating centered finite-difference coefficients in z-direction
        if (omega_wrt(1) .or. omega_wrt(2) .or. liutex_wrt &
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

    end subroutine s_initialize_derived_variables_module

    !>  This subroutine receives as input the specific heat ratio
        !!      function, gamma_sf, and derives from it the specific heat
        !!      ratio. The latter is stored in the derived flow quantity
        !!      storage variable, q_sf.
        !!  @param q_sf Specific heat ratio
    subroutine s_derive_specific_heat_ratio(q_sf)

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        integer :: i, j, k !< Generic loop iterators

        ! Computing specific heat ratio from specific heat ratio function
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end
                    q_sf(i, j, k) = 1._wp + 1._wp/gamma_sf(i, j, k)
                end do
            end do
        end do

    end subroutine s_derive_specific_heat_ratio

    !>  This subroutine admits as inputs the specific heat ratio
        !!      function and the liquid stiffness function, gamma_sf and
        !!      pi_inf_sf, respectively. These are used to calculate the
        !!      values of the liquid stiffness, which are stored in the
        !!      derived flow quantity storage variable, q_sf.
        !!  @param q_sf Liquid stiffness
    subroutine s_derive_liquid_stiffness(q_sf)

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        integer :: i, j, k !< Generic loop iterators

        ! Calculating the values of the liquid stiffness from those of the
        ! specific heat ratio function and the liquid stiffness function
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end
                    q_sf(i, j, k) = pi_inf_sf(i, j, k)/(gamma_sf(i, j, k) + 1._wp)
                end do
            end do
        end do

    end subroutine s_derive_liquid_stiffness

    !> This subroutine admits as inputs the primitive variables,
        !!      the density, the specific heat ratio function and liquid
        !!      stiffness function. It then computes from those variables
        !!      the values of the speed of sound, which are stored in the
        !!      derived flow quantity storage variable, q_sf.
        !! @param q_prim_vf Primitive variables
        !! @param q_sf Speed of sound
    subroutine s_derive_sound_speed(q_prim_vf, q_sf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        integer :: i, j, k !< Generic loop iterators

        ! Fluid bulk modulus for alternate sound speed
        real(wp) :: blkmod1, blkmod2

        ! Computing speed of sound values from those of pressure, density,
        ! specific heat ratio function and the liquid stiffness function
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end

                    ! Compute mixture sound speed
                    if (alt_soundspeed .neqv. .true.) then
                        q_sf(i, j, k) = (((gamma_sf(i, j, k) + 1._wp)* &
                                          q_prim_vf(E_idx)%sf(i, j, k) + &
                                          pi_inf_sf(i, j, k))/(gamma_sf(i, j, k)* &
                                                               rho_sf(i, j, k)))
                    else
                        blkmod1 = ((gammas(1) + 1._wp)*q_prim_vf(E_idx)%sf(i, j, k) + &
                                   pi_infs(1))/gammas(1)
                        blkmod2 = ((gammas(2) + 1._wp)*q_prim_vf(E_idx)%sf(i, j, k) + &
                                   pi_infs(2))/gammas(2)
                        q_sf(i, j, k) = (1._wp/(rho_sf(i, j, k)*(q_prim_vf(adv_idx%beg)%sf(i, j, k)/blkmod1 + &
                                                                 (1._wp - q_prim_vf(adv_idx%beg)%sf(i, j, k))/blkmod2)))
                    end if

                    if (mixture_err .and. q_sf(i, j, k) < 0._wp) then
                        q_sf(i, j, k) = 1.e-16_wp
                    else
                        q_sf(i, j, k) = sqrt(q_sf(i, j, k))
                    end if
                end do
            end do
        end do

    end subroutine s_derive_sound_speed

    !>  This subroutine derives the flux_limiter at cell boundary
        !!      i+1/2. This is an approximation because the velocity used
        !!      to determine the upwind direction is the velocity at the
        !!      cell center i instead of the contact velocity at the cell
        !!      boundary from the Riemann solver.
        !!  @param i Component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Flux limiter
    subroutine s_derive_flux_limiter(i, q_prim_vf, q_sf)

        integer, intent(in) :: i

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp), dimension(-offset_x%beg:m + offset_x%end, &
                            -offset_y%beg:n + offset_y%end, &
                            -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        real(wp) :: top, bottom, slope !< Flux limiter calcs
        integer :: j, k, l !< Generic loop iterators

        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end
                    if (i == 1) then
                        if (q_prim_vf(cont_idx%end + i)%sf(j, k, l) >= 0._wp) then
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
                        if (q_prim_vf(cont_idx%end + i)%sf(j, k, l) >= 0._wp) then
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
                        if (q_prim_vf(cont_idx%end + i)%sf(j, k, l) >= 0._wp) then
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

                    if (abs(top) < 1.e-8_wp) top = 0._wp
                    if (abs(bottom) < 1.e-8_wp) bottom = 0._wp

                    if (f_approx_equal(top, bottom)) then
                        slope = 1._wp
                        !       ELSEIF((top == 0._wp .AND. bottom /= 0._wp) &
                        !               .OR.            &
                        !           (bottom == 0._wp .AND. top /= 0._wp)) THEN
                        !           slope = 0._wp
                    else
                        slope = (top*bottom)/(bottom**2._wp + 1.e-16_wp)
                    end if

                    ! Flux limiter function
                    if (flux_lim == 1) then ! MINMOD (MM)
                        q_sf(j, k, l) = max(0._wp, min(1._wp, slope))
                    elseif (flux_lim == 2) then ! MUSCL (MC)
                        q_sf(j, k, l) = max(0._wp, min(2._wp*slope, 5.e-1_wp*(1._wp + slope), 2._wp))
                    elseif (flux_lim == 3) then ! OSPRE (OP)
                        q_sf(j, k, l) = (15.e-1_wp*(slope**2._wp + slope))/(slope**2._wp + slope + 1._wp)
                    elseif (flux_lim == 4) then ! SUPERBEE (SB)
                        q_sf(j, k, l) = max(0._wp, min(1._wp, 2._wp*slope), min(slope, 2._wp))
                    elseif (flux_lim == 5) then ! SWEBY (SW) (beta = 1.5)
                        q_sf(j, k, l) = max(0._wp, min(15.e-1_wp*slope, 1._wp), min(slope, 15.e-1_wp))
                    elseif (flux_lim == 6) then ! VAN ALBADA (VA)
                        q_sf(j, k, l) = (slope**2._wp + slope)/(slope**2._wp + 1._wp)
                    elseif (flux_lim == 7) then ! VAN LEER (VL)
                        q_sf(j, k, l) = (abs(slope) + slope)/(1._wp + abs(slope))
                    end if
                end do
            end do
        end do
    end subroutine s_derive_flux_limiter

    !>  Computes the solution to the linear system Ax=b w/ sol = x
        !!  @param A Input matrix
        !!  @param b right-hane-side
        !!  @param sol Solution
        !!  @param ndim Problem size
    subroutine s_solve_linear_system(A, b, sol, ndim)

        integer, intent(in) :: ndim
        real(wp), dimension(ndim, ndim), intent(inout) :: A
        real(wp), dimension(ndim), intent(inout) :: b
        real(wp), dimension(ndim), intent(out) :: sol

        !EXTERNAL DGESV

        integer :: i, j, k

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

    end subroutine s_solve_linear_system

    !>  This subroutine receives as inputs the indicator of the
        !!      component of the vorticity that should be outputted and
        !!      the primitive variables. From those inputs, it proceeds
        !!      to calculate values of the desired vorticity component,
        !!      which are subsequently stored in derived flow quantity
        !!      storage variable, q_sf.
        !!  @param i Vorticity component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Vorticity component
    subroutine s_derive_vorticity_component(i, q_prim_vf, q_sf)

        integer, intent(in) :: i

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        integer :: j, k, l, r !< Generic loop iterators

        ! Computing the vorticity component in the x-coordinate direction
        if (i == 1) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k, l) = 0._wp

                        do r = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                q_sf(j, k, l) = &
                                    q_sf(j, k, l) + 1._wp/y_cc(k)* &
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

                        q_sf(j, k, l) = 0._wp

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

                        q_sf(j, k, l) = 0._wp

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

    end subroutine s_derive_vorticity_component

    !> This subroutine gets as inputs the primitive variables. From those
        !!      inputs, it proceeds to calculate the value of the Q_M
        !!      function, which are subsequently stored in the derived flow
        !!      quantity storage variable, q_sf.
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Q_M
    subroutine s_derive_qm(q_prim_vf, q_sf)
        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        real(wp), &
            dimension(1:3, 1:3) :: q_jacobian_sf, S, S2, O, O2

        real(wp) :: trS, Q, IIS
        integer :: j, k, l, r, jj, kk !< Generic loop iterators

        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end

                    ! Get velocity gradient tensor
                    q_jacobian_sf(:, :) = 0._wp

                    do r = -fd_number, fd_number
                        do jj = 1, 3
                            ! d()/dx
                            q_jacobian_sf(jj, 1) = &
                                q_jacobian_sf(jj, 1) + &
                                fd_coeff_x(r, j)* &
                                q_prim_vf(mom_idx%beg + jj - 1)%sf(r + j, k, l)
                            ! d()/dy
                            q_jacobian_sf(jj, 2) = &
                                q_jacobian_sf(jj, 2) + &
                                fd_coeff_y(r, k)* &
                                q_prim_vf(mom_idx%beg + jj - 1)%sf(j, r + k, l)
                            ! d()/dz
                            q_jacobian_sf(jj, 3) = &
                                q_jacobian_sf(jj, 3) + &
                                fd_coeff_z(r, l)* &
                                q_prim_vf(mom_idx%beg + jj - 1)%sf(j, k, r + l)
                        end do
                    end do

                    ! Decompose J into asymmetric matrix, S, and a skew-symmetric matrix, O
                    do jj = 1, 3
                        do kk = 1, 3
                            S(jj, kk) = 0.5_wp* &
                                        (q_jacobian_sf(jj, kk) + q_jacobian_sf(kk, jj))
                            O(jj, kk) = 0.5_wp* &
                                        (q_jacobian_sf(jj, kk) - q_jacobian_sf(kk, jj))
                        end do
                    end do

                    ! Compute S2 = S*S'
                    do jj = 1, 3
                        do kk = 1, 3
                            O2(jj, kk) = O(jj, 1)*O(kk, 1) + &
                                         O(jj, 2)*O(kk, 2) + &
                                         O(jj, 3)*O(kk, 3)
                            S2(jj, kk) = S(jj, 1)*S(kk, 1) + &
                                         S(jj, 2)*S(kk, 2) + &
                                         S(jj, 3)*S(kk, 3)
                        end do
                    end do

                    ! Compute Q
                    Q = 0.5_wp*((O2(1, 1) + O2(2, 2) + O2(3, 3)) - &
                                (S2(1, 1) + S2(2, 2) + S2(3, 3)))
                    trS = S(1, 1) + S(2, 2) + S(3, 3)
                    IIS = 0.5_wp*((S(1, 1) + S(2, 2) + S(3, 3))**2 - &
                                  (S2(1, 1) + S2(2, 2) + S2(3, 3)))
                    q_sf(j, k, l) = Q + IIS

                end do
            end do
        end do

    end subroutine s_derive_qm

    !> This subroutine gets as inputs the primitive variables. From those
        !!      inputs, it proceeds to calculate the Liutex vector and its
        !!      magnitude based on Xu et al. (2019).
        !!  @param q_prim_vf Primitive variables
        !!  @param liutex_mag Liutex magnitude
        !!  @param liutex_axis Liutex axis
    impure subroutine s_derive_liutex(q_prim_vf, liutex_mag, liutex_axis)
        integer, parameter :: nm = 3
        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(out) :: liutex_mag !< Liutex magnitude

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end, nm), &
            intent(out) :: liutex_axis !< Liutex rigid rotation axis

        character, parameter :: ivl = 'N' !< compute left eigenvectors
        character, parameter :: ivr = 'V' !< compute right eigenvectors
        real(wp), dimension(nm, nm) :: vgt !< velocity gradient tensor
        real(wp), dimension(nm) :: lr, li !< real and imaginary parts of eigenvalues
        real(wp), dimension(nm, nm) :: vl, vr !< left and right eigenvectors
        integer, parameter :: lwork = 4*nm !< size of work array (4*nm recommended)
        real(wp), dimension(lwork) :: work !< work array
        integer :: info

        real(wp), dimension(nm) :: eigvec !< real eigenvector
        real(wp) :: eigvec_mag !< magnitude of real eigenvector
        real(wp) :: omega_proj !< projection of vorticity on real eigenvector
        real(wp) :: lci !< imaginary part of complex eigenvalue
        real(wp) :: alpha

        integer :: j, k, l, r, i !< Generic loop iterators
        integer :: idx

        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end

                    ! Get velocity gradient tensor (VGT)
                    vgt(:, :) = 0._wp

                    do r = -fd_number, fd_number
                        do i = 1, 3
                            ! d()/dx
                            vgt(i, 1) = &
                                vgt(i, 1) + &
                                fd_coeff_x(r, j)* &
                                q_prim_vf(mom_idx%beg + i - 1)%sf(r + j, k, l)
                            ! d()/dy
                            vgt(i, 2) = &
                                vgt(i, 2) + &
                                fd_coeff_y(r, k)* &
                                q_prim_vf(mom_idx%beg + i - 1)%sf(j, r + k, l)
                            ! d()/dz
                            vgt(i, 3) = &
                                vgt(i, 3) + &
                                fd_coeff_z(r, l)* &
                                q_prim_vf(mom_idx%beg + i - 1)%sf(j, k, r + l)
                        end do
                    end do

                    ! Call appropriate LAPACK routine based on precision
#ifdef MFC_SINGLE_PRECISION
                    call sgeev(ivl, ivr, nm, vgt, nm, lr, li, vl, nm, vr, nm, work, lwork, info)
#else
                    call dgeev(ivl, ivr, nm, vgt, nm, lr, li, vl, nm, vr, nm, work, lwork, info)
#endif

                    ! Find real eigenvector
                    idx = 1
                    do r = 2, 3
                        if (abs(li(r)) < abs(li(idx))) then
                            idx = r
                        end if
                    end do
                    eigvec = vr(:, idx)

                    ! Normalize real eigenvector if it is effectively non-zero
                    eigvec_mag = sqrt(eigvec(1)**2._wp &
                                      + eigvec(2)**2._wp &
                                      + eigvec(3)**2._wp)
                    if (eigvec_mag > sgm_eps) then
                        eigvec = eigvec/eigvec_mag
                    else
                        eigvec = 0._wp
                    end if

                    ! Compute vorticity projected on the eigenvector
                    omega_proj = (vgt(3, 2) - vgt(2, 3))*eigvec(1) &
                                 + (vgt(1, 3) - vgt(3, 1))*eigvec(2) &
                                 + (vgt(2, 1) - vgt(1, 2))*eigvec(3)

                    ! As eigenvector can have +/- signs, we can choose the sign
                    ! so that omega_proj is positive
                    if (omega_proj < 0._wp) then
                        eigvec = -eigvec
                        omega_proj = -omega_proj
                    end if

                    ! Find imaginary part of complex eigenvalue
                    lci = li(mod(idx, 3) + 1)

                    ! Compute Liutex magnitude
                    alpha = omega_proj**2._wp - 4._wp*lci**2._wp ! (2*alpha)^2
                    if (alpha > 0._wp) then
                        liutex_mag(j, k, l) = omega_proj - sqrt(alpha)
                    else
                        liutex_mag(j, k, l) = omega_proj
                    end if

                    ! Compute Liutex axis
                    liutex_axis(j, k, l, 1) = eigvec(1)
                    liutex_axis(j, k, l, 2) = eigvec(2)
                    liutex_axis(j, k, l, 3) = eigvec(3)

                end do
            end do
        end do

    end subroutine s_derive_liutex

    !>  This subroutine gets as inputs the conservative variables
        !!      and density. From those inputs, it proceeds to calculate
        !!      the values of the numerical Schlieren function, which are
        !!      subsequently stored in the derived flow quantity storage
        !!      variable, q_sf.
        !!  @param q_cons_vf Conservative variables
        !!  @param q_sf Numerical Schlieren function
    impure subroutine s_derive_numerical_schlieren_function(q_cons_vf, q_sf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_cons_vf

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        real(wp) :: drho_dx, drho_dy, drho_dz !<
            !! Spatial derivatives of the density in the x-, y- and z-directions

        real(wp), dimension(2) :: gm_rho_max !<
            !! Maximum value of the gradient magnitude (gm) of the density field
            !! in entire computational domain and not just the local sub-domain.
            !! The first position in the variable contains the maximum value and
            !! the second contains the rank of the processor on which it occurred.

        integer :: i, j, k, l !< Generic loop iterators

        ! Computing Gradient Magnitude of Density

        ! Contributions from the x- and y-coordinate directions
        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end

                    drho_dx = 0._wp
                    drho_dy = 0._wp

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

                        drho_dz = 0._wp

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

        ! Determining the local maximum of the gradient magnitude of density
        ! and bookkeeping the result, along with rank of the local processor
        gm_rho_max = (/maxval(gm_rho_sf), real(proc_rank, wp)/)

        ! Comparing the local maximum gradient magnitude of the density on
        ! this processor to the those computed on the remaining processors.
        ! This allows for the global maximum to be computed and the rank of
        ! the processor on which it has occurred to be recorded.
        if (num_procs > 1) call s_mpi_reduce_maxloc(gm_rho_max)

        ! Computing Numerical Schlieren Function

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

                        q_sf(j, k, l) = 0._wp

                        do i = 1, adv_idx%end - E_idx
                            q_sf(j, k, l) = &
                                q_sf(j, k, l) - schlieren_alpha(i)* &
                                q_cons_vf(i + E_idx)%sf(j, k, l)* &
                                gm_rho_sf(j, k, l)/gm_rho_max(1)
                        end do
                    end do
                end do
            end do
        end if

        ! Up until now, only the inside of the exponential of the numerical
        ! Schlieren function has been evaluated and stored. Then, to finish
        ! the computation, the exponential of the inside quantity is taken.
        q_sf = exp(q_sf)

    end subroutine s_derive_numerical_schlieren_function

    !>  Deallocation procedures for the module
    impure subroutine s_finalize_derived_variables_module

        ! Deallocating the variable containing the gradient magnitude of the
        ! density field provided that the numerical Schlieren function was
        ! was outputted during the post-process
        if (schlieren_wrt) deallocate (gm_rho_sf)

        ! Deallocating the variables that might have been used to bookkeep
        ! the finite-difference coefficients in the x-, y- and z-directions
        if (allocated(fd_coeff_x)) deallocate (fd_coeff_x)
        if (allocated(fd_coeff_y)) deallocate (fd_coeff_y)
        if (allocated(fd_coeff_z)) deallocate (fd_coeff_z)

    end subroutine s_finalize_derived_variables_module

end module m_derived_variables
