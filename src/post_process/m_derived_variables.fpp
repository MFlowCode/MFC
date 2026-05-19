!>
!! @file
!! @brief Contains module m_derived_variables

!> @brief Computes derived flow quantities (sound speed, vorticity, Schlieren, etc.) from conservative and primitive variables

module m_derived_variables

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic
    use m_variables_conversion

    implicit none

    private; public :: s_initialize_derived_variables_module, s_derive_specific_heat_ratio, s_derive_liquid_stiffness, &
        & s_derive_sound_speed, s_derive_flux_limiter, s_derive_vorticity_component, s_derive_qm, s_derive_liutex, &
        & s_derive_numerical_schlieren_function, s_compute_speed_of_sound, s_finalize_derived_variables_module

    real(wp), allocatable, dimension(:,:,:) :: gm_rho_sf  !< Density gradient magnitude for numerical Schlieren
    !> @name Finite-difference (fd) coefficients in x-, y- and z-coordinate directions. Note that because sufficient boundary
    !! information is available for all the active coordinate directions, the centered family of the finite-difference schemes is
    !! used.
    !> @{
    real(wp), allocatable, dimension(:,:), public :: fd_coeff_x
    real(wp), allocatable, dimension(:,:), public :: fd_coeff_y
    real(wp), allocatable, dimension(:,:), public :: fd_coeff_z
    !> @}

contains

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_derived_variables_module

        ! Allocate density gradient magnitude if Schlieren output requested
        if (schlieren_wrt) then
            allocate (gm_rho_sf(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end))
        end if

        ! Allocate FD coefficients (up to 4th order; higher orders need extension)

        if (omega_wrt(2) .or. omega_wrt(3) .or. schlieren_wrt .or. liutex_wrt) then
            allocate (fd_coeff_x(-fd_number:fd_number,-offset_x%beg:m + offset_x%end))
        end if

        if (omega_wrt(1) .or. omega_wrt(3) .or. liutex_wrt .or. (n > 0 .and. schlieren_wrt)) then
            allocate (fd_coeff_y(-fd_number:fd_number,-offset_y%beg:n + offset_y%end))
        end if

        if (omega_wrt(1) .or. omega_wrt(2) .or. liutex_wrt .or. (p > 0 .and. schlieren_wrt)) then
            allocate (fd_coeff_z(-fd_number:fd_number,-offset_z%beg:p + offset_z%end))
        end if

    end subroutine s_initialize_derived_variables_module

    !> Derive the specific heat ratio from the specific heat ratio function gamma_sf. The latter is stored in the derived flow
    !! quantity storage variable, q_sf.
    subroutine s_derive_specific_heat_ratio(q_sf)

        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end), &
             & intent(inout) :: q_sf

        integer :: i, j, k
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end
                    q_sf(i, j, k) = 1._wp + 1._wp/gamma_sf(i, j, k)
                end do
            end do
        end do

    end subroutine s_derive_specific_heat_ratio

    !> Compute the liquid stiffness from the specific heat ratio function gamma_sf and the liquid stiffness function pi_inf_sf,
    !! respectively. These are used to calculate the values of the liquid stiffness, which are stored in the derived flow quantity
    !! storage variable, q_sf.
    subroutine s_derive_liquid_stiffness(q_sf)

        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end), &
             & intent(inout) :: q_sf

        integer :: i, j, k
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end
                    q_sf(i, j, k) = pi_inf_sf(i, j, k)/(gamma_sf(i, j, k) + 1._wp)
                end do
            end do
        end do

    end subroutine s_derive_liquid_stiffness

    !> Compute the speed of sound from the primitive variables, density, specific heat ratio function, and liquid stiffness
    !! function. It then computes from those variables the values of the speed of sound, which are stored in the derived flow
    !! quantity storage variable, q_sf.
    subroutine s_derive_sound_speed(q_prim_vf, q_sf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end), &
             & intent(inout) :: q_sf

        integer  :: i, j, k
        real(wp) :: blkmod1, blkmod2

        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end
                    if (alt_soundspeed .neqv. .true.) then
                        q_sf(i, j, k) = (((gamma_sf(i, j, k) + 1._wp)*q_prim_vf(eqn_idx%E)%sf(i, j, k) + pi_inf_sf(i, j, &
                             & k))/(gamma_sf(i, j, k)*rho_sf(i, j, k)))
                    else
                        blkmod1 = ((gammas(1) + 1._wp)*q_prim_vf(eqn_idx%E)%sf(i, j, k) + pi_infs(1))/gammas(1)
                        blkmod2 = ((gammas(2) + 1._wp)*q_prim_vf(eqn_idx%E)%sf(i, j, k) + pi_infs(2))/gammas(2)
                        q_sf(i, j, k) = (1._wp/(rho_sf(i, j, k)*(q_prim_vf(eqn_idx%adv%beg)%sf(i, j, &
                             & k)/blkmod1 + (1._wp - q_prim_vf(eqn_idx%adv%beg)%sf(i, j, k))/blkmod2)))
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

    !> Derive the flux limiter at cell boundary i+1/2. This is an approximation because the velocity used to determine the upwind
    !! direction is the velocity at the cell center i instead of the contact velocity at the cell boundary from the Riemann solver.
    subroutine s_derive_flux_limiter(i, q_prim_vf, q_sf)

        integer, intent(in)                                 :: i
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end), &
             & intent(inout) :: q_sf

        real(wp) :: top, bottom, slope
        integer  :: j, k, l
        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end
                    if (i == 1) then
                        if (q_prim_vf(eqn_idx%cont%end + i)%sf(j, k, l) >= 0._wp) then
                            top = q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l) - q_prim_vf(eqn_idx%adv%beg)%sf(j - 1, k, l)
                            bottom = q_prim_vf(eqn_idx%adv%beg)%sf(j + 1, k, l) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l)
                        else
                            top = q_prim_vf(eqn_idx%adv%beg)%sf(j + 2, k, l) - q_prim_vf(eqn_idx%adv%beg)%sf(j + 1, k, l)
                            bottom = q_prim_vf(eqn_idx%adv%beg)%sf(j + 1, k, l) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l)
                        end if
                    else if (i == 2) then
                        if (q_prim_vf(eqn_idx%cont%end + i)%sf(j, k, l) >= 0._wp) then
                            top = q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k - 1, l)
                            bottom = q_prim_vf(eqn_idx%adv%beg)%sf(j, k + 1, l) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l)
                        else
                            top = q_prim_vf(eqn_idx%adv%beg)%sf(j, k + 2, l) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k + 1, l)
                            bottom = q_prim_vf(eqn_idx%adv%beg)%sf(j, k + 1, l) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l)
                        end if
                    else
                        if (q_prim_vf(eqn_idx%cont%end + i)%sf(j, k, l) >= 0._wp) then
                            top = q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l - 1)
                            bottom = q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l + 1) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l)
                        else
                            top = q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l + 2) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l + 1)
                            bottom = q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l + 1) - q_prim_vf(eqn_idx%adv%beg)%sf(j, k, l)
                        end if
                    end if

                    if (abs(top) < 1.e-8_wp) top = 0._wp
                    if (abs(bottom) < 1.e-8_wp) bottom = 0._wp

                    if (f_approx_equal(top, bottom)) then
                        slope = 1._wp
                    else
                        slope = (top*bottom)/(bottom**2._wp + 1.e-16_wp)
                    end if

                    if (flux_lim == 1) then  ! MINMOD (MM)
                        q_sf(j, k, l) = max(0._wp, min(1._wp, slope))
                    else if (flux_lim == 2) then  ! MUSCL (MC)
                        q_sf(j, k, l) = max(0._wp, min(2._wp*slope, 5.e-1_wp*(1._wp + slope), 2._wp))
                    else if (flux_lim == 3) then  ! OSPRE (OP)
                        q_sf(j, k, l) = (15.e-1_wp*(slope**2._wp + slope))/(slope**2._wp + slope + 1._wp)
                    else if (flux_lim == 4) then  ! SUPERBEE (SB)
                        q_sf(j, k, l) = max(0._wp, min(1._wp, 2._wp*slope), min(slope, 2._wp))
                    else if (flux_lim == 5) then  ! SWEBY (SW) (beta = 1.5)
                        q_sf(j, k, l) = max(0._wp, min(15.e-1_wp*slope, 1._wp), min(slope, 15.e-1_wp))
                    else if (flux_lim == 6) then  ! VAN ALBADA (VA)
                        q_sf(j, k, l) = (slope**2._wp + slope)/(slope**2._wp + 1._wp)
                    else if (flux_lim == 7) then  ! VAN LEER (VL)
                        q_sf(j, k, l) = (abs(slope) + slope)/(1._wp + abs(slope))
                    end if
                end do
            end do
        end do

    end subroutine s_derive_flux_limiter

    !> Compute the specified component of the vorticity from the primitive variables. From those inputs, it proceeds to calculate
    !! values of the desired vorticity component, which are subsequently stored in derived flow quantity storage variable, q_sf.
    subroutine s_derive_vorticity_component(i, q_prim_vf, q_sf)

        integer, intent(in)                                 :: i
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end), &
             & intent(inout) :: q_sf

        integer :: j, k, l, r
        if (i == 1) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end
                        q_sf(j, k, l) = 0._wp

                        do r = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                q_sf(j, k, l) = q_sf(j, k, l) + 1._wp/y_cc(k)*(fd_coeff_y(r, &
                                     & k)*y_cc(r + k)*q_prim_vf(eqn_idx%mom%end)%sf(j, r + k, l) - fd_coeff_z(r, &
                                     & l)*q_prim_vf(eqn_idx%mom%beg + 1)%sf(j, k, r + l))
                            else
                                q_sf(j, k, l) = q_sf(j, k, l) + fd_coeff_y(r, k)*q_prim_vf(eqn_idx%mom%end)%sf(j, r + k, &
                                     & l) - fd_coeff_z(r, l)*q_prim_vf(eqn_idx%mom%beg + 1)%sf(j, k, r + l)
                            end if
                        end do
                    end do
                end do
            end do
        else if (i == 2) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end
                        q_sf(j, k, l) = 0._wp

                        do r = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                q_sf(j, k, l) = q_sf(j, k, l) + fd_coeff_z(r, l)/y_cc(k)*q_prim_vf(eqn_idx%mom%beg)%sf(j, k, &
                                     & r + l) - fd_coeff_x(r, j)*q_prim_vf(eqn_idx%mom%end)%sf(r + j, k, l)
                            else
                                q_sf(j, k, l) = q_sf(j, k, l) + fd_coeff_z(r, l)*q_prim_vf(eqn_idx%mom%beg)%sf(j, k, &
                                     & r + l) - fd_coeff_x(r, j)*q_prim_vf(eqn_idx%mom%end)%sf(r + j, k, l)
                            end if
                        end do
                    end do
                end do
            end do
        else
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end
                        q_sf(j, k, l) = 0._wp

                        do r = -fd_number, fd_number
                            q_sf(j, k, l) = q_sf(j, k, l) + fd_coeff_x(r, j)*q_prim_vf(eqn_idx%mom%beg + 1)%sf(r + j, k, &
                                 & l) - fd_coeff_y(r, k)*q_prim_vf(eqn_idx%mom%beg)%sf(j, r + k, l)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_derive_vorticity_component

    !> Compute the Q_M criterion from the primitive variables. The Q_M function, which are subsequently stored in the derived flow
    !! quantity storage variable, q_sf.
    subroutine s_derive_qm(q_prim_vf, q_sf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end), &
             & intent(inout) :: q_sf

        real(wp), dimension(1:3,1:3) :: q_jacobian_sf, S, S2, O, O2
        real(wp)                     :: trS, Q, IIS
        integer                      :: j, k, l, r, jj, kk
        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end
                    ! Get velocity gradient tensor
                    q_jacobian_sf(:,:) = 0._wp

                    do r = -fd_number, fd_number
                        do jj = 1, 3
                            ! d()/dx
                            q_jacobian_sf(jj, 1) = q_jacobian_sf(jj, 1) + fd_coeff_x(r, &
                                          & j)*q_prim_vf(eqn_idx%mom%beg + jj - 1)%sf(r + j, k, l)
                            ! d()/dy
                            q_jacobian_sf(jj, 2) = q_jacobian_sf(jj, 2) + fd_coeff_y(r, &
                                          & k)*q_prim_vf(eqn_idx%mom%beg + jj - 1)%sf(j, r + k, l)
                            ! d()/dz
                            q_jacobian_sf(jj, 3) = q_jacobian_sf(jj, 3) + fd_coeff_z(r, &
                                          & l)*q_prim_vf(eqn_idx%mom%beg + jj - 1)%sf(j, k, r + l)
                        end do
                    end do

                    ! Decompose velocity gradient into symmetric strain-rate S and skew-symmetric rotation-rate O
                    do jj = 1, 3
                        do kk = 1, 3
                            S(jj, kk) = 0.5_wp*(q_jacobian_sf(jj, kk) + q_jacobian_sf(kk, jj))
                            O(jj, kk) = 0.5_wp*(q_jacobian_sf(jj, kk) - q_jacobian_sf(kk, jj))
                        end do
                    end do

                    do jj = 1, 3
                        do kk = 1, 3
                            O2(jj, kk) = O(jj, 1)*O(kk, 1) + O(jj, 2)*O(kk, 2) + O(jj, 3)*O(kk, 3)
                            S2(jj, kk) = S(jj, 1)*S(kk, 1) + S(jj, 2)*S(kk, 2) + S(jj, 3)*S(kk, 3)
                        end do
                    end do

                    ! Q-criterion: Q = (||O||^2 - ||S||^2)/2, Hunt et al. CTR (1988)
                    Q = 0.5_wp*((O2(1, 1) + O2(2, 2) + O2(3, 3)) - (S2(1, 1) + S2(2, 2) + S2(3, 3)))
                    trS = S(1, 1) + S(2, 2) + S(3, 3)
                    ! Second invariant of strain-rate tensor
                    IIS = 0.5_wp*((S(1, 1) + S(2, 2) + S(3, 3))**2 - (S2(1, 1) + S2(2, 2) + S2(3, 3)))
                    q_sf(j, k, l) = Q + IIS
                end do
            end do
        end do

    end subroutine s_derive_qm

    !> Compute the Liutex vector and its magnitude based on Xu et al. (2019).
    impure subroutine s_derive_liutex(q_prim_vf, liutex_mag, liutex_axis)

        ! Liutex vortex identification via real eigenvector of velocity gradient, Xu et al. PoF (2019)

        integer, parameter                                  :: nm = 3
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        !> Liutex magnitude

        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end), &
             & intent(out) :: liutex_mag
        !> Liutex rigid rotation axis
        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end,nm), &
             & intent(out) :: liutex_axis
        character, parameter        :: ivl = 'N'     !< compute left eigenvectors
        character, parameter        :: ivr = 'V'     !< compute right eigenvectors
        real(wp), dimension(nm, nm) :: vgt           !< velocity gradient tensor
        real(wp), dimension(nm)     :: lr, li        !< real and imaginary parts of eigenvalues
        real(wp), dimension(nm, nm) :: vl, vr        !< left and right eigenvectors
        integer, parameter          :: lwork = 4*nm  !< size of work array (4*nm recommended)
        real(wp), dimension(lwork)  :: work          !< work array
        integer                     :: info
        real(wp), dimension(nm)     :: eigvec        !< real eigenvector
        real(wp)                    :: eigvec_mag    !< magnitude of real eigenvector
        real(wp)                    :: omega_proj    !< projection of vorticity on real eigenvector
        real(wp)                    :: lci           !< imaginary part of complex eigenvalue
        real(wp)                    :: alpha
        integer                     :: j, k, l, r, i
        integer                     :: idx

        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end
                    ! Get velocity gradient tensor (VGT)
                    vgt(:,:) = 0._wp

                    do r = -fd_number, fd_number
                        do i = 1, 3
                            ! d()/dx
                            vgt(i, 1) = vgt(i, 1) + fd_coeff_x(r, j)*q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(r + j, k, l)
                            ! d()/dy
                            vgt(i, 2) = vgt(i, 2) + fd_coeff_y(r, k)*q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, r + k, l)
                            ! d()/dz
                            vgt(i, 3) = vgt(i, 3) + fd_coeff_z(r, l)*q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, r + l)
                        end do
                    end do

                    ! Call appropriate LAPACK routine based on precision
#ifdef MFC_SINGLE_PRECISION
                    call sgeev(ivl, ivr, nm, vgt, nm, lr, li, vl, nm, vr, nm, work, lwork, info)
#else
                    call dgeev(ivl, ivr, nm, vgt, nm, lr, li, vl, nm, vr, nm, work, lwork, info)
#endif

                    ! Find eigenvector with smallest imaginary eigenvalue (real eigenvector of VGT)
                    idx = 1
                    do r = 2, 3
                        if (abs(li(r)) < abs(li(idx))) then
                            idx = r
                        end if
                    end do
                    eigvec = vr(:,idx)

                    ! Normalize real eigenvector if it is effectively non-zero
                    eigvec_mag = sqrt(eigvec(1)**2._wp + eigvec(2)**2._wp + eigvec(3)**2._wp)
                    if (eigvec_mag > sgm_eps) then
                        eigvec = eigvec/eigvec_mag
                    else
                        eigvec = 0._wp
                    end if

                    ! Compute vorticity projected on the eigenvector
                    omega_proj = (vgt(3, 2) - vgt(2, 3))*eigvec(1) + (vgt(1, 3) - vgt(3, 1))*eigvec(2) + (vgt(2, 1) - vgt(1, &
                                  & 2))*eigvec(3)

                    ! As eigenvector can have +/- signs, we can choose the sign so that omega_proj is positive
                    if (omega_proj < 0._wp) then
                        eigvec = -eigvec
                        omega_proj = -omega_proj
                    end if

                    ! Imaginary eigenvalue of the complex conjugate pair (cyclic index selection)
                    lci = li(mod(idx, 3) + 1)

                    ! Discriminant: determines whether rotation dominates strain
                    alpha = omega_proj**2._wp - 4._wp*lci**2._wp
                    ! Liutex magnitude = omega_proj - sqrt(discriminant) when rotation dominates
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

    !> Compute the values of the numerical Schlieren function, which are subsequently stored in the derived flow quantity storage
    !! variable, q_sf.
    impure subroutine s_derive_numerical_schlieren_function(q_cons_vf, q_sf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf

        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end), &
             & intent(inout) :: q_sf

        real(wp)               :: drho_dx, drho_dy, drho_dz  !< Spatial derivatives of the density in the x-, y- and z-directions
        real(wp), dimension(2) :: gm_rho_max                 !< Global (max gradient magnitude, rank) pair for density
        integer                :: i, j, k, l

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

        if (p > 0) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end
                        drho_dz = 0._wp

                        do i = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                drho_dz = drho_dz + fd_coeff_z(i, l)/y_cc(k)*rho_sf(j, k, i + l)
                            else
                                drho_dz = drho_dz + fd_coeff_z(i, l)*rho_sf(j, k, i + l)
                            end if
                        end do

                        gm_rho_sf(j, k, l) = gm_rho_sf(j, k, l) + drho_dz*drho_dz
                    end do
                end do
            end do
        end if

        gm_rho_sf = sqrt(gm_rho_sf)

        gm_rho_max = (/maxval(gm_rho_sf), real(proc_rank, wp)/)

        if (num_procs > 1) call s_mpi_reduce_maxloc(gm_rho_max)

        ! The form of the numerical Schlieren function depends on the choice of the multicomponent flow model. For the gamma/pi_inf
        ! model, the exponential of the negative, normalized, gradient magnitude of the density is computed. For the volume fraction
        ! model, the amplitude of the exponential's inside is also modulated with respect to the identity of the fluid in which the
        ! function is evaluated. For more information, refer to Marquina and Mulet (2003).

        if (model_eqns == 1) then  ! Gamma/pi_inf model
            q_sf = -gm_rho_sf/gm_rho_max(1)
        else  ! Volume fraction model
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end
                        q_sf(j, k, l) = 0._wp

                        do i = 1, eqn_idx%adv%end - eqn_idx%E
                            q_sf(j, k, l) = q_sf(j, k, l) - schlieren_alpha(i)*q_cons_vf(i + eqn_idx%E)%sf(j, k, l)*gm_rho_sf(j, &
                                 & k, l)/gm_rho_max(1)
                        end do
                    end do
                end do
            end do
        end if

        ! Up until now, only the inside of the exponential of the numerical Schlieren function has been evaluated and stored. Then,
        ! to finish the computation, the exponential of the inside quantity is taken.
        q_sf = exp(q_sf)

    end subroutine s_derive_numerical_schlieren_function

    !> Deallocation procedures for the module
    impure subroutine s_finalize_derived_variables_module

        ! Deallocating the variable containing the gradient magnitude of the density field provided that the numerical Schlieren
        ! function was was outputted during the post-process
        if (schlieren_wrt) deallocate (gm_rho_sf)

        ! Deallocating the variables that might have been used to bookkeep the finite-difference coefficients in the x-, y- and
        ! z-directions
        if (allocated(fd_coeff_x)) deallocate (fd_coeff_x)
        if (allocated(fd_coeff_y)) deallocate (fd_coeff_y)
        if (allocated(fd_coeff_z)) deallocate (fd_coeff_z)

    end subroutine s_finalize_derived_variables_module

end module m_derived_variables
