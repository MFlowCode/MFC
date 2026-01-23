!>
!! @file m_perturbation.fpp
!! @brief Contains module m_perturbation

!> @brief This module contains subroutines that compute perturbations to the
!!              initial mean flow fields.
module m_perturbation

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_boundary_common   ! Boundary conditions module

    use m_helper

    use m_simplex_noise

    use ieee_arithmetic

    implicit none

    real(wp), allocatable, dimension(:, :, :, :) :: q_prim_temp

contains

    impure subroutine s_initialize_perturbation_module()

        if (elliptic_smoothing) then
            allocate (q_prim_temp(0:m, 0:n, 0:p, 1:sys_size))
        end if

    end subroutine s_initialize_perturbation_module

    impure subroutine s_perturb_sphere(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k, l !< generic loop operators

        real(wp) :: perturb_alpha

        real(wp) :: rand_real
        call random_seed()

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    call random_number(rand_real)

                    perturb_alpha = q_prim_vf(E_idx + perturb_sph_fluid)%sf(i, j, k)

                    ! Perturb partial density fields to match perturbed volume fraction fields
                    !    IF ((perturb_alpha >= 25e-2_wp) .AND. (perturb_alpha <= 75e-2_wp)) THEN
                    if ((.not. f_approx_equal(perturb_alpha, 0._wp)) .and. (.not. f_approx_equal(perturb_alpha, 1._wp))) then

                        ! Derive new partial densities
                        do l = 1, num_fluids
                            q_prim_vf(l)%sf(i, j, k) = q_prim_vf(E_idx + l)%sf(i, j, k)*fluid_rho(l)
                        end do

                    end if
                end do
            end do
        end do

    end subroutine s_perturb_sphere

    impure subroutine s_perturb_surrounding_flow(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k !<  generic loop iterators

        real(wp) :: perturb_alpha
        real(wp) :: rand_real
        call random_seed()

        ! Perturb partial density or velocity of surrounding flow by some random small amount of noise
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    perturb_alpha = q_prim_vf(E_idx + perturb_flow_fluid)%sf(i, j, k)
                    call random_number(rand_real)
                    rand_real = rand_real*perturb_flow_mag
                    q_prim_vf(mom_idx%beg)%sf(i, j, k) = (1._wp + rand_real)*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    q_prim_vf(mom_idx%end)%sf(i, j, k) = rand_real*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    if (bubbles_euler) then
                        q_prim_vf(alf_idx)%sf(i, j, k) = (1._wp + rand_real)*q_prim_vf(alf_idx)%sf(i, j, k)
                    end if
                end do
            end do
        end do
    end subroutine s_perturb_surrounding_flow

    impure subroutine s_elliptic_smoothing(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type
        integer :: i, j, k, l, q

        do q = 1, elliptic_smoothing_iters

            ! Communication of buffer regions and apply boundary conditions
            call s_populate_variables_buffers(bc_type, q_prim_vf, pb%sf, mv%sf)

            ! Perform smoothing and store in temp array
            if (n == 0) then
                do j = 0, m
                    do i = 1, sys_size
                        q_prim_temp(j, 0, 0, i) = (1._wp/4._wp)* &
                                                  (q_prim_vf(i)%sf(j + 1, 0, 0) + q_prim_vf(i)%sf(j - 1, 0, 0) + &
                                                   2._wp*q_prim_vf(i)%sf(j, 0, 0))
                    end do
                end do
            else if (p == 0) then
                do k = 0, n
                    do j = 0, m
                        do i = 1, sys_size
                            q_prim_temp(j, k, 0, i) = (1._wp/8._wp)* &
                                                      (q_prim_vf(i)%sf(j + 1, k, 0) + q_prim_vf(i)%sf(j - 1, k, 0) + &
                                                       q_prim_vf(i)%sf(j, k + 1, 0) + q_prim_vf(i)%sf(j, k - 1, 0) + &
                                                       4._wp*q_prim_vf(i)%sf(j, k, 0))
                        end do
                    end do
                end do
            else
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do i = 1, sys_size
                                q_prim_temp(j, k, l, i) = (1._wp/12._wp)* &
                                                          (q_prim_vf(i)%sf(j + 1, k, l) + q_prim_vf(i)%sf(j - 1, k, l) + &
                                                           q_prim_vf(i)%sf(j, k + 1, l) + q_prim_vf(i)%sf(j, k - 1, l) + &
                                                           q_prim_vf(i)%sf(j, k, l + 1) + q_prim_vf(i)%sf(j, k, l - 1) + &
                                                           6._wp*q_prim_vf(i)%sf(j, k, l))
                            end do
                        end do
                    end do
                end do
            end if

            ! Copy smoothed data back to array of scalar fields
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        do i = 1, sys_size
                            q_prim_vf(i)%sf(j, k, l) = q_prim_temp(j, k, l, i)
                        end do
                    end do
                end do
            end do
        end do

    end subroutine s_elliptic_smoothing

    subroutine s_perturb_simplex(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp) :: mag, freq, scale, vel_rsm
        real(wp), dimension(:, :), allocatable :: ofs
        integer :: nOffsets
        real(wp) :: xl, yl, zl

        integer :: i, j, k, l, q

        nOffsets = max(num_dims, num_fluids)

        allocate (ofs(nOffsets, num_dims))

        ! Store offsets
        do i = 1, num_dims
            do j = 1, num_dims
                ofs(j, i) = simplex_params%perturb_vel_offset(j, i)
            end do
        end do

        ! Perturb velocities
        do i = 1, num_dims
            if (simplex_params%perturb_vel(i)) then
                freq = simplex_params%perturb_vel_freq(i)
                scale = simplex_params%perturb_vel_scale(i)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            xl = freq*(x_cc(j) + ofs(i, 1))
                            yl = freq*(y_cc(k) + ofs(i, 2))
                            if (num_dims == 2) then
                                mag = f_simplex2d(xl, yl)
                            elseif (num_dims == 3) then
                                zl = freq*(z_cc(l) + ofs(i, 3))
                                mag = f_simplex3d(xl, yl, zl)
                            end if

                            vel_rsm = 0._wp
                            do q = 1, num_dims
                                vel_rsm = vel_rsm + q_prim_vf(momxb + q - 1)%sf(j, k, l)**2._wp
                            end do
                            vel_rsm = sqrt(vel_rsm)

                            q_prim_vf(momxb + i - 1)%sf(j, k, l) = q_prim_vf(momxb + i - 1)%sf(j, k, l) + &
                                                                   vel_rsm*scale*mag
                        end do
                    end do
                end do
            end if
        end do

        ! Store offsets
        do i = 1, num_dims
            do j = 1, num_fluids
                ofs(j, i) = simplex_params%perturb_dens_offset(j, i)
            end do
        end do

        ! Perturb densities
        do i = 1, num_fluids
            if (simplex_params%perturb_dens(i)) then
                freq = simplex_params%perturb_dens_freq(i)
                scale = simplex_params%perturb_dens_scale(i)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            xl = freq*(x_cc(j) + ofs(i, 1))
                            yl = freq*(y_cc(k) + ofs(i, 2))
                            if (num_dims == 2) then
                                mag = f_simplex2d(xl, yl)
                            elseif (num_dims == 3) then
                                zl = freq*(z_cc(l) + ofs(i, 3))
                                mag = f_simplex3d(xl, yl, zl)
                            end if
                            q_prim_vf(contxb + i - 1)%sf(j, k, l) = q_prim_vf(contxb + i - 1)%sf(j, k, l) + &
                                                                    q_prim_vf(contxb + i - 1)%sf(j, k, l)*scale*mag
                        end do
                    end do
                end do
            end if
        end do

        deallocate (ofs)

    end subroutine s_perturb_simplex

    !>  This subroutine computes velocity perturbations for a temporal mixing
        !!              layer with a hyperbolic tangent mean streamwise velocity
        !!              profile, using an inverter version of the spectrum-based
        !!              synthetic turbulence generation method proposed by
        !!              Guo et al. (2023, JFM).
    subroutine s_perturb_mixlayer(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(mixlayer_perturb_nk) :: k, Ek
        real(wp), dimension(3, 3) :: Rij, Lmat
        real(wp), dimension(3) :: velfluc, sig_tmp, sig, khat, xi
        real(wp) :: dk, alpha, Eksum, q, uu0, phi
        integer :: i, j, l, r, ierr

        ! Initialize parameters
        dk = 1._wp/mixlayer_perturb_nk

        ! Compute prescribed energy spectra
        Eksum = 0._wp
        do i = 1, mixlayer_perturb_nk
            k(i) = dk*i
            Ek(i) = (k(i)/mixlayer_perturb_k0)**4._wp*exp(-2._wp*(k(i)/mixlayer_perturb_k0)**2._wp)
            Eksum = Eksum + Ek(i)
        end do

        ! Main loop
        do r = 0, n
            ! Compute prescribed Reynolds stress tensor with about half
            ! magnitude of its self-similar value
            Rij(:, :) = 0._wp
            uu0 = patch_icpp(1)%vel(1)**2._wp &
                  *(1._wp - tanh(y_cc(r)*mixlayer_vel_coef)**2._wp)
            Rij(1, 1) = 0.05_wp*uu0
            Rij(2, 2) = 0.03_wp*uu0
            Rij(3, 3) = 0.03_wp*uu0
            Rij(1, 2) = -0.02_wp*uu0
            Rij(2, 1) = Rij(1, 2)

            ! Cholesky decomposition for inhomogeneity and anisotropy
            Lmat = 0._wp
            Lmat(1, 1) = sqrt(Rij(1, 1))
            if (abs(Lmat(1, 1)) < sgm_eps) Lmat(1, 1) = sgm_eps
            Lmat(2, 1) = Rij(2, 1)/Lmat(1, 1)
            Lmat(2, 2) = sqrt(Rij(2, 2) - Lmat(2, 1)**2._wp)
            if (abs(Lmat(2, 2)) < sgm_eps) Lmat(2, 2) = sgm_eps
            Lmat(3, 1) = Rij(3, 1)/Lmat(1, 1)
            Lmat(3, 2) = (Rij(3, 2) - Lmat(3, 1)*Lmat(2, 1))/Lmat(2, 2)
            Lmat(3, 3) = sqrt(Rij(3, 3) - Lmat(3, 1)**2._wp - Lmat(3, 2)**2._wp)

            ! Compute perturbation for each Fourier component
            do i = 1, mixlayer_perturb_nk
                ! Generate random numbers for unit wavevector khat,
                ! random unit vector xi, and random mode phase phi
                if (proc_rank == 0) then
                    call s_generate_random_perturbation(khat, xi, phi, i, y_cc(r))
                end if

#ifdef MFC_MPI
                call MPI_BCAST(khat, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
                call MPI_BCAST(xi, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
                call MPI_BCAST(phi, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

                ! Compute mode direction by two-time cross product
                sig_tmp = f_cross(xi, khat)
                sig_tmp = sig_tmp/sqrt(sum(sig_tmp**2._wp))
                sig = f_cross(khat, sig_tmp)

                ! Compute perturbation for each grid
                do l = 0, p
                    do j = 0, m
                        q = sqrt(Ek(i)/Eksum)
                        alpha = k(i)*(khat(1)*x_cc(j) + khat(2)*y_cc(r) + khat(3)*z_cc(l)) + 2._wp*pi*phi
                        velfluc = 2._wp*q*sig*cos(alpha)
                        velfluc = matmul(Lmat, velfluc)
                        q_prim_vf(momxb)%sf(j, r, l) = q_prim_vf(momxb)%sf(j, r, l) + velfluc(1)
                        q_prim_vf(momxb + 1)%sf(j, r, l) = q_prim_vf(momxb + 1)%sf(j, r, l) + velfluc(2)
                        q_prim_vf(momxb + 2)%sf(j, r, l) = q_prim_vf(momxb + 2)%sf(j, r, l) + velfluc(3)
                    end do
                end do
            end do
        end do

    end subroutine s_perturb_mixlayer

    subroutine s_generate_random_perturbation(khat, xi, phi, ik, yloc)
        integer, intent(in) :: ik
        real(wp), intent(in) :: yloc
        real(wp), dimension(3), intent(out) :: khat, xi
        real(wp), intent(out) :: phi
        real(wp) :: theta, eta
        integer :: seed, kfac, yfac

        kfac = ik*amplifier
        yfac = nint((sin(yloc) + 1._wp)*amplifier)
        seed = nint(0.5_wp*modmul(kfac) + 0.5_wp*modmul(yfac))

        call s_prng(theta, seed)
        call s_prng(eta, seed)
        khat = f_unit_vector(theta, eta)

        call s_prng(theta, seed)
        call s_prng(eta, seed)
        xi = f_unit_vector(theta, eta)

        call s_prng(phi, seed)

    end subroutine s_generate_random_perturbation

    ! Generate a random unit vector (spherical distribution)
    function f_unit_vector(theta, eta) result(vec)
        real(wp), intent(in) :: theta, eta
        real(wp) :: zeta, xi
        real(wp), dimension(3) :: vec

        xi = 2._wp*pi*theta
        zeta = acos(2._wp*eta - 1._wp)
        vec(1) = sin(zeta)*cos(xi)
        vec(2) = sin(zeta)*sin(xi)
        vec(3) = cos(zeta)

    end function f_unit_vector

    !>  This function generates a pseudo-random number between 0 and 1 based on
    !!  linear congruential generator.
    subroutine s_prng(var, seed)
        integer, intent(inout) :: seed
        real(wp), intent(out) :: var
        integer :: i

        seed = mod(modmul(seed), modulus)
        var = seed/real(modulus, wp)

    end subroutine s_prng

    function modmul(a) result(val)
        integer, intent(in) :: a
        integer :: val
        real(wp) :: x, y

        x = (multiplier/real(modulus, wp))*a + (increment/real(modulus, wp))
        y = nint((x - floor(x))*decimal_trim)/decimal_trim
        val = nint(y*modulus)

    end function modmul

    impure subroutine s_finalize_perturbation_module()

        if (elliptic_smoothing) then
            deallocate (q_prim_temp)
        end if

    end subroutine s_finalize_perturbation_module

end module m_perturbation
