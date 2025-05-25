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
    
    use ieee_arithmetic

    implicit none

    real(wp), allocatable, dimension(:, :, :, :) :: q_prim_temp

contains

    subroutine s_initialize_perturbation_module()

        if (elliptic_smoothing) then
            allocate (q_prim_temp(0:m, 0:n, 0:p, 1:sys_size))
        end if

    end subroutine s_initialize_perturbation_module

    subroutine s_perturb_sphere(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k, l !< generic loop operators

        real(wp) :: perturb_alpha
        real(wp) :: alpha_unadv
        real(wp) :: rand_real
        call random_seed()

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    call random_number(rand_real)

                    perturb_alpha = q_prim_vf(E_idx + perturb_sph_fluid)%sf(i, j, k)

                    ! Perturb partial density fields to match perturbed volume fraction fields
                    !    IF ((perturb_alpha >= 25e-2_wp) .AND. (perturb_alpha <= 75e-2_wp)) THEN
                    if ((perturb_alpha /= 0._wp) .and. (perturb_alpha /= 1._wp)) then

                        ! Derive new partial densities
                        do l = 1, num_fluids
                            q_prim_vf(l)%sf(i, j, k) = q_prim_vf(E_idx + l)%sf(i, j, k)*fluid_rho(l)
                        end do

                    end if
                end do
            end do
        end do

    end subroutine s_perturb_sphere

    subroutine s_perturb_surrounding_flow(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k, l !<  generic loop iterators

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

    subroutine s_elliptic_smoothing(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type
        integer :: i, j, k, l, q

        do q = 1, elliptic_smoothing_iters

            ! Communication of buffer regions and apply boundary conditions
            call s_populate_variables_buffers(q_prim_vf, pb%sf, mv%sf, bc_type)

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

    !>  This subroutine computes velocity perturbations for a temporal mixing
        !!              layer with a hyperbolic tangent mean streamwise velocity
        !!              profile, using an inverter version of the spectrum-based
        !!              synthetic turbulence generation method proposed by
        !!              Guo et al. (2023, JFM).
    subroutine s_perturb_mixlayer(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(mixlayer_perturb_nk) :: k, Ek
        real(wp), dimension(3,3) :: Rij, Lmat
        real(wp), dimension(3) :: velfluc, sig_tmp, sig
        real(wp) :: dk, k0, alpha, Eksum, q, uu0
        real(wp), dimension(3,0:n,mixlayer_perturb_nk) :: khat, xi
        real(wp), dimension(0:n,mixlayer_perturb_nk) :: phi
        integer :: i, j, l, r, ierr

        ! Initialize random number        
        call s_initialize_perturb_mixlayer(khat,xi,phi)

        ! Initialize parameters
        dk = 1._wp/mixlayer_perturb_nk
        k0 = 0.4446_wp  ! Most unstable mode obtained from linear stability
                        ! analysis. See Michalke (1964, JFM) for details

        ! Compute pre-determined energy spectra
        Eksum = 0_wp
        do i = 1, mixlayer_perturb_nk
            k(i) = dk*i
            Ek(i) = (k(i)/k0)**4._wp*exp(-2._wp*(k(i)/k0)**2._wp)
            Eksum = Eksum + Ek(i)
        end do

        ! Main loop
        do r = 0, n
            ! Compute pre-determined Reynolds stress tensor with about half
            ! magnitude of its self-similar value
            Rij(:, :) = 0_wp
            uu0 = patch_icpp(1)%vel(1)**2._wp &
                  *(1_wp - tanh(y_cc(r)*mixlayer_vel_coef)**2._wp)
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
            Lmat(2, 2) = sqrt(Rij(2, 2) - Lmat(2, 1)**2_wp)
            if (abs(Lmat(2, 2)) < sgm_eps) Lmat(2, 2) = sgm_eps
            Lmat(3, 1) = Rij(3, 1)/Lmat(1, 1)
            Lmat(3, 2) = (Rij(3, 2) - Lmat(3, 1)*Lmat(2, 1))/Lmat(2, 2)
            Lmat(3, 3) = sqrt(Rij(3, 3) - Lmat(3, 1)**2_wp - Lmat(3, 2)**2_wp)

            ! Compute perturbation for each Fourier component
            do i = 1, mixlayer_perturb_nk
!                 ! Generate random numbers for unit wavevector khat,
!                 ! random unit vector xi, and random mode phase phi
!                 if (proc_rank == 0) then
!                     khat = f_random_unit_vector()
!                     xi = f_random_unit_vector()
!                     phi = f_prng()
!                 end if
        
! #ifdef MFC_MPI
!                 call MPI_BCAST(khat, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
!                 call MPI_BCAST(xi, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
!                 call MPI_BCAST(phi, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
! #endif

                ! Compute mode direction by two-time cross product
                sig_tmp = f_cross(xi(:,r,i), khat(:,r,i))
                sig_tmp = sig_tmp/sqrt(sum(sig_tmp**2))
                sig = f_cross(khat(:,r,i), sig_tmp)

                ! Compute perturbation for each grid
                do l = 0, p
                    do j = 0, m
                        q = sqrt(Ek(i)/Eksum)
                        alpha = k(i)*(khat(1,r,i)*x_cc(j) + khat(2,r,i)*y_cc(r) + khat(3,r,i)*z_cc(l)) + 2_wp*pi*phi(r,i)
                        velfluc = 2_wp*q*sig*cos(alpha)
                        velfluc = matmul(Lmat, velfluc)
                        q_prim_vf(momxb)%sf(j, r, l) = q_prim_vf(momxb)%sf(j, r, l) + velfluc(1)
                        q_prim_vf(momxb + 1)%sf(j, r, l) = q_prim_vf(momxb + 1)%sf(j, r, l) + velfluc(2)
                        q_prim_vf(momxb + 2)%sf(j, r, l) = q_prim_vf(momxb + 2)%sf(j, r, l) + velfluc(3)
                    end do
                end do
            end do
        end do

        print *, proc_rank, j, r, l, q_prim_vf(momxb + 2)%sf(j, r, l)

    end subroutine s_perturb_mixlayer


    subroutine s_initialize_perturb_mixlayer(khat, xi, phi)

        real(wp), dimension(3,0:n,mixlayer_perturb_nk), intent(out) :: khat, xi
        real(wp), dimension(0:n,mixlayer_perturb_nk), intent(out) :: phi
        integer :: i, j, k, l

        if (proc_rank == 0) then
            do j = 0, n
                do i = 1, mixlayer_perturb_nk
                    khat(:,j,i) = f_random_unit_vector()
                    xi(:,j,i) = f_random_unit_vector()
                    phi(j,i) = f_prng()
                end do
            end do
        end if

    end subroutine s_initialize_perturb_mixlayer

    ! Generate a random unit vector (spherical distribution)
    function f_random_unit_vector() result(vec)
        real(wp) :: vec(3)
        real(wp) :: theta, phi

        theta = f_prng()
        phi = f_prng()
        theta = 2.0_wp*pi*theta
        phi = acos(2.0_wp*phi - 1.0_wp)
        vec(1) = sin(phi)*cos(theta)
        vec(2) = sin(phi)*sin(theta)
        vec(3) = cos(phi)

    end function f_random_unit_vector

    !>  This function generates a pseudo-random number between 0 and 1 based on
    !!  linear congruential generator. 
    function f_prng() result(val)
        
        real(wp) :: val

        prng_seed = mod(multiplier*prng_seed + increment, modulus)
        val = prng_seed / modulus

    end function f_prng

    subroutine s_finalize_perturbation_module()

        if (elliptic_smoothing) then
            deallocate (q_prim_temp)
        end if

    end subroutine s_finalize_perturbation_module

end module m_perturbation
