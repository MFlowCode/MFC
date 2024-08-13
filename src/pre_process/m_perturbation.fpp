!>
!! @file m_perturbation.fpp
!! @brief Contains module m_perturbation

!> @brief This module contains subroutines that compute perturbations to the
!!              initial mean flow fields.
module m_perturbation

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_eigen_solver          ! Subroutines to solve eigenvalue problem for
    ! complex general matrix

    use ieee_arithmetic

    ! ==========================================================================

    implicit none

    integer :: mixlayer_nvar ! Number of variables in linear stability analysis solver for mixing layer
    integer, allocatable, dimension(:) :: mixlayer_var ! Index of variables in linear stability analysis solver
    integer :: nbp ! Number of grid cell boundary points in y-direction
    integer :: mixlayer_bc_fd ! Order of finite difference applied at the boundaries of mixing layer
    integer :: n_bc_skip ! Number of points skipped in the linear stability analysis due to the boundary condition

contains

    subroutine s_initialize_perturbation_module()

        if (mixlayer_perturb) then
            mixlayer_bc_fd = 2
            nbp = n + 2
            if (model_eqns == 2 .and. num_fluids == 1) then
                n_bc_skip = mixlayer_bc_fd*2
                mixlayer_nvar = 5 ! 1 continuity + 3 momentum + 1 energy
                allocate (mixlayer_var(mixlayer_nvar))

                mixlayer_var(1) = contxb
                mixlayer_var(2) = momxb
                mixlayer_var(3) = momxb + 1
                mixlayer_var(4) = momxb + 2
                mixlayer_var(5) = momxb + 3
            end if
        end if

    end subroutine s_initialize_perturbation_module

    subroutine s_perturb_sphere(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k, l !< generic loop operators

        real(kind(0d0)) :: perturb_alpha
        real(kind(0d0)) :: alpha_unadv
        real(kind(0d0)) :: rand_real
        call random_seed()

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    call random_number(rand_real)

                    perturb_alpha = q_prim_vf(E_idx + perturb_sph_fluid)%sf(i, j, k)

                    ! Perturb partial density fields to match perturbed volume fraction fields
                    !    IF ((perturb_alpha >= 25d-2) .AND. (perturb_alpha <= 75d-2)) THEN
                    if ((perturb_alpha /= 0d0) .and. (perturb_alpha /= 1d0)) then

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

        real(kind(0d0)) :: perturb_alpha
        real(kind(0d0)) :: rand_real
        call random_seed()

        ! Perturb partial density or velocity of surrounding flow by some random small amount of noise
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    perturb_alpha = q_prim_vf(E_idx + perturb_flow_fluid)%sf(i, j, k)
                    call random_number(rand_real)
                    rand_real = rand_real*perturb_flow_mag
                    q_prim_vf(mom_idx%beg)%sf(i, j, k) = (1.d0 + rand_real)*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    q_prim_vf(mom_idx%end)%sf(i, j, k) = rand_real*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    if (bubbles) then
                        q_prim_vf(alf_idx)%sf(i, j, k) = (1.d0 + rand_real)*q_prim_vf(alf_idx)%sf(i, j, k)
                    end if
                end do
            end do
        end do

    end subroutine s_perturb_surrounding_flow

    !>  This subroutine computes velocity perturbations for a temporal mixing
        !!              layer with hypertangent mean streamwise velocity profile
        !!              obtained from linear stability analysis. For a 2D case,
        !!              instability waves with spatial wavenumbers, (4,0), (2,0),
        !!              and (1,0) are superposed. For a 3D waves, (4,4), (4,-4),
        !!              (2,2), (2,-2), (1,1), (1,-1) areadded on top of 2D waves.
    subroutine s_superposition_instability_wave(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(mixlayer_nvar, 0:m, 0:n, 0:p) :: wave, wave1, wave2, wave_tmp
        real(kind(0d0)) :: uratio, Ldomain
        integer :: i, j, k, q

        uratio = 1d0/patch_icpp(1)%vel(1)
        Ldomain = mixlayer_domain*patch_icpp(1)%length_y

        wave = 0d0
        wave1 = 0d0
        wave2 = 0d0

        ! Compute 2D waves
        call s_instability_wave(2*pi*4.0/Ldomain, 0d0, wave_tmp, 0d0)
        wave1 = wave1 + wave_tmp
        call s_instability_wave(2*pi*2.0/Ldomain, 0d0, wave_tmp, 0d0)
        wave1 = wave1 + wave_tmp
        call s_instability_wave(2*pi*1.0/Ldomain, 0d0, wave_tmp, 0d0)
        wave1 = wave1 + wave_tmp
        wave = wave1*0.05

        if (p > 0) then
            ! Compute 3D waves with phase shifts.
            call s_instability_wave(2*pi*4.0/Ldomain, 2*pi*4.0/Ldomain, wave_tmp, 2*pi*11d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*2.0/Ldomain, 2*pi*2.0/Ldomain, wave_tmp, 2*pi*13d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*1.0/Ldomain, 2*pi*1.0/Ldomain, wave_tmp, 2*pi*17d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*4.0/Ldomain, -2*pi*4.0/Ldomain, wave_tmp, 2*pi*19d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*2.0/Ldomain, -2*pi*2.0/Ldomain, wave_tmp, 2*pi*23d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*1.0/Ldomain, -2*pi*1.0/Ldomain, wave_tmp, 2*pi*29d0/31d0)
            wave2 = wave2 + wave_tmp
            wave = wave + 0.15*wave2
        end if

        ! Superpose velocity perturbuations (instability waves) to the velocity field
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    q_prim_vf(contxb)%sf(i, j, k) = q_prim_vf(contxb)%sf(i, j, k) + wave(mixlayer_var(1), i, j, k) ! rho
                    q_prim_vf(momxb)%sf(i, j, k) = q_prim_vf(momxb)%sf(i, j, k) + wave(mixlayer_var(2), i, j, k)/uratio ! u
                    q_prim_vf(momxb + 1)%sf(i, j, k) = q_prim_vf(momxb + 1)%sf(i, j, k) + wave(mixlayer_var(3), i, j, k)/uratio ! v
                    if (p > 0) then
                        q_prim_vf(momxb + 2)%sf(i, j, k) = q_prim_vf(momxb + 2)%sf(i, j, k) + wave(mixlayer_var(4), i, j, k)/uratio ! w
                    end if
                    q_prim_vf(E_idx)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k) + wave(mixlayer_var(5), i, j, k)/uratio**2 ! p

                    if (bubbles .and. (.not. qbmm)) then
                        do q = 1, nb
                            call s_compute_equilibrium_state(q_prim_vf(E_idx)%sf(i, j, k), R0(q), q_prim_vf(bub_idx%rs(q))%sf(i, j, k))
                        end do
                    end if
                end do
            end do
        end do

    end subroutine s_superposition_instability_wave ! ----------------------

    !>  This subroutine computes equilibrium bubble radius of the perturbed pressure field
    subroutine s_compute_equilibrium_state(fP, fR0, fR)
        real(kind(0d0)), intent(in) :: fP, fR0
        real(kind(0d0)), intent(inout) :: fR
        real(kind(0d0)) :: f0, f1
        real(kind(0d0)) :: gam_b
        integer :: ii, jj

        gam_b = 1d0 + 1d0/fluid_pp(num_fluids + 1)%gamma

        ! Loop
        ii = 1
        do while (.true.)

            f0 = (Ca + 2d0/Web)*(fR0/fR)**(3d0*gam_b) - 2d0/(Web*fR) + 1d0 - Ca - fP
            f1 = -3d0*gam_b*(Ca + 2d0/Web)*(fR0/fR)**(3d0*gam_b + 1d0) + 2d0/(Web*fR**2d0)

            if (abs(f0) <= 1e-10) then
                ! Converged
                exit
            else
                ! Update radius
                fR = fR - f0/f1
            end if

            ! Failed case
            if (ieee_is_nan(f0) .or. &
                ieee_is_nan(f1) .or. &
                ii > 1000 .or. &
                fR < 0d0) then

                print *, "Failed to compute equilibrium radius"

                fR = fR0
                exit
            end if

            ii = ii + 1
        end do

    end subroutine s_compute_equilibrium_state

    !>  This subroutine computes instability waves for a given set of spatial
        !!              wavenumbers (alpha, beta) in x and z directions.
        !!              The eigenvalue problem is derived from the linearized
        !!              Euler equations with parallel mean flow assumption
        !!              (See Sandham 1989 PhD thesis for details).
    subroutine s_instability_wave(alpha, beta, wave, shift)
        real(kind(0d0)), intent(in) :: alpha, beta !<  spatial wavenumbers
        real(kind(0d0)), dimension(mixlayer_nvar, 0:m, 0:n, 0:p), intent(inout) :: wave !< instability wave
        real(kind(0d0)) :: shift !< phase shift
        real(kind(0d0)), dimension(0:nbp - 1) :: u_mean !<  mean density and velocity profiles
        real(kind(0d0)) :: rho_mean, p_mean !< mean density and pressure
        real(kind(0d0)), dimension(0:nbp - 1, 0:nbp - 1) :: d !< differential operator in y dir
        real(kind(0d0)) :: gam, pi_inf, mach, c1, adv
        real(kind(0d0)) :: xratio, uratio
        integer :: i, j !<  generic loop iterators

        xratio = mixlayer_vel_coef
        uratio = 1d0/patch_icpp(1)%vel(1)

        ! Set fluid flow properties
        if (bubbles) then
            adv = patch_icpp(1)%alpha(num_fluids)
        else
            adv = 0d0
        end if
        gam = 1d0 + 1d0/fluid_pp(1)%gamma
        pi_inf = fluid_pp(1)%pi_inf*(gam - 1d0)/gam*uratio**2
        rho_mean = patch_icpp(1)%alpha_rho(1)
        p_mean = patch_icpp(1)%pres*uratio**2
        c1 = sqrt((gam*(p_mean + pi_inf))/(rho_mean*(1d0 - adv)))
        mach = 1d0/c1

        ! Assign mean profiles
        do j = 0, n + 1
            u_mean(j) = tanh(y_cb(j - 1)*xratio)
        end do

        ! Compute differential operator in y-dir
        ! based on 2nd order central difference
        d = 0d0
        d(0, 0) = -1d0/((y_cb(0) - y_cb(-1))*xratio)
        d(0, 1) = 1d0/((y_cb(0) - y_cb(-1))*xratio)
        do j = 1, n
            d(j, j - 1) = -1d0/((y_cb(j) - y_cb(j - 2))*xratio)
            d(j, j + 1) = 1d0/((y_cb(j) - y_cb(j - 2))*xratio)
        end do
        d(n + 1, n) = -1d0/((y_cb(n) - y_cb(n - 1))*xratio)
        d(n + 1, n + 1) = 1d0/((y_cb(n) - y_cb(n - 1))*xratio)

        ! Compute
        call s_solve_linear_system(alpha, beta, u_mean, rho_mean, p_mean, d, gam, pi_inf, mach, wave, shift)

    end subroutine s_instability_wave

    !> This subroutine solves linear system from linear stability analysis and
        !!              generate instability waves for the given set of spatial
        !!              wave numbers and phase shift.
    subroutine s_solve_linear_system(alpha, beta, u_mean, rho_mean, p_mean, d, gam, pi_inf, mach, wave, shift)
        real(kind(0d0)), intent(in) :: alpha, beta !<  spatial wavenumbers
        real(kind(0d0)), dimension(0:nbp - 1), intent(in) :: u_mean !<  mean velocity profiles
        real(kind(0d0)), intent(in) :: rho_mean, p_mean !< mean density and pressure
        real(kind(0d0)), dimension(0:nbp - 1, 0:nbp - 1), intent(in) :: d !< differential operator in y dir
        real(kind(0d0)), intent(in) :: gam, pi_inf, mach, shift
        real(kind(0d0)), dimension(mixlayer_nvar, 0:m, 0:n, 0:p), intent(inout) :: wave

        real(kind(0d0)), dimension(0:nbp - 1) :: drho_mean, du_mean !< y-derivatives of mean profiles
        real(kind(0d0)), dimension(0:mixlayer_nvar*nbp - 1, 0:mixlayer_nvar*nbp - 1) :: ar, ai    !< matrices for eigenvalue problem
        real(kind(0d0)), dimension(0:mixlayer_nvar*nbp - 1, 0:mixlayer_nvar*nbp - 1) :: br, bi, ci !< matrices for eigenvalue problem
        real(kind(0d0)), dimension(0:mixlayer_nvar*n - n_bc_skip - 1, 0:mixlayer_nvar*n - n_bc_skip - 1) :: hr, hi    !< matrices for eigenvalue problem

        real(kind(0d0)), dimension(0:mixlayer_nvar*n - n_bc_skip - 1, 0:mixlayer_nvar*n - n_bc_skip - 1) :: zr, zi !< eigenvectors
        real(kind(0d0)), dimension(0:mixlayer_nvar*n - n_bc_skip - 1) :: wr, wi !< eigenvalues
        real(kind(0d0)), dimension(0:mixlayer_nvar*n - n_bc_skip - 1) :: fv1, fv2, fv3 !< temporary memory

        integer :: ierr
        integer :: i, j, k, l !<  generic loop iterators
        integer :: ii, jj !< block matrix indices

        ! Compute y-derivatives of rho and u
        do j = 0, nbp - 1
            drho_mean(j) = 0
            du_mean(j) = 0
            do k = 0, nbp - 1
                drho_mean(j) = 0d0
                du_mean(j) = du_mean(j) + d(j, k)*u_mean(k)
            end do
        end do

        ! Compute B and C, then A = B + C. Here, A is the matrix for the linear
        ! systems of equation (i.e. we are going to solve x for Ax = lambda x).
        ! Here, B includes components of A without differential operator, and
        ! C includes components of A with differential operator.
        br = 0d0
        bi = 0d0
        ci = 0d0
        do j = 0, nbp - 1
            ii = mixlayer_var(1); jj = mixlayer_var(1); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = alpha*u_mean(j); 
            ii = mixlayer_var(1); jj = mixlayer_var(2); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = alpha*rho_mean; 
            ii = mixlayer_var(1); jj = mixlayer_var(3); bi((ii - 1)*nbp + j, (jj - 1)*nbp + j) = -drho_mean(j); 
            ii = mixlayer_var(1); jj = mixlayer_var(4); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = beta*rho_mean; 
            ii = mixlayer_var(2); jj = mixlayer_var(2); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = alpha*u_mean(j); 
            ii = mixlayer_var(2); jj = mixlayer_var(3); bi((ii - 1)*nbp + j, (jj - 1)*nbp + j) = -du_mean(j); 
            ii = mixlayer_var(2); jj = mixlayer_var(5); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = alpha/rho_mean; 
            ii = mixlayer_var(3); jj = mixlayer_var(3); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = alpha*u_mean(j); 
            ii = mixlayer_var(4); jj = mixlayer_var(4); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = alpha*u_mean(j); 
            ii = mixlayer_var(4); jj = mixlayer_var(5); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = beta/rho_mean; 
            ii = mixlayer_var(5); jj = mixlayer_var(2); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = gam*(p_mean + pi_inf)*alpha; 
            ii = mixlayer_var(5); jj = mixlayer_var(4); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = gam*(p_mean + pi_inf)*beta; 
            ii = mixlayer_var(5); jj = mixlayer_var(5); br((ii - 1)*nbp + j, (jj - 1)*nbp + j) = alpha*u_mean(j); 
            do k = 0, n + 1
                ii = mixlayer_var(1); jj = mixlayer_var(3); ci((ii - 1)*nbp + j, (jj - 1)*nbp + k) = -rho_mean*d(j, k); 
                ii = mixlayer_var(3); jj = mixlayer_var(5); ci((ii - 1)*nbp + j, (jj - 1)*nbp + k) = -d(j, k)/rho_mean; 
                ii = mixlayer_var(5); jj = mixlayer_var(3); ci((ii - 1)*nbp + j, (jj - 1)*nbp + k) = -gam*(p_mean + pi_inf)*d(j, k); 
            end do
        end do
        ar = br
        ai = bi + ci

        ! Apply BC to ar and ai matrices
        if (bc_y%beg == -6 .and. bc_y%end == -6) then
            ! Nonreflecting subsonic buffer BC
            call s_instability_nonreflecting_subsonic_buffer_bc(ar, ai, hr, hi, rho_mean, mach)
        end if

        ! Compute eigenvalues and eigenvectors
        call cg(mixlayer_nvar*n - n_bc_skip, mixlayer_nvar*n - n_bc_skip, hr, hi, wr, wi, zr, zi, fv1, fv2, fv3, ierr)

        ! Generate instability wave
        call s_generate_wave(wr, wi, zr, zi, rho_mean, mach, alpha, beta, wave, shift)

    end subroutine s_solve_linear_system

    !> This subroutine applies non-reflecting subsonic buffer boundary condition
        !!              to the linear system of equations (i.e. matrix A).
    subroutine s_instability_nonreflecting_subsonic_buffer_bc(ar, ai, hr, hi, rho_mean, mach)
        real(kind(0d0)), dimension(0:mixlayer_nvar*nbp - 1, 0:mixlayer_nvar*nbp - 1), intent(inout) :: ar, ai    !< matrices for eigenvalue problem
        real(kind(0d0)), dimension(0:mixlayer_nvar*n - n_bc_skip - 1, 0:mixlayer_nvar*n - n_bc_skip - 1), intent(out) :: hr, hi    !< matrices for eigenvalue problem
        real(kind(0d0)), intent(in) :: rho_mean !<  mean density profiles
        real(kind(0d0)), intent(in) :: mach
        real(kind(0d0)), dimension(0:mixlayer_nvar*n - 1, 0:mixlayer_nvar*n - 1) :: fr, fi    !< matrices for eigenvalue problem
        real(kind(0d0)), dimension(0:mixlayer_nvar*n - n_bc_skip - 1, 0:mixlayer_nvar*n - 1) :: gr, gi    !< matrices for eigenvalue problem
        integer :: i, j, k, l, ii, jj

        ! Condition 1: v = 0 at BC - no action required here

        ! Condition 2: du/dy = 0 at BC
        do j = 0, mixlayer_nvar*nbp - 1
            ! at y_domain%beg
            ii = mixlayer_var(1)*nbp
            ar(j, ii + 1) = ar(j, ii + 1) + ar(j, ii)
            ai(j, ii + 1) = ai(j, ii + 1) + ai(j, ii)
            ! at y_domain%end
            ii = mixlayer_var(1)*nbp + nbp - 1
            ar(j, ii - 1) = ar(j, ii - 1) + ar(j, ii)
            ai(j, ii - 1) = ai(j, ii - 1) + ai(j, ii)
        end do

        ! Condition 3: dw/dy = 0 at BC
        do j = 0, mixlayer_nvar*nbp - 1
            ! at y_domain%beg
            ii = (mixlayer_var(3))*nbp
            ar(j, ii + 1) = ar(j, ii + 1) + ar(j, ii)
            ai(j, ii + 1) = ai(j, ii + 1) + ai(j, ii)
            ! at y_domain%end
            ii = (mixlayer_var(3))*nbp + nbp - 1
            ar(j, ii - 1) = ar(j, ii - 1) + ar(j, ii)
            ai(j, ii - 1) = ai(j, ii - 1) + ai(j, ii)
        end do

        ! Condition 4: dp/dy +- rho c dv/dy = 0 at BC
        do j = 0, mixlayer_nvar*nbp - 1
            ! at y_domain%beg
            ii = mixlayer_var(4)*nbp
            ar(j, ii + 1) = ar(j, ii + 1) + ar(j, ii)
            ai(j, ii + 1) = ai(j, ii + 1) + ai(j, ii)
            jj = mixlayer_var(2)*nbp
            ar(j, jj + 1) = ar(j, jj + 1) + ar(j, ii)*rho_mean/mach
            ai(j, jj + 1) = ai(j, jj + 1) + ai(j, ii)*rho_mean/mach
            ! at y_domain%end
            ii = mixlayer_var(4)*nbp + nbp - 1
            ar(j, ii - 1) = ar(j, ii - 1) + ar(j, ii)
            ai(j, ii - 1) = ai(j, ii - 1) + ai(j, ii)
            jj = mixlayer_var(2)*nbp + nbp - 1
            ar(j, jj - 1) = ar(j, jj - 1) - ar(j, ii)*rho_mean/mach
            ai(j, jj - 1) = ai(j, jj - 1) - ai(j, ii)*rho_mean/mach
        end do

        ! Condition 5: c^2 drho/dy +- dp/dy = 0 at BC
        do j = 0, mixlayer_nvar*nbp - 1
            ! at y_domain%beg
            ii = 0
            ar(j, ii + 1) = ar(j, ii + 1) + ar(j, ii)
            ai(j, ii + 1) = ai(j, ii + 1) + ai(j, ii)
            jj = mixlayer_var(2)*nbp
            ar(j, jj + 1) = ar(j, jj + 1) + ar(j, ii)*rho_mean*mach
            ai(j, jj + 1) = ai(j, jj + 1) + ai(j, ii)*rho_mean*mach
            ! at y_domain%end
            ii = nbp - 1
            ar(j, ii - 1) = ar(j, ii - 1) + ar(j, ii)
            ai(j, ii - 1) = ai(j, ii - 1) + ai(j, ii)
            jj = mixlayer_var(2)*nbp + nbp - 1
            ar(j, jj - 1) = ar(j, jj - 1) - ar(j, ii)*rho_mean*mach
            ai(j, jj - 1) = ai(j, jj - 1) - ai(j, ii)*rho_mean*mach
        end do

        ! Remove unnecessary rows of the matrix A (rho, u, v, w, p at the boundaries)
        fr = 0d0
        fi = 0d0
        do ii = 1, mixlayer_nvar
            do jj = 1, mixlayer_nvar
                do k = 0, n - 1
                    do l = 0, n - 1
                        fr((ii - 1)*n + k, (jj - 1)*n + l) = ar((ii - 1)*nbp + k + 1, (jj - 1)*nbp + l + 1)
                        fi((ii - 1)*n + k, (jj - 1)*n + l) = ai((ii - 1)*nbp + k + 1, (jj - 1)*nbp + l + 1)
                    end do
                end do
            end do
        end do

        gr = 0d0
        gi = 0d0
        do ii = 1, mixlayer_nvar
            do j = 0, mixlayer_nvar*n - 1
                if (ii <= mixlayer_var(2)) then
                    do k = 0, n - 1
                        gr((ii - 1)*n + k, j) = fr((ii - 1)*n + k, j)
                        gi((ii - 1)*n + k, j) = fi((ii - 1)*n + k, j)
                    end do
                elseif (ii == mixlayer_var(3)) then
                    do k = 0, n - n_bc_skip - 1
                        gr((ii - 1)*n + k, j) = fr((ii - 1)*n + k + mixlayer_bc_fd, j)
                        gi((ii - 1)*n + k, j) = fi((ii - 1)*n + k + mixlayer_bc_fd, j)
                    end do
                else
                    do k = 0, n - 1
                        gr((ii - 1)*n - n_bc_skip + k, j) = fr((ii - 1)*n + k, j)
                        gi((ii - 1)*n - n_bc_skip + k, j) = fi((ii - 1)*n + k, j)
                    end do
                end if
            end do
        end do

        hr = 0d0
        hi = 0d0
        do i = 0, mixlayer_nvar*n - n_bc_skip - 1
            do jj = 1, mixlayer_nvar
                if (jj <= mixlayer_var(2)) then
                    do k = 0, n - 1
                        hr(i, (jj - 1)*n + k) = gr(i, (jj - 1)*n + k)
                        hi(i, (jj - 1)*n + k) = gi(i, (jj - 1)*n + k)
                    end do
                elseif (jj == mixlayer_var(3)) then
                    do k = 0, n - n_bc_skip - 1
                        hr(i, (jj - 1)*n + k) = gr(i, (jj - 1)*n + k + mixlayer_bc_fd)
                        hi(i, (jj - 1)*n + k) = gi(i, (jj - 1)*n + k + mixlayer_bc_fd)
                    end do
                else
                    do k = 0, n - 1
                        hr(i, (jj - 1)*n - n_bc_skip + k) = gr(i, (jj - 1)*n + k)
                        hi(i, (jj - 1)*n - n_bc_skip + k) = gi(i, (jj - 1)*n + k)
                    end do
                end if
            end do
        end do

    end subroutine s_instability_nonreflecting_subsonic_buffer_bc

    !>  This subroutine generates an instability wave using the most unstable
        !!              eigenvalue and corresponding eigenvector among the
        !!              given set of eigenvalues and eigenvectors.
    subroutine s_generate_wave(wr, wi, zr, zi, rho_mean, mach, alpha, beta, wave, shift)
        real(kind(0d0)), dimension(0:mixlayer_nvar*n - n_bc_skip - 1), intent(in) :: wr, wi !< eigenvalues
        real(kind(0d0)), dimension(0:mixlayer_nvar*n - n_bc_skip - 1, 0:mixlayer_nvar*n - n_bc_skip - 1), intent(in) :: zr, zi !< eigenvectors
        real(kind(0d0)), intent(in) :: rho_mean
        real(kind(0d0)), dimension(mixlayer_nvar, 0:m, 0:n, 0:p), intent(inout) :: wave
        real(kind(0d0)), intent(in) :: alpha, beta, mach, shift
        real(kind(0d0)), dimension(0:mixlayer_nvar*n - n_bc_skip - 1) :: vr, vi, vnr, vni !< most unstable eigenvector
        real(kind(0d0)), dimension(0:mixlayer_nvar*nbp - 1) :: xbr, xbi !< eigenvectors
        real(kind(0d0)), dimension(0:mixlayer_nvar*(nbp - 1) - 1) :: xcr, xci !< eigenvectors
        real(kind(0d0)) :: ang, norm
        real(kind(0d0)) :: tr, ti, cr, ci !< temporary memory
        real(kind(0d0)) :: xratio
        integer idx
        integer i, j, k

        xratio = mixlayer_vel_coef

        ! Find the most unstable eigenvalue and corresponding eigenvector
        k = 0
        do i = 1, mixlayer_nvar*n - n_bc_skip - 1
            if (wi(i) > wi(k)) then
                k = i
            end if
        end do
        vr = zr(:, k)
        vi = zi(:, k)

        ! Normalize the eigenvector by its component with the largest modulus.
        norm = 0d0
        do i = 0, mixlayer_nvar*n - n_bc_skip - 1
            if (dsqrt(vr(i)**2 + vi(i)**2) > norm) then
                idx = i
                norm = dsqrt(vr(i)**2 + vi(i)**2)
            end if
        end do

        tr = vr(idx)
        ti = vi(idx)
        do i = 0, mixlayer_nvar*n - n_bc_skip - 1
            call cdiv(vr(i), vi(i), tr, ti, cr, ci)
            vnr(i) = cr
            vni(i) = ci
        end do

        ! Reassign missing values at boundaries based on the boundary condition
        xbr = 0d0
        xbi = 0d0
        do i = 1, mixlayer_nvar
            if (i <= mixlayer_var(2)) then
                do k = 0, n - 1
                    xbr((i - 1)*nbp + k + 1) = vnr((i - 1)*n + k)
                    xbi((i - 1)*nbp + k + 1) = vni((i - 1)*n + k)
                end do
            elseif (i == mixlayer_var(3)) then
                do k = 0, n - n_bc_skip - 1
                    xbr((i - 1)*nbp + mixlayer_bc_fd + k + 1) = vnr((i - 1)*n + k)
                    xbi((i - 1)*nbp + mixlayer_bc_fd + k + 1) = vni((i - 1)*n + k)
                end do
            else
                do k = 0, n - 1
                    xbr((i - 1)*nbp + k + 1) = vnr((i - 1)*n - n_bc_skip + k)
                    xbi((i - 1)*nbp + k + 1) = vni((i - 1)*n - n_bc_skip + k)
                end do
            end if
        end do

        ! rho at boundaries
        xbr(0) = xbr(1) + xbr(mixlayer_var(2)*nbp + 1)*rho_mean*mach
        xbi(0) = xbi(1) + xbi(mixlayer_var(2)*nbp + 1)*rho_mean*mach
        xbr(nbp - 1) = xbr(n) - xbr(mixlayer_var(2)*nbp + n)*rho_mean*mach
        xbi(nbp - 1) = xbi(n) - xbi(mixlayer_var(2)*nbp + n)*rho_mean*mach

        ! u at boundaries
        xbr(mixlayer_var(1)*nbp) = xbr(mixlayer_var(1)*nbp + 1)
        xbi(mixlayer_var(1)*nbp) = xbi(mixlayer_var(1)*nbp + 1)
        xbr(mixlayer_var(1)*nbp + nbp - 1) = xbr(mixlayer_var(1)*nbp + n)
        xbi(mixlayer_var(1)*nbp + nbp - 1) = xbi(mixlayer_var(1)*nbp + n)

        ! w at boundaries
        xbr((mixlayer_var(3))*nbp + 0) = xbr((mixlayer_var(3))*nbp + 1)
        xbi((mixlayer_var(3))*nbp + 0) = xbi((mixlayer_var(3))*nbp + 1)
        xbr((mixlayer_var(3))*nbp + nbp - 1) = xbr((mixlayer_var(3))*nbp + n)
        xbi((mixlayer_var(3))*nbp + nbp - 1) = xbi((mixlayer_var(3))*nbp + n)

        ! p at boundaries
        xbr(mixlayer_var(4)*nbp + 0) = xbr(mixlayer_var(4)*nbp + 1) + xbr(mixlayer_var(2)*nbp + 1)*rho_mean/mach
        xbi(mixlayer_var(4)*nbp + 0) = xbi(mixlayer_var(4)*nbp + 1) + xbi(mixlayer_var(2)*nbp + 1)*rho_mean/mach
        xbr(mixlayer_var(4)*nbp + nbp - 1) = xbr(mixlayer_var(4)*nbp + n) - xbr(mixlayer_var(2)*nbp + n)*rho_mean/mach
        xbi(mixlayer_var(4)*nbp + nbp - 1) = xbi(mixlayer_var(4)*nbp + n) - xbi(mixlayer_var(2)*nbp + n)*rho_mean/mach

        ! Compute average to get cell-centered values
        xcr = 0d0
        xci = 0d0
        do i = 1, mixlayer_nvar
            do k = 0, n
                xcr((i - 1)*(nbp - 1) + k) = 5d-1*(xbr((i - 1)*nbp + k) + xbr((i - 1)*nbp + k + 1))
                xci((i - 1)*(nbp - 1) + k) = 5d-1*(xbi((i - 1)*nbp + k) + xbi((i - 1)*nbp + k + 1))
            end do
        end do

        ! Generate instability waves in x- and z-directions with phase shifts
        ! wave = Re(eigfunc * exp(i*(alpha*x + beta*z)))
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (beta == 0) then
                        ang = alpha*(x_cc(i)*xratio)
                    else
                        ang = alpha*(x_cc(i)*xratio) + beta*(z_cc(k)*xratio) + shift
                    end if
                    wave(mixlayer_var(1), i, j, k) = xcr(j)*cos(ang) - xci(j)*sin(ang) ! rho
                    wave(mixlayer_var(2), i, j, k) = xcr(mixlayer_var(1)*(nbp - 1) + j)*cos(ang) - xci(mixlayer_var(1)*(nbp - 1) + j)*sin(ang) ! u
                    wave(mixlayer_var(3), i, j, k) = xcr(mixlayer_var(2)*(nbp - 1) + j)*cos(ang) - xci(mixlayer_var(2)*(nbp - 1) + j)*sin(ang) ! v
                    wave(mixlayer_var(4), i, j, k) = xcr(mixlayer_var(3)*(nbp - 1) + j)*cos(ang) - xci(mixlayer_var(3)*(nbp - 1) + j)*sin(ang) ! w
                    wave(mixlayer_var(5), i, j, k) = xcr(mixlayer_var(4)*(nbp - 1) + j)*cos(ang) - xci(mixlayer_var(4)*(nbp - 1) + j)*sin(ang) ! p
                end do
            end do
        end do

    end subroutine s_generate_wave

end module m_perturbation
