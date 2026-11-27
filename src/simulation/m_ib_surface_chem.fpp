!>
!! @file
!! @brief Contains module m_ib_surface_chem

#:include 'macros.fpp'

!> @brief Reacting-surface state of an immersed boundary: the Newton solve for the surface composition and temperature, and the
!! bounded ghost-state blend. Takes no ghost-point arguments; m_ibm couples it to the ghost points.
module m_ib_surface_chem

    use m_precision_select, only: wp
    use m_constants, only: T_surface_min, T_surface_max
    use m_thermochem, only: num_species, gas_constant, get_mixture_molecular_weight, get_mixture_thermal_conductivity_mixavg, &
        & get_species_mass_diffusivities_mixavg
    use m_surface_thermochem, only: get_surface_net_production_rates, get_surface_reaction_heat_flux

    implicit none

    private; public :: s_solve_surface, s_blend_ghost_temperature, s_blend_ghost_species

    ! The species Jacobian is a forward difference, so the Newton update carries O(sqrt(epsilon)) noise and the iteration turns
    ! linear near the root. Cantera's solveSP runs the same 1e-7 forward-difference step against a 1e-4 relative tolerance;
    ! 1e-6 here stays well inside that. energy_tol is relative to the reaction heat max(|q|, 1).
    real(wp), parameter :: species_tol = 1.e-6_wp, energy_tol = 1.e-6_wp

contains

    !> Species flux residual at a heterogeneous reacting surface.
    subroutine s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_s, d
        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species), W_species(num_species)
        real(wp), intent(out) :: R_species(num_species), omega_s(num_species), mdot_s
        real(wp)              :: mw_IP, mw_s, rho_s, sum_BG
        real(wp)              :: Xs_IP(num_species), Xs_s(num_species)
        real(wp)              :: D_s(num_species), B_s(num_species), G_s(num_species)
        integer               :: k

        call get_mixture_molecular_weight(Ys_IP, mw_IP)
        call get_mixture_molecular_weight(Ys_s, mw_s)
        rho_s = pres*mw_s/(gas_constant*T_s)

        do k = 1, num_species
            Xs_IP(k) = Ys_IP(k)*mw_IP/W_species(k)
            Xs_s(k) = Ys_s(k)*mw_s/W_species(k)
        end do

        call get_species_mass_diffusivities_mixavg(pres, T_s, Ys_s, D_s)
        call get_surface_net_production_rates(rho_s, T_s, Ys_s, omega_s)

        mdot_s = 0._wp
        do k = 1, num_species
            mdot_s = mdot_s + W_species(k)*omega_s(k)
        end do

        sum_BG = 0._wp
        do k = 1, num_species
            B_s(k) = rho_s*D_s(k)*W_species(k)/mw_s
            G_s(k) = (Xs_IP(k) - Xs_s(k))/d
            sum_BG = sum_BG + B_s(k)*G_s(k)
        end do

        do k = 1, num_species
            R_species(k) = -B_s(k)*G_s(k) + Ys_s(k)*(sum_BG + mdot_s) - W_species(k)*omega_s(k)
        end do

    end subroutine s_surface_species_residual

    !> Surface energy residual: gas-side conduction balances heterogeneous reaction heat. Radiation and solid-side conduction are
    !! omitted.
    subroutine s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy, q_rxn)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_s, d
        real(wp), intent(in)  :: Ys_s(num_species)
        real(wp), intent(out) :: R_energy, q_rxn
        real(wp)              :: mw_s, rho_s, k_s

        call get_mixture_molecular_weight(Ys_s, mw_s)
        rho_s = pres*mw_s/(gas_constant*T_s)

        call get_mixture_thermal_conductivity_mixavg(T_s, Ys_s, k_s)
        call get_surface_reaction_heat_flux(rho_s, T_s, Ys_s, q_rxn)

        R_energy = k_s*(T_s - T_IP)/d - q_rxn

    end subroutine s_surface_energy_residual

    !> Scaled Newton residual of the species balances at fixed T_s. k_bath selects which species balance the sum constraint
    !! displaces; see s_pick_bath_species for why it must be the most abundant one.
    subroutine s_surface_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, k_bath, flux_scale, R, omega_s, mdot_s)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_s, d, flux_scale
        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species), W_species(num_species)
        integer, intent(in)   :: k_bath
        real(wp), intent(out) :: R(num_species), omega_s(num_species), mdot_s
        real(wp)              :: R_species(num_species)
        integer               :: k

        call s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)

        do k = 1, num_species
            R(k) = R_species(k)/flux_scale
        end do
        R(k_bath) = sum(Ys_s) - 1._wp

    end subroutine s_surface_residual

    !> Ghost state for an immersed surface, from the image-point state and the surface state.
    !!
    !! The natural closure is the linear mirror phi_g = 2*phi_s - phi_IP, which is what makes the midpoint of the ghost/image
    !! pair reproduce the surface value. Its -1 coefficient on phi_IP is unconditional, though, so it drives strictly positive
    !! quantities negative whenever the surface value sits far below the image-point value: a cold wall in hot gas gives a
    !! negative ghost temperature and hence a negative ghost density, and a species the surface consumes faster than half the
    !! free-stream value gives a negative mass fraction. Rather than clamp the result -- which would break sum(Y) = 1 and hide
    !! the excursion -- the ghost state is the largest convex blend of the mirror back toward the surface value that stays
    !! physical:
    !!
    !!     phi_g(theta) = phi_s + theta*(phi_s - phi_IP),   theta in [0, 1]
    !!
    !! theta = 1 is the full second-order mirror; theta = 0 is the first-order Dirichlet ghost phi_g = phi_s, which is
    !! Gibou et al. (JCP 176:205, 2002) Eq. 17 and is what they likewise fall back to where the linear form is ill-behaved.
    !! Each species and the temperature get their own theta (s_blend_ghost_species keeps sum(Y) = 1 through the bath species).
    !! Sharing one would let a trace radical the surface consumes to ~1e-9 (O and OH at a burning carbon wall) drag every other
    !! mirror down with it, turning a 2100 K ghost into 1205 K over a species whose own excursion was 1e-9. The cost of limiting
    !! is small: the mirror already delivers only first-order wall flux (Ezra et al., Int. J. Heat Mass Transfer, 2025),
    !! and a boundary closure one order below the interior scheme retains the interior convergence rate (Gustafsson,
    !! Math. Comp. 29:396, 1975).
    !!
    !! This bounds the ghost state; it does not make the boundary strictly conservative. Baskaya et al. (Computers & Fluids
    !! 270:106134, 2024) trace the same cold-wall, large-gradient regime in ablation to ghost-cell mass conservation error
    !! that surfaces as spurious blowing, and resolve it only by moving to a flux-based cut-cell boundary -- a different
    !! discretization from this one, not a tuning of it.
    subroutine s_blend_ghost_temperature(T_IP, T_s, T_lo, T_hi, theta_max, T_g)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: T_IP, T_s, T_lo, T_hi
        real(wp), intent(in)  :: theta_max  !< Largest mirror weight the caller allows (f_face_theta), at most 1
        real(wp), intent(out) :: T_g
        ! Stop short of the bound so roundoff cannot cross it.
        real(wp), parameter :: blend_safety = 0.9_wp
        real(wp)            :: theta_T

        ! Largest theta <= theta_max keeping T_g in [T_lo, T_hi]. The bounds are compared before dividing, so +-huge() cannot
        ! overflow.
        theta_T = theta_max
        if (T_s < T_IP .and. blend_safety*(T_s - T_lo) < theta_T*(T_IP - T_s)) then
            theta_T = blend_safety*(T_s - T_lo)/(T_IP - T_s)
        else if (T_s > T_IP .and. blend_safety*(T_hi - T_s) < theta_T*(T_s - T_IP)) then
            theta_T = blend_safety*(T_hi - T_s)/(T_s - T_IP)
        end if
        theta_T = max(theta_T, 0._wp)

        T_g = T_s + theta_T*(T_s - T_IP)

    end subroutine s_blend_ghost_temperature

    !> Species half of the ghost-state blend: each species takes its own positivity-limited theta and the bath species (the most
    !! abundant at the surface) absorbs the difference, so sum(Y_g) = sum(Y_s). One theta shared by every species collapsed to
    !! ~0.006 at a burning carbon wall, where the surface consumes the trace radicals (O, OH), and made every species ghost first
    !! order. The shared theta is the fallback if the bath species would go negative.
    subroutine s_blend_ghost_species(Ys_IP, Ys_s, Ys_g)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species)
        real(wp), intent(out) :: Ys_g(num_species)
        real(wp), parameter   :: blend_safety = 0.9_wp
        real(wp)              :: theta_Y
        integer               :: k, k_bath

        call s_pick_bath_species(Ys_s, k_bath)
        Ys_g(k_bath) = Ys_s(k_bath)
        do k = 1, num_species
            if (k == k_bath) cycle
            theta_Y = 1._wp
            if (Ys_s(k) < Ys_IP(k)) theta_Y = min(1._wp, blend_safety*Ys_s(k)/(Ys_IP(k) - Ys_s(k)))
            Ys_g(k) = Ys_s(k) + theta_Y*(Ys_s(k) - Ys_IP(k))
            Ys_g(k_bath) = Ys_g(k_bath) + Ys_s(k) - Ys_g(k)
        end do
        if (Ys_g(k_bath) >= 0._wp) return

        theta_Y = 1._wp
        do k = 1, num_species
            if (Ys_s(k) < Ys_IP(k)) theta_Y = min(theta_Y, blend_safety*Ys_s(k)/(Ys_IP(k) - Ys_s(k)))
        end do
        do k = 1, num_species
            Ys_g(k) = Ys_s(k) + theta_Y*(Ys_s(k) - Ys_IP(k))
        end do

    end subroutine s_blend_ghost_species

    !> Index of the most abundant species, i.e. the balance the sum constraint displaces.
    !!
    !! Closing sum(Y) = 1 by dropping one species' flux balance is standard (Surface CHEMKIN), but the dropped species absorbs
    !! the roundoff of every other balance, so it must be the bath gas. Dropping a fixed index instead -- the last species in
    !! the mechanism -- lands that error on whatever the mechanism happens to list last, which for the reduced GRI mechanism
    !! shipped with the reacting-surface example is H2O2, a trace radical whose own surface balance then goes unenforced.
    !! Cantera's solveSP re-scans for the largest species each iteration (evalSurfLarge) for this reason.
    subroutine s_pick_bath_species(Ys_s, k_bath)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: Ys_s(num_species)
        integer, intent(out) :: k_bath
        integer              :: k

        k_bath = 1
        do k = 2, num_species
            if (Ys_s(k) > Ys_s(k_bath)) k_bath = k
        end do

    end subroutine s_pick_bath_species

    !> Solve for a reacting surface. thermal_bc=0: zero-normal-gradient T; 1: prescribed T; 2: energy balance.
    subroutine s_solve_surface(pres, T_IP, T_wall, d, Ys_IP, W_species, thermal_bc, Ys_s, T_s, mdot_s, converged)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_wall, d
        real(wp), intent(in)  :: Ys_IP(num_species), W_species(num_species)
        integer, intent(in)   :: thermal_bc
        real(wp), intent(out) :: Ys_s(num_species), T_s, mdot_s
        logical, intent(out)  :: converged
        real(wp)              :: omega_s(num_species)

        converged = .false.
        Ys_s = Ys_IP
        mdot_s = 0._wp

        select case (thermal_bc)
        case (0)
            T_s = T_IP
        case (1)
            T_s = T_wall
        case (2)
            call s_solve_surface_energy(pres, T_IP, d, Ys_IP, W_species, Ys_s, T_s, mdot_s, converged)
            return
        case default
            T_s = T_IP
            return
        end select

        call s_solve_surface_species(pres, T_s, d, Ys_IP, W_species, species_tol, Ys_s, omega_s, mdot_s, converged)

    end subroutine s_solve_surface

    !> Damped Newton solve of the species flux balances at a fixed surface temperature, from the initial guess in Ys_s. The residual
    !! is scaled by its value at the image-point composition, so a warm start does not tighten the tolerance. It iterates to tol,
    !! and a stall (no descent step, singular Jacobian or max_iter) still counts as converged below species_tol.
    subroutine s_solve_surface_species(pres, T_s, d, Ys_IP, W_species, tol, Ys_s, omega_s, mdot_s, converged)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)    :: pres, T_s, d, tol
        real(wp), intent(in)    :: Ys_IP(num_species), W_species(num_species)
        real(wp), intent(inout) :: Ys_s(num_species)
        real(wp), intent(out)   :: omega_s(num_species), mdot_s
        logical, intent(out)    :: converged
        integer, parameter      :: max_iter = 30, max_backtrack = 20
        real(wp), parameter     :: fd_eps_Y = 1.e-7_wp
        real(wp), parameter     :: Y_tol = 100._wp*epsilon(1._wp)
        real(wp)                :: A(num_species, num_species), rhs(num_species), delta(num_species)
        real(wp)                :: R(num_species), R_pert(num_species), R_trial(num_species), R_species(num_species)
        real(wp)                :: omega_pert(num_species), omega_trial(num_species)
        real(wp)                :: Ys_pert(num_species), Ys_trial(num_species)
        real(wp)                :: mdot_pert, mdot_trial, flux_scale, dx, lambda, norm_R, norm_trial
        logical                 :: linear_success, accepted
        integer                 :: iter, j, iback, k_bath, k_bath_prev

        converged = .false.

        call s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_IP, W_species, R_species, omega_s, mdot_s)
        flux_scale = max(maxval(abs(R_species)), 1.e-12_wp)

        call s_pick_bath_species(Ys_s, k_bath)
        call s_surface_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, k_bath, flux_scale, R, omega_s, mdot_s)
        norm_R = maxval(abs(R))
        if (norm_R < tol) then
            converged = .true.
            return
        end if
        ! The rates are not evaluated outside the surface temperature range
        if (T_s < T_surface_min .or. T_s > T_surface_max) return

        do iter = 1, max_iter
            do j = 1, num_species
                Ys_pert = Ys_s
                dx = fd_eps_Y*max(abs(Ys_s(j)), 1._wp)
                if (Ys_s(j) + dx > 1._wp) dx = -dx
                Ys_pert(j) = Ys_pert(j) + dx

                call s_surface_residual(pres, T_s, d, Ys_IP, Ys_pert, W_species, k_bath, flux_scale, R_pert, omega_pert, mdot_pert)
                A(:,j) = (R_pert - R)/dx
            end do

            rhs = -R
            call s_solve_surface_linear_system(A, rhs, delta, linear_success)
            if (.not. linear_success) exit

            lambda = 1._wp
            accepted = .false.
            do iback = 1, max_backtrack
                Ys_trial = Ys_s + lambda*delta

                if (minval(Ys_trial) < -Y_tol .or. maxval(Ys_trial) > 1._wp + Y_tol) then
                    lambda = 0.5_wp*lambda
                    cycle
                end if

                where (Ys_trial < 0._wp) Ys_trial = 0._wp
                where (Ys_trial > 1._wp) Ys_trial = 1._wp

                call s_surface_residual(pres, T_s, d, Ys_IP, Ys_trial, W_species, k_bath, flux_scale, R_trial, omega_trial, &
                                        & mdot_trial)
                norm_trial = maxval(abs(R_trial))
                if (norm_trial < norm_R) then
                    accepted = .true.
                    exit
                end if
                lambda = 0.5_wp*lambda
            end do

            if (.not. accepted) exit

            Ys_s = Ys_trial
            R = R_trial
            omega_s = omega_trial
            mdot_s = mdot_trial
            norm_R = norm_trial

            ! The composition has moved, so the species the sum constraint displaces may no longer be the most abundant one.
            ! Re-closing changes what R means, so the residual is rebuilt on the rare steps where the choice actually changes.
            k_bath_prev = k_bath
            call s_pick_bath_species(Ys_s, k_bath)
            if (k_bath /= k_bath_prev) then
                call s_surface_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, k_bath, flux_scale, R, omega_s, mdot_s)
                norm_R = maxval(abs(R))
            end if

            if (norm_R < tol) then
                converged = .true.
                return
            end if
        end do

        converged = norm_R < species_tol

    end subroutine s_solve_surface_species

    !> Energy balance F(T_s) = k (T_s - T_IP)/d - q(T_s) = 0 for thermal_bc = 2, with the species balance solved at each T_s.
    !!
    !! F(T_s) can have up to three roots. For burning carbon these are the kinetic root near T_IP, an unstable middle root, and
    !! the diffusion-limited root far above the gas. The root returned is the first sign change of F met when walking from T_IP
    !! in the direction of -F(T_IP). It is the steady state that a zero-heat-capacity wall starting at the gas temperature
    !! relaxes to (dT_s/dt ~ -F), so it is always a stable root. It is the kinetic root while one exists, and the
    !! diffusion-limited root once the surface has ignited. It needs no state from earlier steps. Surface thermal history and
    !! extinction hysteresis are thermal_bc = 3's job. A coupled Newton in (Y, T_s) cannot be used here: after ignition
    !! dF/dT_s < 0 at T_IP, so its T_s step points away from the root, and a line search on max|R| rejects every step toward
    !! it, because |F| has to grow before it falls.
    !!
    !! The walk takes steps of at most step_frac*T_s, shortened to twice the secant estimate of the root so that the steps
    !! close in on it. Once F changes sign, Illinois regula falsi finishes, with a bisection step whenever two steps in a row
    !! fail to halve the bracket. A kinetic root closer than one step to its unstable partner, that is within ~5% of the
    !! ignition turning point, can be stepped over. No root before T_surface_min or T_surface_max is reported as not
    !! converged.
    subroutine s_solve_surface_energy(pres, T_IP, d, Ys_IP, W_species, Ys_s, T_s, mdot_s, converged)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, d
        real(wp), intent(in)  :: Ys_IP(num_species), W_species(num_species)
        real(wp), intent(out) :: Ys_s(num_species), T_s, mdot_s
        logical, intent(out)  :: converged
        integer, parameter    :: max_walk = 100, max_bracket = 60
        real(wp), parameter   :: step_frac = 0.05_wp
        real(wp)              :: Ys_b(num_species), Ys_c(num_species)
        real(wp)              :: T_a, F_a, T_b, F_b, T_c, F_c, T_p, F_p, mdot_b, mdot_c, scale_c, dir, step, dT_secant, width
        logical               :: ok, bracketed
        integer               :: n, nslow

        converged = .false.

        ! Start: T_a = T_IP, clamped into [T_surface_min, T_surface_max]
        T_a = min(max(T_IP, T_surface_min), T_surface_max)
        Ys_s = Ys_IP
        call s_surface_energy_state(pres, T_IP, T_a, d, Ys_IP, W_species, Ys_s, F_a, scale_c, mdot_s, ok)
        T_s = T_a
        if (.not. ok) return
        if (abs(F_a) < energy_tol*scale_c) then
            converged = .true.
            return
        end if

        ! Walk from T_IP until F changes sign
        dir = merge(1._wp, -1._wp, F_a < 0._wp)
        T_p = T_a
        F_p = F_a
        Ys_b = Ys_s
        bracketed = .false.
        do n = 1, max_walk
            step = step_frac*T_a
            if (F_a /= F_p) then
                dT_secant = -F_a*(T_a - T_p)/(F_a - F_p)
                if (dT_secant*dir > 0._wp) step = min(step, 2._wp*abs(dT_secant))
            end if
            T_b = min(max(T_a + dir*step, T_surface_min), T_surface_max)
            call s_surface_energy_state(pres, T_IP, T_b, d, Ys_IP, W_species, Ys_b, F_b, scale_c, mdot_b, ok)
            if (.not. ok) return
            if (abs(F_b) < energy_tol*scale_c) then
                call s_accept_surface(T_b, Ys_b, mdot_b, T_s, Ys_s, mdot_s, converged)
                return
            end if
            if ((F_b < 0._wp) .neqv. (F_a < 0._wp)) then
                bracketed = .true.
                exit
            end if
            if (T_b == T_surface_min .or. T_b == T_surface_max) return
            T_p = T_a
            F_p = F_a
            T_a = T_b
            F_a = F_b
        end do
        if (.not. bracketed) return

        ! Illinois regula falsi on [T_a, T_b], whose ends keep opposite signs
        nslow = 0
        width = abs(T_b - T_a)
        Ys_c = Ys_b
        do n = 1, max_bracket
            if (nslow < 2) then
                T_c = T_b - F_b*(T_b - T_a)/(F_b - F_a)
            else
                T_c = 0.5_wp*(T_a + T_b)
            end if
            call s_surface_energy_state(pres, T_IP, T_c, d, Ys_IP, W_species, Ys_c, F_c, scale_c, mdot_c, ok)
            if (.not. ok) return
            if (abs(F_c) < energy_tol*scale_c) then
                call s_accept_surface(T_c, Ys_c, mdot_c, T_s, Ys_s, mdot_s, converged)
                return
            end if
            if ((F_c < 0._wp) .eqv. (F_b < 0._wp)) then
                F_a = 0.5_wp*F_a
            else
                T_a = T_b
                F_a = F_b
            end if
            T_b = T_c
            F_b = F_c
            nslow = merge(nslow + 1, 0, abs(T_b - T_a) > 0.5_wp*width)
            width = abs(T_b - T_a)
            if (width <= 100._wp*epsilon(1._wp)*T_b) then
                call s_accept_surface(T_c, Ys_c, mdot_c, T_s, Ys_s, mdot_s, converged)
                return
            end if
        end do

    end subroutine s_solve_surface_energy

    !> Species solve at T_s (warm-started from Ys_s) and the energy residual F there, with its scale max(|q|, 1).
    subroutine s_surface_energy_state(pres, T_IP, T_s, d, Ys_IP, W_species, Ys_s, F, F_scale, mdot_s, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)    :: pres, T_IP, T_s, d
        real(wp), intent(in)    :: Ys_IP(num_species), W_species(num_species)
        real(wp), intent(inout) :: Ys_s(num_species)
        real(wp), intent(out)   :: F, F_scale, mdot_s
        logical, intent(out)    :: ok
        ! F inherits the species solve's error as noise of about tol*flux_scale/(diffusive flux). At a hot surface the
        ! image-point scale is ~1e4 times the diffusive flux, so species_tol would leave F with 1e-5 noise and stall the root
        ! solve. The species residual is cheap to drive further, to where its own roundoff stops it.
        real(wp), parameter :: nested_tol = 1.e-12_wp
        real(wp)            :: omega_s(num_species), q_rxn

        call s_solve_surface_species(pres, T_s, d, Ys_IP, W_species, nested_tol, Ys_s, omega_s, mdot_s, ok)
        F = 0._wp
        F_scale = 1._wp
        if (.not. ok) return
        call s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, F, q_rxn)
        F_scale = max(abs(q_rxn), 1._wp)

    end subroutine s_surface_energy_state

    !> Copy an accepted energy-balance state to the outputs.
    subroutine s_accept_surface(T, Ys, mdot, T_s, Ys_s, mdot_s, converged)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: T, Ys(num_species), mdot
        real(wp), intent(out) :: T_s, Ys_s(num_species), mdot_s
        logical, intent(out)  :: converged

        T_s = T
        Ys_s = Ys
        mdot_s = mdot
        converged = .true.

    end subroutine s_accept_surface

    !> Small dense linear solve with partial pivoting for the local surface Newton system.
    subroutine s_solve_surface_linear_system(A, b, x, success)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(inout) :: A(num_species, num_species)
        real(wp), intent(inout) :: b(num_species)
        real(wp), intent(out)   :: x(num_species)
        logical, intent(out)    :: success
        real(wp)                :: factor, pivot_value, tmp, row_tmp(num_species)
        integer                 :: i, j, k, pivot

        success = .true.
        x = 0._wp

        do k = 1, num_species - 1
            pivot = k
            pivot_value = abs(A(k, k))
            do i = k + 1, num_species
                if (abs(A(i, k)) > pivot_value) then
                    pivot = i
                    pivot_value = abs(A(i, k))
                end if
            end do
            if (pivot_value <= epsilon(1._wp)) then
                success = .false.
                return
            end if

            if (pivot /= k) then
                row_tmp(1:num_species) = A(k,1:num_species)
                A(k,1:num_species) = A(pivot,1:num_species)
                A(pivot,1:num_species) = row_tmp(1:num_species)
                tmp = b(k)
                b(k) = b(pivot)
                b(pivot) = tmp
            end if

            do i = k + 1, num_species
                factor = A(i, k)/A(k, k)
                A(i, k) = 0._wp
                do j = k + 1, num_species
                    A(i, j) = A(i, j) - factor*A(k, j)
                end do
                b(i) = b(i) - factor*b(k)
            end do
        end do

        if (abs(A(num_species, num_species)) <= epsilon(1._wp)) then
            success = .false.
            return
        end if

        x(num_species) = b(num_species)/A(num_species, num_species)
        do i = num_species - 1, 1, -1
            tmp = b(i)
            do j = i + 1, num_species
                tmp = tmp - A(i, j)*x(j)
            end do
            if (abs(A(i, i)) <= epsilon(1._wp)) then
                success = .false.
                return
            end if
            x(i) = tmp/A(i, i)
        end do

    end subroutine s_solve_surface_linear_system

end module m_ib_surface_chem
