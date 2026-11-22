!>
!! @file
!! @brief Contains module m_ib_wall_column

#:include 'macros.fpp'

!> @brief Sub-grid wall conduction column of an isothermal immersed wall (patch_ib%wall_column): state and 1D numerics.
!!
!! A linear ghost mirror carries at most about 2 k dT/h, while a wall layer thinner than a cell carries k dT/(sqrt(pi) delta):
!! a cold particle in hot gas takes up too little heat early on. Each first-layer ghost of such a wall carries a 1D low-Mach
!! conduction column along the wall normal, from the wall (T_w) to s_m = max(3h, L) (T from the gas), with the body's
!! curvature. Its step-mean wall flux reaches the gas as one equivalent conduction temperature per ghost (m_ibm).
module m_ib_wall_column

    use m_derived_types
    use m_global_parameters

    implicit none

    private
    public :: wc_n, wc_ntab, wc_guard, wc_len_frac, wc_T, wc_w, wc_g, wc_ext, wc_q, wc_cut, wc_on, wc_gp, wc_slot, n_wc, &
        & s_wc_step, f_wc_ell, f_wc_node, f_wc_target

    integer, parameter  :: wc_n = 96              !< Column cells
    integer, parameter  :: wc_nsub = 16           !< SDIRK2 substeps per flow step
    integer, parameter  :: wc_ntab = 8            !< Temperatures of the per-step property tables
    real(wp), parameter :: wc_ratio = 1.08_wp     !< Cell stretching from the wall
    real(wp), parameter :: wc_guard = 4._wp       !< Strength of the one-sided guard
    real(wp), parameter :: wc_len_frac = 0.15_wp  !< Default column length / radius

    !> Per column ghost (the first-layer ghosts of wall_column bodies, n_wc of them; wc_gp(c) is the ghost's index in ghost_points,
    !! wc_slot its inverse during setup): column temperatures at the cell centres (wc_T), step-mean wall heat flux into the body per
    !! unit wall area (wc_q), guard heat booked to the body (wc_cut) and state (wc_on: -1 none, 0 the mirror this step, 1 active).
    !! Per face (2 per direction): area times l-coordinate weight (wc_w), area/spacing (wc_g, 0 off a gas face) and the column
    !! temperature averaged over the gas cell (wc_ext).
    real(wp), allocatable, dimension(:,:) :: wc_T, wc_w, wc_g, wc_ext
    real(wp), allocatable, dimension(:)   :: wc_q, wc_cut
    integer, allocatable, dimension(:)    :: wc_on, wc_gp, wc_slot
    integer                               :: n_wc = 0
    $:GPU_DECLARE(create='[wc_T, wc_w, wc_g, wc_ext, wc_q, wc_cut, wc_on, wc_gp, wc_slot]')

contains

    !> Wall distance of column node k: 0 the wall, 1..wc_n the cell centres, wc_n + 1 the top s_m. Faces are geometric from the
    !! wall, f_k = s_m (r^k - 1)/(r^N - 1).
    function f_wc_node(k, sm) result(s)

        $:GPU_ROUTINE(function_name='f_wc_node', parallelism='[seq]', cray_inline=True)
        integer, intent(in)  :: k
        real(wp), intent(in) :: sm
        real(wp)             :: s

        if (k <= 0) then
            s = 0._wp
        else if (k > wc_n) then
            s = sm
        else
            s = 0.5_wp*(f_wc_face(k - 1, sm) + f_wc_face(k, sm))
        end if

    end function f_wc_node

    function f_wc_face(k, sm) result(s)

        $:GPU_ROUTINE(function_name='f_wc_face', parallelism='[seq]', cray_inline=True)
        integer, intent(in)  :: k
        real(wp), intent(in) :: sm
        real(wp)             :: s

        s = sm*(wc_ratio**k - 1._wp)/(wc_ratio**wc_n - 1._wp)

    end function f_wc_face

    !> Area of the surface at wall distance s per unit wall area, (1 + kappa s)^(d - 1)
    function f_wc_area(kap, s) result(a)

        $:GPU_ROUTINE(function_name='f_wc_area', parallelism='[seq]', cray_inline=True)
        real(wp), intent(in) :: kap, s
        real(wp)             :: a

        a = (1._wp + kap*s)**(num_dims - 1)

    end function f_wc_area

    !> Steady radial conduction coordinate, int ds/A: ln(1 + kappa s)/kappa (2D), s/(1 + kappa s) (3D)
    function f_wc_ell(kap, s) result(l)

        $:GPU_ROUTINE(function_name='f_wc_ell', parallelism='[seq]', cray_inline=True)
        real(wp), intent(in) :: kap, s
        real(wp)             :: l

        if (num_dims == 2) then
            l = log(1._wp + kap*s)/kap
        else
            l = s/(1._wp + kap*s)
        end if

    end function f_wc_ell

    !> Column temperature at wall distance s >= 0 (linear between nodes, the last cell's value beyond it)
    function f_wc_interp(Tc, Tw, sm, s) result(T)

        $:GPU_ROUTINE(function_name='f_wc_interp', parallelism='[seq]', cray_inline=True)
        real(wp), dimension(wc_n), intent(in) :: Tc
        real(wp), intent(in)                  :: Tw, sm, s
        real(wp)                              :: T, sc, t0, x0, x1
        integer                               :: k

        sc = min(max(s, 0._wp), f_wc_node(wc_n, sm))
        k = int(log(1._wp + sc*(wc_ratio**wc_n - 1._wp)/sm)/log(wc_ratio)) + 1  ! cell holding sc
        k = min(max(k, 1), wc_n)
        if (sc >= f_wc_node(k, sm)) k = k + 1
        k = min(max(k - 1, 0), wc_n - 1)  ! lower node
        x0 = f_wc_node(k, sm); x1 = f_wc_node(k + 1, sm)
        t0 = Tw
        if (k > 0) t0 = Tc(k)
        T = t0 + (Tc(k + 1) - t0)*(sc - x0)/(x1 - x0)

    end function f_wc_interp

    !> Guard target at a point of signed wall distance s: the column value, and behind the wall its mirror-odd extension 2 T_w -
    !! T(-s), exact for a linear profile, so a smooth steady cell's average is its centre value
    function f_wc_target(Tc, Tw, sm, s) result(T)

        $:GPU_ROUTINE(function_name='f_wc_target', parallelism='[seq]', cray_inline=True)
        real(wp), dimension(wc_n), intent(in) :: Tc
        real(wp), intent(in)                  :: Tw, sm, s
        real(wp)                              :: T

        T = f_wc_interp(Tc, Tw, sm, abs(s))
        if (s < 0._wp) T = 2._wp*Tw - T

    end function f_wc_target

    !> Linear interpolation in a table on wc_ntab uniform temperatures from Tlo, constant beyond its ends
    function f_wc_tab(tab, T, Tlo, dTs) result(v)

        $:GPU_ROUTINE(function_name='f_wc_tab', parallelism='[seq]', cray_inline=True)
        real(wp), dimension(wc_ntab), intent(in) :: tab
        real(wp), intent(in)                     :: T, Tlo, dTs
        real(wp)                                 :: v, x
        integer                                  :: j

        x = min(max((T - Tlo)/dTs, 0._wp), real(wc_ntab - 1, wp))
        j = min(int(x), wc_ntab - 2)
        v = tab(j + 1) + (x - real(j, wp))*(tab(j + 2) - tab(j + 1))

    end function f_wc_tab

    !> Phi = int lambda dT of the piecewise-linear table, lambda constant beyond its ends
    function f_wc_phi(T, Tlo, dTs, lam_t, phi_t) result(phi)

        $:GPU_ROUTINE(function_name='f_wc_phi', parallelism='[seq]', cray_inline=True)
        real(wp), intent(in)                     :: T, Tlo, dTs
        real(wp), dimension(wc_ntab), intent(in) :: lam_t, phi_t
        real(wp)                                 :: phi, x, y
        integer                                  :: j

        x = min(max((T - Tlo)/dTs, 0._wp), real(wc_ntab - 1, wp))
        j = min(int(x), wc_ntab - 2)
        y = (x - real(j, wp))*dTs
        phi = phi_t(j + 1) + lam_t(j + 1)*y + 0.5_wp*(lam_t(j + 2) - lam_t(j + 1))/dTs*y*y
        phi = phi + f_wc_tab(lam_t, T, Tlo, dTs)*(T - Tlo - x*dTs)

    end function f_wc_phi

    !> Kirchhoff mean conductivity (Phi(Tb) - Phi(Ta))/(Tb - Ta)
    function f_wc_kbar(Ta, Tb, Tlo, dTs, lam_t, phi_t) result(k)

        $:GPU_ROUTINE(function_name='f_wc_kbar', parallelism='[seq]', cray_inline=True)
        real(wp), intent(in)                     :: Ta, Tb, Tlo, dTs
        real(wp), dimension(wc_ntab), intent(in) :: lam_t, phi_t
        real(wp)                                 :: k

        if (abs(Tb - Ta) > 1.e-9_wp*(abs(Ta) + abs(Tb))) then
            k = (f_wc_phi(Tb, Tlo, dTs, lam_t, phi_t) - f_wc_phi(Ta, Tlo, dTs, lam_t, phi_t))/(Tb - Ta)
        else
            k = f_wc_tab(lam_t, 0.5_wp*(Ta + Tb), Tlo, dTs)
        end if

    end function f_wc_kbar

    !> Solve (rcv - c L) x = d in place, L tridiagonal with sub lo, diagonal -(lo + up) and super up (Thomas)
    subroutine s_wc_solve(lo, up, rcv, c, d)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(wc_n), intent(in)    :: lo, up, rcv
        real(wp), intent(in)                     :: c
        real(wp), dimension(wc_n), intent(inout) :: d
        real(wp), dimension(wc_n)                :: w
        real(wp)                                 :: den
        integer                                  :: k

        den = rcv(1) + c*(lo(1) + up(1))
        w(1) = -c*up(1)/den
        d(1) = d(1)/den
        do k = 2, wc_n
            den = rcv(k) + c*(lo(k) + up(k)) + c*lo(k)*w(k - 1)
            w(k) = -c*up(k)/den
            d(k) = (d(k) + c*lo(k)*d(k - 1))/den
        end do
        do k = wc_n - 1, 1, -1
            d(k) = d(k) - w(k)*d(k + 1)
        end do

    end subroutine s_wc_solve

    !> Advance one column over dt and return its step-mean wall heat flux into the body q, per unit wall area. Low-Mach energy
    !! equation of gas at rest at an isothermal wall, frozen composition, rho = rcoef/T (rcoef = p W/R_u):
    !!
    !!     rho c_p (dT/dt + u dT/ds) = (1/A) d/ds (A lambda dT/ds),   d(A u)/ds = d(A lambda dT/ds)/ds / (rcoef c_p),   u(0) = 0
    !!
    !! with T(0) = Tw and T(sm) = Tm. Finite volumes on the stretched cells; advection central where the cell Peclet number is
    !! below 1, else upwind. In time, wc_nsub substeps of the L-stable, stiffly accurate 2-stage SDIRK (gamma = 1 - 1/sqrt 2),
    !! coefficients lagged per substep, u from a diffusion-only implicit predictor. Both stages are implicit, so a step start does
    !! not ring into the wall flux. lam_t and cp_t tabulate lambda(T) and c_p(T) from Tlo in steps dTs.
    subroutine s_wc_step(Tc, Tw, Tm, rcoef, kap, sm, dt, Tlo, dTs, lam_t, cp_t, q)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(wc_n), intent(inout) :: Tc
        real(wp), intent(in)                     :: Tw, Tm, rcoef, kap, sm, dt, Tlo, dTs
        real(wp), dimension(wc_ntab), intent(in) :: lam_t, cp_t
        real(wp), intent(out)                    :: q
        real(wp), dimension(wc_ntab)             :: phi_t
        real(wp), dimension(0:wc_n)              :: a
        real(wp), dimension(wc_n)                :: rcv, lo, up, Y1, Y3
        real(wp)                                 :: h, gam, TL, TR, D0, D1, Au0, Au1, u, dl, dr, Ac, V, ep, F
        integer                                  :: it, k

        phi_t(1) = 0._wp
        do k = 2, wc_ntab
            phi_t(k) = phi_t(k - 1) + 0.5_wp*(lam_t(k - 1) + lam_t(k))*dTs
        end do

        h = dt/real(wc_nsub, wp)
        gam = 1._wp - 1._wp/sqrt(2._wp)
        q = 0._wp
        do it = 1, wc_nsub
            ! Face conductances A lambda/ds, wall to top
            do k = 0, wc_n
                TL = Tw
                if (k > 0) TL = Tc(k)
                TR = Tm
                if (k < wc_n) TR = Tc(k + 1)
                a(k) = f_wc_area(kap, f_wc_face(k, sm))*f_wc_kbar(TL, TR, Tlo, dTs, lam_t, phi_t)/(f_wc_node(k + 1, &
                  & sm) - f_wc_node(k, sm))
            end do
            do k = 1, wc_n
                rcv(k) = rcoef/Tc(k)*f_wc_tab(cp_t, Tc(k), Tlo, dTs)*(f_wc_face(k, sm) - f_wc_face(k - 1, sm))*f_wc_area(kap, &
                    & f_wc_node(k, sm))
                lo(k) = a(k - 1)
                up(k) = a(k)
                Y3(k) = rcv(k)*Tc(k)
            end do

            ! Diffusion-only predictor (lagging u from the substep start lets a step start's singular wall flux sweep the column)
            Y3(1) = Y3(1) + h*a(0)*Tw
            Y3(wc_n) = Y3(wc_n) + h*a(wc_n)*Tm
            call s_wc_solve(lo, up, rcv, h, Y3)

            ! Contraction velocity A u from the predicted conduction, and the advective couplings
            D0 = a(0)*(Y3(1) - Tw)
            Au0 = 0._wp
            do k = 1, wc_n
                TR = Tm
                if (k < wc_n) TR = Y3(k + 1)
                D1 = a(k)*(TR - Y3(k))
                Au1 = Au0 + (D1 - D0)/(rcoef*f_wc_tab(cp_t, Tc(k), Tlo, dTs))
                Ac = f_wc_area(kap, f_wc_node(k, sm))
                u = 0.5_wp*(Au0 + Au1)/Ac
                dl = f_wc_node(k, sm) - f_wc_node(k - 1, sm)
                dr = f_wc_node(k + 1, sm) - f_wc_node(k, sm)
                V = (f_wc_face(k, sm) - f_wc_face(k - 1, sm))*Ac
                if (abs(u)*rcv(k)*Ac*(dl + dr) < V*(a(k - 1) + a(k))) then
                    ep = -rcv(k)*u/(dl + dr)
                    lo(k) = a(k - 1) - ep
                else
                    ep = rcv(k)*max(-u, 0._wp)/dr
                    lo(k) = a(k - 1) + rcv(k)*max(u, 0._wp)/dl
                end if
                up(k) = a(k) + ep
                D0 = D1
                Au0 = Au1
            end do

            ! rcv dY/dt = L Y + b, b the wall and top values: two implicit stages
            do k = 1, wc_n
                Y1(k) = rcv(k)*Tc(k)
            end do
            Y1(1) = Y1(1) + gam*h*lo(1)*Tw
            Y1(wc_n) = Y1(wc_n) + gam*h*up(wc_n)*Tm
            call s_wc_solve(lo, up, rcv, gam*h, Y1)
            do k = 1, wc_n
                F = -(lo(k) + up(k))*Y1(k)
                if (k > 1) then
                    F = F + lo(k)*Y1(k - 1)
                else
                    F = F + lo(1)*Tw
                end if
                if (k < wc_n) then
                    F = F + up(k)*Y1(k + 1)
                else
                    F = F + up(wc_n)*Tm
                end if
                Y3(k) = rcv(k)*Tc(k) + (1._wp - gam)*h*F
            end do
            Y3(1) = Y3(1) + gam*h*lo(1)*Tw
            Y3(wc_n) = Y3(wc_n) + gam*h*up(wc_n)*Tm
            call s_wc_solve(lo, up, rcv, gam*h, Y3)

            q = q + a(0)*((1._wp - gam)*Y1(1) + gam*Y3(1) - Tw)
            do k = 1, wc_n
                Tc(k) = Y3(k)
            end do
        end do
        q = q/real(wc_nsub, wp)

    end subroutine s_wc_step

end module m_ib_wall_column
