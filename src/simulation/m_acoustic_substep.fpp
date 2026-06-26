!>
!! @file m_acoustic_substep.fpp
!! @brief Contains module m_acoustic_substep

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Cheap forward-backward acoustic substep kernel for the split-explicit low-Mach time integrator (acoustic_substepping).
!!
!! In split mode the HLLC solver returns only the slow (advective) flux:
!! contact-upwinded momentum advection and the volume-fraction flux; the mass
!! and total-energy face fluxes are zeroed there. This module advances the
!! conserved state by `n_micro` cheap second-order forward-backward microsteps
!! that transport mass and total energy and perform the acoustic (pressure)
!! work, with the slow forcing held frozen.
!!
!! Per microstep (5-equation model, total-energy formulation, step size dtau):
!!   1. narrow halo exchange of the conserved field
!!   2. forward : advance alpha_k*rho_k and rho*E by the second-order centered
!!                divergence of (alpha_k*rho_k * u) and ((rho*E + p) * u), using
!!                the current (frozen-this-microstep) velocity u, plus
!!                dtau * rhs_slow
!!   3. EOS     : recompute the cell pressure p from the updated conserved state
!!   4a. divergence: build div_sf = centered div(rho*u) over the interior + one ghost layer
!!   4b. backward: advance rho*u by  -dtau*grad(p_new) + dtau*rhs_slow_mom
!!                 + acoustic_div_damp*[Dx_i^2 * d/dx_i(div(rho*u))]   (grad-div divergence damping)
!!
!! The divergence damping is a grad-div smoother that suppresses ONLY compressive
!! (nonzero-divergence) acoustic noise. The damping increment is the CENTERED gradient of
!! the CENTERED momentum divergence div_sf, so the two share one discrete stencil: the
!! resulting operator is exactly -v v^T in Fourier space (rank one), giving a zero eigenvalue
!! on every discretely divergence-free mode -> vortical flow is left untouched to machine
!! precision (a per-component Laplacian would instead dissipate it). acoustic_div_damp is
!! DIMENSIONLESS: the Dx_i^2 normalization cancels the grid spacing (raw differences on a
!! uniform grid), NOT scaled by dtau and NOT divided by Dx^2. The discrete operator's most
!! negative eigenvalue is -num_dims, so the explicit-smoother stability bound is
!! acoustic_div_damp <= 2/num_dims per microstep (default 0.1 is well inside this for 1-3D).
!!
!! The mixture coefficients (gamma, pi_inf, qv) come from the cell-local volume
!! fractions and partial densities. For num_fluids>1 the volume fractions evolve by
!! pure slow (material) advection (dtau*rhs_slow_vf(%adv) each microstep); their
!! acoustic div(u) coupling is O(M^2) and dropped in the low-Mach limit. The EOS
!! pressure is therefore cell-local and cheap.
module m_acoustic_substep

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_boundary_common
    use m_weno, only: s_weno

    implicit none

    private
    public :: s_initialize_acoustic_substep_module, s_acoustic_substep, s_finalize_acoustic_substep_module

    !> Cell pressure scratch (working precision), spanning the buffered domain.
    real(wp), allocatable, dimension(:,:,:) :: p_sf
    $:GPU_DECLARE(create='[p_sf]')

    !> Momentum-divergence scratch div(rho*u) (working precision), spanning the buffered domain. Built from the halo-valid momentum
    !! each microstep, so its centered gradient (the grad-div damping) shares the divergence stencil and exactly annihilates
    !! divergence-free vortical modes.
    real(wp), allocatable, dimension(:,:,:) :: div_sf
    $:GPU_DECLARE(create='[div_sf]')

    !> Start-of-microstep snapshot of the transported conserved vars (partial densities 1:num_fluids, then total energy in slot
    !! num_fluids+1), spanning the buffered domain. The forward sweep reads its centered flux stencil ONLY from this frozen state,
    !! so the flux divergence telescopes and species mass is conserved for any num_fluids (no in-place read hazard).
    real(wp), allocatable, dimension(:,:,:,:) :: q_snap
    $:GPU_DECLARE(create='[q_snap]')

    !> Primitive acoustic sub-vector reconstructed at faces: slot 1 = pressure, slot 2 = direction-normal velocity, slot 3 = mixture
    !! density. These are exactly what the flagged-face acoustic flux consumes (sound speed is derived from p and rho, not
    !! reconstructed separately).
    integer, parameter :: n_acoustic = 3
    integer, parameter :: acou_p = 1, acou_u = 2, acou_rho = 3

    !> Primitive acoustic sub-vector over the buffered domain, rebuilt each microstep as the WENO reconstruction input (pressure and
    !! density once, normal velocity refreshed per direction).
    type(scalar_field), allocatable, dimension(:) :: q_acoustic_vf
    $:GPU_DECLARE(create='[q_acoustic_vf]')

    !> WENO-reconstructed left/right acoustic face states (slots: pressure, normal velocity, mixture density), reused per direction;
    !! the reconstructed states feed the flagged-face acoustic flux (A4).
    real(wp), allocatable, dimension(:,:,:,:) :: acL_rs_vf, acR_rs_vf
    $:GPU_DECLARE(create='[acL_rs_vf, acR_rs_vf]')

    !> Per-face discontinuity mask (1 = flagged -> robust WENO+HLLC-acoustic tier, 0 = smooth -> cheap centered tier), last index =
    !! face-normal direction (1=x, 2=y, 3=z), addressed at cell (j,k,l) for the face between (j,k,l) and its +dir neighbour. Produced
    !! here (A3) and consumed by the flux tier selection (A4); integer is the cheapest GPU-portable form of a per-face boolean.
    integer, allocatable, dimension(:,:,:,:) :: acoustic_flag
    $:GPU_DECLARE(create='[acoustic_flag]')

    !> Dimensionless smoothness factor for the WENO-indicator face test (no physical scale). The WENO smoothness indicators are
    !! squared cell differences; a face is flagged when its across-face indicator (p_{j+1}-p_j)^2 exceeds the larger neighbour
    !! indicator by more than this factor -- i.e. the gradient across the face is more than twice the steepest neighbouring
    !! gradient, which a smooth resolved field cannot produce (slopes vary slowly) but a discontinuity does by orders of magnitude.
    !! 2 is the natural dyadic relative factor and matches the beta-ratio scale at which WENO's nonlinear weights begin to
    !! discriminate stencils; it sits in the wide gap above the worst smooth value measured (~1.05) and far below any discontinuity.
    real(wp), parameter :: weno_smooth_factor = 2._wp

    !> Volume-fraction jump tolerance flagging a material interface. alpha in [0,1] is dimensionless, so this is a pure jump in a
    !! [0,1] field: far above roundoff, far below any real interface jump (O(0.1-1)). For num_fluids==1 alpha==1 everywhere, so the
    !! jump is ~0 and this flag stays false.
    real(wp), parameter :: eps_alpha = 1.e-6_wp

contains

    !> Allocate scratch storage for the acoustic substep kernel.
    impure subroutine s_initialize_acoustic_substep_module

        integer :: i  !< Acoustic sub-vector slot iterator

        @:ALLOCATE(p_sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ALLOCATE(div_sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ALLOCATE(q_snap(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:num_fluids + 1))

        @:ALLOCATE(q_acoustic_vf(1:n_acoustic))
        do i = 1, n_acoustic
            @:ALLOCATE(q_acoustic_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                       & idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_acoustic_vf(i))
        end do
        @:ALLOCATE(acL_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:n_acoustic))
        @:ALLOCATE(acR_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:n_acoustic))
        @:ALLOCATE(acoustic_flag(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:num_dims))

    end subroutine s_initialize_acoustic_substep_module

    !> Advance the conserved field in place by `n_micro` forward-backward acoustic microsteps.
    !! @param q_cons_vf    Conserved field (advanced in place).
    !! @param rhs_slow_vf  Frozen slow (advective) forcing, added dtau-scaled each microstep.
    !! @param bc_type      Boundary-condition type fields for the halo exchange.
    !! @param dtau         Microstep size (= dt_stage / n_substeps_stage).
    !! @param n_micro      Number of microsteps to take.
    impure subroutine s_acoustic_substep(q_cons_vf, rhs_slow_vf, bc_type, dtau, n_micro)

        type(scalar_field), dimension(sys_size), intent(inout)     :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in)        :: rhs_slow_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp), intent(in)                                       :: dtau
        integer, intent(in)                                        :: n_micro
        real(wp)                                                   :: rho  !< Mixture density max(sum_k alpha_k*rho_k, sgm_eps)
        !> Frozen mixture stiffened-gas coefficients (sum_k alpha_k*{gammas,pi_infs}, sum_k alpha_rho_k*qvs)
        real(wp) :: gamma, pinf, qv
        real(wp) :: rho_n                               !< Neighbour-cell mixture density (velocity denominator)
        real(wp) :: dyn_p                               !< Dynamic pressure 0.5*rho*|u|^2
        real(wp) :: div                                 !< Centered flux divergence accumulator
        real(wp) :: dpdx                                !< Centered pressure-gradient component
        real(wp) :: damp                                !< Grad-div divergence-damping stencil accumulator
        real(wp) :: inv_2dx, inv_2dy, inv_2dz           !< 1/(2*dx) etc. for centered first differences
        real(wp) :: u_xp, u_xm, u_yp, u_ym, u_zp, u_zm  !< Stencil-neighbour cell-center velocities (computed once per cell)
        integer  :: q                                   !< Microstep iterator
        integer  :: f                                   !< Fluid iterator (mixture sums)
        integer  :: i, j, k, l                          !< Generic / spatial iterators
        integer  :: momxb, momxe, contb                 !< First/last momentum and first continuity equation indices
        integer  :: esnap                               !< Total-energy slot in the q_snap snapshot (= num_fluids + 1)
        integer  :: db_x, de_x, db_y, de_y, db_z, de_z  !< div_sf loop bounds (interior + one ghost layer)
        integer  :: advb, adve                          !< First/last volume-fraction (advection) equation indices
        integer  :: joff, koff, loff                    !< +1 face-neighbour offset for the current direction (one nonzero)
        real(wp) :: eps_w                               !< Local copy of weno_eps (loop-invariant module scalar)
        real(wp) :: d_m, d_0, d_p                       !< Acoustic-pressure first differences across the face stencil
        real(wp) :: beta_m, beta_0, beta_p              !< WENO smoothness indicators (squared differences) of the stencil
        real(wp) :: beta_face, beta_nbr_max             !< Across-face indicator and larger neighbour-stencil indicator (+eps_w)
        logical  :: flag_disc                           !< Combined per-face discontinuity flag

        momxb = eqn_idx%mom%beg
        momxe = eqn_idx%mom%end
        contb = eqn_idx%cont%beg
        esnap = num_fluids + 1
        advb = eqn_idx%adv%beg
        adve = eqn_idx%adv%end
        eps_w = weno_eps

        ! div_sf is needed one ghost layer beyond the interior so its centered gradient can be taken
        ! at interior cells; momentum carries >= 2 ghost layers, so this layer is always halo-valid.
        db_x = idwint(1)%beg - 1; de_x = idwint(1)%end + 1
        if (n > 0) then; db_y = idwint(2)%beg - 1; de_y = idwint(2)%end + 1; else; db_y = 0; de_y = 0; end if
        if (p > 0) then; db_z = idwint(3)%beg - 1; de_z = idwint(3)%end + 1; else; db_z = 0; de_z = 0; end if

        do q = 1, n_micro
            ! 1. Narrow halo exchange of the conserved field (alpha_k*rho_k, rho*u, rho*E
            !    are all carried in q_cons_vf; the centered stencil needs one ghost cell).
            call s_populate_variables_buffers(bc_type, q_cons_vf)

            ! 2a. Cell pressure from the CURRENT conserved state (used by the energy flux), via the
            !     direct mixture stiffened-gas EOS (model_eqns=2): inverting
            !     rho*E = gamma*p + pi_inf + 0.5*rho*|u|^2 + qv with frozen mixture coefficients
            !     rho=sum_k alpha_rho_k, gamma=sum_k alpha_k*gammas(k), pi_inf=sum_k alpha_k*pi_infs(k),
            !     qv=sum_k alpha_rho_k*qvs(k) (matches s_convert_species_to_mixture_variables_acc).
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, f, rho, gamma, pinf, qv, dyn_p]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        rho = 0._wp; gamma = 0._wp; pinf = 0._wp; qv = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do f = 1, num_fluids
                            rho = rho + q_cons_vf(eqn_idx%cont%beg + f - 1)%sf(j, k, l)
                            qv = qv + q_cons_vf(eqn_idx%cont%beg + f - 1)%sf(j, k, l)*qvs(f)
                            gamma = gamma + q_cons_vf(eqn_idx%adv%beg + f - 1)%sf(j, k, l)*gammas(f)
                            pinf = pinf + q_cons_vf(eqn_idx%adv%beg + f - 1)%sf(j, k, l)*pi_infs(f)
                        end do
                        rho = max(rho, sgm_eps)

                        dyn_p = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = momxb, momxe
                            dyn_p = dyn_p + 0.5_wp*q_cons_vf(i)%sf(j, k, l)*q_cons_vf(i)%sf(j, k, l)/rho
                        end do

                        p_sf(j, k, l) = (q_cons_vf(eqn_idx%E)%sf(j, k, l) - dyn_p - pinf - qv)/gamma

                        ! Freeze the transported conserved state for the forward sweep's flux stencil. Reading the
                        ! divergence from this snapshot (never an in-place-updated neighbour) makes it telescope, so
                        ! species mass and energy are conserved for any num_fluids.
                        $:GPU_LOOP(parallelism='[seq]')
                        do f = 1, num_fluids
                            q_snap(j, k, l, f) = q_cons_vf(eqn_idx%cont%beg + f - 1)%sf(j, k, l)
                        end do
                        q_snap(j, k, l, esnap) = q_cons_vf(eqn_idx%E)%sf(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! p_sf is filled over the full buffered range (idwbuff) from the already-exchanged
            ! q_cons_vf, so the centered energy flux at interior edges sees a consistent p in
            ! the one ghost cell it touches. No separate pressure halo exchange is needed.

            ! 2b. FORWARD: advance alpha_k*rho_k and rho*E by the centered divergence of
            !     (alpha_k*rho_k * u) and ((rho*E + p) * u), plus dtau*rhs_slow. Velocity u
            !     is held frozen (computed from the current momentum) for this forward sweep.
            inv_2dx = 0._wp; inv_2dy = 0._wp; inv_2dz = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, f, rho_n, div, inv_2dx, inv_2dy, inv_2dz, u_xp, u_xm, u_yp, &
                                & u_ym, u_zp, u_zm]')
            do l = idwint(3)%beg, idwint(3)%end
                do k = idwint(2)%beg, idwint(2)%end
                    do j = idwint(1)%beg, idwint(1)%end
                        inv_2dx = 1._wp/(2._wp*dx(j))

                        ! Stencil-neighbour velocities u = mom/rho_mix, formed once per cell from the FROZEN snapshot
                        ! mixture density (sum_k alpha_k*rho_k) and reused by the mass and energy fluxes. Momentum is
                        ! untouched by the forward sweep, so mom == its frozen value; reading rho from the snapshot
                        ! makes u_j a single-valued function of the frozen state, which is what lets the flux telescope.
                        rho_n = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do f = 1, num_fluids
                            rho_n = rho_n + q_snap(j + 1, k, l, f)
                        end do
                        u_xp = q_cons_vf(momxb)%sf(j + 1, k, l)/max(rho_n, sgm_eps)
                        rho_n = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do f = 1, num_fluids
                            rho_n = rho_n + q_snap(j - 1, k, l, f)
                        end do
                        u_xm = q_cons_vf(momxb)%sf(j - 1, k, l)/max(rho_n, sgm_eps)
                        if (n > 0) then
                            rho_n = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do f = 1, num_fluids
                                rho_n = rho_n + q_snap(j, k + 1, l, f)
                            end do
                            u_yp = q_cons_vf(momxb + 1)%sf(j, k + 1, l)/max(rho_n, sgm_eps)
                            rho_n = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do f = 1, num_fluids
                                rho_n = rho_n + q_snap(j, k - 1, l, f)
                            end do
                            u_ym = q_cons_vf(momxb + 1)%sf(j, k - 1, l)/max(rho_n, sgm_eps)
                        end if
                        if (p > 0) then
                            rho_n = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do f = 1, num_fluids
                                rho_n = rho_n + q_snap(j, k, l + 1, f)
                            end do
                            u_zp = q_cons_vf(momxe)%sf(j, k, l + 1)/max(rho_n, sgm_eps)
                            rho_n = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do f = 1, num_fluids
                                rho_n = rho_n + q_snap(j, k, l - 1, f)
                            end do
                            u_zm = q_cons_vf(momxe)%sf(j, k, l - 1)/max(rho_n, sgm_eps)
                        end if

                        ! Mass equations: d/dt(alpha_k*rho_k) = -div(alpha_k*rho_k * u) + slow.
                        ! Flux stencil reads the frozen snapshot; apply onto the live field from the frozen center value.
                        $:GPU_LOOP(parallelism='[seq]')
                        do f = 1, num_fluids
                            div = inv_2dx*(u_xp*q_snap(j + 1, k, l, f) - u_xm*q_snap(j - 1, k, l, f))
                            if (n > 0) then
                                inv_2dy = 1._wp/(2._wp*dy(k))
                                div = div + inv_2dy*(u_yp*q_snap(j, k + 1, l, f) - u_ym*q_snap(j, k - 1, l, f))
                            end if
                            if (p > 0) then
                                inv_2dz = 1._wp/(2._wp*dz(l))
                                div = div + inv_2dz*(u_zp*q_snap(j, k, l + 1, f) - u_zm*q_snap(j, k, l - 1, f))
                            end if
                            i = eqn_idx%cont%beg + f - 1
                            q_cons_vf(i)%sf(j, k, l) = q_snap(j, k, l, f) + dtau*(-div + rhs_slow_vf(i)%sf(j, k, l))
                        end do

                        ! Total-energy equation: d/dt(rho*E) = -div((rho*E + p) * u) + slow. p_sf is the frozen
                        ! cell pressure (built from the same pre-sweep state as q_snap), so the energy flux is consistent.
                        div = inv_2dx*((q_snap(j + 1, k, l, esnap) + p_sf(j + 1, k, l))*u_xp - (q_snap(j - 1, k, l, &
                                       & esnap) + p_sf(j - 1, k, l))*u_xm)
                        if (n > 0) then
                            inv_2dy = 1._wp/(2._wp*dy(k))
                            div = div + inv_2dy*((q_snap(j, k + 1, l, esnap) + p_sf(j, k + 1, l))*u_yp - (q_snap(j, k - 1, l, &
                                                 & esnap) + p_sf(j, k - 1, l))*u_ym)
                        end if
                        if (p > 0) then
                            inv_2dz = 1._wp/(2._wp*dz(l))
                            div = div + inv_2dz*((q_snap(j, k, l + 1, esnap) + p_sf(j, k, l + 1))*u_zp - (q_snap(j, k, l - 1, &
                                                 & esnap) + p_sf(j, k, l - 1))*u_zm)
                        end if
                        q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_snap(j, k, l, esnap) + dtau*(-div + rhs_slow_vf(eqn_idx%E)%sf(j, k, l))

                        ! Volume fractions: pure slow (material) advection. The frozen slow forcing
                        ! rhs_slow_vf(%adv) already holds the contact-upwinded HLLC volume-fraction
                        ! flux, integrating to dt_stage*rhs_slow over the stage. The acoustic div(u)
                        ! coupling to alpha is O(M^2) and intentionally dropped in the low-Mach limit.
                        if (num_fluids > 1) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) + dtau*rhs_slow_vf(i)%sf(j, k, l)
                            end do
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! 3. Halo exchange + EOS on the UPDATED state, so the backward pressure gradient
            !    sees the freshly advanced mass and energy.
            call s_populate_variables_buffers(bc_type, q_cons_vf)

            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, f, rho, gamma, pinf, qv, dyn_p]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        rho = 0._wp; gamma = 0._wp; pinf = 0._wp; qv = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do f = 1, num_fluids
                            rho = rho + q_cons_vf(eqn_idx%cont%beg + f - 1)%sf(j, k, l)
                            qv = qv + q_cons_vf(eqn_idx%cont%beg + f - 1)%sf(j, k, l)*qvs(f)
                            gamma = gamma + q_cons_vf(eqn_idx%adv%beg + f - 1)%sf(j, k, l)*gammas(f)
                            pinf = pinf + q_cons_vf(eqn_idx%adv%beg + f - 1)%sf(j, k, l)*pi_infs(f)
                        end do
                        rho = max(rho, sgm_eps)

                        dyn_p = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = momxb, momxe
                            dyn_p = dyn_p + 0.5_wp*q_cons_vf(i)%sf(j, k, l)*q_cons_vf(i)%sf(j, k, l)/rho
                        end do

                        p_sf(j, k, l) = (q_cons_vf(eqn_idx%E)%sf(j, k, l) - dyn_p - pinf - qv)/gamma
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! 3b. Two-tier acoustic reconstruction: build the primitive acoustic sub-vector (pressure,
            !     direction-normal velocity, mixture density) over the buffered domain and WENO-reconstruct
            !     its L/R face states per direction. p_sf and the conserved-field halos are valid here.
            !     Pressure and density are direction-independent (filled once); the normal velocity is
            !     refreshed before each directional reconstruction. The reconstructed states are wired but
            !     not yet consumed by the flux (A4) -- so this leaves the smooth/unflagged answer unchanged.
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, f, rho]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        rho = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do f = 1, num_fluids
                            rho = rho + q_cons_vf(eqn_idx%cont%beg + f - 1)%sf(j, k, l)
                        end do
                        q_acoustic_vf(acou_p)%sf(j, k, l) = p_sf(j, k, l)
                        q_acoustic_vf(acou_rho)%sf(j, k, l) = max(rho, sgm_eps)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            do i = 1, num_dims
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            q_acoustic_vf(acou_u)%sf(j, k, l) = q_cons_vf(momxb + i - 1)%sf(j, k, &
                                          & l)/q_acoustic_vf(acou_rho)%sf(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_reconstruct_acoustic_boundary_values(q_acoustic_vf, acL_rs_vf, acR_rs_vf, i)

                ! 3c. Per-face discontinuity criterion for direction i: flag the robust acoustic tier (A4) when EITHER the
                !     reconstructed acoustic pressure is non-smooth across the face OR the volume fraction varies across it.
                !     The shock test reuses WENO's own smoothness indicators (no new physical threshold): the WENO beta_* are
                !     squared cell-difference stencils, so beta_0 = (p_{j+1}-p_j)^2 measures the smoothness of the stencil
                !     STRADDLING the face and beta_m/beta_p the smoothness of the left/right neighbour stencils. A face is smooth
                !     when its across-face indicator does not exceed the larger neighbour indicator by more than weno_smooth_factor
                !     (using the LARGER neighbour as the local scale keeps the test finite at smooth extrema / flat shoulders where
                !     one one-sided difference vanishes -- the failure mode that sinks a min(beta) test). At a discontinuity the
                !     across-face difference dwarfs the smooth-side differences, exactly the imbalance WENO's nonlinear weights act on.
                joff = 0; koff = 0; loff = 0
                if (i == 1) then; joff = 1; else if (i == 2) then; koff = 1; else; loff = 1; end if
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, f, d_m, d_0, d_p, beta_m, beta_0, beta_p, beta_face, beta_nbr_max, &
                                    & flag_disc]')
                do l = idwint(3)%beg, idwint(3)%end
                    do k = idwint(2)%beg, idwint(2)%end
                        do j = idwint(1)%beg, idwint(1)%end
                            d_m = q_acoustic_vf(acou_p)%sf(j, k, l) - q_acoustic_vf(acou_p)%sf(j - joff, k - koff, l - loff)
                            d_0 = q_acoustic_vf(acou_p)%sf(j + joff, k + koff, l + loff) - q_acoustic_vf(acou_p)%sf(j, k, l)
                            d_p = q_acoustic_vf(acou_p)%sf(j + 2*joff, k + 2*koff, l + 2*loff) &
                                  - q_acoustic_vf(acou_p)%sf(j + joff, k + koff, l + loff)
                            beta_m = d_m*d_m; beta_0 = d_0*d_0; beta_p = d_p*d_p
                            beta_face = beta_0                       ! across-face smoothness indicator
                            beta_nbr_max = max(beta_m, beta_p) + eps_w  ! larger neighbour-stencil indicator (finite local scale)
                            flag_disc = (beta_face > weno_smooth_factor*beta_nbr_max)

                            ! Volume-fraction (material-interface) variation across the face, any fluid.
                            $:GPU_LOOP(parallelism='[seq]')
                            do f = advb, adve
                                if (abs(q_cons_vf(f)%sf(j + joff, k + koff, l + loff) - q_cons_vf(f)%sf(j, k, l)) > eps_alpha) then
                                    flag_disc = .true.
                                end if
                            end do

                            if (flag_disc) then
                                acoustic_flag(j, k, l, i) = 1
                            else
                                acoustic_flag(j, k, l, i) = 0
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do

            ! 4a. DIVERGENCE: build div_sf = centered div(rho*u) over the interior + one ghost
            !     layer. The momentum halo is current from the exchange above and was untouched
            !     by the forward sweep, so all neighbours are valid.
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, div, inv_2dx, inv_2dy, inv_2dz]')
            do l = db_z, de_z
                do k = db_y, de_y
                    do j = db_x, de_x
                        inv_2dx = 1._wp/(2._wp*dx(j))
                        div = inv_2dx*(q_cons_vf(momxb)%sf(j + 1, k, l) - q_cons_vf(momxb)%sf(j - 1, k, l))
                        if (n > 0) then
                            inv_2dy = 1._wp/(2._wp*dy(k))
                            div = div + inv_2dy*(q_cons_vf(momxb + 1)%sf(j, k + 1, l) - q_cons_vf(momxb + 1)%sf(j, k - 1, l))
                        end if
                        if (p > 0) then
                            inv_2dz = 1._wp/(2._wp*dz(l))
                            div = div + inv_2dz*(q_cons_vf(momxe)%sf(j, k, l + 1) - q_cons_vf(momxe)%sf(j, k, l - 1))
                        end if
                        div_sf(j, k, l) = div
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! 4b. BACKWARD: advance rho*u by  -dtau*grad(p_new) + dtau*rhs_slow_mom
            !     + acoustic_div_damp*[Dx_i^2 * d/dx_i(div(rho*u))]   (grad-div divergence damping).
            !     The damping is the centered gradient of div_sf scaled by Dx_i^2 so that
            !     acoustic_div_damp is dimensionless; sharing the divergence stencil makes the
            !     operator annihilate divergence-free vortical modes exactly.
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, dpdx, damp, inv_2dx, inv_2dy, inv_2dz]')
            do l = idwint(3)%beg, idwint(3)%end
                do k = idwint(2)%beg, idwint(2)%end
                    do j = idwint(1)%beg, idwint(1)%end
                        inv_2dx = 1._wp/(2._wp*dx(j))

                        ! x-momentum
                        dpdx = inv_2dx*(p_sf(j + 1, k, l) - p_sf(j - 1, k, l))
                        damp = 0.5_wp*dx(j)*(div_sf(j + 1, k, l) - div_sf(j - 1, k, l))
                        q_cons_vf(momxb)%sf(j, k, l) = q_cons_vf(momxb)%sf(j, k, l) + dtau*(-dpdx + rhs_slow_vf(momxb)%sf(j, k, &
                                  & l)) + acoustic_div_damp*damp

                        ! y-momentum
                        if (n > 0) then
                            inv_2dy = 1._wp/(2._wp*dy(k))
                            dpdx = inv_2dy*(p_sf(j, k + 1, l) - p_sf(j, k - 1, l))
                            damp = 0.5_wp*dy(k)*(div_sf(j, k + 1, l) - div_sf(j, k - 1, l))
                            q_cons_vf(momxb + 1)%sf(j, k, l) = q_cons_vf(momxb + 1)%sf(j, k, &
                                      & l) + dtau*(-dpdx + rhs_slow_vf(momxb + 1)%sf(j, k, l)) + acoustic_div_damp*damp
                        end if

                        ! z-momentum
                        if (p > 0) then
                            inv_2dz = 1._wp/(2._wp*dz(l))
                            dpdx = inv_2dz*(p_sf(j, k, l + 1) - p_sf(j, k, l - 1))
                            damp = 0.5_wp*dz(l)*(div_sf(j, k, l + 1) - div_sf(j, k, l - 1))
                            q_cons_vf(momxe)%sf(j, k, l) = q_cons_vf(momxe)%sf(j, k, l) + dtau*(-dpdx + rhs_slow_vf(momxe)%sf(j, &
                                      & k, l)) + acoustic_div_damp*damp
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_acoustic_substep

    !> WENO-reconstruct the acoustic sub-vector's L/R face states for one direction. Mirrors m_rhs's
    !! s_reconstruct_cell_boundary_values: pick the per-direction reconstruction bounds (interior trimmed by weno_polyn on the
    !! normal axis) and dispatch to s_weno.
    !! @param v_vf      Primitive acoustic sub-vector (pressure, normal velocity, mixture density).
    !! @param vL/vR     Reconstructed left/right face states (slot index matches v_vf).
    !! @param norm_dir  Reconstruction direction (1=x, 2=y, 3=z).
    subroutine s_reconstruct_acoustic_boundary_values(v_vf, vL, vR, norm_dir)

        type(scalar_field), dimension(1:n_acoustic), intent(in)                                :: v_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL, vR
        integer, intent(in)                                                                    :: norm_dir
        type(int_bounds_info)                                                                  :: is1, is2, is3
        integer                                                                                :: recon_dir

        if (norm_dir == 1) then
            is1 = idwbuff(1); is2 = idwbuff(2); is3 = idwbuff(3); recon_dir = 1
        else if (norm_dir == 2) then
            is1 = idwbuff(2); is2 = idwbuff(1); is3 = idwbuff(3); recon_dir = 2
        else
            is1 = idwbuff(3); is2 = idwbuff(2); is3 = idwbuff(1); recon_dir = 3
        end if
        is1%beg = is1%beg + weno_polyn; is1%end = is1%end - weno_polyn

        call s_weno(v_vf, vL, vR, recon_dir, is1, is2, is3)

    end subroutine s_reconstruct_acoustic_boundary_values

    !> Deallocate scratch storage.
    impure subroutine s_finalize_acoustic_substep_module

        integer :: i  !< Acoustic sub-vector slot iterator

        @:DEALLOCATE(p_sf)
        @:DEALLOCATE(div_sf)
        @:DEALLOCATE(q_snap)

        do i = 1, n_acoustic
            @:DEALLOCATE(q_acoustic_vf(i)%sf)
        end do
        @:DEALLOCATE(q_acoustic_vf)
        @:DEALLOCATE(acL_rs_vf)
        @:DEALLOCATE(acR_rs_vf)
        @:DEALLOCATE(acoustic_flag)

    end subroutine s_finalize_acoustic_substep_module

end module m_acoustic_substep
