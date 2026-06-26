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
    use m_riemann_solver_hllc, only: s_acoustic_face_flux, s_convective_face_flux

    implicit none

    private
    public :: s_initialize_acoustic_substep_module, s_acoustic_substep, s_finalize_acoustic_substep_module, frozen_conv_mom_face

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

    !> Full primitive face-reconstruction vector over the buffered domain, slots 1..eqn_idx%E in the standard primitive layout
    !! (partial densities at %cont, ALL velocity components at %mom, pressure at %E) -- mirrors m_rhs's q_prim_qp%vf(1:E) so the
    !! flagged-face convective HLLC flux (next task) consumes states in the same slot order. Rebuilt only on flagged microsteps.
    type(scalar_field), allocatable, dimension(:) :: q_prim_vf
    $:GPU_DECLARE(create='[q_prim_vf]')

    !> WENO-reconstructed left/right FULL primitive face states (slots 1..eqn_idx%E, layout matching q_prim_vf and m_rhs's
    !! qL_rsx_vf/qR_rsx_vf), reused per direction; wired here for the convective HLLC flux but not yet consumed (B2). Reconstruction
    !! is gated by any_flagged so it is skipped entirely on smooth flow.
    real(wp), allocatable, dimension(:,:,:,:) :: qL_rs_vf, qR_rs_vf
    $:GPU_DECLARE(create='[qL_rs_vf, qR_rs_vf]')

    !> Per-face discontinuity mask (1 = flagged -> robust WENO+HLLC-acoustic tier, 0 = smooth -> cheap centered tier), last index =
    !! face-normal direction (1=x, 2=y, 3=z), addressed at cell (j,k,l) for the face between (j,k,l) and its +dir neighbour.
    !! Produced here (A3) and consumed by the flux tier selection (A4); integer is the cheapest GPU-portable form of a per-face
    !! boolean.
    integer, allocatable, dimension(:,:,:,:) :: acoustic_flag
    $:GPU_DECLARE(create='[acoustic_flag]')

    !> Robust-tier face-flux deltas at flagged faces, stored per face-normal direction so the per-direction WENO buffers
    !! (acL_rs_vf/acR_rs_vf/qL_rs_vf/qR_rs_vf, overwritten each direction) can be reused. dEflux = full live HLLC energy face flux
    !! (convective E*u + acoustic pressure work) MINUS the full centered energy face flux (the energy delta consumed by the forward
    !! sweep); dmassflux(:,:,:,dir,f) = live HLLC partial-density (alpha*rho_k) face flux MINUS the centered convective mass flux,
    !! per fluid (the mass delta consumed by the forward sweep). flux_mom_rob = robust normal-momentum (acoustic star-pressure) flux
    !! consumed by the backward sweep; the convective normal+transverse momentum transport stays in the frozen slow forcing
    !! rhs_slow_vf (already contact-upwinded HLLC, frozen at stage entry). All are addressed at cell (j,k,l) for the +dir face and
    !! are exactly zero at unflagged (smooth) faces, so the centered tier is left bit-identical there.
    real(wp), allocatable, dimension(:,:,:,:) :: dEflux, flux_mom_rob
    $:GPU_DECLARE(create='[dEflux, flux_mom_rob]')
    real(wp), allocatable, dimension(:,:,:,:,:) :: dmassflux
    $:GPU_DECLARE(create='[dmassflux]')

    !> Live momentum machinery at flagged faces (per face-normal direction, per velocity component). frozen_conv_mom_face
    !! = the EXACT per-face slow momentum flux the slow (advective) HLLC path computed at the RK stage state U^(s-1) --
    !! whose divergence IS the momentum part of the frozen slow forcing rhs_slow_vf(mom). NOTE: despite the "conv" in the name
    !! this is the FULL slow momentum flux INCLUDING the star-pressure p_Star on the normal component, not a convective-only
    !! rho*u_n*u_d flux; storing p_Star is LOAD-BEARING (it cancels the slow pressure so the acoustic tier flux_mom_rob re-adds
    !! it live without double-counting -- a "convective-only" cleanup would silently break this).
    !! It is PUBLIC and populated by m_rhs (s_compute_rhs) directly from flux_n(dir)%vf(mom) each stage, in GLOBAL
    !! momentum-component layout (component d = flux_n(dir)%vf(mom%beg+d-1)), addressed at cell (j,k,l) for the +dir face -- the same
    !! convention as dmomflux below. Storing the slow path's actual flux (not a reconstruction from U^n) makes the cancellation
    !! exact on EVERY RK stage. dmomflux = live HLLC convective momentum flux (re-evaluated each microstep, mapped from flux_mom's
    !! rotated slots: slot 1 = face-normal, 2..nv = transverse) MINUS frozen_conv_mom_face, the per-face delta consumed by the
    !! backward sweep so the NET convective momentum becomes live full-HLLC (the frozen part cancels the slow-forcing convective
    !! momentum) while viscous/source momentum stays frozen-slow. dmomflux is exactly zero at smooth faces so the centered tier is
    !! left bit-identical there.
    real(wp), allocatable, dimension(:,:,:,:,:) :: frozen_conv_mom_face, dmomflux
    $:GPU_DECLARE(create='[frozen_conv_mom_face, dmomflux]')

    !> Dimensionless max-indicator-dominance ratio above which a face's pressure stencil is judged discontinuous. The face test
    !! compares the LARGEST WENO smoothness indicator beta_max = max(beta_m, beta_0, beta_p) (squared cell differences of the
    !! acoustic pressure; Borges et al. 2008) to the relative pressure-squared floor beta_floor = (eps_rel*p_scale)^2 (see eps_rel),
    !! never to weno_eps, so the ratio is dimensionless, pressure-scale-invariant, and independent of the user's WENO epsilon. The
    !! ratio equals (d_max/p_scale)^2 / eps_rel^2: a face is flagged once the largest per-cell pressure jump exceeds
    !! sqrt(disc_ratio_thresh)*eps_rel of the local pressure scale. This is the max-indicator form (NOT the WENO-Z spread tau_z =
    !! |beta_p - beta_m|, which is exactly 0 at a face-centered 1-cell jump where beta_m ~ beta_p ~ 0 and so MISSES a sharp shock):
    !! the dominance form fires on a sharp jump (a prior experiment measured ~8.1e11 at the shock) while staying ~0 on smooth flow.
    !! At disc_ratio_thresh=1e8 the flag fires when the per-cell pressure jump exceeds 1% of the local pressure scale: far above the
    !! smooth golden 7F8ED027 (M~0.01 acoustic pulse, dp/p ~ 1e-4 -> ratio ~1e4, zero flags, bit-identical) and far below a genuine
    !! shock (~1e11). The relative floor makes the whole test invariant to the pressure normalization.
    real(wp), parameter :: disc_ratio_thresh = 1.e8_wp

    !> Dimensionless relative-roundoff fraction for the face test. A face is a candidate only if its across-face pressure gradient
    !! exceeds eps_rel of the local pressure scale (max |p| over the 4-point stencil); below that the gradient is floating-point
    !! roundoff, not a resolved feature. It sets the relative median floor (eps_rel*p_scale)^2 and the candidate guard, making both
    !! scale as pressure^2 so the whole test is invariant to the pressure normalization. This is a roundoff-rejection threshold, NOT
    !! a tuned shock sensor.
    real(wp), parameter :: eps_rel = 1.e-6_wp

    !> Volume-fraction jump tolerance flagging a material interface. alpha in [0,1] is dimensionless, so this is a pure jump in a
    !! [0,1] field: far above roundoff, far below any real interface jump (O(0.1-1)). For num_fluids==1 alpha==1 everywhere, so the
    !! jump is ~0 and this flag stays false.
    real(wp), parameter :: eps_alpha = 1.e-6_wp

contains

    !> Allocate scratch storage for the acoustic substep kernel.
    impure subroutine s_initialize_acoustic_substep_module

        integer :: i      !< Acoustic sub-vector / primitive slot iterator
        integer :: nfull  !< Number of reconstructed full-primitive variables (= eqn_idx%E: partial densities, velocities, pressure)

        nfull = eqn_idx%E

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

        @:ALLOCATE(q_prim_vf(1:nfull))
        do i = 1, nfull
            @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_prim_vf(i))
        end do
        @:ALLOCATE(qL_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:nfull))
        @:ALLOCATE(qR_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:nfull))

        @:ALLOCATE(acoustic_flag(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:num_dims))
        @:ALLOCATE(dEflux(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:num_dims))
        @:ALLOCATE(flux_mom_rob(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:num_dims))
        @:ALLOCATE(dmassflux(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:num_dims, 1:num_fluids))
        @:ALLOCATE(frozen_conv_mom_face(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                   & idwbuff(3)%beg:idwbuff(3)%end, 1:num_dims, 1:num_vels))
        @:ALLOCATE(dmomflux(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:num_dims, 1:num_vels))

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
        real(wp) :: rho_n  !< Neighbour-cell mixture density (velocity denominator)
        real(wp) :: dyn_p  !< Dynamic pressure 0.5*rho*|u|^2
        real(wp) :: div  !< Centered flux divergence accumulator
        real(wp) :: dpdx  !< Centered pressure-gradient component
        real(wp) :: damp  !< Grad-div divergence-damping stencil accumulator
        real(wp) :: inv_2dx, inv_2dy, inv_2dz  !< 1/(2*dx) etc. for centered first differences
        real(wp) :: u_xp, u_xm, u_yp, u_ym, u_zp, u_zm  !< Stencil-neighbour cell-center velocities (computed once per cell)
        integer  :: q  !< Microstep iterator
        integer  :: f  !< Fluid iterator (mixture sums)
        integer  :: i, j, k, l  !< Generic / spatial iterators
        integer  :: momxb, momxe, contb  !< First/last momentum and first continuity equation indices
        integer  :: esnap  !< Total-energy slot in the q_snap snapshot (= num_fluids + 1)
        integer  :: db_x, de_x, db_y, de_y, db_z, de_z  !< div_sf loop bounds (interior + one ghost layer)
        integer  :: advb, adve  !< First/last volume-fraction (advection) equation indices
        integer  :: joff, koff, loff  !< +1 face-neighbour offset for the current direction (one nonzero)
        real(wp) :: p_m4, p_0c, p_p1, p_p2  !< Acoustic pressures over the 4-point face stencil (j-off, j, j+off, j+2off)
        real(wp) :: p_scale, beta_floor  !< Local pressure scale max|p| and the relative floor (eps_rel*p_scale)^2
        real(wp) :: d_m, d_0, d_p  !< Acoustic-pressure first differences across the face stencil
        real(wp) :: beta_m, beta_0, beta_p  !< WENO smoothness indicators (squared differences) of the stencil
        real(wp) :: beta_max  !< Largest stencil indicator max(beta_m, beta_0, beta_p) for the max-indicator-dominance test
        logical  :: flag_disc  !< Combined per-face discontinuity flag
        real(wp) :: gam_a, pin_a, gam_b, pin_b  !< Left/right cell mixture stiffened-gas coefficients for the face EOS
        real(wp) :: pres_a, rho_a, c_a  !< Left reconstructed face state (pressure, mixture density) + sound speed
        real(wp) :: pres_b, rho_b, c_b  !< Right reconstructed face state + sound speed
        real(wp) :: qv_a, qv_b, ke_a, ke_b, E_a, E_b  !< Left/right qv-sum, kinetic energy, and rebuilt total-energy density (EOS)
        real(wp) :: flmom, flE  !< HLLC acoustic face fluxes (normal momentum, pressure work)
        real(wp) :: flux_E_c  !< HLLC convective energy face flux (E*u contact-upwinded)
        real(wp) :: se_face, se_conv  !< Centered pressure-work and centered convective-energy face fluxes the robust tier replaces
        real(wp) :: u_a, u_b, cm  !< Cell-center normal velocities at the two face cells and the centered per-fluid mass face flux
        integer  :: nv, ic, d  !< Velocity-component count, transverse-fill counter, component iterator

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            !> Reconstructed L/R partial densities + HLLC convective mass flux
            real(wp), dimension(3) :: rhoYks_a, rhoYks_b, flux_cont
            !> L/R velocities (slot 1 = normal) + HLLC convective mom flux
            real(wp), dimension(3) :: vel_a_arr, vel_b_arr, flux_mom_c
        #:else
            !> Reconstructed L/R partial densities + conv mass flux
            real(wp), dimension(num_fluids) :: rhoYks_a, rhoYks_b, flux_cont
            real(wp), dimension(num_dims)   :: vel_a_arr, vel_b_arr, flux_mom_c  !< L/R velocities (slot 1 = normal) + conv mom flux
        #:endif
        real(wp) :: sm_r, sm_l, gd_r, gd_l, dfm_r, dfm_l  !< Backward momentum robust-correction face pressures/damping/deltas
        integer  :: flg_r, flg_l  !< Right/left face discontinuity flags (backward momentum correction)
        integer  :: fb_x, fb_y, fb_z  !< Flag/flux loop low bounds (one ghost below interior on the normal axis)
        integer  :: any_flag_max  !< Per-microstep MAX reduction of acoustic_flag (>0 iff any face flagged); gates reconstruction
        logical  :: any_flagged  !< .true. iff some interior(+ghost) face is flagged this microstep -> reconstruction must run

        momxb = eqn_idx%mom%beg
        momxe = eqn_idx%mom%end
        contb = eqn_idx%cont%beg
        esnap = num_fluids + 1
        advb = eqn_idx%adv%beg
        adve = eqn_idx%adv%end
        nv = momxe - momxb + 1  ! Number of velocity/momentum components (= num_dims for model_eqns=2)

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

            ! 2b'. Two-tier acoustic substep, from the FROZEN pre-forward state so it pairs with the snapshot the forward sweep
            !      transports (pre-forward pressure p_sf and frozen momentum). Build the cell-centered acoustic pressure / mixture
            !      density, flag the discontinuous faces per direction (cell values only, no reconstruction), then -- ONLY if some
            !      face is flagged this microstep -- WENO-reconstruct the acoustic sub-vector AND the full primitive state and
            !      evaluate the flagged-face flux. dEflux carries the robust-minus-centered pressure-work delta consumed by the
            !      forward energy sweep; flux_mom_rob carries the robust normal-momentum flux consumed by the backward sweep. Both
            !      are exactly zero at smooth faces, so the centered tier is left bit-identical there; and on smooth flow no face is
            !      flagged, so the (unused) reconstruction is skipped entirely -- bit-identical AND cheaper.
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

            ! FLAG PASS (FROZEN per RK stage -- computed once, on the first microstep q==1, from the stage-entry state, then held
            ! fixed for q>1): cheaper than per-microstep, consistent with the frozen slow flux, and the discontinuity moves < 1 cell
            ! per stage (advective CFL). Per direction it marks discontinuous faces from CELL values only (no reconstruction) and
            ! zeros the robust-tier face deltas (dEflux/flux_mom_rob/dmassflux), while MAX-reducing any_flag_max (per-RANK
            ! reduction:
            ! processor-boundary faces flag identically via the halo, so no MPI collective is needed). The loop runs one cell below
            ! the interior on the normal axis so every interior cell's lower (-dir) face also has a stored flag/flux. The shock test
            ! is the MAX-INDICATOR DOMINANCE form on the acoustic-pressure stencil: beta_m/beta_0/beta_p are squared cell
            ! differences
            ! (pressure^2; Borges et al. 2008) and a face is flagged when the LARGEST, beta_max = max(beta_m, beta_0, beta_p),
            ! exceeds disc_ratio_thresh times the relative roundoff floor beta_floor = (eps_rel*p_scale)^2. beta_max and beta_floor
            ! are both pressure^2, so the ratio is invariant to the pressure scale and free of weno_eps; it stays ~0 where the
            ! pressure is smooth and explodes across a genuine jump -- crucially it FIRES on a face-centered 1-cell jump (where the
            ! WENO-Z spread tau_z=|beta_p-beta_m| would be ~0 and miss the shock). A volume-fraction jump also flags the face.
            ! Zeroing
            ! the deltas once here is sufficient: the FLUX PASS overwrites only flagged faces each microstep, smooth faces stay 0.
            if (q == 1) then
                any_flag_max = 0
                do i = 1, num_dims
                    joff = 0; koff = 0; loff = 0
                    if (i == 1) then; joff = 1; else if (i == 2) then; koff = 1; else; loff = 1; end if
                    fb_x = idwint(1)%beg; fb_y = idwint(2)%beg; fb_z = idwint(3)%beg
                    if (i == 1) then; fb_x = idwint(1)%beg - 1; else if (i == 2) then; fb_y = idwint(2)%beg &
                        & - 1; else; fb_z = idwint(3)%beg - 1; end if
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, f, d, p_m4, p_0c, p_p1, p_p2, p_scale, beta_floor, d_m, &
                                        & d_0, d_p, beta_m, beta_0, beta_p, beta_max, flag_disc]', reduction='[[any_flag_max]]', &
                                        & reductionOp='[MAX]', copy='[any_flag_max]')
                    do l = fb_z, idwint(3)%end
                        do k = fb_y, idwint(2)%end
                            do j = fb_x, idwint(1)%end
                                p_m4 = q_acoustic_vf(acou_p)%sf(j - joff, k - koff, l - loff)
                                p_0c = q_acoustic_vf(acou_p)%sf(j, k, l)
                                p_p1 = q_acoustic_vf(acou_p)%sf(j + joff, k + koff, l + loff)
                                p_p2 = q_acoustic_vf(acou_p)%sf(j + 2*joff, k + 2*koff, l + 2*loff)
                                d_m = p_0c - p_m4; d_0 = p_p1 - p_0c; d_p = p_p2 - p_p1
                                beta_m = d_m*d_m; beta_0 = d_0*d_0; beta_p = d_p*d_p
                                p_scale = max(abs(p_m4), abs(p_0c), abs(p_p1), abs(p_p2))
                                beta_floor = (eps_rel*p_scale)*(eps_rel*p_scale)
                                beta_max = max(beta_m, beta_0, beta_p)
                                flag_disc = beta_max > disc_ratio_thresh*beta_floor

                                ! Volume-fraction (material-interface) variation across the face, any fluid.
                                $:GPU_LOOP(parallelism='[seq]')
                                do f = advb, adve
                                    if (abs(q_cons_vf(f)%sf(j + joff, k + koff, l + loff) - q_cons_vf(f)%sf(j, k, &
                                        & l)) > eps_alpha) then
                                        flag_disc = .true.
                                    end if
                                end do

                                if (flag_disc) then
                                    acoustic_flag(j, k, l, i) = 1
                                else
                                    acoustic_flag(j, k, l, i) = 0
                                end if
                                dEflux(j, k, l, i) = 0._wp
                                flux_mom_rob(j, k, l, i) = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do f = 1, num_fluids
                                    dmassflux(j, k, l, i, f) = 0._wp
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do d = 1, num_dims
                                    dmomflux(j, k, l, i, d) = 0._wp
                                end do
                                any_flag_max = max(any_flag_max, acoustic_flag(j, k, l, i))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
                any_flagged = any_flag_max > 0
            end if

            ! FLUX PASS (only when the stage has any flagged face; runs LIVE each microstep so the deltas track the evolving shock):
            ! build the full primitive vector, reconstruct both the acoustic sub-vector and the full primitive state per direction,
            ! and at each flagged face assemble the FULL live HLLC face flux F_full = s_convective_face_flux + s_acoustic_face_flux
            ! into the robust-tier deltas (smooth faces stay at the zeros set in the FLAG PASS). The net flagged-face update per
            ! microstep is pure live full HLLC: for MASS and ENERGY the slow flux is exactly zero (the HLLC slow flux zeroes the
            ! mass
            ! and total-energy faces), so the delta replaces the full CENTERED flux with the live HLLC flux (dmassflux for each
            ! alpha*rho_k, dEflux for convective E*u + acoustic pressure work). For MOMENTUM the convective normal+transverse
            ! transport already rides the frozen slow forcing (rhs_slow_vf, contact-upwinded HLLC frozen at stage entry), so only
            ! the
            ! CENTERED acoustic momentum is replaced by the live HLLC acoustic star-pressure flux flux_mom_rob -- the frozen-slow
            ! convective momentum is the inherent split lag (< 1 cell/stage), not retransported here.
            if (any_flagged) then
                ! Full primitive vector over the buffered domain in m_rhs slot order: partial densities at %cont (= conserved for
                ! the
                ! 5-equation model), ALL velocity components at %mom (mom/rho), pressure at %E. Rebuilt each flagged microstep.
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, f]')
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do f = 1, num_fluids
                                q_prim_vf(eqn_idx%cont%beg + f - 1)%sf(j, k, l) = q_cons_vf(eqn_idx%cont%beg + f - 1)%sf(j, k, l)
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                q_prim_vf(momxb + i - 1)%sf(j, k, l) = q_cons_vf(momxb + i - 1)%sf(j, k, &
                                          & l)/q_acoustic_vf(acou_rho)%sf(j, k, l)
                            end do
                            q_prim_vf(eqn_idx%E)%sf(j, k, l) = p_sf(j, k, l)
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

                    call s_reconstruct_boundary_values(q_acoustic_vf, acL_rs_vf, acR_rs_vf, i)
                    call s_reconstruct_boundary_values(q_prim_vf, qL_rs_vf, qR_rs_vf, i)

                    joff = 0; koff = 0; loff = 0
                    if (i == 1) then; joff = 1; else if (i == 2) then; koff = 1; else; loff = 1; end if
                    fb_x = idwint(1)%beg; fb_y = idwint(2)%beg; fb_z = idwint(3)%beg
                    if (i == 1) then; fb_x = idwint(1)%beg - 1; else if (i == 2) then; fb_y = idwint(2)%beg &
                        & - 1; else; fb_z = idwint(3)%beg - 1; end if
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, f, d, ic, gam_a, pin_a, gam_b, pin_b, pres_a, rho_a, c_a, &
                                        & pres_b, rho_b, c_b, qv_a, qv_b, ke_a, ke_b, E_a, E_b, flmom, flE, flux_E_c, se_face, &
                                        & se_conv, u_a, u_b, cm, rhoYks_a, rhoYks_b, vel_a_arr, vel_b_arr, flux_cont, flux_mom_c]')
                    do l = fb_z, idwint(3)%end
                        do k = fb_y, idwint(2)%end
                            do j = fb_x, idwint(1)%end
                                if (acoustic_flag(j, k, l, i) == 1) then
                                    ! Mixture stiffened-gas coefficients of the two cells straddling the face (frozen volume
                                    ! fractions); sound speed and total-energy density rebuilt from the SAME EOS as the pressure
                                    ! recompute: c^2 = (1/gamma + 1)(p + pi_inf/(gamma+1))/rho and
                                    ! rho*E = gamma*p + pi_inf + qv + 0.5*rho*|u|^2.
                                    gam_a = 0._wp; pin_a = 0._wp; gam_b = 0._wp; pin_b = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do f = 1, num_fluids
                                        gam_a = gam_a + q_cons_vf(advb + f - 1)%sf(j, k, l)*gammas(f)
                                        pin_a = pin_a + q_cons_vf(advb + f - 1)%sf(j, k, l)*pi_infs(f)
                                        gam_b = gam_b + q_cons_vf(advb + f - 1)%sf(j + joff, k + koff, l + loff)*gammas(f)
                                        pin_b = pin_b + q_cons_vf(advb + f - 1)%sf(j + joff, k + koff, l + loff)*pi_infs(f)
                                    end do

                                    ! Reconstructed L/R partial densities (qL/qR_rs_vf %cont slots) and the qv-coefficient sums.
                                    rho_a = 0._wp; rho_b = 0._wp; qv_a = 0._wp; qv_b = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do f = 1, num_fluids
                                        rhoYks_a(f) = qL_rs_vf(j, k, l, contb + f - 1)
                                        rhoYks_b(f) = qR_rs_vf(j + joff, k + koff, l + loff, contb + f - 1)
                                        rho_a = rho_a + rhoYks_a(f); rho_b = rho_b + rhoYks_b(f)
                                        qv_a = qv_a + rhoYks_a(f)*qvs(f); qv_b = qv_b + rhoYks_b(f)*qvs(f)
                                    end do
                                    rho_a = max(rho_a, sgm_eps); rho_b = max(rho_b, sgm_eps)

                                    ! Velocity arrays with slot 1 = face-normal (component i), 2..nv = transverse (any order; only
                                    ! flux_cont and flux_E -- which read slot 1 -- are consumed, but the full kinetic energy below
                                    ! sums all components).
                                    vel_a_arr(1) = qL_rs_vf(j, k, l, momxb + i - 1)
                                    vel_b_arr(1) = qR_rs_vf(j + joff, k + koff, l + loff, momxb + i - 1)
                                    ic = 1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do d = 1, num_dims
                                        if (d /= i) then
                                            ic = ic + 1
                                            vel_a_arr(ic) = qL_rs_vf(j, k, l, momxb + d - 1)
                                            vel_b_arr(ic) = qR_rs_vf(j + joff, k + koff, l + loff, momxb + d - 1)
                                        end if
                                    end do

                                    pres_a = qL_rs_vf(j, k, l, eqn_idx%E)
                                    pres_b = qR_rs_vf(j + joff, k + koff, l + loff, eqn_idx%E)
                                    c_a = sqrt(max((1._wp/gam_a + 1._wp)*(pres_a + pin_a/(gam_a + 1._wp))/rho_a, sgm_eps))
                                    c_b = sqrt(max((1._wp/gam_b + 1._wp)*(pres_b + pin_b/(gam_b + 1._wp))/rho_b, sgm_eps))
                                    ke_a = 0._wp; ke_b = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do d = 1, nv
                                        ke_a = ke_a + 0.5_wp*rho_a*vel_a_arr(d)*vel_a_arr(d)
                                        ke_b = ke_b + 0.5_wp*rho_b*vel_b_arr(d)*vel_b_arr(d)
                                    end do
                                    E_a = gam_a*pres_a + pin_a + qv_a + ke_a
                                    E_b = gam_b*pres_b + pin_b + qv_b + ke_b

                                    ! Full live HLLC face flux = convective (mass alpha*rho_k*u, energy E*u) + acoustic (normal-mom
                                    ! star pressure, energy pressure work), on a single consistent reconstructed L/R state.
                                    call s_acoustic_face_flux(pres_a, vel_a_arr(1), rho_a, c_a, pres_b, vel_b_arr(1), rho_b, c_b, &
                                                              & flmom, flE)
                                    call s_convective_face_flux(num_fluids, nv, rhoYks_a, vel_a_arr, rho_a, c_a, pres_a, E_a, &
                                                                & rhoYks_b, vel_b_arr, rho_b, c_b, pres_b, E_b, flux_cont, &
                                                                & flux_mom_c, flux_E_c)

                                    ! Centered face fluxes the forward sweep applies, from the same frozen snapshot (q_snap),
                                    ! pre-forward pressure and cell normal velocity u = mom_i/rho the sweep uses:
                                    ! convective mass 0.5*(arho_j u_j + arho_{j+1} u_{j+1}), convective energy 0.5*(E_j u_j +
                                    ! E_{j+1} u_{j+1}), pressure work 0.5*(p_j u_j + p_{j+1} u_{j+1}). Each delta replaces that
                                    ! centered face flux with the live HLLC flux; smooth faces (zeroed in the FLAG PASS) stay
                                    ! bit-identical, and each delta is single-valued per face so mass and energy telescope.
                                    u_a = q_acoustic_vf(acou_u)%sf(j, k, l)
                                    u_b = q_acoustic_vf(acou_u)%sf(j + joff, k + koff, l + loff)
                                    se_face = 0.5_wp*(q_acoustic_vf(acou_p)%sf(j, k, l)*u_a + q_acoustic_vf(acou_p)%sf(j + joff, &
                                                      & k + koff, l + loff)*u_b)
                                    se_conv = 0.5_wp*(q_snap(j, k, l, esnap)*u_a + q_snap(j + joff, k + koff, l + loff, esnap)*u_b)
                                    flux_mom_rob(j, k, l, i) = flmom
                                    dEflux(j, k, l, i) = (flux_E_c + flE) - (se_conv + se_face)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do f = 1, num_fluids
                                        cm = 0.5_wp*(q_snap(j, k, l, f)*u_a + q_snap(j + joff, k + koff, l + loff, f)*u_b)
                                        dmassflux(j, k, l, i, f) = flux_cont(f) - cm
                                    end do

                                    ! Live convective-momentum delta = live HLLC convective momentum flux (flux_mom_c) MINUS the
                                    ! frozen slow per-face flux frozen_conv_mom_face (populated by m_rhs from the slow path's
                                    ! flux_n(mom) at the RK stage state U^(s-1)), mapping the rotated live flux (slot 1 = face-normal
                                    ! component i, slots 2..nv = transverse) back to the global velocity component. Subtracting the
                                    ! slow path's ACTUAL flux makes the frozen part cancel the convective momentum carried by
                                    ! rhs_slow_vf(mom) EXACTLY on every RK stage, so the net convective momentum is live full HLLC.
                                    dmomflux(j, k, l, i, i) = flux_mom_c(1) - frozen_conv_mom_face(j, k, l, i, i)
                                    ic = 1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do d = 1, num_dims
                                        if (d /= i) then
                                            ic = ic + 1
                                            dmomflux(j, k, l, i, d) = flux_mom_c(ic) - frozen_conv_mom_face(j, k, l, i, d)
                                        end if
                                    end do
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if

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

                        ! Robust-tier mass correction: -dtau*div(dmassflux), where dmassflux(:,:,:,dir,f) = live HLLC alpha*rho_k
                        ! face flux - centered convective mass flux at flagged faces and 0 at smooth faces. This swaps the centered
                        ! convective alpha*rho_k transport just applied for the upwinded HLLC one (the slow mass flux is zero, so
                        ! there is nothing else to subtract); the divergence is single-valued per face so species mass stays
                        ! conserved, and at smooth faces it adds exactly zero -> the centered mass update is bit-identical.
                        $:GPU_LOOP(parallelism='[seq]')
                        do f = 1, num_fluids
                            i = eqn_idx%cont%beg + f - 1
                            q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) - dtau*(dmassflux(j, k, l, 1, &
                                      & f) - dmassflux(j - 1, k, l, 1, f))/dx(j)
                            if (n > 0) q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) - dtau*(dmassflux(j, k, l, 2, &
                                & f) - dmassflux(j, k - 1, l, 2, f))/dy(k)
                            if (p > 0) q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) - dtau*(dmassflux(j, k, l, 3, &
                                & f) - dmassflux(j, k, l - 1, 3, f))/dz(l)
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

                        ! Robust-tier energy correction: -dtau*div(dEflux), where dEflux = full live HLLC energy flux (convective
                        ! E*u
                        ! + acoustic pressure work) - full centered energy flux at flagged faces and 0 at smooth faces. This swaps
                        ! BOTH the convective rho*E transport and the pressure work just applied (the slow energy flux is zero) for
                        ! the upwinded HLLC energy flux; the divergence is single-valued per face (same value seen by both adjacent
                        ! cells) so energy stays conserved, and at smooth faces it adds exactly zero -> bit-identical.
                        q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, l) - dtau*(dEflux(j, k, l, &
                                  & 1) - dEflux(j - 1, k, l, 1))/dx(j)
                        if (n > 0) q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, l) - dtau*(dEflux(j, k, l, &
                            & 2) - dEflux(j, k - 1, l, 2))/dy(k)
                        if (p > 0) q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, l) - dtau*(dEflux(j, k, l, &
                            & 3) - dEflux(j, k, l - 1, 3))/dz(l)

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
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, d, dpdx, damp, inv_2dx, inv_2dy, inv_2dz, flg_r, flg_l, sm_r, &
                                & sm_l, gd_r, gd_l, dfm_r, dfm_l]')
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

                        ! Robust-tier normal-momentum correction. At a flagged face replace the centered acoustic momentum flux
                        ! (dtau*p_centered minus the grad-div damping flux gd) by the HLLC star-pressure flux dtau*flux_mom_rob and
                        ! drop the damping. Written as a single-valued face-flux delta dfm = flag*(dtau*flux_mom_rob - dtau*p_avg +
                        ! gd):
                        ! the centered pressure gradient (Sm) and damping (gd) here are the exact face-flux equivalents of the
                        ! centered backward sweep above (using the post-forward p_sf / div_sf it used), so subtracting them removes
                        ! the centered contribution of that face and the -div(dfm) update stays momentum-conservative. flag and
                        ! flux_mom_rob are zero at smooth faces, so this adds exactly zero there (centered tier bit-identical).
                        flg_r = acoustic_flag(j, k, l, 1); flg_l = acoustic_flag(j - 1, k, l, 1)
                        sm_r = 0.5_wp*(p_sf(j, k, l) + p_sf(j + 1, k, l)); sm_l = 0.5_wp*(p_sf(j - 1, k, l) + p_sf(j, k, l))
                        gd_r = acoustic_div_damp*0.5_wp*dx(j)*dx(j)*(div_sf(j, k, l) + div_sf(j + 1, k, l))
                        gd_l = acoustic_div_damp*0.5_wp*dx(j - 1)*dx(j - 1)*(div_sf(j - 1, k, l) + div_sf(j, k, l))
                        dfm_r = real(flg_r, wp)*(dtau*flux_mom_rob(j, k, l, 1) - dtau*sm_r + gd_r)
                        dfm_l = real(flg_l, wp)*(dtau*flux_mom_rob(j - 1, k, l, 1) - dtau*sm_l + gd_l)
                        q_cons_vf(momxb)%sf(j, k, l) = q_cons_vf(momxb)%sf(j, k, l) - (dfm_r - dfm_l)/dx(j)

                        if (n > 0) then
                            flg_r = acoustic_flag(j, k, l, 2); flg_l = acoustic_flag(j, k - 1, l, 2)
                            sm_r = 0.5_wp*(p_sf(j, k, l) + p_sf(j, k + 1, l)); sm_l = 0.5_wp*(p_sf(j, k - 1, l) + p_sf(j, k, l))
                            gd_r = acoustic_div_damp*0.5_wp*dy(k)*dy(k)*(div_sf(j, k, l) + div_sf(j, k + 1, l))
                            gd_l = acoustic_div_damp*0.5_wp*dy(k - 1)*dy(k - 1)*(div_sf(j, k - 1, l) + div_sf(j, k, l))
                            dfm_r = real(flg_r, wp)*(dtau*flux_mom_rob(j, k, l, 2) - dtau*sm_r + gd_r)
                            dfm_l = real(flg_l, wp)*(dtau*flux_mom_rob(j, k - 1, l, 2) - dtau*sm_l + gd_l)
                            q_cons_vf(momxb + 1)%sf(j, k, l) = q_cons_vf(momxb + 1)%sf(j, k, l) - (dfm_r - dfm_l)/dy(k)
                        end if

                        if (p > 0) then
                            flg_r = acoustic_flag(j, k, l, 3); flg_l = acoustic_flag(j, k, l - 1, 3)
                            sm_r = 0.5_wp*(p_sf(j, k, l) + p_sf(j, k, l + 1)); sm_l = 0.5_wp*(p_sf(j, k, l - 1) + p_sf(j, k, l))
                            gd_r = acoustic_div_damp*0.5_wp*dz(l)*dz(l)*(div_sf(j, k, l) + div_sf(j, k, l + 1))
                            gd_l = acoustic_div_damp*0.5_wp*dz(l - 1)*dz(l - 1)*(div_sf(j, k, l - 1) + div_sf(j, k, l))
                            dfm_r = real(flg_r, wp)*(dtau*flux_mom_rob(j, k, l, 3) - dtau*sm_r + gd_r)
                            dfm_l = real(flg_l, wp)*(dtau*flux_mom_rob(j, k, l - 1, 3) - dtau*sm_l + gd_l)
                            q_cons_vf(momxe)%sf(j, k, l) = q_cons_vf(momxe)%sf(j, k, l) - (dfm_r - dfm_l)/dz(l)
                        end if

                        ! Robust-tier convective-momentum correction: -dtau*div(dmomflux), where dmomflux(:,:,:,dir,d) = live HLLC
                        ! convective momentum flux (normal + transverse) MINUS the frozen stage-entry HLLC convective momentum flux
                        ! at flagged faces (0 at smooth faces). The frozen part cancels the convective momentum already carried by
                        ! the frozen slow forcing rhs_slow_vf(mom), so the NET convective momentum becomes live full HLLC while the
                        ! viscous/source momentum stays frozen-slow (correct). The divergence is single-valued per face (cell j
                        ! right
                        ! face = cell j+1 left face) so momentum telescopes/conserves; at smooth faces it adds exactly zero so the
                        ! centered tier (pressure gradient + grad-div damping) is left bit-identical.
                        $:GPU_LOOP(parallelism='[seq]')
                        do d = 1, num_dims
                            i = momxb + d - 1
                            q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) - dtau*(dmomflux(j, k, l, 1, d) - dmomflux(j - 1, &
                                      & k, l, 1, d))/dx(j)
                            if (n > 0) q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) - dtau*(dmomflux(j, k, l, 2, &
                                & d) - dmomflux(j, k - 1, l, 2, d))/dy(k)
                            if (p > 0) q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) - dtau*(dmomflux(j, k, l, 3, &
                                & d) - dmomflux(j, k, l - 1, 3, d))/dz(l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_acoustic_substep

    !> WENO-reconstruct a primitive vector's L/R face states for one direction. Mirrors m_rhs's s_reconstruct_cell_boundary_values:
    !! pick the per-direction reconstruction bounds (interior trimmed by weno_polyn on the normal axis) and dispatch to s_weno
    !! (which reconstructs ubound(v_vf,1) variables). Used for both the acoustic sub-vector and the full primitive vector.
    !! @param v_vf      Primitive vector (acoustic sub-vector, or full state in m_rhs slot order).
    !! @param vL/vR     Reconstructed left/right face states (slot index matches v_vf).
    !! @param norm_dir  Reconstruction direction (1=x, 2=y, 3=z).
    subroutine s_reconstruct_boundary_values(v_vf, vL, vR, norm_dir)

        type(scalar_field), dimension(1:), intent(in)                                          :: v_vf
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

    end subroutine s_reconstruct_boundary_values

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

        do i = 1, size(q_prim_vf)
            @:DEALLOCATE(q_prim_vf(i)%sf)
        end do
        @:DEALLOCATE(q_prim_vf)
        @:DEALLOCATE(qL_rs_vf)
        @:DEALLOCATE(qR_rs_vf)

        @:DEALLOCATE(acoustic_flag)
        @:DEALLOCATE(dEflux)
        @:DEALLOCATE(flux_mom_rob)
        @:DEALLOCATE(dmassflux)
        @:DEALLOCATE(frozen_conv_mom_face)
        @:DEALLOCATE(dmomflux)

    end subroutine s_finalize_acoustic_substep_module

end module m_acoustic_substep
