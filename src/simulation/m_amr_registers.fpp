!>
!!@file
!!@brief Contains module m_amr_registers

#:include 'macros.fpp'

!> @brief AMR flux registers: per-RK-stage refluxing at the coarse/fine block boundary (SP4). Depends only on m_derived_types +
!! m_global_parameters so both m_rhs (capture) and m_time_steppers (apply) can use it without cycles. NOTE: m_amr uses m_rhs which
!! uses m_amr_registers, so adding "use m_amr" here would create a compilation cycle. Region info is therefore read from
!! amr_region_lo/hi and amr_isect_lo/hi (m_global_parameters), which s_set_amr_fine_geometry keeps mirroring across regrids. creg
!! uses 0-based transverse indexing relative to the rank's block INTERSECTION (= the block at np=1); freg uses 0-based LOCAL fine
!! indexing (aligned: fine children of isect cell t are 2*t and 2*t+1). All arrays are preallocated at max size so regrid requires
!! no reallocation.
!!
!! Multi-fluid (5-eq HLLC, the amr-gated path): the volume-fraction ADVECTIVE flux alpha_i*u_star travels through flux_n (the
!! "VOLUME FRACTION FLUX" block of m_riemann_solver_hllc, same form as the mass flux), so the uniform 1:sys_size capture below
!! refluxes the per-fluid masses, momentum, energy, AND alpha's advective part with no extra registers. The non-conservative
!! remainder (the +alpha*d(u_star)/dx compression term m_rhs assembles from flux_src_n = u_star) is deliberately NOT captured:
!! alpha is genuinely non-conservative, so forcing flux-matching on u_star would be wrong; coarse/fine volume-fraction
!! consistency is instead maintained by mpp_lim's clamp+renormalize (required by the checker for amr with num_fluids > 1).
!!
!! Viscous (SP11): the viscous stress/work face fluxes travel through flux_src_n for the momentum and energy equations
!! (m_rhs s_compute_additional_physics_rhs: rhs += (flux_src_n(j-1) - flux_src_n(j))/dx, identical face indexing and sign to the
!! advective flux_n). They are captured into the SAME registers (added on top of the advective flux for mom..E) so the c/f reflux
!! matches the TOTAL advective+viscous flux; energy conservation therefore includes the viscous work. Fine-ghost velocity gradients
!! at the c/f boundary come from the conservative-linear cons prolongation (no special gradient reconstruction) - like the alpha
!! K-term, that inconsistency is bounded, and conservation is enforced by the flux-register matching.
!!
!! Chemistry species diffusion (SP17): the mixture-averaged species mass fluxes travel through flux_src_n for the species equations,
!! and the thermal-conduction + enthalpy energy flux through the energy equation - the same face-difference assembly as viscous.
!! They are captured into the SAME registers (species always; energy only when NOT viscous, since a viscous run already captures the
!! combined flux_src_n(E)) so the c/f reflux matches the total advective+diffusive flux and species/element/energy conservation holds
!! across the block boundary. Fine-ghost species gradients come from the species-closure cons prolongation - bounded like viscous.
module m_amr_registers

    use m_derived_types
    use m_global_parameters

    implicit none

    private; public :: s_initialize_amr_registers, s_amr_capture_boundary_flux, s_amr_apply_reflux, s_amr_zero_fine_registers, &
        & s_amr_apply_reflux_state, s_finalize_amr_registers, s_amr_reflux_face_flags, freg, creg

    !> SSP-RK3 effective flux weights: q^{n+1} = q^n + dt*(L(q^n)/6 + L(q^(1))/6 + 2*L(q^(2))/3).
    real(wp), parameter :: rk3_w(3) = [1._wp/6._wp, 1._wp/6._wp, 2._wp/3._wp]

    !> Registers for the two block faces normal to one direction: (1:sys_size, transverse-1, transverse-2, 1:amr_max_blocks). The
    !! trailing dimension is the block slot (indexed by amr_cur); each slot's registers are captured/applied independently.
    type t_face_reg
        real(wp), allocatable :: lo(:,:,:,:)
        real(wp), allocatable :: hi(:,:,:,:)
    end type t_face_reg

    type(t_face_reg) :: creg(3)  !< coarse flux at the block boundary faces (relative 0-based transverse)
    type(t_face_reg) :: freg(3)  !< fine flux at the covering fine faces (0-based fine transverse)
    $:GPU_DECLARE(create='[creg, freg]')

contains

    !> Reflux-face participation for THIS rank: own_lo(d)/own_hi(d) = it owns the coarse cell layer just OUTSIDE the block's
    !! low/high face in dim d (where the coarse capture and both reflux applies run; at an interior face the same rank also holds
    !! the inside cells) - i.e. the outside layer lies in its subdomain in dim d and the block's transverse range overlaps its
    !! subdomain. Fine-level distribution: participation is derived from the REPLICATED block range vs this rank's coarse subdomain
    !! (NOT amr_isect, which is owner-only under whole-block ownership); tlo/thi return the GLOBAL transverse overlap
    !! [max(region_lo, sidx) : min(region_hi, sidx+ext)] per dim, so capture and apply share a block-relative frame aligned with the
    !! owner's freg. All true / full-block at np=1. Also returns sidx/ext (collapsed dims pinned to 0). Reads the COARSE grid state
    !! in m/n/p.
    impure subroutine s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)

        integer, intent(out) :: sidx(3), ext(3)
        logical, intent(out) :: own_lo(3), own_hi(3)
        integer, intent(out) :: tlo(3), thi(3)
        logical              :: tv(3), tvd
        integer              :: d, t

        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        ! global transverse overlap of the block with this rank's coarse subdomain (collapsed dims pin to 0)
        do d = 1, 3
            tlo(d) = max(amr_region_lo(d), sidx(d))
            thi(d) = min(amr_region_hi(d), sidx(d) + ext(d))
        end do
        tv(1) = tlo(1) <= thi(1)
        tv(2) = (n_glb == 0) .or. tlo(2) <= thi(2)
        tv(3) = (p_glb == 0) .or. tlo(3) <= thi(3)
        own_lo = .false.; own_hi = .false.
        do d = 1, num_dims
            tvd = .true.
            do t = 1, num_dims
                if (t /= d) tvd = tvd .and. tv(t)
            end do
            own_lo(d) = tvd .and. amr_region_lo(d) - 1 >= sidx(d) .and. amr_region_lo(d) - 1 <= sidx(d) + ext(d)
            own_hi(d) = tvd .and. amr_region_hi(d) + 1 >= sidx(d) .and. amr_region_hi(d) + 1 <= sidx(d) + ext(d)
            ! max_grid_size tiling: a face shared with an adjacent sub-block is fine-fine, NOT a c/f boundary - exclude it from
            ! reflux (its outside cell is inside the neighbour block; refluxing there would corrupt that cell mid-step). The
            ! block-to-block fine-fine halo already makes the shared flux match. (No seams without tiling, so np=1/untiled: no-op.)
            if (own_lo(d) .and. f_amr_face_is_seam(d, -1)) own_lo(d) = .false.
            if (own_hi(d) .and. f_amr_face_is_seam(d, 1)) own_hi(d) = .false.
        end do

    end subroutine s_amr_reflux_face_flags

    !> True iff the current block's face on `side` (+1 high / -1 low) in dim d is shared with an ADJACENT sub-block (max_grid_size
    !! tiling) - i.e. some other block's opposite face is exactly one cell away with matching transverse extents. Such a seam is
    !! fine-fine, not a c/f boundary. Reads the replicated block list (amr_region_*_all) - no tiling means no match.
    pure logical function f_amr_face_is_seam(d, side) result(seam)

        integer, intent(in) :: d, side
        integer             :: y, t
        logical             :: match

        seam = .false.
        do y = 1, amr_num_blocks
            if (y == amr_cur) cycle
            if (side == 1) then
                if (amr_region_lo_all(d, y) /= amr_region_hi(d) + 1) cycle
            else
                if (amr_region_hi_all(d, y) /= amr_region_lo(d) - 1) cycle
            end if
            match = .true.
            do t = 1, num_dims
                if (t /= d) match = match .and. amr_region_lo_all(t, y) == amr_region_lo(t) .and. amr_region_hi_all(t, &
                    & y) == amr_region_hi(t)
            end do
            if (match) then; seam = .true.; return; end if
        end do

    end function f_amr_face_is_seam

    impure subroutine s_initialize_amr_registers(maxc_fit)

        integer, intent(in) :: maxc_fit(3)  !< amr_maxc_fit from m_amr (min-over-ranks local half-extent = max block a rank owns)
        integer             :: maxc1, maxc2, maxc3, max_f1, max_f2, max_f3

        if (.not. amr) return
        ! Registers on ALL ranks: regrid moves the block faces, so any rank can become a participant (fine cells for freg; outside
        ! face layer for creg capture + apply and for RECEIVING freg from the block owner). Fine-level distribution: freg is
        ! captured for the WHOLE block and indexed block-relative by every applier, so registers must span a whole block. The
        ! largest block a rank can own is amr_maxc_fit (the scratch-constraint cap), so size to it - matches m_amr's fine arrays
        ! and right-sizes the face registers to ~1/num_procs^(d-1) the global-half memory at scale.
        maxc1 = maxc_fit(1)
        maxc2 = 1; maxc3 = 1
        if (n_glb > 0) maxc2 = maxc_fit(2)
        if (p_glb > 0) maxc3 = maxc_fit(3)
        max_f1 = 2*maxc1 - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = 2*maxc2 - 1
        if (p_glb > 0) max_f3 = 2*maxc3 - 1
        ! creg: relative 0-based transverse (0:maxc_t-1); freg: 0-based fine (0:max_f_t).
        ! Device-resident (@:ALLOCATE): capture and both applies run as kernels; no host copies are read.
        @:ALLOCATE(creg(1)%lo(1:sys_size,0:maxc2 - 1,0:maxc3 - 1,1:amr_max_blocks), creg(1)%hi(1:sys_size,0:maxc2 - 1, &
                   & 0:maxc3 - 1,1:amr_max_blocks))
        @:ALLOCATE(freg(1)%lo(1:sys_size,0:max_f2,0:max_f3,1:amr_max_blocks), freg(1)%hi(1:sys_size,0:max_f2,0:max_f3, &
                   & 1:amr_max_blocks))
        if (n_glb > 0) then
            @:ALLOCATE(creg(2)%lo(1:sys_size,0:maxc1 - 1,0:maxc3 - 1,1:amr_max_blocks), creg(2)%hi(1:sys_size,0:maxc1 - 1, &
                       & 0:maxc3 - 1,1:amr_max_blocks))
            @:ALLOCATE(freg(2)%lo(1:sys_size,0:max_f1,0:max_f3,1:amr_max_blocks), freg(2)%hi(1:sys_size,0:max_f1,0:max_f3, &
                       & 1:amr_max_blocks))
        end if
        if (p_glb > 0) then
            @:ALLOCATE(creg(3)%lo(1:sys_size,0:maxc1 - 1,0:maxc2 - 1,1:amr_max_blocks), creg(3)%hi(1:sys_size,0:maxc1 - 1, &
                       & 0:maxc2 - 1,1:amr_max_blocks))
            @:ALLOCATE(freg(3)%lo(1:sys_size,0:max_f1,0:max_f2,1:amr_max_blocks), freg(3)%hi(1:sys_size,0:max_f1,0:max_f2, &
                       & 1:amr_max_blocks))
        end if

    end subroutine s_initialize_amr_registers

    !> Capture the c/f boundary-face fluxes for direction id from the just-finalized flux array. Runs INSIDE s_compute_rhs: coarse
    !! call (amr_in_fine_advance false, coarse globals) fills creg at the block boundary faces; fine call (flag true, globals
    !! swapped to the fine block) fills freg at fine faces -1 and m/n/p. creg uses relative 0-based transverse; freg uses 0-based
    !! fine.
    impure subroutine s_amr_capture_boundary_flux(id, flux_dir, flux_src, stage)

        integer, intent(in)            :: id
        type(vector_field), intent(in) :: flux_dir
        type(vector_field), intent(in) :: flux_src
        integer, intent(in)            :: stage
        integer                        :: eq, t1, t2, jlo, jhi, t1_lo, t1_hi, t2_lo, t2_hi, o1, o2, islot, save_cur
        integer                        :: sidx(3), ext(3), tlo(3), thi(3)
        logical                        :: own_lo(3), own_hi(3), cap_lo, cap_hi
        real(wp)                       :: coef
        logical                        :: accum

        if (.not. amr) return
        if (igr) return  ! stage-1 IGR coupling is restriction-only: the fused IGR flux kernels do not expose face fluxes to capture
        if (amr_in_fine_advance .and. .not. amr_rank_owns_block) return
        islot = amr_cur  ! working block slot (local => captured by value in the device kernels below)
        ! flux data was just written by device kernels; the face reads below run as device kernels too
        if (amr_subcycle) then
            if (amr_in_fine_advance) then
                coef = 0.5_wp*rk3_w(stage); accum = .true.  ! zeroed by s_amr_zero_fine_registers before substep 1
            else
                coef = rk3_w(stage); accum = (stage > 1)  ! stage 1 overwrites = implicit zero per coarse step
            end if
        else
            coef = 1._wp; accum = .false.  ! overwrite each stage - default behavior, byte-identical
        end if
        if (amr_in_fine_advance) then
            ! fine branch: globals swapped; jlo=-1, jhi=current fine extent in direction id
            select case (id)
            case (1); jlo = -1; jhi = m; t1_hi = n; t2_hi = p
            case (2); jlo = -1; jhi = n; t1_hi = m; t2_hi = p
            case (3); jlo = -1; jhi = p; t1_hi = m; t2_hi = n
            end select
            $:GPU_PARALLEL_LOOP(collapse=3)
            do t2 = 0, t2_hi
                do t1 = 0, t1_hi
                    do eq = 1, sys_size
                        select case (id)
                        case (1)
                            if (accum) then
                                freg(1)%lo(eq, t1, t2, islot) = freg(1)%lo(eq, t1, t2, islot) + coef*real(flux_dir%vf(eq)%sf(jlo, &
                                     & t1, t2), wp)
                                freg(1)%hi(eq, t1, t2, islot) = freg(1)%hi(eq, t1, t2, islot) + coef*real(flux_dir%vf(eq)%sf(jhi, &
                                     & t1, t2), wp)
                            else
                                freg(1)%lo(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(jlo, t1, t2), wp)
                                freg(1)%hi(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(jhi, t1, t2), wp)
                            end if
                        case (2)
                            if (accum) then
                                freg(2)%lo(eq, t1, t2, islot) = freg(2)%lo(eq, t1, t2, islot) + coef*real(flux_dir%vf(eq)%sf(t1, &
                                     & jlo, t2), wp)
                                freg(2)%hi(eq, t1, t2, islot) = freg(2)%hi(eq, t1, t2, islot) + coef*real(flux_dir%vf(eq)%sf(t1, &
                                     & jhi, t2), wp)
                            else
                                freg(2)%lo(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(t1, jlo, t2), wp)
                                freg(2)%hi(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(t1, jhi, t2), wp)
                            end if
                        case (3)
                            if (accum) then
                                freg(3)%lo(eq, t1, t2, islot) = freg(3)%lo(eq, t1, t2, islot) + coef*real(flux_dir%vf(eq)%sf(t1, &
                                     & t2, jlo), wp)
                                freg(3)%hi(eq, t1, t2, islot) = freg(3)%hi(eq, t1, t2, islot) + coef*real(flux_dir%vf(eq)%sf(t1, &
                                     & t2, jhi), wp)
                            else
                                freg(3)%lo(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(t1, t2, jlo), wp)
                                freg(3)%hi(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(t1, t2, jhi), wp)
                            end if
                        end select
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            ! total-flux matching: add the viscous momentum/energy face fluxes (flux_src) into the same fine
            ! registers so the c/f reflux sees advective+viscous. Base coef/accum are applied above; always
            ! accumulate here. Inviscid path skips this entirely (registers stay byte-identical).
            if (viscous) then
                $:GPU_PARALLEL_LOOP(collapse=3)
                do t2 = 0, t2_hi
                    do t1 = 0, t1_hi
                        do eq = eqn_idx%mom%beg, eqn_idx%E
                            select case (id)
                            case (1)
                                freg(1)%lo(eq, t1, t2, islot) = freg(1)%lo(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(jlo, &
                                     & t1, t2), wp)
                                freg(1)%hi(eq, t1, t2, islot) = freg(1)%hi(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(jhi, &
                                     & t1, t2), wp)
                            case (2)
                                freg(2)%lo(eq, t1, t2, islot) = freg(2)%lo(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(t1, &
                                     & jlo, t2), wp)
                                freg(2)%hi(eq, t1, t2, islot) = freg(2)%hi(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(t1, &
                                     & jhi, t2), wp)
                            case (3)
                                freg(3)%lo(eq, t1, t2, islot) = freg(3)%lo(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(t1, &
                                     & t2, jlo), wp)
                                freg(3)%hi(eq, t1, t2, islot) = freg(3)%hi(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(t1, &
                                     & t2, jhi), wp)
                            end select
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            ! total-flux matching (chemistry species diffusion): the mixture-averaged species mass fluxes
            ! travel through flux_src_n for the species equations; the thermal-conduction + enthalpy energy
            ! flux travels through the energy equation and is captured here only when NOT viscous (the
            ! viscous block above already captured flux_src_n(E), which holds viscous+diffusion combined).
            if (chemistry .and. chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(collapse=2)
                do t2 = 0, t2_hi
                    do t1 = 0, t1_hi
                        $:GPU_LOOP(parallelism='[seq]')
                        do eq = eqn_idx%species%beg, eqn_idx%species%end
                            select case (id)
                            case (1)
                                freg(1)%lo(eq, t1, t2, islot) = freg(1)%lo(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(jlo, &
                                     & t1, t2), wp)
                                freg(1)%hi(eq, t1, t2, islot) = freg(1)%hi(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(jhi, &
                                     & t1, t2), wp)
                            case (2)
                                freg(2)%lo(eq, t1, t2, islot) = freg(2)%lo(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(t1, &
                                     & jlo, t2), wp)
                                freg(2)%hi(eq, t1, t2, islot) = freg(2)%hi(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(t1, &
                                     & jhi, t2), wp)
                            case (3)
                                freg(3)%lo(eq, t1, t2, islot) = freg(3)%lo(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(t1, &
                                     & t2, jlo), wp)
                                freg(3)%hi(eq, t1, t2, islot) = freg(3)%hi(eq, t1, t2, islot) + coef*real(flux_src%vf(eq)%sf(t1, &
                                     & t2, jhi), wp)
                            end select
                        end do
                        if (.not. viscous) then
                            select case (id)
                            case (1)
                                freg(1)%lo(eqn_idx%E, t1, t2, islot) = freg(1)%lo(eqn_idx%E, t1, t2, &
                                     & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(jlo, t1, t2), wp)
                                freg(1)%hi(eqn_idx%E, t1, t2, islot) = freg(1)%hi(eqn_idx%E, t1, t2, &
                                     & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(jhi, t1, t2), wp)
                            case (2)
                                freg(2)%lo(eqn_idx%E, t1, t2, islot) = freg(2)%lo(eqn_idx%E, t1, t2, &
                                     & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(t1, jlo, t2), wp)
                                freg(2)%hi(eqn_idx%E, t1, t2, islot) = freg(2)%hi(eqn_idx%E, t1, t2, &
                                     & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(t1, jhi, t2), wp)
                            case (3)
                                freg(3)%lo(eqn_idx%E, t1, t2, islot) = freg(3)%lo(eqn_idx%E, t1, t2, &
                                     & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(t1, t2, jlo), wp)
                                freg(3)%hi(eqn_idx%E, t1, t2, islot) = freg(3)%hi(eqn_idx%E, t1, t2, &
                                     & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(t1, t2, jhi), wp)
                            end select
                        end if
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        else
            ! coarse branch: a face's capture runs on the rank owning the coarse cells just OUTSIDE it (its
            ! flux_n covers that face; at a rank-interior face the same rank also holds the inside cells).
            ! jlo/jhi = LOCAL flux indices of the block's low/high faces; t1/t2 = 0-based transverse indices
            ! relative to this rank's block INTERSECTION (o1/o2 = local transverse origins), aligned with the
            ! fine registers: the fine children of isect-relative cell t are faces 2*t and 2*t+1. At np=1 the
            ! intersection is the block and both flags hold, recovering the single-rank behavior exactly.
            ! ONE coarse s_compute_rhs pass fills EVERY active block's registers: revisit each slot's region+intersection in turn.
            save_cur = amr_cur
            do islot = 1, amr_num_blocks
                call s_amr_select_slot(islot)
                call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
                cap_lo = own_lo(id); cap_hi = own_hi(id)
                if (cap_lo .or. cap_hi) then
                    ! block-relative transverse frame (0-based from region_lo, aligned with the owner's freg): this rank fills
                    ! creg over its owned overlap [tlo-region_lo : thi-region_lo]; o1/o2 map that back to LOCAL flux indices.
                    select case (id)
                    case (1); jlo = amr_region_lo(1) - 1 - sidx(1); jhi = amr_region_hi(1) - sidx(1)
                        t1_lo = tlo(2) - amr_region_lo(2); t1_hi = thi(2) - amr_region_lo(2); o1 = amr_region_lo(2) - sidx(2)
                        t2_lo = tlo(3) - amr_region_lo(3); t2_hi = thi(3) - amr_region_lo(3); o2 = amr_region_lo(3) - sidx(3)
                    case (2); jlo = amr_region_lo(2) - 1 - sidx(2); jhi = amr_region_hi(2) - sidx(2)
                        t1_lo = tlo(1) - amr_region_lo(1); t1_hi = thi(1) - amr_region_lo(1); o1 = amr_region_lo(1) - sidx(1)
                        t2_lo = tlo(3) - amr_region_lo(3); t2_hi = thi(3) - amr_region_lo(3); o2 = amr_region_lo(3) - sidx(3)
                    case (3); jlo = amr_region_lo(3) - 1 - sidx(3); jhi = amr_region_hi(3) - sidx(3)
                        t1_lo = tlo(1) - amr_region_lo(1); t1_hi = thi(1) - amr_region_lo(1); o1 = amr_region_lo(1) - sidx(1)
                        t2_lo = tlo(2) - amr_region_lo(2); t2_hi = thi(2) - amr_region_lo(2); o2 = amr_region_lo(2) - sidx(2)
                    end select
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do t2 = t2_lo, t2_hi
                        do t1 = t1_lo, t1_hi
                            do eq = 1, sys_size
                                select case (id)
                                case (1)
                                    if (cap_lo) then
                                        if (accum) then
                                            creg(1)%lo(eq, t1, t2, islot) = creg(1)%lo(eq, t1, t2, &
                                                 & islot) + coef*real(flux_dir%vf(eq)%sf(jlo, o1 + t1, o2 + t2), wp)
                                        else
                                            creg(1)%lo(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(jlo, o1 + t1, o2 + t2), wp)
                                        end if
                                    end if
                                    if (cap_hi) then
                                        if (accum) then
                                            creg(1)%hi(eq, t1, t2, islot) = creg(1)%hi(eq, t1, t2, &
                                                 & islot) + coef*real(flux_dir%vf(eq)%sf(jhi, o1 + t1, o2 + t2), wp)
                                        else
                                            creg(1)%hi(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(jhi, o1 + t1, o2 + t2), wp)
                                        end if
                                    end if
                                case (2)
                                    if (cap_lo) then
                                        if (accum) then
                                            creg(2)%lo(eq, t1, t2, islot) = creg(2)%lo(eq, t1, t2, &
                                                 & islot) + coef*real(flux_dir%vf(eq)%sf(o1 + t1, jlo, o2 + t2), wp)
                                        else
                                            creg(2)%lo(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(o1 + t1, jlo, o2 + t2), wp)
                                        end if
                                    end if
                                    if (cap_hi) then
                                        if (accum) then
                                            creg(2)%hi(eq, t1, t2, islot) = creg(2)%hi(eq, t1, t2, &
                                                 & islot) + coef*real(flux_dir%vf(eq)%sf(o1 + t1, jhi, o2 + t2), wp)
                                        else
                                            creg(2)%hi(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(o1 + t1, jhi, o2 + t2), wp)
                                        end if
                                    end if
                                case (3)
                                    if (cap_lo) then
                                        if (accum) then
                                            creg(3)%lo(eq, t1, t2, islot) = creg(3)%lo(eq, t1, t2, &
                                                 & islot) + coef*real(flux_dir%vf(eq)%sf(o1 + t1, o2 + t2, jlo), wp)
                                        else
                                            creg(3)%lo(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(o1 + t1, o2 + t2, jlo), wp)
                                        end if
                                    end if
                                    if (cap_hi) then
                                        if (accum) then
                                            creg(3)%hi(eq, t1, t2, islot) = creg(3)%hi(eq, t1, t2, &
                                                 & islot) + coef*real(flux_dir%vf(eq)%sf(o1 + t1, o2 + t2, jhi), wp)
                                        else
                                            creg(3)%hi(eq, t1, t2, islot) = coef*real(flux_dir%vf(eq)%sf(o1 + t1, o2 + t2, jhi), wp)
                                        end if
                                    end if
                                end select
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    ! total-flux matching (coarse side): add viscous momentum/energy face fluxes into creg, same
                    ! face gating and transverse offsets as the base capture; always accumulate.
                    if (viscous) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do t2 = t2_lo, t2_hi
                            do t1 = t1_lo, t1_hi
                                do eq = eqn_idx%mom%beg, eqn_idx%E
                                    select case (id)
                                    case (1)
                                        if (cap_lo) creg(1)%lo(eq, t1, t2, islot) = creg(1)%lo(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(jlo, o1 + t1, o2 + t2), wp)
                                        if (cap_hi) creg(1)%hi(eq, t1, t2, islot) = creg(1)%hi(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(jhi, o1 + t1, o2 + t2), wp)
                                    case (2)
                                        if (cap_lo) creg(2)%lo(eq, t1, t2, islot) = creg(2)%lo(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(o1 + t1, jlo, o2 + t2), wp)
                                        if (cap_hi) creg(2)%hi(eq, t1, t2, islot) = creg(2)%hi(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(o1 + t1, jhi, o2 + t2), wp)
                                    case (3)
                                        if (cap_lo) creg(3)%lo(eq, t1, t2, islot) = creg(3)%lo(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(o1 + t1, o2 + t2, jlo), wp)
                                        if (cap_hi) creg(3)%hi(eq, t1, t2, islot) = creg(3)%hi(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(o1 + t1, o2 + t2, jhi), wp)
                                    end select
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                    ! total-flux matching (coarse side, chemistry species diffusion): species mass fluxes
                    ! into creg (and the energy flux only when NOT viscous, as on the fine side); same face
                    ! gating and transverse offsets as the base capture; always accumulate.
                    if (chemistry .and. chem_params%diffusion) then
                        $:GPU_PARALLEL_LOOP(collapse=2)
                        do t2 = t2_lo, t2_hi
                            do t1 = t1_lo, t1_hi
                                $:GPU_LOOP(parallelism='[seq]')
                                do eq = eqn_idx%species%beg, eqn_idx%species%end
                                    select case (id)
                                    case (1)
                                        if (cap_lo) creg(1)%lo(eq, t1, t2, islot) = creg(1)%lo(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(jlo, o1 + t1, o2 + t2), wp)
                                        if (cap_hi) creg(1)%hi(eq, t1, t2, islot) = creg(1)%hi(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(jhi, o1 + t1, o2 + t2), wp)
                                    case (2)
                                        if (cap_lo) creg(2)%lo(eq, t1, t2, islot) = creg(2)%lo(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(o1 + t1, jlo, o2 + t2), wp)
                                        if (cap_hi) creg(2)%hi(eq, t1, t2, islot) = creg(2)%hi(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(o1 + t1, jhi, o2 + t2), wp)
                                    case (3)
                                        if (cap_lo) creg(3)%lo(eq, t1, t2, islot) = creg(3)%lo(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(o1 + t1, o2 + t2, jlo), wp)
                                        if (cap_hi) creg(3)%hi(eq, t1, t2, islot) = creg(3)%hi(eq, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eq)%sf(o1 + t1, o2 + t2, jhi), wp)
                                    end select
                                end do
                                if (.not. viscous) then
                                    select case (id)
                                    case (1)
                                        if (cap_lo) creg(1)%lo(eqn_idx%E, t1, t2, islot) = creg(1)%lo(eqn_idx%E, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(jlo, o1 + t1, o2 + t2), wp)
                                        if (cap_hi) creg(1)%hi(eqn_idx%E, t1, t2, islot) = creg(1)%hi(eqn_idx%E, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(jhi, o1 + t1, o2 + t2), wp)
                                    case (2)
                                        if (cap_lo) creg(2)%lo(eqn_idx%E, t1, t2, islot) = creg(2)%lo(eqn_idx%E, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(o1 + t1, jlo, o2 + t2), wp)
                                        if (cap_hi) creg(2)%hi(eqn_idx%E, t1, t2, islot) = creg(2)%hi(eqn_idx%E, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(o1 + t1, jhi, o2 + t2), wp)
                                    case (3)
                                        if (cap_lo) creg(3)%lo(eqn_idx%E, t1, t2, islot) = creg(3)%lo(eqn_idx%E, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(o1 + t1, o2 + t2, jlo), wp)
                                        if (cap_hi) creg(3)%hi(eqn_idx%E, t1, t2, islot) = creg(3)%hi(eqn_idx%E, t1, t2, &
                                            & islot) + coef*real(flux_src%vf(eqn_idx%E)%sf(o1 + t1, o2 + t2, jhi), wp)
                                    end select
                                end if
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                end if  ! cap_lo .or. cap_hi
            end do
            call s_amr_select_slot(save_cur)
        end if

    end subroutine s_amr_capture_boundary_flux

    !> Correct the coarse rhs in the first cell OUTSIDE each block face so the coarse update sees the (child-averaged) fine flux at
    !! every c/f face. Signs follow rhs = (flux_left - flux_right)/dx: low face is the outside cell's RIGHT face => rhs += (F_coarse
    !! - Fbar_fine)/dx; high face is the outside cell's LEFT face => rhs += (Fbar_fine - F_coarse)/dx. Cells INSIDE the block need
    !! no correction (end-of-step restriction overwrites them). c1/c2 are relative 0-based coarse transverse indices.
    impure subroutine s_amr_apply_reflux(rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: eq, c1, c2, f10, f20, dd1, dd2, nch, islot
        integer                                                :: bl1, bh1, bl2, bh2, bl3, bh3, ol1, ol2, ol3, oh1, oh2, oh3
        integer                                                :: tl1, tl2, tl3, dd1_hi, dd2_hi, sidx(3), ext(3), tlo(3), thi(3)
        logical                                                :: d2, d3, own_lo(3), own_hi(3), has_lo, has_hi
        real(wp)                                               :: fblo, fbhi, mlo, mhi

        if (.not. amr) return
        if (igr) return  ! stage-1 IGR: restriction-only coupling (no captured fluxes)
        islot = amr_cur  ! working block slot (local => captured by value in the device kernels below)
        ! per-face participation: each face's correction runs on the rank owning its OUTSIDE cell layer (all faces at np=1);
        ! the owner's freg is broadcast to every participant by s_mpi_bcast_amr_reflux_faces before this is called
        call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
        if (.not. (any(own_lo) .or. any(own_hi))) return
        ! device kernels: the coarse rhs stays device-resident for the coarse RK update kernel
        d2 = n_glb > 0; d3 = p_glb > 0
        ! block-relative transverse frame (aligned with the owner's freg): loop c over each dim's owned overlap [bl:bh] =
        ! [tlo-region_lo : thi-region_lo]; creg/freg use c directly, LOCAL cell = tl + c with tl = region_lo - sidx.
        bl1 = tlo(1) - amr_region_lo(1); bh1 = thi(1) - amr_region_lo(1)
        bl2 = tlo(2) - amr_region_lo(2); bh2 = thi(2) - amr_region_lo(2)
        bl3 = tlo(3) - amr_region_lo(3); bh3 = thi(3) - amr_region_lo(3)
        ol1 = amr_region_lo(1) - 1 - sidx(1); oh1 = amr_region_hi(1) + 1 - sidx(1)
        ol2 = amr_region_lo(2) - 1 - sidx(2); oh2 = amr_region_hi(2) + 1 - sidx(2)
        ol3 = amr_region_lo(3) - 1 - sidx(3); oh3 = amr_region_hi(3) + 1 - sidx(3)
        tl1 = amr_region_lo(1) - sidx(1); tl2 = amr_region_lo(2) - sidx(2); tl3 = amr_region_lo(3) - sidx(3)
        ! x-faces: transverse dims (y, z); children in each active transverse dim
        has_lo = own_lo(1); has_hi = own_hi(1)
        if (has_lo .or. has_hi) then
            nch = 1
            if (n_glb > 0) nch = nch*2
            if (p_glb > 0) nch = nch*2
            dd1_hi = merge(1, 0, n_glb > 0); dd2_hi = merge(1, 0, p_glb > 0)
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dx(ol1)
            if (has_hi) mhi = dx(oh1)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = bl3, bh3
                    do c1 = bl2, bh2
                        f20 = 0; if (d3) f20 = 2*c2
                        f10 = 0; if (d2) f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, dd2_hi
                            do dd1 = 0, dd1_hi
                                fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20 + dd2, islot)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) rhs_vf(eq)%sf(ol1, tl2 + c1, tl3 + c2) = rhs_vf(eq)%sf(ol1, tl2 + c1, &
                            & tl3 + c2) + (creg(1)%lo(eq, c1, c2, islot) - fblo)/mlo
                        if (has_hi) rhs_vf(eq)%sf(oh1, tl2 + c1, tl3 + c2) = rhs_vf(eq)%sf(oh1, tl2 + c1, &
                            & tl3 + c2) + (fbhi - creg(1)%hi(eq, c1, c2, islot))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        ! y-faces (n_glb > 0): transverse dims (x, z); x is always active (2 children)
        has_lo = own_lo(2); has_hi = own_hi(2)
        if (n_glb > 0 .and. (has_lo .or. has_hi)) then
            nch = 2
            if (p_glb > 0) nch = nch*2
            dd2_hi = merge(1, 0, p_glb > 0)
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dy(ol2)
            if (has_hi) mhi = dy(oh2)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = bl3, bh3
                    do c1 = bl1, bh1
                        f20 = 0; if (d3) f20 = 2*c2
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, dd2_hi
                            do dd1 = 0, 1
                                fblo = fblo + freg(2)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                fbhi = fbhi + freg(2)%hi(eq, f10 + dd1, f20 + dd2, islot)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) rhs_vf(eq)%sf(tl1 + c1, ol2, tl3 + c2) = rhs_vf(eq)%sf(tl1 + c1, ol2, &
                            & tl3 + c2) + (creg(2)%lo(eq, c1, c2, islot) - fblo)/mlo
                        if (has_hi) rhs_vf(eq)%sf(tl1 + c1, oh2, tl3 + c2) = rhs_vf(eq)%sf(tl1 + c1, oh2, &
                            & tl3 + c2) + (fbhi - creg(2)%hi(eq, c1, c2, islot))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        ! z-faces (p_glb > 0): transverse dims (x, y); both always active in 3D (4 children)
        has_lo = own_lo(3); has_hi = own_hi(3)
        if (p_glb > 0 .and. (has_lo .or. has_hi)) then
            nch = 4
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dz(ol3)
            if (has_hi) mhi = dz(oh3)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = bl2, bh2
                    do c1 = bl1, bh1
                        f20 = 2*c2
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, 1
                            do dd1 = 0, 1
                                fblo = fblo + freg(3)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                fbhi = fbhi + freg(3)%hi(eq, f10 + dd1, f20 + dd2, islot)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) rhs_vf(eq)%sf(tl1 + c1, tl2 + c2, ol3) = rhs_vf(eq)%sf(tl1 + c1, tl2 + c2, &
                            & ol3) + (creg(3)%lo(eq, c1, c2, islot) - fblo)/mlo
                        if (has_hi) rhs_vf(eq)%sf(tl1 + c1, tl2 + c2, oh3) = rhs_vf(eq)%sf(tl1 + c1, tl2 + c2, &
                            & oh3) + (fbhi - creg(3)%hi(eq, c1, c2, islot))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_amr_apply_reflux

    !> Zero the fine registers (called by the subcycle driver before substep 1 - stage-1 overwrite cannot work across two substeps).
    impure subroutine s_amr_zero_fine_registers()

        integer :: d, eq, t1, t2, t1_hi, t2_hi, islot

        if (.not. amr) return
        if (igr) return  ! stage-1 IGR: restriction-only coupling (no captured fluxes)
        if (.not. amr_rank_owns_block) return
        islot = amr_cur  ! working block slot (local => captured by value in the device kernels below)
        do d = 1, 3
            if (allocated(freg(d)%lo)) then
                t1_hi = ubound(freg(d)%lo, 2); t2_hi = ubound(freg(d)%lo, 3)
                $:GPU_PARALLEL_LOOP(collapse=3)
                do t2 = 0, t2_hi
                    do t1 = 0, t1_hi
                        do eq = 1, sys_size
                            freg(d)%lo(eq, t1, t2, islot) = 0._wp
                            freg(d)%hi(eq, t1, t2, islot) = 0._wp
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end do

    end subroutine s_amr_zero_fine_registers

    !> Berger-Colella state correction (subcycle mode only): after restriction, correct the first coarse cell OUTSIDE each block
    !! face with the time-accumulated flux mismatch: low face: q += dt*(F_c_eff - Fbar_f_eff)/dx ; high face: q += dt*(Fbar_f_eff -
    !! F_c_eff)/dx. Registers hold EFFECTIVE (rk3_w-weighted, substep-averaged) fluxes in subcycle mode.
    impure subroutine s_amr_apply_reflux_state(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        integer                                                :: eq, c1, c2, f10, f20, dd1, dd2, nch, islot
        integer                                                :: bl1, bh1, bl2, bh2, bl3, bh3, ol1, ol2, ol3, oh1, oh2, oh3
        integer                                                :: tl1, tl2, tl3, dd1_hi, dd2_hi, sidx(3), ext(3), tlo(3), thi(3)
        logical                                                :: d2, d3, own_lo(3), own_hi(3), has_lo, has_hi
        real(wp)                                               :: fblo, fbhi, mlo, mhi, dtl

        if (.not. amr) return
        if (igr) return  ! stage-1 IGR: restriction-only coupling (no captured fluxes)
        islot = amr_cur  ! working block slot (local => captured by value in the device kernels below)
        ! per-face participation and block-relative index conventions: see s_amr_apply_reflux
        call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
        if (.not. (any(own_lo) .or. any(own_hi))) return
        ! device kernels: the restricted coarse state stays device-resident
        d2 = n_glb > 0; d3 = p_glb > 0
        dtl = dt
        bl1 = tlo(1) - amr_region_lo(1); bh1 = thi(1) - amr_region_lo(1)
        bl2 = tlo(2) - amr_region_lo(2); bh2 = thi(2) - amr_region_lo(2)
        bl3 = tlo(3) - amr_region_lo(3); bh3 = thi(3) - amr_region_lo(3)
        ol1 = amr_region_lo(1) - 1 - sidx(1); oh1 = amr_region_hi(1) + 1 - sidx(1)
        ol2 = amr_region_lo(2) - 1 - sidx(2); oh2 = amr_region_hi(2) + 1 - sidx(2)
        ol3 = amr_region_lo(3) - 1 - sidx(3); oh3 = amr_region_hi(3) + 1 - sidx(3)
        tl1 = amr_region_lo(1) - sidx(1); tl2 = amr_region_lo(2) - sidx(2); tl3 = amr_region_lo(3) - sidx(3)
        has_lo = own_lo(1); has_hi = own_hi(1)
        if (has_lo .or. has_hi) then
            nch = 1
            if (n_glb > 0) nch = nch*2
            if (p_glb > 0) nch = nch*2
            dd1_hi = merge(1, 0, n_glb > 0); dd2_hi = merge(1, 0, p_glb > 0)
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dx(ol1)
            if (has_hi) mhi = dx(oh1)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = bl3, bh3
                    do c1 = bl2, bh2
                        f20 = 0; if (d3) f20 = 2*c2
                        f10 = 0; if (d2) f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, dd2_hi
                            do dd1 = 0, dd1_hi
                                fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20 + dd2, islot)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) q_cons(eq)%sf(ol1, tl2 + c1, tl3 + c2) = q_cons(eq)%sf(ol1, tl2 + c1, &
                            & tl3 + c2) + dtl*(creg(1)%lo(eq, c1, c2, islot) - fblo)/mlo
                        if (has_hi) q_cons(eq)%sf(oh1, tl2 + c1, tl3 + c2) = q_cons(eq)%sf(oh1, tl2 + c1, &
                            & tl3 + c2) + dtl*(fbhi - creg(1)%hi(eq, c1, c2, islot))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        has_lo = own_lo(2); has_hi = own_hi(2)
        if (n_glb > 0 .and. (has_lo .or. has_hi)) then
            nch = 2
            if (p_glb > 0) nch = nch*2
            dd2_hi = merge(1, 0, p_glb > 0)
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dy(ol2)
            if (has_hi) mhi = dy(oh2)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = bl3, bh3
                    do c1 = bl1, bh1
                        f20 = 0; if (d3) f20 = 2*c2
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, dd2_hi
                            do dd1 = 0, 1
                                fblo = fblo + freg(2)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                fbhi = fbhi + freg(2)%hi(eq, f10 + dd1, f20 + dd2, islot)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) q_cons(eq)%sf(tl1 + c1, ol2, tl3 + c2) = q_cons(eq)%sf(tl1 + c1, ol2, &
                            & tl3 + c2) + dtl*(creg(2)%lo(eq, c1, c2, islot) - fblo)/mlo
                        if (has_hi) q_cons(eq)%sf(tl1 + c1, oh2, tl3 + c2) = q_cons(eq)%sf(tl1 + c1, oh2, &
                            & tl3 + c2) + dtl*(fbhi - creg(2)%hi(eq, c1, c2, islot))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        has_lo = own_lo(3); has_hi = own_hi(3)
        if (p_glb > 0 .and. (has_lo .or. has_hi)) then
            nch = 4
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dz(ol3)
            if (has_hi) mhi = dz(oh3)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = bl2, bh2
                    do c1 = bl1, bh1
                        f20 = 2*c2
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, 1
                            do dd1 = 0, 1
                                fblo = fblo + freg(3)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                fbhi = fbhi + freg(3)%hi(eq, f10 + dd1, f20 + dd2, islot)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) q_cons(eq)%sf(tl1 + c1, tl2 + c2, ol3) = q_cons(eq)%sf(tl1 + c1, tl2 + c2, &
                            & ol3) + dtl*(creg(3)%lo(eq, c1, c2, islot) - fblo)/mlo
                        if (has_hi) q_cons(eq)%sf(tl1 + c1, tl2 + c2, oh3) = q_cons(eq)%sf(tl1 + c1, tl2 + c2, &
                            & oh3) + dtl*(fbhi - creg(3)%hi(eq, c1, c2, islot))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_amr_apply_reflux_state

    impure subroutine s_finalize_amr_registers()

        integer :: d

        if (.not. amr) return
        do d = 1, 3
            if (allocated(creg(d)%lo)) then
                @:DEALLOCATE(creg(d)%lo, creg(d)%hi)
            end if
            if (allocated(freg(d)%lo)) then
                @:DEALLOCATE(freg(d)%lo, freg(d)%hi)
            end if
        end do

    end subroutine s_finalize_amr_registers

end module m_amr_registers
