!>
!!@file
!!@brief Contains module m_amr_registers

#:include 'macros.fpp'

!> @brief AMR flux registers: per-RK-stage refluxing at the coarse/fine block boundary (SP4). Depends only on m_derived_types +
!! m_global_parameters so both m_rhs (capture) and m_time_steppers (apply) can use it without cycles. "use m_amr" would cycle (m_amr
!! -> m_rhs -> m_amr_registers), so region info is read from amr_region_lo/hi and amr_isect_lo/hi (m_global_parameters), mirrored
!! across regrids by s_set_amr_fine_geometry. creg uses 0-based transverse indexing relative to the rank's block INTERSECTION (= the
!! block at np=1); freg uses 0-based LOCAL fine (fine children of isect cell t are 2*t and 2*t+1). All arrays are preallocated at
!! max size, so regrid needs no reallocation.
!!
!! Multi-fluid (5-eq HLLC, amr-gated path): the volume-fraction ADVECTIVE flux alpha_i*u_star travels through flux_rsx_vf (the
!! "VOLUME FRACTION FLUX" block of m_riemann_solver_hllc, same form as the mass flux), so the uniform 1:sys_size capture below
!! refluxes per-fluid masses, momentum, energy, AND alpha's advective part with no extra registers. The non-conservative remainder
!! (the +alpha*d(u_star)/dx compression term m_rhs assembles from flux_src_n = u_star) is deliberately NOT captured: alpha is
!! genuinely non-conservative, so flux-matching u_star would be wrong; coarse/fine volume-fraction consistency is instead held by
!! mpp_lim's clamp+renormalize (required by the checker for amr with num_fluids > 1).
!!
!! Viscous (SP11): the viscous stress/work face fluxes travel through flux_src_n for mom and energy (m_rhs
!! s_compute_additional_physics_rhs: rhs += (flux_src_n(j-1) - flux_src_n(j))/dx, identical face indexing and sign to advective
!! flux_rsx_vf). Captured into the SAME registers (added on top of advective flux for mom..E) so the c/f reflux matches the TOTAL
!! advective+viscous flux; energy conservation thus includes viscous work. Fine-ghost velocity gradients at the c/f boundary come
!! from the conservative-linear cons prolongation (no special gradient reconstruction) - like the alpha K-term, that inconsistency
!! is bounded, and conservation is enforced by the flux-register matching.
!!
!! Chemistry species diffusion (SP17): the mixture-averaged species mass fluxes travel through flux_src_n for the species
!! equations, and the thermal-conduction + enthalpy energy flux through the energy equation - same face-difference assembly as
!! viscous. Captured into the SAME registers (species always; energy only when NOT viscous, since a viscous run already captures
!! the combined flux_src_n(E)) so the c/f reflux matches the total advective+diffusive flux and species/element/energy conservation
!! holds across the block boundary. Fine-ghost species gradients come from the species-closure cons prolongation - bounded like
!! viscous.
module m_amr_registers

    use m_derived_types
    use m_global_parameters
    use m_riemann_state, only: flux_rsx_vf, flux_src_rsx_vf

    implicit none

    private; public :: s_initialize_amr_registers, s_amr_capture_boundary_flux, s_amr_apply_reflux, s_amr_zero_fine_registers, &
        & s_amr_apply_reflux_state, s_finalize_amr_registers, s_amr_reflux_face_flags, s_amr_reflux_apply_faces, &
        & s_amr_parent_foot, s_amr_reg_prepare, freg, creg, f_amr_face_is_seam

    !> SSP-RK3 effective flux weights: q^{n+1} = q^n + dt*(L(q^n)/6 + L(q^(1))/6 + 2*L(q^(2))/3).
    real(wp), parameter :: rk3_w(3) = [1._wp/6._wp, 1._wp/6._wp, 2._wp/3._wp]

    !> Registers for the two block faces normal to one direction: (1:sys_size, transverse-1, transverse-2, 1:amr_max_blocks). The
    !! trailing dimension is the block slot (indexed by amr_cur); each slot is captured/applied independently.
    type t_face_reg
        real(wp), allocatable :: lo(:,:,:,:)
        real(wp), allocatable :: hi(:,:,:,:)
    end type t_face_reg

    type(t_face_reg) :: creg(3)  !< coarse flux at block boundary faces (relative 0-based transverse)
    type(t_face_reg) :: freg(3)  !< fine flux at covering fine faces (0-based fine transverse)
    $:GPU_DECLARE(create='[creg, freg]')

    !> Slot capacity the registers are currently sized for, and the transverse extents they were built with.
    !!
    !! These used to be dimensioned 1:amr_max_blocks outright, which is a pure memory tax: amr_max_blocks is a SAFETY CAP, not a
    !! block count. At 400^3 with ~243 live blocks against the shipped cap of 8192 that is a 34x over-provision costing a MEASURED
    !! 1.41 MB of device memory per unused slot, i.e. 10.8 GiB per GCD. They now grow geometrically like the flat store
    !! (s_amr_st_reserve), so the cap bounds correctness only and memory follows the actual block count.
    integer, parameter :: amr_reg_floor = 64  !< initial slot capacity; growth doubles from here
    integer            :: amr_reg_cap = 0
    integer            :: rc1, rc2, rc3       !< coarse transverse extents (creg is 0:rc*-1)
    integer            :: rf1, rf2, rf3       !< fine transverse extents (freg is 0:rf*)
    !> Mesh-epoch/tripwire keys of the last participation-map build (s_amr_reg_prepare; mirror of the seam-pair cache keys).
    integer(8) :: amr_reg_epoch_built = -1_8
    integer    :: amr_reg_nblk_built = -1

    !> Per-slot geometry scratch for the batched creg capture kernels (1:amr_max_blocks): host-filled from per-slot flags/overlap,
    !! then GPU_UPDATE'd so ONE kernel iterates the slot dimension instead of O(blocks) tiny launches. bactive gates the slot;
    !! bt1lo/bt1hi/bt2lo/bt2hi are the per-slot transverse window (a slot outside the rectangular max caps is cycled); bjlo/bjhi are
    !! the normal-face flux indices; bo1/bo2 the transverse origins; bclo/bchi the per-face capture gates.
    integer, allocatable :: bjlo(:), bjhi(:), bo1(:), bo2(:), bt1lo(:), bt1hi(:), bt2lo(:), bt2hi(:)
    logical, allocatable :: bclo(:), bchi(:), bactive(:)
    $:GPU_DECLARE(create='[bjlo, bjhi, bo1, bo2, bt1lo, bt1hi, bt2lo, bt2hi, bclo, bchi, bactive]')

    !> Per-slot geometry scratch for the batched reflux APPLY kernels (mirror of the capture batching above): a_act gates the slot,
    !! a_lo/a_hi the per-face applies, a_ol/a_oh the outside coarse cell's local index in the face dim, a_t1/t2/t3 the local
    !! transverse origins, a_b1l..a_b2h the transverse windows (block-relative, freg/creg-aligned), a_mlo/a_mhi the outside-cell
    !! widths with the cyl area factors folded in. Filled per direction on the host, GPU_UPDATE'd, consumed by ONE kernel per face
    !! direction instead of O(blocks) tiny launches.
    integer, allocatable  :: a_ol(:), a_oh(:), a_t1(:), a_t2(:), a_t3(:), a_b1l(:), a_b1h(:), a_b2l(:), a_b2h(:)
    logical, allocatable  :: a_lo(:), a_hi(:), a_act(:)
    real(wp), allocatable :: a_mlo(:), a_mhi(:)
    $:GPU_DECLARE(create='[a_ol, a_oh, a_t1, a_t2, a_t3, a_b1l, a_b1h, a_b2l, a_b2h, a_lo, a_hi, a_act, a_mlo, a_mhi]')

contains

    #:def REG_GROW(A, L2, U2, L3, U3)
        if (allocated(${A}$)) then
            if (oldcap > amr_reg_grow_dev_cap) then
                ! near-limit fallback (mirror of s_amr_st_reserve): the device staging below transiently holds
                ! old + tmp = 2*oldcap slots on the device and growth fires at the memory high-water mark, so above
                ! the threshold keep the host round trip - slow, but its device peak is max(old, new).
                $:GPU_UPDATE(host='[' + A + ']')
                allocate (rtmp(1:sys_size,${L2}$:${U2}$,${L3}$:${U3}$,1:oldcap))
                rtmp = ${A}$(:,:,:,1:oldcap)
                @:DEALLOCATE(${A}$)
                @:ALLOCATE(${A}$(1:sys_size, ${L2}$:${U2}$, ${L3}$:${U3}$, 1:newcap))
                ${A}$ = 0._wp
                ${A}$(:,:,:,1:oldcap) = rtmp
                deallocate (rtmp)
                $:GPU_UPDATE(device='[' + A + ']')
            else
                ! stage the live slots on the DEVICE (rtmp_d is device-mapped by @:ALLOCATE); no PCIe traffic
                @:ALLOCATE(rtmp_d(1:sys_size, ${L2}$:${U2}$, ${L3}$:${U3}$, 1:oldcap))
                $:GPU_PARALLEL_LOOP(collapse=4)
                do c4 = 1, oldcap
                    do t2 = ${L3}$, ${U3}$
                        do t1 = ${L2}$, ${U2}$
                            do eq = 1, sys_size
                                rtmp_d(eq, t1, t2, c4) = ${A}$(eq, t1, t2, c4)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                @:DEALLOCATE(${A}$)
                @:ALLOCATE(${A}$(1:sys_size, ${L2}$:${U2}$, ${L3}$:${U3}$, 1:newcap))
                $:GPU_PARALLEL_LOOP(collapse=4)
                do c4 = 1, oldcap
                    do t2 = ${L3}$, ${U3}$
                        do t1 = ${L2}$, ${U2}$
                            do eq = 1, sys_size
                                ${A}$(eq, t1, t2, c4) = rtmp_d(eq, t1, t2, c4)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                @:DEALLOCATE(rtmp_d)
                $:GPU_PARALLEL_LOOP(collapse=4)
                do c4 = oldcap + 1, newcap
                    do t2 = ${L3}$, ${U3}$
                        do t1 = ${L2}$, ${U2}$
                            do eq = 1, sys_size
                                ${A}$(eq, t1, t2, c4) = 0._wp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if
    #:enddef

    !> Grow the reflux registers to cover at least nslot slots, doubling and never shrinking.
    !!
    !! Keyed on amr_reg_n, the DENSE participation-local count (s_amr_reg_prepare): register slots are dense indices from
    !! amr_reg_of, so capacity follows what this rank owns or participates in, not the global block count. The capture and
    !! apply kernels sweep slot = 1..amr_reg_n with bactive/a_act gating the filled subset.
    !!
    !! Contents are PRESERVED across growth, mirroring s_amr_st_reserve. Every caller of s_amr_alloc_slot today is a between-step
    !! operation (regrid, restart, slot reconcile, L0 tile build) and stage 1 overwrites the registers anyway, so discarding would
    !! probably be safe - but freg accumulates across RK stages and across subcycle substeps, so "probably" is not the right
    !! standard for a silent conservation error. Growth stages through a DEVICE temporary like the store (the registers are
    !! device-authoritative; every host consumer pulls its slot to the host immediately before reading, so the host
    !! mirror coming out of a device-path growth UNDEFINED is the store's contract, not a new one). Above the transient threshold
    !! the old host round trip remains as the OOM-safe path.
    impure subroutine s_amr_reg_reserve(nslot)

        integer, intent(in) :: nslot
        integer             :: oldcap, newcap
        integer             :: eq, t1, t2, c4
        !> device staging transiently doubles one register array's footprint; register slots are faces (~1.4 MB across all 12 arrays
        !! vs tens of MB for a store column), so the byte transient of 512 register slots matches the store's 32-column threshold.
        !! Above it, fall back to the host round trip (device peak max(old, new)).
        integer, parameter    :: amr_reg_grow_dev_cap = 512
        real(wp), allocatable :: rtmp(:,:,:,:), rtmp_d(:,:,:,:)

        if (nslot <= amr_reg_cap) return
        oldcap = amr_reg_cap
        newcap = min(amr_max_blocks, max(2*oldcap, nslot))

        @:REG_GROW(creg(1)%lo, 0, rc2 - 1, 0, rc3 - 1)
        @:REG_GROW(creg(1)%hi, 0, rc2 - 1, 0, rc3 - 1)
        @:REG_GROW(freg(1)%lo, 0, rf2, 0, rf3)
        @:REG_GROW(freg(1)%hi, 0, rf2, 0, rf3)
        @:REG_GROW(creg(2)%lo, 0, rc1 - 1, 0, rc3 - 1)
        @:REG_GROW(creg(2)%hi, 0, rc1 - 1, 0, rc3 - 1)
        @:REG_GROW(freg(2)%lo, 0, rf1, 0, rf3)
        @:REG_GROW(freg(2)%hi, 0, rf1, 0, rf3)
        @:REG_GROW(creg(3)%lo, 0, rc1 - 1, 0, rc2 - 1)
        @:REG_GROW(creg(3)%hi, 0, rc1 - 1, 0, rc2 - 1)
        @:REG_GROW(freg(3)%lo, 0, rf1, 0, rf2)
        @:REG_GROW(freg(3)%hi, 0, rf1, 0, rf2)

        amr_reg_cap = newcap

    end subroutine s_amr_reg_reserve

    !> Build/refresh the participation-local register index (amr_reg_of/amr_reg_n, m_global_parameters) and size the registers to
    !! the DENSE count. Lazily keyed on the mesh epoch (every regrid, migration, restart, and slot renumbering bumps it; a
    !! block-count change is the tripwire). A global slot g maps iff this rank (a) owns g, (b) owns g's parent (the parent-side
    !! child-creg capture, the freg receives, and the reflux-to-parent apply all index the CHILD's slot on the parent's owner), or
    !! (c) reflux-face-participates in g per s_amr_reflux_face_flags - the coarse capture, the L0/L1 apply, and the reflux face-wave
    !! receives are gated by exactly these flags, so the map cannot under-cover them. The register footprint was the O(GLOBAL boxes)
    !! device term that broke np32 weak scaling (~1.4 MB/slot across the 12 arrays, reserved toward amr_num_blocks); it is now
    !! O(owned + participation halo). The mapped range is ZEROED after a rebuild: dense slots alias across rebuilds (block g's new
    !! slot may hold another block's stale flux), and zeroing keeps the standing garbage-until-captured contract deterministic.
    !! Contents at a rebuild are dead by construction - the epoch only moves between steps, and every consumer overwrites (stage-1)
    !! or zeroes (s_amr_zero_fine_registers) before its first read of a step.
    impure subroutine s_amr_reg_prepare()

        integer :: g, kc, dch, save_cur, d, t, eq, t1, t2, t1_hi, t2_hi, islot
        integer :: sidx(3), ext(3)
        logical :: tv(3), tvd, need, is_child

        if (.not. amr) return
        if (amr_reg_epoch_built == amr_mesh_epoch .and. amr_reg_nblk_built == amr_num_blocks) return
        save_cur = amr_cur
        amr_reg_of = 0
        amr_reg_n = 0
        do g = 1, amr_num_blocks
            need = amr_owns_all(g)
            if (.not. need .and. amr_block_level(g) <= 1) then
                ! (c) reflux-face participation WITHOUT the fine-fine seam clip - the formula of
                ! f_amr_reflux_participates (m_amr) evaluated for THIS rank, keep lockstep. The subcycle p2p exchange
                ! gates its whole-slot receives on that UNCLIPPED predicate, so a rank whose only participating faces
                ! are tiling seams still posts into the block's register slot and must be mapped. The seam-clipped
                ! s_amr_reflux_face_flags fills (coarse capture, L0/L1 apply, face-wave) are a strict subset.
                call s_amr_select_slot(g)
                sidx = 0; ext = 0
                sidx(1) = start_idx(1); ext(1) = m
                if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
                if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
                tv(1) = amr_region_lo(1) <= sidx(1) + ext(1) .and. amr_region_hi(1) >= sidx(1)
                tv(2) = (n_glb == 0) .or. (amr_region_lo(2) <= sidx(2) + ext(2) .and. amr_region_hi(2) >= sidx(2))
                tv(3) = (p_glb == 0) .or. (amr_region_lo(3) <= sidx(3) + ext(3) .and. amr_region_hi(3) >= sidx(3))
                do d = 1, num_dims
                    tvd = .true.
                    do t = 1, num_dims
                        if (t /= d) tvd = tvd .and. tv(t)
                    end do
                    if (tvd .and. amr_region_lo(d) - 1 >= sidx(d) .and. amr_region_lo(d) - 1 <= sidx(d) + ext(d)) need = .true.
                    if (tvd .and. amr_region_hi(d) + 1 >= sidx(d) .and. amr_region_hi(d) + 1 <= sidx(d) + ext(d)) need = .true.
                end do
            end if
            if (need) then
                amr_reg_n = amr_reg_n + 1
                amr_reg_of(g) = amr_reg_n
            end if
        end do
        ! (b) children of owned blocks - the inline child test of the fine-branch capture below; keep lockstep
        do g = 1, amr_num_blocks
            if (.not. amr_owns_all(g)) cycle
            do kc = 1, amr_num_blocks
                if (amr_reg_of(kc) /= 0) cycle
                if (amr_block_level(kc) /= amr_block_level(g) + 1) cycle
                is_child = .true.
                do dch = 1, 3
                    is_child = is_child .and. amr_region_lo_all(dch, kc) <= amr_region_hi_all(dch, &
                        & g) .and. amr_region_hi_all(dch, kc) >= amr_region_lo_all(dch, g)
                end do
                if (.not. is_child) cycle
                amr_reg_n = amr_reg_n + 1
                amr_reg_of(kc) = amr_reg_n
            end do
        end do
        call s_amr_select_slot(save_cur)  ! also refreshes amr_reg_cur under the new map
        amr_reg_epoch_built = amr_mesh_epoch
        amr_reg_nblk_built = amr_num_blocks
        call s_amr_reg_reserve(amr_reg_n)
        do d = 1, 3
            if (allocated(freg(d)%lo) .and. amr_reg_n > 0) then
                t1_hi = ubound(freg(d)%lo, 2); t2_hi = ubound(freg(d)%lo, 3)
                $:GPU_PARALLEL_LOOP(collapse=4)
                do islot = 1, amr_reg_n
                    do t2 = 0, t2_hi
                        do t1 = 0, t1_hi
                            do eq = 1, sys_size
                                freg(d)%lo(eq, t1, t2, islot) = 0._wp
                                freg(d)%hi(eq, t1, t2, islot) = 0._wp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (allocated(creg(d)%lo) .and. amr_reg_n > 0) then
                t1_hi = ubound(creg(d)%lo, 2); t2_hi = ubound(creg(d)%lo, 3)
                $:GPU_PARALLEL_LOOP(collapse=4)
                do islot = 1, amr_reg_n
                    do t2 = 0, t2_hi
                        do t1 = 0, t1_hi
                            do eq = 1, sys_size
                                creg(d)%lo(eq, t1, t2, islot) = 0._wp
                                creg(d)%hi(eq, t1, t2, islot) = 0._wp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end do

    end subroutine s_amr_reg_prepare

    !> Reflux-face participation for THIS rank: own_lo(d)/own_hi(d) = it owns the coarse cell layer just OUTSIDE the block's
    !! low/high face in dim d (where the coarse capture and both reflux applies run; at an interior face the same rank also holds
    !! the inside cells) - i.e. the outside layer lies in its subdomain in dim d and the block's transverse range overlaps it.
    !! Fine-level distribution: participation derives from the REPLICATED block range vs this rank's coarse subdomain (NOT
    !! amr_isect, which is owner-only under whole-block ownership); tlo/thi return the GLOBAL transverse overlap [max(region_lo,
    !! sidx) : min(region_hi, sidx+ext)] per dim, so capture and apply share a block-relative frame aligned with the owner's freg.
    !! All true / full-block at np=1. Also returns sidx/ext (collapsed dims pinned to 0). Reads the COARSE grid m/n/p.
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
        ! global transverse overlap of block with this rank's coarse subdomain (collapsed dims pin to 0)
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
            ! block-to-block fine-fine halo already matches the shared flux. (No seams without tiling, so np=1/untiled: no-op.)
            if (own_lo(d) .and. f_amr_face_is_seam(d, -1)) own_lo(d) = .false.
            if (own_hi(d) .and. f_amr_face_is_seam(d, 1)) own_hi(d) = .false.
        end do

    end subroutine s_amr_reflux_face_flags

    !> True iff the current block's face on `side` (+1 high / -1 low) in dim d is shared with an ADJACENT sub-block (max_grid_size
    !! tiling) - i.e. another block's opposite face is exactly one cell away with matching transverse extents. Such a seam is
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
        ! Registers on ALL ranks: regrid moves block faces, so any rank can participate (fine cells for freg; outside-face layer
        ! for creg capture/apply and for receiving freg from the block owner). Fine-level distribution: freg is captured for the
        ! WHOLE block and indexed block-relative by every applier, so registers must span a whole block. The largest block a rank
        ! can own is amr_maxc_fit (the scratch-constraint cap), so size to it - matches m_amr's fine arrays and right-sizes the face
        ! registers to ~1/num_procs^(d-1) the global-half memory at scale.
        maxc1 = maxc_fit(1)
        maxc2 = 1; maxc3 = 1
        if (n_glb > 0) maxc2 = maxc_fit(2)
        if (p_glb > 0) maxc3 = maxc_fit(3)
        max_f1 = amr_ref_ratio*maxc1 - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = amr_ref_ratio*maxc2 - 1
        if (p_glb > 0) max_f3 = amr_ref_ratio*maxc3 - 1
        ! creg: relative 0-based transverse (0:maxc_t-1); freg: 0-based fine (0:max_f_t).
        ! Device-resident (@:ALLOCATE): capture and both applies run as kernels; no host copies read.
        ! Stash the transverse extents so s_amr_reg_reserve can rebuild the same shapes when the slot dimension grows.
        rc1 = maxc1; rc2 = maxc2; rc3 = maxc3
        rf1 = max_f1; rf2 = max_f2; rf3 = max_f3
        ! Start at a small slot capacity and grow on demand; do NOT size to amr_max_blocks (see amr_reg_cap above).
        amr_reg_cap = min(amr_max_blocks, amr_reg_floor)
        @:ALLOCATE(creg(1)%lo(1:sys_size,0:rc2 - 1,0:rc3 - 1,1:amr_reg_cap), creg(1)%hi(1:sys_size,0:rc2 - 1, 0:rc3 - 1, &
                   & 1:amr_reg_cap))
        @:ALLOCATE(freg(1)%lo(1:sys_size,0:rf2,0:rf3,1:amr_reg_cap), freg(1)%hi(1:sys_size,0:rf2,0:rf3, 1:amr_reg_cap))
        if (n_glb > 0) then
            @:ALLOCATE(creg(2)%lo(1:sys_size,0:rc1 - 1,0:rc3 - 1,1:amr_reg_cap), creg(2)%hi(1:sys_size,0:rc1 - 1, 0:rc3 - 1, &
                       & 1:amr_reg_cap))
            @:ALLOCATE(freg(2)%lo(1:sys_size,0:rf1,0:rf3,1:amr_reg_cap), freg(2)%hi(1:sys_size,0:rf1,0:rf3, 1:amr_reg_cap))
        end if
        if (p_glb > 0) then
            @:ALLOCATE(creg(3)%lo(1:sys_size,0:rc1 - 1,0:rc2 - 1,1:amr_reg_cap), creg(3)%hi(1:sys_size,0:rc1 - 1, 0:rc2 - 1, &
                       & 1:amr_reg_cap))
            @:ALLOCATE(freg(3)%lo(1:sys_size,0:rf1,0:rf2,1:amr_reg_cap), freg(3)%hi(1:sys_size,0:rf1,0:rf2, 1:amr_reg_cap))
        end if
        ! per-slot geometry scratch for the batched capture kernels (device-resident: host-filled, GPU_UPDATE'd before each call)
        @:ALLOCATE(bjlo(1:amr_max_blocks), bjhi(1:amr_max_blocks), bo1(1:amr_max_blocks), bo2(1:amr_max_blocks))
        @:ALLOCATE(bt1lo(1:amr_max_blocks), bt1hi(1:amr_max_blocks), bt2lo(1:amr_max_blocks), bt2hi(1:amr_max_blocks))
        @:ALLOCATE(bclo(1:amr_max_blocks), bchi(1:amr_max_blocks), bactive(1:amr_max_blocks))
        @:ALLOCATE(a_ol(1:amr_max_blocks), a_oh(1:amr_max_blocks), a_t1(1:amr_max_blocks), a_t2(1:amr_max_blocks))
        @:ALLOCATE(a_t3(1:amr_max_blocks), a_b1l(1:amr_max_blocks), a_b1h(1:amr_max_blocks), a_b2l(1:amr_max_blocks))
        @:ALLOCATE(a_b2h(1:amr_max_blocks), a_lo(1:amr_max_blocks), a_hi(1:amr_max_blocks), a_act(1:amr_max_blocks))
        @:ALLOCATE(a_mlo(1:amr_max_blocks), a_mhi(1:amr_max_blocks))
        ! participation-local register index (host-only ints; the register REALS are what the dense map shrinks)
        allocate (amr_reg_of(1:amr_max_blocks))
        amr_reg_of = 0; amr_reg_n = 0; amr_reg_cur = 0
        amr_reg_epoch_built = -1_8; amr_reg_nblk_built = -1

    end subroutine s_initialize_amr_registers

    !> Parent-fine footprint of block k inside its parent pblk, from REPLICATED metadata only, so every rank computes the same box
    !! (a rank needs it for a block it does NOT own, whose own amr_isect_lo/hi is the empty non-owner footprint). Mirrors the
    !! level>=2 branch of s_set_amr_fine_geometry exactly; rr is the global amr_ref_ratio because a level>=2 block's parent is never
    !! an L0 tile (the only slot with a per-slot ratio of 1). Lives here rather than in m_amr so the child-creg capture below and
    !! m_amr's P2P gather/restrict/reflux share one copy of the formula ("use m_amr" would cycle).
    pure subroutine s_amr_parent_foot(k, pblk, plo, phi)

        integer, intent(in)  :: k, pblk
        integer, intent(out) :: plo(3), phi(3)
        integer              :: d, rr

        rr = amr_ref_ratio
        do d = 1, 3
            plo(d) = rr*(amr_region_lo_all(d, k) - amr_region_lo_all(d, pblk))
            phi(d) = rr*(amr_region_hi_all(d, k) - amr_region_lo_all(d, pblk)) + (rr - 1)
        end do
        if (n_glb == 0) then; plo(2) = 0; phi(2) = 0; end if
        if (p_glb == 0) then; plo(3) = 0; phi(3) = 0; end if

    end subroutine s_amr_parent_foot

    !> Shared creg boundary-flux capture (dense eq range), BATCHED over the slot dimension: for each active slot in [1:nb],
    !! creg(id)%lo/hi(eq, t1, t2, slot) [+=/=] cf * flux(face, bo1(slot)+t1, bo2(slot)+t2) for eq in [eqb:eqe], over the per-slot
    !! transverse window [bt1lo:bt1hi] x [bt2lo:bt2hi]. acc=.true. accumulates, .false. overwrites (the merge picks the old value or
    !! 0 with no arithmetic, so a stage-1 overwrite reads no uninitialized creg). bclo/bchi gate the low/high face (unowned coarse
    !! faces off; child faces always on). The device kernel collapses (slot, t2, t1, eq) over the rectangular caps
    !! [0:maxt2]x[0:maxt1] (max over slots) and cycles inactive slots / out-of-window cells - one launch replaces O(blocks) per-slot
    !! launches. Per-slot geometry (bjlo etc.) is host-filled and GPU_UPDATE'd by the caller. Used for the advective (flat=T,
    !! eqb=1..sys_size) and viscous (flux_src, eqb=mom..E) captures on BOTH the coarse-self and child sides.
    impure subroutine s_amr_capture_creg_dense_batch(nb, id, advective, cf, acc, maxt1, maxt2, eqb, eqe)

        integer, intent(in) :: nb, id, maxt1, maxt2, eqb, eqe
        !> Which flat Riemann buffer to read: T = flux_rsx_vf (advective), F = flux_src_rsx_vf (viscous). Both are plain module
        !! arrays now, so this routine takes NO field dummies at all - the whole point of the flattening.
        logical, intent(in)  :: advective
        real(wp), intent(in) :: cf
        logical, intent(in)  :: acc
        integer              :: eq, t1, t2, slot, i1, i2, i3, j1, j2, j3
        real(wp)             :: v_lo, v_hi

        $:GPU_PARALLEL_LOOP(collapse=4, private='[i1, i2, i3, j1, j2, j3, v_lo, v_hi]')
        do slot = 1, nb
            do t2 = 0, maxt2
                do t1 = 0, maxt1
                    do eq = eqb, eqe
                        if (.not. bactive(slot)) cycle
                        if (t1 < bt1lo(slot) .or. t1 > bt1hi(slot) .or. t2 < bt2lo(slot) .or. t2 > bt2hi(slot)) cycle
                        select case (id)
                        case (1)
                            i1 = bjlo(slot); i2 = bo1(slot) + t1; i3 = bo2(slot) + t2
                            j1 = bjhi(slot); j2 = i2; j3 = i3
                        case (2)
                            i1 = bo1(slot) + t1; i2 = bjlo(slot); i3 = bo2(slot) + t2
                            j1 = i1; j2 = bjhi(slot); j3 = i3
                        case (3)
                            i1 = bo1(slot) + t1; i2 = bo2(slot) + t2; i3 = bjlo(slot)
                            j1 = i1; j2 = i2; j3 = bjhi(slot)
                        end select
                        ! The flux reads MUST stay inside the bclo/bchi guards. A slot goes active when EITHER face is owned
                        ! (s_amr_capture_boundary_flux: cap_lo .or. cap_hi), and the UNOWNED face's index is still computed - it
                        ! then points a whole block width outside this rank's subdomain (jlo down to -amr_max_grid_size). Reading
                        ! it unguarded is an out-of-bounds device access against flux_rsx_vf's tight (-1:m_alloc) bounds. It hides
                        ! at np=1, where the intersection IS the block and both flags hold, so goldens do not catch it.
                        if (bclo(slot)) then
                            if (advective) then
                                v_lo = flux_rsx_vf(i1, i2, i3, eq)
                            else
                                v_lo = flux_src_rsx_vf(i1, i2, i3, eq)
                            end if
                            select case (id)
                            case (1); creg(1)%lo(eq, t1, t2, slot) = merge(creg(1)%lo(eq, t1, t2, slot), 0._wp, acc) + cf*v_lo
                            case (2); creg(2)%lo(eq, t1, t2, slot) = merge(creg(2)%lo(eq, t1, t2, slot), 0._wp, acc) + cf*v_lo
                            case (3); creg(3)%lo(eq, t1, t2, slot) = merge(creg(3)%lo(eq, t1, t2, slot), 0._wp, acc) + cf*v_lo
                            end select
                        end if
                        if (bchi(slot)) then
                            if (advective) then
                                v_hi = flux_rsx_vf(j1, j2, j3, eq)
                            else
                                v_hi = flux_src_rsx_vf(j1, j2, j3, eq)
                            end if
                            select case (id)
                            case (1); creg(1)%hi(eq, t1, t2, slot) = merge(creg(1)%hi(eq, t1, t2, slot), 0._wp, acc) + cf*v_hi
                            case (2); creg(2)%hi(eq, t1, t2, slot) = merge(creg(2)%hi(eq, t1, t2, slot), 0._wp, acc) + cf*v_hi
                            case (3); creg(3)%hi(eq, t1, t2, slot) = merge(creg(3)%hi(eq, t1, t2, slot), 0._wp, acc) + cf*v_hi
                            end select
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_capture_creg_dense_batch

    !> Shared creg boundary-flux capture (chemistry species diffusion), BATCHED over the slot dimension: always-accumulate the
    !! species mass fluxes, plus the energy flux only when NOT viscous (the viscous pass already captured flux_src(E)). Species use
    !! a seq inner loop (runtime range). The device kernel collapses (slot, t2, t1) over the rectangular caps [0:maxt2]x[0:maxt1]
    !! (max over slots) and cycles inactive slots / out-of-window cells. Per-slot geometry is host-filled + GPU_UPDATE'd by the
    !! caller. Used for the chem capture on BOTH the coarse-self and child sides. TWIN of the chemistry freg capture in
    !! s_amr_capture_boundary_flux (fine branch): same species-always + energy-only-when-not-viscous policy - keep lockstep.
    impure subroutine s_amr_capture_creg_chem_batch(nb, id, cf, maxt1, maxt2)

        integer, intent(in)  :: nb, id, maxt1, maxt2
        real(wp), intent(in) :: cf
        integer              :: eq, t1, t2, slot

        $:GPU_PARALLEL_LOOP(collapse=3)
        do slot = 1, nb
            do t2 = 0, maxt2
                do t1 = 0, maxt1
                    if (.not. bactive(slot)) cycle
                    if (t1 < bt1lo(slot) .or. t1 > bt1hi(slot) .or. t2 < bt2lo(slot) .or. t2 > bt2hi(slot)) cycle
                    $:GPU_LOOP(parallelism='[seq]')
                    do eq = eqn_idx%species%beg, eqn_idx%species%end
                        select case (id)
                        case (1)
                            if (bclo(slot)) creg(1)%lo(eq, t1, t2, slot) = creg(1)%lo(eq, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bjlo(slot), bo1(slot) + t1, bo2(slot) + t2, eq)
                            if (bchi(slot)) creg(1)%hi(eq, t1, t2, slot) = creg(1)%hi(eq, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bjhi(slot), bo1(slot) + t1, bo2(slot) + t2, eq)
                        case (2)
                            if (bclo(slot)) creg(2)%lo(eq, t1, t2, slot) = creg(2)%lo(eq, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bo1(slot) + t1, bjlo(slot), bo2(slot) + t2, eq)
                            if (bchi(slot)) creg(2)%hi(eq, t1, t2, slot) = creg(2)%hi(eq, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bo1(slot) + t1, bjhi(slot), bo2(slot) + t2, eq)
                        case (3)
                            if (bclo(slot)) creg(3)%lo(eq, t1, t2, slot) = creg(3)%lo(eq, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bo1(slot) + t1, bo2(slot) + t2, bjlo(slot), eq)
                            if (bchi(slot)) creg(3)%hi(eq, t1, t2, slot) = creg(3)%hi(eq, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bo1(slot) + t1, bo2(slot) + t2, bjhi(slot), eq)
                        end select
                    end do
                    if (.not. viscous) then
                        select case (id)
                        case (1)
                            if (bclo(slot)) creg(1)%lo(eqn_idx%E, t1, t2, slot) = creg(1)%lo(eqn_idx%E, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bjlo(slot), bo1(slot) + t1, bo2(slot) + t2, eqn_idx%E)
                            if (bchi(slot)) creg(1)%hi(eqn_idx%E, t1, t2, slot) = creg(1)%hi(eqn_idx%E, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bjhi(slot), bo1(slot) + t1, bo2(slot) + t2, eqn_idx%E)
                        case (2)
                            if (bclo(slot)) creg(2)%lo(eqn_idx%E, t1, t2, slot) = creg(2)%lo(eqn_idx%E, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bo1(slot) + t1, bjlo(slot), bo2(slot) + t2, eqn_idx%E)
                            if (bchi(slot)) creg(2)%hi(eqn_idx%E, t1, t2, slot) = creg(2)%hi(eqn_idx%E, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bo1(slot) + t1, bjhi(slot), bo2(slot) + t2, eqn_idx%E)
                        case (3)
                            if (bclo(slot)) creg(3)%lo(eqn_idx%E, t1, t2, slot) = creg(3)%lo(eqn_idx%E, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bo1(slot) + t1, bo2(slot) + t2, bjlo(slot), eqn_idx%E)
                            if (bchi(slot)) creg(3)%hi(eqn_idx%E, t1, t2, slot) = creg(3)%hi(eqn_idx%E, t1, t2, &
                                & slot) + cf*flux_src_rsx_vf(bo1(slot) + t1, bo2(slot) + t2, bjhi(slot), eqn_idx%E)
                        end select
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_capture_creg_chem_batch

    !> Capture the c/f boundary-face fluxes for direction id from the just-finalized flux array. Runs INSIDE s_compute_rhs: coarse
    !! call (amr_in_fine_advance false, coarse globals) fills creg at block boundary faces; fine call (flag true, globals swapped to
    !! the fine block) fills freg at fine faces -1 and m/n/p. creg uses relative 0-based transverse; freg uses 0-based fine.
    impure subroutine s_amr_capture_boundary_flux(id, stage)

        integer, intent(in) :: id
        integer, intent(in) :: stage
        integer             :: eq, t1, t2, jlo, jhi, t1_lo, t1_hi, t2_lo, t2_hi, o1, o2, islot, save_cur, sreg
        integer             :: sidx(3), ext(3), tlo(3), thi(3), cflo(3), cfhi(3), kc, dch, maxt1, maxt2
        logical             :: own_lo(3), own_hi(3), cap_lo, cap_hi
        real(wp)            :: coef, ccoef
        logical             :: accum, cacc, is_child

        if (.not. amr) return
        ! Refresh the participation map + register capacity on a topology change; no-op (two integer compares) otherwise.
        call s_amr_reg_prepare()
        if (igr) return  ! stage-1 IGR coupling is restriction-only: the fused IGR flux kernels do not expose face fluxes to capture
        if (amr_in_fine_advance .and. .not. amr_rank_owns_block) return
        ! a level-0 L0 tile advancing through the fine path is COARSE, not a fine block: skip the freg self-capture and the
        ! parent-of-level-1 child-creg loop (which would overwrite the real fine block's creg in the tile-swapped frame). Its creg
        ! comes from the dedicated L0 coarse RHS (amr_in_fine_advance=F). Pure-AMR has no level-0 slots so this never fires.
        if (amr_in_fine_advance .and. amr_block_level(amr_cur) == 0) return
        islot = amr_reg_cur  ! working block's DENSE register slot (local => captured by value in the device kernels below)
        ! flux data was just written by device kernels; the face reads below run as device kernels too
        if (amr_subcycle) then
            if (amr_in_fine_advance) then
                coef = 0.5_wp*rk3_w(stage); accum = .true.  ! zeroed by s_amr_zero_fine_registers before substep 1
            else
                coef = rk3_w(stage); accum = (stage > 1)  ! stage 1 overwrites = implicit zero per coarse step
            end if
        else if (amr_in_fine_advance .and. amr_block_level(amr_cur) >= 2) then
            ! lock-step L2->L1 reflux: parent is already RK-updated by reflux time, so freg must hold the rk3_w-weighted
            ! step-integral flux for the once-per-step STATE correction (stage 1 overwrites = implicit zero, cf. coarse creg).
            coef = rk3_w(stage); accum = (stage > 1)
        else
            coef = 1._wp; accum = .false.  ! overwrite each stage - default, byte-identical
        end if
        if (amr_in_fine_advance) then
            ! fine branch: globals swapped; jlo=-1, jhi=current fine extent in direction id.
            ! TWIN of the creg capture: the advective / viscous (flux_src mom..E) / chemistry (flux_src species always, energy only
            ! when NOT viscous) captures below stay lockstep with s_amr_capture_creg_dense_batch + s_amr_capture_creg_chem_batch,
            ! which encode the identical policy on the coarse side. The "energy only when not viscous" rule lives in FOUR places -
            ! here (freg viscous + chemistry blocks) and both creg batch helpers - change one, change all, or the c/f reflux
            ! subtracts mismatched coarse/fine fluxes (a conservation leak no single-level test catches).
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
                                freg(1)%lo(eq, t1, t2, islot) = freg(1)%lo(eq, t1, t2, islot) + coef*flux_rsx_vf(jlo, t1, t2, eq)
                                freg(1)%hi(eq, t1, t2, islot) = freg(1)%hi(eq, t1, t2, islot) + coef*flux_rsx_vf(jhi, t1, t2, eq)
                            else
                                freg(1)%lo(eq, t1, t2, islot) = coef*flux_rsx_vf(jlo, t1, t2, eq)
                                freg(1)%hi(eq, t1, t2, islot) = coef*flux_rsx_vf(jhi, t1, t2, eq)
                            end if
                        case (2)
                            if (accum) then
                                freg(2)%lo(eq, t1, t2, islot) = freg(2)%lo(eq, t1, t2, islot) + coef*flux_rsx_vf(t1, jlo, t2, eq)
                                freg(2)%hi(eq, t1, t2, islot) = freg(2)%hi(eq, t1, t2, islot) + coef*flux_rsx_vf(t1, jhi, t2, eq)
                            else
                                freg(2)%lo(eq, t1, t2, islot) = coef*flux_rsx_vf(t1, jlo, t2, eq)
                                freg(2)%hi(eq, t1, t2, islot) = coef*flux_rsx_vf(t1, jhi, t2, eq)
                            end if
                        case (3)
                            if (accum) then
                                freg(3)%lo(eq, t1, t2, islot) = freg(3)%lo(eq, t1, t2, islot) + coef*flux_rsx_vf(t1, t2, jlo, eq)
                                freg(3)%hi(eq, t1, t2, islot) = freg(3)%hi(eq, t1, t2, islot) + coef*flux_rsx_vf(t1, t2, jhi, eq)
                            else
                                freg(3)%lo(eq, t1, t2, islot) = coef*flux_rsx_vf(t1, t2, jlo, eq)
                                freg(3)%hi(eq, t1, t2, islot) = coef*flux_rsx_vf(t1, t2, jhi, eq)
                            end if
                        end select
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            ! total-flux matching: add the viscous mom/energy face fluxes (flux_src) into the same fine registers so the c/f reflux
            ! sees advective+viscous. Base coef applied above; always accumulate here. Inviscid path skips this (registers stay
            ! byte-identical).
            if (viscous) then
                $:GPU_PARALLEL_LOOP(collapse=3)
                do t2 = 0, t2_hi
                    do t1 = 0, t1_hi
                        do eq = eqn_idx%mom%beg, eqn_idx%E
                            select case (id)
                            case (1)
                                freg(1)%lo(eq, t1, t2, islot) = freg(1)%lo(eq, t1, t2, islot) + coef*flux_src_rsx_vf(jlo, t1, t2, &
                                     & eq)
                                freg(1)%hi(eq, t1, t2, islot) = freg(1)%hi(eq, t1, t2, islot) + coef*flux_src_rsx_vf(jhi, t1, t2, &
                                     & eq)
                            case (2)
                                freg(2)%lo(eq, t1, t2, islot) = freg(2)%lo(eq, t1, t2, islot) + coef*flux_src_rsx_vf(t1, jlo, t2, &
                                     & eq)
                                freg(2)%hi(eq, t1, t2, islot) = freg(2)%hi(eq, t1, t2, islot) + coef*flux_src_rsx_vf(t1, jhi, t2, &
                                     & eq)
                            case (3)
                                freg(3)%lo(eq, t1, t2, islot) = freg(3)%lo(eq, t1, t2, islot) + coef*flux_src_rsx_vf(t1, t2, jlo, &
                                     & eq)
                                freg(3)%hi(eq, t1, t2, islot) = freg(3)%hi(eq, t1, t2, islot) + coef*flux_src_rsx_vf(t1, t2, jhi, &
                                     & eq)
                            end select
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            ! total-flux matching (chemistry species diffusion): the mixture-averaged species mass fluxes travel through
            ! flux_src_rsx_vf
            ! for the species equations; the thermal-conduction + enthalpy energy flux travels through the energy equation, captured
            ! here only when NOT viscous (the viscous block above already captured flux_src_rsx_vf(E), which holds
            ! viscous+diffusion).
            if (chemistry .and. chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(collapse=2)
                do t2 = 0, t2_hi
                    do t1 = 0, t1_hi
                        $:GPU_LOOP(parallelism='[seq]')
                        do eq = eqn_idx%species%beg, eqn_idx%species%end
                            select case (id)
                            case (1)
                                freg(1)%lo(eq, t1, t2, islot) = freg(1)%lo(eq, t1, t2, islot) + coef*flux_src_rsx_vf(jlo, t1, t2, &
                                     & eq)
                                freg(1)%hi(eq, t1, t2, islot) = freg(1)%hi(eq, t1, t2, islot) + coef*flux_src_rsx_vf(jhi, t1, t2, &
                                     & eq)
                            case (2)
                                freg(2)%lo(eq, t1, t2, islot) = freg(2)%lo(eq, t1, t2, islot) + coef*flux_src_rsx_vf(t1, jlo, t2, &
                                     & eq)
                                freg(2)%hi(eq, t1, t2, islot) = freg(2)%hi(eq, t1, t2, islot) + coef*flux_src_rsx_vf(t1, jhi, t2, &
                                     & eq)
                            case (3)
                                freg(3)%lo(eq, t1, t2, islot) = freg(3)%lo(eq, t1, t2, islot) + coef*flux_src_rsx_vf(t1, t2, jlo, &
                                     & eq)
                                freg(3)%hi(eq, t1, t2, islot) = freg(3)%hi(eq, t1, t2, islot) + coef*flux_src_rsx_vf(t1, t2, jhi, &
                                     & eq)
                            end select
                        end do
                        if (.not. viscous) then
                            select case (id)
                            case (1)
                                freg(1)%lo(eqn_idx%E, t1, t2, islot) = freg(1)%lo(eqn_idx%E, t1, t2, &
                                     & islot) + coef*flux_src_rsx_vf(jlo, t1, t2, eqn_idx%E)
                                freg(1)%hi(eqn_idx%E, t1, t2, islot) = freg(1)%hi(eqn_idx%E, t1, t2, &
                                     & islot) + coef*flux_src_rsx_vf(jhi, t1, t2, eqn_idx%E)
                            case (2)
                                freg(2)%lo(eqn_idx%E, t1, t2, islot) = freg(2)%lo(eqn_idx%E, t1, t2, &
                                     & islot) + coef*flux_src_rsx_vf(t1, jlo, t2, eqn_idx%E)
                                freg(2)%hi(eqn_idx%E, t1, t2, islot) = freg(2)%hi(eqn_idx%E, t1, t2, &
                                     & islot) + coef*flux_src_rsx_vf(t1, jhi, t2, eqn_idx%E)
                            case (3)
                                freg(3)%lo(eqn_idx%E, t1, t2, islot) = freg(3)%lo(eqn_idx%E, t1, t2, &
                                     & islot) + coef*flux_src_rsx_vf(t1, t2, jlo, eqn_idx%E)
                                freg(3)%hi(eqn_idx%E, t1, t2, islot) = freg(3)%hi(eqn_idx%E, t1, t2, &
                                     & islot) + coef*flux_src_rsx_vf(t1, t2, jhi, eqn_idx%E)
                            end select
                        end if
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            ! multi-level lock-step: this fine block (amr_cur) is the COARSE side (parent) of its level+1 children. Capture creg for
            ! each child from THIS block's fine flux at the child's footprint faces - the child's amr_isect_lo/hi is already in this
            ! parent's fine frame, so it indexes flux_rsx_vf directly (face jlo=isect_lo-1, jhi=isect_hi; transverse origin o1/o2).
            ! creg holds the rk3_w-weighted step-integral flux for the once-per-step STATE reflux into this parent
            ! (s_amr_reflux_to_parent). Captures the TOTAL flux - advective (flux_rsx_vf), then viscous (flux_src, mom..E), then
            ! chemistry species+energy - mirroring the coarse-self branch below, so viscous/chemistry multi-level conserves (no
            ! checker gate). creg is the PARENT's OWN flux, so the parent owner captures it for EVERY child of this block -
            ! including children owned by another rank, which supply only the matching freg (s_amr_p2p_freg_to_parent). Framing
            ! therefore comes from s_amr_parent_foot (replicated metadata), NOT amr_isect_*_all(:,kc), which is the empty sentinel
            ! for a child this rank does not own. Under tower co-location every child IS owned, so this captures the identical set.
            ! Fill per-slot (per-child) geometry, then one batched kernel per capture category. Each child's creg lives at its
            ! DENSE register slot (sreg = amr_reg_of(kc); a child of an owned block is always mapped - s_amr_reg_prepare clause
            ! (b) is this loop's twin); both faces always owned (the parent spans the whole child footprint), t1lo=t2lo=0.
            ccoef = rk3_w(stage); cacc = (stage > 1)
            bactive = .false.
            maxt1 = 0; maxt2 = 0
            do kc = 1, amr_num_blocks
                if (amr_block_level(kc) /= amr_block_level(amr_cur) + 1) cycle
                is_child = .true.
                do dch = 1, 3
                    is_child = is_child .and. amr_region_lo_all(dch, kc) <= amr_region_hi_all(dch, &
                        & amr_cur) .and. amr_region_hi_all(dch, kc) >= amr_region_lo_all(dch, amr_cur)
                end do
                if (.not. is_child) cycle
                call s_amr_parent_foot(kc, amr_cur, cflo, cfhi)
                select case (id)
                case (1); jlo = cflo(1) - 1; jhi = cfhi(1)
                    o1 = cflo(2); t1_hi = cfhi(2) - cflo(2)
                    o2 = cflo(3); t2_hi = cfhi(3) - cflo(3)
                case (2); jlo = cflo(2) - 1; jhi = cfhi(2)
                    o1 = cflo(1); t1_hi = cfhi(1) - cflo(1)
                    o2 = cflo(3); t2_hi = cfhi(3) - cflo(3)
                case (3); jlo = cflo(3) - 1; jhi = cfhi(3)
                    o1 = cflo(1); t1_hi = cfhi(1) - cflo(1)
                    o2 = cflo(2); t2_hi = cfhi(2) - cflo(2)
                end select
                sreg = amr_reg_of(kc)
                bactive(sreg) = .true.; bclo(sreg) = .true.; bchi(sreg) = .true.
                bjlo(sreg) = jlo; bjhi(sreg) = jhi; bo1(sreg) = o1; bo2(sreg) = o2
                bt1lo(sreg) = 0; bt1hi(sreg) = t1_hi; bt2lo(sreg) = 0; bt2hi(sreg) = t2_hi
                maxt1 = max(maxt1, t1_hi); maxt2 = max(maxt2, t2_hi)
            end do
            if (any(bactive(1:amr_reg_n))) then
                $:GPU_UPDATE(device='[bjlo, bjhi, bo1, bo2, bt1lo, bt1hi, bt2lo, bt2hi, bclo, bchi, bactive]')
                ! shared capture into each CHILD's creg (parent-fine frame): advective, then total-flux viscous, then chemistry
                call s_amr_capture_creg_dense_batch(amr_reg_n, id, .true., ccoef, cacc, maxt1, maxt2, 1, sys_size)
                if (viscous) call s_amr_capture_creg_dense_batch(amr_reg_n, id, .false., ccoef, .true., maxt1, maxt2, &
                    & eqn_idx%mom%beg, eqn_idx%E)
                if (chemistry .and. chem_params%diffusion) call s_amr_capture_creg_chem_batch(amr_reg_n, id, ccoef, maxt1, maxt2)
            end if
        else
            ! coarse branch: a face's capture runs on the rank owning the coarse cells just OUTSIDE it (its flux_rsx_vf covers that
            ! face;
            ! at a rank-interior face the same rank also holds the inside cells). jlo/jhi = LOCAL flux indices of the block's
            ! low/high faces; t1/t2 = 0-based transverse indices relative to this rank's block INTERSECTION (o1/o2 = local
            ! transverse origins), aligned with the fine registers: fine children of isect-relative cell t are faces 2*t and 2*t+1.
            ! At np=1 the intersection is the block and both flags hold, recovering single-rank behavior exactly.
            ! ONE coarse s_compute_rhs pass fills EVERY active block's registers: revisit each slot's region+intersection in turn.
            save_cur = amr_cur
            bactive = .false.
            maxt1 = 0; maxt2 = 0
            do islot = 1, amr_num_blocks
                ! a level>=2 block's coarse side is its PARENT (creg captured in the fine branch), not L0
                if (amr_block_level(islot) >= 2) cycle
                call s_amr_select_slot(islot)
                call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
                cap_lo = own_lo(id); cap_hi = own_hi(id)
                if (cap_lo .or. cap_hi) then
                    ! block-relative transverse frame (0-based from region_lo, aligned with the owner's freg): this rank fills creg
                    ! over its owned overlap [tlo-region_lo : thi-region_lo]; o1/o2 map that back to LOCAL flux indices.
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
                    ! cap_lo/cap_hi is s_amr_reg_prepare's clause (c) verbatim, so amr_reg_of(islot) is always mapped here
                    sreg = amr_reg_of(islot)
                    bactive(sreg) = .true.; bclo(sreg) = cap_lo; bchi(sreg) = cap_hi
                    bjlo(sreg) = jlo; bjhi(sreg) = jhi; bo1(sreg) = o1; bo2(sreg) = o2
                    bt1lo(sreg) = t1_lo; bt1hi(sreg) = t1_hi; bt2lo(sreg) = t2_lo; bt2hi(sreg) = t2_hi
                    maxt1 = max(maxt1, t1_hi); maxt2 = max(maxt2, t2_hi)
                end if  ! cap_lo .or. cap_hi
            end do
            call s_amr_select_slot(save_cur)
            if (any(bactive(1:amr_reg_n))) then
                $:GPU_UPDATE(device='[bjlo, bjhi, bo1, bo2, bt1lo, bt1hi, bt2lo, bt2hi, bclo, bchi, bactive]')
                ! shared capture into each coarse block's creg (region/sidx frame, per-face ownership gating): advective, then
                ! total-flux viscous, then chemistry species+energy
                call s_amr_capture_creg_dense_batch(amr_reg_n, id, .true., coef, accum, maxt1, maxt2, 1, sys_size)
                if (viscous) call s_amr_capture_creg_dense_batch(amr_reg_n, id, .false., coef, .true., maxt1, maxt2, &
                    & eqn_idx%mom%beg, eqn_idx%E)
                if (chemistry .and. chem_params%diffusion) call s_amr_capture_creg_chem_batch(amr_reg_n, id, coef, maxt1, maxt2)
            end if
        end if

    end subroutine s_amr_capture_boundary_flux

    !> Correct the coarse rhs in the first cell OUTSIDE each block face so the coarse update sees the (child-averaged) fine flux at
    !! every c/f face. Signs follow rhs = (flux_left - flux_right)/dx: low face is the outside cell's RIGHT face => rhs += (F_coarse
    !! - Fbar_fine)/dx; high face is the outside cell's LEFT face => rhs += (Fbar_fine - F_coarse)/dx. Cells INSIDE the block need
    !! no correction (end-of-step restriction overwrites them). c1/c2 are relative 0-based coarse transverse indices.
    impure subroutine s_amr_apply_reflux(rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: eq, c1, c2, c1w, c2w, k, save_cur, nact, gmax1, gmax2
        integer                                                :: f10, f20, dd1, dd2, nch, rr, dd1_hi, dd2_hi, sreg
        integer                                                :: bl1, bh1, bl2, bh2, bl3, bh3
        integer                                                :: i2, i3, sidx(3), ext(3), tlo(3), thi(3)
        logical                                                :: d2, d3, own_lo(3), own_hi(3)
        real(wp)                                               :: fblo, fbhi, wsum, rf

        if (.not. amr) return
        ! Refresh the participation map + register capacity on a topology change; no-op (two integer compares) otherwise.
        call s_amr_reg_prepare()
        if (igr) return  ! stage-1 IGR: restriction-only coupling (no captured fluxes)
        rr = amr_ref_ratio
        d2 = n_glb > 0; d3 = p_glb > 0
        save_cur = amr_cur

        ! BATCHED over the level-1 blocks, one kernel per face direction (mirror of the capture-side batching,
        ! s_amr_capture_creg_dense_batch): the per-box form launched up to 3 tiny face kernels per block per step, and the
        ! per-launch overhead - not the arithmetic - dominated the reflux-apply phase. The per-(face, eq, cell) arithmetic
        ! and child-sum order below are IDENTICAL to the per-box form, and block corrections are disjoint (the merge
        ! invariant keeps blocks >= buff_size apart), so the outputs are bit-identical. Host precompute walks the slots
        ! with the SAME select_slot + s_amr_reflux_face_flags the per-box form used; the a_* descriptors are pushed once
        ! per direction.

        ! x-faces: transverse dims (y, z); children in each active transverse dim
        nch = 1
        if (n_glb > 0) nch = nch*rr
        if (p_glb > 0) nch = nch*rr
        dd1_hi = merge(rr - 1, 0, n_glb > 0); dd2_hi = merge(rr - 1, 0, p_glb > 0)
        nact = 0; gmax1 = 0; gmax2 = 0
        a_act = .false.
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
            if (.not. (own_lo(1) .or. own_hi(1))) cycle
            bl2 = tlo(2) - amr_region_lo(2); bh2 = thi(2) - amr_region_lo(2)
            bl3 = tlo(3) - amr_region_lo(3); bh3 = thi(3) - amr_region_lo(3)
            ! own_lo/own_hi is s_amr_reg_prepare's clause (c) verbatim, so amr_reg_of(k) is always mapped here
            sreg = amr_reg_of(k)
            a_act(sreg) = .true.; a_lo(sreg) = own_lo(1); a_hi(sreg) = own_hi(1)
            a_ol(sreg) = amr_region_lo(1) - 1 - sidx(1); a_oh(sreg) = amr_region_hi(1) + 1 - sidx(1)
            a_t2(sreg) = amr_region_lo(2) - sidx(2); a_t3(sreg) = amr_region_lo(3) - sidx(3)
            a_b1l(sreg) = bl2; a_b1h(sreg) = bh2; a_b2l(sreg) = bl3; a_b2h(sreg) = bh3
            a_mlo(sreg) = 1._wp; a_mhi(sreg) = 1._wp
            if (own_lo(1)) a_mlo(sreg) = dx(a_ol(sreg))
            if (own_hi(1)) a_mhi(sreg) = dx(a_oh(sreg))
            nact = nact + 1
            gmax1 = max(gmax1, bh2 - bl2); gmax2 = max(gmax2, bh3 - bl3)
        end do
        call s_amr_select_slot(save_cur)
        if (nact > 0) then
            $:GPU_UPDATE(device='[a_ol, a_oh, a_t2, a_t3, a_b1l, a_b1h, a_b2l, a_b2h, a_lo, a_hi, a_act, a_mlo, a_mhi]')
            if (cyl_coord) then
                ! axisymmetric x-face (axial): the rr covering fine faces stack in the RADIAL (transverse) direction at
                ! DIFFERENT radii, so Fbar_fine must be area-weighted by fine-face radius (fine y_cc rebuilt from the coarse
                ! y_cb of transverse cell tl2+c1). Outside-cell axial divergence has no radial factor (axial face area ~
                ! cell volume ~ y_cc, cancels), so the width stays dx.
                $:GPU_PARALLEL_LOOP(collapse=4, private='[c1, c2, f10, f20, dd1, fblo, fbhi, wsum, rf, i2, i3]')
                do k = 1, amr_reg_n
                    do c2w = 0, gmax2
                        do c1w = 0, gmax1
                            do eq = 1, sys_size
                                if (.not. a_act(k)) cycle
                                c1 = a_b1l(k) + c1w; c2 = a_b2l(k) + c2w
                                if (c1 > a_b1h(k) .or. c2 > a_b2h(k)) cycle
                                f20 = 0
                                f10 = rr*c1
                                fblo = 0._wp; fbhi = 0._wp; wsum = 0._wp
                                do dd1 = 0, dd1_hi
                                    rf = y_cb(a_t2(k) + c1 - 1) + (real(dd1, &
                                              & wp) + 0.5_wp)*(y_cb(a_t2(k) + c1) - y_cb(a_t2(k) + c1 - 1))/real(rr, wp)
                                    fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20, k)*rf
                                    fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20, k)*rf
                                    wsum = wsum + rf
                                end do
                                fblo = fblo/wsum; fbhi = fbhi/wsum
                                i2 = a_t2(k) + c1; i3 = a_t3(k) + c2
                                if (a_lo(k)) rhs_vf(eq)%sf(a_ol(k), i2, i3) = rhs_vf(eq)%sf(a_ol(k), i2, i3) + (creg(1)%lo(eq, &
                                    & c1, c2, k) - fblo)/a_mlo(k)
                                if (a_hi(k)) rhs_vf(eq)%sf(a_oh(k), i2, i3) = rhs_vf(eq)%sf(a_oh(k), i2, &
                                    & i3) + (fbhi - creg(1)%hi(eq, c1, c2, k))/a_mhi(k)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else
                $:GPU_PARALLEL_LOOP(collapse=4, private='[c1, c2, f10, f20, dd1, dd2, fblo, fbhi, i2, i3]')
                do k = 1, amr_reg_n
                    do c2w = 0, gmax2
                        do c1w = 0, gmax1
                            do eq = 1, sys_size
                                if (.not. a_act(k)) cycle
                                c1 = a_b1l(k) + c1w; c2 = a_b2l(k) + c2w
                                if (c1 > a_b1h(k) .or. c2 > a_b2h(k)) cycle
                                f20 = 0; if (d3) f20 = rr*c2
                                f10 = 0; if (d2) f10 = rr*c1
                                fblo = 0._wp; fbhi = 0._wp
                                do dd2 = 0, dd2_hi
                                    do dd1 = 0, dd1_hi
                                        fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20 + dd2, k)
                                        fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20 + dd2, k)
                                    end do
                                end do
                                fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                                i2 = a_t2(k) + c1; i3 = a_t3(k) + c2
                                if (a_lo(k)) rhs_vf(eq)%sf(a_ol(k), i2, i3) = rhs_vf(eq)%sf(a_ol(k), i2, i3) + (creg(1)%lo(eq, &
                                    & c1, c2, k) - fblo)/a_mlo(k)
                                if (a_hi(k)) rhs_vf(eq)%sf(a_oh(k), i2, i3) = rhs_vf(eq)%sf(a_oh(k), i2, &
                                    & i3) + (fbhi - creg(1)%hi(eq, c1, c2, k))/a_mhi(k)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

        ! y-faces (n_glb > 0): transverse dims (x, z); x is always active (2 children)
        if (n_glb > 0) then
            nch = rr
            if (p_glb > 0) nch = nch*rr
            dd2_hi = merge(rr - 1, 0, p_glb > 0)
            nact = 0; gmax1 = 0; gmax2 = 0
            a_act = .false.
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= 1) cycle
                call s_amr_select_slot(k)
                call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
                if (.not. (own_lo(2) .or. own_hi(2))) cycle
                bl1 = tlo(1) - amr_region_lo(1); bh1 = thi(1) - amr_region_lo(1)
                bl3 = tlo(3) - amr_region_lo(3); bh3 = thi(3) - amr_region_lo(3)
                sreg = amr_reg_of(k)
                a_act(sreg) = .true.; a_lo(sreg) = own_lo(2); a_hi(sreg) = own_hi(2)
                a_ol(sreg) = amr_region_lo(2) - 1 - sidx(2); a_oh(sreg) = amr_region_hi(2) + 1 - sidx(2)
                a_t1(sreg) = amr_region_lo(1) - sidx(1); a_t3(sreg) = amr_region_lo(3) - sidx(3)
                a_b1l(sreg) = bl1; a_b1h(sreg) = bh1; a_b2l(sreg) = bl3; a_b2h(sreg) = bh3
                a_mlo(sreg) = 1._wp; a_mhi(sreg) = 1._wp
                if (own_lo(2)) a_mlo(sreg) = dy(a_ol(sreg))
                if (own_hi(2)) a_mhi(sreg) = dy(a_oh(sreg))
                ! cyl_coord (axisymmetric): the radial c/f flux correction is area-weighted - low/high face carries radius
                ! y_cb, outside cell volume carries y_cc, so fold r_face/r_cell into the width (kernel divides by it).
                if (cyl_coord) then
                    if (own_lo(2)) a_mlo(sreg) = a_mlo(sreg)*y_cc(a_ol(sreg))/y_cb(a_ol(sreg))
                    if (own_hi(2)) a_mhi(sreg) = a_mhi(sreg)*y_cc(a_oh(sreg))/y_cb(a_oh(sreg) - 1)
                end if
                nact = nact + 1
                gmax1 = max(gmax1, bh1 - bl1); gmax2 = max(gmax2, bh3 - bl3)
            end do
            call s_amr_select_slot(save_cur)
            if (nact > 0) then
                $:GPU_UPDATE(device='[a_ol, a_oh, a_t1, a_t3, a_b1l, a_b1h, a_b2l, a_b2h, a_lo, a_hi, a_act, a_mlo, a_mhi]')
                $:GPU_PARALLEL_LOOP(collapse=4, private='[c1, c2, f10, f20, dd1, dd2, fblo, fbhi, i2, i3]')
                do k = 1, amr_reg_n
                    do c2w = 0, gmax2
                        do c1w = 0, gmax1
                            do eq = 1, sys_size
                                if (.not. a_act(k)) cycle
                                c1 = a_b1l(k) + c1w; c2 = a_b2l(k) + c2w
                                if (c1 > a_b1h(k) .or. c2 > a_b2h(k)) cycle
                                f20 = 0; if (d3) f20 = rr*c2
                                f10 = rr*c1
                                fblo = 0._wp; fbhi = 0._wp
                                do dd2 = 0, dd2_hi
                                    do dd1 = 0, rr - 1
                                        fblo = fblo + freg(2)%lo(eq, f10 + dd1, f20 + dd2, k)
                                        fbhi = fbhi + freg(2)%hi(eq, f10 + dd1, f20 + dd2, k)
                                    end do
                                end do
                                fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                                i2 = a_t1(k) + c1; i3 = a_t3(k) + c2
                                if (a_lo(k)) rhs_vf(eq)%sf(i2, a_ol(k), i3) = rhs_vf(eq)%sf(i2, a_ol(k), i3) + (creg(2)%lo(eq, &
                                    & c1, c2, k) - fblo)/a_mlo(k)
                                if (a_hi(k)) rhs_vf(eq)%sf(i2, a_oh(k), i3) = rhs_vf(eq)%sf(i2, a_oh(k), &
                                    & i3) + (fbhi - creg(2)%hi(eq, c1, c2, k))/a_mhi(k)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

        ! z-faces (p_glb > 0): transverse dims (x, y); both always active in 3D (4 children)
        if (p_glb > 0) then
            nch = rr*rr
            nact = 0; gmax1 = 0; gmax2 = 0
            a_act = .false.
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= 1) cycle
                call s_amr_select_slot(k)
                call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
                if (.not. (own_lo(3) .or. own_hi(3))) cycle
                bl1 = tlo(1) - amr_region_lo(1); bh1 = thi(1) - amr_region_lo(1)
                bl2 = tlo(2) - amr_region_lo(2); bh2 = thi(2) - amr_region_lo(2)
                sreg = amr_reg_of(k)
                a_act(sreg) = .true.; a_lo(sreg) = own_lo(3); a_hi(sreg) = own_hi(3)
                a_ol(sreg) = amr_region_lo(3) - 1 - sidx(3); a_oh(sreg) = amr_region_hi(3) + 1 - sidx(3)
                a_t1(sreg) = amr_region_lo(1) - sidx(1); a_t2(sreg) = amr_region_lo(2) - sidx(2)
                a_b1l(sreg) = bl1; a_b1h(sreg) = bh1; a_b2l(sreg) = bl2; a_b2h(sreg) = bh2
                a_mlo(sreg) = 1._wp; a_mhi(sreg) = 1._wp
                if (own_lo(3)) a_mlo(sreg) = dz(a_ol(sreg))
                if (own_hi(3)) a_mhi(sreg) = dz(a_oh(sreg))
                nact = nact + 1
                gmax1 = max(gmax1, bh1 - bl1); gmax2 = max(gmax2, bh2 - bl2)
            end do
            call s_amr_select_slot(save_cur)
            if (nact > 0) then
                $:GPU_UPDATE(device='[a_ol, a_oh, a_t1, a_t2, a_b1l, a_b1h, a_b2l, a_b2h, a_lo, a_hi, a_act, a_mlo, a_mhi]')
                $:GPU_PARALLEL_LOOP(collapse=4, private='[c1, c2, f10, f20, dd1, dd2, fblo, fbhi, i2, i3]')
                do k = 1, amr_reg_n
                    do c2w = 0, gmax2
                        do c1w = 0, gmax1
                            do eq = 1, sys_size
                                if (.not. a_act(k)) cycle
                                c1 = a_b1l(k) + c1w; c2 = a_b2l(k) + c2w
                                if (c1 > a_b1h(k) .or. c2 > a_b2h(k)) cycle
                                f20 = rr*c2
                                f10 = rr*c1
                                fblo = 0._wp; fbhi = 0._wp
                                do dd2 = 0, rr - 1
                                    do dd1 = 0, rr - 1
                                        fblo = fblo + freg(3)%lo(eq, f10 + dd1, f20 + dd2, k)
                                        fbhi = fbhi + freg(3)%hi(eq, f10 + dd1, f20 + dd2, k)
                                    end do
                                end do
                                fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                                i2 = a_t1(k) + c1; i3 = a_t2(k) + c2
                                if (a_lo(k)) rhs_vf(eq)%sf(i2, i3, a_ol(k)) = rhs_vf(eq)%sf(i2, i3, a_ol(k)) + (creg(3)%lo(eq, &
                                    & c1, c2, k) - fblo)/a_mlo(k)
                                if (a_hi(k)) rhs_vf(eq)%sf(i2, i3, a_oh(k)) = rhs_vf(eq)%sf(i2, i3, &
                                    & a_oh(k)) + (fbhi - creg(3)%hi(eq, c1, c2, k))/a_mhi(k)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

    end subroutine s_amr_apply_reflux

    !> Zero the fine registers (called by the subcycle driver before substep 1 - stage-1 overwrite cannot work across two substeps).
    impure subroutine s_amr_zero_fine_registers()

        integer :: d, eq, t1, t2, t1_hi, t2_hi, islot

        if (.not. amr) return
        ! Refresh the participation map + register capacity on a topology change; no-op (two integer compares) otherwise.
        call s_amr_reg_prepare()
        if (igr) return  ! stage-1 IGR: restriction-only coupling (no captured fluxes)
        if (.not. amr_rank_owns_block) return
        islot = amr_reg_cur  ! working block's DENSE register slot (local => captured by value in the device kernels below)
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
        integer :: d, sidx(3), ext(3), tlo(3), thi(3), olo(3), ohi(3), glo(3), ghi(3), woff(3)
        logical :: own_lo(3), own_hi(3)
        real(wp) :: w_lo(3), w_hi(3), mlo(3), mhi(3)

        if (.not. amr) return
        ! Refresh the participation map + register capacity on a topology change; no-op (two integer compares) otherwise.
        call s_amr_reg_prepare()
        if (igr) return  ! stage-1 IGR: restriction-only coupling (no captured fluxes)
        call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
        if (.not. (any(own_lo) .or. any(own_hi))) return
        ! L0/L1 (coarse) frame for the shared kernel: outside cell = region boundary +/-1 in local (sidx-offset) coords; creg-local
        ! loop range is the owned transverse overlap [tlo:thi] block-relative; ownership -> unit face weights; cell widths from the
        ! global coarse grid (amr_ref_ratio = 2, dt = coarse step).
        olo = 0; ohi = 0; glo = 0; ghi = 0; woff = 0; w_lo = 0._wp; w_hi = 0._wp; mlo = 1._wp; mhi = 1._wp
        do d = 1, num_dims
            olo(d) = amr_region_lo(d) - 1 - sidx(d); ohi(d) = amr_region_hi(d) + 1 - sidx(d)
            glo(d) = tlo(d) - amr_region_lo(d); ghi(d) = thi(d) - amr_region_lo(d)
            woff(d) = amr_region_lo(d) - sidx(d)
            if (own_lo(d)) w_lo(d) = 1._wp
            if (own_hi(d)) w_hi(d) = 1._wp
        end do
        if (own_lo(1)) mlo(1) = dx(olo(1))
        if (own_hi(1)) mhi(1) = dx(ohi(1))
        if (n_glb > 0) then
            if (own_lo(2)) mlo(2) = dy(olo(2))
            if (own_hi(2)) mhi(2) = dy(ohi(2))
            ! cyl_coord (axisymmetric): area-weight the radial c/f correction by r_face/r_cell (mirror of s_amr_apply_reflux). L0/L1
            ! coarse frame -> global y_cb/y_cc for the owned outside cell. r_+ = y_cb(olo(2)) low; r_- = y_cb(ohi(2)-1) high.
            if (cyl_coord) then
                if (own_lo(2)) mlo(2) = mlo(2)*y_cc(olo(2))/y_cb(olo(2))
                if (own_hi(2)) mhi(2) = mhi(2)*y_cc(ohi(2))/y_cb(ohi(2) - 1)
            end if
        end if
        if (p_glb > 0) then
            if (own_lo(3)) mlo(3) = dz(olo(3))
            if (own_hi(3)) mhi(3) = dz(ohi(3))
        end if
        call s_amr_reflux_apply_faces(q_cons, amr_reg_cur, 2, dt, olo, ohi, glo, ghi, woff, w_lo, w_hi, mlo, mhi)

    end subroutine s_amr_apply_reflux_state

    !> Shared Berger-Colella STATE reflux kernel: apply q(outside) += w*dtl*(F_coarse - Fbar_fine)/m on the low face and +=
    !! w*dtl*(Fbar_fine - F_coarse)/m on the high face for each active dim, where F_coarse is creg and Fbar_fine averages freg over
    !! the rr**(ndim-1) covering fine faces. Used by BOTH s_amr_apply_reflux_state (L0/L1, coarse/sidx frame, unit weights from
    !! ownership, rr=2) and s_amr_reflux_to_parent (L2->L1, parent-fine frame, sibling-seam weights, rr=amr_ref_ratio). All framing
    !! is caller-passed so the flux-correction math is single-sourced: islot - DENSE register slot (amr_reg_cur); rr - refinement
    !! ratio (fine faces per coarse face per transverse dim); dtl - reflux dt; olo/ohi(d) - outside coarse-cell index just
    !! below/above the block face in dim d; glo/ghi(d) - creg-local loop range in dim d (transverse for the two faces d' /= d);
    !! woff(d) - transverse write origin so the cell index is woff(d) + g; w_lo/w_hi(d) - per-face weight (0 skips the write:
    !! unowned face at np>1, or a fine-fine sibling-tile seam); mlo/mhi(d) - outside-cell width for the low/high face
    !! (invalid/unused where weight is 0). A zero weight SKIPS the write (not multiply-by-0) because the outside index may be out of
    !! bounds on an unowned face.
    impure subroutine s_amr_reflux_apply_faces(q, islot, rr, dtl, olo, ohi, glo, ghi, woff, w_lo, w_hi, mlo, mhi)

        type(scalar_field), dimension(sys_size), intent(inout) :: q
        integer, intent(in) :: islot, rr, olo(3), ohi(3), glo(3), ghi(3), woff(3)
        real(wp), intent(in) :: dtl, w_lo(3), w_hi(3), mlo(3), mhi(3)
        integer :: eq, g1, g2, f10, f20, dd1, dd2, nch, dd1_hi, dd2_hi, ol, oh, w2, w3, w1, gl1, gh1, gl2, gh2, gl3, gh3
        real(wp) :: fblo, fbhi, wl, wh, ml, mh, wsum, rf

        ! loop bounds hoisted to scalars: array-element bounds (glo(d)/ghi(d)) drive the collapsed inner loop and would force the
        ! host arrays present on the device (an ACC present error)

        gl1 = glo(1); gh1 = ghi(1); gl2 = glo(2); gh2 = ghi(2); gl3 = glo(3); gh3 = ghi(3)

        ! x-faces: transverse (y, z)
        if (w_lo(1) /= 0._wp .or. w_hi(1) /= 0._wp) then
            nch = 1; if (n_glb > 0) nch = nch*rr; if (p_glb > 0) nch = nch*rr
            dd1_hi = merge(rr - 1, 0, n_glb > 0); dd2_hi = merge(rr - 1, 0, p_glb > 0)
            ol = olo(1); oh = ohi(1); w2 = woff(2); w3 = woff(3); wl = w_lo(1); wh = w_hi(1); ml = mlo(1); mh = mhi(1)
            if (cyl_coord) then
                ! axisymmetric x-face: area-weight Fbar_fine by fine-face radius (rebuilt from the coarse y_cb of transverse cell
                ! w2+g1) - the rr covering fine faces sit at different radii. cyl reaches here only single-level (L0 frame), so
                ! global y_cb is the correct coarse grid.
                $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi, wsum, rf]')
                do eq = 1, sys_size
                    do g2 = gl3, gh3
                        do g1 = gl2, gh2
                            f20 = 0
                            f10 = rr*g1
                            fblo = 0._wp; fbhi = 0._wp; wsum = 0._wp
                            do dd1 = 0, dd1_hi
                                rf = y_cb(w2 + g1 - 1) + (real(dd1, wp) + 0.5_wp)*(y_cb(w2 + g1) - y_cb(w2 + g1 - 1))/real(rr, wp)
                                fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20, islot)*rf
                                fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20, islot)*rf
                                wsum = wsum + rf
                            end do
                            fblo = fblo/wsum; fbhi = fbhi/wsum
                            if (wl /= 0._wp) q(eq)%sf(ol, w2 + g1, w3 + g2) = q(eq)%sf(ol, w2 + g1, &
                                & w3 + g2) + wl*dtl*(creg(1)%lo(eq, g1, g2, islot) - fblo)/ml
                            if (wh /= 0._wp) q(eq)%sf(oh, w2 + g1, w3 + g2) = q(eq)%sf(oh, w2 + g1, &
                                & w3 + g2) + wh*dtl*(fbhi - creg(1)%hi(eq, g1, g2, islot))/mh
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else
                $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
                do eq = 1, sys_size
                    do g2 = gl3, gh3
                        do g1 = gl2, gh2
                            f20 = 0; if (p_glb > 0) f20 = rr*g2
                            f10 = 0; if (n_glb > 0) f10 = rr*g1
                            fblo = 0._wp; fbhi = 0._wp
                            do dd2 = 0, dd2_hi
                                do dd1 = 0, dd1_hi
                                    fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                    fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20 + dd2, islot)
                                end do
                            end do
                            fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                            if (wl /= 0._wp) q(eq)%sf(ol, w2 + g1, w3 + g2) = q(eq)%sf(ol, w2 + g1, &
                                & w3 + g2) + wl*dtl*(creg(1)%lo(eq, g1, g2, islot) - fblo)/ml
                            if (wh /= 0._wp) q(eq)%sf(oh, w2 + g1, w3 + g2) = q(eq)%sf(oh, w2 + g1, &
                                & w3 + g2) + wh*dtl*(fbhi - creg(1)%hi(eq, g1, g2, islot))/mh
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if
        ! y-faces (n_glb > 0): transverse (x, z); x always active
        if (n_glb > 0 .and. (w_lo(2) /= 0._wp .or. w_hi(2) /= 0._wp)) then
            nch = rr; if (p_glb > 0) nch = nch*rr
            dd2_hi = merge(rr - 1, 0, p_glb > 0)
            ol = olo(2); oh = ohi(2); w1 = woff(1); w3 = woff(3); wl = w_lo(2); wh = w_hi(2); ml = mlo(2); mh = mhi(2)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do g2 = gl3, gh3
                    do g1 = gl1, gh1
                        f20 = 0; if (p_glb > 0) f20 = rr*g2
                        f10 = rr*g1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, dd2_hi
                            do dd1 = 0, rr - 1
                                fblo = fblo + freg(2)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                fbhi = fbhi + freg(2)%hi(eq, f10 + dd1, f20 + dd2, islot)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (wl /= 0._wp) q(eq)%sf(w1 + g1, ol, w3 + g2) = q(eq)%sf(w1 + g1, ol, w3 + g2) + wl*dtl*(creg(2)%lo(eq, &
                            & g1, g2, islot) - fblo)/ml
                        if (wh /= 0._wp) q(eq)%sf(w1 + g1, oh, w3 + g2) = q(eq)%sf(w1 + g1, oh, &
                            & w3 + g2) + wh*dtl*(fbhi - creg(2)%hi(eq, g1, g2, islot))/mh
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        ! z-faces (p_glb > 0): transverse (x, y); both active in 3D
        if (p_glb > 0 .and. (w_lo(3) /= 0._wp .or. w_hi(3) /= 0._wp)) then
            nch = rr*rr
            ol = olo(3); oh = ohi(3); w1 = woff(1); w2 = woff(2); wl = w_lo(3); wh = w_hi(3); ml = mlo(3); mh = mhi(3)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do g2 = gl2, gh2
                    do g1 = gl1, gh1
                        f20 = rr*g2
                        f10 = rr*g1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, rr - 1
                            do dd1 = 0, rr - 1
                                fblo = fblo + freg(3)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                fbhi = fbhi + freg(3)%hi(eq, f10 + dd1, f20 + dd2, islot)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (wl /= 0._wp) q(eq)%sf(w1 + g1, w2 + g2, ol) = q(eq)%sf(w1 + g1, w2 + g2, ol) + wl*dtl*(creg(3)%lo(eq, &
                            & g1, g2, islot) - fblo)/ml
                        if (wh /= 0._wp) q(eq)%sf(w1 + g1, w2 + g2, oh) = q(eq)%sf(w1 + g1, w2 + g2, &
                            & oh) + wh*dtl*(fbhi - creg(3)%hi(eq, g1, g2, islot))/mh
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_amr_reflux_apply_faces

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
        @:DEALLOCATE(bjlo, bjhi, bo1, bo2, bt1lo, bt1hi, bt2lo, bt2hi, bclo, bchi, bactive)
        @:DEALLOCATE(a_ol, a_oh, a_t1, a_t2, a_t3, a_b1l, a_b1h, a_b2l, a_b2h, a_lo, a_hi, a_act, a_mlo, a_mhi)
        if (allocated(amr_reg_of)) deallocate (amr_reg_of)
        amr_reg_n = 0; amr_reg_cur = 0

    end subroutine s_finalize_amr_registers

end module m_amr_registers
