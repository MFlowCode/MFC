!>
!!@file
!!@brief Contains module m_amr_registers

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR).
#! The module arrays these kernels name are the batch tables a_*/b* (allocated unconditionally at
#! reserve), flux_rsx_vf/flux_src_rsx_vf (m_riemann_solvers; allocated whenever these kernels can
#! launch, both under .not. igr) and the local rtmp_d, which @:ALLOCATE puts on the device.
#! A kernel naming an unallocated array aborts, so any new kernel here must only name arrays
#! allocated on every path that can reach it.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief AMR flux registers: per-RK-stage refluxing at the coarse/fine block boundary. Depends only on m_derived_types +
!! m_global_parameters so both m_rhs (capture) and m_time_steppers (apply) can use it without cycles: region info is read from
!! amr_region_lo/hi and amr_isect_lo/hi (m_global_parameters), mirrored across regrids by s_set_amr_fine_geometry. creg uses 0-based
!! transverse indexing relative to the rank's block intersection; freg uses 0-based local fine (fine children of isect cell t are
!! 2*t and 2*t+1). The slot dimension grows on demand (s_amr_reg_reserve); the transverse extents are fixed at init.
!!
!! What is captured: the uniform 1:sys_size capture of flux_rsx_vf refluxes per-fluid masses, momentum, energy and the
!! volume-fraction advective part alpha_i*u_star (5-eq HLLC). The non-conservative +alpha*d(u_star)/dx remainder is deliberately not
!! captured (alpha is genuinely non-conservative); coarse/fine volume-fraction consistency is held by mpp_lim instead. The viscous
!! stress/work fluxes (mom..E) and the chemistry species-diffusion and conduction fluxes (species; energy only when not viscous)
!! travel through flux_src_n with the same face indexing and sign, and are added into the same registers so the reflux matches the
!! total flux. Fine-ghost gradients at the c/f boundary come from the cons prolongation; that inconsistency is bounded, and
!! conservation is enforced by the flux-register matching.
module m_amr_registers

    use m_derived_types
    use m_global_parameters
    use m_riemann_state, only: flux_rsx_vf, flux_src_rsx_vf

    implicit none

    private; public :: s_initialize_amr_registers, s_amr_capture_boundary_flux, s_amr_apply_reflux, s_amr_zero_fine_registers, &
        & s_finalize_amr_registers, s_amr_reflux_face_flags, s_amr_reflux_faces, s_amr_reflux_apply_faces, s_amr_parent_foot, &
        & s_amr_reg_prepare, freg, creg, f_amr_face_is_seam

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
    !! The registers are not dimensioned 1:amr_max_blocks: amr_max_blocks is a safety cap, not a block count, and each unused
    !! slot costs device memory. They grow geometrically like the flat store (s_amr_st_reserve), so the cap bounds correctness
    !! only and memory follows the actual block count.
    integer, parameter :: amr_reg_floor = 64  !< initial slot capacity; growth doubles from here
    integer            :: amr_reg_cap = 0
    integer            :: rc(3), rf(3)        !< coarse / fine transverse extents (creg is 0:rc-1, freg is 0:rf)
    !> Mesh-epoch/tripwire keys of the last participation-map build (s_amr_reg_prepare; mirror of the seam-pair cache keys).
    integer(8) :: amr_reg_epoch_built = -1_8
    integer    :: amr_reg_nblk_built = -1

    !> Per-slot geometry scratch for the batched capture kernels (1:amr_max_blocks): host-filled by s_amr_capture_slot, then
    !! GPU_UPDATE'd so one kernel iterates the slot dimension instead of O(blocks) tiny launches. bactive gates the slot;
    !! bt1lo/bt1hi/bt2lo/bt2hi are the per-slot transverse window (a slot outside the rectangular caps bmax1/bmax2 is cycled);
    !! bjlo/bjhi are the normal-face flux indices; bo1/bo2 the transverse origins; bclo/bchi the per-face capture gates.
    integer, allocatable :: bjlo(:), bjhi(:), bo1(:), bo2(:), bt1lo(:), bt1hi(:), bt2lo(:), bt2hi(:)
    logical, allocatable :: bclo(:), bchi(:), bactive(:)
    integer              :: bmax1 = 0, bmax2 = 0
    $:GPU_DECLARE(create='[bjlo, bjhi, bo1, bo2, bt1lo, bt1hi, bt2lo, bt2hi, bclo, bchi, bactive]')

    !> Per-slot geometry scratch for the batched reflux apply kernels (mirror of the capture batching above): a_act gates the slot,
    !! a_lo/a_hi the per-face applies, a_ol/a_oh the outside coarse cell's local index in the face dim, a_ta/a_tb the local
    !! transverse origins, a_b1l..a_b2h the transverse windows (block-relative, freg/creg-aligned), a_mlo/a_mhi the outside-cell
    !! widths. Filled per direction on the host, GPU_UPDATE'd, consumed by one kernel per face direction instead of O(blocks) tiny
    !! launches.
    integer, allocatable  :: a_ol(:), a_oh(:), a_ta(:), a_tb(:), a_b1l(:), a_b1h(:), a_b2l(:), a_b2h(:)
    logical, allocatable  :: a_lo(:), a_hi(:), a_act(:)
    real(wp), allocatable :: a_mlo(:), a_mhi(:)
    $:GPU_DECLARE(create='[a_ol, a_oh, a_ta, a_tb, a_b1l, a_b1h, a_b2l, a_b2h, a_lo, a_hi, a_act, a_mlo, a_mhi]')

contains

    #:def REG_GROW(A, L2, U2, L3, U3)
        if (allocated(${A}$)) then
            if (oldcap > amr_reg_grow_dev_cap) then
                ! near-limit fallback (mirror of s_amr_st_reserve): the device staging below transiently holds
                ! old + tmp = 2*oldcap slots on the device and growth fires at the memory high-water mark, so above
                ! the threshold keep the host round trip: slow, but its device peak is max(old, new).
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
                ! stage the live slots on the device (rtmp_d is device-mapped by @:ALLOCATE); no PCIe traffic
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
    !! Keyed on amr_reg_n, the dense participation-local count (s_amr_reg_prepare): register slots are dense indices from
    !! amr_reg_of, so capacity follows what this rank owns or participates in, not the global block count. The capture and
    !! apply kernels sweep slot = 1..amr_reg_n with bactive/a_act gating the filled subset.
    !!
    !! Contents are preserved across growth, mirroring s_amr_st_reserve: every caller of s_amr_alloc_slot is a between-step
    !! operation (regrid, restart, slot reconcile, L0 tile build) and stage 1 overwrites the registers, but freg accumulates
    !! across RK stages, and a silent conservation error is not worth the saving. Growth stages
    !! through a device temporary like the store (the registers are device-authoritative; every host consumer pulls its slot to
    !! the host immediately before reading, so the host mirror coming out of a device-path growth is undefined, the same contract
    !! as the store). Above the transient threshold the host round trip is the OOM-safe path.
    impure subroutine s_amr_reg_reserve(nslot)

        integer, intent(in) :: nslot
        integer             :: oldcap, newcap
        integer             :: eq, t1, t2, c4
        !> device staging transiently doubles one register array's footprint; register slots are faces (small compared with a store
        !! column), so a slot-count guard is enough here. Above it, fall back to the host round trip (device peak max(old, new)).
        integer, parameter    :: amr_reg_grow_dev_cap = 512
        real(wp), allocatable :: rtmp(:,:,:,:), rtmp_d(:,:,:,:)

        if (nslot <= amr_reg_cap) return
        oldcap = amr_reg_cap
        newcap = min(amr_max_blocks, max(2*oldcap, nslot))

        #:for D, A, B in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            @:REG_GROW(creg(${D}$)%lo, 0, rc(${A}$) - 1, 0, rc(${B}$) - 1)
            @:REG_GROW(creg(${D}$)%hi, 0, rc(${A}$) - 1, 0, rc(${B}$) - 1)
            @:REG_GROW(freg(${D}$)%lo, 0, rf(${A}$), 0, rf(${B}$))
            @:REG_GROW(freg(${D}$)%hi, 0, rf(${A}$), 0, rf(${B}$))
        #:endfor

        amr_reg_cap = newcap

    end subroutine s_amr_reg_reserve

    !> Build/refresh the participation-local register index (amr_reg_of/amr_reg_n, m_global_parameters) and size the registers to
    !! the dense count. Lazily keyed on the mesh epoch (every regrid, migration, restart, and slot renumbering bumps it; a
    !! block-count change is the tripwire). A global slot g maps iff this rank (a) owns g, (b) owns g's parent (the parent-side
    !! child-creg capture, the freg receives, and the reflux-to-parent apply all index the child's slot on the parent's owner), or
    !! (c) reflux-face-participates in g per s_amr_reflux_face_flags (the coarse capture, the L0/L1 apply, and the reflux face-wave
    !! receives are gated by exactly these flags, so the map cannot under-cover them). This keeps the register footprint O(owned +
    !! participation halo) rather than O(global boxes). The mapped range is zeroed after a rebuild: dense slots alias across
    !! rebuilds (block g's new slot may hold another block's stale flux), and zeroing keeps the garbage-until-captured contract
    !! deterministic. Contents at a rebuild are dead by construction: the epoch only moves between steps, and every consumer
    !! overwrites (stage-1) or zeroes (s_amr_zero_fine_registers) before its first read of a step.
    impure subroutine s_amr_reg_prepare()

        integer :: g, kc, save_cur
        logical :: own_lo(3), own_hi(3), need

        if (.not. amr) return
        if (amr_reg_epoch_built == amr_mesh_epoch .and. amr_reg_nblk_built == amr_num_blocks) return
        save_cur = amr_cur
        amr_reg_of = 0
        amr_reg_n = 0
        do g = 1, amr_num_blocks
            need = amr_owns_all(g)
            if (.not. need .and. amr_block_level(g) <= 1) then
                ! (c) reflux-face participation without the seam clip: the clipped fills (coarse capture, L0/L1 apply, face
                ! wave) are a strict subset, so every rank that posts into the block's register slot is mapped
                call s_amr_select_slot(g)
                call s_amr_reflux_faces(amr_sidx, amr_ext, .false., own_lo, own_hi)
                need = any(own_lo .or. own_hi)
            end if
            if (need) then
                amr_reg_n = amr_reg_n + 1
                amr_reg_of(g) = amr_reg_n
            end if
        end do
        ! (b) children of owned blocks (f_amr_is_child, the fine-branch capture's test)
        do g = 1, amr_num_blocks
            if (.not. amr_owns_all(g)) cycle
            do kc = 1, amr_num_blocks
                if (amr_reg_of(kc) /= 0 .or. .not. f_amr_is_child(kc, g)) cycle
                amr_reg_n = amr_reg_n + 1
                amr_reg_of(kc) = amr_reg_n
            end do
        end do
        call s_amr_select_slot(save_cur)  ! also refreshes amr_reg_cur under the new map
        amr_reg_epoch_built = amr_mesh_epoch
        amr_reg_nblk_built = amr_num_blocks
        call s_amr_reg_reserve(amr_reg_n)
        call s_amr_zero_registers(.true., 1, amr_reg_n)
        call s_amr_zero_registers(.false., 1, amr_reg_n)

    end subroutine s_amr_reg_prepare

    !> The faces of the current block (amr_region_lo/hi, set on every rank by s_amr_select_slot) whose outside coarse layer lies in
    !! the coarse subdomain [sidx, sidx + ext] with transverse overlap: the faces that subdomain's rank refluxes. With clip,
    !! fine-fine seam faces are dropped (max_grid_size tiling: the outside cell is an adjacent block's interior, which the fine-fine
    !! halo already matches; refluxing there would corrupt that cell mid-step). Every reflux participant test derives from this one
    !! predicate, so senders and receivers agree by construction.
    pure subroutine s_amr_reflux_faces(sidx, ext, clip, s_lo, s_hi)

        integer, intent(in)  :: sidx(3), ext(3)
        logical, intent(in)  :: clip
        logical, intent(out) :: s_lo(3), s_hi(3)
        logical              :: tv(3), tvd
        integer              :: d

        tv = (amr_region_lo <= sidx + ext .and. amr_region_hi >= sidx) .or. .not. amr_dim
        s_lo = .false.; s_hi = .false.
        do d = 1, num_dims
            tvd = all(tv .or. [1, 2, 3] == d)
            s_lo(d) = tvd .and. amr_region_lo(d) - 1 >= sidx(d) .and. amr_region_lo(d) - 1 <= sidx(d) + ext(d)
            s_hi(d) = tvd .and. amr_region_hi(d) + 1 >= sidx(d) .and. amr_region_hi(d) + 1 <= sidx(d) + ext(d)
            if (clip) then
                if (s_lo(d) .and. f_amr_face_is_seam(d, -1)) s_lo(d) = .false.
                if (s_hi(d) .and. f_amr_face_is_seam(d, 1)) s_hi(d) = .false.
            end if
        end do

    end subroutine s_amr_reflux_faces

    !> This rank's reflux faces (seam-clipped) with its subdomain and the global transverse overlap [tlo, thi] per dim, so capture
    !! and apply share a block-relative frame aligned with the owner's freg. All true / full-block at np=1.
    impure subroutine s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)

        integer, intent(out) :: sidx(3), ext(3)
        logical, intent(out) :: own_lo(3), own_hi(3)
        integer, intent(out) :: tlo(3), thi(3)

        sidx = amr_sidx; ext = amr_ext
        tlo = max(amr_region_lo, sidx); thi = min(amr_region_hi, sidx + ext)
        call s_amr_reflux_faces(sidx, ext, .true., own_lo, own_hi)

    end subroutine s_amr_reflux_face_flags

    !> True iff the current block's face on `side` (+1 high / -1 low) in dim d is shared with an adjacent sub-block (max_grid_size
    !! tiling), i.e. another block's opposite face is exactly one cell away with matching transverse extents. Such a seam is
    !! fine-fine, not a c/f boundary. Reads the replicated block list (amr_region_*_all); no tiling means no match.
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

        if (.not. amr) return
        ! Registers on all ranks: regrid moves block faces, so any rank can participate (fine cells for freg; outside-face layer
        ! for creg capture/apply and for receiving freg from the block owner). freg is captured for the whole block and indexed
        ! block-relative by every applier, so registers must span a whole block. The largest block a rank can own is amr_maxc_fit
        ! (the scratch-constraint cap), so size to it; this matches m_amr's fine arrays.
        ! creg: relative 0-based transverse (0:maxc_t-1); freg: 0-based fine (0:max_f_t). Device-resident (@:ALLOCATE): capture
        ! and both applies run as kernels; no host copies read. The transverse extents are stashed so s_amr_reg_reserve can
        ! rebuild the same shapes when the slot dimension grows.
        rc = merge(maxc_fit, 1, amr_dim)
        rf = merge(amr_ref_ratio*rc - 1, 0, amr_dim)
        ! Start at a small slot capacity and grow on demand; do not size to amr_max_blocks (see amr_reg_cap above).
        amr_reg_cap = min(amr_max_blocks, amr_reg_floor)
        #:for D, A, B in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            if (amr_dim(${D}$)) then
                @:ALLOCATE(creg(${D}$)%lo(1:sys_size, 0:rc(${A}$) - 1, 0:rc(${B}$) - 1, 1:amr_reg_cap))
                @:ALLOCATE(creg(${D}$)%hi(1:sys_size, 0:rc(${A}$) - 1, 0:rc(${B}$) - 1, 1:amr_reg_cap))
                @:ALLOCATE(freg(${D}$)%lo(1:sys_size, 0:rf(${A}$), 0:rf(${B}$), 1:amr_reg_cap))
                @:ALLOCATE(freg(${D}$)%hi(1:sys_size, 0:rf(${A}$), 0:rf(${B}$), 1:amr_reg_cap))
            end if
        #:endfor
        ! per-slot geometry scratch for the batched capture kernels (device-resident: host-filled, GPU_UPDATE'd before each call)
        @:ALLOCATE(bjlo(1:amr_max_blocks), bjhi(1:amr_max_blocks), bo1(1:amr_max_blocks), bo2(1:amr_max_blocks))
        @:ALLOCATE(bt1lo(1:amr_max_blocks), bt1hi(1:amr_max_blocks), bt2lo(1:amr_max_blocks), bt2hi(1:amr_max_blocks))
        @:ALLOCATE(bclo(1:amr_max_blocks), bchi(1:amr_max_blocks), bactive(1:amr_max_blocks))
        bactive = .false.
        @:ALLOCATE(a_ol(1:amr_max_blocks), a_oh(1:amr_max_blocks), a_ta(1:amr_max_blocks), a_tb(1:amr_max_blocks))
        @:ALLOCATE(a_b1l(1:amr_max_blocks), a_b1h(1:amr_max_blocks), a_b2l(1:amr_max_blocks))
        @:ALLOCATE(a_b2h(1:amr_max_blocks), a_lo(1:amr_max_blocks), a_hi(1:amr_max_blocks), a_act(1:amr_max_blocks))
        @:ALLOCATE(a_mlo(1:amr_max_blocks), a_mhi(1:amr_max_blocks))
        ! participation-local register index (host-only ints; the register reals are what the dense map shrinks)
        allocate (amr_reg_of(1:amr_max_blocks))
        amr_reg_of = 0; amr_reg_n = 0; amr_reg_cur = 0
        amr_reg_epoch_built = -1_8; amr_reg_nblk_built = -1

    end subroutine s_initialize_amr_registers

    !> Block kc is a level+1 child of block g: one level finer with an overlapping region box (proper nesting makes overlap
    !! containment).
    pure logical function f_amr_is_child(kc, g) result(c)

        integer, intent(in) :: kc, g

        c = amr_block_level(kc) == amr_block_level(g) + 1 .and. all(amr_region_lo_all(:,kc) <= amr_region_hi_all(:, &
                            & g) .and. amr_region_hi_all(:,kc) >= amr_region_lo_all(:,g))

    end function f_amr_is_child

    !> Parent-fine footprint of block k inside its parent pblk, from replicated metadata only, so every rank computes the same box
    !! (a rank needs it for a block it does not own, whose own amr_isect_lo/hi is the empty non-owner footprint). Mirrors the
    !! level>=2 branch of s_set_amr_fine_geometry exactly; rr is the global amr_ref_ratio because a level>=2 block's parent is never
    !! an L0 tile (the only slot with a per-slot ratio of 1). Lives here rather than in m_amr so the child-creg capture below and
    !! m_amr's P2P gather/restrict/reflux share one copy of the formula ("use m_amr" would cycle).
    pure subroutine s_amr_parent_foot(k, pblk, plo, phi)

        integer, intent(in)  :: k, pblk
        integer, intent(out) :: plo(3), phi(3)

        plo = merge(amr_ref_ratio*(amr_region_lo_all(:,k) - amr_region_lo_all(:,pblk)), 0, amr_dim)
        phi = merge(amr_ref_ratio*(amr_region_hi_all(:,k) - amr_region_lo_all(:,pblk)) + amr_ref_ratio - 1, 0, amr_dim)

    end subroutine s_amr_parent_foot

    !> Batched member ibm of the current advance: select its slot and return its offset ko in the batch slab and its own extents
    !! ext. Outside a batch (amr_bat_n = 0) this is the one-block path, ko = 0.
    subroutine s_amr_bat_member(ibm, ko, ext)

        integer, intent(in)  :: ibm
        integer, intent(out) :: ko(3), ext(3)

        ko = 0
        if (amr_bat_n > 0) then
            call s_amr_select_slot(amr_bat_blk(ibm))
            ko(amr_bat_sd) = (ibm - 1)*amr_bat_w
            ext = amr_bat_mext(:,ibm)
        else
            ext = [m, n, p]
        end if

    end subroutine s_amr_bat_member

    !> Fill capture slot sreg with the face pair of one block along direction id: blo/bhi is the block's box in local flux indices
    !! (faces at blo(id)-1 and bhi(id)), wlo/whi the transverse window to capture (the block's own extent, or on the coarse side
    !! this rank's owned overlap), clo/chi the per-face gates. Transverse indices are 0-based from blo, aligned across the coarse
    !! and fine registers: fine children of block-relative coarse cell t are faces 2*t and 2*t+1.
    subroutine s_amr_capture_slot(sreg, id, blo, bhi, wlo, whi, clo, chi)

        integer, intent(in) :: sreg, id, blo(3), bhi(3), wlo(3), whi(3)
        logical, intent(in) :: clo, chi
        integer             :: ta, tb

        ta = merge(2, 1, id == 1); tb = merge(2, 3, id == 3)
        bactive(sreg) = .true.; bclo(sreg) = clo; bchi(sreg) = chi
        bjlo(sreg) = blo(id) - 1; bjhi(sreg) = bhi(id); bo1(sreg) = blo(ta); bo2(sreg) = blo(tb)
        bt1lo(sreg) = wlo(ta) - blo(ta); bt1hi(sreg) = whi(ta) - blo(ta)
        bt2lo(sreg) = wlo(tb) - blo(tb); bt2hi(sreg) = whi(tb) - blo(tb)
        bmax1 = max(bmax1, bt1hi(sreg)); bmax2 = max(bmax2, bt2hi(sreg))

    end subroutine s_amr_capture_slot

    !> Boundary-flux capture, batched over the register slots: for each active slot in [1:nb], reg(id)%lo/hi(eq, t1, t2, slot)
    !! [+=/=] cf * flux(face, bo1(slot)+t1, bo2(slot)+t2) for eq in [eqb:eqe], over the per-slot transverse window [bt1lo:bt1hi] x
    !! [bt2lo:bt2hi]; reg is freg (fine) or creg. acc=.true. accumulates, .false. overwrites (the merge picks the old value or 0
    !! with no arithmetic, so a stage-1 overwrite reads no uninitialized register). bclo/bchi gate the low/high face. The device
    !! kernel collapses (slot, t2, t1, eq) over the rectangular caps [0:maxt2]x[0:maxt1] (max over slots) and cycles inactive slots
    !! / out-of-window cells, so one launch replaces O(blocks) per-slot launches.
    impure subroutine s_amr_capture_batch(nb, id, fine, advective, cf, acc, maxt1, maxt2, eqb, eqe)

        integer, intent(in) :: nb, id, maxt1, maxt2, eqb, eqe
        !> Which flat Riemann buffer to read: T = flux_rsx_vf (advective), F = flux_src_rsx_vf (viscous, chemistry). Both are plain
        !! module arrays, so this routine takes no field dummies.
        logical, intent(in)  :: fine, advective, acc
        real(wp), intent(in) :: cf
        integer              :: eq, t1, t2, slot, i1, i2, i3, j1, j2, j3
        real(wp)             :: v

        $:GPU_PARALLEL_LOOP(collapse=4, private='[i1, i2, i3, j1, j2, j3, v]')
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
                        ! The flux reads must stay inside the bclo/bchi guards. A coarse slot goes active when either face is
                        ! owned, and the unowned face's index is still computed; it then points a whole block width outside this
                        ! rank's subdomain (jlo down to -amr_max_grid_size). Reading it unguarded is an out-of-bounds device access
                        ! against flux_rsx_vf's tight (-1:m_alloc) bounds. It hides at np=1, where the intersection is the block
                        ! and both flags hold, so single-rank tests do not catch it.
                        if (bclo(slot)) then
                            if (advective) then
                                v = flux_rsx_vf(i1, i2, i3, eq)
                            else
                                v = flux_src_rsx_vf(i1, i2, i3, eq)
                            end if
                            if (fine) then
                                freg(id)%lo(eq, t1, t2, slot) = merge(freg(id)%lo(eq, t1, t2, slot), 0._wp, acc) + cf*v
                            else
                                creg(id)%lo(eq, t1, t2, slot) = merge(creg(id)%lo(eq, t1, t2, slot), 0._wp, acc) + cf*v
                            end if
                        end if
                        if (bchi(slot)) then
                            if (advective) then
                                v = flux_rsx_vf(j1, j2, j3, eq)
                            else
                                v = flux_src_rsx_vf(j1, j2, j3, eq)
                            end if
                            if (fine) then
                                freg(id)%hi(eq, t1, t2, slot) = merge(freg(id)%hi(eq, t1, t2, slot), 0._wp, acc) + cf*v
                            else
                                creg(id)%hi(eq, t1, t2, slot) = merge(creg(id)%hi(eq, t1, t2, slot), 0._wp, acc) + cf*v
                            end if
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_capture_batch

    !> Total-flux capture of the filled slots (s_amr_capture_slot) into freg (fine) or creg, one batched kernel per category:
    !! advective (all equations); viscous (flux_src, mom..E); chemistry diffusion (flux_src species, plus the thermal-conduction +
    !! enthalpy energy flux only when not viscous, since the viscous pass already captured flux_src(E), which holds
    !! viscous+diffusion). The coarse and fine sides share this one routine, so the c/f reflux always subtracts matching fluxes.
    impure subroutine s_amr_capture_total(nb, id, fine, cf, acc)

        integer, intent(in)  :: nb, id
        logical, intent(in)  :: fine, acc
        real(wp), intent(in) :: cf

        if (.not. any(bactive(1:nb))) return
        $:GPU_UPDATE(device='[bjlo, bjhi, bo1, bo2, bt1lo, bt1hi, bt2lo, bt2hi, bclo, bchi, bactive]')
        call s_amr_capture_batch(nb, id, fine, .true., cf, acc, bmax1, bmax2, 1, sys_size)
        if (viscous) call s_amr_capture_batch(nb, id, fine, .false., cf, .true., bmax1, bmax2, eqn_idx%mom%beg, eqn_idx%E)
        if (chemistry .and. chem_params%diffusion) then
            call s_amr_capture_batch(nb, id, fine, .false., cf, .true., bmax1, bmax2, eqn_idx%species%beg, eqn_idx%species%end)
            if (.not. viscous) call s_amr_capture_batch(nb, id, fine, .false., cf, .true., bmax1, bmax2, eqn_idx%E, eqn_idx%E)
        end if
        bactive = .false.; bmax1 = 0; bmax2 = 0

    end subroutine s_amr_capture_total

    !> Capture the c/f boundary-face fluxes for direction id from the just-finalized flux array. Runs inside s_compute_rhs: coarse
    !! call (amr_in_fine_advance false, coarse globals) fills creg at block boundary faces; fine call (flag true, globals swapped to
    !! the fine block) fills freg at fine faces -1 and m/n/p. creg uses relative 0-based transverse; freg uses 0-based fine.
    impure subroutine s_amr_capture_boundary_flux(id, stage)

        integer, intent(in) :: id, stage
        integer             :: islot, save_cur, kc, ibm, ko(3), ext(3), sidx(3), tlo(3), thi(3), cflo(3), cfhi(3)
        logical             :: own_lo(3), own_hi(3), accum
        real(wp)            :: coef

        if (.not. amr) return
        ! Refresh the participation map + register capacity on a topology change; no-op (two integer compares) otherwise.
        call s_amr_reg_prepare()
        if (igr) return  ! stage-1 IGR coupling is restriction-only: the fused IGR flux kernels do not expose face fluxes to capture
        if (amr_in_fine_advance .and. .not. amr_rank_owns_block) return
        ! a level-0 L0 tile advancing through the fine path is coarse, not a fine block: skip the freg self-capture and the
        ! parent-of-level-1 child-creg loop (which would overwrite the real fine block's creg in the tile-swapped frame). Its creg
        ! comes from the dedicated L0 coarse RHS (amr_in_fine_advance=F). Pure-AMR has no level-0 slots so this never fires.
        if (amr_in_fine_advance .and. amr_block_level(amr_cur) == 0) return
        ! flux data was just written by device kernels; the face reads below run as device kernels too
        save_cur = amr_cur
        if (amr_in_fine_advance) then
            ! fine branch: globals swapped; faces -1 and the fine extent in direction id. Batched advance (amr_bat_n > 0): the slab
            ! holds amr_bat_n same-extent blocks, member ibm at offset ko along amr_bat_sd; each member's faces go into its own
            ! register slot. L2->L1 reflux: the parent is already RK-updated by reflux time, so freg must hold the rk3_w-weighted
            ! step-integral flux for the once-per-step state correction (stage 1 overwrites = implicit zero, cf. creg below).
            if (amr_block_level(amr_cur) >= 2) then
                coef = rk3_w(stage); accum = (stage > 1)
            else
                coef = 1._wp; accum = .false.  ! overwrite each stage (default)
            end if
            do ibm = 1, max(1, amr_bat_n)
                call s_amr_bat_member(ibm, ko, ext)
                call s_amr_capture_slot(amr_reg_cur, id, ko, ko + ext, ko, ko + ext, .true., .true.)
            end do
            call s_amr_capture_total(amr_reg_n, id, .true., coef, accum)
            ! multi-level lock-step: this fine block is the coarse side (parent) of its level+1 children. Capture creg for each
            ! child from this block's fine flux at the child's footprint faces; the footprint is in this parent's fine frame, so
            ! it indexes flux_rsx_vf directly. creg holds the rk3_w-weighted step-integral flux for the once-per-step state reflux
            ! into this parent (s_amr_reflux_to_parent). creg is the parent's own flux, so the parent owner captures it for every
            ! child of this block, including children owned by another rank, which supply only the matching freg
            ! (s_amr_restrict_wave). Framing therefore comes from s_amr_parent_foot (replicated metadata), not
            ! amr_isect_*_all(:,kc), which is the empty sentinel for a child this rank does not own. Each child's creg lives at
            ! its dense register slot (a child of an owned block is always mapped, s_amr_reg_prepare clause (b) is this loop's
            ! twin); both faces always owned (the parent spans the whole child footprint).
            do ibm = 1, max(1, amr_bat_n)
                call s_amr_bat_member(ibm, ko, ext)
                do kc = 1, amr_num_blocks
                    if (.not. f_amr_is_child(kc, amr_cur)) cycle
                    call s_amr_parent_foot(kc, amr_cur, cflo, cfhi)
                    call s_amr_capture_slot(amr_reg_of(kc), id, cflo + ko, cfhi + ko, cflo + ko, cfhi + ko, .true., .true.)
                end do
            end do
            if (amr_bat_n > 0) call s_amr_select_slot(save_cur)
            call s_amr_capture_total(amr_reg_n, id, .false., rk3_w(stage), stage > 1)
        else
            ! coarse branch: a face's capture runs on the rank owning the coarse cells just outside it (its flux_rsx_vf covers
            ! that face; at a rank-interior face the same rank also holds the inside cells). This rank fills creg over its owned
            ! transverse overlap [tlo:thi] in the block-relative frame. At np=1 the intersection is the block and both flags hold,
            ! recovering single-rank behavior exactly. One coarse s_compute_rhs pass fills every active block's registers: revisit
            ! each slot's region+intersection in turn. own_lo/own_hi is s_amr_reg_prepare's clause (c) verbatim, so
            ! amr_reg_of(islot) is always mapped here.
            do islot = 1, amr_num_blocks
                ! a level>=2 block's coarse side is its parent (creg captured in the fine branch), not L0
                if (amr_block_level(islot) >= 2) cycle
                call s_amr_select_slot(islot)
                call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
                if (own_lo(id) .or. own_hi(id)) call s_amr_capture_slot(amr_reg_of(islot), id, amr_region_lo - sidx, &
                    & amr_region_hi - sidx, tlo - sidx, thi - sidx, own_lo(id), own_hi(id))
            end do
            call s_amr_select_slot(save_cur)
            call s_amr_capture_total(amr_reg_n, id, .false., 1._wp, .false.)
        end if

    end subroutine s_amr_capture_boundary_flux

    !> Correct the coarse rhs in the first cell outside each block face so the coarse update sees the (child-averaged) fine flux at
    !! every c/f face. Signs follow rhs = (flux_left - flux_right)/dx: low face is the outside cell's right face => rhs += (F_coarse
    !! - Fbar_fine)/dx; high face is the outside cell's left face => rhs += (Fbar_fine - F_coarse)/dx. Cells inside the block need
    !! no correction (end-of-step restriction overwrites them). c1/c2 are relative 0-based coarse transverse indices.
    impure subroutine s_amr_apply_reflux(rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: eq, c1, c2, c1w, c2w, k, save_cur, nact, gmax1, gmax2
        integer                                                :: f10, f20, dd1, dd2, nch, rr, dd1_hi, dd2_hi, sreg
        integer                                                :: bla, bha, blb, bhb
        integer                                                :: i2, i3, sidx(3), ext(3), tlo(3), thi(3)
        logical                                                :: dta, dtb, own_lo(3), own_hi(3)
        real(wp)                                               :: fblo, fbhi

        if (.not. amr) return
        ! Refresh the participation map + register capacity on a topology change; no-op (two integer compares) otherwise.
        call s_amr_reg_prepare()
        if (igr) return  ! stage-1 IGR: restriction-only coupling (no captured fluxes)
        rr = amr_ref_ratio
        save_cur = amr_cur

        ! Batched over the level-1 blocks, one kernel per face direction (mirror of the capture-side batching,
        ! s_amr_capture_batch), since per-launch overhead rather than arithmetic dominates a per-block form.
        ! Block corrections are disjoint (the merge invariant keeps blocks >= buff_size apart), so the batched kernel is
        ! equivalent to a per-block loop. Host precompute walks the slots with select_slot + s_amr_reflux_face_flags; the
        ! a_* descriptors are pushed once per direction. Per direction d the transverse dims are (ta, tb) and the fine-face
        ! register holds rr children per active transverse dim.
        #:for D, TA, TB, DX, IDX in [(1, 2, 3, 'dx', 'a_ol(k), i2, i3'), (2, 1, 3, 'dy', 'i2, a_ol(k), i3'), (3, 1, 2, 'dz', &
                                      & 'i2, i3, a_ol(k)')]
            if (amr_dim(${D}$)) then
                dta = amr_dim(${TA}$); dtb = amr_dim(${TB}$)
                nch = rr**count([dta, dtb])
                dd1_hi = merge(rr - 1, 0, dta); dd2_hi = merge(rr - 1, 0, dtb)
                nact = 0; gmax1 = 0; gmax2 = 0
                a_act = .false.
                do k = 1, amr_num_blocks
                    if (amr_block_level(k) /= 1) cycle
                    call s_amr_select_slot(k)
                    call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
                    if (.not. (own_lo(${D}$) .or. own_hi(${D}$))) cycle
                    bla = tlo(${TA}$) - amr_region_lo(${TA}$); bha = thi(${TA}$) - amr_region_lo(${TA}$)
                    blb = tlo(${TB}$) - amr_region_lo(${TB}$); bhb = thi(${TB}$) - amr_region_lo(${TB}$)
                    sreg = amr_reg_of(k)
                    a_act(sreg) = .true.; a_lo(sreg) = own_lo(${D}$); a_hi(sreg) = own_hi(${D}$)
                    a_ol(sreg) = amr_region_lo(${D}$) - 1 - sidx(${D}$); a_oh(sreg) = amr_region_hi(${D}$) + 1 - sidx(${D}$)
                    a_ta(sreg) = amr_region_lo(${TA}$) - sidx(${TA}$); a_tb(sreg) = amr_region_lo(${TB}$) - sidx(${TB}$)
                    a_b1l(sreg) = bla; a_b1h(sreg) = bha; a_b2l(sreg) = blb; a_b2h(sreg) = bhb
                    a_mlo(sreg) = 1._wp; a_mhi(sreg) = 1._wp
                    if (own_lo(${D}$)) a_mlo(sreg) = ${DX}$(a_ol(sreg))
                    if (own_hi(${D}$)) a_mhi(sreg) = ${DX}$(a_oh(sreg))
                    nact = nact + 1
                    gmax1 = max(gmax1, bha - bla); gmax2 = max(gmax2, bhb - blb)
                end do
                call s_amr_select_slot(save_cur)
                if (nact > 0) then
                    $:GPU_UPDATE(device='[a_ol, a_oh, a_ta, a_tb, a_b1l, a_b1h, a_b2l, a_b2h, a_lo, a_hi, a_act, a_mlo, a_mhi]')
                    $:GPU_PARALLEL_LOOP(collapse=4, private='[c1, c2, f10, f20, dd1, dd2, fblo, fbhi, i2, i3]')
                    do k = 1, amr_reg_n
                        do c2w = 0, gmax2
                            do c1w = 0, gmax1
                                do eq = 1, sys_size
                                    if (.not. a_act(k)) cycle
                                    c1 = a_b1l(k) + c1w; c2 = a_b2l(k) + c2w
                                    if (c1 > a_b1h(k) .or. c2 > a_b2h(k)) cycle
                                    f10 = 0; if (dta) f10 = rr*c1
                                    f20 = 0; if (dtb) f20 = rr*c2
                                    fblo = 0._wp; fbhi = 0._wp
                                    do dd2 = 0, dd2_hi
                                        do dd1 = 0, dd1_hi
                                            fblo = fblo + freg(${D}$)%lo(eq, f10 + dd1, f20 + dd2, k)
                                            fbhi = fbhi + freg(${D}$)%hi(eq, f10 + dd1, f20 + dd2, k)
                                        end do
                                    end do
                                    fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                                    i2 = a_ta(k) + c1; i3 = a_tb(k) + c2
                                    if (a_lo(k)) rhs_vf(eq)%sf(${IDX}$) = rhs_vf(eq)%sf(${IDX}$) + (creg(${D}$)%lo(eq, c1, c2, &
                                        & k) - fblo)/a_mlo(k)
                                    if (a_hi(k)) rhs_vf(eq)%sf(${IDX.replace('a_ol', 'a_oh')}$) &
                                        & = rhs_vf(eq)%sf(${IDX.replace('a_ol', 'a_oh')}$) + (fbhi - creg(${D}$)%hi(eq, c1, c2, &
                                        & k))/a_mhi(k)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        #:endfor

    end subroutine s_amr_apply_reflux

    !> Zero the working block's fine registers.
    impure subroutine s_amr_zero_fine_registers()

        if (.not. amr) return
        ! Refresh the participation map + register capacity on a topology change; no-op (two integer compares) otherwise.
        call s_amr_reg_prepare()
        if (igr) return  ! stage-1 IGR: restriction-only coupling (no captured fluxes)
        if (.not. amr_rank_owns_block) return
        call s_amr_zero_registers(.true., amr_reg_cur, amr_reg_cur)

    end subroutine s_amr_zero_fine_registers

    !> Zero register slots slo..shi of freg (fine) or creg in every active direction.
    impure subroutine s_amr_zero_registers(fine, slo, shi)

        logical, intent(in) :: fine
        integer, intent(in) :: slo, shi
        integer             :: d, ta, tb, eq, t1, t2, t1_hi, t2_hi, islot

        if (shi < slo) return
        do d = 1, num_dims
            ta = merge(2, 1, d == 1); tb = merge(2, 3, d == 3)
            t1_hi = merge(rf(ta), rc(ta) - 1, fine); t2_hi = merge(rf(tb), rc(tb) - 1, fine)
            $:GPU_PARALLEL_LOOP(collapse=4)
            do islot = slo, shi
                do t2 = 0, t2_hi
                    do t1 = 0, t1_hi
                        do eq = 1, sys_size
                            if (fine) then
                                freg(d)%lo(eq, t1, t2, islot) = 0._wp
                                freg(d)%hi(eq, t1, t2, islot) = 0._wp
                            else
                                creg(d)%lo(eq, t1, t2, islot) = 0._wp
                                creg(d)%hi(eq, t1, t2, islot) = 0._wp
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_amr_zero_registers

    !> Shared Berger-Colella state reflux kernel: apply q(outside) += w*dtl*(F_coarse - Fbar_fine)/m on the low face and +=
    !! w*dtl*(Fbar_fine - F_coarse)/m on the high face for each active dim, where F_coarse is creg and Fbar_fine averages freg over
    !! the rr**(ndim-1) covering fine faces. Used by both s_amr_apply_reflux_state (L0/L1, coarse/sidx frame, unit weights from
    !! ownership, rr=2) and s_amr_reflux_to_parent (L2->L1, parent-fine frame, sibling-seam weights, rr=amr_ref_ratio). All framing
    !! is caller-passed so the flux-correction math is single-sourced: islot = dense register slot (amr_reg_cur); rr = refinement
    !! ratio (fine faces per coarse face per transverse dim); dtl = reflux dt; olo/ohi(d) = outside coarse-cell index just
    !! below/above the block face in dim d; glo/ghi(d) = creg-local loop range in dim d (transverse for the two faces d' /= d);
    !! woff(d) = transverse write origin so the cell index is woff(d) + g; w_lo/w_hi(d) = per-face weight (0 skips the write:
    !! unowned face at np>1, or a fine-fine sibling-tile seam); mlo/mhi(d) = outside-cell width for the low/high face
    !! (invalid/unused where weight is 0). A zero weight skips the write (not multiply-by-0) because the outside index may be out of
    !! bounds on an unowned face.
    impure subroutine s_amr_reflux_apply_faces(q, islot, rr, dtl, olo, ohi, glo, ghi, woff, w_lo, w_hi, mlo, mhi)

        type(scalar_field), dimension(sys_size), intent(inout) :: q
        integer, intent(in) :: islot, rr, olo(3), ohi(3), glo(3), ghi(3), woff(3)
        real(wp), intent(in) :: dtl, w_lo(3), w_hi(3), mlo(3), mhi(3)
        integer :: eq, g1, g2, f10, f20, dd1, dd2, nch, dd1_hi, dd2_hi, ol, oh, wa, wb, gla, gha, glb, ghb
        logical :: dta, dtb
        real(wp) :: fblo, fbhi, wl, wh, ml, mh

        ! Per face direction d the transverse dims are (ta, tb). Loop bounds and framing are hoisted to scalars: array-element
        ! bounds (glo(d)/ghi(d)) would force the host arrays present on the device (an ACC present error).

        #:for D, TA, TB, IDX in [(1, 2, 3, 'ol, wa + g1, wb + g2'), (2, 1, 3, 'wa + g1, ol, wb + g2'), (3, 1, 2, &
                                  & 'wa + g1, wb + g2, ol')]
            if (amr_dim(${D}$) .and. (w_lo(${D}$) /= 0._wp .or. w_hi(${D}$) /= 0._wp)) then
                dta = amr_dim(${TA}$); dtb = amr_dim(${TB}$)
                nch = rr**count([dta, dtb]); dd1_hi = merge(rr - 1, 0, dta); dd2_hi = merge(rr - 1, 0, dtb)
                gla = glo(${TA}$); gha = ghi(${TA}$); glb = glo(${TB}$); ghb = ghi(${TB}$); wa = woff(${TA}$); wb = woff(${TB}$)
                ol = olo(${D}$); oh = ohi(${D}$); wl = w_lo(${D}$); wh = w_hi(${D}$); ml = mlo(${D}$); mh = mhi(${D}$)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
                do eq = 1, sys_size
                    do g2 = glb, ghb
                        do g1 = gla, gha
                            f10 = 0; if (dta) f10 = rr*g1
                            f20 = 0; if (dtb) f20 = rr*g2
                            fblo = 0._wp; fbhi = 0._wp
                            do dd2 = 0, dd2_hi
                                do dd1 = 0, dd1_hi
                                    fblo = fblo + freg(${D}$)%lo(eq, f10 + dd1, f20 + dd2, islot)
                                    fbhi = fbhi + freg(${D}$)%hi(eq, f10 + dd1, f20 + dd2, islot)
                                end do
                            end do
                            fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                            if (wl /= 0._wp) q(eq)%sf(${IDX}$) = q(eq)%sf(${IDX}$) + wl*dtl*(creg(${D}$)%lo(eq, g1, g2, &
                                & islot) - fblo)/ml
                            if (wh /= 0._wp) q(eq)%sf(${IDX.replace('ol', 'oh')}$) = q(eq)%sf(${IDX.replace('ol', 'oh')}$) &
                                & + wh*dtl*(fbhi - creg(${D}$)%hi(eq, g1, g2, islot))/mh
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_amr_reflux_apply_faces

    impure subroutine s_finalize_amr_registers()

        integer :: d

        if (.not. amr) return
        do d = 1, num_dims
            @:DEALLOCATE(creg(d)%lo, creg(d)%hi, freg(d)%lo, freg(d)%hi)
        end do
        @:DEALLOCATE(bjlo, bjhi, bo1, bo2, bt1lo, bt1hi, bt2lo, bt2hi, bclo, bchi, bactive)
        @:DEALLOCATE(a_ol, a_oh, a_ta, a_tb, a_b1l, a_b1h, a_b2l, a_b2h, a_lo, a_hi, a_act, a_mlo, a_mhi)
        if (allocated(amr_reg_of)) deallocate (amr_reg_of)
        amr_reg_n = 0; amr_reg_cur = 0

    end subroutine s_finalize_amr_registers

end module m_amr_registers
