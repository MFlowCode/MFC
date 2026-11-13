!>
!!@file
!!@brief Contains module m_amr_registers

#! AMD OpenMP lane: kernels here may only name arrays allocated on every path that reaches them.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief AMR flux registers: per-RK-stage refluxing at the coarse/fine block boundary. Depends only on m_derived_types and
!! m_global_parameters so m_rhs (capture) and m_time_steppers (apply) can use it without cycles. creg is indexed 0-based relative to
!! the rank's block intersection; freg 0-based local fine (children of coarse cell t are 2*t, 2*t+1). The slot dimension grows on
!! demand (s_amr_reg_reserve). The captured flux is the 1:sys_size advective flux of flux_rsx_vf plus the viscous, species-diffusion
!! and conduction fluxes of flux_src_rsx_vf; the non-conservative alpha*div(u) remainder is not captured (held by mpp_lim).
module m_amr_registers

    use m_derived_types
    use m_global_parameters
    use m_riemann_state, only: flux_rsx_vf, flux_src_rsx_vf

    implicit none

    private; public :: s_initialize_amr_registers, s_amr_capture_boundary_flux, s_amr_apply_reflux, s_finalize_amr_registers, &
        & f_amr_block_near, s_amr_reflux_face_flags, s_amr_reflux_faces, s_amr_reflux_apply_faces, s_amr_parent_foot, &
        & s_amr_reg_prepare, freg, creg, f_amr_face_is_seam

    !> SSP-RK3 effective flux weights: q^{n+1} = q^n + dt*(L(q^n)/6 + L(q^(1))/6 + 2*L(q^(2))/3).
    real(wp), parameter :: rk3_w(3) = [1._wp/6._wp, 1._wp/6._wp, 2._wp/3._wp]

    !> Registers for the two block faces normal to one direction: (1:sys_size, transverse-1, transverse-2, slot).
    type t_face_reg
        real(wp), allocatable :: lo(:,:,:,:)
        real(wp), allocatable :: hi(:,:,:,:)
    end type t_face_reg

    type(t_face_reg) :: creg(3)  !< coarse flux at block boundary faces (relative 0-based transverse)
    type(t_face_reg) :: freg(3)  !< fine flux at covering fine faces (0-based fine transverse)
    $:GPU_DECLARE(create='[creg, freg]')

    !> Slot capacity the registers are sized for (grows geometrically; amr_max_blocks is only a cap) and the transverse extents.
    integer, parameter :: amr_reg_floor = 64  !< initial slot capacity; growth doubles from here
    integer            :: amr_reg_cap = 0
    integer            :: rc(3), rf(3)        !< coarse / fine transverse extents (creg is 0:rc-1, freg is 0:rf)
    !> Mesh-epoch and block-count keys of the last s_amr_reg_prepare build.
    integer(8) :: amr_reg_epoch_built = -1_8
    integer    :: amr_reg_nblk_built = -1

    !> Per-slot capture tables, host-filled then GPU_UPDATE'd once. Rows: 1:2 normal-face flux indices (jlo, jhi), 3:4 transverse
    !! origins, 5:8 transverse window (t1lo, t1hi, t2lo, t2hi), 9:10 low/high face gates, 11 slot gate. Set 0 is the fine side's
    !! per-call table; sets 1:3 are the coarse per-direction tables (s_amr_coarse_tables). bmax(:, set) holds each set's caps.
    integer, allocatable :: bcap(:,:,:)
    integer              :: bmax(2,0:3) = 0
    $:GPU_DECLARE(create='[bcap]')

    !> Per-slot reflux apply tables, built with the coarse capture tables. aint rows: 1:2 outside coarse-cell index in the face dim
    !! (low, high), 3:4 transverse origins, 5:8 transverse windows, 9:10 per-face applies, 11 slot gate; areal rows 1:2 outside-cell
    !! widths. a_n/a_m hold each direction's active count and caps.
    integer, allocatable  :: aint(:,:,:)
    real(wp), allocatable :: areal(:,:,:)
    integer               :: a_n(3) = 0, a_m(2, 3) = 0
    $:GPU_DECLARE(create='[aint, areal]')
    !> Mesh-epoch and block-count keys of the last coarse-table build.
    integer(8) :: amr_ctab_epoch_built = -1_8
    integer    :: amr_ctab_nblk_built = -1

contains
    !> Grow the registers to at least nslot slots (doubling, never shrinking), preserving contents: freg accumulates across RK
    !! stages. Growth stages through a device temporary like the store, so the host mirror is undefined afterwards.
    impure subroutine s_amr_reg_reserve(nslot)

        integer, intent(in) :: nslot
        integer             :: oldcap, newcap
        integer             :: eq, t1, t2, c4
        !> slot count above which growth uses the host round trip instead of device staging
        integer, parameter    :: amr_reg_grow_dev_cap = 512
        real(wp), allocatable :: rtmp(:,:,:,:), rtmp_d(:,:,:,:)

        if (nslot <= amr_reg_cap) return
        oldcap = amr_reg_cap
        newcap = min(amr_max_blocks, max(2*oldcap, nslot))

        #:for D, A, B in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            #:for R, OFF in [('c', ' - 1'), ('f', '')]
                #:for S in ['lo', 'hi']
                    #:set ARR = R + 'reg(' + str(D) + ')%' + S
                    #:set U2 = 'r' + R + '(' + str(A) + ')' + OFF
                    #:set U3 = 'r' + R + '(' + str(B) + ')' + OFF
                    if (allocated(${ARR}$)) then
                        if (oldcap > amr_reg_grow_dev_cap) then
                            ! above the threshold stage through the host: device staging would transiently hold 2*oldcap slots
                            $:GPU_UPDATE(host='[' + ARR + ']')
                            allocate (rtmp(1:sys_size,0:${U2}$,0:${U3}$,1:oldcap))
                            rtmp = ${ARR}$(:,:,:,1:oldcap)
                            @:DEALLOCATE(${ARR}$)
                            @:ALLOCATE(${ARR}$(1:sys_size, 0:${U2}$, 0:${U3}$, 1:newcap))
                            ${ARR}$ = 0._wp
                            ${ARR}$(:,:,:,1:oldcap) = rtmp
                            deallocate (rtmp)
                            $:GPU_UPDATE(device='[' + ARR + ']')
                        else
                            ! stage the live slots on the device
                            @:ALLOCATE(rtmp_d(1:sys_size, 0:${U2}$, 0:${U3}$, 1:oldcap))
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do c4 = 1, oldcap
                                do t2 = 0, ${U3}$
                                    do t1 = 0, ${U2}$
                                        do eq = 1, sys_size
                                            rtmp_d(eq, t1, t2, c4) = ${ARR}$(eq, t1, t2, c4)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                            @:DEALLOCATE(${ARR}$)
                            @:ALLOCATE(${ARR}$(1:sys_size, 0:${U2}$, 0:${U3}$, 1:newcap))
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do c4 = 1, oldcap
                                do t2 = 0, ${U3}$
                                    do t1 = 0, ${U2}$
                                        do eq = 1, sys_size
                                            ${ARR}$(eq, t1, t2, c4) = rtmp_d(eq, t1, t2, c4)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                            @:DEALLOCATE(rtmp_d)
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do c4 = oldcap + 1, newcap
                                do t2 = 0, ${U3}$
                                    do t1 = 0, ${U2}$
                                        do eq = 1, sys_size
                                            ${ARR}$(eq, t1, t2, c4) = 0._wp
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    end if
                #:endfor
            #:endfor
        #:endfor

        amr_reg_cap = newcap

    end subroutine s_amr_reg_reserve

    !> Build the dense participation-local register index (amr_reg_of/amr_reg_n) and size and zero the registers. Keyed on the mesh
    !! epoch (block-count change is the tripwire). Block g maps iff this rank owns g or g's parent, or reflux-face-participates in g
    !! per s_amr_reflux_faces, so the map covers every capture, receive and apply that indexes g's slot. Dense slots alias across
    !! rebuilds, hence the zeroing; contents are dead at a rebuild because stage 1 overwrites them.
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
                ! unclipped, so the seam-clipped participants are a subset
                call s_amr_select_slot(g)
                call s_amr_reflux_faces(amr_sidx, amr_ext, .false., own_lo, own_hi)
                need = any(own_lo .or. own_hi)
            end if
            if (need) then
                amr_reg_n = amr_reg_n + 1
                amr_reg_of(g) = amr_reg_n
            end if
        end do
        ! children of owned blocks (the fine-branch capture's test)
        do g = 1, amr_num_blocks
            if (.not. amr_owns_all(g)) cycle
            do kc = 1, amr_num_blocks
                if (amr_reg_of(kc) /= 0 .or. .not. f_amr_is_child(kc, g)) cycle
                amr_reg_n = amr_reg_n + 1
                amr_reg_of(kc) = amr_reg_n
            end do
        end do
        call s_amr_select_slot(save_cur)  ! refreshes amr_reg_cur under the new map
        amr_reg_epoch_built = amr_mesh_epoch
        amr_reg_nblk_built = amr_num_blocks
        call s_amr_reg_reserve(amr_reg_n)
        call s_amr_zero_registers(.true., 1, amr_reg_n)
        call s_amr_zero_registers(.false., 1, amr_reg_n)

    end subroutine s_amr_reg_prepare

    !> Faces of the current block whose outside coarse layer lies in the coarse subdomain [sidx, sidx + ext] with transverse overlap
    !! (the faces that rank refluxes). With clip, fine-fine seam faces are dropped. Every reflux participant test derives from this.
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
            ! nested: the O(blocks) seam search runs only for faces that survived the box tests
            if (clip) then
                if (s_lo(d)) then
                    if (f_amr_face_is_seam(d, -1)) s_lo(d) = .false.
                end if
                if (s_hi(d)) then
                    if (f_amr_face_is_seam(d, 1)) s_hi(d) = .false.
                end if
            end if
        end do

    end subroutine s_amr_reflux_faces

    !> False only if block k's region widened by its outside face layer misses the coarse subdomain [sidx, sidx + ext].
    pure logical function f_amr_block_near(k, sidx, ext) result(near)

        integer, intent(in) :: k, sidx(3), ext(3)

        near = all((amr_region_lo_all(:,k) - 1 <= sidx + ext .and. amr_region_hi_all(:,k) + 1 >= sidx) .or. .not. amr_dim)

    end function f_amr_block_near

    !> This rank's seam-clipped reflux faces and the global transverse overlap [tlo, thi], in the block-relative frame of the
    !! owner's freg.
    impure subroutine s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)

        integer, intent(out) :: sidx(3), ext(3)
        logical, intent(out) :: own_lo(3), own_hi(3)
        integer, intent(out) :: tlo(3), thi(3)

        sidx = amr_sidx; ext = amr_ext
        tlo = max(amr_region_lo, sidx); thi = min(amr_region_hi, sidx + ext)
        call s_amr_reflux_faces(sidx, ext, .true., own_lo, own_hi)

    end subroutine s_amr_reflux_face_flags

    !> True iff the current block's face on `side` (+1 high / -1 low) in dim d abuts an adjacent block with matching transverse
    !! extents (a fine-fine seam, not a c/f boundary).
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

        integer, intent(in) :: maxc_fit(3)  !< amr_maxc_fit from m_amr (largest block a rank can own)

        if (.not. amr) return
        ! All ranks hold registers (regrid moves faces) spanning a whole block, up to amr_maxc_fit; creg is 0:rc-1, freg 0:rf.
        rc = merge(maxc_fit, 1, amr_dim)
        rf = merge(amr_ref_ratio*rc - 1, 0, amr_dim)
        amr_reg_cap = min(amr_max_blocks, amr_reg_floor)
        #:for D, A, B in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            if (amr_dim(${D}$)) then
                @:ALLOCATE(creg(${D}$)%lo(1:sys_size, 0:rc(${A}$) - 1, 0:rc(${B}$) - 1, 1:amr_reg_cap))
                @:ALLOCATE(creg(${D}$)%hi(1:sys_size, 0:rc(${A}$) - 1, 0:rc(${B}$) - 1, 1:amr_reg_cap))
                @:ALLOCATE(freg(${D}$)%lo(1:sys_size, 0:rf(${A}$), 0:rf(${B}$), 1:amr_reg_cap))
                @:ALLOCATE(freg(${D}$)%hi(1:sys_size, 0:rf(${A}$), 0:rf(${B}$), 1:amr_reg_cap))
            end if
        #:endfor
        @:ALLOCATE(bcap(1:11, 1:amr_max_blocks, 0:3))
        bcap = 0; bmax = 0
        @:ALLOCATE(aint(1:11, 1:amr_max_blocks, 1:3), areal(1:2, 1:amr_max_blocks, 1:3))
        aint = 0; areal = 1._wp; amr_ctab_epoch_built = -1_8; amr_ctab_nblk_built = -1
        allocate (amr_reg_of(1:amr_max_blocks))
        amr_reg_of = 0; amr_reg_n = 0; amr_reg_cur = 0
        amr_reg_epoch_built = -1_8; amr_reg_nblk_built = -1

    end subroutine s_initialize_amr_registers

    !> Block kc is a level+1 child of block g: one level finer with an overlapping region box.
    pure logical function f_amr_is_child(kc, g) result(c)

        integer, intent(in) :: kc, g

        c = amr_block_level(kc) == amr_block_level(g) + 1 .and. all(amr_region_lo_all(:,kc) <= amr_region_hi_all(:, &
                            & g) .and. amr_region_hi_all(:,kc) >= amr_region_lo_all(:,g))

    end function f_amr_is_child

    !> Parent-fine footprint of block k inside its parent pblk from replicated metadata, so every rank gets the same box (the
    !! block's own amr_isect_lo/hi is empty on non-owners). Lives here because "use m_amr" would cycle.
    pure subroutine s_amr_parent_foot(k, pblk, plo, phi)

        integer, intent(in)  :: k, pblk
        integer, intent(out) :: plo(3), phi(3)

        plo = merge(amr_ref_ratio*(amr_region_lo_all(:,k) - amr_region_lo_all(:,pblk)), 0, amr_dim)
        phi = merge(amr_ref_ratio*(amr_region_hi_all(:,k) - amr_region_lo_all(:,pblk)) + amr_ref_ratio - 1, 0, amr_dim)

    end subroutine s_amr_parent_foot

    !> Select batch member ibm's slot; returns its offset ko in the batch slab and its extents ext (ko = 0 outside a batch).
    subroutine s_amr_bat_member(ibm, ko, ext)

        integer, intent(in)  :: ibm
        integer, intent(out) :: ko(3), ext(3)

        ko = 0
        if (amr_bat_n > 0) then
            call s_amr_select_slot(amr_bat_blk(ibm))
            ko = amr_bat_off(:,ibm)
            ext = amr_bat_mext(:,ibm)
        else
            ext = [m, n, p]
        end if

    end subroutine s_amr_bat_member

    !> Fill capture table `set`, slot sreg, for one block's face pair along id: blo/bhi is the block box in local flux indices
    !! (faces at blo(id)-1 and bhi(id)), wlo/whi the transverse window, clo/chi the per-face gates. Transverse indices are 0-based
    !! from blo, aligned across the coarse and fine registers.
    subroutine s_amr_capture_slot(set, sreg, id, blo, bhi, wlo, whi, clo, chi)

        integer, intent(in) :: set, sreg, id, blo(3), bhi(3), wlo(3), whi(3)
        logical, intent(in) :: clo, chi
        integer             :: ta, tb

        ta = merge(2, 1, id == 1); tb = merge(2, 3, id == 3)
        bcap(:,sreg, set) = [blo(id) - 1, bhi(id), blo(ta), blo(tb), wlo(ta) - blo(ta), whi(ta) - blo(ta), wlo(tb) - blo(tb), &
             & whi(tb) - blo(tb), merge(1, 0, clo), merge(1, 0, chi), 1]
        bmax(:,set) = max(bmax(:,set), bcap(6:8:2,sreg, set))

    end subroutine s_amr_capture_slot

    !> Batched boundary-flux capture: reg(id)%lo/hi(eq, t1, t2, slot) [+=/=] cf * flux over each slot's window (bcap rows 5:8) for
    !! eq in [eqb:eqe]; reg is freg (fine) or creg. acc selects accumulate vs overwrite; advective selects flux_rsx_vf vs
    !! flux_src_rsx_vf. Inactive slots and out-of-window cells are cycled.
    impure subroutine s_amr_capture_batch(set, nb, id, fine, advective, cf, acc, maxt1, maxt2, eqb, eqe)

        integer, intent(in)  :: set, nb, id, maxt1, maxt2, eqb, eqe
        logical, intent(in)  :: fine, advective, acc
        real(wp), intent(in) :: cf
        integer              :: eq, t1, t2, slot, i1, i2, i3, j1, j2, j3
        real(wp)             :: v

        $:GPU_PARALLEL_LOOP(collapse=4, private='[i1, i2, i3, j1, j2, j3, v]')
        do slot = 1, nb
            do t2 = 0, maxt2
                do t1 = 0, maxt1
                    do eq = eqb, eqe
                        if (bcap(11, slot, set) == 0) cycle
                        if (t1 < bcap(5, slot, set) .or. t1 > bcap(6, slot, set) .or. t2 < bcap(7, slot, set) .or. t2 > bcap(8, &
                            & slot, set)) cycle
                        select case (id)
                        case (1)
                            i1 = bcap(1, slot, set); i2 = bcap(3, slot, set) + t1; i3 = bcap(4, slot, set) + t2
                            j1 = bcap(2, slot, set); j2 = i2; j3 = i3
                        case (2)
                            i1 = bcap(3, slot, set) + t1; i2 = bcap(1, slot, set); i3 = bcap(4, slot, set) + t2
                            j1 = i1; j2 = bcap(2, slot, set); j3 = i3
                        case (3)
                            i1 = bcap(3, slot, set) + t1; i2 = bcap(4, slot, set) + t2; i3 = bcap(1, slot, set)
                            j1 = i1; j2 = i2; j3 = bcap(2, slot, set)
                        end select
                        ! Reads stay inside the face gates: an unowned face's index can lie a block width outside this rank's
                        ! flux_rsx_vf bounds (hidden at np=1, where both faces are owned).
                        #:for SIDE, GATE, IX in [('lo', 9, 'i'), ('hi', 10, 'j')]
                            if (bcap(${GATE}$, slot, set) == 1) then
                                if (advective) then
                                    v = flux_rsx_vf(${IX}$1, ${IX}$2, ${IX}$3, eq)
                                else
                                    v = flux_src_rsx_vf(${IX}$1, ${IX}$2, ${IX}$3, eq)
                                end if
                                if (fine) then
                                    freg(id)%${SIDE}$(eq, t1, t2, slot) = merge(freg(id)%${SIDE}$(eq, t1, t2, slot), 0._wp, &
                                         & acc) + cf*v
                                else
                                    creg(id)%${SIDE}$(eq, t1, t2, slot) = merge(creg(id)%${SIDE}$(eq, t1, t2, slot), 0._wp, &
                                         & acc) + cf*v
                                end if
                            end if
                        #:endfor
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_capture_batch

    !> Capture the filled slots (s_amr_capture_slot) into freg (fine) or creg: advective, viscous (mom..E), and chemistry diffusion
    !! (species, plus the energy flux only when not viscous, since the viscous pass already captured it). Coarse and fine share this
    !! routine so the reflux subtracts matching fluxes. Set 0 is pushed here and cleared after.
    impure subroutine s_amr_capture_total(set, nb, id, fine, cf, acc)

        integer, intent(in)  :: set, nb, id
        logical, intent(in)  :: fine, acc
        real(wp), intent(in) :: cf
        integer              :: m1, m2

        if (.not. any(bcap(11,1:nb,set) == 1)) return
        if (set == 0) then
            $:GPU_UPDATE(device='[bcap(:, 1:nb, 0)]')
        end if
        m1 = bmax(1, set); m2 = bmax(2, set)
        call s_amr_capture_batch(set, nb, id, fine, .true., cf, acc, m1, m2, 1, sys_size)
        if (viscous) call s_amr_capture_batch(set, nb, id, fine, .false., cf, .true., m1, m2, eqn_idx%mom%beg, eqn_idx%E)
        if (chemistry .and. chem_params%diffusion) then
            call s_amr_capture_batch(set, nb, id, fine, .false., cf, .true., m1, m2, eqn_idx%species%beg, eqn_idx%species%end)
            if (.not. viscous) call s_amr_capture_batch(set, nb, id, fine, .false., cf, .true., m1, m2, eqn_idx%E, eqn_idx%E)
        end if
        if (set == 0) then
            bcap(11,:,0) = 0; bmax(:,0) = 0
        end if

    end subroutine s_amr_capture_total

    !> Capture c/f boundary-face fluxes for direction id from the finalized flux array, inside s_compute_rhs. The coarse call fills
    !! creg at block faces; the fine call (amr_in_fine_advance, globals swapped to the fine block) fills freg at the fine faces.
    impure subroutine s_amr_capture_boundary_flux(id, stage)

        integer, intent(in) :: id, stage
        integer             :: save_cur, kc, ibm, ko(3), ext(3), cflo(3), cfhi(3)
        logical             :: accum
        real(wp)            :: coef

        if (.not. amr) return
        call s_amr_reg_prepare()
        if (igr) return  ! IGR coupling is restriction-only: the fused flux kernels expose no face fluxes
        if (amr_in_fine_advance .and. .not. amr_rank_owns_block) return
        save_cur = amr_cur
        if (amr_in_fine_advance) then
            ! fine branch. Batched advance (amr_bat_n > 0): member ibm sits at offset ko in the slab and captures into its own slot.
            ! Level >= 2 and subcycled blocks accumulate the rk3_w-weighted (and dt-weighted) step integral for the once-per-step
            ! state reflux; otherwise each stage overwrites.
            if (amr_sub_on) then
                coef = rk3_w(stage)*dt; accum = .not. (amr_sub_first .and. stage == 1)
            else if (amr_block_level(amr_cur) >= 2) then
                coef = rk3_w(stage); accum = (stage > 1)
            else
                coef = 1._wp; accum = .false.
            end if
            do ibm = 1, max(1, amr_bat_n)
                call s_amr_bat_member(ibm, ko, ext)
                call s_amr_capture_slot(0, amr_reg_cur, id, ko, ko + ext, ko, ko + ext, .true., .true.)
            end do
            call s_amr_capture_total(0, amr_reg_n, id, .true., coef, accum)
            ! This block is the coarse side of its level+1 children: capture each child's creg at the child's footprint faces
            ! (s_amr_parent_foot, since amr_isect_*_all is empty for non-owned children), at the child's dense slot. Both faces are
            ! always owned.
            do ibm = 1, max(1, amr_bat_n)
                if (amr_block_level(save_cur) >= amr_max_level) exit  ! finest level has no children
                call s_amr_bat_member(ibm, ko, ext)
                do kc = 1, amr_num_blocks
                    if (.not. f_amr_is_child(kc, amr_cur)) cycle
                    call s_amr_parent_foot(kc, amr_cur, cflo, cfhi)
                    call s_amr_capture_slot(0, amr_reg_of(kc), id, cflo + ko, cfhi + ko, cflo + ko, cfhi + ko, .true., .true.)
                end do
            end do
            if (amr_bat_n > 0) call s_amr_select_slot(save_cur)
            ! children's coarse side: the parent's step integral
            call s_amr_capture_total(0, amr_reg_n, id, .false., rk3_w(stage)*merge(dt, 1._wp, amr_sub_on), stage > 1)
        else
            call s_amr_coarse_tables()
            ! subcycled: the coarse step's dt-weighted flux integral
            if (amr_sub_on) then
                call s_amr_capture_total(id, amr_reg_n, id, .false., rk3_w(stage)*dt, stage > 1)
            else
                call s_amr_capture_total(id, amr_reg_n, id, .false., 1._wp, .false.)
            end if
        end if

    end subroutine s_amr_capture_boundary_flux

    !> Build the coarse per-direction capture tables (bcap sets 1:3) and reflux apply tables (aint/areal) once per topology change
    !! and push them to the device. Each face is captured by the rank owning the coarse cells just outside it, over its owned
    !! transverse overlap [tlo:thi] in the block-relative frame. Called only on the coarse grid.
    impure subroutine s_amr_coarse_tables()

        integer :: k, d, save_cur, sreg, ta, tb, sidx(3), ext(3), tlo(3), thi(3)
        logical :: own_lo(3), own_hi(3)

        if (amr_ctab_epoch_built == amr_mesh_epoch .and. amr_ctab_nblk_built == amr_num_blocks) return
        save_cur = amr_cur
        bcap(11,:,1:3) = 0; bmax(:,1:3) = 0; aint(11,:,:) = 0; a_n = 0; a_m = 0
        do k = 1, amr_num_blocks
            ! a level >= 2 block's coarse side is its parent (captured in the fine branch)
            if (amr_block_level(k) >= 2) cycle
            if (.not. f_amr_block_near(k, amr_sidx, amr_ext)) cycle
            call s_amr_select_slot(k)
            call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi, tlo, thi)
            sreg = amr_reg_of(k)
            do d = 1, num_dims
                if (.not. (own_lo(d) .or. own_hi(d))) cycle
                call s_amr_capture_slot(d, sreg, d, amr_region_lo - sidx, amr_region_hi - sidx, tlo - sidx, thi - sidx, &
                                        & own_lo(d), own_hi(d))
                if (amr_block_level(k) /= 1) cycle
                ta = merge(2, 1, d == 1); tb = merge(2, 3, d == 3)
                aint(:,sreg, d) = [amr_region_lo(d) - 1 - sidx(d), amr_region_hi(d) + 1 - sidx(d), amr_region_lo(ta) - sidx(ta), &
                     & amr_region_lo(tb) - sidx(tb), tlo(ta) - amr_region_lo(ta), thi(ta) - amr_region_lo(ta), &
                     & tlo(tb) - amr_region_lo(tb), thi(tb) - amr_region_lo(tb), merge(1, 0, own_lo(d)), merge(1, 0, own_hi(d)), 1]
                areal(:,sreg, d) = 1._wp
                #:for D, DX in [(1, 'dx'), (2, 'dy'), (3, 'dz')]
                    if (d == ${D}$) then
                        if (own_lo(d)) areal(1, sreg, d) = ${DX}$(aint(1, sreg, d))
                        if (own_hi(d)) areal(2, sreg, d) = ${DX}$(aint(2, sreg, d))
                    end if
                #:endfor
                a_n(d) = a_n(d) + 1
                a_m(:,d) = max(a_m(:,d), aint(6:8:2,sreg, d) - aint(5:7:2,sreg, d))
            end do
        end do
        call s_amr_select_slot(save_cur)
        $:GPU_UPDATE(device='[bcap, aint, areal]')
        amr_ctab_epoch_built = amr_mesh_epoch
        amr_ctab_nblk_built = amr_num_blocks

    end subroutine s_amr_coarse_tables

    !> Correct the coarse rhs in the first cell outside each block face so the update sees the child-averaged fine flux: the low
    !! face is the outside cell's right face (rhs += (F_coarse - Fbar_fine)/dx), the high face its left face (rhs += (Fbar_fine -
    !! F_coarse)/dx). c1/c2 are relative 0-based coarse transverse indices. CCE 19.0.0: two or more of the per-direction loop nests
    !! in one procedure fail (ftn-7991 on offload, bad code on CPU), so each lives in its own procedure
    !! (s_amr_apply_reflux_d1/d2/d3); keep them separate.
    impure subroutine s_amr_apply_reflux(rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        if (.not. amr) return
        call s_amr_reg_prepare()
        if (igr) return

        #:for D in [1, 2, 3]
            if (amr_dim(${D}$)) call s_amr_apply_reflux_d${D}$(rhs_vf)
        #:endfor

    end subroutine s_amr_apply_reflux

    ! Batched over level-1 blocks, one kernel per face direction (tables from s_amr_coarse_tables). Same-side corrections of
    ! distinct blocks are disjoint, so the batched kernel is race-free. Per direction d the transverse dims are (ta, tb).
    #:for D, TA, TB, IDX in [(1, 2, 3, 'aint(1, k, 1), i2, i3'), (2, 1, 3, 'i2, aint(1, k, 2), i3'), (3, 1, 2, &
                              & 'i2, i3, aint(1, k, 3)')]
        !> Reflux correction for the faces normal to direction ${D}$. Called only when amr_dim(${D}$).
        impure subroutine s_amr_apply_reflux_d${D}$(rhs_vf)

            type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
            integer                                                :: eq, c1, c2, c1w, c2w, k, gmax1, gmax2
            integer                                                :: f10, f20, dd1, dd2, nch, rr, dd1_hi, dd2_hi, side
            integer                                                :: i2, i3
            logical                                                :: dta, dtb
            real(wp)                                               :: fb

            call s_amr_coarse_tables()
            if (a_n(${D}$) == 0) return
            rr = amr_ref_ratio
            dta = amr_dim(${TA}$); dtb = amr_dim(${TB}$)
            nch = rr**count([dta, dtb])
            dd1_hi = merge(rr - 1, 0, dta); dd2_hi = merge(rr - 1, 0, dtb)
            gmax1 = a_m(1, ${D}$); gmax2 = a_m(2, ${D}$)
            ! separate launches: one block's high outside cell can be another's low outside cell
            do side = 0, 1
                $:GPU_PARALLEL_LOOP(collapse=4, private='[c1, c2, f10, f20, dd1, dd2, fb, i2, i3]')
                do k = 1, amr_reg_n
                    do c2w = 0, gmax2
                        do c1w = 0, gmax1
                            do eq = 1, sys_size
                                if (aint(11, k, ${D}$) == 0 .or. aint(9 + side, k, ${D}$) == 0) cycle
                                c1 = aint(5, k, ${D}$) + c1w; c2 = aint(7, k, ${D}$) + c2w
                                if (c1 > aint(6, k, ${D}$) .or. c2 > aint(8, k, ${D}$)) cycle
                                i2 = aint(3, k, ${D}$) + c1; i3 = aint(4, k, ${D}$) + c2
                                f10 = 0; if (dta) f10 = rr*c1
                                f20 = 0; if (dtb) f20 = rr*c2
                                fb = 0._wp
                                do dd2 = 0, dd2_hi
                                    do dd1 = 0, dd1_hi
                                        if (side == 0) then
                                            fb = fb + freg(${D}$)%lo(eq, f10 + dd1, f20 + dd2, k)
                                        else
                                            fb = fb + freg(${D}$)%hi(eq, f10 + dd1, f20 + dd2, k)
                                        end if
                                    end do
                                end do
                                fb = fb/real(nch, wp)
                                if (side == 0) then
                                    rhs_vf(eq)%sf(${IDX}$) = rhs_vf(eq)%sf(${IDX}$) + (creg(${D}$)%lo(eq, c1, c2, &
                                           & k) - fb)/areal(1, k, ${D}$)
                                else
                                    rhs_vf(eq)%sf(${IDX.replace('aint(1,', 'aint(2,')}$) &
                                           & = rhs_vf(eq)%sf(${IDX.replace('aint(1,', 'aint(2,')}$) + (fb - creg(${D}$)%hi(eq, &
                                           & c1, c2, k))/areal(2, k, ${D}$)
                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do

        end subroutine s_amr_apply_reflux_d${D}$
    #:endfor

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

    !> Berger-Colella state reflux: q(outside) += w*dtl*(F_coarse - Fbar_fine)/m on the low face and the negative on the high face,
    !! per active dim, with F_coarse = creg and Fbar_fine = freg averaged over the rr**(ndim-1) covering fine faces. Framing is
    !! caller-passed: islot dense register slot; olo/ohi outside coarse-cell index below/above the face; glo/ghi creg-local loop
    !! range; woff transverse write origin; w_lo/w_hi per-face weight; mlo/mhi outside-cell width. A zero weight skips the write
    !! because the outside index may be out of bounds on an unowned face.
    impure subroutine s_amr_reflux_apply_faces(q, islot, rr, dtl, olo, ohi, glo, ghi, woff, w_lo, w_hi, mlo, mhi)

        type(scalar_field), dimension(sys_size), intent(inout) :: q
        integer, intent(in) :: islot, rr, olo(3), ohi(3), glo(3), ghi(3), woff(3)
        real(wp), intent(in) :: dtl, w_lo(3), w_hi(3), mlo(3), mhi(3)
        integer :: eq, g1, g2, f10, f20, dd1, dd2, nch, dd1_hi, dd2_hi, ol, oh, wa, wb, gla, gha, glb, ghb
        logical :: dta, dtb
        real(wp) :: fblo, fbhi, wl, wh, ml, mh

        ! Loop bounds are hoisted to scalars: array-element bounds would force the host arrays present on the device (ACC present
        ! error).

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
        @:DEALLOCATE(bcap, aint, areal)
        if (allocated(amr_reg_of)) deallocate (amr_reg_of)
        amr_reg_n = 0; amr_reg_cur = 0

    end subroutine s_finalize_amr_registers

end module m_amr_registers
