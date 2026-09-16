!>
!!@file
!!@brief Contains module m_amr_transfer

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Prolongation, restriction and Berger-Colella reflux, including the restrict/reflux/freg waves.
module m_amr_transfer

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_box, only: f_morton  ! shared 3D Morton key (single-sourced with m_sfc_partition)
    use m_global_parameters
    use m_constants, only: num_fluids_max, model_eqns_6eq, mapCells, K_ib, K_pc, BC_GHOST_EXTRAP
    use m_pressure_relaxation, only: s_pressure_relaxation_procedure
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_min, s_mpi_allreduce_integer_max, s_mpi_allreduce_sum, s_mpi_allreduce_min, &
        & s_mpi_allreduce_max, s_mpi_allreduce_integer_sum, s_mpi_sendrecv_variables_buffers, s_mpi_allreduce_array_max
    use m_rhs, only: s_compute_rhs, q_prim_qp
    use m_variables_conversion, only: s_convert_species_to_mixture_variables_kernel, s_compute_pressure, enforce_density_floor_vc
    use m_phase_change, only: s_infinite_relaxation_k, pc_iter_count
    use m_amr_registers, only: s_amr_zero_fine_registers, s_amr_reflux_apply_faces, s_amr_parent_foot, freg, creg, &
        & s_amr_reg_prepare, f_amr_face_is_seam
    use m_rank_timing, only: s_rank_time_tic, s_rank_time_toc
    use m_phase_timing
    use m_amr_xchg_audit  ! per-call-site accounting of every AMR p2p transfer (s_xa_rec + XA_* site ids)
    use m_ibm, only: s_ibm_alloc_fine, s_ibm_setup_fine, s_ibm_swap_to_fine, s_ibm_restore_from_fine, s_ibm_correct_state, &
        & s_ibm_load_fine_markers, s_update_mib, moving_immersed_boundary_flag, num_gps, ib_markers
    use m_hypoelastic, only: s_hypoelastic_update_fd_coeffs
    use m_weno, only: s_compute_weno_coefficients
    use m_active_box, only: ab_active
    use m_bubbles_EL, only: s_lag_cloud_bbox_local
    use m_igr, only: jac, jac_old
    use m_amr_state
    use m_amr_distribution
    use m_amr_store
    use m_amr_exchange
    use m_amr_frame

    implicit none

    private
    public :: s_amr_freg_wave, s_amr_p2p_reflux_faces, s_amr_reduce_xchg_flag, s_amr_reflux_faces_wave, s_amr_reflux_to_parent, &
        & s_amr_restrict_to_parent, s_amr_restrict_wave, s_interpolate_coarse_to_fine, s_populate_amr_fine, &
        & s_restrict_fine_to_coarse, s_set_amr_fine_geometry

contains

    !> True iff rank r is a reflux applier for the current block: it owns the coarse cell layer just outside some block face and its
    !! subdomain overlaps the block transversely. Mirrors s_amr_reflux_face_flags, but parameterized by r's subdomain from the
    !! computed decomposition (s_amr_rank_decomp, so the block owner can decide which ranks to send freg to, and each rank agrees on
    !! whether it receives). Uses amr_region_lo/hi (the current block, set on every rank by s_amr_select_slot). Deliberately no
    !! f_amr_face_is_seam clip (unlike the flags); the participation-map build (s_amr_reg_prepare, m_amr_registers) copies this
    !! unclipped formula for its clause (c), because both exchange paths gate their freg receives on it. Keep them in lockstep.
    pure logical function f_amr_reflux_participates(r) result(part)

        integer, intent(in) :: r
        integer             :: sidx(3), ext(3), d, t
        logical             :: tv(3), tvd

        call s_amr_rank_decomp(r, sidx, ext)
        tv(1) = amr_region_lo(1) <= sidx(1) + ext(1) .and. amr_region_hi(1) >= sidx(1)
        tv(2) = (n_glb == 0) .or. (amr_region_lo(2) <= sidx(2) + ext(2) .and. amr_region_hi(2) >= sidx(2))
        tv(3) = (p_glb == 0) .or. (amr_region_lo(3) <= sidx(3) + ext(3) .and. amr_region_hi(3) >= sidx(3))
        part = .false.
        do d = 1, num_dims
            tvd = .true.
            do t = 1, num_dims
                if (t /= d) tvd = tvd .and. tv(t)
            end do
            if (tvd .and. amr_region_lo(d) - 1 >= sidx(d) .and. amr_region_lo(d) - 1 <= sidx(d) + ext(d)) part = .true.
            if (tvd .and. amr_region_hi(d) + 1 >= sidx(d) .and. amr_region_hi(d) + 1 <= sidx(d) + ext(d)) part = .true.
        end do

    end function f_amr_reflux_participates

    !> Per-face refinement of f_amr_reflux_participates: the faces of the current block that rank r actually applies (it owns the
    !! outside coarse layer with transverse overlap, minus fine-fine seam faces), mirroring s_amr_reflux_face_flags term for term
    !! (same ownership formula, same f_amr_face_is_seam exclusion). The reflux-faces wave ships exactly these: sender and every
    !! receiver derive the identical set from replicated data, so a face a rank never applies never rides the wire.
    pure subroutine s_amr_reflux_faces_for(r, s_lo, s_hi)

        integer, intent(in)  :: r
        logical, intent(out) :: s_lo(3), s_hi(3)
        integer              :: sidx(3), ext(3), d, t
        logical              :: tv(3), tvd

        call s_amr_rank_decomp(r, sidx, ext)
        tv(1) = amr_region_lo(1) <= sidx(1) + ext(1) .and. amr_region_hi(1) >= sidx(1)
        tv(2) = (n_glb == 0) .or. (amr_region_lo(2) <= sidx(2) + ext(2) .and. amr_region_hi(2) >= sidx(2))
        tv(3) = (p_glb == 0) .or. (amr_region_lo(3) <= sidx(3) + ext(3) .and. amr_region_hi(3) >= sidx(3))
        s_lo = .false.; s_hi = .false.
        do d = 1, num_dims
            tvd = .true.
            do t = 1, num_dims
                if (t /= d) tvd = tvd .and. tv(t)
            end do
            s_lo(d) = tvd .and. amr_region_lo(d) - 1 >= sidx(d) .and. amr_region_lo(d) - 1 <= sidx(d) + ext(d) &
                 & .and. .not. f_amr_face_is_seam(d, -1)
            s_hi(d) = tvd .and. amr_region_hi(d) + 1 >= sidx(d) .and. amr_region_hi(d) + 1 <= sidx(d) + ext(d) &
                 & .and. .not. f_amr_face_is_seam(d, 1)
        end do

    end subroutine s_amr_reflux_faces_for

    !> Fine-level distribution: deliver the current block's fine flux registers freg (captured by the owner during the fine advance)
    !! to exactly the coarse-outside-owners that apply the reflux, point-to-point. The owner sends its whole freg slot
    !! (block-relative; each applier reads its own transverse slice) to every participant; non-owner participants receive it.
    !! Device-resident: owner stages its slot to host, receivers push it back. No-op without MPI/at np=1.
    impure subroutine s_amr_p2p_reflux_faces()

#ifdef MFC_MPI
        integer              :: owner, r, ierr, nreq, cnt, idx, ncand
        integer              :: cand(num_procs), glo(3), ghi(3)
        integer, allocatable :: reqs(:)

        if (.not. amr) return
        if (num_procs == 1) return
        owner = amr_block_owner(amr_cur)
        if (proc_rank == owner) then
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    $:GPU_UPDATE(host='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur), freg(' + str(D) &
                                 & + ')%hi(:, :, :, amr_reg_cur)]')
                end if
            #:endfor
            ! participating ranks by O(overlap) inversion (region grown by 1) filtered by the participation predicate, rather
            ! than an O(P) rank scan. Ascending (owner-excluded at use), so the ISENDs match the receivers' order.
            glo = 0; ghi = 0
            glo(1) = amr_region_lo(1) - 1; ghi(1) = amr_region_hi(1) + 1
            if (n_glb > 0) then; glo(2) = amr_region_lo(2) - 1; ghi(2) = amr_region_hi(2) + 1; end if
            if (p_glb > 0) then; glo(3) = amr_region_lo(3) - 1; ghi(3) = amr_region_hi(3) + 1; end if
            call s_amr_ranks_overlapping(glo, ghi, cand, ncand)
            nreq = 0
            do idx = 1, ncand
                r = cand(idx)
                if (r /= owner .and. f_amr_reflux_participates(r)) nreq = nreq + 1
            end do
            if (nreq > 0) then
                allocate (reqs(2*num_dims*nreq))
                nreq = 0
                do idx = 1, ncand
                    r = cand(idx)
                    if (r == owner .or. .not. f_amr_reflux_participates(r)) cycle
                    #:for D in [1, 2, 3]
                        if (${D}$ <= num_dims) then
                            ! not slot amr_reg_cur: it is 0 on the L0-tiles path
                            cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                            nreq = nreq + 1
                            call s_xa_rec(XA_F5_FACE_SND, 1, cnt, ${2*D}$)
                            call MPI_ISEND(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, r, ${2*D}$, MPI_COMM_WORLD, reqs(nreq), &
                                           & ierr)
                            nreq = nreq + 1
                            call s_xa_rec(XA_F5_FACE_SND, 1, cnt, ${2*D + 1}$)
                            call MPI_ISEND(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, r, ${2*D + 1}$, MPI_COMM_WORLD, &
                                           & reqs(nreq), ierr)
                        end if
                    #:endfor
                end do
                call s_phase_tic(PH_RFWAIT)
                call s_wait_tic()
                call MPI_WAITALL(nreq, reqs, MPI_STATUSES_IGNORE, ierr)
                call s_wait_toc(WT_REFLUX)
                call s_phase_toc(PH_RFWAIT)
                deallocate (reqs)
            end if
        else if (f_amr_reflux_participates(proc_rank)) then
            ! Post all 2*num_dims receives, then one wait, so the faces are not serialised against the owner's send
            ! order. The slice (:,:,:,amr_reg_cur) is contiguous (the dense register slot is the last dimension) and the
            ! owner side ISENDs the identical shape, so there is no temporary-buffer hazard.
            call s_phase_tic(PH_RFRECV)
            allocate (reqs(2*num_dims))
            nreq = 0
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    ! not slot amr_reg_cur: it is 0 on the L0-tiles path
                    cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                    nreq = nreq + 1
                    call s_xa_rec(XA_F5_FACE_RCV, 2, cnt, ${2*D}$)
                    call MPI_IRECV(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, owner, ${2*D}$, MPI_COMM_WORLD, reqs(nreq), ierr)
                    nreq = nreq + 1
                    call s_xa_rec(XA_F5_FACE_RCV, 2, cnt, ${2*D + 1}$)
                    call MPI_IRECV(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, owner, ${2*D + 1}$, MPI_COMM_WORLD, reqs(nreq), &
                                   & ierr)
                end if
            #:endfor
            call s_wait_tic()
            call MPI_WAITALL(nreq, reqs, MPI_STATUSES_IGNORE, ierr)
            call s_wait_toc(WT_REFLUX)
            deallocate (reqs)
            ! Device update only after the wait: the buffers hold nothing valid until then.
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur), freg(' + str(D) &
                                 & + ')%hi(:, :, :, amr_reg_cur)]')
                end if
            #:endfor
            call s_phase_toc(PH_RFRECV)
        end if
#endif

    end subroutine s_amr_p2p_reflux_faces

    !> The per-stage level-1 reflux-face exchange as one wave. Every rank walks the level-1 slots ascending: all receives post first
    !! (zero-copy, directly into the freg host mirrors; each box owns a register slot, so no pool is needed), then the owners stage
    !! + multicast, then one waitall, then the receivers push to device. Every message carries its own keyed tag (band 0, per-pair
    !! seq in the (ascending box, ascending dim, lo-then-hi) plan order both sides derive from the same predicates), so matching
    !! does not depend on posting order or MPI non-overtaking. Under MFC_DEBUG the identity headers travel as separate 8-word
    !! companion messages, one per (box, peer) group ahead of its payloads (a prefix cannot ride a zero-copy payload); they are
    !! never recorded in [amr-xa], so the family words stay exactly comparable. The register arrays are sized up front: the apply
    !! can reallocate them, so nothing may post into freg before s_amr_reg_prepare.
    impure subroutine s_amr_reflux_faces_wave()

#ifdef MFC_MPI
        use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        integer  :: k, r, ierr, nreq, cnt, idx, ncand, tq, nhr, nhs, j, kk2, sq
        integer  :: cand(num_procs), glo(3), ghi(3)
        logical  :: s_lo(3), s_hi(3), u_lo(3), u_hi(3)
        logical  :: cl(3, num_procs), ch(3, num_procs)
        real(wp) :: nanv

        if (num_procs == 1) return
        call s_amr_reg_prepare()
        call s_amr_m1_wave_open(0)
        nanv = ieee_value(0._wp, ieee_quiet_nan)
        nreq = 0; nhr = 0; nhs = 0
        ! participates => the raw region +/-1 touches my interior slab => the region +/-amr_cpat_mar (>= 2) intersects my
        ! coarse range (a superset of the slab) => the block is in amr_l1p. The exact predicates below keep the survivor set and
        ! its ascending order identical to a full block scan.
        call s_amr_refresh_lists()
        if (XA_NH > 0) then
            ! A posted request reads its header slot until the wave drains, so both header pools are sized once here: growing
            ! them inside the post loops reallocates under in-flight sends and receives (the receiver then checks a zeroed header).
            call s_amr_refresh_my_blocks()
            call s_amr_fw_szr(amr_fw_rq, XA_NH*max(amr_n_l1p, 1), amr_fw_dev)
            call s_amr_fw_szr(amr_fw_sq, XA_NH*max(amr_n_my*num_procs, 1), amr_fw_dev)
        end if
        do kk2 = 1, amr_n_l1p
            k = amr_l1p_blk(kk2)
            call s_amr_select_slot(k)
            if (amr_block_owner(k) == proc_rank) cycle
            if (.not. f_amr_reflux_participates(proc_rank)) cycle
            ! face-selective multicast: receive exactly the faces this rank applies (s_amr_reflux_faces_for mirrors the
            ! apply's own_lo/own_hi + seam gates); the owner derives the same set per participant, so the pairing is
            ! exact with no metadata exchange. Debug arm: unreceived faces are NaN-flooded so any hidden reader aborts.
            call s_amr_reflux_faces_for(proc_rank, s_lo, s_hi)
            ! record the block: the apply pass below iterates this list rather than rescanning every block to re-derive the
            ! same order. Built unconditionally (not only under the XA_NH > 0 audit).
            nhr = nhr + 1
            call s_amr_fw_szi(amr_fw_rblk, nhr)
            amr_fw_rblk(nhr) = k
            if (XA_NH > 0) then
                @:ASSERT(size(amr_fw_rq) >= XA_NH*nhr, "amr_fw_rq header pool sized below the wave's receive count")
                nreq = nreq + 1
                call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                amr_fw_reqw(nreq) = XA_NH
                tq = f_amr_m1_tag(0, f_amr_m1_seq(amr_block_owner(k), 2))
                call MPI_IRECV(amr_fw_rq(XA_NH*(nhr - 1) + 1), XA_NH, mpi_p, amr_block_owner(k), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    ! not slot amr_reg_cur: it is 0 on the L0-tiles path
                    cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                    if (s_lo(${D}$)) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = cnt
                        sq = f_amr_m1_seq(amr_block_owner(k), 2); tq = f_amr_m1_tag(0, sq)
                        call s_xa_rec(XA_F5W_FACE_RCV, 2, cnt, tq, peer=amr_block_owner(k), key=k*8 + ${D}$*2, seq=sq)
                        call MPI_IRECV(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, amr_block_owner(k), tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
#ifdef MFC_DEBUG
                        ! amr_reg_cur is 0 when this rank holds no register for the block (s_amr_select_slot's unmapped
                        ! sentinel). There is then no buffer to poison, and nothing that could read one. Only this debug
                        ! branch can meet that case: the receives above are posted under s_lo/s_hi, which are false for a
                        ! block this rank has no register for. Indexing freg with 0 here is an out-of-bounds write.
                    else if (amr_reg_cur > 0) then
                        freg(${D}$)%lo(:,:,:,amr_reg_cur) = nanv
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
#endif
                    end if
                    if (s_hi(${D}$)) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = cnt
                        sq = f_amr_m1_seq(amr_block_owner(k), 2); tq = f_amr_m1_tag(0, sq)
                        call s_xa_rec(XA_F5W_FACE_RCV, 2, cnt, tq, peer=amr_block_owner(k), key=k*8 + ${D}$*2 + 1, seq=sq)
                        call MPI_IRECV(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, amr_block_owner(k), tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
#ifdef MFC_DEBUG
                    else if (amr_reg_cur > 0) then
                        freg(${D}$)%hi(:,:,:,amr_reg_cur) = nanv
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
#endif
                    end if
                end if
            #:endfor
        end do
        call s_amr_refresh_my_blocks()
        do kk2 = 1, amr_n_my  ! owned list; level filter kept
            k = amr_my_blk(kk2)
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            if (amr_block_owner(k) /= proc_rank) cycle  ! belt-and-braces
            glo = 0; ghi = 0
            glo(1) = amr_region_lo(1) - 1; ghi(1) = amr_region_hi(1) + 1
            if (n_glb > 0) then; glo(2) = amr_region_lo(2) - 1; ghi(2) = amr_region_hi(2) + 1; end if
            if (p_glb > 0) then; glo(3) = amr_region_lo(3) - 1; ghi(3) = amr_region_hi(3) + 1; end if
            call s_amr_ranks_overlapping(glo, ghi, cand, ncand)
            ! face-selective multicast: each participant's ship set is its apply set (s_amr_reflux_faces_for), derived
            ! here per candidate; the device->host pull covers only the union of shipped faces.
            u_lo = .false.; u_hi = .false.
            do idx = 1, ncand
                r = cand(idx)
                cl(:,idx) = .false.; ch(:,idx) = .false.
                if (r == proc_rank .or. .not. f_amr_reflux_participates(r)) cycle
                call s_amr_reflux_faces_for(r, s_lo, s_hi)
                cl(:,idx) = s_lo; ch(:,idx) = s_hi
                u_lo = u_lo .or. s_lo; u_hi = u_hi .or. s_hi
            end do
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    if (u_lo(${D}$)) then
                        $:GPU_UPDATE(host='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
                    end if
                    if (u_hi(${D}$)) then
                        $:GPU_UPDATE(host='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
                    end if
                end if
            #:endfor
            do idx = 1, ncand
                r = cand(idx)
                if (r == proc_rank .or. .not. f_amr_reflux_participates(r)) cycle
                if (XA_NH > 0) then
                    nhs = nhs + 1
                    @:ASSERT(size(amr_fw_sq) >= XA_NH*nhs, "amr_fw_sq header pool sized below the wave's send count")
                    call s_xa_hdr_pack(amr_fw_sq(XA_NH*(nhs - 1) + 1:XA_NH*nhs), XA_F5W_FACE_SND, k, [0, 0, 0], [0, 0, 0])
                    nreq = nreq + 1
                    call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                    amr_fw_reqw(nreq) = -1
                    tq = f_amr_m1_tag(0, f_amr_m1_seq(r, 1))
                    call MPI_ISEND(amr_fw_sq(XA_NH*(nhs - 1) + 1), XA_NH, mpi_p, r, tq, MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                end if
                #:for D in [1, 2, 3]
                    if (${D}$ <= num_dims) then
                        ! not slot amr_reg_cur: it is 0 on the L0-tiles path
                        cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                        if (cl(${D}$, idx)) then
                            nreq = nreq + 1
                            call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                            amr_fw_reqw(nreq) = -1
                            sq = f_amr_m1_seq(r, 1); tq = f_amr_m1_tag(0, sq)
                            call s_xa_rec(XA_F5W_FACE_SND, 1, cnt, tq, peer=r, key=k*8 + ${D}$*2, seq=sq)
                            call MPI_ISEND(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, r, tq, MPI_COMM_WORLD, &
                                           & amr_fw_req(nreq), ierr)
                        end if
                        if (ch(${D}$, idx)) then
                            nreq = nreq + 1
                            call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                            amr_fw_reqw(nreq) = -1
                            sq = f_amr_m1_seq(r, 1); tq = f_amr_m1_tag(0, sq)
                            call s_xa_rec(XA_F5W_FACE_SND, 1, cnt, tq, peer=r, key=k*8 + ${D}$*2 + 1, seq=sq)
                            call MPI_ISEND(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, r, tq, MPI_COMM_WORLD, &
                                           & amr_fw_req(nreq), ierr)
                        end if
                    end if
                #:endfor
            end do
        end do
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call s_phase_tic(PH_RFWAIT)
                call s_wait_tic()
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                call s_wait_toc(WT_REFLUX)
                call s_phase_toc(PH_RFWAIT)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "reflux-faces wave: received message length differs from the plan")
                end do
            end block
#else
            call s_phase_tic(PH_RFWAIT)
            call s_wait_tic()
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
            call s_wait_toc(WT_REFLUX)
            call s_phase_toc(PH_RFWAIT)
#endif
        end if
        call s_phase_tic(PH_RFRECV)
        ! the post pass recorded exactly the blocks this rank receives, in this order, so iterate that list
        do j = 1, nhr
            k = amr_fw_rblk(j)
            call s_amr_select_slot(k)
            if (XA_NH > 0) then
                call s_xa_hdr_check(amr_fw_rq(XA_NH*(j - 1) + 1:XA_NH*j), XA_F5W_FACE_SND, k, [0, 0, 0], [0, 0, 0])
            end if
            ! push only the received faces; an unreceived face keeps its device content (never applied here, and
            ! NaN-poisoned in debug, so any hidden reader aborts)
            call s_amr_reflux_faces_for(proc_rank, s_lo, s_hi)
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    if (s_lo(${D}$)) then
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
                    end if
                    if (s_hi(${D}$)) then
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
                    end if
                end if
            #:endfor
        end do
        call s_phase_toc(PH_RFRECV)
#endif

    end subroutine s_amr_reflux_faces_wave

    !> The split-ownership level>=2 freg exchange as one wave, run once before the reflux fold (the registers are final after the
    !! advance, and the applies keep their per-box reverse-order position). Same zero-copy, companion-header design as the faces
    !! wave; keyed tags on band 1 (the faces wave is band 0) keep the two disjoint. The subcycle path keeps its per-box exchange
    !! inside s_amr_reflux_to_parent (do_xchg).
    impure subroutine s_amr_freg_wave()

#ifdef MFC_MPI
        use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        integer  :: k, ierr, nreq, cnt, pblk, cowner, powner, tq, nhr, nhs, j, kk2, sq
        real(wp) :: w_lo(3), w_hi(3), nanv

        if (num_procs == 1) return
        call s_amr_reg_prepare()
        call s_amr_refresh_lists()
        call s_amr_m1_wave_open(1)
        nanv = ieee_value(0._wp, ieee_quiet_nan)
        nreq = 0; nhr = 0; nhs = 0
        if (XA_NH > 0) then
            ! header pools sized once before any post (see the faces wave)
            call s_amr_refresh_my_blocks()
            call s_amr_fw_szr(amr_fw_rq, XA_NH*max(amr_n_fch, 1), amr_fw_dev)
            call s_amr_fw_szr(amr_fw_sq, XA_NH*max(amr_n_my*num_procs, 1), amr_fw_dev)
        end if
        ! amr_fch_blk is this loop's survivor set (level >= 2, my parent, foreign child), ascending; the exact tests stay as
        ! belt-and-braces
        do kk2 = 1, amr_n_fch
            k = amr_fch_blk(kk2)
            call s_amr_select_slot(k)
            pblk = amr_parent_blk(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. powner /= proc_rank) cycle
            ! seam clip: a face weighted 0 by the sibling-seam rule is never consumed by the parent-side reflux apply
            ! (s_amr_reflux_to_parent multiplies it away), so it never ships; both sides derive the identical skip from
            ! s_amr_sibling_face_weights on replicated metadata. Debug arm: skipped-face mirrors are NaN-flooded so any
            ! other consumer of an unshipped face aborts within the step.
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            ! same as the faces wave: record the block so the apply pass need not rescan the block list
            nhr = nhr + 1
            call s_amr_fw_szi(amr_fw_rblk, nhr)
            amr_fw_rblk(nhr) = k
            if (XA_NH > 0) then
                @:ASSERT(size(amr_fw_rq) >= XA_NH*nhr, "amr_fw_rq header pool sized below the wave's receive count")
                nreq = nreq + 1
                call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                amr_fw_reqw(nreq) = XA_NH
                tq = f_amr_m1_tag(1, f_amr_m1_seq(cowner, 2))
                call MPI_IRECV(amr_fw_rq(XA_NH*(nhr - 1) + 1), XA_NH, mpi_p, cowner, tq, MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
            end if
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    ! not slot amr_reg_cur: it is 0 on the L0-tiles path
                    cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                    if (w_lo(${D}$) > 0._wp) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = cnt
                        sq = f_amr_m1_seq(cowner, 2); tq = f_amr_m1_tag(1, sq)
                        call s_xa_rec(XA_F5W_FREG_RCV, 2, cnt, tq, peer=cowner, key=k*8 + ${D}$*2 + 0, seq=sq)
                        call MPI_IRECV(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, cowner, tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
#ifdef MFC_DEBUG
                    else if (amr_reg_cur > 0) then
                        freg(${D}$)%lo(:,:,:,amr_reg_cur) = nanv
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
#endif
                    end if
                    if (w_hi(${D}$) > 0._wp) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = cnt
                        sq = f_amr_m1_seq(cowner, 2); tq = f_amr_m1_tag(1, sq)
                        call s_xa_rec(XA_F5W_FREG_RCV, 2, cnt, tq, peer=cowner, key=k*8 + ${D}$*2 + 1, seq=sq)
                        call MPI_IRECV(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, cowner, tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
#ifdef MFC_DEBUG
                    else if (amr_reg_cur > 0) then
                        freg(${D}$)%hi(:,:,:,amr_reg_cur) = nanv
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
#endif
                    end if
                end if
            #:endfor
        end do
        call s_amr_refresh_my_blocks()
        do kk2 = 1, amr_n_my  ! owned list; level filter kept
            k = amr_my_blk(kk2)
            if (amr_block_level(k) < 2) cycle
            call s_amr_select_slot(k)
            pblk = amr_parent_blk(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. cowner /= proc_rank) cycle
            ! seam clip, send side: the identical weight derivation as the recv walk (replicated metadata), so the
            ! posted sends pair the posted recvs exactly. Skipped faces also skip their device->host pulls.
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    if (w_lo(${D}$) > 0._wp) then
                        $:GPU_UPDATE(host='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
                    end if
                    if (w_hi(${D}$) > 0._wp) then
                        $:GPU_UPDATE(host='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
                    end if
                end if
            #:endfor
            if (XA_NH > 0) then
                nhs = nhs + 1
                @:ASSERT(size(amr_fw_sq) >= XA_NH*nhs, "amr_fw_sq header pool sized below the wave's send count")
                call s_xa_hdr_pack(amr_fw_sq(XA_NH*(nhs - 1) + 1:XA_NH*nhs), XA_F5W_FREG_SND, k, [0, 0, 0], [0, 0, 0])
                nreq = nreq + 1
                call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                amr_fw_reqw(nreq) = -1
                tq = f_amr_m1_tag(1, f_amr_m1_seq(powner, 1))
                call MPI_ISEND(amr_fw_sq(XA_NH*(nhs - 1) + 1), XA_NH, mpi_p, powner, tq, MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
            end if
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    ! not slot amr_reg_cur: it is 0 on the L0-tiles path
                    cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                    if (w_lo(${D}$) > 0._wp) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = -1
                        sq = f_amr_m1_seq(powner, 1); tq = f_amr_m1_tag(1, sq)
                        call s_xa_rec(XA_F5W_FREG_SND, 1, cnt, tq, peer=powner, key=k*8 + ${D}$*2 + 0, seq=sq)
                        call MPI_ISEND(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, powner, tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
                    end if
                    if (w_hi(${D}$) > 0._wp) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = -1
                        sq = f_amr_m1_seq(powner, 1); tq = f_amr_m1_tag(1, sq)
                        call s_xa_rec(XA_F5W_FREG_SND, 1, cnt, tq, peer=powner, key=k*8 + ${D}$*2 + 1, seq=sq)
                        call MPI_ISEND(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, powner, tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
                    end if
                end if
            #:endfor
        end do
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call s_wait_tic()
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                call s_wait_toc(WT_RESTR)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "freg wave: received message length differs from the plan")
                end do
            end block
#else
            call s_wait_tic()
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
            call s_wait_toc(WT_RESTR)
#endif
        end if
        ! iterate the recorded receive list, as the faces wave does
        do j = 1, nhr
            k = amr_fw_rblk(j)
            call s_amr_select_slot(k)
            pblk = amr_parent_blk(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (XA_NH > 0) then
                call s_xa_hdr_check(amr_fw_rq(XA_NH*(j - 1) + 1:XA_NH*j), XA_F5W_FREG_SND, k, [0, 0, 0], [0, 0, 0])
            end if
            ! push only the faces that shipped; a skipped face keeps its device content (dead under weight 0, and
            ! NaN-poisoned in debug, so any other reader aborts)
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    if (w_lo(${D}$) > 0._wp) then
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
                    end if
                    if (w_hi(${D}$) > 0._wp) then
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
                    end if
                end if
            #:endfor
        end do
#endif

    end subroutine s_amr_freg_wave

    !> Set the fine level's geometry (region, intersection, extents, bounds, coordinates) for the box lo:hi. Arrays are preallocated
    !! at max size; this only updates metadata and refills coords. Collective: all ranks must call together (init and regrid do); it
    !! also refreshes the allreduced amr_xchg_coarse_ghosts flag for the new box. Invariant: a level-l block's fine extent is
    !! amr_ref_ratio**l * (coarse-region width) - 1, not amr_ref_ratio*width. (amr_ref_ratio*width holds only for the level-1
    !! initial block; nested boxes compound by amr_ref_ratio per level.) Every fine-extent computation (here, the restart-reader
    !! check, the load-weight, the fmul) uses amr_ref_ratio**level; assuming amr_ref_ratio*width rejects level>=2 blocks as corrupt.
    impure subroutine s_set_amr_fine_geometry(lo, hi)

        integer, intent(in) :: lo(3), hi(3)
        integer             :: sidx(3), ext(3), nmar, bad_loc, pblk

        amr_slots(amr_cur)%region%lo = lo; amr_slots(amr_cur)%region%hi = hi
        amr_region_lo = lo; amr_region_hi = hi  ! global mirror for m_amr_registers (no use-cycle)
        amr_region_lo_all(:,amr_cur) = lo; amr_region_hi_all(:,amr_cur) = hi

        ! Fine-level distribution: a block is owned whole by amr_block_owner(k). The owner holds fine cells for the entire block;
        ! every other rank holds none. amr_isect_lo/hi records the block's coarse footprint (= the whole block on the owner); it
        ! drives the coarse<->fine gather/scatter (which coarse cells the owner pulls in / pushes back). At np=1 the owner is rank 0
        ! and the footprint is the whole domain-resident block.
        amr_rank_owns_block = (amr_block_owner(amr_cur) == proc_rank)
        pblk = 0
        if (amr_rank_owns_block) then
            amr_isect_lo = lo; amr_isect_hi = hi
            if (amr_block_level(amr_cur) >= 2) then
                ! multi-level: express the coarse footprint in the parent block's fine-cell frame (a level-l block's coarse side
                ! is level l-1). parent-fine index of L0 cell c is rr*(c - R1.lo) where rr is the parent's amr_ref_ratio; the
                ! block spans rr fine cells per parent-covered L0 cell. m below then gets amr_ref_ratio*(footprint) cells, as for
                ! a level-1 block over L0. amr_cg / the prolong read this frame. rr is the global amr_ref_ratio, not
                ! amr_slots(pblk)%amr_ref_ratio: that field is written by s_amr_alloc_slot, which a rank owning this block but not
                ! its parent never calls for pblk, so it would read undefined. The two agree wherever both are defined (only an L0
                ! tile carries a per-slot ratio of 1, and a level>=2 block's parent is never an L0 tile). This is the same
                ! footprint s_amr_parent_foot derives from replicated metadata.
                pblk = f_amr_parent_block(amr_cur)
                call s_amr_parent_foot(amr_cur, pblk, amr_isect_lo, amr_isect_hi)
            end if
        else
            amr_isect_lo = 1; amr_isect_hi = 0  ! empty footprint
            if (n_glb > 0) then; amr_isect_lo(2) = 1; amr_isect_hi(2) = 0; end if
            if (p_glb > 0) then; amr_isect_lo(3) = 1; amr_isect_hi(3) = 0; end if
        end if
        amr_isect_lo_all(:,amr_cur) = amr_isect_lo; amr_isect_hi_all(:,amr_cur) = amr_isect_hi
        amr_owns_all(amr_cur) = amr_rank_owns_block
        ! fine extents cover the whole block on the owner; -1 (empty) on non-owners
        amr_slots(amr_cur)%m = amr_ref_ratio*max(amr_isect_hi(1) - amr_isect_lo(1) + 1, 0) - 1
        amr_slots(amr_cur)%n = 0; amr_slots(amr_cur)%p = 0
        if (n_glb > 0) amr_slots(amr_cur)%n = amr_ref_ratio*max(amr_isect_hi(2) - amr_isect_lo(2) + 1, 0) - 1
        if (p_glb > 0) amr_slots(amr_cur)%p = amr_ref_ratio*max(amr_isect_hi(3) - amr_isect_lo(3) + 1, 0) - 1
        amr_slots(amr_cur)%idwbuff(1)%beg = -buff_size; amr_slots(amr_cur)%idwbuff(1)%end = amr_slots(amr_cur)%m + buff_size
        amr_slots(amr_cur)%idwbuff(2)%beg = 0; amr_slots(amr_cur)%idwbuff(2)%end = 0
        amr_slots(amr_cur)%idwbuff(3)%beg = 0; amr_slots(amr_cur)%idwbuff(3)%end = 0
        if (n_glb > 0) then
            amr_slots(amr_cur)%idwbuff(2)%beg = -buff_size; amr_slots(amr_cur)%idwbuff(2)%end = amr_slots(amr_cur)%n + buff_size
        end if
        if (p_glb > 0) then
            amr_slots(amr_cur)%idwbuff(3)%beg = -buff_size; amr_slots(amr_cur)%idwbuff(3)%end = amr_slots(amr_cur)%p + buff_size
        end if
        ! coord building only on ranks with fine cells (others never read their coord arrays)
        if (amr_rank_owns_block) then
            ! Every level builds the same way: replay the ancestor chain from the global L0 boundaries. The owner may hold no part
            ! of the coarse slice it refines, and (level>=2) may not own the parent at all, so neither the local coarse coords nor
            ! the parent's slot can be read here. At level 1 the chain is one step.
            call s_amr_build_block_coords(amr_cur, amr_gxcb, amr_slots(amr_cur)%x_cb, amr_slots(amr_cur)%x_cc, &
                                          & amr_slots(amr_cur)%dx, 1)
            if (n_glb > 0) call s_amr_build_block_coords(amr_cur, amr_gycb, amr_slots(amr_cur)%y_cb, amr_slots(amr_cur)%y_cc, &
                & amr_slots(amr_cur)%dy, 2)
            if (p_glb > 0) call s_amr_build_block_coords(amr_cur, amr_gzcb, amr_slots(amr_cur)%z_cb, amr_slots(amr_cur)%z_cc, &
                & amr_slots(amr_cur)%dz, 3)
        end if

        ! Fine ghost prolongation reads up to nmar coarse cells past each face of the intersection; if that stencil leaves any
        ! rank's interior (block near/at/across a rank boundary), the coarse cons ghosts it reads must be halo-exchanged before
        ! every fill (the solver populates only prim ghosts). All ranks agree on the flag, so the pairwise exchanges are called
        ! consistently.
        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        nmar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        bad_loc = 0
        if (amr_rank_owns_block) then
            if (amr_isect_lo(1) - sidx(1) < nmar .or. sidx(1) + ext(1) - amr_isect_hi(1) < nmar) bad_loc = 1
            if (n_glb > 0 .and. (amr_isect_lo(2) - sidx(2) < nmar .or. sidx(2) + ext(2) - amr_isect_hi(2) < nmar)) bad_loc = 1
            if (p_glb > 0 .and. (amr_isect_lo(3) - sidx(3) < nmar .or. sidx(3) + ext(3) - amr_isect_hi(3) < nmar)) bad_loc = 1
        end if
        ! Accumulate, do not reduce: the caller closes the scan with s_amr_reduce_xchg_flag. Every caller loops over blocks and
        ! wants "does any block need the exchange" (the OR over blocks, not the last block's answer), in one collective rather
        ! than one per block.
        amr_xchg_bad = max(amr_xchg_bad, bad_loc)

    end subroutine s_set_amr_fine_geometry

    !> Close a geometry scan: one allreduce of the accumulated flag, then reset so the next scan starts clean. Must be called after
    !! every s_set_amr_fine_geometry loop (or single call); the flag it sets is read by the fine advance
    !! (s_amr_exchange_coarse_cons_halo).
    impure subroutine s_amr_reduce_xchg_flag()

        integer :: bad_glb

#ifdef MFC_MPI
        integer  :: ierr
        real(wp) :: t0, t1, tmin, tmax

        t0 = f_amr_wtime()
#endif
        call s_mpi_allreduce_integer_max(amr_xchg_bad, bad_glb)
#ifdef MFC_MPI
        t1 = f_amr_wtime()
        if (rank_time_wrt) then
            call MPI_ALLREDUCE(t0, tmin, 1, mpi_p, MPI_MIN, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(t0, tmax, 1, mpi_p, MPI_MAX, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) print '(A,ES10.3,A,ES10.3)', '[amr-rb] xchg_skew ', tmax - tmin, ' xchg_coll ', t1 - t0
        end if
#endif
        amr_xchg_coarse_ghosts = bad_glb == 1
        amr_xchg_bad = 0

    end subroutine s_amr_reduce_xchg_flag

    !> Conservative-linear prolongation for a single variable pair. Reads coarse interior/ghost from qc; writes fine interior to qf.
    !! Minmod-limited slopes.
    impure subroutine s_prolong_one_var(qc, loc, ivar, pos, inject)

        type(scalar_field), intent(in) :: qc
        integer, intent(in)            :: loc, ivar  !< flat-store slot and variable of the fine target
        logical, optional, intent(in)  :: pos        !< floor the child at bub_pos_frac*u0 (bubble radius-moment realizability)
        logical, optional, intent(in)  :: inject     !< piecewise-constant (child = u0): QBMM moment realizability preservation
        integer                        :: fi, fj, fk, ci, cj, ck, ox, oy, oz, rrat, mm, nn, pp, il1, il2, il3
        real(wp)                       :: u0, sx, sy, sz, xix, xiy, xiz, child, bpf
        logical                        :: floor_pos, pw_const, d2, d3

        floor_pos = .false.; if (present(pos)) floor_pos = pos
        pw_const = .false.; if (present(inject)) pw_const = inject

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): amr_isect_lo is global and equals
        ! region_lo on the owner, so amr_isect_lo + f/rr - amr_cpat_off = nmar + f/rr is the patch-local coarse index.
        ! Device kernel: reads the patch's device mirror (pushed once per prolong dispatch by s_interpolate_coarse_to_fine) and
        ! writes the fine slot in place. CPU builds compile this to the identical plain loop.
        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        bpf = bub_pos_frac
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, child]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    u0 = real(qc%sf(ci, cj, ck), wp)
                    sx = minmod(real(qc%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc%sf(ci - 1, cj, ck), wp))
                    sy = 0._wp
                    if (d2) sy = minmod(real(qc%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc%sf(ci, cj - 1, ck), wp))
                    sz = 0._wp
                    if (d3) sz = minmod(real(qc%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc%sf(ci, cj, ck - 1), wp))
                    if (pw_const) then
                        sx = 0._wp; sy = 0._wp; sz = 0._wp
                    end if
                    child = u0 + sx*xix + sy*xiy + sz*xiz
                    if (floor_pos) child = max(child, bpf*u0)
                    amr_cons_st(fi, fj, fk, ivar, loc) = child
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_one_var

    !> Conservative-linear prolongation: fill amr_fine interior from coarse (level-0), minmod-limited. Symmetric child offsets
    !! (+/-1/4 of a coarse cell) => the amr_ref_ratio^d children average to the coarse value. Multi-fluid volume fractions take the
    !! sum-preserving closure path instead (single-fluid runs never branch, so their prolongation is untouched). Twin
    !! s_amr_prolong_pbmv (q<->pb/mv): pb/mv sibling of this prolongation (piecewise-constant there); keep the child-offset frame
    !! and volume-fraction closure in lockstep.
    impure subroutine s_interpolate_coarse_to_fine()

        integer :: i, bstride

        ! the prolong kernels read the gathered patch's device mirror; the level-1 patch is host-filled by the gather
        ! unpack, so push it once per dispatch (for a level>=2 block the patch was device-produced and this re-push of the
        ! pulled bytes is redundant but harmless)

        do i = 1, sys_size
            $:GPU_UPDATE(device='[amr_cg(i)%sf]')
        end do
        bstride = 1
        if (bubbles_euler) bstride = (eqn_idx%bub%end - eqn_idx%bub%beg + 1)/nb
        do i = 1, sys_size
            ! Lagrangian bubbles: alphas sum to the local liquid fraction beta (not 1), so the sum-to-one closure would corrupt
            ! the EL state; each alpha prolongs plainly instead
            if (num_fluids > 1 .and. (.not. bubbles_lagrange) .and. i >= eqn_idx%adv%beg .and. i <= eqn_idx%adv%end) cycle
            if (chemistry .and. i >= eqn_idx%species%beg .and. i <= eqn_idx%species%end) cycle  ! sum/positivity closure below
            ! QBMM carries a bivariate 6-moment set per R0 bin whose CHyQMOM inversion requires realizability (variance c20 =
            ! m20/m00 - (m10/m00)^2 > 0); per-component minmod prolongation can break that joint constraint, so the whole bub
            ! block is injected piecewise-constant (each child inherits the coarse cell's realizable moment set exactly). Non-QBMM
            ! Euler-Euler bubbles instead floor their positive moments (radius nR, non-polytropic partial pressure npb / vapor
            ! mass nmv); the signed velocity moment nV (offset 1 in each bin's stride) prolongs freely.
            call s_prolong_one_var(amr_cg(i), amr_loc_of(amr_cur), i, &
                                   & pos=bubbles_euler .and. .not. qbmm .and. i >= eqn_idx%bub%beg .and. i <= eqn_idx%bub%end &
                                   & .and. mod(i - eqn_idx%bub%beg, bstride) /= 1, &
                                   & inject=qbmm .and. i >= eqn_idx%bub%beg .and. i <= eqn_idx%bub%end)
        end do
        if (num_fluids > 1 .and. (.not. bubbles_lagrange)) call s_prolong_alphas_closure(amr_cg, amr_loc_of(amr_cur))
        if (chemistry) call s_prolong_species_closure(amr_cg, amr_loc_of(amr_cur))

    end subroutine s_interpolate_coarse_to_fine

    !> Sum-preserving volume-fraction prolongation (num_fluids > 1): fluids adv%beg..adv%end-1 are interpolated with minmod slopes
    !! under a shared per-cell limiter switch (a sign change for any fluid in a dim zeroes that dim's slope for all fluids, so the
    !! closure fluid's effective slope is limited consistently) and clamped to [0,1]; the last fluid closes alpha_n = 1 -
    !! sum(others), so sum(alpha) = 1 on the fine level by construction. For two fluids the closure is also in [0,1]; for >2 fluids
    !! any residual closure undershoot is handled by mpp_lim (required by the checker). Same fine/coarse index mapping as
    !! s_prolong_one_var.
    impure subroutine s_prolong_alphas_closure(qc, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: loc
        integer                                             :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        integer                                             :: rrat, mm, nn, pp, il1, il2, il3, advb, adve
        real(wp)                                            :: xix, xiy, xiz, u0, sx, sy, sz, av, asum
        logical                                             :: shx, shy, shz, d2, d3

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): patch-frame offset.
        ! Device kernel; the shared limiter switch (s_alpha_shared_switch) is inlined verbatim.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, asum, shx, shy, shz, i]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ! shared per-cell limiter switch (inlined s_alpha_shared_switch): per dim, slopes stay on only if no
                    ! fluid's centered differences change sign there (symmetric in the fluids, incl. the closure fluid)
                    shx = .true.; shy = d2; shz = d3
                    do i = advb, adve
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        if ((real(qc(i)%sf(ci + 1, cj, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci - 1, cj, ck), &
                            & wp)) <= 0._wp) shx = .false.
                        if (d2) then
                            if ((real(qc(i)%sf(ci, cj + 1, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj - 1, ck), &
                                & wp)) <= 0._wp) shy = .false.
                        end if
                        if (d3) then
                            if ((real(qc(i)%sf(ci, cj, ck + 1), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj, ck - 1), &
                                & wp)) <= 0._wp) shz = .false.
                        end if
                    end do
                    asum = 0._wp
                    do i = advb, adve - 1
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = 0._wp
                        if (shx) sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (d2 .and. shy) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, cj - 1, &
                            & ck), wp))
                        sz = 0._wp
                        if (d3 .and. shz) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                        amr_cons_st(fi, fj, fk, i, loc) = av
                        asum = asum + av
                    end do
                    amr_cons_st(fi, fj, fk, adve, loc) = 1._wp - asum
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_alphas_closure

    !> Species mass-fraction prolongation closure (chemistry): each partial density rho*Y_k is minmod-prolonged and clamped
    !! non-negative, then all species are rescaled so sum_k(rho*Y_k) equals the (already prolonged) continuity density at the fine
    !! cell. This keeps the fine species realizable (Y_k >= 0, and sum(Y_k) = 1 exactly under the cons->prim recovery rho = sum
    !! rho*Y_k) and consistent with the continuity variable the reaction source reads. Same index mapping as s_prolong_one_var; cont
    !! is prolonged in the main loop before this runs.
    impure subroutine s_prolong_species_closure(qc, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: loc
        integer                                             :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        integer                                             :: rrat, mm, nn, pp, il1, il2, il3, spb, spe, cte
        real(wp)                                            :: xix, xiy, xiz, u0, sx, sy, sz, av, rsum, rscale
        logical                                             :: d2, d3

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): patch-frame offset.
        ! Device kernel; the rescale re-reads only this thread's own cell, so the loop nest is safe.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        spb = eqn_idx%species%beg; spe = eqn_idx%species%end; cte = eqn_idx%cont%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, rsum, rscale, i]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    rsum = 0._wp
                    do i = spb, spe
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (d2) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, cj - 1, ck), wp))
                        sz = 0._wp
                        if (d3) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, ck - 1), wp))
                        av = max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp)
                        amr_cons_st(fi, fj, fk, i, loc) = av
                        rsum = rsum + av
                    end do
                    rscale = real(amr_cons_st(fi, fj, fk, cte, loc), wp)/max(rsum, 1.e-30_wp)
                    do i = spb, spe
                        amr_cons_st(fi, fj, fk, i, loc) = real(amr_cons_st(fi, fj, fk, i, loc), wp)*rscale
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_species_closure

    !> Disblock prolongation. Guard: no-op unless amr.
    impure subroutine s_populate_amr_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        integer                                                :: islot

        if (.not. amr) return
        ! Prolong every block (max_grid_size tiling can make several) from its gathered coarse patch. The P2P gather pulls each
        ! patch's inter-rank coarse cells from neighbour interiors, so no coarse-ghost halo exchange is needed; host q_cons_base
        ! holds
        ! the ICs here (this runs before s_initialize_gpu_vars). All ranks call the gather (P2P); only owners prolong.
        do islot = f_l0_slot(1), amr_num_blocks
            call s_amr_select_slot(islot)
            call s_amr_gather_coarse_patch(q_cons_base, .false.)
            call s_amr_gather_send_flush()  ! this site has blocking semantics
            ! non-polytropic QBMM: gather the coarse pb/mv patch too (all ranks, P2P; owners prolong from it below)
            if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, .false.)
            if (amr_rank_owns_block) then
                ! the prolong is a device kernel (writes the slot in place); no push, since a host->device push here would
                ! clobber the device result with the stale host mirror
                call s_interpolate_coarse_to_fine()
                ! non-polytropic QBMM: seed the block's quadrature side-state from the coarse fields
                if (qbmm .and. .not. polytropic) call s_amr_prolong_pbmv()
            end if
        end do
        if (amr_max_level >= 2) call s_amr_build_static_multilevel(q_cons_base)
        call s_amr_select_slot(f_l0_slot(1))

    end subroutine s_populate_amr_fine

    !> Build the static multi-level hierarchy (amr_regrid_int = 0): nest exactly one level-2 block inside level-1 block 1 by a fixed
    !! geometric inset (a regrid would place it by sensor-on-fine instead), prolong the parent state into it, and keep it persistent
    !! so the advance driver steps it every timestep. The restrict/reflux identity it relies on is protected by the static
    !! multi-level goldens and the runtime conservation-defect probe.
    impure subroutine s_amr_build_static_multilevel(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        integer                                                :: L2, n1, par, inset(3)

        if (amr_max_level < 2) return
        n1 = amr_num_blocks
        if (n1 < 1) return
        ! the static hierarchy nests exactly one level-2 block; without pool room it would silently refine only to level 1 (an
        ! under-resolved but "successful" run). n1 (the level-1 tile count) is only known here, not at checker time, so abort at the
        ! point of failure. Replicated inputs -> every rank takes the same branch (collective-safe).
        if (n1 + 1 > l0_slot_off + amr_max_fine) call s_mpi_abort('amr static multi-level (amr_max_level > 1, ' &
            & // 'amr_regrid_int = 0): amr_max_blocks is too small to nest the level-2 block (need >= level-1 block count + 1); ' &
            & // 'increase amr_max_blocks')
        L2 = n1 + 1
        ! parent is the first fine block, f_l0_slot(1), not slot 1, which under coexist is the first level-0 tile. Insetting a
        ! tile instead would put the level-2 box in the wrong place and size it off the tile (a plausible-looking box that
        ! silently corrupts the run).
        par = f_l0_slot(1)
        inset = 0
        inset(1) = max((amr_region_hi_all(1, par) - amr_region_lo_all(1, par) + 1)/4, amr_cpat_mar)
        if (n_glb > 0) inset(2) = max((amr_region_hi_all(2, par) - amr_region_lo_all(2, par) + 1)/4, amr_cpat_mar)
        if (p_glb > 0) inset(3) = max((amr_region_hi_all(3, par) - amr_region_lo_all(3, par) + 1)/4, amr_cpat_mar)
        amr_region_lo_all(:,L2) = amr_region_lo_all(:,par) + inset
        amr_region_hi_all(:,L2) = amr_region_hi_all(:,par) - inset
        ! Guard the fixed-inset box against configs this single-block static builder cannot represent; the dynamic regrid path has
        ! the analogous checks (proper-nesting skip + amr_maxc_fit/2 clamp), but the static path bypasses them. Replicated inputs ->
        ! every rank takes the same branch (collective-safe). (a) a level-1 block smaller than 2*inset inverts the box; (b) a
        ! level-2
        ! L0-extent > amr_maxc_fit/2 makes its parent-fine transverse extent (2*L0) overrun the creg register (allocated
        ! 0:amr_maxc_fit-1), a silent out-of-bounds device write in the L2->L1 reflux capture.
        if (amr_region_lo_all(1, L2) > amr_region_hi_all(1, L2) .or. (n_glb > 0 .and. amr_region_lo_all(2, &
            & L2) > amr_region_hi_all(2, L2)) .or. (p_glb > 0 .and. amr_region_lo_all(3, L2) > amr_region_hi_all(3, &
            & L2))) call s_mpi_abort('amr static multi-level: level-1 block 1 is too small to nest a level-2 block (the fixed ' &
            & // 'inset inverts the box); enlarge the base amr block or reduce amr_cpat_mar')
        if (amr_ref_ratio*(amr_region_hi_all(1, L2) - amr_region_lo_all(1, &
            & L2) + 1) > amr_maxc_fit(1) .or. (n_glb > 0 .and. amr_ref_ratio*(amr_region_hi_all(2, L2) - amr_region_lo_all(2, &
            & L2) + 1) > amr_maxc_fit(2)) .or. (p_glb > 0 .and. amr_ref_ratio*(amr_region_hi_all(3, L2) - amr_region_lo_all(3, &
            & L2) + 1) > amr_maxc_fit(3))) &
            & call s_mpi_abort('amr static multi-level: the nested level-2 block exceeds the per-rank scratch cap ' &
            & // '(2*L0-extent > amr_maxc_fit); static multi-level does not tile the level-2 block - use a smaller base amr ' &
            & // 'block or the dynamic regrid path (amr_regrid_int > 0)')
        amr_block_level(L2) = 2
        amr_block_owner(L2) = amr_block_owner(par); amr_myblk_dirty = .true.
        amr_num_blocks = L2; amr_num_levels = 2
        call s_amr_reconcile_slots()
        amr_cur = L2
        call s_set_amr_fine_geometry(amr_region_lo_all(:,L2), amr_region_hi_all(:,L2))
        call s_amr_reduce_xchg_flag()
        call s_amr_gather_coarse_patch(q_cons_base, .false.)  ! q_coarse ignored for level>=2 (reads the parent block); pass the
        call s_amr_gather_send_flush()  ! this site has blocking semantics
        ! always-allocated base field, not amr_slots(1) (the parent slot is unallocated on a non-owner rank at np>1)
        if (amr_rank_owns_block) then
            ! the prolong is a device kernel: the persistent L2 block's device q_cons is valued in place (a host->device push
            ! here would clobber the device result)
            call s_interpolate_coarse_to_fine()
        end if
        ! persistent L2 block: keep the level-2 block in the active set (amr_num_blocks = L2, amr_num_levels = 2) so the advance
        ! driver steps it across timesteps; no free/revert.
        ! restore amr_cg + the patch frame (amr_cpat_off) to the first fine block: the L2 gather above overwrote them with the
        ! parent-fine frame, and the normal single-block conservation check that follows reads that block's frame. f_l0_slot(1),
        ! not slot 1: under coexist slot 1 is a level-0 tile, and selecting it here would leave the grid globals describing tile
        ! geometry for the rest of init, so s_initialize_weno_module (m_start_up, called after this) would size its device-mapped
        ! coefficient tables off the wrong bounds.
        call s_amr_select_slot(f_l0_slot(1))
        call s_amr_gather_coarse_patch(q_cons_base, .false.)
        call s_amr_gather_send_flush()  ! this site has blocking semantics

    end subroutine s_amr_build_static_multilevel

    !> Volume-weighted restriction: each covered coarse cell = volume-weighted average of its amr_ref_ratio^d fine children (equal
    !! weight on Cartesian grids where children share a volume; radius-weighted by fine y_cc on cyl_coord, where cell volume ~
    !! radius: amr_rvw, single-sourced so device and host paths agree bit-for-bit). Writes the caller's coarse target: in production
    !! the level-0 state q_cons_ts(1)%vf (the deliberate fold-back of fine data each step, plus coarse pb/mv for non-polytropic
    !! QBMM); init-time diagnostics pass a scratch buffer instead. Device kernel.
    impure subroutine s_restrict_fine_to_coarse(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        integer :: nchild, rr, dj_hi, dk_hi, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr
        integer :: rlo(3), rhi(3), ilo(3), ihi(3), bl(3), bh(3)
        real(wp), allocatable :: sbuf(:,:), rbuf(:)
        integer, allocatable :: reqs(:), drank(:)

        if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_tic()

        ! multi-level: a level>=2 block folds back into its parent block's fine array (the coarse side of level l is level l-1),
        ! not the L0 coarse_tgt. Same restriction kernel, targeted at the parent in the parent-fine frame. When child and parent
        ! sit on different ranks the fold is a P2P pair, so both participants must enter or the receiver never posts.
        if (amr_block_level(amr_cur) >= 2) then
            if (amr_rank_owns_block .or. amr_block_owner(f_amr_parent_block(amr_cur)) == proc_rank) then
                call s_amr_restrict_to_parent()
            end if
            if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_toc()
            return
        end if

        ! whole-block-per-rank fold-back: the block owner restricts its fine block to coarse averages over the covered cells
        ! [region_lo:region_hi] and scatters them point-to-point to the coarse-cell owners: the owner overwrites the covered cells
        ! it holds locally and sends each other coarse-owner exactly its covered slice (all sys_size in one message). Covered
        ! cells are in-domain (no ghosts), so each is owned by exactly one interior owner. At np=1 the owner owns every covered
        ! cell, sends nothing, and overwrites locally with the same child-sum.
        rr = amr_slots(amr_cur)%amr_ref_ratio
        nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        rlo = 0; rhi = 0
        rlo(1) = amr_region_lo_all(1, amr_cur); rhi(1) = amr_region_hi_all(1, amr_cur)
        if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, amr_cur); rhi(2) = amr_region_hi_all(2, amr_cur); end if
        if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, amr_cur); rhi(3) = amr_region_hi_all(3, amr_cur); end if
        owner = amr_block_owner(amr_cur)
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        maxsz = sys_size*(rhi(1) - rlo(1) + 1)*(rhi(2) - rlo(2) + 1)*(rhi(3) - rlo(3) + 1)

        ! cyl_coord: fine radial volume weights = this block's fine cell-center radii, pushed to device for the restriction kernels.
        ! Only the owner restricts (device overwrite + device scatter pack), so only the owner needs them; single-sourced from y_cc
        ! so the owner-local and scattered child-averages are bit-identical.
        if (cyl_coord .and. proc_rank == owner) then
            amr_rvw(0:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%y_cc(0:amr_slots(amr_cur)%n)
            $:GPU_UPDATE(device='[amr_rvw]')
        end if

        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()

        if (proc_rank == owner) then
            ! overwrite the covered cells this rank owns, then send each other coarse-owner its covered slice
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            if (num_procs == 1) then
                ! np=1 device-native fold-back: restrict the fine block (device) into the coarse (device) over the covered cells
                ! only, with no host round-trip. Never push the whole coarse array back to the device here: that would clobber
                ! the device-advanced non-covered coarse cells with the stale host copy, a GPU-only divergence (invisible on CPU
                ! where host==device) that IGR/MHD/acoustic amplify. The owner holds every covered cell at np=1.
                if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_overwrite_device_sf(coarse_tgt, &
                    & amr_loc_of(amr_cur), bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
                if (qbmm .and. .not. polytropic .and. amr_rank_owns_block) call s_restrict_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, &
                    & amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf)
                if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_toc()
                return
            end if
            ! owner-local covered cells: restrict fine(device) -> coarse(device) touching only those cells (no whole-coarse device
            ! push, which would clobber the device-advanced non-covered coarse cells; same hazard as at np=1)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_overwrite_device_sf(coarse_tgt, &
                & amr_loc_of(amr_cur), bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
            ! cached destination list (every listed rank's interior overlaps the region by construction)
            nsrc = 0
            do idx = 1, amr_ovl_scatter_n(amr_cur)
                if (amr_ovl_scatter(idx, amr_cur) /= owner) nsrc = nsrc + 1
            end do
            if (nsrc > 0) then
                allocate (sbuf(maxsz, nsrc), reqs(nsrc), drank(nsrc))
                nsrc = 0
                do idx = 1, amr_ovl_scatter_n(amr_cur)
                    r = amr_ovl_scatter(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_interior(r, ilo, ihi)
                    call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
                    nsrc = nsrc + 1; drank(nsrc) = r
                    boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    ! pack this destination's covered slice on the device: restrict averages straight into the wire buffer (same
                    ! child-sum order and wp values as the device overwrite above), with no full-field host pull
                    call s_amr_restrict_pack_device(amr_loc_of(amr_cur), bl, bh, rlo, rr, dj_hi, dk_hi, nchild, sbuf(1:boxsz,nsrc))
#ifdef MFC_MPI
                    call s_xa_rec(XA_F7A_SND, 1, boxsz, amr_cur)
                    call MPI_ISEND(sbuf(1, nsrc), boxsz, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
#ifdef MFC_MPI
                call s_wait_tic()
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
                call s_wait_toc(WT_RESTR)
#endif
                deallocate (sbuf, reqs, drank)
            end if
        else
            ! coarse-owner: if I hold covered cells, receive my slice from the owner and overwrite my local coarse
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                allocate (rbuf(boxsz))
#ifdef MFC_MPI
                call s_xa_rec(XA_F7A_RCV, 2, boxsz, amr_cur)
                call s_wait_tic()
                call MPI_RECV(rbuf, boxsz, mpi_p, owner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                call s_wait_toc(WT_RESTR)
#endif
                ! Device unpack of the covered box, writing only those cells (a whole-array push would clobber the device-advanced
                ! non-covered coarse cells with this rank's stale host copy). This must not be a host unpack followed by a
                ! strided GPU_UPDATE(device=) of the box: a non-contiguous 3-D array section in an OpenMP target update is copied
                ! by AMD flang as size(section) contiguous elements, so only the first row lands on the cells it names and the
                ! remainder silently overwrites neighbouring cells with stale host data (only at np >= 2 with a block whose owner
                ! holds none of its covered cells). The wire layout (ci fastest, then cj, ck, i) is exactly
                ! s_l0_pack_unpack_block's, so it unpacks s_amr_restrict_pack_device's buffer as-is.
                call s_l0_pack_unpack_block_sf(coarse_tgt, bl(1) - o1, bl(2) - o2, bl(3) - o3, bh(1) - bl(1), bh(2) - bl(2), &
                                               & bh(3) - bl(3), rbuf, .false.)
                deallocate (rbuf)
            end if
        end if

        ! non-polytropic QBMM: distributed pb/mv fold-back; the owner restricts the covered cells it holds and scatters each other
        ! coarse-owner its slice (mirror of the q_cons scatter above). Called on all ranks so the P2P send/recv pair up (np=1
        ! handled
        ! by the direct s_restrict_pbmv in the num_procs==1 branch above, which returns before reaching here)
        if (qbmm .and. .not. polytropic) call s_amr_scatter_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf)
        if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_toc()

    end subroutine s_restrict_fine_to_coarse

    !> Multi-level restriction: fold the current level>=2 block's fine averages back into its parent block's fine array over the
    !! covered cells. Same child-sum kernel as the L0 fold-back, targeted at the parent in the parent-fine frame (amr_isect already
    !! parent-fine; offset 0 = the parent's local fine indexing). Local when child and parent share an owner; otherwise a P2P pair.
    impure subroutine s_amr_restrict_to_parent()

        integer               :: pblk, rr, nchild, dj_hi, dk_hi, cowner, powner, boxsz, ierr
        integer               :: plo(3), phi(3)
        real(wp), allocatable :: xbuf(:)

        ! Reads amr_rvw's device copy without a GPU_UPDATE, which is safe only while cyl_coord + amr_max_level > 1 is
        ! checker-gated (this path never runs under cyl_coord). If that gate lifts, refresh amr_rvw here first (see its
        ! declaration).

        pblk = f_amr_parent_block(amr_cur)
        cowner = amr_block_owner(amr_cur); powner = amr_block_owner(pblk)
        if (proc_rank /= cowner .and. proc_rank /= powner) return  ! not a participant

        ! Same replicated-metadata box as the gather, so the folding child and the receiving parent agree without a handshake.
        call s_amr_parent_foot(amr_cur, pblk, plo, phi)
        if (plo(1) > phi(1) .or. plo(2) > phi(2) .or. plo(3) > phi(3)) return  ! empty footprint

        rr = amr_ref_ratio
        nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)

        if (powner == cowner) then
            ! co-located (np=1, or a co-located tower): fold straight into the parent.
            call s_amr_restrict_overwrite_device_st(amr_loc_of(pblk), amr_loc_of(amr_cur), plo, phi, 0, 0, 0, plo, rr, dj_hi, &
                                                    & dk_hi, nchild)
            return
        end if

#ifdef MFC_MPI
        ! Split ownership: the child restricts locally and ships coarse cells (rr**num_dims fewer values than shipping its fine
        ! block), which restriction being an overwrite (not an accumulate) makes correct. Reuses the L0<->L1 scatter's pack and
        ! unpack; their wire layout (ci fastest, then cj, ck, i) is compatible.
        boxsz = sys_size*(phi(1) - plo(1) + 1)*(phi(2) - plo(2) + 1)*(phi(3) - plo(3) + 1)
        allocate (xbuf(boxsz))
        if (proc_rank == cowner) then
            call s_amr_restrict_pack_device(amr_loc_of(amr_cur), plo, phi, plo, rr, dj_hi, dk_hi, nchild, xbuf)
            call s_xa_rec(XA_F7B_SND, 1, boxsz, amr_cur)
            call MPI_SEND(xbuf, boxsz, mpi_p, powner, amr_cur, MPI_COMM_WORLD, ierr)
        else
            call s_xa_rec(XA_F7B_RCV, 2, boxsz, amr_cur)
            call s_wait_tic()
            call MPI_RECV(xbuf, boxsz, mpi_p, cowner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            call s_wait_toc(WT_RESTR)
            ! Device unpack of just the covered box; never a host unpack plus a strided GPU_UPDATE (see the L0 scatter's note: AMD
            ! flang copies a non-contiguous 3-D section as contiguous elements and silently corrupts neighbouring cells).
            call s_l0_pack_unpack_block_st(amr_loc_of(pblk), plo(1), plo(2), plo(3), phi(1) - plo(1), phi(2) - plo(2), &
                                           & phi(3) - plo(3), xbuf, .false.)
        end if
        deallocate (xbuf)
#endif

    end subroutine s_amr_restrict_to_parent

    !> The lock-step restrict fold as plan-based waves, on the lock-step np>1 path: per level (finest first) the split-ownership
    !! child->parent folds are one aggregated message per (child-owner, parent-owner) peer (s_amr_restrict_parent_wave), followed by
    !! the per-box reflux-to-parent applies (their freg was already exchanged by s_amr_freg_wave); then the level-1 -> L0
    !! covered-cell scatter is one aggregated message per (owner, coarse-owner) peer (s_amr_restrict_l1_wave). Both sides derive
    !! identical transfer lists from replicated metadata (the region-box x rank-interior intersections), so the wire layout needs no
    !! handshake and the F7 family words are exact. Level order (finest first) preserves child-before-parent folding; within a level
    !! the covered targets are disjoint, and the restrict/reflux interleave is order-free because sibling-shared faces carry weight
    !! 0 (s_amr_sibling_face_weights). Subcycle and np=1 use the per-box loop.
    impure subroutine s_amr_restrict_wave(coarse_tgt, dt_reflux)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        real(wp), intent(in)                                   :: dt_reflux
        integer                                                :: lev, k, kk, io, ifc, ko, kf

        call s_amr_refresh_lists()
        do lev = amr_max_level, 2, -1
            if (relax) then
                ! s_amr_relax_fine returns unless amr_rank_owns_block (the multi-owner amr_owns_all notion) -> amr_own_blk
                do kk = amr_n_own, 1, -1
                    k = amr_own_blk(kk)
                    if (amr_block_level(k) /= lev) cycle
                    call s_amr_select_slot(k)
                    call s_amr_relax_fine()
                end do
            end if
            call s_phase_tic(PH_RESTR); call s_phase_tic(PH_RSREST)
            call s_amr_restrict_parent_wave(lev)
            call s_phase_toc(PH_RSREST); call s_phase_toc(PH_RESTR)
            ! s_amr_reflux_to_parent returns unless `own_child .or. own_parent`: own_child = amr_rank_owns_block, the
            ! multi-owner amr_owns_all notion (amr_own_blk, not amr_my_blk); own_parent = the level-lev children of my parents,
            ! i.e. amr_fch_blk plus the ones I own myself, already in amr_own_blk. The two lists overlap, and a duplicate visit
            ! would reflux twice, so walk their union: both lists are ascending, and a two-cursor merge from the tails visits each
            ! block once, in descending block order (the parent registers accumulate in that order).
            io = amr_n_own; ifc = amr_n_fch
            do
                do while (io >= 1)
                    if (amr_block_level(amr_own_blk(io)) == lev) exit
                    io = io - 1
                end do
                do while (ifc >= 1)
                    if (amr_block_level(amr_fch_blk(ifc)) == lev) exit
                    ifc = ifc - 1
                end do
                ko = 0; kf = 0
                if (io >= 1) ko = amr_own_blk(io)
                if (ifc >= 1) kf = amr_fch_blk(ifc)
                k = max(ko, kf)
                if (k == 0) exit
                if (ko == k) io = io - 1
                if (kf == k) ifc = ifc - 1
                call s_amr_select_slot(k)
                call s_phase_tic(PH_RESTR); call s_phase_tic(PH_RSRFP)
                call s_amr_reflux_to_parent(dt_reflux, .false.)
                call s_phase_toc(PH_RSRFP); call s_phase_toc(PH_RESTR)
            end do
        end do
        if (relax) then
            do kk = amr_n_own, 1, -1  ! same amr_rank_owns_block predicate as the level>=2 loop above
                k = amr_own_blk(kk)
                if (amr_block_level(k) /= 1) cycle
                call s_amr_select_slot(k)
                call s_amr_relax_fine()
            end do
        end if
        call s_phase_tic(PH_RESTR); call s_phase_tic(PH_RSREST)
        call s_amr_restrict_l1_wave(coarse_tgt)
        call s_phase_toc(PH_RSREST); call s_phase_toc(PH_RESTR)
        ! non-polytropic QBMM pb/mv fold-back keeps its per-box pairing (every rank walks the same reverse order, so the
        ! blocking pairs match exactly as they did inside the per-box fold)
        if (qbmm .and. .not. polytropic) then
            do k = amr_num_blocks, 1, -1
                if (amr_block_level(k) /= 1) cycle
                call s_amr_select_slot(k)
                call s_amr_scatter_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf)
            end do
        end if
        call s_amr_select_slot(1)

    end subroutine s_amr_restrict_wave

    !> One per-level wave of the split-ownership level>=2 child->parent restrict folds (the F7B pairs). Co-located folds run inline
    !! on the child-owner (bit-for-bit the per-box kernel); every cross-rank fold ships its parent-frame covered box in one
    !! aggregated message per (child-owner, parent-owner) peer. Both sides walk the same replicated block list ascending with
    !! per-peer running offsets, so the wire layout agrees with no metadata exchange.
    impure subroutine s_amr_restrict_parent_wave(lev)

        integer, intent(in) :: lev

#ifdef MFC_MPI
        integer :: k, pblk, cowner, powner, rr, nchild, dj_hi, dk_hi, ierr, ip, idx, r, cnt, boff, qbase, nreq, tq, sq, kk
        integer :: plo(3), phi(3)

        rr = amr_ref_ratio
        nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        call s_amr_m1_wave_open(7)
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if
        ! send plan + co-located folds (child-owner side)
        amr_fw_snx = 0; amr_fw_snp = 0
        call s_amr_refresh_my_blocks()
        call s_amr_refresh_lists()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            cowner = amr_block_owner(k)
            pblk = amr_parent_blk(k)
            powner = amr_block_owner(pblk)
            call s_amr_parent_foot(k, pblk, plo, phi)
            if (plo(1) > phi(1) .or. plo(2) > phi(2) .or. plo(3) > phi(3)) cycle
            if (cowner == powner) then
                call s_amr_restrict_overwrite_device_st(amr_loc_of(pblk), amr_loc_of(k), plo, phi, 0, 0, 0, plo, rr, dj_hi, &
                                                        & dk_hi, nchild)
                cycle
            end if
            if (amr_fw_map(powner) == 0) then
                amr_fw_snp = amr_fw_snp + 1
                call s_amr_fw_szi(amr_fw_sprank, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqsz, amr_fw_snp)
                call s_amr_fw_szi(amr_fw_snxp, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqbase, amr_fw_snp)
                amr_fw_map(powner) = amr_fw_snp
                amr_fw_sprank(amr_fw_snp) = powner
            end if
            cnt = sys_size*(phi(1) - plo(1) + 1)*(phi(2) - plo(2) + 1)*(phi(3) - plo(3) + 1)
            amr_fw_snx = amr_fw_snx + 1
            call s_amr_fw_szi(amr_fw_sblk, amr_fw_snx); call s_amr_fw_szi3(amr_fw_sbl, amr_fw_snx)
            call s_amr_fw_szi3(amr_fw_sbh, amr_fw_snx); call s_amr_fw_szi(amr_fw_spi, amr_fw_snx)
            call s_amr_fw_szi(amr_fw_sqo, amr_fw_snx); call s_amr_fw_szi(amr_fw_spo, amr_fw_snx)
            amr_fw_sblk(amr_fw_snx) = k
            amr_fw_sbl(:,amr_fw_snx) = plo; amr_fw_sbh(:,amr_fw_snx) = phi
            amr_fw_spo(amr_fw_snx) = cnt
            amr_fw_spi(amr_fw_snx) = amr_fw_map(powner)
            amr_fw_sqo(amr_fw_snx) = amr_fw_pq(powner) + amr_fw_nx(powner)*XA_NH
            amr_fw_pq(powner) = amr_fw_pq(powner) + cnt
            amr_fw_nx(powner) = amr_fw_nx(powner) + 1
        end do
        qbase = 0
        do ip = 1, amr_fw_snp
            r = amr_fw_sprank(ip)
            amr_fw_snxp(ip) = amr_fw_nx(r)
            amr_fw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_sqbase(ip) = qbase; qbase = qbase + amr_fw_sqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_sq, qbase, amr_fw_dev)
        ! recv plan (parent-owner side): the same replicated walk, so per-peer transfer order matches the sender's.
        ! amr_fch_blk holds exactly these survivors across all levels >= 2, ascending; the level filter narrows to lev
        amr_fw_rnx = 0; amr_fw_rnp = 0
        do kk = 1, amr_n_fch
            k = amr_fch_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            pblk = amr_parent_blk(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. proc_rank /= powner) cycle
            call s_amr_parent_foot(k, pblk, plo, phi)
            if (plo(1) > phi(1) .or. plo(2) > phi(2) .or. plo(3) > phi(3)) cycle
            if (amr_fw_map(cowner) == 0) then
                amr_fw_rnp = amr_fw_rnp + 1
                call s_amr_fw_szi(amr_fw_rprank, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqsz, amr_fw_rnp)
                call s_amr_fw_szi(amr_fw_rnxp, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqbase, amr_fw_rnp)
                amr_fw_map(cowner) = amr_fw_rnp
                amr_fw_rprank(amr_fw_rnp) = cowner
            end if
            cnt = sys_size*(phi(1) - plo(1) + 1)*(phi(2) - plo(2) + 1)*(phi(3) - plo(3) + 1)
            amr_fw_rnx = amr_fw_rnx + 1
            call s_amr_fw_szi(amr_fw_rblk, amr_fw_rnx); call s_amr_fw_szi3(amr_fw_rbl, amr_fw_rnx)
            call s_amr_fw_szi3(amr_fw_rbh, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpi, amr_fw_rnx)
            call s_amr_fw_szi(amr_fw_rqo, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpo, amr_fw_rnx)
            amr_fw_rblk(amr_fw_rnx) = k
            amr_fw_rbl(:,amr_fw_rnx) = plo; amr_fw_rbh(:,amr_fw_rnx) = phi
            amr_fw_rpo(amr_fw_rnx) = cnt
            amr_fw_rpi(amr_fw_rnx) = amr_fw_map(cowner)
            amr_fw_rqo(amr_fw_rnx) = amr_fw_pq(cowner) + amr_fw_nx(cowner)*XA_NH
            amr_fw_pq(cowner) = amr_fw_pq(cowner) + cnt
            amr_fw_nx(cowner) = amr_fw_nx(cowner) + 1
        end do
        qbase = 0
        do ip = 1, amr_fw_rnp
            r = amr_fw_rprank(ip)
            amr_fw_rnxp(ip) = amr_fw_nx(r)
            amr_fw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rqbase(ip) = qbase; qbase = qbase + amr_fw_rqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_rq, qbase, amr_fw_dev)
        nreq = amr_fw_snp + amr_fw_rnp
        if (nreq == 0) return
        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
        nreq = 0
        do ip = 1, amr_fw_rnp
            sq = f_amr_m1_seq(amr_fw_rprank(ip), 2); tq = f_amr_m1_tag(7, sq)
            call s_xa_rec(XA_F7BW_RCV, 2, amr_fw_rqsz(ip) - amr_fw_rnxp(ip)*XA_NH, tq, peer=amr_fw_rprank(ip), &
                          & key=amr_fw_rnxp(ip), seq=sq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rqsz(ip)
            if (amr_fw_dev) then
                #:call GPU_HOST_DATA(use_device_addr='[amr_fw_rq]')
                    call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                #:endcall GPU_HOST_DATA
            else
                call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
        end do
        do idx = 1, amr_fw_snx
            cnt = amr_fw_spo(idx)
            boff = amr_fw_sqbase(amr_fw_spi(idx)) + amr_fw_sqo(idx)
            call s_amr_restrict_pack_device(amr_loc_of(amr_fw_sblk(idx)), amr_fw_sbl(:,idx), amr_fw_sbh(:,idx), amr_fw_sbl(:, &
                                            & idx), rr, dj_hi, dk_hi, nchild, amr_fw_sq(boff + XA_NH + 1:boff + XA_NH + cnt))
            if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F7BW_SND, amr_fw_sblk(idx), amr_fw_sbl(:,idx), &
                & amr_fw_sbh(:,idx))
        end do
        do ip = 1, amr_fw_snp
            sq = f_amr_m1_seq(amr_fw_sprank(ip), 1); tq = f_amr_m1_tag(7, sq)
            call s_xa_rec(XA_F7BW_SND, 1, amr_fw_sqsz(ip) - amr_fw_snxp(ip)*XA_NH, tq, peer=amr_fw_sprank(ip), &
                          & key=amr_fw_snxp(ip), seq=sq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = -1
            if (amr_fw_dev) then
                #:call GPU_HOST_DATA(use_device_addr='[amr_fw_sq]')
                    call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                #:endcall GPU_HOST_DATA
            else
                call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
        end do
#ifdef MFC_DEBUG
        block
            integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
            call s_wait_tic()
            call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
            call s_wait_toc(WT_RESTR)
            do q = 1, nreq
                if (amr_fw_reqw(q) < 0) cycle
                call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                @:ASSERT(gotw == amr_fw_reqw(q), "restrict parent wave: received message length differs from the plan")
            end do
        end block
#else
        call s_wait_tic()
        call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
        call s_wait_toc(WT_RESTR)
#endif
        do idx = 1, amr_fw_rnx
            cnt = amr_fw_rpo(idx)
            boff = amr_fw_rqbase(amr_fw_rpi(idx)) + amr_fw_rqo(idx)
            if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(boff + 1:boff + XA_NH), XA_F7BW_SND, amr_fw_rblk(idx), amr_fw_rbl(:, &
                & idx), amr_fw_rbh(:,idx))
            ! Device unpack of the covered box only (the strided-update flang trap; see s_restrict_fine_to_coarse)
            call s_l0_pack_unpack_block_st(amr_loc_of(amr_parent_blk(amr_fw_rblk(idx))), amr_fw_rbl(1, idx), amr_fw_rbl(2, idx), &
                                           & amr_fw_rbl(3, idx), amr_fw_rbh(1, idx) - amr_fw_rbl(1, idx), amr_fw_rbh(2, &
                                           & idx) - amr_fw_rbl(2, idx), amr_fw_rbh(3, idx) - amr_fw_rbl(3, idx), &
                                           & amr_fw_rq(boff + XA_NH + 1:boff + XA_NH + cnt), .false.)
        end do
#endif

    end subroutine s_amr_restrict_parent_wave

    !> The level-1 -> L0 covered-cell scatter (F7A) as one wave: every owned level-1 block's covered slabs for every listed
    !! coarse-owner ship in one aggregated message per peer; the owner-local covered overwrite and the (cyl_coord) amr_rvw push stay
    !! grouped per block during the pack walk. The receiver plan is my-interior x region(k) over the level-1 blocks I do not own, by
    !! construction (s_amr_ranks_overlapping) exactly the sender's list membership.
    impure subroutine s_amr_restrict_l1_wave(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt

#ifdef MFC_MPI
        integer :: k, owner, rr, nchild, dj_hi, dk_hi, ierr, ip, idx, r, cnt, boff, qbase, nreq, tq, sq, o1, o2, o3, cur, kk
        integer :: rlo(3), rhi(3), ilo(3), ihi(3), milo(3), mihi(3), bl(3), bh(3)

        call s_amr_m1_wave_open(6)
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as the per-box path)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if
        call s_amr_rank_interior(proc_rank, milo, mihi)
        ! send plan (block-owner side): the same (interior x region) covered slabs the per-box path sent, k-grouped
        amr_fw_snx = 0; amr_fw_snp = 0
        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            rlo = 0; rhi = 0
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k); end if
            if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k); end if
            do idx = 1, amr_ovl_scatter_n(k)
                r = amr_ovl_scatter(idx, k)
                if (r == proc_rank) cycle
                call s_amr_rank_interior(r, ilo, ihi)
                call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
                if (bl(1) > bh(1) .or. bl(2) > bh(2) .or. bl(3) > bh(3)) cycle
                if (amr_fw_map(r) == 0) then
                    amr_fw_snp = amr_fw_snp + 1
                    call s_amr_fw_szi(amr_fw_sprank, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqsz, amr_fw_snp)
                    call s_amr_fw_szi(amr_fw_snxp, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqbase, amr_fw_snp)
                    amr_fw_map(r) = amr_fw_snp
                    amr_fw_sprank(amr_fw_snp) = r
                end if
                cnt = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                amr_fw_snx = amr_fw_snx + 1
                call s_amr_fw_szi(amr_fw_sblk, amr_fw_snx); call s_amr_fw_szi3(amr_fw_sbl, amr_fw_snx)
                call s_amr_fw_szi3(amr_fw_sbh, amr_fw_snx); call s_amr_fw_szi(amr_fw_spi, amr_fw_snx)
                call s_amr_fw_szi(amr_fw_sqo, amr_fw_snx); call s_amr_fw_szi(amr_fw_spo, amr_fw_snx)
                amr_fw_sblk(amr_fw_snx) = k
                amr_fw_sbl(:,amr_fw_snx) = bl; amr_fw_sbh(:,amr_fw_snx) = bh
                amr_fw_spo(amr_fw_snx) = cnt
                amr_fw_spi(amr_fw_snx) = amr_fw_map(r)
                amr_fw_sqo(amr_fw_snx) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
                amr_fw_pq(r) = amr_fw_pq(r) + cnt
                amr_fw_nx(r) = amr_fw_nx(r) + 1
            end do
        end do
        qbase = 0
        do ip = 1, amr_fw_snp
            r = amr_fw_sprank(ip)
            amr_fw_snxp(ip) = amr_fw_nx(r)
            amr_fw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_sqbase(ip) = qbase; qbase = qbase + amr_fw_sqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_sq, qbase, amr_fw_dev)
        ! recv plan (coarse-owner side): my interior x region(k) over level-1 blocks I do not own
        amr_fw_rnx = 0; amr_fw_rnp = 0
        ! walk the cached receive list, not every block. The list carries exactly the blocks that pass
        ! level + not-mine + overlap; bl/bh are recomputed because the body needs them.
        call s_amr_refresh_lists()
        do kk = 1, amr_n_l1r
            k = amr_l1r_blk(kk)
            owner = amr_block_owner(k)
            rlo = 0; rhi = 0
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k); end if
            if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k); end if
            call s_amr_box_isect(rlo, rhi, milo, mihi, bl, bh)
            if (amr_fw_map(owner) == 0) then
                amr_fw_rnp = amr_fw_rnp + 1
                call s_amr_fw_szi(amr_fw_rprank, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqsz, amr_fw_rnp)
                call s_amr_fw_szi(amr_fw_rnxp, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqbase, amr_fw_rnp)
                amr_fw_map(owner) = amr_fw_rnp
                amr_fw_rprank(amr_fw_rnp) = owner
            end if
            cnt = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
            amr_fw_rnx = amr_fw_rnx + 1
            call s_amr_fw_szi(amr_fw_rblk, amr_fw_rnx); call s_amr_fw_szi3(amr_fw_rbl, amr_fw_rnx)
            call s_amr_fw_szi3(amr_fw_rbh, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpi, amr_fw_rnx)
            call s_amr_fw_szi(amr_fw_rqo, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpo, amr_fw_rnx)
            amr_fw_rblk(amr_fw_rnx) = k
            amr_fw_rbl(:,amr_fw_rnx) = bl; amr_fw_rbh(:,amr_fw_rnx) = bh
            amr_fw_rpo(amr_fw_rnx) = cnt
            amr_fw_rpi(amr_fw_rnx) = amr_fw_map(owner)
            amr_fw_rqo(amr_fw_rnx) = amr_fw_pq(owner) + amr_fw_nx(owner)*XA_NH
            amr_fw_pq(owner) = amr_fw_pq(owner) + cnt
            amr_fw_nx(owner) = amr_fw_nx(owner) + 1
        end do
        qbase = 0
        do ip = 1, amr_fw_rnp
            r = amr_fw_rprank(ip)
            amr_fw_rnxp(ip) = amr_fw_nx(r)
            amr_fw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rqbase(ip) = qbase; qbase = qbase + amr_fw_rqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_rq, qbase, amr_fw_dev)
        nreq = amr_fw_snp + amr_fw_rnp
        call s_amr_fw_szi(amr_fw_req, max(nreq, 1)); call s_amr_fw_szi(amr_fw_reqw, max(nreq, 1))
        nreq = 0
        do ip = 1, amr_fw_rnp
            sq = f_amr_m1_seq(amr_fw_rprank(ip), 2); tq = f_amr_m1_tag(6, sq)
            call s_xa_rec(XA_F7W_RCV, 2, amr_fw_rqsz(ip) - amr_fw_rnxp(ip)*XA_NH, tq, peer=amr_fw_rprank(ip), &
                          & key=amr_fw_rnxp(ip), seq=sq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rqsz(ip)
            if (amr_fw_dev) then
                #:call GPU_HOST_DATA(use_device_addr='[amr_fw_rq]')
                    call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                #:endcall GPU_HOST_DATA
            else
                call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
        end do
        ! owner-local covered overwrites + device packs, grouped per owned block: amr_rvw is a single device mirror, so a
        ! block's (cyl_coord) radii push must immediately precede that block's overwrite/pack kernels; the transfer list is
        ! k-grouped by construction, so a monotone cursor drains each block's sends inside its group
        cur = 1
        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            rr = amr_slots(k)%amr_ref_ratio
            nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
            dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
            if (cyl_coord) then
                amr_rvw(0:amr_slots(k)%n) = amr_slots(k)%y_cc(0:amr_slots(k)%n)
                $:GPU_UPDATE(device='[amr_rvw]')
            end if
            rlo = 0; rhi = 0
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k); end if
            if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k); end if
            call s_amr_box_isect(rlo, rhi, milo, mihi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_overwrite_device_sf(coarse_tgt, &
                & amr_loc_of(k), bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
            do while (cur <= amr_fw_snx)
                if (amr_fw_sblk(cur) /= k) exit
                cnt = amr_fw_spo(cur)
                boff = amr_fw_sqbase(amr_fw_spi(cur)) + amr_fw_sqo(cur)
                call s_amr_restrict_pack_device(amr_loc_of(k), amr_fw_sbl(:,cur), amr_fw_sbh(:,cur), rlo, rr, dj_hi, dk_hi, &
                                                & nchild, amr_fw_sq(boff + XA_NH + 1:boff + XA_NH + cnt))
                if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F7W_SND, k, amr_fw_sbl(:,cur), &
                    & amr_fw_sbh(:,cur))
                cur = cur + 1
            end do
        end do
        do ip = 1, amr_fw_snp
            sq = f_amr_m1_seq(amr_fw_sprank(ip), 1); tq = f_amr_m1_tag(6, sq)
            call s_xa_rec(XA_F7W_SND, 1, amr_fw_sqsz(ip) - amr_fw_snxp(ip)*XA_NH, tq, peer=amr_fw_sprank(ip), &
                          & key=amr_fw_snxp(ip), seq=sq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = -1
            if (amr_fw_dev) then
                #:call GPU_HOST_DATA(use_device_addr='[amr_fw_sq]')
                    call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                #:endcall GPU_HOST_DATA
            else
                call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
        end do
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call s_wait_tic()
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                call s_wait_toc(WT_RESTR)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "restrict L1 wave: received message length differs from the plan")
                end do
            end block
#else
            call s_wait_tic()
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
            call s_wait_toc(WT_RESTR)
#endif
        end if
        do idx = 1, amr_fw_rnx
            cnt = amr_fw_rpo(idx)
            boff = amr_fw_rqbase(amr_fw_rpi(idx)) + amr_fw_rqo(idx)
            if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(boff + 1:boff + XA_NH), XA_F7W_SND, amr_fw_rblk(idx), amr_fw_rbl(:,idx), &
                & amr_fw_rbh(:,idx))
            ! Device unpack of the covered box only (the strided-update flang trap; see s_restrict_fine_to_coarse)
            call s_l0_pack_unpack_block_sf(coarse_tgt, amr_fw_rbl(1, idx) - o1, amr_fw_rbl(2, idx) - o2, amr_fw_rbl(3, idx) - o3, &
                                           & amr_fw_rbh(1, idx) - amr_fw_rbl(1, idx), amr_fw_rbh(2, idx) - amr_fw_rbl(2, idx), &
                                           & amr_fw_rbh(3, idx) - amr_fw_rbl(3, idx), &
                                           & amr_fw_rq(boff + XA_NH + 1:boff + XA_NH + cnt), .false.)
        end do
#endif

    end subroutine s_amr_restrict_l1_wave

    !> Deliver the current level>=2 block's fine flux registers to its parent block's owner, which holds the matching creg and
    !! applies the correction. One blocking send/recv pair per dimension, mirroring s_amr_p2p_reflux_faces:
    !! freg(d)%lo/hi(:,:,:,slot) is contiguous (trailing slot index fixed, leading dims full), so it goes on the wire with no pack
    !! and its GPU_UPDATE is a contiguous transfer. Several remote children of one parent reuse these tags, which is safe because
    !! MPI does not overtake between a fixed (source, tag, comm) triple and both owners walk the sibling loop in the same replicated
    !! block order. Tag base is disjoint from s_amr_p2p_reflux_faces so an L0/L1 delivery can never be mistaken for a parent
    !! delivery.
    impure subroutine s_amr_p2p_freg_to_parent(pblk)

        integer, intent(in) :: pblk

#ifdef MFC_MPI
        integer :: cowner, powner, cnt, ierr

        cowner = amr_block_owner(amr_cur)
        powner = amr_block_owner(pblk)
        if (proc_rank == cowner) then
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    ! not slot amr_reg_cur: it is 0 on the L0-tiles path
                    cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                    $:GPU_UPDATE(host='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur), freg(' + str(D) &
                                 & + ')%hi(:, :, :, amr_reg_cur)]')
                    call s_xa_rec(XA_F5_FREG_SND, 1, cnt, ${40 + 2*D}$)
                    call MPI_SEND(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, powner, ${40 + 2*D}$, MPI_COMM_WORLD, ierr)
                    call s_xa_rec(XA_F5_FREG_SND, 1, cnt, ${41 + 2*D}$)
                    call MPI_SEND(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, powner, ${41 + 2*D}$, MPI_COMM_WORLD, ierr)
                end if
            #:endfor
        else
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    ! not slot amr_reg_cur: it is 0 on the L0-tiles path
                    cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                    call s_xa_rec(XA_F5_FREG_RCV, 2, cnt, ${40 + 2*D}$)
                    call s_wait_tic()
                    call MPI_RECV(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, cowner, ${40 + 2*D}$, MPI_COMM_WORLD, &
                                  & MPI_STATUS_IGNORE, ierr)
                    call s_xa_rec(XA_F5_FREG_RCV, 2, cnt, ${41 + 2*D}$)
                    call MPI_RECV(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, cowner, ${41 + 2*D}$, MPI_COMM_WORLD, &
                                  & MPI_STATUS_IGNORE, ierr)
                    call s_wait_toc(WT_RESTR)
                    $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur), freg(' + str(D) &
                                 & + ')%hi(:, :, :, amr_reg_cur)]')
                end if
            #:endfor
        end if
#endif

    end subroutine s_amr_p2p_freg_to_parent

    !> Sibling-seam face weights for level>=2 block kb under parent pblk: 0 on a face shared with a same-parent sibling tile
    !! (fine-fine, not c/f: refluxing there double-writes and leaks; the outside parent cell is covered by the sibling's restrict),
    !! 1 otherwise. Replicated metadata only (f_amr_parent_block + f_amr_seam read amr_region_*_all), so every rank derives the same
    !! weights; the reflux apply and the freg wave's wire-skip both call this, so what ships and what is consumed cannot drift
    !! apart.
    impure subroutine s_amr_sibling_face_weights(kb, pblk, w_lo, w_hi)

        integer, intent(in)   :: kb, pblk
        real(wp), intent(out) :: w_lo(3), w_hi(3)
        integer               :: c, y, d

        ! iterate pblk's cached children (same-parent guarantees same level) instead of scanning every block with an
        ! O(global blocks) parent lookup per candidate. Every caller sits in a routine that already refreshed the epoch-keyed
        ! lists; refreshing here would let a loop body reallocate the very list its caller iterates.

        w_lo = 1._wp; w_hi = 1._wp
        if (pblk <= 0) return  ! orphan parent: no siblings (all weights 1); the CSR would index ptr(-1)
        do c = amr_child_ptr(pblk - 1) + 1, amr_child_ptr(pblk)
            y = amr_child_idx(c)
            if (y == kb) cycle
            do d = 1, num_dims
                if (f_amr_seam(kb, y, d)) w_hi(d) = 0._wp  ! sibling just above -> shared high face
                if (f_amr_seam(y, kb, d)) w_lo(d) = 0._wp  ! sibling just below -> shared low face
            end do
        end do

    end subroutine s_amr_sibling_face_weights

    !> Multi-level reflux: apply the Berger-Colella C/F flux correction from the current level>=2 block into its parent block's
    !! cells just outside the block footprint, in the parent-fine frame (mirror of the L0 s_amr_apply_reflux targeted at the parent;
    !! "the coarse" is level l-1). State form: q_parent(outside) += dt*(F_coarse - Fbar_fine)/dxf on the low face and +=
    !! dt*(Fbar_fine - F_coarse)/dxf on the high face, where Fbar_fine is the child-averaged fine register. creg/freg key off this
    !! block's slot. Per-face parent-fine dx (stretched-grid safe).
    !!
    !! The parent's owner applies: it holds the parent field and the parent-side creg (captured over its own advance). Only freg
    !! crosses the wire, and only when the two owners differ. Both participants must reach this routine or the P2P pair deadlocks
    !! (cf. the restrict).
    impure subroutine s_amr_reflux_to_parent(dt_reflux, do_xchg)

        real(wp), intent(in) :: dt_reflux
        !> exchange the split-ownership freg here (the subcycle per-box path); the lock-step driver ships them inside the
        !! restrict-parent wave first
        logical, intent(in) :: do_xchg
        integer             :: pblk, d, olo(3), ohi(3), glo(3), ghi(3), woff(3), plo(3), phi(3)
        real(wp)            :: w_lo(3), w_hi(3), mlo(3), mhi(3)
        logical             :: own_child, own_parent

        call s_amr_refresh_lists()  ! cached parent (f_amr_parent_block is an O(global blocks) scan; this runs per block)
        pblk = amr_parent_blk(amr_cur)
        own_child = amr_rank_owns_block
        own_parent = (amr_block_owner(pblk) == proc_rank)
        if (.not. (own_child .or. own_parent)) return
        if (do_xchg .and. (own_child .neqv. own_parent)) call s_amr_p2p_freg_to_parent(pblk)
        if (.not. own_parent) return
        ! max_grid_size tiling of a level>=2 feature: a face shared with an adjacent sibling tile (same parent) is fine-fine, not
        ! a c/f boundary; its "outside" parent cell is covered by the sibling's restrict, so refluxing there double-writes and
        ! leaks. Skip those faces (weight 0); the fine-fine halo already matched the shared seam flux. No siblings -> all weights
        ! 1 (no-op).
        call s_amr_sibling_face_weights(amr_cur, pblk, w_lo, w_hi)
        ! parent-fine frame for the shared reflux kernel: outside cell = isect boundary +/-1; creg-local loop range 0:extent;
        ! transverse write at the isect origin. Per-face parent-fine cell widths (dx at the low/high outside cell, olo/ohi),
        ! mirroring the L0/L1 s_amr_apply_reflux_state so a stretched parent grid corrects each C/F face with its own width.
        ! Footprint from replicated metadata (s_amr_parent_foot), not amr_isect_lo/hi: on the parent's owner the child's own isect
        ! is
        ! the empty non-owner sentinel whenever the two differ. Identical box while co-located. rr likewise comes from the global
        ! amr_ref_ratio rather than amr_slots(amr_cur), whose slot need not be allocated on this rank.
        call s_amr_parent_foot(amr_cur, pblk, plo, phi)
        olo = 0; ohi = 0; glo = 0; ghi = 0; woff = 0; mlo = 1._wp; mhi = 1._wp
        do d = 1, num_dims
            olo(d) = plo(d) - 1; ohi(d) = phi(d) + 1
            ghi(d) = phi(d) - plo(d)
            woff(d) = plo(d)
        end do
        mlo(1) = amr_slots(pblk)%dx(olo(1)); mhi(1) = amr_slots(pblk)%dx(ohi(1))
        if (n_glb > 0) then; mlo(2) = amr_slots(pblk)%dy(olo(2)); mhi(2) = amr_slots(pblk)%dy(ohi(2)); end if
        if (p_glb > 0) then; mlo(3) = amr_slots(pblk)%dz(olo(3)); mhi(3) = amr_slots(pblk)%dz(ohi(3)); end if
        call s_amr_br_load_faces(amr_loc_of(pblk), olo, ohi, glo, ghi, woff, w_lo, w_hi)
        call s_amr_reflux_apply_faces(amr_cons_br, amr_reg_cur, amr_ref_ratio, dt_reflux, olo, ohi, glo, ghi, woff, w_lo, w_hi, &
                                      & mlo, mhi)
        call s_amr_br_store_faces(amr_loc_of(pblk), olo, ohi, glo, ghi, woff, w_lo, w_hi)

    end subroutine s_amr_reflux_to_parent

    !> Device-native restriction overwrite: restrict the fine block (device) to coarse averages over the covered coarse cells
    !! [bl:bh] global and write coarse_tgt (device) directly, with no host round-trip and only the covered cells touched (a
    !! whole-coarse device push would clobber non-covered cells). Child-sum order: ddk, ddj, then ddi; /nchild; stp cast. The fine
    !! source (the flat store) and coarse_tgt are device-resident. Twin: s_amr_restrict_pack_device runs this same child-sum into a
    !! wire buffer; any change to the loop order, arithmetic, or casts here must be mirrored there byte-identically (owner-local and
    !! scattered coarse cells must match bit-for-bit). Twin (q<->pb/mv) s_amr_restrict_pbmv_box_device runs this same child-sum on
    !! pb/mv; keep the stencil in lockstep. Two targets, one body: the coarse destination is the level-0 monolithic field (`_sf`) or
    !! a parent block in the flat store (`_st`); the fine source is always a block, so it is always the store.
    #:for SFX, CT in [('sf', ''), ('st', 'amr_cons_st')]
        #:set CW = (lambda ix: CT + '(ci - o1, cj - o2, ck - o3, ' + ix + ', ctloc)') if CT else (lambda ix: 'coarse_tgt(' + ix &
                    & + ')%sf(ci - o1, cj - o2, ck - o3)')
        impure subroutine s_amr_restrict_overwrite_device_${SFX}$(${'ctloc' if CT else 'coarse_tgt'}$, loc, bl, bh, o1, o2, o3, &
            & rlo, rr, dj_hi, dk_hi, nchild)

            #:if CT
                integer, intent(in) :: ctloc
            #:else
                type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
            #:endif
            integer, intent(in) :: loc
            integer, intent(in) :: bl(3), bh(3), o1, o2, o3, rlo(3), rr, dj_hi, dk_hi, nchild
            integer             :: i, ci, cj, ck, fi0, fj0, fk0, ddi, ddj, ddk, bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3
            real(wp)            :: acc, wacc, w

            bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
            rl1 = rlo(1); rl2 = rlo(2); rl3 = rlo(3)
            if (cyl_coord) then
                ! axisymmetric volume-weighted fold-back: weight each fine child by its cell-center radius (amr_rvw = fine y_cc, on
                ! device). Same child order as the Cartesian path and the scatter pack, so CPU==GPU and np=1==np>=2.
                $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc, wacc, w]')
                do i = 1, sys_size
                    do ck = bl3, bh3
                        do cj = bl2, bh2
                            do ci = bl1, bh1
                                fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                                acc = 0._wp; wacc = 0._wp
                                do ddk = 0, dk_hi
                                    do ddj = 0, dj_hi
                                        w = amr_rvw(fj0 + ddj)
                                        do ddi = 0, rr - 1
                                            acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, loc), wp)*w
                                            wacc = wacc + w
                                        end do
                                    end do
                                end do
                                ${CW('i')}$ = real(acc/wacc, stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                return
            end if
            $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc]')
            do i = 1, sys_size
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                            acc = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    do ddi = 0, rr - 1
                                        acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, loc), wp)
                                    end do
                                end do
                            end do
                            ${CW('i')}$ = real(acc/real(nchild, wp), stp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_restrict_overwrite_device_${SFX}$
    #:endfor

    !> Device pack of one destination's covered restrict slice (np>=2 scatter): restrict the fine block (device) over the covered
    !! coarse box [bl:bh] global straight into the contiguous wire buffer buf (host, via copyout); only the slice crosses PCIe, not
    !! the full fine field. Same child-sum order and wp values as s_amr_restrict_overwrite_device (no stp cast: the wire carries wp
    !! and the receiver casts), packed with ci fastest, then cj, ck, i, matching the receiver's sequential unpack. Twin:
    !! s_amr_restrict_overwrite_device runs this same child-sum in place; any change to the loop order, arithmetic, or casts here
    !! must be mirrored there byte-identically (owner-local and scattered coarse cells must match bit-for-bit). Twin (q<->pb/mv)
    !! s_amr_restrict_pbmv_pack_device runs this same child-sum into a wire buffer; keep them in lockstep.
    impure subroutine s_amr_restrict_pack_device(loc, bl, bh, rlo, rr, dj_hi, dk_hi, nchild, buf)

        integer, intent(in) :: loc
        integer, intent(in) :: bl(3), bh(3), rlo(3), rr, dj_hi, dk_hi, nchild
        real(wp), intent(inout), contiguous :: buf(:)
        integer :: i, ci, cj, ck, fi0, fj0, fk0, ddi, ddj, ddk, bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3, n1, n2, n3
        real(wp) :: acc, wacc, w

        bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
        rl1 = rlo(1); rl2 = rlo(2); rl3 = rlo(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        if (cyl_coord) then
            ! axisymmetric volume-weighted pack (amr_rvw = fine y_cc, on device); same child order as the overwrite kernel
            $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc, wacc, w]', copyout='[buf]')
            do i = 1, sys_size
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                            acc = 0._wp; wacc = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    w = amr_rvw(fj0 + ddj)
                                    do ddi = 0, rr - 1
                                        acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, loc), wp)*w
                                        wacc = wacc + w
                                    end do
                                end do
                            end do
                            buf(1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*(i - 1)))) = acc/wacc
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            return
        end if
        $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc]', copyout='[buf]')
        do i = 1, sys_size
            do ck = bl3, bh3
                do cj = bl2, bh2
                    do ci = bl1, bh1
                        fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                        acc = 0._wp
                        do ddk = 0, dk_hi
                            do ddj = 0, dj_hi
                                do ddi = 0, rr - 1
                                    acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, loc), wp)
                                end do
                            end do
                        end do
                        buf(1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*(i - 1)))) = acc/real(nchild, wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_restrict_pack_device

end module m_amr_transfer
