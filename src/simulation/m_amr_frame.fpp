!>
!!@file
!!@brief Contains module m_amr_frame

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Block-frame swap (fine grid state in/out of the shared solver) and the pb/mv side-state services.
module m_amr_frame

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

    implicit none

    private
    public :: s_amr_backup_pbmv, s_amr_fine_rk_update_pbmv, s_amr_pressure_relax_fine, s_amr_prolong_pbmv, s_amr_relax_fine, &
        & s_amr_restore_coarse, s_amr_scatter_pbmv, s_amr_swap_to_fine, s_restrict_pbmv

contains

    !> Apply phase-change relaxation (relax) to the current fine block's interior, before restriction. Relaxation is a cell-local,
    !! mass/energy-conserving equilibration (no stencil, no ghosts), so it needs no coarse/fine coupling: it just runs over the fine
    !! interior. Swaps m/n/p to the fine extents so s_infinite_relaxation_k's 0:m,0:n,0:p loop covers this block. Matches the coarse
    !! timing (once per full step; the coarse relax runs once after s_tvd_rk on q_cons_ts(1)) but on the fine solution so the fine
    !! dynamics equilibrate at fine resolution, not only the restricted coarse average.
    impure subroutine s_amr_relax_fine()

        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_amr_br_load(amr_loc_of(amr_cur))
        call s_infinite_relaxation_k(amr_cons_br)
        call s_amr_br_store(amr_loc_of(amr_cur))
        call s_amr_restore_coarse()

    end subroutine s_amr_relax_fine

    !> 6-equation model: apply the per-stage pressure relaxation to the fine block's interior (cell-local equilibration, no
    !! stencil), mirroring the coarse per-stage call. Swaps the grid so the routine's 0:m,0:n,0:p loop covers this block.
    impure subroutine s_amr_pressure_relax_fine()

        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_amr_br_load(amr_loc_of(amr_cur))
        call s_pressure_relaxation_procedure(amr_cons_br)
        call s_amr_br_store(amr_loc_of(amr_cur))
        call s_amr_restore_coarse()

    end subroutine s_amr_pressure_relax_fine

    !> Non-polytropic QBMM: piecewise-constant prolongation of the block's pb/mv interior from the gathered coarse side-state
    !! amr_cg_pb/mv (patch-local frame; the callers run s_amr_gather_coarse_patch_pbmv on all ranks first, so np>=2 couples to the
    !! correct coarse rank). Host loops + device push; the gather is host-current (.not. pull_host). Twin
    !! s_interpolate_coarse_to_fine (pb/mv<->q): pb/mv is piecewise-constant where q_cons is minmod-limited, but the child-offset
    !! frame and realizability/closure policy track it; keep those in lockstep.
    impure subroutine s_amr_prolong_pbmv()

        integer :: fi, fj, fk, q, ib_, ci, cj, ck, rr, lo1, lo2, lo3, ox, oy, oz

        ! Host prolongation (both call paths make the coarse pb/mv host mirrors current first: init writes them on the host,
        ! regrid refreshes them from the device); the device copy of the fine side-state is pushed at the end

        ! coarse pb/mv are read from the gathered block-local patch amr_cg_pb/mv (fine-level distribution): patch-local frame, cell
        ! 0
        ! == amr_cpat_off (matching s_prolong_one_var). The gather is a host loop (.not. pull_host) done by the callers.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rr = amr_slots(amr_cur)%amr_ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = 0, amr_slots(amr_cur)%p
                    ck = 0
                    if (p_glb > 0) ck = lo3 + fk/rr - oz
                    do fj = 0, amr_slots(amr_cur)%n
                        cj = 0
                        if (n_glb > 0) cj = lo2 + fj/rr - oy
                        do fi = 0, amr_slots(amr_cur)%m
                            ci = lo1 + fi/rr - ox
                            amr_slots(amr_cur)%pb_f%sf(fi, fj, fk, q, ib_) = amr_cg_pb(ci, cj, ck, q, ib_)
                            amr_slots(amr_cur)%mv_f%sf(fi, fj, fk, q, ib_) = amr_cg_mv(ci, cj, ck, q, ib_)
                        end do
                    end do
                end do
            end do
        end do
        $:GPU_UPDATE(device='[amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf]')

    end subroutine s_amr_prolong_pbmv

    !> Non-polytropic QBMM: device copy of the block's pb/mv into the step-entry backup (SSP-RK). Twin s_amr_copy_fine_fields
    !! (pb/mv<->q): q_cons sibling of this step-entry backup; keep them in lockstep.
    impure subroutine s_amr_backup_pbmv(pb_s, mv_s, pb_d, mv_d)

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_s, mv_s
        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_d, mv_d
        integer :: fi, fj, fk, q, ib_, b1, e1, b2, e2, b3, e3

        b1 = amr_slots(amr_cur)%idwbuff(1)%beg; e1 = amr_slots(amr_cur)%idwbuff(1)%end
        b2 = amr_slots(amr_cur)%idwbuff(2)%beg; e2 = amr_slots(amr_cur)%idwbuff(2)%end
        b3 = amr_slots(amr_cur)%idwbuff(3)%beg; e3 = amr_slots(amr_cur)%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=5)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = b3, e3
                    do fj = b2, e2
                        do fi = b1, e1
                            pb_d(fi, fj, fk, q, ib_) = pb_s(fi, fj, fk, q, ib_)
                            mv_d(fi, fj, fk, q, ib_) = mv_s(fi, fj, fk, q, ib_)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_backup_pbmv

    !> Non-polytropic QBMM: SSP-RK stage update of the block's pb/mv (device kernel, interior only; mirror of the coarse pb_ts/mv_ts
    !! stage combination in m_time_steppers). Twin s_amr_fine_rk_update + s_tvd_rk (m_time_steppers): same SSP-RK stage combination
    !! on pb/mv; keep all three in lockstep.
    impure subroutine s_amr_fine_rk_update_pbmv(pb_u, mv_u, pb_s, mv_s, rpb, rmv, c1, c2, c3, c4, dtl)

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_u, mv_u
        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_s, mv_s
        real(wp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: rpb, rmv
        real(wp), intent(in) :: c1, c2, c3, c4, dtl
        integer              :: fi, fj, fk, q, ib_, fm, fn, fp

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        $:GPU_PARALLEL_LOOP(collapse=5)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = 0, fp
                    do fj = 0, fn
                        do fi = 0, fm
                            pb_u(fi, fj, fk, q, ib_) = (c1*pb_u(fi, fj, fk, q, ib_) + c2*pb_s(fi, fj, fk, q, &
                                 & ib_) + c3*dtl*rpb(fi, fj, fk, q, ib_))/c4
                            mv_u(fi, fj, fk, q, ib_) = (c1*mv_u(fi, fj, fk, q, ib_) + c2*mv_s(fi, fj, fk, q, &
                                 & ib_) + c3*dtl*rmv(fi, fj, fk, q, ib_))/c4
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_rk_update_pbmv

    !> Non-polytropic QBMM: volume-weighted restriction of the block's pb/mv onto the coarse side-state under the block (device
    !! kernel; same equal-weight child average as the q_cons restrict).
    impure subroutine s_restrict_pbmv(pb_c, mv_c, pb_fin, mv_fin)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_c, mv_c

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_fin, mv_fin
        integer  :: ci, cj, ck, q, ib_, fi0, fj0, fk0, ddi, ddj, ddk, nchild, ox, oy, oz, rr
        integer  :: c1lo, c1hi, c2lo, c2hi, c3lo, c3hi, dj_hi, dk_hi
        real(wp) :: accp, accm

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        rr = amr_slots(amr_cur)%amr_ref_ratio
        nchild = rr
        if (n_glb > 0) nchild = nchild*rr
        if (p_glb > 0) nchild = nchild*rr
        c1lo = amr_isect_lo(1); c1hi = amr_isect_hi(1)
        c2lo = amr_isect_lo(2); c2hi = merge(amr_isect_hi(2), amr_isect_lo(2), n_glb > 0)
        c3lo = amr_isect_lo(3); c3hi = merge(amr_isect_hi(3), amr_isect_lo(3), p_glb > 0)
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        $:GPU_PARALLEL_LOOP(collapse=5, private='[fi0, fj0, fk0, ddi, accp, accm, ddj, ddk]')
        do ib_ = 1, nb
            do q = 1, nnode
                do ck = c3lo, c3hi
                    do cj = c2lo, c2hi
                        do ci = c1lo, c1hi
                            fi0 = (ci - c1lo)*rr; fj0 = (cj - c2lo)*rr; fk0 = (ck - c3lo)*rr
                            accp = 0._wp; accm = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    do ddi = 0, rr - 1
                                        accp = accp + real(pb_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                        accm = accm + real(mv_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                    end do
                                end do
                            end do
                            pb_c(ci - ox, cj - oy, ck - oz, q, ib_) = accp/real(nchild, wp)
                            mv_c(ci - ox, cj - oy, ck - oz, q, ib_) = accm/real(nchild, wp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_restrict_pbmv

    !> Non-polytropic QBMM np>=2: device restriction of pb/mv over the covered box [bl:bh] global into pb_c/mv_c (device), fine
    !! origin (ci-rlo)*rr, local coarse index cell - o. Only the covered cells the owner holds are touched (no whole-array push, the
    !! same GPU-only clobber the q_cons device overwrite avoids). Same child-sum as s_restrict_pbmv. Twin:
    !! s_amr_restrict_pbmv_pack_device runs this same child-sum into a wire buffer; any change to the loop order, arithmetic, or
    !! casts here must be mirrored there byte-identically (local and scattered coarse pb/mv must match bit-for-bit). Twin
    !! (pb/mv<->q) s_amr_restrict_overwrite_device is the q_cons sibling of this child-sum; keep the stencil in lockstep.
    impure subroutine s_amr_restrict_pbmv_box_device(pb_c, mv_c, pb_fin, mv_fin, bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_c, mv_c

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_fin, mv_fin
        integer, intent(in) :: bl(3), bh(3), o1, o2, o3, rlo(3), rr, dj_hi, dk_hi, nchild
        integer             :: ci, cj, ck, q, ib_, fi0, fj0, fk0, ddi, ddj, ddk, bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3
        real(wp)            :: accp, accm

        bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
        rl1 = rlo(1); rl2 = rlo(2); rl3 = rlo(3)
        $:GPU_PARALLEL_LOOP(collapse=5, private='[fi0, fj0, fk0, ddi, accp, accm, ddj, ddk]')
        do ib_ = 1, nb
            do q = 1, nnode
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                            accp = 0._wp; accm = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    do ddi = 0, rr - 1
                                        accp = accp + real(pb_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                        accm = accm + real(mv_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                    end do
                                end do
                            end do
                            pb_c(ci - o1, cj - o2, ck - o3, q, ib_) = real(accp/real(nchild, wp), stp)
                            mv_c(ci - o1, cj - o2, ck - o3, q, ib_) = real(accm/real(nchild, wp), stp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_restrict_pbmv_box_device

    !> Non-polytropic QBMM np>=2 scatter pack: device restriction of pb/mv over the covered box [bl:bh] global straight into the
    !! contiguous wire buffer buf (host, via copyout; pb block then mv block, ci fastest); only the slice crosses PCIe, not the full
    !! fine side-state. Same child-sum as s_amr_restrict_pbmv_box_device; the wire carries wp and the receiver casts to stp. Twin:
    !! s_amr_restrict_pbmv_box_device runs this same child-sum in place; any change to the loop order, arithmetic, or casts here
    !! must be mirrored there byte-identically (local and scattered coarse pb/mv must match bit-for-bit). Twin (pb/mv<->q)
    !! s_amr_restrict_pack_device is the q_cons sibling of this packed child-sum; keep them in lockstep.
    impure subroutine s_amr_restrict_pbmv_pack_device(pb_fin, mv_fin, bl, bh, rlo, rr, dj_hi, dk_hi, nchild, buf)

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_fin, mv_fin
        integer, intent(in) :: bl(3), bh(3), rlo(3), rr, dj_hi, dk_hi, nchild
        real(wp), intent(inout), contiguous :: buf(:)
        integer :: ci, cj, ck, q, ib_, fi0, fj0, fk0, ddi, ddj, ddk, bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3
        integer :: n1, n2, n3, half
        real(wp) :: accp, accm

        bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
        rl1 = rlo(1); rl2 = rlo(2); rl3 = rlo(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        half = n1*n2*n3*nnode*nb
        $:GPU_PARALLEL_LOOP(collapse=5, private='[fi0, fj0, fk0, ddi, ddj, ddk, accp, accm]', copyout='[buf]')
        do ib_ = 1, nb
            do q = 1, nnode
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                            accp = 0._wp; accm = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    do ddi = 0, rr - 1
                                        accp = accp + real(pb_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                        accm = accm + real(mv_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                    end do
                                end do
                            end do
                            buf(1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*((q - 1) + nnode*(ib_ - 1))))) &
                                & = accp/real(nchild, wp)
                            buf(half + 1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*((q - 1) + nnode*(ib_ - 1))))) &
                                & = accm/real(nchild, wp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_restrict_pbmv_pack_device

    !> Device unpack of a received pb/mv covered slice into the coarse fields: exact inverse of s_amr_restrict_pbmv_pack_device's
    !! wire layout (ci fastest, then cj, ck, q, ib_, all of pb followed by all of mv). Unpacking on the device is required, not a
    !! convenience: a host unpack plus a strided GPU_UPDATE(device=) of the covered box is miscopied as a contiguous run; see the
    !! note at the q_cons unpack in s_restrict_fine_to_coarse, of which this is the pb/mv twin.
    impure subroutine s_amr_restrict_pbmv_unpack_device(pb_c, mv_c, bl, bh, o1, o2, o3, buf)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_c, mv_c
        integer, intent(in) :: bl(3), bh(3), o1, o2, o3
        real(wp), intent(inout), contiguous :: buf(:)
        integer :: ci, cj, ck, q, ib_, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3, half

        bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        half = n1*n2*n3*nnode*nb
        $:GPU_PARALLEL_LOOP(collapse=5, copyin='[buf]')
        do ib_ = 1, nb
            do q = 1, nnode
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            pb_c(ci - o1, cj - o2, ck - o3, q, &
                                 & ib_) = real(buf(1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*((q - 1) + nnode*(ib_ &
                                 & - 1))))), stp)
                            mv_c(ci - o1, cj - o2, ck - o3, q, &
                                 & ib_) = real(buf(half + 1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*((q - 1) &
                                 & + nnode*(ib_ - 1))))), stp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_restrict_pbmv_unpack_device

    !> Non-polytropic QBMM: distributed fine->coarse fold-back of the block's pb/mv onto the coarse side-state pb_ts/mv_ts,
    !! mirroring the q_cons scatter in s_restrict_fine_to_coarse. The owner restricts the covered cells it holds (device, owned box)
    !! and sends each other coarse-owner its covered pb/mv slice (device-packed averages, pb block then mv block); that owner
    !! overwrites its local coarse. Called on all ranks at np>=2 (owner/non-owner split inside) so the P2P send/recv pair up; np=1
    !! is handled locally by the direct s_restrict_pbmv call (this routine is reached only when num_procs > 1). Twin
    !! s_restrict_fine_to_coarse (scatter half, pb/mv<->q): same scatter skeleton (owner sends covered coarse slices, each
    !! coarse-owner overwrites its local cells); keep them in lockstep.
    impure subroutine s_amr_scatter_pbmv(pb_fin, mv_fin)

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_fin, mv_fin
        integer               :: nchild, rr, dj_hi, dk_hi, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr, cellsz
        integer               :: rlo(3), rhi(3), ilo(3), ihi(3), bl(3), bh(3)
        real(wp), allocatable :: sbuf(:,:), rbuf(:)
        integer, allocatable  :: reqs(:), drank(:)

        rr = amr_slots(amr_cur)%amr_ref_ratio
        nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        cellsz = 2*nnode*nb
        rlo = 0; rhi = 0
        rlo(1) = amr_region_lo_all(1, amr_cur); rhi(1) = amr_region_hi_all(1, amr_cur)
        if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, amr_cur); rhi(2) = amr_region_hi_all(2, amr_cur); end if
        if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, amr_cur); rhi(3) = amr_region_hi_all(3, amr_cur); end if
        owner = amr_block_owner(amr_cur)
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        maxsz = cellsz*(rhi(1) - rlo(1) + 1)*(rhi(2) - rlo(2) + 1)*(rhi(3) - rlo(3) + 1)

        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()

        if (proc_rank == owner) then
            ! overwrite the covered cells this rank owns (device, owned box), then send each other coarse-owner its covered slice
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_pbmv_box_device(pb_ts(1)%sf, &
                & mv_ts(1)%sf, pb_fin, mv_fin, bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
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
                    boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    ! pack this destination's covered pb/mv slice on the device (restrict averages straight into the wire
                    ! buffer, same child-sum as the device overwrite above), with no full-field fine host pull
                    call s_amr_restrict_pbmv_pack_device(pb_fin, mv_fin, bl, bh, rlo, rr, dj_hi, dk_hi, nchild, sbuf(1:boxsz,nsrc))
#ifdef MFC_MPI
                    call s_xa_rec(XA_F7C_SND, 1, boxsz, amr_cur)
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
            ! coarse-owner: if I hold covered cells, receive my pb/mv slice from the owner and overwrite my local coarse
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                allocate (rbuf(boxsz))
#ifdef MFC_MPI
                call s_xa_rec(XA_F7C_RCV, 2, boxsz, amr_cur)
                call s_wait_tic()
                call MPI_RECV(rbuf, boxsz, mpi_p, owner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                call s_wait_toc(WT_RESTR)
#endif
                ! Device unpack, writing only the covered cells (a whole-array push would clobber device-advanced non-covered
                ! coarse cells with this rank's stale host copy)
                call s_amr_restrict_pbmv_unpack_device(pb_ts(1)%sf, mv_ts(1)%sf, bl, bh, o1, o2, o3, rbuf)
                deallocate (rbuf)
            end if
        end if

    end subroutine s_amr_scatter_pbmv

    !> Swap the global grid state to the fine block. Must be paired with s_amr_restore_coarse.
    impure subroutine s_amr_swap_to_fine()

        ! Saving on a nested swap would overwrite the sw_* bounce buffers with fine state, and the eventual restore would install
        ! fine extents as the coarse grid (silent corruption of everything after). Hence every save below is depth-guarded; the
        ! installs are not, since re-installing the same slot is idempotent.
        amr_swap_depth = amr_swap_depth + 1
        if (amr_swap_depth == 1) then
            sw_m = m; sw_n = n; sw_p = p
            sw_idwint = idwint; sw_idwbuff = idwbuff
        end if
        ! the acoustic source's precomputed spatials are coarse-grid cell indices: applying them on the fine block would inject at
        ! wrong cells (or out of bounds). The support is guaranteed not to overlap the block (checked at startup), so the fine RHS
        ! correctly skips the source.
        if (amr_swap_depth == 1) sw_acoustic_source = acoustic_source
        acoustic_source = .false.
        ! active-box windows are coarse cell indices: applying them on the swapped fine grid would window the wrong cells. Blocks
        ! are
        ! contained in the active window (init check + regrid clamp), so the fine advance legitimately treats its whole block as
        ! active.
        if (amr_swap_depth == 1) sw_ab_active = ab_active
        ab_active = .false.
        $:GPU_UPDATE(device='[ab_active]')
        m = amr_slots(amr_cur)%m; n = amr_slots(amr_cur)%n; p = amr_slots(amr_cur)%p
        idwint(1)%beg = 0; idwint(1)%end = m
        idwint(2)%beg = 0; idwint(2)%end = n
        idwint(3)%beg = 0; idwint(3)%end = p
        idwbuff = amr_slots(amr_cur)%idwbuff
        ! save coarse coords to bounce buffers, then copy fine coords into global arrays
        if (amr_swap_depth == 1) then
            sw_x_cb = x_cb; sw_x_cc = x_cc; sw_dx = dx
            if (n_glb > 0) then; sw_y_cb = y_cb; sw_y_cc = y_cc; sw_dy = dy; end if
            if (p_glb > 0) then; sw_z_cb = z_cb; sw_z_cc = z_cc; sw_dz = dz; end if
        end if
        x_cb(-1:amr_slots(amr_cur)%m) = amr_slots(amr_cur)%x_cb(-1:amr_slots(amr_cur)%m)
        x_cc(0:amr_slots(amr_cur)%m) = amr_slots(amr_cur)%x_cc(0:amr_slots(amr_cur)%m)
        dx(0:amr_slots(amr_cur)%m) = amr_slots(amr_cur)%dx(0:amr_slots(amr_cur)%m)
        if (n_glb > 0) then
            y_cb(-1:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%y_cb(-1:amr_slots(amr_cur)%n)
            y_cc(0:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%y_cc(0:amr_slots(amr_cur)%n)
            dy(0:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%dy(0:amr_slots(amr_cur)%n)
        end if
        if (p_glb > 0) then
            z_cb(-1:amr_slots(amr_cur)%p) = amr_slots(amr_cur)%z_cb(-1:amr_slots(amr_cur)%p)
            z_cc(0:amr_slots(amr_cur)%p) = amr_slots(amr_cur)%z_cc(0:amr_slots(amr_cur)%p)
            dz(0:amr_slots(amr_cur)%p) = amr_slots(amr_cur)%dz(0:amr_slots(amr_cur)%p)
        end if
        ! Extend the fine grid into the ghost shell (s_build_level_coords only fills the interior 0:m). Ghost cells use the exact
        ! parent-cell bisection, the same formula as the interior, with floor division for negative indices. Fine-level
        ! distribution: the owner may not hold the block's coarse coordinate slice locally, so ghost parent boundaries come from the
        ! global boundaries amr_g?cb (cl is a global coarse index, region_lo + floor(jg/rr)), matching the interior build. Blocks
        ! stay buff_size inside the domain, so every ghost parent is an in-domain coarse cell with exact coords.
        block
            integer               :: jg, cl, pblk2, k, rr, pnf
            real(wp), allocatable :: cxb(:), cyb(:), czb(:), tcc(:), tdx(:)
            rr = amr_slots(amr_cur)%amr_ref_ratio
            ! ghost parent boundaries: a level>=2 block's coarse side is its parent's fine grid (indexed in the parent-fine
            ! amr_isect frame, matching the interior s_build_level_coords), not the L0 global boundaries. amr_isect_lo is a
            ! parent-fine index, so indexing amr_g?cb (sized for L0) would read out of bounds (garbage on host, NaN on the device
            ! copy). Source the parent's fine coords for level>=2, the global L0 boundaries for level 1.
            if (amr_block_level(amr_cur) >= 2) then
                ! Rebuild the parent's fine boundaries from replicated metadata; do not read amr_slots(pblk2)%x_cb. That array is
                ! allocated only on the parent's owner, and under per-level distribution this block's owner need not be it; taking
                ! lbound/ubound of an unallocated allocatable is undefined. Same ancestor replay as the interior build, so the
                ! ghost bisection and the interior agree exactly.
                pblk2 = f_amr_parent_block(amr_cur)
                pnf = amr_ref_ratio**amr_block_level(pblk2)*(amr_region_hi_all(1, pblk2) - amr_region_lo_all(1, pblk2) + 1) - 1
                allocate (cxb(-1:pnf), tcc(0:pnf), tdx(0:pnf))
                call s_amr_build_block_coords(pblk2, amr_gxcb, cxb, tcc, tdx, 1)
                deallocate (tcc, tdx)
                if (n_glb > 0) then
                    pnf = amr_ref_ratio**amr_block_level(pblk2)*(amr_region_hi_all(2, pblk2) - amr_region_lo_all(2, pblk2) + 1) - 1
                    allocate (cyb(-1:pnf), tcc(0:pnf), tdx(0:pnf))
                    call s_amr_build_block_coords(pblk2, amr_gycb, cyb, tcc, tdx, 2)
                    deallocate (tcc, tdx)
                end if
                if (p_glb > 0) then
                    pnf = amr_ref_ratio**amr_block_level(pblk2)*(amr_region_hi_all(3, pblk2) - amr_region_lo_all(3, pblk2) + 1) - 1
                    allocate (czb(-1:pnf), tcc(0:pnf), tdx(0:pnf))
                    call s_amr_build_block_coords(pblk2, amr_gzcb, czb, tcc, tdx, 3)
                    deallocate (tcc, tdx)
                end if
            else
                allocate (cxb(lbound(amr_gxcb, 1):ubound(amr_gxcb, 1))); cxb = amr_gxcb
                if (n_glb > 0) then; allocate (cyb(lbound(amr_gycb, 1):ubound(amr_gycb, 1))); cyb = amr_gycb; end if
                if (p_glb > 0) then; allocate (czb(lbound(amr_gzcb, 1):ubound(amr_gzcb, 1))); czb = amr_gzcb; end if
            end if
            do jg = amr_slots(amr_cur)%m + 1, amr_slots(amr_cur)%m + buff_size
                cl = amr_isect_lo(1) + floor(real(jg, wp)/real(rr, wp))
                k = modulo(jg, rr)
                if (k == rr - 1) then
                    x_cb(jg) = cxb(cl)
                else
                    x_cb(jg) = (real(rr - 1 - k, wp)*cxb(cl - 1) + real(k + 1, wp)*cxb(cl))/real(rr, wp)
                end if
                dx(jg) = x_cb(jg) - x_cb(jg - 1); x_cc(jg) = 0.5_wp*(x_cb(jg - 1) + x_cb(jg))
            end do
            ! unified boundary formula (matches the interior subdivision): boundary jg belongs to
            ! parent c = isect_lo + floor(jg/rr); sub-position k=modulo(jg,rr) picks the rr-way split
            do jg = -1 - buff_size, -1
                cl = amr_isect_lo(1) + floor(real(jg, wp)/real(rr, wp))
                k = modulo(jg, rr)
                if (k == rr - 1) then
                    x_cb(jg) = cxb(cl)
                else
                    x_cb(jg) = (real(rr - 1 - k, wp)*cxb(cl - 1) + real(k + 1, wp)*cxb(cl))/real(rr, wp)
                end if
            end do
            do jg = -buff_size, -1
                dx(jg) = x_cb(jg) - x_cb(jg - 1); x_cc(jg) = 0.5_wp*(x_cb(jg - 1) + x_cb(jg))
            end do
            if (n_glb > 0) then
                do jg = amr_slots(amr_cur)%n + 1, amr_slots(amr_cur)%n + buff_size
                    cl = amr_isect_lo(2) + floor(real(jg, wp)/real(rr, wp))
                    k = modulo(jg, rr)
                    if (k == rr - 1) then
                        y_cb(jg) = cyb(cl)
                    else
                        y_cb(jg) = (real(rr - 1 - k, wp)*cyb(cl - 1) + real(k + 1, wp)*cyb(cl))/real(rr, wp)
                    end if
                    dy(jg) = y_cb(jg) - y_cb(jg - 1); y_cc(jg) = 0.5_wp*(y_cb(jg - 1) + y_cb(jg))
                end do
                ! unified boundary formula (matches the interior subdivision): boundary jg belongs to
                ! parent c = isect_lo + floor(jg/rr); sub-position k=modulo(jg,rr) picks the rr-way split
                do jg = -1 - buff_size, -1
                    cl = amr_isect_lo(2) + floor(real(jg, wp)/real(rr, wp))
                    k = modulo(jg, rr)
                    if (k == rr - 1) then
                        y_cb(jg) = cyb(cl)
                    else
                        y_cb(jg) = (real(rr - 1 - k, wp)*cyb(cl - 1) + real(k + 1, wp)*cyb(cl))/real(rr, wp)
                    end if
                end do
                do jg = -buff_size, -1
                    dy(jg) = y_cb(jg) - y_cb(jg - 1); y_cc(jg) = 0.5_wp*(y_cb(jg - 1) + y_cb(jg))
                end do
            end if
            if (p_glb > 0) then
                do jg = amr_slots(amr_cur)%p + 1, amr_slots(amr_cur)%p + buff_size
                    cl = amr_isect_lo(3) + floor(real(jg, wp)/real(rr, wp))
                    k = modulo(jg, rr)
                    if (k == rr - 1) then
                        z_cb(jg) = czb(cl)
                    else
                        z_cb(jg) = (real(rr - 1 - k, wp)*czb(cl - 1) + real(k + 1, wp)*czb(cl))/real(rr, wp)
                    end if
                    dz(jg) = z_cb(jg) - z_cb(jg - 1); z_cc(jg) = 0.5_wp*(z_cb(jg - 1) + z_cb(jg))
                end do
                ! unified boundary formula (matches the interior subdivision): boundary jg belongs to
                ! parent c = isect_lo + floor(jg/rr); sub-position k=modulo(jg,rr) picks the rr-way split
                do jg = -1 - buff_size, -1
                    cl = amr_isect_lo(3) + floor(real(jg, wp)/real(rr, wp))
                    k = modulo(jg, rr)
                    if (k == rr - 1) then
                        z_cb(jg) = czb(cl)
                    else
                        z_cb(jg) = (real(rr - 1 - k, wp)*czb(cl - 1) + real(k + 1, wp)*czb(cl))/real(rr, wp)
                    end if
                end do
                do jg = -buff_size, -1
                    dz(jg) = z_cb(jg) - z_cb(jg - 1); z_cc(jg) = 0.5_wp*(z_cb(jg - 1) + z_cb(jg))
                end do
            end if
        end block
        ! batched advance: the leader's grid is installed above; extend it into the slab of amr_bat_n stacked blocks (stride
        ! amr_bat_w along amr_bat_sd), since the flux divergence reads dx/dy/dz at every slab cell. Cell boundaries (x_cb etc.)
        ! are not replicated: nothing on the batched path reads them (WENO coefficients are not recomputed on a uniform grid).
        if (amr_bat_n > 1) then
            block
                integer :: ibm, o, e
                e = amr_bat_ext(amr_bat_sd)
                do ibm = 2, amr_bat_n
                    o = (ibm - 1)*amr_bat_w
                    select case (amr_bat_sd)
                    case (1)
                        x_cc(o - buff_size:o + e + buff_size) = x_cc(-buff_size:e + buff_size)
                        dx(o - buff_size:o + e + buff_size) = dx(-buff_size:e + buff_size)
                    case (2)
                        y_cc(o - buff_size:o + e + buff_size) = y_cc(-buff_size:e + buff_size)
                        dy(o - buff_size:o + e + buff_size) = dy(-buff_size:e + buff_size)
                    case (3)
                        z_cc(o - buff_size:o + e + buff_size) = z_cc(-buff_size:e + buff_size)
                        dz(o - buff_size:o + e + buff_size) = dz(-buff_size:e + buff_size)
                    end select
                end do
                e = (amr_bat_n - 1)*amr_bat_w + e
                select case (amr_bat_sd)
                case (1); m = e
                case (2); n = e
                case (3); p = e
                end select
                idwint(amr_bat_sd)%end = e; idwbuff(amr_bat_sd)%end = e + buff_size
            end block
        end if
        ! sync the swapped extents/bounds/coordinates to the device: RHS kernels read the device copies of these GPU_DECLARE'd
        ! globals (stale coarse bounds = OOB kernels)
        call s_amr_sync_grid_state_to_device()
        ! hypoelastic stress sources use grid-spacing-dependent FD coefficients: recompute them from the (now fine) grid, else every
        ! fine velocity gradient is halved
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        ! nonuniform coarse grid (stretched, or the axisymmetric axis half-cell): the per-cell WENO coefficients must be rebuilt for
        ! the block's own grid (no-op flag on uniform grids)
        if (amr_weno_coef_recompute) call s_amr_recompute_weno_coefs()

        ! IGR: save the coarse sigma state and seed the fine solve. jac holds this stage's converged coarse sigma (the coarse RHS
        ! ran
        ! first), so its parent values are both the best initial guess and the frozen Dirichlet ghost data for the block-local
        ! Jacobi
        ! solve (the per-iteration BC/halo populate is skipped under amr_in_fine_advance). Piecewise-constant parent injection over
        ! the full buffered fine range.
        if (igr) call s_amr_igr_swap_sigma()

    end subroutine s_amr_swap_to_fine

    !> Restore the global grid state saved by s_amr_swap_to_fine.
    impure subroutine s_amr_restore_coarse(sync_device)

        !> .false. only from the batched stage when another batch follows at once: its swap re-pushes the whole grid state before
        !! any kernel reads it, so the restore-side push is dead work there.
        logical, intent(in), optional :: sync_device

        @:ASSERT(amr_swap_depth > 0, "s_amr_restore_coarse without a matching s_amr_swap_to_fine")
        amr_swap_depth = amr_swap_depth - 1
        ! inner restore: the enclosing swap frame still wants its slot installed, and the sw_* buffers still hold the coarse state
        if (amr_swap_depth > 0) return
        m = sw_m; n = sw_n; p = sw_p
        idwint = sw_idwint; idwbuff = sw_idwbuff
        acoustic_source = sw_acoustic_source
        ab_active = sw_ab_active
        $:GPU_UPDATE(device='[ab_active]')
        ! restore full coarse coords from bounce buffers
        x_cb = sw_x_cb; x_cc = sw_x_cc; dx = sw_dx
        if (n_glb > 0) then; y_cb = sw_y_cb; y_cc = sw_y_cc; dy = sw_dy; end if
        if (p_glb > 0) then; z_cb = sw_z_cb; z_cc = sw_z_cc; dz = sw_dz; end if
        ! sync the restored coarse extents/bounds/coordinates back to the device
        if (.not. present(sync_device)) then
            call s_amr_sync_grid_state_to_device()
        else if (sync_device) then
            call s_amr_sync_grid_state_to_device()
        end if
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        if (amr_weno_coef_recompute) call s_amr_recompute_weno_coefs()
        if (igr) call s_amr_igr_restore_sigma()

    end subroutine s_amr_restore_coarse

    !> Save the coarse jac/jac_old and seed the (already swapped-in) fine block's sigma state by piecewise-constant parent injection
    !! from the saved coarse sigma: interior = initial guess, ghost shell = frozen Dirichlet coupling data for the block-local
    !! solve.
    impure subroutine s_amr_igr_swap_sigma()

        integer :: j, k, l, ci, cj, ck, ibm, g, o1, o2, o3, mb1, me1, mb2, me2, mb3, me3
        integer :: cb1, ce1, cb2, ce2, cb3, ce3, fb1, fe1, fb2, fe2, fb3, fe3
        integer :: lo1, lo2, lo3, ox, oy, oz

        ! bounds/offsets hoisted to scalars: sw_idwbuff (and friends) are host-only module state; referencing them inside the
        ! kernels makes OpenACC's present lookup fail (OpenMP's implicit map(to) tolerates it, so only OpenACC builds crash)

        cb1 = sw_idwbuff(1)%beg; ce1 = sw_idwbuff(1)%end
        cb2 = sw_idwbuff(2)%beg; ce2 = sw_idwbuff(2)%end
        cb3 = sw_idwbuff(3)%beg; ce3 = sw_idwbuff(3)%end
        fb1 = idwbuff(1)%beg; fe1 = idwbuff(1)%end
        fb2 = idwbuff(2)%beg; fe2 = idwbuff(2)%end
        fb3 = idwbuff(3)%beg; fe3 = idwbuff(3)%end
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        ! Save the coarse sigma, outermost swap only. A nested swap must not re-save, or sw_jac would take fine state and both the
        ! seed below and s_amr_igr_restore_sigma would work from it. The seed that follows is not guarded: it reads sw_jac, which
        ! still holds the coarse state, so every nested block seeds from the correct parent.
        if (amr_swap_depth == 1) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = cb3, ce3
                do k = cb2, ce2
                    do j = cb1, ce1
                        sw_jac(j, k, l) = jac(j, k, l)
                        sw_jac_old(j, k, l) = jac_old(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        if (amr_bat_n > 1) then
            ! batched slab: each member's buffered range is seeded from its own parent; the slab bounds above are the leader's
            do ibm = 1, amr_bat_n
                g = amr_bat_blk(ibm)
                lo1 = amr_isect_lo_all(1, g); lo2 = amr_isect_lo_all(2, g); lo3 = amr_isect_lo_all(3, g)
                o1 = 0; o2 = 0; o3 = 0
                select case (amr_bat_sd)
                case (1); o1 = (ibm - 1)*amr_bat_w
                case (2); o2 = (ibm - 1)*amr_bat_w
                case default; o3 = (ibm - 1)*amr_bat_w
                end select
                mb1 = -buff_size; me1 = amr_bat_mext(1, ibm) + buff_size
                mb2 = 0; me2 = 0; mb3 = 0; me3 = 0
                if (n_glb > 0) then; mb2 = -buff_size; me2 = amr_bat_mext(2, ibm) + buff_size; end if
                if (p_glb > 0) then; mb3 = -buff_size; me3 = amr_bat_mext(3, ibm) + buff_size; end if
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, ci, cj, ck]', copyin='[lo1, lo2, lo3, o1, o2, o3, mb1, me1, &
                                    & mb2, me2, mb3, me3]')
                do l = mb3, me3
                    do k = mb2, me2
                        do j = mb1, me1
                            ci = lo1 + floor(real(j, wp)/real(amr_ref_ratio, wp)) - ox
                            cj = 0; ck = 0
                            if (n_glb > 0) cj = lo2 + floor(real(k, wp)/real(amr_ref_ratio, wp)) - oy
                            if (p_glb > 0) ck = lo3 + floor(real(l, wp)/real(amr_ref_ratio, wp)) - oz
                            ci = min(max(ci, cb1), ce1)
                            cj = min(max(cj, cb2), ce2)
                            ck = min(max(ck, cb3), ce3)
                            jac(j + o1, k + o2, l + o3) = sw_jac(ci, cj, ck)
                            jac_old(j + o1, k + o2, l + o3) = sw_jac(ci, cj, ck)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            return
        end if
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, ci, cj, ck]')
        do l = fb3, fe3
            do k = fb2, fe2
                do j = fb1, fe1
                    ci = lo1 + floor(real(j, wp)/real(amr_ref_ratio, wp)) - ox
                    cj = 0; ck = 0
                    if (n_glb > 0) cj = lo2 + floor(real(k, wp)/real(amr_ref_ratio, wp)) - oy
                    if (p_glb > 0) ck = lo3 + floor(real(l, wp)/real(amr_ref_ratio, wp)) - oz
                    ci = min(max(ci, cb1), ce1)
                    cj = min(max(cj, cb2), ce2)
                    ck = min(max(ck, cb3), ce3)
                    jac(j, k, l) = sw_jac(ci, cj, ck)
                    jac_old(j, k, l) = sw_jac(ci, cj, ck)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_igr_swap_sigma

    !> Restore the coarse jac/jac_old saved by s_amr_igr_swap_sigma (bounds already restored).
    impure subroutine s_amr_igr_restore_sigma()

        integer :: j, k, l, b1, e1, b2, e2, b3, e3

        b1 = idwbuff(1)%beg; e1 = idwbuff(1)%end
        b2 = idwbuff(2)%beg; e2 = idwbuff(2)%end
        b3 = idwbuff(3)%beg; e3 = idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = b3, e3
            do k = b2, e2
                do j = b1, e1
                    jac(j, k, l) = sw_jac(j, k, l)
                    jac_old(j, k, l) = sw_jac_old(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_igr_restore_sigma

end module m_amr_frame
