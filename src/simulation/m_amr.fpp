!>
!!@file
!!@brief Contains module m_amr

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A kernel naming an unallocated
#! module array aborts, so conditionally allocated ones (sw_jac/jac: igr) are named only under their allocation condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Block-structured AMR: up to amr_max_blocks refined blocks (2:1 or 4:1 per amr_ref_ratio), optionally nested to
!! amr_max_level, advanced with the shared solver via grid-state swap and conservatively coupled to each block's parent level (ghost
!! prolongation, Berger-Colella flux reflux, restriction) and dynamic regrid.
module m_amr

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_min, s_mpi_allreduce_integer_max
    use m_ibm, only: s_ibm_alloc_fine
    use m_amr_state
    use m_amr_distribution
    use m_amr_wave
    use m_amr_store
    use m_amr_exchange
    use m_amr_seam
    use m_amr_frame
    use m_amr_transfer
    use m_amr_advance

    implicit none

    !> Re-exports every AMR service module.

contains

    !> Build the static refined level-1 block (no-op unless amr); runs after the level-0 grid and time-steppers are set up.
    impure subroutine s_initialize_amr_module()

        amr_dim = [.true., n_glb > 0, p_glb > 0]
        amr_sidx = 0; amr_sidx(1:num_dims) = start_idx
        amr_ext = merge([m, n, p], 0, amr_dim)

        if (.not. amr) return

#ifdef MFC_GPU
        amr_fw_dev = rdma_mpi
#endif
        amr_br_batch = amr_bat_lim

        call s_amr_alloc_pool()
        amr_block_level = 1  ! init default (level-1); regrid re-tags each block's level for nesting
        amr_num_levels = 1
        amr_num_blocks = 1
        amr_cur = 1
        amr_seam_pairs_dirty = .true.; amr_seam_pairs_nblk = -1  ! force a seam-list build on the first fine-fine halo
        amr_mesh_epoch = amr_mesh_epoch + 1

        call s_amr_init_extents()
        call s_amr_init_report()
        call s_amr_scr_init()
        call s_amr_init_swap_buffers()
        call s_amr_build_global_cb(0)  ! the fine-distribution owner rebuilds whole-block fine coordinates from these
        call s_amr_init_coarse_patch()
        ! the coarse decomposition is a structured split computed O(1) per rank by s_amr_rank_decomp; validate it against this rank
        call s_amr_validate_decomp()
        ! per-slot fine-grid IB marker fields, sized like q_cons
        if (ib) call s_ibm_alloc_fine(amr_max_blocks, mbuf_lo(1), mbuf_hi(1), mbuf_lo(2), mbuf_hi(2), mbuf_lo(3), mbuf_hi(3))
        call s_amr_init_first_blocks()
        call s_amr_init_tags()

    end subroutine s_initialize_amr_module

    !> Block caps and preallocation extents. buff_size is unknown at checker time, so the geometric aborts live here.
    impure subroutine s_amr_init_extents()

        integer :: d, bad_loc, bad_glb, fit_d

        call s_amr_compute_isect(amr_block_beg, amr_block_end)

        ! the fine ghost shell and reflux outside cells must stay inside the global domain
        if (amr_block_beg(1) < buff_size .or. amr_block_end(1) > m_glb - buff_size .or. (n_glb > 0 .and. (amr_block_beg(2) &
            & < buff_size .or. amr_block_end(2) > n_glb - buff_size)) .or. (p_glb > 0 .and. (amr_block_beg(3) < buff_size &
            & .or. amr_block_end(3) > p_glb - buff_size))) then
            call s_mpi_abort('amr block must lie at least buff_size cells inside the domain boundaries')
        end if

        ! The fine advance reuses solver scratch sized to the local grid, so an IB block (never tiled) must fit a rank's local
        ! half-extent; non-IB blocks are tiled to <= amr_maxc_fit.
        bad_loc = 0
        if (ib) then
            if (amr_ref_ratio*(amr_block_end(1) - amr_block_beg(1) + 1) - 1 > m) bad_loc = 1
            if (n_glb > 0 .and. amr_ref_ratio*(amr_block_end(2) - amr_block_beg(2) + 1) - 1 > n) bad_loc = 1
            if (p_glb > 0 .and. amr_ref_ratio*(amr_block_end(3) - amr_block_beg(3) + 1) - 1 > p) bad_loc = 1
        end if
        call s_mpi_allreduce_integer_max(bad_loc, bad_glb)
        if (bad_glb == 1) then
            call s_mpi_abort('amr fine extent exceeds a rank local grid (solver scratch is local-sized): an immersed-body block ' &
                             & // 'is owned whole and un-tiled, so it may cover at most about half of any rank subdomain per ' &
                             & // 'dimension; shrink the body region or use fewer ranks')
        end if

        ! max coarse block cells per dim; 1 for collapsed dims
        amr_maxc = merge(([m_glb, n_glb, p_glb] + 1)/amr_ref_ratio, 1, amr_dim)

        ! Regrid box cap: min over ranks of the local half-extent, or the absolute amr_max_grid_size (rank-count-invariant box set,
        ! scratch sized to the cap). Every box is clamped to amr_maxc_fit, which sizes the fine/coord arrays.
        amr_maxc_fit = amr_maxc
        do d = 1, num_dims
            call s_mpi_allreduce_integer_min((amr_ext(d) + 1)/amr_ref_ratio, fit_d)
            if (amr_max_grid_size > 0) then
                amr_maxc_fit(d) = min(amr_maxc(d), amr_max_grid_size)
            else
                amr_maxc_fit(d) = min(amr_maxc(d), fit_d)
            end if
        end do

        ! max fine extents and buffered bounds
        max_f = merge(amr_ref_ratio*amr_maxc_fit - 1, 0, amr_dim)
        call s_amr_set_mbuf()

    end subroutine s_amr_init_extents

    !> Rank-0 advisories and the memory demand. Collective (reduces over ranks).
    impure subroutine s_amr_init_report()

        integer  :: nonuni, nonuni_glb
        real(wp) :: slot_gib, cells, nfam

        ! stacked members reuse the batch leader's coordinates, exact only on bitwise-uniform spacing; warn if it differs at
        ! roundoff

        nonuni = 0
        if (any(dx(0:m) /= dx(0))) nonuni = 1
        if (n_glb > 0) then; if (any(dy(0:n) /= dy(0))) nonuni = 1; end if
        if (p_glb > 0) then; if (any(dz(0:p) /= dz(0))) nonuni = 1; end if
        call s_mpi_allreduce_integer_max(nonuni, nonuni_glb)
        if (proc_rank /= 0) return

        if (nonuni_glb == 1) print '(A)', &
            & ' [amr] NOTE: the grid''s cell spacing is not bitwise uniform: stacked blocks reuse the batch leader''s ' &
            & // 'coordinate arrays, so members after the first see the leader''s spacing'
        ! the SFC map spreads whole blocks, so fewer blocks than ranks leaves some ranks without fine work
        if (num_procs > amr_max_blocks) print '(A,I0,A,I0,A)', ' [amr] WARNING: amr_max_blocks (', amr_max_blocks, &
            & ') < num_procs (', num_procs, &
            & '): the fine level can occupy at most amr_max_blocks ranks - raise amr_max_blocks for better fine-level balance'
        ! fine blocks advance at the coarse dt; the finest-cell CFL is unknown at init, so warn rather than abort
        if (amr_ref_ratio > 2 .or. amr_max_level > 1) print '(A,I0,A)', &
            & ' [amr] WARNING: fine blocks advance at the coarse dt, but the ' &
            & // 'finest cell is amr_ref_ratio**amr_max_level = ', amr_ref_ratio**amr_max_level, &
            & 'x smaller - ensure dt satisfies the FINEST cell CFL (roughly the coarse-stable dt divided by that ' &
            & // 'factor), else the fine block may go unstable'
        ! amr_regrid_int = 0 is static AMR and never adapts
        if (amr_regrid_int == 0) print '(A)', &
            & ' [amr] NOTE: amr_regrid_int = 0 - the block set is STATIC and never adapts. ' &
            & // 'Set amr_regrid_int > 0 (4-8 is a reasonable start) for adaptive refinement.'
        if (amr_max_grid_size == 0 .and. num_procs > 1) print '(A)', &
            & ' [amr] NOTE: amr_max_grid_size = 0 derives the block cap from the ' &
            & // 'decomposition, so it SHRINKS as ranks are added and the box set depends on rank ' &
            & // 'count. Pin it to make the box set rank-invariant.'
        if (amr_regrid_int > 0 .and. amr_regrid_int < 4) print '(A,I0,A)', ' [amr] NOTE: amr_regrid_int = ', amr_regrid_int, &
            & ' regrids often; the tag sweep is per-CELL and flat in box count. ' &
            & // 'Raise it unless the refined feature moves quickly.'

        ! Memory demand: a block costs 2 field families (q_cons, q_cons_stor; q_prim/rhs are pooled) x sys_size arrays on the mbuf
        ! extents. Exceeding device memory presents as a hang (one rank dies in __tgt_target_data_begin_mapper, the rest block in
        ! MPI).
        cells = real(product(mbuf_hi - mbuf_lo + 1), wp)
        nfam = 2._wp
        slot_gib = cells*real(sys_size, wp)*nfam*real(storage_size(1._wp)/8, wp)/1024._wp**3
        print '(A,I0,A,I0,A,ES10.3,A,F8.3,A)', ' [amr] per-block slot: ', nint(cells), ' cells x sys_size x ', nint(nfam), &
            & ' fields = ', cells*real(sys_size, wp)*nfam, ' words (', slot_gib, ' GiB per owned block)'
        print '(A,F9.2,A,I0,A)', ' [amr]   worst case if one rank owned every block: ', slot_gib*real(amr_max_blocks, wp), &
            & ' GiB (amr_max_blocks = ', amr_max_blocks, '). Typical is amr_max_blocks/num_procs blocks per rank.'

    end subroutine s_amr_init_report

    !> Coarse-patch gather buffer (see amr_cg): the largest block's coarse footprint + 2*amr_cpat_mar halo, device-mapped.
    impure subroutine s_amr_init_coarse_patch()

        type(scalar_field), allocatable :: tmp_cg(:)
        integer                         :: i

        amr_cpat_mar = f_amr_cpat_mar()
        amr_cpat_hi = merge(amr_maxc_fit - 1 + 2*amr_cpat_mar, 0, amr_dim)
        ! CCE OpenMP-offload: allocating a bare module-scope scalar_field array aborts (lib-4425); allocate a local, then
        ! move_alloc.
        allocate (tmp_cg(1:sys_size))
        @:ALLOCATE(amr_slab_tab(1:15, 1:64))
        call move_alloc(tmp_cg, amr_cg)
        $:GPU_ENTER_DATA(create='[amr_cg]')
        do i = 1, sys_size
            @:ALLOCATE(amr_cg(i)%sf(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3)))
            amr_cg(i)%sf = 0._stp  ! padding beyond a block's patch is never read; keep it finite for the device copy
            @:ACC_SETUP_SFs(amr_cg(i))
        end do

    end subroutine s_amr_init_coarse_patch

    !> Place the initial block(s) and set their geometry. With dynamic regrid and bodies the block gets the same body-containment
    !! expansion regrid boxes get; non-IB blocks are tiled to <= amr_maxc_fit, IB keeps one block. Slot arrays are allocated lazily
    !! by s_amr_reconcile_slots.
    impure subroutine s_amr_init_first_blocks()

        type(t_box), allocatable :: tiled(:)
        integer                  :: blk_lo(3), blk_hi(3), nt, capt, kk

        blk_lo = amr_block_beg; blk_hi = amr_block_end
        if (ib .and. amr_regrid_int > 0) call s_amr_expand_box_over_bodies(blk_lo, blk_hi)
        allocate (tiled(amr_max_blocks)); nt = 0; capt = 0
        if (ib) then
            nt = 1; tiled(1)%lo = blk_lo; tiled(1)%hi = blk_hi
        else
            call s_amr_tile_box(blk_lo, blk_hi, tiled, nt, amr_max_blocks, capt)
        end if
        amr_num_blocks = nt
        ! set regions first: the owner assignment reads amr_region_*_all and s_set_amr_fine_geometry needs the owners
        do kk = 1, nt
            amr_region_lo_all(:,kk) = tiled(kk)%lo; amr_region_hi_all(:,kk) = tiled(kk)%hi
        end do
        call s_amr_assign_block_owners()
        call s_amr_reconcile_slots()  ! allocate this rank's owned blocks
        do kk = 1, nt
            amr_cur = kk
            call s_set_amr_fine_geometry(tiled(kk)%lo, tiled(kk)%hi)
        end do
        call s_amr_reduce_xchg_flag()
        call s_amr_select_slot(1)  ! refresh the per-block mirrors
        deallocate (tiled)

    end subroutine s_amr_init_first_blocks

    !> Tag bases: narrow-reduction tags sit above the per-box tags; keyed bands start at the next 65536 boundary above them.
    impure subroutine s_amr_init_tags()

#ifdef MFC_MPI
        integer(kind=MPI_ADDRESS_KIND) :: tag_ub
        logical                        :: tag_ub_set
        integer                        :: ierr
#endif

        amr_tag_narrow = amr_max_blocks + 400
        amr_tag_base = ((amr_max_blocks + 800)/65536 + 1)*65536
#ifdef MFC_MPI
        call MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, tag_ub_set, ierr)
        @:ASSERT(tag_ub_set, "MPI_TAG_UB attribute unavailable")
        @:ASSERT(amr_max_blocks + 800 <= tag_ub, &
                 & "AMR tag space exceeds MPI_TAG_UB: amr_max_blocks is too large for this MPI's tag range")
        @:ASSERT(amr_tag_base + amr_n_bands*65536 <= tag_ub, "AMR keyed-tag band space exceeds MPI_TAG_UB")
#endif

    end subroutine s_amr_init_tags

    impure subroutine s_finalize_amr_module()

        integer :: i, islot

        if (.not. amr) return
        do islot = 1, amr_max_blocks
            call s_amr_free_slot(islot)
        end do
        call s_amr_free_pool()
        ! gfortran/ifx abort on deallocating an unallocated array: guard each
        if (allocated(amr_fw_rblk)) deallocate (amr_fw_rblk)
        if (allocated(amr_sw_sq)) deallocate (amr_sw_sq, amr_sw_rq)
        #:for A in ['amr_fw_sq', 'amr_fw_rq']
            if (allocated(${A}$)) then
                if (amr_fw_dev) then
                    $:GPU_EXIT_DATA(delete='[' + A + ']')
                end if
                deallocate (${A}$)
            end if
        #:endfor
        #:for A in ['amr_my_blk', 'amr_l1r_blk', 'amr_l1p_blk', 'amr_fch_blk', 'amr_own_blk', 'amr_parent_blk', &
            'amr_child_ptr', 'amr_child_idx', 'amr_cg_pz']
            if (allocated(${A}$)) deallocate (${A}$)
        #:endfor
        do i = 1, sys_size
            @:DEALLOCATE(amr_cg(i)%sf)
        end do
        @:DEALLOCATE(amr_cg)
        @:DEALLOCATE(amr_slab_tab)
        if (allocated(amr_pc)) then
            @:DEALLOCATE(amr_pc)
        end if
        call s_amr_free_swap_buffers()

    end subroutine s_finalize_amr_module

end module m_amr
