!>
!!@file
!!@brief Contains module m_amr_restart

#:include 'macros.fpp'

!> @brief AMR fine-level restart I/O: writes/reads the fine-level restart file alongside the level-0 restart (serial per-rank
!! unformatted files, or one shared MPI-IO file under parallel_io). Split out of m_amr; block/slot state stays in m_amr (and
!! m_global_parameters).
module m_amr_restart

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field
    use m_global_parameters
    use m_constants, only: amr_restart_blk_hdr_ints, amr_restart_blk_own_ints
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_min
    use m_amr_restart_io
    use m_amr, only: s_amr_reduce_xchg_flag, amr_slots, amr_cons_st, amr_loc_of, amr_seam_pairs_dirty, amr_mesh_epoch, &
        & s_amr_alloc_slot, s_amr_reconcile_slots, s_amr_assign_block_owners, s_set_amr_fine_geometry
    use m_amr_regrid, only: s_amr_check_seam_topology

    implicit none

    private
    public :: s_write_amr_restart, s_read_amr_restart

contains

    !> Write the fine-level restart file for save step t_step alongside the level-0 restart (whose format stays untouched): the
    !! writing rank count, the active-block count, and for each block its box + each rank's intersection-local fine conservative
    !! state. Serial mode: one unformatted file per rank inside its level-0 step directory. Parallel mode: one shared MPI-IO file
    !! (3-int global header [-np, nboxes, sys_size], then per block a 7-int box+level header, a 4-int ownership record [owner+1, m,
    !! n, p; validated on read], then the ranks' fine blocks concatenated in rank order). At the same rank count the decomposition
    !! must match (enforced by the ownership record); a different rank count repartitions on read.
    impure subroutine s_write_amr_restart(t_step)

        integer, intent(in)                  :: t_step
        character(LEN=path_len + 3*name_len) :: file_loc
        integer                              :: i, k, wext(3)

#ifdef MFC_MPI
        integer :: ifile, ierr, cnt, idx, fi, fj, fk, bhdr(amr_restart_blk_hdr_ints), ibytes, sbytes
        integer, allocatable :: myown_all(:)
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: my_off, disp0, ddisp
        integer(kind=MPI_OFFSET_KIND), allocatable :: my_cnt_vec(:), my_off_vec(:), tot_cnt_vec(:)
        logical :: file_exist
        real(stp), allocatable :: buf(:)
#endif

        if (.not. amr) return
        ! host consumer: fine state is device-current during stepping (pull every owned slot)
        do k = 1, amr_num_blocks
            if (amr_owns_all(k)) then
                $:GPU_UPDATE(host='[amr_cons_st(:, :, :, :, amr_loc_of(k))]')
            end if
        end do

        if (.not. parallel_io) then
            ! per-rank file in the step directory just created by the level-0 serial write
            file_loc = f_amr_restart_path(t_step)
            open (2, FILE=trim(file_loc), form='unformatted', STATUS='new')
            write (2) num_procs, amr_num_blocks, sys_size
            do k = 1, amr_num_blocks
                ! per-block header: region box, refinement level (a level-l block's fine extent is amr_ref_ratio**l, not
                ! amr_ref_ratio, of the region - the reader needs the level to rebuild multi-level geometry), extents
                ! -1 extents when this rank does not own the slot: that is the reader's data-presence flag, and a slot this
                ! rank never built (an L0 tile owned elsewhere) has no extents of its own to report.
                wext = -1
                if (amr_owns_all(k)) wext = [amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p]
                write (2) amr_region_lo_all(:,k), amr_region_hi_all(:,k), amr_block_level(k), wext
                if (amr_owns_all(k)) then
                    do i = 1, sys_size
                        write (2) amr_cons_st(0:amr_slots(k)%m,0:amr_slots(k)%n,0:amr_slots(k)%p,i, amr_loc_of(k))
                    end do
                end if
            end do
            close (2)
        else
#ifdef MFC_MPI
            ibytes = storage_size(0)/8; sbytes = storage_size(0._stp)/8
            file_loc = f_amr_restart_path(t_step)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
            ! MPI-IO file handles default to MPI_ERRORS_RETURN: failures silent unless checked
            if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart write: MPI_FILE_OPEN failed for ' // trim(file_loc))
            ! Format v2: a negative rank count marks the compact per-block ownership record (see m_constants). A v1 reader
            ! sees a rank mismatch and aborts with its existing message rather than misparsing the block records.
            if (proc_rank == 0) call MPI_FILE_WRITE_AT(ifile, int(0, MPI_OFFSET_KIND), [-num_procs, amr_num_blocks, sys_size], 3, &
                & MPI_INTEGER, status, ierr)
            disp0 = int(3*ibytes, MPI_OFFSET_KIND)  ! running byte offset past the 3-int global header
            ! hoist per-block metadata collectives: one EXSCAN/ALLREDUCE/ALLGATHER over all blocks
            allocate (my_cnt_vec(amr_num_blocks), my_off_vec(amr_num_blocks), tot_cnt_vec(amr_num_blocks))
            allocate (myown_all(amr_restart_blk_own_ints*amr_num_blocks))
            do k = 1, amr_num_blocks
                cnt = sys_size*(amr_slots(k)%m + 1)*(amr_slots(k)%n + 1)*(amr_slots(k)%p + 1)
                if (.not. amr_owns_all(k)) cnt = 0
                my_cnt_vec(k) = int(cnt, MPI_OFFSET_KIND)
                myown_all(amr_restart_blk_own_ints*(k - 1) + 1:amr_restart_blk_own_ints*k) = 0
                if (amr_owns_all(k)) myown_all(amr_restart_blk_own_ints*(k - 1) + 1:amr_restart_blk_own_ints*k) = [proc_rank + 1, &
                    & amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p]
            end do
            my_off_vec = int(0, MPI_OFFSET_KIND)
            call MPI_EXSCAN(my_cnt_vec, my_off_vec, amr_num_blocks, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) my_off_vec = int(0, MPI_OFFSET_KIND)
            call MPI_ALLREDUCE(my_cnt_vec, tot_cnt_vec, amr_num_blocks, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
            ! Per-block (owner + 1, m, n, p): every rank contributes zeros except the one owner, so a MAX reduction recovers the
            ! record exactly. It carries the same information as the v1 3*num_procs extent vector (readers rebuild it from
            ! their own decomposition and abort on mismatch, catching a different ownership pattern or load_balance split that
            ! would otherwise silently misalign the concatenated per-rank data slices below) in O(blocks) rather than
            ! O(blocks x ranks), in memory and in the file.
            call MPI_ALLREDUCE(MPI_IN_PLACE, myown_all, amr_restart_blk_own_ints*amr_num_blocks, MPI_INTEGER, MPI_MAX, &
                               & MPI_COMM_WORLD, ierr)
            do k = 1, amr_num_blocks
                cnt = int(my_cnt_vec(k), kind(cnt))
                my_off = my_off_vec(k)
                if (proc_rank == 0) then
                    ! amr_restart_blk_hdr_ints-int per-block header: region box (6) + refinement level (a level-l block's fine
                    ! extent is amr_ref_ratio**l, not amr_ref_ratio, of the region - the reader needs the level to rebuild
                    ! multi-level geometry). Header layout is single-sourced in m_constants so both readers stay in lockstep.
                    bhdr(1:3) = amr_region_lo_all(:,k); bhdr(4:6) = amr_region_hi_all(:,k)
                    bhdr(amr_restart_blk_hdr_ints) = amr_block_level(k)
                    call MPI_FILE_WRITE_AT(ifile, disp0, bhdr, amr_restart_blk_hdr_ints, MPI_INTEGER, status, ierr)
                end if
                if (proc_rank == 0) then
                    call MPI_FILE_WRITE_AT(ifile, disp0 + int(amr_restart_blk_hdr_ints*ibytes, MPI_OFFSET_KIND), &
                                           & myown_all(amr_restart_blk_own_ints*(k - 1) + 1), amr_restart_blk_own_ints, &
                                           & MPI_INTEGER, status, ierr)
                end if
                ddisp = disp0 + int((amr_restart_blk_hdr_ints + amr_restart_blk_own_ints)*ibytes, MPI_OFFSET_KIND)
                allocate (buf(max(cnt, 1)))
                idx = 0
                ! cnt == 0 on a non-owning rank, where buf is the 1-element placeholder and this slot's q_cons is not allocated
                ! (lazy owned-only sizing) while its m/n/p metadata IS replicated - so packing would run the full extent loop over
                ! an unallocated slot and overrun buf. The collective WRITE_AT_ALL below still runs on every rank, with cnt = 0.
                if (cnt > 0) then
                    do i = 1, sys_size
                        do fk = 0, amr_slots(k)%p
                            do fj = 0, amr_slots(k)%n
                                do fi = 0, amr_slots(k)%m
                                    idx = idx + 1
                                    buf(idx) = amr_cons_st(fi, fj, fk, i, amr_loc_of(k))
                                end do
                            end do
                        end do
                    end do
                end if
                call MPI_FILE_WRITE_AT_ALL(ifile, ddisp + my_off*int(sbytes, MPI_OFFSET_KIND), buf, cnt*mpi_io_type, mpi_io_p, &
                                           & status, ierr)
                if (ierr /= MPI_SUCCESS) &
                    & call s_mpi_abort('amr restart write: data write failed (disk full/quota?); the file is unusable')
                deallocate (buf)
                disp0 = ddisp + tot_cnt_vec(k)*int(sbytes, MPI_OFFSET_KIND)
            end do
            deallocate (my_cnt_vec, my_off_vec, tot_cnt_vec, myown_all)
            ! close is where buffered MPI-IO data flushes on many stacks - a failure here truncates the file
            call MPI_FILE_CLOSE(ifile, ierr)
            if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart write: MPI_FILE_CLOSE failed; the file may be truncated')
#endif
        end if

    end subroutine s_write_amr_restart

    !> Restore the fine level from the AMR restart file at t_step_start (n_start under cfl_dt), if one exists: for each saved block
    !! rebuild the box via s_set_amr_fine_geometry, then read each rank's intersection-local fine state (exact stp round-trip).
    !! parallel_io repartitions across rank counts (each block is one contiguous region-sized chunk under whole-block ownership,
    !! re-assigned to this run's owners); serial (per-rank files) needs the writing rank count. restored = false on a fresh start,
    !! or (with a one-line warning) on a restart without the file; the caller then re-prolongs from coarse. Collective: all ranks
    !! call it together.
    impure subroutine s_read_amr_restart(restored)

        logical, intent(out)                 :: restored
        character(LEN=path_len + 3*name_len) :: file_loc
        logical                              :: file_exist
        integer                              :: ts, have_loc, have_glb

        restored = .false.
        if (.not. amr) return
        if (cfl_dt) then
            ts = n_start
        else
            ts = t_step_start
        end if
        if (ts == 0) return  ! fresh start: the fine level is prolonged from the pre_process ICs

        file_loc = f_amr_restart_path(ts)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        have_loc = merge(1, 0, file_exist)
        call s_mpi_allreduce_integer_min(have_loc, have_glb)
        if (have_glb == 0) then
            if (proc_rank == 0) then
                print '(A)', &
                    & ' [amr] WARNING: no AMR restart file at this step; the fine level is re-initialized by ' &
                    & // 'prolongation from coarse (fine-level accuracy is lost across this restart)'
            end if
            return
        end if

        if (parallel_io) then
            call s_amr_restart_read_parallel(file_loc)
        else
            call s_amr_restart_read_serial(file_loc)
        end if

        call s_amr_select_slot(1)
        ! restored levels without a regrid: the per-level fill waves iterate 2..amr_num_levels, so leaving it at the
        ! default 1 would silently skip every level>=2 fill until the first regrid recomputes it
        amr_num_levels = max(1, maxval(amr_block_level(1:amr_num_blocks)))
        amr_seam_pairs_dirty = .true.  ! restored a new block set: the cached seam-pair list must be rebuilt
        amr_mesh_epoch = amr_mesh_epoch + 1
        call s_amr_check_seam_topology()  ! abort on seam topologies no halo reconciles (e.g. restart mode-switch)
        restored = .true.

    end subroutine s_read_amr_restart

    !> Install the restored block set: rebuild whole-block owners from the regions, allocate this run's owned slots (freeing any
    !! stale init slots), then each block's geometry under the correct owner.
    impure subroutine s_amr_restart_install_blocks()

        integer :: k

        call s_amr_assign_block_owners()
        call s_amr_reconcile_slots()
        do k = l0_slot_off + 1, amr_num_blocks  ! fine blocks only: an L0 tile is not refined, and s_l0_tiles_init built it
            amr_cur = k
            call s_set_amr_fine_geometry(amr_region_lo_all(:,k), amr_region_hi_all(:,k))
        end do
        call s_amr_reduce_xchg_flag()

    end subroutine s_amr_restart_install_blocks

    !> Serial (per-rank files) read: every block's region, level and, where this rank owned it at write (rm >= 0), its fine state.
    !! Whole-block ownership is decomposition-deterministic, so the file's data-presence flag drives the read and is verified
    !! against this run's ownership once the owners are rebuilt.
    impure subroutine s_amr_restart_read_serial(file_loc)

        character(LEN=*), intent(in) :: file_loc
        integer                      :: i, k, ghdr(3), reg(6), lvl, rm, rn, rp, nblk
        logical, allocatable         :: had_data(:)

        open (2, FILE=trim(file_loc), form='unformatted', ACTION='read', STATUS='old')
        read (2) ghdr
        call s_amr_restart_check_header(ghdr, sys_size, 'amr restart', nblk)
        call s_amr_restart_check_blocks(nblk)
        amr_num_blocks = nblk
        allocate (had_data(amr_num_blocks))
        do k = 1, amr_num_blocks
            read (2) reg, lvl, rm, rn, rp
            call s_amr_restart_check_record(k, reg, lvl)
            had_data(k) = rm >= 0
            if (.not. had_data(k)) cycle
            ! An L0 tile's record: the tile already exists (s_l0_tiles_init) with this rank's own state, and its data here is a
            ! copy of the level-0 field the ordinary restart file restores, which s_l0_copy_coarse_to_tiles reseeds the tiles
            ! from. Consume the record - the stream is positional - and leave the tile alone.
            if (lvl == 0) then
                do i = 1, sys_size
                    read (2) amr_cons_st(0:rm,0:rn,0:rp,i, amr_loc_of(k))
                end do
                cycle
            end if
            ! whole-block owner extents are region-derived per level (a level-l block covers amr_ref_ratio**l fine cells per L0
            ! cell of its region); a stored extent that disagrees is corrupt
            if (rm /= (amr_ref_ratio**lvl)*(reg(4) - reg(1) + 1) - 1 .or. rn /= merge((amr_ref_ratio**lvl)*(reg(5) - reg(2) + 1) &
                & - 1, 0, n_glb > 0) .or. rp /= merge((amr_ref_ratio**lvl)*(reg(6) - reg(3) + 1) - 1, 0, p_glb > 0)) then
                call s_mpi_abort('amr restart: block fine extents disagree with the region (corrupt file)')
            end if
            call s_amr_alloc_slot(k)
            ! zero the whole host column first: the read fills only the interior, but the push covers the full padded column,
            ! and the store grows device-side so the host pad bytes are otherwise undefined
            amr_cons_st(:,:,:,:,amr_loc_of(k)) = 0._stp
            do i = 1, sys_size
                read (2) amr_cons_st(0:rm,0:rn,0:rp,i, amr_loc_of(k))
            end do
            ! Push this block before the next s_amr_alloc_slot: allocating a slot can grow the shared flat store, and
            ! s_amr_st_reserve preserves the growth by pulling device->host first, which would overwrite the block just read with
            ! the not-yet-written device copy (restored fine state becomes NaN). The store is device-authoritative at every alloc
            ! point; a host writer must close that gap itself.
            $:GPU_UPDATE(device='[amr_cons_st(:, :, :, :, amr_loc_of(k))]')
        end do
        close (2)
        call s_amr_restart_install_blocks()
        do k = 1, amr_num_blocks
            if (had_data(k) .neqv. amr_owns_all(k)) then
                call s_mpi_abort('amr restart decomposition mismatch: the file''s block ownership differs from this' &
                                 & // ' run''s (identical decomposition - rank count and load_balance settings - required)')
            end if
        end do
        deallocate (had_data)

    end subroutine s_amr_restart_read_serial

    !> Parallel (MPI-IO) read. The writer's rank count sets only the file layout: whole-block ownership makes each block's fine data
    !! one contiguous region-sized chunk, so any rank count can read it; owners are re-assigned for this run and each new owner
    !! reads its whole blocks. At the writer's rank count the file's ownership record must reproduce this run's decomposition (a
    !! re-derived load_balance split would silently misalign the chunk).
    impure subroutine s_amr_restart_read_parallel(file_loc)

        character(LEN=*), intent(in) :: file_loc

#ifdef MFC_MPI
        integer                             :: i, k, ghdr(3), nblk, ifile, ierr, cnt, idx, fi, fj, fk
        integer                             :: mown(amr_restart_blk_own_ints)
        type(t_amr_restart_catalog)         :: cat
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND)       :: fsz
        real(stp), allocatable              :: buf(:)

        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
        if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart read: MPI_FILE_OPEN failed for ' // trim(file_loc))
        ! MPI-IO errors are silent by default and a read past EOF returns short with an uninitialized tail: the file size is
        ! compared against the catalog's layout below, so a truncated file fails closed instead of restoring garbage
        call MPI_FILE_GET_SIZE(ifile, fsz, ierr)
        call MPI_FILE_READ_AT_ALL(ifile, int(0, MPI_OFFSET_KIND), ghdr, 3, MPI_INTEGER, status, ierr)
        call s_amr_restart_check_header(ghdr, sys_size, 'amr restart', nblk)
        call s_amr_restart_check_blocks(nblk)
        amr_num_blocks = nblk
        call s_amr_restart_read_catalog(ifile, amr_num_blocks, sys_size, cat)
        if (cat%end_disp /= fsz) call s_mpi_abort('amr restart read: file size does not match the expected layout ' &
            & // '(truncated or corrupt amr restart file)')
        do k = 1, amr_num_blocks
            call s_amr_restart_check_record(k, cat%reg(:,k), cat%lvl(k))
        end do
        call s_amr_restart_install_blocks()
        do k = 1, amr_num_blocks
            cnt = 0
            if (amr_owns_all(k)) then
                cnt = sys_size*(amr_slots(k)%m + 1)*(amr_slots(k)%n + 1)*(amr_slots(k)%p + 1)
                mown = [proc_rank + 1, amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p]
                if (abs(ghdr(1)) == num_procs .and. any(cat%own(:,k) /= mown)) then
                    call s_mpi_abort('amr restart: the per-block owner/extent record in the file does not match ' &
                                     & // 'this run''s decomposition; with the same rank count the ownership and ' &
                                     & // '(with load_balance) the weighted splits must match the run that wrote the restart')
                end if
            end if
            allocate (buf(max(cnt, 1)))
            ! collective: every rank calls it, non-owners with count 0
            call MPI_FILE_READ_AT_ALL(ifile, cat%data_disp(k), buf, cnt*mpi_io_type, mpi_io_p, status, ierr)
            ! owner only: zero the whole host column (the unpack fills only the interior, the push covers the padded column)
            if (cnt > 0) amr_cons_st(:,:,:,:,amr_loc_of(k)) = 0._stp
            idx = 0
            do i = 1, sys_size
                do fk = 0, amr_slots(k)%p
                    do fj = 0, amr_slots(k)%n
                        do fi = 0, amr_slots(k)%m
                            idx = idx + 1
                            amr_cons_st(fi, fj, fk, i, amr_loc_of(k)) = buf(idx)
                        end do
                    end do
                end do
            end do
            deallocate (buf)
        end do
        call MPI_FILE_CLOSE(ifile, ierr)
        ! push the restored fine state to the device. Here (not in s_read_amr_restart, which runs after both readers): the
        ! serial reader pushes each block as it reads it, because allocating the next slot can grow the store, and a later
        ! whole-set push would copy back host columns that the device-side growth and the reconcile's compaction have left
        ! undefined - silently corrupting the state this routine just restored (invisible on a CPU build).
        do k = 1, amr_num_blocks
            if (amr_owns_all(k)) then
                $:GPU_UPDATE(device='[amr_cons_st(:, :, :, :, amr_loc_of(k))]')
            end if
        end do
#endif

    end subroutine s_amr_restart_read_parallel

    !> The file's block count must fit this run's metadata.
    impure subroutine s_amr_restart_check_blocks(nblk)

        integer, intent(in) :: nblk

        if (nblk < 1 .or. nblk > amr_max_blocks) then
            call s_mpi_abort('amr restart: the file holds more fine blocks than amr_max_blocks ' &
                             & // 'in this run; restart with amr_max_blocks at least the written block count')
        end if

    end subroutine s_amr_restart_check_blocks

    !> Corrupt/foreign-file guard for block k's record: a box outside the global domain or a level outside 1..amr_max_level would
    !! drive the geometry build and coordinate reads out of bounds silently in release builds. Stores the accepted region and level
    !! (set before the owner/geometry rebuild: s_amr_assign_block_owners and s_set_amr_fine_geometry key off amr_block_level to
    !! place L>=2 blocks under their parent).
    impure subroutine s_amr_restart_check_record(k, reg, lvl)

        integer, intent(in) :: k, reg(6), lvl

        if (reg(1) < 0 .or. reg(4) > m_glb .or. reg(1) > reg(4) .or. (n_glb > 0 .and. (reg(2) < 0 .or. reg(5) > n_glb .or. reg(2) &
            & > reg(5))) .or. (p_glb > 0 .and. (reg(3) < 0 .or. reg(6) > p_glb .or. reg(3) > reg(6)))) then
            call s_mpi_abort('amr restart: corrupt block record (box outside the global domain)')
        end if
        ! level 0 is an L0 TILE. Under coexist the tiles are the pool's fixed prefix [1, l0_slot_off] and amr_num_blocks counts
        ! them, so the writer emits them here; they are rebuilt by s_l0_tiles_init before this read, from the same l0_ntile and
        ! the same decomposition, so the record must agree with the tile that is already there rather than redefine it. A
        ! level-0 record outside the prefix, or a tile whose box moved, means the restart does not belong to this case.
        if (lvl == 0) then
            if (k > l0_slot_off) call s_mpi_abort('amr restart: corrupt block record (level-0 record outside the L0 tile prefix)')
            if (any(amr_region_lo_all(:,k) /= reg(1:3)) .or. any(amr_region_hi_all(:,k) /= reg(4:6))) then
                call s_mpi_abort('amr restart: the L0 tile layout differs from this run''s (identical l0_ntile and rank count ' &
                                 & // 'required to restart a coexist case)')
            end if
            return
        end if
        if (lvl < 1 .or. lvl > amr_max_level) then
            call s_mpi_abort('amr restart: corrupt block record (block level outside 1..amr_max_level)')
        end if
        amr_region_lo_all(:,k) = reg(1:3); amr_region_hi_all(:,k) = reg(4:6)
        amr_block_level(k) = lvl

    end subroutine s_amr_restart_check_record

end module m_amr_restart
