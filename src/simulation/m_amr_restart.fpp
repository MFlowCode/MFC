!>
!!@file
!!@brief Contains module m_amr_restart

#:include 'macros.fpp'

!> @brief AMR fine-level restart I/O: writes/reads the fine-level restart file alongside the level-0 restart (serial per-rank
!! unformatted files, or one shared MPI-IO file under parallel_io).
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
    use m_amr_state, only: max_f
    use m_amr, only: s_amr_reduce_xchg_flag, amr_slots, amr_cons_st, amr_loc_of, amr_seam_pairs_dirty, amr_mesh_epoch, &
        & s_amr_alloc_slot, s_amr_reconcile_slots, s_amr_assign_block_owners, s_set_amr_fine_geometry
    use m_amr_regrid, only: s_amr_check_seam_topology
    use m_amr_distribution, only: s_amr_adopt_block_owners

    implicit none

    private
    public :: s_write_amr_restart, s_read_amr_restart

contains

    !> Serial: one unformatted file per rank in the level-0 step directory just created by the serial write.
    impure subroutine s_amr_restart_write_serial(t_step)

        integer, intent(in)                  :: t_step
        character(LEN=path_len + 3*name_len) :: file_loc
        integer                              :: i, k, wext(3)

        file_loc = f_amr_restart_path(t_step)
        open (2, FILE=trim(file_loc), form='unformatted', STATUS='new')
        write (2) num_procs, amr_num_blocks, sys_size, amr_ref_ratio, amr_max_level, storage_size(0._stp)/8
        do k = 1, amr_num_blocks
            ! per-block header: region box, level (the reader needs it to rebuild multi-level geometry), extents; -1 extents
            ! mark a slot this rank does not own (the reader's data-presence flag)
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

    end subroutine s_amr_restart_write_serial

#ifdef MFC_MPI
    !> Parallel: one shared MPI-IO file.
    impure subroutine s_amr_restart_write_parallel(t_step)

        integer, intent(in) :: t_step
        character(LEN=path_len + 3*name_len) :: file_loc
        integer :: i, k
        integer :: ifile, ierr, cnt, idx, fi, fj, fk, bhdr(amr_restart_blk_hdr_ints), ibytes, sbytes
        integer, allocatable :: myown_all(:)
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: my_off, disp0, ddisp
        integer(kind=MPI_OFFSET_KIND), allocatable :: my_cnt_vec(:), my_off_vec(:), tot_cnt_vec(:)
        logical :: file_exist
        real(stp), allocatable :: buf(:)

        ibytes = storage_size(0)/8; sbytes = storage_size(0._stp)/8
        file_loc = f_amr_restart_path(t_step)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) then
            call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if
        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
        ! MPI-IO errors are silent unless checked
        if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart write: MPI_FILE_OPEN failed for ' // trim(file_loc))
        if (proc_rank == 0) call MPI_FILE_WRITE_AT(ifile, int(0, MPI_OFFSET_KIND), [-num_procs, amr_num_blocks, sys_size, &
            & amr_ref_ratio, amr_max_level, sbytes], amr_restart_ghdr_ints, MPI_INTEGER, status, ierr)
        disp0 = int(amr_restart_ghdr_ints*ibytes, MPI_OFFSET_KIND)  ! running byte offset past the global header
        ! one EXSCAN/ALLREDUCE over all blocks, not per block
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
        ! Per-block (owner + 1, m, n, p): every rank contributes zeros but the owner, so MAX recovers the record. The reader
        ! checks it against its own decomposition, catching an ownership pattern that would misalign the concatenated data.
        call MPI_ALLREDUCE(MPI_IN_PLACE, myown_all, amr_restart_blk_own_ints*amr_num_blocks, MPI_INTEGER, MPI_MAX, &
                           & MPI_COMM_WORLD, ierr)
        do k = 1, amr_num_blocks
            cnt = int(my_cnt_vec(k), kind(cnt))
            my_off = my_off_vec(k)
            if (proc_rank == 0) then
                ! block header (layout single-sourced in m_constants): region box (6) + level
                bhdr(1:3) = amr_region_lo_all(:,k); bhdr(4:6) = amr_region_hi_all(:,k)
                bhdr(amr_restart_blk_hdr_ints) = amr_block_level(k)
                call MPI_FILE_WRITE_AT(ifile, disp0, bhdr, amr_restart_blk_hdr_ints, MPI_INTEGER, status, ierr)
            end if
            if (proc_rank == 0) then
                call MPI_FILE_WRITE_AT(ifile, disp0 + int(amr_restart_blk_hdr_ints*ibytes, MPI_OFFSET_KIND), &
                                       & myown_all(amr_restart_blk_own_ints*(k - 1) + 1), amr_restart_blk_own_ints, MPI_INTEGER, &
                                       & status, ierr)
            end if
            ddisp = disp0 + int((amr_restart_blk_hdr_ints + amr_restart_blk_own_ints)*ibytes, MPI_OFFSET_KIND)
            allocate (buf(max(cnt, 1)))
            idx = 0
            ! non-owners have cnt == 0 and an unallocated q_cons, so skip packing; the collective write below still runs
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
        ! close is where buffered MPI-IO data flushes: a failure here truncates the file
        call MPI_FILE_CLOSE(ifile, ierr)
        if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart write: MPI_FILE_CLOSE failed; the file may be truncated')

    end subroutine s_amr_restart_write_parallel
#endif

    !> Write the fine-level restart file for save step t_step alongside the level-0 restart: the writing rank count, the active
    !! block count, and per block its box + level and each rank's fine conservative state. Serial: one unformatted file per rank in
    !! its level-0 step directory. Parallel: one shared MPI-IO file (header [-np, nboxes, sys_size, ratio, max level, word bytes],
    !! then per block a box+level header, a 4-int ownership record [owner+1, m, n, p], then the ranks' blocks concatenated in rank
    !! order). At the same rank count the ownership record must match on read; a different rank count repartitions.
    impure subroutine s_write_amr_restart(t_step)

        integer, intent(in) :: t_step
        integer             :: k

        if (.not. amr) return
        ! host consumer: fine state is device-current during stepping
        do k = 1, amr_num_blocks
            if (amr_owns_all(k)) then
                $:GPU_UPDATE(host='[amr_cons_st(:, :, :, :, amr_loc_of(k))]')
            end if
        end do

        if (.not. parallel_io) then
            call s_amr_restart_write_serial(t_step)
        else
#ifdef MFC_MPI
            call s_amr_restart_write_parallel(t_step)
#endif
        end if

    end subroutine s_write_amr_restart

    !> Restore the fine level from the AMR restart file at t_step_start (n_start under cfl_dt), if one exists: rebuild each saved
    !! block via s_set_amr_fine_geometry and read the fine state (exact stp round-trip). parallel_io repartitions across rank
    !! counts; serial needs the writing rank count. restored = false on a fresh start (the caller prolongs from coarse); a missing
    !! file aborts. Collective.
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
        if (have_glb == 0) call s_mpi_abort('amr restart: no AMR fine-level restart file ' // trim(file_loc) &
            & // ' (missing on at ' // 'least one rank); the fine level cannot be restored')

        if (parallel_io) then
            call s_amr_restart_read_parallel(file_loc)
        else
            call s_amr_restart_read_serial(file_loc)
        end if

        call s_amr_select_slot(1)
        ! the per-level fill waves iterate 2..amr_num_levels, which would skip every level>=2 fill until the first regrid
        amr_num_levels = max(1, maxval(amr_block_level(1:amr_num_blocks)))
        amr_seam_pairs_dirty = .true.  ! restored a new block set: the cached seam-pair list must be rebuilt
        amr_mesh_epoch = amr_mesh_epoch + 1
        call s_amr_check_seam_topology()  ! abort on seam topologies no halo reconciles (e.g. restart mode-switch)
        restored = .true.

    end subroutine s_read_amr_restart

    !> Install the restored block set: owners from the regions (or the recorded table), this run's owned slots (freeing stale init
    !! slots), then each block's geometry.
    impure subroutine s_amr_restart_install_blocks(own)

        integer, intent(in), optional :: own(:)  !< recorded owner per block, when it is valid for this run
        integer                       :: k

        if (present(own)) then
            call s_amr_adopt_block_owners(own)
        else
            call s_amr_assign_block_owners()
        end if
        call s_amr_reconcile_slots()
        do k = 1, amr_num_blocks
            amr_cur = k
            call s_set_amr_fine_geometry(amr_region_lo_all(:,k), amr_region_hi_all(:,k))
        end do
        call s_amr_reduce_xchg_flag()

    end subroutine s_amr_restart_install_blocks

    !> Serial (per-rank files) read: every block's region, level and, where this rank owned it at write (rm >= 0), its fine state.
    !! The presence flag is verified against this run's ownership once the owners are rebuilt.
    impure subroutine s_amr_restart_read_serial(file_loc)

        character(LEN=*), intent(in) :: file_loc
        integer                      :: i, k, ghdr(amr_restart_ghdr_ints), reg(6), lvl, rm, rn, rp, nblk
        logical, allocatable         :: had_data(:)

        open (2, FILE=trim(file_loc), form='unformatted', ACTION='read', STATUS='old')
        read (2) ghdr
        call s_amr_restart_check_header(ghdr, sys_size, 'amr restart', nblk, amr_ref_ratio, amr_max_level)
        call s_amr_restart_check_blocks(nblk)
        amr_num_blocks = nblk
        allocate (had_data(amr_num_blocks))
        do k = 1, amr_num_blocks
            read (2) reg, lvl, rm, rn, rp
            call s_amr_restart_check_record(k, reg, lvl)
            had_data(k) = rm >= 0
            if (.not. had_data(k)) cycle
            ! owner extents are region-derived per level; a disagreeing stored extent is corrupt
            if (rm /= (amr_ref_ratio**lvl)*(reg(4) - reg(1) + 1) - 1 .or. rn /= merge((amr_ref_ratio**lvl)*(reg(5) - reg(2) + 1) &
                & - 1, 0, n_glb > 0) .or. rp /= merge((amr_ref_ratio**lvl)*(reg(6) - reg(3) + 1) - 1, 0, p_glb > 0)) then
                call s_mpi_abort('amr restart: block fine extents disagree with the region (corrupt file)')
            end if
            call s_amr_alloc_slot(k)
            ! zero the whole host column first: the read fills only the interior but the push covers the padded column
            amr_cons_st(:,:,:,:,amr_loc_of(k)) = 0._stp
            do i = 1, sys_size
                read (2) amr_cons_st(0:rm,0:rn,0:rp,i, amr_loc_of(k))
            end do
            ! Push this block before the next s_amr_alloc_slot: a store growth pulls device->host first and would overwrite this
            ! block with the stale device copy. The store is device-authoritative at every alloc point.
            $:GPU_UPDATE(device='[amr_cons_st(:, :, :, :, amr_loc_of(k))]')
        end do
        close (2)
        call s_amr_restart_install_blocks()
        do k = 1, amr_num_blocks
            if (had_data(k) .neqv. amr_owns_all(k)) then
                call s_mpi_abort('amr restart decomposition mismatch: the file''s block ownership differs from this' &
                                 & // ' run''s (identical decomposition and rank count required)')
            end if
        end do
        deallocate (had_data)

    end subroutine s_amr_restart_read_serial

    !> Parallel (MPI-IO) read. The writer's rank count sets only the file layout: each block's data is one contiguous chunk, so any
    !! rank count can read it; owners are re-assigned and each new owner reads its blocks. At the writer's rank count the file's
    !! ownership record must reproduce this run's decomposition.
    impure subroutine s_amr_restart_read_parallel(file_loc)

        character(LEN=*), intent(in) :: file_loc

#ifdef MFC_MPI
        integer                             :: i, k, ghdr(amr_restart_ghdr_ints), nblk, ifile, ierr, cnt, idx, fi, fj, fk
        integer                             :: mown(amr_restart_blk_own_ints)
        type(t_amr_restart_catalog)         :: cat
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND)       :: fsz
        real(stp), allocatable              :: buf(:)

        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
        if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart read: MPI_FILE_OPEN failed for ' // trim(file_loc))
        ! a read past EOF returns short with an uninitialized tail, so compare the file size to the catalog's layout below
        call MPI_FILE_GET_SIZE(ifile, fsz, ierr)
        call MPI_FILE_READ_AT_ALL(ifile, int(0, MPI_OFFSET_KIND), ghdr, amr_restart_ghdr_ints, MPI_INTEGER, status, ierr)
        call s_amr_restart_check_header(ghdr, sys_size, 'amr restart', nblk, amr_ref_ratio, amr_max_level)
        call s_amr_restart_check_blocks(nblk)
        amr_num_blocks = nblk
        call s_amr_restart_read_catalog(ifile, amr_num_blocks, sys_size, cat)
        if (cat%end_disp /= fsz) call s_mpi_abort('amr restart read: file size does not match the expected layout ' &
            & // '(truncated or corrupt amr restart file)')
        do k = 1, amr_num_blocks
            call s_amr_restart_check_record(k, cat%reg(:,k), cat%lvl(k))
        end do
        ! at the writer's rank count, keep the recorded owners (every fine block must carry one)
        if (abs(ghdr(1)) == num_procs .and. all(cat%own(1,1:amr_num_blocks) > 0)) then
            call s_amr_restart_install_blocks(cat%own(1,1:amr_num_blocks) - 1)
        else
            call s_amr_restart_install_blocks()
        end if
        do k = 1, amr_num_blocks
            cnt = 0
            if (amr_owns_all(k)) then
                cnt = sys_size*(amr_slots(k)%m + 1)*(amr_slots(k)%n + 1)*(amr_slots(k)%p + 1)
                mown = [proc_rank + 1, amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p]
                if (abs(ghdr(1)) == num_procs .and. any(cat%own(:,k) /= mown)) then
                    call s_mpi_abort('amr restart: the per-block owner/extent record in the file does not match ' &
                                     & // 'this run''s decomposition; with the same rank count the ownership must match ' &
                                     & // 'the run that wrote the restart')
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
        ! push the restored state to the device here, not in s_read_amr_restart: the serial reader pushes per block (a store growth
        ! while allocating the next slot would otherwise lose the host data), and a late whole-set push would copy back columns that
        ! growth and compaction left undefined
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

    !> Corrupt/foreign-file guard for block k's record (box outside the domain, level outside 1..amr_max_level, extent over the
    !! block cap). Stores the accepted region and level, which the owner assignment and geometry build key off.
    impure subroutine s_amr_restart_check_record(k, reg, lvl)

        integer, intent(in) :: k, reg(6), lvl

        if (reg(1) < 0 .or. reg(4) > m_glb .or. reg(1) > reg(4) .or. (n_glb > 0 .and. (reg(2) < 0 .or. reg(5) > n_glb .or. reg(2) &
            & > reg(5))) .or. (p_glb > 0 .and. (reg(3) < 0 .or. reg(6) > p_glb .or. reg(3) > reg(6)))) then
            call s_mpi_abort('amr restart: corrupt block record (box outside the global domain)')
        end if
        if (lvl < 1 .or. lvl > amr_max_level) then
            call s_mpi_abort('amr restart: corrupt block record (block level outside 1..amr_max_level)')
        end if
        ! slot coordinate arrays are sized to max_f, which derives from the rank count unless amr_max_grid_size is pinned
        if (any(merge((amr_ref_ratio**lvl)*(reg(4:6) - reg(1:3) + 1) - 1 > max_f, .false., amr_dim))) then
            call s_mpi_abort('amr restart: a block exceeds this rank count''s block cap; restart on no more ranks ' &
                             & // 'than wrote the checkpoint, or set amr_max_grid_size in both runs')
        end if
        amr_region_lo_all(:,k) = reg(1:3); amr_region_hi_all(:,k) = reg(4:6)
        amr_block_level(k) = lvl

    end subroutine s_amr_restart_check_record

end module m_amr_restart
