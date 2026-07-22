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
    use m_constants, only: amr_restart_blk_hdr_ints
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_min
    use m_amr, only: amr_slots, amr_seam_pairs_dirty, s_amr_alloc_slot, s_amr_reconcile_slots, s_amr_assign_block_owners, &
        & s_amr_gather_coarse_patch_pbmv, s_amr_prolong_pbmv, s_set_amr_fine_geometry
    use m_amr_regrid, only: s_amr_check_seam_topology

    implicit none

    private
    public :: s_write_amr_restart, s_read_amr_restart

contains

    !> Write the fine-level restart file for save step t_step alongside the level-0 restart (whose format stays untouched): the
    !! writing rank count, the active-block count, and for EACH block its box + each rank's intersection-local fine conservative
    !! state. Serial mode: one unformatted file per rank inside its level-0 step directory. Parallel mode: one shared MPI-IO file
    !! (3-int global header [np, nboxes, sys_size], then per block a 7-int box+level header, a 3*np-int per-rank fine-extents record
    !! [m,n,p per rank, 0s for non-owners; validated on read], then the ranks' fine blocks concatenated in rank order). Same rank
    !! count + decomposition required to restart (enforced by the extents record).
    impure subroutine s_write_amr_restart(t_step)

        integer, intent(in)                  :: t_step
        character(LEN=path_len + 3*name_len) :: file_loc
        integer                              :: i, k

#ifdef MFC_MPI
        integer :: ifile, ierr, cnt, idx, fi, fj, fk, reg(6), bhdr(amr_restart_blk_hdr_ints), ibytes, sbytes
        integer :: myext(3)
        integer, allocatable :: wext(:), myext_all(:), wext_all(:)
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: my_cnt, my_off, disp0, ddisp
        integer(kind=MPI_OFFSET_KIND), allocatable :: my_cnt_vec(:), my_off_vec(:), tot_cnt_vec(:)
        logical :: file_exist
        real(stp), allocatable :: buf(:)
#endif

        if (.not. amr) return
        ! host consumer: fine state is device-current during stepping (pull every owned slot)
        do k = 1, amr_num_blocks
            if (amr_owns_all(k)) then
                do i = 1, sys_size
                    $:GPU_UPDATE(host='[amr_slots(k)%q_cons(i)%sf]')
                end do
            end if
        end do

        if (.not. parallel_io) then
            ! per-rank file in the step directory just created by the level-0 serial write
            write (file_loc, '(A,I0,A,I0,A)') trim(case_dir) // '/p_all/p', proc_rank, '/', t_step, '/amr_fine.dat'
            open (2, FILE=trim(file_loc), form='unformatted', STATUS='new')
            write (2) num_procs, amr_num_blocks, sys_size
            do k = 1, amr_num_blocks
                ! per-block header: region box, refinement LEVEL (a level-l block's fine extent is amr_ref_ratio**l, not
                ! amr_ref_ratio, of the region - the reader needs the level to rebuild multi-level geometry), extents
                write (2) amr_slots(k)%region%lo, amr_slots(k)%region%hi, amr_block_level(k), amr_slots(k)%m, amr_slots(k)%n, &
                       & amr_slots(k)%p
                if (amr_owns_all(k)) then
                    do i = 1, sys_size
                        write (2) amr_slots(k)%q_cons(i)%sf(0:amr_slots(k)%m,0:amr_slots(k)%n,0:amr_slots(k)%p)
                    end do
                end if
            end do
            close (2)
        else
#ifdef MFC_MPI
            ibytes = storage_size(0)/8; sbytes = storage_size(0._stp)/8
            write (file_loc, '(A,I0,A)') 'amr_', t_step, '.dat'
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
            ! MPI-IO file handles default to MPI_ERRORS_RETURN: failures silent unless checked
            if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart write: MPI_FILE_OPEN failed for ' // trim(file_loc))
            if (proc_rank == 0) call MPI_FILE_WRITE_AT(ifile, int(0, MPI_OFFSET_KIND), [num_procs, amr_num_blocks, sys_size], 3, &
                & MPI_INTEGER, status, ierr)
            disp0 = int(3*ibytes, MPI_OFFSET_KIND)  ! running byte offset past the 3-int global header
            ! hoist per-block metadata collectives: one EXSCAN/ALLREDUCE/ALLGATHER over ALL blocks
            allocate (my_cnt_vec(amr_num_blocks), my_off_vec(amr_num_blocks), tot_cnt_vec(amr_num_blocks))
            allocate (myext_all(3*amr_num_blocks), wext_all(3*num_procs*amr_num_blocks))
            do k = 1, amr_num_blocks
                cnt = sys_size*(amr_slots(k)%m + 1)*(amr_slots(k)%n + 1)*(amr_slots(k)%p + 1)
                if (.not. amr_owns_all(k)) cnt = 0
                my_cnt_vec(k) = int(cnt, MPI_OFFSET_KIND)
                myext_all(3*(k - 1) + 1:3*(k - 1) + 3) = 0
                if (amr_owns_all(k)) myext_all(3*(k - 1) + 1:3*(k - 1) + 3) = [amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p]
            end do
            my_off_vec = int(0, MPI_OFFSET_KIND)
            call MPI_EXSCAN(my_cnt_vec, my_off_vec, amr_num_blocks, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) my_off_vec = int(0, MPI_OFFSET_KIND)
            call MPI_ALLREDUCE(my_cnt_vec, tot_cnt_vec, amr_num_blocks, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
            ! per-rank fine extents (0s for non-owning ranks): readers rebuild this vector from their own decomposition and abort on
            ! mismatch - a different rank count, ownership pattern, or load_balance split would otherwise silently misalign the
            ! concatenated per-rank data slices below
            call MPI_ALLGATHER(myext_all, 3*amr_num_blocks, MPI_INTEGER, wext_all, 3*amr_num_blocks, MPI_INTEGER, MPI_COMM_WORLD, &
                               & ierr)
            if (.not. allocated(wext)) allocate (wext(3*num_procs))
            do k = 1, amr_num_blocks
                cnt = int(my_cnt_vec(k), kind(cnt))
                my_off = my_off_vec(k)
                if (proc_rank == 0) then
                    ! amr_restart_blk_hdr_ints-int per-block header: region box (6) + refinement LEVEL (a level-l block's fine
                    ! extent is amr_ref_ratio**l, not amr_ref_ratio, of the region - the reader needs the level to rebuild
                    ! multi-level geometry). Header layout is single-sourced in m_constants so both readers stay in lockstep.
                    bhdr(1:3) = amr_slots(k)%region%lo; bhdr(4:6) = amr_slots(k)%region%hi
                    bhdr(amr_restart_blk_hdr_ints) = amr_block_level(k)
                    call MPI_FILE_WRITE_AT(ifile, disp0, bhdr, amr_restart_blk_hdr_ints, MPI_INTEGER, status, ierr)
                end if
                ! wext_all layout: rank r's extents for block k at wext_all(3*amr_num_blocks*r + 3*(k-1) + 1 : +3)
                do i = 0, num_procs - 1
                    wext(3*i + 1:3*i + 3) = wext_all(3*amr_num_blocks*i + 3*(k - 1) + 1:3*amr_num_blocks*i + 3*(k - 1) + 3)
                end do
                if (proc_rank == 0) then
                    call MPI_FILE_WRITE_AT(ifile, disp0 + int(amr_restart_blk_hdr_ints*ibytes, MPI_OFFSET_KIND), wext, &
                                           & 3*num_procs, MPI_INTEGER, status, ierr)
                end if
                ddisp = disp0 + int((amr_restart_blk_hdr_ints + 3*num_procs)*ibytes, MPI_OFFSET_KIND)
                allocate (buf(max(cnt, 1)))
                idx = 0
                do i = 1, sys_size
                    do fk = 0, amr_slots(k)%p
                        do fj = 0, amr_slots(k)%n
                            do fi = 0, amr_slots(k)%m
                                idx = idx + 1
                                buf(idx) = amr_slots(k)%q_cons(i)%sf(fi, fj, fk)
                            end do
                        end do
                    end do
                end do
                call MPI_FILE_WRITE_AT_ALL(ifile, ddisp + my_off*int(sbytes, MPI_OFFSET_KIND), buf, cnt*mpi_io_type, mpi_io_p, &
                                           & status, ierr)
                if (ierr /= MPI_SUCCESS) &
                    & call s_mpi_abort('amr restart write: data write failed (disk full/quota?); the file is unusable')
                deallocate (buf)
                disp0 = ddisp + tot_cnt_vec(k)*int(sbytes, MPI_OFFSET_KIND)
            end do
            deallocate (my_cnt_vec, my_off_vec, tot_cnt_vec, myext_all, wext_all)
            ! close is where buffered MPI-IO data flushes on many stacks - a failure here truncates the file
            call MPI_FILE_CLOSE(ifile, ierr)
            if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart write: MPI_FILE_CLOSE failed; the file may be truncated')
#endif
        end if

    end subroutine s_write_amr_restart

    !> Restore the fine level from the AMR restart file at t_step_start (n_start under cfl_dt), if one exists: for each saved block
    !! rebuild the box via s_set_amr_fine_geometry, then read each rank's intersection-local fine state (exact stp round-trip).
    !! parallel_io REPARTITIONS across rank counts (each block is one contiguous region-sized chunk under whole-block ownership,
    !! re-assigned to this run's owners); serial (per-rank files) still needs the writing rank count. restored = false on a fresh
    !! start, or - with a one-line warning - on a legacy restart without the file; the caller then re-prolongs from coarse.
    !! Collective: ALL ranks call it together.
    impure subroutine s_read_amr_restart(restored)

        logical, intent(out)                 :: restored
        character(LEN=path_len + 3*name_len) :: file_loc
        character(LEN=300)                   :: msg
        logical                              :: file_exist
        integer                              :: i, k, ts, have_loc, have_glb, ghdr(3), reg(6), lvl, rm, rn, rp
        logical, allocatable                 :: had_data(:)

#ifdef MFC_MPI
        integer :: ifile, ierr, cnt, idx, fi, fj, fk, ibytes, sbytes, np_old, bhdr(amr_restart_blk_hdr_ints)
        integer :: myext(3)
        integer, allocatable :: wext(:), rext(:), myext_all(:), wext_all(:)
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: my_cnt, my_off, disp0, ddisp, fsz
        integer(kind=MPI_OFFSET_KIND), allocatable :: blk_base(:), my_cnt_vec(:), my_off_vec(:)
        real(stp), allocatable :: buf(:)
#endif

        restored = .false.
        if (.not. amr) return
        if (cfl_dt) then
            ts = n_start
        else
            ts = t_step_start
        end if
        if (ts == 0) return  ! fresh start: the fine level is prolonged from the pre_process ICs

        if (.not. parallel_io) then
            write (file_loc, '(A,I0,A,I0,A)') trim(case_dir) // '/p_all/p', proc_rank, '/', ts, '/amr_fine.dat'
        else
            write (file_loc, '(A,I0,A)') 'amr_', ts, '.dat'
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
        end if
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

        if (.not. parallel_io) then
            open (2, FILE=trim(file_loc), form='unformatted', ACTION='read', STATUS='old')
            read (2) ghdr
            if (ghdr(1) /= num_procs) then
                write (msg, &
                       & '(A,I0,A,I0,A)') 'amr restart rank-count mismatch: the serial (non-parallel_io) AMR restart ' &
                       & // 'file was written with ', ghdr(1), ' ranks but this run has ', num_procs, &
                       & '; restart with the same rank count, or use parallel_io (which repartitions across rank counts)'
                call s_mpi_abort(trim(msg))
            end if
            if (ghdr(3) /= sys_size) then
                write (msg, '(A,I0,A,I0,A)') 'amr restart sys_size mismatch: the AMR restart file has ', ghdr(3), &
                       & ' conserved variables but this run has ', sys_size, &
                       & '; the physics configuration ' &
                       & // '(num_fluids/model_eqns/bubbles/chemistry) must match the run that wrote the restart'
                call s_mpi_abort(trim(msg))
            end if
            if (ghdr(2) < 1 .or. ghdr(2) > amr_max_blocks) then
                call s_mpi_abort('amr restart: the file holds more fine blocks than amr_max_blocks ' &
                                 & // 'in this run; restart with amr_max_blocks at least the written block count')
            end if
            amr_num_blocks = ghdr(2)
            allocate (had_data(amr_num_blocks))
            ! PASS 1: read every block's region + (present iff rm>=0, i.e. this rank owned it at write) the owner's fine state.
            ! Whole-block ownership is decomposition-deterministic, so the file's data-presence flag drives the read here; the owner
            ! map is rebuilt from the regions in pass 2.
            do k = 1, amr_num_blocks
                read (2) reg, lvl, rm, rn, rp
                ! corrupt/foreign-file guard: a box outside the global domain would drive the geometry build and coordinate reads
                ! out of bounds silently in release builds
                if (reg(1) < 0 .or. reg(4) > m_glb .or. reg(1) > reg(4) .or. (n_glb > 0 .and. (reg(2) < 0 .or. reg(5) > n_glb &
                    & .or. reg(2) > reg(5))) .or. (p_glb > 0 .and. (reg(3) < 0 .or. reg(6) > p_glb .or. reg(3) > reg(6)))) then
                    call s_mpi_abort('amr restart: corrupt block record (box outside the global domain)')
                end if
                if (lvl < 1 .or. lvl > amr_max_level) then
                    call s_mpi_abort('amr restart: corrupt block record (block level outside 1..amr_max_level)')
                end if
                amr_region_lo_all(:,k) = reg(1:3); amr_region_hi_all(:,k) = reg(4:6)
                ! set the level BEFORE the owner/geometry rebuild below: s_amr_assign_block_owners and s_set_amr_fine_geometry key
                ! off amr_block_level to place L>=2 blocks under their parent
                amr_block_level(k) = lvl
                had_data(k) = rm >= 0
                if (had_data(k)) then
                    ! whole-block owner extents are region-derived per level (a level-l block covers amr_ref_ratio**l fine cells per
                    ! L0 cell of its region, not amr_ref_ratio); a stored extent that disagrees is corrupt
                    if (rm /= (amr_ref_ratio**lvl)*(reg(4) - reg(1) + 1) - 1 .or. rn /= merge((amr_ref_ratio**lvl)*(reg(5) &
                        & - reg(2) + 1) - 1, 0, n_glb > 0) .or. rp /= merge((amr_ref_ratio**lvl)*(reg(6) - reg(3) + 1) - 1, 0, &
                        & p_glb > 0)) then
                        call s_mpi_abort('amr restart: block fine extents disagree with the region (corrupt file)')
                    end if
                    ! serial (same rank count): had_data == this run's ownership, so this is the owned slot
                    call s_amr_alloc_slot(k)
                    do i = 1, sys_size
                        read (2) amr_slots(k)%q_cons(i)%sf(0:rm,0:rn,0:rp)
                    end do
                end if
            end do
            close (2)
            ! PASS 2: rebuild whole-block owners from the regions, then each block's geometry under the correct owner; verify the
            ! data read (write-owner) matches who owns the block in this run
            call s_amr_assign_block_owners()
            ! free any init slots not in the restart set (had_data slots stay: they are owned)
            call s_amr_reconcile_slots()
            do k = 1, amr_num_blocks
                amr_cur = k
                call s_set_amr_fine_geometry(amr_region_lo_all(:,k), amr_region_hi_all(:,k))
                if (had_data(k) .neqv. amr_owns_all(k)) then
                    call s_mpi_abort('amr restart decomposition mismatch: the file''s block ownership differs from this' &
                                     & // ' run''s (identical decomposition - rank count and load_balance settings - required)')
                end if
            end do
            deallocate (had_data)
        else
#ifdef MFC_MPI
            ibytes = storage_size(0)/8; sbytes = storage_size(0._stp)/8
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            ! MPI-IO errors are silent by default (MPI_ERRORS_RETURN on file handles) and a read past EOF is not even an error - it
            ! returns short with an uninitialized tail. Grab the size up front; the exact expected byte count is compared after the
            ! layout records are consumed below.
            if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart read: MPI_FILE_OPEN failed for ' // trim(file_loc))
            call MPI_FILE_GET_SIZE(ifile, fsz, ierr)
            call MPI_FILE_READ_AT_ALL(ifile, int(0, MPI_OFFSET_KIND), ghdr, 3, MPI_INTEGER, status, ierr)
            ! Repartition-on-restart: the writer's rank count sets only the file layout (the 3*np_old per-block extents record).
            ! Whole-block ownership makes each block's fine data one contiguous region-sized chunk, so ANY new rank count can read
            ! it
            ! - pass 2 re-assigns owners for THIS run and each new owner reads its whole blocks. np_old == num_procs is
            ! byte-identical to the same-rank path (and keeps the layout check).
            np_old = ghdr(1)
            if (ghdr(3) /= sys_size) then
                write (msg, '(A,I0,A,I0,A)') 'amr restart sys_size mismatch: the AMR restart file has ', ghdr(3), &
                       & ' conserved variables but this run has ', sys_size, &
                       & '; the physics configuration ' &
                       & // '(num_fluids/model_eqns/bubbles/chemistry) must match the run that wrote the restart'
                call s_mpi_abort(trim(msg))
            end if
            if (ghdr(2) < 1 .or. ghdr(2) > amr_max_blocks) then
                call s_mpi_abort('amr restart: the file holds more fine blocks than amr_max_blocks ' &
                                 & // 'in this run; restart with amr_max_blocks at least the written block count')
            end if
            amr_num_blocks = ghdr(2)
            allocate (wext(3*np_old), rext(3*num_procs), blk_base(amr_num_blocks))
            ! PASS 1: read every block's region (collective) and lay out the file offsets. Under whole-block ownership the per-block
            ! data size is fixed by the region (one owner holds all sys_size*cells), so all offsets are known before the owner map
            ! is rebuilt in pass 2.
            disp0 = int(3*ibytes, MPI_OFFSET_KIND)
            do k = 1, amr_num_blocks
                call MPI_FILE_READ_AT_ALL(ifile, disp0, bhdr, amr_restart_blk_hdr_ints, MPI_INTEGER, status, ierr)
                reg = bhdr(1:6); lvl = bhdr(amr_restart_blk_hdr_ints)
                ! corrupt/foreign-file guard: a box outside the global domain would drive the geometry build and coordinate reads
                ! out of bounds silently in release builds
                if (reg(1) < 0 .or. reg(4) > m_glb .or. reg(1) > reg(4) .or. (n_glb > 0 .and. (reg(2) < 0 .or. reg(5) > n_glb &
                    & .or. reg(2) > reg(5))) .or. (p_glb > 0 .and. (reg(3) < 0 .or. reg(6) > p_glb .or. reg(3) > reg(6)))) then
                    call s_mpi_abort('amr restart: corrupt block record (box outside the global domain)')
                end if
                if (lvl < 1 .or. lvl > amr_max_level) then
                    call s_mpi_abort('amr restart: corrupt block record (block level outside 1..amr_max_level)')
                end if
                amr_region_lo_all(:,k) = reg(1:3); amr_region_hi_all(:,k) = reg(4:6)
                ! set the level before the owner/geometry rebuild: s_amr_assign_block_owners and s_set_amr_fine_geometry key off
                ! amr_block_level to place L>=2 blocks under their parent
                amr_block_level(k) = lvl
                blk_base(k) = disp0
                ! data size is region-derived per level: a level-l block covers amr_ref_ratio**l fine cells per L0 cell
                cnt = sys_size*((amr_ref_ratio**lvl)*(reg(4) - reg(1) + 1))*merge((amr_ref_ratio**lvl)*(reg(5) - reg(2) + 1), 1, &
                                & n_glb > 0)*merge((amr_ref_ratio**lvl)*(reg(6) - reg(3) + 1), 1, p_glb > 0)
                disp0 = disp0 + int((amr_restart_blk_hdr_ints + 3*np_old)*ibytes, MPI_OFFSET_KIND) + int(cnt, &
                                    & MPI_OFFSET_KIND)*int(sbytes, MPI_OFFSET_KIND)
            end do
            ! PASS 2: rebuild whole-block owners from the regions, then per block build geometry under the correct owner, validate
            ! the writer's layout, and read this rank's owned slice at its offset.
            call s_amr_assign_block_owners()
            ! allocate this run's owned blocks (frees any stale init slots) before the read below
            call s_amr_reconcile_slots()
            do k = 1, amr_num_blocks
                amr_cur = k
                call s_set_amr_fine_geometry(amr_region_lo_all(:,k), amr_region_hi_all(:,k))
            end do
            ! hoist per-block metadata collectives: one ALLGATHER/EXSCAN over ALL blocks
            allocate (my_cnt_vec(amr_num_blocks), my_off_vec(amr_num_blocks))
            allocate (myext_all(3*amr_num_blocks), wext_all(3*num_procs*amr_num_blocks))
            do k = 1, amr_num_blocks
                cnt = sys_size*(amr_slots(k)%m + 1)*(amr_slots(k)%n + 1)*(amr_slots(k)%p + 1)
                if (.not. amr_owns_all(k)) cnt = 0
                my_cnt_vec(k) = int(cnt, MPI_OFFSET_KIND)
                myext_all(3*(k - 1) + 1:3*(k - 1) + 3) = 0
                if (amr_owns_all(k)) myext_all(3*(k - 1) + 1:3*(k - 1) + 3) = [amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p]
            end do
            my_off_vec = int(0, MPI_OFFSET_KIND)
            call MPI_EXSCAN(my_cnt_vec, my_off_vec, amr_num_blocks, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) my_off_vec = int(0, MPI_OFFSET_KIND)
            ! same rank count: validate the writer's per-rank layout against this run's decomposition (a re-derived load_balance
            ! split would silently misalign every rank's slice). Repartitioning (np_old /= num_procs) intentionally uses a DIFFERENT
            ! decomposition, so the layout cannot match - skip the check; whole-block ownership makes each block one contiguous
            ! chunk
            ! the new owner reads wholly, and the file-size check below still fails closed on a truncated/corrupt file.
            if (np_old == num_procs) then
                call MPI_ALLGATHER(myext_all, 3*amr_num_blocks, MPI_INTEGER, wext_all, 3*amr_num_blocks, MPI_INTEGER, &
                                   & MPI_COMM_WORLD, ierr)
            end if
            if (.not. allocated(wext)) allocate (wext(3*np_old))
            if (.not. allocated(rext)) allocate (rext(3*num_procs))
            do k = 1, amr_num_blocks
                cnt = int(my_cnt_vec(k), kind(cnt))
                my_off = my_off_vec(k)
                if (np_old == num_procs) then
                    call MPI_FILE_READ_AT_ALL(ifile, blk_base(k) + int(amr_restart_blk_hdr_ints*ibytes, MPI_OFFSET_KIND), wext, &
                                              & 3*np_old, MPI_INTEGER, status, ierr)
                    do i = 0, num_procs - 1
                        rext(3*i + 1:3*i + 3) = wext_all(3*amr_num_blocks*i + 3*(k - 1) + 1:3*amr_num_blocks*i + 3*(k - 1) + 3)
                    end do
                    if (any(rext /= wext)) then
                        call s_mpi_abort('amr restart: the per-rank fine-block layout in the file does not match ' &
                                         & // 'this run''s decomposition; with the same rank count the ownership and ' &
                                         & // '(with load_balance) the weighted splits must match the run that wrote the restart')
                    end if
                end if
                ddisp = blk_base(k) + int((amr_restart_blk_hdr_ints + 3*np_old)*ibytes, MPI_OFFSET_KIND)
                allocate (buf(max(cnt, 1)))
                call MPI_FILE_READ_AT_ALL(ifile, ddisp + my_off*int(sbytes, MPI_OFFSET_KIND), buf, cnt*mpi_io_type, mpi_io_p, &
                                          & status, ierr)
                idx = 0
                do i = 1, sys_size
                    do fk = 0, amr_slots(k)%p
                        do fj = 0, amr_slots(k)%n
                            do fi = 0, amr_slots(k)%m
                                idx = idx + 1
                                amr_slots(k)%q_cons(i)%sf(fi, fj, fk) = buf(idx)
                            end do
                        end do
                    end do
                end do
                deallocate (buf)
            end do
            deallocate (blk_base, my_cnt_vec, my_off_vec, myext_all, wext_all)
            ! disp0 now equals the exact byte count a complete file must have: a truncated file (crashed writer, filesystem hiccup)
            ! passes every layout check above but returns short reads with garbage tails - fail closed instead of restoring
            ! uninitialized data as the fine level
            if (disp0 /= fsz) then
                call s_mpi_abort('amr restart read: file size does not match the expected layout ' &
                                 & // '(truncated or corrupt amr restart file)')
            end if
            call MPI_FILE_CLOSE(ifile, ierr)
#endif
        end if

        ! push restored fine state to the device (mirrors s_populate_amr_fine's push; host reads above)
        do k = 1, amr_num_blocks
            if (amr_owns_all(k)) then
                do i = 1, sys_size
                    $:GPU_UPDATE(device='[amr_slots(k)%q_cons(i)%sf]')
                end do
            end if
        end do
        ! non-polytropic QBMM: the restart file carries q_cons only; re-prolong each block's side-state from the restored coarse
        ! pb/mv (one-time piecewise-constant smoothing)
        if (qbmm .and. .not. polytropic) then
            do k = 1, amr_num_blocks
                call s_amr_select_slot(k)
                ! gather coarse pb/mv patch on ALL ranks (P2P), then owners re-prolong from it
                call s_amr_gather_coarse_patch_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, .false.)
                if (amr_owns_all(k)) call s_amr_prolong_pbmv()
            end do
        end if
        call s_amr_select_slot(1)
        amr_seam_pairs_dirty = .true.  ! restored a new block set: the cached seam-pair list must be rebuilt
        call s_amr_check_seam_topology()  ! abort on seam topologies no halo reconciles (e.g. restart mode-switch)
        restored = .true.

    end subroutine s_read_amr_restart

end module m_amr_restart
