!>
!!@file
!!@brief Contains module m_amr_restart_io

#:include 'macros.fpp'

!> @brief The AMR fine-level restart file layout shared by m_amr_restart (writer/reader) and post_process's overlay reader
!! (m_data_input): path, fine cell-boundary subdivision, global-header check and the parallel_io per-block catalog. Callers read the
!! data themselves.
!!
!! Layout (parallel_io): amr_restart_ghdr_ints-int header [-num_procs, nblk, nvar, amr_ref_ratio, amr_max_level, bytes per stored
!! word] (the rank count is positive in the serial files), then per block amr_restart_blk_hdr_ints ints [region lo(3), hi(3), level], amr_restart_blk_own_ints ints [owner + 1, m, n, p] and the owner's
!! data chunk nvar*(m+1)*(n+1)*(p+1) stp words (absent without an owner). The serial per-rank file
!! (p_all/p<rank>/<step>/amr_fine.dat) holds the same records as unformatted sequential records, data only for owned blocks.
module m_amr_restart_io

#ifdef MFC_MPI
    use mpi
#endif
    use m_precision_select
    use m_global_parameters_common, only: case_dir, parallel_io, mpiiofs
    use m_mpi_common, only: proc_rank, num_procs, s_mpi_abort
    use m_constants, only: amr_restart_blk_hdr_ints, amr_restart_blk_own_ints, path_len, name_len

    implicit none

    private
    public :: f_amr_restart_path, s_amr_subdivide_cb, s_amr_restart_check_header, amr_restart_ghdr_ints

    integer, parameter :: amr_restart_ghdr_ints = 6  !< global header length (ints)

#ifdef MFC_MPI
    public :: t_amr_restart_catalog, s_amr_restart_read_catalog

    !> Every block's header, ownership record and data-chunk byte offset, read collectively by all ranks.
    type t_amr_restart_catalog
        integer                                    :: nblk = 0
        integer, allocatable                       :: reg(:,:)      !< (6, nblk): region lo(3), hi(3) in global coarse cells
        integer, allocatable                       :: lvl(:)        !< (nblk): refinement level
        integer, allocatable                       :: own(:,:)      !< (amr_restart_blk_own_ints, nblk): owner + 1, m, n, p
        integer(kind=MPI_OFFSET_KIND), allocatable :: data_disp(:)  !< (nblk): byte offset of the block's data chunk
        integer(kind=MPI_OFFSET_KIND)              :: end_disp = 0  !< byte offset past the last chunk = the file's expected size
    end type t_amr_restart_catalog
#endif

contains

    !> The fine-level restart file for time step t_step: one file per rank without parallel_io, one shared file with it.
    function f_amr_restart_path(t_step) result(file_loc)

        integer, intent(in)                  :: t_step
        character(LEN=path_len + 3*name_len) :: file_loc

        if (.not. parallel_io) then
            write (file_loc, '(A,I0,A,I0,A)') trim(case_dir) // '/p_all/p', proc_rank, '/', t_step, '/amr_fine.dat'
        else
            write (file_loc, '(A,I0,A)') 'amr_', t_step, '.dat'
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
        end if

    end function f_amr_restart_path

    !> Fine cell boundaries fcb(-1:nfine) by rr-way subdivision of the coarse cells of pcb starting at coarse index lo; pcb_lb is
    !! lbound(pcb) in the caller (an assumed-shape dummy resets it to 1). Coarse right edges are copied exactly, so a subdivision of
    !! a subdivision reproduces the deeper level bit-for-bit.
    pure subroutine s_amr_subdivide_cb(pcb, pcb_lb, lo, nfine, rr, fcb)

        real(wp), intent(in)    :: pcb(:)
        integer, intent(in)     :: pcb_lb, lo, nfine, rr
        real(wp), intent(inout) :: fcb(-1:)
        integer                 :: fi, c, off, k
        real(wp)                :: xl, xr

        off = 1 - pcb_lb  ! pcb(j) = coarse_cb(j + pcb_lb - 1); coarse_cb(c) = pcb(c + off)
        fcb(-1) = pcb(lo - 1 + off)
        do fi = 0, nfine
            c = lo + fi/rr
            xl = pcb(c - 1 + off); xr = pcb(c + off)
            k = mod(fi, rr)
            if (k == rr - 1) then
                fcb(fi) = xr
            else
                fcb(fi) = (real(rr - 1 - k, wp)*xl + real(k + 1, wp)*xr)/real(rr, wp)
            end if
        end do

    end subroutine s_amr_subdivide_cb

    !> Validate the global header: the rank count (serial files need a match), nvar, the stored word size and, when given, the run's
    !! amr_ref_ratio (equal) and amr_max_level (at least the file's). Returns the block count.
    impure subroutine s_amr_restart_check_header(ghdr, nvar, who, nblk, ref_ratio, max_level)

        integer, intent(in)           :: ghdr(amr_restart_ghdr_ints), nvar
        character(len=*), intent(in)  :: who  !< message prefix ('amr restart' / 'amr post')
        integer, intent(out)          :: nblk
        integer, intent(in), optional :: ref_ratio, max_level
        character(LEN=300)            :: msg

        if (.not. parallel_io .and. ghdr(1) /= num_procs) then
            write (msg, '(A,I0,A,I0,A)') trim(who) // ': the AMR restart file was written with ', ghdr(1), &
                   & ' ranks but this run has ', num_procs, '; restart with the same rank count, or use parallel_io'
            call s_mpi_abort(trim(msg))
        end if
        if (ghdr(3) /= nvar) then
            write (msg, '(A,I0,A,I0,A)') trim(who) // ': the AMR restart file has ', ghdr(3), &
                   & ' conserved variables but this run has ', nvar, &
                   & '; the physics configuration ' // '(num_fluids/model_eqns/bubbles/chemistry) must match the run that wrote it'
            call s_mpi_abort(trim(msg))
        end if
        if (ghdr(6) /= storage_size(0._stp)/8) call s_mpi_abort(trim(who) // ': the AMR restart file was written at a different ' &
            & // 'storage precision (mixed_precision) than this run')
        if (present(ref_ratio)) then
            if (ghdr(4) /= ref_ratio) call s_mpi_abort(trim(who) // ': amr_ref_ratio differs from the run that wrote the file')
        end if
        if (present(max_level)) then
            if (ghdr(5) > max_level) call s_mpi_abort(trim(who) // ': the file holds blocks deeper than this run''s amr_max_level')
        end if
        nblk = ghdr(2)

    end subroutine s_amr_restart_check_header

#ifdef MFC_MPI
    !> Read every block's header and ownership record from an open parallel_io file (rank 0 reads, all ranks receive) and lay out
    !! the data offsets.
    impure subroutine s_amr_restart_read_catalog(ifile, nblk, nvar, cat)

        integer, intent(in)                      :: ifile, nblk, nvar
        type(t_amr_restart_catalog), intent(out) :: cat
        integer, parameter                       :: rec = amr_restart_blk_hdr_ints + amr_restart_blk_own_ints
        integer                                  :: k, ierr, ibytes, sbytes, brec(rec)
        integer, dimension(MPI_STATUS_SIZE)      :: status
        integer(kind=MPI_OFFSET_KIND)            :: disp, chunk

        ibytes = storage_size(0)/8; sbytes = storage_size(0._stp)/8
        cat%nblk = nblk
        allocate (cat%reg(6, nblk), cat%lvl(nblk), cat%own(amr_restart_blk_own_ints, nblk), cat%data_disp(nblk))
        if (proc_rank == 0) then
            disp = int(amr_restart_ghdr_ints*ibytes, MPI_OFFSET_KIND)
            do k = 1, nblk
                call MPI_FILE_READ_AT(ifile, disp, brec, rec, MPI_INTEGER, status, ierr)
                cat%reg(:,k) = brec(1:6); cat%lvl(k) = brec(amr_restart_blk_hdr_ints)
                cat%own(:,k) = brec(amr_restart_blk_hdr_ints + 1:rec)
                cat%data_disp(k) = disp + int(rec*ibytes, MPI_OFFSET_KIND)
                ! a block with no owner contributed no data at all
                chunk = int(0, MPI_OFFSET_KIND)
                if (cat%own(1, k) > 0) chunk = int(nvar, MPI_OFFSET_KIND)*int(cat%own(2, k) + 1, MPI_OFFSET_KIND)*int(cat%own(3, &
                    & k) + 1, MPI_OFFSET_KIND)*int(cat%own(4, k) + 1, MPI_OFFSET_KIND)
                disp = cat%data_disp(k) + chunk*int(sbytes, MPI_OFFSET_KIND)
            end do
            cat%end_disp = disp
        end if
        call MPI_BCAST(cat%reg, 6*nblk, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(cat%lvl, nblk, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(cat%own, amr_restart_blk_own_ints*nblk, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(cat%data_disp, nblk, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(cat%end_disp, 1, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)

    end subroutine s_amr_restart_read_catalog
#endif
end module m_amr_restart_io
