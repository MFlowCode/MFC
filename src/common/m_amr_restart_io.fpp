!>
!!@file
!!@brief Contains module m_amr_restart_io

#:include 'macros.fpp'

!> @brief The AMR fine-level restart file's layout, shared by the simulation's writer/reader (m_amr_restart) and post_process's
!! overlay reader (m_data_input): the file path, the fine cell-boundary subdivision, the global-header check and, under parallel_io,
!! the per-block catalog (header + ownership record + data offset of every block). The data itself is read by the callers, which
!! differ in what they keep (the simulation restores whole blocks into its store; post_process keeps each rank's intersection
!! slice).
!!
!! Layout (parallel_io): 3-int global header [-num_procs, nblk, nvar] (the NEGATIVE rank count marks the current format), then per
!! block: amr_restart_blk_hdr_ints ints [region lo(3), hi(3), level], amr_restart_blk_own_ints ints [owner + 1, m, n, p], and the
!! owner's contiguous data chunk nvar*(m+1)*(n+1)*(p+1) stp words (absent when the block has no owner). The serial per-rank file
!! (p_all/p<rank>/<step>/amr_fine.dat) holds the same records as Fortran unformatted sequential records, data only for owned blocks.
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
    public :: f_amr_restart_path, s_amr_subdivide_cb, s_amr_restart_check_header

#ifdef MFC_MPI
    public :: t_amr_restart_catalog, s_amr_restart_read_catalog

    !> Every block's header and ownership record plus the byte offset of its data chunk, read collectively by all ranks.
    type t_amr_restart_catalog
        integer                                    :: nblk = 0
        integer, allocatable                       :: reg(:,:)      !< (6, nblk): region lo(3), hi(3) in global coarse cells
        integer, allocatable                       :: lvl(:)        !< (nblk): refinement level (0 = an L0 tile, no fine data)
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
    !! lbound(pcb) in the caller (an assumed-shape dummy resets it to 1). Every coarse cell's right edge is copied exactly, so the
    !! subdivision of a subdivision reproduces the deeper level bit-for-bit, and at rr = 2 it is a plain bisection.
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

    !> Validate the 3-int global header against this run: the rank count (serial files are per rank, so it must match; a parallel_io
    !! file may be read at any rank count, but only in the current format) and the conserved-variable count nvar the caller expects.
    !! Returns the block count.
    impure subroutine s_amr_restart_check_header(ghdr, nvar, who, nblk)

        integer, intent(in)          :: ghdr(3), nvar
        character(len=*), intent(in) :: who  !< message prefix ('amr restart' / 'amr post')
        integer, intent(out)         :: nblk
        character(LEN=300)           :: msg

        if (parallel_io) then
            if (ghdr(1) > 0) call s_mpi_abort(trim(who) // ': the AMR restart file is in the retired per-rank-extent ' &
                & // 'layout; regenerate it with the current code')
        else if (ghdr(1) /= num_procs) then
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
        nblk = ghdr(2)

    end subroutine s_amr_restart_check_header

#ifdef MFC_MPI
    !> Read every block's header and ownership record from an open parallel_io file (collective: all ranks call, all get the same
    !! catalog) and lay out the data offsets. nvar is the conserved-variable count the writer used.
    impure subroutine s_amr_restart_read_catalog(ifile, nblk, nvar, cat)

        integer, intent(in)                      :: ifile, nblk, nvar
        type(t_amr_restart_catalog), intent(out) :: cat
        integer                                  :: k, ierr, ibytes, sbytes, bhdr(amr_restart_blk_hdr_ints)
        integer, dimension(MPI_STATUS_SIZE)      :: status
        integer(kind=MPI_OFFSET_KIND)            :: disp, chunk

        ibytes = storage_size(0)/8; sbytes = storage_size(0._stp)/8
        cat%nblk = nblk
        allocate (cat%reg(6, nblk), cat%lvl(nblk), cat%own(amr_restart_blk_own_ints, nblk), cat%data_disp(nblk))
        disp = int(3*ibytes, MPI_OFFSET_KIND)
        do k = 1, nblk
            call MPI_FILE_READ_AT_ALL(ifile, disp, bhdr, amr_restart_blk_hdr_ints, MPI_INTEGER, status, ierr)
            cat%reg(:,k) = bhdr(1:6); cat%lvl(k) = bhdr(amr_restart_blk_hdr_ints)
            call MPI_FILE_READ_AT_ALL(ifile, disp + int(amr_restart_blk_hdr_ints*ibytes, MPI_OFFSET_KIND), cat%own(:,k), &
                                      & amr_restart_blk_own_ints, MPI_INTEGER, status, ierr)
            cat%data_disp(k) = disp + int((amr_restart_blk_hdr_ints + amr_restart_blk_own_ints)*ibytes, MPI_OFFSET_KIND)
            ! a block with no owner contributed no data at all
            chunk = int(0, MPI_OFFSET_KIND)
            if (cat%own(1, k) > 0) chunk = int(nvar, MPI_OFFSET_KIND)*int(cat%own(2, k) + 1, MPI_OFFSET_KIND)*int(cat%own(3, &
                & k) + 1, MPI_OFFSET_KIND)*int(cat%own(4, k) + 1, MPI_OFFSET_KIND)
            disp = cat%data_disp(k) + chunk*int(sbytes, MPI_OFFSET_KIND)
        end do
        cat%end_disp = disp

    end subroutine s_amr_restart_read_catalog
#endif
end module m_amr_restart_io
