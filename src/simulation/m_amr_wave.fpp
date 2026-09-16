!>
!!@file
!!@brief Contains module m_amr_wave

#:include 'macros.fpp'

!> @brief The one place the AMR exchange waves talk to MPI. A wave is a set of point-to-point transfers posted together: all
!! receives, then the packs, then all sends, then one wait. Two layers:
!!
!! - the request layer (s_amr_wave_open / irecv / isend / wait): keyed tags, the exchange audit record per message, the request
!!   list and the debug check that every receive arrived with exactly the planned length. The zero-copy waves (reflux faces, freg)
!!   use it directly on their register arrays.
!! - the plan layer (t_amr_wave_side): the pooled waves (fills, seams, restricts) list their transfers as (peer, block, box, words)
!!   on a send side and a receive side; s_amr_wave_close lays every peer's transfers out contiguously in a pool with a debug
!!   header slot in front of each, and s_amr_wave_post / s_amr_wave_send issue one message per peer. Both ranks of a pair enumerate
!!   the same transfers in the same order from replicated metadata, so the layout agrees with no metadata exchange.
!!
!! Pools are module arrays owned by the caller (device-resident under rdma_mpi: MPI reads and writes them by device address).
module m_amr_wave

#ifdef MFC_MPI
    use mpi
#endif
    use m_precision_select
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort
    use m_amr_xchg_audit
    use m_amr_state, only: amr_fw_dev
    use m_amr_distribution, only: s_amr_m1_wave_open, f_amr_m1_seq, f_amr_m1_tag

    implicit none

    private
    public :: t_amr_wave_side, t_amr_wave, s_amr_wave_open, s_amr_wave_irecv, s_amr_wave_isend, s_amr_wave_irecv_raw, &
        & s_amr_wave_isend_raw, s_amr_wave_wait, f_amr_wave_nreq, s_amr_wave_reset, s_amr_wave_add, s_amr_wave_close, &
        & s_amr_wave_post, s_amr_wave_send, s_amr_wave_slice, s_amr_wave_hdr_pack, s_amr_wave_hdr_check, s_amr_wave_size_int, &
        & s_amr_wave_size_int3, s_amr_wave_size_real

    !> One side (send or receive) of a pooled wave. Transfers are appended in enumeration order; s_amr_wave_close groups them by
    !! peer into one contiguous pool run per peer, each transfer preceded by XA_NH header words.
    type t_amr_wave_side
        integer :: nx = 0     !< transfers
        integer :: np = 0     !< peers
        integer :: words = 0  !< pool words in use after close
        !> per transfer (off: pool offset of its header)
        integer, allocatable :: blk(:), bl(:,:), bh(:,:), cnt(:), peer(:), pi(:), off(:)
        integer, allocatable :: prank(:), pwords(:), pbase(:), pnx(:)  !< per peer: rank, message words, pool base, transfers
    end type t_amr_wave_side

    !> the pooled waves' send/receive sides. The fill and restrict waves run one at a time and share amr_wsend/amr_wrecv (and the
    !! amr_fw_sq/rq pools); the seam wave is posted before the fills and drained after them, so it keeps its own pair.
    type(t_amr_wave_side), public :: amr_wsend, amr_wrecv, amr_wseam_snd, amr_wseam_rcv

    !> An open wave: its tag band and request list; sends carry -1 in reqw, receives their planned word count (checked after the
    !! wait under MFC_DEBUG). The seam wave stays open across the fill waves, so it is a second instance.
    type t_amr_wave
        integer              :: band = 0, nreq = 0
        integer, allocatable :: req(:), reqw(:)
    end type t_amr_wave
    type(t_amr_wave), public, target :: amr_wave, amr_wave_seam
    !> rank-indexed build scratch for s_amr_wave_close (touched entries re-zeroed, so it stays all-zero between builds)
    integer, allocatable :: bmap(:), bnx(:), bwords(:)

contains

    !> Open a wave on tag band b: bumps the band generation, clears the keyed-tag seq counters and the request list.
    impure subroutine s_amr_wave_open(w, b)

        type(t_amr_wave), intent(inout) :: w
        integer, intent(in)             :: b

        w%band = b; w%nreq = 0
        call s_amr_m1_wave_open(b)

    end subroutine s_amr_wave_open

    pure integer function f_amr_wave_nreq(w) result(n)

        type(t_amr_wave), intent(in) :: w

        n = w%nreq

    end function f_amr_wave_nreq

    !> Post one receive of n words from peer into buf (device-resident when dev), recorded at audit site with key.
    impure subroutine s_amr_wave_irecv(w, buf, n, peer, site, key, dev, rec, nrec)

        type(t_amr_wave), intent(inout)     :: w
        real(wp), intent(inout), contiguous :: buf(:)
        integer, intent(in)                 :: n, peer, site, key
        logical, intent(in)                 :: dev
        logical, intent(in), optional       :: rec   !< record in the exchange audit (default; false for header-only messages)
        integer, intent(in), optional       :: nrec  !< payload words to record (default n; a pooled message also carries headers)
        integer                             :: ierr, sq, tq, nr
        logical                             :: do_rec

        do_rec = .true.; if (present(rec)) do_rec = rec
        nr = n; if (present(nrec)) nr = nrec
        sq = f_amr_m1_seq(peer, 2); tq = f_amr_m1_tag(w%band, sq)
        if (do_rec) call s_xa_rec(site, 2, nr, tq, peer=peer, key=key, seq=sq)
        w%nreq = w%nreq + 1
        call s_amr_wave_size_int(w%req, w%nreq); call s_amr_wave_size_int(w%reqw, w%nreq)
        w%reqw(w%nreq) = n
#ifdef MFC_MPI
        if (dev) then
            #:call GPU_HOST_DATA(use_device_addr='[buf]')
                call MPI_IRECV(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
            #:endcall GPU_HOST_DATA
        else
            call MPI_IRECV(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
        end if
#endif

    end subroutine s_amr_wave_irecv

    !> Post one send of n words from buf to peer (see s_amr_wave_irecv).
    impure subroutine s_amr_wave_isend(w, buf, n, peer, site, key, dev, rec, nrec)

        type(t_amr_wave), intent(inout)  :: w
        real(wp), intent(in), contiguous :: buf(:)
        integer, intent(in)              :: n, peer, site, key
        logical, intent(in)              :: dev
        logical, intent(in), optional    :: rec
        integer, intent(in), optional    :: nrec
        integer                          :: ierr, sq, tq, nr
        logical                          :: do_rec

        do_rec = .true.; if (present(rec)) do_rec = rec
        nr = n; if (present(nrec)) nr = nrec
        sq = f_amr_m1_seq(peer, 1); tq = f_amr_m1_tag(w%band, sq)
        if (do_rec) call s_xa_rec(site, 1, nr, tq, peer=peer, key=key, seq=sq)
        w%nreq = w%nreq + 1
        call s_amr_wave_size_int(w%req, w%nreq); call s_amr_wave_size_int(w%reqw, w%nreq)
        w%reqw(w%nreq) = -1
#ifdef MFC_MPI
        if (dev) then
            #:call GPU_HOST_DATA(use_device_addr='[buf]')
                call MPI_ISEND(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
            #:endcall GPU_HOST_DATA
        else
            call MPI_ISEND(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
        end if
#endif

    end subroutine s_amr_wave_isend

    !> Host-array variants for the zero-copy waves, which send register sections (any rank) directly: sequence association on an
    !! assumed-size dummy, no device address.
    impure subroutine s_amr_wave_irecv_raw(w, buf, n, peer, site, key)

        type(t_amr_wave), intent(inout) :: w
        real(wp), intent(inout)         :: buf(*)
        integer, intent(in)             :: n, peer, site, key
        integer                         :: ierr, sq, tq

        sq = f_amr_m1_seq(peer, 2); tq = f_amr_m1_tag(w%band, sq)
        call s_xa_rec(site, 2, n, tq, peer=peer, key=key, seq=sq)
        w%nreq = w%nreq + 1
        call s_amr_wave_size_int(w%req, w%nreq); call s_amr_wave_size_int(w%reqw, w%nreq)
        w%reqw(w%nreq) = n
#ifdef MFC_MPI
        call MPI_IRECV(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
#endif

    end subroutine s_amr_wave_irecv_raw

    impure subroutine s_amr_wave_isend_raw(w, buf, n, peer, site, key)

        type(t_amr_wave), intent(inout) :: w
        real(wp), intent(in)            :: buf(*)
        integer, intent(in)             :: n, peer, site, key
        integer                         :: ierr, sq, tq

        sq = f_amr_m1_seq(peer, 1); tq = f_amr_m1_tag(w%band, sq)
        call s_xa_rec(site, 1, n, tq, peer=peer, key=key, seq=sq)
        w%nreq = w%nreq + 1
        call s_amr_wave_size_int(w%req, w%nreq); call s_amr_wave_size_int(w%reqw, w%nreq)
        w%reqw(w%nreq) = -1
#ifdef MFC_MPI
        call MPI_ISEND(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
#endif

    end subroutine s_amr_wave_isend_raw

    !> Wait for every request of the open wave; under MFC_DEBUG every receive's length must equal its plan (a short message means
    !! the two sides enumerated different transfers).
    impure subroutine s_amr_wave_wait(w)

        type(t_amr_wave), intent(inout) :: w
        integer                         :: ierr

        if (w%nreq == 0) return
#ifdef MFC_MPI
#ifdef MFC_DEBUG
        block
            integer :: st(MPI_STATUS_SIZE, w%nreq), gotw, q
            call MPI_WAITALL(w%nreq, w%req, st, ierr)
            do q = 1, w%nreq
                if (w%reqw(q) < 0) cycle
                call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                @:ASSERT(gotw == w%reqw(q), "amr wave: a received message length differs from the plan")
            end do
        end block
#else
        call MPI_WAITALL(w%nreq, w%req, MPI_STATUSES_IGNORE, ierr)
#endif
#endif
        w%nreq = 0

    end subroutine s_amr_wave_wait

    !> Start enumerating a side's transfers.
    impure subroutine s_amr_wave_reset(s)

        type(t_amr_wave_side), intent(inout) :: s

        s%nx = 0; s%np = 0; s%words = 0

    end subroutine s_amr_wave_reset

    !> Append one transfer of `words` payload words for block blk over box [bl, bh] with peer.
    impure subroutine s_amr_wave_add(s, peer, blk, bl, bh, words)

        type(t_amr_wave_side), intent(inout) :: s
        integer, intent(in)                  :: peer, blk, bl(3), bh(3), words

        s%nx = s%nx + 1
        call s_amr_wave_size_int(s%blk, s%nx); call s_amr_wave_size_int3(s%bl, s%nx); call s_amr_wave_size_int3(s%bh, s%nx)
        call s_amr_wave_size_int(s%cnt, s%nx); call s_amr_wave_size_int(s%peer, s%nx); call s_amr_wave_size_int(s%pi, s%nx)
        call s_amr_wave_size_int(s%off, s%nx)
        s%blk(s%nx) = blk; s%bl(:,s%nx) = bl; s%bh(:,s%nx) = bh; s%cnt(s%nx) = words; s%peer(s%nx) = peer

    end subroutine s_amr_wave_add

    !> Lay the side out in its pool: peers in first-appearance order, each peer's transfers contiguous in appended order with XA_NH
    !! header words in front of every transfer; grows pool (device-resident when dev) to the words in use.
    impure subroutine s_amr_wave_close(s, pool, dev)

        type(t_amr_wave_side), intent(inout) :: s
        real(wp), allocatable, intent(inout) :: pool(:)
        logical, intent(in)                  :: dev
        integer                              :: t, r, ip, base

        if (.not. allocated(bmap)) then
            allocate (bmap(0:num_procs - 1), bnx(0:num_procs - 1), bwords(0:num_procs - 1))
            bmap = 0; bnx = 0; bwords = 0
        end if
        s%np = 0
        do t = 1, s%nx
            r = s%peer(t)
            if (bmap(r) == 0) then
                s%np = s%np + 1
                call s_amr_wave_size_int(s%prank, s%np); call s_amr_wave_size_int(s%pwords, s%np)
                call s_amr_wave_size_int(s%pbase, s%np); call s_amr_wave_size_int(s%pnx, s%np)
                bmap(r) = s%np; s%prank(s%np) = r
            end if
            s%pi(t) = bmap(r)
            s%off(t) = bwords(r) + bnx(r)*XA_NH  ! relative to the peer's base
            bwords(r) = bwords(r) + s%cnt(t)
            bnx(r) = bnx(r) + 1
        end do
        base = 0
        do ip = 1, s%np
            r = s%prank(ip)
            s%pnx(ip) = bnx(r)
            s%pwords(ip) = bwords(r) + bnx(r)*XA_NH
            s%pbase(ip) = base; base = base + s%pwords(ip)
            bmap(r) = 0; bnx(r) = 0; bwords(r) = 0
        end do
        do t = 1, s%nx
            s%off(t) = s%pbase(s%pi(t)) + s%off(t)
        end do
        s%words = base
        call s_amr_wave_size_real(pool, base, dev)

    end subroutine s_amr_wave_close

    !> Pool bounds [lo, hi] of transfer t's payload (its XA_NH header words sit at off+1:off+XA_NH).
    pure subroutine s_amr_wave_slice(s, t, lo, hi)

        type(t_amr_wave_side), intent(in) :: s
        integer, intent(in)               :: t
        integer, intent(out)              :: lo, hi

        lo = s%off(t) + XA_NH + 1; hi = s%off(t) + XA_NH + s%cnt(t)

    end subroutine s_amr_wave_slice

    !> Debug identity header of transfer t (no-op when the audit is off).
    impure subroutine s_amr_wave_hdr_pack(s, pool, t, site)

        type(t_amr_wave_side), intent(in) :: s
        real(wp), intent(inout)           :: pool(:)
        integer, intent(in)               :: t, site

        if (XA_NH > 0) call s_xa_hdr_pack(pool(s%off(t) + 1:s%off(t) + XA_NH), site, s%blk(t), s%bl(:,t), s%bh(:,t))

    end subroutine s_amr_wave_hdr_pack

    impure subroutine s_amr_wave_hdr_check(s, pool, t, site)

        type(t_amr_wave_side), intent(in) :: s
        real(wp), intent(in)              :: pool(:)
        integer, intent(in)               :: t, site

        if (XA_NH > 0) call s_xa_hdr_check(pool(s%off(t) + 1:s%off(t) + XA_NH), site, s%blk(t), s%bl(:,t), s%bh(:,t))

    end subroutine s_amr_wave_hdr_check

    !> One receive per peer into the side's pool runs.
    impure subroutine s_amr_wave_post(w, s, pool, site, dev)

        type(t_amr_wave), intent(inout)     :: w
        type(t_amr_wave_side), intent(in)   :: s
        real(wp), intent(inout), contiguous :: pool(:)
        integer, intent(in)                 :: site
        logical, intent(in)                 :: dev
        integer                             :: ip

        do ip = 1, s%np
            call s_amr_wave_irecv(w, pool(s%pbase(ip) + 1:s%pbase(ip) + s%pwords(ip)), s%pwords(ip), s%prank(ip), site, &
                                  & s%pnx(ip), dev, nrec=s%pwords(ip) - s%pnx(ip)*XA_NH)
        end do

    end subroutine s_amr_wave_post

    !> One send per peer from the side's pool runs (the caller has packed every transfer's slice and header).
    impure subroutine s_amr_wave_send(w, s, pool, site, dev)

        type(t_amr_wave), intent(inout)   :: w
        type(t_amr_wave_side), intent(in) :: s
        real(wp), intent(in), contiguous  :: pool(:)
        integer, intent(in)               :: site
        logical, intent(in)               :: dev
        integer                           :: ip

        do ip = 1, s%np
            call s_amr_wave_isend(w, pool(s%pbase(ip) + 1:s%pbase(ip) + s%pwords(ip)), s%pwords(ip), s%prank(ip), site, &
                                  & s%pnx(ip), dev, nrec=s%pwords(ip) - s%pnx(ip)*XA_NH)
        end do

    end subroutine s_amr_wave_send

    !> High-water sizing for plan tables: a grow preserves the entries already appended.
    impure subroutine s_amr_wave_size_int(a, n)

        integer, allocatable, intent(inout) :: a(:)
        integer, intent(in)                 :: n
        integer, allocatable                :: tmp(:)

        if (.not. allocated(a)) then
            allocate (a(max(n, 64)))
            return
        end if
        if (size(a) >= n) return
        call move_alloc(a, tmp)
        allocate (a(max(n, 2*size(tmp))))
        a(1:size(tmp)) = tmp

    end subroutine s_amr_wave_size_int

    impure subroutine s_amr_wave_size_int3(a, n)

        integer, allocatable, intent(inout) :: a(:,:)
        integer, intent(in)                 :: n
        integer, allocatable                :: tmp(:,:)

        if (.not. allocated(a)) then
            allocate (a(3, max(n, 64)))
            return
        end if
        if (size(a, 2) >= n) return
        call move_alloc(a, tmp)
        allocate (a(3, max(n, 2*size(tmp, 2))))
        a(:,1:size(tmp, 2)) = tmp

    end subroutine s_amr_wave_size_int3

    !> Wire pools, preserving on grow. dev keeps the pool device-resident across (re)allocation: the old image is deleted from the
    !! device before it is freed and the new one created after; contents never survive a wave, so nothing is copied.
    impure subroutine s_amr_wave_size_real(a, n, dev)

        real(wp), allocatable, intent(inout) :: a(:)
        integer, intent(in)                  :: n
        logical, intent(in)                  :: dev
        real(wp), allocatable                :: tmp(:)

        if (.not. allocated(a)) then
            allocate (a(max(n, 64)))
            if (dev) then
                $:GPU_ENTER_DATA(create='[a]')
            end if
            return
        end if
        if (size(a) >= n) return
        if (dev) then
            $:GPU_EXIT_DATA(delete='[a]')
        end if
        call move_alloc(a, tmp)
        allocate (a(max(n, 2*size(tmp))))
        a(1:size(tmp)) = tmp
        if (dev) then
            $:GPU_ENTER_DATA(create='[a]')
        end if

    end subroutine s_amr_wave_size_real

end module m_amr_wave
