!>
!!@file
!!@brief Contains module m_amr_wave

#:include 'macros.fpp'

!> @brief MPI point-to-point waves for AMR exchange: post all receives, pack, post all sends, wait once. Pooled waves list transfers
!! on a send/receive side and s_amr_wave_close lays them out contiguously per peer; both ranks enumerate the same transfers in the
!! same order, so no metadata is exchanged. Pools are caller-owned (device-resident under rdma_mpi).
module m_amr_wave

#ifdef MFC_MPI
    use mpi
#endif
    use m_precision_select
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_amr_state, only: s_amr_size_int, s_amr_size_int3, s_amr_size_real, amr_tsq, amr_tsq_tch, amr_n_tsq, amr_tag_gen, &
        & amr_tag_base

    implicit none

    private
    public :: t_amr_wave_side, t_amr_wave, s_amr_wave_open, s_amr_wave_req, s_amr_wave_req_raw, s_amr_wave_wait, f_amr_wave_nreq, &
        & s_amr_wave_reset, s_amr_wave_add, s_amr_wave_add_slabs, s_amr_wave_close, s_amr_wave_post, s_amr_wave_send, &
        & s_amr_wave_slice

    !> One side (send or receive) of a pooled wave.
    type t_amr_wave_side
        integer :: nx = 0     !< transfers
        integer :: np = 0     !< peers
        integer :: words = 0  !< pool words in use after close
        !> per transfer (off: pool offset)
        integer, allocatable :: blk(:), bl(:,:), bh(:,:), cnt(:), peer(:), pi(:), off(:)
        integer, allocatable :: prank(:), pwords(:), pbase(:), pnx(:)  !< per peer: rank, message words, pool base, transfers
    end type t_amr_wave_side

    !> Fill and restrict waves share amr_wsend/amr_wrecv; the seam wave stays open across the fills, so it has its own pair.
    type(t_amr_wave_side), public :: amr_wsend, amr_wrecv, amr_wseam_snd, amr_wseam_rcv

    !> An open wave: tag band and requests; reqw is -1 for sends, the planned word count for receives.
    type t_amr_wave
        integer              :: band = 0, nreq = 0
        integer, allocatable :: req(:), reqw(:)
    end type t_amr_wave
    type(t_amr_wave), public, target :: amr_wave, amr_wave_seam
    !> Rank-indexed scratch for s_amr_wave_close; all-zero between calls.
    integer, allocatable :: bmap(:), bnx(:), bwords(:)

contains

    !> Open a wave on tag band b.
    impure subroutine s_amr_wave_open(w, b)

        type(t_amr_wave), intent(inout) :: w
        integer, intent(in)             :: b

        w%band = b; w%nreq = 0
        call s_amr_tag_wave_open(b)

    end subroutine s_amr_wave_open

    pure integer function f_amr_wave_nreq(w) result(n)

        type(t_amr_wave), intent(in) :: w

        n = w%nreq

    end function f_amr_wave_nreq

    !> Post one request of n words between buf (device-resident when dev) and peer: dir = 1 sends, 2 receives.
    impure subroutine s_amr_wave_req(w, dir, buf, n, peer, dev)

        type(t_amr_wave), intent(inout)     :: w
        integer, intent(in)                 :: dir, n, peer
        real(wp), intent(inout), contiguous :: buf(:)
        logical, intent(in)                 :: dev
        integer                             :: ierr, tq

        tq = f_amr_wave_tag(w, dir, peer)
        w%reqw(w%nreq) = merge(n, -1, dir == 2)
#ifdef MFC_MPI
        if (dev) then
            #:call GPU_HOST_DATA(use_device_addr='[buf]')
                if (dir == 2) then
                    call MPI_IRECV(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
                else
                    call MPI_ISEND(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
                end if
            #:endcall GPU_HOST_DATA
        else if (dir == 2) then
            call MPI_IRECV(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
        else
            call MPI_ISEND(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
        end if
#endif

    end subroutine s_amr_wave_req

    !> Host-array variant of s_amr_wave_req (assumed-size, any-rank sections, no device address).
    impure subroutine s_amr_wave_req_raw(w, dir, buf, n, peer)

        type(t_amr_wave), intent(inout) :: w
        integer, intent(in)             :: dir, n, peer
        real(wp), intent(inout)         :: buf(*)
        integer                         :: ierr, tq

        tq = f_amr_wave_tag(w, dir, peer)
        w%reqw(w%nreq) = merge(n, -1, dir == 2)
#ifdef MFC_MPI
        if (dir == 2) then
            call MPI_IRECV(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
        else
            call MPI_ISEND(buf, n, mpi_p, peer, tq, MPI_COMM_WORLD, w%req(w%nreq), ierr)
        end if
#endif

    end subroutine s_amr_wave_req_raw

    !> Open a keyed-tag wave: clear the touched per-peer seq counters and bump the band's generation. Call at wave entry on every
    !! rank.
    impure subroutine s_amr_tag_wave_open(band)

        integer, intent(in) :: band
        integer             :: i, d

        if (.not. allocated(amr_tsq)) then
            allocate (amr_tsq(0:num_procs - 1,2), amr_tsq_tch(num_procs, 2))
            amr_tsq = 0; amr_n_tsq = 0
        end if
        do d = 1, 2
            do i = 1, amr_n_tsq(d)
                amr_tsq(amr_tsq_tch(i, d), d) = 0
            end do
            amr_n_tsq(d) = 0
        end do
        amr_tag_gen(band) = mod(amr_tag_gen(band) + 1, 16)

    end subroutine s_amr_tag_wave_open

    !> Next seq on the (peer, dir) channel of the open wave. dir: 1 = send, 2 = recv (separate counters).
    impure integer function f_amr_wave_seq(peer, idir) result(sq)

        integer, intent(in) :: peer, idir

        if (amr_tsq(peer, idir) == 0) then
            amr_n_tsq(idir) = amr_n_tsq(idir) + 1
            amr_tsq_tch(amr_n_tsq(idir), idir) = peer
        end if
        amr_tsq(peer, idir) = amr_tsq(peer, idir) + 1
        sq = amr_tsq(peer, idir)
        @:ASSERT(sq < 4096, "keyed-tag seq overflows its 12-bit field")

    end function f_amr_wave_seq

    pure integer function f_amr_wave_tag_of(band, sq) result(t)

        integer, intent(in) :: band, sq

        t = amr_tag_base + band*65536 + amr_tag_gen(band)*4096 + sq

    end function f_amr_wave_tag_of

    !> Keyed tag of the next request to peer; grows the request list (w%nreq becomes the new request's index).
    impure integer function f_amr_wave_tag(w, dir, peer) result(tq)

        type(t_amr_wave), intent(inout) :: w
        integer, intent(in)             :: dir, peer
        integer                         :: sq

        sq = f_amr_wave_seq(peer, dir); tq = f_amr_wave_tag_of(w%band, sq)
        w%nreq = w%nreq + 1
        call s_amr_size_int(w%req, w%nreq); call s_amr_size_int(w%reqw, w%nreq)

    end function f_amr_wave_tag

    !> Wait for all requests; under MFC_DEBUG each receive length must equal its plan.
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
        call s_amr_size_int(s%blk, s%nx); call s_amr_size_int3(s%bl, s%nx); call s_amr_size_int3(s%bh, s%nx)
        call s_amr_size_int(s%cnt, s%nx); call s_amr_size_int(s%peer, s%nx); call s_amr_size_int(s%pi, s%nx)
        call s_amr_size_int(s%off, s%nx)
        s%blk(s%nx) = blk; s%bl(:,s%nx) = bl; s%bh(:,s%nx) = bh; s%cnt(s%nx) = words; s%peer(s%nx) = peer

    end subroutine s_amr_wave_add

    !> Add the msl slabs [tb:te] of block blk as transfers with peer.
    impure subroutine s_amr_wave_add_slabs(s, peer, blk, msl, tb, te)

        type(t_amr_wave_side), intent(inout) :: s
        integer, intent(in)                  :: peer, blk, msl, tb(3, 6), te(3, 6)
        integer                              :: isl

        do isl = 1, msl
            call s_amr_wave_add(s, peer, blk, tb(:,isl), te(:,isl), sys_size*product(te(:,isl) - tb(:,isl) + 1))
        end do

    end subroutine s_amr_wave_add_slabs

    !> Lay the side out in its pool (peers in first-appearance order, transfers contiguous per peer); grows pool to the words in
    !! use.
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
                call s_amr_size_int(s%prank, s%np); call s_amr_size_int(s%pwords, s%np)
                call s_amr_size_int(s%pbase, s%np); call s_amr_size_int(s%pnx, s%np)
                bmap(r) = s%np; s%prank(s%np) = r
            end if
            s%pi(t) = bmap(r)
            s%off(t) = bwords(r)  ! relative to the peer's base
            bwords(r) = bwords(r) + s%cnt(t)
            bnx(r) = bnx(r) + 1
        end do
        base = 0
        do ip = 1, s%np
            r = s%prank(ip)
            s%pnx(ip) = bnx(r)
            s%pwords(ip) = bwords(r)
            s%pbase(ip) = base; base = base + s%pwords(ip)
            bmap(r) = 0; bnx(r) = 0; bwords(r) = 0
        end do
        do t = 1, s%nx
            s%off(t) = s%pbase(s%pi(t)) + s%off(t)
        end do
        s%words = base
        call s_amr_size_real(pool, base, dev)

    end subroutine s_amr_wave_close

    !> Pool bounds [lo, hi] of transfer t's payload.
    pure subroutine s_amr_wave_slice(s, t, lo, hi)

        type(t_amr_wave_side), intent(in) :: s
        integer, intent(in)               :: t
        integer, intent(out)              :: lo, hi

        lo = s%off(t) + 1; hi = s%off(t) + s%cnt(t)

    end subroutine s_amr_wave_slice

    !> One message per peer (dir = 1 sends, 2 receives) over the side's pool runs.
    impure subroutine s_amr_wave_peers(w, dir, s, pool, dev)

        type(t_amr_wave), intent(inout)     :: w
        integer, intent(in)                 :: dir
        type(t_amr_wave_side), intent(in)   :: s
        real(wp), intent(inout), contiguous :: pool(:)
        logical, intent(in)                 :: dev
        integer                             :: ip

        do ip = 1, s%np
            call s_amr_wave_req(w, dir, pool(s%pbase(ip) + 1:s%pbase(ip) + s%pwords(ip)), s%pwords(ip), s%prank(ip), dev)
        end do

    end subroutine s_amr_wave_peers

    !> One receive per peer into the side's pool runs.
    impure subroutine s_amr_wave_post(w, s, pool, dev)

        type(t_amr_wave), intent(inout)     :: w
        type(t_amr_wave_side), intent(in)   :: s
        real(wp), intent(inout), contiguous :: pool(:)
        logical, intent(in)                 :: dev

        call s_amr_wave_peers(w, 2, s, pool, dev)

    end subroutine s_amr_wave_post

    !> One send per peer (the caller has packed every slice).
    impure subroutine s_amr_wave_send(w, s, pool, dev)

        type(t_amr_wave), intent(inout)     :: w
        type(t_amr_wave_side), intent(in)   :: s
        real(wp), intent(inout), contiguous :: pool(:)
        logical, intent(in)                 :: dev

        call s_amr_wave_peers(w, 1, s, pool, dev)

    end subroutine s_amr_wave_send

end module m_amr_wave
