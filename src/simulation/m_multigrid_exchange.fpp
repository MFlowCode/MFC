!>
!! @file
!! @brief Contains module m_multigrid_exchange

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Ghost-layer exchange of the multigrid levels: periodic seams a rank owns alone are copied in place, and the sides facing
!! other ranks travel in one aggregated message per neighbor rank
module m_multigrid_exchange

#ifdef MFC_MPI
    use mpi
#endif
    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_multigrid_levels

    implicit none

    private; public :: s_initialize_multigrid_exchange, s_mg_exchange, s_finalize_multigrid_exchange

    !> Each exchange sends one message per distinct neighbor rank (mg_nbr), holding all the sides it shares with this one
    real(wp), allocatable, dimension(:) :: mg_sbuf, mg_rbuf
    $:GPU_DECLARE(create='[mg_sbuf, mg_rbuf]')
    !> Side q's action (0 none, 1 pack for another rank, 3 periodic self copy) and its offset in mg_sbuf (1) and mg_rbuf (2) by
    !! level; message i goes to rank mg_nlist(i), from mg_sbeg + 1 for mg_slen values, and arrives at mg_rbeg + 1 for mg_rlen
    integer                             :: mg_nmsg
    integer, dimension(6)               :: mg_smode, mg_nlist
    integer, dimension(6, mg_maxlev, 2) :: mg_boff
    integer, dimension(6, mg_maxlev)    :: mg_sbeg, mg_slen, mg_rbeg, mg_rlen
    $:GPU_DECLARE(create='[mg_smode, mg_boff]')

contains

    !> Buffers and the message layout of every level, fixed by the neighbor ranks: this rank sends its sides in its own side order,
    !! and receives each neighbor's sides in that neighbor's side order, so both ends agree without further communication
    impure subroutine s_initialize_multigrid_exchange()

        integer :: lv, q, qq, i, tot

        tot = 2*((n + 1)*(p + 1) + (m + 1)*(p + 1) + (m + 1)*(n + 1))
        @:ALLOCATE(mg_sbuf(tot), mg_rbuf(tot))
        @:PIN_HOST(mg_sbuf, mg_rbuf)
        mg_nmsg = 0
        do q = 1, 6
            mg_smode(q) = 0
            if (mg_nbr(q) < 0) cycle
            mg_smode(q) = merge(3, 1, mg_nbr(q) == proc_rank)
            if (mg_smode(q) == 3 .or. any(mg_nlist(1:mg_nmsg) == mg_nbr(q))) cycle
            mg_nmsg = mg_nmsg + 1; mg_nlist(mg_nmsg) = mg_nbr(q)
        end do
        mg_boff = 0
        do lv = 1, mg_nlev
            tot = 0
            do i = 1, mg_nmsg
                mg_sbeg(i, lv) = tot
                do q = 1, 6
                    if (mg_nbr(q) /= mg_nlist(i)) cycle
                    mg_boff(q, lv, 1) = tot; tot = tot + f_side_size(lv, q)
                end do
                mg_slen(i, lv) = tot - mg_sbeg(i, lv)
            end do
            tot = 0
            do i = 1, mg_nmsg
                mg_rbeg(i, lv) = tot
                do qq = 1, 6
                    q = f_opposite(qq)
                    if (mg_nbr(q) /= mg_nlist(i)) cycle
                    mg_boff(q, lv, 2) = tot; tot = tot + f_side_size(lv, q)
                end do
                mg_rlen(i, lv) = tot - mg_rbeg(i, lv)
            end do
        end do
        $:GPU_UPDATE(device='[mg_smode, mg_boff]')

    end subroutine s_initialize_multigrid_exchange

    !> Fill the ghost layer of the block of mg_e at offset off, laid out as level lv, from the neighbors: one message per distinct
    !! neighbor rank holding all the sides it shares with this one, and periodic seams a rank owns alone copied in place
    impure subroutine s_mg_exchange(lv, off)

        integer, intent(in)    :: lv, off
        integer, dimension(12) :: req
        integer                :: i, nt, ierr

        call s_mg_sides(lv, off, 1)
        if (mg_nmsg == 0) return
#ifdef MFC_MPI
        if (rdma_mpi) then
            $:GPU_WAIT()
            #:call GPU_HOST_DATA(use_device_addr='[mg_sbuf, mg_rbuf]')
                do i = 1, mg_nmsg
                    call MPI_IRECV(mg_rbuf(mg_rbeg(i, lv) + 1), mg_rlen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, req(i), &
                                   & ierr)
                    call MPI_ISEND(mg_sbuf(mg_sbeg(i, lv) + 1), mg_slen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, &
                                   & req(mg_nmsg + i), ierr)
                end do
                call MPI_WAITALL(2*mg_nmsg, req, MPI_STATUSES_IGNORE, ierr)
            #:endcall GPU_HOST_DATA
        else
            nt = mg_sbeg(mg_nmsg, lv) + mg_slen(mg_nmsg, lv)
            $:GPU_UPDATE(host='[mg_sbuf(1:nt)]')
            do i = 1, mg_nmsg
                call MPI_IRECV(mg_rbuf(mg_rbeg(i, lv) + 1), mg_rlen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, req(i), ierr)
                call MPI_ISEND(mg_sbuf(mg_sbeg(i, lv) + 1), mg_slen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, &
                               & req(mg_nmsg + i), ierr)
            end do
            call MPI_WAITALL(2*mg_nmsg, req, MPI_STATUSES_IGNORE, ierr)
            nt = mg_rbeg(mg_nmsg, lv) + mg_rlen(mg_nmsg, lv)
            $:GPU_UPDATE(device='[mg_rbuf(1:nt)]')
        end if
#endif
        call s_mg_sides(lv, off, 2)

    end subroutine s_mg_exchange

    !> Side q of the opposite end of the same direction
    pure integer function f_opposite(q)

        integer, intent(in) :: q

        f_opposite = q + merge(1, -1, mod(q, 2) == 1)

    end function f_opposite

    !> Cells in the layer of side q on level lv
    pure integer function f_side_size(lv, q)

        integer, intent(in) :: lv, q

        select case ((q + 1)/2)
        case (1); f_side_size = mg_ny(lv)*mg_nz(lv)
        case (2); f_side_size = mg_nx(lv)*mg_nz(lv)
        case default; f_side_size = mg_nx(lv)*mg_ny(lv)
        end select

    end function f_side_size

    !> All six sides of the block at offset off laid out as level lv, in one kernel. Phase 1 copies each periodic seam a rank owns
    !! alone into its ghost layer and packs each side facing another rank into mg_sbuf; phase 2 unpacks mg_rbuf into those ghosts
    impure subroutine s_mg_sides(lv, off, phase)

        integer, intent(in) :: lv, off, phase
        integer             :: nx, ny, nz, ex, ey, base, nmax, q, a, b, md, n1, n2, s1, s2, st, nn, lay, gho, idx, bidx, lvl, ph

        ! Dummies may alias host variables, which a device kernel must not reference

        lvl = lv; ph = phase
        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); ex = nx + 2*mg_gx; ey = ny + 2*mg_gy
        base = off + (mg_gz*ey + mg_gy)*ex + mg_gx + 1
        nmax = max(nx, ny, nz)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[q, a, b, md, n1, n2, s1, s2, st, nn, lay, gho, idx, bidx]')
        do q = 1, 6
            do b = 0, nmax - 1
                do a = 0, nmax - 1
                    ! Phase 1 acts on self seams (3) and packs (1), phase 2 unpacks (becomes 2)
                    md = mg_smode(q)
                    if (ph == 2) md = merge(2, 0, md == 1)
                    if ((q + 1)/2 == 1) then
                        nn = nx; st = 1; n1 = ny; s1 = ex; n2 = nz; s2 = ex*ey
                    else if ((q + 1)/2 == 2) then
                        nn = ny; st = ex; n1 = nx; s1 = 1; n2 = nz; s2 = ex*ey
                    else
                        nn = nz; st = ex*ey; n1 = nx; s1 = 1; n2 = ny; s2 = ex
                    end if
                    if (md > 0 .and. a < n1 .and. b < n2) then
                        lay = merge(0, nn - 1, mod(q, 2) == 1)
                        gho = merge(-1, nn, mod(q, 2) == 1)
                        idx = base + a*s1 + b*s2
                        bidx = mg_boff(q, lvl, merge(1, 2, md == 1)) + b*n1 + a + 1
                        if (md == 1) then
                            mg_sbuf(bidx) = mg_e(idx + lay*st)
                        else if (md == 2) then
                            mg_e(idx + gho*st) = mg_rbuf(bidx)
                        else
                            mg_e(idx + gho*st) = mg_e(idx + (nn - 1 - lay)*st)
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_sides

    impure subroutine s_finalize_multigrid_exchange()

        @:UNPIN_HOST(mg_sbuf, mg_rbuf)
        @:DEALLOCATE(mg_sbuf, mg_rbuf)

    end subroutine s_finalize_multigrid_exchange

end module m_multigrid_exchange
