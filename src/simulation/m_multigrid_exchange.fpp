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

    !> Each exchange sends one message per distinct neighbor rank (mg_nbr_rank), holding all the sides it shares with this one
    real(wp), allocatable, dimension(:) :: mg_send_buf, mg_recv_buf
    $:GPU_DECLARE(create='[mg_send_buf, mg_recv_buf]')
    !> Each side's action (mg_side_none, mg_side_pack or mg_side_self) and its offset in mg_send_buf and mg_recv_buf by level;
    !! message i goes to rank mg_msg_rank(i), from mg_send_beg + 1 for mg_send_len values, and arrives at mg_recv_beg + 1 for
    !! mg_recv_len
    integer                          :: mg_num_msg
    integer, dimension(6)            :: mg_side_mode, mg_msg_rank
    integer, dimension(6, mg_maxlev) :: mg_send_off, mg_recv_off
    integer, dimension(6, mg_maxlev) :: mg_send_beg, mg_send_len, mg_recv_beg, mg_recv_len
    $:GPU_DECLARE(create='[mg_side_mode, mg_send_off, mg_recv_off]')

contains

    !> Buffers and per-level message layout: sides go out in this rank's side order and arrive in the neighbor's, so both agree
    impure subroutine s_initialize_multigrid_exchange()

        integer :: lv, side, side_nb, i, tot

        tot = 2*((n + 1)*(p + 1) + (m + 1)*(p + 1) + (m + 1)*(n + 1))
        @:ALLOCATE(mg_send_buf(tot), mg_recv_buf(tot))
        mg_num_msg = 0
        do side = 1, 6
            mg_side_mode(side) = mg_side_none
            if (mg_nbr_rank(side) < 0) cycle
            mg_side_mode(side) = merge(mg_side_self, mg_side_pack, mg_nbr_rank(side) == proc_rank)
            if (mg_side_mode(side) == mg_side_self .or. any(mg_msg_rank(1:mg_num_msg) == mg_nbr_rank(side))) cycle
            mg_num_msg = mg_num_msg + 1
            mg_msg_rank(mg_num_msg) = mg_nbr_rank(side)
        end do
        mg_send_off = 0
        mg_recv_off = 0
        do lv = 1, mg_num_lev
            tot = 0
            do i = 1, mg_num_msg
                mg_send_beg(i, lv) = tot
                do side = 1, 6
                    if (mg_nbr_rank(side) /= mg_msg_rank(i)) cycle
                    mg_send_off(side, lv) = tot
                    tot = tot + f_side_size(lv, side)
                end do
                mg_send_len(i, lv) = tot - mg_send_beg(i, lv)
            end do
            tot = 0
            do i = 1, mg_num_msg
                mg_recv_beg(i, lv) = tot
                do side_nb = 1, 6
                    side = f_opposite(side_nb)
                    if (mg_nbr_rank(side) /= mg_msg_rank(i)) cycle
                    mg_recv_off(side, lv) = tot
                    tot = tot + f_side_size(lv, side)
                end do
                mg_recv_len(i, lv) = tot - mg_recv_beg(i, lv)
            end do
        end do
        $:GPU_UPDATE(device='[mg_side_mode, mg_send_off, mg_recv_off]')

    end subroutine s_initialize_multigrid_exchange

    !> Fill the ghost layer of v, laid out as level lv, from the neighbors: one message per distinct neighbor rank holding all the
    !! sides it shares with this one, and periodic seams a rank owns alone copied in place
    impure subroutine s_mg_exchange(lv, v)

        integer, intent(in) :: lv

        real(wp), dimension(-mg_ghost_x:mg_nx(lv) - 1 + mg_ghost_x,-mg_ghost_y:mg_ny(lv) - 1 + mg_ghost_y, &
             & -mg_ghost_z:mg_nz(lv) - 1 + mg_ghost_z), intent(inout) :: v
        integer :: nt

        call s_mg_sides(lv, v, unpack=.false.)
        if (mg_num_msg == 0) return
#ifdef MFC_MPI
        if (rdma_mpi) then
            $:GPU_WAIT()
            #:call GPU_HOST_DATA(use_device_addr='[mg_send_buf, mg_recv_buf]')
                call s_mg_sendrecv(lv, mg_send_buf, mg_recv_buf)
            #:endcall GPU_HOST_DATA
        else
            nt = mg_send_beg(mg_num_msg, lv) + mg_send_len(mg_num_msg, lv)
            $:GPU_UPDATE(host='[mg_send_buf(1:nt)]')
            call s_mg_sendrecv(lv, mg_send_buf, mg_recv_buf)
            nt = mg_recv_beg(mg_num_msg, lv) + mg_recv_len(mg_num_msg, lv)
            $:GPU_UPDATE(device='[mg_recv_buf(1:nt)]')
        end if
#endif
        call s_mg_sides(lv, v, unpack=.true.)

    end subroutine s_mg_exchange

    !> Post level lv's messages between send_buf and recv_buf (device addresses under rdma_mpi) and wait for them
    impure subroutine s_mg_sendrecv(lv, send_buf, recv_buf)

        integer, intent(in)                   :: lv
        real(wp), dimension(*), intent(in)    :: send_buf
        real(wp), dimension(*), intent(inout) :: recv_buf
        integer, dimension(12)                :: req
        integer                               :: i, ierr

#ifdef MFC_MPI
        do i = 1, mg_num_msg
            call MPI_IRECV(recv_buf(mg_recv_beg(i, lv) + 1), mg_recv_len(i, lv), mpi_p, mg_msg_rank(i), 7101, MPI_COMM_WORLD, &
                           & req(i), ierr)
            call MPI_ISEND(send_buf(mg_send_beg(i, lv) + 1), mg_send_len(i, lv), mpi_p, mg_msg_rank(i), 7101, MPI_COMM_WORLD, &
                           & req(mg_num_msg + i), ierr)
        end do
        call MPI_WAITALL(2*mg_num_msg, req, MPI_STATUSES_IGNORE, ierr)
#endif

    end subroutine s_mg_sendrecv

    !> The side at the opposite end of the same direction
    pure function f_opposite(side) result(side_opp)

        integer, intent(in) :: side
        integer             :: side_opp

        side_opp = side + merge(1, -1, mod(side, 2) == 1)

    end function f_opposite

    !> Cells in the layer of a side on level lv
    pure function f_side_size(lv, side) result(cells)

        integer, intent(in) :: lv, side
        integer             :: cells

        select case ((side + 1)/2)
        case (1)
            cells = mg_ny(lv)*mg_nz(lv)
        case (2)
            cells = mg_nx(lv)*mg_nz(lv)
        case default
            cells = mg_nx(lv)*mg_ny(lv)
        end select

    end function f_side_size

    !> All six sides of v on level lv in one kernel: pack and copy self seams, or with unpack, fill ghosts from mg_recv_buf
    impure subroutine s_mg_sides(lv, v, unpack)

        integer, intent(in) :: lv
        logical, intent(in) :: unpack

        real(wp), dimension(-mg_ghost_x:mg_nx(lv) - 1 + mg_ghost_x,-mg_ghost_y:mg_ny(lv) - 1 + mg_ghost_y, &
             & -mg_ghost_z:mg_nz(lv) - 1 + mg_ghost_z), intent(inout) :: v
        integer :: nx, ny, nz, nmax, side, dir, a, b, mode, n_t1, n_t2, n_norm, layer, ghost, src, bidx
        integer :: j_g, k_g, l_g, j_s, k_s, l_s
        logical :: low

        nx = mg_nx(lv)
        ny = mg_ny(lv)
        nz = mg_nz(lv)
        nmax = max(nx, ny, nz)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[side, dir, low, a, b, mode, n_t1, n_t2, n_norm, layer, ghost, src, bidx, j_g, &
                            & k_g, l_g, j_s, k_s, l_s]')
        do side = 1, 6
            do b = 0, nmax - 1
                do a = 0, nmax - 1
                    mode = mg_side_mode(side)
                    if (unpack) mode = merge(mg_side_unpack, mg_side_none, mode == mg_side_pack)
                    dir = (side + 1)/2
                    low = mod(side, 2) == 1
                    if (dir == 1) then
                        n_norm = nx
                        n_t1 = ny
                        n_t2 = nz
                    else if (dir == 2) then
                        n_norm = ny
                        n_t1 = nx
                        n_t2 = nz
                    else
                        n_norm = nz
                        n_t1 = nx
                        n_t2 = ny
                    end if
                    if (mode /= mg_side_none .and. a < n_t1 .and. b < n_t2) then
                        layer = merge(0, n_norm - 1, low)
                        ghost = merge(-1, n_norm, low)
                        ! A pack reads the side's own layer, a self seam the opposite one
                        src = merge(layer, n_norm - 1 - layer, mode == mg_side_pack)
                        if (dir == 1) then
                            j_g = ghost
                            k_g = a
                            l_g = b
                            j_s = src
                            k_s = a
                            l_s = b
                        else if (dir == 2) then
                            j_g = a
                            k_g = ghost
                            l_g = b
                            j_s = a
                            k_s = src
                            l_s = b
                        else
                            j_g = a
                            k_g = b
                            l_g = ghost
                            j_s = a
                            k_s = b
                            l_s = src
                        end if
                        bidx = b*n_t1 + a + 1
                        if (mode == mg_side_pack) then
                            mg_send_buf(mg_send_off(side, lv) + bidx) = v(j_s, k_s, l_s)
                        else if (mode == mg_side_unpack) then
                            v(j_g, k_g, l_g) = mg_recv_buf(mg_recv_off(side, lv) + bidx)
                        else
                            v(j_g, k_g, l_g) = v(j_s, k_s, l_s)
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_sides

    impure subroutine s_finalize_multigrid_exchange()

        @:DEALLOCATE(mg_send_buf, mg_recv_buf)

    end subroutine s_finalize_multigrid_exchange

end module m_multigrid_exchange
