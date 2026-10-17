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
    !> Each side's action (0 none, 1 pack for another rank, 3 periodic self copy) and its offset in mg_send_buf (1) and mg_recv_buf
    !! (2) by level; message i goes to rank mg_msg_rank(i), from mg_send_beg + 1 for mg_send_len values, and arrives at mg_recv_beg
    !! + 1 for mg_recv_len
    integer                             :: mg_num_msg
    integer, dimension(6)               :: mg_smode, mg_msg_rank
    integer, dimension(6, mg_maxlev, 2) :: mg_buf_off
    integer, dimension(6, mg_maxlev)    :: mg_send_beg, mg_send_len, mg_recv_beg, mg_recv_len
    $:GPU_DECLARE(create='[mg_smode, mg_buf_off]')

contains

    !> Buffers and the message layout of every level, fixed by the neighbor ranks: this rank sends its sides in its own side order,
    !! and receives each neighbor's sides in that neighbor's side order, so both ends agree without further communication
    impure subroutine s_initialize_multigrid_exchange()

        integer :: lv, side, side_nb, i, tot

        tot = 2*((n + 1)*(p + 1) + (m + 1)*(p + 1) + (m + 1)*(n + 1))
        @:ALLOCATE(mg_send_buf(tot), mg_recv_buf(tot))
        @:PIN_HOST(mg_send_buf, mg_recv_buf)
        mg_num_msg = 0
        do side = 1, 6
            mg_smode(side) = 0
            if (mg_nbr_rank(side) < 0) cycle
            mg_smode(side) = merge(3, 1, mg_nbr_rank(side) == proc_rank)
            if (mg_smode(side) == 3 .or. any(mg_msg_rank(1:mg_num_msg) == mg_nbr_rank(side))) cycle
            mg_num_msg = mg_num_msg + 1
            mg_msg_rank(mg_num_msg) = mg_nbr_rank(side)
        end do
        mg_buf_off = 0
        do lv = 1, mg_num_lev
            tot = 0
            do i = 1, mg_num_msg
                mg_send_beg(i, lv) = tot
                do side = 1, 6
                    if (mg_nbr_rank(side) /= mg_msg_rank(i)) cycle
                    mg_buf_off(side, lv, 1) = tot
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
                    mg_buf_off(side, lv, 2) = tot
                    tot = tot + f_side_size(lv, side)
                end do
                mg_recv_len(i, lv) = tot - mg_recv_beg(i, lv)
            end do
        end do
        $:GPU_UPDATE(device='[mg_smode, mg_buf_off]')

    end subroutine s_initialize_multigrid_exchange

    !> Fill the ghost layer of v, laid out as level lv, from the neighbors: one message per distinct neighbor rank holding all the
    !! sides it shares with this one, and periodic seams a rank owns alone copied in place
    impure subroutine s_mg_exchange(lv, v)

        integer, intent(in) :: lv

        real(wp), dimension(-mg_ghost_x:mg_nx(lv) - 1 + mg_ghost_x,-mg_ghost_y:mg_ny(lv) - 1 + mg_ghost_y, &
             & -mg_ghost_z:mg_nz(lv) - 1 + mg_ghost_z), intent(inout) :: v
        integer, dimension(12) :: req
        integer                :: i, nt, ierr

        call s_mg_sides(lv, v, 1)
        if (mg_num_msg == 0) return
#ifdef MFC_MPI
        if (rdma_mpi) then
            $:GPU_WAIT()
            #:call GPU_HOST_DATA(use_device_addr='[mg_send_buf, mg_recv_buf]')
                do i = 1, mg_num_msg
                    call MPI_IRECV(mg_recv_buf(mg_recv_beg(i, lv) + 1), mg_recv_len(i, lv), mpi_p, mg_msg_rank(i), 7101, &
                                   & MPI_COMM_WORLD, req(i), ierr)
                    call MPI_ISEND(mg_send_buf(mg_send_beg(i, lv) + 1), mg_send_len(i, lv), mpi_p, mg_msg_rank(i), 7101, &
                                   & MPI_COMM_WORLD, req(mg_num_msg + i), ierr)
                end do
                call MPI_WAITALL(2*mg_num_msg, req, MPI_STATUSES_IGNORE, ierr)
            #:endcall GPU_HOST_DATA
        else
            nt = mg_send_beg(mg_num_msg, lv) + mg_send_len(mg_num_msg, lv)
            $:GPU_UPDATE(host='[mg_send_buf(1:nt)]')
            do i = 1, mg_num_msg
                call MPI_IRECV(mg_recv_buf(mg_recv_beg(i, lv) + 1), mg_recv_len(i, lv), mpi_p, mg_msg_rank(i), 7101, &
                               & MPI_COMM_WORLD, req(i), ierr)
                call MPI_ISEND(mg_send_buf(mg_send_beg(i, lv) + 1), mg_send_len(i, lv), mpi_p, mg_msg_rank(i), 7101, &
                               & MPI_COMM_WORLD, req(mg_num_msg + i), ierr)
            end do
            call MPI_WAITALL(2*mg_num_msg, req, MPI_STATUSES_IGNORE, ierr)
            nt = mg_recv_beg(mg_num_msg, lv) + mg_recv_len(mg_num_msg, lv)
            $:GPU_UPDATE(device='[mg_recv_buf(1:nt)]')
        end if
#endif
        call s_mg_sides(lv, v, 2)

    end subroutine s_mg_exchange

    !> The side at the opposite end of the same direction
    pure integer function f_opposite(side)

        integer, intent(in) :: side

        f_opposite = side + merge(1, -1, mod(side, 2) == 1)

    end function f_opposite

    !> Cells in the layer of a side on level lv
    pure integer function f_side_size(lv, side)

        integer, intent(in) :: lv, side

        select case ((side + 1)/2)
        case (1)
            f_side_size = mg_ny(lv)*mg_nz(lv)
        case (2)
            f_side_size = mg_nx(lv)*mg_nz(lv)
        case default
            f_side_size = mg_nx(lv)*mg_ny(lv)
        end select

    end function f_side_size

    !> All six sides of v, laid out as level lv, in one kernel. Phase 1 copies each periodic seam a rank owns alone into its ghost
    !! layer and packs each side facing another rank into mg_send_buf; phase 2 unpacks mg_recv_buf into those ghosts
    impure subroutine s_mg_sides(lv, v, phase)

        integer, intent(in) :: lv, phase

        real(wp), dimension(-mg_ghost_x:mg_nx(lv) - 1 + mg_ghost_x,-mg_ghost_y:mg_ny(lv) - 1 + mg_ghost_y, &
             & -mg_ghost_z:mg_nz(lv) - 1 + mg_ghost_z), intent(inout) :: v
        integer :: nx, ny, nz, nmax, side, a, b, mode, n_t1, n_t2, n_norm, layer, ghost, src, bidx, lv_loc, phase_loc
        integer :: j_g, k_g, l_g, j_s, k_s, l_s

        ! Dummies may alias host variables, which a device kernel must not reference

        lv_loc = lv
        phase_loc = phase
        nx = mg_nx(lv)
        ny = mg_ny(lv)
        nz = mg_nz(lv)
        nmax = max(nx, ny, nz)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[side, a, b, mode, n_t1, n_t2, n_norm, layer, ghost, src, bidx, j_g, k_g, l_g, &
                            & j_s, k_s, l_s]')
        do side = 1, 6
            do b = 0, nmax - 1
                do a = 0, nmax - 1
                    ! Phase 1 acts on self seams (3) and packs (1), phase 2 unpacks (becomes 2)
                    mode = mg_smode(side)
                    if (phase_loc == 2) mode = merge(2, 0, mode == 1)
                    if ((side + 1)/2 == 1) then
                        n_norm = nx
                        n_t1 = ny
                        n_t2 = nz
                    else if ((side + 1)/2 == 2) then
                        n_norm = ny
                        n_t1 = nx
                        n_t2 = nz
                    else
                        n_norm = nz
                        n_t1 = nx
                        n_t2 = ny
                    end if
                    if (mode > 0 .and. a < n_t1 .and. b < n_t2) then
                        layer = merge(0, n_norm - 1, mod(side, 2) == 1)
                        ghost = merge(-1, n_norm, mod(side, 2) == 1)
                        ! A pack reads the side's own layer, a self seam the opposite one
                        src = merge(layer, n_norm - 1 - layer, mode == 1)
                        if ((side + 1)/2 == 1) then
                            j_g = ghost
                            k_g = a
                            l_g = b
                            j_s = src
                            k_s = a
                            l_s = b
                        else if ((side + 1)/2 == 2) then
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
                        bidx = mg_buf_off(side, lv_loc, merge(1, 2, mode == 1)) + b*n_t1 + a + 1
                        if (mode == 1) then
                            mg_send_buf(bidx) = v(j_s, k_s, l_s)
                        else if (mode == 2) then
                            v(j_g, k_g, l_g) = mg_recv_buf(bidx)
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

        @:UNPIN_HOST(mg_send_buf, mg_recv_buf)
        @:DEALLOCATE(mg_send_buf, mg_recv_buf)

    end subroutine s_finalize_multigrid_exchange

end module m_multigrid_exchange
