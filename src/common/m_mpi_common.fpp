
#:include 'macros.fpp'

!> @brief The module serves as a proxy to the parameters and subroutines
!!          available in the MPI implementation's MPI module. Specifically,
!!          the purpose of the proxy is to harness basic MPI commands into
!!          more complicated procedures as to accomplish the communication
!!          goals for the simulation.
module m_mpi_common

#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_helper

    use ieee_arithmetic

    use m_nvtx

    implicit none

    integer, private :: v_size
    $:GPU_DECLARE(create='[v_size]')
    !! Generic flags used to identify and report MPI errors

    real(wp), private, allocatable, dimension(:) :: buff_send !<
    !! This variable is utilized to pack and send the buffer of the cell-average
    !! primitive variables, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    real(wp), private, allocatable, dimension(:) :: buff_recv !<
    !! buff_recv is utilized to receive and unpack the buffer of the cell-
    !! average primitive variables, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    $:GPU_DECLARE(create='[buff_send, buff_recv]')

    integer :: halo_size
    $:GPU_DECLARE(create='[halo_size]')

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    impure subroutine s_initialize_mpi_common_module

#ifdef MFC_MPI
        ! Allocating buff_send/recv and. Please note that for the sake of
        ! simplicity, both variables are provided sufficient storage to hold
        ! the largest buffer in the computational domain.

        if (qbmm .and. .not. polytropic) then
            v_size = sys_size + 2*nb*4
        else
            v_size = sys_size
        end if

        if (n > 0) then
            if (p > 0) then
                halo_size = nint(-1._wp + 1._wp*buff_size*(v_size)* &
                                         & (m + 2*buff_size + 1)* &
                                         & (n + 2*buff_size + 1)* &
                                         & (p + 2*buff_size + 1)/ &
                                         & (cells_bounds%mnp_min + 2*buff_size + 1))
            else
                halo_size = -1 + buff_size*(v_size)* &
                                         & (cells_bounds%mn_max + 2*buff_size + 1)
            end if
        else
            halo_size = -1 + buff_size*(v_size)
        end if

        $:GPU_UPDATE(device='[halo_size, v_size]')

        @:ALLOCATE(buff_send(0:halo_size), buff_recv(0:halo_size))
#endif

    end subroutine s_initialize_mpi_common_module

    !> The subroutine initializes the MPI execution environment
        !!      and queries both the number of processors which will be
        !!      available for the job and the local processor rank.
    impure subroutine s_mpi_initialize

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors
#endif

#ifndef MFC_MPI

        ! Serial run only has 1 processor
        num_procs = 1
        ! Local processor rank is 0
        proc_rank = 0

#else

        ! Initializing the MPI environment
        call MPI_INIT(ierr)

        ! Checking whether the MPI environment has been properly initialized
        if (ierr /= MPI_SUCCESS) then
            print '(A)', 'Unable to initialize MPI environment. Exiting.'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Querying the number of processors available for the job
        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

        ! Querying the rank of the local processor
        call MPI_COMM_RANK(MPI_COMM_WORLD, proc_rank, ierr)

#endif

    end subroutine s_mpi_initialize

    !! @param q_cons_vf Conservative variables
    !! @param ib_markers track if a cell is within the immersed boundary
    !! @param levelset closest distance from every cell to the IB
    !! @param levelset_norm normalized vector from every cell to the closest point to the IB
    !! @param beta Eulerian void fraction from lagrangian bubbles
    impure subroutine s_initialize_mpi_data(q_cons_vf, ib_markers, levelset, levelset_norm, beta)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(integer_field), optional, intent(in) :: ib_markers
        type(levelset_field), optional, intent(IN) :: levelset
        type(levelset_norm_field), optional, intent(IN) :: levelset_norm
        type(scalar_field), intent(in), optional :: beta

        integer, dimension(num_dims) :: sizes_glb, sizes_loc
        integer, dimension(1) :: airfoil_glb, airfoil_loc, airfoil_start

#ifdef MFC_MPI

        ! Generic loop iterator
        integer :: i, j
        integer :: ierr !< Generic flag used to identify and report MPI errors

        !Altered system size for the lagrangian subgrid bubble model
        integer :: alt_sys

        if (present(beta)) then
            alt_sys = sys_size + 1
        else
            alt_sys = sys_size
        end if

        do i = 1, sys_size
            MPI_IO_DATA%var(i)%sf => q_cons_vf(i)%sf(0:m, 0:n, 0:p)
        end do

        if (present(beta)) then
            MPI_IO_DATA%var(alt_sys)%sf => beta%sf(0:m, 0:n, 0:p)
        end if

        !Additional variables pb and mv for non-polytropic qbmm
#ifdef MFC_PRE_PROCESS
        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do j = 1, nnode
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j)%sf => pb%sf(0:m, 0:n, 0:p, j, i)
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j + nb*nnode)%sf => mv%sf(0:m, 0:n, 0:p, j, i)
                end do
            end do
        end if
#endif

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do j = 1, nnode
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j)%sf => pb_ts(1)%sf(0:m, 0:n, 0:p, j, i)
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j + nb*nnode)%sf => mv_ts(1)%sf(0:m, 0:n, 0:p, j, i)
                end do
            end do
        end if
#endif
        ! Define global(g) and local(l) sizes for flow variables
        sizes_glb(1) = m_glb + 1; sizes_loc(1) = m + 1
        if (n > 0) then
            sizes_glb(2) = n_glb + 1; sizes_loc(2) = n + 1
            if (p > 0) then
                sizes_glb(3) = p_glb + 1; sizes_loc(3) = p + 1
            end if
        end if

        ! Define the view for each variable
        do i = 1, alt_sys
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                          MPI_ORDER_FORTRAN, mpi_p, MPI_IO_DATA%view(i), ierr)
            call MPI_TYPE_COMMIT(MPI_IO_DATA%view(i), ierr)
        end do

#ifndef MFC_POST_PROCESS
        if (qbmm .and. .not. polytropic) then
            do i = sys_size + 1, sys_size + 2*nb*4
                call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                              MPI_ORDER_FORTRAN, mpi_p, MPI_IO_DATA%view(i), ierr)
                call MPI_TYPE_COMMIT(MPI_IO_DATA%view(i), ierr)

            end do
        end if
#endif

        if (present(ib_markers)) then

#ifdef MFC_PRE_PROCESS
            MPI_IO_IB_DATA%var%sf => ib_markers%sf
            MPI_IO_levelset_DATA%var%sf => levelset%sf
            MPI_IO_levelsetnorm_DATA%var%sf => levelset_norm%sf
#else
            MPI_IO_IB_DATA%var%sf => ib_markers%sf(0:m, 0:n, 0:p)

#ifndef MFC_POST_PROCESS
            MPI_IO_levelset_DATA%var%sf => levelset%sf(0:m, 0:n, 0:p, 1:num_ibs)
            MPI_IO_levelsetnorm_DATA%var%sf => levelset_norm%sf(0:m, 0:n, 0:p, 1:num_ibs, 1:3)
#endif

#endif
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                          MPI_ORDER_FORTRAN, MPI_INTEGER, MPI_IO_IB_DATA%view, ierr)
            call MPI_TYPE_COMMIT(MPI_IO_IB_DATA%view, ierr)

#ifndef MFC_POST_PROCESS
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                          MPI_ORDER_FORTRAN, mpi_p, MPI_IO_levelset_DATA%view, ierr)
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                          MPI_ORDER_FORTRAN, mpi_p, MPI_IO_levelsetnorm_DATA%view, ierr)

            call MPI_TYPE_COMMIT(MPI_IO_levelset_DATA%view, ierr)
            call MPI_TYPE_COMMIT(MPI_IO_levelsetnorm_DATA%view, ierr)
#endif
        end if

#ifndef MFC_POST_PROCESS
        if (present(ib_markers)) then
            do j = 1, num_ibs
                if (patch_ib(j)%c > 0) then

#ifdef MFC_PRE_PROCESS
                    allocate (MPI_IO_airfoil_IB_DATA%var(1:2*Np))
#endif

                    airfoil_glb(1) = 3*Np*num_procs
                    airfoil_loc(1) = 3*Np
                    airfoil_start(1) = 3*proc_rank*Np

#ifdef MFC_PRE_PROCESS
                    do i = 1, Np
                        MPI_IO_airfoil_IB_DATA%var(i)%x = airfoil_grid_l(i)%x
                        MPI_IO_airfoil_IB_DATA%var(i)%y = airfoil_grid_l(i)%y
                    end do
#endif

                    call MPI_TYPE_CREATE_SUBARRAY(1, airfoil_glb, airfoil_loc, airfoil_start, &
                                                  MPI_ORDER_FORTRAN, mpi_p, MPI_IO_airfoil_IB_DATA%view(1), ierr)
                    call MPI_TYPE_COMMIT(MPI_IO_airfoil_IB_DATA%view(1), ierr)

#ifdef MFC_PRE_PROCESS
                    do i = 1, Np
                        MPI_IO_airfoil_IB_DATA%var(Np + i)%x = airfoil_grid_u(i)%x
                        MPI_IO_airfoil_IB_DATA%var(Np + i)%y = airfoil_grid_u(i)%y
                    end do
#endif
                    call MPI_TYPE_CREATE_SUBARRAY(1, airfoil_glb, airfoil_loc, airfoil_start, &
                                                  MPI_ORDER_FORTRAN, mpi_p, MPI_IO_airfoil_IB_DATA%view(2), ierr)
                    call MPI_TYPE_COMMIT(MPI_IO_airfoil_IB_DATA%view(2), ierr)

                end if
            end do

        end if
#endif

#endif

    end subroutine s_initialize_mpi_data

    impure subroutine s_mpi_gather_data(my_vector, counts, gathered_vector, root)

        integer, intent(in) :: counts          ! Array of vector lengths for each process
        real(wp), intent(in), dimension(counts) :: my_vector   ! Input vector on each process
        integer, intent(in) :: root               ! Rank of the root process
        real(wp), allocatable, intent(out) :: gathered_vector(:) ! Gathered vector on the root process

        integer :: i
        integer :: ierr !< Generic flag used to identify and report MPI errors
        integer, allocatable :: recounts(:), displs(:)

#ifdef MFC_MPI

        allocate (recounts(num_procs))

        call MPI_GATHER(counts, 1, MPI_INTEGER, recounts, 1, MPI_INTEGER, root, &
                        MPI_COMM_WORLD, ierr)

        allocate (displs(size(recounts)))

        displs(1) = 0

        do i = 2, size(recounts)
            displs(i) = displs(i - 1) + recounts(i - 1)
        end do

        allocate (gathered_vector(sum(recounts)))
        call MPI_GATHERV(my_vector, counts, mpi_p, gathered_vector, recounts, displs, mpi_p, &
                         root, MPI_COMM_WORLD, ierr)
#endif
    end subroutine s_mpi_gather_data

    impure subroutine mpi_bcast_time_step_values(proc_time, time_avg)

        real(wp), dimension(0:num_procs - 1), intent(inout) :: proc_time
        real(wp), intent(inout) :: time_avg

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        call MPI_GATHER(time_avg, 1, mpi_p, proc_time(0), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

#endif

    end subroutine mpi_bcast_time_step_values

    impure subroutine s_prohibit_abort(condition, message)
        character(len=*), intent(in) :: condition, message

        print *, ""
        print *, "CASE FILE ERROR"
        print *, "  - Prohibited condition: ", trim(condition)
        if (len_trim(message) > 0) then
            print *, "  - Note: ", trim(message)
        end if
        print *, ""
        call s_mpi_abort(code=CASE_FILE_ERROR_CODE)
    end subroutine s_prohibit_abort

    !>  The goal of this subroutine is to determine the global
        !!      extrema of the stability criteria in the computational
        !!      domain. This is performed by sifting through the local
        !!      extrema of each stability criterion. Note that each of
        !!      the local extrema is from a single process, within its
        !!      assigned section of the computational domain. Finally,
        !!      note that the global extrema values are only bookkeept
        !!      on the rank 0 processor.
        !!  @param icfl_max_loc Local maximum ICFL stability criterion
        !!  @param vcfl_max_loc Local maximum VCFL stability criterion
        !!  @param Rc_min_loc Local minimum Rc stability criterion
        !!  @param icfl_max_glb Global maximum ICFL stability criterion
        !!  @param vcfl_max_glb Global maximum VCFL stability criterion
        !!  @param Rc_min_glb Global minimum Rc stability criterion
    impure subroutine s_mpi_reduce_stability_criteria_extrema(icfl_max_loc, &
                                                              vcfl_max_loc, &
                                                              Rc_min_loc, &
                                                              icfl_max_glb, &
                                                              vcfl_max_glb, &
                                                              Rc_min_glb)

        real(wp), intent(in) :: icfl_max_loc
        real(wp), intent(in) :: vcfl_max_loc
        real(wp), intent(in) :: Rc_min_loc

        real(wp), intent(out) :: icfl_max_glb
        real(wp), intent(out) :: vcfl_max_glb
        real(wp), intent(out) :: Rc_min_glb

#ifdef MFC_SIMULATION
#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Reducing local extrema of ICFL, VCFL, CCFL and Rc numbers to their
        ! global extrema and bookkeeping the results on the rank 0 processor
        call MPI_REDUCE(icfl_max_loc, icfl_max_glb, 1, &
                        mpi_p, MPI_MAX, 0, &
                        MPI_COMM_WORLD, ierr)

        if (viscous) then
            call MPI_REDUCE(vcfl_max_loc, vcfl_max_glb, 1, &
                            mpi_p, MPI_MAX, 0, &
                            MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(Rc_min_loc, Rc_min_glb, 1, &
                            mpi_p, MPI_MIN, 0, &
                            MPI_COMM_WORLD, ierr)
        end if

#else

        icfl_max_glb = icfl_max_loc

        if (viscous) then
            vcfl_max_glb = vcfl_max_loc
            Rc_min_glb = Rc_min_loc
        end if

#endif
#endif

    end subroutine s_mpi_reduce_stability_criteria_extrema

    !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the sum of all
        !!      values. The reduced variable is recorded back onto the
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
    impure subroutine s_mpi_allreduce_sum(var_loc, var_glb)

        real(wp), intent(in) :: var_loc
        real(wp), intent(out) :: var_glb

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, mpi_p, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_sum

    !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the minimum of all
        !!      values. The reduced variable is recorded back onto the
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
    impure subroutine s_mpi_allreduce_min(var_loc, var_glb)

        real(wp), intent(in) :: var_loc
        real(wp), intent(out) :: var_glb

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, mpi_p, &
                           MPI_MIN, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_min

    !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the maximum of all
        !!      values. The reduced variable is recorded back onto the
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
    impure subroutine s_mpi_allreduce_max(var_loc, var_glb)

        real(wp), intent(in) :: var_loc
        real(wp), intent(out) :: var_glb

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, mpi_p, &
                           MPI_MAX, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_max

    !>  The following subroutine takes the inputted variable and
        !!      determines its minimum value on the entire computational
        !!      domain. The result is stored back into inputted variable.
        !!  @param var_loc holds the local value to be reduced among
        !!      all the processors in communicator. On output, the variable holds
        !!      the minimum value, reduced amongst all of the local values.
    impure subroutine s_mpi_reduce_min(var_loc)

        real(wp), intent(inout) :: var_loc

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Temporary storage variable that holds the reduced minimum value
        real(wp) :: var_glb

        ! Performing reduction procedure and eventually storing its result
        ! into the variable that was initially inputted into the subroutine
        call MPI_REDUCE(var_loc, var_glb, 1, mpi_p, &
                        MPI_MIN, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(var_glb, 1, mpi_p, &
                       0, MPI_COMM_WORLD, ierr)

        var_loc = var_glb

#endif

    end subroutine s_mpi_reduce_min

    !>  The following subroutine takes the first element of the
        !!      2-element inputted variable and determines its maximum
        !!      value on the entire computational domain. The result is
        !!      stored back into the first element of the variable while
        !!      the rank of the processor that is in charge of the sub-
        !!      domain containing the maximum is stored into the second
        !!      element of the variable.
        !!  @param var_loc On input, this variable holds the local value and processor rank,
        !!  which are to be reduced among all the processors in communicator.
        !!  On output, this variable holds the maximum value, reduced amongst
        !!  all of the local values, and the process rank to which the value
        !!  belongs.
    impure subroutine s_mpi_reduce_maxloc(var_loc)

        real(wp), dimension(2), intent(inout) :: var_loc

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        real(wp), dimension(2) :: var_glb  !<
            !! Temporary storage variable that holds the reduced maximum value
            !! and the rank of the processor with which the value is associated

        ! Performing reduction procedure and eventually storing its result
        ! into the variable that was initially inputted into the subroutine
        call MPI_REDUCE(var_loc, var_glb, 1, mpi_2p, &
                        MPI_MAXLOC, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(var_glb, 1, mpi_2p, &
                       0, MPI_COMM_WORLD, ierr)

        var_loc = var_glb

#endif

    end subroutine s_mpi_reduce_maxloc

    !> The subroutine terminates the MPI execution environment.
        !! @param prnt error message to be printed
    impure subroutine s_mpi_abort(prnt, code)

        character(len=*), intent(in), optional :: prnt
        integer, intent(in), optional :: code

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors
#endif

        if (present(prnt)) then
            print *, prnt
            call flush (6)

        end if

#ifndef MFC_MPI
        if (present(code)) then
            stop code
        else
            stop 1
        end if
#else
        ! Terminating the MPI environment
        if (present(code)) then
            call MPI_ABORT(MPI_COMM_WORLD, code, ierr)
        else
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if
#endif

    end subroutine s_mpi_abort

    !>Halts all processes until all have reached barrier.
    impure subroutine s_mpi_barrier

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Calling MPI_BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_barrier

    !> The subroutine finalizes the MPI execution environment.
    impure subroutine s_mpi_finalize

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Finalizing the MPI environment
        call MPI_FINALIZE(ierr)

#endif

    end subroutine s_mpi_finalize

    !>  The goal of this procedure is to populate the buffers of
        !!      the cell-average conservative variables by communicating
        !!      with the neighboring processors.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
    subroutine s_mpi_sendrecv_variables_buffers(q_comm, &
                                                mpi_dir, &
                                                pbc_loc, &
                                                nVar, &
                                                pb_in, mv_in)

        type(scalar_field), dimension(1:), intent(inout) :: q_comm
        real(wp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb_in, mv_in
        integer, intent(in) :: mpi_dir, pbc_loc, nVar

        integer :: i, j, k, l, r, q !< Generic loop iterators

        integer :: buffer_counts(1:3), buffer_count

        type(int_bounds_info) :: boundary_conditions(1:3)
        integer :: beg_end(1:2), grid_dims(1:3)
        integer :: dst_proc, src_proc, recv_tag, send_tag

        logical :: beg_end_geq_0, qbmm_comm

        integer :: pack_offset, unpack_offset

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        call nvtxStartRange("RHS-COMM-PACKBUF")

        qbmm_comm = .false.

        if (present(pb_in) .and. present(mv_in) .and. qbmm .and. .not. polytropic) then
            qbmm_comm = .true.
            v_size = nVar + 2*nb*4
            buffer_counts = (/ &
                            buff_size*v_size*(n + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                            /)
        else
            v_size = nVar
            buffer_counts = (/ &
                            buff_size*v_size*(n + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                            /)
        end if

        $:GPU_UPDATE(device='[v_size]')

        buffer_count = buffer_counts(mpi_dir)
        boundary_conditions = (/bc_x, bc_y, bc_z/)
        beg_end = (/boundary_conditions(mpi_dir)%beg, boundary_conditions(mpi_dir)%end/)
        beg_end_geq_0 = beg_end(max(pbc_loc, 0) - pbc_loc + 1) >= 0

        ! Implements:
        ! pbc_loc  bc_x >= 0 -> [send/recv]_tag  [dst/src]_proc
        ! -1 (=0)      0            ->     [1,0]       [0,0]      | 0 0 [1,0] [beg,beg]
        ! -1 (=0)      1            ->     [0,0]       [1,0]      | 0 1 [0,0] [end,beg]
        ! +1 (=1)      0            ->     [0,1]       [1,1]      | 1 0 [0,1] [end,end]
        ! +1 (=1)      1            ->     [1,1]       [0,1]      | 1 1 [1,1] [beg,end]

        send_tag = f_logical_to_int(.not. f_xor(beg_end_geq_0, pbc_loc == 1))
        recv_tag = f_logical_to_int(pbc_loc == 1)

        dst_proc = beg_end(1 + f_logical_to_int(f_xor(pbc_loc == 1, beg_end_geq_0)))
        src_proc = beg_end(1 + f_logical_to_int(pbc_loc == 1))

        grid_dims = (/m, n, p/)

        pack_offset = 0
        if (f_xor(pbc_loc == 1, beg_end_geq_0)) then
            pack_offset = grid_dims(mpi_dir) - buff_size + 1
        end if

        unpack_offset = 0
        if (pbc_loc == 1) then
            unpack_offset = grid_dims(mpi_dir) + buff_size + 1
        end if

        ! Pack Buffer to Send
        #:for mpi_dir in [1, 2, 3]
            if (mpi_dir == ${mpi_dir}$) then
                #:if mpi_dir == 1
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[r]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, buff_size - 1
                                do i = 1, nVar
                                    r = (i - 1) + v_size*(j + buff_size*(k + (n + 1)*l))
                                    buff_send(r) = q_comm(i)%sf(j + pack_offset, k, l)
                                end do
                            end do
                        end do
                    end do

                    if (qbmm_comm) then
                        $:GPU_PARALLEL_LOOP(collapse=4,private='[r]')
                        do l = 0, p
                            do k = 0, n
                                do j = 0, buff_size - 1
                                    do i = nVar + 1, nVar + 4
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                (j + buff_size*(k + (n + 1)*l))
                                            buff_send(r) = pb_in(j + pack_offset, k, l, i - nVar, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do l = 0, p
                            do k = 0, n
                                do j = 0, buff_size - 1
                                    do i = nVar + 1, nVar + 4
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                (j + buff_size*(k + (n + 1)*l))
                                            buff_send(r) = mv_in(j + pack_offset, k, l, i - nVar, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:elif mpi_dir == 2
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[r]')
                    do i = 1, nVar
                        do l = 0, p
                            do k = 0, buff_size - 1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         (k + buff_size*l))
                                    buff_send(r) = q_comm(i)%sf(j, k + pack_offset, l)
                                end do
                            end do
                        end do
                    end do

                    if (qbmm_comm) then
                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do i = nVar + 1, nVar + 4
                            do l = 0, p
                                do k = 0, buff_size - 1
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 (k + buff_size*l))
                                            buff_send(r) = pb_in(j, k + pack_offset, l, i - nVar, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do i = nVar + 1, nVar + 4
                            do l = 0, p
                                do k = 0, buff_size - 1
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 (k + buff_size*l))
                                            buff_send(r) = mv_in(j, k + pack_offset, l, i - nVar, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:else
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[r]')
                    do i = 1, nVar
                        do l = 0, buff_size - 1
                            do k = -buff_size, n + buff_size
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                    buff_send(r) = q_comm(i)%sf(j, k, l + pack_offset)
                                end do
                            end do
                        end do
                    end do

                    if (qbmm_comm) then
                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do i = nVar + 1, nVar + 4
                            do l = 0, buff_size - 1
                                do k = -buff_size, n + buff_size
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                            buff_send(r) = pb_in(j, k, l + pack_offset, i - nVar, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do i = nVar + 1, nVar + 4
                            do l = 0, buff_size - 1
                                do k = -buff_size, n + buff_size
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                            buff_send(r) = mv_in(j, k, l + pack_offset, i - nVar, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:endif
            end if
        #:endfor
        call nvtxEndRange ! Packbuf

        ! Send/Recv
#ifdef MFC_SIMULATION
        #:for rdma_mpi in [False, True]
            if (rdma_mpi .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
                #:if rdma_mpi
                    #:call GPU_HOST_DATA(use_device='[buff_send, buff_recv]')
                        call nvtxStartRange("RHS-COMM-SENDRECV-RDMA")

                        call MPI_SENDRECV( &
                            buff_send, buffer_count, mpi_p, dst_proc, send_tag, &
                            buff_recv, buffer_count, mpi_p, src_proc, recv_tag, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                        call nvtxEndRange ! RHS-MPI-SENDRECV-(NO)-RDMA

                    #:endcall GPU_HOST_DATA
                    $:GPU_WAIT()
                #:else
                    call nvtxStartRange("RHS-COMM-DEV2HOST")
                    $:GPU_UPDATE(host='[buff_send]')
                    call nvtxEndRange
                    call nvtxStartRange("RHS-COMM-SENDRECV-NO-RMDA")

                    call MPI_SENDRECV( &
                        buff_send, buffer_count, mpi_p, dst_proc, send_tag, &
                        buff_recv, buffer_count, mpi_p, src_proc, recv_tag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    call nvtxEndRange ! RHS-MPI-SENDRECV-(NO)-RDMA

                    call nvtxStartRange("RHS-COMM-HOST2DEV")
                    $:GPU_UPDATE(device='[buff_recv]')
                    call nvtxEndRange
                #:endif
            end if
        #:endfor
#else
        call MPI_SENDRECV( &
            buff_send, buffer_count, mpi_p, dst_proc, send_tag, &
            buff_recv, buffer_count, mpi_p, src_proc, recv_tag, &
            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif

        ! Unpack Received Buffer
        call nvtxStartRange("RHS-COMM-UNPACKBUF")
        #:for mpi_dir in [1, 2, 3]
            if (mpi_dir == ${mpi_dir}$) then
                #:if mpi_dir == 1
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[r]')
                    do l = 0, p
                        do k = 0, n
                            do j = -buff_size, -1
                                do i = 1, nVar
                                    r = (i - 1) + v_size* &
                                        (j + buff_size*((k + 1) + (n + 1)*l))
                                    q_comm(i)%sf(j + unpack_offset, k, l) = buff_recv(r)
#if defined(__INTEL_COMPILER)
                                    if (ieee_is_nan(q_comm(i)%sf(j, k, l))) then
                                        print *, "Error", j, k, l, i
                                        error stop "NaN(s) in recv"
                                    end if
#endif
                                end do
                            end do
                        end do
                    end do

                    if (qbmm_comm) then
                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do l = 0, p
                            do k = 0, n
                                do j = -buff_size, -1
                                    do i = nVar + 1, nVar + 4
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                (j + buff_size*((k + 1) + (n + 1)*l))
                                            pb_in(j + unpack_offset, k, l, i - nVar, q) = buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do l = 0, p
                            do k = 0, n
                                do j = -buff_size, -1
                                    do i = nVar + 1, nVar + 4
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                (j + buff_size*((k + 1) + (n + 1)*l))
                                            mv_in(j + unpack_offset, k, l, i - nVar, q) = buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:elif mpi_dir == 2
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[r]')
                    do i = 1, nVar
                        do l = 0, p
                            do k = -buff_size, -1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + buff_size*l))
                                    q_comm(i)%sf(j, k + unpack_offset, l) = buff_recv(r)
#if defined(__INTEL_COMPILER)
                                    if (ieee_is_nan(q_comm(i)%sf(j, k, l))) then
                                        print *, "Error", j, k, l, i
                                        error stop "NaN(s) in recv"
                                    end if
#endif
                                end do
                            end do
                        end do
                    end do

                    if (qbmm_comm) then
                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do i = nVar + 1, nVar + 4
                            do l = 0, p
                                do k = -buff_size, -1
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + buff_size*l))
                                            pb_in(j, k + unpack_offset, l, i - nVar, q) = buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do i = nVar + 1, nVar + 4
                            do l = 0, p
                                do k = -buff_size, -1
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + buff_size*l))
                                            mv_in(j, k + unpack_offset, l, i - nVar, q) = buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:else
                    ! Unpacking buffer from bc_z%beg
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[r]')
                    do i = 1, nVar
                        do l = -buff_size, -1
                            do k = -buff_size, n + buff_size
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + (n + 2*buff_size + 1)* &
                                          (l + buff_size)))
                                    q_comm(i)%sf(j, k, l + unpack_offset) = buff_recv(r)
#if defined(__INTEL_COMPILER)
                                    if (ieee_is_nan(q_comm(i)%sf(j, k, l))) then
                                        print *, "Error", j, k, l, i
                                        error stop "NaN(s) in recv"
                                    end if
#endif
                                end do
                            end do
                        end do
                    end do

                    if (qbmm_comm) then
                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do i = nVar + 1, nVar + 4
                            do l = -buff_size, -1
                                do k = -buff_size, n + buff_size
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + (n + 2*buff_size + 1)* &
                                                  (l + buff_size)))
                                            pb_in(j, k, l + unpack_offset, i - nVar, q) = buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        $:GPU_PARALLEL_LOOP(collapse=5,private='[r]')
                        do i = nVar + 1, nVar + 4
                            do l = -buff_size, -1
                                do k = -buff_size, n + buff_size
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + (n + 2*buff_size + 1)* &
                                                  (l + buff_size)))
                                            mv_in(j, k, l + unpack_offset, i - nVar, q) = buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:endif
            end if
        #:endfor
        call nvtxEndRange
#endif

    end subroutine s_mpi_sendrecv_variables_buffers

    !>  The purpose of this procedure is to optimally decompose
        !!      the computational domain among the available processors.
        !!      This is performed by attempting to award each processor,
        !!      in each of the coordinate directions, approximately the
        !!      same number of cells, and then recomputing the affected
        !!      global parameters.
    subroutine s_mpi_decompose_computational_domain

#ifdef MFC_MPI

        integer :: num_procs_x, num_procs_y, num_procs_z !<
            !! Optimal number of processors in the x-, y- and z-directions

        real(wp) :: tmp_num_procs_x, tmp_num_procs_y, tmp_num_procs_z !<
            !! Non-optimal number of processors in the x-, y- and z-directions

        real(wp) :: fct_min !<
            !! Processor factorization (fct) minimization parameter

        integer :: MPI_COMM_CART !<
            !! Cartesian processor topology communicator

        integer :: rem_cells !<
            !! Remaining number of cells, in a particular coordinate direction,
            !! after the majority is divided up among the available processors

        integer :: recon_order !<
            !! WENO or MUSCL reconstruction order

        integer :: i, j !< Generic loop iterators
        integer :: ierr !< Generic flag used to identify and report MPI errors

        if (recon_type == WENO_TYPE) then
            recon_order = weno_order
        else
            recon_order = muscl_order
        end if

        if (num_procs == 1 .and. parallel_io) then
            do i = 1, num_dims
                start_idx(i) = 0
            end do
            return
        end if

        if (igr) then
            recon_order = igr_order
        else
            recon_order = weno_order
        end if

        ! 3D Cartesian Processor Topology
        if (n > 0) then

            if (p > 0) then

                if (cyl_coord .and. p > 0) then
                    ! Implement pencil processor blocking if using cylindrical coordinates so
                    ! that all cells in azimuthal direction are stored on a single processor.
                    ! This is necessary for efficient application of Fourier filter near axis.

                    ! Initial values of the processor factorization optimization
                    num_procs_x = 1
                    num_procs_y = num_procs
                    num_procs_z = 1
                    ierr = -1

                    ! Computing minimization variable for these initial values
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    tmp_num_procs_z = num_procs_z
                    fct_min = 10._wp*abs((m + 1)/tmp_num_procs_x &
                                         - (n + 1)/tmp_num_procs_y)

                    ! Searching for optimal computational domain distribution
                    do i = 1, num_procs

                        if (mod(num_procs, i) == 0 &
                            .and. &
                            (m + 1)/i >= num_stcls_min*recon_order) then

                            tmp_num_procs_x = i
                            tmp_num_procs_y = num_procs/i

                            if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                               - (n + 1)/tmp_num_procs_y) &
                                .and. &
                                (n + 1)/tmp_num_procs_y &
                                >= &
                                num_stcls_min*recon_order) then

                                num_procs_x = i
                                num_procs_y = num_procs/i
                                fct_min = abs((m + 1)/tmp_num_procs_x &
                                              - (n + 1)/tmp_num_procs_y)
                                ierr = 0

                            end if

                        end if

                    end do

                else

                    ! Initial estimate of optimal processor topology
                    num_procs_x = 1
                    num_procs_y = 1
                    num_procs_z = num_procs
                    ierr = -1

                    ! Benchmarking the quality of this initial guess
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    tmp_num_procs_z = num_procs_z
                    fct_min = 10._wp*abs((m + 1)/tmp_num_procs_x &
                                         - (n + 1)/tmp_num_procs_y) &
                              + 10._wp*abs((n + 1)/tmp_num_procs_y &
                                           - (p + 1)/tmp_num_procs_z)

                    ! Optimization of the initial processor topology
                    do i = 1, num_procs

                        if (mod(num_procs, i) == 0 &
                            .and. &
                            (m + 1)/i >= num_stcls_min*recon_order) then

                            do j = 1, num_procs/i

                                if (mod(num_procs/i, j) == 0 &
                                    .and. &
                                    (n + 1)/j >= num_stcls_min*recon_order) then

                                    tmp_num_procs_x = i
                                    tmp_num_procs_y = j
                                    tmp_num_procs_z = num_procs/(i*j)

                                    if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                                       - (n + 1)/tmp_num_procs_y) &
                                        + abs((n + 1)/tmp_num_procs_y &
                                              - (p + 1)/tmp_num_procs_z) &
                                        .and. &
                                        (p + 1)/tmp_num_procs_z &
                                        >= &
                                        num_stcls_min*recon_order) &
                                        then

                                        num_procs_x = i
                                        num_procs_y = j
                                        num_procs_z = num_procs/(i*j)
                                        fct_min = abs((m + 1)/tmp_num_procs_x &
                                                      - (n + 1)/tmp_num_procs_y) &
                                                  + abs((n + 1)/tmp_num_procs_y &
                                                        - (p + 1)/tmp_num_procs_z)
                                        ierr = 0

                                    end if

                                end if

                            end do

                        end if

                    end do

                end if

                ! Verifying that a valid decomposition of the computational
                ! domain has been established. If not, the simulation exits.
                if (proc_rank == 0 .and. ierr == -1) then
                    call s_mpi_abort('Unsupported combination of values '// &
                                     'of num_procs, m, n, p and '// &
                                     'weno/muscl/igr_order. Exiting.')
                end if

                ! Creating new communicator using the Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 3, (/num_procs_x, &
                                                          num_procs_y, num_procs_z/), &
                                     (/.true., .true., .true./), &
                                     .false., MPI_COMM_CART, ierr)

                ! Finding the Cartesian coordinates of the local process
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 3, &
                                     proc_coords, ierr)
                ! END: 3D Cartesian Processor Topology

                ! Global Parameters for z-direction

                ! Number of remaining cells
                rem_cells = mod(p + 1, num_procs_z)

                ! Optimal number of cells per processor
                p = (p + 1)/num_procs_z - 1

                ! Distributing the remaining cells
                do i = 1, rem_cells
                    if (proc_coords(3) == i - 1) then
                        p = p + 1; exit
                    end if
                end do

                ! Boundary condition at the beginning
                if (proc_coords(3) > 0 .or. (bc_z%beg == BC_PERIODIC .and. num_procs_z > 1)) then
                    proc_coords(3) = proc_coords(3) - 1
                    call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                       bc_z%beg, ierr)
                    proc_coords(3) = proc_coords(3) + 1
                end if

                ! Boundary condition at the end
                if (proc_coords(3) < num_procs_z - 1 .or. (bc_z%end == BC_PERIODIC .and. num_procs_z > 1)) then
                    proc_coords(3) = proc_coords(3) + 1
                    call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                       bc_z%end, ierr)
                    proc_coords(3) = proc_coords(3) - 1
                end if

#ifdef MFC_POST_PROCESS
                ! Ghost zone at the beginning
                if (proc_coords(3) > 0 .and. format == 1) then
                    offset_z%beg = 2
                else
                    offset_z%beg = 0
                end if

                ! Ghost zone at the end
                if (proc_coords(3) < num_procs_z - 1 .and. format == 1) then
                    offset_z%end = 2
                else
                    offset_z%end = 0
                end if
#endif

                ! Beginning and end sub-domain boundary locations
                if (parallel_io) then
                    if (proc_coords(3) < rem_cells) then
                        start_idx(3) = (p + 1)*proc_coords(3)
                    else
                        start_idx(3) = (p + 1)*proc_coords(3) + rem_cells
                    end if
                else
#ifdef MFC_PRE_PROCESS
                    if (old_grid .neqv. .true.) then
                        dz = (z_domain%end - z_domain%beg)/real(p_glb + 1, wp)

                        if (proc_coords(3) < rem_cells) then
                            z_domain%beg = z_domain%beg + dz*real((p + 1)* &
                                                                  proc_coords(3))
                            z_domain%end = z_domain%end - dz*real((p + 1)* &
                                                                  (num_procs_z - proc_coords(3) - 1) &
                                                                  - (num_procs_z - rem_cells))
                        else
                            z_domain%beg = z_domain%beg + dz*real((p + 1)* &
                                                                  proc_coords(3) + rem_cells)
                            z_domain%end = z_domain%end - dz*real((p + 1)* &
                                                                  (num_procs_z - proc_coords(3) - 1))
                        end if
                    end if
#endif
                end if

                ! 2D Cartesian Processor Topology
            else

                ! Initial estimate of optimal processor topology
                num_procs_x = 1
                num_procs_y = num_procs
                ierr = -1

                ! Benchmarking the quality of this initial guess
                tmp_num_procs_x = num_procs_x
                tmp_num_procs_y = num_procs_y
                fct_min = 10._wp*abs((m + 1)/tmp_num_procs_x &
                                     - (n + 1)/tmp_num_procs_y)

                ! Optimization of the initial processor topology
                do i = 1, num_procs

                    if (mod(num_procs, i) == 0 &
                        .and. &
                        (m + 1)/i >= num_stcls_min*recon_order) then

                        tmp_num_procs_x = i
                        tmp_num_procs_y = num_procs/i

                        if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                           - (n + 1)/tmp_num_procs_y) &
                            .and. &
                            (n + 1)/tmp_num_procs_y &
                            >= &
                            num_stcls_min*recon_order) then

                            num_procs_x = i
                            num_procs_y = num_procs/i
                            fct_min = abs((m + 1)/tmp_num_procs_x &
                                          - (n + 1)/tmp_num_procs_y)
                            ierr = 0

                        end if

                    end if

                end do

                ! Verifying that a valid decomposition of the computational
                ! domain has been established. If not, the simulation exits.
                if (proc_rank == 0 .and. ierr == -1) then
                    call s_mpi_abort('Unsupported combination of values '// &
                                     'of num_procs, m, n and '// &
                                     'weno/muscl/igr_order. Exiting.')
                end if

                ! Creating new communicator using the Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 2, (/num_procs_x, &
                                                          num_procs_y/), (/.true., &
                                                                           .true./), .false., MPI_COMM_CART, &
                                     ierr)

                ! Finding the Cartesian coordinates of the local process
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, &
                                     proc_coords, ierr)

            end if
            ! END: 2D Cartesian Processor Topology

            ! Global Parameters for y-direction

            ! Number of remaining cells
            rem_cells = mod(n + 1, num_procs_y)

            ! Optimal number of cells per processor
            n = (n + 1)/num_procs_y - 1

            ! Distributing the remaining cells
            do i = 1, rem_cells
                if (proc_coords(2) == i - 1) then
                    n = n + 1; exit
                end if
            end do

            ! Boundary condition at the beginning
            if (proc_coords(2) > 0 .or. (bc_y%beg == BC_PERIODIC .and. num_procs_y > 1)) then
                proc_coords(2) = proc_coords(2) - 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%beg, ierr)
                proc_coords(2) = proc_coords(2) + 1
            end if

            ! Boundary condition at the end
            if (proc_coords(2) < num_procs_y - 1 .or. (bc_y%end == BC_PERIODIC .and. num_procs_y > 1)) then
                proc_coords(2) = proc_coords(2) + 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%end, ierr)
                proc_coords(2) = proc_coords(2) - 1
            end if

#ifdef MFC_POST_PROCESS
            ! Ghost zone at the beginning
            if (proc_coords(2) > 0 .and. format == 1) then
                offset_y%beg = 2
            else
                offset_y%beg = 0
            end if

            ! Ghost zone at the end
            if (proc_coords(2) < num_procs_y - 1 .and. format == 1) then
                offset_y%end = 2
            else
                offset_y%end = 0
            end if
#endif

            ! Beginning and end sub-domain boundary locations
            if (parallel_io) then
                if (proc_coords(2) < rem_cells) then
                    start_idx(2) = (n + 1)*proc_coords(2)
                else
                    start_idx(2) = (n + 1)*proc_coords(2) + rem_cells
                end if
            else
#ifdef MFC_PRE_PROCESS
                if (old_grid .neqv. .true.) then
                    dy = (y_domain%end - y_domain%beg)/real(n_glb + 1, wp)

                    if (proc_coords(2) < rem_cells) then
                        y_domain%beg = y_domain%beg + dy*real((n + 1)* &
                                                              proc_coords(2))
                        y_domain%end = y_domain%end - dy*real((n + 1)* &
                                                              (num_procs_y - proc_coords(2) - 1) &
                                                              - (num_procs_y - rem_cells))
                    else
                        y_domain%beg = y_domain%beg + dy*real((n + 1)* &
                                                              proc_coords(2) + rem_cells)
                        y_domain%end = y_domain%end - dy*real((n + 1)* &
                                                              (num_procs_y - proc_coords(2) - 1))
                    end if
                end if
#endif
            end if

            ! 1D Cartesian Processor Topology
        else

            ! Optimal processor topology
            num_procs_x = num_procs

            ! Creating new communicator using the Cartesian topology
            call MPI_CART_CREATE(MPI_COMM_WORLD, 1, (/num_procs_x/), &
                                 (/.true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding the Cartesian coordinates of the local process
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 1, &
                                 proc_coords, ierr)

        end if

        ! Global Parameters for x-direction

        ! Number of remaining cells
        rem_cells = mod(m + 1, num_procs_x)

        ! Optimal number of cells per processor
        m = (m + 1)/num_procs_x - 1

        ! Distributing the remaining cells
        do i = 1, rem_cells
            if (proc_coords(1) == i - 1) then
                m = m + 1; exit
            end if
        end do

        call s_update_cell_bounds(cells_bounds, m, n, p)

        ! Boundary condition at the beginning
        if (proc_coords(1) > 0 .or. (bc_x%beg == BC_PERIODIC .and. num_procs_x > 1)) then
            proc_coords(1) = proc_coords(1) - 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%beg, ierr)
            proc_coords(1) = proc_coords(1) + 1
        end if

        ! Boundary condition at the end
        if (proc_coords(1) < num_procs_x - 1 .or. (bc_x%end == BC_PERIODIC .and. num_procs_x > 1)) then
            proc_coords(1) = proc_coords(1) + 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%end, ierr)
            proc_coords(1) = proc_coords(1) - 1
        end if

#ifdef MFC_POST_PROCESS
        ! Ghost zone at the beginning
        if (proc_coords(1) > 0 .and. format == 1 .and. n > 0) then
            offset_x%beg = 2
        else
            offset_x%beg = 0
        end if

        ! Ghost zone at the end
        if (proc_coords(1) < num_procs_x - 1 .and. format == 1 .and. n > 0) then
            offset_x%end = 2
        else
            offset_x%end = 0
        end if
#endif

        ! Beginning and end sub-domain boundary locations
        if (parallel_io) then
            if (proc_coords(1) < rem_cells) then
                start_idx(1) = (m + 1)*proc_coords(1)
            else
                start_idx(1) = (m + 1)*proc_coords(1) + rem_cells
            end if
        else
#ifdef MFC_PRE_PROCESS
            if (old_grid .neqv. .true.) then
                dx = (x_domain%end - x_domain%beg)/real(m_glb + 1, wp)

                if (proc_coords(1) < rem_cells) then
                    x_domain%beg = x_domain%beg + dx*real((m + 1)* &
                                                          proc_coords(1))
                    x_domain%end = x_domain%end - dx*real((m + 1)* &
                                                          (num_procs_x - proc_coords(1) - 1) &
                                                          - (num_procs_x - rem_cells))
                else
                    x_domain%beg = x_domain%beg + dx*real((m + 1)* &
                                                          proc_coords(1) + rem_cells)
                    x_domain%end = x_domain%end - dx*real((m + 1)* &
                                                          (num_procs_x - proc_coords(1) - 1))
                end if
            end if
#endif
        end if
#endif

    end subroutine s_mpi_decompose_computational_domain

    !>  The goal of this procedure is to populate the buffers of
        !!      the grid variables by communicating with the neighboring
        !!      processors. Note that only the buffers of the cell-width
        !!      distributions are handled in such a way. This is because
        !!      the buffers of cell-boundary locations may be calculated
        !!      directly from those of the cell-width distributions.
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
#ifndef MFC_PRE_PROCESS
    subroutine s_mpi_sendrecv_grid_variables_buffers(mpi_dir, pbc_loc)

        integer, intent(in) :: mpi_dir
        integer, intent(in) :: pbc_loc

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! MPI Communication in x-direction
        if (mpi_dir == 1) then

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_x%end >= 0) then      ! PBC at the beginning and end

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        dx(m - buff_size + 1), buff_size, &
                        mpi_p, bc_x%end, 0, &
                        dx(-buff_size), buff_size, &
                        mpi_p, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the beginning only

                    ! Send/receive buffer to/from bc_x%beg/bc_x%beg
                    call MPI_SENDRECV( &
                        dx(0), buff_size, &
                        mpi_p, bc_x%beg, 1, &
                        dx(-buff_size), buff_size, &
                        mpi_p, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            else                        ! PBC at the end

                if (bc_x%beg >= 0) then      ! PBC at the end and beginning

                    ! Send/receive buffer to/from bc_x%beg/bc_x%end
                    call MPI_SENDRECV( &
                        dx(0), buff_size, &
                        mpi_p, bc_x%beg, 1, &
                        dx(m + 1), buff_size, &
                        mpi_p, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the end only

                    ! Send/receive buffer to/from bc_x%end/bc_x%end
                    call MPI_SENDRECV( &
                        dx(m - buff_size + 1), buff_size, &
                        mpi_p, bc_x%end, 0, &
                        dx(m + 1), buff_size, &
                        mpi_p, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            end if
            ! END: MPI Communication in x-direction

            ! MPI Communication in y-direction
        elseif (mpi_dir == 2) then

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_y%end >= 0) then      ! PBC at the beginning and end

                    ! Send/receive buffer to/from bc_y%end/bc_y%beg
                    call MPI_SENDRECV( &
                        dy(n - buff_size + 1), buff_size, &
                        mpi_p, bc_y%end, 0, &
                        dy(-buff_size), buff_size, &
                        mpi_p, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the beginning only

                    ! Send/receive buffer to/from bc_y%beg/bc_y%beg
                    call MPI_SENDRECV( &
                        dy(0), buff_size, &
                        mpi_p, bc_y%beg, 1, &
                        dy(-buff_size), buff_size, &
                        mpi_p, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            else                        ! PBC at the end

                if (bc_y%beg >= 0) then      ! PBC at the end and beginning

                    ! Send/receive buffer to/from bc_y%beg/bc_y%end
                    call MPI_SENDRECV( &
                        dy(0), buff_size, &
                        mpi_p, bc_y%beg, 1, &
                        dy(n + 1), buff_size, &
                        mpi_p, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the end only

                    ! Send/receive buffer to/from bc_y%end/bc_y%end
                    call MPI_SENDRECV( &
                        dy(n - buff_size + 1), buff_size, &
                        mpi_p, bc_y%end, 0, &
                        dy(n + 1), buff_size, &
                        mpi_p, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            end if
            ! END: MPI Communication in y-direction

            ! MPI Communication in z-direction
        else

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_z%end >= 0) then      ! PBC at the beginning and end

                    ! Send/receive buffer to/from bc_z%end/bc_z%beg
                    call MPI_SENDRECV( &
                        dz(p - buff_size + 1), buff_size, &
                        mpi_p, bc_z%end, 0, &
                        dz(-buff_size), buff_size, &
                        mpi_p, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the beginning only

                    ! Send/receive buffer to/from bc_z%beg/bc_z%beg
                    call MPI_SENDRECV( &
                        dz(0), buff_size, &
                        mpi_p, bc_z%beg, 1, &
                        dz(-buff_size), buff_size, &
                        mpi_p, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            else                        ! PBC at the end

                if (bc_z%beg >= 0) then      ! PBC at the end and beginning

                    ! Send/receive buffer to/from bc_z%beg/bc_z%end
                    call MPI_SENDRECV( &
                        dz(0), buff_size, &
                        mpi_p, bc_z%beg, 1, &
                        dz(p + 1), buff_size, &
                        mpi_p, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the end only

                    ! Send/receive buffer to/from bc_z%end/bc_z%end
                    call MPI_SENDRECV( &
                        dz(p - buff_size + 1), buff_size, &
                        mpi_p, bc_z%end, 0, &
                        dz(p + 1), buff_size, &
                        mpi_p, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            end if

        end if
        ! END: MPI Communication in z-direction
#endif

    end subroutine s_mpi_sendrecv_grid_variables_buffers
#endif

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_mpi_common_module

#ifdef MFC_MPI
        deallocate (buff_send, buff_recv)
#endif

    end subroutine s_finalize_mpi_common_module

end module m_mpi_common
