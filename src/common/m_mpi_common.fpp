
#:include 'macros.fpp'

!> @brief The module serves as a proxy to the parameters and subroutines
!!          available in the MPI implementation's MPI module. Specifically,
!!          the purpose of the proxy is to harness basic MPI commands into
!!          more complicated procedures as to accomplish the communication
!!          goals for the simulation.
module m_mpi_common

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_nvtx

    use m_helper

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    ! ==========================================================================

    implicit none

    !> @name Generic flags used to identify and report MPI errors
    !> @{
    integer, private :: ierr
    !> @}

#ifndef MFC_PRE_PROCESS
    real(kind(0d0)), private, allocatable, dimension(:), target :: q_cons_buff_send !<
    !! This variable is utilized to pack and send the buffer of the cell-average
    !! conservative variables, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    real(kind(0d0)), private, allocatable, dimension(:), target :: q_cons_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the cell-
    !! average conservative variables, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    !$acc declare create(q_cons_buff_send, q_cons_buff_recv)

    integer :: v_size

    !$acc declare create(v_size)
#endif

    integer, dimension(3, -1:1) :: neighbor_procs
    !$acc declare create(neighbor_procs)

contains

    !> The subroutine initializes the MPI execution environment
        !!      and queries both the number of processors which will be
        !!      available for the job and the local processor rank.
    subroutine s_mpi_initialize

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
            print '(A)', 'Unable to initialize MPI environment. Exiting ...'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Querying the number of processors available for the job
        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

        ! Querying the rank of the local processor
        call MPI_COMM_RANK(MPI_COMM_WORLD, proc_rank, ierr)

#endif

    end subroutine s_mpi_initialize

    subroutine s_initialize_mpi_common_module()

#if defined(MFC_MPI) && !defined(MFC_PRE_PROCESS)

        ! Allocating q_cons_buff_send/recv and ib_buff_send/recv. Please note that
        ! for the sake of simplicity, both variables are provided sufficient
        ! storage to hold the largest buffer in the computational domain.

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE(q_cons_buff_send(0:-1 + buff_size*(sys_size + 2*nb*4)* &
                                             & (m + 2*buff_size + 1)* &
                                             & (n + 2*buff_size + 1)* &
                                             & (p + 2*buff_size + 1)/ &
                                             & (min(m, n, p) + 2*buff_size + 1)))
                else
                    @:ALLOCATE(q_cons_buff_send(0:-1 + buff_size*(sys_size + 2*nb*4)* &
                                             & (max(m, n) + 2*buff_size + 1)))
                end if
            else
                @:ALLOCATE(q_cons_buff_send(0:-1 + buff_size*(sys_size + 2*nb*4)))
            end if

            @:ALLOCATE(q_cons_buff_recv(0:ubound(q_cons_buff_send, 1)))

            v_size = sys_size + 2*nb*4
        else
#endif
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE(q_cons_buff_send(0:-1 + buff_size*sys_size* &
                                             & (m + 2*buff_size + 1)* &
                                             & (n + 2*buff_size + 1)* &
                                             & (p + 2*buff_size + 1)/ &
                                             & (min(m, n, p) + 2*buff_size + 1)))
                else
                    @:ALLOCATE(q_cons_buff_send(0:-1 + buff_size*sys_size* &
                                             & (max(m, n) + 2*buff_size + 1)))
                end if
            else
                @:ALLOCATE(q_cons_buff_send(0:-1 + buff_size*sys_size))
            end if

            @:ALLOCATE(q_cons_buff_recv(0:ubound(q_cons_buff_send, 1)))

            v_size = sys_size
#if MFC_SIMULATION
        end if
#endif

!$acc update device(v_size)

#endif

    end subroutine s_initialize_mpi_common_module

    subroutine s_finalize_mpi_common_module()

#if defined(MFC_MPI) && !defined(MFC_PRE_PROCESS)
        @:DEALLOCATE(q_cons_buff_send, q_cons_buff_recv)
#endif

    end subroutine s_finalize_mpi_common_module

    !! @param q_cons_vf Conservative variables
    !! @param ib_markers track if a cell is within the immersed boundary
    !! @param levelset closest distance from every cell to the IB
    !! @param levelset_norm normalized vector from every cell to the closest point to the IB
    subroutine s_initialize_mpi_data(q_cons_vf, ib_markers, levelset, levelset_norm)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_cons_vf

        type(integer_field), &
            optional, &
            intent(in) :: ib_markers

        type(levelset_field), &
            optional, &
            intent(IN) :: levelset

        type(levelset_norm_field), &
            optional, &
            intent(IN) :: levelset_norm

        integer, dimension(num_dims) :: sizes_glb, sizes_loc
        integer, dimension(1) :: airfoil_glb, airfoil_loc, airfoil_start

#ifdef MFC_MPI

        ! Generic loop iterator
        integer :: i, j, q, k, l

        do i = 1, sys_size
            MPI_IO_DATA%var(i)%sf => q_cons_vf(i)%sf(0:m, 0:n, 0:p)
        end do

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
        do i = 1, sys_size
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                          MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), ierr)
            call MPI_TYPE_COMMIT(MPI_IO_DATA%view(i), ierr)
        end do

#ifndef MFC_POST_PROCESS
        if (qbmm .and. .not. polytropic) then
            do i = sys_size + 1, sys_size + 2*nb*4
                call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                              MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), ierr)
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
                                          MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_levelset_DATA%view, ierr)
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                          MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_levelsetnorm_DATA%view, ierr)

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
                                                  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_airfoil_IB_DATA%view(1), ierr)
                    call MPI_TYPE_COMMIT(MPI_IO_airfoil_IB_DATA%view(1), ierr)

#ifdef MFC_PRE_PROCESS
                    do i = 1, Np
                        MPI_IO_airfoil_IB_DATA%var(Np + i)%x = airfoil_grid_u(i)%x
                        MPI_IO_airfoil_IB_DATA%var(Np + i)%y = airfoil_grid_u(i)%y
                    end do
#endif
                    call MPI_TYPE_CREATE_SUBARRAY(1, airfoil_glb, airfoil_loc, airfoil_start, &
                                                  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_airfoil_IB_DATA%view(2), ierr)
                    call MPI_TYPE_COMMIT(MPI_IO_airfoil_IB_DATA%view(2), ierr)

                end if
            end do

        end if
#endif

#endif

    end subroutine s_initialize_mpi_data

    subroutine mpi_bcast_time_step_values(proc_time, time_avg)

        real(kind(0d0)), dimension(0:num_procs - 1), intent(inout) :: proc_time
        real(kind(0d0)), intent(inout) :: time_avg

#ifdef MFC_MPI

        call MPI_GATHER(time_avg, 1, MPI_DOUBLE_PRECISION, proc_time(0), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

#endif

    end subroutine mpi_bcast_time_step_values

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
    subroutine s_mpi_reduce_stability_criteria_extrema(icfl_max_loc, &
                                                       vcfl_max_loc, &
                                                       ccfl_max_loc, &
                                                       Rc_min_loc, &
                                                       icfl_max_glb, &
                                                       vcfl_max_glb, &
                                                       ccfl_max_glb, &
                                                       Rc_min_glb)

        real(kind(0d0)), intent(in) :: icfl_max_loc
        real(kind(0d0)), intent(in) :: vcfl_max_loc
        real(kind(0d0)), intent(in) :: ccfl_max_loc
        real(kind(0d0)), intent(in) :: Rc_min_loc

        real(kind(0d0)), intent(out) :: icfl_max_glb
        real(kind(0d0)), intent(out) :: vcfl_max_glb
        real(kind(0d0)), intent(out) :: ccfl_max_glb
        real(kind(0d0)), intent(out) :: Rc_min_glb

#ifdef MFC_SIMULATION
#ifdef MFC_MPI

        ! Reducing local extrema of ICFL, VCFL, CCFL and Rc numbers to their
        ! global extrema and bookkeeping the results on the rank 0 processor
        call MPI_REDUCE(icfl_max_loc, icfl_max_glb, 1, &
                        MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                        MPI_COMM_WORLD, ierr)

        if (viscous) then
            call MPI_REDUCE(vcfl_max_loc, vcfl_max_glb, 1, &
                            MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                            MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(Rc_min_loc, Rc_min_glb, 1, &
                            MPI_DOUBLE_PRECISION, MPI_MIN, 0, &
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
    subroutine s_mpi_allreduce_sum(var_loc, var_glb)

        real(kind(0d0)), intent(in) :: var_loc
        real(kind(0d0)), intent(out) :: var_glb

#ifdef MFC_MPI

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
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
    subroutine s_mpi_allreduce_min(var_loc, var_glb)

        real(kind(0d0)), intent(in) :: var_loc
        real(kind(0d0)), intent(out) :: var_glb

#ifdef MFC_MPI

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
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
    subroutine s_mpi_allreduce_max(var_loc, var_glb)

        real(kind(0d0)), intent(in) :: var_loc
        real(kind(0d0)), intent(out) :: var_glb

#ifdef MFC_MPI

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                           MPI_MAX, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_max

    !>  The following subroutine takes the inputted variable and
        !!      determines its minimum value on the entire computational
        !!      domain. The result is stored back into inputted variable.
        !!  @param var_loc holds the local value to be reduced among
        !!      all the processors in communicator. On output, the variable holds
        !!      the minimum value, reduced amongst all of the local values.
    subroutine s_mpi_reduce_min(var_loc)

        real(kind(0d0)), intent(inout) :: var_loc

#ifdef MFC_MPI

        ! Temporary storage variable that holds the reduced minimum value
        real(kind(0d0)) :: var_glb

        ! Performing reduction procedure and eventually storing its result
        ! into the variable that was initially inputted into the subroutine
        call MPI_REDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                        MPI_MIN, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(var_glb, 1, MPI_DOUBLE_PRECISION, &
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
    subroutine s_mpi_reduce_maxloc(var_loc)

        real(kind(0d0)), dimension(2), intent(inout) :: var_loc

#ifdef MFC_MPI

        real(kind(0d0)), dimension(2) :: var_glb  !<
            !! Temporary storage variable that holds the reduced maximum value
            !! and the rank of the processor with which the value is associated

        ! Performing reduction procedure and eventually storing its result
        ! into the variable that was initially inputted into the subroutine
        call MPI_REDUCE(var_loc, var_glb, 1, MPI_2DOUBLE_PRECISION, &
                        MPI_MAXLOC, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(var_glb, 1, MPI_2DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, ierr)

        var_loc = var_glb

#endif

    end subroutine s_mpi_reduce_maxloc

    !> The subroutine terminates the MPI execution environment.
        !! @param prnt error message to be printed
    subroutine s_mpi_abort(prnt)

        character(len=*), intent(in), optional :: prnt

        if (present(prnt)) then
            print *, prnt
            call flush (6)

        end if

#ifndef MFC_MPI

        stop 1

#else

        ! Terminating the MPI environment
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

#endif

    end subroutine s_mpi_abort

    !>Halts all processes until all have reached barrier.
    subroutine s_mpi_barrier

#ifdef MFC_MPI

        ! Calling MPI_BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_barrier

    !> The subroutine finalizes the MPI execution environment.
    subroutine s_mpi_finalize

#ifdef MFC_MPI

        ! Finalizing the MPI environment
        call MPI_FINALIZE(ierr)

#endif

    end subroutine s_mpi_finalize

#ifndef MFC_PRE_PROCESS
    !>  The goal of this procedure is to populate the buffers of
        !!      the cell-average conservative variables by communicating
        !!      with the neighboring processors.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
    subroutine s_mpi_sendrecv_variables_buffers(q_cons_vf, &
                                                bc_id_sfs, &
#ifdef MFC_SIMULATION
                                                pb, mv, &
#endif
                                                mpi_dir, &
                                                pbc_loc, &
                                                bc_id_has_bc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(t_bc_id_sf), dimension(1:3, -1:1), intent(in) :: bc_id_sfs

#ifdef MFC_SIMULATION
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
#endif
        integer, intent(in) :: mpi_dir, pbc_loc

        logical, dimension(1:3, -1:1, -num_bcs_max:0), intent(in) :: bc_id_has_bc

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        type(bc_patch_parameters) :: bc

        integer :: buffer_counts(1:3)

        real(kind(0d0)), pointer :: p_send, p_recv

        integer :: r

#ifdef MFC_MPI

!$acc update device(v_size)

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            buffer_counts = (/ &
                            buff_size*(sys_size + 2*nb*4)*(n + 1)*(p + 1), &
                            buff_size*(sys_size + 2*nb*4)*(m + 2*buff_size + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                            /)
        else
#endif
            buffer_counts = (/ &
                            buff_size*sys_size*(n + 1)*(p + 1), &
                            buff_size*sys_size*(m + 2*buff_size + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                            /)
#ifdef MFC_SIMULATION
        end if
#endif

        if (bc_id_has_bc(mpi_dir, -pbc_loc, 0)) then
            call nvtxStartRange("RHS-COMM-PACKBUF")

            #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="mpi_dir", loc="-pbc_loc", inner_loops=[("i", 1, "sys_size")])
                q_cons_buff_send(pack_idr) = q_cons_vf(i)%sf(sx, sy, sz)
            #:endblock

#ifdef MFC_SIMULATION
            if (qbmm .and. .not. polytropic) then
                #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="mpi_dir", loc="-pbc_loc", inner_loops=[("i", "sys_size + 1", "sys_size + 4"), ("q", 1, "nb")])
                    q_cons_buff_send(pack_idr + (q - 1)*4) = pb(sx, sy, sz, i - sys_size, q)
                    q_cons_buff_send(pack_idr + (q - 1)*4 + nb*4) = mv(sx, sy, sz, i - sys_size, q)
                #:endblock
            end if
#endif

            call nvtxEndRange ! Packbuf
        end if

        p_send => q_cons_buff_send(0)
        p_recv => q_cons_buff_recv(0)

        #:for rdma_mpi in [False, True]
            if (rdma_mpi .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
                #:if rdma_mpi
                    !$acc data attach(p_send, p_recv)
                    !$acc host_data use_device(p_send, p_recv)

                    call nvtxStartRange("RHS-COMM-SENDRECV-RDMA")
                #:else
                    call nvtxStartRange("RHS-COMM-DEV2HOST")
                    !$acc update host(q_cons_buff_send)
                    call nvtxEndRange

                    call nvtxStartRange("RHS-COMM-SENDRECV-NO-RMDA")
                #:endif

                call MPI_SENDRECV( &
                    p_send, buffer_counts(mpi_dir), MPI_DOUBLE_PRECISION, neighbor_procs(mpi_dir, -pbc_loc), 0, &
                    p_recv, buffer_counts(mpi_dir), MPI_DOUBLE_PRECISION, neighbor_procs(mpi_dir, +pbc_loc), 0, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                call nvtxEndRange ! RHS-MPI-SENDRECV-(NO)-RDMA

                #:if rdma_mpi
                    !$acc end host_data
                    !$acc end data
                    !$acc wait
                #:else
                    call nvtxStartRange("RHS-COMM-HOST2DEV")
                    !$acc update device(q_cons_buff_recv)
                    call nvtxEndRange
                #:endif
            end if
        #:endfor

        if (bc_id_has_bc(mpi_dir, pbc_loc, 0)) then
            call nvtxStartRange("RHS-COMM-UNPACKBUF")

            #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="mpi_dir", loc="pbc_loc", inner_loops=[("i", 1, "sys_size")])
                if (bc_id_sfs(mpi_dir, pbc_loc)%sf(exlhs, eylhs, ezlhs)%type >= 0) then
                    q_cons_vf(i)%sf(x, y, z) = q_cons_buff_recv(pack_idr)
                end if
            #:endblock

#ifdef MFC_SIMULATION
            if (qbmm .and. .not. polytropic) then
                #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="mpi_dir", loc="pbc_loc", inner_loops=[("i", "sys_size + 1", "sys_size + 4"), ("q", 1, "nb")])
                    if (bc_id_sfs(mpi_dir, pbc_loc)%sf(exlhs, eylhs, ezlhs)%type >= 0) then
                        pb(x, y, z, i - sys_size, q) = q_cons_buff_recv(pack_idr + (q - 1)*4)
                        mv(x, y, z, i - sys_size, q) = q_cons_buff_recv(pack_idr + (q - 1)*4 + nb*4)
                    end if
                #:endblock
            end if
#endif

            call nvtxEndRange
        end if
#endif

    end subroutine s_mpi_sendrecv_variables_buffers

    !>  The goal of this procedure is to populate the buffers of
        !!      the grid variables by communicating with the neighboring
        !!      processors. Note that only the buffers of the cell-width
        !!      distributions are handled in such a way. This is because
        !!      the buffers of cell-boundary locations may be calculated
        !!      directly from those of the cell-width distributions.
    subroutine s_mpi_sendrecv_grid_spacing_buffers(bc_id_sfs)

        type(t_bc_id_sf), dimension(1:3, -1:1), intent(in) :: bc_id_sfs

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        type(t_bc_id) :: bc

        integer :: iter_loc

        integer :: send_offset, recv_offset
        integer :: extra_cell_count

        #:for cmp, dir_id, extent in zip(['x', 'y', 'z'], [1, 2, 3], ['m', 'n', 'p'])

            if (${dir_id}$ > num_dims) then
                return
            end if

            do iter_loc = -1, 1, 2

#ifdef MFC_MPI
                if (iter_loc == -1) then
                    send_offset = ${extent}$-buff_size + 1
                    recv_offset = -buff_size
                else
                    send_offset = 0
                    recv_offset = ${extent}$+1
                end if

                call MPI_SENDRECV( &
                    d${cmp}$ (send_offset), buff_size, &
                    MPI_DOUBLE_PRECISION, neighbor_procs(${dir_id}$, -iter_loc), 0, &
                    d${cmp}$ (recv_offset), buff_size, &
                    MPI_DOUBLE_PRECISION, neighbor_procs(${dir_id}$, +iter_loc), 0, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif

                ! Note: This does a little TOO much work (iterating over too many cells)
                !       but this is because dx, dy, dz SHOULD be two dimensional arrays,
                !       not one-dimensional arrays. Indeed, the boundary conditions
                !       can vary along two axes, not just one if the case is 3D.
                #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir=dir_id, loc="iter_loc", gpu=False)
                    bc = bc_id_sfs(${dir_id}$, iter_loc)%sf(exlhs, eylhs, ezlhs)

                    if (bc%type <= -3) then
                        d${cmp}$ (${cmp}$) = d${cmp}$ (e${cmp}$)
                    elseif (bc%type == -2) then
                        d${cmp}$ (${cmp}$) = d${cmp}$ (s${cmp}$)
                    elseif (bc%type == -1) then
                        d${cmp}$ (${cmp}$) = d${cmp}$ (p${cmp}$)
                    end if
                #:endblock

            end do

#ifdef MFC_POST_PROCESS
            extra_cell_count = offset_${cmp}$%beg
#else
            extra_cell_count = buff_size
#endif

            do i = 1, extra_cell_count
                ! Computing the cell-boundary locations buffer, at the beginning of
                ! the coordinate direction, from the cell-width distribution buffer.
                ${cmp}$_cb(-1 - i) = ${cmp}$_cb(-i) - d${cmp}$ (-i)
            end do

#ifdef MFC_POST_PROCESS
            extra_cell_count = offset_${cmp}$%end
#else
            extra_cell_count = buff_size
#endif

            do i = 1, extra_cell_count
                ! Populating the cell-boundary locations buffer, at the end of the
                ! coordinate direction, from buffer of the cell-width distribution.
                ${cmp}$_cb(${extent}$+i) = ${cmp}$_cb(${extent}$+(i - 1)) + d${cmp}$ (${extent}$+i)
            end do

            do i = 1, buff_size
                ! Computing the cell-center locations buffer, at the beginning of
                ! the coordinate direction, from the cell-width distribution buffer.
                ${cmp}$_cc(-i) = ${cmp}$_cc(1 - i) - (d${cmp}$ (1 - i) + d${cmp}$ (-i))/2d0

                ! Populating the cell-center locations buffer, at the end of the
                ! coordinate direction, from buffer of the cell-width distribution.
                ${cmp}$_cc(${extent}$+i) = ${cmp}$_cc(${extent}$+(i - 1)) + (d${cmp}$ (${extent}$+(i - 1)) + d${cmp}$ (${extent}$+i))/2d0
            end do
        #:endfor

    end subroutine s_mpi_sendrecv_grid_spacing_buffers
#endif

    !>  The purpose of this procedure is to optimally decompose
        !!      the computational domain among the available processors.
        !!      This is performed by attempting to award each processor,
        !!      in each of the coordinate directions, approximately the
        !!      same number of cells, and then recomputing the affected
        !!      global parameters.
    subroutine s_mpi_decompose_computational_domain

#ifdef MFC_MPI

        real(kind(0d0)) :: tmp_num_procs_x, tmp_num_procs_y, tmp_num_procs_z !<
                !! Non-optimal number of processors in the x-, y- and z-directions

        real(kind(0d0)) :: fct_min !<
                !! Processor factorization (fct) minimization parameter

        integer :: MPI_COMM_CART !<
                !! Cartesian processor topology communicator

        integer, dimension(1:3) :: rem_cells !<
                !! Remaining number of cells, in a particular coordinate direction,
                !! after the majority is divided up among the available processors

        integer :: i, j !< Generic loop iterators

        integer :: iter_dir

#endif

        proc_coords(:) = 0
        proc_nums(:) = 1

#ifdef MFC_MPI

        if (num_procs == 1 .and. parallel_io) then
            do i = 1, num_dims
                start_idx(i) = 0
            end do
            return
        end if

        ! 3D Cartesian Processor Topology ==================================
        if (n > 0) then

            if (p > 0) then

                if (cyl_coord .and. p > 0) then
                    ! Implement pencil processor blocking if using cylindrical coordinates so
                    ! that all cells in azimuthal direction are stored on a single processor.
                    ! This is necessary for efficient application of Fourier filter near axis.

                    ! Initial values of the processor factorization optimization
                    proc_nums(1) = 1
                    proc_nums(2) = num_procs
                    proc_nums(3) = 1
                    ierr = -1

                    ! Computing minimization variable for these initial values
                    tmp_num_procs_x = proc_nums(1)
                    tmp_num_procs_y = proc_nums(2)
                    tmp_num_procs_z = proc_nums(3)
                    fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                                       - (n + 1)/tmp_num_procs_y)

                    ! Searching for optimal computational domain distribution
                    do i = 1, num_procs

                        if (mod(num_procs, i) == 0 &
                            .and. &
                            (m + 1)/i >= num_stcls_min*weno_order) then

                            tmp_num_procs_x = i
                            tmp_num_procs_y = num_procs/i

                            if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                               - (n + 1)/tmp_num_procs_y) &
                                .and. &
                                (n + 1)/tmp_num_procs_y &
                                >= &
                                num_stcls_min*weno_order) then

                                proc_nums(1) = i
                                proc_nums(2) = num_procs/i
                                fct_min = abs((m + 1)/tmp_num_procs_x &
                                              - (n + 1)/tmp_num_procs_y)
                                ierr = 0

                            end if

                        end if

                    end do

                else

                    ! Initial estimate of optimal processor topology
                    proc_nums(1) = 1
                    proc_nums(2) = 1
                    proc_nums(3) = num_procs
                    ierr = -1

                    ! Benchmarking the quality of this initial guess
                    tmp_num_procs_x = proc_nums(1)
                    tmp_num_procs_y = proc_nums(2)
                    tmp_num_procs_z = proc_nums(3)
                    fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                                       - (n + 1)/tmp_num_procs_y) &
                              + 10d0*abs((n + 1)/tmp_num_procs_y &
                                         - (p + 1)/tmp_num_procs_z)

                    ! Optimization of the initial processor topology
                    do i = 1, num_procs

                        if (mod(num_procs, i) == 0 &
                            .and. &
                            (m + 1)/i >= num_stcls_min*weno_order) then

                            do j = 1, num_procs/i

                                if (mod(num_procs/i, j) == 0 &
                                    .and. &
                                    (n + 1)/j >= num_stcls_min*weno_order) then

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
                                        num_stcls_min*weno_order) &
                                        then

                                        proc_nums(1) = i
                                        proc_nums(2) = j
                                        proc_nums(3) = num_procs/(i*j)
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
                                     'weno_order. Exiting ...')
                end if

                ! Creating new communicator using the Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 3, (/proc_nums(1), &
                                                          proc_nums(2), proc_nums(3)/), &
                                     (/.true., .true., .true./), &
                                     .false., MPI_COMM_CART, ierr)

                ! Finding the Cartesian coordinates of the local process
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 3, &
                                     proc_coords, ierr)
                ! END: 3D Cartesian Processor Topology =============================

                ! Global Parameters for z-direction ================================

                ! Number of remaining cells
                rem_cells(3) = mod(p + 1, proc_nums(3))

#ifdef MFC_PRE_PROCESS
                ! Preliminary uniform cell-width spacing
                if (old_grid .neqv. .true.) then
                    dz = (z_domain%end - z_domain%beg)/real(p + 1, kind(0d0))
                end if
#endif

                ! Optimal number of cells per processor
                p = (p + 1)/proc_nums(3) - 1

                ! Distributing the remaining cells
                do i = 1, rem_cells(3)
                    if (proc_coords(3) == i - 1) then
                        p = p + 1; exit
                    end if
                end do

                ! ==================================================================

                ! 2D Cartesian Processor Topology ==================================
            else

                ! Initial estimate of optimal processor topology
                proc_nums(1) = 1
                proc_nums(2) = num_procs
                ierr = -1

                ! Benchmarking the quality of this initial guess
                tmp_num_procs_x = proc_nums(1)
                tmp_num_procs_y = proc_nums(2)
                fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                                   - (n + 1)/tmp_num_procs_y)

                ! Optimization of the initial processor topology
                do i = 1, num_procs

                    if (mod(num_procs, i) == 0 &
                        .and. &
                        (m + 1)/i >= num_stcls_min*weno_order) then

                        tmp_num_procs_x = i
                        tmp_num_procs_y = num_procs/i

                        if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                           - (n + 1)/tmp_num_procs_y) &
                            .and. &
                            (n + 1)/tmp_num_procs_y &
                            >= &
                            num_stcls_min*weno_order) then

                            proc_nums(1) = i
                            proc_nums(2) = num_procs/i
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
                                     'weno_order. Exiting ...')
                end if

                ! Creating new communicator using the Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 2, (/proc_nums(1), &
                                                          proc_nums(2)/), (/.true., &
                                                                            .true./), .false., MPI_COMM_CART, &
                                     ierr)

                ! Finding the Cartesian coordinates of the local process
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, &
                                     proc_coords, ierr)

            end if
            ! END: 2D Cartesian Processor Topology =============================

            ! Global Parameters for y-direction ================================

            ! Number of remaining cells
            rem_cells(2) = mod(n + 1, proc_nums(2))

#ifdef MFC_PRE_PROCESS
            ! Preliminary uniform cell-width spacing
            if (old_grid .neqv. .true.) then
                dy = (y_domain%end - y_domain%beg)/real(n + 1, kind(0d0))
            end if
#endif

            ! Optimal number of cells per processor
            n = (n + 1)/proc_nums(2) - 1

            ! Distributing the remaining cells
            do i = 1, rem_cells(2)
                if (proc_coords(2) == i - 1) then
                    n = n + 1; exit
                end if
            end do

            ! ==================================================================

            ! 1D Cartesian Processor Topology ==================================
        else

            ! Optimal processor topology
            proc_nums(1) = num_procs

#ifdef MFC_POST_PROCESS
            ! Number of cells in undecomposed computational domain needed
            ! for sub-domain reassembly during formatted data output
            m_root = m
#endif

            ! Creating new communicator using the Cartesian topology
            call MPI_CART_CREATE(MPI_COMM_WORLD, 1, (/proc_nums(1)/), &
                                 (/.true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding the Cartesian coordinates of the local process
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 1, &
                                 proc_coords, ierr)

        end if
        ! ==================================================================

        ! Global Parameters for x-direction ================================

        ! Number of remaining cells
        rem_cells(1) = mod(m + 1, proc_nums(1))

#ifdef MFC_PRE_PROCESS
        ! Preliminary uniform cell-width spacing
        if (old_grid .neqv. .true.) then
            dx = (x_domain%end - x_domain%beg)/real(m + 1, kind(0d0))
        end if
#endif

        ! Optimal number of cells per processor
        m = (m + 1)/proc_nums(1) - 1

        ! Distributing the remaining cells
        do i = 1, rem_cells(1)
            if (proc_coords(1) == i - 1) then
                m = m + 1; exit
            end if
        end do
        ! ==================================================================

        #:for cmp, dir, extent in zip(['z', 'y', 'x'], [3, 2, 1], ['p', 'n', 'm'])

            if (${dir}$ <= num_dims) then

#ifdef MFC_PRE_PROCESS
                ! Beginning and end sub-domain boundary locations
                if (parallel_io .neqv. .true.) then
                    if (old_grid .neqv. .true.) then
                        if (proc_coords(${dir}$) < rem_cells(${dir}$)) then
                            ${cmp}$_domain%beg = ${cmp}$_domain%beg + d${cmp}$*real((${extent}$+1)* &
                                                                                    proc_coords(${dir}$))
                            ${cmp}$_domain%end = ${cmp}$_domain%end - d${cmp}$*real((${extent}$+1)* &
                                                                                    (proc_nums(${dir}$) - proc_coords(${dir}$) - 1) &
                                                                                    - (proc_nums(${dir}$) - rem_cells(${dir}$)))
                        else
                            ${cmp}$_domain%beg = ${cmp}$_domain%beg + d${cmp}$*real((${extent}$+1)* &
                                                                                    proc_coords(${dir}$) + rem_cells(${dir}$))
                            ${cmp}$_domain%end = ${cmp}$_domain%end - d${cmp}$*real((${extent}$+1)* &
                                                                                    (proc_nums(${dir}$) - proc_coords(${dir}$) - 1))
                        end if
                    end if
                else
                    if (proc_coords(${dir}$) < rem_cells(${dir}$)) then
                        start_idx(${dir}$) = (${extent}$+1)*proc_coords(${dir}$)
                    else
                        start_idx(${dir}$) = (${extent}$+1)*proc_coords(${dir}$) + rem_cells(${dir}$)
                    end if
                end if
#endif

#ifdef MFC_POST_PROCESS
                ! Ghost zone at the beginning
                if (proc_coords(${dir}$) > 0 .and. format == 1) then
                    offset_${cmp}$%beg = 2
                else
                    offset_${cmp}$%beg = 0
                end if

                ! Ghost zone at the end
                if (proc_coords(${dir}$) < proc_nums(${dir}$) - 1 .and. format == 1) then
                    offset_${cmp}$%end = 2
                else
                    offset_${cmp}$%end = 0
                end if
#endif

#ifndef MFC_PRE_PROCESS
                if (parallel_io) then
                    if (proc_coords(${dir}$) < rem_cells(${dir}$)) then
                        start_idx(${dir}$) = (${extent}$+1)*proc_coords(${dir}$)
                    else
                        start_idx(${dir}$) = (${extent}$+1)*proc_coords(${dir}$) + rem_cells(${dir}$)
                    end if
                end if
#endif

            end if

        #:endfor

        do iter_dir = 1, num_dims
            proc_coords(iter_dir) = proc_coords(iter_dir) - 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, neighbor_procs(iter_dir, -1), ierr)
            proc_coords(iter_dir) = proc_coords(iter_dir) + 1
            proc_coords(iter_dir) = proc_coords(iter_dir) + 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, neighbor_procs(iter_dir, +1), ierr)
            proc_coords(iter_dir) = proc_coords(iter_dir) - 1
        end do

        if (proc_nums(1) > 1) then
            if (bc_x%beg == -1 .or. proc_coords(1) > 0) then
                bc_x%beg = neighbor_procs(1, -1)
            end if

            if (bc_x%end == +1 .or. proc_coords(1) < proc_nums(1) - 1) then
                bc_x%end = neighbor_procs(1, +1)
            end if
        end if

        if (num_dims >= 2 .and. proc_nums(2) > 1) then
            if (bc_y%beg == -1 .or. proc_coords(2) > 0) then
                bc_y%beg = neighbor_procs(2, -1)
            end if

            if (bc_y%end == +1 .or. proc_coords(2) < proc_nums(2) - 1) then
                bc_y%end = neighbor_procs(2, +1)
            end if
        end if

        if (num_dims >= 3 .and. proc_nums(3) > 1) then
            if (bc_z%beg == -1 .or. proc_coords(3) > 0) then
                bc_z%beg = neighbor_procs(3, -1)
            end if

            if (bc_z%end == +1 .or. proc_coords(3) < proc_nums(3) - 1) then
                bc_z%end = neighbor_procs(3, +1)
            end if
        end if

#endif

    end subroutine s_mpi_decompose_computational_domain

    subroutine s_prohibit_abort(condition, message)
        character(len=*), intent(in) :: condition, message

        print *, ""
        print *, "===================================================================================================="
        print *, "                                          CASE FILE ERROR                                           "
        print *, "----------------------------------------------------------------------------------------------------"
        print *, "Prohibited condition: ", trim(condition)
        if (len_trim(message) > 0) then
            print *, "Note: ", trim(message)
        end if
        print *, "===================================================================================================="
        print *, ""
        call s_mpi_abort
    end subroutine s_prohibit_abort

end module m_mpi_common
