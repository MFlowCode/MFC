!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers

#:include 'macros.fpp'

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles              !< Bubble dynamics routines

    use m_ibm

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_fftw

    use m_nvtx

    use m_time_tmp             !< Lagrangian solver

    use m_particles

    use m_particles_types

    use m_particles_output

    use m_kernel_functions

    use m_mpi_common
    ! ==========================================================================

    implicit none

    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: rhs_pb

    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: rhs_mv

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

    !$acc declare create(q_cons_ts,q_prim_vf,rhs_vf,q_prim_ts, rhs_mv, rhs_pb)

    TYPE(vector_field), ALLOCATABLE, DIMENSION(:) :: rhs_vp_adapt
    !! Adaptive time step, Lagrangian solver

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_time_steppers_module() ! -----------------------

        type(int_bounds_info) :: ix_t, iy_t, iz_t !<
            !! Indical bounds in the x-, y- and z-directions

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        IF(coupledflag) THEN !Euler-Lagrangian solver
            num_ts = 2
        ELSE
            if (time_stepper == 1) then
                num_ts = 1
            elseif (any(time_stepper == (/2, 3/))) then
                num_ts = 2
            end if
        END IF

        ! Setting the indical bounds in the x-, y- and z-directions
        ix_t%beg = -buff_size; ix_t%end = m + buff_size

        if (n > 0) then
            iy_t%beg = -buff_size; iy_t%end = n + buff_size

            if (p > 0) then
                iz_t%beg = -buff_size; iz_t%end = p + buff_size
            else
                iz_t%beg = 0; iz_t%end = 0
            end if
        else
            iy_t%beg = 0; iy_t%end = 0
            iz_t%beg = 0; iz_t%end = 0
        end if

        ! Allocating the cell-average conservative variables
        @:ALLOCATE(q_cons_ts(1:num_ts))

        do i = 1, num_ts
            @:ALLOCATE(q_cons_ts(i)%vf(1:sys_size))
        end do

        do i = 1, num_ts
            do j = 1, sys_size
                @:ALLOCATE(q_cons_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
            end do
        end do

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then
            @:ALLOCATE(q_prim_ts(0:3))

            do i = 0, 3
                @:ALLOCATE(q_prim_ts(i)%vf(1:sys_size))
            end do

            do i = 0, 3
                do j = 1, sys_size
                    @:ALLOCATE(q_prim_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, &
                        iy_t%beg:iy_t%end, &
                        iz_t%beg:iz_t%end))
                end do
            end do
        end if

        ! Allocating the cell-average primitive variables
        @:ALLOCATE(q_prim_vf(1:sys_size))

        do i = 1, adv_idx%end
            @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end))
        end do

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
            end do
        end if

        @:ALLOCATE(pb_ts(1:2))
        !Initialize bubble variables pb and mv at all quadrature nodes for all R0 bins
        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(pb_ts(1)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ALLOCATE(pb_ts(2)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ALLOCATE(rhs_pb(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
        else if (qbmm .and. polytropic) then
            @:ALLOCATE(pb_ts(1)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ALLOCATE(pb_ts(2)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ALLOCATE(rhs_pb(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
        end if

        @:ALLOCATE(mv_ts(1:2))

        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(mv_ts(1)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ALLOCATE(mv_ts(2)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ALLOCATE(rhs_mv(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
        else if (qbmm .and. polytropic) then
            @:ALLOCATE(mv_ts(1)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ALLOCATE(mv_ts(2)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ALLOCATE(rhs_mv(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
        end if

        if (hypoelasticity) then

            do i = stress_idx%beg, stress_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
            end do
        end if

        ! Allocating the cell-average RHS variables
        @:ALLOCATE(rhs_vf(1:sys_size))

        do i = 1, sys_size
            @:ALLOCATE(rhs_vf(i)%sf(0:m, 0:n, 0:p))
        end do

        ! Allocating the cell-average RHS variable for adaptive method, Lagrangian solver
        IF(coupledflag .OR. (solverapproach.EQ.2)) THEN
            ALLOCATE(rhs_vp_adapt(1:6))
            DO i = 1, 6
                ALLOCATE(rhs_vp_adapt(i)%vf(1:sys_size))
            END DO
            DO j = 1, sys_size
                DO i = 1, 6
                    ALLOCATE(rhs_vp_adapt(i)%vf(j)%sf(0:m,0:n,0:p))
                END DO
            END DO
        END IF

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

    end subroutine s_initialize_time_steppers_module ! ---------------------

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    subroutine s_1st_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 1 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step)

#ifdef DEBUG
        print *, 'got rhs'
#endif

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

#ifdef DEBUG
        print *, 'wrote runtime info'
#endif

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_1st_order_tvd_rk ! ------------------------------------

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_2nd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 2 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if
        ! ==================================================================

        ! Stage 2 of 2 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/2d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/2d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/2d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_2nd_order_tvd_rk ! ------------------------------------

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_3rd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q
        real(kind(0d0)) :: ts_error, denom, error_fraction, time_step_factor !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 3 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! ==================================================================

        ! Stage 2 of 3 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            (3d0*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/4d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    (3d0*pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/4d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    (3d0*mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/4d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! ==================================================================

        ! Stage 3 of 3 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + 2d0*q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + 2d0*dt*rhs_vf(i)%sf(j, k, l))/3d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + 2d0*pb_ts(2)%sf(j, k, l, q, i) &
                                     + 2d0*dt*rhs_pb(j, k, l, q, i))/3d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + 2d0*mv_ts(2)%sf(j, k, l, q, i) &
                                     + 2d0*dt*rhs_mv(j, k, l, q, i))/3d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

        call cpu_time(finish)

        time = time + (finish - start)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_3rd_order_tvd_rk ! ------------------------------------

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step) ! ----------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        do i = 1, sys_size
            !$acc update host(q_prim_vf(i)%sf)
        end do

        if (t_step == t_step_start) then
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 1) then
            do i = 1, sys_size
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 2) then
            do i = 1, sys_size
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 3) then
            do i = 1, sys_size
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        else ! All other timesteps
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_ts(2)%vf(i)%sf(:, :, :)
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_ts(1)%vf(i)%sf(:, :, :)
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_ts(0)%vf(i)%sf(:, :, :)
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        end if

    end subroutine s_time_step_cycling ! -----------------------------------

    !> Cash-Karp Runge-Kutta 4th order time-stepping algorithm
        !! @param t_step Current time-step
    SUBROUTINE rkqs (realtime, hnext, hdid, t_step) ! ----------------------

        USE m_mpi_particles
        USE m_particles_types
        USE m_particles_output

        IMPLICIT NONE

        LOGICAL                       :: largestep
        REAL(KIND(0.D0))              :: newtime,errmax,qtime,hdid,hnext,dttarget
        REAL(KIND(0.D0))              :: RKh,htemp,SAFETY=0.9d0,PGROW=-0.2d0, &
                                         PSHRNK=-0.25d0,ERRCON=1.89d-4
        INTEGER                       :: i,j,k
        REAL(KIND(0.D0)), INTENT(IN)  :: realtime
        INTEGER, INTENT(IN)           :: t_step

        qtime = realtime
        dttarget = dt

        !IF(proc_rank == 0) CALL s_write_run_time_information(q_cons_ts(1)%vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        !> Starting adaptive Runge-Kutta
        501 CONTINUE

        RKh=min(hnext,dttarget)
        RKh=max(Rkh,1.0d-12)
        IF (num_procs > 1) THEN
            CALL get_min(RKh)
        ENDIF

        largestep = .FALSE.
        IF (coupledFlag.OR.bubblesources) THEN
            CALL RKparticledyn(qtime,1,q_cons_ts(1)%vf,t_step,q_prim_vf,rhs_vp_adapt(1)%vf,largestep)
            IF (largestep) STOP 'error at the 0 step'
        ELSE
            CALL RKparticledyn(qtime,1,q_cons_ts(1)%vf,t_step,q_prim_vf)
        ENDIF

        !> Take a step
        502   errmax=0.0d0
        CALL rkck(qtime,RKh,errmax,largestep,t_step)
        hdid=RKh
        hnext=RKh

        !> Update values
        qtime = qtime + hdid
        CALL updateRK (q_cons_ts,.TRUE.,q_prim_vf)

        IF (particlestatFlag) CALL particle_stats ()
        IF (.NOT.stillparticlesflag.AND.(num_procs.GT.1)) CALL transfer_particles
        hnext = min(hnext,dt0)

        RETURN

    END SUBROUTINE rkqs ! ------------------------------------------------------------------------


    !> Cash-Karp Runge-Kutta step
    SUBROUTINE rkck(qtime,RKh,errmax,largestep,t_step) ! ------------------------------------------
        !> USES derivs
        !> Given values for n variables y and their derivatives dydx known at x, use the .fth-order
        !> Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
        !> the incremented variables as yout. Also return an estimate of the local truncation error
        !> in yout using the embedded fourth-order method. The user supplies the subroutine
        !> derivs(x,y,dydx), which returns derivatives dydx at x.

        LOGICAL                                      :: largestep
        REAL(KIND(0.D0))                             :: RKh,qtime,errmax
        INTEGER, INTENT(IN)                          :: t_step
        INTEGER                                      :: i, j
        REAL(KIND(0.D0))                             :: A2=0.2d0,A3=0.3d0,A4=0.6d0,A5=1.0d0,A6=0.875d0
        REAL(KIND(0.D0)),DIMENSION(6) ::                                                    &
        RKcoef1=(/0.2d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/),                    &
        RKcoef2=(/3.0d0/40.0d0,9.0d0/40.0d0,0.0d0,0.0d0,0.0d0,0.0d0/),      &
        RKcoef3=(/0.3d0,-0.9d0,1.2d0,0.0d0,0.0d0,0.0d0/),                   &
        RKcoef4=(/-11.0d0/54.0d0,2.5d0,-70.0d0/27.0d0,35.d0/27.d0,0.0d0,0.0d0/),                                     &
        RKcoef5=(/1631.0d0/55296.0d0,175.0d0/512.0d0,575.d0/13824.d0,44275.d0/110592.d0,253.d0/4096.d0,0.0d0/),      &
        RKcoef6=(/37.d0/378.d0,0.0d0,250.d0/621.d0,125.0d0/594.0d0,0.0d0,512.0d0/1771.0d0/),                         &
        RKcoefE=(/37.d0/378.d0-2825.0d0/27648.0d0,0.0d0,250.d0/621.d0-18575.0d0/48384.0d0,                           &
        125.0d0/594.0d0-13525.0d0/55296.0d0,-277.0d0/14336.0d0,512.0d0/1771.0d0-0.25d0/)

        IF (coupledFlag.OR.bubblesources) THEN

            !> First step
            if (proc_rank==0) print*, 'rkqs 1st step at', qtime
            call s_mpi_barrier()
            CALL update(RKh,1,RKcoef1,largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)

            !> Second step
            if (proc_rank==0) print*, 'rkqs 2nd step at', qtime+A2*RKh
            call s_mpi_barrier()
            CALL RKparticledyn(qtime+A2*RKh,2,q_cons_ts(2)%vf,t_step,q_prim_vf,rhs_vp_adapt(2)%vf,largestep)
            CALL update(RKh,2,RKcoef2,largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)

            !> Third step
            if (proc_rank==0) print*, 'rkqs 3rd step at', qtime+A3*RKh
            CALL RKparticledyn(qtime+A3*RKh,3,q_cons_ts(2)%vf,t_step,q_prim_vf,rhs_vp_adapt(3)%vf,largestep)
            CALL update(RKh,3,RKcoef3,largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)

            !> Fourth step
            if (proc_rank==0) print*, 'rkqs 4th step at', qtime+A4*RKh
            CALL RKparticledyn(qtime+A4*RKh,4,q_cons_ts(2)%vf,t_step,q_prim_vf,rhs_vp_adapt(4)%vf,largestep)
            CALL update(RKh,4,RKcoef4,largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)

            !> Fifth step
            if (proc_rank==0) print*, 'rkqs 5th step at', qtime+A5*RKh
            CALL RKparticledyn(qtime+A5*RKh,5,q_cons_ts(2)%vf,t_step,q_prim_vf,rhs_vp_adapt(5)%vf,largestep)
            CALL update(RKh,5,RKcoef5,largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)

            !> Sixth step
            if (proc_rank==0) print*, 'rkqs 6th step at', qtime+A6*RKh
            CALL RKparticledyn(qtime+A6*RKh,6,q_cons_ts(2)%vf,t_step,q_prim_vf,rhs_vp_adapt(6)%vf,largestep)
            CALL update(RKh,6,RKcoef6,largestep, q_cons_ts, rhs_vp_adapt, q_prim_vf)

            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO

            DO i = adv_idx%beg, sys_size
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO

            CALL RKerror (qtime+RKh,RKh,RKcoefE,errmax,largestep,t_step,q_cons_ts,q_prim_vf,rhs_vp_adapt) 

        ELSE

            !> First step
            CALL update(RKh,1,RKcoef1,largestep)
            IF (largestep) RETURN

            !> Second step
            CALL RKparticledyn(qtime+A2*RKh,2,q_cons_ts(1)%vf,t_step,q_prim_vf)
            CALL update(RKh,2,RKcoef2,largestep)
            IF (largestep) RETURN

            !> Third step
            CALL RKparticledyn(qtime+A3*RKh,3,q_cons_ts(1)%vf,t_step,q_prim_vf)
            CALL update(RKh,3,RKcoef3,largestep)
            IF (largestep) RETURN

            !> Fourth step
            CALL RKparticledyn (qtime+A4*RKh,4,q_cons_ts(1)%vf,t_step,q_prim_vf)
            CALL update (RKh,4,RKcoef4,largestep)
            IF (largestep) RETURN

            !> Fifth step
            CALL RKparticledyn (qtime+A5*RKh,5,q_cons_ts(1)%vf,t_step,q_prim_vf)
            CALL update (RKh,5,RKcoef5,largestep)
            IF (largestep) RETURN

            !> Sixth step
            CALL RKparticledyn (qtime+A6*RKh,6,q_cons_ts(1)%vf,t_step,q_prim_vf)
            CALL update (RKh,6,RKcoef6,largestep)
            IF (largestep) RETURN

            CALL RKerror (qtime+RKh,RKh,RKcoefE,errmax,largestep,t_step)

        ENDIF

    END SUBROUTINE rkck ! --------------------------------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_time_steppers_module() ! -------------------------

        integer :: i, j !< Generic loop iterators

        ! Deallocating the cell-average conservative variables
        do i = 1, num_ts

            do j = 1, sys_size
                @:DEALLOCATE(q_cons_ts(i)%vf(j)%sf)
            end do

            @:DEALLOCATE(q_cons_ts(i)%vf)

        end do

        @:DEALLOCATE(q_cons_ts)

        ! Deallocating the cell-average primitive ts variables
        if (probe_wrt) then
            do i = 0, 3
                do j = 1, sys_size
                    @:DEALLOCATE(q_prim_ts(i)%vf(j)%sf)
                end do
                @:DEALLOCATE(q_prim_ts(i)%vf)
            end do
            @:DEALLOCATE(q_prim_ts)
        end if

        ! Deallocating the cell-average primitive variables
        do i = 1, adv_idx%end
            @:DEALLOCATE(q_prim_vf(i)%sf)
        end do

        if (hypoelasticity) then
            do i = stress_idx%beg, stress_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        @:DEALLOCATE(q_prim_vf)

        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size
            @:DEALLOCATE(rhs_vf(i)%sf)
        end do

        @:DEALLOCATE(rhs_vf)

        ! Deallocating the cell-average RHS variable for adaptive method, Lagrangian solver
        IF(coupledflag .OR. (solverapproach.EQ.2)) THEN
            DO i = 1, 6
                DO j = 1, adv_idx%end
                    DEALLOCATE(rhs_vp_adapt(i)%vf(j)%sf)
                END DO
                DEALLOCATE(rhs_vp_adapt(i)%vf)
            END DO
            DEALLOCATE(rhs_vp_adapt)
        END IF

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module ! -----------------------

end module m_time_steppers
