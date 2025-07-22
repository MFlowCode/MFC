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

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hane-side (RHS) evaluation procedures

    use m_pressure_relaxation  !< Pressure relaxation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles_EE           !< Ensemble-averaged bubble dynamics routines

    use m_bubbles_EL           !< Lagrange bubble dynamics routines

    use m_ibm

    use m_hyperelastic

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_boundary_common

    use m_helper

    use m_sim_helpers

    use m_fftw

    use m_nvtx

    use m_thermochem, only: num_species

    use m_body_forces

    implicit none

    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(integer_field), allocatable, dimension(:, :) :: bc_type !<
    !! Boundary condition identifiers

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    real(wp), allocatable, dimension(:, :, :, :, :) :: rhs_pb

    type(scalar_field) :: q_T_sf !<
    !! Cell-average temperature variables at the current time-stage

    real(wp), allocatable, dimension(:, :, :, :, :) :: rhs_mv

    real(wp), allocatable, dimension(:, :, :) :: max_dt

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

    $:GPU_DECLARE(create='[q_cons_ts,q_prim_vf,q_T_sf,rhs_vf,q_prim_ts,rhs_mv,rhs_pb,max_dt]')

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    impure subroutine s_initialize_time_steppers_module

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (time_stepper == 1) then
            num_ts = 1
        elseif (any(time_stepper == (/2, 3/))) then
            num_ts = 2
        end if

        ! Allocating the cell-average conservative variables
        @:ALLOCATE(q_cons_ts(1:num_ts))

        do i = 1, num_ts
            @:ALLOCATE(q_cons_ts(i)%vf(1:sys_size))
        end do

        do i = 1, num_ts
            do j = 1, sys_size
                @:ALLOCATE(q_cons_ts(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
            end do
            @:ACC_SETUP_VFs(q_cons_ts(i))
        end do

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then
            @:ALLOCATE(q_prim_ts(0:3))

            do i = 0, 3
                @:ALLOCATE(q_prim_ts(i)%vf(1:sys_size))
            end do

            do i = 0, 3
                do j = 1, sys_size
                    @:ALLOCATE(q_prim_ts(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                end do
            end do

            do i = 0, 3
                @:ACC_SETUP_VFs(q_prim_ts(i))
            end do
        end if

        ! Allocating the cell-average primitive variables
        @:ALLOCATE(q_prim_vf(1:sys_size))

        if (.not. igr) then
            do i = 1, adv_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do

            if (bubbles_euler) then
                do i = bub_idx%beg, bub_idx%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
                if (adv_n) then
                    @:ALLOCATE(q_prim_vf(n_idx)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(n_idx))
                end if
            end if

            if (mhd) then
                do i = B_idx%beg, B_idx%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (elasticity) then
                do i = stress_idx%beg, stress_idx%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (hyperelasticity) then
                do i = xibeg, xiend + 1
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (cont_damage) then
                @:ALLOCATE(q_prim_vf(damage_idx)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(damage_idx))
            end if

            if (model_eqns == 3) then
                do i = internalEnergies_idx%beg, internalEnergies_idx%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (surface_tension) then
                @:ALLOCATE(q_prim_vf(c_idx)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(c_idx))
            end if

            if (chemistry) then
                do i = chemxb, chemxe
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do

                @:ALLOCATE(q_T_sf%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_T_sf)
            end if
        end if

        @:ALLOCATE(pb_ts(1:2))
        !Initialize bubble variables pb and mv at all quadrature nodes for all R0 bins
        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(pb_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
        else if (qbmm .and. polytropic) then
            @:ALLOCATE(pb_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
        else
            @:ALLOCATE(pb_ts(1)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(0,0,0,0,0))
        end if

        @:ALLOCATE(mv_ts(1:2))

        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(mv_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))

        else if (qbmm .and. polytropic) then
            @:ALLOCATE(mv_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
        else
            @:ALLOCATE(mv_ts(1)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(0,0,0,0,0))
        end if

        ! Allocating the cell-average RHS variables
        @:ALLOCATE(rhs_vf(1:sys_size))

        if (igr) then
            do i = 1, sys_size
                @:ALLOCATE(rhs_vf(i)%sf(-1:m+1,-1:n+1,-1:p+1))
                @:ACC_SETUP_SFs(rhs_vf(i))
            end do
        else
            do i = 1, sys_size
                @:ALLOCATE(rhs_vf(i)%sf(0:m, 0:n, 0:p))
                @:ACC_SETUP_SFs(rhs_vf(i))
            end do
        end if

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

        if (cfl_dt) then
            @:ALLOCATE(max_dt(0:m, 0:n, 0:p))
        end if

        ! Allocating arrays to store the bc types
        @:ALLOCATE(bc_type(1:num_dims,-1:1))

        @:ALLOCATE(bc_type(1,-1)%sf(0:0,0:n,0:p))
        @:ALLOCATE(bc_type(1,1)%sf(0:0,0:n,0:p))
        if (n > 0) then
            @:ALLOCATE(bc_type(2,-1)%sf(-buff_size:m+buff_size,0:0,0:p))
            @:ALLOCATE(bc_type(2,1)%sf(-buff_size:m+buff_size,0:0,0:p))
            if (p > 0) then
                @:ALLOCATE(bc_type(3,-1)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,0:0))
                @:ALLOCATE(bc_type(3,1)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,0:0))
            end if
        end if

        do i = 1, num_dims
            do j = -1, 1, 2
                @:ACC_SETUP_SFs(bc_type(i,j))
            end do
        end do

    end subroutine s_initialize_time_steppers_module

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    impure subroutine s_1st_order_tvd_rk(t_step, time_avg)

        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg

        integer :: i, j, k, l, q !< Generic loop iterator

        ! Stage 1 of 1
        call nvtxStartRange("TIMESTEP")

        call s_compute_rhs(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, 1)

#ifdef DEBUG
        print *, 'got rhs'
#endif

        if (run_time_info) then
            if (igr) then
                call s_write_run_time_information(q_cons_ts(1)%vf, t_step)
            else
                call s_write_run_time_information(q_prim_vf, t_step)
            end if
        end if

#ifdef DEBUG
        print *, 'wrote runtime info'
#endif

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        if (bubbles_lagrange .and. .not. adap_dt) call s_update_lagrange_tdv_rk(stage=1)

        $:GPU_PARALLEL_LOOP(collapse=4)
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
            $:GPU_PARALLEL_LOOP(collapse=5)
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
            $:GPU_PARALLEL_LOOP(collapse=5)
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

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, dt)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

    end subroutine s_1st_order_tvd_rk

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    impure subroutine s_2nd_order_tvd_rk(t_step, time_avg)

        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(wp) :: start, finish

        ! Stage 1 of 2

        call cpu_time(start)

        call nvtxStartRange("TIMESTEP")

        call s_compute_rhs(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, 1)

        if (run_time_info) then
            if (igr) then
                call s_write_run_time_information(q_cons_ts(1)%vf, t_step)
            else
                call s_write_run_time_information(q_prim_vf, t_step)
            end if
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        if (bubbles_lagrange .and. .not. adap_dt) call s_update_lagrange_tdv_rk(stage=1)

        $:GPU_PARALLEL_LOOP(collapse=4)
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
            $:GPU_PARALLEL_LOOP(collapse=5)
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
            $:GPU_PARALLEL_LOOP(collapse=5)
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

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, dt)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! Stage 2 of 2

        call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step, time_avg, 2)

        if (bubbles_lagrange .and. .not. adap_dt) call s_update_lagrange_tdv_rk(stage=2)

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/2._wp
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            $:GPU_PARALLEL_LOOP(collapse=5)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/2._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            $:GPU_PARALLEL_LOOP(collapse=5)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/2._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, 2._wp*dt/3._wp)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

        call cpu_time(finish)

    end subroutine s_2nd_order_tvd_rk

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    impure subroutine s_3rd_order_tvd_rk(t_step, time_avg)

        integer, intent(IN) :: t_step
        real(wp), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q !< Generic loop iterator

        real(wp) :: start, finish

        ! Stage 1 of 3

        if (.not. adap_dt) then
            call cpu_time(start)
            call nvtxStartRange("TIMESTEP")
        end if

        call s_compute_rhs(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, 1)

        if (run_time_info) then
            if (igr) then
                call s_write_run_time_information(q_cons_ts(1)%vf, t_step)
            else
                call s_write_run_time_information(q_prim_vf, t_step)
            end if
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        if (bubbles_lagrange .and. .not. adap_dt) call s_update_lagrange_tdv_rk(stage=1)

        $:GPU_PARALLEL_LOOP(collapse=4)
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
            $:GPU_PARALLEL_LOOP(collapse=5)
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
            $:GPU_PARALLEL_LOOP(collapse=5)
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

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, dt)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! Stage 2 of 3

        call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step, time_avg, 2)

        if (bubbles_lagrange .and. .not. adap_dt) call s_update_lagrange_tdv_rk(stage=2)

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            (3._wp*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/4._wp
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            $:GPU_PARALLEL_LOOP(collapse=5)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    (3._wp*pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/4._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            $:GPU_PARALLEL_LOOP(collapse=5)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    (3._wp*mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/4._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, dt/4._wp)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! Stage 3 of 3
        call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step, time_avg, 3)

        if (bubbles_lagrange .and. .not. adap_dt) call s_update_lagrange_tdv_rk(stage=3)

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + 2._wp*q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + 2._wp*dt*rhs_vf(i)%sf(j, k, l))/3._wp
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            $:GPU_PARALLEL_LOOP(collapse=5)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + 2._wp*pb_ts(2)%sf(j, k, l, q, i) &
                                     + 2._wp*dt*rhs_pb(j, k, l, q, i))/3._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            $:GPU_PARALLEL_LOOP(collapse=5)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + 2._wp*mv_ts(2)%sf(j, k, l, q, i) &
                                     + 2._wp*dt*rhs_mv(j, k, l, q, i))/3._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, 2._wp*dt/3._wp)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if

        call nvtxStartRange("RHS-ELASTIC")
        if (hyperelasticity) call s_hyperelastic_rmt_stress_update(q_cons_ts(1)%vf, q_prim_vf)
        call nvtxEndRange

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        if (.not. adap_dt) then
            call nvtxEndRange
            call cpu_time(finish)

            time = time + (finish - start)
        end if
    end subroutine s_3rd_order_tvd_rk

    !> Strang splitting scheme with 3rd order TVD RK time-stepping algorithm for
        !!      the flux term and adaptive time stepping algorithm for
        !!      the source term
        !! @param t_step Current time-step
    subroutine s_strang_splitting(t_step, time_avg)

        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg

        real(wp) :: start, finish

        call cpu_time(start)

        call nvtxStartRange("TIMESTEP")

        ! Stage 1 of 3
        call s_adaptive_dt_bubble(1)

        ! Stage 2 of 3
        call s_3rd_order_tvd_rk(t_step, time_avg)

        ! Stage 3 of 3
        call s_adaptive_dt_bubble(3)

        call nvtxEndRange

        call cpu_time(finish)

        time = time + (finish - start)

    end subroutine s_strang_splitting

    !> Bubble source part in Strang operator splitting scheme
        !! @param t_step Current time-step
    impure subroutine s_adaptive_dt_bubble(stage)

        integer, intent(in) :: stage

        type(vector_field) :: gm_alpha_qp

        call s_convert_conservative_to_primitive_variables( &
            q_cons_ts(1)%vf, &
            q_T_sf, &
            q_prim_vf, &
            idwint)

        if (bubbles_euler) then

            call s_compute_bubble_EE_source(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, divu)
            call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        elseif (bubbles_lagrange) then

            call s_populate_variables_buffers(bc_type, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            call s_compute_bubble_EL_dynamics(q_prim_vf, stage)
            call s_transfer_data_to_tmp()
            call s_smear_voidfraction()
            if (stage == 3) then
                if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()
                if (lag_params%write_bubbles) then
                    $:GPU_UPDATE(host='[gas_p,gas_mv,intfc_rad,intfc_vel]')
                    call s_write_lag_particles(mytime)
                end if
                call s_write_void_evol(mytime)
            end if

        end if

    end subroutine s_adaptive_dt_bubble

    impure subroutine s_compute_dt()

        real(wp) :: rho        !< Cell-avg. density
        real(wp), dimension(num_vels) :: vel        !< Cell-avg. velocity
        real(wp) :: vel_sum    !< Cell-avg. velocity sum
        real(wp) :: pres       !< Cell-avg. pressure
        real(wp), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
        real(wp) :: gamma      !< Cell-avg. sp. heat ratio
        real(wp) :: pi_inf     !< Cell-avg. liquid stiffness function
        real(wp) :: c          !< Cell-avg. sound speed
        real(wp) :: H          !< Cell-avg. enthalpy
        real(wp), dimension(2) :: Re         !< Cell-avg. Reynolds numbers
        type(vector_field) :: gm_alpha_qp

        real(wp) :: dt_local
        integer :: j, k, l !< Generic loop iterators

        if (.not. igr) then
            call s_convert_conservative_to_primitive_variables( &
                q_cons_ts(1)%vf, &
                q_T_sf, &
                q_prim_vf, &
                idwint)
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[vel, alpha, Re]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (igr) then
                        call s_compute_enthalpy(q_cons_ts(1)%vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, j, k, l)
                    else
                        call s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, j, k, l)
                    end if

                    ! Compute mixture sound speed
                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, alpha, vel_sum, 0._wp, c)

                    call s_compute_dt_from_cfl(vel, c, max_dt, rho, Re, j, k, l)
                end do
            end do
        end do

        #:call GPU_PARALLEL(copyout='[dt_local]', copyin='[max_dt]')
            dt_local = minval(max_dt)
        #:endcall GPU_PARALLEL

        if (num_procs == 1) then
            dt = dt_local
        else
            call s_mpi_allreduce_min(dt_local, dt)
        end if

        $:GPU_UPDATE(device='[dt]')

    end subroutine s_compute_dt

    !> This subroutine applies the body forces source term at each
        !! Runge-Kutta stage
    subroutine s_apply_bodyforces(q_cons_vf, q_prim_vf_in, rhs_vf_in, ldt)

        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(1:sys_size), intent(in) :: q_prim_vf_in
        type(scalar_field), dimension(1:sys_size), intent(inout) :: rhs_vf_in

        real(wp), intent(in) :: ldt !< local dt

        integer :: i, j, k, l

        call nvtxStartRange("RHS-BODYFORCES")
        call s_compute_body_forces_rhs(q_prim_vf_in, q_cons_vf, rhs_vf_in)

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = momxb, E_idx
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) + &
                                                   ldt*rhs_vf_in(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        call nvtxEndRange

    end subroutine s_apply_bodyforces

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step)

        integer, intent(in) :: t_step

        integer :: i, j, k, l !< Generic loop iterator

        if (t_step == t_step_start) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts(3)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
        elseif (t_step == t_step_start + 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts(2)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
        elseif (t_step == t_step_start + 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts(1)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
        elseif (t_step == t_step_start + 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts(0)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
        else ! All other timesteps
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts(3)%vf(i)%sf(j, k, l) = q_prim_ts(2)%vf(i)%sf(j, k, l)
                            q_prim_ts(2)%vf(i)%sf(j, k, l) = q_prim_ts(1)%vf(i)%sf(j, k, l)
                            q_prim_ts(1)%vf(i)%sf(j, k, l) = q_prim_ts(0)%vf(i)%sf(j, k, l)
                            q_prim_ts(0)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_time_step_cycling

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_time_steppers_module

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

        if (.not. igr) then
            ! Deallocating the cell-average primitive variables
            do i = 1, adv_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do

            if (mhd) then
                do i = B_idx%beg, B_idx%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (elasticity) then
                do i = stress_idx%beg, stress_idx%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (hyperelasticity) then
                do i = xibeg, xiend + 1
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (cont_damage) then
                @:DEALLOCATE(q_prim_vf(damage_idx)%sf)
            end if

            if (bubbles_euler) then
                do i = bub_idx%beg, bub_idx%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (model_eqns == 3) then
                do i = internalEnergies_idx%beg, internalEnergies_idx%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if
        end if

        @:DEALLOCATE(q_prim_vf)

        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size
            @:DEALLOCATE(rhs_vf(i)%sf)
        end do

        @:DEALLOCATE(rhs_vf)

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module

end module m_time_steppers
