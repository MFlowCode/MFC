!>
!! @file
!! @brief Contains module m_time_steppers

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Total-variation-diminishing (TVD) Runge--Kutta time integrators (1st-, 2nd-, and 3rd-order SSP)
module m_time_steppers

    use m_phase_timing
    use m_derived_types
    use m_global_parameters
    use m_rhs
    use m_chemistry
    use m_species_bounds
    use m_reactive_burn, only: s_reactive_burn_substep
    use m_pressure_relaxation
    use m_hypoelastic, only: s_enforce_cont_damage_bounds
    use m_data_output
    use m_bubbles_EE
    use m_bubbles_EL
    use m_particles_EL
    use m_ibm
    use m_collisions, only: collisions_active
    use m_mpi_proxy
    use m_boundary_common
    use m_boundary_primitives, only: f_vel_ramp, bc_vel_ramp
    use m_helper
    use m_sim_helpers
    use m_fftw
    use m_nvtx
    use m_thermochem, only: num_species
    use m_body_forces
    use m_stage_hooks, only: s_stage_hooks_pre_ib, s_stage_hooks_post_ib
    use m_derived_variables
    use m_constants, only: model_eqns_6eq, time_stepper_rk1, time_stepper_rk2, time_stepper_rk3, dflt_T_guess
    use m_amr_stage, only: s_amr_stage_begin, s_amr_stage_fine, s_amr_step_fold, s_amr_subcycle_fine, s_amr_burn_fine
    use m_amr_advance, only: s_amr_fine_stage_ib
    use m_amr_registers, only: s_amr_set_rk_table

    implicit none

    type(vector_field), allocatable, dimension(:)    :: q_cons_ts  !< Cell-average conservative variables at each time-stage (TS)
    type(scalar_field), allocatable, dimension(:)    :: q_prim_vf  !< Cell-average primitive variables at the current time-stage
    type(scalar_field), allocatable, dimension(:)    :: rhs_vf     !< Cell-average RHS variables at the current time-stage
    type(integer_field), allocatable, dimension(:,:) :: bc_type    !< Boundary condition identifiers
    !> Cell-average primitive variables at consecutive TIMESTEPS
    type(vector_field), allocatable, dimension(:) :: q_prim_ts1, q_prim_ts2
    real(wp), allocatable, dimension(:,:,:,:,:)   :: rhs_pb
    type(scalar_field)                            :: q_T_sf  !< Cell-average temperature variables at the current time-stage
    real(wp), allocatable, dimension(:,:,:,:,:)   :: rhs_mv
    integer, private                              :: num_ts  !< Number of time stages in the time-stepping scheme
    integer                                       :: stor    !< storage index
    real(wp), allocatable, dimension(:,:)         :: rk_coef
    real(wp), allocatable, dimension(:)           :: rk_c    !< stage times as fractions of dt (rk_c(nstage + 1) = 1)
    integer, private                              :: num_probe_ts
    !> Max dt*(species loss rate) of the explicit reaction source over the run
    real(wp), private :: chem_stiff_max = 0._wp
    logical, private  :: visc_cap_warned = .false.  !< The one-time notice that visc_dt_cap limits the step was printed
    logical, private  :: chem_stiff_warned = .false.

    $:GPU_DECLARE(create='[q_cons_ts, q_prim_vf, q_T_sf, rhs_vf, q_prim_ts1, q_prim_ts2, rhs_mv, rhs_pb, rk_coef, stor, bc_type]')

    !> @cond
#if defined(__NVCOMPILER_GPU_UNIFIED_MEM)
    real(stp), allocatable, dimension(:,:,:,:), pinned, target :: q_cons_ts_pool_host
#elif defined(FRONTIER_UNIFIED)
    real(stp), pointer, contiguous, dimension(:,:,:,:) :: q_cons_ts_pool_host, q_cons_ts_pool_device
    integer(kind=8)                                    :: pool_dims(4), pool_starts(4)
    integer(kind=8)                                    :: pool_size
    type(c_ptr)                                        :: cptr_host, cptr_device
#endif
    !> @endcond

contains

    !> Initialize the time steppers module
    impure subroutine s_initialize_time_steppers_module

#ifdef FRONTIER_UNIFIED
        use hipfort
        use hipfort_hipmalloc
        use hipfort_check
#if defined(MFC_OpenACC)
        use openacc
#endif
#endif
        integer :: i, j  !< Generic loop iterators

        if (time_stepper == time_stepper_rk1) then
            num_ts = 1
        else if (any(time_stepper == (/time_stepper_rk2, time_stepper_rk3/))) then
            num_ts = 2
        end if

        if (probe_wrt) then
            num_probe_ts = 2
        end if

        ! Allocating the cell-average conservative variables
        @:ALLOCATE(q_cons_ts(1:num_ts))
        @:PREFER_GPU(q_cons_ts)

        do i = 1, num_ts
            @:ALLOCATE(q_cons_ts(i)%vf(1:sys_size))
            @:PREFER_GPU(q_cons_ts(i)%vf)
        end do

        !> @cond
#if defined(__NVCOMPILER_GPU_UNIFIED_MEM)
        if (num_ts == 2 .and. nv_uvm_out_of_core) then
            ! host allocation for q_cons_ts(2)%vf(j)%sf for all j
            allocate (q_cons_ts_pool_host(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end, &
                      & idwbuff(3)%beg:idwbuff(3)%end,1:sys_size))
        end if

        do j = 1, sys_size
            ! q_cons_ts(1) lives on the device
            @:ALLOCATE(q_cons_ts(1)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                       & idwbuff(3)%beg:idwbuff(3)%end))
            @:PREFER_GPU(q_cons_ts(1)%vf(j)%sf)
            if (num_ts == 2) then
                if (nv_uvm_out_of_core) then
                    ! q_cons_ts(2) lives on the host
                    q_cons_ts(2)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end, &
                              & idwbuff(3)%beg:idwbuff(3)%end) => q_cons_ts_pool_host(:,:,:,j)
                else
                    @:ALLOCATE(q_cons_ts(2)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:PREFER_GPU(q_cons_ts(2)%vf(j)%sf)
                end if
            end if
        end do

        do i = 1, num_ts
            @:ACC_SETUP_VFs(q_cons_ts(i))
        end do
#elif defined(FRONTIER_UNIFIED)
        ! Allocate to memory regions using hip calls that we will attach pointers to
        do i = 1, 3
            pool_dims(i) = idwbuff(i)%end - idwbuff(i)%beg + 1
            pool_starts(i) = idwbuff(i)%beg
        end do
        pool_dims(4) = sys_size
        pool_starts(4) = 1
#ifdef MFC_MIXED_PRECISION
        pool_size = 1_8*(idwbuff(1)%end - idwbuff(1)%beg + 1)*(idwbuff(2)%end - idwbuff(2)%beg + 1)*(idwbuff(3)%end &
                         & - idwbuff(3)%beg + 1)*sys_size
        call hipCheck(hipMalloc_(cptr_device, pool_size*2_8))
        call c_f_pointer(cptr_device, q_cons_ts_pool_device, shape=pool_dims)
        q_cons_ts_pool_device(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:) => q_cons_ts_pool_device

        call hipCheck(hipMallocManaged_(cptr_host, pool_size*2_8, hipMemAttachGlobal))
        call c_f_pointer(cptr_host, q_cons_ts_pool_host, shape=pool_dims)
        q_cons_ts_pool_host(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:) => q_cons_ts_pool_host
#else
        ! Doing hipMalloc then mapping should be most performant
        call hipCheck(hipMalloc(q_cons_ts_pool_device, dims8=pool_dims, lbounds8=pool_starts))
        ! Without this map CCE will still create a device copy, because it's silly like that
#if defined(MFC_OpenACC)
        call acc_map_data(q_cons_ts_pool_device, c_loc(q_cons_ts_pool_device), c_sizeof(q_cons_ts_pool_device))
#endif
        ! CCE see it can access this and will leave it on the host. It will stay on the host so long as HSA_XNACK=1 NOTE: WE CANNOT
        ! DO ATOMICS INTO THIS MEMORY. We have to change a property to use atomics here Otherwise leaving this as fine-grained will
        ! actually help performance since it can't be cached in GPU L2
        if (num_ts == 2) then
            call hipCheck(hipMallocManaged(q_cons_ts_pool_host, dims8=pool_dims, lbounds8=pool_starts, flags=hipMemAttachGlobal))
#if defined(MFC_OpenMP)
            call hipCheck(hipMemAdvise(c_loc(q_cons_ts_pool_host), c_sizeof(q_cons_ts_pool_host), &
                          & hipMemAdviseSetPreferredLocation, -1))
#endif
        end if
#endif

        do j = 1, sys_size
            ! q_cons_ts(1) lives on the device
            q_cons_ts(1)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end, &
                      & idwbuff(3)%beg:idwbuff(3)%end) => q_cons_ts_pool_device(:,:,:,j)
            if (num_ts == 2) then
                ! q_cons_ts(2) lives on the host
                q_cons_ts(2)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end, &
                          & idwbuff(3)%beg:idwbuff(3)%end) => q_cons_ts_pool_host(:,:,:,j)
            end if
        end do

        do i = 1, num_ts
            @:ACC_SETUP_VFs(q_cons_ts(i))
            do j = 1, sys_size
                $:GPU_UPDATE(device='[q_cons_ts(i)%vf(j)]')
            end do
        end do
#else
        !> @endcond
        do i = 1, num_ts
            do j = 1, sys_size
                @:ALLOCATE(q_cons_ts(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
            @:ACC_SETUP_VFs(q_cons_ts(i))
        end do
        !> @cond
#endif
        !> @endcond

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then
            @:ALLOCATE(q_prim_ts1(1:num_probe_ts))

            do i = 1, num_probe_ts
                @:ALLOCATE(q_prim_ts1(i)%vf(1:sys_size))
            end do

            do i = 1, num_probe_ts
                do j = 1, sys_size
                    ! Zeroed: the acceleration probe reads all four levels from the first step, before they hold a state
                    @:ALLOCATE(q_prim_ts1(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    q_prim_ts1(i)%vf(j)%sf = 0._wp
                end do
                @:ACC_SETUP_VFs(q_prim_ts1(i))
            end do

            @:ALLOCATE(q_prim_ts2(1:num_probe_ts))

            do i = 1, num_probe_ts
                @:ALLOCATE(q_prim_ts2(i)%vf(1:sys_size))
            end do

            do i = 1, num_probe_ts
                do j = 1, sys_size
                    @:ALLOCATE(q_prim_ts2(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    q_prim_ts2(i)%vf(j)%sf = 0._wp
                end do
                @:ACC_SETUP_VFs(q_prim_ts2(i))
            end do
        end if

        ! Allocating the cell-average primitive variables
        @:ALLOCATE(q_prim_vf(1:sys_size))

        if (.not. igr) then
            do i = 1, eqn_idx%adv%end
                @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do

            if (bubbles_euler) then
                do i = eqn_idx%bub%beg, eqn_idx%bub%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
                if (adv_n) then
                    @:ALLOCATE(q_prim_vf(eqn_idx%n)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(eqn_idx%n))
                end if
            end if

            if (mhd) then
                do i = eqn_idx%B%beg, eqn_idx%B%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (hypoelasticity) then
                do i = eqn_idx%stress%beg, eqn_idx%stress%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (cont_damage) then
                @:ALLOCATE(q_prim_vf(eqn_idx%damage)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(eqn_idx%damage))
            end if

            if (hyper_cleaning) then
                @:ALLOCATE(q_prim_vf(eqn_idx%psi)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(eqn_idx%psi))
            end if

            if (model_eqns == model_eqns_6eq) then
                do i = eqn_idx%int_en%beg, eqn_idx%int_en%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (surface_tension) then
                @:ALLOCATE(q_prim_vf(eqn_idx%c)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(eqn_idx%c))
            end if

            if (chemistry) then
                do i = eqn_idx%species%beg, eqn_idx%species%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (chemistry .or. heat_conduction) then
                ! allocation bounds: q_T_sf is the one array here that crosses into the AMR fine advance, so it must hold a slab of
                ! refined blocks as well as the coarse subdomain
                @:ALLOCATE(q_T_sf%sf(idwbuff_alloc(1)%beg:idwbuff_alloc(1)%end, idwbuff_alloc(2)%beg:idwbuff_alloc(2)%end, &
                           & idwbuff_alloc(3)%beg:idwbuff_alloc(3)%end))
                ! the cache is the Newton guess of every conversion; a fine block wider than the coarse subdomain reads it at fine
                ! indices no coarse conversion wrote, so seed the whole allocation
                q_T_sf%sf = dflt_T_guess
                @:ACC_SETUP_SFs(q_T_sf)
                ! @:ALLOCATE does not copy and the ACC_SETUP copyin is Cray-only: push the seed so other offload builds do not
                ! Newton-iterate on uninitialized device memory
                $:GPU_UPDATE(device='[q_T_sf%sf]')
            end if
        end if

        @:ALLOCATE(pb_ts(1:2))
        ! Initialize bubble variables pb and mv at all quadrature nodes for all R0 bins
        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(pb_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                       & 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                       & 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                       & 1:nnode, 1:nb))
        else if (qbmm .and. polytropic) then
            @:ALLOCATE(pb_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                       & idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                       & idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(idwbuff(1)%beg:idwbuff(1)%beg + 1, idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                       & idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
        else
            @:ALLOCATE(pb_ts(1)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(0,0,0,0,0))
        end if

        @:ALLOCATE(mv_ts(1:2))

        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(mv_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                       & 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                       & 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                       & 1:nnode, 1:nb))
        else if (qbmm .and. polytropic) then
            @:ALLOCATE(mv_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                       & idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                       & idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(idwbuff(1)%beg:idwbuff(1)%beg + 1, idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                       & idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
        else
            @:ALLOCATE(mv_ts(1)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(0,0,0,0,0))
        end if

        ! Allocating the cell-average RHS variables
        @:ALLOCATE(rhs_vf(1:sys_size))
        @:PREFER_GPU(rhs_vf)

        if (igr) then
            do i = 1, sys_size
                @:ALLOCATE(rhs_vf(i)%sf(-1:m+1,-1:n+1,-1:p+1))
                @:ACC_SETUP_SFs(rhs_vf(i))
                @:PREFER_GPU(rhs_vf(i)%sf)
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

        ! Allocating arrays to store the bc types
        @:ALLOCATE(bc_type(1:num_dims,1:2))

        @:ALLOCATE(bc_type(1,1)%sf(0:0,0:n,0:p))
        @:ALLOCATE(bc_type(1,2)%sf(0:0,0:n,0:p))
        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            if (n > 0) then
                @:ALLOCATE(bc_type(2,1)%sf(-buff_size:m+buff_size,0:0,0:p))
                @:ALLOCATE(bc_type(2,2)%sf(-buff_size:m+buff_size,0:0,0:p))
                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                    if (p > 0) then
                        @:ALLOCATE(bc_type(3,1)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,0:0))
                        @:ALLOCATE(bc_type(3,2)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,0:0))
                    end if
                #:endif
            end if
        #:endif

        do i = 1, num_dims
            do j = 1, 2
                @:ACC_SETUP_SFs(bc_type(i,j))
            end do
        end do

        if (any(time_stepper == (/time_stepper_rk1, time_stepper_rk2, time_stepper_rk3/))) then
            ! temporary array index for TVD RK
            if (time_stepper == time_stepper_rk1) then
                stor = 1
            else
                stor = 2
            end if

            ! TVD RK coefficients
            @:ALLOCATE(rk_coef(time_stepper, 4))
            if (time_stepper == time_stepper_rk1) then
                rk_coef(1,:) = (/1._wp, 0._wp, 1._wp, 1._wp/)
            else if (time_stepper == time_stepper_rk2) then
                rk_coef(1,:) = (/1._wp, 0._wp, 1._wp, 1._wp/)
                rk_coef(2,:) = (/1._wp, 1._wp, 1._wp, 2._wp/)
            else if (time_stepper == time_stepper_rk3) then
                rk_coef(1,:) = (/1._wp, 0._wp, 1._wp, 1._wp/)
                rk_coef(2,:) = (/1._wp, 3._wp, 1._wp, 4._wp/)
                rk_coef(3,:) = (/2._wp, 1._wp, 2._wp, 3._wp/)
            end if
            $:GPU_UPDATE(device='[rk_coef, stor]')
            ! Stage times as fractions of dt: stage s evaluates its RHS at rk_c(s) and produces the state at rk_c(s + 1)
            allocate (rk_c(time_stepper + 1))
            rk_c(1) = 0._wp
            do i = 1, time_stepper
                rk_c(i + 1) = (rk_coef(i, 1)*rk_c(i) + rk_coef(i, 3))/rk_coef(i, 4)
            end do
            if (amr) call s_amr_set_rk_table(rk_coef, rk_c)
        end if

    end subroutine s_initialize_time_steppers_module

    !> Advance the solution one full step using a TVD Runge-Kutta time integrator
    impure subroutine s_tvd_rk(t_step, time_avg, nstage)

        $:CRAY_DIR('OPTIMIZE (-haggress)')
        integer, intent(in)     :: t_step
        real(wp), intent(inout) :: time_avg
        integer, intent(in)     :: nstage
        integer                 :: i, j, k, l, q, s  !< Generic loop iterator
        integer(kind=8)         :: step_t0, step_t1
        real(wp)                :: a_w               !< weight of this stage's state in the step's result (swept-cell booking)
        integer(kind=8)         :: stage_t0, stage_t1, clock_rate, clock_max
        real(wp)                :: stage_time
        integer, parameter      :: n_warmup = 2      !< time steps excluded before the timing floor (warmup/JIT/first-touch)

        call system_clock(step_t0)
        call nvtxStartRange("TIMESTEP")

        ! Adaptive dt: initial stage
        if (adap_dt) call s_adaptive_dt_bubble(1)
        amr_sub_on = amr .and. amr_subcycle  ! level 0's stages capture its dt-weighted face fluxes

        do s = 1, nstage
            call system_clock(stage_t0)
            ! Explicitly time-dependent terms (sources, forcing, inflow ramps) are evaluated at the stage's own time; the
            ! step-start time in every stage makes them first order. The GRCBC ramp reads t_rhs on the device.
            t_rhs = f_step_time(rk_c(s), t_step, dt)
            $:GPU_UPDATE(device='[t_rhs]')
            ! Dirichlet inflows ramp too; their ghost fill lives in common code, which cannot read t_rhs
            if (any([bc_x%vel_in_ramp, bc_y%vel_in_ramp, bc_z%vel_in_ramp] > 0._wp)) then
                bc_vel_ramp = [f_vel_ramp(bc_x%vel_in_ramp, bc_x%vel_in_t0, bc_x%vel_in_frac0, t_rhs), &
                                          & f_vel_ramp(bc_y%vel_in_ramp, bc_y%vel_in_t0, bc_y%vel_in_frac0, t_rhs), &
                                          & f_vel_ramp(bc_z%vel_in_ramp, bc_z%vel_in_t0, bc_z%vel_in_frac0, t_rhs)]
                $:GPU_UPDATE(device='[bc_vel_ramp]')
            end if
            call s_amr_stage_begin(q_cons_ts(1)%vf)
            call s_phase_tic(PH_COARSE)
            call s_compute_rhs(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, &
                               & t_step, s)
            call s_phase_toc(PH_COARSE)

            if (s == 1) then
                if (run_time_info) then
                    if (igr) then
                        call s_write_run_time_information(q_cons_ts(1)%vf, t_step)
                    end if
                    if (.not. igr) then
                        call s_write_run_time_information(q_prim_vf, t_step)
                    end if
                end if

                if (probe_wrt) then
                    call s_time_step_cycling(t_step)
                    call s_compute_derived_variables(t_step, q_cons_ts(1)%vf, q_prim_ts1, q_prim_ts2)
                end if

                call s_write_ib_history(t_step)

                if (merge(mytime >= t_stop, t_step == t_step_stop, cfl_dt)) then
                    amr_sub_on = .false.
                    return
                end if
            end if

            if (s == 1 .and. bodyForces) call s_update_mean_sources(q_cons_ts(1)%vf, q_T_sf, t_step)

            ! subcycled: level 0 steps alone; the fine levels follow after its step (s_amr_subcycle_fine)
            if (.not. amr_sub_on) call s_amr_stage_fine(s, t_step, rk_coef(s,:), q_cons_ts(1)%vf, q_T_sf, bc_type, rhs_vf, &
                & pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv)
            ! AMR: move the bodies once the fine frames have added their faces to this stage's forces, then run the fine IB
            ! rebuild and correction the fine stage left waiting, at this stage's pose as the coarse ones are; the fine
            ! containment check reads the pose and rotation on the host
            if (amr .and. moving_immersed_boundary_flag) then
                call s_advance_ib_pose(s, t_step)
                call s_update_ib_rotation_matrices(num_ibs)
                $:GPU_UPDATE(host='[patch_ib(1:num_ibs)]')
                call s_amr_fine_stage_ib(pb_ts(1)%sf, mv_ts(1)%sf)
            end if
            if (bubbles_lagrange .and. .not. adap_dt) then
                call s_update_lagrange_tdv_rk(q_prim_vf, bc_type, stage=s, coef=rk_coef(s,:))
                if (s == time_stepper) call s_write_lag_step_output()
            end if
            if (particles_lagrange) then
                call s_write_lag_particle_evol_if_due(s)
                call s_update_lagrange_particles_tdv_rk(bc_type, stage=s, coef=rk_coef(s,:))
                if (s == time_stepper) call s_write_lag_step_output()
            end if
            ! TWIN of s_amr_fine_rk_update_batch (m_amr_advance): change the stage algebra here and it must follow
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            if (s == 1 .and. nstage > 1) then
                                q_cons_ts(stor)%vf(i)%sf(j, k, l) = q_cons_ts(1)%vf(i)%sf(j, k, l)
                            end if
                            if (igr) then
                                q_cons_ts(1)%vf(i)%sf(j, k, l) = (rk_coef(s, 1)*q_cons_ts(1)%vf(i)%sf(j, k, l) + rk_coef(s, &
                                          & 2)*q_cons_ts(stor)%vf(i)%sf(j, k, l) + rk_coef(s, 3)*rhs_vf(i)%sf(j, k, &
                                          & l))/rk_coef(s, 4)
                            else
                                q_cons_ts(1)%vf(i)%sf(j, k, l) = (rk_coef(s, 1)*q_cons_ts(1)%vf(i)%sf(j, k, l) + rk_coef(s, &
                                          & 2)*q_cons_ts(stor)%vf(i)%sf(j, k, l) + rk_coef(s, 3)*dt*rhs_vf(i)%sf(j, k, &
                                          & l))/rk_coef(s, 4)
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            ! Evolve pb and mv for non-polytropic qbmm
            if (qbmm .and. (.not. polytropic)) then
                $:GPU_PARALLEL_LOOP(collapse=5)
                do i = 1, nb
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                do q = 1, nnode
                                    if (s == 1 .and. nstage > 1) then
                                        pb_ts(stor)%sf(j, k, l, q, i) = pb_ts(1)%sf(j, k, l, q, i)
                                        mv_ts(stor)%sf(j, k, l, q, i) = mv_ts(1)%sf(j, k, l, q, i)
                                    end if
                                    pb_ts(1)%sf(j, k, l, q, i) = (rk_coef(s, 1)*pb_ts(1)%sf(j, k, l, q, i) + rk_coef(s, &
                                          & 2)*pb_ts(stor)%sf(j, k, l, q, i) + rk_coef(s, 3)*dt*rhs_pb(j, k, l, q, i))/rk_coef(s, 4)
                                    mv_ts(1)%sf(j, k, l, q, i) = (rk_coef(s, 1)*mv_ts(1)%sf(j, k, l, q, i) + rk_coef(s, &
                                          & 2)*mv_ts(stor)%sf(j, k, l, q, i) + rk_coef(s, 3)*dt*rhs_mv(j, k, l, q, i))/rk_coef(s, 4)
                                end do
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            call s_stage_hooks_pre_ib(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, rk_coef(s, 3)*dt/rk_coef(s, 4), .false., q_T_sf)

            if (ib) then
                ! check if any IBMS are moving, and if so, update the markers, ghost points, levelsets, and levelset norms
                if (moving_immersed_boundary_flag) then
                    if (.not. amr) call s_advance_ib_pose(s, t_step)
                    call s_update_mib(num_ibs)
                    ! the cells this placement covered or uncovered, before the ghost correction
                    if (ib_exchange .or. ib_uncover_fill) then
                        a_w = 1._wp
                        do q = s + 1, nstage
                            a_w = a_w*rk_coef(q, 1)/rk_coef(q, 4)
                        end do
                        if (ib_uncover_fill) call s_fill_uncovered_cells(q_cons_ts(1)%vf, q_cons_ts(merge(stor, 1, &
                            & nstage > 1))%vf, q_prim_vf)
                        call s_add_swept_cells(q_cons_ts(1)%vf, q_cons_ts(merge(stor, 1, nstage > 1))%vf, a_w, dt)
                    end if
                end if

                ! update the ghost fluid properties point values based on IB state
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf, 0)

                ! Wall columns advance once per step, from the first stage's RHS state, and set each stage's conduction temperature
                if (ib_wall_col) then
                    if (s == 1) call s_ibm_advance_wall_columns(q_cons_ts(1)%vf, q_prim_vf, q_T_sf, dt)
                    call s_wall_column_tcond(q_cons_ts(1)%vf, q_prim_vf, q_T_sf)
                end if
            end if

            call s_stage_hooks_post_ib(q_cons_ts(1)%vf, 0)

            ! Grind: minimum wall-clock time of a full RK stage (compute + halo H2D/D2H +
            ! update + IBM correction, aside from I/O) over steady-state stages. Wall clock
            ! (not cpu_time, which counts MPI spin-wait) and the minimum (jitter only adds
            ! time) make it reproducible; the min over stages drops the I/O-bearing s==1 stage.
            call system_clock(stage_t1, clock_rate, clock_max)
            stage_time = f_clock_seconds(stage_t0, stage_t1, clock_rate, clock_max)
            if (t_step - t_step_start < n_warmup) then
                time_avg = 0._wp
            else if (time_avg <= 0._wp) then
                time_avg = stage_time
            else
                time_avg = min(time_avg, stage_time)
            end if
        end do

        ! Operator-split reaction: integrate the stiff chemistry ODE per cell after the flow update,
        ! with sub-stepping, instead of adding the reaction source to the flow RHS (chem_params%reaction_substeps > 0).
        if (chemistry .and. chem_params%reactions .and. chem_params%reaction_substeps > 0) then
            call nvtxStartRange("CHEM-REACTION-SUBSTEP")
            call s_chemistry_reaction_substep(q_cons_ts(1)%vf, q_T_sf, dt, idwint)
            call nvtxEndRange
        end if

        if (chemistry) call s_report_species_cleanup(t_step)
        call s_report_chem_stiffness(t_step)

        ! Operator-split condensed-phase burn: integrate the progress variable per cell after the flow
        ! update, with sub-stepping, instead of adding the source to the flow RHS (rburn%substeps > 0).
        if (reactive_burn .and. rburn%substeps > 0) then
            call nvtxStartRange("BURN-SUBSTEP")
            call s_reactive_burn_substep(q_cons_ts(1)%vf, dt, idwint)
            call nvtxEndRange
        end if

        if (amr_sub_on) call s_amr_subcycle_fine(q_cons_ts(1)%vf, q_cons_ts(stor)%vf, rk_coef, q_T_sf, bc_type, pb_ts(1)%sf, &
            & rhs_pb, mv_ts(1)%sf, rhs_mv, t_step)
        amr_sub_on = .false.
        if (amr .and. reactive_burn .and. rburn%substeps > 0) call s_amr_burn_fine()
        if (amr) call s_amr_step_fold(q_cons_ts(1)%vf)

        if (ib) then
            ! Before the handoff: ib_heat and ib_shell are indexed by the patch_ib slots the stages used
            if (ib_exchange) call s_update_ib_temperatures(dt, q_cons_ts(1)%vf)
            if (ib_mass_relax > 0) call s_conserve_ib_mass(q_cons_ts(1)%vf)
            if (moving_immersed_boundary_flag) then
                call s_wrap_periodic_ibs()  ! wraps the positions of IBs to the local proc
                call s_handoff_ib_ownership()  ! recomputes which ranks own which IBs and communicate to neighbors
            else if (ib_state_wrt .or. ib_force_wrt) then
                call s_compute_ib_forces()
            end if
        end if

        ! Adaptive dt: final stage
        if (adap_dt) call s_adaptive_dt_bubble(3)

        call nvtxEndRange
        ! wall clock: cpu_time also counts host-thread and MPI spin time
        call system_clock(step_t1, clock_rate, clock_max)
        wall_time = f_clock_seconds(step_t0, step_t1, clock_rate, clock_max)

        if (t_step - t_step_start >= 2) then
            wall_time_avg = (wall_time + (t_step - t_step_start - 2)*wall_time_avg)/(t_step - t_step_start - 1)
        else
            wall_time_avg = 0._wp
        end if

    end subroutine s_tvd_rk

    !> Seconds between two system_clock counts, correcting for the counter wrapping past cmax.
    pure function f_clock_seconds(t0, t1, rate, cmax) result(sec)

        integer(kind=8), intent(in) :: t0, t1, rate, cmax
        real(wp)                    :: sec

        sec = real(merge(t1 - t0, t1 - t0 + cmax + 1_8, t1 >= t0), wp)/real(rate, wp)

    end function f_clock_seconds

    !> Bubble source part in Strang operator splitting scheme
    impure subroutine s_adaptive_dt_bubble(stage)

        integer, intent(in) :: stage

        ! The interior only: before the first RHS the conservative ghosts are unset (a zero bubble number density
        ! divides by zero there). The Euler source reads no ghost cell; Lagrange refills the primitive ones below.

        call s_convert_conservative_to_primitive_variables(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, idwint)

        if (bubbles_euler) then
            call s_compute_bubble_EE_source(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, divu)
            call s_comp_alpha_from_n(q_cons_ts(1)%vf)
        else if (bubbles_lagrange) then
            call s_populate_variables_buffers(bc_type, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf, q_T_sf)
            call s_compute_bubble_EL_dynamics(q_prim_vf, bc_type)
            if (stage == 3) then
                if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()
                call s_write_lag_step_output()
            end if
        end if

    end subroutine s_adaptive_dt_bubble

    !> Compute the global time step size from CFL stability constraints across all cells
    !> Sets dt from the CFL criteria, each limited by the finest cell over the coarse cell under AMR. With dt_cand it only returns
    !! the global candidates at cfl_target (ICFL, VCFL, CCFL, TCFL, collision cap, WCOL) and leaves dt alone.
    impure subroutine s_compute_dt(dt_cand)

        real(wp), dimension(6), intent(out), optional :: dt_cand
        real(wp), dimension(6)                        :: dt_candidates_glb
        real(wp)                                      :: dt_prev
        real(wp)                                      :: rho  !< Cell-avg. density

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: vel               !< Cell-avg. velocity
            real(wp), dimension(3) :: alpha, alpha_rho  !< Cell-avg. volume fraction, partial density
        #:else
            real(wp), dimension(num_vels)   :: vel               !< Cell-avg. velocity
            real(wp), dimension(num_fluids) :: alpha, alpha_rho  !< Cell-avg. volume fraction, partial density
        #:endif
        real(wp) :: vel_sum                          !< Cell-avg. velocity sum
        real(wp) :: pres                             !< Cell-avg. pressure
        real(wp) :: gamma                            !< Cell-avg. sp. heat ratio
        real(wp) :: pi_inf                           !< Cell-avg. liquid stiffness function
        real(wp) :: qv                               !< Cell-avg. fluid reference energy
        real(wp) :: c                                !< Cell-avg. sound speed
        real(wp), dimension(2) :: Re                 !< Cell-avg. Reynolds numbers
        real(wp) :: Dth                              !< Cell thermal diffusivity (s_compute_cell_diffusivity)
        real(wp), dimension(4) :: max_dt             !< Cell dt candidates (inviscid, viscous, capillary, thermal)
        real(wp) :: icfl_dt_local, vcfl_dt_local, ccfl_dt_local, tcfl_dt_local, coll_dt_local
        real(wp), dimension(6) :: dt_candidates_loc  !< Rank-local dt candidates (ICFL, VCFL, CCFL, TCFL, collision cap, WCOL)
        real(wp) :: rf, pw1, pw2, pw3                !< amr_cfl_rf of the cell, amr_cfl_pow (scalars for the kernel)
        real(wp), dimension(num_species) :: Ys       !< Cell mass fractions, for chem_params%transport_dt
        real(wp), dimension(2) :: Re_chem            !< Transport-limited Reynolds numbers, for chem_params%transport_dt
        real(wp) :: T_cell                           !< Cell temperature, passed as a scalar (array elements are misaddressed)
        logical :: is_fluid_cell                     !< Cell lies outside every immersed boundary
        integer :: bad_cell                          !< 1 + linear index of a fluid cell with no valid dt, 0 if none
        integer :: j, k, l                           !< Generic loop iterators
        integer :: fl                                !< Fluid loop iterator

        if (.not. igr) then
            call s_convert_conservative_to_primitive_variables(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, idwint)
        end if

        dt_prev = dt
        icfl_dt_local = huge(1.0_wp)
        vcfl_dt_local = huge(1.0_wp)
        ccfl_dt_local = huge(1.0_wp)
        tcfl_dt_local = huge(1.0_wp)
        coll_dt_local = huge(1.0_wp)
        bad_cell = 0
        pw1 = amr_cfl_pow(1); pw2 = amr_cfl_pow(2); pw3 = amr_cfl_pow(3)
        ! One loop body, two kernels: without amr it is the plain CFL loop, so the non-AMR dt does not depend on the AMR terms
        #:for ON_AMR in [False, True]
            if (amr .eqv. ${'.true.' if ON_AMR else '.false.'}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, firstprivate=('[pw1, pw2, pw3]' if ON_AMR else None), private='[vel, alpha, &
                                    & alpha_rho, Re, Dth, rho, vel_sum, pres, gamma, pi_inf, c, qv, fl, max_dt, is_fluid_cell, &
                                    & Ys, Re_chem, T_cell' + (', rf' if ON_AMR else '') + ']', reduction='[[icfl_dt_local, &
                                    & vcfl_dt_local, ccfl_dt_local, tcfl_dt_local], [bad_cell]]', reductionOp='[min, max]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            ! Cells inside an immersed boundary hold ghost-derived, non-physical state and must not set the global
                            ! dt.
                            is_fluid_cell = .true.
                            if (ib) is_fluid_cell = (ib_markers%sf(j, k, l) == 0)

                            if (is_fluid_cell) then
                                if (igr) then
                                    call s_compute_cell_state(q_cons_ts(1)%vf, pres, rho, gamma, pi_inf, Re, alpha, alpha_rho, &
                                                              & vel, vel_sum, qv, j, k, l)
                                else
                                    call s_compute_cell_state(q_prim_vf, pres, rho, gamma, pi_inf, Re, alpha, alpha_rho, vel, &
                                                              & vel_sum, qv, j, k, l)
                                end if

                                ! Compute mixture sound speed
                                call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, alpha, c, alpha_rho)

                                if (any_non_newtonian) then
                                    Re(1) = 0._wp
                                    do fl = 1, num_fluids
                                        if (is_non_newtonian(fl)) then
                                            Re(1) = Re(1) + alpha(fl)*hb_mu_max(fl)
                                        else
                                            Re(1) = Re(1) + alpha(fl)*fluid_inv_re(fl)
                                        end if
                                    end do
                                    Re(1) = 1._wp/max(Re(1), sgm_eps)
                                end if

                                if (chemistry .and. chem_params%transport_dt) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do fl = 1, num_species
                                        Ys(fl) = q_prim_vf(eqn_idx%species%beg + fl - 1)%sf(j, k, l)
                                    end do
                                    T_cell = q_T_sf%sf(j, k, l)
                                    call s_compute_transport_dt_re(pres, rho, T_cell, Ys, Re_chem)
                                    ! Keep the bounds Cantera does not replace: the HB mu_max in Re(1) and bulk viscosity in Re(2).
                                    ! Without viscous, Re holds dflt_real, not a bound
                                    if (viscous) then
                                        if (any_non_newtonian) Re_chem(1) = min(Re_chem(1), Re(1))
                                        Re_chem(2) = min(Re_chem(2), Re(2))
                                    end if
                                    Re = Re_chem
                                end if

                                call s_compute_cell_diffusivity(alpha, alpha_rho, Dth)
                                call s_compute_dt_from_cfl(vel, c, max_dt, rho, Re, alpha, alpha_rho, Dth, j, k, l)
                                #:if ON_AMR
                                    rf = amr_cfl_rf(j, k, l)
                                    max_dt(1) = max_dt(1)/rf**pw1
                                    max_dt(2) = max_dt(2)/rf**pw2
                                    max_dt(3) = max_dt(3)/rf**pw3
                                    max_dt(4) = max_dt(4)/rf**pw2
                                #:endif

                                ! min() drops a NaN candidate, so a cell with a NaN or non-positive dt would otherwise vanish from
                                ! the reduction and dt would come from the cells still valid.
                                if (.not. (max_dt(1) > 0._wp)) bad_cell = max(bad_cell, 1 + j + (m + 1)*(k + (n + 1)*l))

                                icfl_dt_local = min(icfl_dt_local, max_dt(1))
                                vcfl_dt_local = min(vcfl_dt_local, max_dt(2))
                                ccfl_dt_local = min(ccfl_dt_local, max_dt(3))
                                tcfl_dt_local = min(tcfl_dt_local, max_dt(4))
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        if (bad_cell > 0) call s_abort_on_bad_dt_cell(bad_cell - 1)

        ! restrict the time step so an ongoing collision spans at least collision_temporal_resolution time steps; the collision
        ! flag is rank-local, so the cap enters as a candidate before the global elementwise min propagates it to all ranks
        if (collision_model > 0 .and. collision_temporal_resolution > 0) then
            if (collisions_active) coll_dt_local = collision_time/real(collision_temporal_resolution, wp)
            collisions_active = .false.
        end if

        dt_candidates_loc(1) = icfl_dt_local
        dt_candidates_loc(2) = vcfl_dt_local
        dt_candidates_loc(3) = ccfl_dt_local
        dt_candidates_loc(4) = tcfl_dt_local
        dt_candidates_loc(5) = coll_dt_local
        dt_candidates_loc(6) = huge(1.0_wp)
        if (ib_wall_col) dt_candidates_loc(6) = f_wall_column_dt(q_prim_vf, q_T_sf)

        if (num_procs == 1) then
            dt_candidates_glb = dt_candidates_loc
        else
            call s_mpi_allreduce_min_vec(dt_candidates_loc, dt_candidates_glb)
        end if

        if (present(dt_cand)) then
            dt_cand = dt_candidates_glb
            return
        end if

        dt = minval(dt_candidates_glb)
        dt_limiter = dt_limiter_names(minloc(dt_candidates_glb, dim=1))
        if (dt_limiter == 'VCFL' .and. cfl_target > visc_dt_cap .and. .not. visc_cap_warned) then
            visc_cap_warned = .true.
            if (proc_rank == 0) print '(A,ES11.4,A,F5.2,A,I0,A,F5.2,A)', ' WARNING: viscous step limited to dt = ', dt, &
                & ' (VCFL = ', visc_dt_cap, ', stability limit of ' // trim(visc_dt_stencil) // ' in ', num_dims, &
                & 'D); cfl_target = ', cfl_target, ' exceeds it'
        end if

        ! limit how much the time step can grow relative to the previous step
        if (ramp_ratio > 0._wp .and. dt_prev > 0._wp .and. ramp_ratio*dt_prev < dt) then
            dt = ramp_ratio*dt_prev
            dt_limiter = 'RAMP'
        end if

        $:GPU_UPDATE(device='[dt]')

    end subroutine s_compute_dt

    !> AMR with a fixed dt: abort if dt exceeds the finest level's inviscid limit (ICFL = 1, as s_write_run_time_information).
    !! reach: test every cell as if refined to amr_max_level, which a regrid may do later, so the failure is at start-up.
    impure subroutine s_amr_check_fixed_dt(reach)

        logical, intent(in), optional :: reach
        real(wp), dimension(6)        :: dt_cand
        real(wp), allocatable         :: rf_save(:,:,:)
        logical                       :: all_fine
        character(len=256)            :: msg

        ! the cell state and signal speed behind the candidates are the non-relativistic ones

        if (cfl_dt .or. .not. amr .or. relativity) return
        all_fine = .false.; if (present(reach)) all_fine = reach
        if (all_fine) then
            rf_save = amr_cfl_rf; amr_cfl_rf = real(amr_ref_ratio**amr_max_level, wp)
            $:GPU_UPDATE(device='[amr_cfl_rf]')
        end if
        call s_compute_dt_limits(dt_cand)
        if (all_fine) then
            amr_cfl_rf = rf_save
            $:GPU_UPDATE(device='[amr_cfl_rf]')
        end if
        if (dt <= dt_cand(1)) return
        if (all_fine) then
            write (msg, '(A,ES11.4,A,I0,A,ES11.4,A)') 'amr: dt = ', dt, &
                   & ' exceeds the inviscid CFL limit a regrid would set by refining the fastest cells to level ', amr_max_level, &
                   & '; the largest stable dt (ICFL = 1) there is ', dt_cand(1), '. Reduce dt, or set cfl_adap_dt or amr_subcycle.'
        else
            write (msg, '(A,ES11.4,A,ES11.4,A)') 'amr: dt = ', dt, &
                   & ' exceeds the inviscid CFL limit of the finest cells; the largest ' // 'stable dt (ICFL = 1) is ', &
                   & dt_cand(1), '. Reduce dt, or set cfl_adap_dt or amr_subcycle.'
        end if
        call s_mpi_abort(trim(msg))

    end subroutine s_amr_check_fixed_dt

    !> The dt candidates of s_compute_dt at cfl_target = 1: the stability limits of the current state
    impure subroutine s_compute_dt_limits(dt_cand)

        real(wp), dimension(6), intent(out) :: dt_cand
        real(wp)                            :: cfl_save

        cfl_save = cfl_target; cfl_target = 1._wp
        $:GPU_UPDATE(device='[cfl_target]')
        call s_compute_dt(dt_cand)
        cfl_target = cfl_save
        $:GPU_UPDATE(device='[cfl_target]')

    end subroutine s_compute_dt_limits

    !> Fixed dt with viscosity: warn at start-up if dt exceeds the viscous stability limit (visc_dt_cap) of the stencil. The
    !! chemistry temperature is the Newton guess of the next conversion, so it is restored and the run stays bitwise unchanged.
    impure subroutine s_check_fixed_dt_viscous()

        real(wp), dimension(6) :: dt_cand
        real(stp), allocatable :: T_save(:,:,:)

        if (cfl_dt .or. .not. viscous .or. relativity) return
        if (chemistry) then
            $:GPU_UPDATE(host='[q_T_sf%sf]')
            T_save = q_T_sf%sf
        end if
        call s_compute_dt_limits(dt_cand)
        if (chemistry) then
            q_T_sf%sf = T_save
            $:GPU_UPDATE(device='[q_T_sf%sf]')
        end if
        if (proc_rank == 0 .and. dt > dt_cand(2)) print '(A,ES11.4,A,F7.3,A,I0,A,ES11.4,A)', ' WARNING: dt = ', dt, ' is ', &
            & dt/dt_cand(2), ' times the viscous stability limit of ' // trim(visc_dt_stencil) // ' in ', num_dims, &
            & 'D (stable up to dt = ', dt_cand(2), '); the run may diverge. Reduce dt or set cfl_adap_dt.'

    end subroutine s_check_fixed_dt_viscous

    !> Abort with the location of a fluid cell whose state admits no valid time step (NaN or non-positive inviscid dt).
    impure subroutine s_abort_on_bad_dt_cell(idx)

        integer, intent(in) :: idx  !< Rank-local linear index j + (m+1)*(k + (n+1)*l)
        integer             :: j, k, l
        real(wp)            :: loc(3)
        character(len=256)  :: msg

        j = mod(idx, m + 1)
        k = mod(idx/(m + 1), n + 1)
        l = idx/((m + 1)*(n + 1))
        loc = 0._wp
        loc(1) = x_cc(j)
        if (n > 0) loc(2) = y_cc(k)
        if (p > 0) loc(3) = z_cc(l)
        write (msg, '(A,I0,A,3(I0,1X),A,3(ES12.4,1X),A)') 'No valid time step: rank ', proc_rank, ', local cell (j,k,l) = ', j, &
               & k, l, ', x = ', loc, '. Its state is NaN or unphysical.'
        call s_mpi_abort(trim(msg))

    end subroutine s_abort_on_bad_dt_cell

    !> Advance the immersed boundary positions and velocities to the current Runge-Kutta stage (forces from this stage's RHS)
    subroutine s_advance_ib_pose(s, t_step)

        integer, intent(in)    :: s, t_step
        integer                :: i
        integer                :: gbl_id  ! used for analytic ib patch motion
        real(wp)               :: t_stage  ! time of the state produced by RK stage s (used by prescribed kinematics)
        real(wp), dimension(3) :: vel_new, angular_new

        call nvtxStartRange("PROPAGATE-IMMERSED-BOUNDARIES")

        if (moving_immersed_boundary_flag) call s_compute_ib_forces()

        t_stage = f_step_time(rk_c(s + 1), t_step, dt)

        $:GPU_PARALLEL_LOOP(private='[i, gbl_id, vel_new, angular_new]', copyin='[s, t_stage]')
        do i = 1, num_ibs
            if (s == 1) then
                patch_ib(i)%step_vel = patch_ib(i)%vel
                patch_ib(i)%step_angular_vel = patch_ib(i)%angular_vel
                patch_ib(i)%step_angles = patch_ib(i)%angles
                patch_ib(i)%step_x_centroid = patch_ib(i)%x_centroid
                patch_ib(i)%step_y_centroid = patch_ib(i)%y_centroid
                patch_ib(i)%step_z_centroid = patch_ib(i)%z_centroid
            end if

            ! Compute forces BEFORE the RK velocity blend so the device copy of patch_ib%vel matches the host (pre-blend) when
            ! velocity-dependent collision damping forces are evaluated on the GPU.
            if (patch_ib(i)%moving_ibm > 0) then
                if (patch_ib(i)%kin_model > 0) then
                    call s_prescribed_kinematics(i, t_stage)
                else
                    vel_new = (rk_coef(s, 1)*patch_ib(i)%vel + rk_coef(s, 2)*patch_ib(i)%step_vel)/rk_coef(s, 4)
                    angular_new = (rk_coef(s, 1)*patch_ib(i)%angular_vel + rk_coef(s, 2)*patch_ib(i)%step_angular_vel)/rk_coef(s, 4)

                    if (patch_ib(i)%moving_ibm >= 2) then  ! if we are using two-way coupling, apply force and torque
                        ! update the velocity from the force value
                        vel_new = vel_new + rk_coef(s, 3)*dt*(patch_ib(i)%force/patch_ib(i)%mass)/rk_coef(s, 4)

                        ! update the angular velocity with the torque value add the torque to the angular momentum
                        angular_new = (angular_new*patch_ib(i)%moment) + (rk_coef(s, 3)*dt*patch_ib(i)%torque/rk_coef(s, 4))
                        ! update the moment of inertia to be based on the direction of the angular momentum
                        if (num_dims == 3) call s_compute_moment_of_inertia(patch_ib(i), angular_new, patch_ib(i)%moment)
                        angular_new = angular_new/patch_ib(i)%moment
                    end if

                    ! Update the angle and position from the current-stage velocities, before they are overwritten
                    patch_ib(i)%angles = (rk_coef(s, 1)*patch_ib(i)%angles + rk_coef(s, 2)*patch_ib(i)%step_angles + rk_coef(s, &
                             & 3)*patch_ib(i)%angular_vel*dt)/rk_coef(s, 4)
                    patch_ib(i)%x_centroid = (rk_coef(s, 1)*patch_ib(i)%x_centroid + rk_coef(s, &
                             & 2)*patch_ib(i)%step_x_centroid + rk_coef(s, 3)*patch_ib(i)%vel(1)*dt)/rk_coef(s, 4)
                    patch_ib(i)%y_centroid = (rk_coef(s, 1)*patch_ib(i)%y_centroid + rk_coef(s, &
                             & 2)*patch_ib(i)%step_y_centroid + rk_coef(s, 3)*patch_ib(i)%vel(2)*dt)/rk_coef(s, 4)
                    patch_ib(i)%z_centroid = (rk_coef(s, 1)*patch_ib(i)%z_centroid + rk_coef(s, &
                             & 2)*patch_ib(i)%step_z_centroid + rk_coef(s, 3)*patch_ib(i)%vel(3)*dt)/rk_coef(s, 4)

                    patch_ib(i)%vel = vel_new
                    patch_ib(i)%angular_vel = angular_new

                    ! analytic velocities for 1-way coupling, at t_stage: the time of the state this stage produces
                    if (patch_ib(i)%moving_ibm == 1) then
                        @:mib_analytical()
                    end if
                end if
            end if

            ! Wall-point acceleration for the ghost pressure: the change over the step, taken at its first stage. Not at step 0:
            ! step_vel is the unset initial velocity there, so a prescribed velocity nonzero at t = 0 would read as v/dt
            if (s == 1 .and. t_step > 0 .and. patch_ib(i)%moving_ibm > 0) then
                patch_ib(i)%accel = (patch_ib(i)%vel - patch_ib(i)%step_vel)/dt
                patch_ib(i)%angular_accel = (patch_ib(i)%angular_vel - patch_ib(i)%step_angular_vel)/dt
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_advance_ib_pose

    !> Save the temporary q_prim_vf vector into q_prim_ts for use in p_main
    subroutine s_time_step_cycling(t_step)

        integer, intent(in) :: t_step
        integer             :: i, j, k, l  !< Generic loop iterator

        if (t_step == t_step_start) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts2(2)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (t_step == t_step_start + 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts2(1)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (t_step == t_step_start + 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts1(2)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (t_step == t_step_start + 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts1(1)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else  ! All other timesteps
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts2(2)%vf(i)%sf(j, k, l) = q_prim_ts2(1)%vf(i)%sf(j, k, l)
                            q_prim_ts2(1)%vf(i)%sf(j, k, l) = q_prim_ts1(2)%vf(i)%sf(j, k, l)
                            q_prim_ts1(2)%vf(i)%sf(j, k, l) = q_prim_ts1(1)%vf(i)%sf(j, k, l)
                            q_prim_ts1(1)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_time_step_cycling

    !> Explicit reaction source (reaction_substeps = 0): warn once, with step and time, when dt*(max species loss rate), the
    !! stiffness estimate of the reaction Jacobian, exceeds the real-axis stability limit of the RK scheme.
    impure subroutine s_report_chem_stiffness(t_step)

        integer, intent(in) :: t_step
        real(wp)            :: x_loc, x_glb

        if (.not. f_chem_stiffness_tracked()) return
        x_loc = dt*chem_loss_max
        chem_loss_max = 0._wp
        chem_stiff_max = max(chem_stiff_max, x_loc)
        if (chem_stiff_warned) return
        call s_mpi_allreduce_max(x_loc, x_glb)
        if (x_glb <= f_chem_stiffness_limit()) return
        chem_stiff_warned = .true.
        if (proc_rank == 0) print '(A,I0,A,ES10.3,A,ES10.3,A,F4.2,A)', &
            & 'WARNING: explicit chemistry is outside its stability limit at step ', t_step, ' (t = ', mytime, &
            & ' s): dt*|lambda_max| >= ', x_glb, ' vs ', f_chem_stiffness_limit(), &
            & '. Set chem_params%reaction_substeps > 0 or reduce dt.'

    end subroutine s_report_chem_stiffness

    !> End-of-run line with the largest dt*|lambda_max| of the explicit reaction source seen on any rank
    impure subroutine s_report_chem_stiffness_summary()

        real(wp) :: x_glb

        if (.not. f_chem_stiffness_tracked()) return
        call s_mpi_allreduce_max(chem_stiff_max, x_glb)
        if (proc_rank == 0) print '(A,ES10.3,A,F4.2,A)', 'Explicit chemistry: max dt*|lambda_max| = ', x_glb, &
            & ' (stability limit ', f_chem_stiffness_limit(), ')'

    end subroutine s_report_chem_stiffness_summary

    !> True when run_time_info is on and the reaction source is advanced explicitly inside the RK stages
    logical function f_chem_stiffness_tracked()

        f_chem_stiffness_tracked = run_time_info .and. chemistry .and. chem_params%reactions &
            & .and. chem_params%reaction_substeps == 0

    end function f_chem_stiffness_tracked

    !> Real-axis stability limit of the SSP-RK scheme (RK3 2.51, RK1 and RK2 2)
    real(wp) function f_chem_stiffness_limit()
        f_chem_stiffness_limit = merge(2.51_wp, 2._wp, time_stepper == time_stepper_rk3)

    end function f_chem_stiffness_limit

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_time_steppers_module

#ifdef FRONTIER_UNIFIED
        use hipfort
        use hipfort_hipmalloc
        use hipfort_check
#endif
        integer :: i, j  !< Generic loop iterators
        ! Deallocating the cell-average conservative variables
#if defined(__NVCOMPILER_GPU_UNIFIED_MEM)
        do j = 1, sys_size
            @:DEALLOCATE(q_cons_ts(1)%vf(j)%sf)
            if (num_ts == 2) then
                if (nv_uvm_out_of_core) then
                    nullify (q_cons_ts(2)%vf(j)%sf)
                else
                    @:DEALLOCATE(q_cons_ts(2)%vf(j)%sf)
                end if
            end if
        end do
        if (num_ts == 2 .and. nv_uvm_out_of_core) then
            deallocate (q_cons_ts_pool_host)
        end if
#elif defined(FRONTIER_UNIFIED)
        do i = 1, num_ts
            do j = 1, sys_size
                nullify (q_cons_ts(i)%vf(j)%sf)
            end do
        end do
#ifdef MFC_MIXED_PRECISION
        call hipCheck(hipHostFree_(c_loc(q_cons_ts_pool_host)))
        nullify (q_cons_ts_pool_host)
        call hipCheck(hipFree_(c_loc(q_cons_ts_pool_device)))
        nullify (q_cons_ts_pool_device)
#else
        call hipCheck(hipHostFree(q_cons_ts_pool_host))
        call hipCheck(hipFree(q_cons_ts_pool_device))
#endif
#else
        do i = 1, num_ts
            do j = 1, sys_size
                @:DEALLOCATE(q_cons_ts(i)%vf(j)%sf)
            end do
        end do
#endif
        do i = 1, num_ts
            @:DEALLOCATE(q_cons_ts(i)%vf)
        end do

        @:DEALLOCATE(q_cons_ts)

        ! Deallocating the cell-average primitive ts variables
        if (probe_wrt) then
            do i = 1, num_probe_ts
                do j = 1, sys_size
                    @:DEALLOCATE(q_prim_ts1(i)%vf(j)%sf,q_prim_ts2(i)%vf(j)%sf )
                end do
                @:DEALLOCATE(q_prim_ts1(i)%vf, q_prim_ts2(i)%vf)
            end do
            @:DEALLOCATE(q_prim_ts1, q_prim_ts2)
        end if

        if (.not. igr) then
            ! Deallocating the cell-average primitive variables
            do i = 1, eqn_idx%adv%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do

            if (mhd) then
                do i = eqn_idx%B%beg, eqn_idx%B%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (hypoelasticity) then
                do i = eqn_idx%stress%beg, eqn_idx%stress%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (cont_damage) then
                @:DEALLOCATE(q_prim_vf(eqn_idx%damage)%sf)
            end if

            if (hyper_cleaning) then
                @:DEALLOCATE(q_prim_vf(eqn_idx%psi)%sf)
            end if

            if (bubbles_euler) then
                do i = eqn_idx%bub%beg, eqn_idx%bub%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (model_eqns == model_eqns_6eq) then
                do i = eqn_idx%int_en%beg, eqn_idx%int_en%end
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

        call s_report_chem_stiffness_summary()

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

        if (chemistry .or. heat_conduction) then
            @:DEALLOCATE(q_T_sf%sf)
        end if
        @:DEALLOCATE(pb_ts(1)%sf)
        @:DEALLOCATE(pb_ts(2)%sf)
        @:DEALLOCATE(rhs_pb)
        @:DEALLOCATE(pb_ts)
        @:DEALLOCATE(mv_ts(1)%sf)
        @:DEALLOCATE(mv_ts(2)%sf)
        @:DEALLOCATE(rhs_mv)
        @:DEALLOCATE(mv_ts)
        do i = 1, num_dims
            @:DEALLOCATE(bc_type(i,1)%sf)
            @:DEALLOCATE(bc_type(i,2)%sf)
        end do
        @:DEALLOCATE(bc_type)
        if (any(time_stepper == (/1, 2, 3/))) then
            @:DEALLOCATE(rk_coef)
            deallocate (rk_c)
        end if

    end subroutine s_finalize_time_steppers_module

end module m_time_steppers
