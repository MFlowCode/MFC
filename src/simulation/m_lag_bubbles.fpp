!>
!! @file m_lag_bubbles.fpp
!! @brief Contains module m_lag_bubbles

#:include 'macros.fpp'

!> @brief This module is used to add the lagrangian subgrid bubble model
module m_lag_bubbles

    ! Dependencies =============================================================

    use m_global_parameters     !< Definitions of the global parameters

    use m_derived_types         !< Definitions of the derived types

    use m_rhs                   !< Right-hand-side (RHS) evaluation procedures

    use m_data_output           !< Run-time info & solution data output procedures

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_kernel_functions      !< Definitions of the kernel functions

    use m_variables_conversion

    use m_compile_specific

    use m_boundary_conditions

    ! ==========================================================================

    implicit none

#ifdef CRAY_ACC_WAR
    !> @name Lagrangian bubble variables
    !> @{

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :), lag_RKcoef)
    !$acc declare link(lag_RKcoef)

    @:CRAY_DECLARE_GLOBAL(integer, dimension(:, :), lag_bub_id)
    !$acc declare link(lag_bub_id)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), lag_bub_R0)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), lag_bub_Rmax)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), lag_bub_Rmin)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), lag_bub_gas_mg)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), lag_bub_gas_betaT)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), lag_bub_gas_betaC)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), lag_bub_dphidt)
    !$acc declare link(lag_bub_R0, lag_bub_Rmax, lag_bub_Rmin, lag_bub_gas_mg, lag_bub_gas_betaT, lag_bub_gas_betaC, lag_bub_dphidt)
    @:CRAY_DECLARE_GLOBAL(logical, dimension(:) :: lag_bub_equilibrium)
    !$acc declare link(lag_bub_equilibrium)

    !(nBub, 1 -> actual val or 2 -> temp val)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:), lag_bub_gas_p)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:), lag_bub_gas_mv)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:), lag_bub_interface_rad)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:), lag_bub_interface_vel)
    !$acc declare link(lag_bub_gas_p, lag_bub_gas_mv, lag_bub_interface_rad, lag_bub_interface_vel)

    !(nBub, 1-> x or 2->y or 3 ->z, 1 -> actual or 2 -> temporal val)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:,:), lag_bub_motion_pos)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:,:), lag_bub_motion_posPrev)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:,:), lag_bub_motion_vel)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:,:), lag_bub_motion_s)
    !$acc declare link(lag_bub_motion_pos, lag_bub_motion_posPrev, lag_bub_motion_vel, lag_bub_motion_s)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:), lag_bub_interface_draddt)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:), lag_bub_interface_dveldt)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:), lag_bub_gas_dpdt)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:), lag_bub_gas_dmvdt)
    !$acc declare link(lag_bub_interface_draddt, lag_bub_interface_dveldt, lag_bub_gas_dpdt, lag_bub_gas_dmvdt)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:,:), lag_bub_motion_dposdt)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:,:,:), lag_bub_motion_dveldt)
    !$acc declare link(lag_bub_motion_dposdt, lag_bub_motion_dveldt)
#else
    real(kind(0d0)), allocatable, dimension(:, :) :: lag_RKcoef    !< RK 4th-5th time stepper coefficients
    !$acc declare create(lag_RKcoef)

    integer, allocatable, dimension(:, :) :: lag_bub_id  !< Global and local IDs
    !$acc declare create(lag_bub_id)

    !(nBub)
    real(kind(0d0)), allocatable, dimension(:) :: lag_bub_R0    !< Initial bubble radius
    real(kind(0d0)), allocatable, dimension(:) :: lag_bub_Rmax  !< Maximum radius
    real(kind(0d0)), allocatable, dimension(:) :: lag_bub_Rmin  !< Minimum radius
    real(kind(0d0)), allocatable, dimension(:) :: lag_bub_gas_mg        !< bubble's gas mass
    real(kind(0d0)), allocatable, dimension(:) :: lag_bub_gas_betaT     !< heatflux model (Preston et al., 2007)
    real(kind(0d0)), allocatable, dimension(:) :: lag_bub_gas_betaC     !< massflux model (Preston et al., 2007)
    real(kind(0d0)), allocatable, dimension(:) :: lag_bub_dphidt        !< subgrid velocity potential (Maeda & Colonius, 2018)
    !$acc declare create(lag_bub_R0, lag_bub_Rmax, lag_bub_Rmin, lag_bub_gas_mg, lag_bub_gas_betaT, lag_bub_gas_betaC, lag_bub_dphidt)

    logical, allocatable, dimension(:) :: lag_bub_equilibrium   !< Equilibrium flag
    !$acc declare create(lag_bub_equilibrium)

    !(nBub, 1 -> actual val or 2 -> temp val)
    real(kind(0d0)), allocatable, dimension(:, :) :: lag_bub_gas_p       !< Pressure in the bubble
    real(kind(0d0)), allocatable, dimension(:, :) :: lag_bub_gas_mv      !< Vapor mass in the bubble
    real(kind(0d0)), allocatable, dimension(:, :) :: lag_bub_interface_rad   !< Bubble radius
    real(kind(0d0)), allocatable, dimension(:, :) :: lag_bub_interface_vel   !< Velocity of the bubble interface
    !$acc declare create(lag_bub_gas_p, lag_bub_gas_mv, lag_bub_interface_rad, lag_bub_interface_vel)

    !(nBub, 1-> x or 2->y or 3 ->z, 1 -> actual or 2 -> temporal val)
    real(kind(0d0)), allocatable, dimension(:, :, :) :: lag_bub_motion_pos        !< Bubble's position
    real(kind(0d0)), allocatable, dimension(:, :, :) :: lag_bub_motion_posPrev    !< Bubble's previous position
    real(kind(0d0)), allocatable, dimension(:, :, :) :: lag_bub_motion_vel        !< Bubble's velocity
    real(kind(0d0)), allocatable, dimension(:, :, :) :: lag_bub_motion_s          !< Bubble's computational cell position in real format
    !$acc declare create(lag_bub_motion_pos, lag_bub_motion_posPrev, lag_bub_motion_vel, lag_bub_motion_s)

    real(kind(0d0)), allocatable, dimension(:, :) :: lag_bub_interface_draddt    !< Time derivative of bubble's radius
    real(kind(0d0)), allocatable, dimension(:, :) :: lag_bub_interface_dveldt    !< Time derivative of bubble's interface velocity
    real(kind(0d0)), allocatable, dimension(:, :) :: lag_bub_gas_dpdt            !< Time derivative of gas pressure
    real(kind(0d0)), allocatable, dimension(:, :) :: lag_bub_gas_dmvdt           !< Time derivative of the vapor mass in the bubble
    !$acc declare create(lag_bub_interface_draddt, lag_bub_interface_dveldt, lag_bub_gas_dpdt, lag_bub_gas_dmvdt)
    real(kind(0d0)), allocatable, dimension(:, :, :) :: lag_bub_motion_dposdt     !< Time derivative of the bubble's position
    real(kind(0d0)), allocatable, dimension(:, :, :) :: lag_bub_motion_dveldt     !< Time derivative of the bubble's velocity
    !$acc declare create(lag_bub_motion_dposdt, lag_bub_motion_dveldt)
#endif

    type(vector_field) :: q_particle    !< Projection of the lagrangian particles in the Eulerian framework
    !$acc declare create(q_particle)

    integer :: q_particle_idx           !< Size of the q_particle vector field
    integer :: nBubs                    !< Number of bubbles in the local domain
    !$acc declare create(nBubs)

    real(kind(0.d0)) :: Rmax_glb, Rmin_glb                  !< Maximum and minimum bubbe size in the local domain
    real(kind(0.d0)) :: lag_voidmax, lag_voidavg, lag_vol  !< Global void fraction indicators
    !$acc declare create(Rmax_glb, Rmin_glb, lag_voidmax, lag_voidavg, lag_vol)

contains

    !> Initializes the lagrangian subgrid bubble solver
        !! @param q_cons_vf Initial conservative variables
    subroutine s_initialize_lag_bubbles(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        ! RKCK 4th/5th time stepper coefficients (Cash J. and Karp A., 1990)
        real(kind(0.d0)), dimension(6) :: &
            RKcoef1 = (/0.2d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/), &
            RKcoef2 = (/3.0d0/40.0d0, 9.0d0/40.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/), &
            RKcoef3 = (/0.3d0, -0.9d0, 1.2d0, 0.0d0, 0.0d0, 0.0d0/), &
            RKcoef4 = (/-11.0d0/54.0d0, 2.5d0, -70.0d0/27.0d0, 35.d0/27.d0, 0.0d0, 0.0d0/), &
            RKcoef5 = (/1631.0d0/55296.0d0, 175.0d0/512.0d0, 575.d0/13824.d0, 44275.d0/110592.d0, 253.d0/4096.d0, 0.0d0/), &
            RKcoef6 = (/37.d0/378.d0, 0.0d0, 250.d0/621.d0, 125.0d0/594.0d0, 0.0d0, 512.0d0/1771.0d0/), &
            RKcoefE = (/37.d0/378.d0 - 2825.0d0/27648.0d0, 0.0d0, 250.d0/621.d0 - 18575.0d0/48384.0d0, &
                        125.0d0/594.0d0 - 13525.0d0/55296.0d0, -277.0d0/14336.0d0, 512.0d0/1771.0d0 - 0.25d0/)

        integer :: i

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================
        !$acc update device(ix, iy, iz)

        ! Allocate space for the Eulerian fields needed to map the effect of the bubbles
        ! comp 1: one minus the voidfraction (1-beta)
        ! comp 2: Temporal derivative of the void fraction, dbetadt
        ! comp 3 - q_particle_idx : Extra-allocated variables in those cases where source terms are required

        if (lag_solver_approach == 1) then ! One-way coupling
            q_particle_idx = 3
        elseif (lag_solver_approach == 2) then ! Two-way coupling
            q_particle_idx = 4
            if (lag_cluster_type >= 4) q_particle_idx = 10 !Subgrid noise model for 2D approximation
        else
            call s_mpi_abort('Please check the lag_solver_approach input')
        end if

        @:ALLOCATE(q_particle%vf(1:q_particle_idx))

        do i = 1, q_particle_idx
            @:ALLOCATE(q_particle%vf(i)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end do

        @:ACC_SETUP_VFs(q_particle)

        ! Allocating space for lagrangian variables
        @:ALLOCATE_GLOBAL(lag_bub_id(1:lag_nBubs_glb, 1:2))
        @:ALLOCATE_GLOBAL(lag_bub_R0(1:lag_nBubs_glb))
        @:ALLOCATE_GLOBAL(lag_bub_Rmax(1:lag_nBubs_glb))
        @:ALLOCATE_GLOBAL(lag_bub_Rmin(1:lag_nBubs_glb))
        @:ALLOCATE_GLOBAL(lag_bub_equilibrium(1:lag_nBubs_glb))
        @:ALLOCATE_GLOBAL(lag_bub_gas_mg(1:lag_nBubs_glb))
        @:ALLOCATE_GLOBAL(lag_bub_gas_betaT(1:lag_nBubs_glb))
        @:ALLOCATE_GLOBAL(lag_bub_gas_betaC(1:lag_nBubs_glb))
        @:ALLOCATE_GLOBAL(lag_bub_dphidt(1:lag_nBubs_glb))
        @:ALLOCATE_GLOBAL(lag_bub_gas_p(1:lag_nBubs_glb, 1:2))
        @:ALLOCATE_GLOBAL(lag_bub_gas_mv(1:lag_nBubs_glb, 1:2))
        @:ALLOCATE_GLOBAL(lag_bub_interface_rad(1:lag_nBubs_glb, 1:2))
        @:ALLOCATE_GLOBAL(lag_bub_interface_vel(1:lag_nBubs_glb, 1:2))
        @:ALLOCATE_GLOBAL(lag_bub_motion_pos(1:lag_nBubs_glb, 1:3, 1:2))
        @:ALLOCATE_GLOBAL(lag_bub_motion_posPrev(1:lag_nBubs_glb, 1:3, 1:2))
        @:ALLOCATE_GLOBAL(lag_bub_motion_vel(1:lag_nBubs_glb, 1:3, 1:2))
        @:ALLOCATE_GLOBAL(lag_bub_motion_s(1:lag_nBubs_glb, 1:3, 1:2))
        @:ALLOCATE_GLOBAL(lag_bub_interface_draddt(1:lag_nBubs_glb, 1:6))
        @:ALLOCATE_GLOBAL(lag_bub_interface_dveldt(1:lag_nBubs_glb, 1:6))
        @:ALLOCATE_GLOBAL(lag_bub_gas_dpdt(1:lag_nBubs_glb, 1:6))
        @:ALLOCATE_GLOBAL(lag_bub_gas_dmvdt(1:lag_nBubs_glb, 1:6))
        @:ALLOCATE_GLOBAL(lag_bub_motion_dposdt(1:lag_nBubs_glb, 1:3, 1:6))
        @:ALLOCATE_GLOBAL(lag_bub_motion_dveldt(1:lag_nBubs_glb, 1:3, 1:6))

        !< Allocate space for the RKCK 4th/5th time stepper coefficients
        @:ALLOCATE_GLOBAL(lag_RKcoef(1:7, 1:6))
        do i = 1, 6
            lag_RKcoef(1, i) = RKcoef1(i)
            lag_RKcoef(2, i) = RKcoef2(i)
            lag_RKcoef(3, i) = RKcoef3(i)
            lag_RKcoef(4, i) = RKcoef4(i)
            lag_RKcoef(5, i) = RKcoef5(i)
            lag_RKcoef(6, i) = RKcoef6(i)
            lag_RKcoef(7, i) = RKcoefE(i)
        end do

        ! Reading bubbles
        call s_read_input_lag_bubbles(q_cons_vf)

    end subroutine s_initialize_lag_bubbles

    !> The purpose of this procedure is to read the input file with the bubbles' information
        !! @param q_cons_vf Conservative variables
    subroutine s_read_input_lag_bubbles(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        real(kind(0.d0)), dimension(8) :: inputparticle
        real(kind(0.d0)) :: qtime
        integer :: id, nparticles, save_count
        integer :: i, j, k, l
        logical :: file_exist, indomain

        character(LEN=path_len + 2*name_len) :: path_D_dir !<

        dt_max = dt !Initial and largest dt

        ! Initialize number of particles
        nparticles = 0
        id = 0

        ! Read the input lag_bubble file or restart point
        if (cfl_dt) then
            save_count = n_start
            qtime = n_start*t_save
        else
            save_count = t_step_start
            qtime = t_step_start*dt
        end if

        if (save_count == 0) then
            if (proc_rank == 0) print *, 'Reading lag_bubbles input file at save_count: ', save_count
            inquire (file='input/lag_bubbles.dat', exist=file_exist)
            if (file_exist) then
                open (unit=85, file='input/lag_bubbles.dat', form='formatted')
101             read (85, *, end=102) (inputparticle(i), i=1, 8)
                indomain = particle_in_domain(inputparticle(1:3))
                id = id + 1
                if (id > lag_nBubs_glb .and. proc_rank == 0) call s_mpi_abort('Actual number of bubbles is larger than lag_nBubs_glb')
                if (indomain) then
                    nparticles = nparticles + 1
                    call s_add_lag_bubble(inputparticle, q_cons_vf, nparticles)
                    lag_bub_id(nparticles, 1) = id  !global ID
                    lag_bub_id(nparticles, 2) = nparticles  !local ID
                    nBubs = nparticles ! local number of bubbles
                end if
                goto 101
102             continue
            else
                stop "if you include lag_bubbles, you have to initialize them in input/lag_bubbles.dat"
            end if
        else
            if (proc_rank == 0) print *, 'Restarting lag_bubbles at save_count: ', save_count
            call s_add_lag_bubble_restart(nparticles, save_count)
        end if

        print *, " LAG BUBBLES RUNNING, in proc", proc_rank, "number:", nparticles, "/", id

        Rmax_glb = min(dflt_real, -dflt_real)
        Rmin_glb = max(dflt_real, -dflt_real)

        !$acc update device(lag_bubbles, lag_solver_approach, lag_bubble_model, lag_cluster_type, lag_pressure_corrector, &
        !$acc lag_adap_dt, lag_smooth_type, lag_heatTransfer_model, lag_massTransfer_model, lag_write_bubbles, lag_write_bubble_stats, &
        !$acc lag_nBubs_glb, lag_epsilonb, lag_rkck_tolerance, lag_charwidth, lag_valmaxvoid, csonhost, vischost, Thost, gammagas, &
        !$acc gammavapor, pvap, cpgas, cpvapor, kgas, kvapor, Rgas, Rvap, diffcoefvap, sigmabubble)

        !$acc update device(lag_bub_id, lag_bub_R0, lag_bub_Rmax, lag_bub_Rmin, lag_bub_gas_mg, lag_bub_gas_betaT, lag_bub_gas_betaC, &
        !$acc lag_bub_dphidt, lag_bub_equilibrium, lag_bub_gas_p, lag_bub_gas_mv, lag_bub_interface_rad, lag_bub_interface_vel, &
        !$acc lag_bub_motion_pos, lag_bub_motion_posPrev, lag_bub_motion_vel, lag_bub_motion_s, lag_bub_interface_draddt, &
        !$acc lag_bub_interface_dveldt, lag_bub_gas_dpdt, lag_bub_gas_dmvdt, lag_bub_motion_dposdt, lag_bub_motion_dveldt, Rmax_glb, &
        !$acc Rmin_glb, nBubs, lag_RKcoef)

        !$acc update device(dx, dy, dz, x_cb, x_cc, y_cb, y_cc, z_cb, z_cc)

        !Populate temporal variables
        call s_transfer_data_to_tmp()

        if (lag_write_bubbles) call s_write_lag_particles(qtime)

        call s_smear_voidfraction()

        if (save_count == 0) then
            ! Create ./D directory
            write (path_D_dir, '(A,I0,A,I0)') trim(case_dir)//'/D'
            call my_inquire(path_D_dir, file_exist)
            if (.not. file_exist) call s_create_directory(trim(path_D_dir))

            call s_write_restart_lag_bubbles(save_count)
            call s_write_void_evol(mytime)
        end if

    end subroutine s_read_input_lag_bubbles

    !> The purpose of this procedure is to add information of the bubbles when starting fresh
        !! @param inputparticle Particle information
        !! @param q_cons_vf Conservative variables
        !! @param nparticles Local ID of the particle
    subroutine s_add_lag_bubble(inputparticle, q_cons_vf, nparticles)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(kind(0.d0)), dimension(8), intent(in) :: inputparticle
        integer, intent(in) :: nparticles
        integer :: i

        real(kind(0.d0)) :: pliq, volparticle, concvap, totalmass, kparticle, cpparticle
        real(kind(0.d0)) :: omegaN, PeG, PeT, rhol, cson, pcrit, qv, gamma, pi_inf, dynP
        integer, dimension(3) :: cell
        real(kind(0.d0)), dimension(2) :: Re
        real(kind(0.d0)) :: massflag, heatflag

        massflag = 0.0d0
        heatflag = 0.0d0
        if (lag_massTransfer_model) massflag = 1.0d0
        if (lag_heatTransfer_model) heatflag = 1.0d0

        lag_bub_R0(nparticles) = inputparticle(7)
        lag_bub_Rmax(nparticles) = min(dflt_real, -dflt_real)
        lag_bub_Rmin(nparticles) = max(dflt_real, -dflt_real)
        lag_bub_equilibrium(nparticles) = .false.
        lag_bub_dphidt(nparticles) = 0.0d0
        lag_bub_interface_rad(nparticles, 1) = inputparticle(7)
        lag_bub_interface_vel(nparticles, 1) = inputparticle(8)
        lag_bub_motion_pos(nparticles, 1:3, 1) = inputparticle(1:3)
        lag_bub_motion_posPrev(nparticles, 1:3, 1) = lag_bub_motion_pos(nparticles, 1:3, 1)
        lag_bub_motion_vel(nparticles, 1:3, 1) = inputparticle(4:6)

        if (cyl_coord .and. p == 0) then
            lag_bub_motion_pos(nparticles, 2, 1) = dsqrt(lag_bub_motion_pos(nparticles, 2, 1)**2d0 + &
                                                         lag_bub_motion_pos(nparticles, 3, 1)**2d0)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            lag_bub_motion_pos(nparticles, 3, 1) = atan2(inputparticle(3), inputparticle(2))
            lag_bub_motion_posPrev(nparticles, 1:3, 1) = lag_bub_motion_pos(nparticles, 1:3, 1)
        end if

        cell = -buff_size
        call s_locate_cell(lag_bub_motion_pos(nparticles, 1:3, 1), cell, lag_bub_motion_s(nparticles, 1:3, 1))

        ! If particle is in the ghost cells, find the closest non-ghost cell
        cell(1) = min(max(cell(1), 0), m)
        cell(2) = min(max(cell(2), 0), n)
        if (p > 0) cell(3) = min(max(cell(3), 0), p)

        call s_convert_to_mixture_variables(q_cons_vf, cell(1), cell(2), cell(3), rhol, gamma, pi_inf, qv, Re)

        dynP = 0.0d0
        do i = 1, num_dims
            dynP = dynP + 0.5*q_cons_vf(num_fluids + i)%sf(cell(1), cell(2), cell(3))**2/rhol
        end do
        pliq = (q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)) - dynP - pi_inf)/gamma

        if (pliq < 0) print *, "Negative pressure", proc_rank, &
            q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)), pi_inf, gamma, pliq, cell, dynP

        ! Initial particle pressure
        lag_bub_gas_p(nparticles, 1) = pliq + 2.0*sigmabubble/lag_bub_R0(nparticles)
        if (sigmabubble /= 0.0d0) then
            pcrit = pvap - 4.0d0*sigmabubble/(3.d0*sqrt(3.0d0*lag_bub_gas_p(nparticles, 1)*lag_bub_R0(nparticles)**3/(2.0d0*sigmabubble)))
            pref = lag_bub_gas_p(nparticles, 1)
        else
            pcrit = 0.0d0
        end if

        ! Initial particle mass
        volparticle = 4.0d0/3.0d0*pi*lag_bub_R0(nparticles)**3 ! volume
        lag_bub_gas_mv(nparticles, 1) = pvap*volparticle*(1.0d0/(Rvap*Thost))*(massflag) ! vapermass
        lag_bub_gas_mg(nparticles) = (lag_bub_gas_p(nparticles, 1) - pvap*(massflag))*volparticle*(1.0d0/(Rgas*Thost)) ! gasmass
        if (lag_bub_gas_mg(nparticles) <= 0.0d0) stop 'the initial mass of gas inside the bubble is negative. Check your initial conditions'
        totalmass = lag_bub_gas_mg(nparticles) + lag_bub_gas_mv(nparticles, 1) ! totalmass

        ! Bubble natural frequency
        concvap = lag_bub_gas_mv(nparticles, 1)/(lag_bub_gas_mv(nparticles, 1) + lag_bub_gas_mg(nparticles))
        omegaN = (3.0d0*(lag_bub_gas_p(nparticles, 1) - pvap*(massflag)) + 4.0d0*sigmabubble/lag_bub_R0(nparticles))/rhol
        if (pvap*(massflag) > lag_bub_gas_p(nparticles, 1)) then
            print *, 'Not allowed: bubble initially located in a region with pressure below the vapor pressure'
            print *, 'location:', lag_bub_motion_pos(nparticles, 1:3, 1)
            stop
        end if
        omegaN = dsqrt(omegaN/lag_bub_R0(nparticles)**2)

        cpparticle = concvap*cpvapor + (1.0d0 - concvap)*cpgas
        kparticle = concvap*kvapor + (1.0d0 - concvap)*kgas

        ! Mass and heat transfer coefficients (based on Preston 2007)
        PeT = totalmass/volparticle*cpparticle*lag_bub_R0(nparticles)**2*omegaN/kparticle
        lag_bub_gas_betaT(nparticles) = f_transfercoeff(PeT, 1.0d0)*(heatflag)
        PeG = lag_bub_R0(nparticles)**2*omegaN/diffcoefvap
        lag_bub_gas_betaC(nparticles) = f_transfercoeff(PeG, 1.0d0)*(massflag)

        ! Terms to work out directly the heat flux in getfluxes
        lag_bub_gas_betaT(nparticles) = lag_bub_gas_betaT(nparticles)*kparticle

        if (lag_bub_gas_mg(nparticles) <= 0.0d0) stop "Error encontered: Negative gas mass in the bubble, &
&                                                                            check if the bubble is inside the domain."

    end subroutine s_add_lag_bubble

    !> The purpose of this procedure is to add information of the bubbles from a restart point.
        !! @param nparticles Local ID of the particle
        !! @param save_count File identifier
    subroutine s_add_lag_bubble_restart(nparticles, save_count)

        integer :: nparticles, save_count

        character(LEN=path_len + 2*name_len) :: file_loc

#ifdef MFC_MPI
        real(kind(0.d0)), dimension(20) :: inputvals
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: disp
        integer :: view

        integer, dimension(3) :: cell
        logical :: indomain, particle_file, file_exist

        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer :: ifile, ireq, ierr, data_size, tot_data, id
        integer :: i

        write (file_loc, '(a,i0,a)') 'lag_bubbles_mpi_io_', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (file=trim(file_loc), exist=file_exist)

        if (file_exist) then
            if (proc_rank == 0) then
                open (9, file=trim(file_loc), form='unformatted', status='unknown')
                read (9) tot_data, mytime, dt
                close (9)
                print *, 'Reading lag_bubbles_mpi_io: ', tot_data, mytime, dt
            end if
        else
            print '(a)', trim(file_loc)//' is missing. exiting ...'
            call s_mpi_abort
        end if

        call MPI_BCAST(tot_data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mytime, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        gsizes(1) = tot_data
        gsizes(2) = 21
        lsizes(1) = tot_data
        lsizes(2) = 21
        start_idx_part(1) = 0
        start_idx_part(2) = 0

        call MPI_type_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                      MPI_ORDER_FORTRAN, MPI_doUBLE_PRECISION, view, ierr)
        call MPI_type_COMMIT(view, ierr)

        ! Open the file to write all flow variables
        write (file_loc, '(a,i0,a)') 'lag_bubbles_', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (file=trim(file_loc), exist=particle_file)

        if (particle_file) then
            call MPI_FILE_open(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)
            disp = 0d0
            call MPI_FILE_SET_VIEW(ifile, disp, MPI_doUBLE_PRECISION, view, &
                                   'native', mpi_info_null, ierr)
            allocate (MPI_IO_DATA_lag_bubbles(tot_data, 1:21))
            call MPI_FILE_read_ALL(ifile, MPI_IO_DATA_lag_bubbles, 21*tot_data, &
                                   MPI_doUBLE_PRECISION, status, ierr)
            do i = 1, tot_data
                id = int(MPI_IO_DATA_lag_bubbles(i, 1))
                inputvals(1:20) = MPI_IO_DATA_lag_bubbles(i, 2:21)
                indomain = particle_in_domain(inputvals(1:3))
                if (indomain .and. (id > 0)) then
                    nparticles = nparticles + 1
                    nBubs = nparticles ! local number of bubbles
                    lag_bub_id(nparticles, 1) = id  !global ID
                    lag_bub_id(nparticles, 2) = nparticles  !local ID
                    lag_bub_motion_pos(nparticles, 1:3, 1) = inputvals(1:3)
                    lag_bub_motion_posPrev(nparticles, 1:3, 1) = inputvals(4:6)
                    lag_bub_motion_vel(nparticles, 1:3, 1) = inputvals(7:9)
                    lag_bub_interface_rad(nparticles, 1) = inputvals(10)
                    lag_bub_interface_vel(nparticles, 1) = inputvals(11)
                    lag_bub_R0(nparticles) = inputvals(12)
                    lag_bub_Rmax(nparticles) = inputvals(13)
                    lag_bub_Rmin(nparticles) = inputvals(14)
                    lag_bub_dphidt(nparticles) = inputvals(15)
                    lag_bub_gas_p(nparticles, 1) = inputvals(16)
                    lag_bub_gas_mv(nparticles, 1) = inputvals(17)
                    lag_bub_gas_mg(nparticles) = inputvals(18)
                    lag_bub_gas_betaT(nparticles) = inputvals(19)
                    lag_bub_gas_betaC(nparticles) = inputvals(20)
                    lag_bub_equilibrium(nparticles) = .false.
                    cell = -buff_size
                    call s_locate_cell(lag_bub_motion_pos(nparticles, 1:3, 1), cell, lag_bub_motion_s(nparticles, 1:3, 1))
                end if
            end do
            deallocate (MPI_IO_DATA_lag_bubbles)
        end if
        call MPI_FILE_CLOSE(ifile, ierr)
#endif

    end subroutine s_add_lag_bubble_restart

    !>  Contains the two-way and one-way Euler-Lagrange coupled algorithm, including the bubble dynamics subroutines.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
        !! @param rhs_vf Calculated change of conservative variables
        !! @param t_step Current global time step
        !! @param qtime Current time from the RKCK stepper
        !! @param step Current step from the RKCK stepper
    subroutine s_compute_el_coupled_solver(q_cons_vf, q_prim_vf, rhs_vf, pb, rhs_pb, mv, rhs_mv, t_step, time_avg, qtime, step)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, rhs_pb
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: mv, rhs_mv
        integer, intent(in) :: step, t_step
        real(kind(0.d0)), intent(inout) :: time_avg
        real(kind(0.d0)), intent(in) :: qtime

        real(kind(0.d0)) :: gammaparticle, vaporflux, heatflux
        integer :: i, j, k, l
        real(kind(0.d0)) :: preterm1, term2, paux, pint, Romega, term1_fac, Rb

        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size
        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        !$acc update device(ix, iy, iz)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        rhs_vf(i)%sf(j, k, l) = 0.0d0
                    end do
                end do
            end do
        end do

        time_tmp = qtime
        !$acc update device (time_tmp)

        !< Compute eulerian background
        call s_compute_rhs(q_cons_vf, q_prim_vf, rhs_vf, pb, rhs_pb, mv, rhs_mv, t_step, time_avg)

        !< Compute lagrangian bubbles
        call s_populate_variables_buffers(q_prim_vf, pb, mv)

        if (lag_pressure_corrector) then ! Subgrid p_inf model from Maeda and Colonius (2018).

            ! Calculate velocity potentials (valid for one bubble per cell)
            !$acc parallel loop gang vector default(present) private(k)
            do k = 1, nBubs
                paux = f_pressure_inf(k, q_prim_vf, 2, preterm1, term2, Romega)
                pint = f_pressureliq_int(k)
                pint = pint + 0.5d0*lag_bub_interface_vel(k, 2)**2

                if (lag_cluster_type == 2) then
                    lag_bub_dphidt(k) = (paux - pint) + term2
                    ! Accounting for the potential induced by the bubble averaged over the control volume
                    ! Note that this is based on the incompressible flow assumption near the bubble.
                    Rb = lag_bub_interface_rad(k, 2)
                    term1_fac = 3.0d0/2.0d0*(Rb*(Romega**2d0 - Rb**2d0))/(Romega**3d0 - Rb**3d0)
                    lag_bub_dphidt(k) = lag_bub_dphidt(k)/(1 - term1_fac)
                end if
            end do
        end if

        ! Gaseous core evolution
        !$acc parallel loop gang vector default(present) private(k) copyin(step)
        do k = 1, nBubs
            if (.not. lag_bub_equilibrium(k)) then
                call s_compute_interface_fluxes(k, vaporflux, heatflux, gammaparticle)
                lag_bub_gas_dpdt(k, step) = -3.0d0*gammaparticle/lag_bub_interface_rad(k, 2)* &
                                            (lag_bub_gas_p(k, 2)*lag_bub_interface_vel(k, 2) - &
                                             heatflux - (Rvap*Thost)*vaporflux)
                lag_bub_gas_dmvdt(k, step) = 4.0d0*pi*lag_bub_interface_rad(k, 2)**2*vaporflux
            else
                lag_bub_gas_dpdt(k, step) = 0.0d0
                lag_bub_gas_dmvdt(k, step) = 0.0d0
            end if
        end do

        ! Radial motion model
        if (lag_bubble_model == 1) then
            !$acc parallel loop gang vector default(present) private(k) copyin(step)
            do k = 1, nBubs
                if (.not. lag_bub_equilibrium(k)) then
                    call s_compute_KM(k, step, q_prim_vf)
                else
                    lag_bub_interface_dveldt(k, step) = 0.0d0
                end if
                lag_bub_interface_draddt(k, step) = lag_bub_interface_vel(k, 2)
            end do
        else
            !$acc parallel loop gang vector default(present) private(k) copyin(step)
            do k = 1, nBubs
                lag_bub_interface_dveldt(k, step) = 0.0d0
                lag_bub_interface_draddt(k, step) = 0.0d0
            end do
        end if

        ! Bubbles remain in a fixed position
        !$acc parallel loop collapse(2) gang vector default(present) private(k) copyin(step)
        do k = 1, nBubs
            do l = 1, 3
                lag_bub_motion_dposdt(k, l, step) = 0.0d0
                lag_bub_motion_dveldt(k, l, step) = 0.0d0
            end do
        end do

        !< Euler-Lagrange coupling
        if (lag_solver_approach == 2) call s_add_sources(q_cons_vf, q_prim_vf, rhs_vf)

    end subroutine s_compute_el_coupled_solver

    !>  The purpose of this subroutine is to smear the effect of the bubbles in the Eulerian framework
    subroutine s_smear_voidfraction()

        integer :: i, j, k, l

        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
        if (n > 0) iy%beg = -buff_size; 
        if (p > 0) iz%beg = -buff_size; 
        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        !$acc update device(ix, iy, iz)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, q_particle_idx
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        q_particle%vf(i)%sf(j, k, l) = 0.0d0
                    end do
                end do
            end do
        end do

        call s_smoothfunction(nBubs, lag_bub_interface_rad, lag_bub_interface_vel, &
                              lag_bub_motion_s, lag_bub_motion_pos, q_particle)

        !Store 1-beta
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end
                    q_particle%vf(1)%sf(j, k, l) = 1.0d0 - q_particle%vf(1)%sf(j, k, l)
                end do
            end do
        end do

        ! Limiting void fraction given max value
        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    q_particle%vf(1)%sf(i, j, k) = max(q_particle%vf(1)%sf(i, j, k), 1.d0 - lag_valmaxvoid)
                end do
            end do
        end do

    end subroutine s_smear_voidfraction

    !>  The purpose of this subroutine is to add the bubbles source terms following the formulation of Maeda and Colonius (2018)
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Conservative variables
        !! @param rhs_vf Time derivative of the conservative variables
    subroutine s_add_sources(q_cons_vf, q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k, l

        if (lag_cluster_type >= 4) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        do l = 1, E_idx
                            if (q_particle%vf(1)%sf(i, j, k) > (1.0d0 - lag_valmaxvoid)) then
                                rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + &
                                                        q_cons_vf(l)%sf(i, j, k)*(q_particle%vf(2)%sf(i, j, k) + &
                                                                                  q_particle%vf(5)%sf(i, j, k))
                            end if
                        end do
                    end do
                end do
            end do
        else
            !$acc parallel loop collapse(4) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        do l = 1, E_idx
                            if (q_particle%vf(1)%sf(i, j, k) > (1.0d0 - lag_valmaxvoid)) then
                                rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + &
                                                        q_cons_vf(l)%sf(i, j, k)/q_particle%vf(1)%sf(i, j, k)* &
                                                        q_particle%vf(2)%sf(i, j, k)
                            end if
                        end do
                    end do
                end do
            end do
        end if

        do l = 1, num_dims

            call s_gradient_dir(q_prim_vf(E_idx), q_particle%vf(3), l)

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        if (q_particle%vf(1)%sf(i, j, k) > (1.0d0 - lag_valmaxvoid)) then
                            rhs_vf(num_fluids + l)%sf(i, j, k) = rhs_vf(num_fluids + l)%sf(i, j, k) - &
                                                                 (1.0d0 - q_particle%vf(1)%sf(i, j, k))/ &
                                                                 q_particle%vf(1)%sf(i, j, k)* &
                                                                 q_particle%vf(3)%sf(i, j, k)
                        end if
                    end do
                end do
            end do

            !source in energy
            !$acc parallel loop collapse(3) gang vector default(present)
            do k = iz%beg, iz%end
                do j = iy%beg, iy%end
                    do i = ix%beg, ix%end
                        q_particle%vf(3)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k)*q_prim_vf(num_fluids + l)%sf(i, j, k)
                    end do
                end do
            end do

            call s_gradient_dir(q_particle%vf(3), q_particle%vf(4), l)

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        if (q_particle%vf(1)%sf(i, j, k) > (1.0d0 - lag_valmaxvoid)) then
                            rhs_vf(E_idx)%sf(i, j, k) = rhs_vf(E_idx)%sf(i, j, k) - &
                                                        q_particle%vf(4)%sf(i, j, k)*(1.0d0 - q_particle%vf(1)%sf(i, j, k))/ &
                                                        q_particle%vf(1)%sf(i, j, k)
                        end if
                    end do
                end do
            end do
        end do

    end subroutine s_add_sources

    !>  This subroutine computes the Keller-Miksis equation
        !! @param nparticles Particle identifier
        !! @param step Current time-stage in the 4th/5th order Runge—Kutta–Cash–Karp time-stepping algorithm
        !! @param q_prim_vf Primitive variables
    subroutine s_compute_KM(nparticles, step, q_prim_vf)
        !$acc routine seq

        integer, intent(in) :: nparticles, step
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(kind(0.d0)) :: pliqint, pbubble, deltaP, pinf, termI, aux1
        real(kind(0.d0)) :: aux2, velint, rhol, cson, E, H, qv, gamma, pi_inf
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho, alpha
        integer, dimension(3) :: cell
        real(kind(0.0d0)), dimension(3) :: scoord

        real(kind(0.d0)), dimension(num_dims) :: vel
        real(kind(0.d0)), dimension(2) :: Re
        integer :: i

        pliqint = f_pressureliq_int(nparticles) ! pressure at the bubble wall

        scoord = lag_bub_motion_s(nparticles, 1:3, 2)
        call s_get_cell(scoord, cell)
        pinf = f_pressure_inf(nparticles, q_prim_vf, 1, aux1, aux2) ! getting p_inf

        do i = 1, num_fluids
            alpha_rho(i) = q_prim_vf(advxb + i - 1)%sf(cell(1), cell(2), cell(3))*q_prim_vf(i)%sf(cell(1), cell(2), cell(3))
            alpha(i) = q_prim_vf(advxb + i - 1)%sf(cell(1), cell(2), cell(3))
        end do

        call s_convert_species_to_mixture_variables_acc(rhol, gamma, pi_inf, qv, alpha, alpha_rho, Re, cell(1), cell(2), cell(3))

        ! Computing speed of sound
        do i = 1, num_dims
            vel(i) = q_prim_vf(i + num_fluids)%sf(cell(1), cell(2), cell(3))
        end do
        E = gamma*pinf + pi_inf + 0.5d0*rhol*dot_product(vel, vel)
        H = (E + pinf)/rhol
        cson = sqrt((H - 0.5d0*dot_product(vel, vel))/gamma)

        deltaP = pliqint - pinf
        termI = 0.0d0
        velint = lag_bub_interface_vel(nparticles, 2) - lag_bub_gas_dmvdt(nparticles, step)/(4.0d0*pi*lag_bub_interface_rad(nparticles, 2)**2*rhol)

        lag_bub_interface_dveldt(nparticles, step) = ((1.0d0 + velint/cson)*deltaP/rhol &
                                                      + termI &
                                                      + lag_bub_gas_dpdt(nparticles, step)*lag_bub_interface_rad(nparticles, 2)/rhol/cson &
                                                      - velint**2*3.0d0/2.0d0*(1.0d0 - velint/3.0d0/cson)) &
                                                     /(lag_bub_interface_rad(nparticles, 2)*(1.0d0 - velint/cson))

    end subroutine s_compute_KM

    !>  This subroutine computes the fluxes at the bubble's interface
        !! @param nparticles Particle identifier
        !! @param vaporflux Mass flux
        !! @param heatflux Heat flux
        !! @param gammabubble Specific heat of the vapor-gas mixture in the bubble
    subroutine s_compute_interface_fluxes(nparticles, vaporflux, heatflux, gammabubble)
        !$acc routine seq
        integer, intent(in) :: nparticles
        real(kind(0.d0)), intent(out) :: vaporflux, heatflux, gammabubble
        real(kind(0.d0)) :: concvapint, bubbleTemp, volbubble, avgconc, Rmixt, rhogas

        concvapint = 0.d0
        if (lag_massTransfer_model) then
            concvapint = (Rvap/Rgas)*(lag_bub_gas_p(nparticles, 2)/pvap - 1.0d0)
            concvapint = 1.0d0/(1.0d0 + concvapint)
        end if

        bubbleTemp = lag_bub_gas_mg(nparticles)*Rgas + lag_bub_gas_mv(nparticles, 2)*Rvap
        volbubble = 4.0d0/3.0d0*pi*lag_bub_interface_rad(nparticles, 2)**3
        bubbleTemp = lag_bub_gas_p(nparticles, 2)*volbubble/bubbleTemp

        gammabubble = concvapint*gammavapor + (1.0d0 - concvapint)*gammagas
        heatflux = -(gammabubble - 1.0d0)/gammabubble*lag_bub_gas_betaT(nparticles)*(bubbleTemp - Thost)/lag_bub_interface_rad(nparticles, 2)

        avgconc = lag_bub_gas_mv(nparticles, 2)/(lag_bub_gas_mg(nparticles) + lag_bub_gas_mv(nparticles, 2))
        Rmixt = concvapint*Rvap + (1.0d0 - concvapint)*Rgas

        concvapint = min(concvapint, 0.99d0)
        vaporflux = (1.0d0 - concvapint)*lag_bub_interface_rad(nparticles, 2)
        rhogas = (lag_bub_gas_mg(nparticles) + lag_bub_gas_mv(nparticles, 2))/(4.0d0/3.0d0*pi*lag_bub_interface_rad(nparticles, 2)**3)
        vaporflux = -diffcoefvap*lag_bub_gas_betaC(nparticles)*(avgconc - concvapint)*rhogas/vaporflux

    end subroutine s_compute_interface_fluxes

    !> This function calculates the heat and mass transfer coefficients following the formulation proposed by Preston et al. (2007)
        !! @param Pe Peclet number of heat or mass transfer
        !! @param omegaN Characteristic frequency
    function f_transfercoeff(Pe, omegaN)

        real(kind(0.d0)) :: f_transfercoeff, Pe, omegaN
        complex :: transferfunc, auxc

        transferfunc = csqrt(cmplx(0.0d0, Pe*omegaN))
        auxc = (cexp(-cmplx(2.0d0, 0.0d0)*transferfunc) + 1.0d0)/(-cexp(-cmplx(2.0d0, 0.0d0)*transferfunc) + 1.0d0)
        transferfunc = transferfunc*auxc - 1.0d0
        transferfunc = 1.0d0/transferfunc - 3.0d0/cmplx(0.0d0, Pe*omegaN)
        f_transfercoeff = real(1.0d0/transferfunc)

    end function f_transfercoeff

    !> This function calculates the pressure at the bubble wall (bubble-liquid interface).
        !! @param nparticles Particle identifier
    function f_pressureliq_int(nparticles)
        !$acc routine seq
        integer, intent(in) :: nparticles
        real(kind(0.d0)) :: f_pressureliq_int

        real(kind(0.d0)) :: radius, bubblevel, pbubble

        pbubble = lag_bub_gas_p(nparticles, 2)
        radius = lag_bub_interface_rad(nparticles, 2)
        bubblevel = lag_bub_interface_vel(nparticles, 2)

        f_pressureliq_int = pbubble - 2.0d0*sigmabubble/radius - 4.0d0*vischost*bubblevel/radius

    end function f_pressureliq_int

    !> The purpose of this procedure is obtain the pressure that drives the bubble oscillations p_inf
        !! @param nparticles Particle identifier
        !! @param q_prim_vf  Primitive variables
        !! @param ptype 1: p at infinity, 2: averaged P at the bubble location
    function f_pressure_inf(nparticles, q_prim_vf, ptype, preterm1, term2, Romega)
        !$acc routine seq
        integer, intent(in) :: nparticles
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: ptype
        real(kind(0.d0)), intent(out), optional :: preterm1, term2, Romega
        real(kind(0.d0)) :: f_pressure_inf

        real(kind(0.d0)), dimension(3) :: scoord
        real(kind(0.d0)), dimension(3) :: distance, center
        real(kind(0.d0)) :: jac, dij, dpotjdt, dc, vol, aux, &
                            volgas, term1, Rbeq, denom, stddsv, &
                            charvol, charpres, charvol2, charpres2
        integer, dimension(3) :: cell, cellaux
        integer, dimension(3) :: epsilonbaux

        integer :: i, j, k, dir
        integer :: cellx, celly, cellz
        logical :: celloutside

        scoord = lag_bub_motion_s(nparticles, 1:3, 2)
        f_pressure_inf = 0.0d0
        call s_get_cell(scoord, cell)

        if ((lag_cluster_type == 1)) then
            !getting p_cell in terms of only the current cell by interpolation

            call s_get_char_vol(cell(1), cell(2), cell(3), vol)

            ! Getting the cell volulme as Omega
            call s_get_char_dist(cell(1), cell(2), cell(3), stddsv)

            !p_cell (interpolated)
            f_pressure_inf = f_interpolate(scoord, q_prim_vf(E_idx))

            !R_Omega
            dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)

        else if (lag_cluster_type >= 2) then
            ! Bubble dynamic closure from Maeda and Colonius (2018)

            ! Range of cells included in Omega
            if (lag_smooth_type == 1) epsilonbaux(:) = 3

            charvol = 0.d0
            charpres = 0.d0
            charvol2 = 0.d0
            charpres2 = 0.d0
            vol = 0.d0
            if (num_dims == 3) then
                k = -epsilonbaux(3)
            else
                k = 0
            end if
            i = -epsilonbaux(1); j = -epsilonbaux(2)

3001        if ((i <= epsilonbaux(1)) .and. (j <= epsilonbaux(2))) then
                celloutside = .false.
                cellaux(1) = cell(1) + i
                cellaux(2) = cell(2) + j
                cellaux(3) = cell(3) + k

                !Check ghost part in x-direction
                if (cellaux(1) < -buff_size) then
                    celloutside = .true.
                    i = i + 1
                end if

                !Check ghost part in y-direction
                if (cellaux(2) < -buff_size) then
                    celloutside = .true.
                    j = j + 1
                end if
                if (cyl_coord .and. (num_dims /= 3)) then
                    if (y_cc(cellaux(2)) < 0d0) then
                        celloutside = .true.
                        j = j + 1
                    end if
                end if

                !Check ghost part in z-direction
                if (num_dims == 3) then
                    if (cellaux(3) < -buff_size) then
                        celloutside = .true.
                        k = k + 1
                    end if
                end if

                if (cellaux(1) > m + buff_size) celloutside = .true.
                if (cellaux(2) > n + buff_size) celloutside = .true.
                if (cellaux(3) > p + buff_size) celloutside = .true.

                if (.not. celloutside) then
                    call s_get_char_vol(cellaux(1), cellaux(2), cellaux(3), vol)
                    charvol = charvol + vol
                    charpres = charpres + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3))*vol
                    charvol2 = charvol2 + vol*q_particle%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                    charpres2 = charpres2 + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                *vol*q_particle%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                end if

                if (j < epsilonbaux(2)) then
                    j = j + 1
                    goto 3001
                end if

3002            j = -epsilonbaux(2)
                i = i + 1
                goto 3001

            end if

3003        if ((num_dims == 3) .and. (k < epsilonbaux(3))) then
                k = k + 1
                i = -epsilonbaux(1); j = -epsilonbaux(2)
                goto 3001
            end if
            f_pressure_inf = charpres2/charvol2
            vol = charvol
            dc = (3.0d0*abs(vol)/(4.0d0*pi))**(1.0d0/3.0d0)
        else

            stop "Check cluterflag. Exiting ..."

        end if

        if (lag_pressure_corrector .and. present(preterm1)) then

            !Valid if only one bubble exists per cell
            volgas = lag_bub_interface_rad(nparticles, 2)**3
            denom = lag_bub_interface_rad(nparticles, 2)**2
            term1 = lag_bub_dphidt(nparticles)*lag_bub_interface_rad(nparticles, 2)**2
            term2 = lag_bub_interface_vel(nparticles, 2)*lag_bub_interface_rad(nparticles, 2)**2

            Rbeq = volgas**(1.0d0/3.0d0) !surrogate bubble radius
            aux = dc**3 - Rbeq**3
            term2 = term2/denom
            term2 = 3.0d0/2.0d0*term2**2*Rbeq**3*(1.0d0 - Rbeq/dc)/aux
            preterm1 = 3.0d0/2.0d0*Rbeq*(dc**2 - Rbeq**2)/(aux*denom)

            !Control volume radius
            if (present(Romega)) Romega = dc

            ! Getting p_inf
            if (ptype == 1) then
                f_pressure_inf = f_pressure_inf + preterm1*term1 + term2
            end if

        end if

    end function f_pressure_inf

    !>  This subroutine updates the Euler-Lagrange temporal variables before entering to the next time-stage in the time stepper.
        !! @param RKstep Current time step in the adaptative stepper
        !! @param q_cons_ts Conservative variables
        !! @param rhs_ts Time derivatives of the conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_update_tmp_rkck(RKstep, q_cons_ts, rhs_ts, q_prim_vf)

        integer, intent(in) :: RKstep
        type(vector_field), dimension(:), intent(inout) :: q_cons_ts
        type(vector_field), dimension(:), intent(inout) :: rhs_ts
        type(scalar_field), dimension(:), intent(inout) :: q_prim_vf

        integer :: i, j, k, l, q, ierr
        logical :: change, indomain
        integer, dimension(3) :: cell
        real(kind(0.d0)) radiusOld, velOld
        logical :: aux

        call s_transfer_data_to_tmp()

        !$acc parallel loop gang vector default(present) reduction(IOR: lag_largestep) private(k) copyin(RKstep)
        do k = 1, nBubs

            radiusOld = lag_bub_interface_rad(k, 2)
            velOld = lag_bub_interface_vel(k, 2)

            !$acc loop seq
            do i = 1, RKstep
                lag_bub_interface_rad(k, 2) = lag_bub_interface_rad(k, 2) + dt*lag_RKcoef(RKstep, i)*lag_bub_interface_draddt(k, i)
                lag_bub_interface_vel(k, 2) = lag_bub_interface_vel(k, 2) + dt*lag_RKcoef(RKstep, i)*lag_bub_interface_dveldt(k, i)
                lag_bub_motion_pos(k, 1:3, 2) = lag_bub_motion_pos(k, 1:3, 2) + dt*lag_RKcoef(RKstep, i)*lag_bub_motion_dposdt(k, 1:3, i)
                lag_bub_motion_vel(k, 1:3, 2) = lag_bub_motion_vel(k, 1:3, 2) + dt*lag_RKcoef(RKstep, i)*lag_bub_motion_dveldt(k, 1:3, i)
                lag_bub_gas_p(k, 2) = lag_bub_gas_p(k, 2) + dt*lag_RKcoef(RKstep, i)*lag_bub_gas_dpdt(k, i)
                lag_bub_gas_mv(k, 2) = lag_bub_gas_mv(k, 2) + dt*lag_RKcoef(RKstep, i)*lag_bub_gas_dmvdt(k, i)
            end do

            if ((lag_bub_interface_rad(k, 2) <= 0.0d0) .or. &                    ! no negative radius
                (lag_bub_motion_pos(k, 1, 2) /= lag_bub_motion_pos(k, 1, 2))) then   ! finite bubble location

                print *, 'Large time step. Radius from:', radiusOld, ' to :', lag_bub_interface_rad(k, 2), &
                    '; and velocity from', velOld, ' to :', lag_bub_interface_vel(k, 2)
                lag_largestep = .true.
                if (dt < 2.d-15) then
                    print *, 'WARNING large step', dt, '. Removing bubble:', lag_bub_id(k, 1)
                    call s_remove_lag_bubble(k)
                end if
                goto 711
            end if

        end do

711     continue
!$acc update host(lag_largestep)
#ifdef MFC_MPI
        if (num_procs > 1) then
            call MPI_ALLREDUCE(lag_largestep, aux, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
            lag_largestep = aux
        end if
#endif

        if (lag_largestep) return

        ! Update background fluid variables
        !$acc parallel loop collapse(4) gang vector default(present)
        do l = 1, sys_size
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        q_cons_ts(2)%vf(l)%sf(i, j, k) = &
                            q_cons_ts(1)%vf(l)%sf(i, j, k)
                    end do
                end do
            end do
        end do

        !$acc parallel loop collapse(4) gang vector default(present) copyin(RKstep)
        do l = 1, sys_size
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        !$acc loop seq
                        do q = 1, RKstep
                            q_cons_ts(2)%vf(l)%sf(i, j, k) = &
                                q_cons_ts(2)%vf(l)%sf(i, j, k) + &
                                dt*lag_RKcoef(RKstep, q)*rhs_ts(q)%vf(l)%sf(i, j, k)
                        end do
                    end do
                end do
            end do
        end do

    end subroutine s_update_tmp_rkck

    !>  This subroutine calculates the maximum error between the 4th and 5th order Runge-Kutta-Cash-Karp solutions
        !!      for the same time step size. If the errors are smaller than a tolerance, then the algorithm employs
        !!      the 5th order solution, while if not, both eulerian/lagrangian variables are re-calculated with a
        !!      smaller time step size.
    subroutine s_calculate_rkck_truncation_error()

        real(kind(0.d0)) :: erraux, errb
        integer :: i, j, k, l, l1, nb

        !$acc parallel loop gang vector default(present) reduction(MAX: lag_errmax) private(k)
        do k = 1, nBubs
            errb = 0.0d0
            if (.not. lag_bub_equilibrium(k)) then
                !Bubble radius error
                erraux = 0.0d0
                !$acc loop seq
                do i = 1, 6
                    erraux = erraux + lag_RKcoef(7, i)*lag_bub_interface_draddt(k, i)
                end do
                errb = max(errb, abs(erraux)*dt/lag_bub_R0(k))

                !Interface velocity error
                erraux = 0.0d0
                !$acc loop seq
                do i = 1, 6
                    erraux = erraux + lag_RKcoef(7, i)*lag_bub_interface_dveldt(k, i)
                end do
                errb = max(errb, abs(erraux)*dt)

                !Bubble velocity error
                !$acc loop seq
                do j = 1, 3
                    erraux = 0.0d0
                    !$acc loop seq
                    do i = 1, 6
                        erraux = erraux + lag_RKcoef(7, i)*lag_bub_motion_dposdt(k, j, i)
                    end do
                    errb = max(errb, abs(erraux)*dt/(abs(lag_bub_motion_vel(k, j, 2)) + 1.0d-4))
                end do
            end if
            lag_errmax = max(lag_errmax, errb)
        end do

        !$acc update host(lag_errmax)

    end subroutine s_calculate_rkck_truncation_error

    !>  This subroutine updates the conservative fields and the lagrangian variables after accepting the performed time step.
        !! @param q_cons_ts Conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_update_rkck(q_cons_ts, q_prim_vf)

        type(vector_field), dimension(:), intent(inout) :: q_cons_ts
        type(scalar_field), dimension(:), intent(in) :: q_prim_vf

        integer :: i, j, k, l

        !$acc parallel loop gang vector default(present) private(k)
        do k = 1, nBubs
            !Accept time step (actual vars = temporal vars)
            lag_bub_motion_pos(k, 1:3, 1) = lag_bub_motion_pos(k, 1:3, 2)
            lag_bub_motion_vel(k, 1:3, 1) = lag_bub_motion_vel(k, 1:3, 2)
            lag_bub_interface_rad(k, 1) = lag_bub_interface_rad(k, 2)
            lag_bub_interface_vel(k, 1) = lag_bub_interface_vel(k, 2)
            lag_bub_gas_p(k, 1) = lag_bub_gas_p(k, 2)
            lag_bub_gas_mv(k, 1) = lag_bub_gas_mv(k, 2)
        end do

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    do l = 1, sys_size
                        q_cons_ts(1)%vf(l)%sf(i, j, k) = q_cons_ts(2)%vf(l)%sf(i, j, k)
                    end do
                end do
            end do
        end do

        call s_smear_voidfraction()

    end subroutine s_update_rkck

    !> This function performs a bilinear interpolation.
          !! @param coord Interpolation coordinates
          !! @param q Input scalar field
    function f_interpolate(coord, q)
        !$acc routine seq
        type(scalar_field), intent(in) :: q
        real(kind(0.d0)), dimension(3), intent(in) :: coord

        real(kind(0.d0)) :: f_interpolate, tmp
        real(kind(0.d0)), dimension(3) :: psi !local coordinates
        integer, dimension(3) :: cell

        call s_get_psi(coord, psi, cell)

        if (p == 0) then  !2D
            tmp = q%sf(cell(1), cell(2), cell(3))*(1.0d0 - psi(1))*(1.0d0 - psi(2))
            tmp = tmp + q%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1.0d0 - psi(2))
            tmp = tmp + q%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)
            tmp = tmp + q%sf(cell(1), cell(2) + 1, cell(3))*(1.0d0 - psi(1))*psi(2)
        else              !3D
            tmp = q%sf(cell(1), cell(2), cell(3))*(1.0d0 - psi(1))*(1.0d0 - psi(2))*(1.0d0 - psi(3))
            tmp = tmp + q%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1.0d0 - psi(2))*(1.0d0 - psi(3))
            tmp = tmp + q%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)*(1.0d0 - psi(3))
            tmp = tmp + q%sf(cell(1), cell(2) + 1, cell(3))*(1.0d0 - psi(1))*psi(2)*(1.0d0 - psi(3))
            tmp = tmp + q%sf(cell(1), cell(2), cell(3) + 1)*(1.0d0 - psi(1))*(1.0d0 - psi(2))*psi(3)
            tmp = tmp + q%sf(cell(1) + 1, cell(2), cell(3) + 1)*psi(1)*(1.0d0 - psi(2))*psi(3)
            tmp = tmp + q%sf(cell(1) + 1, cell(2) + 1, cell(3) + 1)*psi(1)*psi(2)*psi(3)
            tmp = tmp + q%sf(cell(1), cell(2) + 1, cell(3) + 1)*(1.0d0 - psi(1))*psi(2)*psi(3)
        end if

        f_interpolate = tmp

    end function f_interpolate

    !> This subroutine returns the computational coordinate of the cell for the given position.
          !! @param pos Input coordinates
          !! @param cell Computational coordinate of the cell
          !! @param scoord Calculated particle coordinates
    subroutine s_locate_cell(pos, cell, scoord)

        real(kind(0.d0)), dimension(3) :: pos
        real(kind(0.d0)), dimension(3), optional :: scoord
        integer, dimension(3) :: cell
        integer :: i, j, k

        do while (pos(1) < x_cb(cell(1) - 1))
            cell(1) = cell(1) - 1
        end do

        do while (pos(1) > x_cb(cell(1)))
            cell(1) = cell(1) + 1
        end do

        do while (pos(2) < y_cb(cell(2) - 1))
            cell(2) = cell(2) - 1
        end do

        do while (pos(2) > y_cb(cell(2)))
            cell(2) = cell(2) + 1
        end do

        if (p > 0) then
            do while (pos(3) < z_cb(cell(3) - 1))
                cell(3) = cell(3) - 1
            end do
            do while (pos(3) > z_cb(cell(3)))
                cell(3) = cell(3) + 1
            end do
        end if

        ! The numbering of the cell of which left boundary is the domain boundary is 0.
        ! if comp.coord of the pos is s, the real coordinate of s is
        ! (the coordinate of the left boundary of the Floor(s)-th cell)
        ! + (s-(int(s))*(cell-width).
        ! In other words,  the coordinate of the center of the cell is x_cc(cell).

        !coordinates in computational space
        if (present(scoord)) then
            scoord(1) = cell(1) + (pos(1) - x_cb(cell(1) - 1))/dx(cell(1))
            scoord(2) = cell(2) + (pos(2) - y_cb(cell(2) - 1))/dy(cell(2))
            scoord(3) = 0.0d0
            if (p > 0) scoord(3) = cell(3) + (pos(3) - z_cb(cell(3) - 1))/dz(cell(3))
            call s_get_cell(scoord, cell)
        end if

    end subroutine s_locate_cell

    !> This subroutine transfer data into the temporal variables in the lagrangian solver
          !! @param particle_data  Particle data
    subroutine s_transfer_data_to_tmp()

        integer :: k

        !$acc parallel loop gang vector default(present) private(k)
        do k = 1, nBubs
            lag_bub_gas_p(k, 2) = lag_bub_gas_p(k, 1)
            lag_bub_gas_mv(k, 2) = lag_bub_gas_mv(k, 1)
            lag_bub_interface_rad(k, 2) = lag_bub_interface_rad(k, 1)
            lag_bub_interface_vel(k, 2) = lag_bub_interface_vel(k, 1)
            lag_bub_motion_pos(k, 1:3, 2) = lag_bub_motion_pos(k, 1:3, 1)
            lag_bub_motion_posPrev(k, 1:3, 2) = lag_bub_motion_posPrev(k, 1:3, 1)
            lag_bub_motion_vel(k, 1:3, 2) = lag_bub_motion_vel(k, 1:3, 1)
            lag_bub_motion_s(k, 1:3, 2) = lag_bub_motion_s(k, 1:3, 1)
        end do

    end subroutine s_transfer_data_to_tmp

    !> The purpose of this procedure is to determine if the global coordinates of the bubbles
        !!      are present in the current MPI processor (including ghost cells).
        !! @param pos_part Spatial coordinates of the bubble
    function particle_in_domain(pos_part)

        logical :: particle_in_domain
        real(kind(0.d0)), dimension(3) :: pos_part

        ! 2D
        if (p == 0 .and. cyl_coord .neqv. .true.) then
            ! Defining a virtual z-axis that has the same dimensions as y-axis
            ! defined in the input file
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
                                  (pos_part(3) < lag_charwidth/2d0) .and. (pos_part(3) >= -lag_charwidth/2d0))
        else
            ! cyl_coord
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (abs(pos_part(2)) < y_cb(n + buff_size)) .and. (abs(pos_part(2)) >= max(y_cb(-buff_size - 1), 0d0)))
        end if

        ! 3D
        if (p > 0) then
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
                                  (pos_part(3) < z_cb(p + buff_size)) .and. (pos_part(3) >= z_cb(-buff_size - 1)))
        end if

        ! For symmetric boundary condition
        if (bc_x%beg == -2) then
            particle_in_domain = (particle_in_domain .and. (pos_part(1) >= x_cb(-1)))
        end if
        if (bc_x%end == -2) then
            particle_in_domain = (particle_in_domain .and. (pos_part(1) < x_cb(m)))
        end if
        if (bc_y%beg == -2 .and. (.not. cyl_coord)) then
            particle_in_domain = (particle_in_domain .and. (pos_part(2) >= y_cb(-1)))
        end if
        if (bc_y%end == -2 .and. (.not. cyl_coord)) then
            particle_in_domain = (particle_in_domain .and. (pos_part(2) < y_cb(n)))
        end if

        if (p > 0) then
            if (bc_z%beg == -2) then
                particle_in_domain = (particle_in_domain .and. (pos_part(3) >= z_cb(-1)))
            end if
            if (bc_z%end == -2) then
                particle_in_domain = (particle_in_domain .and. (pos_part(3) < z_cb(p)))
            end if
        end if

    end function particle_in_domain

    !> The purpose of this procedure is to determine if the lagrangian bubble is located in the
        !!       physical domain. The ghost cells are not part of the physical domain.
        !! @param pos_part Spatial coordinates of the bubble
    function particle_in_domain_physical(pos_part)

        logical :: particle_in_domain_physical
        real(kind(0.d0)), dimension(3) :: pos_part

        particle_in_domain_physical = ((pos_part(1) < x_cb(m)) .and. (pos_part(1) >= x_cb(-1)) .and. &
                                       (pos_part(2) < y_cb(n)) .and. (pos_part(2) >= y_cb(-1)))

        if (p > 0) then
            particle_in_domain_physical = (particle_in_domain_physical .and. (pos_part(3) < z_cb(p)) .and. (pos_part(3) >= z_cb(-1)))
        end if

    end function particle_in_domain_physical

    !> The purpose of this procedure is to calculate the gradient of a scalar field along the x, y and z directions
        !!      following a second-order central difference considering uneven widths
        !! @param q Input scalar field
        !! @param dq Output gradient of q
        !! @param dir Gradient spatial direction
    subroutine s_gradient_dir(q, dq, dir)

        type(scalar_field), intent(inout) :: q
        type(scalar_field), intent(inout) :: dq
        integer, intent(in) :: dir

        integer :: i, j, k, l, lmax
        real(kind(0.d0)) :: aux1, aux2

        if (dir == 1) then
            ! Gradient in x dir.
            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        aux1 = dx(i) + dx(i - 1)
                        aux2 = dx(i) + dx(i + 1)
                        dq%sf(i, j, k) = q%sf(i, j, k)*(dx(i + 1) - dx(i - 1)) &
                                         + q%sf(i + 1, j, k)*aux1 &
                                         - q%sf(i - 1, j, k)*aux2
                        dq%sf(i, j, k) = dq%sf(i, j, k)/(aux1*aux2)
                    end do
                end do
            end do
        else
            if (dir == 2) then
                ! Gradient in y dir.
                !$acc parallel loop collapse(3) gang vector default(present)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            aux1 = dy(j) + dy(j - 1)
                            aux2 = dy(j) + dy(j + 1)
                            dq%sf(i, j, k) = q%sf(i, j, k)*(dy(j + 1) - dy(j - 1)) &
                                             + q%sf(i, j + 1, k)*aux1 &
                                             - q%sf(i, j - 1, k)*aux2
                            dq%sf(i, j, k) = dq%sf(i, j, k)/(aux1*aux2)
                        end do
                    end do
                end do
            else
                ! Gradient in z dir.
                !$acc parallel loop collapse(3) gang vector default(present)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            aux1 = dz(k) + dz(k - 1)
                            aux2 = dz(k) + dz(k + 1)
                            dq%sf(i, j, k) = q%sf(i, j, k)*(dz(k + 1) - dz(k - 1)) &
                                             + q%sf(i, j, k + 1)*aux1 &
                                             - q%sf(i, j, k - 1)*aux2
                            dq%sf(i, j, k) = dq%sf(i, j, k)/(aux1*aux2)
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_gradient_dir

    !> Subroutine that writes on each time step the changes of the lagrangian bubbles.
        !!  @param q_time Current time
    subroutine s_write_lag_particles(qtime)

        real(kind(0.d0)) :: qtime, pinf
        integer :: i, j, k
        integer, dimension(3) :: cell

        character(LEN=path_len + 2*name_len) :: file_loc

        write (file_loc, '(A,I0,A)') 'lag_bubble_evol_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)

        if (qtime == 0.d0) then
            open (11, FILE=trim(file_loc), FORM='formatted', position='rewind')
            write (11, *) 'currentTime, particleID, x, y, z, ', &
                'coreVaporMass, coreVaporConcentration, radius, interfaceVelocity, ', &
                'corePressure'
        else
            open (11, FILE=trim(file_loc), FORM='formatted', position='append')
        end if

        ! Cycle through list
        do k = 1, nBubs
            write (11, '(6X,f12.6,I24.8,8e24.8)') &
                qtime, &
                lag_bub_id(k, 1), &
                lag_bub_motion_pos(k, 1, 1), &
                lag_bub_motion_pos(k, 2, 1), &
                lag_bub_motion_pos(k, 3, 1), &
                lag_bub_gas_mv(k, 1), &
                lag_bub_gas_mv(k, 1)/(lag_bub_gas_mv(k, 1) + lag_bub_gas_mg(k)), &
                lag_bub_interface_rad(k, 1), &
                lag_bub_interface_vel(k, 1), &
                lag_bub_gas_p(k, 1)
        end do

        close (11)

    end subroutine s_write_lag_particles

    !>  Subroutine that writes some useful statistics related to the volume fraction
            !!       of the particles (void fraction) in the computatioational domain
            !!       on each time step.
            !!  @param q_time Current time
    subroutine s_write_void_evol(qtime)

        real(kind(0.d0)) :: qtime, volcell, voltot
        real(kind(0.d0)) :: voidmax_glb, voidavg_glb, vol_glb
        integer :: i, j, k
        integer, dimension(3) :: cell
        logical :: prevfile

        character(LEN=path_len + 2*name_len) :: file_loc

        if (proc_rank == 0) then
            write (file_loc, '(A)') 'voidfraction.dat'
            file_loc = trim(case_dir)//'/D/'//trim(file_loc)
            if (qtime == 0.d0) then
                open (12, FILE=trim(file_loc), FORM='formatted', position='rewind')
                !write (12, *) 'currentTime, averageVoidFraction, ', &
                !    'maximumVoidFraction, totalParticlesVolume'
                !write (12, *) 'The averageVoidFraction value does ', &
                !    'not reflect the real void fraction in the cloud since the ', &
                !    'cells which do not have bubbles are not accounted'
            else
                open (12, FILE=trim(file_loc), FORM='formatted', position='append')
            end if
        end if

        lag_voidmax = 0.0d0
        lag_voidavg = 0.0d0
        lag_vol = 0.0d0
        !$acc update device(lag_voidmax, lag_voidavg, lag_vol)

        !$acc parallel loop collapse(3) gang vector default(present) reduction(+:lag_vol,lag_voidavg) &
        !$acc reduction(MAX:lag_voidmax) private(cell)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    lag_voidmax = max(lag_voidmax, 1.0d0 - q_particle%vf(1)%sf(i, j, k))
                    cell(1) = i
                    cell(2) = j
                    cell(3) = k
                    call s_get_char_vol(cell(1), cell(2), cell(3), volcell)
                    if ((1.0d0 - q_particle%vf(1)%sf(i, j, k)) > 5.0d-11) then
                        lag_voidavg = lag_voidavg + (1.0d0 - q_particle%vf(1)%sf(i, j, k))*volcell
                        lag_vol = lag_vol + volcell
                    end if
                end do
            end do
        end do

!$acc update host(lag_voidmax, lag_voidavg, lag_vol)
#ifdef MFC_MPI
        if (num_procs > 1) then
            call s_mpi_allreduce_max(lag_voidmax, voidmax_glb)
            lag_voidmax = voidmax_glb
            call s_mpi_allreduce_sum(lag_vol, vol_glb)
            lag_vol = vol_glb
            call s_mpi_allreduce_sum(lag_voidavg, voidavg_glb)
            lag_voidavg = voidavg_glb
        end if
#endif
        voltot = lag_voidavg
        ! This voidavg value does not reflect the real void fraction in the cloud
        ! since the cell which does not have bubbles are not accounted
        if (lag_vol > 0.0d0) lag_voidavg = lag_voidavg/lag_vol

        if (proc_rank == 0) then
            write (12, '(6X,4e24.6)') &
                qtime, &
                lag_voidavg, &
                voltot, &
                lag_voidmax
            close (12)
        end if

    end subroutine s_write_void_evol

    !>  Subroutine that writes the restarting files for the particles in the lagrangian solver.
        !!  @param t_step Current time step
    subroutine s_write_restart_lag_bubbles(t_step)

        ! Generic string used to store the address of a particular file
        character(LEN=path_len + 2*name_len) :: file_loc

        ! Generic logical used for purpose of asserting whether a particular
        ! directory is or is not located in the designated location
        logical :: dir_check
        logical :: file_exist

        integer :: i, k, t_step
        integer :: nparticles, tot_part, tot_part_wrtn, npart_wrtn

#ifdef MFC_MPI
        ! For Parallel I/O
        integer :: ifile, ireq, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer :: view
        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer, dimension(num_procs) :: part_order, part_ord_mpi

        nparticles = 0.0d0
        if (nBubs /= 0) then
            do k = 1, nBubs
                if (particle_in_domain_physical(lag_bub_motion_pos(k, 1:3, 1))) then
                    nparticles = nparticles + 1
                end if
            end do
        end if

        ! Total number of particles
        call MPI_ALLREDUCE(nparticles, tot_part, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        ! Total number of particles written so far
        call MPI_ALLREDUCE(npart_wrtn, tot_part_wrtn, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        lsizes(1) = max(1, nparticles)
        lsizes(2) = 21

        ! if the partcle number is zero, put 1 since MPI cannot deal with writing
        ! zero particle
        part_order(:) = 1
        part_order(proc_rank + 1) = max(1, nparticles)

        call MPI_ALLREDUCE(part_order, part_ord_mpi, num_procs, MPI_integer, &
                           MPI_MAX, MPI_COMM_WORLD, ierr)

        gsizes(1) = sum(part_ord_mpi(1:num_procs))
        gsizes(2) = 21

        start_idx_part(1) = sum(part_ord_mpi(1:proc_rank + 1)) - part_ord_mpi(proc_rank + 1)
        start_idx_part(2) = 0

        write (file_loc, '(A,I0,A)') 'lag_bubbles_mpi_io_', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) then
            call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if

        ! Writing down the total number of particles
        if (proc_rank == 0) then
            open (9, FILE=trim(file_loc), FORM='unformatted', STATUS='unknown')
            write (9) gsizes(1), mytime, dt
            close (9)
        end if

        call MPI_type_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                      MPI_ORDER_FORTRAN, MPI_doUBLE_PRECISION, view, ierr)
        call MPI_type_COMMIT(view, ierr)

        allocate (MPI_IO_DATA_lag_bubbles(1:max(1, nparticles), 1:21))

        ! Open the file to write all flow variables
        write (file_loc, '(A,I0,A)') 'lag_bubbles_', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) then
            call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if

        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                           mpi_info_int, ifile, ierr)

        disp = 0d0

        call MPI_FILE_SET_VIEW(ifile, disp, MPI_doUBLE_PRECISION, view, &
                               'native', mpi_info_null, ierr)

        ! Cycle through list
        i = 1

        if (nparticles == 0) then
            MPI_IO_DATA_lag_bubbles(1, 1:21) = 0d0
        else

            do k = 1, nBubs

                if (particle_in_domain_physical(lag_bub_motion_pos(k, 1:3, 1))) then

                    MPI_IO_DATA_lag_bubbles(i, 1) = real(lag_bub_id(k, 1))
                    MPI_IO_DATA_lag_bubbles(i, 2:4) = lag_bub_motion_pos(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 5:7) = lag_bub_motion_posPrev(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 8:10) = lag_bub_motion_vel(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 11) = lag_bub_interface_rad(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 12) = lag_bub_interface_vel(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 13) = lag_bub_R0(k)
                    MPI_IO_DATA_lag_bubbles(i, 14) = lag_bub_Rmax(k)
                    MPI_IO_DATA_lag_bubbles(i, 15) = lag_bub_Rmin(k)
                    MPI_IO_DATA_lag_bubbles(i, 16) = lag_bub_dphidt(k)
                    MPI_IO_DATA_lag_bubbles(i, 17) = lag_bub_gas_p(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 18) = lag_bub_gas_mv(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 19) = lag_bub_gas_mg(k)
                    MPI_IO_DATA_lag_bubbles(i, 20) = lag_bub_gas_betaT(k)
                    MPI_IO_DATA_lag_bubbles(i, 21) = lag_bub_gas_betaC(k)

                    i = i + 1

                end if

            end do

        end if

        call MPI_FILE_write_ALL(ifile, MPI_IO_DATA_lag_bubbles, 21*max(1, nparticles), &
                                MPI_doUBLE_PRECISION, status, ierr)

        call MPI_FILE_CLOSE(ifile, ierr)

        deallocate (MPI_IO_DATA_lag_bubbles)

#endif

    end subroutine s_write_restart_lag_bubbles

    !>  This procedure calculates the maximum and minimum radius of each bubble.
    subroutine s_calculate_lag_bubble_stats()

        integer :: k

        !$acc parallel loop gang vector default(present) reduction(MAX: Rmax_glb) reduction(MAX: Rmin_glb) private(k)
        do k = 1, nBubs
            Rmax_glb = max(Rmax_glb, lag_bub_interface_rad(k, 1)/lag_bub_R0(k))
            Rmin_glb = min(Rmin_glb, lag_bub_interface_rad(k, 1)/lag_bub_R0(k))
            lag_bub_Rmax(k) = max(lag_bub_Rmax(k), lag_bub_interface_rad(k, 1)/lag_bub_R0(k))
            lag_bub_Rmin(k) = min(lag_bub_Rmin(k), lag_bub_interface_rad(k, 1)/lag_bub_R0(k))
        end do

    end subroutine s_calculate_lag_bubble_stats

    !>  Subroutine that writes the maximum and minimum radius of each bubble.
    subroutine s_write_lag_bubble_stats()

        integer :: i, k
        character(LEN=path_len + 2*name_len) :: file_loc

        write (file_loc, '(A,I0,A)') 'stats_lag_bubbles_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)

        !$acc update host(Rmax_glb, Rmin_glb)

        open (13, FILE=trim(file_loc), FORM='formatted', position='rewind')
        write (13, *) 'proc_rank, particleID, x, y, z, Rmax_glb, Rmin_glb'

        do k = 1, nBubs
            write (13, '(6X,2I24.8,5e24.8)') &
                proc_rank, &
                lag_bub_id(k, 1), &
                lag_bub_motion_pos(k, 1, 1), &
                lag_bub_motion_pos(k, 2, 1), &
                lag_bub_motion_pos(k, 3, 1), &
                lag_bub_Rmax(k), &
                lag_bub_Rmin(k)
        end do

        close (13)

    end subroutine s_write_lag_bubble_stats

    !> The purpose of this subroutine is to remove one specific particle
          !! @param nparticles Particle id
    subroutine s_remove_lag_bubble(nparticles)
        !$acc routine seq
        integer, intent(in) :: nparticles

        integer :: i

        do i = nparticles, nBubs - 1
            lag_bub_id(i, 1) = lag_bub_id(i + 1, 1)
            lag_bub_R0(i) = lag_bub_R0(i + 1)
            lag_bub_Rmax(i) = lag_bub_Rmax(i + 1)
            lag_bub_Rmin(i) = lag_bub_Rmin(i + 1)
            lag_bub_equilibrium(i) = lag_bub_equilibrium(i + 1)
            lag_bub_gas_mg(i) = lag_bub_gas_mg(i + 1)
            lag_bub_gas_betaT(i) = lag_bub_gas_betaT(i + 1)
            lag_bub_gas_betaC(i) = lag_bub_gas_betaC(i + 1)
            lag_bub_dphidt(i) = lag_bub_dphidt(i + 1)
            lag_bub_gas_p(i, 1:2) = lag_bub_gas_p(i + 1, 1:2)
            lag_bub_gas_mv(i, 1:2) = lag_bub_gas_mv(i + 1, 1:2)
            lag_bub_interface_rad(i, 1:2) = lag_bub_interface_rad(i + 1, 1:2)
            lag_bub_interface_vel(i, 1:2) = lag_bub_interface_vel(i + 1, 1:2)
            lag_bub_motion_pos(i, 1:3, 1:2) = lag_bub_motion_pos(i + 1, 1:3, 1:2)
            lag_bub_motion_posPrev(i, 1:3, 1:2) = lag_bub_motion_posPrev(i + 1, 1:3, 1:2)
            lag_bub_motion_vel(i, 1:3, 1:2) = lag_bub_motion_vel(i + 1, 1:3, 1:2)
            lag_bub_motion_s(i, 1:3, 1:2) = lag_bub_motion_s(i + 1, 1:3, 1:2)
            lag_bub_interface_draddt(i, 1:6) = lag_bub_interface_draddt(i + 1, 1:6)
            lag_bub_interface_dveldt(i, 1:6) = lag_bub_interface_dveldt(i + 1, 1:6)
            lag_bub_gas_dpdt(i, 1:6) = lag_bub_gas_dpdt(i + 1, 1:6)
            lag_bub_gas_dmvdt(i, 1:6) = lag_bub_gas_dmvdt(i + 1, 1:6)
        end do

        nBubs = nBubs - 1

    end subroutine s_remove_lag_bubble

    !> The purpose of this subroutine is to deallocate variables
    subroutine s_finalize_lagrangian_solver()

        integer :: i, j, k, imax

        do i = 1, q_particle_idx
            @:DEALLOCATE(q_particle%vf(i)%sf)
        end do
        @:DEALLOCATE(q_particle%vf)

        !Deallocating space
        @:DEALLOCATE_GLOBAL(lag_RKcoef)
        @:DEALLOCATE_GLOBAL(lag_bub_id)
        @:DEALLOCATE_GLOBAL(lag_bub_R0)
        @:DEALLOCATE_GLOBAL(lag_bub_Rmax)
        @:DEALLOCATE_GLOBAL(lag_bub_Rmin)
        @:DEALLOCATE_GLOBAL(lag_bub_equilibrium)
        @:DEALLOCATE_GLOBAL(lag_bub_gas_mg)
        @:DEALLOCATE_GLOBAL(lag_bub_gas_betaT)
        @:DEALLOCATE_GLOBAL(lag_bub_gas_betaC)
        @:DEALLOCATE_GLOBAL(lag_bub_dphidt)
        @:DEALLOCATE_GLOBAL(lag_bub_gas_p)
        @:DEALLOCATE_GLOBAL(lag_bub_gas_mv)
        @:DEALLOCATE_GLOBAL(lag_bub_interface_rad)
        @:DEALLOCATE_GLOBAL(lag_bub_interface_vel)
        @:DEALLOCATE_GLOBAL(lag_bub_motion_pos)
        @:DEALLOCATE_GLOBAL(lag_bub_motion_posPrev)
        @:DEALLOCATE_GLOBAL(lag_bub_motion_vel)
        @:DEALLOCATE_GLOBAL(lag_bub_motion_s)
        @:DEALLOCATE_GLOBAL(lag_bub_interface_draddt)
        @:DEALLOCATE_GLOBAL(lag_bub_interface_dveldt)
        @:DEALLOCATE_GLOBAL(lag_bub_gas_dpdt)
        @:DEALLOCATE_GLOBAL(lag_bub_gas_dmvdt)
        @:DEALLOCATE_GLOBAL(lag_bub_motion_dposdt)
        @:DEALLOCATE_GLOBAL(lag_bub_motion_dveldt)

    end subroutine s_finalize_lagrangian_solver

end module m_lag_bubbles
