!>
!! @file m_bubbles_EL.fpp
!! @brief Contains module m_bubbles_EL

#:include 'macros.fpp'

!> @brief This module is used to to compute the volume-averaged bubble model
module m_bubbles_EL

    use m_global_parameters             !< Definitions of the global parameters

    use m_mpi_proxy                     !< Message passing interface (MPI) module proxy

    use m_bubbles_EL_kernels            !< Definitions of the kernel functions

    use m_bubbles                       !< General bubble dynamics procedures

    use m_variables_conversion          !< State variables type conversion procedures

    use m_compile_specific

    use m_boundary_conditions

    use m_sim_helpers

    use m_helper

    implicit none

    !(nBub)
    integer, allocatable, dimension(:, :) :: lag_id                 !< Global and local IDs
    real(wp), allocatable, dimension(:) :: bub_R0            !< Initial bubble radius
    real(wp), allocatable, dimension(:) :: Rmax_stats        !< Maximum radius
    real(wp), allocatable, dimension(:) :: Rmin_stats        !< Minimum radius
    real(wp), allocatable, dimension(:) :: gas_mg            !< Bubble's gas mass
    real(wp), allocatable, dimension(:) :: gas_betaT         !< heatflux model (Preston et al., 2007)
    real(wp), allocatable, dimension(:) :: gas_betaC         !< massflux model (Preston et al., 2007)
    real(wp), allocatable, dimension(:) :: bub_dphidt        !< subgrid velocity potential (Maeda & Colonius, 2018)
    !(nBub, 1 -> actual val or 2 -> temp val)
    real(wp), allocatable, dimension(:, :) :: gas_p          !< Pressure in the bubble
    real(wp), allocatable, dimension(:, :) :: gas_mv         !< Vapor mass in the bubble
    real(wp), allocatable, dimension(:, :) :: intfc_rad      !< Bubble radius
    real(wp), allocatable, dimension(:, :) :: intfc_vel      !< Velocity of the bubble interface
    !(nBub, 1-> x or 2->y or 3 ->z, 1 -> actual or 2 -> temporal val)
    real(wp), allocatable, dimension(:, :, :) :: mtn_pos     !< Bubble's position
    real(wp), allocatable, dimension(:, :, :) :: mtn_posPrev !< Bubble's previous position
    real(wp), allocatable, dimension(:, :, :) :: mtn_vel     !< Bubble's velocity
    real(wp), allocatable, dimension(:, :, :) :: mtn_s       !< Bubble's computational cell position in real format
    !(nBub, 1-> x or 2->y or 3 ->z, time-stage)
    real(wp), allocatable, dimension(:, :) :: intfc_draddt   !< Time derivative of bubble's radius
    real(wp), allocatable, dimension(:, :) :: intfc_dveldt   !< Time derivative of bubble's interface velocity
    real(wp), allocatable, dimension(:, :) :: gas_dpdt       !< Time derivative of gas pressure
    real(wp), allocatable, dimension(:, :) :: gas_dmvdt      !< Time derivative of the vapor mass in the bubble
    real(wp), allocatable, dimension(:, :, :) :: mtn_dposdt  !< Time derivative of the bubble's position
    real(wp), allocatable, dimension(:, :, :) :: mtn_dveldt  !< Time derivative of the bubble's velocity

    !$acc declare create(lag_id, bub_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, gas_betaC, bub_dphidt,       &
    !$acc gas_p, gas_mv, intfc_rad, intfc_vel, mtn_pos, mtn_posPrev, mtn_vel, mtn_s, intfc_draddt, intfc_dveldt, &
    !$acc gas_dpdt, gas_dmvdt, mtn_dposdt, mtn_dveldt)

    real(wp), allocatable, dimension(:, :) :: lag_RKCKcoef   !< RKCK 4th-5th time stepper coefficients
    integer, private :: lag_num_ts                                  !<  Number of time stages in the time-stepping scheme

    !$acc declare create(lag_RKCKcoef, lag_num_ts)

    integer :: nBubs                            !< Number of bubbles in the local domain
    real(wp) :: Rmax_glb, Rmin_glb       !< Maximum and minimum bubbe size in the local domain
    type(vector_field) :: q_beta                !< Projection of the lagrangian particles in the Eulerian framework
    integer :: q_beta_idx                       !< Size of the q_beta vector field

    !$acc declare create(nBubs, Rmax_glb, Rmin_glb, q_beta, q_beta_idx)

contains

    !> Initializes the lagrangian subgrid bubble solver
        !! @param q_cons_vf Initial conservative variables
    subroutine s_initialize_bubbles_EL_module(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        integer :: nBubs_glb, i

        ! Setting number of time-stages for selected time-stepping scheme
        lag_num_ts = time_stepper
        if (time_stepper == 4) lag_num_ts = num_ts_rkck

        ! Allocate space for the Eulerian fields needed to map the effect of the bubbles
        if (lag_params%solver_approach == 1) then
            ! One-way coupling
            q_beta_idx = 3
        elseif (lag_params%solver_approach == 2) then
            ! Two-way coupling
            q_beta_idx = 4
            if (lag_params%cluster_type >= 4) then
                !Subgrid noise model for 2D approximation
                q_beta_idx = 6
            end if
        else
            call s_mpi_abort('Please check the lag_params%solver_approach input')
        end if

        !$acc update device(lag_num_ts, q_beta_idx)

        @:ALLOCATE(q_beta%vf(1:q_beta_idx))

        do i = 1, q_beta_idx
            @:ALLOCATE(q_beta%vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
        end do

        @:ACC_SETUP_VFs(q_beta)

        ! Allocating space for lagrangian variables
        nBubs_glb = lag_params%nBubs_glb

        @:ALLOCATE(lag_id(1:nBubs_glb, 1:2))
        @:ALLOCATE(bub_R0(1:nBubs_glb))
        @:ALLOCATE(Rmax_stats(1:nBubs_glb))
        @:ALLOCATE(Rmin_stats(1:nBubs_glb))
        @:ALLOCATE(gas_mg(1:nBubs_glb))
        @:ALLOCATE(gas_betaT(1:nBubs_glb))
        @:ALLOCATE(gas_betaC(1:nBubs_glb))
        @:ALLOCATE(bub_dphidt(1:nBubs_glb))
        @:ALLOCATE(gas_p(1:nBubs_glb, 1:2))
        @:ALLOCATE(gas_mv(1:nBubs_glb, 1:2))
        @:ALLOCATE(intfc_rad(1:nBubs_glb, 1:2))
        @:ALLOCATE(intfc_vel(1:nBubs_glb, 1:2))
        @:ALLOCATE(mtn_pos(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_posPrev(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_vel(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_s(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(intfc_draddt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(intfc_dveldt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(gas_dpdt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(gas_dmvdt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(mtn_dposdt(1:nBubs_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(mtn_dveldt(1:nBubs_glb, 1:3, 1:lag_num_ts))

        if (time_stepper == 4) then
            !< Allocate space for the RKCK 4th/5th time stepper coefficients
            @:ALLOCATE(lag_RKCKcoef(1:lag_num_ts+1, 1:lag_num_ts))
            do i = 1, lag_num_ts
                ! Populate RKCK coefficients (from constants)
                lag_RKCKcoef(1, i) = rkck_coef1(i)
                lag_RKCKcoef(2, i) = rkck_coef2(i)
                lag_RKCKcoef(3, i) = rkck_coef3(i)
                lag_RKCKcoef(4, i) = rkck_coef4(i)
                lag_RKCKcoef(5, i) = rkck_coef5(i)
                lag_RKCKcoef(6, i) = rkck_coef6(i)
                lag_RKCKcoef(7, i) = rkck_coefE(i)
            end do
        end if

        ! Starting bubbles
        call s_start_lagrange_inputs()
        call s_read_input_bubbles(q_cons_vf)

    end subroutine s_initialize_bubbles_EL_module

    !> The purpose of this procedure is to start lagrange bubble parameters applying nondimensionalization if needed
    subroutine s_start_lagrange_inputs()

        integer :: id_bubbles, id_host
        real(wp) :: rho0, c0, T0, x0, p0

        id_bubbles = num_fluids
        id_host = num_fluids - 1

        !Reference values
        rho0 = lag_params%rho0
        c0 = lag_params%c0
        T0 = lag_params%T0
        x0 = lag_params%x0
        p0 = rho0*c0*c0

        !Update inputs
        Tw = lag_params%Thost/T0
        pv = fluid_pp(id_host)%pv/p0
        gamma_v = fluid_pp(id_bubbles)%gamma_v
        gamma_n = fluid_pp(id_host)%gamma_v
        k_vl = fluid_pp(id_bubbles)%k_v*(T0/(x0*rho0*c0*c0*c0))
        k_nl = fluid_pp(id_host)%k_v*(T0/(x0*rho0*c0*c0*c0))
        cp_v = fluid_pp(id_bubbles)%cp_v*(T0/(c0*c0))
        cp_n = fluid_pp(id_host)%cp_v*(T0/(c0*c0))
        R_v = (R_uni/fluid_pp(id_bubbles)%M_v)*(T0/(c0*c0))
        R_n = (R_uni/fluid_pp(id_host)%M_v)*(T0/(c0*c0))
        lag_params%diffcoefvap = lag_params%diffcoefvap/(x0*c0)
        ss = fluid_pp(id_host)%ss/(rho0*x0*c0*c0)
        mul0 = fluid_pp(id_host)%mul0/(rho0*x0*c0)

        ! Parameters used in bubble_model
        Web = 1._wp/ss
        Re_inv = mul0

        ! Need improvements to accept polytropic gas compression, isothermal and adiabatic thermal models, and
        ! the Gilmore and RP bubble models.
        polytropic = .false.    ! Forcing no polytropic model
        thermal = 3             ! Forcing constant transfer coefficient model based on Preston et al., 2007
        ! If Keller-Miksis model is not selected, then no radial motion

    end subroutine s_start_lagrange_inputs

    !> The purpose of this procedure is to obtain the initial bubbles' information
        !! @param q_cons_vf Conservative variables
    subroutine s_read_input_bubbles(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        real(wp), dimension(8) :: inputBubble
        real(wp) :: qtime
        integer :: id, bub_id, save_count
        integer :: i, ios
        logical :: file_exist, indomain

        character(LEN=path_len + 2*name_len) :: path_D_dir !<

        ! Initialize number of particles
        bub_id = 0
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
            if (proc_rank == 0) print *, 'Reading lagrange bubbles input file.'
            inquire (file='input/lag_bubbles.dat', exist=file_exist)
            if (file_exist) then
                open (94, file='input/lag_bubbles.dat', form='formatted', iostat=ios)
                do while (ios == 0)
                    read (94, *, iostat=ios) (inputBubble(i), i=1, 8)
                    if (ios /= 0) cycle
                    indomain = particle_in_domain(inputBubble(1:3))
                    id = id + 1
                    if (id > lag_params%nBubs_glb .and. proc_rank == 0) then
                        call s_mpi_abort('Current number of bubbles is larger than nBubs_glb')
                    end if
                    if (indomain) then
                        bub_id = bub_id + 1
                        call s_add_bubbles(inputBubble, q_cons_vf, bub_id)
                        lag_id(bub_id, 1) = id      !global ID
                        lag_id(bub_id, 2) = bub_id  !local ID
                        nBubs = bub_id              ! local number of bubbles
                    end if
                end do
                close (94)
            else
                stop "if you include lagrange bubbles, you have to initialize them in input/lag_bubbles.dat"
            end if
        else
            if (proc_rank == 0) print *, 'Restarting lagrange bubbles at save_count: ', save_count
            call s_restart_bubbles(bub_id, save_count)
        end if

        print *, " Lagrange bubbles running, in proc", proc_rank, "number:", bub_id, "/", id

        !$acc update device(bubbles_lagrange, lag_params)

        !$acc update device(lag_id, bub_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, gas_betaC,   &
        !$acc bub_dphidt, gas_p, gas_mv, intfc_rad, intfc_vel, mtn_pos, mtn_posPrev, mtn_vel,       &
        !$acc mtn_s, intfc_draddt, intfc_dveldt, gas_dpdt, gas_dmvdt, mtn_dposdt, mtn_dveldt, nBubs)

        Rmax_glb = min(dflt_real, -dflt_real)
        Rmin_glb = max(dflt_real, -dflt_real)
        !$acc update device(Rmax_glb, Rmin_glb)

        if (time_stepper == 4) then
            !Initial and largest dt - rkck stepper
            dt_max = dt
            !$acc update device(lag_RKCKcoef)
        end if

        !$acc update device(dx, dy, dz, x_cb, x_cc, y_cb, y_cc, z_cb, z_cc)

        !Populate temporal variables
        call s_transfer_data_to_tmp()
        call s_smear_voidfraction()

        if (lag_params%write_bubbles) call s_write_lag_particles(qtime)

        if (save_count == 0) then
            ! Create ./D directory
            write (path_D_dir, '(A,I0,A,I0)') trim(case_dir)//'/D'
            call my_inquire(path_D_dir, file_exist)
            if (.not. file_exist) call s_create_directory(trim(path_D_dir))
            call s_write_restart_lag_bubbles(save_count) ! Needed for post_processing
            call s_write_void_evol(qtime)
        end if

    end subroutine s_read_input_bubbles

    !> The purpose of this procedure is to obtain the information of the bubbles when starting fresh
        !! @param inputBubble Bubble information
        !! @param q_cons_vf Conservative variables
        !! @param bub_id Local id of the bubble
    subroutine s_add_bubbles(inputBubble, q_cons_vf, bub_id)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), dimension(8), intent(in) :: inputBubble
        integer, intent(in) :: bub_id
        integer :: i

        real(wp) :: pliq, volparticle, concvap, totalmass, kparticle, cpparticle
        real(wp) :: omegaN, PeG, PeT, rhol, pcrit, qv, gamma, pi_inf, dynP
        integer, dimension(3) :: cell
        real(wp), dimension(2) :: Re
        real(wp) :: massflag, heatflag, Re_trans, Im_trans

        massflag = 0._wp
        heatflag = 0._wp
        if (lag_params%massTransfer_model) massflag = 1._wp
        if (lag_params%heatTransfer_model) heatflag = 1._wp

        bub_R0(bub_id) = inputBubble(7)
        Rmax_stats(bub_id) = min(dflt_real, -dflt_real)
        Rmin_stats(bub_id) = max(dflt_real, -dflt_real)
        bub_dphidt(bub_id) = 0._wp
        intfc_rad(bub_id, 1) = inputBubble(7)
        intfc_vel(bub_id, 1) = inputBubble(8)
        mtn_pos(bub_id, 1:3, 1) = inputBubble(1:3)
        mtn_posPrev(bub_id, 1:3, 1) = mtn_pos(bub_id, 1:3, 1)
        mtn_vel(bub_id, 1:3, 1) = inputBubble(4:6)

        if (cyl_coord .and. p == 0) then
            mtn_pos(bub_id, 2, 1) = sqrt(mtn_pos(bub_id, 2, 1)**2._wp + &
                                         mtn_pos(bub_id, 3, 1)**2._wp)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            mtn_pos(bub_id, 3, 1) = atan2(inputBubble(3), inputBubble(2))
            mtn_posPrev(bub_id, 1:3, 1) = mtn_pos(bub_id, 1:3, 1)
        end if

        cell = -buff_size
        call s_locate_cell(mtn_pos(bub_id, 1:3, 1), cell, mtn_s(bub_id, 1:3, 1))

        ! If particle is in the ghost cells, find the closest non-ghost cell
        cell(1) = min(max(cell(1), 0), m)
        cell(2) = min(max(cell(2), 0), n)
        if (p > 0) cell(3) = min(max(cell(3), 0), p)
        call s_convert_to_mixture_variables(q_cons_vf, cell(1), cell(2), cell(3), &
                                            rhol, gamma, pi_inf, qv, Re)
        dynP = 0._wp
        do i = 1, num_dims
            dynP = dynP + 0.5_wp*q_cons_vf(contxe + i)%sf(cell(1), cell(2), cell(3))**2/rhol
        end do
        pliq = (q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)) - dynP - pi_inf)/gamma
        if (pliq < 0) print *, "Negative pressure", proc_rank, &
            q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)), pi_inf, gamma, pliq, cell, dynP

        ! Initial particle pressure
        gas_p(bub_id, 1) = pliq + 2._wp*(1._wp/Web)/bub_R0(bub_id)
        if ((1._wp/Web) /= 0._wp) then
            pcrit = pv - 4._wp*(1._wp/Web)/(3._wp*sqrt(3._wp*gas_p(bub_id, 1)*bub_R0(bub_id)**3._wp/(2._wp*(1._wp/Web))))
            pref = gas_p(bub_id, 1)
        else
            pcrit = 0._wp
        end if

        ! Initial particle mass
        volparticle = 4._wp/3._wp*pi*bub_R0(bub_id)**3._wp ! volume
        gas_mv(bub_id, 1) = pv*volparticle*(1._wp/(R_v*Tw))*(massflag) ! vapermass
        gas_mg(bub_id) = (gas_p(bub_id, 1) - pv*(massflag))*volparticle*(1._wp/(R_n*Tw)) ! gasmass
        if (gas_mg(bub_id) <= 0._wp) then
            stop 'the initial mass of gas inside the bubble is negative. Check your initial conditions'
        end if
        totalmass = gas_mg(bub_id) + gas_mv(bub_id, 1) ! totalmass

        ! Bubble natural frequency
        concvap = gas_mv(bub_id, 1)/(gas_mv(bub_id, 1) + gas_mg(bub_id))
        omegaN = (3._wp*(gas_p(bub_id, 1) - pv*(massflag)) + 4._wp*(1._wp/Web)/bub_R0(bub_id))/rhol
        if (pv*(massflag) > gas_p(bub_id, 1)) then
            print *, 'Not allowed: bubble initially located in a region with pressure below the vapor pressure'
            print *, 'location:', mtn_pos(bub_id, 1:3, 1)
            stop
        end if
        omegaN = sqrt(omegaN/bub_R0(bub_id)**2._wp)

        cpparticle = concvap*cp_v + (1._wp - concvap)*cp_n
        kparticle = concvap*k_vl + (1._wp - concvap)*k_nl

        ! Mass and heat transfer coefficients (based on Preston 2007)
        PeT = totalmass/volparticle*cpparticle*bub_R0(bub_id)**2._wp*omegaN/kparticle
        call s_transcoeff(1._wp, PeT, Re_trans, Im_trans)
        gas_betaT(bub_id) = Re_trans*(heatflag)*kparticle

        PeG = bub_R0(bub_id)**2._wp*omegaN/lag_params%diffcoefvap
        call s_transcoeff(1._wp, PeG, Re_trans, Im_trans)
        gas_betaC(bub_id) = Re_trans*(massflag)*lag_params%diffcoefvap

        if (gas_mg(bub_id) <= 0._wp) stop "Negative gas mass in the bubble, check if the bubble is in the domain."

    end subroutine s_add_bubbles

    !> The purpose of this procedure is to obtain the information of the bubbles from a restart point.
        !! @param bub_id Local ID of the particle
        !! @param save_count File identifier
    subroutine s_restart_bubbles(bub_id, save_count)

        integer, intent(inout) :: bub_id, save_count

        character(LEN=path_len + 2*name_len) :: file_loc

#ifdef MFC_MPI
        real(wp), dimension(20) :: inputvals
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: disp
        integer :: view

        integer, dimension(3) :: cell
        logical :: indomain, particle_file, file_exist

        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer :: ifile, ierr, tot_data, id
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
            print '(a)', trim(file_loc)//' is missing. exiting.'
            call s_mpi_abort
        end if

        call MPI_BCAST(tot_data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mytime, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

        gsizes(1) = tot_data
        gsizes(2) = 21
        lsizes(1) = tot_data
        lsizes(2) = 21
        start_idx_part(1) = 0
        start_idx_part(2) = 0

        call MPI_type_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                      MPI_ORDER_FORTRAN, mpi_p, view, ierr)
        call MPI_type_COMMIT(view, ierr)

        ! Open the file to write all flow variables
        write (file_loc, '(a,i0,a)') 'lag_bubbles_', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (file=trim(file_loc), exist=particle_file)

        if (particle_file) then
            call MPI_FILE_open(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)
            disp = 0._wp
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, &
                                   'native', mpi_info_null, ierr)
            allocate (MPI_IO_DATA_lag_bubbles(tot_data, 1:21))
            call MPI_FILE_read_ALL(ifile, MPI_IO_DATA_lag_bubbles, 21*tot_data, &
                                   mpi_p, status, ierr)
            do i = 1, tot_data
                id = int(MPI_IO_DATA_lag_bubbles(i, 1))
                inputvals(1:20) = MPI_IO_DATA_lag_bubbles(i, 2:21)
                indomain = particle_in_domain(inputvals(1:3))
                if (indomain .and. (id > 0)) then
                    bub_id = bub_id + 1
                    nBubs = bub_id                  ! local number of bubbles
                    lag_id(bub_id, 1) = id          ! global ID
                    lag_id(bub_id, 2) = bub_id      ! local ID
                    mtn_pos(bub_id, 1:3, 1) = inputvals(1:3)
                    mtn_posPrev(bub_id, 1:3, 1) = inputvals(4:6)
                    mtn_vel(bub_id, 1:3, 1) = inputvals(7:9)
                    intfc_rad(bub_id, 1) = inputvals(10)
                    intfc_vel(bub_id, 1) = inputvals(11)
                    bub_R0(bub_id) = inputvals(12)
                    Rmax_stats(bub_id) = inputvals(13)
                    Rmin_stats(bub_id) = inputvals(14)
                    bub_dphidt(bub_id) = inputvals(15)
                    gas_p(bub_id, 1) = inputvals(16)
                    gas_mv(bub_id, 1) = inputvals(17)
                    gas_mg(bub_id) = inputvals(18)
                    gas_betaT(bub_id) = inputvals(19)
                    gas_betaC(bub_id) = inputvals(20)
                    cell = -buff_size
                    call s_locate_cell(mtn_pos(bub_id, 1:3, 1), cell, mtn_s(bub_id, 1:3, 1))
                end if
            end do
            deallocate (MPI_IO_DATA_lag_bubbles)
        end if
        call MPI_FILE_CLOSE(ifile, ierr)
#endif

    end subroutine s_restart_bubbles

    !>  Contains the two-way and one-way Euler-Lagrange coupled algorithm, including the bubble dynamics subroutines.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
        !! @param rhs_vf Calculated change of conservative variables
        !! @param stage Current stage in the time-stepper algorithm
    subroutine s_compute_EL_coupled_solver(q_cons_vf, q_prim_vf, rhs_vf, stage)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer, intent(in) :: stage

        real(wp) :: vaporflux, pliqint
        real(wp) :: preterm1, term2, paux, pint, Romega, term1_fac, Rb
        real(wp) :: conc_v, R_m, gamma_m, fpb, fmass_n, fmass_v
        real(wp) :: fR, fV, fbeta_c, fbeta_t, fR0, fpbdt
        real(wp) :: pinf, aux1, aux2, velint, cson, rhol
        real(wp) :: gamma, pi_inf, qv
        real(wp), dimension(contxe) :: myalpha_rho, myalpha
        real(wp), dimension(2) :: Re
        integer, dimension(3) :: cell

        integer :: i, k, l

        call nvtxStartRange("LAGRANGE-BUBBLE-DYNAMICS")

        !< BUBBLE DYNAMICS

        ! Subgrid p_inf model based on Maeda and Colonius (2018).
        if (lag_params%pressure_corrector) then
            ! Calculate velocity potentials (valid for one bubble per cell)
            !$acc parallel loop gang vector default(present) private(k, cell)
            do k = 1, nBubs
                call s_get_pinf(k, q_prim_vf, 2, paux, cell, preterm1, term2, Romega)
                fR0 = bub_R0(k)
                fR = intfc_rad(k, 2)
                fV = intfc_vel(k, 2)
                fpb = gas_p(k, 2)
                pint = f_cpbw_KM(fR0, fR, fV, fpb)
                pint = pint + 0.5_wp*fV**2._wp
                if (lag_params%cluster_type == 2) then
                    bub_dphidt(k) = (paux - pint) + term2
                    ! Accounting for the potential induced by the bubble averaged over the control volume
                    ! Note that this is based on the incompressible flow assumption near the bubble.
                    Rb = intfc_rad(k, 2)
                    term1_fac = 3._wp/2._wp*(Rb*(Romega**2._wp - Rb**2._wp))/(Romega**3._wp - Rb**3._wp)
                    bub_dphidt(k) = bub_dphidt(k)/(1._wp - term1_fac)
                end if
            end do
        end if

        ! Radial motion model
        if (bubble_model == 2) then
            !$acc parallel loop gang vector default(present) private(k, myalpha_rho, myalpha, Re, cell) copyin(stage)
            do k = 1, nBubs
                ! Keller-Miksis model

                ! Current bubble state
                fpb = gas_p(k, 2)
                fmass_n = gas_mg(k)
                fmass_v = gas_mv(k, 2)
                fR = intfc_rad(k, 2)
                fV = intfc_vel(k, 2)
                fbeta_c = gas_betaC(k)
                fbeta_t = gas_betaT(k)
                fR0 = bub_R0(k)

                ! Mixture properties in the bubble
                conc_v = 0._wp
                if (lag_params%massTransfer_model) conc_v = 1._wp/(1._wp + (R_v/R_n)*(fpb/pv - 1._wp))
                R_m = (fmass_n*R_n + fmass_v*R_v)
                gamma_m = conc_v*gamma_v + (1._wp - conc_v)*gamma_n

                ! Vapor and heat fluxes
                vaporflux = f_vflux(fR, fV, fmass_v, k, fmass_n, fbeta_c, conc_v)
                fpbdt = f_bpres_dot(vaporflux, fR, fV, fpb, fmass_v, k, fbeta_t, R_m, gamma_m, conc_v)
                gas_dmvdt(k, stage) = 4._wp*pi*fR**2._wp*vaporflux

                ! Pressure at the bubble wall
                pliqint = f_cpbw_KM(fR0, fR, fV, fpb)

                ! Obtaining driving pressure
                call s_get_pinf(k, q_prim_vf, 1, pinf, cell, aux1, aux2)

                ! Obtain liquid density and computing speed of sound from pinf
                !$acc loop seq
                do i = 1, contxe
                    myalpha_rho(i) = q_prim_vf(advxb + i - 1)%sf(cell(1), cell(2), cell(3))* &
                                     q_prim_vf(i)%sf(cell(1), cell(2), cell(3))
                    myalpha(i) = q_prim_vf(advxb + i - 1)%sf(cell(1), cell(2), cell(3))
                end do
                call s_convert_species_to_mixture_variables_acc(rhol, gamma, pi_inf, qv, myalpha, &
                                                                myalpha_rho, Re, cell(1), cell(2), cell(3))
                call s_compute_cson_from_pinf(k, q_prim_vf, pinf, cell, rhol, gamma, pi_inf, cson)

                ! Velocity correction due to massflux
                velint = fV - gas_dmvdt(k, stage)/(4._wp*pi*fR**2._wp*rhol)

                ! Interphase acceleration and update vars
                intfc_dveldt(k, stage) = f_rddot_KM(fpbdt, pinf, pliqint, rhol, fR, velint, fR0, cson)
                gas_dpdt(k, stage) = fpbdt
                intfc_draddt(k, stage) = fV

            end do
        else
            if (proc_rank == 0) print *, 'WARNING: Lagrange bubbles work with Keller Miksis model!', &
                ' Deactivating radial motion.'
            !$acc parallel loop gang vector default(present) private(k) copyin(stage)
            do k = 1, nBubs
                intfc_dveldt(k, stage) = 0._wp
                intfc_draddt(k, stage) = 0._wp
            end do
        end if

        ! Bubbles remain in a fixed position
        !$acc parallel loop collapse(2) gang vector default(present) private(k) copyin(stage)
        do k = 1, nBubs
            do l = 1, 3
                mtn_dposdt(k, l, stage) = 0._wp
                mtn_dveldt(k, l, stage) = 0._wp
            end do
        end do

        call nvtxEndRange

        !< EULER-LAGRANGE COUPLING
        call s_smear_voidfraction()
        if (lag_params%solver_approach == 2) call s_add_rhs_sources(q_cons_vf, q_prim_vf, rhs_vf)

    end subroutine s_compute_EL_coupled_solver

    !>  This procedure computes the speed of sound from a given driving pressure
        !! @param bub_id Bubble id
        !! @param q_prim_vf Primitive variables
        !! @param pinf Driving pressure
        !! @param cell Bubble cell
        !! @param rhol Liquid density
        !! @param gamma Liquid specific heat ratio
        !! @param pi_inf Liquid stiffness
        !! @param cson Calculated speed of sound
    subroutine s_compute_cson_from_pinf(bub_id, q_prim_vf, pinf, cell, rhol, gamma, pi_inf, cson)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_compute_cson_from_pinf
#else
        !$acc routine seq
#endif
        integer, intent(in) :: bub_id
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), intent(in) :: pinf, rhol, gamma, pi_inf
        integer, dimension(3), intent(in) :: cell
        real(wp), intent(out) :: cson

        real(wp) :: E, H
        real(wp), dimension(num_dims) :: vel
        integer :: i

        !$acc loop seq
        do i = 1, num_dims
            vel(i) = q_prim_vf(i + contxe)%sf(cell(1), cell(2), cell(3))
        end do
        E = gamma*pinf + pi_inf + 0.5_wp*rhol*dot_product(vel, vel)
        H = (E + pinf)/rhol
        cson = sqrt((H - 0.5_wp*dot_product(vel, vel))/gamma)

    end subroutine s_compute_cson_from_pinf

    !>  The purpose of this subroutine is to smear the effect of the bubbles in the Eulerian framework
    subroutine s_smear_voidfraction()

        integer :: i, j, k, l

        call nvtxStartRange("BUBBLES-LAGRANGE-KERNELS")

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, q_beta_idx
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        q_beta%vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do

        call s_smoothfunction(nBubs, intfc_rad, intfc_vel, &
                              mtn_s, mtn_pos, q_beta)

        !Store 1-beta
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    q_beta%vf(1)%sf(j, k, l) = 1._wp - q_beta%vf(1)%sf(j, k, l)
                    ! Limiting void fraction given max value
                    q_beta%vf(1)%sf(j, k, l) = max(q_beta%vf(1)%sf(j, k, l), &
                                                   1._wp - lag_params%valmaxvoid)
                end do
            end do
        end do

        call nvtxEndRange

    end subroutine s_smear_voidfraction

    !>  The purpose of this subroutine is to obtain the bubble source terms based on Maeda and Colonius (2018)
        !!      and add them to the RHS scalar field.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Conservative variables
        !! @param rhs_vf Time derivative of the conservative variables
    subroutine s_add_rhs_sources(q_cons_vf, q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf

        integer :: i, j, k, l

        if (lag_params%cluster_type >= 4) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        do l = 1, E_idx
                            if (q_beta%vf(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + &
                                                        q_cons_vf(l)%sf(i, j, k)*(q_beta%vf(2)%sf(i, j, k) + &
                                                                                  q_beta%vf(5)%sf(i, j, k))
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
                            if (q_beta%vf(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + &
                                                        q_cons_vf(l)%sf(i, j, k)/q_beta%vf(1)%sf(i, j, k)* &
                                                        q_beta%vf(2)%sf(i, j, k)
                            end if
                        end do
                    end do
                end do
            end do
        end if

        do l = 1, num_dims

            call s_gradient_dir(q_prim_vf(E_idx), q_beta%vf(3), l)

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        if (q_beta%vf(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                            rhs_vf(contxe + l)%sf(i, j, k) = rhs_vf(contxe + l)%sf(i, j, k) - &
                                                             (1._wp - q_beta%vf(1)%sf(i, j, k))/ &
                                                             q_beta%vf(1)%sf(i, j, k)* &
                                                             q_beta%vf(3)%sf(i, j, k)
                        end if
                    end do
                end do
            end do

            !source in energy
            !$acc parallel loop collapse(3) gang vector default(present)
            do k = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(2)%beg, idwbuff(2)%end
                    do i = idwbuff(1)%beg, idwbuff(1)%end
                        q_beta%vf(3)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k)*q_prim_vf(contxe + l)%sf(i, j, k)
                    end do
                end do
            end do

            call s_gradient_dir(q_beta%vf(3), q_beta%vf(4), l)

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        if (q_beta%vf(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                            rhs_vf(E_idx)%sf(i, j, k) = rhs_vf(E_idx)%sf(i, j, k) - &
                                                        q_beta%vf(4)%sf(i, j, k)*(1._wp - q_beta%vf(1)%sf(i, j, k))/ &
                                                        q_beta%vf(1)%sf(i, j, k)
                        end if
                    end do
                end do
            end do
        end do

    end subroutine s_add_rhs_sources

    !> The purpose of this procedure is obtain the bubble driving pressure p_inf
        !! @param bub_id Particle identifier
        !! @param q_prim_vf  Primitive variables
        !! @param ptype 1: p at infinity, 2: averaged P at the bubble location
        !! @param f_pinfl Driving pressure
        !! @param cell Bubble cell
        !! @param Romega Control volume radius
    subroutine s_get_pinf(bub_id, q_prim_vf, ptype, f_pinfl, cell, preterm1, term2, Romega)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_get_pinf
#else
        !$acc routine seq
#endif
        integer, intent(in) :: bub_id, ptype
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), intent(out) :: f_pinfl
        integer, dimension(3), intent(out) :: cell
        real(wp), intent(out), optional :: preterm1, term2, Romega

        real(wp), dimension(3) :: scoord, psi
        real(wp) :: dc, vol, aux
        real(wp) :: volgas, term1, Rbeq, denom
        real(wp) :: charvol, charpres, charvol2, charpres2
        integer, dimension(3) :: cellaux
        integer :: i, j, k
        integer :: mapCells_pinf, smearGrid, smearGridz
        logical :: celloutside

        scoord = mtn_s(bub_id, 1:3, 2)
        f_pinfl = 0._wp

        !< Find current bubble cell
        cell(:) = int(scoord(:))
        !$acc loop seq
        do i = 1, num_dims
            if (scoord(i) < 0._wp) cell(i) = cell(i) - 1
        end do

        if ((lag_params%cluster_type == 1)) then
            !< Getting p_cell in terms of only the current cell by interpolation

            !< Getting the cell volulme as Omega
            if (p > 0) then
                vol = dx(cell(1))*dy(cell(2))*dz(cell(3))
            else
                if (cyl_coord) then
                    vol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
                else
                    vol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
                end if
            end if

            !< Obtain bilinear interpolation coefficients, based on the current location of the bubble.
            psi(1) = (scoord(1) - real(cell(1)))*dx(cell(1)) + x_cb(cell(1) - 1)
            if (cell(1) == (m + buff_size)) then
                cell(1) = cell(1) - 1
                psi(1) = 1._wp
            else if (cell(1) == (-buff_size)) then
                psi(1) = 0._wp
            else
                if (psi(1) < x_cc(cell(1))) cell(1) = cell(1) - 1
                psi(1) = abs((psi(1) - x_cc(cell(1)))/(x_cc(cell(1) + 1) - x_cc(cell(1))))
            end if

            psi(2) = (scoord(2) - real(cell(2)))*dy(cell(2)) + y_cb(cell(2) - 1)
            if (cell(2) == (n + buff_size)) then
                cell(2) = cell(2) - 1
                psi(2) = 1._wp
            else if (cell(2) == (-buff_size)) then
                psi(2) = 0._wp
            else
                if (psi(2) < y_cc(cell(2))) cell(2) = cell(2) - 1
                psi(2) = abs((psi(2) - y_cc(cell(2)))/(y_cc(cell(2) + 1) - y_cc(cell(2))))
            end if

            if (p > 0) then
                psi(3) = (scoord(3) - real(cell(3)))*dz(cell(3)) + z_cb(cell(3) - 1)
                if (cell(3) == (p + buff_size)) then
                    cell(3) = cell(3) - 1
                    psi(3) = 1._wp
                else if (cell(3) == (-buff_size)) then
                    psi(3) = 0._wp
                else
                    if (psi(3) < z_cc(cell(3))) cell(3) = cell(3) - 1
                    psi(3) = abs((psi(3) - z_cc(cell(3)))/(z_cc(cell(3) + 1) - z_cc(cell(3))))
                end if
            else
                psi(3) = 0._wp
            end if

            !< Perform bilinear interpolation
            if (p == 0) then  !2D
                f_pinfl = q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3))*(1._wp - psi(1))*(1._wp - psi(2))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1._wp - psi(2))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3))*(1._wp - psi(1))*psi(2)
            else              !3D
                f_pinfl = q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3))*(1._wp - psi(1))*(1._wp - psi(2))*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1._wp - psi(2))*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3))*(1._wp - psi(1))*psi(2)*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3) + 1)*(1._wp - psi(1))*(1._wp - psi(2))*psi(3)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3) + 1)*psi(1)*(1._wp - psi(2))*psi(3)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2) + 1, cell(3) + 1)*psi(1)*psi(2)*psi(3)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3) + 1)*(1._wp - psi(1))*psi(2)*psi(3)
            end if

            !R_Omega
            dc = (3._wp*vol/(4._wp*pi))**(1._wp/3._wp)

        else if (lag_params%cluster_type >= 2) then
            ! Bubble dynamic closure from Maeda and Colonius (2018)

            ! Range of cells included in Omega
            if (lag_params%smooth_type == 1) then
                mapCells_pinf = mapCells
            else
                stop "lag_params%cluster_type: 2 requires lag_params%smooth_type: 1."
            end if

            ! Include the cell that contains the bubble (mapCells+1+mapCells)
            smearGrid = mapCells_pinf - (-mapCells_pinf) + 1
            smearGridz = smearGrid
            if (p == 0) smearGridz = 1

            charvol = 0._wp
            charpres = 0._wp
            charvol2 = 0._wp
            charpres2 = 0._wp
            vol = 0._wp

            !$acc loop seq
            do i = 1, smearGrid
                !$acc loop seq
                do j = 1, smearGrid
                    !$acc loop seq
                    do k = 1, smearGridz
                        cellaux(1) = cell(1) + i - (mapCells + 1)
                        cellaux(2) = cell(2) + j - (mapCells + 1)
                        cellaux(3) = cell(3) + k - (mapCells + 1)
                        if (p == 0) cellaux(3) = 0

                        !< check if the current cell is outside the computational domain or not (including ghost cells)
                        celloutside = .false.
                        if (num_dims == 2) then
                            if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                                celloutside = .true.
                            end if
                            if (cyl_coord .and. y_cc(cellaux(2)) < 0._wp) then
                                celloutside = .true.
                            end if
                            if ((cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
                                celloutside = .true.
                            end if
                        else
                            if ((cellaux(3) < -buff_size) .or. (cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                                celloutside = .true.
                            end if

                            if ((cellaux(3) > p + buff_size) .or. (cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
                                celloutside = .true.
                            end if
                        end if
                        if (.not. celloutside) then
                            if (cyl_coord .and. (p == 0) .and. (y_cc(cellaux(2)) < 0._wp)) then
                                celloutside = .true.
                            end if
                        end if

                        if (.not. celloutside) then
                            !< Obtaining the cell volulme
                            if (p > 0) then
                                vol = dx(cellaux(1))*dy(cellaux(2))*dz(cellaux(3))
                            else
                                if (cyl_coord) then
                                    vol = dx(cellaux(1))*dy(cellaux(2))*y_cc(cellaux(2))*2._wp*pi
                                else
                                    vol = dx(cellaux(1))*dy(cellaux(2))*lag_params%charwidth
                                end if
                            end if
                            !< Update values
                            charvol = charvol + vol
                            charpres = charpres + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3))*vol
                            charvol2 = charvol2 + vol*q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                            charpres2 = charpres2 + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                        *vol*q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                        end if

                    end do
                end do
            end do

            f_pinfl = charpres2/charvol2
            vol = charvol
            dc = (3._wp*abs(vol)/(4._wp*pi))**(1._wp/3._wp)
        else

            stop "Check cluterflag. Exiting."

        end if

        if (lag_params%pressure_corrector) then

            !Valid if only one bubble exists per cell
            volgas = intfc_rad(bub_id, 2)**3._wp
            denom = intfc_rad(bub_id, 2)**2._wp
            term1 = bub_dphidt(bub_id)*intfc_rad(bub_id, 2)**2._wp
            term2 = intfc_vel(bub_id, 2)*intfc_rad(bub_id, 2)**2._wp

            Rbeq = volgas**(1._wp/3._wp) !surrogate bubble radius
            aux = dc**3._wp - Rbeq**3._wp
            term2 = term2/denom
            term2 = 3._wp/2._wp*term2**2._wp*Rbeq**3._wp*(1._wp - Rbeq/dc)/aux
            preterm1 = 3._wp/2._wp*Rbeq*(dc**2._wp - Rbeq**2._wp)/(aux*denom)

            !Control volume radius
            if (ptype == 2) Romega = dc

            ! Getting p_inf
            if (ptype == 1) then
                f_pinfl = f_pinfl + preterm1*term1 + term2
            end if

        end if

    end subroutine s_get_pinf

    !>  This subroutine updates the Lagrange variables using the tvd RK time steppers.
        !!      The time derivative of the bubble variables must be stored at every stage to avoid precision errors.
        !! @param stage Current tvd RK stage
    subroutine s_update_lagrange_tdv_rk(stage)

        integer, intent(in) :: stage

        integer :: k

        if (time_stepper == 1) then ! 1st order TVD RK
            !$acc parallel loop gang vector default(present) private(k)
            do k = 1, nBubs
                !u{1} = u{n} +  dt * RHS{n}
                intfc_rad(k, 1) = intfc_rad(k, 1) + dt*intfc_draddt(k, 1)
                intfc_vel(k, 1) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                gas_p(k, 1) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                gas_mv(k, 1) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                if (intfc_rad(k, 1) <= 0._wp) stop "Negative bubble radius encountered, please reduce dt"
            end do

            call s_transfer_data_to_tmp()
            call s_write_void_evol(mytime)
            if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()

            if (lag_params%write_bubbles) then
                !$acc update host(gas_p, gas_mv, intfc_rad, intfc_vel)
                call s_write_lag_particles(mytime)
            end if

        elseif (time_stepper == 2) then ! 2nd order TVD RK
            if (stage == 1) then
                !$acc parallel loop gang vector default(present) private(k)
                do k = 1, nBubs
                    !u{1} = u{n} +  dt * RHS{n}
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*intfc_draddt(k, 1)
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                    mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                    mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                    gas_p(k, 2) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                    if (intfc_rad(k, 2) <= 0._wp) stop "Negative bubble radius encountered, please reduce dt"
                end do

            elseif (stage == 2) then
                !$acc parallel loop gang vector default(present) private(k)
                do k = 1, nBubs
                    !u{1} = u{n} + (1/2) * dt * (RHS{n} + RHS{1})
                    intfc_rad(k, 1) = intfc_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/2._wp
                    intfc_vel(k, 1) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/2._wp
                    mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/2._wp
                    mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/2._wp
                    gas_p(k, 1) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/2._wp
                    gas_mv(k, 1) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/2._wp
                    if (intfc_rad(k, 1) <= 0._wp) stop "Negative bubble radius encountered, please reduce dt"
                end do

                call s_transfer_data_to_tmp()
                call s_write_void_evol(mytime)
                if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()

                if (lag_params%write_bubbles) then
                    !$acc update host(gas_p, gas_mv, intfc_rad, intfc_vel)
                    call s_write_lag_particles(mytime)
                end if

            end if

        elseif (time_stepper == 3) then ! 3rd order TVD RK
            if (stage == 1) then
                !$acc parallel loop gang vector default(present) private(k)
                do k = 1, nBubs
                    !u{1} = u{n} +  dt * RHS{n}
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*intfc_draddt(k, 1)
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                    mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                    mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                    gas_p(k, 2) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                    if (intfc_rad(k, 2) <= 0._wp) stop "Negative bubble radius encountered, please reduce dt"
                end do

            elseif (stage == 2) then
                !$acc parallel loop gang vector default(present) private(k)
                do k = 1, nBubs
                    !u{2} = u{n} + (1/4) * dt * [RHS{n} + RHS{1}]
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/4._wp
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/4._wp
                    mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/4._wp
                    mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/4._wp
                    gas_p(k, 2) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/4._wp
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/4._wp
                    if (intfc_rad(k, 2) <= 0._wp) stop "Negative bubble radius encountered, please reduce dt"
                end do
            elseif (stage == 3) then
                !$acc parallel loop gang vector default(present) private(k)
                do k = 1, nBubs
                    !u{n+1} = u{n} + (2/3) * dt * [(1/4)* RHS{n} + (1/4)* RHS{1} + RHS{2}]
                    intfc_rad(k, 1) = intfc_rad(k, 1) + (2._wp/3._wp)*dt*(intfc_draddt(k, 1)/4._wp + intfc_draddt(k, 2)/4._wp + intfc_draddt(k, 3))
                    intfc_vel(k, 1) = intfc_vel(k, 1) + (2._wp/3._wp)*dt*(intfc_dveldt(k, 1)/4._wp + intfc_dveldt(k, 2)/4._wp + intfc_dveldt(k, 3))
                    mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + (2._wp/3._wp)*dt*(mtn_dposdt(k, 1:3, 1)/4._wp + mtn_dposdt(k, 1:3, 2)/4._wp + mtn_dposdt(k, 1:3, 3))
                    mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + (2._wp/3._wp)*dt*(mtn_dveldt(k, 1:3, 1)/4._wp + mtn_dveldt(k, 1:3, 2)/4._wp + mtn_dveldt(k, 1:3, 3))
                    gas_p(k, 1) = gas_p(k, 1) + (2._wp/3._wp)*dt*(gas_dpdt(k, 1)/4._wp + gas_dpdt(k, 2)/4._wp + gas_dpdt(k, 3))
                    gas_mv(k, 1) = gas_mv(k, 1) + (2._wp/3._wp)*dt*(gas_dmvdt(k, 1)/4._wp + gas_dmvdt(k, 2)/4._wp + gas_dmvdt(k, 3))
                    if (intfc_rad(k, 1) <= 0._wp) stop "Negative bubble radius encountered, please reduce dt"
                end do

                call s_transfer_data_to_tmp()
                call s_write_void_evol(mytime)
                if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()

                if (lag_params%write_bubbles) then
                    !$acc update host(gas_p, gas_mv, intfc_rad, intfc_vel)
                    call s_write_lag_particles(mytime)
                end if

            end if

        end if

    end subroutine s_update_lagrange_tdv_rk

    !>  This subroutine updates the Euler-Lagrange temporal variables before entering to the next time-stage in the RKCK stepper.
        !! @param RKstep Current time step in the RKCK adaptive stepper
        !! @param q_cons_ts Conservative variables
        !! @param rhs_ts Time derivatives of the conservative variables
        !! @param lag_largestep Negative radius flag
    subroutine s_update_tmp_rkck(RKstep, q_cons_ts, rhs_ts, lag_largestep)

        integer, intent(in) :: RKstep
        type(vector_field), dimension(:), intent(inout) :: q_cons_ts
        type(vector_field), dimension(:), intent(inout) :: rhs_ts
        real(wp), intent(out) :: lag_largestep

        integer :: i, j, k, l, q
        real(wp) :: radiusOld, velOld, aux_glb
        integer :: remove_id

        call s_transfer_data_to_tmp()

        lag_largestep = 0._wp
        remove_id = 0
        !$acc parallel loop gang vector default(present) reduction(+:lag_largestep) &
        !$acc reduction(MAX: remove_id) copyin(RKstep) copy(lag_largestep, remove_id)
        do k = 1, nBubs

            radiusOld = intfc_rad(k, 2)
            velOld = intfc_vel(k, 2)

            !$acc loop seq
            do i = 1, RKstep
                intfc_rad(k, 2) = intfc_rad(k, 2) + dt*lag_RKCKcoef(RKstep, i)*intfc_draddt(k, i)
                intfc_vel(k, 2) = intfc_vel(k, 2) + dt*lag_RKCKcoef(RKstep, i)*intfc_dveldt(k, i)
                mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 2) + dt*lag_RKCKcoef(RKstep, i)*mtn_dposdt(k, 1:3, i)
                mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 2) + dt*lag_RKCKcoef(RKstep, i)*mtn_dveldt(k, 1:3, i)
                gas_p(k, 2) = gas_p(k, 2) + dt*lag_RKCKcoef(RKstep, i)*gas_dpdt(k, i)
                gas_mv(k, 2) = gas_mv(k, 2) + dt*lag_RKCKcoef(RKstep, i)*gas_dmvdt(k, i)
            end do

            if ((intfc_rad(k, 2) <= 0._wp) .or. &                       ! no negative radius
                (mtn_pos(k, 1, 2) /= mtn_pos(k, 1, 2))) then  ! finite bubble location
                print *, 'Negative bubble radius encountered'
                lag_largestep = lag_largestep + 1._wp
                if (dt < 2._wp*verysmall_dt) then
                    remove_id = max(remove_id, k)
                end if
            end if

        end do

        if (remove_id /= 0) call s_remove_lag_bubble(remove_id)

#ifdef MFC_MPI
        if (num_procs > 1) then
            call s_mpi_allreduce_sum(lag_largestep, aux_glb)
            lag_largestep = aux_glb
        end if
#endif

        if (lag_largestep > 0._wp) return

        ! Update background fluid variables
        !$acc parallel loop collapse(4) gang vector default(present) copyin(RKstep)
        do l = 1, sys_size
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        q_cons_ts(2)%vf(l)%sf(i, j, k) = &
                            q_cons_ts(1)%vf(l)%sf(i, j, k)
                        !$acc loop seq
                        do q = 1, RKstep
                            q_cons_ts(2)%vf(l)%sf(i, j, k) = &
                                q_cons_ts(2)%vf(l)%sf(i, j, k) + &
                                dt*lag_RKCKcoef(RKstep, q)*rhs_ts(q)%vf(l)%sf(i, j, k)
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
        !! @param rkck_errmax Truncation error
    subroutine s_calculate_rkck_truncation_error(rkck_errmax)

        real(wp), intent(out) :: rkck_errmax

        real(wp) :: erraux, errb
        integer :: i, j, k

        rkck_errmax = 0._wp
        !$acc parallel loop gang vector default(present) reduction(MAX: rkck_errmax) copy(rkck_errmax)
        do k = 1, nBubs
            errb = 0._wp

            !Bubble radius error
            erraux = 0._wp
            !$acc loop seq
            do i = 1, lag_num_ts
                erraux = erraux + lag_RKCKcoef(7, i)*intfc_draddt(k, i)
            end do
            errb = max(errb, abs(erraux)*dt/bub_R0(k))

            !Interface velocity error
            erraux = 0._wp
            !$acc loop seq
            do i = 1, lag_num_ts
                erraux = erraux + lag_RKCKcoef(7, i)*intfc_dveldt(k, i)
            end do
            errb = max(errb, abs(erraux)*dt)

            !Bubble velocity error
            !$acc loop seq
            do j = 1, 3
                erraux = 0._wp
                !$acc loop seq
                do i = 1, lag_num_ts
                    erraux = erraux + lag_RKCKcoef(7, i)*mtn_dposdt(k, j, i)
                end do
                errb = max(errb, abs(erraux)*dt/(abs(mtn_vel(k, j, 2)) + 1.0d-4))
            end do
            rkck_errmax = max(rkck_errmax, errb)
        end do

    end subroutine s_calculate_rkck_truncation_error

    !>  This subroutine updates the conservative fields and the lagrangian variables after accepting the performed time step.
        !! @param q_cons_ts Conservative variables
    subroutine s_update_rkck(q_cons_ts)

        type(vector_field), dimension(:), intent(inout) :: q_cons_ts

        integer :: i, j, k, l

        !$acc parallel loop gang vector default(present) private(k)
        do k = 1, nBubs
            !Accept time step (actual vars = temporal vars)
            mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 2)
            mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 2)
            intfc_rad(k, 1) = intfc_rad(k, 2)
            intfc_vel(k, 1) = intfc_vel(k, 2)
            gas_p(k, 1) = gas_p(k, 2)
            gas_mv(k, 1) = gas_mv(k, 2)
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

    end subroutine s_update_rkck

    !>  This subroutine computes the next time step in the adaptive RKCK stepper in the CPU.
        !! @param lag_largestep Negative radius flag
        !! @param restart_rkck_step Restart the current time step
        !! @param rkck_errmax Truncation error
    subroutine s_compute_rkck_dt(lag_largestep, restart_rkck_step, rkck_errmax)

        real(wp), intent(in) :: lag_largestep
        logical, intent(out) :: restart_rkck_step
        real(wp), intent(inout), optional :: rkck_errmax

        real(wp) :: htemp, aux_glb

        restart_rkck_step = .false.

        if (lag_largestep > 0._wp) then ! Encountered negative radius, so reduce dt and restart time step

            if (rkck_adap_dt) then
                if (dt > verysmall_dt) then
                    restart_rkck_step = .true.
                    dt = SHRNKDT*dt
                    if (num_procs > 1) then
                        call s_mpi_allreduce_min(dt, aux_glb)
                        dt = aux_glb
                    end if
                    !$acc update device(dt)
                    if (proc_rank == 0) print *, '>>>>> WARNING: Reducing dt and restarting time step, now dt: ', dt
                else
                    call s_mpi_abort('Time step smaller than 1e-14')
                end if
            else
                call s_mpi_abort('Time step too large, please reduce dt or enable rkck_adap_dt')
            end if

            return

        end if

        if (rkck_adap_dt) then ! Checking truncation error

            rkck_errmax = min(rkck_errmax, 1._wp)
            if (num_procs > 1) then
                call s_mpi_allreduce_max(rkck_errmax, aux_glb)
                rkck_errmax = aux_glb
            end if
            rkck_errmax = rkck_errmax/rkck_tolerance ! Scale relative to user required tolerance.

            if ((rkck_errmax > 1._wp)) then   ! Truncation error too large, reduce dt and restart time step
                restart_rkck_step = .true.
                htemp = SAFETY*dt*(rkck_errmax**PSHRNK)
                dt = sign(max(abs(htemp), (1._wp - SAFETY)*abs(dt)), dt)  ! No more than a factor of 10.
            else                            ! Step succeeded. Compute size of next step.
                if (rkck_errmax > ERRCON) then
                    dt = SAFETY*dt*(rkck_errmax**PGROW) ! No more than a factor of 5 increase.
                else
                    dt = (1._wp/SHRNKDT)*dt ! Truncation error too small (< 1.89e-4), increase time step
                end if
            end if

            !dt precision accuracy is 16 digits
            dt = (ceiling(dt*RNDDEC)*RNDDEC + ceiling(dt*(RNDDEC**2._wp) - ceiling(dt*RNDDEC)*RNDDEC))/(RNDDEC**2._wp)
            dt = min(dt, dt_max)

            if (dt < 0._wp) call s_mpi_abort('dt must not be negative')
            if (num_procs > 1) then
                call s_mpi_allreduce_min(dt, aux_glb)
                dt = aux_glb
            end if
            !$acc update device(dt)

            if (restart_rkck_step) then
                if (proc_rank == 0) print '("WARNING: Truncation error found. Restaring time step, and now dt = "ES8.6"")', &
                    dt
            end if

        end if

    end subroutine s_compute_rkck_dt

    !> This subroutine returns the computational coordinate of the cell for the given position.
          !! @param pos Input coordinates
          !! @param cell Computational coordinate of the cell
          !! @param scoord Calculated particle coordinates
    subroutine s_locate_cell(pos, cell, scoord)

        real(wp), dimension(3), intent(in) :: pos
        real(wp), dimension(3), intent(out) :: scoord
        integer, dimension(3), intent(inout) :: cell

        integer :: i

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
        scoord(1) = cell(1) + (pos(1) - x_cb(cell(1) - 1))/dx(cell(1))
        scoord(2) = cell(2) + (pos(2) - y_cb(cell(2) - 1))/dy(cell(2))
        scoord(3) = 0._wp
        if (p > 0) scoord(3) = cell(3) + (pos(3) - z_cb(cell(3) - 1))/dz(cell(3))
        cell(:) = int(scoord(:))
        do i = 1, num_dims
            if (scoord(i) < 0._wp) cell(i) = cell(i) - 1
        end do

    end subroutine s_locate_cell

    !> This subroutine transfer data into the temporal variables.
    subroutine s_transfer_data_to_tmp()

        integer :: k

        !$acc parallel loop gang vector default(present) private(k)
        do k = 1, nBubs
            gas_p(k, 2) = gas_p(k, 1)
            gas_mv(k, 2) = gas_mv(k, 1)
            intfc_rad(k, 2) = intfc_rad(k, 1)
            intfc_vel(k, 2) = intfc_vel(k, 1)
            mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1)
            mtn_posPrev(k, 1:3, 2) = mtn_posPrev(k, 1:3, 1)
            mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1)
            mtn_s(k, 1:3, 2) = mtn_s(k, 1:3, 1)
        end do

    end subroutine s_transfer_data_to_tmp

    !> The purpose of this procedure is to determine if the global coordinates of the bubbles
        !!      are present in the current MPI processor (including ghost cells).
        !! @param pos_part Spatial coordinates of the bubble
    function particle_in_domain(pos_part)

        logical :: particle_in_domain
        real(wp), dimension(3), intent(in) :: pos_part

        ! 2D
        if (p == 0 .and. cyl_coord .neqv. .true.) then
            ! Defining a virtual z-axis that has the same dimensions as y-axis
            ! defined in the input file
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
                                  (pos_part(3) < lag_params%charwidth/2._wp) .and. (pos_part(3) >= -lag_params%charwidth/2._wp))
        else
            ! cyl_coord
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (abs(pos_part(2)) < y_cb(n + buff_size)) .and. (abs(pos_part(2)) >= max(y_cb(-buff_size - 1), 0._wp)))
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
        real(wp), dimension(3), intent(in) :: pos_part

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

        integer :: i, j, k

        if (dir == 1) then
            ! Gradient in x dir.
            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        dq%sf(i, j, k) = q%sf(i, j, k)*(dx(i + 1) - dx(i - 1)) &
                                         + q%sf(i + 1, j, k)*(dx(i) + dx(i - 1)) &
                                         - q%sf(i - 1, j, k)*(dx(i) + dx(i + 1))
                        dq%sf(i, j, k) = dq%sf(i, j, k)/ &
                                         ((dx(i) + dx(i - 1))*(dx(i) + dx(i + 1)))
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
                            dq%sf(i, j, k) = q%sf(i, j, k)*(dy(j + 1) - dy(j - 1)) &
                                             + q%sf(i, j + 1, k)*(dy(j) + dy(j - 1)) &
                                             - q%sf(i, j - 1, k)*(dy(j) + dy(j + 1))
                            dq%sf(i, j, k) = dq%sf(i, j, k)/ &
                                             ((dy(j) + dy(j - 1))*(dy(j) + dy(j + 1)))
                        end do
                    end do
                end do
            else
                ! Gradient in z dir.
                !$acc parallel loop collapse(3) gang vector default(present)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            dq%sf(i, j, k) = q%sf(i, j, k)*(dz(k + 1) - dz(k - 1)) &
                                             + q%sf(i, j, k + 1)*(dz(k) + dz(k - 1)) &
                                             - q%sf(i, j, k - 1)*(dz(k) + dz(k + 1))
                            dq%sf(i, j, k) = dq%sf(i, j, k)/ &
                                             ((dz(k) + dz(k - 1))*(dz(k) + dz(k + 1)))
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_gradient_dir

    !> Subroutine that writes on each time step the changes of the lagrangian bubbles.
        !!  @param q_time Current time
    subroutine s_write_lag_particles(qtime)

        real(wp), intent(in) :: qtime
        integer :: k

        character(LEN=path_len + 2*name_len) :: file_loc

        write (file_loc, '(A,I0,A)') 'lag_bubble_evol_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)

        if (qtime == 0._wp) then
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
                lag_id(k, 1), &
                mtn_pos(k, 1, 1), &
                mtn_pos(k, 2, 1), &
                mtn_pos(k, 3, 1), &
                gas_mv(k, 1), &
                gas_mv(k, 1)/(gas_mv(k, 1) + gas_mg(k)), &
                intfc_rad(k, 1), &
                intfc_vel(k, 1), &
                gas_p(k, 1)
        end do

        close (11)

    end subroutine s_write_lag_particles

    !>  Subroutine that writes some useful statistics related to the volume fraction
            !!       of the particles (void fraction) in the computatioational domain
            !!       on each time step.
            !!  @param q_time Current time
    subroutine s_write_void_evol(qtime)

        real(wp), intent(in) :: qtime
        real(wp) :: volcell, voltot
        real(wp) :: lag_void_max, lag_void_avg, lag_vol
        real(wp) :: void_max_glb, void_avg_glb, vol_glb

        integer :: i, j, k
        integer, dimension(3) :: cell

        character(LEN=path_len + 2*name_len) :: file_loc

        if (proc_rank == 0) then
            write (file_loc, '(A)') 'voidfraction.dat'
            file_loc = trim(case_dir)//'/D/'//trim(file_loc)
            if (qtime == 0._wp) then
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

        lag_void_max = 0._wp
        lag_void_avg = 0._wp
        lag_vol = 0._wp
        !$acc parallel loop collapse(3) gang vector default(present) reduction(+:lag_vol,lag_void_avg) &
        !$acc reduction(MAX:lag_void_max) copy(lag_vol, lag_void_avg, lag_void_max)
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    lag_void_max = max(lag_void_max, 1._wp - q_beta%vf(1)%sf(i, j, k))
                    call s_get_char_vol(i, j, k, volcell)
                    if ((1._wp - q_beta%vf(1)%sf(i, j, k)) > 5.0d-11) then
                        lag_void_avg = lag_void_avg + (1._wp - q_beta%vf(1)%sf(i, j, k))*volcell
                        lag_vol = lag_vol + volcell
                    end if
                end do
            end do
        end do

#ifdef MFC_MPI
        if (num_procs > 1) then
            call s_mpi_allreduce_max(lag_void_max, void_max_glb)
            lag_void_max = void_max_glb
            call s_mpi_allreduce_sum(lag_vol, vol_glb)
            lag_vol = vol_glb
            call s_mpi_allreduce_sum(lag_void_avg, void_avg_glb)
            lag_void_avg = void_avg_glb
        end if
#endif
        voltot = lag_void_avg
        ! This voidavg value does not reflect the real void fraction in the cloud
        ! since the cell which does not have bubbles are not accounted
        if (lag_vol > 0._wp) lag_void_avg = lag_void_avg/lag_vol

        if (proc_rank == 0) then
            write (12, '(6X,4e24.8)') &
                qtime, &
                lag_void_avg, &
                lag_void_max, &
                voltot
            close (12)
        end if

    end subroutine s_write_void_evol

    !>  Subroutine that writes the restarting files for the particles in the lagrangian solver.
        !!  @param t_step Current time step
    subroutine s_write_restart_lag_bubbles(t_step)

        ! Generic string used to store the address of a particular file
        integer, intent(in) :: t_step

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist
        integer :: bub_id, tot_part, tot_part_wrtn, npart_wrtn
        integer :: i, k

#ifdef MFC_MPI
        ! For Parallel I/O
        integer :: ifile, ierr
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer :: view
        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer, dimension(num_procs) :: part_order, part_ord_mpi

        bub_id = 0._wp
        if (nBubs /= 0) then
            do k = 1, nBubs
                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then
                    bub_id = bub_id + 1
                end if
            end do
        end if

        if (.not. parallel_io) return

        ! Total number of particles
        call MPI_ALLREDUCE(bub_id, tot_part, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        ! Total number of particles written so far
        call MPI_ALLREDUCE(npart_wrtn, tot_part_wrtn, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        lsizes(1) = max(1, bub_id)
        lsizes(2) = 21

        ! if the partcle number is zero, put 1 since MPI cannot deal with writing
        ! zero particle
        part_order(:) = 1
        part_order(proc_rank + 1) = max(1, bub_id)

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
                                      MPI_ORDER_FORTRAN, mpi_p, view, ierr)
        call MPI_type_COMMIT(view, ierr)

        allocate (MPI_IO_DATA_lag_bubbles(1:max(1, bub_id), 1:21))

        ! Open the file to write all flow variables
        write (file_loc, '(A,I0,A)') 'lag_bubbles_', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) then
            call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        end if

        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                           mpi_info_int, ifile, ierr)

        disp = 0._wp

        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, &
                               'native', mpi_info_null, ierr)

        ! Cycle through list
        i = 1

        if (bub_id == 0) then
            MPI_IO_DATA_lag_bubbles(1, 1:21) = 0._wp
        else

            do k = 1, nBubs

                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then

                    MPI_IO_DATA_lag_bubbles(i, 1) = real(lag_id(k, 1))
                    MPI_IO_DATA_lag_bubbles(i, 2:4) = mtn_pos(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 5:7) = mtn_posPrev(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 8:10) = mtn_vel(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 11) = intfc_rad(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 12) = intfc_vel(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 13) = bub_R0(k)
                    MPI_IO_DATA_lag_bubbles(i, 14) = Rmax_stats(k)
                    MPI_IO_DATA_lag_bubbles(i, 15) = Rmin_stats(k)
                    MPI_IO_DATA_lag_bubbles(i, 16) = bub_dphidt(k)
                    MPI_IO_DATA_lag_bubbles(i, 17) = gas_p(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 18) = gas_mv(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 19) = gas_mg(k)
                    MPI_IO_DATA_lag_bubbles(i, 20) = gas_betaT(k)
                    MPI_IO_DATA_lag_bubbles(i, 21) = gas_betaC(k)

                    i = i + 1

                end if

            end do

        end if

        call MPI_FILE_write_ALL(ifile, MPI_IO_DATA_lag_bubbles, 21*max(1, bub_id), &
                                mpi_p, status, ierr)

        call MPI_FILE_CLOSE(ifile, ierr)

        deallocate (MPI_IO_DATA_lag_bubbles)

#endif

    end subroutine s_write_restart_lag_bubbles

    !>  This procedure calculates the maximum and minimum radius of each bubble.
    subroutine s_calculate_lag_bubble_stats()

        integer :: k

        !$acc parallel loop gang vector default(present) reduction(MAX:Rmax_glb) &
        !$acc reduction(MIN: Rmin_glb) copy(Rmax_glb, Rmin_glb)
        do k = 1, nBubs
            Rmax_glb = max(Rmax_glb, intfc_rad(k, 1)/bub_R0(k))
            Rmin_glb = min(Rmin_glb, intfc_rad(k, 1)/bub_R0(k))
            Rmax_stats(k) = max(Rmax_stats(k), intfc_rad(k, 1)/bub_R0(k))
            Rmin_stats(k) = min(Rmin_stats(k), intfc_rad(k, 1)/bub_R0(k))
        end do

    end subroutine s_calculate_lag_bubble_stats

    !>  Subroutine that writes the maximum and minimum radius of each bubble.
    subroutine s_write_lag_bubble_stats()

        integer :: k
        character(LEN=path_len + 2*name_len) :: file_loc

        write (file_loc, '(A,I0,A)') 'stats_lag_bubbles_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)

        !$acc update host(Rmax_glb, Rmin_glb)

        open (13, FILE=trim(file_loc), FORM='formatted', position='rewind')
        write (13, *) 'proc_rank, particleID, x, y, z, Rmax_glb, Rmin_glb'

        do k = 1, nBubs
            write (13, '(6X,2I24.8,5e24.8)') &
                proc_rank, &
                lag_id(k, 1), &
                mtn_pos(k, 1, 1), &
                mtn_pos(k, 2, 1), &
                mtn_pos(k, 3, 1), &
                Rmax_stats(k), &
                Rmin_stats(k)
        end do

        close (13)

    end subroutine s_write_lag_bubble_stats

    !> The purpose of this subroutine is to remove one specific particle if dt is too small.
          !! @param bub_id Particle id
    subroutine s_remove_lag_bubble(bub_id)

        integer, intent(in) :: bub_id

        integer :: i

        !$acc loop seq
        do i = bub_id, nBubs - 1
            lag_id(i, 1) = lag_id(i + 1, 1)
            bub_R0(i) = bub_R0(i + 1)
            Rmax_stats(i) = Rmax_stats(i + 1)
            Rmin_stats(i) = Rmin_stats(i + 1)
            gas_mg(i) = gas_mg(i + 1)
            gas_betaT(i) = gas_betaT(i + 1)
            gas_betaC(i) = gas_betaC(i + 1)
            bub_dphidt(i) = bub_dphidt(i + 1)
            gas_p(i, 1:2) = gas_p(i + 1, 1:2)
            gas_mv(i, 1:2) = gas_mv(i + 1, 1:2)
            intfc_rad(i, 1:2) = intfc_rad(i + 1, 1:2)
            intfc_vel(i, 1:2) = intfc_vel(i + 1, 1:2)
            mtn_pos(i, 1:3, 1:2) = mtn_pos(i + 1, 1:3, 1:2)
            mtn_posPrev(i, 1:3, 1:2) = mtn_posPrev(i + 1, 1:3, 1:2)
            mtn_vel(i, 1:3, 1:2) = mtn_vel(i + 1, 1:3, 1:2)
            mtn_s(i, 1:3, 1:2) = mtn_s(i + 1, 1:3, 1:2)
            intfc_draddt(i, 1:lag_num_ts) = intfc_draddt(i + 1, 1:lag_num_ts)
            intfc_dveldt(i, 1:lag_num_ts) = intfc_dveldt(i + 1, 1:lag_num_ts)
            gas_dpdt(i, 1:lag_num_ts) = gas_dpdt(i + 1, 1:lag_num_ts)
            gas_dmvdt(i, 1:lag_num_ts) = gas_dmvdt(i + 1, 1:lag_num_ts)
        end do

        nBubs = nBubs - 1
        !$acc update device(nBubs)

    end subroutine s_remove_lag_bubble

    !> The purpose of this subroutine is to deallocate variables
    subroutine s_finalize_lagrangian_solver()

        integer :: i

        do i = 1, q_beta_idx
            @:DEALLOCATE(q_beta%vf(i)%sf)
        end do
        @:DEALLOCATE(q_beta%vf)

        !Deallocating space
        if (time_stepper == 4) then
            @:DEALLOCATE(lag_RKCKcoef)
        end if
        @:DEALLOCATE(lag_id)
        @:DEALLOCATE(bub_R0)
        @:DEALLOCATE(Rmax_stats)
        @:DEALLOCATE(Rmin_stats)
        @:DEALLOCATE(gas_mg)
        @:DEALLOCATE(gas_betaT)
        @:DEALLOCATE(gas_betaC)
        @:DEALLOCATE(bub_dphidt)
        @:DEALLOCATE(gas_p)
        @:DEALLOCATE(gas_mv)
        @:DEALLOCATE(intfc_rad)
        @:DEALLOCATE(intfc_vel)
        @:DEALLOCATE(mtn_pos)
        @:DEALLOCATE(mtn_posPrev)
        @:DEALLOCATE(mtn_vel)
        @:DEALLOCATE(mtn_s)
        @:DEALLOCATE(intfc_draddt)
        @:DEALLOCATE(intfc_dveldt)
        @:DEALLOCATE(gas_dpdt)
        @:DEALLOCATE(gas_dmvdt)
        @:DEALLOCATE(mtn_dposdt)
        @:DEALLOCATE(mtn_dveldt)

    end subroutine s_finalize_lagrangian_solver

end module m_bubbles_EL
