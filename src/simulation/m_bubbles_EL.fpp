!>
!! @file m_bubbles_EL.fpp
!! @brief Contains module m_bubbles_EL

#:include 'macros.fpp'

!> @brief This module is used to add the lagrangian subgrid bubble model
module m_bubbles_EL

    ! Dependencies =============================================================

    use m_global_parameters             !< Definitions of the global parameters

    use m_derived_types                 !< Definitions of the derived types

    use m_rhs                           !< Right-hand-side (RHS) evaluation procedures

    use m_data_output                   !< Run-time info & solution data output procedures

    use m_mpi_proxy                     !< Message passing interface (MPI) module proxy

    use m_bubbles_EL_kernels            !< Definitions of the kernel functions

    use m_variables_conversion

    use m_compile_specific

    use m_boundary_conditions

    ! ==========================================================================

    implicit none

    real(kind(0d0)), allocatable, dimension(:, :) :: lag_RKCKcoef    !< RKCK 4th-5th time stepper coefficients
    !$acc declare create(lag_RKCKcoef)

    integer, allocatable, dimension(:, :) :: lag_id  !< Global and local IDs
    !$acc declare create(lag_id)

    !(nBub)
    real(kind(0d0)), allocatable, dimension(:) :: bub_R0    !< Initial bubble radius
    real(kind(0d0)), allocatable, dimension(:) :: Rmax_stats  !< Maximum radius
    real(kind(0d0)), allocatable, dimension(:) :: Rmin_stats  !< Minimum radius
    real(kind(0d0)), allocatable, dimension(:) :: gas_mg        !< Bubble's gas mass
    real(kind(0d0)), allocatable, dimension(:) :: gas_betaT     !< heatflux model (Preston et al., 2007)
    real(kind(0d0)), allocatable, dimension(:) :: gas_betaC     !< massflux model (Preston et al., 2007)
    real(kind(0d0)), allocatable, dimension(:) :: bub_dphidt        !< subgrid velocity potential (Maeda & Colonius, 2018)
    !$acc declare create(bub_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, gas_betaC, bub_dphidt)

    !(nBub, 1 -> actual val or 2 -> temp val)
    real(kind(0d0)), allocatable, dimension(:, :) :: gas_p       !< Pressure in the bubble
    real(kind(0d0)), allocatable, dimension(:, :) :: gas_mv      !< Vapor mass in the bubble
    real(kind(0d0)), allocatable, dimension(:, :) :: intfc_rad   !< Bubble radius
    real(kind(0d0)), allocatable, dimension(:, :) :: intfc_vel   !< Velocity of the bubble interface
    !$acc declare create(gas_p, gas_mv, intfc_rad, intfc_vel)

    !(nBub, 1-> x or 2->y or 3 ->z, 1 -> actual or 2 -> temporal val)
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mtn_pos        !< Bubble's position
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mtn_posPrev    !< Bubble's previous position
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mtn_vel        !< Bubble's velocity
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mtn_s          !< Bubble's computational cell position in real format
    !$acc declare create(mtn_pos, mtn_posPrev, mtn_vel, mtn_s)

    real(kind(0d0)), allocatable, dimension(:, :) :: intfc_draddt    !< Time derivative of bubble's radius
    real(kind(0d0)), allocatable, dimension(:, :) :: intfc_dveldt    !< Time derivative of bubble's interface velocity
    real(kind(0d0)), allocatable, dimension(:, :) :: gas_dpdt            !< Time derivative of gas pressure
    real(kind(0d0)), allocatable, dimension(:, :) :: gas_dmvdt           !< Time derivative of the vapor mass in the bubble
    !$acc declare create(intfc_draddt, intfc_dveldt, gas_dpdt, gas_dmvdt)
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mtn_dposdt     !< Time derivative of the bubble's position
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mtn_dveldt     !< Time derivative of the bubble's velocity
    !$acc declare create(mtn_dposdt, mtn_dveldt)

    type(vector_field) :: q_beta    !< Projection of the lagrangian particles in the Eulerian framework
    integer :: q_beta_idx           !< Size of the q_beta vector field
    !$acc declare create(q_beta, q_beta_idx)

    integer :: nBubs                !< Number of bubbles in the local domain
    !$acc declare create(nBubs)

    real(kind(0.d0)) :: Rmax_glb, Rmin_glb      !< Maximum and minimum bubbe size in the local domain
    !$acc declare create(Rmax_glb, Rmin_glb)

    integer, private :: lag_num_ts              !<  Number of time stages in the time-stepping scheme
    !$acc declare create(lag_num_ts)

contains

    !> Initializes the lagrangian subgrid bubble solver
        !! @param q_cons_vf Initial conservative variables
    subroutine s_initialize_lag_bubbles(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        integer :: nBubs_glb, i

        ! Setting number of time-stages for selected time-stepping scheme
        lag_num_ts = time_stepper
        if (time_stepper == 4) lag_num_ts = num_ts_rkck

        ! Allocate space for the Eulerian fields needed to map the effect of the bubbles
        if (lag_params%solver_approach == 1) then ! One-way coupling
            q_beta_idx = 3
        elseif (lag_params%solver_approach == 2) then ! Two-way coupling
            q_beta_idx = 4
            if (lag_params%cluster_type >= 4) q_beta_idx = 6 !Subgrid noise model for 2D approximation
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
                lag_RKCKcoef(1, i) = RKCKcoef1(i)
                lag_RKCKcoef(2, i) = RKCKcoef2(i)
                lag_RKCKcoef(3, i) = RKCKcoef3(i)
                lag_RKCKcoef(4, i) = RKCKcoef4(i)
                lag_RKCKcoef(5, i) = RKCKcoef5(i)
                lag_RKCKcoef(6, i) = RKCKcoef6(i)
                lag_RKCKcoef(7, i) = RKCKcoefE(i)
            end do
        end if

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
        integer :: i, j, k, l, ios
        logical :: file_exist, indomain

        character(LEN=path_len + 2*name_len) :: path_D_dir !<

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
            if (proc_rank == 0) print *, 'Reading lagrange bubbles input file.'
            inquire (file='input/lag_bubbles.dat', exist=file_exist)
            if (file_exist) then
                open (94, file='input/lag_bubbles.dat', form='formatted', iostat=ios)
                do while (ios == 0)
                    read (94, *, iostat=ios) (inputparticle(i), i=1, 8)
                    if (ios /= 0) cycle
                    indomain = particle_in_domain(inputparticle(1:3))
                    id = id + 1
                    if (id > lag_params%nBubs_glb .and. proc_rank == 0) call s_mpi_abort('Current number of bubbles is larger than nBubs_glb')
                    if (indomain) then
                        nparticles = nparticles + 1
                        call s_add_lag_bubble(inputparticle, q_cons_vf, nparticles)
                        lag_id(nparticles, 1) = id  !global ID
                        lag_id(nparticles, 2) = nparticles  !local ID
                        nBubs = nparticles ! local number of bubbles
                    end if
                end do
                close (94)
            else
                stop "if you include lagrange bubbles, you have to initialize them in input/lag_bubbles.dat"
            end if
        else
            if (proc_rank == 0) print *, 'Restarting lagrange bubbles at save_count: ', save_count
            call s_add_lag_bubble_restart(nparticles, save_count)
        end if

        print *, " Lagrange bubbles running, in proc", proc_rank, "number:", nparticles, "/", id

        !$acc update device(bubbles_lagrange, lag_params)

        !$acc update device(lag_id, bub_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, gas_betaC,    &
        !$acc bub_dphidt, gas_p, gas_mv, intfc_rad, intfc_vel, mtn_pos, mtn_posPrev, &
        !$acc mtn_vel, mtn_s, intfc_draddt, intfc_dveldt, gas_dpdt, gas_dmvdt,       &
        !$acc mtn_dposdt, mtn_dveldt, nBubs)

        Rmax_glb = min(dflt_real, -dflt_real)
        Rmin_glb = max(dflt_real, -dflt_real)
        !$acc update device(Rmax_glb, Rmin_glb)

        if (time_stepper == 4) then
            dt_max = dt     !Initial and largest dt - rkck stepper
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
        if (lag_params%massTransfer_model) massflag = 1.0d0
        if (lag_params%heatTransfer_model) heatflag = 1.0d0

        bub_R0(nparticles) = inputparticle(7)
        Rmax_stats(nparticles) = min(dflt_real, -dflt_real)
        Rmin_stats(nparticles) = max(dflt_real, -dflt_real)
        bub_dphidt(nparticles) = 0.0d0
        intfc_rad(nparticles, 1) = inputparticle(7)
        intfc_vel(nparticles, 1) = inputparticle(8)
        mtn_pos(nparticles, 1:3, 1) = inputparticle(1:3)
        mtn_posPrev(nparticles, 1:3, 1) = mtn_pos(nparticles, 1:3, 1)
        mtn_vel(nparticles, 1:3, 1) = inputparticle(4:6)

        if (cyl_coord .and. p == 0) then
            mtn_pos(nparticles, 2, 1) = dsqrt(mtn_pos(nparticles, 2, 1)**2d0 + &
                                              mtn_pos(nparticles, 3, 1)**2d0)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            mtn_pos(nparticles, 3, 1) = atan2(inputparticle(3), inputparticle(2))
            mtn_posPrev(nparticles, 1:3, 1) = mtn_pos(nparticles, 1:3, 1)
        end if

        cell = -buff_size
        call s_locate_cell(mtn_pos(nparticles, 1:3, 1), cell, mtn_s(nparticles, 1:3, 1))

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
        gas_p(nparticles, 1) = pliq + 2.0*lag_params%sigmabubble/bub_R0(nparticles)
        if (lag_params%sigmabubble /= 0.0d0) then
            pcrit = lag_params%pvap - 4.0d0*lag_params%sigmabubble/(3.d0*sqrt(3.0d0*gas_p(nparticles, 1)*bub_R0(nparticles)**3/(2.0d0*lag_params%sigmabubble)))
            pref = gas_p(nparticles, 1)
        else
            pcrit = 0.0d0
        end if

        ! Initial particle mass
        volparticle = 4.0d0/3.0d0*pi*bub_R0(nparticles)**3 ! volume
        gas_mv(nparticles, 1) = lag_params%pvap*volparticle*(1.0d0/(lag_params%Rvapor*lag_params%Thost))*(massflag) ! vapermass
        gas_mg(nparticles) = (gas_p(nparticles, 1) - lag_params%pvap*(massflag))*volparticle*(1.0d0/(lag_params%Rgas*lag_params%Thost)) ! gasmass
        if (gas_mg(nparticles) <= 0.0d0) stop 'the initial mass of gas inside the bubble is negative. Check your initial conditions'
        totalmass = gas_mg(nparticles) + gas_mv(nparticles, 1) ! totalmass

        ! Bubble natural frequency
        concvap = gas_mv(nparticles, 1)/(gas_mv(nparticles, 1) + gas_mg(nparticles))
        omegaN = (3.0d0*(gas_p(nparticles, 1) - lag_params%pvap*(massflag)) + 4.0d0*lag_params%sigmabubble/bub_R0(nparticles))/rhol
        if (lag_params%pvap*(massflag) > gas_p(nparticles, 1)) then
            print *, 'Not allowed: bubble initially located in a region with pressure below the vapor pressure'
            print *, 'location:', mtn_pos(nparticles, 1:3, 1)
            stop
        end if
        omegaN = dsqrt(omegaN/bub_R0(nparticles)**2)

        cpparticle = concvap*lag_params%cpvapor + (1.0d0 - concvap)*lag_params%cpgas
        kparticle = concvap*lag_params%kvapor + (1.0d0 - concvap)*lag_params%kgas

        ! Mass and heat transfer coefficients (based on Preston 2007)
        PeT = totalmass/volparticle*cpparticle*bub_R0(nparticles)**2*omegaN/kparticle
        gas_betaT(nparticles) = f_transfercoeff(PeT, 1.0d0)*(heatflag)
        PeG = bub_R0(nparticles)**2*omegaN/lag_params%diffcoefvap
        gas_betaC(nparticles) = f_transfercoeff(PeG, 1.0d0)*(massflag)

        ! Terms to work out directly the heat flux in getfluxes
        gas_betaT(nparticles) = gas_betaT(nparticles)*kparticle

        if (gas_mg(nparticles) <= 0.0d0) stop "Negative gas mass in the bubble, check if the bubble is in the domain."

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
                    lag_id(nparticles, 1) = id  !global ID
                    lag_id(nparticles, 2) = nparticles  !local ID
                    mtn_pos(nparticles, 1:3, 1) = inputvals(1:3)
                    mtn_posPrev(nparticles, 1:3, 1) = inputvals(4:6)
                    mtn_vel(nparticles, 1:3, 1) = inputvals(7:9)
                    intfc_rad(nparticles, 1) = inputvals(10)
                    intfc_vel(nparticles, 1) = inputvals(11)
                    bub_R0(nparticles) = inputvals(12)
                    Rmax_stats(nparticles) = inputvals(13)
                    Rmin_stats(nparticles) = inputvals(14)
                    bub_dphidt(nparticles) = inputvals(15)
                    gas_p(nparticles, 1) = inputvals(16)
                    gas_mv(nparticles, 1) = inputvals(17)
                    gas_mg(nparticles) = inputvals(18)
                    gas_betaT(nparticles) = inputvals(19)
                    gas_betaC(nparticles) = inputvals(20)
                    cell = -buff_size
                    call s_locate_cell(mtn_pos(nparticles, 1:3, 1), cell, mtn_s(nparticles, 1:3, 1))
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
        !! @param stage Current stage in the time-stepper algorithm
    subroutine s_compute_el_coupled_solver(q_cons_vf, q_prim_vf, rhs_vf, stage)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer, intent(in) :: stage

        real(kind(0.d0)) :: gammaparticle, vaporflux, heatflux
        integer :: i, j, k, l
        real(kind(0.d0)) :: preterm1, term2, paux, pint, Romega, term1_fac, Rb

        call nvtxStartRange("DYNAMICS-LAGRANGE-BUBBLES")

        !< BUBBLE DYNAMICS

        ! Subgrid p_inf model from Maeda and Colonius (2018).
        if (lag_params%pressure_corrector) then
            ! Calculate velocity potentials (valid for one bubble per cell)
            !$acc parallel loop gang vector default(present) private(k)
            do k = 1, nBubs
                paux = f_pressure_inf(k, q_prim_vf, 2, preterm1, term2, Romega)
                pint = f_pressureliq_int(k)
                pint = pint + 0.5d0*intfc_vel(k, 2)**2

                if (lag_params%cluster_type == 2) then
                    bub_dphidt(k) = (paux - pint) + term2
                    ! Accounting for the potential induced by the bubble averaged over the control volume
                    ! Note that this is based on the incompressible flow assumption near the bubble.
                    Rb = intfc_rad(k, 2)
                    term1_fac = 3.0d0/2.0d0*(Rb*(Romega**2d0 - Rb**2d0))/(Romega**3d0 - Rb**3d0)
                    bub_dphidt(k) = bub_dphidt(k)/(1 - term1_fac)
                end if
            end do
        end if

        ! Gaseous core evolution
        !$acc parallel loop gang vector default(present) private(k) copyin(stage)
        do k = 1, nBubs
            call s_compute_interface_fluxes(k, vaporflux, heatflux, gammaparticle)
            gas_dpdt(k, stage) = -3.0d0*gammaparticle/intfc_rad(k, 2)* &
                                 (gas_p(k, 2)*intfc_vel(k, 2) - &
                                  heatflux - (lag_params%Rvapor*lag_params%Thost)*vaporflux)
            gas_dmvdt(k, stage) = 4.0d0*pi*intfc_rad(k, 2)**2*vaporflux
        end do

        ! Radial motion model
        if (lag_params%bubble_model == 1) then
            !$acc parallel loop gang vector default(present) private(k) copyin(stage)
            do k = 1, nBubs
                call s_compute_KM(k, stage, q_prim_vf)
                intfc_draddt(k, stage) = intfc_vel(k, 2)
            end do
        else
            !$acc parallel loop gang vector default(present) private(k) copyin(stage)
            do k = 1, nBubs
                intfc_dveldt(k, stage) = 0.0d0
                intfc_draddt(k, stage) = 0.0d0
            end do
        end if

        ! Bubbles remain in a fixed position
        !$acc parallel loop collapse(2) gang vector default(present) private(k) copyin(stage)
        do k = 1, nBubs
            do l = 1, 3
                mtn_dposdt(k, l, stage) = 0.0d0
                mtn_dveldt(k, l, stage) = 0.0d0
            end do
        end do

        call nvtxEndRange

        !< EULER-LAGRANGE COUPLING
        call s_smear_voidfraction()
        if (lag_params%solver_approach == 2) call s_add_sources(q_cons_vf, q_prim_vf, rhs_vf)

    end subroutine s_compute_el_coupled_solver

    !>  The purpose of this subroutine is to smear the effect of the bubbles in the Eulerian framework
    subroutine s_smear_voidfraction()

        integer :: i, j, k, l

        call nvtxStartRange("EL-KERNELS")

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, q_beta_idx
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        q_beta%vf(i)%sf(j, k, l) = 0.0d0
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
                    q_beta%vf(1)%sf(j, k, l) = 1.0d0 - q_beta%vf(1)%sf(j, k, l)
                    ! Limiting void fraction given max value
                    q_beta%vf(1)%sf(j, k, l) = max(q_beta%vf(1)%sf(j, k, l), &
                                                   1.d0 - lag_params%valmaxvoid)
                end do
            end do
        end do

        call nvtxEndRange

    end subroutine s_smear_voidfraction

    !>  The purpose of this subroutine is to add the bubbles source terms based on Maeda and Colonius (2018)
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Conservative variables
        !! @param rhs_vf Time derivative of the conservative variables
    subroutine s_add_sources(q_cons_vf, q_prim_vf, rhs_vf)

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
                            if (q_beta%vf(1)%sf(i, j, k) > (1.0d0 - lag_params%valmaxvoid)) then
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
                            if (q_beta%vf(1)%sf(i, j, k) > (1.0d0 - lag_params%valmaxvoid)) then
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
                        if (q_beta%vf(1)%sf(i, j, k) > (1.0d0 - lag_params%valmaxvoid)) then
                            rhs_vf(num_fluids + l)%sf(i, j, k) = rhs_vf(num_fluids + l)%sf(i, j, k) - &
                                                                 (1.0d0 - q_beta%vf(1)%sf(i, j, k))/ &
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
                        q_beta%vf(3)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k)*q_prim_vf(num_fluids + l)%sf(i, j, k)
                    end do
                end do
            end do

            call s_gradient_dir(q_beta%vf(3), q_beta%vf(4), l)

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        if (q_beta%vf(1)%sf(i, j, k) > (1.0d0 - lag_params%valmaxvoid)) then
                            rhs_vf(E_idx)%sf(i, j, k) = rhs_vf(E_idx)%sf(i, j, k) - &
                                                        q_beta%vf(4)%sf(i, j, k)*(1.0d0 - q_beta%vf(1)%sf(i, j, k))/ &
                                                        q_beta%vf(1)%sf(i, j, k)
                        end if
                    end do
                end do
            end do
        end do

    end subroutine s_add_sources

    !>  This subroutine computes the Keller-Miksis equation
        !! @param nparticles Particle identifier
        !! @param step Current time-stage in RKCK stepper, otherwise step = 1
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

        scoord = mtn_s(nparticles, 1:3, 2)
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
        velint = intfc_vel(nparticles, 2) - gas_dmvdt(nparticles, step)/(4.0d0*pi*intfc_rad(nparticles, 2)**2*rhol)

        intfc_dveldt(nparticles, step) = ((1.0d0 + velint/cson)*deltaP/rhol &
                                          + termI &
                                          + gas_dpdt(nparticles, step)*intfc_rad(nparticles, 2)/rhol/cson &
                                          - velint**2*3.0d0/2.0d0*(1.0d0 - velint/3.0d0/cson)) &
                                         /(intfc_rad(nparticles, 2)*(1.0d0 - velint/cson))

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
        if (lag_params%massTransfer_model) then
            concvapint = (lag_params%Rvapor/lag_params%Rgas)*(gas_p(nparticles, 2)/lag_params%pvap - 1.0d0)
            concvapint = 1.0d0/(1.0d0 + concvapint)
        end if

        bubbleTemp = gas_mg(nparticles)*lag_params%Rgas + gas_mv(nparticles, 2)*lag_params%Rvapor
        volbubble = 4.0d0/3.0d0*pi*intfc_rad(nparticles, 2)**3
        bubbleTemp = gas_p(nparticles, 2)*volbubble/bubbleTemp

        gammabubble = concvapint*lag_params%gammavapor + (1.0d0 - concvapint)*lag_params%gammagas
        heatflux = -(gammabubble - 1.0d0)/gammabubble*gas_betaT(nparticles)*(bubbleTemp - lag_params%Thost)/intfc_rad(nparticles, 2)

        avgconc = gas_mv(nparticles, 2)/(gas_mg(nparticles) + gas_mv(nparticles, 2))
        Rmixt = concvapint*lag_params%Rvapor + (1.0d0 - concvapint)*lag_params%Rgas

        concvapint = min(concvapint, 0.99d0)
        vaporflux = (1.0d0 - concvapint)*intfc_rad(nparticles, 2)
        rhogas = (gas_mg(nparticles) + gas_mv(nparticles, 2))/(4.0d0/3.0d0*pi*intfc_rad(nparticles, 2)**3)
        vaporflux = -lag_params%diffcoefvap*gas_betaC(nparticles)*(avgconc - concvapint)*rhogas/vaporflux

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

        pbubble = gas_p(nparticles, 2)
        radius = intfc_rad(nparticles, 2)
        bubblevel = intfc_vel(nparticles, 2)

        f_pressureliq_int = pbubble - 2.0d0*lag_params%sigmabubble/radius - 4.0d0*lag_params%vischost*bubblevel/radius

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
        real(kind(0.d0)) :: dc, vol, aux, dist
        real(kind(0.d0)) :: volgas, term1, Rbeq, denom
        real(kind(0.d0)) :: charvol, charpres, charvol2, charpres2
        integer, dimension(3) :: cell, cellaux

        integer :: i, j, k, dir
        integer :: cellx, celly, cellz
        integer :: mapCells_pinf, smearGrid, smearGridz
        logical :: celloutside

        scoord = mtn_s(nparticles, 1:3, 2)
        f_pressure_inf = 0.0d0
        call s_get_cell(scoord, cell)

        if ((lag_params%cluster_type == 1)) then
            !getting p_cell in terms of only the current cell by interpolation

            call s_get_char_vol(cell(1), cell(2), cell(3), vol)

            ! Getting the cell volulme as Omega
            call s_get_char_dist(cell(1), cell(2), cell(3), dist)

            !p_cell (interpolated)
            f_pressure_inf = f_interpolate(scoord, q_prim_vf(E_idx))

            !R_Omega
            dc = (3.0d0*vol/(4.0d0*pi))**(1.0d0/3.0d0)

        else if (lag_params%cluster_type >= 2) then
            ! Bubble dynamic closure from Maeda and Colonius (2018)

            ! Range of cells included in Omega
            if (lag_params%smooth_type == 1) then
                mapCells_pinf = mapCells
            else
                stop "lag_params%cluster_type: 2 requires lag_params%smooth_type: 1."
            end if

            smearGrid = mapCells_pinf - (-mapCells_pinf) + 1 ! Include the cell that contains the bubble (mapCells+1+mapCells)
            smearGridz = smearGrid
            if (p == 0) smearGridz = 1

            charvol = 0.d0
            charpres = 0.d0
            charvol2 = 0.d0
            charpres2 = 0.d0
            vol = 0.d0

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

                        call s_check_celloutside(cellaux, celloutside)
                        if (.not. celloutside) then
                            if (cyl_coord .and. (p == 0) .and. (y_cc(cellaux(2)) < 0d0)) then
                                celloutside = .true.
                            end if
                        end if

                        if (.not. celloutside) then
                            call s_get_char_vol(cellaux(1), cellaux(2), cellaux(3), vol)
                            charvol = charvol + vol
                            charpres = charpres + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3))*vol
                            charvol2 = charvol2 + vol*q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                            charpres2 = charpres2 + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                        *vol*q_beta%vf(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                        end if

                    end do
                end do
            end do

            f_pressure_inf = charpres2/charvol2
            vol = charvol
            dc = (3.0d0*abs(vol)/(4.0d0*pi))**(1.0d0/3.0d0)
        else

            stop "Check cluterflag. Exiting ..."

        end if

        if (lag_params%pressure_corrector .and. present(preterm1)) then

            !Valid if only one bubble exists per cell
            volgas = intfc_rad(nparticles, 2)**3
            denom = intfc_rad(nparticles, 2)**2
            term1 = bub_dphidt(nparticles)*intfc_rad(nparticles, 2)**2
            term2 = intfc_vel(nparticles, 2)*intfc_rad(nparticles, 2)**2

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

    !>  This subroutine updates the Lagrange variables in the tvd RK time steppers.
        !! @param stage Current tvd RK stage
    subroutine s_update_lag_tdv_rk(stage)

        integer, intent(in) :: stage

        integer :: k, l, q

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
                if (intfc_rad(k, 1) <= 0.0d0) stop "Negative bubble radius encountered, please reduce dt"
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
                    if (intfc_rad(k, 2) <= 0.0d0) stop "Negative bubble radius encountered, please reduce dt"
                end do

            elseif (stage == 2) then
                !$acc parallel loop gang vector default(present) private(k)
                do k = 1, nBubs
                    !u{1} = u{n} + (1/2) * dt * (RHS{n} + RHS{1})
                    intfc_rad(k, 1) = intfc_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/2d0
                    intfc_vel(k, 1) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/2d0
                    mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/2d0
                    mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/2d0
                    gas_p(k, 1) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/2d0
                    gas_mv(k, 1) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/2d0
                    if (intfc_rad(k, 1) <= 0.0d0) stop "Negative bubble radius encountered, please reduce dt"
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
                    if (intfc_rad(k, 2) <= 0.0d0) stop "Negative bubble radius encountered, please reduce dt"
                end do

            elseif (stage == 2) then
                !$acc parallel loop gang vector default(present) private(k)
                do k = 1, nBubs
                    !u{2} = u{n} + (1/4) * dt * [RHS{n} + RHS{1}]
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/4d0
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/4d0
                    mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/4d0
                    mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/4d0
                    gas_p(k, 2) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/4d0
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/4d0
                    if (intfc_rad(k, 2) <= 0.0d0) stop "Negative bubble radius encountered, please reduce dt"
                end do
            elseif (stage == 3) then
                !$acc parallel loop gang vector default(present) private(k)
                do k = 1, nBubs
                    !u{n+1} = u{n} + (2/3) * dt * [(1/4)* RHS{n} + (1/4)* RHS{1} + RHS{2}]
                    intfc_rad(k, 1) = intfc_rad(k, 1) + (2d0/3d0)*dt*(intfc_draddt(k, 1)/4d0 + intfc_draddt(k, 2)/4d0 + intfc_draddt(k, 3))
                    intfc_vel(k, 1) = intfc_vel(k, 1) + (2d0/3d0)*dt*(intfc_dveldt(k, 1)/4d0 + intfc_dveldt(k, 2)/4d0 + intfc_dveldt(k, 3))
                    mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + (2d0/3d0)*dt*(mtn_dposdt(k, 1:3, 1)/4d0 + mtn_dposdt(k, 1:3, 2)/4d0 + mtn_dposdt(k, 1:3, 3))
                    mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + (2d0/3d0)*dt*(mtn_dveldt(k, 1:3, 1)/4d0 + mtn_dveldt(k, 1:3, 2)/4d0 + mtn_dveldt(k, 1:3, 3))
                    gas_p(k, 1) = gas_p(k, 1) + (2d0/3d0)*dt*(gas_dpdt(k, 1)/4d0 + gas_dpdt(k, 2)/4d0 + gas_dpdt(k, 3))
                    gas_mv(k, 1) = gas_mv(k, 1) + (2d0/3d0)*dt*(gas_dmvdt(k, 1)/4d0 + gas_dmvdt(k, 2)/4d0 + gas_dmvdt(k, 3))
                    if (intfc_rad(k, 1) <= 0.0d0) stop "Negative bubble radius encountered, please reduce dt"
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

    end subroutine s_update_lag_tdv_rk

    !>  This subroutine updates the Euler-Lagrange temporal variables before entering to the next time-stage in the RKCK stepper.
        !! @param RKstep Current time step in the adaptative stepper
        !! @param q_cons_ts Conservative variables
        !! @param rhs_ts Time derivatives of the conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_update_tmp_rkck(RKstep, q_cons_ts, rhs_ts, q_prim_vf, lag_largestep)

        integer, intent(in) :: RKstep
        type(vector_field), dimension(:), intent(inout) :: q_cons_ts
        type(vector_field), dimension(:), intent(inout) :: rhs_ts
        type(scalar_field), dimension(:), intent(inout) :: q_prim_vf
        real(kind(0.d0)), intent(out) :: lag_largestep

        integer :: i, j, k, l, q
        integer, dimension(3) :: cell
        real(kind(0.d0)) :: radiusOld, velOld, aux_glb
        integer :: remove_id

        call s_transfer_data_to_tmp()

        lag_largestep = 0.0d0
        remove_id = 0
        !$acc parallel loop gang vector default(present) reduction(+: lag_largestep, remove_id) private(k) copyin(RKstep)
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

            if ((intfc_rad(k, 2) <= 0.0d0) .or. &                       ! no negative radius
                (mtn_pos(k, 1, 2) /= mtn_pos(k, 1, 2))) then  ! finite bubble location
                print *, 'Negative bubble radius encountered'
                lag_largestep = lag_largestep + 1.0d0
                if (dt < 2.0d0*verysmall_dt) then
                    remove_id = k
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

        if (lag_largestep > 0.0d0) return

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
    subroutine s_calculate_rkck_truncation_error(rkck_errmax)

        real(kind(0.d0)), intent(out) :: rkck_errmax

        real(kind(0.d0)) :: erraux, errb
        integer :: i, j, k, l, l1, nb

        rkck_errmax = 0.0d0
        !$acc parallel loop gang vector default(present) reduction(MAX: rkck_errmax) private(k)
        do k = 1, nBubs
            errb = 0.0d0
            !Bubble radius error
            erraux = 0.0d0
            !$acc loop seq
            do i = 1, lag_num_ts
                erraux = erraux + lag_RKCKcoef(7, i)*intfc_draddt(k, i)
            end do
            errb = max(errb, abs(erraux)*dt/bub_R0(k))

            !Interface velocity error
            erraux = 0.0d0
            !$acc loop seq
            do i = 1, lag_num_ts
                erraux = erraux + lag_RKCKcoef(7, i)*intfc_dveldt(k, i)
            end do
            errb = max(errb, abs(erraux)*dt)

            !Bubble velocity error
            !$acc loop seq
            do j = 1, 3
                erraux = 0.0d0
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
        !! @param q_prim_vf Primitive variables
    subroutine s_update_rkck(q_cons_ts, q_prim_vf)

        type(vector_field), dimension(:), intent(inout) :: q_cons_ts
        type(scalar_field), dimension(:), intent(in) :: q_prim_vf

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

    !>  This subroutine computes the next time step in the adaptive RKCK stepper.
        !! @param lag_largestep Indicates that a negative radius was calculated
        !! @param restart_rkck_step Flag to restart the current time step
        !! @param rkck_errmax Maximum error between the 4th and 5th RKCK results.
    subroutine s_compute_rkck_dt(lag_largestep, restart_rkck_step, rkck_errmax)

        real(kind(0.d0)), intent(in) :: lag_largestep
        logical, intent(out) :: restart_rkck_step
        real(kind(0.d0)), intent(inout), optional :: rkck_errmax

        real(kind(0.d0)) :: htemp, aux_glb

        restart_rkck_step = .false.

        if (lag_largestep > 0.0d0) then ! Encountered negative radius, so reduce dt and restart time step

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

            rkck_errmax = min(rkck_errmax, 1.0d0)
            if (num_procs > 1) then
                call s_mpi_allreduce_max(rkck_errmax, aux_glb)
                rkck_errmax = aux_glb
            end if
            rkck_errmax = rkck_errmax/lag_params%rkck_tolerance ! Scale relative to user required tolerance.

            if ((rkck_errmax > 1.0d0)) then   ! Truncation error too large, reduce dt and restart time step
                restart_rkck_step = .true.
                htemp = SAFETY*dt*((floor(rkck_errmax*RNDDEC)/RNDDEC)**PSHRNK)
                dt = sign(max(abs(htemp), SHRNKDTMAX*abs(dt)), dt)  ! No more than a factor of 10.
                if (proc_rank == 0) print *, '>>>>> WARNING: Truncation error found. Reducing dt and restaring time step, now dt: ', dt
            else                            ! Step succeeded. Compute size of next step.
                if (rkck_errmax > ERRCON) then
                    dt = SAFETY*dt*((floor(rkck_errmax*RNDDEC)/RNDDEC)**PGROW) ! No more than a factor of 5 increase.
                else
                    dt = (1.0d0/SHRNKDT)*dt ! Truncation error too small (< 1.89e-4), increase time step
                end if
            end if

            dt = min(dt, dt_max)
            if (num_procs > 1) then
                call s_mpi_allreduce_min(dt, aux_glb)
                dt = aux_glb
            end if
            !$acc update device(dt)

        end if

    end subroutine s_compute_rkck_dt

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
        real(kind(0.d0)), dimension(3) :: pos_part

        ! 2D
        if (p == 0 .and. cyl_coord .neqv. .true.) then
            ! Defining a virtual z-axis that has the same dimensions as y-axis
            ! defined in the input file
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
                                  (pos_part(3) < lag_params%charwidth/2d0) .and. (pos_part(3) >= -lag_params%charwidth/2d0))
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

        real(kind(0.d0)) :: qtime, volcell, voltot
        real(kind(0.d0)) :: lag_voidmax, lag_voidavg, lag_vol
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
        !$acc parallel loop collapse(3) gang vector default(present) reduction(+:lag_vol,lag_voidavg) &
        !$acc reduction(MAX:lag_voidmax) private(cell)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    lag_voidmax = max(lag_voidmax, 1.0d0 - q_beta%vf(1)%sf(i, j, k))
                    cell(1) = i
                    cell(2) = j
                    cell(3) = k
                    call s_get_char_vol(cell(1), cell(2), cell(3), volcell)
                    if ((1.0d0 - q_beta%vf(1)%sf(i, j, k)) > 5.0d-11) then
                        lag_voidavg = lag_voidavg + (1.0d0 - q_beta%vf(1)%sf(i, j, k))*volcell
                        lag_vol = lag_vol + volcell
                    end if
                end do
            end do
        end do

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
            write (12, '(6X,4e24.8)') &
                qtime, &
                lag_voidavg, &
                lag_voidmax, &
                voltot
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
                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then
                    nparticles = nparticles + 1
                end if
            end do
        end if

        if (.not. parallel_io) return

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
            Rmax_glb = max(Rmax_glb, intfc_rad(k, 1)/bub_R0(k))
            Rmin_glb = min(Rmin_glb, intfc_rad(k, 1)/bub_R0(k))
            Rmax_stats(k) = max(Rmax_stats(k), intfc_rad(k, 1)/bub_R0(k))
            Rmin_stats(k) = min(Rmin_stats(k), intfc_rad(k, 1)/bub_R0(k))
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
          !! @param nparticles Particle id
    subroutine s_remove_lag_bubble(nparticles)

        integer, intent(in) :: nparticles

        integer :: i

        !$acc loop seq
        do i = nparticles, nBubs - 1
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

        integer :: i, j, k, imax

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
