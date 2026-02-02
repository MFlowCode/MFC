!>
!! @file m_particles_EL.fpp
!! @brief Contains module m_particles_EL

#:include 'macros.fpp'

!> @brief This module is used to to compute the volume-averaged particle model
module m_particles_EL

    use m_global_parameters             !< Definitions of the global parameters

    use m_mpi_proxy                     !< Message passing interface (MPI) module proxy

    use m_particles_EL_kernels            !< Definitions of the kernel functions

    use m_particles                       !< General bubble dynamics procedures

    use m_variables_conversion          !< State variables type conversion procedures

    use m_compile_specific

    use m_boundary_common

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_sim_helpers

    use m_helper

    use m_mpi_common

    use m_ibm

    implicit none

    !(nBub)
    integer, allocatable, dimension(:, :) :: lag_id                 !< Global and local IDs
    real(wp), allocatable, dimension(:) :: particle_R0            !< Initial particle radius
    real(wp), allocatable, dimension(:) :: Rmax_stats        !< Maximum radius
    real(wp), allocatable, dimension(:) :: Rmin_stats        !< Minimum radius
    $:GPU_DECLARE(create='[lag_id, particle_R0, Rmax_stats, Rmin_stats]')

    real(wp), allocatable, dimension(:) :: particle_mass            !< Particle Mass
    real(wp), allocatable, dimension(:) :: gas_betaT         !< heatflux model (Preston et al., 2007)
    real(wp), allocatable, dimension(:) :: gas_betaC         !< massflux model (Preston et al., 2007)
    ! real(wp), allocatable, dimension(:) :: bub_dphidt        !< subgrid velocity potential (Maeda & Colonius, 2018)
    $:GPU_DECLARE(create='[particle_mass, gas_betaT, gas_betaC]')
    ! bub_dphidt removed

    !(nPart, 1 -> actual val or 2 -> temp val)
    ! real(wp), allocatable, dimension(:, :) :: gas_p          !< Pressure in the bubble
    ! real(wp), allocatable, dimension(:, :) :: gas_mv         !< Vapor mass in the bubble
    real(wp), allocatable, dimension(:, :) :: particle_rad      !< Particle radius
    ! real(wp), allocatable, dimension(:, :) :: intfc_vel      !< Velocity of the bubble interface
    ! $:GPU_DECLARE(create='[gas_p, gas_mv, intfc_rad, intfc_vel]')
    $:GPU_DECLARE(create='[particle_rad]')
    !New for particles
    !(nPart, 1-> x or 2->y or 3 ->z, 1 -> actual or 2 -> temporal val)
    real(wp), allocatable, dimension(:, :, :) :: mtn_pos     !< Particle's position
    real(wp), allocatable, dimension(:, :, :) :: mtn_posPrev !< Particle's previous position
    real(wp), allocatable, dimension(:, :, :) :: mtn_vel     !< Particle's velocity
    real(wp), allocatable, dimension(:, :, :) :: mtn_velPrev !< Particle's previous velocity
    real(wp), allocatable, dimension(:, :, :) :: mtn_s       !< Particle's computational cell position in real format
    $:GPU_DECLARE(create='[mtn_pos, mtn_posPrev, mtn_vel, mtn_velPrev, mtn_s]')
    !(nPart, 1-> x or 2->y or 3 ->z, time-stage)
    real(wp), allocatable, dimension(:, :) :: intfc_draddt   !< Time derivative of bubble's radius
    ! real(wp), allocatable, dimension(:, :) :: intfc_dveldt   !< Time derivative of bubble's interface velocity
    ! real(wp), allocatable, dimension(:, :) :: gas_dpdt       !< Time derivative of gas pressure
    ! real(wp), allocatable, dimension(:, :) :: gas_dmvdt      !< Time derivative of the vapor mass in the bubble
    real(wp), allocatable, dimension(:, :, :) :: mtn_dposdt  !< Time derivative of the particle's position
    real(wp), allocatable, dimension(:, :, :) :: mtn_dveldt  !< Time derivative of the particle's velocity
    ! $:GPU_DECLARE(create='[intfc_draddt, intfc_dveldt, gas_dpdt, gas_dmvdt, mtn_dposdt, mtn_dveldt]')
    $:GPU_DECLARE(create='[intfc_draddt, mtn_dposdt, mtn_dveldt]')

    integer, private :: lag_num_ts                                  !<  Number of time stages in the time-stepping scheme
    $:GPU_DECLARE(create='[lag_num_ts]')

    real(wp) :: Rmax_glb, Rmin_glb       !< Maximum and minimum bubbe size in the local domain
    !< Projection of the lagrangian particles in the Eulerian framework
    type(scalar_field), dimension(:), allocatable :: q_particles
    integer :: q_particles_idx                     !< Size of the q vector field for particle cell (q)uantities

    $:GPU_DECLARE(create='[Rmax_glb,Rmin_glb,q_particles,q_particles_idx]')

    !Particle Source terms for fluid coupling
    real(wp), allocatable, dimension (:,:) :: f_p !< force on each particle
    $:GPU_DECLARE(create='[f_p]')

    integer, parameter :: LAG_EVOL_ID = 11 ! File id for lag_bubbles_evol_*.dat
    integer, parameter :: LAG_STATS_ID = 12 ! File id for stats_lag_bubbles_*.dat
    integer, parameter :: LAG_VOID_ID = 13 ! File id for voidfraction.dat

    integer, allocatable, dimension(:) :: keep_bubble
    integer, allocatable, dimension(:, :) :: wrap_bubble_loc, wrap_bubble_dir
    $:GPU_DECLARE(create='[keep_bubble]')
    $:GPU_DECLARE(create='[wrap_bubble_loc, wrap_bubble_dir]')

contains

    !> Initializes the lagrangian subgrid particle solver
        !! @param q_cons_vf Initial conservative variables
    impure subroutine s_initialize_particles_EL_module(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        integer :: nParticles_glb, i

        ! Setting number of time-stages for selected time-stepping scheme
        lag_num_ts = time_stepper

        ! Allocate space for the Eulerian fields needed to map the effect of the particles
        if (lag_params%solver_approach == 1) then
            ! One-way coupling
            q_particles_idx = 1 !For tracking volume fraction
        elseif (lag_params%solver_approach == 2) then
            ! Two-way coupling
            q_particles_idx = 4 !For tracking volume fraction, x-mom, y-mom, and energy sources
            if (p > 0) then
                q_particles_idx = 5 !For tracking volume fraction, x-mom, y-mom, z-mom, and energy sources
            end if
        else
            call s_mpi_abort('Please check the lag_params%solver_approach input')
        end if

        pcomm_coords(1)%beg = x_cb(mapcells)
        pcomm_coords(1)%end = x_cb(m - mapcells - 1)
        $:GPU_UPDATE(device='[pcomm_coords(1)]')
        if (n > 0) then
            pcomm_coords(2)%beg = y_cb(mapcells)
            pcomm_coords(2)%end = y_cb(n - mapcells - 1)
            $:GPU_UPDATE(device='[pcomm_coords(2)]')
            if (p > 0) then
                pcomm_coords(3)%beg = z_cb(mapCells)
                pcomm_coords(3)%end = z_cb(p - mapCells - 1)
                $:GPU_UPDATE(device='[pcomm_coords(3)]')
            end if
        end if

        $:GPU_UPDATE(device='[lag_num_ts, q_particles_idx]')

        @:ALLOCATE(q_particles(1:q_particles_idx))

        do i = 1, q_particles_idx
            @:ALLOCATE(q_particles(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_particles(i))
        end do

        ! Allocating space for lagrangian variables
        nParticles_glb = lag_params%nParticles_glb

        @:ALLOCATE(lag_id(1:nParticles_glb, 1:2))
        @:ALLOCATE(particle_R0(1:nParticles_glb))
        @:ALLOCATE(Rmax_stats(1:nParticles_glb))
        @:ALLOCATE(Rmin_stats(1:nParticles_glb))
        @:ALLOCATE(particle_mass(1:nParticles_glb))
        @:ALLOCATE(gas_betaT(1:nParticles_glb))
        @:ALLOCATE(gas_betaC(1:nParticles_glb))
        ! @:ALLOCATE(bub_dphidt(1:nParticles_glb))
        ! @:ALLOCATE(gas_p(1:nParticles_glb, 1:2))
        ! @:ALLOCATE(gas_mv(1:nParticles_glb, 1:2))
        @:ALLOCATE(particle_rad(1:nParticles_glb, 1:2))
        ! @:ALLOCATE(intfc_vel(1:nParticles_glb, 1:2))
        @:ALLOCATE(mtn_pos(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_posPrev(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_velPrev(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_vel(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_s(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(intfc_draddt(1:nParticles_glb, 1:lag_num_ts))
        ! @:ALLOCATE(intfc_dveldt(1:nParticles_glb, 1:lag_num_ts))
        ! @:ALLOCATE(gas_dpdt(1:nParticles_glb, 1:lag_num_ts))
        ! @:ALLOCATE(gas_dmvdt(1:nParticles_glb, 1:lag_num_ts))
        @:ALLOCATE(mtn_dposdt(1:nParticles_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(mtn_dveldt(1:nParticles_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(f_p(1:nParticles_glb, 1:3)) 

        @:ALLOCATE(keep_bubble(1:nParticles_glb))
        @:ALLOCATE(wrap_bubble_loc(1:nParticles_glb, 1:num_dims), wrap_bubble_dir(1:nParticles_glb, 1:num_dims))

        if (adap_dt .and. f_is_default(adap_dt_tol)) adap_dt_tol = dflt_adap_dt_tol

        if (num_procs > 1) call s_initialize_solid_particles_mpi(lag_num_ts)

        ! Starting bubbles
        if (lag_params%write_void_evol) call s_open_void_evol
        if (lag_params%write_bubbles) call s_open_lag_bubble_evol()
        if (lag_params%write_bubbles_stats) call s_open_lag_particle_stats()

        ! if (lag_params%vel_model > 0) then
        !     moving_lag_particles = .true.
        !     lag_pressure_force = lag_params%pressure_force
        !     lag_gravity_force = lag_params%gravity_force
        !     lag_vel_model = lag_params%vel_model
        !     lag_drag_model = lag_params%drag_model
        ! end if
        moving_lag_particles = .true.

        $:GPU_UPDATE(device='[moving_lag_particles, lag_pressure_force, &
            & lag_gravity_force, lag_vel_model, lag_drag_model]')

        call s_read_input_particles(q_cons_vf)

    end subroutine s_initialize_particles_EL_module

    !> The purpose of this procedure is to obtain the initial bubbles' information
        !! @param q_cons_vf Conservative variables
    impure subroutine s_read_input_particles(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        real(wp), dimension(8) :: inputParticle
        real(wp) :: qtime
        integer :: id, particle_id, save_count
        integer :: i, ios
        logical :: file_exist, indomain
        integer, dimension(3) :: cell

        character(LEN=path_len + 2*name_len) :: path_D_dir !<

        ! Initialize number of particles
        particle_id = 0
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
            if (proc_rank == 0) print *, 'Reading lagrange particles input file.'
            call my_inquire(trim(lag_params%input_path), file_exist)
            if (file_exist) then
                open (94, file=trim(lag_params%input_path), form='formatted', iostat=ios)
                do while (ios == 0)
                    read (94, *, iostat=ios) (inputParticle(i), i=1, 8)
                    if (ios /= 0) cycle
                    indomain = particle_in_domain_physical(inputParticle(1:3))
                    id = id + 1
                    if (id > lag_params%nParticles_glb .and. proc_rank == 0) then
                        call s_mpi_abort("Current number of particles is larger than nParticles_glb")
                    end if
                    if (indomain) then
                        particle_id = particle_id + 1
                        call s_add_particles(inputParticle, q_cons_vf, particle_id)
                        lag_id(particle_id, 1) = id      !global ID
                        lag_id(particle_id, 2) = particle_id  !local ID
                        n_el_particles_loc = particle_id              ! local number of bubbles
                    end if
                end do
                close (94)
            else
                call s_mpi_abort("Initialize the lagrange particles in "//trim(lag_params%input_path))
            end if
        else
            if (proc_rank == 0) print *, 'Restarting lagrange particles at save_count: ', save_count
            call s_restart_bubbles(particle_id, save_count)
        end if

        print *, " Lagrange parrticles running, in proc", proc_rank, "number:", particle_id, "/", id

        if (num_procs > 1) then
            call s_mpi_reduce_int_sum(n_el_particles_loc, n_el_particles_glb)
        else
            n_el_particles_glb = n_el_particles_loc
        end if

        if (proc_rank == 0) then
            if (n_el_particles_glb == 0) call s_mpi_abort('No particles in the domain. Check '//trim(lag_params%input_path))
        end if

        if (num_procs > 1) then
            call s_add_particles_to_transfer_list(n_el_particles_loc, mtn_pos(:, :, 1))
            ! call s_mpi_sendrecv_particles(bub_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, &
            !                               gas_betaC, bub_dphidt, lag_id, gas_p, gas_mv, &
            !                               intfc_rad, intfc_vel, mtn_pos, mtn_posPrev, mtn_vel, &
            !                               mtn_s, intfc_draddt, intfc_dveldt, gas_dpdt, &
            !                               gas_dmvdt, mtn_dposdt, mtn_dveldt, lag_num_ts, n_el_particles_loc, &
            !                               dest=1)

            call s_mpi_sendrecv_solid_particles(particle_R0, Rmax_stats, Rmin_stats, particle_mass, gas_betaT, &
                                                gas_betaC, lag_id, &
                                                particle_rad, mtn_pos, mtn_posPrev, mtn_vel, mtn_velPrev, &
                                                mtn_s, intfc_draddt, mtn_dposdt, mtn_dveldt, lag_num_ts, n_el_particles_loc, &
                                                dest=1)

        end if

        $:GPU_UPDATE(device='[bubbles_lagrange, lag_params]')

        ! $:GPU_UPDATE(device='[lag_id,bub_R0,Rmax_stats,Rmin_stats,gas_mg, &
        !     & gas_betaT,gas_betaC,bub_dphidt,gas_p,gas_mv, &
        !     & intfc_rad,intfc_vel,mtn_pos,mtn_posPrev,mtn_vel, &
        !     & mtn_s,intfc_draddt,intfc_dveldt,gas_dpdt,gas_dmvdt, &
        !     & mtn_dposdt,mtn_dveldt,n_el_particles_loc]')

        $:GPU_UPDATE(device='[lag_id,particle_R0,Rmax_stats,Rmin_stats,particle_mass, &
            & gas_betaT,gas_betaC, &
            & particle_rad,mtn_pos,mtn_posPrev,mtn_vel,mtn_velPrev, &
            & mtn_s,intfc_draddt, &
            & mtn_dposdt,mtn_dveldt,n_el_particles_loc]')

        ! !$:GPU_PARALLEL_LOOP(private='[cell]')
        do i = 1, n_el_particles_loc
            cell = fd_number - buff_size
            call s_locate_cell(mtn_pos(i, 1:3, 1), cell, mtn_s(i, 1:3, 1))
        end do

        Rmax_glb = min(dflt_real, -dflt_real)
        Rmin_glb = max(dflt_real, -dflt_real)
        $:GPU_UPDATE(device='[Rmax_glb, Rmin_glb]')

        $:GPU_UPDATE(device='[dx,dy,dz,x_cb,x_cc,y_cb,y_cc,z_cb,z_cc]')

        !Populate temporal variables
        call s_transfer_data_to_tmp_particles()
        call s_smear_particle_sources()

        if (save_count == 0) then
            ! Create ./D directory
            if (proc_rank == 0) then
                write (path_D_dir, '(A,I0,A,I0)') trim(case_dir)//'/D'
                call my_inquire(trim(path_D_dir), file_exist)
                if (.not. file_exist) call s_create_directory(trim(path_D_dir))
            end if
            call s_mpi_barrier()
            call s_write_restart_lag_particles(save_count) ! Needed for post_processing
            if (lag_params%write_void_evol) call s_write_void_evol_particles(qtime)
        end if

        if (lag_params%write_bubbles) call s_write_lag_particle_evol(qtime)

    end subroutine s_read_input_particles

    !> The purpose of this procedure is to obtain the information of the particles when starting fresh
        !! @param inputPart Particle information
        !! @param q_cons_vf Conservative variables
        !! @param part_id Local id of the particle
    impure subroutine s_add_particles(inputPart, q_cons_vf, part_id)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), dimension(8), intent(in) :: inputPart
        integer, intent(in) :: part_id
        integer :: i

        real(wp) :: pliq, volparticle, concvap, totalmass, kparticle, cpparticle
        real(wp) :: omegaN_local, PeG, PeT, rhol, pcrit, qv, gamma, pi_inf, dynP
        integer, dimension(3) :: cell
        real(wp), dimension(2) :: Re
        real(wp) :: massflag, heatflag, Re_trans, Im_trans

        massflag = 0._wp
        heatflag = 0._wp
        if (lag_params%massTransfer_model) massflag = 1._wp
        if (lag_params%heatTransfer_model) heatflag = 1._wp

        particle_R0(part_id) = inputPart(7)
        Rmax_stats(part_id) = min(dflt_real, -dflt_real)
        Rmin_stats(part_id) = max(dflt_real, -dflt_real)
        ! bub_dphidt(part_id) = 0._wp
        particle_rad(part_id, 1) = inputPart(7)
        ! intfc_vel(part_id, 1) = inputPart(8)
        mtn_pos(part_id, 1:3, 1) = inputPart(1:3)
        mtn_posPrev(part_id, 1:3, 1) = mtn_pos(part_id, 1:3, 1)
        mtn_vel(part_id, 1:3, 1) = inputPart(4:6)
        mtn_velPrev(part_id, 1:3, 1) = inputPart(4:6)

        !Initialize Particle Sources
        f_p(part_id,1:3) = 0._wp

        if (cyl_coord .and. p == 0) then
            mtn_pos(part_id, 2, 1) = sqrt(mtn_pos(part_id, 2, 1)**2._wp + &
                                          mtn_pos(part_id, 3, 1)**2._wp)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            mtn_pos(part_id, 3, 1) = atan2(inputPart(3), inputPart(2))
            mtn_posPrev(part_id, 1:3, 1) = mtn_pos(part_id, 1:3, 1)
        end if

        cell = fd_number - buff_size
        call s_locate_cell(mtn_pos(part_id, 1:3, 1), cell, mtn_s(part_id, 1:3, 1))

        ! Check if the bubble is located in the ghost cell of a symmetric, or wall boundary
        if ((any(bc_x%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(1) < 0) .or. &
            (any(bc_x%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(1) > m) .or. &
            (any(bc_y%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(2) < 0) .or. &
            (any(bc_y%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(2) > n)) then
            call s_mpi_abort("Lagrange bubble is in the ghost cells of a symmetric or wall boundary.")
        end if

        if (p > 0) then
            if ((any(bc_z%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(3) < 0) .or. &
                (any(bc_z%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(3) > p)) then
                call s_mpi_abort("Lagrange bubble is in the ghost cells of a symmetric or wall boundary.")
            end if
        end if

        ! call s_convert_to_mixture_variables(q_cons_vf, cell(1), cell(2), cell(3), &
        !                                     rhol, gamma, pi_inf, qv, Re)
        ! dynP = 0._wp
        ! do i = 1, num_dims
        !     dynP = dynP + 0.5_wp*q_cons_vf(contxe + i)%sf(cell(1), cell(2), cell(3))**2/rhol
        ! end do
        ! pliq = (q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)) - dynP - pi_inf)/gamma
        ! if (pliq < 0) print *, "Negative pressure", proc_rank, &
        !     q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)), pi_inf, gamma, pliq, cell, dynP

        ! Initial particle mass
        volparticle = 4._wp/3._wp*pi*particle_R0(part_id)**3 ! volume
        ! gas_mv(particle_id, 1) = pv*volparticle*(1._wp/(R_v*Tw))*(massflag) ! vapermass
        particle_mass(part_id) = volparticle*rho0ref_particle !(gas_p(particle_id, 1) - pv*(massflag))*volparticle*(1._wp/(R_n*Tw)) ! gasmass
        if (particle_mass(part_id) <= 0._wp) then
            call s_mpi_abort("The initial particle mass is negative. Check the initial conditions.")
        end if
        totalmass = particle_mass(part_id) !+ gas_mv(particle_id, 1) ! totalmass


    end subroutine s_add_particles

    !> The purpose of this procedure is to obtain the information of the bubbles from a restart point.
        !! @param part_id Local ID of the particle
        !! @param save_count File identifier
    impure subroutine s_restart_bubbles(part_id, save_count)

        integer, intent(inout) :: part_id, save_count

        character(LEN=path_len + 2*name_len) :: file_loc
        real(wp) :: file_time, file_dt
        integer :: file_num_procs, file_tot_part, tot_part

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

        integer, dimension(:), allocatable :: proc_particle_counts
        real(wp), dimension(1:1, 1:lag_io_vars) :: dummy
        dummy = 0._wp

        ! Construct file path
        write (file_loc, '(A,I0,A)') 'lag_bubbles_', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)

        ! Check if file exists
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort('Restart file '//trim(file_loc)//' does not exist!')
        end if

        if (.not. parallel_io) return

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)

            call MPI_FILE_READ(ifile, file_tot_part, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_READ(ifile, file_time, 1, mpi_p, status, ierr)
            call MPI_FILE_READ(ifile, file_dt, 1, mpi_p, status, ierr)
            call MPI_FILE_READ(ifile, file_num_procs, 1, MPI_INTEGER, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        call MPI_BCAST(file_tot_part, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(file_time, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(file_dt, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(file_num_procs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        allocate (proc_particle_counts(file_num_procs))

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)

            ! Skip to processor counts position
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs), &
                       MPI_OFFSET_KIND)
            call MPI_FILE_SEEK(ifile, disp, MPI_SEEK_SET, ierr)
            call MPI_FILE_READ(ifile, proc_particle_counts, file_num_procs, MPI_INTEGER, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        call MPI_BCAST(proc_particle_counts, file_num_procs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! Set time variables from file
        mytime = file_time
        dt = file_dt

        part_id = proc_particle_counts(proc_rank + 1)

        start_idx_part(1) = 0
        do i = 1, proc_rank
            start_idx_part(1) = start_idx_part(1) + proc_particle_counts(i)
        end do

        start_idx_part(2) = 0
        lsizes(1) = part_id
        lsizes(2) = lag_io_vars

        gsizes(1) = file_tot_part
        gsizes(2) = lag_io_vars

        if (part_id > 0) then

            allocate (MPI_IO_DATA_lag_bubbles(part_id, 1:lag_io_vars))

            call MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                          MPI_ORDER_FORTRAN, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)

            ! Skip extended header
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs) + &
                       file_num_procs*sizeof(proc_particle_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA_lag_bubbles, &
                                   lag_io_vars*part_id, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
            call MPI_TYPE_FREE(view, ierr)

            n_el_particles_loc = part_id

            do i = 1, part_id
                lag_id(i, 1) = int(MPI_IO_DATA_lag_bubbles(i, 1))
                mtn_pos(i, 1:3, 1) = MPI_IO_DATA_lag_bubbles(i, 2:4)
                mtn_posPrev(i, 1:3, 1) = MPI_IO_DATA_lag_bubbles(i, 5:7)
                mtn_vel(i, 1:3, 1) = MPI_IO_DATA_lag_bubbles(i, 8:10)
                particle_rad(i, 1) = MPI_IO_DATA_lag_bubbles(i, 11)
                ! intfc_vel(i, 1) = MPI_IO_DATA_lag_bubbles(i, 12)
                particle_R0(i) = MPI_IO_DATA_lag_bubbles(i, 13)
                Rmax_stats(i) = MPI_IO_DATA_lag_bubbles(i, 14)
                Rmin_stats(i) = MPI_IO_DATA_lag_bubbles(i, 15)
                ! bub_dphidt(i) = MPI_IO_DATA_lag_bubbles(i, 16)
                ! gas_p(i, 1) = MPI_IO_DATA_lag_bubbles(i, 17)
                ! gas_mv(i, 1) = MPI_IO_DATA_lag_bubbles(i, 18)
                particle_mass(i) = MPI_IO_DATA_lag_bubbles(i, 19)
                gas_betaT(i) = MPI_IO_DATA_lag_bubbles(i, 20)
                gas_betaC(i) = MPI_IO_DATA_lag_bubbles(i, 21)
                cell = -buff_size
                call s_locate_cell(mtn_pos(i, 1:3, 1), cell, mtn_s(i, 1:3, 1))
            end do

            deallocate (MPI_IO_DATA_lag_bubbles)

        else
            n_el_particles_loc = 0

            call MPI_TYPE_CONTIGUOUS(0, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)

            ! Skip extended header
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs) + &
                       file_num_procs*sizeof(proc_particle_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_READ_ALL(ifile, dummy, 0, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
            call MPI_TYPE_FREE(view, ierr)
        end if

        if (proc_rank == 0) then
            write (*, '(A,I0,A,I0)') 'Read ', file_tot_part, ' particles from restart file at t_step = ', save_count
            write (*, '(A,E15.7,A,E15.7)') 'Restart time = ', mytime, ', dt = ', dt
        end if

        deallocate (proc_particle_counts)
#endif

    end subroutine s_restart_bubbles

    !>  Contains the particle dynamics subroutines.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
        !! @param rhs_vf Calculated change of conservative variables
        !! @param t_step Current time step
        !! @param stage Current stage in the time-stepper algorithm
    subroutine s_compute_particle_EL_dynamics(q_prim_vf, stage)
#ifdef MFC_OpenMP
        !DIR$ OPTIMIZE (-O1)
#endif
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in) :: stage

        integer, dimension(3) :: cell

        real(wp) :: myMass, myR, myBeta_c, myBeta_t, myR0, myRe, rho_fluid, myVolumeFrac

        real(wp), dimension(3) :: myVel, myPos

        integer :: k, l

        call nvtxStartRange("LAGRANGE-PARTICLE-DYNAMICS")

        $:GPU_PARALLEL_LOOP(private='[k,l,cell,myPb,myMass,myR,myBeta_c,myBeta_t,myR0,myPos,myVel,myRe,rho_fluid,f_p]', copyin='[stage]')
        do k = 1, n_el_particles_loc

            cell = -buff_size
            call s_locate_cell(mtn_pos(k, 1:3, 2), cell, mtn_s(k, 1:3, 2))

            ! Current particle state
            myMass = particle_mass(k)
            myR = particle_rad(k, 2)
            myBeta_c = gas_betaC(k)
            myBeta_t = gas_betaT(k)
            myR0 = particle_R0(k)
            myPos = mtn_pos(k, :, 2)
            myVel = mtn_vel(k, :, 2)
            myRe = 1.48e-5 !fluid_pp(1)%Re(1) !Need a viscosity model for when modeling inviscid eulerian fluid
            myVolumeFrac = 1._wp - q_particles(1)%sf(cell(1), cell(2), cell(3))
            rho_fluid = q_prim_vf(1)%sf(cell(1), cell(2), cell(3))

            do l = 1, num_dims
                f_p(k,l) = f_get_particle_force(myPos(l), myR, myVel, myMass, myRe, rho_fluid, myVolumeFrac, cell, l, q_prim_vf)
                mtn_dposdt(k, l, stage) = myVel(l)
                mtn_dveldt(k, l, stage) = f_p(k,l)/myMass
                intfc_draddt(k, stage) = 0._wp
            end do
            ! $:GPU_ATOMIC(atomic='update')

        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_compute_particle_EL_dynamics

    !>  The purpose of this subroutine is to obtain the bubble source terms based on Maeda and Colonius (2018)
        !!      and add them to the RHS scalar field.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Conservative variables
        !! @param rhs_vf Time derivative of the conservative variables
    subroutine s_compute_particles_EL_source(q_cons_vf, q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        integer :: i, j, k, l

        call s_smear_particle_sources()

        !> Apply particle sources to the Eulerian RHS
        $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
        do k = idwint(3)%beg, idwint(3)%end
            do j = idwint(2)%beg, idwint(2)%end
                do i = idwint(1)%beg, idwint(1)%end
                    if (q_particles(1)%sf(i,j,k) > (1._wp - lag_params%valmaxvoid)) then

                      rhs_vf(momxb)%sf(i,j,k) = rhs_vf(momxb)%sf(i,j,k) + q_particles(2)%sf(i,j,k)
                      rhs_vf(momxb+1)%sf(i,j,k) = rhs_vf(momxb+1)%sf(i,j,k) + q_particles(3)%sf(i,j,k)

                      if (num_dims == 3) then
                        rhs_vf(momxb+2)%sf(i,j,k) = rhs_vf(momxb+2)%sf(i,j,k) + q_particles(4)%sf(i,j,k)
                        ! Energy source
                        rhs_vf(E_idx)%sf(i,j,k) = rhs_vf(E_idx)%sf(i,j,k) &
                                              + q_particles(2)%sf(i,j,k) * q_prim_vf(momxb)%sf(i,j,k) &
                                              + q_particles(3)%sf(i,j,k) * q_prim_vf(momxb+1)%sf(i,j,k) &
                                              + q_particles(4)%sf(i,j,k) * q_prim_vf(momxb+2)%sf(i,j,k) &
                                              + q_particles(5)%sf(i,j,k)
                      else
                        ! Energy source
                        rhs_vf(E_idx)%sf(i,j,k) = rhs_vf(E_idx)%sf(i,j,k) &
                                              + q_particles(2)%sf(i,j,k) * q_prim_vf(momxb)%sf(i,j,k) &
                                              + q_particles(3)%sf(i,j,k) * q_prim_vf(momxb+1)%sf(i,j,k) &
                                              + q_particles(4)%sf(i,j,k)
                      end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_particles_EL_source
    

    !>  The purpose of this subroutine is to smear the effect of the particles in the Eulerian framework
    subroutine s_smear_particle_sources()

        integer :: i, j, k, l

        call nvtxStartRange("BUBBLES-LAGRANGE-KERNELS")
        $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
        do i = 1, q_particles_idx
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        q_particles(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_smoothfunction(n_el_particles_loc, particle_rad, &
                              mtn_s, mtn_pos, mtn_vel, mtn_velPrev, particle_mass, q_particles, f_p)

        !Store 1-beta
        $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    q_particles(1)%sf(j, k, l) = 1._wp - q_particles(1)%sf(j, k, l)
                    ! Limiting void fraction given max value
                    q_particles(1)%sf(j, k, l) = max(q_particles(1)%sf(j, k, l), &
                                                1._wp - lag_params%valmaxvoid)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call nvtxEndRange

        call nvtxStartRange("BUBBLES-LAGRANGE-BETA-COMM")
        if (num_procs > 1) then
            #:for DIRC, DIRI in [('x', 1), ('y', 2), ('z', 3)]
                #:for LOCC, LOCI in [('beg', -1), ('end', 1)]
                    if (bc_${DIRC}$%${LOCC}$ >= 0) then
                        call s_mpi_sendrecv_variables_buffers(q_particles, ${DIRI}$, ${LOCI}$, 2)
                    end if
                #:endfor
            #:endfor
        end if
        call nvtxEndRange

    end subroutine s_smear_particle_sources

    !>  This subroutine updates the Lagrange variables using the tvd RK time steppers.
        !!      The time derivative of the particle variables must be stored at every stage to avoid precision errors.
        !! @param stage Current tvd RK stage
    impure subroutine s_update_lagrange_particles_tdv_rk(q_prim_vf, stage)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: stage

        integer :: k
        moving_lag_particles = .true.

        if (time_stepper == 1) then ! 1st order TVD RK

            $:GPU_PARALLEL_LOOP(private='[k]')
            do k = 1, n_el_particles_loc
                !u{1} = u{n} +  dt * RHS{n}
                particle_rad(k, 1) = particle_rad(k, 1) + dt*intfc_draddt(k, 1)
                ! intfc_vel(k, 1) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                ! gas_p(k, 1) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                ! gas_mv(k, 1) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                if (moving_lag_particles) then
                    mtn_posPrev(k, 1:3, 1) = mtn_pos(k, 1:3, 1)
                    mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                    mtn_velPrev(k, 1:3, 1) = mtn_vel(k, 1:3, 1)
                    mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()

            call s_transfer_data_to_tmp_particles()
            if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf)
            if (lag_params%write_void_evol) call s_write_void_evol_particles(mytime)
            if (lag_params%write_bubbles_stats) call s_calculate_lag_particle_stats()
            if (lag_params%write_bubbles) then
                ! $:GPU_UPDATE(host='[gas_p,gas_mv,particle_rad,intfc_vel]')
                $:GPU_UPDATE(host='[particle_rad]')
                call s_write_lag_particle_evol(mytime)
            end if

        elseif (time_stepper == 2) then ! 2nd order TVD RK
            if (stage == 1) then

                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    !u{1} = u{n} +  dt * RHS{n}
                    particle_rad(k, 2) = particle_rad(k, 1) + dt*intfc_draddt(k, 1)
                    ! intfc_vel(k, 2) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                    ! gas_p(k, 2) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                    ! gas_mv(k, 2) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                    if (moving_lag_particles) then
                        mtn_posPrev(k, 1:3, 2) = mtn_pos(k, 1:3, 1)
                        mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                        mtn_velPrev(k, 1:3, 2) = mtn_vel(k, 1:3, 1)
                        mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf)

            elseif (stage == 2) then

                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    !u{1} = u{n} + (1/2) * dt * (RHS{n} + RHS{1})
                    particle_rad(k, 1) = particle_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/2._wp
                    ! intfc_vel(k, 1) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/2._wp
                    ! gas_p(k, 1) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/2._wp
                    ! gas_mv(k, 1) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/2._wp
                    if (moving_lag_particles) then
                        mtn_posPrev(k, 1:3, 1) = mtn_pos(k, 1:3, 2)
                        mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/2._wp
                        mtn_velPrev(k, 1:3, 1) = mtn_vel(k, 1:3, 2)
                        mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/2._wp
                    end if

                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_transfer_data_to_tmp_particles()
                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf)
                if (lag_params%write_void_evol) call s_write_void_evol_particles(mytime)
                if (lag_params%write_bubbles_stats) call s_calculate_lag_particle_stats()
                if (lag_params%write_bubbles) then
                    $:GPU_UPDATE(host='[particle_rad]')
                    call s_write_lag_particle_evol(mytime)
                end if

            end if

        elseif (time_stepper == 3) then ! 3rd order TVD RK
            if (stage == 1) then

                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    !u{1} = u{n} +  dt * RHS{n}
                    particle_rad(k, 2) = particle_rad(k, 1) + dt*intfc_draddt(k, 1)
                    ! intfc_vel(k, 2) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                    ! gas_p(k, 2) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                    ! gas_mv(k, 2) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                    if (moving_lag_particles) then
                        mtn_posPrev(k, 1:3, 2) = mtn_pos(k, 1:3, 1)
                        mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                        mtn_velPrev(k, 1:3, 2) = mtn_vel(k, 1:3, 1)
                        mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf)

            elseif (stage == 2) then

                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    !u{2} = u{n} + (1/4) * dt * [RHS{n} + RHS{1}]
                    particle_rad(k, 2) = particle_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/4._wp
                    ! intfc_vel(k, 2) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/4._wp
                    ! gas_p(k, 2) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/4._wp
                    ! gas_mv(k, 2) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/4._wp
                    if (moving_lag_particles) then
                        mtn_posPrev(k, 1:3, 2) = mtn_pos(k, 1:3, 2)
                        mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/4._wp
                        mtn_velPrev(k, 1:3, 2) = mtn_vel(k, 1:3, 2)
                        mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/4._wp
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf)

            elseif (stage == 3) then

                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    !u{n+1} = u{n} + (2/3) * dt * [(1/4)* RHS{n} + (1/4)* RHS{1} + RHS{2}]
                    particle_rad(k, 1) = particle_rad(k, 1) + (2._wp/3._wp)*dt*(intfc_draddt(k, 1)/4._wp + intfc_draddt(k, 2)/4._wp + intfc_draddt(k, 3))
                    ! intfc_vel(k, 1) = intfc_vel(k, 1) + (2._wp/3._wp)*dt*(intfc_dveldt(k, 1)/4._wp + intfc_dveldt(k, 2)/4._wp + intfc_dveldt(k, 3))
                    ! gas_p(k, 1) = gas_p(k, 1) + (2._wp/3._wp)*dt*(gas_dpdt(k, 1)/4._wp + gas_dpdt(k, 2)/4._wp + gas_dpdt(k, 3))
                    ! gas_mv(k, 1) = gas_mv(k, 1) + (2._wp/3._wp)*dt*(gas_dmvdt(k, 1)/4._wp + gas_dmvdt(k, 2)/4._wp + gas_dmvdt(k, 3))
                    if (moving_lag_particles) then
                        mtn_posPrev(k, 1:3, 1) = mtn_pos(k, 1:3, 2)
                        mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + (2._wp/3._wp)*dt*(mtn_dposdt(k, 1:3, 1)/4._wp + mtn_dposdt(k, 1:3, 2)/4._wp + mtn_dposdt(k, 1:3, 3))
                        mtn_velPrev(k, 1:3, 1) = mtn_vel(k, 1:3, 2)
                        mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + (2._wp/3._wp)*dt*(mtn_dveldt(k, 1:3, 1)/4._wp + mtn_dveldt(k, 1:3, 2)/4._wp + mtn_dveldt(k, 1:3, 3))
                    end if

                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_transfer_data_to_tmp_particles()
                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf)
                if (lag_params%write_void_evol) call s_write_void_evol_particles(mytime)
                if (lag_params%write_bubbles_stats) call s_calculate_lag_particle_stats()
                if (lag_params%write_bubbles) then
                    $:GPU_UPDATE(host='[particle_mass,particle_rad]')
                    call s_write_lag_particle_evol(mytime)
                end if

            end if

        end if

    end subroutine s_update_lagrange_particles_tdv_rk

    !> This subroutine enforces reflective and wall boundary conditions for EL particles
        !! @param dest Destination for the bubble position update
    impure subroutine s_enforce_EL_particles_boundary_conditions(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer :: k, i, q
        integer :: patch_id, newBubs, new_idx
        real(wp) :: offset
        integer, dimension(3) :: cell

        $:GPU_PARALLEL_LOOP(private='[k, cell]')
        do k = 1, n_el_particles_loc
            keep_bubble(k) = 1
            wrap_bubble_loc(k, :) = 0
            wrap_bubble_dir(k, :) = 0

            ! Relocate bubbles at solid boundaries and delete bubbles that leave
            ! buffer regions
            if (any(bc_x%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                .and. mtn_pos(k, 1, 2) < x_cb(-1) + particle_rad(k, 2)) then
                mtn_pos(k, 1, 2) = x_cb(-1) + particle_rad(k, 2)
            elseif (any(bc_x%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                    .and. mtn_pos(k, 1, 2) > x_cb(m) - particle_rad(k, 2)) then
                mtn_pos(k, 1, 2) = x_cb(m) - particle_rad(k, 2)
            elseif (bc_x%beg == BC_PERIODIC .and. mtn_pos(k, 1, 2) < pcomm_coords(1)%beg .and. &
                    mtn_posPrev(k, 1, 2) > pcomm_coords(1)%beg) then
                wrap_bubble_dir(k, 1) = 1
                wrap_bubble_loc(k, 1) = -1
            elseif (bc_x%end == BC_PERIODIC .and. mtn_pos(k, 1, 2) > pcomm_coords(1)%end .and. &
                    mtn_posPrev(k, 1, 2) < pcomm_coords(1)%end) then
                wrap_bubble_dir(k, 1) = 1
                wrap_bubble_loc(k, 1) = 1
            elseif (mtn_pos(k, 1, 2) >= x_cb(m + mapcells + 1)) then
                keep_bubble(k) = 0
            elseif (mtn_pos(k, 1, 2) < x_cb(-mapcells - 2)) then
                keep_bubble(k) = 0
            end if

            if (any(bc_y%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                .and. mtn_pos(k, 2, 2) < y_cb(-1) + particle_rad(k, 2)) then
                mtn_pos(k, 2, 2) = y_cb(-1) + particle_rad(k, 2)
            else if (any(bc_y%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                     .and. mtn_pos(k, 2, 2) > y_cb(n) - particle_rad(k, 2)) then
                mtn_pos(k, 2, 2) = y_cb(n) - particle_rad(k, 2)
            elseif (bc_y%beg == BC_PERIODIC .and. mtn_pos(k, 2, 2) < pcomm_coords(2)%beg .and. &
                    mtn_posPrev(k, 2, 2) > pcomm_coords(2)%beg) then
                wrap_bubble_dir(k, 2) = 1
                wrap_bubble_loc(k, 2) = -1
            elseif (bc_y%end == BC_PERIODIC .and. mtn_pos(k, 2, 2) > pcomm_coords(2)%end .and. &
                    mtn_posPrev(k, 2, 2) < pcomm_coords(2)%end) then
                wrap_bubble_dir(k, 2) = 1
                wrap_bubble_loc(k, 2) = 1
            elseif (mtn_pos(k, 2, 2) >= y_cb(n + mapcells + 1)) then
                keep_bubble(k) = 0
            elseif (mtn_pos(k, 2, 2) < y_cb(-mapcells - 2)) then
                keep_bubble(k) = 0
            end if

            if (p > 0) then
                if (any(bc_z%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                    .and. mtn_pos(k, 3, 2) < z_cb(-1) + particle_rad(k, 2)) then
                    mtn_pos(k, 3, 2) = z_cb(-1) + particle_rad(k, 2)
                else if (any(bc_z%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                         .and. mtn_pos(k, 3, 2) > z_cb(p) - particle_rad(k, 2)) then
                    mtn_pos(k, 3, 2) = z_cb(p) - particle_rad(k, 2)
                elseif (bc_z%beg == BC_PERIODIC .and. mtn_pos(k, 3, 2) < pcomm_coords(3)%beg .and. &
                        mtn_posPrev(k, 3, 2) > pcomm_coords(3)%beg) then
                    wrap_bubble_dir(k, 3) = 1
                    wrap_bubble_loc(k, 3) = -1
                elseif (bc_z%end == BC_PERIODIC .and. mtn_pos(k, 3, 2) > pcomm_coords(3)%end .and. &
                        mtn_posPrev(k, 3, 2) < pcomm_coords(3)%end) then
                    wrap_bubble_dir(k, 3) = 1
                    wrap_bubble_loc(k, 3) = 1
                elseif (mtn_pos(k, 3, 2) >= z_cb(p + mapCells + 1)) then
                    keep_bubble(k) = 0
                elseif (mtn_pos(k, 3, 2) < z_cb(-mapCells - 2)) then
                    keep_bubble(k) = 0
                end if
            end if

            if (keep_bubble(k) == 1) then
                ! Remove bubbles that are no longer in a liquid
                cell = fd_number - buff_size
                call s_locate_cell(mtn_pos(k, 1:3, 2), cell, mtn_s(k, 1:3, 2))

                if (q_prim_vf(advxb)%sf(cell(1), cell(2), cell(3)) < (1._wp - lag_params%valmaxvoid)) then
                    keep_bubble(k) = 0
                end if

                ! Move bubbles back to surface of IB
                if (ib) then
                    cell = fd_number - buff_size
                    call s_locate_cell(mtn_pos(k, 1:3, 2), cell, mtn_s(k, 1:3, 2))

                    if (ib_markers%sf(cell(1), cell(2), cell(3)) /= 0) then
                        patch_id = ib_markers%sf(cell(1), cell(2), cell(3))

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_dims
                            mtn_pos(k, i, 2) = mtn_pos(k, i, 2) - &
                                               levelset_norm%sf(cell(1), cell(2), cell(3), patch_id, i) &
                                               *levelset%sf(cell(1), cell(2), cell(3), patch_id)
                        end do

                        cell = fd_number - buff_size
                        call s_locate_cell(mtn_pos(k, 1:3, 2), cell, mtn_s(k, 1:3, 2))
                    end if
                end if
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxStartRange("LAG-BC")
        call nvtxStartRange("LAG-BC-DEV2HOST")
        $:GPU_UPDATE(host='[particle_R0, Rmax_stats, Rmin_stats, particle_mass, gas_betaT, &
            & gas_betaC, lag_id, particle_rad, &
            & mtn_pos, mtn_posPrev, mtn_vel, mtn_velPrev, mtn_s, intfc_draddt, &
            & mtn_dposdt, mtn_dveldt, keep_bubble, n_el_particles_loc, &
            & wrap_bubble_dir, wrap_bubble_loc]')
        call nvtxEndRange
        
        if (n_el_particles_loc > 0) then

            newBubs = 0
            do k = 1, n_el_particles_loc
                if (keep_bubble(k) == 1) then
                    newBubs = newBubs + 1
                    if (newBubs /= k) then
                        call s_copy_lag_particle(newBubs, k)
                        wrap_bubble_dir(newBubs,:) = wrap_bubble_dir(k,:)
                        wrap_bubble_loc(newBubs,:) = wrap_bubble_loc(k,:)
                    end if
                end if
            end do

            n_el_particles_loc = newBubs

            ! Handle periodic wrapping of bubbles on same processor
            newBubs = 0
            do k = 1, n_el_particles_loc
                if (any(wrap_bubble_dir(k, :) == 1)) then
                    newBubs = newBubs + 1
                    new_idx = n_el_particles_loc + newBubs
                    call s_copy_lag_particle(new_idx, k)
                    do i = 1, num_dims
                        if (wrap_bubble_dir(k, i) == 1) then
                            offset = glb_bounds(i)%end - glb_bounds(i)%beg
                            if (wrap_bubble_loc(k, i) == 1) then
                                do q = 1, 2
                                    mtn_pos(new_idx, i, q) = mtn_pos(new_idx, i, q) - offset
                                    mtn_posPrev(new_idx, i, q) = mtn_posPrev(new_idx, i, q) - offset
                                end do
                            else if (wrap_bubble_loc(k, i) == -1) then
                                do q = 1, 2
                                    mtn_pos(new_idx, i, q) = mtn_pos(new_idx, i, q) + offset
                                    mtn_posPrev(new_idx, i, q) = mtn_posPrev(new_idx, i, q) + offset
                                end do
                            end if
                        end if
                    end do
                end if
            end do
            n_el_particles_loc = n_el_particles_loc + newBubs
        end if

        ! Handle MPI transfer of particle going to another processor's local domain
        if (num_procs > 1) then
            call nvtxStartRange("LAG-BC-TRANSFER-LIST")
            call s_add_particles_to_transfer_list(n_el_particles_loc, mtn_pos(:, :, 2), mtn_posPrev(:, :, 2))
            call nvtxEndRange

            call nvtxStartRange("LAG-BC-SENDRECV")
            ! call s_mpi_sendrecv_particles(particle_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, &
            !                               gas_betaC, bub_dphidt, lag_id, gas_p, gas_mv, &
            !                               particle_rad, intfc_vel, mtn_pos, mtn_posPrev, mtn_vel, &
            !                               mtn_s, intfc_draddt, intfc_dveldt, gas_dpdt, &
            !                               gas_dmvdt, mtn_dposdt, mtn_dveldt, lag_num_ts, n_el_particles_loc, &
            !                               2)

            call s_mpi_sendrecv_solid_particles(particle_R0, Rmax_stats, Rmin_stats, particle_mass, gas_betaT, &
                                                gas_betaC, lag_id, &
                                                particle_rad, mtn_pos, mtn_posPrev, mtn_vel, mtn_velPrev, &
                                                mtn_s, intfc_draddt, mtn_dposdt, mtn_dveldt, lag_num_ts, n_el_particles_loc, &
                                                2)
            call nvtxEndRange
        end if

        call nvtxStartRange("LAG-BC-HOST2DEV")
        ! $:GPU_UPDATE(device='[particle_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, &
        !     & gas_betaC, bub_dphidt, lag_id, gas_p, gas_mv, particle_rad, intfc_vel, &
        !     & mtn_pos, mtn_posPrev, mtn_vel, mtn_s, intfc_draddt, intfc_dveldt, &
        !     & gas_dpdt, gas_dmvdt, mtn_dposdt, mtn_dveldt, n_el_particles_loc]')
        $:GPU_UPDATE(host='[particle_R0, Rmax_stats, Rmin_stats, particle_mass, gas_betaT, &
            & gas_betaC, lag_id, particle_rad, &
            & mtn_pos, mtn_posPrev, mtn_vel, mtn_velPrev, mtn_s, intfc_draddt, &
            & mtn_dposdt, mtn_dveldt, n_el_particles_loc]')
        call nvtxEndRange
        call nvtxEndRange

        ! $:GPU_PARALLEL_LOOP(private='[cell]')
        ! do k = 1, n_el_particles_loc
        !     cell = fd_number - buff_size
        !     call s_locate_cell(mtn_pos(k, 1:3, 2), cell, mtn_s(k, 1:3, 2))
        ! end do

        ! Update void fraction and communicate buffers
        call s_smear_particle_sources()

    end subroutine s_enforce_EL_particles_boundary_conditions

    !> This subroutine returns the computational coordinate of the cell for the given position.
          !! @param pos Input coordinates
          !! @param cell Computational coordinate of the cell
          !! @param scoord Calculated particle coordinates
    subroutine s_locate_cell(pos, cell, scoord)
        $:GPU_ROUTINE(function_name='s_locate_cell',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), dimension(3), intent(in) :: pos
        real(wp), dimension(3), intent(out) :: scoord
        integer, dimension(3), intent(inout) :: cell

        integer :: i

        do while (pos(1) < x_cb(cell(1) - 1))
            cell(1) = cell(1) - 1
        end do

        do while (pos(1) >= x_cb(cell(1)))
            cell(1) = cell(1) + 1
        end do

        do while (pos(2) < y_cb(cell(2) - 1))
            cell(2) = cell(2) - 1
        end do

        do while (pos(2) >= y_cb(cell(2)))
            cell(2) = cell(2) + 1
        end do

        if (p > 0) then
            do while (pos(3) < z_cb(cell(3) - 1))
                cell(3) = cell(3) - 1
            end do
            do while (pos(3) >= z_cb(cell(3)))
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
    impure subroutine s_transfer_data_to_tmp_particles()

        integer :: k

        $:GPU_PARALLEL_LOOP(private='[k]')
        do k = 1, n_el_particles_loc
            ! gas_p(k, 2) = gas_p(k, 1)
            ! gas_mv(k, 2) = gas_mv(k, 1)
            particle_rad(k, 2) = particle_rad(k, 1)
            ! intfc_vel(k, 2) = intfc_vel(k, 1)
            mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1)
            mtn_posPrev(k, 1:3, 2) = mtn_posPrev(k, 1:3, 1)
            mtn_velPrev(k, 1:3, 2) = mtn_velPrev(k, 1:3, 1)
            mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1)
            mtn_s(k, 1:3, 2) = mtn_s(k, 1:3, 1)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_transfer_data_to_tmp_particles

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
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size - fd_number)) .and. &
                                  (pos_part(1) >= x_cb(fd_number - buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size - fd_number)) .and. &
                                  (pos_part(2) >= y_cb(fd_number - buff_size - 1)) .and. &
                                  (pos_part(3) < lag_params%charwidth/2._wp) .and. (pos_part(3) > -lag_params%charwidth/2._wp))
        else
            ! cyl_coord
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size - fd_number)) .and. &
                                  (pos_part(1) >= x_cb(fd_number - buff_size - 1)) .and. &
                                  (abs(pos_part(2)) < y_cb(n + buff_size - fd_number)) .and. &
                                  (abs(pos_part(2)) >= max(y_cb(fd_number - buff_size - 1), 0._wp)))
        end if

        ! 3D
        if (p > 1) then
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size - fd_number)) .and. &
                                  (pos_part(1) >= x_cb(fd_number - buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size - fd_number)) .and. &
                                  (pos_part(2) >= y_cb(fd_number - buff_size - 1)) .and. &
                                  (pos_part(3) < z_cb(p + buff_size - fd_number)) .and. &
                                  (pos_part(3) >= z_cb(fd_number - buff_size - 1)))
        end if

        ! For symmetric and wall boundary condition
        if (any(bc_x%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then
            particle_in_domain = (particle_in_domain .and. (pos_part(1) >= x_cb(-1)))
        end if
        if (any(bc_x%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then
            particle_in_domain = (particle_in_domain .and. (pos_part(1) < x_cb(m)))
        end if
        if (any(bc_y%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. (.not. cyl_coord)) then
            particle_in_domain = (particle_in_domain .and. (pos_part(2) >= y_cb(-1)))
        end if
        if (any(bc_y%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. (.not. cyl_coord)) then
            particle_in_domain = (particle_in_domain .and. (pos_part(2) < y_cb(n)))
        end if
        if (p > 0) then
            if (any(bc_z%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then
                particle_in_domain = (particle_in_domain .and. (pos_part(3) >= z_cb(-1)))
            end if
            if (any(bc_z%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then
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

        real(stp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:), intent(inout) :: q, dq
        integer, intent(in) :: dir

        integer :: i, j, k

        if (dir == 1) then
            ! Gradient in x dir.
            $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        dq(i, j, k) = q(i, j, k)*(dx(i + 1) - dx(i - 1)) &
                                      + q(i + 1, j, k)*(dx(i) + dx(i - 1)) &
                                      - q(i - 1, j, k)*(dx(i) + dx(i + 1))
                        dq(i, j, k) = dq(i, j, k)/ &
                                      ((dx(i) + dx(i - 1))*(dx(i) + dx(i + 1)))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (dir == 2) then
            ! Gradient in y dir.
            $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        dq(i, j, k) = q(i, j, k)*(dy(j + 1) - dy(j - 1)) &
                                      + q(i, j + 1, k)*(dy(j) + dy(j - 1)) &
                                      - q(i, j - 1, k)*(dy(j) + dy(j + 1))
                        dq(i, j, k) = dq(i, j, k)/ &
                                      ((dy(j) + dy(j - 1))*(dy(j) + dy(j + 1)))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (dir == 3) then
            ! Gradient in z dir.
            $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        dq(i, j, k) = q(i, j, k)*(dz(k + 1) - dz(k - 1)) &
                                      + q(i, j, k + 1)*(dz(k) + dz(k - 1)) &
                                      - q(i, j, k - 1)*(dz(k) + dz(k + 1))
                        dq(i, j, k) = dq(i, j, k)/ &
                                      ((dz(k) + dz(k - 1))*(dz(k) + dz(k + 1)))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_gradient_dir

    impure subroutine s_open_lag_bubble_evol

        character(LEN=path_len + 2*name_len) :: file_loc
        logical file_exist
        character(LEN=25) :: FMT

        write (file_loc, '(A,I0,A)') 'lag_bubble_evol_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)
        call my_inquire(trim(file_loc), file_exist)

        if (precision == 1) then
            FMT = "(A16,A14,8A16)"
        else
            FMT = "(A24,A14,8A24)"
        end if

        if (.not. file_exist) then
            open (LAG_EVOL_ID, FILE=trim(file_loc), FORM='formatted', position='rewind')
            write (LAG_EVOL_ID, FMT) 'currentTime', 'particleID', 'x', 'y', 'z', &
                'coreVaporMass', 'coreVaporConcentration', 'radius', 'interfaceVelocity', &
                'corePressure'
        else
            open (LAG_EVOL_ID, FILE=trim(file_loc), FORM='formatted', position='append')
        end if

    end subroutine s_open_lag_bubble_evol

    !> Subroutine that writes on each time step the changes of the lagrangian bubbles.
        !!  @param q_time Current time
    impure subroutine s_write_lag_particle_evol(qtime)

        real(wp), intent(in) :: qtime
        integer :: k, ios
        character(LEN=25) :: FMT

        character(LEN=path_len + 2*name_len) :: file_loc, path
        logical :: file_exist

        if (precision == 1) then
            FMT = "(F16.8,I14,8F16.8)"
        else
            FMT = "(F24.16,I14,8F24.16)"
        end if

        ! Cycle through list
        do k = 1, n_el_particles_loc
            write (LAG_EVOL_ID, FMT) &
                qtime, &
                lag_id(k, 1), &
                mtn_pos(k, 1, 1), &
                mtn_pos(k, 2, 1), &
                mtn_pos(k, 3, 1), &
                ! gas_mv(k, 1), &
                ! gas_mv(k, 1)/(gas_mv(k, 1) + gas_mg(k)), &
                particle_rad(k, 1)
            ! intfc_vel(k, 1), &
            ! gas_p(k, 1)
        end do

    end subroutine s_write_lag_particle_evol

    impure subroutine s_close_lag_particle_evol

        close (LAG_EVOL_ID)

    end subroutine s_close_lag_particle_evol

    subroutine s_open_void_evol

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        if (proc_rank == 0) then
            write (file_loc, '(A)') 'voidfraction.dat'
            file_loc = trim(case_dir)//'/D/'//trim(file_loc)
            call my_inquire(trim(file_loc), file_exist)
            if (.not. file_exist) then
                open (LAG_VOID_ID, FILE=trim(file_loc), FORM='formatted', position='rewind')
                !write (12, *) 'currentTime, averageVoidFraction, ', &
                !    'maximumVoidFraction, totalParticlesVolume'
                !write (12, *) 'The averageVoidFraction value does ', &
                !    'not reflect the real void fraction in the cloud since the ', &
                !    'cells which do not have bubbles are not accounted'
            else
                open (LAG_VOID_ID, FILE=trim(file_loc), FORM='formatted', position='append')
            end if
        end if

    end subroutine s_open_void_evol

    !>  Subroutine that writes some useful statistics related to the volume fraction
            !!       of the particles (void fraction) in the computatioational domain
            !!       on each time step.
            !!  @param q_time Current time
    impure subroutine s_write_void_evol_particles(qtime)

        real(wp), intent(in) :: qtime
        real(wp) :: volcell, voltot
        real(wp) :: lag_void_max, lag_void_avg, lag_vol
        real(wp) :: void_max_glb, void_avg_glb, vol_glb

        integer :: i, j, k

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        lag_void_max = 0._wp
        lag_void_avg = 0._wp
        lag_vol = 0._wp
        $:GPU_PARALLEL_LOOP(private='[volcell]', collapse=3, reduction='[[lag_vol, lag_void_avg], [lag_void_max]]', reductionOp='[+, MAX]', copy='[lag_vol, lag_void_avg, lag_void_max]')
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    lag_void_max = max(lag_void_max, 1._wp - q_particles(1)%sf(i, j, k))
                    call s_get_char_vol(i, j, k, volcell)
                    if ((1._wp - q_particles(1)%sf(i, j, k)) > 5.0d-11) then
                        lag_void_avg = lag_void_avg + (1._wp - q_particles(1)%sf(i, j, k))*volcell
                        lag_vol = lag_vol + volcell
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

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
            write (LAG_VOID_ID, '(6X,4e24.8)') &
                qtime, &
                lag_void_avg, &
                lag_void_max, &
                voltot
        end if

    end subroutine s_write_void_evol_particles

    subroutine s_close_void_evol

        if (proc_rank == 0) close (LAG_VOID_ID)

    end subroutine s_close_void_evol

    !>  Subroutine that writes the restarting files for the particles in the lagrangian solver.
        !!  @param t_step Current time step
    impure subroutine s_write_restart_lag_particles(t_step)

        ! Generic string used to store the address of a particular file
        integer, intent(in) :: t_step

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist
        integer :: part_id, tot_part
        integer :: i, k

#ifdef MFC_MPI
        ! For Parallel I/O
        integer :: ifile, ierr
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer :: view
        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer, dimension(num_procs) :: part_order, part_ord_mpi
        integer, dimension(num_procs) :: proc_particle_counts
        real(wp), dimension(1:1, 1:lag_io_vars) :: dummy
        dummy = 0._wp

        part_id = 0._wp
        if (n_el_particles_loc /= 0) then
            do k = 1, n_el_particles_loc
                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then
                    part_id = part_id + 1
                end if
            end do
        end if

        if (.not. parallel_io) return

        lsizes(1) = part_id
        lsizes(2) = lag_io_vars

        ! Total number of particles
        call MPI_ALLREDUCE(part_id, tot_part, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        call MPI_ALLGATHER(part_id, 1, MPI_INTEGER, proc_particle_counts, 1, MPI_INTEGER, &
                           MPI_COMM_WORLD, ierr)

        ! Calculate starting index for this processor's particles
        call MPI_EXSCAN(lsizes(1), start_idx_part(1), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (proc_rank == 0) start_idx_part(1) = 0
        start_idx_part(2) = 0

        gsizes(1) = tot_part
        gsizes(2) = lag_io_vars

        write (file_loc, '(A,I0,A)') 'lag_bubbles_', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)

        ! Clean up existing file
        if (proc_rank == 0) then
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, &
                               ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Write header using MPI I/O for consistency
            call MPI_FILE_WRITE(ifile, tot_part, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_WRITE(ifile, mytime, 1, mpi_p, status, ierr)
            call MPI_FILE_WRITE(ifile, dt, 1, mpi_p, status, ierr)
            call MPI_FILE_WRITE(ifile, num_procs, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_WRITE(ifile, proc_particle_counts, num_procs, MPI_INTEGER, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if (part_id > 0) then
            allocate (MPI_IO_DATA_lag_bubbles(max(1, part_id), 1:lag_io_vars))

            i = 1
            do k = 1, n_el_particles_loc
                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then
                    MPI_IO_DATA_lag_bubbles(i, 1) = real(lag_id(k, 1))
                    MPI_IO_DATA_lag_bubbles(i, 2:4) = mtn_pos(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 5:7) = mtn_posPrev(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 8:10) = mtn_vel(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 11) = particle_rad(k, 1)
                    ! MPI_IO_DATA_lag_bubbles(i, 12) = intfc_vel(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 13) = particle_R0(k)
                    MPI_IO_DATA_lag_bubbles(i, 14) = Rmax_stats(k)
                    MPI_IO_DATA_lag_bubbles(i, 15) = Rmin_stats(k)
                    ! MPI_IO_DATA_lag_bubbles(i, 16) = bub_dphidt(k)
                    ! MPI_IO_DATA_lag_bubbles(i, 17) = gas_p(k, 1)
                    ! MPI_IO_DATA_lag_bubbles(i, 18) = gas_mv(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 19) = particle_mass(k)
                    MPI_IO_DATA_lag_bubbles(i, 20) = gas_betaT(k)
                    MPI_IO_DATA_lag_bubbles(i, 21) = gas_betaC(k)
                    i = i + 1
                end if
            end do

            call MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                          MPI_ORDER_FORTRAN, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, &
                               ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Skip header (written by rank 0)
            disp = int(sizeof(tot_part) + 2*sizeof(mytime) + sizeof(num_procs) + &
                       num_procs*sizeof(proc_particle_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA_lag_bubbles, &
                                    lag_io_vars*part_id, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)

            deallocate (MPI_IO_DATA_lag_bubbles)

        else
            call MPI_TYPE_CONTIGUOUS(0, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, &
                               ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Skip header (written by rank 0)
            disp = int(sizeof(tot_part) + 2*sizeof(mytime) + sizeof(num_procs) + &
                       num_procs*sizeof(proc_particle_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_WRITE_ALL(ifile, dummy, 0, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

#endif

    end subroutine s_write_restart_lag_particles

    !>  This procedure calculates the maximum and minimum radius of each bubble.
    subroutine s_calculate_lag_particle_stats()

        integer :: k

        $:GPU_PARALLEL_LOOP(private='[k]', reduction='[[Rmax_glb], [Rmin_glb]]', &
            & reductionOp='[MAX, MIN]', copy='[Rmax_glb,Rmin_glb]')
        do k = 1, n_el_particles_loc
            Rmax_glb = max(Rmax_glb, particle_rad(k, 1)/particle_R0(k))
            Rmin_glb = min(Rmin_glb, particle_rad(k, 1)/particle_R0(k))
            Rmax_stats(k) = max(Rmax_stats(k), particle_rad(k, 1)/particle_R0(k))
            Rmin_stats(k) = min(Rmin_stats(k), particle_rad(k, 1)/particle_R0(k))
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_calculate_lag_particle_stats

    impure subroutine s_open_lag_particle_stats()

        character(LEN=path_len + 2*name_len) :: file_loc
        character(LEN=20) :: FMT
        logical :: file_exist

        write (file_loc, '(A,I0,A)') 'stats_lag_bubbles_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)
        call my_inquire(trim(file_loc), file_exist)

        if (precision == 1) then
            FMT = "(A10,A14,5A16)"
        else
            FMT = "(A10,A14,5A24)"
        end if

        if (.not. file_exist) then
            open (LAG_STATS_ID, FILE=trim(file_loc), FORM='formatted', position='rewind')
            write (LAG_STATS_ID, *) 'proc_rank, particleID, x, y, z, Rmax_glb, Rmin_glb'
        else
            open (LAG_STATS_ID, FILE=trim(file_loc), FORM='formatted', position='append')
        end if

    end subroutine s_open_lag_particle_stats

    !>  Subroutine that writes the maximum and minimum radius of each bubble.
    impure subroutine s_write_lag_particle_stats()

        integer :: k
        character(LEN=path_len + 2*name_len) :: file_loc
        character(LEN=20) :: FMT

        $:GPU_UPDATE(host='[Rmax_glb,Rmin_glb]')

        if (precision == 1) then
            FMT = "(I10,I14,5F16.8)"
        else
            FMT = "(I10,I14,5F24.16)"
        end if

        do k = 1, n_el_particles_loc
            write (LAG_STATS_ID, FMT) &
                proc_rank, &
                lag_id(k, 1), &
                mtn_pos(k, 1, 1), &
                mtn_pos(k, 2, 1), &
                mtn_pos(k, 3, 1), &
                Rmax_stats(k), &
                Rmin_stats(k)
        end do

    end subroutine s_write_lag_particle_stats

    subroutine s_close_lag_particle_stats

        close (LAG_STATS_ID)

    end subroutine s_close_lag_particle_stats

    !> The purpose of this subroutine is to remove one specific particle if dt is too small.
          !! @param part_id Particle id
    impure subroutine s_copy_lag_particle(dest, src)

        integer, intent(in) :: src, dest

        particle_R0(dest) = particle_R0(src)
        Rmax_stats(dest) = Rmax_stats(src)
        Rmin_stats(dest) = Rmin_stats(src)
        particle_mass(dest) = particle_mass(src)
        gas_betaT(dest) = gas_betaT(src)
        gas_betaC(dest) = gas_betaC(src)
        ! bub_dphidt(dest) = bub_dphidt(src)
        lag_id(dest, 1) = lag_id(src, 1)
        ! gas_p(dest, 1:2) = gas_p(src, 1:2)
        ! gas_mv(dest, 1:2) = gas_mv(src, 1:2)
        particle_rad(dest, 1:2) = particle_rad(src, 1:2)
        ! intfc_vel(dest, 1:2) = intfc_vel(src, 1:2)
        mtn_vel(dest, 1:3, 1:2) = mtn_vel(src, 1:3, 1:2)
        mtn_s(dest, 1:3, 1:2) = mtn_s(src, 1:3, 1:2)
        mtn_pos(dest, 1:3, 1:2) = mtn_pos(src, 1:3, 1:2)
        mtn_posPrev(dest, 1:3, 1:2) = mtn_posPrev(src, 1:3, 1:2)
        mtn_velPrev(dest, 1:3, 1:2) = mtn_velPrev(src, 1:3, 1:2)
        intfc_draddt(dest, 1:lag_num_ts) = intfc_draddt(src, 1:lag_num_ts)
        ! intfc_dveldt(dest, 1:lag_num_ts) = intfc_dveldt(src, 1:lag_num_ts)
        ! gas_dpdt(dest, 1:lag_num_ts) = gas_dpdt(src, 1:lag_num_ts)
        ! gas_dmvdt(dest, 1:lag_num_ts) = gas_dmvdt(src, 1:lag_num_ts)
        mtn_dposdt(dest, 1:3, 1:lag_num_ts) = mtn_dposdt(src, 1:3, 1:lag_num_ts)
        mtn_dveldt(dest, 1:3, 1:lag_num_ts) = mtn_dveldt(src, 1:3, 1:lag_num_ts)

    end subroutine s_copy_lag_particle

    !> The purpose of this subroutine is to deallocate variables
    impure subroutine s_finalize_particle_lagrangian_solver()

        integer :: i

        if (lag_params%write_void_evol) call s_close_void_evol
        if (lag_params%write_bubbles) call s_close_lag_particle_evol()
        if (lag_params%write_bubbles_stats) call s_close_lag_particle_stats()

        do i = 1, q_particles_idx
            @:DEALLOCATE(q_particles(i)%sf)
        end do
        @:DEALLOCATE(q_particles)

        !Deallocating space
        @:DEALLOCATE(lag_id)
        @:DEALLOCATE(particle_R0)
        @:DEALLOCATE(Rmax_stats)
        @:DEALLOCATE(Rmin_stats)
        @:DEALLOCATE(particle_mass)
        @:DEALLOCATE(gas_betaT)
        @:DEALLOCATE(gas_betaC)
        ! @:DEALLOCATE(bub_dphidt)
        ! @:DEALLOCATE(gas_p)
        ! @:DEALLOCATE(gas_mv)
        @:DEALLOCATE(particle_rad)
        ! @:DEALLOCATE(intfc_vel)
        @:DEALLOCATE(mtn_pos)
        @:DEALLOCATE(mtn_posPrev)
        @:DEALLOCATE(mtn_velPrev)
        @:DEALLOCATE(mtn_vel)
        @:DEALLOCATE(mtn_s)
        @:DEALLOCATE(intfc_draddt)
        ! @:DEALLOCATE(intfc_dveldt)
        ! @:DEALLOCATE(gas_dpdt)
        ! @:DEALLOCATE(gas_dmvdt)
        @:DEALLOCATE(mtn_dposdt)
        @:DEALLOCATE(mtn_dveldt)
        @:DEALLOCATE(f_p)

    end subroutine s_finalize_particle_lagrangian_solver

end module m_particles_EL
