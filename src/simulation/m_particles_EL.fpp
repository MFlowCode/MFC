!>
!! @file m_particles_EL.fpp
!! @brief Contains module m_particles_EL

#:include 'macros.fpp'

!> @brief This module is used to to compute the volume-averaged particle model
module m_particles_EL

    use m_global_parameters             !< Definitions of the global parameters

    use m_mpi_proxy                     !< Message passing interface (MPI) module proxy

    use m_particles_EL_kernels            !< Definitions of the kernel functions

    use m_particles                       !< General particle procedures

    use m_variables_conversion          !< State variables type conversion procedures

    use m_compile_specific

    use m_boundary_common

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_sim_helpers

    use m_helper

    use m_mpi_common

    use m_ibm

    use m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables

    use m_muscl                !< Monotonic Upstream-centered (MUSCL)
                               !! schemes for conservation laws

    implicit none

    real(wp), save :: next_write_time = 0._wp

    integer, allocatable, dimension(:, :) :: lag_part_id                 !< Global and local IDs
    real(wp), allocatable, dimension(:) :: particle_R0            !< Initial particle radius
    real(wp), allocatable, dimension(:) :: Rmax_stats_part        !< Maximum radius
    real(wp), allocatable, dimension(:) :: Rmin_stats_part        !< Minimum radius
    $:GPU_DECLARE(create='[lag_part_id, particle_R0, Rmax_stats_part, Rmin_stats_part]')

    real(wp), allocatable, dimension(:) :: particle_mass            !< Particle Mass
    $:GPU_DECLARE(create='[particle_mass]')

    integer, allocatable, dimension(:) :: linked_list            !< particle cell linked list
    $:GPU_DECLARE(create='[linked_list]')

    !(nPart, 1 -> actual val or 2 -> temp val)
    real(wp), allocatable, dimension(:, :) :: particle_rad      !< Particle radius
    $:GPU_DECLARE(create='[particle_rad]')

    !(nPart, 1-> x or 2->y or 3 ->z, 1 -> actual or 2 -> temporal val)
    real(wp), allocatable, dimension(:, :, :) :: particle_pos     !< Particle's position
    real(wp), allocatable, dimension(:, :, :) :: particle_posPrev !< Particle's previous position
    real(wp), allocatable, dimension(:, :, :) :: particle_vel     !< Particle's velocity
    real(wp), allocatable, dimension(:, :, :) :: particle_s       !< Particle's computational cell position in real format
    $:GPU_DECLARE(create='[particle_pos, particle_posPrev, particle_vel, particle_s]')
    !(nPart, 1-> x or 2->y or 3 ->z, time-stage)
    real(wp), allocatable, dimension(:, :) :: particle_draddt   !< Time derivative of particle's radius
    real(wp), allocatable, dimension(:, :, :) :: particle_dposdt  !< Time derivative of the particle's position
    real(wp), allocatable, dimension(:, :, :) :: particle_dveldt  !< Time derivative of the particle's velocity
    $:GPU_DECLARE(create='[particle_draddt, particle_dposdt, particle_dveldt]')

    integer, private :: lag_num_ts                                  !<  Number of time stages in the time-stepping scheme
    $:GPU_DECLARE(create='[lag_num_ts]')

    real(wp) :: Rmax_glb, Rmin_glb       !< Maximum and minimum bubbe size in the local domain
    !< Projection of the lagrangian particles in the Eulerian framework
    type(scalar_field), dimension(:), allocatable :: q_particles
    integer :: q_particles_idx                     !< Size of the q vector field for particle cell (q)uantities
    integer, parameter :: alphaf_id = 1
    integer, parameter :: alphaupx_id = 2 !< x particle momentum index
    integer, parameter :: alphaupy_id = 3 !< y particle momentum index
    integer, parameter :: alphaupz_id = 4 !< z particle momentum index
    integer, parameter :: Smx_id = 5
    integer, parameter :: Smy_id = 6
    integer, parameter :: Smz_id = 7
    integer, parameter :: SE_id = 8

    type(scalar_field), dimension(:), allocatable :: field_vars !< For cell quantities (field gradients, etc.)
    integer, parameter :: dPx_id = 1 !< Spatial pressure gradient in x, y, and z
    integer, parameter :: dPy_id = 2
    integer, parameter :: dPz_id = 3
    integer, parameter :: dRhox_id = 4 !< Spatial density gradient in x, y, and z
    integer, parameter :: dRhoy_id = 5
    integer, parameter :: dRhoz_id = 6
    integer, parameter :: dRhoux_id = 7 !< Spatial momentum flux gradient in x, y, and z
    integer, parameter :: dRhouy_id = 8
    integer, parameter :: dRhouz_id = 9
    integer, parameter :: dalphafx_id = 10 !< Spatial volume fraction gradient in x, y, and z
    integer, parameter :: dalphafy_id = 11
    integer, parameter :: dalphafz_id = 12
    integer, parameter :: dalphap_upx_id = 13 !< Spatial particle momentum gradient in x, y, and z
    integer, parameter :: dalphap_upy_id = 14
    integer, parameter :: dalphap_upz_id = 15
    integer, parameter :: nField_vars = 15

    type(scalar_field), dimension(:), allocatable :: weights_x_interp !< For precomputing weights
    type(scalar_field), dimension(:), allocatable :: weights_y_interp !< For precomputing weights
    type(scalar_field), dimension(:), allocatable :: weights_z_interp !< For precomputing weights
    integer :: nWeights_interp

    type(scalar_field), dimension(:), allocatable :: weights_x_grad !< For precomputing weights
    type(scalar_field), dimension(:), allocatable :: weights_y_grad !< For precomputing weights
    type(scalar_field), dimension(:), allocatable :: weights_z_grad !< For precomputing weights
    integer :: nWeights_grad

    type(integer_field), dimension(:), allocatable :: particle_head !< particle heads at each cell

    $:GPU_DECLARE(create='[Rmax_glb,Rmin_glb,q_particles,q_particles_idx,field_vars]')
    $:GPU_DECLARE(create='[weights_x_interp,weights_y_interp,weights_z_interp,nWeights_interp]')
    $:GPU_DECLARE(create='[weights_x_grad,weights_y_grad,weights_z_grad,nWeights_grad]')

    real(wp), allocatable, dimension(:, :, :, :) :: myL_rsx_vf, myL_rsy_vf, myL_rsz_vf, myR_rsx_vf, myR_rsy_vf, myR_rsz_vf
    $:GPU_DECLARE(create='[myL_rsx_vf,myL_rsy_vf,myL_rsz_vf,myR_rsx_vf,myR_rsy_vf,myR_rsz_vf]')
    integer, parameter :: r_alphaf_id = 1 !< Reconstructed fluid volume fraction
    integer, parameter :: r_alphaupx_id = 2 !< Reconstructed particle velocity in x, y, and z
    integer, parameter :: r_alphaupy_id = 3
    integer, parameter :: r_alphaupz_id = 4
    integer, parameter :: r_rhouf_id = 5 !< Reconstructed fluid momentum
    integer, parameter :: nRecon = 5

    !Particle Source terms for fluid coupling
    real(wp), allocatable, dimension(:, :) :: f_p !< force on each particle
    $:GPU_DECLARE(create='[f_p]')

    integer, parameter :: LAG_EVOL_ID = 11 ! File id for lag_bubbles_evol_*.dat
    integer, parameter :: LAG_STATS_ID = 12 ! File id for stats_lag_bubbles_*.dat
    integer, parameter :: LAG_VOID_ID = 13 ! File id for voidfraction.dat

    integer, allocatable, dimension(:) :: keep_bubble
    integer, allocatable, dimension(:, :) :: wrap_bubble_loc, wrap_bubble_dir
    $:GPU_DECLARE(create='[keep_bubble]')
    $:GPU_DECLARE(create='[wrap_bubble_loc, wrap_bubble_dir]')

    integer :: ncc

    type(int_bounds_info) :: is1_loc, is2_loc, is3_loc
    $:GPU_DECLARE(create='[is1_loc,is2_loc,is3_loc]')

    type(int_bounds_info) :: iv_loc !< Vector field indical bounds
    $:GPU_DECLARE(create='[iv_loc]')

contains

    !> Initializes the lagrangian subgrid particle solver
        !! @param q_cons_vf Initial conservative variables
    impure subroutine s_initialize_particles_EL_module(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        integer :: nParticles_glb, i, j, k, nf, l, npts

        !PRIM TO CONS VARIABLES
        real(wp) :: dyn_pres, pi_inf, qv, gamma, pres, T
        real(wp) :: rhou, alpharhou, rho_f, alpharho
        real(wp), dimension(3) :: fluid_vel
        real(wp) :: rhoYks(1:num_species)
        integer :: save_count
        real(wp) :: qtime

        if (cfl_dt) then
            save_count = n_start
            qtime = n_start*t_save
        else
            save_count = t_step_start
            qtime = t_step_start*dt
        end if

        pi_inf = 0._wp
        qv = 0._wp
        gamma = gammas(1)

        ! Setting number of time-stages for selected time-stepping scheme
        lag_num_ts = time_stepper

        ! Allocate space for the Eulerian fields needed to map the effect of the particles
        if (lag_params%solver_approach == 1) then
            ! One-way coupling
            q_particles_idx = 1 !For tracking volume fraction
        elseif (lag_params%solver_approach == 2) then
            !Two-way coupling
            q_particles_idx = 8 !For tracking volume fraction(1), x-mom(2), y-mom(3), z-mom(4), and energy(5) sources, and alpha_p u_p (x(6),y(7),z(8))
        else
            call s_mpi_abort('Please check the lag_params%solver_approach input')
        end if

        nWeights_interp = lag_params%interpolation_order + 1
        nWeights_grad = fd_order + 1

        ncc = 3 !max(1,fd_order,lag_params%interpolation_order,buff_size)

        pcomm_coords(1)%beg = x_cb(-1 + ncc)
        pcomm_coords(1)%end = x_cb(m - ncc)
        $:GPU_UPDATE(device='[pcomm_coords(1)]')
        if (n > 0) then
            pcomm_coords(2)%beg = y_cb(-1 + ncc)
            pcomm_coords(2)%end = y_cb(n - ncc)
            $:GPU_UPDATE(device='[pcomm_coords(2)]')
            if (p > 0) then
                pcomm_coords(3)%beg = z_cb(-1 + ncc)
                pcomm_coords(3)%end = z_cb(p - ncc)
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

        @:ALLOCATE(particle_head(1))
        @:ALLOCATE(particle_head(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
            idwbuff(2)%beg:idwbuff(2)%end, &
            idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(particle_head(1))

        @:ALLOCATE(field_vars(1:nField_vars))
        do i = 1, nField_vars
            @:ALLOCATE(field_vars(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(field_vars(i))
        end do

        @:ALLOCATE(myL_rsx_vf(idwbuff(1)%beg:idwbuff(1)%end, &
            idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:nRecon))
        @:ALLOCATE(myR_rsx_vf(idwbuff(1)%beg:idwbuff(1)%end, &
            idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:nRecon))

        @:ALLOCATE(myL_rsy_vf(idwbuff(2)%beg:idwbuff(2)%end, &
            idwbuff(1)%beg:idwbuff(1)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:nRecon))
        @:ALLOCATE(myR_rsy_vf(idwbuff(2)%beg:idwbuff(2)%end, &
            idwbuff(1)%beg:idwbuff(1)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:nRecon))

        if (p > 0) then
            @:ALLOCATE(myL_rsz_vf(idwbuff(3)%beg:idwbuff(3)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, 1:nRecon))
            @:ALLOCATE(myR_rsz_vf(idwbuff(3)%beg:idwbuff(3)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, 1:nRecon))
        else
            @:ALLOCATE(myL_rsz_vf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:nRecon))
            @:ALLOCATE(myR_rsz_vf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:nRecon))
        end if

        @:ALLOCATE(weights_x_interp(1:nWeights_interp))
        do i = 1, nWeights_interp
            @:ALLOCATE(weights_x_interp(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,1:1,1:1))
            @:ACC_SETUP_SFs(weights_x_interp(i))
        end do

        @:ALLOCATE(weights_y_interp(1:nWeights_interp))
        do i = 1, nWeights_interp
            @:ALLOCATE(weights_y_interp(i)%sf(idwbuff(2)%beg:idwbuff(2)%end,1:1,1:1))
            @:ACC_SETUP_SFs(weights_y_interp(i))
        end do

        @:ALLOCATE(weights_z_interp(1:nWeights_interp))
        do i = 1, nWeights_interp
            @:ALLOCATE(weights_z_interp(i)%sf(idwbuff(3)%beg:idwbuff(3)%end,1:1,1:1))
            @:ACC_SETUP_SFs(weights_z_interp(i))
        end do

        @:ALLOCATE(weights_x_grad(1:nWeights_grad))
        do i = 1, nWeights_grad
            @:ALLOCATE(weights_x_grad(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,1:1,1:1))
            @:ACC_SETUP_SFs(weights_x_grad(i))
        end do

        @:ALLOCATE(weights_y_grad(1:nWeights_grad))
        do i = 1, nWeights_grad
            @:ALLOCATE(weights_y_grad(i)%sf(idwbuff(2)%beg:idwbuff(2)%end,1:1,1:1))
            @:ACC_SETUP_SFs(weights_y_grad(i))
        end do

        @:ALLOCATE(weights_z_grad(1:nWeights_grad))
        do i = 1, nWeights_grad
            @:ALLOCATE(weights_z_grad(i)%sf(idwbuff(3)%beg:idwbuff(3)%end,1:1,1:1))
            @:ACC_SETUP_SFs(weights_z_grad(i))
        end do

        ! Allocating space for lagrangian variables
        nParticles_glb = lag_params%nParticles_glb

        @:ALLOCATE(lag_part_id(1:nParticles_glb, 1:2))
        @:ALLOCATE(particle_R0(1:nParticles_glb))
        @:ALLOCATE(Rmax_stats_part(1:nParticles_glb))
        @:ALLOCATE(Rmin_stats_part(1:nParticles_glb))
        @:ALLOCATE(particle_mass(1:nParticles_glb))
        @:ALLOCATE(particle_rad(1:nParticles_glb, 1:2))
        @:ALLOCATE(particle_pos(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_posPrev(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_vel(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_s(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_draddt(1:nParticles_glb, 1:lag_num_ts))
        @:ALLOCATE(particle_dposdt(1:nParticles_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(particle_dveldt(1:nParticles_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(f_p(1:nParticles_glb, 1:3))
        @:ALLOCATE(linked_list(1:nParticles_glb))

        @:ALLOCATE(keep_bubble(1:nParticles_glb))
        @:ALLOCATE(wrap_bubble_loc(1:nParticles_glb, 1:num_dims), wrap_bubble_dir(1:nParticles_glb, 1:num_dims))

        if (adap_dt .and. f_is_default(adap_dt_tol)) adap_dt_tol = dflt_adap_dt_tol

        if (num_procs > 1) call s_initialize_solid_particles_mpi(lag_num_ts)

        ! Starting particles
        if (lag_params%write_void_evol) call s_open_void_evol
        if (lag_params%write_bubbles) call s_open_lag_bubble_evol()
        if (lag_params%write_bubbles_stats) call s_open_lag_particle_stats()

        if (lag_params%vel_model > 0) then
            moving_lag_particles = .true.
            lag_pressure_force = lag_params%pressure_force
            lag_gravity_force = lag_params%gravity_force
            lag_vel_model = lag_params%vel_model
            lag_drag_model = lag_params%drag_model
        end if

        $:GPU_UPDATE(device='[moving_lag_particles, lag_pressure_force, &
            & lag_gravity_force, lag_vel_model, lag_drag_model]')

        call s_read_input_particles(q_cons_vf)

        npts = (nWeights_interp - 1)/2
        call s_compute_barycentric_weights(npts) !For interpolation

        npts = (nWeights_grad - 1)/2
        call s_compute_fornberg_fd_weights(npts) !For finite differences

        call s_enforce_EL_particles_boundary_conditions(q_cons_vf, 0)

        if (lag_params%solver_approach == 2) then
            if (save_count == 0) then

                !> Correcting initial conditions so they account for particles
                $:GPU_PARALLEL_LOOP(private='[i,j,k,dyn_pres,fluid_vel,rho_f,alpharho,rhou,alpharhou]', collapse=3, copyin = '[pi_inf, qv, gamma, rhoYks]')
                do k = idwint(3)%beg, idwint(3)%end
                    do j = idwint(2)%beg, idwint(2)%end
                        do i = idwint(1)%beg, idwint(1)%end
                        !!!!!!!!! Mass
                            do l = 1, num_fluids !num_fluid is just 1 right now
                                rho_f = q_cons_vf(l)%sf(i, j, k)
                                alpharho = q_particles(alphaf_id)%sf(i, j, k)*rho_f
                                q_cons_vf(l)%sf(i, j, k) = alpharho
                            end do

                        !!!!!!!!! Momentum
                            dyn_pres = 0._wp
                            do l = momxb, momxe
                                fluid_vel(l - momxb + 1) = q_cons_vf(l)%sf(i, j, k)/rho_f
                                rhou = q_cons_vf(l)%sf(i, j, k)
                                alpharhou = q_particles(alphaf_id)%sf(i, j, k)*rhou
                                q_cons_vf(l)%sf(i, j, k) = alpharhou
                                dyn_pres = dyn_pres + q_cons_vf(l)%sf(i, j, k)* &
                                           fluid_vel(l - momxb + 1)/2._wp
                            end do

                        !!!!!!!!!Energy
                            call s_compute_pressure(q_cons_vf(E_idx)%sf(i, j, k), &
                                                    q_cons_vf(alf_idx)%sf(i, j, k), &
                                                    dyn_pres, pi_inf, gamma, alpharho, &
                                                    qv, rhoYks, pres, T)

                            q_cons_vf(E_idx)%sf(i, j, k) = &
                                gamma*pres + dyn_pres + pi_inf + qv !Updating energy in cons
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

            end if
        end if

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
                        lag_part_id(particle_id, 1) = id      !global ID
                        lag_part_id(particle_id, 2) = particle_id  !local ID
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

        $:GPU_UPDATE(device='[particles_lagrange, lag_params]')

        $:GPU_UPDATE(device='[lag_part_id,particle_R0,Rmax_stats_part,Rmin_stats_part,particle_mass,f_p, &
            & particle_rad,particle_pos,particle_posPrev,particle_vel, &
            & particle_s,particle_draddt, &
            & particle_dposdt,particle_dveldt,n_el_particles_loc]')

        Rmax_glb = min(dflt_real, -dflt_real)
        Rmin_glb = max(dflt_real, -dflt_real)
        $:GPU_UPDATE(device='[Rmax_glb, Rmin_glb]')

        $:GPU_UPDATE(device='[dx,dy,dz,x_cb,x_cc,y_cb,y_cc,z_cb,z_cc]')

        !Populate temporal variables
        call s_transfer_data_to_tmp_particles()

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
        Rmax_stats_part(part_id) = min(dflt_real, -dflt_real)
        Rmin_stats_part(part_id) = max(dflt_real, -dflt_real)
        particle_rad(part_id, 1) = inputPart(7)
        particle_pos(part_id, 1:3, 1) = inputPart(1:3)
        particle_posPrev(part_id, 1:3, 1) = particle_pos(part_id, 1:3, 1)
        particle_vel(part_id, 1:3, 1) = inputPart(4:6)

        !Initialize Particle Sources
        f_p(part_id, 1:3) = 0._wp

        if (cyl_coord .and. p == 0) then
            particle_pos(part_id, 2, 1) = sqrt(particle_pos(part_id, 2, 1)**2._wp + &
                                               particle_pos(part_id, 3, 1)**2._wp)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            particle_pos(part_id, 3, 1) = atan2(inputPart(3), inputPart(2))
            particle_posPrev(part_id, 1:3, 1) = particle_pos(part_id, 1:3, 1)
        end if

        cell = fd_number - buff_size
        call s_locate_cell(particle_pos(part_id, 1:3, 1), cell, particle_s(part_id, 1:3, 1))

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

        ! Initial particle mass
        volparticle = 4._wp/3._wp*pi*particle_R0(part_id)**3 ! volume
        ! gas_mv(particle_id, 1) = pv*volparticle*(1._wp/(R_v*Tw))*(massflag) ! vapermass
        particle_mass(part_id) = volparticle*rho0ref_particle !(gas_p(particle_id, 1) - pv*(massflag))*volparticle*(1._wp/(R_n*Tw)) ! gasmass
        if (particle_mass(part_id) <= 0._wp) then
            call s_mpi_abort("The initial particle mass is negative. Check the initial conditions.")
        end if
        ! totalmass = particle_mass(part_id) !+ gas_mv(particle_id, 1) ! totalmass

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
                lag_part_id(i, 1) = int(MPI_IO_DATA_lag_bubbles(i, 1))
                particle_pos(i, 1:3, 1) = MPI_IO_DATA_lag_bubbles(i, 2:4)
                particle_posPrev(i, 1:3, 1) = MPI_IO_DATA_lag_bubbles(i, 5:7)
                particle_vel(i, 1:3, 1) = MPI_IO_DATA_lag_bubbles(i, 8:10)
                particle_rad(i, 1) = MPI_IO_DATA_lag_bubbles(i, 11)
                ! intfc_vel(i, 1) = MPI_IO_DATA_lag_bubbles(i, 12)
                particle_R0(i) = MPI_IO_DATA_lag_bubbles(i, 13)
                Rmax_stats_part(i) = MPI_IO_DATA_lag_bubbles(i, 14)
                Rmin_stats_part(i) = MPI_IO_DATA_lag_bubbles(i, 15)
                ! bub_dphidt(i) = MPI_IO_DATA_lag_bubbles(i, 16)
                ! gas_p(i, 1) = MPI_IO_DATA_lag_bubbles(i, 17)
                ! gas_mv(i, 1) = MPI_IO_DATA_lag_bubbles(i, 18)
                particle_mass(i) = MPI_IO_DATA_lag_bubbles(i, 19)
                ! gas_betaT(i) = MPI_IO_DATA_lag_bubbles(i, 20)
                ! gas_betaC(i) = MPI_IO_DATA_lag_bubbles(i, 21)
                cell = -buff_size
                call s_locate_cell(particle_pos(i, 1:3, 1), cell, particle_s(i, 1:3, 1))
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
    subroutine s_compute_particle_EL_dynamics(q_prim_vf, stage, q_cons_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z)
#ifdef MFC_OpenMP
        !DIR$ OPTIMIZE (-O1)
#endif
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in) :: stage
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vL_x, vL_y, vL_z
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vR_x, vR_y, vR_z

        integer, dimension(3) :: cell
        real(wp) :: myMass, myR, myBeta_c, myBeta_t, myR0, myRe, myVolumeFrac, myGamma, rmass_add
        real(wp), dimension(3) :: myVel, myPos, force_vec

        integer :: k, l, i, j

        $:GPU_PARALLEL_LOOP(private='[k]')
        do k = 1, n_el_particles_loc
            f_p(k, :) = 0._wp
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (lag_params%pressure_force .or. lag_params%added_mass_model > 0) then
            do l = 1, num_dims
                if (l == 1) then
                    call s_gradient_field(vL_x, vR_x, field_vars(dPx_id)%sf, l, E_idx)
                elseif (l == 2) then
                    call s_gradient_field(vL_y, vR_y, field_vars(dPy_id)%sf, l, E_idx)
                elseif (l == 3) then
                    call s_gradient_field(vL_z, vR_z, field_vars(dPz_id)%sf, l, E_idx)
                end if
            end do
        end if

        if (lag_params%added_mass_model > 0) then

            $:GPU_PARALLEL_LOOP(private='[i,j,k]')
            do i = idwbuff(1)%beg, idwbuff(1)%end
                do j = idwbuff(2)%beg, idwbuff(2)%end
                    do k = idwbuff(3)%beg, idwbuff(3)%end

                        myL_rsx_vf(i, j, k, r_rhouf_id) = vL_x(i, j, k, 1)*vL_x(i, j, k, momxb)
                        myL_rsy_vf(j, i, k, r_rhouf_id) = vL_y(j, i, k, 1)*vL_y(j, i, k, momxb + 1)
                        myR_rsx_vf(i, j, k, r_rhouf_id) = vR_x(i, j, k, 1)*vR_x(i, j, k, momxb)
                        myR_rsy_vf(j, i, k, r_rhouf_id) = vR_y(j, i, k, 1)*vR_y(j, i, k, momxb + 1)

                        if (num_dims == 3) then
                            myL_rsz_vf(k, j, i, r_rhouf_id) = vL_z(k, j, i, 1)*vL_z(k, j, i, momxb + 2)
                            myR_rsz_vf(k, j, i, r_rhouf_id) = vR_z(k, j, i, 1)*vR_z(k, j, i, momxb + 2)
                        else
                            myL_rsz_vf(i, j, k, r_rhouf_id) = vL_z(i, j, k, 1)*vL_z(i, j, k, momxb + 2)
                            myR_rsz_vf(i, j, k, r_rhouf_id) = vR_z(i, j, k, 1)*vR_z(i, j, k, momxb + 2)
                        end if

                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            do l = 1, num_dims
                if (l == 1) then
                    call s_gradient_field(vL_x, vR_x, field_vars(dRhox_id)%sf, l, 1)
                    call s_gradient_field(myL_rsx_vf, myR_rsx_vf, field_vars(dRhoux_id)%sf, l, r_rhouf_id)
                elseif (l == 2) then
                    call s_gradient_field(vL_y, vR_y, field_vars(dRhoy_id)%sf, l, 1)
                    call s_gradient_field(myL_rsy_vf, myR_rsy_vf, field_vars(dRhouy_id)%sf, l, r_rhouf_id)
                elseif (l == 3) then
                    call s_gradient_field(vL_z, vR_z, field_vars(dRhoz_id)%sf, l, 1)
                    call s_gradient_field(myL_rsz_vf, myR_rsz_vf, field_vars(dRhouz_id)%sf, l, r_rhouf_id)
                end if
            end do

        end if

        myGamma = (1._wp/fluid_pp(1)%gamma) + 1._wp
        myRe = 1.48e-5 !fluid_pp(1)%Re(1) !Need a viscosity model for when modeling inviscid eulerian fluid

        call nvtxStartRange("LAGRANGE-PARTICLE-DYNAMICS")

        if (lag_params%collision_force) then
            !First compute collision forces
            call s_compute_particle_EL_collisions(stage)
        end if

        !Then compute drag/pressure/added mass forces and convert to particle acceleration
        $:GPU_PARALLEL_LOOP(private='[k,l,cell,myMass,myR,myR0,myPos,myVel,myVolumeFrac,force_vec,rmass_add]',&
        & copyin='[stage, myGamma, myRe]')
        do k = 1, n_el_particles_loc

            cell = fd_number - buff_size
            call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))

            ! Current particle state
            myMass = particle_mass(k)
            myR = particle_rad(k, 2)
            myR0 = particle_R0(k)
            myPos = particle_pos(k, :, 2)
            myVel = particle_vel(k, :, 2)
            myVolumeFrac = 1._wp - q_particles(alphaf_id)%sf(cell(1), cell(2), cell(3))

            particle_dposdt(k, :, stage) = 0._wp
            particle_dveldt(k, :, stage) = 0._wp
            particle_draddt(k, stage) = 0._wp

            call s_get_particle_force(myPos, myR, myVel, myMass, myRe, myGamma, myVolumeFrac, cell, &
                                      q_prim_vf, field_vars, weights_x_interp, weights_y_interp, weights_z_interp, &
                                      force_vec, rmass_add)

            f_p(k, :) = f_p(k, :) + force_vec(:)

            do l = 1, num_dims
                particle_dposdt(k, l, stage) = myVel(l)
                particle_dveldt(k, l, stage) = f_p(k, l)/(myMass + rmass_add)
                particle_draddt(k, stage) = 0._wp
            end do

        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

        call nvtxStartRange("DYNAMICS-DEV2HOST")
        $:GPU_UPDATE(host='[f_p, particle_dposdt, particle_dveldt, particle_draddt]')
        call nvtxEndRange

        ! call s_eliminate_ghost_particles()
        ! call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage)

    end subroutine s_compute_particle_EL_dynamics

    !>  Contains the particle collision force computation.
    subroutine s_compute_particle_EL_collisions(stage)

        integer, intent(in) :: stage
        integer, dimension(3) :: cell
        integer, dimension(3) :: cellaux
        integer :: k, l, q, ip, jp, kp, ii, jj, kk
        logical :: celloutside
        real(wp) :: pidtksp2, ksp, nu1, nu2, Rp1, Rp2, E1, E2, Estar, cor, rmag, Rstar, dij, eta_n, kappa_n, mp1, mp2, dt_loc
        real(wp), dimension(3) :: xp1, xp2, vp1, vp2, v_rel, rpij, nij, vnij, Fnpp_ij
        integer :: kpz
        integer :: newBubs

        if (lag_num_ts == 1) then
            dt_loc = dt
        elseif (lag_num_ts == 2) then
            if (stage == 1) then
                dt_loc = dt
            elseif (stage == 2) then
                dt_loc = dt/2._wp
            end if
        elseif (lag_num_ts == 3) then
            if (stage == 1) then
                dt_loc = dt
            elseif (stage == 2) then
                dt_loc = dt/4._wp
            elseif (stage == 3) then
                dt_loc = (2._wp/3._wp)*dt
            end if
        end if

        kpz = 0
        if (num_dims == 3) kpz = 1

        ksp = 10._wp
        nu1 = 0.35
        nu2 = 0.35
        E1 = 1e9
        E2 = 1e9
        cor = 0.7

        pidtksp2 = (pi**2)/((dt_loc*ksp)**2)

        Estar = 1/(((1 - nu1**2)/E1) + ((1 - nu2**2)/E2))
        Estar = (4._wp/3._wp)*Estar

        call nvtxStartRange("LAGRANGE-PARTICLE-LINKED-LIST")

        $:GPU_PARALLEL_LOOP(private='[k,cell]')
        do k = 1, n_el_particles_loc

            cell = fd_number - buff_size
            call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))

            !!!!! Particle linked list building
            $:GPU_ATOMIC(atomic='capture')
            linked_list(k) = particle_head(1)%sf(cell(1), cell(2), cell(3))
            particle_head(1)%sf(cell(1), cell(2), cell(3)) = k
            $:END_GPU_ATOMIC()

        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

        call nvtxStartRange("LAGRANGE-PARTICLE-COLLISIONS")

        $:GPU_PARALLEL_LOOP(private='[k,cell,ip,jp,kp,Rp1,xp1,mp1,vp1,kk,jj,ii,cellaux,q,Rp2,xp2,mp2,vp2,v_rel,Rstar,rpij,rmag,nij,vnij,dij,kappa_n,eta_n,Fnpp_ij]',&
        & copyin='[ksp,nu1,nu2,E1,E2,cor,pidtksp2,Estar,kpz,dt_loc]')
        do k = 1, n_el_particles_loc

            cell = fd_number - buff_size
            call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))

            ip = cell(1)
            jp = cell(2)
            kp = cell(3)

            Rp1 = particle_rad(k, 2)
            xp1 = particle_pos(k, :, 2)
            mp1 = particle_mass(k)
            vp1 = particle_vel(k, :, 2)

            do kk = kp - kpz, kp + kpz
                do jj = jp - 1, jp + 1
                    do ii = ip - 1, ip + 1

                        cellaux(1) = ii
                        cellaux(2) = jj
                        cellaux(3) = kk

                        call s_check_celloutside_wbuff(cellaux, celloutside)

                        if (.not. celloutside) then
                            q = particle_head(1)%sf(ii, jj, kk)
                            ! Traverse linked list in that cell
                            do while (q /= -1)

                                if (q > k) then

                                    Rp2 = particle_rad(q, 2)
                                    xp2 = particle_pos(q, :, 2)
                                    mp2 = particle_mass(q)
                                    vp2 = particle_vel(q, :, 2)
                                    v_rel = vp2 - vp1

                                    Rstar = (Rp1*Rp2)/(Rp1 + Rp2)
                                    rpij = xp2 - xp1
                                    rmag = sqrt(rpij(1)**2 + rpij(2)**2 + rpij(3)**2)
                                    nij = rpij/rmag
                                    vnij = dot_product(v_rel, nij)*nij
                                    dij = (Rp1 + Rp2) - rmag

                                    if (dij > 0._wp) then

                                        kappa_n = min((pidtksp2*mp1), (Estar*sqrt(Rstar)*sqrt(abs(dij))))

                                        eta_n = ((-2._wp*sqrt(kappa_n)*log(cor))/sqrt((log(cor))**2 + pi**2))*(1._wp/sqrt((1/mp1) + (1/mp2)))

                                        Fnpp_ij = -kappa_n*dij*nij - eta_n*vnij

                                        f_p(k, :) = f_p(k, :) + Fnpp_ij
                                        f_p(q, :) = f_p(q, :) - Fnpp_ij

                                        ! write(*,*) "pidtksp2 * mp1 = ", pidtksp2 * mp1
                                        ! write(*,*) "(Estar * sqrt(Rstar) * sqrt(abs(dij))) = ", (Estar * sqrt(Rstar) * sqrt(abs(dij)))
                                        ! write(*,*) "Colliding particles"
                                        ! write(*,*) "Fnpp_ij = ", Fnpp_ij
                                        ! write(*,*) "k = ", k
                                        ! write(*,*) "f_p(k,:) = ", f_p(k,:)
                                        ! write(*,*) "q = ", q
                                        ! write(*,*) "f_p(q,:) = ", f_p(q,:)
                                    end if

                                end if

                                q = linked_list(q)

                            end do
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_compute_particle_EL_collisions

    !>  Contains deletion of ghost particles
    subroutine s_eliminate_ghost_particles()

        integer :: newBubs, k
        !Remove ghost particles immediately after done using them

        $:GPU_PARALLEL_LOOP(private='[k]')
        do k = 1, n_el_particles_loc
            keep_bubble(k) = 1

            if (particle_pos(k, 1, 2) >= x_cb(m)) then
                keep_bubble(k) = 0
            elseif (particle_pos(k, 1, 2) < x_cb(-1)) then
                keep_bubble(k) = 0
            end if

            if (particle_pos(k, 2, 2) >= y_cb(n)) then
                keep_bubble(k) = 0
            elseif (particle_pos(k, 2, 2) < y_cb(-1)) then
                keep_bubble(k) = 0
            end if

            if (p > 0) then
                if (particle_pos(k, 3, 2) >= z_cb(p)) then
                    keep_bubble(k) = 0
                elseif (particle_pos(k, 3, 2) < z_cb(-1)) then
                    keep_bubble(k) = 0
                end if
            end if

        end do
        $:END_GPU_PARALLEL_LOOP()

        if (n_el_particles_loc > 0) then
            call nvtxStartRange("LAG-BC")
            call nvtxStartRange("LAG-BC-DEV2HOST")
            $:GPU_UPDATE(host='[particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
                & lag_part_id, particle_rad, &
                & particle_pos, particle_posPrev, particle_vel, particle_s, particle_draddt, &
                & particle_dposdt, particle_dveldt, keep_bubble, n_el_particles_loc, &
                & wrap_bubble_dir, wrap_bubble_loc]')
            call nvtxEndRange

            newBubs = 0
            do k = 1, n_el_particles_loc
                if (keep_bubble(k) == 1) then
                    newBubs = newBubs + 1
                    if (newBubs /= k) then
                        call s_copy_lag_particle(newBubs, k)
                    end if
                end if
            end do

            n_el_particles_loc = newBubs

            call nvtxStartRange("LAG-BC-HOST2DEV")
            $:GPU_UPDATE(device='[particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
            & lag_part_id, particle_rad, &
            & particle_pos, particle_posPrev, particle_vel, particle_s, particle_draddt, &
            & particle_dposdt, particle_dveldt, n_el_particles_loc]')
            call nvtxEndRange
        end if

    end subroutine s_eliminate_ghost_particles

    !> The purpose of this subroutine is to check if the current cell is outside the computational domain or not (including ghost cells).
            !! @param cellaux Tested cell to smear the particle effect in.
            !! @param celloutside If true, then cellaux is outside the computational domain.
    subroutine s_check_celloutside_wbuff(cellaux, celloutside)
        $:GPU_ROUTINE(function_name='s_check_celloutside_wbuff',parallelism='[seq]', &
            & cray_inline=True)

        integer, dimension(3), intent(inout) :: cellaux
        logical, intent(out) :: celloutside

        celloutside = .false.

        if (num_dims == 2) then
            if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                celloutside = .true.
            end if

            if ((cellaux(1) > m + buff_size) .or. (cellaux(2) > n + buff_size)) then
                celloutside = .true.
            end if
        else
            if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size) .or. (cellaux(3) < -buff_size)) then
                celloutside = .true.
            end if

            if ((cellaux(1) > m + buff_size) .or. (cellaux(2) > n + buff_size) .or. (cellaux(3) > p + buff_size)) then
                celloutside = .true.
            end if
        end if

    end subroutine s_check_celloutside_wbuff

    !>  The purpose of this subroutine is to obtain the bubble source terms based on Maeda and Colonius (2018)
        !!      and add them to the RHS scalar field.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
        !! @param rhs_vf Time derivative of the conservative variables
    subroutine s_compute_particles_EL_source(q_cons_vf, q_prim_vf, rhs_vf, stage)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        integer :: i, j, k, l, nf, stage
        real(wp) :: dalphapdt, alpha_f, udot_gradalpha

        call s_smear_particle_sources() !To fill in q_particles with volume fraction and source term contributions

        do l = 1, num_dims
            iv_loc%beg = r_alphaf_id; iv_loc%end = r_alphaupz_id
            call s_reconstruct_cell_boundary_values_loc( &
                q_particles(iv_loc%beg:iv_loc%end), &
                myL_rsx_vf, myL_rsy_vf, myL_rsz_vf, &
                myR_rsx_vf, myR_rsy_vf, myR_rsz_vf, &
                l)

        end do

        do l = 1, num_dims
            if (l == 1) then
                call s_gradient_field(myL_rsx_vf, myR_rsx_vf, field_vars(dalphafx_id)%sf, l, r_alphaf_id)
                call s_gradient_field(myL_rsx_vf, myR_rsx_vf, field_vars(dalphap_upx_id)%sf, l, r_alphaupx_id)
            elseif (l == 2) then
                call s_gradient_field(myL_rsy_vf, myR_rsy_vf, field_vars(dalphafy_id)%sf, l, r_alphaf_id)
                call s_gradient_field(myL_rsy_vf, myR_rsy_vf, field_vars(dalphap_upy_id)%sf, l, r_alphaupy_id)
            elseif (l == 3) then
                call s_gradient_field(myL_rsz_vf, myR_rsz_vf, field_vars(dalphafz_id)%sf, l, r_alphaf_id)
                call s_gradient_field(myL_rsz_vf, myR_rsz_vf, field_vars(dalphap_upz_id)%sf, l, r_alphaupz_id)
            end if
        end do

        !> Apply particle sources to the Eulerian RHS
        $:GPU_PARALLEL_LOOP(private='[i,j,k,alpha_f,dalphapdt,udot_gradalpha]', collapse=3)
        do k = idwint(3)%beg, idwint(3)%end
            do j = idwint(2)%beg, idwint(2)%end
                do i = idwint(1)%beg, idwint(1)%end
                    if (q_particles(alphaf_id)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then

                        alpha_f = q_particles(alphaf_id)%sf(i, j, k)

                        dalphapdt = 0
                        do l = 1, num_dims
                            dalphapdt = dalphapdt + field_vars(dalphap_upx_id + l - 1)%sf(i, j, k)
                        end do
                        dalphapdt = -dalphapdt
                        !Add any contribution to dalphapdt from particles growing or shrinking

                        udot_gradalpha = 0._wp
                        do l = 1, num_dims
                            udot_gradalpha = udot_gradalpha + q_prim_vf(momxb + l - 1)%sf(i, j, k)*field_vars(dalphafx_id + l - 1)%sf(i, j, k)
                        end do

                        !> Step 1: Source terms for volume fraction corrections
                        !cons_var/alpha_f * (dalpha_p/dt - u dot grad(alpha_f))
                        do l = 1, E_idx
                            rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + &
                                                    (q_cons_vf(l)%sf(i, j, k)/alpha_f)*(dalphapdt - udot_gradalpha)
                        end do

                        !momentum term -1/alpha_f * (p*grad(alpha_f) - Tau^v dot grad(alpha_f)) !Viscous term not implemented
                        do l = 1, num_dims
                            rhs_vf(momxb + l - 1)%sf(i, j, k) = rhs_vf(momxb + l - 1)%sf(i, j, k) - &
                                                                ((1._wp/alpha_f)* &
                                                                 (q_prim_vf(E_idx)%sf(i, j, k)*field_vars(dalphafx_id + l - 1)%sf(i, j, k)))
                        end do

                        !energy term -1/alpha_f * (p*u dot grad(alpha_f) - (Tau^v dot u) dot grad(alpha_f)) !Viscous term not implemented
                        rhs_vf(E_idx)%sf(i, j, k) = rhs_vf(E_idx)%sf(i, j, k) - &
                                                    ((1._wp/alpha_f)* &
                                                     (q_prim_vf(E_idx)%sf(i, j, k)*udot_gradalpha))

                        !> Step 2: Add the drag/pressure/added mass forces to the fluid

                        rhs_vf(momxb)%sf(i, j, k) = rhs_vf(momxb)%sf(i, j, k) + q_particles(Smx_id)%sf(i, j, k)*(1._wp/alpha_f)
                        rhs_vf(momxb + 1)%sf(i, j, k) = rhs_vf(momxb + 1)%sf(i, j, k) + q_particles(Smy_id)%sf(i, j, k)*(1._wp/alpha_f)

                        ! Energy source
                        rhs_vf(E_idx)%sf(i, j, k) = rhs_vf(E_idx)%sf(i, j, k) + &
                                                    (q_particles(Smx_id)%sf(i, j, k)*q_prim_vf(momxb)%sf(i, j, k) &
                                                     + q_particles(Smy_id)%sf(i, j, k)*q_prim_vf(momxb + 1)%sf(i, j, k) &
                                                     + q_particles(SE_id)%sf(i, j, k))*(1._wp/alpha_f)

                        if (num_dims == 3) then
                            rhs_vf(momxb + 2)%sf(i, j, k) = rhs_vf(momxb + 2)%sf(i, j, k) + q_particles(Smz_id)%sf(i, j, k)*(1._wp/alpha_f)
                            ! Energy source
                            rhs_vf(E_idx)%sf(i, j, k) = rhs_vf(E_idx)%sf(i, j, k) + &
                                                        (q_particles(Smz_id)%sf(i, j, k)*q_prim_vf(momxb + 2)%sf(i, j, k))*(1._wp/alpha_f)
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

        $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    particle_head(1)%sf(j, k, l) = -1
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
        do i = 1, nField_vars
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        field_vars(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

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
                              particle_s, particle_pos, particle_vel, q_particles, f_p, ncc)

        !Store 1-q_particles(1)
        $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    q_particles(alphaf_id)%sf(j, k, l) = 1._wp - q_particles(alphaf_id)%sf(j, k, l)
                    ! Limiting void fraction given max value
                    q_particles(alphaf_id)%sf(j, k, l) = max(q_particles(alphaf_id)%sf(j, k, l), &
                                                             1._wp - lag_params%valmaxvoid)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call nvtxEndRange

    end subroutine s_smear_particle_sources

    !>  This subroutine updates the Lagrange variables using the tvd RK time steppers.
        !!      The time derivative of the particle variables must be stored at every stage to avoid precision errors.
        !! @param stage Current tvd RK stage
    impure subroutine s_update_lagrange_particles_tdv_rk(q_prim_vf, stage)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: stage

        integer :: k

        if (time_stepper == 1) then ! 1st order TVD RK

            $:GPU_PARALLEL_LOOP(private='[k]')
            do k = 1, n_el_particles_loc
                !u{1} = u{n} +  dt * RHS{n}
                particle_rad(k, 1) = particle_rad(k, 1) + dt*particle_draddt(k, 1)
                if (moving_lag_particles) then
                    particle_posPrev(k, 1:3, 1) = particle_pos(k, 1:3, 1)
                    particle_pos(k, 1:3, 1) = particle_pos(k, 1:3, 1) + dt*particle_dposdt(k, 1:3, 1)
                    particle_vel(k, 1:3, 1) = particle_vel(k, 1:3, 1) + dt*particle_dveldt(k, 1:3, 1)
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()

            call s_transfer_data_to_tmp_particles()
            if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage)
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
                    particle_rad(k, 2) = particle_rad(k, 1) + dt*particle_draddt(k, 1)
                    if (moving_lag_particles) then
                        particle_posPrev(k, 1:3, 2) = particle_pos(k, 1:3, 1)
                        particle_pos(k, 1:3, 2) = particle_pos(k, 1:3, 1) + dt*particle_dposdt(k, 1:3, 1)
                        particle_vel(k, 1:3, 2) = particle_vel(k, 1:3, 1) + dt*particle_dveldt(k, 1:3, 1)
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage)

            elseif (stage == 2) then

                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    !u{1} = u{n} + (1/2) * dt * (RHS{n} + RHS{1})
                    particle_rad(k, 1) = particle_rad(k, 1) + dt*(particle_draddt(k, 1) + particle_draddt(k, 2))/2._wp
                    if (moving_lag_particles) then
                        particle_posPrev(k, 1:3, 1) = particle_pos(k, 1:3, 2)
                        particle_pos(k, 1:3, 1) = particle_pos(k, 1:3, 1) + dt*(particle_dposdt(k, 1:3, 1) + particle_dposdt(k, 1:3, 2))/2._wp
                        particle_vel(k, 1:3, 1) = particle_vel(k, 1:3, 1) + dt*(particle_dveldt(k, 1:3, 1) + particle_dveldt(k, 1:3, 2))/2._wp
                    end if

                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_transfer_data_to_tmp_particles()
                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage)
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
                    particle_rad(k, 2) = particle_rad(k, 1) + dt*particle_draddt(k, 1)
                    if (moving_lag_particles) then
                        particle_posPrev(k, 1:3, 2) = particle_pos(k, 1:3, 1)
                        particle_pos(k, 1:3, 2) = particle_pos(k, 1:3, 1) + dt*particle_dposdt(k, 1:3, 1)
                        particle_vel(k, 1:3, 2) = particle_vel(k, 1:3, 1) + dt*particle_dveldt(k, 1:3, 1)
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage)

            elseif (stage == 2) then

                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    !u{2} = u{n} + (1/4) * dt * [RHS{n} + RHS{1}]
                    particle_rad(k, 2) = particle_rad(k, 1) + dt*(particle_draddt(k, 1) + particle_draddt(k, 2))/4._wp
                    if (moving_lag_particles) then
                        particle_posPrev(k, 1:3, 2) = particle_pos(k, 1:3, 2)
                        particle_pos(k, 1:3, 2) = particle_pos(k, 1:3, 1) + dt*(particle_dposdt(k, 1:3, 1) + particle_dposdt(k, 1:3, 2))/4._wp
                        particle_vel(k, 1:3, 2) = particle_vel(k, 1:3, 1) + dt*(particle_dveldt(k, 1:3, 1) + particle_dveldt(k, 1:3, 2))/4._wp
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage)

            elseif (stage == 3) then

                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    !u{n+1} = u{n} + (2/3) * dt * [(1/4)* RHS{n} + (1/4)* RHS{1} + RHS{2}]
                    particle_rad(k, 1) = particle_rad(k, 1) + (2._wp/3._wp)*dt*(particle_draddt(k, 1)/4._wp + particle_draddt(k, 2)/4._wp + particle_draddt(k, 3))
                    if (moving_lag_particles) then
                        particle_posPrev(k, 1:3, 1) = particle_pos(k, 1:3, 2)
                        particle_pos(k, 1:3, 1) = particle_pos(k, 1:3, 1) + (2._wp/3._wp)*dt*(particle_dposdt(k, 1:3, 1)/4._wp + particle_dposdt(k, 1:3, 2)/4._wp + particle_dposdt(k, 1:3, 3))
                        particle_vel(k, 1:3, 1) = particle_vel(k, 1:3, 1) + (2._wp/3._wp)*dt*(particle_dveldt(k, 1:3, 1)/4._wp + particle_dveldt(k, 1:3, 2)/4._wp + particle_dveldt(k, 1:3, 3))
                    end if

                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_transfer_data_to_tmp_particles()
                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage)
                if (lag_params%write_void_evol) call s_write_void_evol_particles(mytime)
                if (lag_params%write_bubbles_stats) call s_calculate_lag_particle_stats()
                if (lag_params%write_bubbles .and. mytime >= next_write_time) then
                    $:GPU_UPDATE(host='[particle_mass,particle_rad]')
                    call s_write_lag_particle_evol(mytime)
                    next_write_time = next_write_time + t_save
                end if

            end if

        end if

    end subroutine s_update_lagrange_particles_tdv_rk

    !> This subroutine enforces reflective and wall boundary conditions for EL particles
        !! @param dest Destination for the bubble position update
    impure subroutine s_enforce_EL_particles_boundary_conditions(q_prim_vf, nstage)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer :: k, i, q
        integer :: patch_id, newBubs, new_idx
        real(wp) :: offset
        integer, dimension(3) :: cell
        integer, intent(in) :: nstage
        logical :: inc_ghost = .true.

        call nvtxStartRange("LAG-BC")
        call nvtxStartRange("LAG-BC-DEV2HOST")
        $:GPU_UPDATE(host='[particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
            & lag_part_id, particle_rad, &
            & particle_pos, particle_posPrev, particle_vel, particle_s, particle_draddt, &
            & particle_dposdt, particle_dveldt, keep_bubble, n_el_particles_loc, &
            & wrap_bubble_dir, wrap_bubble_loc]')
        call nvtxEndRange

        ! Handle MPI transfer of particles going to another processor's local domain
        if (num_procs > 1) then
            call nvtxStartRange("LAG-BC-TRANSFER-LIST")
            call s_add_particles_to_transfer_list(n_el_particles_loc, particle_pos(:, :, 2), particle_posPrev(:, :, 2), inc_ghost)
            call nvtxEndRange

            call nvtxStartRange("LAG-BC-SENDRECV")
            call s_mpi_sendrecv_solid_particles(particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
                                                lag_part_id, &
                                                particle_rad, particle_pos, particle_posPrev, particle_vel, &
                                                particle_s, particle_draddt, particle_dposdt, particle_dveldt, lag_num_ts, n_el_particles_loc, &
                                                2)
            call nvtxEndRange
        end if

        call nvtxStartRange("LAG-BC-HOST2DEV")
        $:GPU_UPDATE(device='[particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
            & lag_part_id, particle_rad, &
            & particle_pos, particle_posPrev, particle_vel, particle_s, particle_draddt, &
            & particle_dposdt, particle_dveldt, n_el_particles_loc]')
        call nvtxEndRange

        $:GPU_PARALLEL_LOOP(private='[k, cell]',copyin='[nstage]')
        do k = 1, n_el_particles_loc
            keep_bubble(k) = 1
            wrap_bubble_loc(k, :) = 0
            wrap_bubble_dir(k, :) = 0

            ! Relocate particles at solid boundaries and delete particles that leave
            ! buffer regions
            if (any(bc_x%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                .and. particle_pos(k, 1, 2) < x_cb(-1) + particle_rad(k, 2)) then
                particle_pos(k, 1, 2) = x_cb(-1) + particle_rad(k, 2)
                if (nstage == lag_num_ts) then
                    particle_pos(k, 1, 1) = particle_pos(k, 1, 2)
                end if
            elseif (any(bc_x%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                    .and. particle_pos(k, 1, 2) > x_cb(m) - particle_rad(k, 2)) then
                particle_pos(k, 1, 2) = x_cb(m) - particle_rad(k, 2)
                if (nstage == lag_num_ts) then
                    particle_pos(k, 1, 1) = particle_pos(k, 1, 2)
                end if
            elseif (bc_x%beg == BC_PERIODIC .and. particle_pos(k, 1, 2) < pcomm_coords(1)%beg .and. &
                    particle_posPrev(k, 1, 2) >= pcomm_coords(1)%beg) then
                wrap_bubble_dir(k, 1) = 1
                wrap_bubble_loc(k, 1) = -1
            elseif (bc_x%end == BC_PERIODIC .and. particle_pos(k, 1, 2) > pcomm_coords(1)%end .and. &
                    particle_posPrev(k, 1, 2) <= pcomm_coords(1)%end) then
                wrap_bubble_dir(k, 1) = 1
                wrap_bubble_loc(k, 1) = 1
            elseif (particle_pos(k, 1, 2) >= x_cb(m + ncc) .or. particle_pos(k, 1, 2) > glb_bounds(1)%end) then
                keep_bubble(k) = 0
            elseif (particle_pos(k, 1, 2) < x_cb(-1 - ncc) .or. particle_pos(k, 1, 2) < glb_bounds(1)%beg) then
                keep_bubble(k) = 0
            end if

            if (any(bc_y%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                .and. particle_pos(k, 2, 2) < y_cb(-1) + particle_rad(k, 2)) then
                particle_pos(k, 2, 2) = y_cb(-1) + particle_rad(k, 2)
                if (nstage == lag_num_ts) then
                    particle_pos(k, 2, 1) = particle_pos(k, 2, 2)
                end if
            else if (any(bc_y%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                     .and. particle_pos(k, 2, 2) > y_cb(n) - particle_rad(k, 2)) then
                particle_pos(k, 2, 2) = y_cb(n) - particle_rad(k, 2)
                if (nstage == lag_num_ts) then
                    particle_pos(k, 2, 1) = particle_pos(k, 2, 2)
                end if
            elseif (bc_y%beg == BC_PERIODIC .and. particle_pos(k, 2, 2) < pcomm_coords(2)%beg .and. &
                    particle_posPrev(k, 2, 2) >= pcomm_coords(2)%beg) then
                wrap_bubble_dir(k, 2) = 1
                wrap_bubble_loc(k, 2) = -1
            elseif (bc_y%end == BC_PERIODIC .and. particle_pos(k, 2, 2) > pcomm_coords(2)%end .and. &
                    particle_posPrev(k, 2, 2) <= pcomm_coords(2)%end) then
                wrap_bubble_dir(k, 2) = 1
                wrap_bubble_loc(k, 2) = 1
            elseif (particle_pos(k, 2, 2) >= y_cb(n + ncc) .or. particle_pos(k, 2, 2) > glb_bounds(2)%end) then
                keep_bubble(k) = 0
            elseif (particle_pos(k, 2, 2) < y_cb(-1 - ncc) .or. particle_pos(k, 2, 2) < glb_bounds(2)%beg) then
                keep_bubble(k) = 0
            end if

            if (p > 0) then
                if (any(bc_z%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                    .and. particle_pos(k, 3, 2) < z_cb(-1) + particle_rad(k, 2)) then
                    particle_pos(k, 3, 2) = z_cb(-1) + particle_rad(k, 2)
                    if (nstage == lag_num_ts) then
                        particle_pos(k, 3, 1) = particle_pos(k, 3, 2)
                    end if
                else if (any(bc_z%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                         .and. particle_pos(k, 3, 2) > z_cb(p) - particle_rad(k, 2)) then
                    particle_pos(k, 3, 2) = z_cb(p) - particle_rad(k, 2)
                    if (nstage == lag_num_ts) then
                        particle_pos(k, 3, 1) = particle_pos(k, 3, 2)
                    end if
                elseif (bc_z%beg == BC_PERIODIC .and. particle_pos(k, 3, 2) < pcomm_coords(3)%beg .and. &
                        particle_posPrev(k, 3, 2) >= pcomm_coords(3)%beg) then
                    wrap_bubble_dir(k, 3) = 1
                    wrap_bubble_loc(k, 3) = -1
                elseif (bc_z%end == BC_PERIODIC .and. particle_pos(k, 3, 2) > pcomm_coords(3)%end .and. &
                        particle_posPrev(k, 3, 2) <= pcomm_coords(3)%end) then
                    wrap_bubble_dir(k, 3) = 1
                    wrap_bubble_loc(k, 3) = 1
                elseif (particle_pos(k, 3, 2) >= z_cb(p + ncc) .or. particle_pos(k, 3, 2) > glb_bounds(3)%end) then
                    keep_bubble(k) = 0
                elseif (particle_pos(k, 3, 2) < z_cb(-1 - ncc) .or. particle_pos(k, 3, 2) < glb_bounds(3)%beg) then
                    keep_bubble(k) = 0
                end if
            end if

            if (keep_bubble(k) == 1) then
                ! Remove bubbles that are no longer in a liquid
                cell = fd_number - buff_size
                call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))

                ! if (q_prim_vf(advxb)%sf(cell(1), cell(2), cell(3)) < (1._wp - lag_params%valmaxvoid)) then
                !     keep_bubble(k) = 0
                ! end if

                ! Move bubbles back to surface of IB
                if (ib) then
                    cell = fd_number - buff_size
                    call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))

                    if (ib_markers%sf(cell(1), cell(2), cell(3)) /= 0) then
                        patch_id = ib_markers%sf(cell(1), cell(2), cell(3))

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_dims
                            particle_pos(k, i, 2) = particle_pos(k, i, 2) - &
                                                    levelset_norm%sf(cell(1), cell(2), cell(3), patch_id, i) &
                                                    *levelset%sf(cell(1), cell(2), cell(3), patch_id)
                        end do

                        cell = fd_number - buff_size
                        call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))
                    end if
                end if
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (n_el_particles_loc > 0) then
            call nvtxStartRange("LAG-BC")
            call nvtxStartRange("LAG-BC-DEV2HOST")
            $:GPU_UPDATE(host='[particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
                & lag_part_id, particle_rad, &
                & particle_pos, particle_posPrev, particle_vel, particle_s, particle_draddt, &
                & particle_dposdt, particle_dveldt, keep_bubble, n_el_particles_loc, &
                & wrap_bubble_dir, wrap_bubble_loc]')
            call nvtxEndRange

            newBubs = 0
            do k = 1, n_el_particles_loc
                if (keep_bubble(k) == 1) then
                    newBubs = newBubs + 1
                    if (newBubs /= k) then
                        call s_copy_lag_particle(newBubs, k)
                        wrap_bubble_dir(newBubs, :) = wrap_bubble_dir(k, :)
                        wrap_bubble_loc(newBubs, :) = wrap_bubble_loc(k, :)
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
                                    particle_pos(new_idx, i, q) = particle_pos(new_idx, i, q) - offset
                                    particle_posPrev(new_idx, i, q) = particle_posPrev(new_idx, i, q) - offset
                                end do
                            else if (wrap_bubble_loc(k, i) == -1) then
                                do q = 1, 2
                                    particle_pos(new_idx, i, q) = particle_pos(new_idx, i, q) + offset
                                    particle_posPrev(new_idx, i, q) = particle_posPrev(new_idx, i, q) + offset
                                end do
                            end if
                        end if
                    end do
                end if
            end do
            call nvtxStartRange("LAG-BC-HOST2DEV")
            $:GPU_UPDATE(device='[particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
            & lag_part_id, particle_rad, &
            & particle_pos, particle_posPrev, particle_vel, particle_s, particle_draddt, &
            & particle_dposdt, particle_dveldt, n_el_particles_loc]')
            call nvtxEndRange
        end if

        $:GPU_PARALLEL_LOOP(private='[cell]')
        do k = 1, n_el_particles_loc
            cell = fd_number - buff_size
            call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))
        end do

        ! Update void fraction and communicate buffers
        call s_smear_particle_sources()

        call nvtxEndRange ! LAG-BC

    end subroutine s_enforce_EL_particles_boundary_conditions

    !>  The purpose of this subroutine is to WENO-reconstruct the
        !!      left and the right cell-boundary values, including values
        !!      at the Gaussian quadrature points, from the cell-averaged
        !!      variables.
        !!  @param v_vf Cell-average variables
        !!  @param vL_qp Left WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param vR_qp Right WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_reconstruct_cell_boundary_values_loc(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, &
                                                      norm_dir)

        type(scalar_field), dimension(iv_loc%beg:iv_loc%end), intent(in) :: v_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vL_x, vL_y, vL_z
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vR_x, vR_y, vR_z
        integer, intent(in) :: norm_dir

        integer :: recon_dir !< Coordinate direction of the reconstruction

        integer :: i, j, k, l

        #:for SCHEME, TYPE in [('weno','WENO_TYPE'), ('muscl','MUSCL_TYPE')]
            if (recon_type == ${TYPE}$ .or. dummy) then
                ! Reconstruction in s1-direction
                if (norm_dir == 1) then
                    is1_loc = idwbuff(1); is2_loc = idwbuff(2); is3_loc = idwbuff(3)
                    recon_dir = 1; is1_loc%beg = is1_loc%beg + ${SCHEME}$_polyn
                    is1_loc%end = is1_loc%end - ${SCHEME}$_polyn

                elseif (norm_dir == 2) then
                    is1_loc = idwbuff(2); is2_loc = idwbuff(1); is3_loc = idwbuff(3)
                    recon_dir = 2; is1_loc%beg = is1_loc%beg + ${SCHEME}$_polyn
                    is1_loc%end = is1_loc%end - ${SCHEME}$_polyn

                else
                    is1_loc = idwbuff(3); is2_loc = idwbuff(2); is3_loc = idwbuff(1)
                    recon_dir = 3; is1_loc%beg = is1_loc%beg + ${SCHEME}$_polyn
                    is1_loc%end = is1_loc%end - ${SCHEME}$_polyn
                end if

                if (n > 0) then
                    if (p > 0) then
                        call s_${SCHEME}$ (v_vf(iv_loc%beg:iv_loc%end), &
                                           vL_x(:, :, :, iv_loc%beg:iv_loc%end), vL_y(:, :, :, iv_loc%beg:iv_loc%end), vL_z(:, :, :, iv_loc%beg:iv_loc%end), vR_x(:, :, :, iv_loc%beg:iv_loc%end), vR_y(:, :, :, iv_loc%beg:iv_loc%end), vR_z(:, :, :, iv_loc%beg:iv_loc%end), &
                                           recon_dir, &
                                           is1_loc, is2_loc, is3_loc)
                    else
                        call s_${SCHEME}$ (v_vf(iv_loc%beg:iv_loc%end), &
                                           vL_x(:, :, :, iv_loc%beg:iv_loc%end), vL_y(:, :, :, iv_loc%beg:iv_loc%end), vL_z(:, :, :, :), vR_x(:, :, :, iv_loc%beg:iv_loc%end), vR_y(:, :, :, iv_loc%beg:iv_loc%end), vR_z(:, :, :, :), &
                                           recon_dir, &
                                           is1_loc, is2_loc, is3_loc)
                    end if
                else

                    call s_${SCHEME}$ (v_vf(iv_loc%beg:iv_loc%end), &
                                       vL_x(:, :, :, iv_loc%beg:iv_loc%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv_loc%beg:iv_loc%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                                       recon_dir, &
                                       is1_loc, is2_loc, is3_loc)
                end if
            end if
        #:endfor
    end subroutine s_reconstruct_cell_boundary_values_loc

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
            particle_rad(k, 2) = particle_rad(k, 1)
            particle_pos(k, 1:3, 2) = particle_pos(k, 1:3, 1)
            particle_posPrev(k, 1:3, 2) = particle_posPrev(k, 1:3, 1)
            particle_vel(k, 1:3, 2) = particle_vel(k, 1:3, 1)
            particle_s(k, 1:3, 2) = particle_s(k, 1:3, 1)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_transfer_data_to_tmp_particles

    !> This subroutine transfer data into the actual variables.
    impure subroutine s_transfer_data_to_actual_particles()

        integer :: k

        $:GPU_PARALLEL_LOOP(private='[k]')
        do k = 1, n_el_particles_loc
            particle_rad(k, 1) = particle_rad(k, 2)
            particle_pos(k, 1:3, 1) = particle_pos(k, 1:3, 2)
            particle_posPrev(k, 1:3, 1) = particle_posPrev(k, 1:3, 2)
            particle_vel(k, 1:3, 1) = particle_vel(k, 1:3, 2)
            particle_s(k, 1:3, 1) = particle_s(k, 1:3, 2)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_transfer_data_to_actual_particles

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

    !> The purpose of this procedure is to calculate the gradient from reconstructed states along the x, y and z
        !! @param q Input scalar field
        !! @param dq Output gradient of q
        !! @param dir Gradient spatial direction
    subroutine s_gradient_field(vL_field, vR_field, dq, dir, field_var)

        ! real(stp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:), intent(in) :: q
        real(stp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:), intent(out) :: dq
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vL_field
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vR_field
        integer, intent(in) :: dir, field_var

        integer :: i, j, k
        real(wp) :: mydx

        if (dir == 1) then

            $:GPU_PARALLEL_LOOP(private='[i,j,k,mydx]', collapse=3,copyin='[dir, field_var]')
            do k = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(2)%beg, idwbuff(2)%end
                    do i = idwbuff(1)%beg, idwbuff(1)%end
                        mydx = dx(i)
                        dq(i, j, k) = (vR_field(i, j, k, field_var) - vL_field(i, j, k, field_var))/mydx
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        elseif (dir == 2) then

            $:GPU_PARALLEL_LOOP(private='[i,j,k,mydx]', collapse=3,copyin='[dir, field_var]')
            do k = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    do i = idwbuff(2)%beg, idwbuff(2)%end
                        mydx = dy(i)
                        dq(j, i, k) = (vR_field(i, j, k, field_var) - vL_field(i, j, k, field_var))/mydx
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        elseif (dir == 3) then

            $:GPU_PARALLEL_LOOP(private='[i,j,k,mydx]', collapse=3,copyin='[dir, field_var]')
            do k = idwbuff(1)%beg, idwbuff(1)%end
                do j = idwbuff(2)%beg, idwbuff(2)%end
                    do i = idwbuff(3)%beg, idwbuff(3)%end
                        mydx = dz(i)
                        dq(k, j, i) = (vR_field(i, j, k, field_var) - vL_field(i, j, k, field_var))/mydx
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end if

    end subroutine s_gradient_field

    !> The purpose of this procedure is to compute the Fornberg finite difference weights for derivatives (only done once at start time)
    impure subroutine s_compute_fornberg_fd_weights(npts)

        integer, intent(in) :: npts
        integer :: i, j, k, a, m_order
        integer :: s_idx
        real(wp) :: x0, y0, z0
        real(wp) :: x_stencil(nWeights_grad)
        real(wp) :: c(nWeights_grad, 0:1)

        m_order = 1   ! first derivative

        $:GPU_PARALLEL_LOOP(private='[i,a,x_stencil,c,s_idx,x0]', copyin='[npts,m_order]')
        do i = idwbuff(1)%beg + npts, idwbuff(1)%end - npts
            do a = -npts, npts
                s_idx = a + npts + 1
                x_stencil(s_idx) = x_cc(i + a)
            end do
            x0 = x_cc(i)

            call s_fornberg_weights(x0, x_stencil, nWeights_grad, m_order, c)

            do a = -npts, npts
                s_idx = a + npts + 1
                weights_x_grad(s_idx)%sf(i, 1, 1) = c(s_idx, 1)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[j,a,x_stencil,c,s_idx,y0]', copyin='[npts,m_order]')
        do j = idwbuff(2)%beg + npts, idwbuff(2)%end - npts
            do a = -npts, npts
                s_idx = a + npts + 1
                x_stencil(s_idx) = y_cc(j + a)
            end do
            y0 = y_cc(j)

            call s_fornberg_weights(y0, x_stencil, nWeights_grad, m_order, c)

            do a = -npts, npts
                s_idx = a + npts + 1
                weights_y_grad(s_idx)%sf(j, 1, 1) = c(s_idx, 1)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (num_dims == 3) then

            $:GPU_PARALLEL_LOOP(private='[k,a,x_stencil,c,s_idx,z0]', copyin='[npts,m_order]')
            do k = idwbuff(3)%beg + npts, idwbuff(3)%end - npts
                do a = -npts, npts
                    s_idx = a + npts + 1
                    x_stencil(s_idx) = z_cc(k + a)
                end do
                z0 = z_cc(k)

                call s_fornberg_weights(z0, x_stencil, nWeights_grad, m_order, c)

                do a = -npts, npts
                    s_idx = a + npts + 1
                    weights_z_grad(s_idx)%sf(k, 1, 1) = c(s_idx, 1)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end if

    end subroutine s_compute_fornberg_fd_weights

    !> The purpose of this procedure is to compute the Fornberg finite difference weights on a local stencil
    subroutine s_fornberg_weights(x0, stencil, npts, m_order, coeffs)
        ! $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: npts       ! number of stencil points
        integer, intent(in) :: m_order       ! highest derivative order
        real(wp), intent(in) :: x0     ! evaluation point
        real(wp), intent(in) :: stencil(npts)   ! stencil coordinates
        real(wp), intent(out) :: coeffs(npts, 0:m_order)

        integer :: i, j, k, mn
        real(wp) :: c1, c2, c3, c4, c5

        coeffs = 0.0_wp
        c1 = 1.0_wp
        c4 = stencil(1) - x0
        coeffs(1, 0) = 1.0_wp

        do i = 2, npts
            mn = min(i - 1, m_order)
            c2 = 1.0_wp
            c5 = c4
            c4 = stencil(i) - x0

            do j = 1, i - 1
                c3 = stencil(i) - stencil(j)
                c2 = c2*c3

                if (j == i - 1) then
                    do k = mn, 1, -1
                        coeffs(i, k) = c1*(k*coeffs(i - 1, k - 1) - c5*coeffs(i - 1, k))/c2
                    end do
                    coeffs(i, 0) = -c1*c5*coeffs(i - 1, 0)/c2
                end if

                do k = mn, 1, -1
                    coeffs(j, k) = (c4*coeffs(j, k) - k*coeffs(j, k - 1))/c3
                end do
                coeffs(j, 0) = c4*coeffs(j, 0)/c3
            end do

            c1 = c2
        end do

    end subroutine s_fornberg_weights

    !> The purpose of this procedure is to compute the barycentric weights for interpolation (only done once at start time)
    impure subroutine s_compute_barycentric_weights(npts)

        integer, intent(in) :: npts
        integer :: i, j, k, l, a, b
        real(wp) :: prod_x, prod_y, prod_z, dx_loc, dy_loc, dz_loc

        $:GPU_PARALLEL_LOOP(private='[i,a,b,prod_x,dx_loc]', copyin = '[npts]')
        do i = idwbuff(1)%beg + npts, idwbuff(1)%end - npts
            do a = -npts, npts
                prod_x = 1._wp
                do b = -npts, npts
                    if (a /= b) then
                        dx_loc = x_cc(i + a) - x_cc(i + b)
                        prod_x = prod_x*dx_loc
                    end if
                end do
                weights_x_interp(a + npts + 1)%sf(i, 1, 1) = 1._wp/prod_x
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[j,a,b,prod_y,dy_loc]', copyin = '[npts]')
        do j = idwbuff(2)%beg + npts, idwbuff(2)%end - npts
            do a = -npts, npts
                prod_y = 1._wp
                do b = -npts, npts
                    if (a /= b) then
                        dy_loc = y_cc(j + a) - y_cc(j + b)
                        prod_y = prod_y*dy_loc
                    end if
                end do
                weights_y_interp(a + npts + 1)%sf(j, 1, 1) = 1._wp/prod_y
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (num_dims == 3) then
            $:GPU_PARALLEL_LOOP(private='[k,a,b,prod_z,dz_loc]', copyin = '[npts]')
            do k = idwbuff(3)%beg + npts, idwbuff(3)%end - npts
                do a = -npts, npts
                    prod_z = 1._wp
                    do b = -npts, npts
                        if (a /= b) then
                            dz_loc = z_cc(k + a) - z_cc(k + b)
                            prod_z = prod_z*dz_loc
                        end if
                    end do
                    weights_z_interp(a + npts + 1)%sf(k, 1, 1) = 1._wp/prod_z
                end do
            end do
        end if

    end subroutine s_compute_barycentric_weights

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
                'Vx', 'Vy', 'Vz', 'Fp_x', 'Fp_y', 'Fp_z', &
                'radius'
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
            ! FMT = "(F16.8,I14,8F16.8)"
            FMT = "(F16.8,I14,10F16.8)"
        else
            ! FMT = "(F24.16,I14,8F24.16)"
            FMT = "(F24.16,I14,10F24.16)"
        end if

        ! Cycle through list
        do k = 1, n_el_particles_loc
            write (LAG_EVOL_ID, FMT) &
                qtime, &
                lag_part_id(k, 1), &
                particle_pos(k, 1, 1), &
                particle_pos(k, 2, 1), &
                particle_pos(k, 3, 1), &
                particle_vel(k, 1, 1), &
                particle_vel(k, 2, 1), &
                particle_vel(k, 3, 1), &
                f_p(k, 1), &
                f_p(k, 2), &
                f_p(k, 3), &
                particle_rad(k, 1)
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
                    lag_void_max = max(lag_void_max, 1._wp - q_particles(alphaf_id)%sf(i, j, k))
                    call s_get_char_vol(i, j, k, volcell)
                    if ((1._wp - q_particles(alphaf_id)%sf(i, j, k)) > 5.0d-11) then
                        lag_void_avg = lag_void_avg + (1._wp - q_particles(alphaf_id)%sf(i, j, k))*volcell
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
                if (particle_in_domain_physical(particle_pos(k, 1:3, 1))) then
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
                if (particle_in_domain_physical(particle_pos(k, 1:3, 1))) then
                    MPI_IO_DATA_lag_bubbles(i, 1) = real(lag_part_id(k, 1))
                    MPI_IO_DATA_lag_bubbles(i, 2:4) = particle_pos(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 5:7) = particle_posPrev(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 8:10) = particle_vel(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 11) = particle_rad(k, 1)
                    ! MPI_IO_DATA_lag_bubbles(i, 12) = intfc_vel(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 13) = particle_R0(k)
                    MPI_IO_DATA_lag_bubbles(i, 14) = Rmax_stats_part(k)
                    MPI_IO_DATA_lag_bubbles(i, 15) = Rmin_stats_part(k)
                    ! MPI_IO_DATA_lag_bubbles(i, 16) = bub_dphidt(k)
                    ! MPI_IO_DATA_lag_bubbles(i, 17) = gas_p(k, 1)
                    ! MPI_IO_DATA_lag_bubbles(i, 18) = gas_mv(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 19) = particle_mass(k)
                    ! MPI_IO_DATA_lag_bubbles(i, 20) = gas_betaT(k)
                    ! MPI_IO_DATA_lag_bubbles(i, 21) = gas_betaC(k)
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
            Rmax_stats_part(k) = max(Rmax_stats_part(k), particle_rad(k, 1)/particle_R0(k))
            Rmin_stats_part(k) = min(Rmin_stats_part(k), particle_rad(k, 1)/particle_R0(k))
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
                lag_part_id(k, 1), &
                particle_pos(k, 1, 1), &
                particle_pos(k, 2, 1), &
                particle_pos(k, 3, 1), &
                Rmax_stats_part(k), &
                Rmin_stats_part(k)
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
        Rmax_stats_part(dest) = Rmax_stats_part(src)
        Rmin_stats_part(dest) = Rmin_stats_part(src)
        particle_mass(dest) = particle_mass(src)
        lag_part_id(dest, 1) = lag_part_id(src, 1)
        particle_rad(dest, 1:2) = particle_rad(src, 1:2)
        particle_vel(dest, 1:3, 1:2) = particle_vel(src, 1:3, 1:2)
        particle_s(dest, 1:3, 1:2) = particle_s(src, 1:3, 1:2)
        particle_pos(dest, 1:3, 1:2) = particle_pos(src, 1:3, 1:2)
        particle_posPrev(dest, 1:3, 1:2) = particle_posPrev(src, 1:3, 1:2)
        particle_draddt(dest, 1:lag_num_ts) = particle_draddt(src, 1:lag_num_ts)
        f_p(dest, 1:3) = f_p(src, 1:3)
        particle_dposdt(dest, 1:3, 1:lag_num_ts) = particle_dposdt(src, 1:3, 1:lag_num_ts)
        particle_dveldt(dest, 1:3, 1:lag_num_ts) = particle_dveldt(src, 1:3, 1:lag_num_ts)

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

        do i = 1, nField_vars
            @:DEALLOCATE(field_vars(i)%sf)
        end do
        @:DEALLOCATE(field_vars)

        do i = 1, nWeights_interp
            @:DEALLOCATE(weights_x_interp(i)%sf)
        end do
        @:DEALLOCATE(weights_x_interp)

        do i = 1, nWeights_interp
            @:DEALLOCATE(weights_y_interp(i)%sf)
        end do
        @:DEALLOCATE(weights_y_interp)

        do i = 1, nWeights_interp
            @:DEALLOCATE(weights_z_interp(i)%sf)
        end do
        @:DEALLOCATE(weights_z_interp)

        do i = 1, nWeights_grad
            @:DEALLOCATE(weights_x_grad(i)%sf)
        end do
        @:DEALLOCATE(weights_x_grad)

        do i = 1, nWeights_grad
            @:DEALLOCATE(weights_y_grad(i)%sf)
        end do
        @:DEALLOCATE(weights_y_grad)

        do i = 1, nWeights_grad
            @:DEALLOCATE(weights_z_grad(i)%sf)
        end do
        @:DEALLOCATE(weights_z_grad)

        !Deallocating space
        @:DEALLOCATE(lag_part_id)
        @:DEALLOCATE(particle_R0)
        @:DEALLOCATE(Rmax_stats_part)
        @:DEALLOCATE(Rmin_stats_part)
        @:DEALLOCATE(particle_mass)
        @:DEALLOCATE(particle_rad)
        @:DEALLOCATE(particle_pos)
        @:DEALLOCATE(particle_posPrev)
        @:DEALLOCATE(particle_vel)
        @:DEALLOCATE(particle_s)
        @:DEALLOCATE(particle_draddt)
        @:DEALLOCATE(particle_dposdt)
        @:DEALLOCATE(particle_dveldt)
        @:DEALLOCATE(f_p)

    end subroutine s_finalize_particle_lagrangian_solver

end module m_particles_EL
