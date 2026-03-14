!>
!! @file m_particles_EL.fpp
!! @brief Contains module m_particles_EL

#:include 'macros.fpp'

!> @brief This module is used to to compute the volume-averaged particle model
module m_particles_EL

    use m_global_parameters             !< Definitions of the global parameters

    use m_mpi_proxy                     !< Message passing interface (MPI) module proxy

    use m_particles_EL_kernels            !< Definitions of the kernel functions

    use m_variables_conversion          !< State variables type conversion procedures

    use m_compile_specific

    use m_boundary_common

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_sim_helpers

    use m_helper

    use m_mpi_common

    use m_ibm

    implicit none

    real(wp), save :: next_write_time = 0._wp

    integer, allocatable, dimension(:, :) :: lag_part_id                 !< Global and local IDs
    integer, allocatable, dimension(:) :: gid_to_local
    real(wp), allocatable, dimension(:) :: particle_R0            !< Initial particle radius
    real(wp), allocatable, dimension(:) :: Rmax_stats_part        !< Maximum radius
    real(wp), allocatable, dimension(:) :: Rmin_stats_part        !< Minimum radius
    $:GPU_DECLARE(create='[lag_part_id, gid_to_local, particle_R0, Rmax_stats_part, Rmin_stats_part]')

    real(wp), allocatable, dimension(:) :: particle_mass            !< Particle Mass
    $:GPU_DECLARE(create='[particle_mass]')
    real(wp), allocatable, dimension(:) :: p_AM            !< Particle Added Mass
    $:GPU_DECLARE(create='[p_AM]')

    integer, allocatable, dimension(:) :: p_owner_rank            !< Particle Added Mass
    $:GPU_DECLARE(create='[p_owner_rank]')

    integer, allocatable, dimension(:) :: linked_list            !< particle cell linked list
    $:GPU_DECLARE(create='[linked_list]')

    integer, allocatable, dimension(:, :, :) :: particle_head !< particle heads at each cell
    $:GPU_DECLARE(create='[particle_head]')

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
    type(scalar_field), dimension(:), allocatable :: kahan_comp !< Kahan compensation for q_beta accumulation
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
    integer, parameter :: drhox_id = 4 !< Spatial density gradient in x, y, and z
    integer, parameter :: drhoy_id = 5
    integer, parameter :: drhoz_id = 6
    integer, parameter :: dufx_id = 7 !< Spatial velocity gradient in x, y, and z
    integer, parameter :: dufy_id = 8
    integer, parameter :: dufz_id = 9
    integer, parameter :: dalphafx_id = 10 !< Spatial fluid volume fraction gradient in x, y, and z
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

    $:GPU_DECLARE(create='[Rmax_glb,Rmin_glb,q_particles,kahan_comp,q_particles_idx,field_vars]')
    $:GPU_DECLARE(create='[weights_x_interp,weights_y_interp,weights_z_interp,nWeights_interp]')
    $:GPU_DECLARE(create='[weights_x_grad,weights_y_grad,weights_z_grad,nWeights_grad]')

    !Particle Source terms for fluid coupling
    real(wp), allocatable, dimension(:, :) :: f_p !< force on each particle
    $:GPU_DECLARE(create='[f_p]')

    real(wp), allocatable, dimension(:) :: gSum !< gaussian sum for each particle
    $:GPU_DECLARE(create='[gSum]')

    integer, allocatable :: force_recv_ids(:) !< ids of collision forces received from other ranks
    real(wp), allocatable :: force_recv_vals(:) !< collision forces received from other ranks
    $:GPU_DECLARE(create='[force_recv_ids, force_recv_vals]')

    integer, parameter :: LAG_EVOL_ID = 11 ! File id for lag_bubbles_evol_*.dat
    integer, parameter :: LAG_STATS_ID = 12 ! File id for stats_lag_bubbles_*.dat
    integer, parameter :: LAG_VOID_ID = 13 ! File id for voidfraction.dat

    integer, allocatable, dimension(:) :: keep_bubble
    integer, allocatable, dimension(:, :) :: wrap_bubble_loc, wrap_bubble_dir
    $:GPU_DECLARE(create='[keep_bubble]')
    $:GPU_DECLARE(create='[wrap_bubble_loc, wrap_bubble_dir]')

    integer :: error_flag !Error flag for collisions
    $:GPU_DECLARE(create='[error_flag]')

    integer, parameter :: ncc = 1 !< Number of collisions cells at boundaries
    real(wp) :: eps_overlap = 1.e-12

contains

    !> Initializes the lagrangian subgrid particle solver
        !! @param q_cons_vf Initial conservative variables
    impure subroutine s_initialize_particles_EL_module(q_cons_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type

        integer :: nParticles_glb, i, j, k, nf, l, npts

        !PRIM TO CONS VARIABLES
        real(wp) :: dyn_pres, pi_inf, qv, gamma, pres, T
        real(wp) :: rhou, alpharhou, rho_f, alpharho
        real(wp), dimension(3) :: fluid_vel
        real(wp) :: rhoYks(1:num_species)
        integer :: save_count
        real(wp) :: qtime

        real(wp) :: myR, func_sum
        real(wp), dimension(3) :: myPos, myVel, myForce
        integer, dimension(3) :: cell
        logical :: only_beta = .true.

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

        pcomm_coords(1)%beg = x_cb(-1)
        pcomm_coords(1)%end = x_cb(m)
        $:GPU_UPDATE(device='[pcomm_coords(1)]')
        if (n > 0) then
            pcomm_coords(2)%beg = y_cb(-1)
            pcomm_coords(2)%end = y_cb(n)
            $:GPU_UPDATE(device='[pcomm_coords(2)]')
            if (p > 0) then
                pcomm_coords(3)%beg = z_cb(-1)
                pcomm_coords(3)%end = z_cb(p)
                $:GPU_UPDATE(device='[pcomm_coords(3)]')
            end if
        end if

        pcomm_coords_ghost(1)%beg = x_cb(-1 + ncc)
        pcomm_coords_ghost(1)%end = x_cb(m - ncc)
        $:GPU_UPDATE(device='[pcomm_coords_ghost(1)]')
        if (n > 0) then
            pcomm_coords_ghost(2)%beg = y_cb(-1 + ncc)
            pcomm_coords_ghost(2)%end = y_cb(n - ncc)
            $:GPU_UPDATE(device='[pcomm_coords_ghost(2)]')
            if (p > 0) then
                pcomm_coords_ghost(3)%beg = z_cb(-1 + ncc)
                pcomm_coords_ghost(3)%end = z_cb(p - ncc)
                $:GPU_UPDATE(device='[pcomm_coords_ghost(3)]')
            end if
        end if

        $:GPU_UPDATE(device='[lag_num_ts, q_particles_idx]')

        @:ALLOCATE(q_particles(1:q_particles_idx))
        @:ALLOCATE(kahan_comp(1:q_particles_idx))
        do i = 1, q_particles_idx
            @:ALLOCATE(q_particles(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_particles(i))
            @:ALLOCATE(kahan_comp(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(kahan_comp(i))
        end do

        @:ALLOCATE(field_vars(1:nField_vars))
        do i = 1, nField_vars
            @:ALLOCATE(field_vars(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(field_vars(i))
        end do

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
        @:ALLOCATE(gid_to_local(1:nParticles_glb))
        @:ALLOCATE(particle_R0(1:nParticles_glb))
        @:ALLOCATE(Rmax_stats_part(1:nParticles_glb))
        @:ALLOCATE(Rmin_stats_part(1:nParticles_glb))
        @:ALLOCATE(particle_mass(1:nParticles_glb))
        @:ALLOCATE(p_AM(1:nParticles_glb))
        @:ALLOCATE(p_owner_rank(1:nParticles_glb))
        @:ALLOCATE(particle_rad(1:nParticles_glb, 1:2))
        @:ALLOCATE(particle_pos(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_posPrev(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_vel(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_s(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_draddt(1:nParticles_glb, 1:lag_num_ts))
        @:ALLOCATE(particle_dposdt(1:nParticles_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(particle_dveldt(1:nParticles_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(f_p(1:nParticles_glb, 1:3))
        @:ALLOCATE(gSum(1:nParticles_glb))

        @:ALLOCATE(linked_list(1:nParticles_glb))

        @:ALLOCATE(particle_head(idwbuff(1)%beg:idwbuff(1)%end, &
            idwbuff(2)%beg:idwbuff(2)%end, &
            idwbuff(3)%beg:idwbuff(3)%end))

        @:ALLOCATE(force_recv_ids(1:lag_params%nParticles_glb))
        @:ALLOCATE(force_recv_vals(1:3*lag_params%nParticles_glb))

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

        ! Allocate cell list arrays for atomic-free Gaussian smearing
        @:ALLOCATE(cell_list_start(0:m, 0:n, 0:p))
        @:ALLOCATE(cell_list_count(0:m, 0:n, 0:p))
        @:ALLOCATE(cell_list_idx(1:lag_params%nParticles_glb))

        call s_read_input_particles(q_cons_vf, bc_type)

        call s_reset_cell_vars()

        $:GPU_PARALLEL_LOOP(private='[k,cell,myR,myPos,myVel,myForce,func_sum]',copyin='[only_beta]')
        do k = 1, n_el_particles_loc

            cell = fd_number - buff_size
            call s_locate_cell(particle_pos(k, 1:3, 1), cell, particle_s(k, 1:3, 1))

            myR = particle_R0(k)
            myPos = particle_pos(k, 1:3, 1)
            myVel = particle_vel(k, 1:3, 1)
            myForce = f_p(k, :)
            !Compute the total gaussian contribution for each particle for normalization
            call s_compute_gaussian_contribution(myR, myPos, cell, func_sum)
            gSum(k) = func_sum

            call s_gaussian_atomic(myR, myVel, myPos, myForce, func_sum, cell, q_particles, only_beta)

        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_finalize_beta_field(bc_type, only_beta)

        npts = (nWeights_interp - 1)/2
        call s_compute_barycentric_weights(npts) !For interpolation

        npts = (nWeights_grad - 1)/2
        call s_compute_fornberg_fd_weights(npts) !For finite differences

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
    impure subroutine s_read_input_particles(q_cons_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type

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
                        n_el_particles_loc = particle_id              ! local number of particles
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

        $:GPU_UPDATE(device='[lag_part_id,particle_R0,Rmax_stats_part,Rmin_stats_part,particle_mass, &
            & f_p,p_AM,p_owner_rank,gid_to_local, &
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
        real(wp) :: massflag, heatflag, Re_trans, Im_trans, myR, func_sum
        real(wp), dimension(3) :: myPos, myVel

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
        p_AM(part_id) = 0._wp
        p_owner_rank(part_id) = proc_rank
        gid_to_local(part_id) = -1

        if (cyl_coord .and. p == 0) then
            particle_pos(part_id, 2, 1) = sqrt(particle_pos(part_id, 2, 1)**2._wp + &
                                               particle_pos(part_id, 3, 1)**2._wp)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            particle_pos(part_id, 3, 1) = atan2(inputPart(3), inputPart(2))
            particle_posPrev(part_id, 1:3, 1) = particle_pos(part_id, 1:3, 1)
        end if

        cell = fd_number - buff_size
        call s_locate_cell(particle_pos(part_id, 1:3, 1), cell, particle_s(part_id, 1:3, 1))

        ! Check if the particle is located in the ghost cell of a symmetric, or wall boundary
        if ((any(bc_x%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(1) < 0) .or. &
            (any(bc_x%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(1) > m) .or. &
            (any(bc_y%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(2) < 0) .or. &
            (any(bc_y%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(2) > n)) then
            call s_mpi_abort("Lagrange particle is in the ghost cells of a symmetric or wall boundary.")
        end if

        if (p > 0) then
            if ((any(bc_z%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(3) < 0) .or. &
                (any(bc_z%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(3) > p)) then
                call s_mpi_abort("Lagrange particle is in the ghost cells of a symmetric or wall boundary.")
            end if
        end if

        ! Initial particle mass
        volparticle = 4._wp/3._wp*pi*particle_R0(part_id)**3 ! volume
        particle_mass(part_id) = volparticle*rho0ref_particle ! mass
        if (particle_mass(part_id) <= 0._wp) then
            call s_mpi_abort("The initial particle mass is negative or zero. Check the particle file.")
        end if

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
    subroutine s_compute_particle_EL_dynamics(q_prim_vf, bc_type, stage, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, rhs_vf)
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type
        type(scalar_field), dimension(sys_size), intent(in) :: rhs_vf
        integer, intent(in) :: stage
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vL_x, vL_y, vL_z
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vR_x, vR_y, vR_z

        integer, dimension(3) :: cell, cellijk
        real(wp) :: myMass, myR, myBeta_c, myBeta_t, myR0, myRe, mydrhodt, myVolumeFrac, myGamma, rmass_add, func_sum
        real(wp), dimension(3) :: myVel, myPos, force_vec, s_cell
        logical :: only_beta = .false.

        integer :: k, l, i, j

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

            do l = 1, num_dims
                if (l == 1) then
                    call s_gradient_field(vL_x, vR_x, field_vars(drhox_id)%sf, l, 1)
                    ! call s_gradient_field(vL_x, vR_x, field_vars(dufx_id)%sf, l, momxb)
                elseif (l == 2) then
                    call s_gradient_field(vL_y, vR_y, field_vars(drhoy_id)%sf, l, 1)
                    ! call s_gradient_field(vL_y, vR_y, field_vars(dufy_id)%sf, l, momxb+1)
                elseif (l == 3) then
                    call s_gradient_field(vL_z, vR_z, field_vars(drhoz_id)%sf, l, 1)
                    ! call s_gradient_field(vL_z, vR_z, field_vars(dufz_id)%sf, l, momxb+2)
                end if
            end do

        end if

        myGamma = (1._wp/fluid_pp(1)%gamma) + 1._wp
        myRe = 1.845e-5_wp !fluid_pp(1)%Re(1) !Need a viscosity model for when modeling inviscid eulerian fluid !< Dynamic viscosity

        call nvtxStartRange("LAGRANGE-PARTICLE-DYNAMICS")

        !> Compute Fluid-Particle Forces (drag/pressure/added mass) and convert to particle acceleration
        $:GPU_PARALLEL_LOOP(private='[i,k,l,cell,s_cell,myMass,myR,myR0,myPos,myVel,myVolumeFrac,force_vec,rmass_add,func_sum,mydrhodt]',&
        & copyin='[stage, myGamma, myRe, only_beta]')
        do k = 1, n_el_particles_loc

            f_p(k, :) = 0._wp
            p_owner_rank(k) = proc_rank

            s_cell = particle_s(k, 1:3, 2)
            cell = int(s_cell(:))
            do i = 1, num_dims
                if (s_cell(i) < 0._wp) cell(i) = cell(i) - 1
            end do

            ! Current particle state
            myMass = particle_mass(k)
            myR = particle_rad(k, 2)
            myR0 = particle_R0(k)
            myPos = particle_pos(k, :, 2)
            myVel = particle_vel(k, :, 2)
            myVolumeFrac = 1._wp - q_particles(alphaf_id)%sf(cell(1), cell(2), cell(3))
            mydrhodt = rhs_vf(1)%sf(cell(1), cell(2), cell(3))

            particle_dposdt(k, :, stage) = 0._wp
            particle_dveldt(k, :, stage) = 0._wp
            particle_draddt(k, stage) = 0._wp

            call s_get_particle_force(myPos, myR, myVel, myMass, myRe, myGamma, myVolumeFrac, mydrhodt, cell, &
                                      q_prim_vf, field_vars, weights_x_interp, weights_y_interp, weights_z_interp, &
                                      force_vec, rmass_add)

            p_AM(k) = rMass_add
            f_p(k, :) = f_p(k, :) + force_vec(:)

            if (.not. lag_params%collision_force) then
                myMass = particle_mass(k) + p_AM(k)
                myVel = particle_vel(k, :, 2)
                do l = 1, num_dims
                    particle_dposdt(k, l, stage) = myVel(l)
                    particle_dveldt(k, l, stage) = f_p(k, l)/myMass
                    particle_draddt(k, stage) = 0._wp
                end do
            end if

            if (lag_params%solver_approach == 2) then
                func_sum = gSum(k)
                call s_gaussian_atomic(myR, myVel, myPos, force_vec, func_sum, cell, q_particles, only_beta)
            end if

        end do
        $:END_GPU_PARALLEL_LOOP()

        if (lag_params%solver_approach == 2) then
            call s_finalize_beta_field(bc_type, only_beta)
        end if

        call nvtxStartRange("LAGRANGE-PARTICLE-COLLISIONS")
        if (lag_params%collision_force) then
            !> Compute Particle-Particle collision forces
            call s_compute_particle_EL_collisions(stage, bc_type)

            $:GPU_PARALLEL_LOOP(private='[k,l,myMass,myVel]')
            do k = 1, n_el_particles_loc
                myMass = particle_mass(k) + p_AM(k)
                myVel = particle_vel(k, :, 2)
                do l = 1, num_dims
                    particle_dposdt(k, l, stage) = myVel(l)
                    particle_dveldt(k, l, stage) = f_p(k, l)/myMass
                    particle_draddt(k, stage) = 0._wp
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end if
        call nvtxEndRange

        call nvtxEndRange

    end subroutine s_compute_particle_EL_dynamics

    !>  Contains the particle collision force computation.
    subroutine s_compute_particle_EL_collisions(stage, bc_type)

        integer, intent(in) :: stage
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type
        integer, dimension(3) :: cell
        real(wp), dimension(3) :: s_cell
        integer, dimension(3) :: cellaux
        integer :: i, k, l, q, ip, jp, kp, ii, jj, kk
        logical :: celloutside
        real(wp) :: pidtksp2, ksp, nu1, nu2, Rp1, Rp2, E1, E2, Estar, cor, rmag, Rstar, dij, eta_n, kappa_n, mp1, mp2, dt_loc
        real(wp), dimension(3) :: xp1, xp2, vp1, vp2, v_rel, rpij, nij, vnij, Fnpp_ij, force_vec
        integer :: kpz
        integer :: total_recv
        integer :: glb_id, count
        integer :: n_el_particles_loc_before_ghost

        if (num_procs > 1) then
            n_el_particles_loc_before_ghost = n_el_particles_loc
            call s_reset_force_buffers()
            call s_add_ghost_particles()
        end if

        ! if (lag_num_ts == 1) then
        !     dt_loc = dt
        ! elseif (lag_num_ts == 2) then
        !     if (stage == 1) then
        !         dt_loc = dt
        !     elseif (stage == 2) then
        !         dt_loc = dt/2._wp
        !     end if
        ! elseif (lag_num_ts == 3) then
        !     if (stage == 1) then
        !         dt_loc = dt
        !     elseif (stage == 2) then
        !         dt_loc = dt/4._wp
        !     elseif (stage == 3) then
        !         dt_loc = (2._wp/3._wp)*dt
        !     end if
        ! end if

        kpz = 0
        if (num_dims == 3) kpz = 1

        ksp = 10._wp
        nu1 = 0.35_wp
        nu2 = 0.35_wp
        E1 = 1.e9_wp
        E2 = 1.e9_wp
        cor = 0.7_wp

        pidtksp2 = (pi**2)/((dt*ksp)**2)

        Estar = 1._wp/(((1._wp - nu1**2)/E1) + ((1._wp - nu2**2)/E2))
        Estar = (4._wp/3._wp)*Estar

        call s_reset_linked_list()

        call nvtxStartRange("LAGRANGE-PARTICLE-COLLISIONS")
        error_flag = 0
        $:GPU_UPDATE(device='[error_flag]')

        $:GPU_PARALLEL_LOOP(private='[i,k,cell,ip,jp,kp,Rp1,xp1,mp1,vp1,kk,jj,ii,cellaux,q,Rp2,xp2,mp2,vp2,v_rel,Rstar,rpij,rmag,nij,vnij,dij,kappa_n,eta_n,Fnpp_ij,force_vec,s_cell,celloutside,count]',&
        & copyin='[ksp,nu1,nu2,E1,E2,cor,pidtksp2,Estar,kpz]')
        do k = 1, n_el_particles_loc

            if (.not. particle_in_domain_physical(particle_pos(k, 1:3, 2))) then
                cycle
            end if

            s_cell = particle_s(k, 1:3, 2)
            cell = int(s_cell(:))
            do i = 1, num_dims
                if (s_cell(i) < 0._wp) cell(i) = cell(i) - 1
            end do

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
                            q = particle_head(ii, jj, kk)
                            ! Traverse linked list in that cell

                            count = 0
                            do while (q /= -1)

                                count = count + 1
                                if (count > n_el_particles_loc) then
                                    $:GPU_ATOMIC(atomic='write')
                                    error_flag = 1
                                    exit
                                end if

                                if (lag_part_id(q, 1) > lag_part_id(k, 1)) then

                                    Rp2 = particle_rad(q, 2)
                                    xp2 = particle_pos(q, :, 2)
                                    mp2 = particle_mass(q)
                                    vp2 = particle_vel(q, :, 2)
                                    v_rel = vp2 - vp1

                                    Rstar = (Rp1*Rp2)/(Rp1 + Rp2)
                                    rpij = xp2 - xp1
                                    rmag = sqrt(rpij(1)**2 + rpij(2)**2 + rpij(3)**2)
                                    rmag = max(rmag, eps_overlap)
                                    nij = rpij/rmag
                                    vnij = dot_product(v_rel, nij)*nij
                                    dij = (Rp1 + Rp2) - rmag

                                    if (dij > 0._wp) then

                                        kappa_n = min((pidtksp2*mp1), (pidtksp2*mp2), (Estar*sqrt(Rstar)*sqrt(abs(dij))))

                                        eta_n = ((-2._wp*sqrt(kappa_n)*log(cor))/sqrt((log(cor))**2 + pi**2))*(1._wp/sqrt((1._wp/mp1) + (1._wp/mp2)))

                                        Fnpp_ij = -kappa_n*dij*nij - eta_n*vnij

                                        f_p(k, :) = f_p(k, :) + Fnpp_ij

                                        if (p_owner_rank(q) == proc_rank) then
                                            ! f_p(q, :) = f_p(q, :) - Fnpp_ij

                                            $:GPU_ATOMIC(atomic='update')
                                            f_p(q, 1) = f_p(q, 1) - Fnpp_ij(1)

                                            $:GPU_ATOMIC(atomic='update')
                                            f_p(q, 2) = f_p(q, 2) - Fnpp_ij(2)

                                            $:GPU_ATOMIC(atomic='update')
                                            f_p(q, 3) = f_p(q, 3) - Fnpp_ij(3)

                                        else
                                            call s_add_force_to_send_buffer(p_owner_rank(q), lag_part_id(q, 1), -Fnpp_ij)
                                        end if
                                    end if
                                end if

                                q = linked_list(q)

                            end do
                        end if

                    end do
                end do
            end do

            !>Check each local particle for wall collisions

            call s_compute_wall_collisions(xp1, vp1, Rp1, mp1, Estar, pidtksp2, cor, force_vec)
            f_p(k, :) = f_p(k, :) + force_vec

        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

        $:GPU_UPDATE(host='[error_flag]')
        if (error_flag == 1) then
            call s_mpi_abort("Linked list infinite loop detected")
        end if

        if (num_procs > 1) then

            n_el_particles_loc = n_el_particles_loc_before_ghost
            $:GPU_UPDATE(device='[n_el_particles_loc]')

            total_recv = 0
            force_recv_ids = 0
            force_recv_vals = 0.

            call s_transfer_collision_forces(total_recv, force_recv_ids, force_recv_vals)

            $:GPU_UPDATE(device = '[force_recv_ids,force_recv_vals]')

            $:GPU_PARALLEL_LOOP(private='[i,k]',copyin = '[total_recv]')
            do i = 1, total_recv
                k = gid_to_local(force_recv_ids(i))
                if (k > 0) then
                    $:GPU_ATOMIC(atomic='update')
                    f_p(k, 1) = f_p(k, 1) + force_recv_vals(3*(i - 1) + 1)

                    $:GPU_ATOMIC(atomic='update')
                    f_p(k, 2) = f_p(k, 2) + force_recv_vals(3*(i - 1) + 2)

                    $:GPU_ATOMIC(atomic='update')
                    f_p(k, 3) = f_p(k, 3) + force_recv_vals(3*(i - 1) + 3)

                end if
            end do
            $:END_GPU_PARALLEL_LOOP()

        end if

    end subroutine s_compute_particle_EL_collisions

    !> This subroutine checks for particles at solid walls to compute a collision force
    subroutine s_compute_wall_collisions(pos, vel, rad, mass, Es, pidtksp, core, wcol_force)
        $:GPU_ROUTINE(function_name='s_compute_wall_collisions',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), dimension(3), intent(in) :: pos, vel
        real(wp), intent(in) :: rad, mass, Es, pidtksp, core
        real(wp), dimension(3), intent(inout) :: wcol_force

        real(wp) :: dij

        wcol_force = 0._wp

        ! Check for particles at solid boundaries
        if (any(bc_x%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then
            dij = rad - (pos(1) - x_cb(-1))

            if (dij > 0._wp) then
                call s_compute_wall_collision_force(dij, vel, rad, mass, Es, pidtksp, core, 1, 1._wp, wcol_force)
            end if
        end if

        if (any(bc_x%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then

            dij = rad - (x_cb(m) - pos(1))

            if (dij > 0._wp) then
                call s_compute_wall_collision_force(dij, vel, rad, mass, Es, pidtksp, core, 1, -1._wp, wcol_force)
            end if
        end if

        if (any(bc_y%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then

            dij = rad - (pos(2) - y_cb(-1))

            if (dij > 0._wp) then
                call s_compute_wall_collision_force(dij, vel, rad, mass, Es, pidtksp, core, 2, 1._wp, wcol_force)
            end if
        end if

        if (any(bc_y%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then

            dij = rad - (y_cb(n) - pos(2))

            if (dij > 0._wp) then
                call s_compute_wall_collision_force(dij, vel, rad, mass, Es, pidtksp, core, 2, -1._wp, wcol_force)
            end if
        end if

        if (p > 0) then
            if (any(bc_z%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then

                dij = rad - (pos(3) - z_cb(-1))

                if (dij > 0._wp) then
                    call s_compute_wall_collision_force(dij, vel, rad, mass, Es, pidtksp, core, 3, 1._wp, wcol_force)
                end if
            end if

            if (any(bc_z%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then

                dij = rad - (z_cb(p) - pos(3))

                if (dij > 0._wp) then
                    call s_compute_wall_collision_force(dij, vel, rad, mass, Es, pidtksp, core, 3, -1._wp, wcol_force)
                end if
            end if

        end if

    end subroutine s_compute_wall_collisions

    !> This subroutine computes the collision force with a solid wall
    subroutine s_compute_wall_collision_force(dij, vel, rad, mass, Es, pidtksp, core, dir, normal, wcol_force)
        $:GPU_ROUTINE(function_name='s_compute_wall_collision_force',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), dimension(3), intent(in) :: vel
        real(wp), intent(in) :: dij, rad, mass, Es, pidtksp, core, normal
        integer, intent(in) :: dir
        real(wp), dimension(3), intent(inout) :: wcol_force

        real(wp), dimension(3) :: nij, v_rel, vnij
        real(wp) :: kappa_n, eta_n

        ! Normal points away from wall (into domain)
        nij = 0._wp
        nij(dir) = normal

        ! Relative velocity (wall has zero velocity)
        v_rel = vel
        vnij = dot_product(v_rel, nij)*nij

        ! Wall has infinite mass so use mp1 only
        kappa_n = min((pidtksp*mass), (Es*sqrt(rad)*sqrt(abs(dij))))

        eta_n = ((-2._wp*sqrt(kappa_n)*log(core))/sqrt((log(core))**2 + pi**2)) &
                *(1._wp/sqrt(1._wp/mass))

        wcol_force = wcol_force + (kappa_n*dij*nij - eta_n*vnij)

    end subroutine s_compute_wall_collision_force

    !> This subroutine adds temporary ghost particles for collision purposes
    subroutine s_add_ghost_particles()

        integer :: k, i, q
        integer :: patch_id, newBubs
        integer, dimension(3) :: cell
        logical :: inc_ghost = .true.

        call nvtxStartRange("LAG-GHOSTADD")
        call nvtxStartRange("LAG-GHOSTADD-DEV2HOST")
        $:GPU_UPDATE(host='[p_owner_rank, particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
            & lag_part_id, particle_rad, &
            & particle_pos, particle_posPrev, particle_vel, particle_s, particle_draddt, &
            & particle_dposdt, particle_dveldt, n_el_particles_loc, &
            & wrap_bubble_dir, wrap_bubble_loc]')
        call nvtxEndRange

        ! Handle MPI transfer of particles going to another processor's local domain
        if (num_procs > 1) then
            call nvtxStartRange("LAG-GHOSTADD-TRANSFER-LIST")
            call s_add_particles_to_transfer_list(n_el_particles_loc, particle_pos(:, :, 2), particle_posPrev(:, :, 2), inc_ghost)
            call nvtxEndRange

            call nvtxStartRange("LAG-GHOSTADD-SENDRECV")
            call s_mpi_sendrecv_solid_particles(p_owner_rank, particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
                                                lag_part_id, &
                                                particle_rad, particle_pos, particle_posPrev, particle_vel, &
                                                particle_s, particle_draddt, particle_dposdt, particle_dveldt, lag_num_ts, n_el_particles_loc, &
                                                2)
            call nvtxEndRange
        end if

        call nvtxStartRange("LAG-GHOSTADD-HOST2DEV")
        $:GPU_UPDATE(device='[p_owner_rank,particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
            & lag_part_id, particle_rad, &
            & particle_pos, particle_posPrev, particle_vel, particle_s, particle_draddt, &
            & particle_dposdt, particle_dveldt, n_el_particles_loc]')
        call nvtxEndRange

        call nvtxEndRange ! LAG-GHOSTADD

        $:GPU_PARALLEL_LOOP(private='[k,cell]')
        do k = 1, n_el_particles_loc

            cell = fd_number - buff_size
            call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))

        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_add_ghost_particles

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

        !Spatial derivative of the fluid volume fraction and eulerian particle momentum fields.
        do l = 1, num_dims
            call s_gradient_dir_fornberg(q_particles(alphaf_id)%sf, field_vars(dalphafx_id + l - 1)%sf, l)
            call s_gradient_dir_fornberg(q_particles(alphaupx_id + l - 1)%sf, field_vars(dalphap_upx_id + l - 1)%sf, l)
        end do

        !> Apply particle sources to the Eulerian RHS
        $:GPU_PARALLEL_LOOP(private='[i,j,k,alpha_f,dalphapdt,udot_gradalpha]', collapse=3)
        do k = idwint(3)%beg, idwint(3)%end
            do j = idwint(2)%beg, idwint(2)%end
                do i = idwint(1)%beg, idwint(1)%end
                    if (q_particles(alphaf_id)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then

                        alpha_f = q_particles(alphaf_id)%sf(i, j, k)

                        dalphapdt = 0._wp
                        udot_gradalpha = 0._wp
                        do l = 1, num_dims
                            dalphapdt = dalphapdt + field_vars(dalphap_upx_id + l - 1)%sf(i, j, k)
                            udot_gradalpha = udot_gradalpha + q_prim_vf(momxb + l - 1)%sf(i, j, k)*field_vars(dalphafx_id + l - 1)%sf(i, j, k)
                        end do
                        dalphapdt = -dalphapdt
                        !Add any contribution to dalphapdt from particles growing or shrinking

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
    subroutine s_reset_linked_list()

        integer :: j, k, l

        $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    particle_head(j, k, l) = -1
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[k]')
        do k = 1, n_el_particles_loc
            linked_list(k) = -1
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_build_linked_list()

    end subroutine s_reset_linked_list

    !>  The purpose of this subroutine is to smear the effect of the particles in the Eulerian framework
    subroutine s_reset_cell_vars()

        integer :: i, j, k, l

        $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
        do i = 1, max(nField_vars, q_particles_idx)  ! outermost is largest of the i-like dims
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        ! Zero field_vars if i <= nField_vars
                        if (i <= nField_vars) field_vars(i)%sf(j, k, l) = 0._wp
                        ! Zero q_particles if i <= q_particles_idx
                        if (i <= q_particles_idx) then
                            q_particles(i)%sf(j, k, l) = 0._wp
                            kahan_comp(i)%sf(j, k, l) = 0._wp
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_reset_cell_vars

    subroutine s_finalize_beta_field(bc_type, onlyBeta)

        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type
        integer :: j, k, l
        logical, intent(in) :: onlyBeta

        call nvtxStartRange("PARTICLES-LAGRANGE-BETA-COMM")
        if (onlyBeta) then
            call s_populate_beta_buffers(q_particles, bc_type, 1)
        else
            call s_populate_beta_buffers(q_particles, bc_type, q_particles_idx)
        end if
        call nvtxEndRange

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

    end subroutine s_finalize_beta_field

    subroutine s_build_linked_list()
        integer :: k, glb_id, i
        integer, dimension(3) :: cell
        real(wp), dimension(3) :: s_cell
        logical :: celloutside

        $:GPU_PARALLEL_LOOP(private='[i,k,cell,s_cell,glb_id,celloutside]')
        do k = 1, n_el_particles_loc

            glb_id = lag_part_id(k, 1)
            gid_to_local(glb_id) = k

            s_cell = particle_s(k, 1:3, 2)
            cell = int(s_cell(:))
            do i = 1, num_dims
                if (s_cell(i) < 0._wp) cell(i) = cell(i) - 1
            end do

            call s_check_celloutside_wbuff(cell, celloutside)

            if (.not. celloutside) then
              !!!!! Particle linked list building
                $:GPU_ATOMIC(atomic='capture')
                linked_list(k) = particle_head(cell(1), cell(2), cell(3))
                particle_head(cell(1), cell(2), cell(3)) = k
                $:END_GPU_ATOMIC_CAPTURE()
            end if

        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_build_linked_list

    !>  This subroutine updates the Lagrange variables using the tvd RK time steppers.
        !!      The time derivative of the particle variables must be stored at every stage to avoid precision errors.
        !! @param stage Current tvd RK stage
    impure subroutine s_update_lagrange_particles_tdv_rk(q_prim_vf, bc_type, stage)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type
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
            if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage, bc_type)
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

                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage, bc_type)

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
                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage, bc_type)
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

                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage, bc_type)

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

                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage, bc_type)

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
                if (moving_lag_particles) call s_enforce_EL_particles_boundary_conditions(q_prim_vf, stage, bc_type)
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
    impure subroutine s_enforce_EL_particles_boundary_conditions(q_prim_vf, nstage, bc_type)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type
        integer, intent(in) :: nstage
        real(wp) :: offset

        integer :: k, i, q
        integer :: patch_id, newBubs, new_idx
        integer, dimension(3) :: cell
        logical :: inc_ghost = .false.
        real(wp) :: myR, func_sum
        real(wp), dimension(3) :: myPos, myVel, myForce
        logical :: only_beta = .true.

        call nvtxStartRange("LAG-BC")
        call nvtxStartRange("LAG-BC-DEV2HOST")
        $:GPU_UPDATE(host='[p_owner_rank, particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
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
            call s_mpi_sendrecv_solid_particles(p_owner_rank, particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
                                                lag_part_id, &
                                                particle_rad, particle_pos, particle_posPrev, particle_vel, &
                                                particle_s, particle_draddt, particle_dposdt, particle_dveldt, lag_num_ts, n_el_particles_loc, &
                                                2)
            call nvtxEndRange
        end if

        call nvtxStartRange("LAG-BC-HOST2DEV")
        $:GPU_UPDATE(device='[p_owner_rank, particle_R0, Rmax_stats_part, Rmin_stats_part, particle_mass, f_p, &
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
                .and. particle_pos(k, 1, 2) < x_cb(-1) + eps_overlap*particle_rad(k, 2)) then
                particle_pos(k, 1, 2) = x_cb(-1) + eps_overlap*particle_rad(k, 2)
                if (nstage == lag_num_ts) then
                    particle_pos(k, 1, 1) = particle_pos(k, 1, 2)
                end if
            elseif (any(bc_x%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                    .and. particle_pos(k, 1, 2) > x_cb(m) - eps_overlap*particle_rad(k, 2)) then
                particle_pos(k, 1, 2) = x_cb(m) - eps_overlap*particle_rad(k, 2)
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
            elseif (particle_pos(k, 1, 2) >= x_cb(m)) then
                keep_bubble(k) = 0
            elseif (particle_pos(k, 1, 2) < x_cb(-1)) then
                keep_bubble(k) = 0
            end if

            if (any(bc_y%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                .and. particle_pos(k, 2, 2) < y_cb(-1) + eps_overlap*particle_rad(k, 2)) then
                particle_pos(k, 2, 2) = y_cb(-1) + eps_overlap*particle_rad(k, 2)
                if (nstage == lag_num_ts) then
                    particle_pos(k, 2, 1) = particle_pos(k, 2, 2)
                end if
            else if (any(bc_y%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                     .and. particle_pos(k, 2, 2) > y_cb(n) - eps_overlap*particle_rad(k, 2)) then
                particle_pos(k, 2, 2) = y_cb(n) - eps_overlap*particle_rad(k, 2)
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
            elseif (particle_pos(k, 2, 2) >= y_cb(n)) then
                keep_bubble(k) = 0
            elseif (particle_pos(k, 2, 2) < y_cb(-1)) then
                keep_bubble(k) = 0
            end if

            if (p > 0) then
                if (any(bc_z%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                    .and. particle_pos(k, 3, 2) < z_cb(-1) + eps_overlap*particle_rad(k, 2)) then
                    particle_pos(k, 3, 2) = z_cb(-1) + eps_overlap*particle_rad(k, 2)
                    if (nstage == lag_num_ts) then
                        particle_pos(k, 3, 1) = particle_pos(k, 3, 2)
                    end if
                else if (any(bc_z%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) &
                         .and. particle_pos(k, 3, 2) > z_cb(p) - eps_overlap*particle_rad(k, 2)) then
                    particle_pos(k, 3, 2) = z_cb(p) - eps_overlap*particle_rad(k, 2)
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
                elseif (particle_pos(k, 3, 2) >= z_cb(p)) then
                    keep_bubble(k) = 0
                elseif (particle_pos(k, 3, 2) < z_cb(-1)) then
                    keep_bubble(k) = 0
                end if
            end if

            if (keep_bubble(k) == 1) then
                ! Remove bubbles that are no longer in a liquid
                cell = fd_number - buff_size
                call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))

                if (q_prim_vf(advxb)%sf(cell(1), cell(2), cell(3)) < (1._wp - lag_params%valmaxvoid)) then
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

        call s_reset_cell_vars()

        $:GPU_PARALLEL_LOOP(private='[cell,myR,myPos,myVel,myForce,func_sum]',copyin='[only_beta]')
        do k = 1, n_el_particles_loc
            myR = particle_rad(k, 2)
            myPos = particle_pos(k, 1:3, 2)
            myVel = particle_vel(k, 1:3, 2)
            myForce = f_p(k, :)

            cell = fd_number - buff_size
            call s_locate_cell(particle_pos(k, 1:3, 2), cell, particle_s(k, 1:3, 2))

            !Compute the total gaussian contribution for each particle for normalization
            call s_compute_gaussian_contribution(myR, myPos, cell, func_sum)
            gSum(k) = func_sum

            call s_gaussian_atomic(myR, myVel, myPos, myForce, func_sum, cell, q_particles, only_beta)

        end do

        ! Update void fraction and communicate buffers
        call s_finalize_beta_field(bc_type, only_beta)

        call nvtxEndRange ! LAG-BC

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
            particle_rad(k, 2) = particle_rad(k, 1)
            particle_pos(k, 1:3, 2) = particle_pos(k, 1:3, 1)
            particle_posPrev(k, 1:3, 2) = particle_posPrev(k, 1:3, 1)
            particle_vel(k, 1:3, 2) = particle_vel(k, 1:3, 1)
            particle_s(k, 1:3, 2) = particle_s(k, 1:3, 1)
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

    !> The purpose of this procedure is to calculate the gradient from reconstructed states along the x, y and z
        !! @param vL_field left edge reconstructed values
        !! @param vR_field right edge reconstructed values
        !! @param dq Output gradient of q
        !! @param dir Gradient spatial direction
        !! @param field_var variable index for reconstructed states
    subroutine s_gradient_field(vL_field, vR_field, dq, dir, field_var)

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

    !> The purpose of this procedure is to calculate the gradient of a scalar field along the x, y and z directions using Fornberg's method
    !! @param q Input scalar field
    !! @param dq Output gradient of q
    !! @param dir Gradient spatial direction
    subroutine s_gradient_dir_fornberg(q, dq, dir)

        real(stp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:), intent(in) :: q
        real(stp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:), intent(out) :: dq
        integer, intent(in) :: dir

        integer :: i, j, k, a, npts, s_idx

        npts = (nWeights_grad - 1)/2

        if (dir == 1) then
            $:GPU_PARALLEL_LOOP(private='[i,j,k,s_idx,a]', collapse=3,copyin='[npts]')
            do k = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(2)%beg, idwbuff(2)%end
                    do i = idwbuff(1)%beg + 2, idwbuff(1)%end - 2
                        dq(i, j, k) = 0._wp
                        do a = -npts, npts
                            s_idx = a + npts + 1
                            dq(i, j, k) = dq(i, j, k) + weights_x_grad(s_idx)%sf(i, 1, 1)*q(i + a, j, k)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (dir == 2) then
            $:GPU_PARALLEL_LOOP(private='[i,j,k,s_idx,a]', collapse=3,copyin='[npts]')
            do k = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do i = idwbuff(1)%beg, idwbuff(1)%end
                        dq(i, j, k) = 0._wp
                        do a = -npts, npts
                            s_idx = a + npts + 1
                            dq(i, j, k) = dq(i, j, k) + weights_y_grad(s_idx)%sf(j, 1, 1)*q(i, j + a, k)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (dir == 3) then
            $:GPU_PARALLEL_LOOP(private='[i,j,k,s_idx,a]', collapse=3,copyin='[npts]')
            do k = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do j = idwbuff(2)%beg, idwbuff(2)%end
                    do i = idwbuff(1)%beg, idwbuff(1)%end
                        dq(i, j, k) = 0._wp
                        do a = -npts, npts
                            s_idx = a + npts + 1
                            dq(i, j, k) = dq(i, j, k) + weights_z_grad(s_idx)%sf(k, 1, 1)*q(i, j, k + a)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_gradient_dir_fornberg

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
            @:DEALLOCATE(kahan_comp(i)%sf)
        end do
        @:DEALLOCATE(q_particles)
        @:DEALLOCATE(kahan_comp)

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

        !Deallocating space
        @:DEALLOCATE(lag_part_id)
        @:DEALLOCATE(gid_to_local)
        @:DEALLOCATE(particle_R0)
        @:DEALLOCATE(Rmax_stats_part)
        @:DEALLOCATE(Rmin_stats_part)
        @:DEALLOCATE(particle_mass)
        @:DEALLOCATE(p_AM)
        @:DEALLOCATE(p_owner_rank)
        @:DEALLOCATE(particle_rad)
        @:DEALLOCATE(particle_pos)
        @:DEALLOCATE(particle_posPrev)
        @:DEALLOCATE(particle_vel)
        @:DEALLOCATE(particle_s)
        @:DEALLOCATE(particle_draddt)
        @:DEALLOCATE(particle_dposdt)
        @:DEALLOCATE(particle_dveldt)
        @:DEALLOCATE(f_p)
        @:DEALLOCATE(gSum)

        @:DEALLOCATE(force_recv_ids)
        @:DEALLOCATE(force_recv_vals)

        @:DEALLOCATE(keep_bubble)
        @:DEALLOCATE(wrap_bubble_loc, wrap_bubble_dir)

        @:DEALLOCATE(linked_list)
        @:DEALLOCATE(particle_head)

        ! Deallocate cell list arrays
        @:DEALLOCATE(cell_list_start)
        @:DEALLOCATE(cell_list_count)
        @:DEALLOCATE(cell_list_idx)

    end subroutine s_finalize_particle_lagrangian_solver

end module m_particles_EL
