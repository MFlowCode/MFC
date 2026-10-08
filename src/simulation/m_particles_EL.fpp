!>
!! @file m_particles_EL.fpp
!! @brief Contains module m_particles_EL

#:include 'macros.fpp'

!> @brief Euler-Lagrange solid particle solver with two-way coupling.
!!
!! Tracks non-deformable solid particles in compressible flow using Gaussian volume-averaging (Maeda & Colonius, J. Computational
!! Physics, 361, 2018). Supports multiple drag correlations, pressure gradient and added mass forces.
!! Derived from the m_bubbles_EL module. Kernel functions are in m_particles_EL_kernels.
module m_particles_EL

    use m_global_parameters     !< Definitions of the global parameters
    use m_mpi_proxy             !< Message passing interface (MPI) module proxy
    use m_particles_EL_kernels  !< Definitions of the kernel functions
    use m_euler_lagrange        !< Utilities shared with the Lagrangian bubbles
    use m_variables_conversion  !< State variables type conversion procedures
    use m_eos
    use m_compile_specific
    use m_boundary_common
    use m_helper_basic          !< Functions to compare floating point numbers
    use m_sim_helpers
    use m_helper
    use m_mpi_common
    use m_ibm
    use m_chemistry
    use ieee_arithmetic

    implicit none

    private
    public :: s_initialize_particles_EL_module, s_finalize_particle_lagrangian_solver, s_compute_particle_EL_dynamics, &
        & s_compute_particle_gradients, s_compute_particles_EL_source, s_update_lagrange_particles_tdv_rk, &
        & s_write_restart_lag_particles, s_write_lag_particle_evol, s_sync_particles_for_save, q_particles, alphaf_id

    real(wp)                             :: next_write_time
    integer, allocatable, dimension(:,:) :: lag_part_id  !< Global and local IDs
    $:GPU_DECLARE(create='[lag_part_id]')

    real(wp), allocatable, dimension(:) :: particle_mass  !< Particle Mass
    $:GPU_DECLARE(create='[particle_mass]')
    real(wp), allocatable, dimension(:) :: p_AM  !< Particle Added Mass
    $:GPU_DECLARE(create='[p_AM]')

    integer(seed_kind), allocatable, dimension(:) :: particle_seed  !< State of the particle's fluctuation-model RNG
    $:GPU_DECLARE(create='[particle_seed]')

    integer, allocatable, dimension(:) :: p_owner_rank  !< MPI rank that owns this particle
    $:GPU_DECLARE(create='[p_owner_rank]')

    ! Particle state arrays use dimensions (nParticles_glb, component, stage): component: 1=x, 2=y, 3=z for position/velocity stage:
    ! 1=committed state at current time level, 2=intermediate RK stage value

    real(wp), allocatable, dimension(:) :: particle_rad  !< Particle radius (constant: particles are rigid)
    $:GPU_DECLARE(create='[particle_rad]')

    ! (nPart, 1-> x or 2->y or 3 ->z, 1 -> actual or 2 -> temporal val)
    real(wp), allocatable, dimension(:,:,:) :: particle_pos      !< Particle's position
    real(wp), allocatable, dimension(:,:,:) :: particle_posPrev  !< Particle's previous position
    real(wp), allocatable, dimension(:,:,:) :: particle_vel      !< Particle's velocity
    real(wp), allocatable, dimension(:,:,:) :: particle_s        !< Particle's computational cell position in real format
    $:GPU_DECLARE(create='[particle_pos, particle_posPrev, particle_vel, particle_s]')
    ! (nPart, 1-> x or 2->y or 3 ->z, time-stage)
    real(wp), allocatable, dimension(:,:,:) :: particle_dposdt  !< Time derivative of the particle's position
    real(wp), allocatable, dimension(:,:,:) :: particle_dveldt  !< Time derivative of the particle's velocity
    $:GPU_DECLARE(create='[particle_dposdt, particle_dveldt]')

    integer, private :: lag_num_ts  !< Number of time stages in the time-stepping scheme
    $:GPU_DECLARE(create='[lag_num_ts]')

    !> Eulerian projection of particle data (volume fraction, momentum, sources)
    type(scalar_field), dimension(:), allocatable :: q_particles
    type(scalar_field), dimension(:), allocatable :: kahan_comp       !< Kahan compensation for q_particles accumulation
    integer                                       :: q_particles_idx  !< Size of the q vector field for particle cell (q)uantities

    !> Interpolated Eulerian field gradients at particle locations
    type(scalar_field), dimension(:), allocatable :: field_vars        !< For cell quantities (field gradients, etc.)
    type(scalar_field), dimension(:), allocatable :: rhs_old           !< For previous rhs values
    type(scalar_field), dimension(:), allocatable :: weights_x_interp  !< For precomputing weights
    type(scalar_field), dimension(:), allocatable :: weights_y_interp  !< For precomputing weights
    type(scalar_field), dimension(:), allocatable :: weights_z_interp  !< For precomputing weights
    integer                                       :: nWeights_interp
    type(scalar_field), dimension(:), allocatable :: weights_x_grad    !< For precomputing weights
    type(scalar_field), dimension(:), allocatable :: weights_y_grad    !< For precomputing weights
    type(scalar_field), dimension(:), allocatable :: weights_z_grad    !< For precomputing weights
    integer                                       :: nWeights_grad

    $:GPU_DECLARE(create='[q_particles, kahan_comp, q_particles_idx, field_vars, rhs_old]')
    $:GPU_DECLARE(create='[weights_x_interp, weights_y_interp, weights_z_interp, nWeights_interp]')
    $:GPU_DECLARE(create='[weights_x_grad, weights_y_grad, weights_z_grad, nWeights_grad]')

    ! Particle Source terms for fluid coupling
    real(wp), allocatable, dimension(:,:) :: f_p  !< force on each particle
    $:GPU_DECLARE(create='[f_p]')

    real(wp), allocatable, dimension(:,:) :: fqs_fluct  !< QS fluctuation force on each particle
    $:GPU_DECLARE(create='[fqs_fluct]')

    real(wp), allocatable, dimension(:) :: gSum  !< gaussian sum for each particle
    $:GPU_DECLARE(create='[gSum]')

    integer, allocatable, dimension(:) :: keep_particle
    $:GPU_DECLARE(create='[keep_particle]')

    real(wp)                              :: eps_overlap = 1.e-12
    real(wp), allocatable, dimension(:,:) :: fluid_vel_at_particle  !< fluid velocity at each particle
    real(wp), allocatable, dimension(:)   :: density_at_particle    !< density at each particle
    real(wp), allocatable, dimension(:)   :: pres_at_particle       !< fluid pressure at each particle
    $:GPU_DECLARE(create='[fluid_vel_at_particle, density_at_particle, pres_at_particle]')
    integer, allocatable, dimension(:)    :: force_status   !< First non-finite force term of each particle (0 if none)
    real(wp), allocatable, dimension(:,:) :: force_re_mach  !< Particle Reynolds and Mach numbers, reported with force_status
    $:GPU_DECLARE(create='[force_status, force_re_mach]')

contains

    !> Allocate the particle and projected-field arrays, read the particles (input file or restart), project them onto the grid, and
    !! precompute the interpolation and finite-difference weights.
    impure subroutine s_initialize_particles_EL_module(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer                                                    :: nParticles_glb, i, j, k, l, npts
        integer                                                    :: save_count
        real(wp)                                                   :: qtime
        integer                                                    :: ind_end_loc

        call s_get_lag_restart_point(save_count, qtime)

        ! The evolution file is written every t_save; a restart resumes at the first multiple of t_save not before qtime
        next_write_time = 0._wp
        if (save_count > 0 .and. t_save > 0._wp) next_write_time = t_save*ceiling(qtime/t_save)

        ! Setting number of time-stages for selected time-stepping scheme
        lag_num_ts = time_stepper

        ! Allocate space for the Eulerian fields needed to map the effect of the particles
        if (particle_params%solver_approach == 1) then
            ! One-way coupling
            q_particles_idx = 7  ! For tracking volume fraction, alpha_p u_p (x(2),y(3),z(4)), alpha_p u_p^2 (x(5),y(6),z(7))
        else if (particle_params%solver_approach == 2) then
            ! Two-way coupling
            ! For tracking volume fraction(1), alpha_p u_p (x(2),y(3),z(4)), alpha_p u_p^2 (x(5),y(6),z(7)), x-mom(8), y-mom(9),
            ! z-mom(10), and energy(11) sources
            q_particles_idx = 11
        else
            call s_mpi_abort('Please check the particle_params%solver_approach input')
        end if

        nWeights_interp = particle_params%interpolation_order + 1
        nWeights_grad = fd_order + 1

        call s_set_lag_comm_coords()

        $:GPU_UPDATE(device='[lag_num_ts, q_particles_idx]')

        @:ALLOCATE(q_particles(1:q_particles_idx))
        @:ALLOCATE(kahan_comp(1:q_particles_idx))
        do i = 1, q_particles_idx
            @:ALLOCATE(q_particles(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                       & idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_particles(i))
            @:ALLOCATE(kahan_comp(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                       & idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(kahan_comp(i))
        end do

        @:ALLOCATE(field_vars(1:nField_vars))
        do i = 1, nField_vars
            @:ALLOCATE(field_vars(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                       & idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(field_vars(i))
        end do

        ! Fluid density and momentum RHS for the added mass; fields allocated only when it is on
        @:ALLOCATE(rhs_old(1:eqn_idx%mom%end))
        if (particle_params%added_mass_force > 0) then
            do i = 1, eqn_idx%mom%end
                @:ALLOCATE(rhs_old(i)%sf(idwint(1)%beg:idwint(1)%end, idwint(2)%beg:idwint(2)%end, idwint(3)%beg:idwint(3)%end))
                @:ACC_SETUP_SFs(rhs_old(i))
            end do
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
        nParticles_glb = particle_params%nParticles_glb

        @:ALLOCATE(lag_part_id(1:nParticles_glb, 1:2))
        @:ALLOCATE(particle_mass(1:nParticles_glb))
        @:ALLOCATE(particle_seed(1:nParticles_glb))
        @:ALLOCATE(p_AM(1:nParticles_glb))
        @:ALLOCATE(p_owner_rank(1:nParticles_glb))
        @:ALLOCATE(particle_rad(1:nParticles_glb))
        @:ALLOCATE(particle_pos(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_posPrev(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_vel(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_s(1:nParticles_glb, 1:3, 1:2))
        @:ALLOCATE(particle_dposdt(1:nParticles_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(particle_dveldt(1:nParticles_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(f_p(1:nParticles_glb, 1:3))
        @:ALLOCATE(fqs_fluct(1:nParticles_glb, 1:3))
        @:ALLOCATE(gSum(1:nParticles_glb))

        @:ALLOCATE(keep_particle(1:nParticles_glb))

        @:ALLOCATE(fluid_vel_at_particle(1:nParticles_glb, 1:3))
        @:ALLOCATE(density_at_particle(1:nParticles_glb))
        @:ALLOCATE(pres_at_particle(1:nParticles_glb))
        @:ALLOCATE(force_status(1:nParticles_glb), force_re_mach(1:nParticles_glb, 1:2))

        if (adap_dt .and. f_is_default(adap_dt_tol)) adap_dt_tol = dflt_adap_dt_tol

        if (num_procs > 1) call s_initialize_solid_particles_mpi(lag_num_ts)

        call s_create_D_dir()

        ! Starting particles
        if (particle_params%write_void_evol) call s_open_void_evol
        if (particle_params%write_particles) then
            call s_open_lag_evol('particle', merge(17, 25, precision == 1), [character(len=11)::'currentTime','particleID', 'x', &
                                 & 'y', 'z', 'Vx', 'Vy', 'Vz', 'Fp_x', 'Fp_y', 'Fp_z', 'radius', 'vFx_at_p', 'vFy_at_p', &
                                 & 'vFz_at_p', 'rhoF_at_p'])
        end if

        call s_read_input_particles()

        call s_initialize_particle_kernels()

        if (particle_params%qs_fluct_force) then
            ind_end_loc = alphaup2z_id
        else if (particle_params%solver_approach == 2) then
            ind_end_loc = alphaupz_id
        else
            ind_end_loc = alphaf_id
        end if

        call s_smear_field_contributions(bc_type, alphaf_id, ind_end_loc, .true.)

        npts = (nWeights_interp - 1)/2
        call s_compute_barycentric_weights(npts)  ! For interpolation

        npts = (nWeights_grad - 1)/2
        call s_compute_fornberg_fd_weights(npts)  ! For finite differences

        if (particle_params%added_mass_force > 0) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
            do k = idwint(3)%beg, idwint(3)%end
                do j = idwint(2)%beg, idwint(2)%end
                    do i = idwint(1)%beg, idwint(1)%end
                        do l = 1, eqn_idx%mom%end
                            rhs_old(l)%sf(i, j, k) = 0._wp
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! Void fraction evolution at t = 0, now that the particles are smeared onto the grid
        if (save_count == 0 .and. particle_params%write_void_evol) call s_write_void_evol(qtime, q_particles(alphaf_id)%sf, &
            & particle_params%charwidth)

    end subroutine s_initialize_particles_EL_module

    !> Read the initial particles from particle_params%input_path (x, y, z, u, v, w, radius per line), or from the restart file when
    !! restarting, keeping those inside this rank's domain.
    impure subroutine s_read_input_particles()

        real(wp), allocatable, dimension(:,:) :: rows  !< x, y, z, u, v, w, radius per particle
        integer, allocatable, dimension(:)    :: ids
        integer                               :: k, particle_id, n_read, save_count
        real(wp)                              :: qtime

        call s_get_lag_restart_point(save_count, qtime)
        particle_id = 0

        if (save_count == 0) then
            if (proc_rank == 0) print *, 'Reading lagrange particles input file.'
            call s_read_lag_input(particle_params%input_path, 7, particle_params%nParticles_glb, rows, ids, particle_id, n_read)
            do k = 1, particle_id
                call s_add_particles(rows(k,:), k, ids(k))
                lag_part_id(k, 1) = ids(k)  ! global ID
                lag_part_id(k, 2) = k  ! local ID
            end do
            n_el_particles_loc = particle_id
        else
            if (proc_rank == 0) print *, 'Restarting lagrange particles at save_count: ', save_count
            call s_restart_particles(particle_id, save_count)
        end if

        call s_count_lag_glb(n_el_particles_loc, n_el_particles_glb, particle_params%input_path)
        if (proc_rank == 0) print '(A,I0)', ' Lagrangian particles in the domain: ', n_el_particles_glb

        $:GPU_UPDATE(device='[particles_lagrange, particle_params]')

        $:GPU_UPDATE(device='[lag_part_id, particle_mass, particle_seed, f_p, fqs_fluct, p_AM, p_owner_rank, particle_rad, &
                     & particle_pos, particle_posPrev, particle_vel, particle_s, particle_dposdt, particle_dveldt, n_el_particles_loc]')

        $:GPU_UPDATE(device='[dx, dy, dz, x_cb, x_cc, y_cb, y_cc, z_cb, z_cc]')

        ! Populate temporal variables
        call s_transfer_data_to_tmp_particles()

    end subroutine s_read_input_particles

    !> Store one particle from the input file: position, velocity, radius, mass from particle_pp%rho0ref_particle, and its
    !! fluctuation-model seed.
    impure subroutine s_add_particles(inputPart, part_id, glb_part_id)

        real(wp), dimension(7), intent(in) :: inputPart
        integer, intent(in)                :: part_id, glb_part_id
        real(wp)                           :: volparticle
        integer, dimension(3)              :: cell

        particle_rad(part_id) = inputPart(7)
        particle_pos(part_id,1:3,1) = inputPart(1:3)
        particle_posPrev(part_id,1:3,1) = particle_pos(part_id,1:3,1)
        if (.not. particle_params%stationary) then
            particle_vel(part_id,1:3,1) = inputPart(4:6)
        else
            particle_vel(part_id,1:3,1) = 0._wp
        end if

        ! Initialize Particle Sources
        f_p(part_id,1:3) = 0._wp
        fqs_fluct(part_id,1:3) = 0._wp
        p_AM(part_id) = 0._wp
        p_owner_rank(part_id) = proc_rank

        if (cyl_coord .and. p == 0) then
            particle_pos(part_id, 2, 1) = sqrt(particle_pos(part_id, 2, 1)**2._wp + particle_pos(part_id, 3, 1)**2._wp)
            ! Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            particle_pos(part_id, 3, 1) = atan2(inputPart(3), inputPart(2))
            particle_posPrev(part_id,1:3,1) = particle_pos(part_id,1:3,1)
        end if

        cell = fd_number - buff_size
        call s_locate_cell(particle_pos(part_id,1:3,1), cell, particle_s(part_id,1:3,1))

        ! Check if the particle is located in the ghost cell of a symmetric, or wall boundary
        if ((any(bc_x%beg == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, &
            & BC_NO_SLIP_WALL/)) .and. cell(1) < 0) .or. (any(bc_x%end == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, &
            & BC_NO_SLIP_WALL/)) .and. cell(1) > m) .or. (any(bc_y%beg == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, &
            & BC_NO_SLIP_WALL/)) .and. cell(2) < 0) .or. (any(bc_y%end == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, &
            & BC_NO_SLIP_WALL/)) .and. cell(2) > n)) then
            call s_mpi_abort("Lagrange particle is in the ghost cells of a symmetric or wall boundary.")
        end if

        if (p > 0) then
            if ((any(bc_z%beg == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, &
                & BC_NO_SLIP_WALL/)) .and. cell(3) < 0) .or. (any(bc_z%end == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, &
                & BC_NO_SLIP_WALL/)) .and. cell(3) > p)) then
                call s_mpi_abort("Lagrange particle is in the ghost cells of a symmetric or wall boundary.")
            end if
        end if

        ! Initial particle mass
        volparticle = 4._wp/3._wp*pi*particle_rad(part_id)**3  ! volume
        particle_mass(part_id) = volparticle*particle_pp%rho0ref_particle  ! mass
        if (particle_mass(part_id) <= 0._wp) then
            call s_mpi_abort("The initial particle mass is negative or zero. Check the particle file.")
        end if

        particle_seed(part_id) = glb_part_id  ! s_prng_splitmix32 gives unrelated streams for different seeds

    end subroutine s_add_particles

    !> Read this rank's particles from the restart file written at save_count.
    impure subroutine s_restart_particles(part_id, save_count)

        integer, intent(inout)                :: part_id, save_count
        real(wp), allocatable, dimension(:,:) :: io_data
        integer, dimension(3)                 :: cell
        integer                               :: i

        call s_read_lag_restart('particles', save_count, io_data, part_id)
        if (.not. allocated(io_data)) return

        n_el_particles_loc = part_id
        do i = 1, part_id
            lag_part_id(i, 1) = int(io_data(i, 1))
            particle_pos(i,1:3,1) = io_data(i,2:4)
            particle_posPrev(i,1:3,1) = io_data(i,5:7)
            particle_vel(i,1:3,1) = io_data(i,8:10)
            particle_rad(i) = io_data(i, 11)
            ! rhs_old (fluid acceleration for the added mass) is not saved: the first stage after a restart uses zero
            fqs_fluct(i,1:3) = io_data(i,16:18)
            particle_mass(i) = io_data(i, 19)
            particle_seed(i) = ior(ishft(int(nint(io_data(i, 21)), seed_kind), 16), int(nint(io_data(i, 20)), seed_kind))
            ! Restart files without a saved seed hold 0: seed from the ID as on a fresh start
            if (particle_seed(i) == 0) particle_seed(i) = lag_part_id(i, 1)
            cell = -buff_size
            call s_locate_cell(particle_pos(i,1:3,1), cell, particle_s(i,1:3,1))
        end do
        deallocate (io_data)

    end subroutine s_restart_particles

    !> Interpolate the fluid state to each particle, compute the drag, pressure-gradient, added-mass and fluctuating forces, and
    !! store the particle RHS for this RK stage. With two-way coupling, also project the forces onto the grid.
    subroutine s_compute_particle_EL_dynamics(q_prim_vf, bc_type, stage)

        type(scalar_field), dimension(sys_size), intent(in)        :: q_prim_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer, intent(in)                                        :: stage
        integer, dimension(3)                                      :: cell
        real(wp)                                                   :: myMass, myR, myRe, myGamma, rmass_add, myFluidRho, myPres
        real(wp)                                                   :: qv, pi_inf, gamma, vel_sum, c, pres, rho
        real(wp), dimension(2)                                     :: Re

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: vel        !< Cell-avg. velocity
            real(wp), dimension(3) :: alpha      !< Cell-avg. volume fraction
            real(wp), dimension(3) :: alpha_rho  !< Cell-avg. partial density
        #:else
            real(wp), dimension(num_vels)   :: vel        !< Cell-avg. velocity
            real(wp), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
            real(wp), dimension(num_fluids) :: alpha_rho  !< Cell-avg. partial density
        #:endif

        real(wp), dimension(3) :: myVel, myPos, force_vec, s_cell, my_fqs_fluct, new_fqs_fluct, myFluidVel
        integer(seed_kind)     :: mySeed, new_seed
        integer                :: k, l, i, i_c, j_c, k_c
        integer                :: my_status, max_status
        real(wp)               :: my_re_p, my_mach_p

        call nvtxStartRange("LAGRANGE-PARTICLE-DYNAMICS")

        ! Compute Fluid-Particle Forces (drag/pressure/added mass) and convert to particle acceleration
        $:GPU_PARALLEL_LOOP(private='[i, k, l, cell, s_cell, myMass, myR, myPos, myVel, mySeed, my_fqs_fluct, new_fqs_fluct, &
                            & force_vec, rmass_add, new_seed, myFluidVel, myFluidRho, myPres, qv, pi_inf, gamma, myGamma, &
                            & vel_sum, vel, alpha, alpha_rho, Re, myRe, pres, rho, i_c, j_c, k_c, c, my_status, my_re_p, &
                            & my_mach_p]', copyin='[stage]')
        do k = 1, n_el_particles_loc
            f_p(k,:) = 0._wp
            p_owner_rank(k) = proc_rank

            s_cell = particle_s(k,1:3,2)
            cell = int(s_cell(:))
            do i = 1, num_dims
                if (s_cell(i) < 0._wp) cell(i) = cell(i) - 1
            end do

            ! Current particle state
            myMass = particle_mass(k)
            myR = particle_rad(k)
            myPos = particle_pos(k,:,2)
            myVel = particle_vel(k,:,2)

            mySeed = particle_seed(k)
            my_fqs_fluct = fqs_fluct(k,:)

            particle_dposdt(k,:,stage) = 0._wp
            particle_dveldt(k,:,stage) = 0._wp
            density_at_particle(k) = 0._wp

            call s_interp_fluid_properties(myPos, cell, q_prim_vf, weights_x_interp, weights_y_interp, weights_z_interp, &
                                           & myFluidVel, myFluidRho, myPres)

            fluid_vel_at_particle(k,:) = myFluidVel
            density_at_particle(k) = myFluidRho
            pres_at_particle(k) = myPres

            ! Mixture properties of the host cell (scalar indices: device routine with a seq loop)
            i_c = cell(1); j_c = cell(2); k_c = cell(3)
            call s_compute_cell_state(q_prim_vf, pres, rho, gamma, pi_inf, Re, alpha, alpha_rho, vel, vel_sum, qv, i_c, j_c, k_c)

            ! Compute mixture sound speed
            call s_compute_speed_of_sound(myPres, myFluidRho, gamma, pi_inf, alpha, c)

            myGamma = (1._wp/gamma) + 1._wp

            call s_get_drag_viscosity(q_prim_vf, myPres, myFluidrho, pi_inf, alpha, Re, myRe, cell(1), cell(2), cell(3))

            call s_get_particle_force(myPos, myR, myVel, myRe, myGamma, mySeed, my_fqs_fluct, stage == 1, cell, q_particles, &
                                      & field_vars, rhs_old, weights_x_interp, weights_y_interp, weights_z_interp, force_vec, &
                                      & rmass_add, new_seed, new_fqs_fluct, myFluidVel, myFluidRho, c, my_status, my_re_p, &
                                      & my_mach_p)
            force_status(k) = my_status
            force_re_mach(k, 1) = my_re_p
            force_re_mach(k, 2) = my_mach_p

            p_AM(k) = rMass_add
            f_p(k,:) = f_p(k,:) + force_vec(:)

            if (particle_params%qs_fluct_force) then
                particle_seed(k) = new_seed
                fqs_fluct(k,:) = new_fqs_fluct
            end if

            myMass = particle_mass(k) + p_AM(k)
            myVel = particle_vel(k,:,2)
            do l = 1, num_dims
                particle_dposdt(k, l, stage) = myVel(l)
                particle_dveldt(k, l, stage) = f_p(k, l)/myMass
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Abort on any non-finite particle force instead of letting it propagate
        max_status = 0
        $:GPU_PARALLEL_LOOP(private='[k]', reduction='[[max_status]]', reductionOp='[MAX]', copy='[max_status]')
        do k = 1, n_el_particles_loc
            max_status = max(max_status, force_status(k))
        end do
        $:END_GPU_PARALLEL_LOOP()
        if (max_status > 0) call s_abort_nonfinite_force()

        if (particle_params%solver_approach == 2) then
            call s_smear_field_contributions(bc_type, Smx_id, SE_id, .false.)
        end if

        call nvtxEndRange

    end subroutine s_compute_particle_EL_dynamics

    !> Report every particle whose force is not finite (position, velocities, fluid state, Re, Ma and the force term that failed
    !! first) and abort the run.
    impure subroutine s_abort_nonfinite_force()

        integer                      :: k, n_bad
        integer, parameter           :: max_reported = 10
        character(len=16), parameter :: qs_names(3) = [character(len=16)::'Gidaspow','Parmar', 'Osnes']
        character(len=16), parameter :: qs_routines(3) = [character(len=16)::'QS_Gidaspow','QS_Parmar', 'QS_Osnes']

        $:GPU_UPDATE(host='[force_status, force_re_mach, lag_part_id, particle_pos, particle_vel, fluid_vel_at_particle, &
                     & density_at_particle, pres_at_particle, f_p]')

        n_bad = count(force_status(1:n_el_particles_loc) > 0)
        print '(A,I0,A,I0,A)', ' Non-finite particle force on rank ', proc_rank, ' (', n_bad, ' particles):'
        do k = 1, n_el_particles_loc
            if (force_status(k) == 0) cycle
            print '(A,I0,2A)', '   particle ', lag_part_id(k, 1), ': non-finite ', trim(force_term_names(force_status(k)))
            select case (force_status(k))
            case (1)
                print '(5A)', '     where to look     ', trim(qs_names(particle_params%qs_force)), ' correlation (', &
                    & trim(qs_routines(particle_params%qs_force)), &
                    & ' in m_particles_EL_kernels.fpp): fluid density, slip ' &
                    & // 'velocity, drag viscosity (mu_ref, suth), Re_p, Ma_p'
            case (2)
                print '(A)', &
                    & '     where to look     pressure gradient at the particle (s_gradient_field, from the ' &
                    & // 'reconstructed face states), fluid pressure'
            case (3)
                print '(A)', &
                    & '     where to look     added mass in s_get_particle_force: fluid density and its gradient, fluid ' &
                    & // 'acceleration from rhs_old, Ma_p'
            case (4)
                print '(A)', &
                    & '     where to look     s_compute_qs_fluctuations: granular temperature, particle volume ' &
                    & // 'fraction, fluctuating force of the previous step'
            end select
            print '(A,3ES14.6)', '     particle force    ', f_p(k,1:3)
            print '(A,3ES14.6)', '     position          ', particle_pos(k,1:3,2)
            print '(A,3ES14.6)', '     particle velocity ', particle_vel(k,1:3,2)
            print '(A,3ES14.6)', '     fluid velocity    ', fluid_vel_at_particle(k,1:3)
            print '(A,2ES14.6)', '     fluid rho, p      ', density_at_particle(k), pres_at_particle(k)
            print '(A,2ES14.6)', '     Re_p, Ma_p        ', force_re_mach(k,1:2)
            if (count(force_status(1:k) > 0) == max_reported) exit
        end do

        call s_mpi_abort('Non-finite particle force; see the report above. Check dt, the fluid state near the particle, ' &
                         & // 'and the particle input.')

    end subroutine s_abort_nonfinite_force

    !> Fluid viscosity for the drag correlations: from the mixture Reynolds number when viscous, from the mixture model with
    !! chemistry, otherwise mu_ref with an optional Sutherland correction.
    subroutine s_get_drag_viscosity(q_prim_vf, fluid_pres, fluid_rho, pi_inf_mix, alpha, Re_mix, mu, i, j, k)

        $:GPU_ROUTINE(function_name='s_get_drag_viscosity',parallelism='[seq]', cray_inline=True)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), intent(in)                                :: fluid_pres, fluid_rho, pi_inf_mix
        real(wp), dimension(num_fluids), intent(in)         :: alpha
        real(wp), dimension(2), intent(in)                  :: Re_mix
        real(wp), intent(out)                               :: mu
        integer, intent(in)                                 :: i, j, k
        real(wp), dimension(num_species)                    :: Ys
        real(wp), parameter                                 :: tref = 273._wp  !< Sutherland reference temperature
        real(wp)                                            :: mu_f
        real(wp)                                            :: fluid_temp, mix_mol_weight, cv_mix
        integer                                             :: l, d

        if (chemistry) then
            do d = 1, num_species
                Ys(d) = q_prim_vf(eqn_idx%species%beg + d - 1)%sf(i, j, k)
            end do
            call get_mixture_molecular_weight(Ys, mix_mol_weight)
            fluid_temp = fluid_pres*mix_mol_weight/(gas_constant*fluid_rho)
            call get_mixture_viscosity_mixavg(fluid_temp, Ys, mu)
        else
            if (viscous) then
                mu = 1._wp/Re_mix(1)
            else
                mu = 0._wp
                cv_mix = 0._wp
                do l = 1, num_fluids
                    cv_mix = cv_mix + alpha(l)*cvs(l)/gammas(l)  ! R = cv*(gamma - 1), with gammas = 1/(gamma - 1)
                end do
                fluid_temp = (fluid_pres + pi_inf_mix)/(fluid_rho*cv_mix)
                do l = 1, num_fluids
                    ! Constant viscosity unless a Sutherland constant is given
                    mu_f = particle_params%mu_ref(l)
                    if (particle_params%suth(l) > 0._wp) then
                        mu_f = mu_f*sqrt(fluid_temp/tref)*(1._wp + particle_params%suth(l)/tref)/(1._wp + particle_params%suth(l) &
                                         & /fluid_temp)
                    end if
                    mu = mu + alpha(l)*mu_f
                end do
            end if
        end if

    end subroutine s_get_drag_viscosity

    !> Project the particle quantities ind_start:ind_end onto the grid with the Gaussian kernel (Maeda and Colonius, 2018), exchange
    !! them across rank boundaries, and turn the projected particle volume fraction into the fluid volume fraction. recompute_gSum
    !! relocates the particles and recomputes the kernel normalization.
    subroutine s_smear_field_contributions(bc_type, ind_start, ind_end, recompute_gSum)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer, intent(in)                                        :: ind_start, ind_end
        logical, intent(in)                                        :: recompute_gSum
        integer                                                    :: i, j, k, l, nVar
        real(wp)                                                   :: myR, func_sum
        real(wp), dimension(3)                                     :: myVel, myPos, s_cell, myForce
        integer, dimension(3)                                      :: cell
        integer, dimension(:), allocatable                         :: vars_send

        nVar = ind_end - ind_start + 1
        allocate (vars_send(nVar))

        do i = 1, nVar
            vars_send(i) = ind_start + (i - 1)
        end do

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
        do i = ind_start, ind_end
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        if (i <= q_particles_idx) then
                            q_particles(i)%sf(j, k, l) = 0._wp
                            kahan_comp(i)%sf(j, k, l) = 0._wp
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[i, k, cell, s_cell, myR, myPos, myVel, myForce, func_sum]')
        do k = 1, n_el_particles_loc
            myR = particle_rad(k)
            myPos = particle_pos(k,1:3,2)
            myVel = particle_vel(k,1:3,2)
            ! The gas receives the particle momentum change m du/dt; f_p excludes the -m_a du/dt solved implicitly
            myForce = f_p(k,:)*particle_mass(k)/(particle_mass(k) + p_AM(k))

            if (recompute_gSum) then
                cell = fd_number - buff_size
                call s_locate_cell(particle_pos(k,1:3,2), cell, particle_s(k,1:3,2))

                ! Compute the total gaussian contribution for each particle for normalization
                call s_compute_gaussian_contribution(myPos, cell, func_sum)
                gSum(k) = func_sum
            else
                s_cell = particle_s(k,1:3,2)
                cell = int(s_cell(:))
                do i = 1, num_dims
                    if (s_cell(i) < 0._wp) cell(i) = cell(i) - 1
                end do
                func_sum = gSum(k)
            end if

            call s_gaussian_atomic(myR, myVel, myPos, myForce, func_sum, cell, q_particles, ind_start, ind_end)
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxStartRange("PARTICLES-LAGRANGE-BETA-COMM")
        call s_populate_beta_buffers(q_particles, kahan_comp, bc_type, nVar, vars_send)
        call nvtxEndRange

        if (alphaf_id >= ind_start .and. alphaf_id <= ind_end) then
            ! Store 1-q_particles(1)
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        q_particles(alphaf_id)%sf(j, k, l) = 1._wp - q_particles(alphaf_id)%sf(j, k, l)
                        ! Limiting void fraction given max value
                        q_particles(alphaf_id)%sf(j, k, l) = max(q_particles(alphaf_id)%sf(j, k, l), &
                                    & 1._wp - particle_params%valmaxvoid)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_smear_field_contributions

    !> Add the two-way coupling source terms to the fluid RHS: volume-fraction terms and the projected particle forces (Maeda and
    !! Colonius, 2018).
    subroutine s_compute_particles_EL_source(q_cons_vf, q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: i, j, k, l
        real(wp)                                               :: dalphapdt, alpha_f, udot_gradalpha

        ! Spatial derivative of the fluid volume fraction and eulerian particle momentum fields.

        do l = 1, num_dims
            call s_gradient_dir_fornberg(q_particles(alphaf_id)%sf, field_vars(dalphafx_id + l - 1)%sf, l)
            call s_gradient_dir_fornberg(q_particles(alphaupx_id + l - 1)%sf, field_vars(dalphap_upx_id + l - 1)%sf, l)
        end do

        ! Pressure terms -(alpha_p/alpha_f) dp/dx_l and -(alpha_p/alpha_f) d(p u_l)/dx_l, the form the bubble solver uses: with
        ! the full particle force deposited, a cloud at rest in uniform pressure stays at rest
        do l = 1, num_dims
            call s_gradient_dir_fornberg(q_prim_vf(eqn_idx%E)%sf, field_vars(dsrc_tmp_id)%sf, l)
            call s_add_pressure_source(rhs_vf(eqn_idx%mom%beg + l - 1)%sf)

            $:GPU_PARALLEL_LOOP(private='[i, j, k]', collapse=3)
            do k = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(2)%beg, idwbuff(2)%end
                    do i = idwbuff(1)%beg, idwbuff(1)%end
                        field_vars(src_tmp_id)%sf(i, j, k) = q_prim_vf(eqn_idx%E)%sf(i, j, &
                                   & k)*q_prim_vf(eqn_idx%mom%beg + l - 1)%sf(i, j, k)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_gradient_dir_fornberg(field_vars(src_tmp_id)%sf, field_vars(dsrc_tmp_id)%sf, l)
            call s_add_pressure_source(rhs_vf(eqn_idx%E)%sf)
        end do

        ! Apply particle sources to the Eulerian RHS
        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, alpha_f, dalphapdt, udot_gradalpha]', collapse=3)
        do k = idwint(3)%beg, idwint(3)%end
            do j = idwint(2)%beg, idwint(2)%end
                do i = idwint(1)%beg, idwint(1)%end
                    if (q_particles(alphaf_id)%sf(i, j, k) > (1._wp - particle_params%valmaxvoid)) then
                        alpha_f = q_particles(alphaf_id)%sf(i, j, k)

                        dalphapdt = 0._wp
                        udot_gradalpha = 0._wp
                        do l = 1, num_dims
                            dalphapdt = dalphapdt + field_vars(dalphap_upx_id + l - 1)%sf(i, j, k)
                            udot_gradalpha = udot_gradalpha + q_prim_vf(eqn_idx%mom%beg + l - 1)%sf(i, j, &
                                & k)*field_vars(dalphafx_id + l - 1)%sf(i, j, k)
                        end do
                        dalphapdt = -dalphapdt
                        ! Add any contribution to dalphapdt from particles growing or shrinking

                        ! Step 1: Source terms for volume fraction corrections
                        ! cons_var/alpha_f * (dalpha_p/dt - u dot grad(alpha_f))
                        do l = 1, eqn_idx%E
                            rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + (q_cons_vf(l)%sf(i, j, &
                                   & k)/alpha_f)*(dalphapdt - udot_gradalpha)
                        end do

                        ! Step 2: Interphase momentum and energy exchange (minus the particle's momentum change and work)
                        do l = 1, num_dims
                            rhs_vf(eqn_idx%mom%beg + l - 1)%sf(i, j, k) = rhs_vf(eqn_idx%mom%beg + l - 1)%sf(i, j, &
                                   & k) + q_particles(Smx_id + l - 1)%sf(i, j, k)/alpha_f
                        end do
                        rhs_vf(eqn_idx%E)%sf(i, j, k) = rhs_vf(eqn_idx%E)%sf(i, j, k) + q_particles(SE_id)%sf(i, j, k)/alpha_f

                        if (particle_params%added_mass_force > 0) then
                            do l = 1, eqn_idx%mom%end
                                rhs_old(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k)
                            end do
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_particles_EL_source

    !> rhs -= (alpha_p/alpha_f)*field_vars(dsrc_tmp_id) in the cells where the particle sources apply.
    subroutine s_add_pressure_source(rhs)

        real(stp), dimension(idwint(1)%beg:,idwint(2)%beg:,idwint(3)%beg:), intent(inout) :: rhs
        integer                                                                           :: i, j, k
        real(wp)                                                                          :: alpha_f

        $:GPU_PARALLEL_LOOP(private='[i, j, k, alpha_f]', collapse=3)
        do k = idwint(3)%beg, idwint(3)%end
            do j = idwint(2)%beg, idwint(2)%end
                do i = idwint(1)%beg, idwint(1)%end
                    alpha_f = q_particles(alphaf_id)%sf(i, j, k)
                    if (alpha_f > (1._wp - particle_params%valmaxvoid)) then
                        rhs(i, j, k) = rhs(i, j, k) - real((1._wp - alpha_f)/alpha_f*field_vars(dsrc_tmp_id)%sf(i, j, k), kind=stp)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_add_pressure_source

    !> Advance particle positions and velocities by one stage of the TVD Runge-Kutta scheme, then apply boundary conditions and rank
    !! handover and write the requested output.
    impure subroutine s_update_lagrange_particles_tdv_rk(bc_type, stage)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer, intent(in)                                        :: stage
        integer                                                    :: k

        if (particle_params%write_particles .and. mytime == 0._wp .and. stage == 1) then
            call s_write_lag_particle_evol(mytime)
            next_write_time = next_write_time + t_save
        end if

        if (particle_params%write_particles .and. mytime >= next_write_time .and. stage == lag_num_ts) then
            call s_write_lag_particle_evol(mytime)
            next_write_time = next_write_time + t_save
        end if

        if (time_stepper == 1) then  ! 1st order TVD RK
            $:GPU_PARALLEL_LOOP(private='[k]')
            do k = 1, n_el_particles_loc
                ! u{1} = u{n} +  dt * RHS{n}
                if (.not. particle_params%stationary) then
                    particle_posPrev(k,1:3,1) = particle_pos(k,1:3,1)
                    particle_pos(k,1:3,1) = particle_pos(k,1:3,1) + dt*particle_dposdt(k,1:3,1)
                    particle_vel(k,1:3,1) = particle_vel(k,1:3,1) + dt*particle_dveldt(k,1:3,1)
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()

            call s_transfer_data_to_tmp_particles()

            if (.not. particle_params%stationary) call s_enforce_EL_particles_boundary_conditions(stage, bc_type)
            if (particle_params%write_void_evol) call s_write_void_evol(mytime, q_particles(alphaf_id)%sf, &
                & particle_params%charwidth)
        else if (time_stepper == 2) then  ! 2nd order TVD RK
            if (stage == 1) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    ! u{1} = u{n} +  dt * RHS{n}
                    if (.not. particle_params%stationary) then
                        particle_posPrev(k,1:3,2) = particle_pos(k,1:3,1)
                        particle_pos(k,1:3,2) = particle_pos(k,1:3,1) + dt*particle_dposdt(k,1:3,1)
                        particle_vel(k,1:3,2) = particle_vel(k,1:3,1) + dt*particle_dveldt(k,1:3,1)
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (.not. particle_params%stationary) call s_enforce_EL_particles_boundary_conditions(stage, bc_type)
            else if (stage == 2) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    ! u{1} = u{n} + (1/2) * dt * (RHS{n} + RHS{1})
                    if (.not. particle_params%stationary) then
                        particle_posPrev(k,1:3,1) = particle_pos(k,1:3,2)
                        particle_pos(k,1:3,1) = particle_pos(k,1:3,1) + dt*(particle_dposdt(k,1:3,1) + particle_dposdt(k,1:3, &
                                     & 2))/2._wp
                        particle_vel(k,1:3,1) = particle_vel(k,1:3,1) + dt*(particle_dveldt(k,1:3,1) + particle_dveldt(k,1:3, &
                                     & 2))/2._wp
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_transfer_data_to_tmp_particles()

                if (.not. particle_params%stationary) call s_enforce_EL_particles_boundary_conditions(stage, bc_type)
                if (particle_params%write_void_evol) call s_write_void_evol(mytime, q_particles(alphaf_id)%sf, &
                    & particle_params%charwidth)
            end if
        else if (time_stepper == 3) then  ! 3rd order TVD RK
            if (stage == 1) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    ! u{1} = u{n} +  dt * RHS{n}
                    if (.not. particle_params%stationary) then
                        particle_posPrev(k,1:3,2) = particle_pos(k,1:3,1)
                        particle_pos(k,1:3,2) = particle_pos(k,1:3,1) + dt*particle_dposdt(k,1:3,1)
                        particle_vel(k,1:3,2) = particle_vel(k,1:3,1) + dt*particle_dveldt(k,1:3,1)
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (.not. particle_params%stationary) call s_enforce_EL_particles_boundary_conditions(stage, bc_type)
            else if (stage == 2) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    ! u{2} = u{n} + (1/4) * dt * [RHS{n} + RHS{1}]
                    if (.not. particle_params%stationary) then
                        particle_posPrev(k,1:3,2) = particle_pos(k,1:3,2)
                        particle_pos(k,1:3,2) = particle_pos(k,1:3,1) + dt*(particle_dposdt(k,1:3,1) + particle_dposdt(k,1:3, &
                                     & 2))/4._wp
                        particle_vel(k,1:3,2) = particle_vel(k,1:3,1) + dt*(particle_dveldt(k,1:3,1) + particle_dveldt(k,1:3, &
                                     & 2))/4._wp
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (.not. particle_params%stationary) call s_enforce_EL_particles_boundary_conditions(stage, bc_type)
            else if (stage == 3) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, n_el_particles_loc
                    ! u{n+1} = u{n} + (2/3) * dt * [(1/4)* RHS{n} + (1/4)* RHS{1} + RHS{2}]
                    if (.not. particle_params%stationary) then
                        particle_posPrev(k,1:3,1) = particle_pos(k,1:3,2)
                        particle_pos(k,1:3,1) = particle_pos(k,1:3,1) + (2._wp/3._wp)*dt*(particle_dposdt(k,1:3, &
                                     & 1)/4._wp + particle_dposdt(k,1:3,2)/4._wp + particle_dposdt(k,1:3,3))
                        particle_vel(k,1:3,1) = particle_vel(k,1:3,1) + (2._wp/3._wp)*dt*(particle_dveldt(k,1:3, &
                                     & 1)/4._wp + particle_dveldt(k,1:3,2)/4._wp + particle_dveldt(k,1:3,3))
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_transfer_data_to_tmp_particles()

                if (.not. particle_params%stationary) call s_enforce_EL_particles_boundary_conditions(stage, bc_type)
                if (particle_params%write_void_evol) call s_write_void_evol(mytime, q_particles(alphaf_id)%sf, &
                    & particle_params%charwidth)
            end if
        end if

    end subroutine s_update_lagrange_particles_tdv_rk

    !> Enforce boundary conditions on Lagrangian particles. Phases: (1) GPU->host transfer, (2) MPI particle exchange with
    !! neighbors, (3) host->GPU transfer, (4) per-particle BC (clamp at walls, remove outside the domain), (5) compaction to remove
    !! deleted particles, (6) re-smear onto Eulerian grid.
    impure subroutine s_enforce_EL_particles_boundary_conditions(nstage, bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer, intent(in)                                        :: nstage
        integer                                                    :: k
        integer                                                    :: newParts
        integer                                                    :: ind_end_loc

        call nvtxStartRange("LAG-BC")
        call nvtxStartRange("LAG-BC-DEV2HOST")
        $:GPU_UPDATE(host='[p_owner_rank, particle_mass, particle_seed, f_p, fqs_fluct, lag_part_id, particle_rad, particle_pos, &
                     & particle_posPrev, particle_vel, particle_s, particle_dposdt, particle_dveldt, keep_particle, n_el_particles_loc]')
        call nvtxEndRange

        ! Handle MPI transfer of particles going to another processor's local domain
        if (num_procs > 1) then
            call nvtxStartRange("LAG-BC-TRANSFER-LIST")
            call s_add_particles_to_transfer_list(n_el_particles_loc, particle_pos(:,:,2), particle_posPrev(:,:,2))
            call nvtxEndRange

            call nvtxStartRange("LAG-BC-SENDRECV")
            call s_mpi_sendrecv_solid_particles(p_owner_rank, particle_mass, particle_seed, f_p, fqs_fluct, lag_part_id, &
                                                & particle_rad, particle_pos, particle_posPrev, particle_vel, particle_s, &
                                                & particle_dposdt, particle_dveldt, lag_num_ts, n_el_particles_loc, 2)
            call nvtxEndRange
        end if

        call nvtxStartRange("LAG-BC-HOST2DEV")
        $:GPU_UPDATE(device='[p_owner_rank, particle_mass, particle_seed, f_p, fqs_fluct, lag_part_id, particle_rad, &
                     & particle_pos, particle_posPrev, particle_vel, particle_s, particle_dposdt, particle_dveldt, n_el_particles_loc]')
        call nvtxEndRange

        $:GPU_PARALLEL_LOOP(private='[k]',copyin='[nstage]')
        do k = 1, n_el_particles_loc
            keep_particle(k) = 1

            ! Relocate particles at solid boundaries and delete particles that leave buffer regions
            if (any(bc_x%beg == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. particle_pos(k, 1, &
                & 2) < x_cb(-1) + eps_overlap) then
                particle_pos(k, 1, 2) = x_cb(-1) + eps_overlap
                if (nstage == lag_num_ts) then
                    particle_pos(k, 1, 1) = particle_pos(k, 1, 2)
                end if
            else if (any(bc_x%end == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. particle_pos(k, 1, &
                     & 2) > x_cb(m) - eps_overlap) then
                particle_pos(k, 1, 2) = x_cb(m) - eps_overlap
                if (nstage == lag_num_ts) then
                    particle_pos(k, 1, 1) = particle_pos(k, 1, 2)
                end if
            else if (particle_pos(k, 1, 2) >= x_cb(m)) then
                keep_particle(k) = 0
            else if (particle_pos(k, 1, 2) < x_cb(-1)) then
                keep_particle(k) = 0
            end if

            if (any(bc_y%beg == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. particle_pos(k, 2, &
                & 2) < y_cb(-1) + eps_overlap) then
                particle_pos(k, 2, 2) = y_cb(-1) + eps_overlap
                if (nstage == lag_num_ts) then
                    particle_pos(k, 2, 1) = particle_pos(k, 2, 2)
                end if
            else if (any(bc_y%end == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. particle_pos(k, 2, &
                     & 2) > y_cb(n) - eps_overlap) then
                particle_pos(k, 2, 2) = y_cb(n) - eps_overlap
                if (nstage == lag_num_ts) then
                    particle_pos(k, 2, 1) = particle_pos(k, 2, 2)
                end if
            else if (particle_pos(k, 2, 2) >= y_cb(n)) then
                keep_particle(k) = 0
            else if (particle_pos(k, 2, 2) < y_cb(-1)) then
                keep_particle(k) = 0
            end if

            if (p > 0) then
                if (any(bc_z%beg == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. particle_pos(k, 3, &
                    & 2) < z_cb(-1) + eps_overlap) then
                    particle_pos(k, 3, 2) = z_cb(-1) + eps_overlap
                    if (nstage == lag_num_ts) then
                        particle_pos(k, 3, 1) = particle_pos(k, 3, 2)
                    end if
                else if (any(bc_z%end == (/BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. particle_pos(k, 3, &
                         & 2) > z_cb(p) - eps_overlap) then
                    particle_pos(k, 3, 2) = z_cb(p) - eps_overlap
                    if (nstage == lag_num_ts) then
                        particle_pos(k, 3, 1) = particle_pos(k, 3, 2)
                    end if
                else if (particle_pos(k, 3, 2) >= z_cb(p)) then
                    keep_particle(k) = 0
                else if (particle_pos(k, 3, 2) < z_cb(-1)) then
                    keep_particle(k) = 0
                end if
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (n_el_particles_loc > 0) then
            call nvtxStartRange("LAG-BC")
            call nvtxStartRange("LAG-BC-DEV2HOST")
            $:GPU_UPDATE(host='[p_owner_rank, particle_mass, particle_seed, f_p, fqs_fluct, lag_part_id, particle_rad, &
                         & particle_pos, particle_posPrev, particle_vel, particle_s, particle_dposdt, particle_dveldt, &
                         & keep_particle, n_el_particles_loc]')
            call nvtxEndRange

            newParts = 0
            do k = 1, n_el_particles_loc
                if (keep_particle(k) == 1) then
                    newParts = newParts + 1
                    if (newParts /= k) then
                        call s_copy_lag_particle(newParts, k)
                    end if
                end if
            end do

            n_el_particles_loc = newParts

            call nvtxStartRange("LAG-BC-HOST2DEV")
            $:GPU_UPDATE(device='[p_owner_rank, particle_mass, particle_seed, f_p, fqs_fluct, lag_part_id, particle_rad, &
                         & particle_pos, particle_posPrev, particle_vel, particle_s, particle_dposdt, particle_dveldt, n_el_particles_loc]')
            call nvtxEndRange
        end if

        if (particle_params%qs_fluct_force) then
            ind_end_loc = alphaup2z_id
        else if (particle_params%solver_approach == 2) then
            ind_end_loc = alphaupz_id
        else
            ind_end_loc = alphaf_id
        end if

        call s_smear_field_contributions(bc_type, alphaf_id, ind_end_loc, .true.)

        call nvtxEndRange  ! LAG-BC

    end subroutine s_enforce_EL_particles_boundary_conditions

    !> Copy the committed particle state (time level 1) into the stage state (time level 2).
    impure subroutine s_transfer_data_to_tmp_particles()

        integer :: k

        $:GPU_PARALLEL_LOOP(private='[k]')
        do k = 1, n_el_particles_loc
            particle_pos(k,1:3,2) = particle_pos(k,1:3,1)
            particle_posPrev(k,1:3,2) = particle_posPrev(k,1:3,1)
            particle_vel(k,1:3,2) = particle_vel(k,1:3,1)
            particle_s(k,1:3,2) = particle_s(k,1:3,1)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_transfer_data_to_tmp_particles

    !> Store the particle field gradients along dir (pressure, velocity, and density with added mass) from this sweep's
    !! reconstructed face states. Called inside the RHS direction loop, so the face states need no extra copy.
    subroutine s_compute_particle_gradients(qL, qR, dir)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: qL, qR
        integer, intent(in)                                                                 :: dir
        integer                                                                             :: i

        if (.not. (particle_params%pressure_gradient_force .or. particle_params%added_mass_force > 0)) return

        call s_gradient_field(qL, qR, field_vars(dPx_id + dir - 1)%sf, dir, eqn_idx%E, eqn_idx%E)
        do i = 1, num_dims
            call s_gradient_field(qL, qR, field_vars(duidxj_id(i, dir))%sf, dir, eqn_idx%mom%beg + i - 1, eqn_idx%mom%beg + i - 1)
        end do
        if (particle_params%added_mass_force > 0) then
            ! Mixture density: the sum of the partial densities
            call s_gradient_field(qL, qR, field_vars(drhox_id + dir - 1)%sf, dir, eqn_idx%cont%beg, eqn_idx%cont%end)
        end if

    end subroutine s_compute_particle_gradients

    !> Cell-centered derivative along dir of the sum of components fbeg:fend: the difference of the right and left reconstructed
    !! face states over the cell width. The face states are stored in (x, y, z) order for every sweep direction, so dq shares their
    !! (i, j, k) indexing.
    subroutine s_gradient_field(vL_field, vR_field, dq, dir, fbeg, fend)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:), intent(out)  :: dq
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: vL_field, vR_field
        integer, intent(in)                                                                 :: dir, fbeg, fend
        integer                                                                             :: i, j, k, f
        real(wp)                                                                            :: mydx, dsum

        $:GPU_PARALLEL_LOOP(private='[i, j, k, f, mydx, dsum]', collapse=3, copyin='[dir, fbeg, fend]')
        do k = idwbuff(3)%beg, idwbuff(3)%end
            do j = idwbuff(2)%beg, idwbuff(2)%end
                do i = idwbuff(1)%beg, idwbuff(1)%end
                    if (dir == 1) then
                        mydx = dx(i)
                    else if (dir == 2) then
                        mydx = dy(j)
                    else
                        mydx = dz(k)
                    end if
                    dsum = 0._wp
                    do f = fbeg, fend
                        dsum = dsum + (vR_field(i, j, k, f) - vL_field(i, j, k, f))
                    end do
                    dq(i, j, k) = dsum/mydx
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_gradient_field

    !> Derivative of a cell-centered field along dir with the precomputed Fornberg finite-difference weights.
    subroutine s_gradient_dir_fornberg(q, dq, dir)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:), intent(in)  :: q
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:), intent(out) :: dq
        integer, intent(in)                                                                :: dir
        integer                                                                            :: i, j, k, a, npts, s_idx

        npts = (nWeights_grad - 1)/2

        if (dir == 1) then
            $:GPU_PARALLEL_LOOP(private='[i, j, k, s_idx, a]', collapse=3,copyin='[npts]')
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
        else if (dir == 2) then
            $:GPU_PARALLEL_LOOP(private='[i, j, k, s_idx, a]', collapse=3,copyin='[npts]')
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
        else if (dir == 3) then
            $:GPU_PARALLEL_LOOP(private='[i, j, k, s_idx, a]', collapse=3,copyin='[npts]')
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

    !> Fornberg first-derivative weights at every cell center in each direction, computed once on the host and copied to the device.
    impure subroutine s_compute_fornberg_fd_weights(npts)

        integer, intent(in) :: npts
        integer             :: i, j, k, a, m_order
        integer             :: s_idx
        real(wp)            :: x0, y0, z0
        real(wp)            :: x_stencil(nWeights_grad)
        real(wp)            :: c(nWeights_grad,0:1)

        m_order = 1  ! first derivative

        ! Computed once on the host: the run-time-sized stencil arrays must not be GPU-private
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

        if (num_dims == 3) then
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
        end if

        do i = 1, nWeights_grad
            $:GPU_UPDATE(device='[weights_x_grad(i)%sf, weights_y_grad(i)%sf, weights_z_grad(i)%sf]')
        end do

    end subroutine s_compute_fornberg_fd_weights

    !> Fornberg (1988) finite-difference weights at x0 on the given stencil, for derivatives up to m_order.
    subroutine s_fornberg_weights(x0, stencil, npts, m_order, coeffs)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)   :: npts  ! number of stencil points
        integer, intent(in)   :: m_order  ! highest derivative order
        real(wp), intent(in)  :: x0  ! evaluation point
        real(wp), intent(in)  :: stencil(npts)  ! stencil coordinates
        real(wp), intent(out) :: coeffs(npts,0:m_order)
        integer               :: i, j, k, mn
        real(wp)              :: c1, c2, c3, c4, c5

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

    !> Barycentric Lagrange interpolation weights at every cell in each direction, computed once.
    impure subroutine s_compute_barycentric_weights(npts)

        integer, intent(in) :: npts
        integer             :: i, j, k, a, b
        real(wp)            :: prod_x, prod_y, prod_z, dx_loc, dy_loc, dz_loc

        $:GPU_PARALLEL_LOOP(private='[i, a, b, prod_x, dx_loc]', copyin = '[npts]')
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

        $:GPU_PARALLEL_LOOP(private='[j, a, b, prod_y, dy_loc]', copyin = '[npts]')
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
            $:GPU_PARALLEL_LOOP(private='[k, a, b, prod_z, dz_loc]', copyin = '[npts]')
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
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_barycentric_weights

    !> Append each local particle's position, velocity, force, radius and fluid state at the particle to the evolution file.
    impure subroutine s_write_lag_particle_evol(qtime)

        real(wp), intent(in) :: qtime
        integer              :: k
        character(LEN=32)    :: FMT

        if (precision == 1) then
            FMT = "(ES17.8E3,I14,14ES17.8E3)"
        else
            FMT = "(ES25.16E3,I14,14ES25.16E3)"
        end if

        $:GPU_UPDATE(host='[lag_part_id, particle_pos, particle_vel, f_p, particle_rad, fluid_vel_at_particle, density_at_particle]')

        ! Cycle through list
        do k = 1, n_el_particles_loc
            write (LAG_EVOL_ID, FMT) qtime, lag_part_id(k, 1), particle_pos(k, 1, 1), particle_pos(k, 2, 1), particle_pos(k, 3, &
                   & 1), particle_vel(k, 1, 1), particle_vel(k, 2, 1), particle_vel(k, 3, 1), f_p(k, 1), f_p(k, 2), f_p(k, 3), &
                   & particle_rad(k), fluid_vel_at_particle(k, 1), fluid_vel_at_particle(k, 2), fluid_vel_at_particle(k, 3), &
                   & density_at_particle(k)
        end do

    end subroutine s_write_lag_particle_evol

    !> Copy the particle state and projected volume fraction to the host before a save, and abort if any particle position or
    !! velocity is NaN.
    impure subroutine s_sync_particles_for_save()

        integer :: k

        $:GPU_UPDATE(host='[lag_part_id, particle_pos, particle_posPrev, particle_vel, particle_rad, particle_mass, &
                     & particle_seed, fqs_fluct]')
        $:GPU_UPDATE(host='[q_particles(alphaf_id)%sf]')

        do k = 1, n_el_particles_loc
            if (any(ieee_is_nan(particle_pos(k,:,1))) .or. any(ieee_is_nan(particle_vel(k,:,1)))) then
                call s_mpi_abort("Particle position or velocity is NaN: check the particle input file, or reduce dt.")
            end if
        end do

    end subroutine s_sync_particles_for_save

    !> Write the particles in this rank's physical domain to restart_data/lustre_lag_particles_<t_step>.dat (parallel_io only), in
    !! the bubble column layout so post_process reads both.
    impure subroutine s_write_restart_lag_particles(t_step)

        integer, intent(in)                   :: t_step
        integer                               :: i, k, n_loc
        integer, dimension(3)                 :: cell
        real(wp), dimension(3)                :: s_loc
        real(wp), allocatable, dimension(:,:) :: io_data

        if (.not. parallel_io) return

        n_loc = count([(particle_in_domain_physical(particle_pos(k,1:3,1)), k=1, n_el_particles_loc)])
        allocate (io_data(max(1, n_loc),1:lag_io_vars))
        io_data = 0._wp  ! Unused columns are written as zero

        i = 0
        do k = 1, n_el_particles_loc
            if (.not. particle_in_domain_physical(particle_pos(k,1:3,1))) cycle
            i = i + 1
            io_data(i, 1) = real(lag_part_id(k, 1))
            io_data(i,2:4) = particle_pos(k,1:3,1)
            io_data(i,5:7) = particle_posPrev(k,1:3,1)
            io_data(i,8:10) = particle_vel(k,1:3,1)
            io_data(i, 11) = particle_rad(k)
            ! Particle volume fraction in the host cell, located from the current position
            cell = fd_number - buff_size
            call s_locate_cell(particle_pos(k,1:3,1), cell, s_loc)
            io_data(i, 12) = 1._wp - q_particles(alphaf_id)%sf(cell(1), cell(2), cell(3))
            ! Bubble column layout: R0 is the (constant) radius and the radius-ratio extremes are 1
            io_data(i, 13) = particle_rad(k)
            io_data(i,14:15) = 1._wp
            io_data(i,16:18) = fqs_fluct(k,1:3)  ! Drag-fluctuation state
            io_data(i, 19) = particle_mass(k)
            ! Seed as two 16-bit halves, exact in single and double precision
            io_data(i, 20) = real(ibits(particle_seed(k), 0, 16), wp)
            io_data(i, 21) = real(ibits(particle_seed(k), 16, 16), wp)
        end do

        call s_write_lag_restart('particles', t_step, io_data, n_loc)
        deallocate (io_data)

    end subroutine s_write_restart_lag_particles

    !> Copy particle src's state into slot dest, to compact the particle arrays after removals.
    impure subroutine s_copy_lag_particle(dest, src)

        integer, intent(in) :: src, dest

        p_owner_rank(dest) = p_owner_rank(src)
        particle_mass(dest) = particle_mass(src)
        particle_seed(dest) = particle_seed(src)
        lag_part_id(dest, 1) = lag_part_id(src, 1)
        particle_rad(dest) = particle_rad(src)
        particle_vel(dest,1:3,1:2) = particle_vel(src,1:3,1:2)
        particle_s(dest,1:3,1:2) = particle_s(src,1:3,1:2)
        particle_pos(dest,1:3,1:2) = particle_pos(src,1:3,1:2)
        particle_posPrev(dest,1:3,1:2) = particle_posPrev(src,1:3,1:2)
        f_p(dest,1:3) = f_p(src,1:3)
        fqs_fluct(dest,1:3) = fqs_fluct(src,1:3)
        particle_dposdt(dest,1:3,1:lag_num_ts) = particle_dposdt(src,1:3,1:lag_num_ts)
        particle_dveldt(dest,1:3,1:lag_num_ts) = particle_dveldt(src,1:3,1:lag_num_ts)

    end subroutine s_copy_lag_particle

    !> Close the particle output files and deallocate the particle arrays.
    impure subroutine s_finalize_particle_lagrangian_solver()

        integer :: i

        if (particle_params%write_void_evol) call s_close_void_evol
        if (particle_params%write_particles) call s_close_lag_evol()

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

        if (particle_params%added_mass_force > 0) then
            do i = 1, eqn_idx%mom%end
                @:DEALLOCATE(rhs_old(i)%sf)
            end do
        end if
        @:DEALLOCATE(rhs_old)

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

        ! Deallocating space
        @:DEALLOCATE(lag_part_id)
        @:DEALLOCATE(particle_mass)
        @:DEALLOCATE(particle_seed)
        @:DEALLOCATE(p_AM)
        @:DEALLOCATE(p_owner_rank)
        @:DEALLOCATE(particle_rad)
        @:DEALLOCATE(particle_pos)
        @:DEALLOCATE(particle_posPrev)
        @:DEALLOCATE(particle_vel)
        @:DEALLOCATE(particle_s)
        @:DEALLOCATE(particle_dposdt)
        @:DEALLOCATE(particle_dveldt)
        @:DEALLOCATE(f_p)
        @:DEALLOCATE(fqs_fluct)
        @:DEALLOCATE(gSum)

        @:DEALLOCATE(keep_particle)

        @:DEALLOCATE(fluid_vel_at_particle)
        @:DEALLOCATE(density_at_particle)
        @:DEALLOCATE(pres_at_particle)
        @:DEALLOCATE(force_status, force_re_mach)

    end subroutine s_finalize_particle_lagrangian_solver

end module m_particles_EL
