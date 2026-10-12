!>
!! @file
!! @brief Contains module m_mpi_proxy

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief MPI halo exchange, domain decomposition, and buffer packing/unpacking for the simulation solver
module m_mpi_proxy

#ifdef MFC_MPI
    use mpi  !< Message passing interface (MPI) module
#endif

    use m_helper_basic
    use m_helper
    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_nvtx
    use ieee_arithmetic

    implicit none

    integer, private, allocatable, dimension(:) :: ib_buff_send  !< IB marker send buffer for halo exchange
    integer, private, allocatable, dimension(:) :: ib_buff_recv  !< IB marker receive buffer for halo exchange
    integer                                     :: i_halo_size
    $:GPU_DECLARE(create='[i_halo_size]')

    integer, dimension(-1:1,-1:1,-1:1)          :: p_send_counts, p_recv_counts
    integer, dimension(:,:,:,:), allocatable    :: p_send_ids
    character(len=1), dimension(:), allocatable :: p_send_buff, p_recv_buff
    integer                                     :: p_buff_size, p_var_size
    !! Lagrangian entity communication variables
    integer, parameter :: MAX_NEIGHBORS = 27
    integer            :: send_requests(MAX_NEIGHBORS), recv_requests(MAX_NEIGHBORS)
    integer            :: neighbor_list(MAX_NEIGHBORS, 3)
    integer            :: n_neighbors
    $:GPU_DECLARE(create='[p_send_counts]')
    integer, allocatable, dimension(:) :: lag_send_ids    !< Local index of the entity packed into each send slot
    integer, allocatable, dimension(:) :: lag_recv_ranks  !< Rank each receive slot came from

contains

    !> Initialize the MPI proxy module
    subroutine s_initialize_mpi_proxy_module()

#ifdef MFC_MPI
        if (ib) then
            if (n > 0) then
                if (p > 0) then
                    i_halo_size = -1 + buff_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1)*(p + 2*buff_size + 1) &
                                                  & /(cells_bounds%mnp_min + 2*buff_size + 1)
                else
                    i_halo_size = -1 + buff_size*(cells_bounds%mn_max + 2*buff_size + 1)
                end if
            else
                i_halo_size = -1 + buff_size
            end if

            $:GPU_UPDATE(device='[i_halo_size]')
            @:ALLOCATE(ib_buff_send(0:i_halo_size), ib_buff_recv(0:i_halo_size))
        end if
#endif

    end subroutine s_initialize_mpi_proxy_module

    !> Size the MPI slot of one Lagrangian entity (an integer ID, n_real reals and n_int8 64-bit integers) and allocate the
    !! communication buffers for up to n_glb entities.
    subroutine s_initialize_lag_entity_mpi(n_real, n_int8, n_glb)

        integer, intent(in) :: n_real, n_int8, n_glb
        integer             :: real_size, int_size, int8_size
        integer             :: ierr  !< Generic flag used to identify and report MPI errors

#ifdef MFC_MPI
        call MPI_Pack_size(1, mpi_p, MPI_COMM_WORLD, real_size, ierr)
        call MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, int_size, ierr)
        call MPI_Pack_size(1, MPI_INTEGER8, MPI_COMM_WORLD, int8_size, ierr)
        p_var_size = n_real*real_size + int_size + n_int8*int8_size
        ! Start with room for an even share of the entities; s_begin_lag_entity_exchange grows it on demand
        p_buff_size = int(min(max(int(n_glb/num_procs, 8), 1_8)*int(p_var_size, 8), int(huge(p_buff_size), 8)))
        call s_allocate_lag_entity_comm(n_glb)
        @:ALLOCATE(lag_send_ids(1:n_glb), lag_recv_ranks(1:n_glb))
#endif

    end subroutine s_initialize_lag_entity_mpi

    !> Allocate the particle send/receive buffers (sized by p_buff_size) and build the list of neighbor ranks.
    subroutine s_allocate_lag_entity_comm(n_glb)

        integer, intent(in) :: n_glb  !< Global number of Lagrangian particles or bubbles
        integer             :: i, j, k

        @:ALLOCATE(p_send_buff(0:p_buff_size), p_recv_buff(0:p_buff_size))
        @:ALLOCATE(p_send_ids(nidx(1)%beg:nidx(1)%end, nidx(2)%beg:nidx(2)%end, nidx(3)%beg:nidx(3)%end, 0:n_glb))

        n_neighbors = 0
        do k = nidx(3)%beg, nidx(3)%end
            do j = nidx(2)%beg, nidx(2)%end
                do i = nidx(1)%beg, nidx(1)%end
                    if (abs(i) + abs(j) + abs(k) /= 0) then
                        n_neighbors = n_neighbors + 1
                        neighbor_list(n_neighbors, 1) = i
                        neighbor_list(n_neighbors, 2) = j
                        neighbor_list(n_neighbors, 3) = k
                    end if
                end do
            end do
        end do

    end subroutine s_allocate_lag_entity_comm

    !> Since only the processor with rank 0 reads and verifies the consistency of user inputs, these are initially not available to
    !! the other processors. Then, the purpose of this subroutine is to distribute the user inputs to the remaining processors in
    !! the communicator.
    impure subroutine s_mpi_bcast_user_inputs()

#ifdef MFC_MPI
        integer :: i, j  !< Generic loop iterator
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        ! Generated: case_dir, namelist scalars (INT/LOG/REAL), CASE_OPT guard, fluid_pp loop,
        !            bub_pp guard, lag_params guard, chem_params guard
        #:include 'generated_bcast.fpp'

        ! manual: m_glb, n_glb, p_glb (computed in s_read_input_file, not namelist-bound)
        call MPI_BCAST(m_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(n_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(p_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! manual: bc_x/y/z member broadcasts (struct members not in NAMELIST_VARS)
        #:for VAR in [ 'bc_x%beg', 'bc_x%end', 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end']
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'bc_x%grcbc_in', 'bc_x%grcbc_out', 'bc_x%grcbc_vel_out',  &
            & 'bc_y%grcbc_in', 'bc_y%grcbc_out', 'bc_y%grcbc_vel_out',            &
            & 'bc_z%grcbc_in', 'bc_z%grcbc_out', 'bc_z%grcbc_vel_out',            &
            & 'bc_x%isothermal_in', 'bc_y%isothermal_in', 'bc_z%isothermal_in',   &
            & 'bc_x%isothermal_out', 'bc_y%isothermal_out', 'bc_z%isothermal_out']
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'bc_x%vb1','bc_x%vb2','bc_x%vb3','bc_x%ve1','bc_x%ve2','bc_x%ve3', &
            & 'bc_y%vb1','bc_y%vb2','bc_y%vb3','bc_y%ve1','bc_y%ve2','bc_y%ve3',           &
            & 'bc_z%vb1','bc_z%vb2','bc_z%vb3','bc_z%ve1','bc_z%ve2','bc_z%ve3',           &
            & 'bc_x%pres_in','bc_x%pres_out','bc_y%pres_in','bc_y%pres_out',               &
            & 'bc_z%pres_in','bc_z%pres_out',                                               &
            & 'bc_x%Twall_in', 'bc_x%Twall_out', 'bc_y%Twall_in', 'bc_y%Twall_out',       &
            & 'bc_z%Twall_in', 'bc_z%Twall_out',                                           &
            & 'bc_x%vel_in_ramp', 'bc_x%vel_in_t0', 'bc_x%vel_in_frac0',                   &
            & 'bc_y%vel_in_ramp', 'bc_y%vel_in_t0', 'bc_y%vel_in_frac0',                   &
            & 'bc_z%vel_in_ramp', 'bc_z%vel_in_t0', 'bc_z%vel_in_frac0']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, 3
            #:for VAR in [ 'bc_x%vel_in', 'bc_x%vel_out', 'bc_y%vel_in', 'bc_y%vel_out',  &
                & 'bc_z%vel_in', 'bc_z%vel_out' ]
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: cfl_dt (runtime-computed logical), bc_io (BC-file existence)
        call MPI_BCAST(cfl_dt, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(bc_io, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        ! manual: shear_stress, bulk_stress (derived from Re_size post-init on all ranks),
        !         bodyForces (derived from bf_x/y/z)
        call MPI_BCAST(shear_stress, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(bulk_stress, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(bodyForces, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        ! manual: bc_x per-fluid inflow arrays (loop over num_fluids_max)
        do i = 1, num_fluids_max
            #:for VAR in ['bc_x%alpha_rho_in','bc_x%alpha_in','bc_y%alpha_rho_in','bc_y%alpha_in','bc_z%alpha_rho_in', &
                & 'bc_z%alpha_in']
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: patch_ib (sim member subset differs from pre; uses count=3, adds mass/moving_ibm)
        do i = 1, num_ibs
            #:for VAR in [ 'radius', 'length_x', 'length_y', 'length_z', &
                & 'x_centroid', 'y_centroid', 'z_centroid', 'slip', 'mass', 'v_blow', &
                & 'Twall', 'burn_rate_exp', 'burn_rate_pref']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            #:for VAR in ['vel', 'angular_vel', 'angles', 'kin_hinge', 'kin_offset']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            #:for VAR in ['kin_phi0', 'kin_theta0', 'kin_theta_mean', 'kin_freq', 'kin_phase', 'kin_t0', 'kin_ramp', &
                & 'kin_pitch_rate', 'kin_smooth']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%kin_model, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%moving_ibm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%airfoil_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%inj_species, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%thermal_bc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%surface_reaction, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        end do

        ! manual: ib_airfoil (kept manual alongside patch_ib)
        do i = 1, num_ib_airfoils_max
            #:for VAR in ['c', 'p', 't', 'm']
                call MPI_BCAST(ib_airfoil(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: stl_models loop (num_stl_models scalar is generated; grouped array members)
        do i = 1, num_stl_models_max
            call MPI_BCAST(stl_models(i)%model_filepath, len(stl_models(i)%model_filepath), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(stl_models(i)%model_threshold, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:for VAR in ['model_translate', 'model_scale']
                call MPI_BCAST(stl_models(i)%${VAR}$, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: particle_cloud (runtime loop to num_particle_clouds; irregular member subset)
        do i = 1, num_particle_clouds
            #:for VAR in ['x_centroid', 'y_centroid', 'z_centroid', 'length_x', 'length_y', 'length_z', &
                & 'radius', 'mass', 'min_spacing', 'shell_inner_radius', 'shell_outer_radius']
                call MPI_BCAST(particle_cloud(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(particle_cloud(i)%num_particles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%moving_ibm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%cloud_geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%shell_axis, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%packing_method, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%periodic, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        end do

        ! manual: acoustic/probe (combined loop; complex acoustic member set)
        do j = 1, num_probes_max
            do i = 1, 3
                call MPI_BCAST(acoustic(j)%loc(i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            end do

            call MPI_BCAST(acoustic(j)%dipole, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'pulse', 'support', 'num_elements', 'element_on', 'bb_num_freq' ]
                call MPI_BCAST(acoustic(j)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'mag', 'length', 'height', &
                'wavelength', 'frequency', 'gauss_sigma_dist', 'gauss_sigma_time', &
                'npulse', 'dir', 'delay', 'foc_length', 'aperture', &
                'element_spacing_angle', 'element_polygon_ratio', 'rotate_angle', &
                'bb_bandwidth', 'bb_lowest_freq' ]
                call MPI_BCAST(acoustic(j)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'x','y','z' ]
                call MPI_BCAST(probe(j)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: spatial-support body-force derived-type members (the bf_spatial_support toggle is broadcast by
        ! generated_bcast.fpp)
        #:for VAR in ['amp','x_centroid','y_centroid','sigma','conv_vel']
            call MPI_BCAST(spatial_bf%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor
        call MPI_BCAST(spatial_bf%freq, 8, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(spatial_bf%phase, 8, mpi_p, 0, MPI_COMM_WORLD, ierr)

        ! manual: synthetic turbulence namelist arrays (registered as indexed
        ! variants only; scalars are broadcast by generated_bcast.fpp)
        call MPI_BCAST(synth_n_waves_per_shell, num_synth_shells_max, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(synth_k_shell, num_synth_shells_max, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(synth_amp_shell, num_synth_shells_max, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(turb_pos, num_turb_sources_max*3, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(synth_L, num_turb_sources_max*3, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_bcast_user_inputs

    !> Adds particles to the transfer list for the MPI communication.
    !! @param nBub Current LOCAL number of bubbles
    !! @param pos Current position of each bubble
    !! @param posPrev Previous position of each bubble (optional, not used
    !!                for communication of initial condition)
    impure subroutine s_add_lag_entities_to_transfer_list(nBub, pos, posPrev)

        integer, intent(in)                  :: nBub
        real(wp), dimension(:,:), intent(in) :: pos, posPrev
        integer                              :: bubID
        integer                              :: i, j, k
        integer                              :: dx, dy, dz

        do k = nidx(3)%beg, nidx(3)%end
            do j = nidx(2)%beg, nidx(2)%end
                do i = nidx(1)%beg, nidx(1)%end
                    p_send_counts(i, j, k) = 0
                end do
            end do
        end do

        do k = 1, nbub
            dx = 0; dy = 0; dz = 0
            if (f_crosses_boundary(k, 1, -1, pos, posPrev)) then
                dx = -1
            else if (f_crosses_boundary(k, 1, 1, pos, posPrev)) then
                dx = 1
            end if
            if (n > 0) then
                if (f_crosses_boundary(k, 2, -1, pos, posPrev)) then
                    dy = -1
                else if (f_crosses_boundary(k, 2, 1, pos, posPrev)) then
                    dy = 1
                end if
            end if
            if (p > 0) then
                if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                    dz = -1
                else if (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                    dz = 1
                end if
            end if
            if (abs(dx) + abs(dy) + abs(dz) /= 0) then
                call s_add_lag_entity_to_direction(k, dx, dy, dz)
            end if
        end do

    contains

        logical function f_crosses_boundary(particle_id, dir, loc, pos, posPrev)

            integer, intent(in)                            :: particle_id, dir, loc
            real(wp), dimension(:,:), intent(in)           :: pos
            real(wp), dimension(:,:), optional, intent(in) :: posPrev

            if (loc == -1) then  ! Beginning of the domain
                if (nidx(dir)%beg == 0) then
                    f_crosses_boundary = .false.
                    return
                end if

                f_crosses_boundary = (posPrev(particle_id, dir) >= pcomm_coords(dir)%beg .and. pos(particle_id, &
                                      & dir) < pcomm_coords(dir)%beg)
            else if (loc == 1) then  ! End of the domain
                if (nidx(dir)%end == 0) then
                    f_crosses_boundary = .false.
                    return
                end if

                f_crosses_boundary = (posPrev(particle_id, dir) <= pcomm_coords(dir)%end .and. pos(particle_id, &
                                      & dir) > pcomm_coords(dir)%end)
            end if

        end function f_crosses_boundary

        subroutine s_add_lag_entity_to_direction(particle_id, dir_x, dir_y, dir_z)

            integer, intent(in) :: particle_id, dir_x, dir_y, dir_z

            p_send_ids(dir_x, dir_y, dir_z, p_send_counts(dir_x, dir_y, dir_z)) = particle_id
            p_send_counts(dir_x, dir_y, dir_z) = p_send_counts(dir_x, dir_y, dir_z) + 1

        end subroutine s_add_lag_entity_to_direction

    end subroutine s_add_lag_entities_to_transfer_list

    !> Exchange the number of particles to send/receive with every neighbor rank (p_send_counts -> p_recv_counts).
    impure subroutine s_exchange_lag_entity_counts()

        integer :: i, j, k, l, send_count, recv_count, ierr

#ifdef MFC_MPI
        send_count = 0
        recv_count = 0

        ! Post all receives first
        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)
            recv_count = recv_count + 1
            call MPI_Irecv(p_recv_counts(i, j, k), 1, MPI_INTEGER, neighbor_ranks(i, j, k), neighbor_tag(i, j, k), &
                           & MPI_COMM_WORLD, recv_requests(recv_count), ierr)
        end do

        ! Post all sends
        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)
            send_count = send_count + 1
            call MPI_Isend(p_send_counts(i, j, k), 1, MPI_INTEGER, neighbor_ranks(i, j, k), neighbor_tag(-i, -j, -k), &
                           & MPI_COMM_WORLD, send_requests(send_count), ierr)
        end do

        if (recv_count > 0) call MPI_Waitall(recv_count, recv_requests(1:recv_count), MPI_STATUSES_IGNORE, ierr)
        if (send_count > 0) call MPI_Waitall(send_count, send_requests(1:send_count), MPI_STATUSES_IGNORE, ierr)
#endif

    end subroutine s_exchange_lag_entity_counts

    !> Exchange the transfer counts with the neighbors, grow the buffers to fit, and list the entities to send in slot order
    !! (n_send) and the rank of each of the n_recv entities to receive.
    impure subroutine s_begin_lag_entity_exchange(n_send, n_recv)

        integer, intent(out) :: n_send, n_recv
        integer              :: i, j, k, l, q
        integer(kind=8)      :: max_bytes
        character(len=256)   :: mpi_dbg_msg

        n_send = 0
        n_recv = 0
#ifdef MFC_MPI
        call s_exchange_lag_entity_counts()

        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)
            do q = 0, p_send_counts(i, j, k) - 1
                n_send = n_send + 1
                lag_send_ids(n_send) = p_send_ids(i, j, k, q)
            end do
            lag_recv_ranks(n_recv + 1:n_recv + p_recv_counts(i, j, k)) = neighbor_ranks(i, j, k)
            n_recv = n_recv + p_recv_counts(i, j, k)
        end do

        ! Every message is part of these totals, so its size and offset fit in 32 bits and in the buffers
        max_bytes = int(max(n_send, n_recv), 8)*int(p_var_size, 8)
        if (max_bytes > int(huge(p_buff_size), 8)) then
            write (mpi_dbg_msg, '(a,i0,a,i0)') 'Lagrangian MPI buffer exceeds 32-bit MPI count limit on rank ', proc_rank, &
                   & ': size=', max_bytes
            call s_mpi_abort(trim(mpi_dbg_msg))
        end if
        if (max_bytes > p_buff_size) then
            p_buff_size = int(max_bytes)
            @:DEALLOCATE(p_send_buff, p_recv_buff)
            @:ALLOCATE(p_send_buff(0:p_buff_size), p_recv_buff(0:p_buff_size))
        end if
#endif

    end subroutine s_begin_lag_entity_exchange

    !> Send the packed slots to the neighbors and receive theirs, in the slot order of s_begin_lag_entity_exchange.
    impure subroutine s_exchange_lag_entities()

        integer :: i, j, k, l, ierr
        integer :: send_count, recv_count, send_slot, recv_slot

#ifdef MFC_MPI
        send_count = 0
        recv_count = 0
        send_slot = 0
        recv_slot = 0
        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)
            if (p_recv_counts(i, j, k) > 0) then
                recv_count = recv_count + 1
                call MPI_Irecv(p_recv_buff(recv_slot*p_var_size), p_recv_counts(i, j, k)*p_var_size, MPI_PACKED, &
                               & neighbor_ranks(i, j, k), neighbor_tag(i, j, k), MPI_COMM_WORLD, recv_requests(recv_count), ierr)
                recv_slot = recv_slot + p_recv_counts(i, j, k)
            end if
            if (p_send_counts(i, j, k) > 0) then
                send_count = send_count + 1
                call MPI_Isend(p_send_buff(send_slot*p_var_size), p_send_counts(i, j, k)*p_var_size, MPI_PACKED, &
                               & neighbor_ranks(i, j, k), neighbor_tag(-i, -j, -k), MPI_COMM_WORLD, send_requests(send_count), ierr)
                send_slot = send_slot + p_send_counts(i, j, k)
            end if
        end do

        call MPI_Waitall(recv_count, recv_requests(1:recv_count), MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall(send_count, send_requests(1:send_count), MPI_STATUSES_IGNORE, ierr)
#endif

    end subroutine s_exchange_lag_entities

    !> Return a unique tag for each neighbor based on its position relative to the current process.
    integer function neighbor_tag(i, j, k) result(tag)

        integer, intent(in) :: i, j, k

        tag = (k + 1)*9 + (j + 1)*3 + (i + 1)

    end function neighbor_tag

    subroutine s_wrap_lag_entity_positions(pos, posPrev, nbubs, dest)

        real(wp), dimension(:,:,:) :: pos, posPrev
        integer                    :: nbubs, dest
        integer                    :: i, q
        real(wp)                   :: offset

        do i = 1, nbubs
            if (periodic_bc(1)) then
                offset = glb_bounds(1)%end - glb_bounds(1)%beg
                if (pos(i, 1, dest) > x_cb(m + buff_size)) then
                    do q = 1, 2
                        pos(i, 1, q) = pos(i, 1, q) - offset
                        posPrev(i, 1, q) = posPrev(i, 1, q) - offset
                    end do
                end if
                if (pos(i, 1, dest) < x_cb(-1 - buff_size)) then
                    do q = 1, 2
                        pos(i, 1, q) = pos(i, 1, q) + offset
                        posPrev(i, 1, q) = posPrev(i, 1, q) + offset
                    end do
                end if
            end if

            if (periodic_bc(2)) then
                offset = glb_bounds(2)%end - glb_bounds(2)%beg
                if (pos(i, 2, dest) > y_cb(n + buff_size)) then
                    do q = 1, 2
                        pos(i, 2, q) = pos(i, 2, q) - offset
                        posPrev(i, 2, q) = posPrev(i, 2, q) - offset
                    end do
                end if
                if (pos(i, 2, dest) < y_cb(-buff_size - 1)) then
                    do q = 1, 2
                        pos(i, 2, q) = pos(i, 2, q) + offset
                        posPrev(i, 2, q) = posPrev(i, 2, q) + offset
                    end do
                end if
            end if

            if (periodic_bc(3)) then
                offset = glb_bounds(3)%end - glb_bounds(3)%beg
                if (pos(i, 3, dest) > z_cb(p + buff_size)) then
                    do q = 1, 2
                        pos(i, 3, q) = pos(i, 3, q) - offset
                        posPrev(i, 3, q) = posPrev(i, 3, q) - offset
                    end do
                end if
                if (pos(i, 3, dest) < z_cb(-1 - buff_size)) then
                    do q = 1, 2
                        pos(i, 3, q) = pos(i, 3, q) + offset
                        posPrev(i, 3, q) = posPrev(i, 3, q) + offset
                    end do
                end if
            end if
        end do

    end subroutine s_wrap_lag_entity_positions

    !> Broadcast random phase numbers from rank 0 to all MPI processes
    impure subroutine s_mpi_send_random_number(phi_rn, num_freq)

        integer, intent(in)                            :: num_freq
        real(wp), intent(inout), dimension(1:num_freq) :: phi_rn

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors
        call MPI_BCAST(phi_rn, num_freq, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_send_random_number

    !> Finalize the MPI proxy module
    subroutine s_finalize_mpi_proxy_module()

#ifdef MFC_MPI
        if (ib) then
            @:DEALLOCATE(ib_buff_send, ib_buff_recv)
        end if

        if (allocated(p_send_buff)) then
            @:DEALLOCATE(p_send_buff, p_recv_buff)
            @:DEALLOCATE(p_send_ids)
        end if
#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
