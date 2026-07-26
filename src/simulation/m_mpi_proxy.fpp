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
    !! EL Bubbles communication variables
    integer, parameter :: MAX_NEIGHBORS = 27
    integer            :: send_requests(MAX_NEIGHBORS), recv_requests(MAX_NEIGHBORS)
    integer            :: recv_offsets(MAX_NEIGHBORS)
    integer            :: neighbor_list(MAX_NEIGHBORS, 3)
    integer            :: n_neighbors
    $:GPU_DECLARE(create='[p_send_counts]')

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

    !! Initialize the MPI buffers and variables required for the particle communication.
    subroutine s_initialize_particles_mpi(lag_num_ts)

        integer, intent(in) :: lag_num_ts
        integer             :: i, j, k
        integer             :: real_size, int_size, nReal
        integer             :: ierr  !< Generic flag used to identify and report MPI errors

#ifdef MFC_MPI
        call MPI_Pack_size(1, mpi_p, MPI_COMM_WORLD, real_size, ierr)
        call MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, int_size, ierr)
        nReal = 7 + 16*2 + 10*lag_num_ts
        p_var_size = nReal*real_size + int_size
        p_buff_size = lag_params%nBubs_glb*p_var_size
        @:ALLOCATE(p_send_buff(0:p_buff_size), p_recv_buff(0:p_buff_size))
        @:ALLOCATE(p_send_ids(nidx(1)%beg:nidx(1)%end, nidx(2)%beg:nidx(2)%end, nidx(3)%beg:nidx(3)%end, 0:lag_params%nBubs_glb))
        ! First, collect all neighbor information
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
#endif

    end subroutine s_initialize_particles_mpi

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
            & 'bc_z%Twall_in', 'bc_z%Twall_out']
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
                & 'burn_rate_exp', 'burn_rate_pref']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            #:for VAR in ['vel', 'angular_vel', 'angles']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%moving_ibm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%airfoil_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%inj_species, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
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
                & 'radius', 'mass', 'min_spacing']
                call MPI_BCAST(particle_cloud(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(particle_cloud(i)%num_particles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%moving_ibm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%packing_method, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        end do

        ! manual: acoustic/probe/integral (combined loop; complex acoustic member set)
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

            #:for VAR in [ 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax' ]
                call MPI_BCAST(integral(j)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
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
    impure subroutine s_add_particles_to_transfer_list(nBub, pos, posPrev)

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
                call s_add_particle_to_direction(k, dx, dy, dz)
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

        subroutine s_add_particle_to_direction(particle_id, dir_x, dir_y, dir_z)

            integer, intent(in) :: particle_id, dir_x, dir_y, dir_z

            p_send_ids(dir_x, dir_y, dir_z, p_send_counts(dir_x, dir_y, dir_z)) = particle_id
            p_send_counts(dir_x, dir_y, dir_z) = p_send_counts(dir_x, dir_y, dir_z) + 1

        end subroutine s_add_particle_to_direction

    end subroutine s_add_particles_to_transfer_list

    !> Perform the MPI communication for lagrangian particles/bubbles.
    impure subroutine s_mpi_sendrecv_particles(bub_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, gas_betaC, bub_dphidt, lag_id, &
        & gas_p, gas_mv, rad, rvel, pos, posPrev, vel, scoord, drad, drvel, dgasp, dgasmv, dpos, dvel, lag_num_ts, nbubs, dest)

        real(wp), dimension(:)     :: bub_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, gas_betaC, bub_dphidt
        integer, dimension(:,:)    :: lag_id
        real(wp), dimension(:,:)   :: gas_p, gas_mv, rad, rvel, drad, drvel, dgasp, dgasmv
        real(wp), dimension(:,:,:) :: pos, posPrev, vel, scoord, dpos, dvel
        integer                    :: position, bub_id, lag_num_ts, tag, partner, send_tag, recv_tag, nbubs, p_recv_size, dest
        integer                    :: i, j, k, l, q, r
        integer                    :: req_send, req_recv, ierr  !< Generic flag used to identify and report MPI errors
        integer                    :: send_count, send_offset, recv_count, recv_offset

#ifdef MFC_MPI
        ! Phase 1: Exchange particle counts using non-blocking communication
        send_count = 0
        recv_count = 0

        ! Post all receives first
        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)
            partner = neighbor_ranks(i, j, k)
            recv_tag = neighbor_tag(i, j, k)

            recv_count = recv_count + 1
            call MPI_Irecv(p_recv_counts(i, j, k), 1, MPI_INTEGER, partner, recv_tag, MPI_COMM_WORLD, recv_requests(recv_count), &
                           & ierr)
        end do

        ! Post all sends
        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)
            partner = neighbor_ranks(i, j, k)
            send_tag = neighbor_tag(-i, -j, -k)

            send_count = send_count + 1
            call MPI_Isend(p_send_counts(i, j, k), 1, MPI_INTEGER, partner, send_tag, MPI_COMM_WORLD, send_requests(send_count), &
                           & ierr)
        end do

        ! Wait for all count exchanges to complete
        if (recv_count > 0) then
            call MPI_Waitall(recv_count, recv_requests(1:recv_count), MPI_STATUSES_IGNORE, ierr)
        end if
        if (send_count > 0) then
            call MPI_Waitall(send_count, send_requests(1:send_count), MPI_STATUSES_IGNORE, ierr)
        end if

        ! Phase 2: Exchange particle data using non-blocking communication
        send_count = 0
        recv_count = 0

        ! Post all receives for particle data first
        recv_offset = 1
        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)

            if (p_recv_counts(i, j, k) > 0) then
                partner = neighbor_ranks(i, j, k)
                p_recv_size = p_recv_counts(i, j, k)*p_var_size
                recv_tag = neighbor_tag(i, j, k)

                recv_count = recv_count + 1
                call MPI_Irecv(p_recv_buff(recv_offset), p_recv_size, MPI_PACKED, partner, recv_tag, MPI_COMM_WORLD, &
                               & recv_requests(recv_count), ierr)
                recv_offsets(l) = recv_offset
                recv_offset = recv_offset + p_recv_size
            end if
        end do

        ! Pack and send particle data
        send_offset = 0
        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)

            if (p_send_counts(i, j, k) > 0 .and. abs(i) + abs(j) + abs(k) /= 0) then
                partner = neighbor_ranks(i, j, k)
                send_tag = neighbor_tag(-i, -j, -k)

                ! Pack data for sending
                position = 0
                do q = 0, p_send_counts(i, j, k) - 1
                    bub_id = p_send_ids(i, j, k, q)

                    call MPI_Pack(lag_id(bub_id, 1), 1, MPI_INTEGER, p_send_buff(send_offset), p_buff_size, position, &
                                  & MPI_COMM_WORLD, ierr)
                    call MPI_Pack(bub_R0(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Rmax_stats(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                  & ierr)
                    call MPI_Pack(Rmin_stats(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                  & ierr)
                    call MPI_Pack(gas_mg(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(gas_betaT(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                  & ierr)
                    call MPI_Pack(gas_betaC(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                  & ierr)
                    call MPI_Pack(bub_dphidt(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                  & ierr)
                    do r = 1, 2
                        call MPI_Pack(gas_p(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, &
                                      & MPI_COMM_WORLD, ierr)
                        call MPI_Pack(gas_mv(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, &
                                      & MPI_COMM_WORLD, ierr)
                        call MPI_Pack(rad(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                      & ierr)
                        call MPI_Pack(rvel(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                      & ierr)
                        call MPI_Pack(pos(bub_id,:,r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                      & ierr)
                        call MPI_Pack(posPrev(bub_id,:,r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, &
                                      & MPI_COMM_WORLD, ierr)
                        call MPI_Pack(vel(bub_id,:,r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                      & ierr)
                        call MPI_Pack(scoord(bub_id,:,r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, &
                                      & MPI_COMM_WORLD, ierr)
                    end do
                    do r = 1, lag_num_ts
                        call MPI_Pack(drad(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, &
                                      & ierr)
                        call MPI_Pack(drvel(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, &
                                      & MPI_COMM_WORLD, ierr)
                        call MPI_Pack(dgasp(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, &
                                      & MPI_COMM_WORLD, ierr)
                        call MPI_Pack(dgasmv(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, &
                                      & MPI_COMM_WORLD, ierr)
                        call MPI_Pack(dpos(bub_id,:,r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, &
                                      & MPI_COMM_WORLD, ierr)
                        call MPI_Pack(dvel(bub_id,:,r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, &
                                      & MPI_COMM_WORLD, ierr)
                    end do
                end do

                send_count = send_count + 1
                call MPI_Isend(p_send_buff(send_offset), position, MPI_PACKED, partner, send_tag, MPI_COMM_WORLD, &
                               & send_requests(send_count), ierr)
                send_offset = send_offset + position
            end if
        end do

        ! Wait for all recvs for contiguous data to complete
        call MPI_Waitall(recv_count, recv_requests(1:recv_count), MPI_STATUSES_IGNORE, ierr)

        ! Process received data as it arrives
        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)

            if (p_recv_counts(i, j, k) > 0 .and. abs(i) + abs(j) + abs(k) /= 0) then
                p_recv_size = p_recv_counts(i, j, k)*p_var_size
                recv_offset = recv_offsets(l)

                position = 0
                ! Unpack received data
                do q = 0, p_recv_counts(i, j, k) - 1
                    nbubs = nbubs + 1
                    bub_id = nbubs
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, lag_id(bub_id, 1), 1, MPI_INTEGER, &
                                    & MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, bub_R0(bub_id), 1, mpi_p, MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, Rmax_stats(bub_id), 1, mpi_p, &
                                    & MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, Rmin_stats(bub_id), 1, mpi_p, &
                                    & MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_mg(bub_id), 1, mpi_p, MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_betaT(bub_id), 1, mpi_p, MPI_COMM_WORLD, &
                                    & ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_betaC(bub_id), 1, mpi_p, MPI_COMM_WORLD, &
                                    & ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, bub_dphidt(bub_id), 1, mpi_p, &
                                    & MPI_COMM_WORLD, ierr)
                    do r = 1, 2
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_p(bub_id, r), 1, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_mv(bub_id, r), 1, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, rad(bub_id, r), 1, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, rvel(bub_id, r), 1, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, pos(bub_id,:,r), 3, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, posPrev(bub_id,:,r), 3, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, vel(bub_id,:,r), 3, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, scoord(bub_id,:,r), 3, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                    end do
                    do r = 1, lag_num_ts
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, drad(bub_id, r), 1, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, drvel(bub_id, r), 1, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, dgasp(bub_id, r), 1, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, dgasmv(bub_id, r), 1, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, dpos(bub_id,:,r), 3, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, dvel(bub_id,:,r), 3, mpi_p, &
                                        & MPI_COMM_WORLD, ierr)
                    end do
                    lag_id(bub_id, 2) = bub_id
                end do
                recv_offset = recv_offset + p_recv_size
            end if
        end do

        ! Wait for all sends to complete
        if (send_count > 0) then
            call MPI_Waitall(send_count, send_requests(1:send_count), MPI_STATUSES_IGNORE, ierr)
        end if
#endif

        if (any(periodic_bc)) then
            call s_wrap_particle_positions(pos, posPrev, nbubs, dest)
        end if

    end subroutine s_mpi_sendrecv_particles

    !> Return a unique tag for each neighbor based on its position relative to the current process.
    integer function neighbor_tag(i, j, k) result(tag)

        integer, intent(in) :: i, j, k

        tag = (k + 1)*9 + (j + 1)*3 + (i + 1)

    end function neighbor_tag

    subroutine s_wrap_particle_positions(pos, posPrev, nbubs, dest)

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

    end subroutine s_wrap_particle_positions

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
