!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief The module serves as a proxy to the parameters and subroutines
!!          available in the MPI implementation's MPI module. Specifically,
!!          the purpose of the proxy is to harness basic MPI commands into
!!          more complicated procedures as to accomplish the communication
!!          goals for the simulation.
module m_mpi_proxy

#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_common

    use m_nvtx

    use ieee_arithmetic

    implicit none

    integer, private, allocatable, dimension(:) :: ib_buff_send !<
    !! This variable is utilized to pack and send the buffer of the immersed
    !! boundary markers, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    integer, private, allocatable, dimension(:) :: ib_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the
    !! immersed boundary markers, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    integer :: i_halo_size
    $:GPU_DECLARE(create='[i_halo_size]')

    integer, dimension(-1:1, -1:1, -1:1) :: p_send_counts, p_recv_counts
    integer, dimension(:, :, :, :), allocatable :: p_send_ids
    character(len=1), dimension(:), allocatable :: p_send_buff, p_recv_buff
    integer :: p_buff_size, p_var_size
    !! EL Bubbles communication variables
    integer, parameter :: MAX_NEIGHBORS = 27
    integer :: send_requests(MAX_NEIGHBORS), recv_requests(MAX_NEIGHBORS)
    integer :: recv_offsets(MAX_NEIGHBORS)
    integer :: neighbor_list(MAX_NEIGHBORS, 3)
    integer :: n_neighbors
    $:GPU_DECLARE(create='[p_send_counts]')

contains

    subroutine s_initialize_mpi_proxy_module()

#ifdef MFC_MPI
        if (ib) then
            if (n > 0) then
                if (p > 0) then
                    i_halo_size = -1 + buff_size* &
                                            & (m + 2*buff_size + 1)* &
                                            & (n + 2*buff_size + 1)* &
                                            & (p + 2*buff_size + 1)/ &
                                            & (cells_bounds%mnp_min + 2*buff_size + 1)
                else
                    i_halo_size = -1 + buff_size* &
                                            & (cells_bounds%mn_max + 2*buff_size + 1)
                end if
            else
                i_halo_size = -1 + buff_size
            end if

            $:GPU_UPDATE(device='[i_halo_size]')
            @:ALLOCATE(ib_buff_send(0:i_halo_size), ib_buff_recv(0:i_halo_size))
        end if
#endif

    end subroutine s_initialize_mpi_proxy_module

    !! This subroutine initializes the MPI buffers and variables
        !! required for the particle communication.
        !! @param lag_num_ts Number of stages in time-stepping scheme
    subroutine s_initialize_particles_mpi(lag_num_ts)

        integer :: i, j, k
        integer :: real_size, int_size, nReal, lag_num_ts
        integer :: ierr !< Generic flag used to identify and report MPI errors

#ifdef MFC_MPI
        call MPI_Pack_size(1, mpi_p, MPI_COMM_WORLD, real_size, ierr)
        call MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, int_size, ierr)
        nReal = 7 + 16*2 + 10*lag_num_ts
        p_var_size = 20*(nReal*real_size + int_size)
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

    !>  Since only the processor with rank 0 reads and verifies
        !!      the consistency of user inputs, these are initially not
        !!      available to the other processors. Then, the purpose of
        !!      this subroutine is to distribute the user inputs to the
        !!      remaining processors in the communicator.
    impure subroutine s_mpi_bcast_user_inputs()

#ifdef MFC_MPI

        integer :: i, j !< Generic loop iterator
        integer :: ierr !< Generic flag used to identify and report MPI errors

        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in ['k_x', 'k_y', 'k_z', 'w_x', 'w_y', 'w_z', 'p_x', 'p_y', &
            & 'p_z', 'g_x', 'g_y', 'g_z']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in ['t_step_old', 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',  &
            & 't_step_start','t_step_stop','t_step_save','t_step_print',       &
            & 'model_eqns','time_stepper', 'riemann_solver', 'low_Mach',       &
            & 'wave_speeds', 'avg_state', 'precision', 'bc_x%beg', 'bc_x%end', &
            & 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end',  'fd_order',     &
            & 'num_probes', 'num_integrals', 'bubble_model', 'thermal',        &
            & 'num_source', 'relax_model', 'num_ibs', 'n_start',    &
            & 'num_bc_patches', 'num_igr_iters', 'num_igr_warm_start_iters', &
            & 'adap_dt_max_iters' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'run_time_info','cyl_coord', 'mpp_lim',     &
            &  'mp_weno', 'rdma_mpi', 'powell', 'cont_damage', 'bc_io', &
            & 'weno_Re_flux', 'alt_soundspeed', 'null_weights', 'mixture_err',   &
            & 'parallel_io', 'hypoelasticity', 'bubbles_euler', 'polytropic',    &
            & 'polydisperse', 'qbmm', 'acoustic_source', 'probe_wrt', 'integral_wrt',   &
            & 'prim_vars_wrt', 'weno_avg', 'file_per_process', 'relax',          &
            & 'adv_n', 'adap_dt', 'ib', 'bodyForces', 'bf_x', 'bf_y', 'bf_z',    &
            & 'bc_x%grcbc_in', 'bc_x%grcbc_out', 'bc_x%grcbc_vel_out',          &
            & 'bc_y%grcbc_in', 'bc_y%grcbc_out', 'bc_y%grcbc_vel_out',          &
            & 'bc_z%grcbc_in', 'bc_z%grcbc_out', 'bc_z%grcbc_vel_out',          &
            & 'cfl_adap_dt', 'cfl_const_dt', 'cfl_dt', 'surface_tension',       &
            & 'shear_stress', 'bulk_stress', 'bubbles_lagrange',                &
            & 'hyperelasticity', 'down_sample', 'int_comp','fft_wrt' ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        if (chemistry) then
            #:for VAR in [ 'diffusion', 'reactions' ]
                call MPI_BCAST(chem_params%${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'gamma_method' ]
                call MPI_BCAST(chem_params%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end if

        if (bubbles_lagrange) then
            #:for VAR in [ 'heatTransfer_model', 'massTransfer_model', 'pressure_corrector', &
                & 'write_bubbles', 'write_bubbles_stats', 'write_void_evol', 'pressure_force', &
                & 'gravity_force']
                call MPI_BCAST(lag_params%${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in ['solver_approach', 'cluster_type', 'smooth_type', 'nBubs_glb', 'vel_model', &
                & 'drag_model']
                call MPI_BCAST(lag_params%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'c0', 'rho0', 'T0', 'x0', 'diffcoefvap', 'epsilonb','charwidth', &
                & 'valmaxvoid', 'Thost']
                call MPI_BCAST(lag_params%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(lag_params%input_path, len(lag_params%input_path), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        end if

        #:for VAR in [ 'dt','weno_eps','teno_CT','pref','rhoref','R0ref','Web','Ca', 'sigma', &
            & 'Re_inv', 'poly_sigma', 'palpha_eps', 'ptgalpha_eps', 'pi_fac',    &
            & 'bc_x%vb1','bc_x%vb2','bc_x%vb3','bc_x%ve1','bc_x%ve2','bc_x%ve2', &
            & 'bc_y%vb1','bc_y%vb2','bc_y%vb3','bc_y%ve1','bc_y%ve2','bc_y%ve3', &
            & 'bc_z%vb1','bc_z%vb2','bc_z%vb3','bc_z%ve1','bc_z%ve2','bc_z%ve3', &
            & 'bc_x%pres_in','bc_x%pres_out','bc_y%pres_in','bc_y%pres_out', 'bc_z%pres_in','bc_z%pres_out', &
            & 't_stop', 't_save', 'cfl_target', 'Bx0', 'alf_factor',  &
            & 'tau_star', 'cont_damage_s', 'alpha_bar', 'adap_dt_tol', &
            & 'ic_eps', 'ic_beta' ]
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, 3
            #:for VAR in [ 'bc_x%vel_in', 'bc_x%vel_out', 'bc_y%vel_in', 'bc_y%vel_out',  &
                & 'bc_z%vel_in', 'bc_z%vel_out']
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        #:if not MFC_CASE_OPTIMIZATION
            call MPI_BCAST(mapped_weno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wenoz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(teno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(num_fluids, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wenoz_q, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(mhd, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(relativity, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr_pres_lim, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr_iter_solver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(viscous, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(recon_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(muscl_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(muscl_lim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endif

        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf','mul0','ss','pv','gamma_v','M_v',  &
                & 'mu_v','k_v', 'cp_v','G', 'cv', 'qv', 'qvp' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(fluid_pp(i)%Re(1), 2, mpi_p, 0, MPI_COMM_WORLD, ierr)
        end do

        do i = 1, num_fluids_max
            #:for VAR in ['bc_x%alpha_rho_in','bc_x%alpha_in','bc_y%alpha_rho_in','bc_y%alpha_in','bc_z%alpha_rho_in','bc_z%alpha_in']
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        do i = 1, num_ibs
            #:for VAR in [ 'radius', 'length_x', 'length_y', &
                & 'x_centroid', 'y_centroid', 'c', 'm', 'p', 't', 'theta', 'slip']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            #:for VAR in ['vel', 'angular_vel', 'angles']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%moving_ibm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        end do

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

        ! NVIDIA UVM variables
        call MPI_BCAST(nv_uvm_out_of_core, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nv_uvm_igr_temps_on_gpu, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nv_uvm_pref_gpu, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        ! Extra BC Variable
        call MPI_BCAST(periodic_bc, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_bcast_user_inputs

    subroutine s_mpi_sendrecv_ib_buffers(ib_markers, mpi_dir, pbc_loc)

        type(integer_field), intent(inout) :: ib_markers

        integer, intent(in) :: mpi_dir, pbc_loc

        integer :: i, j, k, l, r, q !< Generic loop iterators

        integer :: buffer_counts(1:3), buffer_count

        type(int_bounds_info) :: boundary_conditions(1:3)
        integer :: beg_end(1:2), grid_dims(1:3)
        integer :: dst_proc, src_proc, recv_tag, send_tag

        logical :: beg_end_geq_0, qbmm_comm

        integer :: pack_offset, unpack_offset

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        call nvtxStartRange("IB-MARKER-COMM-PACKBUF")

        buffer_counts = (/ &
                        buff_size*(n + 1)*(p + 1), &
                        buff_size*(m + 2*buff_size + 1)*(p + 1), &
                        buff_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                        /)

        buffer_count = buffer_counts(mpi_dir)
        boundary_conditions = (/bc_x, bc_y, bc_z/)
        beg_end = (/boundary_conditions(mpi_dir)%beg, boundary_conditions(mpi_dir)%end/)
        beg_end_geq_0 = beg_end(max(pbc_loc, 0) - pbc_loc + 1) >= 0

        ! Implements:
        ! pbc_loc  bc_x >= 0 -> [send/recv]_tag  [dst/src]_proc
        ! -1 (=0)      0            ->     [1,0]       [0,0]      | 0 0 [1,0] [beg,beg]
        ! -1 (=0)      1            ->     [0,0]       [1,0]      | 0 1 [0,0] [end,beg]
        ! +1 (=1)      0            ->     [0,1]       [1,1]      | 1 0 [0,1] [end,end]
        ! +1 (=1)      1            ->     [1,1]       [0,1]      | 1 1 [1,1] [beg,end]

        send_tag = f_logical_to_int(.not. f_xor(beg_end_geq_0, pbc_loc == 1))
        recv_tag = f_logical_to_int(pbc_loc == 1)

        dst_proc = beg_end(1 + f_logical_to_int(f_xor(pbc_loc == 1, beg_end_geq_0)))
        src_proc = beg_end(1 + f_logical_to_int(pbc_loc == 1))

        grid_dims = (/m, n, p/)

        pack_offset = 0
        if (f_xor(pbc_loc == 1, beg_end_geq_0)) then
            pack_offset = grid_dims(mpi_dir) - buff_size + 1
        end if

        unpack_offset = 0
        if (pbc_loc == 1) then
            unpack_offset = grid_dims(mpi_dir) + buff_size + 1
        end if

        ! Pack Buffer to Send
        #:for mpi_dir in [1, 2, 3]
            if (mpi_dir == ${mpi_dir}$) then
                #:if mpi_dir == 1
                    #:call GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                        do l = 0, p
                            do k = 0, n
                                do j = 0, buff_size - 1
                                    r = (j + buff_size*(k + (n + 1)*l))
                                    ib_buff_send(r) = ib_markers%sf(j + pack_offset, k, l)
                                end do
                            end do
                        end do
                    #:endcall GPU_PARALLEL_LOOP
                #:elif mpi_dir == 2
                    #:call GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                        do l = 0, p
                            do k = 0, buff_size - 1
                                do j = -buff_size, m + buff_size
                                    r = ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         (k + buff_size*l))
                                    ib_buff_send(r) = ib_markers%sf(j, k + pack_offset, l)
                                end do
                            end do
                        end do
                    #:endcall GPU_PARALLEL_LOOP
                #:else
                    #:call GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                        do l = 0, buff_size - 1
                            do k = -buff_size, n + buff_size
                                do j = -buff_size, m + buff_size
                                    r = ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                    ib_buff_send(r) = ib_markers%sf(j, k, l + pack_offset)
                                end do
                            end do
                        end do
                    #:endcall GPU_PARALLEL_LOOP
                #:endif
            end if
        #:endfor
        call nvtxEndRange ! Packbuf

        #:for rdma_mpi in [False, True]
            if (rdma_mpi .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
                #:if rdma_mpi
                    #:call GPU_HOST_DATA(use_device_addr='[ib_buff_send, ib_buff_recv]')

                        call nvtxStartRange("IB-MARKER-SENDRECV-RDMA")
                        call MPI_SENDRECV( &
                            ib_buff_send, buffer_count, MPI_INTEGER, dst_proc, send_tag, &
                            ib_buff_recv, buffer_count, MPI_INTEGER, src_proc, recv_tag, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                        call nvtxEndRange

                    #:endcall GPU_HOST_DATA
                    $:GPU_WAIT()
                #:else
                    call nvtxStartRange("IB-MARKER-DEV2HOST")
                    $:GPU_UPDATE(host='[ib_buff_send]')
                    call nvtxEndRange

                    call nvtxStartRange("IB-MARKER-SENDRECV-NO-RMDA")
                    call MPI_SENDRECV( &
                        ib_buff_send, buffer_count, MPI_INTEGER, dst_proc, send_tag, &
                        ib_buff_recv, buffer_count, MPI_INTEGER, src_proc, recv_tag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                    call nvtxEndRange

                    call nvtxStartRange("IB-MARKER-HOST2DEV")
                    $:GPU_UPDATE(device='[ib_buff_recv]')
                    call nvtxEndRange
                #:endif
            end if
        #:endfor

        ! Unpack Received Buffer
        call nvtxStartRange("IB-MARKER-COMM-UNPACKBUF")
        #:for mpi_dir in [1, 2, 3]
            if (mpi_dir == ${mpi_dir}$) then
                #:if mpi_dir == 1
                    #:call GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                        do l = 0, p
                            do k = 0, n
                                do j = -buff_size, -1
                                    r = (j + buff_size*((k + 1) + (n + 1)*l))
                                    ib_markers%sf(j + unpack_offset, k, l) = ib_buff_recv(r)
                                end do
                            end do
                        end do
                    #:endcall GPU_PARALLEL_LOOP
                #:elif mpi_dir == 2
                    #:call GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                        do l = 0, p
                            do k = -buff_size, -1
                                do j = -buff_size, m + buff_size
                                    r = ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + buff_size*l))
                                    ib_markers%sf(j, k + unpack_offset, l) = ib_buff_recv(r)
                                end do
                            end do
                        end do
                    #:endcall GPU_PARALLEL_LOOP
                #:else
                    ! Unpacking buffer from bc_z%beg
                    #:call GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                        do l = -buff_size, -1
                            do k = -buff_size, n + buff_size
                                do j = -buff_size, m + buff_size
                                    r = ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + (n + 2*buff_size + 1)* &
                                          (l + buff_size)))
                                    ib_markers%sf(j, k, l + unpack_offset) = ib_buff_recv(r)
                                end do
                            end do
                        end do
                    #:endcall GPU_PARALLEL_LOOP
                #:endif
            end if
        #:endfor
        call nvtxEndRange
#endif

    end subroutine s_mpi_sendrecv_ib_buffers

    !> This subroutine adds particles to the transfer list for the MPI
        !! communication.
        !! @param nBub Current LOCAL number of bubbles
        !! @param pos Current position of each bubble
        !! @param posPrev Previous position of each bubble (optional, not used
        !!                for communication of initial condition)
    impure subroutine s_add_particles_to_transfer_list(nBub, pos, posPrev)

        real(wp), dimension(:, :) :: pos
        real(wp), dimension(:, :), optional :: posPrev
        integer :: bubID, nbub
        integer :: i, j, k

        do k = nidx(3)%beg, nidx(3)%end
            do j = nidx(2)%beg, nidx(2)%end
                do i = nidx(1)%beg, nidx(1)%end
                    p_send_counts(i, j, k) = 0
                end do
            end do
        end do

        do k = 1, nbub
            if (f_crosses_boundary(k, 1, -1, pos, posPrev)) then
                call s_add_particle_to_direction(k, -1, 0, 0)
                if (n > 0) then
                    if (f_crosses_boundary(k, 2, -1, pos, posPrev)) then
                        call s_add_particle_to_direction(k, -1, -1, 0)
                        call s_add_particle_to_direction(k, 0, -1, 0)
                        if (p > 0) then
                            if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, -1, -1, -1)
                                call s_add_particle_to_direction(k, 0, -1, -1)
                                call s_add_particle_to_direction(k, -1, 0, -1)
                                call s_add_particle_to_direction(k, 0, 0, -1)
                            elseif (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, -1, -1, 1)
                                call s_add_particle_to_direction(k, 0, -1, 1)
                                call s_add_particle_to_direction(k, -1, 0, 1)
                                call s_add_particle_to_direction(k, 0, 0, 1)
                            end if
                        end if
                    elseif (f_crosses_boundary(k, 2, 1, pos, posPrev)) then
                        call s_add_particle_to_direction(k, -1, 1, 0)
                        call s_add_particle_to_direction(k, 0, 1, 0)
                        if (p > 0) then
                            if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, -1, 1, -1)
                                call s_add_particle_to_direction(k, 0, 1, -1)
                                call s_add_particle_to_direction(k, -1, 0, -1)
                                call s_add_particle_to_direction(k, 0, 0, -1)
                            elseif (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, -1, 1, 1)
                                call s_add_particle_to_direction(k, 0, 1, 1)
                                call s_add_particle_to_direction(k, -1, 0, 1)
                                call s_add_particle_to_direction(k, 0, 0, 1)
                            end if
                        end if
                    else
                        if (p > 0) then
                            if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, -1, 0, -1)
                                call s_add_particle_to_direction(k, 0, 0, -1)
                            elseif (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, -1, 0, 1)
                                call s_add_particle_to_direction(k, 0, 0, 1)
                            end if
                        end if
                    end if
                end if
            elseif (f_crosses_boundary(k, 1, 1, pos, posPrev)) then
                call s_add_particle_to_direction(k, 1, 0, 0)
                if (n > 0) then
                    if (f_crosses_boundary(k, 2, -1, pos, posPrev)) then
                        call s_add_particle_to_direction(k, 1, -1, 0)
                        call s_add_particle_to_direction(k, 0, -1, 0)
                        if (p > 0) then
                            if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, 1, -1, -1)
                                call s_add_particle_to_direction(k, 0, -1, -1)
                                call s_add_particle_to_direction(k, 1, 0, -1)
                                call s_add_particle_to_direction(k, 0, 0, -1)
                            elseif (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, 1, -1, 1)
                                call s_add_particle_to_direction(k, 0, -1, 1)
                                call s_add_particle_to_direction(k, 1, 0, 1)
                                call s_add_particle_to_direction(k, 0, 0, 1)
                            end if
                        end if
                    elseif (f_crosses_boundary(k, 2, 1, pos, posPrev)) then
                        call s_add_particle_to_direction(k, 1, 1, 0)
                        call s_add_particle_to_direction(k, 0, 1, 0)
                        if (p > 0) then
                            if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, 1, 1, -1)
                                call s_add_particle_to_direction(k, 0, 1, -1)
                                call s_add_particle_to_direction(k, 1, 0, -1)
                                call s_add_particle_to_direction(k, 0, 0, -1)
                            elseif (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, 1, 1, 1)
                                call s_add_particle_to_direction(k, 0, 1, 1)
                                call s_add_particle_to_direction(k, 1, 0, 1)
                                call s_add_particle_to_direction(k, 0, 0, 1)
                            end if
                        end if
                    else
                        if (p > 0) then
                            if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, 1, 0, -1)
                                call s_add_particle_to_direction(k, 0, 0, -1)
                            elseif (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                                call s_add_particle_to_direction(k, 1, 0, 1)
                                call s_add_particle_to_direction(k, 0, 0, 1)
                            end if
                        end if
                    end if
                end if
            elseif (f_crosses_boundary(k, 2, -1, pos, posPrev)) then
                call s_add_particle_to_direction(k, 0, -1, 0)
                if (p > 0) then
                    if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                        call s_add_particle_to_direction(k, 0, -1, -1)
                        call s_add_particle_to_direction(k, 0, 0, -1)
                    elseif (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                        call s_add_particle_to_direction(k, 0, -1, 1)
                        call s_add_particle_to_direction(k, 0, 0, 1)
                    end if
                end if
            elseif (f_crosses_boundary(k, 2, 1, pos, posPrev)) then
                call s_add_particle_to_direction(k, 0, 1, 0)
                if (p > 0) then
                    if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                        call s_add_particle_to_direction(k, 0, 1, -1)
                        call s_add_particle_to_direction(k, 0, 0, -1)
                    elseif (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                        call s_add_particle_to_direction(k, 0, 1, 1)
                        call s_add_particle_to_direction(k, 0, 0, 1)
                    end if
                end if
            elseif (p > 0) then
                if (f_crosses_boundary(k, 3, -1, pos, posPrev)) then
                    call s_add_particle_to_direction(k, 0, 0, -1)
                elseif (f_crosses_boundary(k, 3, 1, pos, posPrev)) then
                    call s_add_particle_to_direction(k, 0, 0, 1)
                end if
            end if

        end do

    contains

        logical function f_crosses_boundary(particle_id, dir, loc, pos, posPrev)

            integer, intent(in) :: particle_id, dir, loc
            real(wp), dimension(:, :), intent(in) :: pos
            real(wp), dimension(:, :), optional, intent(in) :: posPrev

            if (loc == -1) then ! Beginning of the domain
                if (nidx(dir)%beg == 0) then
                    f_crosses_boundary = .false.
                    return
                end if

                if (present(posPrev)) then
                    f_crosses_boundary = (posPrev(particle_id, dir) > pcomm_coords(dir)%beg .and. &
                                          pos(particle_id, dir) < pcomm_coords(dir)%beg)
                else
                    f_crosses_boundary = (pos(particle_id, dir) < pcomm_coords(dir)%beg)
                end if
            elseif (loc == 1) then ! End of the domain
                if (nidx(dir)%end == 0) then
                    f_crosses_boundary = .false.
                    return
                end if

                if (present(posPrev)) then
                    f_crosses_boundary = (posPrev(particle_id, dir) < pcomm_coords(dir)%end .and. &
                                          pos(particle_id, dir) > pcomm_coords(dir)%end)
                else
                    f_crosses_boundary = (pos(particle_id, dir) > pcomm_coords(dir)%end)
                end if
            end if

        end function f_crosses_boundary

        subroutine s_add_particle_to_direction(particle_id, dir_x, dir_y, dir_z)

            integer, intent(in) :: particle_id, dir_x, dir_y, dir_z

            p_send_ids(dir_x, dir_y, dir_z, p_send_counts(dir_x, dir_y, dir_z)) = particle_id
            p_send_counts(dir_x, dir_y, dir_z) = p_send_counts(dir_x, dir_y, dir_z) + 1

        end subroutine s_add_particle_to_direction

    end subroutine s_add_particles_to_transfer_list

    !> This subroutine performs the MPI communication for lagrangian particles/
        !! bubbles.
        !! @param bub_R0 Initial radius of each bubble
        !! @param Rmax_stats Maximum radius of each bubble
        !! @param Rmin_stats Minimum radius of each bubble
        !! @param gas_mg Mass of gas in each bubble
        !! @param gas_betaT Heat flux model coefficient for each bubble
        !! @param gas_betaC mass flux model coefficient for each bubble
        !! @param bub_dphidt Subgrid velocity potential for each bubble
        !! @param lag_id Global and local ID of each bubble
        !! @param gas_p Pressure of the gas in each bubble
        !! @param gas_mv Mass of vapor in each bubble
        !! @param rad Radius of each bubble
        !! @param rvel Radial velocity of each bubble
        !! @param pos Position of each bubble
        !! @param posPrev Previous position of each bubble
        !! @param vel Velocity of each bubble
        !! @param scoord Cell index in real format of each bubble
        !! @param drad Radial velocity of each bubble
        !! @param drvel Radial acceleration of each bubble
        !! @param dgasp Time derivative of gas pressure in each bubble
        !! @param dgasmv Time derivative of vapor mass in each bubble
        !! @param dpos Time derivative of position of each bubble
        !! @param dvel Time derivative of velocity of each bubble
        !! @param lag_num_ts Number of stages in time-stepping scheme
        !! @param nBubs Local number of bubbles
    impure subroutine s_mpi_sendrecv_particles(bub_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, &
                                               gas_betaC, bub_dphidt, lag_id, gas_p, gas_mv, rad, &
                                               rvel, pos, posPrev, vel, scoord, drad, drvel, dgasp, &
                                               dgasmv, dpos, dvel, lag_num_ts, nbubs, dest)

        real(wp), dimension(:) :: bub_R0, Rmax_stats, Rmin_stats, gas_mg, gas_betaT, gas_betaC, bub_dphidt
        integer, dimension(:, :) :: lag_id
        real(wp), dimension(:, :) :: gas_p, gas_mv, rad, rvel, drad, drvel, dgasp, dgasmv
        real(wp), dimension(:, :, :) :: pos, posPrev, vel, scoord, dpos, dvel
        integer :: position, bub_id, lag_num_ts, tag, partner, send_tag, recv_tag, nbubs, p_recv_size, dest

        integer :: i, j, k, l, q, r
        integer :: req_send, req_recv, ierr !< Generic flag used to identify and report MPI errors
        integer :: send_count, send_offset, recv_count, recv_offset

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
            call MPI_Irecv(p_recv_counts(i, j, k), 1, MPI_INTEGER, partner, recv_tag, &
                           MPI_COMM_WORLD, recv_requests(recv_count), ierr)
        end do

        ! Post all sends
        do l = 1, n_neighbors
            i = neighbor_list(l, 1)
            j = neighbor_list(l, 2)
            k = neighbor_list(l, 3)
            partner = neighbor_ranks(i, j, k)
            send_tag = neighbor_tag(-i, -j, -k)

            send_count = send_count + 1
            call MPI_Isend(p_send_counts(i, j, k), 1, MPI_INTEGER, partner, send_tag, &
                           MPI_COMM_WORLD, send_requests(send_count), ierr)
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
                call MPI_Irecv(p_recv_buff(recv_offset), p_recv_size, MPI_PACKED, partner, recv_tag, &
                               MPI_COMM_WORLD, recv_requests(recv_count), ierr)
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

            if (p_send_counts(i, j, k) > 0) then
                partner = neighbor_ranks(i, j, k)
                send_tag = neighbor_tag(-i, -j, -k)

                ! Pack data for sending
                position = 0
                do q = 0, p_send_counts(i, j, k) - 1
                    bub_id = p_send_ids(i, j, k, q)
                    call MPI_Pack(lag_id(bub_id, 1), 1, MPI_INTEGER, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(bub_R0(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Rmax_stats(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(Rmin_stats(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(gas_mg(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(gas_betaT(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(gas_betaC(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    call MPI_Pack(bub_dphidt(bub_id), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    do r = 1, 2
                        call MPI_Pack(gas_p(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(gas_mv(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(rad(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(rvel(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(pos(bub_id, :, r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(posPrev(bub_id, :, r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(vel(bub_id, :, r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(scoord(bub_id, :, r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    end do
                    do r = 1, lag_num_ts
                        call MPI_Pack(drad(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(drvel(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(dgasp(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(dgasmv(bub_id, r), 1, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(dpos(bub_id, :, r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                        call MPI_Pack(dvel(bub_id, :, r), 3, mpi_p, p_send_buff(send_offset), p_buff_size, position, MPI_COMM_WORLD, ierr)
                    end do
                end do

                send_count = send_count + 1
                call MPI_Isend(p_send_buff(send_offset), position, MPI_PACKED, partner, send_tag, &
                               MPI_COMM_WORLD, send_requests(send_count), ierr)
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
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, lag_id(bub_id, 1), 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, bub_R0(bub_id), 1, mpi_p, MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, Rmax_stats(bub_id), 1, mpi_p, MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, Rmin_stats(bub_id), 1, mpi_p, MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_mg(bub_id), 1, mpi_p, MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_betaT(bub_id), 1, mpi_p, MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_betaC(bub_id), 1, mpi_p, MPI_COMM_WORLD, ierr)
                    call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, bub_dphidt(bub_id), 1, mpi_p, MPI_COMM_WORLD, ierr)
                    do r = 1, 2
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_p(bub_id, r), 1, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, gas_mv(bub_id, r), 1, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, rad(bub_id, r), 1, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, rvel(bub_id, r), 1, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, pos(bub_id, :, r), 3, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, posPrev(bub_id, :, r), 3, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, vel(bub_id, :, r), 3, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, scoord(bub_id, :, r), 3, mpi_p, MPI_COMM_WORLD, ierr)
                    end do
                    do r = 1, lag_num_ts
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, drad(bub_id, r), 1, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, drvel(bub_id, r), 1, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, dgasp(bub_id, r), 1, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, dgasmv(bub_id, r), 1, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, dpos(bub_id, :, r), 3, mpi_p, MPI_COMM_WORLD, ierr)
                        call MPI_Unpack(p_recv_buff(recv_offset), p_recv_size, position, dvel(bub_id, :, r), 3, mpi_p, MPI_COMM_WORLD, ierr)
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

    !! This function returns a unique tag for each neighbor based on its position
        !! relative to the current process.
        !! @param i, j, k Indices of the neighbor in the range [-1, 1]
        !! @return tag Unique integer tag for the neighbor
    integer function neighbor_tag(i, j, k) result(tag)

        integer, intent(in) :: i, j, k

        tag = (k + 1)*9 + (j + 1)*3 + (i + 1)

    end function neighbor_tag

    subroutine s_wrap_particle_positions(pos, posPrev, nbubs, dest)

        real(wp), dimension(:, :, :) :: pos, posPrev
        integer :: nbubs, dest
        integer :: i, q
        real :: offset

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

    impure subroutine s_mpi_send_random_number(phi_rn, num_freq)
        integer, intent(in) :: num_freq
        real(wp), intent(inout), dimension(1:num_freq) :: phi_rn

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors
        call MPI_BCAST(phi_rn, num_freq, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_send_random_number

    subroutine s_finalize_mpi_proxy_module()

#ifdef MFC_MPI
        if (ib) then
            @:DEALLOCATE(ib_buff_send, ib_buff_recv)
        end if
#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
