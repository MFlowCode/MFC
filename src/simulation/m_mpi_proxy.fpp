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
            & 'adap_dt_max_iters', 'int_comp', 'collision_model' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'run_time_info','cyl_coord', 'mpp_lim',     &
            &  'mp_weno', 'rdma_mpi', 'cont_damage', 'bc_io', &
            & 'weno_Re_flux', 'alt_soundspeed', 'null_weights', 'mixture_err',   &
            & 'parallel_io', 'hypoelasticity', 'bubbles_euler', 'polytropic',    &
            & 'polydisperse', 'qbmm', 'acoustic_source', 'probe_wrt', 'integral_wrt',   &
            & 'prim_vars_wrt', 'weno_avg', 'file_per_process', 'relax',          &
            & 'adv_n', 'adap_dt', 'ib', 'bodyForces', 'bf_x', 'bf_y', 'bf_z',    &
            & 'bc_x%grcbc_in', 'bc_x%grcbc_out', 'bc_x%grcbc_vel_out',          &
            & 'bc_y%grcbc_in', 'bc_y%grcbc_out', 'bc_y%grcbc_vel_out',          &
            & 'bc_z%grcbc_in', 'bc_z%grcbc_out', 'bc_z%grcbc_vel_out',          &
            & 'bc_x%isothermal_in', 'bc_y%isothermal_in', 'bc_z%isothermal_in',        &
            & 'bc_x%isothermal_out', 'bc_y%isothermal_out', 'bc_z%isothermal_out', &
            & 'cfl_adap_dt', 'cfl_const_dt', 'cfl_dt', 'surface_tension',       &
            & 'shear_stress', 'bulk_stress', 'bubbles_lagrange',                &
            & 'hyperelasticity', 'down_sample', 'fft_wrt', &
            & 'hyper_cleaning', 'ib_state_wrt']
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        if (chemistry) then
            #:for VAR in [ 'diffusion', 'reactions' ]
                call MPI_BCAST(chem_params%${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'gamma_method', 'transport_model' ]
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
                & 'drag_model', 'charNz']
                call MPI_BCAST(lag_params%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in ['epsilonb','charwidth','valmaxvoid']
                call MPI_BCAST(lag_params%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(lag_params%input_path, len(lag_params%input_path), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        end if

        #:for VAR in [ 'dt','weno_eps','teno_CT','pref','rhoref','R0ref','Web','Ca', 'sigma', &
            & 'Re_inv', 'poly_sigma', 'palpha_eps', 'ptgalpha_eps', 'pi_fac',    &
            & 'bc_x%vb1','bc_x%vb2','bc_x%vb3','bc_x%ve1','bc_x%ve2','bc_x%ve3', &
            & 'bc_y%vb1','bc_y%vb2','bc_y%vb3','bc_y%ve1','bc_y%ve2','bc_y%ve3', &
            & 'bc_z%vb1','bc_z%vb2','bc_z%vb3','bc_z%ve1','bc_z%ve2','bc_z%ve3', &
            & 'bc_x%pres_in','bc_x%pres_out','bc_y%pres_in','bc_y%pres_out', 'bc_z%pres_in','bc_z%pres_out', &
            & 'bc_x%Twall_in', 'bc_x%Twall_out', 'bc_y%Twall_in', 'bc_y%Twall_out',  &
            & 'bc_z%Twall_in', 'bc_z%Twall_out', &
            & 't_stop', 't_save', 'cfl_target', 'Bx0', 'alf_factor',  &
            & 'tau_star', 'cont_damage_s', 'alpha_bar', 'adap_dt_tol', &
            & 'ic_eps', 'ic_beta', 'hyper_cleaning_speed', &
            & 'hyper_cleaning_tau', 'coefficient_of_restitution', 'collision_time', &
            & 'ib_coefficient_of_friction' ]
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, 3
            #:for VAR in [ 'bc_x%vel_in', 'bc_x%vel_out', 'bc_y%vel_in', 'bc_y%vel_out',  &
                & 'bc_z%vel_in', 'bc_z%vel_out' ]
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
            #:for VAR in [ 'gamma','pi_inf','G','cv','qv','qvp' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(fluid_pp(i)%Re(1), 2, mpi_p, 0, MPI_COMM_WORLD, ierr)
        end do

        if (bubbles_euler .or. bubbles_lagrange) then
            #:for VAR in [ 'R0ref','p0ref','rho0ref','T0ref', &
                'ss','pv','vd','mu_l','mu_v','mu_g','gam_v','gam_g',&
                'M_v','M_g','k_v','k_g','cp_v','cp_g','R_v','R_g']
                call MPI_BCAST(bub_pp%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end if

        do i = 1, num_fluids_max
            #:for VAR in ['bc_x%alpha_rho_in','bc_x%alpha_in','bc_y%alpha_rho_in','bc_y%alpha_in','bc_z%alpha_rho_in', &
                & 'bc_z%alpha_in']
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        do i = 1, num_ibs
            #:for VAR in [ 'radius', 'length_x', 'length_y', 'length_z', &
                & 'x_centroid', 'y_centroid', 'z_centroid', 'c', 'm', 'p', 't', 'theta', 'slip', 'mass', &
                & 'model_threshold']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            #:for VAR in ['vel', 'angular_vel', 'angles', 'model_translate', 'model_scale']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%moving_ibm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_spc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_filepath, len(patch_ib(i)%model_filepath), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
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

        ! Particle MPI buffers are only allocated when num_procs > 1
        if (allocated(p_send_buff)) then
            @:DEALLOCATE(p_send_buff, p_recv_buff)
            @:DEALLOCATE(p_send_ids)
        end if
#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
