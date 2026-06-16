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
            & 'x_domain%beg', 'x_domain%end', 'y_domain%beg', 'y_domain%end',              &
            & 'z_domain%beg', 'z_domain%end',                                               &
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
                & 'x_centroid', 'y_centroid', 'z_centroid', 'slip', 'mass', 'Twall']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            #:for VAR in ['vel', 'angular_vel', 'angles']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%isothermal, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%moving_ibm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%airfoil_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
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
#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
