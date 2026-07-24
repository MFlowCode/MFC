!>
!! @file
!! @brief Contains module m_mpi_proxy

!> @brief Broadcasts user inputs and decomposes the domain across MPI ranks for pre-processing
module m_mpi_proxy

#ifdef MFC_MPI
    use mpi
#endif

    use m_helper
    use m_derived_types
    use m_global_parameters
    use m_mpi_common

    implicit none

contains
    !> Since only processor with rank 0 is in charge of reading and checking the consistency of the user provided inputs, these are
    !! not available to the remaining processors. This subroutine is then in charge of broadcasting the required information.
    impure subroutine s_mpi_bcast_user_inputs

#ifdef MFC_MPI
        integer :: i, j
        integer :: ierr

        ! Generated: case_dir, namelist scalars (INT/LOG/REAL), fluid_rho, fluid_pp loop, bub_pp
        #:include 'generated_bcast.fpp'

        ! manual: m/n/p_glb (computed from m/n/p post-read, not namelist-bound)
        call MPI_BCAST(m_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(n_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(p_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! manual: bc_x/y/z member broadcasts (struct members not in NAMELIST_VARS)
        #:for VAR in ['bc_x%isothermal_in', 'bc_y%isothermal_in', 'bc_z%isothermal_in',   &
            & 'bc_x%isothermal_out', 'bc_y%isothermal_out', 'bc_z%isothermal_out']
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        ! manual: domain bounds and wall temperatures (REAL struct members, not in NAMELIST_VARS)
        #:for VAR in [ 'x_domain%beg', 'x_domain%end', 'y_domain%beg',         &
            & 'y_domain%end', 'z_domain%beg', 'z_domain%end',                   &
            & 'bc_x%Twall_in', 'bc_x%Twall_out',                                &
            & 'bc_y%Twall_in', 'bc_y%Twall_out', 'bc_z%Twall_in',              &
            & 'bc_z%Twall_out']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        ! manual: BC type codes (INTEGER struct members)
        #:for VAR in [ 'bc_x%beg', 'bc_x%end', 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end']
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        ! wall-velocity members consumed by s_slip_wall/s_no_slip_wall on all ranks
        #:for DIM in ['x', 'y', 'z']
            #:for DIR in [1, 2, 3]
                call MPI_BCAST(bc_${DIM}$%vb${DIR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
                call MPI_BCAST(bc_${DIM}$%ve${DIR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        #:endfor

        ! manual: cfl_dt (runtime-computed logical), bc_io (BC-file existence)
        call MPI_BCAST(cfl_dt, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(bc_io, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        ! manual: patch_bc (complex array members broadcast with size())
        do i = 1, num_bc_patches_max
            #:for VAR in ['geometry', 'type', 'dir', 'loc']
                call MPI_BCAST(patch_bc(i)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_bc(i)%radius, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in ['centroid', 'length']
                call MPI_BCAST(patch_bc(i)%${VAR}$, size(patch_bc(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: patch_icpp (complex members: alter_patch, sph_har_coeff, size() arrays)
        do i = 1, num_patches_max
            #:for VAR in [ 'geometry', 'smooth_patch_id', 'hcid']
                call MPI_BCAST(patch_icpp(i)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%smoothen, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%non_axis_sym, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alter_patch(0), num_patches_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'x_centroid', 'y_centroid', 'z_centroid',           &
                & 'length_x', 'length_y', 'length_z', 'radius', 'epsilon',     &
                & 'beta', 'smooth_coeff', 'rho', 'p0', 'm0', 'r0', 'v0',       &
                & 'pres', 'gamma', 'pi_inf', 'cv', 'qv', 'qvp',        &
                & 'cf_val', 'Bx', 'By', 'Bz']
                call MPI_BCAST(patch_icpp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ '2', '3', '4', '5', '6', '7', '8', '9']
                call MPI_BCAST(patch_icpp(i)%a(${VAR}$), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'normal', 'radii', 'vel', 'tau_e', 'alpha_rho', 'alpha', &
                'fourier_cos', 'fourier_sin' ]
                call MPI_BCAST(patch_icpp(i)%${VAR}$, size(patch_icpp(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%sph_har_coeff, size(patch_icpp(i)%sph_har_coeff), mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%modal_clip_r_to_min, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%modal_r_min, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%modal_use_exp_form, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            call MPI_BCAST(patch_icpp(i)%model_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            if (chemistry) then
                call MPI_BCAST(patch_icpp(i)%Y, size(patch_icpp(i)%Y), mpi_p, 0, MPI_COMM_WORLD, ierr)
            end if
        end do

        ! manual: patch_ib (per-target member subset; pre uses size() for vel/angular_vel/angles)
        do i = 1, num_ibs
            #:for VAR in ['vel', 'angular_vel', 'angles']
                call MPI_BCAST(patch_ib(i)%${VAR}$, size(patch_ib(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'x_centroid', 'y_centroid', 'z_centroid', &
                & 'length_x', 'length_y', 'length_z', 'radius', 'twall']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%airfoil_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%slip, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%isothermal, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
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

        ! Hardcoded-patch variables (#1290; not namelist-bound, set on rank 0)
        call MPI_BCAST(interface_file, len(interface_file), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(normFac, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(normMag, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(g0_ic, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(p0_ic, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

        ! manual: simplex_params density perturbation (nested i-j loops; irregular structure)
        do i = 1, num_fluids_max
            call MPI_BCAST(simplex_params%perturb_dens(i), 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(simplex_params%perturb_dens_freq(i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(simplex_params%perturb_dens_scale(i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

            do j = 1, 3
                call MPI_BCAST(simplex_params%perturb_dens_offset(i, j), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            end do
        end do

        ! manual: simplex_params velocity perturbation
        do i = 1, 3
            call MPI_BCAST(simplex_params%perturb_vel(i), 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(simplex_params%perturb_vel_freq(i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(simplex_params%perturb_vel_scale(i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

            do j = 1, 3
                call MPI_BCAST(simplex_params%perturb_vel_offset(i, j), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            end do
        end do
#endif

    end subroutine s_mpi_bcast_user_inputs

end module m_mpi_proxy
