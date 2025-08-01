!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy

!> @brief This module serves as a proxy to the parameters and subroutines
!!              available in the MPI implementation's MPI module. Specifically,
!!              the role of the proxy is to harness basic MPI commands into more
!!              complex procedures as to achieve the required pre-processing
!!              communication goals.
module m_mpi_proxy

#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_helper

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_common

    implicit none

contains
    !> Since only processor with rank 0 is in charge of reading
            !!       and checking the consistency of the user provided inputs,
            !!       these are not available to the remaining processors. This
            !!       subroutine is then in charge of broadcasting the required
            !!       information.
    impure subroutine s_mpi_bcast_user_inputs

#ifdef MFC_MPI

        ! Generic loop iterator
        integer :: i
        ! Generic flag used to identify and report MPI errors
        integer :: ierr

        ! Logistics
        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in ['t_step_old', 't_step_start', 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',  &
            & 'loops_x', 'loops_y', 'loops_z', 'model_eqns', 'num_fluids',     &
            & 'weno_order', 'precision', 'perturb_flow_fluid', &
            & 'perturb_sph_fluid', 'num_patches', 'thermal', 'nb', 'dist_type',&
            & 'relax_model', 'num_ibs', 'n_start', 'elliptic_smoothing_iters', &
            & 'num_bc_patches', 'mixlayer_perturb_nk', 'recon_type', 'muscl_order']
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'old_grid','old_ic','stretch_x','stretch_y','stretch_z',&
            & 'cyl_coord','mpp_lim','hypoelasticity', 'relax', 'parallel_io',  &
            & 'perturb_flow', 'perturb_sph', 'mixlayer_vel_profile',           &
            & 'mixlayer_perturb', 'bubbles_euler', 'polytropic', 'polydisperse',&
            & 'qbmm', 'file_per_process', 'adv_n', 'ib' , 'cfl_adap_dt',       &
            & 'cfl_const_dt', 'cfl_dt', 'surface_tension',                     &
            & 'hyperelasticity', 'pre_stress', 'elliptic_smoothing', 'viscous',&
            & 'bubbles_lagrange', 'bc_io', 'mhd', 'relativity', 'cont_damage', &
            & 'igr' ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor
        call MPI_BCAST(fluid_rho(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in [ 'x_domain%beg', 'x_domain%end', 'y_domain%beg',         &
            & 'y_domain%end', 'z_domain%beg', 'z_domain%end', 'a_x', 'a_y',    &
            & 'a_z', 'x_a', 'x_b', 'y_a', 'y_b', 'z_a', 'z_b', 'bc_x%beg',     &
            & 'bc_x%end', 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end',      &
            & 'perturb_flow_mag', 'pref', 'rhoref', 'poly_sigma', 'R0ref',     &
            & 'Web', 'Ca', 'Re_inv', 'sigR', 'sigV', 'rhoRV', 'palpha_eps',    &
            & 'ptgalpha_eps', 'sigma', 'pi_fac', 'mixlayer_vel_coef', 'Bx0',   &
            & 'mixlayer_perturb_k0']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, num_bc_patches_max
            #:for VAR in ['geometry', 'type', 'dir', 'loc']
                call MPI_BCAST(patch_bc(i)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_bc(i)%radius, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in ['centroid', 'length']
                call MPI_BCAST(patch_bc(i)%${VAR}$, size(patch_bc(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        do i = 1, num_patches_max
            #:for VAR in [ 'geometry', 'smooth_patch_id']
                call MPI_BCAST(patch_icpp(i)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%smoothen, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%non_axis_sym, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alter_patch(0), num_patches_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'x_centroid', 'y_centroid', 'z_centroid',           &
                & 'length_x', 'length_y', 'length_z', 'radius', 'epsilon',     &
                & 'beta', 'smooth_coeff', 'rho', 'p0', 'm0', 'r0', 'v0',       &
                & 'pres', 'gamma', 'pi_inf', 'hcid', 'cv', 'qv', 'qvp',        &
                & 'model_threshold', 'cf_val', 'Bx', 'By', 'Bz']
                call MPI_BCAST(patch_icpp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ '2', '3', '4', '5', '6', '7', '8', '9']
                call MPI_BCAST(patch_icpp(i)%a(${VAR}$), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%model_filepath, len(patch_icpp(i)%model_filepath), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'model_translate', 'model_scale', 'model_rotate', &
                'normal', 'radii', 'vel', 'tau_e', 'alpha_rho', 'alpha' ]
                call MPI_BCAST(patch_icpp(i)%${VAR}$, size(patch_icpp(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%model_spc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            if (chemistry) then
                call MPI_BCAST(patch_icpp(i)%Y, size(patch_icpp(i)%Y), mpi_p, 0, MPI_COMM_WORLD, ierr)
            end if
            ! Broadcast IB variables
            call MPI_BCAST(patch_ib(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_filepath, len(patch_ib(i)%model_filepath), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_threshold, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_spc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'x_centroid', 'y_centroid', 'z_centroid',           &
                & 'length_x', 'length_y', 'length_z', 'radius', 'c', 'p', 't', 'm', 'theta']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%slip, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'model_translate', 'model_scale', 'model_rotate']
                call MPI_BCAST(patch_ib(i)%${VAR}$, size(patch_ib(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf','mul0','ss','pv','gamma_v','M_v',  &
                & 'mu_v','k_v', 'G', 'cv', 'qv', 'qvp' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do
#endif

    end subroutine s_mpi_bcast_user_inputs

end module m_mpi_proxy
