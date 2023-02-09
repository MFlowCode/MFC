# 1 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy

!> @brief This module serves as a proxy to the parameters and subroutines
!!              available in the MPI implementation's MPI module. Specifically,
!!              the role of the proxy is to harness basic MPI commands into more
!!              complex procedures as to achieve the required pre-processing
!!              communication goals.
module m_mpi_proxy

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_common
    ! ==========================================================================

    implicit none

    integer, private :: err_code, ierr !<
    !! Generic flags used to identify and report MPI errors

contains


    !> Since only processor with rank 0 is in charge of reading
        !!       and checking the consistency of the user provided inputs,
        !!       these are not available to the remaining processors. This
        !!       subroutine is then in charge of broadcasting the required
        !!       information.
    subroutine s_mpi_bcast_user_inputs() ! ---------------------------------

#ifdef MFC_MPI

        ! Generic loop iterator
        integer :: i

        ! Logistics
        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(t_step_old, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(p, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(m_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(n_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(p_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(loops_x, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(loops_y, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(loops_z, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(model_eqns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(num_fluids, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(precision, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(perturb_flow_fluid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(perturb_sph_fluid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(num_patches, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(thermal, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(dist_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 52 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(R0_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
# 54 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"

# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(old_grid, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(old_ic, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(stretch_x, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(stretch_y, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(stretch_z, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(cyl_coord, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(adv_alphan, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(mpp_lim, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(hypoelasticity, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(parallel_io, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(perturb_flow, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(perturb_sph, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(bubbles, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(polytropic, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(polydisperse, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 59 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(qbmm, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
# 61 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
        call MPI_BCAST(fluid_rho(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)


# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(x_domain%beg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(x_domain%end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(y_domain%beg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(y_domain%end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(z_domain%beg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(z_domain%end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(a_x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(a_y, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(a_z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(x_a, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(x_b, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(y_a, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(y_b, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(z_a, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(z_b, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(bc_x%beg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(bc_x%end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(bc_y%beg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(bc_y%end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(bc_z%beg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(bc_z%end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(pref, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(rhoref, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(poly_sigma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(R0ref, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(Web, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(Ca, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(Re_inv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(sigR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(sigV, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 70 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(rhoRV, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 72 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"

        do i = 1, num_patches_max
            call MPI_BCAST(patch_icpp(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%smooth_patch_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            call MPI_BCAST(patch_icpp(i)%smoothen, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alter_patch(0), num_patches_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%x_centroid, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%y_centroid, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%z_centroid, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%length_x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%length_y, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%length_z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%radius, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%epsilon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%beta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%smooth_coeff, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%rho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%p0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%m0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%r0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%v0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%pres, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%gamma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 84 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(patch_icpp(i)%pi_inf, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 86 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
            call MPI_BCAST(patch_icpp(i)%normal(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%radii(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%vel(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%tau_e(1), 6, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alpha_rho(1), num_fluids_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alpha(1), num_fluids_max - 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        end do

        ! Fluids physical parameters
        do i = 1, num_fluids_max
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%gamma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%pi_inf, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%mul0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%ss, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%pv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%gamma_v, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%M_v, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%mu_v, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%k_v, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 98 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
                call MPI_BCAST(fluid_pp(i)%G, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/pre_process/m_mpi_proxy.fpp"
        end do
#endif

    end subroutine s_mpi_bcast_user_inputs ! -------------------------------


    !> Description: This subroutine takes care of efficiently distributing
        !!              the computational domain among the available processors
        !!             as well as recomputing some of the global parameters so
        !!              that they reflect the configuration of sub-domain that is
        !!              overseen by the local processor.
    subroutine s_mpi_decompose_computational_domain() ! --------------------

#ifdef MFC_MPI

        ! # of processors in the x-, y- and z-coordinate directions
        integer :: num_procs_x, num_procs_y, num_procs_z

        ! Temporary # of processors in x-, y- and z-coordinate directions
        ! used during the processor factorization optimization procedure
        real(kind(0d0)) :: tmp_num_procs_x, tmp_num_procs_y, tmp_num_procs_z

        ! Processor factorization (fct) minimization parameter
        real(kind(0d0)) :: fct_min

        ! Cartesian processor topology communicator
        integer :: MPI_COMM_CART

        ! Number of remaining cells for a particular coordinate direction
        ! after the bulk has evenly been distributed among the available
        ! processors for that coordinate direction
        integer :: rem_cells

        ! Generic loop iterators
        integer :: i, j

        if (num_procs == 1 .and. parallel_io) then
            do i = 1, num_dims
                start_idx(i) = 0
            end do
            return
        end if

        ! Performing the computational domain decomposition. The procedure
        ! is optimized by ensuring that each processor contains a close to
        ! equivalent piece of the computational domain. Note that explicit
        ! type-casting is omitted here for code legibility purposes.

        ! Generating 3D Cartesian Processor Topology =======================

        if (n > 0) then

            if (p > 0) then

                if (cyl_coord .and. p > 0) then
                    ! Implement pencil processor blocking if using cylindrical coordinates so
                    ! that all cells in azimuthal direction are stored on a single processor.
                    ! This is necessary for efficient application of Fourier filter near axis.

                    ! Initial values of the processor factorization optimization
                    num_procs_x = 1
                    num_procs_y = num_procs
                    num_procs_z = 1
                    ierr = -1

                    ! Computing minimization variable for these initial values
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    tmp_num_procs_z = num_procs_z
                    fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                                       - (n + 1)/tmp_num_procs_y)

                    ! Searching for optimal computational domain distribution
                    do i = 1, num_procs

                        if (mod(num_procs, i) == 0 &
                            .and. &
                            (m + 1)/i >= num_stcls_min*weno_order) then

                            tmp_num_procs_x = i
                            tmp_num_procs_y = num_procs/i

                            if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                               - (n + 1)/tmp_num_procs_y) &
                                .and. &
                                (n + 1)/tmp_num_procs_y &
                                >= &
                                num_stcls_min*weno_order) then

                                num_procs_x = i
                                num_procs_y = num_procs/i
                                fct_min = abs((m + 1)/tmp_num_procs_x &
                                              - (n + 1)/tmp_num_procs_y)
                                ierr = 0

                            end if

                        end if

                    end do

                else

                    ! Initial values of the processor factorization optimization
                    num_procs_x = 1
                    num_procs_y = 1
                    num_procs_z = num_procs
                    ierr = -1

                    ! Computing minimization variable for these initial values
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    tmp_num_procs_z = num_procs_z
                    fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                                       - (n + 1)/tmp_num_procs_y) &
                              + 10d0*abs((n + 1)/tmp_num_procs_y &
                                         - (p + 1)/tmp_num_procs_z)

                    ! Searching for optimal computational domain distribution
                    do i = 1, num_procs

                        if (mod(num_procs, i) == 0 &
                            .and. &
                            (m + 1)/i >= num_stcls_min*weno_order) then

                            do j = 1, (num_procs/i)

                                if (mod(num_procs/i, j) == 0 &
                                    .and. &
                                    (n + 1)/j >= num_stcls_min*weno_order) then

                                    tmp_num_procs_x = i
                                    tmp_num_procs_y = j
                                    tmp_num_procs_z = num_procs/(i*j)

                                    if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                                       - (n + 1)/tmp_num_procs_y) &
                                        + abs((n + 1)/tmp_num_procs_y &
                                              - (p + 1)/tmp_num_procs_z) &
                                        .and. &
                                        (p + 1)/tmp_num_procs_z &
                                        >= &
                                        num_stcls_min*weno_order) &
                                        then

                                        num_procs_x = i
                                        num_procs_y = j
                                        num_procs_z = num_procs/(i*j)
                                        fct_min = abs((m + 1)/tmp_num_procs_x &
                                                      - (n + 1)/tmp_num_procs_y) &
                                                  + abs((n + 1)/tmp_num_procs_y &
                                                        - (p + 1)/tmp_num_procs_z)
                                        ierr = 0

                                    end if

                                end if

                            end do

                        end if

                    end do

                end if

                ! Checking whether the decomposition of the computational
                ! domain was successful
                if (proc_rank == 0 .and. ierr == -1) then
                    print '(A)', 'Unable to decompose computational '// &
                        'domain for selected number of '// &
                        'processors. Exiting ...'
                    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                end if

                ! Creating a new communicator using Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 3, (/num_procs_x, &
                                                          num_procs_y, num_procs_z/), &
                                     (/.true., .true., .true./), &
                                     .false., MPI_COMM_CART, ierr)

                ! Finding corresponding Cartesian coordinates of the local
                ! processor rank in newly declared cartesian communicator
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 3, &
                                     proc_coords, ierr)

                ! END: Generating 3D Cartesian Processor Topology ==================

                ! Sub-domain Global Parameters in z-direction ======================

                ! Number of remaining cells after majority is distributed
                rem_cells = mod(p + 1, num_procs_z)

                ! Preliminary uniform cell-width spacing
                if (old_grid .neqv. .true.) then
                    dz = (z_domain%end - z_domain%beg)/real(p + 1, kind(0d0))
                end if

                ! Optimal number of cells per processor
                p = (p + 1)/num_procs_z - 1

                ! Distributing any remaining cells
                do i = 1, rem_cells
                    if (proc_coords(3) == i - 1) then
                        p = p + 1
                        exit
                    end if
                end do

                ! Beginning and end sub-domain boundary locations
                if (parallel_io .neqv. .true.) then
                    if (old_grid .neqv. .true.) then
                        if (proc_coords(3) < rem_cells) then
                            z_domain%beg = z_domain%beg + dz*real((p + 1)* &
                                                                  proc_coords(3))
                            z_domain%end = z_domain%end - dz*real((p + 1)* &
                                                                  (num_procs_z - proc_coords(3) - 1) &
                                                                  - (num_procs_z - rem_cells))
                        else
                            z_domain%beg = z_domain%beg + dz*real((p + 1)* &
                                                                  proc_coords(3) + rem_cells)
                            z_domain%end = z_domain%end - dz*real((p + 1)* &
                                                                  (num_procs_z - proc_coords(3) - 1))
                        end if
                    end if
                else
                    if (proc_coords(3) < rem_cells) then
                        start_idx(3) = (p + 1)*proc_coords(3)
                    else
                        start_idx(3) = (p + 1)*proc_coords(3) + rem_cells
                    end if
                end if

                ! ==================================================================

                ! Generating 2D Cartesian Processor Topology =======================

            else

                ! Initial values of the processor factorization optimization
                num_procs_x = 1
                num_procs_y = num_procs
                ierr = -1

                ! Computing minimization variable for these initial values
                tmp_num_procs_x = num_procs_x
                tmp_num_procs_y = num_procs_y
                fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                                   - (n + 1)/tmp_num_procs_y)

                ! Searching for optimal computational domain distribution
                do i = 1, num_procs

                    if (mod(num_procs, i) == 0 &
                        .and. &
                        (m + 1)/i >= num_stcls_min*weno_order) then

                        tmp_num_procs_x = i
                        tmp_num_procs_y = num_procs/i

                        if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                           - (n + 1)/tmp_num_procs_y) &
                            .and. &
                            (n + 1)/tmp_num_procs_y &
                            >= &
                            num_stcls_min*weno_order) then

                            num_procs_x = i
                            num_procs_y = num_procs/i
                            fct_min = abs((m + 1)/tmp_num_procs_x &
                                          - (n + 1)/tmp_num_procs_y)
                            ierr = 0

                        end if

                    end if

                end do

                ! Checking whether the decomposition of the computational
                ! domain was successful
                if (proc_rank == 0 .and. ierr == -1) then
                    print '(A)', 'Unable to decompose computational '// &
                        'domain for selected number of '// &
                        'processors. Exiting ...'
                    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                end if

                ! Creating a new communicator using Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 2, (/num_procs_x, &
                                                          num_procs_y/), (/.true., &
                                                                           .true./), .false., MPI_COMM_CART, &
                                     ierr)

                ! Finding corresponding Cartesian coordinates of the local
                ! processor rank in newly declared cartesian communicator
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, &
                                     proc_coords, ierr)

            end if

            ! END: Generating 2D Cartesian Processor Topology ==================

            ! Sub-domain Global Parameters in y-direction ======================

            ! Number of remaining cells after majority has been distributed
            rem_cells = mod(n + 1, num_procs_y)

            ! Preliminary uniform cell-width spacing
            if (old_grid .neqv. .true.) then
                dy = (y_domain%end - y_domain%beg)/real(n + 1, kind(0d0))
            end if

            ! Optimal number of cells per processor
            n = (n + 1)/num_procs_y - 1

            ! Distributing any remaining cells
            do i = 1, rem_cells
                if (proc_coords(2) == i - 1) then
                    n = n + 1
                    exit
                end if
            end do

            ! Beginning and end sub-domain boundary locations
            if (parallel_io .neqv. .true.) then
                if (old_grid .neqv. .true.) then
                    if (proc_coords(2) < rem_cells) then
                        y_domain%beg = y_domain%beg + dy*real((n + 1)* &
                                                              proc_coords(2))
                        y_domain%end = y_domain%end - dy*real((n + 1)* &
                                                              (num_procs_y - proc_coords(2) - 1) &
                                                              - (num_procs_y - rem_cells))
                    else
                        y_domain%beg = y_domain%beg + dy*real((n + 1)* &
                                                              proc_coords(2) + rem_cells)
                        y_domain%end = y_domain%end - dy*real((n + 1)* &
                                                              (num_procs_y - proc_coords(2) - 1))
                    end if
                end if
            else
                if (proc_coords(2) < rem_cells) then
                    start_idx(2) = (n + 1)*proc_coords(2)
                else
                    start_idx(2) = (n + 1)*proc_coords(2) + rem_cells
                end if
            end if

            ! ==================================================================

            ! Generating 1D Cartesian Processor Topology =======================

        else

            ! Number of processors in the coordinate direction is equal to
            ! the total number of processors available
            num_procs_x = num_procs

            ! Creating a new communicator using Cartesian topology
            call MPI_CART_CREATE(MPI_COMM_WORLD, 1, (/num_procs_x/), &
                                 (/.true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding the corresponding Cartesian coordinates of the local
            ! processor rank in the newly declared cartesian communicator
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 1, &
                                 proc_coords, ierr)

        end if

        ! ==================================================================

        ! Sub-domain Global Parameters in x-direction ======================

        ! Number of remaining cells after majority has been distributed
        rem_cells = mod(m + 1, num_procs_x)

        ! Preliminary uniform cell-width spacing
        if (old_grid .neqv. .true.) then
            dx = (x_domain%end - x_domain%beg)/real(m + 1, kind(0d0))
        end if

        ! Optimal number of cells per processor
        m = (m + 1)/num_procs_x - 1

        ! Distributing any remaining cells
        do i = 1, rem_cells
            if (proc_coords(1) == i - 1) then
                m = m + 1
                exit
            end if
        end do

        ! Beginning and end sub-domain boundary locations
        if (parallel_io .neqv. .true.) then
            if (old_grid .neqv. .true.) then
                if (proc_coords(1) < rem_cells) then
                    x_domain%beg = x_domain%beg + dx*real((m + 1)* &
                                                          proc_coords(1))
                    x_domain%end = x_domain%end - dx*real((m + 1)* &
                                                          (num_procs_x - proc_coords(1) - 1) &
                                                          - (num_procs_x - rem_cells))
                else
                    x_domain%beg = x_domain%beg + dx*real((m + 1)* &
                                                          proc_coords(1) + rem_cells)
                    x_domain%end = x_domain%end - dx*real((m + 1)* &
                                                          (num_procs_x - proc_coords(1) - 1))
                end if
            end if
        else
            if (proc_coords(1) < rem_cells) then
                start_idx(1) = (m + 1)*proc_coords(1)
            else
                start_idx(1) = (m + 1)*proc_coords(1) + rem_cells
            end if
        end if

        ! ==================================================================

#endif

    end subroutine s_mpi_decompose_computational_domain ! ------------------

end module m_mpi_proxy
