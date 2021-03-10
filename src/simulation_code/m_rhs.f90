!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /
!!    / /  / / __/ / /___
!!   /_/  /_/_/    \____/
!!
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names
!!  are listed in the copyright file included with this source
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published
!!  by the Free Software Foundation, either version 3 of the license
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file m_rhs.f90
!! @brief Contains module m_rhs
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The module contains the subroutines used to calculate the right-
!!              hand-side (RHS) in the quasi-conservative, shock- and interface-
!!              capturing finite-volume framework for the multicomponent Navier-
!!              Stokes equations supplemented by appropriate advection equations
!!              used to capture the material interfaces. The system of equations
!!              is closed by the stiffened gas equation of state, as well as any
!!              required mixture relationships. Capillarity effects are included
!!              and are modeled by the means of a volume force acting across the
!!              diffuse material interface region. The implementation details of
!!              surface tension may be found in Perigaud and Saurel (2005). Note
!!              that both viscous and surface tension effects are only available
!!              in the volume fraction model.
module m_rhs

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables

    use m_riemann_solvers      !< Exact and approximate Riemann problem solvers

    use m_cbc                  !< Characteristic boundary conditions (CBC)

    use m_bubbles              !< Bubble dynamic routines

    use m_qbmm                 !< Moment inversion
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_rhs_module, &
 s_compute_rhs, &
 s_pressure_relaxation_procedure, &
 s_populate_variables_buffers, &
 s_finalize_rhs_module, &
 s_get_crv, &
 s_get_viscous

    type(vector_field), allocatable, dimension(:, :, :) :: q_cons_qp !<
    !! This variable contains the WENO-reconstructed values of the cell-average
    !! conservative variables, which are located in q_cons_vf, at cell-interior
    !! Gaussian quadrature points (QP).

    type(vector_field), allocatable, dimension(:, :, :) :: q_prim_qp !<
    !! The primitive variables at cell-interior Gaussian quadrature points. These
    !! are calculated from the conservative variables and gradient magnitude (GM)
    !! of the volume fractions, q_cons_qp and gm_alpha_qp, respectively.

    !> @name The left (L) and the right (R) WENO-reconstructed cell-boundary values,
    !! including cell-boundary Gaussian quadrature points, of the cell-average
    !! conservative variables. The latter are stored in the variable q_cons_qp
    !! (NDQP - normal direction quadrature points).
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: qL_cons_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: qR_cons_ndqp
    !> @}

    !> @name The left and right WENO-reconstructed cell-boundary values, that include
    !! cell-boundary Gaussian quadrature points, of the cell-averaged primitive
    !! variables. The latter are stored in the variable q_prim_qp.
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: qL_prim_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: qR_prim_ndqp
    !> @}

    !> @name The first-order spatial derivatives of the primitive variables at cell-
    !! interior Guassian quadrature points. These are WENO-reconstructed from
    !! their respective cell-average values, obtained through the application
    !! of the divergence theorem on the integral-average cell-boundary values
    !! of the primitive variables, located in qK_prim_ndqp, where K = L or R.
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: dq_prim_dx_qp
    type(vector_field), allocatable, dimension(:, :, :) :: dq_prim_dy_qp
    type(vector_field), allocatable, dimension(:, :, :) :: dq_prim_dz_qp
    type(vector_field), allocatable, dimension(:, :, :) :: gm_vel_qp
    !> @}

    !> @name The left and right WENO-reconstructed cell-boundary values of the cell-
    !! average first-order spatial derivatives of the primitive variables. The
    !! cell-average of the first-order spatial derivatives may be found in the
    !! variables dq_prim_ds_qp, where s = x, y or z.
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: dqL_prim_dx_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: dqL_prim_dy_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: dqL_prim_dz_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: dqR_prim_dx_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: dqR_prim_dy_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: dqR_prim_dz_ndqp
    !> @}

    type(vector_field), allocatable, dimension(:, :, :) :: gm_alpha_qp  !<
    !! The gradient magnitude of the volume fractions at cell-interior Gaussian
    !! quadrature points. gm_alpha_qp is calculated from individual first-order
    !! spatial derivatives located in dq_prim_ds_qp.

    !> @name The left and right WENO-reconstructed cell-boundary values of the cell-
    !! average gradient magnitude of volume fractions, located in gm_alpha_qp.
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: gm_alphaL_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: gm_alphaR_ndqp
    !> @}

    !> @name The left and right cell-boundary values of the vector components of the
    !! unit normal (UN) of the volume fractions. These are directly calculated
    !! from the first-order spatial derivatives, stored in dqK_prim_ds_ndqp.
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: un_alphaL_x_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: un_alphaL_y_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: un_alphaL_z_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: un_alphaR_x_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: un_alphaR_y_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: un_alphaR_z_ndqp
    !> @}

    type(scalar_field), allocatable, dimension(:) :: kappa_vf !<
    !! The cell-average curvature of volume fractions. kappa_vf is calculated by
    !! the divergence theorem using the integral-average cell-boundary values of
    !! of the components of the unit normals, which are located in the variables
    !! un_alphaK_s_ndqp.

    !> @name The left and right WENO-reconstructed cell-boundary values of the cell-
    !! average curvature of the volume fractions. The cell-average quantities
    !! are located in the variable kappa_vf.
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: kappaL_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: kappaR_ndqp
    !> @}

    !> @name The cell-boundary values of the fluxes (src - source, gsrc - geometrical
    !! source). These are computed by applying the chosen Riemann problem solver
    !! on the left and right cell-boundary values of the primitive variables,
    !! qK_prim_ndqp, the first-order spatial derivatives, dqK_prim_ds_ndqp, as
    !! well as the curvature of volume fractions, kappaK_ndqp.
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: flux_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: flux_src_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: flux_gsrc_ndqp
    !> @}

    !> @name Additional vector and scalar fields needed in alternate computation
    !! of the cell-average curvature of volume fractions.
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: dgm_alpha_dx_qp
    type(vector_field), allocatable, dimension(:, :, :) :: dgm_alpha_dy_qp
    type(vector_field), allocatable, dimension(:, :, :) :: dgm_alpha_dz_qp
    type(vector_field), allocatable, dimension(:, :, :) :: un_alpha_x_qp
    type(vector_field), allocatable, dimension(:, :, :) :: un_alpha_y_qp
    type(vector_field), allocatable, dimension(:, :, :) :: un_alpha_z_qp

    type(scalar_field), allocatable, dimension(:) :: laplacian_vf
    type(scalar_field), allocatable, dimension(:) :: alt_kappa_vf
    !> @}

    type(scalar_field), allocatable, dimension(:) :: reg_src_vf !<
    !! Additional field for regularization terms

    !> @name Additional field for capillary source terms
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :) :: We_mtm_src
    real(kind(0d0)), allocatable, dimension(:, :, :) :: We_nrg_src
    type(scalar_field), allocatable, dimension(:) :: tau_We_vf
    type(scalar_field), allocatable, dimension(:) :: tau_Re_vf
    !> @}

    !> @name Additional variables for applying a flux limiter to the advection equation
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :) :: flux_lim_func
    type(vector_field), allocatable, dimension(:, :, :) :: lo_flux_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: lo_flux_src_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: lo_flux_gsrc_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: hi_flux_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: hi_flux_src_ndqp
    type(vector_field), allocatable, dimension(:, :, :) :: hi_flux_gsrc_ndqp
    !> @}

    !> @name The indical bounds in the coordinate directions of the Gaussian numerical
    !! integration that identify the locations of the quadrature points that are
    !! necessary to evaluate cell-interior integral-averages.
    !> @{
    type(bounds_info) :: ieta, iksi, itau
    !> @}

    !> @name The indical bounds in the coordinate directions of the Gaussian numerical
    !! integration that identify the locations of the quadrature points that are
    !! necessary to evaluate cell-boundary integral-averages.
    !> @{
    type(bounds_info) :: ichi, ipsi
    !> @}

    type(bounds_info) :: iv !< Vector field indical bounds

    !> @name Indical bounds in the x-, y- and z-directions
    !> @{
    type(bounds_info) :: ix, iy, iz
    !> @}

    !> @name Bubble dynamic source terms
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :) :: bub_adv_src
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: bub_r_src, bub_v_src, bub_p_src, bub_m_src
    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: bub_mom_src
    ! REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:,:) :: mom_sp
    ! REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: mom_3d

    type(scalar_field) :: divu !< matrix for div(u)
    !> @}

    !> @name Monopole source terms
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mono_mass_src, mono_e_src
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: mono_mom_src
    !> @}

    !> @name Saved fluxes for testing
    !> @{
    type(vector_field), allocatable, dimension(:, :, :) :: myflux_vf, myflux_src_vf
    type(scalar_field) :: alf_sum
    !> @}

    character(50) :: file_path !< Local file path for saving debug files

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_rhs_module() ! ---------------------------------

        integer :: i, j, k, l !< Generic loop iterators

        ! Configuring Cell-Interior Quadrature Points ======================
        ieta%beg = 0; iksi%beg = 0; itau%beg = 0

        if (commute_err) then
            ieta%beg = -1; if (n > 0) iksi%beg = -1; if (p > 0) itau%beg = -1
        end if

        ieta%end = -ieta%beg; iksi%end = -iksi%beg; itau%end = -itau%beg
        ! ==================================================================

        ! Configuring Cell-Boundary Quadrature Points ======================
        ichi%beg = 0; ipsi%beg = 0

        if (split_err) then
            ichi%beg = -1; if (p > 0) ipsi%beg = -1
        end if

        ichi%end = -ichi%beg; ipsi%end = -ipsi%beg
        ! ==================================================================

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================

        if (We_size > 0 .and. We_rhs_flux) then
            allocate (tau_We_vf(1:sys_size))
            do i = 1, num_dims
                allocate (tau_We_vf(cont_idx%end + i)%sf(ix%beg:ix%end, &
                                                         iy%beg:iy%end, &
                                                         iz%beg:iz%end))
            end do
            allocate (tau_We_vf(E_idx)%sf(ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
        end if

        if (any(Re_size > 0) .and. cyl_coord) then
            allocate (tau_Re_vf(1:sys_size))
            do i = 1, num_dims
                allocate (tau_Re_vf(cont_idx%end + i)%sf(ix%beg:ix%end, &
                                                         iy%beg:iy%end, &
                                                         iz%beg:iz%end))
            end do
            allocate (tau_Re_vf(E_idx)%sf(ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
        end if

        ! Allocation/Association of q_cons_qp and q_prim_qp ================
        allocate (q_cons_qp(ieta%beg:ieta%end, &
                            iksi%beg:iksi%end, &
                            itau%beg:itau%end))
        allocate (q_prim_qp(ieta%beg:ieta%end, &
                            iksi%beg:iksi%end, &
                            itau%beg:itau%end))

        do k = itau%beg, itau%end
            do j = iksi%beg, iksi%end
                do i = ieta%beg, ieta%end

                    allocate (q_cons_qp(i, j, k)%vf(1:sys_size))
                    allocate (q_prim_qp(i, j, k)%vf(1:sys_size))

                    if (i /= 0 .and. abs(j) >= abs(k)) then

                        do l = 1, sys_size
                            allocate (q_cons_qp(i, j, k)%vf(l)%sf(ix%beg:ix%end, &
                                                                  iy%beg:iy%end, &
                                                                  iz%beg:iz%end))
                            q_prim_qp(i, j, k)%vf(l)%sf => &
                                q_cons_qp(i, j, k)%vf(l)%sf
                        end do

                        do l = mom_idx%beg, E_idx
                            allocate (q_prim_qp(i, j, k)%vf(l)%sf(ix%beg:ix%end, &
                                                                  iy%beg:iy%end, &
                                                                  iz%beg:iz%end))
                        end do
                        if (model_eqns == 3) then
                            do l = internalEnergies_idx%beg, internalEnergies_idx%end
                                allocate (q_prim_qp(i, j, k)%vf(l)%sf(ix%beg:ix%end, &
                                                                      iy%beg:iy%end, &
                                                                      iz%beg:iz%end))
                            end do
                        end if

                        if (bubbles) then
                            do l = bub_idx%beg, bub_idx%end
                                allocate (q_prim_qp(i, j, k)%vf(l)%sf(ix%beg:ix%end, &
                                                                      iy%beg:iy%end, &
                                                                      iz%beg:iz%end))
                            end do
                        end if

                        if (hypoelasticity) then
                            do l = stress_idx%beg, stress_idx%end
                                allocate (q_prim_qp(i, j, k)%vf(l)%sf(ix%beg:ix%end, &
                                                                      iy%beg:iy%end, &
                                                                      iz%beg:iz%end))
                            end do
                        end if
                    end if

                end do
            end do
        end do
        ! ==================================================================

        if (qbmm) then
            allocate (mom_sp(1:nmomsp), mom_3d(0:2, 0:2, nb))
            do i = 0, 2; do j = 0, 2; do k = 1, nb
                    allocate (mom_3d(i, j, k)%sf( &
                              ix%beg:ix%end, &
                              iy%beg:iy%end, &
                              iz%beg:iz%end))
                end do; end do; end do
            do i = 1, nmomsp
                allocate (mom_sp(i)%sf( &
                          ix%beg:ix%end, &
                          iy%beg:iy%end, &
                          iz%beg:iz%end))
            end do
        end if

        ! Allocation/Association of qK_cons_ndqp and qK_prim_ndqp ==========
        allocate (qL_cons_ndqp(1:num_dims, &
                               ichi%beg:ichi%end, &
                               ipsi%beg:ipsi%end))
        allocate (qR_cons_ndqp(1:num_dims, &
                               ichi%beg:ichi%end, &
                               ipsi%beg:ipsi%end))
        allocate (qL_prim_ndqp(1:num_dims, &
                               ichi%beg:ichi%end, &
                               ipsi%beg:ipsi%end))
        allocate (qR_prim_ndqp(1:num_dims, &
                               ichi%beg:ichi%end, &
                               ipsi%beg:ipsi%end))

        allocate (myflux_vf(1:num_dims, &
                            ichi%beg:ichi%end, &
                            ipsi%beg:ipsi%end))
        allocate (myflux_src_vf(1:num_dims, &
                                ichi%beg:ichi%end, &
                                ipsi%beg:ipsi%end))

        allocate (alf_sum%sf( &
                  ix%beg:ix%end, &
                  iy%beg:iy%end, &
                  iz%beg:iz%end))

        do k = ipsi%beg, ipsi%end
            do j = ichi%beg, ichi%end
                do i = 1, num_dims

                    allocate (qL_cons_ndqp(i, j, k)%vf(1:sys_size))
                    allocate (qR_cons_ndqp(i, j, k)%vf(1:sys_size))
                    allocate (qL_prim_ndqp(i, j, k)%vf(1:sys_size))
                    allocate (qR_prim_ndqp(i, j, k)%vf(1:sys_size))

                    allocate (myflux_vf(i, j, k)%vf(1:sys_size))
                    allocate (myflux_src_vf(i, j, k)%vf(1:sys_size))

                    if (abs(j) >= abs(k)) then
                        do l = 1, sys_size
                            allocate (myflux_vf(i, j, k)%vf(l)%sf( &
                                      ix%beg:ix%end, &
                                      iy%beg:iy%end, &
                                      iz%beg:iz%end))
                            allocate (myflux_src_vf(i, j, k)%vf(l)%sf( &
                                      ix%beg:ix%end, &
                                      iy%beg:iy%end, &
                                      iz%beg:iz%end))
                        end do

                        if (i == 1) then

                            do l = 1, cont_idx%end
                                allocate (qL_cons_ndqp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (qR_cons_ndqp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                            end do

                            if (weno_vars == 1) then
                                do l = mom_idx%beg, E_idx
                                    allocate (qL_cons_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (qR_cons_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                            end if

                            do l = mom_idx%beg, E_idx
                                allocate (qL_prim_ndqp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (qR_prim_ndqp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                            end do
                            if (model_eqns == 3) then
                                do l = internalEnergies_idx%beg, internalEnergies_idx%end
                                    allocate (qL_prim_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (qR_prim_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                            end if

                            do l = adv_idx%beg, sys_size
                                allocate (qL_cons_ndqp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (qR_cons_ndqp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                            end do

                            if (bubbles) then
                                do l = bub_idx%beg, bub_idx%end
                                    allocate (qL_prim_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (qR_prim_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                            end if

                            if (hypoelasticity) then
                                do l = stress_idx%beg, stress_idx%end
                                    allocate (qL_prim_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (qR_prim_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                                do l = mom_idx%beg, mom_idx%end
                                    allocate (qL_prim_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (qR_prim_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                            end if
                        else

                            do l = 1, sys_size
                                qL_cons_ndqp(i, j, k)%vf(l)%sf => &
                                    qL_cons_ndqp(1, j, k)%vf(l)%sf
                                qR_cons_ndqp(i, j, k)%vf(l)%sf => &
                                    qR_cons_ndqp(1, j, k)%vf(l)%sf
                                qL_prim_ndqp(i, j, k)%vf(l)%sf => &
                                    qL_prim_ndqp(1, j, k)%vf(l)%sf
                                qR_prim_ndqp(i, j, k)%vf(l)%sf => &
                                    qR_prim_ndqp(1, j, k)%vf(l)%sf
                            end do

                            if ((char_decomp .neqv. .true.) &
                                .and. &
                                any(Re_size > 0)) then
                                if (weno_vars == 1) then
                                    do l = 1, mom_idx%end
                                        allocate (qL_cons_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                        allocate (qR_cons_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                    end do
                                else
                                    do l = mom_idx%beg, mom_idx%end
                                        allocate (qL_prim_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                        allocate (qR_prim_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                    end do
                                    if (model_eqns == 3) then
                                        do l = internalEnergies_idx%beg, internalEnergies_idx%end
                                            allocate (qL_prim_ndqp(i, j, k)%vf(l)%sf( &
                                                      ix%beg:ix%end, &
                                                      iy%beg:iy%end, &
                                                      iz%beg:iz%end))
                                            allocate (qR_prim_ndqp(i, j, k)%vf(l)%sf( &
                                                      ix%beg:ix%end, &
                                                      iy%beg:iy%end, &
                                                      iz%beg:iz%end))
                                        end do
                                    end if
                                end if
                            end if

                            do l = 1, crv_size
                                allocate (qL_cons_ndqp(i, j, k)%vf( &
                                          E_idx + crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (qR_cons_ndqp(i, j, k)%vf( &
                                          E_idx + crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                            end do

                        end if

                        do l = 1, cont_idx%end
                            qL_prim_ndqp(i, j, k)%vf(l)%sf => &
                                qL_cons_ndqp(i, j, k)%vf(l)%sf
                            qR_prim_ndqp(i, j, k)%vf(l)%sf => &
                                qR_cons_ndqp(i, j, k)%vf(l)%sf
                        end do

                        !DO l = adv_idx%beg, adv_idx%end
                        !    qL_prim_ndqp(i,j,k)%vf(l)%sf => &
                        !    qL_cons_ndqp(i,j,k)%vf(l)%sf
                        !    qR_prim_ndqp(i,j,k)%vf(l)%sf => &
                        !    qR_cons_ndqp(i,j,k)%vf(l)%sf
                        !END DO

                        if (adv_alphan) then
                            do l = adv_idx%beg, adv_idx%end
                                qL_prim_ndqp(i, j, k)%vf(l)%sf => &
                                    qL_cons_ndqp(i, j, k)%vf(l)%sf
                                qR_prim_ndqp(i, j, k)%vf(l)%sf => &
                                    qR_cons_ndqp(i, j, k)%vf(l)%sf
                            end do
                        else
                            do l = adv_idx%beg, adv_idx%end + 1
                                qL_prim_ndqp(i, j, k)%vf(l)%sf => &
                                    qL_cons_ndqp(i, j, k)%vf(l)%sf
                                qR_prim_ndqp(i, j, k)%vf(l)%sf => &
                                    qR_cons_ndqp(i, j, k)%vf(l)%sf
                            end do
                        end if

                    end if

                end do
            end do
        end do
        ! END: Allocation/Association of qK_cons_ndqp and qK_prim_ndqp =====

        ! Allocation of dq_prim_ds_qp ======================================

        if (any(Re_size > 0) .or. We_size > 0 .or. tvd_rhs_flux .or. hypoelasticity) then

            allocate (dq_prim_dx_qp(ieta%beg:ieta%end, &
                                    iksi%beg:iksi%end, &
                                    itau%beg:itau%end))
            allocate (dq_prim_dy_qp(ieta%beg:ieta%end, &
                                    iksi%beg:iksi%end, &
                                    itau%beg:itau%end))
            allocate (dq_prim_dz_qp(ieta%beg:ieta%end, &
                                    iksi%beg:iksi%end, &
                                    itau%beg:itau%end))
            allocate (gm_vel_qp(ieta%beg:ieta%end, &
                                iksi%beg:iksi%end, &
                                itau%beg:itau%end))

            do k = itau%beg, itau%end
                do j = iksi%beg, iksi%end
                    do i = ieta%beg, ieta%end

                        allocate (dq_prim_dx_qp(i, j, k)%vf(1:sys_size))
                        allocate (dq_prim_dy_qp(i, j, k)%vf(1:sys_size))
                        allocate (dq_prim_dz_qp(i, j, k)%vf(1:sys_size))
                        allocate (gm_vel_qp(i, j, k)%vf(1:sys_size))

                        if (hypoelasticity .or. all((/i, j, k/) == 0) .and. any(Re_size > 0)) then

                            do l = mom_idx%beg, mom_idx%end
                                allocate (dq_prim_dx_qp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (gm_vel_qp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                            end do

                            if (n > 0) then

                                do l = mom_idx%beg, mom_idx%end
                                    allocate (dq_prim_dy_qp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do

                                if (p > 0) then
                                    do l = mom_idx%beg, mom_idx%end
                                        allocate (dq_prim_dz_qp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                    end do
                                end if

                            end if

                        end if

                        if ((i /= 0 .and. abs(j) >= abs(k)) &
                            .or. &
                            all((/i, j, k/) == 0)) then

                            do l = 1, crv_size
                                allocate (dq_prim_dx_qp(i, j, k)%vf( &
                                          E_idx + crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                if (n > 0) then
                                    allocate (dq_prim_dy_qp(i, j, k)%vf( &
                                              E_idx + crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end if
                            end do

                            if (p > 0) then
                                do l = 1, crv_size
                                    allocate (dq_prim_dz_qp(i, j, k)%vf( &
                                              E_idx + crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                            end if

                        end if

                    end do
                end do
            end do

        end if
        ! END: Allocation of dq_prim_ds_qp =================================

        ! Allocation/Association of dqK_prim_ds_ndqp =======================
        allocate (dqL_prim_dx_ndqp(1:num_dims, &
                                   ichi%beg:ichi%end, &
                                   ipsi%beg:ipsi%end))
        allocate (dqL_prim_dy_ndqp(1:num_dims, &
                                   ichi%beg:ichi%end, &
                                   ipsi%beg:ipsi%end))
        allocate (dqL_prim_dz_ndqp(1:num_dims, &
                                   ichi%beg:ichi%end, &
                                   ipsi%beg:ipsi%end))
        allocate (dqR_prim_dx_ndqp(1:num_dims, &
                                   ichi%beg:ichi%end, &
                                   ipsi%beg:ipsi%end))
        allocate (dqR_prim_dy_ndqp(1:num_dims, &
                                   ichi%beg:ichi%end, &
                                   ipsi%beg:ipsi%end))
        allocate (dqR_prim_dz_ndqp(1:num_dims, &
                                   ichi%beg:ichi%end, &
                                   ipsi%beg:ipsi%end))

        if (any(Re_size > 0) .or. We_size > 0 .or. hypoelasticity) then
            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end
                    do i = 1, num_dims

                        allocate (dqL_prim_dx_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (dqL_prim_dy_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (dqL_prim_dz_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (dqR_prim_dx_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (dqR_prim_dy_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (dqR_prim_dz_ndqp(i, j, k)%vf(1:sys_size))

                        if (abs(j) >= abs(k)) then

                            if (any(Re_size > 0) .or. hypoelasticity) then

                                do l = mom_idx%beg, mom_idx%end
                                    allocate (dqL_prim_dx_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (dqR_prim_dx_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do

                                if (n > 0) then
                                    do l = mom_idx%beg, mom_idx%end
                                        allocate (dqL_prim_dy_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                        allocate (dqR_prim_dy_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                    end do
                                end if

                                if (p > 0) then
                                    do l = mom_idx%beg, mom_idx%end
                                        allocate (dqL_prim_dz_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                        allocate (dqR_prim_dz_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                    end do
                                end if

                            end if

                            do l = 1, crv_size
                                allocate (dqL_prim_dx_ndqp(i, j, k)%vf( &
                                          E_idx + crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (dqR_prim_dx_ndqp(i, j, k)%vf( &
                                          E_idx + crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                if (n > 0) then
                                    allocate (dqL_prim_dy_ndqp(i, j, k)%vf( &
                                              E_idx + crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (dqR_prim_dy_ndqp(i, j, k)%vf( &
                                              E_idx + crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    if (p > 0) then
                                        allocate (dqL_prim_dz_ndqp(i, j, k)%vf( &
                                                  E_idx + crv_idx(l))%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                        allocate (dqR_prim_dz_ndqp(i, j, k)%vf( &
                                                  E_idx + crv_idx(l))%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                    end if
                                end if
                            end do

                        end if

                    end do
                end do
            end do
        end if
        ! END: Allocation/Association of dqK_prim_ds_ndqp ==================

        ! Allocation of gm_alpha_qp ========================================
        allocate (gm_alpha_qp(ieta%beg:ieta%end, &
                              iksi%beg:iksi%end, &
                              itau%beg:itau%end))

        if (We_size > 0) then
            do k = itau%beg, itau%end
                do j = iksi%beg, iksi%end
                    do i = ieta%beg, ieta%end

                        allocate (gm_alpha_qp(i, j, k)%vf(1:num_fluids))

                        if (any(sum(abs((/i, j, k/))) == (/0, num_dims/))) then
                            do l = 1, crv_size
                                allocate (gm_alpha_qp(i, j, k)%vf(crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                            end do
                        end if

                    end do
                end do
            end do
        end if
        ! ==================================================================

        ! Allocation of gm_alphaK_ndqp =====================================
        allocate (gm_alphaL_ndqp(1:num_dims, &
                                 ichi%beg:ichi%end, &
                                 ipsi%beg:ipsi%end))
        allocate (gm_alphaR_ndqp(1:num_dims, &
                                 ichi%beg:ichi%end, &
                                 ipsi%beg:ipsi%end))

        if (We_size > 0) then
            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end
                    do i = 1, num_dims

                        allocate (gm_alphaL_ndqp(i, j, k)%vf(1:num_fluids))
                        allocate (gm_alphaR_ndqp(i, j, k)%vf(1:num_fluids))

                        if (abs(j) >= abs(k)) then
                            do l = 1, crv_size
                                allocate (gm_alphaL_ndqp(i, j, k)%vf( &
                                          crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (gm_alphaR_ndqp(i, j, k)%vf( &
                                          crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                            end do
                        end if

                    end do
                end do
            end do
        end if
        ! ==================================================================

        ! Allocation of un_alphaK_s_ndqp ===================================
        if (We_size > 0) then

            allocate (un_alphaL_x_ndqp(1:num_dims, &
                                       ichi%beg:ichi%end, &
                                       ipsi%beg:ipsi%end))
            allocate (un_alphaL_y_ndqp(1:num_dims, &
                                       ichi%beg:ichi%end, &
                                       ipsi%beg:ipsi%end))
            allocate (un_alphaL_z_ndqp(1:num_dims, &
                                       ichi%beg:ichi%end, &
                                       ipsi%beg:ipsi%end))
            allocate (un_alphaR_x_ndqp(1:num_dims, &
                                       ichi%beg:ichi%end, &
                                       ipsi%beg:ipsi%end))
            allocate (un_alphaR_y_ndqp(1:num_dims, &
                                       ichi%beg:ichi%end, &
                                       ipsi%beg:ipsi%end))
            allocate (un_alphaR_z_ndqp(1:num_dims, &
                                       ichi%beg:ichi%end, &
                                       ipsi%beg:ipsi%end))

            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end

                    do i = 1, num_dims
                        allocate (un_alphaL_x_ndqp(i, j, k)%vf(1:num_fluids))
                        allocate (un_alphaL_y_ndqp(i, j, k)%vf(1:num_fluids))
                        allocate (un_alphaL_z_ndqp(i, j, k)%vf(1:num_fluids))
                        allocate (un_alphaR_x_ndqp(i, j, k)%vf(1:num_fluids))
                        allocate (un_alphaR_y_ndqp(i, j, k)%vf(1:num_fluids))
                        allocate (un_alphaR_z_ndqp(i, j, k)%vf(1:num_fluids))
                    end do

                    if (abs(j) >= abs(k)) then

                        do l = 1, crv_size
                            allocate (un_alphaL_x_ndqp(1, j, k)%vf(crv_idx(l))%sf( &
                                      ix%beg:ix%end, &
                                      iy%beg:iy%end, &
                                      iz%beg:iz%end))
                            allocate (un_alphaR_x_ndqp(1, j, k)%vf(crv_idx(l))%sf( &
                                      ix%beg:ix%end, &
                                      iy%beg:iy%end, &
                                      iz%beg:iz%end))
                            if (n > 0) then
                                allocate (un_alphaL_y_ndqp(2, j, k)%vf(crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (un_alphaR_y_ndqp(2, j, k)%vf(crv_idx(l))%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                if (p > 0) then
                                    allocate (un_alphaL_z_ndqp(3, j, k)%vf(crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (un_alphaR_z_ndqp(3, j, k)%vf(crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end if
                            end if
                        end do
                    end if

                end do
            end do

        end if
        ! END: Allocation of un_alphaK_s_ndqp ==============================

        ! Allocation of kappa_vf ===========================================
        if (We_size > 0) then

            allocate (kappa_vf(1:num_fluids))

            do i = 1, crv_size
                allocate (kappa_vf(crv_idx(i))%sf(ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
            end do

        end if
        ! ==================================================================

        ! Allocation/Association of kappaK_ndqp ============================
        allocate (kappaL_ndqp(1:num_dims, &
                              ichi%beg:ichi%end, &
                              ipsi%beg:ipsi%end))
        allocate (kappaR_ndqp(1:num_dims, &
                              ichi%beg:ichi%end, &
                              ipsi%beg:ipsi%end))

        if (We_size > 0) then
            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end
                    do i = 1, num_dims

                        allocate (kappaL_ndqp(i, j, k)%vf(1:num_fluids))
                        allocate (kappaR_ndqp(i, j, k)%vf(1:num_fluids))

                        if (abs(j) >= abs(k)) then
                            if (i == 1) then
                                do l = 1, crv_size
                                    allocate (kappaL_ndqp(i, j, k)%vf( &
                                              crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (kappaR_ndqp(i, j, k)%vf( &
                                              crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                            else
                                do l = 1, crv_size
                                    kappaL_ndqp(i, j, k)%vf(crv_idx(l))%sf => &
                                        kappaL_ndqp(1, j, k)%vf(crv_idx(l))%sf
                                    kappaR_ndqp(i, j, k)%vf(crv_idx(l))%sf => &
                                        kappaR_ndqp(1, j, k)%vf(crv_idx(l))%sf
                                end do
                            end if
                        end if

                    end do
                end do
            end do
        end if
        ! ==================================================================

        if (alt_crv) then
            ! Allocation of dgm_alpha_ds_qp & un_alpha_s_qp ================
            allocate (dgm_alpha_dx_qp(ieta%beg:ieta%end, &
                                      iksi%beg:iksi%end, &
                                      itau%beg:itau%end))
            allocate (dgm_alpha_dy_qp(ieta%beg:ieta%end, &
                                      iksi%beg:iksi%end, &
                                      itau%beg:itau%end))
            allocate (dgm_alpha_dz_qp(ieta%beg:ieta%end, &
                                      iksi%beg:iksi%end, &
                                      itau%beg:itau%end))
            allocate (un_alpha_x_qp(ieta%beg:ieta%end, &
                                    iksi%beg:iksi%end, &
                                    itau%beg:itau%end))
            allocate (un_alpha_y_qp(ieta%beg:ieta%end, &
                                    iksi%beg:iksi%end, &
                                    itau%beg:itau%end))
            allocate (un_alpha_z_qp(ieta%beg:ieta%end, &
                                    iksi%beg:iksi%end, &
                                    itau%beg:itau%end))

            if (We_size > 0) then
                do k = itau%beg, itau%end
                    do j = iksi%beg, iksi%end
                        do i = ieta%beg, ieta%end
                            allocate (dgm_alpha_dx_qp(i, j, k)%vf(1:num_fluids))
                            allocate (dgm_alpha_dy_qp(i, j, k)%vf(1:num_fluids))
                            allocate (dgm_alpha_dz_qp(i, j, k)%vf(1:num_fluids))
                            allocate (un_alpha_x_qp(i, j, k)%vf(1:num_fluids))
                            allocate (un_alpha_y_qp(i, j, k)%vf(1:num_fluids))
                            allocate (un_alpha_z_qp(i, j, k)%vf(1:num_fluids))

                            if (any(sum(abs((/i, j, k/))) == (/0, num_dims/))) then
                                do l = 1, crv_size
                                    allocate (dgm_alpha_dx_qp(i, j, k)%vf(crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (dgm_alpha_dy_qp(i, j, k)%vf(crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (dgm_alpha_dz_qp(i, j, k)%vf(crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (un_alpha_x_qp(i, j, k)%vf(crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (un_alpha_y_qp(i, j, k)%vf(crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (un_alpha_z_qp(i, j, k)%vf(crv_idx(l))%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                            end if

                        end do
                    end do
                end do
            end if
        end if
        ! ==============================================================

        ! Allocation of laplacian_vf and alt_kappa_vf
        if (We_size > 0) then
            if (alt_crv) then
                allocate (laplacian_vf(1:num_fluids))
                allocate (alt_kappa_vf(1:num_fluids))

                do i = 1, crv_size
                    allocate (laplacian_vf(crv_idx(i))%sf(ix%beg:ix%end, &
                                                          iy%beg:iy%end, &
                                                          iz%beg:iz%end))
                    allocate (alt_kappa_vf(crv_idx(i))%sf(ix%beg:ix%end, &
                                                          iy%beg:iy%end, &
                                                          iz%beg:iz%end))
                end do
            end if

            if (We_src) then
                allocate (We_mtm_src(0:m, 0:n, 0:p))
                allocate (We_nrg_src(0:m, 0:n, 0:p))
            end if

        end if

        ! Allocation of regularization terms
        if (regularization) then
            allocate (reg_src_vf(1:sys_size))
            do i = 1, sys_size
                allocate (reg_src_vf(i)%sf(0:m, 0:n, 0:p))
            end do
        end if

        if (bubbles) then
            allocate (bub_adv_src(0:m, 0:n, 0:p))
            if (qbmm) then
                allocate (bub_mom_src(1:nb, 1:nmom, 0:m, 0:n, 0:p))
            else
                allocate (bub_r_src(1:nb, 0:m, 0:n, 0:p))
                allocate (bub_v_src(1:nb, 0:m, 0:n, 0:p))
                allocate (bub_p_src(1:nb, 0:m, 0:n, 0:p))
                allocate (bub_m_src(1:nb, 0:m, 0:n, 0:p))
            end if
        end if

        if (monopole) then
            allocate (mono_mass_src(0:m, 0:n, 0:p))
            allocate (mono_mom_src(1:num_dims, 0:m, 0:n, 0:p))
            allocate (mono_E_src(0:m, 0:n, 0:p))
        end if

        allocate (divu%sf( &
                  ix%beg:ix%end, &
                  iy%beg:iy%end, &
                  iz%beg:iz%end))

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -1; if (n > 0) iy%beg = -1; if (p > 0) iz%beg = -1

        ix%end = m; iy%end = n; iz%end = p
        ! ==================================================================

        ! Allocation/Association of flux_ndqp, flux_src_ndqp, and flux_gsrc_ndqp ===
        allocate (flux_ndqp(1:num_dims, &
                            ichi%beg:ichi%end, &
                            ipsi%beg:ipsi%end))
        allocate (flux_src_ndqp(1:num_dims, &
                                ichi%beg:ichi%end, &
                                ipsi%beg:ipsi%end))
        allocate (flux_gsrc_ndqp(1:num_dims, &
                                 ichi%beg:ichi%end, &
                                 ipsi%beg:ipsi%end))

        ! Allocation of variables for flux limiter
        if (tvd_rhs_flux) then
            allocate (flux_lim_func(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))

            allocate (lo_flux_ndqp(1:num_dims, &
                                   ichi%beg:ichi%end, &
                                   ipsi%beg:ipsi%end))
            allocate (lo_flux_src_ndqp(1:num_dims, &
                                       ichi%beg:ichi%end, &
                                       ipsi%beg:ipsi%end))
            allocate (lo_flux_gsrc_ndqp(1:num_dims, &
                                        ichi%beg:ichi%end, &
                                        ipsi%beg:ipsi%end))
            allocate (hi_flux_ndqp(1:num_dims, &
                                   ichi%beg:ichi%end, &
                                   ipsi%beg:ipsi%end))
            allocate (hi_flux_src_ndqp(1:num_dims, &
                                       ichi%beg:ichi%end, &
                                       ipsi%beg:ipsi%end))
            allocate (hi_flux_gsrc_ndqp(1:num_dims, &
                                        ichi%beg:ichi%end, &
                                        ipsi%beg:ipsi%end))
        end if

        do k = ipsi%beg, ipsi%end
            do j = ichi%beg, ichi%end
                do i = 1, num_dims

                    allocate (flux_ndqp(i, j, k)%vf(1:sys_size))
                    allocate (flux_src_ndqp(i, j, k)%vf(1:sys_size))
                    allocate (flux_gsrc_ndqp(i, j, k)%vf(1:sys_size))

                    if (abs(j) >= abs(k)) then
                        if (i == 1) then

                            do l = 1, sys_size
                                allocate (flux_ndqp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (flux_gsrc_ndqp(i, j, k)%vf(l)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                            end do

                            if (any(Re_size > 0) .or. We_size > 0 .or. hypoelasticity) then
                                do l = mom_idx%beg, E_idx
                                    allocate (flux_src_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                            end if

                            allocate (flux_src_ndqp(i, j, k)%vf(adv_idx%beg)%sf( &
                                      ix%beg:ix%end, &
                                      iy%beg:iy%end, &
                                      iz%beg:iz%end))
                            if (riemann_solver == 1) then
                                do l = adv_idx%beg + 1, adv_idx%end
                                    allocate (flux_src_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do
                            else
                                !IF ( (num_fluids > 1) .AND. (bubbles .NEQV. .TRUE.)) THEN
                                do l = adv_idx%beg + 1, adv_idx%end
                                    flux_src_ndqp(i, j, k)%vf(l)%sf => &
                                        flux_src_ndqp(i, j, k)%vf(adv_idx%beg)%sf
                                end do
                                !END IF
                            end if

                        else

                            do l = 1, sys_size
                                flux_ndqp(i, j, k)%vf(l)%sf => &
                                    flux_ndqp(1, j, k)%vf(l)%sf
                                flux_src_ndqp(i, j, k)%vf(l)%sf => &
                                    flux_src_ndqp(1, j, k)%vf(l)%sf
                                flux_gsrc_ndqp(i, j, k)%vf(l)%sf => &
                                    flux_gsrc_ndqp(1, j, k)%vf(l)%sf
                            end do

                        end if
                    end if

                end do
            end do
        end do

        if (tvd_rhs_flux) then
            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end
                    do i = 1, num_dims

                        allocate (lo_flux_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (lo_flux_src_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (lo_flux_gsrc_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (hi_flux_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (hi_flux_src_ndqp(i, j, k)%vf(1:sys_size))
                        allocate (hi_flux_gsrc_ndqp(i, j, k)%vf(1:sys_size))

                        if (abs(j) >= abs(k)) then
                            if (i == 1) then
                                do l = 1, sys_size
                                    allocate (lo_flux_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (hi_flux_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (lo_flux_gsrc_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                    allocate (hi_flux_gsrc_ndqp(i, j, k)%vf(l)%sf( &
                                              ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end))
                                end do

                                if (any(Re_size > 0) .or. We_size > 0 .or. hypoelasticity) then
                                    do l = mom_idx%beg, E_idx
                                        allocate (lo_flux_src_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                        allocate (hi_flux_src_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                    end do
                                end if

                                allocate (lo_flux_src_ndqp(i, j, k)%vf(adv_idx%beg)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                allocate (hi_flux_src_ndqp(i, j, k)%vf(adv_idx%beg)%sf( &
                                          ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
                                if (riemann_solver == 1) then
                                    do l = adv_idx%beg + 1, adv_idx%end
                                        allocate (lo_flux_src_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                        allocate (hi_flux_src_ndqp(i, j, k)%vf(l)%sf( &
                                                  ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end))
                                    end do
                                else
                                    do l = adv_idx%beg + 1, adv_idx%end
                                        lo_flux_src_ndqp(i, j, k)%vf(l)%sf => &
                                            lo_flux_src_ndqp(i, j, k)%vf(adv_idx%beg)%sf
                                        hi_flux_src_ndqp(i, j, k)%vf(l)%sf => &
                                            hi_flux_src_ndqp(i, j, k)%vf(adv_idx%beg)%sf
                                    end do
                                end if
                            else
                                do l = 1, sys_size
                                    lo_flux_ndqp(i, j, k)%vf(l)%sf => &
                                        lo_flux_ndqp(1, j, k)%vf(l)%sf
                                    hi_flux_ndqp(i, j, k)%vf(l)%sf => &
                                        hi_flux_ndqp(1, j, k)%vf(l)%sf
                                    lo_flux_src_ndqp(i, j, k)%vf(l)%sf => &
                                        lo_flux_src_ndqp(1, j, k)%vf(l)%sf
                                    hi_flux_src_ndqp(i, j, k)%vf(l)%sf => &
                                        hi_flux_src_ndqp(1, j, k)%vf(l)%sf
                                    lo_flux_gsrc_ndqp(i, j, k)%vf(l)%sf => &
                                        lo_flux_gsrc_ndqp(1, j, k)%vf(l)%sf
                                    hi_flux_gsrc_ndqp(i, j, k)%vf(l)%sf => &
                                        hi_flux_gsrc_ndqp(1, j, k)%vf(l)%sf
                                end do
                            end if
                        end if
                    end do
                end do
            end do
        end if
        ! END: Allocation/Association of flux_ndqp, flux_src_ndqp, and flux_gsrc_ndqp ===

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem
        if (riemann_solver == 1) then
            s_riemann_solver => s_hll_riemann_solver
        elseif (riemann_solver == 2) then
            s_riemann_solver => s_hllc_riemann_solver
        else
            s_riemann_solver => s_exact_riemann_solver
        end if

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables
        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables
        else if (bubbles) then          ! Volume fraction for bubbles
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables_bubbles
        else                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if

    end subroutine s_initialize_rhs_module ! -------------------------------

    !> The purpose of this procedure is to employ the inputted
        !!      cell-average conservative variables in order to compute
        !!      the cell-average RHS variables of the semidiscrete form
        !!      of the governing equations by utilizing the appropriate
        !!      Riemann solver.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param rhs_vf Cell-average RHS variables
        !!  @param t_step Current time-step
    subroutine s_compute_rhs(q_cons_vf, q_prim_vf, rhs_vf, t_step) ! -------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf
        integer, intent(IN) :: t_step

        real(kind(0d0)) :: top, bottom  !< Numerator and denominator when evaluating flux limiter function

        real(kind(0d0)), dimension(0:m, 0:n, 0:p) :: blkmod1, blkmod2, alpha1, alpha2, Kterm !<
            !! terms  for K div(u)

            !! For shear modulus term calculation in stress eq when hypoelastic = true :
        real(kind(0d0)), dimension(num_fluids) :: alpha_K, alpha_rho_K
        real(kind(0d0)) :: rho_K, G_K
        real(kind(0d0)), dimension(0:m, 0:n, 0:p) :: rho_K_field, G_K_field
        real(kind(0d0)) :: gamma_K, pi_inf_K
        real(kind(0d0)), dimension(2) :: Re_K
        real(kind(0d0)), dimension(1:num_fluids, 1:num_fluids) :: We_K

        integer :: i, j, k, l, r, ii !< Generic loop iterators

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================

        if (DEBUG) print *, 'Start rhs'

        ! Association/Population of Working Variables ======================
        do i = 1, sys_size
            q_cons_qp(0, 0, 0)%vf(i)%sf => q_cons_vf(i)%sf
            q_prim_qp(0, 0, 0)%vf(i)%sf => q_prim_vf(i)%sf
        end do

        ! print*, 'cons vars: '
        ! do i = 1,sys_size
        !     print*, 'cons var: ', i, q_cons_qp(0,0,0)%vf(i)%sf(1,0,0)
        ! end do

        call s_populate_conservative_variables_buffers()

        if (DEBUG) print *, 'pop cons vars'
        if ((model_eqns == 2 .or. model_eqns == 3) .and. (adv_alphan .neqv. .true.)) then
            q_cons_qp(0, 0, 0)%vf(sys_size)%sf = 1d0

            do i = adv_idx%beg, adv_idx%end
                q_cons_qp(0, 0, 0)%vf(sys_size)%sf = &
                    q_cons_qp(0, 0, 0)%vf(sys_size)%sf - &
                    q_cons_qp(0, 0, 0)%vf(i)%sf
            end do
        end if

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_cons_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_cons_vf(i)%sf = q_cons_vf(i)%sf*(1.d0 - q_cons_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! ==================================================================

        ! Computing Volume Fraction Gradients and Curvatures ===============
        if (We_size > 0 .and. lsq_deriv) then
            if (p > 0) then
                call s_compute_lsq_gradient_curvature(q_prim_qp(0, 0, 0)%vf, &
                                                      dq_prim_dx_qp(0, 0, 0)%vf, &
                                                      dq_prim_dy_qp(0, 0, 0)%vf, &
                                                      dq_prim_dz_qp(0, 0, 0)%vf, &
                                                      gm_alpha_qp(0, 0, 0)%vf, &
                                                      kappa_vf)
            else
                call s_compute_lsq_gradient_curvature(q_prim_qp(0, 0, 0)%vf, &
                                                      dq_prim_dx_qp(0, 0, 0)%vf, &
                                                      dq_prim_dy_qp(0, 0, 0)%vf, &
                                                      dq_prim_dy_qp(0, 0, 0)%vf, &
                                                      gm_alpha_qp(0, 0, 0)%vf, &
                                                      kappa_vf)
            end if
        end if

        if (crv_size > 0) call s_get_crv(q_cons_vf, q_prim_vf, rhs_vf)

        ! Converting Conservative to Primitive Variables ===================
        iv%beg = 1; iv%end = adv_idx%end

        call s_reconstruct_cell_interior_values(q_cons_qp)

        if ((model_eqns == 2 .or. model_eqns == 3) &
            .and. &
            (adv_alphan .neqv. .true.) &
            .and. &
            commute_err) then
            do k = itau%beg, itau%end, 2
                do j = iksi%beg, iksi%end, 2
                    do i = ieta%beg, ieta%end, 2

                        q_cons_qp(i, j, k)%vf(adv_idx%end)%sf = 1d0

                        do l = adv_idx%beg, adv_idx%end
                            q_cons_qp(i, j, k)%vf(adv_idx%end)%sf = &
                                q_cons_qp(i, j, k)%vf(adv_idx%end)%sf - &
                                q_cons_qp(i, j, k)%vf(l)%sf
                        end do

                    end do
                end do
            end do
        end if

        !convert conservative variables to primitive
        !   (except first and last, \alpha \rho and \alpha)
        !we do this, though there is no actual loop, just (0,0,0)
        do k = itau%beg, itau%end, 2
            do j = iksi%beg, iksi%end, 2
                do i = ieta%beg, ieta%end, 2
                    call s_convert_conservative_to_primitive_variables( &
                        q_cons_qp(i, j, k)%vf, &
                        q_prim_qp(i, j, k)%vf, &
                        gm_alpha_qp(i, j, k)%vf, &
                        ix, iy, iz)
                end do
            end do
        end do

        if (DEBUG) print *, 'conv to prim vars'

        iv%beg = mom_idx%beg; iv%end = E_idx

        call s_average_cell_interior_values(q_prim_qp)

        if (DEBUG) print *, 'got cell interior values'
        if (t_step == t_step_stop) return
        ! ==================================================================

        ! Computing Velocity Gradients =

        if (any(Re_size > 0) .or. hypoelasticity) call s_get_viscous(q_cons_vf, q_prim_vf, rhs_vf)

        if (DEBUG) print *, 'Before qbmm'
        ! compute required moments
        if (qbmm) call s_mom_inv(q_prim_vf, mom_sp, mom_3d, ix, iy, iz)

        ! Dimensional Splitting Loop =======================================
        do i = 1, num_dims

            ! Configuring Coordinate Direction Indexes ======================
            ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

            if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

            ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
            ! ===============================================================

            ! Reconstructing Primitive/Conservative Variables ===============
            if (char_decomp .or. (all(Re_size == 0) .and. We_size == 0)) then

                iv%beg = 1; 
                if (adv_alphan) then
                    iv%end = adv_idx%end
                    if (bubbles) iv%end = sys_size
                    if (hypoelasticity) iv%end = sys_size
                else
                    iv%end = adv_idx%end + 1
                end if

                !reconstruct either primitive or conservative vars
                if (weno_vars == 1) then
                    call s_reconstruct_cell_boundary_values( &
                        q_cons_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                        qL_cons_ndqp(i, :, :), &
                        qR_cons_ndqp(i, :, :), &
                        weno_vars, i)
                else
                    call s_reconstruct_cell_boundary_values( &
                        q_prim_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                        qL_prim_ndqp(i, :, :), &
                        qR_prim_ndqp(i, :, :), &
                        weno_vars, i)
                end if

            else
                ! ===============================================================

                ! Reconstructing Continuity Variables ===========================
                if (weno_vars == 2 .or. all(Re_size == 0)) then

                    iv%beg = cont_idx%beg; iv%end = cont_idx%end

                    call s_reconstruct_cell_boundary_values( &
                        q_cons_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                        qL_cons_ndqp(i, :, :), &
                        qR_cons_ndqp(i, :, :), &
                        dflt_int, i)

                end if
                ! ===============================================================

                ! Reconstructing Momentum/Velocity Variables ====================
                if (all(Re_size == 0)) then

                    iv%beg = mom_idx%beg; iv%end = mom_idx%end

                    if (weno_vars == 1) then
                        call s_reconstruct_cell_boundary_values( &
                            q_cons_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                            qL_cons_ndqp(i, :, :), &
                            qR_cons_ndqp(i, :, :), &
                            dflt_int, i)
                    else
                        call s_reconstruct_cell_boundary_values( &
                            q_prim_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                            qL_prim_ndqp(i, :, :), &
                            qR_prim_ndqp(i, :, :), &
                            dflt_int, i)
                    end if

                end if
                ! ===============================================================

                ! Reconstructing Partial or Mixture Energy/Pressure Variables ===
                iv%beg = E_idx; iv%end = iv%beg

                if (weno_vars == 1) then
                    call s_reconstruct_cell_boundary_values( &
                        q_cons_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                        qL_cons_ndqp(i, :, :), &
                        qR_cons_ndqp(i, :, :), &
                        dflt_int, i)
                else
                    call s_reconstruct_cell_boundary_values( &
                        q_prim_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                        qL_prim_ndqp(i, :, :), &
                        qR_prim_ndqp(i, :, :), &
                        dflt_int, i)
                end if
                ! ===============================================================

                ! Reconstructing Volume Fraction Variables ======================
                if (We_size == 0) then

                    iv%beg = adv_idx%beg; iv%end = adv_idx%end

                    call s_reconstruct_cell_boundary_values( &
                        q_cons_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                        qL_cons_ndqp(i, :, :), &
                        qR_cons_ndqp(i, :, :), &
                        dflt_int, i)

                else

                    do l = adv_idx%beg, adv_idx%end
                        if (all(crv_idx /= l - E_idx)) then

                            iv%beg = l; iv%end = iv%beg

                            call s_reconstruct_cell_boundary_values( &
                                q_cons_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                                qL_cons_ndqp(i, :, :), &
                                qR_cons_ndqp(i, :, :), &
                                dflt_int, i)

                        end if
                    end do

                end if
                ! ===============================================================

            end if

            if ((model_eqns == 2 .or. model_eqns == 3) .and. (adv_alphan .neqv. .true.)) then
                do k = ipsi%beg, ipsi%end, 2
                    do j = ichi%beg, ichi%end, 2

                        qL_cons_ndqp(i, j, k)%vf(sys_size)%sf = 1d0
                        qR_cons_ndqp(i, j, k)%vf(sys_size)%sf = 1d0

                        do l = adv_idx%beg, adv_idx%end

                            qL_cons_ndqp(i, j, k)%vf(sys_size)%sf = &
                                qL_cons_ndqp(i, j, k)%vf(sys_size)%sf - &
                                qL_cons_ndqp(i, j, k)%vf(l)%sf

                            qR_cons_ndqp(i, j, k)%vf(sys_size)%sf = &
                                qR_cons_ndqp(i, j, k)%vf(sys_size)%sf - &
                                qR_cons_ndqp(i, j, k)%vf(l)%sf
                        end do

                    end do
                end do
            end if
            ! END: Reconstructing Volume Fraction Variables =================

            ! Converting Conservative to Primitive Variables ================
            if (weno_vars == 1) then
                do k = ipsi%beg, ipsi%end, 2
                    do j = ichi%beg, ichi%end, 2
                        call s_convert_conservative_to_primitive_variables( &
                            qL_cons_ndqp(i, j, k)%vf, &
                            qL_prim_ndqp(i, j, k)%vf, &
                            gm_alphaL_ndqp(i, j, k)%vf, &
                            ix, iy, iz)
                        call s_convert_conservative_to_primitive_variables( &
                            qR_cons_ndqp(i, j, k)%vf, &
                            qR_prim_ndqp(i, j, k)%vf, &
                            gm_alphaR_ndqp(i, j, k)%vf, &
                            ix, iy, iz)
                    end do
                end do
            end if
            ! ===============================================================

            ! Reconstructing First-Order Spatial Derivatives of Velocity ====
            if (any(Re_size > 0) .or. hypoelasticity) then

                iv%beg = mom_idx%beg; iv%end = mom_idx%end

                if (weno_Re_flux .or. hypoelasticity) then

                    call s_reconstruct_cell_boundary_values( &
                        dq_prim_dx_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                        dqL_prim_dx_ndqp(i, :, :), &
                        dqR_prim_dx_ndqp(i, :, :), &
                        dflt_int, i)

                    if (n > 0) then

                        call s_reconstruct_cell_boundary_values( &
                            dq_prim_dy_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                            dqL_prim_dy_ndqp(i, :, :), &
                            dqR_prim_dy_ndqp(i, :, :), &
                            dflt_int, i)
                        if (p > 0) then
                            call s_reconstruct_cell_boundary_values( &
                                dq_prim_dz_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                                dqL_prim_dz_ndqp(i, :, :), &
                                dqR_prim_dz_ndqp(i, :, :), &
                                dflt_int, i)
                        end if

                    end if

                end if

            end if
            ! ===============================================================

            ! Reconstructing Curvatures =====================================
            do l = 1, crv_size

                iv%beg = crv_idx(l); iv%end = iv%beg

                if (alt_crv .neqv. .true.) then
                    call s_reconstruct_cell_boundary_values( &
                        kappa_vf(iv%beg:iv%end), &
                        kappaL_ndqp(i, :, :), &
                        kappaR_ndqp(i, :, :), &
                        dflt_int, i)
                else
                    call s_reconstruct_cell_boundary_values( &
                        alt_kappa_vf(iv%beg:iv%end), &
                        kappaL_ndqp(i, :, :), &
                        kappaR_ndqp(i, :, :), &
                        dflt_int, i)
                end if

            end do
            ! ===============================================================

            ! Configuring Coordinate Direction Indexes ======================
            if (i == 1) then
                ix%beg = -1; iy%beg = 0; iz%beg = 0
            elseif (i == 2) then
                ix%beg = 0; iy%beg = -1; iz%beg = 0
            else
                ix%beg = 0; iy%beg = 0; iz%beg = -1
            end if

            ix%end = m; iy%end = n; iz%end = p
            ! ===============================================================

            if (tvd_rhs_flux) then
                ! Compute Riemann solver flux and source flux
                do k = ipsi%beg, ipsi%end, 2
                    do j = ichi%beg, ichi%end, 2
                        call s_riemann_solver(qR_prim_ndqp(i, j, k)%vf, &
                                              dqR_prim_dx_ndqp(i, j, k)%vf, &
                                              dqR_prim_dy_ndqp(i, j, k)%vf, &
                                              dqR_prim_dz_ndqp(i, j, k)%vf, &
                                              gm_alphaR_ndqp(i, j, k)%vf, &
                                              kappaR_ndqp(i, j, k)%vf, &
                                              qL_prim_ndqp(i, j, k)%vf, &
                                              dqL_prim_dx_ndqp(i, j, k)%vf, &
                                              dqL_prim_dy_ndqp(i, j, k)%vf, &
                                              dqL_prim_dz_ndqp(i, j, k)%vf, &
                                              gm_alphaL_ndqp(i, j, k)%vf, &
                                              kappaL_ndqp(i, j, k)%vf, &
                                              q_prim_qp(0, 0, 0)%vf, &
                                              hi_flux_ndqp(i, j, k)%vf, &
                                              hi_flux_src_ndqp(i, j, k)%vf, &
                                              hi_flux_gsrc_ndqp(i, j, k)%vf, &
                                              i, ix, iy, iz)
                    end do
                end do

                do l = itau%beg, itau%end, 2
                    do k = iksi%beg, iksi%end, 2
                        do j = ieta%beg, ieta%end, 2
                            call s_riemann_solver(q_prim_qp(j, k, l)%vf, &
                                                  dq_prim_dx_qp(j, k, l)%vf, &
                                                  dq_prim_dy_qp(j, k, l)%vf, &
                                                  dq_prim_dz_qp(j, k, l)%vf, &
                                                  gm_alpha_qp(j, k, l)%vf, &
                                                  kappa_vf, &
                                                  q_prim_qp(j, k, l)%vf, &
                                                  dq_prim_dx_qp(j, k, l)%vf, &
                                                  dq_prim_dy_qp(j, k, l)%vf, &
                                                  dq_prim_dz_qp(j, k, l)%vf, &
                                                  gm_alpha_qp(j, k, l)%vf, &
                                                  kappa_vf, &
                                                  q_prim_qp(0, 0, 0)%vf, &
                                                  lo_flux_ndqp(i, k, l)%vf, &
                                                  lo_flux_src_ndqp(i, k, l)%vf, &
                                                  lo_flux_gsrc_ndqp(i, k, l)%vf, &
                                                  i, ix, iy, iz)
                        end do
                    end do
                end do
            else
                ! Computing Riemann Solver Flux and Source Flux =================
                if (DEBUG) print *, 'about to call s_riemann_solver'
                do k = ipsi%beg, ipsi%end, 2
                    do j = ichi%beg, ichi%end, 2
                        call s_riemann_solver(qR_prim_ndqp(i, j, k)%vf, &
                                              dqR_prim_dx_ndqp(i, j, k)%vf, &
                                              dqR_prim_dy_ndqp(i, j, k)%vf, &
                                              dqR_prim_dz_ndqp(i, j, k)%vf, &
                                              gm_alphaR_ndqp(i, j, k)%vf, &
                                              kappaR_ndqp(i, j, k)%vf, &
                                              qL_prim_ndqp(i, j, k)%vf, &
                                              dqL_prim_dx_ndqp(i, j, k)%vf, &
                                              dqL_prim_dy_ndqp(i, j, k)%vf, &
                                              dqL_prim_dz_ndqp(i, j, k)%vf, &
                                              gm_alphaL_ndqp(i, j, k)%vf, &
                                              kappaL_ndqp(i, j, k)%vf, &
                                              q_prim_qp(0, 0, 0)%vf, &
                                              flux_ndqp(i, j, k)%vf, &
                                              flux_src_ndqp(i, j, k)%vf, &
                                              flux_gsrc_ndqp(i, j, k)%vf, &
                                              i, ix, iy, iz)
                    end do
                end do

                ! do j = 1,sys_size
                !     print*, 'fluxes ', flux_ndqp(i,0,0)%vf(j)%sf(:,0,0)
                ! end do
                ! call s_mpi_abort()

                iv%beg = 1; iv%end = adv_idx%end

                call s_average_cell_boundary_values(flux_ndqp(i, :, :))

                if (any(Re_size > 0) .or. We_size > 0 .or. hypoelasticity) then
                    iv%beg = mom_idx%beg
                else
                    iv%beg = adv_idx%beg
                end if

                if (riemann_solver /= 1) iv%end = adv_idx%beg

                call s_average_cell_boundary_values(flux_src_ndqp(i, :, :))
                call s_average_cell_boundary_values(flux_gsrc_ndqp(i, :, :))
                ! ===============================================================
            end if

            if (alt_soundspeed .or. regularization) then
                do j = 0, m
                    do k = 0, n
                        do l = 0, p
                            blkmod1(j, k, l) = ((fluid_pp(1)%gamma + 1d0)*q_prim_qp(0, 0, 0)%vf(E_idx)%sf(j, k, l) + &
                                                fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                            blkmod2(j, k, l) = ((fluid_pp(2)%gamma + 1d0)*q_prim_qp(0, 0, 0)%vf(E_idx)%sf(j, k, l) + &
                                                fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                            alpha1(j, k, l) = q_cons_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j, k, l)

                            if (bubbles) then
                                alpha2(j, k, l) = q_cons_qp(0, 0, 0)%vf(alf_idx - 1)%sf(j, k, l)
                            else
                                alpha2(j, k, l) = q_cons_qp(0, 0, 0)%vf(adv_idx%end)%sf(j, k, l)
                            end if

                            Kterm(j, k, l) = alpha1(j, k, l)*alpha2(j, k, l)*(blkmod2(j, k, l) - blkmod1(j, k, l))/ &
                                             (alpha1(j, k, l)*blkmod2(j, k, l) + alpha2(j, k, l)*blkmod1(j, k, l))
                        end do
                    end do
                end do
            end if

            ! RHS Contribution in x-direction ===============================
            if (i == 1) then

                ! Compute upwind slope and flux limiter function value if TVD
                ! flux limiter is chosen
                if (tvd_rhs_flux) call s_get_tvd_flux(q_cons_vf, q_prim_vf, rhs_vf, i)

                ! Applying characteristic boundary conditions
                if (bc_x%beg <= -5) then
                    call s_cbc(q_prim_qp(0, 0, 0)%vf, flux_ndqp(i, 0, 0)%vf, &
                               flux_src_ndqp(i, 0, 0)%vf, i, -1, ix, iy, iz)
                end if

                if (bc_x%end <= -5) then
                    call s_cbc(q_prim_qp(0, 0, 0)%vf, flux_ndqp(i, 0, 0)%vf, &
                               flux_src_ndqp(i, 0, 0)%vf, i, 1, ix, iy, iz)
                end if

                ! Applying the Riemann fluxes
                do j = 1, sys_size
                    do k = 0, m
                        rhs_vf(j)%sf(k, :, :) = 1d0/dx(k)* &
                                                (flux_ndqp(i, 0, 0)%vf(j)%sf(k - 1, 0:n, 0:p) &
                                                 - flux_ndqp(i, 0, 0)%vf(j)%sf(k, 0:n, 0:p))
                    end do
                end do

                ! Applying source terms to the RHS of the advection equations
                if (riemann_solver == 1) then
                    !HLL, no K \div(u) so this just adds (subtracts?)
                    ! \alpha_i \div(u) to RHS of \alpha_i transport equation
                    do j = adv_idx%beg, adv_idx%end
                        do k = 0, m
                            rhs_vf(j)%sf(k, :, :) = &
                                rhs_vf(j)%sf(k, :, :) + 1d0/dx(k)* &
                                q_prim_qp(0, 0, 0)%vf(cont_idx%end + i)%sf(k, 0:n, 0:p)* &
                                (flux_src_ndqp(i, 0, 0)%vf(j)%sf(k - 1, 0:n, 0:p) &
                                 - flux_src_ndqp(i, 0, 0)%vf(j)%sf(k, 0:n, 0:p))

                            ! HLLC Version: different in strang ways. sign on div u seems correct here, not above
                            !     rhs_vf(j)%sf(k,:,:) = &
                            !     rhs_vf(j)%sf(k,:,:) + 1d0/dx(k) * &
                            !     q_cons_qp(0,0,0)%vf(j)%sf(k,0:n,0:p) * &
                            !     ( flux_src_ndqp(i,0,0)%vf(j)%sf( k ,0:n,0:p) &
                            !     - flux_src_ndqp(i,0,0)%vf(j)%sf(k-1,0:n,0:p) )
                        end do
                    end do
                else
                    do j = adv_idx%beg, adv_idx%end
                        if (alt_soundspeed .or. regularization) then
                            if (adv_alphan .and. (j == adv_idx%end) .and. (bubbles .neqv. .true.)) then
                                !adv_idx%end, -k div(u)
                                do k = 0, m
                                    rhs_vf(j)%sf(k, :, :) = &
                                        rhs_vf(j)%sf(k, :, :) + 1d0/dx(k)* &
                                        (q_cons_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p) - Kterm(k, :, :))* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(k, 0:n, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(k - 1, 0:n, 0:p))
                                end do
                            else if (adv_alphan .and. (j == adv_idx%beg) .and. (bubbles .neqv. .true.)) then
                                !adv_idx%beg, +k div(u)
                                do k = 0, m
                                    rhs_vf(j)%sf(k, :, :) = &
                                        rhs_vf(j)%sf(k, :, :) + 1d0/dx(k)* &
                                        (q_cons_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p) + Kterm(k, :, :))* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(k, 0:n, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(k - 1, 0:n, 0:p))
                                end do
                            else if (adv_alphan .and. (j == adv_idx%beg) .and. bubbles) then
                                !liquid part, +k div(u)
                                do k = 0, m
                                    rhs_vf(j)%sf(k, :, :) = &
                                        rhs_vf(j)%sf(k, :, :) + 1d0/dx(k)* &
                                        (q_cons_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p) + Kterm(k, :, :))* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(k, 0:n, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(k - 1, 0:n, 0:p))
                                end do
                            else if (adv_alphan .and. (j == alf_idx - 1) .and. bubbles) then
                                !resolved gas, -k div(u)
                                do k = 0, m
                                    rhs_vf(j)%sf(k, :, :) = &
                                        rhs_vf(j)%sf(k, :, :) + 1d0/dx(k)* &
                                        (q_cons_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p) - Kterm(k, :, :))* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(k, 0:n, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(k - 1, 0:n, 0:p))
                                end do
                            else if (adv_alphan .and. (j == alf_idx) .and. bubbles) then
                                !bubble part, no k div(u)
                                do k = 0, m
                                    rhs_vf(j)%sf(k, :, :) = &
                                        rhs_vf(j)%sf(k, :, :) + 1d0/dx(k)* &
                                        q_cons_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(k, 0:n, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(k - 1, 0:n, 0:p))
                                end do
                            end if
                        else
                            !no k \div u, just adds other part of the transport equation
                            do k = 0, m
                                rhs_vf(j)%sf(k, :, :) = &
                                    rhs_vf(j)%sf(k, :, :) + 1d0/dx(k)* &
                                    q_cons_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)* &
                                    (flux_src_ndqp(i, 0, 0)%vf(j)%sf(k, 0:n, 0:p) &
                                     - flux_src_ndqp(i, 0, 0)%vf(j)%sf(k - 1, 0:n, 0:p))
                            end do
                        end if
                    end do
                end if

                if (DEBUG) print *, 'pre-QBMM rhs'
                do j = 1, sys_size
                    ! DO k = 0,m
                    !     IF ( ABS(rhs_vf(j)%sf(k,0,0)) > 1.d-12 ) THEN
                    !         PRINT*, 'large RHS pre QBMM:', rhs_vf(j)%sf(k,0,0)
                    !         CALL s_mpi_abort()
                    !     END IF
                    ! END DO
                    ! print*, 'rhs = ', rhs_vf(j)%sf(1,0,0)
                end do

                if (bubbles) then
                    if (qbmm) then
                        ! advection source
                        rhs_vf(alf_idx)%sf(0:m, 0:n, 0:p) = rhs_vf(alf_idx)%sf(0:m, 0:n, 0:p) + mom_sp(2)%sf(0:m, 0:n, 0:p)
                        ! bubble sources
                        j = bub_idx%beg
                        do k = 1, nb
                            rhs_vf(j)%sf(0:m, 0:n, 0:p) = rhs_vf(j)%sf(0:m, 0:n, 0:p) + mom_3d(0, 0, k)%sf(0:m, 0:n, 0:p)
                            rhs_vf(j + 1)%sf(0:m, 0:n, 0:p) = rhs_vf(j + 1)%sf(0:m, 0:n, 0:p) + mom_3d(1, 0, k)%sf(0:m, 0:n, 0:p)
                            rhs_vf(j + 2)%sf(0:m, 0:n, 0:p) = rhs_vf(j + 2)%sf(0:m, 0:n, 0:p) + mom_3d(0, 1, k)%sf(0:m, 0:n, 0:p)
                            rhs_vf(j + 3)%sf(0:m, 0:n, 0:p) = rhs_vf(j + 3)%sf(0:m, 0:n, 0:p) + mom_3d(2, 0, k)%sf(0:m, 0:n, 0:p)
                            rhs_vf(j + 4)%sf(0:m, 0:n, 0:p) = rhs_vf(j + 4)%sf(0:m, 0:n, 0:p) + mom_3d(1, 1, k)%sf(0:m, 0:n, 0:p)
                            rhs_vf(j + 5)%sf(0:m, 0:n, 0:p) = rhs_vf(j + 5)%sf(0:m, 0:n, 0:p) + mom_3d(0, 2, k)%sf(0:m, 0:n, 0:p)
                            j = j + 6
                        end do
                    else
                        call s_get_divergence(i, q_prim_vf, divu)
                        call s_compute_bubble_source(i, q_prim_vf, q_cons_vf, divu, &
                                                     bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src)

                        rhs_vf(alf_idx)%sf(:, :, :) = rhs_vf(alf_idx)%sf(:, :, :) + bub_adv_src(:, :, :)
                        if (num_fluids > 1) rhs_vf(adv_idx%beg)%sf(:, :, :) = rhs_vf(adv_idx%beg)%sf(:, :, :) - bub_adv_src(:, :, :)

                        do k = 1, nb
                            rhs_vf(bub_idx%rs(k))%sf(:, :, :) = rhs_vf(bub_idx%rs(k))%sf(:, :, :) + bub_r_src(k, :, :, :)
                            rhs_vf(bub_idx%vs(k))%sf(:, :, :) = rhs_vf(bub_idx%vs(k))%sf(:, :, :) + bub_v_src(k, :, :, :)
                            if (polytropic .neqv. .true.) then
                                rhs_vf(bub_idx%ps(k))%sf(:, :, :) = rhs_vf(bub_idx%ps(k))%sf(:, :, :) + bub_p_src(k, :, :, :)
                                rhs_vf(bub_idx%ms(k))%sf(:, :, :) = rhs_vf(bub_idx%ms(k))%sf(:, :, :) + bub_m_src(k, :, :, :)
                            end if
                        end do
                    end if
                end if

                if (DEBUG) print *, 'after bub sources'
                do j = alf_idx, sys_size
                    ! do j = 1, sys_size
                    ! print*, 'rhs = ', rhs_vf(j)%sf(1,0,0)
                end do

                ! do j = 1,sys_size
                ! do k = 1,m
                !     IF ( ABS(rhs_vf(j)%sf(k,0,0) - rhs_vf(j)%sf(k-1,0,0)) > 1.d-14) THEN
                !         print*, 'detected discontinuity in rhs at equation ', j
                !         print*, 'rhs: ', rhs_vf(j)%sf(:,0,0)
                !         call s_mpi_abort()
                !     END IF
                ! end do
                ! end do

                if (monopole) then
                    mono_mass_src = 0d0; mono_mom_src = 0d0; mono_e_src = 0d0; 
                    do j = 1, num_mono
                        call s_get_monopole(i, q_prim_vf, t_step, mono(j))
                    end do
                    do k = cont_idx%beg, cont_idx%end
                        rhs_vf(k)%sf(:, :, :) = rhs_vf(k)%sf(:, :, :) + mono_mass_src(:, :, :)
                    end do
                    do k = mom_idx%beg, mom_idx%end
                        rhs_vf(k)%sf(:, :, :) = rhs_vf(k)%sf(:, :, :) + mono_mom_src(k - cont_idx%end, :, :, :)
                    end do
                    rhs_vf(E_idx)%sf(:, :, :) = rhs_vf(E_idx)%sf(:, :, :) + mono_e_src(:, :, :)
                end if

                ! Hypoelastic rhs terms
                if (hypoelasticity) then

                    ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

                    if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

                    ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg

                    iv%beg = mom_idx%beg; iv%end = mom_idx%end

                    call s_apply_scalar_divergence_theorem(qL_prim_ndqp(1, 0, 0)%vf(iv%beg:iv%end), &
                                                           qR_prim_ndqp(1, 0, 0)%vf(iv%beg:iv%end), &
                                                           dq_prim_dx_qp(0, 0, 0)%vf(iv%beg:iv%end), 1)

                    call s_reconstruct_cell_interior_values(dq_prim_dx_qp)

                    ! Building shear modulus and viscosity mixture variable fields (dimension m*n*p)
                    do j = 0, m
                        do k = 0, n
                            do l = 0, p
                                call s_convert_to_mixture_variables(q_prim_qp(0, 0, 0)%vf, rho_K, gamma_K, &
                                                                    pi_inf_K, Re_K, We_K, j, k, l, &
                                                                    G_K, fluid_pp(:)%G)
                                rho_K_field(j, k, l) = rho_K
                                G_K_field(j, k, l) = G_K
                            end do
                        end do
                    end do

                    do k = 0, m

                        j = stress_idx%beg

                        rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                (q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                 2.0*G_K_field(k, 0:n, 0:p)*(2.0/3.0)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p))

                        if (n > 0) then

                            call s_apply_scalar_divergence_theorem(qL_prim_ndqp(1, 0, 0)%vf(iv%beg:iv%end), &
                                                                   qR_prim_ndqp(1, 0, 0)%vf(iv%beg:iv%end), &
                                                                   dq_prim_dy_qp(0, 0, 0)%vf(iv%beg:iv%end), 2)

                            call s_reconstruct_cell_interior_values(dq_prim_dy_qp)

                            j = stress_idx%beg
                            rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                    (q_prim_qp(0, 0, 0)%vf(j + 1)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j + 1)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) - &
                                                     2.0*G_K_field(k, 0:n, 0:p)*(1.0/3.0)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p))

                            j = stress_idx%beg + 1
                            rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                    (q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j - 1)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j + 1)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                     2.0*G_K_field(k, 0:n, 0:p)*(1.0/2.0)*(dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                                                           dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p)))

                            j = stress_idx%beg + 2
                            rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                    (q_prim_qp(0, 0, 0)%vf(j - 1)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j - 1)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                     2.0*G_K_field(k, 0:n, 0:p)*(dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) - (1.0/3.0)* &
                                                                                 (dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                                                  dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p))))

                            if (p > 0) then

                                call s_apply_scalar_divergence_theorem(qL_prim_ndqp(1, 0, 0)%vf(iv%beg:iv%end), &
                                                                       qR_prim_ndqp(1, 0, 0)%vf(iv%beg:iv%end), &
                                                                       dq_prim_dz_qp(0, 0, 0)%vf(iv%beg:iv%end), 3)

                                call s_reconstruct_cell_interior_values(dq_prim_dz_qp)

                                j = stress_idx%beg
                                rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                        (q_prim_qp(0, 0, 0)%vf(j + 3)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j + 3)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         2.0*G_K_field(k, 0:n, 0:p)*(1.0/3.0)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p))

                                j = stress_idx%beg + 1
                                rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                        (q_prim_qp(0, 0, 0)%vf(j + 3)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j + 2)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p))

                                j = stress_idx%beg + 2
                                rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                        (q_prim_qp(0, 0, 0)%vf(j + 2)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j + 2)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         2.0*G_K_field(k, 0:n, 0:p)*(1.0/3.0)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p))

                                j = stress_idx%beg + 3
                                rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                        (q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j - 3)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j + 1)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j - 2)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j + 2)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) + &
                                                         2.0*G_K_field(k, 0:n, 0:p)*(1.0/2.0)*(dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                                                               dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p)))

                                j = stress_idx%beg + 4
                                rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                        (q_prim_qp(0, 0, 0)%vf(j - 1)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j - 3)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j - 2)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j + 1)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) + &
                                                         2.0*G_K_field(k, 0:n, 0:p)*(1.0/2.0)*(dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                                                               dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p)))

                                j = stress_idx%end
                                rhs_vf(j)%sf(k, :, :) = rhs_vf(j)%sf(k, :, :) + rho_K_field(k, 0:n, 0:p)* &
                                                        (q_prim_qp(0, 0, 0)%vf(j - 2)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j - 2)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j - 1)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j - 1)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) + &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - &
                                                         q_prim_qp(0, 0, 0)%vf(j)%sf(k, 0:n, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) + &
                                                         2.0*G_K_field(k, 0:n, 0:p)*(dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p) - (1.0/3.0)* &
                                                                                     (dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(k, 0:n, 0:p) + &
                                                                                      dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(k, 0:n, 0:p) + &
                                                                                      dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(k, 0:n, 0:p))))

                            end if
                        end if
                    end do
                    !DEALLOCATE(var%sf,grad_x%sf,grad_y%sf,grad_z%sf,norm%sf)
                end if

                ! Applying source terms to the RHS of the internal energy equations
                if (model_eqns == 3) then
                    do j = 1, num_fluids
                        do k = 0, m
                            rhs_vf(j + internalEnergies_idx%beg - 1)%sf(k, :, :) = &
                                rhs_vf(j + internalEnergies_idx%beg - 1)%sf(k, :, :) - 1d0/dx(k)* &
                                q_cons_qp(0, 0, 0)%vf(j + adv_idx%beg - 1)%sf(k, 0:n, 0:p)* &
                                q_prim_qp(0, 0, 0)%vf(E_idx)%sf(k, 0:n, 0:p)* &
                                (flux_src_ndqp(i, 0, 0)%vf(adv_idx%beg)%sf(k, 0:n, 0:p) - &
                                 flux_src_ndqp(i, 0, 0)%vf(adv_idx%beg)%sf(k - 1, 0:n, 0:p))
                        end do
                    end do
                end if

                ! Applying the viscous and capillary source fluxes from the Riemann solver
                if (any(Re_size > 0) .or. (We_size > 0 .and. We_riemann_flux)) then
                    do j = mom_idx%beg, E_idx
                        do k = 0, m
                            rhs_vf(j)%sf(k, :, :) = &
                                rhs_vf(j)%sf(k, :, :) + 1d0/dx(k)* &
                                (flux_src_ndqp(i, 0, 0)%vf(j)%sf(k - 1, 0:n, 0:p) &
                                 - flux_src_ndqp(i, 0, 0)%vf(j)%sf(k, 0:n, 0:p))
                        end do
                    end do
                end if

                ! Applying the capillary stress tensor calculated outside the Riemann solver
                if (We_size > 0 .and. We_rhs_flux) then
                    if (p > 0) then
                        call s_compute_capillary_stress_tensor(i, q_prim_qp(0, 0, 0)%vf, &
                                                               dq_prim_dx_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dz_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               gm_alpha_qp(0, 0, 0)%vf)
                    else
                        call s_compute_capillary_stress_tensor(i, q_prim_qp(0, 0, 0)%vf, &
                                                               dq_prim_dx_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               gm_alpha_qp(0, 0, 0)%vf)
                    end if
                    ! Divergence of capillary stress tensor computed using finite differences
                    do j = mom_idx%beg, E_idx
                        do k = 0, m
                            rhs_vf(j)%sf(k, :, :) = &
                                rhs_vf(j)%sf(k, :, :) - 1d0/(x_cc(k + 1) - x_cc(k - 1))* &
                                (tau_We_vf(j)%sf(k + 1, 0:n, 0:p) &
                                 - tau_We_vf(j)%sf(k - 1, 0:n, 0:p))
                        end do
                    end do
                end if

                ! Applying capillary source terms in momentum and energy equations
                if (We_size > 0 .and. We_src) then
                    call s_compute_capillary_source(i, q_prim_vf)
                    rhs_vf(cont_idx%end + i)%sf(:, :, :) = rhs_vf(cont_idx%end + i)%sf(:, :, :) + We_mtm_src(:, :, :)
                    rhs_vf(E_idx)%sf(:, :, :) = rhs_vf(E_idx)%sf(:, :, :) + We_nrg_src(:, :, :)
                end if

                ! ===============================================================

                ! RHS Contribution in y-direction ===============================
            elseif (i == 2) then
                if (DEBUG) print *, 'get dir 2'
                ! Compute upwind slope and flux limiter function value if TVD
                ! flux limiter is chosen

                if (tvd_rhs_flux) call s_get_tvd_flux(q_cons_vf, q_prim_vf, rhs_vf, i)

                ! Applying characteristic boundary conditions
                if (bc_y%beg <= -5 .and. bc_y%beg /= -13) then
                    call s_cbc(q_prim_qp(0, 0, 0)%vf, flux_ndqp(i, 0, 0)%vf, &
                               flux_src_ndqp(i, 0, 0)%vf, i, -1, ix, iy, iz)
                end if

                if (bc_y%end <= -5) then
                    call s_cbc(q_prim_qp(0, 0, 0)%vf, flux_ndqp(i, 0, 0)%vf, &
                               flux_src_ndqp(i, 0, 0)%vf, i, 1, ix, iy, iz)
                end if

                ! Applying the Riemann fluxes
                do j = 1, sys_size
                    do k = 0, n
                        rhs_vf(j)%sf(:, k, :) = &
                            rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                            (flux_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p) &
                             - flux_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p))
                    end do
                end do

                ! Applying source terms to the RHS of the advection equations
                if (riemann_solver == 1) then
                    do j = adv_idx%beg, adv_idx%end
                        do k = 0, n
                            rhs_vf(j)%sf(:, k, :) = &
                                rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                                q_prim_qp(0, 0, 0)%vf(cont_idx%end + i)%sf(0:m, k, 0:p)* &
                                (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p) &
                                 - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p))
                        end do
                    end do
                else
                    do j = adv_idx%beg, adv_idx%end
                        do k = 0, n
                            if (alt_soundspeed .or. regularization) then
                                if (adv_alphan .and. (j == adv_idx%beg) .and. bubbles) then
                                    !liquid part, +k div(u)
                                    rhs_vf(j)%sf(:, k, :) = &
                                        rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                                        (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p) + Kterm(:, k, :))* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p))
                                else if (adv_alphan .and. (j == alf_idx - 1) .and. bubbles) then
                                    !resolved gas, -k div(u)
                                    rhs_vf(j)%sf(:, k, :) = &
                                        rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                                        (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p) - Kterm(:, k, :))* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p))
                                else if (adv_alphan .and. (j == alf_idx) .and. bubbles) then
                                    !bubble part, no k div(u)
                                    rhs_vf(j)%sf(:, k, :) = &
                                        rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                                        q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p))
                                else if (adv_alphan .and. (j == adv_idx%end)) then
                                    rhs_vf(j)%sf(:, k, :) = &
                                        rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                                        (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p) - Kterm(:, k, :))* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p))
                                    if (cyl_coord) then
                                        rhs_vf(j)%sf(:, k, :) = &
                                            rhs_vf(j)%sf(:, k, :) - Kterm(:, k, :)/2d0/y_cc(k)* &
                                            (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p) &
                                             + flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p))
                                    end if
                                else
                                    rhs_vf(j)%sf(:, k, :) = &
                                        rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                                        (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p) + Kterm(:, k, :))* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p))
                                    if (cyl_coord) then
                                        rhs_vf(j)%sf(:, k, :) = &
                                            rhs_vf(j)%sf(:, k, :) + Kterm(:, k, :)/2d0/y_cc(k)* &
                                            (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p) &
                                             + flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p))
                                    end if
                                end if
                            else
                                rhs_vf(j)%sf(:, k, :) = &
                                    rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                                    q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)* &
                                    (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p) &
                                     - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p))
                            end if
                        end do
                    end do
                end if

                if (bubbles) then
                    call s_get_divergence(i, q_prim_vf, divu)
                    call s_compute_bubble_source(i, q_prim_vf, q_cons_vf, divu, &
                                                 bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src)

                    rhs_vf(alf_idx)%sf(:, :, :) = rhs_vf(alf_idx)%sf(:, :, :) + bub_adv_src(:, :, :)
                    if (num_fluids > 1) rhs_vf(adv_idx%beg)%sf(:, :, :) = rhs_vf(adv_idx%beg)%sf(:, :, :) - bub_adv_src(:, :, :)

                    do k = 1, nb
                        rhs_vf(bub_idx%rs(k))%sf(:, :, :) = rhs_vf(bub_idx%rs(k))%sf(:, :, :) + bub_r_src(k, :, :, :)
                        rhs_vf(bub_idx%vs(k))%sf(:, :, :) = rhs_vf(bub_idx%vs(k))%sf(:, :, :) + bub_v_src(k, :, :, :)
                        if (polytropic .neqv. .true.) then
                            rhs_vf(bub_idx%ps(k))%sf(:, :, :) = rhs_vf(bub_idx%ps(k))%sf(:, :, :) + bub_p_src(k, :, :, :)
                            rhs_vf(bub_idx%ms(k))%sf(:, :, :) = rhs_vf(bub_idx%ms(k))%sf(:, :, :) + bub_m_src(k, :, :, :)
                        end if
                    end do
                end if

                if (monopole) then
                    mono_mass_src = 0d0; mono_mom_src = 0d0; mono_e_src = 0d0; 
                    do j = 1, num_mono
                        call s_get_monopole(i, q_prim_vf, t_step, mono(j))
                    end do

                    do k = cont_idx%beg, cont_idx%end
                        rhs_vf(k)%sf(:, :, :) = rhs_vf(k)%sf(:, :, :) + mono_mass_src(:, :, :)
                    end do

                    do k = mom_idx%beg, mom_idx%end
                        rhs_vf(k)%sf(:, :, :) = rhs_vf(k)%sf(:, :, :) + mono_mom_src(k - cont_idx%end, :, :, :)
                    end do
                    rhs_vf(E_idx)%sf(:, :, :) = rhs_vf(E_idx)%sf(:, :, :) + mono_e_src(:, :, :)
                end if

                ! Hypoelastic rhs terms
                if (hypoelasticity) then
                    do k = 0, n
                        j = stress_idx%beg

                        rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                (q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                 2.0*G_K_field(0:m, k, 0:p)*(2.0/3.0)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p))

                        ! We must have n > 0 given i == 2
                        j = stress_idx%beg
                        rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                (q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) - &
                                                 2.0*G_K_field(0:m, k, 0:p)*(1.0/3.0)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p))

                        j = stress_idx%beg + 1
                        rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                (q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                 2.0*G_K_field(0:m, k, 0:p)*(1.0/2.0)*(dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                                                       dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p)))

                        j = stress_idx%beg + 2
                        rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                (q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                 2.0*G_K_field(0:m, k, 0:p)*(dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) - (1.0/3.0)* &
                                                                             (dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                                              dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p))))

                        if (p > 0) then
                            j = stress_idx%beg
                            rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                    (q_prim_qp(0, 0, 0)%vf(j + 3)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j + 3)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     2.0*G_K_field(0:m, k, 0:p)*(1.0/3.0)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p))

                            j = stress_idx%beg + 1
                            rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                    (q_prim_qp(0, 0, 0)%vf(j + 3)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j + 2)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p))

                            j = stress_idx%beg + 2
                            rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                    (q_prim_qp(0, 0, 0)%vf(j + 2)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j + 2)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     2.0*G_K_field(0:m, k, 0:p)*(1.0/3.0)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p))

                            j = stress_idx%beg + 3
                            rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                    (q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j - 3)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j - 2)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j + 2)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) + &
                                                     2.0*G_K_field(0:m, k, 0:p)*(1.0/2.0)*(dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                                                           dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p)))

                            j = stress_idx%beg + 4
                            rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                    (q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j - 3)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j - 2)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) + &
                                                     2.0*G_K_field(0:m, k, 0:p)*(1.0/2.0)*(dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                                                           dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p)))

                            j = stress_idx%end
                            rhs_vf(j)%sf(:, k, :) = rhs_vf(j)%sf(:, k, :) + rho_K_field(:, k, :)* &
                                                    (q_prim_qp(0, 0, 0)%vf(j - 2)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j - 2)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) + &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - &
                                                     q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, k, 0:p)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) + &
                                                     2.0*G_K_field(0:m, k, 0:p)*(dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p) - (1.0/3.0)* &
                                                                                 (dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, k, 0:p) + &
                                                                                  dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, k, 0:p) + &
                                                                                  dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, k, 0:p))))

                        end if
                    end do
                end if

                ! Applying source terms to the RHS of the internal energy equations
                if (model_eqns == 3) then
                    do j = 1, num_fluids
                        do k = 0, n
                            rhs_vf(j + internalEnergies_idx%beg - 1)%sf(:, k, :) = &
                                rhs_vf(j + internalEnergies_idx%beg - 1)%sf(:, k, :) - 1d0/dy(k)* &
                                q_cons_qp(0, 0, 0)%vf(j + adv_idx%beg - 1)%sf(0:m, k, 0:p)* &
                                q_prim_qp(0, 0, 0)%vf(E_idx)%sf(0:m, k, 0:p)* &
                                (flux_src_ndqp(i, 0, 0)%vf(adv_idx%beg)%sf(0:m, k, 0:p) - &
                                 flux_src_ndqp(i, 0, 0)%vf(adv_idx%beg)%sf(0:m, k - 1, 0:p))
                        end do
                    end do

                    ! Applying the additional geometrical inviscid Riemann
                    ! source fluxes for the internal energy equations
                    ! using the average of velocities at cell boundaries
                    if (cyl_coord) then
                        do j = 1, num_fluids
                            do k = 0, n
                                rhs_vf(j + internalEnergies_idx%beg - 1)%sf(:, k, :) = &
                                    rhs_vf(j + internalEnergies_idx%beg - 1)%sf(:, k, :) - 5d-1/y_cc(k)* &
                                    q_cons_qp(0, 0, 0)%vf(j + adv_idx%beg - 1)%sf(0:m, k, 0:p)* &
                                    q_prim_qp(0, 0, 0)%vf(E_idx)%sf(0:m, k, 0:p)* &
                                    (flux_src_ndqp(i, 0, 0)%vf(adv_idx%beg)%sf(0:m, k, 0:p) + &
                                     flux_src_ndqp(i, 0, 0)%vf(adv_idx%beg)%sf(0:m, k - 1, 0:p))
                            end do
                        end do
                    end if
                end if

                ! Applying the geometrical inviscid Riemann source fluxes calculated as average
                ! of values at cell boundaries
                if (cyl_coord) then
                    do j = 1, sys_size
                        do k = 0, n
                            rhs_vf(j)%sf(:, k, :) = &
                                rhs_vf(j)%sf(:, k, :) - 5d-1/y_cc(k)* &
                                (flux_gsrc_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p) &
                                 + flux_gsrc_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p))
                        end do
                    end do
                end if

                ! Applying the viscous and capillary source fluxes from the Riemann solver
                if (any(Re_size > 0) .or. (We_size > 0 .and. We_riemann_flux) .or. hypoelasticity) then
                    do j = mom_idx%beg, E_idx
                        if (cyl_coord .and. ((bc_y%beg == -2) .or. (bc_y%beg == -13))) then
                            if (p > 0) then
                                call s_compute_viscous_stress_tensor(q_prim_qp(0, 0, 0)%vf, &
                                                                     dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg:mom_idx%end), &
                                                                     dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg:mom_idx%end), &
                                                                     dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg:mom_idx%end))
                            else
                                call s_compute_viscous_stress_tensor(q_prim_qp(0, 0, 0)%vf, &
                                                                     dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg:mom_idx%end), &
                                                                     dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg:mom_idx%end), &
                                                                     dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg:mom_idx%end))
                            end if
                            do k = 1, n
                                rhs_vf(j)%sf(:, k, :) = &
                                    rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                                    (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p) &
                                     - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p))
                            end do
                            rhs_vf(j)%sf(:, 0, :) = &
                                rhs_vf(j)%sf(:, 0, :) + 1d0/(y_cc(1) - y_cc(-1))* &
                                (tau_Re_vf(j)%sf(0:m, -1, 0:p) &
                                 - tau_Re_vf(j)%sf(0:m, 1, 0:p))
                        else
                            do k = 0, n
                                rhs_vf(j)%sf(:, k, :) = &
                                    rhs_vf(j)%sf(:, k, :) + 1d0/dy(k)* &
                                    (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p) &
                                     - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p))
                            end do
                        end if
                    end do
                    ! Applying the geometrical viscous Riemann source fluxes calculated as average
                    ! of values at cell boundaries
                    if (cyl_coord) then
                        do j = mom_idx%beg, E_idx
                            if ((bc_y%beg == -2) .or. (bc_y%beg == -13)) then
                                do k = 1, n
                                    rhs_vf(j)%sf(:, k, :) = &
                                        rhs_vf(j)%sf(:, k, :) - 5d-1/y_cc(k)* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p) &
                                         + flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p))
                                end do
                                rhs_vf(j)%sf(:, 0, :) = &
                                    rhs_vf(j)%sf(:, 0, :) - 1d0/y_cc(0)* &
                                    tau_Re_vf(j)%sf(0:m, 0, 0:p)
                            else
                                do k = 0, n
                                    rhs_vf(j)%sf(:, k, :) = &
                                        rhs_vf(j)%sf(:, k, :) - 5d-1/y_cc(k)* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k - 1, 0:p) &
                                         + flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, k, 0:p))
                                end do
                            end if
                        end do
                    end if
                end if

                ! Applying the capillary stress tensor calculated outside the Riemann solver
                if (We_size > 0 .and. We_rhs_flux) then
                    if (p > 0) then
                        call s_compute_capillary_stress_tensor(i, q_prim_qp(0, 0, 0)%vf, &
                                                               dq_prim_dx_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dz_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               gm_alpha_qp(0, 0, 0)%vf)
                    else
                        call s_compute_capillary_stress_tensor(i, q_prim_qp(0, 0, 0)%vf, &
                                                               dq_prim_dx_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               gm_alpha_qp(0, 0, 0)%vf)
                    end if
                    ! Divergence of capillary stress tensor computed using finite differences
                    do j = mom_idx%beg, E_idx
                        do k = 0, n
                            rhs_vf(j)%sf(:, k, :) = &
                                rhs_vf(j)%sf(:, k, :) - 1d0/(y_cc(k + 1) - y_cc(k - 1))* &
                                (tau_We_vf(j)%sf(0:m, k + 1, 0:p) &
                                 - tau_We_vf(j)%sf(0:m, k - 1, 0:p))
                        end do
                    end do
                    ! Applying the capillary contribution to geometric source term
                    if (cyl_coord) then
                        do j = mom_idx%beg, E_idx
                            do k = 0, n
                                rhs_vf(j)%sf(:, k, :) = &
                                    rhs_vf(j)%sf(:, k, :) - tau_We_vf(j)%sf(0:m, k, 0:p)/y_cc(k)
                            end do
                        end do
                    end if
                end if

                ! Applying capillary source terms in momentum and energy equations
                if (We_size > 0 .and. We_src) then
                    call s_compute_capillary_source(i, q_prim_vf)
                    rhs_vf(cont_idx%end + i)%sf(:, :, :) = rhs_vf(cont_idx%end + i)%sf(:, :, :) + We_mtm_src(:, :, :)
                    rhs_vf(E_idx)%sf(:, :, :) = rhs_vf(E_idx)%sf(:, :, :) + We_nrg_src(:, :, :)
                end if

                ! Applying interface sharpening regularization source terms
                if (regularization .and. num_dims == 2) then
                    call s_compute_regularization_source(i, q_prim_vf)
                    do j = cont_idx%beg, adv_idx%end
                        rhs_vf(j)%sf(:, :, :) = rhs_vf(j)%sf(:, :, :) + reg_src_vf(j)%sf(:, :, :)
                    end do
                end if

                ! ===============================================================

                ! RHS Contribution in z-direction ===============================
            else
                if (DEBUG) print *, 'dir = 3'
                ! Compute upwind slope and flux limiter function value if TVD
                ! flux limiter is chosen
                if (tvd_rhs_flux) call s_get_tvd_flux(q_cons_vf, q_prim_vf, rhs_vf, i)

                ! Applying characteristic boundary conditions
                if (bc_z%beg <= -5) then
                    call s_cbc(q_prim_qp(0, 0, 0)%vf, flux_ndqp(i, 0, 0)%vf, &
                               flux_src_ndqp(i, 0, 0)%vf, i, -1, ix, iy, iz)
                end if

                if (bc_z%end <= -5) then
                    call s_cbc(q_prim_qp(0, 0, 0)%vf, flux_ndqp(i, 0, 0)%vf, &
                               flux_src_ndqp(i, 0, 0)%vf, i, 1, ix, iy, iz)
                end if

                ! Applying the Riemann fluxes
                do j = 1, sys_size
                    if (grid_geometry == 3) then
                        do l = 0, n
                            do k = 0, p
                                rhs_vf(j)%sf(:, l, k) = &
                                    rhs_vf(j)%sf(:, l, k) + 1d0/dz(k)/y_cc(l)* &
                                    (flux_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k - 1) &
                                     - flux_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k))
                            end do
                        end do
                    else
                        do k = 0, p
                            rhs_vf(j)%sf(:, :, k) = &
                                rhs_vf(j)%sf(:, :, k) + 1d0/dz(k)* &
                                (flux_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k - 1) &
                                 - flux_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k))
                        end do
                    end if
                end do

                ! Applying source terms to the RHS of the advection equations
                if (riemann_solver == 1) then
                    do j = adv_idx%beg, adv_idx%end
                        if (grid_geometry == 3) then
                            do l = 0, n
                                do k = 0, p
                                    rhs_vf(j)%sf(:, l, k) = &
                                        rhs_vf(j)%sf(:, l, k) + 1d0/dz(k)/y_cc(l)* &
                                        q_prim_qp(0, 0, 0)%vf(cont_idx%end + i)%sf(0:m, l, k)* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k - 1) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k))
                                end do
                            end do
                        else
                            do k = 0, p
                                rhs_vf(j)%sf(:, :, k) = &
                                    rhs_vf(j)%sf(:, :, k) + 1d0/dz(k)* &
                                    q_prim_qp(0, 0, 0)%vf(cont_idx%end + i)%sf(0:m, 0:n, k)* &
                                    (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k - 1) &
                                     - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k))
                            end do
                        end if
                    end do
                else
                    do j = adv_idx%beg, adv_idx%end
                        if (grid_geometry == 3) then
                            do l = 0, n
                                do k = 0, p
                                    if (alt_soundspeed .or. regularization) then
                                        if (adv_alphan .and. j == adv_idx%end) then
                                            rhs_vf(j)%sf(:, l, k) = &
                                                rhs_vf(j)%sf(:, l, k) + 1d0/dz(k)/y_cc(l)* &
                                                (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, l, k) - Kterm(:, l, k))* &
                                                (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k) &
                                                 - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k - 1))
                                        else
                                            rhs_vf(j)%sf(:, l, k) = &
                                                rhs_vf(j)%sf(:, l, k) + 1d0/dz(k)/y_cc(l)* &
                                                (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, l, k) + Kterm(:, l, k))* &
                                                (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k) &
                                                 - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k - 1))
                                        end if
                                    else
                                        rhs_vf(j)%sf(:, l, k) = &
                                            rhs_vf(j)%sf(:, l, k) + 1d0/dz(k)/y_cc(l)* &
                                            q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, l, k)* &
                                            (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k) &
                                             - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k - 1))
                                    end if
                                end do
                            end do
                        else
                            do k = 0, p
                                if (alt_soundspeed .or. regularization) then
                                    if (adv_alphan .and. (j == adv_idx%beg) .and. bubbles) then
                                        !liquid part, +k div(u)
                                        rhs_vf(j)%sf(:, :, k) = &
                                            rhs_vf(j)%sf(:, :, k) + 1d0/dz(k)* &
                                            (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k) + Kterm(:, :, k))* &
                                            (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k) &
                                             - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k - 1))

                                    else if (adv_alphan .and. (j == alf_idx - 1) .and. bubbles) then
                                        !resolved gas, -k div(u)
                                        rhs_vf(j)%sf(:, :, k) = &
                                            rhs_vf(j)%sf(:, :, k) + 1d0/dz(k)* &
                                            (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k) - Kterm(:, :, k))* &
                                            (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k) &
                                             - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k - 1))
                                    else if (adv_alphan .and. (j == alf_idx) .and. bubbles) then
                                        !bubble part, no k div(u)
                                        rhs_vf(j)%sf(:, :, k) = &
                                            rhs_vf(j)%sf(:, :, k) + 1d0/dz(k)* &
                                            q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)* &
                                            (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k) &
                                             - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k - 1))
                                    else if (adv_alphan .and. j == adv_idx%end) then
                                        rhs_vf(j)%sf(:, :, k) = &
                                            rhs_vf(j)%sf(:, :, k) + 1d0/dz(k)* &
                                            (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k) - Kterm(:, :, k))* &
                                            (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k) &
                                             - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k - 1))
                                    else
                                        rhs_vf(j)%sf(:, :, k) = &
                                            rhs_vf(j)%sf(:, :, k) + 1d0/dz(k)* &
                                            (q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k) + Kterm(:, :, k))* &
                                            (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k) &
                                             - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k - 1))
                                    end if
                                else
                                    rhs_vf(j)%sf(:, :, k) = &
                                        rhs_vf(j)%sf(:, :, k) + 1d0/dz(k)* &
                                        q_cons_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)* &
                                        (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k) &
                                         - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k - 1))
                                end if
                            end do
                        end if
                    end do
                end if

                if (bubbles) then
                    call s_get_divergence(i, q_prim_vf, divu)
                    call s_compute_bubble_source(i, q_prim_vf, q_cons_vf, divu, &
                                                 bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src)

                    rhs_vf(alf_idx)%sf(:, :, :) = rhs_vf(alf_idx)%sf(:, :, :) + bub_adv_src(:, :, :)
                    if (num_fluids > 1) rhs_vf(adv_idx%beg)%sf(:, :, :) = rhs_vf(adv_idx%beg)%sf(:, :, :) - bub_adv_src(:, :, :)

                    do k = 1, nb
                        rhs_vf(bub_idx%rs(k))%sf(:, :, :) = rhs_vf(bub_idx%rs(k))%sf(:, :, :) + bub_r_src(k, :, :, :)
                        rhs_vf(bub_idx%vs(k))%sf(:, :, :) = rhs_vf(bub_idx%vs(k))%sf(:, :, :) + bub_v_src(k, :, :, :)
                        if (polytropic .neqv. .true.) then
                            rhs_vf(bub_idx%ps(k))%sf(:, :, :) = rhs_vf(bub_idx%ps(k))%sf(:, :, :) + bub_p_src(k, :, :, :)
                            rhs_vf(bub_idx%ms(k))%sf(:, :, :) = rhs_vf(bub_idx%ms(k))%sf(:, :, :) + bub_m_src(k, :, :, :)
                        end if
                    end do
                end if

                if (monopole) then
                    mono_mass_src = 0d0; mono_mom_src = 0d0; mono_e_src = 0d0; 
                    do j = 1, num_mono
                        call s_get_monopole(i, q_prim_vf, t_step, mono(j))
                    end do
                    do k = cont_idx%beg, cont_idx%end
                        rhs_vf(k)%sf(:, :, :) = rhs_vf(k)%sf(:, :, :) + mono_mass_src(:, :, :)
                    end do
                    do k = mom_idx%beg, mom_idx%end
                        rhs_vf(k)%sf(:, :, :) = rhs_vf(k)%sf(:, :, :) + mono_mom_src(k - cont_idx%end, :, :, :)
                    end do
                    rhs_vf(E_idx)%sf(:, :, :) = rhs_vf(E_idx)%sf(:, :, :) + mono_e_src(:, :, :)
                end if

                ! Hypoelastic rhs terms
                if (hypoelasticity) then
                    do k = 0, p
                        j = stress_idx%beg
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 2.0*G_K_field(0:m, 0:n, k)*(2.0/3.0)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k))

                        ! We must have n > 0
                        j = stress_idx%beg
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) - &
                                                 2.0*G_K_field(0:m, 0:n, k)*(1.0/3.0)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k))

                        j = stress_idx%beg + 1
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 2.0*G_K_field(0:m, 0:n, k)*(1.0/2.0)*(dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                                                       dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k)))

                        j = stress_idx%beg + 2
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 2.0*G_K_field(0:m, 0:n, k)*(dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) - (1.0/3.0)* &
                                                                             (dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                                              dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k))))

                        ! We must have p > 0
                        j = stress_idx%beg
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j + 3)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 3)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 2.0*G_K_field(0:m, 0:n, k)*(1.0/3.0)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k))

                        j = stress_idx%beg + 1
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j + 3)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 2)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k))

                        j = stress_idx%beg + 2
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j + 2)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 2)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 2.0*G_K_field(0:m, 0:n, k)*(1.0/3.0)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k))

                        j = stress_idx%beg + 3
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 3)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 2)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 2)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) + &
                                                 2.0*G_K_field(0:m, 0:n, k)*(1.0/2.0)*(dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                                                       dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k)))

                        j = stress_idx%beg + 4
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 3)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 2)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j + 1)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) + &
                                                 2.0*G_K_field(0:m, 0:n, k)*(1.0/2.0)*(dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                                                       dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k)))

                        j = stress_idx%end
                        rhs_vf(j)%sf(:, :, k) = rhs_vf(j)%sf(:, :, k) + rho_K_field(:, :, k)* &
                                                (q_prim_qp(0, 0, 0)%vf(j - 2)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 2)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j - 1)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) + &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - &
                                                 q_prim_qp(0, 0, 0)%vf(j)%sf(0:m, 0:n, k)*dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) + &
                                                 2.0*G_K_field(0:m, 0:n, k)*(dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k) - (1.0/3.0)* &
                                                                             (dq_prim_dx_qp(0, 0, 0)%vf(mom_idx%beg)%sf(0:m, 0:n, k) + &
                                                                              dq_prim_dy_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k) + &
                                                                              dq_prim_dz_qp(0, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k))))

                    end do
                end if

                ! Applying source terms to the RHS of the internal energy equations
                if (model_eqns == 3) then
                    do j = 1, num_fluids
                        do k = 0, p
                            rhs_vf(j + internalEnergies_idx%beg - 1)%sf(:, :, k) = &
                                rhs_vf(j + internalEnergies_idx%beg - 1)%sf(:, :, k) - 1d0/dz(k)* &
                                q_cons_qp(0, 0, 0)%vf(j + adv_idx%beg - 1)%sf(0:m, 0:n, k)* &
                                q_prim_qp(0, 0, 0)%vf(E_idx)%sf(0:m, 0:n, k)* &
                                (flux_src_ndqp(i, 0, 0)%vf(adv_idx%beg)%sf(0:m, 0:n, k) - &
                                 flux_src_ndqp(i, 0, 0)%vf(adv_idx%beg)%sf(0:m, 0:n, k - 1))
                        end do
                    end do
                end if

                ! Applying the geometrical inviscid Riemann source fluxes calculated as average
                ! of values at cell boundaries
                if (grid_geometry == 3) then
                    do j = 1, sys_size
                        do l = 0, n
                            do k = 0, p
                                rhs_vf(j)%sf(:, l, k) = &
                                    rhs_vf(j)%sf(:, l, k) - 5d-1/y_cc(l)* &
                                    (flux_gsrc_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k - 1) &
                                     + flux_gsrc_ndqp(i, 0, 0)%vf(j)%sf(0:m, l, k))
                            end do
                        end do
                    end do
                end if

                ! Applying the viscous and capillary source fluxes from the Riemann solver
                if (any(Re_size > 0) .or. (We_size > 0 .and. We_riemann_flux) .or. hypoelasticity) then
                    do j = mom_idx%beg, E_idx
                        do k = 0, p
                            rhs_vf(j)%sf(:, :, k) = &
                                rhs_vf(j)%sf(:, :, k) + 1d0/dz(k)* &
                                (flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k - 1) &
                                 - flux_src_ndqp(i, 0, 0)%vf(j)%sf(0:m, 0:n, k))
                        end do
                    end do
                    ! Modifying momentum components of geometric source term
                    if (grid_geometry == 3) then
                        do k = 0, p
                            rhs_vf(mom_idx%beg + 1)%sf(:, :, k) = &
                                rhs_vf(mom_idx%beg + 1)%sf(:, :, k) + 5d-1* &
                                (flux_src_ndqp(i, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k - 1) &
                                 + flux_src_ndqp(i, 0, 0)%vf(mom_idx%end)%sf(0:m, 0:n, k))

                            rhs_vf(mom_idx%end)%sf(:, :, k) = &
                                rhs_vf(mom_idx%end)%sf(:, :, k) - 5d-1* &
                                (flux_src_ndqp(i, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k - 1) &
                                 + flux_src_ndqp(i, 0, 0)%vf(mom_idx%beg + 1)%sf(0:m, 0:n, k))
                        end do
                    end if
                end if

                ! Applying the capillary stress tensor calculated outside the Riemann solver
                if (We_size > 0 .and. We_rhs_flux) then
                    if (p > 0) then
                        call s_compute_capillary_stress_tensor(i, q_prim_qp(0, 0, 0)%vf, &
                                                               dq_prim_dx_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dz_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               gm_alpha_qp(0, 0, 0)%vf)
                    else
                        call s_compute_capillary_stress_tensor(i, q_prim_qp(0, 0, 0)%vf, &
                                                               dq_prim_dx_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg:adv_idx%end), &
                                                               gm_alpha_qp(0, 0, 0)%vf)
                    end if
                    ! Divergence of capillary stress tensor computed using finite differences
                    do j = mom_idx%beg, E_idx
                        if (grid_geometry == 3) then
                            do l = 0, n
                                do k = 0, p
                                    rhs_vf(j)%sf(:, l, k) = &
                                        rhs_vf(j)%sf(:, l, k) - 1d0/(z_cc(k + 1) - z_cc(k - 1))/y_cc(l)* &
                                        (tau_We_vf(j)%sf(0:m, l, k + 1) &
                                         - tau_We_vf(j)%sf(0:m, l, k - 1))
                                end do
                            end do
                        else
                            do k = 0, p
                                rhs_vf(j)%sf(:, :, k) = &
                                    rhs_vf(j)%sf(:, :, k) - 1d0/(z_cc(k + 1) - z_cc(k - 1))* &
                                    (tau_We_vf(j)%sf(0:m, 0:n, k + 1) &
                                     - tau_We_vf(j)%sf(0:m, 0:n, k - 1))
                            end do
                        end if
                    end do
                end if

                ! Applying capillary source terms in momentum and energy equations
                if (We_size > 0 .and. We_src) then
                    call s_compute_capillary_source(i, q_prim_vf)
                    rhs_vf(cont_idx%end + i)%sf(:, :, :) = rhs_vf(cont_idx%end + i)%sf(:, :, :) + We_mtm_src(:, :, :)
                    rhs_vf(E_idx)%sf(:, :, :) = rhs_vf(E_idx)%sf(:, :, :) + We_nrg_src(:, :, :)
                end if

                ! Applying interface sharpening regularization source terms
                if (regularization .and. num_dims == 3) then
                    call s_compute_regularization_source(i, q_prim_vf)
                    do j = cont_idx%beg, adv_idx%end
                        rhs_vf(j)%sf(:, :, :) = rhs_vf(j)%sf(:, :, :) + reg_src_vf(j)%sf(:, :, :)
                    end do
                end if

            end if
            ! ===============================================================

        end do
        ! END: Dimensional Splitting Loop ==================================

        ! Disassociation of Working Variables ==============================
        do i = 1, sys_size
            nullify (q_cons_qp(0, 0, 0)%vf(i)%sf, q_prim_qp(0, 0, 0)%vf(i)%sf)
        end do
        ! ==================================================================

    end subroutine s_compute_rhs ! -----------------------------------------

    !>  This subroutine takes the consecutive changes in volume fraction
        !!  at a cell boundary and computes the desired slope and flux limiter
        !!  function value at the cell boundary
        !!  @param top The top flux limiter
        !!  @param bottom The bottom flux limiter
        !!  @param flux_lim_func The flux limiter
    subroutine s_compute_flux_lim(top, bottom, flux_lim_func) ! ---------------

        real(kind(0d0)), intent(INOUT) :: top, bottom
        real(kind(0d0)) :: slope
        real(kind(0d0)), intent(OUT) :: flux_lim_func

        ! Limit the flux limiter to be applied only where the change in
        ! volume fraction is greater than machine precision so that insignificant
        ! fluctuations do not trip the limiter
        if (abs(top) < 1d-8) top = 0d0
        if (abs(bottom) < 1d-8) bottom = 0d0
        ! If top = bottom, then cell is in a smooth region of the flow
        ! and the high order flux should be used. Also ensures that areas
        ! of no change in volume fraction (0/0) use the high order flux
        if (top == bottom) then
            slope = 1d0
        else
            slope = (top*bottom)/max(bottom**2d0, sgm_eps)
        end if

        ! Flux limiter function
        if (flux_lim == 1) then ! MINMOD (MM)
            flux_lim_func = max(0d0, min(1d0, slope))
        elseif (flux_lim == 2) then ! MUSCL (MC)
            flux_lim_func = max(0d0, min(2d0*slope, 5d-1*(1d0 + slope), 2d0))
        elseif (flux_lim == 3) then ! OSPRE (OP)
            flux_lim_func = (15d-1*(slope**2d0 + slope))/(slope**2d0 + slope + 1d0)
        elseif (flux_lim == 4) then ! SUPERBEE (SB)
            flux_lim_func = max(0d0, min(1d0, 2d0*slope), min(slope, 2d0))
        elseif (flux_lim == 5) then ! SWEBY (SW) (beta = 1.5)
            flux_lim_func = max(0d0, min(15d-1*slope, 1d0), min(slope, 15d-1))
        elseif (flux_lim == 6) then ! VAN ALBADA (VA)
            flux_lim_func = (slope**2d0 + slope)/(slope**2d0 + 1d0)
        elseif (flux_lim == 7) then ! VAN LEER (VL)
            flux_lim_func = (abs(slope) + slope)/(1d0 + abs(slope))
        end if
    end subroutine s_compute_flux_lim ! ------------------------------------

    !> The purpose of this subroutine is to compute the viscous
        !!      stress tensor for the cells directly next to the axis in
        !!      cylindrical coordinates. This is necessary to avoid the
        !!      1/r singularity that arises at the cell boundary coinciding
        !!      with the axis, i.e., y_cb(-1) = 0.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param grad_x_vf Cell-average primitive variable derivatives, x-dir
        !!  @param grad_y_vf Cell-average primitive variable derivatives, y-dir
        !!  @param grad_z_vf Cell-average primitive variable derivatives, z-dir
    subroutine s_compute_viscous_stress_tensor(q_prim_vf, grad_x_vf, grad_y_vf, grad_z_vf) ! ---

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(num_dims), intent(IN) :: grad_x_vf, grad_y_vf, grad_z_vf

        real(kind(0d0)) :: rho_visc, gamma_visc, pi_inf_visc  !< Mixture variables
        real(kind(0d0)), dimension(2) :: Re_visc
        real(kind(0d0)), dimension(1:num_fluids, 1:num_fluids) :: We_visc

        real(kind(0d0)), dimension(num_dims, num_dims) :: tau_Re !< Capillary stress tensor components

        type(bounds_info) :: ix, iy, iz

        integer :: i, j, k, l !< Generic loop iterator

        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size
        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg

        do i = mom_idx%beg, E_idx
            tau_Re_vf(i)%sf = 0d0
        end do

        if (Re_size(1) > 0) then    ! Shear stresses
            do l = iz%beg, iz%end
                do k = -1, 1
                    do j = ix%beg, ix%end

                        call s_convert_to_mixture_variables(q_prim_vf, rho_visc, &
                                                            gamma_visc, pi_inf_visc, &
                                                            Re_visc, We_visc, j, k, l)

                        tau_Re(2, 1) = (grad_y_vf(1)%sf(j, k, l) + &
                                        grad_x_vf(2)%sf(j, k, l))/ &
                                       Re_visc(1)

                        tau_Re(2, 2) = (4d0*grad_y_vf(2)%sf(j, k, l) &
                                        - 2d0*grad_x_vf(1)%sf(j, k, l) &
                                        - 2d0*q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)/y_cc(k))/ &
                                       (3d0*Re_visc(1))

                        do i = 1, 2
                            tau_Re_vf(cont_idx%end + i)%sf(j, k, l) = &
                                tau_Re_vf(cont_idx%end + i)%sf(j, k, l) - &
                                tau_Re(2, i)

                            tau_Re_vf(E_idx)%sf(j, k, l) = &
                                tau_Re_vf(E_idx)%sf(j, k, l) - &
                                q_prim_vf(cont_idx%end + i)%sf(j, k, l)*tau_Re(2, i)
                        end do

                    end do
                end do
            end do
        end if

        if (Re_size(2) > 0) then    ! Bulk stresses
            do l = iz%beg, iz%end
                do k = -1, 1
                    do j = ix%beg, ix%end

                        call s_convert_to_mixture_variables(q_prim_vf, rho_visc, &
                                                            gamma_visc, pi_inf_visc, &
                                                            Re_visc, We_visc, j, k, l)

                        tau_Re(2, 2) = (grad_x_vf(1)%sf(j, k, l) + &
                                        grad_y_vf(2)%sf(j, k, l) + &
                                        q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)/y_cc(k))/ &
                                       Re_visc(2)

                        tau_Re_vf(mom_idx%beg + 1)%sf(j, k, l) = &
                            tau_Re_vf(mom_idx%beg + 1)%sf(j, k, l) - &
                            tau_Re(2, 2)

                        tau_Re_vf(E_idx)%sf(j, k, l) = &
                            tau_Re_vf(E_idx)%sf(j, k, l) - &
                            q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*tau_Re(2, 2)

                    end do
                end do
            end do
        end if

        if (p == 0) return

        if (Re_size(1) > 0) then    ! Shear stresses
            do l = iz%beg, iz%end
                do k = -1, 1
                    do j = ix%beg, ix%end

                        call s_convert_to_mixture_variables(q_prim_vf, rho_visc, &
                                                            gamma_visc, pi_inf_visc, &
                                                            Re_visc, We_visc, j, k, l)

                        tau_Re(2, 2) = -(2d0/3d0)*grad_z_vf(3)%sf(j, k, l)/y_cc(k)/ &
                                       Re_visc(1)

                        tau_Re(2, 3) = ((grad_z_vf(2)%sf(j, k, l) - &
                                         q_prim_vf(mom_idx%end)%sf(j, k, l))/ &
                                        y_cc(k) + grad_y_vf(3)%sf(j, k, l))/ &
                                       Re_visc(1)

                        do i = 2, 3
                            tau_Re_vf(cont_idx%end + i)%sf(j, k, l) = &
                                tau_Re_vf(cont_idx%end + i)%sf(j, k, l) - &
                                tau_Re(2, i)

                            tau_Re_vf(E_idx)%sf(j, k, l) = &
                                tau_Re_vf(E_idx)%sf(j, k, l) - &
                                q_prim_vf(cont_idx%end + i)%sf(j, k, l)*tau_Re(2, i)
                        end do

                    end do
                end do
            end do
        end if

        if (Re_size(2) > 0) then    ! Bulk stresses
            do l = iz%beg, iz%end
                do k = -1, 1
                    do j = ix%beg, ix%end

                        tau_Re(2, 2) = grad_z_vf(3)%sf(j, k, l)/y_cc(k)/ &
                                       Re_visc(2)

                        tau_Re_vf(mom_idx%beg + 1)%sf(j, k, l) = &
                            tau_Re_vf(mom_idx%beg + 1)%sf(j, k, l) - &
                            tau_Re(2, 2)

                        tau_Re_vf(E_idx)%sf(j, k, l) = &
                            tau_Re_vf(E_idx)%sf(j, k, l) - &
                            q_prim_vf(mom_idx%beg + 1)%sf(j, k, l)*tau_Re(2, 2)

                    end do
                end do
            end do
        end if

    end subroutine s_compute_viscous_stress_tensor ! ----------------------------------------

    !>  The purpose of this procedure is to calculate the capillary
        !!      stress tensor at the cell centers in the non-conservative
        !!      formulation of the equations.
        !!  @param norm_dir Dimensional split index
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param grad_x_vf Cell-average primitive variable derivatives, x-dir
        !!  @param grad_y_vf Cell-average primitive variable derivatives, y-dir
        !!  @param grad_z_vf Cell-average primitive variable derivatives, z-dir
        !!  @param norm_vf Normal-dir derivative
    subroutine s_compute_capillary_stress_tensor(norm_dir, q_prim_vf, grad_x_vf, grad_y_vf, grad_z_vf, norm_vf) ! -------

        integer, intent(IN) :: norm_dir
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(num_fluids), intent(IN) :: norm_vf, grad_x_vf, grad_y_vf, grad_z_vf

        real(kind(0d0)) :: rho_cap, gamma_cap, pi_inf_cap  !< Mixture variables
        real(kind(0d0)), dimension(2) :: Re_cap
        real(kind(0d0)), dimension(1:num_fluids, 1:num_fluids) :: We_cap

        real(kind(0d0)), dimension(num_dims, num_dims) :: tau_We !< Capillary stress tensor components

        type(bounds_info) :: ix, iy, iz

        integer :: i, j, k, l, r !< Generic loop iterator

        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size
        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg

        do i = mom_idx%beg, E_idx
            tau_We_vf(i)%sf = 0d0
        end do

        ! Capillary stresses in x-direction =============================
        if (norm_dir == 1) then
            do i = 1, We_size
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            call s_convert_to_mixture_variables(q_prim_vf, rho_cap, &
                                                                gamma_cap, pi_inf_cap, &
                                                                Re_cap, We_cap, j, k, l)

                            tau_We(1, 1) = 1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           (norm_vf(We_idx(i, 2))%sf(j, k, l) - &
                                            grad_x_vf(We_idx(i, 2))%sf(j, k, l)* &
                                            grad_x_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                            max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps))

                            tau_We(1, 2) = -1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           grad_x_vf(We_idx(i, 2))%sf(j, k, l)* &
                                           grad_y_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                           max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps)

                            do r = 1, 2
                                tau_We_vf(cont_idx%end + r)%sf(j, k, l) = &
                                    tau_We_vf(cont_idx%end + r)%sf(j, k, l) - &
                                    tau_We(1, r)
                            end do

                            tau_We(1, 1) = -1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           grad_x_vf(We_idx(i, 2))%sf(j, k, l)* &
                                           grad_x_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                           max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps)

                            do r = 1, 2
                                tau_We_vf(E_idx)%sf(j, k, l) = &
                                    tau_We_vf(E_idx)%sf(j, k, l) - &
                                    q_prim_vf(cont_idx%end + r)%sf(j, k, l)*tau_We(1, r)
                            end do

                        end do
                    end do
                end do
            end do

            if (p > 0) then
                do i = 1, We_size
                    do l = iz%beg, iz%end
                        do k = iy%beg, iy%end
                            do j = ix%beg, ix%end

                                call s_convert_to_mixture_variables(q_prim_vf, rho_cap, &
                                                                    gamma_cap, pi_inf_cap, &
                                                                    Re_cap, We_cap, j, k, l)

                                tau_We(1, 3) = -1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                               grad_x_vf(We_idx(i, 2))%sf(j, k, l)* &
                                               grad_z_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                               max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps)

                                tau_We_vf(mom_idx%end)%sf(j, k, l) = &
                                    tau_We_vf(mom_idx%end)%sf(j, k, l) - &
                                    tau_We(1, 3)

                                tau_We_vf(E_idx)%sf(j, k, l) = &
                                    tau_We_vf(E_idx)%sf(j, k, l) - &
                                    q_prim_vf(mom_idx%end)%sf(j, k, l)*tau_We(1, 3)

                            end do
                        end do
                    end do
                end do
            end if

            ! End capillary stresses in x-direction =============================

            ! Capillary stresses in y-direction =============================
        elseif (norm_dir == 2) then
            do i = 1, We_size
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            call s_convert_to_mixture_variables(q_prim_vf, rho_cap, &
                                                                gamma_cap, pi_inf_cap, &
                                                                Re_cap, We_cap, j, k, l)

                            tau_We(2, 1) = -1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           grad_y_vf(We_idx(i, 2))%sf(j, k, l)* &
                                           grad_x_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                           max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps)

                            tau_We(2, 2) = 1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           (norm_vf(We_idx(i, 2))%sf(j, k, l) - &
                                            grad_y_vf(We_idx(i, 2))%sf(j, k, l)* &
                                            grad_y_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                            max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps))

                            do r = 1, 2
                                tau_We_vf(cont_idx%end + r)%sf(j, k, l) = &
                                    tau_We_vf(cont_idx%end + r)%sf(j, k, l) - &
                                    tau_We(2, r)
                            end do

                            tau_We(2, 2) = -1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           grad_y_vf(We_idx(i, 2))%sf(j, k, l)* &
                                           grad_y_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                           max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps)

                            do r = 1, 2
                                tau_We_vf(E_idx)%sf(j, k, l) = &
                                    tau_We_vf(E_idx)%sf(j, k, l) - &
                                    q_prim_vf(cont_idx%end + r)%sf(j, k, l)*tau_We(2, r)
                            end do

                        end do
                    end do
                end do
            end do

            if (p > 0) then
                do i = 1, We_size
                    do l = iz%beg, iz%end
                        do k = iy%beg, iy%end
                            do j = ix%beg, ix%end

                                call s_convert_to_mixture_variables(q_prim_vf, rho_cap, &
                                                                    gamma_cap, pi_inf_cap, &
                                                                    Re_cap, We_cap, j, k, l)

                                tau_We(2, 3) = -1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                               grad_y_vf(We_idx(i, 2))%sf(j, k, l)* &
                                               grad_z_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                               max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps)

                                tau_We_vf(mom_idx%end)%sf(j, k, l) = &
                                    tau_We_vf(mom_idx%end)%sf(j, k, l) - &
                                    tau_We(2, 3)

                                tau_We_vf(E_idx)%sf(j, k, l) = &
                                    tau_We_vf(E_idx)%sf(j, k, l) - &
                                    q_prim_vf(mom_idx%end)%sf(j, k, l)*tau_We(2, 3)

                            end do
                        end do
                    end do
                end do
            end if

            ! End capillary stresses in y-direction =============================

            ! Capillary stresses in z-direction =============================
        else
            do i = 1, We_size
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            tau_We(3, 1) = -1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           grad_z_vf(We_idx(i, 2))%sf(j, k, l)* &
                                           grad_x_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                           max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps)

                            tau_We(3, 2) = -1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           grad_z_vf(We_idx(i, 2))%sf(j, k, l)* &
                                           grad_y_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                           max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps)

                            tau_We(3, 3) = 1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           (norm_vf(We_idx(i, 2))%sf(j, k, l) - &
                                            grad_z_vf(We_idx(i, 2))%sf(j, k, l)* &
                                            grad_z_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                            max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps))

                            do r = 1, 3
                                tau_We_vf(cont_idx%end + r)%sf(j, k, l) = &
                                    tau_We_vf(cont_idx%end + r)%sf(j, k, l) - &
                                    tau_We(3, r)
                            end do

                            tau_We(3, 3) = -1d0/We_cap(We_idx(i, 1), We_idx(i, 2))* &
                                           grad_z_vf(We_idx(i, 2))%sf(j, k, l)* &
                                           grad_z_vf(We_idx(i, 2))%sf(j, k, l)/ &
                                           max(norm_vf(We_idx(i, 2))%sf(j, k, l), sgm_eps)

                            do r = 1, 3
                                tau_We_vf(E_idx)%sf(j, k, l) = &
                                    tau_We_vf(E_idx)%sf(j, k, l) - &
                                    q_prim_vf(cont_idx%end + r)%sf(j, k, l)*tau_We(3, r)
                            end do

                        end do
                    end do
                end do
            end do

        end if
        ! End capillary stresses in z-direction =============================

    end subroutine s_compute_capillary_stress_tensor ! ---------------------------------

    !>  The purpose of this procedure is to compute the source term
        !!      that accounts for capillary effects in the momentum equations
        !!  @param i Dimensional split index
        !!  @param q_prim_vf Cell-average primitive variables
    subroutine s_compute_capillary_source(i, q_prim_vf) ! ------------------------------

        integer, intent(IN) :: i
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        real(kind(0d0)) :: rho_cap, gamma_cap, pi_inf_cap !< Mixture variables
        real(kind(0d0)), dimension(2) :: Re_cap
        real(kind(0d0)), dimension(1:num_fluids, 1:num_fluids) :: We_cap
        real(kind(0d0)), dimension(1:num_fluids) :: kappa_cap

        real(kind(0d0)) :: F_sv

        integer :: j, k, l, r !< Generic loop iterators

        do j = 0, m
            do k = 0, n
                do l = 0, p
                    do r = 1, crv_size
                        if (alt_crv .neqv. .true.) then
                            kappa_cap(crv_idx(r)) = kappa_vf(crv_idx(r))%sf(j, k, l)
                        else
                            kappa_cap(crv_idx(r)) = alt_kappa_vf(crv_idx(r))%sf(j, k, l)
                        end if
                    end do
                    call s_convert_to_mixture_variables(q_prim_vf, rho_cap, &
                                                        gamma_cap, pi_inf_cap, &
                                                        Re_cap, We_cap, j, k, l)
                    We_mtm_src(j, k, l) = 0d0
                    We_nrg_src(j, k, l) = 0d0

                    do r = 1, We_size
                        if (i == 1) then
                            F_sv = -1d0/We_cap(We_idx(r, 1), We_idx(r, 2))* &
                                   kappa_cap(We_idx(r, 2))* &
                                   dq_prim_dx_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j, k, l)
                        elseif (i == 2) then
                            F_sv = -1d0/We_cap(We_idx(r, 1), We_idx(r, 2))* &
                                   kappa_cap(We_idx(r, 2))* &
                                   dq_prim_dy_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j, k, l)
                        elseif (i == 3) then
                            F_sv = -1d0/We_cap(We_idx(r, 1), We_idx(r, 2))* &
                                   kappa_cap(We_idx(r, 2))* &
                                   dq_prim_dz_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j, k, l)
                        end if

                        We_mtm_src(j, k, l) = We_mtm_src(j, k, l) + F_sv
                        We_nrg_src(j, k, l) = We_nrg_src(j, k, l) + F_sv* &
                                              q_prim_vf(cont_idx%end + i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
    end subroutine s_compute_capillary_source ! ----------------------------------

    !> Gets the divergence term for k div(U)
    !> @param idir Coordinate direction
    !> @param q_prim_vf Primitive variables
    !> @param mydivu Output divergence term div(U)
    subroutine s_get_divergence(idir, q_prim_vf, mydivu)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), intent(inout) :: mydivu
        integer, intent(IN) :: idir
        integer :: j, k, l !< Generic loop iterators

        !contribute to divergence computation \div(u)
        if (idir == 1) mydivu%sf(:, :, :) = 0d0

        do j = 0, m
            do k = 0, n
                do l = 0, p
                    if (idir == 1) then
                        mydivu%sf(j, k, l) = 0.5d0/dx(j)*(q_prim_vf(cont_idx%end + idir)%sf(j + 1, k, l) - &
                                                          q_prim_vf(cont_idx%end + idir)%sf(j - 1, k, l))
                    else if (idir == 2) then
                        mydivu%sf(j, k, l) = mydivu%sf(j, k, l) + 0.5d0/dy(k)*(q_prim_vf(cont_idx%end + idir)%sf(j, k + 1, l) - &
                                                                               q_prim_vf(cont_idx%end + idir)%sf(j, k - 1, l))
                    else if (idir == 3) then
                        mydivu%sf(j, k, l) = mydivu%sf(j, k, l) + 0.5d0/dz(l)*(q_prim_vf(cont_idx%end + idir)%sf(j, k, l + 1) - &
                                                                               q_prim_vf(cont_idx%end + idir)%sf(j, k, l - 1))
                    end if
                end do
            end do
        end do

    end subroutine s_get_divergence

    !> The purpose of this procedure is to compute the source term
        !! that are needed for generating one-way acoustic waves
        !! @param idir Coordinate direction
        !! @param q_prim_vf Primitive variables
        !! @param t_step Current time-step
        !! @param mymono Monopole parameters
    subroutine s_get_monopole(idir, q_prim_vf, t_step, mymono) ! ------------------------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(mono_parameters), intent(IN) :: mymono
        integer, intent(IN) :: idir, t_step

        integer :: ndirs, j, k, l

        real(kind(0d0)) :: mytime, sound, n_tait, B_tait
        real(kind(0d0)) :: s2, myRho, const_sos

        real(kind(0d0)), dimension(2) :: Re
        real(kind(0d0)), dimension(num_fluids, &
                                   num_fluids) :: We

        ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3

        if (idir == ndirs) then
            mytime = t_step*dt
            if (proc_rank == 0) print *, 'time', mytime, 'delay', mymono%delay, dflt_real
            if ((mytime < mymono%delay) .and. mymono%delay /= dflt_real) return

            do j = 0, m; do k = 0, n; do l = 0, p
                    call s_convert_to_mixture_variables(q_prim_vf, myRho, n_tait, B_tait, Re, We, j, k, l)
                    n_tait = 1.d0/n_tait + 1.d0 !make this the usual little 'gamma'

                    sound = n_tait*(q_prim_vf(E_idx)%sf(j, k, l) + ((n_tait - 1d0)/n_tait)*B_tait)/myRho
                    sound = dsqrt(sound)

                    const_sos = dsqrt(n_tait)

                    s2 = f_g(mytime, sound, const_sos, mymono)*f_delta(j, k, l, mymono%loc, mymono%length, mymono)

                    mono_mass_src(j, k, l) = mono_mass_src(j, k, l) + s2/sound
                    if (n == 0) then

                        ! 1D
                        if (mymono%dir < -0.1d0) then
                            !left-going wave
                            mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) - s2
                        else
                            !right-going wave
                            mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2
                        end if
                    else if (p == 0) then
                        ! IF ( (j==1) .AND. (k==1) .AND. proc_rank == 0) &
                        !    PRINT*, '====== Monopole magnitude: ', f_g(mytime,sound,const_sos,mymono)

                        if (mymono%dir .ne. dflt_real) then
                            ! 2d
                            !mono_mom_src(1,j,k,l) = s2
                            !mono_mom_src(2,j,k,l) = s2
                            mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(mymono%dir)
                            mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*sin(mymono%dir)
                        end if
                    else
                        ! 3D
                        if (mymono%dir .ne. dflt_real) then
                            mono_mom_src(1, j, k, l) = mono_mom_src(1, j, k, l) + s2*cos(mymono%dir)
                            mono_mom_src(2, j, k, l) = mono_mom_src(2, j, k, l) + s2*sin(mymono%dir)
                        end if
                    end if

                    if (model_eqns .ne. 4) then
                        mono_E_src(j, k, l) = mono_E_src(j, k, l) + s2*sound/(n_tait - 1.d0)
                    end if
                end do; end do; end do
        end if

    end subroutine s_get_monopole

    !> This function gives the temporally varying amplitude of the pulse
        !! @param mytime Simulation time
        !! @param sos Sound speed
        !! @param mysos Alternative speed of sound for testing
        !! @param mymono Monopole parameterrs
    function f_g(mytime, sos, mysos, mymono)

        real(kind(0d0)), intent(IN) :: mytime, sos, mysos
        type(mono_parameters), intent(IN) :: mymono
        real(kind(0d0)) :: period, t0, sigt, pa
        real(kind(0d0)) :: offset
        real(kind(0d0)) :: f_g

        offset = 0d0
        if (mymono%delay /= dflt_real) offset = mymono%delay

        if (mymono%pulse == 1) then
            ! Sine wave
            period = mymono%length/sos
            f_g = 0d0
            if (mytime <= (mymono%npulse*period + offset)) then
                f_g = mymono%mag*sin((mytime + offset)*2.d0*pi/period)
            end if
        else if (mymono%pulse == 2) then
            ! Gaussian pulse
            sigt = mymono%length/sos/7.d0
            t0 = 3.5d0*sigt
            f_g = mymono%mag/(dsqrt(2.d0*pi)*sigt)* &
                  dexp(-0.5d0*((mytime - t0)**2.d0)/(sigt**2.d0))
        else if (mymono%pulse == 3) then
            ! Square wave
            sigt = mymono%length/sos
            t0 = 0d0; f_g = 0d0
            if (mytime > t0 .and. mytime < sigt) then
                f_g = mymono%mag
            end if
        else
            print '(A)', 'No pulse type detected. Exiting ...'
            call s_mpi_abort()
        end if

    end function f_g

    !> This function give the spatial support of the acoustic source
        !! @param j First coordinate-direction location index
        !! @param k Second coordinate-direction location index
        !! @param l Third coordinate-direction location index
        !! @param mono_loc Nominal source term location
        !! @param mono_leng Length of source term in space
        !! @param mymono Monopole parameters
    function f_delta(j, k, l, mono_loc, mono_leng, mymono)

        real(kind(0d0)), dimension(3), intent(IN) :: mono_loc
        type(mono_parameters), intent(IN) :: mymono
        real(kind(0d0)), intent(IN) :: mono_leng
        integer, intent(in) :: j, k, l

        integer :: q
        real(kind(0d0)) :: h, hx, hy, hz
        real(kind(0d0)) :: hxnew, hynew
        real(kind(0d0)) :: sig
        real(kind(0d0)) :: f_delta

        if (n == 0) then
            sig = dx(j)
            sig = sig*2.5d0
        else if (p == 0) then
            sig = maxval((/dx(j), dy(k)/))
            sig = sig*2.5d0
        else
            sig = maxval((/dx(j), dy(k), dz(l)/))
            sig = sig*2.5d0
        end if

        if (n == 0) then      !1D
            if (mymono%support == 1) then
                ! 1D delta function
                hx = abs(mono_loc(1) - x_cc(j))

                f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                          dexp(-0.5d0*(hx/(sig/2.d0))**2.d0)
            else if (mymono%support == 0) then
                ! Support for all x
                f_delta = 1.d0
            end if
        else if (p == 0) then !2D
            hx = mono_loc(1) - x_cc(j)
            hy = mono_loc(2) - y_cc(k)
            if (mymono%support == 1) then
                ! 2D delta function
                sig = mono_leng/20.d0
                h = dsqrt(hx**2.d0 + hy**2.d0)

                f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                          dexp(-0.5d0*((h/(sig/2.d0))**2.d0))
            else if (mymono%support == 2) then
                !only support for y \pm some value
                if (abs(hy) < mymono%length) then
                    f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                              dexp(-0.5d0*(hx/(sig/2.d0))**2.d0)
                else
                    f_delta = 0d0
                end if
            else if (mymono%support == 3) then
                ! Only support along some line
                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)

                ! Rotate actual point by -theta
                hxnew = cos(mymono%dir)*hx + sin(mymono%dir)*hy
                hynew = -1.d0*sin(mymono%dir)*hx + cos(mymono%dir)*hy
                if (abs(hynew) < mymono%loc(3)/2.d0) then
                    f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2.d0))**2.d0)
                else
                    f_delta = 0d0
                end if
            else if (mymono%support == 4) then
                ! Support for all y
                f_delta = 1.d0/(dsqrt(2.d0*pi)*sig)* &
                          dexp(-0.5d0*(hx/sig)**2.d0)
            end if
        else !3D
            if (mymono%support == 3) then
                ! Only support along some patch

                hx = x_cc(j) - mono_loc(1)
                hy = y_cc(k) - mono_loc(2)
                hz = z_cc(l) - mono_loc(3)

                ! Rotate actual point by -theta
                hxnew = cos(mymono%dir)*hx + sin(mymono%dir)*hy
                hynew = -1.d0*sin(mymono%dir)*hx + cos(mymono%dir)*hy

                if (abs(hynew) < mymono%length/2. .and. &
                    abs(hz) < mymono%length/2.) then
                    f_delta = 1.d0/(dsqrt(2.d0*pi)*sig/2.d0)* &
                              dexp(-0.5d0*(hxnew/(sig/2.d0))**2.d0)
                else
                    f_delta = 0d0
                end if
            else
                print '(a)', 'Monopole support not properly defined'
                call s_mpi_abort()
            end if
        end if

    end function f_delta

    !>  The purpose of this procedure is to compute the interface
        !!      sharpening regularization source terms. Only applicable
        !!      for 2-fluid system!
        !!  @param i Dimensional split index
        !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_compute_regularization_source(i, q_prim_vf) ! -----------------

        integer, intent(IN) :: i
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), allocatable :: var
        type(scalar_field), allocatable :: grad_x, grad_y, grad_z
        type(scalar_field), allocatable :: alpharho_grad_x, alpharho_grad_y, alpharho_grad_z
        type(scalar_field), allocatable :: norm
        type(scalar_field), allocatable :: un_alpha_x, un_alpha_y, un_alpha_z

        real(kind(0d0)), dimension(0:m, 0:n, 0:p) :: Lheaviside, U0, velmag
        real(kind(0d0)) :: U0_loc, U0_glb
        real(kind(0d0)), dimension(0:m, 0:n, 0:p) :: Rnohat, R1hat, R2hat
        real(kind(0d0)), dimension(num_dims) :: vel

        type(bounds_info) :: ix, iy, iz

        integer :: j, k, l, r !< Generic loop iterators

        ix%beg = -buff_size; iy%beg = -buff_size
        ix%end = m + buff_size; iy%end = n + buff_size
        if (p > 0) then
            iz%beg = -buff_size; iz%end = p + buff_size
        else
            iz%beg = 0; iz%end = 0
        end if
        allocate (var%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (grad_x%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (grad_y%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (grad_z%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (alpharho_grad_x%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (alpharho_grad_y%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (alpharho_grad_z%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (norm%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (un_alpha_x%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (un_alpha_y%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (un_alpha_z%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))

        do j = 0, m
            do k = 0, n
                do l = 0, p
                    if ((q_prim_vf(adv_idx%beg)%sf(j, k, l) > 1d-6) &
                        .and. &
                        (q_prim_vf(adv_idx%beg)%sf(j, k, l) < (1d0 - 1d-6))) then
                        Lheaviside(j, k, l) = 1d0
                    else
                        Lheaviside(j, k, l) = 0d0
                    end if

                    do r = 1, num_dims
                        vel(r) = q_prim_vf(cont_idx%end + r)%sf(j, k, l)
                    end do

                    velmag(j, k, l) = sqrt(dot_product(vel, vel))

                    U0(j, k, l) = 4d0*q_prim_vf(adv_idx%beg)%sf(j, k, l)* &
                                  (1d0 - q_prim_vf(adv_idx%beg)%sf(j, k, l))* &
                                  velmag(j, k, l)
                end do
            end do
        end do

        U0_loc = maxval(U0)
        if (num_procs > 1) then
            call s_mpi_allreduce_max(U0_loc, U0_glb)
        else
            U0_glb = U0_loc
        end if

        var%sf(:, :, :) = q_prim_vf(adv_idx%beg)%sf(:, :, :)
        call s_compute_fd_gradient(var, grad_x, grad_y, grad_z, norm)
        un_alpha_x%sf(:, :, :) = grad_x%sf(:, :, :)/max(norm%sf(:, :, :), sgm_eps)
        un_alpha_y%sf(:, :, :) = grad_y%sf(:, :, :)/max(norm%sf(:, :, :), sgm_eps)
        un_alpha_z%sf(:, :, :) = grad_z%sf(:, :, :)/max(norm%sf(:, :, :), sgm_eps)

        do j = ix%beg, ix%end
            do k = iy%beg, iy%end
                do l = iz%beg, iz%end
                    var%sf(j, k, l) = reg_eps*norm%sf(j, k, l) - q_prim_vf(adv_idx%beg)%sf(j, k, l)* &
                                      (1d0 - q_prim_vf(adv_idx%beg)%sf(j, k, l))
                end do
            end do
        end do
        call s_compute_fd_gradient(var, grad_x, grad_y, grad_z, norm)
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    if (p > 0) then
                        Rnohat(j, k, l) = Lheaviside(j, k, l)*U0_glb* &
                                          (un_alpha_x%sf(j, k, l)*grad_x%sf(j, k, l) + &
                                           un_alpha_y%sf(j, k, l)*grad_y%sf(j, k, l) + &
                                           un_alpha_z%sf(j, k, l)*grad_z%sf(j, k, l))
                    else
                        Rnohat(j, k, l) = Lheaviside(j, k, l)*U0_glb* &
                                          (un_alpha_x%sf(j, k, l)*grad_x%sf(j, k, l) + &
                                           un_alpha_y%sf(j, k, l)*grad_y%sf(j, k, l))
                    end if
                end do
            end do
        end do

        do r = cont_idx%beg, cont_idx%end
            var%sf(:, :, :) = q_prim_vf(r)%sf(:, :, :)
            call s_compute_fd_gradient(var, alpharho_grad_x, alpharho_grad_y, alpharho_grad_z, norm)
            do j = ix%beg, ix%end
                do k = iy%beg, iy%end
                    do l = iz%beg, iz%end
                        if (p > 0) then
                            var%sf(j, k, l) = reg_eps* &
                                              (un_alpha_x%sf(j, k, l)*alpharho_grad_x%sf(j, k, l) + &
                                               un_alpha_y%sf(j, k, l)*alpharho_grad_y%sf(j, k, l) + &
                                               un_alpha_z%sf(j, k, l)*alpharho_grad_z%sf(j, k, l))
                        else
                            var%sf(j, k, l) = reg_eps* &
                                              (un_alpha_x%sf(j, k, l)*alpharho_grad_x%sf(j, k, l) + &
                                               un_alpha_y%sf(j, k, l)*alpharho_grad_y%sf(j, k, l))
                        end if
                    end do
                end do
            end do
            call s_compute_fd_gradient(var, grad_x, grad_y, grad_z, norm)
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        if (p > 0) then
                            var%sf(j, k, l) = Lheaviside(j, k, l)*U0_glb* &
                                              (un_alpha_x%sf(j, k, l)*(grad_x%sf(j, k, l) - &
                                                                       (1d0 - 2d0*q_prim_vf(adv_idx%beg)%sf(j, k, l))*alpharho_grad_x%sf(j, k, l)) + &
                                               un_alpha_y%sf(j, k, l)*(grad_y%sf(j, k, l) - &
                                                                       (1d0 - 2d0*q_prim_vf(adv_idx%beg)%sf(j, k, l))*alpharho_grad_y%sf(j, k, l)) + &
                                               un_alpha_z%sf(j, k, l)*(grad_z%sf(j, k, l) - &
                                                                       (1d0 - 2d0*q_prim_vf(adv_idx%beg)%sf(j, k, l))*alpharho_grad_z%sf(j, k, l)))
                        else
                            var%sf(j, k, l) = Lheaviside(j, k, l)*U0_glb* &
                                              (un_alpha_x%sf(j, k, l)*(grad_x%sf(j, k, l) - &
                                                                       (1d0 - 2d0*q_prim_vf(adv_idx%beg)%sf(j, k, l))*alpharho_grad_x%sf(j, k, l)) + &
                                               un_alpha_y%sf(j, k, l)*(grad_y%sf(j, k, l) - &
                                                                       (1d0 - 2d0*q_prim_vf(adv_idx%beg)%sf(j, k, l))*alpharho_grad_y%sf(j, k, l)))
                        end if
                    end do
                end do
            end do
            if (r == cont_idx%beg) then
                R1hat(:, :, :) = var%sf(0:m, 0:n, 0:p)
            elseif (r == cont_idx%end) then
                R2hat(:, :, :) = var%sf(0:m, 0:n, 0:p)
            end if
        end do

        reg_src_vf(cont_idx%beg)%sf(:, :, :) = R1hat(:, :, :)
        reg_src_vf(cont_idx%end)%sf(:, :, :) = R2hat(:, :, :)
        do r = mom_idx%beg, mom_idx%end
            reg_src_vf(r)%sf(:, :, :) = q_prim_vf(r)%sf(:, :, :)*(R1hat(:, :, :) + R2hat(:, :, :))
        end do
        reg_src_vf(E_idx)%sf(:, :, :) = 5d-1*velmag(:, :, :)**2d0*(R1hat(:, :, :) + R2hat(:, :, :)) + &
                                        (q_prim_vf(E_idx)%sf(:, :, :)*(fluid_pp(1)%gamma - fluid_pp(2)%gamma) + &
                                         fluid_pp(1)%pi_inf - fluid_pp(2)%pi_inf)*Rnohat(:, :, :)
        reg_src_vf(adv_idx%beg)%sf(:, :, :) = Rnohat(:, :, :)
        if (adv_alphan) then
            reg_src_vf(adv_idx%end)%sf(:, :, :) = -Rnohat(:, :, :)
        end if

        deallocate (var%sf, grad_x%sf, grad_y%sf, grad_z%sf, norm%sf)
        deallocate (un_alpha_x%sf, un_alpha_y%sf, un_alpha_z%sf)
        deallocate (alpharho_grad_x%sf, alpharho_grad_y%sf, alpharho_grad_z%sf)

    end subroutine s_compute_regularization_source ! ----------------------------------

    !>  Computes the scalar gradient fields via finite differences
        !!  @param var Variable to compute derivative of
        !!  @param grad_x First coordinate direction component of the derivative
        !!  @param grad_y Second coordinate direction component of the derivative
        !!  @param grad_z Third coordinate direction component of the derivative
        !!  @param norm Norm of the gradient vector
    subroutine s_compute_fd_gradient(var, grad_x, grad_y, grad_z, norm)

        type(scalar_field), intent(IN) :: var
        type(scalar_field), intent(INOUT) :: grad_x
        type(scalar_field), intent(INOUT) :: grad_y
        type(scalar_field), intent(INOUT) :: grad_z
        type(scalar_field), intent(INOUT) :: norm

        type(bounds_info) :: ix, iy, iz

        integer :: j, k, l !< Generic loop iterators

        ix%beg = -buff_size; ix%end = m + buff_size; 
        if (n > 0) then
            iy%beg = -buff_size; iy%end = n + buff_size
            if (p > 0) then
                iz%beg = -buff_size; iz%end = p + buff_size
            else
                iz%beg = -1; iz%end = 1
            end if
        else
            iy%beg = -1; iy%end = 1
        end if

        do j = ix%beg + 1, ix%end - 1
            do k = iy%beg + 1, iy%end - 1
                do l = iz%beg + 1, iz%end - 1
                    grad_x%sf(j, k, l) = (var%sf(j + 1, k, l) - var%sf(j - 1, k, l))/(x_cc(j + 1) - x_cc(j - 1))
                    if (n > 0) then
                        grad_y%sf(j, k, l) = (var%sf(j, k + 1, l) - var%sf(j, k - 1, l))/(y_cc(k + 1) - y_cc(k - 1))
                        if (p > 0) then
                            grad_z%sf(j, k, l) = (var%sf(j, k, l + 1) - var%sf(j, k, l - 1))/(z_cc(l + 1) - z_cc(l - 1))
                        end if
                    end if
                end do
            end do
        end do
        grad_x%sf(ix%beg, :, :) = (-3d0*var%sf(ix%beg, :, :) + 4d0*var%sf(ix%beg + 1, :, :) - var%sf(ix%beg + 2, :, :))/ &
                                  (x_cc(ix%beg + 2) - x_cc(ix%beg))
        grad_x%sf(ix%end, :, :) = (3d0*var%sf(ix%end, :, :) - 4d0*var%sf(ix%end - 1, :, :) + var%sf(ix%end - 2, :, :))/ &
                                  (x_cc(ix%end) - x_cc(ix%end - 2))
        if (n > 0) then
            grad_y%sf(:, iy%beg, :) = (-3d0*var%sf(:, iy%beg, :) + 4d0*var%sf(:, iy%beg + 1, :) - var%sf(:, iy%beg + 2, :))/ &
                                      (y_cc(iy%beg + 2) - y_cc(iy%beg))
            grad_y%sf(:, iy%end, :) = (3d0*var%sf(:, iy%end, :) - 4d0*var%sf(:, iy%end - 1, :) + var%sf(:, iy%end - 2, :))/ &
                                      (y_cc(iy%end) - y_cc(iy%end - 2))
            if (p > 0) then
                grad_z%sf(:, :, iz%beg) = (-3d0*var%sf(:, :, iz%beg) + 4d0*var%sf(:, :, iz%beg + 1) - var%sf(:, :, iz%beg + 2))/ &
                                          (z_cc(iz%beg + 2) - z_cc(iz%beg))
                grad_z%sf(:, :, iz%end) = (3d0*var%sf(:, :, iz%end) - 4d0*var%sf(:, :, iz%end - 1) + var%sf(:, :, iz%end - 2))/ &
                                          (z_cc(iz%end) - z_cc(iz%end - 2))
            end if
        end if

        if (bc_x%beg <= -3) then
            grad_x%sf(0, :, :) = (-3d0*var%sf(0, :, :) + 4d0*var%sf(1, :, :) - var%sf(2, :, :))/ &
                                 (x_cc(2) - x_cc(0))
        end if
        if (bc_x%end <= -3) then
            grad_x%sf(m, :, :) = (3d0*var%sf(m, :, :) - 4d0*var%sf(m - 1, :, :) + var%sf(m - 2, :, :))/ &
                                 (x_cc(m) - x_cc(m - 2))
        end if
        if (n > 0) then
            if (bc_y%beg <= -3 .and. bc_y%beg /= -13) then
                grad_y%sf(:, 0, :) = (-3d0*var%sf(:, 0, :) + 4d0*var%sf(:, 1, :) - var%sf(:, 2, :))/ &
                                     (y_cc(2) - y_cc(0))
            end if
            if (bc_y%end <= -3) then
                grad_y%sf(:, n, :) = (3d0*var%sf(:, n, :) - 4d0*var%sf(:, n - 1, :) + var%sf(:, n - 2, :))/ &
                                     (y_cc(n) - y_cc(n - 2))
            end if
            if (p > 0) then
                if (bc_z%beg <= -3) then
                    grad_z%sf(:, :, 0) = (-3d0*var%sf(:, :, 0) + 4d0*var%sf(:, :, 1) - var%sf(:, :, 2))/ &
                                         (z_cc(2) - z_cc(0))
                end if
                if (bc_z%end <= -3) then
                    grad_z%sf(:, :, p) = (3d0*var%sf(:, :, p) - 4d0*var%sf(:, :, p - 1) + var%sf(:, :, p - 2))/ &
                                         (z_cc(p) - z_cc(p - 2))
                end if
            end if
        end if

        if (p == 0) then
            iz%beg = 0; iz%end = 0
            if (n == 0) then
                iy%beg = 0; iy%end = 0
            end if
        end if
        do j = ix%beg, ix%end
            do k = iy%beg, iy%end
                do l = iz%beg, iz%end
                    if (p > 0) then
                        norm%sf(j, k, l) = sqrt(grad_x%sf(j, k, l)**2d0 + &
                                                grad_y%sf(j, k, l)**2d0 + &
                                                grad_z%sf(j, k, l)**2d0)
                    elseif (n > 0) then
                        norm%sf(j, k, l) = sqrt(grad_x%sf(j, k, l)**2d0 + &
                                                grad_y%sf(j, k, l)**2d0)
                    else
                        norm%sf(j, k, l) = grad_x%sf(j, k, l)
                    end if
                end do
            end do
        end do

    end subroutine s_compute_fd_gradient ! --------------------------------------

    !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_pressure_relaxation_procedure(q_cons_vf) ! ----------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf

        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        real(kind(0d0))                                   ::  pres_relax
        real(kind(0d0)), dimension(num_fluids)            :: pres_K_init
        real(kind(0d0))                                   ::      f_pres
        real(kind(0d0))                                   ::     df_pres
        real(kind(0d0)), dimension(num_fluids)            ::     rho_K_s
        real(kind(0d0))                                   ::   sum_alpha
        real(kind(0d0))                                   ::         rho
        real(kind(0d0))                                   ::    dyn_pres
        real(kind(0d0))                                   ::        E_We
        real(kind(0d0))                                   ::       gamma
        real(kind(0d0))                                   ::      pi_inf
        real(kind(0d0)), dimension(num_fluids)            ::   gamma_min
        real(kind(0d0)), dimension(num_fluids)            ::    pres_inf
        real(kind(0d0)), dimension(2)                     ::          Re
        real(kind(0d0)), dimension(num_fluids, num_fluids) ::          We
        !> @}

        integer :: i, j, k, l, iter !< Generic loop iterators
        integer :: relax !< Relaxation procedure determination variable

        do i = 1, num_fluids
            gamma_min(i) = 1d0/fluid_pp(i)%gamma + 1d0
            pres_inf(i) = fluid_pp(i)%pi_inf/(1d0 + fluid_pp(i)%gamma)
        end do

        do j = 0, m
            do k = 0, n
                do l = 0, p

                    ! Numerical correction of the volume fractions
                    if (mpp_lim) then
                        sum_alpha = 0d0
                        do i = 1, num_fluids
                            if ((q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) .lt. 0d0) .or. &
                                (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) .lt. 0d0)) then
                                q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) = 0d0
                            end if

                            if (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) .gt. 1d0) &
                                q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) = 1d0
                            sum_alpha = sum_alpha + q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                        end do
                        do i = 1, num_fluids
                            q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) = q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)/sum_alpha
                        end do
                    end if

                    ! Pressures relaxation procedure ===================================

                    ! Is the pressure relaxation procedure necessary?
                    relax = 1
                    do i = 1, num_fluids
                        if (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) .gt. (1d0 - sgm_eps)) relax = 0
                    end do

                    if (relax == 1) then
                        ! Initial state
                        pres_relax = 0d0
                        do i = 1, num_fluids
                            if (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) .gt. sgm_eps) then
                                pres_K_init(i) = &
                                    (q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l)/ &
                                     q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) &
                                     - fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma

                                if (pres_K_init(i) .le. -(1d0 - 1d-8)*pres_inf(i) + 1d-8) &
                                    pres_K_init(i) = -(1d0 - 1d-8)*pres_inf(i) + 1d-8
                            else
                                pres_K_init(i) = 0d0
                            end if
                            pres_relax = pres_relax + q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)*pres_K_init(i)
                        end do

                        ! Iterative process for relaxed pressure determination
                        iter = 0
                        f_pres = 1d-9
                        df_pres = 1d9
                        do i = 1, num_fluids
                            rho_K_s(i) = 0d0
                        end do

                        do while (DABS(f_pres) .gt. 1d-10)
                            pres_relax = pres_relax - f_pres/df_pres

                            ! Convergence
                            iter = iter + 1
                            if (iter == 50) then
                                print '(A)', 'Pressure relaxation procedure failed to converge to a solution. Exiting ...'
                                call s_mpi_abort()
                            end if

                            ! Physical pressure
                            do i = 1, num_fluids
                                if (pres_relax .le. -(1d0 - 1d-8)*pres_inf(i) + 1d-8) &
                                    pres_relax = -(1d0 - 1d-8)*pres_inf(i) + 1d0
                            end do

                            ! Newton-Raphson method
                            f_pres = -1d0
                            df_pres = 0d0
                            do i = 1, num_fluids
                                if (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) .gt. sgm_eps) then
                                    rho_K_s(i) = q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)/ &
                                                 max(q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l), sgm_eps) &
                                                 *((pres_relax + pres_inf(i))/(pres_K_init(i) + &
                                                                               pres_inf(i)))**(1d0/gamma_min(i))

                                    f_pres = f_pres + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) &
                                             /rho_K_s(i)

                                    df_pres = df_pres - q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) &
                                              /(gamma_min(i)*rho_K_s(i)*(pres_relax + pres_inf(i)))
                                end if
                            end do

                        end do

                        ! Cell update of the volume fraction
                        do i = 1, num_fluids
                            if (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) .gt. sgm_eps) &
                                q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) = q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) &
                                                                             /rho_K_s(i)
                        end do
                    end if

                    ! ==================================================================

                    ! Mixture-total-energy correction ==================================

                    ! The mixture-total-energy correction of the mixture pressure P is not necessary here
                    ! because the primitive variables are directly recovered later on by the conservative
                    ! variables (see s_convert_conservative_to_primitive_variables called in s_compute_rhs).
                    ! However, the internal-energy equations should be reset with the corresponding mixture
                    ! pressure from the correction. This step is carried out below.

                    call s_convert_to_mixture_variables(q_cons_vf, rho, &
                                                        gamma, pi_inf, &
                                                        Re, We, j, k, l)

                    dyn_pres = 0d0
                    do i = mom_idx%beg, mom_idx%end
                        dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j, k, l)* &
                                   q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
                    end do

                    E_We = 0d0

                    pres_relax = (q_cons_vf(E_idx)%sf(j, k, l) - dyn_pres - pi_inf - E_We)/gamma

                    do i = 1, num_fluids
                        q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) = &
                            q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)* &
                            (fluid_pp(i)%gamma*pres_relax + fluid_pp(i)%pi_inf)
                    end do
                    ! ==================================================================
                end do
            end do
        end do

    end subroutine s_pressure_relaxation_procedure ! -----------------------

    !>  This subroutine compute the TVD flux function
        !!  @param q_cons_vf Cell-averaged conservative variables
        !!  @param q_prim_vf Cell-averaged primitive variables
        !!  @param rhs_vf Cell-averaged RHS variables
        !!  @param i Dimensional splitting index
    subroutine s_get_tvd_flux(q_cons_vf, q_prim_vf, rhs_vf, i) ! -------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf
        integer, intent(in) :: i

        integer :: j, k, l, r !< Generic loop iterators

        real(kind(0d0)) :: top, bottom

        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end
                    ! Upwind direction for slope chosen based on contact velocity
                    ! from Riemann solver (stored in adv_idx of flux_src). Slope
                    ! is calculated using the first volume fraction.
                    if (hi_flux_src_ndqp(i, 0, 0)%vf(adv_idx%beg)%sf(j, k, l) >= 0d0) then
                        top = q_prim_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j, k, l) - &
                              q_prim_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j - 1, k, l)
                        bottom = q_prim_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j + 1, k, l) - &
                                 q_prim_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j, k, l)
                    else
                        top = q_prim_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j + 2, k, l) - &
                              q_prim_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j + 1, k, l)
                        bottom = q_prim_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j + 1, k, l) - &
                                 q_prim_qp(0, 0, 0)%vf(adv_idx%beg)%sf(j, k, l)
                    end if

                    call s_compute_flux_lim(top, bottom, flux_lim_func(j, k, l))
                end do
            end do
        end do

        ! Compile final flux and flux source
        do k = ipsi%beg, ipsi%end, 2
            do j = ichi%beg, ichi%end, 2
                do l = 1, sys_size
                    flux_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) = &
                        lo_flux_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) + &
                        flux_lim_func(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end)* &
                        (hi_flux_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) - &
                         lo_flux_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
                end do

                if (any(Re_size > 0) .or. We_size > 0 .or. hypoelasticity) then
                    do l = mom_idx%beg, E_idx
                        flux_src_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) = &
                            lo_flux_src_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) + &
                            flux_lim_func(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end)* &
                            (hi_flux_src_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) - &
                             lo_flux_src_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
                    end do
                end if
                do l = adv_idx%beg, adv_idx%end
                    flux_src_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) = &
                        lo_flux_src_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) + &
                        flux_lim_func(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end)* &
                        (hi_flux_src_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) - &
                         lo_flux_src_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
                end do
                if (cyl_coord) then
                    do l = 1, sys_size
                        flux_gsrc_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) = &
                            lo_flux_gsrc_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) + &
                            flux_lim_func(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end)* &
                            (hi_flux_gsrc_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end) - &
                             lo_flux_gsrc_ndqp(i, j, k)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
                    end do
                end if
            end do
        end do

        iv%beg = 1; iv%end = adv_idx%end

        call s_average_cell_boundary_values(flux_ndqp(i, :, :))

        if (any(Re_size > 0) .or. We_size > 0 .or. hypoelasticity) then
            iv%beg = mom_idx%beg
        else
            iv%beg = adv_idx%beg
        end if

        if (riemann_solver /= 1) iv%end = adv_idx%beg

        call s_average_cell_boundary_values(flux_src_ndqp(i, :, :))
        call s_average_cell_boundary_values(flux_gsrc_ndqp(i, :, :))
    end subroutine s_get_tvd_flux

    !>  Computes viscous terms
        !!  @param q_cons_vf Cell-averaged conservative variables
        !!  @param q_prim_vf Cell-averaged primitive variables
        !!  @param rhs_vf Cell-averaged RHS variables
    subroutine s_get_viscous(q_cons_vf, q_prim_vf, rhs_vf) ! -------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf

        integer :: i, j, k, l, r !< Generic loop iterators

        do i = 1, num_dims
            ! WENO reconstruct variables to cell boundaries
            if (weno_vars == 1) then

                iv%beg = 1; iv%end = mom_idx%end

                call s_reconstruct_cell_boundary_values( &
                    q_cons_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                    qL_cons_ndqp(i, :, :), &
                    qR_cons_ndqp(i, :, :), &
                    dflt_int, i)

                do k = ipsi%beg, ipsi%end, 2
                    do j = ichi%beg, ichi%end, 2
                        do l = mom_idx%beg, mom_idx%end

                            qL_prim_ndqp(i, j, k)%vf(l)%sf = sgm_eps
                            qR_prim_ndqp(i, j, k)%vf(l)%sf = sgm_eps

                            do r = 1, cont_idx%end
                                qL_prim_ndqp(i, j, k)%vf(l)%sf = &
                                    qL_prim_ndqp(i, j, k)%vf(l)%sf + &
                                    qL_cons_ndqp(i, j, k)%vf(r)%sf
                                qR_prim_ndqp(i, j, k)%vf(l)%sf = &
                                    qR_prim_ndqp(i, j, k)%vf(l)%sf + &
                                    qR_cons_ndqp(i, j, k)%vf(r)%sf
                            end do

                            qL_prim_ndqp(i, j, k)%vf(l)%sf = &
                                qL_cons_ndqp(i, j, k)%vf(l)%sf/ &
                                qL_prim_ndqp(i, j, k)%vf(l)%sf
                            qR_prim_ndqp(i, j, k)%vf(l)%sf = &
                                qR_cons_ndqp(i, j, k)%vf(l)%sf/ &
                                qR_prim_ndqp(i, j, k)%vf(l)%sf

                        end do
                    end do
                end do

            else

                iv%beg = mom_idx%beg; iv%end = mom_idx%end

                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                    qL_prim_ndqp(i, :, :), &
                    qR_prim_ndqp(i, :, :), &
                    dflt_int, i)

            end if

            iv%beg = mom_idx%beg; iv%end = mom_idx%end

            call s_average_cell_boundary_values(qL_prim_ndqp(i, :, :))
            call s_average_cell_boundary_values(qR_prim_ndqp(i, :, :))
        end do

        if (weno_Re_flux) then
            ! Compute velocity gradient at cell centers using scalar
            ! divergence theorem
            do i = 1, num_dims

                if (i == 1) then
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim_ndqp(i, 0, 0)%vf(iv%beg:iv%end), &
                        qR_prim_ndqp(i, 0, 0)%vf(iv%beg:iv%end), &
                        dq_prim_dx_qp(0, 0, 0)%vf(iv%beg:iv%end), i)
                elseif (i == 2) then
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim_ndqp(i, 0, 0)%vf(iv%beg:iv%end), &
                        qR_prim_ndqp(i, 0, 0)%vf(iv%beg:iv%end), &
                        dq_prim_dy_qp(0, 0, 0)%vf(iv%beg:iv%end), i)
                else
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim_ndqp(i, 0, 0)%vf(iv%beg:iv%end), &
                        qR_prim_ndqp(i, 0, 0)%vf(iv%beg:iv%end), &
                        dq_prim_dz_qp(0, 0, 0)%vf(iv%beg:iv%end), i)
                end if

            end do

        else ! Compute velocity gradient at cell centers using finite differences

            iv%beg = mom_idx%beg; iv%end = mom_idx%end

            do k = iv%beg, iv%end

                do j = ix%beg + 1, ix%end
                    dqL_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(j, :, :) = &
                        (q_prim_qp(0, 0, 0)%vf(k)%sf(j, :, :) - &
                         q_prim_qp(0, 0, 0)%vf(k)%sf(j - 1, :, :))/ &
                        (x_cc(j) - x_cc(j - 1))
                end do

                do j = ix%beg, ix%end - 1
                    dqR_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(j, :, :) = &
                        (q_prim_qp(0, 0, 0)%vf(k)%sf(j + 1, :, :) - &
                         q_prim_qp(0, 0, 0)%vf(k)%sf(j, :, :))/ &
                        (x_cc(j + 1) - x_cc(j))
                end do

                if (n > 0) then
                    do j = iy%beg + 1, iy%end
                        dqL_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, j, :) = &
                            (q_prim_qp(0, 0, 0)%vf(k)%sf(:, j, :) - &
                             q_prim_qp(0, 0, 0)%vf(k)%sf(:, j - 1, :))/ &
                            (y_cc(j) - y_cc(j - 1))
                    end do
                    do j = iy%beg, iy%end - 1
                        dqR_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, j, :) = &
                            (q_prim_qp(0, 0, 0)%vf(k)%sf(:, j + 1, :) - &
                             q_prim_qp(0, 0, 0)%vf(k)%sf(:, j, :))/ &
                            (y_cc(j + 1) - y_cc(j))
                    end do
                    do j = iy%beg + 1, iy%end
                        dqL_prim_dx_ndqp(2, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j, :) = &
                            (dqL_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j, :) + &
                             dqR_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j, :) + &
                             dqL_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j - 1, :) + &
                             dqR_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j - 1, :))
                    end do
                    do j = iy%beg, iy%end - 1
                        dqR_prim_dx_ndqp(2, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j, :) = &
                            (dqL_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j + 1, :) + &
                             dqR_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j + 1, :) + &
                             dqL_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j, :) + &
                             dqR_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, j, :))
                    end do
                    do j = ix%beg + 1, ix%end
                        dqL_prim_dy_ndqp(1, 0, 0)%vf(k)%sf(j, iy%beg + 1:iy%end - 1, :) = &
                            (dqL_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(j, iy%beg + 1:iy%end - 1, :) + &
                             dqR_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(j, iy%beg + 1:iy%end - 1, :) + &
                             dqL_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(j - 1, iy%beg + 1:iy%end - 1, :) + &
                             dqR_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(j - 1, iy%beg + 1:iy%end - 1, :))
                    end do
                    do j = ix%beg, ix%end - 1
                        dqR_prim_dy_ndqp(1, 0, 0)%vf(k)%sf(j, iy%beg + 1:iy%end - 1, :) = &
                            (dqL_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(j + 1, iy%beg + 1:iy%end - 1, :) + &
                             dqR_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(j + 1, iy%beg + 1:iy%end - 1, :) + &
                             dqL_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(j, iy%beg + 1:iy%end - 1, :) + &
                             dqR_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(j, iy%beg + 1:iy%end - 1, :))
                    end do
                    dqL_prim_dx_ndqp(2, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, iy%beg + 1:iy%end, :) = 25d-2* &
                                                                                                      dqL_prim_dx_ndqp(2, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, iy%beg + 1:iy%end, :)
                    dqR_prim_dx_ndqp(2, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, iy%beg:iy%end - 1, :) = 25d-2* &
                                                                                                      dqR_prim_dx_ndqp(2, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, iy%beg:iy%end - 1, :)
                    dqL_prim_dy_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end, iy%beg + 1:iy%end - 1, :) = 25d-2* &
                                                                                                      dqL_prim_dy_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end, iy%beg + 1:iy%end - 1, :)
                    dqR_prim_dy_ndqp(1, 0, 0)%vf(k)%sf(ix%beg:ix%end - 1, iy%beg + 1:iy%end - 1, :) = 25d-2* &
                                                                                                      dqR_prim_dy_ndqp(1, 0, 0)%vf(k)%sf(ix%beg:ix%end - 1, iy%beg + 1:iy%end - 1, :)

                    if (p > 0) then

                        do j = iz%beg + 1, iz%end
                            dqL_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, :, j) = &
                                (q_prim_qp(0, 0, 0)%vf(k)%sf(:, :, j) - &
                                 q_prim_qp(0, 0, 0)%vf(k)%sf(:, :, j - 1))/ &
                                (z_cc(j) - z_cc(j - 1))
                        end do
                        do j = iz%beg, iz%end - 1
                            dqR_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, :, j) = &
                                (q_prim_qp(0, 0, 0)%vf(k)%sf(:, :, j + 1) - &
                                 q_prim_qp(0, 0, 0)%vf(k)%sf(:, :, j))/ &
                                (z_cc(j + 1) - z_cc(j))
                        end do
                        do j = ix%beg + 1, ix%end
                            dqL_prim_dz_ndqp(1, 0, 0)%vf(k)%sf(j, :, iz%beg + 1:iz%end - 1) = &
                                (dqL_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(j, :, iz%beg + 1:iz%end - 1) + &
                                 dqR_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(j, :, iz%beg + 1:iz%end - 1) + &
                                 dqL_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(j - 1, :, iz%beg + 1:iz%end - 1) + &
                                 dqR_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(j - 1, :, iz%beg + 1:iz%end - 1))
                        end do
                        do j = ix%beg, ix%end - 1
                            dqR_prim_dz_ndqp(1, 0, 0)%vf(k)%sf(j, :, iz%beg + 1:iz%end - 1) = &
                                (dqL_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(j + 1, :, iz%beg + 1:iz%end - 1) + &
                                 dqR_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(j + 1, :, iz%beg + 1:iz%end - 1) + &
                                 dqL_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(j, :, iz%beg + 1:iz%end - 1) + &
                                 dqR_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(j, :, iz%beg + 1:iz%end - 1))
                        end do
                        do j = iy%beg + 1, iy%end
                            dqL_prim_dz_ndqp(2, 0, 0)%vf(k)%sf(:, j, iz%beg + 1:iz%end - 1) = &
                                (dqL_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, j, iz%beg + 1:iz%end - 1) + &
                                 dqR_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, j, iz%beg + 1:iz%end - 1) + &
                                 dqL_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, j - 1, iz%beg + 1:iz%end - 1) + &
                                 dqR_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, j - 1, iz%beg + 1:iz%end - 1))
                        end do
                        do j = iy%beg, iy%end - 1
                            dqR_prim_dz_ndqp(2, 0, 0)%vf(k)%sf(:, j, iz%beg + 1:iz%end - 1) = &
                                (dqL_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, j + 1, iz%beg + 1:iz%end - 1) + &
                                 dqR_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, j + 1, iz%beg + 1:iz%end - 1) + &
                                 dqL_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, j, iz%beg + 1:iz%end - 1) + &
                                 dqR_prim_dz_ndqp(3, 0, 0)%vf(k)%sf(:, j, iz%beg + 1:iz%end - 1))
                        end do
                        do j = iz%beg + 1, iz%end
                            dqL_prim_dy_ndqp(3, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j) = &
                                (dqL_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j) + &
                                 dqR_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j) + &
                                 dqL_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j - 1) + &
                                 dqR_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j - 1))
                        end do
                        do j = iz%beg, iz%end - 1
                            dqR_prim_dy_ndqp(3, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j) = &
                                (dqL_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j + 1) + &
                                 dqR_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j + 1) + &
                                 dqL_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j) + &
                                 dqR_prim_dy_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, j))
                        end do
                        do j = iz%beg + 1, iz%end
                            dqL_prim_dx_ndqp(3, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j) = &
                                (dqL_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j) + &
                                 dqR_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j) + &
                                 dqL_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j - 1) + &
                                 dqR_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j - 1))
                        end do
                        do j = iz%beg, iz%end - 1
                            dqR_prim_dx_ndqp(3, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j) = &
                                (dqL_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j + 1) + &
                                 dqR_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j + 1) + &
                                 dqL_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j) + &
                                 dqR_prim_dx_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, j))
                        end do

                        dqL_prim_dz_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end, :, iz%beg + 1:iz%end - 1) = 25d-2* &
                                                                                                          dqL_prim_dz_ndqp(1, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end, :, iz%beg + 1:iz%end - 1)
                        dqR_prim_dz_ndqp(1, 0, 0)%vf(k)%sf(ix%beg:ix%end - 1, :, iz%beg + 1:iz%end - 1) = 25d-2* &
                                                                                                          dqR_prim_dz_ndqp(1, 0, 0)%vf(k)%sf(ix%beg:ix%end - 1, :, iz%beg + 1:iz%end - 1)
                        dqL_prim_dz_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end, iz%beg + 1:iz%end - 1) = 25d-2* &
                                                                                                          dqL_prim_dz_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end, iz%beg + 1:iz%end - 1)
                        dqR_prim_dz_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg:iy%end - 1, iz%beg + 1:iz%end - 1) = 25d-2* &
                                                                                                          dqR_prim_dz_ndqp(2, 0, 0)%vf(k)%sf(:, iy%beg:iy%end - 1, iz%beg + 1:iz%end - 1)
                        dqL_prim_dy_ndqp(3, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, iz%beg + 1:iz%end) = 25d-2* &
                                                                                                          dqL_prim_dy_ndqp(3, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, iz%beg + 1:iz%end)
                        dqR_prim_dy_ndqp(3, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, iz%beg:iz%end - 1) = 25d-2* &
                                                                                                          dqR_prim_dy_ndqp(3, 0, 0)%vf(k)%sf(:, iy%beg + 1:iy%end - 1, iz%beg:iz%end - 1)
                        dqL_prim_dx_ndqp(3, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, iz%beg + 1:iz%end) = 25d-2* &
                                                                                                          dqL_prim_dx_ndqp(3, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, iz%beg + 1:iz%end)
                        dqR_prim_dx_ndqp(3, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, iz%beg:iz%end - 1) = 25d-2* &
                                                                                                          dqR_prim_dx_ndqp(3, 0, 0)%vf(k)%sf(ix%beg + 1:ix%end - 1, :, iz%beg:iz%end - 1)

                        call s_compute_fd_gradient(q_prim_qp(0, 0, 0)%vf(k), &
                                                   dq_prim_dx_qp(0, 0, 0)%vf(k), &
                                                   dq_prim_dy_qp(0, 0, 0)%vf(k), &
                                                   dq_prim_dz_qp(0, 0, 0)%vf(k), &
                                                   gm_vel_qp(0, 0, 0)%vf(k))

                    else

                        call s_compute_fd_gradient(q_prim_qp(0, 0, 0)%vf(k), &
                                                   dq_prim_dx_qp(0, 0, 0)%vf(k), &
                                                   dq_prim_dy_qp(0, 0, 0)%vf(k), &
                                                   dq_prim_dy_qp(0, 0, 0)%vf(k), &
                                                   gm_vel_qp(0, 0, 0)%vf(k))

                    end if

                else
                    call s_compute_fd_gradient(q_prim_qp(0, 0, 0)%vf(k), &
                                               dq_prim_dx_qp(0, 0, 0)%vf(k), &
                                               dq_prim_dx_qp(0, 0, 0)%vf(k), &
                                               dq_prim_dx_qp(0, 0, 0)%vf(k), &
                                               gm_vel_qp(0, 0, 0)%vf(k))

                end if

            end do

        end if

    end subroutine s_get_viscous

    !>  Computes curvatures for surface tension terms
        !!  @param q_cons_vf Cell-averaged conservative variables
        !!  @param q_prim_vf Cell-averaged primitive variables
        !!  @param rhs_vf Cell-averaged RHS variables
    subroutine s_get_crv(q_cons_vf, q_prim_vf, rhs_vf) ! -------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf

        integer :: i, j, k, l, r !< Generic loop iterators

        do l = 1, crv_size

            iv%beg = E_idx + crv_idx(l); iv%end = iv%beg

            ! ==================================================================
            do i = 1, num_dims

                ! Reconstruct the volume fraction at the cell boundaries
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                    qL_prim_ndqp(i, :, :), &
                    qR_prim_ndqp(i, :, :), &
                    dflt_int, i)
                ! Average the quadrature points at the cell boundaries to
                ! obtain higher-order estimate of the volume fraction
                call s_average_cell_boundary_values(qL_prim_ndqp(i, :, :))
                call s_average_cell_boundary_values(qR_prim_ndqp(i, :, :))

            end do
            ! ==================================================================

            if (lsq_deriv) then
                call s_reconstruct_cell_interior_values(dq_prim_dx_qp)
                if (n > 0) then
                    call s_reconstruct_cell_interior_values(dq_prim_dy_qp)
                    if (p > 0) then
                        call s_reconstruct_cell_interior_values(dq_prim_dz_qp)
                    end if
                end if
            else
                ! Apply the scalar divergence theorem to compute the volume
                ! fraction gradient at the cell center
                call s_apply_scalar_divergence_theorem( &
                    qL_prim_ndqp(1, 0, 0)%vf(iv%beg:iv%end), &
                    qR_prim_ndqp(1, 0, 0)%vf(iv%beg:iv%end), &
                    dq_prim_dx_qp(0, 0, 0)%vf(iv%beg:iv%end), 1)

                call s_reconstruct_cell_interior_values(dq_prim_dx_qp)

                if (n > 0) then
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim_ndqp(2, 0, 0)%vf(iv%beg:iv%end), &
                        qR_prim_ndqp(2, 0, 0)%vf(iv%beg:iv%end), &
                        dq_prim_dy_qp(0, 0, 0)%vf(iv%beg:iv%end), 2)

                    call s_reconstruct_cell_interior_values(dq_prim_dy_qp)

                    if (p > 0) then
                        call s_apply_scalar_divergence_theorem( &
                            qL_prim_ndqp(3, 0, 0)%vf(iv%beg:iv%end), &
                            qR_prim_ndqp(3, 0, 0)%vf(iv%beg:iv%end), &
                            dq_prim_dz_qp(0, 0, 0)%vf(iv%beg:iv%end), 3)

                        call s_reconstruct_cell_interior_values(dq_prim_dz_qp)
                    end if
                end if
            end if

            ! ==================================================================

            ! Compute the gradient magnitude at the cell center
            do k = itau%beg, itau%end, 2
                do j = iksi%beg, iksi%end, 2
                    do i = ieta%beg, ieta%end, 2
                        call s_compute_gradient_magnitude( &
                            dq_prim_dx_qp(i, j, k)%vf(iv%beg:iv%end), &
                            dq_prim_dy_qp(i, j, k)%vf(iv%beg:iv%end), &
                            dq_prim_dz_qp(i, j, k)%vf(iv%beg:iv%end), &
                            gm_alpha_qp(i, j, k)%vf(crv_idx(l): &
                                                    crv_idx(l)))
                    end do
                end do
            end do

            ! ==================================================================

            ! Reconstruct the volume fraction gradient at the cell
            ! boundaries
            do i = 1, num_dims

                call s_reconstruct_cell_boundary_values( &
                    dq_prim_dx_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                    dqL_prim_dx_ndqp(i, :, :), &
                    dqR_prim_dx_ndqp(i, :, :), &
                    dflt_int, i)
                if (n > 0) then
                    call s_reconstruct_cell_boundary_values( &
                        dq_prim_dy_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                        dqL_prim_dy_ndqp(i, :, :), &
                        dqR_prim_dy_ndqp(i, :, :), &
                        dflt_int, i)
                    if (p > 0) then
                        call s_reconstruct_cell_boundary_values( &
                            dq_prim_dz_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                            dqL_prim_dz_ndqp(i, :, :), &
                            dqR_prim_dz_ndqp(i, :, :), &
                            dflt_int, i)
                    end if
                end if

            end do
            ! ==================================================================

            ! Compute the gradient magnitude at the cell boundaries
            do k = ipsi%beg, ipsi%end, 2
                do j = ichi%beg, ichi%end, 2
                    do i = 1, num_dims
                        call s_compute_gradient_magnitude( &
                            dqL_prim_dx_ndqp(i, j, k)%vf(iv%beg:iv%end), &
                            dqL_prim_dy_ndqp(i, j, k)%vf(iv%beg:iv%end), &
                            dqL_prim_dz_ndqp(i, j, k)%vf(iv%beg:iv%end), &
                            gm_alphaL_ndqp(i, j, k)%vf(crv_idx(l): &
                                                       crv_idx(l)))
                        call s_compute_gradient_magnitude( &
                            dqR_prim_dx_ndqp(i, j, k)%vf(iv%beg:iv%end), &
                            dqR_prim_dy_ndqp(i, j, k)%vf(iv%beg:iv%end), &
                            dqR_prim_dz_ndqp(i, j, k)%vf(iv%beg:iv%end), &
                            gm_alphaR_ndqp(i, j, k)%vf(crv_idx(l): &
                                                       crv_idx(l)))
                    end do
                end do
            end do
            ! ==================================================================
            if (alt_crv) then
                ! ==================================================================

                ! Compute the normalized volume fraction gradient at the
                ! cell center
                do k = itau%beg, itau%end, 2
                    do j = iksi%beg, iksi%end, 2
                        do i = ieta%beg, ieta%end, 2
                            un_alpha_x_qp(i, j, k)%vf(crv_idx(l))%sf = &
                                dq_prim_dx_qp(i, j, k)%vf(E_idx + crv_idx(l))%sf/ &
                                max(gm_alpha_qp(i, j, k)%vf(crv_idx(l))%sf, sgm_eps)
                            if (n > 0) then
                                un_alpha_y_qp(i, j, k)%vf(crv_idx(l))%sf = &
                                    dq_prim_dy_qp(i, j, k)%vf(E_idx + crv_idx(l))%sf/ &
                                    max(gm_alpha_qp(i, j, k)%vf(crv_idx(l))%sf, sgm_eps)
                                if (p > 0) then
                                    un_alpha_z_qp(i, j, k)%vf(crv_idx(l))%sf = &
                                        dq_prim_dz_qp(i, j, k)%vf(E_idx + crv_idx(l))%sf/ &
                                        max(gm_alpha_qp(i, j, k)%vf(crv_idx(l))%sf, sgm_eps)
                                end if
                            end if
                        end do
                    end do
                end do
                ! ==================================================================

                ! Apply the scalar divergence theorem to compute the gradient
                ! of the gradient magnitude at the cell center
                call s_apply_scalar_divergence_theorem( &
                    gm_alphaL_ndqp(1, 0, 0)%vf(crv_idx(l):crv_idx(l)), &
                    gm_alphaR_ndqp(1, 0, 0)%vf(crv_idx(l):crv_idx(l)), &
                    dgm_alpha_dx_qp(0, 0, 0)%vf(crv_idx(l):crv_idx(l)), 1)
                ! CALL s_reconstruct_cell_interior_values(dgm_alpha_dx_qp)

                if (n > 0) then
                    call s_apply_scalar_divergence_theorem( &
                        gm_alphaL_ndqp(2, 0, 0)%vf(crv_idx(l):crv_idx(l)), &
                        gm_alphaR_ndqp(2, 0, 0)%vf(crv_idx(l):crv_idx(l)), &
                        dgm_alpha_dy_qp(0, 0, 0)%vf(crv_idx(l):crv_idx(l)), 2)
                    ! CALL s_reconstruct_cell_interior_values(dgm_alpha_dy_qp)
                    if (p > 0) then
                        call s_apply_scalar_divergence_theorem( &
                            gm_alphaL_ndqp(3, 0, 0)%vf(crv_idx(l):crv_idx(l)), &
                            gm_alphaR_ndqp(3, 0, 0)%vf(crv_idx(l):crv_idx(l)), &
                            dgm_alpha_dz_qp(0, 0, 0)%vf(crv_idx(l):crv_idx(l)), 3)
                        ! CALL s_reconstruct_cell_interior_values(dgm_alpha_dz_qp)
                    end if
                end if
                ! ==================================================================

                !Apply the vector divergence theorem to compute the
                ! Laplacian of the volume fraction at the cell center
                call s_apply_vector_divergence_theorem( &
                    dqL_prim_dx_ndqp, &
                    dqL_prim_dy_ndqp, &
                    dqL_prim_dz_ndqp, &
                    dqR_prim_dx_ndqp, &
                    dqR_prim_dy_ndqp, &
                    dqR_prim_dz_ndqp, &
                    laplacian_vf(crv_idx(l):crv_idx(l)))
                ! ==================================================================

                ! Compute alternate volume fraction curvature
                alt_kappa_vf(crv_idx(l))%sf = (un_alpha_x_qp(0, 0, 0)%vf(crv_idx(l))%sf* &
                                               dgm_alpha_dx_qp(0, 0, 0)%vf(crv_idx(l))%sf + &
                                               un_alpha_y_qp(0, 0, 0)%vf(crv_idx(l))%sf* &
                                               dgm_alpha_dy_qp(0, 0, 0)%vf(crv_idx(l))%sf + &
                                               un_alpha_z_qp(0, 0, 0)%vf(crv_idx(l))%sf* &
                                               dgm_alpha_dz_qp(0, 0, 0)%vf(crv_idx(l))%sf - &
                                               laplacian_vf(crv_idx(l))%sf)/ &
                                              max(gm_alpha_qp(0, 0, 0)%vf(crv_idx(l))%sf, sgm_eps)
                ! ==================================================================
            else
                ! ==================================================================

                ! Compute the normalized volume fraction gradient at the cell
                ! boundaries
                do k = ipsi%beg, ipsi%end, 2
                    do j = ichi%beg, ichi%end, 2
                        un_alphaL_x_ndqp(1, j, k)%vf(crv_idx(l))%sf = &
                            dqL_prim_dx_ndqp(1, j, k)%vf(E_idx + crv_idx(l))%sf/ &
                            max(gm_alphaL_ndqp(1, j, k)%vf(crv_idx(l))%sf, sgm_eps)
                        un_alphaR_x_ndqp(1, j, k)%vf(crv_idx(l))%sf = &
                            dqR_prim_dx_ndqp(1, j, k)%vf(E_idx + crv_idx(l))%sf/ &
                            max(gm_alphaR_ndqp(1, j, k)%vf(crv_idx(l))%sf, sgm_eps)
                        if (n > 0) then
                            un_alphaL_y_ndqp(2, j, k)%vf(crv_idx(l))%sf = &
                                dqL_prim_dy_ndqp(2, j, k)%vf(E_idx + crv_idx(l))%sf/ &
                                max(gm_alphaL_ndqp(2, j, k)%vf(crv_idx(l))%sf, sgm_eps)
                            un_alphaR_y_ndqp(2, j, k)%vf(crv_idx(l))%sf = &
                                dqR_prim_dy_ndqp(2, j, k)%vf(E_idx + crv_idx(l))%sf/ &
                                max(gm_alphaR_ndqp(2, j, k)%vf(crv_idx(l))%sf, sgm_eps)
                            if (p > 0) then
                                un_alphaL_z_ndqp(3, j, k)%vf(crv_idx(l))%sf = &
                                    dqL_prim_dz_ndqp(3, j, k)%vf(E_idx + crv_idx(l))%sf/ &
                                    max(gm_alphaL_ndqp(3, j, k)%vf(crv_idx(l))%sf, sgm_eps)
                                un_alphaR_z_ndqp(3, j, k)%vf(crv_idx(l))%sf = &
                                    dqR_prim_dz_ndqp(3, j, k)%vf(E_idx + crv_idx(l))%sf/ &
                                    max(gm_alphaR_ndqp(3, j, k)%vf(crv_idx(l))%sf, sgm_eps)
                            end if
                        end if
                    end do
                end do
                ! ==================================================================

                iv%beg = crv_idx(l); iv%end = iv%beg
                ! Average the quadrature points at the cell boundaries to
                ! obtain higher-order estimate of the normalized volume fraction
                ! gradient
                call s_average_cell_boundary_values(un_alphaL_x_ndqp(1, :, :))
                call s_average_cell_boundary_values(un_alphaR_x_ndqp(1, :, :))
                if (n > 0) then
                    call s_average_cell_boundary_values(un_alphaL_y_ndqp(2, :, :))
                    call s_average_cell_boundary_values(un_alphaR_y_ndqp(2, :, :))
                    if (p > 0) then
                        call s_average_cell_boundary_values(un_alphaL_z_ndqp(3, :, :))
                        call s_average_cell_boundary_values(un_alphaR_z_ndqp(3, :, :))
                    end if
                end if
                ! ==================================================================
                if (lsq_deriv .neqv. .true.) then
                    ! Apply the vector divergence theorem to compute the divergence
                    ! of the normalized volume fraction gradient
                    call s_apply_vector_divergence_theorem( &
                        un_alphaL_x_ndqp, &
                        un_alphaL_y_ndqp, &
                        un_alphaL_z_ndqp, &
                        un_alphaR_x_ndqp, &
                        un_alphaR_y_ndqp, &
                        un_alphaR_z_ndqp, &
                        kappa_vf(iv%beg:iv%end))

                    kappa_vf(crv_idx(l))%sf = -kappa_vf(crv_idx(l))%sf
                end if
                ! ==================================================================
            end if
        end do
        ! END: Computing Volume Fraction Gradients and Curvatures ==========

    end subroutine s_get_crv

    !>  The purpose of this procedure is to populate the buffers
        !!      of the conservative variables, depending on the selected
        !!      boundary conditions.
        !!  @param v_vf Scalar field for which buffers are populated
    subroutine s_populate_variables_buffers(v_vf) ! ---------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: v_vf

        integer :: i, j, k !< Generic loop iterators

        ! Population of Buffers in x-direction =============================
        if (bc_x%beg <= -3) then         ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(-j, 0:n, 0:p) = &
                        v_vf(i)%sf(0, 0:n, 0:p)
                end do
            end do

        elseif (bc_x%beg == -2) then     ! Symmetry BC at beginning

            do j = 1, buff_size

                do i = 1, cont_idx%end
                    v_vf(i)%sf(-j, 0:n, 0:p) = &
                        v_vf(i)%sf(j - 1, 0:n, 0:p)
                end do

                v_vf(mom_idx%beg)%sf(-j, 0:n, 0:p) = &
                    -v_vf(mom_idx%beg)%sf(j - 1, 0:n, 0:p)

                do i = mom_idx%beg + 1, sys_size
                    v_vf(i)%sf(-j, 0:n, 0:p) = &
                        v_vf(i)%sf(j - 1, 0:n, 0:p)
                end do

            end do

        elseif (bc_x%beg == -1) then     ! Periodic BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(-j, 0:n, 0:p) = &
                        v_vf(i)%sf(m - (j - 1), 0:n, 0:p)
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                v_vf, 1, -1)

        end if

        if (bc_x%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(m + j, 0:n, 0:p) = &
                        v_vf(i)%sf(m, 0:n, 0:p)
                end do
            end do

        elseif (bc_x%end == -2) then     ! Symmetry BC at end

            do j = 1, buff_size

                do i = 1, cont_idx%end
                    v_vf(i)%sf(m + j, 0:n, 0:p) = &
                        v_vf(i)%sf(m - (j - 1), 0:n, 0:p)
                end do

                v_vf(mom_idx%beg)%sf(m + j, 0:n, 0:p) = &
                    -v_vf(mom_idx%beg)%sf(m - (j - 1), 0:n, 0:p)

                do i = mom_idx%beg + 1, sys_size
                    v_vf(i)%sf(m + j, 0:n, 0:p) = &
                        v_vf(i)%sf(m - (j - 1), 0:n, 0:p)
                end do

            end do

        elseif (bc_x%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(m + j, 0:n, 0:p) = &
                        v_vf(i)%sf(j - 1, 0:n, 0:p)
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                v_vf, 1, 1)

        end if

        ! END: Population of Buffers in x-direction ========================

        ! Population of Buffers in y-direction =============================

        if (n == 0) then

            return

        elseif (bc_y%beg <= -3 .and. bc_y%beg /= -13) then     ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(:, -j, 0:p) = &
                        v_vf(i)%sf(:, 0, 0:p)
                end do
            end do

        elseif (bc_y%beg == -13) then    ! Axis BC at beginning

            do j = 1, buff_size
                do k = 0, p
                    if (z_cc(k) < pi) then
                        do i = 1, mom_idx%beg
                            v_vf(i)%sf(:, -j, k) = &
                                v_vf(i)%sf(:, j - 1, k + ((p + 1)/2))
                        end do

                        v_vf(mom_idx%beg + 1)%sf(:, -j, k) = &
                            -v_vf(mom_idx%beg + 1)%sf(:, j - 1, k + ((p + 1)/2))

                        v_vf(mom_idx%end)%sf(:, -j, k) = &
                            -v_vf(mom_idx%end)%sf(:, j - 1, k + ((p + 1)/2))

                        do i = E_idx, sys_size
                            v_vf(i)%sf(:, -j, k) = &
                                v_vf(i)%sf(:, j - 1, k + ((p + 1)/2))
                        end do
                    else
                        do i = 1, mom_idx%beg
                            v_vf(i)%sf(:, -j, k) = &
                                v_vf(i)%sf(:, j - 1, k - ((p + 1)/2))
                        end do

                        v_vf(mom_idx%beg + 1)%sf(:, -j, k) = &
                            -v_vf(mom_idx%beg + 1)%sf(:, j - 1, k - ((p + 1)/2))

                        v_vf(mom_idx%end)%sf(:, -j, k) = &
                            -v_vf(mom_idx%end)%sf(:, j - 1, k - ((p + 1)/2))

                        do i = E_idx, sys_size
                            v_vf(i)%sf(:, -j, k) = &
                                v_vf(i)%sf(:, j - 1, k - ((p + 1)/2))
                        end do
                    end if
                end do
            end do

        elseif (bc_y%beg == -2) then     ! Symmetry BC at beginning

            do j = 1, buff_size

                do i = 1, mom_idx%beg
                    v_vf(i)%sf(:, -j, 0:p) = &
                        v_vf(i)%sf(:, j - 1, 0:p)
                end do

                v_vf(mom_idx%beg + 1)%sf(:, -j, 0:p) = &
                    -v_vf(mom_idx%beg + 1)%sf(:, j - 1, 0:p)

                do i = mom_idx%beg + 2, sys_size
                    v_vf(i)%sf(:, -j, 0:p) = &
                        v_vf(i)%sf(:, j - 1, 0:p)
                end do

            end do

        elseif (bc_y%beg == -1) then     ! Periodic BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(:, -j, 0:p) = &
                        v_vf(i)%sf(:, n - (j - 1), 0:p)
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                v_vf, 2, -1)

        end if

        if (bc_y%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(:, n + j, 0:p) = &
                        v_vf(i)%sf(:, n, 0:p)
                end do
            end do

        elseif (bc_y%end == -2) then     ! Symmetry BC at end

            do j = 1, buff_size

                do i = 1, mom_idx%beg
                    v_vf(i)%sf(:, n + j, 0:p) = &
                        v_vf(i)%sf(:, n - (j - 1), 0:p)
                end do

                v_vf(mom_idx%beg + 1)%sf(:, n + j, 0:p) = &
                    -v_vf(mom_idx%beg + 1)%sf(:, n - (j - 1), 0:p)

                do i = mom_idx%beg + 2, sys_size
                    v_vf(i)%sf(:, n + j, 0:p) = &
                        v_vf(i)%sf(:, n - (j - 1), 0:p)
                end do

            end do

        elseif (bc_y%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(:, n + j, 0:p) = &
                        v_vf(i)%sf(:, j - 1, 0:p)
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                v_vf, 2, 1)

        end if

        ! END: Population of Buffers in y-direction ========================

        ! Population of Buffers in z-direction =============================

        if (p == 0) then

            return

        elseif (bc_z%beg <= -3) then     ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(:, :, -j) = &
                        v_vf(i)%sf(:, :, 0)
                end do
            end do

        elseif (bc_z%beg == -2) then     ! Symmetry BC at beginning

            do j = 1, buff_size

                do i = 1, mom_idx%beg + 1
                    v_vf(i)%sf(:, :, -j) = &
                        v_vf(i)%sf(:, :, j - 1)
                end do

                v_vf(mom_idx%end)%sf(:, :, -j) = &
                    -v_vf(mom_idx%end)%sf(:, :, j - 1)

                do i = E_idx, sys_size
                    v_vf(i)%sf(:, :, -j) = &
                        v_vf(i)%sf(:, :, j - 1)
                end do

            end do

        elseif (bc_z%beg == -1) then     ! Periodic BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(:, :, -j) = &
                        v_vf(i)%sf(:, :, p - (j - 1))
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                v_vf, 3, -1)

        end if

        if (bc_z%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(:, :, p + j) = &
                        v_vf(i)%sf(:, :, p)
                end do
            end do

        elseif (bc_z%end == -2) then     ! Symmetry BC at end

            do j = 1, buff_size

                do i = 1, mom_idx%beg + 1
                    v_vf(i)%sf(:, :, p + j) = &
                        v_vf(i)%sf(:, :, p - (j - 1))
                end do

                v_vf(mom_idx%end)%sf(:, :, p + j) = &
                    -v_vf(mom_idx%end)%sf(:, :, p - (j - 1))

                do i = E_idx, sys_size
                    v_vf(i)%sf(:, :, p + j) = &
                        v_vf(i)%sf(:, :, p - (j - 1))
                end do

            end do

        elseif (bc_z%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(:, :, p + j) = &
                        v_vf(i)%sf(:, :, j - 1)
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                v_vf, 3, 1)

        end if

        ! END: Population of Buffers in z-direction ========================

    end subroutine s_populate_variables_buffers ! -------------

    !>  The purpose of this procedure is to populate the buffers
        !!      of the conservative variables, depending on the selected
        !!      boundary conditions.
    subroutine s_populate_conservative_variables_buffers() ! ---------------

        integer :: i, j, k, l, r !< Generic loop iterators

        ! Population of Buffers in x-direction =============================

        if (bc_x%beg <= -3) then         ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(-j, 0:n, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(0, 0:n, 0:p)
                end do
            end do

        elseif (bc_x%beg == -2) then     ! Symmetry BC at beginning

            do j = 1, buff_size

                do i = 1, cont_idx%end
                    q_cons_qp(0, 0, 0)%vf(i)%sf(-j, 0:n, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(j - 1, 0:n, 0:p)
                end do

                q_cons_qp(0, 0, 0)%vf(mom_idx%beg)%sf(-j, 0:n, 0:p) = &
                    -q_cons_qp(0, 0, 0)%vf(mom_idx%beg)%sf(j - 1, 0:n, 0:p)

                do i = mom_idx%beg + 1, sys_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(-j, 0:n, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(j - 1, 0:n, 0:p)
                end do

            end do

        elseif (bc_x%beg == -1) then     ! Periodic BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(-j, 0:n, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(m - (j - 1), 0:n, 0:p)
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp(0, 0, 0)%vf, 1, -1)

        end if

        if (bc_x%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(m + j, 0:n, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(m, 0:n, 0:p)
                end do
            end do

        elseif (bc_x%end == -2) then     ! Symmetry BC at end

            do j = 1, buff_size

                do i = 1, cont_idx%end
                    q_cons_qp(0, 0, 0)%vf(i)%sf(m + j, 0:n, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(m - (j - 1), 0:n, 0:p)
                end do

                q_cons_qp(0, 0, 0)%vf(mom_idx%beg)%sf(m + j, 0:n, 0:p) = &
                    -q_cons_qp(0, 0, 0)%vf(mom_idx%beg)%sf(m - (j - 1), 0:n, 0:p)

                do i = mom_idx%beg + 1, sys_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(m + j, 0:n, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(m - (j - 1), 0:n, 0:p)
                end do

            end do

        elseif (bc_x%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(m + j, 0:n, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(j - 1, 0:n, 0:p)
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp(0, 0, 0)%vf, 1, 1)

        end if

        if (commute_err) then
            do i = ieta%beg, ieta%end, 2
                do j = iksi%beg, iksi%end, 2
                    do k = itau%beg, itau%end, 2
                        do l = 1, sys_size
                            do r = 1, buff_size
                                q_cons_qp(i, j, k)%vf(l)%sf(-r, 0:n, 0:p) = &
                                    q_cons_qp(0, 0, 0)%vf(l)%sf(-r, 0:n, 0:p)
                                q_cons_qp(i, j, k)%vf(l)%sf(m + r, 0:n, 0:p) = &
                                    q_cons_qp(0, 0, 0)%vf(l)%sf(m + r, 0:n, 0:p)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        ! END: Population of Buffers in x-direction ========================

        ! Population of Buffers in y-direction =============================

        if (n == 0) then

            return

        elseif (bc_y%beg <= -3 .and. bc_y%beg /= -13) then     ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, -j, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, 0, 0:p)
                end do
            end do

        elseif (bc_y%beg == -13) then    ! Axis BC at beginning

            do j = 1, buff_size
                do k = 0, p
                    if (z_cc(k) < pi) then
                        do i = 1, mom_idx%beg
                            q_cons_qp(0, 0, 0)%vf(i)%sf(:, -j, k) = &
                                q_cons_qp(0, 0, 0)%vf(i)%sf(:, j - 1, k + ((p + 1)/2))
                        end do

                        q_cons_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(:, -j, k) = &
                            -q_cons_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(:, j - 1, k + ((p + 1)/2))

                        q_cons_qp(0, 0, 0)%vf(mom_idx%end)%sf(:, -j, k) = &
                            -q_cons_qp(0, 0, 0)%vf(mom_idx%end)%sf(:, j - 1, k + ((p + 1)/2))

                        do i = E_idx, sys_size
                            q_cons_qp(0, 0, 0)%vf(i)%sf(:, -j, k) = &
                                q_cons_qp(0, 0, 0)%vf(i)%sf(:, j - 1, k + ((p + 1)/2))
                        end do
                    else
                        do i = 1, mom_idx%beg
                            q_cons_qp(0, 0, 0)%vf(i)%sf(:, -j, k) = &
                                q_cons_qp(0, 0, 0)%vf(i)%sf(:, j - 1, k - ((p + 1)/2))
                        end do

                        q_cons_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(:, -j, k) = &
                            -q_cons_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(:, j - 1, k - ((p + 1)/2))

                        q_cons_qp(0, 0, 0)%vf(mom_idx%end)%sf(:, -j, k) = &
                            -q_cons_qp(0, 0, 0)%vf(mom_idx%end)%sf(:, j - 1, k - ((p + 1)/2))

                        do i = E_idx, sys_size
                            q_cons_qp(0, 0, 0)%vf(i)%sf(:, -j, k) = &
                                q_cons_qp(0, 0, 0)%vf(i)%sf(:, j - 1, k - ((p + 1)/2))
                        end do
                    end if
                end do
            end do

        elseif (bc_y%beg == -2) then     ! Symmetry BC at beginning

            do j = 1, buff_size

                do i = 1, mom_idx%beg
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, -j, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, j - 1, 0:p)
                end do

                q_cons_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(:, -j, 0:p) = &
                    -q_cons_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(:, j - 1, 0:p)

                do i = mom_idx%beg + 2, sys_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, -j, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, j - 1, 0:p)
                end do

            end do

        elseif (bc_y%beg == -1) then     ! Periodic BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, -j, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, n - (j - 1), 0:p)
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp(0, 0, 0)%vf, 2, -1)

        end if

        if (bc_y%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, n + j, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, n, 0:p)
                end do
            end do

        elseif (bc_y%end == -2) then     ! Symmetry BC at end

            do j = 1, buff_size

                do i = 1, mom_idx%beg
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, n + j, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, n - (j - 1), 0:p)
                end do

                q_cons_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(:, n + j, 0:p) = &
                    -q_cons_qp(0, 0, 0)%vf(mom_idx%beg + 1)%sf(:, n - (j - 1), 0:p)

                do i = mom_idx%beg + 2, sys_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, n + j, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, n - (j - 1), 0:p)
                end do

            end do

        elseif (bc_y%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, n + j, 0:p) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, j - 1, 0:p)
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp(0, 0, 0)%vf, 2, 1)

        end if

        if (commute_err) then
            do i = ieta%beg, ieta%end, 2
                do j = iksi%beg, iksi%end, 2
                    do k = itau%beg, itau%end, 2
                        do l = 1, sys_size
                            do r = 1, buff_size
                                q_cons_qp(i, j, k)%vf(l)%sf(:, -r, 0:p) = &
                                    q_cons_qp(0, 0, 0)%vf(l)%sf(:, -r, 0:p)
                                q_cons_qp(i, j, k)%vf(l)%sf(:, n + r, 0:p) = &
                                    q_cons_qp(0, 0, 0)%vf(l)%sf(:, n + r, 0:p)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        ! END: Population of Buffers in y-direction ========================

        ! Population of Buffers in z-direction =============================

        if (p == 0) then

            return

        elseif (bc_z%beg <= -3) then     ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, -j) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, 0)
                end do
            end do

        elseif (bc_z%beg == -2) then     ! Symmetry BC at beginning

            do j = 1, buff_size

                do i = 1, mom_idx%beg + 1
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, -j) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, j - 1)
                end do

                q_cons_qp(0, 0, 0)%vf(mom_idx%end)%sf(:, :, -j) = &
                    -q_cons_qp(0, 0, 0)%vf(mom_idx%end)%sf(:, :, j - 1)

                do i = E_idx, sys_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, -j) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, j - 1)
                end do

            end do

        elseif (bc_z%beg == -1) then     ! Periodic BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, -j) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, p - (j - 1))
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp(0, 0, 0)%vf, 3, -1)

        end if

        if (bc_z%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, p + j) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, p)
                end do
            end do

        elseif (bc_z%end == -2) then     ! Symmetry BC at end

            do j = 1, buff_size

                do i = 1, mom_idx%beg + 1
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, p + j) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, p - (j - 1))
                end do

                q_cons_qp(0, 0, 0)%vf(mom_idx%end)%sf(:, :, p + j) = &
                    -q_cons_qp(0, 0, 0)%vf(mom_idx%end)%sf(:, :, p - (j - 1))

                do i = E_idx, sys_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, p + j) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, p - (j - 1))
                end do

            end do

        elseif (bc_z%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, p + j) = &
                        q_cons_qp(0, 0, 0)%vf(i)%sf(:, :, j - 1)
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp(0, 0, 0)%vf, 3, 1)

        end if

        if (commute_err) then
            do i = ieta%beg, ieta%end, 2
                do j = iksi%beg, iksi%end, 2
                    do k = itau%beg, itau%end, 2
                        do l = 1, sys_size
                            do r = 1, buff_size
                                q_cons_qp(i, j, k)%vf(l)%sf(:, :, -r) = &
                                    q_cons_qp(0, 0, 0)%vf(l)%sf(:, :, -r)
                                q_cons_qp(i, j, k)%vf(l)%sf(:, :, p + r) = &
                                    q_cons_qp(0, 0, 0)%vf(l)%sf(:, :, p + r)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        ! END: Population of Buffers in z-direction ========================

    end subroutine s_populate_conservative_variables_buffers ! -------------

    !>  The goal of this procedure is to WENO-reconstruct the
        !!      inputted cell-averaged variables at the cell-interior
        !!      Gaussian quadrature points.
        !!  @param v_qp Inputted cell-averaged variables and their WENO-reconstructed
        !!          values at the chosen cell-interior Gaussian quadrature points
    subroutine s_reconstruct_cell_interior_values(v_qp) ! ------------------

        type(vector_field), &
            dimension(ieta%beg:ieta%end, &
                      iksi%beg:iksi%end, &
                      itau%beg:itau%end), &
            intent(INOUT) :: v_qp

        type(bounds_info) :: is1, is2, is3 !< Indical bounds in the s1-, s2- and s3-directions

        ! Reconstruction in x-direction ====================================
        if (commute_err .neqv. .true.) return

        is1 = ix; is2 = iy; is3 = iz

        is1%beg = is1%beg + weno_polyn; is1%end = is1%end - weno_polyn

        call s_weno(v_qp(0, 0, 0)%vf(iv%beg:iv%end), &
                    v_qp(-1, 0, 0)%vf(iv%beg:iv%end), &
                    v_qp(1, 0, 0)%vf(iv%beg:iv%end), &
                    dflt_int, 1, 1, 2, is1, is2, is3)
        ! ==================================================================

        ! Reconstruction in y-direction ====================================
        if (n == 0) return

        is2%beg = is2%beg + weno_polyn; is2%end = is2%end - weno_polyn

        call s_weno(v_qp(-1, 0, 0)%vf(iv%beg:iv%end), &
                    v_qp(-1, -1, 0)%vf(iv%beg:iv%end), &
                    v_qp(-1, 1, 0)%vf(iv%beg:iv%end), &
                    dflt_int, 1, 2, 2, is1, is2, is3)
        call s_weno(v_qp(1, 0, 0)%vf(iv%beg:iv%end), &
                    v_qp(1, -1, 0)%vf(iv%beg:iv%end), &
                    v_qp(1, 1, 0)%vf(iv%beg:iv%end), &
                    dflt_int, 1, 2, 2, is1, is2, is3)
        ! ==================================================================

        ! Reconstruction in z-direction ====================================
        if (p == 0) return

        is3%beg = is3%beg + weno_polyn; is3%end = is3%end - weno_polyn

        call s_weno(v_qp(-1, -1, 0)%vf(iv%beg:iv%end), &
                    v_qp(-1, -1, -1)%vf(iv%beg:iv%end), &
                    v_qp(-1, -1, 1)%vf(iv%beg:iv%end), &
                    dflt_int, 1, 3, 2, is1, is2, is3)
        call s_weno(v_qp(-1, 1, 0)%vf(iv%beg:iv%end), &
                    v_qp(-1, 1, -1)%vf(iv%beg:iv%end), &
                    v_qp(-1, 1, 1)%vf(iv%beg:iv%end), &
                    dflt_int, 1, 3, 2, is1, is2, is3)
        call s_weno(v_qp(1, -1, 0)%vf(iv%beg:iv%end), &
                    v_qp(1, -1, -1)%vf(iv%beg:iv%end), &
                    v_qp(1, -1, 1)%vf(iv%beg:iv%end), &
                    dflt_int, 1, 3, 2, is1, is2, is3)
        call s_weno(v_qp(1, 1, 0)%vf(iv%beg:iv%end), &
                    v_qp(1, 1, -1)%vf(iv%beg:iv%end), &
                    v_qp(1, 1, 1)%vf(iv%beg:iv%end), &
                    dflt_int, 1, 3, 2, is1, is2, is3)
        ! ==================================================================

    end subroutine s_reconstruct_cell_interior_values ! --------------------

    !>  The goal of this subroutine is to numerically approximate
        !!      the cell-interior integral-average of the given variables
        !!      by taking the arithmetic mean of their WENO-reconstructed
        !!      values at the cell-interior Gaussian quadrature points.
        !!  @param v_qp The inputted WENO-reconstructed values of cell-averaged variables
        !!      at cell-interior Gaussian quadrature points and their numerically
        !!      approximated cell-interior integral-average.
    subroutine s_average_cell_interior_values(v_qp) ! ----------------------

        type(vector_field), &
            dimension(ieta%beg:ieta%end, &
                      iksi%beg:iksi%end, &
                      itau%beg:itau%end), &
            intent(INOUT) :: v_qp

        integer :: i, j, k, l !< Generic loop iterators

        if (commute_err .neqv. .true.) return

        do l = iv%beg, iv%end

            v_qp(0, 0, 0)%vf(l)%sf = 0d0

            do k = itau%beg, itau%end, 2
                do j = iksi%beg, iksi%end, 2
                    do i = ieta%beg, ieta%end, 2
                        v_qp(0, 0, 0)%vf(l)%sf = v_qp(0, 0, 0)%vf(l)%sf &
                                                 + v_qp(i, j, k)%vf(l)%sf
                    end do
                end do
            end do

            v_qp(0, 0, 0)%vf(l)%sf = v_qp(0, 0, 0)%vf(l)%sf/(2**num_dims)

        end do

    end subroutine s_average_cell_interior_values ! ------------------------

    !>  The purpose of this subroutine is to WENO-reconstruct the
        !!      left and the right cell-boundary values, including values
        !!      at the Gaussian quadrature points, from the cell-averaged
        !!      variables.
        !!  @param v_vf Cell-average variables
        !!  @param vL_qp Left WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param vR_qp Right WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param cd_vars Characteristic decomposition state variables type
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_reconstruct_cell_boundary_values(v_vf, vL_qp, vR_qp, & ! -
                                                  cd_vars, norm_dir)

        type(scalar_field), dimension(iv%beg:iv%end), intent(IN) :: v_vf

        type(vector_field), &
            dimension(ichi%beg:ichi%end, &
                      ipsi%beg:ipsi%end), &
            intent(INOUT) :: vL_qp, vR_qp

        integer, intent(IN) :: cd_vars
        integer, intent(IN) :: norm_dir

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        type(bounds_info) :: is1, is2, is3 !< Indical bounds in the s1-, s2- and s3-directions

        ! Reconstruction in s1-direction ===================================
        is1 = ix; is2 = iy; is3 = iz

        if (norm_dir == 1) then
            weno_dir = 1; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn
        elseif (norm_dir == 2) then
            weno_dir = 2; is2%beg = is2%beg + weno_polyn
            is2%end = is2%end - weno_polyn
        else
            weno_dir = 3; is3%beg = is3%beg + weno_polyn
            is3%end = is3%end - weno_polyn
        end if

        call s_weno(v_vf(iv%beg:iv%end), &
                    vL_qp(0, 0)%vf(iv%beg:iv%end), &
                    vR_qp(0, 0)%vf(iv%beg:iv%end), &
                    cd_vars, norm_dir, weno_dir, 1, &
                    is1, is2, is3)
        ! ==================================================================

        ! Reconstruction in s2-direction ===================================
        if (split_err .neqv. .true.) return

        if (norm_dir /= 1) then
            weno_dir = 1; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn
        else
            weno_dir = 2; is2%beg = is2%beg + weno_polyn
            is2%end = is2%end - weno_polyn
        end if

        call s_weno(vL_qp(0, 0)%vf(iv%beg:iv%end), &
                    vL_qp(-1, 0)%vf(iv%beg:iv%end), &
                    vL_qp(1, 0)%vf(iv%beg:iv%end), &
                    dflt_int, norm_dir, weno_dir, 2, &
                    is1, is2, is3)
        call s_weno(vR_qp(0, 0)%vf(iv%beg:iv%end), &
                    vR_qp(-1, 0)%vf(iv%beg:iv%end), &
                    vR_qp(1, 0)%vf(iv%beg:iv%end), &
                    dflt_int, norm_dir, weno_dir, 2, &
                    is1, is2, is3)
        ! ==================================================================

        ! Reconstruction in s3-direction ===================================
        if (p == 0) return

        if (norm_dir /= 3) then
            weno_dir = 3; is3%beg = is3%beg + weno_polyn
            is3%end = is3%end - weno_polyn
        else
            weno_dir = 2; is2%beg = is2%beg + weno_polyn
            is2%end = is2%end - weno_polyn
        end if

        call s_weno(vL_qp(-1, 0)%vf(iv%beg:iv%end), &
                    vL_qp(-1, -1)%vf(iv%beg:iv%end), &
                    vL_qp(-1, 1)%vf(iv%beg:iv%end), &
                    dflt_int, norm_dir, weno_dir, 2, &
                    is1, is2, is3)
        call s_weno(vL_qp(1, 0)%vf(iv%beg:iv%end), &
                    vL_qp(1, -1)%vf(iv%beg:iv%end), &
                    vL_qp(1, 1)%vf(iv%beg:iv%end), &
                    dflt_int, norm_dir, weno_dir, 2, &
                    is1, is2, is3)
        call s_weno(vR_qp(-1, 0)%vf(iv%beg:iv%end), &
                    vR_qp(-1, -1)%vf(iv%beg:iv%end), &
                    vR_qp(-1, 1)%vf(iv%beg:iv%end), &
                    dflt_int, norm_dir, weno_dir, 2, &
                    is1, is2, is3)
        call s_weno(vR_qp(1, 0)%vf(iv%beg:iv%end), &
                    vR_qp(1, -1)%vf(iv%beg:iv%end), &
                    vR_qp(1, 1)%vf(iv%beg:iv%end), &
                    dflt_int, norm_dir, weno_dir, 2, &
                    is1, is2, is3)
        ! ==================================================================

    end subroutine s_reconstruct_cell_boundary_values ! --------------------

    !>  The goal of the procedure is to numerically approximate
        !!      the left or right cell-boundary integral-average of the
        !!      given variables by getting the arithmetic mean of their
        !!      WENO-reconstructed values at the cell-boundary Gaussian
        !!      quadrature points.
        !!  @param vK_qp The inputted WENO-reconstructed values of cell-averaged variables
        !!      at the left or right cell-boundary Gaussian quadrature points and
        !!      their numerically approximated cell-boundary integral-average.
    subroutine s_average_cell_boundary_values(vK_qp) ! ---------------------

        type(vector_field), &
            dimension(ichi%beg:ichi%end, &
                      ipsi%beg:ipsi%end), &
            intent(INOUT) :: vK_qp

        integer :: i, j, k !< Generic loop iterators

        if (split_err .neqv. .true.) return

        do k = iv%beg, iv%end

            vK_qp(0, 0)%vf(k)%sf = 0d0

            do j = ipsi%beg, ipsi%end, 2
                do i = ichi%beg, ichi%end, 2
                    vK_qp(0, 0)%vf(k)%sf = vK_qp(0, 0)%vf(k)%sf &
                                           + vK_qp(i, j)%vf(k)%sf
                end do
            end do

            vK_qp(0, 0)%vf(k)%sf = vK_qp(0, 0)%vf(k)%sf/(2**(num_dims - 1))

        end do

    end subroutine s_average_cell_boundary_values ! ------------------------

    !>  The purpose of this subroutine is to employ the inputted
        !!      left and right cell-boundary integral-averaged variables
        !!      to compute the relevant cell-average first-order spatial
        !!      derivatives in the x-, y- or z-direction by means of the
        !!      scalar divergence theorem.
        !!  @param vL_vf Left cell-boundary integral averages
        !!  @param vR_vf Right cell-boundary integral averages
        !!  @param dv_ds_vf Cell-average first-order spatial derivatives
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_apply_scalar_divergence_theorem(vL_vf, vR_vf, & ! --------
                                                 dv_ds_vf, &
                                                 norm_dir)

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(IN) :: vL_vf, vR_vf

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(INOUT) :: dv_ds_vf

        integer, intent(IN) :: norm_dir

        integer :: i, j, k, l !< Generic loop iterators

        ! First-Order Spatial Derivatives in x-direction ===================
        if (norm_dir == 1) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.
            do i = iv%beg, iv%end
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg + 1, ix%end - 1
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1d0/((1d0 + wa_flg)*dx(j)) &
                                *(wa_flg*vL_vf(i)%sf(j + 1, k, l) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j - 1, k, l))
                        end do
                    end do
                end do
            end do

            ! Modifying the way the scalar divergence theorem is utilized
            ! near beginning and end Riemann state extrapolation boundary
            ! conditions, respectively, as to omit any contributions from
            ! outside the physical computational domain in the calculated
            ! cell-averaged first-order spatial derivatives
            if (weno_avg .neqv. .true.) return

            if (bc_x%beg == -4) then     ! Riemann state extrap. BC at beg.
                do i = iv%beg, iv%end
                    do l = iz%beg, iz%end
                        do k = iy%beg, iy%end
                            dv_ds_vf(i)%sf(0, k, l) = 1d0/dx(0) &
                                                      *(vR_vf(i)%sf(0, k, l) &
                                                        - vL_vf(i)%sf(0, k, l))
                        end do
                    end do
                end do
            end if

            if (bc_x%end == -4) then     ! Riemann state extrap. BC at end
                do i = iv%beg, iv%end
                    do l = iz%beg, iz%end
                        do k = iy%beg, iy%end
                            dv_ds_vf(i)%sf(m, k, l) = 1d0/dx(m) &
                                                      *(vR_vf(i)%sf(m, k, l) &
                                                        - vL_vf(i)%sf(m, k, l))
                        end do
                    end do
                end do
            end if
            ! END: First-Order Spatial Derivatives in x-direction ==============

            ! First-Order Spatial Derivatives in y-direction ===================
        elseif (norm_dir == 2) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.
            do i = iv%beg, iv%end
                do l = iz%beg, iz%end
                    do k = iy%beg + 1, iy%end - 1
                        do j = ix%beg, ix%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1d0/((1d0 + wa_flg)*dy(k)) &
                                *(wa_flg*vL_vf(i)%sf(j, k + 1, l) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j, k - 1, l))
                        end do
                    end do
                end do
            end do

            ! Modifying the way the scalar divergence theorem is utilized
            ! near beginning and end Riemann state extrapolation boundary
            ! conditions, respectively, as to omit any contributions from
            ! outside the physical computational domain in the calculated
            ! cell-averaged first-order spatial derivatives
            if (weno_avg .neqv. .true.) return

            if (bc_y%beg == -4) then     ! Riemann state extrap. BC at beg.
                do i = iv%beg, iv%end
                    do l = iz%beg, iz%end
                        do j = ix%beg, ix%end
                            dv_ds_vf(i)%sf(j, 0, l) = 1d0/dy(0) &
                                                      *(vR_vf(i)%sf(j, 0, l) &
                                                        - vL_vf(i)%sf(j, 0, l))
                        end do
                    end do
                end do
            end if

            if (bc_y%end == -4) then     ! Riemann state extrap. BC at end
                do i = iv%beg, iv%end
                    do l = iz%beg, iz%end
                        do j = ix%beg, ix%end
                            dv_ds_vf(i)%sf(j, n, l) = 1d0/dy(n) &
                                                      *(vR_vf(i)%sf(j, n, l) &
                                                        - vL_vf(i)%sf(j, n, l))
                        end do
                    end do
                end do
            end if
            ! END: First-Order Spatial Derivatives in y-direction ==============

            ! First-Order Spatial Derivatives in z-direction ===================
        else

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.
            do i = iv%beg, iv%end
                do l = iz%beg + 1, iz%end - 1
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1d0/((1d0 + wa_flg)*dz(l)) &
                                *(wa_flg*vL_vf(i)%sf(j, k, l + 1) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j, k, l - 1))
                        end do
                    end do
                end do
            end do

            ! Modifying the way the scalar divergence theorem is utilized
            ! near beginning and end Riemann state extrapolation boundary
            ! conditions, respectively, as to omit any contributions from
            ! outside the physical computational domain in the calculated
            ! cell-averaged first-order spatial derivatives
            if (weno_avg .neqv. .true.) return

            if (bc_z%beg == -4) then     ! Riemann state extrap. BC at beg.
                do i = iv%beg, iv%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            dv_ds_vf(i)%sf(j, k, 0) = 1d0/dz(0) &
                                                      *(vR_vf(i)%sf(j, k, 0) &
                                                        - vL_vf(i)%sf(j, k, 0))
                        end do
                    end do
                end do
            end if

            if (bc_z%end == -4) then     ! Riemann state extrap. BC at end
                do i = iv%beg, iv%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            dv_ds_vf(i)%sf(j, k, p) = 1d0/dz(p) &
                                                      *(vR_vf(i)%sf(j, k, p) &
                                                        - vL_vf(i)%sf(j, k, p))
                        end do
                    end do
                end do
            end if

        end if
        ! END: First-Order Spatial Derivatives in z-direction ==============

    end subroutine s_apply_scalar_divergence_theorem ! ---------------------

    !>  The goal of this procedure is to utilize the inputted
        !!      left and right cell-boundary integral-averaged vector
        !!      components in the x-, y-, and z-directions to compute
        !!      the vector divergence by using the divergence theorem.
        !!  @param vL_x_ndqp Left cell-boundary integral-average x-dir component
        !!  @param vL_y_ndqp Left cell-boundary integral-average y-dir component
        !!  @param vL_z_ndqp Left cell-boundary integral-average z-dir component
        !!  @param vR_x_ndqp Right cell-boundary integral-average x-dir component
        !!  @param vR_y_ndqp Right cell-boundary integral-average y-dir component
        !!  @param vR_z_ndqp Right cell-boundary integral-average z-dir component
        !!  @param div_v_vf Cell-average divergence
    subroutine s_apply_vector_divergence_theorem( & ! ----------------
        vL_x_ndqp, vL_y_ndqp, vL_z_ndqp, &
        vR_x_ndqp, vR_y_ndqp, vR_z_ndqp, &
        div_v_vf)

        type(vector_field), &
            dimension(1:num_dims, &
                      ichi%beg:ichi%end, &
                      ipsi%beg:ipsi%end), &
            intent(IN) :: vL_x_ndqp, vR_x_ndqp, &
                          vL_y_ndqp, vR_y_ndqp, &
                          vL_z_ndqp, vR_z_ndqp

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(INOUT) :: div_v_vf

        integer :: i, j, k, l !< Generic loop iterators

        ! First-Order Spatial Derivatives in x-direction ===================

        ! General application of the vector divergence theorem which uses
        ! the left and right cell-boundary integral-averages, inside each
        ! cell, or an arithmetic mean of these two at the cell-boundaries,
        ! in order to obtain cell-average first-order spatial derivatives
        ! inside the cell
        do i = iv%beg, iv%end
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg + 1, ix%end - 1
                        div_v_vf(i)%sf(j, k, l) = 1d0/((1d0 + wa_flg)*dx(j)) &
                                                  *(wa_flg*vL_x_ndqp(1, 0, 0)%vf(i)%sf(j + 1, k, l) &
                                                    + vR_x_ndqp(1, 0, 0)%vf(i)%sf(j, k, l) &
                                                    - vL_x_ndqp(1, 0, 0)%vf(i)%sf(j, k, l) &
                                                    - wa_flg*vR_x_ndqp(1, 0, 0)%vf(i)%sf(j - 1, k, l))
                    end do
                end do
            end do
        end do

        ! Customizing the way the vector divergence theorem is applied near
        ! beginning and end Riemann state extrapolation boundary conditions,
        ! respectively, as to avoid contributions from outside the physical
        ! computational domain in the computation of the cell-average first-
        ! order spatial derivatives
        if (weno_avg) then

            if (bc_x%beg == -4) then     ! Riemann state extrap. BC at beg.
                do i = iv%beg, iv%end
                    do l = iz%beg, iz%end
                        do k = iy%beg, iy%end
                            div_v_vf(i)%sf(0, k, l) = &
                                div_v_vf(i)%sf(0, k, l) - 5d-1/dx(0) &
                                *(vL_x_ndqp(1, 0, 0)%vf(i)%sf(1, k, l) &
                                  - vR_x_ndqp(1, 0, 0)%vf(i)%sf(0, k, l) &
                                  + vL_x_ndqp(1, 0, 0)%vf(i)%sf(0, k, l) &
                                  - vR_x_ndqp(1, 0, 0)%vf(i)%sf(-1, k, l))
                        end do
                    end do
                end do
            end if

            if (bc_x%end == -4) then     ! Riemann state extrap. BC at end
                do i = iv%beg, iv%end
                    do l = iz%beg, iz%end
                        do k = iy%beg, iy%end
                            div_v_vf(i)%sf(m, k, l) = &
                                div_v_vf(i)%sf(m, k, l) - 5d-1/dx(m) &
                                *(vL_x_ndqp(1, 0, 0)%vf(i)%sf(m + 1, k, l) &
                                  - vR_x_ndqp(1, 0, 0)%vf(i)%sf(m, k, l) &
                                  + vL_x_ndqp(1, 0, 0)%vf(i)%sf(m, k, l) &
                                  - vR_x_ndqp(1, 0, 0)%vf(i)%sf(m - 1, k, l))
                        end do
                    end do
                end do
            end if

        end if

        ! END: First-Order Spatial Derivatives in x-direction ==============

        ! First-Order Spatial Derivatives in y-direction ===================

        ! General application of the vector divergence theorem which uses
        ! the left and right cell-boundary integral-averages, inside each
        ! cell, or an arithmetic mean of these two at the cell-boundaries,
        ! in order to obtain cell-average first-order spatial derivatives
        ! inside the cell
        if (n == 0) return

        do i = iv%beg, iv%end
            do l = iz%beg, iz%end
                do k = iy%beg + 1, iy%end - 1
                    do j = ix%beg, ix%end
                        div_v_vf(i)%sf(j, k, l) = div_v_vf(i)%sf(j, k, l) &
                                                  + 1d0/((1d0 + wa_flg)*dy(k)) &
                                                  *(wa_flg*vL_y_ndqp(2, 0, 0)%vf(i)%sf(j, k + 1, l) &
                                                    + vR_y_ndqp(2, 0, 0)%vf(i)%sf(j, k, l) &
                                                    - vL_y_ndqp(2, 0, 0)%vf(i)%sf(j, k, l) &
                                                    - wa_flg*vR_y_ndqp(2, 0, 0)%vf(i)%sf(j, k - 1, l))
                    end do
                end do
            end do
        end do

        ! Customizing the way the vector divergence theorem is applied near
        ! beginning and end Riemann state extrapolation boundary conditions,
        ! respectively, as to avoid contributions from outside the physical
        ! computational domain in the computation of the cell-average first-
        ! order spatial derivatives
        if (weno_avg) then

            if (bc_y%beg == -4) then     ! Riemann state extrap. BC at beg.
                do i = iv%beg, iv%end
                    do l = iz%beg, iz%end
                        do j = ix%beg, ix%end
                            div_v_vf(i)%sf(j, 0, l) = &
                                div_v_vf(i)%sf(j, 0, l) - 5d-1/dy(0) &
                                *(vL_y_ndqp(2, 0, 0)%vf(i)%sf(j, 1, l) &
                                  - vR_y_ndqp(2, 0, 0)%vf(i)%sf(j, 0, l) &
                                  + vL_y_ndqp(2, 0, 0)%vf(i)%sf(j, 0, l) &
                                  - vR_y_ndqp(2, 0, 0)%vf(i)%sf(j, -1, l))
                        end do
                    end do
                end do
            end if

            if (bc_y%end == -4) then     ! Riemann state extrap. BC at end
                do i = iv%beg, iv%end
                    do l = iz%beg, iz%end
                        do j = ix%beg, ix%end
                            div_v_vf(i)%sf(j, n, l) = &
                                div_v_vf(i)%sf(j, n, l) - 5d-1/dy(n) &
                                *(vL_y_ndqp(2, 0, 0)%vf(i)%sf(j, n + 1, l) &
                                  - vR_y_ndqp(2, 0, 0)%vf(i)%sf(j, n, l) &
                                  + vL_y_ndqp(2, 0, 0)%vf(i)%sf(j, n, l) &
                                  - vR_y_ndqp(2, 0, 0)%vf(i)%sf(j, n - 1, l))
                        end do
                    end do
                end do
            end if

        end if

        ! END: First-Order Spatial Derivatives in y-direction ==============

        ! First-Order Spatial Derivatives in z-direction ===================

        ! General application of the vector divergence theorem which uses
        ! the left and right cell-boundary integral-averages, inside each
        ! cell, or an arithmetic mean of these two at the cell-boundaries,
        ! in order to obtain cell-average first-order spatial derivatives
        ! inside the cell
        if (p == 0) return

        do i = iv%beg, iv%end
            do l = iz%beg + 1, iz%end - 1
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        div_v_vf(i)%sf(j, k, l) = div_v_vf(i)%sf(j, k, l) &
                                                  + 1d0/((1d0 + wa_flg)*dz(l)) &
                                                  *(wa_flg*vL_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, l + 1) &
                                                    + vR_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, l) &
                                                    - vL_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, l) &
                                                    - wa_flg*vR_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, l - 1))
                    end do
                end do
            end do
        end do

        ! Customizing the way the vector divergence theorem is applied near
        ! beginning and end Riemann state extrapolation boundary conditions,
        ! respectively, as to avoid contributions from outside the physical
        ! computational domain in the computation of the cell-average first-
        ! order spatial derivatives
        if (weno_avg) then

            if (bc_z%beg == -4) then     ! Riemann state extrap. BC at beg.
                do i = iv%beg, iv%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            div_v_vf(i)%sf(j, k, 0) = &
                                div_v_vf(i)%sf(j, k, 0) - 5d-1/dz(0) &
                                *(vL_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, 1) &
                                  - vR_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, 0) &
                                  + vL_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, 0) &
                                  - vR_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, -1))
                        end do
                    end do
                end do
            end if

            if (bc_z%end == -4) then     ! Riemann state extrap. BC at end
                do i = iv%beg, iv%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            div_v_vf(i)%sf(j, k, p) = &
                                div_v_vf(i)%sf(j, k, p) - 5d-1/dz(p) &
                                *(vL_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, p + 1) &
                                  - vR_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, p) &
                                  + vL_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, p) &
                                  - vR_z_ndqp(3, 0, 0)%vf(i)%sf(j, k, p - 1))
                        end do
                    end do
                end do
            end if

        end if

        ! END: First-Order Spatial Derivatives in z-direction ==============

    end subroutine s_apply_vector_divergence_theorem ! ---------------------

    !>  The purpose of the procedure is to utilize the inputted
        !!      cell-averaged first-order spatial derivatives in the x-,
        !!      y- and z-directions to calculate the gradient magnitude.
        !!  @param dv_dx_vf Cell-average first-order spatial derivatives in the x-dir
        !!  @param dv_dy_vf Cell-average first-order spatial derivatives in the y-dir
        !!  @param dv_dz_vf Cell-average first-order spatial derivatives in the z-dir
        !!  @param gm_v_vf  Gradient magnitude
    subroutine s_compute_gradient_magnitude(dv_dx_vf, & ! -----------------
                                            dv_dy_vf, &
                                            dv_dz_vf, &
                                            gm_v_vf)

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(IN) :: dv_dx_vf, &
                          dv_dy_vf, &
                          dv_dz_vf

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(INOUT) :: gm_v_vf

        integer :: i, j, k, l !< Generic loop iterators

        ! Scalar Product Contribution in x-direction =======================
        do i = iv%beg, iv%end
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        gm_v_vf(i)%sf(j, k, l) = dv_dx_vf(i)%sf(j, k, l) &
                                                 *dv_dx_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        ! ==================================================================

        ! Scalar Product Contribution in y-direction =======================
        if (n > 0) then

            do i = iv%beg, iv%end
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            gm_v_vf(i)%sf(j, k, l) = gm_v_vf(i)%sf(j, k, l) &
                                                     + dv_dy_vf(i)%sf(j, k, l) &
                                                     *dv_dy_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            ! ==================================================================

            ! Scalar Product Contribution in z-direction =======================
            if (p > 0) then
                do i = iv%beg, iv%end
                    do l = iz%beg, iz%end
                        do k = iy%beg, iy%end
                            do j = ix%beg, ix%end
                                gm_v_vf(i)%sf(j, k, l) = gm_v_vf(i)%sf(j, k, l) &
                                                         + dv_dz_vf(i)%sf(j, k, l) &
                                                         *dv_dz_vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
            end if

        end if
        ! ==================================================================

        ! Square Root of the Scalar Product ================================
        do i = iv%beg, iv%end
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        gm_v_vf(i)%sf(j, k, l) = sqrt(gm_v_vf(i)%sf(j, k, l))
                    end do
                end do
            end do
        end do
        ! ==================================================================

    end subroutine s_compute_gradient_magnitude ! --------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_rhs_module() ! -----------------------------------

        integer :: i, j, k, l !< Generic loop iterators

        ! Deallocation/Disassociation of q_cons_qp and q_prim_qp ===========
        do k = itau%beg, itau%end
            do j = iksi%beg, iksi%end
                do i = ieta%beg, ieta%end

                    if (i /= 0 .and. abs(j) >= abs(k)) then

                        do l = 1, sys_size
                            deallocate (q_cons_qp(i, j, k)%vf(l)%sf)
                        end do

                        do l = mom_idx%beg, E_idx
                            deallocate (q_prim_qp(i, j, k)%vf(l)%sf)
                        end do

                        if (model_eqns == 3) then
                            do l = internalEnergies_idx%beg, internalEnergies_idx%end
                                deallocate (q_prim_qp(i, j, k)%vf(l)%sf)
                            end do
                        end if

                        do l = 1, sys_size
                            nullify (q_prim_qp(i, j, k)%vf(l)%sf)
                        end do

                    end if

                    deallocate (q_cons_qp(i, j, k)%vf, q_prim_qp(i, j, k)%vf)

                end do
            end do
        end do

        deallocate (q_cons_qp, q_prim_qp)
        ! ==================================================================

        ! Deallocation/Disassociation of qK_cons_ndqp and qK_prim_ndqp =====
        do k = ipsi%beg, ipsi%end
            do j = ichi%beg, ichi%end
                do i = num_dims, 1, -1

                    if (abs(j) >= abs(k)) then

                        do l = 1, cont_idx%end
                            nullify (qL_prim_ndqp(i, j, k)%vf(l)%sf)
                            nullify (qR_prim_ndqp(i, j, k)%vf(l)%sf)
                        end do

                        do l = adv_idx%beg, adv_idx%end
                            nullify (qL_prim_ndqp(i, j, k)%vf(l)%sf)
                            nullify (qR_prim_ndqp(i, j, k)%vf(l)%sf)
                        end do

                        if (i /= 1) then

                            if ((char_decomp .neqv. .true.) &
                                .and. &
                                any(Re_size > 0)) then
                                if (weno_vars == 1) then
                                    do l = 1, mom_idx%end
                                        deallocate (qL_cons_ndqp(i, j, k)%vf(l)%sf)
                                        deallocate (qR_cons_ndqp(i, j, k)%vf(l)%sf)
                                    end do
                                else
                                    do l = mom_idx%beg, mom_idx%end
                                        deallocate (qL_prim_ndqp(i, j, k)%vf(l)%sf)
                                        deallocate (qR_prim_ndqp(i, j, k)%vf(l)%sf)
                                    end do
                                    if (model_eqns == 3) then
                                        do l = internalEnergies_idx%beg, internalEnergies_idx%end
                                            deallocate (qL_prim_ndqp(i, j, k)%vf(l)%sf)
                                            deallocate (qR_prim_ndqp(i, j, k)%vf(l)%sf)
                                        end do
                                    end if
                                end if
                            end if

                            do l = 1, crv_size
                                deallocate (qL_cons_ndqp(i, j, k)%vf( &
                                            E_idx + crv_idx(l))%sf)
                                deallocate (qR_cons_ndqp(i, j, k)%vf( &
                                            E_idx + crv_idx(l))%sf)
                            end do

                            do l = 1, sys_size
                                nullify (qL_cons_ndqp(i, j, k)%vf(l)%sf)
                                nullify (qR_cons_ndqp(i, j, k)%vf(l)%sf)
                                nullify (qL_prim_ndqp(i, j, k)%vf(l)%sf)
                                nullify (qR_prim_ndqp(i, j, k)%vf(l)%sf)
                            end do

                        else

                            do l = 1, cont_idx%end
                                deallocate (qL_cons_ndqp(i, j, k)%vf(l)%sf)
                                deallocate (qR_cons_ndqp(i, j, k)%vf(l)%sf)
                            end do

                            if (weno_vars == 1) then
                                do l = mom_idx%beg, E_idx
                                    deallocate (qL_cons_ndqp(i, j, k)%vf(l)%sf)
                                    deallocate (qR_cons_ndqp(i, j, k)%vf(l)%sf)
                                end do
                            end if

                            do l = mom_idx%beg, E_idx
                                deallocate (qL_prim_ndqp(i, j, k)%vf(l)%sf)
                                deallocate (qR_prim_ndqp(i, j, k)%vf(l)%sf)
                            end do
                            if (model_eqns == 3) then
                                do l = internalEnergies_idx%beg, internalEnergies_idx%end
                                    deallocate (qL_prim_ndqp(i, j, k)%vf(l)%sf)
                                    deallocate (qR_prim_ndqp(i, j, k)%vf(l)%sf)
                                end do
                            end if

                            do l = adv_idx%beg, adv_idx%end
                                deallocate (qL_cons_ndqp(i, j, k)%vf(l)%sf)
                                deallocate (qR_cons_ndqp(i, j, k)%vf(l)%sf)
                            end do

                        end if

                    end if

                    deallocate (qL_cons_ndqp(i, j, k)%vf, qL_prim_ndqp(i, j, k)%vf)
                    deallocate (qR_cons_ndqp(i, j, k)%vf, qR_prim_ndqp(i, j, k)%vf)

                end do
            end do
        end do

        deallocate (qL_cons_ndqp, qR_cons_ndqp, qL_prim_ndqp, qR_prim_ndqp)
        ! END: Deallocation/Disassociation of qK_cons_ndqp and qK_prim_ndqp

        ! Deallocation of dq_prim_ds_qp ====================================
        if (any(Re_size > 0) .or. We_size > 0 .or. hypoelasticity) then

            do k = itau%beg, itau%end
                do j = iksi%beg, iksi%end
                    do i = ieta%beg, ieta%end

                        if (all((/i, j, k/) == 0) .and. (any(Re_size > 0) .or. hypoelasticity)) then

                            do l = mom_idx%beg, mom_idx%end
                                deallocate (dq_prim_dx_qp(i, j, k)%vf(l)%sf)
                                deallocate (gm_vel_qp(i, j, k)%vf(l)%sf)
                            end do

                            if (n > 0) then

                                do l = mom_idx%beg, mom_idx%end
                                    deallocate (dq_prim_dy_qp(i, j, k)%vf(l)%sf)
                                end do

                                if (p > 0) then
                                    do l = mom_idx%beg, mom_idx%end
                                        deallocate (dq_prim_dz_qp(i, j, k)%vf(l)%sf)
                                    end do
                                end if

                            end if

                        end if

                        if ((i /= 0 .and. abs(j) >= abs(k)) &
                            .or. &
                            all((/i, j, k/) == 0)) then
                            do l = 1, crv_size
                                deallocate (dq_prim_dx_qp(i, j, k)%vf( &
                                            E_idx + crv_idx(l))%sf)
                                if (n > 0) then
                                    deallocate (dq_prim_dy_qp(i, j, k)%vf( &
                                                E_idx + crv_idx(l))%sf)
                                    if (p > 0) then
                                        deallocate (dq_prim_dz_qp(i, j, k)%vf( &
                                                    E_idx + crv_idx(l))%sf)
                                    end if
                                end if
                            end do
                        end if

                        deallocate (dq_prim_dx_qp(i, j, k)%vf)
                        deallocate (dq_prim_dy_qp(i, j, k)%vf)
                        deallocate (dq_prim_dz_qp(i, j, k)%vf)
                        deallocate (gm_vel_qp(i, j, k)%vf)

                    end do
                end do
            end do

            deallocate (dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp, gm_vel_qp)

        end if
        ! END: Deallocation of dq_prim_ds_qp ===============================

        ! Deallocation/Disassociation of dqK_prim_ds_ndqp ==================
        if (any(Re_size > 0) .or. We_size > 0) then
            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end
                    do i = num_dims, 1, -1

                        if (abs(j) >= abs(k)) then

                            if (any(Re_size > 0)) then

                                do l = mom_idx%beg, mom_idx%end
                                    deallocate (dqL_prim_dx_ndqp(i, j, k)%vf(l)%sf)
                                    deallocate (dqR_prim_dx_ndqp(i, j, k)%vf(l)%sf)
                                end do

                                if (n > 0) then
                                    do l = mom_idx%beg, mom_idx%end
                                        deallocate (dqL_prim_dy_ndqp(i, j, k)%vf(l)%sf)
                                        deallocate (dqR_prim_dy_ndqp(i, j, k)%vf(l)%sf)
                                    end do
                                end if

                                if (p > 0) then
                                    do l = mom_idx%beg, mom_idx%end
                                        deallocate (dqL_prim_dz_ndqp(i, j, k)%vf(l)%sf)
                                        deallocate (dqR_prim_dz_ndqp(i, j, k)%vf(l)%sf)
                                    end do
                                end if

                            end if

                            do l = 1, crv_size
                                deallocate (dqL_prim_dx_ndqp(i, j, k)%vf( &
                                            E_idx + crv_idx(l))%sf)
                                deallocate (dqR_prim_dx_ndqp(i, j, k)%vf( &
                                            E_idx + crv_idx(l))%sf)
                                if (n > 0) then
                                    deallocate (dqL_prim_dy_ndqp(i, j, k)%vf( &
                                                E_idx + crv_idx(l))%sf)
                                    deallocate (dqR_prim_dy_ndqp(i, j, k)%vf( &
                                                E_idx + crv_idx(l))%sf)
                                end if
                            end do

                            if (p > 0) then
                                do l = 1, crv_size
                                    deallocate (dqL_prim_dz_ndqp(i, j, k)%vf( &
                                                E_idx + crv_idx(l))%sf)
                                    deallocate (dqR_prim_dz_ndqp(i, j, k)%vf( &
                                                E_idx + crv_idx(l))%sf)
                                end do
                            end if

                        end if

                        deallocate (dqL_prim_dx_ndqp(i, j, k)%vf)
                        deallocate (dqL_prim_dy_ndqp(i, j, k)%vf)
                        deallocate (dqL_prim_dz_ndqp(i, j, k)%vf)
                        deallocate (dqR_prim_dx_ndqp(i, j, k)%vf)
                        deallocate (dqR_prim_dy_ndqp(i, j, k)%vf)
                        deallocate (dqR_prim_dz_ndqp(i, j, k)%vf)

                    end do
                end do
            end do
        end if

        deallocate (dqL_prim_dx_ndqp, dqL_prim_dy_ndqp, dqL_prim_dz_ndqp)
        deallocate (dqR_prim_dx_ndqp, dqR_prim_dy_ndqp, dqR_prim_dz_ndqp)
        ! END: Deallocation/Disassociation of dqK_prim_ds_ndqp =============

        ! Deallocation of gm_alpha_qp ======================================
        if (We_size > 0) then
            do k = itau%beg, itau%end
                do j = iksi%beg, iksi%end
                    do i = ieta%beg, ieta%end

                        if (any(sum(abs((/i, j, k/))) == (/0, num_dims/))) then
                            do l = 1, crv_size
                                deallocate (gm_alpha_qp(i, j, k)%vf(crv_idx(l))%sf)
                            end do
                        end if

                        deallocate (gm_alpha_qp(i, j, k)%vf)

                    end do
                end do
            end do
        end if

        deallocate (gm_alpha_qp)
        ! ==================================================================

        ! Deallocation of gm_alphaK_ndqp ===================================
        if (We_size > 0) then
            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end
                    do i = 1, num_dims

                        if (any(sum(abs((/j, k/))) == (/0, num_dims/))) then
                            do l = 1, crv_size
                                deallocate (gm_alphaL_ndqp(i, j, k)%vf( &
                                            crv_idx(l))%sf)
                                deallocate (gm_alphaR_ndqp(i, j, k)%vf( &
                                            crv_idx(l))%sf)
                            end do
                        end if

                        deallocate (gm_alphaL_ndqp(i, j, k)%vf)
                        deallocate (gm_alphaR_ndqp(i, j, k)%vf)

                    end do
                end do
            end do
        end if

        deallocate (gm_alphaL_ndqp, gm_alphaR_ndqp)
        ! ==================================================================

        ! Deallocation of un_alphaK_s_ndqp =================================
        if (We_size > 0) then

            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end

                    if (abs(j) >= abs(k)) then

                        do l = 1, crv_size
                            deallocate (un_alphaL_x_ndqp(1, j, k)%vf(crv_idx(l))%sf)
                            deallocate (un_alphaR_x_ndqp(1, j, k)%vf(crv_idx(l))%sf)
                            if (n > 0) then
                                deallocate (un_alphaL_y_ndqp(2, j, k)%vf(crv_idx(l))%sf)
                                deallocate (un_alphaR_y_ndqp(2, j, k)%vf(crv_idx(l))%sf)
                            end if
                        end do

                        if (p > 0) then
                            do l = 1, crv_size
                                deallocate (un_alphaL_z_ndqp(3, j, k)%vf( &
                                            crv_idx(l))%sf)
                                deallocate (un_alphaR_z_ndqp(3, j, k)%vf( &
                                            crv_idx(l))%sf)
                            end do
                        end if

                    end if

                    do i = 1, num_dims
                        deallocate (un_alphaL_x_ndqp(i, j, k)%vf)
                        deallocate (un_alphaL_y_ndqp(i, j, k)%vf)
                        deallocate (un_alphaL_z_ndqp(i, j, k)%vf)
                        deallocate (un_alphaR_x_ndqp(i, j, k)%vf)
                        deallocate (un_alphaR_y_ndqp(i, j, k)%vf)
                        deallocate (un_alphaR_z_ndqp(i, j, k)%vf)
                    end do

                end do
            end do

            deallocate (un_alphaL_x_ndqp, un_alphaL_y_ndqp, un_alphaL_z_ndqp)
            deallocate (un_alphaR_x_ndqp, un_alphaR_y_ndqp, un_alphaR_z_ndqp)

        end if
        ! ==================================================================

        ! Deallocation of kappa_vf =========================================
        if (We_size > 0) then

            do i = 1, crv_size
                deallocate (kappa_vf(crv_idx(i))%sf)
            end do

            deallocate (kappa_vf)

        end if
        ! ==================================================================

        ! Deallocation/Disassociation of kappaK_ndqp =======================
        if (We_size > 0) then
            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end
                    do i = num_dims, 1, -1

                        if (abs(j) >= abs(k)) then
                            if (i /= 1) then
                                do l = 1, crv_size
                                    nullify (kappaL_ndqp(i, j, k)%vf(crv_idx(l))%sf)
                                    nullify (kappaR_ndqp(i, j, k)%vf(crv_idx(l))%sf)
                                end do
                            else
                                do l = 1, crv_size
                                    deallocate (kappaL_ndqp(i, j, k)%vf( &
                                                crv_idx(l))%sf)
                                    deallocate (kappaR_ndqp(i, j, k)%vf( &
                                                crv_idx(l))%sf)
                                end do
                            end if
                        end if

                        deallocate (kappaL_ndqp(i, j, k)%vf, kappaR_ndqp(i, j, k)%vf)

                    end do
                end do
            end do
        end if

        deallocate (kappaL_ndqp, kappaR_ndqp)
        ! ==================================================================

        if (alt_crv) then
            ! Deallocation of dgm_alpha_ds_qp & un_alpha_s_qp
            if (We_size > 0) then
                do k = itau%beg, itau%end
                    do j = iksi%beg, iksi%end
                        do i = ieta%beg, ieta%end

                            if (any(sum(abs((/i, j, k/))) == (/0, num_dims/))) then
                                do l = 1, crv_size
                                    deallocate (dgm_alpha_dx_qp(i, j, k)%vf(crv_idx(l))%sf)
                                    deallocate (dgm_alpha_dy_qp(i, j, k)%vf(crv_idx(l))%sf)
                                    deallocate (dgm_alpha_dz_qp(i, j, k)%vf(crv_idx(l))%sf)
                                    deallocate (un_alpha_x_qp(i, j, k)%vf(crv_idx(l))%sf)
                                    deallocate (un_alpha_y_qp(i, j, k)%vf(crv_idx(l))%sf)
                                    deallocate (un_alpha_z_qp(i, j, k)%vf(crv_idx(l))%sf)
                                end do
                            end if

                            deallocate (dgm_alpha_dx_qp(i, j, k)%vf)
                            deallocate (dgm_alpha_dy_qp(i, j, k)%vf)
                            deallocate (dgm_alpha_dz_qp(i, j, k)%vf)
                            deallocate (un_alpha_x_qp(i, j, k)%vf)
                            deallocate (un_alpha_y_qp(i, j, k)%vf)
                            deallocate (un_alpha_z_qp(i, j, k)%vf)

                        end do
                    end do
                end do
            end if

            deallocate (dgm_alpha_dx_qp)
            deallocate (dgm_alpha_dy_qp)
            deallocate (dgm_alpha_dz_qp)
            deallocate (un_alpha_x_qp)
            deallocate (un_alpha_y_qp)
            deallocate (un_alpha_z_qp)
        end if

        ! Deallocation of alt_kappa_vf & laplacian_vf
        if (We_size > 0) then
            if (alt_crv) then
                do i = 1, crv_size
                    deallocate (alt_kappa_vf(crv_idx(i))%sf)
                    deallocate (laplacian_vf(crv_idx(i))%sf)
                end do

                deallocate (alt_kappa_vf)
                deallocate (laplacian_vf)
            end if

            if (We_src) then
                deallocate (We_mtm_src)
                deallocate (We_nrg_src)
            end if

        end if

        if (We_size > 0 .and. We_rhs_flux) then
            do i = 1, num_dims
                deallocate (tau_We_vf(cont_idx%end + i)%sf)
            end do
            deallocate (tau_We_vf(E_idx)%sf)

            deallocate (tau_We_vf)
        end if

        if (any(Re_size > 0) .and. cyl_coord) then
            do i = 1, num_dims
                deallocate (tau_Re_vf(cont_idx%end + i)%sf)
            end do
            deallocate (tau_Re_vf(E_idx)%sf)

            deallocate (tau_Re_vf)
        end if

        ! Deallocation of reg_src_vf
        if (regularization) then
            do i = 1, sys_size
                deallocate (reg_src_vf(i)%sf)
            end do
            deallocate (reg_src_vf)
        end if

        if (tvd_rhs_flux) then
            ! Deallocation of flux limiter variables
            deallocate (flux_lim_func)

            do k = ipsi%beg, ipsi%end
                do j = ichi%beg, ichi%end
                    do i = num_dims, 1, -1
                        if (abs(j) >= abs(k)) then
                            if (i /= 1) then
                                do l = 1, sys_size
                                    nullify (lo_flux_ndqp(i, j, k)%vf(l)%sf)
                                    nullify (hi_flux_ndqp(i, j, k)%vf(l)%sf)
                                    nullify (lo_flux_src_ndqp(i, j, k)%vf(l)%sf)
                                    nullify (hi_flux_src_ndqp(i, j, k)%vf(l)%sf)
                                    nullify (lo_flux_gsrc_ndqp(i, j, k)%vf(l)%sf)
                                    nullify (hi_flux_gsrc_ndqp(i, j, k)%vf(l)%sf)
                                end do
                            else
                                do l = 1, sys_size
                                    deallocate (lo_flux_ndqp(i, j, k)%vf(l)%sf)
                                    deallocate (hi_flux_ndqp(i, j, k)%vf(l)%sf)
                                    deallocate (lo_flux_gsrc_ndqp(i, j, k)%vf(l)%sf)
                                    deallocate (hi_flux_gsrc_ndqp(i, j, k)%vf(l)%sf)
                                end do

                                if (any(Re_size > 0) .or. We_size > 0) then
                                    do l = mom_idx%beg, E_idx
                                        deallocate (lo_flux_src_ndqp(i, j, k)%vf(l)%sf)
                                        deallocate (hi_flux_src_ndqp(i, j, k)%vf(l)%sf)
                                    end do
                                end if

                                if (riemann_solver == 1) then
                                    do l = adv_idx%beg + 1, adv_idx%end
                                        deallocate (lo_flux_src_ndqp(i, j, k)%vf(l)%sf)
                                        deallocate (hi_flux_src_ndqp(i, j, k)%vf(l)%sf)
                                    end do
                                else
                                    do l = adv_idx%beg + 1, adv_idx%end
                                        nullify (lo_flux_src_ndqp(i, j, k)%vf(l)%sf)
                                        nullify (hi_flux_src_ndqp(i, j, k)%vf(l)%sf)
                                    end do
                                end if

                                deallocate (lo_flux_src_ndqp(i, j, k)%vf(adv_idx%beg)%sf)
                                deallocate (hi_flux_src_ndqp(i, j, k)%vf(adv_idx%beg)%sf)
                            end if
                        end if

                        deallocate (lo_flux_ndqp(i, j, k)%vf, lo_flux_src_ndqp(i, j, k)%vf, lo_flux_gsrc_ndqp(i, j, k)%vf)
                        deallocate (hi_flux_ndqp(i, j, k)%vf, hi_flux_src_ndqp(i, j, k)%vf, hi_flux_gsrc_ndqp(i, j, k)%vf)
                    end do
                end do
            end do

            deallocate (lo_flux_ndqp, lo_flux_src_ndqp, lo_flux_gsrc_ndqp)
            deallocate (hi_flux_ndqp, hi_flux_src_ndqp, hi_flux_gsrc_ndqp)
        end if

        ! Deallocation/Disassociation of flux_ndqp, flux_src_ndqp, and flux_gsrc_ndqp ====
        do k = ipsi%beg, ipsi%end
            do j = ichi%beg, ichi%end
                do i = num_dims, 1, -1

                    if (abs(j) >= abs(k)) then
                        if (i /= 1) then

                            do l = 1, sys_size
                                nullify (flux_ndqp(i, j, k)%vf(l)%sf)
                                nullify (flux_src_ndqp(i, j, k)%vf(l)%sf)
                                nullify (flux_gsrc_ndqp(i, j, k)%vf(l)%sf)
                            end do

                        else

                            do l = 1, sys_size
                                deallocate (flux_ndqp(i, j, k)%vf(l)%sf)
                                deallocate (flux_gsrc_ndqp(i, j, k)%vf(l)%sf)
                            end do

                            if (any(Re_size > 0) .or. We_size > 0) then
                                do l = mom_idx%beg, E_idx
                                    deallocate (flux_src_ndqp(i, j, k)%vf(l)%sf)
                                end do
                            end if

                            if (riemann_solver == 1) then
                                do l = adv_idx%beg + 1, adv_idx%end
                                    deallocate (flux_src_ndqp(i, j, k)%vf(l)%sf)
                                end do
                            else
                                do l = adv_idx%beg + 1, adv_idx%end
                                    nullify (flux_src_ndqp(i, j, k)%vf(l)%sf)
                                end do
                            end if

                            deallocate (flux_src_ndqp(i, j, k)%vf(adv_idx%beg)%sf)

                        end if
                    end if

                    deallocate (flux_ndqp(i, j, k)%vf, flux_src_ndqp(i, j, k)%vf, flux_gsrc_ndqp(i, j, k)%vf)

                end do
            end do
        end do

        deallocate (flux_ndqp, flux_src_ndqp, flux_gsrc_ndqp)

        ! END: Deallocation/Disassociation of flux_ndqp, flux_src_ndqp, and flux_gsrc_ndqp  ===

        ! Disassociating procedural pointer to the subroutine which was
        ! utilized to calculate the solution of a given Riemann problem
        s_riemann_solver => null()

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_rhs_module ! ---------------------------------

end module m_rhs
