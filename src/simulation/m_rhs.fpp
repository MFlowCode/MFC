!>
!! @file m_rhs.f90
!! @brief Contains module m_rhs

#:include 'macros.fpp'

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

    use m_hypoelastic

    use m_monopole

    use m_viscous
    
    use m_nvtx
    
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_rhs_module, &
 s_compute_rhs, &
 s_pressure_relaxation_procedure, &
 s_finalize_rhs_module


    type(vector_field) :: q_cons_qp !<
    !! This variable contains the WENO-reconstructed values of the cell-average
    !! conservative variables, which are located in q_cons_vf, at cell-interior
    !! Gaussian quadrature points (QP).

    type(vector_field) :: q_prim_qp !<
    !! The primitive variables at cell-interior Gaussian quadrature points. These
    !! are calculated from the conservative variables and gradient magnitude (GM)
    !! of the volume fractions, q_cons_qp and gm_alpha_qp, respectively.

    !> @name The first-order spatial derivatives of the primitive variables at cell-
    !! interior Guassian quadrature points. These are WENO-reconstructed from
    !! their respective cell-average values, obtained through the application
    !! of the divergence theorem on the integral-average cell-boundary values
    !! of the primitive variables, located in qK_prim_n, where K = L or R.
    !> @{
    type(vector_field) :: dq_prim_dx_qp
    type(vector_field) :: dq_prim_dy_qp
    type(vector_field) :: dq_prim_dz_qp
    !> @}

    !> @name The left and right WENO-reconstructed cell-boundary values of the cell-
    !! average first-order spatial derivatives of the primitive variables. The
    !! cell-average of the first-order spatial derivatives may be found in the
    !! variables dq_prim_ds_qp, where s = x, y or z.
    !> @{
    type(vector_field), allocatable, dimension(:) :: dqL_prim_dx_n
    type(vector_field), allocatable, dimension(:) :: dqL_prim_dy_n
    type(vector_field), allocatable, dimension(:) :: dqL_prim_dz_n
    type(vector_field), allocatable, dimension(:) :: dqR_prim_dx_n
    type(vector_field), allocatable, dimension(:) :: dqR_prim_dy_n
    type(vector_field), allocatable, dimension(:) :: dqR_prim_dz_n
    !> @}

    type(vector_field) :: gm_alpha_qp  !<
    !! The gradient magnitude of the volume fractions at cell-interior Gaussian
    !! quadrature points. gm_alpha_qp is calculated from individual first-order
    !! spatial derivatives located in dq_prim_ds_qp.

    !> @name The left and right WENO-reconstructed cell-boundary values of the cell-
    !! average gradient magnitude of volume fractions, located in gm_alpha_qp.
    !> @{
    type(vector_field), allocatable, dimension(:) :: gm_alphaL_n
    type(vector_field), allocatable, dimension(:) :: gm_alphaR_n
    !> @}

    !> @name The cell-boundary values of the fluxes (src - source, gsrc - geometrical
    !! source). These are computed by applying the chosen Riemann problem solver
    !! .on the left and right cell-boundary values of the primitive variables
    !> @{
    type(vector_field), allocatable, dimension(:) :: flux_n
    type(vector_field), allocatable, dimension(:) :: flux_src_n
    type(vector_field), allocatable, dimension(:) :: flux_gsrc_n
    !> @}

    !> @name Additional field for capillary source terms
    !> @{
    type(scalar_field), allocatable, dimension(:) :: tau_Re_vf
    !> @}

    type(vector_field), allocatable, dimension(:) :: qL_prim, qR_prim

    type(int_bounds_info) :: iv !< Vector field indical bounds

    !> @name Indical bounds in the x-, y- and z-directions
    !> @{
    type(int_bounds_info) :: ix, iy, iz
    !> @}

    type(int_bounds_info) :: is1, is2, is3

    type(int_bounds_info) :: ixt, iyt, izt

    !> @name Bubble dynamic source terms
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :) :: bub_adv_src
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: bub_r_src, bub_v_src, bub_p_src, bub_m_src
    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: bub_mom_src

    type(scalar_field) :: divu !< matrix for div(u)
    !> @}

    !> @name Monopole source terms
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :) :: mono_mass_src, mono_e_src
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: mono_mom_src
    !> @}

    !> @name Saved fluxes for testing
    !> @{
    type(scalar_field) :: alf_sum
    !> @}

    real(kind(0d0)), allocatable, dimension(:, :, :) :: blkmod1, blkmod2, alpha1, alpha2, Kterm
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, qR_rsx_vf, qR_rsy_vf, qR_rsz_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf

    real(kind(0d0)), allocatable, dimension(:) :: gamma_min, pres_inf
    !$acc declare create(gamma_min, pres_inf)

    real(kind(0d0)), allocatable, dimension(:, :) :: Res
    !$acc declare create(Res)

!$acc declare create(q_cons_qp,q_prim_qp,  &
!$acc   dq_prim_dx_qp,dq_prim_dy_qp,dq_prim_dz_qp,dqL_prim_dx_n,dqL_prim_dy_n, &
!$acc   dqL_prim_dz_n,dqR_prim_dx_n,dqR_prim_dy_n,dqR_prim_dz_n,gm_alpha_qp,       &
!$acc   gm_alphaL_n,gm_alphaR_n,flux_n,flux_src_n,flux_gsrc_n,       &
!$acc   tau_Re_vf,qL_prim, qR_prim, iv,ix, iy, iz,is1,is2,is3,bub_adv_src,bub_r_src,bub_v_src, bub_p_src, bub_m_src, &
!$acc   bub_mom_src,alf_sum, &
!$acc   blkmod1, blkmod2, alpha1, alpha2, Kterm, divu, qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
!$acc   dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf, &
!$acc   ixt, iyt, izt)

    real(kind(0d0)), allocatable, dimension(:, :, :) :: nbub !< Bubble number density
!$acc declare create(nbub)

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_rhs_module() ! ---------------------------------

        integer :: i, j, k, l, d !< Generic loop iterators

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================

!$acc update device(ix, iy, iz)

        if (any(Re_size > 0) .and. cyl_coord) then
            @:ALLOCATE(tau_Re_vf(1:sys_size))
            do i = 1, num_dims
                @:ALLOCATE(tau_Re_vf(cont_idx%end + i)%sf(ix%beg:ix%end, &
                                                       &  iy%beg:iy%end, &
                                                       &  iz%beg:iz%end))
            end do
            @:ALLOCATE(tau_Re_vf(E_idx)%sf(ix%beg:ix%end, &
                                         & iy%beg:iy%end, &
                                         & iz%beg:iz%end))
        end if
        
        ixt = ix; iyt = iy; izt = iz

        @:ALLOCATE(q_cons_qp%vf(1:sys_size))
        @:ALLOCATE(q_prim_qp%vf(1:sys_size))

        do l = 1, sys_size
            @:ALLOCATE(q_cons_qp%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end do

        do l = mom_idx%beg, E_idx
            @:ALLOCATE(q_prim_qp%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end do

        do l = adv_idx%end + 1, sys_size
            @:ALLOCATE(q_prim_qp%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end do

        do l = 1, cont_idx%end
            q_prim_qp%vf(l)%sf => &
                q_cons_qp%vf(l)%sf
            !$acc enter data attach(q_prim_qp%vf(l)%sf)
        end do

        do l = adv_idx%beg, adv_idx%end
            q_prim_qp%vf(l)%sf => &
                q_cons_qp%vf(l)%sf
            !$acc enter data attach(q_prim_qp%vf(l)%sf)
        end do

        ! ==================================================================

        if (qbmm) then
            @:ALLOCATE(mom_sp(1:nmomsp), mom_3d(0:2, 0:2, nb))

            do i = 0, 2
                do j = 0, 2
                    do k = 1, nb
                        @:ALLOCATE(mom_3d(i, j, k)%sf( &
                                      & ix%beg:ix%end, &
                                      & iy%beg:iy%end, &
                                      & iz%beg:iz%end))
                    end do
                end do
            end do
            do i = 1, nmomsp
                @:ALLOCATE(mom_sp(i)%sf( &
                        & ix%beg:ix%end, &
                        & iy%beg:iy%end, &
                        & iz%beg:iz%end))
            end do
        end if

        ! Allocation/Association of qK_cons_n and qK_prim_n ==========
        @:ALLOCATE(qL_prim(1:num_dims))
        @:ALLOCATE(qR_prim(1:num_dims))

        do i = 1, num_dims
            @:ALLOCATE(qL_prim(i)%vf(1:sys_size))
            @:ALLOCATE(qR_prim(i)%vf(1:sys_size))
        end do

        if (weno_Re_flux) then

            do i = 1, num_dims
                do l = mom_idx%beg, mom_idx%end
                    @:ALLOCATE(qL_prim(i)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
                    @:ALLOCATE(qR_prim(i)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
                end do
            end do
        end if

        if (mpp_lim .and. bubbles) then
            @:ALLOCATE(alf_sum%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end if
        ! END: Allocation/Association of qK_cons_n and qK_prim_n =====

        @:ALLOCATE(qL_rsx_vf(ix%beg:ix%end, &
                                 iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
        @:ALLOCATE(qR_rsx_vf(ix%beg:ix%end, &
                                 iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))

        if (n > 0) then

            @:ALLOCATE(qL_rsy_vf(iy%beg:iy%end, &
                                     ix%beg:ix%end, iz%beg:iz%end, 1:sys_size))
            @:ALLOCATE(qR_rsy_vf(iy%beg:iy%end, &
                                     ix%beg:ix%end, iz%beg:iz%end, 1:sys_size))
        else
            @:ALLOCATE(qL_rsy_vf(ix%beg:ix%end, &
                                     iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
            @:ALLOCATE(qR_rsy_vf(ix%beg:ix%end, &
                                     iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
        end if

        if (p > 0) then
            @:ALLOCATE(qL_rsz_vf(iz%beg:iz%end, &
                                     iy%beg:iy%end, ix%beg:ix%end, 1:sys_size))
            @:ALLOCATE(qR_rsz_vf(iz%beg:iz%end, &
                                     iy%beg:iy%end, ix%beg:ix%end, 1:sys_size))
        else
            @:ALLOCATE(qL_rsz_vf(ix%beg:ix%end, &
                                     iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
            @:ALLOCATE(qR_rsz_vf(ix%beg:ix%end, &
                                     iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))

        end if
        ! Allocation of dq_prim_ds_qp ======================================

        if (any(Re_size > 0)) then

            @:ALLOCATE(dq_prim_dx_qp%vf(1:sys_size))
            @:ALLOCATE(dq_prim_dy_qp%vf(1:sys_size))
            @:ALLOCATE(dq_prim_dz_qp%vf(1:sys_size))
            
            if (any(Re_size > 0)) then

                do l = mom_idx%beg, mom_idx%end
                    @:ALLOCATE(dq_prim_dx_qp%vf(l)%sf( &
                              & ix%beg:ix%end, &
                              & iy%beg:iy%end, &
                              & iz%beg:iz%end))
                end do

                if (n > 0) then

                    do l = mom_idx%beg, mom_idx%end
                        @:ALLOCATE(dq_prim_dy_qp%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                    end do

                    if (p > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:ALLOCATE(dq_prim_dz_qp%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                        end do
                    end if

                end if

            end if

        end if
        ! END: Allocation of dq_prim_ds_qp =================================

        ! Allocation/Association of dqK_prim_ds_n =======================
        @:ALLOCATE(dqL_prim_dx_n(1:num_dims))
        @:ALLOCATE(dqL_prim_dy_n(1:num_dims))
        @:ALLOCATE(dqL_prim_dz_n(1:num_dims))
        @:ALLOCATE(dqR_prim_dx_n(1:num_dims))
        @:ALLOCATE(dqR_prim_dy_n(1:num_dims))
        @:ALLOCATE(dqR_prim_dz_n(1:num_dims))

        if (any(Re_size > 0)) then
            do i = 1, num_dims
                @:ALLOCATE(dqL_prim_dx_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqL_prim_dy_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqL_prim_dz_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqR_prim_dx_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqR_prim_dy_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqR_prim_dz_n(i)%vf(1:sys_size))
                
                if (any(Re_size > 0)) then

                    do l = mom_idx%beg, mom_idx%end
                        @:ALLOCATE(dqL_prim_dx_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                        @:ALLOCATE(dqR_prim_dx_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                    end do

                    if (n > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:ALLOCATE(dqL_prim_dy_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                            @:ALLOCATE(dqR_prim_dy_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                        end do
                    end if

                    if (p > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:ALLOCATE(dqL_prim_dz_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                            @:ALLOCATE(dqR_prim_dz_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                        end do
                    end if

                end if

            end do
        end if
        ! END: Allocation/Association of d K_prim_ds_n ==================

        if (any(Re_size > 0)) then
            if (weno_Re_flux) then
                @:ALLOCATE(dqL_rsx_vf(ix%beg:ix%end, &
                                          iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                @:ALLOCATE(dqR_rsx_vf(ix%beg:ix%end, &
                                          iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))

                if (n > 0) then

                    @:ALLOCATE(dqL_rsy_vf(iy%beg:iy%end, &
                                              ix%beg:ix%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE(dqR_rsy_vf(iy%beg:iy%end, &
                                              ix%beg:ix%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                else
                    @:ALLOCATE(dqL_rsy_vf(ix%beg:ix%end, &
                                              iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE(dqR_rsy_vf(ix%beg:ix%end, &
                                              iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))

                end if

                if (p > 0) then
                    @:ALLOCATE(dqL_rsz_vf(iz%beg:iz%end, &
                                              iy%beg:iy%end, ix%beg:ix%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE(dqR_rsz_vf(iz%beg:iz%end, &
                                              iy%beg:iy%end, ix%beg:ix%end, mom_idx%beg:mom_idx%end))
                else
                    @:ALLOCATE(dqL_rsz_vf(ix%beg:ix%end, &
                                              iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE(dqR_rsz_vf(ix%beg:ix%end, &
                                              iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))

                end if
            end if
        end if
        ! ==================================================================

        ! Allocation of gm_alphaK_n =====================================
        @:ALLOCATE(gm_alphaL_n(1:num_dims))
        @:ALLOCATE(gm_alphaR_n(1:num_dims))
        ! ==================================================================

        if (bubbles) then
            @:ALLOCATE(bub_adv_src(0:m, 0:n, 0:p))
            if (qbmm) then
                @:ALLOCATE(bub_mom_src(1:nmom, 0:m, 0:n, 0:p, 1:nb))
            else
                @:ALLOCATE(bub_r_src(0:m, 0:n, 0:p, 1:nb))
                @:ALLOCATE(bub_v_src(0:m, 0:n, 0:p, 1:nb))
                @:ALLOCATE(bub_p_src(0:m, 0:n, 0:p, 1:nb))
                @:ALLOCATE(bub_m_src(0:m, 0:n, 0:p, 1:nb))
            end if
        end if

        if (monopole) then
           @:ALLOCATE(mono_mass_src(0:m, 0:n, 0:p))
           @:ALLOCATE(mono_mom_src(1:num_dims, 0:m, 0:n, 0:p))
           @:ALLOCATE(mono_E_src(0:m, 0:n, 0:p))
        end if

        @:ALLOCATE(divu%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))

        ! Allocation/Association of flux_n, flux_src_n, and flux_gsrc_n ===
        @:ALLOCATE(flux_n(1:num_dims))
        @:ALLOCATE(flux_src_n(1:num_dims))
        @:ALLOCATE(flux_gsrc_n(1:num_dims))

        do i = 1, num_dims

            @:ALLOCATE(flux_n(i)%vf(1:sys_size))
            @:ALLOCATE(flux_src_n(i)%vf(1:sys_size))
            @:ALLOCATE(flux_gsrc_n(i)%vf(1:sys_size))

            if (i == 1) then

                do l = 1, sys_size
                    @:ALLOCATE(flux_n(i)%vf(l)%sf( &
                             & ix%beg:ix%end, &
                             & iy%beg:iy%end, &
                             & iz%beg:iz%end))
                    @:ALLOCATE(flux_gsrc_n(i)%vf(l)%sf( &
                            & ix%beg:ix%end, &
                            & iy%beg:iy%end, &
                            & iz%beg:iz%end))
                end do

                if (any(Re_size > 0)) then
                    do l = mom_idx%beg, E_idx
                        @:ALLOCATE(flux_src_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                    end do
                end if

                @:ALLOCATE(flux_src_n(i)%vf(adv_idx%beg)%sf( &
                         & ix%beg:ix%end, &
                         & iy%beg:iy%end, &
                         & iz%beg:iz%end))

                if (riemann_solver == 1) then
                    do l = adv_idx%beg + 1, adv_idx%end
                        @:ALLOCATE(flux_src_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                    end do
                else
                    do l = adv_idx%beg + 1, adv_idx%end
                        flux_src_n(i)%vf(l)%sf => &
                            flux_src_n(i)%vf(adv_idx%beg)%sf
                        !$acc enter data attach(flux_src_n(i)%vf(l)%sf(ix%beg:ix%end,iy%beg:iy%end,iz%beg:iz%end))
                    end do
                end if

            else
                do l = 1, sys_size
                    @:ALLOCATE(flux_gsrc_n(i)%vf(l)%sf( &
                              ix%beg:ix%end, &
                              iy%beg:iy%end, &
                              iz%beg:iz%end))
                end do
                do l = 1, sys_size
                    flux_n(i)%vf(l)%sf => &
                        flux_n(1)%vf(l)%sf
                    flux_src_n(i)%vf(l)%sf => &
                        flux_src_n(1)%vf(l)%sf

                    !$acc enter data attach(flux_n(i)%vf(l)%sf,flux_src_n(i)%vf(l)%sf)
                end do

            end if
        end do
        ! END: Allocation/Association of flux_n, flux_src_n, and flux_gsrc_n ===

        if (alt_soundspeed) then
            @:ALLOCATE(blkmod1(0:m, 0:n, 0:p), blkmod2(0:m, 0:n, 0:p), alpha1(0:m, 0:n, 0:p), alpha2(0:m, 0:n, 0:p), Kterm(0:m, 0:n, 0:p))
        end if

        @:ALLOCATE(gamma_min(1:num_fluids), pres_inf(1:num_fluids))

        do i = 1, num_fluids
            gamma_min(i) = 1d0/fluid_pp(i)%gamma + 1d0
            pres_inf(i) = fluid_pp(i)%pi_inf/(1d0 + fluid_pp(i)%gamma)
        end do
!$acc update device(gamma_min, pres_inf)

        if (any(Re_size > 0)) then
            @:ALLOCATE(Res(1:2, 1:maxval(Re_size)))
        end if

        if (any(Re_size > 0)) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
!$acc update device(Res, Re_idx, Re_size)
        end if

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem
        if (riemann_solver == 1) then
            s_riemann_solver => s_hll_riemann_solver
        elseif (riemann_solver == 2) then
            s_riemann_solver => s_hllc_riemann_solver
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


!$acc parallel loop collapse(5) gang vector default(present)
        do i = 1, sys_size
            do l = startz, p - startz
                do k = starty, n - starty
                    do j = startx, m - startx
                        do d = 1, num_dims
                            flux_gsrc_n(d)%vf(i)%sf(j, k, l) = 0d0;
                        end do
                    end do
                end do
            end do
        end do

        if (bubbles) then
            @:ALLOCATE(nbub(0:m, 0:n, 0:p))
        end if

    end subroutine s_initialize_rhs_module ! -------------------------------

    subroutine s_compute_rhs(q_cons_vf, q_prim_vf, rhs_vf, t_step) ! -------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf
        integer, intent(IN) :: t_step
        
        real(kind(0d0)) :: top, bottom  !< Numerator and denominator when evaluating flux limiter function
        real(kind(0d0)), dimension(num_fluids) :: myalpha_rho, myalpha

        real(kind(0d0)) :: tmp1, tmp2, tmp3, tmp4, &
                           c_gas, c_liquid, &
                           Cpbw, Cpinf, Cpinf_dot, &
                           myH, myHdot, rddot, alf_gas

        real(kind(0d0)) :: pb, mv, vflux, pldot, pbdot

        real(kind(0d0)) :: n_tait, B_tait, angle, angle_z

        real(kind(0d0)), dimension(nb) :: Rtmp, Vtmp
        real(kind(0d0)) :: myR, myV, alf, myP, myRho, R2Vav
        integer :: ndirs

        real(kind(0d0)) :: mytime, sound
        real(kind(0d0)) :: start, finish
        real(kind(0d0)) :: s2, const_sos, s1

        integer :: i, j, k, l, r, q, ii, id !< Generic loop iterators
        integer :: term_index

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================

        !$acc update device(ix, iy, iz)

        ! Association/Population of Working Variables ======================
        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        q_cons_qp%vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        call nvtxStartRange("RHS-MPI")
        call s_populate_conservative_variables_buffers()
        call nvtxEndRange

        
        ! ==================================================================

        ! Converting Conservative to Primitive Variables ==================

        if (mpp_lim .and. bubbles) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        alf_sum%sf(j, k, l) = 0d0
                        !$acc loop seq
                        do i = advxb, advxe - 1
                            alf_sum%sf(j, k, l) = alf_sum%sf(j, k, l) + q_cons_qp%vf(i)%sf(j, k, l)
                        end do
                        !$acc loop seq
                        do i = advxb, advxe - 1
                            q_cons_qp%vf(i)%sf(j, k, l) = q_cons_qp%vf(i)%sf(j, k, l)*(1.d0 - q_cons_qp%vf(alf_idx)%sf(j, k, l)) &
                                                          /alf_sum%sf(j, k, l)
                        end do
                    end do
                end do
            end do

        end if

        call nvtxStartRange("RHS-CONVERT")
        call s_convert_conservative_to_primitive_variables( &
            q_cons_qp%vf, &
            q_prim_qp%vf, &
            gm_alpha_qp%vf, &
            ix, iy, iz)
        call nvtxEndRange


        
        if (t_step == t_step_stop) return
        ! ==================================================================

        if (qbmm) call s_mom_inv(q_prim_qp%vf, mom_sp, mom_3d, ix, iy, iz)

        call nvtxStartRange("Viscous")
        if (any(Re_size > 0)) call s_get_viscous(qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                                            dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, &
                                            qL_prim, &
                                            qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                                            dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n, &
                                            qR_prim, &
                                            q_prim_qp, &
                                            dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp, &
                                            ix, iy, iz)
        call nvtxEndRange()
        
        ! Dimensional Splitting Loop =======================================

        do id = 1, num_dims

            ! Configuring Coordinate Direction Indexes ======================
            ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

            if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

            ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
            ! ===============================================================
            ! Reconstructing Primitive/Conservative Variables ===============
            
            if (all(Re_size == 0)) then
                    iv%beg = 1; iv%end = sys_size
                !call nvtxStartRange("RHS-WENO")
                call nvtxStartRange("RHS-WENO")
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(1:sys_size), &
                    qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                    qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                    id)
                call nvtxEndRange
            else
                call nvtxStartRange("RHS-WENO")
                iv%beg = 1; iv%end = contxe
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(iv%beg:iv%end), &
                    qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                    qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                    id)

                iv%beg = E_idx; iv%end = E_idx
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(iv%beg:iv%end), &
                    qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                    qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                    id)

                iv%beg = advxb; iv%end = advxe
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(iv%beg:iv%end), &
                    qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                    qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                    id)

                iv%beg = mom_idx%beg; iv%end = mom_idx%end
                if (weno_Re_flux) then
                    call s_reconstruct_cell_boundary_values_visc_deriv( &
                        dq_prim_dx_qp%vf(iv%beg:iv%end), &
                        dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, &
                        dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf, &
                        id, dqL_prim_dx_n(id)%vf(iv%beg:iv%end), dqR_prim_dx_n(id)%vf(iv%beg:iv%end), &
                        ix, iy, iz)
                    if (n > 0) then
                        call s_reconstruct_cell_boundary_values_visc_deriv( &
                            dq_prim_dy_qp%vf(iv%beg:iv%end), &
                            dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, &
                            dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf, &
                            id, dqL_prim_dy_n(id)%vf(iv%beg:iv%end), dqR_prim_dy_n(id)%vf(iv%beg:iv%end), &
                            ix, iy, iz)
                        if (p > 0) then
                            call s_reconstruct_cell_boundary_values_visc_deriv( &
                                dq_prim_dz_qp%vf(iv%beg:iv%end), &
                                dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, &
                                dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf, &
                                id, dqL_prim_dz_n(id)%vf(iv%beg:iv%end), dqR_prim_dz_n(id)%vf(iv%beg:iv%end), &
                                ix, iy, iz)
                        end if
                    end if
                end if
                call nvtxEndRange
            end if

            ! Configuring Coordinate Direction Indexes ======================
            if (id == 1) then
                ix%beg = -1; iy%beg = 0; iz%beg = 0
            elseif (id == 2) then
                ix%beg = 0; iy%beg = -1; iz%beg = 0
            else
                ix%beg = 0; iy%beg = 0; iz%beg = -1
            end if
            ix%end = m; iy%end = n; iz%end = p
            ! ===============================================================
            call nvtxStartRange("RHS-Riemann")

            ! Computing Riemann Solver Flux and Source Flux =================

            call s_riemann_solver(qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                                  dqR_prim_dx_n(id)%vf, &
                                  dqR_prim_dy_n(id)%vf, &
                                  dqR_prim_dz_n(id)%vf, &
                                  qR_prim(id)%vf, &
                                  qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                                  dqL_prim_dx_n(id)%vf, &
                                  dqL_prim_dy_n(id)%vf, &
                                  dqL_prim_dz_n(id)%vf, &
                                  qL_prim(id)%vf, &
                                  q_prim_qp%vf, &
                                  flux_n(id)%vf, &
                                  flux_src_n(id)%vf, &
                                  flux_gsrc_n(id)%vf, &
                                  id, ix, iy, iz)
            call nvtxEndRange

            ! ===============================================================


            if (alt_soundspeed) then
!$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            blkmod1(j, k, l) = ((gammas(1) + 1d0)*q_prim_qp%vf(E_idx)%sf(j, k, l) + &
                                                pi_infs(1))/gammas(1)
                            blkmod2(j, k, l) = ((gammas(2) + 1d0)*q_prim_qp%vf(E_idx)%sf(j, k, l) + &
                                                pi_infs(2))/gammas(2)
                            alpha1(j, k, l) = q_cons_qp%vf(advxb)%sf(j, k, l)

                            if (bubbles) then
                                alpha2(j, k, l) = q_cons_qp%vf(alf_idx - 1)%sf(j, k, l)
                            else
                                alpha2(j, k, l) = q_cons_qp%vf(advxe)%sf(j, k, l)
                            end if

                            Kterm(j, k, l) = alpha1(j, k, l)*alpha2(j, k, l)*(blkmod2(j, k, l) - blkmod1(j, k, l))/ &
                                             (alpha1(j, k, l)*blkmod2(j, k, l) + alpha2(j, k, l)*blkmod1(j, k, l))
                        end do
                    end do
                end do
            end if

            call nvtxStartRange("RHS_Flux_Add")
            if (id == 1) then

                if (bc_x%beg <= -5) then
                    call s_cbc(q_prim_qp%vf, flux_n(id)%vf, &
                               flux_src_n(id)%vf, id, -1, ix, iy, iz)
                end if

                if (bc_x%end <= -5) then
                    call s_cbc(q_prim_qp%vf, flux_n(id)%vf, &
                               flux_src_n(id)%vf, id, 1, ix, iy, iz)
                end if

                !$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, sys_size
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                rhs_vf(j)%sf(k, l, q) = 1d0/dx(k)* &
                                                        (flux_n(1)%vf(j)%sf(k - 1, l, q) &
                                                         - flux_n(1)%vf(j)%sf(k, l, q))
                            end do
                        end do
                    end do
                end do

                if (riemann_solver == 1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do j = advxb, advxe
                        do q = 0, p
                            do l = 0, n
                                do k = 0, m
                                    rhs_vf(j)%sf(k, l, q) = &
                                        rhs_vf(j)%sf(k, l, q) + 1d0/dx(k)* &
                                        q_prim_qp%vf(contxe + id)%sf(k, l, q)* &
                                        (flux_src_n(1)%vf(j)%sf(k - 1, l, q) &
                                         - flux_src_n(1)%vf(j)%sf(k, l, q))
                                end do
                            end do
                        end do
                    end do
                else
                    if (alt_soundspeed) then
                        do j = advxb, advxe
                            if ((j == advxe) .and. (bubbles .neqv. .true.)) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do q = 0, p
                                    do l = 0, n
                                        do k = 0, m
                                            rhs_vf(j)%sf(k, l, q) = &
                                                rhs_vf(j)%sf(k, l, q) + 1d0/dx(k)* &
                                                (q_cons_qp%vf(j)%sf(k, l, q) - Kterm(k, l, q))* &
                                                (flux_src_n(1)%vf(j)%sf(k, l, q) &
                                                 - flux_src_n(1)%vf(j)%sf(k - 1, l, q))
                                        end do
                                    end do
                                end do
                            else if ((j == advxb) .and. (bubbles .neqv. .true.)) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do q = 0, p
                                    do l = 0, n
                                        do k = 0, m
                                            rhs_vf(j)%sf(k, l, q) = &
                                                rhs_vf(j)%sf(k, l, q) + 1d0/dx(k)* &
                                                (q_cons_qp%vf(j)%sf(k, l, q) + Kterm(k, l, q))* &
                                                (flux_src_n(1)%vf(j)%sf(k, l, q) &
                                                 - flux_src_n(1)%vf(j)%sf(k - 1, l, q))
                                        end do
                                    end do
                                end do
                            end if
                        end do
                    else
                        !$acc parallel loop collapse(4) gang vector default(present)
                        do j = advxb, advxe
                            do q = 0, p
                                do l = 0, n
                                    do k = 0, m
                                        rhs_vf(j)%sf(k, l, q) = &
                                            rhs_vf(j)%sf(k, l, q) + 1d0/dx(k)* &
                                            q_cons_qp%vf(j)%sf(k, l, q)* &
                                            (flux_src_n(1)%vf(j)%sf(k, l, q) &
                                             - flux_src_n(1)%vf(j)%sf(k - 1, l, q))
                                    end do
                                end do
                            end do
                        end do
                    end if
                end if

                if (bubbles) then
                    if (qbmm) then

                        !$acc parallel loop collapse(3) gang vector default(present)
                        do l = 0, p
                            do q = 0, n
                                do i = 0, m

                                    rhs_vf(alf_idx)%sf(i, q, l) = rhs_vf(alf_idx)%sf(i, q, l) + mom_sp(2)%sf(i, q, l)
                                    j = bubxb
                                    !$acc loop seq
                                    do k = 1, nb
                                        rhs_vf(j)%sf(i, q, l) = &
                                            rhs_vf(j)%sf(i, q, l) + mom_3d(0, 0, k)%sf(i, q, l)
                                        rhs_vf(j + 1)%sf(i, q, l) = &
                                            rhs_vf(j + 1)%sf(i, q, l) + mom_3d(1, 0, k)%sf(i, q, l)
                                        rhs_vf(j + 2)%sf(i, q, l) = &
                                            rhs_vf(j + 2)%sf(i, q, l) + mom_3d(0, 1, k)%sf(i, q, l)
                                        rhs_vf(j + 3)%sf(i, q, l) = &
                                            rhs_vf(j + 3)%sf(i, q, l) + mom_3d(2, 0, k)%sf(i, q, l)
                                        rhs_vf(j + 4)%sf(i, q, l) = &
                                            rhs_vf(j + 4)%sf(i, q, l) + mom_3d(1, 1, k)%sf(i, q, l)
                                        rhs_vf(j + 5)%sf(i, q, l) = &
                                            rhs_vf(j + 5)%sf(i, q, l) + mom_3d(0, 2, k)%sf(i, q, l)
                                        j = j + 6
                                    end do
                                end do
                            end do
                        end do
                    else
!$acc parallel loop collapse(3) gang vector default(present)
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    divu%sf(j, k, l) = 0d0
                                    divu%sf(j, k, l) = &
                                        5d-1/dx(j)*(q_prim_qp%vf(contxe + id)%sf(j + 1, k, l) - &
                                                    q_prim_qp%vf(contxe + id)%sf(j - 1, k, l))

                                end do
                            end do
                        end do

                        ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3
                        if (id == ndirs) then                        
                            call s_compute_bubble_source(bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src, divu, nbub, &
                                                 q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), t_step, id, rhs_vf)
                        end if
                    end if
                end if    

                if (monopole) then
                    ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3
                    if (id == ndirs) then 
                        call s_monopole_calculations(mono_mass_src, mono_mom_src, mono_e_src, &
                                             q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), t_step, id, &
                                             rhs_vf)
                    end if
                end if

                if (model_eqns == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                        rhs_vf(i + intxb - 1)%sf(j, k, l) - 1d0/dx(j)* &
                                        q_cons_qp%vf(i + advxb - 1)%sf(j, k, l)* &
                                        q_prim_qp%vf(E_idx)%sf(j, k, l)* &
                                        (flux_src_n(1)%vf(advxb)%sf(j, k, l) - &
                                         flux_src_n(1)%vf(advxb)%sf(j - 1, k, l))
                                end do
                            end do
                        end do
                    end do
                end if

                if (any(Re_size > 0)) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                !$acc loop seq
                                do i = momxb, E_idx
                                    rhs_vf(i)%sf(j, k, l) = &
                                        rhs_vf(i)%sf(j, k, l) + 1d0/dx(j)* &
                                        (flux_src_n(1)%vf(i)%sf(j - 1, k, l) &
                                         - flux_src_n(1)%vf(i)%sf(j, k, l))
                                end do
                            end do
                        end do
                    end do
                end if

            elseif (id == 2) then
                ! RHS Contribution in y-direction ===============================
                ! Applying the Riemann fluxes

                if (bc_y%beg <= -5 .and. bc_y%beg /= -13) then
                    call s_cbc(q_prim_qp%vf, flux_n(id)%vf, &
                               flux_src_n(id)%vf, id, -1, ix, iy, iz)
                end if

                if (bc_y%end <= -5) then
                    call s_cbc(q_prim_qp%vf, flux_n(id)%vf, &
                               flux_src_n(id)%vf, id, 1, ix, iy, iz)
                end if

                !$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do q = 0, m
                                rhs_vf(j)%sf(q, k, l) = &
                                    rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                    (flux_n(2)%vf(j)%sf(q, k - 1, l) &
                                     - flux_n(2)%vf(j)%sf(q, k, l))
                            end do
                        end do
                    end do
                end do
                ! Applying source terms to the RHS of the advection equations

                if (riemann_solver == 1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do j = advxb, advxe
                        do l = 0, p
                            do k = 0, n
                                do q = 0, m
                                    rhs_vf(j)%sf(q, k, l) = &
                                        rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                        q_prim_qp%vf(contxe + id)%sf(q, k, l)* &
                                        (flux_src_n(2)%vf(j)%sf(q, k - 1, l) &
                                         - flux_src_n(2)%vf(j)%sf(q, k, l))
                                end do
                            end do
                        end do
                    end do
                else

                    if (alt_soundspeed) then
                        do j = advxb, advxe
                            if ((j == advxe) .and. (bubbles .neqv. .true.)) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do l = 0, p
                                    do k = 0, n
                                        do q = 0, m
                                            rhs_vf(j)%sf(q, k, l) = &
                                                rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                                (q_cons_qp%vf(j)%sf(q, k, l) - Kterm(q, k, l))* &
                                                (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                 - flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                        end do
                                    end do
                                end do
                                if (cyl_coord) then
                                    !$acc parallel loop collapse(3) gang vector default(present)
                                    do l = 0, p
                                        do k = 0, n
                                            do q = 0, m
                                                rhs_vf(j)%sf(q, k, l) = &
                                                    rhs_vf(j)%sf(q, k, l) - &
                                                    (Kterm(q, k, l)/2d0/y_cc(k))* &
                                                    (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                     + flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                            end do
                                        end do
                                    end do
                                end if
                            else if ((j == advxb) .and. (bubbles .neqv. .true.)) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do l = 0, p
                                    do k = 0, n
                                        do q = 0, m
                                            rhs_vf(j)%sf(q, k, l) = &
                                                rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                                (q_cons_qp%vf(j)%sf(q, k, l) + Kterm(q, k, l))* &
                                                (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                 - flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                        end do
                                    end do
                                end do
                                if (cyl_coord) then
                                    !$acc parallel loop collapse(3) gang vector default(present)
                                    do l = 0, p
                                        do k = 0, n
                                            do q = 0, m
                                                rhs_vf(j)%sf(q, k, l) = &
                                                    rhs_vf(j)%sf(q, k, l) + &
                                                    (Kterm(q, k, l)/2d0/y_cc(k))* &
                                                    (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                     + flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                            end do
                                        end do
                                    end do
                                end if
                            end if
                        end do
                    else
!$acc parallel loop collapse(4) gang vector default(present)
                        do j = advxb, advxe
                            do l = 0, p
                                do k = 0, n
                                    do q = 0, m
                                        rhs_vf(j)%sf(q, k, l) = &
                                            rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                            q_cons_qp%vf(j)%sf(q, k, l)* &
                                            (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                             - flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                    end do
                                end do
                            end do
                        end do
                    end if
                end if

                if (bubbles .and. (.not. qbmm)) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                divu%sf(j, k, l) = divu%sf(j, k, l) + &
                                                   5d-1/dy(k)*(q_prim_qp%vf(contxe + id)%sf(j, k + 1, l) - &
                                                               q_prim_qp%vf(contxe + id)%sf(j, k - 1, l))

                            end do
                        end do
                    end do

                    ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3
                    if (id == ndirs) then                        
                        call s_compute_bubble_source(bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src, divu, nbub, &
                                             q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), t_step, id, rhs_vf)
                    end if

                end if



                if (monopole) then
                    ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3
                    if (id == ndirs) then 
                        call s_monopole_calculations(mono_mass_src, mono_mom_src, mono_e_src, &
                                             q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), t_step, id, &
                                             rhs_vf)
                    end if
                end if

                if (model_eqns == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                        rhs_vf(i + intxb - 1)%sf(j, k, l) - 1d0/dy(k)* &
                                        q_cons_qp%vf(i + advxb - 1)%sf(j, k, l)* &
                                        q_prim_qp%vf(E_idx)%sf(j, k, l)* &
                                        (flux_src_n(2)%vf(advxb)%sf(j, k, l) - &
                                         flux_src_n(2)%vf(advxb)%sf(j, k - 1, l))
                                end do
                            end do
                        end do
                    end do

                    if (cyl_coord) then
                        !$acc parallel loop collapse(4) gang vector default(present)
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    do i = 1, num_fluids
                                        rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                            rhs_vf(i + intxb - 1)%sf(j, k, l) - 5d-1/y_cc(k)* &
                                            q_cons_qp%vf(i + advxb - 1)%sf(j, k, l)* &
                                            q_prim_qp%vf(E_idx)%sf(j, k, l)* &
                                            (flux_src_n(2)%vf(advxb)%sf(j, k, l) + &
                                             flux_src_n(2)%vf(advxb)%sf(j, k - 1, l))
                                    end do
                                end do
                            end do
                        end do
                    end if
                end if

                if (cyl_coord) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do j = 1, sys_size
                        do l = 0, p
                            do k = 0, n
                                do q = 0, m
                                    rhs_vf(j)%sf(q, k, l) = &
                                        rhs_vf(j)%sf(q, k, l) - 5d-1/y_cc(k)* &
                                        (flux_gsrc_n(2)%vf(j)%sf(q, k - 1, l) &
                                         + flux_gsrc_n(2)%vf(j)%sf(q, k, l))
                                end do
                            end do
                        end do
                    end do
                end if

                if (any(Re_size > 0)) then
                    if (cyl_coord .and. ((bc_y%beg == -2) .or. (bc_y%beg == -13))) then
                        if (p > 0) then
                            call s_compute_viscous_stress_tensor(q_prim_qp%vf, &
                                                                 dq_prim_dx_qp%vf(mom_idx%beg:mom_idx%end), &
                                                                 dq_prim_dy_qp%vf(mom_idx%beg:mom_idx%end), &
                                                                 dq_prim_dz_qp%vf(mom_idx%beg:mom_idx%end), &
                                                                 tau_Re_vf, &
                                                                 ixt, iyt, izt)
                        else
                            call s_compute_viscous_stress_tensor(q_prim_qp%vf, &
                                                                 dq_prim_dx_qp%vf(mom_idx%beg:mom_idx%end), &
                                                                 dq_prim_dy_qp%vf(mom_idx%beg:mom_idx%end), &
                                                                 dq_prim_dy_qp%vf(mom_idx%beg:mom_idx%end), &
                                                                 tau_Re_vf, &
                                                                 ixt, iyt, izt)
                        end if

                        !$acc parallel loop collapse(3) gang vector default(present)
                        do l = 0, p
                            do k = 1, n
                                do j = 0, m
                                    !$acc loop seq
                                    do i = momxb, E_idx
                                        rhs_vf(i)%sf(j, k, l) = &
                                            rhs_vf(i)%sf(j, k, l) + 1d0/dy(k)* &
                                            (flux_src_n(2)%vf(i)%sf(j, k - 1, l) &
                                             - flux_src_n(2)%vf(i)%sf(j, k, l))
                                    end do
                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(2) gang vector default(present)
                        do l = 0, p
                            do j = 0, m
                                !$acc loop seq
                                do i = momxb, E_idx
                                    rhs_vf(i)%sf(j, 0, l) = &
                                        rhs_vf(i)%sf(j, 0, l) + 1d0/(y_cc(1) - y_cc(-1))* &
                                        (tau_Re_vf(i)%sf(j, -1, l) &
                                         - tau_Re_vf(i)%sf(j, 1, l))
                                end do
                            end do
                        end do
                    else
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    !$acc loop seq
                                    do i = momxb, E_idx
                                        rhs_vf(i)%sf(j, k, l) = &
                                            rhs_vf(i)%sf(j, k, l) + 1d0/dy(k)* &
                                            (flux_src_n(2)%vf(i)%sf(j, k - 1, l) &
                                             - flux_src_n(2)%vf(i)%sf(j, k, l))
                                    end do
                                end do
                            end do
                        end do
                    end if
                    ! Applying the geometrical viscous Riemann source fluxes calculated as average
                    ! of values at cell boundaries
                    if (cyl_coord) then
                        if ((bc_y%beg == -2) .or. (bc_y%beg == -13)) then

                            !$acc parallel loop collapse(3) gang vector default(present)
                            do l = 0, p
                                do k = 1, n
                                    do j = 0, m
                                        !$acc loop seq
                                        do i = momxb, E_idx
                                            rhs_vf(i)%sf(j, k, l) = &
                                                rhs_vf(i)%sf(j, k, l) - 5d-1/y_cc(k)* &
                                                (flux_src_n(2)%vf(i)%sf(j, k - 1, l) &
                                                 + flux_src_n(2)%vf(i)%sf(j, k, l))
                                        end do
                                    end do
                                end do
                            end do

                            !$acc parallel loop collapse(2) gang vector default(present)
                            do l = 0, p
                                do j = 0, m
                                    !$acc loop seq
                                    do i = momxb, E_idx
                                        rhs_vf(i)%sf(j, 0, l) = &
                                            rhs_vf(i)%sf(j, 0, l) - 1d0/y_cc(0)* &
                                            tau_Re_vf(i)%sf(j, 0, l)
                                    end do
                                end do
                            end do

                        else

                            !$acc parallel loop collapse(3) gang vector default(present)
                            do l = 0, p
                                do k = 0, n
                                    do j = 0, m
                                        !$acc loop seq
                                        do i = momxb, E_idx
                                            rhs_vf(i)%sf(j, k, l) = &
                                                rhs_vf(i)%sf(j, k, l) - 5d-1/y_cc(k)* &
                                                (flux_src_n(2)%vf(i)%sf(j, k - 1, l) &
                                                 + flux_src_n(2)%vf(i)%sf(j, k, l))
                                        end do
                                    end do
                                end do
                            end do

                        end if
                    end if
                end if

            elseif (id == 3) then
                ! RHS Contribution in z-direction ===============================

                ! Applying the Riemann fluxes

                if (bc_z%beg <= -5) then
                    call s_cbc(q_prim_qp%vf, flux_n(id)%vf, &
                               flux_src_n(id)%vf, id, -1, ix, iy, iz)
                end if

                if (bc_z%end <= -5) then
                    call s_cbc(q_prim_qp%vf, flux_n(id)%vf, &
                               flux_src_n(id)%vf, id, 1, ix, iy, iz)
                end if

                if (grid_geometry == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do j = 1, sys_size
                        do k = 0, p
                            do q = 0, n
                                do l = 0, m
                                    rhs_vf(j)%sf(l, q, k) = &
                                        rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)/y_cc(q)* &
                                        q_prim_qp%vf(contxe + id)%sf(l, q, k)* &
                                        (flux_n(3)%vf(j)%sf(l, q, k - 1) &
                                         - flux_n(3)%vf(j)%sf(l, q, k))
                                end do
                            end do
                        end do
                    end do

                    if (riemann_solver == 1) then
                        do j = advxb, advxe
                            do k = 0, p
                                do q = 0, n
                                    do l = 0, m
                                        rhs_vf(j)%sf(l, q, k) = &
                                            rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)/y_cc(q)* &
                                            q_prim_qp%vf(contxe + id)%sf(l, q, k)* &
                                            (flux_src_n(3)%vf(j)%sf(l, q, k - 1) &
                                             - flux_src_n(3)%vf(j)%sf(l, q, k))
                                    end do
                                end do
                            end do
                        end do
                    else

                        if (alt_soundspeed) then
                            do j = advxb, advxe
                                if ((j == advxe) .and. (bubbles .neqv. .true.)) then
                                    !$acc parallel loop collapse(3) gang vector default(present)
                                    do k = 0, p
                                        do q = 0, n
                                            do l = 0, m
                                                rhs_vf(j)%sf(l, q, k) = &
                                                    rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)/y_cc(q)* &
                                                    (q_cons_qp%vf(j)%sf(l, q, k) - Kterm(l, q, k))* &
                                                    (flux_src_n(3)%vf(j)%sf(l, q, k) &
                                                     - flux_src_n(3)%vf(j)%sf(l, q, k - 1))
                                            end do
                                        end do
                                    end do
                                else if ((j == advxb) .and. (bubbles .neqv. .true.)) then
                                    !$acc parallel loop collapse(3) gang vector default(present)
                                    do k = 0, p
                                        do q = 0, n
                                            do l = 0, m
                                                rhs_vf(j)%sf(l, q, k) = &
                                                    rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)/y_cc(q)* &
                                                    (q_cons_qp%vf(j)%sf(l, q, k) + Kterm(l, q, k))* &
                                                    (flux_src_n(3)%vf(j)%sf(l, q, k) &
                                                     - flux_src_n(3)%vf(j)%sf(l, q, k - 1))
                                            end do
                                        end do
                                    end do
                                end if
                            end do
                        else
                            !$acc parallel loop collapse(4) gang vector default(present)
                            do j = advxb, advxe
                                do k = 0, p
                                    do q = 0, n
                                        do l = 0, m
                                            rhs_vf(j)%sf(l, q, k) = &
                                                rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)/y_cc(q)* &
                                                q_cons_qp%vf(j)%sf(l, q, k)* &
                                                (flux_src_n(3)%vf(j)%sf(l, q, k) &
                                                 - flux_src_n(3)%vf(j)%sf(l, q, k - 1))
                                        end do
                                    end do
                                end do
                            end do
                        end if
                    end if

                    !$acc parallel loop collapse(4) gang vector default(present)
                    do j = 1, sys_size
                        do k = 0, p
                            do q = 0, n
                                do l = 0, m
                                    rhs_vf(j)%sf(l, q, k) = &
                                        rhs_vf(j)%sf(l, q, k) - 5d-1/y_cc(q)* &
                                        (flux_gsrc_n(3)%vf(j)%sf(l, q, k - 1) &
                                         - flux_gsrc_n(3)%vf(j)%sf(l, q, k))
                                end do
                            end do
                        end do
                    end do

                else
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do j = 1, sys_size
                        do k = 0, p
                            do q = 0, n
                                do l = 0, m
                                    rhs_vf(j)%sf(l, q, k) = &
                                        rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                        (flux_n(3)%vf(j)%sf(l, q, k - 1) &
                                         - flux_n(3)%vf(j)%sf(l, q, k))
                                end do
                            end do
                        end do
                    end do

                    if (riemann_solver == 1) then
                        !$acc parallel loop collapse(4) gang vector default(present)
                        do j = advxb, advxe
                            do k = 0, p
                                do q = 0, n
                                    do l = 0, m
                                        rhs_vf(j)%sf(l, q, k) = &
                                            rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                            q_prim_qp%vf(contxe + id)%sf(l, q, k)* &
                                            (flux_src_n(3)%vf(j)%sf(l, q, k - 1) &
                                             - flux_src_n(3)%vf(j)%sf(l, q, k))
                                    end do
                                end do
                            end do
                        end do
                    else

                        if (alt_soundspeed) then
                            do j = advxb, advxe
                                if ((j == advxe) .and. (bubbles .neqv. .true.)) then
                                    !$acc parallel loop collapse(3) gang vector default(present)
                                    do k = 0, p
                                        do q = 0, n
                                            do l = 0, m
                                                rhs_vf(j)%sf(l, q, k) = &
                                                    rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                                    (q_cons_qp%vf(j)%sf(l, q, k) - Kterm(l, q, k))* &
                                                    (flux_src_n(3)%vf(j)%sf(l, q, k) &
                                                     - flux_src_n(3)%vf(j)%sf(l, q, k - 1))
                                            end do
                                        end do
                                    end do
                                else if ((j == advxb) .and. (bubbles .neqv. .true.)) then
                                    !$acc parallel loop collapse(3) gang vector default(present)
                                    do k = 0, p
                                        do q = 0, n
                                            do l = 0, m
                                                rhs_vf(j)%sf(l, q, k) = &
                                                    rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                                    (q_cons_qp%vf(j)%sf(l, q, k) + Kterm(l, q, k))* &
                                                    (flux_src_n(3)%vf(j)%sf(l, q, k) &
                                                     - flux_src_n(3)%vf(j)%sf(l, q, k - 1))
                                            end do
                                        end do
                                    end do
                                end if
                            end do
                        else
                            !$acc parallel loop collapse(4) gang vector default(present)
                            do j = advxb, advxe
                                do k = 0, p
                                    do q = 0, n
                                        do l = 0, m
                                            rhs_vf(j)%sf(l, q, k) = &
                                                rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                                q_cons_qp%vf(j)%sf(l, q, k)* &
                                                (flux_src_n(3)%vf(j)%sf(l, q, k) &
                                                 - flux_src_n(3)%vf(j)%sf(l, q, k - 1))
                                        end do
                                    end do
                                end do
                            end do
                        end if
                    end if
                end if

                call nvtxStartRange("bubbles")
                if (bubbles .and. (.not. qbmm)) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                divu%sf(j, k, l) = divu%sf(j, k, l) + &
                                                   5d-1/dz(l)*(q_prim_qp%vf(contxe + id)%sf(j, k, l + 1) - &
                                                               q_prim_qp%vf(contxe + id)%sf(j, k, l - 1))

                            end do
                        end do
                    end do

                    ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3
                    if (id == ndirs) then                        
                        call s_compute_bubble_source(bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src, divu, nbub, &
                                             q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), t_step, id, rhs_vf)
                    end if

                end if
                call nvtxEndRange()

                call nvtxStartRange("Monopole")

                if (monopole) then
                    ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3
                    if (id == ndirs) then 
                        call s_monopole_calculations(mono_mass_src, mono_mom_src, mono_e_src, &
                                             q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), t_step, id, &
                                             rhs_vf)
                    end if
                end if

                call nvtxEndRange()

                if (model_eqns == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                        rhs_vf(i + intxb - 1)%sf(j, k, l) - 1d0/dz(l)* &
                                        q_cons_qp%vf(i + advxb - 1)%sf(j, k, l)* &
                                        q_prim_qp%vf(E_idx)%sf(j, k, l)* &
                                        (flux_src_n(3)%vf(advxb)%sf(j, k, l) - &
                                         flux_src_n(3)%vf(advxb)%sf(j, k, l - 1))
                                end do
                            end do
                        end do
                    end do
                end if

                if (any(Re_size > 0)) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                !$acc loop seq
                                do i = momxb, E_idx
                                    rhs_vf(i)%sf(j, k, l) = &
                                        rhs_vf(i)%sf(j, k, l) + 1d0/dz(l)* &
                                        (flux_src_n(3)%vf(i)%sf(j, k, l - 1) &
                                         - flux_src_n(3)%vf(i)%sf(j, k, l))
                                end do
                            end do
                        end do
                    end do

                    if (grid_geometry == 3) then
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    rhs_vf(momxb + 1)%sf(j, k, l) = &
                                        rhs_vf(momxb + 1)%sf(j, k, l) + 5d-1* &
                                        (flux_src_n(3)%vf(momxe)%sf(j, k, l - 1) &
                                         + flux_src_n(3)%vf(momxe)%sf(j, k, l))

                                    rhs_vf(momxe)%sf(j, k, l) = &
                                        rhs_vf(momxe)%sf(j, k, l) - 5d-1* &
                                        (flux_src_n(3)%vf(momxb + 1)%sf(j, k, l - 1) &
                                         + flux_src_n(3)%vf(momxb + 1)%sf(j, k, l))
                                end do
                            end do
                        end do
                    end if
                end if

            end if  ! id loop
            call nvtxEndRange

            ! RHS additions for hypoelasticity
            call nvtxStartRange("RHS_Hypoelasticity")

            if (hypoelasticity) then

                call s_compute_hypoelastic_rhs(id, q_prim_qp%vf, rhs_vf)

            end if
            call nvtxEndRange
        end do
        ! END: Dimensional Splitting Loop =================================

        if (run_time_info .or. probe_wrt) then

            ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
            if (n > 0) iy%beg = -buff_size; 
            if (p > 0) iz%beg = -buff_size; 
            ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
            !$acc update device(ix, iy, iz)

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            q_prim_vf(i)%sf(j, k, l) = q_prim_qp%vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do

        end if

        ! ==================================================================

    end subroutine s_compute_rhs ! -----------------------------------------


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
        real(kind(0d0)) :: pres_relax
        real(kind(0d0)), dimension(num_fluids) :: pres_K_init
        real(kind(0d0)) :: f_pres
        real(kind(0d0)) :: df_pres
        real(kind(0d0)), dimension(num_fluids) :: rho_K_s
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho
        real(kind(0d0)), dimension(num_fluids) :: alpha
        real(kind(0d0)) :: sum_alpha
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)), dimension(2) :: Re

        integer :: i, j, k, l, q, iter !< Generic loop iterators
        integer :: relax !< Relaxation procedure determination variable

        !$acc parallel loop collapse(3) gang vector private(pres_K_init, rho_K_s, alpha_rho, alpha, Re, pres_relax)
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    ! Numerical correction of the volume fractions
                    if (mpp_lim) then
                        sum_alpha = 0d0

                        !$acc loop seq
                        do i = 1, num_fluids
                            if ((q_cons_vf(i + contxb - 1)%sf(j, k, l) < 0d0) .or. &
                                (q_cons_vf(i + advxb - 1)%sf(j, k, l) < 0d0)) then
                                q_cons_vf(i + contxb - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + intxb - 1)%sf(j, k, l) = 0d0
                            end if

                            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > 1d0) &
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 1d0
                            sum_alpha = sum_alpha + q_cons_vf(i + advxb - 1)%sf(j, k, l)
                        end do

                        !$acc loop seq
                        do i = 1, num_fluids
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + advxb - 1)%sf(j, k, l)/sum_alpha
                        end do
                    end if

                    ! Pressures relaxation procedure ===================================

                    ! Is the pressure relaxation procedure necessary?
                    relax = 1

                    !$acc loop seq
                    do i = 1, num_fluids
                        if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > (1d0 - sgm_eps)) relax = 0
                    end do

                    if (relax == 1) then
                        ! Initial state
                        pres_relax = 0d0

                        !$acc loop seq
                        do i = 1, num_fluids
                            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > sgm_eps) then
                                pres_K_init(i) = &
                                    (q_cons_vf(i + intxb - 1)%sf(j, k, l)/ &
                                     q_cons_vf(i + advxb - 1)%sf(j, k, l) &
                                     - pi_infs(i))/gammas(i)

                                if (pres_K_init(i) <= -(1d0 - 1d-8)*pres_inf(i) + 1d-8) &
                                    pres_K_init(i) = -(1d0 - 1d-8)*pres_inf(i) + 1d-8
                            else
                                pres_K_init(i) = 0d0
                            end if
                            pres_relax = pres_relax + q_cons_vf(i + advxb - 1)%sf(j, k, l)*pres_K_init(i)
                        end do

                        ! Iterative process for relaxed pressure determination
                        f_pres = 1d-9
                        df_pres = 1d9

                        !$acc loop seq
                        do i = 1, num_fluids
                            rho_K_s(i) = 0d0
                        end do

                        !$acc loop seq
                        do iter = 0, 49

                            if (DABS(f_pres) > 1d-10) then
                                pres_relax = pres_relax - f_pres/df_pres

                                ! Physical pressure
                                do i = 1, num_fluids
                                    if (pres_relax <= -(1d0 - 1d-8)*pres_inf(i) + 1d-8) &
                                        pres_relax = -(1d0 - 1d-8)*pres_inf(i) + 1d0
                                end do

                                ! Newton-Raphson method
                                f_pres = -1d0
                                df_pres = 0d0

                                !$acc loop seq
                                do i = 1, num_fluids
                                    if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > sgm_eps) then
                                        rho_K_s(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/ &
                                                     max(q_cons_vf(i + advxb - 1)%sf(j, k, l), sgm_eps) &
                                                     *((pres_relax + pres_inf(i))/(pres_K_init(i) + &
                                                                                   pres_inf(i)))**(1d0/gamma_min(i))

                                        f_pres = f_pres + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                                                 /rho_K_s(i)

                                        df_pres = df_pres - q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                                                  /(gamma_min(i)*rho_K_s(i)*(pres_relax + pres_inf(i)))
                                    end if
                                end do
                            end if

                        end do

                        ! Cell update of the volume fraction
                        !$acc loop seq
                        do i = 1, num_fluids
                            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > sgm_eps) &
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l) &
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

                    !$acc loop seq
                    do i = 1, num_fluids
                        alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
                        alpha(i) = q_cons_vf(E_idx + i)%sf(j, k, l)
                    end do

                    if (bubbles) then
                        rho = 0d0
                        gamma = 0d0
                        pi_inf = 0d0

                        if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                            !$acc loop seq
                            do i = 1, num_fluids
                                rho = rho + alpha_rho(i)
                                gamma = gamma + alpha(i)*gammas(i)
                                pi_inf = pi_inf + alpha(i)*pi_infs(i)
                            end do
                        else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                            !$acc loop seq
                            do i = 1, num_fluids - 1
                                rho = rho + alpha_rho(i)
                                gamma = gamma + alpha(i)*gammas(i)
                                pi_inf = pi_inf + alpha(i)*pi_infs(i)
                            end do
                        else
                            rho = alpha_rho(1)
                            gamma = gammas(1)
                            pi_inf = pi_infs(1)
                        end if
                    else
                        rho = 0d0
                        gamma = 0d0
                        pi_inf = 0d0

                        sum_alpha = 0d0

                        if (mpp_lim) then
                            !$acc loop seq
                            do i = 1, num_fluids
                                alpha_rho(i) = max(0d0, alpha_rho(i))
                                alpha(i) = min(max(0d0, alpha(i)), 1d0)
                                sum_alpha = sum_alpha + alpha(i)
                            end do

                            alpha = alpha/max(sum_alpha, sgm_eps)

                        end if

                        !$acc loop seq
                        do i = 1, num_fluids
                            rho = rho + alpha_rho(i)
                            gamma = gamma + alpha(i)*gammas(i)
                            pi_inf = pi_inf + alpha(i)*pi_infs(i)
                        end do

                        if (any(Re_size > 0)) then
                            !$acc loop seq
                            do i = 1, 2
                                Re(i) = dflt_real

                                if (Re_size(i) > 0) Re(i) = 0d0
                                !$acc loop seq
                                do q = 1, Re_size(i)
                                    Re(i) = alpha(Re_idx(i, q))/Res(i, q) &
                                            + Re(i)
                                end do

                                Re(i) = 1d0/max(Re(i), sgm_eps)

                            end do
                        end if
                    end if

                    dyn_pres = 0d0

                    !$acc loop seq
                    do i = momxb, momxe
                        dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j, k, l)* &
                                   q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
                    end do

                    pres_relax = (q_cons_vf(E_idx)%sf(j, k, l) - dyn_pres - pi_inf)/gamma

                    !$acc loop seq
                    do i = 1, num_fluids
                        q_cons_vf(i + intxb - 1)%sf(j, k, l) = &
                            q_cons_vf(i + advxb - 1)%sf(j, k, l)* &
                            (gammas(i)*pres_relax + pi_infs(i))
                    end do
                    ! ==================================================================
                end do
            end do
        end do

    end subroutine s_pressure_relaxation_procedure ! -----------------------


    !>  The purpose of this procedure is to populate the buffers
        !!      of the conservative variables, depending on the selected
        !!      boundary conditions.
    subroutine s_populate_conservative_variables_buffers() ! ---------------

        integer :: i, j, k, l, r !< Generic loop iterators

        ! Population of Buffers in x-direction =============================

        if (bc_x%beg <= -3) then         ! Ghost-cell extrap. BC at beginning

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_qp%vf(i)%sf(-j, k, l) = &
                                q_cons_qp%vf(i)%sf(0, k, l)
                        end do
                    end do
                end do
            end do

        elseif (bc_x%beg == -2) then     ! Symmetry BC at beginning

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 1, buff_size
                        !$acc loop seq
                        do i = 1, contxe
                            q_cons_qp%vf(i)%sf(-j, k, l) = &
                                q_cons_qp%vf(i)%sf(j - 1, k, l)
                        end do

                        q_cons_qp%vf(momxb)%sf(-j, k, l) = &
                            - q_cons_qp%vf(momxb)%sf(j - 1, k, l)

                        !$acc loop seq
                        do i = momxb + 1, sys_size
                            q_cons_qp%vf(i)%sf(-j, k, l) = &
                                q_cons_qp%vf(i)%sf(j - 1, k, l)
                        end do
                    end do
                end do
            end do

        elseif (bc_x%beg == -1) then     ! Periodic BC at beginning

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_qp%vf(i)%sf(-j, k, l) = &
                                q_cons_qp%vf(i)%sf(m - (j - 1), k, l)
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 1, -1)

        end if

        if (bc_x%end <= -3) then         ! Ghost-cell extrap. BC at end

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_qp%vf(i)%sf(m + j, k, l) = &
                                q_cons_qp%vf(i)%sf(m, k, l)
                        end do
                    end do
                end do
            end do

        elseif (bc_x%end == -2) then     ! Symmetry BC at end

            !$acc parallel loop collapse(3) default(present)
            do l = 0, p
                do k = 0, n
                    do j = 1, buff_size

                        !$acc loop seq
                        do i = 1, contxe
                            q_cons_qp%vf(i)%sf(m + j, k, l) = &
                                q_cons_qp%vf(i)%sf(m - (j - 1), k, l)
                        end do

                        q_cons_qp%vf(momxb)%sf(m + j, k, l) = &
                            -q_cons_qp%vf(momxb)%sf(m - (j - 1), k, l)

                        !$acc loop seq
                        do i = momxb + 1, sys_size
                            q_cons_qp%vf(i)%sf(m + j, k, l) = &
                                q_cons_qp%vf(i)%sf(m - (j - 1), k, l)
                        end do

                    end do
                end do
            end do

        elseif (bc_x%end == -1) then     ! Periodic BC at end

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            q_cons_qp%vf(i)%sf(m + j, k, l) = &
                                q_cons_qp%vf(i)%sf(j - 1, k, l)
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 1, 1)

        end if

        ! END: Population of Buffers in x-direction ========================

        ! Population of Buffers in y-direction =============================

        if (n == 0) then

            return

        elseif (bc_y%beg <= -3 .and. bc_y%beg /= -13) then     ! Ghost-cell extrap. BC at beginning

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(l, -j, k) = &
                                q_cons_qp%vf(i)%sf(l, 0, k)
                        end do
                    end do
                end do
            end do

        elseif (bc_y%beg == -13) then    ! Axis BC at beginning

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
                        if (z_cc(k) < pi) then
                            !$acc loop seq
                            do i = 1, momxb
                                q_cons_qp%vf(i)%sf(l, -j, k) = &
                                    q_cons_qp%vf(i)%sf(l, j - 1, k + ((p + 1)/2))
                            end do

                            q_cons_qp%vf(momxb + 1)%sf(l, -j, k) = &
                                -q_cons_qp%vf(momxb + 1)%sf(l, j - 1, k + ((p + 1)/2))

                            q_cons_qp%vf(momxe)%sf(l, -j, k) = &
                                -q_cons_qp%vf(momxe)%sf(l, j - 1, k + ((p + 1)/2))

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_cons_qp%vf(i)%sf(l, -j, k) = &
                                    q_cons_qp%vf(i)%sf(l, j - 1, k + ((p + 1)/2))
                            end do
                        else
                            !$acc loop seq
                            do i = 1, momxb
                                q_cons_qp%vf(i)%sf(l, -j, k) = &
                                    q_cons_qp%vf(i)%sf(l, j - 1, k - ((p + 1)/2))
                            end do

                            q_cons_qp%vf(momxb + 1)%sf(l, -j, k) = &
                                -q_cons_qp%vf(momxb + 1)%sf(l, j - 1, k - ((p + 1)/2))

                            q_cons_qp%vf(momxe)%sf(l, -j, k) = &
                                -q_cons_qp%vf(momxe)%sf(l, j - 1, k - ((p + 1)/2))

                            !$acc loop seq
                            do i = E_idx, sys_size
                                q_cons_qp%vf(i)%sf(l, -j, k) = &
                                    q_cons_qp%vf(i)%sf(l, j - 1, k - ((p + 1)/2))
                            end do
                        end if
                    end do
                end do
            end do

        elseif (bc_y%beg == -2) then     ! Symmetry BC at beginning
            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
                        !$acc loop seq
                        do i = 1, momxb
                            q_cons_qp%vf(i)%sf(l, -j, k) = &
                                q_cons_qp%vf(i)%sf(l, j - 1, k)
                        end do

                        q_cons_qp%vf(momxb + 1)%sf(l, -j, k) = &
                            - q_cons_qp%vf(momxb + 1)%sf(l, j - 1, k)

                        !$acc loop seq
                        do i = momxb + 2, sys_size
                            q_cons_qp%vf(i)%sf(l, -j, k) = &
                                q_cons_qp%vf(i)%sf(l, j - 1, k)
                        end do
                    end do
                end do
            end do

        elseif (bc_y%beg == -1) then     ! Periodic BC at beginning
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(l, -j, k) = &
                                q_cons_qp%vf(i)%sf(l, n - (j - 1), k)
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 2, -1)

        end if

        if (bc_y%end <= -3) then         ! Ghost-cell extrap. BC at end
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(l, n + j, k) = &
                                q_cons_qp%vf(i)%sf(l, n, k)
                        end do
                    end do
                end do
            end do

        elseif (bc_y%end == -2) then     ! Symmetry BC at end

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 0, p
                do j = 1, buff_size
                    do l = -buff_size, m + buff_size
                        !$acc loop seq
                        do i = 1, momxb
                            q_cons_qp%vf(i)%sf(l, n + j, k) = &
                                q_cons_qp%vf(i)%sf(l, n - (j - 1), k)
                        end do

                        q_cons_qp%vf(momxb + 1)%sf(l, n + j, k) = &
                            -q_cons_qp%vf(momxb + 1)%sf(l, n - (j - 1), k)

                        !$acc loop seq
                        do i = momxb + 2, sys_size
                            q_cons_qp%vf(i)%sf(l, n + j, k) = &
                                q_cons_qp%vf(i)%sf(l, n - (j - 1), k)
                        end do
                    end do
                end do
            end do

        elseif (bc_y%end == -1) then     ! Periodic BC at end

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(l, n + j, k) = &
                                q_cons_qp%vf(i)%sf(l, j - 1, k)
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 2, 1)

        end if

        ! END: Population of Buffers in y-direction ========================

        ! Population of Buffers in z-direction =============================

        if (p == 0) then

            return

        elseif (bc_z%beg <= -3) then     ! Ghost-cell extrap. BC at beginning

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(k, l, -j) = &
                                q_cons_qp%vf(i)%sf(k, l, 0)
                        end do
                    end do
                end do
            end do

        elseif (bc_z%beg == -2) then     ! Symmetry BC at beginning

            !$acc parallel loop collapse(3) gang vector default(present)
            do j = 1, buff_size
                do l = -buff_size, n + buff_size
                    do k = -buff_size, m + buff_size
                        !$acc loop seq
                        do i = 1, momxb + 1
                            q_cons_qp%vf(i)%sf(k, l, -j) = &
                                q_cons_qp%vf(i)%sf(k, l, j - 1)
                        end do

                        q_cons_qp%vf(momxe)%sf(k, l, -j) = &
                            -q_cons_qp%vf(momxe)%sf(k, l, j - 1)

                        !$acc loop seq
                        do i = E_idx, sys_size
                            q_cons_qp%vf(i)%sf(k, l, -j) = &
                                q_cons_qp%vf(i)%sf(k, l, j - 1)
                        end do
                    end do
                end do
            end do

        elseif (bc_z%beg == -1) then     ! Periodic BC at beginning
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(k, l, -j) = &
                                q_cons_qp%vf(i)%sf(k, l, p - (j - 1))
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 3, -1)

        end if

        if (bc_z%end <= -3) then         ! Ghost-cell extrap. BC at end
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(k, l, p + j) = &
                                q_cons_qp%vf(i)%sf(k, l, p)
                        end do
                    end do
                end do
            end do

        elseif (bc_z%end == -2) then     ! Symmetry BC at end
            !$acc parallel loop collapse(3) gang vector default(present)
            do j = 1, buff_size
                do l = -buff_size, n + buff_size
                    do k = -buff_size, m + buff_size
                        !$acc loop seq
                        do i = 1, momxb + 1
                            q_cons_qp%vf(i)%sf(k, l, p + j) = &
                                q_cons_qp%vf(i)%sf(k, l, p - (j - 1))
                        end do

                        q_cons_qp%vf(momxe)%sf(k, l, p + j) = &
                            -q_cons_qp%vf(momxe)%sf(k, l, p - (j - 1))

                        !$acc loop seq
                        do i = E_idx, sys_size
                            q_cons_qp%vf(i)%sf(k, l, p + j) = &
                                q_cons_qp%vf(i)%sf(k, l, p - (j - 1))
                        end do
                    end do
                end do
            end do

        elseif (bc_z%end == -1) then     ! Periodic BC at end
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            q_cons_qp%vf(i)%sf(k, l, p + j) = &
                                q_cons_qp%vf(i)%sf(k, l, j - 1)
                        end do
                    end do
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                q_cons_qp%vf, 3, 1)

        end if

        ! END: Population of Buffers in z-direction ========================

    end subroutine s_populate_conservative_variables_buffers ! -------------

    !>  The purpose of this subroutine is to WENO-reconstruct the
        !!      left and the right cell-boundary values, including values
        !!      at the Gaussian quadrature points, from the cell-averaged
        !!      variables.
        !!  @param v_vf Cell-average variables
        !!  @param vL_qp Left WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param vR_qp Right WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_reconstruct_cell_boundary_values(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, & ! -
                                                      norm_dir)

        type(scalar_field), dimension(iv%beg:iv%end), intent(IN) :: v_vf

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z

        integer, intent(IN) :: norm_dir

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l
        ! Reconstruction in s1-direction ===================================

        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            weno_dir = 1; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
            weno_dir = 2; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        else
            is1 = iz; is2 = iy; is3 = ix
            weno_dir = 3; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        end if

        if (n > 0) then
            if (p > 0) then

                call s_weno(v_vf(iv%beg:iv%end), &
                    vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                                norm_dir, weno_dir, &
                                is1, is2, is3)
            else
                call s_weno(v_vf(iv%beg:iv%end), &
                    vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                                norm_dir, weno_dir, &
                                is1, is2, is3)
            end if
        else

            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                            norm_dir, weno_dir, &
                            is1, is2, is3)
        end if

        ! ==================================================================
    end subroutine s_reconstruct_cell_boundary_values ! --------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_rhs_module() ! -----------------------------------

        integer :: i, j, k, l !< Generic loop iterators

        do j = cont_idx%beg, cont_idx%end
            !$acc exit data detach(q_prim_qp%vf(j)%sf)
            nullify (q_prim_qp%vf(j)%sf)
        end do

        do j = adv_idx%beg, adv_idx%end
            !$acc exit data detach(q_prim_qp%vf(j)%sf)
            nullify (q_prim_qp%vf(j)%sf)
        end do

        do j = mom_idx%beg, E_idx
            @:DEALLOCATE(q_cons_qp%vf(j)%sf)
            @:DEALLOCATE(q_prim_qp%vf(j)%sf)
        end do

        @:DEALLOCATE(q_cons_qp%vf, q_prim_qp%vf)
        @:DEALLOCATE(qL_rsx_vf, qR_rsx_vf)

        if (n > 0) then
            @:DEALLOCATE(qL_rsy_vf, qR_rsy_vf)
        end if

        if (p > 0) then
            @:DEALLOCATE(qL_rsz_vf, qR_rsz_vf)
        end if

        if (weno_Re_flux) then
            @:DEALLOCATE(dqL_rsx_vf, dqR_rsx_vf)

            if (n > 0) then
                @:DEALLOCATE(dqL_rsy_vf, dqR_rsy_vf)
            end if

            if (p > 0) then
                @:DEALLOCATE(dqL_rsz_vf, dqR_rsz_vf)
            end if
        end if

        if (mpp_lim .and. bubbles) then
            !deallocate(alf_sum%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
            !$acc exit data delete(alf_sum%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end if

        if (any(Re_size > 0)) then
            do l = mom_idx%beg, mom_idx%end
                @:DEALLOCATE(dq_prim_dx_qp%vf(l)%sf)
            end do

            if (n > 0) then

                do l = mom_idx%beg, mom_idx%end
                    @:DEALLOCATE(dq_prim_dy_qp%vf(l)%sf)
                end do

                if (p > 0) then
                    do l = mom_idx%beg, mom_idx%end
                        @:DEALLOCATE(dq_prim_dz_qp%vf(l)%sf)
                    end do
                end if

            end if

            @:DEALLOCATE(dq_prim_dx_qp%vf)
            @:DEALLOCATE(dq_prim_dy_qp%vf)
            @:DEALLOCATE(dq_prim_dz_qp%vf)
        end if

        if (any(Re_size > 0)) then
            do i = num_dims, 1, -1
                if (any(Re_size > 0)) then

                    do l = mom_idx%beg, mom_idx%end
                        @:DEALLOCATE(dqL_prim_dx_n(i)%vf(l)%sf)
                        @:DEALLOCATE(dqR_prim_dx_n(i)%vf(l)%sf)
                    end do

                    if (n > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:DEALLOCATE(dqL_prim_dy_n(i)%vf(l)%sf)
                            @:DEALLOCATE(dqR_prim_dy_n(i)%vf(l)%sf)
                        end do
                    end if

                    if (p > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:DEALLOCATE(dqL_prim_dz_n(i)%vf(l)%sf)
                            @:DEALLOCATE(dqR_prim_dz_n(i)%vf(l)%sf)
                        end do
                    end if

                end if

                @:DEALLOCATE(dqL_prim_dx_n(i)%vf)
                @:DEALLOCATE(dqL_prim_dy_n(i)%vf)
                @:DEALLOCATE(dqL_prim_dz_n(i)%vf)
                @:DEALLOCATE(dqR_prim_dx_n(i)%vf)
                @:DEALLOCATE(dqR_prim_dy_n(i)%vf)
                @:DEALLOCATE(dqR_prim_dz_n(i)%vf)
            end do
        end if

        @:DEALLOCATE(dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n)
        @:DEALLOCATE(dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n)

        if (any(Re_size > 0) .and. cyl_coord) then
            do i = 1, num_dims
                @:DEALLOCATE(tau_Re_vf(cont_idx%end + i)%sf)
            end do
            @:DEALLOCATE(tau_Re_vf(E_idx)%sf)
            @:DEALLOCATE(tau_Re_vf)
        end if

        do i = num_dims, 1, -1
            if (i /= 1) then
                do l = 1, sys_size
                    nullify (flux_n(i)%vf(l)%sf)
                    nullify (flux_src_n(i)%vf(l)%sf)
                    @:DEALLOCATE(flux_gsrc_n(i)%vf(l)%sf)
                end do
            else
                do l = 1, sys_size
                    @:DEALLOCATE(flux_n(i)%vf(l)%sf)
                    @:DEALLOCATE(flux_gsrc_n(i)%vf(l)%sf)
                end do

                if (any(Re_size > 0)) then
                    do l = mom_idx%beg, E_idx
                        @:DEALLOCATE(flux_src_n(i)%vf(l)%sf)
                    end do
                end if

                if (riemann_solver == 1) then
                    do l = adv_idx%beg + 1, adv_idx%end
                        @:DEALLOCATE(flux_src_n(i)%vf(l)%sf)
                    end do
                else
                    do l = adv_idx%beg + 1, adv_idx%end
                        nullify (flux_src_n(i)%vf(l)%sf)
                    end do
                end if

                @:DEALLOCATE(flux_src_n(i)%vf(adv_idx%beg)%sf)
            end if

            @:DEALLOCATE(flux_n(i)%vf, flux_src_n(i)%vf, flux_gsrc_n(i)%vf)
        end do

        @:DEALLOCATE(flux_n, flux_src_n, flux_gsrc_n)

        s_riemann_solver => null()
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_rhs_module ! ---------------------------------

end module m_rhs
