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
!! @file m_cbc.f90
!! @brief Contains module m_cbc
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The module features a large database of characteristic boundary
!!              conditions (CBC) for the Euler system of equations. This system
!!              is augmented by the appropriate advection equations utilized to
!!              capture the material interfaces. The closure is achieved by the
!!              stiffened equation of state and mixture relations. At this time,
!!              the following CBC are available:
!!                           1) Slip Wall
!!                           2) Nonreflecting Subsonic Buffer
!!                           3) Nonreflecting Subsonic Inflow
!!                           4) Nonreflecting Subsonic Outflow
!!                           5) Force-Free Subsonic Outflow
!!                           6) Constant Pressure Subsonic Outflow
!!                           7) Supersonic Inflow
!!                           8) Supersonic Outflow
!!              Please refer to Thompson (1987, 1990) for detailed descriptions.
module m_cbc

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_cbc_module, s_cbc, s_finalize_cbc_module

    abstract interface ! =======================================================

        !> Abstract interface to the procedures that are utilized to calculate
        !! the L variables. For additional information refer to the following:
        !!            1) s_compute_slip_wall_L
        !!            2) s_compute_nonreflecting_subsonic_buffer_L
        !!            3) s_compute_nonreflecting_subsonic_inflow_L
        !!            4) s_compute_nonreflecting_subsonic_outflow_L
        !!            5) s_compute_force_free_subsonic_outflow_L
        !!            6) s_compute_constant_pressure_subsonic_outflow_L
        !!            7) s_compute_supersonic_inflow_L
        !!            8) s_compute_supersonic_outflow_L
        !! @param dflt_int Default null integer
        subroutine s_compute_abstract_L(dflt_int)

            integer, intent(IN) :: dflt_int

        end subroutine s_compute_abstract_L

    end interface ! ============================================================

    type(scalar_field), allocatable, dimension(:) :: q_prim_rs_vf !<
    !! The cell-average primitive variables. They are obtained by reshaping (RS)
    !! q_prim_vf in the coordinate direction normal to the domain boundary along
    !! which the CBC is applied.

    type(scalar_field), allocatable, dimension(:) :: F_rs_vf, F_src_rs_vf !<
    !! Cell-average fluxes (src - source). These are directly determined from the
    !! cell-average primitive variables, q_prims_rs_vf, and not a Riemann solver.

    type(scalar_field), allocatable, dimension(:) :: flux_rs_vf, flux_src_rs_vf !<
    !! The cell-boundary-average of the fluxes. They are initially determined by
    !! reshaping flux_vf and flux_src_vf in a coordinate direction normal to the
    !! domain boundary along which CBC is applied. flux_rs_vf and flux_src_rs_vf
    !! are subsequently modified based on the selected CBC.

    real(kind(0d0)), allocatable, dimension(:)   :: alpha_rho   !< Cell averaged partial densiy
    real(kind(0d0))                              :: rho         !< Cell averaged density
    real(kind(0d0)), allocatable, dimension(:)   :: vel         !< Cell averaged velocity
    real(kind(0d0))                              :: pres        !< Cell averaged pressure
    real(kind(0d0))                              :: E           !< Cell averaged energy
    real(kind(0d0))                              :: H           !< Cell averaged enthalpy
    real(kind(0d0)), allocatable, dimension(:)   :: adv         !< Cell averaged advected variables
    real(kind(0d0)), allocatable, dimension(:)   :: mf          !< Cell averaged mass fraction
    real(kind(0d0))                              :: gamma       !< Cell averaged specific heat ratio
    real(kind(0d0))                              :: pi_inf      !< Cell averaged liquid stiffness
    real(kind(0d0))                              :: c           !< Cell averaged speed of sound
    real(kind(0d0)), dimension(2)   :: Re          !< Cell averaged Reynolds numbers
    real(kind(0d0)), allocatable, dimension(:, :) :: We          !< Cell averaged Weber numbers

    real(kind(0d0)), allocatable, dimension(:) :: dalpha_rho_ds !< Spatial derivatives in s-dir of partial density
    real(kind(0d0)), allocatable, dimension(:) ::       dvel_ds !< Spatial derivatives in s-dir of velocity
    real(kind(0d0))                            ::      dpres_ds !< Spatial derivatives in s-dir of pressure
    real(kind(0d0)), allocatable, dimension(:) ::       dadv_ds !< Spatial derivatives in s-dir of advection variables
    !! Note that these are only obtained in those cells on the domain boundary along which the
    !! CBC is applied by employing finite differences (FD) on the cell-average primitive variables, q_prim_rs_vf.

    real(kind(0d0)), dimension(3) :: lambda !< Eigenvalues (see Thompson 1987,1990)
    real(kind(0d0)), allocatable, dimension(:) :: L      !< L matrix (see Thompson 1987,1990)

    real(kind(0d0)), allocatable, dimension(:) :: ds !< Cell-width distribution in the s-direction

    ! CBC Coefficients =========================================================
    real(kind(0d0)), target, allocatable, dimension(:, :) :: fd_coef_x !< Finite diff. coefficients x-dir
    real(kind(0d0)), target, allocatable, dimension(:, :) :: fd_coef_y !< Finite diff. coefficients y-dir
    real(kind(0d0)), target, allocatable, dimension(:, :) :: fd_coef_z !< Finite diff. coefficients z-dir
    !! The first dimension identifies the location of a coefficient in the FD
    !! formula, while the last dimension denotes the location of the CBC.

    real(kind(0d0)), pointer, dimension(:, :) :: fd_coef => null()

    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: pi_coef_x !< Polynominal interpolant coefficients in x-dir
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: pi_coef_y !< Polynominal interpolant coefficients in y-dir
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: pi_coef_z !< Polynominal interpolant coefficients in z-dir
    !! The first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the location of the CBC.

    real(kind(0d0)), pointer, dimension(:, :, :) :: pi_coef => null()
    ! ==========================================================================

    procedure(s_compute_abstract_L), pointer :: s_compute_L => null() !<
    !! Pointer to procedure used to calculate L variables, based on choice of CBC

    type(bounds_info) :: is1, is2, is3 !< Indical bounds in the s1-, s2- and s3-directions

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_cbc_module() ! ---------------------------------

        if (all((/bc_x%beg, bc_x%end/) > -5) &
            .and. &
            (n > 0 .and. all((/bc_y%beg, bc_y%end/) > -5)) &
            .and. &
            (p > 0 .and. all((/bc_z%beg, bc_z%end/) > -5))) return

        ! Allocating the cell-average primitive variables
        allocate (q_prim_rs_vf(1:sys_size))

        ! Allocating the cell-average and cell-boundary-average fluxes
        allocate (F_rs_vf(1:sys_size), F_src_rs_vf(1:sys_size))
        allocate (flux_rs_vf(1:sys_size), flux_src_rs_vf(1:sys_size))

        ! Allocating the cell-average partial densities, the velocity, the
        ! advected variables, the mass fractions, as well as Weber numbers
        allocate (alpha_rho(1:cont_idx%end))
        allocate (vel(1:num_dims))
        allocate (adv(1:adv_idx%end - E_idx))
        allocate (mf(1:cont_idx%end))

        allocate (We(1:num_fluids, 1:num_fluids))

        ! Allocating the first-order spatial derivatives in the s-direction
        ! of the partial densities, the velocity and the advected variables
        allocate (dalpha_rho_ds(1:cont_idx%end))
        allocate (dvel_ds(1:num_dims))
        allocate (dadv_ds(1:adv_idx%end - E_idx))

        ! Allocating L, see Thompson (1987, 1990)
        allocate (L(1:adv_idx%end))

        ! Allocating the cell-width distribution in the s-direction
        allocate (ds(0:buff_size))

        ! Allocating/Computing CBC Coefficients in x-direction =============
        if (all((/bc_x%beg, bc_x%end/) <= -5)) then

            allocate (fd_coef_x(0:buff_size, -1:1))

            if (weno_order > 1) then
                allocate (pi_coef_x(0:weno_polyn - 1, 0:weno_order - 3, -1:1))
            end if

            call s_compute_cbc_coefficients(1, -1)
            call s_compute_cbc_coefficients(1, 1)

        elseif (bc_x%beg <= -5) then

            allocate (fd_coef_x(0:buff_size, -1:-1))

            if (weno_order > 1) then
                allocate (pi_coef_x(0:weno_polyn - 1, 0:weno_order - 3, -1:-1))
            end if

            call s_compute_cbc_coefficients(1, -1)

        elseif (bc_x%end <= -5) then

            allocate (fd_coef_x(0:buff_size, 1:1))

            if (weno_order > 1) then
                allocate (pi_coef_x(0:weno_polyn - 1, 0:weno_order - 3, 1:1))
            end if

            call s_compute_cbc_coefficients(1, 1)

        end if
        ! ==================================================================

        ! Allocating/Computing CBC Coefficients in y-direction =============
        if (n > 0) then

            if (all((/bc_y%beg, bc_y%end/) <= -5)) then

                allocate (fd_coef_y(0:buff_size, -1:1))

                if (weno_order > 1) then
                    allocate (pi_coef_y(0:weno_polyn - 1, 0:weno_order - 3, -1:1))
                end if

                call s_compute_cbc_coefficients(2, -1)
                call s_compute_cbc_coefficients(2, 1)

            elseif (bc_y%beg <= -5) then

                allocate (fd_coef_y(0:buff_size, -1:-1))

                if (weno_order > 1) then
                    allocate (pi_coef_y(0:weno_polyn - 1, 0:weno_order - 3, -1:-1))
                end if

                call s_compute_cbc_coefficients(2, -1)

            elseif (bc_y%end <= -5) then

                allocate (fd_coef_y(0:buff_size, 1:1))

                if (weno_order > 1) then
                    allocate (pi_coef_y(0:weno_polyn - 1, 0:weno_order - 3, 1:1))
                end if

                call s_compute_cbc_coefficients(2, 1)

            end if

        end if
        ! ==================================================================

        ! Allocating/Computing CBC Coefficients in z-direction =============
        if (p > 0) then

            if (all((/bc_z%beg, bc_z%end/) <= -5)) then

                allocate (fd_coef_z(0:buff_size, -1:1))

                if (weno_order > 1) then
                    allocate (pi_coef_z(0:weno_polyn - 1, 0:weno_order - 3, -1:1))
                end if

                call s_compute_cbc_coefficients(3, -1)
                call s_compute_cbc_coefficients(3, 1)

            elseif (bc_z%beg <= -5) then

                allocate (fd_coef_z(0:buff_size, -1:-1))

                if (weno_order > 1) then
                    allocate (pi_coef_z(0:weno_polyn - 1, 0:weno_order - 3, -1:-1))
                end if

                call s_compute_cbc_coefficients(3, -1)

            elseif (bc_z%end <= -5) then

                allocate (fd_coef_z(0:buff_size, 1:1))

                if (weno_order > 1) then
                    allocate (pi_coef_z(0:weno_polyn - 1, 0:weno_order - 3, 1:1))
                end if

                call s_compute_cbc_coefficients(3, 1)

            end if

        end if
        ! ==================================================================

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables
        elseif (bubbles) then                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables_bubbles
        else                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if

    end subroutine s_initialize_cbc_module ! -------------------------------

    subroutine s_compute_cbc_coefficients(cbc_dir, cbc_loc) ! --------------
        ! Description: The purpose of this subroutine is to compute the grid
        !              dependent FD and PI coefficients, or CBC coefficients,
        !              provided the CBC coordinate direction and location.

        ! CBC coordinate direction and location
        integer, intent(IN) :: cbc_dir, cbc_loc

        ! Cell-boundary locations in the s-direction
        real(kind(0d0)), dimension(0:buff_size + 1) :: s_cb

        ! Generic loop iterator
        integer :: i

        ! Associating CBC coefficients pointers
        call s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc)

        ! Determining the cell-boundary locations in the s-direction
        s_cb(0) = 0d0

        do i = 0, buff_size
            s_cb(i + 1) = s_cb(i) + ds(i)
        end do

        ! Computing CBC1 Coefficients ======================================
        if (weno_order == 1) then

            fd_coef(:, cbc_loc) = 0d0
            fd_coef(0, cbc_loc) = -2d0/(ds(0) + ds(1))
            fd_coef(1, cbc_loc) = -fd_coef(0, cbc_loc)

            ! ==================================================================

            ! Computing CBC2 Coefficients ======================================
        elseif (weno_order == 3) then

            fd_coef(:, cbc_loc) = 0d0
            fd_coef(0, cbc_loc) = -6d0/(3d0*ds(0) + 2d0*ds(1) - ds(2))
            fd_coef(1, cbc_loc) = -4d0*fd_coef(0, cbc_loc)/3d0
            fd_coef(2, cbc_loc) = fd_coef(0, cbc_loc)/3d0

            pi_coef(0, 0, cbc_loc) = (s_cb(0) - s_cb(1))/(s_cb(0) - s_cb(2))

            ! ==================================================================

            ! Computing CBC4 Coefficients ======================================
        else

            fd_coef(:, cbc_loc) = 0d0
            fd_coef(0, cbc_loc) = -50d0/(25d0*ds(0) + 2d0*ds(1) &
                                         - 1d1*ds(2) + 1d1*ds(3) &
                                         - 3d0*ds(4))
            fd_coef(1, cbc_loc) = -48d0*fd_coef(0, cbc_loc)/25d0
            fd_coef(2, cbc_loc) = 36d0*fd_coef(0, cbc_loc)/25d0
            fd_coef(3, cbc_loc) = -16d0*fd_coef(0, cbc_loc)/25d0
            fd_coef(4, cbc_loc) = 3d0*fd_coef(0, cbc_loc)/25d0

            pi_coef(0, 0, cbc_loc) = &
                ((s_cb(0) - s_cb(1))*(s_cb(1) - s_cb(2))* &
                 (s_cb(1) - s_cb(3)))/((s_cb(1) - s_cb(4))* &
                                       (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(2)))
            pi_coef(0, 1, cbc_loc) = &
                ((s_cb(1) - s_cb(0))*(s_cb(1) - s_cb(2))* &
                 ((s_cb(1) - s_cb(3))*(s_cb(1) - s_cb(3)) - &
                  (s_cb(0) - s_cb(4))*((s_cb(3) - s_cb(1)) + &
                                       (s_cb(4) - s_cb(1)))))/ &
                ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                 (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
            pi_coef(0, 2, cbc_loc) = &
                (s_cb(1) - s_cb(0))*((s_cb(1) - s_cb(2))* &
                                     (s_cb(1) - s_cb(3)) + ((s_cb(0) - s_cb(2)) + &
                                                            (s_cb(1) - s_cb(3)))*(s_cb(0) - s_cb(4)))/ &
                ((s_cb(2) - s_cb(0))*(s_cb(0) - s_cb(3))* &
                 (s_cb(0) - s_cb(4)))
            pi_coef(1, 0, cbc_loc) = &
                ((s_cb(0) - s_cb(2))*(s_cb(2) - s_cb(1))* &
                 (s_cb(2) - s_cb(3)))/((s_cb(2) - s_cb(4))* &
                                       (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(1)))
            pi_coef(1, 1, cbc_loc) = &
                ((s_cb(0) - s_cb(2))*(s_cb(1) - s_cb(2))* &
                 ((s_cb(1) - s_cb(3))*(s_cb(2) - s_cb(3)) + &
                  (s_cb(0) - s_cb(4))*((s_cb(1) - s_cb(3)) + &
                                       (s_cb(2) - s_cb(4)))))/ &
                ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                 (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
            pi_coef(1, 2, cbc_loc) = &
                ((s_cb(1) - s_cb(2))*(s_cb(2) - s_cb(3))* &
                 (s_cb(2) - s_cb(4)))/((s_cb(0) - s_cb(2))* &
                                       (s_cb(0) - s_cb(3))*(s_cb(0) - s_cb(4)))

        end if
        ! END: Computing CBC4 Coefficients =================================

        ! Nullifying CBC coefficients
        nullify (fd_coef, pi_coef)

    end subroutine s_compute_cbc_coefficients ! ----------------------------

        !!  The goal of the procedure is to associate the FD and PI
        !!      coefficients, or CBC coefficients, with the appropriate
        !!      targets, based on the coordinate direction and location
        !!      of the CBC.
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
    subroutine s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc) ! ---

        integer, intent(IN) :: cbc_dir, cbc_loc

        integer :: i !< Generic loop iterator

        ! Associating CBC Coefficients in x-direction ======================
        if (cbc_dir == 1) then

            fd_coef => fd_coef_x; if (weno_order > 1) pi_coef => pi_coef_x

            if (cbc_loc == -1) then
                do i = 0, buff_size
                    ds(i) = dx(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dx(m - i)
                end do
            end if
            ! ==================================================================

            ! Associating CBC Coefficients in y-direction ======================
        elseif (cbc_dir == 2) then

            fd_coef => fd_coef_y; if (weno_order > 1) pi_coef => pi_coef_y

            if (cbc_loc == -1) then
                do i = 0, buff_size
                    ds(i) = dy(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dy(n - i)
                end do
            end if
            ! ==================================================================

            ! Associating CBC Coefficients in z-direction ======================
        else

            fd_coef => fd_coef_z; if (weno_order > 1) pi_coef => pi_coef_z

            if (cbc_loc == -1) then
                do i = 0, buff_size
                    ds(i) = dz(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dz(p - i)
                end do
            end if

        end if
        ! ==================================================================

    end subroutine s_associate_cbc_coefficients_pointers ! -----------------

    !>  The following is the implementation of the CBC based on
        !!      the work of Thompson (1987, 1990) on hyperbolic systems.
        !!      The CBC is indirectly applied in the computation of the
        !!      right-hand-side (RHS) near the relevant domain boundary
        !!      through the modification of the fluxes.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_cbc(q_prim_vf, flux_vf, flux_src_vf, & ! -----------------
                     cbc_dir, cbc_loc, &
                     ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf

        integer, intent(IN) :: cbc_dir, cbc_loc

        type(bounds_info), intent(IN) :: ix, iy, iz

        ! First-order time derivatives of the partial densities, density,
        ! velocity, pressure, advection variables, and the specific heat
        ! ratio and liquid stiffness functions
        real(kind(0d0)), dimension(cont_idx%end)   :: dalpha_rho_dt
        real(kind(0d0))                            ::       drho_dt
        real(kind(0d0)), dimension(num_dims)       ::       dvel_dt
        real(kind(0d0))                            ::      dpres_dt
        real(kind(0d0)), dimension(adv_idx%end - E_idx) ::       dadv_dt
        real(kind(0d0))                            ::     dgamma_dt
        real(kind(0d0))                            ::    dpi_inf_dt

        integer :: i, j, k, r !< Generic loop iterators

        real(kind(0d0)) :: blkmod1, blkmod2 !< Fluid bulk modulus for Wood mixture sound speed

        ! Reshaping of inputted data and association of the FD and PI
        ! coefficients, or CBC coefficients, respectively, hinging on
        ! selected CBC coordinate direction
        call s_initialize_cbc(q_prim_vf, flux_vf, flux_src_vf, &
                              cbc_dir, cbc_loc, &
                              ix, iy, iz)

        call s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc)

        ! PI2 of flux_rs_vf and flux_src_rs_vf at j = 1/2 ==================
        if (weno_order == 3) then

            call s_convert_primitive_to_flux_variables(q_prim_rs_vf, &
                                                       F_rs_vf, &
                                                       F_src_rs_vf, &
                                                       is1, is2, is3)

            do i = 1, adv_idx%end
                flux_rs_vf(i)%sf(0, :, :) = F_rs_vf(i)%sf(0, :, :) &
                                            + pi_coef(0, 0, cbc_loc)* &
                                            (F_rs_vf(i)%sf(1, :, :) - &
                                             F_rs_vf(i)%sf(0, :, :))
            end do

            do i = adv_idx%beg, adv_idx%end
                flux_src_rs_vf(i)%sf(0, :, :) = F_src_rs_vf(i)%sf(0, :, :) + &
                                                (F_src_rs_vf(i)%sf(1, :, :) - &
                                                 F_src_rs_vf(i)%sf(0, :, :)) &
                                                *pi_coef(0, 0, cbc_loc)
            end do
            ! ==================================================================

            ! PI4 of flux_rs_vf and flux_src_rs_vf at j = 1/2, 3/2 =============
        elseif (weno_order == 5) then

            call s_convert_primitive_to_flux_variables(q_prim_rs_vf, &
                                                       F_rs_vf, &
                                                       F_src_rs_vf, &
                                                       is1, is2, is3)

            do i = 1, adv_idx%end
                do j = 0, 1
                    flux_rs_vf(i)%sf(j, :, :) = F_rs_vf(i)%sf(j, :, :) &
                                                + pi_coef(j, 0, cbc_loc)* &
                                                (F_rs_vf(i)%sf(3, :, :) - &
                                                 F_rs_vf(i)%sf(2, :, :)) &
                                                + pi_coef(j, 1, cbc_loc)* &
                                                (F_rs_vf(i)%sf(2, :, :) - &
                                                 F_rs_vf(i)%sf(1, :, :)) &
                                                + pi_coef(j, 2, cbc_loc)* &
                                                (F_rs_vf(i)%sf(1, :, :) - &
                                                 F_rs_vf(i)%sf(0, :, :))
                end do
            end do

            do i = adv_idx%beg, adv_idx%end
                do j = 0, 1
                    flux_src_rs_vf(i)%sf(j, :, :) = F_src_rs_vf(i)%sf(j, :, :) + &
                                                    (F_src_rs_vf(i)%sf(3, :, :) - &
                                                     F_src_rs_vf(i)%sf(2, :, :)) &
                                                    *pi_coef(j, 0, cbc_loc) + &
                                                    (F_src_rs_vf(i)%sf(2, :, :) - &
                                                     F_src_rs_vf(i)%sf(1, :, :)) &
                                                    *pi_coef(j, 1, cbc_loc) + &
                                                    (F_src_rs_vf(i)%sf(1, :, :) - &
                                                     F_src_rs_vf(i)%sf(0, :, :)) &
                                                    *pi_coef(j, 2, cbc_loc)
                end do
            end do

        end if
        ! ==================================================================

        ! FD2 or FD4 of RHS at j = 0 =======================================
        do r = is3%beg, is3%end
            do k = is2%beg, is2%end

                ! Transferring the Primitive Variables =======================
                do i = 1, cont_idx%end
                    alpha_rho(i) = q_prim_rs_vf(i)%sf(0, k, r)
                end do

                do i = 1, num_dims
                    vel(i) = q_prim_rs_vf(cont_idx%end + i)%sf(0, k, r)
                end do

                pres = q_prim_rs_vf(E_idx)%sf(0, k, r)

                call s_convert_to_mixture_variables(q_prim_rs_vf, &
                                                    rho, gamma, &
                                                    pi_inf, Re, &
                                                    We, 0, k, r)

                E = gamma*pres + pi_inf + 5d-1*rho*sum(vel**2d0)

                H = (E + pres)/rho

                do i = 1, adv_idx%end - E_idx
                    adv(i) = q_prim_rs_vf(E_idx + i)%sf(0, k, r)
                end do

                mf = alpha_rho/rho

                ! Compute mixture sound speed
                if (alt_soundspeed .or. regularization) then
                    blkmod1 = ((fluid_pp(1)%gamma + 1d0)*pres + &
                               fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                    blkmod2 = ((fluid_pp(2)%gamma + 1d0)*pres + &
                               fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                    c = (1d0/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))
                elseif (model_eqns == 3) then
                    c = 0d0
                    do i = 1, num_fluids
                        c = c + q_prim_rs_vf(i + adv_idx%beg - 1)%sf(0, k, r)*(1d0/fluid_pp(i)%gamma + 1d0)* &
                            (pres + fluid_pp(i)%pi_inf/(fluid_pp(i)%gamma + 1d0))
                    end do
                    c = c/rho
                else
                    c = ((H - 5d-1*sum(vel**2d0))/gamma)
                end if

                c = sqrt(c)

!                  IF (mixture_err .AND. c < 0d0) THEN
!                    c = sgm_eps
!                  ELSE
!                    c = SQRT(c)
!                  END IF

                ! ============================================================

                ! First-Order Spatial Derivatives of Primitive Variables =====
                dalpha_rho_ds = 0d0
                dvel_ds = 0d0
                dpres_ds = 0d0
                dadv_ds = 0d0

                do j = 0, buff_size

                    do i = 1, cont_idx%end
                        dalpha_rho_ds(i) = q_prim_rs_vf(i)%sf(j, k, r)* &
                                           fd_coef(j, cbc_loc) + &
                                           dalpha_rho_ds(i)
                    end do

                    do i = 1, num_dims
                        dvel_ds(i) = q_prim_rs_vf(cont_idx%end + i)%sf(j, k, r)* &
                                     fd_coef(j, cbc_loc) + &
                                     dvel_ds(i)
                    end do

                    dpres_ds = q_prim_rs_vf(E_idx)%sf(j, k, r)* &
                               fd_coef(j, cbc_loc) + &
                               dpres_ds

                    do i = 1, adv_idx%end - E_idx
                        dadv_ds(i) = q_prim_rs_vf(E_idx + i)%sf(j, k, r)* &
                                     fd_coef(j, cbc_loc) + &
                                     dadv_ds(i)
                    end do

                end do
                ! ============================================================

                ! First-Order Temporal Derivatives of Primitive Variables ====
                lambda(1) = vel(dir_idx(1)) - c
                lambda(2) = vel(dir_idx(1))
                lambda(3) = vel(dir_idx(1)) + c

                call s_compute_L(dflt_int)

                ! Be careful about the cylindrical coordinate!
                if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
                    dpres_dt = -5d-1*(L(adv_idx%end) + L(1)) + rho*c*c*vel(dir_idx(1)) &
                               /y_cc(n)
                else
                    dpres_dt = -5d-1*(L(adv_idx%end) + L(1))
                end if

                do i = 1, cont_idx%end
                    dalpha_rho_dt(i) = &
                        -(L(i + 1) - mf(i)*dpres_dt)/(c*c)
                end do

                do i = 1, num_dims
                    dvel_dt(dir_idx(i)) = dir_flg(dir_idx(i))* &
                                          (L(1) - L(adv_idx%end))/(2d0*rho*c) + &
                                          (dir_flg(dir_idx(i)) - 1d0)* &
                                          L(mom_idx%beg + i)
                end do

                ! The treatment of void fraction source is unclear
                if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
                    do i = 1, adv_idx%end - E_idx
                        dadv_dt(i) = -L(mom_idx%end + i) !+ adv(i) * vel(dir_idx(1))/y_cc(n)
                    end do
                else
                    do i = 1, adv_idx%end - E_idx
                        dadv_dt(i) = -L(mom_idx%end + i)
                    end do
                end if

                drho_dt = 0d0; dgamma_dt = 0d0; dpi_inf_dt = 0d0

                if (model_eqns == 1) then
                    drho_dt = dalpha_rho_dt(1)
                    dgamma_dt = dadv_dt(1)
                    dpi_inf_dt = dadv_dt(2)
                else
                    do i = 1, num_fluids
                        drho_dt = drho_dt + dalpha_rho_dt(i)
                        dgamma_dt = dgamma_dt + dadv_dt(i)*fluid_pp(i)%gamma
                        dpi_inf_dt = dpi_inf_dt + dadv_dt(i)*fluid_pp(i)%pi_inf
                    end do
                end if
                ! ============================================================

                ! flux_rs_vf and flux_src_rs_vf at j = -1/2 ==================
                do i = 1, cont_idx%end
                    flux_rs_vf(i)%sf(-1, k, r) = flux_rs_vf(i)%sf(0, k, r) &
                                                 + ds(0)*dalpha_rho_dt(i)
                end do

                do i = mom_idx%beg, mom_idx%end
                    flux_rs_vf(i)%sf(-1, k, r) = flux_rs_vf(i)%sf(0, k, r) &
                                                 + ds(0)*(vel(i - cont_idx%end)*drho_dt &
                                                          + rho*dvel_dt(i - cont_idx%end))
                end do

                flux_rs_vf(E_idx)%sf(-1, k, r) = flux_rs_vf(E_idx)%sf(0, k, r) &
                                                 + ds(0)*(pres*dgamma_dt &
                                                          + gamma*dpres_dt &
                                                          + dpi_inf_dt &
                                                          + rho*sum(vel*dvel_dt) &
                                                          + 5d-1*drho_dt*sum(vel**2d0))

                if (riemann_solver == 1) then

                    do i = adv_idx%beg, adv_idx%end
                        flux_rs_vf(i)%sf(-1, k, r) = 0d0
                    end do

                    do i = adv_idx%beg, adv_idx%end
                        flux_src_rs_vf(i)%sf(-1, k, r) = &
                            1d0/max(abs(vel(dir_idx(1))), sgm_eps) &
                            *sign(1d0, vel(dir_idx(1))) &
                            *(flux_rs_vf(i)%sf(0, k, r) &
                              + vel(dir_idx(1)) &
                              *flux_src_rs_vf(i)%sf(0, k, r) &
                              + ds(0)*dadv_dt(i - E_idx))
                    end do

                else

                    do i = adv_idx%beg, adv_idx%end
                        flux_rs_vf(i)%sf(-1, k, r) = flux_rs_vf(i)%sf(0, k, r) - &
                                                     adv(i - E_idx)*flux_src_rs_vf(i)%sf(0, k, r) + &
                                                     ds(0)*dadv_dt(i - E_idx)
                    end do

                    do i = adv_idx%beg, adv_idx%end
                        flux_src_rs_vf(i)%sf(-1, k, r) = 0d0
                    end do

                end if
                ! END: flux_rs_vf and flux_src_rs_vf at j = -1/2 =============

            end do
        end do
        ! END: FD2 or FD4 of RHS at j = 0 ==================================

        ! The reshaping of outputted data and disssociation of the FD and PI
        ! coefficients, or CBC coefficients, respectively, based on selected
        ! CBC coordinate direction.
        call s_finalize_cbc(flux_vf, flux_src_vf, &
                            cbc_dir, cbc_loc, &
                            ix, iy, iz)

        nullify (fd_coef, pi_coef)

    end subroutine s_cbc ! -------------------------------------------------

    !>  The L variables for the slip wall CBC, see pg. 451 of
        !!      Thompson (1990). At the slip wall (frictionless wall),
        !!      the normal component of velocity is zero at all times,
        !!      while the transverse velocities may be nonzero.
        !!  @param dflt_int Default null integer
    subroutine s_compute_slip_wall_L(dflt_int) ! -----------------------------------

        integer, intent(IN) :: dflt_int

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        L(2:adv_idx%end - 1) = 0d0

        L(adv_idx%end) = L(1)

    end subroutine s_compute_slip_wall_L ! ---------------------------------

    !>  The L variables for the nonreflecting subsonic buffer CBC
        !!      see pg. 13 of Thompson (1987). The nonreflecting subsonic
        !!      buffer reduces the amplitude of any reflections caused by
        !!      outgoing waves.
        !!  @param dflt_int Default null integer
    subroutine s_compute_nonreflecting_subsonic_buffer_L(dflt_int) ! ---------------

        integer, intent(IN) :: dflt_int

        integer :: i !< Generic loop iterator

        L(1) = (5d-1 - 5d-1*sign(1d0, lambda(1)))*lambda(1) &
               *(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, mom_idx%beg
            L(i) = (5d-1 - 5d-1*sign(1d0, lambda(2)))*lambda(2) &
                   *(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = mom_idx%beg + 1, mom_idx%end
            L(i) = (5d-1 - 5d-1*sign(1d0, lambda(2)))*lambda(2) &
                   *(dvel_ds(dir_idx(i - cont_idx%end)))
        end do

        do i = E_idx, adv_idx%end - 1
            L(i) = (5d-1 - 5d-1*sign(1d0, lambda(2)))*lambda(2) &
                   *(dadv_ds(i - mom_idx%end))
        end do

        L(adv_idx%end) = (5d-1 - 5d-1*sign(1d0, lambda(3)))*lambda(3) &
                         *(dpres_ds + rho*c*dvel_ds(dir_idx(1)))

    end subroutine s_compute_nonreflecting_subsonic_buffer_L ! -------------

    !>  The L variables for the nonreflecting subsonic inflow CBC
        !!      see pg. 455, Thompson (1990). This nonreflecting subsonic
        !!      CBC assumes an incoming flow and reduces the amplitude of
        !!      any reflections caused by outgoing waves.
        !! @param dflt_int Default null integer
    subroutine s_compute_nonreflecting_subsonic_inflow_L(dflt_int) ! ---------------

        integer, intent(IN) :: dflt_int

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        L(2:adv_idx%end) = 0d0

    end subroutine s_compute_nonreflecting_subsonic_inflow_L ! -------------

    !>  The L variables for the nonreflecting subsonic outflow
        !!      CBC see pg. 454 of Thompson (1990). This nonreflecting
        !!      subsonic CBC presumes an outgoing flow and reduces the
        !!      amplitude of any reflections caused by outgoing waves.
        !! @param dflt_int Default null integer
    subroutine s_compute_nonreflecting_subsonic_outflow_L(dflt_int) ! --------------

        integer, intent(IN) :: dflt_int

        integer :: i !> Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, mom_idx%beg
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = mom_idx%beg + 1, mom_idx%end
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - cont_idx%end)))
        end do

        do i = E_idx, adv_idx%end - 1
            L(i) = lambda(2)*(dadv_ds(i - mom_idx%end))
        end do

        ! bubble index
        L(adv_idx%end) = 0d0

    end subroutine s_compute_nonreflecting_subsonic_outflow_L ! ------------

    !>  The L variables for the force-free subsonic outflow CBC,
        !!      see pg. 454 of Thompson (1990). The force-free subsonic
        !!      outflow sets to zero the sum of all of the forces which
        !!      are acting on a fluid element for the normal coordinate
        !!      direction to the boundary. As a result, a fluid element
        !!      at the boundary is simply advected outward at the fluid
        !!      velocity.
        !! @param dflt_int Default null integer
    subroutine s_compute_force_free_subsonic_outflow_L(dflt_int) ! -----------------

        integer, intent(IN) :: dflt_int

        integer :: i !> Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, mom_idx%beg
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = mom_idx%beg + 1, mom_idx%end
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - cont_idx%end)))
        end do

        do i = E_idx, adv_idx%end - 1
            L(i) = lambda(2)*(dadv_ds(i - mom_idx%end))
        end do

        L(adv_idx%end) = L(1) + 2d0*rho*c*lambda(2)*dvel_ds(dir_idx(1))

    end subroutine s_compute_force_free_subsonic_outflow_L ! ---------------

    !>  L variables for the constant pressure subsonic outflow
        !!      CBC see pg. 455 Thompson (1990). The constant pressure
        !!      subsonic outflow maintains a fixed pressure at the CBC
        !!      boundary in absence of any transverse effects.
        !! @param dflt_int Default null integer
    subroutine s_compute_constant_pressure_subsonic_outflow_L(dflt_int) ! ----------

        integer, intent(IN) :: dflt_int

        integer :: i !> Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, mom_idx%beg
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = mom_idx%beg + 1, mom_idx%end
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - cont_idx%end)))
        end do

        do i = E_idx, adv_idx%end - 1
            L(i) = lambda(2)*(dadv_ds(i - mom_idx%end))
        end do

        L(adv_idx%end) = -L(1)

    end subroutine s_compute_constant_pressure_subsonic_outflow_L ! --------

    !>  L variables for the supersonic inflow CBC, see pg. 453
        !!      Thompson (1990). The supersonic inflow CBC is a steady
        !!      state, or nearly a steady state, CBC in which only the
        !!      transverse terms may generate a time dependence at the
        !!      inflow boundary.
        !! @param dflt_int Default null integer
    subroutine s_compute_supersonic_inflow_L(dflt_int) ! ---------------------------

        integer, intent(IN) :: dflt_int

        L = 0d0

    end subroutine s_compute_supersonic_inflow_L ! -------------------------

    !>  L variables for the supersonic outflow CBC, see pg. 453
        !!      of Thompson (1990). For the supersonic outflow CBC, the
        !!      flow evolution at the boundary is determined completely
        !!      by the interior data.
        !! @param dflt_int Default null integer
    subroutine s_compute_supersonic_outflow_L(dflt_int) ! --------------------------

        integer, intent(IN) :: dflt_int

        integer :: i !< Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, mom_idx%beg
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = mom_idx%beg + 1, mom_idx%end
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - cont_idx%end)))
        end do

        do i = E_idx, adv_idx%end - 1
            L(i) = lambda(2)*(dadv_ds(i - mom_idx%end))
        end do

        L(adv_idx%end) = lambda(3)*(dpres_ds + rho*c*dvel_ds(dir_idx(1)))

    end subroutine s_compute_supersonic_outflow_L ! ------------------------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are required for the setup of the
        !!      selected CBC.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_initialize_cbc(q_prim_vf, flux_vf, flux_src_vf, & ! ------
                                cbc_dir, cbc_loc, &
                                ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: flux_vf, flux_src_vf

        integer, intent(IN) :: cbc_dir, cbc_loc
        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: dj !< Indical shift based on CBC location

        integer :: i, j, k, r !< Generic loop iterators

        ! Configuring the coordinate direction indexes and flags
        if (cbc_dir == 1) then
            is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1d0, 0d0, 0d0/)
        elseif (cbc_dir == 2) then
            is1%beg = 0; is1%end = buff_size; is2 = ix; is3 = iz
            dir_idx = (/2, 1, 3/); dir_flg = (/0d0, 1d0, 0d0/)
        else
            is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = ix
            dir_idx = (/3, 1, 2/); dir_flg = (/0d0, 0d0, 1d0/)
        end if

        ! Determining the indicial shift based on CBC location
        dj = max(0, cbc_loc)

        ! Allocation/Association of Primitive and Flux Variables ===========
        do i = 1, sys_size
            allocate (q_prim_rs_vf(i)%sf(0:buff_size, &
                                         is2%beg:is2%end, &
                                         is3%beg:is3%end))
        end do

        if (weno_order > 1) then

            do i = 1, adv_idx%end
                allocate (F_rs_vf(i)%sf(0:buff_size, &
                                        is2%beg:is2%end, &
                                        is3%beg:is3%end))
            end do

            allocate (F_src_rs_vf(adv_idx%beg)%sf(0:buff_size, &
                                                  is2%beg:is2%end, &
                                                  is3%beg:is3%end))

            if (riemann_solver == 1) then
                do i = adv_idx%beg + 1, adv_idx%end
                    allocate (F_src_rs_vf(i)%sf(0:buff_size, &
                                                is2%beg:is2%end, &
                                                is3%beg:is3%end))
                end do
            else
                do i = adv_idx%beg + 1, adv_idx%end
                    F_src_rs_vf(i)%sf => F_src_rs_vf(adv_idx%beg)%sf
                end do
            end if

        end if

        do i = 1, adv_idx%end
            allocate (flux_rs_vf(i)%sf(-1:buff_size, &
                                       is2%beg:is2%end, &
                                       is3%beg:is3%end))
        end do

        allocate (flux_src_rs_vf(adv_idx%beg)%sf(-1:buff_size, &
                                                 is2%beg:is2%end, &
                                                 is3%beg:is3%end))

        if (riemann_solver == 1) then
            do i = adv_idx%beg + 1, adv_idx%end
                allocate (flux_src_rs_vf(i)%sf(-1:buff_size, &
                                               is2%beg:is2%end, &
                                               is3%beg:is3%end))
            end do
        else
            do i = adv_idx%beg + 1, adv_idx%end
                flux_src_rs_vf(i)%sf => flux_src_rs_vf(adv_idx%beg)%sf
            end do
        end if
        ! END: Allocation/Association of Primitive and Flux Variables ======

        ! Reshaping Inputted Data in x-direction ===========================
        if (cbc_dir == 1) then

            do i = 1, sys_size
                do r = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = 0, buff_size
                            q_prim_rs_vf(i)%sf(j, k, r) = &
                                q_prim_vf(i)%sf(dj*(m - 2*j) + j, k, r)
                        end do
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = 0, buff_size
                        q_prim_rs_vf(mom_idx%beg)%sf(j, k, r) = &
                            q_prim_vf(mom_idx%beg)%sf(dj*(m - 2*j) + j, k, r)* &
                            sign(1d0, -real(cbc_loc, kind(0d0)))
                    end do
                end do
            end do

            do i = 1, adv_idx%end
                do r = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = -1, buff_size
                            flux_rs_vf(i)%sf(j, k, r) = &
                                flux_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = -1, buff_size
                        flux_rs_vf(mom_idx%beg)%sf(j, k, r) = &
                            flux_vf(mom_idx%beg)%sf(dj*((m - 1) - 2*j) + j, k, r)
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = -1, buff_size
                        flux_src_rs_vf(adv_idx%beg)%sf(j, k, r) = &
                            flux_src_vf(adv_idx%beg)%sf(dj*((m - 1) - 2*j) + j, k, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                do i = adv_idx%beg + 1, adv_idx%end
                    do r = iz%beg, iz%end
                        do k = iy%beg, iy%end
                            do j = -1, buff_size
                                flux_src_rs_vf(i)%sf(j, k, r) = &
                                    flux_src_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r)
                            end do
                        end do
                    end do
                end do
            else
                do r = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = -1, buff_size
                            flux_src_rs_vf(adv_idx%beg)%sf(j, k, r) = &
                                flux_src_vf(adv_idx%beg)%sf(dj*((m - 1) - 2*j) + j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Inputted Data in x-direction ======================

            ! Reshaping Inputted Data in y-direction ===========================
        elseif (cbc_dir == 2) then

            do i = 1, sys_size
                do r = iz%beg, iz%end
                    do k = ix%beg, ix%end
                        do j = 0, buff_size
                            q_prim_rs_vf(i)%sf(j, k, r) = &
                                q_prim_vf(i)%sf(k, dj*(n - 2*j) + j, r)
                        end do
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = ix%beg, ix%end
                    do j = 0, buff_size
                        q_prim_rs_vf(mom_idx%beg + 1)%sf(j, k, r) = &
                            q_prim_vf(mom_idx%beg + 1)%sf(k, dj*(n - 2*j) + j, r)* &
                            sign(1d0, -real(cbc_loc, kind(0d0)))
                    end do
                end do
            end do

            do i = 1, adv_idx%end
                do r = iz%beg, iz%end
                    do k = ix%beg, ix%end
                        do j = -1, buff_size
                            flux_rs_vf(i)%sf(j, k, r) = &
                                flux_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = ix%beg, ix%end
                    do j = -1, buff_size
                        flux_rs_vf(mom_idx%beg + 1)%sf(j, k, r) = &
                            flux_vf(mom_idx%beg + 1)%sf(k, dj*((n - 1) - 2*j) + j, r)
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = ix%beg, ix%end
                    do j = -1, buff_size
                        flux_src_rs_vf(adv_idx%beg)%sf(j, k, r) = &
                            flux_src_vf(adv_idx%beg)%sf(k, dj*((n - 1) - 2*j) + j, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                do i = adv_idx%beg + 1, adv_idx%end
                    do r = iz%beg, iz%end
                        do k = ix%beg, ix%end
                            do j = -1, buff_size
                                flux_src_rs_vf(i)%sf(j, k, r) = &
                                    flux_src_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r)
                            end do
                        end do
                    end do
                end do
            else
                do r = iz%beg, iz%end
                    do k = ix%beg, ix%end
                        do j = -1, buff_size
                            flux_src_rs_vf(adv_idx%beg)%sf(j, k, r) = &
                                flux_src_vf(adv_idx%beg)%sf(k, dj*((n - 1) - 2*j) + j, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Inputted Data in y-direction ======================

            ! Reshaping Inputted Data in z-direction ===========================
        else

            do i = 1, sys_size
                do r = ix%beg, ix%end
                    do k = iy%beg, iy%end
                        do j = 0, buff_size
                            q_prim_rs_vf(i)%sf(j, k, r) = &
                                q_prim_vf(i)%sf(r, k, dj*(p - 2*j) + j)
                        end do
                    end do
                end do
            end do

            do r = ix%beg, ix%end
                do k = iy%beg, iy%end
                    do j = 0, buff_size
                        q_prim_rs_vf(mom_idx%end)%sf(j, k, r) = &
                            q_prim_vf(mom_idx%end)%sf(r, k, dj*(p - 2*j) + j)* &
                            sign(1d0, -real(cbc_loc, kind(0d0)))
                    end do
                end do
            end do

            do i = 1, adv_idx%end
                do r = ix%beg, ix%end
                    do k = iy%beg, iy%end
                        do j = -1, buff_size
                            flux_rs_vf(i)%sf(j, k, r) = &
                                flux_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            do r = ix%beg, ix%end
                do k = iy%beg, iy%end
                    do j = -1, buff_size
                        flux_rs_vf(mom_idx%end)%sf(j, k, r) = &
                            flux_vf(mom_idx%end)%sf(r, k, dj*((p - 1) - 2*j) + j)
                    end do
                end do
            end do

            do r = ix%beg, ix%end
                do k = iy%beg, iy%end
                    do j = -1, buff_size
                        flux_src_rs_vf(adv_idx%beg)%sf(j, k, r) = &
                            flux_src_vf(adv_idx%beg)%sf(r, k, dj*((p - 1) - 2*j) + j)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                do i = adv_idx%beg + 1, adv_idx%end
                    do r = ix%beg, ix%end
                        do k = iy%beg, iy%end
                            do j = -1, buff_size
                                flux_src_rs_vf(i)%sf(j, k, r) = &
                                    flux_src_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j)
                            end do
                        end do
                    end do
                end do
            else
                do r = ix%beg, ix%end
                    do k = iy%beg, iy%end
                        do j = -1, buff_size
                            flux_src_rs_vf(adv_idx%beg)%sf(j, k, r) = &
                                flux_src_vf(adv_idx%beg)%sf(r, k, dj*((p - 1) - 2*j) + j)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if

        end if
        ! END: Reshaping Inputted Data in z-direction ======================

        ! Association of the procedural pointer to the appropriate procedure
        ! that will be utilized in the evaluation of L variables for the CBC

        ! ==================================================================
        if ((cbc_dir == 1 .and. cbc_loc == -1 .and. bc_x%beg == -5) &
            .or. (cbc_dir == 1 .and. cbc_loc == 1 .and. bc_x%end == -5) &
            .or. (cbc_dir == 2 .and. cbc_loc == -1 .and. bc_y%beg == -5) &
            .or. (cbc_dir == 2 .and. cbc_loc == 1 .and. bc_y%end == -5) &
            .or. (cbc_dir == 3 .and. cbc_loc == -1 .and. bc_z%beg == -5) &
            .or. (cbc_dir == 3 .and. cbc_loc == 1 .and. bc_z%end == -5)) &
            then

            s_compute_L => s_compute_slip_wall_L

        elseif ((cbc_dir == 1 .and. cbc_loc == -1 .and. bc_x%beg == -6) &
                .or. (cbc_dir == 1 .and. cbc_loc == 1 .and. bc_x%end == -6) &
                .or. (cbc_dir == 2 .and. cbc_loc == -1 .and. bc_y%beg == -6) &
                .or. (cbc_dir == 2 .and. cbc_loc == 1 .and. bc_y%end == -6) &
                .or. (cbc_dir == 3 .and. cbc_loc == -1 .and. bc_z%beg == -6) &
                .or. (cbc_dir == 3 .and. cbc_loc == 1 .and. bc_z%end == -6)) &
            then

            s_compute_L => s_compute_nonreflecting_subsonic_buffer_L

        elseif ((cbc_dir == 1 .and. cbc_loc == -1 .and. bc_x%beg == -7) &
                .or. (cbc_dir == 1 .and. cbc_loc == 1 .and. bc_x%end == -7) &
                .or. (cbc_dir == 2 .and. cbc_loc == -1 .and. bc_y%beg == -7) &
                .or. (cbc_dir == 2 .and. cbc_loc == 1 .and. bc_y%end == -7) &
                .or. (cbc_dir == 3 .and. cbc_loc == -1 .and. bc_z%beg == -7) &
                .or. (cbc_dir == 3 .and. cbc_loc == 1 .and. bc_z%end == -7)) &
            then

            s_compute_L => s_compute_nonreflecting_subsonic_inflow_L

        elseif ((cbc_dir == 1 .and. cbc_loc == -1 .and. bc_x%beg == -8) &
                .or. (cbc_dir == 1 .and. cbc_loc == 1 .and. bc_x%end == -8) &
                .or. (cbc_dir == 2 .and. cbc_loc == -1 .and. bc_y%beg == -8) &
                .or. (cbc_dir == 2 .and. cbc_loc == 1 .and. bc_y%end == -8) &
                .or. (cbc_dir == 3 .and. cbc_loc == -1 .and. bc_z%beg == -8) &
                .or. (cbc_dir == 3 .and. cbc_loc == 1 .and. bc_z%end == -8)) &
            then

            s_compute_L => s_compute_nonreflecting_subsonic_outflow_L

        elseif ((cbc_dir == 1 .and. cbc_loc == -1 .and. bc_x%beg == -9) &
                .or. (cbc_dir == 1 .and. cbc_loc == 1 .and. bc_x%end == -9) &
                .or. (cbc_dir == 2 .and. cbc_loc == -1 .and. bc_y%beg == -9) &
                .or. (cbc_dir == 2 .and. cbc_loc == 1 .and. bc_y%end == -9) &
                .or. (cbc_dir == 3 .and. cbc_loc == -1 .and. bc_z%beg == -9) &
                .or. (cbc_dir == 3 .and. cbc_loc == 1 .and. bc_z%end == -9)) &
            then

            s_compute_L => s_compute_force_free_subsonic_outflow_L

        elseif ((cbc_dir == 1 .and. cbc_loc == -1 .and. bc_x%beg == -10) &
                .or. (cbc_dir == 1 .and. cbc_loc == 1 .and. bc_x%end == -10) &
                .or. (cbc_dir == 2 .and. cbc_loc == -1 .and. bc_y%beg == -10) &
                .or. (cbc_dir == 2 .and. cbc_loc == 1 .and. bc_y%end == -10) &
                .or. (cbc_dir == 3 .and. cbc_loc == -1 .and. bc_z%beg == -10) &
                .or. (cbc_dir == 3 .and. cbc_loc == 1 .and. bc_z%end == -10)) &
            then

            s_compute_L => s_compute_constant_pressure_subsonic_outflow_L

        elseif ((cbc_dir == 1 .and. cbc_loc == -1 .and. bc_x%beg == -11) &
                .or. (cbc_dir == 1 .and. cbc_loc == 1 .and. bc_x%end == -11) &
                .or. (cbc_dir == 2 .and. cbc_loc == -1 .and. bc_y%beg == -11) &
                .or. (cbc_dir == 2 .and. cbc_loc == 1 .and. bc_y%end == -11) &
                .or. (cbc_dir == 3 .and. cbc_loc == -1 .and. bc_z%beg == -11) &
                .or. (cbc_dir == 3 .and. cbc_loc == 1 .and. bc_z%end == -11)) &
            then

            s_compute_L => s_compute_supersonic_inflow_L

        else

            s_compute_L => s_compute_supersonic_outflow_L

        end if
        ! ==================================================================

    end subroutine s_initialize_cbc ! --------------------------------------

    !>  Deallocation and/or the disassociation procedures that
        !!      are necessary in order to finalize the CBC application
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_finalize_cbc(flux_vf, flux_src_vf, & ! -------------------
                              cbc_dir, cbc_loc, &
                              ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf

        integer, intent(IN) :: cbc_dir, cbc_loc
        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: dj !< Indical shift based on CBC location

        integer :: i, j, k, r !< Generic loop iterators

        ! Determining the indicial shift based on CBC location
        dj = max(0, cbc_loc)

        ! Reshaping Outputted Data in x-direction ==========================
        if (cbc_dir == 1) then

            do i = 1, adv_idx%end
                do r = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                flux_rs_vf(i)%sf(j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = -1, buff_size
                        flux_vf(mom_idx%beg)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                            flux_rs_vf(mom_idx%beg)%sf(j, k, r)
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = -1, buff_size
                        flux_src_vf(adv_idx%beg)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                            flux_src_rs_vf(adv_idx%beg)%sf(j, k, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                do i = adv_idx%beg + 1, adv_idx%end
                    do r = iz%beg, iz%end
                        do k = iy%beg, iy%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                    flux_src_rs_vf(i)%sf(j, k, r)
                            end do
                        end do
                    end do
                end do
            else
                do r = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = -1, buff_size
                            flux_src_vf(adv_idx%beg)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                flux_src_rs_vf(adv_idx%beg)%sf(j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Outputted Data in x-direction =====================

            ! Reshaping Outputted Data in y-direction ==========================
        elseif (cbc_dir == 2) then

            do i = 1, adv_idx%end
                do r = iz%beg, iz%end
                    do k = ix%beg, ix%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                flux_rs_vf(i)%sf(j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = ix%beg, ix%end
                    do j = -1, buff_size
                        flux_vf(mom_idx%beg + 1)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                            flux_rs_vf(mom_idx%beg + 1)%sf(j, k, r)
                    end do
                end do
            end do

            do r = iz%beg, iz%end
                do k = ix%beg, ix%end
                    do j = -1, buff_size
                        flux_src_vf(adv_idx%beg)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                            flux_src_rs_vf(adv_idx%beg)%sf(j, k, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                do i = adv_idx%beg + 1, adv_idx%end
                    do r = iz%beg, iz%end
                        do k = ix%beg, ix%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                    flux_src_rs_vf(i)%sf(j, k, r)
                            end do
                        end do
                    end do
                end do
            else
                do r = iz%beg, iz%end
                    do k = ix%beg, ix%end
                        do j = -1, buff_size
                            flux_src_vf(adv_idx%beg)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                flux_src_rs_vf(adv_idx%beg)%sf(j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Outputted Data in y-direction =====================

            ! Reshaping Outputted Data in z-direction ==========================
        else

            do i = 1, adv_idx%end
                do r = ix%beg, ix%end
                    do k = iy%beg, iy%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                flux_rs_vf(i)%sf(j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            do r = ix%beg, ix%end
                do k = iy%beg, iy%end
                    do j = -1, buff_size
                        flux_vf(mom_idx%end)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                            flux_rs_vf(mom_idx%end)%sf(j, k, r)
                    end do
                end do
            end do

            do r = ix%beg, ix%end
                do k = iy%beg, iy%end
                    do j = -1, buff_size
                        flux_src_vf(adv_idx%beg)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                            flux_src_rs_vf(adv_idx%beg)%sf(j, k, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                do i = adv_idx%beg + 1, adv_idx%end
                    do r = ix%beg, ix%end
                        do k = iy%beg, iy%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                    flux_src_rs_vf(i)%sf(j, k, r)
                            end do
                        end do
                    end do
                end do
            else
                do r = ix%beg, ix%end
                    do k = iy%beg, iy%end
                        do j = -1, buff_size
                            flux_src_vf(adv_idx%beg)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                flux_src_rs_vf(adv_idx%beg)%sf(j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if

        end if
        ! END: Reshaping Outputted Data in z-direction =====================

        ! Deallocation/Disassociation of Primitive and Flux Variables ======
        do i = 1, sys_size
            deallocate (q_prim_rs_vf(i)%sf)
        end do

        if (weno_order > 1) then

            do i = 1, adv_idx%end
                deallocate (F_rs_vf(i)%sf)
            end do

            deallocate (F_src_rs_vf(adv_idx%beg)%sf)

            if (riemann_solver == 1) then
                do i = adv_idx%beg + 1, adv_idx%end
                    deallocate (F_src_rs_vf(i)%sf)
                end do
            else
                do i = adv_idx%beg + 1, adv_idx%end
                    nullify (F_src_rs_vf(i)%sf)
                end do
            end if

        end if

        do i = 1, adv_idx%end
            deallocate (flux_rs_vf(i)%sf)
        end do

        deallocate (flux_src_rs_vf(adv_idx%beg)%sf)

        if (riemann_solver == 1) then
            do i = adv_idx%beg + 1, adv_idx%end
                deallocate (flux_src_rs_vf(i)%sf)
            end do
        else
            do i = adv_idx%beg + 1, adv_idx%end
                nullify (flux_src_rs_vf(i)%sf)
            end do
        end if
        ! ==================================================================

        ! Nullifying procedural pointer used in evaluation of L for the CBC
        s_compute_L => null()

    end subroutine s_finalize_cbc ! ----------------------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_cbc_module() ! -----------------------------------

        if (all((/bc_x%beg, bc_x%end/) > -5) &
            .and. &
            (n > 0 .and. all((/bc_y%beg, bc_y%end/) > -5)) &
            .and. &
            (p > 0 .and. all((/bc_z%beg, bc_z%end/) > -5))) return

        ! Deallocating the cell-average primitive variables
        deallocate (q_prim_rs_vf)

        ! Deallocating the cell-average and cell-boundary-average fluxes
        deallocate (F_rs_vf, F_src_rs_vf)
        deallocate (flux_rs_vf, flux_src_rs_vf)

        ! Deallocating the cell-average partial densities, the velocity, the
        ! advection variables, the mass fractions and also the Weber numbers
        deallocate (alpha_rho, vel, adv, mf, We)

        ! Deallocating the first-order spatial derivatives, in s-direction,
        ! of the partial densities, the velocity and the advected variables
        deallocate (dalpha_rho_ds, dvel_ds, dadv_ds)

        ! Deallocating L, see Thompson (1987, 1990)
        deallocate (L)

        ! Deallocating the cell-width distribution in the s-direction
        deallocate (ds)

        ! Deallocating CBC Coefficients in x-direction =====================
        if (any((/bc_x%beg, bc_x%end/) <= -5)) then
            deallocate (fd_coef_x); if (weno_order > 1) deallocate (pi_coef_x)
        end if
        ! ==================================================================

        ! Deallocating CBC Coefficients in y-direction =====================
        if (n > 0 .and. any((/bc_y%beg, bc_y%end/) <= -5)) then
            deallocate (fd_coef_y); if (weno_order > 1) deallocate (pi_coef_y)
        end if
        ! ==================================================================

        ! Deallocating CBC Coefficients in z-direction =====================
        if (p > 0 .and. any((/bc_z%beg, bc_z%end/) <= -5)) then
            deallocate (fd_coef_z); if (weno_order > 1) deallocate (pi_coef_z)
        end if
        ! ==================================================================

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_cbc_module ! ---------------------------------

end module m_cbc
