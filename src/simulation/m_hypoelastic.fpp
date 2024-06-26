!>
!! @file m_hypoelastic.f90
!! @brief Contains module m_hypoelastic

#:include 'macros.fpp'

!> @brief This module is used to compute source terms for hypoelastic model
module m_hypoelastic

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_hypoelastic_module, &
 s_compute_hypoelastic_rhs

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), Gs)
    !$acc declare link(Gs)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), du_dx, du_dy, du_dz)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), dv_dx, dv_dy, dv_dz)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), dw_dx, dw_dy, dw_dz)
    !$acc declare link(du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), rho_K_field, G_K_field)
    !$acc declare link(rho_K_field, G_K_field)
#else
    real(kind(0d0)), allocatable, dimension(:) :: Gs
    !$acc declare create(Gs)

    real(kind(0d0)), allocatable, dimension(:, :, :) :: du_dx, du_dy, du_dz
    real(kind(0d0)), allocatable, dimension(:, :, :) :: dv_dx, dv_dy, dv_dz
    real(kind(0d0)), allocatable, dimension(:, :, :) :: dw_dx, dw_dy, dw_dz
    !$acc declare create(du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz)

    real(kind(0d0)), allocatable, dimension(:, :, :) :: rho_K_field, G_K_field
    !$acc declare create(rho_K_field, G_K_field)

#endif

contains

    subroutine s_initialize_hypoelastic_module

        integer :: i

        @:ALLOCATE_GLOBAL(Gs(1:num_fluids))
        @:ALLOCATE_GLOBAL(rho_K_field(0:m,0:n,0:p), G_K_field(0:m,0:n,0:p))
        @:ALLOCATE_GLOBAL(du_dx(0:m,0:n,0:p))
        if (n > 0) then
            @:ALLOCATE_GLOBAL(du_dy(0:m,0:n,0:p), dv_dx(0:m,0:n,0:p), dv_dy(0:m,0:n,0:p))
            if (p > 0) then
                @:ALLOCATE_GLOBAL(du_dz(0:m,0:n,0:p), dv_dz(0:m,0:n,0:p))
                @:ALLOCATE_GLOBAL(dw_dx(0:m,0:n,0:p), dw_dy(0:m,0:n,0:p), dw_dz(0:m,0:n,0:p))
            end if
        end if

        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        !$acc update device(Gs)

    end subroutine s_initialize_hypoelastic_module

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the elastic stress equations
        !!  @param idir Dimension splitting index
        !!  @param q_prim_vf Primitive variables
        !!  @param rhs_vf rhs variables
    subroutine s_compute_hypoelastic_rhs(idir, q_prim_vf, rhs_vf)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(kind(0d0)) :: rho_K, G_K

        integer :: i, k, l, q !< Loop variables
        integer :: ndirs  !< Number of coordinate directions

        ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3

        if (idir == 1) then
            ! calculate velocity gradients + rho_K and G_K
            ! TODO: re-organize these loops one by one for GPU efficiency if possible?

            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        du_dx(k, l, q) = &
                            (q_prim_vf(momxb)%sf(k - 2, l, q) &
                             - 8d0*q_prim_vf(momxb)%sf(k - 1, l, q) &
                             + 8d0*q_prim_vf(momxb)%sf(k + 1, l, q) &
                             - q_prim_vf(momxb)%sf(k + 2, l, q)) &
                            /(12d0*dx(k))
                    end do
                end do
            end do

            if (ndirs > 1) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            du_dy(k, l, q) = &
                                (q_prim_vf(momxb)%sf(k, l - 2, q) &
                                 - 8d0*q_prim_vf(momxb)%sf(k, l - 1, q) &
                                 + 8d0*q_prim_vf(momxb)%sf(k, l + 1, q) &
                                 - q_prim_vf(momxb)%sf(k, l + 2, q)) &
                                /(12d0*dy(l))
                            dv_dx(k, l, q) = &
                                (q_prim_vf(momxb + 1)%sf(k - 2, l, q) &
                                 - 8d0*q_prim_vf(momxb + 1)%sf(k - 1, l, q) &
                                 + 8d0*q_prim_vf(momxb + 1)%sf(k + 1, l, q) &
                                 - q_prim_vf(momxb + 1)%sf(k + 2, l, q)) &
                                /(12d0*dx(k))
                            dv_dy(k, l, q) = &
                                (q_prim_vf(momxb + 1)%sf(k, l - 2, q) &
                                 - 8d0*q_prim_vf(momxb + 1)%sf(k, l - 1, q) &
                                 + 8d0*q_prim_vf(momxb + 1)%sf(k, l + 1, q) &
                                 - q_prim_vf(momxb + 1)%sf(k, l + 2, q)) &
                                /(12d0*dy(l))
                        end do
                    end do
                end do

                ! 3D
                if (ndirs == 3) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                du_dz(k, l, q) = &
                                    (q_prim_vf(momxb)%sf(k, l, q - 2) &
                                     - 8d0*q_prim_vf(momxb)%sf(k, l, q - 1) &
                                     + 8d0*q_prim_vf(momxb)%sf(k, l, q + 1) &
                                     - q_prim_vf(momxb)%sf(k, l, q + 2)) &
                                    /(12d0*dz(q))
                                dv_dz(k, l, q) = &
                                    (q_prim_vf(momxb + 1)%sf(k, l, q - 2) &
                                     - 8d0*q_prim_vf(momxb + 1)%sf(k, l, q - 1) &
                                     + 8d0*q_prim_vf(momxb + 1)%sf(k, l, q + 1) &
                                     - q_prim_vf(momxb + 1)%sf(k, l, q + 2)) &
                                    /(12d0*dz(q))
                                dw_dx(k, l, q) = &
                                    (q_prim_vf(momxe)%sf(k - 2, l, q) &
                                     - 8d0*q_prim_vf(momxe)%sf(k - 1, l, q) &
                                     + 8d0*q_prim_vf(momxe)%sf(k + 1, l, q) &
                                     - q_prim_vf(momxe)%sf(k + 2, l, q)) &
                                    /(12d0*dx(k))
                                dw_dy(k, l, q) = &
                                    (q_prim_vf(momxe)%sf(k, l - 2, q) &
                                     - 8d0*q_prim_vf(momxe)%sf(k, l - 1, q) &
                                     + 8d0*q_prim_vf(momxe)%sf(k, l + 1, q) &
                                     - q_prim_vf(momxe)%sf(k, l + 2, q)) &
                                    /(12d0*dy(l))
                                dw_dz(k, l, q) = &
                                    (q_prim_vf(momxe)%sf(k, l, q - 2) &
                                     - 8d0*q_prim_vf(momxe)%sf(k, l, q - 1) &
                                     + 8d0*q_prim_vf(momxe)%sf(k, l, q + 1) &
                                     - q_prim_vf(momxe)%sf(k, l, q + 2)) &
                                    /(12d0*dz(q))
                            end do
                        end do
                    end do
                end if
            end if

            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rho_K = 0d0; G_K = 0d0
                        do i = 1, num_fluids
                            rho_K = rho_K + q_prim_vf(i)%sf(k, l, q) !alpha_rho_K(1)
                            G_K = G_K + q_prim_vf(advxb - 1 + i)%sf(k, l, q)*Gs(i)  !alpha_K(1) * Gs(1)
                        end do
                        rho_K_field(k, l, q) = rho_K
                        G_K_field(k, l, q) = G_K

                        !TODO: take this out if not needed
                        if (G_K < 1000) then
                            G_K_field(k, l, q) = 0
                        end if
                    end do
                end do
            end do

            ! apply rhs source term to elastic stress equation
            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(strxb)%sf(k, l, q) = &
                            rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l, q)* &
                            ((4d0*G_K_field(k, l, q)/3d0) + &
                             q_prim_vf(strxb)%sf(k, l, q))* &
                            du_dx(k, l, q)
                    end do
                end do
            end do

        elseif (idir == 2) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(strxb)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                    (q_prim_vf(strxb + 1)%sf(k, l, q)*du_dy(k, l, q) + &
                                                     q_prim_vf(strxb + 1)%sf(k, l, q)*du_dy(k, l, q) - &
                                                     q_prim_vf(strxb)%sf(k, l, q)*dv_dy(k, l, q) - &
                                                     2d0*G_K_field(k, l, q)*(1d0/3d0)*dv_dy(k, l, q))

                        rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 1)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb)%sf(k, l, q)*dv_dx(k, l, q) - &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*du_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dy(k, l, q) - &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         2d0*G_K_field(k, l, q)*(1d0/2d0)*(du_dy(k, l, q) + &
                                                                                           dv_dx(k, l, q)))

                        rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dx(k, l, q) - &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dv_dy(k, l, q) - &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         2d0*G_K_field(k, l, q)*(dv_dy(k, l, q) - (1d0/3d0)* &
                                                                                 (du_dx(k, l, q) + &
                                                                                  dv_dy(k, l, q))))
                    end do
                end do
            end do

        elseif (idir == 3) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(strxb)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                    (q_prim_vf(strxb + 3)%sf(k, l, q)*du_dz(k, l, q) + &
                                                     q_prim_vf(strxb + 3)%sf(k, l, q)*du_dz(k, l, q) - &
                                                     q_prim_vf(strxb)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                     2d0*G_K_field(k, l, q)*(1d0/3d0)*dw_dz(k, l, q))

                        rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 4)%sf(k, l, q)*du_dz(k, l, q) + &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dv_dz(k, l, q) - &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dw_dz(k, l, q))

                        rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dz(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dz(k, l, q) - &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                         2d0*G_K_field(k, l, q)*(1d0/3d0)*dw_dz(k, l, q))

                        rhs_vf(strxb + 3)%sf(k, l, q) = rhs_vf(strxb + 3)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 3)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb)%sf(k, l, q)*dw_dx(k, l, q) - &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*du_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dw_dy(k, l, q) - &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 5)%sf(k, l, q)*du_dz(k, l, q) + &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dw_dz(k, l, q) + &
                                                         2d0*G_K_field(k, l, q)*(1d0/2d0)*(du_dz(k, l, q) + &
                                                                                           dw_dx(k, l, q)))

                        rhs_vf(strxb + 4)%sf(k, l, q) = rhs_vf(strxb + 4)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 3)%sf(k, l, q)*dv_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dw_dx(k, l, q) - &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dw_dy(k, l, q) - &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 5)%sf(k, l, q)*dv_dz(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dw_dz(k, l, q) + &
                                                         2d0*G_K_field(k, l, q)*(1d0/2d0)*(dv_dz(k, l, q) + &
                                                                                           dw_dy(k, l, q)))

                        rhs_vf(strxe)%sf(k, l, q) = rhs_vf(strxe)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                    (q_prim_vf(strxe - 2)%sf(k, l, q)*dw_dx(k, l, q) + &
                                                     q_prim_vf(strxe - 2)%sf(k, l, q)*dw_dx(k, l, q) - &
                                                     q_prim_vf(strxe)%sf(k, l, q)*du_dx(k, l, q) + &
                                                     q_prim_vf(strxe - 1)%sf(k, l, q)*dw_dy(k, l, q) + &
                                                     q_prim_vf(strxe - 1)%sf(k, l, q)*dw_dy(k, l, q) - &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dw_dz(k, l, q) + &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dw_dz(k, l, q) + &
                                                     2d0*G_K_field(k, l, q)*(dw_dz(k, l, q) - (1d0/3d0)* &
                                                                             (du_dx(k, l, q) + &
                                                                              dv_dy(k, l, q) + &
                                                                              dw_dz(k, l, q))))
                    end do
                end do
            end do
        end if

    end subroutine s_compute_hypoelastic_rhs

end module m_hypoelastic
