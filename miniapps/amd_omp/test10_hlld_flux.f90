! test10_hlld_flux.f90
!
! Tests the HLLD F flux computation using multiple PRIVATE DT variables
! (vec3 for B and vel, arr7 for U and F, scalar_lr for rho/pres/E/pTot)
! inside an OpenMP target offload parallel loop.
!
! Uses Dai-Woodward MHD test values scaled by cell index.
! Checks that all 7 flux components are computed correctly on device.
program test10_hlld_flux
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 10000

    type :: vec3
        real(wp) :: L(3), R(3)
    end type
    type :: arr7
        real(wp) :: L(7), R(7)
    end type
    type :: scalar_lr
        real(wp) :: L, R
    end type

    ! Dai-Woodward left state (dimensionless)
    real(wp), parameter :: rho0 = 1.08_wp
    real(wp), parameter :: vx0 = 1.2_wp, vy0 = 0.01_wp, vz0 = 0.5_wp
    real(wp), parameter :: Bx0 = 0.5641895835_wp
    real(wp), parameter :: By0 = 1.0149412604_wp
    real(wp), parameter :: Bz0 = 0.5641895835_wp
    real(wp), parameter :: pres0 = 0.95_wp
    real(wp), parameter :: gam  = 5.0_wp/3.0_wp

    real(wp) :: F_out(N, 7), F_expected(N, 7)
    type(vec3) :: vel_priv, B_priv
    type(arr7) :: U_priv, F_priv
    type(scalar_lr) :: rho_priv, pm_priv, E_priv, pT_priv
    real(wp) :: rho_in(N), vx_in(N), vy_in(N), vz_in(N)
    real(wp) :: Bx_in(N), By_in(N), Bz_in(N), pres_in(N)
    integer :: i, j, nerr
    real(wp) :: pm, E, pTot

    do i = 1, N
        rho_in(i)  = rho0 * (1._wp + 0.001_wp*real(i, wp))
        vx_in(i)   = vx0
        vy_in(i)   = vy0 * real(i, wp)
        vz_in(i)   = vz0
        Bx_in(i)   = Bx0
        By_in(i)   = By0 * (1._wp + 0.001_wp*real(i, wp))
        Bz_in(i)   = Bz0
        pres_in(i) = pres0 * (1._wp + 0.001_wp*real(i, wp))

        pm   = 0.5_wp*(Bx_in(i)**2 + By_in(i)**2 + Bz_in(i)**2)
        E    = pres_in(i)/(gam - 1._wp) + 0.5_wp*rho_in(i)*(vx_in(i)**2 + vy_in(i)**2 + vz_in(i)**2) + pm
        pTot = pres_in(i) + pm

        F_expected(i, 1) = rho_in(i)*vx_in(i)
        F_expected(i, 2) = rho_in(i)*vx_in(i)**2 - Bx_in(i)**2 + pTot
        F_expected(i, 3) = rho_in(i)*vx_in(i)*vy_in(i) - Bx_in(i)*By_in(i)
        F_expected(i, 4) = rho_in(i)*vx_in(i)*vz_in(i) - Bx_in(i)*Bz_in(i)
        F_expected(i, 5) = vx_in(i)*By_in(i) - vy_in(i)*Bx_in(i)
        F_expected(i, 6) = vx_in(i)*Bz_in(i) - vz_in(i)*Bx_in(i)
        F_expected(i, 7) = (E + pTot)*vx_in(i) - Bx_in(i)*(vx_in(i)*Bx_in(i) + vy_in(i)*By_in(i) + vz_in(i)*Bz_in(i))
    end do
    F_out = 0._wp

    !$omp target teams distribute parallel do &
    !$omp   map(to:rho_in,vx_in,vy_in,vz_in,Bx_in,By_in,Bz_in,pres_in) &
    !$omp   map(from:F_out) &
    !$omp   private(vel_priv, B_priv, U_priv, F_priv, rho_priv, pm_priv, E_priv, pT_priv, j)
    do i = 1, N
        rho_priv%L = rho_in(i)
        vel_priv%L(1) = vx_in(i);  vel_priv%L(2) = vy_in(i);  vel_priv%L(3) = vz_in(i)
        B_priv%L(1)   = Bx_in(i);  B_priv%L(2)   = By_in(i);  B_priv%L(3)   = Bz_in(i)

        pm_priv%L = 0.5_wp*(B_priv%L(1)**2 + B_priv%L(2)**2 + B_priv%L(3)**2)
        E_priv%L  = pres_in(i)/(gam - 1._wp) + &
                    0.5_wp*rho_priv%L*(vel_priv%L(1)**2 + vel_priv%L(2)**2 + vel_priv%L(3)**2) + &
                    pm_priv%L
        pT_priv%L = pres_in(i) + pm_priv%L

        U_priv%L(1) = rho_priv%L
        U_priv%L(2) = rho_priv%L*vel_priv%L(1)
        U_priv%L(3) = rho_priv%L*vel_priv%L(2)
        U_priv%L(4) = rho_priv%L*vel_priv%L(3)
        U_priv%L(5) = B_priv%L(2)
        U_priv%L(6) = B_priv%L(3)
        U_priv%L(7) = E_priv%L

        F_priv%L(1) = U_priv%L(2)
        F_priv%L(2) = U_priv%L(2)*vel_priv%L(1) - B_priv%L(1)*B_priv%L(1) + pT_priv%L
        F_priv%L(3) = U_priv%L(2)*vel_priv%L(2) - B_priv%L(1)*B_priv%L(2)
        F_priv%L(4) = U_priv%L(2)*vel_priv%L(3) - B_priv%L(1)*B_priv%L(3)
        F_priv%L(5) = vel_priv%L(1)*B_priv%L(2) - vel_priv%L(2)*B_priv%L(1)
        F_priv%L(6) = vel_priv%L(1)*B_priv%L(3) - vel_priv%L(3)*B_priv%L(1)
        F_priv%L(7) = (E_priv%L + pT_priv%L)*vel_priv%L(1) - &
                      B_priv%L(1)*(vel_priv%L(1)*B_priv%L(1) + vel_priv%L(2)*B_priv%L(2) + vel_priv%L(3)*B_priv%L(3))

        do j = 1, 7
            F_out(i, j) = F_priv%L(j)
        end do
    end do
    !$omp end target teams distribute parallel do

    nerr = 0
    do i = 1, N
        do j = 1, 7
            if (abs(F_out(i, j) - F_expected(i, j)) > 1.e-8_wp*max(abs(F_expected(i,j)), 1._wp)) nerr = nerr + 1
        end do
    end do

    if (nerr == 0) then
        print *, "PASS test10: HLLD F flux computed correctly with private DT vars"
    else
        print *, "FAIL test10:", nerr, "errors -- HLLD flux DT computation"
    end if
end program test10_hlld_flux
