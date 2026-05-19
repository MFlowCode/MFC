! test11_hlld_full.f90
!
! Full HLLD Riemann solver for a single interface with hardcoded
! Dai-Woodward left/right states. Computes U, U_star, F, F_star,
! double-star states, selects F_hlld via wave-speed logic.
!
! Reference solution computed on host; GPU result compared.
! All DT variable patterns match the actual MFC HLLD code.

program test11_hlld_full
    implicit none
    integer, parameter :: wp = 8

    type :: vec3;    real(wp) :: L(3), R(3); end type
    type :: arr7;    real(wp) :: L(7), R(7); end type
    type :: slr;     real(wp) :: L, R;       end type

    ! Dai-Woodward initial states (Miyoshi & Kusano 2005, Table 1)
    real(wp), parameter :: rhoL = 1.08_wp,  rhoR = 1.0_wp
    real(wp), parameter :: vxL  = 1.2_wp,   vxR  = 0.0_wp
    real(wp), parameter :: vyL  = 0.01_wp,  vyR  = 0.0_wp
    real(wp), parameter :: vzL  = 0.5_wp,   vzR  = 0.0_wp
    real(wp), parameter :: BxL  = 0.5641895835_wp, BxR = 0.5641895835_wp
    real(wp), parameter :: ByL  = 1.0149412604_wp, ByR = 1.1283791670_wp
    real(wp), parameter :: BzL  = 0.5641895835_wp, BzR = 0.5641895835_wp
    real(wp), parameter :: presL = 0.95_wp, presR = 1.0_wp
    real(wp), parameter :: gam  = 5.0_wp/3.0_wp

    real(wp) :: F_gpu(7), F_ref(7)
    type(vec3) :: vel, B
    type(arr7) :: U, U_star, U_double, F, F_star
    type(slr)  :: rho, pres, E, H, gamma, pi_inf, qv, vel_rms
    type(slr)  :: c, c_fast, pres_mag, pTot, s, rho_star, E_star, sqrt_rho_star
    type(slr)  :: v_star, w_star, E_double_lr, p_star
    real(wp)   :: F_hlld(7)
    real(wp)   :: s_M, s_starL, s_starR, denom_ds, sign_Bx
    real(wp)   :: v_double, w_double, By_double, Bz_double, E_double
    real(wp)   :: B2, term, disc
    integer    :: i

    F_gpu = 0._wp

    !$omp target map(from:F_gpu) &
    !$omp   private(vel, B, U, U_star, U_double, F, F_star, &
    !$omp           rho, pres, E, H, gamma, pi_inf, qv, vel_rms, &
    !$omp           c, c_fast, pres_mag, pTot, s, rho_star, E_star, sqrt_rho_star, &
    !$omp           v_star, w_star, E_double_lr, p_star, &
    !$omp           s_M, s_starL, s_starR, denom_ds, sign_Bx, &
    !$omp           v_double, w_double, By_double, Bz_double, E_double, &
    !$omp           B2, term, disc, F_hlld, i)
        ! Load states
        rho%L = rhoL;   rho%R = rhoR
        pres%L = presL; pres%R = presR
        vel%L(1) = vxL; vel%L(2) = vyL; vel%L(3) = vzL
        vel%R(1) = vxR; vel%R(2) = vyR; vel%R(3) = vzR
        B%L(1) = BxL;   B%L(2) = ByL;   B%L(3) = BzL
        B%R(1) = BxR;   B%R(2) = ByR;   B%R(3) = BzR

        gamma%L = gam;  gamma%R = gam
        pi_inf%L = 0._wp; pi_inf%R = 0._wp
        qv%L = 0._wp;   qv%R = 0._wp

        vel_rms%L = vel%L(1)**2 + vel%L(2)**2 + vel%L(3)**2
        vel_rms%R = vel%R(1)**2 + vel%R(2)**2 + vel%R(3)**2

        pres_mag%L = 0.5_wp*(B%L(1)**2 + B%L(2)**2 + B%L(3)**2)
        pres_mag%R = 0.5_wp*(B%R(1)**2 + B%R(2)**2 + B%R(3)**2)

        E%L = gamma%L*pres%L/(gamma%L - 1._wp) + 0.5_wp*rho%L*vel_rms%L + pres_mag%L
        E%R = gamma%R*pres%R/(gamma%R - 1._wp) + 0.5_wp*rho%R*vel_rms%R + pres_mag%R

        pTot%L = pres%L + pres_mag%L
        pTot%R = pres%R + pres_mag%R

        ! Sound speed (ideal gas, no pi_inf)
        c%L = sqrt(gamma%L*pres%L/rho%L)
        c%R = sqrt(gamma%R*pres%R/rho%R)

        ! Fast magnetosonic speed
        B2 = B%L(1)**2 + B%L(2)**2 + B%L(3)**2
        term = c%L**2 + B2/rho%L
        disc = term**2 - 4._wp*c%L**2*(B%L(1)**2/rho%L)
        disc = max(disc, 0._wp)
        c_fast%L = sqrt(0.5_wp*(term + sqrt(disc)))
        B2 = B%R(1)**2 + B%R(2)**2 + B%R(3)**2
        term = c%R**2 + B2/rho%R
        disc = term**2 - 4._wp*c%R**2*(B%R(1)**2/rho%R)
        disc = max(disc, 0._wp)
        c_fast%R = sqrt(0.5_wp*(term + sqrt(disc)))

        s%L = min(vel%L(1) - c_fast%L, vel%R(1) - c_fast%R)
        s%R = max(vel%R(1) + c_fast%R, vel%L(1) + c_fast%L)

        s_M = (((s%R - vel%R(1))*rho%R*vel%R(1) - (s%L - vel%L(1))*rho%L*vel%L(1) - pTot%R + pTot%L) &
               /((s%R - vel%R(1))*rho%R - (s%L - vel%L(1))*rho%L))

        p_star%L = pTot%L + rho%L*(s%L - vel%L(1))*(s_M - vel%L(1))
        p_star%R = pTot%R + rho%R*(s%R - vel%R(1))*(s_M - vel%R(1))

        rho_star%L = rho%L*(s%L - vel%L(1))/(s%L - s_M)
        rho_star%R = rho%R*(s%R - vel%R(1))/(s%R - s_M)

        E_star%L = ((s%L - vel%L(1))*E%L - pTot%L*vel%L(1) + p_star%L*s_M)/(s%L - s_M)
        E_star%R = ((s%R - vel%R(1))*E%R - pTot%R*vel%R(1) + p_star%R*s_M)/(s%R - s_M)

        ! Build U and U_star
        U%L(1) = rho%L;           U%R(1) = rho%R
        U%L(2) = rho%L*vel%L(1);  U%R(2) = rho%R*vel%R(1)
        U%L(3) = rho%L*vel%L(2);  U%R(3) = rho%R*vel%R(2)
        U%L(4) = rho%L*vel%L(3);  U%R(4) = rho%R*vel%R(3)
        U%L(5) = B%L(2);           U%R(5) = B%R(2)
        U%L(6) = B%L(3);           U%R(6) = B%R(3)
        U%L(7) = E%L;              U%R(7) = E%R

        U_star%L(1) = rho_star%L;           U_star%R(1) = rho_star%R
        U_star%L(2) = rho_star%L*s_M;       U_star%R(2) = rho_star%R*s_M
        U_star%L(3) = rho_star%L*vel%L(2);  U_star%R(3) = rho_star%R*vel%R(2)
        U_star%L(4) = rho_star%L*vel%L(3);  U_star%R(4) = rho_star%R*vel%R(3)
        U_star%L(5) = B%L(2);               U_star%R(5) = B%R(2)
        U_star%L(6) = B%L(3);               U_star%R(6) = B%R(3)
        U_star%L(7) = E_star%L;             U_star%R(7) = E_star%R

        ! Build F
        F%L(1) = U%L(2)
        F%L(2) = U%L(2)*vel%L(1) - B%L(1)*B%L(1) + pTot%L
        F%L(3) = U%L(2)*vel%L(2) - B%L(1)*B%L(2)
        F%L(4) = U%L(2)*vel%L(3) - B%L(1)*B%L(3)
        F%L(5) = vel%L(1)*B%L(2) - vel%L(2)*B%L(1)
        F%L(6) = vel%L(1)*B%L(3) - vel%L(3)*B%L(1)
        F%L(7) = (E%L + pTot%L)*vel%L(1) - B%L(1)*(vel%L(1)*B%L(1) + vel%L(2)*B%L(2) + vel%L(3)*B%L(3))
        F%R(1) = U%R(2)
        F%R(2) = U%R(2)*vel%R(1) - B%R(1)*B%R(1) + pTot%R
        F%R(3) = U%R(2)*vel%R(2) - B%R(1)*B%R(2)
        F%R(4) = U%R(2)*vel%R(3) - B%R(1)*B%R(3)
        F%R(5) = vel%R(1)*B%R(2) - vel%R(2)*B%R(1)
        F%R(6) = vel%R(1)*B%R(3) - vel%R(3)*B%R(1)
        F%R(7) = (E%R + pTot%R)*vel%R(1) - B%R(1)*(vel%R(1)*B%R(1) + vel%R(2)*B%R(2) + vel%R(3)*B%R(3))

        do i = 1, 7
            F_star%L(i) = F%L(i) + s%L*(U_star%L(i) - U%L(i))
            F_star%R(i) = F%R(i) + s%R*(U_star%R(i) - U%R(i))
        end do

        s_starL = s_M - abs(B%L(1))/sqrt(rho_star%L)
        s_starR = s_M + abs(B%L(1))/sqrt(rho_star%R)

        sqrt_rho_star%L = sqrt(rho_star%L)
        sqrt_rho_star%R = sqrt(rho_star%R)
        denom_ds = sqrt_rho_star%L + sqrt_rho_star%R
        sign_Bx = sign(1._wp, B%L(1))

        v_star%L = vel%L(2); v_star%R = vel%R(2)
        w_star%L = vel%L(3); w_star%R = vel%R(3)

        v_double = (sqrt_rho_star%L*v_star%L + sqrt_rho_star%R*v_star%R + (B%R(2) - B%L(2))*sign_Bx)/denom_ds
        w_double = (sqrt_rho_star%L*w_star%L + sqrt_rho_star%R*w_star%R + (B%R(3) - B%L(3))*sign_Bx)/denom_ds
        By_double = (sqrt_rho_star%L*B%R(2) + sqrt_rho_star%R*B%L(2) &
                     + sqrt_rho_star%L*sqrt_rho_star%R*(v_star%R - v_star%L)*sign_Bx)/denom_ds
        Bz_double = (sqrt_rho_star%L*B%R(3) + sqrt_rho_star%R*B%L(3) &
                     + sqrt_rho_star%L*sqrt_rho_star%R*(w_star%R - w_star%L)*sign_Bx)/denom_ds
        E_double_lr%L = E_star%L - sqrt_rho_star%L*((v_star%L*B%L(2) + w_star%L*B%L(3)) &
                                                     - (v_double*By_double + w_double*Bz_double))*sign_Bx
        E_double_lr%R = E_star%R + sqrt_rho_star%R*((v_star%R*B%R(2) + w_star%R*B%R(3)) &
                                                     - (v_double*By_double + w_double*Bz_double))*sign_Bx
        E_double = 0.5_wp*(E_double_lr%L + E_double_lr%R)

        U_double%L(1) = rho_star%L;          U_double%R(1) = rho_star%R
        U_double%L(2) = rho_star%L*s_M;      U_double%R(2) = rho_star%R*s_M
        U_double%L(3) = rho_star%L*v_double; U_double%R(3) = rho_star%R*v_double
        U_double%L(4) = rho_star%L*w_double; U_double%R(4) = rho_star%R*w_double
        U_double%L(5) = By_double;            U_double%R(5) = By_double
        U_double%L(6) = Bz_double;            U_double%R(6) = Bz_double
        U_double%L(7) = E_double;             U_double%R(7) = E_double

        if (0._wp <= s%L) then
            do i = 1, 7; F_hlld(i) = F%L(i); end do
        else if (0._wp <= s_starL) then
            do i = 1, 7; F_hlld(i) = F_star%L(i); end do
        else if (0._wp <= s_M) then
            do i = 1, 7; F_hlld(i) = F_star%L(i) + s_starL*(U_double%L(i) - U_star%L(i)); end do
        else if (0._wp <= s_starR) then
            do i = 1, 7; F_hlld(i) = F_star%R(i) + s_starR*(U_double%R(i) - U_star%R(i)); end do
        else if (0._wp <= s%R) then
            do i = 1, 7; F_hlld(i) = F_star%R(i); end do
        else
            do i = 1, 7; F_hlld(i) = F%R(i); end do
        end if

        do i = 1, 7
            F_gpu(i) = F_hlld(i)
        end do
    !$omp end target

    ! Reference: run same computation on CPU
    call compute_hlld_ref(F_ref)

    write(*,'(a,7f12.6)') "GPU: ", F_gpu
    write(*,'(a,7f12.6)') "CPU: ", F_ref
    if (all(abs(F_gpu - F_ref) < 1.e-8_wp*max(maxval(abs(F_ref)), 1._wp))) then
        write(*,'(a)') "PASS test11: full HLLD flux matches reference"
    else
        write(*,'(a)') "FAIL test11: HLLD flux mismatch"
        do i = 1, 7
            if (abs(F_gpu(i) - F_ref(i)) > 1.e-8_wp*max(abs(F_ref(i)), 1._wp)) &
                write(*,'(a,i2,a,f20.15,a,f20.15)') "  F(", i, ") gpu=", F_gpu(i), " ref=", F_ref(i)
        end do
    end if

contains

    subroutine compute_hlld_ref(Fout)
        real(wp), intent(out) :: Fout(7)
        type(vec3) :: vel, B
        type(arr7) :: U, U_star, U_double, F, F_star
        type(slr)  :: rho, pres, E, gamma, pres_mag, pTot, s, rho_star, E_star, sqrt_rho_star
        type(slr)  :: v_star, w_star, E_double_lr, p_star, c, c_fast
        real(wp)   :: F_hlld(7), s_M, s_starL, s_starR, denom_ds, sign_Bx
        real(wp)   :: v_double, w_double, By_double, Bz_double, E_double
        real(wp)   :: B2, term, disc, vel_rms_L, vel_rms_R
        integer    :: i

        rho%L = rhoL;   rho%R = rhoR
        pres%L = presL; pres%R = presR
        vel%L(1) = vxL; vel%L(2) = vyL; vel%L(3) = vzL
        vel%R(1) = vxR; vel%R(2) = vyR; vel%R(3) = vzR
        B%L(1) = BxL;   B%L(2) = ByL;   B%L(3) = BzL
        B%R(1) = BxR;   B%R(2) = ByR;   B%R(3) = BzR
        gamma%L = gam;  gamma%R = gam

        vel_rms_L = vxL**2 + vyL**2 + vzL**2
        vel_rms_R = vxR**2 + vyR**2 + vzR**2
        pres_mag%L = 0.5_wp*(BxL**2 + ByL**2 + BzL**2)
        pres_mag%R = 0.5_wp*(BxR**2 + ByR**2 + BzR**2)
        E%L = gam*presL/(gam - 1._wp) + 0.5_wp*rhoL*vel_rms_L + pres_mag%L
        E%R = gam*presR/(gam - 1._wp) + 0.5_wp*rhoR*vel_rms_R + pres_mag%R
        pTot%L = presL + pres_mag%L
        pTot%R = presR + pres_mag%R
        c%L = sqrt(gam*presL/rhoL); c%R = sqrt(gam*presR/rhoR)
        B2 = BxL**2 + ByL**2 + BzL**2; term = c%L**2 + B2/rhoL
        disc = max(term**2 - 4._wp*c%L**2*(BxL**2/rhoL), 0._wp); c_fast%L = sqrt(0.5_wp*(term + sqrt(disc)))
        B2 = BxR**2 + ByR**2 + BzR**2; term = c%R**2 + B2/rhoR
        disc = max(term**2 - 4._wp*c%R**2*(BxR**2/rhoR), 0._wp); c_fast%R = sqrt(0.5_wp*(term + sqrt(disc)))
        s%L = min(vxL - c_fast%L, vxR - c_fast%R)
        s%R = max(vxR + c_fast%R, vxL + c_fast%L)
        s_M = (((s%R - vxR)*rhoR*vxR - (s%L - vxL)*rhoL*vxL - pTot%R + pTot%L) &
               /((s%R - vxR)*rhoR - (s%L - vxL)*rhoL))
        p_star%L = pTot%L + rhoL*(s%L - vxL)*(s_M - vxL)
        rho_star%L = rhoL*(s%L - vxL)/(s%L - s_M)
        rho_star%R = rhoR*(s%R - vxR)/(s%R - s_M)
        E_star%L = ((s%L - vxL)*E%L - pTot%L*vxL + p_star%L*s_M)/(s%L - s_M)
        p_star%R = pTot%R + rhoR*(s%R - vxR)*(s_M - vxR)
        E_star%R = ((s%R - vxR)*E%R - pTot%R*vxR + p_star%R*s_M)/(s%R - s_M)

        U%L = [rhoL, rhoL*vxL, rhoL*vyL, rhoL*vzL, ByL, BzL, E%L]
        U%R = [rhoR, rhoR*vxR, rhoR*vyR, rhoR*vzR, ByR, BzR, E%R]
        U_star%L = [rho_star%L, rho_star%L*s_M, rho_star%L*vyL, rho_star%L*vzL, ByL, BzL, E_star%L]
        U_star%R = [rho_star%R, rho_star%R*s_M, rho_star%R*vyR, rho_star%R*vzR, ByR, BzR, E_star%R]
        F%L = [U%L(2), U%L(2)*vxL - BxL**2 + pTot%L, U%L(2)*vyL - BxL*ByL, U%L(2)*vzL - BxL*BzL, &
               vxL*ByL - vyL*BxL, vxL*BzL - vzL*BxL, (E%L+pTot%L)*vxL - BxL*(vxL*BxL+vyL*ByL+vzL*BzL)]
        F%R = [U%R(2), U%R(2)*vxR - BxR**2 + pTot%R, U%R(2)*vyR - BxR*ByR, U%R(2)*vzR - BxR*BzR, &
               vxR*ByR - vyR*BxR, vxR*BzR - vzR*BxR, (E%R+pTot%R)*vxR - BxR*(vxR*BxR+vyR*ByR+vzR*BzR)]
        F_star%L = F%L + s%L*(U_star%L - U%L)
        F_star%R = F%R + s%R*(U_star%R - U%R)
        s_starL = s_M - abs(BxL)/sqrt(rho_star%L)
        s_starR = s_M + abs(BxL)/sqrt(rho_star%R)
        sqrt_rho_star%L = sqrt(rho_star%L); sqrt_rho_star%R = sqrt(rho_star%R)
        denom_ds = sqrt_rho_star%L + sqrt_rho_star%R; sign_Bx = sign(1._wp, BxL)
        v_star%L = vyL; v_star%R = vyR; w_star%L = vzL; w_star%R = vzR
        v_double = (sqrt_rho_star%L*vyL + sqrt_rho_star%R*vyR + (ByR - ByL)*sign_Bx)/denom_ds
        w_double = (sqrt_rho_star%L*vzL + sqrt_rho_star%R*vzR + (BzR - BzL)*sign_Bx)/denom_ds
        By_double = (sqrt_rho_star%L*ByR + sqrt_rho_star%R*ByL + sqrt_rho_star%L*sqrt_rho_star%R*(vyR - vyL)*sign_Bx)/denom_ds
        Bz_double = (sqrt_rho_star%L*BzR + sqrt_rho_star%R*BzL + sqrt_rho_star%L*sqrt_rho_star%R*(vzR - vzL)*sign_Bx)/denom_ds
        E_double_lr%L = E_star%L - sqrt_rho_star%L*((vyL*ByL + vzL*BzL) - (v_double*By_double + w_double*Bz_double))*sign_Bx
        E_double_lr%R = E_star%R + sqrt_rho_star%R*((vyR*ByR + vzR*BzR) - (v_double*By_double + w_double*Bz_double))*sign_Bx
        E_double = 0.5_wp*(E_double_lr%L + E_double_lr%R)
        U_double%L = [rho_star%L, rho_star%L*s_M, rho_star%L*v_double, rho_star%L*w_double, By_double, Bz_double, E_double]
        U_double%R = [rho_star%R, rho_star%R*s_M, rho_star%R*v_double, rho_star%R*w_double, By_double, Bz_double, E_double]

        if (0._wp <= s%L) then; Fout = F%L
        else if (0._wp <= s_starL) then; Fout = F_star%L
        else if (0._wp <= s_M) then; Fout = F_star%L + s_starL*(U_double%L - U_star%L)
        else if (0._wp <= s_starR) then; Fout = F_star%R + s_starR*(U_double%R - U_star%R)
        else if (0._wp <= s%R) then; Fout = F_star%R
        else; Fout = F%R
        end if
    end subroutine

end program test11_hlld_full
