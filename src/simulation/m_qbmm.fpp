!>
!! @file m_qbmm.f90
!! @brief Contains module m_qbmm

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief This module is used to compute moment inversion via qbmm
module m_qbmm

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_qbmm_module, s_mom_inv, s_coeff

    real(wp), allocatable, dimension(:, :, :, :, :) :: momrhs

    #:if MFC_CASE_OPTIMIZATION
        integer, parameter :: nterms = ${nterms}$
    #:else
        integer :: nterms
    #:endif

    type(int_bounds_info) :: is1, is2, is3

    integer, allocatable, dimension(:)    :: bubrs
    integer, allocatable, dimension(:, :) :: bubmoms

!$acc declare create(momrhs, nterms, is1, is2, is3)
!$acc declare create(bubrs, bubmoms)

contains

    subroutine s_initialize_qbmm_module()

        integer :: i1, i2, q, i, j

        #:if not MFC_CASE_OPTIMIZATION

            if (bubble_model == 2) then
                ! Keller-Miksis without viscosity/surface tension
                nterms = 12
            else if (bubble_model == 3) then
                ! Rayleigh-Plesset with viscosity/surface tension
                nterms = 6
            end if

            !$acc update device(nterms)

        #:endif

        @:ALLOCATE(momrhs(3, 0:2, 0:2, nterms, nb))
        momrhs = 0._wp

        ! Assigns the required RHS moments for moment transport equations
        ! The rhs%(:,3) is only to be used for R0 quadrature, not for computing X/Y indices
        do q = 1, nb
            do i1 = 0, 2; do i2 = 0, 2
                    if ((i1 + i2) <= 2) then
                        if (bubble_model == 3) then
                            momrhs(1, i1, i2, 1, q) = -1._wp + i1
                            momrhs(2, i1, i2, 1, q) = -1._wp + i2
                            momrhs(3, i1, i2, 1, q) = 0._wp

                            momrhs(1, i1, i2, 2, q) = -1._wp + i1
                            momrhs(2, i1, i2, 2, q) = 1._wp + i2
                            momrhs(3, i1, i2, 2, q) = 0._wp

                            momrhs(1, i1, i2, 3, q) = -1._wp + i1 - 3._wp*gam
                            momrhs(2, i1, i2, 3, q) = -1._wp + i2
                            momrhs(3, i1, i2, 3, q) = 3._wp*gam

                            momrhs(1, i1, i2, 4, q) = -1._wp + i1
                            momrhs(2, i1, i2, 4, q) = 1._wp + i2
                            momrhs(3, i1, i2, 4, q) = 0._wp

                            if (Re_inv /= dflt_real) then
                                ! add viscosity
                                momrhs(1, i1, i2, 5, q) = -2._wp + i1
                                momrhs(2, i1, i2, 5, q) = i2
                                momrhs(3, i1, i2, 5, q) = 0._wp
                            end if

                            if (Web /= dflt_real) then
                                ! add surface tension
                                momrhs(1, i1, i2, 6, q) = -2._wp + i1
                                momrhs(2, i1, i2, 6, q) = -1._wp + i2
                                momrhs(3, i1, i2, 6, q) = 0._wp
                            end if
                        else if (bubble_model == 2) then
                            ! KM with approximation of 1/(1-V/C) = 1+V/C
                            momrhs(1, i1, i2, 1, q) = -1._wp + i1
                            momrhs(2, i1, i2, 1, q) = 1._wp + i2
                            momrhs(3, i1, i2, 1, q) = 0._wp

                            momrhs(1, i1, i2, 2, q) = -1._wp + i1
                            momrhs(2, i1, i2, 2, q) = 2._wp + i2
                            momrhs(3, i1, i2, 2, q) = 0._wp

                            momrhs(1, i1, i2, 3, q) = -1._wp + i1
                            momrhs(2, i1, i2, 3, q) = 3._wp + i2
                            momrhs(3, i1, i2, 3, q) = 0._wp

                            momrhs(1, i1, i2, 4, q) = -1._wp + i1
                            momrhs(2, i1, i2, 4, q) = -1._wp + i2
                            momrhs(3, i1, i2, 4, q) = 0._wp

                            momrhs(1, i1, i2, 5, q) = -1._wp + i1
                            momrhs(2, i1, i2, 5, q) = i2
                            momrhs(3, i1, i2, 5, q) = 0._wp

                            momrhs(1, i1, i2, 6, q) = -1._wp + i1
                            momrhs(2, i1, i2, 6, q) = 1._wp + i2
                            momrhs(3, i1, i2, 6, q) = 0._wp

                            momrhs(1, i1, i2, 7, q) = -1._wp + i1 - 3._wp*gam
                            momrhs(2, i1, i2, 7, q) = -1._wp + i2
                            momrhs(3, i1, i2, 7, q) = 3._wp*gam

                            momrhs(1, i1, i2, 8, q) = -1._wp + i1 - 3._wp*gam
                            momrhs(2, i1, i2, 8, q) = i2
                            momrhs(3, i1, i2, 8, q) = 3._wp*gam

                            momrhs(1, i1, i2, 9, q) = -1._wp + i1 - 3._wp*gam
                            momrhs(2, i1, i2, 9, q) = 1._wp + i2
                            momrhs(3, i1, i2, 9, q) = 3._wp*gam

                            momrhs(1, i1, i2, 10, q) = -1._wp + i1 - 3._wp*gam
                            momrhs(2, i1, i2, 10, q) = i2
                            momrhs(3, i1, i2, 10, q) = 3._wp*gam

                            momrhs(1, i1, i2, 11, q) = -1._wp + i1 - 3._wp*gam
                            momrhs(2, i1, i2, 11, q) = 1._wp + i2
                            momrhs(3, i1, i2, 11, q) = 3._wp*gam

                            momrhs(1, i1, i2, 12, q) = -1._wp + i1
                            momrhs(2, i1, i2, 12, q) = 1._wp + i2
                            momrhs(3, i1, i2, 12, q) = 0._wp
                        end if
                    end if
                end do; end do
        end do

        !$acc update device(momrhs)

        @:ALLOCATE(bubrs(1:nb))
        @:ALLOCATE(bubmoms(1:nb, 1:nmom))

        do i = 1, nb
            bubrs(i) = bub_idx%rs(i)
        end do
!$acc update device(bubrs)

        do j = 1, nmom
            do i = 1, nb
                bubmoms(i, j) = bub_idx%moms(i, j)
            end do
        end do
!$acc update device(bubmoms)

    end subroutine s_initialize_qbmm_module

    subroutine s_coeff(pres, rho, c, coeffs)
!$acc routine seq
        real(wp), intent(IN) :: pres, rho, c
        real(wp), dimension(nterms, 0:2, 0:2), intent(OUT) :: coeffs
        integer :: i1, i2

        coeffs = 0._wp

        do i2 = 0, 2; do i1 = 0, 2
                if ((i1 + i2) <= 2) then
                    if (bubble_model == 3) then
                        ! RPE
                        coeffs(1, i1, i2) = -1._wp*i2*pres/rho
                        coeffs(2, i1, i2) = -3._wp*i2/2._wp
                        coeffs(3, i1, i2) = i2/rho
                        coeffs(4, i1, i2) = i1
                        if (Re_inv /= dflt_real) coeffs(5, i1, i2) = -4._wp*i2*Re_inv/rho
                        if (Web /= dflt_real) coeffs(6, i1, i2) = -2._wp*i2/Web/rho
                    else if (bubble_model == 2) then
                        ! KM with approximation of 1/(1-V/C) = 1+V/C
                        coeffs(1, i1, i2) = -3._wp*i2/2._wp
                        coeffs(2, i1, i2) = -i2/c
                        coeffs(3, i1, i2) = i2/(2._wp*c*c)
                        coeffs(4, i1, i2) = -i2*pres/rho
                        coeffs(5, i1, i2) = -2._wp*i2*pres/(c*rho)
                        coeffs(6, i1, i2) = -i2*pres/(c*c*rho)
                        coeffs(7, i1, i2) = i2/rho
                        coeffs(8, i1, i2) = 2._wp*i2/(c*rho)
                        coeffs(9, i1, i2) = i2/(c*c*rho)
                        coeffs(10, i1, i2) = -3._wp*i2*gam/(c*rho)
                        coeffs(11, i1, i2) = -3._wp*i2*gam/(c*c*rho)
                        coeffs(12, i1, i2) = i1
                    end if
                end if
            end do; end do

    end subroutine s_coeff

    subroutine s_mom_inv(q_prim_vf, momsp, moms3d, ix, iy, iz)

        type(scalar_field), dimension(:), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(:), intent(INOUT) :: momsp
        type(scalar_field), dimension(0:, 0:, :), intent(INOUT) :: moms3d
        type(int_bounds_info), intent(IN) :: ix, iy, iz

        real(wp), dimension(nmom) :: moms
        real(wp), dimension(nb) :: Rvec
        real(wp), dimension(nnode, nb) :: wght, abscX, abscY
        real(wp), dimension(nterms, 0:2, 0:2) :: mom3d_terms, coeff
        real(wp) :: pres, rho, nbub, c, alf, R3, momsum
        real(wp) :: start, finish
        real(wp) :: n_tait, B_tait

        integer :: j, k, l, q, r, s !< Loop variables
        integer :: id1, id2, id3
        integer :: i1, i2

        is1 = ix; is2 = iy; is3 = iz

        !$acc update device(is1, is2, is3)


!$acc parallel loop collapse(3) gang vector default(present) private(moms, wght, abscX, abscY, coeff)
        do id3 = is3%beg, is3%end
            do id2 = is2%beg, is2%end
                do id1 = is1%beg, is1%end

                    alf = q_prim_vf(alf_idx)%sf(id1, id2, id3)
                    pres = q_prim_vf(E_idx)%sf(id1, id2, id3)
                    rho = q_prim_vf(contxb)%sf(id1, id2, id3)
                    if (bubble_model == 2) then
                        n_tait = gammas(1)
                        n_tait = 1._wp/n_tait + 1._wp !make this the usual little 'gamma'
                        B_tait = pi_infs(1)
                        B_tait = B_tait*(n_tait-1)/n_tait ! make this the usual pi_inf
                        c = n_tait*(pres + B_tait)/(rho*(1._wp - alf))
                        if (c > 0._wp) then
                            c = sqrt(c)
                        else
                            c = sgm_eps
                        end if
                    end if

                    call s_coeff(pres, rho, c, coeff)

                    ! SHB: Manually adjusted pressure here for no-coupling case
                    ! pres = 1._wp/0.3_wp

                    if (alf > small_alf) then

                        R3 = 0._wp

                        !$acc loop seq
                        do q = 1, nb
                            R3 = R3 + weight(q)*q_prim_vf(bubrs(q))%sf(id1, id2, id3)**3._wp
                        end do

                        nbub = (3._wp/(4._wp*pi))*alf/R3

                        !$acc loop seq
                        do q = 1, nb
                            !$acc loop seq
                            do r = 1, nmom
                                moms(r) = q_prim_vf(bubmoms(q, r))%sf(id1, id2, id3)
                            end do

                           

                            call s_chyqmom(moms, wght(:, q), abscX(:, q), abscY(:, q))


                            !$acc loop seq
                            do i2 = 0, 2
                                !$acc loop seq
                                do i1 = 0, 2
                                    if ((i1 + i2) <= 2) then
                                        momsum = 0._wp
                                        !$acc loop seq
                                        do j = 1, nterms           
                                            momsum = momsum  + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q)) &
                                                            *f_quad2D(abscX(:, q), abscY(:, q), wght(:, q), momrhs(:, i1, i2, j, q))
                                        end do
                                        moms3d(i1, i2, q)%sf(id1, id2, id3) = nbub * momsum

                                    end if
                                end do
                            end do

                            
                        end do

                        momsp(1)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3._wp, 0._wp, 0._wp)
                        momsp(2)%sf(id1, id2, id3) = 4._wp*pi*nbub*f_quad(abscX, abscY, wght, 2._wp, 1._wp, 0._wp)
                        momsp(3)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3._wp, 2._wp, 0._wp)
                        if (abs(gam - 1._wp) <= (1._wp * (10._wp ** -(4)))) then
                            ! Gam \approx 1, don't risk imaginary quadrature
                            momsp(4)%sf(id1, id2, id3) = 1._wp
                        else
                            momsp(4)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3._wp*(1._wp - gam), 0._wp, 3._wp*gam)
                        end if

                    
                    else
                        !$acc loop seq
                        do q = 1, nb
                            !$acc loop seq
                            do i1 = 0, 2
                                !$acc loop seq
                                do i2 = 0, 2
                                    moms3d(i1, i2, q)%sf(id1, id2, id3) = 0._wp
                                end do
                            end do
                        end do

                        momsp(1)%sf(id1, id2, id3) = 0._wp
                        momsp(2)%sf(id1, id2, id3) = 0._wp
                        momsp(3)%sf(id1, id2, id3) = 0._wp
                        momsp(4)%sf(id1, id2, id3) = 0._wp

                    end if

                end do
            end do
        end do


    end subroutine s_mom_inv

    subroutine s_chyqmom(momin, wght, abscX, abscY)
!$acc routine seq
        real(wp), dimension(nnode), intent(INOUT) :: wght, abscX, abscY
        real(wp), dimension(nmom), intent(IN) :: momin

        real(wp), dimension(0:2, 0:2) :: moms
        real(wp), dimension(3) :: M1, M3
        real(wp), dimension(2) :: myrho, myrho3, up, up3, Vf
        real(wp) :: bu, bv, d20, d11, d02, c20, c11, c02
        real(wp) :: mu2avg, mu2, vp21, vp22, rho21, rho22

        moms(0, 0) = momin(1)
        moms(1, 0) = momin(2)
        moms(0, 1) = momin(3)
        moms(2, 0) = momin(4)
        moms(1, 1) = momin(5)
        moms(0, 2) = momin(6)

        bu = moms(1, 0)/moms(0, 0)
        bv = moms(0, 1)/moms(0, 0)
        d20 = moms(2, 0)/moms(0, 0)
        d11 = moms(1, 1)/moms(0, 0)
        d02 = moms(0, 2)/moms(0, 0)
        c20 = d20 - bu**2._wp; 
        c11 = d11 - bu*bv; 
        c02 = d02 - bv**2._wp; 
        M1 = (/1._wp, 0._wp, c20/)
        call s_hyqmom(myrho, up, M1)
        Vf = c11*up/c20

        mu2avg = c02 - sum(myrho(:)*(Vf(:)**2._wp))
        mu2avg = maxval((/mu2avg, 0._wp/))
        mu2 = mu2avg
        M3 = (/1._wp, 0._wp, mu2/)
        call s_hyqmom(myrho3, up3, M3)

        vp21 = up3(1)
        vp22 = up3(2)
        rho21 = myrho3(1)
        rho22 = myrho3(2)

        wght(1) = myrho(1)*rho21
        wght(2) = myrho(1)*rho22
        wght(3) = myrho(2)*rho21
        wght(4) = myrho(2)*rho22
        wght = moms(0, 0)*wght

        abscX(1) = up(1)
        abscX(2) = up(1)
        abscX(3) = up(2)
        abscX(4) = up(2)
        abscX = bu + abscX

        abscY(1) = Vf(1) + vp21
        abscY(2) = Vf(1) + vp22
        abscY(3) = Vf(2) + vp21
        abscY(4) = Vf(2) + vp22
        abscY = bv + abscY

    end subroutine s_chyqmom

    subroutine s_hyqmom(frho, fup, fmom)
        !$acc routine seq
        real(wp), dimension(2), intent(INOUT) :: frho, fup
        real(wp), dimension(3), intent(IN) :: fmom
        real(wp) :: bu, d2, c2

        bu = fmom(2)/fmom(1)
        d2 = fmom(3)/fmom(1)
        c2 = d2 - bu**2._wp
        frho(1) = fmom(1)/2._wp; 
        frho(2) = fmom(1)/2._wp; 
        c2 = maxval((/c2, verysmall/))
        fup(1) = bu - sqrt(c2)
        fup(2) = bu + sqrt(c2)

    end subroutine s_hyqmom

    function f_quad(abscX, abscY, wght, q, r, s)
        !$acc routine seq
        real(wp), dimension(nnode, nb), intent(IN) :: abscX, abscY, wght
        real(wp), intent(IN) :: q, r, s
        real(wp) :: f_quad_RV, f_quad
        integer :: i

        f_quad = 0._wp
        do i = 1, nb
            f_quad_RV = sum(wght(:, i)*(abscX(:, i)**q)*(abscY(:, i)**r))
            f_quad = f_quad + weight(i)*(R0(i)**s)*f_quad_RV
        end do

    end function f_quad

    function f_quad2D(abscX, abscY, wght, pow)
        !$acc routine seq
        real(wp), dimension(nnode), intent(IN) :: abscX, abscY, wght
        real(wp), dimension(3), intent(IN) :: pow
        real(wp) :: f_quad2D

        f_quad2D = sum(wght(:)*(abscX(:)**pow(1))*(abscY(:)**pow(2)))
    end function f_quad2D

end module m_qbmm
