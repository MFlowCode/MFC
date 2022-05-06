!>
!! @file m_qbmm.f90
!! @brief Contains module m_qbmm
!! @author S. Bryngelson
!! @version 1.0
!! @date MAY 28, 2020

!> @brief This module is used to compute moment inversion via qbmm
module m_qbmm

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_qbmm_module, s_mom_inv, s_coeff

    real(kind(0d0)), parameter :: verysmall = 1.d-12
    real(kind(0d0)), allocatable, dimension(:, :, :, :, :)  :: momrhs
    integer         :: nterms !< Number of rhs terms in each moment transport equations
    type(bounds_info) :: is1, is2, is3


    real(kind(0d0)) :: momxb, momxe
    real(kind(0d0)) :: contxb, contxe
    real(kind(0d0)) :: bubxb, bubxe
    real(kind(0d0)) :: advxb, advxe
    real(kind(0d0)),allocatable, dimension(:) :: gammas, pi_infs, bubrs
    real(kind(0d0)), allocatable, dimension(:, :) :: bubmoms

!$acc declare create(verysmall, momrhs, nterms, is1, is2, is3)
!$acc declare create(momxb, momxe, bubxb, bubxe, contxb, contxe, advxb, advxe, gammas, pi_infs, bubrs, bubmoms)


contains

    subroutine s_initialize_qbmm_module()

        integer :: i1, i2, q, i, j

        if (bubble_model == 2) then
            ! Keller-Miksis without viscosity/surface tension
            nterms = 12
        else if (bubble_model == 3) then
            ! Rayleigh-Plesset with viscosity/surface tension
            nterms = 6
        end if

        !$acc update device(nterms)
        
        allocate (momrhs(3, 0:2, 0:2, nterms, nb))
        momrhs = 0d0

        ! Assigns the required RHS moments for moment transport equations
        ! The rhs%(:,3) is only to be used for R0 quadrature, not for computing X/Y indices
        do q = 1, nb
            do i1 = 0, 2; do i2 = 0, 2
                    if ((i1 + i2) <= 2) then
                        if (bubble_model == 3) then
                            momrhs(1, i1, i2,  1, q) = -1.d0 + i1
                            momrhs(2, i1, i2,  1, q) = -1.d0 + i2
                            momrhs(3, i1, i2,  1, q) = 0d0

                            momrhs(1, i1, i2,  2, q) = -1.d0 + i1
                            momrhs(2, i1, i2,  2, q) = 1.d0 + i2
                            momrhs(3, i1, i2,  2, q) = 0d0

                            momrhs(1, i1, i2,  3, q) = -1.d0 + i1 - 3.d0*gam
                            momrhs(2, i1, i2,  3, q) = -1.d0 + i2
                            momrhs(3, i1, i2,  3, q) = 3.d0*gam

                            momrhs(1, i1, i2,  4, q) = -1.d0 + i1
                            momrhs(2, i1, i2,  4, q) = 1.d0 + i2
                            momrhs(3, i1, i2,  4, q) = 0d0

                            if (Re_inv .ne. dflt_real) then
                                ! add viscosity
                                momrhs(1, i1, i2,  5, q) = -2.d0 + i1
                                momrhs(2, i1, i2,  5, q) = i2
                                momrhs(3, i1, i2,  5, q) = 0d0
                            end if

                            if (Web .ne. dflt_real) then
                                ! add surface tension
                                momrhs(1, i1, i2,  6, q) = -2.d0 + i1
                                momrhs(2, i1, i2,  6, q) = -1.d0 + i2
                                momrhs(3, i1, i2,  6, q) = 0d0
                            end if
                        else if (bubble_model == 2) then
                            ! KM with approximation of 1/(1-V/C) = 1+V/C
                            momrhs(1, i1, i2,  1, q) = -1d0 + i1
                            momrhs(2, i1, i2,  1, q) = 1d0 + i2
                            momrhs(3, i1, i2,  1, q) = 0d0

                            momrhs(1, i1, i2,  2, q) = -1d0 + i1
                            momrhs(2, i1, i2,  2, q) = 2d0 + i2
                            momrhs(3, i1, i2,  2, q) = 0d0

                            momrhs(1, i1, i2,  3, q) = -1d0 + i1
                            momrhs(2, i1, i2,  3, q) = 3d0 + i2
                            momrhs(3, i1, i2,  3, q) = 0d0

                            momrhs(1, i1, i2,  4, q) = -1d0 + i1
                            momrhs(2, i1, i2,  4, q) = -1d0 + i2
                            momrhs(3, i1, i2,  4, q) = 0d0

                            momrhs(1, i1, i2,  5, q) = -1d0 + i1
                            momrhs(2, i1, i2,  5, q) = i2
                            momrhs(3, i1, i2,  5, q) = 0d0

                            momrhs(1, i1, i2,  6, q) = -1d0 + i1
                            momrhs(2, i1, i2,  6, q) = 1d0 + i2
                            momrhs(3, i1, i2,  6, q) = 0d0

                            momrhs(1, i1, i2,  7, q) = -1d0 + i1 - 3d0*gam
                            momrhs(2, i1, i2,  7, q) = -1d0 + i2
                            momrhs(3, i1, i2,  7, q) = 3d0*gam

                            momrhs(1, i1, i2,  8, q) = -1d0 + i1 - 3d0*gam
                            momrhs(2, i1, i2,  8, q) = i2
                            momrhs(3, i1, i2,  8, q) = 3d0*gam

                            momrhs(1, i1, i2,  9, q) = -1d0 + i1 - 3d0*gam
                            momrhs(2, i1, i2,  9, q) = 1d0 + i2
                            momrhs(3, i1, i2,  9, q) = 3d0*gam

                            momrhs(1, i1, i2,  10, q) = -1d0 + i1 - 3d0*gam
                            momrhs(2, i1, i2,  10, q) = i2
                            momrhs(3, i1, i2,  10, q) = 3d0*gam

                            momrhs(1, i1, i2,  11, q) = -1d0 + i1 - 3d0*gam
                            momrhs(2, i1, i2,  11, q) = 1d0 + i2
                            momrhs(3, i1, i2,  11, q) = 3d0*gam

                            momrhs(1, i1, i2,  12, q) = -1d0 + i1
                            momrhs(2, i1, i2,  12, q) = 1d0 + i2
                            momrhs(3, i1, i2,  12, q) = 0d0
                        end if
                    end if
                end do; end do
        end do

        !$acc update device(momrhs)

        momxb = mom_idx%beg; momxe = mom_idx%end
        bubxb = bub_idx%beg; bubxe = bub_idx%end
        advxb = adv_idx%beg; advxe = adv_idx%end
        contxb = cont_idx%beg; contxe = cont_idx%end
!$acc update device(momxb, momxe, bubxb, bubxe, advxb, advxe, contxb, contxe)

        allocate(gammas(1:num_fluids))
        allocate(pi_infs(1:num_fluids))

        allocate(bubrs(1:nb))

        allocate(bubmoms(1:nb, 1:nmom))

        do i = 1, num_fluids
            gammas(i) = fluid_pp(i)%gamma
            pi_infs(i) = fluid_pp(i)%pi_inf
        end do
!$acc update device(gammas, pi_infs)

        do i = 1, nb
            bubrs(i) = bub_idx%rs(i)
        end do
!$acc update device(bubrs)

        do j = 1, nmom
            do i = 1, nb
                bubmoms(i,j) = bub_idx%moms(i,j)
            end do
        end do
!$acc update device(bubmoms)


    end subroutine s_initialize_qbmm_module

    subroutine s_coeff(pres, rho, c, coeffs)
!$acc routine seq
        real(kind(0.d0)), intent(IN) :: pres, rho, c
        real(kind(0.d0)), dimension(:, 0:, 0:), intent(OUT) :: coeffs
        integer :: i1, i2

        coeffs = 0d0

        do i1 = 0, 2; do i2 = 0, 2
                if ((i1 + i2) <= 2) then
                    if (bubble_model == 3) then
                        ! RPE
                        coeffs(1, i1, i2) = -1d0*i2*pres/rho
                        coeffs(2, i1, i2) = -3d0*i2/2d0
                        coeffs(3, i1, i2) = i2/rho
                        coeffs(4, i1, i2) = i1
                        if (Re_inv /= dflt_real) coeffs(5, i1, i2) = -4d0*i2*Re_inv/rho
                        if (Web /= dflt_real) coeffs(6, i1, i2) = -2d0*i2/Web/rho
                    else if (bubble_model == 2) then
                        ! KM with approximation of 1/(1-V/C) = 1+V/C
                        coeffs(1, i1, i2) = -3d0*i2/2d0
                        coeffs(2, i1, i2) = -i2/c
                        coeffs(3, i1, i2) = i2/(2d0*c*c)
                        coeffs(4, i1, i2) = -i2*pres/rho
                        coeffs(5, i1, i2) = -2d0*i2*pres/(c*rho)
                        coeffs(6, i1, i2) = -i2*pres/(c*c*rho)
                        coeffs(7, i1, i2) = i2/rho
                        coeffs(8, i1, i2) = 2d0*i2/(c*rho)
                        coeffs(9, i1, i2) = i2/(c*c*rho)
                        coeffs(10, i1, i2) = -3d0*i2*gam/(c*rho)
                        coeffs(11, i1, i2) = -3d0*i2*gam/(c*c*rho)
                        coeffs(12, i1, i2) = i1
                    end if
                end if
            end do; end do

    end subroutine s_coeff

    subroutine s_mom_inv(q_prim_vf, momsp, moms3d, ix, iy, iz)

        type(scalar_field), dimension(:), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(:), intent(INOUT) :: momsp
        type(scalar_field), dimension(0:, 0:, :), intent(INOUT) :: moms3d
        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)), dimension(nmom) :: moms
        real(kind(0d0)), dimension(nb) :: Rvec
        real(kind(0d0)), dimension(nnode, nb) :: wght, abscX, abscY
        real(kind(0d0)), dimension(nterms, 0:2, 0:2) :: mom3d_terms, coeff
        real(kind(0d0)) :: pres, rho, nbub, c, alf
        real(kind(0d0)) :: n_tait, B_tait
        real(kind(0d0)), dimension(3) :: momrhs_vec

        integer :: j, k, l, q, r, s !< Loop variables
        integer :: id1, id2, id3
        integer :: i1, i2

        is1 = ix; is2 = iy; is3 = iz

        !$acc update device(is1, is2, is3)


!$acc parallel loop collapse(3) gang vector default(present) private(moms, Rvec, wght, abscX, abscY, mom3d_terms, coeff, momrhs_vec)
        do id3 = is3%beg, is3%end 
            do id2 = is2%beg, is2%end 
                do id1 = is1%beg, is1%end

                alf = q_prim_vf(alf_idx)%sf(id1, id2, id3)
                pres = q_prim_vf(E_idx)%sf(id1, id2, id3)
                rho = q_prim_vf(contxb)%sf(id1, id2, id3)
                if (bubble_model == 2) then
                    n_tait = gammas(1)
                    n_tait = 1.d0/n_tait + 1.d0 !make this the usual little 'gamma'
                    B_tait = pi_infs(1)
                    c = n_tait*(pres + B_tait)/(rho*(1.d0 - alf))
                    if (c > 0.d0) then
                        c = DSQRT(c)
                    else
                        c = sgm_eps
                    end if
                end if

                call s_coeff(pres, rho, c, coeff)



                ! SHB: Manually adjusted pressure here for no-coupling case
                ! pres = 1d0/0.3d0

                if (alf > small_alf) then

                    !$acc loop seq
                    do q = 1, nb
                        Rvec(q) = q_prim_vf(bubrs(q))%sf(id1, id2, id3)
                    end do

                    call s_comp_n_from_prim(alf, Rvec, nbub)

                    !$acc loop seq
                    do q = 1, nb
                        !$acc loop seq
                        do r = 1, nmom
                            moms(r) = q_prim_vf(bubmoms(q, r))%sf(id1, id2, id3)
                        end do

                        ! IF(id1==0) THEN
                        !     PRINT*, 'pres: ', pres
                        !     PRINT*, 'nb : ', nbub
                        !     PRINT*, 'alf: ', alf
                        !     DO s = 1,nmom
                        !         PRINT*, 'mom: ', moms(s)
                        !     END DO
                        ! END IF

                        call s_chyqmom(moms, wght(:, q), abscX(:, q), abscY(:, q))



                        !$acc loop seq
                        do j = 1, nterms
                            !$acc loop seq
                            do i2 = 0, 2 
                                !$acc loop seq
                                do i1 = 0, 2
                                    if ((i1 + i2) <= 2) then

                                        mom3d_terms(j, i1, i2) = coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2,  j, q)) &
                                                                 *f_quad2D(abscX(:, q), abscY(:, q), wght(:, q), momrhs(:, i1, i2, j, q))
                                    end if
                                end do 
                            end do
                        end do

                        !$acc loop seq
                        do i1 = 0, 2
                            !$acc loop seq
                             do i2 = 0, 2
                                if ((i1 + i2) <= 2) then
                                    moms3d(i1, i2, q)%sf(id1, id2, id3) = nbub*sum(mom3d_terms(:, i1, i2))
                                    ! IF (moms3d(i1,i2,q)%sf(id1,id2,id3) .NE. moms3d(i1,i2,q)%sf(id1,id2,id3)) THEN
                                    !     PRINT*, 'nan in mom3d', i1,i2,id1
                                    !     PRINT*, 'nbu: ', nbub
                                    !     PRINT*, 'alf: ', alf
                                    !     PRINT*, 'moms: ', moms(:)
                                    !     CALL s_mpi_abort()
                                    ! END IF
                                end if
                            end do
                        end do
                    end do

                    momsp(1)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3d0, 0d0, 0d0)
                    momsp(2)%sf(id1, id2, id3) = 4.d0*pi*nbub*f_quad(abscX, abscY, wght, 2d0, 1d0, 0d0)
                    momsp(3)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3d0, 2d0, 0d0)
                    if (abs(gam - 1.d0) <= 1.d-4) then
                        ! Gam \approx 1, don't risk imaginary quadrature
                        momsp(4)%sf(id1, id2, id3) = 1.d0
                    else
                        momsp(4)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3d0*(1d0 - gam), 0d0, 3d0*gam)
                    end if

                    !!$acc loop seq
                    !do i1 = 1, 4
                       ! if (momsp(i1)%sf(id1, id2, id3) /= momsp(i1)%sf(id1, id2, id3)) then
                       !     print *, 'NaN in sp moment', i1, 'location', id1, id2, id3
                       !     print *, 'Rs', Rvec(:)
                       !     print *, 'alpha', alf
                       !     print *, 'nbub', nbub
                       !     print *, 'abscX', abscX(:, :)
                       !     print *, 'abscY', abscY(:, :)
                       !     print *, 'wght', wght(:, :)
                        !    call s_mpi_abort()
                        !end if
                    !end do
                else
                    !$acc loop seq
                    do q = 1, nb
                        !$acc loop seq
                        do i1 = 0, 2
                            !$acc loop seq 
                            do i2 = 0, 2
                                moms3d(i1, i2, q)%sf(id1, id2, id3) = 0d0
                            end do 
                        end do
                    end do

                    momsp(1)%sf(id1, id2, id3) = 0d0
                    momsp(2)%sf(id1, id2, id3) = 0d0
                    momsp(3)%sf(id1, id2, id3) = 0d0
                    momsp(4)%sf(id1, id2, id3) = 0d0

                end if

            end do
         end do
      end do

    end subroutine s_mom_inv

    subroutine s_chyqmom(momin, wght, abscX, abscY)
!$acc routine seq
        real(kind(0d0)), dimension(:), intent(INOUT) :: wght, abscX, abscY
        real(kind(0d0)), dimension(:), intent(IN) :: momin

        real(kind(0d0)), dimension(0:nmom, 0:nmom) :: moms
        real(kind(0d0)), dimension(3) :: M1, M3
        real(kind(0d0)), dimension(2) :: myrho, myrho3, up, up3, Vf
        real(kind(0d0)) :: bu, bv, d20, d11, d02, c20, c11, c02
        real(kind(0d0)) :: mu2avg, mu2, vp21, vp22, rho21, rho22

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
        c20 = d20 - bu**2d0; 
        c11 = d11 - bu*bv; 
        c02 = d02 - bv**2d0; 
        M1 = (/1d0, 0d0, c20/)
        call s_hyqmom(myrho, up, M1)
        Vf = c11*up/c20

        mu2avg = c02 - sum(myrho(:)*(Vf(:)**2d0))
        mu2avg = maxval((/mu2avg, 0d0/))
        mu2 = mu2avg
        M3 = (/1d0, 0d0, mu2/)
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
        real(kind(0d0)), dimension(2), intent(INOUT) :: frho, fup
        real(kind(0d0)), dimension(3), intent(IN) :: fmom
        real(kind(0d0)) :: bu, d2, c2

        bu = fmom(2)/fmom(1)
        d2 = fmom(3)/fmom(1)
        c2 = d2 - bu**2d0
        frho(1) = fmom(1)/2d0; 
        frho(2) = fmom(1)/2d0; 
        c2 = maxval((/c2, verysmall/))
        fup(1) = bu - DSQRT(c2)
        fup(2) = bu + DSQRT(c2)

    end subroutine s_hyqmom

    function f_quad(abscX, abscY, wght, q, r, s)
    !$acc routine seq
        real(kind(0.d0)), dimension(:, :), intent(IN) :: abscX, abscY, wght
        real(kind(0.d0)), intent(IN) :: q, r, s
        real(kind(0.d0)) :: f_quad_RV, f_quad
        integer :: i

        f_quad = 0d0
        do i = 1, nb
            f_quad_RV = sum(wght(:, i)*(abscX(:, i)**q)*(abscY(:, i)**r))
            f_quad = f_quad + weight(i)*(R0(i)**s)*f_quad_RV
        end do

    end function f_quad

    function f_quad2D(abscX, abscY, wght, pow)
    !$acc routine seq
        real(kind(0.d0)), dimension(:), intent(IN) :: abscX, abscY, wght
        real(kind(0.d0)), dimension(:), intent(IN) :: pow
        real(kind(0.d0)) :: f_quad2D

        f_quad2D = sum(wght(:)*(abscX(:)**pow(1))*(abscY(:)**pow(2)))
    end function f_quad2D

end module m_qbmm
