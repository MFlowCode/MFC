!>
!! @file
!! @brief Contains module m_thinc

#:include 'macros.fpp'

!> @brief THINC and MTHINC interface compression for volume fraction sharpening. THINC (int_comp=1): 1D directional interface
!! compression applied after MUSCL/WENO reconstruction. Uses hyperbolic tangent profile per dimension. MTHINC (int_comp=2):
!! Multi-dimensional THINC that reconstructs a tanh profile oriented along the interface normal (computed from the gradient of
!! alpha), then integrates that profile over each cell face using Gaussian quadrature. A Newton iteration enforces the conservation
!! constraint (cell-averaged alpha is preserved). Reference: B. Xie and F. Xiao, "Toward efficient and accurate interface capturing
!! on arbitrary hybrid unstructured grids: The THINC method with quadratic surface representation and Gaussian quadrature," Journal
!! of Computational Physics, vol. 349, pp. 415-440, 2017.
module m_thinc

    use m_derived_types
    use m_global_parameters
    use m_helper

#ifdef MFC_OpenACC
    use openacc
#endif

    implicit none

    private; public :: s_initialize_thinc_module, s_thinc_compression, s_compute_mthinc_normals, s_finalize_thinc_module

    !> 3-point Gauss-Legendre quadrature on [-1/2, 1/2]
    !> Node locations: +-sqrt(3/5)/2, 0
    real(wp) :: gq3_pts(3) = [-5e-1_wp*0.7745966692414834_wp, 0._wp, 5e-1_wp*0.7745966692414834_wp]
    !> Weights: 5/18, 8/18, 5/18
    real(wp) :: gq3_wts(3) = [5._wp/18._wp, 8._wp/18._wp, 5._wp/18._wp]
    $:GPU_DECLARE(copyin='[gq3_pts, gq3_wts]')

    real(wp), allocatable, dimension(:,:,:,:) :: mthinc_nhat !> Unit normal vector
    real(wp), allocatable, dimension(:,:,:)   :: mthinc_d    !> Interface position parameter
    $:GPU_DECLARE(create='[mthinc_nhat, mthinc_d]')

contains

    !> @brief Stable difference: log_cosh(a+h) - log_cosh(a-h) = 2*atanh(tanh(a)*tanh(h)). Avoids catastrophic cancellation when h
    !! is small relative to a.
    function f_log_cosh_diff(a, h) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: a, h
        real(wp)             :: res, t

        t = tanh(a)*tanh(h)
        res = 2._wp*atanh(sign(min(abs(t), 1._wp - epsilon(1._wp)), t))

    end function f_log_cosh_diff

    !> @brief Analytical 1-D integral of the THINC function
    function f_thinc_integral_1d(a, b) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: a, b
        real(wp)             :: res

        if (abs(b) < verysmall) then
            res = 5e-1_wp*(1._wp + tanh(a))
        else
            res = 5e-1_wp + f_log_cosh_diff(a, 5e-1_wp*b)/(2._wp*b)
        end if

    end function f_thinc_integral_1d

    !> @brief Volume integral of H(xi) = 0.5*(1 + tanh(beta*(n.xi + d))) over the cell [-1/2, 1/2]^ndim
    function f_mthinc_volume_integral(n1, n2, n3, d, beta, ndim) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: n1, n2, n3, d, beta
        integer, intent(in)  :: ndim
        real(wp)             :: res, a
        integer              :: q1, q2

        if (ndim == 1) then
            res = f_thinc_integral_1d(beta*d, beta*n1)
        else if (ndim == 2) then
            res = 0._wp
            do q1 = 1, 3
                a = beta*(n2*gq3_pts(q1) + d)
                res = res + gq3_wts(q1)*f_thinc_integral_1d(a, beta*n1)
            end do
        else
            res = 0._wp
            do q1 = 1, 3
                do q2 = 1, 3
                    a = beta*(n2*gq3_pts(q1) + n3*gq3_pts(q2) + d)
                    res = res + gq3_wts(q1)*gq3_wts(q2)*f_thinc_integral_1d(a, beta*n1)
                end do
            end do
        end if

    end function f_mthinc_volume_integral

    !> @brief Derivative dV/dd of the volume integral (for Newton iteration)
    function f_mthinc_volume_integral_dd(n1, n2, n3, d, beta, ndim) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: n1, n2, n3, d, beta
        integer, intent(in)  :: ndim
        real(wp)             :: res, th
        integer              :: q1, q2, q3

        res = 0._wp
        if (ndim == 1) then
            do q1 = 1, 3
                th = tanh(beta*(n1*gq3_pts(q1) + d))
                res = res + gq3_wts(q1)*(1._wp - th*th)
            end do
        else if (ndim == 2) then
            do q1 = 1, 3
                do q2 = 1, 3
                    th = tanh(beta*(n1*gq3_pts(q1) + n2*gq3_pts(q2) + d))
                    res = res + gq3_wts(q1)*gq3_wts(q2)*(1._wp - th*th)
                end do
            end do
        else
            do q1 = 1, 3
                do q2 = 1, 3
                    do q3 = 1, 3
                        th = tanh(beta*(n1*gq3_pts(q1) + n2*gq3_pts(q2) + n3*gq3_pts(q3) + d))
                        res = res + gq3_wts(q1)*gq3_wts(q2)*gq3_wts(q3)*(1._wp - th*th)
                    end do
                end do
            end do
        end if
        res = 5e-1_wp*beta*res

    end function f_mthinc_volume_integral_dd

    !> @brief Solve for the interface-position parameter d
    function f_mthinc_solve_d(n1, n2, n3, beta, alpha_cell, ndim) result(d)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: n1, n2, n3, beta, alpha_cell
        integer, intent(in)  :: ndim
        real(wp)             :: d, V, residual, dV
        integer              :: iter

        d = 0._wp
        do iter = 1, 30
            V = f_mthinc_volume_integral(n1, n2, n3, d, beta, ndim)
            residual = V - alpha_cell
            if (abs(residual) < verysmall) exit
            dV = f_mthinc_volume_integral_dd(n1, n2, n3, d, beta, ndim)
            if (abs(dV) < verysmall) exit
            d = d - residual/dV
        end do

    end function f_mthinc_solve_d

    !> @brief Face-averaged THINC function at a cell face. face_dir: 1=x, 2=y, 3=z. face_pos: -0.5 (low) or +0.5 (high). The face
    !! coordinate in face_dir is fixed; remaining directions are integrated over [-1/2, 1/2] (analytically along one, Gauss for
    !! others).
    function f_mthinc_face_average(n1, n2, n3, d, beta, face_dir, face_pos, ndim) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: n1, n2, n3, d, beta, face_pos
        integer, intent(in)  :: face_dir, ndim
        real(wp)             :: res, a, n_face, n_t0, n_t1
        integer              :: q

        if (ndim == 1) then
            res = 5e-1_wp*(1._wp + tanh(beta*(n1*face_pos + d)))
        else if (ndim == 2) then
            ! One transverse direction
            if (face_dir == 1) then
                n_face = n1; n_t0 = n2
            else
                n_face = n2; n_t0 = n1
            end if
            a = beta*(n_face*face_pos + d)
            res = f_thinc_integral_1d(a, beta*n_t0)
        else
            ! Two transverse directions
            if (face_dir == 1) then
                n_face = n1; n_t0 = n2; n_t1 = n3
            else if (face_dir == 2) then
                n_face = n2; n_t0 = n1; n_t1 = n3
            else
                n_face = n3; n_t0 = n1; n_t1 = n2
            end if
            res = 0._wp
            do q = 1, 3
                a = beta*(n_face*face_pos + n_t0*gq3_pts(q) + d)
                res = res + gq3_wts(q)*f_thinc_integral_1d(a, beta*n_t1)
            end do
        end if

    end function f_mthinc_face_average

    subroutine s_initialize_thinc_module()

        if (int_comp == 2) then
            @:ALLOCATE(mthinc_nhat(1:3, idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                       & idwbuff(3)%beg:idwbuff(3)%end))
            @:ALLOCATE(mthinc_d(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        end if

    end subroutine s_initialize_thinc_module

    !> @brief Compute the unit normal and interface-position parameter d at each interface cell
    subroutine s_compute_mthinc_normals(v_vf)

        type(scalar_field), dimension(:), intent(in) :: v_vf
        integer                                      :: j, k, l
        real(wp)                                     :: nr_x, nr_y, nr_z, nmag, nmax, ac
        type(int_bounds_info), dimension(3)          :: id_norm

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    mthinc_nhat(1, j, k, l) = 0._wp
                    mthinc_nhat(2, j, k, l) = 0._wp
                    mthinc_nhat(3, j, k, l) = 0._wp
                    mthinc_d(j, k, l) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        id_norm(1)%beg = idwbuff(1)%beg + 1; id_norm(1)%end = idwbuff(1)%end - 1
        id_norm(2)%beg = 0; id_norm(2)%end = 0
        id_norm(3)%beg = 0; id_norm(3)%end = 0
        if (n > 0) then
            id_norm(2)%beg = idwbuff(2)%beg + 1; id_norm(2)%end = idwbuff(2)%end - 1
        end if
        if (p > 0) then
            id_norm(3)%beg = idwbuff(3)%beg + 1; id_norm(3)%end = idwbuff(3)%end - 1
        end if

        ! Compute unit normal and solve for d at interior
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, nr_x, nr_y, nr_z, nmag, nmax, ac]', copyin='[id_norm]')
        do l = id_norm(3)%beg, id_norm(3)%end
            do k = id_norm(2)%beg, id_norm(2)%end
                do j = id_norm(1)%beg, id_norm(1)%end
                    ac = v_vf(eqn_idx%adv%beg)%sf(j, k, l)

                    if (ac >= ic_eps .and. ac <= 1._wp - ic_eps) then
                        nr_x = (v_vf(eqn_idx%adv%beg)%sf(j + 1, k, l) - v_vf(eqn_idx%adv%beg)%sf(j - 1, k, &
                                & l))*(x_cb(j) - x_cb(j - 1))/(x_cc(j + 1) - x_cc(j - 1))

                        nr_y = 0._wp
                        if (n > 0) then
                            nr_y = (v_vf(eqn_idx%adv%beg)%sf(j, k + 1, l) - v_vf(eqn_idx%adv%beg)%sf(j, k - 1, &
                                    & l))*(y_cb(k) - y_cb(k - 1))/(y_cc(k + 1) - y_cc(k - 1))
                        end if

                        nr_z = 0._wp
                        if (p > 0) then
                            nr_z = (v_vf(eqn_idx%adv%beg)%sf(j, k, l + 1) - v_vf(eqn_idx%adv%beg)%sf(j, k, &
                                    & l - 1))*(z_cb(l) - z_cb(l - 1))/(z_cc(l + 1) - z_cc(l - 1))
                        end if

                        nmag = sqrt(nr_x*nr_x + nr_y*nr_y + nr_z*nr_z)

                        if (nmag > verysmall) then
                            nr_x = nr_x/nmag
                            nr_y = nr_y/nmag
                            nr_z = nr_z/nmag

                            ! Snap near-grid-aligned normals to exact alignment
                            nmax = max(abs(nr_x), abs(nr_y), abs(nr_z))
                            if (abs(nr_x) < mthinc_align_tol*nmax) nr_x = 0._wp
                            if (abs(nr_y) < mthinc_align_tol*nmax) nr_y = 0._wp
                            if (abs(nr_z) < mthinc_align_tol*nmax) nr_z = 0._wp
                            nmag = sqrt(nr_x*nr_x + nr_y*nr_y + nr_z*nr_z)
                            if (nmag > verysmall) then
                                nr_x = nr_x/nmag
                                nr_y = nr_y/nmag
                                nr_z = nr_z/nmag

                                mthinc_nhat(1, j, k, l) = nr_x
                                mthinc_nhat(2, j, k, l) = nr_y
                                mthinc_nhat(3, j, k, l) = nr_z

                                mthinc_d(j, k, l) = f_mthinc_solve_d(nr_x, nr_y, nr_z, ic_beta, ac, num_dims)
                            end if
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_mthinc_normals

    !> @brief Applies THINC (int_comp=1) or MTHINC (int_comp=2) interface compression to sharpen volume-fraction and density
    !! reconstructions at material interfaces. Called after WENO/MUSCL reconstruction per direction.
    subroutine s_thinc_compression(v_rs_ws, vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, recon_dir, &
                                   & is1_d, is2_d, is3_d)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in)    :: v_rs_ws
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL_rs_vf_x, vL_rs_vf_y, &
             & vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z
        integer, intent(in)               :: recon_dir
        type(int_bounds_info), intent(in) :: is1_d, is2_d, is3_d
        integer                           :: j, k, l, ix, iy, iz
        real(wp)                          :: aCL, aCR, aC, aTHINC, qmin, qmax, A, B, C
        real(wp)                          :: sgn, moncon, beta_eff
        real(wp)                          :: nh1, nh2, nh3, d_local, rho1, rho2
        real(wp)                          :: rho_b, rho_e

        #:for REC_DIR, XYZ, CC_PRI in [(1, 'x', 'x_cc'), (2, 'y', 'y_cc'), (3, 'z', 'z_cc')]
            if (recon_dir == ${REC_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, ix, iy, iz, aCL, aC, aCR, aTHINC, moncon, sgn, qmin, qmax, A, &
                                    & B, C, beta_eff, nh1, nh2, nh3, d_local, rho1, rho2, rho_b, rho_e]')
                do l = is3_d%beg, is3_d%end
                    do k = is2_d%beg, is2_d%end
                        do j = is1_d%beg, is1_d%end
                            aCL = v_rs_ws(j - 1, k, l, eqn_idx%adv%beg)
                            aC = v_rs_ws(j, k, l, eqn_idx%adv%beg)
                            aCR = v_rs_ws(j + 1, k, l, eqn_idx%adv%beg)

                            if (aC >= ic_eps .and. aC <= 1._wp - ic_eps) then
                                if (int_comp == 2 .and. n > 0) then  ! MTHINC
                                    ! Map reshaped (j,k,l) to physical (ix,iy,iz)
                                    #:if REC_DIR == 1
                                        ix = j; iy = k; iz = l
                                    #:elif REC_DIR == 2
                                        ix = k; iy = j; iz = l
                                    #:else
                                        ix = l; iy = k; iz = j
                                    #:endif

                                    nh1 = mthinc_nhat(1, ix, iy, iz)
                                    nh2 = mthinc_nhat(2, ix, iy, iz)
                                    nh3 = mthinc_nhat(3, ix, iy, iz)
                                    d_local = mthinc_d(ix, iy, iz)

                                    ! Skip if no valid normal was computed
                                    if (nh1*nh1 + nh2*nh2 + nh3*nh3 > 5e-1_wp) then
                                        rho1 = v_rs_ws(j, k, l, eqn_idx%cont%beg)/aC
                                        rho2 = v_rs_ws(j, k, l, eqn_idx%cont%end)/(1._wp - aC)

                                        ! Left face
                                        aTHINC = f_mthinc_face_average(nh1, nh2, nh3, d_local, ic_beta, ${REC_DIR}$, -5e-1_wp, &
                                                                       & num_dims)
                                        if (aTHINC < ic_eps) aTHINC = ic_eps
                                        if (aTHINC > 1._wp - ic_eps) aTHINC = 1._wp - ic_eps
                                        vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%beg) = rho1*aTHINC
                                        vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%end) = rho2*(1._wp - aTHINC)
                                        vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%beg) = aTHINC
                                        vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%end) = 1._wp - aTHINC

                                        ! Right face
                                        aTHINC = f_mthinc_face_average(nh1, nh2, nh3, d_local, ic_beta, ${REC_DIR}$, 5e-1_wp, &
                                                                       & num_dims)
                                        if (aTHINC < ic_eps) aTHINC = ic_eps
                                        if (aTHINC > 1._wp - ic_eps) aTHINC = 1._wp - ic_eps
                                        vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%beg) = rho1*aTHINC
                                        vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%end) = rho2*(1._wp - aTHINC)
                                        vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%beg) = aTHINC
                                        vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%end) = 1._wp - aTHINC
                                    end if
                                else  ! THINC
                                    moncon = (aCR - aC)*(aC - aCL)

                                    if (moncon > moncon_cutoff) then
                                        if (aCR - aCL > 0._wp) then
                                            sgn = 1._wp
                                        else
                                            sgn = -1._wp
                                        end if

                                        beta_eff = ic_beta

                                        qmin = min(aCR, aCL)
                                        qmax = max(aCR, aCL) - qmin

                                        C = (aC - qmin + sgm_eps)/(qmax + sgm_eps)
                                        B = exp(sgn*beta_eff*(2._wp*C - 1._wp))
                                        A = (B/cosh(beta_eff) - 1._wp)/tanh(beta_eff)

                                        rho_b = v_rs_ws(j, k, l, eqn_idx%cont%beg)/aC
                                        rho_e = v_rs_ws(j, k, l, eqn_idx%cont%end)/(1._wp - aC)

                                        ! Left face
                                        aTHINC = qmin + 5e-1_wp*qmax*(1._wp + sgn*A)
                                        if (aTHINC < ic_eps) aTHINC = ic_eps
                                        if (aTHINC > 1._wp - ic_eps) aTHINC = 1._wp - ic_eps
                                        vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%beg) = rho_b*aTHINC
                                        vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%end) = rho_e*(1._wp - aTHINC)
                                        vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%beg) = aTHINC
                                        vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%end) = 1._wp - aTHINC

                                        ! Right face
                                        aTHINC = qmin + 5e-1_wp*qmax*(1._wp + sgn*(tanh(beta_eff) + A)/(1._wp + A*tanh(beta_eff)))
                                        if (aTHINC < ic_eps) aTHINC = ic_eps
                                        if (aTHINC > 1._wp - ic_eps) aTHINC = 1._wp - ic_eps
                                        vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%beg) = rho_b*aTHINC
                                        vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%end) = rho_e*(1._wp - aTHINC)
                                        vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%beg) = aTHINC
                                        vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%end) = 1._wp - aTHINC
                                    end if
                                end if
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_thinc_compression

    subroutine s_finalize_thinc_module()

        if (int_comp == 2) then
            @:DEALLOCATE(mthinc_nhat)
            @:DEALLOCATE(mthinc_d)
        end if

    end subroutine s_finalize_thinc_module

end module m_thinc
