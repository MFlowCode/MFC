!> @brief Stable difference: log_cosh(a+h) - log_cosh(a-h) = 2*atanh(tanh(a)*tanh(h)). Avoids catastrophic cancellation when h is
!! small relative to a.
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

!> @brief Face-averaged THINC function at a cell face. face_dir: 1=x, 2=y, 3=z. face_pos: -0.5 (low) or +0.5 (high). The face
!! coordinate in face_dir is fixed; remaining directions are integrated over [-1/2, 1/2] (analytically along one, Gauss for others).
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
