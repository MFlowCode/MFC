module m_simplex_noise

    use m_constants

    use m_precision_select

    implicit none

    private; public :: f_simplex3d, &
 f_simplex2d

    integer, parameter :: p_vec(0:511) = [ &
                          151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, &
                          69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, &
                          252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, &
                          171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, &
                          122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, &
                          161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, &
                          159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, &
                          147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, &
                          183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, &
                          129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, &
                          228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, &
                          239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, &
                          4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, &
                          61, 156, 180, &
                          151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, &
                          69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, &
                          252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, &
                          171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, &
                          122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, &
                          161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, &
                          159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, &
                          147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, &
                          183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, &
                          129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, &
                          228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, &
                          239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, &
                          4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, &
                          61, 156, 180]

    real(wp), parameter :: grad3(12, 3) = reshape([ &
                                                  1._wp, 1._wp, 0._wp, &
                                                  -1._wp, 1._wp, 0._wp, &
                                                  1._wp, -1._wp, 0._wp, &
                                                  -1._wp, -1._wp, 0._wp, &
                                                  1._wp, 0._wp, 1._wp, &
                                                  -1._wp, 0._wp, 1._wp, &
                                                  1._wp, 0._wp, -1._wp, &
                                                  -1._wp, 0._wp, -1._wp, &
                                                  0._wp, 1._wp, 1._wp, &
                                                  0._wp, -1._wp, 1._wp, &
                                                  0._wp, 1._wp, -1._wp, &
                                                  0._wp, -1._wp, -1._wp], shape=[12, 3])

    real(wp), parameter :: grad2(10, 2) = reshape([ &
                                                  1._wp, 1._wp, &
                                                  -1._wp, 1._wp, &
                                                  1._wp, -1._wp, &
                                                  -1._wp, -1._wp, &
                                                  1._wp, 0._wp, &
                                                  -1._wp, 0._wp, &
                                                  0._wp, 1._wp, &
                                                  0._wp, -1._wp, &
                                                  1._wp, 1._wp, &
                                                  -1._wp, 1._wp], shape=[10, 2])

contains

    function f_simplex3d(xin, yin, zin) result(n)

        real(wp), intent(in) :: xin, yin, zin
        real(wp) :: n
        real(wp) :: n0, n1, n2, n3
        real(wp) :: f3, g3
        real(wp) :: x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3
        integer :: i, j, k, i1, j1, k1, i2, j2, k2
        integer :: ii, jj, kk, gi0, gi1, gi2, gi3
        real(wp) :: s, t, r, t0, t1, t2, t3
        real(wp) :: g(3)
        real(wp) :: x, y, z

        f3 = 1._wp/3._wp
        g3 = 1._wp/6._wp

        s = (xin + yin + zin)*f3
        i = floor(xin + s)
        j = floor(yin + s)
        k = floor(zin + s)

        t = (i + j + k)*g3

        x0 = xin - (i - t)
        y0 = yin - (j - t)
        z0 = zin - (k - t)

        if (x0 >= y0) then
            if (y0 >= z0) then
                i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0
            else if (x0 >= z0) then
                i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1
            else
                i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1
            end if
        else
            if (y0 < z0) then
                i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1
            else if (x0 < z0) then
                i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1
            else
                i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0
            end if
        end if

        x1 = x0 - i1 + g3
        y1 = y0 - j1 + g3
        z1 = z0 - k1 + g3
        x2 = x0 - i2 + 2._wp*g3
        y2 = y0 - j2 + 2._wp*g3
        z2 = z0 - k2 + 2._wp*g3
        x3 = x0 - 1._wp + 3._wp*g3
        y3 = y0 - 1._wp + 3._wp*g3
        z3 = z0 - 1._wp + 3._wp*g3

        ii = iand(i, 255)
        jj = iand(j, 255)
        kk = iand(k, 255)

        gi0 = mod(p_vec(ii + p_vec(jj + p_vec(kk) + 1) + 1), 12) + 1
        gi1 = mod(p_vec(ii + i1 + p_vec(jj + j1 + p_vec(kk + k1) + 1) + 1), 12) + 1
        gi2 = mod(p_vec(ii + i2 + p_vec(jj + j2 + p_vec(kk + k2) + 1) + 1), 12) + 1
        gi3 = mod(p_vec(ii + 1 + p_vec(jj + 1 + p_vec(kk + 1) + 1) + 1), 12) + 1

        t0 = 0.5_wp - x0*x0 - y0*y0 - z0*z0
        if (t0 < 0._wp) then
            n0 = 0._wp
        else
            t0 = t0*t0
            n0 = t0*t0*dot_product(grad3(gi0, :), [x0, y0, z0])
        end if

        t1 = 0.5_wp - x1*x1 - y1*y1 - z1*z1
        if (t1 < 0._wp) then
            n1 = 0._wp
        else
            t1 = t1*t1
            n1 = t1*t1*dot_product(grad3(gi1, :), [x1, y1, z1])
        end if

        t2 = 0.5_wp - x2*x2 - y2*y2 - z2*z2
        if (t2 < 0._wp) then
            n2 = 0._wp
        else
            t2 = t2*t2
            n2 = t2*t2*dot_product(grad3(gi2, :), [x2, y2, z2])
        end if

        t3 = 0.5_wp - x3*x3 - y3*y3 - z3*z3
        if (t3 < 0._wp) then
            n3 = 0._wp
        else
            t3 = t3*t3
            n3 = t3*t3*dot_product(grad3(gi3, :), [x3, y3, z3])
        end if

        n = 32._wp*(n0 + n1 + n2 + n3)

    end function f_simplex3d

    function f_simplex2d(xin, yin) result(n)

        real(wp), intent(in) :: xin, yin
        real(wp) :: n
        real(wp), parameter :: F2 = 0.5_wp*(sqrt(3._wp) - 1._wp)
        real(wp), parameter :: G2 = (3._wp - sqrt(3._wp))/6._wp
        integer :: i, j, ii, jj, gi0, gi1, gi2
        real(wp) :: s, t, x0, y0, x1, y1, x2, y2
        real(wp) :: t0, t1, t2, n0, n1, n2
        integer :: i1, j1

        s = (xin + yin)*F2
        i = floor(xin + s)
        j = floor(yin + s)

        t = real(i + j, 8)*G2

        x0 = xin - (i - t)
        y0 = yin - (j - t)

        if (x0 > y0) then
            i1 = 1; j1 = 0
        else
            i1 = 0; j1 = 1
        end if

        x1 = x0 - i1 + G2
        y1 = y0 - j1 + G2
        x2 = x0 - 1._wp + 2._wp*G2
        y2 = y0 - 1._wp + 2._wp*G2

        ii = mod(i, 255)
        jj = mod(j, 255)

        gi0 = mod(p_vec(ii + p_vec(jj)), 10) + 1
        gi1 = mod(p_vec(ii + i1 + p_vec(jj + j1)), 10) + 1
        gi2 = mod(p_vec(ii + 1 + p_vec(jj + 1)), 10) + 1

        t0 = 0.5_wp - x0*x0 - y0*y0
        if (t0 < 0._wp) then
            n0 = 0._wp
        else
            t0 = t0*t0
            n0 = t0*t0*dot2(gi0, x0, y0)
        end if

        t1 = 0.5_wp - x1*x1 - y1*y1
        if (t1 < 0._wp) then
            n1 = 0._wp
        else
            t1 = t1*t1
            n1 = t1*t1*dot2(gi1, x1, y1)
        end if

        t2 = 0.5_wp - x2*x2 - y2*y2
        if (t2 < 0._wp) then
            n2 = 0._wp
        else
            t2 = t2*t2
            n2 = t2*t2*dot2(gi2, x2, y2)
        end if

        n = 70._wp*(n0 + n1 + n2)

    end function f_simplex2d

    function dot2(g, x, y) result(dot)

        integer, intent(in) :: g
        real(wp), intent(in) :: x, y
        real(wp) :: dot
        dot = grad2(g + 1, 1)*x + grad2(g + 1, 2)*y

    end function

end module m_simplex_noise
