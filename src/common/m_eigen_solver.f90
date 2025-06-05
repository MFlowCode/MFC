!>
!! @file m_eigen_solver.f90
!! @brief Contains module m_eigen_solver

!> @brief The purpose of the module is to solve an eigenvalue problem
!!              for a complex general matrix. Subroutines are imported
!!              from EISPACK (https://netlib.org/eispack/) with minor
!!              modifications for compatibility.
module m_eigen_solver

    use m_precision_select

    implicit none

    private; 
    public :: cg, cbal, corth, comqr2, csroot, cdiv, pythag

contains

    !>  This subroutine calls the recommended sequence of subroutines from the
        !!              eigensystem subroutine package (eispack) to find the
        !!              eigenvalues and eigenvectors (if desired) of a complex
        !!              general matrix.
        !! @param nm the row dimension of the two-dimensional array parameters
        !! @param nl the order of the matrix a=(ar,ai)
        !! @param ar the real part of the complex general matrix
        !! @param ai the imaginary part of the complex general matrix
        !! @param wr the real part of the eigenvalues
        !! @param wi the imaginary part of the eigenvalues
        !! @param zr the real part of the eigenvectors
        !! @param zi the imaginary part of the eigenvectors
        !! @param fv1 temporary storage array
        !! @param fv2 temporary storage array
        !! @param fv3 temporary storage array
        !! @param ierr an error completion code
    pure subroutine cg(nm, nl, ar, ai, wr, wi, zr, zi, fv1, fv2, fv3, ierr)
        integer, intent(in) :: nm, nl
        real(wp), dimension(nm, nl), intent(inout) :: ar, ai
        real(wp), dimension(nl), intent(out) :: wr, wi
        real(wp), dimension(nm, nl), intent(out) :: zr, zi
        real(wp), dimension(nl), intent(out) :: fv1, fv2, fv3
        integer, intent(out) :: ierr

        integer :: is1, is2

        if (nl <= nm) go to 10
        ierr = 10*nl
        go to 50

10      call cbal(nm, nl, ar, ai, is1, is2, fv1)
        call corth(nm, nl, is1, is2, ar, ai, fv2, fv3)
        call comqr2(nm, nl, is1, is2, fv2, fv3, ar, ai, wr, wi, zr, zi, ierr)
        if (ierr /= 0) go to 50
        call cbabk2(nm, nl, is1, is2, fv1, nl, zr, zi)

50      return
    end subroutine cg

    !>  This subroutine is a translation of the algol procedure cbalance,
        !!              which is a complex version of balance, num. math. 13,
        !!              293-304(1969) by parlett and reinsch. handbook for auto.
        !!              comp., vol.ii-linear algebra, 315-326(1971).
        !!              This subroutine balances a complex matrix and isolates
        !!              eigenvalues whenever possible.
        !! @param nm the row dimension of the two-dimensional array parameters
        !! @param nl the order of the matrix
        !! @param ar the real part of the complex matrix to be balanced
        !! @param ai the imaginary part of the complex matrix to be balanced
        !! @param low one of two integers such that ar(i,j) and ai(i,j)
        !!            are equal to zero if
        !!            (1) i is greater than j and
        !!            (2) j=1, ,low-1 or i=igh+1,  ,nl.
        !! @param igh one of two integers such that ar(i,j) and ai(i,j)
        !!            are equal to zero if
        !!            (1) i is greater than j and
        !!            (2) j=1, ,low-1 or i=igh+1, ,nl.
        !! @param scale the information determining the permutations and scaling
        !!              factors used.
    pure subroutine cbal(nm, nl, ar, ai, low, igh, scale)
        integer, intent(in) :: nm, nl
        real(wp), dimension(nm, nl), intent(inout) :: ar, ai
        integer, intent(out) :: low, igh
        real(wp), dimension(nl), intent(out) :: scale

        integer :: i, j, k, l, ml, jj, iexc
        real(wp) :: c, f, g, r, s, b2, radix
        logical :: noconv

        radix = 16.0_wp

        b2 = radix*radix
        k = 1
        l = nl
        go to 100
! in-line procedure for row and column exchange
20      scale(ml) = j
        if (j == ml) go to 50

        do 30 i = 1, l
            f = ar(i, j)
            ar(i, j) = ar(i, ml)
            ar(i, ml) = f
            f = ai(i, j)
            ai(i, j) = ai(i, ml)
            ai(i, ml) = f
30      end do

        do 40 i = k, nl
            f = ar(j, i)
            ar(j, i) = ar(ml, i)
            ar(ml, i) = f
            f = ai(j, i)
            ai(j, i) = ai(ml, i)
            ai(ml, i) = f
40      end do

50      go to(80, 130), iexc
! search for rows isolating an eigenvalue and push them down
80      if (l == 1) go to 280
        l = l - 1
! for j=l step -1 until 1 do
100     do 120 jj = 1, l
            j = l + 1 - jj

            do 110 i = 1, l
                if (i == j) go to 110
                if (ar(j, i) /= 0.0_wp .or. ai(j, i) /= 0.0_wp) go to 120
110         end do

            ml = l
            iexc = 1
            go to 20
120     end do

        go to 140
! search for columns isolating an eigenvalue and push them left
130     k = k + 1

140     do 170 j = k, l

            do 150 i = k, l
                if (i == j) go to 150
                if (ar(i, j) /= 0.0_wp .or. ai(i, j) /= 0.0_wp) go to 170
150         end do

            ml = k
            iexc = 2
            go to 20
170     end do
! now balance the submatrix in rows k to l
        do 180 i = k, l
            scale(i) = 1.0_wp
180     end do
! iterative loop for norm reduction
190     noconv = .false.

        do 270 i = k, l
            c = 0.0_wp
            r = 0.0_wp

            do 200 j = k, l
                if (j == i) go to 200
                c = c + abs(ar(j, i)) + abs(ai(j, i))
                r = r + abs(ar(i, j)) + abs(ai(i, j))
200         end do
!     guard against zero c or r due to underflow
            if (c == 0.0_wp .or. r == 0.0_wp) go to 270
            g = r/radix
            f = 1.0_wp
            s = c + r
210         if (c >= g) go to 220
            f = f*radix
            c = c*b2
            go to 210
220         g = r*radix
230         if (c < g) go to 240
            f = f/radix
            c = c/b2
            go to 230
!     now balance
240         if ((c + r)/f >= 0.95_wp*s) go to 270
            g = 1.0_wp/f
            scale(i) = scale(i)*f
            noconv = .true.

            do 250 j = k, nl
                ar(i, j) = ar(i, j)*g
                ai(i, j) = ai(i, j)*g
250         end do

            do 260 j = 1, l
                ar(j, i) = ar(j, i)*f
                ai(j, i) = ai(j, i)*f
260         end do

270     end do

        if (noconv) go to 190

280     low = k
        igh = l
        return
    end subroutine cbal

    !>  This subroutine is a translation of a complex analogue of the algol
        !!              procedure orthes, num. math. 12, 349-368(1968) by martin
        !!              and wilkinson. handbook for auto. comp., vol.ii-linear
        !!              algebra, 339-358(1971). Given a complex general matrix,
        !!              this subroutine reduces a submatrix situated in rows and
        !!              columns low through igh to upper hessenberg form by
        !!              unitary similarity transformations.
        !! @param nm the row dimension of the two-dimensional array parameters
        !! @param nl the order of the matrix
        !! @param ar the real part of the complex matrix
        !! @param ai the imaginary part of the complex matrix
        !! @param low an integer determined by the balancing subroutine cbal.
        !!            if  cbal  has not been used, set low=1.
        !! @param igh an integer determined by the balancing subroutine cbal.
        !!            if  cbal  has not been used, set igh=nl.
        !! @param ortr further information about the transformations
        !! @param orti further information about the transformations
    pure subroutine corth(nm, nl, low, igh, ar, ai, ortr, orti)
        integer, intent(in) :: nm, nl, low, igh
        real(wp), dimension(nm, nl), intent(inout) :: ar, ai
        real(wp), dimension(igh), intent(out) :: ortr, orti

        integer :: i, j, ml, ii, jj, la, mp, kp1, mll
        real(wp) :: f, g, h, fi, fr, scale

        mll = 6

        la = igh - 1
        kp1 = low + 1
        if (la < kp1) go to 200

        do 180 ml = kp1, la
            h = 0.0_wp
            ortr(ml) = 0.0_wp
            orti(ml) = 0.0_wp
            scale = 0.0_wp
!     scale column (algol tol then not needed)
            do 90 i = ml, igh
                scale = scale + abs(ar(i, ml - 1)) + abs(ai(i, ml - 1))
90          end do
            if (scale == 0._wp) go to 180
            mp = ml + igh
!     for i=igh step -1 until ml do
            do 100 ii = ml, igh
                i = mp - ii
                ortr(i) = ar(i, ml - 1)/scale
                orti(i) = ai(i, ml - 1)/scale
                h = h + ortr(i)*ortr(i) + orti(i)*orti(i)
100         end do

            g = sqrt(h)
            call pythag(ortr(ml), orti(ml), f)
            if (f == 0._wp) go to 103
            h = h + f*g
            g = g/f
            ortr(ml) = (1.0_wp + g)*ortr(ml)
            orti(ml) = (1.0_wp + g)*orti(ml)
            go to 105

103         ortr(ml) = g
            ar(ml, ml - 1) = scale
!     form (i-(u*ut)/h) * a
105         do 130 j = ml, nl
                fr = 0.0_wp
                fi = 0.0_wp
!     for i=igh step -1 until ml do
                do 110 ii = ml, igh
                    i = mp - ii
                    fr = fr + ortr(i)*ar(i, j) + orti(i)*ai(i, j)
                    fi = fi + ortr(i)*ai(i, j) - orti(i)*ar(i, j)
110             end do

                fr = fr/h
                fi = fi/h

                do 120 i = ml, igh
                    ar(i, j) = ar(i, j) - fr*ortr(i) + fi*orti(i)
                    ai(i, j) = ai(i, j) - fr*orti(i) - fi*ortr(i)
120             end do

130         end do

!     form (i-(u*ut)/h)*a*(i-(u*ut)/h)
            do 160 i = 1, igh
                fr = 0.0_wp
                fi = 0.0_wp

!     for j=igh step -1 until ml do
                do 140 jj = ml, igh
                    j = mp - jj
                    fr = fr + ortr(j)*ar(i, j) - orti(j)*ai(i, j)
                    fi = fi + ortr(j)*ai(i, j) + orti(j)*ar(i, j)
140             end do

                fr = fr/h
                fi = fi/h

                do 150 j = ml, igh
                    ar(i, j) = ar(i, j) - fr*ortr(j) - fi*orti(j)
                    ai(i, j) = ai(i, j) + fr*orti(j) - fi*ortr(j)
150             end do

160         end do

            ortr(ml) = scale*ortr(ml)
            orti(ml) = scale*orti(ml)
            ar(ml, ml - 1) = -g*ar(ml, ml - 1)
            ai(ml, ml - 1) = -g*ai(ml, ml - 1)
180     end do

200     return
    end subroutine corth

    !>  This subroutine is a translation of a unitary analogue of the algol
        !!              procedure  comlr2, num. math. 16, 181-204(1970) by
        !!              peters and wilkinson. handbook for auto. comp.,
        !!              vol.ii-linear algebra, 372-395(1971). The unitary
        !!              analogue substitutes the qr algorithm of francis (comp.
        !!              jour. 4, 332-345(1962)) for the lr algorithm.
        !!              This subroutine finds the eigenvalues and eigenvectors
        !!              of a complex upper hessenberg matrix by the qr method.
        !!              The eigenvectors of a complex general matrix can
        !!              also be found if  corth  has been used to reduce this
        !!              general matrix to hessenberg form.
        !! @param nm the row dimension of the two-dimensional array parameters
        !! @param nl the order of the matrix
        !! @param low an integer determined by the balancing subroutine cbal.
        !!            if  cbal  has not been used, set low=1.
        !! @param igh an integer determined by the balancing subroutine cbal.
        !!            if  cbal  has not been used, set igh=nl.
        !! @param ortr information about the unitary transformations used in the
        !!             reduction by  corth
        !! @param orti information about the unitary transformations used in the
        !!             reduction by  corth
        !! @param hr the real part of the complex upper hessenberg matrix.
        !! @param hi the imaginary part of the complex upper hessenberg matrix.
        !! @param wr the real part of the eigenvalues
        !! @param wi the imaginary part of the eigenvalues
        !! @param zr the real part of the eigenvectors
        !! @param zi the imaginary part of the eigenvectors
        !! @param ierr an error completion code
    pure subroutine comqr2(nm, nl, low, igh, ortr, orti, hr, hi, wr, wi, zr, zi, ierr)
        integer, intent(in) :: nm, nl, low, igh
        real(wp), dimension(nm, nl), intent(inout) :: hr, hi
        real(wp), dimension(nl), intent(out) :: wr, wi
        real(wp), dimension(nm, nl), intent(out) :: zr, zi
        real(wp), dimension(igh), intent(inout) :: ortr, orti
        integer, intent(out) :: ierr

        integer :: i, j, k, l, ml, en, ii, jj, ll, nn, ip1, itn, its, lp1, enm1, iend
        real(wp) :: si, sr, ti, tr, xi, xr, xxi, xxr, yi, yr, zzi, zzr, &
                    norm, tst1, tst2, c
!
        ierr = 0
!     initialize eigenvector matrix
        do 101 j = 1, nl
!
            do 100 i = 1, nl
                zr(i, j) = 0.0_wp
                zi(i, j) = 0.0_wp
100         end do
            zr(j, j) = 1.0_wp
101     end do
!     form the matrix of accumulated transformations
!                from the information left by corth
        iend = igh - low - 1
        if (iend < 0) go to 180
        if (iend == 0) go to 150
        if (iend > 0) go to 105
!     for i=igh-1 step -1 until low+1 do
105     do 140 ii = 1, iend
            i = igh - ii
            if (abs(ortr(i)) == 0._wp .and. abs(orti(i)) == 0._wp) go to 140
            if (abs(hr(i, i - 1)) == 0._wp .and. abs(hi(i, i - 1)) == 0._wp) go to 140
!     norm below is negative of h formed in corth
            norm = hr(i, i - 1)*ortr(i) + hi(i, i - 1)*orti(i)
            ip1 = i + 1

            do 110 k = ip1, igh
                ortr(k) = hr(k, i - 1)
                orti(k) = hi(k, i - 1)
110         end do

            do 130 j = i, igh
                sr = 0.0_wp
                si = 0.0_wp

                do 115 k = i, igh
                    sr = sr + ortr(k)*zr(k, j) + orti(k)*zi(k, j)
                    si = si + ortr(k)*zi(k, j) - orti(k)*zr(k, j)
115             end do

                sr = sr/norm
                si = si/norm

                do 120 k = i, igh
                    zr(k, j) = zr(k, j) + sr*ortr(k) - si*orti(k)
                    zi(k, j) = zi(k, j) + sr*orti(k) + si*ortr(k)
120             end do

130         end do

140     end do
!     create real subdiagonal elements
150     l = low + 1

        do 170 i = l, igh
            ll = min0(i + 1, igh)
            if (abs(hi(i, i - 1)) == 0._wp) go to 170
            call pythag(hr(i, i - 1), hi(i, i - 1), norm)
            yr = hr(i, i - 1)/norm
            yi = hi(i, i - 1)/norm
            hr(i, i - 1) = norm
            hi(i, i - 1) = 0.0_wp

            do 155 j = i, nl
                si = yr*hi(i, j) - yi*hr(i, j)
                hr(i, j) = yr*hr(i, j) + yi*hi(i, j)
                hi(i, j) = si
155         end do

            do 160 j = 1, ll
                si = yr*hi(j, i) + yi*hr(j, i)
                hr(j, i) = yr*hr(j, i) - yi*hi(j, i)
                hi(j, i) = si
160         end do

            do 165 j = low, igh
                si = yr*zi(j, i) + yi*zr(j, i)
                zr(j, i) = yr*zr(j, i) - yi*zi(j, i)
                zi(j, i) = si
165         end do
170     end do
!     store roots isolated by cbal
180     do 200 i = 1, nl
            if (i >= low .and. i <= igh) go to 200
            wr(i) = hr(i, i)
            wi(i) = hi(i, i)
200     end do

        en = igh
        tr = 0.0_wp
        ti = 0.0_wp
        itn = 30*nl
!     search for next eigenvalue
220     if (en < low) go to 680
        its = 0
        enm1 = en - 1
!     look for single small sub-diagonal element
!                for l=en step -1 until low do
240     do 260 ll = low, en
            l = en + low - ll
            if (l == low) go to 300
            tst1 = abs(hr(l - 1, l - 1)) + abs(hi(l - 1, l - 1)) &
                   + abs(hr(l, l)) + abs(hi(l, l))
            tst2 = tst1 + abs(hr(l, l - 1))
            if (tst2 == tst1) go to 300
260     end do
!     form shift
300     if (l == en) go to 660
        if (itn == 0) go to 1000
        if (its == 10 .or. its == 20) go to 320
        sr = hr(en, en)
        si = hi(en, en)
        xr = hr(enm1, en)*hr(en, enm1)
        xi = hi(enm1, en)*hr(en, enm1)
        if (xr == 0.0_wp .and. xi == 0.0_wp) go to 340
        yr = (hr(enm1, enm1) - sr)/2.0_wp
        yi = (hi(enm1, enm1) - si)/2.0_wp
        call csroot(yr**2 - yi**2 + xr, 2.0_wp*yr*yi + xi, zzr, zzi)
        if (yr*zzr + yi*zzi >= 0.0_wp) go to 310
        zzr = -zzr
        zzi = -zzi
310     call cdiv(xr, xi, yr + zzr, yi + zzi, xxr, xxi)
        sr = sr - xxr
        si = si - xxi
        go to 340
!     form exceptional shift
320     sr = abs(hr(en, enm1)) + abs(hr(enm1, en - 2))
        si = 0.0_wp

340     do 360 i = low, en
            hr(i, i) = hr(i, i) - sr
            hi(i, i) = hi(i, i) - si
360     end do

        tr = tr + sr
        ti = ti + si
        its = its + 1
        itn = itn - 1
!     reduce to triangle (rows)
        lp1 = l + 1

        do 500 i = lp1, en
            sr = hr(i, i - 1)
            hr(i, i - 1) = 0.0_wp
            call pythag(hr(i - 1, i - 1), hi(i - 1, i - 1), c)
            call pythag(c, sr, norm)
            xr = hr(i - 1, i - 1)/norm
            wr(i - 1) = xr
            xi = hi(i - 1, i - 1)/norm
            wi(i - 1) = xi
            hr(i - 1, i - 1) = norm
            hi(i - 1, i - 1) = 0.0_wp
            hi(i, i - 1) = sr/norm

            do 490 j = i, nl
                yr = hr(i - 1, j)
                yi = hi(i - 1, j)
                zzr = hr(i, j)
                zzi = hi(i, j)
                hr(i - 1, j) = xr*yr + xi*yi + hi(i, i - 1)*zzr
                hi(i - 1, j) = xr*yi - xi*yr + hi(i, i - 1)*zzi
                hr(i, j) = xr*zzr - xi*zzi - hi(i, i - 1)*yr
                hi(i, j) = xr*zzi + xi*zzr - hi(i, i - 1)*yi
490         end do

500     end do

        si = hi(en, en)
        if (abs(si) == 0._wp) go to 540
        call pythag(hr(en, en), si, norm)
        sr = hr(en, en)/norm
        si = si/norm
        hr(en, en) = norm
        hi(en, en) = 0.0_wp
        if (en == nl) go to 540
        ip1 = en + 1

        do 520 j = ip1, nl
            yr = hr(en, j)
            yi = hi(en, j)
            hr(en, j) = sr*yr + si*yi
            hi(en, j) = sr*yi - si*yr
520     end do
!     inverse operation (columns)
540     do 600 j = lp1, en
            xr = wr(j - 1)
            xi = wi(j - 1)

            do 580 i = 1, j
                yr = hr(i, j - 1)
                yi = 0.0_wp
                zzr = hr(i, j)
                zzi = hi(i, j)
                if (i == j) go to 560
                yi = hi(i, j - 1)
                hi(i, j - 1) = xr*yi + xi*yr + hi(j, j - 1)*zzi
560             hr(i, j - 1) = xr*yr - xi*yi + hi(j, j - 1)*zzr
                hr(i, j) = xr*zzr + xi*zzi - hi(j, j - 1)*yr
                hi(i, j) = xr*zzi - xi*zzr - hi(j, j - 1)*yi
580         end do

            do 590 i = low, igh
                yr = zr(i, j - 1)
                yi = zi(i, j - 1)
                zzr = zr(i, j)
                zzi = zi(i, j)
                zr(i, j - 1) = xr*yr - xi*yi + hi(j, j - 1)*zzr
                zi(i, j - 1) = xr*yi + xi*yr + hi(j, j - 1)*zzi
                zr(i, j) = xr*zzr + xi*zzi - hi(j, j - 1)*yr
                zi(i, j) = xr*zzi - xi*zzr - hi(j, j - 1)*yi
590         end do
600     end do

        if (abs(si) == 0._wp) go to 240

        do 630 i = 1, en
            yr = hr(i, en)
            yi = hi(i, en)
            hr(i, en) = sr*yr - si*yi
            hi(i, en) = sr*yi + si*yr
630     end do

        do 640 i = low, igh
            yr = zr(i, en)
            yi = zi(i, en)
            zr(i, en) = sr*yr - si*yi
            zi(i, en) = sr*yi + si*yr
640     end do

        go to 240
!     a root found
660     hr(en, en) = hr(en, en) + tr
        wr(en) = hr(en, en)
        hi(en, en) = hi(en, en) + ti
        wi(en) = hi(en, en)
        en = enm1
        go to 220
!     all roots found.  backsubstitute to find
!                vectors of upper triangular form
680     norm = 0.0_wp

        do i = 1, nl
            do j = i, nl
                tr = abs(hr(i, j)) + abs(hi(i, j))
                if (tr > norm) norm = tr
            end do
        end do

        if (nl == 1 .or. norm == 0._wp) go to 1001
!     for en=nl step -1 until 2 do
        do 800 nn = 2, nl
            en = nl + 2 - nn
            xr = wr(en)
            xi = wi(en)
            hr(en, en) = 1.0_wp
            hi(en, en) = 0.0_wp
            enm1 = en - 1
!     for i=en-1 step -1 until 1 do
            do 780 ii = 1, enm1
                i = en - ii
                zzr = 0.0_wp
                zzi = 0.0_wp
                ip1 = i + 1

                do 740 j = ip1, en
                    zzr = zzr + hr(i, j)*hr(j, en) - hi(i, j)*hi(j, en)
                    zzi = zzi + hr(i, j)*hi(j, en) + hi(i, j)*hr(j, en)
740             end do

                yr = xr - wr(i)
                yi = xi - wi(i)
                if (yr /= 0.0_wp .or. yi /= 0.0_wp) go to 765
                tst1 = norm
                yr = tst1
760             yr = 0.01_wp*yr
                tst2 = norm + yr
                if (tst2 > tst1) go to 760
765             continue
                call cdiv(zzr, zzi, yr, yi, hr(i, en), hi(i, en))
!     overflow control
                tr = abs(hr(i, en)) + abs(hi(i, en))
                if (tr == 0.0_wp) go to 780
                tst1 = tr
                tst2 = tst1 + 1.0_wp/tst1
                if (tst2 > tst1) go to 780
                do 770 j = i, en
                    hr(j, en) = hr(j, en)/tr
                    hi(j, en) = hi(j, en)/tr
770             end do

780         end do

800     end do
!     end backsubstitution
!     vectors of isolated roots
        do 840 i = 1, nl
            if (i >= low .and. i <= igh) go to 840

            do 820 j = I, nl
                zr(i, j) = hr(i, j)
                zi(i, j) = hi(i, j)
820         end do

840     end do
!     multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=nl step -1 until low do
        do jj = low, nl
            j = nl + low - jj
            ml = min0(j, igh)

            do i = low, igh
                zzr = 0.0_wp
                zzi = 0.0_wp

                do 860 k = low, ml
                    zzr = zzr + zr(i, k)*hr(k, j) - zi(i, k)*hi(k, j)
                    zzi = zzi + zr(i, k)*hi(k, j) + zi(i, k)*hr(k, j)
860             end do

                zr(i, j) = zzr
                zi(i, j) = zzi
            end do
        end do

        go to 1001
!     set error. all eigenvalues have not
!                converged after 30*nl iterations
1000    ierr = en
1001    return
    end subroutine comqr2

    !>  This subroutine is a translation of the algol procedure cbabk2, which is
        !!              a complex version of balbak, num. math. 13,
        !!              293-304(1969) by parlett and reinsch. handbook for auto.
        !!              comp., vol.ii-linear algebra, 315-326(1971).
        !!              This subroutine forms the eigenvectors of a complex
        !!              general matrix by back transforming those of the
        !!              correspondingbalanced matrix determined by cbal.
        !! @param nm the row dimension of the two-dimensional array parameters
        !! @param nl the order of the matrix
        !! @param ar the real part of the complex matrix to be balanced
        !! @param ai the imaginary part of the complex matrix to be balanced
        !! @param low an integer determined by the balancing subroutine cbal
        !! @param igh an integer determined by the balancing subroutine cbal
        !! @param scale the information determining the permutations and scaling
        !!              factors used.
        !! @param ml the number of eigenvectors to be back transformed
        !! @param zr the real part of the eigenvectors to be back transformed in
        !!           their first ml columns
        !! @param zi the imaginary part of the eigenvectors to be back
        !!           transformed in their first ml columns
    pure subroutine cbabk2(nm, nl, low, igh, scale, ml, zr, zi)
        integer, intent(in) :: nm, nl, low, igh
        real(wp), intent(in) :: scale(nl)
        integer, intent(in) :: ml
        real(wp), intent(inout) :: zr(nm, ml), zi(nm, ml)

        integer :: i, j, k, ii
        real(wp) :: s

        if (ml == 0) go to 200
        if (igh == low) go to 120

        do 110 i = low, igh
            s = scale(i)
!     left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.0_wp/scale(i).
            do 100 j = 1, ml
                zr(i, j) = zr(i, j)*s
                zi(i, j) = zi(i, j)*s
100         end do

110     end do
!     for i=low-1 step -1 until 1,
!     igh+1 step 1 until nl do
120     do 140 ii = 1, nl
            i = ii
            if (i >= low .and. i <= igh) go to 140
            if (i < low) i = low - ii
            k = scale(i)
            if (k == i) go to 140

            do 130 j = 1, ml
                s = zr(i, j)
                zr(i, j) = zr(k, j)
                zr(k, j) = s
                s = zi(i, j)
                zi(i, j) = zi(k, j)
                zi(k, j) = s
130         end do

140     end do

200     return
    end subroutine cbabk2

    pure elemental subroutine csroot(xr, xi, yr, yi)
        real(wp), intent(in) :: xr, xi
        real(wp), intent(out) :: yr, yi

!     (yr,yi) = complex sqrt(xr,xi)
!     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)

        real(wp) :: s, tr, ti, c
        tr = xr
        ti = xi
        call pythag(tr, ti, c)
        s = sqrt(0.5_wp*(c + abs(tr)))
        if (tr >= 0.0_wp) yr = s
        if (ti < 0.0_wp) s = -s
        if (tr <= 0.0_wp) yi = s
        if (tr < 0.0_wp) yr = 0.5_wp*(ti/yi)
        if (tr > 0.0_wp) yi = 0.5_wp*(ti/yr)
        return
    end subroutine csroot

    pure elemental subroutine cdiv(ar, ai, br, bi, cr, ci)
        real(wp), intent(in) :: ar, ai, br, bi
        real(wp), intent(out) :: cr, ci
        real(wp) :: s, ars, ais, brs, bis

        s = abs(br) + abs(bi)
        ars = ar/s
        ais = ai/s
        brs = br/s
        bis = bi/s
        s = brs**2._wp + bis**2._wp
        cr = (ars*brs + ais*bis)/s
        ci = (ais*brs - ars*bis)/s
        return
    end subroutine cdiv

    pure elemental subroutine pythag(a, b, c)
        real(wp), intent(in) :: a, b
        real(wp), intent(out) :: c

!     finds sqrt(a**2+b**2) without overflow or destructive underflow

        real(wp) :: p, r, s, t, u
        p = max(abs(a), abs(b))
        if (p == 0.0_wp) go to 20
        r = (min(abs(a), abs(b))/p)**2
10      continue
        t = 4.0_wp + r
        if (t == 4.0_wp) go to 20
        s = r/t
        u = 1.0_wp + 2.0_wp*s
        p = u*p
        r = (s/u)**2._wp*r
        go to 10
20      c = p
        return
    end subroutine pythag

end module m_eigen_solver
