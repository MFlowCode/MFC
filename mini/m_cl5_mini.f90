! Minimal reproducer for ifx 2025.1.1 SPIR-V ICE #5633
!
! ICE trigger: multiple !$omp target teams loop regions in one module, each
! calling a different !$omp declare target subroutine in the SAME module.
! Subroutines access allocatable module arrays and a derived type with
! array-of-arrays fields (interp_coeffs(2,2,2)) plus character field.
!
! Compile (expect ICE):
!   ifx -free -fiopenmp -fopenmp-targets=spir64 -O3 -c m_cl5_mini.f90
! Workaround (passes):
!   ifx -free -fiopenmp -fopenmp-targets=spir64 -O3 -fno-inline -c m_cl5_mini.f90

module m_cl5

    implicit none
    private

    integer, parameter :: wp = kind(1.0d0)

    type :: patch_t
        integer  :: geometry
        real(wp) :: radius, cx, cy, cz
        real(wp) :: lx, ly, lz
        real(wp), dimension(1:3,1:3) :: rot, rot_inv
        real(wp), dimension(1:3)     :: offset
        character(LEN=200) :: label    ! char field in device-mapped type
    end type patch_t

    type :: gp_t
        integer,  dimension(3)       :: loc
        real(wp), dimension(2, 2, 2) :: interp_coeffs  ! 3D array in struct
        integer                      :: pid
        real(wp)                     :: levelset
        real(wp), dimension(3)       :: levelset_norm
        integer                      :: xp, yp, zp
    end type gp_t

    type(patch_t), dimension(8), target :: patches
    real(wp), allocatable :: xc(:), yc(:), zc(:)
    !$omp declare target(patches, xc, yc, zc)

    public :: gp_t, s_init, s_finalize, s_dispatch

contains

    subroutine s_init(n)
        integer, intent(in) :: n
        integer :: i
        allocate(xc(n), yc(n), zc(n))
        do i = 1, n
            xc(i) = real(i,wp)*0.1_wp
            yc(i) = real(i,wp)*0.1_wp
            zc(i) = real(i,wp)*0.1_wp
        end do
        do i = 1, 8
            patches(i)%geometry = i
            patches(i)%radius   = 0.5_wp; patches(i)%cx = 0.5_wp
            patches(i)%cy       = 0.5_wp; patches(i)%cz = 0.5_wp
            patches(i)%lx       = 0.4_wp; patches(i)%ly = 0.4_wp; patches(i)%lz = 0.4_wp
            patches(i)%rot      = 0.0_wp
            patches(i)%rot(1,1) = 1.0_wp; patches(i)%rot(2,2) = 1.0_wp; patches(i)%rot(3,3) = 1.0_wp
            patches(i)%rot_inv  = patches(i)%rot
            patches(i)%offset   = 0.0_wp
            patches(i)%label    = 'none'
        end do
        !$omp target enter data map(to: patches, xc, yc, zc)
    end subroutine s_init

    subroutine s_finalize()
        !$omp target exit data map(delete: patches, xc, yc, zc)
        deallocate(xc, yc, zc)
    end subroutine s_finalize

    subroutine s_dispatch(gps, n)
        type(gp_t), intent(inout) :: gps(:)
        integer,    intent(in)    :: n
        integer :: i

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 1) call s_geo1(gps(i))
        end do

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 2) call s_geo2(gps(i))
        end do

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 3) call s_geo3(gps(i))
        end do

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 4) call s_geo4(gps(i))
        end do

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 5) call s_geo5(gps(i))
        end do

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 6) call s_geo6(gps(i))
        end do

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 7) call s_geo7(gps(i))
        end do

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 8) call s_geo8(gps(i))
        end do

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 9) call s_geo9(gps(i))
        end do

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            if (patches(gps(i)%pid)%geometry == 10) call s_geo10(gps(i))
        end do

    end subroutine s_dispatch

    subroutine s_geo1(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j, k
        real(wp), dimension(3) :: dv
        real(wp) :: d
        id = gp%pid; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dv(1) = xc(i) - patches(id)%cx - real(gp%xp,wp)*0.1_wp
        dv(2) = yc(j) - patches(id)%cy - real(gp%yp,wp)*0.1_wp
        dv(3) = zc(k) - patches(id)%cz - real(gp%zp,wp)*0.1_wp
        d = sqrt(sum(dv**2))
        gp%levelset = d - patches(id)%radius
        if (d > 0._wp) then; gp%levelset_norm = dv/d; else; gp%levelset_norm = 0._wp; end if
    end subroutine s_geo1

    subroutine s_geo2(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j, k
        real(wp) :: dx, dy, dz
        id = gp%pid; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dx = abs(xc(i) - patches(id)%cx) - patches(id)%lx/2._wp
        dy = abs(yc(j) - patches(id)%cy) - patches(id)%ly/2._wp
        dz = abs(zc(k) - patches(id)%cz) - patches(id)%lz/2._wp
        gp%levelset = max(dx, max(dy, dz))
        gp%levelset_norm(1) = sign(1._wp, xc(i) - patches(id)%cx)
        gp%levelset_norm(2) = sign(1._wp, yc(j) - patches(id)%cy)
        gp%levelset_norm(3) = sign(1._wp, zc(k) - patches(id)%cz)
    end subroutine s_geo2

    subroutine s_geo3(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j, k
        real(wp), dimension(3) :: dv, ldv
        real(wp) :: r2d
        id = gp%pid; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dv(1) = xc(i) - patches(id)%cx; dv(2) = yc(j) - patches(id)%cy; dv(3) = zc(k) - patches(id)%cz
        ldv = matmul(patches(id)%rot_inv, dv)
        r2d = sqrt(ldv(1)**2 + ldv(2)**2)
        gp%levelset = r2d - patches(id)%radius
        if (r2d > 0._wp) then
            gp%levelset_norm(1) = ldv(1)/r2d; gp%levelset_norm(2) = ldv(2)/r2d
        else
            gp%levelset_norm(1) = 0._wp; gp%levelset_norm(2) = 0._wp
        end if
        gp%levelset_norm(3) = 0._wp
    end subroutine s_geo3

    subroutine s_geo4(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j, k
        real(wp), dimension(3) :: dv, ldv
        real(wp) :: d
        id = gp%pid; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dv(1) = xc(i) - patches(id)%cx; dv(2) = yc(j) - patches(id)%cy; dv(3) = zc(k) - patches(id)%cz
        ldv = matmul(patches(id)%rot_inv, dv) - patches(id)%offset
        d = sqrt(sum(ldv**2))
        gp%levelset = d - patches(id)%radius
        gp%levelset_norm = ldv / max(d, 1.0e-12_wp)
    end subroutine s_geo4

    subroutine s_geo5(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j
        real(wp) :: ax, ay, d
        real(wp), dimension(3) :: dv
        id = gp%pid; i = gp%loc(1); j = gp%loc(2)
        ax = patches(id)%lx/2._wp; ay = patches(id)%ly/2._wp
        dv(1) = (xc(i) - patches(id)%cx)/ax
        dv(2) = (yc(j) - patches(id)%cy)/ay
        dv(3) = 0._wp
        d = sqrt(dv(1)**2 + dv(2)**2)
        gp%levelset = d - 1._wp
        if (d > 0._wp) then; gp%levelset_norm = dv/d; else; gp%levelset_norm = 0._wp; end if
    end subroutine s_geo5

    subroutine s_geo6(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j
        real(wp) :: dx, dy
        id = gp%pid; i = gp%loc(1); j = gp%loc(2)
        dx = abs(xc(i) - patches(id)%cx) - patches(id)%lx/2._wp
        dy = abs(yc(j) - patches(id)%cy) - patches(id)%ly/2._wp
        gp%levelset = max(dx, dy)
        gp%levelset_norm(1) = sign(1._wp, xc(i) - patches(id)%cx)
        gp%levelset_norm(2) = sign(1._wp, yc(j) - patches(id)%cy)
        gp%levelset_norm(3) = 0._wp
    end subroutine s_geo6

    subroutine s_geo7(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j
        real(wp), dimension(3) :: dv, ldv
        real(wp) :: d
        id = gp%pid; i = gp%loc(1); j = gp%loc(2)
        dv(1) = xc(i) - patches(id)%cx; dv(2) = yc(j) - patches(id)%cy; dv(3) = 0._wp
        ldv = matmul(patches(id)%rot_inv, dv) - patches(id)%offset
        d = sqrt(sum(ldv(1:2)**2))
        gp%levelset = d - patches(id)%radius
        gp%levelset_norm = ldv / max(d, 1.0e-12_wp)
    end subroutine s_geo7

    subroutine s_geo8(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j
        real(wp) :: d1, d2, d3
        real(wp), dimension(3) :: dv
        id = gp%pid; i = gp%loc(1); j = gp%loc(2)
        dv(1) = xc(i) - patches(id)%cx; dv(2) = yc(j) - patches(id)%cy; dv(3) = 0._wp
        d1 = dv(2) + patches(id)%ly/2._wp
        d2 = -dv(2) + dv(1)*sqrt(3._wp)/3._wp + patches(id)%ly/3._wp
        d3 = -dv(2) - dv(1)*sqrt(3._wp)/3._wp + patches(id)%ly/3._wp
        gp%levelset = min(max(-d1,-d2), max(-d1,-d3))
        gp%levelset_norm(1) = dv(1)/max(sqrt(sum(dv(1:2)**2)),1e-12_wp)
        gp%levelset_norm(2) = dv(2)/max(sqrt(sum(dv(1:2)**2)),1e-12_wp)
        gp%levelset_norm(3) = 0._wp
    end subroutine s_geo8

    subroutine s_geo9(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j
        real(wp), dimension(3) :: dv
        real(wp) :: d
        id = gp%pid; i = gp%loc(1); j = gp%loc(2)
        dv(1) = xc(i) - patches(id)%cx + patches(id)%lx*patches(id)%offset(1)
        dv(2) = yc(j) - patches(id)%cy + patches(id)%ly*patches(id)%offset(2)
        dv(3) = 0._wp
        d = sqrt(sum(dv**2))
        gp%levelset = d - patches(id)%radius
        gp%levelset_norm = dv / max(d, 1.0e-12_wp)
    end subroutine s_geo9

    subroutine s_geo10(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j, k
        real(wp), dimension(3) :: dv, ldv
        real(wp) :: d
        id = gp%pid; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dv(1) = xc(i) - patches(id)%cx; dv(2) = yc(j) - patches(id)%cy; dv(3) = zc(k) - patches(id)%cz
        ldv = matmul(patches(id)%rot, dv)
        d = sqrt(sum(ldv**2))
        gp%levelset = d - patches(id)%radius
        gp%levelset_norm = matmul(patches(id)%rot_inv, ldv) / max(d, 1.0e-12_wp)
    end subroutine s_geo10

end module m_cl5

program test_cl5
    use m_cl5
    implicit none
    integer, parameter :: N = 16
    type(gp_t) :: gps(N)
    integer :: i
    call s_init(N)
    do i = 1, N
        gps(i)%levelset      = 0.0d0
        gps(i)%levelset_norm = 0.0d0
        gps(i)%interp_coeffs = 0.0d0
        gps(i)%loc           = [mod(i,N)+1, mod(i*2,N)+1, mod(i*3,N)+1]
        gps(i)%pid           = mod(i-1,8)+1
        gps(i)%xp            = 0; gps(i)%yp = 0; gps(i)%zp = 0
    end do
    call s_dispatch(gps, N)
    write(*,*) 'levelset(1) =', gps(1)%levelset
    call s_finalize()
end program test_cl5
