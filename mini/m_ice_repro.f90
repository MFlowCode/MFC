! Minimal reproducer for ifx 2025.1.1 SPIR-V ICE #5633
!
! Trigger: matmul() called with a derived-type member (dimension(3,3) array)
!          inside an !$omp declare target subroutine, called from an
!          !$omp target teams loop. Module-level allocatable arrays are also
!          accessed from the declare-target routine.
!
! ifx version: 2025.1.1 20250418 (iimpi/2025a)
! Hardware:    Intel GPU Max 1100 (Ponte Vecchio XT)
!
! Compile (ICE -- error #5633: segmentation violation in SPIR-V backend):
!   ifx -free -fiopenmp -fopenmp-targets=spir64 -O3 -c m_ice_repro.f90
!   ifx -free -fiopenmp -fopenmp-targets=spir64 -O2 -c m_ice_repro.f90
!   ifx -free -fiopenmp -fopenmp-targets=spir64 -O1 -c m_ice_repro.f90
!
! Workaround (OK -- inlining suppressed so matmul is called via link, not inline):
!   ifx -free -fiopenmp -fopenmp-targets=spir64 -O3 -fno-inline -c m_ice_repro.f90

module m_ice

    implicit none
    private

    integer, parameter :: wp = kind(1.0d0)

    type :: patch_t
        real(wp), dimension(3, 3) :: mat  ! 3x3 member -- key field
        real(wp) :: cx, cy, cz, radius
    end type patch_t

    type :: gp_t
        integer  :: pid
        integer, dimension(3) :: loc
        real(wp) :: levelset
        real(wp), dimension(3) :: norm
    end type gp_t

    type(patch_t), dimension(4) :: patches
    real(wp), allocatable :: xc(:), yc(:), zc(:)
    !$omp declare target(patches, xc, yc, zc)

    public :: gp_t, s_init, s_finalize, s_dispatch

contains

    subroutine s_init(n)
        integer, intent(in) :: n
        integer :: i
        allocate(xc(n), yc(n), zc(n))
        do i = 1, n
            xc(i) = real(i, wp)*0.1_wp
            yc(i) = real(i, wp)*0.1_wp
            zc(i) = real(i, wp)*0.1_wp
        end do
        do i = 1, 4
            patches(i)%cx = 0.5_wp; patches(i)%cy = 0.5_wp; patches(i)%cz = 0.5_wp
            patches(i)%radius = 0.25_wp
            patches(i)%mat      = 0.0_wp
            patches(i)%mat(1,1) = 1.0_wp
            patches(i)%mat(2,2) = 1.0_wp
            patches(i)%mat(3,3) = 1.0_wp
        end do
        !$omp target enter data map(to: patches, xc, yc, zc)
    end subroutine s_init

    subroutine s_finalize()
        !$omp target exit data map(delete: patches, xc, yc, zc)
        deallocate(xc, yc, zc)
    end subroutine s_finalize

    ! Single !$omp target teams loop calling a single declare-target sub.
    subroutine s_dispatch(gps, n)
        type(gp_t), intent(inout) :: gps(:)
        integer,    intent(in)    :: n
        integer :: i

        !$omp target teams loop private(i) map(tofrom:gps)
        do i = 1, n
            call s_apply(gps(i))
        end do

    end subroutine s_dispatch

    ! !$omp declare target subroutine in the SAME module.
    ! ICE fires when the inliner pulls this into the target loop at -O1/-O2/-O3.
    subroutine s_apply(gp)
        !$omp declare target
        type(gp_t), intent(inout) :: gp
        integer :: id, i, j, k
        real(wp), dimension(3) :: dv, ldv
        real(wp) :: d
        id = gp%pid; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dv(1) = xc(i) - patches(id)%cx
        dv(2) = yc(j) - patches(id)%cy
        dv(3) = zc(k) - patches(id)%cz
        ! matmul with a derived-type member (patches(id)%mat) is the key trigger
        ldv = matmul(patches(id)%mat, dv)
        d = sqrt(sum(ldv**2))
        gp%levelset = d - patches(id)%radius
        if (d > 0.0_wp) then
            gp%norm = ldv/d
        else
            gp%norm = 0.0_wp
        end if
    end subroutine s_apply

end module m_ice

program test_ice
    use m_ice
    implicit none
    integer, parameter :: N = 16
    type(gp_t) :: gps(N)
    integer :: i
    call s_init(N)
    do i = 1, N
        gps(i)%pid      = mod(i-1,4)+1
        gps(i)%loc      = [mod(i,N)+1, mod(i*2,N)+1, mod(i*3,N)+1]
        gps(i)%levelset = 0.0d0
        gps(i)%norm     = 0.0d0
    end do
    call s_dispatch(gps, N)
    write(*,*) 'levelset(1) =', gps(1)%levelset
    call s_finalize()
end program test_ice
