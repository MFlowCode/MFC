! Minimal reproducer for ifx 2025.1.1 SPIR-V ICE #5633
!
! Pattern: single module with 10+ !$omp target teams loop regions, each calling
! a different !$omp declare target subroutine defined in the SAME module.
! The subroutines access a derived type that contains a character field
! (ib_patch_parameters%model_filepath) and 3x3 real arrays (rotation_matrix).
!
! Compile (ICE expected at -O3, pass expected at -O3 -fno-inline):
!   mpif90 -free -fiopenmp -fopenmp-targets=spir64 -fp-model=precise \
!          -march=native -mno-avx512fp16 -O3 -c m_cl4_mini.f90
!
! Fix:
!   mpif90 ... -O3 -fno-inline -c m_cl4_mini.f90

module m_cl4_mini

    implicit none
    private

    integer, parameter :: wp = kind(1.0d0)
    integer, parameter :: PATHLEN = 200
    integer, parameter :: MAX_PATCHES = 10

    ! Replicated from m_derived_types.fpp :: ib_patch_parameters
    type :: ib_patch_parameters
        integer  :: geometry
        real(wp) :: x_centroid, y_centroid, z_centroid
        real(wp) :: step_x_centroid, step_y_centroid, step_z_centroid
        real(wp), dimension(1:3)     :: centroid_offset
        real(wp), dimension(1:3)     :: angles, step_angles
        real(wp), dimension(1:3,1:3) :: rotation_matrix
        real(wp), dimension(1:3,1:3) :: rotation_matrix_inverse
        real(wp) :: c, p, t, m
        real(wp) :: length_x, length_y, length_z
        real(wp) :: radius, theta
        logical  :: slip
        character(LEN=PATHLEN) :: model_filepath  ! <-- char field; triggers SPIR-V bug
        real(wp), dimension(1:3) :: model_translate
        real(wp), dimension(1:3) :: model_scale
        real(wp), dimension(1:3) :: model_rotate
        integer  :: model_spc
        real(wp) :: model_threshold
        integer  :: moving_ibm
        real(wp) :: mass, moment
        real(wp), dimension(1:3) :: force, torque
        real(wp), dimension(1:3) :: vel, step_vel
        real(wp), dimension(1:3) :: angular_vel, step_angular_vel
    end type ib_patch_parameters

    ! Replicated from m_derived_types.fpp :: ghost_point
    type :: ghost_point
        integer, dimension(3)        :: loc
        real(wp), dimension(3)       :: ip_loc
        integer, dimension(3)        :: ip_grid
        real(wp), dimension(2, 2, 2) :: interp_coeffs  ! 3-D array in struct
        integer                      :: ib_patch_id
        real(wp)                     :: levelset
        real(wp), dimension(1:3)     :: levelset_norm
        logical                      :: slip
        integer, dimension(3)        :: DB
        integer                      :: x_periodicity, y_periodicity, z_periodicity
    end type ghost_point

    type :: bounds_info
        real(wp) :: beg, end
    end type bounds_info

    ! Module-level variables accessed by the declare-target subroutines
    type(ib_patch_parameters), dimension(MAX_PATCHES) :: patch_ib
    real(wp), allocatable :: x_cc(:), y_cc(:), z_cc(:)
    type(bounds_info) :: x_domain, y_domain, z_domain
    integer :: num_ibs

    !$omp declare target(patch_ib, x_cc, y_cc, z_cc, x_domain, y_domain, z_domain, num_ibs)

    public :: ghost_point, s_cl4_dispatch, s_cl4_init, s_cl4_finalize

contains

    ! ------------------------------------------------------------------ init --

    subroutine s_cl4_init(nx, nz)
        integer, intent(in) :: nx, nz
        integer :: i
        allocate(x_cc(nx), y_cc(nx), z_cc(nz))
        do i = 1, nx
            x_cc(i) = real(i, wp) * 0.1_wp
            y_cc(i) = real(i, wp) * 0.1_wp
        end do
        do i = 1, nz
            z_cc(i) = real(i, wp) * 0.1_wp
        end do
        x_domain%beg = 0._wp; x_domain%end = real(nx, wp)*0.1_wp
        y_domain%beg = 0._wp; y_domain%end = real(nx, wp)*0.1_wp
        z_domain%beg = 0._wp; z_domain%end = real(nz, wp)*0.1_wp
        num_ibs = 1
        do i = 1, MAX_PATCHES
            patch_ib(i)%geometry   = i
            patch_ib(i)%x_centroid = 0.5_wp
            patch_ib(i)%y_centroid = 0.5_wp
            patch_ib(i)%z_centroid = 0.5_wp
            patch_ib(i)%radius     = 0.2_wp
            patch_ib(i)%length_x   = 0.4_wp
            patch_ib(i)%length_y   = 0.4_wp
            patch_ib(i)%length_z   = 0.4_wp
            patch_ib(i)%rotation_matrix         = 0._wp
            patch_ib(i)%rotation_matrix(1,1)    = 1._wp
            patch_ib(i)%rotation_matrix(2,2)    = 1._wp
            patch_ib(i)%rotation_matrix(3,3)    = 1._wp
            patch_ib(i)%rotation_matrix_inverse = patch_ib(i)%rotation_matrix
            patch_ib(i)%centroid_offset         = 0._wp
            patch_ib(i)%model_filepath          = 'none'
        end do
        !$omp target enter data map(to: patch_ib, x_cc, y_cc, z_cc, x_domain, y_domain, z_domain, num_ibs)
    end subroutine s_cl4_init

    subroutine s_cl4_finalize()
        !$omp target exit data map(delete: patch_ib, x_cc, y_cc, z_cc, x_domain, y_domain, z_domain, num_ibs)
        deallocate(x_cc, y_cc, z_cc)
    end subroutine s_cl4_finalize

    ! ------------------------------------------------------- dispatch (ICE site) --

    !> 10 separate target teams loop regions, each calling a different
    !! !$omp declare target subroutine defined in THIS module.
    subroutine s_cl4_dispatch(gps, num_gps)
        type(ghost_point), dimension(:), intent(inout) :: gps
        integer, intent(in) :: num_gps
        integer :: i, patch_id

        ! 3D geometry loops (5 loops, geometries 8-12)
        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 8) call s_geo_sphere(gps(i))
        end do

        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 9) call s_geo_cuboid(gps(i))
        end do

        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 10) call s_geo_cylinder(gps(i))
        end do

        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 11) call s_geo_3d_airfoil(gps(i))
        end do

        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 12) call s_geo_model(gps(i))
        end do

        ! 2D geometry loops (5 loops, geometries 2-6)
        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 2) call s_geo_circle(gps(i))
        end do

        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 3) call s_geo_rectangle(gps(i))
        end do

        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 4) call s_geo_airfoil(gps(i))
        end do

        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 5) call s_geo_ellipse(gps(i))
        end do

        !$omp target teams loop private(i,patch_id) map(tofrom:gps) map(to:patch_ib,num_ibs)
        do i = 1, num_gps
            patch_id = gps(i)%ib_patch_id
            if (patch_ib(patch_id)%geometry == 6) call s_geo_triangle(gps(i))
        end do

    end subroutine s_cl4_dispatch

    ! ----------------------------------------- declare-target geometry routines --
    ! All defined in the SAME module; each has !$omp declare target.
    ! The ICE fires when the inliner pulls these into the target teams loop bodies.

    subroutine s_geo_sphere(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp) :: dist
        real(wp), dimension(3) :: dv
        integer :: id, i, j, k
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dv(1) = x_cc(i) - patch_ib(id)%x_centroid - real(gp%x_periodicity,wp)*(x_domain%end-x_domain%beg)
        dv(2) = y_cc(j) - patch_ib(id)%y_centroid - real(gp%y_periodicity,wp)*(y_domain%end-y_domain%beg)
        dv(3) = z_cc(k) - patch_ib(id)%z_centroid - real(gp%z_periodicity,wp)*(z_domain%end-z_domain%beg)
        dist = sqrt(sum(dv**2))
        gp%levelset = dist - patch_ib(id)%radius
        if (dist > 0._wp) then
            gp%levelset_norm = dv/dist
        else
            gp%levelset_norm = 0._wp
        end if
    end subroutine s_geo_sphere

    subroutine s_geo_cuboid(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp) :: dx, dy, dz, lx, ly, lz
        integer :: id, i, j, k
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        lx = patch_ib(id)%length_x; ly = patch_ib(id)%length_y; lz = patch_ib(id)%length_z
        dx = abs(x_cc(i) - patch_ib(id)%x_centroid) - lx/2._wp
        dy = abs(y_cc(j) - patch_ib(id)%y_centroid) - ly/2._wp
        dz = abs(z_cc(k) - patch_ib(id)%z_centroid) - lz/2._wp
        gp%levelset = max(dx, max(dy, dz))
        gp%levelset_norm(1) = sign(1._wp, x_cc(i) - patch_ib(id)%x_centroid)
        gp%levelset_norm(2) = sign(1._wp, y_cc(j) - patch_ib(id)%y_centroid)
        gp%levelset_norm(3) = sign(1._wp, z_cc(k) - patch_ib(id)%z_centroid)
    end subroutine s_geo_cuboid

    subroutine s_geo_cylinder(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp) :: r2d, dist
        real(wp), dimension(3) :: dv, rot_dv
        integer :: id, i, j, k
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dv(1) = x_cc(i) - patch_ib(id)%x_centroid
        dv(2) = y_cc(j) - patch_ib(id)%y_centroid
        dv(3) = z_cc(k) - patch_ib(id)%z_centroid
        rot_dv = matmul(patch_ib(id)%rotation_matrix_inverse, dv)
        r2d = sqrt(rot_dv(1)**2 + rot_dv(2)**2)
        dist = r2d - patch_ib(id)%radius
        gp%levelset = dist
        if (r2d > 0._wp) then
            gp%levelset_norm(1) = rot_dv(1)/r2d
            gp%levelset_norm(2) = rot_dv(2)/r2d
        else
            gp%levelset_norm(1) = 0._wp
            gp%levelset_norm(2) = 0._wp
        end if
        gp%levelset_norm(3) = 0._wp
    end subroutine s_geo_cylinder

    subroutine s_geo_3d_airfoil(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp), dimension(3) :: dv, local
        integer :: id, i, j, k
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dv(1) = x_cc(i) - patch_ib(id)%x_centroid
        dv(2) = y_cc(j) - patch_ib(id)%y_centroid
        dv(3) = z_cc(k) - patch_ib(id)%z_centroid
        local = matmul(patch_ib(id)%rotation_matrix_inverse, dv) - patch_ib(id)%centroid_offset
        gp%levelset = sqrt(sum(local**2)) - patch_ib(id)%radius
        gp%levelset_norm = local / max(sqrt(sum(local**2)), 1e-12_wp)
    end subroutine s_geo_3d_airfoil

    subroutine s_geo_model(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp), dimension(3) :: dv
        integer :: id, i, j, k
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2); k = gp%loc(3)
        dv(1) = x_cc(i) - patch_ib(id)%x_centroid + patch_ib(id)%model_translate(1)*patch_ib(id)%model_scale(1)
        dv(2) = y_cc(j) - patch_ib(id)%y_centroid + patch_ib(id)%model_translate(2)*patch_ib(id)%model_scale(2)
        dv(3) = z_cc(k) - patch_ib(id)%z_centroid + patch_ib(id)%model_translate(3)*patch_ib(id)%model_scale(3)
        gp%levelset = sqrt(sum(dv**2)) - patch_ib(id)%radius
        gp%levelset_norm = dv / max(sqrt(sum(dv**2)), 1e-12_wp)
    end subroutine s_geo_model

    subroutine s_geo_circle(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp) :: dist
        real(wp), dimension(3) :: dv
        integer :: id, i, j
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2)
        dv(1) = x_cc(i) - patch_ib(id)%x_centroid - real(gp%x_periodicity,wp)*(x_domain%end-x_domain%beg)
        dv(2) = y_cc(j) - patch_ib(id)%y_centroid - real(gp%y_periodicity,wp)*(y_domain%end-y_domain%beg)
        dv(3) = 0._wp
        dist = sqrt(sum(dv**2))
        gp%levelset = dist - patch_ib(id)%radius
        if (dist > 0._wp) then
            gp%levelset_norm = dv/dist
        else
            gp%levelset_norm = 0._wp
        end if
    end subroutine s_geo_circle

    subroutine s_geo_rectangle(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp) :: dx, dy
        integer :: id, i, j
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2)
        dx = abs(x_cc(i) - patch_ib(id)%x_centroid) - patch_ib(id)%length_x/2._wp
        dy = abs(y_cc(j) - patch_ib(id)%y_centroid) - patch_ib(id)%length_y/2._wp
        gp%levelset = max(dx, dy)
        gp%levelset_norm(1) = sign(1._wp, x_cc(i) - patch_ib(id)%x_centroid)
        gp%levelset_norm(2) = sign(1._wp, y_cc(j) - patch_ib(id)%y_centroid)
        gp%levelset_norm(3) = 0._wp
    end subroutine s_geo_rectangle

    subroutine s_geo_airfoil(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp), dimension(3) :: dv, local
        integer :: id, i, j
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2)
        dv(1) = x_cc(i) - patch_ib(id)%x_centroid
        dv(2) = y_cc(j) - patch_ib(id)%y_centroid
        dv(3) = 0._wp
        local = matmul(patch_ib(id)%rotation_matrix_inverse, dv) - patch_ib(id)%centroid_offset
        gp%levelset = sqrt(sum(local**2)) - patch_ib(id)%radius
        gp%levelset_norm = local / max(sqrt(sum(local**2)), 1e-12_wp)
    end subroutine s_geo_airfoil

    subroutine s_geo_ellipse(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp) :: ax, ay, dist
        real(wp), dimension(3) :: dv
        integer :: id, i, j
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2)
        ax = patch_ib(id)%length_x / 2._wp
        ay = patch_ib(id)%length_y / 2._wp
        dv(1) = (x_cc(i) - patch_ib(id)%x_centroid) / ax
        dv(2) = (y_cc(j) - patch_ib(id)%y_centroid) / ay
        dv(3) = 0._wp
        dist = sqrt(dv(1)**2 + dv(2)**2)
        gp%levelset = dist - 1._wp
        if (dist > 0._wp) then
            gp%levelset_norm = dv/dist
        else
            gp%levelset_norm = 0._wp
        end if
    end subroutine s_geo_ellipse

    subroutine s_geo_triangle(gp)
        !$omp declare target
        type(ghost_point), intent(inout) :: gp
        real(wp) :: d1, d2, d3
        real(wp), dimension(3) :: dv
        integer :: id, i, j
        id = gp%ib_patch_id; i = gp%loc(1); j = gp%loc(2)
        dv(1) = x_cc(i) - patch_ib(id)%x_centroid
        dv(2) = y_cc(j) - patch_ib(id)%y_centroid
        dv(3) = 0._wp
        d1 = dv(2) + patch_ib(id)%length_y/2._wp
        d2 = -dv(2) + dv(1)*sqrt(3._wp)/3._wp + patch_ib(id)%length_y/3._wp
        d3 = -dv(2) - dv(1)*sqrt(3._wp)/3._wp + patch_ib(id)%length_y/3._wp
        gp%levelset = min(max(-d1,-d2), max(-d1,-d3))
        gp%levelset_norm(1) = dv(1) / max(sqrt(sum(dv(1:2)**2)), 1e-12_wp)
        gp%levelset_norm(2) = dv(2) / max(sqrt(sum(dv(1:2)**2)), 1e-12_wp)
        gp%levelset_norm(3) = 0._wp
    end subroutine s_geo_triangle

end module m_cl4_mini

program test_cl4
    use m_cl4_mini
    implicit none

    integer, parameter :: wp = kind(1.0d0)
    integer, parameter :: N = 64, NGPS = 32
    type(ghost_point), allocatable :: gps(:)
    integer :: i

    call s_cl4_init(N, N)

    allocate(gps(NGPS))
    do i = 1, NGPS
        gps(i)%loc           = [mod(i,N)+1, mod(i*2,N)+1, mod(i*3,N)+1]
        gps(i)%ib_patch_id   = 1
        gps(i)%x_periodicity = 0
        gps(i)%y_periodicity = 0
        gps(i)%z_periodicity = 0
        gps(i)%levelset      = 0._wp
        gps(i)%levelset_norm = 0._wp
        gps(i)%interp_coeffs = 0._wp
        gps(i)%slip          = .false.
    end do

    call s_cl4_dispatch(gps, NGPS)

    write(*,'(a,f12.6)') 'levelset(1) = ', gps(1)%levelset
    deallocate(gps)
    call s_cl4_finalize()
    write(*,'(a)') 'OK'

end program test_cl4
