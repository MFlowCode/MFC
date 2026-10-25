!>
!! @file m_particle_cloud.fpp
!! @brief Generates particle clouds: converts particle_cloud specifications into particle centres, then into individual
!!        sphere/circle IB patches (written to the initial IB state file) or Euler-Lagrange bubbles (written to the bubble input
!!        file).

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Generates particle clouds in two stages. Each cloud geometry (box, a patch shape, hemisphere shell, file) first produces
!! only particle centres (particle_cloud_coords); particle_type then decides what they become: IB clouds become ib_patch_parameters
!! entries, EL clouds become rows of the Euler-Lagrange bubble input file. Under file_per_process it runs on every rank: each rank
!! computes the same deterministic placement (so no MPI broadcast of particle positions is needed) and keeps only the IBs
!! f_local_rank_owns_location says are its own. EL clouds are generated on rank 0 alone, since simulation reads one shared bubble
!! file. Otherwise only rank 0 runs it and keeps every particle.
module m_particle_cloud

    use m_global_parameters
    use m_constants
    use m_mpi_common
    use m_helper
    use m_model
    use m_patch_geometries

    implicit none

    private

    public :: s_generate_particle_clouds

    !> Particle centres of one cloud, independent of the particle type they become
    type particle_cloud_coords
        integer                               :: num_particles  !< Number of centres stored
        real(wp), allocatable, dimension(:,:) :: pos            !< Centres, (1:3, 1:capacity)
        real(wp), allocatable, dimension(:)   :: radius         !< Radius of each particle
    end type particle_cloud_coords

    !> The patch-shaped region of the cloud being generated, as an IB patch built from its particle_cloud entry, and the world-space
    !! box that holds it. Set by s_set_cloud_region.
    type(ib_patch_parameters) :: cloud_region
    real(wp), dimension(3)    :: region_bbox_min, region_bbox_max

contains

    !> Generate every particle cloud. IB clouds fill particle_cloud_ibs (only this rank's under file_per_process, see module docs);
    !! only the first num_particle_cloud_ibs entries are written, the remainder is left uninitialized. EL clouds fill the first
    !! num_lag_bubbles columns of lag_bubbles with bubble input rows on rank 0.
    impure subroutine s_generate_particle_clouds(glb_bounds, particle_cloud_ibs, num_particle_cloud_ibs, lag_bubbles, &
        & num_lag_bubbles)

        type(bounds_info), dimension(3), intent(in)                       :: glb_bounds
        type(ib_patch_parameters), allocatable, intent(out), dimension(:) :: particle_cloud_ibs
        integer, intent(out)                                              :: num_particle_cloud_ibs
        real(wp), allocatable, intent(out), dimension(:,:)                :: lag_bubbles
        integer, intent(out)                                              :: num_lag_bubbles
        type(particle_cloud_coords)                                       :: coords
        integer                                                           :: cloud_idx, glbl_idx
        real(wp)                                                          :: t_start, t_end

        allocate (particle_cloud_ibs(min(num_ib_patches_max_namelist, f_num_cloud_particles(particle_type_ib))))
        allocate (lag_bubbles(num_lag_bubble_input_vars, f_num_cloud_particles(particle_type_el)))
        num_particle_cloud_ibs = 0
        num_lag_bubbles = 0
        if (num_particle_clouds == 0) return

        call cpu_time(t_start)
        call s_initialize_cloud_shapes()

        glbl_idx = num_ibs
        do cloud_idx = 1, num_particle_clouds
            if (particle_cloud(cloud_idx)%particle_type == particle_type_el .and. proc_rank /= 0) cycle

            call s_get_cloud_coords(cloud_idx, coords)
            select case (particle_cloud(cloud_idx)%particle_type)
            case (particle_type_ib)
                call s_coords_to_cloud_ibs(cloud_idx, coords, glb_bounds, glbl_idx, particle_cloud_ibs, num_particle_cloud_ibs)
            case (particle_type_el)
                call s_coords_to_lag_bubbles(coords, lag_bubbles, num_lag_bubbles)
            case default
                call s_mpi_abort("Particle cloud particle_type is not a known particle type of MFC. Exiting.")
            end select
            deallocate (coords%pos, coords%radius)
        end do

        call cpu_time(t_end)
        if (proc_rank == 0) print '(a,i0,a,i0,a,f0.3,a)', 'Particle clouds placed ', glbl_idx - num_ibs, ' IBs and ', &
            & num_lag_bubbles, ' EL bubbles in ', t_end - t_start, ' seconds.'

    end subroutine s_generate_particle_clouds

    !> Total requested particle count over the clouds of one particle_type.
    pure function f_num_cloud_particles(particle_type) result(num_particles)

        integer, intent(in) :: particle_type
        integer             :: num_particles

        num_particles = sum(particle_cloud(1:num_particle_clouds)%num_particles, &
                            & mask=particle_cloud(1:num_particle_clouds)%particle_type == particle_type)

    end function f_num_cloud_particles

    !> Fills coords with every particle centre of cloud_idx, in the deterministic order global ids are assigned in. Rejection
    !! sampling handles every sampled geometry (the per-candidate geometry sampling lives in s_sample_cloud_candidate). Lattice
    !! packing fills a box or a patch shape (hemisphere-shell + lattice is rejected in case_validator.py), and file clouds ignore
    !! packing_method.
    impure subroutine s_get_cloud_coords(cloud_idx, coords)

        integer, intent(in)                      :: cloud_idx
        type(particle_cloud_coords), intent(out) :: coords
        integer                                  :: alloc_stat

        allocate (coords%pos(3, particle_cloud(cloud_idx)%num_particles), coords%radius(particle_cloud(cloud_idx)%num_particles), &
                  & stat=alloc_stat)
        if (alloc_stat /= 0) then
            call s_mpi_abort("Error :: Ran out of CPU memory trying to allocate particle cloud coordinates. " &
                             & // "Current system resources cannot place the specified number of particles.")
        end if
        coords%num_particles = 0

        if (particle_cloud(cloud_idx)%cloud_geometry == cloud_geometry_file) then
            call s_read_particle_cloud_file(cloud_idx, coords)
            return
        end if
        if (f_is_patch_cloud(cloud_idx)) call s_set_cloud_region(cloud_idx)

        select case (particle_cloud(cloud_idx)%packing_method)
        case (1)  ! rejection (random) packing method
            call s_particle_cloud_rejection_pack(cloud_idx, coords)
        case (2)  ! lattice packing method
            if (f_is_patch_cloud(cloud_idx)) then
                call s_particle_cloud_lattice_region(cloud_idx, coords)
            else
                call s_particle_cloud_lattice(cloud_idx, coords)
            end if
        case default
            call s_mpi_abort("Particle cloud packing method is not a known packing method of MFC. Exiting.")
        end select

    end subroutine s_get_cloud_coords

    !> Does cloud_idx fill a patch shape (an IB patch geometry number) rather than a box, hemisphere shell or file?
    pure function f_is_patch_cloud(cloud_idx) result(is_patch)

        integer, intent(in) :: cloud_idx
        logical             :: is_patch

        is_patch = all(particle_cloud(cloud_idx)%cloud_geometry /= [cloud_geometry_box, cloud_geometry_hemisphere_shell, &
                       & cloud_geometry_file])

    end function f_is_patch_cloud

    !> Loads what patch-shaped clouds test against: STL models, once and shared with ICPP patches, and airfoil surfaces, resolved to
    !! about 10 surface points per particle diameter.
    impure subroutine s_initialize_cloud_shapes()

        integer :: cloud_idx

        do cloud_idx = 1, num_particle_clouds
            select case (particle_cloud(cloud_idx)%cloud_geometry)
            case (cloud_geometry_stl_2d, cloud_geometry_stl_3d)
                call s_instantiate_STL_models()
            case (cloud_geometry_airfoil_2d, cloud_geometry_airfoil_3d)
                call s_initialize_airfoil_grid(particle_cloud(cloud_idx)%airfoil_id, 2._wp*particle_cloud(cloud_idx)%radius)
            end select
        end do

    end subroutine s_initialize_cloud_shapes

    !> Builds cloud_region, the IB patch with cloud_idx's shape (its geometry number, centroid, lengths, region_radius, angles,
    !! airfoil and model), and its world-space bounding box. The airfoil leading edge or STL model origin sits at the centroid.
    subroutine s_set_cloud_region(cloud_idx)

        integer, intent(in) :: cloud_idx

        cloud_region%geometry = particle_cloud(cloud_idx)%cloud_geometry
        cloud_region%x_centroid = particle_cloud(cloud_idx)%x_centroid
        cloud_region%y_centroid = particle_cloud(cloud_idx)%y_centroid
        cloud_region%z_centroid = particle_cloud(cloud_idx)%z_centroid
        cloud_region%length_x = particle_cloud(cloud_idx)%length_x
        cloud_region%length_y = particle_cloud(cloud_idx)%length_y
        cloud_region%length_z = particle_cloud(cloud_idx)%length_z
        cloud_region%radius = particle_cloud(cloud_idx)%region_radius
        cloud_region%angles = particle_cloud(cloud_idx)%angles
        cloud_region%airfoil_id = particle_cloud(cloud_idx)%airfoil_id
        cloud_region%model_id = particle_cloud(cloud_idx)%model_id
        cloud_region%centroid_offset = 0._wp
        call s_update_patch_rotation_matrix(cloud_region)

        call s_get_patch_bounding_box(cloud_region, [cloud_region%x_centroid, cloud_region%y_centroid, cloud_region%z_centroid], &
                                      & region_bbox_min, region_bbox_max)
        if (num_dims < 3) then
            region_bbox_min(3) = cloud_region%z_centroid
            region_bbox_max(3) = cloud_region%z_centroid
        end if

    end subroutine s_set_cloud_region

    !> Is pos a valid centre for the current patch-shaped cloud: inside its region, tested in the region's own frame exactly as an
    !! IB patch tests a cell, and inside the computational domain, which clips any part of the shape that lies outside it?
    function f_is_inside_cloud_region(pos) result(is_valid)

        real(wp), dimension(3), intent(in) :: pos
        logical                            :: is_valid
        real(wp), dimension(3)             :: xyz_local

        is_valid = pos(1) >= x_domain_glb%beg .and. pos(1) <= x_domain_glb%end .and. pos(2) >= y_domain_glb%beg .and. pos(2) &
                       & <= y_domain_glb%end
        if (num_dims == 3) is_valid = is_valid .and. pos(3) >= z_domain_glb%beg .and. pos(3) <= z_domain_glb%end
        if (.not. is_valid) return

        xyz_local = pos - [cloud_region%x_centroid, cloud_region%y_centroid, cloud_region%z_centroid]
        xyz_local = matmul(cloud_region%rotation_matrix_inverse, xyz_local)
        is_valid = f_is_inside_patch(cloud_region, xyz_local)

    end function f_is_inside_cloud_region

    !> Appends one particle, its centre and radius, to coords.
    subroutine s_add_cloud_coord(coords, centroid, radius)

        type(particle_cloud_coords), intent(inout) :: coords
        real(wp), dimension(3), intent(in)         :: centroid
        real(wp), intent(in)                       :: radius

        coords%num_particles = coords%num_particles + 1
        coords%pos(:,coords%num_particles) = centroid
        coords%radius(coords%num_particles) = radius

    end subroutine s_add_cloud_coord

    !> Rejection-samples particle centres into a box or hemisphere-shell region with a minimum centre-to-centre spacing. coords
    !! doubles as the placed-particle store the spatial hash indexes for overlap rejection. Only the per-candidate geometry sampling
    !! differs between box and hemisphere shell; it is delegated to s_sample_cloud_candidate. Candidates are uniform over the region
    !! (PDF = 1); a cloud with an analytic pdf additionally keeps each candidate with probability pdf/pdf_max, and a cloud without
    !! one draws no extra random numbers, so its placement is unchanged. Particle radii come from s_draw_cloud_radii, largest first;
    !! two particles must be at least the sum of their radii plus min_spacing apart.
    subroutine s_particle_cloud_rejection_pack(cloud_idx, coords)

        integer, intent(in)                        :: cloud_idx
        type(particle_cloud_coords), intent(inout) :: coords
        integer                                    :: n_target, seed
        integer(8)                                 :: n_attempts, max_attempts
        real(wp)                                   :: min_dist, rx, ry, rz, pdf, pdf_max, acceptance
        logical                                    :: overlaps, reject, periodic_pack, has_pdf
        integer                                    :: hash_size, slot
        integer                                    :: bx, by, bz, nx_bins, ny_bins, nz_bins
        integer, allocatable                       :: hash_head(:), chain_next(:)
        real(wp)                                   :: xmin, ymin, zmin, length_x, length_y, length_z, r_next
        real(wp), allocatable                      :: radii(:)

        n_target = particle_cloud(cloud_idx)%num_particles
        allocate (radii(n_target))
        call s_draw_cloud_radii(cloud_idx, radii)
        ! The largest pair separation sizes the spatial-hash bins
        min_dist = 2._wp*radii(1) + particle_cloud(cloud_idx)%min_spacing
        periodic_pack = particle_cloud(cloud_idx)%cloud_geometry == cloud_geometry_box .and. particle_cloud(cloud_idx)%periodic == 1
        length_x = particle_cloud(cloud_idx)%length_x
        length_y = particle_cloud(cloud_idx)%length_y
        length_z = particle_cloud(cloud_idx)%length_z
        xmin = particle_cloud(cloud_idx)%x_centroid - 0.5_wp*length_x
        ymin = particle_cloud(cloud_idx)%y_centroid - 0.5_wp*length_y
        zmin = particle_cloud(cloud_idx)%z_centroid - 0.5_wp*length_z
        ! Periodic wrapping bins the box; otherwise bins are absolute and the counts are unused (and lengths may be unset).
        nx_bins = 1
        ny_bins = 1
        nz_bins = 1
        if (periodic_pack) then
            nx_bins = max(1, ceiling(length_x/min_dist))
            ny_bins = max(1, ceiling(length_y/min_dist))
            if (num_dims == 3) nz_bins = max(1, ceiling(length_z/min_dist))
        end if

        n_attempts = 0
        seed = particle_cloud(cloud_idx)%seed
        if (seed == 0) seed = 1 + cloud_idx*1013904223

        call s_cloud_pdf(cloud_idx, 0._wp, 0._wp, 0._wp, has_pdf, pdf)
        acceptance = 1._wp
        if (has_pdf) call s_estimate_cloud_pdf(cloud_idx, seed, pdf_max, acceptance)
        max_attempts = int(min(real(n_target, wp)*1000._wp/acceptance, real(huge(max_attempts), wp)), 8)

        ! Hash table: 4x overprovisioned for ~25% load factor, minimum 16 buckets. chain_next(i) links placed particle i to the
        ! previous occupant of its bucket.
        hash_size = max(16, 4*n_target)
        allocate (hash_head(hash_size))
        allocate (chain_next(n_target))
        hash_head = -1
        chain_next = -1

        do while (coords%num_particles < n_target .and. n_attempts < max_attempts)
            n_attempts = n_attempts + 1

            r_next = radii(coords%num_particles + 1)
            call s_sample_cloud_candidate(cloud_idx, seed, r_next, rx, ry, rz, reject)
            if (reject) cycle

            if (has_pdf) then
                call s_cloud_pdf(cloud_idx, rx, ry, rz, has_pdf, pdf)
                @:PROHIBIT(.not. (pdf >= 0._wp), "particle_cloud pdf is negative or NaN inside its region")
                if (pdf > pdf_max) then
                    if (proc_rank == 0) print '(a,i0,a)', 'WARNING: particle_cloud(', cloud_idx, &
                        & ')%pdf exceeds its sampled maximum; raising the bound, which slightly biases earlier placements.'
                    pdf_max = 1.1_wp*pdf
                end if
                if (f_xorshift(seed)*pdf_max > pdf) cycle
            end if

            call s_check_cloud_particle_overlap(rx, ry, rz, r_next, coords%pos, coords%radius, hash_head, chain_next, hash_size, &
                                                & min_dist, particle_cloud(cloud_idx)%min_spacing, periodic_pack, xmin, ymin, &
                                                & zmin, length_x, length_y, length_z, nx_bins, ny_bins, nz_bins, overlaps, bx, &
                                                & by, bz)

            if (.not. overlaps) then
                call s_add_cloud_coord(coords, [rx, ry, rz], r_next)

                ! Insert into hash grid as head of bucket chain
                slot = f_bin_hash(bx, by, bz, hash_size)
                chain_next(coords%num_particles) = hash_head(slot)
                hash_head(slot) = coords%num_particles
            end if
        end do

        if (coords%num_particles < n_target) then
            call s_mpi_abort("Error :: Failed to place all particles in particle bed")
        end if

        deallocate (hash_head, chain_next, radii)

    end subroutine s_particle_cloud_rejection_pack

    !> Draws the radius of every particle in cloud_idx, sorted largest first so the biggest find room while the region is empty. One
    !! size by default; size classes take their number fractions exactly; a lognormal in ln(r/radius) with standard deviation
    !! radius_sigma is drawn from its own stream (Box-Muller on f_xorshift) and redrawn outside [radius_min, radius_max].
    subroutine s_draw_cloud_radii(cloud_idx, radii)

        integer, intent(in)                        :: cloud_idx
        real(wp), dimension(:), intent(out)        :: radii
        integer                                    :: n, k, i, seed, n_done, n_class
        integer, dimension(num_radius_classes_max) :: order
        real(wp)                                   :: r_min, r_max, r, cumulative, total

        n = size(radii)
        radii = particle_cloud(cloud_idx)%radius
        call s_cloud_radius_bounds(particle_cloud(cloud_idx), r_min, r_max)

        if (.not. f_is_default(particle_cloud(cloud_idx)%radius_sigma)) then
            seed = ieor(particle_cloud(cloud_idx)%seed, 2135587861)
            if (seed == 0) seed = 1
            do i = 1, n
                do
                    r = particle_cloud(cloud_idx)%radius*exp(particle_cloud(cloud_idx)%radius_sigma*sqrt(-2._wp*log(max(f_xorshift(seed), &
                                       & tiny(1._wp))))*cos(2._wp*pi*f_xorshift(seed)))
                    if (r >= r_min .and. r <= r_max) exit
                end do
                radii(i) = r
            end do
            call s_sort_descending(radii)
        else if (any(particle_cloud(cloud_idx)%radius_fraction > 0._wp)) then
            ! Classes largest first; cumulative rounding makes the counts sum to n exactly
            total = sum(particle_cloud(cloud_idx)%radius_fraction)
            order = [(k, k=1, num_radius_classes_max)]
            do k = 1, num_radius_classes_max - 1
                do i = k + 1, num_radius_classes_max
                    if (particle_cloud(cloud_idx)%radius_ratio(order(i)) > particle_cloud(cloud_idx)%radius_ratio(order(k))) then
                        order([k, i]) = order([i, k])
                    end if
                end do
            end do
            n_done = 0
            cumulative = 0._wp
            do k = 1, num_radius_classes_max
                cumulative = cumulative + particle_cloud(cloud_idx)%radius_fraction(order(k))/total
                n_class = nint(cumulative*real(n, wp)) - n_done
                radii(n_done + 1:n_done + n_class) &
                      & = particle_cloud(cloud_idx)%radius*particle_cloud(cloud_idx)%radius_ratio(order(k))
                n_done = n_done + n_class
            end do
        end if

    end subroutine s_draw_cloud_radii

    !> Sorts a in place, largest first (heapsort).
    subroutine s_sort_descending(a)

        real(wp), dimension(:), intent(inout) :: a
        integer                               :: n, i, root, child
        real(wp)                              :: tmp

        n = size(a)
        do i = n/2, 1, -1
            call s_sift(i, n)
        end do
        do i = n, 2, -1
            tmp = a(1)
            a(1) = a(i)
            a(i) = tmp
            call s_sift(1, i - 1)
        end do
        a = a(n:1:-1)

    contains

        !> Restores the max-heap below start within a(1:last).
        subroutine s_sift(start, last)

            integer, intent(in) :: start, last

            root = start
            do while (2*root <= last)
                child = 2*root
                if (child < last) then
                    if (a(child + 1) > a(child)) child = child + 1
                end if
                if (a(root) >= a(child)) return
                tmp = a(root)
                a(root) = a(child)
                a(child) = tmp
                root = child
            end do

        end subroutine s_sift

    end subroutine s_sort_descending

    !> Evaluates cloud_idx's analytic placement PDF, particle_cloud(i)%pdf, at (x, y, z). The case-file expression is compiled in
    !! through case.fpp (toolchain/mfc/case.py); has_pdf is false and pdf is 1 for a cloud without one.
    subroutine s_cloud_pdf(cloud_idx, x, y, z, has_pdf, pdf)

        integer, intent(in)   :: cloud_idx
        real(wp), intent(in)  :: x, y, z
        logical, intent(out)  :: has_pdf
        real(wp), intent(out) :: pdf

        has_pdf = .false.
        pdf = 1._wp
        #:if defined('particle_cloud_pdf')
            @:particle_cloud_pdf()
        #:endif

    end subroutine s_cloud_pdf

    !> Estimates the maximum of cloud_idx's PDF over its region, padded by 10%, and the fraction of uniform candidates rejection
    !! will keep, from num_pdf_samples uniform candidates. They come from their own stream, derived from seed, so the placement
    !! stream itself is untouched.
    subroutine s_estimate_cloud_pdf(cloud_idx, seed, pdf_max, acceptance)

        integer, intent(in)   :: cloud_idx, seed
        real(wp), intent(out) :: pdf_max, acceptance
        integer, parameter    :: num_pdf_samples = 100000
        integer               :: n, pdf_seed
        real(wp)              :: rx, ry, rz, pdf, pdf_sum
        logical               :: reject, has_pdf

        pdf_seed = ieor(seed, 1540483477)
        if (pdf_seed == 0) pdf_seed = 1
        n = 0
        pdf_sum = 0._wp
        pdf_max = 0._wp
        do while (n < num_pdf_samples)
            call s_sample_cloud_candidate(cloud_idx, pdf_seed, particle_cloud(cloud_idx)%radius, rx, ry, rz, reject)
            if (reject) cycle
            call s_cloud_pdf(cloud_idx, rx, ry, rz, has_pdf, pdf)
            @:PROHIBIT(.not. (pdf >= 0._wp), "particle_cloud pdf is negative or NaN inside its region")
            n = n + 1
            pdf_sum = pdf_sum + pdf
            pdf_max = max(pdf_max, pdf)
        end do

        @:PROHIBIT(pdf_max <= 0._wp, "particle_cloud pdf is zero everywhere inside its region")
        pdf_max = 1.1_wp*pdf_max
        acceptance = pdf_sum/(real(n, wp)*pdf_max)

    end subroutine s_estimate_cloud_pdf

    !> Draws one rejection-sampling candidate centre (rx, ry, rz) for cloud_idx, advancing seed in place. For box geometry the
    !! candidate is uniform in the box and never rejected. For a hemisphere shell the candidate is uniform in the shell volume - 2D
    !! uses theta uniform on [0, pi] with the sqrt radial CDF; 3D uses uniform phi, uniform cos(polar) on [0, 1], and the cube-root
    !! radial CDF - and reject is set when it lands within one particle radius of the flat face (the plane through the centroid
    !! perpendicular to shell_axis: 1=x, 2=y, 3=z; 2D has no z-axis so any value other than 1 falls back to y, matching the
    !! pre-shell_axis default), a hard geometric cut applied after sampling that preserves uniformity over the remaining region. For
    !! a patch shape the candidate is uniform in the region's bounding box and rejected outside the shape or the domain.
    subroutine s_sample_cloud_candidate(cloud_idx, seed, particle_radius, rx, ry, rz, reject)

        integer, intent(in)    :: cloud_idx
        integer, intent(inout) :: seed
        real(wp), intent(in)   :: particle_radius
        real(wp), intent(out)  :: rx, ry, rz
        logical, intent(out)   :: reject
        real(wp)               :: xmin, xmax, ymin, ymax, zmin, zmax
        real(wp)               :: theta, phi, r_shell, rho, u, zdir, r_inner, r_outer

        reject = .false.

        select case (particle_cloud(cloud_idx)%cloud_geometry)
        case (cloud_geometry_box)
            xmin = particle_cloud(cloud_idx)%x_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_x
            xmax = particle_cloud(cloud_idx)%x_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_x
            ymin = particle_cloud(cloud_idx)%y_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_y
            ymax = particle_cloud(cloud_idx)%y_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_y
            zmin = particle_cloud(cloud_idx)%z_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_z
            zmax = particle_cloud(cloud_idx)%z_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_z

            rx = xmin + f_xorshift(seed)*(xmax - xmin)
            ry = ymin + f_xorshift(seed)*(ymax - ymin)
            if (num_dims < 3) then
                rz = particle_cloud(cloud_idx)%z_centroid
            else
                rz = zmin + f_xorshift(seed)*(zmax - zmin)
            end if
        case (cloud_geometry_hemisphere_shell)
            r_inner = particle_cloud(cloud_idx)%shell_inner_radius + particle_radius
            r_outer = particle_cloud(cloud_idx)%shell_outer_radius - particle_radius

            if (num_dims < 3) then
                theta = pi*f_xorshift(seed)
                u = f_xorshift(seed)
                r_shell = sqrt((r_outer**2 - r_inner**2)*u + r_inner**2)
                rz = particle_cloud(cloud_idx)%z_centroid
                if (particle_cloud(cloud_idx)%shell_axis == 1) then
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*sin(theta)
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*cos(theta)
                    if (rx < particle_cloud(cloud_idx)%x_centroid + particle_radius) reject = .true.
                else
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*cos(theta)
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*sin(theta)
                    if (ry < particle_cloud(cloud_idx)%y_centroid + particle_radius) reject = .true.
                end if
            else
                phi = 2._wp*pi*f_xorshift(seed)
                zdir = f_xorshift(seed)
                rho = sqrt(max(0._wp, 1._wp - zdir**2))
                u = f_xorshift(seed)
                r_shell = ((r_outer**3 - r_inner**3)*u + r_inner**3)**(1._wp/3._wp)
                select case (particle_cloud(cloud_idx)%shell_axis)
                case (1)  ! opens toward +x
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*zdir
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*rho*cos(phi)
                    rz = particle_cloud(cloud_idx)%z_centroid + r_shell*rho*sin(phi)
                    if (rx < particle_cloud(cloud_idx)%x_centroid + particle_radius) reject = .true.
                case (2)  ! opens toward +y
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*zdir
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*rho*cos(phi)
                    rz = particle_cloud(cloud_idx)%z_centroid + r_shell*rho*sin(phi)
                    if (ry < particle_cloud(cloud_idx)%y_centroid + particle_radius) reject = .true.
                case default  ! 3: opens toward +z
                    rz = particle_cloud(cloud_idx)%z_centroid + r_shell*zdir
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*rho*cos(phi)
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*rho*sin(phi)
                    if (rz < particle_cloud(cloud_idx)%z_centroid + particle_radius) reject = .true.
                end select
            end if
        case default  ! patch shape: uniform in its bounding box, kept inside the shape
            rx = region_bbox_min(1) + f_xorshift(seed)*(region_bbox_max(1) - region_bbox_min(1))
            ry = region_bbox_min(2) + f_xorshift(seed)*(region_bbox_max(2) - region_bbox_min(2))
            rz = region_bbox_min(3)
            if (num_dims == 3) rz = region_bbox_min(3) + f_xorshift(seed)*(region_bbox_max(3) - region_bbox_min(3))
            reject = .not. f_is_inside_cloud_region([rx, ry, rz])
        end select

    end subroutine s_sample_cloud_candidate

    !> Places particles on the optimally dense lattice for the cloud region: a triangular lattice in 2D, a face-centered cubic
    !! lattice in 3D. The lattice spacing is set by the particle density (num_particles over the region area/volume); if that
    !! spacing falls below the required centre-to-centre distance (2*radius + min_spacing), the region is too dense and the run is
    !! aborted.
    subroutine s_particle_cloud_lattice(cloud_idx, coords)

        integer, intent(in)                        :: cloud_idx
        type(particle_cloud_coords), intent(inout) :: coords
        integer                                    :: n_target, row, col, ncx, ncy, ix, jy, kz, b
        real(wp)                                   :: xmin, xmax, ymin, ymax, zmin, zmax, min_dist
        real(wp)                                   :: spacing, row_dy, cell, x0, px, py
        real(wp), dimension(4)                     :: bx_off, by_off, bz_off

        xmin = particle_cloud(cloud_idx)%x_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_x
        xmax = particle_cloud(cloud_idx)%x_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_x
        ymin = particle_cloud(cloud_idx)%y_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_y
        ymax = particle_cloud(cloud_idx)%y_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_y
        zmin = particle_cloud(cloud_idx)%z_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_z
        zmax = particle_cloud(cloud_idx)%z_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_z

        min_dist = 2._wp*particle_cloud(cloud_idx)%radius + particle_cloud(cloud_idx)%min_spacing
        n_target = particle_cloud(cloud_idx)%num_particles

        if (num_dims < 3) then
            ! Triangular lattice: area per particle = (sqrt(3)/2)*spacing**2.
            spacing = sqrt(2._wp*(xmax - xmin)*(ymax - ymin)/(sqrt(3._wp)*real(n_target, wp)))
        else
            ! Face-centered cubic lattice: volume per particle = spacing**3/sqrt(2).
            spacing = (sqrt(2._wp)*(xmax - xmin)*(ymax - ymin)*(zmax - zmin)/real(n_target, wp))**(1._wp/3._wp)
        end if

        if (spacing < min_dist) then
            call s_mpi_abort("Error :: Particle cloud is too dense for lattice packing; " &
                             & // "reduce num_particles or min_spacing, or enlarge the cloud region")
        end if

        if (num_dims < 3) then
            ! Triangular lattice: rows pitched by spacing*sqrt(3)/2, odd rows shifted by half a spacing.
            row_dy = spacing*sqrt(3._wp)/2._wp
            row = 0
            do while (coords%num_particles < n_target)
                py = ymin + real(row, wp)*row_dy
                x0 = xmin
                if (mod(row, 2) == 1) x0 = xmin + 0.5_wp*spacing
                col = 0
                px = x0
                do while (px <= xmax .and. coords%num_particles < n_target)
                    call s_add_cloud_coord(coords, [px, py, particle_cloud(cloud_idx)%z_centroid], particle_cloud(cloud_idx)%radius)
                    col = col + 1
                    px = x0 + real(col, wp)*spacing
                end do
                row = row + 1
            end do
        else
            ! Face-centered cubic lattice via the conventional cubic cell (side = spacing*sqrt(2)) and its four basis points.
            cell = spacing*sqrt(2._wp)
            bx_off = [0._wp, 0.5_wp, 0.5_wp, 0._wp]*cell
            by_off = [0._wp, 0.5_wp, 0._wp, 0.5_wp]*cell
            bz_off = [0._wp, 0._wp, 0.5_wp, 0.5_wp]*cell
            ncx = max(1, ceiling((xmax - xmin)/cell))
            ncy = max(1, ceiling((ymax - ymin)/cell))
            kz = 0
            do while (coords%num_particles < n_target)
                do jy = 0, ncy - 1
                    do ix = 0, ncx - 1
                        do b = 1, 4
                            if (coords%num_particles >= n_target) exit
                            call s_add_cloud_coord(coords, [xmin + real(ix, wp)*cell + bx_off(b), ymin + real(jy, &
                                                   & wp)*cell + by_off(b), zmin + real(kz, wp)*cell + bz_off(b)], &
                                                   & particle_cloud(cloud_idx)%radius)
                        end do
                    end do
                end do
                kz = kz + 1
            end do
        end if

    end subroutine s_particle_cloud_lattice

    !> Places num_particles centres on a triangular (2D) or face-centered cubic (3D) lattice laid over the bounding box of the
    !! patch-shaped region, keeping the sites f_is_inside_cloud_region accepts and, with an analytic pdf, thinning them by it. The
    !! spacing is the largest at which that still leaves num_particles sites: halved from the box size until enough sites fit, then
    !! bisected. Surplus sites at that spacing are dropped evenly through the lattice order.
    subroutine s_particle_cloud_lattice_region(cloud_idx, coords)

        integer, intent(in)                        :: cloud_idx
        type(particle_cloud_coords), intent(inout) :: coords
        integer                                    :: n_target, n_lo, n_mid, iter
        real(wp)                                   :: min_dist, lo, hi, mid, pdf_max, acceptance, pdf
        logical                                    :: has_pdf

        n_target = particle_cloud(cloud_idx)%num_particles
        min_dist = 2._wp*particle_cloud(cloud_idx)%radius + particle_cloud(cloud_idx)%min_spacing
        call s_cloud_pdf(cloud_idx, 0._wp, 0._wp, 0._wp, has_pdf, pdf)
        pdf_max = 1._wp
        if (has_pdf) call s_estimate_cloud_pdf(cloud_idx, particle_cloud(cloud_idx)%seed, pdf_max, acceptance)

        hi = maxval(region_bbox_max(1:num_dims) - region_bbox_min(1:num_dims))
        lo = hi
        n_lo = f_region_lattice_sites(cloud_idx, lo, has_pdf, pdf_max)
        do while (n_lo < n_target)
            @:PROHIBIT(lo <= min_dist, &
                       & "Particle cloud is too dense for lattice packing in its region; " &
                       & // "reduce num_particles or min_spacing, or enlarge the region")
            hi = lo
            lo = max(0.5_wp*lo, min_dist)
            n_lo = f_region_lattice_sites(cloud_idx, lo, has_pdf, pdf_max)
        end do

        do iter = 1, 60
            if (hi - lo <= 1.e-10_wp*hi) exit
            mid = 0.5_wp*(lo + hi)
            n_mid = f_region_lattice_sites(cloud_idx, mid, has_pdf, pdf_max)
            if (n_mid >= n_target) then
                lo = mid
                n_lo = n_mid
            else
                hi = mid
            end if
        end do

        n_mid = f_region_lattice_sites(cloud_idx, lo, has_pdf, pdf_max, n_lo, coords)

    end subroutine s_particle_cloud_lattice_region

    !> Counts the lattice sites at the given spacing that the current patch-shaped region keeps (see
    !! s_particle_cloud_lattice_region). Given n_sites, that count, and coords, it also appends num_particles of them to coords,
    !! keeping site k when floor(k*num_particles/n_sites) advances, which spreads the dropped surplus evenly. With a pdf, sites are
    !! thinned from a stream reset on every call, so the count and the fill see the same sites.
    function f_region_lattice_sites(cloud_idx, spacing, has_pdf, pdf_max, n_sites, coords) result(n_valid)

        integer, intent(in)                                  :: cloud_idx
        real(wp), intent(in)                                 :: spacing, pdf_max
        logical, intent(in)                                  :: has_pdf
        integer, intent(in), optional                        :: n_sites
        type(particle_cloud_coords), intent(inout), optional :: coords
        integer                                              :: n_valid
        integer                                              :: seed, nx, ny, nz, ix, jy, kz, b, n_basis
        real(wp)                                             :: row_dy, cell, pdf
        real(wp), dimension(3)                               :: site
        real(wp), dimension(3, 4)                            :: basis
        logical                                              :: has_pdf_at

        seed = ieor(particle_cloud(cloud_idx)%seed, 1013904223)
        if (seed == 0) seed = 1
        n_valid = 0

        if (num_dims < 3) then
            ! Triangular lattice: rows pitched by spacing*sqrt(3)/2, odd rows shifted by half a spacing.
            row_dy = spacing*sqrt(3._wp)/2._wp
            nx = floor((region_bbox_max(1) - region_bbox_min(1))/spacing) + 1
            ny = floor((region_bbox_max(2) - region_bbox_min(2))/row_dy) + 1
            nz = 1
            n_basis = 1
            basis(:,1) = 0._wp
        else
            ! Face-centered cubic lattice via the conventional cubic cell (side = spacing*sqrt(2)) and its four basis points.
            cell = spacing*sqrt(2._wp)
            nx = floor((region_bbox_max(1) - region_bbox_min(1))/cell) + 1
            ny = floor((region_bbox_max(2) - region_bbox_min(2))/cell) + 1
            nz = floor((region_bbox_max(3) - region_bbox_min(3))/cell) + 1
            n_basis = 4
            basis = reshape([0._wp, 0._wp, 0._wp, 0.5_wp, 0.5_wp, 0._wp, 0.5_wp, 0._wp, 0.5_wp, 0._wp, 0.5_wp, 0.5_wp]*cell, [3, 4])
        end if

        do kz = 0, nz - 1
            do jy = 0, ny - 1
                do ix = 0, nx - 1
                    do b = 1, n_basis
                        if (num_dims < 3) then
                            site = region_bbox_min + [real(ix, wp)*spacing + merge(0.5_wp*spacing, 0._wp, mod(jy, 2) == 1), &
                                                           & real(jy, wp)*row_dy, 0._wp]
                        else
                            site = region_bbox_min + real([ix, jy, kz], wp)*cell + basis(:,b)
                        end if
                        if (any(site(1:num_dims) > region_bbox_max(1:num_dims))) cycle
                        if (.not. f_is_inside_cloud_region(site)) cycle
                        if (has_pdf) then
                            call s_cloud_pdf(cloud_idx, site(1), site(2), site(3), has_pdf_at, pdf)
                            if (f_xorshift(seed)*pdf_max > pdf) cycle
                        end if

                        n_valid = n_valid + 1
                        if (present(coords)) then
                            if ((int(n_valid, 8)*particle_cloud(cloud_idx)%num_particles)/n_sites > (int(n_valid - 1, &
                                & 8)*particle_cloud(cloud_idx)%num_particles)/n_sites) then
                                call s_add_cloud_coord(coords, site, particle_cloud(cloud_idx)%radius)
                            end if
                        end if
                    end do
                end do
            end do
        end do

    end function f_region_lattice_sites

    !> Reads cloud_idx's particle centres from particle_cloud(cloud_idx)%filepath, one whitespace-separated "x y z" row per
    !! particle; blank and '#' lines are skipped. 2D reads only x and y and places the particle at z_centroid. The row count must
    !! equal num_particles, which simulation uses to number the cloud's IBs.
    impure subroutine s_read_particle_cloud_file(cloud_idx, coords)

        integer, intent(in)                        :: cloud_idx
        type(particle_cloud_coords), intent(inout) :: coords
        character(LEN=pathlen_max)                 :: line
        real(wp), dimension(3)                     :: centroid
        integer                                    :: file_unit, ios

        open (newunit=file_unit, file=trim(particle_cloud(cloud_idx)%filepath), form='formatted', status='old', action='read', &
              & iostat=ios)
        if (ios /= 0) call s_mpi_abort("Cannot open particle cloud file: " // trim(particle_cloud(cloud_idx)%filepath))

        centroid(3) = particle_cloud(cloud_idx)%z_centroid
        do
            read (file_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle

            read (line, *, iostat=ios) centroid(1:min(num_dims, 3))
            if (ios /= 0) call s_mpi_abort("Malformed row in particle cloud file: " // trim(particle_cloud(cloud_idx)%filepath))
            if (coords%num_particles == particle_cloud(cloud_idx)%num_particles) then
                call s_mpi_abort("Particle cloud file holds more rows than num_particles: " &
                                 & // trim(particle_cloud(cloud_idx)%filepath))
            end if
            call s_add_cloud_coord(coords, centroid, particle_cloud(cloud_idx)%radius)
        end do
        close (file_unit)

        if (coords%num_particles /= particle_cloud(cloud_idx)%num_particles) then
            call s_mpi_abort("Particle cloud file holds fewer rows than num_particles: " &
                             & // trim(particle_cloud(cloud_idx)%filepath))
        end if

    end subroutine s_read_particle_cloud_file

    !> Appends one IB patch per centre in coords to particle_cloud_ibs. Every centre takes the next global id, so ids agree across
    !! ranks, but under file_per_process only the centres this rank owns are kept.
    subroutine s_coords_to_cloud_ibs(cloud_idx, coords, glb_bounds, glbl_idx, particle_cloud_ibs, num_particle_cloud_ibs)

        integer, intent(in)                                    :: cloud_idx
        type(particle_cloud_coords), intent(in)                :: coords
        type(bounds_info), dimension(3), intent(in)            :: glb_bounds
        integer, intent(inout)                                 :: glbl_idx, num_particle_cloud_ibs
        type(ib_patch_parameters), intent(inout), dimension(:) :: particle_cloud_ibs
        integer                                                :: i

        do i = 1, coords%num_particles
            glbl_idx = glbl_idx + 1
            if (file_per_process) then
                if (.not. f_local_rank_owns_location(coords%pos(:,i), glb_bounds)) cycle
            end if
            call s_add_cloud_particle(num_particle_cloud_ibs, glbl_idx, coords%pos(:,i), coords%radius(i), particle_cloud_ibs)
        end do

    end subroutine s_coords_to_cloud_ibs

    !> Writes a single placed particle into particle_cloud_ibs at the next free slot, advancing ib_idx, tagged with its
    !! already-assigned, absolute global patch id via glbl_idx. Only the fields s_write_ib_state_0_file writes are set; simulation
    !! fills every other property in s_assign_particle_cloud_ib_defaults (src/simulation/m_start_up.fpp).
    subroutine s_add_cloud_particle(ib_idx, glbl_idx, centroid, radius, particle_cloud_ibs)

        integer, intent(in)                                    :: glbl_idx
        integer, intent(inout)                                 :: ib_idx
        real(wp), dimension(3), intent(in)                     :: centroid
        real(wp), intent(in)                                   :: radius
        type(ib_patch_parameters), intent(inout), dimension(:) :: particle_cloud_ibs

        ib_idx = ib_idx + 1
        @:PROHIBIT(ib_idx > size(particle_cloud_ibs), &
                   & "Too many particle-cloud IBs on one rank. Modify case file or increase num_ib_patches_max_namelist.")

        particle_cloud_ibs(ib_idx)%gbl_patch_id = glbl_idx
        particle_cloud_ibs(ib_idx)%x_centroid = centroid(1)
        particle_cloud_ibs(ib_idx)%y_centroid = centroid(2)
        particle_cloud_ibs(ib_idx)%z_centroid = centroid(3)
        particle_cloud_ibs(ib_idx)%vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%angular_vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%angles(:) = 0._wp
        particle_cloud_ibs(ib_idx)%radius = radius

    end subroutine s_add_cloud_particle

    !> Appends one Euler-Lagrange bubble input row per particle in coords to lag_bubbles: at rest, with its own radius as its
    !! equilibrium and initial radius.
    subroutine s_coords_to_lag_bubbles(coords, lag_bubbles, num_lag_bubbles)

        type(particle_cloud_coords), intent(in) :: coords
        real(wp), intent(inout), dimension(:,:) :: lag_bubbles
        integer, intent(inout)                  :: num_lag_bubbles
        integer                                 :: i

        do i = 1, coords%num_particles
            num_lag_bubbles = num_lag_bubbles + 1
            lag_bubbles(:,num_lag_bubbles) = 0._wp
            lag_bubbles(1:3,num_lag_bubbles) = coords%pos(:,i)
            lag_bubbles(7, num_lag_bubbles) = coords%radius(i)
        end do

    end subroutine s_coords_to_lag_bubbles

    !> Convert a candidate particle centre to spatial-hash bin coordinates.
    subroutine s_get_cloud_bin(px, py, pz, min_dist, periodic_pack, xmin, ymin, zmin, nx_bins, ny_bins, nz_bins, bx, by, bz)

        real(wp), intent(in) :: px, py, pz, min_dist
        logical, intent(in)  :: periodic_pack
        real(wp), intent(in) :: xmin, ymin, zmin
        integer, intent(in)  :: nx_bins, ny_bins, nz_bins
        integer, intent(out) :: bx, by, bz

        if (periodic_pack) then
            bx = modulo(int(floor((px - xmin)/min_dist)), nx_bins)
            by = modulo(int(floor((py - ymin)/min_dist)), ny_bins)
            if (num_dims < 3) then
                bz = 0
            else
                bz = modulo(int(floor((pz - zmin)/min_dist)), nz_bins)
            end if
        else
            bx = int(floor(px/min_dist))
            by = int(floor(py/min_dist))
            if (num_dims < 3) then
                bz = 0
            else
                bz = int(floor(pz/min_dist))
            end if
        end if

    end subroutine s_get_cloud_bin

    !> Check whether a candidate particle of radius p_radius overlaps any already-placed particle (closer than the sum of their
    !! radii plus gap) in neighbouring spatial-hash bins, which are min_dist, the largest pair separation, wide. Scans the 3x3(x3)
    !! bin neighborhood - O(1) average via hash lookup - and also hands back the candidate's own bin so the caller can insert it
    !! without recomputing.
    subroutine s_check_cloud_particle_overlap(px, py, pz, p_radius, placed, placed_radius, hash_head, chain_next, hash_size, &
        & min_dist, gap, periodic_pack, xmin, ymin, zmin, length_x, length_y, length_z, nx_bins, ny_bins, nz_bins, overlaps, bx, &
        & by, bz)

        real(wp), intent(in)                 :: px, py, pz, p_radius, min_dist, gap
        real(wp), intent(in), dimension(:,:) :: placed
        real(wp), intent(in), dimension(:)   :: placed_radius
        integer, intent(in), dimension(:)    :: hash_head, chain_next
        integer, intent(in)                  :: hash_size
        logical, intent(in)                  :: periodic_pack
        real(wp), intent(in)                 :: xmin, ymin, zmin, length_x, length_y, length_z
        integer, intent(in)                  :: nx_bins, ny_bins, nz_bins
        logical, intent(out)                 :: overlaps
        integer, intent(out)                 :: bx, by, bz
        integer                              :: nbx, nby, nbz, slot
        integer                              :: dx_b, dy_b, dz_b, dz_lo, dz_hi, j
        real(wp)                             :: dist_sq, dx, dy, dz

        call s_get_cloud_bin(px, py, pz, min_dist, periodic_pack, xmin, ymin, zmin, nx_bins, ny_bins, nz_bins, bx, by, bz)

        dz_lo = -1
        dz_hi = 1
        if (num_dims < 3) then
            dz_lo = 0
            dz_hi = 0
        end if

        overlaps = .false.

        do dx_b = -1, 1
            do dy_b = -1, 1
                do dz_b = dz_lo, dz_hi
                    nbx = bx + dx_b
                    nby = by + dy_b
                    nbz = bz + dz_b
                    if (periodic_pack) then
                        nbx = modulo(nbx, nx_bins)
                        nby = modulo(nby, ny_bins)
                        if (num_dims == 3) nbz = modulo(nbz, nz_bins)
                    end if
                    slot = f_bin_hash(nbx, nby, nbz, hash_size)
                    j = hash_head(slot)
                    do while (j > 0)
                        dx = abs(px - placed(1, j))
                        dy = abs(py - placed(2, j))
                        if (periodic_pack) then
                            dx = min(dx, length_x - dx)
                            dy = min(dy, length_y - dy)
                        end if
                        if (num_dims < 3) then
                            dist_sq = dx**2 + dy**2
                        else
                            dz = abs(pz - placed(3, j))
                            if (periodic_pack) dz = min(dz, length_z - dz)
                            dist_sq = dx**2 + dy**2 + dz**2
                        end if
                        if (dist_sq < (p_radius + placed_radius(j) + gap)**2) then
                            overlaps = .true.
                            return
                        end if
                        j = chain_next(j)
                    end do
                end do
            end do
        end do

    end subroutine s_check_cloud_particle_overlap

    !> Xorshift PRNG. Advances seed in-place and returns a value in [0, 1).
    function f_xorshift(seed) result(rval)

        integer, intent(inout) :: seed
        real(wp)               :: rval

        seed = ieor(seed, ishft(seed, 13))
        seed = ieor(seed, ishft(seed, -17))
        seed = ieor(seed, ishft(seed, 5))

        ! Mask off the sign bit rather than abs(): at seed = -huge-1, abs() overflows and returns the value
        ! unchanged, yielding rval = -1 and, in the shell sampler, a NaN centre that then passes every overlap
        ! comparison (all comparisons against NaN are false) and gets placed.
        rval = real(iand(seed, huge(seed)), wp)/real(huge(seed), wp)

    end function f_xorshift

    !> Hash bin coordinates to a 1-indexed slot in [1, hash_size]. Uses large prime multipliers to spread bins across buckets. Hash
    !! collisions are benign: the distance check catches false neighbours.
    function f_bin_hash(bx, by, bz, hash_size) result(slot)

        integer, intent(in) :: bx, by, bz, hash_size
        integer             :: slot
        integer(8)          :: key

        key = ieor(ieor(int(bx, 8)*73856093_8, int(by, 8)*19349663_8), int(bz, 8)*83492791_8)
        slot = int(mod(abs(key), int(hash_size, 8))) + 1

    end function f_bin_hash

end module m_particle_cloud
