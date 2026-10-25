!>
!! @file m_particle_cloud.fpp
!! @brief Generates particle clouds as immersed boundary patches or Euler-Lagrange bubbles

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Places each cloud's particle centres, then turns them into IB patches or Euler-Lagrange bubbles. Placement is
!! deterministic, so under file_per_process every rank places every IB and keeps its own; EL clouds are placed on rank 0 only.
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

    !> Patch-shaped region of the current cloud and its world-space bounding box (the box itself for a box cloud)
    type(ib_patch_parameters) :: cloud_region
    real(wp), dimension(3)    :: region_bbox_min, region_bbox_max

contains

    !> Generates every cloud into particle_cloud_ibs (IB) and lag_bubbles (EL); only the first num_particle_cloud_ibs and
    !! num_lag_bubbles entries are set
    impure subroutine s_generate_particle_clouds(glb_bounds, particle_cloud_ibs, num_particle_cloud_ibs, lag_bubbles, &
        & num_lag_bubbles)

        type(bounds_info), dimension(3), intent(in)                       :: glb_bounds
        type(ib_patch_parameters), allocatable, intent(out), dimension(:) :: particle_cloud_ibs
        integer, intent(out)                                              :: num_particle_cloud_ibs
        real(wp), allocatable, intent(out), dimension(:,:)                :: lag_bubbles
        integer, intent(out)                                              :: num_lag_bubbles
        type(particle_cloud_coords)                                       :: coords
        integer                                                           :: cloud_idx, global_index
        real(wp)                                                          :: start_time, end_time

        allocate (particle_cloud_ibs(min(num_ib_patches_max_namelist, f_num_cloud_particles(particle_type_ib))))
        allocate (lag_bubbles(num_lag_bubble_input_vars, f_num_cloud_particles(particle_type_el)))
        num_particle_cloud_ibs = 0
        num_lag_bubbles = 0
        if (num_particle_clouds == 0) return

        call cpu_time(start_time)
        call s_initialize_cloud_shapes()

        global_index = num_ibs
        do cloud_idx = 1, num_particle_clouds
            if (particle_cloud(cloud_idx)%particle_type == particle_type_el .and. proc_rank /= 0) cycle

            call s_get_cloud_coords(cloud_idx, coords)
            select case (particle_cloud(cloud_idx)%particle_type)
            case (particle_type_ib)
                call s_coords_to_cloud_ibs(coords, glb_bounds, global_index, particle_cloud_ibs, num_particle_cloud_ibs)
            case (particle_type_el)
                call s_coords_to_lag_bubbles(coords, lag_bubbles, num_lag_bubbles)
            case default
                call s_mpi_abort("Particle cloud particle_type is not a known particle type of MFC. Exiting.")
            end select
            deallocate (coords%pos, coords%radius)
        end do

        call cpu_time(end_time)
        if (proc_rank == 0) print '(a,i0,a,i0,a,f0.3,a)', 'Particle clouds placed ', global_index - num_ibs, ' IBs and ', &
            & num_lag_bubbles, ' EL bubbles in ', end_time - start_time, ' seconds.'

    end subroutine s_generate_particle_clouds

    !> Fills coords with cloud_idx's particle centres, in the order global ids are assigned
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
        if (particle_cloud(cloud_idx)%cloud_geometry /= cloud_geometry_hemisphere_shell) call s_set_cloud_region(cloud_idx)

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

    !> Loads the STL models and airfoil surfaces that patch-shaped clouds test against
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

    !> Builds cloud_region, an IB patch with cloud_idx's shape, and its world-space bounding box
    subroutine s_set_cloud_region(cloud_idx)

        integer, intent(in)    :: cloud_idx
        real(wp), dimension(3) :: centroid, lengths

        centroid = [particle_cloud(cloud_idx)%x_centroid, particle_cloud(cloud_idx)%y_centroid, &
                                   & particle_cloud(cloud_idx)%z_centroid]
        lengths = [particle_cloud(cloud_idx)%length_x, particle_cloud(cloud_idx)%length_y, particle_cloud(cloud_idx)%length_z]
        if (particle_cloud(cloud_idx)%cloud_geometry == cloud_geometry_box) then
            region_bbox_min = centroid - 0.5_wp*lengths
            region_bbox_max = centroid + 0.5_wp*lengths
        else
            cloud_region%geometry = particle_cloud(cloud_idx)%cloud_geometry
            cloud_region%x_centroid = centroid(1)
            cloud_region%y_centroid = centroid(2)
            cloud_region%z_centroid = centroid(3)
            cloud_region%length_x = lengths(1)
            cloud_region%length_y = lengths(2)
            cloud_region%length_z = lengths(3)
            cloud_region%radius = particle_cloud(cloud_idx)%region_radius
            cloud_region%angles = particle_cloud(cloud_idx)%angles
            cloud_region%airfoil_id = particle_cloud(cloud_idx)%airfoil_id
            cloud_region%model_id = particle_cloud(cloud_idx)%model_id
            cloud_region%centroid_offset = 0._wp
            call s_update_patch_rotation_matrix(cloud_region)
            call s_get_patch_bounding_box(cloud_region, centroid, region_bbox_min, region_bbox_max)
        end if
        if (num_dims < 3) then
            region_bbox_min(3) = centroid(3)
            region_bbox_max(3) = centroid(3)
        end if

    end subroutine s_set_cloud_region

    !> Is position inside both the current cloud region and the computational domain?
    function f_is_inside_cloud_region(position) result(is_valid)

        real(wp), dimension(3), intent(in) :: position
        logical                            :: is_valid
        real(wp), dimension(3)             :: xyz_local

        is_valid = position(1) >= x_domain_glb%beg .and. position(1) <= x_domain_glb%end .and. position(2) >= y_domain_glb%beg &
                            & .and. position(2) <= y_domain_glb%end
        if (num_dims == 3) is_valid = is_valid .and. position(3) >= z_domain_glb%beg .and. position(3) <= z_domain_glb%end
        if (.not. is_valid) return

        xyz_local = position - [cloud_region%x_centroid, cloud_region%y_centroid, cloud_region%z_centroid]
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

    !> Rejection-samples cloud_idx's centres, largest first, each pair at least their radii plus min_spacing apart
    subroutine s_particle_cloud_rejection_pack(cloud_idx, coords)

        integer, intent(in)                        :: cloud_idx
        type(particle_cloud_coords), intent(inout) :: coords
        integer                                    :: num_requested, seed
        integer(8)                                 :: num_attempts, max_attempts
        real(wp)                                   :: bin_width, candidate_x, candidate_y, candidate_z, pdf, pdf_max, acceptance
        logical                                    :: overlaps, reject, periodic_pack, has_pdf
        integer                                    :: hash_size, slot
        integer                                    :: bin_x, bin_y, bin_z, nx_bins, ny_bins, nz_bins
        integer, allocatable                       :: hash_head(:), chain_next(:)
        real(wp)                                   :: xmin, ymin, zmin, length_x, length_y, length_z, next_radius
        real(wp), allocatable                      :: radii(:)

        num_requested = particle_cloud(cloud_idx)%num_particles
        allocate (radii(num_requested))
        call s_draw_cloud_radii(cloud_idx, radii)
        ! largest pair separation sets the hash bin width
        bin_width = 2._wp*radii(1) + particle_cloud(cloud_idx)%min_spacing
        periodic_pack = particle_cloud(cloud_idx)%cloud_geometry == cloud_geometry_box .and. particle_cloud(cloud_idx)%periodic == 1
        length_x = particle_cloud(cloud_idx)%length_x
        length_y = particle_cloud(cloud_idx)%length_y
        length_z = particle_cloud(cloud_idx)%length_z
        xmin = particle_cloud(cloud_idx)%x_centroid - 0.5_wp*length_x
        ymin = particle_cloud(cloud_idx)%y_centroid - 0.5_wp*length_y
        zmin = particle_cloud(cloud_idx)%z_centroid - 0.5_wp*length_z
        ! only periodic packing uses the bin counts
        nx_bins = 1
        ny_bins = 1
        nz_bins = 1
        if (periodic_pack) then
            nx_bins = max(1, ceiling(length_x/bin_width))
            ny_bins = max(1, ceiling(length_y/bin_width))
            if (num_dims == 3) nz_bins = max(1, ceiling(length_z/bin_width))
        end if

        num_attempts = 0
        seed = particle_cloud(cloud_idx)%seed
        if (seed == 0) seed = 1 + cloud_idx*1013904223

        call s_estimate_cloud_pdf(cloud_idx, seed, has_pdf, pdf_max, acceptance)
        max_attempts = int(min(real(num_requested, wp)*1000._wp/acceptance, real(huge(max_attempts), wp)), 8)

        ! Hash table: 4x overprovisioned for ~25% load factor, minimum 16 buckets. chain_next(i) links placed particle i to the
        ! previous occupant of its bucket.
        hash_size = max(16, 4*num_requested)
        allocate (hash_head(hash_size))
        allocate (chain_next(num_requested))
        hash_head = -1
        chain_next = -1

        ! random sequential addition, Widom (1966)
        do while (coords%num_particles < num_requested .and. num_attempts < max_attempts)
            num_attempts = num_attempts + 1

            next_radius = radii(coords%num_particles + 1)
            call s_sample_cloud_candidate(cloud_idx, seed, next_radius, candidate_x, candidate_y, candidate_z, reject)
            if (reject) cycle

            ! rejection sampling against the padded pdf maximum, von Neumann (1951)
            if (has_pdf) then
                call s_cloud_pdf(cloud_idx, candidate_x, candidate_y, candidate_z, has_pdf, pdf)
                @:PROHIBIT(.not. (pdf >= 0._wp), "particle_cloud pdf is negative or NaN inside its region")
                if (pdf > pdf_max) then
                    if (proc_rank == 0) print '(a,i0,a)', 'WARNING: particle_cloud(', cloud_idx, &
                        & ')%pdf exceeds its sampled maximum; raising the bound, which slightly biases earlier placements.'
                    pdf_max = 1.1_wp*pdf
                end if
                if (f_xorshift(seed)*pdf_max > pdf) cycle
            end if

            call s_check_cloud_particle_overlap(candidate_x, candidate_y, candidate_z, next_radius, coords%pos, coords%radius, &
                                                & hash_head, chain_next, hash_size, bin_width, &
                                                & particle_cloud(cloud_idx)%min_spacing, periodic_pack, xmin, ymin, zmin, &
                                                & length_x, length_y, length_z, nx_bins, ny_bins, nz_bins, overlaps, bin_x, &
                                                & bin_y, bin_z)

            if (.not. overlaps) then
                call s_add_cloud_coord(coords, [candidate_x, candidate_y, candidate_z], next_radius)

                ! Insert into hash grid as head of bucket chain
                slot = f_bin_hash(bin_x, bin_y, bin_z, hash_size)
                chain_next(coords%num_particles) = hash_head(slot)
                hash_head(slot) = coords%num_particles
            end if
        end do

        if (coords%num_particles < num_requested) then
            call s_mpi_abort("Error :: Failed to place all particles in particle bed")
        end if

        deallocate (hash_head, chain_next, radii)

    end subroutine s_particle_cloud_rejection_pack

    !> Draws cloud_idx's particle radii, largest first: one size, exact size-class fractions, or a clipped lognormal
    subroutine s_draw_cloud_radii(cloud_idx, radii)

        integer, intent(in)                        :: cloud_idx
        real(wp), dimension(:), intent(out)        :: radii
        integer                                    :: num_particles, k, i, seed, num_assigned, num_in_class
        integer, dimension(num_radius_classes_max) :: class_order
        real(wp)                                   :: min_radius, max_radius, drawn_radius, cumulative_fraction, total_fraction
        real(wp)                                   :: radial_draw, angle_draw

        num_particles = size(radii)
        radii = particle_cloud(cloud_idx)%radius
        call s_cloud_radius_bounds(particle_cloud(cloud_idx), min_radius, max_radius)

        if (.not. f_is_default(particle_cloud(cloud_idx)%radius_sigma)) then
            ! own seed stream, so the placement stream is untouched
            seed = ieor(particle_cloud(cloud_idx)%seed, 2135587861)
            if (seed == 0) seed = 1
            do i = 1, num_particles
                ! lognormal radius from a Box-Muller normal deviate, Box & Muller (1958)
                do
                    radial_draw = f_xorshift(seed)
                    angle_draw = f_xorshift(seed)
                    drawn_radius = particle_cloud(cloud_idx)%radius*exp(particle_cloud(cloud_idx)%radius_sigma*sqrt(-2._wp*log(max(radial_draw, &
                                                  & tiny(1._wp))))*cos(2._wp*pi*angle_draw))
                    if (drawn_radius >= min_radius .and. drawn_radius <= max_radius) exit
                end do
                radii(i) = drawn_radius
            end do
            call s_sort_descending(radii)
        else if (any(particle_cloud(cloud_idx)%radius_fraction > 0._wp)) then
            ! largest class first; cumulative rounding makes the counts sum to num_particles
            total_fraction = sum(particle_cloud(cloud_idx)%radius_fraction)
            class_order = [(k, k=1, num_radius_classes_max)]
            do k = 1, num_radius_classes_max - 1
                do i = k + 1, num_radius_classes_max
                    if (particle_cloud(cloud_idx)%radius_ratio(class_order(i)) &
                        & > particle_cloud(cloud_idx)%radius_ratio(class_order(k))) then
                        class_order([k, i]) = class_order([i, k])
                    end if
                end do
            end do
            num_assigned = 0
            cumulative_fraction = 0._wp
            do k = 1, num_radius_classes_max
                cumulative_fraction = cumulative_fraction + particle_cloud(cloud_idx)%radius_fraction(class_order(k))/total_fraction
                num_in_class = nint(cumulative_fraction*real(num_particles, wp)) - num_assigned
                radii(num_assigned + 1:num_assigned + num_in_class) &
                      & = particle_cloud(cloud_idx)%radius*particle_cloud(cloud_idx)%radius_ratio(class_order(k))
                num_assigned = num_assigned + num_in_class
            end do
        end if

    end subroutine s_draw_cloud_radii

    !> Sorts values in place, largest first (heapsort, Williams 1964).
    subroutine s_sort_descending(values)

        real(wp), dimension(:), intent(inout) :: values
        integer                               :: num_values, i, root, child
        real(wp)                              :: temporary

        num_values = size(values)
        do i = num_values/2, 1, -1
            call s_sift(i, num_values)
        end do
        do i = num_values, 2, -1
            temporary = values(1)
            values(1) = values(i)
            values(i) = temporary
            call s_sift(1, i - 1)
        end do
        values = values(num_values:1:-1)

    contains

        !> Restores the max-heap below start within values(1:last).
        subroutine s_sift(start, last)

            integer, intent(in) :: start, last

            root = start
            do while (2*root <= last)
                child = 2*root
                if (child < last) then
                    if (values(child + 1) > values(child)) child = child + 1
                end if
                if (values(root) >= values(child)) return
                temporary = values(root)
                values(root) = values(child)
                values(child) = temporary
                root = child
            end do

        end subroutine s_sift

    end subroutine s_sort_descending

    !> Evaluates cloud_idx's analytic placement pdf at (x, y, z), compiled in through case.fpp; pdf is 1 without one
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

    !> Estimates the padded pdf maximum and the acceptance rate over the region; both are 1 without a pdf
    subroutine s_estimate_cloud_pdf(cloud_idx, seed, has_pdf, pdf_max, acceptance)

        integer, intent(in)   :: cloud_idx, seed
        logical, intent(out)  :: has_pdf
        real(wp), intent(out) :: pdf_max, acceptance
        integer, parameter    :: num_pdf_samples = 100000
        integer               :: num_samples, pdf_seed
        real(wp)              :: candidate_x, candidate_y, candidate_z, pdf, pdf_sum
        logical               :: reject

        pdf_max = 1._wp
        acceptance = 1._wp
        call s_cloud_pdf(cloud_idx, 0._wp, 0._wp, 0._wp, has_pdf, pdf)
        if (.not. has_pdf) return

        ! own seed stream, so the placement stream is untouched
        pdf_seed = ieor(seed, 1540483477)
        if (pdf_seed == 0) pdf_seed = 1
        num_samples = 0
        pdf_sum = 0._wp
        pdf_max = 0._wp
        do while (num_samples < num_pdf_samples)
            call s_sample_cloud_candidate(cloud_idx, pdf_seed, particle_cloud(cloud_idx)%radius, candidate_x, candidate_y, &
                                          & candidate_z, reject)
            if (reject) cycle
            call s_cloud_pdf(cloud_idx, candidate_x, candidate_y, candidate_z, has_pdf, pdf)
            @:PROHIBIT(.not. (pdf >= 0._wp), "particle_cloud pdf is negative or NaN inside its region")
            num_samples = num_samples + 1
            pdf_sum = pdf_sum + pdf
            pdf_max = max(pdf_max, pdf)
        end do

        @:PROHIBIT(pdf_max <= 0._wp, "particle_cloud pdf is zero everywhere inside its region")
        pdf_max = 1.1_wp*pdf_max
        acceptance = pdf_sum/(real(num_samples, wp)*pdf_max)

    end subroutine s_estimate_cloud_pdf

    !> Draws a uniform candidate centre for cloud_idx, rejecting it near a shell's flat face or outside a patch shape. Shells use
    !! the sqrt (2D) or cube-root (3D) radial CDF; in 2D any shell_axis but 1 opens toward +y.
    subroutine s_sample_cloud_candidate(cloud_idx, seed, particle_radius, candidate_x, candidate_y, candidate_z, reject)

        integer, intent(in)    :: cloud_idx
        integer, intent(inout) :: seed
        real(wp), intent(in)   :: particle_radius
        real(wp), intent(out)  :: candidate_x, candidate_y, candidate_z
        logical, intent(out)   :: reject
        real(wp)               :: theta, phi, r_shell, rho, u, zdir, r_inner, r_outer

        reject = .false.

        select case (particle_cloud(cloud_idx)%cloud_geometry)
        case (cloud_geometry_hemisphere_shell)
            r_inner = particle_cloud(cloud_idx)%shell_inner_radius + particle_radius
            r_outer = particle_cloud(cloud_idx)%shell_outer_radius - particle_radius

            if (num_dims < 3) then
                theta = pi*f_xorshift(seed)
                u = f_xorshift(seed)
                r_shell = sqrt((r_outer**2 - r_inner**2)*u + r_inner**2)
                candidate_z = particle_cloud(cloud_idx)%z_centroid
                if (particle_cloud(cloud_idx)%shell_axis == 1) then
                    candidate_x = particle_cloud(cloud_idx)%x_centroid + r_shell*sin(theta)
                    candidate_y = particle_cloud(cloud_idx)%y_centroid + r_shell*cos(theta)
                    if (candidate_x < particle_cloud(cloud_idx)%x_centroid + particle_radius) reject = .true.
                else
                    candidate_x = particle_cloud(cloud_idx)%x_centroid + r_shell*cos(theta)
                    candidate_y = particle_cloud(cloud_idx)%y_centroid + r_shell*sin(theta)
                    if (candidate_y < particle_cloud(cloud_idx)%y_centroid + particle_radius) reject = .true.
                end if
            else
                phi = 2._wp*pi*f_xorshift(seed)
                zdir = f_xorshift(seed)
                rho = sqrt(max(0._wp, 1._wp - zdir**2))
                u = f_xorshift(seed)
                r_shell = ((r_outer**3 - r_inner**3)*u + r_inner**3)**(1._wp/3._wp)
                select case (particle_cloud(cloud_idx)%shell_axis)
                case (1)  ! opens toward +x
                    candidate_x = particle_cloud(cloud_idx)%x_centroid + r_shell*zdir
                    candidate_y = particle_cloud(cloud_idx)%y_centroid + r_shell*rho*cos(phi)
                    candidate_z = particle_cloud(cloud_idx)%z_centroid + r_shell*rho*sin(phi)
                    if (candidate_x < particle_cloud(cloud_idx)%x_centroid + particle_radius) reject = .true.
                case (2)  ! opens toward +y
                    candidate_y = particle_cloud(cloud_idx)%y_centroid + r_shell*zdir
                    candidate_x = particle_cloud(cloud_idx)%x_centroid + r_shell*rho*cos(phi)
                    candidate_z = particle_cloud(cloud_idx)%z_centroid + r_shell*rho*sin(phi)
                    if (candidate_y < particle_cloud(cloud_idx)%y_centroid + particle_radius) reject = .true.
                case default  ! 3: opens toward +z
                    candidate_z = particle_cloud(cloud_idx)%z_centroid + r_shell*zdir
                    candidate_x = particle_cloud(cloud_idx)%x_centroid + r_shell*rho*cos(phi)
                    candidate_y = particle_cloud(cloud_idx)%y_centroid + r_shell*rho*sin(phi)
                    if (candidate_z < particle_cloud(cloud_idx)%z_centroid + particle_radius) reject = .true.
                end select
            end if
        case default  ! box or patch shape, uniform in its bounding box
            candidate_x = region_bbox_min(1) + f_xorshift(seed)*(region_bbox_max(1) - region_bbox_min(1))
            candidate_y = region_bbox_min(2) + f_xorshift(seed)*(region_bbox_max(2) - region_bbox_min(2))
            candidate_z = region_bbox_min(3)
            if (num_dims == 3) candidate_z = region_bbox_min(3) + f_xorshift(seed)*(region_bbox_max(3) - region_bbox_min(3))
            if (f_is_patch_cloud(cloud_idx)) reject = .not. f_is_inside_cloud_region([candidate_x, candidate_y, candidate_z])
        end select

    end subroutine s_sample_cloud_candidate

    !> Places particles on the optimally dense lattice for the cloud region: a triangular lattice in 2D, a face-centered cubic
    !! lattice in 3D. The lattice spacing is set by the particle density (num_particles over the region area/volume); if that
    !! spacing falls below the required centre-to-centre distance (2*radius + min_spacing), the region is too dense and the run is
    !! aborted.
    subroutine s_particle_cloud_lattice(cloud_idx, coords)

        integer, intent(in)                        :: cloud_idx
        type(particle_cloud_coords), intent(inout) :: coords
        integer                                    :: num_requested, row, col, ncx, ncy, ix, jy, kz, b
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
        num_requested = particle_cloud(cloud_idx)%num_particles

        if (num_dims < 3) then
            ! Triangular lattice: area per particle = (sqrt(3)/2)*spacing**2.
            spacing = sqrt(2._wp*(xmax - xmin)*(ymax - ymin)/(sqrt(3._wp)*real(num_requested, wp)))
        else
            ! Face-centered cubic lattice: volume per particle = spacing**3/sqrt(2).
            spacing = (sqrt(2._wp)*(xmax - xmin)*(ymax - ymin)*(zmax - zmin)/real(num_requested, wp))**(1._wp/3._wp)
        end if

        if (spacing < min_dist) then
            call s_mpi_abort("Error :: Particle cloud is too dense for lattice packing; " &
                             & // "reduce num_particles or min_spacing, or enlarge the cloud region")
        end if

        if (num_dims < 3) then
            ! Triangular lattice: rows pitched by spacing*sqrt(3)/2, odd rows shifted by half a spacing.
            row_dy = spacing*sqrt(3._wp)/2._wp
            row = 0
            do while (coords%num_particles < num_requested)
                py = ymin + real(row, wp)*row_dy
                x0 = xmin
                if (mod(row, 2) == 1) x0 = xmin + 0.5_wp*spacing
                col = 0
                px = x0
                do while (px <= xmax .and. coords%num_particles < num_requested)
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
            do while (coords%num_particles < num_requested)
                do jy = 0, ncy - 1
                    do ix = 0, ncx - 1
                        do b = 1, 4
                            if (coords%num_particles >= num_requested) exit
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

    !> Places num_particles centres on the densest lattice inside a patch shape, bisecting the spacing to where the site count
    !! crosses num_particles (the count is not monotone, so a larger spacing may exist); surplus sites are dropped evenly.
    subroutine s_particle_cloud_lattice_region(cloud_idx, coords)

        integer, intent(in)                        :: cloud_idx
        type(particle_cloud_coords), intent(inout) :: coords
        integer                                    :: num_requested, num_sites_lower, num_sites_middle, i
        real(wp)                                   :: min_distance, spacing_lower, spacing_upper, spacing_middle
        real(wp)                                   :: pdf_max, acceptance
        logical                                    :: has_pdf

        num_requested = particle_cloud(cloud_idx)%num_particles
        min_distance = 2._wp*particle_cloud(cloud_idx)%radius + particle_cloud(cloud_idx)%min_spacing
        call s_estimate_cloud_pdf(cloud_idx, particle_cloud(cloud_idx)%seed, has_pdf, pdf_max, acceptance)

        ! halve the spacing from the region size until enough sites fit
        spacing_upper = maxval(region_bbox_max(1:num_dims) - region_bbox_min(1:num_dims))
        spacing_lower = spacing_upper
        num_sites_lower = f_region_lattice_sites(cloud_idx, spacing_lower, has_pdf, pdf_max)
        do while (num_sites_lower < num_requested)
            @:PROHIBIT(spacing_lower <= min_distance, &
                       & "Particle cloud is too dense for lattice packing in its region; " &
                       & // "reduce num_particles or min_spacing, or enlarge the region")
            spacing_upper = spacing_lower
            spacing_lower = max(0.5_wp*spacing_lower, min_distance)
            num_sites_lower = f_region_lattice_sites(cloud_idx, spacing_lower, has_pdf, pdf_max)
        end do

        ! bisect between too few and enough sites
        do i = 1, 60
            if (spacing_upper - spacing_lower <= 1.e-10_wp*spacing_upper) exit
            spacing_middle = 0.5_wp*(spacing_lower + spacing_upper)
            num_sites_middle = f_region_lattice_sites(cloud_idx, spacing_middle, has_pdf, pdf_max)
            if (num_sites_middle >= num_requested) then
                spacing_lower = spacing_middle
                num_sites_lower = num_sites_middle
            else
                spacing_upper = spacing_middle
            end if
        end do

        num_sites_middle = f_region_lattice_sites(cloud_idx, spacing_lower, has_pdf, pdf_max, num_sites_lower, coords)

    end subroutine s_particle_cloud_lattice_region

    !> Counts the current region's lattice sites at spacing; given num_sites and coords, also appends num_particles of them
    function f_region_lattice_sites(cloud_idx, spacing, has_pdf, pdf_max, num_sites, coords) result(num_valid_sites)

        integer, intent(in)                                  :: cloud_idx
        real(wp), intent(in)                                 :: spacing, pdf_max
        logical, intent(in)                                  :: has_pdf
        integer, intent(in), optional                        :: num_sites
        type(particle_cloud_coords), intent(inout), optional :: coords
        integer                                              :: num_valid_sites
        integer                                              :: seed, num_columns, num_rows, num_layers, i, j, k, l
        integer                                              :: num_basis_points
        real(wp)                                             :: row_spacing, cell_length, pdf
        real(wp), dimension(3)                               :: site
        real(wp), dimension(3, 4)                            :: basis
        logical                                              :: has_site_pdf

        ! same pdf stream every call, so the count and the fill thin the same sites

        seed = ieor(particle_cloud(cloud_idx)%seed, 1013904223)
        if (seed == 0) seed = 1
        num_valid_sites = 0

        ! densest lattice packings of disks and spheres, Conway & Sloane (1999) ch. 1
        if (num_dims < 3) then
            ! Triangular lattice: rows pitched by spacing*sqrt(3)/2, odd rows shifted by half a spacing.
            row_spacing = spacing*sqrt(3._wp)/2._wp
            num_columns = floor((region_bbox_max(1) - region_bbox_min(1))/spacing) + 1
            num_rows = floor((region_bbox_max(2) - region_bbox_min(2))/row_spacing) + 1
            num_layers = 1
            num_basis_points = 1
            basis(:,1) = 0._wp
        else
            ! Face-centered cubic lattice via the conventional cubic cell (side = spacing*sqrt(2)) and its four basis points.
            cell_length = spacing*sqrt(2._wp)
            num_columns = floor((region_bbox_max(1) - region_bbox_min(1))/cell_length) + 1
            num_rows = floor((region_bbox_max(2) - region_bbox_min(2))/cell_length) + 1
            num_layers = floor((region_bbox_max(3) - region_bbox_min(3))/cell_length) + 1
            num_basis_points = 4
            basis = reshape([0._wp, 0._wp, 0._wp, 0.5_wp, 0.5_wp, 0._wp, 0.5_wp, 0._wp, 0.5_wp, 0._wp, 0.5_wp, &
                            & 0.5_wp]*cell_length, [3, 4])
        end if

        do k = 0, num_layers - 1
            do j = 0, num_rows - 1
                do i = 0, num_columns - 1
                    do l = 1, num_basis_points
                        if (num_dims < 3) then
                            site = region_bbox_min + [real(i, wp)*spacing + merge(0.5_wp*spacing, 0._wp, mod(j, 2) == 1), real(j, &
                                                           & wp)*row_spacing, 0._wp]
                        else
                            site = region_bbox_min + real([i, j, k], wp)*cell_length + basis(:,l)
                        end if
                        if (any(site(1:num_dims) > region_bbox_max(1:num_dims))) cycle
                        if (.not. f_is_inside_cloud_region(site)) cycle
                        if (has_pdf) then
                            call s_cloud_pdf(cloud_idx, site(1), site(2), site(3), has_site_pdf, pdf)
                            if (f_xorshift(seed)*pdf_max > pdf) cycle
                        end if

                        num_valid_sites = num_valid_sites + 1
                        if (present(coords)) then
                            ! keep the site each time valid_sites*num_particles/num_sites steps up, dropping the surplus evenly
                            if ((int(num_valid_sites, &
                                & 8)*particle_cloud(cloud_idx)%num_particles)/num_sites > (int(num_valid_sites - 1, &
                                & 8)*particle_cloud(cloud_idx)%num_particles)/num_sites) then
                                call s_add_cloud_coord(coords, site, particle_cloud(cloud_idx)%radius)
                            end if
                        end if
                    end do
                end do
            end do
        end do

    end function f_region_lattice_sites

    !> Reads cloud_idx's "x y z" centre rows from its filepath, skipping blank and '#' lines; 2D reads x y and uses z_centroid
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

    !> Appends coords to particle_cloud_ibs with the next global ids, keeping only this rank's under file_per_process. Only the
    !! fields the IB state file holds are set; simulation sets the rest in s_assign_particle_cloud_ib_defaults.
    subroutine s_coords_to_cloud_ibs(coords, glb_bounds, global_index, particle_cloud_ibs, num_particle_cloud_ibs)

        type(particle_cloud_coords), intent(in)                :: coords
        type(bounds_info), dimension(3), intent(in)            :: glb_bounds
        integer, intent(inout)                                 :: global_index, num_particle_cloud_ibs
        type(ib_patch_parameters), intent(inout), dimension(:) :: particle_cloud_ibs
        integer                                                :: i

        do i = 1, coords%num_particles
            global_index = global_index + 1
            if (file_per_process) then
                if (.not. f_local_rank_owns_location(coords%pos(:,i), glb_bounds)) cycle
            end if
            num_particle_cloud_ibs = num_particle_cloud_ibs + 1
            @:PROHIBIT(num_particle_cloud_ibs > size(particle_cloud_ibs), &
                       & "Too many particle-cloud IBs on one rank. Modify case file or increase num_ib_patches_max_namelist.")
            particle_cloud_ibs(num_particle_cloud_ibs)%gbl_patch_id = global_index
            particle_cloud_ibs(num_particle_cloud_ibs)%x_centroid = coords%pos(1, i)
            particle_cloud_ibs(num_particle_cloud_ibs)%y_centroid = coords%pos(2, i)
            particle_cloud_ibs(num_particle_cloud_ibs)%z_centroid = coords%pos(3, i)
            particle_cloud_ibs(num_particle_cloud_ibs)%vel(:) = 0._wp
            particle_cloud_ibs(num_particle_cloud_ibs)%angular_vel(:) = 0._wp
            particle_cloud_ibs(num_particle_cloud_ibs)%angles(:) = 0._wp
            particle_cloud_ibs(num_particle_cloud_ibs)%radius = coords%radius(i)
        end do

    end subroutine s_coords_to_cloud_ibs

    !> Appends coords to lag_bubbles as Euler-Lagrange bubbles at rest, each at its own equilibrium radius
    subroutine s_coords_to_lag_bubbles(coords, lag_bubbles, num_lag_bubbles)

        type(particle_cloud_coords), intent(in) :: coords
        real(wp), intent(inout), dimension(:,:) :: lag_bubbles
        integer, intent(inout)                  :: num_lag_bubbles
        integer                                 :: i

        do i = 1, coords%num_particles
            num_lag_bubbles = num_lag_bubbles + 1
            lag_bubbles(:,num_lag_bubbles) = 0._wp
            lag_bubbles(1:3,num_lag_bubbles) = coords%pos(:,i)
            lag_bubbles(7, num_lag_bubbles) = coords%radius(i)  ! equilibrium radius R0
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

    !> Does the candidate overlap a placed particle in the 3x3(x3) neighbouring hash bins? Also returns the candidate's own bin.
    subroutine s_check_cloud_particle_overlap(candidate_x, candidate_y, candidate_z, candidate_radius, placed, placed_radius, &
        & hash_head, chain_next, hash_size, bin_width, min_spacing, periodic_pack, xmin, ymin, zmin, length_x, length_y, &
        & length_z, nx_bins, ny_bins, nz_bins, overlaps, bin_x, bin_y, bin_z)

        real(wp), intent(in)                 :: candidate_x, candidate_y, candidate_z, candidate_radius, bin_width, min_spacing
        real(wp), intent(in), dimension(:,:) :: placed
        real(wp), intent(in), dimension(:)   :: placed_radius
        integer, intent(in), dimension(:)    :: hash_head, chain_next
        integer, intent(in)                  :: hash_size
        logical, intent(in)                  :: periodic_pack
        real(wp), intent(in)                 :: xmin, ymin, zmin, length_x, length_y, length_z
        integer, intent(in)                  :: nx_bins, ny_bins, nz_bins
        logical, intent(out)                 :: overlaps
        integer, intent(out)                 :: bin_x, bin_y, bin_z
        integer                              :: neighbor_bin_x, neighbor_bin_y, neighbor_bin_z, slot
        integer                              :: i, j, k, k_lower, k_upper, placed_index
        real(wp)                             :: distance_squared, dx, dy, dz

        call s_get_cloud_bin(candidate_x, candidate_y, candidate_z, bin_width, periodic_pack, xmin, ymin, zmin, nx_bins, ny_bins, &
                             & nz_bins, bin_x, bin_y, bin_z)

        k_lower = -1
        k_upper = 1
        if (num_dims < 3) then
            k_lower = 0
            k_upper = 0
        end if

        overlaps = .false.

        do i = -1, 1
            do j = -1, 1
                do k = k_lower, k_upper
                    neighbor_bin_x = bin_x + i
                    neighbor_bin_y = bin_y + j
                    neighbor_bin_z = bin_z + k
                    if (periodic_pack) then
                        neighbor_bin_x = modulo(neighbor_bin_x, nx_bins)
                        neighbor_bin_y = modulo(neighbor_bin_y, ny_bins)
                        if (num_dims == 3) neighbor_bin_z = modulo(neighbor_bin_z, nz_bins)
                    end if
                    slot = f_bin_hash(neighbor_bin_x, neighbor_bin_y, neighbor_bin_z, hash_size)
                    placed_index = hash_head(slot)
                    do while (placed_index > 0)
                        dx = abs(candidate_x - placed(1, placed_index))
                        dy = abs(candidate_y - placed(2, placed_index))
                        if (periodic_pack) then
                            dx = min(dx, length_x - dx)
                            dy = min(dy, length_y - dy)
                        end if
                        if (num_dims < 3) then
                            distance_squared = dx**2 + dy**2
                        else
                            dz = abs(candidate_z - placed(3, placed_index))
                            if (periodic_pack) dz = min(dz, length_z - dz)
                            distance_squared = dx**2 + dy**2 + dz**2
                        end if
                        if (distance_squared < (candidate_radius + placed_radius(placed_index) + min_spacing)**2) then
                            overlaps = .true.
                            return
                        end if
                        placed_index = chain_next(placed_index)
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
