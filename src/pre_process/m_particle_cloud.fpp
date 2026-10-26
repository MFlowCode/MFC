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

    !> Lagrangian points of one cloud, independent of the particle type they become; rows of a particle cloud file
    type lagrangian_points
        integer                               :: num_particles       !< Number of points stored
        real(wp), allocatable, dimension(:,:) :: pos                 !< Centres, (1:3, 1:capacity)
        real(wp), allocatable, dimension(:,:) :: vel                 !< Initial velocities, (1:3, 1:capacity)
        real(wp), allocatable, dimension(:)   :: radius              !< Radius of each particle
        real(wp), allocatable, dimension(:)   :: interface_velocity  !< Initial bubble wall velocity; zero for IBs
    end type lagrangian_points

    !> Patch-shaped region of the current cloud, its world-space bounding box, and the axis a hemisphere shell opens toward
    type(ib_patch_parameters) :: cloud_region
    real(wp), dimension(3)    :: region_bbox_min, region_bbox_max
    integer                   :: region_shell_axis

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
        type(lagrangian_points)                                           :: points
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

            call s_get_cloud_points(cloud_idx, points)
            select case (particle_cloud(cloud_idx)%particle_type)
            case (particle_type_ib)
                call s_points_to_cloud_ibs(points, glb_bounds, global_index, particle_cloud_ibs, num_particle_cloud_ibs)
            case (particle_type_el)
                call s_points_to_lag_bubbles(points, lag_bubbles, num_lag_bubbles)
            case default
                call s_mpi_abort("Particle cloud particle_type is not a known particle type of MFC. Exiting.")
            end select
            deallocate (points%pos, points%vel, points%radius, points%interface_velocity)
        end do

        call cpu_time(end_time)
        if (proc_rank == 0) print '(a,i0,a,i0,a,f0.3,a)', 'Particle clouds placed ', global_index - num_ibs, ' IBs and ', &
            & num_lag_bubbles, ' EL bubbles in ', end_time - start_time, ' seconds.'

    end subroutine s_generate_particle_clouds

    !> Fills points with cloud_idx's particles, in the order global ids are assigned
    impure subroutine s_get_cloud_points(cloud_idx, points)

        integer, intent(in)                  :: cloud_idx
        type(lagrangian_points), intent(out) :: points
        integer                              :: num_requested, alloc_stat

        num_requested = particle_cloud(cloud_idx)%num_particles
        allocate (points%pos(3, num_requested), points%vel(3, num_requested), points%radius(num_requested), &
                  & points%interface_velocity(num_requested), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call s_mpi_abort("Error :: Ran out of CPU memory trying to allocate particle cloud coordinates. " &
                             & // "Current system resources cannot place the specified number of particles.")
        end if
        points%num_particles = 0

        if (particle_cloud(cloud_idx)%cloud_geometry == cloud_geometry_file) then
            call s_read_particle_cloud_file(cloud_idx, points)
            return
        end if
        call s_set_cloud_region(cloud_idx)

        select case (particle_cloud(cloud_idx)%packing_method)
        case (1)  ! rejection (random) packing method
            call s_particle_cloud_rejection_pack(cloud_idx, points)
        case (2)  ! lattice packing method
            call s_particle_cloud_lattice(cloud_idx, points)
        case default
            call s_mpi_abort("Particle cloud packing method is not a known packing method of MFC. Exiting.")
        end select

    end subroutine s_get_cloud_points

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

    !> Sets the world-space bounding box of cloud_idx's region, and cloud_region for a patch shape
    subroutine s_set_cloud_region(cloud_idx)

        integer, intent(in)    :: cloud_idx
        real(wp), dimension(3) :: centroid, lengths

        centroid = [particle_cloud(cloud_idx)%x_centroid, particle_cloud(cloud_idx)%y_centroid, &
                                   & particle_cloud(cloud_idx)%z_centroid]
        lengths = [particle_cloud(cloud_idx)%length_x, particle_cloud(cloud_idx)%length_y, particle_cloud(cloud_idx)%length_z]

        select case (particle_cloud(cloud_idx)%cloud_geometry)
        case (cloud_geometry_box)
            region_bbox_min = centroid - 0.5_wp*lengths
            region_bbox_max = centroid + 0.5_wp*lengths
        case (cloud_geometry_hemisphere_shell)
            ! in 2D, any shell_axis but 1 opens toward +y
            region_shell_axis = particle_cloud(cloud_idx)%shell_axis
            if (num_dims < 3 .and. region_shell_axis /= 1) region_shell_axis = 2
            region_bbox_min = centroid - particle_cloud(cloud_idx)%shell_outer_radius
            region_bbox_max = centroid + particle_cloud(cloud_idx)%shell_outer_radius
            region_bbox_min(region_shell_axis) = centroid(region_shell_axis)
        case default  ! patch shape
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
        end select

        ! 2D clouds lie in the z_centroid plane
        if (num_dims < 3) then
            region_bbox_min(3) = centroid(3)
            region_bbox_max(3) = centroid(3)
        end if

    end subroutine s_set_cloud_region

    !> Is a particle centred at position inside both the computational domain and cloud_idx's region?
    function f_is_inside_cloud_region(cloud_idx, position) result(is_valid)

        integer, intent(in)                :: cloud_idx
        real(wp), dimension(3), intent(in) :: position
        logical                            :: is_valid
        real(wp), dimension(3)             :: xyz_local
        real(wp)                           :: radial_distance

        is_valid = position(1) >= x_domain_glb%beg .and. position(1) <= x_domain_glb%end .and. position(2) >= y_domain_glb%beg &
                            & .and. position(2) <= y_domain_glb%end
        if (num_dims == 3) is_valid = is_valid .and. position(3) >= z_domain_glb%beg .and. position(3) <= z_domain_glb%end
        if (.not. is_valid) return

        xyz_local = position - [particle_cloud(cloud_idx)%x_centroid, particle_cloud(cloud_idx)%y_centroid, &
                                               & particle_cloud(cloud_idx)%z_centroid]
        select case (particle_cloud(cloud_idx)%cloud_geometry)
        case (cloud_geometry_box)  ! the bounding box is the region
        case (cloud_geometry_hemisphere_shell)  ! between the shell radii, on the open side of the flat face
            radial_distance = norm2(xyz_local(1:num_dims))
            is_valid = xyz_local(region_shell_axis) >= 0._wp
            is_valid = is_valid .and. radial_distance >= particle_cloud(cloud_idx)%shell_inner_radius
            is_valid = is_valid .and. radial_distance <= particle_cloud(cloud_idx)%shell_outer_radius
        case default  ! patch shape
            is_valid = f_is_inside_patch(cloud_region, matmul(cloud_region%rotation_matrix_inverse, xyz_local))
        end select

    end function f_is_inside_cloud_region

    !> Appends one particle to points, at rest unless a velocity or interface velocity is given
    subroutine s_add_lagrangian_point(points, centroid, radius, velocity, interface_velocity)

        type(lagrangian_points), intent(inout)       :: points
        real(wp), dimension(3), intent(in)           :: centroid
        real(wp), intent(in)                         :: radius
        real(wp), dimension(3), intent(in), optional :: velocity
        real(wp), intent(in), optional               :: interface_velocity

        points%num_particles = points%num_particles + 1
        points%pos(:,points%num_particles) = centroid
        points%radius(points%num_particles) = radius
        points%vel(:,points%num_particles) = 0._wp
        points%interface_velocity(points%num_particles) = 0._wp
        if (present(velocity)) points%vel(:,points%num_particles) = velocity
        if (present(interface_velocity)) points%interface_velocity(points%num_particles) = interface_velocity

    end subroutine s_add_lagrangian_point

    !> Rejection-samples cloud_idx's centres, largest first, each pair at least their radii plus min_spacing apart
    subroutine s_particle_cloud_rejection_pack(cloud_idx, points)

        integer, intent(in)                    :: cloud_idx
        type(lagrangian_points), intent(inout) :: points
        integer                                :: num_requested, seed, initial_seed
        integer(8)                             :: num_attempts, max_attempts, estimate_attempt
        real(wp)                               :: bin_width, candidate_x, candidate_y, candidate_z, pdf, pdf_max, acceptance
        logical                                :: overlaps, reject, periodic_pack, has_pdf
        integer                                :: hash_size, slot
        integer                                :: bin_x, bin_y, bin_z, nx_bins, ny_bins, nz_bins
        integer, allocatable                   :: hash_head(:), chain_next(:)
        real(wp)                               :: xmin, ymin, zmin, length_x, length_y, length_z, next_radius
        real(wp), allocatable                  :: radii(:)

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
        initial_seed = seed

        ! the cap is at least 1000 draws per particle, so without a pdf the costly estimate waits for them
        call s_cloud_pdf(cloud_idx, 0._wp, 0._wp, 0._wp, has_pdf, pdf)
        estimate_attempt = 1000_8*num_requested
        if (has_pdf) estimate_attempt = 0
        max_attempts = huge(max_attempts)

        ! Hash table: 4x overprovisioned for ~25% load factor, minimum 16 buckets. chain_next(i) links placed particle i to the
        ! previous occupant of its bucket.
        hash_size = max(16, 4*num_requested)
        allocate (hash_head(hash_size))
        allocate (chain_next(num_requested))
        hash_head = -1
        chain_next = -1

        ! random sequential addition, Widom (1966)
        do while (points%num_particles < num_requested .and. num_attempts < max_attempts)
            ! estimate the pdf maximum and the attempt cap
            if (num_attempts == estimate_attempt) then
                call s_estimate_cloud_pdf(cloud_idx, initial_seed, has_pdf, pdf_max, acceptance)
                max_attempts = int(min(real(num_requested, wp)*1000._wp/acceptance, real(huge(max_attempts), wp)), 8)
                if (num_attempts == max_attempts) exit
            end if
            num_attempts = num_attempts + 1

            next_radius = radii(points%num_particles + 1)
            call s_sample_cloud_candidate(cloud_idx, seed, candidate_x, candidate_y, candidate_z, reject)
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

            call s_check_cloud_particle_overlap(candidate_x, candidate_y, candidate_z, next_radius, points%pos, points%radius, &
                                                & hash_head, chain_next, hash_size, bin_width, &
                                                & particle_cloud(cloud_idx)%min_spacing, periodic_pack, xmin, ymin, zmin, &
                                                & length_x, length_y, length_z, nx_bins, ny_bins, nz_bins, overlaps, bin_x, &
                                                & bin_y, bin_z)

            if (.not. overlaps) then
                call s_add_lagrangian_point(points, [candidate_x, candidate_y, candidate_z], next_radius)

                ! Insert into hash grid as head of bucket chain
                slot = f_bin_hash(bin_x, bin_y, bin_z, hash_size)
                chain_next(points%num_particles) = hash_head(slot)
                hash_head(slot) = points%num_particles
            end if
        end do

        if (points%num_particles < num_requested) then
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
        real(wp)                                   :: radial_draw, angle_draw, log_radius_offset

        num_particles = size(radii)
        radii = particle_cloud(cloud_idx)%radius

        if (.not. f_is_default(particle_cloud(cloud_idx)%radius_sigma)) then
            ! clip to radius_min and radius_max, by default radius*exp(-/+3 radius_sigma)
            min_radius = merge(particle_cloud(cloud_idx)%radius*exp(-3._wp*particle_cloud(cloud_idx)%radius_sigma), &
                               & particle_cloud(cloud_idx)%radius_min, f_is_default(particle_cloud(cloud_idx)%radius_min))
            max_radius = merge(particle_cloud(cloud_idx)%radius*exp(3._wp*particle_cloud(cloud_idx)%radius_sigma), &
                               & particle_cloud(cloud_idx)%radius_max, f_is_default(particle_cloud(cloud_idx)%radius_max))
            ! own seed stream, so the placement stream is untouched
            seed = ieor(particle_cloud(cloud_idx)%seed, 2135587861)
            if (seed == 0) seed = 1
            do i = 1, num_particles
                ! lognormal radius from a Box-Muller normal deviate, Box & Muller (1958)
                do
                    radial_draw = f_xorshift(seed)
                    angle_draw = f_xorshift(seed)
                    log_radius_offset = particle_cloud(cloud_idx)%radius_sigma*sqrt(-2._wp*log(max(radial_draw, &
                                                       & tiny(1._wp))))*cos(2._wp*pi*angle_draw)
                    drawn_radius = particle_cloud(cloud_idx)%radius*exp(log_radius_offset)
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

    !> Estimates the share of bounding-box draws a candidate survives (region, then pdf) and the padded pdf maximum (1 without one)
    subroutine s_estimate_cloud_pdf(cloud_idx, seed, has_pdf, pdf_max, acceptance)

        integer, intent(in)   :: cloud_idx, seed
        logical, intent(out)  :: has_pdf
        real(wp), intent(out) :: pdf_max, acceptance
        integer, parameter    :: num_pdf_samples = 100000
        integer(8), parameter :: max_draws = 100_8*num_pdf_samples
        integer               :: num_samples, pdf_seed
        integer(8)            :: num_draws
        real(wp)              :: candidate_x, candidate_y, candidate_z, pdf, pdf_sum
        logical               :: reject

        pdf_max = 1._wp
        call s_cloud_pdf(cloud_idx, 0._wp, 0._wp, 0._wp, has_pdf, pdf)
        if (has_pdf) pdf_max = 0._wp

        ! own seed stream, so the placement stream is untouched
        pdf_seed = ieor(seed, 1540483477)
        if (pdf_seed == 0) pdf_seed = 1
        num_samples = 0
        num_draws = 0
        pdf_sum = 0._wp
        do while (num_samples < num_pdf_samples .and. num_draws < max_draws)
            num_draws = num_draws + 1
            call s_sample_cloud_candidate(cloud_idx, pdf_seed, candidate_x, candidate_y, candidate_z, reject)
            if (reject) cycle
            num_samples = num_samples + 1
            if (.not. has_pdf) cycle

            call s_cloud_pdf(cloud_idx, candidate_x, candidate_y, candidate_z, has_pdf, pdf)
            @:PROHIBIT(.not. (pdf >= 0._wp), "particle_cloud pdf is negative or NaN inside its region")
            pdf_sum = pdf_sum + pdf
            pdf_max = max(pdf_max, pdf)
        end do

        @:PROHIBIT(num_samples == 0, "No particle fits inside a particle_cloud region; enlarge the region or shrink radius")
        acceptance = real(num_samples, wp)/real(num_draws, wp)
        if (.not. has_pdf) return

        @:PROHIBIT(pdf_max <= 0._wp, "particle_cloud pdf is zero everywhere inside its region")
        pdf_max = 1.1_wp*pdf_max
        acceptance = acceptance*pdf_sum/(real(num_samples, wp)*pdf_max)

    end subroutine s_estimate_cloud_pdf

    !> Draws a uniform candidate centre in cloud_idx's bounding box, rejecting it outside the region, von Neumann (1951)
    subroutine s_sample_cloud_candidate(cloud_idx, seed, candidate_x, candidate_y, candidate_z, reject)

        integer, intent(in)    :: cloud_idx
        integer, intent(inout) :: seed
        real(wp), intent(out)  :: candidate_x, candidate_y, candidate_z
        logical, intent(out)   :: reject

        candidate_x = region_bbox_min(1) + f_xorshift(seed)*(region_bbox_max(1) - region_bbox_min(1))
        candidate_y = region_bbox_min(2) + f_xorshift(seed)*(region_bbox_max(2) - region_bbox_min(2))
        candidate_z = region_bbox_min(3)
        if (num_dims == 3) candidate_z = region_bbox_min(3) + f_xorshift(seed)*(region_bbox_max(3) - region_bbox_min(3))
        reject = .not. f_is_inside_cloud_region(cloud_idx, [candidate_x, candidate_y, candidate_z])

    end subroutine s_sample_cloud_candidate

    !> Places num_particles centres on the densest lattice inside the cloud region, bisecting the spacing to where the site count
    !! crosses num_particles (the count is not monotone, so a larger spacing may exist); surplus sites are dropped evenly.
    subroutine s_particle_cloud_lattice(cloud_idx, points)

        integer, intent(in)                    :: cloud_idx
        type(lagrangian_points), intent(inout) :: points
        integer                                :: num_requested, num_sites_lower, num_sites_middle, i
        real(wp)                               :: min_distance, spacing_lower, spacing_upper, spacing_middle
        real(wp)                               :: pdf_max, acceptance
        logical                                :: has_pdf

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

        num_sites_middle = f_region_lattice_sites(cloud_idx, spacing_lower, has_pdf, pdf_max, num_sites_lower, points)

    end subroutine s_particle_cloud_lattice

    !> Counts the current region's lattice sites at spacing; given num_sites and points, also appends num_particles of them
    function f_region_lattice_sites(cloud_idx, spacing, has_pdf, pdf_max, num_sites, points) result(num_valid_sites)

        integer, intent(in)                              :: cloud_idx
        real(wp), intent(in)                             :: spacing, pdf_max
        logical, intent(in)                              :: has_pdf
        integer, intent(in), optional                    :: num_sites
        type(lagrangian_points), intent(inout), optional :: points
        integer                                          :: num_valid_sites
        integer                                          :: seed, num_columns, num_rows, num_layers, i, j, k, l
        integer                                          :: num_basis_points
        real(wp)                                         :: row_spacing, cell_length, pdf
        real(wp), dimension(3)                           :: site
        real(wp), dimension(3, 4)                        :: basis
        integer(8)                                       :: num_kept, num_kept_previous
        logical                                          :: has_site_pdf

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
                        if (.not. f_is_inside_cloud_region(cloud_idx, site)) cycle
                        if (has_pdf) then
                            call s_cloud_pdf(cloud_idx, site(1), site(2), site(3), has_site_pdf, pdf)
                            if (f_xorshift(seed)*pdf_max > pdf) cycle
                        end if

                        num_valid_sites = num_valid_sites + 1
                        if (present(points)) then
                            ! keep the site each time valid_sites*num_particles/num_sites steps up, dropping the surplus evenly
                            num_kept = int(num_valid_sites, 8)*particle_cloud(cloud_idx)%num_particles/num_sites
                            num_kept_previous = int(num_valid_sites - 1, 8)*particle_cloud(cloud_idx)%num_particles/num_sites
                            if (num_kept > num_kept_previous) call s_add_lagrangian_point(points, site, &
                                & particle_cloud(cloud_idx)%radius)
                        end if
                    end do
                end do
            end do
        end do

    end function f_region_lattice_sites

    !> Reads cloud_idx's "x y z u v w radius" rows, plus "interface_velocity" for EL bubbles; 2D ignores z, '#' lines are skipped
    impure subroutine s_read_particle_cloud_file(cloud_idx, points)

        integer, intent(in)                    :: cloud_idx
        type(lagrangian_points), intent(inout) :: points
        character(LEN=pathlen_max)             :: line
        real(wp), dimension(9)                 :: row
        integer                                :: file_unit, ios, num_values, max_values

        open (newunit=file_unit, file=trim(particle_cloud(cloud_idx)%filepath), form='formatted', status='old', action='read', &
              & iostat=ios)
        if (ios /= 0) call s_mpi_abort("Cannot open particle cloud file: " // trim(particle_cloud(cloud_idx)%filepath))

        ! the input/lag_bubbles.dat layout; only bubbles carry an interface velocity
        max_values = merge(8, 7, particle_cloud(cloud_idx)%particle_type == particle_type_el)
        do
            read (file_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle

            ! count the row's values by reading the most that fit
            do num_values = max_values + 1, 7, -1
                read (line, *, iostat=ios) row(1:num_values)
                if (ios == 0) exit
            end do
            if (ios /= 0 .or. num_values > max_values) then
                call s_mpi_abort("Malformed row in particle cloud file (expected x y z u v w radius, plus " &
                                 & // "interface_velocity for EL bubbles): " // trim(particle_cloud(cloud_idx)%filepath))
            end if

            ! a missing interface velocity starts at rest
            if (num_values == 7) row(8) = 0._wp
            ! 2D bubbles keep z, which axisymmetric runs fold into their radial position
            if (num_dims < 3 .and. max_values == 7) row(3) = particle_cloud(cloud_idx)%z_centroid

            if (.not. (row(7) > 0._wp)) call s_mpi_abort("Particle cloud file radius must be > 0: " &
                & // trim(particle_cloud(cloud_idx)%filepath))
            if (points%num_particles == particle_cloud(cloud_idx)%num_particles) then
                call s_mpi_abort("Particle cloud file holds more rows than num_particles: " &
                                 & // trim(particle_cloud(cloud_idx)%filepath))
            end if

            ! position, radius, velocity, interface velocity
            call s_add_lagrangian_point(points, row(1:3), row(7), row(4:6), row(8))
        end do
        close (file_unit)

        if (points%num_particles /= particle_cloud(cloud_idx)%num_particles) then
            call s_mpi_abort("Particle cloud file holds fewer rows than num_particles: " &
                             & // trim(particle_cloud(cloud_idx)%filepath))
        end if

    end subroutine s_read_particle_cloud_file

    !> Appends points to particle_cloud_ibs with the next global ids, keeping only this rank's under file_per_process. Only the
    !! fields the IB state file holds are set; simulation sets the rest in s_assign_particle_cloud_ib_defaults.
    subroutine s_points_to_cloud_ibs(points, glb_bounds, global_index, particle_cloud_ibs, num_particle_cloud_ibs)

        type(lagrangian_points), intent(in)                    :: points
        type(bounds_info), dimension(3), intent(in)            :: glb_bounds
        integer, intent(inout)                                 :: global_index, num_particle_cloud_ibs
        type(ib_patch_parameters), intent(inout), dimension(:) :: particle_cloud_ibs
        integer                                                :: i

        do i = 1, points%num_particles
            global_index = global_index + 1
            if (file_per_process) then
                if (.not. f_local_rank_owns_location(points%pos(:,i), glb_bounds)) cycle
            end if
            num_particle_cloud_ibs = num_particle_cloud_ibs + 1
            @:PROHIBIT(num_particle_cloud_ibs > size(particle_cloud_ibs), &
                       & "Too many particle-cloud IBs on one rank. Modify case file or increase num_ib_patches_max_namelist.")
            particle_cloud_ibs(num_particle_cloud_ibs)%gbl_patch_id = global_index
            particle_cloud_ibs(num_particle_cloud_ibs)%x_centroid = points%pos(1, i)
            particle_cloud_ibs(num_particle_cloud_ibs)%y_centroid = points%pos(2, i)
            particle_cloud_ibs(num_particle_cloud_ibs)%z_centroid = points%pos(3, i)
            particle_cloud_ibs(num_particle_cloud_ibs)%vel(:) = points%vel(:,i)
            particle_cloud_ibs(num_particle_cloud_ibs)%angular_vel(:) = 0._wp
            particle_cloud_ibs(num_particle_cloud_ibs)%angles(:) = 0._wp
            particle_cloud_ibs(num_particle_cloud_ibs)%radius = points%radius(i)
        end do

    end subroutine s_points_to_cloud_ibs

    !> Appends points to lag_bubbles as Euler-Lagrange bubbles, each starting at its equilibrium radius
    subroutine s_points_to_lag_bubbles(points, lag_bubbles, num_lag_bubbles)

        type(lagrangian_points), intent(in)     :: points
        real(wp), intent(inout), dimension(:,:) :: lag_bubbles
        integer, intent(inout)                  :: num_lag_bubbles
        integer                                 :: i

        do i = 1, points%num_particles
            num_lag_bubbles = num_lag_bubbles + 1
            ! position, velocity, equilibrium radius R0, interface velocity
            lag_bubbles(:,num_lag_bubbles) = [points%pos(:,i), points%vel(:,i), points%radius(i), points%interface_velocity(i)]
        end do

    end subroutine s_points_to_lag_bubbles

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

        ! mask the sign bit; abs() overflows at -huge-1
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
