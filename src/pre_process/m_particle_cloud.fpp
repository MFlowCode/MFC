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

    !> Patch-shaped region of the current cloud, its world-space bounding box, and the axis a hemisphere shell opens toward
    type(ib_patch_parameters) :: cloud_region
    real(wp), dimension(3)    :: region_bbox_min, region_bbox_max
    integer                   :: region_shell_axis

    !> Spatial hash of the placed particles for rejection packing; hash_chain_next(i) links particle i to the previous one in its
    !! bucket
    integer, allocatable, dimension(:) :: hash_head, hash_chain_next
    integer                            :: hash_size
    integer, dimension(3)              :: hash_num_bins
    real(wp)                           :: hash_bin_width, hash_min_spacing
    logical                            :: hash_is_periodic

    !> Lattice of the current cloud: cell pitch per axis, cells per axis, and basis points per cell
    real(wp), dimension(3)    :: lattice_pitch
    real(wp), dimension(3, 4) :: lattice_basis
    integer, dimension(3)     :: lattice_num_cells
    integer                   :: lattice_num_basis

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
        integer                              :: num_requested, seed, alloc_stat
        integer(kind=8)                      :: num_attempts, max_attempts
        real(wp)                             :: particle_radius, pdf_max, estimated_pdf_max
        real(wp), dimension(3)               :: position
        real(wp), allocatable, dimension(:)  :: radii
        logical                              :: is_accepted

        num_requested = particle_cloud(cloud_idx)%num_particles
        allocate (points%pos(3, num_requested), points%vel(3, num_requested), points%radius(num_requested), &
                  & points%interface_velocity(num_requested), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call s_mpi_abort("Error :: Ran out of CPU memory trying to allocate particle cloud coordinates. " &
                             & // "Current system resources cannot place the specified number of particles.")
        end if
        points%num_particles = 0

        ! use a separate generation tree for file-based importing
        if (particle_cloud(cloud_idx)%packing_method == packing_method_file) then
            call s_read_particle_cloud_file(cloud_idx, points)
            return
        end if

        call s_set_cloud_region(cloud_idx)
        allocate (radii(num_requested))
        call s_draw_cloud_radii(cloud_idx, radii)
        seed = f_cloud_seed(cloud_idx)
        call s_initialize_cloud_packing(cloud_idx, radii(1), seed, pdf_max, max_attempts)
        estimated_pdf_max = pdf_max

        num_attempts = 0
        do while (points%num_particles < num_requested .and. num_attempts < max_attempts)
            num_attempts = num_attempts + 1

            ! radii are drawn up front, largest first, so rejected attempts do not bias the size distribution
            particle_radius = radii(points%num_particles + 1)

            ! get a point sampled via the packing method
            select case (particle_cloud(cloud_idx)%packing_method)
            case (packing_method_rejection)
                ! a uniform point in the bounding box that is clear of the placed particles
                call s_get_rejection_packing_point(seed, position)
                if (f_overlaps_placed_particle(position, particle_radius, points)) cycle
            case (packing_method_lattice)
                ! the num_attempts-th site of the lattice spanning the bounding box
                position = f_get_lattice_packing_point(num_attempts)
            end select

            ! check that the point fits in the cloud region
            if (.not. f_is_inside_cloud_region(cloud_idx, position)) cycle

            ! check the probability distribution; the default pdf accepts every point
            call s_check_cloud_pdf(cloud_idx, position, seed, pdf_max, is_accepted)
            if (.not. is_accepted) cycle

            ! this is a valid position
            call s_add_lagrangian_point(points, position, particle_radius)
            if (allocated(hash_head)) call s_insert_into_cloud_hash(points)
        end do

        if (points%num_particles < num_requested) call s_mpi_abort("Error :: Failed to place all particles in particle cloud")
        if (pdf_max > estimated_pdf_max .and. proc_rank == 0) print '(a,i0,a)', 'WARNING: particle_cloud(', cloud_idx, &
            & ')%pdf exceeded its sampled maximum; the raised bound slightly biases earlier placements.'
        if (allocated(hash_head)) deallocate (hash_head, hash_chain_next)

    end subroutine s_get_cloud_points

    !> Placement seed of cloud_idx; a zero seed is replaced, since xorshift stays at zero
    function f_cloud_seed(cloud_idx) result(seed)

        integer, intent(in) :: cloud_idx
        integer             :: seed

        seed = particle_cloud(cloud_idx)%seed
        if (seed == 0) seed = int(mod(1_8 + int(cloud_idx, 8)*1013904223_8, int(huge(seed), 8)))

    end function f_cloud_seed

    !> Estimates the pdf maximum, then sets up the packing method and caps the placement attempts
    subroutine s_initialize_cloud_packing(cloud_idx, largest_radius, seed, pdf_max, max_attempts)

        integer, intent(in)          :: cloud_idx, seed
        real(wp), intent(in)         :: largest_radius
        real(wp), intent(out)        :: pdf_max
        integer(kind=8), intent(out) :: max_attempts
        real(wp)                     :: acceptance

        call s_estimate_cloud_pdf(cloud_idx, seed, pdf_max, acceptance)

        select case (particle_cloud(cloud_idx)%packing_method)
        case (packing_method_rejection)
            ! 1000 surviving draws per particle
            max_attempts = int(min(real(particle_cloud(cloud_idx)%num_particles, wp)*1000._wp/acceptance, &
                               & real(huge(max_attempts), wp)), 8)
            call s_initialize_cloud_hash(cloud_idx, largest_radius)
        case (packing_method_lattice)
            ! every lattice site is tried once
            call s_initialize_cloud_lattice(cloud_idx, seed, pdf_max)
            max_attempts = f_num_lattice_sites()
        case default
            call s_mpi_abort("Particle cloud packing method is not a known packing method of MFC. Exiting.")
        end select

    end subroutine s_initialize_cloud_packing

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

    !> Is a particle centred at position inside the computational domain and cloud_idx's bounding box and region?
    function f_is_inside_cloud_region(cloud_idx, position) result(is_valid)

        integer, intent(in)                :: cloud_idx
        real(wp), dimension(3), intent(in) :: position
        logical                            :: is_valid
        real(wp), dimension(3)             :: xyz_local
        real(wp)                           :: radial_distance

        is_valid = position(1) >= x_domain_glb%beg .and. position(1) <= x_domain_glb%end .and. position(2) >= y_domain_glb%beg &
                            & .and. position(2) <= y_domain_glb%end
        if (num_dims == 3) is_valid = is_valid .and. position(3) >= z_domain_glb%beg .and. position(3) <= z_domain_glb%end
        is_valid = is_valid .and. all(position(1:num_dims) >= region_bbox_min(1:num_dims)) .and. all(position(1:num_dims) &
                                      & <= region_bbox_max(1:num_dims))
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

    !> Draws a uniform point in the current cloud's bounding box for random sequential addition, Widom (1966)
    subroutine s_get_rejection_packing_point(seed, position)

        integer, intent(inout)              :: seed
        real(wp), dimension(3), intent(out) :: position
        integer                             :: axis
        real(wp)                            :: random_draw

        position = region_bbox_min
        do axis = 1, num_dims
            call s_draw_xorshift(seed, random_draw)
            position(axis) = region_bbox_min(axis) + random_draw*(region_bbox_max(axis) - region_bbox_min(axis))
        end do

    end subroutine s_get_rejection_packing_point

    !> Sizes the spatial hash so neighbouring bins span the largest pair separation; periodic box clouds wrap their bins
    subroutine s_initialize_cloud_hash(cloud_idx, largest_radius)

        integer, intent(in)  :: cloud_idx
        real(wp), intent(in) :: largest_radius

        hash_min_spacing = particle_cloud(cloud_idx)%min_spacing
        hash_bin_width = 2._wp*largest_radius + hash_min_spacing
        hash_is_periodic = particle_cloud(cloud_idx)%cloud_geometry == cloud_geometry_box &
                                          & .and. particle_cloud(cloud_idx)%periodic == 1
        hash_num_bins = 1
        if (hash_is_periodic) then
            hash_num_bins(1:num_dims) = max(1, ceiling((region_bbox_max(1:num_dims) - region_bbox_min(1:num_dims))/hash_bin_width))
        end if

        ! 4x overprovisioned for ~25% load factor, minimum 16 buckets
        hash_size = max(16, 4*particle_cloud(cloud_idx)%num_particles)
        allocate (hash_head(hash_size), hash_chain_next(particle_cloud(cloud_idx)%num_particles))
        hash_head = -1
        hash_chain_next = -1

    end subroutine s_initialize_cloud_hash

    !> Spatial-hash bin of position, wrapped for periodic clouds
    function f_get_cloud_bin(position) result(bin)

        real(wp), dimension(3), intent(in) :: position
        integer, dimension(3)              :: bin

        bin = 0
        bin(1:num_dims) = floor((position(1:num_dims) - region_bbox_min(1:num_dims))/hash_bin_width)
        if (hash_is_periodic) bin = modulo(bin, hash_num_bins)

    end function f_get_cloud_bin

    !> Does a particle of radius at position overlap a placed particle in the 3x3(x3) neighbouring hash bins?
    function f_overlaps_placed_particle(position, radius, points) result(overlaps)

        real(wp), dimension(3), intent(in)  :: position
        real(wp), intent(in)                :: radius
        type(lagrangian_points), intent(in) :: points
        logical                             :: overlaps
        integer, dimension(3)               :: bin, neighbor_bin, neighbor_reach
        integer                             :: i, j, k, placed_index
        real(wp), dimension(3)              :: separation, periodic_lengths

        bin = f_get_cloud_bin(position)
        periodic_lengths = region_bbox_max - region_bbox_min
        neighbor_reach = 0
        neighbor_reach(1:num_dims) = 1
        overlaps = .false.

        do k = -neighbor_reach(3), neighbor_reach(3)
            do j = -1, 1
                do i = -1, 1
                    neighbor_bin = bin + [i, j, k]
                    if (hash_is_periodic) neighbor_bin = modulo(neighbor_bin, hash_num_bins)
                    placed_index = hash_head(f_bin_hash(neighbor_bin))
                    do while (placed_index > 0)
                        separation = abs(position - points%pos(:,placed_index))
                        if (hash_is_periodic) separation = min(separation, periodic_lengths - separation)
                        if (sum(separation(1:num_dims)**2) < (radius + points%radius(placed_index) + hash_min_spacing)**2) then
                            overlaps = .true.
                            return
                        end if
                        placed_index = hash_chain_next(placed_index)
                    end do
                end do
            end do
        end do

    end function f_overlaps_placed_particle

    !> Links the newest particle of points in as the head of its hash bucket
    subroutine s_insert_into_cloud_hash(points)

        type(lagrangian_points), intent(in) :: points
        integer                             :: slot

        slot = f_bin_hash(f_get_cloud_bin(points%pos(:,points%num_particles)))
        hash_chain_next(points%num_particles) = hash_head(slot)
        hash_head(slot) = points%num_particles

    end subroutine s_insert_into_cloud_hash

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
                    call s_draw_xorshift(seed, radial_draw)
                    call s_draw_xorshift(seed, angle_draw)
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

    !> Rejection-samples position against cloud_idx's pdf, von Neumann (1951); a pdf above pdf_max raises it
    subroutine s_check_cloud_pdf(cloud_idx, position, seed, pdf_max, is_accepted)

        integer, intent(in)                :: cloud_idx
        real(wp), dimension(3), intent(in) :: position
        integer, intent(inout)             :: seed
        real(wp), intent(inout)            :: pdf_max
        logical, intent(out)               :: is_accepted
        logical                            :: has_pdf
        real(wp)                           :: pdf, random_draw

        is_accepted = .true.
        call s_cloud_pdf(cloud_idx, position(1), position(2), position(3), has_pdf, pdf)
        if (.not. has_pdf) return

        @:PROHIBIT(.not. (pdf >= 0._wp), "particle_cloud pdf is negative or NaN inside its region")
        if (pdf > pdf_max) pdf_max = 1.1_wp*pdf
        call s_draw_xorshift(seed, random_draw)
        is_accepted = random_draw*pdf_max <= pdf

    end subroutine s_check_cloud_pdf

    !> Estimates the share of bounding-box draws a candidate survives (region, then pdf) and the padded pdf maximum (1 without one)
    subroutine s_estimate_cloud_pdf(cloud_idx, seed, pdf_max, acceptance)

        integer, intent(in)        :: cloud_idx, seed
        real(wp), intent(out)      :: pdf_max, acceptance
        integer, parameter         :: num_pdf_samples = 100000
        integer(kind=8), parameter :: max_draws = 100_8*num_pdf_samples
        integer                    :: num_samples, pdf_seed
        integer(kind=8)            :: num_draws
        real(wp)                   :: pdf, pdf_sum
        real(wp), dimension(3)     :: position
        logical                    :: has_pdf

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
            call s_get_rejection_packing_point(pdf_seed, position)
            if (.not. f_is_inside_cloud_region(cloud_idx, position)) cycle
            num_samples = num_samples + 1
            if (.not. has_pdf) cycle

            call s_cloud_pdf(cloud_idx, position(1), position(2), position(3), has_pdf, pdf)
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

    !> Sets the largest lattice spacing whose region holds at least num_particles sites after pdf thinning, bisecting where the site
    !! count crosses num_particles (the count is not monotone, so a larger spacing may exist)
    subroutine s_initialize_cloud_lattice(cloud_idx, seed, pdf_max)

        integer, intent(in)  :: cloud_idx, seed
        real(wp), intent(in) :: pdf_max
        integer              :: num_requested, i
        real(wp)             :: min_distance, spacing_lower, spacing_upper, spacing_middle

        num_requested = particle_cloud(cloud_idx)%num_particles
        min_distance = 2._wp*particle_cloud(cloud_idx)%radius + particle_cloud(cloud_idx)%min_spacing

        ! halve the spacing from the region size until enough sites fit
        spacing_upper = maxval(region_bbox_max(1:num_dims) - region_bbox_min(1:num_dims))
        spacing_lower = spacing_upper
        do while (f_count_lattice_sites(cloud_idx, spacing_lower, seed, pdf_max) < num_requested)
            @:PROHIBIT(spacing_lower <= min_distance, &
                       & "Particle cloud is too dense for lattice packing in its region; " &
                       & // "reduce num_particles or min_spacing, or enlarge the region")
            spacing_upper = spacing_lower
            spacing_lower = max(0.5_wp*spacing_lower, min_distance)
        end do

        ! bisect between too few and enough sites
        do i = 1, 60
            if (spacing_upper - spacing_lower <= 1.e-10_wp*spacing_upper) exit
            spacing_middle = 0.5_wp*(spacing_lower + spacing_upper)
            if (f_count_lattice_sites(cloud_idx, spacing_middle, seed, pdf_max) >= num_requested) then
                spacing_lower = spacing_middle
            else
                spacing_upper = spacing_middle
            end if
        end do

        call s_set_cloud_lattice(spacing_lower)

    end subroutine s_initialize_cloud_lattice

    !> Counts the lattice sites at spacing that the placement loop would accept, replaying its region and pdf checks
    function f_count_lattice_sites(cloud_idx, spacing, seed, pdf_max) result(num_valid_sites)

        integer, intent(in)    :: cloud_idx, seed
        real(wp), intent(in)   :: spacing, pdf_max
        integer                :: num_valid_sites, site_seed
        integer(kind=8)        :: site_index
        real(wp)               :: site_pdf_max
        real(wp), dimension(3) :: position
        logical                :: is_accepted

        call s_set_cloud_lattice(spacing)
        site_seed = seed
        site_pdf_max = pdf_max
        num_valid_sites = 0
        do site_index = 1, f_num_lattice_sites()
            position = f_get_lattice_packing_point(site_index)
            if (.not. f_is_inside_cloud_region(cloud_idx, position)) cycle
            call s_check_cloud_pdf(cloud_idx, position, site_seed, site_pdf_max, is_accepted)
            if (.not. is_accepted) cycle
            num_valid_sites = num_valid_sites + 1
        end do

    end function f_count_lattice_sites

    !> Lays the densest lattice of nearest-neighbour spacing over the bounding box, Conway & Sloane (1999) ch. 1
    subroutine s_set_cloud_lattice(spacing)

        real(wp), intent(in) :: spacing

        lattice_basis = 0._wp
        if (num_dims < 3) then
            ! triangular: rows pitched by spacing*sqrt(3)/2, odd rows shifted by half a spacing
            lattice_pitch = [spacing, spacing*sqrt(3._wp)/2._wp, 0._wp]
            lattice_num_basis = 1
        else
            ! face-centred cubic: a cubic cell of side spacing*sqrt(2) with four basis points
            lattice_pitch = spacing*sqrt(2._wp)
            lattice_num_basis = 4
            lattice_basis = reshape([0._wp, 0._wp, 0._wp, 0.5_wp, 0.5_wp, 0._wp, 0.5_wp, 0._wp, 0.5_wp, 0._wp, 0.5_wp, &
                                    & 0.5_wp]*lattice_pitch(1), [3, 4])
        end if

        lattice_num_cells = 1
        lattice_num_cells(1:num_dims) = floor((region_bbox_max(1:num_dims) - region_bbox_min(1:num_dims)) &
                          & /lattice_pitch(1:num_dims)) + 1

    end subroutine s_set_cloud_lattice

    !> Number of sites in the current lattice
    function f_num_lattice_sites() result(num_sites)

        integer(kind=8) :: num_sites

        num_sites = int(lattice_num_basis, 8)*product(int(lattice_num_cells, 8))

    end function f_num_lattice_sites

    !> Returns the site_index-th site of the current lattice, counting basis points fastest, then x, y and z cells
    function f_get_lattice_packing_point(site_index) result(position)

        integer(kind=8), intent(in) :: site_index
        real(wp), dimension(3)      :: position
        integer(kind=8)             :: remaining_index
        integer, dimension(3)       :: cell
        integer                     :: basis_index

        remaining_index = site_index - 1
        basis_index = int(mod(remaining_index, int(lattice_num_basis, 8))) + 1
        remaining_index = remaining_index/lattice_num_basis
        cell(1) = int(mod(remaining_index, int(lattice_num_cells(1), 8)))
        remaining_index = remaining_index/lattice_num_cells(1)
        cell(2) = int(mod(remaining_index, int(lattice_num_cells(2), 8)))
        cell(3) = int(remaining_index/lattice_num_cells(2))

        position = region_bbox_min + real(cell, wp)*lattice_pitch + lattice_basis(:,basis_index)
        if (num_dims < 3 .and. mod(cell(2), 2) == 1) position(1) = position(1) + 0.5_wp*lattice_pitch(1)

    end function f_get_lattice_packing_point

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

    !> Xorshift PRNG. Advances seed in-place and draws random_draw in [0, 1).
    subroutine s_draw_xorshift(seed, random_draw)

        integer, intent(inout) :: seed
        real(wp), intent(out)  :: random_draw

        seed = ieor(seed, ishft(seed, 13))
        seed = ieor(seed, ishft(seed, -17))
        seed = ieor(seed, ishft(seed, 5))

        ! mask the sign bit; abs() overflows at -huge-1
        random_draw = real(iand(seed, huge(seed)), wp)/real(huge(seed), wp)

    end subroutine s_draw_xorshift

    !> Hash bin coordinates to a 1-indexed slot in [1, hash_size]. Uses large prime multipliers to spread bins across buckets. Hash
    !! collisions are benign: the distance check catches false neighbours.
    function f_bin_hash(bin) result(slot)

        integer, dimension(3), intent(in) :: bin
        integer                           :: slot
        integer(kind=8)                   :: key

        key = ieor(ieor(int(bin(1), 8)*73856093_8, int(bin(2), 8)*19349663_8), int(bin(3), 8)*83492791_8)
        slot = int(mod(abs(key), int(hash_size, 8))) + 1

    end function f_bin_hash

end module m_particle_cloud
