!>
!! @file
!! @brief Contains module m_collisions

#:include 'macros.fpp'

!> @brief Ghost-node immersed boundary method: locates ghost/image points, computes interpolation coefficients, and corrects the flow state
module m_collisions

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_constants

    use m_compute_levelset

    use m_ib_patches

    use m_model

    implicit none

    private; public :: s_apply_collision_forces;

    ! arrays for holding collision lookups indices into the ib_overlap_distance array
    integer, dimension(num_ibs) :: col_idx_lookup_1
    integer, dimension(num_ibs) :: col_idx_lookup_2

    ! overlap distances for computing collisions
    integer, dimension(num_ibs * (num_ibs-1) / 2) :: collision_lookup
    real(wp), dimension(num_ibs, (num_dims*2)) :: wall_overlap_distances


contains

    subroutine s_apply_collision_forces(ghost_points, ib_markers)

      type(ghost_point), dimension(num_gps), intent(in) :: ghost_points
      type(integer_field), intent(in) :: ib_markers

      integer :: num_considered_collisions

      ! return if no collisions
      if (collision_model == 0) return

      $:GPU_UPDATE(host='[ghost_points]') ! TODO :: TEMPORARY UNTIL GPU SUPPORT ENABLED. REMOVE LATER

      ! get is distance used in the force calculation with each IB and each wall 
      call s_detect_ib_wall_collisions(ghost_points, wall_overlap_distances)
      call s_detect_ib_ib_collisions(ghost_points, ib_markers, collision_lookup, num_considered_collisions)

      select case (collision_model)
          case(1) ! soft sphere model
              call s_apply_wall_collision_forces_soft_sphere(wall_overlap_distances)
              call s_appply_ib_ib_collision_forces_soft_sphere(collision_lookup, num_considered_collisions)
      end select

    end subroutine s_apply_collision_forces

    !> @brief applyies collision forces to IBs assuming a soft-sphere collision model (all IBs are circles or spheres)
    subroutine s_apply_wall_collision_forces_soft_sphere(collision_lookup, num_considered_collisions)

      integer, intent(in) :: num_considered_collisions
      integer, dimension(num_considered_collisions), intent(in) :: collision_lookup

      end subroutine s_apply_wall_collision_forces_soft_sphere

    !> @brief applyies collision forces to IBs assuming a soft-sphere collision model (all IBs are circles or spheres)
    subroutine s_apply_wall_collision_forces_soft_sphere(wall_overlap_distances)

        real(wp), dimension(num_ibs, (num_dims*2)), intent(in) :: wall_overlap_distances

        integer :: patch_id, i
        real(wp) :: spring_stiffness, damping_parameter, coefficient_of_friction
        real(wp), dimension(3) :: normal_force, tangental_force, normal_vector, normal_velocity, tangental_vector

        spring_stiffness = 1._wp
        damping_parameter = 1._wp
        coefficient_of_friction = 1._wp

        do patch_id = 1, num_ibs
            do i = 1, num_dims*2
                ! only compute force contributions if there was an overlap
                if (f_approx_equal(wall_overlap_distances(patch_id, i), 0._wp)) cycle

                select case (i)
                    case(1) ! x domain left
                        normal_vector = [-1._wp, 0._wp, 0._wp]
                    case(2) ! x domain right
                        normal_vector = [1._wp, 0._wp, 0._wp]
                    case(3) ! y domain bottom
                        normal_vector = [0._wp, -1._wp, 0._wp]
                    case(4) ! y domain top
                        normal_vector = [0._wp, 1._wp, 0._wp]
                    case(5) ! z domain back
                        normal_vector = [0._wp, 0._wp, 1._wp]
                    case(6) ! z domain front
                        normal_vector = [0._wp, 0._wp, -1._wp]
                end select

                normal_velocity = dot_product(patch_ib(patch_id)%vel, normal_vector)*normal_vector
                tangental_vector = patch_ib(patch_id)%vel - normal_velocity
                tangental_vector = tangental_vector / norm2(tangental_vector)
                normal_force = -spring_stiffness * wall_overlap_distances(patch_id, i)  * normal_vector - damping_parameter * normal_velocity
                tangental_force = -coefficient_of_friction * norm2(normal_force) * tangental_vector
                patch_ib(patch_id)%forces = patch_ib(patch_id)%forces + normal_force + tangental_force
            end do
        end do


    end subroutine s_apply_wall_collision_forces_soft_sphere()

    !> uses ghost-point/image-point information to determine if it is possible if two IBs are colliding, effectively an optimized nearest neighbor search
    subroutine s_detect_ib_ib_collisions(gps, ib_markers, ib_collision_distances, collision_lookup, num_considered_collisions)

      type(ghost_point), dimension(num_gps), intent(in) :: gps
      type(integer_field), intent(in) :: ib_markers
      integer, dimension(num_ibs * (num_ibs-1) / 2), intent(out) :: collision_lookup
      integer, intent(out) :: num_considered_collisions

      integer :: i, j, k, col_idx_1, col_idx_2
      integer gp_idx, gp_patch_id, ip_patch_id
      integer :: max_pairs, pair_idx, out_idx
      logical :: already_found

      ! Temporary array to hold all detected pairs (with potential duplicates)
      integer, dimension(num_gps, 2) :: raw_pairs
      integer :: num_raw

      max_pairs = num_ibs * (num_ibs - 1) / 2
      num_raw = 0

      do gp_idx = 1, num_ibs

          gp_patch_id = gps(i)%ib_patch_id
          i = gps(i)%loc(1)
          h = gps(i)%loc(2)
          k = gps(i)%loc(3)
          ip_patch_id = ib_markers%sf(i, j, k)

          ! Pass 1: Collect all candidate pairs (may contain duplicates)
          if (ib_patch_id < ip_patch_id) then
              num_raw = num_raw + 1
              ! Store with smaller ID first for consistent ordering
              raw_pairs(num_raw, 1) = gp_patch_id
              raw_pairs(num_raw, 2) = ip_patch_id
          end if
      end do

      ! Pass 2: Coalesce into unique pairs
      num_considered_collisions = 0
      collision_lookup = 0
      do pair_idx = 1, num_raw
        already_found = .false.
        do out_idx = 1, num_considered_collisions
            if (collision_lookup(out_idx, 1) == raw_pairs(pair_idx, 1) .and. &
                collision_lookup(out_idx, 2) == raw_pairs(pair_idx, 2)) then
                already_found = .true.
                exit
            end if
        end do

        if (.not. already_found) then
            num_considered_collisions = num_considered_collisions + 1
            collision_lookup(num_considered_collisions, 1) = raw_pairs(pair_idx, 1)
            collision_lookup(num_considered_collisions, 2) = raw_pairs(pair_idx, 2)
        end if
    end do

    end subroutine s_detect_ib_ib_collisions

    !> @brief uses boundary conditions and particle lcoations to check for wall conditions
    subroutine s_detect_ib_wall_collisions(gps, overlap_distances)

        type(ghost_point), dimension(num_gps), intent(in) :: gps
        real(wp), dimension(num_ibs, (num_dims*2) + 1), intent(out) :: overlap_distances

        integer :: gp_idx, i, j, k, patch_id
        real(wp) :: edge_location, overlap_distance

        collision_distances = 0._wp

        ! iterate over all ghost points to detect the one that is most-overlapping in each direction
        do gp_idx = 1, num_gps
            i = gps(gp_idx)%loc(1)
            j = gps(gp_idx)%loc(2)
            if (p /= 0) k = gps(gp_idx)%loc(3)
            patch_id = gps(gp_idx)%ib_patch_id

            ! check if the boundaries are either of the two conditions we should compute collisions with
            if (bc_x%beg == BC_SLIP_WALL .or. bc_x%beg == BC_NO_SLIP_WALL) then
                ! get the location of the true IB surface towards the domain boundary
                edge_location = x_cc(i) - abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(1))
                ! check if that edge actually extends out of the comutational domain
                if (edge_location < x_domain%beg) then
                    overlap_distance = x_domain%beg - edge_location ! the distance that the IB extends out of the domain
                    overlap_distances(patch_id, 1) = max(overlap_distances(patch_id, 1), overlap_distance) ! Save this distance if it is the new maximum distance
                end if
            end if

            if (bc_x%end == BC_SLIP_WALL .or. bc_x%end == BC_NO_SLIP_WALL) then
                edge_location = x_cc(i) + abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(1))
                if (edge_location > x_domain%end) then
                    overlap_distance = edge_location - x_domain%end
                    overlap_distances(patch_id, 2) = max(overlap_distances(patch_id, 2), overlap_distance)
                end if
            end if

            if (bc_y%beg == BC_SLIP_WALL .or. bc_y%beg == BC_NO_SLIP_WALL) then
                edge_location = y_cc(j) - abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(2))
                if (edge_location < y_domain%beg) then
                    overlap_distance = y_domain%beg - edge_location
                    overlap_distances(patch_id, 3) = max(overlap_distances(patch_id, 3), overlap_distance)
                end if
            end if

            if (bc_y%end == BC_SLIP_WALL .or. bc_y%end == BC_NO_SLIP_WALL) then
                edge_location = y_cc(j) + abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(2))
                if (edge_location > y_domain%end) then
                    overlap_distance = edge_location - y_domain%end
                    overlap_distances(patch_id, 4) = max(overlap_distances(patch_id, 4), overlap_distance)
                end if
            end if

            if (p > 0) then
                if (bc_z%beg == BC_SLIP_WALL .or. bc_z%beg == BC_NO_SLIP_WALL) then
                    edge_location = z_cc(k) - abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(3))
                    if (edge_location < z_domain%beg) then
                        overlap_distance = z_domain%beg - edge_location
                        overlap_distances(patch_id, 6) = max(overlap_distances(patch_id, 5), overlap_distance)
                    end if
                end if

                if (bc_z%end == BC_SLIP_WALL .or. bc_z%end == BC_NO_SLIP_WALL) then
                    edge_location = z_cc(k) + abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(3))
                    if (edge_location > z_domain%end) then
                        overlap_distance = edge_location - z_domain%end
                        overlap_distances(patch_id, 6) = max(overlap_distances(patch_id, 6), overlap_distance)
                    end if
                end if
            end if
        end do

    end subroutine s_detect_ib_wall_collisions

end module m_collisions
