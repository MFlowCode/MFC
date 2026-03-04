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

    private; public :: s_apply_collision_forces, s_initialize_collisions_module;

    ! overlap distances for computing collisions
    integer, dimension(num_ibs * (num_ibs-1) / 2, 2) :: collision_lookup
    real(wp), dimension(num_ibs, (num_dims*2)) :: wall_overlap_distances

    real(wp) :: spring_stiffness, damping_parameter

contains

    subroutine s_initialize_collisions_module()

        real(wp) :: e, m

        e = coefficient_of_restitution

        damping_parameter = -2._wp * log(e) / (pi**2 + log(e)**2)
        spring_stiffness = ( collision_time**2 * (pi**2 + log(e)**2) )**-1

    end subroutine s_initialize_collisions_module

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
    subroutine s_apply_ib_collision_forces_soft_sphere(collision_lookup, num_considered_collisions)

        integer, intent(in) :: num_considered_collisions
        integer, dimension(num_considered_collisions), intent(in) :: collision_lookup

        integer :: i, pid1, pid2 ! iterators and patch IDs
        real(wp) :: overlap_distance
        real(wp), dimension(3) :: collision_location, normal_vector, centroid_1, centroid_2
        real(wp) :: k, eta, effective_mass ! the spring stiffness and damping coefficient and mass of a specific interaction

        ! Iterate over all collisions detected
        do i = 1, num_considered_collisions
            pid1 = collision_lookup(i, 1)
            pid2 = collision_lookup(i, 2)

            centroid_1 = [patch_ib(pid1)%x_centroid, patch_ib(pid1)%y_centroid, 0._wp]
            centroid_2 = [patch_ib(pid2)%x_centroid, patch_ib(pid2)%y_centroid, 0._wp]
            if (num_dims == 3) then 
                centroid_1(3) = patch_ib(pid1)%z_centroid
                centroid_2(3) = patch_ib(pid2)%z_centroid
            end if

            normal_vector = centroid_2 - centroid_1
            overlap_distance = patch_ib(pid1)%radius + patch_ib(pid1)%radius - norm2(normal_vector)
            if (overlap_distance > 0._wp) then ! if the two patches are close enough to collide
                normal_vector = normal_vector / norm2(normal_vector)
                collision_location = centroid_1 + normal_vector * (patch_ib(pid1)%radius - 0.5_wp*overlap_distance)
                if (f_local_rank_owns_collision(collision_location)) then

                    ! compute constants of the collision
                    effective_mass = ( patch_ib(pid1)%mass**-1 + patch_ib(pid2)%mass**-1 )**-1
                    k = spring_stiffness * effective_mass
                    eta = damping_parameter * sqrt(effective_mass * k)

                    ! standard soft-sphere collision model
                    normal_velocity = dot_product(patch_ib(pid1)%vel, normal_vector)*normal_vector
                    tangental_vector = patch_ib(pid1)%vel - normal_velocity
                    tangental_vector = tangental_vector / norm2(tangental_vector)
                    normal_force = -k * overlap_distance  * normal_vector - eta * normal_velocity
                    tangental_force = -coefficient_of_friction * norm2(normal_force) * tangental_vector
                    call s_cross_product(normal_vector * patch_ib(pid1)%radius, tangental_force, torque)
                    patch_ib(pid1)%force = patch_ib(pid1)%forces + normal_force + tangental_force
                    patch_ib(pid1)%torque = patch_ib(pid1)%torque + torque    

                    ! apoply equal and opposite force/torque to second sphere
                    atch_ib(pid1)%force = patch_ib(pid1)%forces - normal_force - tangental_force
                    patch_ib(pid1)%torque = patch_ib(pid1)%torque - torque   
                end if
            end if
        end do

    end subroutine s_apply_ib_collision_forces_soft_sphere

    !> @brief applyies collision forces to IBs assuming a soft-sphere collision model (all IBs are circles or spheres)
    subroutine s_apply_wall_collision_forces_soft_sphere(wall_overlap_distances)

        real(wp), dimension(num_ibs, (num_dims*2)), intent(in) :: wall_overlap_distances

        integer :: patch_id, i
        real(wp) :: spring_stiffness, damping_parameter, coefficient_of_friction
        real(wp), dimension(3) :: normal_force, tangental_force, normal_vector, normal_velocity, tangental_vector, collision_location, torque
        real(wp) :: k, eta ! the spring stiffness and damping coefficient for a specific IB

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

                ! ensure the local rank owns that collision before proceeding
                collision_location = [patch_ib(patch_id)%x_centroid, patch_ib(patch_id)%y_centroid, 0._wp]
                if (num_dims == 3) collision_location(3) = patch_ib(patch_id)%z_centroid
                collision_location = collision_location + normal_vector * patch_ib(patch_id)%radius
                if (f_local_rank_owns_collision(collision_location)) then

                    k = spring_stiffness * patch_ib(patch_id)%mass
                    eta = damping_parameter * sqrt(patch_ib(patch_id)%mass * k)

                    ! standard soft-sphere collision  with the wall
                    normal_velocity = dot_product(patch_ib(patch_id)%vel, normal_vector)*normal_vector
                    tangental_vector = patch_ib(patch_id)%vel - normal_velocity
                    tangental_vector = tangental_vector / norm2(tangental_vector)
                    normal_force = -k * wall_overlap_distances(patch_id, i)  * normal_vector - eta * normal_velocity
                    tangental_force = -coefficient_of_friction * norm2(normal_force) * tangental_vector
                    call s_cross_product(normal_vector * patch_ib(patch_id)%radius, tangental_force, torque)
                    patch_ib(patch_id)%force = patch_ib(patch_id)%forces + normal_force + tangental_force
                    patch_ib(patch_id)%torque = patch_ib(patch_id)%torque + torque

                end if
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

    !> @brief function checks if this local MPI processor owns this specific collision
    function f_local_rank_owns_collision(collision_location) result(owns_collision)

      real(wp), dimension(3), intent(in) :: collision_location
      logical, intent(out) :: owns_collision

      #:if defined('MFC_MPI')
          if (num_procs == 1) then
              owns_collision = .true.
          else
              ! catch the edge case where th collision lies just outside the computational domain
              #:for X, ID in [('x', 1), ('y', 2), ('z', 3)]
                  if (num_dims >= ${ID}$) then
                      if (bc_${X}$%beg /= BC_PERIODIC) then
                          if (collision_location(${ID}$) < ${X}$_domain%beg &
                            .or. ${X}$_domain%end < collision_location(${ID}$)) then
                                owns_collision = .true.
                                return
                            end if
                      end if
                  end if
              #:endfor

              ! the object that contains the collision lcoation owns the collisions
              owns_collision = x_cb(-1) <= owns_collision(1) .and. owns_collision(1) < x_cb(m)
              owns_collision = owns_collision .and. y_cb(-1) <= owns_collision(2) .and. owns_collision(2) < y_cb(n)
              if (p /= 0) owns_collision = owns_collision .and. z_cb(-1) <= owns_collision(3) .and. owns_collision(3) < z_cb(p)
          end if
      #:else
          owns_collision = .true.
      #:endif

    end function f_local_rank_owns_collision

end module m_collisions
