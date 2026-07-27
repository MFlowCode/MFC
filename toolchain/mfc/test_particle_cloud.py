import math

import pytest

from mfc.case_validator import CaseConstraintError, CaseValidator

MASK32 = 0xFFFFFFFF
INT32_SIGN = 0x80000000
INT32_MOD = 0x100000000
INT32_MAX = 2147483647.0


def _int32(value):
    value &= MASK32
    if value & INT32_SIGN:
        value -= INT32_MOD
    return value


def _xorshift(seed):
    seed = _int32(seed ^ _int32(seed << 13))
    seed = _int32(seed ^ _int32((seed & MASK32) >> 17))
    seed = _int32(seed ^ _int32(seed << 5))
    return seed, abs(float(seed)) / INT32_MAX


def _hash_cell(point, cell_size):
    return tuple(math.floor(coord / cell_size) for coord in point)


def _is_overlapping(candidate, placed, cells, cell_size, min_distance):
    cell = _hash_cell(candidate, cell_size)
    min_distance_squared = min_distance * min_distance

    for dz in range(-1, 2):
        for dy in range(-1, 2):
            for dx in range(-1, 2):
                neighbor = (cell[0] + dx, cell[1] + dy, cell[2] + dz)
                for particle_idx in cells.get(neighbor, []):
                    distance_squared = sum((candidate[axis] - placed[particle_idx][axis]) ** 2 for axis in range(3))
                    if distance_squared < min_distance_squared:
                        return True

    return False


def _pack_hemi_shell_3d(count, radius, min_spacing, inner_radius, outer_radius, seed_value):
    centroid = (0.0, 0.0, 0.0)
    center_inner = inner_radius + radius
    center_outer = outer_radius - radius
    min_center_distance = 2.0 * radius + min_spacing
    cell_size = min_center_distance
    cells = {}
    placed = []
    seed = seed_value
    attempts = 0

    while len(placed) < count and attempts < count * 1000:
        attempts += 1
        seed, random_phi = _xorshift(seed)
        seed, random_zdir = _xorshift(seed)
        seed, random_radius = _xorshift(seed)

        phi = 2.0 * math.pi * random_phi
        zdir = random_zdir
        rho = math.sqrt(max(0.0, 1.0 - zdir * zdir))
        shell_radius = ((center_outer**3 - center_inner**3) * random_radius + center_inner**3) ** (1.0 / 3.0)
        candidate = (
            centroid[0] + shell_radius * rho * math.cos(phi),
            centroid[1] + shell_radius * rho * math.sin(phi),
            centroid[2] + shell_radius * zdir,
        )

        if candidate[2] < centroid[2] + radius:
            continue
        if _is_overlapping(candidate, placed, cells, cell_size, min_center_distance):
            continue

        cells.setdefault(_hash_cell(candidate, cell_size), []).append(len(placed))
        placed.append(candidate)

    if len(placed) != count:
        raise RuntimeError(f"Only placed {len(placed)} particles after {attempts} attempts")

    return placed


def _valid_cloud_params():
    return {
        "m": 31,
        "n": 31,
        "p": 31,
        "model_eqns": 2,
        "num_fluids": 1,
        "num_patches": 1,
        "t_step_start": 0,
        "t_step_stop": 1,
        "t_step_save": 1,
        "dt": 1.0e-6,
        "x_domain%beg": -1.0,
        "x_domain%end": 1.0,
        "y_domain%beg": -1.0,
        "y_domain%end": 1.0,
        "z_domain%beg": -1.0,
        "z_domain%end": 1.0,
        "bc_x%beg": -1,
        "bc_x%end": -1,
        "bc_y%beg": -1,
        "bc_y%end": -1,
        "bc_z%beg": -1,
        "bc_z%end": -1,
        "riemann_solver": 2,
        "wave_speeds": 1,
        "avg_state": 2,
        "fd_order": 2,
        "ib": "T",
        "num_ibs": 0,
        "num_particle_clouds": 1,
        "particle_cloud(1)%geometry": 2,
        "particle_cloud(1)%packing_method": 1,
        "particle_cloud(1)%num_particles": 8,
        "particle_cloud(1)%radius": 0.04,
        "particle_cloud(1)%mass": 1.0,
        "particle_cloud(1)%min_spacing": 0.01,
        "particle_cloud(1)%shell_inner_radius": 0.1,
        "particle_cloud(1)%shell_outer_radius": 0.3,
    }


def test_hemi_shell_random_packing_respects_clearance_and_overlap():
    radius = 0.04
    min_spacing = 0.01
    inner_radius = 0.1
    outer_radius = 0.4
    particles = _pack_hemi_shell_3d(
        count=16,
        radius=radius,
        min_spacing=min_spacing,
        inner_radius=inner_radius,
        outer_radius=outer_radius,
        seed_value=12345,
    )

    min_center_distance = 2.0 * radius + min_spacing
    for particle in particles:
        radial_distance = math.sqrt(sum(coord * coord for coord in particle))
        assert radial_distance >= inner_radius + radius - 1.0e-12
        assert radial_distance <= outer_radius - radius + 1.0e-12
        assert particle[2] >= radius - 1.0e-12

    for i, particle_i in enumerate(particles):
        for particle_j in particles[i + 1 :]:
            distance = math.sqrt(sum((particle_i[axis] - particle_j[axis]) ** 2 for axis in range(3)))
            assert distance >= min_center_distance - 1.0e-12


def test_hemi_shell_validator_accepts_feasible_shell():
    CaseValidator(_valid_cloud_params()).validate("simulation")


def test_hemi_shell_validator_rejects_negative_inner_radius():
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_inner_radius": -0.01}
    with pytest.raises(CaseConstraintError, match="shell_inner_radius >= 0"):
        CaseValidator(params).validate("simulation")


def test_hemi_shell_validator_rejects_shell_thinner_than_particle_diameter():
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_outer_radius": 0.17}
    with pytest.raises(CaseConstraintError, match="shell_outer_radius > shell_inner_radius"):
        CaseValidator(params).validate("simulation")
