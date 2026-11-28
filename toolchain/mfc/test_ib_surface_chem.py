"""Compile the immersed-boundary reacting-surface solver (m_ib_surface_chem) and check it against Cantera, without a flow solve."""

import subprocess

import cantera as ct
import numpy as np
import pytest

from mfc.params.generators.fortran_gen import generate_constants_fpp
from mfc.test_surface_chemistry_codegen import GAS, cantera_reference, carbon_surface
from mfc.test_thermochem import ROOT, compile_kernel, fypp
from mfc.thermochem import generate_fortran, generate_surface_fortran

DRIVER = """
program reference
    use m_thermochem
    use m_ib_surface_chem
    implicit none
    integer :: ierr, mode, bc
    real(KIND) :: pres, t_ip, t_wall, d, t_s, mdot, t_lo, t_hi, theta_max, t_g
    real(KIND) :: y_ip(num_species), y_s(num_species), y_g(num_species)
    logical :: converged
    do
        read(*,*,iostat=ierr) mode
        if (ierr /= 0) exit
        if (mode == 1) then
            read(*,*) pres, t_ip, t_wall, d, bc, y_ip
            call s_solve_surface(pres, t_ip, t_wall, d, y_ip, molecular_weights, bc, y_s, t_s, mdot, converged)
            write(*,'(I2,1X,*(ES25.16E3,1X))') merge(1, 0, converged), t_s, mdot, y_s
        else
            read(*,*) t_ip, t_s, t_lo, t_hi, theta_max, y_ip, y_s
            call s_blend_ghost_temperature(t_ip, t_s, t_lo, t_hi, theta_max, t_g)
            call s_blend_ghost_species(y_ip, y_s, y_g)
            write(*,'(*(ES25.16E3,1X))') t_g, y_g
        end if
    end do
end program
"""

PRES, D = ct.one_atm, 1.0e-4  # 1 atm, a 0.1 mm image-point distance


def run(tmp_path, gas, surface, lines):
    """Compile m_thermochem, m_surface_thermochem, m_constants and m_ib_surface_chem with DRIVER and feed it lines."""
    (tmp_path / "generated_constants.fpp").write_text(generate_constants_fpp())
    parts = [
        fypp(tmp_path, "m_surface_thermochem", generate_surface_fortran(gas, surface)),
        fypp(tmp_path, "m_constants", (ROOT / "src/common/m_constants.fpp").read_text()),
        fypp(tmp_path, "m_ib_surface_chem", (ROOT / "src/simulation/m_ib_surface_chem.fpp").read_text()),
    ]
    driver = "".join(p.read_text() for p in parts) + DRIVER
    executable = compile_kernel(tmp_path, gas, source=generate_fortran(gas), driver_source=driver)
    out = subprocess.run([str(executable)], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True).stdout
    return [np.fromstring(line, sep=" ") for line in out.splitlines()]


def solve_line(t_ip, t_wall, bc, y_ip):
    return f"1\n{PRES} {t_ip} {t_wall} {D} {bc} " + " ".join(map(str, y_ip))


def gas_state(gas, t, x):
    gas.TPX = t, PRES, x
    return gas.Y.copy()


def flux_residual(gas, surface, t_ip, y_ip, t_s, y_s):
    """The species flux balance of s_surface_species_residual and the energy balance, evaluated with Cantera, and their scales."""
    w = gas.molecular_weights
    gas.TPY = t_s, PRES, y_s
    rho, mw, diff, k_s = gas.density, gas.mean_molecular_weight, gas.mix_diff_coeffs, gas.thermal_conductivity
    x_s = gas.X.copy()
    gas.TPY = t_ip, PRES, y_ip
    rate = cantera_reference(surface, t_s, PRES, y_s)
    omega, q = rate[:-1], rate[-1]
    bg = rho * diff * w / mw * (gas.X - x_s) / D
    mdot = w @ omega
    species = -bg + y_s * (bg.sum() + mdot) - w * omega
    return species, np.max(np.abs(w * omega)), k_s * (t_s - t_ip) / D - q, abs(q), mdot


@pytest.fixture(scope="module")
def carbon():
    return ct.Solution(GAS), carbon_surface()


@pytest.mark.parametrize("bc,t_wall", [(1, 1300.0), (2, 0.0)])
def test_newton_solve_satisfies_the_surface_balance(tmp_path, carbon, bc, t_wall):
    """Prescribed wall temperature (1) and the surface energy balance (2) on burning carbon in moist air at 1200 K: the returned
    state closes the flux balance computed by Cantera, keeps sum(Y) = 1 and its gasification rate is Cantera's."""
    gas, surface = carbon
    y_ip = gas_state(gas, 1200.0, "O2:1,N2:3.76,H2O:0.3")
    (row,) = run(tmp_path, gas, surface, [solve_line(1200.0, t_wall, bc, y_ip)])
    assert row[0] == 1, "the Newton solve did not converge"
    t_s, mdot, y_s = row[1], row[2], row[3:]
    if bc == 1:
        assert t_s == t_wall
    assert np.all(y_s >= 0) and abs(y_s.sum() - 1) < 1e-6
    species, scale, energy, q_scale, mdot_ct = flux_residual(gas, surface, 1200.0, y_ip, t_s, y_s)
    assert np.max(np.abs(species)) < 1e-5 * scale
    assert mdot > 0 and mdot == pytest.approx(mdot_ct, rel=1e-10)
    if bc == 2:
        assert t_s > 1210.0  # the reaction heat lifts the surface above the gas (1223 K)
        assert abs(energy) < 1e-5 * max(q_scale, 1.0)


def test_inert_surface_returns_the_gas_state(tmp_path):
    """Without a surface mechanism every rate is zero: the surface state is the image-point gas, solved at once, no gasification."""
    gas = ct.Solution(GAS)
    y_ip = gas_state(gas, 900.0, "O2:1,N2:3.76")
    (row,) = run(tmp_path, gas, None, [solve_line(900.0, 0.0, 0, y_ip)])
    assert row[0] == 1 and row[1] == 900.0 and row[2] == 0.0
    np.testing.assert_array_equal(row[3:], y_ip)


def test_unknown_thermal_bc_is_not_converged(tmp_path, carbon):
    """A thermal_bc the solve does not know leaves the surface inert and reports no convergence, which the caller counts."""
    gas, surface = carbon
    y_ip = gas_state(gas, 1500.0, "O2:1,N2:3.76")
    (row,) = run(tmp_path, gas, surface, [solve_line(1500.0, 0.0, 7, y_ip)])
    assert row[0] == 0 and row[1] == 1500.0 and row[2] == 0.0


def test_ghost_blend_stays_physical(tmp_path, carbon):
    """For random image-point and surface states the ghost temperature stays in [T_lo, T_hi] and the ghost species stay
    non-negative with the surface's sum."""
    gas, surface = carbon
    rng = np.random.default_rng(11)
    cases, lines = [], []
    for _ in range(40):
        y_ip, y_s = (rng.dirichlet(np.full(gas.n_species, 0.3)) for _ in range(2))
        t_ip, t_s = rng.uniform(250.0, 3000.0, 2)
        t_lo, t_hi = 200.0, 3500.0
        cases.append((t_lo, t_hi, y_s))
        lines.append(f"2\n{t_ip} {t_s} {t_lo} {t_hi} {rng.uniform(0, 1)} " + " ".join(map(str, [*y_ip, *y_s])))
    for (t_lo, t_hi, y_s), row in zip(cases, run(tmp_path, gas, surface, lines)):
        t_g, y_g = row[0], row[1:]
        assert t_lo <= t_g <= t_hi
        assert np.all(y_g >= 0)
        assert y_g.sum() == pytest.approx(y_s.sum(), abs=1e-14)
