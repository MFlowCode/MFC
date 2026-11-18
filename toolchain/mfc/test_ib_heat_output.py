"""Parsing and heat bookkeeping of D/ib_heat.dat (ib_heat_wrt), through the lumped-wall example's reader."""

import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("heat_balance", ROOT / "examples/2D_ibm_lumped_wall/heat_balance.py")
heat_balance = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(heat_balance)

REC_LEN = 10 + 5 * 18 + 1  # I10, five 1X,ES17.9E3 fields, newline: docs/documentation/case.md


def _es(x):
    """x as Fortran ES17.9E3 writes it: three-digit exponent, right-justified in 17 characters."""
    mantissa, exponent = f"{x:.9E}".split("E")
    return f"{mantissa}E{exponent[0]}{abs(int(exponent)):03d}".rjust(17)


def _write(path, bodies):
    """bodies[k] is the list of (time, mdot, mdot_face, qdot, Twall) rows of body k + 1; rows interleave bodies as the writer does."""
    lines = []
    for row in zip(*bodies):
        for gid, fields in enumerate(row, start=1):
            lines.append(f"{gid:10d}" + "".join(" " + _es(x) for x in fields) + "\n")
    assert all(len(ln) == REC_LEN for ln in lines)
    Path(path).write_text("".join(lines))


def test_record_width_matches_the_documented_layout():
    assert len(_es(1.0)) == 17
    assert len(_es(-1.234e-300)) == 17
    assert len(_es(6.02e123)) == 17


def test_reader_orders_rows_by_step_and_bodies_by_id(tmp_path):
    a = [(1.0e-6 * k, 0.0, 0.5, -2.0 * k, 900.0 - k) for k in range(1, 4)]
    b = [(1.0e-6 * k, 1.0e-3, 0.5, 3.0 * k, 300.0 + k) for k in range(1, 4)]
    _write(tmp_path / "ib_heat.dat", [a, b])
    rows = heat_balance.read_ib_heat(tmp_path / "ib_heat.dat", num_ibs=2)
    assert rows.shape == (3, 2, 6)
    np.testing.assert_allclose(rows[:, 0, 0], 1)
    np.testing.assert_allclose(rows[:, 1, 0], 2)
    np.testing.assert_allclose(rows[:, 0, 5], [899.0, 898.0, 897.0])
    np.testing.assert_allclose(rows[:, 1, 4], [3.0, 6.0, 9.0])


def test_a_single_row_file_still_has_three_axes(tmp_path):
    _write(tmp_path / "ib_heat.dat", [[(1.0e-6, 0.0, 0.0, 1.0, 900.0)]])
    assert heat_balance.read_ib_heat(tmp_path / "ib_heat.dat").shape == (1, 1, 6)


def test_balance_residual_vanishes_for_the_explicit_update():
    """T_k = T_{k-1} + Q_k dt_k / C, with the step dt the difference of the row times, is what s_update_ib_temperatures does."""
    heat_capacity = 3.0e-3
    time = np.cumsum(np.r_[0.0, 1.0e-7 * (1 + np.arange(20) % 3)])
    qdot = -np.linspace(5.0, 1.0, time.size)
    twall = np.empty_like(time)
    twall[0] = 900.0
    for k in range(1, time.size):
        twall[k] = twall[k - 1] + qdot[k] * (time[k] - time[k - 1]) / heat_capacity
    np.testing.assert_allclose(heat_balance.balance_residual(time, qdot, twall, heat_capacity), 0.0, atol=1e-12)


def test_balance_residual_sees_a_missing_heat_term():
    heat_capacity = 3.0e-3
    time = np.arange(5) * 1.0e-7
    qdot = np.full(5, -2.0)
    twall = 900.0 + np.arange(5) * (-2.0 * 1.0e-7 / heat_capacity)
    twall[3:] += 0.5
    residual = heat_balance.balance_residual(time, qdot, twall, heat_capacity)
    assert abs(residual[2]) > 1e-4 and np.allclose(np.delete(residual, 2), 0.0, atol=1e-12)
