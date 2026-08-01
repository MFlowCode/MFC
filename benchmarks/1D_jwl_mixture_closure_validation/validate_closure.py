#!/usr/bin/env python3
"""Mixture-closure validation for the JWL weighted-composition model.

Two studies, both scored against a converged pressure-temperature (PT)
equilibrium reference computed independently of the closure:

  lhs         State-space study on the Latin-hypercube protocol of Jackson's
              JWL EOS notes: N states with p log-uniform in [1e4, 1e8] Pa,
              T uniform in [300, 5000] K, Y uniform in [1e-6, 1 - 1e-6],
              seed 12345, constructed ON the PT-equilibrium manifold. The
              shipped weighted-composition closure and the legacy Rocflu
              density/energy-ramp closure it replaced (RFLU_ModJWL behavior,
              Garno et al., Phys. Rev. Fluids 5, 123201, 2020) are evaluated
              at the same (rho, e, Y) and their p/T/c errors tabulated.

  simulation  Reads the 1D mixing-shocktube output of case.py in this
              directory and scores the actual mixture cells the solver
              produced: (a) transcription fidelity, this script's closure
              pressure against MFC's own output pressure, and (b) closure
              accuracy against the PT-equilibrium reference over the mixed
              band 0 < Y < 1.

The PT-equilibrium reference solves the exact scalar reduction of the
two-phase equilibrium system: with both constituents Grueneisen-caloric
(e_k = e_cold,k(rho_k) + cv_k T), temperature equality plus energy
conservation are linear in T, volume additivity is explicit, and only
pressure equality is nonlinear:

    T(rho_p)     = (e - Y e_cold(rho_p)) / (Y cv_j + (1-Y) cv_a)
    rho_a(rho_p) = (1-Y) / (1/rho - Y/rho_p)
    R(rho_p)     = p_JWL(rho_p, T) - p_air(rho_a, T) = 0

solved by a bracketed scan plus bisection and Newton polish, with the
physicality guard T > 0. The equilibrium sound speed follows from the
implicit function theorem. All formulas are transcribed from the validated
standalone study backing this PR's closure (jwl_standalone, 2026).

Exit code 0 iff every PASS gate holds.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

# TNT products + air, matching case.py and toolchain/mfc/jwl_products.py.
A = 3.712e11  # Pa
B = 3.231e9  # Pa
R1 = 4.15
R2 = 0.95
OMEGA0 = 0.30
RHO0 = 1630.0  # kg/m^3
E0 = 1.0089e10  # J/m^3 (legacy ramp reference only; the shipped closure has no E0 dependence)
CV_J = 613.5  # J/(kg K)
AIR_GAMMA = 0.4  # Grueneisen Gamma = gamma_phys - 1 = 1/fluid_pp%gamma
AIR_RHO0 = 1.225  # kg/m^3 (legacy ramp reference only)
AIR_E0 = 2.5575e5  # J/kg   (legacy ramp reference only)
CV_A = 717.5  # J/(kg K)
SGM_EPS = 1.0e-16

# Jackson LHS protocol.
SEED = 12345
P_MIN, P_MAX = 1.0e4, 1.0e8
T_MIN, T_MAX = 300.0, 5000.0
Y_MIN, Y_MAX = 1.0e-6, 1.0 - 1.0e-6


# --- shipped closure: transcription of s_jwl_weighted_composition_state_er ---


def cb_coeffs(Y):
    w = Y * CV_J / (Y * CV_J + (1.0 - Y) * CV_A)
    An = w * A
    Bn = w * B
    omega = AIR_GAMMA + w * (OMEGA0 - AIR_GAMMA)
    cv = Y * CV_J + (1.0 - Y) * CV_A
    return An, Bn, omega, cv


def cb_state(rho, e, Y):
    """Shipped weighted-composition closure: (rho, e, Y) -> (p, T, c2, c2_floor)."""
    rho_s = np.maximum(rho, SGM_EPS)
    Y_s = np.clip(Y, 0.0, 1.0)
    An, Bn, omega, cv = cb_coeffs(Y_s)
    V = RHO0 / rho_s
    exp1 = np.exp(-R1 * V)
    exp2 = np.exp(-R2 * V)
    coef1 = (1.0 - omega / (R1 * V)) * exp1
    coef2 = (1.0 - omega / (R2 * V)) * exp2
    p = An * coef1 + Bn * coef2 + omega * rho_s * e
    T = (p - An * exp1 - Bn * exp2) / (omega * cv * rho_s)
    c2 = exp1 * An * (R1 * RHO0 / rho_s**2 - omega / rho_s - omega / (R1 * RHO0)) + exp2 * Bn * (R2 * RHO0 / rho_s**2 - omega / rho_s - omega / (R2 * RHO0)) + omega * (e + p / rho_s)
    p = np.maximum(p, SGM_EPS)
    T = np.maximum(T, SGM_EPS)
    c2_floor = AIR_GAMMA * np.maximum(p, SGM_EPS) / rho_s
    return p, T, c2, c2_floor


def cb_energy_pr(rho, p, Y):
    """Shipped analytic inverse: (rho, p, Y) -> e (s_jwl_weighted_composition_energy_pr)."""
    rho_s = np.maximum(rho, SGM_EPS)
    An, Bn, omega, _ = cb_coeffs(np.clip(Y, 0.0, 1.0))
    V = RHO0 / rho_s
    C1 = (1.0 - omega / (R1 * V)) * np.exp(-R1 * V)
    C2 = (1.0 - omega / (R2 * V)) * np.exp(-R2 * V)
    return np.maximum((p - An * C1 - Bn * C2) / np.maximum(omega * rho_s, SGM_EPS), 0.0)


# --- legacy comparator: Rocflu density/energy-ramp closure (pre-PR MFC heritage) ---


def jwl_core(rho, e, An, Bn, omega, mA, mB, momega, cv):
    """Blended JWL kernel with coefficient derivatives in c2 (legacy closures only)."""
    rho_s = np.maximum(rho, SGM_EPS)
    V = RHO0 / rho_s
    exp1 = np.exp(-R1 * V)
    exp2 = np.exp(-R2 * V)
    coef1 = (1.0 - omega / (R1 * V)) * exp1
    coef2 = (1.0 - omega / (R2 * V)) * exp2
    p = An * coef1 + Bn * coef2 + omega * rho_s * e
    T = (p - An * exp1 - Bn * exp2) / (omega * cv * rho_s)
    c2 = (
        exp1 * An * (R1 * RHO0 / rho_s**2 - omega / rho_s - omega / (R1 * RHO0) - rho_s * momega / (R1 * RHO0))
        + mA * p * coef1 / rho_s**2
        + exp2 * Bn * (R2 * RHO0 / rho_s**2 - omega / rho_s - omega / (R2 * RHO0) - rho_s * momega / (R2 * RHO0))
        + mB * p * coef2 / rho_s**2
        + omega * (e + p / rho_s)
        + momega * rho_s * e
    )
    return p, T, c2


def legacy_rocflu_state(rho, e, Y):
    """Legacy Rocflu closure: hard Y switches, (rho, e)-ramped coefficients."""
    rho = np.asarray(rho, dtype=float)
    e = np.asarray(e, dtype=float)
    Y = np.asarray(Y, dtype=float)

    # Interior blend: A/B ramp linearly in e between air_e0 and e_j; omega and cv
    # ramp linearly in rho between air_rho0 and rho0 (Y plays no role, matching
    # RFLU_ModJWL.F90).
    ej = E0 / RHO0
    phi = np.clip((e - AIR_E0) / (ej - AIR_E0), 0.0, 1.0)
    ramp_on = (e > AIR_E0) & (e < ej)
    mA = np.where(ramp_on, A / (ej - AIR_E0), 0.0)
    mB = np.where(ramp_on, B / (ej - AIR_E0), 0.0)
    slope_w = (OMEGA0 - AIR_GAMMA) / (RHO0 - AIR_RHO0)
    omega = np.clip(AIR_GAMMA + slope_w * (rho - AIR_RHO0), AIR_GAMMA, OMEGA0)
    momega = np.where((rho >= AIR_RHO0) & (rho < RHO0), slope_w, 0.0)
    slope_cv = (CV_J - CV_A) / (RHO0 - AIR_RHO0)
    cv = np.clip(CV_A + slope_cv * (rho - AIR_RHO0), min(CV_A, CV_J), max(CV_A, CV_J))
    p, T, c2 = jwl_core(rho, e, A * phi, B * phi, omega, mA, mB, momega, cv)

    # Hard endpoint overrides.
    air = Y <= 0.01
    p_air = AIR_GAMMA * rho * e
    p = np.where(air, p_air, p)
    T = np.where(air, e / CV_A, T)
    c2 = np.where(air, (AIR_GAMMA + 1.0) * p_air / rho, c2)
    prod = Y > 0.99
    if np.any(prod):
        pp, Tp, c2p = jwl_core(rho, e, A, B, OMEGA0, 0.0, 0.0, 0.0, CV_J)
        p = np.where(prod, pp, p)
        T = np.where(prod, Tp, T)
        c2 = np.where(prod, c2p, c2)
    return p, T, c2


# --- PT-equilibrium reference ---


def e_cold(rho_p):
    """JWL cold-curve energy, chosen so T = (e - e_cold)/cv matches the pure JWL law."""
    V = RHO0 / np.maximum(rho_p, SGM_EPS)
    return (A / R1 * np.exp(-R1 * V) + B / R2 * np.exp(-R2 * V)) / RHO0


def products_p(rho_p, T):
    V = RHO0 / np.maximum(rho_p, SGM_EPS)
    return A * np.exp(-R1 * V) + B * np.exp(-R2 * V) + OMEGA0 * rho_p * CV_J * T


def scalar_residual(rho, e, Y, rho_p):
    """R(rho_p) = p_products - p_air and its analytic derivative, plus phase states."""
    cvm = Y * CV_J + (1.0 - Y) * CV_A
    u = 1.0 / rho - Y / rho_p  # (1-Y)/rho_a
    ok = (u > SGM_EPS) & (rho_p > 0.0)
    u = np.where(ok, u, 1.0)
    rho_a = (1.0 - Y) / u
    drho_a = -(Y / rho_p**2) * rho_a**2 / (1.0 - Y)

    V = RHO0 / rho_p
    e1 = np.exp(-R1 * V)
    e2 = np.exp(-R2 * V)
    pref = A * e1 + B * e2
    dpref = (A * R1 * e1 + B * R2 * e2) * RHO0 / rho_p**2

    T = (e - Y * e_cold(rho_p)) / cvm
    dT = -(Y / cvm) * pref / rho_p**2  # d e_cold/d rho_p = pref/rho_p^2 (exact JWL identity)

    p_p = pref + OMEGA0 * rho_p * CV_J * T
    p_a = AIR_GAMMA * rho_a * CV_A * T
    R = np.where(ok, p_p - p_a, np.inf)
    dR = dpref + OMEGA0 * CV_J * (T + rho_p * dT) - AIR_GAMMA * CV_A * (drho_a * T + rho_a * dT)
    return R, dR, T, p_p, p_a, rho_a, ok


def pt_solve(rho, e, Y, nscan=200, nbis=80, npolish=4):
    """Converged PT equilibrium: (rho, e, Y) -> (p, T, rho_p, ok), vectorized.

    Bracketed scan over ln(rho_p - Y rho), bisection, Newton polish; roots with
    T <= 0 are rejected (spurious deep-tension family).
    """
    rho = np.asarray(rho, dtype=float)
    n = rho.size
    lo = np.log(Y * rho * 1.0e-12)  # rho_p - Y*rho spans many decades
    hi = np.log(np.maximum(10.0 * RHO0 - Y * rho, Y * rho))
    xs = lo[None, :] + (hi - lo)[None, :] * np.linspace(0.0, 1.0, nscan)[:, None]

    Rs = np.empty((nscan, n))
    Ts = np.empty((nscan, n))
    for k in range(nscan):  # keep memory bounded; each row is vectorized
        Rk, _, Tk, _, _, _, okk = scalar_residual(rho, e, Y, Y * rho + np.exp(xs[k]))
        Rs[k] = np.where(okk, Rk, np.inf)
        Ts[k] = Tk

    # First sign change with a physical (T > 0) left endpoint.
    valid = np.isfinite(Rs[:-1]) & np.isfinite(Rs[1:]) & (Rs[:-1] * Rs[1:] <= 0.0) & (Ts[:-1] > 0.0)
    has = valid.any(axis=0)
    kfirst = np.where(has, valid.argmax(axis=0), 0)
    idx = np.arange(n)
    xlo = xs[kfirst, idx]
    xhi = xs[kfirst + 1, idx]
    Rlo, _, _, _, _, _, _ = scalar_residual(rho, e, Y, Y * rho + np.exp(xlo))

    for _ in range(nbis):
        xm = 0.5 * (xlo + xhi)
        Rm, _, _, _, _, _, _ = scalar_residual(rho, e, Y, Y * rho + np.exp(xm))
        left = Rlo * Rm <= 0.0
        xhi = np.where(left, xm, xhi)
        xlo = np.where(left, xlo, xm)
        Rlo = np.where(left, Rlo, Rm)

    x = 0.5 * (xlo + xhi)
    for _ in range(npolish):
        Rm, dRm, _, _, _, _, okm = scalar_residual(rho, e, Y, Y * rho + np.exp(x))
        dx = Rm / np.where(np.abs(dRm * np.exp(x)) > 0.0, dRm * np.exp(x), 1.0)
        xn = x - dx
        good = okm & (xn > xlo) & (xn < xhi)
        x = np.where(good, xn, x)

    rho_p = Y * rho + np.exp(x)
    R, _, T, p_p, p_a, _, okr = scalar_residual(rho, e, Y, rho_p)
    p = (Y * rho / rho_p) * p_p + (1.0 - Y * rho / rho_p) * p_a
    scale = np.maximum(np.maximum(np.abs(p_p), np.abs(p_a)), 1.0e5)
    ok = has & okr & (T > 0.0) & (np.abs(R) < 1.0e-6 * scale)
    return p, T, rho_p, ok


def equilibrium_c2(rho, e, Y, rho_p):
    """Equilibrium sound speed at the converged root, by the implicit function theorem."""
    R, F_rp, T, p_p, p_a, rho_a, ok = scalar_residual(rho, e, Y, rho_p)
    cvm = Y * CV_J + (1.0 - Y) * CV_A

    V = RHO0 / rho_p
    e1 = np.exp(-R1 * V)
    e2 = np.exp(-R2 * V)
    pref = A * e1 + B * e2
    dpref = (A * R1 * e1 + B * R2 * e2) * RHO0 / rho_p**2

    dT_rp = -(Y / cvm) * pref / rho_p**2
    dT_e = 1.0 / cvm
    dra_rp = -(Y / rho_p**2) * rho_a**2 / (1.0 - Y)
    dra_rho = rho_a**2 / ((1.0 - Y) * rho**2)

    dpp_rp = dpref + OMEGA0 * CV_J * T  # partial in rho_p at fixed T
    dpp_T = OMEGA0 * rho_p * CV_J

    F_rho = -AIR_GAMMA * CV_A * dra_rho * T
    F_e = dpp_T * dT_e - AIR_GAMMA * CV_A * rho_a * dT_e

    F_rp = np.where(np.abs(F_rp) > 1.0e-30, F_rp, np.nan)
    drp_drho = -F_rho / F_rp
    drp_de = -F_e / F_rp

    dp_drho = dpp_rp * drp_drho + dpp_T * dT_rp * drp_drho
    dp_de = dpp_rp * drp_de + dpp_T * (dT_rp * drp_de + dT_e)

    pmix = (Y * rho / rho_p) * p_p + (1.0 - Y * rho / rho_p) * p_a
    c2 = dp_drho + (pmix / rho**2) * dp_de
    return c2, ok & np.isfinite(c2)


# --- LHS state-space study ---


def lhs_manifold_states(nsamp, seed=SEED):
    """Draw (p, T, Y) by Latin hypercube and construct (rho, e, Y) on the PT manifold."""
    rng = np.random.default_rng(seed)

    def stratified():
        return (rng.permutation(nsamp) + rng.random(nsamp)) / nsamp

    p_ex = 10.0 ** (np.log10(P_MIN) + stratified() * (np.log10(P_MAX) - np.log10(P_MIN)))
    T_ex = T_MIN + stratified() * (T_MAX - T_MIN)
    Y = Y_MIN + stratified() * (Y_MAX - Y_MIN)

    # Products density from p(rho_p, T) = p_ex: log-space bisection (monotone branch).
    xlo = np.full(nsamp, np.log(1.0e-4))
    xhi = np.full(nsamp, np.log(5.0 * RHO0))
    flo = products_p(np.exp(xlo), T_ex) - p_ex
    fhi = products_p(np.exp(xhi), T_ex) - p_ex
    ok = flo * fhi <= 0.0
    for _ in range(100):
        xm = 0.5 * (xlo + xhi)
        fm = products_p(np.exp(xm), T_ex) - p_ex
        left = flo * fm <= 0.0
        xhi = np.where(left, xm, xhi)
        xlo = np.where(left, xlo, xm)
        flo = np.where(left, flo, fm)
    rho_p = np.exp(0.5 * (xlo + xhi))
    ok &= np.abs(products_p(rho_p, T_ex) - p_ex) <= 1.0e-6 * p_ex

    rho_a = p_ex / (AIR_GAMMA * CV_A * T_ex)
    e_p = e_cold(rho_p) + CV_J * T_ex
    e_a = CV_A * T_ex
    rho = 1.0 / (Y / rho_p + (1.0 - Y) / rho_a)
    e = Y * e_p + (1.0 - Y) * e_a
    return p_ex, T_ex, Y, rho, e, ok


def pct_row(v):
    v = np.sort(v)
    n = len(v)
    grab = lambda q: v[max(0, min(n - 1, int(q * n + 0.5) - 1))]
    return grab(0.50), grab(0.95), grab(0.99), v[-1]


def run_lhs(nsamp, csv_path=None):
    print("=" * 74)
    print(" JWL mixture-closure state-space study (Jackson LHS protocol)")
    print(f" N = {nsamp}, seed {SEED}, TNT products + air (CW baseline regime)")
    print(f" p in [1e4, {P_MAX:.0e}] Pa log-uniform, T in [300, 5000] K, Y in (0, 1)")
    print(" P_max is the equilibrium mixing ceiling, not the compressed CJ regime;")
    print(" the compressed tail is exercised by the in-simulation reflected band.")
    print("=" * 74)

    p_ex, T_ex, Y, rho, e, ok = lhs_manifold_states(nsamp)

    p_ref, T_ref, rho_p, ok_ref = pt_solve(rho, e, Y)
    c2_ref, ok_c = equilibrium_c2(rho, e, Y, rho_p)
    # Manifold consistency: the reference must recover the constructed state.
    ok &= ok_ref & ok_c & (c2_ref > 0.0) & (np.abs(p_ref - p_ex) <= 1.0e-4 * np.abs(p_ex))
    n_valid = int(ok.sum())
    print(f" valid states: {n_valid}   censored: {nsamp - n_valid}")

    rho, e, Y = rho[ok], e[ok], Y[ok]
    p_ref, T_ref, c_ref = p_ref[ok], T_ref[ok], np.sqrt(c2_ref[ok])

    models = {}
    p_cb, T_cb, c2_cb, c2f = cb_state(rho, e, Y)
    models["weighted-composition (shipped)"] = (p_cb, T_cb, np.sqrt(np.maximum(c2_cb, c2f)))
    p_lg, T_lg, c2_lg = legacy_rocflu_state(rho, e, Y)
    models["legacy Rocflu ramp (replaced)"] = (p_lg, T_lg, np.sqrt(np.maximum(c2_lg, SGM_EPS)))

    errs = {}
    for name, (p, T, c) in models.items():
        errs[name] = (
            100.0 * np.abs(p - p_ref) / np.abs(p_ref),
            100.0 * np.abs(T - T_ref) / np.abs(T_ref),
            100.0 * np.abs(c - c_ref) / c_ref,
        )

    for qi, qname in enumerate(("p", "T", "c")):
        print("-" * 74)
        print(f" model                              |{qname}| err %:   median      p95      p99      max")
        for name, err in errs.items():
            med, p95, p99, mx = pct_row(err[qi])
            print(f"   {name:34s} {med:12.2e} {p95:8.1e} {p99:8.1e} {mx:8.1e}")
    print("-" * 74)

    # Round trip through the shipped analytic inverse: e -> p -> e.
    e_rt = cb_energy_pr(rho, p_cb, Y)
    rt = np.max(np.abs(e_rt - e) / np.abs(e))
    print(f" shipped inverse round trip max |e' - e|/e: {rt:.2e}")

    if csv_path:
        with open(csv_path, "w") as f:
            f.write("model,quantity,error_pct\n")
            for name, err in errs.items():
                key = "shipped" if "shipped" in name else "legacy"
                for qname, v in zip(("p", "T", "c"), err):
                    for x in v:
                        f.write(f"{key},{qname},{x:.6e}\n")
        print(f" wrote {csv_path}")

    med_p, p95_p, p99_p, max_p = pct_row(errs["weighted-composition (shipped)"][0])
    med_l = pct_row(errs["legacy Rocflu ramp (replaced)"][0])[0]
    gates = [
        ("shipped median |p| err <= 1e-8 %", med_p <= 1.0e-8),
        ("shipped p95    |p| err <= 1 %", p95_p <= 1.0),
        ("shipped p99    |p| err <= 10 %", p99_p <= 10.0),
        ("shipped inverse round trip <= 1e-8", rt <= 1.0e-8),
        ("median improvement over legacy >= 1e6 x", med_l / max(med_p, 1e-300) >= 1.0e6),
        ("valid-state fraction >= 95 %", n_valid >= 0.95 * nsamp),
    ]
    return report_gates(gates)


# --- simulation-output study ---


def read_1d_field(case_dir, kind, var, step):
    """Concatenate D/<kind>.<var>.<proc>.<step>.dat across ranks, sorted by x."""
    pat = os.path.join(case_dir, "D", f"{kind}.{var}.*.{step:06d}.dat")
    files = sorted(glob.glob(pat))
    if not files:
        raise FileNotFoundError(pat)
    data = np.vstack([np.loadtxt(f) for f in files])
    order = np.argsort(data[:, 0])
    return data[order, 0], data[order, 1]


def run_simulation(case_dir):
    steps = sorted({int(m.group(1)) for f in glob.glob(os.path.join(case_dir, "D", "cons.1.*.dat")) for m in [re.search(r"\.(\d{6})\.dat$", f)] if m})
    if not steps:
        print(f" no simulation output under {case_dir}/D; run the case first:")
        print(f"   ./mfc.sh run {case_dir}/case.py -n 2")
        return 2
    step = steps[-1]
    print("=" * 74)
    print(f" JWL mixture-closure in-simulation study: {case_dir}, step {step}")
    print("=" * 74)

    x, ar1 = read_1d_field(case_dir, "cons", 1, step)
    _, ar2 = read_1d_field(case_dir, "cons", 2, step)
    _, mom = read_1d_field(case_dir, "cons", 3, step)
    _, En = read_1d_field(case_dir, "cons", 4, step)
    _, p_mfc = read_1d_field(case_dir, "prim", 4, step)

    rho = ar1 + ar2
    Y = np.clip(ar1 / rho, 0.0, 1.0)
    e = (En - 0.5 * mom**2 / rho) / rho

    # (a) transcription fidelity: this script's closure against the solver's output p.
    p_cb, _, _, _ = cb_state(rho, e, Y)
    fid = np.max(np.abs(p_cb - p_mfc) / np.abs(p_mfc))
    print(f" cells: {len(x)}   closure transcription max |p - p_MFC|/p_MFC: {fid:.2e}")

    # (b) closure accuracy against PT equilibrium over the mixed band.
    band = (Y > 1.0e-9) & (Y < 1.0 - 1.0e-9)
    nb = int(band.sum())
    print(f" mixed-band cells (1e-9 < Y < 1 - 1e-9): {nb}")
    if nb == 0:
        print(" no mixed cells to score; increase t_step_stop or resolution")
        return 2
    rb, eb, Yb, pb = rho[band], e[band], Y[band], p_cb[band]
    p_ref, T_ref, rho_p, ok = pt_solve(rb, eb, Yb)
    nok = int(ok.sum())
    print(f" PT reference converged on {nok}/{nb} mixed cells")
    err = 100.0 * np.abs(pb[ok] - p_ref[ok]) / np.abs(p_ref[ok])
    med, p95, p99, mx = pct_row(err)
    print(f" closure |p| err % vs PT equilibrium:   median {med:.2e}   p95 {p95:.2e}   p99 {p99:.2e}   max {mx:.2e}")
    print("-" * 74)

    gates = [
        ("closure transcription fidelity <= 1e-10", fid <= 1.0e-10),
        ("PT reference converged on >= 95 % of mixed cells", nok >= 0.95 * nb),
        ("mixed-band median |p| err <= 1 %", med <= 1.0),
        ("mixed-band p99    |p| err <= 10 %", p99 <= 10.0),
    ]
    return report_gates(gates)


def report_gates(gates):
    print(" PASS/FAIL gates:")
    ok_all = True
    for name, ok in gates:
        print(f"   [{'PASS' if ok else 'FAIL'}] {name}")
        ok_all &= ok
    print("=" * 74)
    return 0 if ok_all else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--study", choices=("lhs", "simulation", "all"), default="all")
    ap.add_argument("--nsamples", type=int, default=10000, help="LHS sample count")
    ap.add_argument("--case-dir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--csv", default=None, help="optional CSV output path for the LHS study")
    ap.add_argument("--mfc", default=None, help=argparse.SUPPRESS)  # toolchain compatibility
    args = ap.parse_args()

    rc = 0
    if args.study in ("lhs", "all"):
        rc |= run_lhs(args.nsamples, args.csv)
    if args.study in ("simulation", "all"):
        rc |= run_simulation(args.case_dir)
    return rc


if __name__ == "__main__":
    sys.exit(main())
