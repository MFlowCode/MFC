"""Well-known JWL detonation-product parameter sets for MFC case files.

Import these in a case.py instead of hand-entering the JWL coefficient block
(and the transcription errors that invites). Each explosive is a PRODUCTS set: the JWL
coefficients A, B [Pa], R1, R2, omega [-], the reference density rho0 [kg/m^3],
the specific detonation energy Q [J/kg] (MFC derives E0 = rho0*Q), and the
products heat capacity cv [J/(kg K)]. Coefficients are from the LLNL Explosives
Handbook (Dobratz & Crawford, UCRL-52997) unless noted; TNT and LX-10 match
MFC's validated in-repo cases.

An inconsistent set fails fast at startup via MFC's init closure scan, so a
mistaken value aborts rather than silently biasing the physics. The scan checks
consistency and invertibility, not calibration, so the cited reference CJ state
in each comment stays a useful sanity check. The products cv is weakly sensitive
(pure-products pressure and sound speed are cv-free; cv sets only the diagnostic
temperature and the mixture cv blend), and a material-specific value comes from a
thermochemical detonation-products calculation, not from the JWL fit -- so where
none is published a principled placeholder of 1000 J/(kg K) is used (TNT and LX-10
carry their in-repo/paper values). Set a specific one if you have it.

Example (TNT charge in air, two fluids):

    from mfc.jwl_products import TNT, AIR, jwl_fluid, ambient_fluid, znd_delta_e

    params = {
        ...,
        **jwl_fluid(1, TNT, AIR),      # products fluid + its ambient references
        **ambient_fluid(2, AIR),       # the co-existing ideal-gas ambient
    }
    # For a resolved ZND detonation (jwl_reactive), add the reactant offset:
    params.update(jwl_fluid(1, TNT, AIR, delta_e=znd_delta_e(TNT)))
"""

import math

# --- detonation products (JWL parameter sets) ---
# A, B [Pa]; R1, R2, omega [-]; rho0 [kg/m^3]; Q [J/kg]; cv [J/(kg K)].
# The trailing comment gives the source and the reference CJ pressure / speed.

TNT = {  # MFC in-repo TNT (E0 = 1.0089e10 J/m^3); P_CJ ~ 21 GPa, D ~ 6900 m/s
    "name": "TNT",
    "A": 3.712e11,
    "B": 3.231e9,
    "R1": 4.15,
    "R2": 0.95,
    "omega": 0.30,
    "rho0": 1630.0,
    "Q": 1.0089e10 / 1630.0,
    "cv": 613.5,
}
LX10 = {  # Garno et al., J. Appl. Phys. 128, 195903 (2020), Table I; P_CJ ~ 37 GPa, D ~ 8800 m/s
    "name": "LX-10",
    "A": 880.2e9,
    "B": 17.437e9,
    "R1": 4.60,
    "R2": 1.20,
    "omega": 0.30,
    "rho0": 1860.0,
    "Q": 5.59e6,
    "cv": 1000.0,
}
PETN = {  # Dobratz & Crawford, UCRL-52997; P_CJ ~ 33.5 GPa, D ~ 8300 m/s
    "name": "PETN",
    "A": 6.171e11,
    "B": 1.6926e10,
    "R1": 4.40,
    "R2": 1.20,
    "omega": 0.25,
    "rho0": 1770.0,
    "Q": 1.01e10 / 1770.0,
    "cv": 1000.0,
}
PETN_KUHL = {  # Kuhl, Bell, Beckner, Khasainov, UCRL-PROC-225822 (2006) Appendix;
    # rho0 is KUHL'S OWN reference/charge density (1.0 g/cc pressed booster, not
    # full PETN TMD 1.77 g/cc) -- A/B are fit to that isochor, so use them together.
    # cv derived from R_DP = 28.76 g/mol via cv = R/omega; Q from the reported heat
    # of detonation 1423 cal/g; T_CJ ~ 4600 K (no separate D_CJ/P_CJ reported).
    "name": "PETN (Kuhl)",
    "A": 5.8e11,
    "B": 9.3e9,
    "R1": 7.0,
    "R2": 1.7,
    "omega": 0.246,
    "rho0": 1000.0,
    "Q": 5.9538e6,
    "cv": 1174.9,
}
COMP_B = {  # Dobratz & Crawford, UCRL-52997 (Composition B-3); P_CJ ~ 29.5 GPa, D ~ 7980 m/s
    "name": "Composition B",
    "A": 5.242e11,
    "B": 7.678e9,
    "R1": 4.20,
    "R2": 1.10,
    "omega": 0.34,
    "rho0": 1717.0,
    "Q": 8.5e9 / 1717.0,
    "cv": 1000.0,
}
HMX = {  # Dobratz & Crawford, UCRL-52997; P_CJ ~ 42 GPa, D ~ 9110 m/s
    "name": "HMX",
    "A": 7.783e11,
    "B": 7.071e9,
    "R1": 4.20,
    "R2": 1.00,
    "omega": 0.30,
    "rho0": 1891.0,
    "Q": 10.5e9 / 1891.0,
    "cv": 1000.0,
}
PBX9404 = {  # Dobratz & Crawford, UCRL-52997; P_CJ ~ 37 GPa, D ~ 8800 m/s
    "name": "PBX-9404",
    "A": 8.524e11,
    "B": 1.802e10,
    "R1": 4.60,
    "R2": 1.30,
    "omega": 0.38,
    "rho0": 1840.0,
    "Q": 10.2e9 / 1840.0,
    "cv": 1000.0,
}

PRODUCTS = {p["name"]: p for p in (TNT, LX10, PETN, PETN_KUHL, COMP_B, HMX, PBX9404)}

# --- ambient gas (the co-existing non-JWL fluid) ---
# gamma is MFC's convention (1/(gamma_phys - 1)); e0 [J/kg] is the energy at p0, rho0.
AIR = {"name": "air", "rho0": 1.225, "e0": 2.5575e5, "gamma": 2.5, "cv": 717.5}


def jwl_fluid(idx, products, ambient=AIR, delta_e=None):
    """fluid_pp(idx)%... entries for a JWL products fluid and its ambient references.

    delta_e (<= 0) installs the Garno reactant energy offset for a resolved ZND
    detonation; use znd_delta_e(products) to compute it.
    """
    px = f"fluid_pp({idx})%"
    d = {
        px + "eos": 2,
        px + "gamma": 2.5,  # unused when an ambient fluid is present; fallback ambient Grueneisen otherwise
        px + "pi_inf": 0.0,
        px + "cv": products["cv"],
        px + "jwl_A": products["A"],
        px + "jwl_B": products["B"],
        px + "jwl_R1": products["R1"],
        px + "jwl_R2": products["R2"],
        px + "jwl_omega": products["omega"],
        px + "jwl_rho0": products["rho0"],
        px + "jwl_Q": products["Q"],
        px + "jwl_air_rho0": ambient["rho0"],
        px + "jwl_air_e0": ambient["e0"],
    }
    if delta_e is not None:
        d[px + "jwl_delta_e"] = delta_e
    return d


def ambient_fluid(idx, ambient=AIR):
    """fluid_pp(idx)%... entries for the co-existing ideal-gas ambient fluid."""
    px = f"fluid_pp({idx})%"
    return {px + "eos": 1, px + "gamma": ambient["gamma"], px + "pi_inf": 0.0, px + "cv": ambient["cv"]}


def znd_delta_e(products, p_amb=101325.0):
    """Garno et al. (2020) Eq. 17 reactant energy offset [J/kg], <= 0.

    Sets unreacted (lambda = 0) explosive at (rho0, e = 0) to ambient pressure so
    it sits on a stiffer Hugoniot than the products, giving jwl_reactive a
    resolved ZND structure (von Neumann spike -> CJ).
    """
    a, b, r1, r2, w, rho0 = (products[k] for k in ("A", "B", "R1", "R2", "omega", "rho0"))
    f1 = a * (1.0 - w / r1) * math.exp(-r1) + b * (1.0 - w / r2) * math.exp(-r2)
    return (p_amb - f1) / (w * rho0)
