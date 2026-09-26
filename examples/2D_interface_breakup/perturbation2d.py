"""
perturbation2d.py -- multimode 2D interfacial perturbation generator for breakup studies.

2D analogue of perturbation3d.py. Generates a single-valued interface displacement
eta(x) on a periodic horizontal domain as a superposition of 1D Fourier modes with
random phases and a prescribed power spectrum. This is the standard "multimode"
initialization used in Rayleigh-Taylor / Faraday / primary atomization studies
(e.g. Dimonte et al., Phys. Fluids 2004; Thornber et al., Phys. Fluids 2010).

The interface is   y = y_mean + eta(x),   with

    eta(x) = sum_modes a_m cos(k_m x + phi_m)

Each mode is an integer number of wavelengths across Lx, so eta is exactly
periodic -- matching periodic (bc = -1) lateral boundaries in MFC. Modes are
selected by their *physical* wavenumber k in the band
[2*pi/lambda_max, 2*pi/lambda_min]. Because the band is set by physical
wavelengths rather than integer mode counts, the seeded spectrum -- and the
shortest seeded wave that sets the cells-per-wave resolution -- is independent of
the box dimension Lx.

Statistical quantities that define the perturbation (same meaning as the 3D code)
eta_rms      target RMS amplitude  sqrt(<eta^2>)            [length]
lambda_min,lambda_max
             physical wavelength band of the seeded modes [length]; the seeded
             k spans [2*pi/lambda_max, 2*pi/lambda_min], independent of the box
p            spectral slope:  modal energy  a_m^2 ~ k_m^p
                p =  0  white / flat ("broadband"); -2 red; +2 blue
seed         RNG seed -> reproducible phases (and amplitudes if randomized)
randomize_amp  Rayleigh-random modal amplitudes (Gaussian random field) if True

Diagnostics returned
eta_rms_realized   sqrt(mean(eta^2))            -- check vs requested eta_rms
eta_mean           mean(eta)                    -- ~0 by construction
eta_max,eta_min    extrema (peak-to-valley)
slope_rms          sqrt(<(deta/dx)^2>)          -- interface steepness
lambda_dom         dominant wavelength = 2*pi / k_peak
lambda_int         integral length = 2*pi / k_bar (energy-weighted mean k)
n_modes            number of modes retained in the band
"""

import numpy as np


def generate_perturbation_2d(
    Lx,
    Nx,
    eta_rms,
    lambda_min,
    lambda_max,
    p=0.0,
    seed=0,
    randomize_amp=True,
    x0=0.0,
):
    """Return (x, eta, diag) for a multimode periodic interface eta(x).

    Parameters
    Lx           horizontal domain length (must match the MFC x-domain width)
    Nx           sample count (use the MFC cell count m+1 for a 1:1 map)
    eta_rms      target RMS amplitude of the interface
    lambda_min,lambda_max
                 physical wavelength band of the seeded modes (same length units
                 as Lx). Modes are selected by physical k in
                 [2*pi/lambda_max, 2*pi/lambda_min], so the seeded wavelengths are
                 a property of the physics -- NOT of the box. Changing Lx no
                 longer rescales the shortest seeded wave (and hence the
                 cells-per-wave resolution). Modes remain integer-commensurate
                 with Lx for exact periodicity.
    p            spectral slope, energy a_m^2 ~ k_m^p
    seed         RNG seed for reproducibility
    randomize_amp  Rayleigh-random amplitudes (Gaussian field) if True
    x0           domain origin; samples are cell centers on [x0, x0+Lx), so
                 they line up with MFC's x_cc absolute coords

    Returns
    x    : (Nx,) sample coordinates, cell centers on [x0, x0 + Lx)
    eta  : (Nx,) interface displacement
    diag : dict of statistical diagnostics + modal arrays
    """
    rng = np.random.default_rng(seed)

    # Cell-centered sample points on the periodic horizontal domain
    x = x0 + (np.arange(Nx) + 0.5) * (Lx / Nx)

    # Physical wavenumber band: select modes whose k falls in
    # [k_min, k_max] = [2*pi/lambda_max, 2*pi/lambda_min].  This is anchored to
    # physical lengths, so the seeded band does not move with the box size.
    k_min = 2.0 * np.pi / lambda_max
    k_max = 2.0 * np.pi / lambda_min

    # Enumerate integer modes n >= 1 -- so eta is exactly periodic on Lx -- and
    # keep those whose physical wavenumber lies in the band. The integer range
    # is wide enough that the largest representable k reaches k_max.
    nx_max = int(np.ceil(Lx / lambda_min))
    n_arr = np.arange(1, nx_max + 1, dtype=float)
    k_all = 2.0 * np.pi * n_arr / Lx
    keep = (k_all >= k_min) & (k_all <= k_max)
    kmag = k_all[keep]
    if kmag.size == 0:
        raise ValueError("no integer modes commensurate with Lx fall in the " "wavelength band [lambda_min, lambda_max]; widen the band or the box")

    # Target spectral envelope: a_m^2 ~ k^p  ->  a_m ~ k^(p/2)
    envelope = kmag ** (p / 2.0)

    # Random phases; optional Rayleigh amplitudes for a true Gaussian random field
    phi = rng.uniform(0.0, 2.0 * np.pi, size=kmag.size)
    if randomize_amp:
        a = envelope * np.sqrt(-2.0 * np.log(rng.uniform(0.0, 1.0, size=kmag.size)))
    else:
        a = envelope.copy()

    # Normalize so that <eta^2> = eta_rms^2.
    # Distinct cosine modes over the periodic box: <eta^2> = 0.5 * sum(a_m^2).
    rms_raw = np.sqrt(0.5 * np.sum(a**2))
    a *= eta_rms / rms_raw

    # Synthesize the interface on the x grid
    phase = np.outer(x, kmag) + phi[None, :]  # (Nx, n_modes)
    eta = np.sum(a[None, :] * np.cos(phase), axis=1)
    deta_dx = np.sum(-a[None, :] * kmag[None, :] * np.sin(phase), axis=1)

    # Diagnostics
    Pk = 0.5 * a**2  # modal energy

    k_peak = kmag[np.argmax(Pk)]
    k_bar = np.sum(kmag * Pk) / np.sum(Pk)

    diag = {
        "eta_rms_target": eta_rms,
        "eta_rms_realized": float(np.sqrt(np.mean(eta**2))),
        "eta_mean": float(np.mean(eta)),
        "eta_max": float(np.max(eta)),
        "eta_min": float(np.min(eta)),
        "slope_rms": float(np.sqrt(np.mean(deta_dx**2))),
        "lambda_dom": float(2.0 * np.pi / k_peak),
        # Realized band (shortest/longest seeded wave) from the actual modes
        "lambda_min": float(2.0 * np.pi / kmag.max()),
        "lambda_max": float(2.0 * np.pi / kmag.min()),
        "lambda_int": float(2.0 * np.pi / k_bar),
        "n_modes": int(kmag.size),
        # modal arrays
        "k": kmag,
        "a": a,
        "phi": phi,
        "Pk": Pk,
    }
    return x, eta, diag


def cells_per_lambda(lambda_min, dx):
    """Cells resolving the *shortest* seeded wave: lambda_min / dx.

    With lambda_min anchored to a physical length (not the box) and an isotropic
    grid (dx == dy), this depends only on the grid spacing -- it does not
    drift with the domain aspect ratio Lx/Ly.  Keep >= ~8-12 for WENO5.
    """
    return lambda_min / dx


def write_interface_2d(x, eta, y_mean, path="interface_profile.dat"):
    """Write the *absolute* interface height y_int(x) = y_mean + eta(x).

    This is the file hcid 209 (2dHardcodedIC.fpp) reads. Format: one row per
    sample, two columns
        x   y_interface
    on a UNIFORM x grid (hcid 209 infers the spacing from the first two rows).
    There is no header; the Fortran reader counts rows to EOF.
    """
    y_int = y_mean + eta
    with open(path, "w") as fh:
        np.savetxt(fh, np.column_stack([x, y_int]))


_DIAG_KEYS = [
    "eta_rms_target",
    "eta_rms_realized",
    "eta_mean",
    "eta_max",
    "eta_min",
    "slope_rms",
    "lambda_dom",
    "lambda_min",
    "lambda_max",
    "lambda_int",
    "n_modes",
]


def print_diagnostics(diag, file=None):
    print("2D interface perturbation diagnostics", file=file)
    print("-" * 40, file=file)
    for key in _DIAG_KEYS:
        print(f"  {key:<18s} {diag[key]:.6g}", file=file)


def save_diagnostics(diag, path="diagnostics.txt", extra=None):
    """Write the perturbation diagnostics (and optional extra key/values) to a file."""
    with open(path, "w") as fh:
        fh.write("2D interface perturbation diagnostics\n")
        fh.write("-" * 40 + "\n")
        for key in _DIAG_KEYS:
            fh.write(f"  {key:<18s} {diag[key]:.6g}\n")
        if extra:
            fh.write("\ncase parameters\n")
            fh.write("-" * 40 + "\n")
            for key, val in extra.items():
                if isinstance(val, (int, float)):
                    fh.write(f"  {key:<18s} {val:.6g}\n")
                else:
                    fh.write(f"  {key:<18s} {val}\n")
    return path


if __name__ == "__main__":
    # Demo: broadband 2D perturbation.
    Lx = 2.0
    Nx = 256
    eta_rms = 0.005 * Lx
    H = 1.0
    lambda_min, lambda_max = H / 8.0, H  # physical band, anchored to H
    p = 0.0

    x, eta, diag = generate_perturbation_2d(
        Lx,
        Nx,
        eta_rms,
        lambda_min=lambda_min,
        lambda_max=lambda_max,
        p=p,
        seed=42,
        randomize_amp=True,
    )

    print_diagnostics(diag)
    print(f"  cells/shortest-wave  {cells_per_lambda(diag['lambda_min'], Lx / Nx):.1f}")
    write_interface_2d(x, eta, y_mean=1.0)
    print("\nwrote interface_profile.dat")

    # Optional plot (headless-safe)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.5, 3.0))
        ax.plot(x, eta)
        ax.set_xlabel("x")
        ax.set_ylabel("eta")
        ax.set_title(f"eta(x) (rms={diag['eta_rms_realized']:.3g})")
        fig.tight_layout()
        fig.savefig("interface_profile.png", dpi=150)
        print("wrote interface_profile.png")
    except Exception as exc:  # pragma: no cover
        print(f"(plot skipped: {exc})")
