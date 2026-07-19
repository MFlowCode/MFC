# JWL equation of state (simple form)

This build adds a Jones-Wilkins-Lee (JWL) equation of state for detonation
products to MFC's five-equation model (`model_eqns = 2`), together with a simple
two-material mixture closure. It is the first tier of a staged contribution: one
JWL products fluid plus one ambient gas, a linear mass-fraction mixture rule, and
no reactions. The reactive and composition-weighted layers are listed under
Roadmap below and land as later pull requests.

## Model

Pure JWL products follow the Mie-Grueneisen form referenced to an isentrope,
with `V = rho0/rho`:

```
p = A (1 - w/(R1 V)) exp(-R1 V) + B (1 - w/(R2 V)) exp(-R2 V) + w rho e
```

`A`, `B` (Pa), `R1`, `R2`, and `w = omega` are cylinder-test fits. A cell holding
a products mass fraction `Y = alpha_rho_products / rho` mixes products and ambient
by blending the coefficients linearly in `Y`:

```
An = Y A,   Bn = Y B,   omega = air_gamma + Y (omega0 - air_gamma),
cv  = Y cv_products + (1 - Y) cv_air
```

The blend is exact at the endpoints: `Y = 1` recovers the pure JWL law and
`Y = 0` recovers the ambient law. A stiffened-gas ambient (for example water)
adds `pi_hat = (1 - Y) pi_inf`; an ideal-gas ambient (`pi_inf = 0`) recovers the
ideal closure bit-identically. Because every coefficient depends on `Y` alone and
never on `rho` or `e`, the pressure stays linear in `e`, so the `(rho, p, Y) -> e`
inverse is closed form and the sound speed is the exact Grueneisen derivative.
At start-up `s_jwl_verify_closure` sweeps the `(rho, e, Y)` envelope and aborts on
any non-positive sound speed or failed pressure/energy round trip. The full
derivation lives in `src/common/m_jwl.fpp`.

## Case setup

Select the EOS per fluid with `fluid_pp(i)%eos`: `1` stiffened gas (default), `2`
JWL. At most one fluid may be JWL, and it requires `model_eqns = 2`. A JWL fluid
is defined by:

| Parameter | Meaning | Requirement |
| :--- | :--- | :--- |
| `jwl_A`, `jwl_B`, `jwl_R1`, `jwl_R2`, `jwl_omega` | JWL products coefficients | required |
| `jwl_rho0` | products reference density | required |
| `jwl_Q` or `jwl_E0` | detonation energy, specific (J/kg) or volumetric (J/m^3); given `jwl_Q`, MFC sets `jwl_E0 = jwl_rho0 * jwl_Q` | one of the two |
| `jwl_air_rho0` | density of the co-existing ambient gas | required |
| `jwl_air_e0` or `jwl_air_p0` | ambient specific internal energy or pressure | one of the two |
| `cv` | positive heat capacity, on both the JWL fluid and the ambient fluid | required |

The ambient Grueneisen coefficient `omega0_air = 1/fluid_pp(air)%gamma` comes from
the non-JWL fluid (or the JWL fluid's own `gamma` when it is alone). Recall that
`fluid_pp%gamma` stores `1/(gamma - 1)`, so air is `2.5`, not `1.4`.

JWL is prohibited with `wave_speeds = 2`, characteristic (CBC) boundaries,
`alt_soundspeed`, elasticity, `igr`, `bubbles_euler`, `mhd`, and `chemistry`;
each of those paths evaluates a stiffened-gas pressure that would bypass the JWL
closure. Immersed boundaries (`ib`) and Lagrangian bubbles (`bubbles_lagrange`)
are supported.

## Examples

- `examples/1D_jwl_single_material_shocktube` — single JWL fluid, two-pressure shock tube.
- `examples/1D_jwl_mixture_test` — products and air across the linear closure.
- `examples/1D_jwl_air_interface_advection` — advected products/air interface.
- `examples/2D_jwl_detonation` — cylindrical products charge bursting into air.

## Roadmap (Tier 2, later pull requests)

- Composition (heat-capacity) weighted closure `w = Y cv_j / (Y cv_j + (1 - Y) cv_a)`
  in place of the linear `w = Y` blend.
- Reactant/product energy offset (`jwl_delta_e`) for resolved ZND detonation structure.
- Reaction sources: program burn, products-air afterburn, and JWL++ pressure-driven reactive burn.
- Reactive-mode immersed-boundary and Lagrangian-bubble coupling.
