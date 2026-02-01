@page parameters Case Parameters Reference

# Case Parameters Reference

> **Auto-generated** from parameter registry
> 
> Regenerate with: `./mfc.sh generate --json-schema`

## Overview

MFC supports **3,400** case parameters organized into families.

**Quick search:** Use `./mfc.sh params <query>` to search parameters from the command line.

## Parameter Families

| Family | Count | Description |
|--------|-------|-------------|
| [patch_icpp](#patch_icpp) | 1954 | Initial condition patch parameters |
| [patch_ib](#patch_ib) | 370 | Immersed boundary patch parameters |
| [fluid_pp](#fluid_pp) | 170 | Fluid material properties |
| [patch_bc](#patch_bc) | 110 | Boundary condition patch parameters |
| [acoustic](#acoustic) | 104 | Acoustic source parameters |
| [simplex_params](#simplex_params) | 78 | Simplex noise perturbation parameters |
| [bc_x](#bc_x) | 39 | X-direction boundary conditions |
| [bc_y](#bc_y) | 39 | Y-direction boundary conditions |
| [bc_z](#bc_z) | 39 | Z-direction boundary conditions |
| [integral](#integral) | 30 | Integral region parameters |
| [probe](#probe) | 30 | Probe/monitoring point parameters |
| [bub_pp](#bub_pp) | 20 | Bubble property parameters |
| [lag_params](#lag_params) | 17 | Lagrangian particle parameters |
| [alpha_rho_e_wrt](#alpha_rho_e_wrt) | 10 | Partial density-energy output flags |
| [alpha_rho_wrt](#alpha_rho_wrt) | 10 | Partial density output flags |
| [alpha_wrt](#alpha_wrt) | 10 | Volume fraction output flags |
| [fluid_rho](#fluid_rho) | 10 | Fluid reference densities |
| [kappa_wrt](#kappa_wrt) | 10 | Curvature output flags |
| [schlieren_alpha](#schlieren_alpha) | 10 | Numerical schlieren coefficients |
| [chem_params](#chem_params) | 4 | Chemistry model parameters |
| [flux_wrt](#flux_wrt) | 3 | Flux output flags |
| [mom_wrt](#mom_wrt) | 3 | Momentum output flags |
| [omega_wrt](#omega_wrt) | 3 | Vorticity output flags |
| [vel_wrt](#vel_wrt) | 3 | Velocity output flags |
| [x_domain](#x_domain) | 2 | X-direction domain parameters |
| [x_output](#x_output) | 2 | X-direction output region |
| [y_domain](#y_domain) | 2 | Y-direction domain parameters |
| [y_output](#y_output) | 2 | Y-direction output region |
| [z_domain](#z_domain) | 2 | Z-direction domain parameters |
| [z_output](#z_output) | 2 | Z-direction output region |
| [general](#general) | 312 | Core simulation parameters (grid, time, model, etc.) |

---

## patch_icpp

*Initial condition patch parameters*

**1954 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `patch_icpp(N)%%Bx` | `patch_icpp(1)%%Bx` | X-component of magnetic field for patch 1 |
| `patch_icpp(N)%%By` | `patch_icpp(1)%%By` | Y-component of magnetic field for patch 1 |
| `patch_icpp(N)%%Bz` | `patch_icpp(1)%%Bz` | Z-component of magnetic field for patch 1 |
| `patch_icpp(N)%%Y(M)` | `patch_icpp(1)%%Y(1)` | Mass fraction of species 1 for patch 1 |
| `patch_icpp(N)%%a(M)` | `patch_icpp(1)%%a(2)` | Shape coefficient a(2) for patch 1 |
| `patch_icpp(N)%%alpha` | `patch_icpp(1)%%alpha` | Volume fraction for patch 1 |
| `patch_icpp(N)%%alpha(M)` | `patch_icpp(1)%%alpha(1)` | Volume fraction of fluid 1 for patch 1 |
| `patch_icpp(N)%%alpha_rho` | `patch_icpp(1)%%alpha_rho` | Partial density for patch 1 |
| `patch_icpp(N)%%alpha_rho(M)` | `patch_icpp(1)%%alpha_rho(1)` | Partial density of fluid 1 for patch 1 |
| `patch_icpp(N)%%alter_patch` | `patch_icpp(10)%%alter_patch` | Enable patch alteration for patch 10 |
| `patch_icpp(N)%%alter_patch(M)` | `patch_icpp(10)%%alter_patch(1)` | Alter patch 1 with patch 10 |
| `patch_icpp(N)%%beta` | `patch_icpp(1)%%beta` | Shape parameter beta for patch 1 |
| `patch_icpp(N)%%cf_val` | `patch_icpp(1)%%cf_val` | Color function value for patch 1 |
| `patch_icpp(N)%%cv` | `patch_icpp(1)%%cv` | Specific heat at constant volume for patch 1 |
| `patch_icpp(N)%%epsilon` | `patch_icpp(1)%%epsilon` | Interface thickness for patch 1 |
| `patch_icpp(N)%%gamma` | `patch_icpp(1)%%gamma` | Specific heat ratio for patch 1 |
| `patch_icpp(N)%%geometry` | `patch_icpp(1)%%geometry` | Geometry type for initial condition patch 1 |
| `patch_icpp(N)%%hcid` | `patch_icpp(1)%%hcid` | Hard-coded patch ID for patch 1 |
| `patch_icpp(N)%%length_x` | `patch_icpp(1)%%length_x` | X-dimension length for patch 1 |
| `patch_icpp(N)%%length_y` | `patch_icpp(1)%%length_y` | Y-dimension length for patch 1 |
| `patch_icpp(N)%%length_z` | `patch_icpp(1)%%length_z` | Z-dimension length for patch 1 |
| `patch_icpp(N)%%m0` | `patch_icpp(1)%%m0` | Initial bubble mass for patch 1 |
| `patch_icpp(N)%%model_filepath` | `patch_icpp(1)%%model_filepath` | STL model file path for patch 1 |
| `patch_icpp(N)%%model_rotate(M)` | `patch_icpp(1)%%model_rotate(1)` | Model rotation component 1 for patch 1 |
| `patch_icpp(N)%%model_scale(M)` | `patch_icpp(1)%%model_scale(1)` | Model scale component 1 for patch 1 |
| `patch_icpp(N)%%model_spc` | `patch_icpp(1)%%model_spc` | Model spacing for patch 1 |
| `patch_icpp(N)%%model_threshold` | `patch_icpp(1)%%model_threshold` | Model threshold for patch 1 |
| `patch_icpp(N)%%model_translate(M)` | `patch_icpp(1)%%model_translate(1)` | Model translation component 1 for patch 1 |
| `patch_icpp(N)%%non_axis_sym` | `patch_icpp(1)%%non_axis_sym` | Non-axisymmetric parameter for patch 1 |
| `patch_icpp(N)%%normal` | `patch_icpp(1)%%normal` | Normal direction for patch 1 |
| `patch_icpp(N)%%normal(M)` | `patch_icpp(1)%%normal(1)` | Normal component 1 for patch 1 |
| `patch_icpp(N)%%p0` | `patch_icpp(1)%%p0` | Initial bubble pressure for patch 1 |
| `patch_icpp(N)%%pi_inf` | `patch_icpp(1)%%pi_inf` | Stiffness pressure for patch 1 |
| `patch_icpp(N)%%pres` | `patch_icpp(1)%%pres` | Pressure for patch 1 |
| `patch_icpp(N)%%qv` | `patch_icpp(1)%%qv` | Heat of formation for patch 1 |
| `patch_icpp(N)%%qvp` | `patch_icpp(1)%%qvp` | Heat of formation prime for patch 1 |
| `patch_icpp(N)%%r0` | `patch_icpp(1)%%r0` | Initial bubble radius for patch 1 |
| `patch_icpp(N)%%radii` | `patch_icpp(1)%%radii` | Radii for patch 1 |
| `patch_icpp(N)%%radii(M)` | `patch_icpp(1)%%radii(1)` | Radius component 1 for patch 1 |
| `patch_icpp(N)%%radius` | `patch_icpp(1)%%radius` | Radius for patch 1 |
| `patch_icpp(N)%%rho` | `patch_icpp(1)%%rho` | Density for patch 1 |
| `patch_icpp(N)%%smooth_coeff` | `patch_icpp(1)%%smooth_coeff` | Smoothing coefficient for patch 1 |
| `patch_icpp(N)%%smooth_patch_id` | `patch_icpp(1)%%smooth_patch_id` | Patch ID to smooth against for patch 1 |
| `patch_icpp(N)%%smoothen` | `patch_icpp(1)%%smoothen` | Enable smoothing for patch 1 |
| `patch_icpp(N)%%tau_e(M)` | `patch_icpp(1)%%tau_e(1)` | Elastic stress component 1 for patch 1 |
| `patch_icpp(N)%%v0` | `patch_icpp(1)%%v0` | Initial bubble velocity for patch 1 |
| `patch_icpp(N)%%vel` | `patch_icpp(1)%%vel` | Velocity magnitude for patch 1 |
| `patch_icpp(N)%%vel(M)` | `patch_icpp(1)%%vel(1)` | Velocity component 1 for patch 1 |
| `patch_icpp(N)%%x_centroid` | `patch_icpp(1)%%x_centroid` | X-coordinate of centroid for patch 1 |
| `patch_icpp(N)%%y_centroid` | `patch_icpp(1)%%y_centroid` | Y-coordinate of centroid for patch 1 |
| `patch_icpp(N)%%z_centroid` | `patch_icpp(1)%%z_centroid` | Z-coordinate of centroid for patch 1 |

---

## patch_ib

*Immersed boundary patch parameters*

**370 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `patch_ib(N)%%angles(M)` | `patch_ib(1)%%angles(1)` | Orientation angle 1 for IB patch 1 |
| `patch_ib(N)%%angular_vel(M)` | `patch_ib(1)%%angular_vel(1)` | Angular velocity component 1 for IB patch 1 |
| `patch_ib(N)%%c` | `patch_ib(1)%%c` | Shape parameter c for IB patch 1 |
| `patch_ib(N)%%geometry` | `patch_ib(1)%%geometry` | Geometry type for immersed boundary 1 |
| `patch_ib(N)%%length_x` | `patch_ib(1)%%length_x` | X-dimension length for IB patch 1 |
| `patch_ib(N)%%length_y` | `patch_ib(1)%%length_y` | Y-dimension length for IB patch 1 |
| `patch_ib(N)%%length_z` | `patch_ib(1)%%length_z` | Z-dimension length for IB patch 1 |
| `patch_ib(N)%%m` | `patch_ib(1)%%m` | Shape parameter m for IB patch 1 |
| `patch_ib(N)%%mass` | `patch_ib(1)%%mass` | Mass for IB patch 1 |
| `patch_ib(N)%%model_filepath` | `patch_ib(1)%%model_filepath` | STL model file path for IB patch 1 |
| `patch_ib(N)%%model_rotate(M)` | `patch_ib(1)%%model_rotate(1)` | Model rotation component 1 for IB patch 1 |
| `patch_ib(N)%%model_scale(M)` | `patch_ib(1)%%model_scale(1)` | Model scale component 1 for IB patch 1 |
| `patch_ib(N)%%model_spc` | `patch_ib(1)%%model_spc` | Model spacing for IB patch 1 |
| `patch_ib(N)%%model_threshold` | `patch_ib(1)%%model_threshold` | Model threshold for IB patch 1 |
| `patch_ib(N)%%model_translate(M)` | `patch_ib(1)%%model_translate(1)` | Model translation component 1 for IB patch 1 |
| `patch_ib(N)%%moving_ibm` | `patch_ib(1)%%moving_ibm` | Enable moving boundary for IB patch 1 |
| `patch_ib(N)%%p` | `patch_ib(1)%%p` | Shape parameter p for IB patch 1 |
| `patch_ib(N)%%radius` | `patch_ib(1)%%radius` | Radius for IB patch 1 |
| `patch_ib(N)%%slip` | `patch_ib(1)%%slip` | Enable slip condition for IB patch 1 |
| `patch_ib(N)%%t` | `patch_ib(1)%%t` | Shape parameter t for IB patch 1 |
| `patch_ib(N)%%theta` | `patch_ib(1)%%theta` | Theta angle for IB patch 1 |
| `patch_ib(N)%%vel(M)` | `patch_ib(1)%%vel(1)` | Velocity component 1 for IB patch 1 |
| `patch_ib(N)%%x_centroid` | `patch_ib(1)%%x_centroid` | X-coordinate of centroid for IB patch 1 |
| `patch_ib(N)%%y_centroid` | `patch_ib(1)%%y_centroid` | Y-coordinate of centroid for IB patch 1 |
| `patch_ib(N)%%z_centroid` | `patch_ib(1)%%z_centroid` | Z-coordinate of centroid for IB patch 1 |

---

## fluid_pp

*Fluid material properties*

**170 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `fluid_pp(N)%%D_v` | `fluid_pp(1)%%D_v` |  |
| `fluid_pp(N)%%G` | `fluid_pp(1)%%G` | Shear modulus for fluid 1 |
| `fluid_pp(N)%%M_v` | `fluid_pp(1)%%M_v` |  |
| `fluid_pp(N)%%Re(M)` | `fluid_pp(1)%%Re(1)` | Reynolds number component 1 for fluid 1 |
| `fluid_pp(N)%%cp_v` | `fluid_pp(1)%%cp_v` |  |
| `fluid_pp(N)%%cv` | `fluid_pp(1)%%cv` | Specific heat at constant volume for fluid 1 |
| `fluid_pp(N)%%gamma` | `fluid_pp(1)%%gamma` | Specific heat ratio for fluid 1 |
| `fluid_pp(N)%%gamma_v` | `fluid_pp(1)%%gamma_v` |  |
| `fluid_pp(N)%%k_v` | `fluid_pp(1)%%k_v` |  |
| `fluid_pp(N)%%mu_v` | `fluid_pp(1)%%mu_v` |  |
| `fluid_pp(N)%%mul0` | `fluid_pp(1)%%mul0` |  |
| `fluid_pp(N)%%pi_inf` | `fluid_pp(1)%%pi_inf` | Stiffness pressure for fluid 1 |
| `fluid_pp(N)%%pv` | `fluid_pp(1)%%pv` |  |
| `fluid_pp(N)%%qv` | `fluid_pp(1)%%qv` | Heat of formation for fluid 1 |
| `fluid_pp(N)%%qvp` | `fluid_pp(1)%%qvp` | Heat of formation prime for fluid 1 |
| `fluid_pp(N)%%ss` | `fluid_pp(1)%%ss` |  |

---

## patch_bc

*Boundary condition patch parameters*

**110 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `patch_bc(N)%%centroid(M)` | `patch_bc(1)%%centroid(1)` | Centroid component 1 for BC patch 1 |
| `patch_bc(N)%%dir` | `patch_bc(1)%%dir` | Direction for BC patch 1 |
| `patch_bc(N)%%geometry` | `patch_bc(1)%%geometry` | Geometry type for BC patch 1 |
| `patch_bc(N)%%length(M)` | `patch_bc(1)%%length(1)` | Length component 1 for BC patch 1 |
| `patch_bc(N)%%loc` | `patch_bc(1)%%loc` | Location for BC patch 1 |
| `patch_bc(N)%%radius` | `patch_bc(1)%%radius` | Radius for BC patch 1 |
| `patch_bc(N)%%type` | `patch_bc(1)%%type` | BC type for patch 1 |

---

## acoustic

*Acoustic source parameters*

**104 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `acoustic(N)%%aperture` | `acoustic(1)%%aperture` | Aperture for acoustic source 1 |
| `acoustic(N)%%bb_bandwidth` | `acoustic(1)%%bb_bandwidth` | Broadband bandwidth for acoustic source 1 |
| `acoustic(N)%%bb_lowest_freq` | `acoustic(1)%%bb_lowest_freq` | Lowest broadband frequency for source 1 |
| `acoustic(N)%%bb_num_freq` | `acoustic(1)%%bb_num_freq` | Number of broadband frequencies for source 1 |
| `acoustic(N)%%delay` | `acoustic(1)%%delay` | Delay for acoustic source 1 |
| `acoustic(N)%%dipole` | `acoustic(1)%%dipole` | Enable dipole for acoustic source 1 |
| `acoustic(N)%%dir` | `acoustic(1)%%dir` | Direction for acoustic source 1 |
| `acoustic(N)%%element_on` | `acoustic(1)%%element_on` | Active element index for acoustic source 1 |
| `acoustic(N)%%element_polygon_ratio` | `acoustic(1)%%element_polygon_ratio` | Element polygon ratio for acoustic source 1 |
| `acoustic(N)%%element_spacing_angle` | `acoustic(1)%%element_spacing_angle` | Element spacing angle for acoustic source 1 |
| `acoustic(N)%%foc_length` | `acoustic(1)%%foc_length` | Focal length for acoustic source 1 |
| `acoustic(N)%%frequency` | `acoustic(1)%%frequency` | Frequency for acoustic source 1 |
| `acoustic(N)%%gauss_sigma_dist` | `acoustic(1)%%gauss_sigma_dist` | Gaussian spatial width for acoustic source 1 |
| `acoustic(N)%%gauss_sigma_time` | `acoustic(1)%%gauss_sigma_time` | Gaussian temporal width for acoustic source 1 |
| `acoustic(N)%%height` | `acoustic(1)%%height` | Height for acoustic source 1 |
| `acoustic(N)%%length` | `acoustic(1)%%length` | Length for acoustic source 1 |
| `acoustic(N)%%loc(M)` | `acoustic(1)%%loc(1)` | Location component 1 for acoustic source 1 |
| `acoustic(N)%%mag` | `acoustic(1)%%mag` | Magnitude for acoustic source 1 |
| `acoustic(N)%%npulse` | `acoustic(1)%%npulse` | Number of pulses for acoustic source 1 |
| `acoustic(N)%%num_elements` | `acoustic(1)%%num_elements` | Number of array elements for acoustic source 1 |
| `acoustic(N)%%pulse` | `acoustic(1)%%pulse` | Pulse type for acoustic source 1 |
| `acoustic(N)%%rotate_angle` | `acoustic(1)%%rotate_angle` | Rotation angle for acoustic source 1 |
| `acoustic(N)%%support` | `acoustic(1)%%support` | Support type for acoustic source 1 |
| `acoustic(N)%%wavelength` | `acoustic(1)%%wavelength` | Wavelength for acoustic source 1 |

---

## simplex_params

*Simplex noise perturbation parameters*

**78 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `simplex_params%%perturb_dens(N)` | `simplex_params%%perturb_dens(1)` | Enable density perturbation for fluid 1 |
| `simplex_params%%perturb_dens_freq(N)` | `simplex_params%%perturb_dens_freq(1)` | Density perturbation frequency for fluid 1 |
| `simplex_params%%perturb_dens_offset(N, M)` | `simplex_params%%perturb_dens_offset(1, 1)` | Density perturbation offset (1) for fluid 1 |
| `simplex_params%%perturb_dens_scale(N)` | `simplex_params%%perturb_dens_scale(1)` | Density perturbation scale for fluid 1 |
| `simplex_params%%perturb_vel(N)` | `simplex_params%%perturb_vel(1)` | Enable velocity perturbation for direction 1 |
| `simplex_params%%perturb_vel_freq(N)` | `simplex_params%%perturb_vel_freq(1)` | Velocity perturbation frequency for direction 1 |
| `simplex_params%%perturb_vel_offset(N, M)` | `simplex_params%%perturb_vel_offset(1,1)` | Velocity perturbation offset (1) for direction 1 |
| `simplex_params%%perturb_vel_scale(N)` | `simplex_params%%perturb_vel_scale(1)` | Velocity perturbation scale for direction 1 |

---

## bc_x

*X-direction boundary conditions*

**39 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `bc_x%%alpha_in(N)` | `bc_x%%alpha_in(1)` | Inlet volume fraction of fluid 1 at x-boundary |
| `bc_x%%alpha_rho_in(N)` | `bc_x%%alpha_rho_in(1)` | Inlet partial density of fluid 1 at x-boundary |
| `bc_x%%beg` | `bc_x%%beg` | Boundary condition at x-begin (-1=periodic, -2=reflective... |
| `bc_x%%end` | `bc_x%%end` | Boundary condition at x-end |
| `bc_x%%grcbc_in` | `bc_x%%grcbc_in` | Enable GRCBC at x-inlet |
| `bc_x%%grcbc_out` | `bc_x%%grcbc_out` | Enable GRCBC at x-outlet |
| `bc_x%%grcbc_vel_out` | `bc_x%%grcbc_vel_out` | Enable GRCBC velocity at x-outlet |
| `bc_x%%pres_in` | `bc_x%%pres_in` | Inlet pressure at x-boundary |
| `bc_x%%pres_out` | `bc_x%%pres_out` | Outlet pressure at x-boundary |
| `bc_x%%vb1` | `bc_x%%vb1` | Boundary velocity component 1 at x-begin |
| `bc_x%%vb2` | `bc_x%%vb2` | Boundary velocity component 2 at x-begin |
| `bc_x%%vb3` | `bc_x%%vb3` | Boundary velocity component 3 at x-begin |
| `bc_x%%ve1` | `bc_x%%ve1` | Boundary velocity component 1 at x-end |
| `bc_x%%ve2` | `bc_x%%ve2` | Boundary velocity component 2 at x-end |
| `bc_x%%ve3` | `bc_x%%ve3` | Boundary velocity component 3 at x-end |
| `bc_x%%vel_in(N)` | `bc_x%%vel_in(1)` | Inlet velocity component 1 at x-boundary |
| `bc_x%%vel_out(N)` | `bc_x%%vel_out(1)` | Outlet velocity component 1 at x-boundary |

---

## bc_y

*Y-direction boundary conditions*

**39 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `bc_y%%alpha_in(N)` | `bc_y%%alpha_in(1)` | Inlet volume fraction of fluid 1 at y-boundary |
| `bc_y%%alpha_rho_in(N)` | `bc_y%%alpha_rho_in(1)` | Inlet partial density of fluid 1 at y-boundary |
| `bc_y%%beg` | `bc_y%%beg` | Boundary condition at y-begin |
| `bc_y%%end` | `bc_y%%end` | Boundary condition at y-end |
| `bc_y%%grcbc_in` | `bc_y%%grcbc_in` | Enable GRCBC at y-inlet |
| `bc_y%%grcbc_out` | `bc_y%%grcbc_out` | Enable GRCBC at y-outlet |
| `bc_y%%grcbc_vel_out` | `bc_y%%grcbc_vel_out` | Enable GRCBC velocity at y-outlet |
| `bc_y%%pres_in` | `bc_y%%pres_in` | Inlet pressure at y-boundary |
| `bc_y%%pres_out` | `bc_y%%pres_out` | Outlet pressure at y-boundary |
| `bc_y%%vb1` | `bc_y%%vb1` | Boundary velocity component 1 at y-begin |
| `bc_y%%vb2` | `bc_y%%vb2` | Boundary velocity component 2 at y-begin |
| `bc_y%%vb3` | `bc_y%%vb3` | Boundary velocity component 3 at y-begin |
| `bc_y%%ve1` | `bc_y%%ve1` | Boundary velocity component 1 at y-end |
| `bc_y%%ve2` | `bc_y%%ve2` | Boundary velocity component 2 at y-end |
| `bc_y%%ve3` | `bc_y%%ve3` | Boundary velocity component 3 at y-end |
| `bc_y%%vel_in(N)` | `bc_y%%vel_in(1)` | Inlet velocity component 1 at y-boundary |
| `bc_y%%vel_out(N)` | `bc_y%%vel_out(1)` | Outlet velocity component 1 at y-boundary |

---

## bc_z

*Z-direction boundary conditions*

**39 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `bc_z%%alpha_in(N)` | `bc_z%%alpha_in(1)` | Inlet volume fraction of fluid 1 at z-boundary |
| `bc_z%%alpha_rho_in(N)` | `bc_z%%alpha_rho_in(1)` | Inlet partial density of fluid 1 at z-boundary |
| `bc_z%%beg` | `bc_z%%beg` | Boundary condition at z-begin |
| `bc_z%%end` | `bc_z%%end` | Boundary condition at z-end |
| `bc_z%%grcbc_in` | `bc_z%%grcbc_in` | Enable GRCBC at z-inlet |
| `bc_z%%grcbc_out` | `bc_z%%grcbc_out` | Enable GRCBC at z-outlet |
| `bc_z%%grcbc_vel_out` | `bc_z%%grcbc_vel_out` | Enable GRCBC velocity at z-outlet |
| `bc_z%%pres_in` | `bc_z%%pres_in` | Inlet pressure at z-boundary |
| `bc_z%%pres_out` | `bc_z%%pres_out` | Outlet pressure at z-boundary |
| `bc_z%%vb1` | `bc_z%%vb1` | Boundary velocity component 1 at z-begin |
| `bc_z%%vb2` | `bc_z%%vb2` | Boundary velocity component 2 at z-begin |
| `bc_z%%vb3` | `bc_z%%vb3` | Boundary velocity component 3 at z-begin |
| `bc_z%%ve1` | `bc_z%%ve1` | Boundary velocity component 1 at z-end |
| `bc_z%%ve2` | `bc_z%%ve2` | Boundary velocity component 2 at z-end |
| `bc_z%%ve3` | `bc_z%%ve3` | Boundary velocity component 3 at z-end |
| `bc_z%%vel_in(N)` | `bc_z%%vel_in(1)` | Inlet velocity component 1 at z-boundary |
| `bc_z%%vel_out(N)` | `bc_z%%vel_out(1)` | Outlet velocity component 1 at z-boundary |

---

## integral

*Integral region parameters*

**30 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `integral(N)%%xmax` | `integral(1)%%xmax` | X-max of integral region 1 |
| `integral(N)%%xmin` | `integral(1)%%xmin` | X-min of integral region 1 |
| `integral(N)%%ymax` | `integral(1)%%ymax` | Y-max of integral region 1 |
| `integral(N)%%ymin` | `integral(1)%%ymin` | Y-min of integral region 1 |
| `integral(N)%%zmax` | `integral(1)%%zmax` | Z-max of integral region 1 |
| `integral(N)%%zmin` | `integral(1)%%zmin` | Z-min of integral region 1 |

---

## probe

*Probe/monitoring point parameters*

**30 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `probe(N)%%x` | `probe(1)%%x` | X-coordinate of probe 1 |
| `probe(N)%%y` | `probe(1)%%y` | Y-coordinate of probe 1 |
| `probe(N)%%z` | `probe(1)%%z` | Z-coordinate of probe 1 |

---

## bub_pp

*Bubble property parameters*

**20 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `bub_pp%%M_g` | Real | Bubble parameter: M_g |
| `bub_pp%%M_v` | Real | Bubble parameter: M_v |
| `bub_pp%%R0ref` | Real | Bubble parameter: R0ref |
| `bub_pp%%R_g` | Real | Bubble parameter: R_g |
| `bub_pp%%R_v` | Real | Bubble parameter: R_v |
| `bub_pp%%T0ref` | Real | Bubble parameter: T0ref |
| `bub_pp%%cp_g` | Real | Bubble parameter: cp_g |
| `bub_pp%%cp_v` | Real | Bubble parameter: cp_v |
| `bub_pp%%gam_g` | Real | Bubble parameter: gam_g |
| `bub_pp%%gam_v` | Real | Bubble parameter: gam_v |
| `bub_pp%%k_g` | Real | Bubble parameter: k_g |
| `bub_pp%%k_v` | Real | Bubble parameter: k_v |
| `bub_pp%%mu_g` | Real | Bubble parameter: mu_g |
| `bub_pp%%mu_l` | Real | Bubble parameter: mu_l |
| `bub_pp%%mu_v` | Real | Bubble parameter: mu_v |
| `bub_pp%%p0ref` | Real | Bubble parameter: p0ref |
| `bub_pp%%pv` | Real | Bubble parameter: pv |
| `bub_pp%%rho0ref` | Real | Bubble parameter: rho0ref |
| `bub_pp%%ss` | Real | Bubble parameter: ss |
| `bub_pp%%vd` | Real | Bubble parameter: vd |

---

## lag_params

*Lagrangian particle parameters*

**17 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `lag_params%%T0` | Real | Lagrangian tracking parameter: T0 |
| `lag_params%%Thost` | Real | Lagrangian tracking parameter: Thost |
| `lag_params%%c0` | Real | Lagrangian tracking parameter: c0 |
| `lag_params%%charwidth` | Real | Lagrangian tracking parameter: charwidth |
| `lag_params%%cluster_type` | Integer | Lagrangian tracking parameter: cluster_type |
| `lag_params%%epsilonb` | Real | Lagrangian tracking parameter: epsilonb |
| `lag_params%%heatTransfer_model` | Logical (T/F) | Lagrangian tracking parameter: heatTransfer_model |
| `lag_params%%massTransfer_model` | Logical (T/F) | Lagrangian tracking parameter: massTransfer_model |
| `lag_params%%nBubs_glb` | Integer | Lagrangian tracking parameter: nBubs_glb |
| `lag_params%%pressure_corrector` | Logical (T/F) | Lagrangian tracking parameter: pressure_corrector |
| `lag_params%%rho0` | Real | Lagrangian tracking parameter: rho0 |
| `lag_params%%smooth_type` | Integer | Lagrangian tracking parameter: smooth_type |
| `lag_params%%solver_approach` | Integer | Lagrangian tracking parameter: solver_approach |
| `lag_params%%valmaxvoid` | Real | Lagrangian tracking parameter: valmaxvoid |
| `lag_params%%write_bubbles` | Logical (T/F) | Lagrangian tracking parameter: write_bubbles |
| `lag_params%%write_bubbles_stats` | Logical (T/F) | Lagrangian tracking parameter: write_bubbles_stats |
| `lag_params%%x0` | Real | Lagrangian tracking parameter: x0 |

---

## alpha_rho_e_wrt

*Partial density-energy output flags*

**10 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `alpha_rho_e_wrt(1)` | Logical (T/F) | Write partial energy for fluid 1 |
| `alpha_rho_e_wrt(10)` | Logical (T/F) | Write partial energy for fluid 10 |
| `alpha_rho_e_wrt(2)` | Logical (T/F) | Write partial energy for fluid 2 |
| `alpha_rho_e_wrt(3)` | Logical (T/F) | Write partial energy for fluid 3 |
| `alpha_rho_e_wrt(4)` | Logical (T/F) | Write partial energy for fluid 4 |
| `alpha_rho_e_wrt(5)` | Logical (T/F) | Write partial energy for fluid 5 |
| `alpha_rho_e_wrt(6)` | Logical (T/F) | Write partial energy for fluid 6 |
| `alpha_rho_e_wrt(7)` | Logical (T/F) | Write partial energy for fluid 7 |
| `alpha_rho_e_wrt(8)` | Logical (T/F) | Write partial energy for fluid 8 |
| `alpha_rho_e_wrt(9)` | Logical (T/F) | Write partial energy for fluid 9 |

---

## alpha_rho_wrt

*Partial density output flags*

**10 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `alpha_rho_wrt(1)` | Logical (T/F) | Write partial density for fluid 1 |
| `alpha_rho_wrt(10)` | Logical (T/F) | Write partial density for fluid 10 |
| `alpha_rho_wrt(2)` | Logical (T/F) | Write partial density for fluid 2 |
| `alpha_rho_wrt(3)` | Logical (T/F) | Write partial density for fluid 3 |
| `alpha_rho_wrt(4)` | Logical (T/F) | Write partial density for fluid 4 |
| `alpha_rho_wrt(5)` | Logical (T/F) | Write partial density for fluid 5 |
| `alpha_rho_wrt(6)` | Logical (T/F) | Write partial density for fluid 6 |
| `alpha_rho_wrt(7)` | Logical (T/F) | Write partial density for fluid 7 |
| `alpha_rho_wrt(8)` | Logical (T/F) | Write partial density for fluid 8 |
| `alpha_rho_wrt(9)` | Logical (T/F) | Write partial density for fluid 9 |

---

## alpha_wrt

*Volume fraction output flags*

**10 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `alpha_wrt(1)` | Logical (T/F) | Write volume fraction for fluid 1 |
| `alpha_wrt(10)` | Logical (T/F) | Write volume fraction for fluid 10 |
| `alpha_wrt(2)` | Logical (T/F) | Write volume fraction for fluid 2 |
| `alpha_wrt(3)` | Logical (T/F) | Write volume fraction for fluid 3 |
| `alpha_wrt(4)` | Logical (T/F) | Write volume fraction for fluid 4 |
| `alpha_wrt(5)` | Logical (T/F) | Write volume fraction for fluid 5 |
| `alpha_wrt(6)` | Logical (T/F) | Write volume fraction for fluid 6 |
| `alpha_wrt(7)` | Logical (T/F) | Write volume fraction for fluid 7 |
| `alpha_wrt(8)` | Logical (T/F) | Write volume fraction for fluid 8 |
| `alpha_wrt(9)` | Logical (T/F) | Write volume fraction for fluid 9 |

---

## fluid_rho

*Fluid reference densities*

**10 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `fluid_rho(1)` | Real | Reference density for fluid 1 |
| `fluid_rho(10)` | Real | Reference density for fluid 10 |
| `fluid_rho(2)` | Real | Reference density for fluid 2 |
| `fluid_rho(3)` | Real | Reference density for fluid 3 |
| `fluid_rho(4)` | Real | Reference density for fluid 4 |
| `fluid_rho(5)` | Real | Reference density for fluid 5 |
| `fluid_rho(6)` | Real | Reference density for fluid 6 |
| `fluid_rho(7)` | Real | Reference density for fluid 7 |
| `fluid_rho(8)` | Real | Reference density for fluid 8 |
| `fluid_rho(9)` | Real | Reference density for fluid 9 |

---

## kappa_wrt

*Curvature output flags*

**10 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `kappa_wrt(1)` | Logical (T/F) | Write curvature for fluid 1 |
| `kappa_wrt(10)` | Logical (T/F) | Write curvature for fluid 10 |
| `kappa_wrt(2)` | Logical (T/F) | Write curvature for fluid 2 |
| `kappa_wrt(3)` | Logical (T/F) | Write curvature for fluid 3 |
| `kappa_wrt(4)` | Logical (T/F) | Write curvature for fluid 4 |
| `kappa_wrt(5)` | Logical (T/F) | Write curvature for fluid 5 |
| `kappa_wrt(6)` | Logical (T/F) | Write curvature for fluid 6 |
| `kappa_wrt(7)` | Logical (T/F) | Write curvature for fluid 7 |
| `kappa_wrt(8)` | Logical (T/F) | Write curvature for fluid 8 |
| `kappa_wrt(9)` | Logical (T/F) | Write curvature for fluid 9 |

---

## schlieren_alpha

*Numerical schlieren coefficients*

**10 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `schlieren_alpha(1)` | Real | Schlieren coefficient for fluid 1 |
| `schlieren_alpha(10)` | Real | Schlieren coefficient for fluid 10 |
| `schlieren_alpha(2)` | Real | Schlieren coefficient for fluid 2 |
| `schlieren_alpha(3)` | Real | Schlieren coefficient for fluid 3 |
| `schlieren_alpha(4)` | Real | Schlieren coefficient for fluid 4 |
| `schlieren_alpha(5)` | Real | Schlieren coefficient for fluid 5 |
| `schlieren_alpha(6)` | Real | Schlieren coefficient for fluid 6 |
| `schlieren_alpha(7)` | Real | Schlieren coefficient for fluid 7 |
| `schlieren_alpha(8)` | Real | Schlieren coefficient for fluid 8 |
| `schlieren_alpha(9)` | Real | Schlieren coefficient for fluid 9 |

---

## chem_params

*Chemistry model parameters*

**4 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `chem_params%%diffusion` | Logical (T/F) | Chemistry parameter: diffusion |
| `chem_params%%gamma_method` | Integer | Chemistry parameter: gamma_method |
| `chem_params%%reactions` | Logical (T/F) | Chemistry parameter: reactions |
| `chem_params%%transport_model` | Integer | Chemistry parameter: transport_model |

---

## flux_wrt

*Flux output flags*

**3 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `flux_wrt(1)` | Logical (T/F) | Write flux component 1 |
| `flux_wrt(2)` | Logical (T/F) | Write flux component 2 |
| `flux_wrt(3)` | Logical (T/F) | Write flux component 3 |

---

## mom_wrt

*Momentum output flags*

**3 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `mom_wrt(1)` | Logical (T/F) | Write momentum component 1 |
| `mom_wrt(2)` | Logical (T/F) | Write momentum component 2 |
| `mom_wrt(3)` | Logical (T/F) | Write momentum component 3 |

---

## omega_wrt

*Vorticity output flags*

**3 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `omega_wrt(1)` | Logical (T/F) | Write vorticity component 1 |
| `omega_wrt(2)` | Logical (T/F) | Write vorticity component 2 |
| `omega_wrt(3)` | Logical (T/F) | Write vorticity component 3 |

---

## vel_wrt

*Velocity output flags*

**3 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `vel_wrt(1)` | Logical (T/F) | Write velocity component 1 |
| `vel_wrt(2)` | Logical (T/F) | Write velocity component 2 |
| `vel_wrt(3)` | Logical (T/F) | Write velocity component 3 |

---

## x_domain

*X-direction domain parameters*

**2 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `x_domain%%beg` | Real | Beginning of the x-direction domain |
| `x_domain%%end` | Real | End of the x-direction domain |

---

## x_output

*X-direction output region*

**2 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `x_output%%beg` | Real | X-direction output beg parameter |
| `x_output%%end` | Real | X-direction output end parameter |

---

## y_domain

*Y-direction domain parameters*

**2 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `y_domain%%beg` | Real | Beginning of the y-direction domain |
| `y_domain%%end` | Real | End of the y-direction domain |

---

## y_output

*Y-direction output region*

**2 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `y_output%%beg` | Real | Y-direction output beg parameter |
| `y_output%%end` | Real | Y-direction output end parameter |

---

## z_domain

*Z-direction domain parameters*

**2 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `z_domain%%beg` | Real | Beginning of the z-direction domain |
| `z_domain%%end` | Real | End of the z-direction domain |

---

## z_output

*Z-direction output region*

**2 parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `z_output%%beg` | Real | Z-direction output beg parameter |
| `z_output%%end` | Real | Z-direction output end parameter |

---

## general

*Core simulation parameters (grid, time, model, etc.)*

**312 parameters**

### Patterns

| Pattern | Example | Description |
|---------|---------|-------------|
| `Bx0` | `Bx0` | Background magnetic field in x-direction |
| `Ca` | `Ca` |  |
| `E_wrt` | `E_wrt` | Write energy field |
| `R0ref` | `R0ref` | Reference bubble radius |
| `Re_inv` | `Re_inv` |  |
| `Web` | `Web` |  |
| `a_x` | `a_x` | Rate of grid stretching in the x-direction |
| `a_y` | `a_y` | Rate of grid stretching in the y-direction |
| `a_z` | `a_z` | Rate of grid stretching in the z-direction |
| `acoustic_source` | `acoustic_source` | Enable acoustic source terms |
| `adap_dt` | `adap_dt` | Enable adaptive time stepping |
| `adap_dt_max_iters` | `adap_dt_max_iters` | Maximum iterations for adaptive time stepping |
| `adap_dt_tol` | `adap_dt_tol` | Tolerance for adaptive time stepping |
| `adv_n` | `adv_n` | Enable advection of number density |
| `alf_factor` | `alf_factor` | Artificial viscosity factor |
| `alpha_bar` | `alpha_bar` | Average volume fraction |
| `alpha_rho_wrt` | `alpha_rho_wrt` | Write partial density field |
| `alpha_wrt` | `alpha_wrt` | Write volume fraction field |
| `alt_soundspeed` | `alt_soundspeed` | Use alternative sound speed formulation |
| `avg_state` | `avg_state` | Average state for Riemann solver (1=Roe, 2=arithmetic) |
| `bf_x` | `bf_x` | Enable body force in x-direction |
| `bf_y` | `bf_y` | Enable body force in y-direction |
| `bf_z` | `bf_z` | Enable body force in z-direction |
| `bubble_model` | `bubble_model` | Bubble dynamics model (1=Gilmore, 2=Keller-Miksis, 3=Rayl... |
| `bubbles_euler` | `bubbles_euler` | Enable Euler-Euler bubble model |
| `bubbles_lagrange` | `bubbles_lagrange` | Enable Lagrangian bubble tracking |
| `c_wrt` | `c_wrt` | Write sound speed field |
| `cantera_file` | `cantera_file` | Cantera mechanism file for chemistry |
| `case_dir` | `case_dir` | Case directory path |
| `cf_wrt` | `cf_wrt` | Write color function field |
| `cfl_adap_dt` | `cfl_adap_dt` | Enable adaptive time stepping based on CFL |
| `cfl_const_dt` | `cfl_const_dt` | Use constant CFL for time stepping |
| `cfl_dt` | `cfl_dt` | Enable CFL-based time stepping |
| `cfl_max` | `cfl_max` | Maximum allowed CFL number |
| `cfl_target` | `cfl_target` | Target CFL number for adaptive time stepping |
| `chem_wrt_T` | `chem_wrt_T` | Write temperature field for chemistry |
| `chem_wrt_Y(N)` | `chem_wrt_Y(1)` | Write mass fraction of species 1 |
| `chemistry` | `chemistry` | Enable chemical reactions |
| `cons_vars_wrt` | `cons_vars_wrt` | Write conservative variables |
| `cont_damage` | `cont_damage` | Enable continuum damage model |
| `cont_damage_s` | `cont_damage_s` | Continuum damage shape parameter |
| `cyl_coord` | `cyl_coord` | Enable cylindrical coordinates (2D: axisymmetric, 3D: cyl... |
| `dist_type` | `dist_type` | Distribution type for polydisperse bubbles |
| `down_sample` | `down_sample` | Enable output downsampling |
| `dt` | `dt` | Time step size |
| `elliptic_smoothing` | `elliptic_smoothing` | Enable elliptic smoothing |
| `elliptic_smoothing_iters` | `elliptic_smoothing_iters` | Number of elliptic smoothing iterations |
| `fd_order` | `fd_order` | Finite difference order for gradients |
| `fft_wrt` | `fft_wrt` | Enable FFT output |
| `file_per_process` | `file_per_process` | Write separate file per MPI process |
| `fluid_rho` | `fluid_rho` | Reference fluid density |
| `flux_lim` | `flux_lim` | Flux limiter type |
| `flux_wrt` | `flux_wrt` | Write flux to output |
| `format` | `format` | Output format (1=Silo, 2=binary) |
| `g_x` | `g_x` | Body force parameter g in x-direction |
| `g_y` | `g_y` | Body force parameter g in y-direction |
| `g_z` | `g_z` | Body force parameter g in z-direction |
| `gamma_wrt` | `gamma_wrt` | Write gamma field |
| `heat_ratio_wrt` | `heat_ratio_wrt` | Write heat capacity ratio field |
| `hyperelasticity` | `hyperelasticity` | Enable hyperelastic model |
| `hypoelasticity` | `hypoelasticity` | Enable hypoelastic model |
| `ib` | `ib` | Enable immersed boundary method |
| `ic_beta` | `ic_beta` | Interface compression beta |
| `ic_eps` | `ic_eps` | Interface compression epsilon |
| `igr` | `igr` | Enable implicit gradient reconstruction |
| `igr_iter_solver` | `igr_iter_solver` | IGR iterative solver type |
| `igr_order` | `igr_order` | Implicit gradient reconstruction order |
| `igr_pres_lim` | `igr_pres_lim` | Enable IGR pressure limiting |
| `int_comp` | `int_comp` | Enable interface compression |
| `integral_wrt` | `integral_wrt` | Write integral data |
| `k_x` | `k_x` | Body force parameter k in x-direction |
| `k_y` | `k_y` | Body force parameter k in y-direction |
| `k_z` | `k_z` | Body force parameter k in z-direction |
| `kappa_wrt` | `kappa_wrt` | Write curvature field |
| `lag_betaC_wrt` | `lag_betaC_wrt` | Write Lagrangian betaC field |
| `lag_betaT_wrt` | `lag_betaT_wrt` | Write Lagrangian betaT field |
| `lag_db_wrt` | `lag_db_wrt` | Write Lagrangian db field |
| `lag_dphidt_wrt` | `lag_dphidt_wrt` | Write Lagrangian dphidt field |
| `lag_header` | `lag_header` | Enable Lagrangian output header |
| `lag_id_wrt` | `lag_id_wrt` | Write Lagrangian id field |
| `lag_mg_wrt` | `lag_mg_wrt` | Write Lagrangian mg field |
| `lag_mv_wrt` | `lag_mv_wrt` | Write Lagrangian mv field |
| `lag_pos_prev_wrt` | `lag_pos_prev_wrt` | Write Lagrangian pos_prev field |
| `lag_pos_wrt` | `lag_pos_wrt` | Write Lagrangian pos field |
| `lag_pres_wrt` | `lag_pres_wrt` | Write Lagrangian pres field |
| `lag_r0_wrt` | `lag_r0_wrt` | Write Lagrangian r0 field |
| `lag_rad_wrt` | `lag_rad_wrt` | Write Lagrangian rad field |
| `lag_rmax_wrt` | `lag_rmax_wrt` | Write Lagrangian rmax field |
| `lag_rmin_wrt` | `lag_rmin_wrt` | Write Lagrangian rmin field |
| `lag_rvel_wrt` | `lag_rvel_wrt` | Write Lagrangian rvel field |
| `lag_txt_wrt` | `lag_txt_wrt` | Write Lagrangian txt field |
| `lag_vel_wrt` | `lag_vel_wrt` | Write Lagrangian vel field |
| `liutex_wrt` | `liutex_wrt` | Write Liutex vortex field |
| `loops_x` | `loops_x` | Number of times to apply grid stretching in x |
| `loops_y` | `loops_y` | Number of times to apply grid stretching in y |
| `loops_z` | `loops_z` | Number of times to apply grid stretching in z |
| `low_Mach` | `low_Mach` | Low Mach number correction |
| `m` | `m` | Number of grid cells in the x-direction |
| `mapped_weno` | `mapped_weno` | Enable mapped WENO scheme |
| `mhd` | `mhd` | Enable magnetohydrodynamics |
| `mixlayer_domain` | `mixlayer_domain` | Mixing layer domain size |
| `mixlayer_perturb` | `mixlayer_perturb` | Enable mixing layer perturbation |
| `mixlayer_perturb_k0` | `mixlayer_perturb_k0` | Base wavenumber for mixing layer perturbation |
| `mixlayer_perturb_nk` | `mixlayer_perturb_nk` | Number of perturbation modes for mixing layer |
| `mixlayer_vel_coef` | `mixlayer_vel_coef` | Velocity coefficient for mixing layer |
| `mixlayer_vel_profile` | `mixlayer_vel_profile` | Enable mixing layer velocity profile |
| `mixture_err` | `mixture_err` | Enable mixture error checking |
| `model_eqns` | `model_eqns` | Model equations (1=gamma-law, 2=5-eq, 3=6-eq, 4=4-eq) |
| `mom_wrt` | `mom_wrt` | Write mom to output |
| `mp_weno` | `mp_weno` | Enable monotonicity-preserving WENO |
| `mpp_lim` | `mpp_lim` | Enable mixture pressure positivity limiter |
| `muscl_lim` | `muscl_lim` | MUSCL limiter type |
| `muscl_order` | `muscl_order` | Order of MUSCL reconstruction |
| `n` | `n` | Number of grid cells in the y-direction |
| `n_start` | `n_start` | Starting time step index |
| `n_start_old` | `n_start_old` | Starting index from previous simulation |
| `nb` | `nb` | Number of bubble bins for polydisperse model |
| `null_weights` | `null_weights` | Allow null WENO weights |
| `num_bc_patches` | `num_bc_patches` | Number of boundary condition patches |
| `num_fluids` | `num_fluids` | Number of fluid components |
| `num_ibs` | `num_ibs` | Number of immersed boundary patches |
| `num_igr_iters` | `num_igr_iters` | Number of IGR iterations |
| `num_igr_warm_start_iters` | `num_igr_warm_start_iters` | Number of IGR warm-start iterations |
| `num_integrals` | `num_integrals` | Number of integral regions |
| `num_patches` | `num_patches` | Number of initial condition patches |
| `num_probes` | `num_probes` | Number of probe points |
| `num_source` | `num_source` | Number of acoustic sources |
| `nv_uvm_igr_temps_on_gpu` | `nv_uvm_igr_temps_on_gpu` | Store IGR temporaries on GPU |
| `nv_uvm_out_of_core` | `nv_uvm_out_of_core` | Enable NVIDIA UVM out-of-core |
| `nv_uvm_pref_gpu` | `nv_uvm_pref_gpu` | Prefer GPU for NVIDIA UVM |
| `old_grid` | `old_grid` | Use grid from previous simulation |
| `old_ic` | `old_ic` | Use initial conditions from previous simulation |
| `omega_wrt` | `omega_wrt` | Write vorticity field |
| `output_partial_domain` | `output_partial_domain` | Enable partial domain output |
| `p` | `p` | Number of grid cells in the z-direction |
| `p_x` | `p_x` | Body force parameter p in x-direction |
| `p_y` | `p_y` | Body force parameter p in y-direction |
| `p_z` | `p_z` | Body force parameter p in z-direction |
| `palpha_eps` | `palpha_eps` | Volume fraction epsilon for pressure relaxation |
| `parallel_io` | `parallel_io` | Enable parallel I/O |
| `perturb_flow` | `perturb_flow` | Enable flow perturbation |
| `perturb_flow_fluid` | `perturb_flow_fluid` | Fluid index for flow perturbation |
| `perturb_flow_mag` | `perturb_flow_mag` | Magnitude of flow perturbation |
| `perturb_sph` | `perturb_sph` | Enable spherical perturbation |
| `perturb_sph_fluid` | `perturb_sph_fluid` | Fluid index for spherical perturbation |
| `pi_fac` | `pi_fac` | Pi infinity factor |
| `pi_inf_wrt` | `pi_inf_wrt` | Write pi_inf field |
| `poly_sigma` | `poly_sigma` | Polydisperse distribution standard deviation |
| `polydisperse` | `polydisperse` | Enable polydisperse bubble distribution |
| `polytropic` | `polytropic` | Enable polytropic gas behavior for bubbles |
| `powell` | `powell` | Enable Powell source terms for MHD |
| `pre_stress` | `pre_stress` | Enable pre-stress initialization |
| `precision` | `precision` | Output precision (1=single, 2=double) |
| `pref` | `pref` | Reference pressure |
| `pres_inf_wrt` | `pres_inf_wrt` | Write reference pressure field |
| `pres_wrt` | `pres_wrt` | Write pressure field |
| `prim_vars_wrt` | `prim_vars_wrt` | Write primitive variables |
| `probe_wrt` | `probe_wrt` | Write probe data |
| `ptgalpha_eps` | `ptgalpha_eps` | Volume fraction epsilon for PTG relaxation |
| `qbmm` | `qbmm` | Enable quadrature-based moment method |
| `qm_wrt` | `qm_wrt` | Write Q-criterion field |
| `rdma_mpi` | `rdma_mpi` | Enable RDMA for MPI communication (GPUs) |
| `recon_type` | `recon_type` | Reconstruction type (1=WENO, 2=MUSCL) |
| `relativity` | `relativity` | Enable special relativity |
| `relax` | `relax` | Enable relaxation terms |
| `relax_model` | `relax_model` | Relaxation model type |
| `rhoRV` | `rhoRV` | Bubble radius-velocity correlation |
| `rho_wrt` | `rho_wrt` | Write density field |
| `rhoref` | `rhoref` | Reference density |
| `riemann_solver` | `riemann_solver` | Riemann solver (1=HLL, 2=HLLC, 3=exact) |
| `run_time_info` | `run_time_info` | Output run-time information |
| `schlieren_alpha` | `schlieren_alpha` | Schlieren alpha coefficient |
| `schlieren_wrt` | `schlieren_wrt` | Write schlieren images |
| `sigR` | `sigR` | Bubble radius standard deviation |
| `sigV` | `sigV` | Bubble velocity standard deviation |
| `sigma` | `sigma` | Surface tension coefficient |
| `sim_data` | `sim_data` | Enable simulation data output |
| `simplex_perturb` | `simplex_perturb` | Enable simplex noise perturbation |
| `stretch_x` | `stretch_x` | Enable grid stretching in the x-direction |
| `stretch_y` | `stretch_y` | Enable grid stretching in the y-direction |
| `stretch_z` | `stretch_z` | Enable grid stretching in the z-direction |
| `surface_tension` | `surface_tension` | Enable surface tension effects |
| `t_save` | `t_save` | Time interval for saving data |
| `t_step_old` | `t_step_old` | Time step to restart from |
| `t_step_print` | `t_step_print` | Time step interval for printing info |
| `t_step_save` | `t_step_save` | Time step interval for saving data |
| `t_step_start` | `t_step_start` | Starting time step index |
| `t_step_stop` | `t_step_stop` | Ending time step index |
| `t_stop` | `t_stop` | Simulation stop time |
| `t_tol` | `t_tol` | Time tolerance |
| `tau_star` | `tau_star` | Non-dimensional relaxation time |
| `teno` | `teno` | Enable TENO scheme |
| `teno_CT` | `teno_CT` | TENO cutoff parameter |
| `thermal` | `thermal` | Thermal model selection |
| `time_stepper` | `time_stepper` | Time integration scheme (1=Euler, 2=TVD-RK2, 3=TVD-RK3) |
| `vel_wrt` | `vel_wrt` | Write velocity field |
| `viscous` | `viscous` | Enable viscous effects |
| `w_x` | `w_x` | Body force parameter w in x-direction |
| `w_y` | `w_y` | Body force parameter w in y-direction |
| `w_z` | `w_z` | Body force parameter w in z-direction |
| `wave_speeds` | `wave_speeds` | Wave speed estimates (1=direct, 2=pressure) |
| `weno_Re_flux` | `weno_Re_flux` | Enable WENO for viscous fluxes |
| `weno_avg` | `weno_avg` | Enable WENO averaging |
| `weno_eps` | `weno_eps` | WENO epsilon parameter for smoothness |
| `weno_order` | `weno_order` | Order of WENO reconstruction (1, 3, 5, or 7) |
| `wenoz` | `wenoz` | Enable WENO-Z scheme |
| `wenoz_q` | `wenoz_q` | WENO-Z power parameter |
| `x_a` | `x_a` | Start of stretching in negative x-direction |
| `x_b` | `x_b` | Start of stretching in positive x-direction |
| `y_a` | `y_a` | Start of stretching in negative y-direction |
| `y_b` | `y_b` | Start of stretching in positive y-direction |
| `z_a` | `z_a` | Start of stretching in negative z-direction |
| `z_b` | `z_b` | Start of stretching in positive z-direction |

---

## Command Line Reference

Search parameters using the CLI:

```bash
# Search for parameters
./mfc.sh params weno

# Show parameter descriptions
./mfc.sh params weno -d

# List all families
./mfc.sh params -f

# Filter by type
./mfc.sh params -t real weno
```
