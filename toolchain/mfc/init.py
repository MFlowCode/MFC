"""MFC Case Template Generator - Create new case files from templates."""

import os
import shutil

from .printer import cons
from .common import MFC_EXAMPLE_DIRPATH, MFCException
from .state import ARG


# Built-in minimal templates
BUILTIN_TEMPLATES = {
    '1D_minimal': '''\
#!/usr/bin/env python3
"""
1D Minimal Case Template
------------------------
A minimal 1D shock tube case to get started with MFC.

Usage:
    ./mfc.sh run case.py
"""
import math
import json

# =============================================================================
# SIMULATION PARAMETERS - Modify these for your case
# =============================================================================

# Grid resolution
Nx = 399                    # Number of cells in x-direction

# Domain size
x_start = 0.0               # Domain start
x_end = 1.0                 # Domain end

# Time stepping
t_end = 0.1                 # End time
Nt = 1000                   # Number of time steps

# Initial conditions for left state (patch 1)
rho_L = 1.0                 # Density
vel_L = 0.0                 # Velocity
pres_L = 1.0                # Pressure

# Initial conditions for right state (patch 2)
rho_R = 0.125               # Density
vel_R = 0.0                 # Velocity
pres_R = 0.1                # Pressure

# Fluid properties
gamma = 1.4                 # Ratio of specific heats

# =============================================================================
# DERIVED QUANTITIES - Usually don't need to modify
# =============================================================================
dx = (x_end - x_start) / (Nx + 1)
dt = t_end / Nt

# =============================================================================
# CASE DICTIONARY - MFC configuration
# =============================================================================
print(json.dumps({
    # Logistics
    "run_time_info": "T",

    # Computational Domain
    "x_domain%beg": x_start,
    "x_domain%end": x_end,
    "m": Nx,
    "n": 0,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": max(1, Nt // 10),

    # Simulation Algorithm
    "num_patches": 2,
    "model_eqns": 2,            # 5-equation model
    "num_fluids": 1,
    "time_stepper": 3,          # TVD RK3
    "weno_order": 5,            # WENO5
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": 2,        # HLLC
    "wave_speeds": 1,
    "avg_state": 2,

    # Boundary Conditions (-3 = extrapolation)
    "bc_x%beg": -3,
    "bc_x%end": -3,

    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",

    # Patch 1: Left state
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": (x_start + x_end) / 4,
    "patch_icpp(1)%length_x": (x_end - x_start) / 2,
    "patch_icpp(1)%vel(1)": vel_L,
    "patch_icpp(1)%pres": pres_L,
    "patch_icpp(1)%alpha_rho(1)": rho_L,
    "patch_icpp(1)%alpha(1)": 1.0,

    # Patch 2: Right state
    "patch_icpp(2)%geometry": 1,
    "patch_icpp(2)%x_centroid": 3 * (x_start + x_end) / 4,
    "patch_icpp(2)%length_x": (x_end - x_start) / 2,
    "patch_icpp(2)%vel(1)": vel_R,
    "patch_icpp(2)%pres": pres_R,
    "patch_icpp(2)%alpha_rho(1)": rho_R,
    "patch_icpp(2)%alpha(1)": 1.0,

    # Fluid Properties
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}))
''',

    '2D_minimal': '''\
#!/usr/bin/env python3
"""
2D Minimal Case Template
------------------------
A minimal 2D case with a circular perturbation.

Usage:
    ./mfc.sh run case.py
"""
import math
import json

# =============================================================================
# SIMULATION PARAMETERS - Modify these for your case
# =============================================================================

# Grid resolution
Nx = 99                     # Cells in x-direction
Ny = 99                     # Cells in y-direction

# Domain size
x_start, x_end = 0.0, 1.0
y_start, y_end = 0.0, 1.0

# Time stepping
dt = 1.0e-6
Nt = 1000

# Background state
rho_bg = 1.0
vel_x_bg = 0.0
vel_y_bg = 0.0
pres_bg = 1.0e5

# Perturbation (circular region)
x_center = 0.5
y_center = 0.5
radius = 0.1
rho_pert = 2.0
pres_pert = 2.0e5

# Fluid properties
gamma = 1.4

# =============================================================================
# CASE DICTIONARY - MFC configuration
# =============================================================================
print(json.dumps({
    # Logistics
    "run_time_info": "T",

    # Computational Domain
    "x_domain%beg": x_start,
    "x_domain%end": x_end,
    "y_domain%beg": y_start,
    "y_domain%end": y_end,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": max(1, Nt // 10),

    # Simulation Algorithm
    "num_patches": 2,
    "model_eqns": 2,
    "num_fluids": 1,
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,

    # Boundary Conditions
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,

    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",

    # Patch 1: Background
    "patch_icpp(1)%geometry": 3,  # Rectangle
    "patch_icpp(1)%x_centroid": (x_start + x_end) / 2,
    "patch_icpp(1)%y_centroid": (y_start + y_end) / 2,
    "patch_icpp(1)%length_x": x_end - x_start,
    "patch_icpp(1)%length_y": y_end - y_start,
    "patch_icpp(1)%vel(1)": vel_x_bg,
    "patch_icpp(1)%vel(2)": vel_y_bg,
    "patch_icpp(1)%pres": pres_bg,
    "patch_icpp(1)%alpha_rho(1)": rho_bg,
    "patch_icpp(1)%alpha(1)": 1.0,

    # Patch 2: Circular perturbation
    "patch_icpp(2)%geometry": 2,  # Circle
    "patch_icpp(2)%x_centroid": x_center,
    "patch_icpp(2)%y_centroid": y_center,
    "patch_icpp(2)%radius": radius,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": vel_x_bg,
    "patch_icpp(2)%vel(2)": vel_y_bg,
    "patch_icpp(2)%pres": pres_pert,
    "patch_icpp(2)%alpha_rho(1)": rho_pert,
    "patch_icpp(2)%alpha(1)": 1.0,

    # Fluid Properties
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}))
''',

    '3D_minimal': '''\
#!/usr/bin/env python3
"""
3D Minimal Case Template
------------------------
A minimal 3D case with a spherical perturbation.

Usage:
    ./mfc.sh run case.py -N 2 -n 4  # Run on 2 nodes with 4 tasks each
"""
import math
import json

# =============================================================================
# SIMULATION PARAMETERS - Modify these for your case
# =============================================================================

# Grid resolution (keep low for testing, increase for production)
Nx = 49                     # Cells in x-direction
Ny = 49                     # Cells in y-direction
Nz = 49                     # Cells in z-direction

# Domain size
x_start, x_end = 0.0, 1.0
y_start, y_end = 0.0, 1.0
z_start, z_end = 0.0, 1.0

# Time stepping
dt = 1.0e-6
Nt = 100                    # Keep low for testing

# Background state
rho_bg = 1.0
pres_bg = 1.0e5

# Spherical perturbation
x_center = 0.5
y_center = 0.5
z_center = 0.5
radius = 0.1
rho_pert = 2.0
pres_pert = 2.0e5

# Fluid properties
gamma = 1.4

# =============================================================================
# CASE DICTIONARY - MFC configuration
# =============================================================================
print(json.dumps({
    # Logistics
    "run_time_info": "T",

    # Computational Domain
    "x_domain%beg": x_start,
    "x_domain%end": x_end,
    "y_domain%beg": y_start,
    "y_domain%end": y_end,
    "z_domain%beg": z_start,
    "z_domain%end": z_end,
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": max(1, Nt // 10),

    # Simulation Algorithm
    "num_patches": 2,
    "model_eqns": 2,
    "num_fluids": 1,
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,

    # Boundary Conditions
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "bc_z%beg": -3,
    "bc_z%end": -3,

    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",

    # Patch 1: Background (cube)
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": (x_start + x_end) / 2,
    "patch_icpp(1)%y_centroid": (y_start + y_end) / 2,
    "patch_icpp(1)%z_centroid": (z_start + z_end) / 2,
    "patch_icpp(1)%length_x": x_end - x_start,
    "patch_icpp(1)%length_y": y_end - y_start,
    "patch_icpp(1)%length_z": z_end - z_start,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": pres_bg,
    "patch_icpp(1)%alpha_rho(1)": rho_bg,
    "patch_icpp(1)%alpha(1)": 1.0,

    # Patch 2: Spherical perturbation
    "patch_icpp(2)%geometry": 8,  # Sphere
    "patch_icpp(2)%x_centroid": x_center,
    "patch_icpp(2)%y_centroid": y_center,
    "patch_icpp(2)%z_centroid": z_center,
    "patch_icpp(2)%radius": radius,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": pres_pert,
    "patch_icpp(2)%alpha_rho(1)": rho_pert,
    "patch_icpp(2)%alpha(1)": 1.0,

    # Fluid Properties
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}))
''',
}


def get_available_templates():
    """Get list of available templates (built-in + examples)."""
    templates = list(BUILTIN_TEMPLATES.keys())

    # Add examples as templates
    if os.path.isdir(MFC_EXAMPLE_DIRPATH):
        for name in sorted(os.listdir(MFC_EXAMPLE_DIRPATH)):
            example_path = os.path.join(MFC_EXAMPLE_DIRPATH, name)
            if os.path.isdir(example_path) and os.path.isfile(os.path.join(example_path, 'case.py')):
                templates.append(f"example:{name}")

    return templates


def list_templates():
    """Print available templates."""
    cons.print("[bold]Available Templates[/bold]\n")

    cons.print("  [bold cyan]Built-in Templates:[/bold cyan]")
    for name in sorted(BUILTIN_TEMPLATES.keys()):
        desc = {
            '1D_minimal': 'Minimal 1D shock tube case',
            '2D_minimal': 'Minimal 2D case with circular perturbation',
            '3D_minimal': 'Minimal 3D case with spherical perturbation',
        }.get(name, '')
        cons.print(f"    [green]{name:20s}[/green] {desc}")

    cons.print()
    cons.print("  [bold cyan]From Examples:[/bold cyan]")

    if os.path.isdir(MFC_EXAMPLE_DIRPATH):
        examples = []
        for name in sorted(os.listdir(MFC_EXAMPLE_DIRPATH)):
            example_path = os.path.join(MFC_EXAMPLE_DIRPATH, name)
            if os.path.isdir(example_path) and os.path.isfile(os.path.join(example_path, 'case.py')):
                examples.append(name)

        # Group by dimension
        for dim in ['0D', '1D', '2D', '3D']:
            dim_examples = [e for e in examples if e.startswith(dim)]
            if dim_examples:
                cons.print(f"    [dim]{dim}:[/dim] {', '.join(dim_examples[:5])}", end='')
                if len(dim_examples) > 5:
                    cons.print(f" [dim]... (+{len(dim_examples) - 5} more)[/dim]")
                else:
                    cons.print()

    cons.print()
    cons.print("  [bold]Usage:[/bold]")
    cons.print("    ./mfc.sh new my_case                        # Use default 1D template")
    cons.print("    ./mfc.sh new my_case --template 2D_minimal  # Use 2D template")
    cons.print("    ./mfc.sh new my_case --template example:1D_sodshocktube  # Copy from example")
    cons.print()


def create_case(name: str, template: str):
    """Create a new case from a template."""
    # Determine output directory
    output_dir = os.path.abspath(name)

    if os.path.exists(output_dir):
        raise MFCException(f"Directory already exists: {output_dir}")

    # Check if it's a built-in template
    if template in BUILTIN_TEMPLATES:
        os.makedirs(output_dir, exist_ok=True)
        case_path = os.path.join(output_dir, 'case.py')

        with open(case_path, 'w') as f:
            f.write(BUILTIN_TEMPLATES[template])

        os.chmod(case_path, 0o755)  # Make executable

        cons.print(f"[bold green]Created[/bold green] {output_dir}/")
        cons.print(f"  Using template: [cyan]{template}[/cyan]")
        cons.print()
        cons.print("  [bold]Next steps:[/bold]")
        cons.print(f"    1. Edit [cyan]{name}/case.py[/cyan] to configure your simulation")
        cons.print(f"    2. Run: [cyan]./mfc.sh run {name}/case.py[/cyan]")
        cons.print()

    # Check if it's an example template
    elif template.startswith('example:'):
        example_name = template[8:]  # Remove 'example:' prefix
        example_path = os.path.join(MFC_EXAMPLE_DIRPATH, example_name)

        if not os.path.isdir(example_path):
            raise MFCException(f"Example not found: {example_name}")

        # Copy the example directory
        shutil.copytree(example_path, output_dir)

        cons.print(f"[bold green]Created[/bold green] {output_dir}/")
        cons.print(f"  Copied from example: [cyan]{example_name}[/cyan]")
        cons.print()
        cons.print("  [bold]Next steps:[/bold]")
        cons.print(f"    1. Review and modify [cyan]{name}/case.py[/cyan]")
        cons.print(f"    2. Run: [cyan]./mfc.sh run {name}/case.py[/cyan]")
        cons.print()

    else:
        available = ', '.join(list(BUILTIN_TEMPLATES.keys())[:3])
        raise MFCException(
            f"Unknown template: {template}\n"
            f"Available built-in templates: {available}\n"
            f"Or use 'example:<name>' to copy from examples.\n"
            f"Run './mfc.sh new --list' to see all available templates."
        )


def init():
    """Main entry point for the init command."""
    if ARG("list"):
        list_templates()
        return

    name = ARG("name")
    template = ARG("template")

    if not name:
        # Show full help like ./mfc.sh new -h
        # pylint: disable=import-outside-toplevel
        import sys
        from .user_guide import print_command_help
        from .cli.commands import MFC_CLI_SCHEMA
        from .cli.argparse_gen import generate_parser
        from .state import MFCConfig

        print_command_help("new", show_argparse=False)
        _, subparser_map = generate_parser(MFC_CLI_SCHEMA, MFCConfig())
        subparser_map["new"].print_help()
        sys.stdout.flush()
        sys.stderr.write("\n./mfc.sh new: error: the following arguments are required: NAME\n")
        sys.exit(2)

    create_case(name, template)
