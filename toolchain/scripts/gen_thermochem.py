#!/usr/bin/env python3
"""
Generate m_thermochem.f90 module using Pyrometheus and Cantera.

This script is called by CMake when MFC_CHEMISTRY=ON to generate
the Fortran thermochemistry module from a Cantera mechanism file.
"""
import argparse
import sys

try:
    import cantera as ct
except ImportError as exc:
    print("ERROR: cantera Python package is required for chemistry.", file=sys.stderr)
    print("Install it with: pip install cantera", file=sys.stderr)
    raise ImportError("cantera is required for chemistry code generation") from exc

try:
    import pyrometheus as pyro
except ImportError as exc:
    print("ERROR: pyrometheus Python package is required for chemistry.", file=sys.stderr)
    print("Install it with: pip install pyrometheus", file=sys.stderr)
    raise ImportError("pyrometheus is required for chemistry code generation") from exc

def main():
    """Generate m_thermochem.f90 using Pyrometheus and Cantera."""
    parser = argparse.ArgumentParser(
        description="Generate Fortran thermochemistry module via Pyrometheus"
    )
    parser.add_argument("--module", required=True, help="Module name (e.g., m_thermochem)")
    parser.add_argument("--out", required=True, help="Output Fortran file path")
    parser.add_argument(
        "--mech",
        default="h2o2.yaml",
        help="Cantera mechanism file (YAML). Defaults to h2o2.yaml"
    )
    args = parser.parse_args()

    # Load the Cantera solution
    # If no mechanism file is specified or it's empty, use Cantera's default h2o2.yaml
    mech_file = args.mech if args.mech else "h2o2.yaml"

    try:
        sol = ct.Solution(mech_file)
    except Exception as exc:
        print(f"ERROR: Failed to load Cantera mechanism '{mech_file}': {exc}", file=sys.stderr)
        sys.exit(1)

    # Generate Fortran code using Pyrometheus
    try:
        code = pyro.FortranCodeGenerator().generate(args.module, sol)
    except Exception as exc:
        print(f"ERROR: Failed to generate Fortran code: {exc}", file=sys.stderr)
        sys.exit(1)

    # Write to output file
    try:
        with open(args.out, 'w', encoding='utf-8') as f:
            f.write(code)
        print(f"Successfully generated {args.out} from {mech_file}")
    except Exception as exc:
        print(f"ERROR: Failed to write output file '{args.out}': {exc}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
