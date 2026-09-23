"""MFC-owned, Fortran-only thermochemistry generation."""

__all__ = ["generate_fortran", "generate_surface_fortran"]


def generate_fortran(solution, module_name="m_thermochem", scalar_type="real(dp)", offload=None):
    """Load the expression machinery only when generating Fortran."""
    from .fortran import generate_fortran as generate

    return generate(solution, module_name, scalar_type, offload)


def generate_surface_fortran(gas, surface=None, module_name="m_surface_thermochem", scalar_type="real(dp)", offload=None):
    """Load the expression machinery only when generating the surface module."""
    from .surface import generate_surface_fortran as generate

    return generate(gas, surface, module_name, scalar_type, offload)
