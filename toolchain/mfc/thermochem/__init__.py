"""MFC-owned, Fortran-only thermochemistry generation."""

__all__ = ["generate_fortran"]


def generate_fortran(solution, module_name="m_thermochem"):
    """Load the expression machinery only when generating Fortran."""
    from .fortran import generate_fortran as generate

    return generate(solution, module_name)
