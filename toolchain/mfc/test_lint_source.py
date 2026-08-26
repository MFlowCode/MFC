"""Tests for the manual registry-bound broadcast lint in lint_source.py."""

from mfc.lint_source import (
    _extract_bcast_roots,
    check_double_precision,
    check_integer_wp,
    check_manual_registry_bcasts,
)

BCAST_TAIL = ", 1, mpi_p, 0, MPI_COMM_WORLD, ierr)"


def test_extract_direct_call_root():
    lines = [f"        call MPI_BCAST(m_glb{BCAST_TAIL}"]
    assert _extract_bcast_roots(lines) == [(1, "m_glb")]


def test_extract_indexed_argument_keeps_root():
    lines = [f"        call MPI_BCAST(phi_rn(1){BCAST_TAIL}"]
    assert _extract_bcast_roots(lines) == [(1, "phi_rn")]


def test_extract_fypp_list_with_continuation():
    lines = [
        "        #:for VAR in [ 'weno_eps', 'teno_CT', &",
        "            & 'poly_sigma']",
        "            call MPI_BCAST(${VAR}$" + BCAST_TAIL,
        "        #:endfor",
    ]
    assert _extract_bcast_roots(lines) == [(3, "weno_eps"), (3, "teno_CT"), (3, "poly_sigma")]


def test_struct_members_and_loop_indices_skipped():
    lines = [
        "        #:for VAR in [ 'bc_x%beg', 'bc_x%end']",
        "            call MPI_BCAST(${VAR}$" + BCAST_TAIL,
        "        #:endfor",
        "        #:for DIM in ['x', 'y', 'z']",
        "            #:for DIR in [1, 2, 3]",
        "                call MPI_BCAST(bc_${DIM}$%vb${DIR}$" + BCAST_TAIL,
        "            #:endfor",
        "        #:endfor",
        f"        call MPI_BCAST(patch_ib(i)%geometry{BCAST_TAIL}",
        "        #:for VAR in ['c', 'p', 't', 'm']",
        "            call MPI_BCAST(ib_airfoil(i)%${VAR}$" + BCAST_TAIL,
        "        #:endfor",
    ]
    assert _extract_bcast_roots(lines) == []


def _write_src(tmp_path, rel: str, body: str):
    """Write a source file under tmp_path/src/<rel> for a lint check to scan."""
    path = tmp_path / "src" / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")


def _write_proxy(tmp_path, target_dir: str, body: str):
    _write_src(tmp_path, f"{target_dir}/m_mpi_proxy.fpp", body)


def test_double_precision_flags_signed_d_exponent(tmp_path):
    _write_src(tmp_path, "simulation/m_x.fpp", "        if (x > 5.0d-11) then\n")
    errors = check_double_precision(tmp_path)
    assert len(errors) == 1
    assert "5.0d-11" in errors[0]


def test_double_precision_flags_whole_mantissa(tmp_path):
    """A multi-digit mantissa is caught, and the message quotes the whole literal."""
    body = "\n".join(
        [
            "        p0 = 101325.d0",
            "        x = 1013.25d3",
            "        y = .5d0",
            "",
        ]
    )
    _write_src(tmp_path, "simulation/m_w.fpp", body)
    errors = check_double_precision(tmp_path)
    assert len(errors) == 3
    assert "101325.d0" in errors[0]
    assert "1013.25d3" in errors[1]
    assert ".5d0" in errors[2]


def test_double_precision_ignores_d_edit_descriptor(tmp_path):
    """A 'D' edit descriptor with a repeat count is not a literal."""
    _write_src(tmp_path, "simulation/m_v.fpp", "        write (*, '(1D12.4)') x\n")
    assert check_double_precision(tmp_path) == []


def test_double_precision_clean_cases(tmp_path):
    body = "\n".join(
        [
            "    integer, dimension(2) :: cart2d12_coords, cart2d13_coords",  # identifier, not a literal
            "        call MPI_CART_COORDS(comm, rank, 2, cart2d12_coords, ierr)",
            "        x = 5.0e-11_wp",  # correct working-precision literal
            "",
        ]
    )
    _write_src(tmp_path, "simulation/m_y.fpp", body)
    assert check_double_precision(tmp_path) == []


def test_integer_wp_flagged(tmp_path):
    _write_src(tmp_path, "simulation/m_z.fpp", "        integer(wp) :: i, j, k, l\n")
    errors = check_integer_wp(tmp_path)
    assert len(errors) == 1
    assert "integer(wp)" in errors[0]


def test_integer_wp_flags_explicit_kind_keyword(tmp_path):
    """integer(kind=wp) is the same defect, spelled the way real(kind=wp) is."""
    _write_src(tmp_path, "simulation/m_k.fpp", "        integer(kind=wp) :: i, j, k, l\n")
    errors = check_integer_wp(tmp_path)
    assert len(errors) == 1
    assert "integer(kind=wp)" in errors[0]


def test_integer_wp_clean(tmp_path):
    body = "        integer :: i\n        real(wp) :: x\n        integer(KIND=MPI_OFFSET_KIND) :: disp\n"
    _write_src(tmp_path, "simulation/m_ok.fpp", body)
    assert check_integer_wp(tmp_path) == []


def test_manual_broadcast_of_registry_scalar_is_flagged(tmp_path):
    _write_proxy(tmp_path, "simulation", f"        call MPI_BCAST(poly_sigma{BCAST_TAIL}\n")

    errors = check_manual_registry_bcasts(tmp_path)
    assert len(errors) == 1
    assert "manual MPI_BCAST of registry-bound scalar 'poly_sigma'" in errors[0]
    assert "m_mpi_proxy.fpp:1" in errors[0]


def test_manual_residue_is_clean(tmp_path):
    body = "\n".join(
        [
            f"        call MPI_BCAST(m_glb{BCAST_TAIL}",
            f"        call MPI_BCAST(cfl_dt{BCAST_TAIL}",
            "        #:for VAR in [ 'bc_x%beg', 'bc_x%end']",
            "            call MPI_BCAST(${VAR}$" + BCAST_TAIL,
            "        #:endfor",
            f"        call MPI_BCAST(patch_icpp(i)%geometry{BCAST_TAIL}",
            "",
        ]
    )
    _write_proxy(tmp_path, "simulation", body)

    assert check_manual_registry_bcasts(tmp_path) == []
