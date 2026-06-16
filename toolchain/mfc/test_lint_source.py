"""Tests for the manual registry-bound broadcast lint in lint_source.py."""

from mfc.lint_source import _extract_bcast_roots, check_manual_registry_bcasts

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
        "            & 'pref']",
        "            call MPI_BCAST(${VAR}$" + BCAST_TAIL,
        "        #:endfor",
    ]
    assert _extract_bcast_roots(lines) == [(3, "weno_eps"), (3, "teno_CT"), (3, "pref")]


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


def _write_proxy(tmp_path, target_dir: str, body: str):
    proxy_dir = tmp_path / "src" / target_dir
    proxy_dir.mkdir(parents=True)
    (proxy_dir / "m_mpi_proxy.fpp").write_text(body, encoding="utf-8")


def test_manual_broadcast_of_registry_scalar_is_flagged(tmp_path):
    _write_proxy(tmp_path, "simulation", f"        call MPI_BCAST(rhoref{BCAST_TAIL}\n")

    errors = check_manual_registry_bcasts(tmp_path)
    assert len(errors) == 1
    assert "manual MPI_BCAST of registry-bound scalar 'rhoref'" in errors[0]
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
