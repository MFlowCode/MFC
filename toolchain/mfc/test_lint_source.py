"""Tests for the manual registry-bound broadcast lint in lint_source.py."""

from mfc.lint_source import (
    _extract_bcast_roots,
    check_device_routine_element_args,
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


_LOOPED = """    subroutine s_curve(rho, i, p)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho
        integer, intent(in) :: i
        real(wp), intent(out) :: p
        integer :: it
        $:GPU_LOOP(parallelism='[seq]')
        do it = 1, 8
            p = p + rho
        end do
    end subroutine s_curve
    subroutine s_wrap(rho, i, p)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho
        integer, intent(in) :: i
        real(wp), intent(out) :: p
        call s_curve(rho, i, p)
    end subroutine s_wrap
    subroutine s_plain(rho, i, p)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho
        integer, intent(in) :: i
        real(wp), intent(out) :: p
        p = rho
    end subroutine s_plain
"""


def _KERNEL(body: str) -> str:
    return "        $:GPU_PARALLEL_LOOP(collapse=3)\n        " + body + "\n        $:END_GPU_PARALLEL_LOOP()\n"


def test_host_call_sites_are_not_flagged(tmp_path):
    _write_src(tmp_path, "simulation/m_x.fpp", _LOOPED + "        call s_curve(q(1)%sf(j, k, l), 1, out(k, l, q))\n")
    assert check_device_routine_element_args(tmp_path) == []


def test_element_into_looped_device_routine_is_flagged(tmp_path):
    _write_src(tmp_path, "simulation/m_x.fpp", _LOOPED + _KERNEL("call s_curve(q(1)%sf(j, k, l), 1, out(k, l, q))"))
    errors = check_device_routine_element_args(tmp_path)
    assert len(errors) == 2
    assert "q(1)%sf(j, k, l)" in errors[0] and "out(k, l, q)" in errors[1]


def test_element_reaches_the_loop_through_a_caller(tmp_path):
    _write_src(tmp_path, "simulation/m_x.fpp", _LOOPED + _KERNEL("call s_wrap(pres, 1, blkmod(k, &\n            & l, q))"))
    errors = check_device_routine_element_args(tmp_path)
    assert len(errors) == 1 and "s_wrap" in errors[0]


def test_scalars_expressions_and_loopless_routines_pass(tmp_path):
    _write_src(
        tmp_path,
        "simulation/m_x.fpp",
        _LOOPED
        + _KERNEL("call s_curve(alpha_rho(i)/max(alpha(i), sgm_eps), i, p_i)\n        call s_plain(q(1)%sf(j, k, l), 1, out(k, l, q))\n        call s_curve(real(q(1)%sf(j, k, l), wp), 1, p_i)"),
    )
    assert check_device_routine_element_args(tmp_path) == []


def test_loop_inside_a_device_function_counts_and_propagates(tmp_path):
    src = """    function f_looped(x, i) result(y)
        $:GPU_ROUTINE(function_name='f_looped', parallelism='[seq]')
        real(wp), intent(in) :: x
        integer, intent(in) :: i
        real(wp) :: y
        integer :: it
        y = x
        $:GPU_LOOP(parallelism='[seq]')
        do it = 1, 8
            y = y + 1._wp
        end do
    end function f_looped
    subroutine s_via_function(x, i, y)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: x
        integer, intent(in) :: i
        real(wp), intent(out) :: y
        y = f_looped(x, i)
    end subroutine s_via_function
"""
    calls = "out(k, l, q) = f_looped(q(1)%sf(k, l, q), 1)\n        call s_via_function(q(1)%sf(k, l, q), 1, tmp)\n        tmp = f_looped(p_scalar, 1)"
    _write_src(tmp_path, "simulation/m_x.fpp", src + _KERNEL(calls))
    errors = check_device_routine_element_args(tmp_path)
    assert [e.split("`")[3] for e in errors] == ["f_looped", "s_via_function"]


def test_constructor_commas_and_unprefixed_functions_and_contained_scoping(tmp_path):
    src = """    function g_looped(x) result(y)
        $:GPU_ROUTINE(function_name='g_looped', parallelism='[seq]')
        real(wp), intent(in) :: x
        real(wp) :: y
        integer :: it
        y = x
        $:GPU_LOOP(parallelism='[seq]')
        do it = 1, 8
            y = y + 1._wp
        end do
    end function g_looped
    subroutine s_outer(a, b)
        real(wp), intent(in) :: a
        real(wp), intent(out) :: b
        b = a
    contains
        subroutine s_inner(x, y)
            $:GPU_ROUTINE(parallelism='[seq]')
            real(wp), intent(in) :: x
            real(wp), intent(out) :: y
            y = g_looped(x)
        end subroutine s_inner
    end subroutine s_outer
"""
    calls = "b = g_looped(q(1)%sf(k, l, q))\\n        c = g_looped(sum([v(1), v(2)]))\\n        call s_outer(q(1)%sf(k, l, q), tmp)"
    _write_src(tmp_path, "simulation/m_x.fpp", src + _KERNEL(calls))
    errors = check_device_routine_element_args(tmp_path)
    # the unprefixed function is found by name; the array constructor is not split into a fake element;
    # s_outer is not tainted by the loop that only its contained s_inner reaches (and is not a device routine)
    assert [e.split("`")[1] for e in errors] == ["q(1)%sf(k, l, q)"]
