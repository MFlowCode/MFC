import os
import shutil
import stat
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

from mfc.common import delete_directory, get_py_program_output


def test_case_files_run_under_the_toolchain_interpreter():
    """A case file must run under sys.executable, not under whatever python3 PATH finds.

    coverage-health.yml invokes build/venv/bin/python3 directly instead of activating the
    venv, so PATH's python3 was the system interpreter and every chemistry case failed with
    ModuleNotFoundError: No module named 'cantera' -- in a job whose venv had cantera. The
    case file reports the interpreter that ran it; it must be this one.
    """
    with tempfile.TemporaryDirectory() as d:
        case = Path(d) / "case.py"
        case.write_text("import sys\nprint(sys.executable)\n")
        out, err = get_py_program_output(str(case), [])
    assert err == 0
    assert out.strip() == sys.executable


def test_a_hostile_python3_on_path_is_not_consulted():
    """Put a python3 on PATH that refuses to run anything; the case file must still load.

    This is the failure the health job hit, made deterministic: PATH's python3 is a
    different, wrong interpreter. Asserting on sys.executable alone would still pass if
    someone reintroduced a PATH lookup on a machine where the two happen to coincide.
    """
    with tempfile.TemporaryDirectory() as d:
        bindir = Path(d) / "bin"
        bindir.mkdir()
        fake = bindir / "python3"
        fake.write_text("#!/bin/sh\necho 'wrong interpreter' >&2\nexit 3\n")
        fake.chmod(fake.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
        case = Path(d) / "case.py"
        case.write_text("print('loaded')\n")
        with patch.dict(os.environ, {"PATH": f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}"}):
            out, err = get_py_program_output(str(case), [])
    assert err == 0, "case file ran under PATH's python3 instead of sys.executable"
    assert out.strip() == "loaded"


def test_delete_directory_tolerates_enoent_from_a_concurrent_unlink():
    """rmtree can raise ENOENT mid-sweep for an entry already gone (#1968); finish the sweep instead."""
    real_rmtree = shutil.rmtree

    def flaky_rmtree(path, *args, **kwargs):
        if not kwargs.get("ignore_errors"):
            os.rmdir(os.path.join(path, "p0"))
            raise FileNotFoundError(2, "No such file or directory", "p0")
        return real_rmtree(path, *args, **kwargs)

    with tempfile.TemporaryDirectory() as d:
        target = Path(d) / "silo_hdf5"
        (target / "p0").mkdir(parents=True)
        (target / "other").mkdir()
        with patch("mfc.common.shutil.rmtree", flaky_rmtree):
            delete_directory(str(target))
        assert not target.exists()
