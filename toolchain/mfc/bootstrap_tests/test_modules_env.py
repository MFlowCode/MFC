"""Tests for the KEY=value loader in toolchain/bootstrap/modules.sh.

`. ./mfc.sh load` exports the assignments written on `toolchain/modules` lines.
Those lines carry two shapes that pull in opposite directions:

  * several assignments on one line -- `CC=nvc CXX=nvc++ FC=nvfortran`
  * a single value that itself contains spaces -- linker/compiler flags

The loader has to keep both intact. The tests below drive the real shell
function out of modules.sh rather than a copy of it, so they fail if the
implementation regresses.
"""

import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

MODULES_SH = Path(__file__).resolve().parents[2] / "bootstrap" / "modules.sh"

pytestmark = pytest.mark.skipif(shutil.which("bash") is None, reason="bash is required to exercise modules.sh")


def _extract_function(name):
    """Return the source text of a shell function defined in modules.sh.

    modules.sh is meant to be `source`d with arguments and runs an interactive
    loader at import time, so the function under test is lifted out on its own.
    """
    text = MODULES_SH.read_text()
    match = re.search(rf"^{name}\(\) \{{.*?^\}}", text, re.MULTILINE | re.DOTALL)
    assert match, f"{name}() not found in {MODULES_SH}"
    return match.group(0)


def _load(entry, preset_env=None):
    """Run one `toolchain/modules` entry through the loader.

    Returns the resulting environment as a dict, so a test can assert on the
    variables the entry was supposed to set.
    """
    preamble = "\n".join(f"export {k}={v}" for k, v in (preset_env or {}).items())
    script = textwrap.dedent(
        """\
        log() {{ :; }}          # modules.sh logging stub
        {preamble}
        {function_src}
        __export_assignments {entry}
        env
        """
    ).format(preamble=preamble, function_src=_extract_function("__export_assignments"), entry=_quote(entry))

    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True)
    env = {}
    for line in proc.stdout.splitlines():
        key, sep, value = line.partition("=")
        if sep:
            env[key] = value
    return env


def _quote(value):
    """Single-quote a string for safe interpolation into the bash script."""
    return "'" + value.replace("'", "'\\''") + "'"


@pytest.mark.parametrize(
    "entry,expected",
    [
        # Every multi-assignment shape currently in toolchain/modules.
        (
            "MFC_CUDA_CC=70,75,80 NVHPC_CUDA_HOME=$CUDA_HOME CC=nvc CXX=nvc++ FC=nvfortran",
            {"MFC_CUDA_CC": "70,75,80", "NVHPC_CUDA_HOME": "/opt/cuda", "CC": "nvc", "CXX": "nvc++", "FC": "nvfortran"},
        ),
        ("CC=nvc CXX=nvc++ FC=nvfortran", {"CC": "nvc", "CXX": "nvc++", "FC": "nvfortran"}),
        ("MPICC=mpiicc MPICXX=mpiicpc MPIFC=mpiifort", {"MPICC": "mpiicc", "MPICXX": "mpiicpc", "MPIFC": "mpiifort"}),
        ('UCX_NET_DEVICES="mlx5_4:1,mlx5_7:1"', {"UCX_NET_DEVICES": "mlx5_4:1,mlx5_7:1"}),
        ('PYTHONPATH=""', {"PYTHONPATH": ""}),
        ("NVHPC_CUDA_HOME=$CUDA_HOME", {"NVHPC_CUDA_HOME": "/opt/cuda"}),
    ],
)
def test_existing_module_lines_still_load(entry, expected):
    """Assignment shapes already in toolchain/modules keep their old meaning."""
    env = _load(entry, preset_env={"CUDA_HOME": "/opt/cuda"})
    for key, value in expected.items():
        assert env.get(key) == value


def test_multi_word_value_survives_intact():
    """A value with spaces is exported whole, not truncated to its first word.

    Written unquoted, `CRAY_CCE_LLD_ARGS` used to keep only the first flag: the
    rest were treated as further names to export and the failure went to stderr,
    which `./mfc.sh load` routinely discards. The build then succeeds with a
    compiler workaround silently missing.
    """
    flags = "-plugin-opt=-mattr=-mai-insts -plugin-opt=-disable-promote-alloca-to-vector"
    env = _load(f"CRAY_CCE_LLD_ARGS={flags}")
    assert env["CRAY_CCE_LLD_ARGS"] == flags


def test_quoted_multi_word_value_survives_intact():
    """Quoting the value in toolchain/modules works too, and drops the quotes."""
    env = _load('CRAY_CCE_LLD_ARGS="-O2 -g"')
    assert env["CRAY_CCE_LLD_ARGS"] == "-O2 -g"


def test_multi_word_value_followed_by_another_assignment():
    """A space-carrying value ends where the next NAME= begins."""
    env = _load("CFLAGS=-O2 -march=native CC=gcc")
    assert env["CFLAGS"] == "-O2 -march=native"
    assert env["CC"] == "gcc"


def test_variable_reference_inside_a_multi_word_value():
    """Values are still expanded once, so they can reference earlier exports."""
    env = _load("LDFLAGS=-L$CUDA_HOME/lib64 -lcudart", preset_env={"CUDA_HOME": "/opt/cuda"})
    assert env["LDFLAGS"] == "-L/opt/cuda/lib64 -lcudart"


def test_glob_characters_in_a_value_are_not_expanded():
    """A value such as `-Wl,*` must not pick up filenames from the cwd."""
    env = _load("MYFLAGS=-Wl,* -O2")
    assert env["MYFLAGS"] == "-Wl,* -O2"


def test_every_assignment_line_in_toolchain_modules_loads_cleanly():
    """No line shipped in toolchain/modules may error out or lose a value."""
    modules_file = MODULES_SH.resolve().parents[1] / "modules"
    entries = []
    for line in modules_file.read_text().splitlines():
        match = re.match(r"^[a-z][a-z0-9]*-(?:all|cpu|gpu)\s+(.*)$", line)
        if match and "=" in match.group(1):
            entries.append(match.group(1))

    assert entries, "expected at least one KEY=value line in toolchain/modules"

    for entry in entries:
        env = _load(entry, preset_env={"CUDA_HOME": "/opt/cuda", "OLCF_AFAR_ROOT": "/sw/afar"})
        for name in re.findall(r"(?:^|\s)([A-Za-z_][A-Za-z0-9_]*)=", entry):
            assert name in env, f"{name} was not exported by: {entry}"
