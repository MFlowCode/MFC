"""Tests for the toolchain/modules line reader in toolchain/bootstrap/modules.sh.

`. ./mfc.sh load` reads each `toolchain/modules` line twice: once to collect the
module names it should `module load`, and once to export the KEY=value
assignments. Those lines carry two shapes that pull in opposite directions:

  * several assignments on one line -- `CC=nvc CXX=nvc++ FC=nvfortran`
  * a single value that itself contains spaces -- linker/compiler flags

Both readers have to agree on where an assignment starts and ends, or a value's
words leak into the module list. The tests below drive the real shell functions
out of modules.sh rather than a copy of them, so they fail if the implementation
regresses.
"""

import re
import shlex
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


def _helpers(*names):
    """Return the source of the named functions plus the glob-guard helpers."""
    return "\n".join(_extract_function(n) for n in ("__noglob_push", "__noglob_pop", *names))


def _load(entry, preset_env=None):
    """Run one `toolchain/modules` entry through the exporter.

    Returns the resulting environment as a dict, so a test can assert on the
    variables the entry was supposed to set.
    """
    preamble = "\n".join(f"export {k}={shlex.quote(v)}" for k, v in (preset_env or {}).items())
    script = textwrap.dedent(
        """\
        log() {{ :; }}          # modules.sh logging stub
        {preamble}
        {function_src}
        __export_assignments {entry}
        env
        """
    ).format(preamble=preamble, function_src=_helpers("__export_assignments"), entry=_quote(entry))

    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True)
    env = {}
    for line in proc.stdout.splitlines():
        key, sep, value = line.partition("=")
        if sep:
            env[key] = value
    return env


def _module_words(entry):
    """Return the module names `. ./mfc.sh load` would pass to `module load`."""
    script = textwrap.dedent(
        """\
        {function_src}
        __module_words {entry}
        """
    ).format(function_src=_helpers("__module_words"), entry=_quote(entry))

    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True)
    return proc.stdout.strip()


def _load_with_stderr(entry):
    """Run one entry through the exporter, returning its environment and stderr."""
    script = textwrap.dedent(
        """\
        log() {{ :; }}          # modules.sh logging stub
        {function_src}
        __export_assignments {entry}
        env
        """
    ).format(function_src=_helpers("__export_assignments"), entry=_quote(entry))

    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True)
    env = {}
    for line in proc.stdout.splitlines():
        key, sep, value = line.partition("=")
        if sep:
            env[key] = value
    return env, proc.stderr.strip()


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


@pytest.mark.parametrize("value", ["-n", "-e", "-E", "-n -e trailing"])
def test_value_that_looks_like_an_echo_option_survives(value):
    """A value starting with -n/-e must be exported, not read as an option.

    Expanding through `echo` made bash's builtin treat these as flags, so
    `MYVAR=-n` exported an empty string. Expansion now runs through the
    positional parameters instead.
    """
    env = _load(f"MYVAR={value}")
    assert env["MYVAR"] == value


def test_module_names_ahead_of_the_first_assignment_are_not_exported():
    """A line may carry module names before its assignments; only the latter are exported.

    The names are not identifiers, so handing them to ``export`` fails. Nothing
    reads that failure -- ``. ./mfc.sh load`` output is routinely redirected --
    which is the same silence the multi-word bug hid behind, so assert on it.
    """
    env, stderr = _load_with_stderr("python cmake CC=nvc CXX=nvc++")

    assert env["CC"] == "nvc"
    assert env["CXX"] == "nvc++"
    assert stderr == ""


def test_return_status_is_zero_when_the_caller_already_set_noglob():
    """`set -f` in the caller must not make the loader report failure."""
    script = textwrap.dedent(
        """\
        log() {{ :; }}
        {function_src}
        set -f
        __export_assignments 'CC=gcc'
        echo "status=$?"
        """
    ).format(function_src=_helpers("__export_assignments"))
    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True)
    assert "status=0" in proc.stdout


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


# The module-name side of the same line. `module load` receives whatever these
# return, and an unknown name aborts `. ./mfc.sh load` outright, so a value word
# leaking into this list is a hard failure rather than a cosmetic one.


@pytest.mark.parametrize(
    "entry,expected",
    [
        ("python cmake/3.22.2", "python cmake/3.22.2"),
        ("nvhpc/22.9 cuda/11.7 openmpi/4.0.5-nvhpc22.9", "nvhpc/22.9 cuda/11.7 openmpi/4.0.5-nvhpc22.9"),
        ("CC=nvc CXX=nvc++ FC=nvfortran", ""),
        ("python cmake CC=nvc CXX=nvc++", "python cmake"),
    ],
)
def test_module_names_are_separated_from_assignments(entry, expected):
    """Module names come first; everything from the first NAME= on is a value."""
    assert _module_words(entry) == expected


@pytest.mark.parametrize(
    "entry",
    [
        # Words of a multi-word value that happen to lack '=' of their own.
        "LDFLAGS=-L$CUDA_HOME/lib64 -lcudart",
        "CFLAGS=-O2 -march=native CC=gcc",
        'CRAY_CCE_LLD_ARGS="-O2 -g"',
        "MYFLAGS=-Wl,* -O2",
    ],
)
def test_value_words_are_never_passed_to_module_load(entry):
    """A value's own words must not be mistaken for module names.

    The classification step used to drop every word containing '=' and keep the
    rest, so `LDFLAGS=-L$CUDA_HOME/lib64 -lcudart` handed `-lcudart` to
    `module load`. That fails and trips `error "Failed to load modules."`,
    aborting the loader before any variable is exported.
    """
    assert _module_words(entry) == ""


def test_module_names_in_toolchain_modules_are_unchanged():
    """Every line shipped today yields exactly the module list it always did."""
    modules_file = MODULES_SH.resolve().parents[1] / "modules"
    checked = 0
    for line in modules_file.read_text().splitlines():
        match = re.match(r"^[a-z][a-z0-9]*-(?:all|cpu|gpu)\s+(.*)$", line)
        if not match:
            continue
        entry = match.group(1)
        # What the pre-existing implementation produced: drop '='-bearing words.
        expected = " ".join(word for word in entry.split() if "=" not in word)
        assert _module_words(entry) == expected, entry
        checked += 1

    assert checked, "expected at least one module line in toolchain/modules"
