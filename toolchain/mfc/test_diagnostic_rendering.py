"""Captured diagnostics must survive being printed.

The console prints through Rich with markup enabled (printer.py), so raw log
text is interpreted as markup. Compiler and MPI output is full of square
brackets -- absolute paths in diagnostics, `[host:pid]` rank prefixes -- and
Rich either raises MarkupError on them or silently eats them as tags.

That turns diagnostic capture into the opposite of its purpose: a crash, or a
message with the identifying parts removed. It is worse inside an MFCException,
because main.py renders that message with markup from inside the handler, so the
failure surfaces as an unhandled traceback instead of the intended report.
"""

from rich.console import Console

from mfc.common import console_safe, get_program_output

HOSTILE = [
    "ftn: error in [/lustre/orion/cfd154/scratch/x.f90] line 3",
    "[node1:12345] MPI abort",
    "nvlink error: undefined reference in [/gpfs/alpine/proj/m_riemann.o]",
]


def render(text):
    """Render the way the CLI does: a Rich console with markup enabled."""
    console = Console(file=open("/dev/null", "w"), record=True, width=200)
    console.print(text, soft_wrap=True)
    return console.export_text()


def test_rich_really_does_mangle_raw_log_text():
    # Guards the premise: if Rich ever stops doing this, console_safe can go.
    raised = False
    try:
        render(HOSTILE[0])
    except Exception:
        raised = True
    assert raised or "[/lustre/orion/cfd154/scratch/x.f90]" not in render(HOSTILE[0])


def test_console_safe_text_renders_without_raising():
    for line in HOSTILE:
        render(console_safe(line))


def test_console_safe_preserves_bracketed_paths():
    out = render(console_safe(HOSTILE[0]))
    assert "[/lustre/orion/cfd154/scratch/x.f90]" in out


def test_console_safe_preserves_mpi_rank_prefixes():
    # "[node1:12345]" is how you tell which rank died; Rich eats it as a tag.
    assert "[node1:12345]" in render(console_safe(HOSTILE[1]))


def test_get_program_output_can_capture_stderr():
    # h5dump writes "unable to open file" to stderr, so without this the
    # diagnostic added to the h5dump failure path is always "(no output)".
    out, code = get_program_output(["bash", "-c", "echo to-stderr >&2; exit 3"], merge_stderr=True)
    assert code == 3
    assert "to-stderr" in out


def test_get_program_output_still_ignores_stderr_by_default():
    # Other callers parse stdout and must not start seeing stderr mixed in.
    out, _ = get_program_output(["bash", "-c", "echo out; echo err >&2"])
    assert "err" not in out
    assert "out" in out
