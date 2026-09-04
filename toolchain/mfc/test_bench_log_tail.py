"""A failing benchmark case must say why in the CI log.

bench.py writes each case's output to build/benchmarks/<hash>/<case>.out and, on
failure, prints only that path -- a path on the cluster that no artifact upload
ever collects. The result is a CI failure whose entire diagnosis is
"failed with exit code 143". Thirteen jobs failed that way over 2026-08-18..31,
eleven of them on Phoenix, and none can be debugged after the fact.
"""

from mfc.common import log_tail


def test_shows_the_end_of_the_log_where_the_error_is(tmp_path):
    log = tmp_path / "case.out"
    log.write_text("start\nmiddle\nNaN(s) in timestep output\n")
    assert "NaN(s) in timestep output" in log_tail(str(log))


def test_keeps_only_the_last_lines_so_a_huge_log_cannot_flood_ci(tmp_path):
    log = tmp_path / "case.out"
    log.write_text("\n".join(f"line {i}" for i in range(5000)))
    body = log_tail(str(log), max_lines=50)
    assert "line 4999" in body
    assert "line 4000" not in body


def test_says_so_when_the_log_was_never_written(tmp_path):
    # A case killed before it opened its log is itself a useful signal, and must
    # not turn into a traceback inside the bench harness.
    body = log_tail(str(tmp_path / "absent.out"))
    assert "absent.out" in body
    assert body.strip() != ""


def test_says_so_when_the_log_is_empty(tmp_path):
    log = tmp_path / "case.out"
    log.write_text("")
    assert log_tail(str(log)).strip() != ""


def test_names_the_log_so_the_full_file_can_still_be_found(tmp_path):
    log = tmp_path / "case.out"
    log.write_text("boom\n")
    assert "case.out" in log_tail(str(log))


def test_reads_only_the_tail_of_a_large_log(tmp_path):
    # Bounded memory matters: this runs on an already-failing path, and solver
    # and benchmark logs reach tens of MB.
    import tracemalloc

    log = tmp_path / "big.out"
    log.write_text("".join(f"line {i} {'x' * 200}\n" for i in range(200_000)))
    assert log.stat().st_size > 40_000_000

    tracemalloc.start()
    body = log_tail(str(log), max_lines=20)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    assert "line 199999" in body
    assert "line 100000" not in body
    assert peak < 5_000_000, f"peak {peak} suggests the whole file was read"
