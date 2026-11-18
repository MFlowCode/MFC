"""Tests for the Lagrangian output files in toolchain/mfc/packer/pack.py.

Each rank writes its own lag_bubble_evol_<rank>.dat and stats_lag_bubbles_<rank>.dat,
starting with a header line. A rank that owns no bubble writes the header and
nothing else, which a one-rank run never produces: with one bubble on two ranks,
the packer failed on that file and the test could not be run.
"""

from mfc.packer import pack


def _case(tmp_path, files):
    for name, text in files.items():
        (tmp_path / "D").mkdir(exist_ok=True)
        (tmp_path / "D" / name).write_text(text)
    return pack.compile(str(tmp_path))


HEADER = " currentTime particleID x y z\n"


def test_a_rank_with_no_bubble_packs_as_an_empty_entry(tmp_path):
    packed, err = _case(tmp_path, {"lag_bubble_evol_0.dat": HEADER, "lag_bubble_evol_1.dat": HEADER + " 0.1 1 0.5 0.5 0.0\n", "stats_lag_bubbles_0.dat": " proc_rank, particleID, x\n"})
    assert err is None, err
    assert packed.find("D/lag_bubble_evol_0.dat").doubles == []
    assert packed.find("D/stats_lag_bubbles_0.dat").doubles == []
    assert packed.find("D/lag_bubble_evol_1.dat").doubles == [0.1, 1.0, 0.5, 0.5, 0.0]


def test_a_malformed_lagrangian_file_is_still_an_error(tmp_path):
    _, err = _case(tmp_path, {"lag_bubble_evol_0.dat": HEADER + " 0.1 nonsense\n"})
    assert err is not None and "lag_bubble_evol_0.dat" in err


def test_an_empty_entry_survives_a_golden_round_trip(tmp_path):
    """The golden is written and read back: an empty entry is written as the path and a trailing space."""
    packed, err = _case(tmp_path, {"lag_bubble_evol_0.dat": HEADER, "lag_bubble_evol_1.dat": HEADER + " 0.1 1 0.5 0.5 0.0\n"})
    assert err is None, err
    golden = tmp_path / "golden.txt"
    golden.write_text("\n".join(str(e) for e in packed.entries.values()))  # the golden lines Pack.save writes
    loaded = pack.load(str(golden))
    assert loaded.find("D/lag_bubble_evol_0.dat").doubles == []
    assert loaded.find("D/lag_bubble_evol_1.dat").doubles == [0.1, 1.0, 0.5, 0.5, 0.0]
