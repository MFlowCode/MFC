"""Tests for the SCC fingerprint core (toolchain/mfc/scc/core.py)."""

import os

import pytest

from mfc.scc import core

x25519 = pytest.importorskip("cryptography.hazmat.primitives.asymmetric.x25519")
serialization = pytest.importorskip("cryptography.hazmat.primitives.serialization")
InvalidTag = pytest.importorskip("cryptography.exceptions").InvalidTag


def _keys():
    priv = x25519.X25519PrivateKey.generate()
    raw = serialization.Encoding.Raw
    return (
        priv.private_bytes(raw, serialization.PrivateFormat.Raw, serialization.NoEncryption()),
        priv.public_key().public_bytes(raw, serialization.PublicFormat.Raw),
    )


def _write(root, rel, text):
    path = os.path.join(root, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def test_line_roundtrip_and_header_binding():
    priv, pub = _keys()
    line = core.encode_line(pub, "build", "simulation", "sha1", {"x": 1})
    assert core.decode_line(priv, line) == ([core.key_id(pub), "build", "simulation", "sha1", "-"], {"x": 1})
    with pytest.raises(InvalidTag):
        core.decode_line(priv, line.replace(" sha1 ", " sha2 ", 1))
    with pytest.raises(InvalidTag):
        core.decode_line(_keys()[0], line)


def test_tree_state_detects_changes(tmp_path):
    root = str(tmp_path)
    _write(root, "src/simulation/m_rhs.fpp", "a")
    _write(root, "toolchain/templates/x.mako", "b")
    os.makedirs(os.path.join(root, os.path.dirname(core.BASELINE_PATH)))
    baseline = core.write_baseline(root, os.path.join(root, core.BASELINE_PATH))

    _write(root, "src/simulation/m_rhs.fpp", "hacked")
    _write(root, "toolchain/templates/x.mako", "tuned")
    _write(root, "src/simulation/m_new.fpp", "new")
    _write(root, "src/simulation/include/case.fpp", "generated")

    tree = core.tree_state(root, ["toolchain/templates/*"])
    status = {c["path"]: c["status"] for c in tree["changes"]}
    assert status["src/simulation/m_rhs.fpp"] == "M" and status["src/simulation/m_new.fpp"] == "A"
    assert sorted(tree["violations"]) == ["src/simulation/m_new.fpp", "src/simulation/m_rhs.fpp"]
    assert core.digest(core.apply_changes(baseline, tree["changes"])) == tree["digest"]


def test_diff_stat():
    assert core.diff_stat(b"a\nb\nc\n", b"a\nB\nc\nd\n") == (2, 1)
    assert core.diff_stat(b"", b"x\n") == (1, 0)
    assert core.diff_stat(b"\xff", b"") is None
