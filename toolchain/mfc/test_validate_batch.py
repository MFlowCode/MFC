"""Batch validation (precheck's example check) uses the default build configuration, whatever the local build was."""

from mfc import state, validate


def test_batch_validation_ignores_the_local_build(monkeypatch):
    seen = []
    monkeypatch.setattr(validate, "_validate_batch_default_build", lambda inputs: seen.append(state.CFG()))
    local = state.MFCConfig(mpi=False, single=True, mixed=True)
    monkeypatch.setattr(state, "gCFG", local)
    validate._validate_batch(["case.py"])
    assert seen == [state.MFCConfig()]
    assert state.CFG() is local
