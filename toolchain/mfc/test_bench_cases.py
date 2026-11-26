"""bench --cases reruns chosen benchmark cases.

The PR template asks an author to rerun a case whose grind speedup is below the
threshold. Running all of bench.yaml again for one case is far slower than
needed.
"""

from pathlib import Path

import pytest
import yaml

from mfc import bench, state
from mfc.cli.argparse_gen import generate_parser
from mfc.cli.commands import MFC_CLI_SCHEMA
from mfc.common import MFC_BENCH_FILEPATH, MFCException

REPO = Path(__file__).resolve().parents[2]
ALL = [case["slug"] for case in yaml.safe_load(Path(MFC_BENCH_FILEPATH).read_text())]


def _parse(*argv):
    parser, _ = generate_parser(MFC_CLI_SCHEMA, state.MFCConfig())
    return vars(parser.parse_args(["bench", "-o", "out.yaml", *argv]))


def _bench(tmp_path, monkeypatch, *argv):
    """Run bench() as the command line would, with the case runs faked.

    Returns the slugs that were run and the slugs in the results file.
    """
    ran = []

    def fake_run(command, **_):
        ran.append(Path(command[2]).parent.name)
        summary = Path(command[command.index("--output-summary") + 1])
        summary.write_text(yaml.safe_dump({"simulation": {"exec": 1.0, "grind": 1.0}}))
        return 0

    monkeypatch.chdir(REPO)
    monkeypatch.setattr(bench, "system", fake_run)
    monkeypatch.setattr(bench, "MFC_BUILD_DIR", str(tmp_path))
    output = tmp_path / "out.yaml"
    args = {**_parse(*argv), "targets": ["simulation"], "output": str(output), "--": []}
    for key, value in args.items():
        monkeypatch.setitem(state.gARG, key, value)

    bench.bench()

    return ran, list(yaml.safe_load(output.read_text())["cases"])


def test_runs_every_case_by_default(tmp_path, monkeypatch):
    assert _bench(tmp_path, monkeypatch) == (ALL, sorted(ALL))


def test_runs_only_the_named_cases(tmp_path, monkeypatch):
    ran, results = _bench(tmp_path, monkeypatch, "--cases", "igr", "ibm")

    # bench.yaml order, whatever order they were named in.
    assert ran == ["ibm", "igr"]
    assert results == ["ibm", "igr"]


def test_an_unknown_case_fails_before_anything_runs(tmp_path, monkeypatch):
    ran = []
    monkeypatch.setattr(bench, "system", lambda *_, **__: ran.append(1))

    with pytest.raises(MFCException) as error:
        monkeypatch.chdir(REPO)
        monkeypatch.setattr(bench, "MFC_BUILD_DIR", str(tmp_path))
        for key, value in {**_parse("--cases", "ibm", "hllc"), "targets": ["simulation"], "--": []}.items():
            monkeypatch.setitem(state.gARG, key, value)
        bench.bench()

    assert "hllc" in str(error.value)
    assert all(slug in str(error.value) for slug in ALL)
    assert ran == []
