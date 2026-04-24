"""
Tests for the archive module.

Covers source collection, archive format round-trips (dir, tar, tar.zst),
the no-op path when --archive is unset, destination-collision tally
fallback, and plan_archive() error cases.
"""

import datetime
import os
import shutil
import subprocess
import tarfile
import tempfile
import types
import unittest
from contextlib import contextmanager
from unittest.mock import patch

from .. import state


def _make_fake_case(dirpath: str):
    from .input import MFCInputFile

    case_py = os.path.join(dirpath, "case.py")
    with open(case_py, "w") as f:
        f.write("# fake case\n")
    with open(os.path.join(dirpath, "simulation.inp"), "w") as f:
        f.write("&user_inputs /\n")
    with open(os.path.join(dirpath, "equations.dat"), "w") as f:
        f.write("eq\n")
    with open(os.path.join(dirpath, "MFC.out"), "w") as f:
        f.write("log\n")
    os.makedirs(os.path.join(dirpath, "D"))
    with open(os.path.join(dirpath, "D", "output.dat"), "w") as f:
        f.write("data\n")

    return MFCInputFile(case_py, dirpath, {})


def _fake_targets():
    return [types.SimpleNamespace(name="simulation")]


def _collect_sources_fn():
    # Module-level dunder names aren't mangled; attribute access inside a class body
    # would be, so we reach through __dict__.
    from . import archive

    return archive.__dict__["__collect_sources"]


class _StateSandbox(unittest.TestCase):
    def setUp(self):
        self._saved_gARG = dict(state.gARG)
        state.gARG.update({"name": "MFC", "output_summary": None})

    def tearDown(self):
        state.gARG.clear()
        state.gARG.update(self._saved_gARG)


class TestCollectSources(_StateSandbox):
    def test_finds_case_namelist_and_artifacts(self):
        collect = _collect_sources_fn()

        with tempfile.TemporaryDirectory() as tmp:
            case = _make_fake_case(tmp)
            sources = collect(case, _fake_targets())

        names = {os.path.basename(s) for s in sources}
        self.assertIn("case.py", names)
        self.assertIn("simulation.inp", names)
        self.assertIn("equations.dat", names)
        self.assertIn("MFC.out", names)
        self.assertIn("D", names)

    def test_skips_missing_artifacts(self):
        collect = _collect_sources_fn()

        with tempfile.TemporaryDirectory() as tmp:
            case = _make_fake_case(tmp)
            sources = collect(case, _fake_targets())

        names = {os.path.basename(s) for s in sources}
        self.assertNotIn("time_data.dat", names)
        self.assertNotIn("restart_data", names)
        self.assertNotIn("p_all", names)


@contextmanager
def _run_archive(fmt: str):
    from . import archive as archive_mod

    src = tempfile.mkdtemp()
    dest_root = tempfile.mkdtemp()
    try:
        state.gARG["archive"] = dest_root
        state.gARG["archive_format"] = fmt
        case = _make_fake_case(src)
        plan = archive_mod.plan_archive(case)
        archive_mod.archive(plan, case, _fake_targets())
        entries = sorted(os.listdir(dest_root))
        assert len(entries) == 1, f"expected one archive entry, got {entries}"
        yield os.path.join(dest_root, entries[0])
    finally:
        shutil.rmtree(src, ignore_errors=True)
        shutil.rmtree(dest_root, ignore_errors=True)


class TestArchiveFormats(_StateSandbox):
    def test_format_dir(self):
        with _run_archive("dir") as path:
            self.assertTrue(os.path.isdir(path))
            self.assertTrue(os.path.isfile(os.path.join(path, "manifest.yaml")))
            self.assertTrue(os.path.isfile(os.path.join(path, "case.py")))
            self.assertTrue(os.path.isfile(os.path.join(path, "simulation.inp")))
            self.assertTrue(os.path.isdir(os.path.join(path, "D")))
            self.assertTrue(os.path.isfile(os.path.join(path, "D", "output.dat")))

    def test_format_tar(self):
        with _run_archive("tar") as path:
            self.assertTrue(path.endswith(".tar"))
            self.assertTrue(tarfile.is_tarfile(path))
            with tarfile.open(path) as tf:
                names = tf.getnames()
            base = os.path.basename(path)[: -len(".tar")]
            self.assertIn(f"{base}/manifest.yaml", names)
            self.assertIn(f"{base}/case.py", names)
            self.assertIn(f"{base}/D/output.dat", names)

    def test_format_tar_zst(self):
        try:
            r = subprocess.run(["tar", "--zstd", "--version"], capture_output=True, check=False, timeout=5)
            if r.returncode != 0:
                self.skipTest("tar --zstd not available")
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.skipTest("tar --zstd not available")

        with _run_archive("tar.zst") as path:
            self.assertTrue(path.endswith(".tar.zst"))
            self.assertGreater(os.path.getsize(path), 0)

            listing = subprocess.run(["tar", "--zstd", "-tf", path], capture_output=True, text=True, check=True)
            base = os.path.basename(path)[: -len(".tar.zst")]
            self.assertIn(f"{base}/manifest.yaml", listing.stdout)
            self.assertIn(f"{base}/case.py", listing.stdout)


class TestArchiveBehavior(_StateSandbox):
    def test_plan_returns_none_when_archive_unset(self):
        from . import archive as archive_mod

        state.gARG["archive"] = None
        state.gARG["archive_format"] = "dir"
        # No case needed since plan_archive returns before touching it.
        self.assertIsNone(archive_mod.plan_archive(case=None))

    def test_plan_bad_format_raises(self):
        from ..common import MFCException
        from . import archive as archive_mod

        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as dest_root:
            state.gARG["archive"] = dest_root
            state.gARG["archive_format"] = "bogus"
            case = _make_fake_case(src)
            with self.assertRaises(MFCException):
                archive_mod.plan_archive(case)

    def test_plan_uses_case_dir_name_as_stem(self):
        from . import archive as archive_mod

        fixed = datetime.datetime(2026, 1, 1, 12, 0, 0)

        with tempfile.TemporaryDirectory() as parent, tempfile.TemporaryDirectory() as dest_root, patch("mfc.run.archive.datetime.datetime") as MockDT:
            MockDT.now.return_value = fixed
            case_dir = os.path.join(parent, "my_cool_case")
            os.makedirs(case_dir)
            state.gARG["archive"] = dest_root
            state.gARG["archive_format"] = "dir"
            case = _make_fake_case(case_dir)

            plan = archive_mod.plan_archive(case)

        self.assertEqual(plan.stem, "my_cool_case-20260101-120000")
        self.assertTrue(plan.dest.endswith("my_cool_case-20260101-120000"))

    def test_plan_collision_gets_tally_suffix(self):
        from . import archive as archive_mod

        fixed = datetime.datetime(2026, 1, 1, 12, 0, 0)

        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as dest_root, patch("mfc.run.archive.datetime.datetime") as MockDT:
            MockDT.now.return_value = fixed
            state.gARG["archive"] = dest_root
            state.gARG["archive_format"] = "dir"
            case = _make_fake_case(src)

            plan1 = archive_mod.plan_archive(case)
            os.makedirs(plan1.dest)  # simulate an existing archive at that path
            plan2 = archive_mod.plan_archive(case)
            os.makedirs(plan2.dest)
            plan3 = archive_mod.plan_archive(case)

        self.assertTrue(plan2.dest.endswith("-2"))
        self.assertTrue(plan3.dest.endswith("-3"))
        self.assertNotEqual(plan1.dest, plan2.dest)
        self.assertNotEqual(plan2.dest, plan3.dest)

    def test_output_summary_outside_case_dir_is_skipped(self):
        collect = _collect_sources_fn()

        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as elsewhere:
            outside = os.path.join(elsewhere, "summary.yaml")
            with open(outside, "w") as f:
                f.write("k: v\n")

            state.gARG["output_summary"] = outside
            case = _make_fake_case(src)
            sources = collect(case, _fake_targets())

        self.assertNotIn(outside, sources)


if __name__ == "__main__":
    unittest.main()
