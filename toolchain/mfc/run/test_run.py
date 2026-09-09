"""Tests for run.run() option validation ordering (issue #1511).

An invalid invocation -- ``--no-mpi`` with more than one rank, ``nodes <= 0``,
a malformed ``--email`` -- must be rejected before run() touches the case
directory. Previously the job script was rendered and written to disk first,
so a rejected run left a stale ``<name>.sh`` behind (and, with ``--clean``, had
already wiped the previous run's outputs).
"""

import types
import unittest
from contextlib import ExitStack
from unittest.mock import Mock, patch

from .. import state
from ..common import MFCException
from . import run as run_mod


def _fake_target(name):
    return types.SimpleNamespace(name=name)


def _fake_case(clean):
    return types.SimpleNamespace(params={}, clean=clean)


class _StateSandbox(unittest.TestCase):
    def setUp(self):
        self._saved_gARG = dict(state.gARG)
        state.gARG.clear()
        state.gARG.update(
            {
                "name": "MFC",
                "input": "case.py",
                "targets": ["simulation"],
                "targets_explicit": True,
                "engine": "interactive",
                "mpi": True,
                "nodes": 1,
                "tasks_per_node": 1,
                "email": "",
                "verbose": 0,
                "clean": True,
                "archive": None,
                "dry_run": True,
                "output_summary": None,
            }
        )

    def tearDown(self):
        state.gARG.clear()
        state.gARG.update(self._saved_gARG)

    def _patched_run(self):
        """Patch every side-effecting collaborator of run() with a Mock.

        Module-level dunder names are not mangled, but attribute access from
        inside a class body would be, so the generators are swapped through
        ``__dict__`` (same trick as test_archive.py).
        """
        stack = ExitStack()
        mocks = {
            "build": stack.enter_context(patch.object(run_mod, "build")),
            "generate_job_script": Mock(),
            "generate_input_files": Mock(),
            "clean": Mock(),
        }
        stack.enter_context(patch.object(run_mod, "get_targets", side_effect=lambda names: [_fake_target(n) for n in names]))
        stack.enter_context(
            patch.dict(
                run_mod.__dict__,
                {
                    "__generate_job_script": mocks["generate_job_script"],
                    "__generate_input_files": mocks["generate_input_files"],
                },
            )
        )
        return stack, mocks


class TestInvalidOptionsAreRejectedBeforeSideEffects(_StateSandbox):
    def _assert_rejected_cleanly(self):
        stack, mocks = self._patched_run()
        with stack:
            with self.assertRaises(MFCException):
                run_mod.run(targets=["simulation"], case=_fake_case(mocks["clean"]))

            self.assertFalse(mocks["generate_job_script"].called, "job script was written before option validation")
            self.assertFalse(mocks["generate_input_files"].called, "input files were written before option validation")
            self.assertFalse(mocks["clean"].called, "case was cleaned before option validation")
            self.assertFalse(mocks["build"].called, "build ran before option validation")

    def test_no_mpi_with_multiple_ranks(self):
        state.gARG.update({"mpi": False, "tasks_per_node": 4})
        self._assert_rejected_cleanly()

    def test_non_positive_nodes(self):
        state.gARG.update({"nodes": 0})
        self._assert_rejected_cleanly()

    def test_non_positive_tasks_per_node(self):
        state.gARG.update({"tasks_per_node": 0})
        self._assert_rejected_cleanly()

    def test_malformed_email(self):
        state.gARG.update({"email": "not-an-address"})
        self._assert_rejected_cleanly()


class TestValidOptionsStillRun(_StateSandbox):
    def test_valid_options_reach_generation(self):
        stack, mocks = self._patched_run()
        with stack:
            run_mod.run(targets=["simulation"], case=_fake_case(mocks["clean"]))

            self.assertTrue(mocks["build"].called)
            self.assertTrue(mocks["clean"].called)
            self.assertTrue(mocks["generate_job_script"].called)
            self.assertTrue(mocks["generate_input_files"].called)
