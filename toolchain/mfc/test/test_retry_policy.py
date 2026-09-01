"""Tests for the test-suite retry decision (issue #1798).

Retries in this harness cost 3x wall clock and, as far as anyone can measure,
rescue almost nothing: bench.py's equivalent rescued 0 of 235 retried cases, and
all 2,795 recorded failed-test entries show `Attempts: 3`.

That second figure is one-sided by construction -- only failures record an
attempt count, so a case that passes on attempt 2 is indistinguishable from one
that passed first try. Measuring rescues is therefore the prerequisite for any
policy change, and is fixed here.

Also fixed: the retry loop never consulted the suite-wide abort flag. On a bad
node, where the 30% failure-rate abort is exactly what should stop the run,
every in-flight case still burned its remaining attempts first.
"""

from mfc.test.test import should_retry


def test_a_failure_is_retried_until_the_attempt_budget_is_spent():
    assert should_retry(attempt=1, max_attempts=3, aborting=False)
    assert should_retry(attempt=2, max_attempts=3, aborting=False)


def test_the_last_attempt_is_not_retried():
    assert not should_retry(attempt=3, max_attempts=3, aborting=False)


def test_a_single_attempt_budget_never_retries():
    assert not should_retry(attempt=1, max_attempts=1, aborting=False)


def test_no_retry_once_the_suite_is_aborting():
    # The abort fires when the failure rate says the environment is broken --
    # a bad node, say. Spending two more attempts per in-flight case is the
    # opposite of the fail-fast the abort exists to provide.
    assert not should_retry(attempt=1, max_attempts=3, aborting=True)
    assert not should_retry(attempt=2, max_attempts=3, aborting=True)


def test_the_summary_reports_how_often_a_retry_actually_helped(capsys):
    from mfc.test.test import _print_test_summary

    _print_test_summary(10, 2, 0, 0, 1.0, [], [], rescued=3)
    assert "recovered by a retry" in capsys.readouterr().out


def test_the_summary_stays_quiet_when_no_retry_helped(capsys):
    from mfc.test.test import _print_test_summary

    _print_test_summary(10, 2, 0, 0, 1.0, [], [], rescued=0)
    assert "recovered by a retry" not in capsys.readouterr().out
