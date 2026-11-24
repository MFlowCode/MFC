"""Tests for --max-tests, which a CI smoke run uses to sample the selected tests."""

from mfc.test.test import cap_evenly


def test_no_cap_and_cap_not_reached_keep_everything():
    cases = list(range(10))
    assert cap_evenly(cases, 0) == (cases, [])
    assert cap_evenly(cases, 10) == (cases, [])
    assert cap_evenly(cases, 99) == (cases, [])


def test_cap_keeps_exactly_cap_cases_spread_over_the_list():
    kept, dropped = cap_evenly(list(range(100)), 10)
    assert kept == [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
    assert sorted(kept + dropped) == list(range(100))


def test_cap_counts_only_cases_that_will_run():
    cases = [("chem", i) if i % 4 == 0 else ("plain", i) for i in range(40)]
    kept, dropped = cap_evenly(cases, 6, lambda c: c[0] == "plain")
    assert sum(c[0] == "plain" for c in kept) == 6
    # The others stay for the later filter to skip, and are not reported as capped.
    assert sum(c[0] == "chem" for c in kept) == 10
    assert all(c[0] == "plain" for c in dropped)


def test_order_is_preserved():
    kept, _ = cap_evenly(list("abcdefghij"), 3)
    assert kept == sorted(kept)
