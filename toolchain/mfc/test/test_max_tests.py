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


def test_cases_marked_keep_survive_and_count_toward_the_cap():
    cases = list(range(100))
    kept, dropped = cap_evenly(cases, 10, keep=lambda c: c in (3, 97))
    assert 3 in kept and 97 in kept
    assert len(kept) == 10 and sorted(kept + dropped) == cases


def test_more_kept_cases_than_the_cap_all_run():
    kept, _ = cap_evenly(list(range(50)), 5, keep=lambda c: c < 8)
    assert kept == list(range(8))
