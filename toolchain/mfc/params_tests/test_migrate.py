"""Tests for case-file migration to named parameter values."""

SAMPLE = """print(json.dumps({
    "riemann_solver": 2,
    "model_eqns": 2,
    "weno_order": 5,
    "num_fluids": 1,
    'time_stepper': 3,
}))
"""


def test_migrate_text_rewrites_named_params():
    from mfc.params.migrate import migrate_text

    out, n = migrate_text(SAMPLE)
    assert '"riemann_solver": "hllc",' in out
    assert '"model_eqns": "5eq",' in out
    assert "'time_stepper': \"rk3\"," in out
    assert '"weno_order": 5,' in out  # not an enumerated vocabulary - untouched
    assert '"num_fluids": 1,' in out  # has min/max but no names - untouched
    assert n == 3


def test_migrate_text_is_idempotent():
    from mfc.params.migrate import migrate_text

    once, _ = migrate_text(SAMPLE)
    twice, n = migrate_text(once)
    assert twice == once
    assert n == 0


def test_migrate_text_leaves_unknown_codes():
    from mfc.params.migrate import migrate_text

    out, n = migrate_text('"riemann_solver": 9,')
    assert out == '"riemann_solver": 9,'
    assert n == 0


def test_migrate_text_rewrites_comments_too():
    # Known limitation of text-based migration: commented-out parameter
    # lines are rewritten as well. Semantics are unaffected (comments do
    # not execute); this test documents the behavior.
    from mfc.params.migrate import migrate_text

    out, n = migrate_text('# "riemann_solver": 2,\n"riemann_solver": 2,\n')
    assert out == '# "riemann_solver": "hllc",\n"riemann_solver": "hllc",\n'
    assert n == 2


def test_migrate_text_handles_inline_comment_without_comma():
    from mfc.params.migrate import migrate_text

    out, n = migrate_text('"time_stepper": 3  # final entry, no trailing comma\n')
    assert out == '"time_stepper": "rk3"  # final entry, no trailing comma\n'
    assert n == 1


def test_migrate_text_handles_end_of_string():
    from mfc.params.migrate import migrate_text

    out, n = migrate_text('"format": 1')
    assert out == '"format": "silo"'
    assert n == 1
