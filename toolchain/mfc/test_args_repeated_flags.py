"""A repeated list-valued flag must not silently discard the earlier values.

`-t` is declared with nargs="+", so it takes a space-separated list. Repeating the
flag does not append: argparse's default "store" action keeps only the last
occurrence, so `-t pre_process -t simulation` builds a script containing only
`simulation`. The queued job then runs simulation against an empty restart_data/
and nothing in the output says a target was dropped.

These tests pin the warning, and pin that the parse result itself is untouched.
"""

import sys
from unittest.mock import patch

from mfc.args import _multi_value_flags, _warn_on_repeated_multi_value_flags, parse
from mfc.cli.argparse_gen import generate_parser
from mfc.cli.commands import MFC_CLI_SCHEMA
from mfc.state import MFCConfig


def _warnings_for(command, argv):
    """Return the warning strings emitted for `argv` under `command`."""
    emitted = []
    with patch("mfc.args.cons") as printer:
        printer.print.side_effect = lambda *a, **k: emitted.append(" ".join(str(x) for x in a))
        _warn_on_repeated_multi_value_flags(command, argv)
    return emitted


def test_repeated_short_targets_flag_warns():
    """The exact invocation from the report: -t pre_process -t simulation."""
    warnings = _warnings_for("run", ["run", "case.py", "-t", "pre_process", "-t", "simulation"])
    assert len(warnings) == 1
    assert "-t" in warnings[0]
    assert "--targets" in warnings[0]
    assert "2 times" in warnings[0]


def test_repeated_long_targets_flag_warns():
    """--targets is the same option and must be counted alongside -t."""
    warnings = _warnings_for("run", ["run", "case.py", "--targets", "pre_process", "--targets", "simulation"])
    assert len(warnings) == 1
    assert "2 times" in warnings[0]


def test_mixed_short_and_long_forms_warn():
    """Mixing the two spellings of one option still drops the first list."""
    warnings = _warnings_for("run", ["run", "case.py", "-t", "pre_process", "--targets", "simulation"])
    assert len(warnings) == 1
    assert "2 times" in warnings[0]


def test_equals_form_is_counted():
    """--targets=simulation is the same occurrence, written differently."""
    warnings = _warnings_for("run", ["run", "case.py", "-t", "pre_process", "--targets=simulation"])
    assert len(warnings) == 1
    assert "2 times" in warnings[0]


def test_three_repeats_reports_the_real_count():
    warnings = _warnings_for("run", ["run", "case.py", "-t", "pre_process", "-t", "simulation", "-t", "post_process"])
    assert len(warnings) == 1
    assert "3 times" in warnings[0]


def test_intended_single_flag_form_is_silent():
    """The documented form -- one flag, several values -- must not warn."""
    assert _warnings_for("run", ["run", "case.py", "-t", "pre_process", "simulation"]) == []


def test_one_target_is_silent():
    """`-t simulation` alone is a legitimate restart invocation."""
    assert _warnings_for("run", ["run", "case.py", "-t", "simulation"]) == []


def test_no_flag_is_silent():
    assert _warnings_for("run", ["run", "case.py"]) == []


def test_repeated_flag_on_build_warns():
    """`targets` reaches build through include_common, not its own argument list."""
    warnings = _warnings_for("build", ["build", "-t", "pre_process", "-t", "simulation"])
    assert len(warnings) == 1


def test_other_list_valued_flags_are_covered():
    """The check is schema-driven, so every nargs="+" flag is covered, not just -t."""
    warnings = _warnings_for("run", ["run", "case.py", "-g", "0", "-g", "1"])
    assert len(warnings) == 1
    assert "--gpus" in warnings[0]


def test_a_value_that_looks_like_the_flag_is_not_miscounted():
    """A test UUID or filename equal to a flag name must not inflate the count."""
    assert _warnings_for("run", ["run", "case.py", "-t", "simulation"]) == []


def test_unknown_command_is_silent():
    """A command not in the schema must not raise on the way to argparse's error."""
    assert _warnings_for("not-a-command", ["not-a-command", "-t", "a", "-t", "b"]) == []


def test_multi_value_flags_finds_targets_for_run():
    """The schema walk must see both the command's own args and its common sets."""
    dests = {dest for _, dest in _multi_value_flags("run")}
    assert "targets" in dests
    assert "gpus" in dests


def test_multi_value_flags_skips_scalar_options():
    """Options without nargs cannot be silently overwritten and must not be listed."""
    dests = {dest for _, dest in _multi_value_flags("test")}
    # `--from` / `--to` are plain scalars; repeating them loses nothing meaningful.
    assert "from" not in dests
    assert "to" not in dests


def test_parse_still_returns_argparses_result_for_a_repeated_flag():
    """The warning is advisory: parsing behaviour is deliberately left unchanged."""
    argv = ["./mfc.sh", "run", "case.py", "-t", "pre_process", "-t", "simulation"]
    with patch.object(sys, "argv", argv), patch("mfc.args.cons"):
        args = parse(MFCConfig())
    assert args["targets"] == ["simulation"]


def test_warning_is_consistent_with_what_argparse_actually_did():
    """Whenever we warn, argparse must genuinely have kept only the last list."""
    parser, _ = generate_parser(MFC_CLI_SCHEMA, MFCConfig())
    argv = ["run", "case.py", "-t", "pre_process", "-t", "simulation"]
    parsed = parser.parse_args(argv)
    assert parsed.targets == ["simulation"]
    assert _warnings_for("run", argv) != []
