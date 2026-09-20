"""Check that no test-suite case needs a build of its own.

A case gets its own hashed install directory whenever anything feeding
`MFCTarget.get_slug` differs from the default build: an analytic initial condition (a
`patch_icpp` expression, codegen'd into a per-case `case.fpp`), a chemistry mechanism,
or any future `case.fpp` input. The test jobs run `./mfc.sh test --no-build`, so a case
whose slug was never pre-built is handed a binary that does not exist:

    error: execve(): .../build/install/gpu-mp-eca8fb2b47/bin/pre_process: No such file
    or directory

What gets pre-built is lane-dependent. Most lanes run `./mfc.sh test --dry-run -a`,
which builds every case-specific variant. The Frontier AMD gpu-omp lane cannot: each
amdflang device link is ~1 h, so `.github/workflows/test.yml` splits its build into two
concurrent jobs, and only these two slugs exist afterwards:

    base -> ./mfc.sh build                                  (the default build)
    chem -> ./mfc.sh test --dry-run -a -o Chemistry         (cases whose trace has a
                                                             "Chemistry" segment)
    eos  -> ./mfc.sh test --dry-run -a -o eos=mie_gruneisen (cases whose trace has an
                                                             "eos=mie_gruneisen" segment)

So the invariant checked here is that every golden test case's slug is one of those: the
default build's, or one shared with a case `--only Chemistry` or `--only eos=mie_gruneisen`
selects. The eos variant exists because eos_state_dependent is a compile-time constant in
every build (Case.get_fpp's _prepend, like chemistry), so a Mie-Gruneisen, JWL or Vinet
fluid moves the slug. All such cases share one slug per target, so selecting on the
Mie-Gruneisen label alone covers the JWL and Vinet ones too. A case that
fails this check is fine on every other lane and red on that one, roughly 1.5 h into
the run, which is why it is worth catching in the lint gate instead.

Fixes, in order of preference:
  * an example: add it to `casesToSkip` in `toolchain/mfc/test/cases.py`, as every
    `*_convergence` and `*_analytical` example already is;
  * a suite case that needs spatial variation: get it geometrically rather than
    analytically (e.g. a cuboid patch, whose bounds are plain numbers);
  * a genuinely new chemistry mechanism: give the case a "Chemistry" trace segment so
    the `chem` pre-build covers it.

If that lane's build commands change, update `_allowed_slugs` to match them.
"""

from __future__ import annotations

import contextlib
import io
import os
import sys
from pathlib import Path

# The targets a golden test can launch. The slug is per target, so all three are checked.
TARGET_NAMES = ("pre_process", "simulation", "post_process")

CHEMISTRY_LABEL = "Chemistry"

# The label the eos pre-build selects on; see the module docstring.
EOS_LABEL = "eos=mie_gruneisen"


def _import_toolchain(repo_root: Path):
    toolchain_dir = str(repo_root / "toolchain")
    if toolchain_dir not in sys.path:
        sys.path.insert(0, toolchain_dir)

    from mfc import state

    # get_fpp reads this argument, which only a real ./mfc.sh invocation sets. The test
    # suite never builds with case optimization, so the default is the honest value.
    state.gARG.setdefault("case_optimization", False)

    from mfc.build import get_target
    from mfc.run import input
    from mfc.test.cases import list_cases

    return get_target, input, list_cases


def _materialize(builder, input_module):
    """Build a slug-able case, muting the case files' own chatter (they print diagnostics).

    `get_slug` resolves a chemistry case's mechanism through `get_cantera_solution`, which
    lives on MFCInputFile rather than on the TestCase the suite hands out, so the params go
    back through MFCInputFile -- rooted at the test's own directory, where the run would
    look for a mechanism shipped beside the case.
    """
    with contextlib.redirect_stdout(io.StringIO()):
        case = builder.to_case()
        return input_module.MFCInputFile("case.py", case.get_dirpath(), case.params)


def _golden_cases(list_cases, input_module):
    """(builder, case) for every golden test, skipping the self-driving convergence runs."""
    for builder in list_cases():
        if builder.kind != "golden":
            continue
        yield builder, _materialize(builder, input_module)


def _slugs(case, targets) -> dict:
    return {target.name: target.get_slug(case) for target in targets}


def _allowed_slugs(cases, targets, input_module) -> dict:
    """Per target, the slugs the split pre-build leaves on disk.

    `./mfc.sh build` with no case file builds `input.load(None, [], {})` -- an empty
    case -- and `--only Chemistry` / `--only eos=mie_gruneisen` select on exact trace
    segments (see _filter_only).
    """
    default_case = input_module.load(None, [], {})
    allowed = {name: {slug} for name, slug in _slugs(default_case, targets).items()}

    for builder, case in cases:
        segments = builder.trace.split(" -> ")
        if CHEMISTRY_LABEL not in segments and EOS_LABEL not in segments:
            continue
        for name, slug in _slugs(case, targets).items():
            allowed[name].add(slug)

    return allowed


@contextlib.contextmanager
def _chdir(path: Path):
    """cases.py names example case files relative to the repo root, and get_slug builds
    paths from the cwd. pytest runs from toolchain/, so anchor both here."""
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def check_no_case_specific_builds(repo_root: Path) -> list:
    get_target, input_module, list_cases = _import_toolchain(repo_root)

    with _chdir(repo_root):
        targets = [get_target(name) for name in TARGET_NAMES]
        cases = list(_golden_cases(list_cases, input_module))
        allowed = _allowed_slugs(cases, targets, input_module)

        errors = []
        for builder, case in cases:
            unbuilt = sorted(name for name, slug in _slugs(case, targets).items() if slug not in allowed[name])
            if unbuilt:
                errors.append(f"  {builder.get_uuid()}  {builder.trace}\n      needs its own build of: {', '.join(unbuilt)}")

    return errors


def main():
    repo_root = Path(__file__).resolve().parents[2]

    errors = check_no_case_specific_builds(repo_root)

    if errors:
        print("Test suite check failed: these cases need a build no CI lane pre-builds.")
        print("Every lane but Frontier AMD gpu-omp would build them on demand; that one")
        print("hands the test a binary that does not exist. See this file's docstring.")
        print("")
        for e in errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
