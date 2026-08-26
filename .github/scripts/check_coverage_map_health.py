"""Fail loudly if the committed coverage map is stale or under-covers. Used by coverage-health.yml."""
import datetime
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "toolchain"))
from mfc.test.coverage import COVERAGE_MAP_PATH, load_map, map_health  # noqa: E402
from mfc.test.cases import list_cases  # noqa: E402  (returns the current test list)

MAX_AGE_DAYS = 10
MIN_FRACTION = 0.80

# coverage-refresh.yml moves this ref to the commit it just rebuilt the map from, on EVERY
# successful refresh -- including the ones whose entries come out identical and therefore
# push no commit at all.
#
# The map's own _meta.git_sha cannot answer the freshness question. It only advances when a
# commit lands, and the refresh deliberately skips the commit when the coverage entries are
# unchanged (a no-op commit on every push is what #1683 removed). So a coverage-relevant
# commit whose coverage happens to be identical pinned git_sha in the past and made this
# check fail every day thereafter, with no way to recover: a refresh that WAS working
# looked broken. The question here is "has a refresh run since the last relevant commit",
# and only the refresh itself can answer it.
#
# It lives outside refs/heads/ and refs/tags/ so that updating it fires no `push` workflow
# trigger. actions/checkout does not fetch it; coverage-health.yml fetches it explicitly.
VERIFIED_REF = "refs/coverage-map/verified"

# Must mirror the `paths:` trigger of coverage-refresh.yml. A superset would raise false
# alarms (a change that never triggers a refresh would look like a missed refresh).
COVERAGE_RELEVANT_PATHS = [":(glob)src/**/*.fpp", "toolchain/mfc/test/cases.py"]


def _git_env():
    """The environment minus every GIT_* variable.

    `cwd` is how callers select the repository here, but git honors an inherited GIT_DIR /
    GIT_INDEX_FILE over the working directory -- and git exports both while running a hook.
    MFC's pre-commit hook runs precheck, whose toolchain lint imports this script and calls
    these functions against throwaway test repos; without the scrub those calls silently
    answer for the developer's real repository instead.
    """
    return {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}


def verified_sha(cwd=None):
    """Commit the last successful refresh verified the map against, or None if unknown.

    None covers both "the ref was never pushed" (a fork, or the window before the first
    refresh after VERIFIED_REF was introduced) and "the fetch did not bring it down". The
    caller must read that as undeterminable and fall back to the wall-clock age rule, not
    as a failure -- an absent ref is not evidence of a broken refresh.
    """
    rev = subprocess.run(["git", "rev-parse", "--verify", "--quiet", f"{VERIFIED_REF}^{{commit}}"], capture_output=True, text=True, check=False, cwd=cwd, env=_git_env())
    return rev.stdout.strip() or None


def verified_after_last_change(git_sha, cwd=None):
    """Did that refresh run at or after the last coverage-relevant commit? None if unknowable.

    This is the direct question the wall-clock age check only approximated: it stays quiet
    over a genuinely quiet repo (no refresh needed, so no heartbeat needed either) and
    fires on the next relevant push after a refresh breaks, rather than 10 days later.
    """
    if not git_sha:
        return None
    last = subprocess.run(["git", "log", "-1", "--format=%H", "--", *COVERAGE_RELEVANT_PATHS], capture_output=True, text=True, check=False, cwd=cwd, env=_git_env())
    if last.returncode != 0 or not last.stdout.strip():
        return None  # shallow clone or no such commit -> fall back to the age rule
    ancestor = subprocess.run(["git", "merge-base", "--is-ancestor", last.stdout.strip(), git_sha], capture_output=True, check=False, cwd=cwd, env=_git_env())
    return {0: True, 1: False}.get(ancestor.returncode)  # anything else -> None (unknown sha, shallow history)


def main():
    entries, meta = load_map(COVERAGE_MAP_PATH)
    if entries is None:
        sys.exit("Coverage map missing or corrupt.")
    # Compute each current test's coverage key. Loading a case runs its case file as a
    # subprocess, so anything that file imports must resolve. Keep skipping the ones that
    # do not rather than crashing -- map_health measures the fraction of *loadable*
    # current tests that are mapped, so a smaller current_keys cannot produce a false
    # "stale" result -- but the skip is a safety net, not the expected path: every case
    # here is meant to load. The 16 chemistry cases that used to land in this list did so
    # because get_py_program_output ran them under PATH's python3 while this job invokes
    # build/venv/bin/python3 directly, leaving the venv's cantera out of reach.
    current_keys = set()
    unloadable = []
    for b in list_cases():
        try:
            current_keys.add(b.to_case().coverage_key())
        except Exception as exc:  # noqa: BLE001 -- a case file that won't import must not crash the health check
            last_line = (str(exc).strip().splitlines() or ["(no message)"])[-1][:140]
            unloadable.append((getattr(b, "trace", repr(b)), last_line))
    if unloadable:
        print(f"Note: {len(unloadable)} case(s) could not be loaded in this lightweight job (excluded from the check):")
        for trace, err in unloadable[:15]:
            print(f"  - {trace}: {err}")
    ok, msg = map_health(
        meta=meta,
        current_keys=current_keys,
        mapped_keys=set(entries),
        now=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        max_age_days=MAX_AGE_DAYS,
        min_fraction=MIN_FRACTION,
        verified_after_last_change=verified_after_last_change(verified_sha()),
    )
    print(msg)
    return 0 if ok else 1


# Guarded so the two git predicates above can be imported and unit-tested against a
# throwaway repository; running the module still performs the full check.
if __name__ == "__main__":
    sys.exit(main())
