"""Toolchain hooks for SCC fingerprints (see core.py) and the `./mfc.sh fingerprint` precheck."""

import os
from typing import Optional

from ..common import MFC_ROOT_DIR, MFCException
from ..printer import cons
from ..state import ARG
from . import core

_tree = None
_reported = False


def _tree_state() -> dict:
    global _tree  # noqa: PLW0603
    if _tree is None:
        _tree = core.tree_state(MFC_ROOT_DIR, core.read_rules(os.path.join(MFC_ROOT_DIR, core.PUBLIC_RULES_PATH)))
    return _tree


def _pubkey() -> bytes:
    try:
        return core.load_pubkey(MFC_ROOT_DIR)
    except (OSError, ValueError) as exc:
        raise MFCException(f"Fingerprint: cannot read the organizer key {core.PUBKEY_PATH} ({exc}).") from exc


def _encode(stage: str, target: str, exe_sha: str, record: dict) -> str:
    global _reported  # noqa: PLW0603
    tree = _tree_state()
    if not _reported:
        _reported = True
        if tree["baseline_sha"] is None:
            cons.print(f"[yellow]Fingerprint: no {core.BASELINE_PATH}; recording the tree without a baseline.[/yellow]")
        elif tree["violations"]:
            cons.print(f"[yellow]Fingerprint: {len(tree['violations'])} change(s) outside the allowed list. Run ./mfc.sh fingerprint for details.[/yellow]")
    try:
        return core.encode_line(_pubkey(), stage, target, exe_sha, {"stage": stage, "tree": tree, **record, **core.context()})
    except ImportError as exc:
        raise MFCException("Fingerprint: the 'cryptography' package is missing; rerun ./mfc.sh to refresh the Python environment.") from exc


def _find_line(log: str, *header: str) -> Optional[str]:
    """Last build record in `log` whose header has `header` from the exe_sha field on."""
    if not os.path.isfile(log):
        return None
    with open(log) as f:
        hits = [line for line in f if line.split(" ")[4 : 4 + len(header)] == list(header)]
    return hits[-1] if hits else None


def _append(path: str, lines) -> None:
    with open(path, "a") as f:
        f.writelines(lines)


def record_build(target, case) -> None:
    """Log the tree that produced the installed binary, once per (binary, tree) pair."""
    exe = core.file_sha(target.get_install_binpath(case))
    log = os.path.join(target.get_install_dirpath(case), core.BUILD_LOG_NAME)
    if exe != "-" and _find_line(log, exe, _tree_state()["digest"]) is None:
        _append(log, [_encode("build", target.name, exe, {"target": target.name, "exe_sha": exe})])


def record_run(targets, case) -> None:
    """Append the matching build records and a run record to the case directory's log."""
    log = os.path.join(os.path.dirname(os.path.abspath(ARG("input"))), core.LOG_NAME)
    existing = ""
    if os.path.isfile(log):
        with open(log) as f:
            existing = f.read()

    exes, copied = {}, []
    for t in (t for t in targets if not t.isDependency):
        exe = core.file_sha(t.get_install_binpath(case))
        line = _find_line(os.path.join(t.get_install_dirpath(case), core.BUILD_LOG_NAME), exe)
        exes[t.name] = {"exe_sha": exe, "build_record": line is not None}
        if line and line not in existing:
            copied.append(line)
    _append(log, copied)

    with open(ARG("input"), "rb") as f:
        case_src = f.read().decode(errors="replace")
    record = {"exes": exes, "prev": core.file_sha(log), "case": {"path": os.path.abspath(ARG("input")), "source": case_src, "params": {k: str(v) for k, v in case.params.items()}}}
    _append(log, [_encode("run", "-", "-", record)])


def _diff_stat(change: dict, base_sha: Optional[str]) -> str:
    """'+added -removed' against the baseline version, recovered from git history by hash."""
    new = b""
    if change["status"] != "D":
        with open(os.path.join(MFC_ROOT_DIR, change["path"]), "rb") as f:
            new = f.read()
    old = b"" if base_sha is None else core.baseline_content(MFC_ROOT_DIR, change["path"], base_sha)
    if old is None:
        return "(no baseline)"
    stat = core.diff_stat(old, new)
    return "(binary)" if stat is None else f"+{stat[0]} -{stat[1]}"


def fingerprint() -> None:
    """`./mfc.sh fingerprint`: check the tree against the public rules."""
    tree = _tree_state()
    if tree["baseline_sha"] is None:
        raise MFCException(f"Fingerprint: {core.BASELINE_PATH} is missing, so changes cannot be determined.")
    _encode("check", "-", "-", {})  # proves the key and the crypto library work

    cons.print(f"[bold]Fingerprint precheck[/bold] against [magenta]{core.PUBLIC_RULES_PATH}[/magenta] ({tree['nfiles']} files)")
    cons.indent()
    violations = set(tree["violations"])
    baseline = core.load_baseline(os.path.join(MFC_ROOT_DIR, core.BASELINE_PATH))
    width = max((len(c["path"]) for c in tree["changes"]), default=0)
    for c in tree["changes"]:
        if c["path"] in violations:
            verdict, stat = "[red]not allowed[/red]", ""
        elif core.matches(c["path"], core.GENERATED):
            verdict, stat = "[dim]generated[/dim]", ""
        else:
            verdict, stat = "[green]allowed[/green]", _diff_stat(c, baseline.get(c["path"]))
        cons.print(f"{c['status']} {c['path']:<{width}}  {stat:<14} {verdict}", highlight=False)
    if not tree["changes"]:
        cons.print("No changes from the baseline.")
    cons.unindent()

    if violations:
        raise MFCException(f"Fingerprint precheck FAILED: {len(violations)} file(s) changed outside the allowed list.")
    cons.print("[bold green]Fingerprint precheck PASSED.[/bold green]")
