"""Fypp-expand MFC sources the way the build does, without a build. Used by `./mfc.sh context`."""

import hashlib
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# Target-independent arguments of the Fypp command in cmake/Fypp.cmake; test_context checks they match.
FYPP_FLAGS = ["-m", "re", "--line-numbering", "--no-folding", "--line-length=999", "--line-numbering-mode=nocontlines"]


def _fypp_exe() -> str:
    """The fypp beside the running interpreter (the toolchain venv), else the one on PATH."""
    beside = Path(sys.executable).parent / "fypp"
    return str(beside) if beside.exists() else (shutil.which("fypp") or "fypp")


def fypp_command(fypp, repo, include_dir, target, src, out, compiler="NVHPC"):
    common = repo / "src" / "common"
    incs = [include_dir / target, repo / "src" / target / "include", common / "include", common, common / "include" / "defaults"]
    defines = [f"MFC_{compiler}", f"MFC_{target.upper()}", f'MFC_COMPILER="{compiler}"', "MFC_CASE_OPTIMIZATION=False", "chemistry=False"]
    return [fypp, *FYPP_FLAGS[:2], *(a for i in incs for a in ("-I", str(i))), *(a for d in defines for a in ("-D", d)), *FYPP_FLAGS[2:], str(src), str(out)]


def _write_includes(repo: Path, include_dir: Path, targets):
    """Write what the build generates before Fypp runs: parameter includes and a default case.fpp."""
    sys.path.insert(0, str(repo / "toolchain"))
    from mfc import state  # pylint: disable=import-outside-toplevel
    from mfc.case import Case  # pylint: disable=import-outside-toplevel
    from mfc.params.generators.fortran_gen import get_generated_files  # pylint: disable=import-outside-toplevel

    for path, content in get_generated_files(include_dir.parent):
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists() or path.read_text() != content:
            path.write_text(content)
    state.gARG.setdefault("case_optimization", False)  # get_fpp reads it; nothing is case-optimized here
    for t in targets:
        (include_dir / t).mkdir(parents=True, exist_ok=True)
        case = include_dir / t / "case.fpp"
        content = Case({}).get_fpp(t, print=False)
        if not case.exists() or case.read_text() != content:
            case.write_text(content)


def _sources(repo: Path, target: str):
    return sorted((repo / "src" / "common").glob("*.fpp")) + sorted((repo / "src" / target).glob("*.fpp"))


def _key(paths) -> str:
    h = hashlib.sha256()
    for p in paths:
        h.update(str(p).encode())
        h.update(p.read_bytes())
    return h.hexdigest()[:16]


def expand(repo, targets, cache, compiler="NVHPC", write_includes=True):
    """Expand each target's sources into cache/<target>/<name>.f90. Returns (files per target, errors)."""
    repo, cache = Path(repo).resolve(), Path(cache).resolve()
    include_dir = cache / "include"
    if write_includes:
        _write_includes(repo, include_dir, targets)
    fypp = _fypp_exe()
    inc_key = _key(sorted(repo.glob("src/*/include/**/*.fpp")) + sorted(include_dir.glob("*/*.fpp")))
    jobs = []
    for t in targets:
        for src in _sources(repo, t):
            # Keyed by source and include content, so unchanged sources are not re-expanded.
            jobs.append((t, src, cache / t / f"{src.name}.{_key([src])}{inc_key}.f90", cache / t / f"{src.name}.f90"))

    def run(job):
        t, src, keyed, final = job
        if not keyed.exists():
            keyed.parent.mkdir(parents=True, exist_ok=True)
            r = subprocess.run(fypp_command(fypp, repo, include_dir, t, src, keyed, compiler), capture_output=True, text=True, check=False)
            if r.returncode:
                keyed.unlink(missing_ok=True)
                final.unlink(missing_ok=True)  # never index the previous expansion against new line numbers
                return f"could not expand {src.relative_to(repo)}: {(r.stderr.strip().splitlines() or ['?'])[-1]}"
        if not final.is_symlink() or final.resolve() != keyed:
            final.unlink(missing_ok=True)
            final.symlink_to(keyed.name)
        return None

    with ThreadPoolExecutor(16) as ex:
        errors = [e for e in ex.map(run, jobs) if e]
    files = {t: [f for tt, _, _, f in jobs if tt == t and f.exists()] for t in targets}
    return files, errors
