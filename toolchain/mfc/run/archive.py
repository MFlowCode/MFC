import dataclasses
import datetime
import os
import shutil
import sys
import tarfile
import tempfile

from ..common import MFCException, does_command_exist, file_dump_yaml, generate_git_tagline, system
from ..printer import cons
from ..state import ARG, CFG
from . import input

ARTIFACT_FILENAMES = [
    "equations.dat",
    "run_time.inf",
    "time_data.dat",
    "io_time_data.dat",
    "fort.1",
    "pre_time_data.dat",
]

ARTIFACT_DIRNAMES = [
    "D",
    "p_all",
    "restart_data",
    "silo_hdf5",
]


def __collect_sources(case: input.MFCInputFile, targets) -> list:
    dirpath = case.dirpath
    name = ARG("name")
    sources = []

    if os.path.isfile(case.filename):
        sources.append(case.filename)

    for target in targets:
        inp = os.path.join(dirpath, f"{target.name}.inp")
        if os.path.isfile(inp):
            sources.append(inp)

    for candidate in ARTIFACT_FILENAMES:
        candidate_path = os.path.join(dirpath, candidate)
        if os.path.isfile(candidate_path):
            sources.append(candidate_path)

    for script_ext in ("sh", "bat"):
        script = os.path.join(dirpath, f"{name}.{script_ext}")
        if os.path.isfile(script):
            sources.append(script)

    out_file = os.path.join(dirpath, f"{name}.out")
    if os.path.isfile(out_file):
        sources.append(out_file)

    summary = ARG("output_summary")
    if summary is not None and os.path.isfile(summary):
        summary_abs = os.path.abspath(summary)
        dirpath_abs = os.path.abspath(dirpath)
        if os.path.commonpath([summary_abs, dirpath_abs]) == dirpath_abs:
            sources.append(summary_abs)
        else:
            cons.print(f"[yellow]Archive: skipping output_summary outside case dir: {summary_abs}[/yellow]")

    for dirname in ARTIFACT_DIRNAMES:
        candidate_dir = os.path.join(dirpath, dirname)
        if os.path.isdir(candidate_dir):
            sources.append(candidate_dir)

    return sources


def __build_manifest(case: input.MFCInputFile, targets, sources: list, archive_path: str, archive_format: str) -> dict:
    dirpath = case.dirpath
    relative_sources = []
    for src in sources:
        try:
            rel = os.path.relpath(src, dirpath)
        except ValueError:
            rel = src
        relative_sources.append(rel)

    return {
        "timestamp": datetime.datetime.now().astimezone().isoformat(),
        "source_case": os.path.abspath(case.filename),
        "source_dir": os.path.abspath(dirpath),
        "invocation": sys.argv[1:],
        "git": generate_git_tagline(),
        "targets": [t.name for t in targets],
        "archive_format": archive_format,
        "archive_path": archive_path,
        "build_lock": dataclasses.asdict(CFG()),
        "contents": sorted(relative_sources),
    }


def __copy_dir(sources: list, case: input.MFCInputFile, dest: str) -> None:
    os.makedirs(dest, exist_ok=True)
    dirpath = case.dirpath

    for src in sources:
        try:
            rel = os.path.relpath(src, dirpath)
        except ValueError:
            rel = os.path.basename(src)

        target_path = os.path.join(dest, rel)
        os.makedirs(os.path.dirname(target_path), exist_ok=True)

        if os.path.isdir(src):
            shutil.copytree(src, target_path, dirs_exist_ok=True, symlinks=True)
        else:
            shutil.copy2(src, target_path)


def __write_tar(sources: list, case: input.MFCInputFile, dest: str, compressed: bool, manifest_file: str) -> None:
    dirpath = case.dirpath
    arcroot = os.path.basename(dest).removesuffix(".tar.zst").removesuffix(".tar")

    def rel_for(path: str) -> str:
        try:
            return os.path.relpath(path, dirpath)
        except ValueError:
            return os.path.basename(path)

    if compressed:
        if not does_command_exist("tar"):
            raise MFCException("Archive: 'tar' binary not found; required for --archive-format tar.zst.")

        with tempfile.TemporaryDirectory() as staging:
            staging_root = os.path.join(staging, arcroot)
            os.makedirs(staging_root, exist_ok=True)

            for src in sources:
                rel = rel_for(src)
                if rel.startswith(".."):
                    raise MFCException(f"Archive: refusing to include source outside case dir: {src}")
                target = os.path.join(staging_root, rel)
                os.makedirs(os.path.dirname(target), exist_ok=True)
                if os.path.isdir(src):
                    shutil.copytree(src, target, symlinks=True)
                else:
                    shutil.copy2(src, target)

            shutil.copy2(manifest_file, os.path.join(staging_root, "manifest.yaml"))

            result = system(
                ["tar", "--zstd", "-cf", dest, "-C", staging, arcroot],
                print_cmd=False,
            )
            if result.returncode != 0:
                raise MFCException(f"Archive: 'tar --zstd' failed with exit code {result.returncode}. Ensure GNU tar >= 1.31.")
        return

    with tarfile.open(dest, "w") as tf:
        for src in sources:
            rel = rel_for(src)
            if rel.startswith(".."):
                raise MFCException(f"Archive: refusing to include source outside case dir: {src}")
            tf.add(src, arcname=os.path.join(arcroot, rel))
        tf.add(manifest_file, arcname=os.path.join(arcroot, "manifest.yaml"))


def archive(case: input.MFCInputFile, targets) -> None:
    archive_root = ARG("archive")
    if archive_root is None:
        return

    archive_format = ARG("archive_format") or "dir"
    name = ARG("name")
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    stem = f"{name}-{timestamp}"

    archive_root = os.path.abspath(os.path.expanduser(archive_root))
    os.makedirs(archive_root, exist_ok=True)

    suffix_map = {"dir": "", "tar": ".tar", "tar.zst": ".tar.zst"}
    if archive_format not in suffix_map:
        raise MFCException(f"Archive: unsupported format '{archive_format}'. Must be one of: {', '.join(suffix_map)}.")
    dest = os.path.join(archive_root, stem + suffix_map[archive_format])

    if os.path.exists(dest):
        raise MFCException(f"Archive: destination already exists: {dest}")

    sources = __collect_sources(case, targets)
    if not sources:
        cons.print("[yellow]Archive: no artifacts found to archive; skipping.[/yellow]")
        return

    cons.print()
    cons.print(f"[bold]Archiving[/bold] to [magenta]{dest}[/magenta] ({archive_format})")

    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as tmp:
        manifest_path = tmp.name

    cons.indent()
    try:
        manifest = __build_manifest(case, targets, sources, dest, archive_format)
        file_dump_yaml(manifest_path, manifest)

        if archive_format == "dir":
            __copy_dir(sources, case, dest)
            shutil.copy2(manifest_path, os.path.join(dest, "manifest.yaml"))
        else:
            __write_tar(sources, case, dest, compressed=(archive_format == "tar.zst"), manifest_file=manifest_path)

        cons.print(f"Wrote [magenta]{len(sources)}[/magenta] artifact(s) + manifest.yaml.")
    finally:
        cons.unindent()
        if os.path.isfile(manifest_path):
            os.unlink(manifest_path)
