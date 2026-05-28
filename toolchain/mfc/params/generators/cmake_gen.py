"""
Generate Fortran parameter .fpp files into the CMake build directory.

Called by CMakeLists.txt at configure time:
  python3 cmake_gen.py <cmake_binary_dir>

Writes generated_namelist.fpp and generated_decls.fpp into
  <cmake_binary_dir>/include/{pre_process,simulation,post_process}/
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from mfc.params.generators.fortran_gen import generate_decls_fpp, generate_namelist_fpp

_TARGETS = [
    ("pre", "pre_process"),
    ("sim", "simulation"),
    ("post", "post_process"),
]


def main(build_dir: Path) -> None:
    for short, full in _TARGETS:
        out_dir = build_dir / "include" / full
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "generated_namelist.fpp").write_text(generate_namelist_fpp(short))
        (out_dir / "generated_decls.fpp").write_text(generate_decls_fpp(short))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <cmake_binary_dir>")
    main(Path(sys.argv[1]))
