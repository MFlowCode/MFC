"""Generate Fortran parameter .fpp files into the CMake build directory.

Called by CMakeLists.txt at configure time:
  python3 cmake_gen.py <cmake_binary_dir>
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from mfc.params.generators.fortran_gen import generate_decls_fpp, generate_namelist_fpp  # noqa: E402

if len(sys.argv) != 2:
    sys.exit(f"Usage: {sys.argv[0]} <cmake_binary_dir>")

build_dir = Path(sys.argv[1])
for short, full in [("pre", "pre_process"), ("sim", "simulation"), ("post", "post_process")]:
    out_dir = build_dir / "include" / full
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "generated_namelist.fpp").write_text(generate_namelist_fpp(short))
    (out_dir / "generated_decls.fpp").write_text(generate_decls_fpp(short))
