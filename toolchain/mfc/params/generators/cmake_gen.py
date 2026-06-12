"""Generate Fortran parameter .fpp files into the CMake build directory.

Invoked by the build-time custom command in cmake/ParamsCodegen.cmake:
  python3 cmake_gen.py <cmake_binary_dir>
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from mfc.params.generators.fortran_gen import get_generated_files  # noqa: E402

if len(sys.argv) != 2:
    sys.exit(f"Usage: {sys.argv[0]} <cmake_binary_dir>")

build_dir = Path(sys.argv[1])
for path, content in get_generated_files(build_dir):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
