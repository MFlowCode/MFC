# typed: strict
# frozen_string_literal: true

# Homebrew formula for MFC (Multiphase Flow Code)
# This formula is automatically deployed to the MFlowCode/homebrew-mfc tap by GitHub Actions workflow
class Mfc < Formula
  desc "Exascale multiphase/multiphysics compressible flow solver"
  homepage "https://mflowcode.github.io/"
  url "https://github.com/MFlowCode/MFC/archive/refs/tags/v5.1.0.tar.gz"
  sha256 "4684bee6a529287f243f8929fb7edb0dfebbb04df7c1806459761c9a6c9261cf"
  license "MIT"
  head "https://github.com/MFlowCode/MFC.git", branch: "master"

  depends_on "cmake" => :build
  depends_on "gcc" => :build

  depends_on "boost"
  depends_on "fftw"
  depends_on "hdf5"
  depends_on "open-mpi"
  depends_on "openblas"
  depends_on "python@3.12"

  # Disable bottles due to Python venv with compiled extensions that can't be relocated
  def pour_bottle?
    false
  end

  def install
    # Create Python virtual environment inside libexec (inside Cellar for proper bottling)
    venv = libexec/"venv"
    system Formula["python@3.12"].opt_bin/"python3.12", "-m", "venv", venv
    system venv/"bin/pip", "install", "--upgrade", "pip", "setuptools", "wheel"

    # Install Cantera from PyPI (required dependency for MFC build)
    system venv/"bin/pip", "install", "cantera==3.1.0"

    # Install MFC Python package and dependencies into venv
    # Keep toolchain in buildpath for now - mfc.sh needs it there
    # Use editable install (-e) to avoid RECORD file issues when venv is symlinked at runtime
    system venv/"bin/pip", "install", "-e", buildpath/"toolchain"

    # Create symlink so mfc.sh uses our pre-installed venv
    mkdir_p "build"
    ln_sf venv, "build/venv"

    # Now build MFC with pre-configured venv
    # Set VIRTUAL_ENV so mfc.sh uses existing venv instead of creating new one
    ENV["VIRTUAL_ENV"] = venv

    # Build MFC using pre-configured venv
    # Must run from buildpath (MFC root directory) where toolchain/ exists
    Dir.chdir(buildpath) do
      system "./mfc.sh", "build", "-t", "pre_process", "simulation", "post_process", "-j", ENV.make_jobs.to_s
    end

    # After build completes, install Python toolchain to prefix
    prefix.install "toolchain"

    # Install only executable files from hashed build dirs
    Dir.glob("build/install/*/bin/*").each do |p|
      next if !File.file?(p) || !File.executable?(p)

      bin.install p
      (bin/File.basename(p)).chmod 0755
    end

    # Install main mfc.sh script
    libexec.install "mfc.sh"

    # Install examples
    prefix.install "examples"

    # Create smart wrapper script that:
    # 1. Works around read-only Cellar issue by using a temp working dir
    # 2. Uses the persistent var-based venv (no copy needed)
    # 3. Minimal copying - only what needs to be writable
    # 4. Resolves input file paths before changing directories
    (bin/"mfc").write <<~EOS
        #!/bin/bash
        set -euo pipefail

        # Unset VIRTUAL_ENV to ensure mfc.sh uses our configured venv
        unset VIRTUAL_ENV || true

        # Save original working directory
        ORIG_DIR="$(pwd)"

        # Process arguments and convert relative paths to absolute paths
        ARGS=()
        for arg in "$@"; do
          # If argument looks like a file path and exists as relative path
          if [[ "$arg" =~ \\.(py|txt|json|yaml|yml)$ ]] && [ -e "${ORIG_DIR}/${arg}" ]; then
            ARGS+=("${ORIG_DIR}/${arg}")
          else
            ARGS+=("$arg")
          fi
        done

        # Create a temporary working directory (Cellar is read-only)
        TMPDIR="$(mktemp -d)"
        trap 'rm -rf "${TMPDIR}"' EXIT
        cd "${TMPDIR}"

        # Copy only mfc.sh (small, fast)
        cp "#{libexec}/mfc.sh" .

        # Create a minimal CMakeLists.txt to prevent the cat error
        # mfc.sh reads this to get version info
        cat > CMakeLists.txt << 'CMAKE_EOF'
      cmake_minimum_required(VERSION 3.18)
      project(MFC VERSION #{version})
      CMAKE_EOF

        # Symlink toolchain (read-only is fine, no copy needed)
        ln -s "#{prefix}/toolchain" toolchain

        # Symlink examples (read-only is fine)
        ln -s "#{prefix}/examples" examples

        # Create build directory structure
        mkdir -p build

        # Symlink the persistent venv (no copy)
        ln -s "#{libexec}/venv" build/venv

        # Copy only pyproject.toml (tiny file, prevents reinstall checks)
        cp "#{prefix}/toolchain/pyproject.toml" build/pyproject.toml

        # Create a sitecustomize.py file that patches MFC paths before any imports
        # This runs automatically at Python startup and ensures paths are correct
        cat > sitecustomize.py << 'SITECUSTOMIZE_EOF'
      import sys
      import os

      # Add toolchain to path
      sys.path.insert(0, "#{prefix}/toolchain")

      # Patch MFC paths before mfc.common is imported anywhere
      _mfc_temp_root = os.getcwd()

      def _patch_mfc_common():
          """Patches mfc.common module to use temp directory for writable files."""
          import mfc.common
          mfc.common.MFC_ROOT_DIR = _mfc_temp_root
          mfc.common.MFC_BUILD_DIR = os.path.join(_mfc_temp_root, "build")
          mfc.common.MFC_LOCK_FILEPATH = os.path.join(mfc.common.MFC_BUILD_DIR, "lock.yaml")
          # Keep toolchain and examples pointing to Homebrew installation
          mfc.common.MFC_TOOLCHAIN_DIR = "#{prefix}/toolchain"
          mfc.common.MFC_EXAMPLE_DIRPATH = "#{prefix}/examples"

      def _patch_mfc_build():
          """Patches MFCTarget to use pre-installed binaries."""
          from mfc.build import MFCTarget

          # Override get_install_binpath to use pre-installed binaries
          _original_get_install_binpath = MFCTarget.get_install_binpath
          def _homebrew_get_install_binpath(self, case):
              return "#{bin}/" + self.name
          MFCTarget.get_install_binpath = _homebrew_get_install_binpath

          # Override is_buildable to skip building main targets
          _original_is_buildable = MFCTarget.is_buildable
          def _homebrew_is_buildable(self):
              if self.name in ["pre_process", "simulation", "post_process", "syscheck"]:
                  return False
              return _original_is_buildable(self)
          MFCTarget.is_buildable = _homebrew_is_buildable

      # Apply patches immediately
      _patch_mfc_common()
      _patch_mfc_build()
      SITECUSTOMIZE_EOF

        # Set PYTHONPATH to include current directory so sitecustomize.py is found
        export PYTHONPATH="${TMPDIR}:#{prefix}/toolchain:${PYTHONPATH:-}"

        # For 'mfc run', add --no-build flag to skip compilation
        if [ "${1-}" = "run" ]; then
          exec ./mfc.sh "${ARGS[@]}" --no-build
        else
          exec ./mfc.sh "${ARGS[@]}"
        fi
    EOS
    (bin/"mfc").chmod 0755
  end

  def post_install
    # Fix executable permissions (Homebrew sometimes overrides them)
    (bin/"mfc").chmod 0755
  end

  def caveats
    <<~EOS
      MFC has been installed successfully!

      To run a case:
        mfc run <case.py>

      Pre-built binaries are also available directly:
        pre_process, simulation, post_process

      Examples are available in:
        #{prefix}/examples

      Example:
        cp #{prefix}/examples/1D_sodshocktube/case.py .
        mfc run case.py

      Note: Cantera 3.1.0 is pre-installed in the MFC virtual environment.
    EOS
  end

  test do
    # Test that all binaries exist and are executable
    %w[pre_process simulation post_process].each do |prog|
      assert_path_exists bin/prog
      assert_predicate bin/prog, :executable?
    end

    # Test that toolchain is installed
    assert_path_exists prefix/"toolchain"

    # Test that venv exists and has required packages
    assert_path_exists libexec/"venv"
    assert_predicate (libexec/"venv/bin/python"), :executable?

    # Test that examples exist
    assert_path_exists prefix/"examples"

    # Test that mfc wrapper works
    system bin/"mfc", "--help"

    # Test running a simple 1D Sod shock tube case from a separate directory
    # This ensures the wrapper script correctly handles relative paths
    testpath_case = testpath/"test_run"
    testpath_case.mkpath

    # Copy case.py from examples to an independent test directory
    cp prefix/"examples/1D_sodshocktube/case.py", testpath_case/"case.py"

    # Run the case from the test directory (this will execute pre_process and simulation)
    # Limit to 1 processor and reduce runtime for testing
    cd testpath_case do
      system bin/"mfc", "run", "case.py", "-j", "1"
    end

    # Verify output files were created in the test directory
    assert_path_exists testpath_case/"silo_hdf5"
    assert_predicate testpath_case/"silo_hdf5", :directory?
  end
end
