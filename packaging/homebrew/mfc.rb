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

  def install
    # Create Python virtual environment (remove existing one first for clean reinstalls)
    venv = var/"mfc/venv"
    mkdir_p venv.parent
    rm_r(venv, force: true)
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
    # 3. Ensures mfc.sh doesn't reinstall packages by copying pyproject.toml
    (bin/"mfc").write <<~EOS
      #!/bin/bash
      set -euo pipefail

      # Unset VIRTUAL_ENV to ensure mfc.sh uses our configured venv
      unset VIRTUAL_ENV || true

      # Create a temporary working directory (Cellar is read-only)
      TMPDIR="$(mktemp -d)"
      trap 'rm -rf "${TMPDIR}"' EXIT

      # Copy mfc.sh to temp dir (it may try to write build artifacts)
      cp "#{libexec}/mfc.sh" "${TMPDIR}/"
      cd "${TMPDIR}"

      # Copy toolchain directory (not symlink) so Python paths resolve correctly
      # This prevents paths from resolving back to read-only Cellar
      cp -R "#{prefix}/toolchain" "toolchain"

      # Patch toolchain to use Homebrew-installed binaries
      # Replace get_install_binpath to return Homebrew bin directory
      cat >> "toolchain/mfc/build.py" << 'PATCH_EOF'

      # Homebrew patch: Override get_install_binpath to use pre-installed binaries
      _original_get_install_binpath = MFCTarget.get_install_binpath
      def _homebrew_get_install_binpath(self, case):
          return "#{bin}/" + self.name
      MFCTarget.get_install_binpath = _homebrew_get_install_binpath

      # Override is_buildable to skip building main targets and syscheck
      _original_is_buildable = MFCTarget.is_buildable
      def _homebrew_is_buildable(self):
          if self.name in ["pre_process", "simulation", "post_process", "syscheck"]:
              return False  # Skip building - use pre-installed binaries
          return _original_is_buildable(self)
      MFCTarget.is_buildable = _homebrew_is_buildable
      PATCH_EOF

      # Copy examples directory (required by mfc.sh Python code)
      cp -R "#{prefix}/examples" "examples"

      # Create build directory and symlink the persistent venv
      mkdir -p "build"
      ln -s "#{var}/mfc/venv" "build/venv"

      # Copy pyproject.toml so mfc.sh thinks dependencies are already installed
      cp "#{prefix}/toolchain/pyproject.toml" "build/pyproject.toml"

      # For 'mfc run', add --no-build flag to skip compilation
      if [ "${1-}" = "run" ]; then
        exec ./mfc.sh "$@" --no-build
      else
        exec ./mfc.sh "$@"
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
    assert_path_exists var/"mfc/venv"
    assert_predicate (var/"mfc/venv/bin/python"), :executable?

    # Test that examples exist
    assert_path_exists prefix/"examples"

    # Test that mfc wrapper works
    system bin/"mfc", "--help"
  end
end
