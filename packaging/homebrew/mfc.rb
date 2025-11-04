# typed: strict
# frozen_string_literal: true

# Homebrew formula for MFC (Multiphase Flow Code)
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
  depends_on "sundials"
  depends_on "yaml-cpp"

  resource "cantera" do
    url "https://github.com/Cantera/cantera.git",
        tag: "v3.1.0"
  end

  def install
    # Create Python virtual environment first (before MFC build)
    venv = libexec/"venv"
    system Formula["python@3.12"].opt_bin/"python3.12", "-m", "venv", venv
    system venv/"bin/pip", "install", "--upgrade", "pip", "setuptools", "wheel"

    # Build and install Cantera 3.1.0 from source BEFORE MFC build
    resource("cantera").stage do
      # Install Cantera build dependencies (including scons)
      system venv/"bin/pip", "install", "cython", "numpy", "ruamel.yaml", "packaging", "scons"

      # Configure Cantera build
      # Run scons with the venv's Python so it can find installed packages
      system venv/"bin/python", "-m", "SCons", "build",
             "python_package=y",
             "f90_interface=n",
             "system_sundials=y",
             "system_yamlcpp=y",
             "system_fmt=n",
             "extra_inc_dirs=#{Formula["sundials"].opt_include}:#{Formula["yaml-cpp"].opt_include}",
             "extra_lib_dirs=#{Formula["sundials"].opt_lib}:#{Formula["yaml-cpp"].opt_lib}",
             "prefix=#{libexec}/cantera",
             "python_cmd=#{venv}/bin/python",
             "-j#{ENV.make_jobs}"

      # Install Cantera
      system venv/"bin/python", "-m", "SCons", "install"

      # Install Cantera Python package into venv
      cd "build/python" do
        system venv/"bin/pip", "install", "--no-build-isolation", "."
      end
    end

    # Install Python toolchain (needed before build)
    prefix.install "toolchain"

    # Install MFC Python package and dependencies into venv
    system venv/"bin/pip", "install", "-e", prefix/"toolchain"

    # Create symlink so mfc.sh uses our pre-installed venv
    mkdir_p "build"
    ln_sf venv, "build/venv"

    # Now build MFC with pre-configured venv
    # Set VIRTUAL_ENV so mfc.sh uses existing venv instead of creating new one
    ENV["VIRTUAL_ENV"] = venv
    ENV["PATH"] = "#{venv}/bin:#{ENV.fetch("PATH", nil)}"

    system "./mfc.sh", "build",
           "-t", "pre_process", "simulation", "post_process",
           "-j", ENV.make_jobs

    # Install binaries
    # MFC installs each binary to a separate hashed subdirectory, find them individually
    %w[pre_process simulation post_process].each do |binary|
      binary_paths = Dir.glob("build/install/*/bin/#{binary}")
      raise "Could not find #{binary}" if binary_paths.empty?

      bin.install binary_paths.first
    end

    # Install mfc.sh script to libexec
    libexec.install "mfc.sh"

    # Install examples
    pkgshare.install "examples"

    # Create a wrapper that sets up a working environment for mfc.sh
    # The wrapper uses a temporary directory since Cellar is read-only and
    # activates the pre-installed Python virtual environment
    (bin/"mfc").write <<~EOS
      #!/bin/bash
      set -e

      # Activate the pre-installed Python virtual environment
      source "#{libexec}/venv/bin/activate"

      # Create a working directory for MFC in user's cache
      MFC_WORK_DIR="${TMPDIR:-/tmp}/mfc-homebrew-$$"
      mkdir -p "$MFC_WORK_DIR"

      # Function to clean up on exit
      cleanup() {
        rm -rf "$MFC_WORK_DIR"
      }
      trap cleanup EXIT

      # Create minimal directory structure that mfc.sh expects
      cd "$MFC_WORK_DIR"
      ln -sf "#{prefix}/toolchain" toolchain
      ln -sf "#{libexec}/mfc.sh" mfc.sh
      ln -sf "#{pkgshare}/examples" examples

      # Link the venv so mfc.sh doesn't try to create its own
      mkdir -p build
      ln -sf "#{libexec}/venv" build/venv

      # Set up environment variables
      export MFC_INSTALL_DIR="#{prefix}"
      export MFC_BIN_DIR="#{bin}"
      export BOOST_INCLUDE="#{Formula["boost"].opt_include}"

      # Run mfc.sh with all arguments
      exec ./mfc.sh "$@"
    EOS
    chmod 0755, bin/"mfc"
  end

  def caveats
    <<~EOS
      MFC has been installed with:
        - mfc command-line tool: #{bin}/mfc
        - pre_process: #{bin}/pre_process
        - simulation:  #{bin}/simulation
        - post_process: #{bin}/post_process

      Examples are available in:
        #{pkgshare}/examples

      To run an example:
        cd #{pkgshare}/examples/1D_sodshocktube
        mfc run case.py

      Documentation: https://mflowcode.github.io/
    EOS
  end

  test do
    # Test that the binaries exist and are executable
    assert_path_exists bin/"mfc"
    assert_predicate bin/"mfc", :executable?
    assert_path_exists bin/"pre_process"
    assert_predicate bin/"pre_process", :executable?
    assert_path_exists bin/"simulation"
    assert_predicate bin/"simulation", :executable?
    assert_path_exists bin/"post_process"
    assert_predicate bin/"post_process", :executable?

    # Verify toolchain and mfc.sh were installed
    assert_path_exists libexec/"mfc.sh"
    assert_path_exists prefix/"toolchain"
    assert_path_exists prefix/"toolchain/mfc"

    # Verify Python venv was created with dependencies
    assert_path_exists libexec/"venv"
    assert_path_exists libexec/"venv/bin/python"
    assert_path_exists libexec/"venv/bin/pip"

    # Verify examples were installed
    assert_path_exists pkgshare/"examples"
    assert_path_exists pkgshare/"examples/1D_sodshocktube/case.py"

    # Test mfc wrapper functionality with pre-installed venv
    system bin/"mfc", "--help"
  end
end
