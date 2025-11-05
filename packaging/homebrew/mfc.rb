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

  def install
    # Create Python virtual environment
    venv = libexec/"venv"
    system Formula["python@3.12"].opt_bin/"python3.12", "-m", "venv", venv
    system venv/"bin/pip", "install", "--upgrade", "pip", "setuptools", "wheel"

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

    # Build MFC using pre-configured venv
    system "./mfc.sh", "build", "-t", "pre_process", "simulation", "post_process", "-j", ENV.make_jobs.to_s

    # Install binaries - they're in hashed subdirectories like build/install/<hash>/bin/*
    Dir.glob("build/install/*/bin/*").each do |binary_path|
      bin.install binary_path
    end

    # Install main mfc.sh script
    libexec.install "mfc.sh"

    # Install toolchain directory (already done above, but make sure it stays in prefix)
    # (already done with prefix.install above)

    # Install examples
    prefix.install "examples"

    # Create smart wrapper script that:
    # 1. Works around read-only Cellar issue
    # 2. Activates venv automatically so cantera/dependencies are available
    (bin/"mfc").write <<~EOS
      #!/bin/bash
      set -e

      # Activate the pre-installed venv so all Python dependencies are available
      # This makes cantera and other packages accessible if users install them in the venv
      source "#{venv}/bin/activate"

      # Create a temporary working directory (Cellar is read-only)
      TMPDIR=$(mktemp -d)
      trap "rm -rf $TMPDIR" EXIT

      # Copy mfc.sh to temp dir (it may try to write build artifacts)
      cp "#{libexec}/mfc.sh" "$TMPDIR/"
      cd "$TMPDIR"

      # Run mfc.sh with all arguments
      exec ./mfc.sh "$@"
    EOS
  end

  def caveats
    <<~EOS
      MFC has been installed successfully!

      To use MFC:
        mfc --help

      Note: For cases requiring chemical kinetics (Cantera), you'll need to install
      Cantera separately in the MFC virtual environment:

        source #{libexec}/venv/bin/activate
        pip install cantera

      Alternatively, use conda:
        conda install -c conda-forge cantera

      Examples are available in:
        #{prefix}/examples
    EOS
  end

  test do
    # Test that all binaries exist and are executable
    %w[pre_process simulation post_process].each do |prog|
      assert_predicate bin/prog, :exist?
      assert_predicate bin/prog, :executable?
    end

    # Test that toolchain is installed
    assert_predicate prefix/"toolchain", :exist?

    # Test that venv exists and has required packages
    assert_predicate libexec/"venv", :exist?
    assert_predicate libexec/"venv/bin/python", :executable?

    # Test that examples exist
    assert_predicate prefix/"examples", :exist?

    # Test that mfc wrapper works
    system bin/"mfc", "--help"
  end
end
