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
  depends_on "python@3.12" => :build

  depends_on "boost"
  depends_on "fftw"
  depends_on "hdf5"
  depends_on "open-mpi"
  depends_on "openblas"

  def install
    # MFC uses a Python wrapper script for building
    # Homebrew's superenv handles compiler setup via gcc dependency
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

    # Install examples
    pkgshare.install "examples"
  end

  def caveats
    <<~EOS
      MFC has been installed with the following binaries:
        - pre_process: #{bin}/pre_process
        - simulation:  #{bin}/simulation
        - post_process: #{bin}/post_process

      Examples are available in:
        #{pkgshare}/examples

      For full development functionality (build, test, etc.),
      clone the repository from: https://github.com/MFlowCode/MFC

      Documentation: https://mflowcode.github.io/
    EOS
  end

  test do
    # Test that the binaries exist and are executable
    assert_path_exists bin/"pre_process"
    assert_predicate bin/"pre_process", :executable?
    assert_path_exists bin/"simulation"
    assert_predicate bin/"simulation", :executable?
    assert_path_exists bin/"post_process"
    assert_predicate bin/"post_process", :executable?

    # Verify examples were installed
    assert_path_exists pkgshare/"examples"
    assert_path_exists pkgshare/"examples/1D_sodshocktube/case.py"
  end
end
