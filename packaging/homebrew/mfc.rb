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
    bin.install "build/install/bin/pre_process"
    bin.install "build/install/bin/simulation"
    bin.install "build/install/bin/post_process"

    # Install mfc.sh script to libexec (for executable scripts)
    libexec.install "mfc.sh"

    # Install Python toolchain
    # The entire toolchain directory is required because mfc.sh depends on:
    # - util.sh for shell utilities
    # - main.py and mfc/ for the Python CLI
    # - bootstrap/ for build/lint/format scripts
    # - templates/ for HPC job submission
    prefix.install "toolchain"
    
    # Install examples
    pkgshare.install "examples"
    
    # Create a wrapper that sets up the environment and calls mfc.sh
    # The wrapper changes to the installation directory because mfc.sh
    # expects to be run from MFC's root (checks for toolchain/util.sh)
    (bin/"mfc").write <<~EOS
      #!/bin/bash
      export BOOST_INCLUDE="#{Formula["boost"].opt_include}"
      cd "#{prefix}" && exec "#{libexec}/mfc.sh" "$@"
    EOS
    chmod 0755, bin/"mfc"
  end

  def caveats
    <<~EOS
      MFC has been installed with:
        - pre_process: #{bin}/pre_process
        - simulation:  #{bin}/simulation
        - post_process: #{bin}/post_process
        - mfc wrapper: #{bin}/mfc

      Examples are available in:
        #{pkgshare}/examples

      To run an example:
        mfc run #{pkgshare}/examples/1D_sodshocktube/case.py

      Documentation: https://mflowcode.github.io/
    EOS
  end

  test do
    # Test that the binaries exist
    assert_predicate bin/"pre_process", :exist?
    assert_predicate bin/"simulation", :exist?
    assert_predicate bin/"post_process", :exist?
    
    # Test mfc wrapper
    system bin/"mfc", "--help"
    
    # Test that binaries can execute
    system bin/"pre_process", "-h"
    system bin/"simulation", "-h"
    
    # Test that mfc.sh is accessible in libexec
    assert_predicate libexec/"mfc.sh", :exist?
    assert_predicate prefix/"toolchain", :exist?
    
    # Test running a simple example case to verify full toolchain
    system bin/"mfc", "run", "#{pkgshare}/examples/1D_sodshocktube/case.py"
  end
end

