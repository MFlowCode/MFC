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
    # Set up environment for MFC
    ENV["BOOST_INCLUDE"] = "#{Formula["boost"].opt_include}"
    ENV["FC"] = "gfortran"
    ENV["CC"] = "gcc"
    ENV["CXX"] = "g++"

    # MFC uses a Python wrapper script for building
    system "./mfc.sh", "build",
           "-t", "pre_process", "simulation", "post_process",
           "-j", ENV.make_jobs

    # Install binaries
    bin.install "build/install/bin/pre_process"
    bin.install "build/install/bin/simulation"
    bin.install "build/install/bin/post_process"

    # Install mfc.sh script to prefix
    prefix.install "mfc.sh"

    # Install Python toolchain
    prefix.install "toolchain"
    
    # Install examples
    pkgshare.install "examples"
    
    # Create a wrapper that sets up the environment and calls mfc.sh
    (bin/"mfc").write <<~EOS
      #!/bin/bash
      export BOOST_INCLUDE="#{Formula["boost"].opt_include}"
      exec "#{prefix}/mfc.sh" "$@"
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
    
    # Test that mfc.sh is accessible
    assert_predicate prefix/"mfc.sh", :exist?
    assert_predicate prefix/"toolchain", :exist?
  end
end

