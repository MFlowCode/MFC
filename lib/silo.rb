class Silo < Formula
  desc "LLNL mesh and field Silo I/O library. Allows creating databases for VisIt."
  homepage "https://wci.llnl.gov/simulation/computer-codes/silo"
  url "https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz"
  sha256 "4b901dfc1eb4656e83419a6fde15a2f6c6a31df84edfad7f1dc296e01b20140e"
  revision 4

  bottle do
    sha256 "624de8ce662b6ce4f51cd100faf7c212ea7ad745ee66965c6cb0c045dddc830d" => :sierra
    sha256 "0c11bc3e37dfd10984a1475a0f6c15714c58e955a28395e5548f29d2dd8d4e0b" => :el_capitan
    sha256 "392947bac24a93a7a5e70ddeb6825ada3cc97de9c7a1cb7a2b0d9806a91b5ac0" => :yosemite
    sha256 "b4cc0e28f1e29c8c2cf51b82787ec6bb1bf80310c775ad796f48a704e06e9d4c" => :x86_64_linux
  end

  option "with-static", "Build as static instead of dynamic library"
  option "without-lite-headers", "Do not install PDB lite headers"

  depends_on :x11 => :optional
  depends_on "gcc" if OS.mac? # for gfortran
  depends_on "readline"
  depends_on "hdf5" => :recommended

  def install
    args = ["--prefix=#{prefix}"]
    args << "--enable-optimization"
    args << "--with-zlib"
    args << "--enable-install-lite-headers" if build.with? "lite-headers"
    args << "--enable-shared" if build.without? "static"
    args << "--enable-x" if build.with? "x11"
    args << "--with-hdf5=#{Formula["hdf5"].opt_prefix}" if build.with? "hdf5"

    ENV.append "LDFLAGS", "-L#{Formula["readline"].opt_lib} -lreadline"
    system "./configure", *args
    system "make", "install"
    if build.with? "hdf5"
      rm lib/"libsiloh5.settings"
    else
      rm lib/"libsilo.settings"
    end
  end

  test do
    (testpath/"test.c").write <<~EOS
        #include <silo.h>

        int main(void)
        {
            DBfile *silodb;
            silodb = DBCreate("test.silo", DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
            if (!silodb)
                return 1;

            DBClose(silodb);
            return 0;
        }
        EOS
    test_args = ["test.c", "-I#{opt_include}", "-L#{opt_lib}", "-o", "test"]
    if build.with? "hdf5"
      test_args << "-lsiloh5"
      test_args << "-L#{Formula["hdf5"].opt_lib}"
      test_args << "-lhdf5"
    else
      test_args << "-lsilo"
    end
    test_args << "-lm"
    system ENV.cc, *test_args
    system "./test"
  end
end
