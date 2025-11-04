# README Update for Homebrew Installation

## Suggested Addition to MFC README.md

Add this section after the existing installation instructions:

---

### Quick Install (macOS/Linux)

#### Homebrew (macOS)

```bash
brew tap sbryngelson/mfc
brew install mfc
```

That's it! MFC is now installed with all dependencies. The binaries (`pre_process`, `simulation`, `post_process`) are automatically added to your PATH.

**Verify installation:**
```bash
which simulation pre_process post_process
simulation  # Should output: "./simulation.inp is missing. Exiting."
```

---

## Alternative Addition (if you want it more prominent)

Add a new "Try MFC" table row to the existing table in README.md:

```markdown
| Path | Command |
| --- | --- |
| **Codespaces** (fastest) | Click the "Codespaces" badge above to launch in 1 click |
| **Homebrew** (macOS) | `brew tap sbryngelson/mfc && brew install mfc` |
| **Local build** | `./mfc.sh build -j $(nproc) && ./mfc.sh test -j $(nproc)` |
```

---

## Documentation Update

Add a new page to the MFC documentation:

### Installation via Homebrew (macOS)

MFC can be easily installed on macOS using Homebrew.

#### Prerequisites
- macOS (Intel or Apple Silicon)
- Homebrew installed ([brew.sh](https://brew.sh))

#### Installation

```bash
# Add the MFC tap
brew tap sbryngelson/mfc

# Install MFC
brew install mfc
```

This will automatically install all required dependencies including:
- CMake
- GCC (for gfortran)
- Python 3.12
- Boost
- FFTW
- HDF5
- Open MPI
- OpenBLAS

The installation takes approximately 15-20 minutes as MFC is built from source.

#### Verification

After installation, verify that MFC is working:

```bash
# Check that binaries are in your PATH
which pre_process simulation post_process

# Test simulation binary
simulation
# Expected output: "./simulation.inp is missing. Exiting."
```

#### Using MFC

The MFC binaries are now available system-wide:

```bash
pre_process   # Generate initial conditions
simulation    # Run simulations
post_process  # Post-process results
```

Examples are included at:
```bash
/opt/homebrew/Cellar/mfc/5.1.0/share/mfc/examples/
```

#### Updating

To update to a newer version of MFC:

```bash
brew update
brew upgrade mfc
```

#### Uninstalling

```bash
brew uninstall mfc
brew untap sbryngelson/mfc
```

---

## GitHub Discussions Announcement

Suggested post for GitHub Discussions:

### 🎉 MFC is now available via Homebrew!

We're excited to announce that MFC can now be installed on macOS with a single command!

```bash
brew tap sbryngelson/mfc
brew install mfc
```

This makes it incredibly easy to get started with MFC on macOS. The Homebrew formula:
- ✅ Installs all dependencies automatically
- ✅ Builds MFC from source (optimized for your machine)
- ✅ Adds binaries to your PATH
- ✅ Includes 124 example cases
- ✅ Takes ~15-20 minutes to complete

Once installed, just run:
```bash
simulation
pre_process
post_process
```

**Tap repository**: https://github.com/sbryngelson/homebrew-mfc

This significantly lowers the barrier to entry for new users and makes MFC more accessible to the macOS community!

---

## Tweet/Social Media

Suggested post:

🚀 MFC is now on Homebrew! 

Install our exascale CFD solver on macOS with one command:

```
brew tap sbryngelson/mfc
brew install mfc
```

✅ All dependencies handled
✅ Optimized build
✅ 124 examples included

Lower barriers, more science! 🔬

#CFD #HPC #OpenSource #Homebrew #macOS

---

## Badges for README

Consider adding a Homebrew badge to the README:

```markdown
[![Homebrew](https://img.shields.io/badge/Homebrew-sbryngelson%2Fmfc-blue)](https://github.com/sbryngelson/homebrew-mfc)
```

Or

```markdown
[![Install with Homebrew](https://img.shields.io/badge/Install%20with-Homebrew-orange)](https://github.com/sbryngelson/homebrew-mfc)
```

---

## Follow-up Actions

1. ✅ Update MFC README.md with Homebrew installation instructions
2. ✅ Add Homebrew section to installation documentation
3. ✅ Post announcement in GitHub Discussions
4. ✅ Share on Twitter/social media
5. ✅ Update any installation tutorials/videos
6. ✅ Add Homebrew badge to README
7. ⏳ Consider submitting to Homebrew core (requires 75+ GitHub stars on the tap)





