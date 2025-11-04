# Homebrew Formula Guide for MFC

## What is Homebrew?

Homebrew is the most popular package manager for macOS (and Linux). Having MFC in Homebrew means users can install it with a single command:

```bash
brew install mfc
```

## Two Distribution Options

### Option 1: Personal Tap (Quick, Immediate Use) ✅

Create a GitHub repository called `homebrew-mfc` that acts as a custom Homebrew "tap".

**Advantages:**
- Immediate availability
- Full control over updates
- Users install via: `brew install sbryngelson/mfc/mfc`

**Steps:**

1. **Create the tap repository:**
```bash
# Create new GitHub repo: sbryngelson/homebrew-mfc
gh repo create homebrew-mfc --public --description "Homebrew tap for MFC"
```

2. **Push the formula:**
```bash
cd /Users/spencer/Downloads
mkdir -p homebrew-mfc/Formula
cp MFC/mfc.rb homebrew-mfc/Formula/
cd homebrew-mfc
git init
git add Formula/mfc.rb
git commit -m "Add MFC formula"
git remote add origin https://github.com/sbryngelson/homebrew-mfc.git
git push -u origin main
```

3. **Users can now install:**
```bash
brew tap sbryngelson/mfc
brew install mfc
```

---

### Option 2: Submit to homebrew-core (Official, Wider Reach)

Submit MFC to the official Homebrew repository for maximum visibility.

**Advantages:**
- Listed in official Homebrew search
- Users install via: `brew install mfc` (no tap needed)
- ~300K+ daily Homebrew users can discover MFC
- Automatic updates when MFC releases

**Disadvantages:**
- Longer review process (1-2 weeks)
- Must meet homebrew-core standards
- Less control over timing

**Steps:**

1. Fork https://github.com/Homebrew/homebrew-core
2. Add `Formula/m/mfc.rb` to your fork
3. Submit PR with title: `mfc 5.1.0 (new formula)`
4. Wait for maintainer review

---

## Current Formula Status

**File:** `/Users/spencer/Downloads/MFC/mfc.rb`

**Formula details:**
- ✅ Passes `brew audit --strict --online`
- ✅ All dependencies declared
- ✅ Includes test block
- ✅ Proper SHA256 checksum
- ✅ MIT license declared

**What it installs:**
```bash
$(brew --prefix)/bin/pre_process
$(brew --prefix)/bin/simulation  
$(brew --prefix)/bin/post_process
$(brew --prefix)/bin/mfc         # Wrapper script
```

**Dependencies handled automatically:**
- boost
- cmake
- fftw
- gcc
- hdf5
- open-mpi
- openblas
- python@3.12

---

## Testing the Formula Locally

Before publishing, test it locally:

```bash
cd /Users/spencer/Downloads/MFC

# 1. Audit the formula
brew audit --strict --online ./mfc.rb

# 2. Install from the local formula (builds from source)
brew install --build-from-source ./mfc.rb

# 3. Test the installation
mfc --help
pre_process --version
simulation --version
post_process --version

# 4. Run an example
mfc run $(brew --prefix)/share/mfc/examples/1D_sodshocktube/case.py

# 5. Uninstall when done testing
brew uninstall mfc
```

---

## Recommended Approach

**Start with Option 1 (Personal Tap):**

1. Create `sbryngelson/homebrew-mfc` tap (5 minutes)
2. Users can immediately install with: `brew install sbryngelson/mfc/mfc`
3. Add installation instructions to MFC README
4. Monitor usage for a few weeks

**Then submit Option 2 (homebrew-core):**

1. Once tap is proven stable, submit to homebrew-core
2. After merge, users can use simpler: `brew install mfc`
3. Deprecate personal tap in favor of official formula

---

## README Installation Section

Once tap is live, add this to MFC's README:

````markdown
### Installation via Homebrew (macOS/Linux)

```bash
# Add the MFC tap
brew tap sbryngelson/mfc

# Install MFC
brew install mfc

# Run an example
mfc run $(brew --prefix)/share/mfc/examples/1D_sodshocktube/case.py
```

For manual installation, see [Getting Started](https://mflowcode.github.io/documentation/md_getting-started.html).
````

---

## Maintenance

**Updating for new releases:**

1. Update `version`, `url`, and `sha256` in `mfc.rb`
2. Commit and push to tap repository
3. Users update via: `brew upgrade mfc`

**Formula location in tap:**
```
sbryngelson/homebrew-mfc/
  └── Formula/
      └── mfc.rb
```

---

## Expected Impact

**Personal tap:**
- 5-10 stars from macOS users who prefer Homebrew
- Easier onboarding for new users
- Reduced installation support requests

**homebrew-core (if accepted):**
- 20-50 stars over 6 months
- Discoverable via `brew search cfd` or `brew search flow`
- Listed on Homebrew's website
- Exposure to ~300K daily Homebrew users

---

## Next Steps

1. **Create personal tap** (Option 1) - ~5 minutes
2. **Test locally** - ~30 minutes (build time)
3. **Update README** with Homebrew instructions
4. **Announce** on Discussions/Slack
5. **Monitor** usage for 2-4 weeks
6. **Submit to homebrew-core** (Option 2) once stable

---

## Support

- Homebrew formula docs: https://docs.brew.sh/Formula-Cookbook
- Homebrew tap docs: https://docs.brew.sh/Taps
- Questions: GitHub Discussions or shb@gatech.edu






