# Homebrew Tap Deployment Checklist

Use this checklist to complete the Homebrew tap deployment setup.

## ✅ Completed (in this PR)

- [x] Updated `.github/workflows/deploy-tap.yml` with improved CD strategy
- [x] Added formula validation (`brew audit` and `brew style`)
- [x] Added support for version tag bumps
- [x] Added tap repository bootstrap logic
- [x] Changed runner to macOS (required for Homebrew tools)
- [x] Updated `HOMEBREW.md` with deployment details
- [x] Created `SETUP_INSTRUCTIONS.md` with step-by-step guide
- [x] Created `TAP_REPO_SETUP.md` with tap repository configuration
- [x] Created `PR_SUMMARY.md` with overview of changes

## 🔲 Next Steps (after PR merge)

### 1. Add Repository Secret

**In the main repository** (sbryngelson/MFC or MFlowCode/MFC after PR merge):

- [ ] Go to GitHub → (your repository) → Settings → Secrets and variables → Actions
- [ ] Click "New repository secret"
- [ ] Name: `TAP_REPO_TOKEN`
- [ ] Value: Create a Personal Access Token:
  1. GitHub → Settings → Developer settings → Personal access tokens → Tokens (classic)
  2. Click "Generate new token (classic)"
  3. Scopes: Check `repo`
  4. Copy the generated token
  5. Paste into the secret value
- [ ] Click "Add secret"

**Detailed instructions**: See `SETUP_INSTRUCTIONS.md` section 1-2

### 2. Initialize Tap Repository

**At https://github.com/MFlowCode/homebrew-mfc**:

Choose one:
- [ ] **Option A**: Click "Add a README" in GitHub UI to create first commit
- [ ] **Option B**: Let the workflow bootstrap automatically on first run

### 3. Configure Tap Repository

**At https://github.com/MFlowCode/homebrew-mfc/settings/actions**:

- [ ] Under "Workflow permissions":
  - [ ] Select "Read and write permissions"
  - [ ] Check "Allow GitHub Actions to create and approve pull requests"
- [ ] Click "Save"

### 4. Add Bottle Workflow to Tap

**Create `.github/workflows/bottle.yml` in the tap repository**:

- [ ] Copy the workflow from `TAP_REPO_SETUP.md` (section "Bottling Workflow")
- [ ] Commit to the tap repository's `main` branch

Alternatively, let the first formula sync complete, then add the workflow for subsequent updates.

### 5. Test the Deployment

**Test formula sync**:
- [ ] Edit `packaging/homebrew/mfc.rb` (make a trivial change)
- [ ] Push to your branch
- [ ] Check Actions tab - workflow should run and sync to tap
- [ ] Verify formula appears at https://github.com/MFlowCode/homebrew-mfc/blob/main/Formula/mfc.rb

**Test version bump** (optional):
- [ ] Create a test tag: `git tag v5.1.0-test && git push origin v5.1.0-test`
- [ ] Check that workflow updates URL and sha256
- [ ] Delete test tag if not needed

### 6. Verify Installation

**Test end-to-end installation**:
```bash
# Add tap
brew tap MFlowCode/mfc

# Install MFC
brew install MFlowCode/mfc/mfc

# Verify
mfc --help
which pre_process simulation post_process
```

## 📋 Reference Documentation

- **Setup Guide**: `SETUP_INSTRUCTIONS.md` - Complete setup walkthrough
- **Tap Configuration**: `TAP_REPO_SETUP.md` - Tap repository setup and bottle workflow
- **PR Overview**: `PR_SUMMARY.md` - Summary of changes and rationale
- **User Guide**: `HOMEBREW.md` - End-user installation and usage documentation

## 🔧 Troubleshooting

### Workflow fails with "Permission denied"
- Verify `TAP_REPO_TOKEN` secret is set in the main repository
- Check that the token has `repo` scope
- Ensure the token hasn't expired

### Formula fails audit
- Read the error message in the Actions log
- Common fixes:
  - Update sha256 checksum if source tarball changed
  - Fix Ruby syntax errors
  - Add missing dependencies

### Tap repository not found
- Initialize the tap repository (Step 2 above)
- Check repository name is exactly `MFlowCode/homebrew-mfc`

### Bottles don't build
- Verify Actions are enabled in tap repository (Step 3 above)
- Check that `.github/workflows/bottle.yml` exists in tap
- Review Actions logs in the tap repository

## 🎉 Success Criteria

You'll know everything is working when:
- ✅ Formula syncs automatically when you push changes
- ✅ Version tags trigger automatic formula bumps
- ✅ Bottles build for macOS 13 (x86_64) and macOS 14 (arm64)
- ✅ Users can install with: `brew install MFlowCode/mfc/mfc`
- ✅ Installation is fast (uses pre-built bottles, not compilation)

## 📞 Need Help?

See the documentation files in `packaging/homebrew/`:
- Start with `SETUP_INSTRUCTIONS.md` for step-by-step guidance
- Check `TAP_REPO_SETUP.md` for tap-specific configuration
- Review `PR_SUMMARY.md` for implementation details

