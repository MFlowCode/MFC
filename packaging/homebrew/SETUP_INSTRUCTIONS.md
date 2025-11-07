# Homebrew Tap Deployment Setup Instructions

This guide walks through the complete setup for automated Homebrew tap deployment.

## 1. Create Personal Access Token (PAT)

You need a token to allow the main MFC repository to push to the tap repository.

### Steps:
1. Go to GitHub → Settings (top-right avatar) → Developer settings
2. Click "Personal access tokens" → "Tokens (classic)"
3. Click "Generate new token (classic)"
4. Configure the token:
   - **Name**: `MFC Homebrew Tap Deploy`
   - **Expiration**: Choose based on your preference (90 days, 1 year, or no expiration)
   - **Scopes**: Check only `repo` (Full control of private repositories)
5. Click "Generate token"
6. **Copy the token** - you won't be able to see it again!

## 2. Add Secret to Main Repository

Add the token as a secret in the repository where the workflow runs.

### For your fork (sbryngelson/MFC):
1. Go to your fork: https://github.com/sbryngelson/MFC
2. Click Settings → Secrets and variables → Actions
3. Click "New repository secret"
4. Configure:
   - **Name**: `TAP_REPO_TOKEN` (must match exactly)
   - **Secret**: Paste the token you copied in step 1
5. Click "Add secret"

### For the upstream repository (MFlowCode/MFC):
The same process, but the maintainers of MFlowCode/MFC will need to add the secret there when the PR is merged.

## 3. Initialize the Tap Repository

The tap repository needs to have at least one commit for the workflow to clone it.

### Option A: Via GitHub UI (easiest)
1. Go to https://github.com/MFlowCode/homebrew-mfc
2. Click "Add a README" or "creating a new file"
3. Commit the file to create the initial commit on `main`

### Option B: Automatic Bootstrap
The deploy workflow will automatically initialize an empty tap repository on first run, so you can skip this step if you prefer.

## 4. Configure Tap Repository Settings

1. Go to https://github.com/MFlowCode/homebrew-mfc/settings/actions
2. Under "Workflow permissions":
   - Select "Read and write permissions"
   - Check "Allow GitHub Actions to create and approve pull requests"
3. Click "Save"

## 5. Add Bottling Workflow to Tap (Optional but Recommended)

Create `.github/workflows/bottle.yml` in the tap repository with the content from `TAP_REPO_SETUP.md`.

This enables automatic bottle building for macOS users, so they don't have to compile from source.

## 6. Test the Setup

### Test formula sync:
1. Make a small change to `packaging/homebrew/mfc.rb` (e.g., update a comment)
2. Push to your `homebrew-formula` branch
3. Check Actions tab - the "Deploy Homebrew Tap" workflow should run
4. Verify the formula appears in https://github.com/MFlowCode/homebrew-mfc/blob/main/Formula/mfc.rb

### Test version bump:
1. Create and push a tag: `git tag v5.1.1 && git push origin v5.1.1`
2. The workflow should update the formula URL and sha256 automatically
3. The tap's bottle workflow will build bottles for the new version

## Troubleshooting

### "Permission denied" when pushing to tap
- Verify the `TAP_REPO_TOKEN` secret is set correctly
- Verify the token has `repo` scope
- Check that the token hasn't expired

### Workflow doesn't trigger
- Verify you're pushing to a tracked branch (`main`, `master`, or `homebrew-formula`)
- Check that changes include `packaging/homebrew/mfc.rb`
- For tags, ensure the tag format is `v*.*.*` (e.g., `v5.1.0`)

### "Repository not found" or clone fails
- Initialize the tap repository with at least one commit (see Step 3)
- Verify the tap repository name is exactly `MFlowCode/homebrew-mfc`

### Formula fails audit
- Check the error message in the Actions log
- Common issues:
  - Invalid URL or sha256 checksum
  - Missing or incorrect license
  - Syntax errors in the Ruby formula

## Security Notes

- The `TAP_REPO_TOKEN` should be a Personal Access Token, not an OAuth token
- Use the minimum scope needed (`repo` for private repos, or `public_repo` for public-only)
- Consider using a fine-grained PAT (beta) with access only to the tap repository
- Regularly rotate tokens (set expiration and regenerate before expiry)
- Never commit tokens to the repository - always use GitHub Secrets

## What Happens After Setup

Once configured, the system works automatically:

1. **Formula changes**: Edit `packaging/homebrew/mfc.rb` → push → automatically synced to tap
2. **New releases**: Push tag `v5.2.0` → formula auto-updates → bottles auto-build
3. **User installation**: `brew install MFlowCode/mfc/mfc` → uses pre-built bottles (fast!)

No manual intervention needed for releases or updates.

