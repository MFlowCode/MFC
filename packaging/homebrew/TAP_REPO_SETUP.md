# Homebrew Tap Repository Setup

This document describes the setup needed in the `MFlowCode/homebrew-mfc` tap repository.

## Initial Setup

The tap repository needs minimal setup:

1. **Initialize the repository** (if empty):
   - Add a README.md via GitHub UI or push an initial commit
   - The deploy workflow from the main repo will bootstrap automatically if needed

2. **Enable Actions with write permissions**:
   - Go to: Settings → Actions → General → Workflow permissions
   - Select: "Read and write permissions"
   - Check: "Allow GitHub Actions to create and approve pull requests"
   - Save

## Bottling Workflow

Create `.github/workflows/bottle.yml` in the tap repository with this content:

```yaml
name: Bottle

on:
  pull_request:
    paths:
      - 'Formula/*.rb'
  push:
    branches: [ main ]
    paths:
      - 'Formula/*.rb'

permissions:
  contents: write

jobs:
  test-bot:
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [macos-14, macos-13]
    steps:
      - uses: actions/checkout@v4
      
      - uses: Homebrew/actions/setup-homebrew@master
      
      - name: Build and test (PR)
        if: github.event_name == 'pull_request'
        run: |
          brew test-bot --only-formulae --tap=$GITHUB_REPOSITORY --formula=mfc
      
      - name: Build, bottle, and publish (push to main)
        if: github.event_name == 'push'
        env:
          HOMEBREW_GITHUB_API_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: |
          brew test-bot --only-formulae --tap=$GITHUB_REPOSITORY --formula=mfc --publish
```

## How It Works

### Formula Updates (from main MFC repo)
1. Changes to `packaging/homebrew/mfc.rb` trigger the deploy workflow
2. The workflow audits and syncs the formula to the tap
3. On push to `main`, the bottle workflow runs

### Bottling (in tap repo)
1. **On PR**: `brew test-bot` builds and tests bottles for both architectures
2. **On merge to main**: 
   - `brew test-bot --publish` builds bottles
   - Uploads bottles to tap's GitHub Releases
   - Commits the `bottle do` block back to the formula

### Version Bumps (from main MFC repo)
1. Push a semver tag (e.g., `v5.2.0`) to the main MFC repo
2. The deploy workflow:
   - Updates the formula URL and sha256
   - Audits the updated formula
   - Pushes to the tap
3. The tap's bottle workflow automatically builds bottles for the new version

## No Secrets Needed

The tap repository doesn't need any secrets:
- `GITHUB_TOKEN` is sufficient for `test-bot` to commit and upload bottles
- The `TAP_REPO_TOKEN` secret only lives in the main MFC repository

## Testing

Users can install from the tap with:

```bash
# Add the tap
brew tap MFlowCode/mfc

# Install MFC
brew install MFlowCode/mfc/mfc

# Or in one command
brew install MFlowCode/mfc/mfc
```

Once bottles are published, Homebrew will automatically use them instead of building from source.

## Maintenance

- The formula is automatically synced from the main repo
- Bottles are automatically built on each formula update
- No manual intervention needed for releases

