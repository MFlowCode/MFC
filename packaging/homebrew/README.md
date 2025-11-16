# MFC Homebrew Tap

Official Homebrew tap for [MFC (Multiphase Flow Code)](https://github.com/MFlowCode/MFC) – an exascale multiphase/multiphysics compressible flow solver.

## Installation

```bash
brew install mflowcode/mfc/mfc
```

## Quick Start

Run the 1D Sod shock tube example:

```bash
mkdir -p ~/mfc_example && cd ~/mfc_example
cp $(brew --prefix mfc)/examples/1D_sodshocktube/case.py .
mfc case.py -n 2
```

## What's Included

- **Prebuilt binaries**: `pre_process`, `simulation`, `post_process`
- **Python toolchain**: Pre-configured virtual environment with all dependencies
- **Examples**: Located at `$(brew --prefix mfc)/examples/`
- **`mfc` wrapper**: Simplified command-line interface for running cases

## Usage

```bash
mfc <case.py> -n <processes>
```

Use `-n X` to set the number of MPI processes.

**Note**: The Homebrew wrapper supports only running cases. For developer commands (`build`, `test`, `clean`, etc.), [clone the repository](https://github.com/MFlowCode/MFC) and use `./mfc.sh`.

## Documentation

- [MFC Documentation](https://mflowcode.github.io/)
- [Getting Started Guide](https://mflowcode.github.io/documentation/md_getting-started.html)
- [Running Cases](https://mflowcode.github.io/documentation/md_running.html)

## About This Tap

This tap is automatically maintained by the [MFC repository](https://github.com/MFlowCode/MFC) via GitHub Actions. The formula is synced from [`packaging/homebrew/mfc.rb`](https://github.com/MFlowCode/MFC/blob/master/packaging/homebrew/mfc.rb).

For issues or contributions related to the Homebrew formula, please visit the [main MFC repository](https://github.com/MFlowCode/MFC).

## Uninstall

```bash
brew uninstall mfc
brew untap mflowcode/mfc
```

## License

MFC is licensed under the [MIT License](https://github.com/MFlowCode/MFC/blob/master/LICENSE).

