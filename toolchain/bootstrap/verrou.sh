#!/bin/bash
#
# Opt-in installer for Verrou (the Valgrind FP-perturbation tool used by
# `./mfc.sh fp-stability`). Verrou is NOT a Python/pip package — it is a fork of
# Valgrind that must be compiled from source (~20 min), so this is a deliberate,
# explicit step rather than something `fp-stability` does silently.
#
#   bash toolchain/bootstrap/verrou.sh            # build into $HOME/.local/verrou
#   VERROU_HOME=/path bash toolchain/bootstrap/verrou.sh
#   bash toolchain/bootstrap/verrou.sh --force    # rebuild even if present
#
# Versions are pinned to match the fp-stability CI workflow.

set -euo pipefail

VALGRIND_VERSION="3.26.0"
VERROU_COMMIT="a58d434"
PREFIX="${VERROU_HOME:-$HOME/.local/verrou}"
FORCE="${1:-}"

echo "==> Verrou bootstrap (Valgrind ${VALGRIND_VERSION} + edf-hpc/verrou@${VERROU_COMMIT}) -> ${PREFIX}"

# Idempotent: skip if already installed and working.
if [ "$FORCE" != "--force" ] && [ -x "${PREFIX}/bin/valgrind" ] && "${PREFIX}/bin/valgrind" --tool=verrou --version >/dev/null 2>&1; then
    echo "==> Verrou already installed at ${PREFIX} (use --force to rebuild). Nothing to do."
    exit 0
fi

# Platform: Valgrind has no working modern-macOS support; Linux only.
if [ "$(uname -s)" != "Linux" ]; then
    echo "ERROR: Verrou requires Linux (Valgrind does not support modern macOS, incl. Apple Silicon)." >&2
    exit 1
fi
case "$(uname -m)" in
    x86_64) ;;
    aarch64|arm64)
        echo "WARNING: $(uname -m) detected. Valgrind builds here, but Verrou's FP backends are" >&2
        echo "         best-validated on x86_64 — treat results as experimental on this arch." >&2
        ;;
    *)
        echo "WARNING: unrecognised arch $(uname -m); the build may fail. Proceeding anyway." >&2
        ;;
esac

# Build dependencies.
missing=""
for tool in tar git make patch autoconf automake; do
    command -v "$tool" >/dev/null 2>&1 || missing="$missing $tool"
done
command -v cc >/dev/null 2>&1 || command -v gcc >/dev/null 2>&1 || missing="$missing gcc"
command -v wget >/dev/null 2>&1 || command -v curl >/dev/null 2>&1 || missing="$missing wget/curl"
if [ -n "$missing" ]; then
    echo "ERROR: missing build dependencies:$missing" >&2
    echo "       Install them (e.g. apt: build-essential automake autoconf libtool; or load HPC modules) and retry." >&2
    exit 1
fi

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT
cd "$workdir"

tarball="valgrind-${VALGRIND_VERSION}.tar.bz2"
url="https://sourceware.org/pub/valgrind/${tarball}"
echo "==> Downloading ${tarball}"
if command -v wget >/dev/null 2>&1; then
    wget -q "$url"
else
    curl -fsSL -o "$tarball" "$url"
fi
tar xf "$tarball"

echo "==> Cloning Verrou @ ${VERROU_COMMIT}"
git clone --quiet https://github.com/edf-hpc/verrou.git
git -C verrou checkout --quiet "$VERROU_COMMIT"

# Merge Verrou into the Valgrind tree and apply its patch.
cp -r verrou "valgrind-${VALGRIND_VERSION}/verrou"
cd "valgrind-${VALGRIND_VERSION}"
cat verrou/valgrind.*diff | patch -p1

echo "==> Building (this takes ~20 min)"
./autogen.sh
./configure --enable-only64bit --prefix="$PREFIX"
make -j"$(nproc)"
make install

echo "==> Verifying"
"${PREFIX}/bin/valgrind" --tool=verrou --version
echo "==> Done. Verrou installed at ${PREFIX}"
echo "    Run:  ./mfc.sh fp-stability   (or set VERROU_HOME=${PREFIX} if you used a custom prefix)"
