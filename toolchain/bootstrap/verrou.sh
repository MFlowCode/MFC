#!/bin/bash
#
# Opt-in installer for Verrou (the Valgrind FP-perturbation tool used by
# `./mfc.sh fp-stability`). Verrou is NOT a Python/pip package — it is a fork of
# Valgrind. By default this downloads a prebuilt, hash-verified artifact (seconds);
# if none is available for this tag/arch it falls back to a source build (~20 min).
# fp-stability auto-runs this on first use when Verrou is absent (printing what it
# does); it is also safe to run by hand. A failed install aborts, never a silent skip.
#
#   bash toolchain/bootstrap/verrou.sh            # install into $HOME/.local/verrou
#   VERROU_HOME=/path bash toolchain/bootstrap/verrou.sh
#   bash toolchain/bootstrap/verrou.sh --force    # reinstall even if present
#   VERROU_BUILD_FROM_SOURCE=1 bash toolchain/bootstrap/verrou.sh   # skip the prebuilt
#
# Versions are pinned to match the fp-stability CI workflow.

set -euo pipefail

VALGRIND_VERSION="3.26.0"
VERROU_COMMIT="a58d434"
# Prebuilt artifacts (built once per arch) live in a small companion repo. The tag
# pins to the (valgrind, verrou) pair above — bump all three together.
VERROU_DIST_REPO="${VERROU_DIST_REPO:-sbryngelson/verrou-dist}"
VERROU_DIST_TAG="${VERROU_DIST_TAG:-v1}"
PREFIX="${VERROU_HOME:-$HOME/.local/verrou}"
FORCE="${1:-}"

echo "==> Verrou bootstrap (Valgrind ${VALGRIND_VERSION} + edf-hpc/verrou@${VERROU_COMMIT}) -> ${PREFIX}"

# Idempotent: skip if already installed and working. Source env.sh first if present
# (a prebuilt tree needs VALGRIND_LIB to run; a source build works either way).
if [ "$FORCE" != "--force" ] && [ -x "${PREFIX}/bin/valgrind" ] \
   && ( [ -f "${PREFIX}/env.sh" ] && . "${PREFIX}/env.sh"; "${PREFIX}/bin/valgrind" --tool=verrou --version >/dev/null 2>&1 ); then
    echo "==> Verrou already installed at ${PREFIX} (use --force to rebuild). Nothing to do."
    exit 0
fi

# Platform: Valgrind has no working modern-macOS support; Linux only.
if [ "$(uname -s)" != "Linux" ]; then
    echo "ERROR: Verrou requires Linux (Valgrind does not support modern macOS, incl. Apple Silicon)." >&2
    exit 1
fi
arch_tag=""
case "$(uname -m)" in
    x86_64) arch_tag="x86_64" ;;
    aarch64|arm64)
        arch_tag="aarch64"
        echo "WARNING: $(uname -m) detected. Valgrind builds here, but Verrou's FP backends are" >&2
        echo "         best-validated on x86_64 — treat results as experimental on this arch." >&2
        ;;
    *)
        echo "WARNING: unrecognised arch $(uname -m); the build may fail. Proceeding anyway." >&2
        ;;
esac

# Fast path: download a prebuilt, hash-verified artifact and source its relocatable
# env.sh, instead of building from source. Any failure (no asset for this arch/tag,
# missing zstd/sha256sum, checksum mismatch, won't run) falls through to the build.
try_prebuilt() {
    [ -n "$arch_tag" ] || return 1
    [ "${VERROU_BUILD_FROM_SOURCE:-}" = "1" ] && return 1
    command -v sha256sum >/dev/null 2>&1 || return 1
    tar --zstd --help >/dev/null 2>&1 || command -v zstd >/dev/null 2>&1 || return 1
    command -v curl >/dev/null 2>&1 || command -v wget >/dev/null 2>&1 || return 1

    local asset base dl
    asset="verrou-${VERROU_COMMIT}-valgrind-${VALGRIND_VERSION}-linux-${arch_tag}.tar.zst"
    base="https://github.com/${VERROU_DIST_REPO}/releases/download/${VERROU_DIST_TAG}/${asset}"
    dl="$(mktemp -d)"

    echo "==> Trying prebuilt ${VERROU_DIST_REPO}@${VERROU_DIST_TAG} (${asset})"
    _fetch() {  # url dest
        if command -v curl >/dev/null 2>&1; then curl -fsSL -o "$2" "$1"; else wget -q -O "$2" "$1"; fi
    }
    if ! _fetch "$base" "$dl/$asset" || ! _fetch "$base.sha256" "$dl/$asset.sha256"; then
        echo "==> No prebuilt for this tag/arch — building from source instead."
        rm -rf "$dl"; return 1
    fi
    if ! ( cd "$dl" && sha256sum -c "$asset.sha256" >/dev/null 2>&1 ); then
        echo "WARNING: prebuilt checksum mismatch — building from source instead." >&2
        rm -rf "$dl"; return 1
    fi

    # Extract + verify in a staging dir, then swap into $PREFIX atomically. set -e
    # is suppressed inside a function used as an `if` condition, so check each step
    # explicitly — otherwise a failed extract would fall through and the source
    # build would install on top of a half-written tree (or a stale one on --force).
    local stage="$dl/stage"
    mkdir -p "$stage"
    if tar --zstd --help >/dev/null 2>&1; then
        tar -C "$stage" --zstd -xf "$dl/$asset" || { echo "WARNING: prebuilt extract failed — building from source instead." >&2; rm -rf "$dl"; return 1; }
    else
        zstd -dc "$dl/$asset" | tar -C "$stage" -xf - || { echo "WARNING: prebuilt extract failed — building from source instead." >&2; rm -rf "$dl"; return 1; }
    fi

    # Valgrind bakes its build prefix into the binary; the artifact's env.sh sets
    # VALGRIND_LIB relative to the tree so the relocated install works. Verify the
    # staged tree runs before committing it.
    if ! ( . "${stage}/env.sh" && "${stage}/bin/valgrind" --tool=verrou --version >/dev/null 2>&1 ); then
        echo "WARNING: prebuilt did not run — building from source instead." >&2
        rm -rf "$dl"; return 1
    fi

    # Commit only now: replace any existing $PREFIX atomically.
    mkdir -p "$(dirname "$PREFIX")"
    rm -rf "$PREFIX"
    if ! mv "$stage" "$PREFIX"; then
        echo "WARNING: could not install prebuilt to ${PREFIX} — building from source instead." >&2
        rm -rf "$dl"; return 1
    fi
    rm -rf "$dl"
    return 0
}

if try_prebuilt; then
    echo "==> Verifying"
    ( . "${PREFIX}/env.sh" && "${PREFIX}/bin/valgrind" --tool=verrou --version )
    echo "==> Done (prebuilt). Verrou installed at ${PREFIX}"
    echo "    Run:  ./mfc.sh fp-stability   (or set VERROU_HOME=${PREFIX} if you used a custom prefix)"
    exit 0
fi

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
