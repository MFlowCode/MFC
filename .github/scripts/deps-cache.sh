#!/bin/bash
# Provides deps_cache_init / deps_cache_restore / deps_cache_save: a cache of the
# built dependencies (LAPACK, hipfort, HDF5, Silo, FFTW) shared by every runner.
# Frontier rebuilt them from source in every job, 14-35 minutes on the login node,
# although they change only with the compiler, the modules or the pinned versions.
#
# An entry holds each dependency's install tree and its install record
# (install_manifest*.txt, which is what makes the toolchain skip it) plus its
# CMakeCache.txt (which the toolchain reads to detect a dependency built by
# another compiler). Build scratch is not kept; an entry is a few MB.
#
# A wrong entry is worse than none: a compute node cannot download sources, so
# the toolchain's rebuild of a stale dependency fails there. So the key covers
# everything the dependency build reads, an entry is published by one atomic
# rename after a successful build, and an entry that is unfinished or has lost
# files is a miss. To drop every entry, delete the cache folder.
#
# Usage: source .github/scripts/deps-cache.sh
#        deps_cache_init <cluster> "<build opts>"   # after the modules are loaded
#        deps_cache_restore                          # before ./mfc.sh build --deps-only
#        deps_cache_save                             # after it succeeds
# MFC_DEPS_CACHE_DIR overrides the folder; "off" disables the cache.

DEPS_CACHE_VERSION=1

# Number of files and their total size: the completeness check of an entry.
_deps_cache_census() {
    ( cd "$1" 2>/dev/null && find . -type f ! -name .complete -printf '%s\n' | awk '{ n++; s += $1 } END { print n + 0, s + 0 }' )
}

_deps_cache_valid() {
    [ -f "$1/.complete" ] && [ "$(_deps_cache_census "$1")" = "$(cat "$1/.complete")" ]
}

deps_cache_init() {
    local c
    deps_cache_dir="${MFC_DEPS_CACHE_DIR:-/lustre/orion/cfd154/proj-shared/runners/.deps-cache}"
    deps_cache_key=$({
        echo "v$DEPS_CACHE_VERSION cluster=$1 opts=$2 arch=$(uname -m)"
        if type module > /dev/null 2>&1; then module -t list 2>&1 | LC_ALL=C sort; fi
        for c in "${FC:-ftn}" "${CC:-cc}" "${CXX:-CC}" cmake; do
            command -v "$c" > /dev/null 2>&1 && "$c" --version 2>&1 | head -n 3
        done
        # the dependency definitions, the find modules they use, and the toolchain's dependency targets
        find toolchain/dependencies toolchain/cmake toolchain/mfc/build.py -type f 2> /dev/null | LC_ALL=C sort | xargs -r sha256sum
    } | sha256sum | cut -d' ' -f1)
    deps_cache_entry="$deps_cache_dir/$deps_cache_key"
    [ "$deps_cache_dir" = off ] && deps_cache_entry=""
    return 0
}

deps_cache_restore() {
    [ -n "${deps_cache_entry:-}" ] || return 0
    if ! _deps_cache_valid "$deps_cache_entry"; then
        echo "deps cache: miss ($deps_cache_key)"
        return 0
    fi
    if [ -e build/staging ] || [ -e build/install ]; then
        echo "deps cache: build/ is not empty; not restoring"
        return 0
    fi
    if mkdir -p build && cp -a "$deps_cache_entry/staging" build/staging && cp -a "$deps_cache_entry/install" build/install; then
        echo "deps cache: hit ($deps_cache_key): $(ls build/staging | tr '\n' ' ')"
    else
        rm -rf build/staging build/install
        echo "deps cache: miss (could not copy $deps_cache_key)"
    fi
    return 0
}

deps_cache_save() {
    [ -n "${deps_cache_entry:-}" ] || return 0
    _deps_cache_valid "$deps_cache_entry" && return 0

    # Dependencies are staged under their bare name; MFC targets under <config>-<hash>.
    local deps=() dir name
    for dir in build/staging/*/; do
        name=$(basename "$dir")
        [[ "$name" =~ ^[a-z0-9_]+$ ]] || continue
        compgen -G "${dir}install_manifest*.txt" > /dev/null || continue
        deps+=("$name")
    done
    if [ ${#deps[@]} -eq 0 ]; then
        echo "deps cache: no installed dependency to save"
        return 0
    fi

    local tmp="$deps_cache_dir/.tmp.$deps_cache_key.$$.$(hostname -s)"
    # Every step exits on failure explicitly: `set -e` is ignored in a subshell
    # that is an `if` condition, and a copy that failed half-way would otherwise be
    # counted, marked complete and published.
    if (
        umask 002
        mkdir -p "$tmp/staging" "$tmp/install" || exit 1
        for name in "${deps[@]}"; do
            mkdir "$tmp/staging/$name" || exit 1
            cp -p build/staging/"$name"/install_manifest*.txt "$tmp/staging/$name/" || exit 1
            [ ! -f "build/staging/$name/CMakeCache.txt" ] || cp -p "build/staging/$name/CMakeCache.txt" "$tmp/staging/$name/" || exit 1
            [ ! -d "build/install/$name" ] || cp -a "build/install/$name" "$tmp/install/" || exit 1
        done
        _deps_cache_census "$tmp" > "$tmp/.complete" || exit 1
        # A damaged entry would otherwise block its key for good.
        if [ -e "$deps_cache_entry" ] && ! _deps_cache_valid "$deps_cache_entry"; then
            mv -T "$deps_cache_entry" "$tmp.damaged" || exit 1
            rm -rf "$tmp.damaged"
        fi
        # rename(2) refuses a non-empty target, so the first job to publish a key wins
        mv -T "$tmp" "$deps_cache_entry" || exit 1
    ) 2> /dev/null; then
        echo "deps cache: saved ($deps_cache_key): ${deps[*]}"
    else
        rm -rf "$tmp" "$tmp.damaged"
        echo "deps cache: not saved ($deps_cache_key)"
    fi
    return 0
}
