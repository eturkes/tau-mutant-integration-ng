#!/usr/bin/env bash
# System (apt) dependencies, determined EMPIRICALLY -- `rv sysdeps` returns [] on trixie
# (its sysreqs DB does not map the distro), so this is the captured install log. Two kinds:
#   - runtime libs the P3M trixie BINARY packages link against, found via:
#       find rv/library -name '*.so' -print0 | xargs -0 -I{} ldd {} | grep "not found"
#   - build toolchain for Bioconductor packages compiled from source (P3M serves Bioc
#     as source on Debian; current direct native deps include limma and edgeR).
# Passwordless sudo in this env. Run BEFORE `rv sync` on a fresh clone. Idempotent.
set -euo pipefail

# P3M CRAN binaries (rproject.toml) are trixie-specific; warn if the distro differs so the
# mismatch is visible (these build tools + libs exist elsewhere, so do not hard-fail).
# shellcheck source=/dev/null
. /etc/os-release 2>/dev/null || true
[ "${VERSION_CODENAME:-}" = "trixie" ] || echo "WARNING: expected Debian trixie, got '${VERSION_CODENAME:-unknown}'." >&2

PKGS=(
  build-essential   # gcc / g++ / make -> compile the source Bioconductor packages
  gfortran          # limma ships Fortran sources
  libglpk40         # libglpk.so.40 runtime -> igraph (a Seurat dependency)
)
sudo apt-get update
sudo apt-get install -y --no-install-recommends "${PKGS[@]}"

# log resolved versions to stdout (the install record; not committed -- apt versions float
# with the trixie point release, and we claim no bitwise repro without Docker).
echo "installed sysdep versions:"
dpkg-query -W -f='  ${Package} ${Version}\n' "${PKGS[@]}"
