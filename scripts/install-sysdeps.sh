#!/usr/bin/env bash
# System (apt) dependencies, determined EMPIRICALLY -- `rv sysdeps` returns [] on trixie
# (its sysreqs DB does not map the distro), so this is the captured install log. Two kinds:
#   - runtime libs the P3M trixie BINARY packages link against, found via:
#       find rv/library -name '*.so' -print0 | xargs -0 -I{} ldd {} | grep "not found"
#   - build toolchain for the Bioconductor packages compiled from source (P3M serves Bioc
#     as source on Debian): limma (Fortran + C), edgeR (C++), BiocParallel (C).
# Passwordless sudo in this env. Run BEFORE `rv sync` on a fresh clone. Idempotent.
set -euo pipefail
PKGS=(
  build-essential   # gcc / g++ / make -> compile the source Bioconductor packages
  gfortran          # limma ships Fortran sources
  libglpk40         # libglpk.so.40 runtime -> igraph (a Seurat dependency)
)
sudo apt-get update
sudo apt-get install -y --no-install-recommends "${PKGS[@]}"
